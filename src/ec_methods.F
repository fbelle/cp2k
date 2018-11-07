!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief test and write routines for Harris energy correction on top of Kim-Gordon subsystem DFT
!>        
!> \author Fabian Belleflamme
!> \par History
!>      - 04.2018 init
! **************************************************************************************************
MODULE ec_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_allocate_matrix_set, dbcsr_copy, dbcsr_create, &
        dbcsr_deallocate_matrix_set, dbcsr_desymmetrize, dbcsr_distribution_type, dbcsr_filter, &
        dbcsr_get_info, dbcsr_init_p, dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_set, &
        dbcsr_trace, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_print, dbcsr_finalize,& 
        dbcsr_copy_into_existing, dbcsr_norm, dbcsr_transposed
   USE ec_environment_types,            ONLY: energy_correction_type
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: allocate_qs_force,&
                                              qs_force_type,&
                                              replicate_qs_force,&
                                              zero_qs_force,&
                                              deallocate_qs_force,&
                                              qs_copy_force
                                           
   !fbelle - ec_mos
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              allocate_mo_set,&
                                              get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_ks_methods,                   ONLY: calculate_w_matrix
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_core_energies,                ONLY: calculate_ptrace
   ! matrix projection
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE iterate_matrix,                  ONLY: invert_Hotelling


#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_methods'

   PUBLIC :: ec_write_forces, ec_calculate_w_matrix, ec_write_bs_forces, ec_w_matrix
   PUBLIC :: ec_mos, ec_ks_w_matrix, ec_w_test_matrix
   PUBLIC :: ec_rho_check
   PUBLIC :: matrix_projection

CONTAINS

! **************************************************************************************************
!> \brief   Simplified from qs_force.F write_forces 
!> \param   qs_env ...
!> \param   qs_force ...
!> \param   ftype 1: all , 2:Harris, 3: rho_elec, default is Harris
!> \date    04.2018
!> \author  fbelle
! **************************************************************************************************
   SUBROUTINE ec_write_forces(qs_env, qs_force, ftype)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force
      INTEGER, INTENT(IN), OPTIONAL                      :: ftype

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_write_forces', routineP = moduleN//':'//routineN

      CHARACTER(LEN=35)                                  :: fmtstr2
      CHARACTER(LEN=48)                                  :: fmtstr1
      INTEGER                                            :: i, iatom, ikind, natom, ndigits
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      REAL(KIND=dp), DIMENSION(3)                        :: grand_total
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER                                            :: my_ftype


      IF (.NOT. ASSOCIATED(qs_force)) THEN
         CALL cp_abort(__LOCATION__, &
                       "The qs_force pointer is not associated "// &
                       "and cannot be printed")
      END IF

      ! get atom mapping
      NULLIFY (atomic_kind_set)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               natom=natom)
      ALLOCATE (atom_of_kind(natom))
      ALLOCATE (kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, &
                               kind_of=kind_of)

      ! Variable precision output of the forces
      ndigits=8

      fmtstr1 = "(/,/,T2,A,/,/,T3,A,T11,A,T23,A,T40,A1,2(  X,A1))"
      WRITE (UNIT=fmtstr1(41:42), FMT="(I2)") ndigits+5

      fmtstr2 = "(/,(T2,I5,4X,I4,T18,A,T34,3F  .  ))"
      WRITE (UNIT=fmtstr2(32:33), FMT="(I2)") ndigits
      WRITE (UNIT=fmtstr2(29:30), FMT="(I2)") ndigits+6

      WRITE (UNIT=6, FMT=fmtstr1) &
         "FORCES [a.u.]", "Atom", "Kind", "Component", "X", "Y", "Z"

      grand_total(:) = 0.0_dp

      IF (PRESENT(ftype)) THEN
         my_ftype = ftype
      ELSE 
         my_ftype = 2
      END IF

      SELECT CASE (my_ftype)
      CASE DEFAULT

         DO iatom = 1, natom
            ikind = kind_of(iatom)
            i = atom_of_kind(iatom)
            WRITE (UNIT=6, FMT=fmtstr2) &
               iatom, ikind, "       overlap", qs_force(ikind)%overlap(1:3, i), &
               iatom, ikind, "       kinetic", qs_force(ikind)%kinetic(1:3, i), &
               iatom, ikind, "       gth_ppl", qs_force(ikind)%gth_ppl(1:3, i), &
               iatom, ikind, "      gth_ppnl", qs_force(ikind)%gth_ppnl(1:3, i), &
               iatom, ikind, "  core_overlap", qs_force(ikind)%core_overlap(1:3, i), &
               iatom, ikind, "      rho_core", qs_force(ikind)%rho_core(1:3, i), &
               iatom, ikind, "      rho_elec", qs_force(ikind)%rho_elec(1:3, i),&
               iatom, ikind, "         total", qs_force(ikind)%total(1:3, i)
            grand_total(1:3) = grand_total(1:3)+qs_force(ikind)%total(1:3, i)
         END DO

      CASE (0) !all

         DO iatom = 1, natom
            ikind = kind_of(iatom)
            i = atom_of_kind(iatom)
            WRITE (UNIT=6, FMT=fmtstr2) &
              iatom, ikind, "       overlap", qs_force(ikind)%overlap(1:3, i), &
              iatom, ikind, "  overlap_admm", qs_force(ikind)%overlap_admm(1:3, i), &
              iatom, ikind, "       kinetic", qs_force(ikind)%kinetic(1:3, i), &
              iatom, ikind, "       gth_ppl", qs_force(ikind)%gth_ppl(1:3, i), &
              iatom, ikind, "      gth_nlcc", qs_force(ikind)%gth_nlcc(1:3, i), &
              iatom, ikind, "      gth_ppnl", qs_force(ikind)%gth_ppnl(1:3, i), &
              iatom, ikind, "  core_overlap", qs_force(ikind)%core_overlap(1:3, i), &
              iatom, ikind, "      rho_core", qs_force(ikind)%rho_core(1:3, i), &
              iatom, ikind, "      rho_elec", qs_force(ikind)%rho_elec(1:3, i),&
              iatom, ikind, "      rho_lri_elec", qs_force(ikind)%rho_lri_elec(1:3, i), &
              iatom, ikind, "      ch_pulay", qs_force(ikind)%ch_pulay(1:3, i), &
              iatom, ikind, "    dispersion", qs_force(ikind)%dispersion(1:3, i), &
              iatom, ikind, "         other", qs_force(ikind)%other(1:3, i), &
              iatom, ikind, "       fock_4c", qs_force(ikind)%fock_4c(1:3, i), &
              iatom, ikind, "     ehrenfest", qs_force(ikind)%ehrenfest(1:3, i), &
              iatom, ikind, "        efield", qs_force(ikind)%efield(1:3, i), &
              iatom, ikind, "           eev", qs_force(ikind)%eev(1:3, i), &
              iatom, ikind, "   mp2_non_sep", qs_force(ikind)%mp2_non_sep(1:3, i), &
              iatom, ikind, "       mp2_sep", qs_force(ikind)%mp2_sep(1:3, i), &
               iatom, ikind, "         total", qs_force(ikind)%total(1:3, i)
            grand_total(1:3) = grand_total(1:3)+qs_force(ikind)%total(1:3, i)
         END DO

      CASE (2) !harris

         DO iatom = 1, natom
            ikind = kind_of(iatom)
            i = atom_of_kind(iatom)
            WRITE (UNIT=6, FMT=fmtstr2) &
              iatom, ikind, "       overlap", qs_force(ikind)%overlap(1:3, i), &
              iatom, ikind, "       kinetic", qs_force(ikind)%kinetic(1:3, i), &
              iatom, ikind, "       gth_ppl", qs_force(ikind)%gth_ppl(1:3, i), &
              iatom, ikind, "      gth_ppnl", qs_force(ikind)%gth_ppnl(1:3, i), &
              iatom, ikind, "  core_overlap", qs_force(ikind)%core_overlap(1:3, i), &
              iatom, ikind, "      rho_core", qs_force(ikind)%rho_core(1:3, i), &
              iatom, ikind, "      rho_elec", qs_force(ikind)%rho_elec(1:3, i),&
              iatom, ikind, "         total", qs_force(ikind)%total(1:3, i)
            grand_total(1:3) = grand_total(1:3)+qs_force(ikind)%total(1:3, i)
         END DO

      CASE (3) !rho_elec

         DO iatom = 1, natom
            ikind = kind_of(iatom)
            i = atom_of_kind(iatom)
            WRITE (UNIT=6, FMT=fmtstr2) &
              iatom, ikind, "      rho_elec", qs_force(ikind)%rho_elec(1:3, i)
         END DO

      END SELECT


      DEALLOCATE (atom_of_kind)
      DEALLOCATE (kind_of)

   END SUBROUTINE ec_write_forces

! **************************************************************************************************
!> \brief Compute matrix_w as needed for the forces
!>        W = P*KS*P
!> \param matrix_w ...
!> \par History
!>       2018.2 created [Fabian Belleflamme]
!> \author Fabian Belleflamme - heavily influenced by calculate_w_matrix in dm_ls_scf 
! **************************************************************************************************
   SUBROUTINE ec_calculate_w_matrix(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_calculate_w_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: scaling
      TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_tmp2, matrix_tmp3
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env = qs_env, dft_control = dft_control)
      
      ! To iterate over matrices
      nspins = dft_control%nspins

      ! allocate ec_env%matrix_w
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_w, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_w(ispin,1)%matrix)
         ! matrix_w same structure as matrix_s
         CALL dbcsr_copy(ec_env%matrix_w(ispin, 1)%matrix, ec_env%matrix_s(1, 1)%matrix , name="W MATRIX")
         CALL dbcsr_set(ec_env%matrix_w(ispin, 1)%matrix, 0.0_dp)
      END DO

      ! initialize temp matrices for multiplication
      DO ispin = 1, nspins
         CALL dbcsr_create(matrix_tmp1, template=ec_env%matrix_s(ispin, 1)%matrix, name="TMP1 MATRIX")
         CALL dbcsr_create(matrix_tmp2, template=ec_env%matrix_s(ispin, 1)%matrix, name="TMP2 MATRIX")
         CALL dbcsr_create(matrix_tmp3, template=ec_env%matrix_s(ispin, 1)%matrix, name="TMP3 MATRIX")
      END DO

         CALL dbcsr_set(matrix_tmp1, 0.0_dp)
         CALL dbcsr_set(matrix_tmp2, 0.0_dp)
         CALL dbcsr_set(matrix_tmp3, 0.0_dp)

         CALL dbcsr_finalize(matrix_tmp1)
         CALL dbcsr_finalize(matrix_tmp2)
         CALL dbcsr_finalize(matrix_tmp3)

      IF (nspins == 1) THEN
         scaling = 0.5_dp
      ELSE
         scaling = 1.0_dp
      ENDIF

      ! in the multiplication used to be ls_scf_env%eps_filter
      DO ispin = 1, nspins
         CALL dbcsr_copy(matrix_tmp3, ec_env%matrix_ks(ispin, 1)%matrix)
         CALL dbcsr_multiply("N", "N", scaling, ec_env%matrix_p(ispin, 1)%matrix, matrix_tmp3, &
                             0.0_dp, matrix_tmp1)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ec_env%matrix_p(ispin,1)%matrix, &
                             0.0_dp, matrix_tmp2)
         CALL dbcsr_copy(ec_env%matrix_w(ispin, 1)%matrix, matrix_tmp2, name="MATRIX W EC")
      END DO
      
      DO ispin = 1, nspins
         CALL dbcsr_release(matrix_tmp1)
         CALL dbcsr_release(matrix_tmp2)
         CALL dbcsr_release(matrix_tmp3)
      END DO
      
      CALL timestop(handle)

   END SUBROUTINE ec_calculate_w_matrix



! **************************************************************************************************
!> \brief Compute matrix_w as needed for the forces
!>        W = P*KS*P
!> \param matrix_w ...
!> \par History
!>       2018.2 created [Fabian Belleflamme]
!> \author Fabian Belleflamme - heavily influenced by calculate_w_matrix in dm_ls_scf 
! **************************************************************************************************
   SUBROUTINE ec_w_matrix(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_w_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: scaling
      TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_tmp2, matrix_tmp3, &
                                                            ref_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(KIND=dp)                                      :: eps_filter

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env = qs_env, dft_control = dft_control)
      
      nspins = dft_control%nspins
      eps_filter = dft_control%qs_control%eps_filter_matrix

      ! allocate ec_env%matrix_w
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_w, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_w(ispin,1)%matrix)
         ! matrix_w same structure as matrix_s
         CALL dbcsr_copy(ec_env%matrix_w(ispin, 1)%matrix, ec_env%matrix_s(1, 1)%matrix)
         CALL dbcsr_set(ec_env%matrix_w(ispin, 1)%matrix, 0.0_dp)
      END DO

      ! init temp matrices
      ref_matrix = ec_env%matrix_s(1, 1)%matrix
      CALL dbcsr_create(matrix_tmp1, template=ref_matrix, matrix_type=dbcsr_type_no_symmetry )
      CALL dbcsr_create(matrix_tmp2, template=ref_matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp3, template=ref_matrix, matrix_type=dbcsr_type_no_symmetry)

      IF (nspins == 1) THEN
         scaling = 0.5_dp
      ELSE
         scaling = 1.0_dp
      ENDIF

      ! perform multiplication W = P*KS*P
      DO ispin = 1, nspins
         CALL dbcsr_copy(matrix_tmp3, ec_env%matrix_ks(ispin, 1)%matrix)
         CALL dbcsr_filter(matrix_tmp3, eps_filter)
         CALL dbcsr_multiply("N", "N", scaling, ec_env%matrix_p(ispin, 1)%matrix, matrix_tmp3, &
                             0.0_dp, matrix_tmp1, filter_eps = eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, ec_env%matrix_p(ispin, 1)%matrix, &
                             0.0_dp, matrix_tmp2, filter_eps = eps_filter)
         CALL dbcsr_copy(ec_env%matrix_w(ispin,1)%matrix, matrix_tmp2, name="EC W MATRIX", keep_sparsity = .TRUE.)
      ENDDO

      CALL dbcsr_release(matrix_tmp1)
      CALL dbcsr_release(matrix_tmp2)
      CALL dbcsr_release(matrix_tmp3)

      CALL timestop(handle)

   END SUBROUTINE ec_w_matrix

! **************************************************************************************************
!> \brief TEST ROUTINE TO CHECK MATRIX W CALCULATION 
!>        W = P*KS*P
!> \param matrix_w ...
!> \par History
!>       2018.2 created [Fabian Belleflamme]
!> \author Fabian Belleflamme 
! **************************************************************************************************
   SUBROUTINE ec_w_test_matrix(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_w_test_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: scaling
      TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_tmp2, matrix_tmp3, &
                                                            ref_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_s, matrix_ks 
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env = qs_env, dft_control = dft_control)
      
      nspins = dft_control%nspins
      eps_filter = dft_control%qs_control%eps_filter_matrix

      
      CALL get_qs_env(qs_env=qs_env, rho = rho, matrix_ks_kp=matrix_ks, matrix_s_kp = matrix_s)
      CALL qs_rho_get(rho, rho_ao_kp = rho_ao )

      ! allocate ec_env%matrix_w
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_w, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_w(ispin,1)%matrix)
         ! matrix_w same structure as matrix_s
         CALL dbcsr_copy(ec_env%matrix_w(ispin, 1)%matrix, matrix_s(1, 1)%matrix)
         CALL dbcsr_set(ec_env%matrix_w(ispin, 1)%matrix, 0.0_dp)
      END DO

      ! init temp matrices
      ref_matrix = matrix_s(1, 1)%matrix
      CALL dbcsr_create(matrix_tmp1, template=ref_matrix, matrix_type=dbcsr_type_no_symmetry,&
         name = "temp1")
      CALL dbcsr_create(matrix_tmp2, template=ref_matrix, matrix_type=dbcsr_type_no_symmetry,&
         name = "temp2")
      CALL dbcsr_create(matrix_tmp3, template=ref_matrix, matrix_type=dbcsr_type_no_symmetry,&
         name = "temp3")

      IF (nspins == 1) THEN
         scaling = 0.5_dp
      ELSE
         scaling = 1.0_dp
      ENDIF

      ! perform multiplication W = P*KS*P
      DO ispin = 1, nspins
         CALL dbcsr_copy(matrix_tmp3, matrix_ks(ispin, 1)%matrix)
         CALL dbcsr_filter(matrix_tmp3, eps_filter)
         CALL dbcsr_multiply("N", "N", scaling, rho_ao(ispin, 1)%matrix, matrix_tmp3, &
                             0.0_dp, matrix_tmp1, filter_eps = eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, rho_ao(ispin, 1)%matrix, &
                             0.0_dp, matrix_tmp2, filter_eps = eps_filter)
     
       ! overwrite sparsity of matrix_tmp2
       CALL dbcsr_copy(ec_env%matrix_w(ispin,1)%matrix, matrix_tmp2, name="EC W MATRIX", keep_sparsity = .TRUE.)
      ENDDO

      CALL dbcsr_print(matrix_tmp3)
      CALL dbcsr_print(matrix_tmp1)
      CALL dbcsr_print(matrix_tmp2)
      CALL dbcsr_print(ec_env%matrix_w(1,1)%matrix)

      CALL dbcsr_release(matrix_tmp1)
      CALL dbcsr_release(matrix_tmp2)
      CALL dbcsr_release(matrix_tmp3)

      CALL timestop(handle)

   END SUBROUTINE ec_w_test_matrix


! **************************************************************************************************
!> \brief   Calculate and write band structure forces  
!> \param   qs_env ...
!> \param   force ...
!> \date    05.2018
!> \author  fbelle
! **************************************************************************************************
   SUBROUTINE ec_write_bs_forces(qs_env, force)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_write_bs_forces', routineP = moduleN//':'//routineN

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of, natom_of_kind
      INTEGER                                            :: natom, ikind, nkind
      TYPE(cp_para_env_type), POINTER                    :: para_env


      ! get atom mapping
      NULLIFY (atomic_kind_set , qs_force)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, natom = natom, para_env = para_env)
      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, &
                               kind_of=kind_of)

      !   *** Allocate the force data structure ***
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (natom_of_kind(nkind))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               natom_of_kind=natom_of_kind)
      CALL allocate_qs_force(qs_force, natom_of_kind)
      DEALLOCATE (natom_of_kind)
      CALL zero_qs_force(qs_force)

      ! copy to qs_force
      CALL qs_copy_force(qs_force, force)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! doesnt work past replicate

      !  *** replicate forces ***
      CALL replicate_qs_force(qs_force, para_env)

      DO ikind = 1, SIZE(force)
         qs_force(ikind)%total(:, :) = qs_force(ikind)%total(:, :)+ &
                                       qs_force(ikind)%core_overlap(:, :)+ &
                                       qs_force(ikind)%gth_ppl(:, :)+ &
                                       qs_force(ikind)%gth_ppnl(:, :)+ &
                                       qs_force(ikind)%kinetic(:, :)+ &
                                       qs_force(ikind)%overlap(:, :)+ &
                                       qs_force(ikind)%rho_core(:, :)+ &
                                       qs_force(ikind)%rho_elec(:, :)
      END DO

      CALL ec_write_forces(qs_env, qs_force, ftype = 2)

      DEALLOCATE (atom_of_kind)
      DEALLOCATE (kind_of)

      CALL deallocate_qs_force(qs_force)

   END SUBROUTINE ec_write_bs_forces


! **************************************************************************************************
!> \brief   Allocation and initialization of mo set
!> \param qs_env ...
!> \param qs_force ...
!> \date    05.2018
!> \author  fbelle
! **************************************************************************************************

   SUBROUTINE ec_mos(qs_env)
   
      TYPE(qs_environment_type), POINTER                 :: qs_env
      
      CHARACTER(len=*), PARAMETER :: routineN = 'ec_mos', routineP = moduleN//':'//routineN
      
!      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      INTEGER                                            :: ispin, nspins
   
      CALL get_qs_env(qs_env, dft_control = dft_control )
      nspins = dft_control%nspins
     



      ! IF LS-DFT : mos not allocated yet 
      IF (dft_control%qs_control%do_ls_scf) THEN
         WRITE(*,*) "ec_mos - ls_scf" 
 
         ALLOCATE (mos(nspins))
         DO ispin = 1, nspins
            NULLIFY (mos(ispin)%mo_set)
!            CALL allocate_mo_set(mo_set=mos(ispin)%mo_set, &
!                                 nao=n_ao, &
!                                 nmo=n_mo(ispin), &
!                                 nelectron=nelectron_spin(ispin), &
!                                 n_el_f=REAL(nelectron_spin(ispin), dp), &
!                                 maxocc=maxocc, &
!                                 flexible_electron_count=dft_control%relax_multiplicity)
!   
!   
!         CALL set_qs_env(qs_env, mos=mos)
   
         END DO
      
      ! If KS-DFT : mos allocated and filled 
      ELSE
         WRITE(*,*) "ec_mos - ks_scf" 
       
      END IF
   
   
   
   
  !    CALL get_qs_env(qs_env=qs_env, &
  !                    dft_control=dft_control, &
  !                    mos=mos, &
  !                    matrix_s_kp=matrix_s, &
  !                    mos_aux_fit=mos_aux_fit, &
  !                    xas_env=xas_env)
  !    CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
  ! 
  !    nmo_mat = dft_control%nspins
  !    IF (dft_control%restricted) nmo_mat = 1 ! right now, there might be more mos than needed derivs
  ! 
  ! !   *** finish initialization of the MOs ***
  !       CPASSERT(ASSOCIATED(mos))
  !       DO ispin = 1, SIZE(mos)
  !          CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b)
  !          IF (.NOT. ASSOCIATED(mo_coeff)) THEN
  !             CALL init_mo_set(mos(ispin)%mo_set, &
  !                              fm_pool=ao_mo_fm_pools(ispin)%pool, &
  !                              name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))// &
  !                              "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))))
  !          END IF
  !          IF (.NOT. ASSOCIATED(mo_coeff_b)) THEN
  !             CALL cp_fm_get_info(mos(ispin)%mo_set%mo_coeff, ncol_global=nmo)
  !             CALL dbcsr_init_p(mos(ispin)%mo_set%mo_coeff_b)
  !             CALL cp_dbcsr_m_by_n_from_row_template(mos(ispin)%mo_set%mo_coeff_b, template=matrix_s(1, 1)%matrix, n=nmo, &
  !                                                    sym=dbcsr_type_no_symmetry)
  !          ENDIF
  !       END DO
   

   END SUBROUTINE ec_mos 

! **************************************************************************************************
!> \brief Compute matrix_w from mos
!>        
!> \param matrix_w ...
!> \par History
!>       2018.5 created [Fabian Belleflamme]
!> \author Fabian Belleflamme 
! **************************************************************************************************
   SUBROUTINE ec_ks_w_matrix(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_ks_w_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      CALL timeset(routineN, handle)

      WRITE(*,*) "ec_ks_w_matrix"

      CALL get_qs_env(qs_env = qs_env, dft_control = dft_control)
      
      ! To iterate over matrices
      nspins = dft_control%nspins
      eps_filter = dft_control%qs_control%eps_filter_matrix

      ! allocate ec_env%matrix_w
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_w, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (ec_env%matrix_w(ispin,1)%matrix)
         ! matrix_w same structure as matrix_s
         CALL dbcsr_copy(ec_env%matrix_w(ispin, 1)%matrix, ec_env%matrix_s(1, 1)%matrix , name="W MATRIX")
         CALL dbcsr_set(ec_env%matrix_w(ispin, 1)%matrix, 0.0_dp)
      END DO

      
      CALL get_qs_env(qs_env, mos=mos)

      DO ispin = 1, nspins
       
         CALL calculate_w_matrix(mos(ispin)%mo_set, ec_env%matrix_w(ispin,1)%matrix)

      END DO


      CALL timestop(handle)

   END SUBROUTINE ec_ks_w_matrix

! **************************************************************************************************
!> \brief init mos environment for energy correction  
!>        
!> \param 
!> \par History
!>       2018.7 created [Fabian Belleflamme]
!> \author Fabian Belleflamme 
! **************************************************************************************************
    SUBROUTINE ec_rho_check(qs_env, matrix_p, matrix_s)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
 
      CHARACTER(len=*), PARAMETER :: routineN = 'ec_rho_check', &
         routineP = moduleN//':'//routineN
 
      INTEGER                                            :: handle, nspins
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp)                                      :: trace!, trace1
      REAL(KIND=dp)                                      :: norm_sym, norm_idempo
      TYPE(dbcsr_type)                                   :: matrix_p_trans, matrix_ps, matrix_psps
      TYPE(dbcsr_type)                                   :: ref_matrix
      TYPE(dbcsr_type)                                   :: matrix_p_desym
      TYPE(dbcsr_type)                                   :: matrix_s_desym
!      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
 
 
      CALL timeset(routineN, handle)
 
      CALL get_qs_env(qs_env = qs_env, &
                      dft_control = dft_control, &
                      para_env = para_env)

      ! setup
      nspins = dft_control%nspins
      ref_matrix = matrix_s(1,1)%matrix

      !--------------------------------------------------
      ! symmetry check
      !--------------------------------------------------
      CALL dbcsr_create(matrix_p_trans, template = ref_matrix, &
                        matrix_type = dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_p_desym, template = ref_matrix, &
                        matrix_type = dbcsr_type_no_symmetry)
    
      CALL dbcsr_desymmetrize(matrix_p(1,1)%matrix, matrix_p_desym)
      !CALL dbcsr_set(matrix_sym, 0.0_dp)
      
      ! P_trans
      CALL dbcsr_transposed(matrix_p_trans, matrix_p_desym)

      ! P_trans - P
      CALL dbcsr_add(matrix_p_desym, matrix_p_trans, alpha_scalar = 1.0_dp, beta_scalar = -1.0_dp)

      ! Norm(P_trans - P)
      CALL dbcsr_norm(matrix = matrix_p_desym, &
                      which_norm = 2, & ! dbcsr_norm_maxabsnorm
                      norm_scalar = norm_sym)

      IF(para_env%mepos==0) THEN
          WRITE(*,*) repeat("-",80)  
          WRITE(*,*) "ec_rho_check - 'Norm(P_trans - P) = 0'", norm_sym
      END IF 
                  
                  
      CALL dbcsr_set(matrix_p_desym, 0.0_dp)
      
      !--------------------------------------------------
      ! idempotency check
      !--------------------------------------------------

      CALL dbcsr_create(matrix_s_desym, template = ref_matrix, &
                        matrix_type = dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_ps, template = ref_matrix, &
                        matrix_type = dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_psps, template = ref_matrix, & 
                        matrix_type = dbcsr_type_no_symmetry)

      CALL dbcsr_set(matrix_s_desym, 0.0_dp)
      CALL dbcsr_set(matrix_ps, 0.0_dp)
      CALL dbcsr_set(matrix_psps, 0.0_dp)

      ! Desym matrix P and S
      CALL dbcsr_desymmetrize(matrix_p(1,1)%matrix, matrix_p_desym)
      CALL dbcsr_desymmetrize(matrix_s(1,1)%matrix, matrix_s_desym)

      ! PS
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p_desym, matrix_s_desym, &
                          0.0_dp, matrix_ps)


!      IF(para_env%mepos==0) THEN
!          WRITE(*,*) "ec_rho_check - matrix_ps"
!          CALL dbcsr_print(matrix_ps)
!      END IF 


      ! PSPS
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_ps, matrix_ps, &
                          0.0_dp, matrix_psps)

      ! PSPS - PS
      CALL dbcsr_add(matrix_psps, matrix_ps, &
                     alpha_scalar = 1.0_dp, beta_scalar = -1.0_dp)

      ! Norm(PSPS - PS)            
      CALL dbcsr_norm(matrix = matrix_psps, &
                      which_norm = 2, & ! dbcsr_norm_maxabsnorm
                      norm_scalar = norm_idempo)

      IF(para_env%mepos==0) THEN
          WRITE(*,*) "ec_rho_check - 'PSPS - PS = 0'", norm_idempo
      END IF 

      CALL dbcsr_release(matrix_ps)
      CALL dbcsr_release(matrix_psps)
      
      !--------------------------------------------------
      ! Particle conservation check
      !--------------------------------------------------
     
      ! Tr[PS]
      CALL calculate_ptrace(matrix_p, matrix_s, trace, nspins)

      ! \int rho_r dr
      !trace1 = pw_integrate_function(rho_r(1)%pw)


      IF(para_env%mepos==0) THEN
          WRITE(*,*) "ec_rho_check - Tr[PS]", trace
          !WRITE(*,*) "ec_rho_check - Int(rho_r)", trace1
         WRITE(*,*) repeat("-",80)  
      END IF 

      ! Clean up
      CALL dbcsr_release(matrix_p_trans)
      CALL dbcsr_release(matrix_p_desym)
      CALL dbcsr_release(matrix_s_desym)


      CALL timestop(handle)
 
    END SUBROUTINE ec_rho_check

! **************************************************************************************************
!> \brief Projection of a dbcsr_matrix from one basis set to another 
!> 
!> \param qs_env ...
!> \param matrix_in      Matrix to project
!> \param matrix_out     projected matrix in output basis
!> \param matrix_s_out   overlap matrix of output basis
!> \param sab_out_vs_in  overlapped neighborlist of in and output basis
!> \param basis_type_in  name of input basis set
!> \param basis_type_out name of output basis set

!> \date    07.2018
!> \author  fbelle
! **************************************************************************************************

   SUBROUTINE matrix_projection(qs_env, matrix_in, matrix_out, matrix_s_out, &
                                sab_out_vs_in, basis_type_in, basis_type_out)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_in
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_out
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_out
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                           ::  sab_out_vs_in
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type_in, basis_type_out

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_rho_update', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin, nspins
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_mixed
      TYPE(dbcsr_type)                                   :: matrix_s_inv, matrix_tmp
      TYPE(dbcsr_type)                                   :: matrix_A


      CALL get_qs_env(qs_env = qs_env, &
                      rho = rho, &
                      dft_control = dft_control, &
                      ks_env = ks_env)
!      CALL get_qs_env(qs_env,para_env=para_env)
      nspins = dft_control%nspins

      !--------------------------------------------------
      ! Density matrix projection ADMM2 (see Guidon2010)
      ! KG BASIS -> HARRIS BASIS
      !--------------------------------------------------

      ! calculate S_mixed
       
      NULLIFY (matrix_s_mixed)
      CALL dbcsr_allocate_matrix_set(matrix_s_mixed, nspins, 1)
    
      DO ispin = 1, nspins
         ALLOCATE (matrix_s_mixed(ispin,1)%matrix)

         CALL build_overlap_matrix(ks_env = ks_env, &
                                   matrixkp_s = matrix_s_mixed, &
                                   matrix_name = "MIXED_OVERLAP", &
                                   basis_type_a = basis_type_out, & ! HARRIS
                                   basis_type_b = basis_type_in, &  ! ORB
                                   sab_nl = sab_out_vs_in)
       END DO 

      ! ---------------------------------
      ! Calculate A = S_out^(-1) * S_mixed
      ! ---------------------------------

      ! S_out^(-1)
      CALL dbcsr_create(matrix_s_inv, template = matrix_s_out(1,1)%matrix, matrix_type = "N")
      CALL invert_Hotelling(matrix_s_inv, matrix_s_out(1,1)%matrix, & 
                            dft_control%qs_control%eps_filter_matrix)
      
      ! A = S_out^(-1) * S_mixed
      CALL dbcsr_create(matrix_A, template = matrix_s_mixed(1,1)%matrix, matrix_type = "N")
      
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_s_inv, matrix_s_mixed(1,1)%matrix, &
                          0.0_dp, matrix_A)
      
      ! Calculate P_out = A * P_in * A^T
      CALL dbcsr_create(matrix_tmp, template = matrix_A)
      DO ispin = 1, nspins
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_A, matrix_in(ispin,1)%matrix, &
                             0.0_dp, matrix_tmp)
         CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, matrix_A, &
                             0.0_dp, matrix_out(ispin,1)%matrix)
      END DO

      CALL dbcsr_release(matrix_s_inv)
      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix_A)
      CALL dbcsr_deallocate_matrix_set(matrix_s_mixed)
    
   END SUBROUTINE matrix_projection

END MODULE ec_methods
