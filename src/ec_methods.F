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
        dbcsr_trace, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_print, dbcsr_finalize
   USE kg_environment_types,            ONLY: energy_correction_type
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   !USE qs_force_types,                  ONLY: qs_force_type
   
   USE qs_force_types,                  ONLY: allocate_qs_force,&
                                              qs_force_type,&
                                              replicate_qs_force,&
                                              zero_qs_force,&
                                              deallocate_qs_force
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_methods'

   PUBLIC :: ec_write_forces, ec_calculate_w_matrix, ec_write_bs_forces

CONTAINS

! **************************************************************************************************
!> \brief   Simplified from qs_force.F write_forces 
!> \param qs_env ...
!> \param qs_force ...
!> \date    04.2018
!> \author  fbelle
! **************************************************************************************************
   SUBROUTINE ec_write_forces(qs_env, qs_force)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: qs_force

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_write_forces2', routineP = moduleN//':'//routineN

      CHARACTER(LEN=35)                                  :: fmtstr2
      CHARACTER(LEN=48)                                  :: fmtstr1
      INTEGER                                            :: i, iatom, ikind, natom, ndigits
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      REAL(KIND=dp), DIMENSION(3)                        :: grand_total
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set


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
      ndigits=6

      fmtstr1 = "(/,/,T2,A,/,/,T3,A,T11,A,T23,A,T40,A1,2(  X,A1))"
      WRITE (UNIT=fmtstr1(41:42), FMT="(I2)") ndigits+5

      fmtstr2 = "(/,(T2,I5,4X,I4,T18,A,T34,3F  .  ))"
      WRITE (UNIT=fmtstr2(32:33), FMT="(I2)") ndigits
      WRITE (UNIT=fmtstr2(29:30), FMT="(I2)") ndigits+6

      WRITE (UNIT=6, FMT=fmtstr1) &
         "FORCES [a.u.]", "Atom", "Kind", "Component", "X", "Y", "Z"

      grand_total(:) = 0.0_dp

      DO iatom = 1, natom
         ikind = kind_of(iatom)
         i = atom_of_kind(iatom)
         WRITE (UNIT=6, FMT=fmtstr2) &
            iatom, ikind, "       overlap", qs_force(ikind)%overlap(1:3, i), &
           !iatom, ikind, "  overlap_admm", qs_force(ikind)%overlap_admm(1:3, i), &
            iatom, ikind, "       kinetic", qs_force(ikind)%kinetic(1:3, i), &
            iatom, ikind, "       gth_ppl", qs_force(ikind)%gth_ppl(1:3, i), &
           !iatom, ikind, "      gth_nlcc", qs_force(ikind)%gth_nlcc(1:3, i), &
            iatom, ikind, "      gth_ppnl", qs_force(ikind)%gth_ppnl(1:3, i), &
            iatom, ikind, "  core_overlap", qs_force(ikind)%core_overlap(1:3, i), &
            iatom, ikind, "      rho_core", qs_force(ikind)%rho_core(1:3, i), &
           iatom, ikind, "      rho_elec", qs_force(ikind)%rho_elec(1:3, i),&
           !iatom, ikind, "      rho_lri_elec", qs_force(ikind)%rho_lri_elec(1:3, i), &
           !iatom, ikind, "      ch_pulay", qs_force(ikind)%ch_pulay(1:3, i), &
           !iatom, ikind, "    dispersion", qs_force(ikind)%dispersion(1:3, i), &
           !iatom, ikind, "         other", qs_force(ikind)%other(1:3, i), &
           !iatom, ikind, "       fock_4c", qs_force(ikind)%fock_4c(1:3, i), &
           !iatom, ikind, "     ehrenfest", qs_force(ikind)%ehrenfest(1:3, i), &
           !iatom, ikind, "        efield", qs_force(ikind)%efield(1:3, i), &
           !iatom, ikind, "           eev", qs_force(ikind)%eev(1:3, i), &
           !iatom, ikind, "   mp2_non_sep", qs_force(ikind)%mp2_non_sep(1:3, i), &
           !iatom, ikind, "       mp2_sep", qs_force(ikind)%mp2_sep(1:3, i), &
            iatom, ikind, "         total", qs_force(ikind)%total(1:3, i)
         grand_total(1:3) = grand_total(1:3)+qs_force(ikind)%total(1:3, i)
      END DO

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
      INTEGER                                            :: i, iatom, ikind, natom, ndigits, nkind
      CHARACTER(LEN=35)                                  :: fmtstr2
      CHARACTER(LEN=48)                                  :: fmtstr1


      ! get atom mapping
      NULLIFY (atomic_kind_set , qs_force)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, natom = natom)
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

      DO ikind = 1, SIZE(force)

         qs_force(ikind)%total(:, :) = qs_force(ikind)%total(:, :)+ &
                                       force(ikind)%core_overlap(:, :)+ &
                                       force(ikind)%gth_ppl(:, :)+ &
                                       force(ikind)%gth_ppnl(:, :)+ &
                                       force(ikind)%kinetic(:, :)+ &
                                       force(ikind)%overlap(:, :)+ &
                                       force(ikind)%rho_elec(:, :)
      END DO

      ndigits=6

      fmtstr1 = "(/,/,T2,A,/,/,T3,A,T11,A,T23,A,T40,A1,2(  X,A1))"
      WRITE (UNIT=fmtstr1(41:42), FMT="(I2)") ndigits+5

      fmtstr2 = "(/,(T2,I5,4X,I4,T18,A,T34,3F  .  ))"
      WRITE (UNIT=fmtstr2(32:33), FMT="(I2)") ndigits
      WRITE (UNIT=fmtstr2(29:30), FMT="(I2)") ndigits+6

      WRITE (UNIT=6, FMT=fmtstr1) &
         "FORCES [a.u.]", "Atom", "Kind", "Component", "X", "Y", "Z"

      DO iatom = 1, natom
         ikind = kind_of(iatom)
         i = atom_of_kind(iatom)
         WRITE (UNIT=6, FMT=fmtstr2) &
            iatom, ikind, "         f_trace", qs_force(ikind)%total(1:3, i)
      END DO

      DEALLOCATE (atom_of_kind)
      DEALLOCATE (kind_of)

      CALL deallocate_qs_force(qs_force)

   END SUBROUTINE ec_write_bs_forces

END MODULE ec_methods
