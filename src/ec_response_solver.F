!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Harris functional response force driver routine
!>
!> 
!> \date 08.2018
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_response_solver
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_copy_into_existing,&
                                              dbcsr_create,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_distribution_type,&
                                              dbcsr_init_p,&
                                              dbcsr_get_info,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_print,&
                                              dbcsr_release,&
                                              dbcsr_transposed,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_norm
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: add_qs_force,&
                                              allocate_qs_force,&
                                              qs_copy_force,&
                                              qs_force_type,&
                                              replicate_qs_force,&
                                              zero_qs_force
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho,&
                                              qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type,&
                                              qs_rho_clear
   USE qs_subsys_types,                 ONLY: qs_subsys_set,&
                                              qs_subsys_type
   USE kg_environment_types,            ONLY: energy_correction_type,&
                                              kg_environment_type
   USE kg_energy_corrections,           ONLY: ec_build_core_hamiltonian, &
                                              ec_build_ks_matrix
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integral_ab,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_write,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type,&
                                              pw_p_type
   USE ec_methods,                      ONLY: ec_write_forces,&
                                              ec_write_bs_forces,& 
                                              ec_calculate_w_matrix,&
                                              ec_w_matrix, &
                                              ec_w2_matrix, &
                                              matrix_projection,&
                                              ec_rho_check
   USE qs_integrate_potential,          ONLY: integrate_v_rspace, &
                                              integrate_v_core_rspace
   USE qs_rho_methods,                  ONLY: duplicate_rho_type                                         
   !fbelle
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE qs_core_energies,                ONLY: calculate_ecore_overlap
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE message_passing,                 ONLY: mp_sum
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_overlap,                      ONLY: build_overlap_matrix,&
                                              build_overlap_force
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE virial_types,                    ONLY: virial_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_core_energies,                ONLY: calculate_ptrace
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_triangular,&
                                              cp_fm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_invert,&
                                              cp_fm_symm,&
                                              cp_fm_upper_to_full
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE kg_tnadd_mat,                    ONLY: build_tnadd_mat
   USE input_constants,                 ONLY: kg_tnadd_atomic

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_response'

! *** Public subroutines ***

   PUBLIC :: ec_response_density
   PUBLIC :: ec_response_forces
   PUBLIC :: ec_nsc_matrix
   PUBLIC :: ec_g2_matrix

CONTAINS

! **************************************************************************************************
!> \brief  
!>          
!>          
!>          
!>          
!> \param qs_env ...
!> \date    08.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_response_density(qs_env)
 
      TYPE(qs_environment_type), POINTER                 :: qs_env
  
      CHARACTER(len=*), PARAMETER :: routineN = 'ec_reponse_density', routineP = moduleN//':'//routineN
      
      INTEGER                                            :: handle
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      INTEGER                                            :: ispin, nspins
      TYPE(dbcsr_type)                                   :: ref_matrix
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: g2_matrix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: g2_matrix_inv
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_l2
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: identity_matrix
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_tmp
      TYPE(dbcsr_type)                                   :: matrix_buf1, matrix_buf2
  
  !!---------------------------------------------------------------------------------------
  !!  Initial general setup
  !!---------------------------------------------------------------------------------------
  
      CALL timeset(routineN, handle)
      NULLIFY (logger)
      logger => cp_get_default_logger()
      
      CALL get_qs_env(qs_env,para_env=para_env)
      IF(para_env%mepos==0) THEN
         WRITE(*,*) "ec_response_density - BEGINNING"
      END IF
  
      ! get qs and kg environment
      CALL get_qs_env(qs_env = qs_env,& 
                      kg_env = kg_env,&
                      dft_control = dft_control,&
                      ks_env = ks_env)
  
      ! shortcut to ec_type
      !  ec_build_neighborlist has been built during kg_energy_correction call
      ec_env => kg_env%ec_env
      nspins = dft_control%nspins
  
      NULLIFY (g2_matrix, g2_matrix_inv, matrix_l2, identity_matrix)
  
      CALL get_qs_env(qs_env, matrix_s_kp = matrix_s)
  
      ! all matrices should be sparse here
      ref_matrix = matrix_s(1,1)%matrix
      !ref_matrix = ec_env%matrix_ks(1,1)%matrix
 
      DO ispin = 1 , nspins
         CALL dbcsr_allocate_matrix_set(g2_matrix, nspins, 1)
         ALLOCATE (g2_matrix(ispin,1)%matrix)
         CALL dbcsr_copy(g2_matrix(ispin,1)%matrix, ref_matrix, & 
                         keep_sparsity = .FALSE., name="G2-MATRIX")
         CALL dbcsr_set(g2_matrix(ispin, 1)%matrix, 0.0_dp)
  
         CALL dbcsr_allocate_matrix_set(g2_matrix_inv, nspins, 1)
         ALLOCATE (g2_matrix_inv(ispin,1)%matrix)
         CALL dbcsr_copy(g2_matrix_inv(ispin,1)%matrix, ref_matrix, & 
                         keep_sparsity = .FALSE., name="G2-INV-MATRIX")
         CALL dbcsr_set(g2_matrix_inv(ispin, 1)%matrix, 0.0_dp)

        !------  
         CALL dbcsr_allocate_matrix_set(matrix_l2, nspins, 1)
         ALLOCATE (matrix_l2(ispin,1)%matrix)
         CALL dbcsr_copy(matrix_l2(ispin,1)%matrix, ref_matrix, & 
                         name="LAMBDA2-MATRIX")
         CALL dbcsr_set(matrix_l2(ispin, 1)%matrix, 0.0_dp)
   
!         ALLOCATE (matrix_tmp(ispin,1)%matrix)
!         CALL dbcsr_create(matrix_tmp(ispin,1)%matrix, &
!                           template = ref_matrix, matrix_type = dbcsr_type_no_symmetry)
!         CALL dbcsr_set(matrix_tmp(ispin, 1)%matrix, 0.0_dp)
!--------

         CALL dbcsr_allocate_matrix_set(identity_matrix, nspins, 1)
         ALLOCATE (identity_matrix(ispin,1)%matrix)
         CALL dbcsr_copy(identity_matrix(ispin,1)%matrix, ref_matrix, & 
                         keep_sparsity = .FALSE., name="LAMBDA2-MATRIX")
         CALL dbcsr_set(identity_matrix(ispin, 1)%matrix, 0.0_dp)
      END DO
  
      ! Calculate NSC matrix
  !    CALL ec_nsc_matrix(qs_env, ec_env, nsc_matrix)
  
      IF(para_env%mepos==0) THEN
        !WRITE(*,*) " ec_response_density - nsc_matrix", ASSOCIATED(nsc_matrix)
        !IF(ASSOCIATED(nsc_matrix)) THEN
        !      CALL dbcsr_print(nsc_matrix(1,1)%matrix)
        !END IF
        WRITE(*,*) " ec_response_density - ec_env%matrix_nsc", ASSOCIATED(ec_env%matrix_nsc)
        IF(ASSOCIATED(ec_env%matrix_nsc)) THEN
              CALL dbcsr_print(ec_env%matrix_nsc(1,1)%matrix)
        END IF
      END IF
  
      ! Calculate g2_matrix
      CALL ec_g2_matrix(qs_env, ec_env, g2_matrix)
  
!      IF(para_env%mepos==0) THEN
!        WRITE(*,*) " ec_response_density - g2_matrix", ASSOCIATED(g2_matrix)
!        IF(ASSOCIATED(g2_matrix)) THEN
!              CALL dbcsr_print(g2_matrix(1,1)%matrix)
!        END IF
!      END IF
  
      ! Calculate g2_matrix^-1
      CALL matrix_inversion(qs_env, matrix_inv = g2_matrix_inv , matrix_in = g2_matrix)

      IF(para_env%mepos==0) THEN
        WRITE(*,*) " ec_response_density - g2_matrix_inv", ASSOCIATED(g2_matrix)
        IF(ASSOCIATED(g2_matrix)) THEN
              CALL dbcsr_print(g2_matrix(1,1)%matrix)
        END IF
      END IF

      !--------------------------------------------------------------------------
  
      ! test matrix inversion
      CALL dbcsr_multiply("N", "N", 1.0_dp, g2_matrix(1,1)%matrix, &
                          g2_matrix_inv(1,1)%matrix, 0.0_dp, identity_matrix(1,1)%matrix)
  
!      IF(para_env%mepos==0) THEN
!        WRITE(*,*) " ec_response_density - identity_matrix", ASSOCIATED(identity_matrix)
!        IF(ASSOCIATED(identity_matrix)) THEN
!              CALL dbcsr_print(identity_matrix(1,1)%matrix)
!        END IF
!      END IF
  
      !--------------------------------------------------------------------------

      IF(para_env%mepos==0) THEN
         WRITE(*,*) " ec_response_density - matrix_l2", ASSOCIATED(matrix_l2)
         IF(ASSOCIATED(matrix_l2)) THEN
               CALL dbcsr_print(matrix_l2(1,1)%matrix)
         END IF
      END IF

      ! Calculate response density
      ! [nsc_matrix] * [g2_matrix]^-1 = lambda2 
      DO ispin = 1 , nspins

         CALL dbcsr_create(matrix_tmp, template = ref_matrix, &
                           matrix_type = dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_buf1, template = ref_matrix, &
                           matrix_type = dbcsr_type_no_symmetry)
         CALL dbcsr_create(matrix_buf2, template = ref_matrix, &
                           matrix_type = dbcsr_type_no_symmetry)
! desym not needed
!         CALL dbcsr_desymmetrize(ec_env%matrix_nsc(ispin, 1)%matrix, matrix_buf1)
!         CALL dbcsr_desymmetrize(g2_matrix_inv(ispin, 1)%matrix, matrix_buf2)

         WRITE(*,*) "matrix_buf1"
         CALL dbcsr_print(matrix_buf1)
         WRITE(*,*) "matrix_buf2"
         CALL dbcsr_print(matrix_buf2)

         CALL dbcsr_multiply("N", "N", 1.0_dp, ec_env%matrix_nsc(ispin,1)%matrix, &
                             g2_matrix_inv(ispin,1)%matrix, 0.0_dp, matrix_tmp)
         !CALL dbcsr_multiply("N", "N", 1.0_dp, ec_env%matrix_nsc(ispin,1)%matrix, &
         !                    g2_matrix_inv(ispin,1)%matrix, 0.0_dp, matrix_l2(ispin,1)%matrix)

!         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_buf1, &
!                             matrix_buf2, 0.0_dp, matrix_tmp)

         WRITE(*,*) "matrix_tmp"
         CALL dbcsr_print(matrix_tmp)

         ! test, "make" a symmetrized response density matrix
         ! 0.5*(matrix + tranposed_matrix)
         CALL matrix_symmetrize(matrix_out = matrix_buf1, &
                                matrix_in = matrix_tmp, &
                                ref_matrix = ref_matrix)
         
         CALL dbcsr_copy(matrix_l2(ispin,1)%matrix, matrix_buf1, name="RESPONSE", keep_sparsity = .TRUE.)
      END DO
 
      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix_buf1)
      CALL dbcsr_release(matrix_buf2)

      IF(para_env%mepos==0) THEN
         WRITE(*,*) " ec_response_density - g2_matrix_inv", ASSOCIATED(g2_matrix_inv)
         IF(ASSOCIATED(g2_matrix_inv)) THEN
               CALL dbcsr_print(g2_matrix_inv(1,1)%matrix)
         END IF
!         WRITE(*,*) " ec_response_density - matrix_tmp", ASSOCIATED(matrix_tmp)
!         IF(ASSOCIATED(matrix_tmp)) THEN
!               CALL dbcsr_print(matrix_tmp(1,1)%matrix)
!         END IF
         WRITE(*,*) " ec_response_density - matrix_l2", ASSOCIATED(matrix_l2)
         IF(ASSOCIATED(matrix_l2)) THEN
               CALL dbcsr_print(matrix_l2(1,1)%matrix)
         END IF
      END IF
 
      ! calculate missing force term
      ! takes lambda2 matrix as de-facto density matrix.
      CALL ec_response_build_core_hamiltonian(qs_env = qs_env, response_matrix = matrix_l2)
      
      CALL ec_elec_response_forces(qs_env = qs_env, ec_env = ec_env, response_matrix = matrix_l2)
      
!      CALL ec_response_forces(qs_env = qs_env, ec_env = ec_env, response_matrix = matrix_l2)
  
  !    IF(para_env%mepos==0) THEN
  !      WRITE(*,*) " ec_response_density - matrix_s - AFTER ec_response_forces", ASSOCIATED(matrix_s)
  !      IF(ASSOCIATED(matrix_s)) THEN
  !            CALL dbcsr_print(matrix_s(1,1)%matrix)
  !      END IF
  !    END IF
  
       ! Deallocate dummy ks matrix                         
       CALL dbcsr_deallocate_matrix_set(g2_matrix)
       CALL dbcsr_deallocate_matrix_set(g2_matrix_inv)
       CALL dbcsr_deallocate_matrix_set(matrix_l2)
       CALL dbcsr_deallocate_matrix_set(identity_matrix)
      
       IF(para_env%mepos==0) THEN
          WRITE(*,*) "ec_response_density - END"
       END IF
      
      CALL timestop(handle)
 
   END SUBROUTINE ec_response_density

! **************************************************************************************************
!> \brief matrix symmetrize: X = 0.5*(X + X^T) + 0.5*(X - X^T) -
!>        reduction to first term of right hand side.
!> \param qs_env ...
!> \param matrix_out ...
!> \param matrix_in ...
!> \param ref_matrix ...
!> \par History
!>      08.2018 created [fbelle]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE matrix_symmetrize(matrix_out, matrix_in, ref_matrix)

      TYPE(dbcsr_type)                                   :: matrix_out, matrix_in, ref_matrix

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_inversion', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: matrix_trans

      CALL timeset(routineN, handle)


      CALL dbcsr_create(matrix_trans, template = ref_matrix, &
                        matrix_type = dbcsr_type_no_symmetry)

      ! Transpose of X
      CALL dbcsr_transposed(matrix_trans, matrix_in)

      ! P_trans - P
!      CALL dbcsr_add(matrix_in, matrix_in, alpha_scalar = 1.0_dp, beta_scalar = -0.5_dp)
!      CALL dbcsr_add(matrix_in, matrix_trans, alpha_scalar = 1.0_dp, beta_scalar = 0.5_dp)
!      CALL dbcsr_copy(matrix_out, matrix_in)

      CALL dbcsr_add(matrix_in, matrix_trans, alpha_scalar = 0.5_dp, beta_scalar = 0.5_dp)
      CALL dbcsr_copy(matrix_out, matrix_in)

      CALL dbcsr_print(matrix_out)

      CALL dbcsr_release(matrix_trans)

      CALL timestop(handle)

   END SUBROUTINE matrix_symmetrize

! **************************************************************************************************
!> \brief matrix inversion via cholesky decomposition 
!> \param qs_env ...
!> \param matrix_inv
!> \param matrix_in
!> \param 
!> \par History
!>      08.2018 created [fbelle]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE matrix_inversion(qs_env, matrix_inv, matrix_in)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_inv, matrix_in

      CHARACTER(LEN=*), PARAMETER :: routineN = 'matrix_inversion', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins, nsize
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_ks, fm_ortho
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type), POINTER                           :: ortho_dbcsr, ref_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, para_env)
      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, para_env=para_env)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      eps_filter = dft_control%qs_control%eps_filter_matrix
      nspins = dft_control%nspins

      CALL dbcsr_get_info(matrix_in(1, 1)%matrix, nfullrows_total=nsize)

      NULLIFY (fm_ortho, fm_ks, fm_struct, ref_matrix)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nsize, &
                               ncol_global=nsize, para_env=para_env)
                            
      CALL cp_fm_create(fm_ortho, fm_struct)
      CALL cp_fm_create(fm_ks, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! factorization
      ref_matrix => matrix_in(1, 1)%matrix
      NULLIFY (ortho_dbcsr)
      CALL dbcsr_init_p(ortho_dbcsr)
      CALL dbcsr_create(ortho_dbcsr, template=ref_matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      DO ispin = 1, nspins
         ref_matrix = matrix_in(ispin, 1)%matrix
         CALL copy_dbcsr_to_fm(ref_matrix, fm_ortho)
         CALL cp_fm_cholesky_decompose(fm_ortho)
         CALL cp_fm_triangular_invert(fm_ortho)
         CALL cp_fm_set_all(fm_ks, 0.0_dp)
         CALL cp_fm_to_fm_triangular(fm_ortho, fm_ks, "U")
         CALL copy_fm_to_dbcsr(fm_ks, ortho_dbcsr)
         CALL cp_fm_set_all(fm_ks, 0.0_dp)

         ! Calculate iniverse
         ! A-1 = (U^-1)(U^-1)^T
         CALL dbcsr_multiply("N", "T", 1.0_dp, ortho_dbcsr, ortho_dbcsr, &
                             0.0_dp, matrix_inv(ispin,1)%matrix, filter_eps = eps_filter)
      END DO

      CALL cp_fm_release(fm_ks)
      CALL cp_fm_release(fm_ortho)
      CALL dbcsr_release(ortho_dbcsr)
      DEALLOCATE (ortho_dbcsr)

      CALL timestop(handle)

   END SUBROUTINE matrix_inversion

! **************************************************************************************************
!> \brief response density core hamiltonian  sparse 
!>        
!> \param qs_env ...
!> \param calculate_forces ...
!> \author 8.2018 fbelle 
! **************************************************************************************************
   SUBROUTINE ec_response_build_core_hamiltonian(qs_env, response_matrix)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      !TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: response_matrix
     
      CHARACTER(len=*), PARAMETER :: routineN = 'ec_response_build_core_hamiltonian', routineP = moduleN//':'//routineN
      
      ! core hamiltonian
      INTEGER                                            :: handle, ispin, nspins, nder, nimages
  !    INTEGER                                            :: i, iatom, ikind, natom
      INTEGER                                            :: natom
      TYPE(cp_logger_type), POINTER                      :: logger
      LOGICAL                                            :: calculate_forces 
      TYPE(kg_environment_type), POINTER                 :: kg_env
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: eps_filter, eps_ppnl
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
  !    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_s, matrix_h, matrix_t, matrix_w, &
                                                            matrix_ks
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap_force
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
  
      INTEGER                                            :: ia, iatom, ikind, natom_kind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: respmat
  
  
      logger => cp_get_default_logger()
      
      CALL timeset(routineN, handle)
     
      NULLIFY(ks_env)
      NULLIFY(sab_orb, sac_ppl, sap_ppnl)
     
      CALL get_qs_env(qs_env, &
                      atomic_kind_set = atomic_kind_set, &
                      kg_env = kg_env, &
                      dft_control = dft_control, &
                      force = force, &
                      para_env = para_env, &
                      ks_env = ks_env, &
                      sab_orb = sab_orb,&
                      sac_ppl = sac_ppl, &
                      sap_ppnl = sap_ppnl)
     
      IF(para_env%mepos==0) THEN
         WRITE(*,*) "ec_response_build_core_hamiltonian - BEGINNING"
      END IF
      
      calculate_forces = .TRUE.
      nder = 1
      nimages = dft_control%nimages
      CPASSERT(nimages == 1)
      nspins = dft_control%nspins
      eps_filter = dft_control%qs_control%eps_filter_matrix
     
     !!---------------------------------------------------------------------------------------
      ! CORE FORCES 
     !!---------------------------------------------------------------------------------------
     
      CALL get_qs_env(qs_env, matrix_s_kp = matrix_s)
      CALL get_qs_env(qs_env, matrix_ks_kp = matrix_ks)
      CALL get_qs_env(qs_env, matrix_h_kp = matrix_h)
      CALL get_qs_env(qs_env, kinetic_kp = matrix_t)
     
      CALL get_qs_env(qs_env=qs_env, natom=natom)
      ALLOCATE (overlap_force(3, natom))
      overlap_force = 0.0_dp
  
      !--------------------------------------------------------------------------
      ! Energy weighed response density matrix
      !--------------------------------------------------------------------------
      NULLIFY(matrix_w)
      
      CALL dbcsr_allocate_matrix_set(matrix_w, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (matrix_w(ispin,1)%matrix)
         CALL dbcsr_create(matrix_w(ispin, 1)%matrix, template = matrix_s(1, 1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_w(ispin, 1)%matrix, sab_orb)
         CALL dbcsr_set(matrix_w(ispin, 1)%matrix, 0.0_dp)
      END DO
      
      CALL ec_w2_matrix(qs_env, matrix_w = matrix_w, matrix_p = response_matrix, matrix_ks = matrix_ks)


      CALL dbcsr_print(matrix_w(1,1)%matrix)
      CALL dbcsr_print(matrix_ks(1,1)%matrix)

      CALL build_overlap_force(ks_env, &
                               force = overlap_force, &
                               basis_type_a = "ORB", &
                               basis_type_b = "ORB", &
                               sab_nl = sab_orb, &
                               matrix_p = matrix_w(1,1)%matrix)
    
      CALL dbcsr_deallocate_matrix_set(matrix_w)
  
      respmat => response_matrix(:, 1)
      ! kinetic forces
      CALL build_kinetic_matrix(ks_env, &
                                matrixkp_t = matrix_t, &
                                !matrix_t = kg_env%tnadd_mat, &
                                matrix_name = "KINETIC ENERGY MATRIX", &
                                basis_type = "ORB", &
                                sab_nl = sab_orb, &
                                eps_filter = eps_filter, &
                                calculate_forces = calculate_forces, &
                                matrixkp_p = response_matrix)
                                !matrix_p = respmat(1)%matrix)
     
      DO ikind = 1, SIZE(atomic_kind_set, 1)
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom_kind)
         DO ia = 1, natom_kind
            iatom = atomic_kind%atom_list(ia)
            force(ikind)%overlap(:, ia) = force(ikind)%overlap(:, ia) + overlap_force(:, iatom)
         END DO
      END DO
  
      DEALLOCATE (overlap_force)
     
      ! compute the ppl contribution to the core hamiltonian
      CALL get_qs_env(qs_env = qs_env, qs_kind_set = qs_kind_set, particle_set = particle_set, &
                      atomic_kind_set = atomic_kind_set, dbcsr_dist = dbcsr_dist)
     
      NULLIFY (cell_to_index, virial)
      use_virial = .FALSE.
  
      ! forces from TNADD matrix
      IF (qs_env%kg_env%tnadd_method == kg_tnadd_atomic) THEN
      CALL build_tnadd_mat(kg_env, response_matrix, force, virial, calculate_forces, use_virial, &
                           qs_kind_set, atomic_kind_set, particle_set, sab_orb, dbcsr_dist)
      END IF
  
      IF (ASSOCIATED(sac_ppl)) THEN
         CALL build_core_ppl(matrix_h, response_matrix, force, & 
                             virial, calculate_forces, use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, & 
                             sab_orb, sac_ppl, nimages, cell_to_index, "ORB")
      END IF
  
      ! compute the ppnl contribution to the core hamiltonian ***
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         CALL build_core_ppnl(matrix_h, response_matrix, force, & 
                              virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, &
                              sab_orb, sap_ppnl, eps_ppnl, nimages, cell_to_index, "ORB")
      END IF


      IF(para_env%mepos==0) THEN
         WRITE(*,*) "ec_response_build_core_hamiltonian - END"
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_response_build_core_hamiltonian

! **************************************************************************************************
!> \brief  
!> \param qs_env ...
!>  
!>  
!> \date    08.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_elec_response_forces(qs_env, ec_env, response_matrix)
 
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: response_matrix
     
      CHARACTER(len=*), PARAMETER :: routineN = 'ec_elec_response_forces', routineP = moduleN//':'//routineN
      
      ! core hamiltonian
      INTEGER                                            :: handle, ispin, nspins, nder, nimages
  !    INTEGER                                            :: i, iatom, ikind, natom
      TYPE(cp_logger_type), POINTER                      :: logger
      LOGICAL                                            :: calculate_forces 
      TYPE(kg_environment_type), POINTER                 :: kg_env
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
  !    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial
      TYPE(cp_para_env_type), POINTER                    :: para_env
      ! electronic
      REAL(dp)                                           :: eexc, ehartree
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, tau_r, v_rspace, v_tau_rspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: input, xc_section
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: old_rho_g
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_ks
  
      logger => cp_get_default_logger()
      
      CALL timeset(routineN, handle)
     
      NULLIFY(ks_env)
      NULLIFY(sab_orb, sac_ppl, sap_ppnl)
     
      CALL get_qs_env(qs_env, &
                      atomic_kind_set = atomic_kind_set, &
                      kg_env = kg_env, &
                      dft_control = dft_control, &
                      force = force, &
                      input = input, &
                      para_env = para_env, &
                      ks_env = ks_env, &
                      sab_orb = sab_orb,&
                      sac_ppl = sac_ppl, &
                      sap_ppnl = sap_ppnl)
     
      IF(para_env%mepos==0) THEN
         WRITE(*,*) "ec_elec_response_forces - BEGINNING"
      END IF
      
      calculate_forces = .TRUE.
      nder = 1
      nimages = dft_control%nimages
      CPASSERT(nimages == 1)
      nspins = dft_control%nspins
      eps_filter = dft_control%qs_control%eps_filter_matrix
     
      NULLIFY(xc_section)
  
      xc_section => section_vals_get_subs_vals(input, "DFT%KG_METHOD%XC")
  
      CALL get_qs_env(qs_env=qs_env, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      CPASSERT(.NOT. use_virial)
  
      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))
  
      NULLIFY (auxbas_pw_pool, poisson_env, pw_pools)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools, poisson_env=poisson_env)
  
      ! Calculate the Hartree potential
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             v_hartree_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
  
      ! get sparse density 
      CALL qs_rho_get(ec_env%old_rho, rho_g = old_rho_g)
  
      ! + v_hartree[rho_in]
      CALL pw_zero(rho_tot_gspace%pw)
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, ec_env%old_rho)
  
      ! Hartree potential
      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree, &
                            v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)
  
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
  
      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      NULLIFY (v_rspace, v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=ec_env%old_rho, xc_section=xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=eexc, just_energy=.FALSE.)
  
      CALL qs_rho_get(ec_env%old_rho, rho_r=rho_r)
      IF (ASSOCIATED(v_tau_rspace)) THEN
         CALL qs_rho_get(ec_env%old_rho, tau_r=tau_r)
      END IF

      CALL get_qs_env(qs_env, matrix_ks_kp = matrix_ks)

      DO ispin = 1, nspins
         ! Add v_hartree + v_xc = v_rspace
         CALL pw_scale(v_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(v_hartree_rspace%pw, v_rspace(ispin)%pw)
  
         CALL integrate_v_rspace(v_rspace = v_rspace(ispin), &
                                 pmat = response_matrix(ispin, 1), &
                                 hmat = matrix_ks(ispin, 1), &
                                 qs_env = qs_env, &
                                 calculate_forces = calculate_forces, &
                                 basis_type = "ORB")
        
         IF (ASSOCIATED(v_tau_rspace)) THEN
            ! integrate over Tau-potential <nabla.a|V|nabla.b>
            CALL pw_scale(v_tau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin), hmat = matrix_ks(ispin, 1), &
                                    qs_env=qs_env, calculate_forces = calculate_forces, &
                                    compute_tau=.TRUE., &
                                    basis_type="ORB")
         END IF
  
      END DO
  
      ! return pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
         IF (ASSOCIATED(v_tau_rspace)) THEN
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
         END IF
      ENDDO
  
!      CALL calculate_ecore_overlap(qs_env, para_env, calculate_forces=calculate_forces)
  
      DEALLOCATE (v_rspace)
  
      CALL timestop(handle)
  
      IF(para_env%mepos==0) THEN
         WRITE(*,*) "ec_elec_response_forces - END"
      END IF
 
   END SUBROUTINE ec_elec_response_forces

! **************************************************************************************************
!> \brief  
!> \param qs_env ...
!>  
!>  
!> \date    08.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_response_forces(qs_env, ec_env, response_matrix)
 
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: response_matrix
     
      CHARACTER(len=*), PARAMETER :: routineN = 'ec_response_forces', routineP = moduleN//':'//routineN
      
      ! core hamiltonian
      INTEGER                                            :: handle, ispin, nspins, nder, nimages
  !    INTEGER                                            :: i, iatom, ikind, natom
      INTEGER                                            :: natom
      TYPE(cp_logger_type), POINTER                      :: logger
      LOGICAL                                            :: calculate_forces 
      TYPE(kg_environment_type), POINTER                 :: kg_env
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: eps_filter, eps_ppnl
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
  !    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_s, matrix_h, matrix_t, matrix_w, &
                                                            matrix_ks
      ! electronic
      REAL(dp)                                           :: eexc, ehartree
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, tau_r, v_rspace, v_tau_rspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: input, xc_section
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: old_rho_g
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: overlap_force
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
  
      INTEGER                                            :: ia, iatom, ikind, natom_kind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: respmat
  
  
      logger => cp_get_default_logger()
      
      CALL timeset(routineN, handle)
     
      NULLIFY(ks_env)
      NULLIFY(sab_orb, sac_ppl, sap_ppnl)
     
      CALL get_qs_env(qs_env, &
                      atomic_kind_set = atomic_kind_set, &
                      kg_env = kg_env, &
                      dft_control = dft_control, &
                      force = force, &
                      input = input, &
                      para_env = para_env, &
                      ks_env = ks_env, &
                      sab_orb = sab_orb,&
                      sac_ppl = sac_ppl, &
                      sap_ppnl = sap_ppnl)
     
      IF(para_env%mepos==0) THEN
         WRITE(*,*) "ec_response_forces - BEGINNING"
      END IF
      
      calculate_forces = .TRUE.
      nder = 1
      nimages = dft_control%nimages
      CPASSERT(nimages == 1)
      nspins = dft_control%nspins
      eps_filter = dft_control%qs_control%eps_filter_matrix
     
     !!---------------------------------------------------------------------------------------
      ! CORE FORCES 
     !!---------------------------------------------------------------------------------------
     
      CALL get_qs_env(qs_env, matrix_s_kp = matrix_s)
      CALL get_qs_env(qs_env, matrix_ks_kp = matrix_ks)
      CALL get_qs_env(qs_env, matrix_h_kp = matrix_h)
      CALL get_qs_env(qs_env, kinetic_kp = matrix_t)
     
      CALL get_qs_env(qs_env=qs_env, natom=natom)
      ALLOCATE (overlap_force(3, natom))
      overlap_force = 0.0_dp
  
      !--------------------------------------------------------------------------
      ! Energy weighed response density matrix
      !--------------------------------------------------------------------------
      NULLIFY(matrix_w)
      
      CALL dbcsr_allocate_matrix_set(matrix_w, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (matrix_w(ispin,1)%matrix)
         !CALL dbcsr_copy(matrix_w(ispin, 1)%matrix, matrix_s(1, 1)%matrix)
         !CALL dbcsr_set(matrix_w(ispin, 1)%matrix, 0.0_dp)
         CALL dbcsr_create(matrix_w(ispin, 1)%matrix, template = matrix_s(1, 1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_w(ispin, 1)%matrix, sab_orb)
         CALL dbcsr_set(matrix_w(ispin, 1)%matrix, 0.0_dp)
      END DO
      
      CALL ec_w2_matrix(qs_env, matrix_w = matrix_w, matrix_p = response_matrix, matrix_ks = matrix_ks)
  
      CALL build_overlap_force(ks_env, &
                               force = overlap_force, &
                               basis_type_a = "ORB", &
                               basis_type_b = "ORB", &
                               sab_nl = sab_orb, &
                               matrix_p = matrix_w(1,1)%matrix)
    
      CALL dbcsr_deallocate_matrix_set(matrix_w)
  
      respmat => response_matrix(:, 1)
  
      ! kinetic forces
      CALL build_kinetic_matrix(ks_env, &
  !                              matrixkp_t = matrix_t, &
                                matrix_t = kg_env%tnadd_mat, &
                                matrix_name = "KINETIC ENERGY MATRIX", &
                                basis_type = "ORB", &
                                sab_nl = sab_orb, &
                                eps_filter = eps_filter, &
                                calculate_forces = calculate_forces, &
                                matrix_p = respmat(1)%matrix)
     
      DO ikind = 1, SIZE(atomic_kind_set, 1)
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom_kind)
         DO ia = 1, natom_kind
            iatom = atomic_kind%atom_list(ia)
            force(ikind)%overlap(:, ia) = force(ikind)%overlap(:, ia) + overlap_force(:, iatom)
         END DO
      END DO
  
      DEALLOCATE (overlap_force)
     
      ! compute the ppl contribution to the core hamiltonian
      CALL get_qs_env(qs_env = qs_env, qs_kind_set = qs_kind_set, particle_set = particle_set, &
                      atomic_kind_set = atomic_kind_set, dbcsr_dist = dbcsr_dist)
     
      NULLIFY (cell_to_index, virial)
      use_virial = .FALSE.
  
      ! forces from TNADD matrix
      IF (qs_env%kg_env%tnadd_method == kg_tnadd_atomic) THEN
  !    CALL build_tnadd_mat(kg_env, response_matrix, force, virial, calculate_forces, use_virial, &
  !                         qs_kind_set, atomic_kind_set, particle_set, sab_orb, dbcsr_dist)
      END IF
  
      IF (ASSOCIATED(sac_ppl)) THEN
         CALL build_core_ppl(matrix_h, response_matrix, force, & 
                             virial, calculate_forces, use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, & 
                             sab_orb, sac_ppl, nimages, cell_to_index, "ORB")
      END IF
  
      ! compute the ppnl contribution to the core hamiltonian ***
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         CALL build_core_ppnl(matrix_h, response_matrix, force, & 
                              virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, &
                              sab_orb, sap_ppnl, eps_ppnl, nimages, cell_to_index, "ORB")
      END IF

!!---------------------------------------------------------------------------------------
    ! ELECTRONIC FORCES 
!!---------------------------------------------------------------------------------------

      NULLIFY(xc_section)
  
      xc_section => section_vals_get_subs_vals(input, "DFT%KG_METHOD%XC")
  
      CALL get_qs_env(qs_env=qs_env, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      CPASSERT(.NOT. use_virial)
  
      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))
  
      NULLIFY (auxbas_pw_pool, poisson_env, pw_pools)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools, poisson_env=poisson_env)
  
      ! Calculate the Hartree potential
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, &
                             in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, &
                             v_hartree_rspace%pw, &
                             use_data=REALDATA3D, &
                             in_space=REALSPACE)
  
      ! get sparse density 
      CALL qs_rho_get(ec_env%old_rho, rho_g = old_rho_g)
  
      ! + v_hartree[rho_in]
      CALL pw_zero(rho_tot_gspace%pw)
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, ec_env%old_rho)
  
      ! Hartree potential
      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree, &
                            v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)
  
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
  
      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      NULLIFY (v_rspace, v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=ec_env%old_rho, xc_section=xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=eexc, just_energy=.FALSE.)
  
      CALL qs_rho_get(ec_env%old_rho, rho_r=rho_r)
      IF (ASSOCIATED(v_tau_rspace)) THEN
         CALL qs_rho_get(ec_env%old_rho, tau_r=tau_r)
      END IF
  
      DO ispin = 1, nspins
         ! Add v_hartree + v_xc = v_rspace
         CALL pw_scale(v_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(v_hartree_rspace%pw, v_rspace(ispin)%pw)
  
         CALL integrate_v_rspace(v_rspace = v_rspace(ispin), &
                                 pmat = response_matrix(ispin, 1), &
                                 hmat = matrix_ks(ispin, 1), &
                                 qs_env = qs_env, &
                                 calculate_forces = calculate_forces, &
                                 basis_type = "ORB")
        
         IF (ASSOCIATED(v_tau_rspace)) THEN
            ! integrate over Tau-potential <nabla.a|V|nabla.b>
            CALL pw_scale(v_tau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin), hmat = matrix_ks(ispin, 1), &
                                    qs_env=qs_env, calculate_forces = calculate_forces, &
                                    compute_tau=.TRUE., &
                                    basis_type="ORB")
         END IF
  
      END DO
  
      ! return pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
         IF (ASSOCIATED(v_tau_rspace)) THEN
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
         END IF
      ENDDO
  
      ! energies
  !    CALL calculate_ecore_overlap(qs_env, para_env, calculate_forces=calculate_forces)
  
      DEALLOCATE (v_rspace)
  
      CALL timestop(handle)
  
      IF(para_env%mepos==0) THEN
         WRITE(*,*) "ec_response_forces - END"
      END IF
 
   END SUBROUTINE ec_response_forces

! **************************************************************************************************
!> \brief  Non-self-consistent forces (Hartree and XC part)
!> \param qs_env ...
!>  adapted for constraint on lagrangian - matrix calculation 
!>  \int (v_hartree[rho_out-rho_in] + (rho_out - rho_in)*\mu'_[rho_in]) \nabla_I rho_in dr

!> \date    08.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_nsc_matrix(qs_env, ec_env, nsc_matrix)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(energy_correction_type), POINTER              :: ec_env
    TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: nsc_matrix

    CHARACTER(len=*), PARAMETER :: routineN = 'ec_nsc_matrix', routineP = moduleN//':'//routineN
    
    TYPE(kg_environment_type), POINTER                 :: kg_env
    INTEGER                                            :: handle
    TYPE(cp_logger_type), POINTER                      :: logger
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
    ! fbelle
    TYPE(pw_env_type), POINTER                         :: pw_env
    TYPE(qs_rho_type), POINTER                         :: rho!, old_rho
    TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                          v_hartree_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, old_rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, old_rho_g
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao, old_rho_ao
    TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
    TYPE(pw_poisson_type), POINTER                     :: poisson_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    INTEGER                                            :: ispin, nspins
    LOGICAL                                            :: calculate_forces 
    TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
    INTEGER                                            :: natom
    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
    ! fbelle
    TYPE(xc_derivative_set_type), POINTER              :: deriv_set
    TYPE(xc_rho_set_type), POINTER                     :: rho_set
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc, v_rspace_nsc
    TYPE(section_vals_type), POINTER                   :: input, xc_section

    logger => cp_get_default_logger()
    
    CALL timeset(routineN, handle)
    
    CALL get_qs_env(qs_env, &
                    kg_env = kg_env, &
                    dft_control = dft_control, &
                    rho = rho, &
                    force = force, &
                    input = input, &
                    para_env = para_env)

    IF(para_env%mepos==0) THEN
       WRITE(*,*) "ec_nsc_matrix - BEGINNING"
    END IF
    
    calculate_forces = .FALSE.
    
    nspins = dft_control%nspins
     
    ! get atomic mappng
    CALL get_qs_env(qs_env, & 
                    atomic_kind_set = atomic_kind_set, &
                    natom = natom)
    ALLOCATE (atom_of_kind(natom), kind_of(natom))
    CALL get_atomic_kind_set(atomic_kind_set, &
                             atom_of_kind = atom_of_kind, &
                             kind_of = kind_of)

!!---------------------------------------------------------------------------------------
    ! *** Calculation of forces due to non-self consistency of energy correction *** ! 
    ! First, setup pw pools
!!---------------------------------------------------------------------------------------

    NULLIFY (pw_env)
    CALL get_qs_env(qs_env, pw_env=pw_env)
    CPASSERT(ASSOCIATED(pw_env))

    NULLIFY (auxbas_pw_pool, poisson_env, pw_pools)
    ! gets the tmp grids
    CALL pw_env_get(pw_env = pw_env, &
                    auxbas_pw_pool = auxbas_pw_pool, &
                    pw_pools = pw_pools, &
                    poisson_env = poisson_env)

    ! Calculate the NSC Hartree potential
    CALL pw_pool_create_pw(pool = auxbas_pw_pool, &
                           pw = v_hartree_gspace%pw, &
                           use_data = COMPLEXDATA1D, &
                           in_space = RECIPROCALSPACE)
    CALL pw_pool_create_pw(pool = auxbas_pw_pool, &
                           pw = rho_tot_gspace%pw, &
                           use_data = COMPLEXDATA1D, &
                           in_space = RECIPROCALSPACE)

    CALL pw_pool_create_pw(pool = auxbas_pw_pool, &
                           pw = v_hartree_rspace%pw, &
                           use_data = REALDATA3D, &
                           in_space = REALSPACE)

    ! get rho_out - after Harris energy correction
    CALL qs_rho_get(rho, & 
                    rho_r = rho_r, &
                    rho_g = rho_g, &
                    rho_ao_kp = rho_ao)

    ! get rho_in - before Harris energy correction
    CALL qs_rho_get(ec_env%old_rho, &
                    rho_r = old_rho_r, &
                    rho_g = old_rho_g, &
                    rho_ao_kp = old_rho_ao )

!!---------------------------------------------------------------------------------------
    ! *** Calculation of forces due to non-self consistency of energy correction *** ! 
    ! First, v_NSC_xc[rho_in]
!!---------------------------------------------------------------------------------------

    ! ** v_xc ** !
    ALLOCATE(v_xc(nspins), v_rspace_nsc(nspins))
    DO ispin = 1, nspins
       NULLIFY(v_xc(nspins)%pw, v_rspace_nsc(nspins)%pw)
       CALL pw_pool_create_pw(pool = auxbas_pw_pool, & 
                              pw = v_xc(ispin)%pw, &
                              use_data = REALDATA3D, &
                              in_space = REALSPACE)
       CALL pw_pool_create_pw(pool = auxbas_pw_pool, & 
                              pw = v_rspace_nsc(ispin)%pw, &
                              use_data = REALDATA3D, &
                              in_space = REALSPACE)
       CALL pw_zero(v_xc(ispin)%pw)
       CALL pw_zero(v_rspace_nsc(nspins)%pw)
    END DO


    NULLIFY(xc_section)

    xc_section => section_vals_get_subs_vals(input, "DFT%KG_METHOD%XC")

    NULLIFY (rho_set, deriv_set)
    ! main ingredient is xc_rho_set_and_dset_create
    CALL xc_prep_2nd_deriv(deriv_set = deriv_set, &        ! containing potentials
                           rho_set = rho_set, &            ! density at which deriv are calced
                           rho_r = old_rho_r, &            ! where derivative evaluated 
                           pw_pool = auxbas_pw_pool, &     ! OK
                           xc_section = ec_env%xc_section)        ! OK 

    CALL xc_calc_2nd_deriv(v_xc = v_xc, & 
                           deriv_set = deriv_set, &
                           rho_set = rho_set, &
                           rho1_set = rho_set, &
                           pw_pool = auxbas_pw_pool, &
                           xc_section = ec_env%xc_section, & 
                           gapw = .FALSE.)

    DO ispin = 1,nspins
      ! \nu'(rho_in)*(rho_out - rho_in)
      v_rspace_nsc(ispin)%pw%cr3d = (rho_r(ispin)%pw%cr3d &
                                    - old_rho_r(ispin)%pw%cr3d) &
                                    * v_xc(ispin)%pw%cr3d
    END DO

    CALL xc_dset_release(deriv_set)
    CALL xc_rho_set_release(rho_set = rho_set, pw_pool = auxbas_pw_pool)

    ! + v_hartree[rho_out-rho_in]

    ! rho_tot_gspace = rho_g - old_rho_g
    CALL pw_zero(rho_tot_gspace%pw)
    DO ispin = 1,nspins
      CALL pw_axpy(rho_g(ispin)%pw,rho_tot_gspace%pw)
      CALL pw_axpy(old_rho_g(ispin)%pw, rho_tot_gspace%pw, alpha = -1._dp)
    END DO

    ! get Hartree potential from rho_tot_gspace
    CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, &
                           vhartree = v_hartree_gspace%pw)
    CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
    CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

    ! Add v_NSC_xc + v_NSC_hartree
    DO ispin = 1,nspins
      CALL pw_scale(v_rspace_nsc(ispin)%pw, v_rspace_nsc(ispin)%pw%pw_grid%dvol)
      CALL pw_axpy(v_hartree_rspace%pw,v_rspace_nsc(ispin)%pw)
    END DO
    
    DO ispin = 1 , nspins
       ! original task_list and ORB Basis
       CALL integrate_v_rspace(v_rspace = v_rspace_nsc(ispin), &
                               pmat = old_rho_ao(ispin, 1), &
                               hmat = nsc_matrix(ispin, 1), &
                               qs_env = qs_env, &
                               calculate_forces = calculate_forces, & 
                               basis_type = "ORB")
       
!       CALL integrate_v_rspace(v_rspace = v_rspace_nsc(ispin), &
!                               pmat = old_rho_ao(ispin, 1), &
!                               hmat = nsc_matrix(ispin, 1), &
!                               qs_env = qs_env, &
!                               calculate_forces = calculate_forces, & 
!                               !basis_type = "ORB")
!                               basis_type = "HARRIS", &
!                               task_list_external = ec_env%task_list)
    END DO

    DO ispin = 1,nspins
       CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
       CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_nsc(ispin)%pw)
    END DO

    DEALLOCATE(v_xc, v_rspace_nsc)

    ! release pw grids and matrices
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)

    CALL timestop(handle)

    IF(para_env%mepos==0) THEN
       WRITE(*,*) "ec_nsc_matrix - END"
    END IF

  END SUBROUTINE ec_nsc_matrix

! **************************************************************************************************
!> \brief  
!> \param qs_env ...
!> \param ec_env ...
!>        
!>        
!>       
!> \date    08.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_g2_matrix(qs_env, ec_env, g2_matrix)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(energy_correction_type), POINTER              :: ec_env
    TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: g2_matrix

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_g2_matrix', &
       routineP = moduleN//':'//routineN

    TYPE(kg_environment_type), POINTER                 :: kg_env
    INTEGER                                            :: handle
    TYPE(cp_logger_type), POINTER                      :: logger
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
    ! fbelle
    TYPE(pw_env_type), POINTER                         :: pw_env
    TYPE(qs_rho_type), POINTER                         :: rho!, old_rho
    TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                          v_hartree_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, old_rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, old_rho_g
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao, old_rho_ao
    TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
    TYPE(pw_poisson_type), POINTER                     :: poisson_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    INTEGER                                            :: ispin, nspins
    LOGICAL                                            :: calculate_forces 
    TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
    INTEGER                                            :: natom
    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
    ! fbelle
    TYPE(xc_derivative_set_type), POINTER              :: deriv_set
    TYPE(xc_rho_set_type), POINTER                     :: rho_set
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc, v_rspace_nsc
    TYPE(section_vals_type), POINTER                   :: input, xc_section

    logger => cp_get_default_logger()
    
    CALL timeset(routineN, handle)
    
    CALL get_qs_env(qs_env, &
                    kg_env = kg_env, &
                    dft_control = dft_control, &
                    rho = rho, &
                    force = force, &
                    input = input, &
                    para_env = para_env)

    IF(para_env%mepos==0) THEN
       WRITE(*,*) "ec_g2_matrix - BEGINNING"
    END IF
    
    calculate_forces = .FALSE.
    
    nspins = dft_control%nspins
     
    ! get atomic mappng
    CALL get_qs_env(qs_env, & 
                    atomic_kind_set = atomic_kind_set, &
                    natom = natom)
    ALLOCATE (atom_of_kind(natom), kind_of(natom))
    CALL get_atomic_kind_set(atomic_kind_set, &
                             atom_of_kind = atom_of_kind, &
                             kind_of = kind_of)

!!---------------------------------------------------------------------------------------
    ! *** Calculation of forces due to non-self consistency of energy correction *** ! 
    ! First, setup pw pools
!!---------------------------------------------------------------------------------------

    NULLIFY (pw_env)
    CALL get_qs_env(qs_env, pw_env=pw_env)
    CPASSERT(ASSOCIATED(pw_env))

    NULLIFY (auxbas_pw_pool, poisson_env, pw_pools)
    ! gets the tmp grids
    CALL pw_env_get(pw_env = pw_env, &
                    auxbas_pw_pool = auxbas_pw_pool, &
                    pw_pools = pw_pools, &
                    poisson_env = poisson_env)

    ! Calculate the NSC Hartree potential
    CALL pw_pool_create_pw(pool = auxbas_pw_pool, &
                           pw = v_hartree_gspace%pw, &
                           use_data = COMPLEXDATA1D, &
                           in_space = RECIPROCALSPACE)
    CALL pw_pool_create_pw(pool = auxbas_pw_pool, &
                           pw = rho_tot_gspace%pw, &
                           use_data = COMPLEXDATA1D, &
                           in_space = RECIPROCALSPACE)

    CALL pw_pool_create_pw(pool = auxbas_pw_pool, &
                           pw = v_hartree_rspace%pw, &
                           use_data = REALDATA3D, &
                           in_space = REALSPACE)

    ! get rho_out - after Harris energy correction
    CALL qs_rho_get(rho, & 
                    rho_r = rho_r, &
                    rho_g = rho_g, &
                    rho_ao_kp = rho_ao)

    ! get rho_in - before Harris energy correction
    CALL qs_rho_get(ec_env%old_rho, &
                    rho_r = old_rho_r, &
                    rho_g = old_rho_g, &
                    rho_ao_kp = old_rho_ao )

!!---------------------------------------------------------------------------------------
    ! First, nu'_xc[rho_in]
!!---------------------------------------------------------------------------------------

    ! ** v_xc ** !
    ALLOCATE(v_xc(nspins), v_rspace_nsc(nspins))
    DO ispin = 1, nspins
       NULLIFY(v_xc(nspins)%pw, v_rspace_nsc(nspins)%pw)
       CALL pw_pool_create_pw(pool = auxbas_pw_pool, & 
                              pw = v_xc(ispin)%pw, &
                              use_data = REALDATA3D, &
                              in_space = REALSPACE)
       CALL pw_pool_create_pw(pool = auxbas_pw_pool, & 
                              pw = v_rspace_nsc(ispin)%pw, &
                              use_data = REALDATA3D, &
                              in_space = REALSPACE)
       CALL pw_zero(v_xc(ispin)%pw)
       CALL pw_zero(v_rspace_nsc(nspins)%pw)
    END DO

    NULLIFY(xc_section)

    NULLIFY (rho_set, deriv_set)
    ! main ingredient is xc_rho_set_and_dset_create
    CALL xc_prep_2nd_deriv(deriv_set = deriv_set, &        ! containing potentials
                           rho_set = rho_set, &            ! density at which deriv are calced
                           rho_r = old_rho_r, &            ! where derivative evaluated 
                           pw_pool = auxbas_pw_pool, &     ! OK
                           xc_section = ec_env%xc_section)        ! OK 

    CALL xc_calc_2nd_deriv(v_xc = v_xc, & 
                           deriv_set = deriv_set, &
                           rho_set = rho_set, &
                           rho1_set = rho_set, &
                           pw_pool = auxbas_pw_pool, &
                           xc_section = ec_env%xc_section, & 
                           gapw = .FALSE.)

    DO ispin = 1,nspins
      ! \nu'(rho_in)*(rho_out - rho_in)
      v_rspace_nsc(ispin)%pw%cr3d = old_rho_r(ispin)%pw%cr3d &
                                    * v_xc(ispin)%pw%cr3d
    END DO

    CALL xc_dset_release(deriv_set)
    CALL xc_rho_set_release(rho_set = rho_set, pw_pool = auxbas_pw_pool)

    ! + v_hartree[rho_in]
    CALL pw_zero(rho_tot_gspace%pw)
    DO ispin = 1,nspins
      CALL pw_axpy(old_rho_g(ispin)%pw, rho_tot_gspace%pw)
    END DO

    ! get Hartree potential from rho_tot_gspace
    CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, &
                           vhartree = v_hartree_gspace%pw)
    CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
    CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

    ! Add v_NSC_xc + v_NSC_hartree
    DO ispin = 1,nspins
      CALL pw_scale(v_rspace_nsc(ispin)%pw, v_rspace_nsc(ispin)%pw%pw_grid%dvol)
      CALL pw_axpy(v_hartree_rspace%pw,v_rspace_nsc(ispin)%pw)
    END DO
    
    DO ispin = 1 , nspins
       ! input density and basis - use sparsity
       CALL integrate_v_rspace(v_rspace = v_rspace_nsc(ispin), &
                               pmat = old_rho_ao(ispin, 1), &
                               hmat = g2_matrix(ispin, 1), &
                               qs_env = qs_env, &
                               calculate_forces = calculate_forces, & 
                               basis_type = "ORB")
    END DO

    DO ispin = 1,nspins
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_nsc(ispin)%pw)
    END DO

    DEALLOCATE(v_xc, v_rspace_nsc)

    ! release pw grids and matrices
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)

    CALL timestop(handle)

    IF(para_env%mepos==0) THEN
       WRITE(*,*) "ec_g2_matrix - END"
    END IF

  END SUBROUTINE ec_g2_matrix

END MODULE ec_response_solver
