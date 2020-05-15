!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Harris functional force driver routine
!>
!>
!> \date 09.2019
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_ao_res_solver
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_transpose
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_write_formatted
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_methods,                   ONLY: dbcsr_get_matrix_type 
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_add_on_diag,&
                                              dbcsr_checksum,&
                                              dbcsr_complete_redistribute,&
                                              dbcsr_copy,&
                                              dbcsr_copy_into_existing,&
                                              dbcsr_create,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_distribution_type,&
                                              dbcsr_dot,&
                                              dbcsr_filter,&
                                              dbcsr_frobenius_norm,&
                                              dbcsr_multiply,&
                                              dbcsr_norm,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_print,&
                                              dbcsr_scale,&
                                              dbcsr_release,&
                                              dbcsr_transposed,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_antisymmetric
   USE ec_energy_corrections,           ONLY: ec_build_core_hamiltonian
   USE ec_env_types,                    ONLY: ec_environment_type
   USE ec_methods,                      ONLY: calculate_kernel, &
                                              ec_projvirt_matrix,&
                                              ec_w2_matrix
   USE ec_response_solver,              ONLY: ec_z_matrix
   USE input_constants,                 ONLY: ls_s_sqrt_proot,&
                                              ls_s_sqrt_ns 
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_types,                 ONLY: linres_control_release,&
                                              linres_control_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_linres_methods,               ONLY: linres_solver
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE dm_ls_scf_curvy,                 ONLY: compute_cg_matrices,& 
                                              commutator_symm,&
                                              transform_matrix_orth
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz,& 
                                              matrix_sqrt_proot
   USE qs_rho_methods,                  ONLY: duplicate_rho_type,&
                                              qs_rho_update_rho,&
                                              qs_rho_rebuild
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type,&
                                              qs_rho_clear
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
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE machine,                         ONLY: m_flush
   USE mathlib,                         ONLY: invmat
   ! 
   USE dm_ls_scf_methods,               ONLY: apply_matrix_preconditioner
      
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_ao_res_solver'

! *** Public subroutines ***

   PUBLIC :: ec_res_ao

CONTAINS

! **************************************************************************************************
!> XXXXXX
!> COPY PASTE from iterate_matrix.F because of type incompatibilities
!> XXXXXX
!> \brief McWeeny purification of a matrix in the non-orthonormal basis
!> \param matrix_p Matrix to purify (needs to be almost idempotent already)
!> \param matrix_s Overlap-Matrix
!> \param threshold Threshold used as filter_eps and convergence criteria
!> \param max_steps Max number of iterations
!> \par History
!>       2013.01 created [Florian Schiffmann]
!>       2014.07 slightly refactored [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE mcweeny(matrix_p, matrix_s, threshold, max_steps)
      TYPE(dbcsr_type)                                   :: matrix_p
      TYPE(dbcsr_type)                                   :: matrix_s
      REAL(KIND=dp)                                      :: threshold
      INTEGER                                            :: max_steps

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mcweeny', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, unit_nr
      REAL(KIND=dp)                                      :: frob_norm, trace
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_ps, matrix_psp, matrix_test

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL dbcsr_create(matrix_ps, template=matrix_p, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_psp, template=matrix_p, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_test, template=matrix_p, matrix_type=dbcsr_type_no_symmetry)

         DO i = 1, max_steps
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p, matrix_s, &
                                0.0_dp, matrix_ps, filter_eps=threshold)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_ps, matrix_p, &
                                0.0_dp, matrix_psp, filter_eps=threshold)
            IF (i == 1) CALL dbcsr_dot(matrix_ps, matrix_ps, trace)

            ! test convergence
            CALL dbcsr_copy(matrix_test, matrix_psp)
            CALL dbcsr_add(matrix_test, matrix_p, 1.0_dp, -1.0_dp)
            frob_norm = dbcsr_frobenius_norm(matrix_test) ! test = PSP - P
            IF (unit_nr > 0) WRITE (unit_nr, '(t3,a,2f16.8)') "McWeeny: Deviation of idempotency", frob_norm
            IF (unit_nr > 0) CALL m_flush(unit_nr)

            ! construct new P
            CALL dbcsr_copy(matrix_p, matrix_psp)
            CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_ps, matrix_psp, &
                                3.0_dp, matrix_p, filter_eps=threshold)

            ! frob_norm < SQRT(trace*threshold)
            IF (frob_norm*frob_norm < trace*threshold) EXIT
         END DO

      CALL dbcsr_release(matrix_ps)
      CALL dbcsr_release(matrix_psp)
      CALL dbcsr_release(matrix_test)
      CALL timestop(handle)

   END SUBROUTINE mcweeny

! ***************************************************************************************************
!> \brief  build preconditioner
!>         Ax = [F, S]_B
!>          
!>          
!>          
!> \param qs_env ...
!> \date    03.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE build_preconditioner(qs_env, matrix_ks, matrix_ps, matrix_cg, matrix_m, eps_filter)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_p_type), DIMENSION(:,:), INTENT(IN), &
                                        POINTER        :: matrix_cg, matrix_ks, matrix_ps
    TYPE(dbcsr_p_type), DIMENSION(:,:), INTENT(OUT), &
                                        POINTER        :: matrix_m
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'build_preconditioner', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle, ispin, nspins
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_b, matrix_bs, matrix_s, matrix_tmp
    TYPE(dft_control_type), POINTER                    :: dft_control

!---------------------------------------------------------------------------------------
!  Initial general setup
!---------------------------------------------------------------------------------------

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    matrix_s_kp = matrix_s,&
                    para_env=para_env)
    nspins = dft_control%nspins
    IF(para_env%mepos==0) WRITE(*,*) "build_preconditioner BEGINNING"

    NULLIFY(matrix_b)
    CALL dbcsr_allocate_matrix_set(matrix_b, nspins, 1)
    DO ispin = 1, nspins
       ALLOCATE (matrix_b(ispin, 1)%matrix)
!       CALL dbcsr_create(matrix_b(ispin, 1)%matrix, name = "[X,P] RSP DNSTY", &
!                         template = matrix_s(1, 1)%matrix)
       CALL dbcsr_copy(matrix_b(ispin, 1)%matrix, matrix_s(1, 1)%matrix)
       CALL dbcsr_set(matrix_b(ispin,1)%matrix, 0.0_dp)
    END DO
    NULLIFY(matrix_tmp)
    CALL dbcsr_allocate_matrix_set(matrix_tmp, nspins, 1)
    DO ispin = 1, nspins
       ALLOCATE (matrix_tmp(ispin, 1)%matrix)
       CALL dbcsr_create(matrix_tmp(ispin, 1)%matrix, name = "[F,B]_S ", &
                         template = matrix_ps(1, 1)%matrix)
       CALL dbcsr_set(matrix_tmp(ispin,1)%matrix, 0.0_dp)
    END DO
    NULLIFY(matrix_bs)
    CALL dbcsr_allocate_matrix_set(matrix_bs, nspins, 1)
    DO ispin = 1, nspins
       ALLOCATE (matrix_bs(ispin, 1)%matrix)
!       CALL dbcsr_create(matrix_bs(ispin, 1)%matrix, name = "B*S work matrix", &
!                         template = matrix_ps(1, 1)%matrix)
       CALL dbcsr_copy(matrix_bs(ispin, 1)%matrix, matrix_ks(1, 1)%matrix)
       CALL dbcsr_set(matrix_bs(ispin,1)%matrix, 0.0_dp)
    END DO

    ! Build intermediate density matrix
    ! B = [cg, P]_S = cg*S*P - P*S*cg
    CALL ec_scomm(qs_env = qs_env, &
                  matrix_a = matrix_cg, &
                  matrix_ps = matrix_ps, &
                  matrix_out = matrix_b, &
                  eps_filter = eps_filter, &
                  alpha = 1.0_dp, &
                  beta = 0.0_dp, &
                  ps_trans = .TRUE.)

    ! B * S
    DO ispin = 1, nspins
       CALL dbcsr_multiply("N","N", 1.0_dp, matrix_b(1, 1)%matrix, matrix_s(1, 1)%matrix, &
                           0.0_dp, matrix_bs(ispin,1)%matrix, filter_eps = eps_filter)
    END DO

    ! Build first part of operator
    ! Ax = [F,S]_B
    ! Ax = F*B*S - S*B*F
    CALL ec_scomm(qs_env, &
                  matrix_a = matrix_ks, &
                  matrix_ps = matrix_bs, &
                  matrix_out = matrix_m, &
                  eps_filter = eps_filter, &
                  alpha = 1.0_dp, &
                  beta = 0.0_dp, &
                  ps_trans = .FALSE.)

!    DO ispin = 1, nspins
!       CALL dbcsr_transposed(matrix_tmp(ispin,1)%matrix, matrix_m(ispin,1)%matrix)
!       CALL dbcsr_copy(matrix_m(ispin,1)%matrix, matrix_tmp(ispin,1)%matrix)
!    END DO

!    ! proj(M) = [E, P]_S = ESP - PSE
!    CALL ec_scomm(qs_env, &
!                  matrix_a = matrix_tmp, &
!                  matrix_b = matrix_s, &
!                  matrix_c = matrix_p, &
!                  matrix_out = matrix_m, &
!                  eps_filter = eps_filter, &
!                  alpha = 1.0_dp, &
!                  beta = 0.0_dp)

!    ! Apply projector on m
!    CALL matrix_proj_t(qs_env = qs_env, &
!                       matrix_p = matrix_p, &
!                       matrix_s = matrix_s, &
!                       matrix_out = matrix_m)

    CALL dbcsr_deallocate_matrix_set(matrix_b)
    CALL dbcsr_deallocate_matrix_set(matrix_bs)
    CALL dbcsr_deallocate_matrix_set(matrix_tmp)

    CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief: Solve AX + B = 0  
!> \param qs_env ...
!>  \param ec_env Harris energy correction environment         
!>  \param p_env Storage of response density
!>          
!> \date    01.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_res_ao(qs_env, ec_env, p_env)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(ec_environment_type), POINTER                 :: ec_env
    TYPE(qs_p_env_type), POINTER                       :: p_env

    CHARACTER(len=*), PARAMETER :: routineN = 'ec_res_ao', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: i, handle, ispin, ncyc, nspins
    TYPE(cp_logger_type), POINTER                      :: logger
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    TYPE(qs_rho_type), POINTER                         :: rho
    TYPE(neighbor_list_set_p_type), DIMENSION(:), &
       POINTER                                         :: sab_orb
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_ks, matrix_m, matrix_p, matrix_s, &
                                                          rho_ao, matrix_nsc, ksmat
    LOGICAL                                            :: converged
    REAL(KIND=dp)                                      :: eps_filter, norm_res
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_Ax, matrix_cg, matrix_cg_z, &
                                                          matrix_res, matrix_rhs, matrix_z, matrix_z0  
    TYPE(dbcsr_type)                                   :: matrix_tmp
    REAL(KIND=dp)                                      :: t1, t2
    INTEGER                                            :: unit_nr
    REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha, beta, new_norm, norm_cA, norm_rr
    INTEGER                                            :: iounit
    REAL(KIND=dp)                                      :: conv_val, old_conv, shift 
    LOGICAL                                            :: do_precond
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_ps
    REAL(KIND=dp), DIMENSION(:), POINTER               :: norm_zmo, norm_zao 

    CALL timeset(routineN, handle)
    NULLIFY (logger)
    logger => cp_get_default_logger()
    IF (logger%para_env%ionode) THEN
       unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
    ELSE
       unit_nr = -1
    ENDIF
    t1 = m_walltime()

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    matrix_ks_kp = ksmat, &
                    matrix_s_kp = matrix_s, &
                    para_env = para_env, &
                    sab_orb = sab_orb, &
                    rho = rho)
    IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao BEGINNING"
    nspins = dft_control%nspins
    eps_filter = 1.0E-30_dp
    iounit = 6

    ! Get the ground state density
    CALL qs_rho_get(rho, rho_ao_kp = rho_ao)

    ALLOCATE(alpha(nspins),beta(nspins), new_norm(nspins), norm_cA(nspins), norm_rr(nspins))
    ALLOCATE(norm_zmo(nspins), norm_zao(nspins))

    ! local matrix P, KS, and NSC
    ! Get local desymmetrized versions of 
    ! the ground-state density matrix (matrix_p)
    ! the Kohn-Sham matrix (matrix_ks)
    ! the RHS property gradient i.e. Harris functional (matrix_nsc)
    NULLIFY(matrix_p, matrix_ks, matrix_nsc)
    CALL dbcsr_allocate_matrix_set(matrix_p, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_ks, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_nsc, nspins, 1)
    DO ispin = 1, nspins
       ALLOCATE (matrix_p(ispin, 1)%matrix)
       ALLOCATE (matrix_ks(ispin, 1)%matrix)
       ALLOCATE (matrix_nsc(ispin, 1)%matrix)
       ! XXX copy afterwards also changes name, unecessary here
       CALL dbcsr_create(matrix_p(ispin, 1)%matrix, name = "P_IN", &
                         template = ksmat(1, 1)%matrix)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_p(ispin, 1)%matrix, sab_orb)
       CALL dbcsr_create(matrix_ks(ispin, 1)%matrix, name = "KS_IN", &
                         template = ksmat(1, 1)%matrix)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks(ispin, 1)%matrix, sab_orb)
       CALL dbcsr_create(matrix_nsc(ispin, 1)%matrix, name = "NSC", &
                         template = ksmat(1, 1)%matrix)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_nsc(ispin, 1)%matrix, sab_orb)

       ! XXX not necessary, taken care of by dbcsr_multiply
       ! XXX in that case, no copy of KS matric required
       ! XXX or of RHS-matrix either
       !CALL dbcsr_desymmetrize(rho_ao(ispin,1)%matrix, matrix_p(ispin,1)%matrix)
       !CALL dbcsr_desymmetrize(ksmat(ispin,1)%matrix, matrix_ks(ispin,1)%matrix)

       CALL dbcsr_copy(matrix_p(ispin,1)%matrix, rho_ao(ispin,1)%matrix)
       CALL dbcsr_copy(matrix_ks(ispin,1)%matrix, ksmat(ispin,1)%matrix)
       CALL dbcsr_copy(matrix_nsc(ispin, 1)%matrix, ec_env%matrix_hz(ispin)%matrix, &
                       name="NSC-matrix")
    END DO

    ! Scale matrix_p by factor 1/2 in closed-shell
    IF (nspins == 1) CALL dbcsr_scale(matrix_p(1,1)%matrix, 0.5_dp)

    ! Ensure idempoteny of ground-state density
    DO ispin = 1, nspins
       CALL mcweeny(matrix_p(ispin,1)%matrix, matrix_s(ispin,1)%matrix, eps_filter, 3)
    END DO

    ! Seting up work matrices : 
    ! matrix_Ax : The Hessian linear transformation trial matrix
    ! matrix_cg : Trial matrix
    ! matrix_cg_z : Approximate solution
    ! matrix_m : Preconditioner
    ! matrix_ps : Precomputed product P*S used in S-commutator 
    ! matrix_res : Residual
    ! matrix_z : Response density
    ! matrix_z0 : Preconditioned trial matrix

    CALL dbcsr_create(matrix_tmp, template=matrix_p(1,1)%matrix, matrix_type=dbcsr_type_no_symmetry)

    NULLIFY(matrix_Ax, matrix_cg, matrix_cg_z, matrix_m, matrix_ps, matrix_res, &
            matrix_rhs, matrix_z, matrix_z0)
    CALL dbcsr_allocate_matrix_set(matrix_Ax, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_cg, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_cg_z, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_m, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_ps, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_res, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_rhs, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_z, nspins, 1)
    CALL dbcsr_allocate_matrix_set(matrix_z0, nspins, 1)
    DO ispin = 1, nspins
       ALLOCATE (matrix_Ax(ispin, 1)%matrix)
       ALLOCATE (matrix_cg(ispin, 1)%matrix)
       ALLOCATE (matrix_cg_z(ispin, 1)%matrix)
       ALLOCATE (matrix_m(ispin, 1)%matrix)
       ALLOCATE (matrix_ps(ispin, 1)%matrix)
       ALLOCATE (matrix_res(ispin, 1)%matrix)
       ALLOCATE (matrix_rhs(ispin, 1)%matrix)
       ALLOCATE (matrix_z(ispin, 1)%matrix)
       ALLOCATE (matrix_z0(ispin, 1)%matrix)
       CALL dbcsr_create(matrix_Ax(ispin, 1)%matrix, name = "linop MATRIX", &
                         template = matrix_s(1, 1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_cg(ispin, 1)%matrix, name = "TRIAL MATRIX", &
                         template = matrix_s(1, 1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_cg_z(ispin, 1)%matrix, name = "MATRIX CG-Z", &
                         template = matrix_s(1, 1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_m(ispin, 1)%matrix, name = "Preconditioner", &
                         template = matrix_s(1, 1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_ps(ispin, 1)%matrix, name = "MATRIX PS", &
                          template = matrix_s(1, 1)%matrix, &
                          matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_res(ispin, 1)%matrix, name = "RESIDUAL", &
                         template = matrix_s(1, 1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_rhs(ispin, 1)%matrix, name = "RHS", &
                         template = matrix_s(1, 1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_z(ispin, 1)%matrix, name = "Z-Matrix", &
                         template = matrix_s(1, 1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_z0(ispin, 1)%matrix, name = "p after precondi-Matrix", &
                         template = matrix_s(1, 1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)

       CALL cp_dbcsr_alloc_block_from_nbl(matrix_Ax(ispin, 1)%matrix, sab_orb)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_cg(ispin, 1)%matrix, sab_orb)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_cg_z(ispin, 1)%matrix, sab_orb)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_m(ispin, 1)%matrix, sab_orb)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_res(ispin, 1)%matrix, sab_orb)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_z0(ispin, 1)%matrix, sab_orb)

       CALL dbcsr_set(matrix_Ax(ispin,1)%matrix, 0.0_dp)
       CALL dbcsr_set(matrix_cg(ispin,1)%matrix, 0.0_dp)
       CALL dbcsr_set(matrix_cg_z(ispin,1)%matrix, 0.0_dp)
       CALL dbcsr_set(matrix_m(ispin,1)%matrix, 0.0_dp)
       CALL dbcsr_set(matrix_ps(ispin,1)%matrix, 0.0_dp)
       CALL dbcsr_set(matrix_res(ispin,1)%matrix, 0.0_dp)
       CALL dbcsr_set(matrix_z(ispin,1)%matrix, 0.0_dp)
       CALL dbcsr_set(matrix_z0(ispin,1)%matrix, 0.0_dp)
    END DO

    ! matrix_ps
    DO ispin = 1, nspins
       CALL dbcsr_multiply("N","N", 1.0_dp, matrix_p(1, 1)%matrix, matrix_s(1, 1)%matrix, &
                           0.0_dp, matrix_ps(ispin,1)%matrix, filter_eps = eps_filter)
    END DO

    !----------------------------------------
    ! Get righ-hand-side (RHS) operator
    !----------------------------------------

    ! E^[1]_Harris = G(\delta P)P*S - S*P*G(\delta P)
    ! ABC-CBA
    ! matrix_rhs antiysymmetric
    CALL ec_scomm(qs_env = qs_env, &
                  matrix_a = matrix_nsc , &
                  matrix_ps = matrix_ps, &
                  matrix_out = matrix_rhs, &
                  eps_filter = eps_filter, &
                  alpha = -2.0_dp, &
                  beta = 0.0_dp, &
                  ps_trans = .FALSE.)

    ! r_0 = b
    CALL dbcsr_copy(matrix_res(ispin,1)%matrix, matrix_rhs(ispin,1)%matrix)

    ! Frist trial matrix is RHS 
!    CALL dbcsr_copy(matrix_cg_z(ispin,1)%matrix, matrix_res(ispin,1)%matrix)
    ! Set first trial matrix to 0
    DO ispin = 1, nspins
       CALL dbcsr_set(matrix_cg_z(ispin,1)%matrix, 0.0_dp)
    END DO

    ! Apply projector on cg_z
!    CALL matrix_proj(qs_env, matrix_ps, matrix_cg_z, eps_filter)

    ! Ax0
    ! Build Hessian transformation Ax0 on first trial matrix
    CALL build_hessian_op(qs_env = qs_env, &
                          matrix_ks = matrix_ks, &   ! IN
                          matrix_p  = matrix_p,  &   ! IN
                          matrix_ps = matrix_ps, &   ! IN
                          matrix_cg = matrix_cg_z, & ! IN trial matrix
                          matrix_Ax = matrix_Ax, &   ! OUT Hessian
                          eps_filter = eps_filter, &
                          iter = 1)

    ! Preconditioner 
    ! First term of Hessian operator (uncoupled term)
    CALL build_preconditioner(qs_env = qs_env, &
                              matrix_ks = matrix_ks, &   ! IN
                              matrix_ps = matrix_ps, &   ! IN
                              matrix_cg = matrix_cg_z,&  ! IN
                              matrix_m = matrix_m, &     ! OUT precond
                              eps_filter = eps_filter)
                              
    ! r_0 = b - Ax0
    DO ispin = 1, nspins
       CALL dbcsr_add(matrix_res(ispin,1)%matrix, matrix_Ax(ispin,1)%matrix, 1.0_dp, -1.0_dp)
    END DO

    ! Apply projector on r
!    CALL matrix_proj_t(qs_env, matrix_ps, matrix_res, eps_filter)

    ! XXX Preconditioner not properly implemented
!    do_precond = .TRUE.
    do_precond = .FALSE.
    ! Preconditioner
    DO ispin = 1, nspins
       IF (do_precond) THEN
          ! z_0 = M * r_0
          ! M ASYM, r_0 ASYM
          ! results in z_0 NO-SYM 
          IF(para_env%mepos==0) WRITE(*,*) "ec_ao_res - APPLY PRECONDITIONER z_0 = M * r_0"

          ! M*matrix_res = Asym*asym = no symmetry
          CALL dbcsr_multiply("N","N", 1.0_dp, matrix_m(ispin,1)%matrix, matrix_res(ispin,1)%matrix,&
                              0.0_dp, matrix_z0(ispin,1)%matrix, filter_eps = eps_filter)

!          WRITE(*,*) "ec_res_ao - M*r_0"
!          CALL dbcsr_print(matrix_z0(1,1)%matrix)
!          CALL dbcsr_print(matrix_z0(1,1)%matrix, matlab_format = .TRUE.)

          ! [M, res] = MR - RM, result asym, but wrong z matrix
!          CALL commutator_symm(matrix_m(ispin,1)%matrix, matrix_res(ispin, 1)%matrix, &
!                               matrix_z0(ispin,1)%matrix, eps_filter, 1.0_dp)

!          IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - M*r_0 - after antisymmetrization"
!          CALL dbcsr_print(matrix_z0(1,1)%matrix, matlab_format = .TRUE.)
       ELSE
          ! z_0 = r_0
          CALL dbcsr_copy(matrix_z0(ispin,1)%matrix, matrix_res(ispin,1)%matrix)
       END IF
    END DO

    ! Init Convergence criteria 
    norm_res = 0.0_dp

    DO ispin = 1, nspins

       ! cg = p_0 = z_0
       CALL dbcsr_copy(matrix_cg(ispin,1)%matrix, matrix_z0(ispin,1)%matrix)

       ! Tr(r_0 * z_0)
       CALL dbcsr_dot(matrix_res(ispin,1)%matrix, matrix_cg(ispin, 1)%matrix, norm_rr(ispin))
       IF(para_env%mepos==0) WRITE(*,*) "Tr[r_0 * z_0]", norm_rr(1)

       IF (norm_rr(ispin) .LT. 0.0_dp) CPABORT("Tr[r_0 * z_0] < 0")
       norm_res = MAX(norm_res, ABS(norm_rr(ispin)))
       IF(para_env%mepos==0) WRITE(*,*) "norm_res", norm_res
    END DO

    alpha(:) = 0.0_dp

    ! header
    IF (iounit > 0) THEN
       IF(para_env%mepos==0) WRITE (iounit, "(/,T3,A,T16,A,T25,A,T38,A,/,T3,A)") &
          "Iteration", "Stepsize", "Convergence", "Time", &
          REPEAT("-", 80)
    ENDIF

    ! Notation in other CG-routines 
    ! Florian  ||| Marcella
    ! norm rr  ||| tr_rz0
    ! norm_res ||| norm_res
    ! norm_cA  ||| tr_pAp
    ! new_norm ||| tr_rz1

    ! Start iteration
    ncyc = 200
    DO i = 1, ncyc
       IF(para_env%mepos==0) WRITE(*,*) "START ITERATION", i

       ! Convergence criteria
       converged=.FALSE.
       ! Max number of iteration reached
       IF (i == ncyc) THEN
          IF (iounit > 0) THEN
             IF(para_env%mepos==0) WRITE (iounit, "(/,T2,A/)") &
                "The linear solver didnt converge! Maximum number of iterations reached."
             CALL m_flush(iounit)
          ENDIF
          EXIT 
       ENDIF
    
       ! 1.0E-6_dp in MO-linear response solver
       IF (norm_res .LT. 1.0E-10_dp) THEN
          converged = .TRUE.
       ENDIF

       t2 = m_walltime()
       IF (iounit > 0) THEN
          IF(para_env%mepos==0) WRITE(*,*) i, MAXVAL(alpha), norm_res, t2-t1
          CALL m_flush(iounit)
       ENDIF

       IF (converged) THEN
          IF (iounit > 0) THEN
             IF(para_env%mepos==0) WRITE (iounit, "(/,T2,A,I4,A,/)") "The linear solver converged in ", i, " iterations."
             CALL m_flush(iounit)
          ENDIF
          EXIT
       ENDIF

       ! Hessian linear transformation on trial matrix p_j (matrix_cg)
       ! Ap_j = [F,S]_B + [G(B),S]_P and B = [p_j,P]_S
       CALL build_hessian_op(qs_env = qs_env, &
                             matrix_ks = matrix_ks, & ! IN
                             matrix_p  = matrix_p,  & ! IN
                             matrix_ps = matrix_ps, & ! IN
                             matrix_cg = matrix_cg, & ! IN trial matrix
                             matrix_Ax = matrix_Ax, & ! OUT (Ax(cg))
                             eps_filter = eps_filter, &
                             iter = i)

       ! Apply projector on Ax
!       CALL matrix_proj_t(qs_env, matrix_ps, matrix_Ax, eps_filter)

       ! Preconditioner 
       CALL build_preconditioner(qs_env = qs_env, &
                                 matrix_ks = matrix_ks, &  ! IN
                                 matrix_ps = matrix_ps, &  ! IN
                                 matrix_cg = matrix_cg,&   ! IN
                                 matrix_m = matrix_m, &
                                 eps_filter = eps_filter)  ! OUT precond

       DO ispin = 1, nspins

          CALL dbcsr_filter(matrix_Ax(ispin,1)%matrix, eps_filter)

          ! tr(Ap_j * p_j)
          CALL dbcsr_dot(matrix_cg(ispin,1)%matrix, matrix_Ax(ispin,1)%matrix, norm_cA(ispin))
          IF(para_env%mepos==0) WRITE(*,*) "norm_cA", norm_cA(ispin)

          IF(norm_cA(ispin) == 0.0_dp) norm_cA(ispin) = 1.0E-10_dp

          ! XXX Needs to be extended
          ! If trace negative, recalculate Hessian for unpreconditioned trial matrix
          IF (norm_cA(ispin) .LT. 0.0_dp) THEN
             CPABORT("tr(Ap_j*p_j) < 0")
          END IF

          ! Determine step-size
          ! alpha = norm_rr/norm_cA
          IF (norm_cA(ispin) .LT. 1.0E-10_dp) THEN
             alpha(ispin) = 1.0_dp
          ELSE
             alpha(ispin) = norm_rr(ispin)/norm_cA(ispin)
          END IF
          IF(para_env%mepos==0) WRITE(*,*) "alpha", alpha(ispin)

          ! x_j+1 = x_j + alpha*p_j
          ! save approximate solution
          CALL dbcsr_add(matrix_cg_z(ispin,1)%matrix, matrix_cg(ispin,1)%matrix, 1.0_dp, alpha(ispin))

          ! r_j+1 = r_j - alpha * Ap_j 
          ! update residual 
          CALL dbcsr_add(matrix_res(ispin,1)%matrix, matrix_Ax(ispin,1)%matrix, 1.0_dp, -alpha(ispin))

          ! XXX Preconditioner not properly implemented
          !do_precond = .TRUE.
          do_precond = .FALSE.
          IF (do_precond) THEN
             ! z0 <- M*res

             ! M*matrix_res = Asym*asym = no symmetry
             CALL dbcsr_multiply("N","N", 1.0_dp, matrix_m(ispin,1)%matrix, &
                                 matrix_res(ispin,1)%matrix,&
                                 1.0_dp, matrix_z0(ispin,1)%matrix, filter_eps = eps_filter)

             ! get antisym
             CALL dbcsr_transposed(matrix_tmp, matrix_z0(ispin,1)%matrix)
             CALL dbcsr_add(matrix_z0(ispin,1)%matrix, matrix_tmp, 0.5_dp, -0.5_dp)
             ! [M, res] = MR - RM, result asym, but wrong z matrix
!             CALL commutator_symm(matrix_m(ispin,1)%matrix, matrix_res(ispin, 1)%matrix, &
!                                  matrix_z0(ispin,1)%matrix, eps_filter, 1.0_dp)

          ELSE
             CALL dbcsr_copy(matrix_z0(ispin,1)%matrix, matrix_res(ispin,1)%matrix)
          END IF
       END DO

       norm_res = 0.0_dp

       DO ispin = 1, nspins
          ! Tr[r_j+1*r_j+1]
          CALL dbcsr_dot(matrix_res(ispin,1)%matrix, matrix_z0(ispin,1)%matrix, new_norm(ispin))
          IF (new_norm(ispin) .LT. 0.0_dp) CPABORT("tr(r_j+1*z_j+1) < 0")
          norm_res = MAX(norm_res, new_norm(ispin))
          IF (para_env%mepos==0) WRITE(*,*) "norm_res", norm_res

          IF (norm_rr(ispin) .LT. 1.0E-30_dp .OR. new_norm(ispin) .LT. 1.0E-30_dp) THEN
             beta(ispin) = 0.0_dp
             converged = .TRUE.
          ELSE
             beta(ispin) = new_norm(ispin)/norm_rr(ispin)
          END IF
          IF (para_env%mepos==0) Write(*,*) "beta", beta(1)

          ! update new trial matrix
          ! p_j+1 = z_j+1 + beta*p_j
          CALL dbcsr_add(matrix_cg(ispin,1)%matrix, matrix_z0(ispin,1)%matrix, beta(ispin), 1.0_dp)
          CALL dbcsr_filter(matrix_cg(ispin,1)%matrix, eps_filter)

          norm_rr(ispin) = new_norm(ispin)

          ! XXX necessary if tr[pAp] negative and needs to be recalculated without preconditioner
          ! tr_rz00(ispin) = norm_rr(ispin)
       END DO
    END DO ! ncyc

    ! Apply projector on cg_z
!    CALL matrix_proj(qs_env, matrix_ps, matrix_cg_z, eps_filter = eps_filter)
!    CALL matrix_proj_general(qs_env, matrix_ps, matrix_cg_z)

    ! Build final Z-Matrix 
    ! Z = [cg_z, P]_S = cg_z*S*P - P*S*cg_z
    !                 = cg_z*(P*S)^T + (cg_z*(PS)^T)^T actual calculation
    CALL ec_scomm(qs_env = qs_env, &
                  matrix_a = matrix_cg_z , &
                  matrix_ps = matrix_ps, &
                  matrix_out = matrix_z, &
                  eps_filter = eps_filter, &
                  alpha = 1.0_dp, &
                  beta = 0.0_dp,&
                  ps_trans = .TRUE.)

    ! Test Z matrix MO vs AO

    ! Z-Matrix from linres_solver (MO)
    norm_zmo(1) = dbcsr_frobenius_norm(p_env%p1(1, 1)%matrix)
    norm_zao(1) = dbcsr_frobenius_norm(matrix_z(1, 1)%matrix)
    
    IF (para_env%mepos==0) WRITE(*,*) "|| Z (MO) ||_F :", norm_zmo(1)
    IF (para_env%mepos==0) WRITE(*,*) "|| Z (AO) ||_F :", norm_zao(1)
    norm_res = norm_zmo(1) - norm_zao(1)
    IF (para_env%mepos==0) WRITE(*,*) "||Z(MO)||_F - ||Z(AO)||_F :", norm_res

    ! build Az
    CALL build_hessian_op(qs_env = qs_env, &
                          matrix_ks = matrix_ks, &   ! IN
                          matrix_p  = matrix_p,  &   ! IN
                          matrix_ps = matrix_ps, &   ! IN
                          matrix_cg = matrix_cg_z, & ! IN converged trial matrix 
                          matrix_Ax = matrix_Ax, &   ! OUT (A(cg_z))
                          eps_filter = eps_filter)

    ! r_z = A(z) - b = 0
    DO ispin = 1, nspins
       CALL dbcsr_add(matrix_Ax(ispin,1)%matrix, matrix_rhs(ispin,1)%matrix, 1.0_dp -1.0_dp)

       CALL dbcsr_dot(matrix_Ax(ispin,1)%matrix, matrix_Ax(ispin,1)%matrix, norm_res)
    END DO
    IF (para_env%mepos==0) WRITE(*,*) "Tr[r_z * r_z]", norm_res

    ! Z_MO - Z_AO
    CALL dbcsr_desymmetrize(ec_env%matrix_z(1, 1)%matrix, matrix_cg(1,1)%matrix)
    CALL dbcsr_add(matrix_cg(1, 1)%matrix, matrix_z(1,1)%matrix, 1.0_dp, -1.0_dp)
    norm_res = dbcsr_frobenius_norm(matrix_cg(1, 1)%matrix)
    IF (para_env%mepos==0) WRITE(*,*) "|| Z_MO - Z_AO ||_F :", norm_res

    ! export z matrix
    DO ispin = 1, nspins
       CALL dbcsr_copy(p_env%p1(ispin), matrix_z(ispin,1)%matrix, keep_sparsity = .TRUE.)
    END DO 

    ! Release matrices
    CALL dbcsr_deallocate_matrix_set(matrix_p)
    CALL dbcsr_deallocate_matrix_set(matrix_ks)
    CALL dbcsr_deallocate_matrix_set(matrix_nsc)
    CALL dbcsr_deallocate_matrix_set(matrix_ps)

    CALL dbcsr_deallocate_matrix_set(matrix_Ax)
    CALL dbcsr_deallocate_matrix_set(matrix_cg)
    CALL dbcsr_deallocate_matrix_set(matrix_cg_z)
    CALL dbcsr_deallocate_matrix_set(matrix_m)
    CALL dbcsr_deallocate_matrix_set(matrix_res)
    CALL dbcsr_deallocate_matrix_set(matrix_rhs)
    CALL dbcsr_deallocate_matrix_set(matrix_z)
    CALL dbcsr_deallocate_matrix_set(matrix_z0)

    DEALLOCATE(alpha, beta, new_norm, norm_cA, norm_rr)
    DEALLOCATE(norm_zmo, norm_zao)

    IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - END"

    CALL timestop(handle)

  END SUBROUTINE ec_res_ao

! ***************************************************************************************************
!> \brief  calculate linear transformation of Hessian matrix on a trial vector matrix matrix_cg 
!>         Ax = [F, S]_B + [G(B), S]_Pin
!>         with response density B = [cg, P]_S,
!>         Kohn-Sham (F) and density matrix (P) of the converged ground-state calculation, 
!>         overlap matrix S, kernel matrix G(B)
!>          
!> \param qs_env ...
!> \matrix_ks Kohn-Sham matrix of ground state calculation
!> \matrix_ps, precomputed product of 
!> \matrix_cg, 
!> \matrix_Ax,
!> \eps_filter, 

!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE build_hessian_op(qs_env, matrix_ks, matrix_ps, matrix_cg, matrix_Ax, &
                              eps_filter, iter)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_p_type), DIMENSION(:,:), INTENT(IN), &
                                        POINTER        :: matrix_cg, matrix_ks, matrix_ps
    TYPE(dbcsr_p_type), DIMENSION(:,:), INTENT(OUT), &
                                        POINTER        :: matrix_Ax
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'build_hessian_op', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle, ispin, nspins
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_rho_type), POINTER                         :: rho, rho_b
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_b, matrix_s, rho_ao_b
    TYPE(dft_control_type), POINTER                    :: dft_control
    REAL(KIND=dp)                                      :: chksum
    REAL(KIND=dp)                                      :: trace
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_bs
    TYPE(dbcsr_type)                                   :: matrix_tmp
    TYPE(neighbor_list_set_p_type), DIMENSION(:), &
       POINTER                                         :: sab_orb

!---------------------------------------------------------------------------------------
!  Initial general setup
!---------------------------------------------------------------------------------------

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    matrix_s_kp = matrix_s,&
                    para_env=para_env,&
                    rho = rho, &
                    sab_orb = sab_orb)
    nspins = dft_control%nspins
    IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op BEGINNING", iter

    NULLIFY(matrix_b)
    CALL dbcsr_allocate_matrix_set(matrix_b, nspins, 1)
    DO ispin = 1, nspins
       ALLOCATE (matrix_b(ispin, 1)%matrix)
       CALL dbcsr_create(matrix_b(ispin, 1)%matrix, name = "[X,P] RSP DNSTY", &
                         template = matrix_s(1, 1)%matrix, &
                         matrix_type = dbcsr_type_no_symmetry)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_b(ispin, 1)%matrix, sab_orb)
!       CALL dbcsr_copy(matrix_b(ispin,1)%matrix, matrix_ks(ispin, 1)%matrix)
       CALL dbcsr_set(matrix_b(ispin,1)%matrix, 0.0_dp)
    END DO
    NULLIFY(matrix_bs)
    CALL dbcsr_allocate_matrix_set(matrix_bs, nspins, 1)
    DO ispin = 1, nspins
       ALLOCATE (matrix_bs(ispin, 1)%matrix)
       CALL dbcsr_create(matrix_bs(ispin, 1)%matrix, name = "B*S work matrix", &
                         template = matrix_ks(1, 1)%matrix, &
                         matrix_type = dbcsr_type_no_symmetry)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_bs(ispin, 1)%matrix, sab_orb)
       CALL dbcsr_set(matrix_bs(ispin,1)%matrix, 0.0_dp)
    END DO
    CALL dbcsr_create(matrix_tmp, template=matrix_s(1,1)%matrix, matrix_type=dbcsr_type_no_symmetry)

!    WRITE(*,*) "build_hessian_op - matrix_b - before B=[cg,P]_S"
!    CALL dbcsr_print(matrix_b(1,1)%matrix)

    ! Build intermediate density matrix
    ! B = [cg, P]_S = cg*S*P - P*S*cg
    !               = cg*(P*S)^T + (cg*(PS)^T)^T is actual calculation
    IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op - calculate B", iter
    CALL ec_scomm(qs_env, &
                  matrix_a = matrix_cg, &
                  matrix_ps = matrix_ps, &
                  matrix_out = matrix_b, &
                  eps_filter = eps_filter, &
                  alpha = 1.0_dp, &
                  beta = 0.0_dp,&
                  ps_trans = .TRUE.)

   !-----------------------
   ! Test symmetrie of B
   CALL dbcsr_transposed(matrix_tmp, matrix_b(1,1)%matrix)
   ! tmp = Z^T - Z = 0 if symmetric
   CALL dbcsr_add(matrix_tmp, matrix_b(1,1)%matrix, 1.0_dp, -1.0_dp)
   WRITE(*,*) "SYMMETRY TEST OF B MATRIX"
!   CALL dbcsr_print(matrix_tmp)
!   CALL dbcsr_print(matrix_tmp, matlab_format = .TRUE.)
   trace = dbcsr_frobenius_norm(matrix_tmp)

   IF (para_env%mepos==0) WRITE(*,*) "|| B^T -  B (AO) ||_F :", trace

   !-----------------------

!    WRITE(*,*) "build_hessian_op - matrix_b - after B=[cg,P]_S"
!    CALL dbcsr_print(matrix_b(1,1)%matrix)
!    CALL dbcsr_print(matrix_b(1,1)%matrix, matlab_format = .TRUE.)

    ! matrix_ls_to_qs 
    ! CALL dbcsr_copy_into_existing(matrix_qs, matrix_ls)

!    WRITE(*,*) "build_hessian_op - matrix_bs - before s copy"
!    CALL dbcsr_print(matrix_bs(1,1)%matrix)

!    ! keep s sparsity
!    DO ispin = 1, nspins
!       CALL dbcsr_copy(matrix_bs(ispin, 1)%matrix, matrix_s(ispin,1)%matrix)
!    END DO
                
!    WRITE(*,*) "build_hessian_op - matrix_bs - after s copy"
!    CALL dbcsr_print(matrix_bs(1,1)%matrix)

!    write(*,*) "build_hessian_op - BEFORE projector"
!    CALL matrix_proj(qs_env, matrix_ps, matrix_b, eps_filter)
!    CALL matrix_proj_general(qs_env, matrix_ps, matrix_b)


!    WRITE(*,*) "build_hessian_op - matrix_b - after projector"
!    CALL dbcsr_print(matrix_b(1,1)%matrix)


    ! matrix_bs 
    ! Analogous to other S-commutators, we precompute the matrix product B*S
    ! to use same efficient routine to calculate S-commutator
    IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op - calculate BS"
    DO ispin = 1, nspins
       CALL dbcsr_multiply("N","N", 1.0_dp, matrix_b(1, 1)%matrix, matrix_s(1, 1)%matrix, &
                           0.0_dp, matrix_bs(ispin,1)%matrix, filter_eps = eps_filter)
    END DO

    ! Ax1 = [F,S]_B = FBS - SBF
    !               = F*BS - (F*BS)^T
    IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op - calculate Ax1"
    CALL ec_scomm(qs_env, &
                  matrix_a = matrix_ks, &
                  matrix_ps = matrix_bs, &
                  matrix_out = matrix_Ax, &
                  eps_filter = eps_filter, &
                  alpha = 1.0_dp, &
                  beta = 0.0_dp,&
                  ps_trans = .FALSE.)

    ! Check 2nd part
    chksum = 0.0_dp
    DO ispin = 1, nspins
       chksum = chksum+dbcsr_checksum(matrix_b(ispin,1)%matrix)
    ENDDO

    IF (para_env%mepos==0) WRITE(*,*) "chksum:", chksum

    ! skip the kernel if the DM is very small
    !IF (chksum .GT. 1.0E-14_dp) THEN
    IF (chksum .GT. 1.0E-30_dp) THEN
       !IF(para_env%mepos==0) WRITE(*,*) "continue with hessian_op2", chksum

       ! Bring B as density on grid
       NULLIFY(rho_b)
       CALL qs_rho_create(rho_b)
       CALL duplicate_rho_type(rho_input = rho, &
                               rho_output = rho_b, &
                               qs_env = qs_env)

       ! Get response density matrix
       CALL qs_rho_get(rho_b, rho_ao_kp = rho_ao_b)

       ! Copy b-density matrix to its rho environment for collocation
       ! going from non-symmetric to symmetric
       DO ispin = 1, nspins
          CALL dbcsr_copy(rho_ao_b(ispin, 1)%matrix, matrix_b(ispin,1)%matrix, keep_sparsity =.TRUE.)
       END DO

!       WRITE(*,*) "build_hessian_op - rho_ao_b"
!       CALL dbcsr_print(rho_ao_b(1,1)%matrix)

       ! update rho_r_b and rho_g_b wrt rho_ao_b
       ! Same basis as reference calculation
       CALL qs_rho_update_rho(rho_b, qs_env)

       ! Calculate kernel
       ! Ax = Ax1           + Ax2
       ! Ax = F*B*S - S*B*F + G(B)*P*S - S*P*G(B)
       !                                IN/OUT     IN
       CALL hessian_op2_nonorth(qs_env, matrix_Ax, matrix_ps, rho_b, eps_filter)

!       WRITE(*,*) "build_hessian_op - matrix_Ax after op2"
!       CALL dbcsr_print(matrix_Ax(1,1)%matrix)

       CALL qs_rho_release(rho_b)

    END IF

    CALL dbcsr_release(matrix_tmp)
    CALL dbcsr_deallocate_matrix_set(matrix_b)
    CALL dbcsr_deallocate_matrix_set(matrix_bs)

    CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief  
!>          
!>          
!>          
!>          
!> \param qs_env ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE hessian_op2_nonorth(qs_env, matrix_Ax, matrix_ps, rho_b, eps_filter)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_Ax, matrix_ps 
    TYPE(qs_rho_type), POINTER                         :: rho_b
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'hessian_op2_nonorth', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle, ispin, nspins
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_rho_type), POINTER                         :: rho
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_s, matrix_G
    TYPE(section_vals_type), POINTER                   :: input, xc_section
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc
    TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
    TYPE(pw_poisson_type), POINTER                     :: poisson_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                          v_hartree_rspace
    TYPE(pw_env_type), POINTER                         :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g_b
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_tmp
    REAL(KIND=dp)                                      :: ehartree
    TYPE(neighbor_list_set_p_type), DIMENSION(:), &
       POINTER                                         :: sab_orb

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    input = input, &
                    matrix_s_kp = matrix_s,&
                    para_env=para_env,&
                    rho = rho, &
                    sab_orb = sab_orb)
    nspins = dft_control%nspins
    IF(para_env%mepos==0) WRITE(*,*) "hessian_op2_nonorth BEGINNING"

    ! Get density matrix
    CALL qs_rho_get(rho, rho_ao_kp = rho_ao)

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

! XC-Kernel
    NULLIFY(v_xc, xc_section)
    xc_section => section_vals_get_subs_vals(input, "DFT%XC")

    ! add xc-kernel
    CALL calculate_kernel(qs_env, &
                          vxc = v_xc, &
                          rho = rho, &
                          rho1 = rho_b, &
                          xc_section = xc_section)
    DO ispin = 1,nspins
       CALL pw_scale(v_xc(ispin)%pw,v_xc(ispin)%pw%pw_grid%dvol)
       ! factor 2 here
       IF (nspins == 1) v_xc(1)%pw%cr3d = 2.0_dp*v_xc(1)%pw%cr3d 
    END DO

! Hartree contribution to Hessian

    ! rho_tot_gspace = rho_g
    ! take B density to build G^{H}[B]
    CALL qs_rho_get(rho_b, rho_g = rho_g_b)
    CALL pw_zero(rho_tot_gspace%pw)
    DO ispin = 1,nspins
       CALL pw_axpy(rho_g_b(ispin)%pw,rho_tot_gspace%pw)
    END DO
    
    ! test
    IF (nspins == 1) CALL pw_scale(rho_tot_gspace%pw, 2.0_dp)

    ! get Hartree potential from rho_tot_gspace
    CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree, &
                           vhartree = v_hartree_gspace%pw)
    CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
    CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

!    IF(para_env%mepos==0) WRITE(*,*) "hessian_op2_non_ortho_and_back - Ehartree", ehartree

    ! Add v_xc + v_H
    DO ispin = 1,nspins
       CALL pw_axpy(v_hartree_rspace%pw,v_xc(ispin)%pw)
    END DO
    ! need to be doen above and separately
    !IF (nspins == 1 ) v_xc(1)%pw%cr3d = 2.0_dp*v_xc(1)%pw%cr3d

    NULLIFY (matrix_G)
    CALL dbcsr_allocate_matrix_set(matrix_G, nspins, 1)

    ! Integrate this
    DO ispin = 1, nspins
       ! Init response kernel matrix
       ALLOCATE (matrix_G(ispin,1)%matrix)
       CALL dbcsr_copy(matrix_G(ispin,1)%matrix, matrix_s(1,1)%matrix, &
                       name="MATRIX Kernel")

       CALL dbcsr_create(matrix_G(ispin,1)%matrix, name="MATRIX Kernel", &
                         template = matrix_s(1,1)%matrix)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_G(ispin, 1)%matrix, sab_orb)
       CALL dbcsr_set(matrix_G(ispin,1)%matrix, 0.0_dp)

       ! matrix G(B)
       ! symmetric matrix coming out here
       ! needs proper setup, otherwise lower triangular blocks are set to zero
       ! when they are actually symmetric to upper triangular blocks.
       CALL integrate_v_rspace(v_rspace = v_xc(ispin), &
                               pmat = rho_ao(ispin,1), &
                               hmat = matrix_G(ispin, 1), &
                               qs_env = qs_env, &
                               calculate_forces = .FALSE., &
                               basis_type = "ORB")
    END DO

!    IF(para_env%mepos==0) WRITE(*,*) "hessian_op2_nonorth - matrix_G"
!    CALL dbcsr_print(matrix_G(1,1)%matrix)
    IF(ASSOCIATED(matrix_p)) THEN
!    IF(para_env%mepos==0) WRITE(*,*) "hessian_op2_nonorth - matrix_p"
    END IF
!    CALL dbcsr_print(matrix_p(1,1)%matrix)
!    IF(para_env%mepos==0) WRITE(*,*) "hessian_op2_nonorth - matrix_Ax1"
!    CALL dbcsr_print(matrix_Ax(1,1)%matrix)

    ! Hessian Ax already contains
    ! Ax1 = F*B*S - S*B*F
    ! Now adding 
    ! Ax = Ax1 + G(B)P*S - S*P*G(B)
    !    = Ax1 + G(B)*PS - (G(B)*PS)^T
    IF(para_env%mepos==0) WRITE(*,*) "hessian_op2_nonorth - ec_scomm for Ax2"
    DO ispin = 1, nspins
       CALL ec_scomm(qs_env = qs_env, &
                     matrix_a = matrix_G, &
                     matrix_ps = matrix_ps, &
                     matrix_out = matrix_Ax, &
                     eps_filter = eps_filter, &
                     alpha = 1.0_dp, &
                     beta = 1.0_dp,&
                     ps_trans = .FALSE.)
    END DO

    ! release pw grids
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
    DO ispin = 1,nspins
       CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
    END DO
    DEALLOCATE(v_xc)

    CALL dbcsr_deallocate_matrix_set(matrix_G)

    IF(para_env%mepos==0) WRITE(*,*) "hessian_op2_nonorth END"
    
    CALL timestop(handle)

  END SUBROUTINE hessian_op2_nonorth

! **************************************************************************************************
!> \brief Projector test P(M) = P_o*M*P_v^T + P_v*M*P_o^T
!>        D_b = [proj(cg),P]_S = proj([cg,P],S)
!>        with P_o = P*S
!>        with P_v = (1-P*S)
!>        preortho equivalent
!> qs_env ...
!> matrix_p ...
!> matrix_s ...
!> matrix_cg trial matrix, w/o projector applied to it 
!>          
!> \param qs_env ...
!> \date    02.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE projector_test(qs_env, matrix_ps)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_p_type), DIMENSION(:,:), &
                                        POINTER        :: matrix_ps

    CHARACTER(len=*), PARAMETER :: routineN = 'projector_test', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    REAL(KIND=dp)                                      :: eps_filter, frob_norm
    TYPE(dbcsr_type)                                   :: matrix_psps
    TYPE(dbcsr_type)                                   :: matrix_pv, matrix_pvpv

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    para_env=para_env)
    eps_filter = 1.0E-7_dp
    IF(para_env%mepos==0) WRITE(*,*) "projector_test BEGINNING"

    CALL dbcsr_create(matrix_psps, template = matrix_ps(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_pv, template = matrix_ps(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_pvpv, template = matrix_ps(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)

    ! matrix psps
    CALL dbcsr_multiply("N","N", 1.0_dp, matrix_ps(1,1)%matrix, matrix_ps(1,1)%matrix, &
                        1.0_dp, matrix_psps, filter_eps = eps_filter)

    ! PSPS - PS
    WRITE(*,*) "PS"
    CALL dbcsr_print(matrix_ps(1,1)%matrix, matlab_format = .TRUE.)
    WRITE(*,*) "PSPS"
    CALL dbcsr_print(matrix_psps, matlab_format = .TRUE.)


    CALL dbcsr_add(matrix_psps, matrix_ps(1,1)%matrix, 1.0_dp, -1.0_dp)

    frob_norm = 0.0_dp
    frob_norm = dbcsr_frobenius_norm(matrix_psps)
    WRITE(*,*) "|| PSPS - PS ||_F :", frob_norm



    ! Pv
    CALL dbcsr_copy(matrix_pv, matrix_ps(1,1)%matrix)
    CALL dbcsr_set(matrix_pv, 0.0_dp)
    CALL dbcsr_add_on_diag(matrix_pv, 1.0_dp)
    CALL dbcsr_add(matrix_pv, matrix_ps(1,1)%matrix, 1.0_dp, -1.0_dp)

    ! matrix pvpv
    CALL dbcsr_multiply("N","N", 1.0_dp, matrix_pv, matrix_pv, &
                        1.0_dp, matrix_pvpv, filter_eps = eps_filter)

    WRITE(*,*) "PV"
    CALL dbcsr_print(matrix_pv, matlab_format = .TRUE.)
    WRITE(*,*) "PVPV"
    CALL dbcsr_print(matrix_pvpv, matlab_format = .TRUE.)

    ! PvPv - Pv
    CALL dbcsr_add(matrix_pvpv, matrix_pv, 1.0_dp, -1.0_dp)

    frob_norm = 0.0_dp
    frob_norm = dbcsr_frobenius_norm(matrix_pvpv)
    WRITE(*,*) "|| PvPv - Pv ||_F :", frob_norm

    CALL dbcsr_release(matrix_psps)
    CALL dbcsr_release(matrix_pv)
    CALL dbcsr_release(matrix_pvpv)

    CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief Projector P(M) = P_o*M*P_v^T + P_v*M*P_o^T
!>        with P_o = P*S
!>        with P_v = (1-P*S)
!>        preortho equivalent
!> qs_env ...
!> matrix_ps
!> matrix_out Matrix to Matrix to which projector is applied.
!>          
!> \param qs_env ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE matrix_proj(qs_env, matrix_ps, matrix_out, eps_filter)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_p_type), DIMENSION(:,:), &
                            INTENT(IN), POINTER        :: matrix_ps
    TYPE(dbcsr_p_type), DIMENSION(:,:), &
                         INTENT(INOUT), POINTER        :: matrix_out
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_proj', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle, ispin, nspins
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_m

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    para_env=para_env)
    nspins = dft_control%nspins
    !IF(para_env%mepos==0) WRITE(*,*) "matrix_proj BEGINNING"

    ! init temp matrices
    CALL dbcsr_create(matrix_m, template = matrix_out(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_tmp1, template = matrix_out(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)

    ! P(M) = PSM - PSM(PS)^T + M(PS)^T - PSM(PS)^T
    ! P(M) = PSM + M(PS)^T - 2*PSM(PS)^T

    ! P(M) = PSM + MSP - 2*PSMSP

!    WRITE(*,*) "matrix_proj - matrix_out BEGINNING"
!    CALL dbcsr_print(matrix_out(1,1)%matrix)

    DO ispin = 1, nspins

       ! M
       CALL dbcsr_copy(matrix_m, matrix_out(ispin,1)%matrix)

       ! cg = PSM 
       CALL dbcsr_multiply("N","N", 1.0_dp, matrix_ps(ispin,1)%matrix, matrix_m, &
                           0.0_dp, matrix_out(ispin,1)%matrix, filter_eps = eps_filter)

!       WRITE(*,*) "cg = PSM"
!       CALL dbcsr_print(matrix_out(1,1)%matrix)

       ! tmp1 = MSP = M(PS)^T
       CALL dbcsr_multiply("N","T", 1.0_dp, matrix_m, matrix_ps(ispin,1)%matrix, &
                           0.0_dp, matrix_tmp1, filter_eps = eps_filter)

!       WRITE(*,*) "cg = MSP"
!       CALL dbcsr_print(matrix_tmp1)

       ! cg = cg + MSP = cg + tmp1
       CALL dbcsr_add(matrix_out(ispin,1)%matrix, matrix_tmp1, 1.0_dp, 1.0_dp)

!       WRITE(*,*) "cg = PSM + MSP"
!       CALL dbcsr_print(matrix_out(1,1)%matrix)

       ! cg = cg - 2*PSMSP = cg - 2*PS*tmp1
       CALL dbcsr_multiply("N","N", -2.0_dp, matrix_ps(ispin,1)%matrix, matrix_tmp1, &
                           0.0_dp, matrix_m, filter_eps = eps_filter)
!       WRITE(*,*) "cg = - 2*PSMSP"
!       CALL dbcsr_print(matrix_m)

       ! cg = cg - 2*PSMSP = cg - 2*PS*tmp1
       CALL dbcsr_multiply("N","N", -2.0_dp, matrix_ps(ispin,1)%matrix, matrix_tmp1, &
                           1.0_dp, matrix_out(ispin,1)%matrix, filter_eps = eps_filter)
       
!       WRITE(*,*) "P(M) = PSM + MSP - 2*PSMSP"
!       CALL dbcsr_print(matrix_out(1,1)%matrix)

    END DO
    
    CALL dbcsr_release(matrix_tmp1)
    CALL dbcsr_release(matrix_m)

    CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief Projector P_T(M) = P_o^T*M*P_v + P_v^T*M*P_o
!>        with P_o = P*S
!>        with P_v = (1-P*S)
!>        postortho equivalent
!> qs_env ...
!> matrix_ps ...
!> matrix_out Matrix to Matrix to which projector is applied.
!>          
!> \param qs_env ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE matrix_proj_t(qs_env, matrix_ps, matrix_out, eps_filter)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_p_type), DIMENSION(:,:), &
                            INTENT(IN), POINTER        :: matrix_ps
    TYPE(dbcsr_p_type), DIMENSION(:,:), &
                         INTENT(INOUT), POINTER        :: matrix_out
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_proj_t', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle, ispin, nspins
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_m

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    para_env=para_env)
    nspins = dft_control%nspins
    !IF(para_env%mepos==0) WRITE(*,*) "matrix_proj_t BEGINNING"

    ! init temp matrices
    CALL dbcsr_create(matrix_m, template = matrix_out(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_tmp1, template = matrix_out(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)

    ! P(M) = MPS + (PS)^T*M - 2*(PS)^T*M*PS
    ! P(M) = MPS + SPM - 2SPMPS

    DO ispin = 1, nspins

       ! M
       CALL dbcsr_copy(matrix_m, matrix_out(ispin,1)%matrix)

       ! cg = MPS 
       CALL dbcsr_multiply("N","N", 1.0_dp, matrix_m, matrix_ps(ispin,1)%matrix, &
                           0.0_dp, matrix_out(ispin,1)%matrix, filter_eps = eps_filter)

!       WRITE(*,*) "cg = MPS"
!       CALL dbcsr_print(matrix_out(1,1)%matrix)

       ! tmp1 = SP*M = (PS)^T*M
       CALL dbcsr_multiply("T","N", 1.0_dp, matrix_ps(ispin,1)%matrix, matrix_m, &
                           0.0_dp, matrix_tmp1, filter_eps = eps_filter)

!       WRITE(*,*) "cg = MPS + SP*M"
!       ! cg = cg + SP*M = cg + tmp1
       CALL dbcsr_add(matrix_out(ispin,1)%matrix, matrix_tmp1, 1.0_dp, 1.0_dp)

!       ! cg = cg - 2*SP*M*PS = cg - 2*tmp1*PS
!       CALL dbcsr_multiply("N","N", -2.0_dp, matrix_ps(ispin,1)%matrix, matrix_tmp1, &
!                           0.0_dp, matrix_m, filter_eps = eps_filter)
!       WRITE(*,*) "- 2*(PS)^T*M*PS"
!       CALL dbcsr_print(matrix_m)

       ! cg = cg - 2*SP*M*PS = cg - 2*tmp1*PS
       CALL dbcsr_multiply("N","N", -2.0_dp, matrix_tmp1, matrix_ps(ispin,1)%matrix, &
                           1.0_dp, matrix_out(ispin,1)%matrix, filter_eps = eps_filter)
       
!       WRITE(*,*) "P(M) = MPS + (PS)^T*M - 2*(PS)^T*M*PS"
!       CALL dbcsr_print(matrix_out(1,1)%matrix)

    END DO
    
    CALL dbcsr_release(matrix_tmp1)
    CALL dbcsr_release(matrix_m)

    CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief Projector P(M) = P_o*M*P_v^T + P_v*M*P_o^T
!>        with P_o = P*S
!>        with P_v = (1-P*S)
!>        preortho equivalent
!> qs_env ...
!> matrix_ps
!> matrix_out Matrix to Matrix to which projector is applied.
!>          
!> \param qs_env ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE matrix_proj_general(qs_env, matrix_ps, matrix_out, eps_filter)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_p_type), DIMENSION(:,:), &
                            INTENT(IN), POINTER        :: matrix_ps
    TYPE(dbcsr_p_type), DIMENSION(:,:), &
                         INTENT(INOUT), POINTER        :: matrix_out
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'matrix_proj_general', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle, ispin, nspins
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_tmp2, matrix_m

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    para_env=para_env)
    nspins = dft_control%nspins
    !IF(para_env%mepos==0) WRITE(*,*) "matrix_proj BEGINNING"

    ! init temp matrices
    CALL dbcsr_create(matrix_m, template = matrix_out(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_tmp2, template = matrix_out(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_tmp1, template = matrix_out(1,1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)

    ! P(M) = PSM - PSM(PS)^T + M(PS)^T - PSM(PS)^T
    ! P(M) = PSM + M(PS)^T - 2*PSM(PS)^T

    ! P(M) = PSM + MSP - 2*PSMSP

!    WRITE(*,*) "matrix_proj - matrix_out BEGINNING"
!    CALL dbcsr_print(matrix_out(1,1)%matrix)

    DO ispin = 1, nspins

       ! M
       CALL dbcsr_copy(matrix_m, matrix_out(ispin,1)%matrix)

       ! cg = PSM 
       CALL dbcsr_multiply("N","N", 1.0_dp, matrix_ps(ispin,1)%matrix, matrix_m, &
                           0.0_dp, matrix_tmp2, filter_eps = eps_filter)

       ! tmp1 = MSP = M(PS)^T
       CALL dbcsr_multiply("N","T", 1.0_dp, matrix_m, matrix_ps(ispin,1)%matrix, &
                           0.0_dp, matrix_tmp1, filter_eps = eps_filter)

!       ! cg = cg + MSP = cg + tmp1
       CALL dbcsr_add(matrix_tmp2, matrix_tmp1, 1.0_dp, 1.0_dp)

       ! cg = cg - 2*PSMSP = cg - 2*PS*tmp1
       CALL dbcsr_multiply("N","N", -2.0_dp, matrix_ps(ispin,1)%matrix, matrix_tmp1, &
                           0.0_dp, matrix_m, filter_eps = eps_filter)

       ! cg = cg - 2*PSMSP = cg - 2*PS*tmp1
       CALL dbcsr_multiply("N","N", -2.0_dp, matrix_ps(ispin,1)%matrix, matrix_tmp1, &
                           1.0_dp, matrix_tmp2, filter_eps = eps_filter)
       
       ! output
       ! Keep symmetry of matrix we applied projector on
       CALL dbcsr_copy(matrix_out(ispin,1)%matrix, matrix_tmp2, keep_sparsity=.TRUE.)

    END DO
    
    CALL dbcsr_release(matrix_tmp1)
    CALL dbcsr_release(matrix_tmp2)
    CALL dbcsr_release(matrix_m)

    CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief S-Commutator 
!>        matrix_out = alpha*[A,C]_B + beta*matrix_out 
!>                   = alpha*(ABC - CBA) + beta*matrix_out
!>
!> \param matrix_a ...
!> \param matrix_b ...
!> \param matrix_c ...
!> \param matrix_out Add S-commutator result onto matrix_out
!> \param alpha Multiplicator for s-commutator 
!> \param alpha Multiplicator for matrix_out 
!> \param ps_trans Use tranpose of PS matrix
!> \par History
!>       2019.9 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_scomm(qs_env, matrix_a, matrix_ps, matrix_out, eps_filter, alpha, beta, ps_trans)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:,:), INTENT(IN), &
         POINTER                                         :: matrix_a, matrix_ps
      TYPE(dbcsr_p_type), DIMENSION(:,:), INTENT(INOUT), &
         POINTER                                         :: matrix_out
      REAL(KIND=dp), INTENT(IN)                          :: eps_filter
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha, beta
      LOGICAL, INTENT(IN), OPTIONAL                      :: ps_trans

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_scomm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      TYPE(dbcsr_type)                                   :: matrix_tmp2, matrix_tmp1, ref_matrix
      TYPE(dbcsr_type)                                   :: matrix_tmp3
      TYPE(dft_control_type), POINTER                    :: dft_control
      REAL(KIND=dp)                                      :: myalpha, mybeta
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL                                            :: my_pst
      REAL(KIND=dp)                                      :: trace

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env = qs_env, &
                      dft_control = dft_control, &
                      para_env = para_env)
      nspins = dft_control%nspins
!      IF (para_env%mepos == 0) WRITE(*,*) "ec_scomm - BEGINNING"

      myalpha = 1.0_dp
      IF (PRESENT(beta)) myalpha = alpha
      mybeta = 0.0_dp
      IF (PRESENT(beta)) mybeta = beta
      my_pst = .FALSE.
      IF (PRESENT(ps_trans)) my_pst = ps_trans

      ! Ref matrix - Get symmetry of input methods
      ref_matrix = matrix_a(1,1)%matrix

      ! init temp matrices
      CALL dbcsr_create(matrix_tmp1, template = ref_matrix, matrix_type = dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp2, template = ref_matrix, matrix_type = dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp3, template = ref_matrix, matrix_type = dbcsr_type_no_symmetry)
      CALL dbcsr_set(matrix_tmp1, 0.0_dp)
      CALL dbcsr_set(matrix_tmp2, 0.0_dp)
      CALL dbcsr_set(matrix_tmp3, 0.0_dp)

      ! ABC - CBA 

      ! CASES : 
      ! all use symmetry
      ! res : G(dp)*P*S - S*P*G(dP)  
      ! B   : cg*S*P - P*S*cg         !! needs PS transposed
      ! Ax1 : F*B*S - S*B*F           !! needs PS to be BS
      ! Ax2 : G(B)*P*S - S*P*G(B)

      ! perform multiplication D = beta*D + alpha*(ABC - CBA)

      ! dbcsr_multiply: C = a*A*B + b*C

      DO ispin = 1, nspins
         IF(.NOT. my_pst) THEN
!            WRITE(*,*) "ec_scomm - just commutator"
            ! We could use symmetry here
            ! D = [A,C]_B = ABC - CBA = ABC - (ABC)^{T}
            ! tmp1 = ABC = < matrix_a > * PS
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_a(ispin, 1)%matrix, matrix_ps(ispin,1)%matrix, &
                                0.0_dp, matrix_tmp1, filter_eps = eps_filter)

!            WRITE(*,*) "ec_scomm - matrix_abc"
!            CALL dbcsr_print(matrix_tmp1, matlab_format = .TRUE.)

            ! tmp_2 = (tmp1)^T
            CALL dbcsr_transposed(matrix_tmp2, matrix_tmp1)
            ! tmp1 = tmp1 - (tmp1)^T
            CALL dbcsr_add(matrix_tmp1, matrix_tmp2, 1.0_dp, -1.0_dp)

!            WRITE(*,*) "ec_scomm - matrix_out before fill"
!            CALL dbcsr_print(matrix_out(1,1)%matrix)

            ! D = beta*D + alpha*(ABC - (ABC)^T)
            ! output antisymmetric matrix
            CALL dbcsr_add(matrix_out(ispin,1)%matrix, matrix_tmp1, mybeta, myalpha)


            ! Test anti-symmetry of output matrix
            CALL dbcsr_set(matrix_tmp2, 0.0_dp)
            CALL dbcsr_copy(matrix_tmp3, matrix_out(1,1)%matrix)
            CALL dbcsr_transposed(matrix_tmp2, matrix_tmp3)
            CALL dbcsr_add(matrix_tmp3, matrix_tmp2, 1.0_dp, 1.0_dp)
            trace = dbcsr_frobenius_norm(matrix_tmp3)
            WRITE(*,*) "ec_scomm - output matrix (anti-symmetric), || out - out^T ||_F  ", trace


         ELSE
            ! Case for B and Z matrix - i.e. matrix_a is ANTI-SYM 
            ! [cg, P]_S = cg*S*P - P*S*cg
            ! tmp1 = A*t

            ! tmp1 = cg*(PS)^T = cg*S*P = A*B*C
            CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_a(ispin, 1)%matrix, matrix_ps(ispin,1)%matrix, &
                                0.0_dp, matrix_tmp1, filter_eps = eps_filter)

            ! tmp_2 = (tmp1)^T
            CALL dbcsr_transposed(matrix_tmp2, matrix_tmp1)
            ! tmp1 = tmp1 - (tmp1)^T
            CALL dbcsr_add(matrix_tmp1, matrix_tmp2, 1.0_dp, 1.0_dp)

            ! Test symmetry of output matrix
            CALL dbcsr_set(matrix_tmp2, 0.0_dp)
            CALL dbcsr_copy(matrix_tmp3, matrix_tmp1)
            CALL dbcsr_transposed(matrix_tmp2, matrix_tmp3)
            CALL dbcsr_add(matrix_tmp3, matrix_tmp2, 1.0_dp, -1.0_dp)
            trace = dbcsr_frobenius_norm(matrix_tmp3)
            WRITE(*,*) "ec_scomm - output matrix (symmetric), || out - out^T ||_F  ", trace

            ! D = ABC + (ABC)^T
            ! Here, beta = 0 and alpha = 1 ALWAYS
            !CALL dbcsr_copy(matrix_out(ispin,1)%matrix, matrix_tmp1, keep_sparsity = .TRUE.)
            CALL dbcsr_copy(matrix_out(ispin,1)%matrix, matrix_tmp1)
         END IF

         CALL dbcsr_filter(matrix_out(ispin,1)%matrix, eps_filter)

      ENDDO

      CALL dbcsr_release(matrix_tmp1)
      CALL dbcsr_release(matrix_tmp2)
      CALL dbcsr_release(matrix_tmp3)

      CALL timestop(handle)

   END SUBROUTINE ec_scomm

! **************************************************************************************************
!> \brief   Calculation of second derivative potential
!> \param   qs_env ...
!> \param   vxc [OUT], contains partially integrated second derivative
!>          taken with respect to rho, evluated in rho and folded with rho1
!> \param
!> \param
!> \date    11.2019
!> \author  fbelle
! **************************************************************************************************
   SUBROUTINE calculate_kernel(qs_env, vxc, rho, rho1, xc_section)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc
      TYPE(qs_rho_type), POINTER                         :: rho, rho1
      TYPE(section_vals_type), POINTER                   :: xc_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_kernel', routineP = moduleN//':'//routineN

      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: handle, ispin, nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      LOGICAL                                            :: lsd
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, rho1_r, rho1_g, tau_pw, vxc_rho
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(xc_rho_set_type), POINTER                     :: rho_set, rho1_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set

      CALL timeset(routineN, handle)
      NULLIFY(logger)
      logger => cp_get_default_logger()

      NULLIFY (auxbas_pw_pool, pw_env, vxc_rho)

      CALL get_qs_env(qs_env, &
                      dft_control = dft_control,&
                      para_env = para_env,&
                      pw_env=pw_env)
      nspins = dft_control%nspins

      CALL pw_env_get(pw_env = pw_env, &
                      auxbas_pw_pool = auxbas_pw_pool)

      ALLOCATE(vxc_rho(nspins))
      DO ispin = 1, nspins
         NULLIFY(vxc_rho(nspins)%pw)
         CALL pw_pool_create_pw(pool = auxbas_pw_pool, &
                                pw = vxc_rho(ispin)%pw, &
                                use_data = REALDATA3D, &
                                in_space = REALSPACE)
         CALL pw_zero(vxc_rho(ispin)%pw)
      END DO

      ! Get grid based density
      NULLIFY(rho_r, rho1_r, rho1_g)
      CALL qs_rho_get(rho, rho_r = rho_r)
      CALL qs_rho_get(rho1, rho_r = rho1_r, rho_g = rho1_g)

      NULLIFY (deriv_set, rho_set, rho1_set, tau_pw)
      ! Init rho1_set, with which to fold second derivative
      bo = rho_r(1)%pw%pw_grid%bounds_local
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      lsd = (nspins == 2)
      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

      ! calculate the rho set used to fold with 2nd derivative
      CALL xc_rho_set_update(rho1_set, rho1_r, rho1_g, tau_pw, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool)

      ! main ingredient is xc_rho_set_and_dset_create
      CALL xc_prep_2nd_deriv(deriv_set = deriv_set, &    ! containing potentials
                             rho_set = rho_set, &        ! density at which derivs are calculated
                             rho_r = rho_r, &            ! place where derivative is evaluated
                             pw_pool = auxbas_pw_pool, & ! pool for grids
                             xc_section = xc_section)

      ! evaluation of 2nd deriv in rho_set density
      ! folding of second deriv with density in rho1_set
      CALL xc_calc_2nd_deriv(v_xc = vxc_rho, &           ! XC-potential
                             deriv_set = deriv_set, &    ! deriv of xc-potential
                             rho_set = rho_set, &        ! density at which deriv are calculated
                             rho1_set = rho1_set, &      ! density with which to fold
                             pw_pool = auxbas_pw_pool, & ! pool for grids
                             xc_section = xc_section, &
                             gapw = .FALSE.)

      ! Release second deriv stuff
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho_set = rho_set, pw_pool = auxbas_pw_pool)
      CALL xc_rho_set_release(rho_set = rho1_set, pw_pool = auxbas_pw_pool)

      ! export vxc
      IF (ASSOCIATED(vxc_rho)) THEN
         ALLOCATE (vxc(nspins))
         DO ispin = 1, nspins
            vxc(ispin)%pw => vxc_rho(ispin)%pw
!            CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
         END DO
      END IF

!      DO ispin = 1, nspins
!         CALL pw_pool_give_back_pw(auxbas_pw_pool, vxc_rho(ispin)%pw)
!      END DO
      DEALLOCATE (vxc_rho)

!      IF(para_env%mepos==0) WRITE(*,*) "calculate_kernel END"

      CALL timestop(handle)

END MODULE ec_ao_res_solver
