!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Harris functional force driver routine
!>
!>
!> \date 09.2019
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_ao_solver
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_norm,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace,&
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
                                              dbcsr_copy,&
                                              dbcsr_copy_into_existing,&
                                              dbcsr_create,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_distribution_type,&
                                              dbcsr_dot,&
                                              dbcsr_filter,&
                                              dbcsr_frobenius_norm,&
                                              dbcsr_get_info,&
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
   USE ec_env_types,                    ONLY: energy_correction_type 
   USE kg_correction,                   ONLY: create_kernel
   USE input_constants,                 ONLY: kg_tnadd_embed,&
                                              ls_s_sqrt_proot,&
                                              ls_s_sqrt_ns,&
                                              ot_precond_full_single_inverse
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kg_correction,                   ONLY: kg_ekin_subset
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_types,                 ONLY: linres_control_create,&
                                              linres_control_release,&
                                              linres_control_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz,& 
                                              matrix_sqrt_proot
   USE qs_linres_methods,               ONLY: p_env_check_i_alloc
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

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_ao_solver'

! *** Public subroutines ***

   PUBLIC :: calculate_2whz_ao_matrix
   PUBLIC :: ec_res_ao
   PUBLIC :: ec_wz_matrix

CONTAINS

! **************************************************************************************************
!> \brief   
!>          
!>  without passing to orthonormal basis         
!>          
!>          
!> \param qs_env ...
!> \date    01.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_res_ao(qs_env, ec_env, p_env, iounit, should_stop)

     TYPE(qs_environment_type), POINTER                 :: qs_env
     TYPE(energy_correction_type), POINTER              :: ec_env
     TYPE(qs_p_env_type), POINTER                       :: p_env
     INTEGER, INTENT(IN)                                :: iounit
     LOGICAL, INTENT(OUT)                               :: should_stop

     CHARACTER(len=*), PARAMETER :: routineN = 'ec_res_ao', routineP = moduleN//':'//routineN
     
     INTEGER                                            :: i, handle, ispin, max_iter_lanczos, &
                                                           nao, nspins, s_sqrt_method, s_sqrt_order
     LOGICAL                                            :: do_precond, restart
     REAL(KIND=dp)                                      :: eps_filter, eps_lanczos, focc, norm_res, &
                                                           threshold, t1, t2
     TYPE(cp_para_env_type), POINTER                    :: para_env
     TYPE(dft_control_type), POINTER                    :: dft_control
     TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_Ax, matrix_cg, matrix_cg_z, &
                                                           matrix_ps, matrix_s, matrix_q, &
                                                           matrix_res, matrix_z, matrix_z0
     REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha, beta, new_norm, norm_cA, &
                                                           norm_rr, tr_rz00
     REAL(KIND=dp)                                      :: norm_zao, norm_zmo
     TYPE(linres_control_type), POINTER                 :: linres_control
     TYPE(section_vals_type), POINTER                   :: solver_section
     ! SQRT crap
     TYPE(dbcsr_type)                                   :: matrix_s_sqrt, matrix_s_sqrt_inv 
     TYPE(dbcsr_type)                                   :: matrix_oldz

     CALL timeset(routineN, handle)

     NULLIFY (alpha, beta, dft_control, linres_control, new_norm, norm_cA, norm_rr, tr_rz00, &
              solver_section, matrix_s, para_env)

     t1 = m_walltime()

     CALL get_qs_env(qs_env,para_env=para_env)
     IF(para_env%mepos==0) WRITE(*,*) repeat("-",80)
     IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao BEGINNING"
     IF(para_env%mepos==0) WRITE(*,*) repeat("-",80)

     CALL get_qs_env(qs_env = qs_env, &
                     dft_control = dft_control, &
                     linres_control = linres_control, &
                     matrix_s = matrix_s)
     nspins = dft_control%nspins

     ! Get number of atomic orbitals, to scale convergence criteria
     CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

     ! Gerneral setup 
     do_precond = .TRUE.

     IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - eps", linres_control%eps
     eps_filter = linres_control%eps_filter
     IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - eps_filter", eps_filter

     ALLOCATE(alpha(nspins),beta(nspins), new_norm(nspins), norm_cA(nspins), norm_rr(nspins))
     ALLOCATE(tr_rz00(nspins))

     NULLIFY(matrix_Ax, matrix_cg, matrix_cg_z, matrix_res, matrix_z, matrix_z0)
     CALL dbcsr_allocate_matrix_set(matrix_Ax, nspins)
     CALL dbcsr_allocate_matrix_set(matrix_cg, nspins)
     CALL dbcsr_allocate_matrix_set(matrix_cg_z, nspins)
     CALL dbcsr_allocate_matrix_set(matrix_res, nspins)
     CALL dbcsr_allocate_matrix_set(matrix_z, nspins)
     CALL dbcsr_allocate_matrix_set(matrix_z0, nspins)
     DO ispin = 1, nspins
        ALLOCATE (matrix_Ax(ispin)%matrix)
        ALLOCATE (matrix_cg(ispin)%matrix)
        ALLOCATE (matrix_cg_z(ispin)%matrix)
        ALLOCATE (matrix_res(ispin)%matrix)
        ALLOCATE (matrix_z(ispin)%matrix)
        ALLOCATE (matrix_z0(ispin)%matrix)
        CALL dbcsr_create(matrix_Ax(ispin)%matrix, name = "linop MATRIX", &
                          template = matrix_s(1)%matrix,&
                          matrix_type=dbcsr_type_no_symmetry)
        CALL dbcsr_create(matrix_cg(ispin)%matrix, name = "TRIAL MATRIX", &
                          template = matrix_s(1)%matrix,&
                          matrix_type=dbcsr_type_no_symmetry)
        CALL dbcsr_create(matrix_cg_z(ispin)%matrix, name = "MATRIX CG-Z", &
                          template = matrix_s(1)%matrix,&
                          matrix_type=dbcsr_type_no_symmetry)
        CALL dbcsr_create(matrix_res(ispin)%matrix, name = "RESIDUE", &
                          template = matrix_s(1)%matrix,&
                          matrix_type=dbcsr_type_no_symmetry)
        CALL dbcsr_create(matrix_z(ispin)%matrix, name = "MATRIX Z", &
                          template = matrix_s(1)%matrix,&
                          matrix_type=dbcsr_type_no_symmetry)
        CALL dbcsr_create(matrix_z0(ispin)%matrix, name = "p after precondi-Matrix", &
                          template = matrix_s(1)%matrix,&
                          matrix_type=dbcsr_type_no_symmetry)
     END DO

     ! Matrix PS and Projector Q
     IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - build projector P and Q"
     NULLIFY(matrix_ps, matrix_q)
     CALL projector_p_and_q(qs_env, matrix_ps, matrix_q, eps_filter)

     ! Test of Z matrix vs MO result
     norm_zmo = dbcsr_frobenius_norm(ec_env%p_env%p1(1)%matrix)
     !CALL dbcsr_desymmetrize(ec_env%p_env%p1(1)%matrix, matrix_oldz)

     ! Preconditioner
     IF (do_precond) THEN

        ! Get info about preconditioner
        solver_section => section_vals_get_subs_vals(qs_env%input, "DFT%ENERGY_CORRECTION%LS_SOLVER")
        CALL section_vals_val_get(solver_section, "S_SQRT_METHOD", i_val=s_sqrt_method)
        CALL section_vals_val_get(solver_section, "S_SQRT_ORDER", i_val=s_sqrt_order)
        CALL section_vals_val_get(solver_section, "EPS_LANCZOS", r_val=eps_lanczos)
        CALL section_vals_val_get(solver_section, "MAX_ITER_LANCZOS", i_val=max_iter_lanczos)
        ! Avoid conversion issues for inversion of overlap matrix
        threshold = MAX(eps_filter, 1.0e-12_dp)

        CALL dbcsr_create(matrix_s_sqrt, template=matrix_s(1)%matrix, &
                               matrix_type=dbcsr_type_no_symmetry)
        CALL dbcsr_create(matrix_s_sqrt_inv, template=matrix_s(1)%matrix, &
                               matrix_type=dbcsr_type_no_symmetry)

        SELECT CASE (s_sqrt_method)
        CASE (ls_s_sqrt_proot)
            CALL matrix_sqrt_proot(matrix_s_sqrt, matrix_s_sqrt_inv, &
                                   matrix_s(1)%matrix, threshold, &
                                   s_sqrt_order, eps_lanczos, max_iter_lanczos, symmetrize=.TRUE.)
        CASE (ls_s_sqrt_ns)
           CALL matrix_sqrt_Newton_Schulz(matrix_s_sqrt, matrix_s_sqrt_inv, &
                                          matrix_s(1)%matrix, threshold, &
                                          s_sqrt_order, eps_lanczos, max_iter_lanczos)
        CASE DEFAULT
           CPABORT("Unknown sqrt method.")
        END SELECT

     END IF

!     !---------------
!     ! Test of hessian routine
!     ! Calculate hessian with Z-Matrix obtained from MO-solver
!     IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - CALL build_hessian_op_from_z"
!     CALL build_hessian_op_from_z(qs_env, matrix_ps, matrix_q, matrix_Ax, eps_filter)
!     !--------------

     !----------------------------------------
     ! Get righ-hand-side (RHS) operator
     !----------------------------------------

     ! Spin factor
     focc = -2.0_dp
     IF (nspins == 1) focc = -4.0_dp

     ! E^[1]_Harris = G(\delta P)P*S - S*P*G(\delta P)
     ! ABC-CBA
     ! RES ANTISYM
     IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - build RHS"
     CALL ec_scomm(qs_env = qs_env, &
                   matrix_a = ec_env%matrix_hz, & 
                   matrix_ps = matrix_ps, &
                   matrix_out = matrix_res, &
                   eps_filter = eps_filter, &
                   alpha = focc, & ! Beware the factor -1.0_dp
                   beta = 0.0_dp, &
                   ps_trans = .FALSE.)

     ! Frist trial matrix 
     DO ispin = 1, nspins
        ! Set first trial matrix to RHS
        IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - Initital guess cg_z0 = b"
        CALL dbcsr_copy(matrix_cg_z(ispin)%matrix, matrix_res(ispin)%matrix)
        ! Set first trial matrix to 0
!        IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - Initital guess cg_z0 = 0"
!        CALL dbcsr_set(matrix_cg_z(ispin)%matrix, 0.0_dp)
     END DO

     !IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - PROJ(cg_z0)"
     CALL projector(qs_env, matrix_ps, matrix_q, matrix_cg_z, eps_filter, .FALSE.)

     ! Ax
     CALL build_hessian_op(qs_env = qs_env, &
                           p_env = p_env,&
                           matrix_ps = matrix_ps, &   ! IN
                           matrix_cg = matrix_cg_z, & ! IN trial matrix
                           matrix_Ax = matrix_Ax, &   ! OUT Hessian
                           eps_filter = eps_filter)
                               
     ! r_0 = b - Ax0
     DO ispin = 1, nspins
        CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -1.0_dp)
     END DO

     ! Apply P(r) = P^T*r*Q + Q^T*r*P 
     !IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - Proj_T(r_0)"
     CALL projector(qs_env, matrix_ps, matrix_q, matrix_res, eps_filter, .TRUE.)

     ! Preconditioner
     IF (do_precond) THEN

        ! z_0 = S^-1/2 * res * S^-1/2
        IF(para_env%mepos==0) WRITE(*,*) "ec_ao_res - PRECONDITIONER z_0 = S^-1/2 * M * S^-1/2"
        DO ispin = 1, nspins
            CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
            CALL apply_matrix_preconditioner(matrix_z0(ispin)%matrix, "forward", &
                                             matrix_s_sqrt, matrix_s_sqrt_inv)
        END DO
     ELSE
        ! z_0 = r_0
        DO ispin = 1, nspins
           CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
        END DO
     END IF

     norm_res = 0.0_dp
     ! Tr(r_0 * z_0)
     DO ispin = 1, nspins

        ! cg = p_0 = z_0
        CALL dbcsr_copy(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix)

        ! Tr(r_0 * z_0)
        CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_cg(ispin)%matrix, norm_rr(ispin))
        IF(para_env%mepos==0) WRITE(*,*) "Tr(r_0 * z_0) (norm_rr)", norm_rr(1)

        IF (norm_rr(ispin) .LT. 0.0_dp) CPABORT("norm_rr < 0")
        norm_res = MAX(norm_res, norm_rr(ispin)/REAL(nao, dp))
        IF(para_env%mepos==0) WRITE(*,*) "norm_res", norm_res
     END DO

     alpha(:) = 0.0_dp

     ! header
     IF (iounit > 0) THEN
        IF(para_env%mepos==0) WRITE (iounit, "(/,T3,A,T16,A,T25,A,T38,A,/,T3,A)") &
           "Iteration", "Stepsize", "Convergence", "Time", &
           REPEAT("-", 80)
     ENDIF

     ! start iteration
     restart = .FALSE.
     iteration: DO i = 1, linres_control%max_iter

        IF(para_env%mepos==0) WRITE(*,*) "START ITERATION", i

        ! Convergence criteria
        linres_control%converged = .FALSE.
        ! default for eps 10E-6 in MO_linres
        IF (norm_res .LT. linres_control%eps) THEN
           linres_control%converged = .TRUE.
        ENDIF

        t2 = m_walltime()
        IF (iounit > 0) THEN
           IF(para_env%mepos==0) WRITE(*,*) i, linres_control%flag, MAXVAL(alpha), norm_res, t2-t1
           CALL m_flush(iounit)
        ENDIF

        IF (linres_control%converged) THEN
           IF (iounit > 0) THEN
              WRITE (iounit, "(/,T2,A,I4,A,/)") "The linear solver converged in ", i, " iterations."
              CALL m_flush(iounit)
           ENDIF
           EXIT iteration
        ELSE IF (should_stop) THEN
           IF (iounit > 0) THEN
              WRITE (iounit, "(/,T2,A,I4,A,/)") "The linear solver did NOT converge! External stop"
              CALL m_flush(iounit)
           END IF
           EXIT iteration
        ENDIF

        ! Ax = [F,S]_B + [G(B),S]_P and B = [X,P]_S
        CALL build_hessian_op(qs_env = qs_env, &
                              p_env = p_env,&
                              matrix_ps = matrix_ps, & ! IN
                              matrix_cg = matrix_cg, & ! IN trial matrix
                              matrix_Ax = matrix_Ax, & ! OUT (Ax(cg))
                              eps_filter = eps_filter)

        ! Apply projector on Ax
        !IF(para_env%mepos==0) WRITE(*,*) "ec_ao_res - Proj_T(Ax)"
        CALL projector(qs_env, matrix_ps, matrix_q, matrix_Ax, eps_filter, .TRUE.)

        DO ispin = 1, nspins

           CALL dbcsr_filter(matrix_Ax(ispin)%matrix, eps_filter)

           ! tr(Ap_j * p_j)
           CALL dbcsr_dot(matrix_cg(ispin)%matrix, matrix_Ax(ispin)%matrix, norm_cA(ispin))
           IF(para_env%mepos==0) WRITE(*,*) "norm_cA", norm_cA(ispin)

        END DO 

        IF (norm_cA(ispin) .LT. 0.0_dp) THEN

           DO ispin = 1, nspins

              ! Recalculate w/o preconditioner
              IF (i > 1) THEN
                 ! p_j = -z_j + beta*p_j
                 CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix, beta(ispin), -1.0_dp)
                 CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_res(ispin)%matrix, new_norm(ispin))
                 beta(ispin) = new_norm(ispin)/tr_rz00(ispin)
                 CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_res(ispin)%matrix, beta(ispin), 1.0_dp)
                 norm_rr(ispin) = new_norm(ispin)
              ELSE
                 CALL dbcsr_copy(matrix_res(ispin)%matrix, matrix_cg(ispin)%matrix)
                 CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_res(ispin)%matrix, norm_rr(ispin))
              END IF
        
           END DO 

           CALL build_hessian_op(qs_env = qs_env, &
                                 p_env = p_env, &
                                 matrix_ps = matrix_ps, &
                                 matrix_cg = matrix_cg, & ! cg
                                 matrix_Ax = matrix_Ax, &
                                 eps_filter = eps_filter)

           ! Matrix projector T
           CALL projector(qs_env, matrix_ps, matrix_q, matrix_Ax, eps_filter, .TRUE.)

           DO ispin = 1, nspins
              CALL dbcsr_dot(matrix_cg(ispin)%matrix, matrix_Ax(ispin)%matrix, norm_cA(ispin))
           END DO 

           CPABORT("tr(Ap_j*p_j) < 0")

        END IF

        DO ispin = 1, nspins
           ! Determine step-size
           ! alpha = norm_rr/norm_cA
           IF (norm_cA(ispin) .LT. linres_control%eps) THEN
              alpha(ispin) = 1.0_dp
           ELSE
              alpha(ispin) = norm_rr(ispin)/norm_cA(ispin)
           END IF
           IF(para_env%mepos==0) WRITE(*,*) "alpha", alpha(ispin)

           ! x_j+1 = x_j + alpha*p_j
           ! save response-denisty of this iteration
           CALL dbcsr_add(matrix_cg_z(ispin)%matrix, matrix_cg(ispin)%matrix, 1.0_dp, alpha(ispin))
           ! r_j+1 = r_j - alpha * Ap_j
!           CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -alpha(ispin))
        END DO

        ! need to recompute the residue
        restart = .FALSE.
        IF (MOD(i, linres_control%restart_every) .EQ. 0) THEN
           !
           ! r_j+1 = b - A * x_j+1
           CALL build_hessian_op(qs_env = qs_env, &
                                 p_env = p_env, &
                                 matrix_ps = matrix_ps, &
                                 matrix_cg = matrix_cg_z, & ! cg
                                 matrix_Ax = matrix_Ax, &
                                 eps_filter = eps_filter)

           CALL ec_scomm(qs_env = qs_env, &
                         matrix_a = ec_env%matrix_hz, & 
                         matrix_ps = matrix_ps, &
                         matrix_out = matrix_res, &
                         eps_filter = eps_filter, &
                         alpha = focc, & ! Beware the factor -1.0_dp
                         beta = 0.0_dp, &
                         ps_trans = .FALSE.)

           DO ispin = 1, nspins
              CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -1.0_dp)
           END DO

           CALL projector(qs_env, matrix_ps, matrix_q, matrix_res, eps_filter, .TRUE.)
           !
           restart = .TRUE.
        ELSE
           ! proj Ap onto the virtual subspace
           CALL projector(qs_env, matrix_ps, matrix_q, matrix_ax, eps_filter, .TRUE.)
           !
           ! r_j+1 = r_j - alpha * Ap_j
           DO ispin = 1, nspins
              CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -alpha(ispin))
           ENDDO
           restart = .FALSE.
        ENDIF

        !------------------------------------
        ! Testing of Z-Matrix
        CALL dbcsr_copy(matrix_z0(1)%matrix, matrix_cg_z(1)%matrix)
        CALL projector(qs_env, matrix_ps, matrix_q, matrix_z0, eps_filter, .FALSE.)
        ! Build intermediate  Z-Matrix 
        CALL ec_scomm(qs_env = qs_env, &
                      matrix_a = matrix_z0, &
                      matrix_ps = matrix_ps, &
                      matrix_out = p_env%p1, &
                      eps_filter = eps_filter, &
                      alpha = 0.5_dp, &
                      beta = 0.0_dp,&
                      ps_trans = .TRUE.)
        norm_zao = dbcsr_frobenius_norm(p_env%p1(1)%matrix)
        IF (para_env%mepos==0) WRITE(*,*) "|| Z (MO) ||_F :", norm_zmo
        IF (para_env%mepos==0) WRITE(*,*) "|| Z (AO) ||_F :", norm_zao
        IF (para_env%mepos==0) WRITE(*,*) "It, conv, ||Z(AO)||_F :", i, norm_res, norm_zao
        !------------------------------------

        IF (do_precond) THEN
           ! z_0 = S^-1/2 * M * S^-1/2
           !IF(para_env%mepos==0) WRITE(*,*) "ec_ao_res - PRECONDITIONER z_0 = S^-1/2 * M * S^-1/2"
           ! M*res = S^-1/2*res*S^-1/2
           DO ispin = 1, nspins
              CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
              ! M*res = S^-1/2*res*S^-1/2
              CALL apply_matrix_preconditioner(matrix_z0(ispin)%matrix, "forward", &
                                               matrix_s_sqrt, matrix_s_sqrt_inv)
           END DO
        ELSE
           DO ispin = 1, nspins
              CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
           END DO
        END IF

        norm_res = 0.0_dp

        DO ispin = 1, nspins
           ! Tr[r_j+1*r_j+1]
           CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_z0(ispin)%matrix, new_norm(ispin))
           IF (new_norm(ispin) .LT. 0.0_dp) CPABORT("tr(r_j+1*z_j+1) < 0")
           norm_res = MAX(norm_res, new_norm(ispin)/REAL(nao,dp))
           IF (para_env%mepos==0) WRITE(*,*) "norm_res", norm_res
           !IF (para_env%mepos==0) WRITE(*,*) "norm_rr", norm_rr(1)
           !IF (para_env%mepos==0) WRITE(*,*) "new_norm", new_norm(1)

           IF (norm_rr(ispin) .LT. linres_control%eps .OR. new_norm(ispin) .LT. linres_control%eps) THEN
              beta(ispin) = 0.0_dp
              linres_control%converged = .TRUE.
           ELSE
              beta(ispin) = new_norm(ispin)/norm_rr(ispin)
           END IF
           IF (para_env%mepos==0) Write(*,*) "beta", beta(1)

           ! update new search vector (matrix)
           ! p_j+1 = z_j+1 + beta*p_j
           CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix, beta(ispin), 1.0_dp)
           CALL dbcsr_filter(matrix_cg(ispin)%matrix, eps_filter)

           ! necessary if tr[pAp] negative and recalculated 
           tr_rz00(ispin) = norm_rr(ispin)
           norm_rr(ispin) = new_norm(ispin)
        END DO
     END DO iteration 

     !IF (para_env%mepos==0) WRITE(*,*) "ec_res_ao - Proj(cg_z)"
     CALL projector(qs_env, matrix_ps, matrix_q, matrix_cg_z, eps_filter, .FALSE.)

     ! Build final Z-Matrix 
     ! Z = [cg_z, P]_S = cg_z*S*P - P*S*cg_z
     !                 = cg_z*(P*S)^T + (cg_z*(PS)^T)^T actual calculation
     IF (para_env%mepos==0) WRITE(*,*) "ec_res_ao - FINAL Z AO"
     CALL ec_scomm(qs_env = qs_env, &
                   matrix_a = matrix_cg_z , &
                   matrix_ps = matrix_ps, &
                   matrix_out = p_env%p1, &
                   eps_filter = eps_filter, &
                   alpha = 0.5_dp, &
                   beta = 0.0_dp,&
                   ps_trans = .TRUE.)

     norm_zao = dbcsr_frobenius_norm(p_env%p1(1)%matrix)
     IF (para_env%mepos==0) WRITE(*,*) "|| Z (MO) ||_F :", norm_zmo
     IF (para_env%mepos==0) WRITE(*,*) "|| Z (AO) ||_F :", norm_zao
     norm_res = norm_zmo - norm_zao
     IF (para_env%mepos==0) WRITE(*,*) "||Z(MO)||_F - ||Z(AO)||_F :", norm_res

     CALL dbcsr_dot(p_env%p1(1)%matrix, matrix_s(1)%matrix, norm_res)
     IF(para_env%mepos==0) WRITE(*,*) "Tr[PzS]", norm_res

     !----------------------------------------------------------
     ! Calculate energy-weighted response density matrix
     ! AO: W_z_in = 0.5*(Z*KS*P + P*KS*Z)
     !----------------------------------------------------------

     ! Wz = 0.25*(Z*KS*P + P*KS*Z)
     ! outputs a matrix in symmetric type
     CALL ec_wz_matrix(qs_env, p_env%p1, p_env%w1, eps_filter)

     ! Release matrices
     CALL dbcsr_deallocate_matrix_set(matrix_Ax)
     CALL dbcsr_deallocate_matrix_set(matrix_cg)
     CALL dbcsr_deallocate_matrix_set(matrix_cg_z)
     CALL dbcsr_deallocate_matrix_set(matrix_ps)
     CALL dbcsr_deallocate_matrix_set(matrix_res)
     CALL dbcsr_deallocate_matrix_set(matrix_q)
     CALL dbcsr_deallocate_matrix_set(matrix_z)
     CALL dbcsr_deallocate_matrix_set(matrix_z0)
     CALL dbcsr_release(matrix_oldz)

     IF (do_precond) THEN
        CALL dbcsr_release(matrix_s_sqrt)
        CALL dbcsr_release(matrix_s_sqrt_inv)
     END IF

     !CALL linres_control_release(linres_control)

     DEALLOCATE(alpha, beta, new_norm, norm_cA, norm_rr)

     IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - END"

     CALL timestop(handle)

  END SUBROUTINE

! **************************************************************************************************
!> \brief Compute matrix_wz as needed for the forces
!>        Wz = 0.5*(Z*KS*P + P*KS*Z) (closed-shell)
!> \param qs_env ...
!> \param matrix_z ...
!> \param matrix_wz The energy weighted response-density matrix
!> \par History
!>       2020.2 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_wz_matrix(qs_env, matrix_z, matrix_wz, eps_filter)

     TYPE(qs_environment_type), POINTER                 :: qs_env
     TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
                                         POINTER        :: matrix_z
     TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
                                         POINTER        :: matrix_wz
     REAL(KIND=dp), INTENT(IN)                          :: eps_filter

     CHARACTER(len=*), PARAMETER :: routineN = 'ec_wz_matrix', &
        routineP = moduleN//':'//routineN

     INTEGER                                            :: handle, ispin, nspins
     REAL(KIND=dp)                                      :: scaling
     TYPE(dbcsr_type)                                   :: matrix_tmp, matrix_tmp2
     TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, matrix_p
     TYPE(dft_control_type), POINTER                    :: dft_control
     TYPE(cp_para_env_type), POINTER                    :: para_env
     TYPE(qs_rho_type), POINTER                         :: rho
     REAL(KIND=dp)                                      :: norm_res

     CALL timeset(routineN, handle)

     CALL get_qs_env(qs_env = qs_env, &
                     dft_control = dft_control, &
                     matrix_ks = matrix_ks, &
                     matrix_s = matrix_s, &
                     rho=rho,&
                     para_env = para_env)
     nspins = dft_control%nspins

     CALL qs_rho_get(rho, rho_ao = matrix_p)

     ! Init temp matrices
     CALL dbcsr_create(matrix_tmp, template = matrix_z(1)%matrix, &
                       matrix_type = dbcsr_type_no_symmetry)
     CALL dbcsr_create(matrix_tmp2, template = matrix_z(1)%matrix, &
                       matrix_type = dbcsr_type_no_symmetry)

     ! If closed-shell, factor 2 is introduced twice 
     ! bc of double appearance of density matrix 
     scaling = 1.0_dp
     IF (nspins == 1) scaling = 0.5_dp

     ! Whz = ZFP + PFZ = Z(FP) + (Z(FP))^T
     DO ispin = 1, nspins

        ! tmp = FP
        CALL dbcsr_multiply("N", "N", scaling, matrix_ks(ispin)%matrix, matrix_p(ispin)%matrix, &
                            0.0_dp, matrix_tmp, filter_eps = eps_filter, retain_sparsity=.FALSE.)

        ! tmp2 = ZFP
        CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_z(ispin)%matrix, matrix_tmp, &
                            0.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity=.FALSE.)

        ! tmp = (ZFP)^T
        CALL dbcsr_transposed(matrix_tmp, matrix_tmp2)

        ! tmp = ZFP + (ZFP)^T
        CALL dbcsr_add(matrix_tmp, matrix_tmp2, 1.0_dp, 1.0_dp)

        ! Whz = ZFP + PFZ
        CALL dbcsr_copy(matrix_wz(ispin)%matrix, matrix_tmp, keep_sparsity = .TRUE.)

        CALL dbcsr_filter(matrix_wz(ispin)%matrix, eps_filter)

     ENDDO

     IF(para_env%mepos==0) WRITE(*,*) "Wz = 0.25*(Z*KS*P + P*KS*Z)"
     norm_res = 0.0_dp
     norm_res = dbcsr_frobenius_norm(matrix_wz(1)%matrix)
     IF(para_env%mepos==0) WRITE(*,*) "||Wz||_F", norm_res
     CALL dbcsr_dot(matrix_s(1)%matrix, matrix_wz(1)%matrix, norm_res)
     IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - Tr[WinSin]", norm_res
     CALL dbcsr_dot(matrix_p(1)%matrix, matrix_wz(1)%matrix, norm_res)
     IF(para_env%mepos==0) WRITE(*,*) "ec_res_ao - Tr[WinPin]", norm_res

     CALL dbcsr_release(matrix_tmp)
     CALL dbcsr_release(matrix_tmp2)

     CALL timestop(handle)

  END SUBROUTINE ec_wz_matrix

! ***************************************************************************************************
!> \brief  calculate lin transformation of Hessian matrix on a trial vector (matrix) matrix_cg 
!>         Ax = [F, S]_B + [G(B), S]_Pin
!>         Here, task is calculation of kernel G(B) 
!>         and second commutator
!>          
!>          
!>          
!> \param qs_env ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE build_hessian_op(qs_env, p_env, matrix_ps, matrix_cg, matrix_Ax, eps_filter)

     TYPE(qs_environment_type), POINTER                 :: qs_env
     TYPE(qs_p_env_type), POINTER                       :: p_env
     TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
                                         POINTER        :: matrix_cg, matrix_ps
     TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
                                         POINTER        :: matrix_Ax
     REAL(KIND=dp), INTENT(IN)                          :: eps_filter

     CHARACTER(len=*), PARAMETER :: routineN = 'build_hessian_op', routineP = moduleN//':'//routineN
     
     INTEGER                                            :: handle, ispin, nspins
     TYPE(cp_para_env_type), POINTER                    :: para_env
     TYPE(qs_rho_type), POINTER                         :: rho
     TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_bs, matrix_ks, &
                                                           matrix_s, rho1_ao
     TYPE(dft_control_type), POINTER                    :: dft_control
     REAL(KIND=dp)                                      :: chksum

     CALL timeset(routineN, handle)

     CALL get_qs_env(qs_env = qs_env, &
                     dft_control = dft_control, &
                     matrix_ks = matrix_ks,&
                     matrix_s = matrix_s,&
                     para_env=para_env,&
                     rho = rho)
     nspins = dft_control%nspins
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op BEGINNING", ASSOCIATED(p_env)

     NULLIFY(matrix_bs)
     CALL dbcsr_allocate_matrix_set(matrix_bs, nspins)
     DO ispin = 1, nspins
        ALLOCATE (matrix_bs(ispin)%matrix)
        CALL dbcsr_create(matrix_bs(ispin)%matrix, name = "B*S work matrix", &
                          template = matrix_s(1)%matrix, &
                          matrix_type = dbcsr_type_no_symmetry)
     END DO

     !IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op - calculate B"
     CALL ec_scomm(qs_env, &
                   matrix_a = matrix_cg, &
                   matrix_ps = matrix_ps, &
                   matrix_out = p_env%p1, &
                   eps_filter = eps_filter, &
                   alpha = 1.0_dp, &
                   beta = 0.0_dp,&
                   ps_trans = .TRUE.)

     ! matrix_bs 
     ! Analogous to other S-commutators, we precompute the matrix product B*S
     ! to use same efficient routine to calculate S-commutator
     !IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op - calculate BS"
     DO ispin = 1, nspins
        CALL dbcsr_multiply("N","N", 1.0_dp, p_env%p1(ispin)%matrix, matrix_s(ispin)%matrix, &
                            0.0_dp, matrix_bs(ispin)%matrix, filter_eps = eps_filter)
     END DO

     ! Ax1 = [F,S]_B = FBS - SBF
     !               = F*BS - (F*BS)^T
     !IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op - calculate Ax1"
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
        chksum = chksum+dbcsr_checksum(p_env%p1(ispin)%matrix)
     ENDDO

     IF (para_env%mepos==0) WRITE(*,*) "chksum:", chksum

     ! skip the kernel if the DM is very small
     IF (chksum .GT. 1.0E-14_dp) THEN
        !IF(para_env%mepos==0) WRITE(*,*) "continue with hessian_op2", chksum

        CALL p_env_check_i_alloc(p_env, qs_env)

        CALL qs_rho_get(p_env%rho1, rho_ao=rho1_ao)

        ! Get response density matrix
        DO ispin = 1, nspins
           CALL dbcsr_copy(rho1_ao(ispin)%matrix, p_env%p1(ispin)%matrix)
        ENDDO

        ! Update grids to matrix_b
        ! Same basis as ground-state calculation
        CALL qs_rho_update_rho(p_env%rho1, qs_env)

        ! Calculate kernel
        ! Ax = Ax1           + Ax2
        ! Ax = F*B*S - S*B*F + G(B)*P*S - S*P*G(B)
        CALL hessian_op2(qs_env, matrix_Ax, matrix_ps, p_env%rho1, eps_filter)


     END IF

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
  SUBROUTINE hessian_op2(qs_env, matrix_Ax, matrix_ps, rho1, eps_filter)

     TYPE(qs_environment_type), POINTER                 :: qs_env
     TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_Ax, matrix_ps 
     TYPE(qs_rho_type), POINTER                         :: rho1
     REAL(KIND=dp), INTENT(IN)                          :: eps_filter

     CHARACTER(len=*), PARAMETER :: routineN = 'hessian_op2', routineP = moduleN//':'//routineN
     
     INTEGER                                            :: handle, ispin, nspins
     REAL(KIND=dp)                                      :: alpha, ekin_mol
     TYPE(cp_para_env_type), POINTER                    :: para_env
     TYPE(qs_rho_type), POINTER                         :: rho
     TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho1_ao
     TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_G
     TYPE(section_vals_type), POINTER                   :: input, xc_section
     TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc
     TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
     TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
     TYPE(pw_poisson_type), POINTER                     :: poisson_env
     TYPE(dft_control_type), POINTER                    :: dft_control
     TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                           v_hartree_rspace
     TYPE(pw_env_type), POINTER                         :: pw_env
     TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho1_g

     CALL timeset(routineN, handle)

     NULLIFY (pw_env)

     CALL get_qs_env(qs_env = qs_env, &
                     dft_control = dft_control, &
                     input = input, &
                     matrix_s = matrix_s,&
                     para_env=para_env,&
                     pw_env = pw_env,&
                     rho = rho)
     nspins = dft_control%nspins
     CPASSERT(ASSOCIATED(pw_env))
     !IF(para_env%mepos==0) WRITE(*,*) "hessian_op2 BEGINNING"

     ! Get density matrix
     CALL qs_rho_get(rho, rho_ao = rho_ao)

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

! XC -Kernel
     NULLIFY(v_xc, xc_section)
     xc_section => section_vals_get_subs_vals(input, "DFT%XC")

     ! add xc-kernel
     CALL create_kernel(qs_env, &
                        vxc = v_xc, &
                        rho = rho, &
                        rho1 = rho1, &
                        xc_section = xc_section)
     DO ispin = 1,nspins
        CALL pw_scale(v_xc(ispin)%pw,v_xc(ispin)%pw%pw_grid%dvol)
     END DO

! Ha rtree contribution to Hessian

     ! rho_tot_gspace = rho_g
     ! take B density to build G^{H}[B]
     CALL qs_rho_get(rho1, rho_g = rho1_g)
     CALL pw_zero(rho_tot_gspace%pw)

     DO ispin = 1,nspins
        CALL pw_axpy(rho1_g(ispin)%pw,rho_tot_gspace%pw)
     END DO

     ! get Hartree potential from rho_tot_gspace
     CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, &
                            vhartree = v_hartree_gspace%pw)
     CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
     CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

     ! Add v_xc + v_H
     DO ispin = 1,nspins
        CALL pw_axpy(v_hartree_rspace%pw,v_xc(ispin)%pw)
     END DO
     IF (nspins == 1) CALL pw_scale(v_xc(1)%pw, 2.0_dp)

     NULLIFY (matrix_G)
     CALL dbcsr_allocate_matrix_set(matrix_G, nspins)

     ! integrate this
     DO ispin = 1, nspins
        ! Init response kernel matrix
        ALLOCATE (matrix_G(ispin)%matrix)
        CALL dbcsr_copy(matrix_G(ispin)%matrix, matrix_s(1)%matrix, &
                        name="MATRIX Kernel")
        CALL dbcsr_set(matrix_G(ispin)%matrix, 0.0_dp)

        ! matrix G(B) - the big price
        ! symettric matrix coming out here
        CALL integrate_v_rspace(v_rspace = v_xc(ispin), &
                                pmat = rho_ao(ispin),&
                                hmat = matrix_G(ispin), &
                                qs_env = qs_env, &
                                calculate_forces = .FALSE., &
                                basis_type = "ORB")
     END DO

     ! Calculate KG correction to kernel
     IF (dft_control%qs_control%do_kg) THEN
        IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN

           alpha = 1.0_dp
           IF (nspins == 1) alpha = 2.0_dp
           
           CPASSERT(dft_control%nimages == 1)
           ekin_mol = 0.0_dp
           CALL qs_rho_get(rho1, rho_ao = rho1_ao)
          
           CALL kg_ekin_subset(qs_env = qs_env, &
                               ks_matrix = matrix_G, &
                               ekin_mol = ekin_mol, &
                               calc_force = .FALSE., &
                               do_kernel = .TRUE., &
                               pmat_ext = rho1_ao, &
                               alpha = alpha)
        END IF
     END IF

     ! Hessian already contains
     ! Ax = F*B*S - S*B*F
     ! Now adding 
     ! Ax = Ax1 + G(B)P*S - S*P*G(B)
     !    = Ax1 + G(B)*PS - (G(B)*PS)^T
     !IF(para_env%mepos==0) WRITE(*,*) "hessian_op2 - ec_scomm for Ax2"
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

!     IF(para_env%mepos==0) WRITE(*,*) "hessian_op2 END"
     
     CALL timestop(handle)

  END SUBROUTINE hessian_op2

! ***************************************************************************************************
!> \brief  calculate lin transformation of Hessian matrix on a trial vector (matrix) matrix_cg 
!>         Ax = [F, S]_B + [G(B), S]_Pin
!>         Here, task is calculation of kernel G(B) 
!>         and second commutator
!>          
!>         !!!!!!!!!!!! CALCULATE HESSIAN WITH MATRIX B GIVEN (AKA MATRIX Z)     
!>          
!> \param qs_env ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE build_hessian_op_from_z(qs_env, matrix_ps, matrix_q, matrix_Ax, eps_filter)

     TYPE(qs_environment_type), POINTER                 :: qs_env
     TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
                                         POINTER        :: matrix_ps, matrix_q
     TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
                                         POINTER        :: matrix_Ax
     REAL(KIND=dp), INTENT(IN)                          :: eps_filter

     CHARACTER(len=*), PARAMETER :: routineN = 'build_hessian_op_from_z', routineP = moduleN//':'//routineN
     
     INTEGER                                            :: handle, ispin, nspins
     TYPE(cp_para_env_type), POINTER                    :: para_env
     TYPE(qs_rho_type), POINTER                         :: rho, rho_b
     TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_b, matrix_bs, matrix_ks, &
                                                           matrix_s, matrix_nsc, matrix_rhs,&
                                                           rho_ao_b
     TYPE(dft_control_type), POINTER                    :: dft_control
     REAL(KIND=dp)                                      :: chksum, focc
     TYPE(energy_correction_type), POINTER              :: ec_env
     TYPE(neighbor_list_set_p_type), DIMENSION(:), &
        POINTER                                         :: sab_orb
     REAL(KIND=dp)                                      :: norm_b, norm_ax, norm_res
!     TYPE(dbcsr_type)                                   :: matrix_tmp1
!     REAL(KIND=dp)                                      :: frob_norm

     CALL timeset(routineN, handle)
     CALL get_qs_env(qs_env = qs_env, &
                     dft_control = dft_control, &
                     ec_env = ec_env, &
                     matrix_ks = matrix_ks,&
                     matrix_s = matrix_s,&
                     para_env=para_env,&
                     rho = rho,&
                     sab_orb = sab_orb)
     nspins = dft_control%nspins
     IF(para_env%mepos==0) WRITE(*,*) repeat("-",80)
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z BEGINNING"
     IF(para_env%mepos==0) WRITE(*,*) repeat("-",80)

     norm_res  = dbcsr_frobenius_norm(ec_env%p_env%p1(1)%matrix)
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z - ||Z(MO)||_F", norm_res

     NULLIFY(matrix_b)
     CALL dbcsr_allocate_matrix_set(matrix_b, nspins)
     DO ispin = 1, nspins
        ALLOCATE (matrix_b(ispin)%matrix)
        CALL dbcsr_create(matrix_b(ispin)%matrix, name = "[X,S]_P RSP DNSTY", &
                          template = matrix_s(1)%matrix)
        CALL cp_dbcsr_alloc_block_from_nbl(matrix_b(ispin)%matrix, sab_orb)
        CALL dbcsr_copy(matrix_b(ispin)%matrix, ec_env%p_env%p1(ispin)%matrix)

        norm_res  = dbcsr_frobenius_norm(matrix_b(1)%matrix)
        IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z - ||B||_F", norm_res

        ! Z matrix stored in ec_env was scaled with 0.5 in closed shell case
        CALL dbcsr_scale(matrix_b(ispin)%matrix, 2.0_dp)
     END DO

     NULLIFY(matrix_bs)
     CALL dbcsr_allocate_matrix_set(matrix_bs, nspins)
     DO ispin = 1, nspins
        ALLOCATE (matrix_bs(ispin)%matrix)
        CALL dbcsr_create(matrix_bs(ispin)%matrix, name = "B*S work matrix", &
                          template = matrix_ks(1)%matrix, &
                          matrix_type = dbcsr_type_no_symmetry)
     END DO
     NULLIFY(matrix_rhs)
     CALL dbcsr_allocate_matrix_set(matrix_rhs, nspins)
     DO ispin = 1, nspins
        ALLOCATE (matrix_rhs(ispin)%matrix)
        CALL dbcsr_create(matrix_rhs(ispin)%matrix, name = "RHS", &
                          template = matrix_s(1)%matrix, &
                          matrix_type = dbcsr_type_no_symmetry)
     END DO
     NULLIFY(matrix_nsc)
     CALL dbcsr_allocate_matrix_set(matrix_nsc, nspins)
     DO ispin = 1, nspins
        ALLOCATE (matrix_nsc(ispin)%matrix)
        CALL dbcsr_create(matrix_nsc(ispin)%matrix, name = "RHS", &
                          template = matrix_s(1)%matrix, &
                          matrix_type = dbcsr_type_no_symmetry)

        CALL dbcsr_copy(matrix_nsc(ispin)%matrix, ec_env%matrix_hz(ispin)%matrix)
     END DO

     ! Spin factor
     focc = -2.0_dp
     IF (nspins == 1) focc = -4.0_dp
     ! E^[1]_Harris = G(\delta P)P*S - S*P*G(\delta P)
     ! ABC-CBA
     ! RES ANTISYM
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z - calculate RHS"
     CALL ec_scomm(qs_env = qs_env, &
                   !matrix_a = matrix_nsc, &
                   matrix_a = ec_env%matrix_hz, &
                   matrix_ps = matrix_ps, &
                   matrix_out = matrix_rhs, &
                   eps_filter = eps_filter, &
                   alpha = -4.0_dp, &
                   beta = 0.0_dp, &
                   ps_trans = .FALSE.)

     norm_res  = dbcsr_frobenius_norm(matrix_rhs(1)%matrix)
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z - ||RHS||_F", norm_res

     ! matrix_bs 
     ! Analogous to other S-commutators, we precompute the matrix product B*S
     ! to use same efficient routine to calculate S-commutator
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z - calculate BS"
     DO ispin = 1, nspins
        CALL dbcsr_multiply("N","N", 1.0_dp, matrix_b(ispin)%matrix, matrix_s(ispin)%matrix, &
                            0.0_dp, matrix_bs(ispin)%matrix, filter_eps = eps_filter, &
                            retain_sparsity = .FALSE.)
     END DO

     norm_res  = dbcsr_frobenius_norm(matrix_bs(1)%matrix)
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z - ||BS||_F", norm_res

     ! Ax1 = [F,S]_B = FBS - SBF
     !               = F*BS - (F*BS)^T
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z - calculate Az1"
     CALL ec_scomm(qs_env, &
                   matrix_a = matrix_ks, &
                   matrix_ps = matrix_bs, &
                   matrix_out = matrix_Ax, &
                   eps_filter = eps_filter, &
                   alpha = 1.0_dp, &
                   beta = 0.0_dp,&
                   ps_trans = .FALSE.)

     norm_res  = dbcsr_frobenius_norm(matrix_Ax(1)%matrix)
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z - ||Az1||_F", norm_res

     ! Check 2nd part
     chksum = 0.0_dp
     DO ispin = 1, nspins
        chksum = chksum+dbcsr_checksum(ec_env%p_env%p1(ispin)%matrix)
     ENDDO

     IF (para_env%mepos==0) WRITE(*,*) "chksum:", chksum

     ! skip the kernel if the DM is very small
     IF (chksum .GT. 1.0E-14_dp) THEN
        !IF(para_env%mepos==0) WRITE(*,*) "continue with hessian_op2", chksum

        ! Bring B as density on grid
        NULLIFY(rho_b)
        CALL qs_rho_create(rho_b)
        CALL duplicate_rho_type(rho_input = rho, &
                                rho_output = rho_b, &
                                qs_env = qs_env)

        ! Get response density matrix
        CALL qs_rho_get(rho_b, rho_ao = rho_ao_b)

        ! Copy b-density matrix to its rho environment for collocation
        ! going from non-symmetric to symmetric
        DO ispin = 1, nspins
           CALL dbcsr_copy(rho_ao_b(ispin)%matrix, matrix_b(ispin)%matrix, keep_sparsity =.TRUE.)
        END DO

        ! update rho_r_b and rho_g_b wrt rho_ao_b
        ! Same basis as reference calculation
        CALL qs_rho_update_rho(rho_b, qs_env)

        !IF (para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z - Ax2" 
        ! Calculate kernel
        ! Ax = Ax1           + Ax2
        ! Ax = F*B*S - S*B*F + G(B)*P*S - S*P*G(B)
        !                                IN/OUT     IN         IN
        CALL hessian_op2(qs_env, matrix_Ax, matrix_ps, rho_b, eps_filter)

        ! Test Z Matrix
!        WRITE(*,*) "build_hessian_op_from_z - projector on MATRIX Az"
!        CALL projector_test(qs_env, matrix_Ax, matrix_p, matrix_ps, eps_filter)
!        WRITE(*,*) "build_hessian_op_from_z - Proj(Ax1+2)"
!        CALL projector_test(qs_env, matrix_Ax, matrix_p, matrix_ps, eps_filter)

        CALL qs_rho_release(rho_b)
     END IF

     !--------------

     ! Spin factor
     focc = -2.0_dp
     IF (nspins == 1) focc = -4.0_dp
     ! E^[1]_Harris = G(\delta P)P*S - S*P*G(\delta P)
     ! ABC-CBA
     ! RES ANTISYM
     DO ispin = 1, nspins
        CALL dbcsr_set(matrix_rhs(ispin)%matrix, 0.0_dp)
     END DO
     CALL ec_scomm(qs_env = qs_env, &
                   matrix_a = matrix_nsc, &
                   matrix_ps = matrix_ps, &
                   matrix_out = matrix_rhs, &
                   eps_filter = eps_filter, &
                   alpha = -4.0_dp, &
                   beta = 0.0_dp, &
                   ps_trans = .FALSE.)

     ! Ax1+2 - b = 0
     norm_res = 0.0_dp 

     IF (para_env%mepos==0) WRITE(*,*) "BEFORE PROJECTED :"
     norm_b  = dbcsr_frobenius_norm(matrix_rhs(1)%matrix)
     norm_ax = dbcsr_frobenius_norm(matrix_Ax(1)%matrix)
     IF (para_env%mepos==0) WRITE(*,*) "|| b_ao  ||_F :", norm_b
     IF (para_env%mepos==0) WRITE(*,*) "|| Az_ao ||_F :", norm_ax
     norm_res = norm_b - norm_ax
     IF (para_env%mepos==0) WRITE(*,*) "||b_ao||_F - ||Ax_ao||_F :", norm_res

     IF (para_env%mepos==0) WRITE(*,*) "PROJECTED :"
     CALL projector(qs_env, matrix_ps, matrix_q, matrix_rhs, eps_filter, .TRUE.)
     CALL projector(qs_env, matrix_ps, matrix_q, matrix_Ax, eps_filter, .TRUE.)

     norm_b  = dbcsr_frobenius_norm(matrix_rhs(1)%matrix)
     norm_ax = dbcsr_frobenius_norm(matrix_Ax(1)%matrix)
     IF (para_env%mepos==0) WRITE(*,*) "|| b_ao  ||_F :", norm_b
     IF (para_env%mepos==0) WRITE(*,*) "|| Az_ao ||_F :", norm_ax
     norm_res = norm_b - norm_ax
     IF (para_env%mepos==0) WRITE(*,*) "||P(b_ao)||_F - ||P(Ax_ao)||_F :", norm_res

     CALL dbcsr_add(matrix_rhs(1)%matrix, matrix_Ax(1)%matrix, 1.0_dp, -1.0_dp)
     norm_res = dbcsr_frobenius_norm(matrix_rhs(1)%matrix)
     IF (para_env%mepos==0) WRITE(*,*) "no projection:"
     IF (para_env%mepos==0) WRITE(*,*) "|| b_ao - Az_ao ||_F :", norm_res
     CALL dbcsr_dot(matrix_rhs(1)%matrix, matrix_rhs(1)%matrix, norm_res)
     IF (para_env%mepos==0) WRITE(*,*) "tr(r*r) :", norm_res
     
     IF (para_env%mepos==0) WRITE(*,*) "Apply projection P(r):"
     CALL projector(qs_env, matrix_ps, matrix_q, matrix_rhs, eps_filter, .TRUE.)
     IF (para_env%mepos==0) WRITE(*,*) "after projection:"
     CALL dbcsr_dot(matrix_rhs(1)%matrix, matrix_rhs(1)%matrix, norm_res)
     IF (para_env%mepos==0) WRITE(*,*) "tr[r*r] :", norm_res
     !--------------

     CALL dbcsr_deallocate_matrix_set(matrix_b)
     CALL dbcsr_deallocate_matrix_set(matrix_bs)
     CALL dbcsr_deallocate_matrix_set(matrix_nsc)
     CALL dbcsr_deallocate_matrix_set(matrix_rhs)

     IF(para_env%mepos==0) WRITE(*,*) repeat("-",80)
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z END"
     IF(para_env%mepos==0) WRITE(*,*) repeat("-",80)

     CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief Projector P(M) = P*M*Q^T + Q*M*P^T, or
!>                  P^T(M) = P^T*M*Q + Q^T*M*P
!>        with P = D*S, Q = (1-D*S), and D=CC^T
!> 
!> \param qs_env ...
!> \param matrix_ps Projector P
!> \param matrix_q Projector Q 
!> \param matrix_io Matrix M to which projector is applied.
!> \param eps_filter ...
!> \param ptrans Decides if P(M) (.FALSE. or P^T(M) (.TRUE.) is applied to M
!>          
!> \date    06.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
!  SUBROUTINE projector(qs_env, matrix_ps, matrix_q, matrix_io, eps_filter, ptrans)
!
!     TYPE(qs_environment_type), POINTER                 :: qs_env
!     TYPE(dbcsr_p_type), DIMENSION(:), &
!                             INTENT(IN), POINTER        :: matrix_ps, matrix_q
!     TYPE(dbcsr_p_type), DIMENSION(:), &
!                          INTENT(INOUT), POINTER        :: matrix_io
!     REAL(KIND=dp), INTENT(IN)                          :: eps_filter
!     LOGICAL, INTENT(IN), OPTIONAL                      :: ptrans
!
!     CHARACTER(len=*), PARAMETER :: routineN = 'projector', routineP = moduleN//':'//routineN
!     
!     INTEGER                                            :: handle, ispin, nspins
!     LOGICAL                                            :: my_ptrans
!     TYPE(cp_para_env_type), POINTER                    :: para_env
!     TYPE(dft_control_type), POINTER                    :: dft_control
!     TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_tmp2
!
!     CALL timeset(routineN, handle)
!
!     CALL get_qs_env(qs_env = qs_env, &
!                     dft_control = dft_control, &
!                     para_env=para_env)
!     nspins = dft_control%nspins
!!     IF(para_env%mepos==0) WRITE(*,*) "projector BEGINNING"
!
!     ! Decides which projector to apply 
!     ! FALSE : Projector onto space of non-redudnant elements
!     ! TRUE : It's transpose 
!     my_ptrans = .FALSE.
!     IF (PRESENT(ptrans)) my_ptrans = ptrans
!  
!     CALL dbcsr_create(matrix_tmp1, template = matrix_ps(1)%matrix, &
!                       matrix_type = dbcsr_type_no_symmetry)
!     CALL dbcsr_create(matrix_tmp2, template = matrix_ps(1)%matrix, &
!                       matrix_type = dbcsr_type_no_symmetry)
!
!     IF (my_ptrans) THEN
!        ! Proj^T(M) = P^T*M*Q + Q^T*M*P
!        ! with P = D*S, Q = (1 - P), and D = CC^T
!        DO ispin = 1, nspins 
!           ! tmp1 = P^T*M
!           CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_ps(ispin)%matrix, matrix_io(ispin)%matrix, &
!                               0.0_dp, matrix_tmp1, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!           ! tmp2 = P^T*M*Q
!           CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, matrix_q(ispin)%matrix, &
!                               0.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!           ! tmp1 = Q^T*M
!           CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_q(ispin)%matrix, matrix_io(ispin)%matrix, &
!                               0.0_dp, matrix_tmp1, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!           ! tmp2 = tmp2 + Q^T*M*P
!           CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, matrix_ps(ispin)%matrix, &
!                               1.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!        END DO
!
!     ELSE
!
!        ! Proj(M) = P*M*Q^T + Q*M*P^T
!        ! with P = D*S, Q = (1 - P), and D = CC^T
!        DO ispin = 1, nspins 
!           ! tmp1 = P*M
!           CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_ps(ispin)%matrix, matrix_io(1)%matrix, &
!                               0.0_dp, matrix_tmp1, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!           ! tmp2 = P*M*Q^T
!           CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp1, matrix_q(ispin)%matrix, &
!                               0.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!           ! tmp1  = Q*M
!           CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_q(ispin)%matrix, matrix_io(ispin)%matrix, &
!                               0.0_dp, matrix_tmp1, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!           ! tmp2 = tmp2 + Q*M*P^T
!           CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp1, matrix_ps(ispin)%matrix, &
!                               1.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!        END DO
!
!     END IF
!
!     ! Export
!     DO ispin = 1, nspins 
!        !CALL dbcsr_copy(matrix_io(ispin)%matrix, matrix_tmp2, keep_sparsity = .FALSE.)
!     END DO
!
!     CALL dbcsr_release(matrix_tmp1)
!     CALL dbcsr_release(matrix_tmp2)
!
!     CALL timestop(handle)
!
!  END SUBROUTINE 


! **************************************************************************************************
!> \brief Projector P(M) = P*M*Q^T + Q*M*P^T, or
!>                  P^T(M) = P^T*M*Q + Q^T*M*P
!>        with P = D*S, Q = (1-D*S), and D=CC^T
!> 
!> \param qs_env ...
!> \param matrix_ps Projector P
!> \param matrix_q Projector Q 
!> \param matrix_io Matrix M to which projector is applied.
!> \param eps_filter ...
!> \param ptrans Decides if P(M) (.FALSE.i) or P^T(M) (.TRUE.) is applied to M
!>          
!> \date    06.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE projector(qs_env, matrix_ps, matrix_q, matrix_io, eps_filter, ptrans)

     TYPE(qs_environment_type), POINTER                 :: qs_env
     TYPE(dbcsr_p_type), DIMENSION(:), &
                             INTENT(IN), POINTER        :: matrix_ps, matrix_q
     TYPE(dbcsr_p_type), DIMENSION(:), &
                          INTENT(INOUT), POINTER        :: matrix_io
     REAL(KIND=dp), INTENT(IN)                          :: eps_filter
     LOGICAL, INTENT(IN), OPTIONAL                      :: ptrans

     CHARACTER(len=*), PARAMETER :: routineN = 'projector', routineP = moduleN//':'//routineN
     
     INTEGER                                            :: handle, ispin, nspins
     LOGICAL                                            :: my_ptrans
     TYPE(cp_para_env_type), POINTER                    :: para_env
     TYPE(dft_control_type), POINTER                    :: dft_control
     TYPE(dbcsr_type)                                   :: matrix_tmp

     CALL timeset(routineN, handle)

     CALL get_qs_env(qs_env = qs_env, &
                     dft_control = dft_control, &
                     para_env=para_env)
     nspins = dft_control%nspins
     !IF(para_env%mepos==0) WRITE(*,*) "projector BEGINNING"

     ! Decides which projector to apply 
     ! FALSE : Projector onto space of non-redudnant elements
     ! TRUE : It's transpose 
     my_ptrans = .FALSE.
     IF (PRESENT(ptrans)) my_ptrans = ptrans
  
     CALL dbcsr_create(matrix_tmp, template = matrix_ps(1)%matrix, &
                       matrix_type = dbcsr_type_no_symmetry)

     DO ispin = 1, nspins 
        ! Proj^T(M) = P^T*M*Q + Q^T*M*P
        ! with P = D*S, Q = (1 - P), and D = CC^T
        IF (my_ptrans) THEN
           ! tmp1 = P^T*M
           CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_ps(ispin)%matrix, matrix_io(ispin)%matrix, &
                               0.0_dp, matrix_tmp, filter_eps = eps_filter) 
           ! tmp2 = P^T*M*Q
           CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, matrix_q(ispin)%matrix, &
                               0.0_dp, matrix_io(ispin)%matrix, filter_eps = eps_filter) 

           ! tmp1 = (P^T*M^T*Q)^T = -(P^T*M*Q)^T

        ELSE

           ! Proj(M) = P*M*Q^T + Q*M*P^T
           ! with P = D*S, Q = (1 - P), and D = CC^T

           ! tmp1 = P*M
           CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_ps(ispin)%matrix, matrix_io(1)%matrix, &
                               0.0_dp, matrix_tmp, filter_eps = eps_filter) 
           ! tmp2 = P*M*Q^T
           CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, matrix_q(ispin)%matrix, &
                               0.0_dp, matrix_io(ispin)%matrix, filter_eps = eps_filter) 
!           ! tmp1  = Q*M
!           CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_q(ispin)%matrix, matrix_io(ispin)%matrix, &
!                               0.0_dp, matrix_tmp1, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!           ! tmp2 = tmp2 + Q*M*P^T
!           CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp1, matrix_ps(ispin)%matrix, &
!                               1.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity = .FALSE.) 

           ! tmp1 = (P*M^T*Q^T)^T = -(P*M*Q^T)^T

        END IF

        CALL dbcsr_transposed(matrix_tmp, matrix_io(ispin)%matrix)
        CALL dbcsr_add(matrix_io(ispin)%matrix, matrix_tmp, 1.0_dp, -1.0_dp)

     END DO

     CALL dbcsr_release(matrix_tmp)

     CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief Projector Q = (1 - DS), with D = CC^T 
!> 
!> \param qs_env ...
!> \param matrix_ps Matrix D*S (D=CC^T), also projector on the occupied, unperturbed ground-state 
!> \param matrix_q Projector Q on the unoccupied, unperturbed ground-state 
!> \param eps_filter ...
!>          
!> \date    06.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE projector_p_and_q(qs_env, matrix_ps, matrix_q, eps_filter)

     TYPE(qs_environment_type), POINTER                 :: qs_env
     TYPE(dbcsr_p_type), DIMENSION(:), &
                             INTENT(INOUT), POINTER     :: matrix_ps, matrix_q
     REAL(KIND=dp), INTENT(IN)                          :: eps_filter

     CHARACTER(len=*), PARAMETER :: routineN = 'projector_p_and_q', routineP = moduleN//':'//routineN
     
     INTEGER                                            :: handle, ispin, nspins
     REAL(KIND=dp)                                      :: focc
     TYPE(cp_para_env_type), POINTER                    :: para_env
     TYPE(dft_control_type), POINTER                    :: dft_control
     TYPE(qs_rho_type), POINTER                         :: rho
     TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao

     CALL timeset(routineN, handle)

     CALL get_qs_env(qs_env = qs_env, &
                     dft_control = dft_control, &
                     matrix_s = matrix_s,&
                     para_env=para_env,&
                     rho = rho)
     nspins = dft_control%nspins
     CALL qs_rho_get(rho, rho_ao = rho_ao)
!     IF(para_env%mepos==0) WRITE(*,*) "projector_p_and_q BEGINNING"
 
     ! Allocate projectors P and Q
     CALL dbcsr_allocate_matrix_set(matrix_ps, nspins)
     CALL dbcsr_allocate_matrix_set(matrix_q, nspins)
     DO ispin = 1, nspins
        ALLOCATE (matrix_ps(ispin)%matrix)
        ALLOCATE (matrix_q(ispin)%matrix)
        CALL dbcsr_create(matrix_ps(ispin)%matrix, name = "MATRIX PS", &
                          template = matrix_s(1)%matrix,&
                          matrix_type=dbcsr_type_no_symmetry)
        CALL dbcsr_create(matrix_q(ispin)%matrix, name = "PROJECTOR Q", &
                          template = matrix_s(1)%matrix,&
                          matrix_type=dbcsr_type_no_symmetry)
     END DO

     ! Double-occupancy in closed-shell calculations: Scale rho_ao by 1/2
     focc = 1.0_dp
     IF (nspins == 1) focc = 0.5_dp

     ! matrix_ps
     ! Precomupte matrix product D*S for S-Commutator 
     ! This is also the projector onto the occupied, unperturbed ground-state
     DO ispin = 1, nspins
        CALL dbcsr_multiply("N","N", focc, rho_ao(ispin)%matrix, matrix_s(ispin)%matrix, &
                            0.0_dp, matrix_ps(ispin)%matrix, filter_eps = eps_filter)
     END DO
        
     ! Q = (1 - PS)
     ! projector onto the un-occupied, unperturbed ground-state
     DO ispin = 1, nspins
        CALL dbcsr_copy(matrix_q(ispin)%matrix, matrix_ps(ispin)%matrix, name = "PROJECTOR Q")
        CALL dbcsr_scale(matrix_q(ispin)%matrix, -1.0_dp)
        CALL dbcsr_add_on_diag(matrix_q(ispin)%matrix, 1.0_dp)
     END DO

     CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief S-Commutator
!>        matrix_out = alpha*[A,C]_B + beta*matrix_out
!>                   = alpha*(ABC - CBA) + beta*matrix_out
!>
!> \param qs_env ... 
!> \param matrix_a Input matrix 
!> \param matrix_ps Precomputed matrix product D*S (D=CC^T)
!> \param matrix_out S-commutator result matrix_out = beta*matrix_out + alpha*[A,C]_B
!>                                                  = beta*matrix_out + alpha*(ABC - CBA)
!> \param eps_filter ...
!> \param alpha Multiplicator for s-commutator 
!> \param beta Multiplicator for input content of matrix_out 
!> \param ps_trans Use tranpose of PS matrix (matrix_out has S-sparsity)
!> \par History
!>       2020.3 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_scomm(qs_env, matrix_a, matrix_ps, matrix_out, eps_filter, alpha, beta, ps_trans)

       TYPE(qs_environment_type), POINTER                 :: qs_env
       TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
          POINTER                                         :: matrix_a, matrix_ps
       TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
          POINTER                                         :: matrix_out
       REAL(KIND=dp), INTENT(IN)                          :: eps_filter
       REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha, beta
       LOGICAL, INTENT(IN), OPTIONAL                      :: ps_trans

       CHARACTER(len=*), PARAMETER :: routineN = 'ec_scomm', &
          routineP = moduleN//':'//routineN

       INTEGER                                            :: handle, ispin, nspins
       TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_tmp2
       TYPE(dft_control_type), POINTER                    :: dft_control
       REAL(KIND=dp)                                      :: myalpha, mybeta
       TYPE(cp_para_env_type), POINTER                    :: para_env
       LOGICAL                                            :: my_pst

       CALL timeset(routineN, handle)

       CALL get_qs_env(qs_env = qs_env, &
                       dft_control = dft_control, &
                       para_env = para_env)
       nspins = dft_control%nspins
       !IF (para_env%mepos == 0) WRITE(*,*) "ec_scomm - BEGINNING"

       my_pst = .FALSE.
       IF (PRESENT(ps_trans)) my_pst = ps_trans
       myalpha = 1.0_dp
       IF (PRESENT(alpha)) myalpha = alpha 
       mybeta = 0.0_dp
       IF (PRESENT(beta)) mybeta = beta

       ! init temp matrices
       CALL dbcsr_create(matrix_tmp1, template=matrix_a(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_tmp2, template=matrix_a(1)%matrix, matrix_type=dbcsr_type_no_symmetry)

       ! ABC - CBA 
       ! Two cases : 
       ! Matrix A symmetric
       !    RHS : G(dp)*P*S - S*P*G(dP)  
       !    Ax1 : F*B*S - S*B*F         !! needs PS to be BS
       !    Ax2 : G(B)*P*S - S*P*G(B)
       ! 
       ! Matrix A skew-symmetric 
       !    B   : cg*S*P - P*S*cg       !! needs PS transposed
       !    Z   : cg_z*S*P - P*S*cg     !! needs PS transposed

       DO ispin = 1, nspins
          IF(.NOT. my_pst) THEN

             ! D = [A,C]_B = ABC - CBA = ABC - (ABC)^{T}
             ! tmp1 = ABC = < matrix_a > * PS
             CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_a(ispin)%matrix, matrix_ps(ispin)%matrix, &
                                 0.0_dp, matrix_tmp1, filter_eps = eps_filter)

             ! tmp_2 = (tmp1)^T
             CALL dbcsr_transposed(matrix_tmp2, matrix_tmp1)

             ! tmp1 = tmp1 - (tmp1)^T
             CALL dbcsr_add(matrix_tmp1, matrix_tmp2, 1.0_dp, -1.0_dp)

             ! D = beta*D + alpha*(ABC - (ABC)^T)
             ! Outputs a skew-symmetric matrix
             CALL dbcsr_add(matrix_out(ispin)%matrix, matrix_tmp1, mybeta, myalpha)

          ELSE

             ! Case for B and Z matrix - i.e. matrix_a is ANTI-SYM 
             ! [cg, P]_S = cg*S*P - P*S*cg
             ! tmp1 = A*t
             ! tmp1 = cg*(PS)^T = cg*S*P = A*B*C
             CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_a(ispin)%matrix, matrix_ps(ispin)%matrix, &
                                 0.0_dp, matrix_tmp1, filter_eps = eps_filter)
             ! tmp_2 = (tmp1)^T
             CALL dbcsr_transposed(matrix_tmp2, matrix_tmp1)
             ! tmp1 = alpha*(tmp1 + (tmp1)^T)
             CALL dbcsr_add(matrix_tmp1, matrix_tmp2, myalpha, myalpha)

             ! D = ABC + (ABC)^T
             ! Outputs a symmetric matrix
             CALL dbcsr_copy(matrix_out(ispin)%matrix, matrix_tmp1, keep_sparsity = .TRUE.)

          END IF

          CALL dbcsr_filter(matrix_out(ispin)%matrix, eps_filter)

       ENDDO

       CALL dbcsr_release(matrix_tmp1)
       CALL dbcsr_release(matrix_tmp2)

       CALL timestop(handle)

   END SUBROUTINE ec_scomm

! **************************************************************************************************
!> \brief Win = 0.5(P*Win*P)
!>
!>
!> \param
!> \par History
!>       2020.2 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE calculate_2whz_ao_matrix(qs_env, matrix_hz, matrix_whz, focc, eps_filter)

       TYPE(qs_environment_type), POINTER                  :: qs_env
       TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
                                            POINTER        :: matrix_hz
       TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
                                            POINTER        :: matrix_whz
       REAL(KIND=dp), INTENT(IN)                           :: focc, eps_filter

       CHARACTER(len=*), PARAMETER :: routineN = 'calculate_2whz_ao_matrix', &
          routineP = moduleN//':'//routineN

       INTEGER                                            :: handle, ispin, nspins
       REAL(KIND=dp)                                      :: falpha, scaling
       TYPE(dbcsr_type)                                   :: matrix_tmp
       TYPE(dft_control_type), POINTER                    :: dft_control
       TYPE(cp_para_env_type), POINTER                    :: para_env
       TYPE(qs_rho_type), POINTER                         :: rho
       TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao
       !REAL(KIND=dp)                                      :: trace

       CALL timeset(routineN, handle)

       CALL get_qs_env(qs_env = qs_env, &
                       dft_control = dft_control, &
                       rho = rho,&
                       para_env=para_env)
       nspins = dft_control%nspins
       CALL qs_rho_get(rho, rho_ao = rho_ao)

       !IF (para_env%mepos==0) WRITE(*,*) "calculate_2whz_ao_matrix - BEGINNING"

       ! init temp matrix
       CALL dbcsr_create(matrix_tmp, template=matrix_hz(1)%matrix, matrix_type=dbcsr_type_no_symmetry)

       ! Operation in MO-solver closed-shell: 
       ! Whz = 0.5*(CC^T*Hz*CC^T)
       ! Here,  P = 2*CC^T
       ! Liear sclaing implies absence of explicit knowledge of molecular orbital coefficients
       !!! Scaling the density matrix at this stage, means dividing the final result by 4
       !!! as two multiplications with P are carried out.
       scaling = 1.0_dp
       IF (nspins == 1) scaling = 0.5_dp

       ! Scaling for closed shell-case on Hz
       falpha = focc*scaling

       DO ispin = 1, nspins

          ! tmp1 = H*CC^T
          CALL dbcsr_multiply("N", "N", scaling, matrix_hz(ispin)%matrix, rho_ao(ispin)%matrix, &
                              0.0_dp, matrix_tmp, filter_eps = eps_filter)
          ! tmp2 = CC^T*tmp1
          ! Whz = Wz + focc*(scaling*P)*Hz*(scaling*P)
          CALL dbcsr_multiply("N", "N", falpha, rho_ao(ispin)%matrix, matrix_tmp, &
                              1.0_dp, matrix_whz(ispin)%matrix, filter_eps = eps_filter, &
                              retain_sparsity = .TRUE.)

       END DO

!       IF(para_env%mepos==0) WRITE(*,*) "Whz = CC^T*Hz*CC^T AO"
!       trace = 0.0_dp
!       trace = dbcsr_frobenius_norm(matrix_whz(1)%matrix)
!       IF (para_env%mepos==0) WRITE(*,*) "calculate_2whz_ao_matrix - ||Whz|||_f", trace

       CALL dbcsr_release(matrix_tmp)

       CALL timestop(handle)

   END SUBROUTINE

END MODULE ec_ao_solver
