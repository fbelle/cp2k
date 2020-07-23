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
MODULE ec_orth_solver
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_external_control,             ONLY: external_control
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
                                              dbcsr_finalize,&
                                              dbcsr_frobenius_norm,&
                                              dbcsr_get_info,&
                                              dbcsr_multiply,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_print,&
                                              dbcsr_scale,&
                                              dbcsr_release,&
                                              dbcsr_trace,&
                                              dbcsr_transposed,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric,&
                                              dbcsr_type_no_symmetry
   USE ec_env_types,                    ONLY: energy_correction_type 
   USE input_constants,                 ONLY: ls_s_sqrt_proot,&
                                              ls_s_sqrt_ns,&
                                              ot_precond_full_single_inverse
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kg_correction,                   ONLY: create_kernel
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_methods,               ONLY: p_env_check_i_alloc
   USE qs_linres_types,                 ONLY: linres_control_create,&
                                              linres_control_release,&
                                              linres_control_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   ! borrowed from curvy steps
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
   USE qs_p_env_methods,                ONLY: p_env_create,&
                                              p_env_psi0_changed
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   ! fbelle
   USE ec_ao_solver,                   ONLY: ec_wz_matrix
      
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_response_solver'

! *** Public subroutines ***

   PUBLIC :: ec_response_ao

CONTAINS

! **************************************************************************************************
!> \brief   
!>          
!>          
!>          
!> \param qs_env ...
!> \date    01.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_response_ao(qs_env, ec_env, p_env, iounit, should_stop)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(energy_correction_type), POINTER              :: ec_env
    TYPE(qs_p_env_type), POINTER                       :: p_env
    INTEGER, INTENT(IN)                                :: iounit
    LOGICAL, INTENT(OUT)                               :: should_stop

    CHARACTER(len=*), PARAMETER :: routineN = 'ec_response_ao', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: i, handle, ispin, max_iter_lanczos, nao, &
                                                          nspins, s_sqrt_method, s_sqrt_order
    LOGICAL                                            :: restart 
    REAL(KIND=dp)                                      :: eps_filter, eps_lanczos, focc, norm_res, &
                                                          t1, t2
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ksmat, matrix_Ax, matrix_cg, matrix_cg_z, &
                                                          matrix_ks, matrix_nsc, matrix_p, matrix_res, &
                                                          matrix_s, matrix_z, matrix_z0, &
                                                          rho_ao
    TYPE(dbcsr_type)                                   :: matrix_s_sqrt, matrix_s_sqrt_inv
    TYPE(qs_rho_type), POINTER                         :: rho
    REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha, beta, new_norm, norm_cA, &
                                                          norm_rr, tr_rz00
    REAL(KIND=dp)                                      :: min_shift, old_conv, shift
    TYPE(dbcsr_type)                                   :: matrix_tmp
    LOGICAL                                            :: do_precond
    TYPE(linres_control_type), POINTER                 :: linres_control
    TYPE(section_vals_type), POINTER                   :: solver_section
    ! line search stuff
    REAL(KIND=dp)                                      :: astep, devi(3), fac, fac1, lin_eq(3, 3), vec(3)
    INTEGER                                            :: info, iter
    REAL(KIND=dp)                                      :: norm_zao, norm_zmo

    CALL timeset(routineN, handle)

    NULLIFY(dft_control, ksmat, matrix_s, linres_control, para_env, rho)

    t1 = m_walltime()

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    linres_control=linres_control, &
                    matrix_ks = ksmat, &
                    matrix_s = matrix_s, &
                    para_env=para_env,&
                    rho = rho)
    nspins = dft_control%nspins
    IF(para_env%mepos==0) WRITE(*,*) "ec_orth_solver - BEGINNING"

    CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=nao)

!    do_precond = .TRUE.
    do_precond = .FALSE.

    !iounit = 6
!    s_sqrt_method = 2
!    s_sqrt_order = 3
!    eps_lanczos = 1.0E-3_dp
!    max_iter_lanczos = 128

    solver_section => section_vals_get_subs_vals(qs_env%input, "DFT%ENERGY_CORRECTION%LS_SOLVER")
    CALL section_vals_val_get(solver_section, "S_SQRT_METHOD", i_val=s_sqrt_method)
    CALL section_vals_val_get(solver_section, "S_SQRT_ORDER", i_val=s_sqrt_order)
    CALL section_vals_val_get(solver_section, "EPS_LANCZOS", r_val=eps_lanczos)
    CALL section_vals_val_get(solver_section, "MAX_ITER_LANCZOS", i_val=max_iter_lanczos)  

    eps_filter = linres_control%eps_filter
    IF(para_env%mepos==0) WRITE(*,*) "ec_orth_solver - eps", linres_control%eps
    IF(para_env%mepos==0) WRITE(*,*) "ec_orth_solver - eps_filter", eps_filter

    IF(para_env%mepos==0) WRITE(*,*) "ec_orth_solver BEGINNING"

    CALL qs_rho_get(rho, rho_ao = rho_ao)

    ALLOCATE(alpha(nspins),beta(nspins), new_norm(nspins), norm_cA(nspins), norm_rr(nspins))
    ALLOCATE(tr_rz00(nspins))

    ! local matrix P, KS, and NSC
    ! to bring into orthogonal basis
    NULLIFY(matrix_p, matrix_ks, matrix_nsc)
    CALL dbcsr_allocate_matrix_set(matrix_p, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_ks, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_nsc, nspins)
    DO ispin = 1, nspins
       ALLOCATE (matrix_p(ispin)%matrix)
       ALLOCATE (matrix_ks(ispin)%matrix)
       ALLOCATE (matrix_nsc(ispin)%matrix)
       CALL dbcsr_create(matrix_p(ispin)%matrix, name = "P_IN ORTHO", &
                         template = ksmat(1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_ks(ispin)%matrix, name = "KS_IN ORTHO", &
                         template = ksmat(1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_nsc(ispin)%matrix, name = "NSC IN ORTHO", &
                         template = ksmat(1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
    
       CALL dbcsr_desymmetrize(rho_ao(ispin)%matrix, matrix_p(ispin)%matrix)
       CALL dbcsr_desymmetrize(ksmat(ispin)%matrix, matrix_ks(ispin)%matrix)
       CALL dbcsr_desymmetrize(ec_env%matrix_hz(ispin)%matrix, matrix_nsc(ispin)%matrix)
    END DO

    ! Scale matrix_p by factor 1/2 in closed-shell
    IF (nspins == 1) CALL dbcsr_scale(matrix_p(1)%matrix, 0.5_dp)

    ! Transform P, KS, and Harris kernel matrix into Orthonormal basis
    CALL dbcsr_create(matrix_s_sqrt, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_s_sqrt_inv, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

    SELECT CASE (s_sqrt_method)
       CASE (ls_s_sqrt_proot)
           CALL matrix_sqrt_proot(matrix_s_sqrt, matrix_s_sqrt_inv, &
                                  matrix_s(1)%matrix, eps_filter, &
                                  s_sqrt_order, eps_lanczos, max_iter_lanczos, symmetrize=.TRUE.)
       CASE (ls_s_sqrt_ns)
          CALL matrix_sqrt_Newton_Schulz(matrix_s_sqrt, matrix_s_sqrt_inv, &
                                         matrix_s(1)%matrix, eps_filter, &
                                         s_sqrt_order, eps_lanczos, max_iter_lanczos)
       CASE DEFAULT
          CPABORT("Unknown sqrt method.")
    END SELECT

    ! Transform into orthogonal basis
    DO ispin = 1, nspins
       CALL transform_m_orth(matrix_p(ispin)%matrix, matrix_s_sqrt, eps_filter)
       CALL transform_m_orth(matrix_ks(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)
       CALL transform_m_orth(matrix_nsc(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)
    END DO

    !----------------------------------------
    ! Create non-symmetric matrices: Ax, cg, res
    ! Create symmetric matrix: z
    !----------------------------------------

    CALL dbcsr_create(matrix_tmp, template=matrix_s(1)%matrix, matrix_type=dbcsr_type_no_symmetry)

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
       CALL dbcsr_create(matrix_z(ispin)%matrix, name = "Z-Matrix", &
                         template = matrix_s(1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_z0(ispin)%matrix, name = "p after precondi-Matrix", &
                         template = matrix_s(1)%matrix,&
                         matrix_type=dbcsr_type_no_symmetry)
    END DO

    !----------------------------------------
    ! Get righ-hand-side operators
    !----------------------------------------

    IF(para_env%mepos==0) WRITE(*,*) "ec_response_ao - Build RHS"
    ! Spin factor
    focc = -2.0_dp
    IF (nspins == 1) focc = -4.0_dp

    ! E^[1]_Harris = -4*G[\delta P]*Pin - Pin*G[\delta P] = -4*[G[\delta P], Pin]
    CALL commutator(matrix_nsc, matrix_p, matrix_res, eps_filter, .FALSE., alpha = focc)

    ! Initial guess cg_Z 
    DO ispin = 1, nspins
       CALL dbcsr_copy(matrix_cg_z(ispin)%matrix, matrix_res(ispin)%matrix)
    END DO

    !---------
    ! Matrix Z from MO
    !---------
!    CALL build_hessian_op_from_z(qs_env, p_env, matrix_ks, matrix_p, matrix_s_sqrt, &
!                                 matrix_s_sqrt_inv, matrix_Ax, matrix_res, eps_filter)
    !---------
    norm_zmo = dbcsr_frobenius_norm(ec_env%p_env%p1(1)%matrix)

    ! Projector on trial matrix
    CALL projector(qs_env, matrix_p, matrix_cg_z, eps_filter, .FALSE.)

    ! Ax0
    CALL build_hessian_op(qs_env = qs_env, &
                          p_env = p_env,&
                          matrix_ks = matrix_ks, &
                          matrix_p = matrix_p, &   ! p
                          matrix_s_sqrt_inv = matrix_s_sqrt_inv, &   
                          matrix_cg = matrix_cg_z, & ! cg
                          matrix_Ax = matrix_Ax, &
                          eps_filter = eps_filter)

    ! r_0 = b - Ax0
    DO ispin = 1, nspins
       CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -1.0_dp)
    END DO

    ! Matrix projector T
    CALL projector(qs_env, matrix_p, matrix_res, eps_filter, .TRUE.)

    ! Preconditioner
    linres_control%flag = ""
    DO ispin = 1, nspins
       IF (do_precond) THEN
          ! z_0 = M * r_0 
          ! M ASYM, r_0 ASYM
          IF(para_env%mepos==0) WRITE(*,*) "z_0 = M * r_0"
          ! M*matrix_res = Asym*asym = no symmetry 
          linres_control%flag = "PCG-AO-ORTHO"
       ELSE
          ! z_0 = r_0
          CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
          linres_control%flag = "CG-AO-ORTHO"
       END IF
    END DO

    norm_res = 0.0_dp

    DO ispin = 1, nspins
       ! cg = p_0 = z_0
       CALL dbcsr_copy(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix)

       ! Tr(r_0 * z_0)
       CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_cg(ispin)%matrix, norm_rr(ispin))
       IF(para_env%mepos==0) WRITE(*,*) "Tr(r_0 * z_0) (norm_rr)", norm_rr(1)

       IF (norm_rr(ispin) .LT. 0.0_dp) CPABORT("norm_rr < 0")
       norm_res = MAX(norm_res, ABS(norm_rr(ispin)/REAL(nao,dp)))
       IF(para_env%mepos==0) WRITE(*,*) "norm_res", norm_res
    END DO

    ! eigenvalue shifting
    min_shift = 0.0_dp
    old_conv = norm_rr(1) 
    shift = MIN(10.0_dp, MAX(min_shift, 0.05_dp*old_conv))
    old_conv= 100.0_dp

    ! header
    IF (iounit > 0) THEN
       IF(para_env%mepos==0) WRITE (iounit, "(/,T3,A,T16,A,T25,A,T38,A,T52,A,/,T3,A)") &
          "Iteration", "Method", "Stepsize", "Convergence", "Time", &
          REPEAT("-", 80)
    ENDIF

    alpha(:) = 0.0_dp
    restart = .FALSE.
    should_stop = .FALSE.

    ! start iteration
    iteration: DO i = 1, linres_control%max_iter 
    
       ! Convergence criteria
       linres_control%converged = .FALSE.
       ! default for eps 10E-6 in MO_linres
       IF (norm_res .LT. linres_control%eps) THEN
          linres_control%converged = .TRUE.
       ENDIF

       t2 = m_walltime()
       IF (iter .EQ. 1 .OR. MOD(iter, 1) .EQ. 0 .OR. linres_control%converged &
           .OR. restart .OR. should_stop) THEN
          IF (iounit > 0) THEN
             WRITE (iounit, "(T5,I5,T18,A3,T28,L1,T38,1E8.2,T48,F16.10,T68,F8.2)") &
                i, linres_control%flag, restart, MAXVAL(alpha), norm_res, t2 - t1
             CALL m_flush(iounit)
          ENDIF
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

       ! Max number of iteration reached
       IF (i == linres_control%max_iter) THEN
          IF (iounit > 0) THEN
             IF(para_env%mepos==0) WRITE (iounit, "(/,T2,A/)") &
                "The linear solver didnt converge! Maximum number of iterations reached."
             CALL m_flush(iounit)
          ENDIF
          linres_control%converged = .FALSE.
       ENDIF

       ! Hessian Ax = [F,B] + [G(B),P]
       CALL build_hessian_op(qs_env = qs_env, &
                             p_env=p_env,&
                             matrix_ks = matrix_ks, &
                             matrix_p = matrix_p, &   ! p
                             matrix_s_sqrt_inv = matrix_s_sqrt_inv, &   
                             matrix_cg = matrix_cg, & ! cg
                             matrix_Ax = matrix_Ax,&
                             eps_filter = eps_filter)

       ! Matrix projector T
       CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter, .TRUE.)

       ! Apply shift
       DO ispin = 1, nspins
          CALL dbcsr_add(matrix_Ax(ispin)%matrix, matrix_cg(ispin)%matrix, 1.0_dp, shift)
       END DO

       DO ispin = 1, nspins

          CALL dbcsr_filter(matrix_Ax(ispin)%matrix, eps_filter)
          ! norm_cA = tr(Ap_j * p_j)
          CALL dbcsr_dot(matrix_cg(ispin)%matrix, matrix_Ax(ispin)%matrix, norm_cA(ispin))
          IF(para_env%mepos==0) WRITE(*,*) "norm_cA", norm_cA(ispin)

          IF (norm_cA(ispin) .LT. 0.0_dp) THEN

             ! Recalculate w/o preconditioner
             IF (i > 1) THEN
                ! p = -z + beta*p
                CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix, beta(ispin), -1.0_dp)
                CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_res(ispin)%matrix, new_norm(ispin))
                beta(ispin) = new_norm(ispin)/tr_rz00(ispin)
                CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_res(ispin)%matrix, beta(ispin), 1.0_dp)
                norm_rr(ispin) = new_norm(ispin)
             ELSE
                CALL dbcsr_copy(matrix_res(ispin)%matrix, matrix_cg(ispin)%matrix)
                CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_res(ispin)%matrix, norm_rr(ispin))
             END IF

             CALL build_hessian_op(qs_env = qs_env, &
                                   p_env = p_env,&
                                   matrix_ks = matrix_ks, &
                                   matrix_p = matrix_p, &   ! p
                                   matrix_s_sqrt_inv = matrix_s_sqrt_inv, &   
                                   matrix_cg = matrix_cg, & ! cg
                                   matrix_Ax = matrix_Ax, &
                                   eps_filter = eps_filter)

             ! Matrix projector T
             CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter, .TRUE.)

             CALL dbcsr_dot(matrix_cg(ispin)%matrix, matrix_Ax(ispin)%matrix, norm_cA(ispin))

             CPABORT("tr(Ap_j*p_j) < 0")

          END IF

       END DO

       DO ispin = 1, nspins
          lin_eq = 0.0_dp
          fac = norm_rr(ispin)/norm_cA(ispin)
          fac1 = fac
          ! Use a 3point line serach and a fit to a quadratic function to determine optimal step size
          DO iter = 1, 3
             CALL dbcsr_copy(matrix_tmp, matrix_res(ispin)%matrix)
             CALL dbcsr_add(matrix_tmp, matrix_Ax(ispin)%matrix, 1.0_dp, -fac)
             devi(iter) = dbcsr_frobenius_norm(matrix_tmp)
             lin_eq(iter, :) = (/fac**2, fac, 1.0_dp/)
             fac = fac1 + fac1*((-1)**iter)*0.5_dp
          END DO
          CALL invmat(lin_eq, info)
          vec = MATMUL(lin_eq, devi)
          astep = -vec(2)/(2.0_dp*vec(1))
          IF(para_env%mepos==0) WRITE(*,*) "astep", astep
          fac = SQRT(norm_rr(ispin)/(norm_cA(ispin)*astep))
    !scale the previous matrices to match the step size
          CALL dbcsr_scale(matrix_Ax(ispin)%matrix, fac)
          CALL dbcsr_scale(matrix_cg(ispin)%matrix, fac)
          norm_cA(ispin) = norm_cA(ispin)*fac**2
       END DO

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
!          CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -alpha(ispin))
       END DO

       ! need to recompute the residue
       restart = .FALSE.
       IF (MOD(iter, linres_control%restart_every) .EQ. 0) THEN
          !
          ! r_j+1 = b - A * x_j+1
          CALL build_hessian_op(qs_env = qs_env, &
                                p_env = p_env, &
                                matrix_ks = matrix_ks, &
                                matrix_p = matrix_p, &
                                matrix_s_sqrt_inv = matrix_s_sqrt_inv, &
                                matrix_cg = matrix_cg_z, & ! cg
                                matrix_Ax = matrix_Ax, &
                                eps_filter = eps_filter)
          ! b 
          CALL commutator(matrix_nsc, matrix_p, matrix_res, eps_filter, .FALSE., alpha = focc)
   
          DO ispin = 1, nspins
             CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -1.0_dp)
          END DO
   
          CALL projector(qs_env, matrix_p, matrix_res, eps_filter, .TRUE.)
          !
          restart = .TRUE.
       ELSE
          ! proj Ap onto the virtual subspace
          CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter, .TRUE.)
          !
          ! r_j+1 = r_j - alpha * Ap_j
          DO ispin = 1, nspins
             CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -alpha(ispin))
          ENDDO
          restart = .FALSE.
       ENDIF

       !-------------
       ! testing of Z Matrix

       ! Z = [cg_z,P]
       CALL dbcsr_copy(matrix_z0(1)%matrix, matrix_cg_z(1)%matrix)
       CALL projector(qs_env, matrix_p, matrix_z0, eps_filter, .FALSE.)
       CALL commutator(matrix_z0, matrix_p, matrix_z, eps_filter, .TRUE., alpha = 0.5_dp)

       ! Transform Z-matrix back into non-orthogonal basis
       CALL transform_m_orth(matrix_z(1)%matrix, matrix_s_sqrt_inv, eps_filter)

       norm_zao = dbcsr_frobenius_norm(matrix_z(1)%matrix)
       IF (para_env%mepos==0) WRITE(*,*) "|| Z (MO) ||_F :", norm_zmo
       IF (para_env%mepos==0) WRITE(*,*) "|| Z (AO) ||_F :", norm_zao
       IF (para_env%mepos==0) WRITE(*,*) "It, conv, ||Z(AO)||_F :", i, norm_res, norm_zao
       !------------

       ! Preconditioner
       linres_control%flag = ""
       IF (do_precond) THEN
          ! z_j+1 = M * r_j+1
          linres_control%flag = "PCG-AO-ORTHO"
       ELSE 
          DO ispin = 1, nspins
             CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
          END DO
          linres_control%flag = "CG-AO-ORTHO"
       END IF

       norm_res = 0.0_dp

       DO ispin = 1, nspins
          ! Tr[r_j+1*z_j+1]
          CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_z0(ispin)%matrix, new_norm(ispin))
          IF (new_norm(ispin) .LT. 0.0_dp) CPABORT("tr(r_j+1*z_j+1) < 0")
          IF (para_env%mepos==0) WRITE(*,*) "new_norm", new_norm(1)
          norm_res = MAX(norm_res, new_norm(ispin)/REAL(nao,dp))

          IF (norm_rr(ispin) .LT. linres_control%eps .OR. new_norm(ispin) .LT. linres_control%eps) THEN
             beta(ispin) = 0.0_dp
             linres_control%converged = .TRUE.
          IF (para_env%mepos==0) WRITE(*,*) "beta = 0.0 CONVERGED"
          ELSE
             beta(ispin) = new_norm(ispin)/norm_rr(ispin)
          END IF
          IF (para_env%mepos==0) WRITE(*,*) "beta", beta(1)

          ! update new search vector (matrix cg)
          ! Here: cg_j+1 = z_j+1 + beta*cg_j
          CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix, beta(ispin), 1.0_dp)
          CALL dbcsr_filter(matrix_cg(ispin)%matrix, eps_filter)

          tr_rz00(ispin) = norm_rr(ispin)
          norm_rr(ispin) = new_norm(ispin)
       END DO

       ! Can we exit the SCF loop?
       CALL external_control(should_stop, "LS_SOLVER", target_time=qs_env%target_time, &
                             start_time=qs_env%start_time)

    END DO iteration 

    ! Matrix projector
    CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter, .FALSE.)

    ! Z = [cg_z,P]
    CALL commutator(matrix_cg_z, matrix_p, matrix_z, eps_filter, .TRUE., alpha = 0.5_dp)

    DO ispin = 1, nspins
       ! Transform Z-matrix back into non-orthogonal basis
       CALL transform_m_orth(matrix_z(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)

       ! Export Z-Matrix 
       CALL dbcsr_copy(p_env%p1(ispin)%matrix, matrix_z(ispin)%matrix, keep_sparsity=.TRUE.)
    END DO

    norm_zao = dbcsr_frobenius_norm(matrix_z(1)%matrix)
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
    CALL ec_wz_matrix(qs_env, p_env%p1, p_env%w1, eps_filter)
 
    ! Release matrices
    CALL dbcsr_release(matrix_tmp)

    CALL dbcsr_release(matrix_s_sqrt)
    CALL dbcsr_release(matrix_s_sqrt_inv)
    
    CALL dbcsr_deallocate_matrix_set(matrix_p)
    CALL dbcsr_deallocate_matrix_set(matrix_ks)
    CALL dbcsr_deallocate_matrix_set(matrix_nsc)
    CALL dbcsr_deallocate_matrix_set(matrix_z)
    CALL dbcsr_deallocate_matrix_set(matrix_Ax)
    CALL dbcsr_deallocate_matrix_set(matrix_res)
    CALL dbcsr_deallocate_matrix_set(matrix_cg)
    CALL dbcsr_deallocate_matrix_set(matrix_cg_z)
    CALL dbcsr_deallocate_matrix_set(matrix_z0)

!    CALL linres_control_release(linres_control)

    DEALLOCATE(alpha, beta, new_norm, norm_cA, norm_rr)
    DEALLOCATE(tr_rz00)

    IF(para_env%mepos==0) WRITE(*,*) "ec_response_ao - END"

    CALL timestop(handle)

  END SUBROUTINE ec_response_ao

! **************************************************************************************************
!> \brief  calculate lin transformation of Hessian matrix on a trial vector (matrix) matrix_cg 
!>         which is stored in response density B = [cg,P] = cg*P - P*cg = cg*P + (cg*P)^T
!>         Ax = [F, B] + [G(B), Pin] in orthonormal basis
!>
!>          
!> \param qs_env ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE build_hessian_op(qs_env, p_env, matrix_ks, matrix_p, matrix_s_sqrt_inv, &
                              matrix_cg, matrix_Ax, eps_filter)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(qs_p_env_type), POINTER                       :: p_env
    TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
                                        POINTER        :: matrix_cg, matrix_ks, matrix_p
    TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s_sqrt_inv 
    TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
                                        POINTER        :: matrix_Ax
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'build_hessian_op', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle, ispin, nspins
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_rho_type), POINTER                         :: rho
    TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_b, rho1_ao
    TYPE(dft_control_type), POINTER                    :: dft_control
    REAL(KIND=dp)                                      :: chksum
!    REAL(KIND=dp)                                      :: norm_res

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    para_env=para_env,&
                    rho = rho)
    nspins = dft_control%nspins
    !IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op BEGINNING"

    NULLIFY(matrix_b)
    CALL dbcsr_allocate_matrix_set(matrix_b, nspins)
    DO ispin = 1, nspins
       ALLOCATE (matrix_b(ispin)%matrix)
       CALL dbcsr_create(matrix_b(ispin)%matrix, name = "[X,P] RSP DNSTY", &
                         template = matrix_p(1)%matrix,&
                         matrix_type = dbcsr_type_no_symmetry)
    END DO

    ! Build intermediate density matrix
    ! cg antisym, P sym => B sym
    ! B = [cg, P] = cg*P - P*cg = cg*P + P^T*cg*T = cg*P + (cg*P)^T
    CALL commutator(matrix_cg, matrix_p, matrix_b, eps_filter, .TRUE.)

    
    ! Build first part of operator
    ! Ax = [F,[X,P]] = [F,B]
    ! Ax = F*B - B*F
    ! F sym, B sym use comm sym => Ax antisym
    CALL commutator(matrix_ks, matrix_b, matrix_Ax, eps_filter, .FALSE.)
    
!    norm_res = dbcsr_frobenius_norm(matrix_Ax(1)%matrix)
!    IF (para_env%mepos==0) WRITE(*,*) "||Ax1||", norm_res
!    norm_res = dbcsr_frobenius_norm(matrix_b(1)%matrix)
!    IF (para_env%mepos==0) WRITE(*,*) "||B||", norm_res

    ! Check 2nd part
    chksum = 0.0_dp
    DO ispin = 1, nspins
       chksum = chksum+dbcsr_checksum(matrix_b(ispin)%matrix)
    ENDDO
    IF (para_env%mepos==0) WRITE(*,*) "chksum:", chksum

    ! skip the kernel if the DM is very small
    IF (chksum .GT. 1.0E-14_dp) THEN

       CALL p_env_check_i_alloc(p_env, qs_env)

       CALL qs_rho_get(p_env%rho1, rho_ao=rho1_ao)

       ! Transform B into NON-ortho basis for collocation
       DO ispin = 1, nspins
          CALL transform_m_orth(matrix_b(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)
          ! Keep sparsity of density matrix 
          CALL dbcsr_copy(rho1_ao(ispin)%matrix, matrix_b(ispin)%matrix, keep_sparsity = .TRUE.)
       END DO

       ! Update grids to trial response density 
       ! Same basis as ground-state calculation
       CALL qs_rho_update_rho(p_env%rho1, qs_env)

       ! Calculate kernel
       ! Ax = F*B - B*F + G(B)*P - P*G(B)
       !                        IN/OUT     IN        IN             IN
       CALL hessian_op2(qs_env, matrix_Ax, matrix_p, matrix_s_sqrt_inv, p_env%rho1, eps_filter)

    END IF

    CALL dbcsr_deallocate_matrix_set(matrix_b)

    CALL timestop(handle)

  END SUBROUTINE 

! **************************************************************************************************
!> \brief  
!>         
!> matrix_Ax 
!> matrix_p Density matrix in orthogonal basis         
!> param rho_b Response density on grid          
!>          
!> \param qs_env ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE hessian_op2(qs_env, matrix_Ax, matrix_p, matrix_s_sqrt_inv, rho1, eps_filter)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
                                        POINTER        :: matrix_Ax 
    TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
                                        POINTER        :: matrix_p 
    TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s_sqrt_inv 
    TYPE(qs_rho_type), INTENT(IN), POINTER             :: rho1
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'hessian_op2', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle, ispin, nspins
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_rho_type), POINTER                         :: rho
    TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, matrix_G, rho_ao
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

    CALL get_qs_env(qs_env,para_env=para_env)
    !IF(para_env%mepos==0) WRITE(*,*) "hessian_op2 BEGINNING"

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    input = input, &
                    matrix_s = matrix_s,&
                    rho = rho)
    nspins = dft_control%nspins

    ! Get non-ortho input density matrix on grid
    CALL qs_rho_get(rho, rho_ao = rho_ao)

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
    CALL create_kernel(qs_env, &
                       vxc = v_xc, &
                       rho = rho, &
                       rho1 = rho1, &
                       xc_section = xc_section)

    DO ispin = 1,nspins
       CALL pw_scale(v_xc(ispin)%pw,v_xc(ispin)%pw%pw_grid%dvol)
    END DO

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

    ! Init response kernel matrix
    NULLIFY (matrix_G)
    CALL dbcsr_allocate_matrix_set(matrix_G, nspins)
    DO ispin = 1, nspins
       ALLOCATE (matrix_G(ispin)%matrix)
       CALL dbcsr_copy(matrix_G(ispin)%matrix, matrix_s(1)%matrix, &
                       name="MATRIX Kernel")

       ! use input density matrix, but in ortho basis
       ! matrix G(B) - the big price
       CALL integrate_v_rspace(v_rspace = v_xc(ispin), &
                               pmat = rho_ao(ispin), &
                               hmat = matrix_G(ispin), &
                               qs_env = qs_env, &
                               calculate_forces = .FALSE., &
                               basis_type = "ORB")
    END DO

    ! Transforming G(B) into ortho
    DO ispin = 1, nspins
       CALL transform_m_orth(matrix_G(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)
    END DO

    ! Hessian already contains  Ax = [F,B] (ORTHO), now adding
    ! Ax = Ax + G(B)P - (G(B)P)^T
    ! G(B) and P sym => Ax2 antisym
    CALL commutator(matrix_G, matrix_p, matrix_Ax, eps_filter, .FALSE., 1.0_dp, 1.0_dp)

    ! release pw grids
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
    DO ispin = 1,nspins
       CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
    END DO
    DEALLOCATE(v_xc)

    !IF(para_env%mepos==0) WRITE(*,*) "hessian_op2 END"

    CALL dbcsr_deallocate_matrix_set(matrix_G)
    
    CALL timestop(handle)

  END SUBROUTINE hessian_op2

! **************************************************************************************************
!> \brief  calculate lin transformation of Hessian matrix on a trial vector (matrix) matrix_cg 
!>         which is stored in response density B = [cg,P] = cg*P - P*cg = cg*P + (cg*P)^T
!>         Ax = [F, B] + [G(B), Pin] in orthonormal basis
!>
!>          
!> \param qs_env ...
!> \date    12.2019
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE build_hessian_op_from_z(qs_env, p_env, matrix_ks, matrix_p, matrix_s_sqrt, &
                                     matrix_s_sqrt_inv, matrix_Ax, matrix_res, eps_filter)

     TYPE(qs_environment_type), POINTER                 :: qs_env
     TYPE(qs_p_env_type), POINTER                       :: p_env
     TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
                                         POINTER        :: matrix_ks, matrix_p
     TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s_sqrt, matrix_s_sqrt_inv
     TYPE(dbcsr_p_type), DIMENSION(:), INTENT(OUT), &
                                         POINTER        :: matrix_Ax
     TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
                                         POINTER        :: matrix_res
     REAL(KIND=dp), INTENT(IN)                          :: eps_filter

     CHARACTER(len=*), PARAMETER :: routineN = 'build_hessian_op', routineP = moduleN//':'//routineN

     INTEGER                                            :: handle, ispin, nspins
     TYPE(cp_para_env_type), POINTER                    :: para_env
     TYPE(qs_rho_type), POINTER                         :: rho, rho_b
     TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_b, rho_ao_b
     TYPE(dft_control_type), POINTER                    :: dft_control
     REAL(KIND=dp)                                      :: chksum
     REAL(KIND=dp)                                      :: norm_ax, norm_b

     CALL timeset(routineN, handle)

     CALL get_qs_env(qs_env,para_env=para_env)
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z BEGINNING"

     CALL get_qs_env(qs_env = qs_env, &
                     dft_control = dft_control, &
                     rho = rho)
     nspins = dft_control%nspins

     NULLIFY(matrix_b)
     CALL dbcsr_allocate_matrix_set(matrix_b, nspins)
     DO ispin = 1, nspins
        ALLOCATE (matrix_b(ispin)%matrix)
        CALL dbcsr_create(matrix_b(ispin)%matrix, name = "[X,P] RSP DNSTY", &
                          template = matrix_p(1)%matrix,&
                          matrix_type = dbcsr_type_no_symmetry)
     END DO

     ! Z nonortho -> Z ortho
     DO ispin = 1, nspins
        CALL dbcsr_desymmetrize(p_env%p1(ispin)%matrix, matrix_b(ispin)%matrix)
        CALL transform_m_orth(matrix_b(ispin)%matrix, matrix_s_sqrt, eps_filter)
        ! Z matrix stored in ec_env was scaled with *0.5 in closed shell case
        CALL dbcsr_scale(matrix_b(ispin)%matrix, 2.0_dp)
     END DO

     ! Build first part of operator
     ! Ax = [F,[X,P]] = [F,B]
     ! Ax = F*B - B*F
     ! F sym, B sym use comm sym => Ax antisym
     CALL commutator(matrix_ks, matrix_b, matrix_Ax, eps_filter, .FALSE.)

     ! Check 2nd part
     chksum = 0.0_dp
     DO ispin = 1, nspins
        chksum = chksum+dbcsr_checksum(matrix_b(ispin)%matrix)
     ENDDO
     IF (para_env%mepos==0) WRITE(*,*) "chksum:", chksum

     ! skip the kernel if the DM is very small
     IF (chksum .GT. 1.0E-14_dp) THEN

        ! Bring B as density on grid
        NULLIFY(rho_b)
        CALL qs_rho_create(rho_b)
        CALL duplicate_rho_type(rho_input = rho, &
                                rho_output = rho_b, &
                                qs_env = qs_env)

        ! Get response density matrix
        CALL qs_rho_get(rho_b, rho_ao = rho_ao_b)

        ! Transform B into NON-ortho basis for collocation
        DO ispin = 1, nspins
           CALL transform_m_orth(matrix_b(ispin)%matrix, matrix_s_sqrt_inv, eps_filter)
           ! Keep sparsity of density matrix 
           CALL dbcsr_copy(rho_ao_b(ispin)%matrix, matrix_b(ispin)%matrix, keep_sparsity = .TRUE.)
        END DO

        ! update rho_r_b and rho_g_b wrt rho_ao_b
        ! Same basis as reference calculation
        CALL qs_rho_update_rho(rho_b, qs_env)

        ! Calculate kernel
        ! Ax = F*B - B*F + G(B)*P - P*G(B)
        !                        IN/OUT     IN        IN             IN
        CALL hessian_op2(qs_env, matrix_Ax, matrix_p, matrix_s_sqrt_inv, rho_b, eps_filter)

        CALL qs_rho_release(rho_b)

     END IF

     CALL projector(qs_env, matrix_p, matrix_res, eps_filter, .TRUE.)
     CALL projector(qs_env, matrix_p, matrix_Ax, eps_filter, .TRUE.)

     norm_b  = dbcsr_frobenius_norm(matrix_res(1)%matrix)
     norm_ax = dbcsr_frobenius_norm(matrix_Ax(1)%matrix)
     IF (para_env%mepos==0) WRITE(*,*) "|| b_ao  ||_F :", norm_b
     IF (para_env%mepos==0) WRITE(*,*) "|| Az_ao ||_F :", norm_ax

     CALL dbcsr_set(matrix_Ax(1)%matrix, 0.0_dp)

     CALL dbcsr_deallocate_matrix_set(matrix_b)
     IF(para_env%mepos==0) WRITE(*,*) "build_hessian_op_from_z ---------"

     CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief computes a commutator exploiting symmetry RES= beta*RES + k*[A,B] = k*[AB-(AB)T]
!> \param a ...
!> \param b ...
!> \param res        Commutator result 
!> \param eps_filter filtering threshold for sparse matrices
!> \param alpha      prefactor alpha in above equation
!> \param beta       prefactor beta in above equation
!> \param anticomm   Calculate anticommutator
!>
!> \par History
!>       2020.07 adapted Fabian Belleflamme  (based on commutator_symm)
! **************************************************************************************************

   SUBROUTINE commutator(a, b, res, eps_filter, anticomm, alpha, beta)

      TYPE(dbcsr_p_type), DIMENSION(:)                   :: a, b, res
      REAL(KIND=dp)                                      :: eps_filter
      LOGICAL                                            :: anticomm
      REAL(KIND=dp), OPTIONAL                            :: alpha, beta

      CHARACTER(LEN=*), PARAMETER :: routineN = 'commutator_symm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: facc, myalpha, mybeta
      TYPE(dbcsr_type)                                   :: work, work2

      CALL timeset(routineN, handle)

      CALL dbcsr_create(work, template=a(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(work2, template=a(1)%matrix, matrix_type=dbcsr_type_no_symmetry)

      myalpha = 1.0_dp
      IF (PRESENT(alpha)) myalpha = alpha
      mybeta = 0.0_dp
      IF (PRESENT(beta)) mybeta = beta
      facc = -1.0_dp
      IF (anticomm) facc = 1.0_dp

      DO ispin = 1, SIZE(a)

         CALL dbcsr_multiply("N", "N", myalpha, a(ispin)%matrix, b(ispin)%matrix, &
                             0.0_dp, work, filter_eps=eps_filter)
         CALL dbcsr_transposed(work2, work)

         ! RES= beta*RES + alpha*{A,B} = beta*RES + alpha*[AB+(AB)T]
         ! RES= beta*RES + alpha*[A,B] = beta*RES + alpha*[AB-(AB)T]
         CALL dbcsr_add(work, work2, 1.0_dp, facc)

         CALL dbcsr_add(res(ispin)%matrix, work, mybeta, 1.0_dp)

      END DO

      CALL dbcsr_release(work)
      CALL dbcsr_release(work2)

      CALL timestop(handle)

   END SUBROUTINE commutator

! **************************************************************************************************
!> \brief transform from ORTH (NON-ORTH) to NON-ORTH (ORTH) 
!>        multiply left&right with (1) S^1/2, or (2) S^1/2_inv
!>        symmetrize: (T) if matri_in symmetric, or (f) if antisymmetric
!>        
!>        
!>          
!> \param qs_env ...
!> \date    01.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE transform_ortho(qs_env, matrix_in, which_trafo, symmetrize)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_type), POINTER, INTENT(INOUT)           :: matrix_in
    INTEGER, INTENT(IN)                                :: which_trafo
    LOGICAL, INTENT(IN)                                :: symmetrize

    CHARACTER(len=*), PARAMETER :: routineN = 'transform_ortho', routineP = moduleN//':'//routineN

    INTEGER                                            :: handle
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    REAL(KIND=dp)                                      :: eps_filter, eps_lanczos
    TYPE(dbcsr_type)                                   :: matrix_s_sqrt, matrix_s_sqrt_inv
    INTEGER                                            :: ls_sqrt_method, max_iter_lanczos, s_sqrt_order
    TYPE(dbcsr_type)                                   :: matrix_tmp, matrix_work
    TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    matrix_s = matrix_s, &
                    para_env=para_env)
    eps_filter = 1.0E-7_dp

    IF(para_env%mepos==0) WRITE(*,*) "transform_ortho BEGINNING"

    !---------------------------------
    CALL dbcsr_create(matrix_s_sqrt, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_s_sqrt_inv, template=matrix_s(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)

    ls_sqrt_method = 2
    s_sqrt_order = 3
    eps_lanczos = 1.0E-3_dp
    max_iter_lanczos = 128
    SELECT CASE (ls_sqrt_method)
    CASE (ls_s_sqrt_proot)
       CALL matrix_sqrt_proot(matrix_s_sqrt, matrix_s_sqrt_inv, &
                              matrix_s(1)%matrix, eps_filter, &
                              s_sqrt_order, &
                              eps_lanczos, max_iter_lanczos, &
                              symmetrize = .TRUE.)
    CASE(ls_s_sqrt_ns)
       CALL matrix_sqrt_Newton_Schulz(matrix_s_sqrt, matrix_s_sqrt_inv, &
                                      matrix_s(1)%matrix, eps_filter, &
                                      s_sqrt_order, &
                                      eps_lanczos, max_iter_lanczos)
    END SELECT

    CALL dbcsr_create(matrix_work, template=matrix_in, matrix_type=dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_tmp, template=matrix_in, matrix_type=dbcsr_type_no_symmetry)

    SELECT CASE (which_trafo)
    CASE(1) ! S^1/2 * M * S^1/2

       ! matrix_trafo = matrix_s_sqrt
       CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_in, matrix_s_sqrt, &
                           0.0_dp, matrix_work, filter_eps=eps_filter)
       CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_s_sqrt, matrix_work, &
                           0.0_dp, matrix_tmp, filter_eps=eps_filter)

    CASE(2) ! S!1/2_inv * M *S^1/2_inv

       ! matrix_trafo = matrix_s_sqrt_inv
       CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_in, matrix_s_sqrt_inv, &
                           0.0_dp, matrix_work, filter_eps=eps_filter)
       CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_s_sqrt_inv, matrix_work, &
                           0.0_dp, matrix_tmp, filter_eps=eps_filter)
    END SELECT


    IF (symmetrize) THEN
       WRITE(*,*) "SYMMETRIZED"
       CALL dbcsr_transposed(matrix_work, matrix_tmp)
       CALL dbcsr_add(matrix_tmp, matrix_work, 0.5_dp, 0.5_dp)
       CALL dbcsr_copy(matrix_in, matrix_tmp)
    ELSE 
       WRITE(*,*) "ANTI-SYMMETRIZED"
       CALL dbcsr_transposed(matrix_work, matrix_tmp)
       CALL dbcsr_add(matrix_tmp, matrix_work, 0.5_dp, -0.5_dp)
       CALL dbcsr_copy(matrix_in, matrix_tmp)
    END IF

    CALL dbcsr_release(matrix_tmp)
    CALL dbcsr_release(matrix_work)
    CALL dbcsr_release(matrix_s_sqrt)
    CALL dbcsr_release(matrix_s_sqrt_inv)

    CALL timestop(handle)

  END SUBROUTINE

! **************************************************************************************************
!> \brief Projector P(M) = P*M*Q^T + Q*M*P^T
!>        with P = D
!>        with Q = (1-D)
!> 
!> \param qs_env ...
!> \param matrix_ps Precomputed prodcut D*S
!> \param matrix_io Matrix to which projector is applied.
!>          
!> \date    06.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
!  SUBROUTINE projector(qs_env, matrix_p, matrix_io, eps_filter, ptrans)
!
!    TYPE(qs_environment_type), POINTER                 :: qs_env
!    TYPE(dbcsr_p_type), DIMENSION(:), &
!                            INTENT(IN), POINTER        :: matrix_p
!    TYPE(dbcsr_p_type), DIMENSION(:), &
!                         INTENT(INOUT), POINTER        :: matrix_io
!    REAL(KIND=dp), INTENT(IN)                          :: eps_filter
!    LOGICAL, INTENT(IN), OPTIONAL                      :: ptrans
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'projector', routineP = moduleN//':'//routineN
!
!    INTEGER                                            :: handle, ispin, nspins
!    LOGICAL                                            :: my_ptrans
!    TYPE(cp_para_env_type), POINTER                    :: para_env
!    TYPE(dft_control_type), POINTER                    :: dft_control
!    TYPE(dbcsr_type)                                   :: matrix_tmp1, matrix_tmp2, matrix_q
!
!    CALL timeset(routineN, handle)
!
!    CALL get_qs_env(qs_env = qs_env, &
!                    dft_control = dft_control, &
!                    para_env=para_env)
!    nspins = dft_control%nspins
!!    IF(para_env%mepos==0) WRITE(*,*) "projector BEGINNING"
!
!    ! Decides which projector to apply 
!    ! FALSE : Projector onto space of non-redudnant elements
!    ! TRUE : It's transpose 
!    my_ptrans = .FALSE.
!    IF (PRESENT(ptrans)) my_ptrans = ptrans
!
!    CALL dbcsr_create(matrix_q, template = matrix_p(1)%matrix, &
!                      matrix_type = dbcsr_type_no_symmetry)
!    CALL dbcsr_create(matrix_tmp1, template = matrix_p(1)%matrix, &
!                      matrix_type = dbcsr_type_no_symmetry)
!    CALL dbcsr_create(matrix_tmp2, template = matrix_p(1)%matrix, &
!                      matrix_type = dbcsr_type_no_symmetry)
!
!    ! Q = (1 - P)
!    CALL dbcsr_copy(matrix_q, matrix_p(1)%matrix)
!    CALL dbcsr_scale(matrix_q, -1.0_dp)
!    CALL dbcsr_add_on_diag(matrix_q, 1.0_dp)
!    CALL dbcsr_finalize(matrix_q)
!
!    IF (my_ptrans) THEN
!       ! Proj^T(M) = P^T*M*Q + Q^T*M*P
!       ! with P = D = CC^T
!       ! and  Q = (1 - P)
!       DO ispin = 1, nspins
!          ! tmp1 = P^T*M
!          CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_p(ispin)%matrix, matrix_io(ispin)%matrix, &
!                              0.0_dp, matrix_tmp1, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!          ! tmp2 = P^T*M*Q
!          CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, matrix_q, &
!                              0.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!
!          ! tmp1 = (P^T*M^T*Q)^T ?=? -(P^T*M*Q)^T
!           CALL dbcsr_transposed(matrix_tmp1, matrix_tmp2)
!           CALL dbcsr_add(matrix_tmp2, matrix_tmp1, 1.0_dp, -1.0_dp)
!
!!          ! tmp1  = Q^T*M
!!          CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_q, matrix_io(ispin)%matrix, &
!!                              0.0_dp, matrix_tmp1, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!!          ! tmp2 = tmp2 + Q^T*M*P
!!          CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp1, matrix_p(ispin)%matrix, &
!!                              1.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!       END DO
!    ELSE
!       ! Proj(M) = P*M*Q^T + Q*M*P^T
!       ! with P = D = CC^T
!       ! and  Q = (1 - P)
!       DO ispin = 1, nspins
!          ! tmp1 = P*M
!          CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin)%matrix, matrix_io(1)%matrix, &
!                              0.0_dp, matrix_tmp1, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!          ! tmp2 = P*M*Q^T
!          CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp1, matrix_q, &
!                              0.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!
!          ! tmp1 = (P*M^T*Q^T)^T ?=? -(P*M*Q^T)^T
!          CALL dbcsr_transposed(matrix_tmp1, matrix_tmp2)
!          CALL dbcsr_add(matrix_tmp2, matrix_tmp1, 1.0_dp, -1.0_dp)
!
!!          ! tmp1  = Q*M
!!          CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_q, matrix_io(ispin)%matrix, &
!!                              0.0_dp, matrix_tmp1, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!!          ! tmp2 = tmp2 + Q*M*P^T
!!          CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp1, matrix_p(ispin)%matrix, &
!!                              1.0_dp, matrix_tmp2, filter_eps = eps_filter, retain_sparsity = .FALSE.) 
!       END DO
!    END IF
!
!    ! Export
!    DO ispin = 1, nspins
!       CALL dbcsr_copy(matrix_io(ispin)%matrix, matrix_tmp2, keep_sparsity = .FALSE.)
!       !CALL dbcsr_complete_redistribute(matrix_tmp2, matrix_io(ispin)%matrix)
!    END DO
!
!    CALL dbcsr_release(matrix_tmp1)
!    CALL dbcsr_release(matrix_tmp2)
!    CALL dbcsr_release(matrix_q)
!
!    CALL timestop(handle)
!
!  END SUBROUTINE

! **************************************************************************************************
!> \brief Projector P(M) = P*M*Q^T + Q*M*P^T
!>        with P = D
!>        with Q = (1-D)
!> 
!> \param qs_env ...
!> \param matrix_ps Precomputed prodcut D*S
!> \param matrix_io Matrix to which projector is applied.
!>          
!> \date    06.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE projector(qs_env, matrix_p, matrix_io, eps_filter, ptrans)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_p_type), DIMENSION(:), &
                            INTENT(IN), POINTER        :: matrix_p
    TYPE(dbcsr_p_type), DIMENSION(:), &
                         INTENT(INOUT), POINTER        :: matrix_io
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter
    LOGICAL, INTENT(IN), OPTIONAL                      :: ptrans

    CHARACTER(len=*), PARAMETER :: routineN = 'projector', routineP = moduleN//':'//routineN

    INTEGER                                            :: handle, ispin, nspins
    LOGICAL                                            :: my_ptrans
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    TYPE(dbcsr_type)                                   :: matrix_tmp, matrix_q

    CALL timeset(routineN, handle)

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    para_env=para_env)
    nspins = dft_control%nspins
!    IF(para_env%mepos==0) WRITE(*,*) "projector BEGINNING"

    ! Decides which projector to apply 
    ! FALSE : Projector onto space of non-redudnant elements
    ! TRUE : It's transpose 
    my_ptrans = .FALSE.
    IF (PRESENT(ptrans)) my_ptrans = ptrans

    CALL dbcsr_create(matrix_q, template = matrix_p(1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_tmp, template = matrix_p(1)%matrix, &
                      matrix_type = dbcsr_type_no_symmetry)

    ! Q = (1 - P)
    CALL dbcsr_copy(matrix_q, matrix_p(1)%matrix)
    CALL dbcsr_scale(matrix_q, -1.0_dp)
    CALL dbcsr_add_on_diag(matrix_q, 1.0_dp)
    CALL dbcsr_finalize(matrix_q)


    DO ispin = 1, nspins

       IF (my_ptrans) THEN

          ! Proj^T(M) = P^T*M*Q + Q^T*M*P
          ! with P = D = CC^T
          ! and  Q = (1 - P)

          ! tmp1 = P^T*M
          CALL dbcsr_multiply("T", "N", 1.0_dp, matrix_p(ispin)%matrix, matrix_io(ispin)%matrix, &
                              0.0_dp, matrix_tmp, filter_eps = eps_filter)
          ! m_io = P^T*M*Q
          CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, matrix_q, &
                              0.0_dp, matrix_io(ispin)%matrix, filter_eps = eps_filter)

          ! tmp1 = (P^T*M^T*Q)^T = -(P^T*M*Q)^T

       ELSE

          ! Proj(M) = P*M*Q^T + Q*M*P^T
          ! with P = D = CC^T
          ! and  Q = (1 - P)

          ! tmp1 = P*M
          CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin)%matrix, matrix_io(1)%matrix, &
                              0.0_dp, matrix_tmp, filter_eps = eps_filter)
          ! tmp2 = P*M*Q^T
          CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_tmp, matrix_q, &
                              0.0_dp, matrix_io(ispin)%matrix, filter_eps = eps_filter)

          ! tmp1 = (P*M^T*Q^T)^T = -(P*M*Q^T)^T

       END IF

       CALL dbcsr_transposed(matrix_tmp, matrix_io(ispin)%matrix)
       CALL dbcsr_add(matrix_io(ispin)%matrix, matrix_tmp, 1.0_dp, -1.0_dp)

    END DO

    CALL dbcsr_release(matrix_tmp)
    CALL dbcsr_release(matrix_q)

    CALL timestop(handle)

  END SUBROUTINE

! **************************************************************************************************
!> \brief performs a tranformation of a matrix back to/into orthonormal basis
!>        in case of P a scaling of 0.5 has to be applied for closed shell case
!> \param matrix       matrix to be transformed
!> \param matrix_trafo transformation matrix
!> \param eps_filter   filtering threshold for sparse matrices
!> \par History
!>       2012.05 created [Florian Schiffmann]
!> \author Florian Schiffmann
!>         
! **************************************************************************************************

   SUBROUTINE transform_m_orth(matrix, matrix_trafo, eps_filter)
      TYPE(dbcsr_type)                                   :: matrix
      TYPE(dbcsr_type)                                   :: matrix_trafo
      REAL(KIND=dp)                                      :: eps_filter

      CHARACTER(LEN=*), PARAMETER :: routineN = 'transform_matrix_orth', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: matrix_tmp, matrix_work

      CALL timeset(routineN, handle)

      CALL dbcsr_create(matrix_work, template=matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp, template=matrix, matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix, matrix_trafo, &
                          0.0_dp, matrix_work, filter_eps=eps_filter)
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_trafo, matrix_work, &
                          0.0_dp, matrix_tmp, filter_eps=eps_filter)
      ! symmetrize results (this is again needed to make sure everything is stable)
      CALL dbcsr_transposed(matrix_work, matrix_tmp)
      CALL dbcsr_add(matrix_tmp, matrix_work, 0.5_dp, 0.5_dp)
      CALL dbcsr_copy(matrix, matrix_tmp)

      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix_work)
      CALL timestop(handle)

   END SUBROUTINE

END MODULE ec_3ao_solver
