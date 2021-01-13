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
MODULE ec_newton
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
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
                                              dbcsr_norm,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_print,&
                                              dbcsr_scale,&
                                              dbcsr_release,&
                                              dbcsr_trace,&
                                              dbcsr_transposed,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric,&
                                              dbcsr_type_no_symmetry,&
                                              dbcsr_type_antisymmetric
   USE dm_ls_scf_qs,                    ONLY: matrix_ls_create,&
                                              matrix_qs_to_ls,&
                                              matrix_ls_to_qs
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
!   USE dm_ls_scf_curvy,                 ONLY: update_p_exp
   USE ec_env_types,                    ONLY: energy_correction_type
   USE input_constants,                 ONLY: kg_tnadd_embed,&
                                              ls_s_sqrt_proot,&
                                              ls_s_sqrt_ns,&
                                              ot_precond_full_single_inverse
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE iterate_matrix,                  ONLY: invert_Hotelling
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_types,                 ONLY: linres_control_create,&
                                              linres_control_release,&
                                              linres_control_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
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
   USE mathconstants,                   ONLY: ifac
   ! fbelle
   USE qs_initial_guess,                ONLY: calculate_atomic_block_dm
      
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_newton'

! *** Public subroutines ***

   PUBLIC :: prepare_newton

CONTAINS

! **************************************************************************************************
!> \brief      Curvy Steps like CG  to solve A(D)X = G(D)
!>             
!>             
!>             
!>          
!> \param qs_env ...
!> \param matrix_ks Harris KS matrix, in 
!> \param qs_env ...
!> \param qs_env ...
!> \param qs_env ...
!>             
!>          
!> \date    01.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
!  SUBROUTINE prepare_newton(qs_env, ksmat, matrix_s, pmat, &
!                            matrix_s_sqrt, matrix_s_sqrt_inv, eps_filter, iounit)
  SUBROUTINE prepare_newton(qs_env, ls_env, iounit)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(ls_scf_env_type), POINTER                     :: ls_env
    INTEGER                                            :: iounit

    CHARACTER(len=*), PARAMETER :: routineN = 'prepare_newton', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: i, handle, ispin, nspins
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    LOGICAL                                            :: converged, restart
    REAL(KIND=dp)                                      :: norm_res
    TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_Ax, matrix_b, matrix_cg, & 
                                                          matrix_cg_z, matrix_res, &
                                                          matrix_tmp, matrix_tmp1, matrix_PKs
!    TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_p 
    REAL(KIND=dp)                                      :: t1, t2
    REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha, beta, new_norm, norm_cA, norm_rr
    INTEGER                                            :: max_iter, nao
    TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
    INTEGER, DIMENSION(2)                              :: nelectron_spin
    TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
    TYPE(section_vals_type), POINTER                   :: ec_section
    CHARACTER(LEN=8)                                   :: flag
    REAL(KIND=dp)                                      :: step_size(2) 
    ! 
    REAL(KIND=dp)                                      :: norm_pao
    TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p_in
    TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_z0
    TYPE(energy_correction_type, POINTER               :: ec_env
    LOGICAL                                            :: do_precond 
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_BCH

    CALL timeset(routineN, handle)

    NULLIFY(dft_control, para_env)
    NULLIFY(atomic_kind_set, qs_kind_set, particle_set)

    t1 = m_walltime()

    CALL get_qs_env(qs_env = qs_env, &
                    dft_control = dft_control, &
                    ec_env = ec_env,&
                    atomic_kind_set = atomic_kind_set, &
                    qs_kind_set = qs_kind_set, &
                    particle_set = particle_set, &
                    nelectron_spin = nelectron_spin, &
                    para_env = para_env)
    nspins = dft_control%nspins
    CALL dbcsr_get_info(ls_env%matrix_ks(1), nfullrows_total=nao)

    IF (iounit > 0) WRITE (iounit,*) repeat("-",80)
    IF (iounit > 0) WRITE (iounit,*) "---- Newton solver ----"
    IF (iounit > 0) WRITE (iounit,*) repeat("-",80)

    ec_section => section_vals_get_subs_vals(qs_env%input, "DFT%ENERGY_CORRECTION")
    !CALL section_vals_val_get(ec_section, "MAX_ITER", i_val=max_iter)
    max_iter = 200
    flag = "HAR"

    IF (iounit > 0) WRITE (iounit,*) "eps", ec_env%eps_default 
    IF (iounit > 0) WRITE (iounit,*) "max_iter", max_iter 

!----------------------------------------------------------
! Init matrices to orthogonalize
!----------------------------------------------------------

    NULLIFY(matrix_p_in)
    CALL dbcsr_allocate_matrix_set(matrix_p_in, nspins)
    DO ispin = 1, nspins
       ALLOCATE (matrix_p_in(ispin)%matrix)
       CALL dbcsr_create(matrix_p_in(ispin)%matrix, name = "P_IN ORTHO", &
                         template = ls_env%matrix_p(1))
    END DO


!----------------------------------------------------------
! EITHER GS OR BLOCK_DM

!    DO ispin = 1, nspins
!       ! get input density matrix ... from ground state calculation
!       !CALL dbcsr_copy(matrix_p(ispin)%matrix, pmat(ispin)) 
!       CALL dbcsr_copy(matrix_p(ispin)%matrix, ls_env%matrix_p(ispin)) 
!    END DO

    ! First guess on Harris density matrix 
    ! in non-orthogonal basis
!    CALL calculate_atomic_block_dm(matrix_p_in, matrix_s, &
!                                   particle_set, atomic_kind_set, qs_kind_set, &
!                                   nspins, nelectron_spin, iounit, para_env)

!----------------------------------------------------------


    norm_res = dbcsr_frobenius_norm(ls_env%matrix_p(1))
    IF (iounit > 0) WRITE (iounit,*)  "||P_out_0||", norm_res

    ! Scale matrix_p by factor 1/2 in closed-shell
    IF (nspins == 1) CALL dbcsr_scale(ls_env%matrix_p(1), 0.5_dp)

    ! Transform initial P_out_0 into orthogonal basis
    DO ispin = 1, nspins
       CALL transform_m_orth(ls_env%matrix_p(ispin), ls_env%matrix_s_sqrt, ls_env%eps_filter)

       CALL transform_m_orth(ls_env%matrix_ks(ispin), ls_env%matrix_s_sqrt_inv, ls_env%eps_filter)
    END DO

    ! Save transformed and scaled input density matrix
    DO ispin = 1, nspins
       CALL dbcsr_copy(matrix_p_in(ispin)%matrix, ls_env%matrix_p(ispin))
    END DO

    norm_res = dbcsr_frobenius_norm(ls_env%matrix_p(1))
    IF (iounit > 0) WRITE (iounit,*)  "||S^1/2*P_out_0*S^1/2||", norm_res

    ALLOCATE(alpha(nspins),beta(nspins), new_norm(nspins), norm_cA(nspins), norm_rr(nspins))
    !----------------------------------------
    ! Create non-symmetric matrices: Ax, cg, res
    !----------------------------------------
    NULLIFY(matrix_Ax, matrix_b, matrix_cg, matrix_cg_z, matrix_res)
    NULLIFY(matrix_tmp, matrix_tmp1)
    NULLIFY(matrix_PKs)
    NULLIFY(matrix_z0)
    CALL dbcsr_allocate_matrix_set(matrix_Ax, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_b, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_cg, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_cg_z, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_res, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_tmp, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_tmp1, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_PKs, nspins)
    CALL dbcsr_allocate_matrix_set(matrix_z0, nspins)
    DO ispin = 1, nspins
       ALLOCATE (matrix_Ax(ispin)%matrix)
       ALLOCATE (matrix_b(ispin)%matrix)
       ALLOCATE (matrix_cg(ispin)%matrix)
       ALLOCATE (matrix_cg_z(ispin)%matrix)
       ALLOCATE (matrix_res(ispin)%matrix)
       ALLOCATE (matrix_tmp(ispin)%matrix)
       ALLOCATE (matrix_tmp1(ispin)%matrix)
       ALLOCATE (matrix_PKs(ispin)%matrix)
       ALLOCATE (matrix_z0(ispin)%matrix)
       CALL dbcsr_create(matrix_Ax(ispin)%matrix, name = "linop MATRIX", &
                         template = ls_env%matrix_ks(1),&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_b(ispin)%matrix, name = "Response MATRIX", &
                         template = ls_env%matrix_ks(1),&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_cg(ispin)%matrix, name = "TRIAL MATRIX", &
                         template = ls_env%matrix_ks(1),&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_cg_z(ispin)%matrix, name = "TRIAL MATRIX", &
                         template = ls_env%matrix_ks(1),&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_res(ispin)%matrix, name = "RESIDUE", &
                         template = ls_env%matrix_ks(1),&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_tmp(ispin)%matrix, name = "TMP MAT", &
                         template = ls_env%matrix_ks(1),&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_tmp1(ispin)%matrix, name = "TMP1 MAT", &
                         template = ls_env%matrix_ks(1),&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_PKs(ispin)%matrix, name = "PKs MAT", &
                         template = ls_env%matrix_ks(1),&
                         matrix_type=dbcsr_type_no_symmetry)
       CALL dbcsr_create(matrix_z0(ispin)%matrix, name = "PKs MAT", &
                         template = ls_env%matrix_ks(1),&
                         matrix_type=dbcsr_type_no_symmetry)
    END DO

!    CALL commutator(matrix_p, ls_env%matrix_ks, matrix_cg, eps_filter, .TRUE., -2.0_dp, 0.0_dp)
!    DO ispin = 1, nspins
!       ! Residual matrix in first step=cg matrix.
!       CALL dbcsr_copy(matrix_res(ispin,1)%matrix, matrix_cg(ispin,1)%matrix)
!       ! Begin the real optimization loop
!       CALL dbcsr_set(matrix_cg_z(ispin,1)%matrix, 0.0_dp)
!    END DO
!
!    ! Precompute -FP-[FP]T which will be used throughout the CG iterations
!    ! b = -1.0_dp * ( (PF) + (PF)^T )
!    CALL commutator(matrix_p, ls_env%matrix_ks, matrix_b, eps_filter, .FALSE., -1.0_dp, 0.0_dp)

    DO ispin = 1, nspins

! Precompute some matrices to save work during iterations
       CALL dbcsr_multiply("N", "N", 1.0_dp, ls_env%matrix_p(ispin), ls_env%matrix_ks(ispin), &
                                0.0_dp, matrix_PKs(ispin)%matrix, filter_eps=ls_env%eps_filter)
       CALL dbcsr_transposed(matrix_b(ispin)%matrix, matrix_PKs(ispin)%matrix)
       CALL dbcsr_copy(matrix_cg(ispin)%matrix, matrix_b(ispin)%matrix)

! Starting CG with guess 0-matrix gives -2*gradient=[Ks*P-(Ks*P)T] for cg_matrix in second step
       CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_PKs(ispin)%matrix, 2.0_dp, -2.0_dp)

! Residual matrix in first step=cg matrix. Keep Pks for later use in CG!
       CALL dbcsr_copy(matrix_res(ispin)%matrix, matrix_cg(ispin)%matrix)

! Precompute -FP-[FP]T which will be used throughout the CG iterations
       CALL dbcsr_add(matrix_b(ispin)%matrix, matrix_PKs(ispin)%matrix, -1.0_dp, -1.0_dp)

       ! Ax0
       CALL commutator_symm(matrix_b(ispin)%matrix, matrix_cg(ispin)%matrix, &
                            matrix_Ax(ispin)%matrix, ls_env%eps_filter, 1.0_dp)
       CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_cg(ispin)%matrix, ls_env%matrix_p(ispin), &
                           0.0_dp, matrix_tmp(ispin)%matrix, filter_eps=ls_env%eps_filter)
       CALL commutator_symm(ls_env%matrix_ks(ispin), matrix_tmp(ispin)%matrix, &
                            matrix_tmp1(ispin)%matrix, ls_env%eps_filter, 2.0_dp)
       CALL dbcsr_add(matrix_Ax(ispin)%matrix, matrix_tmp1(ispin)%matrix, 1.0_dp, 1.0_dp)

    END DO

    ! r0 = b - Ax0
    DO ispin = 1, nspins
       CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, 1.0_dp, -1.0_dp)
    END DO

    DO ispin = 1, nspins
       CALL one_projector(qs_env, ls_env%matrix_p(ispin), matrix_res(ispin)%matrix, ls_env%eps_filter)
    END DO

    do_precond = .FALSE.
    ! Preconditioner
    IF (do_precond) THEN
       ! z_0 = M * r_0
       ! Conjugate gradient returns z_0
       ! IF(para_env%mepos==0) WRITE(*,*) "M * z_0 = r_0"
!       CALL preconditioner(qs_env = qs_env,&
!                           matrix_ks = ls_env%matrix_ks, &
!                           matrix_p = ls_env%matrix_p, &
!                           matrix_rhs = matrix_res, &
!                           matrix_cg_z = matrix_z0,&
!                           eps_filter = els_env%ps_filter,&
!                           iounit = iounit)

!       CALL precond_lvl3(qs_env = qs_env,&
!                         matrix_ks = ls_env%matrix_ks, &
!                         matrix_p = ls_env%matrix_p, &
!                         matrix_rhs = matrix_res, &
!                         matrix_cg_z = matrix_z0,&
!                         eps_filter = els_env%ps_filter,&
!                         iounit = iounit)
    ELSE
       ! z_0 = r_0
       DO ispin = 1, nspins
          CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
       END DO
    END IF

    DO ispin = 1, nspins

       ! cg = p_0 = z_0
       CALL dbcsr_copy(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix)

       ! Tr(r_0 * z_0)
       CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_z0(ispin)%matrix, norm_rr(ispin))
       IF (iounit > 0) WRITE(iounit,*) "Tr(r_0 * z_0) (norm_rr)", norm_rr(1)

       IF (norm_rr(ispin) .LT. 0.0_dp) CPABORT("norm_rr < 0")
       norm_res = MAX(norm_res, ABS(norm_rr(ispin)/REAL(nao,dp)))
       IF(iounit > 0) WRITE(iounit,*) "norm_res", norm_res
    END DO

    ! header
    IF (iounit > 0) THEN
       WRITE (iounit, "(/,T3,A,T16,A,T25,A,T38,A,T52,A,/,T3,A)") &
          "Iteration", "Method", "Stepsize", "Convergence", "Time", &
          REPEAT("-", 80)
    ENDIF

    alpha(:) = 0.0_dp
    converged = .FALSE.
    restart = .FALSE.
    flag= "HAR"

    ! start iteration
    iteration: DO i = 1, max_iter

       ! default for eps 10E-6 in MO_linres
       IF (norm_res .LT. ec_env%eps_default) THEN
          converged = .TRUE.
       ENDIF

       t2 = m_walltime()
       IF (i .EQ. 1 .OR. MOD(i, 1) .EQ. 0 .OR. converged ) THEN
          IF (iounit > 0) THEN
             !WRITE (iounit, "(T5,I5,T18,A3,T28,L1,T38,1E8.2,T48,F16.10,T68,F8.2)") &
             WRITE (iounit, "(T5,I5,T16,A6,T28,L1,T38,1E8.2,T48,F16.10,T68,F8.2)") &
                i, flag, restart, MAXVAL(alpha), norm_res, t2 - t1
!              IF (iounit > 0) WRITE (iounit, "(T2,A18,T20,I3,T25,1E16.10,T45,1E16.10,T68,1E16.10)") &
!                                     "It,alp,norm,time", i, MAXVAL(alpha), norm_res, t2-t1
             CALL m_flush(iounit)
          ENDIF
       ENDIF
       IF (converged) THEN
          IF (iounit > 0) THEN
             WRITE (iounit, "(/,T2,A,I4,A,/)") "The linear solver converged in ", i, " iterations."
             CALL m_flush(iounit)
          ENDIF
          EXIT iteration
       ENDIF

       ! Max number of iteration reached
       IF (i == max_iter) THEN
          IF (iounit > 0) THEN
             WRITE (iounit, "(/,T2,A/)") &
                "The linear solver didnt converge! Maximum number of iterations reached."
             CALL m_flush(iounit)
          ENDIF
          converged = .FALSE.
       ENDIF

!       ! One step to compute: -FPD-DPF-DFP-PFD (not obvious but symmetry allows for some tricks)
!       CALL commutator(matrix_b, matrix_cg, matrix_Ax, ls_env%eps_filter, .TRUE., 1.0_dp, 0.0_dp)
!       ! Compute the missing bits 2*(FDP+PDF) (again use symmetry to compute as a commutator)
!       DO ispin = 1, nspins
!          CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_cg(ispin,1)%matrix, &
!                              matrix_p(ispin,1)%matrix, 0.0_dp, matrix_tmp(ispin,1)%matrix, &
!                              filter_eps=ls_env%eps_filter)
!       END DO
!       CALL commutator(matrix_ks, matrix_tmp, matrix_Ax, ls_env%eps_filter, .TRUE., 2.0_dp, 1.0_dp)

       DO ispin = 1, nspins

          ! Ax
          !------------------
          CALL commutator_symm(matrix_b(ispin)%matrix, matrix_cg(ispin)%matrix, &
                               matrix_Ax(ispin)%matrix, ls_env%eps_filter, 1.0_dp)
          CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_cg(ispin)%matrix, ls_env%matrix_p(ispin), &
                              0.0_dp, matrix_tmp(ispin)%matrix, filter_eps=ls_env%eps_filter)
          CALL commutator_symm(ls_env%matrix_ks(ispin), matrix_tmp(ispin)%matrix, &
                               matrix_tmp1(ispin)%matrix, ls_env%eps_filter, 2.0_dp)
          CALL dbcsr_add(matrix_Ax(ispin)%matrix, matrix_tmp1(ispin)%matrix, 1.0_dp, 1.0_dp)
       END DO
          !------------------

       DO ispin = 1, nspins
          CALL one_projector(qs_env, ls_env%matrix_p(ispin), matrix_Ax(ispin)%matrix, ls_env%eps_filter)

!          CALL dbcsr_add(matrix_Ax(ispin,1)%matrix, matrix_cg(ispin,1)%matrix, 1.0_dp, shift)

          ! norm_cA = tr(Ap_j * p_j)
          CALL dbcsr_dot(matrix_cg(ispin)%matrix, matrix_Ax(ispin)%matrix, norm_cA(ispin))
          IF (iounit > 0) WRITE (iounit,*)  "norm_cA", norm_cA(ispin)

          ! Determine step-size
          ! alpha = norm_rr/norm_cA
          IF (norm_cA(ispin) .LT. ec_env%eps_default) THEN
             alpha(ispin) = 1.0_dp
          ELSE
             alpha(ispin) = norm_rr(ispin)/norm_cA(ispin)
          END IF

          ! x_j+1 = x_j + alpha*p_j
          ! save response-denisty of this iteration
          CALL dbcsr_add(matrix_cg_z(ispin)%matrix, matrix_cg(ispin)%matrix, &
                         1.0_dp, alpha(ispin))

          ! r_j+1 = r_j - alpha * Ap_j
          CALL dbcsr_add(matrix_res(ispin)%matrix, matrix_Ax(ispin)%matrix, &
                         1.0_dp, -alpha(ispin))

          IF (do_precond) THEN
          ELSE
             ! z_j+1 = r_j+1
             CALL dbcsr_copy(matrix_z0(ispin)%matrix, matrix_res(ispin)%matrix)
          END IF

       END DO
      
       norm_res = 0.0_dp
 
       DO ispin = 1, nspins
 
          ! Tr[r_j+1*z_j+1]
          CALL dbcsr_dot(matrix_res(ispin)%matrix, matrix_z0(ispin)%matrix, new_norm(ispin))
          IF (new_norm(ispin) .LT. 0.0_dp) CPABORT("tr(r_j+1*z_j+1) < 0")
          IF (iounit > 0) WRITE (iounit,*)  "new_norm", new_norm(1)
          norm_res = MAX(norm_res, new_norm(ispin)/REAL(nao,dp))
          IF (iounit > 0) WRITE (iounit,*)  "norm_res", norm_res
 
          IF (norm_rr(ispin) .LT. ec_env%eps_default*0.001_dp &
              .OR. new_norm(ispin) .LT. ec_env%eps_default*0.001_dp) THEN
             beta(ispin) = 0.0_dp
             converged = .TRUE.
             IF (iounit > 0) WRITE (iounit,*)  "beta = 0.0 CONVERGED"
          ELSE
             beta(ispin) = new_norm(ispin)/norm_rr(ispin)
          END IF
          IF (iounit > 0) WRITE (iounit,*)  "beta", beta(1)
 
          ! update new search vector (matrix cg)
          ! Here: cg_j+1 = z_j+1 + beta*cg_j
          CALL dbcsr_add(matrix_cg(ispin)%matrix, matrix_z0(ispin)%matrix, beta(ispin), 1.0_dp)
          CALL dbcsr_filter(matrix_cg(ispin)%matrix, ls_env%eps_filter)
 
          norm_rr(ispin) = new_norm(ispin)
       END DO
 
    END DO iteration

    DO ispin = 1, nspins
       CALL dbcsr_set(matrix_Ax(ispin)%matrix, 0.0_dp)
    END DO
 
    DO ispin =1, nspins
       CALL one_projector(qs_env, ls_env%matrix_p(ispin), matrix_cg_z(ispin)%matrix, ls_env%eps_filter)
    END DO

    step_size = 1.0_dp

    norm_pao = dbcsr_frobenius_norm(matrix_cg_z(1)%matrix)
    IF (iounit > 0) WRITE (iounit,*)  "||CGZ||_F :", norm_pao
    IF (iounit > 0) WRITE (iounit,*) "Before update_p_exp"
    norm_pao = dbcsr_frobenius_norm(ls_env%matrix_p(1))
    IF (iounit > 0) WRITE (iounit,*) "||CGZ||_F :", norm_pao
    norm_pao = dbcsr_frobenius_norm(ls_env%matrix_p(1))
    IF (iounit > 0) WRITE (iounit,*) "||P(AO)||_F :", norm_pao

    ! P_IN : curvy_env%matrix_p
    ! P_out : ls_sc_env%matrix_p

!    WRITE(*,*) "p_in"
!    CALL dbcsr_print(p_in(1)%matrix)
!    WRITE(*,*) "p_out"
!    CALL dbcsr_print(p_out(1)%matrix)

    NULLIFY(matrix_BCH)
    CALL dbcsr_allocate_matrix_set(matrix_BCH, nspins, 1)
    DO ispin = 1, nspins
       ALLOCATE (matrix_BCH(ispin, 1)%matrix)
       CALL dbcsr_create(matrix_BCH(ispin, 1)%matrix, name = "BCH dummy", &
                         template = ls_env%matrix_s)
       ! Local matrix KS - ortho transformation destroys symmetric type
       CALL dbcsr_copy(matrix_BCH(ispin, 1)%matrix, ls_env%matrix_s)
       CALL dbcsr_set(matrix_BCH(ispin, 1)%matrix, 0.0_dp)
    END DO

    ! CG_Z -> P_out
    CALL update_p_exp(matrix_p_in = matrix_p_in, &
                      matrix_p_out = ls_env%matrix_p, &
                      matrix_dp = matrix_cg_z, &
                      matrix_BCH = matrix_BCH, &
                      threshold = ls_env%eps_filter, &
                      step_size = step_size)

    IF (iounit > 0) WRITE (iounit,*) "After update_p_exp"
    norm_pao = dbcsr_frobenius_norm(ls_env%matrix_p(1))
    IF (para_env%mepos==0) WRITE(*,*) "||P||_F :", norm_pao
 
    DO ispin = 1, nspins
       ! Transform P_out back into non-orthogonal basis
       CALL transform_m_orth(ls_env%matrix_p(ispin), ls_env%matrix_s_sqrt_inv, ls_env%eps_filter)
 
       ! bring back matrix_ks to non-ortho
       CALL transform_m_orth(ls_env%matrix_ks(ispin), ls_env%matrix_s_sqrt, ls_env%eps_filter)
    END DO

    ! Scale matrix_p by factor 1/2 in closed-shell
    IF (nspins == 1) CALL dbcsr_scale(ls_env%matrix_p(1), 2.0_dp)

    norm_pao = dbcsr_frobenius_norm(ls_env%matrix_p(1))
    IF (iounit > 0) WRITE (iounit,*) "||P(AO)||_F :", norm_pao

    ! Release matrices

    CALL dbcsr_deallocate_matrix_set(matrix_Ax)
    CALL dbcsr_deallocate_matrix_set(matrix_b)
    CALL dbcsr_deallocate_matrix_set(matrix_cg)
    CALL dbcsr_deallocate_matrix_set(matrix_cg_z)
    CALL dbcsr_deallocate_matrix_set(matrix_p_in)
    CALL dbcsr_deallocate_matrix_set(matrix_res)
    CALL dbcsr_deallocate_matrix_set(matrix_tmp)
    CALL dbcsr_deallocate_matrix_set(matrix_tmp1)
    CALL dbcsr_deallocate_matrix_set(matrix_PKs)
    CALL dbcsr_deallocate_matrix_set(matrix_z0)

    CALL dbcsr_deallocate_matrix_set(matrix_BCH)

    DEALLOCATE(alpha, beta, new_norm, norm_cA, norm_rr)

    IF (iounit > 0) WRITE (iounit,*) " ---- newton_solver - END ----"

    CALL timestop(handle)

  END SUBROUTINE prepare_newton 

! **************************************************************************************************
!> \brief computes a commutator exploiting symmetry RES=k*[A,B]=k*[AB-(AB)T]
!> \param a ...
!> \param b ...
!> \param res ...
!> \param eps_filter   filtering threshold for sparse matrices
!> \param prefac      prefactor k in above equation
!> \par History
!>       2012.05 created [Florian Schiffmann]
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE commutator_symm(a, b, res, eps_filter, prefac)
      TYPE(dbcsr_type)                                   :: a, b, res
      REAL(KIND=dp)                                      :: eps_filter, prefac

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'commutator_symm'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: work

      CALL timeset(routineN, handle)

      CALL dbcsr_create(work, template=a, matrix_type=dbcsr_type_no_symmetry)

      CALL dbcsr_multiply("N", "N", prefac, a, b, 0.0_dp, res, filter_eps=eps_filter)
      CALL dbcsr_transposed(work, res)
      CALL dbcsr_add(res, work, 1.0_dp, -1.0_dp)

      CALL dbcsr_release(work)

      CALL timestop(handle)
   END SUBROUTINE commutator_symm

! **************************************************************************************************
!> \brief computes (anti-)commutator exploiting (anti-)symmetry:
!>        A symmetric : RES = beta*RES + k*[A,B] = k*(AB-(AB)^T)
!>        A anti-sym  : RES = beta*RES + k*{A,B} = k*(AB+(AB)^T)
!>
!>
!> \param a          Matrix A
!> \param b          Matrix B
!> \param res        Commutator result
!> \param eps_filter filtering threshold for sparse matrices
!> \param alpha      Scaling of anti-/commutator
!> \param beta       Scaling of contens of inital contens of result matrix
!> \param anticomm   Calculate anticommutator
!>
!> \par History
!>       2020.07 Fabian Belleflamme  (based on commutator_symm)
! **************************************************************************************************
   SUBROUTINE commutator(a, b, res, eps_filter, anticomm, alpha, beta)

      TYPE(dbcsr_p_type), DIMENSION(:,:)                 :: a, b, res
      REAL(KIND=dp)                                      :: eps_filter
      LOGICAL                                            :: anticomm
      REAL(KIND=dp), OPTIONAL                            :: alpha, beta 

      CHARACTER(LEN=*), PARAMETER :: routineN = 'commutator_symm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      REAL(KIND=dp)                                      :: facc, myalpha, mybeta
      TYPE(dbcsr_type)                                   :: work, work2

      CALL timeset(routineN, handle)

      CALL dbcsr_create(work, template=a(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(work2, template=a(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)

      myalpha = 1.0_dp
      IF (PRESENT(alpha)) myalpha = alpha
      mybeta = 0.0_dp
      IF (PRESENT(beta)) mybeta = beta
      facc = -1.0_dp
      IF (anticomm) facc = 1.0_dp

      DO ispin = 1, SIZE(a) 

         CALL dbcsr_multiply("N", "N", myalpha, a(ispin, 1)%matrix, b(ispin, 1)%matrix, &
                             0.0_dp, work, filter_eps=eps_filter)
         CALL dbcsr_transposed(work2, work)

         ! RES= beta*RES + alpha*{A,B} = beta*RES + alpha*[AB+(AB)T]
         ! RES= beta*RES + alpha*[A,B] = beta*RES + alpha*[AB-(AB)T]
         CALL dbcsr_add(work, work2, 1.0_dp, facc)

         CALL dbcsr_add(res(ispin, 1)%matrix, work, mybeta, 1.0_dp)

      END DO 

      CALL dbcsr_release(work)
      CALL dbcsr_release(work2)

      CALL timestop(handle)

   END SUBROUTINE commutator


! **************************************************************************************************
!> \brief Projector P(M) = P*M*Q^T + Q*M*P^T
!>        with P = D
!>        with Q = (1-D)
!> 
!> \param qs_env ...
!> \param matrix_p  Ground-state density in orthonormal basis 
!> \param matrix_io Matrix to which projector is applied.
!>          
!> \date    06.2020
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE one_projector(qs_env, matrix_p, matrix_io, eps_filter)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(dbcsr_type), INTENT(IN), POINTER              :: matrix_p
    TYPE(dbcsr_type), INTENT(INOUT), POINTER           :: matrix_io
    REAL(KIND=dp), INTENT(IN)                          :: eps_filter

    CHARACTER(len=*), PARAMETER :: routineN = 'one_projector', routineP = moduleN//':'//routineN

    INTEGER                                            :: handle
    TYPE(dft_control_type), POINTER                    :: dft_control
    TYPE(dbcsr_type)                                   :: matrix_tmp, matrix_q

    CALL timeset(routineN, handle)

!    IF(para_env%mepos==0) WRITE(*,*) "projector BEGINNING"
    CALL get_qs_env(qs_env, dft_control=dft_control)

    CALL dbcsr_create(matrix_q, template = matrix_p, &
                      matrix_type = dbcsr_type_no_symmetry)
    CALL dbcsr_create(matrix_tmp, template = matrix_p, &
                      matrix_type = dbcsr_type_no_symmetry)

    ! Q = (1 - P)
    CALL dbcsr_copy(matrix_q, matrix_p)
    CALL dbcsr_scale(matrix_q, -1.0_dp)
    CALL dbcsr_add_on_diag(matrix_q, 1.0_dp)
    CALL dbcsr_finalize(matrix_q)

    ! Proj(M) = P*M*Q + Q*M*P
    ! with P = D = CC^T
    ! and  Q = (1 - P)

    ! tmp1 = P*M
    CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p, matrix_io, &
                        0.0_dp, matrix_tmp, filter_eps = eps_filter)
    ! m_io = P*M*Q
    CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, matrix_q, &
                        0.0_dp, matrix_io, filter_eps = eps_filter)

    ! tmp = (P^T*M^T*Q^T)^T = -(P*M*Q)^T
    CALL dbcsr_transposed(matrix_tmp, matrix_io)
    CALL dbcsr_add(matrix_io, matrix_tmp, 1.0_dp, -1.0_dp)

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

! **************************************************************************************************
!> \brief Use the BCH update to get the new idempotent P
!>        Numerics don't allow for perfect idempotency, therefore a mc weeny
!>        step is used to make sure we stay close to the idempotent surface
!> \param matrix_p_in ...
!> \param matrix_p_out ...
!> \param matrix_dp ...
!> \param matrix_BCH ...
!> \param threshold ...
!> \param step_size ...
!> \param BCH_saved ...
!> \param n_bch_hist ...
!> \par History
!>       2012.05 created [Florian Schiffmann]
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE update_p_exp(matrix_p_in, matrix_p_out, matrix_dp, matrix_BCH, threshold, step_size)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p_in
      TYPE(dbcsr_type), DIMENSION(:)                     :: matrix_p_out
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_dp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_BCH
      REAL(KIND=dp)                                      :: threshold, step_size(2)

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_p_exp'

      INTEGER                                            :: handle, i, ispin, nsave, nspin, unit_nr
      LOGICAL                                            :: save_BCH
      REAL(KIND=dp)                                      :: frob_norm, step_fac
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix, matrix_tmp
      INTEGER                                            :: BCH_saved(2), n_bch_hist
      REAL(KIND=dp)                                      :: norm_res 

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF
      WRITE(*,*) "update_p_exp - BEGINNING"

      CALL dbcsr_create(matrix, template=matrix_p_in(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp, template=matrix_p_in(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      nspin = SIZE(matrix_p_in,1)

      n_bch_hist = 0
      BCH_saved = 0

      DO ispin = 1, nspin
         step_fac = 1.0_dp
         frob_norm = 1.0_dp
         nsave = 0

         CALL dbcsr_copy(matrix_tmp, matrix_p_in(ispin)%matrix)
         CALL dbcsr_copy(matrix_p_out(ispin), matrix_p_in(ispin)%matrix)
! If a BCH history is used make good use of it and do a few steps as a copy and scale update of P
! else BCH_saved will be 0 and loop is skipped
WRITE(*,*) "Start reading BCH history"
         DO i = 1, BCH_saved(ispin)
            WRITE(*,*) "iterate over BCH_saved",i
            step_fac = step_fac*step_size(ispin)
            CALL dbcsr_copy(matrix_tmp, matrix_p_out(ispin))
            CALL dbcsr_add(matrix_p_out(ispin), matrix_BCH(ispin, i)%matrix, 1.0_dp, ifac(i)*step_fac)
            CALL dbcsr_add(matrix_tmp, matrix_p_out(ispin), 1.0_dp, -1.0_dp)
            frob_norm = dbcsr_frobenius_norm(matrix_tmp)
            IF (unit_nr .GT. 0) WRITE (unit_nr, "(t3,a,i3,a,f16.8)") "BCH: step", i, " Norm of P_old-Pnew:", frob_norm
            IF (frob_norm .LT. threshold) EXIT
         END DO

WRITE(*,*) "Stop reading BCH history"
WRITE(*,*) "BCH_saved(1): ", BCH_saved(1)
WRITE(*,*) "BCH_saved(2): ", BCH_saved(2)
WRITE(*,*) "Start with BCH loop"
!            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, matrix_dp(1)%matrix, &
!                                0.0_dp, matrix, filter_eps=threshold)

!WRITE(*,*) "matrix"
!CALL dbcsr_print(matrix)
!WRITE(*,*) "matrix_tmp"
!CALL dbcsr_print(matrix_tmp)
!WRITE(*,*) "matrix_dp"
!CALL dbcsr_print(matrix_dp(1)%matrix)

         IF (frob_norm .LT. threshold) CYCLE

! If the copy and scale isn't enough compute a few more BCH steps. 20 seems high but except of the first step it will never be close
         save_BCH = BCH_saved(ispin) == 0 .AND. n_bch_hist .GT. 0
         DO i = BCH_saved(ispin) + 1, 20
WRITE(*,*) "BCH loop - i", i
            step_fac = step_fac*step_size(ispin)
WRITE(*,*) "step_fac", step_fac
WRITE(*,*) "step_size(ispin)", step_size(ispin)
            !allow for a bit of matrix magic here by exploiting matrix and matrix_tmp
            !matrix_tmp is alway the previous order of the BCH series
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_tmp, matrix_dp(ispin)%matrix, &
                                0.0_dp, matrix, filter_eps=threshold)
WRITE(*,*) "BCH loop - after multiply"
            norm_res = dbcsr_frobenius_norm(matrix)
            WRITE(*,*) "update_p_exp, matrix", norm_res

            !(anti)symmetry allows to sum the transposed instead of the full commutator, matrix becomes the latest result

            CALL dbcsr_transposed(matrix_tmp, matrix)
            CALL dbcsr_add(matrix, matrix_tmp, 1.0_dp, 1.0_dp)

            !Finally, add the new BCH order to P, but store the previous one for a convergence check
            CALL dbcsr_copy(matrix_tmp, matrix_p_out(ispin))
            CALL dbcsr_add(matrix_p_out(ispin), matrix, 1.0_dp, ifac(i)*step_fac)
            IF (save_BCH .AND. i .LE. n_bch_hist) THEN
               CALL dbcsr_copy(matrix_BCH(ispin,1)%matrix, matrix)
               nsave = i
            END IF

            CALL dbcsr_add(matrix_tmp, matrix_p_out(ispin), 1.0_dp, -1.0_dp)

            !Stop the BCH-series if two successive P's differ by less the threshold
            frob_norm = dbcsr_frobenius_norm(matrix_tmp)
            IF (unit_nr .GT. 0) WRITE (unit_nr, "(t3,a,i3,a,f16.8)") "BCH: step", i, " Norm of P_old-Pnew:", frob_norm
            IF (frob_norm .LT. threshold) EXIT

            !Copy the latest BCH-matrix on matrix tmp, so we can cycle with all matrices in place
            CALL dbcsr_copy(matrix_tmp, matrix)
            CALL dbcsr_filter(matrix_tmp, threshold)
         END DO
         BCH_saved(ispin) = nsave
         IF (unit_nr .GT. 0) WRITE (unit_nr, "(A)") " "
      END DO

!      CALL purify_mcweeny(matrix_p_out, threshold, 1)
      IF (unit_nr .GT. 0) CALL m_flush(unit_nr)
      CALL dbcsr_release(matrix_tmp)
      CALL dbcsr_release(matrix)
      CALL timestop(handle)

   END SUBROUTINE update_p_exp

! **************************************************************************************************
!> \brief McWeeny purification of a matrix in the orthonormal basis
!> \param matrix_p Matrix to purify (needs to be almost idempotent already)
!> \param threshold Threshold used as filter_eps and convergence criteria
!> \param max_steps Max number of iterations
!> \par History
!>       2013.01 created [Florian Schiffmann]
!>       2014.07 slightly refactored [Ole Schuett]
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE purify_mcweeny(matrix_p, threshold, max_steps)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p
      REAL(KIND=dp)                                      :: threshold
      INTEGER                                            :: max_steps

      CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mcweeny'

      INTEGER                                            :: handle, i, ispin, unit_nr
      REAL(KIND=dp)                                      :: frob_norm, trace
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_type)                                   :: matrix_pp, matrix_tmp

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      CALL dbcsr_create(matrix_pp, template=matrix_p(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(matrix_tmp, template=matrix_p(1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_trace(matrix_p(1)%matrix, trace)

      DO ispin = 1, SIZE(matrix_p)
         DO i = 1, max_steps
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_p(ispin)%matrix, matrix_p(ispin)%matrix, &
                                0.0_dp, matrix_pp, filter_eps=threshold)

            ! test convergence
            CALL dbcsr_copy(matrix_tmp, matrix_pp)
            CALL dbcsr_add(matrix_tmp, matrix_p(ispin)%matrix, 1.0_dp, -1.0_dp)
            frob_norm = dbcsr_frobenius_norm(matrix_tmp) ! tmp = PP - P
            IF (unit_nr > 0) WRITE (unit_nr, '(t3,a,f16.8)') "McWeeny: Deviation of idempotency", frob_norm
            IF (unit_nr > 0) CALL m_flush(unit_nr)

            ! construct new P
            CALL dbcsr_copy(matrix_tmp, matrix_pp)
            CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_pp, matrix_p(ispin)%matrix, &
                                3.0_dp, matrix_tmp, filter_eps=threshold)
            CALL dbcsr_copy(matrix_p(ispin)%matrix, matrix_tmp) ! tmp = 3PP - 2PPP

            ! frob_norm < SQRT(trace*threshold)
            IF (frob_norm*frob_norm < trace*threshold) EXIT
         END DO
      END DO

      CALL dbcsr_release(matrix_pp)
      CALL dbcsr_release(matrix_tmp)
      CALL timestop(handle)

   END SUBROUTINE purify_mcweeny

END MODULE ec_newton
