!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Harris functional force driver routine
!>
!> Force contributions : ec_core_hamiltonian:  kinetic, overlap, gth_ppl, 
!>                                             gth_ppnl
!>                       ec_build_ks_matrix:   core_overlap, rho_core, rho_elec
!>                       ec_nsc_force:         (-1)*rho_elec
!> F_tot = F_kin[n_out] + F_gth_ppl[n_out] + F_gth_ppnl[n_out] 
!>         + F_BS_H[n_in, n_out] + F_BS_XC[n_in, n_out] 
!>         + F_core[rho_core + n_in] + F_NSC[n_in, n_out]
!>
!> w/ F_NSCF_NSC[n_in, n_out] = F_BS_H + F_BS_XC + F_delta_H + F_delta_XC
!>                            = \int (v_hartree[rho_out-rho_in] +
!>                                   (rho_out - rho_in)*\mu'_[rho_in])*\nabla_I(rho_in) dr
!>
!> \date 02.2018
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_force
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix_set,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_print,&
                                              dbcsr_release,&
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
                                              ec_w_matrix
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
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_force'

! *** Public subroutines ***

   PUBLIC :: ec_forces, ec_bs_forces, ec_nsc_force

CONTAINS


! **************************************************************************************************
!> \brief   Calculate the forces in the Harris energy correction.
!>          First, saves rho to ec_env%old_rho
!>          Then, updates rho with new ec_env%matrix_p
!>          Then, calls routines for force calculation
!>          F_tot = F_band + F_II + F_NSC_XC + F_NSC_H
!> \param qs_env ...
!> \date    02.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_forces(qs_env)

    TYPE(qs_environment_type), POINTER                 :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'ec_forces', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle
    LOGICAL                                            :: calculate_forces 
    TYPE(cp_logger_type), POINTER                      :: logger
    TYPE(energy_correction_type), POINTER              :: ec_env
    TYPE(kg_environment_type), POINTER                 :: kg_env
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
    TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
    TYPE(qs_rho_type), POINTER                         :: rho
    TYPE(dft_control_type), POINTER                    :: dft_control
    INTEGER                                            :: ispin, nspins, nkind
!    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r
!    TYPE(pw_p_type), DIMENSION(:), POINTER             :: old_rho_r, old_rho_g
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: old_rho_ao
    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: natom_of_kind
    TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
!fbelle
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: old_rho_r
!    TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_w
    REAL(KIND=dp)                                      :: dia
    REAL(KIND=dp)                                      :: doa
    REAL(KIND=dp)                                      :: delta_rho
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_dif_rspace

    TYPE(pw_env_type), POINTER                         :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
    TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
    REAL(KIND=dp)                                      :: norm_scalar
    TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_delta

!!---------------------------------------------------------------------------------------
!!  Initial general setup
!!---------------------------------------------------------------------------------------

    CALL timeset(routineN, handle)
    NULLIFY (logger)
    logger => cp_get_default_logger()
    
    
    CALL get_qs_env(qs_env,para_env=para_env)
    IF(para_env%mepos==0) THEN
!       WRITE(*,*) "ec_forces - BEGINNING"
    END IF

    ! obviously 
    calculate_forces = .TRUE.
    
    ! get qs and kg environment
    NULLIFY (atomic_kind_set, para_env, rho_ao, old_rho_ao, dft_control)
    CALL get_qs_env(qs_env = qs_env,& 
                    kg_env = kg_env,&
                    atomic_kind_set = atomic_kind_set,&
                    dft_control = dft_control,&
                    force = force,&
                    para_env = para_env,&
                    particle_set = particle_set,&
                    rho = rho)

    ! shortcut to ec_type
    !  ec_build_neighborlist has been built during kg_energy_correction call
    ec_env => kg_env%ec_env
    nspins = dft_control%nspins

!!---------------------------------------------------------------------------------------
!!  storage force environment
!!---------------------------------------------------------------------------------------

    ! Introduce new ec_env%force type to store forces and prevent overwriting
    !    especially to rho_elec
    IF (.NOT. ASSOCIATED(ec_env%force)) THEN
       ! get atom mapping
       nkind = SIZE(atomic_kind_set)
       ALLOCATE (natom_of_kind(nkind))
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                natom_of_kind=natom_of_kind)
       
       CALL allocate_qs_force(ec_env%force, natom_of_kind)
       DEALLOCATE (natom_of_kind)
    END IF
    CALL zero_qs_force(ec_env%force)

!!---------------------------------------------------------------------------------------
!!  storage old density
!!---------------------------------------------------------------------------------------

    ! allocate and initialize the old density
    CALL qs_rho_create(ec_env%old_rho)
    ! copy rho structure, alert, sets density matrix to zero
    CALL duplicate_rho_type(rho_input = rho, &
                            rho_output = ec_env%old_rho, &
                            qs_env = qs_env)

    CALL qs_rho_get(rho, rho_ao_kp = rho_ao )
    CALL qs_rho_get(ec_env%old_rho, rho_ao_kp = old_rho_ao )
    ! copy density matrix manually  
    DO ispin = 1, nspins
       CALL dbcsr_copy(old_rho_ao(ispin,1)%matrix, & 
                         rho_ao(ispin,1)%matrix, keep_sparsity=.FALSE.)
    END DO
    
    ! update rho strucutre with new density matrix of energy correction
    ! copy ec_env%matrix_p to rho%rho_ao(:,:) 
    DO ispin = 1, nspins
       CALL dbcsr_copy(rho_ao(ispin,1)%matrix, ec_env%matrix_p(ispin,1)%matrix, & 
                         keep_sparsity = .FALSE.)
    END DO
    ! updates rho_r and rho_g wrt to new rho%rho_ao
    CALL qs_rho_update_rho(rho, qs_env = qs_env, task_list_external = ec_env%task_list)




!!---------------------------------------------------------------------------------------
!!---------------------------------------------------------------------------------------
!!  Printing of density and integrate over grid
!!---------------------------------------------------------------------------------------
!!---------------------------------------------------------------------------------------

    ! get rho_out - after Harris energy correction
    CALL qs_rho_get(rho, & 
                    rho_r = rho_r)

    ! get rho_in - before Harris energy correction
    CALL qs_rho_get(ec_env%old_rho, &
                    rho_r = old_rho_r)

    ! density integrals over grid: \int rho dr = N_el
    dia = pw_integrate_function(old_rho_r(1)%pw)
    doa = pw_integrate_function(rho_r(1)%pw)

    NULLIFY (pw_env)
    CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
    CPASSERT(ASSOCIATED(pw_env))

    NULLIFY (auxbas_pw_pool, pw_pools)
    ! gets the tmp grids
    CALL pw_env_get(pw_env, &
                    auxbas_pw_pool = auxbas_pw_pool, &
                    pw_pools = pw_pools)

    ALLOCATE(rho_dif_rspace(nspins))
    DO ispin = 1, nspins
       NULLIFY(rho_dif_rspace(nspins)%pw)
       CALL pw_pool_create_pw(pool = auxbas_pw_pool, & 
                              pw = rho_dif_rspace(ispin)%pw, &
                              use_data = REALDATA3D, &
                              in_space = REALSPACE)
       CALL pw_zero(rho_dif_rspace(ispin)%pw)
    END DO


    ! delta density
    DO ispin = 1,nspins
      CALL pw_axpy(rho_r(ispin)%pw,rho_dif_rspace(ispin)%pw)
      CALL pw_axpy(old_rho_r(ispin)%pw, rho_dif_rspace(ispin)%pw, alpha = -1._dp)
    
      CALL pw_scale(rho_dif_rspace(ispin)%pw, rho_dif_rspace(ispin)%pw%pw_grid%dvol)
    END DO

   delta_rho = pw_integrate_function(rho_dif_rspace(1)%pw)

    CALL get_qs_env(qs_env,para_env=para_env)
    IF(para_env%mepos==0) THEN
       WRITE(*,*) " ec_forces - old_rho_r", dia
       WRITE(*,*) " ec_forces - rho_r", doa
       WRITE(*,*) " ec_forces - delta_rho", delta_rho
    END IF

    DO ispin = 1, nspins
       CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_dif_rspace(ispin)%pw)
    END DO
    DEALLOCATE(rho_dif_rspace)



    NULLIFY (matrix_delta)
    
    DO ispin = 1 , nspins
       CALL dbcsr_allocate_matrix_set(matrix_delta, nspins, 1)
       ALLOCATE (matrix_delta(ispin,1)%matrix)

       ! copy rho_ao into delta density matrix
       CALL dbcsr_copy(matrix_delta(ispin,1)%matrix, ec_env%matrix_p(ispin,1)%matrix, & 
                       keep_sparsity = .FALSE., name="DELTA DENSITY-MATRIX")
   
       ! substract old density matrix
       CALL dbcsr_add(matrix_delta(ispin, 1)%matrix, &
                      old_rho_ao(ispin, 1)%matrix, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
       
       CALL dbcsr_norm(matrix = matrix_delta(ispin, 1)%matrix, &
                       which_norm = 2, & ! dbcsr_norm_maxabsnorm
                       norm_scalar = norm_scalar)
       
       IF(para_env%mepos==0) THEN
          WRITE(*,*) "DELTA DENSITY-MATRIX"
          WRITE(*,*) "norm_scalar", norm_scalar
       END IF


    END DO


    ! Deallocate dummy ks matrix                         
    CALL dbcsr_deallocate_matrix_set(matrix_delta)

!!---------------------------------------------------------------------------------------
!!  matrix_w
!!---------------------------------------------------------------------------------------

!    IF(para_env%mepos==0) THEN
!
!       WRITE(*,*) "ec_forces - BEFORE CALCULATE W MATRIX"
!
!       WRITE(*,*) " ec_forces -  ec_env%matrix_p", ASSOCIATED(ec_env%matrix_p)
!       IF(ASSOCIATED(ec_env%matrix_p)) THEN
!             CALL dbcsr_print(ec_env%matrix_p(1,1)%matrix)
!       END IF
!       WRITE(*,*) " ec_forces - rho_ao", ASSOCIATED(rho_ao)
!       IF(ASSOCIATED(rho_ao)) THEN
!             CALL dbcsr_print(rho_ao(1,1)%matrix)
!       END IF
!       WRITE(*,*) " ec_forces - old_rho_ao", ASSOCIATED(old_rho_ao)
!       IF(ASSOCIATED(old_rho_ao)) THEN
!             CALL dbcsr_print(old_rho_ao(1,1)%matrix)
!       END IF
!    END IF


    ! calculate matrix_w = matrix_p*matrix_ks*matrix_p
    !CALL ec_calculate_w_matrix(qs_env, ec_env)
    CALL ec_w_matrix(qs_env, ec_env)

      IF(para_env%mepos==0) THEN
!       WRITE(*,*) " ec_forces - ec_env%matrix_s", ASSOCIATED(ec_env%matrix_s)
!       IF(ASSOCIATED(ec_env%matrix_s)) THEN
!             CALL dbcsr_print(ec_env%matrix_s(1,1)%matrix)
!       END IF
!       WRITE(*,*) " ec_forces - ec_env%matrix_t", ASSOCIATED(ec_env%matrix_t)
!       IF(ASSOCIATED(ec_env%matrix_t)) THEN
!             CALL dbcsr_print(ec_env%matrix_t(1,1)%matrix)
!       END IF
!       WRITE(*,*) " ec_forces - ec_env%matrix_h", ASSOCIATED(ec_env%matrix_h)
!       IF(ASSOCIATED(ec_env%matrix_h)) THEN
!             CALL dbcsr_print(ec_env%matrix_h(1,1)%matrix)
!       END IF
!       WRITE(*,*) " ec_forces -  ec_env%matrix_ks", ASSOCIATED(ec_env%matrix_ks)
!       IF(ASSOCIATED(ec_env%matrix_ks)) THEN
!             CALL dbcsr_print(ec_env%matrix_ks(1,1)%matrix)
!       END IF
!       WRITE(*,*) " ec_forces -  ec_env%matrix_p", ASSOCIATED(ec_env%matrix_p)
!       IF(ASSOCIATED(ec_env%matrix_p)) THEN
!             CALL dbcsr_print(ec_env%matrix_p(1,1)%matrix)
!       END IF
!       WRITE(*,*) " ec_forces - ec_env%matrix_w", ASSOCIATED(ec_env%matrix_w)
!       IF(ASSOCIATED(ec_env%matrix_w)) THEN
!             CALL dbcsr_print(ec_env%matrix_w(1,1)%matrix)
!       END IF
      END IF


!!---------------------------------------------------------------------------------------
!!  FORCES
!!---------------------------------------------------------------------------------------

    ! d/dR[Sum of eigenvalues] = d/dR[trace(matrix_ks[n_in] * n_out)] 
    ! core hamiltonian force contribution: 
    ! kinetic, overlap, gth_ppl, gth_ppnl 
    CALL ec_build_core_hamiltonian(qs_env, ec_env, calculate_forces)


    ! next, grid based forces of band structure: 
    ! rho_elec : F_BS_XC[n_in, n_out] + F_BS_H[n_in, n_out]
    ! rho_core : integrate_v_core_rspace [n_out]
    ! core_overlap : calculate_ecore_overlap 
!    CALL ec_build_ks_matrix(qs_env, ec_env, calculate_forces)
    CALL ec_bs_forces(qs_env, ec_env)

    ! store these forces contributions in ec_env%force
    !                  to            from
    CALL qs_copy_force(ec_env%force, force)

    ! Sweet spot - all BS forces are saved in ec_env%force, and force is about to be purged
    ! can replicate and print force right here, right now
    CALL replicate_qs_force(force, para_env)

    IF(para_env%mepos==0) THEN
       WRITE(*,*) repeat("-",80)  
       WRITE(*,*) "ec_forces - qs_env%...%force - after ec_bs_forces"
       CALL ec_write_forces(qs_env, force, ftype = 2)
!       WRITE(*,*) "ec_forces - ec_env%force - after ec_bs_forces"
!       CALL ec_write_forces(qs_env, ec_env%force)
!       WRITE(*,*) "ec_forces - ec_bs_forces start"
!       CALL ec_write_bs_forces(qs_env, force)
       WRITE(*,*) "ec_forces - ec_bs_forces end"
       WRITE(*,*) repeat("-",80)  
    END IF

    ! empty qs force type
    CALL zero_qs_force(force)
    
    ! forces from band structure and double counting terms
    ! grid based forces rho_elec
    CALL ec_nsc_force(qs_env, ec_env)

    !!---------------------------------------------------------------------------------------
    !! NSC force contributions : Hartree and XC
    !!---------------------------------------------------------------------------------------

    CALL ec_nsc_hartree_force(qs_env, ec_env)

!    IF(para_env%mepos==0) THEN
!       WRITE(*,*) repeat("-",80)  
!       WRITE(*,*) "ec_forces - qs_env%...%force - after ec_nsc_hartree_force"
!       CALL ec_write_forces(qs_env, force)
!       WRITE(*,*) "ec_forces - ec_env%force after - after ec_nsc_hartree_force "
!       CALL ec_write_forces(qs_env, ec_env%force)
!       WRITE(*,*) repeat("-",80)  
!    END IF

    CALL ec_nsc_xc_force(qs_env, ec_env)
    
!    IF(para_env%mepos==0) THEN
!       WRITE(*,*) repeat("-",80)  
!       WRITE(*,*) "ec_forces - qs_env%...%force - after ec_nsc_xc_force"
!       CALL ec_write_forces(qs_env, force)
!       WRITE(*,*) "ec_forces - ec_env%force after - after ec_nsc_xc_force "
!       CALL ec_write_forces(qs_env, ec_env%force)
!       WRITE(*,*) repeat("-",80)  
!    END IF

    !!---------------------------------------------------------------------------------------

    ! copy back forces
    CALL qs_copy_force(force, ec_env%force)

!    IF(para_env%mepos==0) THEN
!       WRITE(*,*) repeat("-",80)  
!       WRITE(*,*) "Forces END of ec_forces"  
!       WRITE(*,*) "ec_forces - qs_env%...%force - END ec_forces "
!       CALL ec_write_forces(qs_env, force)
!       WRITE(*,*) "ec_forces - ec_env%force - END ec_forces "
!       CALL ec_write_forces(qs_env, ec_env%force)
!       WRITE(*,*) repeat("-",80)  
!    END IF

    CALL timestop(handle)

  END SUBROUTINE ec_forces

! **************************************************************************************************
!> \brief  Band structure forces (Hartree and XC part)
!> \param qs_env ...
!> \param ec_env ...
!>        rho_elec : F_BS_XC[n_in, n_out] + F_BS_H[n_in, n_out]
!>        rho_core 
!>        core_overlap
!> \date    05.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_bs_forces(qs_env, ec_env)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(energy_correction_type), POINTER              :: ec_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_bs_forces', &
       routineP = moduleN//':'//routineN

    INTEGER                                            :: handle, ispin, nspins
    REAL(dp)                                           :: eexc
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    TYPE(pw_env_type), POINTER                         :: pw_env
    TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                          v_hartree_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace, v_tau_rspace
    TYPE(pw_poisson_type), POINTER                     :: poisson_env
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
    TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
    TYPE(qs_ks_env_type), POINTER                      :: ks_env
    TYPE(qs_rho_type), POINTER                         :: rho
    TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: my_matrix_ks
    LOGICAL                                            :: calculate_forces 
    TYPE(section_vals_type), POINTER                   :: input, xc_section
    REAL(dp)                                           :: ehartree
    REAL(dp)                                           :: exc
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r

    CALL timeset(routineN, handle)

    calculate_forces = .TRUE.
    
    ! get all information on the electronic density
    NULLIFY (rho, ks_env, input)
    CALL get_qs_env(qs_env = qs_env, &
                    rho = rho, &
                    dft_control = dft_control, &
                    para_env = para_env, &
                    ks_env = ks_env, &
                    input = input)

    nspins = dft_control%nspins

    NULLIFY (pw_env)
    CALL get_qs_env(qs_env = qs_env, pw_env = pw_env)
    CPASSERT(ASSOCIATED(pw_env))

    NULLIFY (auxbas_pw_pool, poisson_env, pw_pools)
    ! gets the tmp grids
    CALL pw_env_get(pw_env, &
                    auxbas_pw_pool = auxbas_pw_pool, &
                    pw_pools = pw_pools, &
                    poisson_env = poisson_env)

    ! Calculate the Hartree potential
    CALL pw_pool_create_pw(auxbas_pw_pool, &
                           v_hartree_gspace%pw, &
                           use_data = COMPLEXDATA1D, &
                           in_space = RECIPROCALSPACE)
    CALL pw_pool_create_pw(auxbas_pw_pool, &
                           rho_tot_gspace%pw, &
                           use_data = COMPLEXDATA1D, &
                           in_space = RECIPROCALSPACE)
    CALL pw_pool_create_pw(auxbas_pw_pool, &
                           v_hartree_rspace%pw, &
                           use_data = REALDATA3D, &
                           in_space = REALSPACE)

    ! input density to build potentials
    !CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, ec_env%old_rho)
    CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

    ! Hartree potential
    CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, vhartree = v_hartree_gspace%pw, ehartree = ehartree)
    CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
    CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)
    
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
   
    ! build xc potential
    NULLIFY (v_rspace, v_tau_rspace)
    NULLIFY (xc_section)

    xc_section => section_vals_get_subs_vals(input, "DFT%KG_METHOD%XC")
    
    CALL qs_vxc_create(ks_env = ks_env, &
                       !rho_struct = ec_env%old_rho, &
                       rho_struct = rho, &
                       xc_section = ec_env%xc_section, &
                       vxc_rho = v_rspace, &
                       vxc_tau = v_tau_rspace, &
                       exc = eexc, &
                       just_energy = .FALSE.)

    
    exc = pw_integrate_function(v_rspace(1)%pw)
    
    CALL qs_rho_get(rho, & 
                    rho_r = rho_r, &
                    rho_g = rho_g)

    IF(para_env%mepos==0) THEN
      WRITE(*,*) repeat("-",80) 
     
      WRITE(*,*) "ec_bs_forces - pw_poisson_solve - ehartree", ehartree
      WRITE(*,*) "ec_bs_forces - qs_vxc_create - eexc", eexc 
      WRITE(*,*) "ec_bs_forces - pw_integrate - exc", exc 

!      WRITE(*,*) "rho_r%pw%in_use", rho_r(1)%pw%in_use
!      WRITE(*,*) "rho_r(1)%pw%cr3d(1,1,1)", rho_r(1)%pw%cr3d(1,1,1)
!      WRITE(*,*) "rho_r(1)%pw%cr3d(1,1,2)", rho_r(1)%pw%cr3d(1,1,2)
!      WRITE(*,*) "rho_r(1)%pw%cr3d(1,1,3)", rho_r(1)%pw%cr3d(1,1,3)
!      
!      WRITE(*,*) "rho_g%pw%in_use", rho_g(1)%pw%in_use
!      WRITE(*,*) "rho_g(1)%pw%cc(1)", rho_g(1)%pw%cc(1)
!      WRITE(*,*) "rho_g(1)%pw%cc(2)", rho_g(1)%pw%cc(2)
!      WRITE(*,*) "rho_g(1)%pw%cc(262144)", rho_g(1)%pw%cc(262144)
!
!
!      WRITE(*,*) "v_rspace%pw%in_use", v_rspace(1)%pw%in_use
!      WRITE(*,*) "v_rspace(1)%pw%cr3d(1,1,1)", v_rspace(1)%pw%cr3d(1,1,1)
!      WRITE(*,*) "v_rspace(1)%pw%cr3d(1,1,2)", v_rspace(1)%pw%cr3d(1,1,2)
!      WRITE(*,*) "v_rspace(1)%pw%cr3d(1,1,3)", v_rspace(1)%pw%cr3d(1,1,3)


      WRITE(*,*) repeat("-",80)  
    END IF




    DO ispin = 1, nspins
       ! Add v_hartree + v_xc = v_rspace
       CALL pw_scale(v_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
       CALL pw_axpy(v_hartree_rspace%pw, v_rspace(ispin)%pw)
       ! integrate over potential <a|V|b>
    
       ! overwrites KS matrix
       NULLIFY (my_matrix_ks)
       CALL dbcsr_allocate_matrix_set(my_matrix_ks, nspins, 1)
       ALLOCATE (my_matrix_ks(ispin,1)%matrix)
       CALL dbcsr_copy(my_matrix_ks(ispin,1)%matrix, ec_env%matrix_ks(ispin,1)%matrix, & 
                       keep_sparsity = .FALSE., name="DUMMY KS-MATRIX")

       CALL integrate_v_rspace(v_rspace = v_rspace(ispin), &
                               pmat = ec_env%matrix_p(ispin, 1), &
                               hmat = my_matrix_ks(ispin, 1), &
                               qs_env = qs_env, &
                               calculate_forces = calculate_forces, &
                               basis_type = "HARRIS", &
                               task_list_external = ec_env%task_list)
    
       ! Deallocate dummy ks matrix                         
       CALL dbcsr_deallocate_matrix_set(my_matrix_ks)
    END DO 

    ! Core-Core repulsion forces (force%rho_core)
    CALL integrate_v_core_rspace(v_hartree_rspace, qs_env)

    ! return pw grids
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)

    DO ispin = 1, nspins
       CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
       DEALLOCATE(v_rspace)
       IF (ASSOCIATED(v_tau_rspace)) THEN
          CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
        DEALLOCATE(v_tau_rspace)
       END IF
    ENDDO
   
    CALL calculate_ecore_overlap(qs_env, para_env, calculate_forces = calculate_forces)

    CALL timestop(handle)

  END SUBROUTINE ec_bs_forces

! **************************************************************************************************
!> \brief  Non-self-consistent forces (Hartree and XC part)
!> \param qs_env ...
!>  \int (v_hartree[rho_out-rho_in] + (rho_out - rho_in)*\mu'_[rho_in]) \nabla_I rho_in dr

!> \date    03.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_nsc_force(qs_env, ec_env)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(energy_correction_type), POINTER              :: ec_env

    CHARACTER(len=*), PARAMETER :: routineN = 'ec_nsc', routineP = moduleN//':'//routineN
    
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
    INTEGER                                            :: i, iatom, ikind, natom
    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
    ! fbelle
    TYPE(xc_derivative_set_type), POINTER              :: deriv_set
    TYPE(xc_rho_set_type), POINTER                     :: rho_set
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc, v_rspace_nsc
    TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: my_matrix_ks
    TYPE(section_vals_type), POINTER                   :: input, xc_section

    logger => cp_get_default_logger()
    
    CALL timeset(routineN//"_forces", handle)
    
    CALL get_qs_env(qs_env, &
                    kg_env = kg_env, &
                    dft_control = dft_control, &
                    rho = rho, &
                    force = force, &
                    input = input, &
                    para_env = para_env)

!    IF(para_env%mepos==0) THEN
!       WRITE(*,*) "ec_nsc_force - BEGINNING"
!    END IF
    
    calculate_forces = .TRUE.
    
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
                           xc_section = xc_section)        ! OK 

    CALL xc_calc_2nd_deriv(v_xc = v_xc, & 
                           deriv_set = deriv_set, &
                           rho_set = rho_set, &
                           rho1_set = rho_set, &
                           pw_pool = auxbas_pw_pool, &
                           xc_section = xc_section, & 
                           gapw = .FALSE.)

    DO ispin = 1,nspins
      !CALL pw_scale(v_xc(ispin)%pw,v_xc(ispin)%pw%pw_grid%dvol)
      ! \nu'(rho_in)*(rho_out - rho_in)
      v_rspace_nsc(ispin)%pw%cr3d = (rho_r(ispin)%pw%cr3d &
                                    - old_rho_r(ispin)%pw%cr3d) &
                                    * v_xc(ispin)%pw%cr3d
    END DO

    CALL xc_dset_release(deriv_set)
    CALL xc_rho_set_release(rho_set = rho_set, pw_pool = auxbas_pw_pool)

    !!---------------------------------------------------------------------------------------
        ! *** Calculation of forces due to non-self consistency of energy correction *** ! 
        ! Second, add v_hartree[rho_out-rho_in] to potential we will integrate 
    !!---------------------------------------------------------------------------------------
   
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
    
    
    !!---------------------------------------------------------------------------------------
        ! *** Calculation of forces due to non-self consistency of energy correction *** ! 
        ! \int v_nsc(r) * \grad_I(n_in) 
        ! old density matrix, thus 'old' task_list
    !!---------------------------------------------------------------------------------------
    !DO ispin = 1 , nspins
    !   ! modifies KS matrix!
    !   CALL integrate_v_rspace(v_rspace = v_rspace_nsc(ispin), &
    !                           pmat = old_rho_ao(ispin, 1), &
    !                           hmat = ec_env%matrix_p(ispin, 1), &
    !                           qs_env = qs_env, &
    !                           calculate_forces = calculate_forces, &
    !                           basis_type = "HARRIS", &
    !                           task_list_external = ec_env%task_list)
    !END DO

    NULLIFY (my_matrix_ks)
    
    DO ispin = 1 , nspins
       ! modifies KS matrix!
       CALL dbcsr_allocate_matrix_set(my_matrix_ks, nspins, 1)
       ALLOCATE (my_matrix_ks(ispin,1)%matrix)

       CALL dbcsr_copy(my_matrix_ks(ispin,1)%matrix, ec_env%matrix_ks(ispin,1)%matrix, & 
                       keep_sparsity = .FALSE., name="DUMMY KS-MATRIX")
   
       CALL integrate_v_rspace(v_rspace = v_rspace_nsc(ispin), &
                               pmat = old_rho_ao(ispin, 1), &
                               hmat = my_matrix_ks(ispin, 1), &
                               qs_env = qs_env, &
                               calculate_forces = calculate_forces)
    
       ! Deallocate dummy ks matrix                         
       CALL dbcsr_deallocate_matrix_set(my_matrix_ks)
    END DO

    !---------------------------------------------------------------------------------------
    ! Need to copy NSC force contrib from qs_env%...%force to stored forces in ec_env%force
    ! acting on rho_elec 
    !---------------------------------------------------------------------------------------

    ! Substract NSC Hartree force contribution from saved forces in ec_env%force
    DO iatom = 1, natom
       ikind = kind_of(iatom)
       i = atom_of_kind(iatom)
       
       ec_env%force(ikind)%rho_elec(1:3, i) = ec_env%force(ikind)%rho_elec(1:3, i) & 
                                                + force(ikind)%rho_elec(1:3, i)
    END DO

    DEALLOCATE (atom_of_kind, kind_of)

    ! gather NSF force accross mpi ranks to print
    DO ikind = 1,SIZE(force)
      CALL mp_sum(force(ikind)%rho_elec, para_env%group)
    END DO

    IF(para_env%mepos==0) THEN
      WRITE(*,*) repeat("-",80)  
      WRITE(*,*) "ec_nsc_forces start - qs_env%...%force - rho_elec"
      CALL ec_write_forces(qs_env, force, ftype = 3)
      WRITE(*,*) "ec_nsc_forces end - qs_env%...%force - rho_elec"
!      WRITE(*,*) "ec_nsc_forces - ec_env%force - rho_elec"
!      CALL ec_write_forces(qs_env, ec_env%force, ftype = 3)
      WRITE(*,*) repeat("-",80)  
    END IF

    ! Delete NSC force contribution from qs_env%...%force 
    CALL zero_qs_force(force)
    
!    IF(para_env%mepos==0) THEN
!       WRITE(*,*) "ec_nsc_force - qs_env%...%force - after copy of rho_elec"
!       CALL ec_write_forces(qs_env, force)
!       WRITE(*,*) "ec_force - ec_env%force - after copy of rho_elec"
!       CALL ec_write_forces(qs_env, ec_env%force)
!    END IF


    DO ispin = 1,nspins
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_nsc(ispin)%pw)
    END DO

    DEALLOCATE(v_xc, v_rspace_nsc)

    ! release pw grids and matrices
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)

!    IF(para_env%mepos==0) WRITE(*,*) "ec_nsc_force - END"
!    CALL ec_write_forces(qs_env, force, calculate_forces=.TRUE.)

    CALL timestop(handle)

  END SUBROUTINE ec_nsc_force

! **************************************************************************************************
!> \brief  Non-self-consistent forces (Hartree part)
!> \param qs_env ...
!>  \int (v_hartree[rho_out-rho_in] \nabla_I rho_in dr

!> \date    06.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_nsc_hartree_force(qs_env, ec_env)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(energy_correction_type), POINTER              :: ec_env

    CHARACTER(len=*), PARAMETER :: routineN = 'ec_nsc_hartree', routineP = moduleN//':'//routineN
    
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
    INTEGER                                            :: i, iatom, ikind, natom
    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
    TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: my_matrix_ks
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace_nsc


    logger => cp_get_default_logger()
    
    CALL timeset(routineN//"_forces", handle)
    
    CALL get_qs_env(qs_env, &
                    dft_control = dft_control, &
                    rho = rho, &
                    force = force, &
                    para_env = para_env)

    IF(para_env%mepos==0) THEN
       WRITE(*,*) "ec_nsc_hartree_force - BEGINNING"
    END IF
    
    calculate_forces = .TRUE.
    
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


    
    ALLOCATE(v_rspace_nsc(nspins))
    DO ispin = 1, nspins
       NULLIFY(v_rspace_nsc(nspins)%pw)
       CALL pw_pool_create_pw(pool = auxbas_pw_pool, & 
                              pw = v_rspace_nsc(ispin)%pw, &
                              use_data = REALDATA3D, &
                              in_space = REALSPACE)
       CALL pw_zero(v_rspace_nsc(nspins)%pw)
    END DO

    !!---------------------------------------------------------------------------------------
        ! *** Calculation of forces due to non-self consistency of energy correction *** ! 
        ! Second, add v_hartree[rho_out-rho_in] to potential we will integrate 
    !!---------------------------------------------------------------------------------------
   
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


    ! Add v_NSC_hartree
    DO ispin = 1,nspins
!      CALL pw_scale(v_rspace_nsc(ispin)%pw, v_rspace_nsc(ispin)%pw%pw_grid%dvol)
      CALL pw_axpy(v_hartree_rspace%pw,v_rspace_nsc(ispin)%pw)
    END DO


    DO ispin = 1 , nspins
       
       ! modifies KS matrix!
       NULLIFY (my_matrix_ks)
       CALL dbcsr_allocate_matrix_set(my_matrix_ks, nspins, 1)
       ALLOCATE (my_matrix_ks(ispin,1)%matrix)
       CALL dbcsr_copy(my_matrix_ks(ispin,1)%matrix, ec_env%matrix_ks(ispin,1)%matrix, & 
                       keep_sparsity = .FALSE., name="DUMMY KS-MATRIX")
   
       CALL integrate_v_rspace(v_rspace = v_rspace_nsc(ispin), &
                               pmat = old_rho_ao(ispin, 1), &
                               hmat = my_matrix_ks(ispin, 1), &
                               qs_env = qs_env, &
                               calculate_forces = calculate_forces)
    
       ! Deallocate dummy ks matrix                         
       CALL dbcsr_deallocate_matrix_set(my_matrix_ks)
    END DO

    !---------------------------------------------------------------------------------------
    ! Need to copy NSC force contrib from qs_env%...%force to stored forces in ec_env%force
    ! acting on rho_elec 
    !---------------------------------------------------------------------------------------

    ! Substract NSC Hartree force contribution from saved forces in ec_env%force
    DO iatom = 1, natom
       ikind = kind_of(iatom)
       i = atom_of_kind(iatom)
       
!       ec_env%force(ikind)%rho_elec(1:3, i) = ec_env%force(ikind)%rho_elec(1:3, i) & 
!                                                + force(ikind)%rho_elec(1:3, i)
    END DO

    DEALLOCATE (atom_of_kind, kind_of)

    ! gather NSF force accross mpi ranks to print
    DO ikind = 1,SIZE(force)
      CALL mp_sum(force(ikind)%rho_elec, para_env%group)
    END DO

    IF(para_env%mepos==0) THEN
      WRITE(*,*) repeat("-",80)  
      WRITE(*,*) "ec_nsc_hartree_forces start - qs_env%...%force - rho_elec"
      CALL ec_write_forces(qs_env, force, ftype = 3)
      WRITE(*,*) "ec_nsc_hartree_forces end - qs_env%...%force - rho_elec"
!      WRITE(*,*) "ec_nsc_forces - ec_env%force - rho_elec"
!      CALL ec_write_forces(qs_env, ec_env%force, ftype = 3)
      WRITE(*,*) repeat("-",80)  
    END IF

    ! Delete NSC force contribution from qs_env%...%force 
    CALL zero_qs_force(force)

    DO ispin = 1,nspins
       CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_nsc(ispin)%pw)
    END DO

    DEALLOCATE(v_rspace_nsc)

    ! release pw grids and matrices
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)

    CALL timestop(handle)

  END SUBROUTINE ec_nsc_hartree_force

! **************************************************************************************************
!> \brief  Non-self-consistent forces (XC part)
!> \param qs_env ...
!>  \int (rho_out - rho_in)*\mu'_[rho_in]) \nabla_I rho_in dr

!> \date    03.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_nsc_xc_force(qs_env, ec_env)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    TYPE(energy_correction_type), POINTER              :: ec_env

    CHARACTER(len=*), PARAMETER :: routineN = 'ec_nsc_xc', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle
    TYPE(cp_logger_type), POINTER                      :: logger
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
    ! fbelle
    TYPE(pw_env_type), POINTER                         :: pw_env
    TYPE(qs_rho_type), POINTER                         :: rho!, old_rho
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, old_rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, old_rho_g
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao, old_rho_ao
    TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
    TYPE(dft_control_type), POINTER                    :: dft_control
    INTEGER                                            :: ispin, nspins
    LOGICAL                                            :: calculate_forces 
    TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
    INTEGER                                            :: i, iatom, ikind, natom
    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
    ! fbelle
    TYPE(xc_derivative_set_type), POINTER              :: deriv_set
    TYPE(xc_rho_set_type), POINTER                     :: rho_set!, rho1_set
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc, v_rspace_nsc
    TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: my_matrix_ks
!    TYPE(section_vals_type), POINTER                   :: input
    ! fbelle
!    TYPE(pw_p_type), DIMENSION(:), POINTER             :: tau
!    TYPE(qs_ks_env_type), POINTER                      :: ks_env
    TYPE(section_vals_type), POINTER                   :: input, xc_section
!    TYPE(section_vals_type), POINTER                   :: xc_fun_section
!    INTEGER, DIMENSION(2, 3)                           :: bo
!    TYPE(xc_rho_cflags_type)                           :: needs
!    LOGICAL                                            :: lsd

    logger => cp_get_default_logger()
    
    CALL timeset(routineN//"_forces", handle)
    
    CALL get_qs_env(qs_env, &
                    dft_control = dft_control, &
                    rho = rho, &
                    force = force, &
                    input = input, &
                    para_env = para_env)

    calculate_forces = .TRUE.
    
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
    CALL get_qs_env(qs_env, pw_env = pw_env)
    CPASSERT(ASSOCIATED(pw_env))

    NULLIFY (auxbas_pw_pool, pw_pools)
    ! gets the tmp grids
    CALL pw_env_get(pw_env = pw_env, &
                    auxbas_pw_pool = auxbas_pw_pool, &
                    pw_pools = pw_pools)

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
                           xc_section = xc_section) ! OK 

    CALL xc_calc_2nd_deriv(v_xc = v_xc, & 
                           deriv_set = deriv_set, &
                           rho_set = rho_set, &
                           rho1_set = rho_set, &
                           pw_pool = auxbas_pw_pool, &
                           xc_section = xc_section, & 
                           gapw = .FALSE.)

    DO ispin = 1,nspins
      !CALL pw_scale(v_xc(ispin)%pw,v_xc(ispin)%pw%pw_grid%dvol)
      ! \nu'(rho_in)*(rho_out - rho_in)
      v_rspace_nsc(ispin)%pw%cr3d = (rho_r(ispin)%pw%cr3d &
                                  - old_rho_r(ispin)%pw%cr3d) &
                                  * v_xc(ispin)%pw%cr3d
    END DO

    CALL xc_dset_release(deriv_set)
    CALL xc_rho_set_release(rho_set = rho_set, pw_pool = auxbas_pw_pool)

    DO ispin = 1,nspins
      CALL pw_scale(v_rspace_nsc(ispin)%pw, v_rspace_nsc(ispin)%pw%pw_grid%dvol)
    END DO
    
    DO ispin = 1 , nspins
       
       ! modifies KS matrix!
       NULLIFY (my_matrix_ks)
       CALL dbcsr_allocate_matrix_set(my_matrix_ks, nspins, 1)
       ALLOCATE (my_matrix_ks(ispin,1)%matrix)

       CALL dbcsr_copy(my_matrix_ks(ispin,1)%matrix, ec_env%matrix_ks(ispin,1)%matrix, & 
                       keep_sparsity = .FALSE., name="DUMMY KS-MATRIX")
   
       CALL integrate_v_rspace(v_rspace = v_rspace_nsc(ispin), &
                               pmat = old_rho_ao(ispin, 1), &
                               hmat = my_matrix_ks(ispin, 1), &
                               qs_env = qs_env, &
                               calculate_forces = calculate_forces)
    
       ! Deallocate dummy ks matrix                         
       CALL dbcsr_deallocate_matrix_set(my_matrix_ks)
    END DO

    !---------------------------------------------------------------------------------------
    ! Need to copy NSC force contrib from qs_env%...%force to stored forces in ec_env%force
    ! acting on rho_elec 
    !---------------------------------------------------------------------------------------

    ! Substract NSC Hartree force contribution from saved forces in ec_env%force
    DO iatom = 1, natom
       ikind = kind_of(iatom)
       i = atom_of_kind(iatom)
       
!       ec_env%force(ikind)%rho_elec(1:3, i) = ec_env%force(ikind)%rho_elec(1:3, i) & 
!                                                + force(ikind)%rho_elec(1:3, i)
    END DO

    DEALLOCATE (atom_of_kind, kind_of)

    ! gather NSF force accross mpi ranks to print
    DO ikind = 1,SIZE(force)
      CALL mp_sum(force(ikind)%rho_elec, para_env%group)
    END DO

    IF(para_env%mepos==0) THEN
      WRITE(*,*) repeat("-",80)  
      WRITE(*,*) "ec_nsc_xc_forces start - qs_env%...%force - rho_elec"
      CALL ec_write_forces(qs_env, force, ftype = 3)
      WRITE(*,*) "ec_nsc_xc_forces end - qs_env%...%force - rho_elec"
!      WRITE(*,*) "ec_nsc_forces - ec_env%force - rho_elec"
!      CALL ec_write_forces(qs_env, ec_env%force, ftype = 3)
      WRITE(*,*) repeat("-",80)  
    END IF

    ! Delete NSC force contribution from qs_env%...%force 
    CALL zero_qs_force(force)

    DO ispin = 1,nspins
       CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
       CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_nsc(ispin)%pw)
    END DO

    DEALLOCATE(v_xc, v_rspace_nsc)

    CALL timestop(handle)

  END SUBROUTINE ec_nsc_xc_force

END MODULE ec_force
