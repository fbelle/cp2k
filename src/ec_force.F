!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

!fbelle
! MODIFY for harris energy (and force correction) just like qs_energy.F and qs_force.F routines
! NEW : contains just force routines necessary

! **************************************************************************************************
!> \brief Quickstep force driver routine
!> \author MK (12.06.2002)

!> Force contributions : ec_core_hamiltonian  kinetic, overlap, gth_ppl, gth_ppnl
!>                       gth_nlcc ... not implemented


! **************************************************************************************************
MODULE ec_force
   USE admm_methods,                    ONLY: calc_aux_mo_derivs_none,&
                                              calc_mixed_overlap_force
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_deallocate_matrix_set,&
                                              dbcsr_finalize,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_print,&
                                              dbcsr_release,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric
   USE dft_plus_u,                      ONLY: plus_u
   USE efield_utils,                    ONLY: calculate_ecore_efield
   USE input_constants,                 ONLY: do_admm_purify_none,&
                                              outer_scf_hirshfeld_constraint
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_environment_type
   USE message_passing,                 ONLY: mp_sum
   USE mulliken,                        ONLY: mulliken_restraint
   USE particle_types,                  ONLY: particle_type
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
   USE qs_dftb_dispersion,              ONLY: calculate_dftb_dispersion
   USE qs_dftb_matrices,                ONLY: build_dftb_matrices
   USE qs_energy,                       ONLY: qs_energies
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_external_potential,           ONLY: external_c_potential,&
                                              external_e_potential
   USE qs_force_types,                  ONLY: allocate_qs_force,&
                                              qs_copy_force,&
                                              qs_force_type,&
                                              replicate_qs_force,&
                                              zero_qs_force
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env,&
                                              calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type,&
                                              qs_rho_clear
   USE qs_scf_post_scf,                 ONLY: qs_scf_compute_properties
   USE qs_subsys_types,                 ONLY: qs_subsys_set,&
                                              qs_subsys_type
   USE ri_environment_methods,          ONLY: build_ri_matrices
   USE rt_propagation_forces,           ONLY: calc_c_mat_force,&
                                              rt_admm_force
   USE se_core_core,                    ONLY: se_core_core_interaction
   USE se_core_matrix,                  ONLY: build_se_core_matrix
   USE virial_types,                    ONLY: virial_type
   !fbelle
   USE kg_environment_types,            ONLY: energy_correction_type,&
                                              kg_environment_type
   USE kg_energy_corrections,           ONLY: ec_build_core_hamiltonian, &
                                              ec_build_ks_matrix,&
                                              ec_calculate_w_matrix
   ! fbelle pw stuff
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integral_ab,&
                                              pw_scale,&
                                              pw_transfer,&
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
                                              pw_p_type
   USE ec_methods,                      ONLY: ec_write_forces,&                                           
                                              ec_write_forces2                                           
   USE qs_ks_utils,                     ONLY: print_densities
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_rho_methods,                  ONLY: duplicate_rho_type                                           
                                           
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_force'

! *** Public subroutines ***

   PUBLIC :: ec_forces, ec_nsc_force

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

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_forces', routineP = moduleN//':'//routineN
    
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
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r
!    TYPE(pw_p_type), DIMENSION(:), POINTER             :: old_rho_r, old_rho_g
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: old_rho_ao
    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: natom_of_kind
    TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
    !fbelle
!    TYPE(qs_force_type), DIMENSION(:), POINTER         :: ec_force
!    INTEGER                                            :: natom
    TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set


!!---------------------------------------------------------------------------------------
!!  Initial general setup
!!---------------------------------------------------------------------------------------

    logger => cp_get_default_logger()
    
    CALL timeset(routineN//"_forces", handle)
    
    CALL get_qs_env(qs_env,para_env=para_env)
    IF(para_env%mepos==0) THEN
       WRITE(*,*) "ec_forces - BEGINNING"
    END IF

    ! obviously 
    calculate_forces = .TRUE.
    
    ! get qs and kg environment
    NULLIFY (para_env, rho_ao, old_rho_ao, dft_control)
    CALL get_qs_env(qs_env=qs_env, kg_env=kg_env, dft_control=dft_control, para_env=para_env)

    !shortcut to ec_type
    ! ec_build_neighborlist has been built during kg_energy_correction call
    ec_env => kg_env%ec_env
    nspins = dft_control%nspins


    ! Introduce new ec_env%force type to store forces and prevent overwriting
    !    especially to rho_elec
    ! get atom mapping
    NULLIFY (atomic_kind_set)
    CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, particle_set=particle_set)

    IF (.NOT. ASSOCIATED(ec_env%force)) THEN
       nkind = SIZE(atomic_kind_set)
       ALLOCATE (natom_of_kind(nkind))
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                natom_of_kind=natom_of_kind)
       CALL allocate_qs_force(ec_env%force, natom_of_kind)
       DEALLOCATE (natom_of_kind)
    END IF
    CALL zero_qs_force(ec_env%force)

    ! store old rho (before energy correction)
    ! after kg_energy_correction rho_structure is not updated
    ! only new density matrix : ec_env%matrix_p 
    CALL get_qs_env(qs_env,rho=rho)
    CALL qs_rho_get(rho, rho_ao_kp=rho_ao, rho_r=rho_r, rho_g=rho_g)
    CALL qs_rho_create(ec_env%old_rho)

    ! copy rho structure, alert, sets density matrix to zero
    ! rho_in  = rho
    ! rho_out = ec_env%old_rho
    CALL duplicate_rho_type(rho, ec_env%old_rho, qs_env)
                   
    ! duplicate_rho_types doesnt copy density matrix ...
    CALL qs_rho_get(ec_env%old_rho, rho_ao_kp=old_rho_ao)
    DO ispin = 1, nspins
       CALL dbcsr_copy(old_rho_ao(ispin,1)%matrix, rho_ao(ispin,1)%matrix, keep_sparsity=.FALSE.)
    END DO

    ! update rho strucutre with new density matrix of energy correction
    ! 'rho_ao_in'  = ec_env%matrix_p 
    ! 'rho_ao_out' = rho%rho_ao(:,:) 
     
    CALL get_qs_env(qs_env, rho=rho)
    CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
    DO ispin = 1, nspins
       CALL dbcsr_copy(rho_ao(ispin,1)%matrix, ec_env%matrix_p(ispin,1)%matrix ,keep_sparsity=.FALSE.)
    END DO
    
    ! updates rho_r and rho_g wrt to new rho%rho_ao
    CALL qs_rho_update_rho(rho, qs_env=qs_env)

    ! calculate matrix_w = matrix_p*matrix_ks*matrix_p
    !IF(para_env%mepos==0) THEN 
    !   WRITE(*,*) "ec_force - ec_calculate_w_matrix"
    !END IF
    CALL ec_calculate_w_matrix(qs_env, ec_env)

    ! *** d/dR[Sum of eigenvalues] = d/dR[trace(H_in * rho_out)] *** !
    !core hamiltonian force contribution: 
    ! kinetic, 
    ! overlap, 
    ! gth_ppl, 
    ! gth_ppnl, 
    CALL ec_build_core_hamiltonian(qs_env, ec_env, calculate_forces)

    ! forces from core overlap
    ! core_overlap
    ! already in ec_build_ks_matrix
    ! CALL calculate_ecore_overlap(qs_env, para_env, calculate_forces=.TRUE.)

    ! next are all grid based forces : 
    ! rho_elec
    ! rho_core
    ! contains CALL calculate_ecore_overlap 

    CALL ec_build_ks_matrix(qs_env, ec_env, calculate_forces)
    
    ! At this point: F_band + F_core have been calculated
    ! Missing: NSC terms F_\delta_H and F_\delta_XC

    ! Store these forces contributions in ec_env%force
    CALL qs_copy_force(ec_env%force, force)

    ! empty qs force type
    CALL zero_qs_force(force)
  
    CALL get_qs_env(qs_env, force=force) 
    WRITE(*,*) "ec_force - qs_env%...%force"
    CALL ec_write_forces2(qs_env, force)

    WRITE(*,*) "ec_force - ec_env%force - after zero_qs-force"
    !CALL get_qs_env(qs_env, force=force)
    CALL ec_write_forces2(qs_env, ec_env%force)

    ! forces from double counting terms
    CALL ec_nsc_force(qs_env)

    WRITE(*,*) "ec_force - qs_env%...%force - after ec_nsc_force"
    CALL ec_write_forces2(qs_env, force)

    WRITE(*,*) "ec_force - ec_env%force - after ec_nsc_force"
    !CALL get_qs_env(qs_env, force=force)
    CALL ec_write_forces2(qs_env, ec_env%force)


!    CALL get_qs_env(qs_env, force=force)
!    WRITE(*,*) "WRITE 1 FUCKING FORCE"
!    WRITE(*,*) "force(1)%overlap(1,1)", force(1)%overlap(1,1)
!    CALL ec_write_forces2(qs_env, force)


    !  *** replicate forces ***
    ! sums up of current forces in total
    CALL replicate_qs_force(force, para_env)


!    WRITE(*,*) "Forces END of ec_forces"  
!     CALL ec_write_forces(qs_env, ec_env%force, calculate_forces=.TRUE.)





    CALL timestop(handle)

  END SUBROUTINE ec_forces


! **************************************************************************************************
!> \brief  Non self consistent forces (Hartree and XC part)
!> \param qs_env ...

!> \date    03.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_nsc_force(qs_env)

    TYPE(qs_environment_type), POINTER                 :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'ec_nsc', routineP = moduleN//':'//routineN
    
    TYPE(energy_correction_type), POINTER              :: ec_env
    TYPE(kg_environment_type), POINTER                 :: kg_env
    INTEGER                                            :: handle
    TYPE(cp_logger_type), POINTER                      :: logger
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
    ! fbelle
    TYPE(pw_env_type), POINTER                         :: pw_env
    TYPE(qs_rho_type), POINTER                         :: rho!, old_rho
    TYPE(pw_p_type)                                    :: rho_delta_gspace, v_hartree_gspace, &
                                                          v_hartree_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, rho_g, old_rho_r, old_rho_g
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao, old_rho_ao
    TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
    TYPE(pw_poisson_type), POINTER                     :: poisson_env
    TYPE(dft_control_type), POINTER                    :: dft_control
    INTEGER                                            :: ispin, nspins
    LOGICAL                                            :: calculate_forces 
    TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_delta

    logger => cp_get_default_logger()
    
    CALL timeset(routineN//"_forces", handle)
    
    CALL get_qs_env(qs_env, kg_env=kg_env, dft_control=dft_control, rho=rho, force=force, para_env=para_env)

    ec_env => kg_env%ec_env

    calculate_forces = .TRUE.
    
    
    WRITE(*,*) "ec_nsc_force - BEGINNING"
!    CALL zero_qs_force(force)
!    CALL ec_write_forces(qs_env, ec_env%force, calculate_forces=.TRUE.)

    nspins = dft_control%nspins
     
    ! *** Calculation of forces due to non-self consistency of energy correction *** ! 
    ! First, v_hartree[rho_out-rho_in]

    NULLIFY (pw_env)
    CALL get_qs_env(qs_env, pw_env=pw_env)
    CPASSERT(ASSOCIATED(pw_env))

    NULLIFY (auxbas_pw_pool, poisson_env, pw_pools)
    ! gets the tmp grids
    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                    pw_pools=pw_pools, poisson_env=poisson_env)

    ! Calculate the NSC Hartree potential
    CALL pw_pool_create_pw(auxbas_pw_pool, &
                           v_hartree_gspace%pw, &
                           use_data=COMPLEXDATA1D, &
                           in_space=RECIPROCALSPACE)
    CALL pw_pool_create_pw(auxbas_pw_pool, &
                           rho_delta_gspace%pw, &
                           use_data=COMPLEXDATA1D, &
                           in_space=RECIPROCALSPACE)

    CALL pw_pool_create_pw(auxbas_pw_pool, &
                           v_hartree_rspace%pw, &
                           use_data=REALDATA3D, &
                           in_space=REALSPACE)

    ! get rho_out - after Harris energy correction
    CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, rho_ao_kp=rho_ao)

    ! get rho_in - before Harris energy correction
    CALL qs_rho_get(ec_env%old_rho, rho_r=old_rho_r, rho_g=old_rho_g, rho_ao_kp=old_rho_ao )

!    WRITE(*,*) "ec_nsc_force - rho_ao "
!    IF(ASSOCIATED(rho_ao)) THEN
!        CALL dbcsr_print(rho_ao(1, 1)%matrix)
!    END IF
!    WRITE(*,*) "ec_nsc_force - old_rho_ao "
!    IF(ASSOCIATED(rho_ao)) THEN
!        CALL dbcsr_print(old_rho_ao(1)%matrix)
!    END IF
   

    ! create delta density matrix
    IF (ASSOCIATED(ec_env%matrix_delta)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_delta)
    CALL dbcsr_allocate_matrix_set(ec_env%matrix_delta, nspins, 1)
    matrix_delta => ec_env%matrix_delta
    DO ispin = 1, nspins
       ALLOCATE (matrix_delta(ispin, 1)%matrix)
       CALL dbcsr_create(matrix_delta(ispin, 1)%matrix, name="DELTA DENSITY MATRIX", &
                         template=ec_env%matrix_p(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
       CALL cp_dbcsr_alloc_block_from_nbl(matrix_delta(ispin, 1)%matrix, ec_env%sab_orb)
       ! matrix_delta = rho_ao - old_rho_ao
       CALL dbcsr_copy(matrix_delta(ispin, 1)%matrix, rho_ao(ispin, 1)%matrix)
       CALL dbcsr_add(matrix_delta(ispin, 1)%matrix, old_rho_ao(ispin, 1)%matrix, alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
    END DO

!    WRITE(*,*) "ec_nsc_force - matrix_delta "
!    IF(ASSOCIATED(ec_env%matrix_delta)) THEN
!        CALL dbcsr_print(ec_env%matrix_delta(1, 1)%matrix)
!    END IF
!
!    WRITE(*,*) "ec_nsc_force - matrix_delta"
!    IF(ASSOCIATED(ec_env%matrix_delta)) THEN
!        CALL dbcsr_print(ec_env%matrix_delta(1, 1)%matrix)
!    END IF

    CALL print_densities(qs_env, rho)
   
!    ! + v_hartree[rho_out-rho_in]
    CALL pw_zero(rho_delta_gspace%pw)
    DO ispin = 1,nspins
      CALL pw_axpy(rho_g(ispin)%pw,rho_delta_gspace%pw)
      CALL pw_axpy(old_rho_g(ispin)%pw,rho_delta_gspace%pw,alpha=-1._dp)
    END DO

    ! get Hartree potential from rho_delta_gspace
    CALL pw_poisson_solve(poisson_env, rho_delta_gspace%pw, &
                           vhartree=v_hartree_gspace%pw)

    CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)

    DO ispin = 1, nspins
       CALL integrate_v_rspace(v_rspace=v_hartree_rspace, &
                               pmat=matrix_delta(ispin,1),&
                               hmat=ec_env%matrix_ks(ispin,1), &
                               qs_env=qs_env, &
                               calculate_forces=calculate_forces)
    END DO
    
    CALL get_qs_env(qs_env, force=force)











    ! release pw grids and matrices
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_delta_gspace%pw)

!    DO ispin = 1, nspins
!       CALL dbcsr_release(matrix_delta(ispin,1)%matrix)
!    END DO

    WRITE(*,*) "ec_nsc_force - END"
!    CALL ec_write_forces(qs_env, force, calculate_forces=.TRUE.)

    CALL timestop(handle)

  END SUBROUTINE ec_nsc_force

END MODULE ec_force
