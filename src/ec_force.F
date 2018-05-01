!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Harris functional force driver routine

!> Force contributions : ec_core_hamiltonian  kinetic, overlap, gth_ppl, gth_ppnl
!>                       gth_nlcc ... not implemented
!> F_tot = F_kin + F_gth_ppl + F_gth_ppnl + F_BS_H[n_out] + F_BS_XC[n_out] + F_core
!>         + F_NSC
!> w/ F_NSC = F_BS_H + F_BS_XC + F_delta_H + F_delta_XC
!> \date 02.2018
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_force
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_fm_types,                     ONLY: cp_fm_type
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
                                              dbcsr_finalize,&
                                              dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_print,&
                                              dbcsr_release,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric,&
                                              dbcsr_type_no_symmetry
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
   USE qs_energy,                       ONLY: qs_energies
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
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
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
                                              ec_calculate_w_matrix
   USE qs_integrate_potential,          ONLY: integrate_v_rspace
   USE qs_rho_methods,                  ONLY: duplicate_rho_type                                         !fbelle
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type
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
    CALL get_qs_env(qs_env=qs_env,& 
                    kg_env=kg_env,&
                    dft_control=dft_control,&
                    force = force,&
                    para_env=para_env)

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
    CALL get_qs_env(qs_env, rho = rho)
    CALL qs_rho_get(rho, &
                    rho_ao_kp = rho_ao, &
                    rho_r = rho_r, &
                    rho_g = rho_g) 
    CALL qs_rho_create(ec_env%old_rho)

    ! copy rho structure, alert, sets density matrix to zero
    CALL duplicate_rho_type(rho_input = rho, &
                            rho_output = ec_env%old_rho, &
                            qs_env = qs_env)

    ! duplicate_rho_types doesnt copy density matrix ...
    CALL qs_rho_get(ec_env%old_rho, rho_ao_kp=old_rho_ao)
    DO ispin = 1, nspins
       CALL dbcsr_copy(old_rho_ao(ispin,1)%matrix, rho_ao(ispin,1)%matrix, keep_sparsity=.FALSE.)
    END DO

    ! update rho strucutre with new density matrix of energy correction
    ! 'rho_ao_in'  = ec_env%matrix_p 
    ! 'rho_ao_out' = rho%rho_ao(:,:) 
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

    ! next are all grid based forces : 
    ! rho_elec - F_BS_XC[n_out] + F_BS_H[n_out]
    ! rho_core
    ! contains CALL calculate_ecore_overlap 
    CALL ec_build_ks_matrix(qs_env, ec_env, calculate_forces)

    ! At this point: F_band + F_core have been calculated
    ! Missing: NSC terms F_NSC + F_rho_core i.e. grid based forces

    ! store these forces contributions in ec_env%force
    !                  to            from
    CALL qs_copy_force(ec_env%force, force)

    ! empty qs force type
    CALL zero_qs_force(force)
  
    ! forces from band structure and double counting terms
    ! grid based forces rho_elec
    CALL ec_nsc_force(qs_env)

    ! copy back forces
    CALL qs_copy_force(force, ec_env%force)

    IF (para_env%mepos==0) THEN
       WRITE(*,*) "Forces END of ec_forces"  
    END IF

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
    TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: density_matrix
    !fbelle
    TYPE(xc_derivative_set_type), POINTER              :: deriv_set
    TYPE(xc_rho_set_type), POINTER                     :: rho_set
    TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_xc, v_rspace_nsc

    logger => cp_get_default_logger()
    
    CALL timeset(routineN//"_forces", handle)
    
    
    CALL get_qs_env(qs_env, &
                    kg_env=kg_env, &
                    dft_control=dft_control, &
                    rho=rho, &
                    force=force, &
                    para_env=para_env)

    IF(para_env%mepos==0) THEN
       WRITE(*,*) "ec_nsc_force - BEGINNING"
    END IF
    
    ec_env => kg_env%ec_env

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

    NULLIFY (rho_set, deriv_set)
    ! main ingredient is xc_rho_set_and_dset_create
    CALL xc_prep_2nd_deriv(deriv_set = deriv_set, &        ! containing potentials
                           rho_set = rho_set, &            ! density at which deriv are calced
                           rho_r = old_rho_r, &      ! where derivative evaluated 
                           pw_pool = auxbas_pw_pool, &     ! OK
                           xc_section = ec_env%xc_section) ! OK 

    CALL xc_calc_2nd_deriv(v_xc = v_xc, & 
                           deriv_set = deriv_set, &
                           rho_set = rho_set, &
                           rho1_set = rho_set, &
                           pw_pool = auxbas_pw_pool, &
                           xc_section = ec_env%xc_section, & 
                           gapw = .FALSE.)

    DO ispin = 1,nspins
      CALL pw_scale(v_xc(ispin)%pw,v_xc(ispin)%pw%pw_grid%dvol)
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
      CALL pw_axpy(old_rho_g(ispin)%pw,rho_tot_gspace%pw,alpha=-1._dp)
    END DO

    ! get Hartree potential from rho_tot_gspace
    CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, &
                           vhartree=v_hartree_gspace%pw)
    CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
    CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

    ! Add v_NSC_xc + v_NSC_hartree
    DO ispin = 1,nspins
      CALL pw_axpy(v_hartree_rspace%pw,v_rspace_nsc(ispin)%pw)
      CALL pw_scale(v_rspace_nsc(ispin)%pw, v_rspace_nsc(ispin)%pw%pw_grid%dvol)
    END DO
    
!    IF(para_env%mepos==0) THEN
!    WRITE(*,*) "ec_nsc_force - qs_env%...%force - before  rho_elec"
!    CALL ec_write_forces(qs_env, force)
!    WRITE(*,*) "ec_force - ec_env%force before rho_elec "
!    CALL ec_write_forces(qs_env, ec_env%force)
!    END IF
    
    ! empty qs force 
    CALL zero_qs_force(force)

    ! Careful, old_rho_ao is sparse, however hmat is not (task_list conflict)
    ! Allocate dummy matrix (density_matrix) to transform old_rho_ao into non-sparse matrix 
    NULLIFY (density_matrix)
    CALL dbcsr_allocate_matrix_set(density_matrix, nspins, 1)
    DO ispin = 1, nspins
       ALLOCATE (density_matrix(ispin,1)%matrix)
       CALL dbcsr_create(density_matrix(ispin, 1)%matrix, name="NON-SPARSE OLD P-MATRIX", &
                         template=rho_ao(1, 1)%matrix)
       CALL dbcsr_copy(density_matrix(1,1)%matrix, rho_ao(1,1)%matrix, keep_sparsity = .FALSE.)
       CALL dbcsr_set(density_matrix(ispin, 1)%matrix, 0.0_dp)
       CALL dbcsr_copy(density_matrix(1,1)%matrix, old_rho_ao(1,1)%matrix, keep_sparsity = .TRUE.)
    END DO

!    IF(para_env%mepos==0) THEN
!    WRITE(*,*) "ec_nsc_force - old_rho_ao"
!    CALL dbcsr_print(old_rho_ao(1,1)%matrix)
!
!    WRITE(*,*) "ec_nsc_force - rho_ao"
!    CALL dbcsr_print(rho_ao(1,1)%matrix)
!    
!    WRITE(*,*) "ec_nsc_force - density_matrix"
!    CALL dbcsr_print(density_matrix(1,1)%matrix)
!    END IF 

    !!---------------------------------------------------------------------------------------
        ! *** Calculation of forces due to non-self consistency of energy correction *** ! 
        ! \int v_nsc(r) * \grad_I(rho_in) 
    !!---------------------------------------------------------------------------------------
    DO ispin = 1 , nspins
       ! modifies KS matrix!
       CALL integrate_v_rspace(v_rspace = v_rspace_nsc(ispin), &
                               pmat = density_matrix(ispin, 1), &
                               hmat = ec_env%matrix_p(ispin, 1), &
                               qs_env = qs_env, &
                               calculate_forces = calculate_forces, &
                               basis_type = "HARRIS", &
                               task_list_external = ec_env%task_list)

!       CALL get_qs_env(qs_env, force = force)
!       WRITE(*,*) "ec_nsc_force - qs_env%...%force - rho_elec"
!       CALL ec_write_forces(qs_env, force)

    END DO


!    IF(para_env%mepos==0) THEN
!    WRITE(*,*) "ec_nsc_force - old_rho_ao - after integrate_v_rspace"
!    CALL dbcsr_print(old_rho_ao(1,1)%matrix)
!    WRITE(*,*) "ec_nsc_force - density_matrix after integrate_v_rspace"
!    CALL dbcsr_print(density_matrix(1,1)%matrix)
!    END IF


!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

    
!    IF(para_env%mepos==0) THEN
!    CALL get_qs_env(qs_env, force=force)
!    WRITE(*,*) "ec_nsc_force - qs_env%...%force - rho_elec - before copy"
!    CALL ec_write_forces(qs_env, force)
!    WRITE(*,*) "ec_force - ec_env%force - rho_elec - before copy"
!    CALL ec_write_forces(qs_env, ec_env%force)
!    END IF

    !---------------------------------------------------------------------------------------
    ! Need to copy NSC Hartree force contrib from qs_env%...%force to stored forces in ec_env%force
    ! acting on rho_elec 
    !---------------------------------------------------------------------------------------

    ! Substract NSC Hartree force contribution from saved forces in ec_env%force
    DO iatom = 1, natom
       ikind = kind_of(iatom)
       i = atom_of_kind(iatom)
       
       ec_env%force(ikind)%rho_elec(1:3, i) = ec_env%force(ikind)%rho_elec(1:3, i) - force(ikind)%rho_elec(1:3, i)
    END DO

    DEALLOCATE (atom_of_kind, kind_of)

    ! Delete NSC force contribution
    CALL zero_qs_force(force)
    
!    IF(para_env%mepos==0) THEN
!    WRITE(*,*) "ec_nsc_force - qs_env%...%force - after copy of rho_elec"
!    CALL ec_write_forces(qs_env, force)
!    WRITE(*,*) "ec_force - ec_env%force - after copy of rho_elec"
!    CALL ec_write_forces(qs_env, ec_env%force)
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
    
    !NULLIFY(density_matrix)
    DO ispin = 1 , nspins
    CALL dbcsr_deallocate_matrix_set(density_matrix)
    END DO

    IF(para_env%mepos==0) WRITE(*,*) "ec_nsc_force - END"
!    CALL ec_write_forces(qs_env, force, calculate_forces=.TRUE.)

    CALL timestop(handle)

  END SUBROUTINE ec_nsc_force

END MODULE ec_force
