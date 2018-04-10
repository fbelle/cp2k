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
                                              dbcsr_deallocate_matrix_set,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_print
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
                                              qs_force_type,&
                                              replicate_qs_force,&
                                              zero_qs_force
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_did_change,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
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
                                              pw_transfer
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
                                           
                                           
                                           
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_force'

! *** Public subroutines ***

   PUBLIC :: ec_forces, ec_nsc_force, ec_write_forces

CONTAINS


! **************************************************************************************************
!> \brief   Calculate the forces in the Harris energy correction.
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
    INTEGER                                            :: ispin, nspins
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_g, rho_r

    logger => cp_get_default_logger()
    
    CALL timeset(routineN//"_forces", handle)
    
    CALL get_qs_env(qs_env,para_env=para_env)
    IF(para_env%mepos==0) THEN
       WRITE(*,*) "ec_forces - BEGINNING"
    END IF

    ! obviously 
    calculate_forces = .TRUE.
    
    ! get qs and kg environment
    CALL get_qs_env(qs_env=qs_env, kg_env=kg_env)
    
    NULLIFY (para_env)
    CALL get_qs_env(qs_env, para_env=para_env)

    !shortcut to ec_type
    ! ec_build_neighborlist has been built during kg_energy_correction call
    ec_env => kg_env%ec_env

!    ! rebuild plane wave environment
!    CALL qs_env_rebuild_pw_env(qs_env)

      !------------------------------------------------------------------------
      !fbelle 
      CALL get_qs_env(qs_env,rho=rho)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
     
      WRITE(*,*) repeat("-",80)  
      WRITE(*,*) "ec_forces - BEFORE  - rho_ao "


      CALL qs_rho_get(rho, rho_ao_kp=rho_ao, rho_r=rho_r, rho_g=rho_g)
      WRITE(*,*) "rho_r%pw%in_use", rho_r(1)%pw%in_use
      WRITE(*,*) "rho_r(1)%pw%cr3d(1,1,1)", rho_r(1)%pw%cr3d(1,1,1)
      WRITE(*,*) "rho_g%pw%in_use", rho_g(1)%pw%in_use
      WRITE(*,*) "rho_g(1)%pw%cc(1)", rho_g(1)%pw%cc(1)
      WRITE(*,*) "rho_g(1)%pw%cc(2)", rho_g(1)%pw%cc(2)
      WRITE(*,*) "rho_g(1)%pw%cc(262144)", rho_g(1)%pw%cc(262144)


      WRITE(*,*) "ASSOCIATED(rho_ao)", ASSOCIATED(rho_ao)
      IF(ASSOCIATED(rho_ao)) THEN
          CALL dbcsr_print(rho_ao(1,1)%matrix)
      END IF
      WRITE(*,*) "ASSOCIATED(ec_env%matrix_p)", ASSOCIATED(ec_env%matrix_p)
      IF(ASSOCIATED(ec_env%matrix_p)) THEN
          CALL dbcsr_print(ec_env%matrix_p(1,1)%matrix)
      END IF
      
      WRITE(*,*) repeat("-",80)  
      !------------------------------------------------------------------------
   
    ! calculate matrix_w = matrix_p*matrix_ks*matrix_p
    IF(para_env%mepos==0) THEN 
       WRITE(*,*) "ec_force - ec_calculate_w_matrix"
    END IF
    CALL ec_calculate_w_matrix(qs_env, ec_env)

    ! *** d/dR[Sum of eigenvalues] = d/dR[trace(H_in * rho_out)] *** !
    !core hamiltonian force contribution: 
    ! kinetic, 
    ! overlap, 
    ! gth_ppl, 
    ! gth_ppnl, 
    !
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

    ! forces from double counting terms
    CALL ec_nsc_force(qs_env)

    CALL get_qs_env(qs_env, force=force)

    !  *** replicate forces ***
    ! sums up of current forces in total
    CALL replicate_qs_force(force, para_env)


    WRITE(*,*) "Forces END of ec_forces"  
    !CALL ec_write_forces(qs_env, calculate_forces)

    
    !update rho_ao with new density matrix of energy correction
    ! maybe move this to end of harris functional (after forces) so that rho_init and rho_final are available
    CALL get_qs_env(qs_env, dft_control=dft_control, rho=rho)
    CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
    nspins = dft_control%nspins
    DO ispin = 1, nspins
       CALL dbcsr_copy(rho_ao(ispin,1)%matrix, ec_env%matrix_p(ispin,1)%matrix ,keep_sparsity=.FALSE.)
    END DO

    
    ! updates rho
    CALL qs_rho_update_rho(rho, qs_env=qs_env)


      !------------------------------------------------------------------------
      !fbelle 
      CALL get_qs_env(qs_env,rho=rho)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
      
      WRITE(*,*) "ec_forces - END  - rho_ao "

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao, rho_r=rho_r, rho_g=rho_g)
      WRITE(*,*) "rho_r%pw%in_use", rho_r(1)%pw%in_use
      WRITE(*,*) "rho_r(1)%pw%cr3d(1,1,1)", rho_r(1)%pw%cr3d(1,1,1)
      WRITE(*,*) "rho_g%pw%in_use", rho_g(1)%pw%in_use
      WRITE(*,*) "rho_g(1)%pw%cc(1)", rho_g(1)%pw%cc(1)
      WRITE(*,*) "rho_g(1)%pw%cc(2)", rho_g(1)%pw%cc(2)
      WRITE(*,*) "rho_g(1)%pw%cc(262144)", rho_g(1)%pw%cc(262144)

      WRITE(*,*) "ASSOCIATED(rho_ao)", ASSOCIATED(rho_ao)
      IF(ASSOCIATED(rho_ao)) THEN
          CALL dbcsr_print(rho_ao(1,1)%matrix)
      END IF
      WRITE(*,*) "ASSOCIATED(ec_env%matrix_p)", ASSOCIATED(ec_env%matrix_p)
      IF(ASSOCIATED(ec_env%matrix_p)) THEN
          CALL dbcsr_print(ec_env%matrix_p(1,1)%matrix)
      END IF
      !------------------------------------------------------------------------


      ! deallocate W Matrix:
      CALL get_qs_env(qs_env=qs_env)
      CALL dbcsr_deallocate_matrix_set(ec_env%matrix_w)


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
    
    INTEGER                                            :: handle
    TYPE(cp_logger_type), POINTER                      :: logger
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
    ! fbelle
    TYPE(pw_env_type), POINTER                         :: pw_env
    TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                          v_hartree_rspace
    TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
    TYPE(pw_poisson_type), POINTER                     :: poisson_env


    logger => cp_get_default_logger()
    
    CALL timeset(routineN//"_forces", handle)
    
    CALL get_qs_env(qs_env, force=force, para_env=para_env)

    WRITE(*,*) "ec_nsc_force - BEGINNING"

     
    ! *** Calculation of forces due to non-self consistency of energy correction *** ! 
    ! First, v_hartree[rho_out-rho_in]

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

    
                        
    ! Get the total density in g-space [ions + electrons]
    ! rho contains old density info
    !CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)
                        
                        
    ! rho_in = rho_ao
    ! rho_out = ec_env%matrix_p 
    !CALL get_qs_env(qs_env,rho=rho)
    !CALL qs_rho_get(rho, rho_ao_kp=rho_ao)





    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
    CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)


    WRITE(*,*) "ec_nsc_force - END"
    !CALL ec_write_forces(qs_env, calculate_forces=.TRUE.)

    CALL timestop(handle)

  END SUBROUTINE ec_nsc_force




! **************************************************************************************************
!> \brief  Write output of forces 
!> \param qs_env ...
!> \date    03.2018
!> \author  Fabian Belleflamme
! **************************************************************************************************
  SUBROUTINE ec_write_forces(qs_env, calculate_forces)

    TYPE(qs_environment_type), POINTER                 :: qs_env
    LOGICAL                                            :: calculate_forces 

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_forces', routineP = moduleN//':'//routineN
    
    INTEGER                                            :: handle
    TYPE(cp_logger_type), POINTER                      :: logger
    TYPE(cp_para_env_type), POINTER                    :: para_env
    TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
    !fbelle
    CHARACTER(LEN=35)                                  :: fmtstr2
    REAL(KIND=dp), DIMENSION(3)                        :: grand_total
    INTEGER                                            :: i, iatom, ikind, natom, ndigits
    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
    TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set

    logger => cp_get_default_logger()
    
    CALL timeset(routineN//"_forces", handle)
    
    CALL get_qs_env(qs_env,para_env=para_env)

    !---------------------------------------------------------------------------------------
    IF(para_env%mepos==0) THEN
    
    IF (calculate_forces) THEN
    
       
       !needs this somehow
       CALL get_qs_env(qs_env=qs_env, force=force)
       
       ! get atom mapping
       NULLIFY (atomic_kind_set)
       CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set)
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, natom=natom)
       ALLOCATE (atom_of_kind(natom), kind_of(natom))
       CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                atom_of_kind=atom_of_kind, &
                                kind_of=kind_of)
       
       ndigits=6
       
       fmtstr2 = "(/,(T2,I5,4X,I4,T18,A,T34,3F  .  ))"
       WRITE (UNIT=fmtstr2(32:33), FMT="(I2)") ndigits
       WRITE (UNIT=fmtstr2(29:30), FMT="(I2)") ndigits+6
       
       !grand_total(:) = 0.0_dp
       
          DO iatom = 1, natom
             ikind = kind_of(iatom)
             i = atom_of_kind(iatom)
             WRITE (UNIT=6, FMT=fmtstr2) &
                iatom, ikind, "       overlap", force(ikind)%overlap(1:3, i), &
                iatom, ikind, "  overlap_admm", force(ikind)%overlap_admm(1:3, i), &
                iatom, ikind, "       kinetic", force(ikind)%kinetic(1:3, i), &
                iatom, ikind, "       gth_ppl", force(ikind)%gth_ppl(1:3, i), &
                iatom, ikind, "      gth_nlcc", force(ikind)%gth_nlcc(1:3, i), &
                iatom, ikind, "      gth_ppnl", force(ikind)%gth_ppnl(1:3, i), &
                iatom, ikind, "  core_overlap", force(ikind)%core_overlap(1:3, i), &
                iatom, ikind, "      rho_core", force(ikind)%rho_core(1:3, i), &
                iatom, ikind, "      rho_elec", force(ikind)%rho_elec(1:3, i), &
                iatom, ikind, "      rho_lri_elec", force(ikind)%rho_lri_elec(1:3, i), &
                iatom, ikind, "      ch_pulay", force(ikind)%ch_pulay(1:3, i), &
                iatom, ikind, "    dispersion", force(ikind)%dispersion(1:3, i), &
                iatom, ikind, "         other", force(ikind)%other(1:3, i), &
                iatom, ikind, "       fock_4c", force(ikind)%fock_4c(1:3, i), &
                iatom, ikind, "     ehrenfest", force(ikind)%ehrenfest(1:3, i), &
                iatom, ikind, "        efield", force(ikind)%efield(1:3, i), &
                iatom, ikind, "           eev", force(ikind)%eev(1:3, i), &
                iatom, ikind, "   mp2_non_sep", force(ikind)%mp2_non_sep(1:3, i), &
                iatom, ikind, "       mp2_sep", force(ikind)%mp2_sep(1:3, i), &
                iatom, ikind, "         total", force(ikind)%total(1:3, i)
             grand_total(1:3) = grand_total(1:3)+force(ikind)%total(1:3, i)
          END DO
    
    END IF
    END IF
    !---------------------------------------------------------------------------------------



    CALL timestop(handle)

  END SUBROUTINE ec_write_forces








END MODULE ec_force
