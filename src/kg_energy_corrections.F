!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for a Harris type energy correction on top of a
!>        Kim-Gordon calculation
!> \par History
!>       03.2014 created
!>       02.2018 Fabian Belleflamme - adding forces
!> \author JGH
! **************************************************************************************************
MODULE kg_energy_corrections
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_triangular_invert
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_triangular,&
                                              cp_fm_type, &
                                              cp_fm_to_fm,&
                                              cp_fm_p_type, &
                                              cp_fm_get_info,&
                                              cp_fm_write_formatted,&
                                              cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_allocate_matrix_set, dbcsr_copy, dbcsr_create, &
        dbcsr_deallocate_matrix_set, dbcsr_desymmetrize, dbcsr_distribution_type, dbcsr_filter, &
        dbcsr_get_info, dbcsr_init_p, dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_set, &
        dbcsr_trace, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric, dbcsr_print, dbcsr_finalize, &
        dbcsr_type_antisymmetric, dbcsr_scale
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE input_constants,                 ONLY: kg_ec_diagonalization,&
                                              kg_ec_functional_harris
   USE input_section_types,             ONLY: section_vals_val_get,&
                                              section_vals_type,&
                                              section_vals_get_subs_vals
   USE kg_environment_types,            ONLY: energy_correction_type,&
                                              kg_environment_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE mao_basis,                       ONLY: mao_generate_basis
   USE molecule_types,                  ONLY: molecule_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_integral_ab,&
                                              pw_integrate_function,&
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
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_dispersion_pairpot,           ONLY: calculate_dispersion_pairpot
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential,          ONLY: integrate_v_rspace, &
                                              integrate_v_core_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                               neighbor_list_set_p_type
   USE qs_neighbor_lists,               ONLY: atom2d_build,&
                                              atom2d_cleanup,&
                                              build_neighbor_lists,&
                                              local_atoms_type,&
                                              pair_radius_setup,&
                                              write_neighbor_lists
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE task_list_methods,               ONLY: generate_qs_task_list
   USE task_list_types,                 ONLY: allocate_task_list,&
                                              deallocate_task_list
   USE virial_types,                    ONLY: virial_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
!fbelle 
   USE qs_energy_utils,                 ONLY: qs_energies_compute_matrix_w
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              allocate_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE dm_ls_scf_qs,                    ONLY: matrix_ls_to_qs
   USE qs_rho_methods,                  ONLY: duplicate_rho_type,&
                                              qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE qs_initial_guess,                ONLY: calculate_atomic_block_dm
   USE ec_methods,                      ONLY: ec_write_forces
   USE qs_ks_types,                     ONLY: get_ks_env
   USE qs_energy_types,                 ONLY: init_qs_energy
   USE qs_core_energies,                ONLY: calculate_ptrace
   USE ec_methods,                      ONLY: ec_write_forces,&
                                              ec_write_bs_forces,& 
                                              ec_calculate_w_matrix,&
                                              ec_mos_init, &
                                              ec_rho_check, &
                                              matrix_projection
   USE input_section_types,             ONLY: section_vals_type
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE iterate_matrix,                  ONLY: invert_Hotelling


#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_energy_corrections'

   PUBLIC :: kg_energy_correction
   ! fbelle
   ! make this subr avail for external use 
   PUBLIC :: ec_build_core_hamiltonian
   PUBLIC :: ec_build_ks_matrix
   PUBLIC :: ec_ks_solver
   PUBLIC :: ec_energy
   PUBLIC :: ec_build_neighborlist
   PUBLIC :: mao_create_matrices
   PUBLIC :: ec_rho_update

CONTAINS

! **************************************************************************************************
!> \brief Initialize ec_env%old_rho 
!> 
!>        This routine needs to perform 3 task:
!>        - alloc & init ec_env%old_rho  
!>        - project rho_ao(KG_BASIS) to old_rho_ao(HARRIS_BASIS) 
!>        - update old_rho_r/-_g wrt old_rho_ao(HARRIS_BASIS) 
!>         
!> \param qs_env ...
!> \param ec_env ...
!> \date    07.2018
!> \author  fbelle
! **************************************************************************************************
   SUBROUTINE ec_rho_update(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
 
      CHARACTER(len=*), PARAMETER :: routineN = 'ec_rho_update3', &
         routineP = moduleN//':'//routineN
 
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: old_rho_ao, rho_ao
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: old_rho_r, old_rho_g
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      REAL(KIND=dp), DIMENSION(:), POINTER               :: old_tot_rho
      INTEGER                                            :: ispin, nspins
      TYPE(cp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp)                                      :: iorr
      REAL(KIND=dp)                                      :: irr
      ! basis set check
      INTEGER                                            :: ikind, nkind
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_kg
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_h
      LOGICAL                                            :: same_basis

      NULLIFY (rho)
      CALL get_qs_env(qs_env = qs_env, &
                      rho = rho, &
                      dft_control = dft_control, &
                      ks_env = ks_env)
      CALL get_qs_env(qs_env,para_env=para_env)
      
      IF(para_env%mepos==0) THEN
         WRITE(*,*) "ec_rho_update - BEGINNING"
      END IF
 
      nspins = dft_control%nspins

      !--------------------------------------------------
      ! init and allocate old_rho - to store input density
      !--------------------------------------------------
 
      ! allocate and initialize grids and density matrix
      CALL qs_rho_rebuild(rho = ec_env%old_rho, & 
                          qs_env = qs_env, &
                          rebuild_ao = .TRUE., &
                          rebuild_grids = .TRUE.)

      ! get sparse input density matrix
      CALL qs_rho_get(rho, rho_ao_kp = rho_ao, rho_r = rho_r )
      ! assign target to store old density
      CALL qs_rho_get(ec_env%old_rho, &
                      rho_r = old_rho_r, &
                      rho_g = old_rho_g, &
                      rho_ao_kp = old_rho_ao, &
                      tot_rho_r = old_tot_rho)

      IF(para_env%mepos==0) THEN
!         WRITE(*,*) "ec_rho_update - BEFORE  "
!         WRITE(*,*) "ec_rho_update3 - rho_ao", ASSOCIATED(rho_ao)
!         IF(ASSOCIATED(rho_ao)) THEN
!               CALL dbcsr_print(rho_ao(1,1)%matrix)
!         END IF
!         WRITE(*,*) "ec_rho_update3 - old_rho_ao", ASSOCIATED(old_rho_ao)
!         IF(ASSOCIATED(old_rho_ao)) THEN
!               CALL dbcsr_print(old_rho_ao(1,1)%matrix)
!         END IF
      END IF

      same_basis = .FALSE.

      CALL get_qs_env(qs_env = qs_env, qs_kind_set = qs_kind_set)
      nkind = SIZE(qs_kind_set)

      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind = qs_kind_set(ikind), basis_set = basis_set_h, basis_type = "HARRIS")
         CALL get_qs_kind(qs_kind = qs_kind_set(ikind), basis_set = basis_set_kg, basis_type = "ORB")

         IF (basis_set_kg%name == basis_set_h%name) THEN
            same_basis = .TRUE.
         ELSE
            same_basis = .FALSE.
         END IF
      END DO

      ! rebuild as non-sparse density matrix
      IF (ASSOCIATED(old_rho_ao)) CALL dbcsr_deallocate_matrix_set(old_rho_ao)
      CALL dbcsr_allocate_matrix_set(old_rho_ao, nspins, 1)
      CALL qs_rho_set(ec_env%old_rho, rho_ao_kp = old_rho_ao)

      DO ispin = 1, nspins
         ALLOCATE (old_rho_ao(ispin,1)%matrix)
         CALL dbcsr_create(matrix = old_rho_ao(ispin,1)%matrix, &
                           template = ec_env%matrix_s(1,1)%matrix, &
                           name = "OLD DENSITY MATRIX", &
                           matrix_type = dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(old_rho_ao(ispin,1)%matrix, ec_env%sab_orb)
         CALL dbcsr_set(old_rho_ao(ispin,1)%matrix, 0.0_dp)
      END DO


      !--------------------------------------------------
      ! Sparse density matrix projection
      ! KG BASIS -> HARRIS BASIS
      ! 1) same basis set size: copy sparse into non-sparse 
      ! 2) different basis set size: ADMM2 (see Guidon2010)
      !--------------------------------------------------
      IF(same_basis) THEN
         
         ! copy sparse density matrix rho_ao into non-sparse old_rho_ao
         DO ispin = 1, nspins
            CALL dbcsr_copy(old_rho_ao(ispin,1)%matrix, & 
                              rho_ao(ispin,1)%matrix, keep_sparsity=.TRUE.)
         END DO

      ELSE

         ! ADMM2 method of basis set tranformation
         CALL matrix_projection(qs_env, &
                                matrix_in = rho_ao, &
                                matrix_out = old_rho_ao, &
                                matrix_s_out = ec_env%matrix_s, &
                                sab_out_vs_in = ec_env%sab_ec_vs_kg, &
                                basis_type_in = "ORB", &
                                basis_type_out = "HARRIS")

      END IF

      DO ispin = 1, nspins
         ! computes the density corresponding to a given density matrix on the grid
         ! old_rho_ao -> old_rho_r & old_rho_g
         CALL calculate_rho_elec(matrix_p = old_rho_ao(ispin,1)%matrix, &
                                 rho = old_rho_r(ispin), &
                                 rho_gspace = old_rho_g(ispin), &
                                 total_rho = old_tot_rho(ispin), &
                                 ks_env = ks_env, &
                                 basis_type = "HARRIS", & 
                                 task_list_external = ec_env%task_list)
      END DO

      IF(para_env%mepos==0) THEN
!         WRITE(*,*) "ec_rho_update - rho_ao", ASSOCIATED(rho_ao)
!         IF(ASSOCIATED(rho_ao)) THEN
!               CALL dbcsr_print(rho_ao(1,1)%matrix)
!         END IF
!         WRITE(*,*) "ec_rho_update - old_rho_ao", ASSOCIATED(old_rho_ao)
!         IF(ASSOCIATED(old_rho_ao)) THEN
!               CALL dbcsr_print(old_rho_ao(1,1)%matrix)
!         END IF
      END IF

      !!---------------------------------------------------------------------------------------
      !!  Printing of density and integrate over grid
      !!---------------------------------------------------------------------------------------

      IF(para_env%mepos==0) THEN
         WRITE(*,*) "ec_rho_update - old_rho_ao after projection"
      END IF

      CALL ec_rho_check(qs_env, matrix_p = old_rho_ao, matrix_s = ec_env%matrix_s)

         ! density integrals over grid: \int rho dr = N_el
         iorr = pw_integrate_function(old_rho_r(1)%pw)
         irr = pw_integrate_function(rho_r(1)%pw)

         IF(para_env%mepos==0) THEN
            WRITE(*,*) " ec_rho_update - old_rho_r", iorr 
            WRITE(*,*) " ec_rho_update - rho_r", irr 
!            IF(ASSOCIATED(old_rho_ao)) THEN
!                  CALL dbcsr_print(old_rho_ao(1,1)%matrix)
!            END IF
         END IF

      !!---------------------------------------------------------------------------------------

      IF(para_env%mepos==0) THEN
         WRITE(*,*) " ec_rho_update - END"
      END IF

   END SUBROUTINE ec_rho_update

! **************************************************************************************************
!> \brief Energy correction to a KG simulation
!>
!> \param qs_env ...
!> \param calculate_forces ...
!> \par History
!>       03.2014 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE kg_energy_correction(qs_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: calculate_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'kg_energy_correction', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, unit_nr
      LOGICAL                                            :: my_calc_forces
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(kg_environment_type), POINTER                 :: kg_env
      ! 
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      my_calc_forces = .FALSE.
      ! fbelle
      ! disable so that CALL ec_build_core_hamiltonian only for energy here
      !IF (PRESENT(calculate_forces)) my_calc_forces = calculate_forces


      CALL get_qs_env(qs_env,para_env=para_env)
      IF(para_env%mepos==0) THEN
         WRITE(*,*) "my_calc_forces", my_calc_forces
         WRITE(*,*) "calculate_forces", calculate_forces 
      END IF

      NULLIFY (ec_env, kg_env)
      CALL get_qs_env(qs_env=qs_env, kg_env=kg_env, dft_control = dft_control)

      ! Check for energy correction
      !IF (kg_env%energy_correction) THEN
      IF (dft_control%qs_control%do_kg_ec) THEN

         ec_env => kg_env%ec_env

         ec_env%etotal = 0.0_dp
         ec_env%eband = 0.0_dp
         ec_env%ehartree = 0.0_dp
         ec_env%exc = 0.0_dp
         ec_env%vhxc = 0.0_dp
         ec_env%edispersion = 0.0_dp

         logger => cp_get_default_logger()
         IF (logger%para_env%mepos == logger%para_env%source) THEN
            unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
         ELSE
            unit_nr = -1
         ENDIF

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
            WRITE (unit_nr, '(T2,A,A,A,A,A)') "!", REPEAT("-", 27), " KG energy correction ", REPEAT("-", 28), "!"
         END IF


         ! build neighbor and task lists
         CALL ec_build_neighborlist(qs_env, ec_env)
         CALL ec_build_core_hamiltonian(qs_env, ec_env, my_calc_forces)

         ! store input density
         CALL ec_rho_update(qs_env, ec_env)

         CALL ec_build_ks_matrix(qs_env, ec_env, my_calc_forces)

         ! MAO basis
         IF (ec_env%mao) THEN
            IF (ASSOCIATED(ec_env%mao_coef)) CALL dbcsr_deallocate_matrix_set(ec_env%mao_coef)
            NULLIFY (ec_env%mao_coef)
            CALL mao_generate_basis(qs_env, ec_env%mao_coef, ref_basis_set="HARRIS", molecular=.TRUE., &
                                    max_iter=ec_env%mao_max_iter, eps_grad=ec_env%mao_eps_grad, unit_nr=unit_nr)
         END IF
         CALL ec_ks_solver(qs_env, ec_env)

         CALL ec_energy(qs_env, ec_env, unit_nr)

!         IF(para_env%mepos==0) THEN
!            WRITE(*,*) "kg_energy_correction - END - ec_env%matrix_ks", ASSOCIATED(ec_env%matrix_ks)
!            IF(ASSOCIATED(ec_env%matrix_ks)) THEN
!                  CALL dbcsr_print(ec_env%matrix_ks(1,1)%matrix)
!            END IF
!            WRITE(*,*) "kg_energy_correction - END - ec_env%matrix_p", ASSOCIATED(ec_env%matrix_p)
!            IF(ASSOCIATED(ec_env%matrix_ks)) THEN
!               CALL dbcsr_print(ec_env%matrix_p(1,1)%matrix)
!            END IF
!         END IF

!         IF (calculate_forces) THEN
!           CALL kg_response_solver(qs_env,ec_env)
!            CALL ec_build_core_hamiltonian(qs_env, ec_env, calculate_forces)
!         END IF
      
!         IF (calculate_forces) THEN
!            WRITE(*,*) "ec_env%matrix_p LEAVING kg_energy_correction"
!            WRITE(*,*) "ASSOCIATED(ec_env%matrix_p)", ASSOCIATED(ec_env%matrix_p)
!            CALL dbcsr_print(ec_env%matrix_p(1,1)%matrix)
!            
!            WRITE(*,*) "qs_env%ks_env%matrix_w  LEAVING kg_energy_correction"
!            CALL get_qs_env(qs_env,ks_env=ks_env,matrix_w=matrix_w)
!            WRITE(*,*) "ASSOCIATED(matrix_w)", ASSOCIATED(matrix_w)
!            CALL dbcsr_print(matrix_w(1)%matrix)
!         END IF


         ! Experimental
         ! copy Harris functional energy into qs_env%energy
         CALL get_qs_env(qs_env, energy = energy)
         CALL init_qs_energy(energy) 

         energy%total = ec_env%etotal


         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE kg_energy_correction

! **************************************************************************************************
!> \brief Construction of the Core Hamiltonian Matrix
!>        Short version of qs_core_hamiltonian
!> \param qs_env ...
!> \param ec_env ...
!> \param calculate_forces ...
!> \author Creation (03.2014,JGH)
! **************************************************************************************************
   SUBROUTINE ec_build_core_hamiltonian(qs_env, ec_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_build_core_hamiltonian', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nder, nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: eps_filter, eps_ppnl
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: kinetic


      CALL get_qs_env(qs_env,para_env=para_env)
      
      IF (calculate_forces) THEN
         CALL timeset(routineN//"_forces", handle)
         IF(para_env%mepos==0) THEN
            WRITE(*,*) "ec_build_core_hamiltonian - FORCES"
         END IF

      ELSE
         CALL timeset(routineN, handle)
      ENDIF




      ! no k-points possible
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      nimages = dft_control%nimages
      CPASSERT(nimages == 1)

      ! check for virial (currently no stress tensor available)
      CALL get_qs_env(qs_env=qs_env, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      CPASSERT(.NOT. use_virial)

      ! get neighbor lists, we need the full sab_orb list from the ec_env
      NULLIFY (sab_orb, sac_ppl, sap_ppnl)
      sab_orb => ec_env%sab_orb
      sac_ppl => ec_env%sac_ppl
      sap_ppnl => ec_env%sap_ppnl

      IF (calculate_forces) THEN
         CALL get_qs_env(qs_env=qs_env, force=force)
      END IF

      ! forces require first derivative of energy
      IF (calculate_forces) THEN
         nder = 1
      ELSE 
         nder = 0   
      END IF
      
      ! Overlap and kinetic energy matrices
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
      eps_filter = dft_control%qs_control%eps_filter_matrix
     
      IF (calculate_forces) THEN

        ! get rho_in - before Harris energy correction
        CALL get_qs_env(qs_env, matrix_s_kp=matrix_s)
        CALL get_qs_env(qs_env, kinetic_kp=kinetic)


        CALL build_overlap_matrix(ks_env, &
                                  nderivative = nder, &
                                  matrixkp_s = ec_env%matrix_s, &
                                  !matrixkp_s = matrix_s, &
                                  matrix_name = "OVERLAP MATRIX", &
                                  basis_type_a = "HARRIS", &
                                  basis_type_b = "HARRIS", &
                                  sab_nl = sab_orb, &
                                  calculate_forces = calculate_forces, &
                                  matrixkp_p = ec_env%matrix_w)

        CALL build_kinetic_matrix(ks_env, &
                                  matrixkp_t = ec_env%matrix_t, &
                                  !matrixkp_t = kinetic, &
                                  matrix_name = "KINETIC ENERGY MATRIX", &
                                  basis_type = "HARRIS", &
                                  sab_nl = sab_orb, &
                                  eps_filter = eps_filter, &
                                  calculate_forces = calculate_forces, &
                                  matrixkp_p = ec_env%matrix_p)
      ELSE
         
         CALL build_overlap_matrix(ks_env, &
                                   nderivative = nder, &
                                   matrixkp_s = ec_env%matrix_s, &
                                   matrix_name=  "OVERLAP MATRIX", &
                                   basis_type_a = "HARRIS", &
                                   basis_type_b = "HARRIS", &
                                   sab_nl = sab_orb)

         CALL build_kinetic_matrix(ks_env, &
                                   matrixkp_t = ec_env%matrix_t, &
                                   matrix_name = "KINETIC ENERGY MATRIX", &
                                   basis_type = "HARRIS", &
                                   sab_nl = sab_orb, &
                                   eps_filter = eps_filter)
      END IF

! ----------------------------------------------------------------------------- 

      IF(para_env%mepos==0) THEN
!         WRITE(*,*) " ec_build_core_hamiltonian ec_env%matrix_s", ASSOCIATED(ec_env%matrix_s)
!         IF(ASSOCIATED(ec_env%matrix_s)) THEN
!               CALL dbcsr_print(ec_env%matrix_s(1,1)%matrix)
!         END IF
!         WRITE(*,*) " ec_build_core_hamiltonian ec_env%matrix_p", ASSOCIATED(ec_env%matrix_p)
!         IF(ASSOCIATED(ec_env%matrix_p)) THEN
!               CALL dbcsr_print(ec_env%matrix_p(1,1)%matrix)
!         END IF
!         WRITE(*,*) " ec_build_core_hamiltonian ec_env%matrix_t", ASSOCIATED(ec_env%matrix_t)
!         IF(ASSOCIATED(ec_env%matrix_t)) THEN
!               CALL dbcsr_print(ec_env%matrix_t(1,1)%matrix)
!         END IF
      END IF

      ! initialize H matrix
      CALL dbcsr_allocate_matrix_set(ec_env%matrix_h, 1, 1)
      ALLOCATE (ec_env%matrix_h(1, 1)%matrix)
      CALL dbcsr_create(ec_env%matrix_h(1, 1)%matrix, template=ec_env%matrix_s(1, 1)%matrix)
      CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_h(1, 1)%matrix, sab_orb)

      ! add kinetic energy
      CALL dbcsr_copy(ec_env%matrix_h(1, 1)%matrix, ec_env%matrix_t(1, 1)%matrix, &
                      keep_sparsity=.TRUE., name="CORE HAMILTONIAN MATRIX")

      ! compute the ppl contribution to the core hamiltonian
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      
      NULLIFY (cell_to_index, virial)
      use_virial = .FALSE.

       IF (ASSOCIATED(sac_ppl)) THEN
          CALL build_core_ppl(ec_env%matrix_h, ec_env%matrix_p, force, & 
                              virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, & 
                              sab_orb, sac_ppl, nimages, cell_to_index, "HARRIS")
       END IF

       ! compute the ppnl contribution to the core hamiltonian ***
       eps_ppnl = dft_control%qs_control%eps_ppnl
       IF (ASSOCIATED(sap_ppnl)) THEN
          CALL build_core_ppnl(ec_env%matrix_h, ec_env%matrix_p, force, & 
                               virial, calculate_forces, use_virial, nder, &
                               qs_kind_set, atomic_kind_set, particle_set, &
                               sab_orb, sap_ppnl, eps_ppnl, nimages, cell_to_index, "HARRIS")
       END IF

      CALL timestop(handle)

   END SUBROUTINE ec_build_core_hamiltonian

! **************************************************************************************************
!> \brief calculate the complete KS matrix
!> \param qs_env ...
!> \param ec_env ...
!> \par History
!>      03.2014 adapted from qs_ks_build_kohn_sham_matrix [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ec_build_ks_matrix(qs_env, ec_env, calculate_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_build_ks_matrix', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: handle, ispin, nspins
      LOGICAL                                            :: use_virial
      REAL(dp)                                           :: eexc, ehartree, eovrl, eself, evhxc
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r, tau_r, v_rspace, v_tau_rspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, virial=virial, dft_control=dft_control, &
                      para_env=para_env, ks_env=ks_env)

      nspins = dft_control%nspins
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      CPASSERT(.NOT. use_virial)

      IF (.NOT. calculate_forces) THEN
         ! Kohn-Sham matrix
         IF (ASSOCIATED(ec_env%matrix_ks)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_ks)

         CALL dbcsr_allocate_matrix_set(ec_env%matrix_ks, nspins, 1)
         DO ispin = 1, nspins
            headline = "KOHN-SHAM MATRIX"
            ALLOCATE (ec_env%matrix_ks(ispin, 1)%matrix)
            CALL dbcsr_create(ec_env%matrix_ks(ispin, 1)%matrix, name=TRIM(headline), &
                              template=ec_env%matrix_s(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
            CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%sab_orb)
            CALL dbcsr_set(ec_env%matrix_ks(ispin, 1)%matrix, 0.0_dp)
         ENDDO
      END IF

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

      ! take reconstructed rho
      rho => ec_env%old_rho

      ! Get the total density in g-space [ions + electrons]
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

      ! Hartree potential
      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree, &
                            v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      NULLIFY (v_rspace, v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=ec_env%xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=eexc, just_energy=.FALSE.)

      evhxc = 0.0_dp

      CALL qs_rho_get(rho, rho_r=rho_r)
      IF (ASSOCIATED(v_tau_rspace)) THEN
         CALL qs_rho_get(rho, tau_r=tau_r)
      END IF

!      DO ispin = 1, nspins
!         evxc = evxc+pw_integral_ab(rho_r(ispin)%pw, v_rspace(ispin)%pw)
!      END DO

      DO ispin = 1, nspins
         ! Add v_hartree + v_xc = v_rspace
         CALL pw_scale(v_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(v_hartree_rspace%pw, v_rspace(ispin)%pw)
         ! integrate over potential <a|V|b>

         CALL integrate_v_rspace(v_rspace = v_rspace(ispin), &
                                 hmat = ec_env%matrix_ks(ispin, 1), &
                                 qs_env = qs_env, &
                                 calculate_forces = .FALSE., &
                                 basis_type = "HARRIS", &
                                 task_list_external = ec_env%task_list)
        
         IF (ASSOCIATED(v_tau_rspace)) THEN
            ! integrate over Tau-potential <nabla.a|V|nabla.b>
            CALL pw_scale(v_tau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            CALL integrate_v_rspace(v_rspace=v_tau_rspace(ispin), hmat=ec_env%matrix_ks(ispin, 1), &
                                    qs_env=qs_env, calculate_forces=calculate_forces, &
                                    compute_tau=.TRUE., &
                                    basis_type="HARRIS", &
                                    task_list_external=ec_env%task_list)
         END IF

         ! calclulate Int(vhxc*rho)dr and Int(vtau*tau)dr
         evhxc = evhxc+pw_integral_ab(rho_r(ispin)%pw, v_rspace(ispin)%pw)/v_rspace(1)%pw%pw_grid%dvol
         IF (ASSOCIATED(v_tau_rspace)) THEN
            evhxc = evhxc+pw_integral_ab(tau_r(ispin)%pw, & 
                v_tau_rspace(ispin)%pw)/v_tau_rspace(ispin)%pw%pw_grid%dvol
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
      CALL calculate_ecore_self(qs_env, E_self_core=eself)
      CALL calculate_ecore_overlap(qs_env, para_env, calculate_forces=calculate_forces, E_overlap_core=eovrl)
  
      CALL get_qs_env(qs_env,para_env=para_env)
      IF(para_env%mepos==0) THEN
          WRITE(*,*) "ec_build_ks_matrix - ehartree", ehartree
          WRITE(*,*) "ec_build_ks_matrix - erxc", eexc
!          !WRITE(*,*) "ehartree_core", ehartree_core
!          !WRITE(*,*) "ehartree_both", ehartree_both
!          WRITE(*,*) "evorl", eovrl
!          WRITE(*,*) "eself", eself
!          !WRITE(*,*) "evxc", evxc
!          WRITE(*,*) "evhxc", evhxc
!          WRITE(*,*) "Juerg correction",  ehartree - evhxc 
      END IF 
      ec_env%exc = eexc
      !ec_env%ehartree = -ehartree+eovrl+eself+ehartree_core
      ec_env%ehartree = ehartree+eovrl+eself
      ec_env%vhxc = evhxc
      !ec_env%vhxc = evxc

      IF (.NOT. calculate_forces) THEN
         ! add the core matrix
         DO ispin = 1, nspins
            CALL dbcsr_add(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%matrix_h(1, 1)%matrix, &
                           alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
         END DO

         ! At this point the ks matrix is up to date, filter it if requested
         DO ispin = 1, nspins
            CALL dbcsr_filter(ec_env%matrix_ks(ispin, 1)%matrix, &
                              dft_control%qs_control%eps_filter_matrix)
         ENDDO
      END IF


      CALL get_qs_env(qs_env,para_env=para_env)
      IF(para_env%mepos==0) THEN
!         WRITE(*,*) "ec_ks_solver - ec_env%matrix_h", ASSOCIATED(ec_env%matrix_h)
!         IF(ASSOCIATED(ec_env%matrix_h)) THEN
!               CALL dbcsr_print(ec_env%matrix_h(1,1)%matrix)
!         END IF
!         WRITE(*,*) "ec_ks_solver - ec_env%matrix_ks", ASSOCIATED(ec_env%matrix_ks)
!         IF(ASSOCIATED(ec_env%matrix_ks)) THEN
!               CALL dbcsr_print(ec_env%matrix_ks(1,1)%matrix)
!         END IF
      END IF

      DEALLOCATE (v_rspace)

      CALL timestop(handle)

   END SUBROUTINE ec_build_ks_matrix

! **************************************************************************************************
!> \brief Solve KS equation for a given matrix
!> \param qs_env ...
!> \param ec_env ...
!> \par History
!>      03.2014 created [JGH]
!> \author JGH
! **************************************************************************************************

   SUBROUTINE ec_ks_solver(qs_env, ec_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_ks_solver', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: handle, ispin, nspins
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, pmat, smat
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      ! create density matrix
      IF (.NOT. ASSOCIATED(ec_env%matrix_p)) THEN
         headline = "DENSITY MATRIX"
         CALL dbcsr_allocate_matrix_set(ec_env%matrix_p, nspins, 1)
         DO ispin = 1, nspins
            ALLOCATE (ec_env%matrix_p(ispin, 1)%matrix)
            CALL dbcsr_create(ec_env%matrix_p(ispin, 1)%matrix, name=TRIM(headline), &
                              template=ec_env%matrix_s(1, 1)%matrix)
            CALL cp_dbcsr_alloc_block_from_nbl(ec_env%matrix_p(ispin, 1)%matrix, ec_env%sab_orb)
         END DO
      END IF

      IF (ec_env%mao) THEN
         CALL mao_create_matrices(ec_env, ksmat, smat, pmat)
      ELSE
         ksmat => ec_env%matrix_ks
         smat => ec_env%matrix_s
         pmat => ec_env%matrix_p
      ENDIF

      SELECT CASE (ec_env%ks_solver)
      CASE (kg_ec_diagonalization)
         CALL ec_diag_solver(qs_env, ksmat, smat, pmat)
      CASE DEFAULT
         CPASSERT(.FALSE.)
      END SELECT

      IF (ec_env%mao) THEN
         CALL mao_release_matrices(ec_env, ksmat, smat, pmat)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE ec_ks_solver

! **************************************************************************************************
!> \brief Create matrices with MAO sizes
!> \param ec_env ...
!> \param ksmat ...
!> \param smat ...
!> \param pmat ...
!> \par History
!>      08.2016 created [JGH]
!> \author JGH
! **************************************************************************************************

   SUBROUTINE mao_create_matrices(ec_env, ksmat, smat, pmat)

      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, smat, pmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mao_create_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_sizes
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      TYPE(dbcsr_type)                                   :: cgmat

      CALL timeset(routineN, handle)

      mao_coef => ec_env%mao_coef

      NULLIFY (ksmat, smat, pmat)
      nspins = SIZE(ec_env%matrix_ks, 1)
      CALL dbcsr_get_info(mao_coef(1)%matrix, col_blk_size=col_blk_sizes, distribution=dbcsr_dist)
      CALL dbcsr_allocate_matrix_set(ksmat, nspins, 1)
      CALL dbcsr_allocate_matrix_set(smat, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (ksmat(ispin, 1)%matrix)
         CALL dbcsr_create(ksmat(ispin, 1)%matrix, dist=dbcsr_dist, name="MAO KS mat", &
                           matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                           col_blk_size=col_blk_sizes, nze=0)
         ALLOCATE (smat(ispin, 1)%matrix)
         CALL dbcsr_create(smat(ispin, 1)%matrix, dist=dbcsr_dist, name="MAO S mat", &
                           matrix_type=dbcsr_type_symmetric, row_blk_size=col_blk_sizes, &
                           col_blk_size=col_blk_sizes, nze=0)
      END DO
      !
      CALL dbcsr_create(cgmat, name="TEMP matrix", template=mao_coef(1)%matrix)
      DO ispin = 1, nspins
         CALL dbcsr_multiply("N", "N", 1.0_dp, ec_env%matrix_s(1, 1)%matrix, mao_coef(ispin)%matrix, &
                             0.0_dp, cgmat)
         CALL dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, smat(ispin, 1)%matrix)
         CALL dbcsr_multiply("N", "N", 1.0_dp, ec_env%matrix_ks(1, 1)%matrix, mao_coef(ispin)%matrix, &
                             0.0_dp, cgmat)
         CALL dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, ksmat(ispin, 1)%matrix)
      END DO
      CALL dbcsr_release(cgmat)

      CALL dbcsr_allocate_matrix_set(pmat, nspins, 1)
      DO ispin = 1, nspins
         ALLOCATE (pmat(ispin, 1)%matrix)
         CALL dbcsr_create(pmat(ispin, 1)%matrix, template=smat(1, 1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(pmat(ispin, 1)%matrix, ec_env%sab_orb)
      END DO

      CALL timestop(handle)

   END SUBROUTINE mao_create_matrices

! **************************************************************************************************
!> \brief Release matrices with MAO sizes
!> \param ec_env ...
!> \param ksmat ...
!> \param smat ...
!> \param pmat ...
!> \par History
!>      08.2016 created [JGH]
!> \author JGH
! **************************************************************************************************

   SUBROUTINE mao_release_matrices(ec_env, ksmat, smat, pmat)

      TYPE(energy_correction_type), POINTER              :: ec_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat, smat, pmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mao_release_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      TYPE(dbcsr_type)                                   :: cgmat

      CALL timeset(routineN, handle)

      mao_coef => ec_env%mao_coef
      nspins = SIZE(mao_coef, 1)

      ! save pmat in full basis format
      CALL dbcsr_create(cgmat, name="TEMP matrix", template=mao_coef(1)%matrix)
      DO ispin = 1, nspins
         CALL dbcsr_multiply("N", "N", 1.0_dp, mao_coef(ispin)%matrix, pmat(ispin, 1)%matrix, 0.0_dp, cgmat)
         CALL dbcsr_multiply("N", "T", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, &
                             ec_env%matrix_p(ispin, 1)%matrix, retain_sparsity=.TRUE.)
      END DO
      CALL dbcsr_release(cgmat)

      CALL dbcsr_deallocate_matrix_set(ksmat)
      CALL dbcsr_deallocate_matrix_set(smat)
      CALL dbcsr_deallocate_matrix_set(pmat)

      CALL timestop(handle)

   END SUBROUTINE mao_release_matrices

! **************************************************************************************************
!> \brief Solve KS equation using diagonalization
!> \param qs_env ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param matrix_p ...
!> \par History
!>      03.2014 created [JGH]
!> \author JGH
! **************************************************************************************************

!   SUBROUTINE ec_diag_solver(qs_env, matrix_ks, matrix_s, matrix_p)
!
!      TYPE(qs_environment_type), POINTER                 :: qs_env
!      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, matrix_p
!
!      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_diag_solver', routineP = moduleN//':'//routineN
!
!      INTEGER                                            :: handle, info, ispin, nmo(2), nsize, &
!                                                            nspins
!      REAL(KIND=dp)                                      :: eps_filter, focc(2)
!      !REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
!      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues
!      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
!      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
!      TYPE(cp_fm_type), POINTER                          :: fm_ks, fm_mo, fm_ortho
!      TYPE(cp_para_env_type), POINTER                    :: para_env
!      TYPE(dbcsr_type), POINTER                          :: buf1_dbcsr, buf2_dbcsr, ortho_dbcsr, &
!                                                            ref_matrix
!      TYPE(dft_control_type), POINTER                    :: dft_control
!
!      INTEGER                                            :: nao, homo
!!      TYPE(cp_blacs_env_type), POINTER                   :: context
!!      TYPE(cp_fm_type), POINTER                          :: ks, p, work
!      TYPE(cp_fm_type), POINTER                          :: weighted_vectors
!      TYPE(mo_set_type), POINTER                         :: mo_set
!      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
!
!      TYPE(energy_correction_type), POINTER              :: ec_env
!      TYPE(scf_control_type), POINTER                    :: scf_control
!      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
!
!      CALL timeset(routineN, handle)
!
!
!      ec_env => qs_env%kg_env%ec_env
!
!!      CALL ec_mos_init(qs_env, mos = ec_env%mos) 
!      CALL get_qs_env(qs_env, scf_control = scf_control,  mos = mos)
!
!      mo_set => mos(1)%mo_set
!
!      ! finish initialization of the MOs 
!      ! because in LS-SCF qs_scf_ensure_mos is not called
!      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
!      DO ispin = 1, SIZE(mos)
!         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=fm_mo)
!         IF (.NOT. ASSOCIATED(fm_mo)) THEN
!            CALL init_mo_set(mos(ispin)%mo_set, &
!                             fm_pool=ao_mo_fm_pools(ispin)%pool, &
!                             name="you get a name when you earn it")
!         END IF
!      END DO
!
!
!
!
!
!
!
!
!      CALL get_mo_set(mo_set = mo_set, &
!                      eigenvalues = eigenvalues, &
!                      mo_coeff=fm_mo)
!
!
!
!
!      NULLIFY (blacs_env, para_env)
!      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, para_env=para_env)
!
!      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
!      eps_filter = dft_control%qs_control%eps_filter_matrix
!      nspins = dft_control%nspins
!
!      nmo = 0
!      CALL get_qs_env(qs_env=qs_env, nelectron_spin=nmo)
!      focc = 1._dp
!      IF (nspins == 1) THEN
!         focc = 2._dp
!         nmo(1) = nmo(1)/2
!      END IF
!
!      CALL dbcsr_get_info(matrix_ks(1, 1)%matrix, nfullrows_total=nsize)
!      !ALLOCATE (eigenvalues(nsize))
!
!      NULLIFY (fm_ortho, fm_ks, fm_mo, fm_struct, ref_matrix)
!      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nsize, &
!                               ncol_global=nsize, para_env=para_env)
!      CALL cp_fm_create(fm_ortho, fm_struct)
!      CALL cp_fm_create(fm_ks, fm_struct)
!      CALL cp_fm_create(fm_mo, fm_struct)
!      CALL cp_fm_struct_release(fm_struct)
!
!      ! factorization
!      ref_matrix => matrix_s(1, 1)%matrix
!      NULLIFY (ortho_dbcsr, buf1_dbcsr, buf2_dbcsr)
!      CALL dbcsr_init_p(ortho_dbcsr)
!      CALL dbcsr_create(ortho_dbcsr, template=ref_matrix, &
!                        matrix_type=dbcsr_type_no_symmetry)
!      CALL dbcsr_init_p(buf1_dbcsr)
!      CALL dbcsr_create(buf1_dbcsr, template=ref_matrix, &
!                        matrix_type=dbcsr_type_no_symmetry)
!      CALL dbcsr_init_p(buf2_dbcsr)
!      CALL dbcsr_create(buf2_dbcsr, template=ref_matrix, &
!                        matrix_type=dbcsr_type_no_symmetry)
!
!      DO ispin = 1, nspins
!
!
!
!         ref_matrix => matrix_s(ispin, 1)%matrix
!         CALL copy_dbcsr_to_fm(ref_matrix, fm_ortho)
!         CALL cp_fm_cholesky_decompose(fm_ortho)
!         CALL cp_fm_triangular_invert(fm_ortho)
!         CALL cp_fm_set_all(fm_ks, 0.0_dp)
!         CALL cp_fm_to_fm_triangular(fm_ortho, fm_ks, "U")
!         CALL copy_fm_to_dbcsr(fm_ks, ortho_dbcsr)
!         CALL cp_fm_set_all(fm_ks, 0.0_dp)
!
!
!
!
!
!
!!         WRITE(*,*) "ref_matrix - matrix_s"
!!         CALL dbcsr_print(ref_matrix)
!!
!!         CALL copy_dbcsr_to_fm(ref_matrix, fm_ortho)
!!         CALL cp_fm_cholesky_decompose(fm_ortho)
!!
!!         CALL copy_fm_to_dbcsr(fm_ortho, ortho_dbcsr)
!!!         WRITE(*,*) "ortho_dbcsr"
!!!         CALL dbcsr_print(ortho_dbcsr)
!!!         CALL dbcsr_multiply("T", "N", 1.0_dp, ortho_dbcsr, ortho_dbcsr, &
!!!                             0.0_dp, buf1_dbcsr, filter_eps=eps_filter)
!!
!!         CALL cp_fm_write_formatted(fm_ortho,6)
!!
!!
!!         CALL cp_fm_triangular_invert(fm_ortho)
!!
!!         CALL cp_fm_write_formatted(fm_ortho,6)
!!
!!         CALL cp_fm_set_all(fm_ks, 0.0_dp)
!!         CALL cp_fm_to_fm_triangular(fm_ortho, fm_ks, "U")
!!         CALL copy_fm_to_dbcsr(fm_ks, ortho_dbcsr)
!!
!!         WRITE(*,*) "ortho_dbcsr"
!!         CALL dbcsr_print(ortho_dbcsr)
!!
!!         CALL dbcsr_set(buf1_dbcsr, 0.0_dp)
!!         CALL dbcsr_multiply("T", "N", 1.0_dp, ortho_dbcsr, ortho_dbcsr, &
!!                             0.0_dp, buf1_dbcsr, filter_eps=eps_filter)
!!
!!         CALL dbcsr_multiply("N", "T", 1.0_dp, ortho_dbcsr, ortho_dbcsr, &
!!                             0.0_dp, buf2_dbcsr, filter_eps=eps_filter)
!!
!!
!!         WRITE(*,*) "buf1_dbcsr"
!!         CALL dbcsr_print(buf1_dbcsr)
!!         WRITE(*,*) "buf2_dbcsr"
!!         CALL dbcsr_print(buf2_dbcsr)
!!
!!         CALL cp_fm_set_all(fm_ks, 0.0_dp)
!
!         ! calculate ZHZ(T)
!         ! calculate Z(T)HZ
!         CALL dbcsr_desymmetrize(matrix_ks(ispin, 1)%matrix, buf1_dbcsr)
!         CALL dbcsr_multiply("N", "N", 1.0_dp, buf1_dbcsr, ortho_dbcsr, &
!                             0.0_dp, buf2_dbcsr, filter_eps=eps_filter)
!         CALL dbcsr_multiply("T", "N", 1.0_dp, ortho_dbcsr, buf2_dbcsr, &
!                             0.0_dp, buf1_dbcsr, filter_eps=eps_filter)
!         ! copy to fm format
!         CALL copy_dbcsr_to_fm(buf1_dbcsr, fm_ks)
!         CALL choose_eigv_solver(fm_ks, fm_mo, eigenvalues, info)
!         CPASSERT(info == 0)
!         ! back transform of mos c = Z(T)*c
!         CALL copy_fm_to_dbcsr(fm_mo, buf1_dbcsr)
!         CALL dbcsr_multiply("N", "N", 1.0_dp, ortho_dbcsr, buf1_dbcsr, &
!                             0.0_dp, buf2_dbcsr, filter_eps=eps_filter)
!         ! density matrix
!         CALL dbcsr_set(matrix_p(ispin, 1)%matrix, 0.0_dp)
!         CALL dbcsr_multiply("N", "T", focc(ispin), buf2_dbcsr, buf2_dbcsr, &
!                             1.0_dp, matrix_p(ispin, 1)%matrix, retain_sparsity=.TRUE., last_k=nmo(ispin))
!
!      END DO
!!---------------------------------
!
!!      ! allocate ec_env%matrix_w
!!      CALL dbcsr_allocate_matrix_set(ec_env%matrix_w, nspins, 1)
!!      DO ispin = 1, nspins
!!         ALLOCATE (ec_env%matrix_w(ispin,1)%matrix)
!!         ! matrix_w same structure as matrix_s
!!         CALL dbcsr_copy(ec_env%matrix_w(ispin, 1)%matrix, ec_env%matrix_s(1, 1)%matrix)
!!         CALL dbcsr_set(ec_env%matrix_w(ispin, 1)%matrix, 0.0_dp)
!!      END DO
!!
!!      NULLIFY (context)
!!      NULLIFY (fm_struct)
!!      NULLIFY (ks)
!!      NULLIFY (p)
!!      NULLIFY (para_env)
!!      NULLIFY (work)
!!
!!      CALL cp_fm_get_info(fm_mo, context=context, nrow_global=nao, para_env=para_env)
!!      CALL cp_fm_struct_create(fm_struct, context=context, nrow_global=nao, &
!!                               ncol_global=nao, para_env=para_env)
!!      CALL cp_fm_create(ks, fm_struct, name="Kohn-Sham matrix")
!!      CALL cp_fm_create(p, fm_struct, name="Density matrix")
!!      CALL cp_fm_create(work, fm_struct, name="Work matrix")
!!      CALL cp_fm_struct_release(fm_struct)
!!      CALL copy_dbcsr_to_fm(matrix_ks(1,1)%matrix, ks)
!!      CALL copy_dbcsr_to_fm(matrix_p(1,1)%matrix, p)
!!      CALL cp_fm_upper_to_full(p, work)
!!      CALL cp_fm_symm("L", "U", nao, nao, 1.0_dp, ks, p, 0.0_dp, work)
!!      CALL cp_gemm("T", "N", nao, nao, nao, 1.0_dp, p, work, 0.0_dp, ks)
!!!      CALL dbcsr_set(ec_env%matrix_w(1,1)%matrix, 0.0_dp)
!!      CALL copy_fm_to_dbcsr(ks, ec_env%matrix_w(1,1)%matrix, keep_sparsity=.TRUE.)
!!
!!      CALL dbcsr_scale(ec_env%matrix_w(1,1)%matrix, 0.5)
!!      
!!
!!      CALL cp_fm_release(work)
!!      CALL cp_fm_release(p)
!!      CALL cp_fm_release(ks)
!!
!!      WRITE(*,*) "matrix_w"
!!      CALL dbcsr_print(ec_env%matrix_w(1,1)%matrix)
!
!      CALL set_mo_occupation(mo_array=mos, &
!                             smear=scf_control%smear)
!
!
!      CALL get_mo_set(mo_set = mo_set, &
!                      nao=nao, &
!                      homo=homo, &
!                      mo_coeff=fm_mo)
!
!      WRITE(*,*) "homo", homo
!      WRITE(*,*) "nao", nao
!
!
!      CALL cp_fm_create(weighted_vectors, mos(1)%mo_set%mo_coeff%matrix_struct, "weighted_vectors")
!      CALL cp_fm_to_fm(fm_mo, weighted_vectors)
!
!!      ! scale every column with the occupation
!!      ALLOCATE (eigocc(mo_set%homo))
!!
!!      DO imo = 1, mo_set%homo
!!         eigocc(imo) = mo_set%eigenvalues(imo)*mo_set%occupation_numbers(imo)
!!      ENDDO
!!      CALL cp_fm_column_scale(weighted_vectors, eigocc)
!!      DEALLOCATE (eigocc)
!!
!!      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=w_matrix, &
!!                                 matrix_v=mo_set%mo_coeff, &
!!                                 matrix_g=weighted_vectors, &
!!                                 ncol=mo_set%homo)
!
!
!!---------------------------------
!
!
!
!
!      CALL cp_fm_release(fm_ks)
!      CALL cp_fm_release(fm_mo)
!      CALL cp_fm_release(fm_ortho)
!      CALL dbcsr_release(ortho_dbcsr)
!      CALL dbcsr_release(buf1_dbcsr)
!      CALL dbcsr_release(buf2_dbcsr)
!      DEALLOCATE (ortho_dbcsr, buf1_dbcsr, buf2_dbcsr)
!      DEALLOCATE (eigenvalues)
!
!      CALL timestop(handle)
!
!   END SUBROUTINE ec_diag_solver




! **************************************************************************************************
!> \brief Solve KS equation using diagonalization
!> \param qs_env ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param matrix_p ...
!> \par History
!>      03.2014 created [JGH]
!> \author JGH
! **************************************************************************************************

   SUBROUTINE ec_diag_solver(qs_env, matrix_ks, matrix_s, matrix_p)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s, matrix_p

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_diag_solver', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, info, ispin, nmo(2), nsize, &
                                                            nspins
      REAL(KIND=dp)                                      :: eps_filter, focc(2)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fm_ks, fm_mo, fm_ortho
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type), POINTER                          :: buf1_dbcsr, buf2_dbcsr, ortho_dbcsr, &
                                                            ref_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, para_env)
      CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env, para_env=para_env)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)
      eps_filter = dft_control%qs_control%eps_filter_matrix
      nspins = dft_control%nspins

      nmo = 0
      CALL get_qs_env(qs_env=qs_env, nelectron_spin=nmo)
      focc = 1._dp
      IF (nspins == 1) THEN
         focc = 2._dp
         nmo(1) = nmo(1)/2
      END IF

      CALL dbcsr_get_info(matrix_ks(1, 1)%matrix, nfullrows_total=nsize)
      ALLOCATE (eigenvalues(nsize))

      NULLIFY (fm_ortho, fm_ks, fm_mo, fm_struct, ref_matrix)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=nsize, &
                               ncol_global=nsize, para_env=para_env)
      CALL cp_fm_create(fm_ortho, fm_struct)
      CALL cp_fm_create(fm_ks, fm_struct)
      CALL cp_fm_create(fm_mo, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! factorization
      ref_matrix => matrix_s(1, 1)%matrix
      NULLIFY (ortho_dbcsr, buf1_dbcsr, buf2_dbcsr)
      CALL dbcsr_init_p(ortho_dbcsr)
      CALL dbcsr_create(ortho_dbcsr, template=ref_matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(buf1_dbcsr)
      CALL dbcsr_create(buf1_dbcsr, template=ref_matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_init_p(buf2_dbcsr)
      CALL dbcsr_create(buf2_dbcsr, template=ref_matrix, &
                        matrix_type=dbcsr_type_no_symmetry)

      DO ispin = 1, nspins
         ref_matrix => matrix_s(ispin, 1)%matrix
         CALL copy_dbcsr_to_fm(ref_matrix, fm_ortho)
         CALL cp_fm_cholesky_decompose(fm_ortho)
         CALL cp_fm_triangular_invert(fm_ortho)
         CALL cp_fm_set_all(fm_ks, 0.0_dp)
         CALL cp_fm_to_fm_triangular(fm_ortho, fm_ks, "U")
         CALL copy_fm_to_dbcsr(fm_ks, ortho_dbcsr)
         CALL cp_fm_set_all(fm_ks, 0.0_dp)

         ! calculate ZHZ(T)
         ! calculate Z(T)HZ
         CALL dbcsr_desymmetrize(matrix_ks(ispin, 1)%matrix, buf1_dbcsr)
         CALL dbcsr_multiply("N", "N", 1.0_dp, buf1_dbcsr, ortho_dbcsr, &
                             0.0_dp, buf2_dbcsr, filter_eps=eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, ortho_dbcsr, buf2_dbcsr, &
                             0.0_dp, buf1_dbcsr, filter_eps=eps_filter)
         ! copy to fm format
         CALL copy_dbcsr_to_fm(buf1_dbcsr, fm_ks)
         CALL choose_eigv_solver(fm_ks, fm_mo, eigenvalues, info)
         CPASSERT(info == 0)
         ! back transform of mos c = Z(T)*c
         CALL copy_fm_to_dbcsr(fm_mo, buf1_dbcsr)
         CALL dbcsr_multiply("N", "N", 1.0_dp, ortho_dbcsr, buf1_dbcsr, &
                             0.0_dp, buf2_dbcsr, filter_eps=eps_filter)
         ! density matrix
         CALL dbcsr_set(matrix_p(ispin, 1)%matrix, 0.0_dp)
         CALL dbcsr_multiply("N", "T", focc(ispin), buf2_dbcsr, buf2_dbcsr, &
                             1.0_dp, matrix_p(ispin, 1)%matrix, retain_sparsity=.TRUE., last_k=nmo(ispin))
      END DO

      CALL cp_fm_release(fm_ks)
      CALL cp_fm_release(fm_mo)
      CALL cp_fm_release(fm_ortho)
      CALL dbcsr_release(ortho_dbcsr)
      CALL dbcsr_release(buf1_dbcsr)
      CALL dbcsr_release(buf2_dbcsr)
      DEALLOCATE (ortho_dbcsr, buf1_dbcsr, buf2_dbcsr)
      DEALLOCATE (eigenvalues)

      CALL timestop(handle)

   END SUBROUTINE ec_diag_solver

! **************************************************************************************************
!> \brief Calculate the energy correction
!> \param qs_env ...
!> \param ec_env ...
!> \param unit_nr ...
!> \author Creation (03.2014,JGH)
! **************************************************************************************************
   SUBROUTINE ec_energy(qs_env, ec_env, unit_nr)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type)                       :: ec_env
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_energy', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nspins
      REAL(KIND=dp)                                      :: eband, energy, trace
      !fbelle
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_h
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_t
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_ks
!      REAL(KIND=dp)                                      :: trace_sp_ii, trace_sp_oo
!      REAL(KIND=dp)                                      :: trace_sp_io, trace_sp_oi
!
!      REAL(KIND=dp)                                      :: ecore_ii, ecore_oo
!      REAL(KIND=dp)                                      :: ecore_io, ecore_oi
!      
!      REAL(KIND=dp)                                      :: ekinetic_ii, ekinetic_oo
!      REAL(KIND=dp)                                      :: ekinetic_io, ekinetic_oi
!      
!      REAL(KIND=dp)                                      :: ebs_ii, ebs_oo
!      REAL(KIND=dp)                                      :: ebs_io, ebs_oi

      CALL timeset(routineN, handle)

      ! dispersion through pairpotentials
      CALL calculate_dispersion_pairpot(qs_env, ec_env%dispersion_env, energy, .FALSE.)
      ec_env%edispersion = ec_env%edispersion+energy

      SELECT CASE (ec_env%energy_functional)
      CASE (kg_ec_functional_harris)

         nspins = SIZE(ec_env%matrix_ks, 1)


      ! fbelle - old ks matrix
!      CALL get_qs_env(qs_env=qs_env, matrix_ks_kp=matrix_ks)     


         eband = 0.0_dp
         DO ispin = 1, nspins
!dbg
            CALL dbcsr_trace(ec_env%matrix_p(ispin, 1)%matrix, ec_env%matrix_s(1, 1)%matrix, trace)
            IF (unit_nr > 0) WRITE (unit_nr, '(T2,A,T16,F16.10)') 'Tr[PS] ', trace
!dbg
! fbelle - original
            CALL dbcsr_trace(ec_env%matrix_ks(ispin, 1)%matrix, ec_env%matrix_p(ispin, 1)%matrix, trace)
            eband = eband+trace
!            CALL calculate_ptrace(matrix_ks, ec_env%matrix_p, trace, nspins)
!            eband = eband+trace
         END DO

         ec_env%eband = eband
         ec_env%etotal = ec_env%eband+ec_env%ehartree+ec_env%exc-ec_env%vhxc+ec_env%edispersion
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "HF Etotal   ", ec_env%etotal
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Eband    ", ec_env%eband
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Ehartree ", ec_env%ehartree
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Exc      ", ec_env%exc
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Evhxc    ", ec_env%vhxc
            WRITE (unit_nr, '(T2,A,T16,F16.10)') "Edisp    ", ec_env%edispersion
         END IF

      CASE DEFAULT

         CPASSERT(.FALSE.)

      END SELECT

!      ! fbelle 
      CALL get_qs_env(qs_env, rho=rho, para_env = para_env)
      CALL qs_rho_get(rho, rho_ao_kp = rho_ao)
      CALL get_qs_env(qs_env=qs_env, matrix_ks_kp = matrix_ks, matrix_s_kp = matrix_s)     
      CALL get_qs_env(qs_env=qs_env, matrix_h_kp = matrix_h)     
      CALL get_qs_env(qs_env=qs_env, kinetic_kp = matrix_t)     


      IF(para_env%mepos==0) THEN
!         WRITE(*,*) " ec_energy matrix_ks", ASSOCIATED(matrix_ks)
!         IF(ASSOCIATED(matrix_ks)) THEN
!               CALL dbcsr_print(matrix_ks(1,1)%matrix)
!         END IF
!         WRITE(*,*) " ec_energy ec_env%matrix_ks", ASSOCIATED(ec_env%matrix_ks)
!         IF(ASSOCIATED(ec_env%matrix_ks)) THEN
!               CALL dbcsr_print(ec_env%matrix_ks(1,1)%matrix)
!         END IF
!         WRITE(*,*) " ec_energy ec_env%matrix_t", ASSOCIATED(ec_env%matrix_t)
!         IF(ASSOCIATED(ec_env%matrix_t)) THEN
!               CALL dbcsr_print(ec_env%matrix_t(1,1)%matrix)
!         END IF
!         WRITE(*,*) " ec_energy ec_env%matrix_p", ASSOCIATED(ec_env%matrix_p)
!         IF(ASSOCIATED(ec_env%matrix_p)) THEN
!               CALL dbcsr_print(ec_env%matrix_p(1,1)%matrix)
!         END IF
!         WRITE(*,*) " ec_energy -  rho_ao", ASSOCIATED(rho_ao)
!         IF(ASSOCIATED(rho_ao)) THEN
!               CALL dbcsr_print(rho_ao(1,1)%matrix)
!         END IF
      END IF
      
      ! core matrix trace
!      CALL calculate_ptrace(ec_env%matrix_h, ec_env%matrix_p, ecore_oo, nspins)
!      CALL calculate_ptrace(ec_env%matrix_h, rho_ao, ecore_oi, nspins)

!      CALL calculate_ptrace(matrix_h, ec_env%matrix_p, ecore_io, nspins)
!      CALL calculate_ptrace(matrix_h, rho_ao, ecore_ii, nspins)

      ! kinetic matrix trace
!      CALL calculate_ptrace(ec_env%matrix_t, ec_env%matrix_p, ekinetic_oo, nspins)
!      CALL calculate_ptrace(ec_env%matrix_t, rho_ao, ekinetic_oi, nspins)
      
!      CALL calculate_ptrace(matrix_t, ec_env%matrix_p, ekinetic_io, nspins)
!      CALL calculate_ptrace(matrix_t, rho_ao, ekinetic_ii, nspins)

      ! overlap
!      CALL calculate_ptrace(ec_env%matrix_s, ec_env%matrix_p, trace_sp_oo, nspins)
!      CALL calculate_ptrace(ec_env%matrix_s, rho_ao, trace_sp_oi, nspins)

!      CALL calculate_ptrace(matrix_s, ec_env%matrix_p, trace_sp_io, nspins)
!      CALL calculate_ptrace(matrix_s, rho_ao, trace_sp_ii, nspins)
      
      ! band structure energy = Tr[KS*P]
!      CALL calculate_ptrace(ec_env%matrix_ks, ec_env%matrix_p, ebs_oo, nspins)
!      CALL calculate_ptrace(ec_env%matrix_ks, rho_ao, ebs_oi, nspins)
      
!      CALL calculate_ptrace(matrix_ks, ec_env%matrix_p, ebs_io , nspins)
!      CALL calculate_ptrace(matrix_ks, rho_ao, ebs_ii, nspins)
      
      CALL get_qs_env(qs_env,para_env=para_env)
      IF(para_env%mepos==0) THEN
         !core
!         WRITE(*,*) "trace[h_out*p_out]",  ecore_oo
!         WRITE(*,*) "trace[h_out*p_in]",   ecore_oi
!         WRITE(*,*) "trace[h_in*p_out]",   ecore_io
!         WRITE(*,*) "trace[h_in*p_in]",    ecore_ii
         
         ! kinetic 
!         WRITE(*,*) "trace[T_out*p_out]",  ekinetic_oo
!         WRITE(*,*) "trace[T_out*p_in]",   ekinetic_oi
!         WRITE(*,*) "trace[T_in*p_out]",   ekinetic_io
!         WRITE(*,*) "trace[T_in*p_in]",    ekinetic_ii
         
!         WRITE(*,*) "trace[ks_out*p_out]", ebs_oo
!         WRITE(*,*) "trace[ks_out*p_in]",  ebs_oi
!         WRITE(*,*) "trace[ks_in*p_out]",  ebs_io
!         WRITE(*,*) "trace[ks_in*p_in]",   ebs_ii

!         WRITE(*,*) "trace[S_in*p_in]",    trace_sp_ii
!         WRITE(*,*) "trace[S_out*p_out]",  trace_sp_oo
!         WRITE(*,*) "trace[S_out*p_in]",   trace_sp_oi
!         WRITE(*,*) "trace[S_in*p_out]",   trace_sp_io
      END IF

      CALL timestop(handle)

   END SUBROUTINE ec_energy

! **************************************************************************************************
!> \brief builds either the full neighborlist or neighborlists of molecular
!> \brief subsets, depending on parameter values
!> \param qs_env ...
!> \param ec_env ...
!> \par History
!>       2012.07 created [Martin Haeufel]
!>       2016.07 Adapted for Harris functional [JGH]
!> \author Martin Haeufel
! **************************************************************************************************
   SUBROUTINE ec_build_neighborlist(qs_env, ec_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(energy_correction_type), POINTER              :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_build_neighborlist', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikind, nkind
      LOGICAL                                            :: gth_potential_present, &
                                                            sgp_potential_present, &
                                                            skip_load_balance_distributed
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: orb_present, ppl_present, ppnl_present
      REAL(dp)                                           :: subcells
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: orb_radius, ppl_radius, ppnl_radius
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: pair_radius
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: distribution_1d
      TYPE(distribution_2d_type), POINTER                :: distribution_2d
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(local_atoms_type), ALLOCATABLE, DIMENSION(:)  :: atom2d
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential
      !fbelle
!      TYPE(section_vals_type), POINTER                   :: neighbor_list_section
!      TYPE(cp_para_env_type), POINTER                    :: para_env
      ! fbelle
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_kg
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: kg_radius
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: kg_present



      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set, gth_potential_present=gth_potential_present, &
                           sgp_potential_present=sgp_potential_present)
      nkind = SIZE(qs_kind_set)
      ALLOCATE (orb_radius(nkind), ppl_radius(nkind), ppnl_radius(nkind))
      ALLOCATE (orb_present(nkind), ppl_present(nkind), ppnl_present(nkind))
      ALLOCATE (pair_radius(nkind, nkind))
      ALLOCATE (atom2d(nkind))

      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      cell=cell, &
                      distribution_2d=distribution_2d, &
                      local_particles=distribution_1d, &
                      particle_set=particle_set, &
                      molecule_set=molecule_set)

      CALL atom2d_build(atom2d, distribution_1d, distribution_2d, atomic_kind_set, &
                        molecule_set, .FALSE., particle_set)

      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom2d(ikind)%list)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set, basis_type="HARRIS")
         IF (ASSOCIATED(basis_set)) THEN
            orb_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=basis_set, kind_radius=orb_radius(ikind))
         ELSE
            orb_present(ikind) = .FALSE.
            orb_radius(ikind) = 0.0_dp
         END IF
         CALL get_qs_kind(qs_kind, gth_potential=gth_potential, sgp_potential=sgp_potential)
         IF (gth_potential_present .OR. sgp_potential_present) THEN
            IF (ASSOCIATED(gth_potential)) THEN
               CALL get_potential(potential=gth_potential, &
                                  ppl_present=ppl_present(ikind), &
                                  ppl_radius=ppl_radius(ikind), &
                                  ppnl_present=ppnl_present(ikind), &
                                  ppnl_radius=ppnl_radius(ikind))
            ELSE IF (ASSOCIATED(sgp_potential)) THEN
               CALL get_potential(potential=sgp_potential, &
                                  ppl_present=ppl_present(ikind), &
                                  ppl_radius=ppl_radius(ikind), &
                                  ppnl_present=ppnl_present(ikind), &
                                  ppnl_radius=ppnl_radius(ikind))
            ELSE
               ppl_present(ikind) = .FALSE.
               ppl_radius(ikind) = 0.0_dp
               ppnl_present(ikind) = .FALSE.
               ppnl_radius(ikind) = 0.0_dp
            END IF
         END IF
      END DO

      CALL section_vals_val_get(qs_env%input, "DFT%SUBCELLS", r_val=subcells)

      ! overlap
      CALL pair_radius_setup(orb_present, orb_present, orb_radius, orb_radius, pair_radius)
      CALL build_neighbor_lists(ec_env%sab_orb, particle_set, atom2d, cell, pair_radius, &
                                subcells=subcells, nlname="sab_orb")
      ! pseudopotential
      IF (gth_potential_present .OR. sgp_potential_present) THEN
         IF (ANY(ppl_present)) THEN
            CALL pair_radius_setup(orb_present, ppl_present, orb_radius, ppl_radius, pair_radius)
            CALL build_neighbor_lists(ec_env%sac_ppl, particle_set, atom2d, cell, pair_radius, &
                                      subcells=subcells, operator_type="ABC", nlname="sac_ppl")
         END IF

         IF (ANY(ppnl_present)) THEN
            CALL pair_radius_setup(orb_present, ppnl_present, orb_radius, ppnl_radius, pair_radius)
            CALL build_neighbor_lists(ec_env%sap_ppnl, particle_set, atom2d, cell, pair_radius, &
                                      subcells=subcells, operator_type="ABBA", nlname="sap_ppnl")
         END IF
      END IF

 
      ! fbelle
      ! Build mixed overlap matrix KG-BASIS vs HARRIS-BASIS (ORB here)
      ALLOCATE (kg_present(nkind))
      ALLOCATE (kg_radius(nkind))
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom2d(ikind)%list)
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_kg, basis_type="ORB")
         IF (ASSOCIATED(basis_set)) THEN
            kg_present(ikind) = .TRUE.
            CALL get_gto_basis_set(gto_basis_set=basis_set_kg, kind_radius=kg_radius(ikind))
         ELSE
            kg_present(ikind) = .FALSE.
            kg_radius(ikind) = 0.0_dp
         END IF
      END DO

      CALL pair_radius_setup(orb_present, kg_present, orb_radius, kg_radius, pair_radius)
      CALL build_neighbor_lists(ec_env%sab_ec_vs_kg, particle_set, atom2d, cell, pair_radius, &
                                symmetric=.FALSE.,  molecular=.FALSE., &
                                subcells=subcells, nlname="sab_ec_vs_kg")



      ! Release work storage
      CALL atom2d_cleanup(atom2d)
      DEALLOCATE (atom2d)
      DEALLOCATE (orb_present, ppl_present, ppnl_present)
      DEALLOCATE (orb_radius, ppl_radius, ppnl_radius)
      DEALLOCATE (pair_radius)
      !fbelle
      DEALLOCATE(kg_present)
      DEALLOCATE (kg_radius)


      ! Task list
      CALL get_qs_env(qs_env, ks_env=ks_env, dft_control=dft_control)
      skip_load_balance_distributed = dft_control%qs_control%skip_load_balance_distributed
      IF (ASSOCIATED(ec_env%task_list)) CALL deallocate_task_list(ec_env%task_list)
      CALL allocate_task_list(ec_env%task_list)
      CALL generate_qs_task_list(ks_env, ec_env%task_list, &
                                 reorder_rs_grid_ranks=.FALSE., soft_valid=.FALSE., &
                                 skip_load_balance_distributed=skip_load_balance_distributed, &
                                 basis_type="HARRIS", sab_orb_external=ec_env%sab_orb)




      ! fbelle 
      ! Print out the neighborlist
      !neighbor_list_section => section_vals_get_subs_vals(qs_env%input, "DFT%KG_METHOD%PRINT%NEIGHBOR_LISTS")

      !CALL get_qs_env(qs_env, para_env = para_env)

      !CALL write_neighbor_lists(ec_env%sab_orb, particle_set, cell, para_env, neighbor_list_section, &
      !                         "/SAB_ORB_FULL", "sab_orb", "HARRIS NEIGHBORLIST")


      CALL timestop(handle)

   END SUBROUTINE ec_build_neighborlist

END MODULE kg_energy_corrections
