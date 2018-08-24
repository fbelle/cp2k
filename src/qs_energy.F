!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Perform a QUICKSTEP wavefunction optimization (single point)
!> \par History
!>      none
!> \author MK (29.10.2002)
! **************************************************************************************************
MODULE qs_energy
   USE almo_scf,                        ONLY: almo_entry_scf
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE dm_ls_scf,                       ONLY: ls_scf
   USE kg_energy_corrections,           ONLY: kg_energy_correction
   USE lri_environment_methods,         ONLY: lri_print_stat
   USE qs_energy_init,                  ONLY: qs_energies_init
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_energy_utils,                 ONLY: qs_energies_compute_matrix_w,&
                                              qs_energies_mp2,&
                                              qs_energies_properties
   USE qs_environment_methods,          ONLY: qs_env_rebuild_pw_env
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
   USE qs_scf,                          ONLY: scf
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize,&
        dbcsr_distribution_type, dbcsr_filter,  dbcsr_get_info,&
        dbcsr_init_p, dbcsr_multiply, dbcsr_p_type, dbcsr_release,&
        dbcsr_set, dbcsr_trace,  dbcsr_type, dbcsr_type_no_symmetry,&
        dbcsr_type_symmetric, dbcsr_print, dbcsr_finalize
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_core_energies,                ONLY: calculate_ptrace
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE ec_methods,                      ONLY: ec_rho_check   
   USE kg_environment_types,            ONLY: energy_correction_type,&
                                              kg_environment_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_energy'

   PUBLIC :: qs_energies

CONTAINS

! **************************************************************************************************
!> \brief   Driver routine for QUICKSTEP single point wavefunction optimization.
!> \param qs_env ...
!> \param consistent_energies ...
!> \param calc_forces ...
!> \date    29.10.2002
!> \par History
!>          - consistent_energies option added (25.08.2005, TdK)
!>          - introduced driver for energy in order to properly decide between
!>            SCF or RTP (fschiff 02.09)
!> \author  MK
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE qs_energies(qs_env, consistent_energies, calc_forces)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: consistent_energies, calc_forces

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_energies', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: my_calc_forces, run_rtp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      !fbelle
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: rho_ao
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_ks
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: kinetic
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER        :: matrix_w
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_h
      !fbelle - trace
!      INTEGER                                            :: nspins
!      REAL(KIND=dp)                                      :: trace_ks_p, ecore, ekinetic, trace_ovrl
!      REAL(KIND=dp)                                      :: trace
      TYPE(kg_environment_type), POINTER                 :: kg_env

      CALL timeset(routineN, handle)

      my_calc_forces = .FALSE.
      IF (PRESENT(calc_forces)) my_calc_forces = calc_forces

      CALL qs_env_rebuild_pw_env(qs_env)

      CALL get_qs_env(qs_env=qs_env, run_rtp=run_rtp)
      IF (.NOT. run_rtp) THEN

         NULLIFY (dft_control, energy)
         CALL qs_energies_init(qs_env, my_calc_forces)
         CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, energy=energy)

         ! *** Perform a SCF run ***
         IF (dft_control%qs_control%do_ls_scf) THEN
            CALL ls_scf(qs_env=qs_env)
         ELSE IF (dft_control%qs_control%do_almo_scf) THEN
            CALL almo_entry_scf(qs_env=qs_env, calc_forces=my_calc_forces)
         ELSE
            CALL scf(qs_env=qs_env)

            ! Compute MP2 energy
            CALL qs_energies_mp2(qs_env, my_calc_forces)

            ! if calculate forces, time to compute the w matrix
            CALL qs_energies_compute_matrix_w(qs_env, my_calc_forces)

         END IF

      !--------------------------------------------------------------------------
 
      CALL get_qs_env(qs_env,para_env=para_env)
      IF(para_env%mepos==0) THEN
        WRITE(*,*) "qs_energies - after KG" 
        
        WRITE(*,*) "core        "  , energy%core
        WRITE(*,*) "core_overlap " , energy%core_overlap
        WRITE(*,*) "core_self   "  , energy%core_self
        WRITE(*,*) "exc         "  , energy%exc
        WRITE(*,*) "hartree     "  , energy%hartree
        WRITE(*,*) "e_hartree,  "  , energy%e_hartree
        WRITE(*,*) "total       "  , energy%total
        WRITE(*,*) "tot_old     "  , energy%tot_old
        WRITE(*,*) "kinetic     "  , energy%kinetic
     END IF
 
      ! fbelle 
      CALL get_qs_env(qs_env, rho=rho, para_env = para_env)
      CALL qs_rho_get(rho, rho_ao_kp = rho_ao)
      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks)     
      CALL get_qs_env(qs_env, matrix_w=matrix_w)
      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s)
      CALL get_qs_env(qs_env, matrix_h_kp=matrix_h)
      CALL get_qs_env(qs_env, kinetic_kp=kinetic)
      CALL get_qs_env(qs_env, kg_env = kg_env)

      CALL ec_rho_check(qs_env, matrix_p = rho_ao, matrix_s = matrix_s)

      IF(para_env%mepos==0) THEN
!         WRITE(*,*) " qs_energies - matrix_ks", ASSOCIATED(matrix_ks)
!         IF(ASSOCIATED(matrix_ks)) THEN
!               CALL dbcsr_print(matrix_ks(1,1)%matrix)
!         END IF
!         WRITE(*,*) " qs_energies - matrix_s", ASSOCIATED(matrix_s)
!         IF(ASSOCIATED(matrix_s)) THEN
!               CALL dbcsr_print(matrix_s(1,1)%matrix)
!         END IF
!         WRITE(*,*) " qs_energies - matrix_h", ASSOCIATED(matrix_h)
!         IF(ASSOCIATED(matrix_h)) THEN
!               CALL dbcsr_print(matrix_h(1,1)%matrix)
!         END IF
!         WRITE(*,*) " qs_energies - rho_ao", ASSOCIATED(rho_ao)
!         IF(ASSOCIATED(rho_ao)) THEN
!               CALL dbcsr_print(rho_ao(1,1)%matrix)
!         END IF
!         WRITE(*,*) " qs_energies - kinetic", ASSOCIATED(kinetic)
!         IF(ASSOCIATED(kinetic)) THEN
!               CALL dbcsr_print(kinetic(1,1)%matrix)
!         END IF
!         WRITE(*,*) " qs_energies - matrix_w", ASSOCIATED(matrix_w)
!         IF(ASSOCIATED(matrix_w)) THEN
!               CALL dbcsr_print(matrix_w(1)%matrix)
!         END IF
!         CALL calculate_ptrace(matrix_ks, rho_ao, trace, dft_control%nspins)
!         WRITE(*,*) "trace[ks_out*p_out]", trace

         WRITE(*,*) " qs_energies - tnadd_mat", ASSOCIATED(kg_env%tnadd_mat)
         IF(ASSOCIATED(kg_env%tnadd_mat)) THEN
         CALL dbcsr_print(kg_env%tnadd_mat(1)%matrix)
         END IF

      END IF

!      ! kinetic matrix trace
!      CALL calculate_ptrace(kinetic, rho_ao, ekinetic, nspins)
!      CALL calculate_ptrace(matrix_s, rho_ao, trace_ovrl, nspins)
!      ! core matrix trace
!      CALL calculate_ptrace(matrix_h, rho_ao, ecore, nspins)
!      ! band structure energy = Tr[KS*P]
!      CALL calculate_ptrace(matrix_ks, rho_ao, trace_ks_p, nspins)
!      
!      CALL get_qs_env(qs_env,para_env=para_env)
!      IF(para_env%mepos==0) THEN
!         WRITE(*,*) "trace[T*p]", ekinetic
!         WRITE(*,*) "trace[S*p]", trace_ovrl
!         WRITE(*,*) "trace[h*p]", ecore
!         WRITE(*,*) "trace[ks*p]", trace_ks_p
!      END IF
 
     !--------------------------------------------------------------------------
     
         IF (dft_control%qs_control%do_kg_ec) THEN
            ! Check for energy correction
            CALL kg_energy_correction(qs_env, calculate_forces=my_calc_forces)
         END IF

      !--------------------------------------------------------------------------
 
      IF(para_env%mepos==0) THEN
!        WRITE(*,*) "qs_energies - after Harris / before qs_ks_update_qs_env" 
!        
!        WRITE(*,*) "core        "  , energy%core
!        WRITE(*,*) "core_overlap " , energy%core_overlap
!        WRITE(*,*) "core_self   "  , energy%core_self
!        WRITE(*,*) "exc         "  , energy%exc
!        WRITE(*,*) "hartree     "  , energy%hartree
!        WRITE(*,*) "e_hartree,  "  , energy%e_hartree
!        WRITE(*,*) "total       "  , energy%total
!        WRITE(*,*) "tot_old     "  , energy%tot_old
!        WRITE(*,*) "kinetic     "  , energy%kinetic
     END IF
     !--------------------------------------------------------------------------


         IF (PRESENT(consistent_energies)) THEN
            ! consistent_energies would update core and KS matrix 
            ! to update energies to qs_forces
            ! since Harris forces are calculated differently we skip this
            !IF (consistent_energies .AND. .NOT. dft_control%qs_control%do_kg_ec ) THEN
            IF (consistent_energies) THEN
!               CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., just_energy=.TRUE.)

               ! add MP2 energy if necessary
               IF (ASSOCIATED(qs_env%mp2_env)) THEN
                  energy%total = energy%total+energy%mp2
               END IF
            END IF
         END IF

      !--------------------------------------------------------------------------
 
      IF(para_env%mepos==0) THEN
!        WRITE(*,*) "qs_energies - after qs_ks_update_qs_env" 
!        
!        WRITE(*,*) "core        "  , energy%core
!        WRITE(*,*) "core_overlap " , energy%core_overlap
!        WRITE(*,*) "core_self   "  , energy%core_self
!        WRITE(*,*) "exc         "  , energy%exc
!        WRITE(*,*) "hartree     "  , energy%hartree
!        WRITE(*,*) "e_hartree,  "  , energy%e_hartree
!        WRITE(*,*) "total       "  , energy%total
!        WRITE(*,*) "tot_old     "  , energy%tot_old
!        WRITE(*,*) "kinetic     "  , energy%kinetic
      END IF
     !--------------------------------------------------------------------------

         CALL qs_energies_properties(qs_env)

         IF (dft_control%qs_control%lrigpw) THEN
            CALL lri_print_stat(qs_env)
         END IF

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_energies

END MODULE qs_energy
