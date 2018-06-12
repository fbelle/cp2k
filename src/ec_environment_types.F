!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for a for a Kim-Gordon-like partitioning into molecular
!>        subunits
!> \par History
!>       2014 created [JGH]
!>       05.2018 moved energy correcction types here
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_environment_types
   USE dbcsr_api,                       ONLY: dbcsr_deallocate_matrix_set,&
                                              dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE qs_dispersion_types,             ONLY: qs_dispersion_release,&
                                              qs_dispersion_type
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE task_list_types,                 ONLY: deallocate_task_list,&
                                              task_list_type
   !fbelle 
   USE qs_rho_types,                    ONLY: qs_rho_create,&
                                              qs_rho_get,&
                                              qs_rho_release,&
                                              qs_rho_type
   USE qs_force_types,                  ONLY: deallocate_qs_force,&
                                              qs_force_type
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              mo_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kg_environment_types'

   PUBLIC :: ec_env_release, energy_correction_type

   TYPE subset_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(task_list_type), POINTER :: task_list
   END TYPE subset_type

! *****************************************************************************
!> \brief Contains information on the energy correction functional for KG
!> \par History
!>       03.2014 created
!> \author JGH
! *****************************************************************************
   TYPE energy_correction_type
      CHARACTER(len=20)                                :: ec_name
      INTEGER                                          :: energy_functional
      INTEGER                                          :: ks_solver
      INTEGER                                          :: factorization
      REAL(KIND=dp)                                    :: eps_default
      ! basis set
      CHARACTER(len=20)                                :: basis
      LOGICAL                                          :: mao
      INTEGER                                          :: mao_max_iter
      REAL(KIND=dp)                                    :: mao_eps_grad
      ! energy components
      REAL(KIND=dp)                                    :: etotal
      REAL(KIND=dp)                                    :: eband, exc, ehartree, vhxc
      REAL(KIND=dp)                                    :: edispersion
      
      ! full neighbor lists and corresponding task list
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(task_list_type), POINTER                    :: task_list
      ! the XC function to be used for the correction, dispersion info
      TYPE(section_vals_type), POINTER                 :: xc_section
      TYPE(qs_dispersion_type), POINTER                :: dispersion_env
      ! matrices in complete basis
      ! KS: Kohn-Sham; H: Core; S: overlap; T: kinetic energy, W: wheighed density;
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_ks
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_h
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_t
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_w
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_delta
      
      ! reduce basis
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER     :: mao_coef
      
      ! store old rho
      TYPE(qs_rho_type), POINTER                    :: old_rho

      ! store force environment 
      TYPE(qs_force_type), DIMENSION(:), POINTER    :: force

      ! store MOS 
      TYPE(mo_set_p_type), DIMENSION(:),POINTER     :: mos

   END TYPE energy_correction_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ec_env ...
! **************************************************************************************************
   SUBROUTINE ec_env_release(ec_env)
      TYPE(energy_correction_type), POINTER                 :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_env_release', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iab

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(ec_env))

      ! neighbor lists
      IF (ASSOCIATED(ec_env%sab_orb)) THEN
         DO iab = 1, SIZE(ec_env%sab_orb)
            CALL deallocate_neighbor_list_set(ec_env%sab_orb(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (ec_env%sab_orb)
      END IF
      IF (ASSOCIATED(ec_env%sac_ppl)) THEN
         DO iab = 1, SIZE(ec_env%sac_ppl)
            CALL deallocate_neighbor_list_set(ec_env%sac_ppl(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (ec_env%sac_ppl)
      END IF
      IF (ASSOCIATED(ec_env%sap_ppnl)) THEN
         DO iab = 1, SIZE(ec_env%sap_ppnl)
            CALL deallocate_neighbor_list_set(ec_env%sap_ppnl(iab)%neighbor_list_set)
         END DO
         DEALLOCATE (ec_env%sap_ppnl)
      END IF
      ! operator matrices
      IF (ASSOCIATED(ec_env%matrix_ks)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_ks)
      IF (ASSOCIATED(ec_env%matrix_h)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_h)
      IF (ASSOCIATED(ec_env%matrix_s)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_s)
      IF (ASSOCIATED(ec_env%matrix_t)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_t)
      IF (ASSOCIATED(ec_env%matrix_p)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_p)
      IF (ASSOCIATED(ec_env%matrix_w)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_w)
      IF (ASSOCIATED(ec_env%matrix_delta)) CALL dbcsr_deallocate_matrix_set(ec_env%matrix_delta)
      IF (ASSOCIATED(ec_env%task_list)) THEN
         CALL deallocate_task_list(ec_env%task_list)
      END IF
      ! reduced basis
      IF (ASSOCIATED(ec_env%mao_coef)) CALL dbcsr_deallocate_matrix_set(ec_env%mao_coef)
      ! dispersion environment
      IF (ASSOCIATED(ec_env%dispersion_env)) THEN
         CALL qs_dispersion_release(ec_env%dispersion_env)
      END IF
      IF (ASSOCIATED(ec_env%old_rho)) THEN
         CALL qs_rho_release(ec_env%old_rho)
      END IF

      ! force environment 
      IF (ASSOCIATED(ec_env%force)) CALL deallocate_qs_force(ec_env%force)

      ! mos
      IF (ASSOCIATED(ec_env%mos)) THEN
         DO i = 1, SIZE(ec_env%mos)
            CALL deallocate_mo_set(ec_env%mos(i)%mo_set)
         END DO
         DEALLOCATE (ec_env%mos)
      END IF

      DEALLOCATE (ec_env)

      CALL timestop(handle)

   END SUBROUTINE ec_env_release

END MODULE ec_environment_types
