!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for an Harris functional energy correction
!> \par History
!>       2019.09 created [Fabian Belleflamme] - Refactoring out of KG 
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_environment_types
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_type
   USE integration_grid_types,          ONLY: deallocate_intgrid,&
                                              integration_grid_type
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_density_release,&
                                              lri_density_type,&
                                              lri_env_release,&
                                              lri_environment_type
   USE molecule_types,                  ONLY: molecule_type
   USE qs_dispersion_types,             ONLY: qs_dispersion_release,&
                                              qs_dispersion_type
   USE qs_grid_atom,                    ONLY: atom_integration_grid_type,&
                                              deallocate_atom_int_grid
   USE qs_neighbor_list_types,          ONLY: deallocate_neighbor_list_set,&
                                              neighbor_list_set_p_type
   USE task_list_types,                 ONLY: deallocate_task_list,&
                                              task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_environment_types'

   PUBLIC :: ec_environment_type, ec_env_release

   TYPE subset_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(task_list_type), POINTER :: task_list
   END TYPE subset_type

! *****************************************************************************
!> \brief Contains information on the energy correction functional
!> \par History
!>       03.2014 created
!>       09.2019 refactoring fbelle 
!> \author JGH
! *****************************************************************************
   TYPE ec_environment_type
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
      ! KS: Kohn-Sham; H: Core; S: overlap; T: kinetic energy;
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_ks
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_h
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_t
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER  :: matrix_p
      ! reduce basis
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER     :: mao_coef
   END TYPE ec_environment_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param ec_env The environment containg the energy correction parameters 
! **************************************************************************************************
   SUBROUTINE ec_env_release(ec_env)
      TYPE(ec_environment_type), POINTER                 :: ec_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'ec_env_release', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iab

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
      IF (ASSOCIATED(ec_env%task_list)) THEN
         CALL deallocate_task_list(ec_env%task_list)
      END IF
      ! reduced basis
      IF (ASSOCIATED(ec_env%mao_coef)) CALL dbcsr_deallocate_matrix_set(ec_env%mao_coef)
      ! dispersion environment
      IF (ASSOCIATED(ec_env%dispersion_env)) THEN
         CALL qs_dispersion_release(ec_env%dispersion_env)
      END IF

      DEALLOCATE (ec_env)

      CALL timestop(handle)

   END SUBROUTINE ec_env_release

END MODULE ec_environment_types
