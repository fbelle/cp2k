!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines used for Harris functional  
!>        Kohn-Sham calculation
!> \par History
!>       10.2020 created
!> \author Fabian Belleflamme
! **************************************************************************************************
MODULE ec_methods
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_row_template
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_init_random,&
                                              cp_fm_get_info,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_sum
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_matrix_pools,                 ONLY: mpools_create,&
                                              mpools_rebuild_fm_pools,&
                                              mpools_release,&
                                              qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_type,&
                                              mo_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ec_methods'

   PUBLIC :: ec_mos_init 

CONTAINS

! **************************************************************************************************
!> \brief init mos environment for energy correction
!>
!> \param qs_env ...
!> \param matrix_s Used as template
!> \param
!>
!> \par History
!>       2020.10 created [Fabian Belleflamme]
!> \author Fabian Belleflamme
! **************************************************************************************************
   SUBROUTINE ec_mos_init(qs_env, matrix_s)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type)                                   :: matrix_s

      CHARACTER(len=*), PARAMETER :: routineN = 'ec_mos_init', &
         routineP = moduleN//':'//routineN

      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      INTEGER                                            :: handle, ispin, nspins, multiplicity, nelectron
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER                                            :: n_ao
      INTEGER, DIMENSION(2)                              :: n_mo, nelectron_spin
      REAL(dp)                                           :: maxocc
      !
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
      INTEGER                                            :: nmo
      TYPE(qs_matrix_pools_type), POINTER                :: my_mpools
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos

      CALL timeset(routineN, handle)

      NULLIFY (ao_mo_fm_pools, dft_control, mos)
      NULLIFY (mo_coeff)

      CALL get_qs_env(qs_env = qs_env, &
                      dft_control = dft_control, &
                      blacs_env=blacs_env, &
                      qs_kind_set = qs_kind_set, &
                      nelectron_spin = nelectron_spin)
      nspins = dft_control%nspins

      !fbelle
      CALL get_qs_env(qs_env,para_env=para_env)

      ! Start setup
      CALL get_qs_kind_set(qs_kind_set, nsgf = n_ao, nelectron=nelectron)

      ! the total number of electrons
      nelectron = nelectron-dft_control%charge
      multiplicity = dft_control%multiplicity

     ! setting maxocc and n_mo
      IF (dft_control%nspins == 1) THEN
         maxocc = 2.0_dp
         nelectron_spin(1) = nelectron
         nelectron_spin(2) = 0
         IF (MODULO(nelectron, 2) == 0) THEN
            n_mo(1) = nelectron/2
         ELSE
            n_mo(1) = INT(nelectron/2._dp)+1
         END IF
         n_mo(2) = 0
      ELSE
         maxocc = 1.0_dp

         ! The simplist spin distribution is written here. Special cases will
         ! need additional user input
         IF (MODULO(nelectron+multiplicity-1, 2) /= 0) THEN
            CPABORT("LSD: try to use a different multiplicity")
         END IF

         nelectron_spin(1) = (nelectron+multiplicity-1)/2
         nelectron_spin(2) = (nelectron-multiplicity+1)/2

         IF (nelectron_spin(2) < 0) THEN
            CPABORT("LSD: too few electrons for this multiplicity")
         END IF

         n_mo(1) = nelectron_spin(1)
         n_mo(2) = nelectron_spin(2)

      END IF

      ! Allocate MO set
      ALLOCATE (mos(nspins))
      DO ispin = 1, nspins
         NULLIFY (mos(ispin)%mo_set)
         CALL allocate_mo_set(mo_set = mos(ispin)%mo_set, &
                              nao = n_ao, &
                              nmo = n_mo(ispin), &
                              nelectron = nelectron_spin(ispin), &
                              n_el_f = REAL(nelectron_spin(ispin), dp), &
                              maxocc = maxocc, &
                              flexible_electron_count = dft_control%relax_multiplicity)
      END DO

      CALL set_qs_env(qs_env, mos=mos)

      ! init pools
      NULLIFY (my_mpools)
      CALL mpools_create(mpools=my_mpools)
      CALL mpools_rebuild_fm_pools(mpools=my_mpools, &
                                   mos=mos, &
                                   blacs_env=blacs_env, &
                                   para_env=para_env)

!   *** finish initialization of the MOs ***

      CPASSERT(ASSOCIATED(mos))
      NULLIFY (mo_coeff, mo_coeff_b)
      DO ispin = 1, SIZE(mos)
         CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b)
         IF (.NOT. ASSOCIATED(mo_coeff)) THEN
            CALL init_mo_set(mos(ispin)%mo_set, &
                             fm_pool=my_mpools%ao_mo_fm_pools(ispin)%pool, &
                             name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))// &
                             "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))))
         END IF

         IF (.NOT. ASSOCIATED(mo_coeff_b)) THEN
            CALL cp_fm_get_info(mos(ispin)%mo_set%mo_coeff, ncol_global=nmo)
            CALL dbcsr_init_p(mos(ispin)%mo_set%mo_coeff_b)
            CALL cp_dbcsr_m_by_n_from_row_template(mos(ispin)%mo_set%mo_coeff_b, &
                                                   template=matrix_s, &
                                                   n = nmo, &
                                                   sym = dbcsr_type_no_symmetry)
         ENDIF
      END DO

      CALL mpools_release(mpools=my_mpools)

      CALL timestop(handle)

   END SUBROUTINE ec_mos_init      

END MODULE ec_methods
