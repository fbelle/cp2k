!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief test and write routines for Harris energy correction on top of Kim-Gordon subsystem DFT
!>        
!> \author Fabian Belleflamme
!> \par History
!>      - 04.2018 init
! **************************************************************************************************
MODULE ec_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type

   
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_methods'

   PUBLIC :: ec_write_forces 

CONTAINS

! **************************************************************************************************
!> \brief  Write ec forces 
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


END MODULE ec_methods
