!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculates the moment integrals <a|r^m|b> and <a|r x d/dr|b>
!> \par History
!>      added angular moments (JGH 11.2012)
!> \author JGH (20.07.2006)
! **************************************************************************************************
MODULE qs_moments
   USE ai_angmom,                       ONLY: angmom
   USE ai_moments,                      ONLY: contract_cossin,&
                                              cossin,&
                                              cossin2,&
                                              moment
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE block_p_types,                   ONLY: block_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_lu_decompose
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_double,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_result_methods,               ONLY: cp_results_erase,&
                                              put_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_distribution_type, &
        dbcsr_get_block_p, dbcsr_p_type, dbcsr_set, dbcsr_trace, dbcsr_type, dbcsr_type_symmetric, &
        dbcsr_print, dbcsr_scale
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_sum
   USE moments_utils,                   ONLY: get_reference_point
   USE orbital_pointers,                ONLY: current_maxl,&
                                              indco,&
                                              ncoset
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: bohr,&
                                              debye
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE rt_propagation_types,            ONLY: get_rtp,&
                                              rt_prop_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_moments'

   ! Public subroutines
   PUBLIC :: build_berry_moment_matrix, build_local_moment_matrix
   PUBLIC :: build_berry_moment_matrix2
   PUBLIC :: build_berry_moment_matrix3
   PUBLIC :: debug_berry_moment_matrix
   PUBLIC :: build_berry_kpoint_matrix
   PUBLIC :: qs_moment_berry_phase, qs_moment_locop
   PUBLIC :: qs_moment_berry_phase2

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param moments ...
!> \param nmoments ...
!> \param ref_point ...
!> \param ref_points ...
! **************************************************************************************************
   SUBROUTINE build_local_moment_matrix(qs_env, moments, nmoments, ref_point, ref_points)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: moments
      INTEGER, INTENT(IN)                                :: nmoments
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: ref_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: ref_points

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_local_moment_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, last_jatom, &
         maxco, maxsgf, natom, ncoa, ncob, nkind, nm, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, npgfa, npgfb, &
                                                            nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: dab, rab2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: mab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rac, rb, rbc, rc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: mint
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      IF (nmoments < 1) RETURN

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, cell, particle_set, sab_orb)

      nm = (6+11*nmoments+6*nmoments**2+nmoments**3)/6-1
      CPASSERT(SIZE(moments) >= nm)

      NULLIFY (qs_kind_set, particle_set, sab_orb, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ! Allocate work storage
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, maxsgf=maxsgf)

      ALLOCATE (mab(maxco, maxco, nm))
      mab(:, :, :) = 0.0_dp

      ALLOCATE (work(maxco, maxsgf))
      work(:, :) = 0.0_dp

      ALLOCATE (mint(nm))
      DO i = 1, nm
         NULLIFY (mint(i)%block)
      END DO

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         IF (inode == 1) last_jatom = 0

         ! this guarentees minimum image convention
         ! anything else would not make sense
         IF (jatom == last_jatom) THEN
            CYCLE
         END IF

         last_jatom = jatom

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         DO i = 1, nm
            NULLIFY (mint(i)%block)
            CALL dbcsr_get_block_p(matrix=moments(i)%matrix, &
                                   row=irow, col=icol, BLOCK=mint(i)%block, found=found)
            mint(i)%block = 0._dp
         END DO

         ! fold atomic position back into unit cell
         IF (PRESENT(ref_points)) THEN
            rc(:) = 0.5_dp*(ref_points(:, iatom)+ref_points(:, jatom))
         ELSE IF (PRESENT(ref_point)) THEN
            rc(:) = ref_point(:)
         ELSE
            rc(:) = 0._dp
         END IF
         ! using PBC here might screw a molecule that fits the box (but e.g. hasn't been shifted by center_molecule)
         ! by folding around the center, such screwing can be avoided for a proper choice of center.
         ra(:) = pbc(particle_set(iatom)%r(:)-rc, cell)+rc
         rb(:) = pbc(particle_set(jatom)%r(:)-rc, cell)+rc
         ! we dont use PBC at this point
         rab(:) = ra(:)-rb(:)
         rac(:) = ra(:)-rc(:)
         rbc(:) = rb(:)-rc(:)
         rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
         dab = SQRT(rab2)

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               ! Calculate the primitive integrals
               CALL moment(la_max(iset), npgfa(iset), zeta(:, iset), &
                           rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), zetb(:, jset), &
                           rpgfb(:, jset), nmoments, rac, rbc, mab)

               ! Contraction step
               DO i = 1, nm

                  CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                             1.0_dp, mab(1, 1, i), SIZE(mab, 1), &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, work(1, 1), SIZE(work, 1))

                  IF (iatom <= jatom) THEN

                     CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                                1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                work(1, 1), SIZE(work, 1), &
                                1.0_dp, mint(i)%block(sgfa, sgfb), &
                                SIZE(mint(i)%block, 1))

                  ELSE

                     CALL dgemm("T", "N", nsgfb(jset), nsgfa(iset), ncoa, &
                                1.0_dp, work(1, 1), SIZE(work, 1), &
                                sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                1.0_dp, mint(i)%block(sgfb, sgfa), &
                                SIZE(mint(i)%block, 1))

                  END IF

               END DO

            END DO
         END DO

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! Release work storage
      DEALLOCATE (mab, basis_set_list)
      DEALLOCATE (work)
      DO i = 1, nm
         NULLIFY (mint(i)%block)
      END DO
      DEALLOCATE (mint)

      CALL timestop(handle)

   END SUBROUTINE build_local_moment_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param magmom ...
!> \param nmoments ...
!> \param ref_point ...
!> \param ref_points ...
! **************************************************************************************************
   SUBROUTINE build_local_magmom_matrix(qs_env, magmom, nmoments, ref_point, ref_points)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magmom
      INTEGER, INTENT(IN)                                :: nmoments
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: ref_point
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: ref_points

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_local_magmom_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, last_jatom, &
         maxco, maxsgf, natom, ncoa, ncob, nkind, nm, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, npgfa, npgfb, &
                                                            nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: dab, rab2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: mab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rac, rb, rbc, rc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(block_p_type), ALLOCATABLE, DIMENSION(:)      :: mint
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      IF (nmoments < 1) RETURN

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, cell, particle_set, sab_orb)
      NULLIFY (matrix_s)

      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)

      ! magnetic dipoles/angular moments only
      nm = 3
      CALL dbcsr_allocate_matrix_set(magmom, nm)

      DO i = 1, nm
         ALLOCATE (magmom(i)%matrix)
         CALL dbcsr_copy(magmom(i)%matrix, matrix_s(1)%matrix, "Moments")
         CALL dbcsr_set(magmom(i)%matrix, 0.0_dp)
      END DO

      NULLIFY (qs_kind_set, particle_set, sab_orb, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ! Allocate work storage
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, maxsgf=maxsgf)

      ALLOCATE (mab(maxco, maxco, nm))
      mab(:, :, :) = 0.0_dp

      ALLOCATE (work(maxco, maxsgf))
      work(:, :) = 0.0_dp

      ALLOCATE (mint(nm))
      DO i = 1, nm
         NULLIFY (mint(i)%block)
      END DO

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         IF (inode == 1) last_jatom = 0

         ! this guarentees minimum image convention
         ! anything else would not make sense
         IF (jatom == last_jatom) THEN
            CYCLE
         END IF

         last_jatom = jatom

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         DO i = 1, nm
            NULLIFY (mint(i)%block)
            CALL dbcsr_get_block_p(matrix=magmom(i)%matrix, &
                                   row=irow, col=icol, BLOCK=mint(i)%block, found=found)
            mint(i)%block = 0._dp
         END DO

         ! fold atomic position back into unit cell
         IF (PRESENT(ref_points)) THEN
            rc(:) = 0.5_dp*(ref_points(:, iatom)+ref_points(:, jatom))
         ELSE IF (PRESENT(ref_point)) THEN
            rc(:) = ref_point(:)
         ELSE
            rc(:) = 0._dp
         END IF
         ! using PBC here might screw a molecule that fits the box (but e.g. hasn't been shifted by center_molecule)
         ! by folding around the center, such screwing can be avoided for a proper choice of center.
         ra(:) = pbc(particle_set(iatom)%r(:)-rc, cell)+rc
         rb(:) = pbc(particle_set(jatom)%r(:)-rc, cell)+rc
         ! we dont use PBC at this point
         rab(:) = ra(:)-rb(:)
         rac(:) = ra(:)-rc(:)
         rbc(:) = rb(:)-rc(:)
         rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
         dab = SQRT(rab2)

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               ! Calculate the primitive integrals
               CALL angmom(la_max(iset), npgfa(iset), zeta(:, iset), &
                           rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), zetb(:, jset), &
                           rpgfb(:, jset), rac, rbc, mab)

               ! Contraction step
               DO i = 1, nm

                  CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                             1.0_dp, mab(1, 1, i), SIZE(mab, 1), &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, work(1, 1), SIZE(work, 1))

                  IF (iatom <= jatom) THEN

                     CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                                1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                work(1, 1), SIZE(work, 1), &
                                1.0_dp, mint(i)%block(sgfa, sgfb), &
                                SIZE(mint(i)%block, 1))

                  ELSE

                     CALL dgemm("T", "N", nsgfb(jset), nsgfa(iset), ncoa, &
                                1.0_dp, work(1, 1), SIZE(work, 1), &
                                sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                1.0_dp, mint(i)%block(sgfb, sgfa), &
                                SIZE(mint(i)%block, 1))

                  END IF

               END DO

            END DO
         END DO

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      ! Release work storage
      DEALLOCATE (mab, basis_set_list)
      DEALLOCATE (work)
      DO i = 1, nm
         NULLIFY (mint(i)%block)
      END DO
      DEALLOCATE (mint)

      CALL timestop(handle)

   END SUBROUTINE build_local_magmom_matrix


! **************************************************************************************************
!> \brief Debugging of moments integral cossin, and its derivatives wrt k
!>        Check analytical integrals using the numerical integrals as refernce
!> \param qs_env ...
!> \param cosmat ...
!> \param sinmat ...
!> \param kvec ...
!> \param sab_orb_external ...
!> \param basis_type ...
!> \param nderivative ... Order of derivative wrt k
! **************************************************************************************************
   SUBROUTINE debug_berry_moment_matrix(qs_env, nmoments)
   
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nmoments

      CHARACTER(LEN=*), PARAMETER :: routineN = 'debug_berry_moment_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j, k, l, nmom, handle
      REAL(dp), DIMENSION(3)                             :: kvec
      REAL(dp)                                           :: std_value
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat

      ! fbelle
      TYPE(dbcsr_type), DIMENSION(:), POINTER            :: dkcos, dksin
      TYPE(dbcsr_type), DIMENSION(:), POINTER            :: d2kcos, d2ksin
      REAL(KIND=dp)                                      :: dx   

      TYPE(dbcsr_type), DIMENSION(:,:), POINTER          :: analyt_cosmat, analyt_sinmat
      TYPE(dbcsr_type), DIMENSION(:,:), POINTER          :: numer_cosmat, numer_sinmat
!      TYPE(dbcsr_type), DIMENSION(:), POINTER            :: numcos, numsin
      REAL (KIND=dp), DIMENSION(3,3)                     :: anactrace, anastrace
      REAL (KIND=dp), DIMENSION(3,3)                     :: numctrace, numstrace

      CALL timeset(routineN, handle)

      ! restrict maximum moment available
      nmom = MIN(nmoments, 2)

      NULLIFY (cell, matrix_s)

      CALL get_qs_env(qs_env, &
                      cell = cell, &
                      matrix_s = matrix_s)

      NULLIFY (cosmat, sinmat)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL dbcsr_set(cosmat, 0.0_dp)
      CALL dbcsr_set(sinmat, 0.0_dp)

      ! Init matrices to store analytically calculated matrices
      NULLIFY (analyt_cosmat, analyt_sinmat)
      ALLOCATE (analyt_cosmat(3,3), analyt_sinmat(3,3))
      ! Init matrices to store numerically calculated matrices
      NULLIFY (numer_cosmat, numer_sinmat)
      ALLOCATE (numer_cosmat(3,3), numer_sinmat(3,3))

      DO i = 1, 3
         DO j = 1, 3
            CALL dbcsr_copy(analyt_cosmat(i,j), matrix_s(1)%matrix, 'COS ANAL')
            CALL dbcsr_copy(analyt_sinmat(i,j), matrix_s(1)%matrix, 'SIN ANAL')
            CALL dbcsr_set(analyt_cosmat(i,j), 0.0_dp)
            CALL dbcsr_set(analyt_sinmat(i,j), 0.0_dp)

            CALL dbcsr_copy(numer_cosmat(i,j), matrix_s(1)%matrix, 'COS NUMER')
            CALL dbcsr_copy(numer_sinmat(i,j), matrix_s(1)%matrix, 'SIN NUMER')
            CALL dbcsr_set(numer_cosmat(i,j), 0.0_dp)
            CALL dbcsr_set(numer_sinmat(i,j), 0.0_dp)
         END DO
      END DO         

!      ! more matrices for numer storage
!      NULLIFY (numcos, numsin)
!      ALLOCATE (numcos(2), numsin(2))
!
!      DO i = 1, 2
!         CALL dbcsr_copy(numcos(i), matrix_s(1)%matrix, 'COS MOM')
!         CALL dbcsr_copy(numsin(i), matrix_s(1)%matrix, 'SIN MOM')
!         CALL dbcsr_set(numcos(i), 0.0_dp)
!         CALL dbcsr_set(numsin(i), 0.0_dp)
!      END DO         

      ! Initialize and allocate matrices to store derivatives
      IF(nmom >= 1) THEN
         NULLIFY(dkcos, dksin)
         ALLOCATE(dkcos(3), dksin(3))
         DO i = 1, 3
            CALL dbcsr_copy(dkcos(i), matrix_s(1)%matrix, 'DKCOS MOM')
            CALL dbcsr_copy(dksin(i), matrix_s(1)%matrix, 'DKSIN MOM')
            CALL dbcsr_set(dkcos(i), 0.0_dp)
            CALL dbcsr_set(dksin(i), 0.0_dp)
         END DO
      END IF
      IF(nmom >= 2) THEN
         NULLIFY(d2kcos, d2ksin)
         ALLOCATE(d2kcos(6), d2ksin(6))
         DO i = 1, 6
            CALL dbcsr_copy(d2kcos(i), matrix_s(1)%matrix, 'D2KCOS MOM')
            CALL dbcsr_copy(d2ksin(i), matrix_s(1)%matrix, 'D2KSIN MOM')
            CALL dbcsr_set(d2kcos(i), 0.0_dp)
            CALL dbcsr_set(d2ksin(i), 0.0_dp)
         END DO
      END IF

      ! Deviation
      dx = 1.0E-3

      !------------------------------------------------------------------------
      ! Analytic call of cossin
      !------------------------------------------------------------------------

      ! 3 distinct kvecs -> each gives 3 first derivs
      !                  -> each gives 6 second derivs

      ! First, evaluate analytic 
      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
            DO i = 1, 3

               kvec(:) = twopi*cell%h_inv(i, :)
               CALL build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec, &
                                               nderivative = 1, dkcos = dkcos, dksin = dksin)
       
               DO j = 1, 3
                  ! copy cosmat sinmat into separate array
                  CALL dbcsr_copy(analyt_cosmat(i, j), dkcos(j)) 
                  CALL dbcsr_copy(analyt_sinmat(i, j), dksin(j)) 

                  ! store only trace
                  CALL dbcsr_trace(analyt_cosmat(i, j), anactrace(i, j))
                  CALL dbcsr_trace(analyt_cosmat(i, j), anastrace(i, j))
               END DO
            END DO

         CASE (2)
            ! Quadrupole
            DO i = 1, 3
               DO j = i, 3
                  kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :))
                  CALL build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec, & 
                                                  nderivative = 2, dkcos = dkcos, dksin = dksin, &
                                                  d2kcos = d2kcos, d2ksin = d2ksin)
                  ! copy cosmat sinmat into separate array
!                  analyt_cosmat(i, j, 1:6) = d2kcos(1:6) 
!                  analyt_sinmat(i, j, 1:6) = d2ksin(1:6) 
               END DO
            END DO
         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO

      WRITE(*,*) "----------------------------------------------"
      WRITE(*,*) "ANALYTIC DERIVATIVE analyt_cosmat(1, 1)"
      WRITE(*,*) "----------------------------------------------"
      ! print analytic matrices
      CALL dbcsr_print(analyt_cosmat(1, 1))
      CALL dbcsr_print(analyt_cosmat(1, 2))
      CALL dbcsr_print(analyt_cosmat(1, 3))

      WRITE(*,*) "TRACE dcos:", anactrace(i,1:3)
      WRITE(*,*) "TRACE dsin:", anastrace(i,1:3)

      WRITE(*,*) "----------------------------------------------"
    
      !------------------------------------------------------------------------
      ! Numerical calculation 
      !------------------------------------------------------------------------

      ! kvec is a 3D object, needs to be perturbed in all 3 directions
      ! just like particle position in forces
      DO i = 1, 3

         ! do each kvec twice: (+) deviation and (-) deviation
         kvec(:) = twopi*cell%h_inv(i, :)

         ! loop over the 3 derivation directions 
         DO k = 1, 3
            ! save inital value of kvec in direction k
            std_value = kvec(k)

            DO j = 1, 2
               CALL dbcsr_set(cosmat, 0.0_dp)
               CALL dbcsr_set(sinmat, 0.0_dp)
               ! Apply deviation on kvecs direction k
               kvec(k) = std_value-(-1.0_dp)**j*dx

               ! Compute cosmat and sinmat
!               cosmat => numcos(j)
!               sinmat => numsin(j)
               CALL build_berry_moment_matrix3(qs_env = qs_env, &
                                               cosmat = cosmat, &
                                               sinmat = sinmat, &
                                               kvec = kvec)
               IF (j = 1) THEN
                  CALL dbcsr_copy(numer_cosmat(i, k), cosmat)
                  CALL dbcsr_copy(numer_sinmat(i, k), sinmat)
                  ! Numcos within j loop
               END IF
               IF (j = 2) THEN
                  CALL dbcsr_copy(numer_cosmat(i,k), cosmat)
                  CALL dbcsr_copy(numer_cosmat(i,k), sinmat)
               END IF

            END DO
            kvec(k) = std_value

            ! Numerical calculation: 
            ! analyt_cosmat ~~ numer_cosmat = (numcos(2) - numcos(1))/dx

!            CALL dbcsr_add(numer_cosmat(i, k), numcos(2), 1.0_dp, -1.0_dp)
!            CALL dbcsr_add(numer_sinmat(i, k), numsin(2), 1.0_dp, -1.0_dp)

            ! dividing by dx
            !CALL dbcsr_scale(numer_cosmat(i, k), 1000.0_dp)
            !CALL dbcsr_scale(numer_sinmat(i, k), 1000.0_dp)

            ! store only trace
!            CALL dbcsr_trace(numer_cosmat(i,k), numctrace(i,k))
!            CALL dbcsr_trace(numer_cosmat(i,k), numstrace(i,k))

         END DO

      END DO

      WRITE(*,*) "----------------------------------------------"
      WRITE(*,*) "NUMERIC DERIVATIVE numer_cosmat(1, 1)"
      WRITE(*,*) "----------------------------------------------"
      CALL dbcsr_print(numer_cosmat(1, 1))
      CALL dbcsr_print(numer_cosmat(1, 2))
      CALL dbcsr_print(numer_cosmat(1, 3))
!      DO k = 1, 3
!         WRITE(*,*) "TRACE dcos:", numctrace(i,1:3)
!         WRITE(*,*) "TRACE dsin:", numstrace(i,1:3)
!      END DO
      WRITE(*,*) "----------------------------------------------"






      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)


      DEALLOCATE(analyt_cosmat)
      DEALLOCATE(analyt_sinmat)
      DEALLOCATE(numer_cosmat)
      DEALLOCATE(numer_sinmat)



!      CALL dbcsr_deallocate_matrix(analyt_cosmat)
!      CALL dbcsr_deallocate_matrix(analyt_sinmat)
!      CALL dbcsr_deallocate_matrix(numer_cosmat)
!      CALL dbcsr_deallocate_matrix(numer_sinmat)
!      CALL dbcsr_deallocate_matrix(numcos)
!      CALL dbcsr_deallocate_matrix(numsin)


!      DO i = 1, 3
!         DO j = 1, 3
!            CALL dbcsr_deallocate_matrix(analyt_cosmat(i,j))
!            CALL dbcsr_deallocate_matrix(analyt_sinmat(i,j))
!           
!            CALL dbcsr_deallocate_matrix(numer_cosmat(i,j))
!            CALL dbcsr_deallocate_matrix(numer_sinmat(i,j))
!         END DO
!      END DO
!      DO i = 1, 2
!            CALL dbcsr_deallocate_matrix(numcos(i))
!            CALL dbcsr_deallocate_matrix(numsin(i))
!      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cosmat ...
!> \param sinmat ...
!> \param kvec ...
!> \param sab_orb_external ...
!> \param basis_type ...
!> \param nderivative ... Order of derivative wrt k
!> \param dkcos dksin ... First derivative wrt k of moment matrix
!> \param d2kcos d2ksin ... Second derivative wrt k of  moment matrix
! **************************************************************************************************
   SUBROUTINE build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec, &
                                         sab_orb_external, basis_type, nderivative, &
                                         dkcos, dksin, d2kcos, d2ksin)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kvec
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_external
      CHARACTER(len=*), OPTIONAL                         :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: nderivative
      TYPE(dbcsr_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: dkcos, dksin, d2kcos, d2ksin 

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_berry_moment_matrix3', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, ldab, ldsa, &
         ldsb, ldwork, natom, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: cblock, sblock, work
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
 
      ! fbelle
      INTEGER                                            :: maxder, nder, i
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: cosmint, sinmint
      REAL(dp), DIMENSION(:, :), POINTER                 :: ckblk1, ckblk2, ckblk3
      REAL(dp), DIMENSION(:, :), POINTER                 :: skblk1, skblk2, skblk3
      REAL(dp), DIMENSION(:, :), POINTER                 :: c2kblk1,c2kblk2,c2kblk3,c2kblk4,c2kblk5,c2kblk6
      REAL(dp), DIMENSION(:, :), POINTER                 :: s2kblk1,s2kblk2,s2kblk3,s2kblk4,s2kblk5,s2kblk6

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, particle_set, sab_orb, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      IF (PRESENT(sab_orb_external)) sab_orb => sab_orb_external

      CALL dbcsr_set(sinmat, 0.0_dp)
      CALL dbcsr_set(cosmat, 0.0_dp)

      IF (PRESENT(dkcos) .AND. PRESENT(dksin)) THEN
         DO i = 1, 3
            CALL dbcsr_set(dkcos(i), 0.0_dp)
            CALL dbcsr_set(dksin(i), 0.0_dp)
         END DO
      END IF
      
      IF (PRESENT(d2kcos) .AND. PRESENT(d2ksin)) THEN
         DO i = 1, 6
            CALL dbcsr_set(d2kcos(i), 0.0_dp)
            CALL dbcsr_set(d2ksin(i), 0.0_dp)
         END DO
      END IF

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=ldwork, basis_type=basis_type)
      ldab = ldwork

      ! derivative size
      IF (PRESENT(nderivative)) THEN
         nder = nderivative
      ELSE
         nder = 0
      END IF
      maxder = ncoset(nder)

      ! Allocate array for (a|eikr|b) intrgral and its derivatives
      ALLOCATE(cosmint(ldab, ldab, maxder))
      ALLOCATE(sinmint(ldab, ldab, maxder))
      ALLOCATE (work(ldwork, ldwork))

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a, basis_type=basis_type)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet
         
         ldsa = SIZE(sphi_a, 1)
         ldsb = SIZE(sphi_b, 1)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         NULLIFY (cblock)
         CALL dbcsr_get_block_p(matrix=cosmat, &
                                row=irow, col=icol, BLOCK=cblock, found=found)
         NULLIFY (sblock)
         CALL dbcsr_get_block_p(matrix=sinmat, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
   
         IF(PRESENT(dkcos) .AND. PRESENT(dksin)) THEN
!            WRITE(*,*) "dkcos - associate blocks"
           
            ! THIS IS AWEFUL - but 'BLOCK' wants a pointer
            NULLIFY (ckblk1,ckblk2,ckblk3)
            NULLIFY (skblk1,skblk2,skblk3)

               CALL dbcsr_get_block_p(matrix=dkcos(1), &
                                      row=irow, col=icol, BLOCK=ckblk1, found=found)
               CALL dbcsr_get_block_p(matrix=dkcos(2), &
                                      row=irow, col=icol, BLOCK=ckblk2, found=found)
               CALL dbcsr_get_block_p(matrix=dkcos(3), &
                                      row=irow, col=icol, BLOCK=ckblk3, found=found)
               CALL dbcsr_get_block_p(matrix=dksin(1), &
                                      row=irow, col=icol, BLOCK=skblk1, found=found)
               CALL dbcsr_get_block_p(matrix=dksin(2), &
                                      row=irow, col=icol, BLOCK=skblk2, found=found)
               CALL dbcsr_get_block_p(matrix=dksin(3), &
                                      row=irow, col=icol, BLOCK=skblk3, found=found)
         END IF         

         IF(PRESENT(d2kcos) .AND. PRESENT(d2ksin)) THEN
!            WRITE(*,*) "d2kcos - associate blocks"
      
            ! THIS IS AWEFUL - but 'BLOCK' wants a pointer
            NULLIFY (c2kblk1,c2kblk2,c2kblk3,c2kblk4,c2kblk5,c2kblk6)
            NULLIFY (s2kblk1,s2kblk2,s2kblk3,s2kblk4,s2kblk5,s2kblk6)

               CALL dbcsr_get_block_p(matrix=d2kcos(1), &
                                      row=irow, col=icol, BLOCK=c2kblk1, found=found)
               CALL dbcsr_get_block_p(matrix=d2kcos(2), &
                                      row=irow, col=icol, BLOCK=c2kblk2, found=found)
               CALL dbcsr_get_block_p(matrix=d2kcos(3), &
                                      row=irow, col=icol, BLOCK=c2kblk3, found=found)
               CALL dbcsr_get_block_p(matrix=d2kcos(4), &
                                      row=irow, col=icol, BLOCK=c2kblk4, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(5), &
                                      row=irow, col=icol, BLOCK=c2kblk5, found=found)
               CALL dbcsr_get_block_p(matrix=d2kcos(6), &
                                      row=irow, col=icol, BLOCK=c2kblk6, found=found)

               CALL dbcsr_get_block_p(matrix=d2ksin(1), &
                                      row=irow, col=icol, BLOCK=s2kblk1, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(2), &
                                      row=irow, col=icol, BLOCK=s2kblk2, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(3), &
                                      row=irow, col=icol, BLOCK=s2kblk3, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(4), &
                                      row=irow, col=icol, BLOCK=s2kblk4, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(5), &
                                      row=irow, col=icol, BLOCK=s2kblk5, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(6), &
                                      row=irow, col=icol, BLOCK=s2kblk6, found=found)
         END IF         

         IF (ASSOCIATED(cblock) .AND. .NOT. ASSOCIATED(sblock) .OR. &
             .NOT. ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN
            CPABORT("")
         ENDIF

         IF (ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN

            ra(:) = pbc(particle_set(iatom)%r(:), cell)
            rb(:) = ra+rab
            dab = SQRT(rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3))

            DO iset = 1, nseta

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               
               DO jset = 1, nsetb

                  IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)

                  ! Calculate the (a|exp(-ikr)|b) integrals and derivatives (a|(-ir)^nder*exp(-ikr)|b)
                  ! Careful, this derivative is wrt k, NOT the ion positions
                  SELECT CASE (nder)
                  CASE(0)
                    CALL cossin2(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1))
        
                  CASE(1)
!                    WRITE(*,*) "cossin2 - First Derivative"
                    CALL cossin2(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1), &
                                dksc = cosmint(:, :, 2:4), dkss = sinmint(:, :, 2:4))

                  CASE(2)
!                    WRITE(*,*) "cossin2 - Second Derivative"
                    CALL cossin2(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab=cosmint(:, :, 1), sinab=sinmint(:, :, 1), &
                                dksc=cosmint(:, :, 2:4), dkss=sinmint(:, :, 2:4), &
                                ddksc=cosmint(:, :, 5:10), ddkss=sinmint(:, :, 5:10))
        
                  CASE DEFAULT
                     CPABORT("Derivatives of (a|exp(-ikr)|b) only up to second order")
                  END SELECT

                 
                  ! Contraction from cartesian to spherical basis functions 
                  !WRITE(*,*) "Contract cosmint(:, :, 1)"
                  CALL contract_cossin(cblock, sblock, &
                                       iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                       jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                       cosmint(:, :, 1), sinmint(:, :, 1), ldab, work, ldwork)

                  IF (nder >= 1) THEN
                    ! Contraction of matrices storing first derivatives
                    !WRITE(*,*) "Contract cosmint(:, :, 2:4)"
                    CALL contract_cossin(ckblk1, skblk1, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 2), sinmint(:, :, 2), ldab, work, ldwork)
                    CALL contract_cossin(ckblk2, skblk2, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 3), sinmint(:, :, 3), ldab, work, ldwork)
                    CALL contract_cossin(ckblk3, skblk3, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 4), sinmint(:, :, 4), ldab, work, ldwork)
                  END IF

                  IF (nder >= 2) THEN
                    ! Contraction of matrices storing second derivatives
                    !WRITE(*,*) "Contract cosmint(:, :, 5:10)"
                    CALL contract_cossin(c2kblk1, s2kblk1, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 5), sinmint(:, :, 5), ldab, work, ldwork)
                    CALL contract_cossin(c2kblk2, s2kblk2, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 6), sinmint(:, :, 6), ldab, work, ldwork)
                    CALL contract_cossin(c2kblk3, s2kblk3, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 7), sinmint(:, :, 7), ldab, work, ldwork)
                    CALL contract_cossin(c2kblk4, s2kblk4, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 8), sinmint(:, :, 8), ldab, work, ldwork)
                    CALL contract_cossin(c2kblk5, s2kblk5, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 9), sinmint(:, :, 9), ldab, work, ldwork)
                    CALL contract_cossin(c2kblk6, s2kblk6, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 10), sinmint(:, :, 10), ldab, work, ldwork)
                  END IF

               END DO
            END DO

         ENDIF

      END DO

!      WRITE(*,*) ""
!      WRITE(*,*) "Integral:"
!      CALL dbcsr_print(cosmat)
!      CALL dbcsr_print(sinmat)
!      CALL dbcsr_print(dkcos(1))
!      CALL dbcsr_print(dkcos(2))
!      CALL dbcsr_print(dkcos(3))

      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (cosmint)
      DEALLOCATE (sinmint)
      DEALLOCATE (work)
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

!      WRITE(*,*) "build_berry_moment_matrix3 - END"

   END SUBROUTINE build_berry_moment_matrix3


! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cosmat ...
!> \param sinmat ...
!> \param kvec ...
!> \param sab_orb_external ...
!> \param basis_type ...
!> \param nderivative ... Order of derivative wrt k
! **************************************************************************************************
   SUBROUTINE build_berry_moment_matrix2(qs_env, cosmat, sinmat, kvec, &
                                         sab_orb_external, basis_type, nderivative, &
                                         dkcos, dksin, dk2cos, dk2sin )

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kvec
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_external
      CHARACTER(len=*), OPTIONAL                         :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: nderivative
      TYPE(dbcsr_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: dkcos, dksin, dk2cos, dk2sin 


      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_berry_moment_matrix2', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, ldab, ldsa, &
         ldsb, ldwork, natom, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: cblock, sblock, work
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
 
      ! fbelle
      INTEGER                                            :: maxder, nder, i, j
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: cosmint, sinmint
      REAL(dp), DIMENSION(:, :, :), POINTER              :: ckblock, skblock
      REAL(dp), DIMENSION(:, :), POINTER                 :: cos_block, sin_block

      WRITE(*,*) "build_berry_moment_matrix2"

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, particle_set, sab_orb, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      IF (PRESENT(sab_orb_external)) sab_orb => sab_orb_external

      CALL dbcsr_set(sinmat, 0.0_dp)
      CALL dbcsr_set(cosmat, 0.0_dp)

      IF (PRESENT(dkcos) .AND. PRESENT(dksin)) THEN

         WRITE(*,*) "SIZE(dkcos, 1)", SIZE(dkcos, 1)

         DO i = 1, 3
            CALL dbcsr_set(dkcos(i), 0.0_dp)
            CALL dbcsr_set(dksin(i), 0.0_dp)
         END DO
      END IF
      
      IF (PRESENT(dk2cos) .AND. PRESENT(dk2sin)) THEN
         DO i = 1, 6
            CALL dbcsr_set(dk2cos(i), 0.0_dp)
            CALL dbcsr_set(dk2sin(i), 0.0_dp)
         END DO
      END IF

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=ldwork, basis_type=basis_type)
      ldab = ldwork

      ! derivative size
      IF (PRESENT(nderivative)) THEN
         nder = nderivative
      ELSE
         nder = 0
      END IF
      maxder = ncoset(nder)

      ! Allocate array for (a|eikr|b) intrgral and its derivatives
      ALLOCATE(cosmint(ldab, ldab, maxder))
      ALLOCATE(sinmint(ldab, ldab, maxder))
      ALLOCATE (work(ldwork, ldwork))

      WRITE(*,*) "ldab", ldab

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a, basis_type=basis_type)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

!         WRITE(*,*) "first_sgfa",  basis_set_a%first_sgf
!         WRITE(*,*) "la_max",  basis_set_a%lmax
!         WRITE(*,*) "la_min",  basis_set_a%lmin
!         WRITE(*,*) "npgfa",  basis_set_a%npgf
!         WRITE(*,*) "nseta",  basis_set_a%nset
!         WRITE(*,*) "nsgfa", basis_set_a%nsgf_set
!         WRITE(*,*) "rpgfa", basis_set_a%pgf_radius
!         WRITE(*,*) "set_radius_a",  basis_set_a%set_radius
!         WRITE(*,*) "zeta", basis_set_a%zet
!         WRITE(*,*) "first_sgfb",  basis_set_b%first_sgf
!         WRITE(*,*) "lb_max",  basis_set_b%lmax
!         WRITE(*,*) "lb_min",  basis_set_b%lmin
!         WRITE(*,*) "npgfb",  basis_set_b%npgf
!         WRITE(*,*) "nsetb",  basis_set_b%nset
!         WRITE(*,*) "nsgfb", basis_set_b%nsgf_set
!         WRITE(*,*) "rpgfb", basis_set_b%pgf_radius
!         WRITE(*,*) "set_radius_b",  basis_set_b%set_radius
!         WRITE(*,*) "zetb", basis_set_b%zet
!         WRITE(*,*) ""


         ldsa = SIZE(sphi_a, 1)
         ldsb = SIZE(sphi_b, 1)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         NULLIFY (cblock)
         CALL dbcsr_get_block_p(matrix=cosmat, &
                                row=irow, col=icol, BLOCK=cblock, found=found)
         NULLIFY (sblock)
         CALL dbcsr_get_block_p(matrix=sinmat, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
   
         IF(ASSOCIATED(dkcos) .AND. ASSOCIATED(dksin)) THEN
            NULLIFY (ckblock)
            NULLIFY (skblock)
            DO i = 1, 3
!               cos_block => ckblock(:, :, i)
!               sin_block => skblock(:, :, i)
               CALL dbcsr_get_block_p(matrix=dkcos(i), &
                                      row=irow, col=icol, BLOCK=cos_block, found=found)
               CALL dbcsr_get_block_p(matrix=dksin(i), &
                                      row=irow, col=icol, BLOCK=sin_block, found=found)
            END DO
         END IF         
                       
                       
                       
                       
         IF (ASSOCIATED(cblock) .AND. .NOT. ASSOCIATED(sblock) .OR. &
             .NOT. ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN
            CPABORT("")
         ENDIF

         IF (ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN

            ra(:) = pbc(particle_set(iatom)%r(:), cell)
            rb(:) = ra+rab
            dab = SQRT(rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3))

            DO iset = 1, nseta

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               
               DO jset = 1, nsetb

                  IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)

                  ! Calculate the (a|exp(-ikr)|b) integrals and derivatives (a|(-ir)^nder*exp(-ikr)|b)
                  ! Careful, this derivative is wrt k, NOT the ion positions
                  SELECT CASE (nder)
                  CASE(0)
                    CALL cossin2(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1))

                    ! Contraction from cartesian to spherical basis functions 
                    CALL contract_cossin(cblock, sblock, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 1), sinmint(:, :, 1), ldab, work, ldwork)
        
                  CASE(1)
                    WRITE(*,*) "cossin2 - First Derivative"
                    CALL cossin2(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1), &
                                dksc = cosmint(:, :, 2:4), dkss = sinmint(:, :, 2:4))

                    ! Contraction from cartesian to spherical basis functions 
                    CALL contract_cossin(cblock, sblock, &
                                         iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                         jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                         cosmint(:, :, 1), sinmint(:, :, 1), ldab, work, ldwork)
                    ! Contraction of matrices storing first derivatives
                    DO i = 1, 3
                       j = i + 1
                       CALL dbcsr_get_block_p(matrix=dkcos(i), &
                                              row=irow, col=icol, BLOCK=cos_block, found=found)
                       CALL dbcsr_get_block_p(matrix=dksin(i), &
                                              row=irow, col=icol, BLOCK=sin_block, found=found)
                       
                       CALL contract_cossin(cos_block, sin_block, &
                                            iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                            jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                            cosmint(:, :, j), sinmint(:, :, j), ldab, work, ldwork)
                    END DO

                  CASE(2)
                    WRITE(*,*) "cossin2 - Second Derivative"
                    CALL cossin2(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab=cosmint(:, :, 1), sinab=sinmint(:, :, 1), &
                                dksc=cosmint(:, :, 2:4), dkss=sinmint(:, :, 2:4), &
                                ddksc=cosmint(:, :, 5:10), ddkss=sinmint(:, :, 5:10))
                  CASE DEFAULT
                     CPABORT("Derivatives of (a|exp(-ikr)|b) only up to second order")
                  END SELECT

!                  CALL contract_cossin(cblock, sblock, &
!                                       iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
!                                       jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
!                                       cosmint(:, :, 1), sinmint(:, :, 1), ldab, work, ldwork)
               END DO
            END DO

         ENDIF

      END DO

      WRITE(*,*) ""
      WRITE(*,*) "Integral:"
      CALL dbcsr_print(cosmat)
      CALL dbcsr_print(sinmat)

      CALL dbcsr_print(dkcos(1))
      CALL dbcsr_print(dkcos(2))
      CALL dbcsr_print(dkcos(3))



      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (cosmint)
      DEALLOCATE (sinmint)
      DEALLOCATE (work)
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

      WRITE(*,*) "build_berry_moment_matrix2 - END"

   END SUBROUTINE build_berry_moment_matrix2


! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cosmat ...
!> \param sinmat ...
!> \param kvec ...
!> \param sab_orb_external ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec, sab_orb_external, basis_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kvec
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_external
      CHARACTER(len=*), OPTIONAL                         :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_berry_moment_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, iatom, icol, ikind, inode, irow, iset, jatom, jkind, jset, ldab, ldsa, &
         ldsb, ldwork, natom, ncoa, ncob, nkind, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: cblock, cosab, sblock, sinab, work
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(cell_type), POINTER                           :: cell
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, particle_set, sab_orb, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      IF (PRESENT(sab_orb_external)) sab_orb => sab_orb_external

      CALL dbcsr_set(sinmat, 0.0_dp)
      CALL dbcsr_set(cosmat, 0.0_dp)

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=ldwork, basis_type=basis_type)
      ldab = ldwork
      ALLOCATE (cosab(ldab, ldab))
      ALLOCATE (sinab(ldab, ldab))
      ALLOCATE (work(ldwork, ldwork))

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)

      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set_a, basis_type=basis_type)
         IF (ASSOCIATED(basis_set_a)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set_a
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

!         WRITE(*,*) "first_sgfa",  basis_set_a%first_sgf
!         WRITE(*,*) "la_max",  basis_set_a%lmax
!         WRITE(*,*) "la_min",  basis_set_a%lmin
!         WRITE(*,*) "npgfa",  basis_set_a%npgf
!         WRITE(*,*) "nseta",  basis_set_a%nset
!         WRITE(*,*) "nsgfa", basis_set_a%nsgf_set
!         WRITE(*,*) "rpgfa", basis_set_a%pgf_radius
!         WRITE(*,*) "set_radius_a",  basis_set_a%set_radius
!         WRITE(*,*) "zeta", basis_set_a%zet
!         WRITE(*,*) "first_sgfb",  basis_set_b%first_sgf
!         WRITE(*,*) "lb_max",  basis_set_b%lmax
!         WRITE(*,*) "lb_min",  basis_set_b%lmin
!         WRITE(*,*) "npgfb",  basis_set_b%npgf
!         WRITE(*,*) "nsetb",  basis_set_b%nset
!         WRITE(*,*) "nsgfb", basis_set_b%nsgf_set
!         WRITE(*,*) "rpgfb", basis_set_b%pgf_radius
!         WRITE(*,*) "set_radius_b",  basis_set_b%set_radius
!         WRITE(*,*) "zetb", basis_set_b%zet
!         WRITE(*,*) ""


         ldsa = SIZE(sphi_a, 1)
         ldsb = SIZE(sphi_b, 1)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         NULLIFY (cblock)
         CALL dbcsr_get_block_p(matrix=cosmat, &
                                row=irow, col=icol, BLOCK=cblock, found=found)
         NULLIFY (sblock)
         CALL dbcsr_get_block_p(matrix=sinmat, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         IF (ASSOCIATED(cblock) .AND. .NOT. ASSOCIATED(sblock) .OR. &
             .NOT. ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN
            CPABORT("")
         ENDIF

         IF (ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN

            ra(:) = pbc(particle_set(iatom)%r(:), cell)
            rb(:) = ra+rab
            dab = SQRT(rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3))

            DO iset = 1, nseta

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)
               
               DO jset = 1, nsetb

                  IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)

                  ! Calculate the primitive integrals
                  CALL cossin(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                              lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                              ra, rb, kvec, cosab, sinab)
                  CALL contract_cossin(cblock, sblock, &
                                       iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                       jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                       cosab, sinab, ldab, work, ldwork)

               END DO
            END DO

         ENDIF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (cosab)
      DEALLOCATE (sinab)
      DEALLOCATE (work)
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE build_berry_moment_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param cosmat ...
!> \param sinmat ...
!> \param kvec ...
! **************************************************************************************************
   SUBROUTINE build_berry_kpoint_matrix(qs_env, cosmat, sinmat, kvec)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: cosmat, sinmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kvec

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_berry_kpoint_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, i, iatom, ic, icol, ikind, inode, irow, iset, jatom, jkind, jset, ldab, &
         ldsa, ldsb, ldwork, natom, ncoa, ncob, nimg, nkind, nseta, nsetb, sgfa, sgfb
      INTEGER, DIMENSION(3)                              :: icell
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgfa, nsgfb, row_blk_sizes
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, use_cell_mapping
      REAL(dp), DIMENSION(:, :), POINTER                 :: cblock, cosab, sblock, sinab, work
      REAL(KIND=dp)                                      :: dab
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, sphi_a, sphi_b, zeta, zetb
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER                  :: basis_set, basis_set_a, basis_set_b
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      dft_control=dft_control)
      nimg = dft_control%nimages
      IF (nimg > 1) THEN
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         use_cell_mapping = .TRUE.
      ELSE
         use_cell_mapping = .FALSE.
      END IF

      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      nkind = SIZE(qs_kind_set)
      natom = SIZE(particle_set)
      ALLOCATE (basis_set_list(nkind))
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind=qs_kind, basis_set=basis_set)
         IF (ASSOCIATED(basis_set)) THEN
            basis_set_list(ikind)%gto_basis_set => basis_set
         ELSE
            NULLIFY (basis_set_list(ikind)%gto_basis_set)
         END IF
      END DO

      ALLOCATE (row_blk_sizes(natom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, &
                            basis=basis_set_list)
      CALL get_ks_env(ks_env, dbcsr_dist=dbcsr_dist)
      ! (re)allocate matrix sets
      CALL dbcsr_allocate_matrix_set(sinmat, 1, nimg)
      CALL dbcsr_allocate_matrix_set(cosmat, 1, nimg)
      DO i = 1, nimg
         ! sin
         ALLOCATE (sinmat(1, i)%matrix)
         CALL dbcsr_create(matrix=sinmat(1, i)%matrix, &
                           name="SINMAT", &
                           dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                           nze=0)
         CALL cp_dbcsr_alloc_block_from_nbl(sinmat(1, i)%matrix, sab_orb)
         CALL dbcsr_set(sinmat(1, i)%matrix, 0.0_dp)
         ! cos
         ALLOCATE (cosmat(1, i)%matrix)
         CALL dbcsr_create(matrix=cosmat(1, i)%matrix, &
                           name="COSMAT", &
                           dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                           nze=0)
         CALL cp_dbcsr_alloc_block_from_nbl(cosmat(1, i)%matrix, sab_orb)
         CALL dbcsr_set(cosmat(1, i)%matrix, 0.0_dp)
      END DO

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=ldwork)
      ldab = ldwork
      ALLOCATE (cosab(ldab, ldab))
      ALLOCATE (sinab(ldab, ldab))
      ALLOCATE (work(ldwork, ldwork))

      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, inode=inode, &
                                iatom=iatom, jatom=jatom, r=rab, cell=icell)
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         ! basis ikind
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         sphi_a => basis_set_a%sphi
         zeta => basis_set_a%zet
         ! basis jkind
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         sphi_b => basis_set_b%sphi
         zetb => basis_set_b%zet

         ldsa = SIZE(sphi_a, 1)
         ldsb = SIZE(sphi_b, 1)

         IF (iatom <= jatom) THEN
            irow = iatom
            icol = jatom
         ELSE
            irow = jatom
            icol = iatom
         END IF

         IF (use_cell_mapping) THEN
            ic = cell_to_index(icell(1), icell(2), icell(3))
            CPASSERT(ic > 0)
         ELSE
            ic = 1
         END IF

         NULLIFY (sblock)
         CALL dbcsr_get_block_p(matrix=sinmat(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         CPASSERT(found)
         NULLIFY (cblock)
         CALL dbcsr_get_block_p(matrix=cosmat(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=cblock, found=found)
         CPASSERT(found)

         ra(:) = pbc(particle_set(iatom)%r(:), cell)
         rb(:) = ra+rab
         dab = SQRT(rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3))

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO jset = 1, nsetb

               IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(1, jset)

               ! Calculate the primitive integrals
               CALL cossin(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                           lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                           ra, rb, kvec, cosab, sinab)
               CALL contract_cossin(cblock, sblock, &
                                    iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                    jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                    cosab, sinab, ldab, work, ldwork)

            END DO
         END DO

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (cosab)
      DEALLOCATE (sinab)
      DEALLOCATE (work)
      DEALLOCATE (basis_set_list)
      DEALLOCATE (row_blk_sizes)

      CALL timestop(handle)

   END SUBROUTINE build_berry_kpoint_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param magnetic ...
!> \param nmoments ...
!> \param reference ...
!> \param ref_point ...
!> \param unit_number ...
! **************************************************************************************************
   SUBROUTINE qs_moment_berry_phase2(qs_env, magnetic, nmoments, reference, ref_point, unit_number)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: magnetic
      INTEGER, INTENT(IN)                                :: nmoments, reference
      REAL(dp), DIMENSION(:), POINTER                    :: ref_point
      INTEGER, INTENT(IN)                                :: unit_number

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_moment_berry_phase2', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)        :: rlab
      CHARACTER(LEN=default_string_length)               :: description
      COMPLEX(dp)                                        :: xphase(3), zdet, zdeta, zi(3), &
                                                            zij(3, 3), zijk(3, 3, 3), &
                                                            zijkl(3, 3, 3, 3), zphase(3), zz
      INTEGER                                            :: handle, i, ia, idim, ikind, ispin, ix, &
                                                            iy, iz, j, k, l, nao, nm, nmo, nmom, &
                                                            nmotot, tmp_dim
      LOGICAL                                            :: floating, ghost, uniform
      REAL(dp)                                           :: charge, ci(3), cij(3, 3), dd, occ, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mmom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rmom
      REAL(dp), DIMENSION(3)                             :: kvec, qq, rcc, ria
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER         :: eigrmat
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_new, opvec
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: op_fm_set
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp

!      REAL(dp)   :: normz1, rediv1, imdiv1
!      REAL(dp)   :: normz2, rediv2, imdiv2
!      COMPLEX(dp)                                        :: z1

      ! fbelle
      COMPLEX(dp)                                        :: d1(4,3), d2(4, 3), d3(4, 3, 3), &
                                                            d4(4, 3, 3, 3)
      TYPE(dbcsr_type), DIMENSION(:), POINTER            :: dkcos, dksin

   
      CPASSERT(ASSOCIATED(qs_env))

      WRITE(*,*) "qs_moment_berry_phase2 - BEGINNING"

      IF (ASSOCIATED(qs_env%ls_scf_env)) THEN
         IF (unit_number > 0) WRITE (unit_number, *) "Periodic moment calculation not implemented in linear scaling code"
         RETURN
      ENDIF

      CALL timeset(routineN, handle)

      ! restrict maximum moment available
      nmom = MIN(nmoments, 2)

      nm = (6+11*nmom+6*nmom**2+nmom**3)/6-1
      ! rmom(:,1)=electronic
      ! rmom(:,2)=nuclear
      ! rmom(:,3)=total
      ALLOCATE (rmom(nm+1, 3))
      ALLOCATE (rlab(nm+1))
      rmom = 0.0_dp
      rlab = ""
      IF (magnetic) THEN
         nm = 3
         ALLOCATE (mmom(nm))
         mmom = 0._dp
      END IF

      WRITE(*,*) "nmom", nmom
      WRITE(*,*) "nm", nm
      WRITE(*,*) "SIZE(rmom,1)", SIZE(rmom,1)
   

     ! fbelle
     ! DEBUGGING OF DERIVATIVES
     CALL debug_berry_moment_matrix(qs_env, nmoments = 1)


      NULLIFY (dft_control, rho, cell, particle_set, results, para_env, &
               local_particles, matrix_s, mos, rho_ao)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      rho=rho, &
                      cell=cell, &
                      results=results, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      para_env=para_env, &
                      local_particles=local_particles, &
                      matrix_s=matrix_s, &
                      mos=mos)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      NULLIFY (cosmat, sinmat)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL dbcsr_set(cosmat, 0.0_dp)
      CALL dbcsr_set(sinmat, 0.0_dp)

      ! Initialize and allocate matrices to store derivatives
      IF(nmom >= 1) THEN
         NULLIFY(dkcos, dksin)
         ALLOCATE(dkcos(3), dksin(3))
         DO i = 1, 3
            CALL dbcsr_copy(dkcos(i), matrix_s(1)%matrix, 'DKCOS MOM')
            CALL dbcsr_copy(dksin(i), matrix_s(1)%matrix, 'DKSIN MOM')
            CALL dbcsr_set(dkcos(i), 0.0_dp)
            CALL dbcsr_set(dksin(i), 0.0_dp)
         END DO
      END IF

!      WRITE(*,*) "cosmat"
!      CALL dbcsr_print(cosmat) 
!      CALL dbcsr_print(dkcos(1))
!      CALL dbcsr_print(dkcos(2))
!      CALL dbcsr_print(dkcos(3))

      ALLOCATE (op_fm_set(2, dft_control%nspins))
      ALLOCATE (opvec(dft_control%nspins))
      ALLOCATE (eigrmat(dft_control%nspins))
      nmotot = 0
      DO ispin = 1, dft_control%nspins
         NULLIFY (tmp_fm_struct, mo_coeff)
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
         nmotot = nmotot+nmo
         CALL cp_fm_create(opvec(ispin)%matrix, mo_coeff%matrix_struct)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nmo, &
                                  ncol_global=nmo, para_env=para_env, context=mo_coeff%matrix_struct%context)
         DO i = 1, SIZE(op_fm_set, 1)
            NULLIFY (op_fm_set(i, ispin)%matrix)
            CALL cp_fm_create(op_fm_set(i, ispin)%matrix, tmp_fm_struct)
         END DO
         CALL cp_cfm_create(eigrmat(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
      END DO

      ! occupation
      DO ispin = 1, dft_control%nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, maxocc=occ, uniform_occupation=uniform)
         IF (.NOT. uniform) THEN
            CPABORT("Berry phase moments for non uniform MOs' occupation numbers not implemented")
         END IF
      END DO

      ! reference point
      CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
      rcc = pbc(rcc, cell)

      ! label
      DO l = 1, nm
         ix = indco(1, l+1)
         iy = indco(2, l+1)
         iz = indco(3, l+1)
         CALL set_label(rlab(l+1), ix, iy, iz)
      END DO

      ! nuclear contribution
      DO ia = 1, SIZE(particle_set)
         atomic_kind => particle_set(ia)%atomic_kind
         CALL get_atomic_kind(atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)
         IF (.NOT. ghost .AND. .NOT. floating) THEN
            rmom(1, 2) = rmom(1, 2)-charge
         ENDIF
      END DO
      ria = twopi*MATMUL(cell%h_inv, rcc)
      zphase = CMPLX(COS(ria), SIN(ria), dp)**rmom(1, 2)

      zi = 0._dp
      zij = 0._dp
      zijk = 0._dp
      zijkl = 0._dp

      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
            WRITE(*,*) "Dipole - Nuclear"
            zi(:) = CMPLX(1._dp, 0._dp, dp)
            DO ia = 1, SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)
               IF (.NOT. ghost .AND. .NOT. floating) THEN
                  ria = particle_set(ia)%r
                  ria = pbc(ria, cell)
                  DO i = 1, 3
                     kvec(:) = twopi*cell%h_inv(i, :)
                     dd = SUM(kvec(:)*ria(:))
                     zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)**charge
                     zi(i) = zi(i)*zdeta
!                     WRITE(*,*) "ia", ia, "zi", zi
                  END DO
               ENDIF
            END DO
            zi = zi*zphase
            ci = AIMAG(LOG(zi))/twopi
            qq = AIMAG(LOG(zi))
            rmom(2:4, 2) = MATMUL(cell%hmat, ci)
         CASE (2)
            ! Quadrupole
!            CPABORT("Berry phase moments bigger than 1 not implemented")
            zij(:, :) = CMPLX(1._dp, 0._dp, dp)
            DO ia = 1, SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge)
               ria = particle_set(ia)%r
               ria = pbc(ria, cell)
               DO i = 1, 3
                  DO j = i, 3
                     kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :))
                     dd = SUM(kvec(:)*ria(:))
                     zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)**charge
                     zij(i, j) = zij(i, j)*zdeta
                     ! upper triangular elements of symmetric rank-2 tensor
                     zij(j, i) = zij(i, j)
                  END DO
               END DO
            END DO
!            DO i = 1, 3
!              WRITE(*,*) "zij(i,:)", zij(i,:) 
!            END DO


!            DO i = 1, 3
!               DO j = 1, 3
!                  zij(i, j) = zij(i, j)*zphase(i)*zphase(j)
!               END DO
!            END DO

            DO i = 1, 3
               DO j = 1, 3
!                  WRITE(*,*) "Indices", i, j
                  zij(i, j) = zij(i, j)*zphase(i)*zphase(j)
!                  WRITE(*,*) "zi(i)", zi(i)
!                  WRITE(*,*) "zi(j)", zi(j)
!                  WRITE(*,*) "zij(i,j)", zij(i,j)
!                  
!                  normz1 = Real(zi(i))*Real(zi(i)) + AIMAG(zi(i))*AIMAG(zi(i))
!                  rediv1 = (Real(zij(i,j))*Real(zi(i)) + AIMAG(zij(i,j))*AIMAG(zi(i)))/normz1 
!                  imdiv1 = (AIMAG(zij(i,j))*Real(zi(i)) - Real(zij(i,j))*AIMAG(zi(i)))/normz1
!                  
!                  WRITE(*,*) "normz1", normz1
!                  WRITE(*,*) "rediv1", rediv1
!                  WRITE(*,*) "imdiv1", imdiv1
!
!                  z1 = zij(i, j)/zi(i)
!                  WRITE(*,*) "z1", z1
!
!                  normz2 = Real(zi(j))*Real(zi(j)) + AIMAG(zi(j))*AIMAG(zi(j))
!                  rediv2 = (rediv1*Real(zi(j)) + imdiv1*AIMAG(zi(j)))/normz2 
!                  imdiv2 = (imdiv1*Real(zi(j)) - imdiv1*AIMAG(zi(j)))/normz2
!
!                  WRITE(*,*) "normz2", normz2
!                  WRITE(*,*) "rediv2", rediv2
!                  WRITE(*,*) "imdiv2", imdiv2
!
                  zz = zij(i, j)/zi(i)/zi(j)
!                  WRITE(*,*) "zz", zz
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
!               WRITE(*,* )  "zij(i, :)",  zij(i, :)
            END DO

            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))

!            DO i = 1, 3
!               WRITE(*,* )  "cij(i, :)",  cij(i, :)
!            END DO
            DO k = 4, 9
!               ! old implementation - does not correctly map
!               ix = indco(1, k+1)
!               iy = indco(2, k+1)
!               iz = indco(3, k+1)
!               WRITE(*,*) "k/ix/iy/iz", k, ix, iy, iz
!               IF (ix == 0) THEN
!                  rmom(k+1, 2) = cij(iy, iz)
!               ELSE IF (iy == 0) THEN
!                  rmom(k+1, 2) = cij(ix, iz)
!               ELSE IF (iz == 0) THEN
!                  rmom(k+1, 2) = cij(ix, iy)
!               END IF
            END DO

         ! Hardcoding of upper triangular matrix elements of cij
         ! into moment vector
         ! I am sorry
         rmom(5,2) = cij(1,1) 
         rmom(6,2) = cij(1,2) 
         rmom(7,2) = cij(1,3) 
         rmom(8,2) = cij(2,2) 
         rmom(9,2) = cij(2,3) 
         rmom(10,2) = cij(3,3) 

         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO

      ! electronic contribution

      ria = twopi*REAL(nmotot, dp)*occ*MATMUL(cell%h_inv, rcc)
      xphase = CMPLX(COS(ria), SIN(ria), dp)

      ! charge
      trace = 0.0_dp
      DO ispin = 1, dft_control%nspins
         CALL dbcsr_trace(rho_ao(ispin)%matrix, matrix_s(1)%matrix, trace)
         rmom(1, 1) = rmom(1, 1)+trace
      END DO

      zi = 0._dp
      zij = 0._dp
      zijk = 0._dp
      zijkl = 0._dp

      ! fraction components
      d1 = 0._dp
      d2 = 0._dp
      d3 = 0._dp
      d4 = 0._dp

      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
            DO i = 1, 3
               kvec(:) = twopi*cell%h_inv(i, :)
!               WRITE(*,*) "kvec(:)", kvec(:)
               CALL build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec, &
                                               nderivative = 1, dkcos = dkcos, dksin = dksin)
               !CALL build_berry_moment_matrix2(qs_env, cosmat, sinmat, kvec, &
               !                                nderivative = 1, dkcos = dkcos, dksin = dksin)
               IF (qs_env%run_rtp) THEN
                  CALL get_qs_env(qs_env, rtp=rtp)
                  CALL get_rtp(rtp, mos_new=mos_new)
                  CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
               ELSE
                  CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
               END IF
               zdet = CMPLX(1._dp, 0._dp, dp)
               DO ispin = 1, dft_control%nspins
                  CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=tmp_dim)
                  DO idim = 1, tmp_dim
                     eigrmat(ispin)%matrix%local_data(:, idim) = &
                        CMPLX(op_fm_set(1, ispin)%matrix%local_data(:, idim), &
                              -op_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
                  END DO
                  CALL cp_cfm_lu_decompose(eigrmat(ispin)%matrix, zdeta)
                  zdet = zdet*zdeta
                  ! multiply twice, account for double occupation
                  IF (dft_control%nspins == 1) THEN
                     zdet = zdet*zdeta
                  ENDIF
               END DO
               zi(i) = zdet
            END DO
            zi = zi*xphase
         CASE (2)
            ! Quadrupole
!            CPABORT("Berry phase moments bigger than 1 not implemented")
            DO i = 1, 3
               DO j = i, 3
                  kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :))
                  CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)
                  IF (qs_env%run_rtp) THEN
                     CALL get_qs_env(qs_env, rtp=rtp)
                     CALL get_rtp(rtp, mos_new=mos_new)
                     CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
                  ELSE
                     CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
                  END IF
                  zdet = CMPLX(1._dp, 0._dp, dp)
                  DO ispin = 1, dft_control%nspins
                     CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=tmp_dim)
                     DO idim = 1, tmp_dim
                        eigrmat(ispin)%matrix%local_data(:, idim) = &
                           ! exp(-ikr) = cos(kr) - i*sin(kr)
                           CMPLX(op_fm_set(1, ispin)%matrix%local_data(:, idim), &
                                 -op_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
                     END DO
                     CALL cp_cfm_lu_decompose(eigrmat(ispin)%matrix, zdeta)
                     zdet = zdet*zdeta
                     IF (dft_control%nspins == 1) THEN
                        zdet = zdet*zdeta
                     ENDIF
                  END DO
                  zij(i, j) = zdet*xphase(i)*xphase(j)
                  zij(j, i) = zdet*xphase(i)*xphase(j)
               END DO
            END DO
         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO
      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole (apply periodic (2 Pi) boundary conditions)
            ci = AIMAG(LOG(zi))
            DO i = 1, 3
               IF (qq(i)+ci(i) > pi) ci(i) = ci(i)-twopi
               IF (qq(i)+ci(i) < -pi) ci(i) = ci(i)+twopi
            END DO
            rmom(2:4, 1) = MATMUL(cell%hmat, ci)/twopi
         CASE (2)
            ! Quadrupole
!            CPABORT("Berry phase moments bigger than 1 not implemented")
            DO i = 1, 3
               DO j = 1, 3
                  zz = zij(i, j)/zi(i)/zi(j)
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
            END DO
            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))
! old code - coordinate mapping issue
!            DO k = 4, 9
!               ix = indco(1, k+1)
!               iy = indco(2, k+1)
!               iz = indco(3, k+1)
!               IF (ix == 0) THEN
!                  rmom(k+1, 1) = cij(iy, iz)
!               ELSE IF (iy == 0) THEN
!                  rmom(k+1, 1) = cij(ix, iz)
!               ELSE IF (iz == 0) THEN
!                  rmom(k+1, 1) = cij(ix, iy)
!               END IF
!            END DO
 
            ! Hardcoding of upper triangular matrix elements of cij
            ! into moment vector
            ! I am sorry
            rmom(5,2) = cij(1,1) 
            rmom(6,2) = cij(1,2) 
            rmom(7,2) = cij(1,3) 
            rmom(8,2) = cij(2,2) 
            rmom(9,2) = cij(2,3) 
            rmom(10,2) = cij(3,3) 
 
         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO

      rmom(:, 3) = rmom(:, 1)+rmom(:, 2)
      description = "[DIPOLE]"
      CALL cp_results_erase(results=results, description=description)
      CALL put_results(results=results, description=description, &
                       values=rmom(2:4, 3))
      IF (magnetic) THEN
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.TRUE., mmom=mmom)
      ELSE
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.TRUE.)
      END IF

      DEALLOCATE (rmom)
      DEALLOCATE (rlab)
      IF (magnetic) THEN
         DEALLOCATE (mmom)
      END IF

      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)

      DEALLOCATE(dkcos)
      DEALLOCATE(dksin)


      DO ispin = 1, dft_control%nspins
         CALL cp_fm_release(opvec(ispin)%matrix)
         CALL cp_cfm_release(eigrmat(ispin)%matrix)
         DO i = 1, SIZE(op_fm_set, 1)
            CALL cp_fm_release(op_fm_set(i, ispin)%matrix)
         END DO
      END DO
      DEALLOCATE (op_fm_set)
      DEALLOCATE (opvec)
      DEALLOCATE (eigrmat)

      CALL timestop(handle)

   END SUBROUTINE qs_moment_berry_phase2

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param magnetic ...
!> \param nmoments ...
!> \param reference ...
!> \param ref_point ...
!> \param unit_number ...
! **************************************************************************************************
   SUBROUTINE qs_moment_berry_phase(qs_env, magnetic, nmoments, reference, ref_point, unit_number)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: magnetic
      INTEGER, INTENT(IN)                                :: nmoments, reference
      REAL(dp), DIMENSION(:), POINTER                    :: ref_point
      INTEGER, INTENT(IN)                                :: unit_number

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_moment_berry_phase', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)        :: rlab
      CHARACTER(LEN=default_string_length)               :: description
      COMPLEX(dp)                                        :: xphase(3), zdet, zdeta, zi(3), &
                                                            zij(3, 3), zijk(3, 3, 3), &
                                                            zijkl(3, 3, 3, 3), zphase(3), zz
      INTEGER                                            :: handle, i, ia, idim, ikind, ispin, ix, &
                                                            iy, iz, j, k, l, nao, nm, nmo, nmom, &
                                                            nmotot, tmp_dim
      LOGICAL                                            :: floating, ghost, uniform
      REAL(dp)                                           :: charge, ci(3), cij(3, 3), dd, occ, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mmom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rmom
      REAL(dp), DIMENSION(3)                             :: kvec, qq, rcc, ria
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER         :: eigrmat
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_new, opvec
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: op_fm_set
      TYPE(cp_fm_struct_type), POINTER                   :: tmp_fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rt_prop_type), POINTER                        :: rtp

      REAL(dp)   :: normz1, rediv1, imdiv1
      REAL(dp)   :: normz2, rediv2, imdiv2
      COMPLEX(dp)                                        :: z1

      CPASSERT(ASSOCIATED(qs_env))

      WRITE(*,*) "qs_moment_berry_phase - BEGINNING"

      IF (ASSOCIATED(qs_env%ls_scf_env)) THEN
         IF (unit_number > 0) WRITE (unit_number, *) "Periodic moment calculation not implemented in linear scaling code"
         RETURN
      ENDIF

      CALL timeset(routineN, handle)

      ! restrict maximum moment available
      nmom = MIN(nmoments, 2)

      nm = (6+11*nmom+6*nmom**2+nmom**3)/6-1
      ! rmom(:,1)=electronic
      ! rmom(:,2)=nuclear
      ! rmom(:,3)=total
      ALLOCATE (rmom(nm+1, 3))
      ALLOCATE (rlab(nm+1))
      rmom = 0.0_dp
      rlab = ""
      IF (magnetic) THEN
         nm = 3
         ALLOCATE (mmom(nm))
         mmom = 0._dp
      END IF

      WRITE(*,*) "nmom", nmom
      WRITE(*,*) "nm", nm
      WRITE(*,*) "SIZE(rmom,1)", SIZE(rmom,1)
    

      NULLIFY (dft_control, rho, cell, particle_set, results, para_env, &
               local_particles, matrix_s, mos, rho_ao)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      rho=rho, &
                      cell=cell, &
                      results=results, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      para_env=para_env, &
                      local_particles=local_particles, &
                      matrix_s=matrix_s, &
                      mos=mos)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      NULLIFY (cosmat, sinmat)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL dbcsr_set(cosmat, 0.0_dp)
      CALL dbcsr_set(sinmat, 0.0_dp)

      ALLOCATE (op_fm_set(2, dft_control%nspins))
      ALLOCATE (opvec(dft_control%nspins))
      ALLOCATE (eigrmat(dft_control%nspins))
      nmotot = 0
      DO ispin = 1, dft_control%nspins
         NULLIFY (tmp_fm_struct, mo_coeff)
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
         nmotot = nmotot+nmo
         CALL cp_fm_create(opvec(ispin)%matrix, mo_coeff%matrix_struct)
         CALL cp_fm_struct_create(tmp_fm_struct, nrow_global=nmo, &
                                  ncol_global=nmo, para_env=para_env, context=mo_coeff%matrix_struct%context)
         DO i = 1, SIZE(op_fm_set, 1)
            NULLIFY (op_fm_set(i, ispin)%matrix)
            CALL cp_fm_create(op_fm_set(i, ispin)%matrix, tmp_fm_struct)
         END DO
         CALL cp_cfm_create(eigrmat(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
      END DO

      ! occupation
      DO ispin = 1, dft_control%nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, maxocc=occ, uniform_occupation=uniform)
         IF (.NOT. uniform) THEN
            CPABORT("Berry phase moments for non uniform MOs' occupation numbers not implemented")
         END IF
      END DO

      ! reference point
      CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
      rcc = pbc(rcc, cell)

      ! label
      DO l = 1, nm
         ix = indco(1, l+1)
         iy = indco(2, l+1)
         iz = indco(3, l+1)
         CALL set_label(rlab(l+1), ix, iy, iz)
      END DO

      ! nuclear contribution
      DO ia = 1, SIZE(particle_set)
         atomic_kind => particle_set(ia)%atomic_kind
         CALL get_atomic_kind(atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)
         IF (.NOT. ghost .AND. .NOT. floating) THEN
            rmom(1, 2) = rmom(1, 2)-charge
         ENDIF
      END DO
      ria = twopi*MATMUL(cell%h_inv, rcc)
      zphase = CMPLX(COS(ria), SIN(ria), dp)**rmom(1, 2)

      WRITE(*,*) "zphase", zphase

      zi = 0._dp
      zij = 0._dp
      zijk = 0._dp
      zijkl = 0._dp

      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
            WRITE(*,*) "Dipole - Nuclear"
            zi(:) = CMPLX(1._dp, 0._dp, dp)
            DO ia = 1, SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)
               IF (.NOT. ghost .AND. .NOT. floating) THEN
                  ria = particle_set(ia)%r
                  ria = pbc(ria, cell)
                  DO i = 1, 3
                     kvec(:) = twopi*cell%h_inv(i, :)
                     dd = SUM(kvec(:)*ria(:))
                     zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)**charge
                     zi(i) = zi(i)*zdeta
!                     WRITE(*,*) "ia", ia, "zi", zi
                  END DO
               ENDIF
            END DO
            zi = zi*zphase
            ci = AIMAG(LOG(zi))/twopi
            qq = AIMAG(LOG(zi))
            rmom(2:4, 2) = MATMUL(cell%hmat, ci)
         CASE (2)
            ! Quadrupole
!            CPABORT("Berry phase moments bigger than 1 not implemented")
            zij(:, :) = CMPLX(1._dp, 0._dp, dp)
            DO ia = 1, SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge)
               ria = particle_set(ia)%r
               ria = pbc(ria, cell)
               DO i = 1, 3
                  DO j = i, 3
                     kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :))
                     dd = SUM(kvec(:)*ria(:))
                     zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)**charge
                     zij(i, j) = zij(i, j)*zdeta
                     ! upper triangular elements of symmetric rank-2 tensor
                     zij(j, i) = zij(i, j)
                  END DO
               END DO
            END DO
!            DO i = 1, 3
!              WRITE(*,*) "zij(i,:)", zij(i,:) 
!            END DO


!            DO i = 1, 3
!               DO j = 1, 3
!                  zij(i, j) = zij(i, j)*zphase(i)*zphase(j)
!               END DO
!            END DO

            DO i = 1, 3
               DO j = 1, 3
                  WRITE(*,*) "Indices", i, j
                  zij(i, j) = zij(i, j)*zphase(i)*zphase(j)
                  WRITE(*,*) "zi(i)", zi(i)
                  WRITE(*,*) "zi(j)", zi(j)
                  WRITE(*,*) "zij(i,j)", zij(i,j)
                  
                  normz1 = Real(zi(i))*Real(zi(i)) + AIMAG(zi(i))*AIMAG(zi(i))
                  rediv1 = (Real(zij(i,j))*Real(zi(i)) + AIMAG(zij(i,j))*AIMAG(zi(i)))/normz1 
                  imdiv1 = (AIMAG(zij(i,j))*Real(zi(i)) - Real(zij(i,j))*AIMAG(zi(i)))/normz1
                  
                  WRITE(*,*) "normz1", normz1
                  WRITE(*,*) "rediv1", rediv1
                  WRITE(*,*) "imdiv1", imdiv1

                  z1 = zij(i, j)/zi(i)
                  WRITE(*,*) "z1", z1

                  normz2 = Real(zi(j))*Real(zi(j)) + AIMAG(zi(j))*AIMAG(zi(j))
                  rediv2 = (rediv1*Real(zi(j)) + imdiv1*AIMAG(zi(j)))/normz2 
                  imdiv2 = (imdiv1*Real(zi(j)) - imdiv1*AIMAG(zi(j)))/normz2

                  WRITE(*,*) "normz2", normz2
                  WRITE(*,*) "rediv2", rediv2
                  WRITE(*,*) "imdiv2", imdiv2

                  zz = zij(i, j)/zi(i)/zi(j)
                  WRITE(*,*) "zz", zz
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
               WRITE(*,* )  "zij(i, :)",  zij(i, :)
            END DO





            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))

            DO i = 1, 3
               WRITE(*,* )  "cij(i, :)",  cij(i, :)
            END DO
            DO k = 4, 9
!               ! old implementation - does not correctly map
!               ix = indco(1, k+1)
!               iy = indco(2, k+1)
!               iz = indco(3, k+1)
!               WRITE(*,*) "k/ix/iy/iz", k, ix, iy, iz
!               IF (ix == 0) THEN
!                  rmom(k+1, 2) = cij(iy, iz)
!               ELSE IF (iy == 0) THEN
!                  rmom(k+1, 2) = cij(ix, iz)
!               ELSE IF (iz == 0) THEN
!                  rmom(k+1, 2) = cij(ix, iy)
!               END IF
            END DO

         ! Hardcoding of upper triangular matrix elements of cij
         ! into moment vector
         ! I am sorry
         rmom(5,2) = cij(1,1) 
         rmom(6,2) = cij(1,2) 
         rmom(7,2) = cij(1,3) 
         rmom(8,2) = cij(2,2) 
         rmom(9,2) = cij(2,3) 
         rmom(10,2) = cij(3,3) 

         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO

      ! electronic contribution

      ria = twopi*REAL(nmotot, dp)*occ*MATMUL(cell%h_inv, rcc)
      xphase = CMPLX(COS(ria), SIN(ria), dp)

      ! charge
      trace = 0.0_dp
      DO ispin = 1, dft_control%nspins
         CALL dbcsr_trace(rho_ao(ispin)%matrix, matrix_s(1)%matrix, trace)
         rmom(1, 1) = rmom(1, 1)+trace
      END DO

      zi = 0._dp
      zij = 0._dp
      zijk = 0._dp
      zijkl = 0._dp

      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
            DO i = 1, 3
               kvec(:) = twopi*cell%h_inv(i, :)
               CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)
               IF (qs_env%run_rtp) THEN
                  CALL get_qs_env(qs_env, rtp=rtp)
                  CALL get_rtp(rtp, mos_new=mos_new)
                  CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
               ELSE
                  CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
               END IF
               zdet = CMPLX(1._dp, 0._dp, dp)
               DO ispin = 1, dft_control%nspins
                  CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=tmp_dim)
                  DO idim = 1, tmp_dim
                     eigrmat(ispin)%matrix%local_data(:, idim) = &
                        CMPLX(op_fm_set(1, ispin)%matrix%local_data(:, idim), &
                              -op_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
                  END DO
                  CALL cp_cfm_lu_decompose(eigrmat(ispin)%matrix, zdeta)
                  zdet = zdet*zdeta
                  IF (dft_control%nspins == 1) THEN
                     zdet = zdet*zdeta
                  ENDIF
               END DO
               zi(i) = zdet
            END DO
            zi = zi*xphase
         CASE (2)
            ! Quadrupole
!            CPABORT("Berry phase moments bigger than 1 not implemented")
            DO i = 1, 3
               DO j = i, 3
                  kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :))
                  CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)
                  IF (qs_env%run_rtp) THEN
                     CALL get_qs_env(qs_env, rtp=rtp)
                     CALL get_rtp(rtp, mos_new=mos_new)
                     CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
                  ELSE
                     CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
                  END IF
                  zdet = CMPLX(1._dp, 0._dp, dp)
                  DO ispin = 1, dft_control%nspins
                     CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=tmp_dim)
                     DO idim = 1, tmp_dim
                        eigrmat(ispin)%matrix%local_data(:, idim) = &
                           CMPLX(op_fm_set(1, ispin)%matrix%local_data(:, idim), &
                                 -op_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
                     END DO
                     CALL cp_cfm_lu_decompose(eigrmat(ispin)%matrix, zdeta)
                     zdet = zdet*zdeta
                     IF (dft_control%nspins == 1) THEN
                        zdet = zdet*zdeta
                     ENDIF
                  END DO
                  zij(i, j) = zdet*xphase(i)*xphase(j)
                  zij(j, i) = zdet*xphase(i)*xphase(j)
               END DO
            END DO
         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO
      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole (apply periodic (2 Pi) boundary conditions)
            ci = AIMAG(LOG(zi))
            DO i = 1, 3
               IF (qq(i)+ci(i) > pi) ci(i) = ci(i)-twopi
               IF (qq(i)+ci(i) < -pi) ci(i) = ci(i)+twopi
            END DO
            rmom(2:4, 1) = MATMUL(cell%hmat, ci)/twopi
         CASE (2)
            ! Quadrupole
!            CPABORT("Berry phase moments bigger than 1 not implemented")
            DO i = 1, 3
               DO j = 1, 3
                  zz = zij(i, j)/zi(i)/zi(j)
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
            END DO
            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))
! old code - coordinate mapping issue
!            DO k = 4, 9
!               ix = indco(1, k+1)
!               iy = indco(2, k+1)
!               iz = indco(3, k+1)
!               IF (ix == 0) THEN
!                  rmom(k+1, 1) = cij(iy, iz)
!               ELSE IF (iy == 0) THEN
!                  rmom(k+1, 1) = cij(ix, iz)
!               ELSE IF (iz == 0) THEN
!                  rmom(k+1, 1) = cij(ix, iy)
!               END IF
!            END DO
 
            ! Hardcoding of upper triangular matrix elements of cij
            ! into moment vector
            ! I am sorry
            rmom(5,2) = cij(1,1) 
            rmom(6,2) = cij(1,2) 
            rmom(7,2) = cij(1,3) 
            rmom(8,2) = cij(2,2) 
            rmom(9,2) = cij(2,3) 
            rmom(10,2) = cij(3,3) 
 
         CASE (3)
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
         CASE DEFAULT
            CPABORT("Berry phase moments bigger than 4 not implemented")
         END SELECT
      END DO

      rmom(:, 3) = rmom(:, 1)+rmom(:, 2)
      description = "[DIPOLE]"
      CALL cp_results_erase(results=results, description=description)
      CALL put_results(results=results, description=description, &
                       values=rmom(2:4, 3))
      IF (magnetic) THEN
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.TRUE., mmom=mmom)
      ELSE
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.TRUE.)
      END IF

      DEALLOCATE (rmom)
      DEALLOCATE (rlab)
      IF (magnetic) THEN
         DEALLOCATE (mmom)
      END IF

      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)

      DO ispin = 1, dft_control%nspins
         CALL cp_fm_release(opvec(ispin)%matrix)
         CALL cp_cfm_release(eigrmat(ispin)%matrix)
         DO i = 1, SIZE(op_fm_set, 1)
            CALL cp_fm_release(op_fm_set(i, ispin)%matrix)
         END DO
      END DO
      DEALLOCATE (op_fm_set)
      DEALLOCATE (opvec)
      DEALLOCATE (eigrmat)

      CALL timestop(handle)

   END SUBROUTINE qs_moment_berry_phase

! **************************************************************************************************
!> \brief ...
!> \param cosmat ...
!> \param sinmat ...
!> \param mos ...
!> \param op_fm_set ...
!> \param opvec ...
! **************************************************************************************************
   SUBROUTINE op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)

      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: op_fm_set
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: opvec

      INTEGER                                            :: i, nao, nmo
      TYPE(cp_fm_type), POINTER                          :: mo_coeff

      DO i = 1, SIZE(op_fm_set, 2) ! spin
         CALL get_mo_set(mo_set=mos(i)%mo_set, nao=nao, mo_coeff=mo_coeff, nmo=nmo)
         CALL cp_dbcsr_sm_fm_multiply(cosmat, mo_coeff, opvec(i)%matrix, ncol=nmo)
         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mo_coeff, opvec(i)%matrix, 0.0_dp, &
                      op_fm_set(1, i)%matrix)
         CALL cp_dbcsr_sm_fm_multiply(sinmat, mo_coeff, opvec(i)%matrix, ncol=nmo)
         CALL cp_gemm("T", "N", nmo, nmo, nao, 1.0_dp, mo_coeff, opvec(i)%matrix, 0.0_dp, &
                      op_fm_set(2, i)%matrix)
      ENDDO

   END SUBROUTINE op_orbbas

! **************************************************************************************************
!> \brief ...
!> \param cosmat ...
!> \param sinmat ...
!> \param mos ...
!> \param op_fm_set ...
!> \param mos_new ...
! **************************************************************************************************
   SUBROUTINE op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)

      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: op_fm_set
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_new

      INTEGER                                            :: i, icol, lcol, nao, newdim, nmo
      LOGICAL                                            :: double_col, double_row
      TYPE(cp_fm_struct_type), POINTER                   :: newstruct, newstruct1
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, work, work1, work2

      DO i = 1, SIZE(op_fm_set, 2) ! spin
         CALL get_mo_set(mo_set=mos(i)%mo_set, nao=nao, mo_coeff=mo_coeff, nmo=nmo)
         CALL cp_fm_get_info(mos_new(2*i)%matrix, ncol_local=lcol, ncol_global=nmo)
         double_col = .TRUE.
         double_row = .FALSE.
         CALL cp_fm_struct_double(newstruct, &
                                  mos_new(2*i)%matrix%matrix_struct, &
                                  mos_new(2*i)%matrix%matrix_struct%context, &
                                  double_col, &
                                  double_row)

         CALL cp_fm_create(work, matrix_struct=newstruct)
         CALL cp_fm_create(work1, matrix_struct=newstruct)
         CALL cp_fm_create(work2, matrix_struct=newstruct)
         CALL cp_fm_get_info(work, ncol_global=newdim)

         CALL cp_fm_set_all(work, 0.0_dp, 0.0_dp)
         DO icol = 1, lcol
            work%local_data(:, icol) = mos_new(2*i-1)%matrix%local_data(:, icol)
            work%local_data(:, icol+lcol) = mos_new(2*i)%matrix%local_data(:, icol)
         END DO

         CALL cp_dbcsr_sm_fm_multiply(cosmat, work, work1, ncol=newdim)
         CALL cp_dbcsr_sm_fm_multiply(sinmat, work, work2, ncol=newdim)

         DO icol = 1, lcol
            work%local_data(:, icol) = work1%local_data(:, icol)-work2%local_data(:, icol+lcol)
            work%local_data(:, icol+lcol) = work1%local_data(:, icol+lcol)+work2%local_data(:, icol)
         END DO

         CALL cp_fm_release(work1)
         CALL cp_fm_release(work2)

         CALL cp_fm_struct_double(newstruct1, &
                                  op_fm_set(1, i)%matrix%matrix_struct, &
                                  op_fm_set(1, i)%matrix%matrix_struct%context, &
                                  double_col, &
                                  double_row)

         CALL cp_fm_create(work1, matrix_struct=newstruct1)

         CALL cp_gemm("T", "N", nmo, newdim, nao, 1.0_dp, mos_new(2*i-1)%matrix, &
                      work, 0.0_dp, work1)

         DO icol = 1, lcol
            op_fm_set(1, i)%matrix%local_data(:, icol) = work1%local_data(:, icol)
            op_fm_set(2, i)%matrix%local_data(:, icol) = work1%local_data(:, icol+lcol)
         END DO

         CALL cp_gemm("T", "N", nmo, newdim, nao, 1.0_dp, mos_new(2*i)%matrix, &
                      work, 0.0_dp, work1)

         DO icol = 1, lcol
            op_fm_set(1, i)%matrix%local_data(:, icol) = &
               op_fm_set(1, i)%matrix%local_data(:, icol)+work1%local_data(:, icol+lcol)
            op_fm_set(2, i)%matrix%local_data(:, icol) = &
               op_fm_set(2, i)%matrix%local_data(:, icol)-work1%local_data(:, icol)
         END DO

         CALL cp_fm_release(work)
         CALL cp_fm_release(work1)
         CALL cp_fm_struct_release(newstruct)
         CALL cp_fm_struct_release(newstruct1)

      ENDDO

   END SUBROUTINE op_orbbas_rtp

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param magnetic ...
!> \param nmoments ...
!> \param reference ...
!> \param ref_point ...
!> \param unit_number ...
! **************************************************************************************************
   SUBROUTINE qs_moment_locop(qs_env, magnetic, nmoments, reference, ref_point, unit_number)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: magnetic
      INTEGER, INTENT(IN)                                :: nmoments, reference
      REAL(dp), DIMENSION(:), POINTER                    :: ref_point
      INTEGER, INTENT(IN)                                :: unit_number

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_moment_locop', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)        :: rlab
      CHARACTER(LEN=default_string_length)               :: description
      INTEGER                                            :: akind, handle, i, ia, iatom, ikind, &
                                                            ispin, ix, iy, iz, l, nm, nmom
      REAL(dp)                                           :: charge, dd, strace, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mmom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rmom
      REAL(dp), DIMENSION(3)                             :: rcc, ria
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: magmom, matrix_s, moments, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      CPASSERT(ASSOCIATED(qs_env))

      CALL timeset(routineN, handle)

      WRITE(*,*) "qs_moment_locop - BEGINNING"

      ! reference point
      CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)

      ! only allow for moments up to maxl set by basis
      nmom = MIN(nmoments, current_maxl)
      ! electronic contribution
      NULLIFY (moments, matrix_s)
      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
      nm = (6+11*nmom+6*nmom**2+nmom**3)/6-1
      CALL dbcsr_allocate_matrix_set(moments, nm)
      DO i = 1, nm
         ALLOCATE (moments(i)%matrix)
         CALL dbcsr_copy(moments(i)%matrix, matrix_s(1)%matrix, "Moments")
         CALL dbcsr_set(moments(i)%matrix, 0.0_dp)
      END DO

      CALL build_local_moment_matrix(qs_env, moments, nmom, ref_point=rcc)

      NULLIFY (dft_control, rho, cell, particle_set, qs_kind_set, results, para_env, matrix_s, rho_ao)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      rho=rho, &
                      cell=cell, &
                      results=results, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      para_env=para_env, &
                      matrix_s=matrix_s)

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      nm = SIZE(moments)
      ALLOCATE (rmom(nm+1, 3))
      ALLOCATE (rlab(nm+1))
      rmom = 0.0_dp
      rlab = ""

      trace = 0.0_dp
      DO ispin = 1, dft_control%nspins
         CALL dbcsr_trace(rho_ao(ispin)%matrix, matrix_s(1)%matrix, trace)
         rmom(1, 1) = rmom(1, 1)+trace
      END DO

      DO i = 1, SIZE(moments)
         strace = 0._dp
         DO ispin = 1, dft_control%nspins
            CALL dbcsr_trace(rho_ao(ispin)%matrix, moments(i)%matrix, trace)
            strace = strace+trace
         END DO
         rmom(i+1, 1) = strace
      END DO

      CALL dbcsr_deallocate_matrix_set(moments)

      ! magnetic moments
      IF (magnetic) THEN
         NULLIFY (magmom)
         CALL build_local_magmom_matrix(qs_env, magmom, nmom, ref_point=rcc)
         nm = SIZE(magmom)
         ALLOCATE (mmom(nm))
         DO i = 1, SIZE(magmom)
            strace = 0._dp
            DO ispin = 1, dft_control%nspins
               CALL dbcsr_trace(rho_ao(ispin)%matrix, magmom(i)%matrix, trace)
               strace = strace+trace
            END DO
            mmom(i) = strace
         END DO
         CALL dbcsr_deallocate_matrix_set(magmom)
      END IF

      ! nuclear contribution
      CALL get_qs_env(qs_env=qs_env, &
                      local_particles=local_particles)
      DO ikind = 1, SIZE(local_particles%n_el)
         DO ia = 1, local_particles%n_el(ikind)
            iatom = local_particles%list(ikind)%array(ia)
            ! fold atomic positions back into unit cell
            ria = pbc(particle_set(iatom)%r-rcc, cell)+rcc
            ria = ria-rcc
            atomic_kind => particle_set(iatom)%atomic_kind
            CALL get_atomic_kind(atomic_kind, kind_number=akind)
            CALL get_qs_kind(qs_kind_set(akind), core_charge=charge)
            rmom(1, 2) = rmom(1, 2)-charge
            DO l = 1, nm
               ix = indco(1, l+1)
               iy = indco(2, l+1)
               iz = indco(3, l+1)
               dd = 1._dp
               IF (ix > 0) dd = dd*ria(1)**ix
               IF (iy > 0) dd = dd*ria(2)**iy
               IF (iz > 0) dd = dd*ria(3)**iz
               rmom(l+1, 2) = rmom(l+1, 2)-charge*dd
               CALL set_label(rlab(l+1), ix, iy, iz)
            END DO
         END DO
      END DO
      CALL mp_sum(rmom(:, 2), para_env%group)
      rmom(:, :) = -rmom(:, :)
      rmom(:, 3) = rmom(:, 1)+rmom(:, 2)

      description = "[DIPOLE]"
      CALL cp_results_erase(results=results, description=description)
      CALL put_results(results=results, description=description, &
                       values=rmom(2:4, 3))

     ! fbelle
     ! WTF there is no quadrupole for local yet     
     CPABORT("Berry phase moments bigger than 1 not implemented")
              
              
      IF (magnetic) THEN
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.FALSE., mmom=mmom)
      ELSE
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.FALSE.)
      END IF

      DEALLOCATE (rmom)
      DEALLOCATE (rlab)
      IF (magnetic) THEN
         DEALLOCATE (mmom)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_moment_locop

! **************************************************************************************************
!> \brief ...
!> \param label ...
!> \param ix ...
!> \param iy ...
!> \param iz ...
! **************************************************************************************************
   SUBROUTINE set_label(label, ix, iy, iz)
      CHARACTER(LEN=*), INTENT(OUT)                      :: label
      INTEGER, INTENT(IN)                                :: ix, iy, iz

      INTEGER                                            :: i

      label = ""
      DO i = 1, ix
         WRITE (label(i:), "(A1)") "X"
      END DO
      DO i = ix+1, ix+iy
         WRITE (label(i:), "(A1)") "Y"
      END DO
      DO i = ix+iy+1, ix+iy+iz
         WRITE (label(i:), "(A1)") "Z"
      END DO

   END SUBROUTINE set_label

! **************************************************************************************************
!> \brief ...
!> \param unit_number ...
!> \param nmom ...
!> \param rmom ...
!> \param rlab ...
!> \param rcc ...
!> \param cell ...
!> \param periodic ...
!> \param mmom ...
! **************************************************************************************************
   SUBROUTINE print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic, mmom)
      INTEGER, INTENT(IN)                                :: unit_number, nmom
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: rmom
      CHARACTER(LEN=8), DIMENSION(:)                     :: rlab
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: rcc
      TYPE(cell_type), POINTER                           :: cell
      LOGICAL                                            :: periodic
      REAL(dp), DIMENSION(:), INTENT(IN), OPTIONAL       :: mmom

      INTEGER                                            :: i, i0, i1, j, l
      REAL(dp)                                           :: dd

      IF (unit_number > 0) THEN
         DO l = 0, nmom
            SELECT CASE (l)
            CASE (0)
               WRITE (unit_number, "(T3,A,T32,3F14.8)") "Reference Point [Bohr]", rcc
               WRITE (unit_number, "(T3,A)") "Charges"
               WRITE (unit_number, "(T5,A,T18,F14.8,T36,A,T42,F14.8,T60,A,T68,F14.8)") &
                  "Electronic=", rmom(1, 1), "Core=", rmom(1, 2), "Total=", rmom(1, 3)
            CASE (1)
               IF (periodic) THEN
                  WRITE (unit_number, "(T3,A)") "Dipole vectors are based on the periodic (Berry phase) operator."
                  WRITE (unit_number, "(T3,A)") "They are defined modulo integer multiples of the cell matrix [Debye]."
                  WRITE (unit_number, "(T3,A,3(F14.8,1X),A)") "[X] [", cell%hmat(1, :)*debye, "] [i]"
                  WRITE (unit_number, "(T3,A,3(F14.8,1X),A)") "[Y]=[", cell%hmat(2, :)*debye, "]*[j]"
                  WRITE (unit_number, "(T3,A,3(F14.8,1X),A)") "[Z] [", cell%hmat(3, :)*debye, "] [k]"
               ELSE
                  WRITE (unit_number, "(T3,A)") "Dipoles are based on the traditional operator."
               ENDIF
               dd = SQRT(SUM(rmom(2:4, 3)**2))*debye
               WRITE (unit_number, "(T3,A)") "Dipole moment [Debye]"
               WRITE (unit_number, "(T5,3(A,A,F14.8,1X),T60,A,T68,F14.8)") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye, i=2, 4), "Total=", dd
            CASE (2)
               WRITE (unit_number, "(T3,A)") "Quadrupole moment [Debye*Angstrom]"
               WRITE (unit_number, "(T17,3(A,A,F14.8,9X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr, i=5, 7)
               WRITE (unit_number, "(T17,3(A,A,F14.8,9X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr, i=8, 10)
            CASE (3)
               WRITE (unit_number, "(T3,A)") "Octapole moment [Debye*Angstrom**2]"
               WRITE (unit_number, "(T7,4(A,A,F14.8,3X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr, i=11, 14)
               WRITE (unit_number, "(T7,4(A,A,F14.8,3X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr, i=15, 18)
               WRITE (unit_number, "(T7,4(A,A,F14.8,3X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr, i=19, 20)
            CASE (4)
               WRITE (unit_number, "(T3,A)") "Hexadecapole moment [Debye*Angstrom**3]"
               WRITE (unit_number, "(T6,4(A,A,F14.8,2X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr/bohr, i=21, 24)
               WRITE (unit_number, "(T6,4(A,A,F14.8,2X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr/bohr, i=25, 28)
               WRITE (unit_number, "(T6,4(A,A,F14.8,2X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr/bohr, i=29, 32)
               WRITE (unit_number, "(T6,4(A,A,F14.8,2X))") &
                  (TRIM(rlab(i)), "=", rmom(i, 3)*debye/bohr/bohr/bohr, i=32, 35)
            CASE DEFAULT
               WRITE (unit_number, "(T3,A,A,I2)") "Higher moment [Debye*Angstrom**(L-1)]", &
                  "  L=", l
               i0 = (6+11*(l-1)+6*(l-1)**2+(l-1)**3)/6
               i1 = (6+11*l+6*l**2+l**3)/6-1
               dd = debye/(bohr)**(l-1)
               DO i = i0, i1, 3
                  WRITE (unit_number, "(T18,3(A,A,F14.8,4X))") &
                     (TRIM(rlab(j+1)), "=", rmom(j+1, 3)*dd, j=i, MIN(i1, i+2))
               END DO
            END SELECT
         END DO
         IF (PRESENT(mmom)) THEN
            IF (nmom >= 1) THEN
               dd = SQRT(SUM(mmom(1:3)**2))
               WRITE (unit_number, "(T3,A)") "Magnetic Dipole Moment [WRONG]"
               WRITE (unit_number, "(T5,3(A,A,F14.8,1X),T60,A,T68,F14.8)") &
                  (TRIM(rlab(i+1)), "=", mmom(i), i=1, 3), "Total=", dd
            END IF
         END IF
      END IF

   END SUBROUTINE print_moments

END MODULE qs_moments

