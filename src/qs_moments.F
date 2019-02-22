!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
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
        dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, &
        dbcsr_distribution_type, dbcsr_get_block_p, dbcsr_multiply, dbcsr_p_type, dbcsr_set, &
        dbcsr_trace, dbcsr_type, dbcsr_type_symmetric
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
!fbelle
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type 
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm, &
                                              copy_fm_to_dbcsr, &
                                              cp_fm_to_dbcsr_row_template
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale, &
                                              cp_fm_triangular_invert, &
                                              cp_fm_trace, &
                                              cp_fm_gemm
   USE qs_mo_methods,                   ONLY: calculate_density_matrix
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_moments'

   ! Public subroutines
   PUBLIC :: build_berry_moment_matrix, build_local_moment_matrix
   PUBLIC :: build_berry_moment_matrix2
   PUBLIC :: build_berry_moment_matrix3
   PUBLIC :: build_berry_moment_matrix4
!   PUBLIC :: debug_berry_moment_matrix
   PUBLIC :: build_berry_kpoint_matrix
   PUBLIC :: qs_moment_berry_phase, qs_moment_locop
   PUBLIC :: qs_moment_berry_phase2
   PUBLIC :: qs_moment_berry_phase3
   PUBLIC :: qs_moment_berry_phase4
   PUBLIC :: numer_berry_moment_matrix

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
         POINTER                                         :: sab_all
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      IF (nmoments < 1) RETURN

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set, cell, particle_set, sab_all)
      NULLIFY (matrix_s)

      CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)

      ! magnetic dipoles/angular moments only
      nm = 3
      CALL dbcsr_allocate_matrix_set(magmom, nm)

      DO i = 1, nm
         ALLOCATE (magmom(i)%matrix)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, magmom(i)%matrix)
         CALL dbcsr_set(magmom(i)%matrix, 0.0_dp)
      END DO

      NULLIFY (qs_kind_set, particle_set, sab_all, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_all=sab_all)

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
      CALL neighbor_list_iterator_create(nl_iterator, sab_all)
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

         irow = iatom
         icol = jatom

         DO i = 1, nm
            NULLIFY (mint(i)%block)
            CALL dbcsr_get_block_p(matrix=magmom(i)%matrix, &
                                   row=irow, col=icol, BLOCK=mint(i)%block, found=found)
            mint(i)%block = 0._dp
            CPASSERT(ASSOCIATED(mint(i)%block))
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

                  CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1), SIZE(work, 1), &
                             1.0_dp, mint(i)%block(sgfa, sgfb), &
                             SIZE(mint(i)%block, 1))

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
!> \brief NUMERICAL CALCULATION OF FIRST DERIVATIVE 
!>        
!> \param qs_env ...
!> \param cosmat ...
!> \param sinmat ...
!> \param kvec ...
!> \param sab_orb_external ...
!> \param basis_type ...
!> \param nderivative ... Order of derivative wrt k
! **************************************************************************************************
   SUBROUTINE numer_berry_moment_matrix(qs_env, kvec, nmoments, dkc, dks)
   
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nmoments
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER, OPTIONAL :: dkc, dks

      CHARACTER(LEN=*), PARAMETER :: routineN = 'numer_berry_moment_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, j, k, nm, nmom
      REAL(dp), DIMENSION(3)                             :: kvec
      REAL(dp)                                           :: std_value
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat

      ! fbelle
      REAL(KIND=dp)                                      :: dx   

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER        :: numer_cosmat, numer_sinmat
      TYPE(cp_para_env_type), POINTER                    :: para_env


      CALL timeset(routineN, handle)

      ! restrict maximum moment available
      nmom = MIN(nmoments, 1)

      nm = (6+11*nmom+6*nmom**2+nmom**3)/6-1

      NULLIFY (cell, matrix_s)

      CALL get_qs_env(qs_env, &
                      cell = cell, &
                       matrix_s = matrix_s)
      CALL get_qs_env(qs_env, para_env=para_env)


      NULLIFY (cosmat, sinmat)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL dbcsr_set(cosmat, 0.0_dp)
      CALL dbcsr_set(sinmat, 0.0_dp)

      ! Init matrices to store numerically calculated matrices
      NULLIFY (numer_cosmat, numer_sinmat)
      CALL dbcsr_allocate_matrix_set(numer_cosmat, 3)
      CALL dbcsr_allocate_matrix_set(numer_sinmat, 3)

      DO j = 1, 3
         ALLOCATE (numer_cosmat(j)%matrix, numer_sinmat(j)%matrix)
         CALL dbcsr_copy(numer_cosmat(j)%matrix, matrix_s(1)%matrix, 'COS NUMER')
         CALL dbcsr_copy(numer_sinmat(j)%matrix, matrix_s(1)%matrix, 'SIN NUMER')
         CALL dbcsr_set(numer_cosmat(j)%matrix, 0.0_dp)
         CALL dbcsr_set(numer_sinmat(j)%matrix, 0.0_dp)
      END DO

      ! Deviation
      dx = 1.0E-7

      !------------------------------------------------------------------------
      ! Numerical calculation 
      !------------------------------------------------------------------------

      IF(para_env%mepos==0) THEN
!         WRITE(*,*) "----------------------------------------------"
!         WRITE(*,*) "NUMERIC DERIVATIVE"
!         WRITE(*,*) "----------------------------------------------"
      END IF

      ! loop over the 3 derivation directions (k is 3D vector)
      DO k = 1, 3
         ! save inital value of kvec in direction k
         std_value = kvec(k)
         CALL dbcsr_set(cosmat, 0.0_dp)
         CALL dbcsr_set(sinmat, 0.0_dp)
         DO j = 1, 2
            ! Apply deviation on kvecs direction k
            kvec(k) = std_value-(-1.0_dp)**j*dx
            ! Compute cosmat and sinmat
            CALL build_berry_moment_matrix(qs_env = qs_env, &
                                           cosmat = cosmat, &
                                           sinmat = sinmat, &
                                           kvec = kvec)

            IF (j == 1) THEN
               CALL dbcsr_copy(numer_cosmat(k)%matrix, cosmat)
               CALL dbcsr_copy(numer_sinmat(k)%matrix, sinmat)
            ELSE
               CALL dbcsr_add(numer_cosmat(k)%matrix, cosmat, 1.0_dp, -1.0_dp)
               CALL dbcsr_add(numer_sinmat(k)%matrix, sinmat, 1.0_dp, -1.0_dp)
            END IF

         END DO
         kvec(k) = std_value

         ! dividing by 0.5*dx
         CALL dbcsr_scale(numer_cosmat(k)%matrix, 0.5_dp/dx)
         CALL dbcsr_scale(numer_sinmat(k)%matrix, 0.5_dp/dx)

         CALL dbcsr_copy(dkc(k)%matrix, numer_cosmat(k)%matrix)
         CALL dbcsr_copy(dks(k)%matrix, numer_sinmat(k)%matrix)

      END DO ! k

      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)

      CALL dbcsr_deallocate_matrix_set(numer_cosmat)
      CALL dbcsr_deallocate_matrix_set(numer_sinmat)

      CALL timestop(handle)

      IF(para_env%mepos==0) THEN
         WRITE(*,*) "numer_berry_moment_matrix - END"      
      END IF

   END SUBROUTINE



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

      INTEGER                                            :: handle, i, j, k, l, nm, nmom
      REAL(dp), DIMENSION(3)                             :: kvec
      REAL(dp)                                           :: std_value
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat

      ! fbelle
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dkcos, dksin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: d2kcos, d2ksin
      REAL(KIND=dp)                                      :: dx   

      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: analyt_cosmat, analyt_sinmat
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: numer_cosmat, numer_sinmat
      REAL (KIND=dp), DIMENSION(3,3)                     :: anactrace, anastrace
      REAL (KIND=dp), DIMENSION(3,3)                     :: numctrace, numstrace
      REAL (KIND=dp), DIMENSION(3,3)                     :: diffcnorm, diffsnorm
      TYPE(cp_para_env_type), POINTER                    :: para_env


      CALL timeset(routineN, handle)

      ! restrict maximum moment available
      nmom = MIN(nmoments, 2)

      nm = (6+11*nmom+6*nmom**2+nmom**3)/6-1

      NULLIFY (cell, matrix_s)

      CALL get_qs_env(qs_env, &
                      cell = cell, &
                       matrix_s = matrix_s)
      CALL get_qs_env(qs_env, para_env=para_env)


      NULLIFY (cosmat, sinmat)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL dbcsr_set(cosmat, 0.0_dp)
      CALL dbcsr_set(sinmat, 0.0_dp)

      ! Init matrices to store analytically calculated matrices
      NULLIFY (analyt_cosmat, analyt_sinmat)
      CALL dbcsr_allocate_matrix_set(analyt_cosmat, 3, 3)
      CALL dbcsr_allocate_matrix_set(analyt_sinmat, 3, 3)
      ! Init matrices to store numerically calculated matrices
      NULLIFY (numer_cosmat, numer_sinmat)
      CALL dbcsr_allocate_matrix_set(numer_cosmat, 3, 3)
      CALL dbcsr_allocate_matrix_set(numer_sinmat, 3, 3)

      DO i = 1, 3
         DO j = 1, 3
            ALLOCATE (analyt_cosmat(i,j)%matrix, analyt_sinmat(i,j)%matrix)
            CALL dbcsr_copy(analyt_cosmat(i,j)%matrix, matrix_s(1)%matrix, 'COS ANAL')
            CALL dbcsr_copy(analyt_sinmat(i,j)%matrix, matrix_s(1)%matrix, 'SIN ANAL')
            CALL dbcsr_set(analyt_cosmat(i,j)%matrix, 0.0_dp)
            CALL dbcsr_set(analyt_sinmat(i,j)%matrix, 0.0_dp)

            ALLOCATE (numer_cosmat(i,j)%matrix, numer_sinmat(i,j)%matrix)
            CALL dbcsr_copy(numer_cosmat(i,j)%matrix, matrix_s(1)%matrix, 'COS NUMER')
            CALL dbcsr_copy(numer_sinmat(i,j)%matrix, matrix_s(1)%matrix, 'SIN NUMER')
            CALL dbcsr_set(numer_cosmat(i,j)%matrix, 0.0_dp)
            CALL dbcsr_set(numer_sinmat(i,j)%matrix, 0.0_dp)
         END DO
      END DO         

      ! Initialize and allocate matrices to store derivatives
      IF(nmom >= 1) THEN

         NULLIFY(dkcos, dksin)
         CALL dbcsr_allocate_matrix_set(dkcos, 3)
         CALL dbcsr_allocate_matrix_set(dksin, 3)
         DO i = 1, 3
            ALLOCATE(dkcos(i)%matrix, dksin(i)%matrix)
            CALL dbcsr_copy(dkcos(i)%matrix, matrix_s(1)%matrix, 'DKCOS MOM')
            CALL dbcsr_copy(dksin(i)%matrix, matrix_s(1)%matrix, 'DKSIN MOM')
            CALL dbcsr_set(dkcos(i)%matrix, 0.0_dp)
            CALL dbcsr_set(dksin(i)%matrix, 0.0_dp)
         END DO
      END IF

      IF(nmom >= 2) THEN
         NULLIFY(d2kcos, d2ksin)
         CALL dbcsr_allocate_matrix_set(dkcos, 6)
         CALL dbcsr_allocate_matrix_set(dksin, 6)
         DO i = 1, 6
            ALLOCATE(d2kcos(i)%matrix, d2ksin(i)%matrix)
            CALL dbcsr_copy(d2kcos(i)%matrix, matrix_s(1)%matrix, 'D2KCOS MOM')
            CALL dbcsr_copy(d2ksin(i)%matrix, matrix_s(1)%matrix, 'D2KSIN MOM')
            CALL dbcsr_set(d2kcos(i)%matrix, 0.0_dp)
            CALL dbcsr_set(d2ksin(i)%matrix, 0.0_dp)
         END DO
      END IF


      ! Deviation
      dx = 1.0E-5

      !------------------------------------------------------------------------
      ! Analytic call of cossin
      !------------------------------------------------------------------------

      IF(para_env%mepos==0) THEN
         WRITE(*,*) "----------------------------------------------"
         WRITE(*,*) "ANALYTIC DERIVATIVE "
         WRITE(*,*) "----------------------------------------------"
      END IF

      ! 3 distinct kvecs -> each gives 3 first derivs
      !                  -> each gives 6 second derivs

      ! First, evaluate analytic 
      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
            DO i = 1, 3
               kvec(:) = twopi*cell%h_inv(i, :)
               CALL build_berry_moment_matrix2(qs_env, cosmat, sinmat, kvec, &
                                               nderivative = 1, dkcos = dkcos, dksin = dksin)
               DO j = 1, 3
                  ! copy cosmat sinmat into separate array
                  CALL dbcsr_copy(analyt_cosmat(i, j)%matrix, dkcos(j)%matrix) 
                  CALL dbcsr_copy(analyt_sinmat(i, j)%matrix, dksin(j)%matrix) 

                  ! store only trace
                  CALL dbcsr_trace(analyt_cosmat(i, j)%matrix, anactrace(i, j))
                  CALL dbcsr_trace(analyt_sinmat(i, j)%matrix, anastrace(i, j))
               END DO
            END DO
         CASE (2)
            ! Quadrupole
            DO i = 1, 3
               DO j = i, 3
                  kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :))
!                  CALL build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec, & 
!                                                  nderivative = 2, dkcos = dkcos, dksin = dksin, &
!                                                  d2kcos = d2kcos, d2ksin = d2ksin)
                  !DO k = 1, 6
                  !   ! copy cosmat sinmat into separate array
                  !   CALL dbcsr_copy(analyt_cosmat(i, j), d2kcos(j)) 
                  !   CALL dbcsr_copy(analyt_sinmat(i, j), d2ksin(j)) 

                  !   ! store only trace
                  !   CALL dbcsr_trace(analyt_cosmat(i, j), anactrace(i, j))
                  !   CALL dbcsr_trace(analyt_sinmat(i, j), anastrace(i, j))
                  !END DO

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

      ! Traces of first derivatives
      IF(para_env%mepos==0) THEN
         WRITE(*,*) "TRACE OF #1 ANALYTIC DERIVATIVE"
         WRITE(*,*) "TRACE dcos:", anactrace(1,1:3)
         WRITE(*,*) "TRACE dcos:", anactrace(2,1:3)
         WRITE(*,*) "TRACE dcos:", anactrace(3,1:3)
         WRITE(*,*) ""
         WRITE(*,*) "TRACE dsin:", anastrace(1,1:3)
         WRITE(*,*) "TRACE dsin:", anastrace(2,1:3)
         WRITE(*,*) "TRACE dsin:", anastrace(3,1:3)
         WRITE(*,*) "----------------------------------------------"
      END IF

      !------------------------------------------------------------------------
      ! Numerical calculation 
      !------------------------------------------------------------------------

      IF(para_env%mepos==0) THEN
         WRITE(*,*) "----------------------------------------------"
         WRITE(*,*) "NUMERIC DERIVATIVE" 
         WRITE(*,*) "----------------------------------------------"
      END IF

      ! kvec is a 3D object, needs to be perturbed in all 3 directions
      ! just like particle position in forces
      DO i = 1, 3
         ! do each kvec twice: (+) deviation and (-) deviation
         kvec(:) = twopi*cell%h_inv(i, :)
         ! loop over the 3 derivation directions (k is 3D vector)
         DO k = 1, 3
            ! save inital value of kvec in direction k
            std_value = kvec(k)
            CALL dbcsr_set(cosmat, 0.0_dp)
            CALL dbcsr_set(sinmat, 0.0_dp)
            DO j = 1, 2
               ! Apply deviation on kvecs direction k
               kvec(k) = std_value-(-1.0_dp)**j*dx
               ! Compute cosmat and sinmat
               CALL build_berry_moment_matrix3(qs_env = qs_env, &
                                               cosmat = cosmat, &
                                               sinmat = sinmat, &
                                               kvec = kvec)

               IF (j == 1) THEN
                  CALL dbcsr_copy(numer_cosmat(i, k)%matrix, cosmat)
                  CALL dbcsr_copy(numer_sinmat(i, k)%matrix, sinmat)
               ELSE
                  CALL dbcsr_add(numer_cosmat(i, k)%matrix, cosmat, 1.0_dp, -1.0_dp)
                  CALL dbcsr_add(numer_sinmat(i, k)%matrix, sinmat, 1.0_dp, -1.0_dp)
               END IF

            END DO
            kvec(k) = std_value

            ! dividing by 0.5*dx
            CALL dbcsr_scale(numer_cosmat(i, k)%matrix, 0.5_dp/dx)
            CALL dbcsr_scale(numer_sinmat(i, k)%matrix, 0.5_dp/dx)

            ! store trace(i,k)
            CALL dbcsr_trace(numer_cosmat(i,k)%matrix, numctrace(i,k))
            CALL dbcsr_trace(numer_sinmat(i,k)%matrix, numstrace(i,k))
         END DO ! k
      END DO

      ! Norm of Differnce: NORM(ANALYT - NUMER)
      DO i = 1, 3
         DO k = 1, 3
            CALL dbcsr_add(analyt_cosmat(i,k)%matrix, numer_cosmat(i,k)%matrix, 1.0_dp, -1.0_dp)
            CALL dbcsr_norm(matrix = analyt_cosmat(i,k)%matrix, which_norm = 1, norm_scalar = diffcnorm(i,k)) 
            CALL dbcsr_add(analyt_sinmat(i,k)%matrix, numer_sinmat(i,k)%matrix, 1.0_dp, -1.0_dp)
            CALL dbcsr_norm(matrix = analyt_sinmat(i,k)%matrix, which_norm = 1, norm_scalar = diffsnorm(i,k)) 
         END DO
      END DO

      IF(para_env%mepos==0) THEN
         WRITE(*,*) "TRACE OF #1 NUMERICAL DERIVATIVE"
         WRITE(*,*) "TRACE dcos:", numctrace(1,1:3)
         WRITE(*,*) "TRACE dcos:", numctrace(2,1:3)
         WRITE(*,*) "TRACE dcos:", numctrace(3,1:3)
         WRITE(*,*) ""
         WRITE(*,*) "TRACE dsin:", numstrace(1,1:3)
         WRITE(*,*) "TRACE dsin:", numstrace(2,1:3)
         WRITE(*,*) "TRACE dsin:", numstrace(3,1:3)
         WRITE(*,*) "----------------------------------------------"
         WRITE(*,*) "MATRIX NORM( ANALYT - NUM )"
         WRITE(*,*) "", diffcnorm(1,:)
         WRITE(*,*) "", diffcnorm(2,:)
         WRITE(*,*) "", diffcnorm(3,:)
         WRITE(*,*) "" 
         WRITE(*,*) "", diffsnorm(1,:)
         WRITE(*,*) "", diffsnorm(2,:)
         WRITE(*,*) "", diffsnorm(3,:)
         WRITE(*,*) "----------------------------------------------"
      END IF

      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)

      IF(nmom >= 1) THEN
         CALL dbcsr_deallocate_matrix_set(dkcos)
         CALL dbcsr_deallocate_matrix_set(dksin)
      END IF

      IF(nmom >= 2) THEN
         CALL dbcsr_deallocate_matrix_set(d2kcos)
         CALL dbcsr_deallocate_matrix_set(d2ksin)
      END IF

      CALL dbcsr_deallocate_matrix_set(analyt_cosmat)
      CALL dbcsr_deallocate_matrix_set(analyt_sinmat)

      CALL dbcsr_deallocate_matrix_set(numer_cosmat)
      CALL dbcsr_deallocate_matrix_set(numer_sinmat)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Calculation of the moment matrix (a|exp(ikr)|b) 
!>        and their derivatives wrt a factor t in front of each k 
!>        to pull down the whole argument of exp(ikr) up to order nderivative
!> \param qs_env The Quickstep environment
!> \param cosmat Real part of (a|exp(ikr)|b) 
!> \param sinmat Imaginary part of (a|exp(ikr)|b) 
!> \param kvec The k-vector for which the integral is evaluated
!> \param sab_orb_external Specify an external neighborlist
!> \param basis_type ...
!> \param nderivative Order of derivative wrt t
!> \param dkcos dksin Real and imaginary part of first derivative wrt t
!> \param d2kcos d2ksin Real and imaginary part of second derivative wrt t
! **************************************************************************************************
   SUBROUTINE build_berry_moment_matrix4(qs_env, cosmat, sinmat, kvec, &
                                         sab_orb_external, basis_type, nderivative, &
                                         dkcos, dksin, d2kcos, d2ksin, d3kcos, d3ksin, d4kcos, d4ksin)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kvec
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_external
      CHARACTER(len=*), OPTIONAL                         :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: nderivative
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: dkcos, dksin, d2kcos, d2ksin 
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: d3kcos, d3ksin, d4kcos, d4ksin 

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_berry_moment_matrix4', &
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
      INTEGER                                            :: maxder, nder, k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: cosmint, sinmint
      REAL(dp), DIMENSION(:, :), POINTER                 :: ckblk, skblk, c2kblk, s2kblk
      REAL(dp), DIMENSION(:, :), POINTER                 :: c3kblk, s3kblk, c4kblk, s4kblk

      CALL timeset(routineN, handle)

      WRITE(*,*) "build_berry_moment_matrix4"

      NULLIFY (qs_kind_set, particle_set, sab_orb, cell)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, cell=cell, &
                      sab_orb=sab_orb)

      IF (PRESENT(sab_orb_external)) sab_orb => sab_orb_external

      CALL dbcsr_set(sinmat, 0.0_dp)
      CALL dbcsr_set(cosmat, 0.0_dp)

      IF (PRESENT(dkcos) .AND. PRESENT(dksin)) THEN
          CALL dbcsr_set(dkcos, 0.0_dp)
          CALL dbcsr_set(dksin, 0.0_dp)
      END IF
      IF (PRESENT(d2kcos) .AND. PRESENT(d2ksin)) THEN
         CALL dbcsr_set(d2kcos, 0.0_dp)
         CALL dbcsr_set(d2ksin, 0.0_dp)
      END IF
      IF (PRESENT(d3kcos) .AND. PRESENT(d3ksin)) THEN
         CALL dbcsr_set(d3kcos, 0.0_dp)
         CALL dbcsr_set(d3ksin, 0.0_dp)
      END IF
      IF (PRESENT(d4kcos) .AND. PRESENT(d4ksin)) THEN
         CALL dbcsr_set(d4kcos, 0.0_dp)
         CALL dbcsr_set(d4ksin, 0.0_dp)
      END IF

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=ldwork, basis_type=basis_type)
      ldab = ldwork

      ! derivative size
      IF (PRESENT(nderivative)) THEN
         nder = nderivative
      ELSE
         nder = 0
      END IF
      maxder = nder + 1 

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
      k = 0
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         k = k + 1
!         WRITE(*,*) "NEIGHBORLIST ITERATION", k
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
   
         ! block deriv 1
         IF(PRESENT(dkcos) .AND. PRESENT(dksin)) THEN
            NULLIFY (ckblk)
            CALL dbcsr_get_block_p(matrix=dkcos, &
                                   row=irow, col=icol, BLOCK=ckblk, found=found)
            NULLIFY (skblk)
            CALL dbcsr_get_block_p(matrix=dksin, &
                                   row=irow, col=icol, BLOCK=skblk, found=found)
         END IF         

         ! block deriv 2
         IF(PRESENT(d2kcos) .AND. PRESENT(d2ksin)) THEN
            NULLIFY (c2kblk)
            CALL dbcsr_get_block_p(matrix=d2kcos, &
                                   row=irow, col=icol, BLOCK=c2kblk, found=found)
            NULLIFY (s2kblk)
            CALL dbcsr_get_block_p(matrix=d2ksin, &
                                   row=irow, col=icol, BLOCK=s2kblk, found=found)
         END IF         
         ! block deriv 3
         IF(PRESENT(d3kcos) .AND. PRESENT(d3ksin)) THEN
            NULLIFY (c3kblk)
            CALL dbcsr_get_block_p(matrix=d3kcos, &
                                   row=irow, col=icol, BLOCK=c3kblk, found=found)
            NULLIFY (s3kblk)
            CALL dbcsr_get_block_p(matrix=d3ksin, &
                                   row=irow, col=icol, BLOCK=s3kblk, found=found)
         END IF         
         ! block deriv 4
         IF(PRESENT(d4kcos) .AND. PRESENT(d4ksin)) THEN
            NULLIFY (c4kblk)
            CALL dbcsr_get_block_p(matrix=d2kcos, &
                                   row=irow, col=icol, BLOCK=c4kblk, found=found)
            NULLIFY (s4kblk)
            CALL dbcsr_get_block_p(matrix=d2ksin, &
                                   row=irow, col=icol, BLOCK=s4kblk, found=found)
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

                  ! Calculate the (a|exp(ikr)|b) integrals and derivatives (a|(ikr)^nder*exp(ikr)|b)
                  ! Careful, this derivative is wrt to a factor t in front of each k (t -> 1)
                  ! to pull down the whole argument of exp(ikr)
                  SELECT CASE (nder)
                  CASE(0)
                     !WRITE(*,*) "build_berry_moment_matrix4 - CASE 0"
                    CALL cossin4(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, nder = maxder, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1))
                  CASE(1)
                     !WRITE(*,*) "build_berry_moment_matrix4 - CASE 1"
                    CALL cossin4(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, nder = maxder, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1), &
                                dksc = cosmint(:, :, 2), dkss = sinmint(:, :, 2))
                  CASE(2)
                    ! WRITE(*,*) "build_berry_moment_matrix4 - CASE 2"
                    CALL cossin4(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, nder = maxder, cosab=cosmint(:, :, 1), sinab=sinmint(:, :, 1), &
                                dksc=cosmint(:, :, 2), dkss=sinmint(:, :, 2), &
                                ddksc=cosmint(:, :, 3), ddkss=sinmint(:, :, 3))
                  CASE(3)
                    CALL cossin4(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, nder = maxder, cosab=cosmint(:, :, 1), sinab=sinmint(:, :, 1), &
                                dksc=cosmint(:, :, 2), dkss=sinmint(:, :, 2), &
                                ddksc=cosmint(:, :, 3), ddkss=sinmint(:, :, 3), &
                                d3ksc=cosmint(:, :, 4), d3kss=sinmint(:, :, 4))
                  CASE(4)
                    CALL cossin4(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, nder = maxder, cosab=cosmint(:, :, 1), sinab=sinmint(:, :, 1), &
                                dksc=cosmint(:, :, 2), dkss=sinmint(:, :, 2), &
                                ddksc=cosmint(:, :, 3), ddkss=sinmint(:, :, 3), & 
                                d3ksc=cosmint(:, :, 4), d3kss=sinmint(:, :, 4), &
                                d4ksc=cosmint(:, :, 5), d4kss=sinmint(:, :, 5))
        
                  CASE DEFAULT
                     CPABORT("Derivatives of (a|(ikr)^nder*exp(ikr)|b) only up to second order")
                  END SELECT
                 
                  ! Contraction from cartesian to spherical basis functions 
                  CALL contract_cossin(cblock, sblock, &
                                       iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                       jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                       cosmint(:, :, 1), sinmint(:, :, 1), ldab, work, ldwork)
                  IF (nder >= 1) THEN
                     ! Contraction of matrices storing first derivatives
                     CALL contract_cossin(ckblk, skblk, &
                                          iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                          jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                          cosmint(:, :, 2), sinmint(:, :, 2), ldab, work, ldwork)
                  END IF
                  IF (nder >= 2) THEN
                     ! Contraction of matrices storing second derivatives
                     CALL contract_cossin(c2kblk, s2kblk, &
                                          iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                          jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                          cosmint(:, :, 3), sinmint(:, :, 3), ldab, work, ldwork)
                  END IF
                  IF (nder >= 3) THEN
                     ! Contraction of matrices storing second derivatives
                     CALL contract_cossin(c3kblk, s3kblk, &
                                          iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                          jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                          cosmint(:, :, 4), sinmint(:, :, 4), ldab, work, ldwork)
                  END IF
                  IF (nder >= 4) THEN
                     ! Contraction of matrices storing second derivatives
                     CALL contract_cossin(c4kblk, s4kblk, &
                                          iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                          jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                          cosmint(:, :, 5), sinmint(:, :, 5), ldab, work, ldwork)
                  END IF

               END DO
            END DO

         ENDIF

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (cosmint)
      DEALLOCATE (sinmint)
      DEALLOCATE (work)
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE build_berry_moment_matrix4


! **************************************************************************************************
!> \brief Calculation of the moment matrix (a|exp(ikr)|b) 
!>        and their derivatives wrt a factor t in front of each k 
!>        to pull down the whole argument of exp(ikr) up to order nderivative
!> \param qs_env The Quickstep environment
!> \param cosmat Real part of (a|exp(ikr)|b) 
!> \param sinmat Imaginary part of (a|exp(ikr)|b) 
!> \param kvec The k-vector for which the integral is evaluated
!> \param reference ...
!> \param ref_point ...
!> \param sab_orb_external Specify an external neighborlist
!> \param basis_type ...
!> \param nderivative Order of derivative wrt t
!> \param dkcos dksin Real and imaginary part of first derivative wrt t
!> \param d2kcos d2ksin Real and imaginary part of second derivative wrt t
! **************************************************************************************************
   SUBROUTINE build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec, reference, ref_point,&
                                         sab_orb_external, basis_type, nderivative, &
                                         dkcos, dksin, d2kcos, d2ksin)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kvec
      INTEGER, INTENT(IN), OPTIONAL                      :: reference
      REAL(dp), DIMENSION(:), POINTER, OPTIONAL          :: ref_point
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_external
      CHARACTER(len=*), OPTIONAL                         :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: nderivative
      TYPE(dbcsr_type), OPTIONAL, POINTER                :: dkcos, dksin, d2kcos, d2ksin 

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
      INTEGER                                            :: maxder, nder, k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: cosmint, sinmint
      REAL(dp), DIMENSION(:, :), POINTER                 :: ckblk, skblk, c2kblk, s2kblk
      REAL(dp), DIMENSION(3)                             :: rcc
!      INTEGER                                            :: i, j 

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
          CALL dbcsr_set(dkcos, 0.0_dp)
          CALL dbcsr_set(dksin, 0.0_dp)
      END IF
      
      IF (PRESENT(d2kcos) .AND. PRESENT(d2ksin)) THEN
         CALL dbcsr_set(d2kcos, 0.0_dp)
         CALL dbcsr_set(d2ksin, 0.0_dp)
      END IF

      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, maxco=ldwork, basis_type=basis_type)
      ldab = ldwork

      ! derivative size
      IF (PRESENT(nderivative)) THEN
         nder = nderivative
      ELSE
         nder = 0
      END IF
      maxder = nder + 1 

      ! Reference point
      IF (PRESENT(reference)) THEN
         CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
         rcc = pbc(rcc, cell)
      ELSE
         rcc = 0.0_dp 
      END IF
      
      ! Allocate array for (a|eikr|b) intrgral and its derivatives
      ALLOCATE(cosmint(ldab, ldab, maxder))
      ALLOCATE(sinmint(ldab, ldab, maxder))
      ALLOCATE (work(ldwork, ldwork))

      WRITE(*,*) "ldab", ldab
      WRITE(*,*) "SIZE(work, 1)", SIZE(work, 1)
      WRITE(*,*) "SIZE(work, 2)", SIZE(work, 2)

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
      k = 0
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         k = k + 1
         WRITE(*,*) "----------------- NEIGHBORLIST ITERATION ------------------", k
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

         WRITE(*,*) "BASIS_SET A"
         WRITE(*,*) "first_sgfa",  basis_set_a%first_sgf
         WRITE(*,*) "la_max",  basis_set_a%lmax
         WRITE(*,*) "la_min",  basis_set_a%lmin
         WRITE(*,*) "npgfa",  basis_set_a%npgf
         WRITE(*,*) "nseta",  basis_set_a%nset
         WRITE(*,*) "nsgfa", basis_set_a%nsgf_set
         WRITE(*,*) "rpgfa", basis_set_a%pgf_radius
         WRITE(*,*) "set_radius_a",  basis_set_a%set_radius
         WRITE(*,*) "zeta", basis_set_a%zet
         WRITE(*,*) ""
         WRITE(*,*) "BASIS_SET B"
         WRITE(*,*) "first_sgfb",  basis_set_b%first_sgf
         WRITE(*,*) "lb_max",  basis_set_b%lmax
         WRITE(*,*) "lb_min",  basis_set_b%lmin
         WRITE(*,*) "npgfb",  basis_set_b%npgf
         WRITE(*,*) "nsetb",  basis_set_b%nset
         WRITE(*,*) "nsgfb", basis_set_b%nsgf_set
         WRITE(*,*) "rpgfb", basis_set_b%pgf_radius
         WRITE(*,*) "set_radius_b",  basis_set_b%set_radius
         WRITE(*,*) "zetb", basis_set_b%zet
         WRITE(*,*) ""

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
            NULLIFY (ckblk)
            CALL dbcsr_get_block_p(matrix=dkcos, &
                                   row=irow, col=icol, BLOCK=ckblk, found=found)
            NULLIFY (skblk)
            CALL dbcsr_get_block_p(matrix=dksin, &
                                   row=irow, col=icol, BLOCK=skblk, found=found)
         END IF         

         IF(PRESENT(d2kcos) .AND. PRESENT(d2ksin)) THEN
            NULLIFY (c2kblk)
            CALL dbcsr_get_block_p(matrix=d2kcos, &
                                   row=irow, col=icol, BLOCK=c2kblk, found=found)
            NULLIFY (s2kblk)
            CALL dbcsr_get_block_p(matrix=d2ksin, &
                                   row=irow, col=icol, BLOCK=s2kblk, found=found)
         END IF         

         IF (ASSOCIATED(cblock) .AND. .NOT. ASSOCIATED(sblock) .OR. &
             .NOT. ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN
            CPABORT("")
         ENDIF

         IF (ASSOCIATED(cblock) .AND. ASSOCIATED(sblock)) THEN

            ra(:) = pbc(particle_set(iatom)%r(:), cell)
            rb(:) = ra+rab
            ! apply reference point
            ! here calculated -r_acc
            ! because r_acc = r_cc - r_a
            ! but routines are to evaluate r_a
            ! Additional multipliction with (-1.0_dp) later.
            IF(PRESENT(reference)) THEN
               ra(:) = ra-rcc
               rb(:) = rb-rcc
            END IF
            dab = SQRT(rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3))

            WRITE(*,*) "START SET ITERATION"
            DO iset = 1, nseta

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(1, iset)

               WRITE(*,*) "iset", iset
!               WRITE(*,*) "ncoa", ncoa
!               WRITE(*,*) "sgfa", sgfa
!               WRITE(*,*) ""

               DO jset = 1, nsetb

                  IF (set_radius_a(iset)+set_radius_b(jset) < dab) CYCLE

                  ncob = npgfb(jset)*ncoset(lb_max(jset))
                  sgfb = first_sgfb(1, jset)

                  WRITE(*,*) "jset", jset
!                  WRITE(*,*) "ncob", ncob
!                  WRITE(*,*) "sgfb", sgfb
!                  WRITE(*,*) ""

                  ! Calculate the (a|exp(ikr)|b) integrals and derivatives (a|(ikr)^nder*exp(ikr)|b)
                  ! Careful, this derivative is wrt to a factor t in front of each k (t -> 1)
                  ! to pull down the whole argument of exp(ikr)
                  SELECT CASE (nder)
                  CASE(0) ! up to 0-th derivative
                    CALL cossin3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1))
                  CASE(1) ! up to first derivative 
                    CALL cossin3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1), &
                                dksc = cosmint(:, :, 2), dkss = sinmint(:, :, 2))

!                    WRITE(*,*) "AFTER COSSIN - CASE 1"
!                    WRITE(*,*) "SIZE(cosmint(:, :, 1),1)", SIZE(cosmint(:, :, 1),1)
!                    WRITE(*,*) "SIZE(cosmint(:, :, 1),2)", SIZE(cosmint(:, :, 1),2)
!                    WRITE(*,*) "SIZE(cosmint(:, :, 2),1)", SIZE(cosmint(:, :, 2),1)
!                    WRITE(*,*) "SIZE(cosmint(:, :, 2),2)", SIZE(cosmint(:, :, 2),2)
!                    
!                    WRITE(*,*) "COS"
!                    DO i = 1, SIZE(cosmint(:,:,1),1)
!!                       DO j = 1, SIZE(cosmint(:,:,1),2)
!!                          WRITE(*,*) "cosmint(i,j,1)", i, j, cosmint(i,j,1)
!!                          WRITE(*,*) "cosmint(i,i,1)", i, i, cosmint(i,i,1)
!!                       END DO 
!                    END DO
!                    WRITE(*,*) "DKCOS"
!                    DO i = 1, SIZE(cosmint(:,:,2),1)
!!                       DO j = 1, SIZE(cosmint(:,:,2),2)
!!                          WRITE(*,*) "cosmint(i,j,2)", i, j, cosmint(i,j,2)
!!                          WRITE(*,*) "cosmint(i,i,2)", i, i, cosmint(i,i,2)
!!                       END DO 
!                    END DO
!                    IF(ldab >= 2) THEN
!                       WRITE(*,*) "cosmint(2, 1, 1)", cosmint(2, 1, 1) 
!                       WRITE(*,*) "cosmint(1, 2, 1)", cosmint(1, 2, 1) 
!                       WRITE(*,*) "cosmint(2, 2, 1)", cosmint(2, 2, 1) 
!                    END IF
!                    WRITE(*,*) "cosmint(1, 1, 2:4)", cosmint(1, 1, 2) 
!                    WRITE(*,*) "sinmint(1, 1, 2:4)", sinmint(1, 1, 2) 
!                    IF(ldab >= 2) THEN
!                      WRITE(*,*) "cosmint(2, 1, 2:4)", cosmint(2, 1, 2:4) 
!                       WRITE(*,*) "cosmint(1, 2, 2:4)", cosmint(1, 2, 2:4) 
!                       WRITE(*,*) "cosmint(2, 2, 2:4)", cosmint(2, 2, 2:4) 
!                    END IF

                  CASE(2) ! up to 2nd derivative
                    CALL cossin3(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab=cosmint(:, :, 1), sinab=sinmint(:, :, 1), &
                                dksc=cosmint(:, :, 2), dkss=sinmint(:, :, 2), &
                                ddksc=cosmint(:, :, 3), ddkss=sinmint(:, :, 3))
        
                  CASE DEFAULT
                     CPABORT("Derivatives of (a|(ikr)^nder*exp(ikr)|b) only up to second order")
                  END SELECT

                  ! Contraction from cartesian to spherical basis functions 
                  CALL contract_cossin(cblock, sblock, &
                                       iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                       jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                       cosmint(1:ldab, 1:ldab, 1), sinmint(1:ldab, 1:ldab, 1), &
                                       ldab, work(1:ldwork,1:ldwork), ldwork)

!                  WRITE(*,*) "AFTER CONTRACTION - CASE 0"
!                  WRITE(*,*) "SIZE(cblock,1)", SIZE(cblock,1)
!                  WRITE(*,*) "SIZE(cblock,2)", SIZE(cblock,2)
!                  WRITE(*,*) "COS"
!                  DO i = 1, SIZE(cblock,1)
!                     DO j = 1, SIZE(cblock,2)
!!                        WRITE(*,*) "cblock(i,j)", i, j, cblock(i,j)
!                     END DO 
!                  END DO

                  IF (nder >= 1) THEN
                     ! Contraction of matrices storing first derivatives

!                     WRITE(*,*) "CONTRACTION OF DKCOS and DKSIN"
!                     WRITE(*,*) "SIZE(sphi_a,1)", SIZE(sphi_a,1)
!                     WRITE(*,*) "SIZE(sphi_a,2)", SIZE(sphi_a,2)
!                     WRITE(*,*) "SIZE(sphi_b,1)", SIZE(sphi_b,1)
!                     WRITE(*,*) "SIZE(sphi_b,2)", SIZE(sphi_b,2)

                     CALL contract_cossin(ckblk, skblk, &
                                          iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                          jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                          cosmint(1:ldab, 1:ldab, 2), sinmint(1:ldab, 1:ldab, 2),&
                                          ldab, work(1:ldwork,1:ldwork), ldwork)

!                     WRITE(*,*) "AFTER CONTRACTION - CASE 1"
!                     WRITE(*,*) "cosmint(1, 1, 2)", cosmint(1, 1, 2) 
!                     WRITE(*,*) "sinmint(1, 1, 2)", sinmint(1, 1, 2) 
!                     WRITE(*,*) "ckblk(1, 1)", ckblk(1, 1) 
!                     WRITE(*,*) "skblk(1, 1)", skblk(1, 1) 

!                     WRITE(*,*) "SIZE(cosmint,1)", SIZE(cosmint,1)
!                     WRITE(*,*) "SIZE(cosmint,2)", SIZE(cosmint,2)
!                     WRITE(*,*) "SIZE(ckblk,1)", SIZE(ckblk,1)
!                     WRITE(*,*) "SIZE(ckblk,2)", SIZE(ckblk,2)
!                     WRITE(*,*) "DKCOS"
!                     DO i = 1, SIZE(ckblk,1)
!                        DO j = 1, SIZE(ckblk,2)
!                           WRITE(*,*) "ckblk(i,j)", i, j, ckblk(i,j)
!                        END DO 
!                     END DO

                  END IF

                  IF (nder >= 2) THEN
                     ! Contraction of matrices storing second derivatives
                     CALL contract_cossin(c2kblk, s2kblk, &
                                          iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                          jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                          cosmint(:, :, 3), sinmint(:, :, 3), ldab, work, ldwork)
                  END IF

               END DO
            END DO

         ENDIF

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (cosmint)
      DEALLOCATE (sinmint)
      DEALLOCATE (work)
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

   END SUBROUTINE build_berry_moment_matrix3

! **************************************************************************************************
!> \brief ... !!!!! VERSION FOR DERIVATIVE WRT k as vector
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
   SUBROUTINE build_berry_moment_matrix2(qs_env, cosmat, sinmat, kvec, &
                                         sab_orb_external, basis_type, nderivative, &
                                         dkcos, dksin, d2kcos, d2ksin)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), POINTER                          :: cosmat, sinmat
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: kvec
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb_external
      CHARACTER(len=*), OPTIONAL                         :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: nderivative
      TYPE(dbcsr_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: dkcos, dksin, d2kcos, d2ksin 

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
            CALL dbcsr_set(dkcos(i)%matrix, 0.0_dp)
            CALL dbcsr_set(dksin(i)%matrix, 0.0_dp)
         END DO
      END IF
      
      IF (PRESENT(d2kcos) .AND. PRESENT(d2ksin)) THEN
         DO i = 1, 6
            CALL dbcsr_set(d2kcos(i)%matrix, 0.0_dp)
            CALL dbcsr_set(d2ksin(i)%matrix, 0.0_dp)
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

               CALL dbcsr_get_block_p(matrix=dkcos(1)%matrix, &
                                      row=irow, col=icol, BLOCK=ckblk1, found=found)
               CALL dbcsr_get_block_p(matrix=dkcos(2)%matrix, &
                                      row=irow, col=icol, BLOCK=ckblk2, found=found)
               CALL dbcsr_get_block_p(matrix=dkcos(3)%matrix, &
                                      row=irow, col=icol, BLOCK=ckblk3, found=found)
                                   
               CALL dbcsr_get_block_p(matrix=dksin(1)%matrix, &
                                      row=irow, col=icol, BLOCK=skblk1, found=found)
               CALL dbcsr_get_block_p(matrix=dksin(2)%matrix, &
                                      row=irow, col=icol, BLOCK=skblk2, found=found)
               CALL dbcsr_get_block_p(matrix=dksin(3)%matrix, &
                                      row=irow, col=icol, BLOCK=skblk3, found=found)
         END IF         

         IF(PRESENT(d2kcos) .AND. PRESENT(d2ksin)) THEN
!            WRITE(*,*) "d2kcos - associate blocks"
      
            ! THIS IS AWEFUL - but 'BLOCK' wants a pointer
            NULLIFY (c2kblk1,c2kblk2,c2kblk3,c2kblk4,c2kblk5,c2kblk6)
            NULLIFY (s2kblk1,s2kblk2,s2kblk3,s2kblk4,s2kblk5,s2kblk6)

               CALL dbcsr_get_block_p(matrix=d2kcos(1)%matrix, &
                                      row=irow, col=icol, BLOCK=c2kblk1, found=found)
               CALL dbcsr_get_block_p(matrix=d2kcos(2)%matrix, &
                                      row=irow, col=icol, BLOCK=c2kblk2, found=found)
               CALL dbcsr_get_block_p(matrix=d2kcos(3)%matrix, &
                                      row=irow, col=icol, BLOCK=c2kblk3, found=found)
               CALL dbcsr_get_block_p(matrix=d2kcos(4)%matrix, &
                                      row=irow, col=icol, BLOCK=c2kblk4, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(5)%matrix, &
                                      row=irow, col=icol, BLOCK=c2kblk5, found=found)
               CALL dbcsr_get_block_p(matrix=d2kcos(6)%matrix, &
                                      row=irow, col=icol, BLOCK=c2kblk6, found=found)

               CALL dbcsr_get_block_p(matrix=d2ksin(1)%matrix, &
                                      row=irow, col=icol, BLOCK=s2kblk1, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(2)%matrix, &
                                      row=irow, col=icol, BLOCK=s2kblk2, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(3)%matrix, &
                                      row=irow, col=icol, BLOCK=s2kblk3, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(4)%matrix, &
                                      row=irow, col=icol, BLOCK=s2kblk4, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(5)%matrix, &
                                      row=irow, col=icol, BLOCK=s2kblk5, found=found)
               CALL dbcsr_get_block_p(matrix=d2ksin(6)%matrix, &
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

                  ! Calculate the (a|exp(ikr)|b) integrals and derivatives (a|(ir)^nder*exp(ikr)|b)
                  ! Careful, this derivative is wrt k, NOT the ion positions
                  SELECT CASE (nder)
                  CASE(0)
                    CALL cossin2(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1))
        
                  CASE(1)
                    CALL cossin2(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab = cosmint(:, :, 1), sinab = sinmint(:, :, 1), &
                                dksc = cosmint(:, :, 2:4), dkss = sinmint(:, :, 2:4))

                  CASE(2)
                    CALL cossin2(la_max(iset), npgfa(iset), zeta(:, iset), rpgfa(:, iset), la_min(iset), &
                                lb_max(jset), npgfb(jset), zetb(:, jset), rpgfb(:, jset), lb_min(jset), &
                                ra, rb, kvec, cosab=cosmint(:, :, 1), sinab=sinmint(:, :, 1), &
                                dksc=cosmint(:, :, 2:4), dkss=sinmint(:, :, 2:4), &
                                ddksc=cosmint(:, :, 5:10), ddkss=sinmint(:, :, 5:10))
        
                  CASE DEFAULT
                     CPABORT("Derivatives of (a|exp(-ikr)|b) only up to second order")
                  END SELECT
                 
                  ! Contraction from cartesian to spherical basis functions 
                  CALL contract_cossin(cblock, sblock, &
                                       iatom, ncoa, nsgfa(iset), sgfa, sphi_a, ldsa, &
                                       jatom, ncob, nsgfb(jset), sgfb, sphi_b, ldsb, &
                                       cosmint(:, :, 1), sinmint(:, :, 1), ldab, work, ldwork)

                  IF (nder >= 1) THEN
                     ! Contraction of matrices storing first derivatives
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

      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (cosmint)
      DEALLOCATE (sinmint)
      DEALLOCATE (work)
      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

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
!> \brief ... based on cleaned phase2, to test higher order derivs
!> \param qs_env ...
!> \param magnetic ...
!> \param nmoments ...
!> \param reference ...
!> \param ref_point ...
!> \param unit_number ...
!> \param method - Method to calculate electronic dipole contribution
!> history fbelle 10/2018 - Added new method for electronic dipole, quadrupole calculation 
! **************************************************************************************************
   SUBROUTINE qs_moment_berry_phase4(qs_env, magnetic, nmoments, reference, ref_point, unit_number)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: magnetic
      INTEGER, INTENT(IN)                                :: nmoments, reference
      REAL(dp), DIMENSION(:), POINTER                    :: ref_point
      INTEGER, INTENT(IN)                                :: unit_number

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_moment_berry_phase4', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)        :: rlab
      CHARACTER(LEN=default_string_length)               :: description
      COMPLEX(dp)                                        :: xphase(3), zdeta, zi(3), &
                                                            zij(3, 3), zijk(3, 3, 3), &
                                                            zijkl(3, 3, 3, 3), zphase(3), zz
      INTEGER                                            :: handle, i, ia, ikind, ispin, ix, &
                                                            iy, iz, j, k, l, nao, nm, nmo, nmom, &
                                                            nmotot
      LOGICAL                                            :: floating, ghost, uniform
      REAL(dp)                                           :: charge, ci(3), cij(3, 3), dd, occ, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mmom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rmom
      REAL(dp), DIMENSION(3)                             :: kvec, qq, rcc, ria
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER         :: eigrmat
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mos_new, opvec
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: op_fm_set, dop_fm_set, d2op_fm_set, d3op_fm_set, d4op_fm_set
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
      !fbelle
      REAL(dp)                                           :: relem, ielem
      COMPLEX(dp)                                        :: det2, det1, et1
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dkcos, dksin

!      TYPE(dbcsr_type), POINTER                          :: dkcos, dksin, d2kcos, d2ksin, d3kcos, d3ksin, d4kcos, d4ksin
      INTEGER                                            :: nelectron
      !fbelle
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:)              :: tr1
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:)            :: tr2
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:,:)          :: tr3
      COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:)        :: tr4
      COMPLEX(dp)                                           :: ztr
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER         :: deigrmat, cumulant_matrix
      INTEGER                                            :: idim, tmp_dim

      CPASSERT(ASSOCIATED(qs_env))

      IF (ASSOCIATED(qs_env%ls_scf_env)) THEN
         IF (unit_number > 0) WRITE (unit_number, *) "Periodic moment calculation not implemented in linear scaling code"
         RETURN
      ENDIF

      CALL timeset(routineN, handle)

      ! restrict maximum moment available
      nmom = MIN(nmoments, 4)
      WRITE(*,*) "nmom", nmom

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

      ! Allocate array to store traces
!      NULLIFY(tr1, tr2, tr3, tr4)
      ALLOCATE(tr1(0:nmom,3))
      ALLOCATE(tr2(0:nmom,3,3))
      ALLOCATE(tr3(0:nmom,3,3,3))
      ALLOCATE(tr4(0:nmom,3,3,3,3))

      tr1 = 0.0_dp
      tr2 = 0.0_dp
      tr3 = 0.0_dp
      tr4 = 0.0_dp

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
                      mos=mos, &
                      nelectron_total = nelectron)

                   WRITE(*,*) "ELECTRONS", nelectron

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
         CALL dbcsr_allocate_matrix_set(dkcos, nmom)
         CALL dbcsr_allocate_matrix_set(dksin, nmom)
         DO i = 1, nmom
            ALLOCATE(dkcos(i)%matrix, dksin(i)%matrix)
            CALL dbcsr_copy(dkcos(i)%matrix, matrix_s(1)%matrix, 'DKCOS MOM')
            CALL dbcsr_copy(dksin(i)%matrix, matrix_s(1)%matrix, 'DKSIN MOM')
            CALL dbcsr_set(dkcos(i)%matrix, 0.0_dp)
            CALL dbcsr_set(dksin(i)%matrix, 0.0_dp)
         END DO
      END IF

      ALLOCATE (op_fm_set(2, dft_control%nspins))
      ALLOCATE (opvec(dft_control%nspins))
      ALLOCATE (eigrmat(dft_control%nspins))

      IF(nmom >= 1) ALLOCATE (dop_fm_set(2, dft_control%nspins))
      IF(nmom >= 2) ALLOCATE (d2op_fm_set(2, dft_control%nspins))
      IF(nmom >= 3) ALLOCATE (d3op_fm_set(2, dft_control%nspins))
      IF(nmom >= 4) ALLOCATE (d4op_fm_set(2, dft_control%nspins))

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
            
            IF(nmom >= 1) THEN
               NULLIFY (dop_fm_set(i, ispin)%matrix)
               CALL cp_fm_create(dop_fm_set(i, ispin)%matrix, tmp_fm_struct)
               ALLOCATE (deigrmat(dft_control%nspins))
               ALLOCATE (cumulant_matrix(dft_control%nspins))
            END IF
            IF(nmom >= 2) THEN
               NULLIFY (d2op_fm_set(i, ispin)%matrix)
               CALL cp_fm_create(d2op_fm_set(i, ispin)%matrix, tmp_fm_struct)
            END IF
            IF(nmom >= 3) THEN
               NULLIFY (d3op_fm_set(i, ispin)%matrix)
               CALL cp_fm_create(d3op_fm_set(i, ispin)%matrix, tmp_fm_struct)
            END IF
            IF(nmom >= 4) THEN
               NULLIFY (d4op_fm_set(i, ispin)%matrix)
               CALL cp_fm_create(d4op_fm_set(i, ispin)%matrix, tmp_fm_struct)
            END IF

         END DO
         CALL cp_cfm_create(eigrmat(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
         CALL cp_cfm_create(deigrmat(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
         CALL cp_cfm_create(cumulant_matrix(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
      END DO

      ! Initialize and allocate matrices to store derivatives
      IF(nmom >= 1) THEN
         ALLOCATE (dop_fm_set(2, dft_control%nspins))

         NULLIFY(dkcos, dksin)
         CALL dbcsr_allocate_matrix_set(dkcos, 3)
         CALL dbcsr_allocate_matrix_set(dksin, 3)
         DO i = 1, 3
            ALLOCATE(dkcos(i)%matrix, dksin(i)%matrix)
            CALL dbcsr_copy(dkcos(i)%matrix, matrix_s(1)%matrix, 'DKCOS MOM')
            CALL dbcsr_copy(dksin(i)%matrix, matrix_s(1)%matrix, 'DKSIN MOM')
            CALL dbcsr_set(dkcos(i)%matrix, 0.0_dp)
            CALL dbcsr_set(dksin(i)%matrix, 0.0_dp)
         END DO
      END IF

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
                     zij(j, i) = zij(i, j)
                  END DO
               END DO
            END DO
            DO i = 1, 3
               DO j = 1, 3
                  zij(i, j) = zij(i, j)*zphase(i)*zphase(j)
                  zz = zij(i, j)/zi(i)/zi(j)
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
            END DO
            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))
            DO k = 4, 9
               ix = indco(1, k+1)
               iy = indco(2, k+1)
               iz = indco(3, k+1)
!               IF (ix == 0) THEN
!                  rmom(k+1, 2) = cij(iy, iz)
!               ELSE IF (iy == 0) THEN
!                  rmom(k+1, 2) = cij(ix, iz)
!               ELSE IF (iz == 0) THEN
!                  rmom(k+1, 2) = cij(ix, iy)
!               END IF
         ! Hardcoding of upper triangular matrix elements of cij
         ! into moment vector
         ! whatever ... 
         rmom(5,2) = cij(1,1) 
         rmom(6,2) = cij(1,2) 
         rmom(7,2) = cij(1,3) 
         rmom(8,2) = cij(2,2) 
         rmom(9,2) = cij(2,3) 
         rmom(10,2) = cij(3,3) 
            END DO
         CASE (3)
            ! Octapole
!            CPABORT("Berry phase moments bigger than 2 not implemented")
         CASE (4)
            ! Hexadecapole
!            CPABORT("Berry phase moments bigger than 3 not implemented")
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

      ztr = 0.0_dp

      IF (nmom >=1 ) THEN
         ! Dipole
         DO i = 1, 3
            WRITE(*,*) "DIPOLE DIRECTION", i
            kvec(:) = twopi*cell%h_inv(i, :)

            ! d/dt (a|e^(ikr)|b)
            CALL build_berry_moment_matrix4(qs_env, cosmat, sinmat, kvec, nderivative = 1, &
                                            dkcos = dkcos(1)%matrix, dksin = dksin(1)%matrix)

            CALL dbcsr_print(cosmat, matlab_format = .TRUE.)
            CALL dbcsr_print(sinmat, matlab_format = .TRUE.)
            CALL dbcsr_print(dkcos(1)%matrix, matlab_format = .TRUE.)
            CALL dbcsr_print(dksin(1)%matrix, matlab_format = .TRUE.)

            ! rotation from AOs onto MOs               
            IF (qs_env%run_rtp) THEN
               CALL get_qs_env(qs_env, rtp=rtp)
               CALL get_rtp(rtp, mos_new=mos_new)
               CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
               CALL op_orbbas_rtp(dkcos(1)%matrix, dksin(1)%matrix, mos, dop_fm_set, mos_new)
            ELSE
               CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
               CALL op_orbbas(dkcos(1)%matrix, dksin(1)%matrix, mos, dop_fm_set, opvec)
            END IF

            ! Electronic dipole moment
            ! Elements :
            ! d11 = (a|(-ikr)*e^(-ikr))|b)
            ! d01 = (a|(e^(-ikr))|b)

            ! get traces
            DO ispin = 1, dft_control%nspins

!               ! Trace deriv real
!               CALL copy_fm_to_dbcsr(dop_fm_set(1, ispin)%matrix, dkcos(1)%matrix)
!               CALL dbcsr_trace(dkcos(1)%matrix, relem)
!               ! Trace deriv img
!               CALL cp_fm_scale(-1.0_dp, dop_fm_set(2, ispin)%matrix)
!               CALL copy_fm_to_dbcsr(dop_fm_set(2, ispin)%matrix, dksin(1)%matrix)
!               CALL dbcsr_trace(dksin(1)%matrix, ielem)
!               det1 = CMPLX(relem, ielem, dp)
!               WRITE(*,*) "det1", det1
!
!               ! real
!               CALL copy_fm_to_dbcsr(op_fm_set(1, ispin)%matrix, cosmat)
!               CALL dbcsr_trace(cosmat, relem)
!               ! img
!               CALL cp_fm_scale(-1.0_dp, op_fm_set(2, ispin)%matrix)
!               CALL copy_fm_to_dbcsr(op_fm_set(2, ispin)%matrix, sinmat)
!               CALL dbcsr_trace(sinmat, ielem)
!
!               ! Tr(d11) / Tr(d01)
!               et1 = CMPLX(relem, ielem, dp)
!               WRITE(*,*) "et1", et1

                CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=tmp_dim)
                DO idim = 1, tmp_dim

                   ! Make complex matrix of d01
                   eigrmat(ispin)%matrix%local_data(:, idim) = &
                      CMPLX(op_fm_set(1, ispin)%matrix%local_data(:, idim), &
                            -op_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
                   ! Make complex matrix of d11
                   deigrmat(ispin)%matrix%local_data(:, idim) = &
                      CMPLX(dop_fm_set(1, ispin)%matrix%local_data(:, idim), &
                            -dop_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
                END DO
                
                ! Invert d01
                CALL cp_cfm_lu_invert(eigrmat(ispin)%matrix)

!>              matrix_c = alpha * op1( matrix_a ) * op2( matrix_b ) + beta*matrix_c.
                CALL cp_cfm_gemm("N","N", idim, idim, idim, &
                                 z_one, eigrmat(ispin)%matrix, deigrmat(ispin)%matrix, &
                                 z_zero, cumulant_matrix(ispin)%matrix)

                ! cfm -> fm
                ! nmo x nmo 
                CALL cp_cfm_to_fm(msource = cumulant_matrix(ispin)%matrix, &
                                  mtargetr = op_fm_set(1, ispin)%matrix, &
                                  mtargeti = op_fm_set(2, ispin)%matrix)

                WRITE(*,*) "nao_op_fm_set - COS - AFTER A^-1 * A"
                CALL cp_fm_write_formatted(dop_fm_set(1, ispin)%matrix, 6)
                WRITE(*,*) "nao_op_fm_set - SIN - AFTER A^-1 * A"
                CALL cp_fm_write_formatted(dop_fm_set(2, ispin)%matrix, 6)

                ! fm -> dbcsr
                CALL copy_fm_to_dbcsr(op_fm_set(1, ispin)%matrix, cosmat)
                CALL copy_fm_to_dbcsr(op_fm_set(2, ispin)%matrix, sinmat)

               ! calculate trace
               CALL dbcsr_trace(cosmat, relem)
               CALL dbcsr_trace(sinmat, ielem)
               et1 =  CMPLX(relem, ielem, dp)

               ! Sum traces for spin
               !tr1(0, i) = tr1(0, i) + et1
               !tr1(1, i) = tr1(1, i) + det1
            END DO

            ! Divide by (-ik)
            !zi(i) = tr1(1, i) / tr1(0, i) / CMPLX(0.0_dp, 1.0_dp*kvec(i), dp)
            zi(i) = et1
            WRITE(*,*) "zi", zi(i)

         END DO
         zi = zi*xphase
         WRITE(*,*) "TOTAL zi", zi

      END IF

      IF (nmom >=2) THEN
         ! Quadrupole
!         CPABORT("Berry phase moments bigger than 1 not implemented")
         DO i = 1, 3
            WRITE(*,*) "QUADRUPOLE"
            WRITE(*,*) "QUADRUPOLE DIRECTION i", i
            kvec(:) = twopi*cell%h_inv(i, :)

            ! d/dt (a|e^(ikr)|b)
            ! For elements d01(i), d11(i), d21(i)
            ! For elements d01(j), d11(j), d21(j)
            CALL build_berry_moment_matrix4(qs_env, cosmat, sinmat, kvec, nderivative = 2, &
                                            dkcos = dkcos(1)%matrix, dksin = dksin(1)%matrix, &
                                            d2kcos = dkcos(2)%matrix, d2ksin = dksin(2)%matrix)

            IF (qs_env%run_rtp) THEN
               CALL get_qs_env(qs_env, rtp=rtp)
               CALL get_rtp(rtp, mos_new=mos_new)
               CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
               CALL op_orbbas_rtp(dkcos(1)%matrix, dksin(1)%matrix, mos, dop_fm_set, mos_new)
               CALL op_orbbas_rtp(dkcos(2)%matrix, dksin(2)%matrix, mos, dop_fm_set, mos_new)
            ELSE
               CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
               CALL op_orbbas(dkcos(1)%matrix, dksin(1)%matrix, mos, dop_fm_set(:,:), opvec)
               CALL op_orbbas(dkcos(2)%matrix, dksin(2)%matrix, mos, d2op_fm_set(:,:), opvec)
            END IF

            ! get traces
            DO ispin = 1, dft_control%nspins

               ! Trace 2nd deriv real
               CALL copy_fm_to_dbcsr(d2op_fm_set(1, ispin)%matrix, dkcos(2)%matrix)
               CALL dbcsr_trace(dkcos(2)%matrix, relem)
               ! Trace 2nd deriv img
               CALL cp_fm_scale(-1.0_dp, d2op_fm_set(2, ispin)%matrix)
               CALL copy_fm_to_dbcsr(d2op_fm_set(2, ispin)%matrix, dksin(2)%matrix)
               CALL dbcsr_trace(dksin(2)%matrix, ielem)
               det2 = CMPLX(relem, ielem, dp)

               ! Trace deriv real
               CALL copy_fm_to_dbcsr(dop_fm_set(1, ispin)%matrix, dkcos(1)%matrix)
               CALL dbcsr_trace(dkcos(1)%matrix, relem)
               ! Trace deriv img
               CALL cp_fm_scale(-1.0_dp, dop_fm_set(2, ispin)%matrix)
               CALL copy_fm_to_dbcsr(dop_fm_set(2, ispin)%matrix, dksin(1)%matrix)
               CALL dbcsr_trace(dksin(1)%matrix, ielem)
               det1 = CMPLX(relem, ielem, dp)
               WRITE(*,*) "det1", det1

               ! real
               CALL copy_fm_to_dbcsr(op_fm_set(1, ispin)%matrix, cosmat)
               CALL dbcsr_trace(cosmat, relem)
               ! img
               CALL cp_fm_scale(-1.0_dp, op_fm_set(2, ispin)%matrix)
               CALL copy_fm_to_dbcsr(op_fm_set(2, ispin)%matrix, sinmat)
               CALL dbcsr_trace(sinmat, ielem)
               et1 = CMPLX(relem, ielem, dp)
               WRITE(*,*) "et1", et1

               ! Sum traces for spin
               tr1(0, i) = tr1(0, i) + et1
               tr1(1, i) = tr1(1, i) + det1
               tr1(2, i) = tr1(2, i) + det2
            END DO

            DO j = i, 3
               WRITE(*,*) "QUADRUPOLE DIRECTION j", j
               kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :))

               CALL build_berry_moment_matrix4(qs_env, cosmat, sinmat, kvec, nderivative = 2, &
                                               dkcos = dkcos(1)%matrix, dksin = dksin(1)%matrix, &
                                               d2kcos = dkcos(2)%matrix, d2ksin = dksin(2)%matrix)

               IF (qs_env%run_rtp) THEN
                  CALL get_qs_env(qs_env, rtp=rtp)
                  CALL get_rtp(rtp, mos_new=mos_new)
                  CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
                  CALL op_orbbas_rtp(dkcos(1)%matrix, dksin(1)%matrix, mos, dop_fm_set, mos_new)
                  CALL op_orbbas_rtp(dkcos(2)%matrix, dksin(2)%matrix, mos, dop_fm_set, mos_new)
               ELSE
                  CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
                  CALL op_orbbas(dkcos(1)%matrix, dksin(1)%matrix, mos, dop_fm_set(:,:), opvec)
                  CALL op_orbbas(dkcos(2)%matrix, dksin(2)%matrix, mos, d2op_fm_set(:,:), opvec)
               END IF

               ! get traces
               DO ispin = 1, dft_control%nspins

                  ! Trace 2nd deriv real
                  CALL copy_fm_to_dbcsr(d2op_fm_set(1, ispin)%matrix, dkcos(2)%matrix)
                  CALL dbcsr_trace(dkcos(2)%matrix, relem)
                  ! Trace 2nd deriv img
                  CALL cp_fm_scale(-1.0_dp, d2op_fm_set(2, ispin)%matrix)
                  CALL copy_fm_to_dbcsr(d2op_fm_set(2, ispin)%matrix, dksin(2)%matrix)
                  CALL dbcsr_trace(dksin(2)%matrix, ielem)
                  det2 = CMPLX(relem, ielem, dp)
                  WRITE(*,*) "det2", det2

                  ! Trace first deriv real
                  CALL copy_fm_to_dbcsr(dop_fm_set(1, ispin)%matrix, dkcos(1)%matrix)
                  CALL dbcsr_trace(dkcos(1)%matrix, relem)
                  ! Trace first deriv img
                  CALL cp_fm_scale(-1.0_dp, dop_fm_set(2, ispin)%matrix)
                  CALL copy_fm_to_dbcsr(dop_fm_set(2, ispin)%matrix, dksin(1)%matrix)
                  CALL dbcsr_trace(dksin(1)%matrix, ielem)
                  det1 = CMPLX(relem, ielem, dp)
                  WRITE(*,*) "det1", det1

                  ! Real  
                  CALL copy_fm_to_dbcsr(op_fm_set(1, ispin)%matrix, cosmat)
                  CALL dbcsr_trace(cosmat, relem)
                  ! img
                  CALL cp_fm_scale(-1.0_dp, op_fm_set(2, ispin)%matrix)
                  CALL copy_fm_to_dbcsr(op_fm_set(2, ispin)%matrix, sinmat)
                  CALL dbcsr_trace(sinmat, ielem)
                  et1 = CMPLX(relem, ielem, dp)
                  WRITE(*,*) "et1", et1

                  ! Sum traces for spin
                  tr2(0, i, j) = tr2(0, i, j) + et1
                  tr2(1, i, j) = tr2(1, i, j) + det1
                  tr2(2, i ,j) = tr2(2, i, j) + det2

               END DO

               ztr = 0.5*(tr2(2, i, j) / tr2(0, i, j) - & 
                     tr2(1, i, j)**2 /tr2(0, i, j)**2 - &
                     tr1(2, i) / tr1(0, i) + &
                     tr1(1, i)**2 / tr1(0, i)**2 - &
                     tr1(2, j) / tr1(0, j) + &
                     tr1(1, j)**2 / tr1(0, j)**2)
                 
               WRITE(*,*) "ztr", ztr           

               zij(i, j) = ztr*xphase(i)*xphase(j)
               zij(j, i) = ztr*xphase(i)*xphase(j)

            END DO
         END DO

      END IF

      IF (nmom >= 3) THEN
            ! Octapole
            CPABORT("Berry phase moments bigger than 2 not implemented")
!            DO i = 1, 3
!               DO j = i, 3
!                  DO k = j, 3
!                     kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :)+cell%h_inv(k, :) )
!
!                  END DO
!               END DO
!            END DO
!
      END IF

      IF (nmom >= 4) THEN
            ! Hexadecapole
            CPABORT("Berry phase moments bigger than 3 not implemented")
!            DO i = 1, 3
!               DO j = i, 3
!                  DO k = j, 3
!                     DO l = k, 3
!                        kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :) &
!                                        +cell%h_inv(k, :)+cell%h_inv(l, :))
!
!                     END DO
!                  END DO
!               END DO
!            END DO
      END IF
 
      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)

            ! Dipole (apply periodic (2 Pi) boundary conditions)
               !ci = AIMAG(LOG(zi))
               !DO i = 1, 3
               !   IF (qq(i)+ci(i) > pi) ci(i) = ci(i)-twopi
               !   IF (qq(i)+ci(i) < -pi) ci(i) = ci(i)+twopi
               !END DO
               !rmom(2:4, 1) = MATMUL(cell%hmat, ci)/twopi
               
               ! Convert from cumulants to moments
               ci = REAL(zi)
               DO i = 1, 3
                  IF (qq(i)+ci(i) > pi) ci(i) = ci(i)-twopi
                  IF (qq(i)+ci(i) < -pi) ci(i) = ci(i)+twopi
               END DO
               
               !rmom(2:4, 1) = ci*nelectron
               rmom(2:4, 1) = 2.0_dp*MATMUL(cell%hmat, ci)/twopi 

            WRITE(*,*) ""
            WRITE(*,*) "nuclear", rmom(2:4, 2)
            WRITE(*,*) "rmom", rmom(2:4, 1)
            WRITE(*,*) ""

         CASE (2)
            ! Quadrupole
!            CPABORT("Berry phase moments bigger than 1 not implemented")
!            DO i = 1, 3
!               DO j = 1, 3
!                  zz = zij(i, j)/zi(i)/zi(j)
!                  cij(i, j) = AIMAG(LOG(zz))/twopi
!               END DO
!            END DO
!            cij = 0.5_dp*cij/twopi/twopi
!            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))

            ! Convert from cumulants to moments
            DO i = 1, 3
               DO j = 1, 3
                  cij(i, j) = REAL(zij(i, j)) + REAL(ci(i))*REAL(ci(j))
               END DO
            END DO
     
            ! Hardcoding of upper triangular matrix elements of cij
            ! into moment vector
            ! whatever ... 
            rmom(5,1) = cij(1,1)*nelectron 
            rmom(6,1) = cij(1,2)*nelectron 
            rmom(7,1) = cij(1,3)*nelectron 
            rmom(8,1) = cij(2,2)*nelectron 
            rmom(9,1) = cij(2,3)*nelectron 
            rmom(10,1) = cij(3,3)*nelectron 

            WRITE(*,*) "XX", rmom(5,1)
            WRITE(*,*) "XY", rmom(6,1)
            WRITE(*,*) "XZ", rmom(7,1)
            WRITE(*,*) "YY", rmom(8,1)
            WRITE(*,*) "ZY", rmom(9,1)
            WRITE(*,*) "ZZ", rmom(10,1)

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

      IF(nmom >= 1) THEN
         CALL dbcsr_deallocate_matrix_set(dkcos)
         CALL dbcsr_deallocate_matrix_set(dksin)
      END IF
     
      DO ispin = 1, dft_control%nspins
         CALL cp_fm_release(opvec(ispin)%matrix)
         CALL cp_cfm_release(eigrmat(ispin)%matrix)
         DO i = 1, SIZE(op_fm_set, 1)
            CALL cp_fm_release(op_fm_set(i, ispin)%matrix)
            IF(nmom >= 1) CALL cp_fm_release(dop_fm_set(i, ispin)%matrix)
            IF(nmom >= 2) CALL cp_fm_release(d2op_fm_set(i, ispin)%matrix)
            IF(nmom >= 3) CALL cp_fm_release(d3op_fm_set(i, ispin)%matrix)
            IF(nmom >= 4) CALL cp_fm_release(d4op_fm_set(i, ispin)%matrix)
         END DO
      END DO
      DEALLOCATE (op_fm_set)

      IF(nmom >= 1) DEALLOCATE (dop_fm_set)
      IF(nmom >= 2) DEALLOCATE (d2op_fm_set)
      IF(nmom >= 3) DEALLOCATE (d3op_fm_set)
      IF(nmom >= 4) DEALLOCATE (d4op_fm_set)

      DEALLOCATE (opvec)
      DEALLOCATE (eigrmat)
      DEALLOCATE (tr1)
      DEALLOCATE (tr2)
      DEALLOCATE (tr3)
      DEALLOCATE (tr4)

      DO ispin = 1, dft_control%nspins
         CALL cp_cfm_release(deigrmat(ispin)%matrix)
         CALL cp_cfm_release(cumulant_matrix(ispin)%matrix)
      END DO

      DEALLOCATE (deigrmat)
      DEALLOCATE (cumulant_matrix)

      CALL timestop(handle)

   END SUBROUTINE qs_moment_berry_phase4

! **************************************************************************************************
!> \brief ...                  THE ONE WITH RATIO BEFORE TRANSFORMATION
!> \param qs_env ...
!> \param magnetic ...
!> \param nmoments ...
!> \param reference ...
!> \param ref_point ...
!> \param unit_number ...
! **************************************************************************************************
   SUBROUTINE qs_moment_berry_phase3(qs_env, magnetic, nmoments, reference, ref_point, unit_number)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: magnetic
      INTEGER, INTENT(IN)                                :: nmoments, reference
      REAL(dp), DIMENSION(:), POINTER                    :: ref_point
      INTEGER, INTENT(IN)                                :: unit_number

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_moment_berry_phase3', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:)        :: rlab
      CHARACTER(LEN=default_string_length)               :: description
      COMPLEX(dp)                                        :: xphase(3), zdet, zdeta, zi(3), &
                                                            zij(3, 3), zijk(3, 3, 3), &
                                                            zijkl(3, 3, 3, 3), zphase(3), zz
      INTEGER                                            :: handle, i, ia, idim, ikind, ispin, ix, &
                                                            iy, iz, j, k, l, nao, nm, nmo, nmom, &
                                                            nmotot
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

      ! fbelle
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: tmp_op_fm_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dkcos, dksin
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: dop_fm_set

!      REAL(dp)                                           :: newci(3)
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: newrmom
      INTEGER                                            :: tmp_dim

      ! fbelle - inversion craze
      REAL(dp)                                           :: tracei(3)
!      REAL(dp)                                           :: strace
!      REAL(dp)                                           :: tracei(3), sintrace, strace
      REAL(dp)                                           :: relem, ielem

      ! nao x nao
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: nao_dop_fm_set, nao_op_fm_set
      TYPE(cp_fm_struct_type), POINTER                   :: nao_fm_struct
      TYPE(cp_cfm_p_type), POINTER                       :: eeigrmat
      TYPE(cp_cfm_p_type), POINTER                       :: deigrmat
      TYPE(cp_cfm_p_type), POINTER                       :: deriveigrmat
      TYPE(cp_cfm_p_type), POINTER                       :: inveigrmat
      TYPE(cp_cfm_p_type), POINTER                       :: cumulant_matrix

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: op_cos_dbcsr, op_sin_dbcsr
      INTEGER                                            :: icol, irow
      INTEGER                                            :: ncol_local, nrow_local
      REAL(KIND=dp), DIMENSION(:), POINTER               :: diag
!      REAL (dp)                                          :: ctrace, sstrace, cdtrace, sdtrace
      COMPLEX(dp)                                        :: ctrace(3)
      COMPLEX(dp)                                        :: etr, et1
      ! new nuclear
      INTEGER                                            :: nelectron
      REAL(dp)                                           :: qn
      COMPLEX(dp), DIMENSION(3)                          :: d01, d11
      REAL(dp), DIMENSION(3)                             :: newqq

      CPASSERT(ASSOCIATED(qs_env))

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
      ALLOCATE (newrmom(nm+1, 3))
      newrmom = 0.0_dp
      rlab = ""
      IF (magnetic) THEN
         nm = 3
         ALLOCATE (mmom(nm))
         mmom = 0._dp
      END IF

     ! fbelle
     ! DEBUGGING OF DERIVATIVES
!     CALL debug_berry_moment_matrix(qs_env, nmoments = 1)

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
                      mos=mos, &
                      nelectron_total = nelectron)

      IF(para_env%mepos==0) THEN
         WRITE(*,*) "qs_moment_berry_phase3 - BEGINNING"
      END IF
      CALL qs_rho_get(rho, rho_ao=rho_ao)

!----------------
! DBCSR NAO x NAO
!----------------

      NULLIFY (cosmat, sinmat)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL dbcsr_set(cosmat, 0.0_dp)
      CALL dbcsr_set(sinmat, 0.0_dp)

      ! Initialize and allocate matrices to store derivatives
      IF(nmom >= 1) THEN
         NULLIFY(dkcos, dksin)
         CALL dbcsr_allocate_matrix_set(dkcos, 3)
         CALL dbcsr_allocate_matrix_set(dksin, 3)
         DO i = 1, 3
            ALLOCATE(dkcos(i)%matrix, dksin(i)%matrix)
            CALL dbcsr_copy(dkcos(i)%matrix, matrix_s(1)%matrix, 'DKCOS MOM')
            CALL dbcsr_copy(dksin(i)%matrix, matrix_s(1)%matrix, 'DKSIN MOM')
            CALL dbcsr_set(dkcos(i)%matrix, 0.0_dp)
            CALL dbcsr_set(dksin(i)%matrix, 0.0_dp)
         END DO
      END IF
      
!----------------

!----------------
! FM NMO x NMO
!----------------

      ALLOCATE (op_fm_set(2, dft_control%nspins))
      ALLOCATE (tmp_op_fm_set(2, dft_control%nspins))
      ALLOCATE (opvec(dft_control%nspins))
      ALLOCATE (eigrmat(dft_control%nspins))
      
      ! Initialize and allocate matrices to store derivatives
      IF(nmom >= 1) THEN
         ALLOCATE (dop_fm_set(2, dft_control%nspins))
      END IF

      ! fbelle - Store projected derivatives 
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
            NULLIFY (dop_fm_set(i, ispin)%matrix)
            NULLIFY (tmp_op_fm_set(i, ispin)%matrix)
            CALL cp_fm_create(op_fm_set(i, ispin)%matrix, tmp_fm_struct)
            CALL cp_fm_create(dop_fm_set(i, ispin)%matrix, tmp_fm_struct)
            CALL cp_fm_create(tmp_op_fm_set(i, ispin)%matrix, tmp_fm_struct)
         END DO

         CALL cp_cfm_create(eigrmat(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
         CALL cp_fm_struct_release(tmp_fm_struct)
         
      END DO

      ALLOCATE (diag(nao))

!----------------

!----------------
! DBCSR NMO x NMO
!----------------

    NULLIFY(op_cos_dbcsr, op_sin_dbcsr)
    CALL dbcsr_allocate_matrix_set(op_cos_dbcsr, dft_control%nspins)
    CALL dbcsr_allocate_matrix_set(op_sin_dbcsr, dft_control%nspins)

    DO ispin = 1, dft_control%nspins

       ALLOCATE(op_cos_dbcsr(ispin)%matrix)
       ALLOCATE(op_sin_dbcsr(ispin)%matrix)
!
!    ! Need to create from scratch with nmo size
       CALL dbcsr_copy(op_cos_dbcsr(ispin)%matrix, matrix_s(1)%matrix, 'COS_DBCSR_NMO')
       CALL dbcsr_copy(op_sin_dbcsr(ispin)%matrix, matrix_s(1)%matrix, 'SIN_DBCSR_NMO')
       CALL dbcsr_set(op_cos_dbcsr(ispin)%matrix, 0.0_dp)
       CALL dbcsr_set(op_sin_dbcsr(ispin)%matrix, 0.0_dp)
    END DO

!----------------

!----------------
! FM NAO x NAO
!----------------

      ! nao x nao
      IF(nmom >= 1) THEN

         NULLIFY (nao_fm_struct)
         ALLOCATE (nao_op_fm_set(2))
         ALLOCATE (nao_dop_fm_set(2))

         ALLOCATE (eeigrmat)
         ALLOCATE (deriveigrmat)
         ALLOCATE (deigrmat)
         ALLOCATE (inveigrmat)
         ALLOCATE (cumulant_matrix)
         CALL cp_fm_struct_create(nao_fm_struct, nrow_global=nao, &
                                  ncol_global=nao, para_env=para_env, context=mo_coeff%matrix_struct%context)
         DO i = 1, SIZE(nao_op_fm_set, 1)
            NULLIFY (nao_op_fm_set(i)%matrix)
            NULLIFY (nao_dop_fm_set(i)%matrix)
            CALL cp_fm_create(nao_op_fm_set(i)%matrix, nao_fm_struct)
            CALL cp_fm_create(nao_dop_fm_set(i)%matrix, nao_fm_struct)
         END DO

         CALL cp_cfm_create(eeigrmat%matrix, nao_op_fm_set(1)%matrix%matrix_struct)
         CALL cp_cfm_create(deriveigrmat%matrix, nao_op_fm_set(1)%matrix%matrix_struct)
         CALL cp_cfm_create(deigrmat%matrix, nao_op_fm_set(1)%matrix%matrix_struct)
         CALL cp_cfm_create(inveigrmat%matrix, nao_op_fm_set(1)%matrix%matrix_struct)
         CALL cp_cfm_create(cumulant_matrix%matrix, nao_op_fm_set(1)%matrix%matrix_struct)

      END IF
      CALL cp_fm_struct_release(nao_fm_struct)

!----------------

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
!            WRITE(*,*) "Dipole - Nuclear"
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
                  END DO
               ENDIF
            END DO
            zi = zi*zphase
            ci = AIMAG(LOG(zi))/twopi
            qq = AIMAG(LOG(zi))
            rmom(2:4, 2) = MATMUL(cell%hmat, ci)


            !----------------------------

               d01(:) = CMPLX(0._dp, 0._dp, dp)
               d11(:) = CMPLX(0._dp, 0._dp, dp)
               zi = 0._dp
               DO ia = 1, SIZE(particle_set)
                  atomic_kind => particle_set(ia)%atomic_kind
                  CALL get_atomic_kind(atomic_kind, kind_number=ikind)
                  CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)

                  ! Normalize charge
                  qn = charge/rmom(1,2)
                  IF (.NOT. ghost .AND. .NOT. floating) THEN
                     ria = particle_set(ia)%r
                     ria = pbc(ria, cell)
                     DO i = 1, 3
                        kvec(:) = twopi*cell%h_inv(i, :)
                        dd = SUM(kvec(:)*ria(:))
                        d01(i) = d01(i) + qn*CMPLX(COS(dd), -1.0_dp*SIN(dd), KIND=dp)
                        d11(i) = d11(i) + qn*ria(i)*CMPLX(COS(dd), -1.0_dp*SIN(dd), KIND=dp)
                     END DO
                  ENDIF
               END DO
               zi(:) = d11(:)/d01(:)
               newqq = REAL(zi, dp)*nelectron
               newrmom(2:4, 2) = newqq 



            !----------------------------








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

            DO i = 1, 3
               DO j = 1, 3
                  zij(i, j) = zij(i, j)*zphase(i)*zphase(j)
                  zz = zij(i, j)/zi(i)/zi(j)
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
            END DO

            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))

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

      tracei = 0._dp
!      strace = 0._dp

      ! MODFIED ONLY TAKES DIAGONAL DERIV
      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! OLD OLD OLD OLD OLD OLD
            ! Dipole
            DO i = 1, 3
               IF(para_env%mepos==0) WRITE(*,*) "DIPOLE DIRECTION i IN OLD", i
               kvec(:) = twopi*cell%h_inv(i, :)
               CALL build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)
               IF (qs_env%run_rtp) THEN
                  CALL get_qs_env(qs_env, rtp=rtp)
                  CALL get_rtp(rtp, mos_new=mos_new)
                  CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
               ELSE
                  CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
               END IF
!               WRITE(*,*) "cosmat" 
!               CALL dbcsr_print(cosmat)
!               WRITE(*,*) "op_fm_set - COS"
!               CALL cp_fm_write_formatted(op_fm_set(1, 1)%matrix, 6)
!               WRITE(*,*) "op_fm_set - SIN"
!               CALL cp_fm_write_formatted(op_fm_set(2, 1)%matrix, 6)
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
            ! OLD OLD OLD OLD OLD OLD

            ! Dipole
            DO i = 1, 3
               IF(para_env%mepos==0) WRITE(*,*) "DIPOLE DIRECTION i", i
               kvec(:) = twopi*cell%h_inv(i, :)
               ! d/dt (a|e^(ikr)|b)
               CALL build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec, nderivative = 1, &
                                               dkcos = dkcos(i)%matrix, dksin = dksin(i)%matrix)
!               ! d/dk (a|e^(ikr)|b)
!               CALL build_berry_moment_matrix2(qs_env, cosmat, sinmat, kvec, &
!                                               nderivative = 1, dkcos = dkcos, dksin = dksin)

               !-------------------------------------------
               ! INVERSE BEFORE TRANSFORMATION
               !-------------------------------------------
             
                ! dbcsr to fm
                ! nao x nao
                CALL copy_dbcsr_to_fm(cosmat, nao_op_fm_set(1)%matrix)
                CALL copy_dbcsr_to_fm(sinmat, nao_op_fm_set(2)%matrix)
                CALL copy_dbcsr_to_fm(dkcos(i)%matrix, nao_dop_fm_set(1)%matrix)
                CALL copy_dbcsr_to_fm(dksin(i)%matrix, nao_dop_fm_set(2)%matrix)

!                WRITE(*,*) "cosmat"
!                !CALL dbcsr_print(cosmat)
!                CALL dbcsr_print(cosmat, matlab_format = .TRUE.)
!                WRITE(*,*) "sinmat"
!                !CALL dbcsr_print(sinmat)
!                CALL dbcsr_print(sinmat, matlab_format = .TRUE.)
!                WRITE(*,*) "dkcos(i)%matrix"
!                !CALL dbcsr_print(dkcos(i)%matrix)
!                CALL dbcsr_print(dkcos(i)%matrix, matlab_format = .TRUE.)
!                WRITE(*,*) "dksin(i)%matrix"
!                CALL dbcsr_print(dksin(i)%matrix, matlab_format = .TRUE.)
!                WRITE(*,*) "nao_op_fm_set - COS"
!                CALL cp_fm_write_formatted(nao_op_fm_set(1)%matrix, 6)
!                WRITE(*,*) "nao_op_fm_set - SIN"
!                CALL cp_fm_write_formatted(nao_op_fm_set(2)%matrix, 6)
!                WRITE(*,*) "nao_dop_fm_set - COS"
!                CALL cp_fm_write_formatted(nao_dop_fm_set(1)%matrix, 6)
!                WRITE(*,*) "nao_dop_fm_set - SIN"
!                CALL cp_fm_write_formatted(nao_dop_fm_set(2)%matrix, 6)

                !------------------------------------------------------------------------------
                !------------------------------------------------------------------------------
                
                et1 = CMPLX(0.0_dp, 0.0_dp, dp)
                DO ispin = 1, dft_control%nspins
                
                   CALL cp_cfm_get_info(deriveigrmat%matrix, ncol_local = ncol_local, &
                                        nrow_local = nrow_local)
                   DO idim = 1, ncol_local
                
                      ! Make complex matrix of d01
                      eeigrmat%matrix%local_data(:, idim) = &
                         CMPLX(nao_op_fm_set(1)%matrix%local_data(:, idim), &
                               -nao_op_fm_set(2)%matrix%local_data(:, idim), dp)
                      ! Make complex matrix of d11
                      deigrmat%matrix%local_data(:, idim) = &
                         CMPLX(nao_dop_fm_set(1)%matrix%local_data(:, idim), &
                               -nao_dop_fm_set(2)%matrix%local_data(:, idim), dp)
                   END DO
                   
                   ! Invert d01
                   CALL cp_cfm_lu_invert(eeigrmat%matrix)
                   ! element by element inversion
                   ! nao x nao
                   DO icol = 1, ncol_local
                      DO irow = 1, nrow_local
                           relem = nao_op_fm_set(1)%matrix%local_data(irow, icol)
                           ielem = -1.0_dp*nao_op_fm_set(2)%matrix%local_data(irow, icol)
!                           ! set cfm matrix element to its reciprocal
!                           IF(relem == 0.0_dp .AND. ielem == 0.0_dp) THEN
!                              inveigrmat%matrix%local_data(irow, icol) = CMPLX(0.0_dp, 0.0_dp, dp)
!                           ELSE
!                              inveigrmat%matrix%local_data(irow, icol) = &
!                                 CMPLX(relem/(relem**2+ielem**2), -1.0_dp*ielem/(relem**2+ielem**2), dp)
!                           END IF
                      END DO
                   END DO

                   CALL cp_cfm_to_fm(msource = deigrmat%matrix, &
                                     mtargetr = nao_op_fm_set(1)%matrix, &
                                     mtargeti = nao_op_fm_set(2)%matrix)
!                   WRITE(*,*) "deriveigrmat - COS deriv"
!                   CALL cp_fm_write_formatted(nao_op_fm_set(1)%matrix, 6)
!                   WRITE(*,*) "deriveigrmat - SIN deriv"
!                   CALL cp_fm_write_formatted(nao_op_fm_set(2)%matrix, 6)
!
!                   CALL cp_cfm_to_fm(msource = eeigrmat%matrix, &
!                                     mtargetr = nao_op_fm_set(1)%matrix, &
!                                     mtargeti = nao_op_fm_set(2)%matrix)
!                   WRITE(*,*) "nao_op_fm_set - COS inv"
!                   CALL cp_fm_write_formatted(nao_op_fm_set(1)%matrix, 6)
!                   WRITE(*,*) "nao_op_fm_set - SIN inv"
!                   CALL cp_fm_write_formatted(nao_op_fm_set(2)%matrix, 6)

                   !---------------------------------------------
                   ! tr(d01^-1 * d11)
                   CALL cp_cfm_trace(matrix_a = eeigrmat%matrix, &
                                     matrix_b = deigrmat%matrix, &
                                     trace = et1)
                
                   etr = z_one*et1
                
                   ! Closed shell
                   IF (dft_control%nspins == 1) THEN
                      etr = etr + et1
                   ENDIF
                   !---------------------------------------------

                   !---------------------------------------------
!                   CALL cp_cfm_get_info(cumulant_matrix%matrix, ncol_local = ncol_local)
!
!                   ! matrix_c = alpha * op1( matrix_a ) * op2( matrix_b ) + beta*matrix_c.
!                   CALL cp_cfm_gemm("N","N", ncol_local, ncol_local, ncol_local, &
!                                    z_one, eeigrmat%matrix, deigrmat%matrix, &
!                                    z_zero, cumulant_matrix%matrix)
!
!!                   CALL  cp_cfm_schur_product(matrix_a = eeigrmat%matrix,&
!!                                              matrix_b = deigrmat%matrix,&
!!                                              matrix_c = cumulant_matrix%matrix)
!
!                   ! cfm -> fm
!                   ! nmo x nmo 
!                   CALL cp_cfm_to_fm(msource = cumulant_matrix%matrix, &
!                                     mtargetr = nao_op_fm_set(1)%matrix, &
!                                     mtargeti = nao_op_fm_set(2)%matrix)
!   
!                   ! fm -> dbcsr
!                   CALL copy_fm_to_dbcsr(nao_op_fm_set(1)%matrix, cosmat)
!                   CALL copy_fm_to_dbcsr(nao_op_fm_set(2)%matrix, sinmat)
!
!                   ! Rotation onto MOs
!                   CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
!
!                   ! fm -> dbcsr
!                   CALL copy_fm_to_dbcsr(op_fm_set(1, 1)%matrix, cosmat)
!                   CALL copy_fm_to_dbcsr(op_fm_set(2, 1)%matrix, sinmat)
!   
!                   ! calculate trace
!                   CALL dbcsr_trace(cosmat, relem)
!                   CALL dbcsr_trace(sinmat, ielem)
!                   et1 =  CMPLX(relem, ielem, dp)
!   
!                   etr = z_one*et1
!               
!                   ! Closed shell
!                   IF (dft_control%nspins == 1) THEN
!                      etr = etr + et1
!                   ENDIF
                  !---------------------------------------------

                END DO
                ctrace(i) = etr / CMPLX(0.0_dp, kvec(i), dp)
            
                !------------------------------------------------------------------------------
                !------------------------------------------------------------------------------

!                ! Store dkcos_fm and dksin_fm in complex matrix
!                ! nao x nao
!                CALL cp_cfm_get_info(deriveigrmat%matrix, ncol_local = ncol_local, nrow_local = nrow_local)
!
!                ! element by element inversion
!                ! nao x nao
!                DO icol = 1, ncol_local
!                   DO irow = 1, nrow_local
!                        relem = nao_op_fm_set(1)%matrix%local_data(irow, icol)
!                        ielem = -1.0_dp*nao_op_fm_set(2)%matrix%local_data(irow, icol)
!                        ! set cfm matrix element to its reciprocal
!                        IF(relem == 0.0_dp .AND. ielem == 0.0_dp) THEN
!                           inveigrmat%matrix%local_data(irow, icol) = CMPLX(0.0_dp, 0.0_dp, dp)
!                        ELSE
!                           inveigrmat%matrix%local_data(irow, icol) = &
!                              CMPLX(relem/(relem**2+ielem**2), -1.0_dp*ielem/(relem**2+ielem**2), dp)
!                        END IF
!                   END DO
!                END DO
!                ! Store first deriv as complex matrix
!                CALL cp_fm_scale(-1.0_dp, nao_dop_fm_set(2)%matrix)
!                CALL cp_fm_to_cfm(msourcer = nao_dop_fm_set(1)%matrix, &
!                                  msourcei = nao_dop_fm_set(2)%matrix, &
!                                  mtarget = deriveigrmat%matrix)
!
!                CALL cp_cfm_to_fm(msource = deriveigrmat%matrix, &
!                                  mtargetr = nao_op_fm_set(1)%matrix, &
!                                  mtargeti = nao_op_fm_set(2)%matrix)
!                WRITE(*,*) "deriveigrmat - COS deriv"
!                CALL cp_fm_write_formatted(nao_op_fm_set(1)%matrix, 6)
!                WRITE(*,*) "deriveigrmat - SIN deriv"
!                CALL cp_fm_write_formatted(nao_op_fm_set(2)%matrix, 6)
!                CALL cp_cfm_to_fm(msource = inveigrmat%matrix, &
!                                  mtargetr = nao_op_fm_set(1)%matrix, &
!                                  mtargeti = nao_op_fm_set(2)%matrix)
!                WRITE(*,*) "nao_op_fm_set - COS inv"
!                CALL cp_fm_write_formatted(nao_op_fm_set(1)%matrix, 6)
!                WRITE(*,*) "nao_op_fm_set - SIN inv"
!                CALL cp_fm_write_formatted(nao_op_fm_set(2)%matrix, 6)
!
!                ! Multiplying Hadamard 
!                CALL cp_cfm_schur_product(inveigrmat%matrix, deriveigrmat%matrix, &
!                                          cumulant_matrix%matrix)
!                ! cfm -> fm
!                ! nao x nao 
!                CALL cp_cfm_to_fm(msource = cumulant_matrix%matrix, &
!                                  mtargetr = nao_op_fm_set(1)%matrix, &
!                                  mtargeti = nao_op_fm_set(2)%matrix)
!
!!                WRITE(*,*) "nao_op_fm_set - COS - AFTER HADAMARD"
!!                CALL cp_fm_write_formatted(nao_op_fm_set(1)%matrix, 6)
!!                WRITE(*,*) "nao_op_fm_set - SIN - AFTER HADAMARD"
!!                CALL cp_fm_write_formatted(nao_op_fm_set(2)%matrix, 6)
!
!                strace = 0.0_dp
!                sintrace = 0.0_dp
!
!                ! fm -> dbcsr
!                CALL copy_fm_to_dbcsr(nao_op_fm_set(1)%matrix, cosmat)
!                CALL copy_fm_to_dbcsr(nao_op_fm_set(2)%matrix, sinmat)
!
!!                CALL dbcsr_multiply("N","N", 1.0_dp, rho_ao(1)%matrix, sinmat, 1.0_dp, &
!!                                    cosmat)
!                
!!                strace = 0.0_dp
!!                sintrace = 0.0_dp
!!                DO ispin = 1, dft_control%nspins
!!                   CALL dbcsr_trace(rho_ao(ispin)%matrix, sinmat, sintrace)
!!
!!                   IF(para_env%mepos==0) WRITE(*,*) "TRACE tr(rho_ao*d11/d01)_nao", sintrace
!!
!!                   CALL dbcsr_trace(sinmat, sintrace)
!!!                   CALL dbcsr_trace(cosmat, sintrace)
!!
!!                   WRITE(*,*) "TRACE tr(d11/d01)_nao", sintrace
!!                   strace = strace + sintrace
!!
!!                   ! Closed shell
!!                   IF (dft_control%nspins == 1) THEN
!!                     strace = strace + sintrace
!!                   ENDIF
!!
!!                END DO
!
!!                WRITE(*,*) "AFTER HADAMARD - COS"
!!                CALL dbcsr_print(cosmat)
!!                WRITE(*,*) "AFTER HADAMARD - SIN"
!!                CALL dbcsr_print(sinmat)
!
!                ! project cosmat and dkcos
!                ! nao x nao -> nmo x nmo
!                ! spin
!                CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
!
!                DO ispin = 1, dft_control%nspins 
!                   ! fm -> dbcsr
!                   ! nmo x nmo
!                   CALL copy_fm_to_dbcsr(op_fm_set(1, ispin)%matrix, op_cos_dbcsr(ispin)%matrix)
!                   CALL copy_fm_to_dbcsr(op_fm_set(2, ispin)%matrix, op_sin_dbcsr(ispin)%matrix)
!
!!                   WRITE(*,*) "op_sin_dbcsr"
!!                   CALL dbcsr_print(op_sin_dbcsr(ispin)%matrix)
!
!                   ! trace
!                   CALL dbcsr_trace(op_sin_dbcsr(ispin)%matrix, sintrace)
!                   WRITE(*,*) "TRACE tr(d11/d01)_nmo", sintrace
!                    
!                   strace = strace + sintrace
!
!                   ! Closed shell
!                   IF (dft_control%nspins == 1) THEN
!                     strace = strace + sintrace
!                   ENDIF
!
!                END DO
!                ctrace(i) = strace / CMPLX(0.0_dp, kvec(i), dp)
            END DO ! i

!            ctrace = ctrace*xphase

           IF(para_env%mepos==0)  WRITE(*,*) "ctrace(:)", ctrace(:)

         CASE (2)
            ! Quadrupole
!            CPABORT("Berry phase moments bigger than 1 not implemented")
            DO i = 1, 3
               DO j = i, 3
                  kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :))
                  ! (a|e^ikr|b)
                  CALL build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec)
                                                 
                  IF (qs_env%run_rtp) THEN
                     CALL get_qs_env(qs_env, rtp=rtp)
                     CALL get_rtp(rtp, mos_new=mos_new)
                     CALL op_orbbas_rtp(cosmat, sinmat, mos, op_fm_set, mos_new)
                  ELSE

                  CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
                  ! First derivs 
                  DO k = 1, 3
                     CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
                  END DO
                  ! Second derivs
                     DO k = 1, 6
                        cosmat => dkcos(k)%matrix
                        sinmat => dksin(k)%matrix
                        !!!!!!!! OVERWRITING OF op_fm_set
                        CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
                     END DO

                  END IF
                  zdet = CMPLX(1._dp, 0._dp, dp)
                  DO ispin = 1, dft_control%nspins
                     CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=tmp_dim)
                     DO idim = 1, tmp_dim
                        eigrmat(ispin)%matrix%local_data(:, idim) = &
                           ! e^ikr -> exp(-ikr) = cos(kr) - i*sin(kr)
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
!------------------------------------------------------------------------------
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

            tracei = REAL(ctrace,dp)
!            tracei = AIMAG(ctrace)

           ! shift new elec based on old nuclear
           IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
           IF(para_env%mepos==0) WRITE(*,*) "Berry nuclear + Greg elec"
           IF(para_env%mepos==0) WRITE(*,*) "tracei", tracei
            DO i = 1, 3
               IF(para_env%mepos==0) WRITE(*,*) "START"
               IF(para_env%mepos==0) WRITE(*,*) "iteration", i
               IF(para_env%mepos==0) WRITE(*,*) "qq(i) =", qq(i)
               IF(para_env%mepos==0) WRITE(*,*) "qq(i)+tracei(i) =", qq(i)+tracei(i)
               DO WHILE  (qq(i)+tracei(i) > pi)
                  IF(para_env%mepos==0) WRITE(*,*) "First Do While", qq(i)+tracei(i)
                  IF(para_env%mepos==0) WRITE(*,*) "tracei(i)", tracei(i)
                  tracei(i) = tracei(i)-twopi
                  IF(para_env%mepos==0) WRITE(*,*) "qq(i)+tracei(i)", qq(i)+tracei(i)
               END DO
               DO WHILE (qq(i)+tracei(i) < -pi)
                  IF(para_env%mepos==0) WRITE(*,*) "Second Do While", qq(i)+tracei(i)
                  IF(para_env%mepos==0) WRITE(*,*) "tracei(i)", tracei(i)
                  tracei(i) = tracei(i)+twopi
                  IF(para_env%mepos==0) WRITE(*,*) "qq(i)+tracei(i)", qq(i)+tracei(i)
               END DO
               IF(para_env%mepos==0) WRITE(*,*) "criteria  after second DO WHILE", qq(i)+tracei(i)
            END DO
            IF(para_env%mepos==0) WRITE(*,*) "tracei AFTER", tracei
!            ! Shift electronic based on new nuclear
!            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
!            IF(para_env%mepos==0) WRITE(*,*) "Greg nuclear + elec"
!            DO i = 1, 3
!               IF(para_env%mepos==0) WRITE(*,*) "iteration", i
!               IF(para_env%mepos==0) WRITE(*,*) "newqq(i) =", newqq(i)
!               IF(para_env%mepos==0) WRITE(*,*) "newqq(i)+tracei2(i) =", newqq(i)+tracei(i)
!               IF(para_env%mepos==0) WRITE(*,*) "START"
!               DO WHILE  (newqq(i)+tracei(i) > pi)
!                  tracei(i) = tracei(i)-twopi
!                  IF(para_env%mepos==0) WRITE(*,*) "newqq(i)+tracei(i)", newqq(i)+tracei(i)
!               END DO
!               IF(para_env%mepos==0) WRITE(*,*) "tracei(i) after -2pi", newqq(i)+tracei(i)
!               DO WHILE (newqq(i)+tracei(i) < -pi)
!                  tracei(i) = tracei(i)+twopi
!                  IF(para_env%mepos==0) WRITE(*,*) "newqq(i)+tracei(i)", newqq(i)+tracei(i)
!               END DO
!               IF(para_env%mepos==0) WRITE(*,*) "criteria  after +2pi", newqq(i)+tracei(i)
!            END DO
!            IF(para_env%mepos==0) WRITE(*,*) "tracei AFTER", tracei
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)

            ! Derivation wrt t pulls down additional k (2pi/L)
!            newrmom(2:4, 1) = MATMUL(cell%hmat, tracei)/twopi
            ! Derivation wrt k
            !newrmom(2:4, 1) = 2.0_dp*tracei
            newrmom(2:4, 1) = tracei
            
            IF(para_env%mepos==0) THEN 
               WRITE(*,*) ""
               WRITE(*,*) "nuclear", rmom(2:4, 2)
               WRITE(*,*) "rmom", rmom(2:4, 1)
               WRITE(*,*) "new nuclear", newrmom(2:4, 2)
               WRITE(*,*) "newrmom", newrmom(2:4, 1)
               WRITE(*,*) ""
               WRITE(*,*) "total", rmom(2:4, 1)+rmom(2:4, 2)
               WRITE(*,*) "newtotal", -1.0_dp*(newrmom(2:4, 1)+newrmom(2:4, 2))

               WRITE(*,*) "DIVIDE: old / new"

               IF (newrmom(2,1) .NE. 0) THEN
                  WRITE(*,*) rmom(2,1) / newrmom(2, 1)
                  WRITE(*,*) "Error X [%]", 100.0_dp*(rmom(2,1) - newrmom(2, 1)) / newrmom(2, 1)
               ELSE
                  WRITE(*,*) 0.0_dp
               END IF
               IF (newrmom(3,1) .NE. 0) THEN
                  WRITE(*,*) rmom(3,1) / newrmom(3,1)
                  WRITE(*,*) "Error Y [%]", 100.0_dp*(rmom(3,1) - newrmom(3, 1)) / newrmom(3, 1)
               ELSE
                  WRITE(*,*) 0.0_dp
               END IF
               IF (newrmom(4,1) .NE. 0) THEN
                  WRITE(*,*) rmom(4,1) / newrmom(4,1)
                  WRITE(*,*) "Error Z [%]", 100.0_dp*(rmom(4,1) - newrmom(4, 1)) / newrmom(4, 1)
               ELSE
                  WRITE(*,*) 0.0_dp
               END IF

!               WRITE(*,*) "cell%h_inv", cell%h_inv
!               WRITE(*,*) "cell%h_mat", cell%hmat
!               WRITE(*,*) "2*cell%h_inv*twopi", 2*twopi*cell%h_inv
!               WRITE(*,*) "cell%hmat/twopi", cell%hmat/twopi
            END IF

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
 
            ! Hardcoding of upper triangular matrix elements of cij
            ! into moment vector
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
!------------------------------------------------------------------------------
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
      DEALLOCATE (newrmom)
      DEALLOCATE (rlab)
      IF (magnetic) THEN
         DEALLOCATE (mmom)
      END IF

      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)

      CALL dbcsr_deallocate_matrix_set(op_cos_dbcsr)
      CALL dbcsr_deallocate_matrix_set(op_sin_dbcsr)

      CALL dbcsr_deallocate_matrix_set(dkcos)
      CALL dbcsr_deallocate_matrix_set(dksin)

      CALL cp_cfm_release(eeigrmat%matrix)
      CALL cp_cfm_release(deigrmat%matrix)
      CALL cp_cfm_release(inveigrmat%matrix)
      CALL cp_cfm_release(cumulant_matrix%matrix)
      CALL cp_cfm_release(deriveigrmat%matrix)

      DO ispin = 1, dft_control%nspins
         CALL cp_fm_release(opvec(ispin)%matrix)
         CALL cp_cfm_release(eigrmat(ispin)%matrix)
         DO i = 1, SIZE(op_fm_set, 1)
            CALL cp_fm_release(op_fm_set(i, ispin)%matrix)
            CALL cp_fm_release(tmp_op_fm_set(i, ispin)%matrix)
            CALL cp_fm_release(dop_fm_set(i, ispin)%matrix)
         END DO
      END DO
      DO i = 1, SIZE(nao_op_fm_set, 1)
         CALL cp_fm_release(nao_op_fm_set(i)%matrix)
         CALL cp_fm_release(nao_dop_fm_set(i)%matrix)
      END DO

      DEALLOCATE (op_fm_set)
      DEALLOCATE (dop_fm_set)
      DEALLOCATE (tmp_op_fm_set)
      DEALLOCATE (opvec)
      DEALLOCATE (eigrmat)
      DEALLOCATE (eeigrmat)
      DEALLOCATE (deigrmat)
      DEALLOCATE (inveigrmat)
      DEALLOCATE (deriveigrmat)
      DEALLOCATE (cumulant_matrix)

      DEALLOCATE (nao_op_fm_set)
      DEALLOCATE (nao_dop_fm_set)


      DEALLOCATE (diag)

      CALL timestop(handle)

   END SUBROUTINE qs_moment_berry_phase3

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param magnetic ...
!> \param nmoments ...
!> \param reference ...
!> \param ref_point ...
!> \param unit_number ... 
! **************************************************************************************************
!   SUBROUTINE qs_moment_berry_num(qs_env, magnetic, nmoments, reference, ref_point, unit_number)
!
!      TYPE(qs_environment_type), POINTER                 :: qs_env
!      LOGICAL, INTENT(IN)                                :: magnetic
!      INTEGER, INTENT(IN)                                :: nmoments, reference
!      REAL(dp), DIMENSION(:), POINTER                    :: ref_point
!      INTEGER, INTENT(IN)                                :: unit_number
!
!      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_moment_berry_num', &
!         routineP = moduleN//':'//routineN
!
!
!
!
!
!   END SUBROUTINE qs_moment_berry_num


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
                                                            iy, iz, j, l, nao, nm, nmo, nmom, &
                                                            nmotot, tmp_dim
      LOGICAL                                            :: floating, ghost, uniform
      REAL(dp)                                           :: charge, ci(3), cij(3, 3), dd, occ, trace
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: mmom
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: rmom, newrmom, newrmom2
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

      ! fbelle
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: dop_fm_set, tmp_op_fm_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dkcos, dksin
      TYPE(cp_cfm_p_type), DIMENSION(:), POINTER         :: deigrmat, cumulant_matrix
      REAL(dp)                                           :: tracei(3), tracei2(3)
!      REAL(dp)                                           :: strace
!      REAL(dp)                                           :: relem, ielem
!      REAL(dp)                                           :: tracei(3)
!      REAL(dp)                                           :: sintrace 
!      INTEGER                                            :: irow, icol, ncol_local, nrow_local
!      INTEGER                                            :: ncol_local!, nrow_local
!      COMPLEX(dp)                                        :: det1
      COMPLEX(dp)                                        :: ctrace(3), ctrace2(3), et1
      INTEGER                                            :: nelectron
      COMPLEX(dp)                                        :: etr
      REAL(dp), DIMENSION(3)                             :: newqq
      COMPLEX(dp), DIMENSION(3)                          :: d01, d11
      COMPLEX(dp)                                        :: ii
      REAL(dp)                                           :: norm, qn

      CPASSERT(ASSOCIATED(qs_env))

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
      ALLOCATE (newrmom(nm+1, 3))
      newrmom = 0.0_dp
      ALLOCATE (newrmom2(nm+1, 3))
      newrmom2 = 0.0_dp
      rlab = ""
      IF (magnetic) THEN
         nm = 3
         ALLOCATE (mmom(nm))
         mmom = 0._dp
      END IF

     ! fbelle
     ! DEBUGGING OF DERIVATIVES
     ! CALL debug_berry_moment_matrix(qs_env, nmoments = 1)

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
                      mos=mos,&
                      nelectron_total = nelectron)

      IF(para_env%mepos==0) THEN
         WRITE(*,*) "qs_moment_berry_phase2 - BEGINNING"
      END IF

      CALL qs_rho_get(rho, rho_ao=rho_ao)

      NULLIFY (cosmat, sinmat)
      ALLOCATE (cosmat, sinmat)
      CALL dbcsr_copy(cosmat, matrix_s(1)%matrix, 'COS MOM')
      CALL dbcsr_copy(sinmat, matrix_s(1)%matrix, 'SIN MOM')
      CALL dbcsr_set(cosmat, 0.0_dp)
      CALL dbcsr_set(sinmat, 0.0_dp)

      ALLOCATE (op_fm_set(2, dft_control%nspins))
      ALLOCATE (tmp_op_fm_set(2, dft_control%nspins))
      ALLOCATE (opvec(dft_control%nspins))
      ALLOCATE (eigrmat(dft_control%nspins))

      ! Initialize and allocate matrices to store derivatives
      IF(nmom >= 1) THEN
         ALLOCATE (dop_fm_set(2, dft_control%nspins))
         ALLOCATE (deigrmat(dft_control%nspins))
         ALLOCATE (cumulant_matrix(dft_control%nspins))

         NULLIFY(dkcos, dksin)
         CALL dbcsr_allocate_matrix_set(dkcos, 3)
         CALL dbcsr_allocate_matrix_set(dksin, 3)
         DO i = 1, 3
            ALLOCATE(dkcos(i)%matrix, dksin(i)%matrix)
            CALL dbcsr_copy(dkcos(i)%matrix, matrix_s(1)%matrix, 'DKCOS MOM')
            CALL dbcsr_copy(dksin(i)%matrix, matrix_s(1)%matrix, 'DKSIN MOM')
            CALL dbcsr_set(dkcos(i)%matrix, 0.0_dp)
            CALL dbcsr_set(dksin(i)%matrix, 0.0_dp)
         END DO
      END IF

      ! fbelle - Store projected derivatives 
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
            NULLIFY (dop_fm_set(i, ispin)%matrix)
            NULLIFY (tmp_op_fm_set(i, ispin)%matrix)
            CALL cp_fm_create(op_fm_set(i, ispin)%matrix, tmp_fm_struct)
            CALL cp_fm_create(dop_fm_set(i, ispin)%matrix, tmp_fm_struct)
            CALL cp_fm_create(tmp_op_fm_set(i, ispin)%matrix, tmp_fm_struct)
         END DO

         CALL cp_cfm_create(eigrmat(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
         CALL cp_cfm_create(deigrmat(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
         CALL cp_cfm_create(cumulant_matrix(ispin)%matrix, op_fm_set(1, ispin)%matrix%matrix_struct)
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

      IF(para_env%mepos==0) WRITE(*,*) "rcc", rcc(:)

      ! label
      DO l = 1, nm
         ix = indco(1, l+1)
         iy = indco(2, l+1)
         iz = indco(3, l+1)
         CALL set_label(rlab(l+1), ix, iy, iz)
      END DO

      ! nuclear contribution

      ! charge
      DO ia = 1, SIZE(particle_set)
         atomic_kind => particle_set(ia)%atomic_kind
         CALL get_atomic_kind(atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)
         IF (.NOT. ghost .AND. .NOT. floating) THEN
            rmom(1, 2) = rmom(1, 2)-charge
            newrmom(1, 2) = newrmom(1, 2)-charge
            newrmom2(1, 2) = newrmom2(1, 2)-charge
         ENDIF
      END DO

      ria = twopi*MATMUL(cell%h_inv, rcc)
      zphase = CMPLX(COS(ria), SIN(ria), dp)**rmom(1, 2)
 
!      IF(para_env%mepos==0) WRITE(*,*) "periodicpos", ria
!      IF(para_env%mepos==0) WRITE(*,*) "Reazphase", REAL(zphase,dp)
!      IF(para_env%mepos==0) WRITE(*,*) "Imgzphase", AIMAG(zphase)

      zi = 0._dp
      zij = 0._dp
      zijk = 0._dp
      zijkl = 0._dp

      ii = CMPLX(0.0_dp, 1.0_dp, KIND=dp)
      norm = SQRT(rmom(1, 2)**2)
      IF(para_env%mepos==0) WRITE(*,*) "NORM", norm 

      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
            ! OLD OLD OLD OLD OLD OLD OLD OLD OLD
!            WRITE(*,*) "Dipole - Nuclear"
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
            IF(para_env%mepos==0) WRITE(*,*) "Nuclear OLD" 
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
            zi(:) = CMPLX(1._dp, 0._dp, dp)
            DO ia = 1, SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)
!               IF(para_env%mepos==0) WRITE(*,*) "ATOM", atomic_kind%name
               IF (.NOT. ghost .AND. .NOT. floating) THEN
                  ria = particle_set(ia)%r
                  ria = pbc(ria, cell)
                  DO i = 1, 3
!                     IF(para_env%mepos==0) WRITE(*,*) "DIRECTION", i
                     kvec(:) = twopi*cell%h_inv(i, :)
                     dd = SUM(kvec(:)*ria(:))
!                     IF(para_env%mepos==0) WRITE(*,*) "dd", dd
                     zdeta = CMPLX(COS(dd), SIN(dd), KIND=dp)**charge
!                     IF(para_env%mepos==0) WRITE(*,*) "zdeta", zdeta
                     zi(i) = zi(i)*zdeta
                  END DO
               ENDIF
            END DO
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 10)
            zi = zi*zphase
            ci = AIMAG(LOG(zi))/twopi
!            IF(para_env%mepos==0) WRITE(*,*) "ci", ci
            qq = AIMAG(LOG(zi))
!            IF(para_env%mepos==0) WRITE(*,*) "qq", qq
            rmom(2:4, 2) = MATMUL(cell%hmat, ci)
            IF(para_env%mepos==0) WRITE(*,*) "Nuclear old", rmom(2:4, 2)
            ! OLD OLD OLD OLD OLD OLD OLD OLD OLD

!            ! NEW NEW NEW NEW 
!            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
!            IF(para_env%mepos==0) WRITE(*,*) "Nuclear New" 
!            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
!            d01(:) = CMPLX(0._dp, 0._dp, dp)
!            d11(:) = CMPLX(0._dp, 0._dp, dp)
!            zi = 0._dp
!            DO ia = 1, SIZE(particle_set)
!               atomic_kind => particle_set(ia)%atomic_kind
!               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
!               CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)
!
!               qn = charge/norm
!!               IF(para_env%mepos==0) WRITE(*,*) "qn", qn
!!               IF(para_env%mepos==0) WRITE(*,*) "ATOM", atomic_kind%name
!               IF (.NOT. ghost .AND. .NOT. floating) THEN
!                  ria = particle_set(ia)%r
!                  ria = pbc(ria, cell)
!                  IF(para_env%mepos==0) WRITE(*,*) "ria", ria
!
!                  DO i = 1, 3
!                     IF(para_env%mepos==0) WRITE(*,*) "DIRECTION", i
!                     kvec(:) = twopi*cell%h_inv(i, :)
!                     IF(para_env%mepos==0) WRITE(*,*) "kvec", kvec(:)
!                     dd = SUM(kvec(:)*ria(:))
!!                     IF(para_env%mepos==0) WRITE(*,*) "dd", dd
!!                     IF(para_env%mepos==0) WRITE(*,*) "d01 (each)", d01(i)
!!                     IF(para_env%mepos==0) WRITE(*,*) "d11 (each)", d11(i)
!
!                     ! ik-derivation
!!                     d01(i) = d01(i) + qn*CMPLX(COS(dd), -1.0_dp*SIN(dd), KIND=dp)
!!                     d11(i) = d11(i) + qn*ria(i)*CMPLX(COS(dd), -1.0_dp*SIN(dd), KIND=dp)
!                     ! t-derivation
!                     d01(i) = d01(i) + qn*CMPLX(COS(dd), -1.0_dp*SIN(dd), KIND=dp)
!                     d11(i) = d11(i) + qn*CMPLX(0.0_dp, -1.0_dp*dd, dp)*&
!                                       CMPLX(COS(dd), -1.0_dp*SIN(dd), KIND=dp)
!                     IF(para_env%mepos==0) WRITE(*,*) "d01", d01(i)
!                     IF(para_env%mepos==0) WRITE(*,*) "d11", d11(i)
!                  END DO
!               ENDIF
!            END DO
!
!            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 10)
!            zi(:) = d11(:)/d01(:)
!            IF(para_env%mepos==0) WRITE(*,*) "zi", zi
!
!            ! t-deriv
!            DO i = 1, 3
!               kvec(:) = twopi*cell%h_inv(i, :)
!               zi(:) = d11(:)/d01(:) / CMPLX(0.0_dp, -1.0_dp*kvec(i), dp)
!            END DO
!            ci = REAL(zi, KIND=dp)*nelectron
!            newqq = REAL(zi, dp)*nelectron
!
!            ! ik-deriv
!!            ci = REAL(zi, KIND=dp)*nelectron
!!            newqq = REAL(zi, dp)*nelectron
!            IF(para_env%mepos==0) WRITE(*,*) "ci", ci
!            IF(para_env%mepos==0) WRITE(*,*) "newqq", newqq
!            
!!            newrmom(2:4, 2) = newqq
!            newrmom(2:4, 2) = ci
!            IF(para_env%mepos==0) WRITE(*,*) "Nuclear new", newrmom(2:4, 2)
!            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
!            ! NEW NEW NEW NEW 

            ! NEW NEW NEW NEW 
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
            IF(para_env%mepos==0) WRITE(*,*) "Nuclear New" 
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
            d01(:) = CMPLX(0._dp, 0._dp, dp)
            d11(:) = CMPLX(0._dp, 0._dp, dp)
            zi = 0._dp
            DO ia = 1, SIZE(particle_set)
               atomic_kind => particle_set(ia)%atomic_kind
               CALL get_atomic_kind(atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), core_charge=charge, ghost=ghost, floating=floating)

               qn = charge/norm
!               IF(para_env%mepos==0) WRITE(*,*) "qn", qn
!               IF(para_env%mepos==0) WRITE(*,*) "ATOM", atomic_kind%name
               IF (.NOT. ghost .AND. .NOT. floating) THEN
                  ria = particle_set(ia)%r
                  ria = pbc(ria, cell)
                  ria = rcc - ria 
                  IF(para_env%mepos==0) WRITE(*,*) "ria", ria

                  DO i = 1, 3
                     IF(para_env%mepos==0) WRITE(*,*) "DIRECTION", i
                     kvec(:) = twopi*cell%h_inv(i, :)
                     IF(para_env%mepos==0) WRITE(*,*) "kvec", kvec(:)
                     dd = SUM(kvec(:)*ria(:))
!                     IF(para_env%mepos==0) WRITE(*,*) "dd", dd
!                     IF(para_env%mepos==0) WRITE(*,*) "d01 (each)", d01(i)
!                     IF(para_env%mepos==0) WRITE(*,*) "d11 (each)", d11(i)

                     ! ik-derivation
!                     d01(i) = d01(i) + qn*CMPLX(COS(dd), -1.0_dp*SIN(dd), KIND=dp)
!                     d11(i) = d11(i) + qn*ria(i)*CMPLX(COS(dd), -1.0_dp*SIN(dd), KIND=dp)
                     ! t-derivation
                     d01(i) = d01(i) + qn*CMPLX(COS(dd), SIN(dd), KIND=dp)
                     d11(i) = d11(i) + qn*CMPLX(0.0_dp, dd, dp)*&
                                       CMPLX(COS(dd), SIN(dd), KIND=dp)
                     IF(para_env%mepos==0) WRITE(*,*) "d01", d01(i)
                     IF(para_env%mepos==0) WRITE(*,*) "d11", d11(i)
                  END DO
               ENDIF
            END DO

            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 10)
            zi(:) = d11(:)/d01(:)
            IF(para_env%mepos==0) WRITE(*,*) "zi", zi

            ! t-deriv
            DO i = 1, 3
               kvec(:) = twopi*cell%h_inv(i, :)
               zi(:) = d11(:)/d01(:) / CMPLX(0.0_dp, -1.0_dp*kvec(i), dp)
            END DO
            ci = REAL(zi, KIND=dp)*nelectron
            newqq = REAL(zi, dp)*nelectron

            ! ik-deriv
!            ci = REAL(zi, KIND=dp)*nelectron
!            newqq = REAL(zi, dp)*nelectron
            IF(para_env%mepos==0) WRITE(*,*) "ci", ci
            IF(para_env%mepos==0) WRITE(*,*) "newqq", newqq
            
!            newrmom(2:4, 2) = newqq
            newrmom(2:4, 2) = ci
            IF(para_env%mepos==0) WRITE(*,*) "Nuclear new", newrmom(2:4, 2)
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
            ! NEW NEW NEW NEW 




         CASE (2)
            ! Quadrupole
            CPABORT("Berry phase moments bigger than 1 not implemented")
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

!      IF(para_env%mepos==0) WRITE(*,*) "Reaxphase", REAL(xphase,dp)
!      IF(para_env%mepos==0) WRITE(*,*) "Imgxphase", AIMAG(xphase)

      ! charge
      trace = 0.0_dp
      DO ispin = 1, dft_control%nspins
         CALL dbcsr_trace(rho_ao(ispin)%matrix, matrix_s(1)%matrix, trace)
         rmom(1, 1) = rmom(1, 1)+trace
         newrmom(1, 1) = newrmom(1, 1)+trace
      END DO

      zi = 0._dp
      zij = 0._dp
      zijk = 0._dp
      zijkl = 0._dp

      IF(para_env%mepos==0) WRITE(*,*) "ELECTRONIC"
      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! OLD OLD OLD OLD OLD OLD
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
            ! OLD OLD OLD OLD OLD OLD

            ! Dipole
            DO i = 1, 3
               IF(para_env%mepos==0) WRITE(*,*) "DIPOLE DIRECTION i", i
               !ORGINIAL
               kvec(:) = twopi*cell%h_inv(i, :)
!               IF(para_env%mepos==0) WRITE(*,*) "kvec(:)", kvec(:)
               ! d/dt (a|e^(ikr)|b)
               CALL build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec, &
                                              reference, ref_point, nderivative = 1, &
                                               dkcos = dkcos(i)%matrix, dksin = dksin(i)%matrix)


!               ! ik-derivative
!               CALL build_berry_moment_matrix2(qs_env, cosmat, sinmat, kvec, nderivative = 1,&
!                                               dkcos = dkcos, dksin = dksin)

!               CALL numer_berry_moment_matrix(qs_env, kvec, nmoments=1, &
!                                              dkc=dkcos, dks=dksin)

               ! modified t-derivation
!               CALL build_berry_moment_matrix4(qs_env, cosmat, sinmat, kvec, nderivative = 1, &
!                                               dkcos = dkcos(i)%matrix, dksin = dksin(i)%matrix)

!               WRITE(*,*) "t-Derivation"
!               IF(para_env%mepos==0) WRITE(*,*) "cosmat"
!               CALL dbcsr_print(cosmat, matlab_format = .TRUE.)
!               IF(para_env%mepos==0) WRITE(*,*) "sinmat"
!               CALL dbcsr_print(sinmat, matlab_format = .TRUE.)
!               IF(para_env%mepos==0) WRITE(*,*) "dkcos"
!               CALL dbcsr_print(dkcos(i)%matrix, matlab_format = .TRUE.)
!               IF(para_env%mepos==0) WRITE(*,*) "dksin"
!               CALL dbcsr_print(dksin(i)%matrix, matlab_format = .TRUE.)
                             
!               ! d/dk (a|e^(ikr)|b)
!               CALL build_berry_moment_matrix4(qs_env, cosmat, sinmat, kvec, nderivative = 1, &
!                                               dkcos = dkcos(i)%matrix, dksin = dksin(i)%matrix)

!               WRITE(*,*) "t-recursion-derivation"                             
!               CALL dbcsr_print(cosmat, matlab_format = .TRUE.)
!               CALL dbcsr_print(sinmat, matlab_format = .TRUE.)
!               CALL dbcsr_print(dkcos(i)%matrix, matlab_format = .TRUE.)
!               CALL dbcsr_print(dksin(i)%matrix, matlab_format = .TRUE.)

               ! rotation onto MOs 
               CALL op_orbbas(cosmat, sinmat, mos, op_fm_set, opvec)
               CALL op_orbbas(dkcos(i)%matrix, dksin(i)%matrix, mos, dop_fm_set, opvec)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

            et1 = CMPLX(0.0_dp, 0.0_dp, dp)
            DO ispin = 1, dft_control%nspins

               CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=tmp_dim)
               DO idim = 1, tmp_dim

                  ! Make complex matrix of d01
                  eigrmat(ispin)%matrix%local_data(:, idim) = &
                     CMPLX(op_fm_set(1, ispin)%matrix%local_data(:, idim), &
                           -op_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
                  ! Make complex matrix of d11
                  deigrmat(ispin)%matrix%local_data(:, idim) = &
                     CMPLX(dop_fm_set(1, ispin)%matrix%local_data(:, idim), &
                           -dop_fm_set(2, ispin)%matrix%local_data(:, idim), dp)
               END DO
               
               ! Invert d01
               CALL cp_cfm_lu_invert(eigrmat(ispin)%matrix)

!               !-------------------------------------------
!               CALL cp_cfm_to_fm(msource = deigrmat(1)%matrix, &
!                                 mtargetr = dop_fm_set(1,1)%matrix, &
!                                 mtargeti = dop_fm_set(2,1)%matrix)
!               IF(para_env%mepos==0)   WRITE(*,*) "deriveigrmat - COS deriv"
!               CALL cp_fm_write_formatted(dop_fm_set(1,1)%matrix, 6)
!               IF(para_env%mepos==0) WRITE(*,*) "deriveigrmat - SIN deriv"
!               CALL cp_fm_write_formatted(dop_fm_set(2,1)%matrix, 6)
!               CALL cp_cfm_to_fm(msource = eigrmat(1)%matrix, &
!                                 mtargetr = op_fm_set(1,1)%matrix, &
!                                 mtargeti = op_fm_set(2,1)%matrix)
!               IF(para_env%mepos==0) WRITE(*,*) "inveigrmat - COS inv"
!               CALL cp_fm_write_formatted(op_fm_set(1,1)%matrix, 6)
!               IF(para_env%mepos==0) WRITE(*,*) "inveigrmat - SIN inv"
!               CALL cp_fm_write_formatted(op_fm_set(2,1)%matrix, 6)
!               !-------------------------------------------

               ! tr(d01^-1 * d11)
               CALL cp_cfm_trace(matrix_a = eigrmat(ispin)%matrix, &
                                 matrix_b = deigrmat(ispin)%matrix, &
                                 trace = et1)

               etr = z_one*et1

               ! Closed shell
               IF (dft_control%nspins == 1) THEN
                  etr = etr + et1
               ENDIF

!               !---------------------------------------------
!               CALL cp_cfm_get_info(cumulant_matrix(1)%matrix, ncol_local = ncol_local)
!
!               ! matrix_c = alpha * op1( matrix_a ) * op2( matrix_b ) + beta*matrix_c.
!!               CALL cp_cfm_gemm("N","N", ncol_local, ncol_local, ncol_local, &
!!                                z_one, eigrmat(1)%matrix, deigrmat(1)%matrix, &
!!                                z_zero, cumulant_matrix(1)%matrix)
!
!               CALL  cp_cfm_schur_product(matrix_a = eigrmat(1)%matrix,&
!                                          matrix_b = deigrmat(1)%matrix,&
!                                          matrix_c = cumulant_matrix(1)%matrix)
!
!               ! cfm -> fm
!               ! nmo x nmo 
!               CALL cp_cfm_to_fm(msource = cumulant_matrix(1)%matrix, &
!                                 mtargetr = op_fm_set(1,1)%matrix, &
!                                 mtargeti = op_fm_set(2,1)%matrix)
!   
!               ! fm -> dbcsr
!               CALL copy_fm_to_dbcsr(op_fm_set(1,1)%matrix, cosmat)
!               CALL copy_fm_to_dbcsr(op_fm_set(2,1)%matrix, sinmat)
!
!               ! calculate trace
!               CALL dbcsr_trace(cosmat, relem)
!               CALL dbcsr_trace(sinmat, ielem)
!               et1 =  CMPLX(relem, ielem, dp)
!   
!               etr = z_one*et1
!               
!               ! Closed shell
!               IF (dft_control%nspins == 1) THEN
!                  etr = etr + et1
!               ENDIF
!               !---------------------------------------------

            END DO
!            ctrace(i) = etr / CMPLX(0.0_dp, kvec(i), dp)
            ! t-derivative
            ctrace(i) = etr 
            ctrace2(i) = etr / CMPLX(0.0_dp, kvec(i), dp)

            ! ik-derivative
!            ctrace(i) = etr * CMPLX(0.0_dp, -1.0_dp*kvec(i), dp)

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!               ! get traces
!!               strace = 0.0_dp
!               DO ispin = 1, dft_control%nspins
!
!                  ! Trace deriv real
!                  CALL dbcsr_set(dkcos(ispin)%matrix, 0.0_dp)
!                  CALL copy_fm_to_dbcsr(dop_fm_set(1, ispin)%matrix, dkcos(ispin)%matrix)
!                  CALL dbcsr_trace(dkcos(ispin)%matrix, relem)
!
!                  ! Trace deriv img
!                  CALL cp_fm_scale(-1.0_dp, dop_fm_set(2, ispin)%matrix)
!                  CALL dbcsr_set(dksin(ispin)%matrix, 0.0_dp)
!                  CALL copy_fm_to_dbcsr(dop_fm_set(2, ispin)%matrix, dksin(ispin)%matrix)
!                  CALL dbcsr_trace(dksin(ispin)%matrix, ielem)
!                  det1 = CMPLX(relem, ielem, dp)
!!                  WRITE(*,*) "det1", det1
!
!                  ! real
!                  CALL dbcsr_set(cosmat, 0.0_dp)
!                  CALL copy_fm_to_dbcsr(op_fm_set(1, ispin)%matrix, cosmat)
!                  CALL dbcsr_trace(cosmat, relem)
!
!                  ! img
!                  CALL cp_fm_scale(-1.0_dp, op_fm_set(2, ispin)%matrix)
!                  CALL dbcsr_set(sinmat, 0.0_dp)
!                  CALL copy_fm_to_dbcsr(op_fm_set(2, ispin)%matrix, sinmat)
!                  CALL dbcsr_trace(sinmat, ielem)
!
!                  ! Tr(d11) / Tr(d01)
!                  et1 = CMPLX(relem, ielem, dp)
!!                  WRITE(*,*) "et1", et1
!               END DO
!
!               ctrace(i) = det1 / et1 / CMPLX(0.0_dp, 1.0_dp*kvec(i), dp)
!               !tracei(i) = AIMAG(det1 / et1)
!               WRITE(*,*) "TRACE", ctrace(i)
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!               ! element by element inversion
!               strace = 0.0_dp
!               DO ispin = 1, dft_control%nspins
!                  CALL cp_cfm_get_info(eigrmat(ispin)%matrix, ncol_local=ncol_local, nrow_local = nrow_local)
!                  DO irow = 1, nrow_local  
!                     DO icol = 1, ncol_local 
!                        relem = op_fm_set(1, ispin)%matrix%local_data(irow, icol)
!                        ielem = -1.0_dp*op_fm_set(2, ispin)%matrix%local_data(irow, icol)
!                        ! set cfm matrix element to its reciprocal
!                        IF(relem == 0.0_dp .AND. ielem == 0.0_dp) THEN
!                           eigrmat(ispin)%matrix%local_data(irow, icol) = CMPLX(0.0_dp, 0.0_dp, dp)
!                        ELSE
!                           eigrmat(ispin)%matrix%local_data(irow, icol) = &
!                              CMPLX(relem/(relem**2+ielem**2), -1.0_dp*ielem/(relem**2+ielem**2), dp)
!                        END IF
!                     END DO
!                  END DO
!                  ! Store first deriv as complex matrix
!                  CALL cp_fm_scale(-1.0_dp, dop_fm_set(2, ispin)%matrix)
!                  CALL cp_fm_to_cfm(msourcer = dop_fm_set(1, ispin)%matrix, &
!                                    msourcei = dop_fm_set(2, ispin)%matrix, &
!                                    mtarget = deigrmat(ispin)%matrix)
!                  ! Multiplying Hadamard 
!                  CALL cp_cfm_schur_product(deigrmat(ispin)%matrix, eigrmat(ispin)%matrix, &
!                                            cumulant_matrix(ispin)%matrix)
!                  ! cfm -> fm
!                  CALL cp_cfm_to_fm(msource = cumulant_matrix(ispin)%matrix, &
!                                    mtargetr = tmp_op_fm_set(1, ispin)%matrix, &
!                                    mtargeti = tmp_op_fm_set(2, ispin)%matrix)
!                  ! fm -> dbcsr
!                  CALL copy_fm_to_dbcsr(tmp_op_fm_set(1, ispin)%matrix, cosmat)
!                  CALL copy_fm_to_dbcsr(tmp_op_fm_set(2, ispin)%matrix, sinmat)
!                  
!                   
!                   ! dbcsr -> trace
!                  CALL dbcsr_trace(sinmat, sintrace)
!                  strace = strace + sintrace
!               END DO
!               ! save trace for i
!               tracei(i) = strace
!------------------------------------------------------------------------------
            END DO ! i
            IF(para_env%mepos==0) WRITE(*,*) "ctrace", ctrace
            IF(para_env%mepos==0) WRITE(*,*) "ctrace2", ctrace2
            ! Beware, missing multiplication with phase??
!            ctrace = ctrace*xphase

         CASE (2)
            ! Quadrupole
!            CPABORT("Berry phase moments bigger than 1 not implemented")
            DO i = 1, 3
               DO j = i, 3
                  kvec(:) = twopi*(cell%h_inv(i, :)+cell%h_inv(j, :))
                  ! (a|e^ikr|b)
                  CALL build_berry_moment_matrix3(qs_env, cosmat, sinmat, kvec, &
                                                 nderivative = 1, dkcos = dkcos(i)%matrix, dksin = dksin(i)%matrix)
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
                           ! e^ikr -> exp(-ikr) = cos(kr) - i*sin(kr)
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
!------------------------------------------------------------------------------
      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole (apply periodic (2 Pi) boundary conditions)
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
            ci = AIMAG(LOG(zi))
            IF(para_env%mepos==0) WRITE(*,*) "ci", ci
            DO i = 1, 3
               IF (qq(i)+ci(i) > pi) ci(i) = ci(i)-twopi
               IF (qq(i)+ci(i) < -pi) ci(i) = ci(i)+twopi
            END DO
            IF(para_env%mepos==0) WRITE(*,*)  "ci after pi", ci
            rmom(2:4, 1) = MATMUL(cell%hmat, ci)/twopi

           ! ik-derivative
!           tracei(:) = REAL(ctrace(:),dp)
           ! t-derivative
           tracei(:) = AIMAG(ctrace(:))
           tracei2(:) = REAL(ctrace2(:),dp)
           
           ! shift new elec based on old nuclear
           IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
           IF(para_env%mepos==0) WRITE(*,*) "Berry nuclear + Greg elec"
           IF(para_env%mepos==0) WRITE(*,*) "tracei", tracei
            DO i = 1, 3
               IF(para_env%mepos==0) WRITE(*,*) "START"
               IF(para_env%mepos==0) WRITE(*,*) "iteration", i
               IF(para_env%mepos==0) WRITE(*,*) "qq(i) =", qq(i)
               IF(para_env%mepos==0) WRITE(*,*) "qq(i)+tracei(i) =", qq(i)+tracei(i)
               DO WHILE  (qq(i)+tracei(i) > pi)
                  IF(para_env%mepos==0) WRITE(*,*) "First Do While", qq(i)+tracei(i)
                  IF(para_env%mepos==0) WRITE(*,*) "tracei(i)", tracei(i)
                  tracei(i) = tracei(i)-twopi
                  IF(para_env%mepos==0) WRITE(*,*) "qq(i)+tracei(i)", qq(i)+tracei(i)
               END DO
               DO WHILE (qq(i)+tracei(i) < -pi)
                  IF(para_env%mepos==0) WRITE(*,*) "Second Do While", qq(i)+tracei(i)
                  IF(para_env%mepos==0) WRITE(*,*) "tracei(i)", tracei(i)
                  tracei(i) = tracei(i)+twopi
                  IF(para_env%mepos==0) WRITE(*,*) "qq(i)+tracei(i)", qq(i)+tracei(i)
               END DO
               IF(para_env%mepos==0) WRITE(*,*) "criteria  after second DO WHILE", qq(i)+tracei(i)
            END DO
            IF(para_env%mepos==0) WRITE(*,*) "tracei AFTER", tracei
           
            ! Shift electronic based on new nuclear
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)
           IF(para_env%mepos==0) WRITE(*,*) "Greg nuclear + elec"
            DO i = 1, 3
               IF(para_env%mepos==0) WRITE(*,*) "iteration", i
               IF(para_env%mepos==0) WRITE(*,*) "newqq(i) =", newqq(i)
               IF(para_env%mepos==0) WRITE(*,*) "newqq(i)+tracei2(i) =", newqq(i)+tracei2(i)
               IF(para_env%mepos==0) WRITE(*,*) "START"
               DO WHILE  (newqq(i)+tracei2(i) > pi)
                  tracei2(i) = tracei2(i)-twopi
                  IF(para_env%mepos==0) WRITE(*,*) "newqq(i)+tracei(i)", newqq(i)+tracei2(i)
               END DO
               IF(para_env%mepos==0) WRITE(*,*) "tracei(i) after -2pi", newqq(i)+tracei2(i)
               DO WHILE (newqq(i)+tracei2(i) < -pi)
                  tracei2(i) = tracei2(i)+twopi
                  IF(para_env%mepos==0) WRITE(*,*) "newqq(i)+tracei(i)", newqq(i)+tracei2(i)
               END DO
               IF(para_env%mepos==0) WRITE(*,*) "criteria  after +2pi", newqq(i)+tracei2(i)
            END DO
            IF(para_env%mepos==0) WRITE(*,*) "tracei2 AFTER", tracei2
            IF(para_env%mepos==0) WRITE(*,*) REPEAT("-", 27)

            ! Derivation wrt t pulls down additional k (2pi/L)

            ! t & ik-derivative
!            newrmom(2:4, 1) = MATMUL(cell%hmat, tracei(:))/twopi
            newrmom(2:4, 1) = MATMUL(cell%hmat, tracei(:))/twopi
            newrmom2(2:4, 1) = tracei2(:)

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
            ! Hardcoding of upper triangular matrix elements of cij
            ! into moment vector
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
!------------------------------------------------------------------------------
      ! Berry phase
      rmom(:, 3) = rmom(:, 1)+rmom(:, 2)
      ! Greg elec and Berry nucl
      newrmom(:, 3) = newrmom(:, 1)+rmom(:, 2)
      ! Greg elec and nucl
      newrmom2(:, 3) = newrmom2(:, 1)+newrmom(:, 2)
      IF(para_env%mepos==0) THEN 
         WRITE(*,*) ""
         WRITE(*,*) "nuclear", rmom(2:4, 2)
         WRITE(*,*) "rmom", rmom(2:4, 1)
         WRITE(*,*) ""
         WRITE(*,*) "newnuclear", newrmom(2:4, 2)
         WRITE(*,*) "newrmom", newrmom(2:4, 1)
         WRITE(*,*) "newrmom2", newrmom2(2:4, 1)
         WRITE(*,*) ""
!         WRITE(*,*) "Electronic contribution OLD [Debye]"
!         WRITE(*,*) rmom(2:4, 1)*debye
!         WRITE(*,*) ""
!         WRITE(*,*) "Electronic contribution NEW [Debye]"
!         WRITE(*,*) newrmom(2:4, 1)*debye
!         WRITE(*,*) ""
         WRITE(*,*) "TOTAL_OLD", rmom(2:4, 2)+rmom(2:4, 1)
         WRITE(*,*) ""
         WRITE(*,*) "TOTAL_oNnE", rmom(2:4, 2)+newrmom(2:4, 1)
         WRITE(*,*) ""
         WRITE(*,*) "TOTAL_NEW", newrmom2(2:4, 3)
         WRITE(*,*) ""

         WRITE(*,*) "DIVIDE: old / oNnE TOTAL (Berry vs 1/2 Greg)"
         IF (newrmom(2,1) .NE. 0) THEN
            WRITE(*,*) rmom(2,3) / newrmom(2, 3)
            WRITE(*,*) "Error X [%]", 100.0_dp*(rmom(2,3) - newrmom(2, 3)) / newrmom(2, 3)
         ELSE
            WRITE(*,*) 0.0_dp
         END IF
         IF (newrmom(3,1) .NE. 0) THEN
            WRITE(*,*) rmom(3,3) / newrmom(3,3)
            WRITE(*,*) "Error Y [%]", 100.0_dp*(rmom(3,3) - newrmom(3, 3)) / newrmom(3, 3)
         ELSE
            WRITE(*,*) 0.0_dp
         END IF
         IF (newrmom(4,1) .NE. 0) THEN
            WRITE(*,*) rmom(4,3) / newrmom(4,3)
            WRITE(*,*) "Error Z [%]", 100.0_dp*(rmom(4,3) - newrmom(4, 3)) / newrmom(4, 3)
         ELSE
            WRITE(*,*) 0.0_dp
         END IF

         WRITE(*,*) "DIVIDE: old / new  TOTAL (Berry vs Greg)"
         IF (newrmom(2,1) .NE. 0) THEN
            WRITE(*,*) rmom(2,3) / newrmom2(2, 3)
            WRITE(*,*) "Error X [%]", 100.0_dp*(rmom(2,3) - newrmom2(2, 3)) / newrmom2(2, 3)
         ELSE
            WRITE(*,*) 0.0_dp
         END IF
         IF (newrmom(3,1) .NE. 0) THEN
            WRITE(*,*) rmom(3,3) / newrmom2(3, 3)
            WRITE(*,*) "Error Y [%]", 100.0_dp*(rmom(3,3) - newrmom2(3, 3)) / newrmom2(3, 3)
         ELSE
            WRITE(*,*) 0.0_dp
         END IF
         IF (newrmom(4,1) .NE. 0) THEN
            WRITE(*,*) rmom(4,3) / newrmom2(4, 3)
            WRITE(*,*) "Error Z [%]", 100.0_dp*(rmom(4,3) - newrmom2(4, 3)) / newrmom2(4, 3) 
         ELSE
            WRITE(*,*) 0.0_dp
         END IF
      END IF

      description = "[DIPOLE]"
      CALL cp_results_erase(results=results, description=description)
      CALL put_results(results=results, description=description, &
                       values=rmom(2:4, 3))
      IF (magnetic) THEN
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.TRUE., mmom=mmom)
      ELSE
         IF(para_env%mepos==0) WRITE(*,*) "Berry phase"
         CALL print_moments(unit_number, nmom, rmom, rlab, rcc, cell, periodic=.TRUE.)
         IF(para_env%mepos==0) WRITE(*,*) "Berry nuclear + Greg elec"
         CALL print_moments(unit_number, nmom, newrmom, rlab, rcc, cell, periodic=.TRUE.)
         IF(para_env%mepos==0) WRITE(*,*) "Greg nuclear + elec"
         CALL print_moments(unit_number, nmom, newrmom2, rlab, rcc, cell, periodic=.TRUE.)
      END IF

      DEALLOCATE (rmom)
      DEALLOCATE (newrmom)
      DEALLOCATE (newrmom2)
      DEALLOCATE (rlab)
      IF (magnetic) THEN
         DEALLOCATE (mmom)
      END IF

      CALL dbcsr_deallocate_matrix(cosmat)
      CALL dbcsr_deallocate_matrix(sinmat)
      CALL dbcsr_deallocate_matrix_set(dkcos)
      CALL dbcsr_deallocate_matrix_set(dksin)

      DO ispin = 1, dft_control%nspins
         CALL cp_fm_release(opvec(ispin)%matrix)
         CALL cp_cfm_release(eigrmat(ispin)%matrix)
         CALL cp_cfm_release(deigrmat(ispin)%matrix)
         CALL cp_cfm_release(cumulant_matrix(ispin)%matrix)

         DO i = 1, SIZE(op_fm_set, 1)
            CALL cp_fm_release(op_fm_set(i, ispin)%matrix)
            CALL cp_fm_release(tmp_op_fm_set(i, ispin)%matrix)
            CALL cp_fm_release(dop_fm_set(i, ispin)%matrix)
         END DO
      END DO

      DEALLOCATE (op_fm_set)
      DEALLOCATE (dop_fm_set)
      DEALLOCATE (tmp_op_fm_set)
      DEALLOCATE (opvec)
      DEALLOCATE (eigrmat)
      DEALLOCATE (deigrmat)
      DEALLOCATE (cumulant_matrix)

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

      CPASSERT(ASSOCIATED(qs_env))

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
      ! rmom(:,1)=total
      ALLOCATE (rmom(nm+1, 3))
      ALLOCATE (rlab(nm+1))
      rmom = 0.0_dp
      rlab = ""
      IF (magnetic) THEN
         nm = 3
         ALLOCATE (mmom(nm))
         mmom = 0._dp
      END IF

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

      zi = 0._dp
      zij = 0._dp
      zijk = 0._dp
      zijkl = 0._dp

      DO l = 1, nmom
         SELECT CASE (l)
         CASE (1)
            ! Dipole
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
                  END DO
               ENDIF
            END DO
            zi = zi*zphase
            ci = AIMAG(LOG(zi))/twopi
            qq = AIMAG(LOG(zi))
            rmom(2:4, 2) = MATMUL(cell%hmat, ci)
         CASE (2)
            ! Quadrupole
            CPABORT("Berry phase moments bigger than 1 not implemented")
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
                     zij(j, i) = zij(i, j)
                  END DO
               END DO
            END DO
            DO i = 1, 3
               DO j = 1, 3
                  zij(i, j) = zij(i, j)*zphase(i)*zphase(j)
                  zz = zij(i, j)/zi(i)/zi(j)
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
            END DO
            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))
            DO k = 4, 9
               ix = indco(1, k+1)
               iy = indco(2, k+1)
               iz = indco(3, k+1)
               IF (ix == 0) THEN
                  rmom(k+1, 2) = cij(iy, iz)
               ELSE IF (iy == 0) THEN
                  rmom(k+1, 2) = cij(ix, iz)
               ELSE IF (iz == 0) THEN
                  rmom(k+1, 2) = cij(ix, iy)
               END IF
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
            CPABORT("Berry phase moments bigger than 1 not implemented")
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
            CPABORT("Berry phase moments bigger than 1 not implemented")
            DO i = 1, 3
               DO j = 1, 3
                  zz = zij(i, j)/zi(i)/zi(j)
                  cij(i, j) = AIMAG(LOG(zz))/twopi
               END DO
            END DO
            cij = 0.5_dp*cij/twopi/twopi
            cij = MATMUL(MATMUL(cell%hmat, cij), TRANSPOSE(cell%hmat))
            DO k = 4, 9
               ix = indco(1, k+1)
               iy = indco(2, k+1)
               iz = indco(3, k+1)
               IF (ix == 0) THEN
                  rmom(k+1, 1) = cij(iy, iz)
               ELSE IF (iy == 0) THEN
                  rmom(k+1, 1) = cij(ix, iz)
               ELSE IF (iz == 0) THEN
                  rmom(k+1, 1) = cij(ix, iy)
               END IF
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

      TYPE(dbcsr_type),                                  :: cosmat, sinmat
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
      TYPE(dbcsr_type), POINTER                          :: rho_magmom_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho

      CPASSERT(ASSOCIATED(qs_env))

      CALL timeset(routineN, handle)

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

      ! magnetic moments
      IF (magnetic) THEN
         NULLIFY (magmom)
         CALL build_local_magmom_matrix(qs_env, magmom, nmom, ref_point=rcc)
         nm = SIZE(magmom)
         ALLOCATE (mmom(nm))

         ! Allocate matrices to store the matrix product to be traced (dbcsr_trace only works for products of
         ! symmetric matrices)
         NULLIFY (rho_magmom_ao)
         ALLOCATE (rho_magmom_ao)
         CALL dbcsr_desymmetrize(matrix_s(1)%matrix, rho_magmom_ao)

         IF (qs_env%run_rtp) THEN
            ! get imaginary part of the density in rho_ao (the real part is not needed since the trace of the product
            ! of a symmetric (REAL(rho_ao)) and an anti-symmetric (L_AO) matrix is zero)
            ! There may be other cases, where the imaginary part of the density is relevant
            NULLIFY (rho_ao)
            CALL qs_rho_get(rho, rho_ao_im=rho_ao)
         ENDIF
         ! if the density is purely real this is an expensive way to calculate zero
         DO i = 1, SIZE(magmom)
            strace = 0._dp
            DO ispin = 1, dft_control%nspins
               CALL dbcsr_set(rho_magmom_ao, 0.0_dp)
               CALL dbcsr_multiply("T", "N", 1.0_dp, rho_ao(ispin)%matrix, magmom(i)%matrix, &
                                   0.0_dp, rho_magmom_ao)
               CALL dbcsr_trace(rho_magmom_ao, trace)
               strace = strace+trace
            END DO
            mmom(i) = strace
         END DO

         CALL dbcsr_deallocate_matrix_set(magmom)
         CALL dbcsr_deallocate_matrix(rho_magmom_ao)
      END IF

      description = "[DIPOLE]"
      CALL cp_results_erase(results=results, description=description)
      CALL put_results(results=results, description=description, &
                       values=rmom(2:4, 3))
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
               WRITE (unit_number, "(T3,A)") "Magnetic Dipole Moment (only orbital contrib.) [a.u.]"
               WRITE (unit_number, "(T5,3(A,A,F14.8,1X),T60,A,T68,F14.8)") &
                  (TRIM(rlab(i+1)), "=", mmom(i), i=1, 3), "Total=", dd
            END IF
         END IF
      END IF

   END SUBROUTINE print_moments

END MODULE qs_moments

