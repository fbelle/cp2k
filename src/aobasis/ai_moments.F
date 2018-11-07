!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the moment integrals over Cartesian Gaussian-type
!>      functions.
!> \par Literature
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      none
!> \author J. Hutter (16.02.2005)
! **************************************************************************************************
MODULE ai_moments

! ax,ay,az  : Angular momentum index numbers of orbital a.
! bx,by,bz  : Angular momentum index numbers of orbital b.
! coset     : Cartesian orbital set pointer.
! dab       : Distance between the atomic centers a and b.
! l{a,b}    : Angular momentum quantum number of shell a or b.
! l{a,b}_max: Maximum angular momentum quantum number of shell a or b.
! l{a,b}_min: Minimum angular momentum quantum number of shell a or b.
! rac       : Distance vector between the atomic center a and reference point c.
! rbc       : Distance vector between the atomic center b and reference point c.
! rpgf{a,b} : Radius of the primitive Gaussian-type function a or b.
! zet{a,b}  : Exponents of the Gaussian-type functions a or b.
! zetp      : Reciprocal of the sum of the exponents of orbital a and b.

   USE ai_derivatives,                  ONLY: adbdr,&
                                              dabdr
   USE ai_overlap,                      ONLY: overlap
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              ncoset
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_moments'

   PUBLIC :: cossin, moment, diffop, diff_momop, contract_cossin, dipole_force
   PUBLIC :: cossin2
   PUBLIC :: cossin3
   PUBLIC :: cossin4

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cos_block ...
!> \param sin_block ...
!> \param iatom ...
!> \param ncoa ...
!> \param nsgfa ...
!> \param sgfa ...
!> \param sphi_a ...
!> \param ldsa ...
!> \param jatom ...
!> \param ncob ...
!> \param nsgfb ...
!> \param sgfb ...
!> \param sphi_b ...
!> \param ldsb ...
!> \param cosab ...
!> \param sinab ...
!> \param ldab ...
!> \param work ...
!> \param ldwork ...
! **************************************************************************************************
   SUBROUTINE contract_cossin(cos_block, sin_block, &
                              iatom, ncoa, nsgfa, sgfa, sphi_a, ldsa, &
                              jatom, ncob, nsgfb, sgfb, sphi_b, ldsb, &
                              cosab, sinab, ldab, work, ldwork)

      REAL(dp), DIMENSION(:, :), POINTER                 :: cos_block, sin_block
      INTEGER, INTENT(IN)                                :: iatom, ncoa, nsgfa, sgfa
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: sphi_a
      INTEGER, INTENT(IN)                                :: ldsa, jatom, ncob, nsgfb, sgfb
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: sphi_b
      INTEGER, INTENT(IN)                                :: ldsb
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: cosab, sinab
      INTEGER, INTENT(IN)                                :: ldab
      REAL(dp), DIMENSION(:, :)                          :: work
      INTEGER, INTENT(IN)                                :: ldwork

! Calculate cosine

      CALL dgemm("N", "N", ncoa, nsgfb, ncob, &
                 1.0_dp, cosab(1, 1), ldab, &
                 sphi_b(1, sgfb), ldsb, &
                 0.0_dp, work(1, 1), ldwork)

      IF (iatom <= jatom) THEN
         CALL dgemm("T", "N", nsgfa, nsgfb, ncoa, &
                    1.0_dp, sphi_a(1, sgfa), ldsa, &
                    work(1, 1), ldwork, &
                    1.0_dp, cos_block(sgfa, sgfb), &
                    SIZE(cos_block, 1))
      ELSE
         CALL dgemm("T", "N", nsgfb, nsgfa, ncoa, &
                    1.0_dp, work(1, 1), ldwork, &
                    sphi_a(1, sgfa), ldsa, &
                    1.0_dp, cos_block(sgfb, sgfa), &
                    SIZE(cos_block, 1))
      END IF

      ! Calculate sine
      CALL dgemm("N", "N", ncoa, nsgfb, ncob, &
                 1.0_dp, sinab(1, 1), ldab, &
                 sphi_b(1, sgfb), ldsb, &
                 0.0_dp, work(1, 1), ldwork)

      IF (iatom <= jatom) THEN
         CALL dgemm("T", "N", nsgfa, nsgfb, ncoa, &
                    1.0_dp, sphi_a(1, sgfa), ldsa, &
                    work(1, 1), ldwork, &
                    1.0_dp, sin_block(sgfa, sgfb), &
                    SIZE(sin_block, 1))
      ELSE
         CALL dgemm("T", "N", nsgfb, nsgfa, ncoa, &
                    1.0_dp, work(1, 1), ldwork, &
                    sphi_a(1, sgfa), ldsa, &
                    1.0_dp, sin_block(sgfb, sgfa), &
                    SIZE(sin_block, 1))
      END IF

   END SUBROUTINE contract_cossin


! **************************************************************************************************
!                           Chris explained Greg derivation, with old recursion derivation 
!> \brief ...                  !!!!!! need new recursion derivation WRT t and then index
!> \param la_max_set ...
!> \param npgfa ...         my deriv -> 3 , Gregs deriv -> 1 , need my structure for deriv, but kvec(i) survives
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min_set ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param rac ...
!> \param rbc ...
!> \param kvec ...
!> \param cosab ...
!> \param sinab ...
!> \param dcosab ...
!> \param dsinab ...
!> \param dsc  First cosine derivative wrt k  
!> \param dss  First sine derivative wrt k 
!> \param ddsc Second cosine derivative wrt k
!> \param ddss Second sine derivative wrt k
! **************************************************************************************************
   SUBROUTINE cossin4(la_max_set, npgfa, zeta, rpgfa, la_min_set, &
                     lb_max, npgfb, zetb, rpgfb, lb_min, &
                     rac, rbc, kvec, cosab, sinab, dcosab, dsinab, &
                     dksc, dkss, ddksc, ddkss)

      INTEGER, INTENT(IN)                                :: la_max_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min_set, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc, kvec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: cosab, sinab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dcosab, dsinab
      !fbelle
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
                  OPTIONAL                               :: dksc, dkss, ddksc, ddkss

      INTEGER :: ax, ay, az, bx, by, bz, cda, cdax, cday, cdaz, coa, coamx, coamy, coamz, coapx, &
         coapy, coapz, cob, da, da_max, dax, day, daz, i, ipgf, j, jpgf, k, la, la_max, la_min, &
         la_start, lb, lb_start, na, nb
      REAL(KIND=dp)                                      :: dab, f0, f1, f2, f3, fax, fay, faz, ftz, &
                                                            fx, fy, fz, k2, kdp, rab2, s, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rap, rbp
      REAL(KIND=dp), DIMENSION(ncoset(la_max_set), &
         ncoset(lb_max), 3)                              :: dscos, dssin
      REAL(KIND=dp), &
         DIMENSION(ncoset(la_max_set+1), ncoset(lb_max)) :: sc, ss
      !fbelle
      REAL(KIND=dp), DIMENSION(3)                        :: r 
      REAL(KIND=dp), DIMENSION(ncoset(la_max_set+1), ncoset(lb_max), 3) :: dsc, dss
      REAL(KIND=dp), DIMENSION(ncoset(la_max_set+1), ncoset(lb_max), 6) :: ddsc, ddss
      REAL(KIND=dp)                                      :: s1, s2, s3, cos1, cos2, cos3, sin1, sin2, sin3


!      WRITE(*,*) "cossin2 - BEGINNING"

      rab = rbc-rac
      rab2 = SUM(rab**2)
      dab = SQRT(rab2)
      k2 = kvec(1)*kvec(1)+kvec(2)*kvec(2)+kvec(3)*kvec(3)

      IF (PRESENT(dcosab)) THEN
         da_max = 1
         la_max = la_max_set+1
         la_min = MAX(0, la_min_set-1)
         dscos = 0.0_dp
         dssin = 0.0_dp
      ELSE
         da_max = 0
         la_max = la_max_set
         la_min = la_min_set
      END IF

      ! initialize all matrix elements to zero
      IF (PRESENT(dcosab)) THEN
         na = ncoset(la_max-1)*npgfa
      ELSE
         na = ncoset(la_max)*npgfa
      END IF
      nb = ncoset(lb_max)*npgfb
      cosab(1:na, 1:nb) = 0.0_dp
      sinab(1:na, 1:nb) = 0.0_dp
      IF (PRESENT(dcosab)) THEN
         dcosab(1:na, 1:nb, :) = 0.0_dp
         dsinab(1:na, 1:nb, :) = 0.0_dp
      END IF

      ! first derivative 
      IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
         DO i = 1, 3
            dksc(1:na, 1:nb, i) = 0.0_dp
            dkss(1:na, 1:nb, i) = 0.0_dp
         END DO
      END IF
      ! second derivative
      IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
         DO i = 1, 6
            ddksc(1:na, 1:nb, i) = 0.0_dp
            ddkss(1:na, 1:nb, i) = 0.0_dp
         END DO
      END IF

!   *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0
      DO ipgf = 1, npgfa

         nb = 0

         DO jpgf = 1, npgfb

            ss = 0.0_dp
            sc = 0.0_dp

            IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
               dsc = 0.0_dp
               dss = 0.0_dp
            END IF
            IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
               ddsc = 0.0_dp
               ddss = 0.0_dp
            END IF

!       *** Screening ***
            IF (rpgfa(ipgf)+rpgfb(jpgf) < dab) THEN
               nb = nb+ncoset(lb_max)
               CYCLE
            END IF

!       *** Calculate some prefactors ***

            zetp = 1.0_dp/(zeta(ipgf)+zetb(jpgf))

            f0 = (pi*zetp)**1.5_dp
            f1 = zetb(jpgf)*zetp
            f2 = 0.5_dp*zetp

            kdp = zetp*DOT_PRODUCT(kvec, zeta(ipgf)*rac+zetb(jpgf)*rbc)

!            r = zetp*(zeta(ipgf)*rac+zetb(jpgf)*rbc)
            
            ! Pulling down k_i*r_i instead of r_i
            r = kvec*zetp*(zeta(ipgf)*rac+zetb(jpgf)*rbc)

!       *** Calculate the basic two-center cos/sin integral [s|cos/sin|s] ***

            s = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*k2*zetp)
            sc(1, 1) = s*COS(kdp)
            ss(1, 1) = s*SIN(kdp)

            !----------------




            !----------------

            !----------------
            s1 = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*kvec(1)*kvec(1)*zetp)
            s2 = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*kvec(2)*kvec(2)*zetp)
            s3 = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*kvec(3)*kvec(3)*zetp)

            cos1 = COS(zetp*(kvec(1)*(zeta(ipgf)*rac(1)+zetb(jpgf)*rbc(1))))
            cos2 = COS(zetp*(kvec(2)*(zeta(ipgf)*rac(2)+zetb(jpgf)*rbc(2))))
            cos3 = COS(zetp*(kvec(3)*(zeta(ipgf)*rac(3)+zetb(jpgf)*rbc(3))))
            sin1 = SIN(zetp*(kvec(1)*(zeta(ipgf)*rac(1)+zetb(jpgf)*rbc(1))))
            sin2 = SIN(zetp*(kvec(2)*(zeta(ipgf)*rac(2)+zetb(jpgf)*rbc(2))))
            sin3 = SIN(zetp*(kvec(3)*(zeta(ipgf)*rac(3)+zetb(jpgf)*rbc(3))))

            IF (PRESENT(dksc)) THEN
               dsc(1,1,1) = -f2*kvec(1)*kvec(1)*s1*cos1 - r(1)*s1*cos1
               dsc(1,1,2) = -f2*kvec(2)*kvec(2)*s2*cos2 - r(2)*s2*cos2 
               dsc(1,1,3) = -f2*kvec(3)*kvec(3)*s3*cos3 - r(3)*s3*cos3
            END IF
            IF (PRESENT(dkss)) THEN
               dss(1,1,1) = -f2*kvec(1)*kvec(1)*s1*sin1 + r(1)*s1*sin1
               dss(1,1,2) = -f2*kvec(2)*kvec(2)*s2*sin2 + r(2)*s2*sin2
               dss(1,1,3) = -f2*kvec(3)*kvec(3)*s3*sin3 + r(3)*s3*sin3
            END IF
            !----------------

!Fabian            
!            IF (PRESENT(dksc)) THEN
!               ! dsc/dk1
!               dsc(1,1,1) = -f2*kvec(1)*sc(1,1) - r(1)*ss(1,1) 
!               ! dsc/dk2
!               dsc(1,1,2) = -f2*kvec(2)*sc(1,1) - r(2)*ss(1,1) 
!               ! dsc/dk3
!               dsc(1,1,3) = -f2*kvec(3)*sc(1,1) - r(3)*ss(1,1)
!            END IF
!            IF (PRESENT(dkss)) THEN
!               ! dss/dk1
!               dss(1,1,1) = -f2*kvec(1)*ss(1,1) + r(1)*sc(1,1) 
!               ! dss/dk2
!               dss(1,1,2) = -f2*kvec(2)*ss(1,1) + r(2)*sc(1,1) 
!               ! dss/dk3
!               dss(1,1,3) = -f2*kvec(3)*ss(1,1) + r(3)*sc(1,1) 
!            END IF

!       *** Recurrence steps: [s|O|s] -> [a|O|b] ***

            IF (la_max > 0) THEN
!         *** Vertical recurrence steps: [s|O|s] -> [a|O|s] ***

               rap(:) = f1*rab(:)

!         *** [p|O|s] = (Pi - Ai)*[s|O|s] +[s|dO|s]  (i = x,y,z) ***

               sc(2, 1) = rap(1)*sc(1, 1)-f2*kvec(1)*ss(1, 1)
               sc(3, 1) = rap(2)*sc(1, 1)-f2*kvec(2)*ss(1, 1)
               sc(4, 1) = rap(3)*sc(1, 1)-f2*kvec(3)*ss(1, 1)
               ss(2, 1) = rap(1)*ss(1, 1)+f2*kvec(1)*sc(1, 1)
               ss(3, 1) = rap(2)*ss(1, 1)+f2*kvec(2)*sc(1, 1)
               ss(4, 1) = rap(3)*ss(1, 1)+f2*kvec(3)*sc(1, 1)

               IF (PRESENT(dksc)) THEN
                  dsc(2,1,1) = rap(1)*dsc(1,1,1) - f2*kvec(1)*ss(1,1) - f2*kvec(1)*dss(1,1,1) 
                  dsc(2,1,2) = rap(1)*dsc(1,1,2) - f2*kvec(1)*dss(1,1,2)
                  dsc(2,1,3) = rap(1)*dsc(1,1,3) - f2*kvec(1)*dss(1,1,3)

                  dsc(3,1,1) = rap(2)*dsc(1,1,1) - f2*kvec(2)*dss(1,1,1) 
                  dsc(3,1,2) = rap(2)*dsc(1,1,2) - f2*kvec(2)*ss(1,1) - f2*kvec(2)*dss(1,1,2)
                  dsc(3,1,3) = rap(2)*dsc(1,1,3) - f2*kvec(2)*dss(1,1,3)

                  dsc(4,1,1) = rap(3)*dsc(1,1,1) - f2*kvec(3)*dss(1,1,1) 
                  dsc(4,1,2) = rap(3)*dsc(1,1,2) - f2*kvec(3)*dss(1,1,2)
                  dsc(4,1,3) = rap(3)*dsc(1,1,3) - f2*kvec(3)*ss(1,1) - f2*kvec(3)*dss(1,1,3)
               END IF

               IF (PRESENT(dkss)) THEN
                  dss(2,1,1) = rap(1)*dss(1,1,1) + f2*kvec(1)*sc(1,1) + f2*kvec(1)*dsc(1,1,1) 
                  dss(2,1,2) = rap(1)*dss(1,1,2) + f2*kvec(1)*dsc(1,1,2)
                  dss(2,1,3) = rap(1)*dss(1,1,3) + f2*kvec(1)*dsc(1,1,3)

                  dss(3,1,1) = rap(2)*dss(1,1,1) + f2*kvec(2)*dsc(1,1,1) 
                  dss(3,1,2) = rap(2)*dss(1,1,2) + f2*kvec(2)*sc(1,1) + f2*kvec(2)*dsc(1,1,2)
                  dss(3,1,3) = rap(2)*dss(1,1,3) + f2*kvec(2)*dsc(1,1,3)

                  dss(4,1,1) = rap(3)*dss(1,1,1) + f2*kvec(3)*dsc(1,1,1) 
                  dss(4,1,2) = rap(3)*dss(1,1,2) + f2*kvec(3)*dsc(1,1,2)
                  dss(4,1,3) = rap(3)*dss(1,1,3) + f2*kvec(3)*sc(1,1) + f2*kvec(3)*dsc(1,1,3)
               END IF

!         *** [a|O|s] = (Pi - Ai)*[a-1i|O|s] + f2*Ni(a-1i)*[a-2i|s] ***
!         ***           + [a-1i|dO|s]                               ***

               DO la = 2, la_max

!           *** Increase the angular momentum component z of function a ***

                  sc(coset(0, 0, la), 1) = rap(3)*sc(coset(0, 0, la-1), 1)+ &
                                           f2*REAL(la-1, dp)*sc(coset(0, 0, la-2), 1)- &
                                           f2*kvec(3)*ss(coset(0, 0, la-1), 1)
                  ss(coset(0, 0, la), 1) = rap(3)*ss(coset(0, 0, la-1), 1)+ &
                                           f2*REAL(la-1, dp)*ss(coset(0, 0, la-2), 1)+ &
                                           f2*kvec(3)*sc(coset(0, 0, la-1), 1)
           
                  IF (PRESENT(dksc)) THEN
                     dsc(coset(0, 0, la), 1, 1) = rap(3)*dsc(coset(0, 0, la-1), 1, 1) + &
                                                  f2*REAL(la-1, dp)*dsc(coset(0, 0, la-2), 1, 1)- &
                                                  f2*kvec(3)*dss(coset(0, 0, la-1), 1, 1)
                     dsc(coset(0, 0, la), 1, 2) = rap(3)*dsc(coset(0, 0, la-1), 1, 2)+ &
                                                  f2*REAL(la-1, dp)*dsc(coset(0, 0, la-2), 1, 2)- &
                                                  f2*kvec(3)*dss(coset(0, 0, la-1), 1, 2)
                     dsc(coset(0, 0, la), 1, 3) = rap(3)*dsc(coset(0, 0, la-1), 1, 3)+ &
                                                  f2*REAL(la-1, dp)*dsc(coset(0, 0, la-2), 1, 3)- &
                                                  f2*kvec(3)*ss(coset(0, 0, la-1), 1)- &
                                                  f2*kvec(3)*dss(coset(0, 0, la-1), 1, 3)
                  END IF
                  IF (PRESENT(dkss)) THEN
                     dss(coset(0, 0, la), 1, 1) = rap(3)*dss(coset(0, 0, la-1), 1, 1)+ &
                                                  f2*REAL(la-1, dp)*dss(coset(0, 0, la-2), 1, 1)+ &
                                                  f2*kvec(3)*dsc(coset(0, 0, la-1), 1, 1)
                     dss(coset(0, 0, la), 1, 2) = rap(3)*dss(coset(0, 0, la-1), 1, 2)+ &
                                                  f2*REAL(la-1, dp)*dss(coset(0, 0, la-2), 1, 2)+ &
                                                  f2*kvec(3)*dsc(coset(0, 0, la-1), 1, 2)
                     dss(coset(0, 0, la), 1, 3) = rap(3)*dss(coset(0, 0, la-1), 1, 3)+ &
                                                  f2*REAL(la-1, dp)*dss(coset(0, 0, la-2), 1, 3)+ &
                                                  f2*kvec(3)*sc(coset(0, 0, la-1), 1)+ &
                                                  f2*kvec(3)*dsc(coset(0, 0, la-1), 1, 3)
                  END IF

!           *** Increase the angular momentum component y of function a ***

                  az = la-1
                  sc(coset(0, 1, az), 1) = rap(2)*sc(coset(0, 0, az), 1)- &
                                           f2*kvec(2)*ss(coset(0, 0, az), 1)
                  ss(coset(0, 1, az), 1) = rap(2)*ss(coset(0, 0, az), 1)+ &
                                           f2*kvec(2)*sc(coset(0, 0, az), 1)

                  IF (PRESENT(dksc)) THEN
                     dsc(coset(0, 1, az), 1, 1) = rap(2)*dsc(coset(0, 0, az), 1, 1)- &
                                                  f2*kvec(2)*dss(coset(0, 0, az), 1, 1) 
                     dsc(coset(0, 1, az), 1, 2) = rap(2)*dsc(coset(0, 0, az), 1, 2)- &
                                                  f2*kvec(2)*ss(coset(0, 0, az), 1)- &
                                                  f2*kvec(2)*dss(coset(0, 0, az), 1, 2) 
                     dsc(coset(0, 1, az), 1, 3) = rap(2)*dsc(coset(0, 0, az), 1, 3)- &
                                                  f2*kvec(2)*dss(coset(0, 0, az), 1, 3) 
                  END IF
                  IF (PRESENT(dkss)) THEN
                     dss(coset(0, 1, az), 1, 1) = rap(2)*dss(coset(0, 0, az), 1, 1)+ &
                                                  f2*kvec(2)*dsc(coset(0, 0, az), 1, 1) 
                     dss(coset(0, 1, az), 1, 2) = rap(2)*dss(coset(0, 0, az), 1, 2)+ &
                                                  f2*kvec(2)*sc(coset(0, 0, az), 1)+ &
                                                  f2*kvec(2)*dsc(coset(0, 0, az), 1, 2) 
                     dss(coset(0, 1, az), 1, 3) = rap(2)*dss(coset(0, 0, az), 1, 3)+ &
                                                  f2*kvec(2)*dsc(coset(0, 0, az), 1, 3) 
                  END IF

                  DO ay = 2, la
                     az = la-ay
                     sc(coset(0, ay, az), 1) = rap(2)*sc(coset(0, ay-1, az), 1)+ &
                                               f2*REAL(ay-1, dp)*sc(coset(0, ay-2, az), 1)- &
                                               f2*kvec(2)*ss(coset(0, ay-1, az), 1)
                     ss(coset(0, ay, az), 1) = rap(2)*ss(coset(0, ay-1, az), 1)+ &
                                               f2*REAL(ay-1, dp)*ss(coset(0, ay-2, az), 1)+ &
                                               f2*kvec(2)*sc(coset(0, ay-1, az), 1)
        
                     IF (PRESENT(dksc)) THEN
                        dsc(coset(0, ay, az), 1, 1) = rap(2)*dsc(coset(0, ay-1, az), 1, 1)+ &
                                                      f2*REAL(ay-1, dp)*dsc(coset(0, ay-2, az), 1, 1)- &
                                                      f2*kvec(2)*dss(coset(0, ay-1, az), 1, 1) 
                        dsc(coset(0, ay, az), 1, 2) = rap(2)*dsc(coset(0, ay-1, az), 1, 2)+ &
                                                      f2*REAL(ay-1, dp)*dsc(coset(0, ay-2, az), 1, 2)- &
                                                      f2*kvec(2)*ss(coset(0, ay-1, az), 1)- &
                                                      f2*kvec(2)*dss(coset(0, ay-1, az), 1, 2) 
                        dsc(coset(0, ay, az), 1, 3) = rap(2)*dsc(coset(0, ay-1, az), 1, 3)+ &
                                                      f2*REAL(ay-1, dp)*dsc(coset(0, ay-2, az), 1, 3)- &
                                                      f2*kvec(2)*dss(coset(0, ay-1, az), 1, 3) 
                     END IF
                     IF (PRESENT(dkss)) THEN
                        dss(coset(0, ay, az), 1, 1) = rap(2)*dss(coset(0, ay-1, az), 1, 1)+ &
                                                      f2*REAL(ay-1, dp)*dss(coset(0, ay-2, az), 1, 1)+ &
                                                      f2*kvec(2)*dsc(coset(0, ay-1, az), 1, 1) 
                        dss(coset(0, ay, az), 1, 2) = rap(2)*dss(coset(0, ay-1, az), 1, 2)+ &
                                                      f2*REAL(ay-1, dp)*dss(coset(0, ay-2, az), 1, 2)+ &
                                                      f2*kvec(2)*sc(coset(0, ay-1, az), 1)+ &
                                                      f2*kvec(2)*dsc(coset(0, ay-1, az), 1, 2) 
                        dss(coset(0, ay, az), 1, 3) = rap(2)*dss(coset(0, ay-1, az), 1, 3)+ &
                                                      f2*REAL(ay-1, dp)*dss(coset(0, ay-2, az), 1, 3)+ &
                                                      f2*kvec(2)*dsc(coset(0, ay-1, az), 1, 3) 
                     END IF
                  END DO

!           *** Increase the angular momentum component x of function a ***

                  DO ay = 0, la-1
                     az = la-1-ay
                     sc(coset(1, ay, az), 1) = rap(1)*sc(coset(0, ay, az), 1)- &
                                               f2*kvec(1)*ss(coset(0, ay, az), 1)
                     ss(coset(1, ay, az), 1) = rap(1)*ss(coset(0, ay, az), 1)+ &
                                               f2*kvec(1)*sc(coset(0, ay, az), 1)
                  
                     IF (PRESENT(dksc)) THEN
                        dsc(coset(1, ay, az), 1, 1) = rap(1)*dsc(coset(0, ay, az), 1, 1)- &
                                                      f2*kvec(1)*ss(coset(0, ay, az), 1)- &
                                                      f2*kvec(1)*dss(coset(0, ay, az), 1, 1)

                        dsc(coset(1, ay, az), 1, 2) = rap(1)*dsc(coset(0, ay, az), 1, 2)- &
                                                      f2*kvec(1)*dss(coset(0, ay, az), 1, 2)

                        dsc(coset(1, ay, az), 1, 3) = rap(1)*dsc(coset(0, ay, az), 1, 3)- &
                                                      f2*kvec(1)*dss(coset(0, ay, az), 1, 3)
                     END IF

                     IF (PRESENT(dkss)) THEN
                        dss(coset(1, ay, az), 1, 1) = rap(1)*dss(coset(0, ay, az), 1, 1)+ &
                                                      f2*kvec(1)*sc(coset(0, ay, az), 1)+ &
                                                      f2*kvec(1)*dsc(coset(0, ay, az), 1, 1)

                        dss(coset(1, ay, az), 1, 2) = rap(1)*dss(coset(0, ay, az), 1, 2)+ &
                                                      f2*kvec(1)*dsc(coset(0, ay, az), 1, 2)

                        dss(coset(1, ay, az), 1, 3) = rap(1)*dss(coset(0, ay, az), 1, 3)+ &
                                                      f2*kvec(1)*dsc(coset(0, ay, az), 1, 3)
                     END IF
                  END DO

                  DO ax = 2, la
                     f3 = f2*REAL(ax-1, dp)
                     DO ay = 0, la-ax
                        az = la-ax-ay
                        sc(coset(ax, ay, az), 1) = rap(1)*sc(coset(ax-1, ay, az), 1)+ &
                                                   f3*sc(coset(ax-2, ay, az), 1)- &
                                                   f2*kvec(1)*ss(coset(ax-1, ay, az), 1)
                        ss(coset(ax, ay, az), 1) = rap(1)*ss(coset(ax-1, ay, az), 1)+ &
                                                   f3*ss(coset(ax-2, ay, az), 1)+ &
                                                   f2*kvec(1)*sc(coset(ax-1, ay, az), 1)

                        IF (PRESENT(dksc)) THEN
                           dsc(coset(ax, ay, az), 1, 1) = rap(1)*dsc(coset(ax-1, ay, az), 1, 1)+ &
                                                          f3*dsc(coset(ax-2, ay, az), 1, 1)- &
                                                          f2*kvec(1)*ss(coset(ax-1, ay, az), 1)- &
                                                          f2*kvec(1)*dss(coset(ax-1, ay, az), 1, 1)
                           dsc(coset(ax, ay, az), 1, 2) = rap(1)*dsc(coset(ax-1, ay, az), 1, 2)+ &
                                                          f3*dsc(coset(ax-2, ay, az), 1, 2)- &
                                                          f2*kvec(1)*dss(coset(ax-1, ay, az), 1, 2)
                           dsc(coset(ax, ay, az), 1, 3) = rap(1)*dsc(coset(ax-1, ay, az), 1, 3)+ &
                                                          f3*dsc(coset(ax-2, ay, az), 1, 3)- &
                                                          f2*kvec(1)*dss(coset(ax-1, ay, az), 1, 3)
                        END IF
                        IF (PRESENT(dkss)) THEN
                           dss(coset(ax, ay, az), 1, 1) = rap(1)*dss(coset(ax-1, ay, az), 1, 1)+ &
                                                          f3*dss(coset(ax-2, ay, az), 1, 1)+ &
                                                          f2*kvec(1)*sc(coset(ax-1, ay, az), 1)+ &
                                                          f2*kvec(1)*dsc(coset(ax-1, ay, az), 1, 1)
                           dss(coset(ax, ay, az), 1, 2) = rap(1)*dss(coset(ax-1, ay, az), 1, 2)+ &
                                                          f3*dss(coset(ax-2, ay, az), 1, 2)+ &
                                                          f2*kvec(1)*dsc(coset(ax-1, ay, az), 1, 2)
                           dss(coset(ax, ay, az), 1, 3) = rap(1)*dss(coset(ax-1, ay, az), 1, 3)+ &
                                                          f3*dss(coset(ax-2, ay, az), 1, 3)+ &
                                                          f2*kvec(1)*dsc(coset(ax-1, ay, az), 1, 3)
                        END IF
                     END DO
                  END DO

               END DO

!         *** Recurrence steps: [a|O|s] -> [a|O|b] ***

               IF (lb_max > 0) THEN

                  DO j = 2, ncoset(lb_max)
                     DO i = 1, ncoset(la_max)
                        sc(i, j) = 0.0_dp
                        ss(i, j) = 0.0_dp

                        IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                           DO k = 1, 3
                              dsc(i, j, k) = 0.0_dp
                              dss(i, j, k) = 0.0_dp
                           END DO
                        END IF
                        IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                           DO k = 1, 6
                              ddsc(i, j, k) = 0.0_dp
                              ddss(i, j, k) = 0.0_dp
                           END DO
                        END IF

                     END DO
                  END DO

!           *** Horizontal recurrence steps ***

                  rbp(:) = rap(:)-rab(:)

!           *** [a|O|p] = [a+1i|O|s] - (Bi - Ai)*[a|O|s] ***

                  IF (lb_max == 1) THEN
                     la_start = la_min
                  ELSE
                     la_start = MAX(0, la_min-1)
                  END IF

                  DO la = la_start, la_max-1
                     DO ax = 0, la
                        DO ay = 0, la-ax
                           az = la-ax-ay
                           sc(coset(ax, ay, az), 2) = sc(coset(ax+1, ay, az), 1)- &
                                                      rab(1)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 3) = sc(coset(ax, ay+1, az), 1)- &
                                                      rab(2)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 4) = sc(coset(ax, ay, az+1), 1)- &
                                                      rab(3)*sc(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = ss(coset(ax+1, ay, az), 1)- &
                                                      rab(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = ss(coset(ax, ay+1, az), 1)- &
                                                      rab(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = ss(coset(ax, ay, az+1), 1)- &
                                                      rab(3)*ss(coset(ax, ay, az), 1)

                        IF (PRESENT(dksc)) THEN
                           DO i = 1, 3
                              dsc(coset(ax, ay, az), 2, i) = dsc(coset(ax+1, ay, az), 1, i)- &
                                                             rab(1)*dsc(coset(ax, ay, az), 1, i)
                              dsc(coset(ax, ay, az), 3, i) = dsc(coset(ax, ay+1, az), 1, i)- &
                                                             rab(2)*dsc(coset(ax, ay, az), 1, i)
                              dsc(coset(ax, ay, az), 4, i) = dsc(coset(ax, ay, az+1), 1, i)- &
                                                             rab(3)*dsc(coset(ax, ay, az), 1, i)
                           END DO
                        END IF
                        
                        IF (PRESENT(dkss)) THEN
                           DO i = 1, 3
                              dss(coset(ax, ay, az), 2, i) = dss(coset(ax+1, ay, az), 1, i)- &
                                                             rab(1)*dss(coset(ax, ay, az), 1, i)
                              dss(coset(ax, ay, az), 3, i) = dss(coset(ax, ay+1, az), 1, i)- &
                                                             rab(2)*dss(coset(ax, ay, az), 1, i)
                              dss(coset(ax, ay, az), 4, i) = dss(coset(ax, ay, az+1), 1, i)- &
                                                             rab(3)*dss(coset(ax, ay, az), 1, i)
                           END DO
                        END IF
                        END DO
                     END DO
                  END DO

!           *** Vertical recurrence step ***

!           *** [a|O|p] = (Pi - Bi)*[a|O|s] + f2*Ni(a)*[a-1i|O|s] ***
!           ***           + [a|dO|s]                              ***

                  DO ax = 0, la_max
                     fx = f2*REAL(ax, dp)
                     DO ay = 0, la_max-ax
                        fy = f2*REAL(ay, dp)
                        az = la_max-ax-ay
                        fz = f2*REAL(az, dp)
                        IF (ax == 0) THEN
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------

                           IF (PRESENT(dksc)) THEN
                              dsc(coset(ax, ay, az), 2, 1) = rbp(1)*dsc(coset(ax, ay, az), 1, 1)- &
                                                             f2*kvec(1)*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 2, 2) = rbp(1)*dsc(coset(ax, ay, az), 1, 2)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 2, 3) = rbp(1)*dsc(coset(ax, ay, az),1 , 3)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(dkss)) THEN
                              dss(coset(ax, ay, az), 2, 1) = rbp(1)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             f2*kvec(1)*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 2, 2) = rbp(1)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 2, 3) = rbp(1)*dss(coset(ax, ay, az),1 , 3)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                        ELSE
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1)+ &
                                                      fx*sc(coset(ax-1, ay, az), 1)- &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1)+ &
                                                      fx*ss(coset(ax-1, ay, az), 1)+ &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------

                           IF (PRESENT(dksc)) THEN
                              dsc(coset(ax, ay, az), 2, 1) = rbp(1)*dsc(coset(ax, ay, az), 1, 1)+ &
                                                             fx*dsc(coset(ax-1, ay, az), 1, 1)- &
                                                             f2*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 2, 2) = rbp(1)*dsc(coset(ax, ay, az), 1, 2)+ &
                                                             fx*kvec(1)*dsc(coset(ax-1, ay, az), 1, 2)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 2, 3) = rbp(1)*dsc(coset(ax, ay, az), 1, 3)+ &
                                                             fx*dsc(coset(ax-1, ay, az), 1, 3)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(dkss)) THEN
                              dss(coset(ax, ay, az), 2, 1) = rbp(1)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             fx*dss(coset(ax-1, ay, az), 1, 1)+ &
                                                             f2*kvec(1)*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 2, 2) = rbp(1)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             fx*dss(coset(ax-1, ay, az), 1, 2)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 2, 3) = rbp(1)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             fx*dss(coset(ax-1, ay, az), 1, 3)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                        END IF
                        IF (ay == 0) THEN
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------
                        
                           IF (PRESENT(dksc)) THEN
                              dsc(coset(ax, ay, az), 3, 1) = rbp(2)*dsc(coset(ax, ay, az), 1, 1)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 3, 2) = rbp(2)*dsc(coset(ax, ay, az), 1, 2)- &
                                                             f2*kvec(2)*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 3, 3) = rbp(2)*dsc(coset(ax, ay, az), 1, 3)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(dkss)) THEN
                              dss(coset(ax, ay, az), 3, 1) = rbp(2)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 3, 2) = rbp(2)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             f2*kvec(2)*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 3, 3) = rbp(2)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                        ELSE
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1)+ &
                                                      fy*sc(coset(ax, ay-1, az), 1)- &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1)+ &
                                                      fy*ss(coset(ax, ay-1, az), 1)+ &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------
                           
                           IF (PRESENT(dksc)) THEN
                              dsc(coset(ax, ay, az), 3, 1) = rbp(2)*dsc(coset(ax, ay, az), 1, 1)+ &
                                                             fy*dsc(coset(ax, ay-1, az), 1, 1)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 3, 2) = rbp(2)*dsc(coset(ax, ay, az), 1, 2)+ &
                                                             fy*dsc(coset(ax, ay-1, az), 1, 2)- &
                                                             f2*kvec(2)*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 3, 3) = rbp(2)*dsc(coset(ax, ay, az), 1, 3)+ &
                                                             fy*dsc(coset(ax, ay-1, az), 1, 3)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(dkss)) THEN
                              dss(coset(ax, ay, az), 3, 1) = rbp(2)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             fy*dss(coset(ax, ay-1, az), 1, 1)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 3, 2) = rbp(2)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             fy*dss(coset(ax, ay-1, az), 1, 2)+ &
                                                             f2*kvec(2)*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 3, 3) = rbp(2)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             fy*dss(coset(ax, ay-1, az), 1, 3)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                        END IF
                        IF (az == 0) THEN
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------

                           IF (PRESENT(dksc)) THEN                        
                              dsc(coset(ax, ay, az), 4, 1) = rbp(3)*dsc(coset(ax, ay, az), 1, 1)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 4, 2) = rbp(3)*dsc(coset(ax, ay, az), 1, 2)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 4, 3) = rbp(3)*dsc(coset(ax, ay, az), 1, 3)- &
                                                             f2*kvec(3)*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(dkss)) THEN                        
                              dss(coset(ax, ay, az), 4, 1) = rbp(3)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 4, 2) = rbp(3)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 4, 3) = rbp(3)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             f2*kvec(3)*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                        ELSE
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1)+ &
                                                      fz*sc(coset(ax, ay, az-1), 1)- &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1)+ &
                                                      fz*ss(coset(ax, ay, az-1), 1)+ &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------
                  
                                          
                           IF (PRESENT(dksc)) THEN                        
                              dsc(coset(ax, ay, az), 4, 1) = rbp(3)*dsc(coset(ax, ay, az), 1, 1)+ &
                                                             fz*dsc(coset(ax, ay, az-1), 1, 1)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 4, 2) = rbp(3)*dsc(coset(ax, ay, az), 1, 2)+ &
                                                             fz*dsc(coset(ax, ay, az-1), 1, 2)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 4, 3) = rbp(3)*dsc(coset(ax, ay, az), 1, 3)+ &
                                                             fz*dsc(coset(ax, ay, az-1), 1, 3)- &
                                                             f2*kvec(3)*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(dkss)) THEN                        
                              dss(coset(ax, ay, az), 4, 1) = rbp(3)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             fz*dss(coset(ax, ay, az-1), 1, 1)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 4, 2) = rbp(3)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             fz*dss(coset(ax, ay, az-1), 1, 2)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 4, 3) = rbp(3)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             fz*dss(coset(ax, ay, az-1), 1, 3)+ &
                                                             f2*kvec(3)*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                        END IF
                     END DO
                  END DO

!           *** Recurrence steps: [a|O|p] -> [a|O|b] ***

                  DO lb = 2, lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|O|b] = [a+1i|O|b-1i] - (Bi - Ai)*[a|O|b-1i] ***

                     IF (lb == lb_max) THEN
                        la_start = la_min
                     ELSE
                        la_start = MAX(0, la_min-1)
                     END IF

                     DO la = la_start, la_max-1
                        DO ax = 0, la
                           DO ay = 0, la-ax
                              az = la-ax-ay

!                   *** Shift of angular momentum component z from a to b ***

                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 sc(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                 rab(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 ss(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                 rab(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))

                              IF (PRESENT(dksc)) THEN
                                 DO i = 1, 3
                                    dsc(coset(ax, ay, az), coset(0, 0, lb), i) = &
                                       dsc(coset(ax, ay, az+1), coset(0, 0, lb-1), i)- &
                                       rab(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), i)
                                 END DO
                              END IF
                              IF (PRESENT(dkss)) THEN
                                 DO i = 1, 3
                                    dss(coset(ax, ay, az), coset(0, 0, lb), i) = &
                                       dss(coset(ax, ay, az+1), coset(0, 0, lb-1), i)- &
                                       rab(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), i)
                                 END DO
                              END IF
!                   *** Shift of angular momentum component y from a to b ***

                              DO by = 1, lb
                                 bz = lb-by
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    sc(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                    rab(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    ss(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                    rab(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))

                                 IF (PRESENT(dksc)) THEN
                                    DO i = 1, 3
                                       dsc(coset(ax, ay, az), coset(0, by, bz), i) = &
                                          dsc(coset(ax, ay+1, az), coset(0, by-1, bz), i)- &
                                          rab(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), i)
                                    END DO
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    DO i = 1, 3
                                       dss(coset(ax, ay, az), coset(0, by, bz), i) = &
                                          dss(coset(ax, ay+1, az), coset(0, by-1, bz), i)- &
                                          rab(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), i)
                                    END DO
                                 END IF

                              END DO

!                   *** Shift of angular momentum component x from a to b ***

                              DO bx = 1, lb
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       sc(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                       rab(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       ss(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                       rab(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))

                                    IF (PRESENT(dksc)) THEN
                                       DO i = 1, 3
                                          dsc(coset(ax, ay, az), coset(bx, by, bz), i) = &
                                             dsc(coset(ax+1, ay, az), coset(bx-1, by, bz), i)- &
                                             rab(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), i)
                                       END DO
                                    END IF
                                    IF (PRESENT(dkss)) THEN
                                       DO i = 1, 3
                                          dss(coset(ax, ay, az), coset(bx, by, bz), i) = &
                                             dss(coset(ax+1, ay, az), coset(bx-1, by, bz), i)- &
                                             rab(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), i)
                                       END DO
                                    END IF
                                 END DO
                              END DO

                           END DO
                        END DO
                     END DO

!             *** Vertical recurrence step ***

!             *** [a|O|b] = (Pi - Bi)*[a|O|b-1i] + f2*Ni(a)*[a-1i|O|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|O|b-2i] + [a|dO|b-1i]            ***

                     DO ax = 0, la_max
                        fx = f2*REAL(ax, dp)
                        DO ay = 0, la_max-ax
                           fy = f2*REAL(ay, dp)
                           az = la_max-ax-ay
                           fz = f2*REAL(az, dp)

!                 *** Increase the angular momentum component z of function b ***

                           f3 = f2*REAL(lb-1, dp)

                           IF (az == 0) THEN
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              !------------------------------------------------
        
                              IF (PRESENT(dksc)) THEN
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                    rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                    f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 1)- &
                                    f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                    rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                    f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 2)- &
                                    f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                    rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                    f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 3)- &
                                    f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                    f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                              END IF
                              IF (PRESENT(dkss)) THEN
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 1)+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 2)+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 3)+ &
                                     f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                              END IF
                           ELSE
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 fz*sc(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 fz*ss(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              !------------------------------------------------

                              IF (PRESENT(dksc)) THEN
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                     rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                     fz*dsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 1)+ &
                                     f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 1)- &
                                     f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                     rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                     fz*dsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 2)+ &
                                     f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 2)- &
                                     f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                     rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                     fz*dsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 3)+ &
                                     f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 3)- &
                                     f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                     f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                              END IF
                              IF (PRESENT(dkss)) THEN
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                     fz*dss(coset(ax, ay, az-1), coset(0, 0, lb-1), 1)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 1)+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                     fz*dss(coset(ax, ay, az-1), coset(0, 0, lb-1), 2)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 2)+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                     fz*dss(coset(ax, ay, az-1), coset(0, 0, lb-1), 3)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 3)+ &
                                     f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                              END IF
                           END IF

!                 *** Increase the angular momentum component y of function b ***

                           IF (ay == 0) THEN
                              bz = lb-1
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz))- &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              !------------------------------------------------

                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                        f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 3)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                        f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)
                                 END IF
                              DO by = 2, lb
                                 bz = lb-by
                                 f3 = f2*REAL(by-1, dp)
                                 !------------------------------------------------
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    f3*sc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    f3*ss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 !------------------------------------------------
                              
                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 1)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 2)- &
                                        f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 3)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 1)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    dss(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 2)+ &
                                        f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    dss(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 3)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                 END IF
                              END DO
                           ELSE
                              bz = lb-1
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 fy*sc(coset(ax, ay-1, az), coset(0, 0, bz))- &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 fy*ss(coset(ax, ay-1, az), coset(0, 0, bz))+ &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              !------------------------------------------------

                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, 0, bz), 1)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, 0, bz), 2)- &
                                        f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, 0, bz), 3)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 3)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, 0, bz), 1)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, 0, bz), 2)+ &
                                        f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, 0, bz), 3)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)
                                 END IF

                              DO by = 2, lb
                                 bz = lb-by
                                 f3 = f2*REAL(by-1, dp)
                                 !------------------------------------------------
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    fy*sc(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                    f3*sc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    fy*ss(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                    f3*ss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 !------------------------------------------------

                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, by-1, bz), 1)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 1)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, by-1, bz), 2)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 2)- &
                                        f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, by-1, bz), 3)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 3)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, by-1, bz), 1)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 1)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    dss(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, by-1, bz), 2)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 2)+ &
                                        f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    dss(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, by-1, bz), 3)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 3)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                 END IF
                              END DO
                           END IF

!                 *** Increase the angular momentum component x of function b ***

                           IF (ax == 0) THEN
                              DO by = 0, lb-1
                                 bz = lb-1-by
                                 !----------------------------------------
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz))- &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz))+ &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                                 !----------------------------------------
                              
                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 1)- &
                                        f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 1)
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 2)- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 2)
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 3)- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 3)
                                 END IF

                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 1)+ &
                                        f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 1)
                                    dss(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 2)+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 2)
                                    dss(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 3)+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 3)
                                 END IF
                              END DO
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx-1, dp)
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    !----------------------------------------
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       f3*sc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       f3*ss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    !----------------------------------------

                                    IF (PRESENT(dksc)) THEN
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 1)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 1)- &
                                           f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 2)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 2)- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 3)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 3)- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                    END IF
                                    IF (PRESENT(dkss)) THEN
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 1)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 1)+ &
                                           f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 2)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 2)+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 3)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 3)+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                    END IF
                                 END DO
                              END DO
                           ELSE
                              DO by = 0, lb-1
                                 bz = lb-1-by
                                 !----------------------------------------
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                    fx*sc(coset(ax-1, ay, az), coset(0, by, bz))- &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz))+ &
                                    fx*ss(coset(ax-1, ay, az), coset(0, by, bz))+ &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                                 !----------------------------------------
                              
                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 1)+ &
                                        fx*dsc(coset(ax-1, ay, az), coset(0, by, bz), 1)- &
                                        f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 1)
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 2)+ &
                                        fx*dsc(coset(ax-1, ay, az), coset(0, by, bz), 2)- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 2)
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 3)+ &
                                        fx*dsc(coset(ax-1, ay, az), coset(0, by, bz), 3)- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 3)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 1)+ &
                                        fx*dss(coset(ax-1, ay, az), coset(0, by, bz), 1)+ &
                                        f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 1)
                                    dss(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 2)+ &
                                        fx*dss(coset(ax-1, ay, az), coset(0, by, bz), 2)+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 2)
                                    dss(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 3)+ &
                                        fx*dss(coset(ax-1, ay, az), coset(0, by, bz), 3)+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 3)
                                 END IF
                              END DO
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx-1, dp)
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    !----------------------------------------
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       fx*sc(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                       f3*sc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       fx*ss(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                       f3*ss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    !----------------------------------------
                                 
                                    IF (PRESENT(dksc)) THEN
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 1)+ &
                                           fx*dsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 1)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 1)- &
                                           f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 2)+ &
                                           fx*dsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 2)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 2)- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 3)+ &
                                           fx*dsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 3)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 3)- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                    END IF
                                    IF (PRESENT(dkss)) THEN
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 1)+ &
                                           fx*dss(coset(ax-1, ay, az), coset(bx-1, by, bz), 1)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 1)+ &
                                           f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 2)+ &
                                           fx*dss(coset(ax-1, ay, az), coset(bx-1, by, bz), 2)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 2)+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 3)+ &
                                           fx*dss(coset(ax-1, ay, az), coset(bx-1, by, bz), 3)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 3)+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                    END IF
                                 
                                 END DO
                              END DO
                           END IF

                        END DO
                     END DO

                  END DO

               END IF

            ELSE

               IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|O|s] -> [s|O|b] ***

                  rbp(:) = (f1-1.0_dp)*rab(:)

!           *** [s|O|p] = (Pi - Bi)*[s|O|s] + [s|dO|s] ***

                  !----------------------------------------
                  sc(1, 2) = rbp(1)*sc(1, 1)-f2*kvec(1)*ss(1, 1)
                  sc(1, 3) = rbp(2)*sc(1, 1)-f2*kvec(2)*ss(1, 1)
                  sc(1, 4) = rbp(3)*sc(1, 1)-f2*kvec(3)*ss(1, 1)
                  ss(1, 2) = rbp(1)*ss(1, 1)+f2*kvec(1)*sc(1, 1)
                  ss(1, 3) = rbp(2)*ss(1, 1)+f2*kvec(2)*sc(1, 1)
                  ss(1, 4) = rbp(3)*ss(1, 1)+f2*kvec(3)*sc(1, 1)
                  !----------------------------------------

               IF (PRESENT(dksc)) THEN
                  dsc(1,2,1) = rbp(1)*dsc(1,1,1) - f2*kvec(1)*ss(1,1) - f2*kvec(1)*dss(1,1,1) 
                  dsc(1,2,2) = rbp(1)*dsc(1,1,2) - f2*kvec(1)*dss(1,1,2)
                  dsc(1,2,3) = rbp(1)*dsc(1,1,3) - f2*kvec(1)*dss(1,1,3)
                  
                  dsc(1,3,1) = rbp(2)*dsc(1,1,1) - f2*kvec(2)*dss(1,1,1) 
                  dsc(1,3,2) = rbp(2)*dsc(1,1,2) - f2*kvec(2)*ss(1,1) - f2*kvec(2)*dss(1,1,2)
                  dsc(1,3,3) = rbp(2)*dsc(1,1,3) - f2*kvec(2)*dss(1,1,3)
                  
                  dsc(1,4,1) = rbp(3)*dsc(1,1,1) - f2*kvec(3)*dss(1,1,1) 
                  dsc(1,4,2) = rbp(3)*dsc(1,1,2) - f2*kvec(3)*dss(1,1,2)
                  dsc(1,4,3) = rbp(3)*dsc(1,1,3) - f2*kvec(3)*ss(1,1) - f2*kvec(3)*dss(1,1,3)
               END IF

               IF (PRESENT(dkss)) THEN
                  dss(1,2,1) = rbp(1)*dss(1,1,1) + f2*kvec(1)*sc(1,1) + f2*kvec(1)*dsc(1,1,1) 
                  dss(1,2,2) = rbp(1)*dss(1,1,2) + f2*kvec(1)*dsc(1,1,2)
                  dss(1,2,3) = rbp(1)*dss(1,1,3) + f2*kvec(1)*dsc(1,1,3)

                  dss(1,3,1) = rbp(2)*dss(1,1,1) + f2*kvec(2)*dsc(1,1,1) 
                  dss(1,3,2) = rbp(2)*dss(1,1,2) + f2*kvec(2)*sc(1,1) + f2*kvec(2)*dsc(1,1,2)
                  dss(1,3,3) = rbp(2)*dss(1,1,3) + f2*kvec(2)*dsc(1,1,3)

                  dss(1,4,1) = rbp(3)*dss(1,1,1) + f2*kvec(3)*dsc(1,1,1) 
                  dss(1,4,2) = rbp(3)*dss(1,1,2) + f2*kvec(3)*dsc(1,1,2)
                  dss(1,4,3) = rbp(3)*dss(1,1,3) + f2*kvec(3)*sc(1,1) + f2*kvec(3)*dsc(1,1,3)
               END IF

!           *** [s|O|b] = (Pi - Bi)*[s|O|b-1i] + f2*Ni(b-1i)*[s|O|b-2i] ***
!           ***           + [s|dO|b-1i]                                 ***

                  DO lb = 2, lb_max

!             *** Increase the angular momentum component z of function b ***

                     !-------------------------------------
                     sc(1, coset(0, 0, lb)) = rbp(3)*sc(1, coset(0, 0, lb-1))+ &
                                              f2*REAL(lb-1, dp)*sc(1, coset(0, 0, lb-2))- &
                                              f2*kvec(3)*ss(1, coset(0, 0, lb-1))
                     ss(1, coset(0, 0, lb)) = rbp(3)*ss(1, coset(0, 0, lb-1))+ &
                                              f2*REAL(lb-1, dp)*ss(1, coset(0, 0, lb-2))+ &
                                              f2*kvec(3)*sc(1, coset(0, 0, lb-1))
                     !-------------------------------------
              
                     IF (PRESENT(dksc)) THEN
                        dsc(1, coset(0, 0, lb), 1) = &
                            rbp(3)*dsc(1, coset(0, 0, lb-1), 1)+ &
                            f2*REAL(lb-1, dp)*dsc(1, coset(0, 0, lb-2), 1)- &
                            f2*kvec(3)*dss(1, coset(0, 0, lb-1), 1)
                        dsc(1, coset(0, 0, lb), 2) = &
                            rbp(3)*dsc(1, coset(0, 0, lb-1), 2)+ &
                            f2*REAL(lb-1, dp)*dsc(1, coset(0, 0, lb-2), 2)- &
                            f2*kvec(3)*dss(1, coset(0, 0, lb-1), 2)
                        dsc(1, coset(0, 0, lb), 3) = &
                            rbp(3)*dsc(1, coset(0, 0, lb-1), 3)+ &
                            f2*REAL(lb-1, dp)*dsc(1, coset(0, 0, lb-2), 3)- &
                            f2*kvec(3)*ss(1, coset(0, 0, lb-1))- &
                            f2*kvec(3)*dss(1, coset(0, 0, lb-1), 3)
                     END IF
                     IF (PRESENT(dkss)) THEN
                        dss(1, coset(0, 0, lb), 1) = &
                            rbp(3)*dss(1, coset(0, 0, lb-1), 1)+ &
                            f2*REAL(lb-1, dp)*dss(1, coset(0, 0, lb-2), 1)+ &
                            f2*kvec(3)*dsc(1, coset(0, 0, lb-1), 1)
                        dss(1, coset(0, 0, lb), 2) = &
                            rbp(3)*dss(1, coset(0, 0, lb-1), 2)+ &
                            f2*REAL(lb-1, dp)*dss(1, coset(0, 0, lb-2), 2)+ &
                            f2*kvec(3)*dsc(1, coset(0, 0, lb-1), 2)
                        dss(1, coset(0, 0, lb), 3) = &
                            rbp(3)*dss(1, coset(0, 0, lb-1), 3)+ &
                            f2*REAL(lb-1, dp)*dss(1, coset(0, 0, lb-2), 3)+ &
                            f2*kvec(3)*sc(1, coset(0, 0, lb-1))- &
                            f2*kvec(3)*dsc(1, coset(0, 0, lb-1), 3)
                     END IF

!             *** Increase the angular momentum component y of function b ***

                     bz = lb-1
                     !-------------------------------------
                     sc(1, coset(0, 1, bz)) = rbp(2)*sc(1, coset(0, 0, bz))- &
                                              f2*kvec(2)*ss(1, coset(0, 0, bz))
                     ss(1, coset(0, 1, bz)) = rbp(2)*ss(1, coset(0, 0, bz))+ &
                                              f2*kvec(2)*sc(1, coset(0, 0, bz))
                     !-------------------------------------

                     IF (PRESENT(dksc)) THEN
                        dsc(1, coset(0, 1, bz), 1) = rbp(2)*dsc(1, coset(0, 0, bz), 1)- &
                                                     f2*kvec(2)*dss(1, coset(0, 0, bz), 1)
                        dsc(1, coset(0, 1, bz), 2) = rbp(2)*dsc(1, coset(0, 0, bz), 2)- &
                                                     f2*kvec(2)*ss(1, coset(0, 0, bz))- &
                                                     f2*kvec(2)*dss(1, coset(0, 0, bz), 2)
                        dsc(1, coset(0, 1, bz), 3) = rbp(2)*dsc(1, coset(0, 0, bz), 3)- &
                                                     f2*kvec(2)*dss(1, coset(0, 0, bz), 3)
                     END IF

                     IF (PRESENT(dkss)) THEN
                        dss(1, coset(0, 1, bz), 1) = rbp(2)*dss(1, coset(0, 0, bz), 1)+ &
                                                     f2*kvec(2)*dsc(1, coset(0, 0, bz), 1)
                        dss(1, coset(0, 1, bz), 2) = rbp(2)*dss(1, coset(0, 0, bz), 2)+ &
                                                     f2*kvec(2)*sc(1, coset(0, 0, bz))+ &
                                                     f2*kvec(2)*dsc(1, coset(0, 0, bz), 2)
                        dss(1, coset(0, 1, bz), 3) = rbp(2)*dss(1, coset(0, 0, bz), 3)+ &
                                                     f2*kvec(2)*dsc(1, coset(0, 0, bz), 3)
                     END IF

                     DO by = 2, lb
                        bz = lb-by
                        !-------------------------------------
                        sc(1, coset(0, by, bz)) = rbp(2)*sc(1, coset(0, by-1, bz))+ &
                                                  f2*REAL(by-1, dp)*sc(1, coset(0, by-2, bz))- &
                                                  f2*kvec(2)*ss(1, coset(0, by-1, bz))
                        ss(1, coset(0, by, bz)) = rbp(2)*ss(1, coset(0, by-1, bz))+ &
                                                  f2*REAL(by-1, dp)*ss(1, coset(0, by-2, bz))+ &
                                                  f2*kvec(2)*sc(1, coset(0, by-1, bz))
                        !-------------------------------------
                    
                        IF (PRESENT(dksc)) THEN
                           dsc(1, coset(0, by, bz), 1) = &
                               rbp(2)*dsc(1, coset(0, by-1, bz), 1)+ &
                               f2*REAL(by-1, dp)*dsc(1, coset(0, by-2, bz), 1)- &
                               f2*kvec(2)*dss(1, coset(0, by-1, bz), 1)
                           dsc(1, coset(0, by, bz), 2) = &
                               rbp(2)*dsc(1, coset(0, by-1, bz), 2)+ &
                               f2*REAL(by-1, dp)*dsc(1, coset(0, by-2, bz), 2)- &
                               f2*kvec(2)*ss(1, coset(0, by-1, bz))- &
                               f2*kvec(2)*dss(1, coset(0, by-1, bz), 2)
                           dsc(1, coset(0, by, bz), 3) = &
                               rbp(2)*dsc(1, coset(0, by-1, bz), 3)+ &
                               f2*REAL(by-1, dp)*dsc(1, coset(0, by-2, bz), 3)- &
                               f2*kvec(2)*dss(1, coset(0, by-1, bz), 3)
                        END IF
                        IF (PRESENT(dkss)) THEN
                           dss(1, coset(0, by, bz), 1) = &
                               rbp(2)*dss(1, coset(0, by-1, bz), 1)+ &
                               f2*REAL(by-1, dp)*dss(1, coset(0, by-2, bz), 1)+ &
                               f2*kvec(2)*dsc(1, coset(0, by-1, bz), 1)
                           dss(1, coset(0, by, bz), 2) = &
                               rbp(2)*dss(1, coset(0, by-1, bz), 2)+ &
                               f2*REAL(by-1, dp)*dss(1, coset(0, by-2, bz), 2)+ &
                               f2*kvec(2)*sc(1, coset(0, by-1, bz))+ &
                               f2*kvec(2)*dsc(1, coset(0, by-1, bz), 2)
                           dss(1, coset(0, by, bz), 3) = &
                               rbp(2)*dss(1, coset(0, by-1, bz), 3)+ &
                               f2*REAL(by-1, dp)*dss(1, coset(0, by-2, bz), 3)+ &
                               f2*kvec(2)*dsc(1, coset(0, by-1, bz), 3)
                        END IF
                    
                     END DO

!             *** Increase the angular momentum component x of function b ***

                     DO by = 0, lb-1
                        bz = lb-1-by
                        !----------------------------------
                        sc(1, coset(1, by, bz)) = rbp(1)*sc(1, coset(0, by, bz))- &
                                                  f2*kvec(1)*ss(1, coset(0, by, bz))
                        ss(1, coset(1, by, bz)) = rbp(1)*ss(1, coset(0, by, bz))+ &
                                                  f2*kvec(1)*sc(1, coset(0, by, bz))
                        !----------------------------------
                                               
                        IF (PRESENT(dksc)) THEN
                           dsc(1, coset(1, by, bz), 1) = &
                               rbp(1)*dsc(1, coset(0, by, bz), 1)- &
                               f2*kvec(1)*ss(1, coset(0, by, bz))- &
                               f2*kvec(1)*dss(1, coset(0, by, bz), 1)
                           dsc(1, coset(1, by, bz), 2) = &
                               rbp(1)*dsc(1, coset(0, by, bz), 2)- &
                               f2*kvec(1)*dss(1, coset(0, by, bz), 2)
                           dsc(1, coset(1, by, bz), 3) = &
                               rbp(1)*dsc(1, coset(0, by, bz), 3)- &
                               f2*kvec(1)*dss(1, coset(0, by, bz), 3)
                        END IF                       
                        IF (PRESENT(dkss)) THEN
                           dss(1, coset(1, by, bz), 1) = &
                               rbp(1)*dss(1, coset(0, by, bz), 1)+ &
                               f2*kvec(1)*sc(1, coset(0, by, bz))+ &
                               f2*kvec(1)*dsc(1, coset(0, by, bz), 1)
                           dss(1, coset(1, by, bz), 2) = &
                               rbp(1)*dss(1, coset(0, by, bz), 2)+ &
                               f2*kvec(1)*dsc(1, coset(0, by, bz), 2)
                           dss(1, coset(1, by, bz), 3) = &
                               rbp(1)*dss(1, coset(0, by, bz), 3)+ &
                               f2*kvec(1)*dsc(1, coset(0, by, bz), 3)
                        END IF                       
                     END DO

                     DO bx = 2, lb
                        f3 = f2*REAL(bx-1, dp)
                        DO by = 0, lb-bx
                           bz = lb-bx-by
                           !----------------------------------
                           sc(1, coset(bx, by, bz)) = rbp(1)*sc(1, coset(bx-1, by, bz))+ &
                                                      f3*sc(1, coset(bx-2, by, bz))- &
                                                      f2*kvec(1)*ss(1, coset(bx-1, by, bz))
                           ss(1, coset(bx, by, bz)) = rbp(1)*ss(1, coset(bx-1, by, bz))+ &
                                                      f3*ss(1, coset(bx-2, by, bz))+ &
                                                      f2*kvec(1)*sc(1, coset(bx-1, by, bz))
                          !----------------------------------
                                                   
                           IF (PRESENT(dksc)) THEN
                              dsc(1, coset(bx, by, bz), 1) = &
                                  rbp(1)*dsc(1, coset(bx-1, by, bz), 1)+ &
                                  f3*dsc(1, coset(bx-2, by, bz), 1)- &
                                  f2*kvec(1)*ss(1, coset(bx-1, by, bz))- &
                                  f2*kvec(1)*dss(1, coset(bx-1, by, bz), 1)
                              dsc(1, coset(bx, by, bz), 2) = &
                                  rbp(1)*dsc(1, coset(bx-1, by, bz), 2)+ &
                                  f3*dsc(1, coset(bx-2, by, bz), 2)- &
                                  f2*kvec(1)*dss(1, coset(bx-1, by, bz), 2)
                              dsc(1, coset(bx, by, bz), 3) = &
                                  rbp(1)*dsc(1, coset(bx-1, by, bz), 3)+ &
                                  f3*dsc(1, coset(bx-2, by, bz), 3)- &
                                  f2*kvec(1)*dss(1, coset(bx-1, by, bz), 3)
                           END IF                       
                           IF (PRESENT(dkss)) THEN
                              dss(1, coset(bx, by, bz), 1) = &
                                  rbp(1)*dss(1, coset(bx-1, by, bz), 1)+ &
                                  f3*dss(1, coset(bx-2, by, bz), 1)+ &
                                  f2*kvec(1)*sc(1, coset(bx-1, by, bz))+ &
                                  f2*kvec(1)*dsc(1, coset(bx-1, by, bz), 1)
                              dss(1, coset(bx, by, bz), 2) = &
                                  rbp(1)*dss(1, coset(bx-1, by, bz), 2)+ &
                                  f3*dss(1, coset(bx-2, by, bz), 2)+ &
                                  f2*kvec(1)*dsc(1, coset(bx-1, by, bz), 2)
                              dss(1, coset(bx, by, bz), 3) = &
                                  rbp(1)*dss(1, coset(bx-1, by, bz), 3)+ &
                                  f3*dss(1, coset(bx-2, by, bz), 3)+ &
                                  f2*kvec(1)*dsc(1, coset(bx-1, by, bz), 3)
                           END IF                       
                        END DO
                     END DO

                  END DO

               END IF

            END IF

            DO j = ncoset(lb_min-1)+1, ncoset(lb_max)
               DO i = ncoset(la_min_set-1)+1, ncoset(la_max_set)
                  cosab(na+i, nb+j) = sc(i, j)
                  sinab(na+i, nb+j) = ss(i, j)

                  IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                     dksc(na+i, nb+j, 1) = dsc(i, j, 1)
                     dksc(na+i, nb+j, 2) = dsc(i, j, 2)
                     dksc(na+i, nb+j, 3) = dsc(i, j, 3)

                     dkss(na+i, nb+j, 1) = dss(i, j, 1)
                     dkss(na+i, nb+j, 2) = dss(i, j, 2)
                     dkss(na+i, nb+j, 3) = dss(i, j, 3)
                  END IF

                  IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                     ddksc(na+i, nb+j, 1) = ddsc(i, j, 1)
                     ddksc(na+i, nb+j, 2) = ddsc(i, j, 2)
                     ddksc(na+i, nb+j, 3) = ddsc(i, j, 3)
                     ddksc(na+i, nb+j, 4) = ddsc(i, j, 4)
                     ddksc(na+i, nb+j, 5) = ddsc(i, j, 5)
                     ddksc(na+i, nb+j, 6) = ddsc(i, j, 6)

                     ddkss(na+i, nb+j, 1) = ddss(i, j, 1)
                     ddkss(na+i, nb+j, 2) = ddss(i, j, 2)
                     ddkss(na+i, nb+j, 3) = ddss(i, j, 3)
                     ddkss(na+i, nb+j, 4) = ddss(i, j, 4)
                     ddkss(na+i, nb+j, 5) = ddss(i, j, 5)
                     ddkss(na+i, nb+j, 6) = ddss(i, j, 6)
                  END IF

               END DO
            END DO

            IF (PRESENT(dcosab)) THEN
               la_start = 0
               lb_start = 0
            ELSE
               la_start = la_min
               lb_start = lb_min
            END IF

            DO da = 0, da_max-1
               ftz = 2.0_dp*zeta(ipgf)
               DO dax = 0, da
                  DO day = 0, da-dax
                     daz = da-dax-day
                     cda = coset(dax, day, daz)-1
                     cdax = coset(dax+1, day, daz)-1
                     cday = coset(dax, day+1, daz)-1
                     cdaz = coset(dax, day, daz+1)-1
                     !*** [da/dAi|O|b] = 2*zeta*[a+1i|O|b] - Ni(a)[a-1i|O|b] ***

                     DO la = la_start, la_max-da-1
                        DO ax = 0, la
                           fax = REAL(ax, dp)
                           DO ay = 0, la-ax
                              fay = REAL(ay, dp)
                              az = la-ax-ay
                              faz = REAL(az, dp)
                              coa = coset(ax, ay, az)
                              coamx = coset(ax-1, ay, az)
                              coamy = coset(ax, ay-1, az)
                              coamz = coset(ax, ay, az-1)
                              coapx = coset(ax+1, ay, az)
                              coapy = coset(ax, ay+1, az)
                              coapz = coset(ax, ay, az+1)
                              DO lb = lb_start, lb_max
                                 DO bx = 0, lb
                                    DO by = 0, lb-bx
                                       bz = lb-bx-by
                                       cob = coset(bx, by, bz)
                                       dscos(coa, cob, cdax) = ftz*sc(coapx, cob)-fax*sc(coamx, cob)
                                       dscos(coa, cob, cday) = ftz*sc(coapy, cob)-fay*sc(coamy, cob)
                                       dscos(coa, cob, cdaz) = ftz*sc(coapz, cob)-faz*sc(coamz, cob)
                                       dssin(coa, cob, cdax) = ftz*ss(coapx, cob)-fax*ss(coamx, cob)
                                       dssin(coa, cob, cday) = ftz*ss(coapy, cob)-fay*ss(coamy, cob)
                                       dssin(coa, cob, cdaz) = ftz*ss(coapz, cob)-faz*ss(coamz, cob)
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO

                  END DO
               END DO
            END DO

            IF (PRESENT(dcosab)) THEN
               DO k = 1, 3
                  DO j = 1, ncoset(lb_max)
                     DO i = 1, ncoset(la_max_set)
                        dcosab(na+i, nb+j, k) = dscos(i, j, k)
                        dsinab(na+i, nb+j, k) = dssin(i, j, k)
                     END DO
                  END DO
               END DO
            END IF

            nb = nb+ncoset(lb_max)

         END DO

         na = na+ncoset(la_max_set)

      END DO

   END SUBROUTINE cossin4

! **************************************************************************************************
!> \brief ...
!> \param la_max_set ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min_set ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param rac ...
!> \param rbc ...
!> \param kvec ...
!> \param cosab ...
!> \param sinab ...
!> \param dcosab ...
!> \param dsinab ...
!> \param dsc  First cosine derivative wrt k  
!> \param dss  First sine derivative wrt k 
!> \param ddsc Second cosine derivative wrt k
!> \param ddss Second sine derivative wrt k
! **************************************************************************************************
   SUBROUTINE cossin3(la_max_set, npgfa, zeta, rpgfa, la_min_set, &
                     lb_max, npgfb, zetb, rpgfb, lb_min, &
                     rac, rbc, kvec, cosab, sinab, dcosab, dsinab, &
                     dksc, dkss, ddksc, ddkss)

      INTEGER, INTENT(IN)                                :: la_max_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min_set, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc, kvec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: cosab, sinab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dcosab, dsinab
      !fbelle
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
                  OPTIONAL                               :: dksc, dkss, ddksc, ddkss

      INTEGER :: ax, ay, az, bx, by, bz, cda, cdax, cday, cdaz, coa, coamx, coamy, coamz, coapx, &
         coapy, coapz, cob, da, da_max, dax, day, daz, i, ipgf, j, jpgf, k, la, la_max, la_min, &
         la_start, lb, lb_start, na, nb
      REAL(KIND=dp)                                      :: dab, f0, f1, f2, f3, fax, fay, faz, ftz, &
                                                            fx, fy, fz, k2, kdp, rab2, s, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rap, rbp
      REAL(KIND=dp), DIMENSION(ncoset(la_max_set), &
         ncoset(lb_max), 3)                              :: dscos, dssin
      REAL(KIND=dp), &
         DIMENSION(ncoset(la_max_set+1), ncoset(lb_max)) :: sc, ss
      !fbelle
      REAL(KIND=dp), DIMENSION(3)                        :: r 
      REAL(KIND=dp), &
         DIMENSION(ncoset(la_max_set+1), ncoset(lb_max)) :: dsc, dss, ddsc, ddss

!      WRITE(*,*) "cossin3 - BEGINNING"

      rab = rbc-rac
      rab2 = SUM(rab**2)
      dab = SQRT(rab2)
      k2 = kvec(1)*kvec(1)+kvec(2)*kvec(2)+kvec(3)*kvec(3)

      IF (PRESENT(dcosab)) THEN
         da_max = 1
         la_max = la_max_set+1
         la_min = MAX(0, la_min_set-1)
         dscos = 0.0_dp
         dssin = 0.0_dp
      ELSE
         da_max = 0
         la_max = la_max_set
         la_min = la_min_set
      END IF

      ! initialize all matrix elements to zero
      IF (PRESENT(dcosab)) THEN
         na = ncoset(la_max-1)*npgfa
      ELSE
         na = ncoset(la_max)*npgfa
      END IF
      nb = ncoset(lb_max)*npgfb
      cosab(1:na, 1:nb) = 0.0_dp
      sinab(1:na, 1:nb) = 0.0_dp
      IF (PRESENT(dcosab)) THEN
         dcosab(1:na, 1:nb, :) = 0.0_dp
         dsinab(1:na, 1:nb, :) = 0.0_dp
      END IF

      ! first derivative 
      IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
         dksc(1:na, 1:nb) = 0.0_dp
         dkss(1:na, 1:nb) = 0.0_dp
      END IF
      ! second derivative
      IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
         ddksc(1:na, 1:nb) = 0.0_dp
         ddkss(1:na, 1:nb) = 0.0_dp
      END IF

!   *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0
      DO ipgf = 1, npgfa

         nb = 0

         DO jpgf = 1, npgfb

            ss = 0.0_dp
            sc = 0.0_dp

            IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
               dsc = 0.0_dp
               dss = 0.0_dp
            END IF
            IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
               ddsc = 0.0_dp
               ddss = 0.0_dp
            END IF

!       *** Screening ***
            IF (rpgfa(ipgf)+rpgfb(jpgf) < dab) THEN
               nb = nb+ncoset(lb_max)
               CYCLE
            END IF

!       *** Calculate some prefactors ***

            zetp = 1.0_dp/(zeta(ipgf)+zetb(jpgf))

            f0 = (pi*zetp)**1.5_dp
            f1 = zetb(jpgf)*zetp
            f2 = 0.5_dp*zetp

            kdp = zetp*DOT_PRODUCT(kvec, zeta(ipgf)*rac+zetb(jpgf)*rbc)

            r = zetp*(zeta(ipgf)*rac+zetb(jpgf)*rbc)

!       *** Calculate the basic two-center cos/sin integral [s|cos/sin|s] ***

            s = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*k2*zetp)
            sc(1, 1) = s*COS(kdp)
            ss(1, 1) = s*SIN(kdp)

            ! We introduce a derivation variable t 
            !
            ! kvec(:) -> t*kvec(:)
            ! kdp -> t*kdp
            ! k2 -> t^2*k2
            !
            ! d/dt (kvec(i)) = kvec(i)
            ! d/dt (s) = 2*k2*s*t
            ! d/dt (cos(t*kdp)) = (-1)*kdp*sin(t*kdp)
    
            ! d/dt (sc) = d/dt (s) * cos(kdp) + s*d/dt cos(kdp)
            !           = -0.25*2*t*k2*zetp*s*cos(kdp) + s*(-1)*kdp*sin(kdp)
            !           = -f2*t*k2*sc(1,1) - kdp*ss(1,1)

            IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
               dsc(1,1) = -f2*k2*sc(1,1) - kdp*ss(1,1) 
               dss(1,1) = -f2*k2*ss(1,1) + kdp*sc(1,1) 
            END IF

           ! d^2/dt^2 (sc) = -f2*sc(1,1) - f2*k2*dsc(1,1) - kdp*dss(1,1)

            IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
               ddsc(1,1) = -f2*k2*sc(1,1) - f2*k2*dsc(1,1) - kdp*dsc(1,1)
               ddss(1,1) = -f2*k2*ss(1,1) - f2*k2*dss(1,1) + kdp*dsc(1,1)
            END IF

!       *** Recurrence steps: [s|O|s] -> [a|O|b] ***

            IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|O|s] -> [a|O|s] ***

               rap(:) = f1*rab(:)

!         *** [p|O|s] = (Pi - Ai)*[s|O|s] +[s|dO|s]  (i = x,y,z) ***

               sc(2, 1) = rap(1)*sc(1, 1)-f2*kvec(1)*ss(1, 1)
               sc(3, 1) = rap(2)*sc(1, 1)-f2*kvec(2)*ss(1, 1)
               sc(4, 1) = rap(3)*sc(1, 1)-f2*kvec(3)*ss(1, 1)
               ss(2, 1) = rap(1)*ss(1, 1)+f2*kvec(1)*sc(1, 1)
               ss(3, 1) = rap(2)*ss(1, 1)+f2*kvec(2)*sc(1, 1)
               ss(4, 1) = rap(3)*ss(1, 1)+f2*kvec(3)*sc(1, 1)

               IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                  dsc(2, 1) = rap(1)*dsc(1, 1)-f2*ss(1, 1)-f2*kvec(1)*dss(1, 1)
                  dsc(3, 1) = rap(2)*dsc(1, 1)-f2*ss(1, 1)-f2*kvec(2)*dss(1, 1)
                  dsc(4, 1) = rap(3)*dsc(1, 1)-f2*ss(1, 1)-f2*kvec(3)*dss(1, 1)
                  dss(2, 1) = rap(1)*dss(1, 1)+f2*sc(1, 1)+f2*kvec(1)*dsc(1, 1)
                  dss(3, 1) = rap(2)*dss(1, 1)+f2*sc(1, 1)+f2*kvec(2)*dsc(1, 1)
                  dss(4, 1) = rap(3)*dss(1, 1)+f2*sc(1, 1)+f2*kvec(3)*dsc(1, 1)
               END IF

               IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                  ddsc(2, 1) = rap(1)*ddsc(1, 1)-2*f2*dss(1, 1)-f2*kvec(1)*ddss(1, 1)
                  ddsc(3, 1) = rap(2)*ddsc(1, 1)-2*f2*dss(1, 1)-f2*kvec(2)*ddss(1, 1)
                  ddsc(4, 1) = rap(3)*ddsc(1, 1)-2*f2*dss(1, 1)-f2*kvec(3)*ddss(1, 1)
                  ddss(2, 1) = rap(1)*ddss(1, 1)+2*f2*dsc(1, 1)+f2*kvec(1)*ddsc(1, 1)
                  ddss(3, 1) = rap(2)*ddss(1, 1)+2*f2*dsc(1, 1)+f2*kvec(2)*ddsc(1, 1)
                  ddss(4, 1) = rap(3)*ddss(1, 1)+2*f2*dsc(1, 1)+f2*kvec(3)*ddsc(1, 1)
               END IF

!         *** [a|O|s] = (Pi - Ai)*[a-1i|O|s] + f2*Ni(a-1i)*[a-2i|s] ***
!         ***           + [a-1i|dO|s]                               ***

               DO la = 2, la_max

!           *** Increase the angular momentum component z of function a ***

                  sc(coset(0, 0, la), 1) = rap(3)*sc(coset(0, 0, la-1), 1)+ &
                                           f2*REAL(la-1, dp)*sc(coset(0, 0, la-2), 1)- &
                                           f2*kvec(3)*ss(coset(0, 0, la-1), 1)
                  ss(coset(0, 0, la), 1) = rap(3)*ss(coset(0, 0, la-1), 1)+ &
                                           f2*REAL(la-1, dp)*ss(coset(0, 0, la-2), 1)+ &
                                           f2*kvec(3)*sc(coset(0, 0, la-1), 1)
           
                  IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                     dsc(coset(0, 0, la), 1) = rap(3)*dsc(coset(0, 0, la-1), 1)+ &
                                               f2*REAL(la-1, dp)*dsc(coset(0, 0, la-2), 1)- &
                                               f2*kvec(3)*ss(coset(0, 0, la-1), 1)- &
                                               f2*kvec(3)*dss(coset(0, 0, la-1), 1)
                     dss(coset(0, 0, la), 1) = rap(3)*dss(coset(0, 0, la-1), 1)+ &
                                               f2*REAL(la-1, dp)*dss(coset(0, 0, la-2), 1)+ &
                                               f2*kvec(3)*sc(coset(0, 0, la-1), 1)+ &
                                               f2*kvec(3)*dsc(coset(0, 0, la-1), 1)
                  END IF
                  
                  IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                     ddsc(coset(0, 0, la), 1) = rap(3)*ddsc(coset(0, 0, la-1), 1)+ &
                                                f2*REAL(la-1, dp)*ddsc(coset(0, 0, la-2), 1)- &
                                                f2*kvec(3)*ss(coset(0, 0, la-1), 1)- &
                                                2*f2*kvec(3)*dss(coset(0, 0, la-1), 1)- &
                                                f2*kvec(3)*ddss(coset(0, 0, la-1), 1)
                     ddss(coset(0, 0, la), 1) = rap(3)*ddss(coset(0, 0, la-1), 1)+ &
                                                f2*REAL(la-1, dp)*ddss(coset(0, 0, la-2), 1)+ &
                                                f2*kvec(3)*dsc(coset(0, 0, la-1), 1)+ &
                                                2*f2*kvec(3)*dsc(coset(0, 0, la-1), 1)+ &
                                                f2*kvec(3)*ddsc(coset(0, 0, la-1), 1)
                  END IF

!           *** Increase the angular momentum component y of function a ***

                  az = la-1
                  sc(coset(0, 1, az), 1) = rap(2)*sc(coset(0, 0, az), 1)- &
                                           f2*kvec(2)*ss(coset(0, 0, az), 1)
                  ss(coset(0, 1, az), 1) = rap(2)*ss(coset(0, 0, az), 1)+ &
                                           f2*kvec(2)*sc(coset(0, 0, az), 1)

                  IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                  
                     dsc(coset(0, 1, az), 1) = rap(2)*dsc(coset(0, 0, az), 1)- &
                                               f2*kvec(2)*ss(coset(0, 0, az), 1)- &
                                               f2*kvec(2)*dss(coset(0, 0, az), 1)
                     dss(coset(0, 1, az), 1) = rap(2)*dss(coset(0, 0, az), 1)+ &
                                               f2*kvec(2)*sc(coset(0, 0, az), 1)+ &
                                               f2*kvec(2)*dsc(coset(0, 0, az), 1)
                  END IF

                  IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                     ddsc(coset(0, 1, az), 1) = rap(2)*ddsc(coset(0, 0, az), 1)- &
                                                f2*kvec(2)*ss(coset(0, 0, az), 1)- &
                                                2*f2*kvec(2)*dss(coset(0, 0, az), 1)- &
                                                f2*kvec(2)*ddss(coset(0, 0, az), 1)
                     ddss(coset(0, 1, az), 1) = rap(2)*ddss(coset(0, 0, az), 1)+ &
                                                f2*kvec(2)*sc(coset(0, 0, az), 1)+ &
                                                2*f2*kvec(2)*dsc(coset(0, 0, az), 1)+ &
                                                f2*kvec(2)*dsc(coset(0, 0, az), 1)
                  END IF

                  DO ay = 2, la
                     az = la-ay
                     sc(coset(0, ay, az), 1) = rap(2)*sc(coset(0, ay-1, az), 1)+ &
                                               f2*REAL(ay-1, dp)*sc(coset(0, ay-2, az), 1)- &
                                               f2*kvec(2)*ss(coset(0, ay-1, az), 1)
                     ss(coset(0, ay, az), 1) = rap(2)*ss(coset(0, ay-1, az), 1)+ &
                                               f2*REAL(ay-1, dp)*ss(coset(0, ay-2, az), 1)+ &
                                               f2*kvec(2)*sc(coset(0, ay-1, az), 1)
        
                     IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                        dsc(coset(0, ay, az), 1) = rap(2)*dsc(coset(0, ay-1, az), 1)+ &
                                                   f2*REAL(ay-1, dp)*dsc(coset(0, ay-2, az), 1)- &
                                                   f2*kvec(2)*ss(coset(0, ay-1, az), 1)- &
                                                   f2*kvec(2)*dss(coset(0, ay-1, az), 1)
                        dss(coset(0, ay, az), 1) = rap(2)*dss(coset(0, ay-1, az), 1)+ &
                                                   f2*REAL(ay-1, dp)*dss(coset(0, ay-2, az), 1)+ &
                                                   f2*kvec(2)*sc(coset(0, ay-1, az), 1)+ &
                                                   f2*kvec(2)*dsc(coset(0, ay-1, az), 1)
                     END IF

                     IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                        ddsc(coset(0, ay, az), 1) = rap(2)*ddsc(coset(0, ay-1, az), 1)+ &
                                                    f2*REAL(ay-1, dp)*ddsc(coset(0, ay-2, az), 1)- &
                                                    f2*kvec(2)*ss(coset(0, ay-1, az), 1)- &
                                                    2*f2*dss(coset(0, ay-1, az), 1)- &
                                                    f2*kvec(2)*ddss(coset(0, ay-1, az), 1)
                        ddss(coset(0, ay, az), 1) = rap(2)*ddss(coset(0, ay-1, az), 1)+ &
                                                    f2*REAL(ay-1, dp)*ddss(coset(0, ay-2, az), 1)+ &
                                                    f2*kvec(2)*sc(coset(0, ay-1, az), 1)+ &
                                                    2*f2*kvec(2)*dsc(coset(0, ay-1, az), 1)+ &
                                                    f2*kvec(2)*ddsc(coset(0, ay-1, az), 1)
                     END IF
                  END DO

!           *** Increase the angular momentum component x of function a ***

                  DO ay = 0, la-1
                     az = la-1-ay
                     sc(coset(1, ay, az), 1) = rap(1)*sc(coset(0, ay, az), 1)- &
                                               f2*kvec(1)*ss(coset(0, ay, az), 1)
                     ss(coset(1, ay, az), 1) = rap(1)*ss(coset(0, ay, az), 1)+ &
                                               f2*kvec(1)*sc(coset(0, ay, az), 1)
                  
                     IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                        dsc(coset(1, ay, az), 1) = rap(1)*dsc(coset(0, ay, az), 1)- &
                                                   f2*kvec(1)*ss(coset(0, ay, az), 1)- &
                                                   f2*kvec(1)*dss(coset(0, ay, az), 1)
                        dss(coset(1, ay, az), 1) = rap(1)*dss(coset(0, ay, az), 1)+ &
                                                   f2*kvec(1)*sc(coset(0, ay, az), 1)+ &
                                                   f2*kvec(1)*dsc(coset(0, ay, az), 1)
                     END IF

                     IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                        ddsc(coset(1, ay, az), 1) = rap(1)*ddsc(coset(0, ay, az), 1)- &
                                                    f2*kvec(1)*dss(coset(0, ay, az), 1)- &
                                                    2*f2*kvec(1)*dss(coset(0, ay, az), 1)- &
                                                    f2*kvec(1)*ddss(coset(0, ay, az), 1)
                        ddss(coset(1, ay, az), 1) = rap(1)*ddss(coset(0, ay, az), 1)+ &
                                                    f2*kvec(1)*sc(coset(0, ay, az), 1)+ &
                                                    2*f2*kvec(1)*dsc(coset(0, ay, az), 1)+ &
                                                    f2*kvec(1)*ddsc(coset(0, ay, az), 1)
                     END IF
                                         
                  END DO

                  DO ax = 2, la
                     f3 = f2*REAL(ax-1, dp)
                     DO ay = 0, la-ax
                        az = la-ax-ay
                        sc(coset(ax, ay, az), 1) = rap(1)*sc(coset(ax-1, ay, az), 1)+ &
                                                   f3*sc(coset(ax-2, ay, az), 1)- &
                                                   f2*kvec(1)*ss(coset(ax-1, ay, az), 1)
                        ss(coset(ax, ay, az), 1) = rap(1)*ss(coset(ax-1, ay, az), 1)+ &
                                                   f3*ss(coset(ax-2, ay, az), 1)+ &
                                                   f2*kvec(1)*sc(coset(ax-1, ay, az), 1)

                        IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                           dsc(coset(ax, ay, az), 1) = rap(1)*dsc(coset(ax-1, ay, az), 1)+ &
                                                       f3*dsc(coset(ax-2, ay, az), 1)- &
                                                       f2*kvec(1)*ss(coset(ax-1, ay, az), 1)- &
                                                       f2*kvec(1)*dss(coset(ax-1, ay, az), 1)
                           dss(coset(ax, ay, az), 1) = rap(1)*dss(coset(ax-1, ay, az), 1)+ &
                                                       f3*dss(coset(ax-2, ay, az), 1)+ &
                                                       f2*kvec(1)*sc(coset(ax-1, ay, az), 1)+ &
                                                       f2*kvec(1)*dsc(coset(ax-1, ay, az), 1)
                        END IF
                        IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                           ddsc(coset(ax, ay, az), 1) = rap(1)*ddsc(coset(ax-1, ay, az), 1)+ &
                                                        f3*ddsc(coset(ax-2, ay, az), 1)- &
                                                        f2*kvec(1)*ss(coset(ax-1, ay, az), 1)- &
                                                        2*f2*kvec(1)*dss(coset(ax-1, ay, az), 1)- &
                                                        f2*kvec(1)*ddss(coset(ax-1, ay, az), 1)
                           ddss(coset(ax, ay, az), 1) = rap(1)*ddss(coset(ax-1, ay, az), 1)+ &
                                                        f3*ddss(coset(ax-2, ay, az), 1)+ &
                                                        f2*kvec(1)*sc(coset(ax-1, ay, az), 1)+ &
                                                        2*f2*kvec(1)*dsc(coset(ax-1, ay, az), 1)+ &
                                                        f2*kvec(1)*ddsc(coset(ax-1, ay, az), 1)
                        END IF

                     END DO
                  END DO

               END DO

!         *** Recurrence steps: [a|O|s] -> [a|O|b] ***

               IF (lb_max > 0) THEN

                  DO j = 2, ncoset(lb_max)
                     DO i = 1, ncoset(la_max)
                        sc(i, j) = 0.0_dp
                        ss(i, j) = 0.0_dp

                        IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                           dsc(i, j) = 0.0_dp
                           dss(i, j) = 0.0_dp
                        END IF
                        IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                           ddsc(i, j) = 0.0_dp
                           ddss(i, j) = 0.0_dp
                        END IF

                     END DO
                  END DO

!           *** Horizontal recurrence steps ***

                  rbp(:) = rap(:)-rab(:)

!           *** [a|O|p] = [a+1i|O|s] - (Bi - Ai)*[a|O|s] ***

                  IF (lb_max == 1) THEN
                     la_start = la_min
                  ELSE
                     la_start = MAX(0, la_min-1)
                  END IF

                  DO la = la_start, la_max-1
                     DO ax = 0, la
                        DO ay = 0, la-ax
                           az = la-ax-ay
                           sc(coset(ax, ay, az), 2) = sc(coset(ax+1, ay, az), 1)- &
                                                      rab(1)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 3) = sc(coset(ax, ay+1, az), 1)- &
                                                      rab(2)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 4) = sc(coset(ax, ay, az+1), 1)- &
                                                      rab(3)*sc(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = ss(coset(ax+1, ay, az), 1)- &
                                                      rab(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = ss(coset(ax, ay+1, az), 1)- &
                                                      rab(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = ss(coset(ax, ay, az+1), 1)- &
                                                      rab(3)*ss(coset(ax, ay, az), 1)

                        IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN

                           DO i = 1, 3
                              dsc(coset(ax, ay, az), 2) = dsc(coset(ax+1, ay, az), 1)- &
                                                          rab(1)*dsc(coset(ax, ay, az), 1)
                              dsc(coset(ax, ay, az), 3) = dsc(coset(ax, ay+1, az), 1)- &
                                                          rab(2)*dsc(coset(ax, ay, az), 1)
                              dsc(coset(ax, ay, az), 4) = dsc(coset(ax, ay, az+1), 1)- &
                                                          rab(3)*dsc(coset(ax, ay, az), 1)
                              dss(coset(ax, ay, az), 2) = dss(coset(ax+1, ay, az), 1)- &
                                                          rab(1)*dss(coset(ax, ay, az), 1)
                              dss(coset(ax, ay, az), 3) = dss(coset(ax, ay+1, az), 1)- &
                                                          rab(2)*dss(coset(ax, ay, az), 1)
                              dss(coset(ax, ay, az), 4) = dss(coset(ax, ay, az+1), 1)- &
                                                          rab(3)*dss(coset(ax, ay, az), 1)
                           END DO
                        END IF
                        IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                           DO i = 1, 6
                              ddsc(coset(ax, ay, az), 2) = ddsc(coset(ax+1, ay, az), 1)- &
                                                           rab(1)*ddsc(coset(ax, ay, az), 1)
                              ddsc(coset(ax, ay, az), 3) = ddsc(coset(ax, ay+1, az), 1)- &
                                                           rab(2)*ddsc(coset(ax, ay, az), 1)
                              ddsc(coset(ax, ay, az), 4) = ddsc(coset(ax, ay, az+1), 1)- &
                                                           rab(3)*ddsc(coset(ax, ay, az), 1)
                              ddss(coset(ax, ay, az), 2) = ddss(coset(ax+1, ay, az), 1)- &
                                                           rab(1)*ddss(coset(ax, ay, az), 1)
                              ddss(coset(ax, ay, az), 3) = ddss(coset(ax, ay+1, az), 1)- &
                                                           rab(2)*ddss(coset(ax, ay, az), 1)
                              ddss(coset(ax, ay, az), 4) = ddss(coset(ax, ay, az+1), 1)- &
                                                           rab(3)*ddss(coset(ax, ay, az), 1)
                           END DO
                        END IF

                        END DO
                     END DO
                  END DO

!           *** Vertical recurrence step ***

!           *** [a|O|p] = (Pi - Bi)*[a|O|s] + f2*Ni(a)*[a-1i|O|s] ***
!           ***           + [a|dO|s]                              ***

                  DO ax = 0, la_max
                     fx = f2*REAL(ax, dp)
                     DO ay = 0, la_max-ax
                        fy = f2*REAL(ay, dp)
                        az = la_max-ax-ay
                        fz = f2*REAL(az, dp)
                        IF (ax == 0) THEN
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)

                           IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN

                              dsc(coset(ax, ay, az), 2) = rbp(1)*dsc(coset(ax, ay, az), 1)- &
                                                          f2*kvec(1)*ss(coset(ax, ay, az), 1)- &
                                                          f2*kvec(1)*dss(coset(ax, ay, az), 1)
                              dss(coset(ax, ay, az), 2) = rbp(1)*dss(coset(ax, ay, az), 1)+ &
                                                          f2*kvec(1)*sc(coset(ax, ay, az), 1)+ &
                                                          f2*kvec(1)*dsc(coset(ax, ay, az), 1)
                           END IF
                           
                           IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN

                              ddsc(coset(ax, ay, az), 2) = rbp(1)*ddsc(coset(ax, ay, az), 1)- &
                                                           f2*kvec(1)*ss(coset(ax, ay, az), 1)- &
                                                           2*f2*kvec(1)*dss(coset(ax, ay, az), 1)- &
                                                           f2*kvec(1)*ddss(coset(ax, ay, az), 1)
                              ddss(coset(ax, ay, az), 2) = rbp(1)*dss(coset(ax, ay, az), 1)+ &
                                                           f2*kvec(1)*sc(coset(ax, ay, az), 1)+ &
                                                           2*f2*kvec(1)*dsc(coset(ax, ay, az), 1)+ &
                                                           f2*kvec(1)*ddsc(coset(ax, ay, az), 1)
                           END IF      
                        ELSE
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1)+ &
                                                      fx*sc(coset(ax-1, ay, az), 1)- &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1)+ &
                                                      fx*ss(coset(ax-1, ay, az), 1)+ &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)

                           IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN

                              dsc(coset(ax, ay, az), 2) = rbp(1)*dsc(coset(ax, ay, az), 1)+ &
                                                          fx*dsc(coset(ax-1, ay, az), 1)- &
                                                          f2*kvec(1)*ss(coset(ax, ay, az), 1)- &
                                                          f2*kvec(1)*dss(coset(ax, ay, az), 1)
                              dss(coset(ax, ay, az), 2) = rbp(1)*dss(coset(ax, ay, az), 1)+ &
                                                          fx*dss(coset(ax-1, ay, az), 1)+ &
                                                          f2*kvec(1)*sc(coset(ax, ay, az), 1)+ &
                                                          f2*kvec(1)*dsc(coset(ax, ay, az), 1)
                           END IF

                           IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                              ddsc(coset(ax, ay, az), 2) = rbp(1)*ddsc(coset(ax, ay, az), 1)+ &
                                                           fx*ddsc(coset(ax-1, ay, az), 1)- &
                                                           f2*kvec(1)*ss(coset(ax, ay, az), 1)- &
                                                           2*f2*kvec(1)*dss(coset(ax, ay, az), 1)- &
                                                           f2*kvec(1)*ddss(coset(ax, ay, az), 1)
                              ddss(coset(ax, ay, az), 2) = rbp(1)*ddss(coset(ax, ay, az), 1)+ &
                                                           fx*ddss(coset(ax-1, ay, az), 1)+ &
                                                           f2*kvec(1)*sc(coset(ax, ay, az), 1)+ &
                                                           2*f2*kvec(1)*dsc(coset(ax, ay, az), 1)+ &
                                                           f2*kvec(1)*ddsc(coset(ax, ay, az), 1)
                           END IF
                        END IF
                        IF (ay == 0) THEN
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                        
                           IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                              dsc(coset(ax, ay, az), 3) = rbp(2)*dsc(coset(ax, ay, az), 1)- &
                                                          f2*kvec(2)*ss(coset(ax, ay, az), 1)- &
                                                          f2*kvec(2)*dss(coset(ax, ay, az), 1)
                              dss(coset(ax, ay, az), 3) = rbp(2)*dss(coset(ax, ay, az), 1)+ &
                                                          f2*kvec(2)*sc(coset(ax, ay, az), 1)+ &
                                                          f2*kvec(2)*dsc(coset(ax, ay, az), 1)
                           END IF
                                                
                           IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                              ddsc(coset(ax, ay, az), 3) = rbp(2)*ddsc(coset(ax, ay, az), 1)- &
                                                           f2*kvec(2)*ss(coset(ax, ay, az), 1)- &
                                                           2*f2*kvec(2)*dss(coset(ax, ay, az), 1)- &
                                                           f2*kvec(2)*ddss(coset(ax, ay, az), 1)
                              ddss(coset(ax, ay, az), 3) = rbp(2)*ddss(coset(ax, ay, az), 1)+ &
                                                           f2*kvec(2)*sc(coset(ax, ay, az), 1)+ &
                                                           2*f2*kvec(2)*dsc(coset(ax, ay, az), 1)+ &
                                                           f2*kvec(2)*ddsc(coset(ax, ay, az), 1)
                           END IF
                                                
                        ELSE
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1)+ &
                                                      fy*sc(coset(ax, ay-1, az), 1)- &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1)+ &
                                                      fy*ss(coset(ax, ay-1, az), 1)+ &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                           
                           IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                              dsc(coset(ax, ay, az), 3) = rbp(2)*dsc(coset(ax, ay, az), 1)+ &
                                                          fy*dsc(coset(ax, ay-1, az), 1)- &
                                                          f2*kvec(2)*ss(coset(ax, ay, az), 1)- &
                                                          f2*kvec(2)*dss(coset(ax, ay, az), 1)
                              dss(coset(ax, ay, az), 3) = rbp(2)*dss(coset(ax, ay, az), 1)+ &
                                                          fy*dss(coset(ax, ay-1, az), 1)+ &
                                                          f2*kvec(2)*sc(coset(ax, ay, az), 1)+ &
                                                          f2*kvec(2)*dsc(coset(ax, ay, az), 1)
                           END IF
                                                
                           IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                              ddsc(coset(ax, ay, az), 3) = rbp(2)*dsc(coset(ax, ay, az), 1)+ &
                                                           fy*dsc(coset(ax, ay-1, az), 1)- &
                                                           f2*kvec(2)*ss(coset(ax, ay, az), 1)- &
                                                           2*f2*kvec(2)*dss(coset(ax, ay, az), 1)- &
                                                           f2*kvec(2)*ddss(coset(ax, ay, az), 1)
                              ddss(coset(ax, ay, az), 3) = rbp(2)*ddss(coset(ax, ay, az), 1)+ &
                                                           fy*ddss(coset(ax, ay-1, az), 1)+ &
                                                           f2*kvec(2)*sc(coset(ax, ay, az), 1)+ &
                                                           2*f2*kvec(2)*dsc(coset(ax, ay, az), 1)+ &
                                                           f2*kvec(2)*ddsc(coset(ax, ay, az), 1)
                           END IF
                                                
                        END IF
                        IF (az == 0) THEN
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)

                           IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                              dsc(coset(ax, ay, az), 4) = rbp(3)*dsc(coset(ax, ay, az), 1)- &
                                                          f2*kvec(3)*ss(coset(ax, ay, az), 1)- &
                                                          f2*kvec(3)*dss(coset(ax, ay, az), 1)
                              dss(coset(ax, ay, az), 4) = rbp(3)*dss(coset(ax, ay, az), 1)+ &
                                                          f2*kvec(3)*sc(coset(ax, ay, az), 1)+ &
                                                          f2*kvec(3)*dsc(coset(ax, ay, az), 1)

                           END IF
                           IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                              ddsc(coset(ax, ay, az), 4) = rbp(3)*ddsc(coset(ax, ay, az), 1)- &
                                                           f2*kvec(3)*ss(coset(ax, ay, az), 1)- &
                                                           2*f2*kvec(3)*dss(coset(ax, ay, az), 1)- &
                                                           f2*kvec(3)*ddss(coset(ax, ay, az), 1)
                              ddss(coset(ax, ay, az), 4) = rbp(3)*ddss(coset(ax, ay, az), 1)+ &
                                                           f2*kvec(3)*sc(coset(ax, ay, az), 1)+ &
                                                           2*f2*kvec(3)*dsc(coset(ax, ay, az), 1)+ &
                                                           f2*kvec(3)*ddsc(coset(ax, ay, az), 1)
                           END IF
                        
                        ELSE
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1)+ &
                                                      fz*sc(coset(ax, ay, az-1), 1)- &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1)+ &
                                                      fz*ss(coset(ax, ay, az-1), 1)+ &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)
                  
                                          
                           IF (PRESENT(dksc)) THEN                        
                              dsc(coset(ax, ay, az), 4) = rbp(3)*dsc(coset(ax, ay, az), 1)+ &
                                                          fz*dsc(coset(ax, ay, az-1), 1)- &
                                                          f2*kvec(3)*ss(coset(ax, ay, az), 1)- &
                                                          f2*kvec(3)*dss(coset(ax, ay, az), 1)
                              dss(coset(ax, ay, az), 4) = rbp(3)*dss(coset(ax, ay, az), 1)+ &
                                                          fz*dss(coset(ax, ay, az-1), 1)+ &
                                                          f2*kvec(3)*sc(coset(ax, ay, az), 1)+ &
                                                          f2*kvec(3)*dsc(coset(ax, ay, az), 1)
                           END IF

                           IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                              ddsc(coset(ax, ay, az), 4) = rbp(3)*ddsc(coset(ax, ay, az), 1)+ &
                                                           fz*ddsc(coset(ax, ay, az-1), 1)- &
                                                           f2*kvec(3)*ss(coset(ax, ay, az), 1)- &
                                                           2*f2*kvec(3)*dss(coset(ax, ay, az), 1)- &
                                                           f2*kvec(3)*ddss(coset(ax, ay, az), 1)
                              ddss(coset(ax, ay, az), 4) = rbp(3)*ddss(coset(ax, ay, az), 1)+ &
                                                           fz*ddss(coset(ax, ay, az-1), 1)+ &
                                                           f2*kvec(3)*sc(coset(ax, ay, az), 1)+ &
                                                           2*f2*kvec(3)*dsc(coset(ax, ay, az), 1)+ &
                                                           f2*kvec(3)*ddsc(coset(ax, ay, az), 1)
                           END IF
                        END IF
                     END DO
                  END DO

!           *** Recurrence steps: [a|O|p] -> [a|O|b] ***

                  DO lb = 2, lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|O|b] = [a+1i|O|b-1i] - (Bi - Ai)*[a|O|b-1i] ***

                     IF (lb == lb_max) THEN
                        la_start = la_min
                     ELSE
                        la_start = MAX(0, la_min-1)
                     END IF

                     DO la = la_start, la_max-1
                        DO ax = 0, la
                           DO ay = 0, la-ax
                              az = la-ax-ay

!                   *** Shift of angular momentum component z from a to b ***

                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 sc(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                 rab(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 ss(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                 rab(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))

                              IF (PRESENT(dksc)) THEN
                                 DO i = 1, 3
                                    dsc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                       dsc(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                       rab(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1))
                                    dss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                       dss(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                       rab(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1))
                                 END DO
                              END IF
                              IF (PRESENT(ddksc)) THEN
                                 DO i = 1, 6
                                    ddsc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                       ddsc(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                       rab(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1))
                                    ddss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                       ddss(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                       rab(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1))
                                 END DO
                              END IF

!                   *** Shift of angular momentum component y from a to b ***

                              DO by = 1, lb
                                 bz = lb-by
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    sc(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                    rab(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    ss(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                    rab(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))

                                 IF (PRESENT(dksc)) THEN
                                    DO i = 1, 3
                                       dsc(coset(ax, ay, az), coset(0, by, bz)) = &
                                          dsc(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                          rab(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz))
                                       dss(coset(ax, ay, az), coset(0, by, bz)) = &
                                          dss(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                          rab(2)*dss(coset(ax, ay, az), coset(0, by-1, bz))
                                    END DO
                                 END IF
                                 IF (PRESENT(ddksc)) THEN
                                    DO i = 1, 6
                                       ddsc(coset(ax, ay, az), coset(0, by, bz)) = &
                                          ddsc(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                          rab(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz))
                                       ddss(coset(ax, ay, az), coset(0, by, bz)) = &
                                          ddss(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                          rab(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz))
                                    END DO
                                 END IF
                              END DO

!                   *** Shift of angular momentum component x from a to b ***

                              DO bx = 1, lb
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       sc(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                       rab(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       ss(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                       rab(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))

                                    IF (PRESENT(dksc)) THEN
                                       DO i = 1, 3
                                          dsc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                              dsc(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                              rab(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz))
                                          dss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                              dss(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                              rab(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz))
                                       END DO
                                    END IF

                                    IF (PRESENT(ddksc)) THEN
                                       DO i = 1, 6
                                          ddsc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                               ddsc(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                               rab(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz))
                                          ddss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                               ddss(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                               rab(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz))
                                       END DO
                                    END IF

                                 END DO
                              END DO

                           END DO
                        END DO
                     END DO

!             *** Vertical recurrence step ***

!             *** [a|O|b] = (Pi - Bi)*[a|O|b-1i] + f2*Ni(a)*[a-1i|O|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|O|b-2i] + [a|dO|b-1i]            ***

                     DO ax = 0, la_max
                        fx = f2*REAL(ax, dp)
                        DO ay = 0, la_max-ax
                           fy = f2*REAL(ay, dp)
                           az = la_max-ax-ay
                           fz = f2*REAL(az, dp)

!                 *** Increase the angular momentum component z of function b ***

                           f3 = f2*REAL(lb-1, dp)

                           IF (az == 0) THEN
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              !------------------------------------------------
        
                              IF (PRESENT(dksc)) THEN
                                 dsc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                     rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                     f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                     f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1))
                                 dss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                     f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1))
                              END IF

                              IF (PRESENT(ddksc)) THEN
                                 ddsc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                      rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                      f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                      f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                      2*f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                      f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1))
                                 ddss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                      rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                      f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                      f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                      2*f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                      f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1))
                              END IF

                           ELSE
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 fz*sc(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 fz*ss(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              !------------------------------------------------

                              IF (PRESENT(dksc)) THEN
                                 dsc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                     rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     fz*dsc(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                     f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                     f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                     f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1))
                                 dss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     fz*dss(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                     f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1))
                              END IF
                              IF (PRESENT(ddksc)) THEN
                                 ddsc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                      rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                      fz*ddsc(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                      f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                      f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                      2*f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                      f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1))
                                 ddss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                      rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                      fz*ddss(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                      f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                      f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                      2*f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                      f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1))
                              END IF

                           END IF

!                 *** Increase the angular momentum component y of function b ***

                           IF (ay == 0) THEN
                              bz = lb-1
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz))- &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              !------------------------------------------------

                              IF (PRESENT(dksc)) THEN
                                 dsc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                     rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz))- &
                                     f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))- &
                                     f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz))
                                 dss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                     rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                     f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                     f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz))
                              END IF

                              IF (PRESENT(ddksc)) THEN
                                 ddsc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                      rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz))- &
                                      f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))- &
                                      2*f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz))- &
                                      f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz))
                                 ddss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                      rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                      f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                      2*f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                      f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz))
                              END IF

                              DO by = 2, lb
                                 bz = lb-by
                                 f3 = f2*REAL(by-1, dp)
                                 !------------------------------------------------
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    f3*sc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    f3*ss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 !------------------------------------------------
                              
                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, by, bz)) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                        f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz))
                                    dss(coset(ax, ay, az), coset(0, by, bz)) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                        f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz)) + &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz))
                                 END IF

                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(0, by, bz)) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                         f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                         2*f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz))
                                    ddss(coset(ax, ay, az), coset(0, by, bz)) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                         f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz)) + &
                                         2*f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz)) + &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz))
                                 END IF
                              
                              END DO
                           ELSE
                              bz = lb-1
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 fy*sc(coset(ax, ay-1, az), coset(0, 0, bz))- &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 fy*ss(coset(ax, ay-1, az), coset(0, 0, bz))+ &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              !------------------------------------------------

                              IF (PRESENT(dksc)) THEN
                                 dsc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                     rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                     fy*dsc(coset(ax, ay-1, az), coset(0, 0, bz))- &
                                     f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))- &
                                     f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz))
                                 dss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                     rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                     fy*dss(coset(ax, ay-1, az), coset(0, 0, bz))+ &
                                     f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                     f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz))
                              END IF

                              IF (PRESENT(ddksc)) THEN
                                 ddsc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                      rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                      fy*ddsc(coset(ax, ay-1, az), coset(0, 0, bz))- &
                                      f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))- &
                                      2*f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz))- &
                                      f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz))
                                 ddss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                      rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                      fy*ddss(coset(ax, ay-1, az), coset(0, 0, bz))+ &
                                      f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                      2*f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                      f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz))
                              END IF

                              DO by = 2, lb
                                 bz = lb-by
                                 f3 = f2*REAL(by-1, dp)
                                 !------------------------------------------------
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    fy*sc(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                    f3*sc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    fy*ss(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                    f3*ss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 !------------------------------------------------

                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, by, bz)) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                        f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz))
                                    dss(coset(ax, ay, az), coset(0, by, bz)) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                        f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz))
                                 END IF

                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(0, by, bz)) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                         f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                         2*f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz))
                                    ddss(coset(ax, ay, az), coset(0, by, bz)) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                         f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                         2*f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz))
                                 END IF

                              END DO
                           END IF

!                 *** Increase the angular momentum component x of function b ***

                           IF (ax == 0) THEN
                              DO by = 0, lb-1
                                 bz = lb-1-by
                                 !----------------------------------------
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz))- &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz))+ &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                                 !----------------------------------------
                              
                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(1, by, bz)) = &
                                        rbp(1)*dsc(coset(ax, ay, az), coset(0, by, bz))- &
                                        f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz))
                                    dss(coset(ax, ay, az), coset(1, by, bz)) = &
                                        rbp(1)*dss(coset(ax, ay, az), coset(0, by, bz))+ &
                                        f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz))
                                 END IF

                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(1, by, bz)) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), coset(0, by, bz))- &
                                         f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))- &
                                         2*f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz))- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz))
                                    ddss(coset(ax, ay, az), coset(1, by, bz)) = &
                                         rbp(1)*ddss(coset(ax, ay, az), coset(0, by, bz))+ &
                                         f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                         2*f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz))+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz))
                                 END IF
                              END DO
                              
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx-1, dp)
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    !----------------------------------------
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       f3*sc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       f3*ss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    !----------------------------------------

                                    IF (PRESENT(dksc)) THEN
                                       dsc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz))+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                           f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz))
                                       dss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz))+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                           f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    END IF
                                    
                                    IF (PRESENT(ddksc)) THEN
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                      coset(bx-1, by, bz))+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                            f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                            2*f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz))
                                       ddss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                      coset(bx-1, by, bz))+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                            f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                            2*f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    END IF
                                 END DO
                              END DO
                           ELSE
                              DO by = 0, lb-1
                                 bz = lb-1-by
                                 !----------------------------------------
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                    fx*sc(coset(ax-1, ay, az), coset(0, by, bz))- &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz))+ &
                                    fx*ss(coset(ax-1, ay, az), coset(0, by, bz))+ &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                                 !----------------------------------------
                              
                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(1, by, bz)) = &
                                        rbp(1)*dsc(coset(ax, ay, az), coset(0, by, bz))+ &
                                        fx*dsc(coset(ax-1, ay, az), coset(0, by, bz))- &
                                        f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz))
                                    dss(coset(ax, ay, az), coset(1, by, bz)) = &
                                        rbp(1)*dss(coset(ax, ay, az), coset(0, by, bz))+ &
                                        fx*dss(coset(ax-1, ay, az), coset(0, by, bz))+ &
                                        f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz))
                                 END IF

                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(1, by, bz)) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), coset(0, by, bz))+ &
                                         fx*ddsc(coset(ax-1, ay, az), coset(0, by, bz))- &
                                         f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))- &
                                         2*f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz))- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz))
                                    ddss(coset(ax, ay, az), coset(1, by, bz)) = &
                                         rbp(1)*ddss(coset(ax, ay, az), coset(0, by, bz))+ &
                                         fx*ddss(coset(ax-1, ay, az), coset(0, by, bz))+ &
                                         2*f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz))+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz))
                                 END IF

                              END DO
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx-1, dp)
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    !----------------------------------------
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       fx*sc(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                       f3*sc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       fx*ss(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                       f3*ss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    !----------------------------------------
                                 
                                    IF (PRESENT(dksc)) THEN
                                       dsc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                      coset(bx-1, by, bz))+ &
                                           fx*dsc(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                           f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz))
                                       dss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                      coset(bx-1, by, bz))+ &
                                           fx*dss(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                           f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    END IF
                                    
                                    IF (PRESENT(ddksc)) THEN
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz))+ &
                                            fx*ddsc(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                            f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                            2*f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz))
                                       ddss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz))+ &
                                            fx*ddss(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                            f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                            2*f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    END IF
                                 
                                 END DO
                              END DO
                           END IF

                        END DO
                     END DO

                  END DO

               END IF

            ELSE

               IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|O|s] -> [s|O|b] ***

                  rbp(:) = (f1-1.0_dp)*rab(:)

!           *** [s|O|p] = (Pi - Bi)*[s|O|s] + [s|dO|s] ***

                  !----------------------------------------
                  sc(1, 2) = rbp(1)*sc(1, 1)-f2*kvec(1)*ss(1, 1)
                  sc(1, 3) = rbp(2)*sc(1, 1)-f2*kvec(2)*ss(1, 1)
                  sc(1, 4) = rbp(3)*sc(1, 1)-f2*kvec(3)*ss(1, 1)
                  ss(1, 2) = rbp(1)*ss(1, 1)+f2*kvec(1)*sc(1, 1)
                  ss(1, 3) = rbp(2)*ss(1, 1)+f2*kvec(2)*sc(1, 1)
                  ss(1, 4) = rbp(3)*ss(1, 1)+f2*kvec(3)*sc(1, 1)
                  !----------------------------------------

               IF (PRESENT(dksc)) THEN
                  dsc(1,2) = rbp(1)*dsc(1,1) - f2*kvec(1)*ss(1,1) - f2*kvec(1)*dss(1,1) 
                  dsc(1,3) = rbp(2)*dsc(1,1) - f2*kvec(2)*ss(1,1) - f2*kvec(2)*dss(1,1)
                  dsc(1,4) = rbp(3)*dsc(1,1) - f2*kvec(3)*ss(1,1) - f2*kvec(3)*dss(1,1)
                  dss(1,2) = rbp(1)*dss(1,1) + f2*kvec(1)*sc(1,1) + f2*kvec(1)*dsc(1,1) 
                  dss(1,3) = rbp(2)*dss(1,1) + f2*kvec(2)*sc(1,1) + f2*kvec(2)*dsc(1,1)
                  dss(1,4) = rbp(3)*dss(1,1) + f2*kvec(3)*sc(1,1) + f2*kvec(3)*dsc(1,1)
               END IF

               IF (PRESENT(ddksc)) THEN
                  ddsc(1,2) = rbp(1)*ddsc(1,1)-f2*kvec(1)*ss(1,1)-2*f2*kvec(1)*dss(1,1)-f2*kvec(1)*ddss(1,1) 
                  ddsc(1,3) = rbp(2)*ddsc(1,1)-f2*kvec(2)*ss(1,1)-2*f2*kvec(2)*dss(1,1)-f2*kvec(1)*ddss(1,1)
                  ddsc(1,4) = rbp(3)*ddsc(1,1)-f2*kvec(3)*ss(1,1)-2*f2*kvec(3)*dss(1,1)-f2*kvec(1)*ddss(1,1)
                  ddss(1,2) = rbp(1)*ddss(1,1)+f2*kvec(1)*sc(1,1)+2*f2*kvec(1)*dsc(1,1)+f2*kvec(1)*ddss(1,1) 
                  ddss(1,3) = rbp(2)*ddss(1,1)+f2*kvec(2)*sc(1,1)+2*f2*kvec(2)*dsc(1,1)+f2*kvec(2)*ddss(1,1)
                  ddss(1,4) = rbp(3)*ddss(1,1)+f2*kvec(3)*sc(1,1)+2*f2*kvec(3)*dsc(1,1)+f2*kvec(3)*ddss(1,1)
               END IF

! fbelle changed until here


!           *** [s|O|b] = (Pi - Bi)*[s|O|b-1i] + f2*Ni(b-1i)*[s|O|b-2i] ***
!           ***           + [s|dO|b-1i]                                 ***

                  DO lb = 2, lb_max

!             *** Increase the angular momentum component z of function b ***

                     !-------------------------------------
                     sc(1, coset(0, 0, lb)) = rbp(3)*sc(1, coset(0, 0, lb-1))+ &
                                              f2*REAL(lb-1, dp)*sc(1, coset(0, 0, lb-2))- &
                                              f2*kvec(3)*ss(1, coset(0, 0, lb-1))
                     ss(1, coset(0, 0, lb)) = rbp(3)*ss(1, coset(0, 0, lb-1))+ &
                                              f2*REAL(lb-1, dp)*ss(1, coset(0, 0, lb-2))+ &
                                              f2*kvec(3)*sc(1, coset(0, 0, lb-1))
                     !-------------------------------------
              
                     IF (PRESENT(dksc)) THEN
                        dsc(1, coset(0, 0, lb)) = rbp(3)*dsc(1, coset(0, 0, lb-1))+ &
                                                  f2*REAL(lb-1, dp)*dsc(1, coset(0, 0, lb-2))- &
                                                  f2*kvec(3)*ss(1, coset(0, 0, lb-1))- &
                                                  f2*kvec(3)*dss(1, coset(0, 0, lb-1))
                        dss(1, coset(0, 0, lb)) = rbp(3)*dss(1, coset(0, 0, lb-1))+ &
                                                  f2*REAL(lb-1, dp)*dss(1, coset(0, 0, lb-2))+ &
                                                  f2*kvec(3)*sc(1, coset(0, 0, lb-1))+ &
                                                  f2*kvec(3)*dsc(1, coset(0, 0, lb-1))
                     END IF

                     IF (PRESENT(ddksc)) THEN
                        ddsc(1, coset(0, 0, lb)) = rbp(3)*ddsc(1, coset(0, 0, lb-1))+ &
                                                   f2*REAL(lb-1, dp)*ddsc(1, coset(0, 0, lb-2))- &
                                                   f2*kvec(3)*ss(1, coset(0, 0, lb-1))- &
                                                   2*f2*kvec(3)*dss(1, coset(0, 0, lb-1))- &
                                                   f2*kvec(3)*ddss(1, coset(0, 0, lb-1))
                        ddss(1, coset(0, 0, lb)) = rbp(3)*ddss(1, coset(0, 0, lb-1))+ &
                                                   f2*REAL(lb-1, dp)*ddss(1, coset(0, 0, lb-2))+ &
                                                   f2*kvec(3)*sc(1, coset(0, 0, lb-1))+ &
                                                   2*f2*kvec(3)*dsc(1, coset(0, 0, lb-1))+ &
                                                   f2*kvec(3)*ddsc(1, coset(0, 0, lb-1))
                     END IF

!             *** Increase the angular momentum component y of function b ***

                     bz = lb-1
                     !-------------------------------------
                     sc(1, coset(0, 1, bz)) = rbp(2)*sc(1, coset(0, 0, bz))- &
                                              f2*kvec(2)*ss(1, coset(0, 0, bz))
                     ss(1, coset(0, 1, bz)) = rbp(2)*ss(1, coset(0, 0, bz))+ &
                                              f2*kvec(2)*sc(1, coset(0, 0, bz))
                     !-------------------------------------

                     IF (PRESENT(dksc)) THEN
                        dsc(1, coset(0, 1, bz)) = rbp(2)*dsc(1, coset(0, 0, bz))- &
                                                  f2*kvec(2)*ss(1, coset(0, 0, bz))- &
                                                  f2*kvec(2)*dss(1, coset(0, 0, bz))
                        dss(1, coset(0, 1, bz)) = rbp(2)*dss(1, coset(0, 0, bz))+ &
                                                  f2*kvec(2)*sc(1, coset(0, 0, bz))+ &
                                                  f2*kvec(2)*dsc(1, coset(0, 0, bz))
                     END IF

                     IF (PRESENT(ddksc)) THEN
                        ddsc(1, coset(0, 1, bz)) = rbp(2)*ddsc(1, coset(0, 0, bz))- &
                                                   f2*kvec(2)*ss(1, coset(0, 0, bz))- &
                                                   2*f2*kvec(2)*dss(1, coset(0, 0, bz))- &
                                                   f2*kvec(2)*ddss(1, coset(0, 0, bz))
                        ddss(1, coset(0, 1, bz)) = rbp(2)*ddss(1, coset(0, 0, bz))+ &
                                                   f2*kvec(2)*sc(1, coset(0, 0, bz))+ &
                                                   2*f2*kvec(2)*dsc(1, coset(0, 0, bz))+ &
                                                   f2*kvec(2)*ddsc(1, coset(0, 0, bz))
                     END IF

                     DO by = 2, lb
                        bz = lb-by
                        !-------------------------------------
                        sc(1, coset(0, by, bz)) = rbp(2)*sc(1, coset(0, by-1, bz))+ &
                                                  f2*REAL(by-1, dp)*sc(1, coset(0, by-2, bz))- &
                                                  f2*kvec(2)*ss(1, coset(0, by-1, bz))
                        ss(1, coset(0, by, bz)) = rbp(2)*ss(1, coset(0, by-1, bz))+ &
                                                  f2*REAL(by-1, dp)*ss(1, coset(0, by-2, bz))+ &
                                                  f2*kvec(2)*sc(1, coset(0, by-1, bz))
                        !-------------------------------------
                    
                        IF (PRESENT(dksc)) THEN
                           dsc(1, coset(0, by, bz)) = rbp(2)*dsc(1, coset(0, by-1, bz))+ &
                                                      f2*REAL(by-1, dp)*dsc(1, coset(0, by-2, bz))- &
                                                      f2*kvec(2)*ss(1, coset(0, by-1, bz))- &
                                                      f2*kvec(2)*dss(1, coset(0, by-1, bz))
                           dss(1, coset(0, by, bz)) = rbp(2)*dss(1, coset(0, by-1, bz))+ &
                                                      f2*REAL(by-1, dp)*dss(1, coset(0, by-2, bz))+ &
                                                      f2*kvec(2)*sc(1, coset(0, by-1, bz))+ &
                                                      f2*kvec(2)*dsc(1, coset(0, by-1, bz))
                        END IF

                        IF (PRESENT(ddksc)) THEN
                           ddsc(1, coset(0, by, bz)) = rbp(2)*ddsc(1, coset(0, by-1, bz))+ &
                                                       f2*REAL(by-1, dp)*ddsc(1, coset(0, by-2, bz))- &
                                                       f2*kvec(2)*ss(1, coset(0, by-1, bz))- &
                                                       2*f2*kvec(2)*dss(1, coset(0, by-1, bz))- &
                                                       f2*kvec(2)*ddss(1, coset(0, by-1, bz))
                           ddss(1, coset(0, by, bz)) = rbp(2)*ddss(1, coset(0, by-1, bz))+ &
                                                       f2*REAL(by-1, dp)*ddss(1, coset(0, by-2, bz))+ &
                                                       f2*kvec(2)*sc(1, coset(0, by-1, bz))+ &
                                                       2*f2*kvec(2)*dsc(1, coset(0, by-1, bz))+ &
                                                       f2*kvec(2)*ddsc(1, coset(0, by-1, bz))
                        END IF
                    
                     END DO

!             *** Increase the angular momentum component x of function b ***

                     DO by = 0, lb-1
                        bz = lb-1-by
                        !----------------------------------
                        sc(1, coset(1, by, bz)) = rbp(1)*sc(1, coset(0, by, bz))- &
                                                  f2*kvec(1)*ss(1, coset(0, by, bz))
                        ss(1, coset(1, by, bz)) = rbp(1)*ss(1, coset(0, by, bz))+ &
                                                  f2*kvec(1)*sc(1, coset(0, by, bz))
                        !----------------------------------
                                               
                        IF (PRESENT(dksc)) THEN
                           dsc(1, coset(1, by, bz)) = rbp(1)*dsc(1, coset(0, by, bz))- &
                                                      f2*kvec(1)*ss(1, coset(0, by, bz))- &
                                                      f2*kvec(1)*dss(1, coset(0, by, bz))
                           dss(1, coset(1, by, bz)) = rbp(1)*dss(1, coset(0, by, bz))+ &
                                                      f2*kvec(1)*sc(1, coset(0, by, bz))+ &
                                                      f2*kvec(1)*dsc(1, coset(0, by, bz))
                        END IF                       

                        IF (PRESENT(ddksc)) THEN
                           ddsc(1, coset(1, by, bz)) = rbp(1)*ddsc(1, coset(0, by, bz))- &
                                                       f2*kvec(1)*ss(1, coset(0, by, bz))- &
                                                       2*f2*kvec(1)*dss(1, coset(0, by, bz))- &
                                                       f2*kvec(1)*ddss(1, coset(0, by, bz))
                           ddss(1, coset(1, by, bz)) = rbp(1)*ddss(1, coset(0, by, bz))+ &
                                                       f2*kvec(1)*sc(1, coset(0, by, bz))+ &
                                                       2*f2*kvec(1)*dsc(1, coset(0, by, bz))+ &
                                                       f2*kvec(1)*ddsc(1, coset(0, by, bz))
                        END IF                       

                     END DO

                     DO bx = 2, lb
                        f3 = f2*REAL(bx-1, dp)
                        DO by = 0, lb-bx
                           bz = lb-bx-by
                           !----------------------------------
                           sc(1, coset(bx, by, bz)) = rbp(1)*sc(1, coset(bx-1, by, bz))+ &
                                                      f3*sc(1, coset(bx-2, by, bz))- &
                                                      f2*kvec(1)*ss(1, coset(bx-1, by, bz))
                           ss(1, coset(bx, by, bz)) = rbp(1)*ss(1, coset(bx-1, by, bz))+ &
                                                      f3*ss(1, coset(bx-2, by, bz))+ &
                                                      f2*kvec(1)*sc(1, coset(bx-1, by, bz))
                          !----------------------------------
                                                   
                           IF (PRESENT(dksc)) THEN
                              dsc(1, coset(bx, by, bz)) = rbp(1)*dsc(1, coset(bx-1, by, bz))+ &
                                                          f3*dsc(1, coset(bx-2, by, bz))- &
                                                          f2*kvec(1)*ss(1, coset(bx-1, by, bz))- &
                                                          f2*kvec(1)*dss(1, coset(bx-1, by, bz))
                              dss(1, coset(bx, by, bz)) = rbp(1)*dss(1, coset(bx-1, by, bz))+ &
                                                          f3*dss(1, coset(bx-2, by, bz))+ &
                                                          f2*kvec(1)*ss(1, coset(bx-1, by, bz))+ &
                                                          f2*kvec(1)*dss(1, coset(bx-1, by, bz))
                           END IF                       
                                                   
                           IF (PRESENT(ddksc)) THEN
                              ddsc(1, coset(bx, by, bz)) = rbp(1)*ddsc(1, coset(bx-1, by, bz))+ &
                                                           f3*ddsc(1, coset(bx-2, by, bz))- &
                                                           f2*kvec(1)*ss(1, coset(bx-1, by, bz))- &
                                                           2*f2*kvec(1)*dss(1, coset(bx-1, by, bz))- &
                                                           f2*kvec(1)*ddss(1, coset(bx-1, by, bz))
                              ddss(1, coset(bx, by, bz)) = rbp(1)*ddss(1, coset(bx-1, by, bz))+ &
                                                           f3*ddss(1, coset(bx-2, by, bz))+ &
                                                           f2*kvec(1)*ss(1, coset(bx-1, by, bz))- &
                                                           2*f2*kvec(1)*dss(1, coset(bx-1, by, bz))- &
                                                           f2*kvec(1)*ddss(1, coset(bx-1, by, bz))
                           END IF                       
                                                   
                        END DO
                     END DO

                  END DO

               END IF

            END IF

            DO j = ncoset(lb_min-1)+1, ncoset(lb_max)
               DO i = ncoset(la_min_set-1)+1, ncoset(la_max_set)
                  cosab(na+i, nb+j) = sc(i, j)
                  sinab(na+i, nb+j) = ss(i, j)

                  IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                     dksc(na+i, nb+j) = dsc(i, j)
                     dkss(na+i, nb+j) = dss(i, j)
                  END IF

                  IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                     ddksc(na+i, nb+j) = ddsc(i, j)
                     ddkss(na+i, nb+j) = ddss(i, j)
                  END IF

               END DO
            END DO

            IF (PRESENT(dcosab)) THEN
               la_start = 0
               lb_start = 0
            ELSE
               la_start = la_min
               lb_start = lb_min
            END IF

            DO da = 0, da_max-1
               ftz = 2.0_dp*zeta(ipgf)
               DO dax = 0, da
                  DO day = 0, da-dax
                     daz = da-dax-day
                     cda = coset(dax, day, daz)-1
                     cdax = coset(dax+1, day, daz)-1
                     cday = coset(dax, day+1, daz)-1
                     cdaz = coset(dax, day, daz+1)-1
                     !*** [da/dAi|O|b] = 2*zeta*[a+1i|O|b] - Ni(a)[a-1i|O|b] ***

                     DO la = la_start, la_max-da-1
                        DO ax = 0, la
                           fax = REAL(ax, dp)
                           DO ay = 0, la-ax
                              fay = REAL(ay, dp)
                              az = la-ax-ay
                              faz = REAL(az, dp)
                              coa = coset(ax, ay, az)
                              coamx = coset(ax-1, ay, az)
                              coamy = coset(ax, ay-1, az)
                              coamz = coset(ax, ay, az-1)
                              coapx = coset(ax+1, ay, az)
                              coapy = coset(ax, ay+1, az)
                              coapz = coset(ax, ay, az+1)
                              DO lb = lb_start, lb_max
                                 DO bx = 0, lb
                                    DO by = 0, lb-bx
                                       bz = lb-bx-by
                                       cob = coset(bx, by, bz)
                                       dscos(coa, cob, cdax) = ftz*sc(coapx, cob)-fax*sc(coamx, cob)
                                       dscos(coa, cob, cday) = ftz*sc(coapy, cob)-fay*sc(coamy, cob)
                                       dscos(coa, cob, cdaz) = ftz*sc(coapz, cob)-faz*sc(coamz, cob)
                                       dssin(coa, cob, cdax) = ftz*ss(coapx, cob)-fax*ss(coamx, cob)
                                       dssin(coa, cob, cday) = ftz*ss(coapy, cob)-fay*ss(coamy, cob)
                                       dssin(coa, cob, cdaz) = ftz*ss(coapz, cob)-faz*ss(coamz, cob)
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO

                  END DO
               END DO
            END DO

            IF (PRESENT(dcosab)) THEN
               DO k = 1, 3
                  DO j = 1, ncoset(lb_max)
                     DO i = 1, ncoset(la_max_set)
                        dcosab(na+i, nb+j, k) = dscos(i, j, k)
                        dsinab(na+i, nb+j, k) = dssin(i, j, k)
                     END DO
                  END DO
               END DO
            END IF

            nb = nb+ncoset(lb_max)

         END DO

         na = na+ncoset(la_max_set)

      END DO

   END SUBROUTINE cossin3

! **************************************************************************************************
!> \brief ...                   FABIAN derivation WRT k, which is wrong, but since 3 derivs, 
!>                              it's ok in terms of component structure
!> \param la_max_set ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min_set ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param rac ...
!> \param rbc ...
!> \param kvec ...
!> \param cosab ...
!> \param sinab ...
!> \param dcosab ...
!> \param dsinab ...
!> \param dsc  First cosine derivative wrt k  
!> \param dss  First sine derivative wrt k 
!> \param ddsc Second cosine derivative wrt k
!> \param ddss Second sine derivative wrt k
! **************************************************************************************************
   SUBROUTINE cossin2(la_max_set, npgfa, zeta, rpgfa, la_min_set, &
                     lb_max, npgfb, zetb, rpgfb, lb_min, &
                     rac, rbc, kvec, cosab, sinab, dcosab, dsinab, &
                     dksc, dkss, ddksc, ddkss)

      INTEGER, INTENT(IN)                                :: la_max_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min_set, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc, kvec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: cosab, sinab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dcosab, dsinab
      !fbelle
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
                  OPTIONAL                               :: dksc, dkss, ddksc, ddkss

      INTEGER :: ax, ay, az, bx, by, bz, cda, cdax, cday, cdaz, coa, coamx, coamy, coamz, coapx, &
         coapy, coapz, cob, da, da_max, dax, day, daz, i, ipgf, j, jpgf, k, la, la_max, la_min, &
         la_start, lb, lb_start, na, nb
      REAL(KIND=dp)                                      :: dab, f0, f1, f2, f3, fax, fay, faz, ftz, &
                                                            fx, fy, fz, k2, kdp, rab2, s, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rap, rbp
      REAL(KIND=dp), DIMENSION(ncoset(la_max_set), &
         ncoset(lb_max), 3)                              :: dscos, dssin
      REAL(KIND=dp), &
         DIMENSION(ncoset(la_max_set+1), ncoset(lb_max)) :: sc, ss
      !fbelle
      REAL(KIND=dp), DIMENSION(3)                        :: r 
      REAL(KIND=dp), DIMENSION(ncoset(la_max_set+1), ncoset(lb_max), 3) :: dsc, dss
      REAL(KIND=dp), DIMENSION(ncoset(la_max_set+1), ncoset(lb_max), 6) :: ddsc, ddss

!      WRITE(*,*) "cossin2 - BEGINNING"

      rab = rbc-rac
      rab2 = SUM(rab**2)
      dab = SQRT(rab2)
      k2 = kvec(1)*kvec(1)+kvec(2)*kvec(2)+kvec(3)*kvec(3)

      IF (PRESENT(dcosab)) THEN
         da_max = 1
         la_max = la_max_set+1
         la_min = MAX(0, la_min_set-1)
         dscos = 0.0_dp
         dssin = 0.0_dp
      ELSE
         da_max = 0
         la_max = la_max_set
         la_min = la_min_set
      END IF

      ! initialize all matrix elements to zero
      IF (PRESENT(dcosab)) THEN
         na = ncoset(la_max-1)*npgfa
      ELSE
         na = ncoset(la_max)*npgfa
      END IF
      nb = ncoset(lb_max)*npgfb
      cosab(1:na, 1:nb) = 0.0_dp
      sinab(1:na, 1:nb) = 0.0_dp
      IF (PRESENT(dcosab)) THEN
         dcosab(1:na, 1:nb, :) = 0.0_dp
         dsinab(1:na, 1:nb, :) = 0.0_dp
      END IF

      ! first derivative 
      IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
         DO i = 1, 3
            dksc(1:na, 1:nb, i) = 0.0_dp
            dkss(1:na, 1:nb, i) = 0.0_dp
         END DO
      END IF
      ! second derivative
      IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
         DO i = 1, 6
            ddksc(1:na, 1:nb, i) = 0.0_dp
            ddkss(1:na, 1:nb, i) = 0.0_dp
         END DO
      END IF

!   *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0
      DO ipgf = 1, npgfa

         nb = 0

         DO jpgf = 1, npgfb

            ss = 0.0_dp
            sc = 0.0_dp

            IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
               dsc = 0.0_dp
               dss = 0.0_dp
            END IF
            IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
               ddsc = 0.0_dp
               ddss = 0.0_dp
            END IF

!       *** Screening ***
            IF (rpgfa(ipgf)+rpgfb(jpgf) < dab) THEN
               nb = nb+ncoset(lb_max)
               CYCLE
            END IF

!       *** Calculate some prefactors ***

            zetp = 1.0_dp/(zeta(ipgf)+zetb(jpgf))

            f0 = (pi*zetp)**1.5_dp
            f1 = zetb(jpgf)*zetp
            f2 = 0.5_dp*zetp

            kdp = zetp*DOT_PRODUCT(kvec, zeta(ipgf)*rac+zetb(jpgf)*rbc)

            r = zetp*(zeta(ipgf)*rac+zetb(jpgf)*rbc)

!       *** Calculate the basic two-center cos/sin integral [s|cos/sin|s] ***

            s = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*k2*zetp)
            sc(1, 1) = s*COS(kdp)
            ss(1, 1) = s*SIN(kdp)

            IF (PRESENT(dksc)) THEN
               ! dsc/dk1
               dsc(1,1,1) = -f2*kvec(1)*sc(1,1) - r(1)*ss(1,1) 
               ! dsc/dk2
               dsc(1,1,2) = -f2*kvec(2)*sc(1,1) - r(2)*ss(1,1) 
               ! dsc/dk3
               dsc(1,1,3) = -f2*kvec(3)*sc(1,1) - r(3)*ss(1,1)
            
            END IF
            IF (PRESENT(dkss)) THEN
               ! dss/dk1
               dss(1,1,1) = -f2*kvec(1)*ss(1,1) + r(1)*sc(1,1) 
               ! dss/dk2
               dss(1,1,2) = -f2*kvec(2)*ss(1,1) + r(2)*sc(1,1) 
               ! dss/dk3
               dss(1,1,3) = -f2*kvec(3)*ss(1,1) + r(3)*sc(1,1) 
            END IF

            IF (PRESENT(ddksc)) THEN
               ! d(dsc(1,1,1)) / dk1
               ddsc(1,1,1) = -f2*sc(1,1) - f2*kvec(1)*dsc(1,1,1) - r(1)*dss(1,1,1)
               ! d(dsc(1,1,1)) / dk2
               ddsc(1,1,2) = -f2*kvec(1)*dsc(1,1,2) - r(1)*dss(1,1,2)
               ! d(dsc(1,1,1)) / dk3
               ddsc(1,1,3) = -f2*kvec(1)*dsc(1,1,3) - r(1)*dss(1,1,3)
               ! d(dsc(1,1,2)) / dk2
               ddsc(1,1,4) = -f2*sc(1,1) - f2*kvec(2)*dsc(1,1,2) - r(2)*dss(1,1,2) 
               ! d(dsc(1,1,2)) / dk3
               ddsc(1,1,5) = -f2*kvec(2)*dsc(1,1,3) - r(2)*dss(1,1,3)
               ! d(dsc(1,1,3)) / dk3
               ddsc(1,1,6) = -f2*sc(1,1) - f2*kvec(3)*dsc(1,1,3) - r(3)*dss(1,1,3)
            END IF

            IF (PRESENT(ddkss)) THEN
               ! d(dss(1,1,1)) / dk1
               ddss(1,1,1) = -f2*ss(1,1) - f2*kvec(1)*dss(1,1,1) + r(1)*dsc(1,1,1)
               ! d(dss(1,1,1)) / dk2
               ddss(1,1,2) = -f2*kvec(1)*dss(1,1,2) + r(1)*dsc(1,1,2)
               ! d(dss(1,1,1)) / dk3
               ddss(1,1,3) = -f2*kvec(1)*dss(1,1,3) + r(1)*dsc(1,1,3)
               ! d(dss(1,1,2)) / dk2
               ddss(1,1,4) = -f2*ss(1,1) - f2*kvec(2)*dss(1,1,2) + r(2)*dsc(1,1,2) 
               ! d(dss(1,1,2)) / dk3
               ddss(1,1,5) = -f2*kvec(2)*dss(1,1,3) + r(2)*dsc(1,1,3)
               ! d(dss(1,1,3)) / dk3
               ddss(1,1,6) = -f2*ss(1,1) - f2*kvec(3)*dss(1,1,3) + r(3)*dsc(1,1,3)
            END IF

!       *** Recurrence steps: [s|O|s] -> [a|O|b] ***

            IF (la_max > 0) THEN
!         *** Vertical recurrence steps: [s|O|s] -> [a|O|s] ***

               rap(:) = f1*rab(:)

!         *** [p|O|s] = (Pi - Ai)*[s|O|s] +[s|dO|s]  (i = x,y,z) ***

               sc(2, 1) = rap(1)*sc(1, 1)-f2*kvec(1)*ss(1, 1)
               sc(3, 1) = rap(2)*sc(1, 1)-f2*kvec(2)*ss(1, 1)
               sc(4, 1) = rap(3)*sc(1, 1)-f2*kvec(3)*ss(1, 1)
               ss(2, 1) = rap(1)*ss(1, 1)+f2*kvec(1)*sc(1, 1)
               ss(3, 1) = rap(2)*ss(1, 1)+f2*kvec(2)*sc(1, 1)
               ss(4, 1) = rap(3)*ss(1, 1)+f2*kvec(3)*sc(1, 1)

               IF (PRESENT(dksc)) THEN
                  dsc(2,1,1) = rap(1)*dsc(1,1,1) - f2*ss(1,1) - f2*kvec(1)*dss(1,1,1) 
                  dsc(2,1,2) = rap(1)*dsc(1,1,2) - f2*kvec(1)*dss(1,1,2)
                  dsc(2,1,3) = rap(1)*dsc(1,1,3) - f2*kvec(1)*dss(1,1,3)

                  dsc(3,1,1) = rap(2)*dsc(1,1,1) - f2*kvec(2)*dss(1,1,1) 
                  dsc(3,1,2) = rap(2)*dsc(1,1,2) - f2*ss(1,1) - f2*kvec(2)*dss(1,1,2)
                  dsc(3,1,3) = rap(2)*dsc(1,1,3) - f2*kvec(2)*dss(1,1,3)

                  dsc(4,1,1) = rap(3)*dsc(1,1,1) - f2*kvec(3)*dss(1,1,1) 
                  dsc(4,1,2) = rap(3)*dsc(1,1,2) - f2*kvec(3)*dss(1,1,2)
                  dsc(4,1,3) = rap(3)*dsc(1,1,3) - f2*ss(1,1) - f2*kvec(3)*dss(1,1,3)
               END IF

               IF (PRESENT(dkss)) THEN
                  dss(2,1,1) = rap(1)*dss(1,1,1) + f2*sc(1,1) + f2*kvec(1)*dsc(1,1,1) 
                  dss(2,1,2) = rap(1)*dss(1,1,2) + f2*kvec(1)*dsc(1,1,2)
                  dss(2,1,3) = rap(1)*dss(1,1,3) + f2*kvec(1)*dsc(1,1,3)

                  dss(3,1,1) = rap(2)*dss(1,1,1) + f2*kvec(2)*dsc(1,1,1) 
                  dss(3,1,2) = rap(2)*dss(1,1,2) + f2*sc(1,1) + f2*kvec(2)*dsc(1,1,2)
                  dss(3,1,3) = rap(2)*dss(1,1,3) + f2*kvec(2)*dsc(1,1,3)

                  dss(4,1,1) = rap(3)*dss(1,1,1) + f2*kvec(3)*dsc(1,1,1) 
                  dss(4,1,2) = rap(3)*dss(1,1,2) + f2*kvec(3)*dsc(1,1,2)
                  dss(4,1,3) = rap(3)*dss(1,1,3) + f2*sc(1,1) + f2*kvec(3)*dsc(1,1,3)
               END IF

               IF (PRESENT(ddksc)) THEN
                  ddsc(2,1,1) = rap(1)*ddsc(1,1,1) - 2*f2*dss(1,1,1) - f2*kvec(1)*ddss(1,1,1) 
                  ddsc(2,1,2) = rap(1)*ddsc(1,1,2) - f2*dss(1,1,2) - f2*kvec(1)*ddss(1,1,2) 
                  ddsc(2,1,3) = rap(1)*ddsc(1,1,3) - f2*dss(1,1,3) - f2*kvec(1)*ddss(1,1,3) 
                  ddsc(2,1,4) = rap(1)*ddsc(1,1,4) - f2*kvec(1)*ddss(1,1,4) 
                  ddsc(2,1,5) = rap(1)*ddsc(1,1,5) - f2*kvec(1)*ddss(1,1,5) 
                  ddsc(2,1,6) = rap(1)*ddsc(1,1,6) - f2*kvec(1)*ddss(1,1,6) 

                  ddsc(3,1,1) = rap(2)*ddsc(1,1,1) - f2*kvec(2)*ddss(1,1,1)
                  ddsc(3,1,2) = rap(2)*ddsc(1,1,2) - f2*dss(1,1,1) - f2*kvec(2)*ddss(1,1,2)
                  ddsc(3,1,3) = rap(2)*ddsc(1,1,3) - f2*kvec(2)*ddss(1,1,3)
                  ddsc(3,1,4) = rap(2)*ddsc(1,1,4) - 2*f2*dss(1,1,2) - f2*kvec(2)*ddss(1,1,4)
                  ddsc(3,1,5) = rap(2)*ddsc(1,1,5) - f2*dss(1,1,3) - f2*kvec(2)*ddss(1,1,5) 
                  ddsc(3,1,6) = rap(2)*ddsc(1,1,6) - f2*kvec(2)*ddss(1,1,6)

                  ddsc(4,1,1) = rap(3)*ddsc(1,1,1) - f2*kvec(3)*ddss(1,1,1)
                  ddsc(4,1,2) = rap(3)*ddsc(1,1,2) - f2*kvec(3)*ddss(1,1,2)
                  ddsc(4,1,3) = rap(3)*ddsc(1,1,3) - f2*dss(1,1,1) - f2*kvec(3)*dss(1,1,3)
                  ddsc(4,1,4) = rap(3)*ddsc(1,1,4) - f2*kvec(3)*ddss(1,1,4)
                  ddsc(4,1,5) = rap(3)*ddsc(1,1,5) - f2*dss(1,1,2) - f2*kvec(3)*ddss(1,1,5)
                  ddsc(4,1,6) = rap(3)*ddsc(1,1,6) - 2*f2*dss(1,1,3) - f2*kvec(3)*ddss(1,1,6)
               END IF

               IF (PRESENT(ddkss)) THEN
                  ddss(2,1,1) = rap(1)*ddss(1,1,1) + 2*f2*dsc(1,1,1) + f2*kvec(1)*ddsc(1,1,1) 
                  ddss(2,1,2) = rap(1)*ddss(1,1,2) + f2*dsc(1,1,2) + f2*kvec(1)*ddsc(1,1,2) 
                  ddss(2,1,3) = rap(1)*ddss(1,1,3) + f2*dsc(1,1,3) + f2*kvec(1)*ddsc(1,1,3) 
                  ddss(2,1,4) = rap(1)*ddss(1,1,4) + f2*kvec(1)*ddsc(1,1,4) 
                  ddss(2,1,5) = rap(1)*ddss(1,1,5) + f2*kvec(1)*ddsc(1,1,5) 
                  ddss(2,1,6) = rap(1)*ddss(1,1,6) + f2*kvec(1)*ddsc(1,1,6) 

                  ddss(3,1,1) = rap(2)*ddss(1,1,1) + f2*kvec(2)*ddsc(1,1,1)
                  ddss(3,1,2) = rap(2)*ddss(1,1,2) + f2*dsc(1,1,1) + f2*kvec(2)*ddsc(1,1,2)
                  ddss(3,1,3) = rap(2)*ddss(1,1,3) + f2*kvec(2)*ddsc(1,1,3)
                  ddss(3,1,4) = rap(2)*ddss(1,1,4) + 2*f2*dsc(1,1,2) + f2*kvec(2)*ddsc(1,1,4)
                  ddss(3,1,5) = rap(2)*ddss(1,1,5) + f2*dsc(1,1,3) + f2*kvec(2)*ddsc(1,1,5) 
                  ddss(3,1,6) = rap(2)*ddss(1,1,6) + f2*kvec(2)*ddsc(1,1,6)

                  ddss(4,1,1) = rap(3)*ddss(1,1,1) + f2*kvec(3)*ddsc(1,1,1)
                  ddss(4,1,2) = rap(3)*ddss(1,1,2) + f2*kvec(3)*ddsc(1,1,2)
                  ddss(4,1,3) = rap(3)*ddss(1,1,3) + f2*dsc(1,1,1) + f2*kvec(3)*ddsc(1,1,3)
                  ddss(4,1,4) = rap(3)*ddss(1,1,4) + f2*kvec(3)*ddsc(1,1,4)
                  ddss(4,1,5) = rap(3)*ddss(1,1,5) + f2*dsc(1,1,2) + f2*kvec(3)*ddsc(1,1,5)
                  ddss(4,1,6) = rap(3)*ddss(1,1,6) + 2*f2*dsc(1,1,3) + f2*kvec(3)*ddsc(1,1,6)
               END IF


!         *** [a|O|s] = (Pi - Ai)*[a-1i|O|s] + f2*Ni(a-1i)*[a-2i|s] ***
!         ***           + [a-1i|dO|s]                               ***

               DO la = 2, la_max

!           *** Increase the angular momentum component z of function a ***

                  sc(coset(0, 0, la), 1) = rap(3)*sc(coset(0, 0, la-1), 1)+ &
                                           f2*REAL(la-1, dp)*sc(coset(0, 0, la-2), 1)- &
                                           f2*kvec(3)*ss(coset(0, 0, la-1), 1)
                  ss(coset(0, 0, la), 1) = rap(3)*ss(coset(0, 0, la-1), 1)+ &
                                           f2*REAL(la-1, dp)*ss(coset(0, 0, la-2), 1)+ &
                                           f2*kvec(3)*sc(coset(0, 0, la-1), 1)
           
                  IF (PRESENT(dksc)) THEN
                     dsc(coset(0, 0, la), 1, 1) = rap(3)*dsc(coset(0, 0, la-1), 1, 1) + &
                                                  f2*REAL(la-1, dp)*dsc(coset(0, 0, la-2), 1, 1)- &
                                                  f2*kvec(3)*dss(coset(0, 0, la-1), 1, 1)
                     dsc(coset(0, 0, la), 1, 2) = rap(3)*dsc(coset(0, 0, la-1), 1, 2)+ &
                                                  f2*REAL(la-1, dp)*dsc(coset(0, 0, la-2), 1, 2)- &
                                                  f2*kvec(3)*dss(coset(0, 0, la-1), 1, 2)
                     dsc(coset(0, 0, la), 1, 3) = rap(3)*dsc(coset(0, 0, la-1), 1, 3)+ &
                                                  f2*REAL(la-1, dp)*dsc(coset(0, 0, la-2), 1, 3)- &
                                                  f2*ss(coset(0, 0, la-1), 1)- &
                                                  f2*kvec(3)*dss(coset(0, 0, la-1), 1, 3)
                  END IF
                  IF  (PRESENT(ddksc)) THEN
                     ddsc(coset(0, 0, la), 1, 1) = rap(3)*ddsc(coset(0, 0, la-1), 1, 1) + &
                                                   f2*REAL(la-1, dp)*ddsc(coset(0, 0, la-2), 1, 1)- & 
                                                   f2*kvec(3)*ddss(coset(0, 0, la-1), 1, 1)
                     ddsc(coset(0, 0, la), 1, 2) = rap(3)*ddsc(coset(0, 0, la-1), 1, 2) + &
                                                   f2*REAL(la-1, dp)*ddsc(coset(0, 0, la-2), 1, 2)- &
                                                   f2*kvec(3)*ddss(coset(0, 0, la-1), 1, 2)
                     ddsc(coset(0, 0, la), 1, 3) = rap(3)*ddsc(coset(0, 0, la-1), 1, 3) + &
                                                   f2*REAL(la-1, dp)*ddsc(coset(0, 0, la-2), 1, 3)- &
                                                   f2*dss(coset(0, 0, la-1), 1, 1)- &
                                                   f2*kvec(3)*ddss(coset(0, 0, la-1), 1, 3)
                     ddsc(coset(0, 0, la), 1, 4) = rap(3)*ddsc(coset(0, 0, la-1), 1, 4) + &
                                                   f2*REAL(la-1, dp)*ddsc(coset(0, 0, la-2), 1, 4)- &
                                                   f2*kvec(3)*ddss(coset(0, 0, la-1), 1, 4)
                     ddsc(coset(0, 0, la), 1, 5) = rap(3)*ddsc(coset(0, 0, la-1), 1, 5) + &
                                                   f2*REAL(la-1, dp)*ddsc(coset(0, 0, la-2), 1, 5)- &
                                                   f2*dss(coset(0, 0, la-1), 1, 2)- &
                                                   f2*kvec(3)*ddss(coset(0, 0, la-1), 1, 5)
                     ddsc(coset(0, 0, la), 1, 6) = rap(3)*ddsc(coset(0, 0, la-1), 1, 6) + &
                                                   f2*REAL(la-1, dp)*ddsc(coset(0, 0, la-2), 1, 6)- &
                                                   2*f2*dss(coset(0, 0, la-1), 1, 3)- &
                                                   f2*kvec(3)*ddss(coset(0, 0, la-1), 1, 6)
                  END IF
                  IF (PRESENT(dkss)) THEN
                     dss(coset(0, 0, la), 1, 1) = rap(3)*dss(coset(0, 0, la-1), 1, 1) + &
                                                  f2*REAL(la-1, dp)*dss(coset(0, 0, la-2), 1, 1)+ &
                                                  f2*kvec(3)*dsc(coset(0, 0, la-1), 1, 1)
                     dss(coset(0, 0, la), 1, 2) = rap(3)*dss(coset(0, 0, la-1), 1, 2)+ &
                                                  f2*REAL(la-1, dp)*dss(coset(0, 0, la-2), 1, 2)+ &
                                                  f2*kvec(3)*dsc(coset(0, 0, la-1), 1, 2)
                     dss(coset(0, 0, la), 1, 3) = rap(3)*dss(coset(0, 0, la-1), 1, 3)+ &
                                                  f2*REAL(la-1, dp)*dss(coset(0, 0, la-2), 1, 3)+ &
                                                  f2*sc(coset(0, 0, la-1), 1)+ &
                                                  f2*kvec(3)*dsc(coset(0, 0, la-1), 1, 3)
                  END IF
                  IF  (PRESENT(ddkss)) THEN
                     ddss(coset(0, 0, la), 1, 1) = rap(3)*ddss(coset(0, 0, la-1), 1, 1) + &
                                                   f2*REAL(la-1, dp)*ddss(coset(0, 0, la-2), 1, 1)+ & 
                                                   f2*kvec(3)*ddsc(coset(0, 0, la-1), 1, 1)
                     ddss(coset(0, 0, la), 1, 2) = rap(3)*ddss(coset(0, 0, la-1), 1, 2) + &
                                                   f2*REAL(la-1, dp)*ddss(coset(0, 0, la-2), 1, 2)+ &
                                                   f2*kvec(3)*ddsc(coset(0, 0, la-1), 1, 2)
                     ddss(coset(0, 0, la), 1, 3) = rap(3)*ddss(coset(0, 0, la-1), 1, 3) + &
                                                   f2*REAL(la-1, dp)*ddss(coset(0, 0, la-2), 1, 3)+ &
                                                   f2*dsc(coset(0, 0, la-1), 1, 1)+ &
                                                   f2*kvec(3)*ddsc(coset(0, 0, la-1), 1, 3)
                     ddss(coset(0, 0, la), 1, 4) = rap(3)*ddss(coset(0, 0, la-1), 1, 4) + &
                                                   f2*REAL(la-1, dp)*ddss(coset(0, 0, la-2), 1, 4)+ &
                                                   f2*kvec(3)*ddsc(coset(0, 0, la-1), 1, 4)
                     ddss(coset(0, 0, la), 1, 5) = rap(3)*ddss(coset(0, 0, la-1), 1, 5) + &
                                                   f2*REAL(la-1, dp)*ddss(coset(0, 0, la-2), 1, 5)+ &
                                                   f2*dsc(coset(0, 0, la-1), 1, 2)+ &
                                                   f2*kvec(3)*ddsc(coset(0, 0, la-1), 1, 5)
                     ddss(coset(0, 0, la), 1, 6) = rap(3)*ddss(coset(0, 0, la-1), 1, 6) + &
                                                   f2*REAL(la-1, dp)*ddss(coset(0, 0, la-2), 1, 6)+ &
                                                   2*f2*dsc(coset(0, 0, la-1), 1, 3)+ &
                                                   f2*kvec(3)*ddsc(coset(0, 0, la-1), 1, 6)
                  END IF

!           *** Increase the angular momentum component y of function a ***

                  az = la-1
                  sc(coset(0, 1, az), 1) = rap(2)*sc(coset(0, 0, az), 1)- &
                                           f2*kvec(2)*ss(coset(0, 0, az), 1)
                  ss(coset(0, 1, az), 1) = rap(2)*ss(coset(0, 0, az), 1)+ &
                                           f2*kvec(2)*sc(coset(0, 0, az), 1)

                  IF (PRESENT(dksc)) THEN
                     dsc(coset(0, 1, az), 1, 1) = rap(2)*dsc(coset(0, 0, az), 1, 1)- &
                                                  f2*kvec(2)*dss(coset(0, 0, az), 1, 1) 
                     dsc(coset(0, 1, az), 1, 2) = rap(2)*dsc(coset(0, 0, az), 1, 2)- &
                                                  f2*ss(coset(0, 0, az), 1)- &
                                                  f2*kvec(2)*dss(coset(0, 0, az), 1, 2) 
                     dsc(coset(0, 1, az), 1, 3) = rap(2)*dsc(coset(0, 0, az), 1, 3)- &
                                                  f2*kvec(2)*dss(coset(0, 0, az), 1, 3) 
                  END IF
                  IF (PRESENT(ddksc)) THEN
                     ddsc(coset(0, 1, az), 1, 1) = rap(2)*ddsc(coset(0, 0, az), 1, 1)- &
                                                   f2*kvec(2)*ddss(coset(0, 0, az), 1, 1) 
                     ddsc(coset(0, 1, az), 1, 2) = rap(2)*ddsc(coset(0, 0, az), 1, 2)- & 
                                                   f2*dss(coset(0, 0, az), 1, 1)- &
                                                   f2*kvec(2)*ddss(coset(0, 0, az), 1, 2)
                     ddsc(coset(0, 1, az), 1, 3) = rap(2)*ddsc(coset(0, 0, az), 1, 3)- &
                                                   f2*kvec(2)*ddss(coset(0, 0, az), 1, 3) 
                     ddsc(coset(0, 1, az), 1, 4) = rap(2)*ddsc(coset(0, 0, az), 1, 4)- &
                                                   2*f2*dss(coset(0, 0, az), 1, 2)- &
                                                   f2*kvec(2)*ddss(coset(0, 0, az), 1, 4) 
                     ddsc(coset(0, 1, az), 1, 5) = rap(2)*ddsc(coset(0, 0, az), 1, 5)- &
                                                   f2*dss(coset(0, 0, az), 1, 3)- &
                                                   f2*kvec(2)*ddss(coset(0, 0, az), 1, 5) 
                     ddsc(coset(0, 1, az), 1, 6) = rap(2)*ddsc(coset(0, 0, az), 1, 6)- &
                                                   f2*kvec(2)*ddss(coset(0, 0, az), 1, 6) 
                  END IF
                  IF (PRESENT(dkss)) THEN
                     dss(coset(0, 1, az), 1, 1) = rap(2)*dss(coset(0, 0, az), 1, 1)+ &
                                                  f2*kvec(2)*dsc(coset(0, 0, az), 1, 1) 
                     dss(coset(0, 1, az), 1, 2) = rap(2)*dss(coset(0, 0, az), 1, 2)+ &
                                                  f2*sc(coset(0, 0, az), 1)+ &
                                                  f2*kvec(2)*dsc(coset(0, 0, az), 1, 2) 
                     dss(coset(0, 1, az), 1, 3) = rap(2)*dss(coset(0, 0, az), 1, 3)+ &
                                                  f2*kvec(2)*dsc(coset(0, 0, az), 1, 3) 
                  END IF
                  IF (PRESENT(ddkss)) THEN
                     ddss(coset(0, 1, az), 1, 1) = rap(2)*ddss(coset(0, 0, az), 1, 1)+ &
                                                   f2*kvec(2)*ddsc(coset(0, 0, az), 1, 1) 
                     ddss(coset(0, 1, az), 1, 2) = rap(2)*ddss(coset(0, 0, az), 1, 2)+ & 
                                                   f2*dsc(coset(0, 0, az), 1, 1)+ &
                                                   f2*kvec(2)*ddsc(coset(0, 0, az), 1, 2)
                     ddss(coset(0, 1, az), 1, 3) = rap(2)*ddss(coset(0, 0, az), 1, 3)+ &
                                                   f2*kvec(2)*ddsc(coset(0, 0, az), 1, 3) 
                     ddss(coset(0, 1, az), 1, 4) = rap(2)*ddss(coset(0, 0, az), 1, 4)+ &
                                                   2*f2*dsc(coset(0, 0, az), 1, 2)+ &
                                                   f2*kvec(2)*ddsc(coset(0, 0, az), 1, 4) 
                     ddss(coset(0, 1, az), 1, 5) = rap(2)*ddss(coset(0, 0, az), 1, 5)+ &
                                                   f2*dsc(coset(0, 0, az), 1, 3)+ &
                                                   f2*kvec(2)*ddsc(coset(0, 0, az), 1, 5) 
                     ddss(coset(0, 1, az), 1, 6) = rap(2)*ddss(coset(0, 0, az), 1, 6)+ &
                                                   f2*kvec(2)*ddsc(coset(0, 0, az), 1, 6) 
                  END IF

                  DO ay = 2, la
                     az = la-ay
                     sc(coset(0, ay, az), 1) = rap(2)*sc(coset(0, ay-1, az), 1)+ &
                                               f2*REAL(ay-1, dp)*sc(coset(0, ay-2, az), 1)- &
                                               f2*kvec(2)*ss(coset(0, ay-1, az), 1)
                     ss(coset(0, ay, az), 1) = rap(2)*ss(coset(0, ay-1, az), 1)+ &
                                               f2*REAL(ay-1, dp)*ss(coset(0, ay-2, az), 1)+ &
                                               f2*kvec(2)*sc(coset(0, ay-1, az), 1)
        
                     IF (PRESENT(dksc)) THEN
                        dsc(coset(0, ay, az), 1, 1) = rap(2)*dsc(coset(0, ay-1, az), 1, 1)+ &
                                                      f2*REAL(ay-1, dp)*dsc(coset(0, ay-2, az), 1, 1)- &
                                                      f2*kvec(2)*dss(coset(0, ay-1, az), 1, 1) 
                        dsc(coset(0, ay, az), 1, 2) = rap(2)*dsc(coset(0, ay-1, az), 1, 2)+ &
                                                      f2*REAL(ay-1, dp)*dsc(coset(0, ay-2, az), 1, 2)- &
                                                      f2*ss(coset(0, ay-1, az), 1)- &
                                                      f2*kvec(2)*dss(coset(0, ay-1, az), 1, 2) 
                        dsc(coset(0, ay, az), 1, 3) = rap(2)*dsc(coset(0, ay-1, az), 1, 3)+ &
                                                      f2*REAL(ay-1, dp)*dsc(coset(0, ay-2, az), 1, 3)- &
                                                      f2*kvec(2)*dss(coset(0, ay-1, az), 1, 3) 
                     END IF
                     IF (PRESENT(dkss)) THEN
                        dss(coset(0, ay, az), 1, 1) = rap(2)*dss(coset(0, ay-1, az), 1, 1)+ &
                                                      f2*REAL(ay-1, dp)*dss(coset(0, ay-2, az), 1, 1)+ &
                                                      f2*kvec(2)*dsc(coset(0, ay-1, az), 1, 1) 
                        dss(coset(0, ay, az), 1, 2) = rap(2)*dss(coset(0, ay-1, az), 1, 2)+ &
                                                      f2*REAL(ay-1, dp)*dss(coset(0, ay-2, az), 1, 2)+ &
                                                      f2*sc(coset(0, ay-1, az), 1)+ &
                                                      f2*kvec(2)*dsc(coset(0, ay-1, az), 1, 2) 
                        dss(coset(0, ay, az), 1, 3) = rap(2)*dss(coset(0, ay-1, az), 1, 3)+ &
                                                      f2*REAL(ay-1, dp)*dss(coset(0, ay-2, az), 1, 3)+ &
                                                      f2*kvec(2)*dsc(coset(0, ay-1, az), 1, 3) 
                     END IF
                     IF (PRESENT(ddksc)) THEN
                        ddsc(coset(0, ay, az), 1, 1) = rap(2)*ddsc(coset(0, ay-1, az), 1, 1)+ &
                                                       f2*REAL(ay-1, dp)*ddsc(coset(0, ay-2, az), 1, 1)- &
                                                       f2*kvec(2)*ddss(coset(0, ay-1, az), 1, 1) 
                        ddsc(coset(0, ay, az), 1, 2) = rap(2)*ddsc(coset(0, ay-1, az), 1, 2)+ & 
                                                       f2*REAL(ay-1, dp)*ddsc(coset(0, ay-2, az), 1, 2)- &
                                                       f2*dss(coset(0, ay-1, az), 1, 1)- &
                                                       f2*kvec(2)*ddss(coset(0, ay-1, az), 1, 2)
                        ddsc(coset(0, ay, az), 1, 3) = rap(2)*ddsc(coset(0, ay-1, az), 1, 3)+ &
                                                       f2*REAL(ay-1, dp)*ddsc(coset(0, ay-2, az), 1, 3)- &
                                                       f2*kvec(2)*ddss(coset(0, ay-1, az), 1, 3) 
                        ddsc(coset(0, ay, az), 1, 4) = rap(2)*ddsc(coset(0, ay-1, az), 1, 4)+ &
                                                       f2*REAL(ay-1, dp)*ddsc(coset(0, ay-2, az), 1, 4)- &
                                                       2*f2*dss(coset(0, ay-1, az), 1, 2)- &
                                                       f2*kvec(2)*ddss(coset(0, 0, az), 1, 4) 
                        ddsc(coset(0, ay, az), 1, 5) = rap(2)*ddsc(coset(0, ay-1, az), 1, 5)+ &
                                                       f2*REAL(ay-1, dp)*ddsc(coset(0, ay-2, az), 1, 5)- &
                                                       f2*dss(coset(0, ay-1, az), 1, 3)- &
                                                       f2*kvec(2)*ddss(coset(0, ay-1, az), 1, 5) 
                        ddsc(coset(0, ay, az), 1, 6) = rap(2)*ddsc(coset(0, ay-1, az), 1, 6)+ &
                                                       f2*REAL(ay-1, dp)*ddsc(coset(0, ay-2, az), 1, 6)- &
                                                       f2*kvec(2)*ddss(coset(0, ay-1, az), 1, 6) 
                     END IF
                     IF (PRESENT(ddkss)) THEN
                        ddss(coset(0, ay, az), 1, 1) = rap(2)*ddss(coset(0, ay-1, az), 1, 1)+ &
                                                       f2*REAL(ay-1, dp)*ddss(coset(0, ay-2, az), 1, 1)+ &
                                                       f2*kvec(2)*ddsc(coset(0, ay-1, az), 1, 1) 
                        ddss(coset(0, ay, az), 1, 2) = rap(2)*ddss(coset(0, ay-1, az), 1, 2)+ & 
                                                       f2*REAL(ay-1, dp)*ddss(coset(0, ay-2, az), 1, 2)+ &
                                                       f2*dsc(coset(0, ay-1, az), 1, 1)+ &
                                                       f2*kvec(2)*ddsc(coset(0, ay-1, az), 1, 2)
                        ddss(coset(0, ay, az), 1, 3) = rap(2)*ddss(coset(0, ay-1, az), 1, 3)+ &
                                                       f2*REAL(ay-1, dp)*ddss(coset(0, ay-2, az), 1, 3)+ &
                                                       f2*kvec(2)*ddsc(coset(0, ay-1, az), 1, 3) 
                        ddss(coset(0, ay, az), 1, 4) = rap(2)*ddss(coset(0, ay-1, az), 1, 4)+ &
                                                       f2*REAL(ay-1, dp)*ddss(coset(0, ay-2, az), 1, 4)+ &
                                                       2*f2*dsc(coset(0, ay-1, az), 1, 2)+ &
                                                       f2*kvec(2)*ddsc(coset(0, ay-1, az), 1, 4) 
                        ddss(coset(0, ay, az), 1, 5) = rap(2)*ddss(coset(0, ay-1, az), 1, 5)+ &
                                                       f2*REAL(ay-1, dp)*ddss(coset(0, ay-2, az), 1, 5)+ &
                                                       f2*dsc(coset(0, ay-1, az), 1, 3)+ &
                                                       f2*kvec(2)*ddsc(coset(0, ay-1, az), 1, 5) 
                        ddss(coset(0, ay, az), 1, 6) = rap(2)*ddss(coset(0, ay-1, az), 1, 6)+ &
                                                       f2*REAL(ay-1, dp)*ddss(coset(0, ay-2, az), 1, 6)+ &
                                                       f2*kvec(2)*ddsc(coset(0, ay-1, az), 1, 6) 
                     END IF
                               
                  END DO

!           *** Increase the angular momentum component x of function a ***

                  DO ay = 0, la-1
                     az = la-1-ay
                     sc(coset(1, ay, az), 1) = rap(1)*sc(coset(0, ay, az), 1)- &
                                               f2*kvec(1)*ss(coset(0, ay, az), 1)
                     ss(coset(1, ay, az), 1) = rap(1)*ss(coset(0, ay, az), 1)+ &
                                               f2*kvec(1)*sc(coset(0, ay, az), 1)
                  
                     IF (PRESENT(dksc)) THEN
                        dsc(coset(1, ay, az), 1, 1) = rap(1)*dsc(coset(0, ay, az), 1, 1)- &
                                                      f2*ss(coset(0, ay, az), 1)- &
                                                      f2*kvec(1)*dss(coset(0, ay, az), 1, 1)

                        dsc(coset(1, ay, az), 1, 2) = rap(1)*dsc(coset(0, ay, az), 1, 2)- &
                                                      f2*kvec(1)*dss(coset(0, ay, az), 1, 2)

                        dsc(coset(1, ay, az), 1, 3) = rap(1)*dsc(coset(0, ay, az), 1, 3)- &
                                                      f2*kvec(1)*dss(coset(0, ay, az), 1, 3)
                     END IF

                     IF (PRESENT(ddksc)) THEN
                        ddsc(coset(1, ay, az), 1, 1) = rap(1)*ddsc(coset(0, ay, az), 1, 1)- &
                                                       2*f2*dss(coset(0, ay, az), 1, 1)- &
                                                       f2*kvec(1)*ddss(coset(0, ay, az), 1, 1)
                        ddsc(coset(1, ay, az), 1, 2) = rap(1)*ddsc(coset(0, ay, az), 1, 2)- &
                                                       f2*dss(coset(0, ay, az), 1, 2)- &
                                                       f2*kvec(1)*ddss(coset(0, ay, az), 1, 2)
                        ddsc(coset(1, ay, az), 1, 3) = rap(1)*ddsc(coset(0, ay, az), 1, 3)- &
                                                       f2*dss(coset(0, ay, az), 1, 3)- &
                                                       f2*kvec(1)*ddss(coset(0, ay, az), 1, 3)
                        ddsc(coset(1, ay, az), 1, 4) = rap(1)*ddsc(coset(0, ay, az), 1, 4)- &
                                                       f2*kvec(1)*ddss(coset(0, ay, az), 1, 4)
                        ddsc(coset(1, ay, az), 1, 5) = rap(1)*ddsc(coset(0, ay, az), 1, 5)- &
                                                       f2*kvec(1)*ddss(coset(0, ay, az), 1, 5)
                        ddsc(coset(1, ay, az), 1, 6) = rap(1)*ddsc(coset(0, ay, az), 1, 6)- &
                                                       f2*kvec(1)*ddss(coset(0, ay, az), 1, 6)
                     END IF

                     IF (PRESENT(dkss)) THEN
                        dss(coset(1, ay, az), 1, 1) = rap(1)*dss(coset(0, ay, az), 1, 1)+ &
                                                      f2*sc(coset(0, ay, az), 1)+ &
                                                      f2*kvec(1)*dsc(coset(0, ay, az), 1, 1)

                        dss(coset(1, ay, az), 1, 2) = rap(1)*dss(coset(0, ay, az), 1, 2)+ &
                                                      f2*kvec(1)*dsc(coset(0, ay, az), 1, 2)

                        dss(coset(1, ay, az), 1, 3) = rap(1)*dss(coset(0, ay, az), 1, 3)+ &
                                                      f2*kvec(1)*dsc(coset(0, ay, az), 1, 3)
                     END IF
                     
                     IF (PRESENT(ddkss)) THEN
                        ddss(coset(1, ay, az), 1, 1) = rap(1)*ddss(coset(0, ay, az), 1, 1)+ &
                                                       2*f2*dsc(coset(0, ay, az), 1, 1)+ &
                                                       f2*kvec(1)*ddsc(coset(0, ay, az), 1, 1)
                        ddss(coset(1, ay, az), 1, 2) = rap(1)*ddss(coset(0, ay, az), 1, 2)+ &
                                                       f2*dsc(coset(0, ay, az), 1, 2)+ &
                                                       f2*kvec(1)*ddsc(coset(0, ay, az), 1, 2)
                        ddss(coset(1, ay, az), 1, 3) = rap(1)*ddss(coset(0, ay, az), 1, 3)+ &
                                                       f2*dsc(coset(0, ay, az), 1, 3)+ &
                                                       f2*kvec(1)*ddsc(coset(0, ay, az), 1, 3)
                        ddss(coset(1, ay, az), 1, 4) = rap(1)*ddss(coset(0, ay, az), 1, 4)+ &
                                                       f2*kvec(1)*ddsc(coset(0, ay, az), 1, 4)
                        ddss(coset(1, ay, az), 1, 5) = rap(1)*ddss(coset(0, ay, az), 1, 5)+ &
                                                       f2*kvec(1)*ddsc(coset(0, ay, az), 1, 5)
                        ddss(coset(1, ay, az), 1, 6) = rap(1)*ddss(coset(0, ay, az), 1, 6)+ &
                                                       f2*kvec(1)*ddsc(coset(0, ay, az), 1, 6)
                     END IF
                                         
                  END DO

                  DO ax = 2, la
                     f3 = f2*REAL(ax-1, dp)
                     DO ay = 0, la-ax
                        az = la-ax-ay
                        sc(coset(ax, ay, az), 1) = rap(1)*sc(coset(ax-1, ay, az), 1)+ &
                                                   f3*sc(coset(ax-2, ay, az), 1)- &
                                                   f2*kvec(1)*ss(coset(ax-1, ay, az), 1)
                        ss(coset(ax, ay, az), 1) = rap(1)*ss(coset(ax-1, ay, az), 1)+ &
                                                   f3*ss(coset(ax-2, ay, az), 1)+ &
                                                   f2*kvec(1)*sc(coset(ax-1, ay, az), 1)

                     IF (PRESENT(dksc)) THEN
                        dsc(coset(ax, ay, az), 1, 1) = rap(1)*dsc(coset(ax-1, ay, az), 1, 1)+ &
                                                       f3*dsc(coset(ax-2, ay, az), 1, 1)- &
                                                       f2*ss(coset(ax-1, ay, az), 1)- &
                                                       f2*kvec(1)*dss(coset(ax-1, ay, az), 1, 1)
                        dsc(coset(ax, ay, az), 1, 2) = rap(1)*dsc(coset(ax-1, ay, az), 1, 2)+ &
                                                       f3*dsc(coset(ax-2, ay, az), 1, 2)- &
                                                       f2*kvec(1)*dss(coset(ax-1, ay, az), 1, 2)
                        dsc(coset(ax, ay, az), 1, 3) = rap(1)*dsc(coset(ax-1, ay, az), 1, 3)+ &
                                                       f3*dsc(coset(ax-2, ay, az), 1, 3)- &
                                                       f2*kvec(1)*dss(coset(ax-1, ay, az), 1, 3)
                     END IF
                     IF (PRESENT(dkss)) THEN
                        dss(coset(ax, ay, az), 1, 1) = rap(1)*dss(coset(ax-1, ay, az), 1, 1)+ &
                                                       f3*dss(coset(ax-2, ay, az), 1, 1)+ &
                                                       f2*sc(coset(ax-1, ay, az), 1)+ &
                                                       f2*kvec(1)*dsc(coset(ax-1, ay, az), 1, 1)
                        dss(coset(ax, ay, az), 1, 2) = rap(1)*dss(coset(ax-1, ay, az), 1, 2)+ &
                                                       f3*dss(coset(ax-2, ay, az), 1, 2)+ &
                                                       f2*kvec(1)*dsc(coset(ax-1, ay, az), 1, 2)
                        dss(coset(ax, ay, az), 1, 3) = rap(1)*dss(coset(ax-1, ay, az), 1, 3)+ &
                                                       f3*dss(coset(ax-2, ay, az), 1, 3)+ &
                                                       f2*kvec(1)*dsc(coset(ax-1, ay, az), 1, 3)
                     END IF
                     IF (PRESENT(ddksc)) THEN
                        ddsc(coset(ax, ay, az), 1, 1) = rap(1)*ddsc(coset(ax-1, ay, az), 1, 1)+ &
                                                        f3*ddsc(coset(ax-2, ay, az), 1, 1)- &
                                                        2*f2*dss(coset(ax-1, ay, az), 1, 1)- &
                                                        f2*kvec(1)*ddss(coset(ax-1, ay, az), 1, 1)
                        ddsc(coset(ax, ay, az), 1, 2) = rap(1)*ddsc(coset(ax-1, ay, az), 1, 2)+ &
                                                        f3*ddsc(coset(ax-2, ay, az), 1, 2)- &
                                                        f2*dss(coset(ax-1, ay, az), 1, 2)- &
                                                        f2*kvec(1)*ddss(coset(ax-1, ay, az), 1, 2)
                        ddsc(coset(ax, ay, az), 1, 3) = rap(1)*ddsc(coset(ax-1, ay, az), 1, 3)+ &
                                                        f3*ddsc(coset(ax-2, ay, az), 1, 3)- &
                                                        f2*dss(coset(ax-1, ay, az), 1, 3)- &
                                                        f2*kvec(1)*ddss(coset(ax-1, ay, az), 1, 3)
                        ddsc(coset(ax, ay, az), 1, 4) = rap(1)*ddsc(coset(ax-1, ay, az), 1, 4)+ &
                                                        f3*ddsc(coset(ax-2, ay, az), 1, 4)- &
                                                        f2*kvec(1)*ddss(coset(ax-1, ay, az), 1, 4)
                        ddsc(coset(ax, ay, az), 1, 5) = rap(1)*ddsc(coset(ax-1, ay, az), 1, 5)+ &
                                                        f3*ddsc(coset(ax-2, ay, az), 1, 5)- &
                                                        f2*kvec(1)*ddss(coset(ax-1, ay, az), 1, 5)
                        ddsc(coset(ax, ay, az), 1, 6) = rap(1)*ddsc(coset(ax-1, ay, az), 1, 6)+ &
                                                        f3*ddsc(coset(ax-2, ay, az), 1, 6)- &
                                                        f2*kvec(1)*ddss(coset(ax-1, ay, az), 1, 6)
                     END IF

                     IF (PRESENT(ddkss)) THEN
                        ddss(coset(ax, ay, az), 1, 1) = rap(1)*ddss(coset(ax-1, ay, az), 1, 1)+ &
                                                        f3*ddss(coset(ax-2, ay, az), 1, 1)+ &
                                                        2*f2*dsc(coset(ax-1, ay, az), 1, 1)+ &
                                                        f2*kvec(1)*ddsc(coset(ax-1, ay, az), 1, 1)
                        ddss(coset(ax, ay, az), 1, 2) = rap(1)*ddss(coset(ax-1, ay, az), 1, 2)+ &
                                                        f3*ddss(coset(ax-2, ay, az), 1, 2)+ &
                                                        f2*dsc(coset(ax-1, ay, az), 1, 2)+ &
                                                        f2*kvec(1)*ddsc(coset(ax-1, ay, az), 1, 2)
                        ddss(coset(ax, ay, az), 1, 3) = rap(1)*ddss(coset(ax-1, ay, az), 1, 3)+ &
                                                        f3*ddss(coset(ax-2, ay, az), 1, 3)+ &
                                                        f2*dsc(coset(ax-1, ay, az), 1, 3)+ &
                                                        f2*kvec(1)*ddsc(coset(ax-1, ay, az), 1, 3)
                        ddss(coset(ax, ay, az), 1, 4) = rap(1)*ddss(coset(ax-1, ay, az), 1, 4)+ &
                                                        f3*ddss(coset(ax-2, ay, az), 1, 4)+ &
                                                        f2*kvec(1)*ddsc(coset(ax-1, ay, az), 1, 4)
                        ddss(coset(ax, ay, az), 1, 5) = rap(1)*ddss(coset(ax-1, ay, az), 1, 5)+ &
                                                        f3*ddss(coset(ax-2, ay, az), 1, 5)+ &
                                                        f2*kvec(1)*ddsc(coset(ax-1, ay, az), 1, 5)
                        ddss(coset(ax, ay, az), 1, 6) = rap(1)*ddss(coset(ax-1, ay, az), 1, 6)+ &
                                                        f3*ddss(coset(ax-2, ay, az), 1, 6)+ &
                                                        f2*kvec(1)*ddsc(coset(ax-1, ay, az), 1, 6)
                     END IF


                     END DO
                  END DO

               END DO

!         *** Recurrence steps: [a|O|s] -> [a|O|b] ***

               IF (lb_max > 0) THEN

                  DO j = 2, ncoset(lb_max)
                     DO i = 1, ncoset(la_max)
                        sc(i, j) = 0.0_dp
                        ss(i, j) = 0.0_dp

                        IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                           DO k = 1, 3
                              dsc(i, j, k) = 0.0_dp
                              dss(i, j, k) = 0.0_dp
                           END DO
                        END IF
                        IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                           DO k = 1, 6
                              ddsc(i, j, k) = 0.0_dp
                              ddss(i, j, k) = 0.0_dp
                           END DO
                        END IF

                     END DO
                  END DO

!           *** Horizontal recurrence steps ***

                  rbp(:) = rap(:)-rab(:)

!           *** [a|O|p] = [a+1i|O|s] - (Bi - Ai)*[a|O|s] ***

                  IF (lb_max == 1) THEN
                     la_start = la_min
                  ELSE
                     la_start = MAX(0, la_min-1)
                  END IF

                  DO la = la_start, la_max-1
                     DO ax = 0, la
                        DO ay = 0, la-ax
                           az = la-ax-ay
                           sc(coset(ax, ay, az), 2) = sc(coset(ax+1, ay, az), 1)- &
                                                      rab(1)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 3) = sc(coset(ax, ay+1, az), 1)- &
                                                      rab(2)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 4) = sc(coset(ax, ay, az+1), 1)- &
                                                      rab(3)*sc(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = ss(coset(ax+1, ay, az), 1)- &
                                                      rab(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = ss(coset(ax, ay+1, az), 1)- &
                                                      rab(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = ss(coset(ax, ay, az+1), 1)- &
                                                      rab(3)*ss(coset(ax, ay, az), 1)

                        IF (PRESENT(dksc)) THEN
                           DO i = 1, 3
                              dsc(coset(ax, ay, az), 2, i) = dsc(coset(ax+1, ay, az), 1, i)- &
                                                             rab(1)*dsc(coset(ax, ay, az), 1, i)
                              dsc(coset(ax, ay, az), 3, i) = dsc(coset(ax, ay+1, az), 1, i)- &
                                                             rab(2)*dsc(coset(ax, ay, az), 1, i)
                              dsc(coset(ax, ay, az), 4, i) = dsc(coset(ax, ay, az+1), 1, i)- &
                                                             rab(3)*dsc(coset(ax, ay, az), 1, i)
                           END DO
                        END IF
                        
                        IF (PRESENT(dkss)) THEN
                           DO i = 1, 3
                              dss(coset(ax, ay, az), 2, i) = dss(coset(ax+1, ay, az), 1, i)- &
                                                             rab(1)*dss(coset(ax, ay, az), 1, i)
                              dss(coset(ax, ay, az), 3, i) = dss(coset(ax, ay+1, az), 1, i)- &
                                                             rab(2)*dss(coset(ax, ay, az), 1, i)
                              dss(coset(ax, ay, az), 4, i) = dss(coset(ax, ay, az+1), 1, i)- &
                                                             rab(3)*dss(coset(ax, ay, az), 1, i)
                           END DO
                        END IF
                        IF (PRESENT(ddksc)) THEN
                           DO i = 1, 6
                              ddsc(coset(ax, ay, az), 2, i) = ddsc(coset(ax+1, ay, az), 1, i)- &
                                                              rab(1)*ddsc(coset(ax, ay, az), 1, i)
                              ddsc(coset(ax, ay, az), 3, i) = ddsc(coset(ax, ay+1, az), 1, i)- &
                                                              rab(2)*ddsc(coset(ax, ay, az), 1, i)
                              ddsc(coset(ax, ay, az), 4, i) = ddsc(coset(ax, ay, az+1), 1, i)- &
                                                              rab(3)*ddsc(coset(ax, ay, az), 1, i)
                           END DO
                        END IF
                        
                        IF (PRESENT(ddkss)) THEN
                           DO i = 1, 6
                              ddss(coset(ax, ay, az), 2, i) = ddss(coset(ax+1, ay, az), 1, i)- &
                                                              rab(1)*ddss(coset(ax, ay, az), 1, i)
                              ddss(coset(ax, ay, az), 3, i) = ddss(coset(ax, ay+1, az), 1, i)- &
                                                              rab(2)*ddss(coset(ax, ay, az), 1, i)
                              ddss(coset(ax, ay, az), 4, i) = ddss(coset(ax, ay, az+1), 1, i)- &
                                                              rab(3)*ddss(coset(ax, ay, az), 1, i)
                           END DO
                        END IF

                        END DO
                     END DO
                  END DO

!           *** Vertical recurrence step ***

!           *** [a|O|p] = (Pi - Bi)*[a|O|s] + f2*Ni(a)*[a-1i|O|s] ***
!           ***           + [a|dO|s]                              ***

                  DO ax = 0, la_max
                     fx = f2*REAL(ax, dp)
                     DO ay = 0, la_max-ax
                        fy = f2*REAL(ay, dp)
                        az = la_max-ax-ay
                        fz = f2*REAL(az, dp)
                        IF (ax == 0) THEN
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------

                           IF (PRESENT(dksc)) THEN
                              dsc(coset(ax, ay, az), 2, 1) = rbp(1)*dsc(coset(ax, ay, az), 1, 1)- &
                                                             f2*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 2, 2) = rbp(1)*dsc(coset(ax, ay, az), 1, 2)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 2, 3) = rbp(1)*dsc(coset(ax, ay, az),1 , 3)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddksc)) THEN
                              ddsc(coset(ax, ay, az), 2, 1) = rbp(1)*ddsc(coset(ax, ay, az), 1, 1)- &
                                                              2*f2*dss(coset(ax, ay, az), 1, 1)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 1)
                              ddsc(coset(ax, ay, az), 2, 2) = rbp(1)*ddsc(coset(ax, ay, az), 1, 2)- &
                                                              f2*dss(coset(ax, ay, az), 1, 2)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 2)
                              ddsc(coset(ax, ay, az), 2, 3) = rbp(1)*ddsc(coset(ax, ay, az), 1, 3)- &
                                                              f2*dss(coset(ax, ay, az), 1, 3)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 3)
                              ddsc(coset(ax, ay, az), 2, 4) = rbp(1)*ddsc(coset(ax, ay, az), 1, 4)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 4)
                              ddsc(coset(ax, ay, az), 2, 5) = rbp(1)*ddsc(coset(ax, ay, az), 1, 5)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 5)
                              ddsc(coset(ax, ay, az), 2, 6) = rbp(1)*ddsc(coset(ax, ay, az), 1, 6)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 6)
                           END IF
                           IF (PRESENT(dkss)) THEN
                              dss(coset(ax, ay, az), 2, 1) = rbp(1)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             f2*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 2, 2) = rbp(1)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 2, 3) = rbp(1)*dss(coset(ax, ay, az),1 , 3)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddkss)) THEN
                              ddss(coset(ax, ay, az), 2, 1) = rbp(1)*ddss(coset(ax, ay, az), 1, 1)+ &
                                                              2*f2*dsc(coset(ax, ay, az), 1, 1)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 1)
                              ddss(coset(ax, ay, az), 2, 2) = rbp(1)*ddss(coset(ax, ay, az), 1, 2)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 2)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 2)
                              ddss(coset(ax, ay, az), 2, 3) = rbp(1)*ddss(coset(ax, ay, az), 1, 3)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 3)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 3)
                              ddss(coset(ax, ay, az), 2, 4) = rbp(1)*ddss(coset(ax, ay, az), 1, 4)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 4)
                              ddss(coset(ax, ay, az), 2, 5) = rbp(1)*ddss(coset(ax, ay, az), 1, 5)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 5)
                              ddss(coset(ax, ay, az), 2, 6) = rbp(1)*ddss(coset(ax, ay, az), 1, 6)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 6)
                           END IF

                        ELSE
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1)+ &
                                                      fx*sc(coset(ax-1, ay, az), 1)- &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1)+ &
                                                      fx*ss(coset(ax-1, ay, az), 1)+ &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------

                           IF (PRESENT(dksc)) THEN
                              dsc(coset(ax, ay, az), 2, 1) = rbp(1)*dsc(coset(ax, ay, az), 1, 1)+ &
                                                             fx*dsc(coset(ax-1, ay, az), 1, 1)- &
                                                             f2*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 2, 2) = rbp(1)*dsc(coset(ax, ay, az), 1, 2)+ &
                                                             fx*dsc(coset(ax-1, ay, az), 1, 2)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 2, 3) = rbp(1)*dsc(coset(ax, ay, az), 1, 3)+ &
                                                             fx*dsc(coset(ax-1, ay, az), 1, 3)- &
                                                             f2*kvec(1)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddksc)) THEN
                              ddsc(coset(ax, ay, az), 2, 1) = rbp(1)*ddsc(coset(ax, ay, az), 1, 1)+ &
                                                              fx*ddsc(coset(ax-1, ay, az), 1, 1)- &
                                                              2*f2*dss(coset(ax, ay, az), 1, 1)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 1)
                              ddsc(coset(ax, ay, az), 2, 2) = rbp(1)*ddsc(coset(ax, ay, az), 1, 2)+ &
                                                              fx*ddsc(coset(ax-1, ay, az), 1, 2)- &
                                                              f2*dss(coset(ax, ay, az), 1, 2)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 2)
                              ddsc(coset(ax, ay, az), 2, 3) = rbp(1)*ddsc(coset(ax, ay, az), 1, 3)+ &
                                                              fx*ddsc(coset(ax-1, ay, az), 1, 3)- &
                                                              f2*dss(coset(ax, ay, az), 1, 3)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 3)
                              ddsc(coset(ax, ay, az), 2, 4) = rbp(1)*ddsc(coset(ax, ay, az), 1, 4)+ &
                                                              fx*ddsc(coset(ax-1, ay, az), 1, 4)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 4)
                              ddsc(coset(ax, ay, az), 2, 5) = rbp(1)*ddsc(coset(ax, ay, az), 1, 5)+ &
                                                              fx*ddsc(coset(ax-1, ay, az), 1, 5)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 5)
                              ddsc(coset(ax, ay, az), 2, 6) = rbp(1)*ddsc(coset(ax, ay, az), 1, 6)+ &
                                                              fx*ddsc(coset(ax-1, ay, az), 1, 6)- &
                                                              f2*kvec(1)*ddss(coset(ax, ay, az), 1, 6)
                           END IF
                           IF (PRESENT(dkss)) THEN
                              dss(coset(ax, ay, az), 2, 1) = rbp(1)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             fx*dss(coset(ax-1, ay, az), 1, 1)+ &
                                                             f2*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 2, 2) = rbp(1)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             fx*dss(coset(ax-1, ay, az), 1, 2)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 2, 3) = rbp(1)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             fx*dss(coset(ax-1, ay, az), 1, 3)+ &
                                                             f2*kvec(1)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddkss)) THEN
                              ddss(coset(ax, ay, az), 2, 1) = rbp(1)*ddss(coset(ax, ay, az), 1, 1)+ &
                                                              fx*ddss(coset(ax-1, ay, az), 1, 1)+ &
                                                              2*f2*dsc(coset(ax, ay, az), 1, 1)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 1)
                              ddss(coset(ax, ay, az), 2, 2) = rbp(1)*ddss(coset(ax, ay, az), 1, 2)+ &
                                                              fx*ddss(coset(ax-1, ay, az), 1, 2)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 2)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 2)
                              ddss(coset(ax, ay, az), 2, 3) = rbp(1)*ddss(coset(ax, ay, az), 1, 3)+ &
                                                              fx*ddss(coset(ax-1, ay, az), 1, 3)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 3)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 3)
                              ddss(coset(ax, ay, az), 2, 4) = rbp(1)*ddss(coset(ax, ay, az), 1, 4)+ &
                                                              fx*ddss(coset(ax-1, ay, az), 1, 4)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 4)
                              ddss(coset(ax, ay, az), 2, 5) = rbp(1)*ddss(coset(ax, ay, az), 1, 5)+ &
                                                              fx*ddss(coset(ax-1, ay, az), 1, 5)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 5)
                              ddss(coset(ax, ay, az), 2, 6) = rbp(1)*ddss(coset(ax, ay, az), 1, 6)+ &
                                                              fx*ddss(coset(ax-1, ay, az), 1, 6)+ &
                                                              f2*kvec(1)*ddsc(coset(ax, ay, az), 1, 6)
                            END IF
                        END IF
                        IF (ay == 0) THEN
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------
                        
                           IF (PRESENT(dksc)) THEN
                              dsc(coset(ax, ay, az), 3, 1) = rbp(2)*dsc(coset(ax, ay, az), 1, 1)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 3, 2) = rbp(2)*dsc(coset(ax, ay, az), 1, 2)- &
                                                             f2*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 3, 3) = rbp(2)*dsc(coset(ax, ay, az), 1, 3)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddksc)) THEN
                              ddsc(coset(ax, ay, az), 3, 1) = rbp(2)*ddsc(coset(ax, ay, az), 1, 1)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 1)
                              ddsc(coset(ax, ay, az), 3, 2) = rbp(2)*ddsc(coset(ax, ay, az), 1, 2)- &
                                                              f2*dss(coset(ax, ay, az), 1, 1)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 2)
                              ddsc(coset(ax, ay, az), 3, 3) = rbp(2)*ddsc(coset(ax, ay, az), 1, 3)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 3)
                              ddsc(coset(ax, ay, az), 3, 4) = rbp(2)*ddsc(coset(ax, ay, az), 1, 4)- &
                                                              2*f2*dss(coset(ax, ay, az), 1, 2)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 4)
                              ddsc(coset(ax, ay, az), 3, 5) = rbp(2)*ddsc(coset(ax, ay, az), 1, 5)- &
                                                              f2*dss(coset(ax, ay, az), 1, 3)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 5)
                              ddsc(coset(ax, ay, az), 3, 6) = rbp(2)*ddsc(coset(ax, ay, az), 1, 6)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 6)
                           END IF
                           IF (PRESENT(dkss)) THEN
                              dss(coset(ax, ay, az), 3, 1) = rbp(2)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 3, 2) = rbp(2)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             f2*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 3, 3) = rbp(2)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddkss)) THEN
                              ddss(coset(ax, ay, az), 3, 1) = rbp(2)*ddss(coset(ax, ay, az), 1, 1)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 1)
                              ddss(coset(ax, ay, az), 3, 2) = rbp(2)*ddss(coset(ax, ay, az), 1, 2)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 1)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 2)
                              ddss(coset(ax, ay, az), 3, 3) = rbp(2)*ddss(coset(ax, ay, az), 1, 3)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 3)
                              ddss(coset(ax, ay, az), 3, 4) = rbp(2)*ddss(coset(ax, ay, az), 1, 4)+ &
                                                              2*f2*dsc(coset(ax, ay, az), 1, 2)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 4)
                              ddss(coset(ax, ay, az), 3, 5) = rbp(2)*ddss(coset(ax, ay, az), 1, 5)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 3)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 5)
                              ddss(coset(ax, ay, az), 3, 6) = rbp(2)*ddss(coset(ax, ay, az), 1, 6)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 6)
                           END IF
                                                
                        ELSE
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1)+ &
                                                      fy*sc(coset(ax, ay-1, az), 1)- &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1)+ &
                                                      fy*ss(coset(ax, ay-1, az), 1)+ &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------
                           
                           IF (PRESENT(dksc)) THEN
                              dsc(coset(ax, ay, az), 3, 1) = rbp(2)*dsc(coset(ax, ay, az), 1, 1)+ &
                                                             fy*dsc(coset(ax, ay-1, az), 1, 1)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 3, 2) = rbp(2)*dsc(coset(ax, ay, az), 1, 2)+ &
                                                             fy*dsc(coset(ax, ay-1, az), 1, 2)- &
                                                             f2*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 3, 3) = rbp(2)*dsc(coset(ax, ay, az), 1, 3)+ &
                                                             fy*dsc(coset(ax, ay-1, az), 1, 3)- &
                                                             f2*kvec(2)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddksc)) THEN
                              ddsc(coset(ax, ay, az), 3, 1) = rbp(2)*ddsc(coset(ax, ay, az), 1, 1)+ &
                                                              fy*ddsc(coset(ax, ay-1, az), 1, 1)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 1)
                              ddsc(coset(ax, ay, az), 3, 2) = rbp(2)*ddsc(coset(ax, ay, az), 1, 2)+ &
                                                              fy*ddsc(coset(ax, ay-1, az), 1, 2)- &
                                                              f2*dss(coset(ax, ay, az), 1, 1)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 2)
                              ddsc(coset(ax, ay, az), 3, 3) = rbp(2)*ddsc(coset(ax, ay, az), 1, 3)+ &
                                                              fy*ddsc(coset(ax, ay-1, az), 1, 3)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 3)
                              ddsc(coset(ax, ay, az), 3, 4) = rbp(2)*ddsc(coset(ax, ay, az), 1, 4)+ &
                                                              fy*ddsc(coset(ax, ay-1, az), 1, 4)- &
                                                              2*f2*dss(coset(ax, ay, az), 1, 2)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 4)
                              ddsc(coset(ax, ay, az), 3, 5) = rbp(2)*ddsc(coset(ax, ay, az), 1, 5)+ &
                                                              fy*ddsc(coset(ax, ay-1, az), 1, 5)- &
                                                              f2*dss(coset(ax, ay, az), 1, 3)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 5)
                              ddsc(coset(ax, ay, az), 3, 6) = rbp(2)*ddsc(coset(ax, ay, az), 1, 6)+ &
                                                              fy*ddsc(coset(ax, ay-1, az), 1, 6)- &
                                                              f2*kvec(2)*ddss(coset(ax, ay, az), 1, 6)
                           END IF
                           IF (PRESENT(dkss)) THEN
                              dss(coset(ax, ay, az), 3, 1) = rbp(2)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             fy*dss(coset(ax, ay-1, az), 1, 1)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 3, 2) = rbp(2)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             fy*dss(coset(ax, ay-1, az), 1, 2)+ &
                                                             f2*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 3, 3) = rbp(2)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             fy*dss(coset(ax, ay-1, az), 1, 3)+ &
                                                             f2*kvec(2)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddkss)) THEN
                              ddss(coset(ax, ay, az), 3, 1) = rbp(2)*ddss(coset(ax, ay, az), 1, 1)+ &
                                                              fy*ddss(coset(ax, ay-1, az), 1, 1)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 1)
                              ddss(coset(ax, ay, az), 3, 2) = rbp(2)*ddss(coset(ax, ay, az), 1, 2)+ &
                                                              fy*ddss(coset(ax, ay-1, az), 1, 2)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 1)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 2)
                              ddss(coset(ax, ay, az), 3, 3) = rbp(2)*ddss(coset(ax, ay, az), 1, 3)+ &
                                                              fy*ddss(coset(ax, ay-1, az), 1, 3)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 3)
                              ddss(coset(ax, ay, az), 3, 4) = rbp(2)*ddss(coset(ax, ay, az), 1, 4)+ &
                                                              fy*ddss(coset(ax, ay-1, az), 1, 4)+ &
                                                              2*f2*dsc(coset(ax, ay, az), 1, 2)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 4)
                              ddss(coset(ax, ay, az), 3, 5) = rbp(2)*ddss(coset(ax, ay, az), 1, 5)+ &
                                                              fy*ddss(coset(ax, ay-1, az), 1, 5)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 3)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 5)
                              ddss(coset(ax, ay, az), 3, 6) = rbp(2)*ddss(coset(ax, ay, az), 1, 6)+ &
                                                              fy*ddss(coset(ax, ay-1, az), 1, 6)+ &
                                                              f2*kvec(2)*ddsc(coset(ax, ay, az), 1, 6)
                           END IF
                                                
                        END IF
                        IF (az == 0) THEN
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------

                           IF (PRESENT(dksc)) THEN                        
                              dsc(coset(ax, ay, az), 4, 1) = rbp(3)*dsc(coset(ax, ay, az), 1, 1)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 4, 2) = rbp(3)*dsc(coset(ax, ay, az), 1, 2)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 4, 3) = rbp(3)*dsc(coset(ax, ay, az), 1, 3)- &
                                                             f2*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddksc)) THEN
                              ddsc(coset(ax, ay, az), 4, 1) = rbp(3)*ddsc(coset(ax, ay, az), 1, 1)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 1)
                              ddsc(coset(ax, ay, az), 4, 2) = rbp(3)*ddsc(coset(ax, ay, az), 1, 2)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 2)
                              ddsc(coset(ax, ay, az), 4, 3) = rbp(3)*ddsc(coset(ax, ay, az), 1, 3)- &
                                                              f2*dss(coset(ax, ay, az), 1, 1)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 3)
                              ddsc(coset(ax, ay, az), 4, 4) = rbp(3)*ddsc(coset(ax, ay, az), 1, 4)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 4)
                              ddsc(coset(ax, ay, az), 4, 5) = rbp(3)*ddsc(coset(ax, ay, az), 1, 5)- &
                                                              f2*dss(coset(ax, ay, az), 1, 2)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 5)
                              ddsc(coset(ax, ay, az), 4, 6) = rbp(3)*ddsc(coset(ax, ay, az), 1, 6)- &
                                                              2*f2*dss(coset(ax, ay, az), 1, 3)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 6)
                           END IF
                           IF (PRESENT(dkss)) THEN                        
                              dss(coset(ax, ay, az), 4, 1) = rbp(3)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 4, 2) = rbp(3)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 4, 3) = rbp(3)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             f2*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddkss)) THEN
                              ddss(coset(ax, ay, az), 4, 1) = rbp(3)*ddss(coset(ax, ay, az), 1, 1)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 1)
                              ddss(coset(ax, ay, az), 4, 2) = rbp(3)*ddss(coset(ax, ay, az), 1, 2)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 2)
                              ddss(coset(ax, ay, az), 4, 3) = rbp(3)*ddss(coset(ax, ay, az), 1, 3)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 1)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 3)
                              ddss(coset(ax, ay, az), 4, 4) = rbp(3)*ddss(coset(ax, ay, az), 1, 4)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 4)
                              ddss(coset(ax, ay, az), 4, 5) = rbp(3)*ddss(coset(ax, ay, az), 1, 5)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 2)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 5)
                              ddss(coset(ax, ay, az), 4, 6) = rbp(3)*ddss(coset(ax, ay, az), 1, 6)+ &
                                                              2*f2*dsc(coset(ax, ay, az), 1, 3)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 6)
                           END IF
                        
                        ELSE
                           !------------------------------------------------
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1)+ &
                                                      fz*sc(coset(ax, ay, az-1), 1)- &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1)+ &
                                                      fz*ss(coset(ax, ay, az-1), 1)+ &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)
                           !------------------------------------------------
                  
                                          
                           IF (PRESENT(dksc)) THEN                        
                              dsc(coset(ax, ay, az), 4, 1) = rbp(3)*dsc(coset(ax, ay, az), 1, 1)+ &
                                                             fz*dsc(coset(ax, ay, az-1), 1, 1)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 1)
                              dsc(coset(ax, ay, az), 4, 2) = rbp(3)*dsc(coset(ax, ay, az), 1, 2)+ &
                                                             fz*dsc(coset(ax, ay, az-1), 1, 2)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 2)
                              dsc(coset(ax, ay, az), 4, 3) = rbp(3)*dsc(coset(ax, ay, az), 1, 3)+ &
                                                             fz*dsc(coset(ax, ay, az-1), 1, 3)- &
                                                             f2*ss(coset(ax, ay, az), 1)- &
                                                             f2*kvec(3)*dss(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddksc)) THEN
                              ddsc(coset(ax, ay, az), 4, 1) = rbp(3)*ddsc(coset(ax, ay, az), 1, 1)+ &
                                                              fz*ddsc(coset(ax, ay, az-1), 1, 1)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 1)
                              ddsc(coset(ax, ay, az), 4, 2) = rbp(3)*ddsc(coset(ax, ay, az), 1, 2)+ &
                                                              fz*ddsc(coset(ax, ay, az-1), 1, 2)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 2)
                              ddsc(coset(ax, ay, az), 4, 3) = rbp(3)*ddsc(coset(ax, ay, az), 1, 3)+ &
                                                              fz*ddsc(coset(ax, ay, az-1), 1, 3)- &
                                                              f2*dss(coset(ax, ay, az), 1, 1)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 3)
                              ddsc(coset(ax, ay, az), 4, 4) = rbp(3)*ddsc(coset(ax, ay, az), 1, 4)+ &
                                                              fz*ddsc(coset(ax, ay, az-1), 1, 4)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 4)
                              ddsc(coset(ax, ay, az), 4, 5) = rbp(3)*ddsc(coset(ax, ay, az), 1, 5)+ &
                                                              fz*ddsc(coset(ax, ay, az-1), 1, 5)- &
                                                              f2*dss(coset(ax, ay, az), 1, 2)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 5)
                              ddsc(coset(ax, ay, az), 4, 6) = rbp(3)*ddsc(coset(ax, ay, az), 1, 6)+ &
                                                              fz*ddsc(coset(ax, ay, az-1), 1, 6)- &
                                                              2*f2*dss(coset(ax, ay, az), 1, 3)- &
                                                              f2*kvec(3)*ddss(coset(ax, ay, az), 1, 6)
                           END IF
                           IF (PRESENT(dkss)) THEN                        
                              dss(coset(ax, ay, az), 4, 1) = rbp(3)*dss(coset(ax, ay, az), 1, 1)+ &
                                                             fz*dss(coset(ax, ay, az-1), 1, 1)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 1)
                              dss(coset(ax, ay, az), 4, 2) = rbp(3)*dss(coset(ax, ay, az), 1, 2)+ &
                                                             fz*dss(coset(ax, ay, az-1), 1, 2)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 2)
                              dss(coset(ax, ay, az), 4, 3) = rbp(3)*dss(coset(ax, ay, az), 1, 3)+ &
                                                             fz*dss(coset(ax, ay, az-1), 1, 3)+ &
                                                             f2*sc(coset(ax, ay, az), 1)+ &
                                                             f2*kvec(3)*dsc(coset(ax, ay, az), 1, 3)
                           END IF
                           IF (PRESENT(ddkss)) THEN
                              ddss(coset(ax, ay, az), 4, 1) = rbp(3)*ddss(coset(ax, ay, az), 1, 1)+ &
                                                              fz*ddss(coset(ax, ay, az-1), 1, 1)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 1)
                              ddss(coset(ax, ay, az), 4, 2) = rbp(3)*ddss(coset(ax, ay, az), 1, 2)+ &
                                                              fz*ddss(coset(ax, ay, az-1), 1, 2)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 2)
                              ddss(coset(ax, ay, az), 4, 3) = rbp(3)*ddss(coset(ax, ay, az), 1, 3)+ &
                                                              fz*ddss(coset(ax, ay, az-1), 1, 3)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 1)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 3)
                              ddss(coset(ax, ay, az), 4, 4) = rbp(3)*ddss(coset(ax, ay, az), 1, 4)+ &
                                                              fz*ddss(coset(ax, ay, az-1), 1, 4)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 4)
                              ddss(coset(ax, ay, az), 4, 5) = rbp(3)*ddss(coset(ax, ay, az), 1, 5)+ &
                                                              fz*ddss(coset(ax, ay, az-1), 1, 5)+ &
                                                              f2*dsc(coset(ax, ay, az), 1, 2)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 5)
                              ddss(coset(ax, ay, az), 4, 6) = rbp(3)*ddss(coset(ax, ay, az), 1, 6)+ &
                                                              fz*ddss(coset(ax, ay, az-1), 1, 6)+ &
                                                              2*f2*dsc(coset(ax, ay, az), 1, 3)+ &
                                                              f2*kvec(3)*ddsc(coset(ax, ay, az), 1, 6)
                           END IF
                        END IF
                     END DO
                  END DO

!           *** Recurrence steps: [a|O|p] -> [a|O|b] ***

                  DO lb = 2, lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|O|b] = [a+1i|O|b-1i] - (Bi - Ai)*[a|O|b-1i] ***

                     IF (lb == lb_max) THEN
                        la_start = la_min
                     ELSE
                        la_start = MAX(0, la_min-1)
                     END IF

                     DO la = la_start, la_max-1
                        DO ax = 0, la
                           DO ay = 0, la-ax
                              az = la-ax-ay

!                   *** Shift of angular momentum component z from a to b ***

                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 sc(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                 rab(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 ss(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                 rab(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))

                              IF (PRESENT(dksc)) THEN
                                 DO i = 1, 3
                                    dsc(coset(ax, ay, az), coset(0, 0, lb), i) = &
                                       dsc(coset(ax, ay, az+1), coset(0, 0, lb-1), i)- &
                                       rab(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), i)
                                 END DO
                              END IF
                              IF (PRESENT(dkss)) THEN
                                 DO i = 1, 3
                                    dss(coset(ax, ay, az), coset(0, 0, lb), i) = &
                                       dss(coset(ax, ay, az+1), coset(0, 0, lb-1), i)- &
                                       rab(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), i)
                                 END DO
                              END IF
                              IF (PRESENT(ddksc)) THEN
                                 DO i = 1, 6
                                    ddsc(coset(ax, ay, az), coset(0, 0, lb), i) = &
                                       ddsc(coset(ax, ay, az+1), coset(0, 0, lb-1), i)- &
                                       rab(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), i)
                                 END DO
                              END IF
                              IF (PRESENT(ddkss)) THEN
                                 DO i = 1, 6
                                    ddss(coset(ax, ay, az), coset(0, 0, lb), i) = &
                                       ddss(coset(ax, ay, az+1), coset(0, 0, lb-1), i)- &
                                       rab(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), i)
                                 END DO
                              END IF

!                   *** Shift of angular momentum component y from a to b ***

                              DO by = 1, lb
                                 bz = lb-by
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    sc(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                    rab(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    ss(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                    rab(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))

                                 IF (PRESENT(dksc)) THEN
                                    DO i = 1, 3
                                       dsc(coset(ax, ay, az), coset(0, by, bz), i) = &
                                          dsc(coset(ax, ay+1, az), coset(0, by-1, bz), i)- &
                                          rab(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), i)
                                    END DO
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    DO i = 1, 3
                                       dss(coset(ax, ay, az), coset(0, by, bz), i) = &
                                          dss(coset(ax, ay+1, az), coset(0, by-1, bz), i)- &
                                          rab(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), i)
                                    END DO
                                 END IF
                                 IF (PRESENT(ddksc)) THEN
                                    DO i = 1, 6
                                       ddsc(coset(ax, ay, az), coset(0, by, bz), i) = &
                                          ddsc(coset(ax, ay+1, az), coset(0, by-1, bz), i)- &
                                          rab(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), i)
                                    END DO
                                 END IF
                                 IF (PRESENT(ddkss)) THEN
                                    DO i = 1, 6
                                       ddss(coset(ax, ay, az), coset(0, by, bz), i) = &
                                          ddss(coset(ax, ay+1, az), coset(0, by-1, bz), i)- &
                                          rab(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), i)
                                    END DO
                                 END IF

                              END DO

!                   *** Shift of angular momentum component x from a to b ***

                              DO bx = 1, lb
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       sc(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                       rab(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       ss(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                       rab(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))

                                    IF (PRESENT(dksc)) THEN
                                       DO i = 1, 3
                                          dsc(coset(ax, ay, az), coset(bx, by, bz), i) = &
                                             dsc(coset(ax+1, ay, az), coset(bx-1, by, bz), i)- &
                                             rab(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), i)
                                       END DO
                                    END IF
                                    IF (PRESENT(dkss)) THEN
                                       DO i = 1, 3
                                          dss(coset(ax, ay, az), coset(bx, by, bz), i) = &
                                             dss(coset(ax+1, ay, az), coset(bx-1, by, bz), i)- &
                                             rab(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), i)
                                       END DO
                                    END IF

                                    IF (PRESENT(ddksc)) THEN
                                       DO i = 1, 6
                                          ddsc(coset(ax, ay, az), coset(bx, by, bz), i) = &
                                             ddsc(coset(ax+1, ay, az), coset(bx-1, by, bz), i)- &
                                             rab(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), i)
                                       END DO
                                    END IF
                                    IF (PRESENT(ddkss)) THEN
                                       DO i = 1, 6
                                          ddss(coset(ax, ay, az), coset(bx, by, bz), i) = &
                                             ddss(coset(ax+1, ay, az), coset(bx-1, by, bz), i)- &
                                             rab(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), i)
                                       END DO
                                    END IF

                                 END DO
                              END DO

                           END DO
                        END DO
                     END DO

!             *** Vertical recurrence step ***

!             *** [a|O|b] = (Pi - Bi)*[a|O|b-1i] + f2*Ni(a)*[a-1i|O|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|O|b-2i] + [a|dO|b-1i]            ***

                     DO ax = 0, la_max
                        fx = f2*REAL(ax, dp)
                        DO ay = 0, la_max-ax
                           fy = f2*REAL(ay, dp)
                           az = la_max-ax-ay
                           fz = f2*REAL(az, dp)

!                 *** Increase the angular momentum component z of function b ***

                           f3 = f2*REAL(lb-1, dp)

                           IF (az == 0) THEN
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              !------------------------------------------------
        
                              IF (PRESENT(dksc)) THEN
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                    rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                    f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 1)- &
                                    f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                    rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                    f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 2)- &
                                    f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                    rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                    f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 3)- &
                                    f2*ss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                    f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                              END IF
                              IF (PRESENT(ddksc)) THEN
                                 ddsc(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                      rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                      f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 1)- &
                                      f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                 ddsc(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                      rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                      f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 2)- &
                                      f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 3)- &
                                       f2*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 4) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 4)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 4)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 4)
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 5) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 5)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 5)- &
                                       f2*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 5)
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 6) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 6)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 6)- &
                                       2*f2*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 6)
                              END IF
                              IF (PRESENT(dkss)) THEN
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 1)+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 2)+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 3)+ &
                                     f2*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                              END IF
                              IF (PRESENT(ddkss)) THEN
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 1)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 2)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 3)+ &
                                       f2*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 4) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 4)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 4)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 4)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 5) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 5)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 5)+ &
                                       f2*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 5)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 6) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 6)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 6)+ &
                                       2*f2*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 6)
                              END IF

                           ELSE
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 fz*sc(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 fz*ss(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              !------------------------------------------------

                              IF (PRESENT(dksc)) THEN
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                     rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                     fz*dsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 1)+ &
                                     f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 1)- &
                                     f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                     rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                     fz*dsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 2)+ &
                                     f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 2)- &
                                     f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                 dsc(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                     rbp(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                     fz*dsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 3)+ &
                                     f3*dsc(coset(ax, ay, az), coset(0, 0, lb-2), 3)- &
                                     f2*ss(coset(ax, ay, az), coset(0, 0, lb-1))- &
                                     f2*kvec(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                              END IF
                              IF (PRESENT(ddksc)) THEN
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                       fz*ddsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 1)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 1)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                       fz*ddsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 2)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 2)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                       fz*ddsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 3)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 3)- &
                                       f2*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 4) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 4)+ &
                                       fz*ddsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 4)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 4)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 4)
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 5) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 5)+ &
                                       fz*ddsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 5)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 5)- &
                                       f2*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 5)
                                  ddsc(coset(ax, ay, az), coset(0, 0, lb), 6) = &
                                       rbp(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 6)+ &
                                       fz*ddsc(coset(ax, ay, az-1), coset(0, 0, lb-1), 6)+ &
                                       f3*ddsc(coset(ax, ay, az), coset(0, 0, lb-2), 6)- &
                                       2*f2*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)- &
                                       f2*kvec(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 6)
                              END IF

                              IF (PRESENT(dkss)) THEN
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                     fz*dss(coset(ax, ay, az-1), coset(0, 0, lb-1), 1)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 1)+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                     fz*dss(coset(ax, ay, az-1), coset(0, 0, lb-1), 2)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 2)+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                 dss(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                     rbp(3)*dss(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                     fz*dss(coset(ax, ay, az-1), coset(0, 0, lb-1), 3)+ &
                                     f3*dss(coset(ax, ay, az), coset(0, 0, lb-2), 3)+ &
                                     f2*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                     f2*kvec(3)*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                              END IF

                              IF (PRESENT(ddkss)) THEN
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 1) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                       fz*ddss(coset(ax, ay, az-1), coset(0, 0, lb-1), 1)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 1)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 2) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                       fz*ddss(coset(ax, ay, az-1), coset(0, 0, lb-1), 2)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 2)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 3) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                       fz*ddss(coset(ax, ay, az-1), coset(0, 0, lb-1), 3)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 3)+ &
                                       f2*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 1)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 4) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 4)+ &
                                       fz*ddss(coset(ax, ay, az-1), coset(0, 0, lb-1), 4)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 4)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 4)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 5) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 5)+ &
                                       fz*ddss(coset(ax, ay, az-1), coset(0, 0, lb-1), 5)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 5)+ &
                                       f2*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 2)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 5)
                                  ddss(coset(ax, ay, az), coset(0, 0, lb), 6) = &
                                       rbp(3)*ddss(coset(ax, ay, az), coset(0, 0, lb-1), 6)+ &
                                       fz*ddss(coset(ax, ay, az-1), coset(0, 0, lb-1), 6)+ &
                                       f3*ddss(coset(ax, ay, az), coset(0, 0, lb-2), 6)+ &
                                       2*f2*dsc(coset(ax, ay, az), coset(0, 0, lb-1), 3)+ &
                                       f2*kvec(3)*ddsc(coset(ax, ay, az), coset(0, 0, lb-1), 6)
                              END IF

                           END IF

!                 *** Increase the angular momentum component y of function b ***

                           IF (ay == 0) THEN
                              bz = lb-1
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz))- &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              !------------------------------------------------

                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                        f2*ss(coset(ax, ay, az), coset(0, 0, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 3)
                                 END IF
                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                         f2*dss(coset(ax, ay, az), coset(0, 0, bz), 1)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 3)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 4) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 4)+ &
                                         2*f2*dss(coset(ax, ay, az), coset(0, 0, bz), 2)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 4)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 5) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 5)+ &
                                         f2*dss(coset(ax, ay, az), coset(0, 0, bz), 3)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 5)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 6) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 6)+ &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 6)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                        f2*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)
                                 END IF

                                 IF (PRESENT(ddkss)) THEN
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 3)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 4) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 4)+ &
                                         2*f2*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 4)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 5) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 5)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 5)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 6) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 6)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 6)
                                 END IF

                              DO by = 2, lb
                                 bz = lb-by
                                 f3 = f2*REAL(by-1, dp)
                                 !------------------------------------------------
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    f3*sc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    f3*ss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 !------------------------------------------------
                              
                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 1)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 2)- &
                                        f2*ss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 3)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                 END IF
                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 1)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 2)- &
                                         f2*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 3)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 4) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 4)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 4)- &
                                         2*f2*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 4)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 5) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 5)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 5)- &
                                         f2*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 5)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 6) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 6)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 6)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 6)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 1)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    dss(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 2)+ &
                                        f2*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    dss(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 3)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                 END IF
                                 IF (PRESENT(ddkss)) THEN
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 1)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 2)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 3)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 4) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 4)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 4)+ &
                                         2*f2*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 4)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 5) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 5)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 5)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 5)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 6) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 6)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 6)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 6)
                                 END IF
                              
                              END DO
                           ELSE
                              bz = lb-1
                              !------------------------------------------------
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 fy*sc(coset(ax, ay-1, az), coset(0, 0, bz))- &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 fy*ss(coset(ax, ay-1, az), coset(0, 0, bz))+ &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              !------------------------------------------------

                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, 0, bz), 1)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, 0, bz), 2)- &
                                        f2*ss(coset(ax, ay, az), coset(0, 0, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    dsc(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, 0, bz), 3)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 3)
                                 END IF
                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, 0, bz), 1)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, 0, bz), 2)- &
                                         f2*dss(coset(ax, ay, az), coset(0, 0, bz), 1)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, 0, bz), 3)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 3)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 4) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 4)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, 0, bz), 4)- &
                                         2*f2*dss(coset(ax, ay, az), coset(0, 0, bz), 2)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 4)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 5) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 5)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, 0, bz), 5)- &
                                         f2*dss(coset(ax, ay, az), coset(0, 0, bz), 3)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 5)
                                    ddsc(coset(ax, ay, az), coset(0, 1, bz), 6) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 6)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, 0, bz), 6)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 6)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, 0, bz), 1)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, 0, bz), 2)+ &
                                        f2*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    dss(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, 0, bz), 3)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)
                                 END IF
                                 IF (PRESENT(ddkss)) THEN
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 1) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, 0, bz), 1)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 1)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 2) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, 0, bz), 2)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, 0, bz), 1)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 2)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 3) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, 0, bz), 3)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 3)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 4) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 4)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, 0, bz), 4)+ &
                                         2*f2*dsc(coset(ax, ay, az), coset(0, 0, bz), 2)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 4)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 5) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 5)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, 0, bz), 5)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, 0, bz), 3)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 5)
                                    ddss(coset(ax, ay, az), coset(0, 1, bz), 6) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, 0, bz), 6)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, 0, bz), 6)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, 0, bz), 6)
                                 END IF

                              DO by = 2, lb
                                 bz = lb-by
                                 f3 = f2*REAL(by-1, dp)
                                 !------------------------------------------------
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    fy*sc(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                    f3*sc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    fy*ss(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                    f3*ss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 !------------------------------------------------

                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, by-1, bz), 1)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 1)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, by-1, bz), 2)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 2)- &
                                        f2*ss(coset(ax, ay, az), coset(0, by-1, bz))- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    dsc(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                        rbp(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                        fy*dsc(coset(ax, ay-1, az), coset(0, by-1, bz), 3)+ &
                                        f3*dsc(coset(ax, ay, az), coset(0, by-2, bz), 3)- &
                                        f2*kvec(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                 END IF
                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, by-1, bz), 1)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 1)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, by-1, bz), 2)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 2)- &
                                         f2*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, by-1, bz), 3)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 3)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 4) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 4)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, by-1, bz), 4)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 4)- &
                                         2*f2*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 4)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 5) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 5)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, by-1, bz), 5)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 5)- &
                                         f2*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 5)
                                    ddsc(coset(ax, ay, az), coset(0, by, bz), 6) = &
                                         rbp(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 6)+ &
                                         fy*ddsc(coset(ax, ay-1, az), coset(0, by-1, bz), 6)+ &
                                         f3*ddsc(coset(ax, ay, az), coset(0, by-2, bz), 6)- &
                                         f2*kvec(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 6)
                                 END IF

                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, by-1, bz), 1)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 1)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    dss(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, by-1, bz), 2)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 2)+ &
                                        f2*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    dss(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                        rbp(2)*dss(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                        fy*dss(coset(ax, ay-1, az), coset(0, by-1, bz), 3)+ &
                                        f3*dss(coset(ax, ay, az), coset(0, by-2, bz), 3)+ &
                                        f2*kvec(2)*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                 END IF
                                 IF (PRESENT(ddkss)) THEN
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 1) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, by-1, bz), 1)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 1)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 1)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 2) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, by-1, bz), 2)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 2)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, by-1, bz), 1)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 2)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 3) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, by-1, bz), 3)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 3)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 3)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 4) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 4)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, by-1, bz), 4)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 4)+ &
                                         2*f2*dsc(coset(ax, ay, az), coset(0, by-1, bz), 2)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 4)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 5) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 5)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, by-1, bz), 5)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 5)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, by-1, bz), 3)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 5)
                                    ddss(coset(ax, ay, az), coset(0, by, bz), 6) = &
                                         rbp(2)*ddss(coset(ax, ay, az), coset(0, by-1, bz), 6)+ &
                                         fy*ddss(coset(ax, ay-1, az), coset(0, by-1, bz), 6)+ &
                                         f3*ddss(coset(ax, ay, az), coset(0, by-2, bz), 6)+ &
                                         f2*kvec(2)*ddsc(coset(ax, ay, az), coset(0, by-1, bz), 6)
                                 END IF

                              END DO
                           END IF

!                 *** Increase the angular momentum component x of function b ***

                           IF (ax == 0) THEN
                              DO by = 0, lb-1
                                 bz = lb-1-by
                                 !----------------------------------------
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz))- &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz))+ &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                                 !----------------------------------------
                              
                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 1)- &
                                        f2*ss(coset(ax, ay, az), coset(0, by, bz))- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 1)
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 2)- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 2)
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 3)- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 3)
                                 END IF

                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 1)- &
                                         2*f2*dss(coset(ax, ay, az), coset(0, by, bz), 1)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 1)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 2)- &
                                         f2*dss(coset(ax, ay, az), coset(0, by, bz), 2)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 2)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 3)- &
                                         f2*dss(coset(ax, ay, az), coset(0, by, bz), 3)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 3)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 4) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                    coset(0, by, bz), 4)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 4)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 5) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 5)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 5)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 6) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 6)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 6)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 1)+ &
                                        f2*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 1)
                                    dss(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 2)+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 2)
                                    dss(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 3)+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 3)
                                 END IF
                                 IF (PRESENT(ddkss)) THEN
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 1)+ &
                                         2*f2*dsc(coset(ax, ay, az), coset(0, by, bz), 1)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 1)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 2)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, by, bz), 2)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 2)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 3)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, by, bz), 3)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 3)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 4) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 4)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 4)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 5) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 5)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 5)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 6) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 6)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 6)
                                 END IF
                              END DO
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx-1, dp)
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    !----------------------------------------
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       f3*sc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       f3*ss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    !----------------------------------------

                                    IF (PRESENT(dksc)) THEN
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 1)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 1)- &
                                           f2*ss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 2)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 2)- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 3)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 3)- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                    END IF
                                    IF (PRESENT(ddksc)) THEN
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 1)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 1)- &
                                            2*f2*dss(coset(ax, ay, az), coset(bx-1, by, bz), 1)- & 
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 2)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 2)- &
                                            f2*dss(coset(ax, ay, az), coset(bx-1, by, bz), 2)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 3)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 3)- &
                                            f2*dss(coset(ax, ay, az), coset(bx-1, by, bz), 3)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 4) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 4)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 4)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 4)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 5) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 5)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 5)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 5)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 6) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 6)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 6)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 6)
                                    END IF
                                    IF (PRESENT(dkss)) THEN
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 1)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 1)+ &
                                           f2*sc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 2)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 2)+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 3)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 3)+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                    END IF
                                    IF (PRESENT(ddkss)) THEN
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 1)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 1)+ &
                                            2*f2*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 1)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 2)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 2)+ &
                                            f2*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 2)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 3)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 3)+ &
                                            f2*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 3)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 4) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 4)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 4)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 4)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 5) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 5)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 5)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 5)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 6) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 6)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 6)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 6)
                                    END IF
                                 
                                 END DO
                              END DO
                           ELSE
                              DO by = 0, lb-1
                                 bz = lb-1-by
                                 !----------------------------------------
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                    fx*sc(coset(ax-1, ay, az), coset(0, by, bz))- &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz))+ &
                                    fx*ss(coset(ax-1, ay, az), coset(0, by, bz))+ &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                                 !----------------------------------------
                              
                                 IF (PRESENT(dksc)) THEN
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 1)+ &
                                        fx*dsc(coset(ax-1, ay, az), coset(0, by, bz), 1)- &
                                        f2*ss(coset(ax, ay, az), coset(0, by, bz))- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 1)
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 2)+ &
                                        fx*dsc(coset(ax-1, ay, az), coset(0, by, bz), 2)- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 2)
                                    dsc(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                        rbp(1)*dsc(coset(ax, ay, az), &
                                                   coset(0, by, bz), 3)+ &
                                        fx*dsc(coset(ax-1, ay, az), coset(0, by, bz), 3)- &
                                        f2*kvec(1)*dss(coset(ax, ay, az), coset(0, by, bz), 3)
                                 END IF
                                 IF (PRESENT(ddksc)) THEN
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 1)+ &
                                         fx*ddsc(coset(ax-1, ay, az), coset(0, by, bz), 1)- &
                                         2*f2*dss(coset(ax, ay, az), coset(0, by, bz), 1)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 1)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 2)+ &
                                         fx*ddsc(coset(ax-1, ay, az), coset(0, by, bz), 2)- &
                                         f2*dss(coset(ax, ay, az), coset(0, by, bz), 2)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 2)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 3)+ &
                                         fx*ddsc(coset(ax-1, ay, az), coset(0, by, bz), 3)- &
                                         f2*dss(coset(ax, ay, az), coset(0, by, bz), 3)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 3)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 4) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                    coset(0, by, bz), 4)+ &
                                         fx*ddsc(coset(ax-1, ay, az), coset(0, by, bz), 4)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 4)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 5) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 5)+ &
                                         fx*ddsc(coset(ax-1, ay, az), coset(0, by, bz), 5)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 5)
                                    ddsc(coset(ax, ay, az), coset(1, by, bz), 6) = &
                                         rbp(1)*ddsc(coset(ax, ay, az), &
                                                     coset(0, by, bz), 6)+ &
                                         fx*ddsc(coset(ax-1, ay, az), coset(0, by, bz), 6)- &
                                         f2*kvec(1)*ddss(coset(ax, ay, az), coset(0, by, bz), 6)
                                 END IF
                                 IF (PRESENT(dkss)) THEN
                                    dss(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 1)+ &
                                        fx*dss(coset(ax-1, ay, az), coset(0, by, bz), 1)+ &
                                        f2*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 1)
                                    dss(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 2)+ &
                                        fx*dss(coset(ax-1, ay, az), coset(0, by, bz), 2)+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 2)
                                    dss(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                        rbp(1)*dss(coset(ax, ay, az), &
                                                   coset(0, by, bz), 3)+ &
                                        fx*dss(coset(ax-1, ay, az), coset(0, by, bz), 3)+ &
                                        f2*kvec(1)*dsc(coset(ax, ay, az), coset(0, by, bz), 3)
                                 END IF
                                 IF (PRESENT(ddkss)) THEN
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 1) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 1)+ &
                                         fx*ddss(coset(ax-1, ay, az), coset(0, by, bz), 1)+ &
                                         2*f2*dsc(coset(ax, ay, az), coset(0, by, bz), 1)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 1)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 2) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 2)+ &
                                         fx*ddss(coset(ax-1, ay, az), coset(0, by, bz), 2)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, by, bz), 2)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 2)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 3) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 3)+ &
                                         fx*ddss(coset(ax-1, ay, az), coset(0, by, bz), 3)+ &
                                         f2*dsc(coset(ax, ay, az), coset(0, by, bz), 3)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 3)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 4) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 4)+ &
                                         fx*ddss(coset(ax-1, ay, az), coset(0, by, bz), 4)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 4)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 5) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 5)+ &
                                         fx*ddss(coset(ax-1, ay, az), coset(0, by, bz), 5)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 5)
                                    ddss(coset(ax, ay, az), coset(1, by, bz), 6) = &
                                         rbp(1)*ddss(coset(ax, ay, az), &
                                                     coset(0, by, bz), 6)+ &
                                         fx*ddss(coset(ax-1, ay, az), coset(0, by, bz), 6)+ &
                                         f2*kvec(1)*ddsc(coset(ax, ay, az), coset(0, by, bz), 6)
                                 END IF
                              END DO
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx-1, dp)
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    !----------------------------------------
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       fx*sc(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                       f3*sc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       fx*ss(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                       f3*ss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    !----------------------------------------
                                 
                                    IF (PRESENT(dksc)) THEN
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 1)+ &
                                           fx*dsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 1)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 1)- &
                                           f2*ss(coset(ax, ay, az), coset(bx-1, by, bz))- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 2)+ &
                                           fx*dsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 2)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 2)- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       dsc(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                           rbp(1)*dsc(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 3)+ &
                                           fx*dsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 3)+ &
                                           f3*dsc(coset(ax, ay, az), coset(bx-2, by, bz), 3)- &
                                           f2*kvec(1)*dss(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                    END IF
                                    IF (PRESENT(ddksc)) THEN
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 1)+ &
                                            fx*ddsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 1)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 1)- &
                                            2*f2*dss(coset(ax, ay, az), coset(bx-1, by, bz), 1)- & 
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 2)+ &
                                            fx*ddsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 2)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 2)- &
                                            f2*dss(coset(ax, ay, az), coset(bx-1, by, bz), 2)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 3)+ &
                                            fx*ddsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 3)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 3)- &
                                            f2*dss(coset(ax, ay, az), coset(bx-1, by, bz), 3)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 4) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 4)+ &
                                            fx*ddsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 4)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 4)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 4)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 5) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 5)+ &
                                            fx*ddsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 5)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 5)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 5)
                                       ddsc(coset(ax, ay, az), coset(bx, by, bz), 6) = &
                                            rbp(1)*ddsc(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 6)+ &
                                            fx*ddsc(coset(ax-1, ay, az), coset(bx-1, by, bz), 6)+ &
                                            f3*ddsc(coset(ax, ay, az), coset(bx-2, by, bz), 6)- &
                                            f2*kvec(1)*ddss(coset(ax, ay, az), coset(bx-1, by, bz), 6)
                                    END IF
                                    IF (PRESENT(dkss)) THEN
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 1)+ &
                                           fx*dss(coset(ax-1, ay, az), coset(bx-1, by, bz), 1)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 1)+ &
                                           f2*sc(coset(ax, ay, az), coset(bx-1, by, bz))+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 2)+ &
                                           fx*dss(coset(ax-1, ay, az), coset(bx-1, by, bz), 2)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 2)+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       dss(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                           rbp(1)*dss(coset(ax, ay, az), &
                                                     coset(bx-1, by, bz), 3)+ &
                                           fx*dss(coset(ax-1, ay, az), coset(bx-1, by, bz), 3)+ &
                                           f3*dss(coset(ax, ay, az), coset(bx-2, by, bz), 3)+ &
                                           f2*kvec(1)*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                    END IF
                                    IF (PRESENT(ddkss)) THEN
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 1) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 1)+ &
                                            fx*ddss(coset(ax-1, ay, az), coset(bx-1, by, bz), 1)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 1)+ &
                                            2*f2*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 1)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 1)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 2) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 2)+ &
                                            fx*ddss(coset(ax-1, ay, az), coset(bx-1, by, bz), 2)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 2)+ &
                                            f2*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 2)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 2)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 3) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                        coset(bx-1, by, bz), 3)+ &
                                            fx*ddss(coset(ax-1, ay, az), coset(bx-1, by, bz), 3)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 3)+ &
                                            f2*dsc(coset(ax, ay, az), coset(bx-1, by, bz), 3)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 3)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 4) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 4)+ &
                                            fx*ddss(coset(ax-1, ay, az), coset(bx-1, by, bz), 4)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 4)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 4)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 5) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 5)+ &
                                            fx*ddss(coset(ax-1, ay, az), coset(bx-1, by, bz), 5)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 5)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 5)
                                       ddss(coset(ax, ay, az), coset(bx, by, bz), 6) = &
                                            rbp(1)*ddss(coset(ax, ay, az), &
                                                       coset(bx-1, by, bz), 6)+ &
                                            fx*ddss(coset(ax-1, ay, az), coset(bx-1, by, bz), 6)+ &
                                            f3*ddss(coset(ax, ay, az), coset(bx-2, by, bz), 6)+ &
                                            f2*kvec(1)*ddsc(coset(ax, ay, az), coset(bx-1, by, bz), 6)
                                    END IF
                                 
                                 END DO
                              END DO
                           END IF

                        END DO
                     END DO

                  END DO

               END IF

            ELSE

               IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|O|s] -> [s|O|b] ***

                  rbp(:) = (f1-1.0_dp)*rab(:)

!           *** [s|O|p] = (Pi - Bi)*[s|O|s] + [s|dO|s] ***

                  !----------------------------------------
                  sc(1, 2) = rbp(1)*sc(1, 1)-f2*kvec(1)*ss(1, 1)
                  sc(1, 3) = rbp(2)*sc(1, 1)-f2*kvec(2)*ss(1, 1)
                  sc(1, 4) = rbp(3)*sc(1, 1)-f2*kvec(3)*ss(1, 1)
                  ss(1, 2) = rbp(1)*ss(1, 1)+f2*kvec(1)*sc(1, 1)
                  ss(1, 3) = rbp(2)*ss(1, 1)+f2*kvec(2)*sc(1, 1)
                  ss(1, 4) = rbp(3)*ss(1, 1)+f2*kvec(3)*sc(1, 1)
                  !----------------------------------------

               IF (PRESENT(dksc)) THEN
                  dsc(1,2,1) = rbp(1)*dsc(1,1,1) - f2*ss(1,1) - f2*kvec(1)*dss(1,1,1) 
                  dsc(1,2,2) = rbp(1)*dsc(1,1,2) - f2*kvec(1)*dss(1,1,2)
                  dsc(1,2,3) = rbp(1)*dsc(1,1,3) - f2*kvec(1)*dss(1,1,3)
                  
                  dsc(1,3,1) = rbp(2)*dsc(1,1,1) - f2*kvec(2)*dss(1,1,1) 
                  dsc(1,3,2) = rbp(2)*dsc(1,1,2) - f2*ss(1,1) - f2*kvec(2)*dss(1,1,2)
                  dsc(1,3,3) = rbp(2)*dsc(1,1,3) - f2*kvec(2)*dss(1,1,3)
                  
                  dsc(1,4,1) = rbp(3)*dsc(1,1,1) - f2*kvec(3)*dss(1,1,1) 
                  dsc(1,4,2) = rbp(3)*dsc(1,1,2) - f2*kvec(3)*dss(1,1,2)
                  dsc(1,4,3) = rbp(3)*dsc(1,1,3) - f2*ss(1,1) - f2*kvec(3)*dss(1,1,3)
               END IF

               IF (PRESENT(dkss)) THEN
                  dss(1,2,1) = rbp(1)*dss(1,1,1) + f2*sc(1,1) + f2*kvec(1)*dsc(1,1,1) 
                  dss(1,2,2) = rbp(1)*dss(1,1,2) + f2*kvec(1)*dsc(1,1,2)
                  dss(1,2,3) = rbp(1)*dss(1,1,3) + f2*kvec(1)*dsc(1,1,3)

                  dss(1,3,1) = rbp(2)*dss(1,1,1) + f2*kvec(2)*dsc(1,1,1) 
                  dss(1,3,2) = rbp(2)*dss(1,1,2) + f2*sc(1,1) + f2*kvec(2)*dsc(1,1,2)
                  dss(1,3,3) = rbp(2)*dss(1,1,3) + f2*kvec(2)*dsc(1,1,3)

                  dss(1,4,1) = rbp(3)*dss(1,1,1) + f2*kvec(3)*dsc(1,1,1) 
                  dss(1,4,2) = rbp(3)*dss(1,1,2) + f2*kvec(3)*dsc(1,1,2)
                  dss(1,4,3) = rbp(3)*dss(1,1,3) + f2*sc(1,1) + f2*kvec(3)*dsc(1,1,3)
               END IF

               IF (PRESENT(ddksc)) THEN
                  ddsc(1,2,1) = rbp(1)*ddsc(1,1,1) - 2*f2*dss(1,1,1) - f2*kvec(1)*ddss(1,1,1) 
                  ddsc(1,2,2) = rbp(1)*ddsc(1,1,2) - f2*dss(1,1,2) - f2*kvec(1)*ddss(1,1,2) 
                  ddsc(1,2,3) = rbp(1)*ddsc(1,1,3) - f2*dss(1,1,3) - f2*kvec(1)*ddss(1,1,3) 
                  ddsc(1,2,4) = rbp(1)*ddsc(1,1,4) - f2*kvec(1)*ddss(1,1,4) 
                  ddsc(1,2,5) = rbp(1)*ddsc(1,1,5) - f2*kvec(1)*ddss(1,1,5) 
                  ddsc(1,2,6) = rbp(1)*ddsc(1,1,6) - f2*kvec(1)*ddss(1,1,6) 
                  
                  ddsc(1,3,1) = rbp(2)*ddsc(1,1,1) - f2*kvec(2)*ddss(1,1,1)
                  ddsc(1,3,2) = rbp(2)*ddsc(1,1,2) - f2*dss(1,1,1) - f2*kvec(2)*ddss(1,1,2)
                  ddsc(1,3,3) = rbp(2)*ddsc(1,1,3) - f2*kvec(2)*ddss(1,1,3)
                  ddsc(1,3,4) = rbp(2)*ddsc(1,1,4) - 2*f2*dss(1,1,2) - f2*kvec(2)*ddss(1,1,4)
                  ddsc(1,3,5) = rbp(2)*ddsc(1,1,5) - f2*dss(1,1,3) - f2*kvec(2)*ddss(1,1,5) 
                  ddsc(1,3,6) = rbp(2)*ddsc(1,1,6) - f2*kvec(2)*ddss(1,1,6)
                  
                  ddsc(1,4,1) = rbp(3)*ddsc(1,1,1) - f2*kvec(3)*ddss(1,1,1)
                  ddsc(1,4,2) = rbp(3)*ddsc(1,1,2) - f2*kvec(3)*ddss(1,1,2)
                  ddsc(1,4,3) = rbp(3)*ddsc(1,1,3) - f2*dss(1,1,1) - f2*kvec(3)*ddss(1,1,3)
                  ddsc(1,4,4) = rbp(3)*ddsc(1,1,4) - f2*kvec(3)*ddss(1,1,4)
                  ddsc(1,4,5) = rbp(3)*ddsc(1,1,5) - f2*dss(1,1,2) - f2*kvec(3)*ddss(1,1,5)
                  ddsc(1,4,6) = rbp(3)*ddsc(1,1,6) - 2*f2*dss(1,1,3) - f2*kvec(3)*ddss(1,1,6)
               END IF

               IF (PRESENT(ddkss)) THEN
                  ddss(1,2,1) = rbp(1)*ddss(1,1,1) + 2*f2*dsc(1,1,1) + f2*kvec(1)*ddsc(1,1,1) 
                  ddss(1,2,2) = rbp(1)*ddss(1,1,2) + f2*dsc(1,1,2) + f2*kvec(1)*ddsc(1,1,2) 
                  ddss(1,2,3) = rbp(1)*ddss(1,1,3) + f2*dsc(1,1,3) + f2*kvec(1)*ddsc(1,1,3) 
                  ddss(1,2,4) = rbp(1)*ddss(1,1,4) + f2*kvec(1)*ddsc(1,1,4) 
                  ddss(1,2,5) = rbp(1)*ddss(1,1,5) + f2*kvec(1)*ddsc(1,1,5) 
                  ddss(1,2,6) = rbp(1)*ddss(1,1,6) + f2*kvec(1)*ddsc(1,1,6) 
                  
                  ddss(1,3,1) = rbp(2)*ddss(1,1,1) + f2*kvec(2)*ddsc(1,1,1)
                  ddss(1,3,2) = rbp(2)*ddss(1,1,2) + f2*dsc(1,1,1) + f2*kvec(2)*ddsc(1,1,2)
                  ddss(1,3,3) = rbp(2)*ddss(1,1,3) + f2*kvec(2)*ddsc(1,1,3)
                  ddss(1,3,4) = rbp(2)*ddss(1,1,4) + 2*f2*dsc(1,1,2) + f2*kvec(2)*ddsc(1,1,4)
                  ddss(1,3,5) = rbp(2)*ddss(1,1,5) + f2*dsc(1,1,3) + f2*kvec(2)*ddsc(1,1,5) 
                  ddss(1,3,6) = rbp(2)*ddss(1,1,6) + f2*kvec(2)*ddsc(1,1,6)
                  
                  ddss(1,4,1) = rbp(3)*ddss(1,1,1) + f2*kvec(3)*ddsc(1,1,1)
                  ddss(1,4,2) = rbp(3)*ddss(1,1,2) + f2*kvec(3)*ddsc(1,1,2)
                  ddss(1,4,3) = rbp(3)*ddss(1,1,3) + f2*dsc(1,1,1) + f2*kvec(3)*ddsc(1,1,3)
                  ddss(1,4,4) = rbp(3)*ddss(1,1,4) + f2*kvec(3)*ddsc(1,1,4)
                  ddss(1,4,5) = rbp(3)*ddss(1,1,5) + f2*dsc(1,1,2) + f2*kvec(3)*ddsc(1,1,5)
                  ddss(1,4,6) = rbp(3)*ddss(1,1,6) + 2*f2*dsc(1,1,3) + f2*kvec(3)*ddsc(1,1,6)
               END IF

!           *** [s|O|b] = (Pi - Bi)*[s|O|b-1i] + f2*Ni(b-1i)*[s|O|b-2i] ***
!           ***           + [s|dO|b-1i]                                 ***

                  DO lb = 2, lb_max

!             *** Increase the angular momentum component z of function b ***

                     !-------------------------------------
                     sc(1, coset(0, 0, lb)) = rbp(3)*sc(1, coset(0, 0, lb-1))+ &
                                              f2*REAL(lb-1, dp)*sc(1, coset(0, 0, lb-2))- &
                                              f2*kvec(3)*ss(1, coset(0, 0, lb-1))
                     ss(1, coset(0, 0, lb)) = rbp(3)*ss(1, coset(0, 0, lb-1))+ &
                                              f2*REAL(lb-1, dp)*ss(1, coset(0, 0, lb-2))+ &
                                              f2*kvec(3)*sc(1, coset(0, 0, lb-1))
                     !-------------------------------------
              
                     IF (PRESENT(dksc)) THEN
                        dsc(1, coset(0, 0, lb), 1) = &
                            rbp(3)*dsc(1, coset(0, 0, lb-1), 1)+ &
                            f2*REAL(lb-1, dp)*dsc(1, coset(0, 0, lb-2), 1)- &
                            f2*kvec(3)*dss(1, coset(0, 0, lb-1), 1)
                        dsc(1, coset(0, 0, lb), 2) = &
                            rbp(3)*dsc(1, coset(0, 0, lb-1), 2)+ &
                            f2*REAL(lb-1, dp)*dsc(1, coset(0, 0, lb-2), 2)- &
                            f2*kvec(3)*dss(1, coset(0, 0, lb-1), 2)
                        dsc(1, coset(0, 0, lb), 3) = &
                            rbp(3)*dsc(1, coset(0, 0, lb-1), 3)+ &
                            f2*REAL(lb-1, dp)*dsc(1, coset(0, 0, lb-2), 3)- &
                            f2*ss(1, coset(0, 0, lb-1))- &
                            f2*kvec(3)*dss(1, coset(0, 0, lb-1), 3)
                     END IF
                     IF (PRESENT(ddksc)) THEN
                        ddsc(1, coset(0, 0, lb), 1) = &
                             rbp(3)*ddsc(1, coset(0, 0, lb-1), 1)+ &
                             f2*REAL(lb-1, dp)*ddsc(1, coset(0, 0, lb-2), 1)- &
                             f2*kvec(3)*ddss(1, coset(0, 0, lb-1), 1)
                        ddsc(1, coset(0, 0, lb), 2) = &
                             rbp(3)*ddsc(1, coset(0, 0, lb-1), 2)+ &
                             f2*REAL(lb-1, dp)*ddsc(1, coset(0, 0, lb-2), 2)- &
                             f2*kvec(3)*ddss(1, coset(0, 0, lb-1), 2)
                        ddsc(1, coset(0, 0, lb), 3) = &
                             rbp(3)*ddsc(1, coset(0, 0, lb-1), 3)+ &
                             f2*REAL(lb-1, dp)*ddsc(1, coset(0, 0, lb-2), 3)- &
                             f2*dss(1, coset(0, 0, lb-1), 1)- &
                             f2*kvec(3)*ddss(1, coset(0, 0, lb-1), 3)
                        ddsc(1, coset(0, 0, lb), 4) = &
                             rbp(3)*ddsc(1, coset(0, 0, lb-1), 4)+ &
                             f2*REAL(lb-1, dp)*ddsc(1, coset(0, 0, lb-2), 4)- &
                             f2*kvec(3)*ddss(1, coset(0, 0, lb-1), 4)
                        ddsc(1, coset(0, 0, lb), 5) = &
                             rbp(3)*ddsc(1, coset(0, 0, lb-1), 5)+ &
                             f2*REAL(lb-1, dp)*ddsc(1, coset(0, 0, lb-2), 5)- &
                             f2*dss(1, coset(0, 0, lb-1), 2)- &
                             f2*kvec(3)*ddss(1, coset(0, 0, lb-1), 5)
                        ddsc(1, coset(0, 0, lb), 6) = &
                             rbp(3)*ddsc(1, coset(0, 0, lb-1), 6)+ &
                             f2*REAL(lb-1, dp)*ddsc(1, coset(0, 0, lb-2), 6)- &
                             2*f2*kvec(3)*dss(1, coset(0, 0, lb-1), 3)- &
                             f2*kvec(3)*ddss(1, coset(0, 0, lb-1), 6)
                     END IF
                     IF (PRESENT(dkss)) THEN
                        dss(1, coset(0, 0, lb), 1) = &
                            rbp(3)*dss(1, coset(0, 0, lb-1), 1)+ &
                            f2*REAL(lb-1, dp)*dss(1, coset(0, 0, lb-2), 1)+ &
                            f2*kvec(3)*dsc(1, coset(0, 0, lb-1), 1)
                        dss(1, coset(0, 0, lb), 2) = &
                            rbp(3)*dss(1, coset(0, 0, lb-1), 2)+ &
                            f2*REAL(lb-1, dp)*dss(1, coset(0, 0, lb-2), 2)+ &
                            f2*kvec(3)*dsc(1, coset(0, 0, lb-1), 2)
                        dss(1, coset(0, 0, lb), 3) = &
                            rbp(3)*dss(1, coset(0, 0, lb-1), 3)+ &
                            f2*REAL(lb-1, dp)*dss(1, coset(0, 0, lb-2), 3)+ &
                            f2*sc(1, coset(0, 0, lb-1))+ &
                            f2*kvec(3)*dsc(1, coset(0, 0, lb-1), 3)
                     END IF
                     IF (PRESENT(ddkss)) THEN
                        ddss(1, coset(0, 0, lb), 1) = &
                             rbp(3)*ddss(1, coset(0, 0, lb-1), 1)+ &
                             f2*REAL(lb-1, dp)*ddss(1, coset(0, 0, lb-2), 1)+ &
                             f2*kvec(3)*ddsc(1, coset(0, 0, lb-1), 1)
                        ddss(1, coset(0, 0, lb), 2) = &
                             rbp(3)*ddss(1, coset(0, 0, lb-1), 2)+ &
                             f2*REAL(lb-1, dp)*ddss(1, coset(0, 0, lb-2), 2)+ &
                             f2*kvec(3)*ddsc(1, coset(0, 0, lb-1), 2)
                        ddss(1, coset(0, 0, lb), 3) = &
                             rbp(3)*ddss(1, coset(0, 0, lb-1), 3)+ &
                             f2*REAL(lb-1, dp)*ddss(1, coset(0, 0, lb-2), 3)+ &
                             f2*dsc(1, coset(0, 0, lb-1), 1)+ &
                             f2*kvec(3)*ddsc(1, coset(0, 0, lb-1), 3)
                        ddss(1, coset(0, 0, lb), 4) = &
                             rbp(3)*ddss(1, coset(0, 0, lb-1), 4)+ &
                             f2*REAL(lb-1, dp)*ddss(1, coset(0, 0, lb-2), 4)+ &
                             f2*kvec(3)*ddsc(1, coset(0, 0, lb-1), 4)
                        ddss(1, coset(0, 0, lb), 5) = &
                             rbp(3)*ddss(1, coset(0, 0, lb-1), 5)+ &
                             f2*REAL(lb-1, dp)*ddss(1, coset(0, 0, lb-2), 5)+ &
                             f2*dsc(1, coset(0, 0, lb-1), 2)+ &
                             f2*kvec(3)*ddsc(1, coset(0, 0, lb-1), 5)
                        ddss(1, coset(0, 0, lb), 6) = &
                             rbp(3)*ddss(1, coset(0, 0, lb-1), 6)+ &
                             f2*REAL(lb-1, dp)*ddss(1, coset(0, 0, lb-2), 6)+ &
                             2*f2*kvec(3)*dsc(1, coset(0, 0, lb-1), 3)+ &
                             f2*kvec(3)*ddsc(1, coset(0, 0, lb-1), 6)
                     END IF

!             *** Increase the angular momentum component y of function b ***

                     bz = lb-1
                     !-------------------------------------
                     sc(1, coset(0, 1, bz)) = rbp(2)*sc(1, coset(0, 0, bz))- &
                                              f2*kvec(2)*ss(1, coset(0, 0, bz))
                     ss(1, coset(0, 1, bz)) = rbp(2)*ss(1, coset(0, 0, bz))+ &
                                              f2*kvec(2)*sc(1, coset(0, 0, bz))
                     !-------------------------------------

                     IF (PRESENT(dksc)) THEN
                        dsc(1, coset(0, 1, bz), 1) = rbp(2)*dsc(1, coset(0, 0, bz), 1)- &
                                                     f2*kvec(2)*dss(1, coset(0, 0, bz), 1)
                        dsc(1, coset(0, 1, bz), 2) = rbp(2)*dsc(1, coset(0, 0, bz), 2)- &
                                                     f2*ss(1, coset(0, 0, bz))- &
                                                     f2*kvec(2)*dss(1, coset(0, 0, bz), 2)
                        dsc(1, coset(0, 1, bz), 3) = rbp(2)*dsc(1, coset(0, 0, bz), 3)- &
                                                     f2*kvec(2)*dss(1, coset(0, 0, bz), 3)
                     END IF

                     IF (PRESENT(ddksc)) THEN
                        ddsc(1, coset(0, 1, bz), 1) = rbp(2)*ddsc(1, coset(0, 0, bz), 1)- &
                                                      f2*kvec(2)*ddss(1, coset(0, 0, bz), 1)
                        ddsc(1, coset(0, 1, bz), 2) = rbp(2)*ddsc(1, coset(0, 0, bz), 2)- &
                                                      f2*dss(1, coset(0, 0, bz), 1)- &
                                                      f2*kvec(2)*ddss(1, coset(0, 0, bz), 2)
                        ddsc(1, coset(0, 1, bz), 3) = rbp(2)*ddsc(1, coset(0, 0, bz), 3)- &
                                                      f2*kvec(2)*ddss(1, coset(0, 0, bz), 3)
                        ddsc(1, coset(0, 1, bz), 4) = rbp(2)*ddsc(1, coset(0, 0, bz), 4)- &
                                                      2*f2*dss(1, coset(0, 0, bz), 2)- &
                                                      f2*kvec(2)*ddss(1, coset(0, 0, bz), 4)
                        ddsc(1, coset(0, 1, bz), 5) = rbp(2)*ddsc(1, coset(0, 0, bz), 5)- &
                                                      f2*dss(1, coset(0, 0, bz), 3)- &
                                                      f2*kvec(2)*ddss(1, coset(0, 0, bz), 5)
                        ddsc(1, coset(0, 1, bz), 6) = rbp(2)*ddsc(1, coset(0, 0, bz), 6)- &
                                                      f2*kvec(2)*ddss(1, coset(0, 0, bz), 6)
                     END IF

                     IF (PRESENT(dkss)) THEN
                        dss(1, coset(0, 1, bz), 1) = rbp(2)*dss(1, coset(0, 0, bz), 1)+ &
                                                     f2*kvec(2)*dsc(1, coset(0, 0, bz), 1)
                        dss(1, coset(0, 1, bz), 2) = rbp(2)*dss(1, coset(0, 0, bz), 2)+ &
                                                     f2*sc(1, coset(0, 0, bz))+ &
                                                     f2*kvec(2)*dsc(1, coset(0, 0, bz), 2)
                        dss(1, coset(0, 1, bz), 3) = rbp(2)*dss(1, coset(0, 0, bz), 3)+ &
                                                     f2*kvec(2)*dsc(1, coset(0, 0, bz), 3)
                     END IF

                     IF (PRESENT(ddkss)) THEN
                        ddss(1, coset(0, 1, bz), 1) = rbp(2)*ddss(1, coset(0, 0, bz), 1)+ &
                                                      f2*kvec(2)*ddsc(1, coset(0, 0, bz), 1)
                        ddss(1, coset(0, 1, bz), 2) = rbp(2)*ddss(1, coset(0, 0, bz), 2)+ &
                                                      f2*dsc(1, coset(0, 0, bz), 1)+ &
                                                      f2*kvec(2)*ddsc(1, coset(0, 0, bz), 2)
                        ddss(1, coset(0, 1, bz), 3) = rbp(2)*ddss(1, coset(0, 0, bz), 3)+ &
                                                      f2*kvec(2)*ddsc(1, coset(0, 0, bz), 3)
                        ddss(1, coset(0, 1, bz), 4) = rbp(2)*ddss(1, coset(0, 0, bz), 4)+ &
                                                      2*f2*dsc(1, coset(0, 0, bz), 2)+ &
                                                      f2*kvec(2)*ddsc(1, coset(0, 0, bz), 4)
                        ddss(1, coset(0, 1, bz), 5) = rbp(2)*ddss(1, coset(0, 0, bz), 5)+ &
                                                      f2*dsc(1, coset(0, 0, bz), 3)+ &
                                                      f2*kvec(2)*ddsc(1, coset(0, 0, bz), 5)
                        ddss(1, coset(0, 1, bz), 6) = rbp(2)*ddss(1, coset(0, 0, bz), 6)+ &
                                                      f2*kvec(2)*ddsc(1, coset(0, 0, bz), 6)
                     END IF

                     DO by = 2, lb
                        bz = lb-by
                        !-------------------------------------
                        sc(1, coset(0, by, bz)) = rbp(2)*sc(1, coset(0, by-1, bz))+ &
                                                  f2*REAL(by-1, dp)*sc(1, coset(0, by-2, bz))- &
                                                  f2*kvec(2)*ss(1, coset(0, by-1, bz))
                        ss(1, coset(0, by, bz)) = rbp(2)*ss(1, coset(0, by-1, bz))+ &
                                                  f2*REAL(by-1, dp)*ss(1, coset(0, by-2, bz))+ &
                                                  f2*kvec(2)*sc(1, coset(0, by-1, bz))
                        !-------------------------------------
                    
                        IF (PRESENT(dksc)) THEN
                           dsc(1, coset(0, by, bz), 1) = &
                               rbp(2)*dsc(1, coset(0, by-1, bz), 1)+ &
                               f2*REAL(by-1, dp)*dsc(1, coset(0, by-2, bz), 1)- &
                               f2*kvec(2)*dss(1, coset(0, by-1, bz), 1)
                           dsc(1, coset(0, by, bz), 2) = &
                               rbp(2)*dsc(1, coset(0, by-1, bz), 2)+ &
                               f2*REAL(by-1, dp)*dsc(1, coset(0, by-2, bz), 2)- &
                               f2*ss(1, coset(0, by-1, bz))- &
                               f2*kvec(2)*dss(1, coset(0, by-1, bz), 2)
                           dsc(1, coset(0, by, bz), 3) = &
                               rbp(2)*dsc(1, coset(0, by-1, bz), 3)+ &
                               f2*REAL(by-1, dp)*dsc(1, coset(0, by-2, bz), 3)- &
                               f2*kvec(2)*dss(1, coset(0, by-1, bz), 3)
                        END IF
                        IF (PRESENT(ddksc)) THEN
                           ddsc(1, coset(0, by, bz), 1) = &
                                rbp(2)*ddsc(1, coset(0, by-1, bz), 1)+ &
                                f2*REAL(by-1, dp)*ddsc(1, coset(0, by-2, bz), 1)- &
                                f2*kvec(2)*ddss(1, coset(0, by-1, bz), 1)
                           ddsc(1, coset(0, by, bz), 2) = &
                                rbp(2)*ddsc(1, coset(0, by-1, bz), 2)+ &
                                f2*REAL(by-1, dp)*ddsc(1, coset(0, by-2, bz), 2)- &
                                f2*dss(1, coset(0, by-1, bz), 1)- &
                                f2*kvec(2)*ddss(1, coset(0, by-1, bz), 2)
                           ddsc(1, coset(0, by, bz), 3) = &
                                rbp(2)*ddsc(1, coset(0, by-1, bz), 3)+ &
                                f2*REAL(by-1, dp)*ddsc(1, coset(0, by-2, bz), 3)- &
                                f2*kvec(2)*ddss(1, coset(0, by-1, bz), 3)
                           ddsc(1, coset(0, by, bz), 4) = &
                                rbp(2)*ddsc(1, coset(0, by-1, bz), 4)+ &
                                f2*REAL(by-1, dp)*ddsc(1, coset(0, by-2, bz), 4)- &
                                2*f2*dss(1, coset(0, by-1, bz), 2)- &
                                f2*kvec(2)*ddss(1, coset(0, by-1, bz), 4)
                           ddsc(1, coset(0, by, bz), 5) = &
                                rbp(2)*ddsc(1, coset(0, by-1, bz), 5)+ &
                                f2*REAL(by-1, dp)*ddsc(1, coset(0, by-2, bz), 5)- &
                                f2*dss(1, coset(0, by-1, bz), 3)- &
                                f2*kvec(2)*ddss(1, coset(0, by-1, bz), 5)
                           ddsc(1, coset(0, by, bz), 6) = &
                                rbp(2)*ddsc(1, coset(0, by-1, bz), 6)+ &
                                f2*REAL(by-1, dp)*ddsc(1, coset(0, by-2, bz), 6)- &
                                f2*kvec(2)*ddss(1, coset(0, by-1, bz), 6)
                        END IF
                        IF (PRESENT(dkss)) THEN
                           dss(1, coset(0, by, bz), 1) = &
                               rbp(2)*dss(1, coset(0, by-1, bz), 1)+ &
                               f2*REAL(by-1, dp)*dss(1, coset(0, by-2, bz), 1)+ &
                               f2*kvec(2)*dsc(1, coset(0, by-1, bz), 1)
                           dss(1, coset(0, by, bz), 2) = &
                               rbp(2)*dss(1, coset(0, by-1, bz), 2)+ &
                               f2*REAL(by-1, dp)*dss(1, coset(0, by-2, bz), 2)+ &
                               f2*sc(1, coset(0, by-1, bz))+ &
                               f2*kvec(2)*dsc(1, coset(0, by-1, bz), 2)
                           dss(1, coset(0, by, bz), 3) = &
                               rbp(2)*dss(1, coset(0, by-1, bz), 3)+ &
                               f2*REAL(by-1, dp)*dss(1, coset(0, by-2, bz), 3)+ &
                               f2*kvec(2)*dsc(1, coset(0, by-1, bz), 3)
                        END IF
                        IF (PRESENT(ddkss)) THEN
                           ddss(1, coset(0, by, bz), 1) = &
                                rbp(2)*ddss(1, coset(0, by-1, bz), 1)+ &
                                f2*REAL(by-1, dp)*ddss(1, coset(0, by-2, bz), 1)+ &
                                f2*kvec(2)*ddsc(1, coset(0, by-1, bz), 1)
                           ddss(1, coset(0, by, bz), 2) = &
                                rbp(2)*ddss(1, coset(0, by-1, bz), 2)+ &
                                f2*REAL(by-1, dp)*ddss(1, coset(0, by-2, bz), 2)+ &
                                f2*dsc(1, coset(0, by-1, bz), 1)+ &
                                f2*kvec(2)*ddsc(1, coset(0, by-1, bz), 2)
                           ddss(1, coset(0, by, bz), 3) = &
                                rbp(2)*ddss(1, coset(0, by-1, bz), 3)+ &
                                f2*REAL(by-1, dp)*ddss(1, coset(0, by-2, bz), 3)+ &
                                f2*kvec(2)*ddsc(1, coset(0, by-1, bz), 3)
                           ddss(1, coset(0, by, bz), 4) = &
                                rbp(2)*ddss(1, coset(0, by-1, bz), 4)+ &
                                f2*REAL(by-1, dp)*ddss(1, coset(0, by-2, bz), 4)+ &
                                2*f2*dsc(1, coset(0, by-1, bz), 2)+ &
                                f2*kvec(2)*ddsc(1, coset(0, by-1, bz), 4)
                           ddss(1, coset(0, by, bz), 5) = &
                                rbp(2)*ddss(1, coset(0, by-1, bz), 5)+ &
                                f2*REAL(by-1, dp)*ddss(1, coset(0, by-2, bz), 5)+ &
                                f2*dsc(1, coset(0, by-1, bz), 3)+ &
                                f2*kvec(2)*ddsc(1, coset(0, by-1, bz), 5)
                           ddss(1, coset(0, by, bz), 6) = &
                                rbp(2)*ddss(1, coset(0, by-1, bz), 6)+ &
                                f2*REAL(by-1, dp)*ddss(1, coset(0, by-2, bz), 6)+ &
                                f2*kvec(2)*ddsc(1, coset(0, by-1, bz), 6)
                        END IF
                    
                     END DO

!             *** Increase the angular momentum component x of function b ***

                     DO by = 0, lb-1
                        bz = lb-1-by
                        !----------------------------------
                        sc(1, coset(1, by, bz)) = rbp(1)*sc(1, coset(0, by, bz))- &
                                                  f2*kvec(1)*ss(1, coset(0, by, bz))
                        ss(1, coset(1, by, bz)) = rbp(1)*ss(1, coset(0, by, bz))+ &
                                                  f2*kvec(1)*sc(1, coset(0, by, bz))
                        !----------------------------------
                                               
                        IF (PRESENT(dksc)) THEN
                           dsc(1, coset(1, by, bz), 1) = &
                               rbp(1)*dsc(1, coset(0, by, bz), 1)- &
                               f2*ss(1, coset(0, by, bz))- &
                               f2*kvec(1)*dss(1, coset(0, by, bz), 1)
                           dsc(1, coset(1, by, bz), 2) = &
                               rbp(1)*dsc(1, coset(0, by, bz), 2)- &
                               f2*kvec(1)*dss(1, coset(0, by, bz), 2)
                           dsc(1, coset(1, by, bz), 3) = &
                               rbp(1)*dsc(1, coset(0, by, bz), 3)- &
                               f2*kvec(1)*dss(1, coset(0, by, bz), 3)
                        END IF                       
                        IF (PRESENT(ddksc)) THEN
                           ddsc(1, coset(1, by, bz), 1) = &
                                rbp(1)*ddsc(1, coset(0, by, bz), 1)- &
                                2*f2*dss(1, coset(0, by, bz), 1)- &
                                f2*kvec(1)*ddss(1, coset(0, by, bz), 1)
                           ddsc(1, coset(1, by, bz), 2) = &
                                rbp(1)*ddsc(1, coset(0, by, bz), 2)- &
                                f2*dss(1, coset(0, by, bz), 2)- &
                                f2*kvec(1)*ddss(1, coset(0, by, bz), 2)
                           ddsc(1, coset(1, by, bz), 3) = &
                                rbp(1)*ddsc(1, coset(0, by, bz), 3)- &
                                f2*dss(1, coset(0, by, bz), 3)- &
                                f2*kvec(1)*ddss(1, coset(0, by, bz), 3)
                           ddsc(1, coset(1, by, bz), 4) = &
                                rbp(1)*ddsc(1, coset(0, by, bz), 4)- &
                                f2*kvec(1)*ddss(1, coset(0, by, bz), 4)
                           ddsc(1, coset(1, by, bz), 5) = &
                                rbp(1)*ddsc(1, coset(0, by, bz), 5)- &
                                f2*kvec(1)*ddss(1, coset(0, by, bz), 5)
                           ddsc(1, coset(1, by, bz), 6) = &
                                rbp(1)*ddsc(1, coset(0, by, bz), 6)- &
                                f2*kvec(1)*ddss(1, coset(0, by, bz), 6)
                        END IF                       
                        IF (PRESENT(dkss)) THEN
                           dss(1, coset(1, by, bz), 1) = &
                               rbp(1)*dss(1, coset(0, by, bz), 1)+ &
                               f2*sc(1, coset(0, by, bz))+ &
                               f2*kvec(1)*dsc(1, coset(0, by, bz), 1)
                           dss(1, coset(1, by, bz), 2) = &
                               rbp(1)*dss(1, coset(0, by, bz), 2)+ &
                               f2*kvec(1)*dsc(1, coset(0, by, bz), 2)
                           dss(1, coset(1, by, bz), 3) = &
                               rbp(1)*dss(1, coset(0, by, bz), 3)+ &
                               f2*kvec(1)*dsc(1, coset(0, by, bz), 3)
                        END IF                       
                        IF (PRESENT(ddkss)) THEN
                           ddss(1, coset(1, by, bz), 1) = &
                                rbp(1)*ddss(1, coset(0, by, bz), 1)+ &
                                2*f2*dsc(1, coset(0, by, bz), 1)+ &
                                f2*kvec(1)*ddsc(1, coset(0, by, bz), 1)
                           ddss(1, coset(1, by, bz), 2) = &
                                rbp(1)*ddss(1, coset(0, by, bz), 2)+ &
                                f2*dsc(1, coset(0, by, bz), 2)+ &
                                f2*kvec(1)*ddsc(1, coset(0, by, bz), 2)
                           ddss(1, coset(1, by, bz), 3) = &
                                rbp(1)*ddss(1, coset(0, by, bz), 3)+ &
                                f2*dsc(1, coset(0, by, bz), 3)+ &
                                f2*kvec(1)*ddsc(1, coset(0, by, bz), 3)
                           ddss(1, coset(1, by, bz), 4) = &
                                rbp(1)*ddss(1, coset(0, by, bz), 4)+ &
                                f2*kvec(1)*ddsc(1, coset(0, by, bz), 4)
                           ddss(1, coset(1, by, bz), 5) = &
                                rbp(1)*ddss(1, coset(0, by, bz), 5)+ &
                                f2*kvec(1)*ddsc(1, coset(0, by, bz), 5)
                           ddss(1, coset(1, by, bz), 6) = &
                                rbp(1)*ddss(1, coset(0, by, bz), 6)+ &
                                f2*kvec(1)*ddsc(1, coset(0, by, bz), 6)
                        END IF                       
                     END DO

                     DO bx = 2, lb
                        f3 = f2*REAL(bx-1, dp)
                        DO by = 0, lb-bx
                           bz = lb-bx-by
                           !----------------------------------
                           sc(1, coset(bx, by, bz)) = rbp(1)*sc(1, coset(bx-1, by, bz))+ &
                                                      f3*sc(1, coset(bx-2, by, bz))- &
                                                      f2*kvec(1)*ss(1, coset(bx-1, by, bz))
                           ss(1, coset(bx, by, bz)) = rbp(1)*ss(1, coset(bx-1, by, bz))+ &
                                                      f3*ss(1, coset(bx-2, by, bz))+ &
                                                      f2*kvec(1)*sc(1, coset(bx-1, by, bz))
                          !----------------------------------
                                                   
                           IF (PRESENT(dksc)) THEN
                              dsc(1, coset(bx, by, bz), 1) = &
                                  rbp(1)*dsc(1, coset(bx-1, by, bz), 1)+ &
                                  f3*dsc(1, coset(bx-2, by, bz), 1)- &
                                  f2*ss(1, coset(bx-1, by, bz))- &
                                  f2*kvec(1)*dss(1, coset(bx-1, by, bz), 1)
                              dsc(1, coset(bx, by, bz), 2) = &
                                  rbp(1)*dsc(1, coset(bx-1, by, bz), 2)+ &
                                  f3*dsc(1, coset(bx-2, by, bz), 2)- &
                                  f2*kvec(1)*dss(1, coset(bx-1, by, bz), 2)
                              dsc(1, coset(bx, by, bz), 3) = &
                                  rbp(1)*dsc(1, coset(bx-1, by, bz), 3)+ &
                                  f3*dsc(1, coset(bx-2, by, bz), 3)- &
                                  f2*kvec(1)*dss(1, coset(bx-1, by, bz), 3)
                           END IF                       
                           IF (PRESENT(ddksc)) THEN
                              ddsc(1, coset(bx, by, bz), 1) = &
                                   rbp(1)*ddsc(1, coset(bx-1, by, bz), 1)+ &
                                   f3*ddsc(1, coset(bx-2, by, bz), 1)- &
                                   2*f2*dss(1, coset(bx-1, by, bz), 1)- &
                                   f2*kvec(1)*ddss(1, coset(bx-1, by, bz), 1)
                              ddsc(1, coset(bx, by, bz), 2) = &
                                   rbp(1)*ddsc(1, coset(bx-1, by, bz), 2)+ &
                                   f3*ddsc(1, coset(bx-2, by, bz), 2)- &
                                   f2*dss(1, coset(bx-1, by, bz), 2)- &
                                   f2*kvec(1)*ddss(1, coset(bx-1, by, bz), 2)
                              ddsc(1, coset(bx, by, bz), 3) = &
                                   rbp(1)*ddsc(1, coset(bx-1, by, bz), 3)+ &
                                   f3*ddsc(1, coset(bx-2, by, bz), 3)- &
                                   f2*dss(1, coset(bx-1, by, bz), 3)- &
                                   f2*kvec(1)*ddss(1, coset(bx-1, by, bz), 3)
                              ddsc(1, coset(bx, by, bz), 4) = &
                                   rbp(1)*ddsc(1, coset(bx-1, by, bz), 4)+ &
                                   f3*ddsc(1, coset(bx-2, by, bz), 4)- &
                                   f2*kvec(1)*ddss(1, coset(bx-1, by, bz), 4)
                              ddsc(1, coset(bx, by, bz), 5) = &
                                   rbp(1)*ddsc(1, coset(bx-1, by, bz), 5)+ &
                                   f3*ddsc(1, coset(bx-2, by, bz), 5)- &
                                   f2*kvec(1)*ddss(1, coset(bx-1, by, bz), 5)
                              ddsc(1, coset(bx, by, bz), 6) = &
                                   rbp(1)*ddsc(1, coset(bx-1, by, bz), 6)+ &
                                   f3*ddsc(1, coset(bx-2, by, bz), 6)- &
                                   f2*kvec(1)*ddss(1, coset(bx-1, by, bz), 6)
                           END IF                       
                           IF (PRESENT(dkss)) THEN
                              dss(1, coset(bx, by, bz), 1) = &
                                  rbp(1)*dss(1, coset(bx-1, by, bz), 1)+ &
                                  f3*dss(1, coset(bx-2, by, bz), 1)+ &
                                  f2*sc(1, coset(bx-1, by, bz))+ &
                                  f2*kvec(1)*dsc(1, coset(bx-1, by, bz), 1)
                              dss(1, coset(bx, by, bz), 2) = &
                                  rbp(1)*dss(1, coset(bx-1, by, bz), 2)+ &
                                  f3*dss(1, coset(bx-2, by, bz), 2)+ &
                                  f2*kvec(1)*dsc(1, coset(bx-1, by, bz), 2)
                              dss(1, coset(bx, by, bz), 3) = &
                                  rbp(1)*dss(1, coset(bx-1, by, bz), 3)+ &
                                  f3*dss(1, coset(bx-2, by, bz), 3)+ &
                                  f2*kvec(1)*dsc(1, coset(bx-1, by, bz), 3)
                           END IF                       
                           IF (PRESENT(ddkss)) THEN
                              ddss(1, coset(bx, by, bz), 1) = &
                                   rbp(1)*ddss(1, coset(bx-1, by, bz), 1)+ &
                                   f3*ddss(1, coset(bx-2, by, bz), 1)+ &
                                   2*f2*dsc(1, coset(bx-1, by, bz), 1)+ &
                                   f2*kvec(1)*ddsc(1, coset(bx-1, by, bz), 1)
                              ddss(1, coset(bx, by, bz), 2) = &
                                   rbp(1)*ddss(1, coset(bx-1, by, bz), 2)+ &
                                   f3*ddss(1, coset(bx-2, by, bz), 2)+ &
                                   f2*dsc(1, coset(bx-1, by, bz), 2)+ &
                                   f2*kvec(1)*ddsc(1, coset(bx-1, by, bz), 2)
                              ddss(1, coset(bx, by, bz), 3) = &
                                   rbp(1)*ddss(1, coset(bx-1, by, bz), 3)+ &
                                   f3*ddss(1, coset(bx-2, by, bz), 3)+ &
                                   f2*dsc(1, coset(bx-1, by, bz), 3)+ &
                                   f2*kvec(1)*ddsc(1, coset(bx-1, by, bz), 3)
                              ddss(1, coset(bx, by, bz), 4) = &
                                   rbp(1)*ddss(1, coset(bx-1, by, bz), 4)+ &
                                   f3*ddss(1, coset(bx-2, by, bz), 4)+ &
                                   f2*kvec(1)*ddsc(1, coset(bx-1, by, bz), 4)
                              ddss(1, coset(bx, by, bz), 5) = &
                                   rbp(1)*ddss(1, coset(bx-1, by, bz), 5)+ &
                                   f3*ddss(1, coset(bx-2, by, bz), 5)+ &
                                   f2*kvec(1)*ddsc(1, coset(bx-1, by, bz), 5)
                              ddss(1, coset(bx, by, bz), 6) = &
                                   rbp(1)*ddss(1, coset(bx-1, by, bz), 6)+ &
                                   f3*ddss(1, coset(bx-2, by, bz), 6)+ &
                                   f2*kvec(1)*ddsc(1, coset(bx-1, by, bz), 6)
                           END IF                       
                                                   
                        END DO
                     END DO

                  END DO

               END IF

            END IF

            DO j = ncoset(lb_min-1)+1, ncoset(lb_max)
               DO i = ncoset(la_min_set-1)+1, ncoset(la_max_set)
                  cosab(na+i, nb+j) = sc(i, j)
                  sinab(na+i, nb+j) = ss(i, j)

                  IF (PRESENT(dksc) .AND. PRESENT(dkss)) THEN
                     dksc(na+i, nb+j, 1) = dsc(i, j, 1)
                     dksc(na+i, nb+j, 2) = dsc(i, j, 2)
                     dksc(na+i, nb+j, 3) = dsc(i, j, 3)

                     dkss(na+i, nb+j, 1) = dss(i, j, 1)
                     dkss(na+i, nb+j, 2) = dss(i, j, 2)
                     dkss(na+i, nb+j, 3) = dss(i, j, 3)
                  END IF

                  IF (PRESENT(ddksc) .AND. PRESENT(ddkss)) THEN
                     ddksc(na+i, nb+j, 1) = ddsc(i, j, 1)
                     ddksc(na+i, nb+j, 2) = ddsc(i, j, 2)
                     ddksc(na+i, nb+j, 3) = ddsc(i, j, 3)
                     ddksc(na+i, nb+j, 4) = ddsc(i, j, 4)
                     ddksc(na+i, nb+j, 5) = ddsc(i, j, 5)
                     ddksc(na+i, nb+j, 6) = ddsc(i, j, 6)

                     ddkss(na+i, nb+j, 1) = ddss(i, j, 1)
                     ddkss(na+i, nb+j, 2) = ddss(i, j, 2)
                     ddkss(na+i, nb+j, 3) = ddss(i, j, 3)
                     ddkss(na+i, nb+j, 4) = ddss(i, j, 4)
                     ddkss(na+i, nb+j, 5) = ddss(i, j, 5)
                     ddkss(na+i, nb+j, 6) = ddss(i, j, 6)
                  END IF

               END DO
            END DO

            IF (PRESENT(dcosab)) THEN
               la_start = 0
               lb_start = 0
            ELSE
               la_start = la_min
               lb_start = lb_min
            END IF

            DO da = 0, da_max-1
               ftz = 2.0_dp*zeta(ipgf)
               DO dax = 0, da
                  DO day = 0, da-dax
                     daz = da-dax-day
                     cda = coset(dax, day, daz)-1
                     cdax = coset(dax+1, day, daz)-1
                     cday = coset(dax, day+1, daz)-1
                     cdaz = coset(dax, day, daz+1)-1
                     !*** [da/dAi|O|b] = 2*zeta*[a+1i|O|b] - Ni(a)[a-1i|O|b] ***

                     DO la = la_start, la_max-da-1
                        DO ax = 0, la
                           fax = REAL(ax, dp)
                           DO ay = 0, la-ax
                              fay = REAL(ay, dp)
                              az = la-ax-ay
                              faz = REAL(az, dp)
                              coa = coset(ax, ay, az)
                              coamx = coset(ax-1, ay, az)
                              coamy = coset(ax, ay-1, az)
                              coamz = coset(ax, ay, az-1)
                              coapx = coset(ax+1, ay, az)
                              coapy = coset(ax, ay+1, az)
                              coapz = coset(ax, ay, az+1)
                              DO lb = lb_start, lb_max
                                 DO bx = 0, lb
                                    DO by = 0, lb-bx
                                       bz = lb-bx-by
                                       cob = coset(bx, by, bz)
                                       dscos(coa, cob, cdax) = ftz*sc(coapx, cob)-fax*sc(coamx, cob)
                                       dscos(coa, cob, cday) = ftz*sc(coapy, cob)-fay*sc(coamy, cob)
                                       dscos(coa, cob, cdaz) = ftz*sc(coapz, cob)-faz*sc(coamz, cob)
                                       dssin(coa, cob, cdax) = ftz*ss(coapx, cob)-fax*ss(coamx, cob)
                                       dssin(coa, cob, cday) = ftz*ss(coapy, cob)-fay*ss(coamy, cob)
                                       dssin(coa, cob, cdaz) = ftz*ss(coapz, cob)-faz*ss(coamz, cob)
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO

                  END DO
               END DO
            END DO

            IF (PRESENT(dcosab)) THEN
               DO k = 1, 3
                  DO j = 1, ncoset(lb_max)
                     DO i = 1, ncoset(la_max_set)
                        dcosab(na+i, nb+j, k) = dscos(i, j, k)
                        dsinab(na+i, nb+j, k) = dssin(i, j, k)
                     END DO
                  END DO
               END DO
            END IF

            nb = nb+ncoset(lb_max)

         END DO

         na = na+ncoset(la_max_set)

      END DO

   END SUBROUTINE cossin2


! **************************************************************************************************
!> \brief ...
!> \param la_max_set ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min_set ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param rac ...
!> \param rbc ...
!> \param kvec ...
!> \param cosab ...
!> \param sinab ...
!> \param dcosab ...
!> \param dsinab ...
! **************************************************************************************************
   SUBROUTINE cossin(la_max_set, npgfa, zeta, rpgfa, la_min_set, &
                     lb_max, npgfb, zetb, rpgfb, lb_min, &
                     rac, rbc, kvec, cosab, sinab, dcosab, dsinab)

      INTEGER, INTENT(IN)                                :: la_max_set, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min_set, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc, kvec
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: cosab, sinab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT), &
         OPTIONAL                                        :: dcosab, dsinab

      INTEGER :: ax, ay, az, bx, by, bz, cda, cdax, cday, cdaz, coa, coamx, coamy, coamz, coapx, &
         coapy, coapz, cob, da, da_max, dax, day, daz, i, ipgf, j, jpgf, k, la, la_max, la_min, &
         la_start, lb, lb_start, na, nb
      REAL(KIND=dp)                                      :: dab, f0, f1, f2, f3, fax, fay, faz, ftz, &
                                                            fx, fy, fz, k2, kdp, rab2, s, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rap, rbp
      REAL(KIND=dp), DIMENSION(ncoset(la_max_set), &
         ncoset(lb_max), 3)                              :: dscos, dssin
      REAL(KIND=dp), &
         DIMENSION(ncoset(la_max_set+1), ncoset(lb_max)) :: sc, ss

      rab = rbc-rac
      rab2 = SUM(rab**2)
      dab = SQRT(rab2)
      k2 = kvec(1)*kvec(1)+kvec(2)*kvec(2)+kvec(3)*kvec(3)

      IF (PRESENT(dcosab)) THEN
         da_max = 1
         la_max = la_max_set+1
         la_min = MAX(0, la_min_set-1)
         dscos = 0.0_dp
         dssin = 0.0_dp
      ELSE
         da_max = 0
         la_max = la_max_set
         la_min = la_min_set
      END IF

      ! initialize all matrix elements to zero
      IF (PRESENT(dcosab)) THEN
         na = ncoset(la_max-1)*npgfa
      ELSE
         na = ncoset(la_max)*npgfa
      END IF
      nb = ncoset(lb_max)*npgfb
      cosab(1:na, 1:nb) = 0.0_dp
      sinab(1:na, 1:nb) = 0.0_dp
      IF (PRESENT(dcosab)) THEN
         dcosab(1:na, 1:nb, :) = 0.0_dp
         dsinab(1:na, 1:nb, :) = 0.0_dp
      END IF
!   *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0
      DO ipgf = 1, npgfa

         nb = 0

         DO jpgf = 1, npgfb

            ss = 0.0_dp
            sc = 0.0_dp

!       *** Screening ***
            IF (rpgfa(ipgf)+rpgfb(jpgf) < dab) THEN
               nb = nb+ncoset(lb_max)
               CYCLE
            END IF

!       *** Calculate some prefactors ***

            zetp = 1.0_dp/(zeta(ipgf)+zetb(jpgf))

            f0 = (pi*zetp)**1.5_dp
            f1 = zetb(jpgf)*zetp
            f2 = 0.5_dp*zetp

            kdp = zetp*DOT_PRODUCT(kvec, zeta(ipgf)*rac+zetb(jpgf)*rbc)

!       *** Calculate the basic two-center cos/sin integral [s|cos/sin|s] ***

            s = f0*EXP(-zeta(ipgf)*f1*rab2)*EXP(-0.25_dp*k2*zetp)
            sc(1, 1) = s*COS(kdp)
            ss(1, 1) = s*SIN(kdp)

!       *** Recurrence steps: [s|O|s] -> [a|O|b] ***

            IF (la_max > 0) THEN

!         *** Vertical recurrence steps: [s|O|s] -> [a|O|s] ***

               rap(:) = f1*rab(:)

!         *** [p|O|s] = (Pi - Ai)*[s|O|s] +[s|dO|s]  (i = x,y,z) ***

               sc(2, 1) = rap(1)*sc(1, 1)-f2*kvec(1)*ss(1, 1)
               sc(3, 1) = rap(2)*sc(1, 1)-f2*kvec(2)*ss(1, 1)
               sc(4, 1) = rap(3)*sc(1, 1)-f2*kvec(3)*ss(1, 1)
               ss(2, 1) = rap(1)*ss(1, 1)+f2*kvec(1)*sc(1, 1)
               ss(3, 1) = rap(2)*ss(1, 1)+f2*kvec(2)*sc(1, 1)
               ss(4, 1) = rap(3)*ss(1, 1)+f2*kvec(3)*sc(1, 1)

!         *** [a|O|s] = (Pi - Ai)*[a-1i|O|s] + f2*Ni(a-1i)*[a-2i|s] ***
!         ***           + [a-1i|dO|s]                               ***

               DO la = 2, la_max

!           *** Increase the angular momentum component z of function a ***

                  sc(coset(0, 0, la), 1) = rap(3)*sc(coset(0, 0, la-1), 1)+ &
                                           f2*REAL(la-1, dp)*sc(coset(0, 0, la-2), 1)- &
                                           f2*kvec(3)*ss(coset(0, 0, la-1), 1)
                  ss(coset(0, 0, la), 1) = rap(3)*ss(coset(0, 0, la-1), 1)+ &
                                           f2*REAL(la-1, dp)*ss(coset(0, 0, la-2), 1)+ &
                                           f2*kvec(3)*sc(coset(0, 0, la-1), 1)

!           *** Increase the angular momentum component y of function a ***

                  az = la-1
                  sc(coset(0, 1, az), 1) = rap(2)*sc(coset(0, 0, az), 1)- &
                                           f2*kvec(2)*ss(coset(0, 0, az), 1)
                  ss(coset(0, 1, az), 1) = rap(2)*ss(coset(0, 0, az), 1)+ &
                                           f2*kvec(2)*sc(coset(0, 0, az), 1)

                  DO ay = 2, la
                     az = la-ay
                     sc(coset(0, ay, az), 1) = rap(2)*sc(coset(0, ay-1, az), 1)+ &
                                               f2*REAL(ay-1, dp)*sc(coset(0, ay-2, az), 1)- &
                                               f2*kvec(2)*ss(coset(0, ay-1, az), 1)
                     ss(coset(0, ay, az), 1) = rap(2)*ss(coset(0, ay-1, az), 1)+ &
                                               f2*REAL(ay-1, dp)*ss(coset(0, ay-2, az), 1)+ &
                                               f2*kvec(2)*sc(coset(0, ay-1, az), 1)
                  END DO

!           *** Increase the angular momentum component x of function a ***

                  DO ay = 0, la-1
                     az = la-1-ay
                     sc(coset(1, ay, az), 1) = rap(1)*sc(coset(0, ay, az), 1)- &
                                               f2*kvec(1)*ss(coset(0, ay, az), 1)
                     ss(coset(1, ay, az), 1) = rap(1)*ss(coset(0, ay, az), 1)+ &
                                               f2*kvec(1)*sc(coset(0, ay, az), 1)
                  END DO

                  DO ax = 2, la
                     f3 = f2*REAL(ax-1, dp)
                     DO ay = 0, la-ax
                        az = la-ax-ay
                        sc(coset(ax, ay, az), 1) = rap(1)*sc(coset(ax-1, ay, az), 1)+ &
                                                   f3*sc(coset(ax-2, ay, az), 1)- &
                                                   f2*kvec(1)*ss(coset(ax-1, ay, az), 1)
                        ss(coset(ax, ay, az), 1) = rap(1)*ss(coset(ax-1, ay, az), 1)+ &
                                                   f3*ss(coset(ax-2, ay, az), 1)+ &
                                                   f2*kvec(1)*sc(coset(ax-1, ay, az), 1)
                     END DO
                  END DO

               END DO

!         *** Recurrence steps: [a|O|s] -> [a|O|b] ***

               IF (lb_max > 0) THEN

                  DO j = 2, ncoset(lb_max)
                     DO i = 1, ncoset(la_max)
                        sc(i, j) = 0.0_dp
                        ss(i, j) = 0.0_dp
                     END DO
                  END DO

!           *** Horizontal recurrence steps ***

                  rbp(:) = rap(:)-rab(:)

!           *** [a|O|p] = [a+1i|O|s] - (Bi - Ai)*[a|O|s] ***

                  IF (lb_max == 1) THEN
                     la_start = la_min
                  ELSE
                     la_start = MAX(0, la_min-1)
                  END IF

                  DO la = la_start, la_max-1
                     DO ax = 0, la
                        DO ay = 0, la-ax
                           az = la-ax-ay
                           sc(coset(ax, ay, az), 2) = sc(coset(ax+1, ay, az), 1)- &
                                                      rab(1)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 3) = sc(coset(ax, ay+1, az), 1)- &
                                                      rab(2)*sc(coset(ax, ay, az), 1)
                           sc(coset(ax, ay, az), 4) = sc(coset(ax, ay, az+1), 1)- &
                                                      rab(3)*sc(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = ss(coset(ax+1, ay, az), 1)- &
                                                      rab(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = ss(coset(ax, ay+1, az), 1)- &
                                                      rab(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = ss(coset(ax, ay, az+1), 1)- &
                                                      rab(3)*ss(coset(ax, ay, az), 1)
                        END DO
                     END DO
                  END DO

!           *** Vertical recurrence step ***

!           *** [a|O|p] = (Pi - Bi)*[a|O|s] + f2*Ni(a)*[a-1i|O|s] ***
!           ***           + [a|dO|s]                              ***

                  DO ax = 0, la_max
                     fx = f2*REAL(ax, dp)
                     DO ay = 0, la_max-ax
                        fy = f2*REAL(ay, dp)
                        az = la_max-ax-ay
                        fz = f2*REAL(az, dp)
                        IF (ax == 0) THEN
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)
                        ELSE
                           sc(coset(ax, ay, az), 2) = rbp(1)*sc(coset(ax, ay, az), 1)+ &
                                                      fx*sc(coset(ax-1, ay, az), 1)- &
                                                      f2*kvec(1)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 2) = rbp(1)*ss(coset(ax, ay, az), 1)+ &
                                                      fx*ss(coset(ax-1, ay, az), 1)+ &
                                                      f2*kvec(1)*sc(coset(ax, ay, az), 1)
                        END IF
                        IF (ay == 0) THEN
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                        ELSE
                           sc(coset(ax, ay, az), 3) = rbp(2)*sc(coset(ax, ay, az), 1)+ &
                                                      fy*sc(coset(ax, ay-1, az), 1)- &
                                                      f2*kvec(2)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 3) = rbp(2)*ss(coset(ax, ay, az), 1)+ &
                                                      fy*ss(coset(ax, ay-1, az), 1)+ &
                                                      f2*kvec(2)*sc(coset(ax, ay, az), 1)
                        END IF
                        IF (az == 0) THEN
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1)- &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1)+ &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)
                        ELSE
                           sc(coset(ax, ay, az), 4) = rbp(3)*sc(coset(ax, ay, az), 1)+ &
                                                      fz*sc(coset(ax, ay, az-1), 1)- &
                                                      f2*kvec(3)*ss(coset(ax, ay, az), 1)
                           ss(coset(ax, ay, az), 4) = rbp(3)*ss(coset(ax, ay, az), 1)+ &
                                                      fz*ss(coset(ax, ay, az-1), 1)+ &
                                                      f2*kvec(3)*sc(coset(ax, ay, az), 1)
                        END IF
                     END DO
                  END DO

!           *** Recurrence steps: [a|O|p] -> [a|O|b] ***

                  DO lb = 2, lb_max

!             *** Horizontal recurrence steps ***

!             *** [a|O|b] = [a+1i|O|b-1i] - (Bi - Ai)*[a|O|b-1i] ***

                     IF (lb == lb_max) THEN
                        la_start = la_min
                     ELSE
                        la_start = MAX(0, la_min-1)
                     END IF

                     DO la = la_start, la_max-1
                        DO ax = 0, la
                           DO ay = 0, la-ax
                              az = la-ax-ay

!                   *** Shift of angular momentum component z from a to b ***

                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 sc(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                 rab(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 ss(coset(ax, ay, az+1), coset(0, 0, lb-1))- &
                                 rab(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))

!                   *** Shift of angular momentum component y from a to b ***

                              DO by = 1, lb
                                 bz = lb-by
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    sc(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                    rab(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    ss(coset(ax, ay+1, az), coset(0, by-1, bz))- &
                                    rab(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))
                              END DO

!                   *** Shift of angular momentum component x from a to b ***

                              DO bx = 1, lb
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       sc(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                       rab(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       ss(coset(ax+1, ay, az), coset(bx-1, by, bz))- &
                                       rab(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))
                                 END DO
                              END DO

                           END DO
                        END DO
                     END DO

!             *** Vertical recurrence step ***

!             *** [a|O|b] = (Pi - Bi)*[a|O|b-1i] + f2*Ni(a)*[a-1i|O|b-1i] + ***
!             ***           f2*Ni(b-1i)*[a|O|b-2i] + [a|dO|b-1i]            ***

                     DO ax = 0, la_max
                        fx = f2*REAL(ax, dp)
                        DO ay = 0, la_max-ax
                           fy = f2*REAL(ay, dp)
                           az = la_max-ax-ay
                           fz = f2*REAL(az, dp)

!                 *** Increase the angular momentum component z of function b ***

                           f3 = f2*REAL(lb-1, dp)

                           IF (az == 0) THEN
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                           ELSE
                              sc(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 fz*sc(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                 f3*sc(coset(ax, ay, az), coset(0, 0, lb-2))- &
                                 f2*kvec(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))
                              ss(coset(ax, ay, az), coset(0, 0, lb)) = &
                                 rbp(3)*ss(coset(ax, ay, az), coset(0, 0, lb-1))+ &
                                 fz*ss(coset(ax, ay, az-1), coset(0, 0, lb-1))+ &
                                 f3*ss(coset(ax, ay, az), coset(0, 0, lb-2))+ &
                                 f2*kvec(3)*sc(coset(ax, ay, az), coset(0, 0, lb-1))
                           END IF

!                 *** Increase the angular momentum component y of function b ***

                           IF (ay == 0) THEN
                              bz = lb-1
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz))- &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              DO by = 2, lb
                                 bz = lb-by
                                 f3 = f2*REAL(by-1, dp)
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    f3*sc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    f3*ss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                              END DO
                           ELSE
                              bz = lb-1
                              sc(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*sc(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 fy*sc(coset(ax, ay-1, az), coset(0, 0, bz))- &
                                 f2*kvec(2)*ss(coset(ax, ay, az), coset(0, 0, bz))
                              ss(coset(ax, ay, az), coset(0, 1, bz)) = &
                                 rbp(2)*ss(coset(ax, ay, az), coset(0, 0, bz))+ &
                                 fy*ss(coset(ax, ay-1, az), coset(0, 0, bz))+ &
                                 f2*kvec(2)*sc(coset(ax, ay, az), coset(0, 0, bz))
                              DO by = 2, lb
                                 bz = lb-by
                                 f3 = f2*REAL(by-1, dp)
                                 sc(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    fy*sc(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                    f3*sc(coset(ax, ay, az), coset(0, by-2, bz))- &
                                    f2*kvec(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))
                                 ss(coset(ax, ay, az), coset(0, by, bz)) = &
                                    rbp(2)*ss(coset(ax, ay, az), coset(0, by-1, bz))+ &
                                    fy*ss(coset(ax, ay-1, az), coset(0, by-1, bz))+ &
                                    f3*ss(coset(ax, ay, az), coset(0, by-2, bz))+ &
                                    f2*kvec(2)*sc(coset(ax, ay, az), coset(0, by-1, bz))
                              END DO
                           END IF

!                 *** Increase the angular momentum component x of function b ***

                           IF (ax == 0) THEN
                              DO by = 0, lb-1
                                 bz = lb-1-by
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz))- &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz))+ &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                              END DO
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx-1, dp)
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       f3*sc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       f3*ss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                 END DO
                              END DO
                           ELSE
                              DO by = 0, lb-1
                                 bz = lb-1-by
                                 sc(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*sc(coset(ax, ay, az), coset(0, by, bz))+ &
                                    fx*sc(coset(ax-1, ay, az), coset(0, by, bz))- &
                                    f2*kvec(1)*ss(coset(ax, ay, az), coset(0, by, bz))
                                 ss(coset(ax, ay, az), coset(1, by, bz)) = &
                                    rbp(1)*ss(coset(ax, ay, az), coset(0, by, bz))+ &
                                    fx*ss(coset(ax-1, ay, az), coset(0, by, bz))+ &
                                    f2*kvec(1)*sc(coset(ax, ay, az), coset(0, by, bz))
                              END DO
                              DO bx = 2, lb
                                 f3 = f2*REAL(bx-1, dp)
                                 DO by = 0, lb-bx
                                    bz = lb-bx-by
                                    sc(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*sc(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       fx*sc(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                       f3*sc(coset(ax, ay, az), coset(bx-2, by, bz))- &
                                       f2*kvec(1)*ss(coset(ax, ay, az), coset(bx-1, by, bz))
                                    ss(coset(ax, ay, az), coset(bx, by, bz)) = &
                                       rbp(1)*ss(coset(ax, ay, az), &
                                                 coset(bx-1, by, bz))+ &
                                       fx*ss(coset(ax-1, ay, az), coset(bx-1, by, bz))+ &
                                       f3*ss(coset(ax, ay, az), coset(bx-2, by, bz))+ &
                                       f2*kvec(1)*sc(coset(ax, ay, az), coset(bx-1, by, bz))
                                 END DO
                              END DO
                           END IF

                        END DO
                     END DO

                  END DO

               END IF

            ELSE

               IF (lb_max > 0) THEN

!           *** Vertical recurrence steps: [s|O|s] -> [s|O|b] ***

                  rbp(:) = (f1-1.0_dp)*rab(:)

!           *** [s|O|p] = (Pi - Bi)*[s|O|s] + [s|dO|s] ***

                  sc(1, 2) = rbp(1)*sc(1, 1)-f2*kvec(1)*ss(1, 1)
                  sc(1, 3) = rbp(2)*sc(1, 1)-f2*kvec(2)*ss(1, 1)
                  sc(1, 4) = rbp(3)*sc(1, 1)-f2*kvec(3)*ss(1, 1)
                  ss(1, 2) = rbp(1)*ss(1, 1)+f2*kvec(1)*sc(1, 1)
                  ss(1, 3) = rbp(2)*ss(1, 1)+f2*kvec(2)*sc(1, 1)
                  ss(1, 4) = rbp(3)*ss(1, 1)+f2*kvec(3)*sc(1, 1)

!           *** [s|O|b] = (Pi - Bi)*[s|O|b-1i] + f2*Ni(b-1i)*[s|O|b-2i] ***
!           ***           + [s|dO|b-1i]                                 ***

                  DO lb = 2, lb_max

!             *** Increase the angular momentum component z of function b ***

                     sc(1, coset(0, 0, lb)) = rbp(3)*sc(1, coset(0, 0, lb-1))+ &
                                              f2*REAL(lb-1, dp)*sc(1, coset(0, 0, lb-2))- &
                                              f2*kvec(3)*ss(1, coset(0, 0, lb-1))
                     ss(1, coset(0, 0, lb)) = rbp(3)*ss(1, coset(0, 0, lb-1))+ &
                                              f2*REAL(lb-1, dp)*ss(1, coset(0, 0, lb-2))+ &
                                              f2*kvec(3)*sc(1, coset(0, 0, lb-1))

!             *** Increase the angular momentum component y of function b ***

                     bz = lb-1
                     sc(1, coset(0, 1, bz)) = rbp(2)*sc(1, coset(0, 0, bz))- &
                                              f2*kvec(2)*ss(1, coset(0, 0, bz))
                     ss(1, coset(0, 1, bz)) = rbp(2)*ss(1, coset(0, 0, bz))+ &
                                              f2*kvec(2)*sc(1, coset(0, 0, bz))

                     DO by = 2, lb
                        bz = lb-by
                        sc(1, coset(0, by, bz)) = rbp(2)*sc(1, coset(0, by-1, bz))+ &
                                                  f2*REAL(by-1, dp)*sc(1, coset(0, by-2, bz))- &
                                                  f2*kvec(2)*ss(1, coset(0, by-1, bz))
                        ss(1, coset(0, by, bz)) = rbp(2)*ss(1, coset(0, by-1, bz))+ &
                                                  f2*REAL(by-1, dp)*ss(1, coset(0, by-2, bz))+ &
                                                  f2*kvec(2)*sc(1, coset(0, by-1, bz))
                     END DO

!             *** Increase the angular momentum component x of function b ***

                     DO by = 0, lb-1
                        bz = lb-1-by
                        sc(1, coset(1, by, bz)) = rbp(1)*sc(1, coset(0, by, bz))- &
                                                  f2*kvec(1)*ss(1, coset(0, by, bz))
                        ss(1, coset(1, by, bz)) = rbp(1)*ss(1, coset(0, by, bz))+ &
                                                  f2*kvec(1)*sc(1, coset(0, by, bz))
                     END DO

                     DO bx = 2, lb
                        f3 = f2*REAL(bx-1, dp)
                        DO by = 0, lb-bx
                           bz = lb-bx-by
                           sc(1, coset(bx, by, bz)) = rbp(1)*sc(1, coset(bx-1, by, bz))+ &
                                                      f3*sc(1, coset(bx-2, by, bz))- &
                                                      f2*kvec(1)*ss(1, coset(bx-1, by, bz))
                           ss(1, coset(bx, by, bz)) = rbp(1)*ss(1, coset(bx-1, by, bz))+ &
                                                      f3*ss(1, coset(bx-2, by, bz))+ &
                                                      f2*kvec(1)*sc(1, coset(bx-1, by, bz))
                        END DO
                     END DO

                  END DO

               END IF

            END IF

            DO j = ncoset(lb_min-1)+1, ncoset(lb_max)
               DO i = ncoset(la_min_set-1)+1, ncoset(la_max_set)
                  cosab(na+i, nb+j) = sc(i, j)
                  sinab(na+i, nb+j) = ss(i, j)
               END DO
            END DO

            IF (PRESENT(dcosab)) THEN
               la_start = 0
               lb_start = 0
            ELSE
               la_start = la_min
               lb_start = lb_min
            END IF

            DO da = 0, da_max-1
               ftz = 2.0_dp*zeta(ipgf)
               DO dax = 0, da
                  DO day = 0, da-dax
                     daz = da-dax-day
                     cda = coset(dax, day, daz)-1
                     cdax = coset(dax+1, day, daz)-1
                     cday = coset(dax, day+1, daz)-1
                     cdaz = coset(dax, day, daz+1)-1
                     !*** [da/dAi|O|b] = 2*zeta*[a+1i|O|b] - Ni(a)[a-1i|O|b] ***

                     DO la = la_start, la_max-da-1
                        DO ax = 0, la
                           fax = REAL(ax, dp)
                           DO ay = 0, la-ax
                              fay = REAL(ay, dp)
                              az = la-ax-ay
                              faz = REAL(az, dp)
                              coa = coset(ax, ay, az)
                              coamx = coset(ax-1, ay, az)
                              coamy = coset(ax, ay-1, az)
                              coamz = coset(ax, ay, az-1)
                              coapx = coset(ax+1, ay, az)
                              coapy = coset(ax, ay+1, az)
                              coapz = coset(ax, ay, az+1)
                              DO lb = lb_start, lb_max
                                 DO bx = 0, lb
                                    DO by = 0, lb-bx
                                       bz = lb-bx-by
                                       cob = coset(bx, by, bz)
                                       dscos(coa, cob, cdax) = ftz*sc(coapx, cob)-fax*sc(coamx, cob)
                                       dscos(coa, cob, cday) = ftz*sc(coapy, cob)-fay*sc(coamy, cob)
                                       dscos(coa, cob, cdaz) = ftz*sc(coapz, cob)-faz*sc(coamz, cob)
                                       dssin(coa, cob, cdax) = ftz*ss(coapx, cob)-fax*ss(coamx, cob)
                                       dssin(coa, cob, cday) = ftz*ss(coapy, cob)-fay*ss(coamy, cob)
                                       dssin(coa, cob, cdaz) = ftz*ss(coapz, cob)-faz*ss(coamz, cob)
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO

                  END DO
               END DO
            END DO

            IF (PRESENT(dcosab)) THEN
               DO k = 1, 3
                  DO j = 1, ncoset(lb_max)
                     DO i = 1, ncoset(la_max_set)
                        dcosab(na+i, nb+j, k) = dscos(i, j, k)
                        dsinab(na+i, nb+j, k) = dssin(i, j, k)
                     END DO
                  END DO
               END DO
            END IF

            nb = nb+ncoset(lb_max)

         END DO

         na = na+ncoset(la_max_set)

      END DO

   END SUBROUTINE cossin

! **************************************************************************************************
!> \brief ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lc_max ...
!> \param rac ...
!> \param rbc ...
!> \param mab ...
! **************************************************************************************************
   SUBROUTINE moment(la_max, npgfa, zeta, rpgfa, la_min, &
                     lb_max, npgfb, zetb, rpgfb, &
                     lc_max, rac, rbc, mab)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lc_max
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: mab

      INTEGER                                            :: ax, ay, az, bx, by, bz, i, ipgf, j, &
                                                            jpgf, k, l, l1, l2, la, la_start, lb, &
                                                            lx, lx1, ly, ly1, lz, lz1, na, nb, ni
      REAL(KIND=dp)                                      :: dab, f0, f1, f2, f2x, f2y, f2z, f3, fx, &
                                                            fy, fz, rab2, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: rab, rap, rbp, rpc
      REAL(KIND=dp), DIMENSION(ncoset(la_max), ncoset(&
         lb_max), ncoset(lc_max))                        :: s

      rab = rbc-rac
      rab2 = SUM(rab**2)
      dab = SQRT(rab2)

!   *** Loop over all pairs of primitive Gaussian-type functions ***

      na = 0

      DO ipgf = 1, npgfa

         nb = 0

         DO jpgf = 1, npgfb

            s = 0.0_dp
!       *** Screening ***

            IF (rpgfa(ipgf)+rpgfb(jpgf) < dab) THEN
               DO k = 1, ncoset(lc_max)-1
                  DO j = nb+1, nb+ncoset(lb_max)
                     DO i = na+1, na+ncoset(la_max)
                        mab(i, j, k) = 0.0_dp
                     END DO
                  END DO
               END DO
               nb = nb+ncoset(lb_max)
               CYCLE
            END IF

!       *** Calculate some prefactors ***

            zetp = 1.0_dp/(zeta(ipgf)+zetb(jpgf))

            f0 = (pi*zetp)**1.5_dp
            f1 = zetb(jpgf)*zetp
            f2 = 0.5_dp*zetp

!       *** Calculate the basic two-center moment integral [s|M|s] ***

            rpc = zetp*(zeta(ipgf)*rac+zetb(jpgf)*rbc)
            s(1, 1, 1) = f0*EXP(-zeta(ipgf)*f1*rab2)
            DO l = 2, ncoset(lc_max)
               lx = indco(1, l)
               ly = indco(2, l)
               lz = indco(3, l)
               l2 = 0
               IF (lz > 0) THEN
                  l1 = coset(lx, ly, lz-1)
                  IF (lz > 1) l2 = coset(lx, ly, lz-2)
                  ni = lz-1
                  i = 3
               ELSE IF (ly > 0) THEN
                  l1 = coset(lx, ly-1, lz)
                  IF (ly > 1) l2 = coset(lx, ly-2, lz)
                  ni = ly-1
                  i = 2
               ELSE IF (lx > 0) THEN
                  l1 = coset(lx-1, ly, lz)
                  IF (lx > 1) l2 = coset(lx-2, ly, lz)
                  ni = lx-1
                  i = 1
               END IF
               s(1, 1, l) = rpc(i)*s(1, 1, l1)
               IF (l2 > 0) s(1, 1, l) = s(1, 1, l)+f2*REAL(ni, dp)*s(1, 1, l2)
            END DO

!       *** Recurrence steps: [s|M|s] -> [a|M|b] ***

            DO l = 1, ncoset(lc_max)

               lx = indco(1, l)
               ly = indco(2, l)
               lz = indco(3, l)
               IF (lx > 0) THEN
                  lx1 = coset(lx-1, ly, lz)
               ELSE
                  lx1 = -1
               END IF
               IF (ly > 0) THEN
                  ly1 = coset(lx, ly-1, lz)
               ELSE
                  ly1 = -1
               END IF
               IF (lz > 0) THEN
                  lz1 = coset(lx, ly, lz-1)
               ELSE
                  lz1 = -1
               END IF
               f2x = f2*REAL(lx, dp)
               f2y = f2*REAL(ly, dp)
               f2z = f2*REAL(lz, dp)

               IF (la_max > 0) THEN

!           *** Vertical recurrence steps: [s|M|s] -> [a|M|s] ***

                  rap(:) = f1*rab(:)

!           *** [p|M|s] = (Pi - Ai)*[s|M|s] + f2*Ni(m-1i)[s|M-1i|s] ***

                  s(2, 1, l) = rap(1)*s(1, 1, l)
                  s(3, 1, l) = rap(2)*s(1, 1, l)
                  s(4, 1, l) = rap(3)*s(1, 1, l)
                  IF (lx1 > 0) s(2, 1, l) = s(2, 1, l)+f2x*s(1, 1, lx1)
                  IF (ly1 > 0) s(3, 1, l) = s(3, 1, l)+f2y*s(1, 1, ly1)
                  IF (lz1 > 0) s(4, 1, l) = s(4, 1, l)+f2z*s(1, 1, lz1)

!           *** [a|M|s] = (Pi - Ai)*[a-1i|M|s] + f2*Ni(a-1i)*[a-2i|M|s] ***
!           ***           + f2*Ni(m-1i)*[a-1i|M-1i|s]                   ***

                  DO la = 2, la_max

!             *** Increase the angular momentum component z of function a ***

                     s(coset(0, 0, la), 1, l) = rap(3)*s(coset(0, 0, la-1), 1, l)+ &
                                                f2*REAL(la-1, dp)*s(coset(0, 0, la-2), 1, l)
                     IF (lz1 > 0) s(coset(0, 0, la), 1, l) = s(coset(0, 0, la), 1, l)+ &
                                                             f2z*s(coset(0, 0, la-1), 1, lz1)

!             *** Increase the angular momentum component y of function a ***

                     az = la-1
                     s(coset(0, 1, az), 1, l) = rap(2)*s(coset(0, 0, az), 1, l)
                     IF (ly1 > 0) s(coset(0, 1, az), 1, l) = s(coset(0, 1, az), 1, l)+ &
                                                             f2y*s(coset(0, 0, az), 1, ly1)

                     DO ay = 2, la
                        az = la-ay
                        s(coset(0, ay, az), 1, l) = rap(2)*s(coset(0, ay-1, az), 1, l)+ &
                                                    f2*REAL(ay-1, dp)*s(coset(0, ay-2, az), 1, l)
                        IF (ly1 > 0) s(coset(0, ay, az), 1, l) = s(coset(0, ay, az), 1, l)+ &
                                                                 f2y*s(coset(0, ay-1, az), 1, ly1)
                     END DO

!             *** Increase the angular momentum component x of function a ***

                     DO ay = 0, la-1
                        az = la-1-ay
                        s(coset(1, ay, az), 1, l) = rap(1)*s(coset(0, ay, az), 1, l)
                        IF (lx1 > 0) s(coset(1, ay, az), 1, l) = s(coset(1, ay, az), 1, l)+ &
                                                                 f2x*s(coset(0, ay, az), 1, lx1)
                     END DO

                     DO ax = 2, la
                        f3 = f2*REAL(ax-1, dp)
                        DO ay = 0, la-ax
                           az = la-ax-ay
                           s(coset(ax, ay, az), 1, l) = rap(1)*s(coset(ax-1, ay, az), 1, l)+ &
                                                        f3*s(coset(ax-2, ay, az), 1, l)
                           IF (lx1 > 0) s(coset(ax, ay, az), 1, l) = s(coset(ax, ay, az), 1, l)+ &
                                                                     f2x*s(coset(ax-1, ay, az), 1, lx1)
                        END DO
                     END DO

                  END DO

!           *** Recurrence steps: [a|M|s] -> [a|M|b] ***

                  IF (lb_max > 0) THEN

                     DO j = 2, ncoset(lb_max)
                        DO i = 1, ncoset(la_max)
                           s(i, j, l) = 0.0_dp
                        END DO
                     END DO

!             *** Horizontal recurrence steps ***

                     rbp(:) = rap(:)-rab(:)

!             *** [a|M|p] = [a+1i|M|s] - (Bi - Ai)*[a|M|s] ***

                     IF (lb_max == 1) THEN
                        la_start = la_min
                     ELSE
                        la_start = MAX(0, la_min-1)
                     END IF

                     DO la = la_start, la_max-1
                        DO ax = 0, la
                           DO ay = 0, la-ax
                              az = la-ax-ay
                              s(coset(ax, ay, az), 2, l) = s(coset(ax+1, ay, az), 1, l)- &
                                                           rab(1)*s(coset(ax, ay, az), 1, l)
                              s(coset(ax, ay, az), 3, l) = s(coset(ax, ay+1, az), 1, l)- &
                                                           rab(2)*s(coset(ax, ay, az), 1, l)
                              s(coset(ax, ay, az), 4, l) = s(coset(ax, ay, az+1), 1, l)- &
                                                           rab(3)*s(coset(ax, ay, az), 1, l)
                           END DO
                        END DO
                     END DO

!             *** Vertical recurrence step ***

!             *** [a|M|p] = (Pi - Bi)*[a|M|s] + f2*Ni(a)*[a-1i|M|s] ***
!             ***           + f2*Ni(m)*[a|M-1i|s]                   ***

                     DO ax = 0, la_max
                        fx = f2*REAL(ax, dp)
                        DO ay = 0, la_max-ax
                           fy = f2*REAL(ay, dp)
                           az = la_max-ax-ay
                           fz = f2*REAL(az, dp)
                           IF (ax == 0) THEN
                              s(coset(ax, ay, az), 2, l) = rbp(1)*s(coset(ax, ay, az), 1, l)
                           ELSE
                              s(coset(ax, ay, az), 2, l) = rbp(1)*s(coset(ax, ay, az), 1, l)+ &
                                                           fx*s(coset(ax-1, ay, az), 1, l)
                           END IF
                           IF (lx1 > 0) s(coset(ax, ay, az), 2, l) = s(coset(ax, ay, az), 2, l)+ &
                                                                     f2x*s(coset(ax, ay, az), 1, lx1)
                           IF (ay == 0) THEN
                              s(coset(ax, ay, az), 3, l) = rbp(2)*s(coset(ax, ay, az), 1, l)
                           ELSE
                              s(coset(ax, ay, az), 3, l) = rbp(2)*s(coset(ax, ay, az), 1, l)+ &
                                                           fy*s(coset(ax, ay-1, az), 1, l)
                           END IF
                           IF (ly1 > 0) s(coset(ax, ay, az), 3, l) = s(coset(ax, ay, az), 3, l)+ &
                                                                     f2y*s(coset(ax, ay, az), 1, ly1)
                           IF (az == 0) THEN
                              s(coset(ax, ay, az), 4, l) = rbp(3)*s(coset(ax, ay, az), 1, l)
                           ELSE
                              s(coset(ax, ay, az), 4, l) = rbp(3)*s(coset(ax, ay, az), 1, l)+ &
                                                           fz*s(coset(ax, ay, az-1), 1, l)
                           END IF
                           IF (lz1 > 0) s(coset(ax, ay, az), 4, l) = s(coset(ax, ay, az), 4, l)+ &
                                                                     f2z*s(coset(ax, ay, az), 1, lz1)
                        END DO
                     END DO

!             *** Recurrence steps: [a|M|p] -> [a|M|b] ***

                     DO lb = 2, lb_max

!               *** Horizontal recurrence steps ***

!               *** [a|M|b] = [a+1i|M|b-1i] - (Bi - Ai)*[a|M|b-1i] ***

                        IF (lb == lb_max) THEN
                           la_start = la_min
                        ELSE
                           la_start = MAX(0, la_min-1)
                        END IF

                        DO la = la_start, la_max-1
                           DO ax = 0, la
                              DO ay = 0, la-ax
                                 az = la-ax-ay

!                     *** Shift of angular momentum component z from a to b ***

                                 s(coset(ax, ay, az), coset(0, 0, lb), l) = &
                                    s(coset(ax, ay, az+1), coset(0, 0, lb-1), l)- &
                                    rab(3)*s(coset(ax, ay, az), coset(0, 0, lb-1), l)

!                     *** Shift of angular momentum component y from a to b ***

                                 DO by = 1, lb
                                    bz = lb-by
                                    s(coset(ax, ay, az), coset(0, by, bz), l) = &
                                       s(coset(ax, ay+1, az), coset(0, by-1, bz), l)- &
                                       rab(2)*s(coset(ax, ay, az), coset(0, by-1, bz), l)
                                 END DO

!                     *** Shift of angular momentum component x from a to b ***

                                 DO bx = 1, lb
                                    DO by = 0, lb-bx
                                       bz = lb-bx-by
                                       s(coset(ax, ay, az), coset(bx, by, bz), l) = &
                                          s(coset(ax+1, ay, az), coset(bx-1, by, bz), l)- &
                                          rab(1)*s(coset(ax, ay, az), coset(bx-1, by, bz), l)
                                    END DO
                                 END DO

                              END DO
                           END DO
                        END DO

!               *** Vertical recurrence step ***

!               *** [a|M|b] = (Pi - Bi)*[a|M|b-1i] + f2*Ni(a)*[a-1i|M|b-1i] + ***
!               ***           f2*Ni(b-1i)*[a|M|b-2i] + f2*Ni(m)[a|M-1i|b-1i]  ***

                        DO ax = 0, la_max
                           fx = f2*REAL(ax, dp)
                           DO ay = 0, la_max-ax
                              fy = f2*REAL(ay, dp)
                              az = la_max-ax-ay
                              fz = f2*REAL(az, dp)

!                   *** Shift of angular momentum component z from a to b ***

                              f3 = f2*REAL(lb-1, dp)

                              IF (az == 0) THEN
                                 s(coset(ax, ay, az), coset(0, 0, lb), l) = &
                                    rbp(3)*s(coset(ax, ay, az), coset(0, 0, lb-1), l)+ &
                                    f3*s(coset(ax, ay, az), coset(0, 0, lb-2), l)
                              ELSE
                                 s(coset(ax, ay, az), coset(0, 0, lb), l) = &
                                    rbp(3)*s(coset(ax, ay, az), coset(0, 0, lb-1), l)+ &
                                    fz*s(coset(ax, ay, az-1), coset(0, 0, lb-1), l)+ &
                                    f3*s(coset(ax, ay, az), coset(0, 0, lb-2), l)
                              END IF
                              IF (lz1 > 0) s(coset(ax, ay, az), coset(0, 0, lb), l) = &
                                 s(coset(ax, ay, az), coset(0, 0, lb), l)+ &
                                 f2z*s(coset(ax, ay, az), coset(0, 0, lb-1), lz1)

!                   *** Shift of angular momentum component y from a to b ***

                              IF (ay == 0) THEN
                                 bz = lb-1
                                 s(coset(ax, ay, az), coset(0, 1, bz), l) = &
                                    rbp(2)*s(coset(ax, ay, az), coset(0, 0, bz), l)
                                 IF (ly1 > 0) s(coset(ax, ay, az), coset(0, 1, bz), l) = &
                                    s(coset(ax, ay, az), coset(0, 1, bz), l)+ &
                                    f2y*s(coset(ax, ay, az), coset(0, 0, bz), ly1)
                                 DO by = 2, lb
                                    bz = lb-by
                                    f3 = f2*REAL(by-1, dp)
                                    s(coset(ax, ay, az), coset(0, by, bz), l) = &
                                       rbp(2)*s(coset(ax, ay, az), coset(0, by-1, bz), l)+ &
                                       f3*s(coset(ax, ay, az), coset(0, by-2, bz), l)
                                    IF (ly1 > 0) s(coset(ax, ay, az), coset(0, by, bz), l) = &
                                       s(coset(ax, ay, az), coset(0, by, bz), l)+ &
                                       f2y*s(coset(ax, ay, az), coset(0, by-1, bz), ly1)
                                 END DO
                              ELSE
                                 bz = lb-1
                                 s(coset(ax, ay, az), coset(0, 1, bz), l) = &
                                    rbp(2)*s(coset(ax, ay, az), coset(0, 0, bz), l)+ &
                                    fy*s(coset(ax, ay-1, az), coset(0, 0, bz), l)
                                 IF (ly1 > 0) s(coset(ax, ay, az), coset(0, 1, bz), l) = &
                                    s(coset(ax, ay, az), coset(0, 1, bz), l)+ &
                                    f2y*s(coset(ax, ay, az), coset(0, 0, bz), ly1)
                                 DO by = 2, lb
                                    bz = lb-by
                                    f3 = f2*REAL(by-1, dp)
                                    s(coset(ax, ay, az), coset(0, by, bz), l) = &
                                       rbp(2)*s(coset(ax, ay, az), coset(0, by-1, bz), l)+ &
                                       fy*s(coset(ax, ay-1, az), coset(0, by-1, bz), l)+ &
                                       f3*s(coset(ax, ay, az), coset(0, by-2, bz), l)
                                    IF (ly1 > 0) s(coset(ax, ay, az), coset(0, by, bz), l) = &
                                       s(coset(ax, ay, az), coset(0, by, bz), l)+ &
                                       f2y*s(coset(ax, ay, az), coset(0, by-1, bz), ly1)
                                 END DO
                              END IF

!                   *** Shift of angular momentum component x from a to b ***

                              IF (ax == 0) THEN
                                 DO by = 0, lb-1
                                    bz = lb-1-by
                                    s(coset(ax, ay, az), coset(1, by, bz), l) = &
                                       rbp(1)*s(coset(ax, ay, az), coset(0, by, bz), l)
                                    IF (lx1 > 0) s(coset(ax, ay, az), coset(1, by, bz), l) = &
                                       s(coset(ax, ay, az), coset(1, by, bz), l)+ &
                                       f2x*s(coset(ax, ay, az), coset(0, by, bz), lx1)
                                 END DO
                                 DO bx = 2, lb
                                    f3 = f2*REAL(bx-1, dp)
                                    DO by = 0, lb-bx
                                       bz = lb-bx-by
                                       s(coset(ax, ay, az), coset(bx, by, bz), l) = &
                                          rbp(1)*s(coset(ax, ay, az), coset(bx-1, by, bz), l)+ &
                                          f3*s(coset(ax, ay, az), coset(bx-2, by, bz), l)
                                       IF (lx1 > 0) s(coset(ax, ay, az), coset(bx, by, bz), l) = &
                                          s(coset(ax, ay, az), coset(bx, by, bz), l)+ &
                                          f2x*s(coset(ax, ay, az), coset(bx-1, by, bz), lx1)
                                    END DO
                                 END DO
                              ELSE
                                 DO by = 0, lb-1
                                    bz = lb-1-by
                                    s(coset(ax, ay, az), coset(1, by, bz), l) = &
                                       rbp(1)*s(coset(ax, ay, az), coset(0, by, bz), l)+ &
                                       fx*s(coset(ax-1, ay, az), coset(0, by, bz), l)
                                    IF (lx1 > 0) s(coset(ax, ay, az), coset(1, by, bz), l) = &
                                       s(coset(ax, ay, az), coset(1, by, bz), l)+ &
                                       f2x*s(coset(ax, ay, az), coset(0, by, bz), lx1)
                                 END DO
                                 DO bx = 2, lb
                                    f3 = f2*REAL(bx-1, dp)
                                    DO by = 0, lb-bx
                                       bz = lb-bx-by
                                       s(coset(ax, ay, az), coset(bx, by, bz), l) = &
                                          rbp(1)*s(coset(ax, ay, az), coset(bx-1, by, bz), l)+ &
                                          fx*s(coset(ax-1, ay, az), coset(bx-1, by, bz), l)+ &
                                          f3*s(coset(ax, ay, az), coset(bx-2, by, bz), l)
                                       IF (lx1 > 0) s(coset(ax, ay, az), coset(bx, by, bz), l) = &
                                          s(coset(ax, ay, az), coset(bx, by, bz), l)+ &
                                          f2x*s(coset(ax, ay, az), coset(bx-1, by, bz), lx1)
                                    END DO
                                 END DO
                              END IF

                           END DO
                        END DO

                     END DO

                  END IF

               ELSE

                  IF (lb_max > 0) THEN

!             *** Vertical recurrence steps: [s|M|s] -> [s|M|b] ***

                     rbp(:) = (f1-1.0_dp)*rab(:)

!             *** [s|M|p] = (Pi - Bi)*[s|M|s] + f2*Ni(m)*[s|M-1i|s] ***

                     s(1, 2, l) = rbp(1)*s(1, 1, l)
                     s(1, 3, l) = rbp(2)*s(1, 1, l)
                     s(1, 4, l) = rbp(3)*s(1, 1, l)
                     IF (lx1 > 0) s(1, 2, l) = s(1, 2, l)+f2x*s(1, 1, lx1)
                     IF (ly1 > 0) s(1, 3, l) = s(1, 3, l)+f2y*s(1, 1, ly1)
                     IF (lz1 > 0) s(1, 4, l) = s(1, 4, l)+f2z*s(1, 1, lz1)

!             *** [s|M|b] = (Pi - Bi)*[s|M|b-1i] + f2*Ni(b-1i)*[s|M|b-2i] ***
!             ***           + f2*Ni(m)*[s|M-1i|b-1i]                      ***

                     DO lb = 2, lb_max

!               *** Increase the angular momentum component z of function b ***

                        s(1, coset(0, 0, lb), l) = rbp(3)*s(1, coset(0, 0, lb-1), l)+ &
                                                   f2*REAL(lb-1, dp)*s(1, coset(0, 0, lb-2), l)
                        IF (lz1 > 0) s(1, coset(0, 0, lb), l) = s(1, coset(0, 0, lb), l)+ &
                                                                f2z*s(1, coset(0, 0, lb-1), lz1)

!               *** Increase the angular momentum component y of function b ***

                        bz = lb-1
                        s(1, coset(0, 1, bz), l) = rbp(2)*s(1, coset(0, 0, bz), l)
                        IF (ly1 > 0) s(1, coset(0, 1, bz), l) = s(1, coset(0, 1, bz), l)+ &
                                                                f2y*s(1, coset(0, 0, bz), ly1)

                        DO by = 2, lb
                           bz = lb-by
                           s(1, coset(0, by, bz), l) = rbp(2)*s(1, coset(0, by-1, bz), l)+ &
                                                       f2*REAL(by-1, dp)*s(1, coset(0, by-2, bz), l)
                           IF (ly1 > 0) s(1, coset(0, by, bz), l) = s(1, coset(0, by, bz), l)+ &
                                                                    f2y*s(1, coset(0, by-1, bz), ly1)
                        END DO

!             *** Increase the angular momentum component x of function b ***

                        DO by = 0, lb-1
                           bz = lb-1-by
                           s(1, coset(1, by, bz), l) = rbp(1)*s(1, coset(0, by, bz), l)
                           IF (lx1 > 0) s(1, coset(1, by, bz), l) = s(1, coset(1, by, bz), l)+ &
                                                                    f2x*s(1, coset(0, by, bz), lx1)
                        END DO

                        DO bx = 2, lb
                           f3 = f2*REAL(bx-1, dp)
                           DO by = 0, lb-bx
                              bz = lb-bx-by
                              s(1, coset(bx, by, bz), l) = rbp(1)*s(1, coset(bx-1, by, bz), l)+ &
                                                           f3*s(1, coset(bx-2, by, bz), l)
                              IF (lx1 > 0) s(1, coset(bx, by, bz), l) = s(1, coset(bx, by, bz), l)+ &
                                                                        f2x*s(1, coset(bx-1, by, bz), lx1)
                           END DO
                        END DO

                     END DO

                  END IF

               END IF

            END DO

            DO k = 2, ncoset(lc_max)
               DO j = 1, ncoset(lb_max)
                  DO i = 1, ncoset(la_max)
                     mab(na+i, nb+j, k-1) = s(i, j, k)
                  END DO
               END DO
            END DO

            nb = nb+ncoset(lb_max)

         END DO

         na = na+ncoset(la_max)

      END DO

   END SUBROUTINE moment

! **************************************************************************************************
!> \brief This returns the derivative of the overlap integrals [a|b], with respect
!>       to the position of the primitive on the  left, i.e.
!>       [da/dR_ai|b] =  2*zeta*[a+1i|b] - Ni(a)[a-1i|b]
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param rab ...
!> \param difab ...
! **************************************************************************************************
   SUBROUTINE diffop(la_max, npgfa, zeta, rpgfa, la_min, &
                     lb_max, npgfb, zetb, rpgfb, lb_min, &
                     rab, difab)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rab
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(OUT)     :: difab

      INTEGER                                            :: lda_min, ldb_min, lds, lmax
      REAL(KIND=dp)                                      :: dab, rab2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: s
      REAL(KIND=dp), DIMENSION(npgfa*ncoset(la_max+1), &
         npgfb*ncoset(lb_max+1))                         :: sab

      rab2 = SUM(rab**2)
      dab = SQRT(rab2)

      lda_min = MAX(0, la_min-1)
      ldb_min = MAX(0, lb_min-1)
      lmax = MAX(la_max+1, lb_max+1)
      lds = ncoset(lmax+1)
      ALLOCATE (s(lds, lds, ncoset(1)))
      sab = 0.0_dp
      s = 0.0_dp
      CALL overlap(la_max+1, lda_min, npgfa, rpgfa, zeta, &
                   lb_max+1, ldb_min, npgfb, rpgfb, zetb, &
                   rab, dab, sab, 0, .FALSE., s, lds)

      CALL dabdr(la_max, npgfa, zeta, rpgfa, la_min, &
                 lb_max, npgfb, rpgfb, lb_min, &
                 dab, sab, difab(:, :, 1), difab(:, :, 2), difab(:, :, 3))

      DEALLOCATE (s)

   END SUBROUTINE diffop

! **************************************************************************************************
!> \brief This returns the derivative of the moment integrals [a|\mu|b], with respect
!>       to the position of the primitive on the  left, i.e.
!>       [da/dR_ai|\mu|b] =  2*zeta*[a+1i|\mu|b] - Ni(a)[a-1i|\mu|b]
!>       order indicates the max order of the moment operator to be calculated
!>       1: dipole
!>       2: quadrupole
!>       ...
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param order ...
!> \param rac ...
!> \param rbc ...
!> \param difmab ...
!> \param mab_ext ...
!> \note
! **************************************************************************************************
   SUBROUTINE diff_momop(la_max, npgfa, zeta, rpgfa, la_min, &
                         lb_max, npgfb, zetb, rpgfb, lb_min, &
                         order, rac, rbc, difmab, mab_ext)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min, order
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(OUT)  :: difmab
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: mab_ext

      INTEGER                                            :: imom, lda, lda_min, ldb, ldb_min, lmax
      REAL(KIND=dp)                                      :: dab, rab(3), rab2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: difmab_tmp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: mab

      rab = rbc-rac
      rab2 = SUM(rab**2)
      dab = SQRT(rab2)

      lda_min = MAX(0, la_min-1)
      ldb_min = MAX(0, lb_min-1)
      lmax = MAX(la_max+1, lb_max+1)
      lda = ncoset(la_max)*npgfa
      ldb = ncoset(lb_max)*npgfb
      ALLOCATE (difmab_tmp(lda, ldb, 3))

      IF (PRESENT(mab_ext)) THEN
         mab => mab_ext
      ELSE
         ALLOCATE (mab(npgfa*ncoset(la_max+1), npgfb*ncoset(lb_max+1), &
                       ncoset(order)-1))
         mab = 0.0_dp
!     *** Calculate the primitive overlap integrals ***
         CALL moment(la_max+1, npgfa, zeta, rpgfa, lda_min, &
                     lb_max+1, npgfb, zetb, rpgfb, &
                     order, rac, rbc, mab)

      END IF
      DO imom = 1, ncoset(order)-1
         difmab_tmp = 0.0_dp
         CALL adbdr(la_max, npgfa, rpgfa, la_min, &
                    lb_max, npgfb, zetb, rpgfb, lb_min, &
                    dab, mab(:, :, imom), difmab_tmp(:, :, 1), &
                    difmab_tmp(:, :, 2), difmab_tmp(:, :, 3))
         difmab(1:lda, 1:ldb, imom, 1) = difmab_tmp(1:lda, 1:ldb, 1)
         difmab(1:lda, 1:ldb, imom, 2) = difmab_tmp(1:lda, 1:ldb, 2)
         difmab(1:lda, 1:ldb, imom, 3) = difmab_tmp(1:lda, 1:ldb, 3)
      END DO

      IF (PRESENT(mab_ext)) THEN
         NULLIFY (mab)
      ELSE
         DEALLOCATE (mab)
      END IF
      DEALLOCATE (difmab_tmp)

   END SUBROUTINE diff_momop

! **************************************************************************************************
!> \brief This returns the derivative of the dipole integrals [a|x|b], with respect
!>       to the position of the primitive on the left and right, i.e.
!>       [da/dR_ai|\mu|b] =  2*zeta*[a+1i|\mu|b] - Ni(a)[a-1i|\mu|b]
!> \param la_max ...
!> \param npgfa ...
!> \param zeta ...
!> \param rpgfa ...
!> \param la_min ...
!> \param lb_max ...
!> \param npgfb ...
!> \param zetb ...
!> \param rpgfb ...
!> \param lb_min ...
!> \param order ...
!> \param rac ...
!> \param rbc ...
!> \param pab ...
!> \param forcea ...
!> \param forceb ...
!> \note
! **************************************************************************************************
   SUBROUTINE dipole_force(la_max, npgfa, zeta, rpgfa, la_min, &
                           lb_max, npgfb, zetb, rpgfb, lb_min, &
                           order, rac, rbc, pab, forcea, forceb)

      INTEGER, INTENT(IN)                                :: la_max, npgfa
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zeta, rpgfa
      INTEGER, INTENT(IN)                                :: la_min, lb_max, npgfb
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zetb, rpgfb
      INTEGER, INTENT(IN)                                :: lb_min, order
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rac, rbc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: forcea, forceb

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dipole_force', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, imom, ipgf, j, jpgf, lda, lda_min, &
                                                            ldb, ldb_min, lmax, na, nb
      REAL(KIND=dp)                                      :: dab, rab(3), rab2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: difmab, mab

      CPASSERT(order == 1)

      rab = rbc-rac
      rab2 = SUM(rab**2)
      dab = SQRT(rab2)

      lda_min = MAX(0, la_min-1)
      ldb_min = MAX(0, lb_min-1)
      lmax = MAX(la_max+1, lb_max+1)
      lda = ncoset(la_max)*npgfa
      ldb = ncoset(lb_max)*npgfb
      ALLOCATE (difmab(lda, ldb, 3))
      ALLOCATE (mab(npgfa*ncoset(la_max+1), npgfb*ncoset(lb_max+1), 3))
      mab = 0.0_dp
      CALL moment(la_max+1, npgfa, zeta, rpgfa, lda_min, &
                  lb_max+1, npgfb, zetb, rpgfb, 1, rac, rbc, mab)

      DO imom = 1, 3
         difmab = 0.0_dp
         CALL adbdr(la_max, npgfa, rpgfa, la_min, lb_max, npgfb, zetb, rpgfb, lb_min, &
                    dab, mab(:, :, imom), difmab(:, :, 1), difmab(:, :, 2), difmab(:, :, 3))
         na = 0
         DO ipgf = 1, npgfa
            nb = 0
            DO jpgf = 1, npgfb
               DO j = nb+ncoset(lb_min-1)+1, nb+ncoset(lb_max)
                  DO i = na+ncoset(la_min-1)+1, na+ncoset(la_max)
                     forceb(imom, 1) = forceb(imom, 1)+pab(i, j)*difmab(i, j, 1)
                     forceb(imom, 2) = forceb(imom, 2)+pab(i, j)*difmab(i, j, 2)
                     forceb(imom, 3) = forceb(imom, 3)+pab(i, j)*difmab(i, j, 3)
                  END DO
               END DO
               nb = nb+ncoset(lb_max)
            END DO
            na = na+ncoset(la_max)
         END DO

         difmab = 0.0_dp
         CALL dabdr(la_max, npgfa, zeta, rpgfa, la_min, lb_max, npgfb, rpgfb, lb_min, &
                    dab, mab(:, :, imom), difmab(:, :, 1), difmab(:, :, 2), difmab(:, :, 3))
         na = 0
         DO ipgf = 1, npgfa
            nb = 0
            DO jpgf = 1, npgfb
               DO j = nb+ncoset(lb_min-1)+1, nb+ncoset(lb_max)
                  DO i = na+ncoset(la_min-1)+1, na+ncoset(la_max)
                     forcea(imom, 1) = forcea(imom, 1)+pab(i, j)*difmab(i, j, 1)
                     forcea(imom, 2) = forcea(imom, 2)+pab(i, j)*difmab(i, j, 2)
                     forcea(imom, 3) = forcea(imom, 3)+pab(i, j)*difmab(i, j, 3)
                  END DO
               END DO
               nb = nb+ncoset(lb_max)
            END DO
            na = na+ncoset(la_max)
         END DO
      END DO

      DEALLOCATE (mab, difmab)

   END SUBROUTINE dipole_force

END MODULE ai_moments
