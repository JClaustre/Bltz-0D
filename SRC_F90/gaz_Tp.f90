!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 16/02/2016
! Objctv: Neutral Gaz temperature calculation.
! note  : boundary conditions : r=0 --> Tg(r) = cnste
!                               r=R --> d/dr(Tg) = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_TPGAZ
  USE F90_KIND  
  USE MOD_PARAM
  USE MOD_DGTSV
  IMPLICIT NONE

  REAL(DOUBLE), PARAMETER :: Tp0 = 600.d0

CONTAINS
  
  !***********************************************************************
  SUBROUTINE TP_Neutral (sys, elec, meta, Tg_p)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), INTENT(IN) :: elec
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(node), INTENT(INOUT) :: Tg_p
    !LOCAL
    INTEGER :: k
    INTEGER :: info = 0 !used for the DGTSV routine
    REAL(DOUBLE) :: Dx, Coef, Coef2, Aux
    REAL(DOUBLE), DIMENSION(node) :: ne, Kpa
    REAL(DOUBLE), DIMENSION(node) :: DL, D, DU, R

    Dx = sys%Ra / real(node-1)
    ne(:) = 0.d0
    DL = 0.d0 ; D = 0.d0; DU = 0.d0 ; R = 0.d0

    DO k = 1, node
       ne(k) = elec%Ni * bessj0(real(2.4048 * real(k-1)*Dx / sys%Ra))
       Kpa(k)= 0.152d0 * (Tg_p(k)/300.d0)**0.71
    END DO

    DO k = 1, node-1
       Coef = 0.666667d0*Clock%Dt / (meta(0)%Ni*kb*Dx*real(k)* Dx**2)
       Coef2 = 2.d0* ne(k) * me * Clock%Dt * meta(0)%Nuel(k) / (mhe*meta(0)%Ni)
       ! Lower boundary condition (Neumann Null)
       D(1) = -1.d0 ; Du(1) = 1.d0  
       R(1) = 0.d0!Dx * meta(0)%Tp*qok

       IF (k == 1) THEN
          Aux = OffDiag(Coef,1,Kpa,Dx)
       ELSE 
          Dl(k-1) = - Aux
          Aux     = OffDiag(Coef,k, Kpa, Dx)
          Du(k)   = - Aux
          D (k)   = 1.0d0 - ( Dl(k-1) + Du(k) ) + Coef2
          R (k)   = Tg_p(k) + Coef2 * elec%Tp*qok 
       END IF
    END DO
    ! Upper Boundary conditions (Dirichlet)
    D(node) = 1.d0 ; Dl(node-1) = 0.d0 
    R(node) = Tp0

    ! Calcul de la solution **************************
    CALL DGTSV (node, 1, DL, D, DU, R, node, info)
    ! ************************************************
    IF (info /= 0)THEN
       print*, "info error:", info
       STOP
    END IF
    ! *******************************************
    Tg_p(:) = R(:)
    meta(0)%Tp = ( sum(R(1:node/2)) / (node/2) ) * koq
    meta(0)%Ni = meta(0)%Prs / (qe * meta(0)%Tp * 7.5006d-3)

  CONTAINS
    FUNCTION OffDiag(Coef, k, Kpa, Dx)
      INTEGER     , INTENT(IN) :: k
      REAL(DOUBLE), INTENT(IN) :: Coef, Dx
      REAL(DOUBLE), DIMENSION(node), INTENT(IN) :: Kpa
      REAL(DOUBLE) :: OffDiag
      ! LOCAL
      REAL(DOUBLE) :: r, Kap
      kap = 0.5d0 * ( Kpa(k+1) + Kpa(k) )
      r = 0.5d0 * ( real(k)*Dx + real(k+1)*Dx )
      OffDiag = Kap * r * Coef
    END FUNCTION OffDiag
  END SUBROUTINE TP_Neutral

  !***********************************************************************
  FUNCTION bessj0(x)
    REAL bessj0,x
    !Returns the Bessel function J0(x) for any real x .
    REAL ax,xx,z
    DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,&
         r5,r6,s1,s2,s3,s4,s5,s6,y !We’ll accumulate polynomials in
                                   !double precision.
    SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
    DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,&
         -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,&
         .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
    DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d0,&
         -11214424.18d0,77392.33017d0,-184.9052456d0/,&
         s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,&
         9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
    if(abs(x).lt.8.)then
       !Direct rational function fit.
       y=x**2
       bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))&
            /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
    else
       !Fitting function (6.5.9).
       ax=abs(x)
       z=8./ax
       y=z**2
       xx=ax-.785398164
       bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))&
            -z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
    endif
    RETURN
  END FUNCTION bessj0
  !***********************************************************************

END MODULE MOD_TPGAZ
