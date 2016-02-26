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
  SUBROUTINE TP_Neutral (sys, elec, meta, OneD)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), INTENT(IN) :: elec
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    TYPE(profil1D), INTENT(INOUT) :: OneD
    !LOCAL
    INTEGER :: k, Nmoy, SizeE, Med
    INTEGER :: info = 0 !used for the DGTSV routine
    REAL(DOUBLE) :: Dx, Coef, Coef2, nuMoy
    REAL(DOUBLE) :: A, B, C, D
    REAL(DOUBLE) :: tm, tp, t2m, t2p
    REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: DL, DI, DU, R

    OneD%nx = SIZE(OneD%Tg)
    ALLOCATE ( DL(OneD%nx), DI(OneD%nx), DU(OneD%nx), R(OneD%nx) )
    Nmoy = int(OneD%nx/2)
    Dx = sys%Ra / real(OneD%nx-1)
    DL = 0.d0 ; DI = 0.d0; DU = 0.d0 ; R = 0.d0

    DO k = 1, OneD%nx
       OneD%ne(k) = elec%Ni * bessj0(real(2.4048 * real(k)*Dx / sys%Ra))
    END DO
    SizeE = size(meta(0)%Nuel(:))
    nuMoy = sum(meta(0)%Nuel(1:SizeE-1)) / (SizeE-1)

    DO k = 1, OneD%nx-1
       Coef = 0.666667d0*Clock%Dt / (meta(0)%Ni*kb*Dx**2)
       Coef2 = 2.d0* MassR * Clock%Dt * nuMoy *OneD%ne(k-1) / meta(0)%Ni
       ! Lower boundary condition (Neumann Null)
       Di(1) = 1.d0 ; Du(1) = -1.d0  
       R(1) = 0.d0
       Med = 1
       IF (k == 1) THEN
          A = Off1(Coef,1,OneD%Tg,Dx, Med)
          B = Off2(1,OneD%Tg,0.71d0)
       ELSE
          
          IF (k-1 == 1) THEN
             Dl(k-1) = - A*( 1.d0 - B )/(Dx*real(k))
             C = A/(Dx*real(k))
             D = B
          ELSE
             Dl(k-1) = - A*( 1.d0 - B )
             C = A ; D = B
          END IF
          
          A = Off1(Coef,k,OneD%Tg,Dx, Med)/(Dx*real(k))
          B = Off2(k,OneD%Tg,0.71d0)
          Du(k) = - A*( 1.d0 + B )

          Di (k) = 1.0d0 + A*(1.d0-B) + C*(1.d0 + D) 
          R (k)  = A*(OneD%Tg(k+1)-OneD%Tg(k)) - C*(OneD%Tg(k)-OneD%Tg(k-1))

          IF (k .LE. OneD%bnd) THEN
             Di (k) = Di (k) + Coef2 
             R (k)  = R (k)  + Coef2 * (elec%Tp*qok - OneD%Tg(k))
          END IF

       END IF
    END DO
    ! Upper Boundary conditions (Dirichlet) 
    ! Carefull /!\ here we solve W = T^k+1 - T^k
    Di(OneD%nx) = 1.d0 ; Dl(OneD%nx) = 0.d0 
    R(OneD%nx) = 0.d0

    ! Calcul de la solution **************************
    CALL DGTSV (OneD%nx, 1, DL, DI, DU, R, OneD%nx, info)
    ! ************************************************
    IF (info /= 0)THEN
       print*, "info error:", info
       STOP
    END IF
    ! *******************************************
    OneD%Tg(:OneD%nx-1) = R(:OneD%nx-1)+OneD%Tg(:OneD%nx-1)
    OneD%Tg(OneD%nx) = Tp0
    OneD%Pg(:) =  meta(0)%Ni * (qe * OneD%Tg(:) * koq * 7.5006d-3)

    meta(0)%Tp  = ( sum(OneD%Tg(1:Nmoy)) / (Nmoy) ) * koq
    meta(0)%Prs = ( sum(OneD%Pg(1:Nmoy)) / (Nmoy) )
    DEALLOCATE (DL, DI, DU, R)

  CONTAINS
    FUNCTION Off1(Coef, k, Tg, Dx, M)
      INTEGER     , INTENT(IN) :: k, M
      REAL(DOUBLE), INTENT(IN) :: Coef, Dx
      REAL(DOUBLE), DIMENSION(:), INTENT(IN) :: Tg
      REAL(DOUBLE) :: Off1
      ! LOCAL
      REAL(DOUBLE) :: r, Kap
      IF (M == 1) kap = 0.156d0   * ((Tg(k+1)+Tg(k))/300.d0)**0.710 / 2.d0
      IF (M == 2) kap = 0.02623d0 * ((Tg(k+1)+Tg(k))/300.d0)**0.788 / 2.d0
      r = 0.5d0 * ( real(k)*Dx + real(k+1)*Dx )
      Off1 = Coef * Kap * r
    END FUNCTION Off1

    FUNCTION Off2(k, Tg, nu)
      INTEGER     , INTENT(IN) :: k
      REAL(DOUBLE), DIMENSION(:), INTENT(IN) :: Tg
      REAL(DOUBLE) :: nu, Off2
      REAL(DOUBLE) :: Tp, Tp2
      Tp = Tg(k+1) - Tg(k) 
      Tp2 = ( Tg(k+1) + Tg(k) )
      Off2 = nu * Tp / Tp2
    END FUNCTION Off2

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
