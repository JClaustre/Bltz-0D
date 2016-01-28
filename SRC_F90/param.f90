!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 15/05/2015
! Objctv: List all physical constants and 
!         parameters needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MOD_PARAM

  USE F90_KIND
  IMPLICIT NONE

  !-----------------------------------------------------------
  TYPE, PUBLIC :: Time
     REAL(DOUBLE) :: Dt=0, SimuTime, Year, DoY, DoY_LT
     REAL(DOUBLE) :: Hours, Minutes, Seconds
     REAL(DOUBLE) :: SumDt, TRstart
     INTEGER :: MaxIter, IterDt, NumIter, Rstart
  END type Time
  !-----------------------------------------------------------
  TYPE, PUBLIC::SysVar
     INTEGER :: nx, P0 ! node number and max-node number
     REAL(DOUBLE) :: Emx, Dx ! grid step
     REAL(DOUBLE) :: Ra, L, volume
     REAL(DOUBLE) :: E, Eef, Freq, Omg, Powr, IPowr
  END type SysVar
  !-----------------------------------------------------------
  TYPE, PUBLIC::Species
     INTEGER           :: Nn, Ns, Nl, N0
     REAL(DOUBLE)      :: Ni, Tp, Prs, En, Deg, Damb
     REAL(DOUBLE)      :: Dfree, mobl, Updens
     CHARACTER(len=10) :: Name
     REAL(DOUBLE), DIMENSION(:), POINTER :: Aij, Nuel
     REAL(DOUBLE), DIMENSION(:), POINTER :: SecRec, SecTot, SecMtM
     REAL(DOUBLE), DIMENSION(:,:), POINTER :: SecIon, SecExc
  END type Species
  !-----------------------------------------------------------
  TYPE, PUBLIC::Diagnos
     REAL(DOUBLE)      :: Tx, EnProd, EnLoss
     CHARACTER(len=10) :: Name
  END type Diagnos
   !-----------------------------------------------------------

  ! *******************************************************************************
  INTEGER, PARAMETER :: Lv=44
  INTEGER, PARAMETER :: NumIon  = 3  ! He+ | He2+ | He2*
  INTEGER, PARAMETER :: NumMeta = 34 ! 1S1 --> 7P1

  TYPE(Time)    :: Clock
  TYPE(SysVar)  :: sys
  TYPE(Species) :: elec
  TYPE(Diagnos), DIMENSION(15) :: diag
  TYPE(Species), DIMENSION(NumIon)    :: ion
  TYPE(Species), DIMENSION(0:NumMeta) :: meta ! (0) --> fundamental state

  REAL(DOUBLE), PARAMETER :: kb  = 1.3807d-23 ! Boltzmann constant (m2 kg s-2 K-1)
  REAL(DOUBLE), PARAMETER :: qe  = 1.602d-19  ! Elementary charge (C)
  REAL(DOUBLE), PARAMETER :: koq = kb/qe      ! conversion to eV
  REAL(DOUBLE), PARAMETER :: qok = qe/kb      ! conversion to K
  REAL(DOUBLE), PARAMETER :: me  = 9.109d-31  ! Electron mass (kg)
  REAL(DOUBLE), PARAMETER :: qome= qe/me
  REAL(DOUBLE), PARAMETER :: mi  = 1.672d-27  ! Proton mass (kg)
  REAL(DOUBLE), PARAMETER :: mhe = 4.002*mi   ! Helium mass (kg)
  REAL(DOUBLE), PARAMETER :: eps = 8.8542d-12 ! Permittivity of free space (F.m-1)
  REAL(DOUBLE), PARAMETER :: Pi  = 4*ATAN(1.d0)
  REAL(DOUBLE), PARAMETER :: gama= dsqrt(2.d0*qome)
  REAL(DOUBLE), PARAMETER :: MassR = 1.3710d-04 ! Mass Ratio (me/mi)
  REAL(DOUBLE), PARAMETER :: Ry  = 13.605692  ! Rydberg energy (eV)
  REAL(DOUBLE), PARAMETER :: Nlosh = 2.6868d+25 ! Loschmidt Number (m-3) => (P=1atm, T=0C)
  !REAL(DOUBLE), PARAMETER :: LnC = 10.d0      ! lnC = ln(Λ) log Coulomb (cf. Fk-Pl)
  CHARACTER(len=1), PARAMETER :: tabul=char(9) ! tabulation 
  ! *******************************************************************************
  REAL(DOUBLE), DIMENSION(:)  , ALLOCATABLE :: F
  REAL(DOUBLE), DIMENSION(:)  , ALLOCATABLE :: U
  REAL(DOUBLE), DIMENSION(2)                :: consv
  REAL(DOUBLE), DIMENSION(NumMeta)          :: Sn   ! Associative rate coeff
  REAL(DOUBLE), DIMENSION(NumMeta,NumMeta)  :: K_ij ! l-change rate coeff
  REAL(DOUBLE), DIMENSION(0:Lv,0:Lv)        :: Fosc ! Oscillator strenght
  REAL(DOUBLE) :: MaxR

CONTAINS

  ! **********************************************************
  FUNCTION IdU(i,Dx)
    REAL(DOUBLE) :: IdU
    INTEGER     , INTENT(IN) :: i
    REAL(DOUBLE), INTENT(IN) :: Dx
    IdU = (dble(i) - 0.5d0) * Dx
  END FUNCTION IdU
  ! **********************************************************

  SUBROUTINE AllocArray(nx)
    INTEGER :: i
    INTEGER, INTENT(IN) :: nx

    DO i = 0, NumMeta
       ALLOCATE ( Meta(i)%SecIon(2 ,nx) ) ; Meta(i)%SecIon(:,:) = 0.d0
       ALLOCATE ( Meta(i)%SecExc(0:NumMeta,nx) ) ; Meta(i)%SecExc(:,:) = 0.d0
       ALLOCATE ( Meta(i)%Aij(0:NumMeta) ) ; Meta(i)%Aij(:) = 0.d0
    END DO

    SELECT CASE (NumIon)
    CASE (3)
       ALLOCATE ( ion(NumIon)%SecIon(2 ,nx) ) ; ion(NumIon)%SecIon(:,:) = 0.d0
       ALLOCATE ( ion(NumIon)%SecExc(1 ,nx) ) ; ion(NumIon)%SecExc(:,:) = 0.d0
    END SELECT
    
    ALLOCATE ( Meta(0)%SecTot(nx) ) ; Meta(0)%SecTot(:) = 0.d0
    ALLOCATE ( Meta(0)%SecMtM(nx) ) ; Meta(0)%SecMtM(:) = 0.d0
    ALLOCATE ( Meta(0)%SecRec(nx) ) ; Meta(0)%SecRec(:) = 0.d0
    ALLOCATE ( Meta(0)%Nuel(nx)  ) ; Meta(0)%Nuel(:)   = 0.d0

    ALLOCATE ( F(nx) ) ; F(:) = 0.d0
    ALLOCATE ( U(nx) ) ; U(:) = 0.d0

  END SUBROUTINE AllocArray
  ! **********************************************************
  SUBROUTINE DelocArray()
    INTEGER :: i
    
    DO i = 0, NumMeta
       DEALLOCATE ( Meta(i)%SecIon, Meta(i)%SecExc, Meta(i)%Aij )
    END DO
    SELECT CASE (NumIon)
    CASE (3) 
       DEALLOCATE ( ion(NumIon)%SecIon, ion(NumIon)%SecExc )
    END SELECT

    DEALLOCATE ( Meta(0)%SecTot, Meta(0)%SecMtM, Meta(0)%SecRec )
    DEALLOCATE ( F, U, Meta(0)%Nuel )
  END SUBROUTINE DelocArray
  ! **********************************************************

  !****************** USEFULL FUNCTION ***********************
  !**** Following Subroutines :
  !**** - Write_Out1D
  !**** - Write_Out2D
  !**** - Write_Out2D
  !**** - PrinTime
  !**** - LoopTime
  !**** - tridag

  ! **********************************************************


  !***********************************************************
  !                    SUBROUTINE Write_Out1D
  !***********************************************************
  SUBROUTINE Write_Out1D( array, FileName )
    !INTENT
    REAL(DOUBLE), DIMENSION(:), INTENT(IN) :: array
    CHARACTER(*), INTENT(IN) :: FileName
    !LOCAL
    INTEGER :: i

    OPEN(UNIT=99,File="./datFile/"//TRIM(ADJUSTL(FileName)),&
         ACTION="WRITE",STATUS="UNKNOWN")
    DO i = lbound(array,1), ubound(array,1)
       write(99,*) i, array(i)
    END DO
    CLOSE(99)
  END SUBROUTINE Write_Out1D
  !***********************************************************
  !                    SUBROUTINE Write_Out2D
  !***********************************************************
  SUBROUTINE Write_Out2D( array, FileName )
    !INTENT
    REAL(DOUBLE), DIMENSION(:,:), INTENT(IN) :: array
    CHARACTER(*), INTENT(IN) :: FileName
    !LOCAL
    INTEGER :: i, j

    OPEN(UNIT=99,File="./datFile/"//TRIM(ADJUSTL(FileName)),&
         ACTION="WRITE",STATUS="UNKNOWN")
    DO i = lbound(array,1), ubound(array,1)
       DO j = lbound(array,2), ubound(array,2)
          write(99,*) i, j, array(i,j)
       END DO
       write(99,*) ""
    END DO
    CLOSE(99)
  END SUBROUTINE Write_Out2D

  !***********************************************************
  !                    SUBROUTINE Write_Out3D
  !***********************************************************
  SUBROUTINE Write_Out3D( array, FileName )
    !INTENT
    REAL(DOUBLE), DIMENSION(:,:,:), INTENT(IN) :: array
    CHARACTER(*), INTENT(IN) :: FileName
    !LOCAL
    INTEGER :: i, j, k

    OPEN(UNIT=99,File="./datFile/"//TRIM(ADJUSTL(FileName)),&
         ACTION="WRITE",STATUS="UNKNOWN")
    DO i = lbound(array,1), ubound(array,1)
       DO j = lbound(array,2), ubound(array,2)
          DO k = lbound(array,3), ubound(array,3)
             write(99,*) i, j, k, array(i,j,k)
          END DO
          write(99,*) ""
       END DO
       write(99,*) ""
    END DO
    CLOSE(99)
  END SUBROUTINE Write_Out3D
  !***********************************************************
  SUBROUTINE PrinTime (t1, t2, rate)
    !INTENT
    INTEGER, INTENT(IN) :: t1, t2, rate
    !LOCAL
    INTEGER :: min, Hrs
    REAL(DOUBLE) :: sec
    Hrs = 0; min = 0; sec = 0.d0

    sec = real(t2-t1) / real(rate)
    if (sec .GE. 60.d0) then
       min = int(sec / 60.d0)
       sec = ((sec / 60.d0) - min ) * 60.d0
       IF (min .GE. 60) THEN
          Hrs = int(min / 60)
          min = int((real(min/60.0d0) - Hrs)*60)
          sec = (real(min/60.d0) - min)*60
       END IF
    END if
    Clock%Hours = Hrs ; Clock%Minutes = min ; Clock%Seconds = sec
    write(*,"(2A,I3,2(A,I2),A)") tabul, "Elapsed CPU Time : ", &
         Hrs, "h", min, ":", int(sec), "s"
  END SUBROUTINE PrinTime
  !***********************************************************
  SUBROUTINE LoopTime(t1, t2, clockR, Nloop)
    INTEGER, INTENT(IN)    :: t1, t2, clockR
    INTEGER, INTENT(INOUT) :: Nloop
    REAL(DOUBLE) :: sec, tot
    sec = 0.d0 ; tot = 0.d0
    sec = real(t2-t1) / real(ClockR)
    Nloop = int( (Clock%SimuTime-Clock%SumDt) /Clock%Dt) + 300

    tot = sec * Nloop /100.d0

    IF (tot .LE. 60.d0) THEN
       write(*,"(2A)", advance="no") tabul, &
            "Hang on to your seat ... calculation time is estimated to "
       write(*,"(F4.1,A,I6,A)") tot, " sec ! (Num Loop = ", Nloop, ")"
    ELSEIF (tot .GT. 60 .and. tot .LE. 5*60.0) THEN
       write(*,"(2A)", advance="no") tabul, &
               "You Have time to check your email ... calculation time is estimated to "
       write(*,"(F5.2,A,I6,A)") tot/60.d0, " min ! (Num Loop = ", Nloop, ")"
    ELSEIF (tot .GT. 5*60.0 .and. tot .LE. 15*60.0) THEN
       write(*,"(2A)", advance="no") tabul, &
            "You can take your coffee ... calculation time is estimated to "
       write(*,"(F5.2,A,I6,A)") tot/60.d0, " min ! (Num Loop = ", Nloop, ")"
    ELSEIF (tot .GT. 15*60.0 .and. tot .LE. 45*60.0) THEN
       write(*,"(2A)", advance="no") tabul, &
            "Think about something else! ... calculation time is estimated to "
       write(*,"(F5.2,A,I6,A)") tot/60.d0, " min ! (Num Loop = ", Nloop, ")"
    ELSEIF (tot .GT. 45*60.0 .and. tot .LE. 3600.) THEN
       write(*,"(2A)", advance="no") tabul, &
            "Maybe you should think about Optimizations ... calculation time is estimated to "
       write(*,"(F5.2,A,I6,A)") tot/60.d0, " min ! (Num Loop = ", Nloop, ")"
    ELSEIF (tot .GT. 3600. ) THEN
       write(*,"(2A)", advance="no") tabul, &
            "One suggestion, look at Parallelization! ... calculation time is estimated to "
       write(*,"(2(I3,A),I8,A)") int(tot/3600.d0), "H", int(((tot/3600.)-int(tot/3600))*60.), "min | (Num Loop = ", Nloop, ")"
    END IF
  END SUBROUTINE LoopTime
  !***********************************************************    
  SUBROUTINE tridag(a,b,c,r,u,n) 
    INTEGER j,n 
    REAL(DOUBLE) :: gam(50000),a(n),b(n),c(n),u(n),r(n), bet
    if (b(1) .EQ. 0d0) Then
       print*, 'ERROR : b(1)=0 in tridag' 
       Stop
    End if
    bet=b(1) 
    u(1)=r(1)/bet 
    do j=2,n 
       gam(j)=c(j-1)/bet 
       bet=b(j)-a(j)*gam(j) 
       if (bet .EQ. 0d0) Then
          print*, 'ERROR : bet=0 in tridag' 
          Stop
       End if
       u(j)=(r(j)-a(j)*u(j-1))/bet 
    end do
    do j=n-1,1,-1 
       u(j)=u(j)-gam(j+1)*u(j+1) 
    end do
  END SUBROUTINE tridag

  !***********************************************************
  FUNCTION arth(first,increment,n)
    REAL(DOUBLE), INTENT(IN) :: first,increment
    INTEGER, INTENT(IN) :: n
    REAL(DOUBLE), DIMENSION(n) :: arth
    INTEGER :: k,k2
    REAL(DOUBLE) :: temp
    if (n > 0) arth(1)=first
    if (n <= 16) then
       do k=2,n
          arth(k)=arth(k-1)+increment
       end do
    else
       do k=2,8
          arth(k)=arth(k-1)+increment
       end do
       temp=increment*8
       k=8
       do
          if (k >= n) exit
          k2=k+k
          arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth

  !***********************************************************
  FUNCTION gammln(xx)
    REAL, INTENT(IN) :: xx
    REAL :: gammln
    !Returns the value ln[Γ( xx )] for xx > 0.
    REAL(DOUBLE) :: tmp,x
    !Internal arithmetic will be done in double precision, a nicety that
    !you can omit if five-figure accuracy is good enough.
    REAL(DOUBLE) :: stp = 2.5066282746310005d0
    REAL(DOUBLE), DIMENSION(6) :: coef 

    coef = (/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,&
    -1.231739572450155d0,0.1208650973866179d-02,-0.5395239384953d-05/)
    !call assert(xx > 0.0, ’gammln_s arg’)
    x=xx
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    gammln=tmp+log(stp*(1.000000000190015d0+&
         sum(coef(:)/arth(x+1.0d0,1.0d0,size(coef))))/x)
  END FUNCTION gammln

  !***********************************************************    
  FUNCTION gser(a,x,gln)
    REAL, INTENT(IN) :: a,x
    REAL, OPTIONAL, INTENT(OUT) :: gln
    REAL :: gser
    INTEGER, PARAMETER :: ITMAX=100
    REAL, PARAMETER :: EPS=1.0d-10
    !Returns the incomplete gamma function P (a, x) evaluated by its
    !series representation as gamser . Also optionally returns ln Γ(a)
    !as gln .
    INTEGER :: n
    REAL :: ap,del,summ
    if (x == 0.0) then
       gser=0.0
       RETURN
    end if
    ap=a
    summ=1.0/a
    del=summ
    do n=1,ITMAX
       ap=ap+1.0
       del=del*x/ap
       summ=summ+del
       if (abs(del) < abs(summ)*EPS) exit
    end do
    if (n > ITMAX) print*, "a too large, ITMAX too small in gser_s"
    if (present(gln)) then
       gln=gammln(a)
       gser=summ*exp(-x+a*log(x)-gln)
    else
       gser=summ*exp(-x+a*log(x)-gammln(a))
    end if
  END FUNCTION gser

  FUNCTION gcf(a,x,gln)
    REAL, INTENT(IN) :: a,x
    REAL, OPTIONAL, INTENT(OUT) :: gln
    REAL :: gcf
    INTEGER, PARAMETER :: ITMAX=100
    REAL(DOUBLE), PARAMETER :: EPS=1.0d-10,FPMIN=1.0d-30
    !Returns the incomplete gamma function Q(a, x) evaluated by its
    !continued fraction repre- sentation as gammcf . Also optionally
    !returns ln Γ(a) as gln .  Parameters: ITMAX is the maximum
    !allowed number of iterations; EPS is the relative accu- racy;
    !FPMIN is a number near the smallest representable floating-point
    !number.
    INTEGER :: i
    REAL :: an,b,c,d,del,h
    if (x == 0.0) then
       gcf=1.0
       RETURN
    end if
    b=x+1.0-a
    !Set up for evaluating continued fraction by modified Lentz’s
    !method (§5.2) with b 0 = 0.
    c=1.0/FPMIN
    d=1.0/b
    h=d
    do i=1,ITMAX
       !Iterate to convergence.
       an=-i*(i-a)
       b=b+2.0
       d=an*d+b
       if (abs(d) < FPMIN) d=FPMIN
       c=b+an/c
       if (abs(c) < FPMIN) c=FPMIN
       d=1.0/d
       del=d*c
       h=h*del
       if (abs(del-1.0) <= EPS) exit
    end do
    if (i > ITMAX) print*, "a too large, ITMAX too small in gcf_s"
    if (present(gln)) then
       gln=gammln(a)
       gcf=exp(-x+a*log(x)-gln)*h
       !Put factors in front.
    else
       gcf=exp(-x+a*log(x)-gammln(a))*h
    end if
  END FUNCTION gcf

  !***********************************************************
  FUNCTION gammp(a,x)
    REAL, INTENT(IN) :: a,x
    REAL :: gammp
    !Returns the incomplete gamma function P (a, x).
    !call assert( x >= 0.0, a > 0.0, ’gammp_s args’)
    if (x<a+1.0) then
       !Use the series representation.
       gammp=gser(a,x)
    else
       !Use the continued fraction representation and take its
       !complement.
       gammp=1.0-gcf(a,x)
    end if
  END FUNCTION gammp

  !***********************************************************
  FUNCTION gammq(a,x)
    REAL, INTENT(IN) :: a,x
    REAL :: gammq
    !Returns the incomplete gamma function Q(a, x) ≡ 1 − P (a, x).
    !call assert( x >= 0.0, a > 0.0, ’gammq_s args’)
    if (x<a+1.0) then
       !Use the series representation and take its complement.
       gammq=1.0-gser(a,x)
    else
       !Use the continued fraction representation.
       gammq=gcf(a,x)
    end if
  END FUNCTION gammq

  !***********************************************************
END MODULE MOD_PARAM
