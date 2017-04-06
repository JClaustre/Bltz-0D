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
     REAL(DOUBLE) :: Dt=0, SimuTime
     REAL(DOUBLE) :: SumDt, TRstart
     REAL(DOUBLE) :: Hours, Minutes, Seconds
     INTEGER      :: MaxIter, IterDt, NumIter, Rstart
  END type Time
  !-----------------------------------------------------------
  TYPE, PUBLIC::SysVar
     INTEGER :: nx, P0, rf! node number and max-node number
     REAL(DOUBLE) :: Emx, Dx ! grid step
     REAL(DOUBLE) :: Ra, L, volume
     REAL(DOUBLE) :: E, Emoy, Freq, Powr, IPowr, Pcent, Pwmoy
  END type SysVar
  !-----------------------------------------------------------
  TYPE, PUBLIC::Species
     INTEGER           :: Nn, Ns, Nl, N0
     REAL(DOUBLE)      :: Ni, Tp, Prs, En, Deg, Damb
     REAL(DOUBLE)      :: Dfree, mobl, Updens, J, NStart
     CHARACTER(len=10) :: Name
     REAL(DOUBLE), DIMENSION(:), POINTER :: Aij, Nuel, Nuei
     REAL(DOUBLE), DIMENSION(:), POINTER :: SecRec, SecTot, SecMtM, SecEI
     REAL(DOUBLE), DIMENSION(:,:), POINTER :: SecIon, SecExc
  END type Species
  !-----------------------------------------------------------
  TYPE, PUBLIC::Excited
     CHARACTER(len=4) :: Name
     REAL(DOUBLE)     :: Dn_o, Ntot
     REAL(DOUBLE)     :: T_relax, Te, Tr, tau_e
     REAL(DOUBLE)     :: polarz
     REAL(DOUBLE), DIMENSION(:), POINTER :: Ni
  END type Excited
  !-----------------------------------------------------------
  TYPE, PUBLIC::laser
     INTEGER      :: plz, OnOff, Ntr
     REAL(DOUBLE) :: Is, sec, Lwave, Stime  
     INTEGER, DIMENSION(9) :: Ck
     INTEGER, DIMENSION(6,6) :: Eij, Fij
     INTEGER, DIMENSION(6) :: lamb
  END type Laser
  !-----------------------------------------------------------
  TYPE, PUBLIC::Diagnos
     REAL(DOUBLE)      :: SumTx, EnProd, EnLoss
     REAL(DOUBLE)      :: InM1, OutM1, InM2, OutM2
     REAL(DOUBLE), DIMENSION(3) :: Tx
     CHARACTER(len=10) :: Name
  END type Diagnos
  !-----------------------------------------------------------
  TYPE, PUBLIC::profil1D
     INTEGER :: nx, bnd
     REAL(DOUBLE) :: SLab, Dx, nuMoy, Tw
     REAL(DOUBLE), DIMENSION(500) :: Tg, Pg, ng, ne, nu
  END type profil1D
  !-----------------------------------------------------------

  ! *******************************************************************************
  INTEGER, PARAMETER :: Lv=44
  INTEGER, PARAMETER :: NumIon  = 3  ! He+ | He2+ | He2*
  INTEGER, PARAMETER :: NumMeta = 34 ! 1S1 --> 7P1
  INTEGER, PARAMETER :: Npop1 = 6    ! Sublevel numbers in 2S3 
  INTEGER, PARAMETER :: Npop2 = 18   ! Sublevel numbers in 2P3

  !CHARACTER(*), PARAMETER :: DirFile = "./datFile/post-Dischrg/760_Torr/0.01_microS/"
  CHARACTER(*), PARAMETER :: DirFile = "./datFile/Default_RunDir_2/"

  TYPE(Time)    :: Clock
  TYPE(SysVar)  :: sys
  TYPE(Species) :: elec
  TYPE(profil1D):: OneD
  TYPE(Laser)   :: lasr
  TYPE(Diagnos), DIMENSION(18) :: diag
  TYPE(Species), DIMENSION(NumIon)    :: ion
  TYPE(Species), DIMENSION(0:NumMeta) :: meta ! (0) --> fundamental state
  TYPE(Excited), DIMENSION(2) :: pop

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
  REAL(DOUBLE), PARAMETER :: PPi = 2.d0*Pi 
  REAL(DOUBLE), PARAMETER :: gama= dsqrt(2.d0*qome)
  REAL(DOUBLE), PARAMETER :: MassR = me/mhe ! Mass Ratio (me/mhe)
  REAL(DOUBLE), PARAMETER :: Ry  = 13.605692  ! Rydberg energy (eV)
  REAL(DOUBLE), PARAMETER :: Nlosh = 2.6868d+25 ! Loschmidt Number (m-3) => (P=1atm, T=0C)
  REAL(DOUBLE), PARAMETER :: Vcel = 2.9979d+08 ! speed of light in vacuum (m/s)
  REAL(DOUBLE), PARAMETER :: Hp   = 6.6261d-34  ! Planck Constant (J s)
  REAL(DOUBLE), PARAMETER :: Hpb  = Hp / (2.d0*Pi) ! H_barre (J s) 
  REAL(DOUBLE), PARAMETER :: fineS= qe*qe / (2*eps*Hp*Vcel) ! Fine structure constant
  !REAL(DOUBLE), PARAMETER :: LnC = 10.d0      ! lnC = ln(Λ) log Coulomb (cf. Fk-Pl)
  CHARACTER(len=1), PARAMETER :: tabul=char(9) ! tabulation 
  ! *******************************************************************************
  REAL(DOUBLE), DIMENSION(:)  , ALLOCATABLE :: F
  REAL(DOUBLE), DIMENSION(:)  , ALLOCATABLE :: U
  REAL(DOUBLE), DIMENSION(2)                :: consv
  REAL(DOUBLE), DIMENSION(34)               :: Sn   ! Associative rate coeff
  REAL(DOUBLE), DIMENSION(NumMeta,NumMeta)  :: K_ij ! l-change rate coeff
  REAL(DOUBLE), DIMENSION(0:Lv,0:Lv)        :: Fosc ! Oscillator strenght
  REAL(DOUBLE), DIMENSION(Npop2,Npop1,3)    :: Tij  ! Transition for each Ck components
  REAL(DOUBLE) :: MaxR                      ! Max rate calculated --> used for adaptative time
  REAL(DOUBLE) :: LnC                       ! lnC = ln(Λ) log Coulomb (cf. Fk-Pl)
  REAL(DOUBLE) :: Vg                        ! Calculate sheath potential for diffusion routine
  REAL(DOUBLE) :: Twnsd_a                   ! First Townsend coefficient
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
    ALLOCATE ( Meta(0)%SecMtM(nx) ) ; Meta(0)%SecMtM(:) = 0.d0 ! Elastic momentum transfer
    ALLOCATE ( Meta(1)%SecMtM(nx) ) ; Meta(1)%SecMtM(:) = 0.d0 ! Effective momentum transfer
    ALLOCATE ( Meta(0)%SecRec(nx) ) ; Meta(0)%SecRec(:) = 0.d0
    ALLOCATE ( Meta(0)%Nuel(nx)  ) ; Meta(0)%Nuel(:)   = 0.d0

    ALLOCATE ( F(nx) ) ; F(:) = 0.d0
    ALLOCATE ( U(nx) ) ; U(:) = 0.d0
    ALLOCATE ( elec%SecEI(nx) ) ; elec%SecEI(:) = 0.d0
    ALLOCATE ( elec%Nuei(nx) )  ; elec%Nuei(:)  = 0.d0

    ALLOCATE ( pop(1)%Ni(6)  ) ; pop(1)%Ni(:) = 0.d0
    ALLOCATE ( pop(2)%Ni(18) ) ; pop(2)%Ni(:) = 0.d0

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

    DEALLOCATE ( Meta(0)%SecTot, Meta(0)%SecMtM, Meta(1)%SecMtM, Meta(0)%SecRec )
    DEALLOCATE ( F, U, Meta(0)%Nuel )
    DEALLOCATE ( pop(1)%Ni, pop(2)%Ni )
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

    OPEN(UNIT=0999,File=TRIM(ADJUSTL(DirFile))//TRIM(ADJUSTL(FileName)),&
         ACTION="WRITE",STATUS="UNKNOWN")
    DO i = lbound(array,1), ubound(array,1)
       write(0999,*) i, array(i)
    END DO
    CLOSE(0999)
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

    OPEN(UNIT=99,File=TRIM(ADJUSTL(DirFile))//TRIM(ADJUSTL(FileName)),&
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

    OPEN(UNIT=99,File=TRIM(ADJUSTL(DirFile))//TRIM(ADJUSTL(FileName)),&
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
       write(*,"(F4.1,A,I6,A)") tot, " sec ! (Num Loop = ", Nloop, ")  "
    ELSEIF (tot .GT. 60 .and. tot .LE. 5*60.0) THEN
       write(*,"(2A)", advance="no") tabul, &
               "You Have time to ... check your email! estimated calculation: "
       write(*,"(F5.2,A,I6,A)") tot/60.d0, " min ! (Num Loop = ", Nloop, ")  "
    ELSEIF (tot .GT. 5*60.0 .and. tot .LE. 15*60.0) THEN
       write(*,"(2A)", advance="no") tabul, &
            "You have time for coffee! estimated calculation: "
       write(*,"(F5.2,A,I6,A)") tot/60.d0, " min ! (Num Loop = ", Nloop, ")  "
    ELSEIF (tot .GT. 15*60.0 .and. tot .LE. 45*60.0) THEN
       write(*,"(2A)", advance="no") tabul, &
            "estimated calculation : "
       write(*,"(F5.2,A,I6,A)") tot/60.d0, " min ! (Num Loop = ", Nloop, ")  "
    ELSEIF (tot .GT. 45*60.0 .and. tot .LE. 3600.) THEN
       write(*,"(2A)", advance="no") tabul, &
            "You have time for ... Optimizations! estimated calculation: "
       write(*,"(F5.2,A,I6,A)") tot/60.d0, " min ! (Num Loop = ", Nloop, ")  "
    ELSEIF (tot .GT. 3600. ) THEN
       write(*,"(2A)", advance="no") tabul, &
            "You have all the time! estimated calculation: "
       write(*,"(2(I3,A),I8,A)") int(tot/3600.d0), "H", int(((tot/3600.)-int(tot/3600))*60.), "min | (Num Loop = ", Nloop, ")  "
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
END MODULE MOD_PARAM
