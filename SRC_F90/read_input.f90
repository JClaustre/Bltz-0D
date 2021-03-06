!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 08/07/2015
! Objctv: Read input parameters, cross-sections and initialize
!         EEDF, densities, etc.
! note  : 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_READ

  USE F90_KIND
  USE MOD_PARAM
  USE MOD_EXCIT
  USE MOD_XCHANGE
  USE MOD_PENNASS
  USE MOD_RECOMB
  USE MOD_IONIZ
  IMPLICIT NONE

CONTAINS
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !***********************************************************************
  !                    SUBROUTINE READ_INPUT
  !***********************************************************************
  !**** Read the input file and init cross-sections
  SUBROUTINE Read_input (sys, ion, elec, Meta, lasr)
    IMPLICIT NONE
    !INTENT
    TYPE(SysVar), INTENT(INOUT) :: sys
    TYPE(Species), INTENT(INOUT) :: elec
    TYPE(Laser)  , INTENT(INOUT) :: lasr
    TYPE(Species), DIMENSION(NumIon), INTENT(INOUT) :: ion
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    !LOCAL 
    REAL(DOUBLE)       :: MassRatio  ! e-/He Gas Mass Ratio
    CHARACTER(len=2)   :: Name       ! Gas Name
    INTEGER, PARAMETER :: Npmax=1001 ! # of points max in file
    INTEGER :: i, j, k, l, idex
    INTEGER :: Npts, Npts0, expand, nul
    INTEGER :: NumA, NumQ, readA, readQ
    REAL(DOUBLE) :: EmaxF, Alpha, Eij, DU0, U0
    INTEGER     , DIMENSION(3) :: WhichQ, WhichA 
    INTEGER     , DIMENSION(43) :: Npts2 
    REAL(DOUBLE), DIMENSION(43) :: EmaxF2
    REAL(DOUBLE), DIMENSION(Npmax) :: SecRead, EnRead, SecML  ! Momentum c-s read
                                                              ! in file
    REAL(DOUBLE), DIMENSION(sys%nx+1)    :: SecMom
    REAL(DOUBLE), DIMENSION(Lv,Npmax) :: Sec !EXCIT/IONIZ C-S READ IN
    !FILE
    INTEGER :: LVG, LIG, NVG, LVVG !  ACTUAL # OF VIBRATIONAL STATES
                                   !  CONSIDERED IN BOLTZMANN
    INTEGER, DIMENSION(0:Lv) :: NVYES ! FOR EACH OF THE EXCITED STATES,
                                      ! NVYES IS 1 IF THIS STATE WILL BE
                                      ! CONSIDERED IN MCR
    REAL(DOUBLE), DIMENSION(0:2,0:Lv) :: Q, A  !Coeff
    INTEGER :: DEL01,XETHETA, IY2,IMOD,MOD_F,MOD, N1,N2,N3
    REAL(DOUBLE) :: PR,BS0,QS0, N4
    CHARACTER(len=10) :: A0
    LOGICAL :: file_exists
    SecMom=0.d0

    !**********************************************************************
    !                     READ DISCHARGE CONDITIONS 
    !**********************************************************************
    Write(*,"(2A)") tabul, "Starting Initialization "
    IF (Clock%Rstart == 0) THEN
       WRITE(*,"(3A)") tabul, 'Directory used : ', DirFile
       WRITE(*,"(2A)",advance="no") tabul, 'Reading Dschrge Condit : [input_he]'
       OPEN (UNIT=90,FILE='./datFile/input_he',STATUS='OLD')
    ELSE IF (Clock%Rstart == 1) THEN
       WRITE(*,"(3A)") tabul, 'Directory used : ', DirFile
       WRITE(*,"(2A)",advance="no") tabul, 'Reading Dschrge Condit : [Rs_input_he]'
       OPEN (UNIT=90,FILE=TRIM(ADJUSTL(DirFile))//'Rstart/Rs_input_he',STATUS='OLD')
    END IF
    READ (90,*) sys%nx
    READ (90,*) Clock%Rstart
    READ (90,*) sys%E
    READ (90,*) meta(0)%Ni, meta(0)%N0
    READ (90,*) sys%Freq, sys%rf
    READ (90,*) meta(0)%Prs
    READ (90,*) meta(0)%Tp, OneD%Tw
    READ (90,*) elec%Tp
    READ (90,*) sys%Ra
    READ (90,*) sys%L
    IF (Clock%Rstart == 0)THEN
       READ (90,*) elec%Ni
    ELSE
       READ (90,*) elec%Ni, ion(1)%Ni, ion(2)%Ni
    END IF
    READ (90,*) sys%Powr, sys%P0
    READ (90,*) Clock%SimuTime
    IF (Clock%Rstart == 0)THEN
       READ (90,*) Clock%MxDt
    ELSE
       READ (90,*) Clock%Dt
       Clock%MxDt = 5.d-09
    END IF
    READ (90,*) sys%Emx
    READ (90,*) Clock%TRstart
    READ (90,*) !/********************************/!
    READ (90,*) lasr%OnOff
    READ (90,*) lasr%Lwave
    READ (90,*) lasr%plz
    READ (90,*) lasr%Is
    READ (90,*) lasr%Sec
    READ (90,*) lasr%Stime
    READ (90,*) lasr%Ntr
    READ (90,*) (lasr%Ck(j), j=1,lasr%Ntr)
    IF (Clock%Rstart == 1) READ (90,*) Clock%SumDt
    CLOSE(90)
    !**********************************************************************
    IF (Clock%Rstart == 0) THEN
       Meta(0)%Tp = Meta(0)%Tp * koq! Tp °K to eV
       elec%Name = "      ELEC"
       Clock%SimuTime = Clock%SimuTime * 1.d-6
       Clock%TRstart  = Clock%TRstart  * 1.d-6
       Clock%Dt = 5.d-12
       !**** Viva USI unit !
       sys%E  = sys%E  * 1.d2
       sys%Freq = sys%Freq * 2.d0 * pi
       meta(0)%Ni = meta(0)%Ni * 1.d+6
       sys%Ra = sys%Ra * 1.d-2
       sys%L  = sys%L  * 1.d-2
       elec%Ni= elec%Ni * 1.d+6
       sys%Powr = sys%Powr * 1d+6
       lasr%Lwave = lasr%Lwave * 1d-9
       lasr%Sec = lasr%Sec * 1d-4
       lasr%Stime = lasr%Stime * 1d-6
    END IF
    !*******************************
    sys%Dx = sys%Emx / dble(sys%nx)

    OPEN(UNIT=50,FILE='./datFile/data_he.cs',ACTION="READ",STATUS="OLD")
    READ(50,*)
    READ(50,'(A2)')Name
    READ(50,*)MassRatio
    READ(50,*)
    !**************************************
    ! **** Momentum transfer cross section
    READ(50,*)
    READ(50,*)
    READ(50,*)npts, EmaxF

    IF(npts .GT. Npmax) THEN
       print*, 'WARNING [INREAD_INPUT()],npts > Npmax... not possible!'
       print*, '<<< Exiting program >>>'
       STOP
    END IF

    READ(50,*)
    l=0

11  READ(50,*) Npts0, expand
    READ(50,*)(SecRead(i), i=1,Npts0)
    IF(l .NE. 0 .AND. SecRead(1) .NE. SECML(L+1))THEN
       PRINT*,'DISCONNECTED DATA SETS IN CROSS-SECTION'
       PRINT*,'NPTS0:',Npts0, ' L+1=', l+1
       PRINT*,'SecRead(1):', SecRead(1),'SECML(L+1):', SECML(l+1)
       STOP
    ENDIF

    IF(expand.EQ.1)THEN
       DO i=1, Npts0-1
          j=l+i
          SECML(j)=SecRead(i)
       END DO
       SECML(j+1)=SecRead(Npts0)
       l=j
    ELSEIF(expand.GT.1)THEN
       DO i=1,Npts0-1
          PR = 0.d0
          DO k=1, expand
             j=expand*(i-1)+k+l
             SECML(j)=((DBLE(expand)-PR)*SecRead(i)+PR*SecRead(i+1))/DBLE(expand)
             PR=PR+1.d0
          END DO
       END DO
       SECML(j+1)=SecRead(NPTS0)
       l=j
    ELSEIF(expand.LT.-1)THEN
       DO i=1,Npts0-1,ABS(expand)
          j=l+1+(i-1)/ABS(expand)
          SECML(j)=SecRead(i)
       END DO
       SECML(j+1)=SecRead(Npts0)
       l=j
    END IF

    IF(j+1.EQ.Npts)THEN
       SECML(j+1)=SecRead(Npts0)
       l=j+1
    ENDIF

    IF(L.LT.Npts)GOTO 11
    !**************************************
    !**** Electronic excited levels
    READ(50,*)
    READ(50,'(I3)')LVG
    READ(50,'(I3)')LIG
    IF(LIG.GT.LVG) PRINT*,'ATTENTION IL Y A',LIG,'IONS,ET SEULEMENT',LVG,'NIVEAUX'
    READ(50,'(I3)')readQ
    IF(readQ.EQ.1) READ(50,*)NumQ,(WhichQ(i),i=1,NumQ)
    READ(50,'(I3)')readA
    IF(readA.EQ.1) READ(50,*)NumA,(WhichA(i),i=1,NumA)
    READ(50,*)
    READ(50,*)
    READ(50,*)
    NVYES(0)=1

    !****  READ SOME IDENTIFIERS FOR EACH STATE :
    !****  LEV:	STATE NAME
    !****  V  :	ENERGY
    !****  MOD_F:1 IF OSCILLATOR STRENGTHS ARE PROVIDED
    !****  NS :	S NUMBER
    !****  NL :	L NUMBER
    !****  MOD:	-FOR GROUND STATE ONLY- # OF CROSS-SECTIONS PROVIDED IN
    !****        TABLES
    !****  NVYES:FOR ALL STATES EXCEPT GROUND. NVYES= 1 IF EXCITATION FROM THIS
    !****        STATE IS CONSIDERED, AND IF THIS LEVEL IS CONSIDERED IN MCR MODEL

    DO idex=0, LVG
       IF(idex.EQ.0)THEN
          READ(50,'(A10,E12.4,5(1X,I5))') meta(idex)%Name,meta(idex)%En,&
               MOD_F,meta(idex)%Ns,meta(idex)%Nl,meta(idex)%Nn,MOD
       ELSEIF(idex .eq. 42+1) THEN
          READ(50,'(A10,E12.4,5(1X,I5))')ion(1)%Name,ion(1)%En,&
               MOD_F,ion(1)%Ns,ion(1)%Nl,ion(1)%Nn,NVYES(idex) 
       ELSE IF(idex .eq. 42+2) THEN
          READ(50,'(A10,E12.4,5(1X,I5))')ion(2)%Name,ion(2)%En,&
               MOD_F,ion(2)%Ns,ion(2)%Nl,ion(2)%Nn,NVYES(idex) 
       ELSE
          READ(50,'(A10,E12.4,5(1X,I5))')A0,N4,MOD_F,N1,N2,N3,NVYES(idex)
          IF (idex .LE. NumMeta) THEN
             meta(idex)%Name = A0 ; meta(idex)%En = N4
             meta(idex)%Ns = N1 ; meta(idex)%Nl = N2 ; meta(idex)%Nn = N3
          END IF
       END IF

       READ(50,*)
       IF(MOD_F.EQ.1)THEN
          READ(50,*)
          READ(50,*)(FOSC(idex,i),i=0,LVG)
          READ(50,*)
       END IF
       IF(readQ.EQ.1)THEN
          READ(50,*)
          READ(50,*)(Q(WhichQ(i),idex),i=1,NumQ)
          READ(50,*)
       END IF
       IF(readA.EQ.1)THEN
          READ(50,*)
          READ(50,*)(A(WhichA(I),idex),i=1,NumA)
          READ(50,*)
       END IF

       IF(idex.NE.0) GOTO 901

       DO IMOD=1,MOD
          READ(50,*)
          READ(50,*)IY2,Npts,EmaxF,nul
          READ(50,*)
          EmaxF2(IY2)=EmaxF
          Npts2(IY2)=Npts
          DO i=1,nul
             Sec(IY2,i)=0.0
          END DO
          l=nul
          IF(l.EQ.Npts)GOTO 900
13        READ(50,*)Npts0, expand
          READ(50,*)(SecRead(i),i=1,Npts0)
          IF(l.NE.nul .AND. SecRead(1).NE.SEC(IY2,l))THEN
             PRINT*,'DISCONNECTED DATA SETS EXCITATION CROSS-SECTION'
             PRINT*,'IY2=',IY2
             PRINT*,'NPTS0:',NPTS0,'L=',l
             PRINT*,'SecRead(1):',SecRead(1),'SEC(IY2,L):',SEC(IY2,l)
          END IF

          IF(expand.EQ.1)THEN
             DO i=1, Npts0-1
                j=l+i
                SEC(IY2,j)=SecRead(i)
             END DO
             SEC(IY2,j+1)=SecRead(Npts0)
             l=j
          ELSEIF(expand.GT.1)THEN
             DO i=1,Npts0-1
                PR=1.d0
                DO k=1,expand
                   j=expand*(i-1)+k+l
                   SEC(IY2,j)=((DBLE(expand)-PR)*SecRead(i)+PR*SecRead(i+1))/DBLE(expand)
                   PR=PR+1.d0
                END DO
             END DO
             l=j
          ELSEIF(expand.LT.-1)THEN
             DO i=1,NPTS0,ABS(expand)
                j=l+(i-1)/ABS(expand)
                SEC(IY2,j)=SecRead(i)
             END DO
             l=j
          END IF

          IF(j+1.EQ.Npts)THEN
             SEC(IY2,j+1)=SecRead(Npts0)
             l=j+1
          END IF
          IF(l.LT.Npts) GOTO 13
          READ(50,*)
900    END DO
901 END DO
    !**************************************

    NVG=0
    DO i=0,LV
       NVG=NVG+NVYES(i)
    ENDDO
    NVG=NVG-1

    READ(50,'(I3)')LVVG
    IF(LVVG.NE.0)THEN
       READ(50,*)DEL01
       READ(50,*)XETHETA
    ENDIF

    READ(50,*)
    READ(50,*)
    READ(50,*)
    READ(50,'(E10.4)')BS0
    READ(50,'(F4.2)')QS0
    READ(50,*)
    CLOSE(50)
    !**************************************
    write(*,"(A)") ' ......... Done'
    !**************************************
    IF (meta(0)%N0 == 1) THEN
       write(*,"(2A,ES10.2)") tabul, "Using Gas density entry Ng (cm-3) : ", meta(0)%Ni*1d-6
       meta(0)%Prs = meta(0)%Ni * qe * meta(0)%Tp *7.5006d-3
       write(*,"(2A,F8.2)") tabul, "Pressure is redefined (Torr) : ", meta(0)%Prs
    ELSEIF (meta(0)%N0 == 0) THEN
       write(*,"(2A,F8.2)") tabul, "Using Gas Pressure entry Prs (Torr) : ", meta(0)%Prs
       meta(0)%Ni = meta(0)%Prs / (qe * meta(0)%Tp * 7.5006d-3)
       write(*,"(2A,ES10.2)") tabul, "Gas density is redefined (cm-3) : ", meta(0)%Ni*1d-6
    END IF
    SELECT CASE (NumIon)                                                    !
    CASE (3)                                                                !
       ion(NumIon)%Name = " HE2-DIMER" ; ion(NumIon)%En = 17.94d0
       write(*,"(4A,F6.2)") tabul, "Init Excimer : ",ion(NumIon)%Name, " | Enrgy (eV) = ", ion(NumIon)%En
    END SELECT

    !**** Cross-Sec Elastic Momentum transfer
    OPEN(UNIT=51,FILE='./datFile/momentum.cs',ACTION="READ",STATUS="OLD")
    do i = 1, 8
       READ(51,*)
    END do
    READ(51,*) Npts 
    READ(51,*) (EnRead(i), i=1,Npts)
    READ(51,*) ; READ(51,*)
    READ(51,*) (SecRead(i), i=1,Npts)
    !**************************************
    !**** Interpolat Cross-Sect Momentum
    DO i=1, sys%nx
       Du0=0.d0
       U0 = IdU(i,sys%Dx)
       DO j = 1, Npts-1
          IF ( U0 == EnRead(j) ) meta(0)%SecMtM(i) = SecRead(j)
          IF ( U0 .gt. EnRead(j) .and. U0 .lt. EnRead(j+1)) Then
             Du0 = EnRead(j+1) - EnRead(j)
             meta(0)%SecMtM(i) = ((EnRead(j+1) - U0)*SecRead(j) )/Du0 &
                  + ((U0 - EnRead(j))*SecRead(j+1) )/Du0
          END IF
       END DO
    END DO
    meta(0)%SecMtM(:) = meta(0)%SecMtM(:) * 1e-20
    do i = 1, 3
       READ(51,*)
    END do
    READ(51,*) Npts 
    READ(51,*) (EnRead(i), i=1,Npts)
    READ(51,*) ; READ(51,*)
    READ(51,*) (SecRead(i), i=1,Npts)
    CLOSE(51)
    !**** Interpolat Cross-Sect Effective Momentum
    DO i=1, sys%nx
       Du0=0.d0
       U0 = IdU(i,sys%Dx)
       DO j = 1, Npts-1
          IF ( U0 == EnRead(j) ) meta(1)%SecMtM(i) = SecRead(j)
          IF ( U0 .gt. EnRead(j) .and. U0 .lt. EnRead(j+1)) Then
             Du0 = EnRead(j+1) - EnRead(j)
             meta(1)%SecMtM(i) = ((EnRead(j+1) - U0)*SecRead(j) )/Du0 &
                  + ((U0 - EnRead(j))*SecRead(j+1) )/Du0
          END IF
       END DO
    END DO
    !**************************************

    !**************************************
    !**** Cross-Sec Electron-Ion Momentum transfer
    DO i = 1, sys%nx
       U0 = IdU(i,sys%Dx)
       LnC = 23.d0 - log((elec%Ni*1d-6)**0.5d0 / elec%Tp**1.5d0)
       elec%SecEI(i) = qe**4  * LnC / (16.d0*pi*eps**2*U0)
    END DO
    !**************************************
    !**** Calcul Degenerescence
    DO i=0,NumMeta
       meta(i)%Deg=DBLE(FLOAT(meta(i)%Ns*(2*meta(i)%Nl+1)))
    END DO
    k=0
    DO i=24,LVG-3,8
       k=k+1
       Alpha=0.d0
       IF (i .LE. NumMeta) THEN
          DO j=1, k
             Alpha=Alpha+2*DBLE(FLOAT(meta(i)%Nl+j))+1.
          END DO
          meta(i)%Deg=meta(i)%Deg+Alpha*DBLE(FLOAT(meta(i)%Ns))
          meta(i+1)%Deg=meta(i+1)%Deg+Alpha*DBLE(FLOAT(meta(i+1)%Ns))
       END IF
    END DO
    !**************************************

    !**** Einstein transition probab : Aij
    DO i = 3, NumMeta
       DO j = 0, i-1
          Eij = meta(i)%En-meta(j)%En
          meta(i)%Aij(j) = 4.339d7 * Eij**2 * Fosc(j,i) * &
               meta(j)%Deg / meta(i)%Deg
          meta(i)%ondemit(j) = Hp * Vcel / (Eij*qe) ! (lambda = h * c / (Ef-Ei)  en metre)
          !IF (i==4.or.i == 3) print*, "Init: ", meta(i)%Name, meta(j)%Name, meta(i)%Aij(j), Fosc(j,i), meta(i)%ondemit(j)
       END DO
    END DO

    !**************************************
    !**** Associative process
    !**** He(n,l,s)+He --> He2+ + e
    CALL Init_Asso(Sn, 0)
    !**************************************
    !**** l-change atomic process
    !**** He(n,l,s)+He <--> He(n,l',s)+He
    CALL Init_lchange (meta, K_ij)
    !**************************************
    !**** Inelastic electron cross section 
    !**** He(n,l,s)+e <--> He(n',l',s')+e
    CALL Init_ExcDxc(sys, meta, Fosc, Q, A)
    !**************************************
    !**** Recombination electron cross sec
    !**** He2+ + e --> He(n,l,s) + He(1S)
    CALL Init_Recomb (sys, meta)
    !**************************************
    !**** Ionization cross sec
    !**** He(n,l,s) + e --> He+ + 2e
    CALL Init_ioniz (sys, meta)
    !**************************************

    !**** READ Probability transitions for laser's pumping
    OPEN(UNIT=51,FILE='./datFile/pumping.cs',ACTION="READ",STATUS="OLD")
    do i = 1, 9
       READ(51,*)
    END do
    do i = 1, 18 ! Read transition elements
       READ(51,*) (Tij(i,j,1), j=1,6)
    END do
    do i = 1, 5
       READ(51,*)
    END do
    do i = 1, 18 ! Read polarization elements
       READ(51,*) (Tij(i,j,2), j=1,6)
    END do
    do i = 1, 3
       READ(51,*)
    END do
    DO i = 1, 18 ! Read Ck transition components
       READ(51,*) (Tij(i,j,3), j=1,6)
    END DO
    do i = 1, 5
       READ(51,*)
    END do
    DO i = 1, 6 ! Read E Matrix
       READ(51,*) (lasr%Eij(i,j), j=1,6)
    END DO
    READ(51,*)
    DO i = 1, 6 ! Read F Matrix
       READ(51,*) (lasr%Fij(i,j), j=1,6)
    END DO
    READ(51,*)
    READ(51,*) (lasr%lamb(j), j=1,6)
!    do i = 1, 9 ! Read Tij coordinates
!       READ(51,*) Ck(i)%Name, Ck(i)%Ti(1), Ck(i)%Ti(2),Ck(i)%Tj(1), Ck(i)%Tj(2), Ck(i)%Ediff
!    END do
    CLOSE(51)
    !**************************************

    !**** Check if directories are OK.
    write(*,"(2A)",advance="no") tabul, "Check directories ... "
    inquire( file=DirFile, exist=file_exists)
    if ( file_exists ) then
       write(*,"(2A)") tabul, "directory exists!"
    else
       !**** workaround: it calls an extern program...
       call execute_command_line ('mkdir -p ' // adjustl(trim( dirFile ) ) )
       call execute_command_line ('mkdir -p ' // adjustl(trim( dirFile ) )//"Rstart/" )
       write(*,"(2A)") "dir doesn't exist! I create it for youuuuuu! --> ", DirFile
    end if
    !**************************************
    CALL date_and_time(Clock%DDay(1),Clock%DDay(2),Clock%DDay(3), Clock%Date)

    write(*,"(2A)") tabul, "Initilization Done"
    write(*,"(A,3(A,I4))") tabul, "Today : ",Clock%Date(3),"/",Clock%Date(2),"/",Clock%Date(1)
    write(*,"(A,3(A,I2))") tabul, "Hour : ", Clock%Date(5),"h ",Clock%Date(6),"min ",Clock%Date(7)
    write(*,"(2A)") tabul, "*************************************************"
  END SUBROUTINE Read_input
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Init (sys, clock, ion, elec, meta, lasr)
    !INTENT
    TYPE(SysVar) , INTENT(INOUT) :: sys
    TYPE(Time)   , INTENT(INOUT) :: clock
    TYPE(Species), INTENT(INOUT) :: elec
    TYPE(Laser)  , INTENT(INOUT) :: lasr
    TYPE(Species), DIMENSION(:) , INTENT(INOUT) :: ion
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    !LOCAL
    INTEGER :: i, j
    REAL(DOUBLE) :: power, Uc, Df, Tp_tmp
    !**** Look for Restart or not
    OPEN  (UNIT=90,FILE='./datFile/input_he',STATUS='OLD')
    READ (90,*) sys%Nx
    READ (90,*) Clock%Rstart
    CLOSE (90)
    !**** Max number of iteration for the simulation : beyond that number... it's insane! 
    clock%MaxIter = int(5E+08)

    IF (Clock%Rstart == 0) clock%SumDt = 0.d0
    !**** Init Grid Variables
    IF (Clock%Rstart == 1) THEN
       OPEN  (UNIT=90,FILE=TRIM(ADJUSTL(DirFile))//'Rstart/Rs_input_he',STATUS='OLD')
       READ  (90,*) sys%Nx
       CLOSE (90)
    END IF
    !**** Alloc Arrays
    CALL AllocArray(sys%nx)
    !**** Read Init
    CALL Read_Input(sys, ion, elec, meta, lasr)
    !**** Init EEDF
    IF (Clock%Rstart == 1) OPEN (UNIT=90,FILE=TRIM(ADJUSTL(DirFile))//'Rstart/EEDF.dat',STATUS='OLD')
    DO i = 1, sys%nx
       U(i)  = IdU(i,sys%Dx)
       meta(0)%Nuel(i) = meta(0)%Ni*meta(0)%SecMtm(i)*gama*dsqrt(U(i))
       OneD%nuMoy = OneD%nuMoy + meta(0)%SecMtm(i)*gama*dsqrt(U(i))

       !**** Maxwllian distribution function\
       IF (Clock%Rstart == 1) THEN
          READ(90,"(I6, ES19.10)") j, F(i)
       ELSE IF (i.LT.sys%nx)THEN
          F(i) = ( 2.d0*elec%Ni / sqrt(pi*elec%Tp**3) ) * exp( -(U(i)/elec%Tp))
       END IF

       consv(1) = consv(1) + F(i)*U(i)**(0.5d0)*sys%Dx
       consv(2) = consv(2) + F(i)*U(i)**(1.5d0)*sys%Dx
       !elec%Nuei(i) = gama*dsqrt(U(i))*elec%SecEI(i)*elec%Ni
    END DO
    
    IF (Clock%Rstart == 1) THEN
       CLOSE (90)
       OPEN (UNIT=99,FILE=TRIM(ADJUSTL(DirFile))//'Rstart/Tg.dat',STATUS='OLD')
       DO  i = 1, size(OneD%Tg)
          READ(99,*) j, OneD%Tg(i)
       END DO
       CLOSE (99)
    END IF
    sys%Volume = pi * sys%Ra**2 * sys%L
    elec%Ni = consv(1)
    write(*,"(2A, F6.2,A)"  ) tabul, "Tpe init : ", 0.66667d0*consv(2)/consv(1), " (eV)"
    write(*,"(2A, ES10.2,A)") tabul,"Ne  init : ", consv(1)*1d-6, " (cm-3) "
    
    !**** Fix the power here function of Elec field ************************!
    IF (sys%P0 .EQ. 1) THEN                                                 !
       DO i = 1, sys%nx - 1                                                 !
          IF (sys%rf == 1) THEN
             Uc = qome * sys%E**2 / (meta(0)%Nuel(i)**2 )                   !
          ELSE
             Uc = qome * sys%E**2 / (meta(0)%Nuel(i)**2 + sys%Freq**2)      !
          END IF
          Df = (F(i+1)-F(i))                                                !
          power = power - Uc * U(i)**(1.5d0) * Df * meta(0)%Nuel(i)*0.6667d0!
       END DO                                                               !
       power = power * qe                                                   !
       sys%Powr = power                                                     !
       write (*,"(2A,ES15.4,A,ES15.6)") tabul, "Power [W/cm3 | W] fixed to : ", &
            sys%Powr*1d-6, " |", sys%Powr*sys%volume                        !
    END IF                                                                  !
    !***********************************************************************!

    !**** 1D profil for gas temperature calculation
    OneD%SLab = 0.0d0 ! (m)
    OneD%nx = size(OneD%Tg)
    OneD%Dx = (sys%Ra + OneD%SLab) / real(OneD%nx-1)
    IF (OneD%SLab .NE. 0.d0) THEN
       OneD%bnd = int(sys%Ra / OneD%Dx)
    ELSE
       OneD%Bnd = OneD%nx
    END IF
    OneD%nuMoy = OneD%nuMoy / sys%nx

    IF (Clock%Rstart == 0)  THEN
       ! Gas temperature in the cylinder (Parabolic profile)
       Tp_tmp = (3.d0/2.d0) * (meta(0)%Tp*qok - OneD%Tw/3.d0)
       DO i = 1, OneD%bnd
          OneD%Tg(i) = Tp_tmp - ((i-1)*OneD%Dx / sys%Ra)**2 * (Tp_tmp-OneD%Tw)
       END DO
       IF (OneD%SLab .NE. 0.d0) THEN
          Do i = OneD%bnd+1, OneD%nx
             OneD%Tg(i) = (i-OneD%bnd+1)*OneD%Dx*(300.d0-OneD%Tw)/OneD%SLab&
                  + OneD%Tw ! Room temperature (K)
          END Do
       END IF
    END IF

    IF (meta(0)%N0 == 1) THEN
       OneD%ng(:) =  meta(0)%Ni
       OneD%Pg(:) =  meta(0)%Ni * (qe * OneD%Tg(:) * koq * 7.5006d-3)
    ELSE
       OneD%Pg(:) =  meta(0)%Prs
       OneD%ng(:) =  meta(0)%Prs / (qe * OneD%Tg(:) * koq * 7.5006d-3)
    END IF

    !**** Init Densities (Ions + excited states) (m-3) *********************!
    IF (Clock%Rstart == 0) THEN                                             !
       ion(2)%Ni = elec%Ni * 0.00d0                                         !
       ion(1)%Ni = elec%Ni * 1.00d0                                         !
       SELECT CASE (NumIon)                                                 !
       CASE (3) ; ion(NumIon)%Ni = 1.0d+0     ! Molecular Excimer           !
       END SELECT                                                           !
       DO i = 1, NumMeta                                                    !
          IF (i.EQ.1) meta(i)%Ni = 1.0d+0     ! Metastable 2S3              !
          IF (i.EQ.2) meta(i)%Ni = 1.0d+0     ! Metastable 2S1              !
          IF (i.EQ.3) meta(i)%Ni = 1.0d+0     ! Radiative state 2P3         !
          IF (i.GE.4) meta(i)%Ni = 1.00+0     ! others excited states       !  
       END DO                                                               !
       !**** Allocate densities for sublevels in 2S3 and 2P3 ***            !
       pop(1)%Ni(:) = meta(1)%Ni/6.d0 ; pop(2)%Ni(:) = meta(3)%Ni/18.d0     !
    ELSE                                                                    !
       OPEN (UNIT=90,FILE=TRIM(ADJUSTL(DirFile))//'Rstart/Density.dat',STATUS='OLD')
       READ(90,*) (meta(i)%Ni, i=1,NumMeta), &                              !
            (pop(1)%Ni(i),i=1,6), (pop(2)%Ni(i), i=1,18)                    !
       SELECT CASE (NumIon)                                                 !
       CASE (3) ; READ(90,*) ion(NumIon)%Ni, &                              !
            (pop(1)%Ni(i),i=1,6), (pop(2)%Ni(i), i=1,18)                    !
       END SELECT                                                           !
       CLOSE (90)                                                           !
       !**** Check electron conservation from Rstart files ***              !
       IF (Clock%Rstart.EQ.1) THEN                                          !
          IF (ion(1)%Ni.GT.ion(2)%Ni) THEN                                  !
             ion(2)%Ni = elec%Ni - ion(1)%Ni                                !
          ELSE                                                              !
             ion(1)%Ni = elec%Ni - ion(2)%Ni                                !
          END IF                                                            !
       END IF                                                               !
    END IF                                                                  !
    !***********************************************************************!
    CALL Write_Out1D( F,  "F_init.dat")

  END SUBROUTINE Init
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Rstart_SaveFiles (sys, clock, ion, elec, meta, pop, F)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Time)   , INTENT(IN) :: clock
    TYPE(Species), INTENT(IN) :: elec
    TYPE(Excited), DIMENSION(2) , INTENT(IN) :: pop
    TYPE(Species), DIMENSION(:) , INTENT(IN) :: ion
    TYPE(Species), DIMENSION(0:), INTENT(IN) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN) :: F
    !LOCAL
    INTEGER :: i

    !**** Save Energy Electron Distribution Function
    OPEN(UNIT=990,File=TRIM(ADJUSTL(DirFile))//"Rstart/EEDF.dat",ACTION="WRITE",STATUS="UNKNOWN")
    DO i = 1, sys%nx
       write(990,"(I6, ES19.10)") i, F(i)
    END DO
    CLOSE(990)
    !**** Save excited states density
    OPEN(UNIT=990,File=TRIM(ADJUSTL(DirFile))//"Rstart/Density.dat",ACTION="WRITE",STATUS="UNKNOWN")
    write(990,"(42ES19.10,6ES19.10,18ES19.10)") (meta(i)%Ni, i=1,NumMeta), &
         (pop(1)%Ni(i), i=1,6), (pop(2)%Ni(i), i=1,18)
    SELECT CASE (NumIon) 
    CASE (3) ; write(990,"(ES19.10,6ES19.10,18ES19.10)") ion(NumIon)%Ni, &
         (pop(1)%Ni(i), i=1,6), (pop(2)%Ni(i), i=1,18)
    END SELECT
    CLOSE(990)
    !**** Save neutral gas temperature profile
    OPEN(UNIT=990,File=TRIM(ADJUSTL(DirFile))//"Rstart/Tg.dat",ACTION="WRITE",STATUS="UNKNOWN")
    DO i = 1, OneD%nx
       write(990,"(I6, ES15.6)") i, OneD%Tg(i)
    END DO
    CLOSE(990)
    !**** Save Parameters
    OPEN(UNIT=990,File=TRIM(ADJUSTL(DirFile))//"Rstart/Rs_input_he",ACTION="WRITE",STATUS="UNKNOWN")
    !*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    WRITE (990,"(I6)") sys%nx
    !*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    WRITE (990,"(I4)")     1
    WRITE (990,"(ES15.6)") sys%E
    WRITE (990,"(ES15.6,I2)") meta(0)%Ni, meta(0)%N0
    WRITE (990,"(ES15.6,I2)") sys%Freq, sys%rf
    WRITE (990,"(ES15.6)") meta(0)%Prs
    WRITE (990,"(2ES15.6)") meta(0)%Tp, OneD%Tw
    WRITE (990,"(ES15.6)") elec%Tp
    WRITE (990,"(ES15.6)") sys%Ra
    WRITE (990,"(ES15.6)") sys%L
    WRITE (990,"(3ES19.10)") elec%Ni, (ion(i)%Ni, i=1,2)
    WRITE (990,"(ES15.6,I2)") sys%Powr, sys%P0
    WRITE (990,"(ES15.6)") Clock%SimuTime
    WRITE (990,"(ES15.6)") Clock%Dt
    WRITE (990,"(ES15.6)") sys%Emx
    WRITE (990,"(ES15.6)") Clock%TRstart
    WRITE (990,"(A)") "!********************************!"
    WRITE (990,"(I4)") lasr%OnOff
    WRITE (990,"(ES15.6)") lasr%Lwave
    WRITE (990,"(I4)") lasr%plz
    WRITE (990,"(ES15.6)") lasr%Is
    WRITE (990,"(ES15.6)") lasr%Sec
    WRITE (990,"(ES15.6)") lasr%Stime
    WRITE (990,"(I4)") lasr%Ntr
    WRITE (990,"(9I4)") (lasr%Ck(i), i=1,lasr%Ntr)
    !*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    WRITE (990,"(ES15.6)") Clock%SumDt
    !*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    CLOSE (990)
    !***********************************************
  END SUBROUTINE Rstart_SaveFiles
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
END MODULE MOD_READ
