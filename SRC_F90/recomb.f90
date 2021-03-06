!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 08/07/2015
! Objctv: Dissociative Recombination processes in He
! note  : Norm : use the rates defined in "Alves et al", and 
!                renormalized the EEDF
!         Recomb: use the cross-section given in "Santos et al"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_RECOMB
  USE F90_KIND  
  USE MOD_PARAM
  IMPLICIT NONE

CONTAINS

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Recomb (sys, meta, Neut, U, Fi, diag, iter)
    !INTENT
    INTEGER      , INTENT(IN) :: iter
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(2) , INTENT(INOUT) :: Neut
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
    TYPE(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag

    !LOCAL
    INTEGER :: i
    REAL(DOUBLE) :: coef, recmb, Dx, ratx!, Rtmp, part
    REAL(DOUBLE) :: energI, energF, U3
    REAL(DOUBLE), DIMENSION(4) :: tx
    tx = (/0.011d0, 0.341d0, 0.645d0, 0.003d0/) ! Santos et al.
    !tx = (/0.037d0, 0.360d0, 0.586d0, 0.017d0/) ! Pedersen et al
    Dx = sys%Dx ; recmb=0.d0 ; Coef = 0.d0 !; Rtmp = 0.d0 ; part = 0.d0
    energI = 0.d0 ; energF = 0.d0

    IF (mod(iter,2000)==0) OPEN(UNIT=919,File=TRIM(ADJUSTL(DirFile))//"Rates_all.dat",&
         ACTION="WRITE",POSITION="APPEND",STATUS='OLD')

    DO i = 1, sys%nx
       U3 = U(i)*U(i)*U(i)
       !**** Initial Energy density
       energI = energI + Fi(i) * sqrt(U3) * sys%Dx
       !****************************************************
       coef = U(i) * meta(0)%SecRec(i) * Fi(i) * gama * ion(2)%Ni
       !Rtmp = Rtmp + ( U(i) * meta(0)%SecRec(i) * Fi(i) * gama * Dx )
       recmb = recmb + coef
       Fi(i) = Fi(i) - Clock%Dt * coef / dsqrt( U(i) )
       !**** Final Energy density
       energF = energF + Fi(i) * sqrt(U3) * sys%Dx
       !****************************************************
    END DO
    ion(2)%Updens  = ion(2)%Updens  - Clock%Dt * recmb * Dx
    !**** Ref. Branching ratio in Santos et al (j.phys D:47 (2014)) 
    Do i = 1, 4
       meta(i)%Updens = meta(i)%Updens + Clock%Dt * (tx(i)*recmb) * Dx
       !**** UpDate neutral density for polarization
       Neut(1)%Updens = Neut(1)%Updens + Clock%Dt * (tx(i)*recmb) * Dx
    END Do
    !**** Energy conservation Diagnostic
    diag(8)%EnLoss = diag(8)%EnLoss + (energI-energF)
    !**** Rate calcul for adaptative time-step
    ratx = recmb * Dx / ion(2)%Ni
    IF (ratx .GT. maxR) maxR = ratx
    !****************
    diag(8)%SumTx =  diag(8)%SumTx + Clock%Dt * recmb * Dx
    !***************** Diagnostic for relative importance of reactions (m-3/s)
    diag(8)%Tx(1) =  recmb * Dx
    !*************** Diagnostic for metastable and 2^3P rates (m-3 s-1)
    diag(8)%InM1 = (tx(1)*recmb) * Dx
    diag(8)%InM2 = (tx(3)*recmb) * Dx
    !****************
    IF (mod(iter,2000)==0) THEN
       WRITE(919,"(A)") "Rates & density (s-1 | m-3 s-1 | m-3) in Disso. Recomb. : e + He2+ -> He* + He(1S)"
       WRITE(919,"(A,3ES12.3)") ion(2)%Name, recmb*Dx/ion(2)%Ni, recmb*Dx, ion(2)%Ni
    END IF
    IF (mod(iter,2000)==0) CLOSE(919)

  END SUBROUTINE Recomb
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Recomb_Norm (sys, meta, U, Fi, diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
    TYPE(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag

    !LOCAL
    INTEGER :: i
    REAL(DOUBLE) :: coef, recmb, Dx
    REAL(DOUBLE) :: energI, energF, U3
    Dx = sys%Dx ; recmb=0.d0 ; Coef = 0.d0
    energI = 0.d0 ; energF = 0.d0

    DO i = 1, sys%nx
       U3 = U(i)*U(i)*U(i)
       !**** Initial Energy density
       energI = energI + Fi(i) * sqrt(U3) * sys%Dx
       !****************************************************
       coef = U(i) * 1.06d-20 * Fi(i) * gama * ion(2)%Ni
       recmb = recmb + coef
       Fi(i) = Fi(i) - Clock%Dt * coef / dsqrt( U(i) )
       !**** Final Energy density
       energF = energF + Fi(i) * sqrt(U3) * sys%Dx
       !****************************************************
    END DO
    ion(2)%Updens  = ion(2)%Updens  - Clock%Dt * recmb * Dx
    meta(1)%Updens = meta(1)%Updens + Clock%Dt * recmb * Dx

    !**** Energy conservation Diagnostic
    diag(8)%EnLoss = diag(8)%EnLoss + (energI-energF)
    diag(8)%SumTx = diag(8)%SumTx + Clock%Dt * recmb * Dx
    !***************** Diagnostic for relative importance of reactions
    diag(8)%Tx(1) =  recmb * Dx
    !****************
  END SUBROUTINE Recomb_Norm

  !***********************************************************************
  SUBROUTINE Recomb_Alves (sys, meta, U, Fi, diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
    TYPE(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    !LOCAL
    INTEGER :: i
    REAL(DOUBLE) :: coef, recmb, part, En1, En2
    En1  = 0.d0 ; En2 = 0.d0 ; part = 0.d0
    recmb = 5.0d-09 * 1.d-6 * (meta(0)%Tp / (elec%Tp)) ! m3 s-1
    coef = recmb * Clock%Dt * elec%Ni * ion(2)%Ni

    do i = 1, sys%nx
       En1 = Fi(i) * U(i)**1.5d0 * sys%Dx
       part = part + Fi(i) * sqrt(U(i)) * sys%Dx
    end do

    elec%Ni = part - coef
    ion(2)%UpDens  = ion(2)%UpDens  - coef
    meta(1)%UpDens = meta(1)%UpDens + coef
    do i = 1, sys%nx
       Fi(i) = Fi(i) * elec%Ni / part
       En2 = Fi(i) * U(i)**1.5d0 * sys%Dx
    END do

    !**** Energy conservation Diagnostic
    diag(8)%EnLoss = diag(8)%EnLoss + abs(En1 - En2)
    diag(8)%SumTx = diag(8)%SumTx + coef

  END SUBROUTINE Recomb_Alves
  !***********************************************************************

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Conv_3Body (meta, ion, Neut, iter)
    !INTENT
    INTEGER      , INTENT(IN) :: iter
    TYPE(Species), DIMENSION(0:), INTENT(InOut) :: meta
    TYPE(Species), DIMENSION(:) , INTENT(InOUT) :: ion
    TYPE(Species), DIMENSION(2) , INTENT(InOUT) :: Neut
    !LOCAL
    REAL(DOUBLE) :: eta, Src, Tp, excim, ratx
    Tp = meta(0)%Tp * qok
    IF (mod(iter,2000)==0) OPEN(UNIT=919,File=TRIM(ADJUSTL(DirFile))//"Rates_all.dat",&
         ACTION="WRITE",POSITION="APPEND",STATUS='OLD')

    !*************************************
    !**** cf. Alves 1992 : η (cm6 s-1)
    !**** He+ + 2He(1S) --> He2+ + He(1S)
    eta = 1.4d-31 * 1.d-12 * (Tp/300.d0)**(-0.6d0) 
    Src = eta * ion(1)%Ni * meta(0)%Ni**2

    IF (ion(1)%Ni .GT. 0.d0) THEN
       ion(1)%Updens = ion(1)%Updens - Clock%Dt * Src
       ion(2)%Updens = ion(2)%Updens + Clock%Dt * Src
       !**** UpDate neutral density for polarization
       Neut(1)%Updens = Neut(1)%Updens + Clock%Dt * Src
       Neut(2)%Updens = Neut(2)%Updens + 2.d0 * Clock%Dt * Src
       !**** Energy conservation Diagnostic
       diag(7)%EnLoss = diag(7)%EnLoss + Clock%Dt * Src * abs(ion(1)%En-ion(2)%En)
       !***************** Diagnostic for relative importance of reactions (m-3/s)
       diag(7)%Tx(1) = Src
       !***************
    END IF
    !**** Rate calcul for adaptative time-step
    ratx = eta * meta(0)%Ni**2
    IF (ratx .GT. maxR) maxR = ratx
    !*************************************
    IF (mod(iter,2000)==0) THEN
       WRITE(919,"(A)") "Rates & density (s-1 | m-3 s-1 | m-3) in He+ + 2He(1S) --> He2+ + He(1S)"
       WRITE(919,"(A, 3ES12.3)") ion(1)%Name, ratx, Src, ion(1)%Ni
    END IF

    !**** cf. Belmonte 2007 : η (cm3 s-1)
    !**** He2+ + He(1S) --> He+ + 2He(1S)
    eta = 1.40d-06 * 1.d-6 * exp(-(ion(1)%En-ion(2)%En)*qok / Tp) / Tp**0.67
    Src = eta * ion(2)%Ni * meta(0)%Ni
    IF (ion(2)%Ni .GT. 0.d0) THEN
       ion(1)%Updens = ion(1)%Updens + Clock%Dt * Src
       ion(2)%Updens = ion(2)%Updens - Clock%Dt * Src
       !**** UpDate neutral density for polarization
       Neut(1)%Updens = Neut(1)%Updens + 2.d0 * Clock%Dt * Src
       Neut(2)%Updens = Neut(2)%Updens + Clock%Dt * Src
       !**** Energy conservation Diagnostic
       diag(7)%EnProd = diag(7)%EnProd + Clock%Dt * Src * abs(ion(1)%En-ion(2)%En)
       !***************** Diagnostic for relative importance of reactions (m-3/s)
       diag(11)%Tx(1) = Src
       !***************** 
    END IF
    !**** Rate calcul for adaptative time-step
    ratx = eta * meta(0)%Ni
    IF (ratx .GT. maxR) maxR = ratx
    !*************************************
    !*************************************
    IF (mod(iter,2000)==0) THEN
       WRITE(919,"(A)") "Rates & density (s-1 | m-3 s-1 | m-3) in He2+ + He(1S) --> He+ + 2He(1S)"
       WRITE(919,"(A,3ES12.3)") ion(2)%Name, eta*meta(0)%Ni, Src, ion(2)%Ni
    END IF
    !**** 3body collisions (Excimer creation) : He2*
    SELECT CASE (NumIon)
    CASE (3)
       !**** He(2P3) + 2He --> He2* + He (cm6 s-1)
       !**** rate from Koymen et al (Chem.Phys.Lett 168 5 1990)
       IF (meta(0)%Tp*qok.GT.300.d0) THEN
          excim = 1.6d-44 * meta(3)%Ni * meta(0)%Ni**2
       ELSE
          excim = (2.5d0 + 267.d0/(meta(0)%Tp*qok))*1d-44 * meta(3)%Ni * meta(0)%Ni**2
       END IF
       meta(3)%UpDens = meta(3)%UpDens - Clock%Dt * excim
       ion(NumIon)%UpDens  = ion(NumIon)%UpDens  + Clock%Dt * excim
       !**** UpDate neutral density for polarization
       Neut(1)%Updens = Neut(1)%Updens + Clock%Dt * excim
       Neut(2)%Updens = Neut(2)%Updens + 2.d0 * Clock%Dt * excim
       !**** Rate calcul for adaptative time-step
       ratx = excim / meta(3)%Ni
       IF (ratx .GT. maxR) maxR = ratx
       IF (mod(iter,2000)==0) THEN
          WRITE(919,"(A)") "Rates & density (s-1 | m-3 s-1 | m-3) in He(2P3) + 2He --> He2* + He"
          WRITE(919,"(A,3ES12.3)") meta(3)%Name, ratx, excim, meta(3)%Ni
       END IF
       !***************** Diagnostic for relative importance of reactions (m-3/s)
       diag(12)%Tx(1) = excim
       !*************** Diagnostic for metastable and 2^3P rates (s-1)
       diag(12)%OutM2 = excim / meta(3)%Ni
       !**** He2* + He --> He(2P3) + 2He (cm3 s-1)
       !**** rate from Belmonte et al (J.Phys.D:Appl.Phys 40 7343 2007)
       excim = 3.6d-20 * ion(3)%Ni * meta(0)%Ni
       meta(3)%UpDens = meta(3)%UpDens + Clock%Dt * excim
       ion(3)%UpDens  = ion(3)%UpDens  - Clock%Dt * excim
       !**** UpDate neutral density for polarization
       Neut(1)%Updens = Neut(1)%Updens + 2.d0 * Clock%Dt * excim
       Neut(2)%Updens = Neut(2)%Updens + Clock%Dt * excim
       !**** Rate calcul for adaptative time-step
       ratx = excim / meta(3)%Ni
       IF (ratx .GT. maxR) maxR = ratx
       IF (mod(iter,2000)==0) THEN
          WRITE(919,"(A)") "Rates & density (s-1 | m-3 s-1 | m-3) in He2* + He --> He(2P3) + 2He"
          WRITE(919,"(A,3ES12.3)") ion(3)%Name,ratx,excim, ion(3)%Ni
       END IF
       !***************** Diagnostic for relative importance of reactions (m-3/s)
       diag(13)%Tx(1) = excim
       !*************** Diagnostic for metastable and 2^3P rates (s-1)
       diag(13)%InM2 = excim
       !**** He(2S3) + 2He --> He2* + He (cm6 s-1)
       IF (meta(0)%Tp*qok.LE.750d0) THEN
          !**** rate from Koymen et al (Chem.Phys.Lett 168 5 1990)
          excim = Tp*(8.7d0*exp(-750.d0/Tp)+0.41d0*exp(-200/Tp))*1d-36*1d-12 &
               * meta(1)%Ni * meta(0)%Ni**2
       ELSE
          excim = 1.5d-46 *meta(1)%Ni * meta(0)%Ni**2
       END IF
       meta(1)%UpDens = meta(1)%UpDens - Clock%Dt * excim
       ion(3)%UpDens  = ion(3)%UpDens  + Clock%Dt * excim
       !**** Rate calcul for adaptative time-step
       ratx = excim / meta(1)%Ni
       IF (ratx .GT. maxR) maxR = ratx
       IF (mod(iter,2000)==0) THEN
          WRITE(919,"(A)") "Rates & density (s-1 | m-3 s-1 | m-3) in He(2S3) + 2He --> He2* + He"
          WRITE(919,"(A,3ES12.3)") meta(1)%Name,ratx,excim,meta(1)%Ni
       END IF
       !***************** Diagnostic for relative importance of reactions (m-3/s)
       diag(14)%Tx(1) = excim
       !*************** Diagnostic for metastable and 2^3S rates (s-1)
       diag(14)%OutM1 = excim / meta(1)%Ni
    END SELECT

    IF (mod(iter,2000)==0) CLOSE(919)
  END SUBROUTINE Conv_3Body

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Init_Recomb (sys, meta)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    !LOCAL
    INTEGER :: i, j, npts
    REAL(DOUBLE) :: Dx, Du, U
    REAL(DOUBLE), DIMENSION(152) :: SecRead, EnRead
    Dx = sys%Dx
    !**************************************
    write(*,"(2A)",advance="no") tabul, 'Reading Cross Section : [recomb_he.cs]'
    !**************************************
    SecRead = 0.d0 ; EnRead = 0.d0
    OPEN(UNIT=51,FILE='./datFile/recomb_he.cs',ACTION="READ",STATUS="OLD")
    do i = 1, 8
       READ(51,*)
    END do
    READ(51,*) Npts
    READ(51,*)(EnRead(i), i=1,Npts)
    READ(51,*) ; READ(51,*)
    READ(51,*)(SecRead(i), i=1,Npts)
    CLOSE(51)
    !**************************************
    write(*,"(A)") ' ......... Done'
    !**************************************

    !**************************************
    !**** Interpolat Cross-Sect recombinat
    DO i=1, sys%nx
       Du=0.d0
       U = IdU(i,Dx)
       meta(0)%SecRec(i) = SecRead(Npts) ! Because the recomb diss
                                         ! stop at 39.64 eV. We take
                                         ! sigma(u>39.64) = sigma(u=39.64)
       DO j = 1, Npts-1
          IF ( U == EnRead(j) ) meta(0)%SecRec(i) = SecRead(j)
          IF ( U .gt. EnRead(j) .and. U .lt. EnRead(j+1)) Then
             Du = EnRead(j+1) - EnRead(j)
             meta(0)%SecRec(i) = ((EnRead(j+1) - U)*SecRead(j) )/Du &
                  + ((U - EnRead(j))*SecRead(j+1) )/Du
          END IF
       END DO
    END DO
    meta(0)%SecRec(:) = meta(0)%SecRec(:) * 1d-20
    meta(0)%SecRec(sys%nx) = 0.d0

  END SUBROUTINE Init_Recomb

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
END MODULE MOD_RECOMB
