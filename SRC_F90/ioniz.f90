!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 08/07/2015
! Objctv: Ionization processes in He
! note  : data's and analytic formula in 
!         Luis Alves et al (doi:10.1088/0022-3727/25/12/007)
!         M Santos et al (doi:10.1088/0022-3727/47/26/265201)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_IONIZ
  USE F90_KIND  
  USE MOD_PARAM
  IMPLICIT NONE

CONTAINS

  !***********************************************************************
  ! Second electron with 1/2 energy of primary one
  SUBROUTINE Ioniz_50 (sys, meta, U, Fi, diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
    !LOCAL
    INTEGER :: i, k, kp, ichi
    REAL(DOUBLE) :: prod, loss, coef, coef1
    REAL(DOUBLE) :: Eij, chi, rchi, Dx, ionz
    Dx = sys%Dx ; prod=0.d0 ; loss = 0.d0

    DO i = 0, NumMeta
       ionz=0.d0
       coef1 = gama * meta(i)%Ni
       Eij = ion(1)%En - meta(i)%En ! ionization threshold
       chi = Eij/Dx ; ichi = int(chi) ; rchi = chi - ichi
       IF (rchi .LT. 0.d0 .OR. Eij .LT. 0.d0) then
          print*, 'probleme rchi<0 in [Ioniz]' ; STOP
       END IF

       DO k = 1, sys%nx
          prod=0.d0 ; loss = 0.d0
          kp = 2*k + ichi
          loss = meta(i)%SecIon(1,k) * Fi(k) * U(k)

          IF ((kp-2 .GT. 0) .and. (kp-2 .LE. sys%nx) ) &
               prod = prod + (1.d0-rchi)*0.5d0 * meta(i)%SecIon(1,kp-2) * Fi(kp-2) * U(kp-2)
          IF (kp-1 .GT. 0 .and. kp-1 .LE. sys%nx) &
               prod = prod + (1.5d0-rchi) * meta(i)%SecIon(1,kp-1) * Fi(kp-1) * U(kp-1)
          IF (kp .LE. sys%nx) &
               prod = prod + 1.5d0 * meta(i)%SecIon(1,kp) * Fi(kp) * U(kp)
          IF (kp+1 .LE. sys%nx) &
               prod = prod + (0.5d0+rchi) * meta(i)%SecIon(1,kp+1) * Fi(kp+1) * U(kp+1)
          IF (kp+2 .LE. sys%nx) &
               prod = prod + rchi*0.5d0 * meta(i)%SecIon(1,kp+2) * Fi(kp+2) * U(kp+2)

          !**** Excited states balance
          IF (k .GE. ichi+1) then
             coef = 1.d0
             if (k == ichi+1) coef = (1.0d0-rchi)
             ionz = ionz + coef * Fi(k) * U(k) * meta(i)%SecIon(1,k)
          END IF
          !**** UpDate EEDF
          Fi(k) = Fi(k) + Clock%Dt * (prod-loss) * coef1 / sqrt(U(k))
       END DO
       diag(2)%EnProd(NumMeta+1) = diag(2)%EnProd(NumMeta+1) + Clock%Dt * ionz * coef1 * Dx*Eij
       diag(2)%DnProd(NumMeta+1) = diag(2)%DnProd(NumMeta+1) + Clock%Dt * ionz * coef1 * Dx
       IF(i.GT.0) diag(2)%DnLoss(i) = diag(2)%DnLoss(i) + Clock%Dt * ionz * coef1 * Dx
       meta(i)%UpDens = meta(i)%UpDens - Clock%Dt * ionz * coef1 * Dx
       ion(1)%Updens  = ion(1)%Updens  + Clock%Dt * ionz * coef1 * Dx
    END DO

  END SUBROUTINE Ioniz_50
  !***********************************************************************

  !***********************************************************************
  ! Second electron with 0 energy
  SUBROUTINE Ioniz_100 (sys, meta, U, Fi, diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
    !LOCAL
    INTEGER :: i, k, kp, ichi, case
    REAL(DOUBLE) :: prod, loss, ratx
    REAL(DOUBLE) :: Eij, chi, rchi, Dx
    REAL(DOUBLE) :: Coef, coef1, cnst, Src
    !SubCYCLING VARIABLES
    REAL(DOUBLE) :: SubDt
    INTEGER :: SubCycl, l
    Dx = sys%Dx ; diag(2)%Tx = 0.d0

    case = 0 ! if 0 then "Francois case" | else "J-P case"
    cnst = dsqrt(2.d0/Dx**3.d0)

    DO i = 0, NumMeta
       coef1 = gama * meta(i)%Ni
       !**** SubCycling for the direct ionization 
       IF (i .NE. 0) THEN
          SubCycl = 1
          SubDt = Clock%Dt
       ELSE
          SubDt = 1.0d-11
          IF (Clock%Dt .GT. SubDt) SubCycl = ceiling(Clock%Dt / SubDt)
          IF (Clock%Dt .LE. SubDt) SubCycl = 1
          IF (Clock%Dt .LE. SubDt) SubDt = Clock%Dt
       END IF
       !*******************************************
       DO l = 1, SubCycl
          Src = 0.d0
          Eij = ion(1)%En - meta(i)%En ! ionization threshold
          IF (case == 0) Eij = Eij + Dx*0.5d0
          chi = Eij/Dx ; ichi = int(chi) ; rchi = chi - ichi
          IF (rchi .LT. 0.d0 .OR. Eij .LT. 0.d0) then
             print*, 'probleme rchi<0 in [Ioniz]', i, meta(i)%En; STOP
          END IF

          DO k = 1, sys%nx
             prod= 0.d0 ; loss= 0.d0 ; Coef = 1.d0
             kp = k + ichi

             IF (k == ichi+1) Coef = (1.d0-rchi)
             loss = Coef * U(k) * meta(i)%SecIon(1,k) * Fi(k)
             IF (kp   .LE. sys%nx) prod = U(kp) * meta(i)%SecIon(1,kp)*Fi(kp) * (1.d0-rchi)
             IF (kp+1 .LE. sys%nx) prod = prod + rchi * U(kp+1) * meta(i)%SecIon(1,kp+1) * Fi(kp+1)
             !**** Excited states balance
             IF (k .GE. ichi+1) Then
                coef = 1.d0
                if (k == ichi+1) coef = (1.d0-rchi) * (1.d0- (rchi/chi) )
                Src = Src + ( coef * U(k) * meta(i)%SecIon(1,k) * Fi(k) )
             END IF
             !**** UpDate Distribution Function
             Fi(k) = Fi(k) + SubDt * (prod-loss) * coef1 / dsqrt(U(k))
          END DO
          IF ( case == 0 ) THEN
             Fi(1) = Fi(1) + SubDt * Src * coef1* cnst * Dx
          ELSE
             Fi(1) = Fi(1) + SubDt * Src * coef1* cnst * Dx * 3.d0 / 2.d0
             Fi(2) = Fi(2) - SubDt * Src * coef1* cnst * Dx / (2.d0 * dsqrt(3.d0))
          END IF
          !**** Diagnostic
          diag(2)%EnProd(NumMeta+1) = diag(2)%EnProd(NumMeta+1) + SubDt * Src * coef1* Dx*(Eij -Dx*0.5d0) 
          diag(2)%DnProd(NumMeta+1) = diag(2)%DnProd(NumMeta+1) + SubDt * Src * coef1* Dx
          IF(i.GT.0) diag(2)%DnLoss(i) = diag(2)%DnLoss(i) + SubDt * Src * coef1* Dx
          meta(i)%Updens = meta(i)%Updens - SubDt * Src * coef1 * Dx
          ion(1)%Ni  = ion(1)%Ni  + SubDt * Src * coef1 * Dx

          ratx = Src * Dx * gama
          if (ratx .GT. maxR) maxR = ratx
          diag(2)%Tx = diag(2)%Tx + ratx * meta(i)%Ni
       END DO
    END DO
  END SUBROUTINE Ioniz_100
  !***********************************************************************
  !***********************************************************************
!  ! Second electron with 0 energy
!  SUBROUTINE Ioniz_OMP (sys, meta, U, Fi, diag)
!    !INTENT
!    TYPE(SysVar) , INTENT(IN) :: sys
!    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
!    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
!    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
!    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
!    !LOCAL
!    INTEGER :: i, k, kp, ichi, case
!    REAL(DOUBLE) :: prod, loss
!    REAL(DOUBLE) :: Eij, chi, rchi, Dx
!    REAL(DOUBLE) :: Coef, coef1, cnst, Src, moy, En1, En2
!    !SubCYCLING VARIABLES
!    REAL(DOUBLE) :: SubDt
!    INTEGER :: SubCycl
!    INTEGER :: l
!    !OMP VARIABLES
!    INTEGER :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS, OMP_GET_NUM_PROCS, NTHREADS, TID, PROCS
!    REAL(DOUBLE), DIMENSION(sys%nx) :: Fo
!
!    Dx = sys%Dx ; moy = 0.d0
!    En1 = 0.d0; En2 = 0.d0
!
!    case = 0 ! if 0 then "Francois case" | else "J-P case"
!    cnst = dsqrt(2.d0/Dx**3.d0)
!
!    PROCS = OMP_GET_NUM_PROCS()
!    CALL OMP_SET_NUM_THREADS( (PROCS/2) )
!    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k, TID, coef1, Eij, chi, ichi, rchi),&
!    !$OMP& PRIVATE(Src, kp, prod, loss, Coef, l, SubDt, SubCycl)
!    !***********************************************************
!    !$OMP DO
!    DO i = 0, 0
!       !**** SubCycling for the direct ionization 
!       SubDt = 1d-11
!       IF (Clock%Dt .GT. SubDt) SubCycl = ceiling(Clock%Dt / SubDt)
!       IF (Clock%Dt .LE. SubDt) SubCycl = 1
!       !*******************************************
!       DO l = 1, SubCycl
!          Src = 0.d0
!          coef1 =  gama * meta(i)%Ni
!          Eij = ion(1)%En - meta(i)%En ! ionization threshold
!          IF (case == 0) Eij = Eij + Dx*0.5d0
!          chi = Eij/Dx ; ichi = int(chi) ; rchi = chi - ichi
!          IF (rchi .LT. 0.d0 .OR. Eij .LT. 0.d0) then
!             print*, 'probleme rchi<0 in [Ioniz]', i, meta(i)%En; STOP
!          END IF
!
!          DO k = 1, sys%nx
!             En1 = En1 + (Fi(k) * U(k)**1.5d0 * Dx)
!             prod= 0.d0 ; loss= 0.d0 ; Coef = 1.d0
!             kp = (k + ichi)
!             IF ( k == (ichi+1) ) Coef = (1.d0-rchi)
!             loss = Coef * U(k) * meta(i)%SecIon(1,k) * Fi(k)
!             IF (kp     .LE. sys%nx) prod = U(kp) * meta(i)%SecIon(1,kp)*Fi(kp) * (1.d0-rchi)
!             IF ((kp+1) .LE. sys%nx) prod = prod + (rchi * U(kp+1) * meta(i)%SecIon(1,(kp+1))*Fi(kp+1) )
!             !**** Excited states balance
!             IF ( k .GE. (ichi+1) ) Then
!                coef = 1.d0
!                if ( k == (ichi+1) ) coef = (1.d0-rchi) * (1.d0- (rchi/chi) )
!                Src = Src + ( coef * U(k) * meta(i)%SecIon(1,k) * Fi(k) )
!             END IF
!             !**** UpDate Distribution Function
!             Fi(k) = Fi(k) + ( SubDt * (prod-loss) * coef1 / dsqrt(U(k)) )
!             En2 = En2 + Fi(k) * U(k)**1.5d0 * Dx
!          END DO
!
!          IF ( case == 0 ) THEN
!             Fi(1) = Fi(1) + ( SubDt * Src * coef1* cnst * Dx )
!          ELSE
!             Fi(1) = Fi(1) + ( SubDt * Src * coef1* cnst * Dx * 3.d0 / 2.d0 )
!             Fi(2) = Fi(2) - ( SubDt * Src * coef1* cnst * Dx / (2.d0 * dsqrt(3.d0)) )
!          END IF
!          diag(2)%EnLoss(1) = diag(2)%EnLoss(1) + abs(En1 - En2)
!          diag(2)%EnProd(NumMeta+1) = diag(2)%EnProd(NumMeta+1) + SubDt * Src * coef1* Dx*(Eij -Dx*0.5d0) 
!          diag(2)%DnProd(NumMeta+1) = diag(2)%DnProd(NumMeta+1) + SubDt * Src * coef1* Dx
!          IF(i.GT.0) diag(2)%DnLoss(i) = diag(2)%DnLoss(i) + SubDt * Src * coef1* Dx
!          meta(i)%Updens = meta(i)%Updens - (SubDt * Src * coef1 * Dx)
!          ion(1)%Updens  = ion(1)%Updens  + (SubDt * Src * coef1 * Dx)
!          moy = moy + Src * Dx * coef1
!       END DO
!    END DO
!    !$OMP END DO
!    diag(2)%Tx = moy
!
!    Fo = Fi
!    !***********************************************************
!    TID = OMP_GET_THREAD_NUM()
!    IF (TID .EQ. 0) NTHREADS = OMP_GET_NUM_THREADS()
!    !***********************************************************
!    !$OMP DO REDUCTION(+:Fi)
!    DO i = 1, NumMeta
!       Src = 0.d0
!       coef1 =  gama * meta(i)%Ni
!       Eij = ion(1)%En - meta(i)%En ! ionization threshold
!       IF (case == 0) Eij = Eij + Dx*0.5d0
!       chi = Eij/Dx ; ichi = int(chi) ; rchi = chi - ichi
!       IF (rchi .LT. 0.d0 .OR. Eij .LT. 0.d0) then
!          print*, 'probleme rchi<0 in [Ioniz]', i, meta(i)%En; STOP
!       END IF
!       
!       DO k = 1, sys%nx
!          En1 = En1 + (Fo(k) * U(k)**1.5d0 * Dx)
!          prod= 0.d0 ; loss= 0.d0 ; Coef = 1.d0
!          kp = (k + ichi)
!          IF ( k == (ichi+1) ) Coef = (1.d0-rchi)
!          loss = Coef * U(k) * meta(i)%SecIon(1,k) * Fo(k)
!          IF (kp     .LE. sys%nx) prod = U(kp) * meta(i)%SecIon(1,kp)*Fo(kp) * (1.d0-rchi)
!          IF ((kp+1) .LE. sys%nx) prod = prod + (rchi * U(kp+1) * meta(i)%SecIon(1,(kp+1))*Fo(kp+1) )
!          !**** Excited states balance
!          IF ( k .GE. (ichi+1) ) Then
!             coef = 1.d0
!             if ( k == (ichi+1) ) coef = (1.d0-rchi) * (1.d0- (rchi/chi) )
!             Src = Src + ( coef * U(k) * meta(i)%SecIon(1,k) * Fo(k) )
!          END IF
!          !**** UpDate Distribution Function
!          Fi(k) = Fi(k) + ( Clock%Dt * (prod-loss) * coef1 / dsqrt(U(k)) )
!          En2 = En2 + Fi(k) * U(k)**1.5d0 * Dx
!       END DO
!
!       IF ( case == 0 ) THEN
!          Fi(1) = Fi(1) + ( Clock%Dt * Src * coef1* cnst * Dx )
!       ELSE
!          Fi(1) = Fi(1) + ( Clock%Dt * Src * coef1* cnst * Dx * 3.d0 / 2.d0 )
!          Fi(2) = Fi(2) - ( Clock%Dt * Src * coef1* cnst * Dx / (2.d0 * dsqrt(3.d0)) )
!       END IF
!       diag(2)%EnLoss(1) = diag(2)%EnLoss(1) + abs(En1 - En2)
!       diag(2)%EnProd(NumMeta+1) = diag(2)%EnProd(NumMeta+1) + Clock%Dt * Src * coef1* Dx*(Eij -Dx*0.5d0) 
!       diag(2)%DnProd(NumMeta+1) = diag(2)%DnProd(NumMeta+1) + Clock%Dt * Src * coef1* Dx
!       IF(i.GT.0) diag(2)%DnLoss(i) = diag(2)%DnLoss(i) + Clock%Dt * Src * coef1* Dx
!       meta(i)%Updens = meta(i)%Updens - (Clock%Dt * Src * coef1 * Dx)
!       ion(1)%Updens  = ion(1)%Updens  + (Clock%Dt * Src * coef1 * Dx)
!       moy = moy + Src * Dx * coef1
!    END DO
!    !$OMP END DO
!    diag(2)%Tx = moy
!    !$OMP END PARALLEL
!  END SUBROUTINE Ioniz_OMP
!  !***********************************************************************

  !***********************************************************************
  ! Second electron with 0 energy from excimer He2*
  SUBROUTINE Ioniz_Excimer100 (sys, ion, U, Fi, diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(:), INTENT(INOUT) :: ion
    Type(Diagnos), DIMENSION(:), INTENT(INOUT) :: diag
    REAL(DOUBLE) , DIMENSION(:), INTENT(IN)    :: U
    REAL(DOUBLE) , DIMENSION(:), INTENT(INOUT) :: Fi
    !LOCAL
    INTEGER :: i, k, kp, ichi, case
    REAL(DOUBLE) :: prod, loss
    REAL(DOUBLE) :: Eij, chi, rchi, Dx
    REAL(DOUBLE) :: Coef, coef1, cnst, Src, moy
    Dx = sys%Dx ; prod=0.d0 ; loss = 0.d0 ; moy = 0.d0

    case = 0 ! if 0 then "Francois case" | else "J-P case"
    cnst = dsqrt(2.d0/Dx**3.d0)

    Src = 0.d0
    coef1 = gama * ion(3)%Ni
    Eij = ion(2)%En - ion(3)%En ! ionization threshold
    IF (case == 0) Eij = Eij + Dx*0.5d0
    chi = Eij/Dx ; ichi = int(chi) ; rchi = chi - ichi
    IF (rchi .LT. 0.d0 .OR. Eij .LT. 0.d0) then
       print*, 'probleme rchi<0 in [Ioniz]', i, meta(i)%En; STOP
    END IF

    DO k = 1, sys%nx
       prod= 0.d0 ; loss= 0.d0 ; Coef = 1.d0
       kp = k + ichi

       IF (k == ichi+1) Coef = (1.d0-rchi)
       loss = Coef * U(k) * ion(3)%SecIon(1,k) * Fi(k)
       IF (kp   .LE. sys%nx) prod = U(kp) * ion(3)%SecIon(1,kp)*Fi(kp) * (1.d0-rchi)
       IF (kp+1 .LE. sys%nx) prod = prod + rchi * U(kp+1) * ion(3)%SecIon(1,kp+1) * Fi(kp+1)
       !**** Excited states balance
       IF (k .GE. ichi+1) Then
          coef = 1.d0
          if (k == ichi+1) coef = (1.d0-rchi) * (1.d0- (rchi/chi) )
          Src = Src + ( coef * U(k) * ion(3)%SecIon(1,k) * Fi(k) )
       END IF
       !**** UpDate Distribution Function
       Fi(k) = Fi(k) + Clock%Dt * (prod-loss) * coef1 / dsqrt(U(k))
    END DO

    IF ( case == 0 ) THEN
       Fi(1) = Fi(1) + Clock%Dt * Src * coef1* cnst * Dx
    ELSE
       Fi(1) = Fi(1) + Clock%Dt * Src * coef1* cnst * Dx * 3.d0 / 2.d0
       Fi(2) = Fi(2) - Clock%Dt * Src * coef1* cnst * Dx / (2.d0 * dsqrt(3.d0))
    END IF
    ion(3)%Updens = ion(3)%Updens - Clock%Dt * Src * coef1 * Dx
    ion(1)%Updens  = ion(1)%Updens  + Clock%Dt * Src * coef1 * Dx
  END SUBROUTINE Ioniz_Excimer100
  !***********************************************************************

  !***********************************************************************
  SUBROUTINE Init_ioniz (sys, meta)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    !LOCAL
    INTEGER :: i, j, k, l, npts
    REAL(DOUBLE) :: Dx, Du, Eij, U
    REAL(DOUBLE), DIMENSION(200) :: SecRead, EnRead
    Dx = sys%Dx

    !**************************************
    write(*,"(2A)",advance="no") tabul, 'Reading Cross Section : [ioniz_he.cs] '
    !**************************************

    OPEN(UNIT=51,FILE='./datFile/ioniz_he.cs',ACTION="READ",STATUS="OLD")
    do i = 1, 8
       READ(51,*)
    END do

    DO l = 0, 2
       SecRead = 0.d0 ; EnRead = 0.d0
       READ(51,*) Npts
       READ(51,*)(EnRead(i), i=1,Npts)
       READ(51,*) ; READ(51,*)
       READ(51,*)(SecRead(i), i=1,Npts)

       !**** Interpolat Cross-Sect Direct Ioniz (He(1S1|2S3|2S1) --> He+)
       DO i=1, sys%nx
          Du=0.d0
          U = IdU(i,Dx)
          DO j = 1, Npts-1
             IF ( U == EnRead(j) ) meta(l)%SecIon(1,i) = 1d-20 * SecRead(j)
             IF ( U .gt. EnRead(j) .and. U .lt. EnRead(j+1)) Then
                Du = EnRead(j+1) - EnRead(j)
                meta(l)%SecIon(1,i) =  1d-20 * ( ((EnRead(j+1) - U)*SecRead(j) )/Du &
                     + ((U - EnRead(j))*SecRead(j+1) )/Du )
             END IF
          END DO
       END DO
       meta(l)%SecIon(1,sys%nx) = 0.d0
       !**************************************
       READ(51,*) ; READ(51,*) ; READ(51,*); READ(51,*)
    END DO
    CLOSE(51)
    !**************************************
    write(*,"(A)") ' ......... Done'
    !**************************************

    !**** He(n>2,l,s) --> He+
    DO i=3, NumMeta
       Eij=ion(1)%En-meta(i)%En
       IF(Eij .LE. 0.d0) GOTO 9905
       DO k=1,sys%Nx
          U=IdU(k,Dx)/Eij
          IF(U.GT.1.d0) meta(i)%SecIon(1,k) = 3.49d0*(Ry/Eij)**2 * &
               ((U-1.d0)/U**2.14d0) * log(1.25d0*U)
          meta(i)%SecIon(1,k) = meta(i)%SecIon(1,k) * 1d-20
       END DO
       meta(i)%SecIon(1,sys%nx) = 0.d0
9905 END DO
    !**** Ionization from Excimer He2*
    !**** He2* + e- --> He2+ + 2e-
    IF (NumIon == 3) THEN
       Eij=ion(2)%En-ion(3)%En
       IF (Eij .GT. 0.d0) THEN
          DO k = 1, sys%nx
             U=IdU(k,Dx)/Eij
             IF(U.GT.1.d0) ion(3)%SecIon(1,k) = 3.49d0*(Ry/Eij)**2 * &
               ((U-1.d0)/U**2.14d0) * log(1.25d0*U)
          ion(3)%SecIon(1,k) = ion(3)%SecIon(1,k) * 1d-20
          END DO
       END IF
    END IF
  END SUBROUTINE Init_ioniz
  !***********************************************************************

END MODULE MOD_IONIZ
