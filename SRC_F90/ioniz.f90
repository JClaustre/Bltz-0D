!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 08/07/2015
! Objctv: Ionization (direct & undirect) processes in Helium
! note  : 50 : Divide the residual energy between the scattered 
!              in & out electron (50/50)
!         100: The scattered-in electron keep all the residual
!              energy from the collision.
!         Subcycles are used (only for the ground state atom) 
!         in order to keep a greater time-step
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
    REAL(DOUBLE) :: Eij, chi, rchi, Dx, ionz, ratx
    !SubCYCLING VARIABLES
    REAL(DOUBLE) :: SubDt
    INTEGER :: SubCycl, l
    diag(2)%Tx = 0.d0
    Dx = sys%Dx ; prod=0.d0 ; loss = 0.d0

    DO i = 0, NumMeta
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
       coef1 = gama * meta(i)%Ni
       Eij = ion(1)%En - meta(i)%En ! ionization threshold
       chi = Eij/Dx ; ichi = int(chi) ; rchi = chi - ichi
       IF (rchi .LT. 0.d0 .OR. Eij .LT. 0.d0) then
          print*, 'probleme rchi<0 in [Ioniz]' ; STOP
       END IF
          
       DO l = 1, SubCycl
          ionz=0.d0

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
             Fi(k) = Fi(k) + SubDt * (prod-loss) * coef1 / sqrt(U(k))
          END DO
          ratx = ionz * Dx * gama
          if (ratx .GT. maxR) maxR = ratx
          diag(2)%Tx = diag(2)%Tx + ratx * meta(i)%Ni
          diag(2)%EnLoss = diag(2)%EnLoss + SubDt * ionz * coef1 * Dx*Eij
          meta(i)%UpDens = meta(i)%UpDens - SubDt * ionz * coef1 * Dx
          ion(1)%Updens  = ion(1)%Updens  + SubDt * ionz * coef1 * Dx
       END DO
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
    INTEGER :: i, k, kp, ichi, cas, nx
    REAL(DOUBLE) :: prod, loss, ratx
    REAL(DOUBLE) :: Eij, chi, rchi, Dx
    REAL(DOUBLE) :: Coef, coef1, cnst, Src
    !SubCYCLING VARIABLES
    REAL(DOUBLE) :: SubDt
    INTEGER :: SubCycl, l
    Dx = sys%Dx ; diag(2)%Tx = 0.d0
    nx = sys%nx

    cas = 1 ! if 0 then "Vidal case" | else "Matte case"
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
          IF (cas == 0) Eij = Eij + Dx*0.5d0
          chi = Eij/Dx ; ichi = int(chi) ; rchi = chi - ichi
          IF (rchi .LT. 0.d0 .OR. Eij .LT. 0.d0) then
             print*, 'probleme rchi<0 in [Ioniz]', i, meta(i)%En; STOP
          END IF

          DO k = 1, nx
             prod= 0.d0 ; loss= 0.d0 ; Coef = 1.d0
             kp = k + ichi

             IF (k == ichi+1) Coef = (1.d0-rchi)
             loss = Coef * U(k) * meta(i)%SecIon(1,k) * Fi(k)
             IF (kp   .LE. nx) prod = U(kp) * meta(i)%SecIon(1,kp)*Fi(kp) * (1.d0-rchi)
             IF (kp+1 .LE. nx) prod = prod + rchi * U(kp+1) * meta(i)%SecIon(1,kp+1) * Fi(kp+1)
             !**** Excited states balance
             IF (k .GE. ichi+1) Then
                coef = 1.d0
                if (k == ichi+1) coef = (1.d0-rchi) * (1.d0- (rchi/chi) )
                Src = Src + ( coef * U(k) * meta(i)%SecIon(1,k) * Fi(k) )
             END IF
             !**** UpDate Distribution Function
             Fi(k) = Fi(k) + SubDt * (prod-loss) * coef1 / dsqrt(U(k))
          END DO
          IF ( cas == 0 ) THEN
             Fi(1) = Fi(1) + SubDt * Src * coef1* cnst * Dx
             !**** Diagnostic
             diag(2)%EnLoss = diag(2)%EnLoss + SubDt * Src * coef1* Dx*(Eij-Dx*0.5d0)
          ELSE
             Fi(1) = Fi(1) + SubDt * Src * coef1* cnst * Dx * 3.d0 / 2.d0
             Fi(2) = Fi(2) - SubDt * Src * coef1* cnst * Dx / (2.d0 * dsqrt(3.d0))
             !**** Diagnostic
             diag(2)%EnLoss = diag(2)%EnLoss + SubDt * Src * coef1 * Dx * Eij
          END IF
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
  ! Second electron with 0 energy from excimer He2*
  SUBROUTINE Ioniz_Dimer100 (sys, ion, U, Fi)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(:), INTENT(INOUT) :: ion
    REAL(DOUBLE) , DIMENSION(:), INTENT(IN)    :: U
    REAL(DOUBLE) , DIMENSION(:), INTENT(INOUT) :: Fi
    !LOCAL
    INTEGER :: k, kp, km, ichi, case, Nion
    REAL(DOUBLE) :: prod, loss, rcmb, ionz
    REAL(DOUBLE) :: Eij, chi, rchi, Dx, br
    REAL(DOUBLE) :: Coef, coef1, coef2, cnst, Si, Sr
    REAL(DOUBLE), DIMENSION(sys%nx) :: Fo
    Dx = sys%Dx ; br = 0.3d0
    SELECT CASE (3)
    CASE (3) 
       Nion = 3
    END SELECT

    case = 0 ! if 0 then "Francois case" | else "J-P case"
    cnst = dsqrt(2.d0/Dx**3.d0)

    Si = 0.d0 ; Sr = 0.d0
    coef1 = gama * ion(Nion)%Ni
    coef2 = gama * ion(2)%Ni
    Eij = ion(2)%En - ion(Nion)%En ! ionization threshold
    IF (case == 0) Eij = Eij + Dx*0.5d0
    chi = Eij/Dx ; ichi = int(chi) ; rchi = chi - ichi

    DO k = 1, sys%nx
       prod= 0.d0 ; loss= 0.d0 ; Coef = 1.d0
       kp = k + ichi
       km = k - ichi
       Fo(k) = Fi(k)
       !**** Ionization process
       IF (k == ichi+1) Coef = (1.d0-rchi)
       loss = Coef * U(k) * ion(Nion)%SecIon(1,k) * Fi(k)
       IF (kp   .LE. sys%nx) prod = U(kp) * ion(Nion)%SecIon(1,kp)*Fi(kp) * (1.d0-rchi)
       IF (kp+1 .LE. sys%nx) prod = prod + rchi * U(kp+1) * ion(Nion)%SecIon(1,kp+1) * Fi(kp+1)
       ionz = (prod - loss) * coef1 / dsqrt(U(k))

       !**** 3-body recombination
       prod= 0.d0 ; loss= 0.d0
       loss = U(k) * ion(Nion)%SecIon(2,k) * Fi(k)
       IF (km   .GT. 0) prod = U(km) * ion(Nion)%SecIon(2,km)* (1.d0-rchi) * Fo(km)
       IF (km-1 .GT. 0) prod = prod + rchi * U(km-1) * ion(Nion)%SecIon(2,km-1) * Fo(km-1)
       rcmb = (prod - loss) * coef2 / dsqrt(U(k))

       !**** Excited states balance
       IF (k .GE. ichi+1) Then
          coef = 1.d0
          if (k == ichi+1) coef = (1.d0-rchi) * (1.d0- (rchi/chi) )
          Si = Si + ( coef * U(k) * ion(Nion)%SecIon(1,k) * Fi(k)*coef1 )
       END IF
       IF (k .LE. (sys%nx-ichi-1) ) Sr = Sr + ( U(k) * ion(Nion)%SecIon(2,k) * Fi(k)*coef2 )

       !**** UpDate Electron Distribution Function
       Fi(k) = Fi(k) + Clock%Dt * (ionz + rcmb)
    END DO

    IF ( case == 0 ) THEN
       Fi(1) = Fi(1) + Clock%Dt * (Si - Sr) * cnst * Dx
       !**** Diagnostic
       diag(13)%EnLoss = diag(13)%EnLoss + Clock%Dt * Si * Dx*(Eij-Dx*0.5d0)
       diag(13)%EnProd = diag(13)%EnProd + Clock%Dt * Sr * Dx*(Eij-Dx*0.5d0)
    ELSE
       Fi(1) = Fi(1) + Clock%Dt * (Si-Sr) * cnst * Dx * 3.d0 / 2.d0
       Fi(2) = Fi(2) - Clock%Dt * (Si-Sr) * cnst * Dx / (2.d0 * dsqrt(3.d0))
       !**** Diagnostic
       diag(13)%EnLoss = diag(13)%EnLoss + Clock%Dt * Si * Dx * Eij
       diag(13)%EnProd = diag(13)%EnProd + Clock%Dt * Sr * Dx * Eij
    END IF
    !**** br == branching ratio
    ion(Nion)%Updens = ion(Nion)%Updens + Clock%Dt * ((1.-br)*Sr-Si) * Dx
    meta(1)%Updens   = meta(1)%Updens   + Clock%Dt * br*Sr * Dx
    ion(1)%Updens    = ion(1)%Updens    + Clock%Dt * (Si-Sr) * Dx
  END SUBROUTINE Ioniz_Dimer100
  !***********************************************************************

  !***********************************************************************
  SUBROUTINE Init_ioniz (sys, meta)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    !LOCAL
    INTEGER :: i, j, k, l, npts, ichi
    REAL(DOUBLE) :: Dx, Du, Eij, U, rchi
    REAL(DOUBLE) :: A, B, C
    REAL(DOUBLE), DIMENSION(200) :: SecRead, EnRead
    REAL(DOUBLE), DIMENSION(6,0:18) :: fit
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

    !**** Fitting Coefficients for ionization in Supp data Santos 
    !**** (J.Phys.D: Appl.Phys 47 2014)
    DO l = 1, 6
       READ(51,*) Npts
       READ(51,*)(fit(l,i), i=0,Npts-1)
       READ(51,*) ; READ(51,*) ; READ(51,*); READ(51,*)
    END DO
    CLOSE(51)
    !**************************************
    write(*,"(A)") ' ......... Done'
    !**************************************

    !**** He(n>2,l,s) --> He+
    DO i=3, NumMeta
       Eij=ion(1)%En-meta(i)%En
       IF(Eij .GT. 0.d0) THEN
          DO k=1,sys%Nx
             U=IdU(k,Dx)/Eij
             IF(U.GT.1.d0) meta(i)%SecIon(1,k) = 3.49d0*(Ry/Eij)**2 * &
                  ((U-1.d0)/U**2.14d0) * log(1.25d0*U)
             meta(i)%SecIon(1,k) = meta(i)%SecIon(1,k) * 1d-20
          END DO
          meta(i)%SecIon(1,sys%nx) = 0.d0
       END IF
    END DO

    !**** Ionization(1)/Recombination(2) from Dimer He2*
    !**** He2* + e- <--> He2+ + 2e-
    SELECT CASE (NumIon)
    CASE (3) 
       Eij=ion(2)%En-ion(NumIon)%En
       A = 9.93844d-15 ; B = 9.8416d-1 ; C = 1.292d-02
       DO k = 1, sys%nx
          U=IdU(k,Dx)
          ion(NumIon)%SecIon(1,k) = 1d-04*( A*log(U/Eij)/(Eij*U) ) * (1.d0-B*exp(C)*exp(-C*U/Eij))
          IF (U .LE. Eij) ion(NumIon)%SecIon(1,k) = 0.d0
       END DO
       !**** Recomb/Ioniz cross-section relation
       ichi = int(Eij/Dx) ; rchi = (Eij/Dx) - ichi
       DO k = 1, sys%nx
          Du=IdU(k,Dx)/Eij
          if(k .LE. sys%nx-ichi) ion(NumIon)%SecIon(2,k) = (sqrt(Pi)/4.d0)*(Du/(Du+1.d0))&
               * ( (1.0d0-rchi) * ion(NumIon)%SecIon(1,k+ichi) )
          if(k .LE. sys%nx-ichi-1) ion(NumIon)%SecIon(2,k) = ion(NumIon)%SecIon(2,k) &
               + (sqrt(Pi)/4.d0)*(Du/(Du+1.d0))* ( rchi * ion(NumIon)%SecIon(1,k+ichi+1) )
       END DO
       ion(NumIon)%SecIon(1,sys%nx) = 0.d0
       ion(NumIon)%SecIon(2,sys%nx) = 0.d0
    END SELECT

  END SUBROUTINE Init_ioniz
  !***********************************************************************

END MODULE MOD_IONIZ
