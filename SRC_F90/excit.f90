!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 08/07/2015
! Objctv: Excitation & SuperElastic processes in He
! note  : data's and empiric formula in 
!         Luis Alves et al (doi:10.1088/0022-3727/25/12/007)
!         M Santos et al (doi:10.1088/0022-3727/47/26/265201)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_EXCIT
  USE F90_KIND  
  USE MOD_PARAM
  IMPLICIT NONE

CONTAINS

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Exc_Begin(sys, meta, U, Fi, diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
    !LOCAL
    INTEGER :: i, j, k, kp, km, ichi
    REAL(DOUBLE) :: n, coef, coef1, coef2
    REAL(DOUBLE) :: C_Exc, C_Dxc, prod, loss
    REAL(DOUBLE) :: chi, rchi, E_ij
    REAL(DOUBLE) :: Sx, Sd, steady, En1, En2, Pn1, Pn2
    REAL(DOUBLE), DIMENSION(sys%nx) :: Fo
    !********************
    n = sys%Dx
    !********************

    !********************************************************
    DO i = 0, NumMeta-1
       coef1 = meta(i)%Ni * gama
       !********************************************************
       DO j = i+1, NumMeta
          coef2 = meta(j)%Ni * gama
          !********************************************************
          Sx = 0.d0 ; Sd = 0.d0
          E_ij = meta(j)%En - meta(i)%En

          IF (E_ij .NE. 0.d0 .AND. meta(i)%SecExc(j,1) .NE. 111.d0) THEN
             chi = E_ij/sys%Dx ; ichi = int(chi) ; rchi = chi - ichi
             IF (rchi .LT. 0.d0 .OR. E_ij .LT. 0.d0) then
                print*, 'probleme rchi<0 in [Inelastic]' ; STOP
             END IF

             DO k = 1, sys%nx
                kp = k + ichi
                km = k - ichi
                Fo(k) = Fi(k)
                !**** Excitation
                prod= 0.d0 ; loss= 0.d0 ; Coef = 1.d0
                IF (k == ichi+1) Coef = (1.d0-rchi)
                loss = Coef * U(k) * meta(i)%SecExc(j,k) * Fi(k)
                IF (kp   .LE. sys%nx) prod = U(kp) * meta(i)%SecExc(j,kp)*Fi(kp) * (1.d0-rchi)
                IF (kp+1 .LE. sys%nx) prod = prod + rchi * U(kp+1) * meta(i)%SecExc(j,kp+1) * Fi(kp+1)
                C_Exc = coef1 * (prod - loss) / sqrt(U(k))
                !**************************** 

                prod= 0.d0 ; loss= 0.d0 ; coef = 1.d0
                !**** De-Excitation
                loss = U(k) * meta(j)%SecExc(i,k) * Fi(k)
                IF (km   .GT. 0) prod = U(km) * meta(j)%SecExc(i,km)* (1.d0-rchi) * Fo(km)
                IF (km-1 .GT. 0) prod = prod + rchi * U(km-1) * meta(j)%SecExc(i,km-1) * Fo(km-1)
                C_Dxc = coef2 * (prod - loss) / sqrt(U(k))
                !**************************** 

                !**** Excited states balance
                IF (k .GE. ichi+1) Then
                   coef = 1.d0
                   if (k == ichi+1) coef = (1.d0-rchi) * (1.d0- (rchi/chi) )
                   Sx = Sx + ( coef * U(k) * meta(i)%SecExc(j,k) * Fi(k)* gama* n )
                END IF
                IF (k .LE. (sys%nx-ichi-1) ) Sd = Sd + ( U(k) * meta(j)%SecExc(i,k) * Fi(k)* gama*n)
                !**************************** 

                !**** UpDate EEDF
                Fi(k) = Fi(k) + Clock%Dt * ( C_Exc + C_Dxc )
                !**************************** 
             END DO
             !**** Diagnostic
             IF (i == 0) THEN
                diag(1)%EnProd(NumMeta+Numion+1) = diag(1)%EnProd(NumMeta+Numion+1) + &
                     Clock%Dt * Sd*meta(j)%Ni* E_ij
                diag(1)%EnLoss(NumMeta+Numion+1) = diag(1)%EnLoss(NumMeta+Numion+1) + &
                     Clock%Dt * Sx*meta(i)%Ni* E_ij
             ELSE
                diag(1)%EnProd(i) = diag(1)%EnProd(i) + Clock%Dt * Sd*meta(j)%Ni * E_ij
                diag(1)%EnLoss(i) = diag(1)%EnLoss(i) + Clock%Dt * Sx*meta(i)%Ni * E_ij
                diag(1)%DnProd(j) = diag(1)%DnProd(j) + Clock%Dt * Sx*meta(i)%Ni
                diag(1)%DnLoss(j) = diag(1)%DnLoss(j) + Clock%Dt * Sd*meta(j)%Ni
             END IF
             !*****************
             !**** UpDate Density
             meta(i)%UpDens = meta(i)%UpDens + Clock%Dt*(Sd*meta(j)%Ni - Sx*meta(i)%Ni)
             meta(j)%UpDens = meta(j)%UpDens + Clock%Dt*(Sx*meta(i)%Ni - Sd*meta(j)%Ni)

             if (Sd .GT. MaxR) MaxR = Sd
             IF (Sx .GT. MaxR) MaxR = Sx
          END IF
       END DO
    END DO
  END SUBROUTINE Exc_Begin

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Exc_Equil(sys, meta, U, Fi, diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
    !LOCAL
    INTEGER :: i, j, k, kp, km, ichi
    REAL(DOUBLE) :: n, coef, coef1, coef2
    REAL(DOUBLE) :: C_Exc, C_Dxc, prod, loss
    REAL(DOUBLE) :: chi, rchi, E_ij
    REAL(DOUBLE) :: Sx, Sd, steady
    REAL(DOUBLE), DIMENSION(sys%nx) :: Fo
    !SubCYCLING VARIABLES
    REAL(DOUBLE) :: SubDt, SubRt
    !Equili VaRIABLES
    REAL(DOUBLE) :: Ni, Nj, Nexpl, Rmx, Rmd
    !********************
    SubDt = Clock%Dt
    SubRt = 2.d-10
    n = sys%Dx
    !********************

    !********************************************************
    DO i = 0, NumMeta-1
       coef1 = meta(i)%Ni * gama
       !********************************************************
       DO j = i+1, NumMeta
          coef2 = meta(j)%Ni * gama
          !********************************************************
          Sx = 0.d0 ; Sd = 0.d0
          E_ij = meta(j)%En - meta(i)%En

          IF (E_ij .NE. 0.d0 .AND. meta(i)%SecExc(j,1) .NE. 111.d0) THEN
             chi = E_ij/sys%Dx ; ichi = int(chi) ; rchi = chi - ichi
             IF (rchi .LT. 0.d0 .OR. E_ij .LT. 0.d0) then
                print*, 'probleme rchi<0 in [Inelastic]' ; STOP
             END IF

             DO k = 1, sys%nx
                !**** Excited states balance
                IF (k .GE. ichi+1) Then
                   coef = 1.d0
                   if (k == ichi+1) coef = (1.d0-rchi) * (1.d0- (rchi/chi) )
                   Sx = Sx + ( coef * U(k) * meta(i)%SecExc(j,k) * Fi(k)* gama* n )
                END IF
                IF (k .LE. (sys%nx-ichi-1) ) Sd = Sd + ( U(k) * meta(j)%SecExc(i,k) * Fi(k)* gama*n)
             END DO

             IF (1./Sx .LE. SubRt .or. 1./Sd .LE. SubRt) THEN
                !**** Equilibre 
                Ni = meta(i)%Ni ; Nj = meta(j)%Ni
                steady = meta(i)%Ni + meta(j)%Ni
                meta(i)%Ni = ( 1.d0 / (1.d0 + Sd/Sx) ) * (steady * Sd ) / Sx
                meta(j)%Ni = ( 1.d0 / (1.d0 + Sx/Sd) ) * (steady * Sx ) / Sd
                Nexpl = 0.d0 ; Nexpl = Ni + SubDt * (Sd*Nj - Sx*Ni)
                Rmx = (meta(i)%Ni-Ni) / (Nexpl-Ni)
                IF ((Nexpl-Ni) .EQ. 0.d0) Rmx = 0.d0
                Nexpl = 0.d0 ; Nexpl = Nj + SubDt * (Sx*Ni - Sd*Nj)
                Rmd = (meta(j)%Ni-Nj) / (Nexpl-Nj)
                IF ((Nexpl-Nj) .EQ. 0.d0) Rmd = 0.d0
             ELSE
                Rmx = 1.d0 ; Rmd = 1.d0
                if (Sd .GT. MaxR) MaxR = Sd
                IF (Sx .GT. MaxR) MaxR = Sx
                meta(j)%UpDens = meta(j)%UpDens + SubDt * (Sx*meta(i)%Ni-Sd*meta(j)%Ni)
                meta(i)%UpDens = meta(i)%UpDens + SubDt * (Sd*meta(j)%Ni-Sx*meta(i)%Ni)
                !*****************
                IF(i.GT.0) diag(1)%DnProd(i) = diag(1)%DnProd(i) + SubDt * Sd*meta(j)%Ni
                IF(i.GT.0) diag(1)%DnLoss(i) = diag(1)%DnLoss(i) + SubDt * Sx*meta(i)%Ni
             END IF

             DO k = 1, sys%nx
                kp = k + ichi
                km = k - ichi
                Fo(k) = Fi(k)
                !**** Excitation
                prod= 0.d0 ; loss= 0.d0 ; Coef = 1.d0
                IF (k == ichi+1) Coef = (1.d0-rchi)
                loss = Coef * U(k) * meta(i)%SecExc(j,k) * Fi(k)
                IF (kp   .LE. sys%nx) prod = U(kp) * meta(i)%SecExc(j,kp)*Fi(kp) * (1.d0-rchi)
                IF (kp+1 .LE. sys%nx) prod = prod + rchi * U(kp+1) * meta(i)%SecExc(j,kp+1) * Fi(kp+1)
                C_Exc = coef1 * (prod - loss) / sqrt(U(k))

                prod= 0.d0 ; loss= 0.d0 ; coef = 1.d0
                !**** De-Excitation
                loss = U(k) * meta(j)%SecExc(i,k) * Fi(k)
                IF (km   .GT. 0) prod = U(km) * meta(j)%SecExc(i,km)* (1.d0-rchi) * Fo(km)
                IF (km-1 .GT. 0) prod = prod + rchi * U(km-1) * meta(j)%SecExc(i,km-1) * Fo(km-1)
                C_Dxc = coef2 * (prod - loss) / sqrt(U(k))

                !**** UpDate EEDF
                Fi(k) = Fi(k) + SubDt * ( C_Exc * Rmx + C_Dxc * Rmd )
                !**************************** 
             END DO
             IF (i == 0) THEN
                diag(1)%EnProd(NumMeta+Numion+1) = diag(1)%EnProd(NumMeta+Numion+1) + &
                     SubDt * Sd*meta(j)%Ni* E_ij * Rmd
                diag(1)%EnLoss(NumMeta+Numion+1) = diag(1)%EnLoss(NumMeta+Numion+1) + &
                     SubDt * Sx*meta(i)%Ni* E_ij * Rmx
             ELSE
                diag(1)%EnProd(i) = diag(1)%EnProd(i) + SubDt * Sd*meta(j)%Ni * E_ij * Rmd
                diag(1)%EnLoss(i) = diag(1)%EnLoss(i) + SubDt * Sx*meta(i)%Ni * E_ij * Rmx
                diag(1)%DnProd(j) = diag(1)%DnProd(j) + SubDt * Sx*meta(i)%Ni
                diag(1)%DnLoss(j) = diag(1)%DnLoss(j) + SubDt * Sd*meta(j)%Ni
             END IF
             !*****************
          END IF
       END DO
    END DO
    !********************
  END SUBROUTINE Exc_Equil

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Exc_Impli(sys, meta, U, Fi, diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
    !LOCAL
    INTEGER :: i, j, k, kp, km, ichi
    REAL(DOUBLE) :: n, coef, coef1, coef2
    REAL(DOUBLE) :: C_Exc, C_Dxc, prod, loss
    REAL(DOUBLE) :: chi, rchi, E_ij
    REAL(DOUBLE) :: Sx, Sd, steady
    REAL(DOUBLE), DIMENSION(sys%nx) :: Fo
    !SubCYCLING VARIABLES
    REAL(DOUBLE) :: SubDt
    !Implicit Density VARIABLES
    REAL(DOUBLE) :: Ni, Nj, Nexpl, Rmx, Rmd
    !********************
    SubDt = Clock%Dt
    n = sys%Dx
    !********************

    !********************************************************
    DO i = 0, NumMeta-1
       coef1 = meta(i)%Ni * gama
       !********************************************************
       DO j = i+1, NumMeta
          coef2 = meta(j)%Ni * gama
          !********************************************************
          Sx = 0.d0 ; Sd = 0.d0
          E_ij = meta(j)%En - meta(i)%En

          IF (E_ij .NE. 0.d0 .AND. meta(i)%SecExc(j,1) .NE. 111.d0) THEN
             chi = E_ij/sys%Dx ; ichi = int(chi) ; rchi = chi - ichi
             IF (rchi .LT. 0.d0 .OR. E_ij .LT. 0.d0) then
                print*, 'probleme rchi<0 in [Inelastic]' ; STOP
             END IF

             !**** Excited states balance
             DO k = 1, sys%nx
                IF (k .GE. ichi+1) Then
                   coef = 1.d0
                   if (k == ichi+1) coef = (1.d0-rchi) * (1.d0- (rchi/chi) )
                   Sx = Sx + ( coef * U(k) * meta(i)%SecExc(j,k) * Fi(k)* gama* n )
                END IF
                IF (k .LE. (sys%nx-ichi-1) ) Sd = Sd + ( U(k) * meta(j)%SecExc(i,k) * Fi(k)* gama*n)
             END DO

             !**** Implicit Density
             Ni = meta(i)%Ni ; Nj = meta(j)%Ni
             meta(i)%Ni = ( 1.d0 / (1.d0 + SubDt*Sx) ) * (Ni + SubDt*meta(j)%Ni*Sd)
             meta(j)%Ni = ( 1.d0 / (1.d0 + SubDt*Sd) ) * (Nj + SubDt*meta(i)%Ni*Sx)
             Nexpl = 0.d0 ; Nexpl = Ni + SubDt * (Sd*Nj - Sx*Ni)
             Rmx = (meta(i)%Ni-Ni) / (Nexpl-Ni)
             IF ((Nexpl-Ni) .EQ. 0.d0 ) Rmx = 0.d0
             Nexpl = 0.d0 ; Nexpl = Nj + SubDt * (Sx*Ni - Sd*Nj)
             Rmd = (meta(j)%Ni-Nj) / (Nexpl-Nj)
             IF ((Nexpl-Nj) .EQ. 0.d0 ) Rmd = 0.d0

             DO k = 1, sys%nx
                kp = k + ichi
                km = k - ichi
                Fo(k) = Fi(k)
                !**** Excitation
                prod= 0.d0 ; loss= 0.d0 ; Coef = 1.d0
                IF (k == ichi+1) Coef = (1.d0-rchi)
                loss = Coef * U(k) * meta(i)%SecExc(j,k) * Fi(k)
                IF (kp   .LE. sys%nx) prod = U(kp) * meta(i)%SecExc(j,kp)*Fi(kp) * (1.d0-rchi)
                IF (kp+1 .LE. sys%nx) prod = prod + rchi * U(kp+1) * meta(i)%SecExc(j,kp+1) * Fi(kp+1)
                C_Exc = coef1 * (prod - loss) / sqrt(U(k))

                prod= 0.d0 ; loss= 0.d0 ; coef = 1.d0
                !**** De-Excitation
                loss = U(k) * meta(j)%SecExc(i,k) * Fi(k)
                IF (km   .GT. 0) prod = U(km) * meta(j)%SecExc(i,km)* (1.d0-rchi) * Fo(km)
                IF (km-1 .GT. 0) prod = prod + rchi * U(km-1) * meta(j)%SecExc(i,km-1) * Fo(km-1)
                C_Dxc = coef2 * (prod - loss) / sqrt(U(k))

                !**** UpDate EEDF
                Fi(k) = Fi(k) + SubDt * ( C_Exc * Rmx + C_Dxc * Rmd )
                !**************************** 
             END DO
             !**** Diagnostic
             IF (i == 0) THEN
                diag(1)%EnProd(NumMeta+Numion+1) = diag(1)%EnProd(NumMeta+Numion+1) + &
                     SubDt * Sd*meta(j)%Ni* E_ij * Rmd
                diag(1)%EnLoss(NumMeta+Numion+1) = diag(1)%EnLoss(NumMeta+Numion+1) + &
                     SubDt * Sx*meta(i)%Ni* E_ij * Rmx
             ELSE
                diag(1)%EnProd(i) = diag(1)%EnProd(i) + SubDt * Sd*meta(j)%Ni * E_ij * Rmd
                diag(1)%EnLoss(i) = diag(1)%EnLoss(i) + SubDt * Sx*meta(i)%Ni * E_ij * Rmx
                diag(1)%DnProd(j) = diag(1)%DnProd(j) + SubDt * Sx*meta(i)%Ni
                diag(1)%DnLoss(j) = diag(1)%DnLoss(j) + SubDt * Sd*meta(j)%Ni
             END IF
             !*****************
          END IF
       END DO
    END DO
    !********************
  END SUBROUTINE Exc_Impli

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Init_ExcDxc(sys, meta, Fosc, Q, A)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(0:,0:), INTENT(IN) :: Fosc !OSCILLATOR STRENGTH
    REAL(DOUBLE) , DIMENSION(0:,0:), INTENT(IN) :: Q, A !Coeff
    !LOCAL
    INTEGER :: i, j, k, ichi
    REAL(DOUBLE) :: rchi, Eij, Dx, Du, coef
    REAL(DOUBLE) :: Rp, G, alpha, pi_alpha2

    Dx = sys%Dx ; pi_alpha2 = 0.879735d0 * 1d-20
    !**********************************************
    !**** Cross Section Calculation ***************
    !**** Excitation/De-excitation process ********
    !**** He(n,l,s)+ e <--> He(n',l',s')+ e *******
    !**********************************************
    DO i=0, NumMeta-1
       DO j=i+1, NumMeta
          Eij=Meta(j)%En-Meta(i)%En
          ichi = int(Eij/Dx) ; rchi = (Eij/Dx) - ichi
          Rp=0.5d0 * (Fosc(i,j)*Ry/Eij)**(-0.7d0)
          !**** 1S1<-->2s3
          IF (i == 0 .and. j == 1) THEN 
             DO k=1, sys%nx
                Du=IdU(k,Dx)/Eij
                meta(i)%SecExc(j,k)= pi_alpha2*(3.76d-2*(Du-1.d0)/Du**3 + &
                     6.65d-1*(Du-1.d0)/Du**5)
                IF(Du.LE.1.d0) meta(i)%SecExc(j,k) = 0.d0
             END DO
          END IF
          !**** 1S1-->n3P | 2S1-->n3P | 1S1,2S1-->n3 | 2S3-->n1P
          IF ( (i == 0 .and. (meta(j)%Nn .GE. 4 .and. meta(j)%Nl == 1 .and. meta(j)%Ns == 3)) .or.&
               (i == 2 .and. (meta(j)%Nn .GE. 2 .and. meta(j)%Nl == 1 .and. meta(j)%Ns == 3)) .or.&
               ((i == 0 .or. i == 2) .and. (meta(j)%Nn .GE. 3 .and. meta(j)%Nl .NE. 1 .and. meta(j)%Ns == 3) ) .or.&
               (i == 1 .and. (meta(j)%Nn .GE. 2 .and. meta(j)%Ns == 1)) ) THEN 
             !print*, meta(i)%Name, meta(j)%Name
             DO k=1, sys%nx
                Du=IdU(k,Dx)/Eij
                meta(i)%SecExc(j,k)= pi_alpha2 * Q(i,j) * (Du**2-1.d0) / Du**A(i,j)
                IF(Du.LE.1.d0) meta(i)%SecExc(j,k) = 0.d0
             END DO
          END IF
          !**** 1S1-->2S1 | 1S1,2S1-->n1 | 2S3-->n3S,D,F,G,H,I
          IF ( (i == 0 .and. j == 2) .or.&
               ((i == 0 .or. i == 2) .and. (meta(j)%Nn .GE. 3 .and. meta(j)%Nl .NE. 1 .and. meta(j)%Ns == 1) ) .or.&
               (i == 1 .and. (meta(j)%Nn .GE. 3 .and. meta(j)%Nl .NE. 1 .and. meta(j)%Ns == 3)) ) THEN 
             !print*, meta(i)%Name, meta(j)%Name
             DO k=1, sys%nx
                Du=IdU(k,Dx)/Eij
                meta(i)%SecExc(j,k)= pi_alpha2 * Q(i,j) * (Du-1.d0) / Du**A(i,j)
                IF(Du.LE.1.d0) meta(i)%SecExc(j,k) = 0.d0
             END DO
          END IF
          !**** 1S1-->n1P | 2S1-->n1P | 2S3-->n3P
          IF ( (i == 0 .and. (meta(j)%Nn .GE. 5 .and. meta(j)%Nl == 1 .and. meta(j)%Ns == 1) ) .or.&
               (i == 2 .and. (meta(j)%Nn .GE. 2 .and. meta(j)%Nl == 1 .and. meta(j)%Ns == 1) ) .or.&
               (i == 1 .and. (meta(j)%Nn .GE. 2 .and. meta(j)%Nl == 1 .and. meta(j)%Ns == 3)) ) THEN 
             !print*, meta(i)%Name, meta(j)%Name
             DO k=1, sys%nx
                Du=IdU(k,Dx)/Eij
                meta(i)%SecExc(j,k)= pi_alpha2 * Q(i,j) * (Ry/Eij)**2 * (Du-1.d0)/Du**2 &
                     * log(1.25d0*Du)
                IF(Du.LE.1.d0) meta(i)%SecExc(j,k) = 0.d0
             END DO
          END IF
          !**** 2P3-->n3S,D | 2P1-->n1S,D
          IF ( (i == 3 .and. (meta(j)%Nn .GE. 3 .and. (meta(j)%Nl == 0 .or. meta(j)%Nl == 2) &
               .and. meta(j)%Ns == 3) ) .or.&
               (i == 4 .and. (meta(j)%Nn .GE. 3 .and. (meta(j)%Nl == 0 .or. meta(j)%Nl == 2) &
               .and. meta(j)%Ns == 1) ) ) THEN 
             Alpha = 1
             IF (i == 3) Alpha = 3
             DO k=1, sys%nx
                Du=IdU(k,Dx)/Eij
                G = (1.d0-exp(-RP*(Du-1.d0)))*log(Du+0.2d0)/Du
                meta(i)%SecExc(j,k)= pi_alpha2 * 4.d0*Alpha * (Ry/Eij)**2 * Fosc(i,j)*G
                IF(Du.LE.1.d0) meta(i)%SecExc(j,k) = 0.d0
             END DO
          END IF
          !**** Tout le reste (n,l,s)-->(n',l±1,s)
          IF ( (meta(i)%Nn .GE. 3 .and. meta(j)%Nl .LE. meta(i)%Nl+1 .and. meta(j)%Ns == meta(i)%Ns) ) THEN 
             !print*, meta(i)%Name, meta(j)%Name
             DO k=1, sys%nx
                Du=IdU(k,Dx)/Eij
                G = (1.d0-exp(-RP*(Du-1.d0)))*log(Du+0.2d0)/Du
                meta(i)%SecExc(j,k)= pi_alpha2 * 4.d0 * (Ry/Eij)**2 * Fosc(i,j)*G
                IF(Du.LE.1.d0) meta(i)%SecExc(j,k) = 0.d0
             END DO
          END IF
       END DO
    END DO
    !**********************************************
    !**** De-excitation cross-Section using *******
    !**** Klein-Rosseland Relation ****************
    !**********************************************
    coef = 0.d0
    DO i=0, NumMeta-1
       DO j=i+1, NumMeta
          Eij=Meta(j)%En-Meta(i)%En
          coef = real(meta(i)%Deg / meta(j)%Deg)
          IF(Eij .LE. 0.d0) GOTO 366
          ichi = int(Eij/Dx) ; rchi = (Eij/Dx) - ichi
          DO k=1,sys%Nx
             Du=IdU(k,Dx)/Eij
             if(k .LE. sys%nx-ichi) meta(j)%SecExc(i,k) = coef*((Du)/(Du+1.d0))&
                  * ( (1.0d0-rchi) * meta(i)%SecExc(j,k+ichi) )
             if(k .LE. sys%nx-ichi-1) meta(j)%SecExc(i,k) = meta(j)%SecExc(i,k) &
                  + coef*((Du)/(Du+1.d0))* ( rchi * meta(i)%SecExc(j,k+ichi+1) )
          END DO
          meta(i)%SecExc(j,sys%nx) = 0.d0
          meta(j)%SecExc(i,sys%nx) = 0.d0
366    END DO
    END DO

    DO i=0, NumMeta-1
       DO j=i+1, NumMeta
          if (sum(meta(i)%SecExc(j,:)) == 0.d0)THEN
             meta(i)%SecExc(j,1) = 111.d0
          END if
       END DO
    END DO

  END SUBROUTINE Init_ExcDxc
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
END MODULE MOD_EXCIT
