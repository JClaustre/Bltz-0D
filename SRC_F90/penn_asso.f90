!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 08/07/2015
! Objctv: Penning and Associative processes in He
! note  : Several associative rates are available 
!         (see below "init_asso" routine) and can be changed in
!         the "read_init.f90" file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_PENNASS
  USE F90_KIND  
  USE MOD_PARAM
  IMPLICIT NONE

CONTAINS

  !***********************************************************************
  SUBROUTINE Penn_Assoc (sys, meta, U, Fi, Diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)    :: U
    TYPE(Diagnos), DIMENSION(:) , INTENT(INOUT) :: Diag
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: Fi
    !LOCAL
    INTEGER :: i, j, k, l, ichi, Nion, Nmeta
    REAL(DOUBLE) :: asso, Penn, beta, Dx, Ndens
    REAL(DOUBLE) :: Eij, chi, rchi, coef1, coef2
    REAL(DOUBLE) :: ratx, Rate
    !*********************************
    beta = 2.9d-9 * 1d-6 * (meta(0)%Tp*qok/300.d0)**(-1.86d0)
    asso=0.d0 ; Penn=0.d0 ; coef1=0.d0 ; coef2=0.d0 ; Rate = 0.d0
    Diag(6)%Tx(1) = 0.d0
    Dx = sys%Dx
    Nmeta = 34
    If (NumMeta.LT.34) Nmeta = NumMeta

    !**** Associative process
    DO i = 5, Nmeta
       Eij = meta(i)%En - ion(2)%En ! associative threshold
       IF (Eij > 0.d0) THEN
          chi = Eij/Dx + 0.5d0 ; ichi = int(chi)
          IF (ichi == 0) ichi = 1
          rchi = ( Eij - U(ichi) ) / Dx
          asso = meta(0)%Ni*meta(i)%Ni * Sn(i)
          
          DO k = 1, sys%nx
             coef1 = (1.d0 - rchi) / (sqrt(U(ichi))*Dx)
             coef2 = rchi / (sqrt(U(ichi+1))*Dx)
             IF (k .NE. ichi  ) coef1 = 0.d0
             IF (k .NE. ichi+1) coef2 = 0.d0
             Fi(k) = Fi(k) + Clock%Dt * (asso * (coef1 + coef2))
          END DO
          !**** Particle balance (explicit --> meta)
          meta(i)%Updens = meta(i)%Updens - Clock%Dt * asso
          meta(0)%Updens = meta(0)%Updens - Clock%Dt * asso
          ion(2)%Updens  = ion(2)%Updens  + Clock%Dt * asso
          !**** Energy conservation Diagnostic
          Diag(6)%EnProd = Diag(6)%EnProd + Clock%Dt * asso * Eij
          !***************
          Diag(6)%SumTx = Diag(6)%SumTx + Clock%Dt * asso

          ratx = Sn(i)*meta(0)%Ni
          IF (ratx .GT. MaxR) MaxR = ratx
          !***************** Diagnostic for relative importance of reactions (m-3/s)
          IF (asso.GT.Rate) THEN
             Rate = asso
             Diag(6)%Tx(2) = real(i)
          END IF
          Diag(6)%Tx(1) = Diag(6)%Tx(1) + asso
       END IF
    END DO

!    !**** Penning process
!    Rate=0.d0 ; Diag(5)%Tx(:)=0.d0
!    diag(5)%OutM1 = 0.d0 ; diag(5)%OutM2 = 0.d0
!    !********************
!    DO i = 1, 3
!       DO j = 1, 3
!          
!          DO l = 1, 2
!             IF (l .EQ. 1) Penn  =  0.3d0* meta(i)%Ni*meta(j)%Ni * beta
!             IF (l .EQ. 2) Penn  =  0.7d0* meta(i)%Ni*meta(j)%Ni * beta
!             Eij = meta(j)%En + meta(i)%En - ion(l)%En ! Penning threshold
!             chi = Eij/Dx + 0.5d0 ; ichi = int(chi)
!             IF (ichi == 0) ichi = 1
!             rchi = ( Eij - U(ichi) ) / Dx
!
!             DO k = 1, sys%nx
!                coef1 = (1.d0 - rchi) / (sqrt(U(ichi))*Dx)
!                coef2 = rchi / (sqrt(U(ichi+1))*Dx)
!                IF (k .NE. ichi  ) coef1 = 0.d0
!                IF (k .NE. ichi+1) coef2 = 0.d0
!                Fi(k) = Fi(k) + Clock%Dt * (Penn * (coef1 + coef2))
!             END DO
!             ion(l)%Updens = ion(l)%Updens + Clock%Dt * Penn
!             !**** Energy conservation Diagnostic
!             diag(5)%EnProd = diag(5)%EnProd + Eij*Clock%Dt * Penn
!             !****************
!             diag(5)%SumTx = diag(5)%SumTx + Clock%Dt * Penn
!             !***************** Diagnostic for relative importance of reactions (m-3/s)
!             IF (Penn.GT.Rate) THEN
!                Rate = Penn
!                Diag(5)%Tx(2) = real(i) ; Diag(5)%Tx(3) = real(j)
!             END IF
!             Diag(5)%Tx(1) = Diag(5)%Tx(1) + Penn
!
!          END DO
!          !**** Update population
!          Penn = meta(i)%Ni*meta(j)%Ni * beta
!          meta(i)%Updens = meta(i)%Updens - Clock%Dt * Penn
!          meta(j)%Updens = meta(j)%Updens - Clock%Dt * Penn
!          ratx = Penn / meta(i)%Ni
!          IF (ratx .GT. MaxR) MaxR = ratx
!          !*************** Diagnostic for metastable and 2^3P rates (s-1)
!          if(i==1) THEN
!             diag(5)%OutM1 = diag(5)%OutM1 + Penn/meta(i)%Ni
!          END if
!          if(j==1) THEN
!             diag(5)%OutM1 = diag(5)%OutM1 + Penn/meta(j)%Ni
!          END if
!          !****************
!          if(i==3) diag(5)%OutM2 = diag(5)%OutM2 + Penn/meta(i)%Ni
!          if(j==3) diag(5)%OutM2 = diag(5)%OutM2 + Penn/meta(j)%Ni
!
!       END DO
!    END DO
!
!    ! Involving Dimer Penning processes 
!    SELECT CASE (NumIon)
!    CASE (3)
!       Rate = 0.d0 ; Diag(4)%Tx(:) = 0.d0
!       diag(4)%OutM1=0.d0 ; diag(4)%OutM2=0.d0
!       Nion = 3
!       ! #1 He* + He2* --> He+ + He + e
!       !               --> He2+ + He + e
!       ! #2 He2* + He2* --> He+ + He + e  (case with i==0)
!       !                --> He2+ + He + e
!       DO i = 0, 3
!          beta = 2.5d-09 * 1d-6
!          IF (i == 0) beta = 1.5d-09* 1d-6
!          If (i == 0) Ndens = ion(Nion)%Ni
!          IF (i .GT. 0) Ndens = meta(i)%Ni
!          DO l = 1, 2
!             IF (l .EQ. 1) Penn  =  0.15d0* Ndens*ion(Nion)%Ni * beta
!             IF (l .EQ. 2) Penn  =  0.85d0* Ndens*ion(Nion)%Ni * beta
!             Eij = ion(Nion)%En + meta(i)%En - ion(l)%En ! Penning threshold
!             IF (i == 0) Eij = 2d0*ion(Nion)%En - ion(l)%En 
!             chi = Eij/Dx + 0.5d0 ; ichi = int(chi)
!             IF (ichi == 0) ichi = 1
!             rchi = ( Eij - U(ichi) ) / Dx
!
!             DO k = 1, sys%nx
!                coef1 = (1.d0 - rchi) / (sqrt(U(ichi))*Dx)
!                coef2 = rchi / (sqrt(U(ichi+1))*Dx)
!                IF (k .NE. ichi  ) coef1 = 0.d0
!                IF (k .NE. ichi+1) coef2 = 0.d0
!                Fi(k) = Fi(k) + Clock%Dt * (Penn * (coef1 + coef2))
!             END DO
!             ion(l)%Updens = ion(l)%Updens + Clock%Dt * Penn
!             !**** Energy conservation Diagnostic
!             diag(4)%EnProd = diag(4)%EnProd + Eij*Clock%Dt * Penn
!             diag(4)%SumTx = diag(4)%SumTx + Clock%Dt * Penn
!             !***************** Diagnostic for relative importance of reactions (m-3/s)
!             IF (Penn.GT.Rate) THEN
!                Rate = Penn             
!                Diag(4)%Tx(2) = real(i)
!             END IF
!             Diag(4)%Tx(1) = Diag(4)%Tx(1) + Penn
!
!             !*************** Diagnostic for metastable and 2^3P rates (cm-3 s-1)
!             IF (i==1) THEN
!                !*************** Diagnostic for metastable and 2^3P rates (s-1)
!                diag(4)%OutM1    = diag(4)%OutM1 + Penn/meta(i)%Ni
!             END IF
!             IF (i==3) THEN
!                diag(4)%OutM2    = diag(4)%OutM2 + Penn/meta(i)%Ni
!             END IF
!             !****************
!          END DO
!          !**** Update population
!          Penn = Ndens * ion(Nion)%Ni * beta
!          IF (i.GT.0) meta(i)%Updens   = meta(i)%Updens   - Clock%Dt * Penn
!          If (i == 0) ion(Nion)%Updens = ion(Nion)%Updens - Clock%Dt * Penn
!          ion(Nion)%Updens = ion(Nion)%Updens - Clock%Dt * Penn
!          ratx = Penn / Ndens
!          IF (ratx .GT. MaxR) MaxR = ratx
!
!       END DO
!
!    END SELECT

  END SUBROUTINE Penn_Assoc
  !***********************************************************************

  SUBROUTINE Init_Asso(Sn, n)
    !INTENT
    INTEGER, INTENT(IN) :: n
    REAL(DOUBLE), DIMENSION(:), INTENT(OUT) :: Sn
    Sn = 0.d0

    Select Case (n)
    Case (0) ! Alves data's
       !********************************************
       ! Table 4. in Alves 92.
       ! Associative ionization rate coef. (cm3 s-1)
       !********************************************
       Sn(8)  = 3.14d-11 * 1.d-6
       Sn(9)  = 1.96d-10 * 1.d-6
       Sn(10) = 1.31d-11 * 1.d-6
       Sn(11) = 2.75d-10 * 1.d-6
       Sn(12) = 2.48d-10 * 1.d-6
       Sn(13) = 3.14d-10 * 1.d-6
       Sn(14) = 3.40d-10 * 1.d-6
       Sn(15) = 5.90d-10 * 1.d-6
       Sn(16) = 3.30d-10 * 1.d-6
       Sn(17) = 3.30d-10 * 1.d-6
       Sn(18) = 1.70d-10 * 1.d-6
       Sn(19) = 3.00d-10 * 1.d-6
       Sn(20) = 2.70d-10 * 1.d-6
       Sn(21) = 3.45d-10 * 1.d-6
       Sn(22) = 3.70d-10 * 1.d-6
       Sn(23) = 6.50d-10 * 1.d-6
       Sn(24) = 7.30d-10 * 1.d-6
       Sn(25) = 7.30d-10 * 1.d-6
       Sn(26) = 1.90d-10 * 1.d-6
       Sn(27) = 3.20d-10 * 1.d-6
       Sn(28) = 2.90d-10 * 1.d-6
       Sn(29) = 3.60d-10 * 1.d-6
       Sn(30) = 3.90d-10 * 1.d-6
       Sn(31) = 6.80d-10 * 1.d-6
       Sn(32) = 8.20d-10 * 1.d-6
       Sn(33) = 8.20d-10 * 1.d-6
       Sn(34) = 2.00d-10 * 1.d-6
       !************************************
    Case (1) ! Santos data's
       Sn(5)  = 2.90d-11 * 1.d-6
       Sn(6)  = 2.90d-10 * 1.d-6
       Sn(7)  = 3.40d-11 * 1.d-6
       Sn(8)  = 3.60d-11 * 1.d-6
       Sn(9)  = 2.30d-10 * 1.d-6
       Sn(10) = 8.60d-12 * 1.d-6
       Sn(11) = 3.80d-11 * 1.d-6
       Sn(12) = 3.40d-11 * 1.d-6
       Sn(13) = 7.80d-11 * 1.d-6
       Sn(14) = 1.50d-10 * 1.d-6
       Sn(15) = 2.50d-10 * 1.d-6
       Sn(16) = 1.10d-10 * 1.d-6
       Sn(17) = 2.00d-10 * 1.d-6
       Sn(18) = 4.20d-11 * 1.d-6
       Sn(19) = 4.10d-11 * 1.d-6
       Sn(20) = 3.70d-11 * 1.d-6
       Sn(21) = 8.60d-11 * 1.d-6
       Sn(22) = 1.63d-10 * 1.d-6
       Sn(23) = 2.75d-10 * 1.d-6
       Sn(24) = 2.16d-10 * 1.d-6
       Sn(25) = 3.92d-10 * 1.d-6
       Sn(26) = 4.70d-11 * 1.d-6
       Sn(27) = 4.40d-11 * 1.d-6
       Sn(28) = 4.00d-11 * 1.d-6
       Sn(29) = 8.90d-11 * 1.d-6
       Sn(30) = 1.72d-10 * 1.d-6
       Sn(31) = 2.88d-10 * 1.d-6
       Sn(32) = 3.07d-10 * 1.d-6
       Sn(33) = 5.58d-10 * 1.d-6
       Sn(34) = 4.90d-11 * 1.d-6 
    Case (2)
       Sn(5)  = 1.00d-13 * 1.d-6
       Sn(6)  = 8.00d-10 * 1.d-6
       Sn(7)  = 2.50d-11 * 1.d-6
       Sn(8)  = 1.00d-13 * 1.d-6
       Sn(9)  = 1.00d-13 * 1.d-6
       Sn(10) = 1.00d-11 * 1.d-6
       Sn(11) = 2.00d-10 * 1.d-6
       Sn(12) = 1.50d-09 * 1.d-6
       Sn(13) = 1.00d-13 * 1.d-6
       Sn(14) = 1.00d-13 * 1.d-6
       Sn(15) = 1.00d-13 * 1.d-6
       Sn(16) = 5.00d-12 * 1.d-6
       Sn(17) = 2.00d-11 * 1.d-6
       Sn(18) = 1.00d-10 * 1.d-6
       Sn(19) = 4.10d-11 * 1.d-6
       Sn(20) = 4.10d-11 * 1.d-6
       Sn(21) = 8.60d-11 * 1.d-6
       Sn(22) = 1.63d-10 * 1.d-6
       Sn(23) = 2.75d-10 * 1.d-6
       Sn(24) = 2.16d-10 * 1.d-6
       Sn(25) = 3.92d-10 * 1.d-6
       Sn(26) = 4.70d-11 * 1.d-6
       Sn(27) = 4.40d-11 * 1.d-6
       Sn(28) = 4.00d-11 * 1.d-6
       Sn(29) = 8.90d-11 * 1.d-6
       Sn(30) = 1.72d-10 * 1.d-6
       Sn(31) = 2.88d-10 * 1.d-6
       Sn(32) = 3.07d-10 * 1.d-6
       Sn(33) = 5.58d-10 * 1.d-6
       Sn(34) = 4.90d-11 * 1.d-6 
    Case Default
       write(*,"(A)") "Choose Alves (0), Santos (1) or Belmonte (2) data's!"
    END Select

  END SUBROUTINE Init_Asso

END MODULE MOD_PENNASS
