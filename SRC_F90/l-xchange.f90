!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 08/07/2015
! Objctv: l-exchange and s-exchange processes in He
! note  : Calculation of an equilibrium solution due to their 
!         high collision frequency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_XCHANGE
  USE F90_KIND  
  USE MOD_PARAM
  IMPLICIT NONE

CONTAINS
  !***********************************************************************
  !**** No time dependent solution ***************************************
  !**** We use : dni/dt = -Kij * Ni*Nhe + Kji * Nii*Nhe
  !**** with dni/dt = 0 and N = Ni + Nii
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE l_change (meta, Kij)
    !INTENT
    TYPE(Species), DIMENSION(0:) , INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:,:), INTENT(IN)  :: Kij
    !LOCAL
    INTEGER :: i, j
    REAL(DOUBLE) :: Coef1, Coef2, Eij
    !********************
    Coef1= 0.d0 ; coef2 = 0.d0
    !**** l-change atomic reations
    DO i = 1 , NumMeta
       DO j = 1, NumMeta

          IF (Kij(i,j) .NE. 0.d0) THEN
             Eij = (meta(i)%En-meta(j)%En)
             IF (Eij .GT. 0.d0) THEN
                Coef1 = Kij(j,i) / Kij(i,j)
                Coef2 = meta(i)%Ni + Meta(j)%Ni

                meta(i)%Ni = (Coef1 * Coef2) / (1.d0 + Coef1)
                meta(j)%Ni = Coef2 / (1.d0 + Coef1)
             END IF
          END IF
       END DO
    END DO

    !**** s-change atomic reations
    DO i = 1, NumMeta
       if (meta(i)%Nl == 3 .and. meta(i)%Ns == 3) THEN
          Coef1 = meta(i)%Ni + meta(i+1)%Ni
          Coef2 = Kij(i+1,i)/Kij(i,i+1)

          Meta(i+1)%Ni = (1.0d0 / (1.0d0 + 1.d0/Coef2)) * Coef1/ Coef2
          Meta(i)%Ni   = (1.0d0 / (1.0d0 + Coef2)) * Coef2 * Coef1
       END if
    END DO
  END SUBROUTINE l_change
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !**** rate (cm-3 s-1) too high!
  SUBROUTINE l_change_old (meta, Kij)
    !INTENT
    TYPE(Species), DIMENSION(0:) , INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:,:), INTENT(IN)  :: Kij
    !LOCAL
    INTEGER :: i, j
    REAL(DOUBLE) :: Coef1, Coef2, Eij
    Coef1= 0.d0 ; coef2 = 0.d0

    !**** l-change atomic reations
    DO i = 1 , NumMeta
       DO j = 1, NumMeta

          IF (Kij(i,j) .NE. 0.d0) THEN
             Eij = (meta(i)%En-meta(j)%En)
             IF (Eij .GT. 0.d0) THEN
                Coef1 = Kij(i,j) * meta(0)%Ni * meta(i)%Ni
                Coef2 = Kij(j,i) * meta(0)%Ni * meta(j)%Ni
                meta(i)%Ni = meta(i)%Ni + Clock%Dt * (Coef2 - Coef1)
                meta(j)%Ni = meta(j)%Ni + Clock%Dt * (Coef1 - Coef2)
             END IF
          END IF

       END DO
    END DO

    !**** s-change atomic reations
    DO i = 1, NumMeta
       if (meta(i)%Nl == 3 .and. meta(i)%Ns == 3) THEN
          Coef1 = meta(i)%Ni + meta(i+1)%Ni
          Coef2 = Kij(i+1,i)/Kij(i,i+1)
          
          meta(i)%Ni   = (1.0d0 / (1.0d0 + Coef2)) * Coef2 * Coef1
          meta(i+1)%Ni = meta(i)%Ni / Coef2
       END if
    END DO

  END SUBROUTINE l_change_old
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Init_lchange (meta, Kij)
    !INTENT
    TYPE(Species), DIMENSION(0:) , INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:,:), INTENT(INOUT) :: Kij
    !LOCAL
    INTEGER :: i
    REAL(DOUBLE) :: Ratio, coeff
    
    !**********************************************
    !**** Cross Section Calculation ***************
    !**** l-change atomic process *****************
    !**** He(n,l,s)+He <--> He(n,l',s)+He *********
    !**** s-change atomic process *****************
    !**** He(n,l,s)+He <--> He(n,l,s')+He *********
    !**********************************************
    Kij=0.d0
    DO i = 1, NumMeta
       coeff= 7.76d-10 * 1d-6 / meta(i)%Deg
       !**** n=3D3-->3P3
       IF (i == 8) THEN
          Ratio = ( (meta(i)%En - meta(i-1)%En) / meta(0)%Tp)
          Kij(i,i-1)= coeff * (Ratio)**(-0.29d0)
          IF (Ratio .GT. 1.d0) Kij(i,i-1) = Kij(i,i-1) * exp(-1.28d0*Ratio)
          Kij(i-1,i) = Kij(i,i-1) * meta(i)%Deg * exp(-Ratio) / meta(i-1)%Deg
          !print*, meta(i)%Name, meta(i-1)%Name, K_ij(i,i-1), K_ij(i-1,i)
       END IF
       !**** n=3P1-->3D1
       IF (i == 10) THEN
          Ratio = ( (meta(i)%En - meta(i-1)%En) / meta(0)%Tp)
          Kij(i,i-1)= coeff * (Ratio)**(-0.29d0)
          IF (Ratio .GT. 1.d0) K_ij(i,i-1) = Kij(i,i-1) * exp(-1.28d0*Ratio)
          Kij(i-1,i) = Kij(i,i-1) * meta(i)%Deg * exp(-Ratio) / meta(i-1)%Deg
       END IF
       if (meta(i)%Nn .GT. 3) THEN
           !**** nP3-->nS3 
          IF (meta(i)%Nl == 1 .and. meta(i)%Ns == 3) THEN
             Ratio = ( (meta(i)%En - meta(i-2)%En) / meta(0)%Tp)
             Kij(i,i-2)= coeff * (Ratio)**(-0.29d0)
             IF (Ratio .GT. 1.d0) Kij(i,i-2) = Kij(i,i-2) * exp(-1.02d0*Ratio)
             Kij(i-2,i) = Kij(i,i-2) * meta(i)%Deg * exp(-Ratio) / meta(i-2)%Deg
          END IF
          !**** nD3-->nS3 
          IF (meta(i)%Nl == 2 .and. meta(i)%Ns == 3) THEN
             Ratio = ( (meta(i)%En - meta(i-3)%En) / meta(0)%Tp)
             Kij(i, i-3)= coeff * (Ratio)**(-0.29d0)
             IF (Ratio .GT. 1.d0) Kij(i,i-3) = Kij(i,i-3) * exp(-1.d0*Ratio)
             Kij(i-3,i) = Kij(i,i-3) * meta(i)%Deg * exp(-Ratio) / meta(i-3)%Deg
          END IF
          !**** nF,G,H,I3-->nP3 
          IF (meta(i)%Nl == 3 .and. meta(i)%Ns == 3) THEN
             Ratio = ( (meta(i)%En - meta(i-3)%En) / meta(0)%Tp)
             Kij(i,i-3)= coeff * (Ratio)**(-0.29d0)
             IF (Ratio .GT. 1.d0) Kij(i,i-3) = Kij(i,i-3) * exp(-.71d0*Ratio)
             Kij(i-3,i) = Kij(i,i-3) * meta(i)%Deg * exp(-Ratio) / meta(i-3)%Deg
          END IF
          !**** nF,G,H,I3-->nS3 
          IF (meta(i)%Nl == 3 .and. meta(i)%Ns == 3) THEN
             Ratio = ( (meta(i)%En - meta(i-5)%En) / meta(0)%Tp)
             Kij(i,i-5)= coeff * (Ratio)**(-0.29d0)
             IF (Ratio .GT. 1.d0) Kij(i,i-5) = Kij(i,i-5) * exp(-1.d0*Ratio)
             Kij(i-5,i) = Kij(i,i-5) * meta(i)%Deg * exp(-Ratio) / meta(i-5)%Deg
          END IF
          !**** nF,G,H,I1-->nS1
          IF (meta(i)%Nl == 3 .and. meta(i)%Ns == 1) THEN
             Ratio = ( (meta(i)%En - meta(i-5)%En) / meta(0)%Tp)
             Kij(i,i-5)= coeff * (Ratio)**(-0.29d0)
             IF (Ratio .GT. 1.d0) Kij(i,i-5) = Kij(i,i-5) * exp(-1.07d0*Ratio)
             Kij(i-5,i) = Kij(i,i-5) * meta(i)%Deg * exp(-Ratio) / meta(i-5)%Deg
          END IF
          !**** nD3-->nP3 | nP1-->nF,G,H,I1
          IF ( (meta(i)%Nl == 2 .and. meta(i)%Ns == 3)  .or. &
               (meta(i)%Nl == 1 .and. meta(i)%Ns == 1) ) THEN
             Ratio = ( (meta(i)%En - meta(i-1)%En) / meta(0)%Tp)
             Kij(i,i-1)= coeff * (Ratio)**(-0.29d0)
             Kij(i-1,i) = Kij(i,i-1) * meta(i)%Deg * exp(-Ratio) / meta(i-1)%Deg
          END IF
          !**** nF,G,H,I3-->nD3 | nF,G,H,I1-->nD1
          IF ( (meta(i)%Nl == 3 .and. meta(i)%Ns == 3)  .or. &
               (meta(i)%Nl == 3 .and. meta(i)%Ns == 1) ) THEN
             Ratio = ( (meta(i)%En - meta(i-2)%En) / meta(0)%Tp)
             Kij(i,i-2)= coeff * (Ratio)**(-0.29d0)
             Kij(i-2,i) = Kij(i,i-2) * meta(i)%Deg * exp(-Ratio) / meta(i-2)%Deg
          END IF
          !**** nP1-->nD1 
          IF (meta(i)%Nl == 1 .and. meta(i)%Ns == 1) THEN
             Ratio = ( (meta(i)%En - meta(i-3)%En) / meta(0)%Tp)
             Kij(i,i-3)= coeff * (Ratio)**(-0.29d0)
             Kij(i-3,i) = Kij(i,i-3) * meta(i)%Deg * exp(-Ratio) / meta(i-3)%Deg
          END IF
          !**** nP1-->nS1 
          IF (meta(i)%Nl == 1 .and. meta(i)%Ns == 1) THEN
             Ratio = ( (meta(i)%En - meta(i-6)%En) / meta(0)%Tp)
             Kij(i,i-6)= coeff * (Ratio)**(-0.29d0)
             IF (Ratio .GT. 1.d0) Kij(i,i-6) = Kij(i,i-6) * exp(-1.05d0*Ratio)
             Kij(i-6,i) = Kij(i,i-6) * meta(i)%Deg * exp(-Ratio) / meta(i-6)%Deg
          END IF
          !**** nD1-->nS1 
          IF (meta(i)%Nl == 2 .and. meta(i)%Ns == 1) THEN
             Ratio = ( (meta(i)%En - meta(i-3)%En) / meta(0)%Tp)
             Kij(i,i-3)= coeff * (Ratio)**(-0.29d0)
             IF (Ratio .GT. 1.d0) Kij(i,i-3) = Kij(i,i-3) * exp(-1.28d0*Ratio)
             Kij(i-3,i) = Kij(i,i-3) * meta(i)%Deg * exp(-Ratio) / meta(i-3)%Deg
          END IF
       END if
    END DO

    !**********************************************
    !**** s-exchange atomic process ***************
    !**** He(n,l,s)+He <--> He(n,l,s')+He *********
    !**** only for n1,3( F,G,H,I ) ****************
    !**********************************************
    DO i = 1, NumMeta
       if (meta(i)%Nl == 3 .and. meta(i)%Ns == 3) THEN
          Kij(i,i+1) = 1.0d-4 * 1d-6
          Kij(i+1,i) = 1.0d-4 * 1d-6 * (meta(i)%Deg / meta(i+1)%Deg) ! * exp(0.d0)
       END if
    END DO

  END SUBROUTINE Init_lchange
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!

END MODULE MOD_XCHANGE
