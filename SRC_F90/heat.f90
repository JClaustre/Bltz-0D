!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Charlotte Boukandou (modified by Jonathan Claustre)
! Date  : 14/07/2015
! Objctv: Elastic collisions, Fokker-Planck and Heating 
!         processes in Helium
! note  : data's and analytic formula in 
!         Luis Alves et al (doi:10.1088/0022-3727/25/12/007)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_CHAUF
  USE F90_KIND
  USE MOD_PARAM
  IMPLICIT NONE

CONTAINS

  !*********** SUBROUTINE FOKKER-PLANCK ***************!
  subroutine FP (sys,elec, f1, U)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), INTENT(IN) :: elec
    REAL(DOUBLE), DIMENSION(:), INTENT(INOUT) :: f1
    REAL(DOUBLE), DIMENSION(:), INTENT(IN)    :: U
    !LOCAL
    INTEGER :: i, l, nx
    REAL(DOUBLE) :: part
    REAL(DOUBLE) :: hdu, tt, nu, v3, cnst, cn, err
    REAL(DOUBLE) :: y00, yy1, yy2, xx, truc, small, LnC
    REAL(DOUBLE), DIMENSION(:) , ALLOCATABLE :: f1_new
    REAL(DOUBLE), DIMENSION(:) , ALLOCATABLE :: II,JJ
    REAL(DOUBLE), DIMENSION(:) , ALLOCATABLE :: utt, usq, A, B, C, D

    !**** Log Coulomb (cf. NRL formulary) (Ne--> cm-3 !)

    LnC = 23.d0 - log((elec%Ni*1d-6)**0.5d0 / elec%Tp**1.5d0)

    nx = sys%nx
    ALLOCATE ( II(0:nx), JJ(0:nx), f1_new(nx) )
    ALLOCATE ( utt(nx), usq(nx), A(nx), B(nx), C(nx), D(nx) )

    hdu = 0.5d0*sys%Dx
    tt  = 2.d0/3.d0
    v3 = ( (2.d0*qome))**(1.5d0)
    nu = (4.d0*pi*elec%Ni*LnC*qe**4) / ( (4.d0*pi*eps*me)**2 * v3)
    
    cnst = -(2.0d0*Clock%Dt * nu) / Sys%Dx
    !********Initialisation:******************************
    part = 0.d0
    do i=1,nx
       utt(i)= 0.5d0*((U(i)+hdu)**(1.5d0) + (U(i)-hdu)**1.5d0)  
       usq(i)= tt * ((U(i)+hdu)**(1.5d0) - (U(i)-hdu)**1.5d0) /sys%Dx
       !*****Calcul des quantites initiales
       part = part + f1(i) * U(i)**(0.5d0) * Sys%Dx
    end do

    !**** On normalise la fonction de distribution
    !**** integ(u½ f(u) du) = 1
    do i=1,nx
       f1(i)=  f1(i) / part
       D(i)= f1(i)
    end do
    II(0) = 0.d0 ; JJ(0) = 0.d0
    l = 0

123 small = 5.d-12
    l = l + 1 
    y00 = 0.d0 ; xx  = 0.d0
    yy1 = 0.d0 ; yy2 = 0.d0

    do i=1,nx
       y00= y00 + f1(i)
    end do

    do i=1,nx
       xx  = xx  + f1(i)*usq(i) 
       yy1 = yy1 + f1(i)*utt(i) 
       yy2 = yy2 + f1(i) 
       truc = (i * sys%Dx)**1.5d0
       II(i)= xx*hdu
       JJ(i)= tt*(yy1 + truc*(y00-yy2))
    end do

    do i=1,nx
       Cn = Cnst / sqrt(U(i))
       A(i)= Cn*(-II(i-1) + JJ(i-1))
       B(i)= 1. + Cn*(II(i) - II(i-1) - JJ(i) - JJ(i-1))
       C(i)= Cn*(II(i) + JJ(i))
       f1_new(i)=f1(i)
    end do

    !----Calcul de f1-------------------------------------------------------------
    call tridag(A,B,C,D,f1,nx)
    !----Critere d arret- Erreur relative-----------------------------------------

    IF (l .GT. 500) THEN
       write(*,"(2A)") tabul, "No CONVERGENCE reach in Fokker-Plank routine !"
       write(*,"(A,2(A,ES15.6))") tabul, "CONVERGENCE error = ", small, " | Calculated error = ", err
       write(*,"(2A)") tabul, "Stop Calculations !"
       Stop
    END IF

    do i=1,nx
       err = abs( (f1_new(i)-f1(i))/f1(i) )
       if(err.gt.small) go to 123 
    end do

    !**** On "denormalise" la fonction de distribution
    do i=1,nx
       f1(i)=  f1(i) * part
    end do

    DEALLOCATE (II, JJ, f1_new)
    DEALLOCATE (utt, usq, A, B, C, D)
  END subroutine FP
  !************** FIN SUBROUTINE FOKKER-PLANCK ****************************
  SUBROUTINE Heating (sys,meta, U,F)
    !INTENT
    TYPE(SysVar), INTENT(INOUT)  :: sys
    REAL(DOUBLE), DIMENSION(:), INTENT(INOUT) :: F
    REAL(DOUBLE), DIMENSION(:), INTENT(IN)    :: U
    TYPE(Species), DIMENSION(0:), INTENT(IN)  :: meta
    !LOCAL
    INTEGER      :: i, nx
    REAL(DOUBLE) :: Dx, part,alpha0, nucm,nucp
    REAL(DOUBLE) :: YY,ZZ,XX, En1, En2
    REAL(DOUBLE) :: power, Uc, Df
    REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: AC1,BC1,CC1,SECMAX
    REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: f0,nuc
    En1 = 0.d0 ; En2 = 0.d0
    nx = sys%nx ; Dx = sys%Dx
    ALLOCATE ( AC1(nx),BC1(nx),CC1(nx),SECMAX(nx) )
    ALLOCATE ( f0(nx),nuc(0:nx) )

    !**** Power calculation
    power = 0.d0
    do i = 1, nx
       nucp  = meta(0)%Nuel(i)
       Uc = qome * nucp / (nucp**2 + sys%Freq**2)
       IF (i .LT. nx) Df = F(i+1) - F(i)
       power = power - (U(i)**(1.5d0) * Uc * Df * 0.6667d0)
    END do
    sys%E = dsqrt ( sys%Powr / (power * qe) )
    !write(*,"(A,4ES15.6)") "Elec Field (V/m) & Input Power (Watt/m-3)", sys%E, Sys%Powr
    !***************************************************

    !****** PARAMETRES COLLISIONS ELASTIQUES*************
    alpha0 = gama**2
    !***************** CALCUL DES QUANTITES INITIALES**************************
    part= 0.d0
    do i=1, nx
       part = part + F(i)*dsqrt(U(i))*Dx    ! nombre de particules initiale
       En1  = En1  + F(i)*U(i)**(1.5d0)*Dx    
    end do
    !**** NORMALISATION DE LA FONCTION DE DISTRIBUTION ***********************
    do i=1,nx
       F(i) = F(i) / part
       f0(i)= F(i)
    end do

    !**** SYSTEME TRIDIAGONALE A RESOUDRE CAS GENEGRAL***************
    do i=1,nx
       XX = (U(i)-0.5d0*Dx)                               
       YY = (U(i)+0.5d0*Dx)
       ZZ = Clock%Dt/(3.d0*Dx*Dx*dsqrt(U(i)))

       IF (i == 1 ) THEN
          nucm = meta(0)%Nuel(i)
          nucp= 0.5d0*(meta(0)%Nuel(i)+meta(0)%Nuel(i+1))
       ELSE IF (i == nx) THEN
          nucp = meta(0)%Nuel(i)
          nucm= 0.5d0*(meta(0)%Nuel(i)+meta(0)%Nuel(i-1))
       ELSE
          nucm= 0.5d0*(meta(0)%Nuel(i)+meta(0)%Nuel(i-1))
          nucp= 0.5d0*(meta(0)%Nuel(i)+meta(0)%Nuel(i+1))
       END IF

       AC1(i)= -ZZ*alpha0 * sys%E**2*(XX**1.5d0)*nucm / (nucm**2 + sys%Freq**2)
       CC1(i)= -ZZ*alpha0 * sys%E**2*(YY**1.5d0)*nucp / (nucp**2 + sys%Freq**2)    
       BC1(i)= 1.d0 + ZZ*alpha0*sys%E**2*( (YY**1.5d0 *nucp / (nucp**2 + sys%Freq**2)) +&
            (XX**1.5d0*nucm / (nucm**2 + sys%Freq**2)) )
    end do
    !*****SOLUTION DU SYSTEME TRIDIAGONALE f1 AU TEMPS k+1-************************
    CALL TRIDAG (AC1,BC1,CC1,f0,F,nx)

    do i=1,nx
       F(i)= F(i) * part
       En2 = En2  + F(i)*U(i)**(1.5d0)*Dx
    end do
    !**** Diagnostic 
    diag(10)%EnProd(1) = diag(10)%EnProd(1) + abs(En2 - En1)
    !***************

    DEALLOCATE ( AC1,BC1,CC1,SECMAX )
    DEALLOCATE ( f0,nuc )
  END SUBROUTINE Heating
  !******************** Fin subroutine chauffage ***************

  !*********** SUBROUTINE Collisions Elastiques ***************!
  SUBROUTINE Elastic (sys,meta,U,F)
    !INTENT
    TYPE(SysVar), INTENT(IN)  :: sys
    REAL(DOUBLE), DIMENSION(:), INTENT(INOUT) :: F
    REAL(DOUBLE), DIMENSION(:), INTENT(IN)    :: U
    TYPE(Species), DIMENSION(0:), INTENT(IN)  :: meta
    !LOCAL
    INTEGER      :: i, nx
    REAL(DOUBLE) :: Dx, XX,YY,ZZ, En, En2
    REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: AEN,BEN,CEN
    REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: f0
    REAL(DOUBLE) :: part,alpha1,alpha2,nucm,nucp
    En = 0.d0 ; En2 = 0.d0
    nx = sys%nx ; Dx = sys%Dx
    ALLOCATE ( AEN(nx),BEN(nx),CEN(nx),f0(nx) )

    !******PARAMETRES COLLISIONS ELASTIQUES*************
    alpha1 = (2.d0 * MassR)   !(2.*me/Mhe)
    alpha2 = (meta(0)%Tp) / Dx 
    !*****************CALCUL DES QUANTITES INITIALES**************************
    part= 0.d0
    do i=1,nx
       part= part + F(i)*u(i)**(0.5d0)*Dx  ! nombre de particules initiale
       En  = En   + F(i)*u(i)**(1.5d0)*Dx
    end do
    !****NORMALISATION DE LA FONCTION DE DISTRIBUTION ***********************
    do i=1,nx
       F(i) = F(i) / part
       f0(i)= F(i)
    end do

    !****************SYSTEME TRIDIAGONALE A RESOUDRE************************************
    do i=1,nx
       XX = (U(i)-0.5d0*Dx)                                
       YY = (U(i)+0.5d0*Dx)
       ZZ = alpha1 * Clock%Dt / ( Dx*dsqrt(U(i)) )
       IF (i == 1 ) THEN
          nucm = meta(0)%Nuel(i)
          nucp = 0.5d0*( meta(0)%Nuel(i)+meta(0)%Nuel(i+1) )
       ELSE IF (i == nx) THEN
          nucp = meta(0)%Nuel(i)
          nucm = 0.5d0*( meta(0)%Nuel(i)+meta(0)%Nuel(i-1) )
       ELSE
          nucp = 0.5d0*( meta(0)%Nuel(i)+meta(0)%Nuel(i+1) )
          nucm = 0.5d0*( meta(0)%Nuel(i)+meta(0)%Nuel(i-1) )
       END IF

       AEN(i)= -ZZ*(XX**1.5d0)*(-0.5d0 + alpha2) * nucm
       CEN(i)= -ZZ*(YY**1.5d0)*( 0.5d0 + alpha2) * nucp
       BEN(i)= 1.d0 + ZZ * ( (YY**1.5d0)*(-0.5d0 + alpha2)*nucp &
            + (XX**1.5d0)*(0.5d0 + alpha2) * nucm)
    end do

    !*****SOLUTION DU SYSTEME TRIDIAGONALE f1 AU TEMPS k+1-***************************************
    CALL TRIDAG (AEN,BEN,CEN,f0,F,nx)

    do i=1,nx
       F(i)= F(i) * part
       En2 = En2  + F(i)*U(i)**(1.5d0)*Dx
    end do
    !**** Diagnostic 
    diag(10)%EnLoss(2) = diag(10)%EnLoss(2) + abs(En - En2)
    if (int(Clock%SumDt/Clock%Dt) .EQ. Clock%NumIter-2)THEN 
       diag(10)%Tx = (En - En2)
       !write(*,"(A,ES15.6)") "Powr Elastic: ", diag(10)%Tx * qe / clock%Dt
    END if
    !***************
    DEALLOCATE ( AEN,BEN,CEN,f0)
  END SUBROUTINE Elastic
END MODULE MOD_CHAUF
