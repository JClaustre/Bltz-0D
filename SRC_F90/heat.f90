!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Charlotte Boukandou (modified by Jonathan Claustre)
! Date  : 14/07/2015
! Objctv: Elastic collisions, Fokker-Planck and Heating 
!         processes in Helium
! note  : Fk-Planck equation --> see paper Vidal & Boukandou
!         doi:10.1016/j.cpc.2017.07.004
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_CHAUF
  USE F90_KIND
  USE MOD_PARAM
  IMPLICIT NONE

  INTEGER, PRIVATE :: Pwk = 0

CONTAINS

  SUBROUTINE POWER_CONTROL (Clock, sys, meta, U, F, Post_D, Cgen)
    !INTENT
    TYPE(Time)   , INTENT(IN)     :: Clock
    TYPE(SysVar) , INTENT(INOUT)  :: sys
    REAL(DOUBLE) , INTENT(IN)     :: Post_D, Cgen
    TYPE(Species), DIMENSION(0:), INTENT(IN)  :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN)  :: U, F
    !LOCAL
    INTEGER      :: i, nx
    REAL(DOUBLE) :: Dx, Fn, nuc, frq
    REAL(DOUBLE) :: power, Uc, Df, GenPwr
    REAL(DOUBLE) :: Delta, Op
    nx = sys%nx ; Dx = sys%Dx ; power = 0.d0
    GenPwr = 1.d-6 ! Time constant to start/end the generator.
    
    IF (Clock%SumDt .LT. Post_D) THEN
       !**** Increase Power
       IF (sys%P0 == 0) THEN
          sys%Powr = sys%IPowr * (1.d0 - exp( -real(Clock%SumDt) / GenPwr) )
          sys%Pcent = sys%Powr

          !**** Power calculation
          power = 0.d0
          do i = 1, nx-1
             nuc  = meta(0)%Ni*meta(0)%SecMtm(i)*gama*dsqrt(U(i))
             Uc = qome * nuc / (nuc**2 + sys%Freq**2)
             IF (i .LT. nx-1) THEN
                Df = F(i+1) - F(i)
             ELSE IF (i.EQ.nx-1) THEN
                !**** linear extrapolation for f(nx)
                Fn = F(nx-2) + (F(nx-1)-F(nx-2))/Dx
                Df = Fn - F(i)
             END IF
             power = power - (U(i)**(1.5d0) * Uc * Df * 0.6667d0)
          END do
          !**** New External Electric Field Calculation 
          sys%E = dsqrt ( sys%Powr / (power * qe) )

       ELSE
          sys%E = sys%Emax * (1.d0 - exp( -real(Clock%SumDt) / GenPwr) )
          sys%Pcent = sys%E
       END IF
       
       !***************************************************
    ELSE
       !**** Decrease External Electric source ***
       IF (Pwk == 0) sys%IPowr = sys%E
       sys%E = sys%IPowr * exp( -real(Pwk*Clock%Dt) / (GenPwr*Cgen))
       IF (sys%E.LT.1d-08) sys%E = 0.d0
       sys%Pcent = sys%E
       Pwk = Pwk+1
    END IF

    !**** Effet de peau 
    IF (sys%P0 == 1) THEN
       Op = 5.64d4 * sqrt(elec%Ni*1d-6)
       Delta = Vcel / Op * sqrt(2.d0*1d9*meta(0)%Prs/sys%freq)

       sys%E = sqrt( sys%Emax**2 * exp(-2d0*sys%Ra/Delta) )
    END IF
    
    !**** RF - electric field *********************
    IF (sys%rf == 1) THEN
       IF (sys%Freq.NE.0.d0) THEN
          sys%E = (sys%E*sqrt(2.d0)) * sin(sys%Freq * Clock%SumDt)
       END IF
       frq = 0.d0
    ELSE
       frq = sys%Freq
    END IF
    !*********************************************
    
    !**** Diagnostic to calculate the absorbed power by the plasma ***
    sys%Emoy = sys%Emoy + dabs(sys%E)
    !**** Fix the power here function of Elec field ************************!
    DO i = 1, sys%nx - 1                                                    !
       nuc  = meta(0)%Ni*meta(0)%SecMtm(i)*gama*dsqrt(U(i))
       Uc = qome * sys%E**2 / (nuc**2 + Frq**2)                             !
       IF (i .LT. nx-1) THEN
          Df = F(i+1) - F(i)
       ELSE IF (i.EQ.nx-1) THEN
          !**** linear extrapolation for f(nx) ***
          Fn = F(nx-2) + (F(nx-1)-F(nx-2))/Dx
          Df = Fn - F(i)
       END IF
       power = power - Uc * U(i)**(1.5d0) * Df * nuc * 0.6667d0             !
    END DO                                                                  !
    power = power * qe                                                      !
    sys%Pwmoy = power
    !***********************************************************************!
  END SUBROUTINE POWER_CONTROL

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Heating (sys,meta, U,F)
    !INTENT
    TYPE(SysVar), INTENT(IN)  :: sys
    REAL(DOUBLE), DIMENSION(:), INTENT(INOUT) :: F
    REAL(DOUBLE), DIMENSION(:), INTENT(IN)    :: U
    TYPE(Species), DIMENSION(0:), INTENT(IN)  :: meta
    !LOCAL
    INTEGER      :: i, nx
    REAL(DOUBLE) :: Dx, part,partf,alpha0, nucm,nucp
    REAL(DOUBLE) :: YY,ZZ,XX, En1, En2, frq, Frac_mol
    REAL(DOUBLE), DIMENSION(sys%nx) :: f0,AC1,BC1,CC1
    En1 = 0.d0 ; En2 = 0.d0
    nx = sys%nx ; Dx = sys%Dx

    !**** Calculation of the Total momentum cross secrtion.
    meta(2)%SecMtm = meta(1)%SecMtm
    Do i = 1, NumMeta
       Frac_mol = ( meta(i)%Ni / meta(0)%Ni )
       meta(2)%SecMtm(:) = meta(2)%SecMtm(:) + meta(0)%SecExc(i,:) * Frac_mol
    END Do

    !****** PARAMETRES COLLISIONS ELASTIQUES*************
    alpha0 = gama*gama
    !***************** CALCUL DES QUANTITES INITIALES**************************
    part= 0.d0 ; partf = 0.d0
    do i=1, nx
       part = part + F(i)*dsqrt(U(i))*Dx    ! nombre de particules initiale
       En1  = En1  + F(i)*U(i)**(1.5d0)*Dx    
    end do
    !**** NORMALISATION DE LA FONCTION DE DISTRIBUTION ***********************
    do i=1,nx
       F(i) = F(i) / part
       f0(i)= F(i)
    end do

    !**** Rf mode ******************************
    IF (sys%rf == 1) THEN
       frq = 0.d0
    ELSE
       frq = sys%Freq
    END IF
    
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

       AC1(i)= -ZZ*alpha0 * sys%E**2*(XX**1.5d0)*nucm / (nucm**2 + frq**2)
       CC1(i)= -ZZ*alpha0 * sys%E**2*(YY**1.5d0)*nucp / (nucp**2 + frq**2)    
       BC1(i)= 1.d0 + ZZ*alpha0*sys%E**2*( (YY**1.5d0 *nucp / (nucp**2 + frq**2)) +&
            (XX**1.5d0*nucm / (nucm**2 + frq**2)) )
    end do

    !*****SOLUTION DU SYSTEME TRIDIAGONALE f1 AU TEMPS k+1-************************
    CALL TRIDAG (AC1,BC1,CC1,f0,F,nx)
    F(nx) = 0.d0
    do i=1,nx
       F(i)= F(i) * part
       partf = partf + F(i)*dsqrt(U(i))*Dx
       En2 = En2  + F(i)*U(i)**(1.5d0)*Dx
    end do

    !**** Diagnostic 
    diag(10)%EnProd = diag(10)%EnProd + dabs(En2 - En1)
    ! *** If Emax is not enough or other... make sure Ne=[Na+Nm] *  
    ion(1)%Ni = partf * ion(1)%Ni / part
    ion(2)%Ni = partf - ion(1)%Ni
    !***************
  END SUBROUTINE Heating
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!

  !*********** SUBROUTINE Collisions Elastiques ***************!
  SUBROUTINE Elastic (sys,meta,U,F)
    !INTENT
    TYPE(SysVar), INTENT(IN)  :: sys
    REAL(DOUBLE), DIMENSION(:), INTENT(INOUT) :: F
    REAL(DOUBLE), DIMENSION(:), INTENT(IN)    :: U
    TYPE(Species), DIMENSION(0:), INTENT(INOUT)  :: meta
    !LOCAL
    INTEGER      :: i, nx
    REAL(DOUBLE) :: Dx, XX,YY,ZZ, En, En2, Ren
    REAL(DOUBLE), DIMENSION(sys%nx) :: AEN,BEN,CEN,f0
    REAL(DOUBLE) :: part,alpha1,alpha2,nucm,nucp
    En = 0.d0 ; En2 = 0.d0 ; Ren = 0.d0
    nx = sys%nx ; Dx = sys%Dx

    !******PARAMETRES COLLISIONS ELASTIQUES*************
    alpha1 = (2.d0 * MassR)   !(2.*me/Mhe)
    alpha2 = (meta(0)%Tp) / Dx
    !*****************CALCUL DES QUANTITES INITIALES**************************
    part= 0.d0
    do i=1,nx
       part= part + F(i)*u(i)**(0.5d0)*Dx  ! nombre de particules initiale
       En  = En   + F(i)*u(i)**(1.5d0)*Dx
       meta(0)%Nuel(i) = meta(0)%Ni*meta(0)%SecMtm(i)*gama*dsqrt(U(i))
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
       !**** Collision Rate e-n :
       Ren = Ren + gama * F(i) * meta(0)%SecMtm(i) * U(i) * Dx * part
    end do

    !*****SOLUTION DU SYSTEME TRIDIAGONALE f1 AU TEMPS k+1-***************************************
    CALL TRIDAG (AEN,BEN,CEN,f0,F,nx)

    do i=1,nx
       F(i)= F(i) * part
       En2 = En2  + F(i)*U(i)**(1.5d0)*Dx
    end do
    !**** Diagnostic 
    diag(11)%EnLoss = diag(11)%EnLoss + dabs(En - En2)
    !**** rate electron-neutral
    diag(18)%Tx(1) = Ren * meta(0)%Ni
    !***************
  END SUBROUTINE Elastic
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!

  !*********** SUBROUTINE FOKKER-PLANCK ***************!
  subroutine FP (sys,elec, f1, U)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Species), INTENT(IN) :: elec
    REAL(DOUBLE), DIMENSION(:), INTENT(INOUT) :: f1
    REAL(DOUBLE), DIMENSION(:), INTENT(IN)    :: U
    !LOCAL
    INTEGER :: i, l, nx
    REAL(DOUBLE) :: part, En, En2, DiagE
    REAL(DOUBLE) :: hdu, tt, nu, v3, cnst, cn, err
    REAL(DOUBLE) :: y00, yy1, yy2, xx, truc, small
    REAL(DOUBLE), DIMENSION(0:sys%nx) :: II,JJ
    REAL(DOUBLE), DIMENSION(sys%nx) :: utt, usq, A, B, C, D, f1_new
    En = 0.d0 ; En2 = 0.d0
    nx = sys%nx
    hdu = 0.5d0*sys%Dx
    tt  = 2.d0/3.d0
    v3 = ( (2.d0*qome))**(1.5d0)
    
    !********Initialisation:******************************
    part = 0.d0
    do i=1,nx
       utt(i)= 0.5d0*((U(i)+hdu)**(1.5d0) + (U(i)-hdu)**1.5d0)  
       usq(i)= tt * ((U(i)+hdu)**(1.5d0) - (U(i)-hdu)**1.5d0) /sys%Dx
       !*****Calcul des quantites initiales
       part = part + f1(i) * U(i)**(0.5d0) * Sys%Dx
       En   = En + f1(i) * U(i)**(1.5d0) * Sys%Dx
    end do
    !**** Log Coulomb (cf. NRL formulary) (Ne--> cm-3 !)
    IF (elec%Tp.LE.10.d0) THEN
       LnC = 23.d0 - log((part*1d-6)**0.5d0 / elec%Tp**1.5d0)
    ELSE IF (elec%Tp.GT.10.d0) THEN
       LnC = 24.d0 - log((part*1d-6)**0.5d0 / elec%Tp )
    END IF

    nu = (4.d0*pi*elec%Ni*LnC*qe**4) / ( (4.d0*pi*eps*me)**2 * v3)
    cnst = -(2.0d0*Clock%Dt * nu) / Sys%Dx
    
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
       write(*,"(2A,ES10.2,A,ES10.2)") tabul, "Stop Calculations ! time: ", Clock%SumDt, " Dt= ", Clock%Dt
       Stop
    END IF

    do i=1,nx
       err = dabs( (f1_new(i)-f1(i))/f1(i) )
       if(err.gt.small) go to 123 
    end do

    !**** On "denormalise" la fonction de distribution
    do i=1,nx
       f1(i) = f1(i) * part
       En2   = En2 + f1(i) * U(i)**(1.5d0) * Sys%Dx
    end do
    DiagE = dabs(1- En2/En)
    IF (DiagE.GT.1.d-08) print*, "Energy In Fokker_planck Routine not well conserved!", DiagE
  END subroutine FP
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
END MODULE MOD_CHAUF
