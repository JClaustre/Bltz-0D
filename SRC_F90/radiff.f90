!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 14/07/2015
! Objctv: Radiative and ambipolar diffusion processes in He
! note  : Norm : Calculate ambipolar diffusion using fixed 
!                mobility
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_RADIFF
  USE F90_KIND  
  USE MOD_PARAM
  IMPLICIT NONE

CONTAINS

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Radiat (sys, meta, Fosc, diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    Type(Diagnos), DIMENSION(:)    , INTENT(INOUT) :: diag
    TYPE(Species), DIMENSION(0:)   , INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(0:,0:), INTENT(IN)    :: Fosc
    !LOCAL
    INTEGER :: i, j
    REAL(DOUBLE) :: Eij, damp, EscapF, emitF
    REAL(DOUBLE) :: Kor, Gcol, Gdop, Gcd

    DO i = 1, NumMeta!18
       DO j = 0, i-1 !10
          IF (meta(i)%Aij(j).NE.0.d0) THEN
             !write(*,"(2A,ES15.6)") meta(i)%Name, meta(j)%Name, meta(i)%Aij(j)
             Eij = meta(i)%En-meta(j)%En
             IF (Eij .GT. 0.d0) THEN
                Kor = 2.876d-10 * 1d-4 * Fosc(j,i) * meta(j)%Ni * sys%Ra /&
                     (dsqrt(meta(0)%Tp*qok)*Eij)
                damp = 0.d0
                IF (Kor .GT. 1.d0) THEN
                   damp = (1.d0 + 3.221d-14* meta(j)%Ni *(1d-6) * meta(i)%Deg / (meta(j)%Deg * Eij**3) )&
                        * (6.6379d-2*Fosc(j,i)*Eij*meta(j)%Deg / (meta(i)%Deg*dsqrt(meta(0)%Tp*qok)) )
                   Gdop = 1.6d0 / ( Kor * dsqrt(pi*log(Kor)) )
                   Gcol = (2.0d0 / pi) * dsqrt( dsqrt(pi)*damp / Kor )
                   Gcd  = 2.0d0 * damp / ( pi * dsqrt(log(Kor)) )
                   EscapF = 0.d0
                   IF (Gcd/Gcol .GT. 8.d0) THEN
                      EscapF = Gcol * erf(Gcd/Gcol)
                      !print*, '>10', i,j,meta(i)%name, meta(j)%Name, meta(i)%Aij(j), Kor, EscapF
                   ELSE
                      EscapF = Gdop / exp(Gcd**2/Gcol**2) + Gcol * erf(Gcd/Gcol)
                      !print*, '<10', i,j,meta(i)%name, meta(j)%Name, meta(i)%Aij(j), Kor, EscapF
                   END IF
                   emitF = meta(i)%Aij(j) * EscapF

                   meta(j)%Updens = meta(j)%Updens + Clock%Dt* emitF * meta(i)%Ni
                   meta(i)%Updens = meta(i)%Updens - Clock%Dt* emitF * meta(i)%Ni
                   !**** Diagnostic
                   diag(3)%EnLoss = diag(3)%EnLoss + Clock%Dt* emitF * meta(i)%Ni * Eij
                   !****************
                END IF
             END IF
          END IF

       END DO
    END DO

  END SUBROUTINE Radiat
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Diffuz(sys,meta,ion,elec,F,U,diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN)    :: sys
    TYPE(Species), INTENT(INOUT) :: elec
    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    TYPE(Species), DIMENSION(:) , INTENT(INOUT) :: ion
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: F
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN) :: U
    !LOCAL
    INTEGER :: i
    REAL(DOUBLE) :: mua, mum, mue, Da, Dm, De, En, En2
    REAL(DOUBLE) :: Coef, Coef2, Sa, Sm, Se, Smeta1, Smeta2, smeta3
    En = 0.d0 ; En2 = 0.d0
    !**** Normalization of EEDF
    !**** Σ f(u)u½ du = 1
    elec%Ni = 0.d0
    DO i = 1, sys%nx
       elec%Ni = elec%Ni + F(i) * sqrt(U(i)) * sys%Dx
       En = En + F(i) * U(i)**(1.5d0) * sys%Dx
    END DO
    F(:) = F(:) / elec%Ni
    !**** Mobility and diffusion Coef calculation 
    !**** Electron
    !**** μ = -(1/3Ng).γ. Σ (df(u)u du/(du.σ))
    !**** D =  (1/3Ng).γ. Σ (f(u)u du/σ)
    Coef = gama * sys%Dx / (3.d0*meta(0)%Ni)
    DO i = 1, sys%nx
       coef2 = 0.d0
       if (meta(0)%SecMtm(i) .ne. 0.d0) Coef2 = Coef * U(i) / meta(0)%SecMtm(i)
       De = De + Coef2*F(i)
       if (i < sys%nx-1) mue = mue - Coef2 * (F(i+1)-F(i)) / sys%Dx
    END DO
    !**** Atomic & Molecular Ion
    mua = 2.68d19*1d2  / (2.96d-3 * dsqrt(meta(0)%Tp*qok) + 3.11d-2) / meta(0)%Ni ! cf. Santos
    !mua = 1.d0 / (2.96d-3 * dsqrt(meta(0)%Tp*qok) + 3.11d-2) ! cf. Belmonte
    Da  = mua * meta(0)%Tp
    mum = 2.68d19*1d2 / meta(0)%Ni
    Dm  = mum * meta(0)%Tp
    elec%Dfree = De ; ion(1)%Dfree = Da ; ion(2)%Dfree = Dm 
    elec%mobl = mue ; ion(1)%mobl = mua ; ion(2)%mobl = mum 
    !**** Ambipolar Diffusion Coefficient (cm² s-¹)
    Coef = (ion(1)%Ni*mua + ion(2)%Ni*mum + elec%Ni*mue)
    !**************************************
    elec%Damb = ion(1)%Ni * (mua*De + mue*Da) + ion(2)%Ni*(mum*De + mue*Dm)
    elec%Damb = elec%Damb / Coef
    !**************************************
    ion(1)%Damb = ion(2)%Ni * (mum*Da - mua*Dm) + elec%Ni*(mua*De + mue*Da)
    ion(1)%Damb = ion(1)%Damb / Coef
    !**************************************
    ion(2)%Damb = ion(1)%Ni * (-mum*Da + mua*Dm) + elec%Ni*(mum*De + mue*Dm)
    ion(2)%Damb = ion(2)%Damb / Coef
    !**** Diffusion Coefficient for He(2S3) | He(2S1) | He2*
    meta(1)%Damb = 8.922d-02 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
    meta(2)%Damb = meta(1)%Damb
    !**************************************
    !**** Sj = Dj.nj / Λ²  (m-3 s-¹)
    Coef2 = (sys%Ra/2.405d0)**2
    Sa = ion(1)%Damb * ion(1)%Ni / Coef2
    Sm = ion(2)%Damb * ion(2)%Ni / Coef2
    Se = elec%Damb * elec%Ni   / Coef2
    Smeta1 = meta(1)%Damb * meta(1)%Ni / Coef2
    Smeta2 = meta(2)%Damb * meta(2)%Ni / Coef2
    !**** particle balance
    elec%Ni    = elec%Ni  - Clock%Dt * Se
    ion(1)%Updens  = ion(1)%Updens  - Clock%Dt * Sa
    ion(2)%Updens  = ion(2)%Updens  - Clock%Dt * Sm
    meta(1)%Updens = meta(1)%Updens - Clock%Dt * Smeta1
    meta(2)%Updens = meta(2)%Updens - Clock%Dt * Smeta2
    SELECT CASE (NumIon) 
    CASE (3)   
       ion(NumIon)%Damb  = 7.102d-02 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
       Smeta3 = ion(NumIon)%Damb  * ion(NumIon)%Ni  / Coef2
       ion(NumIon)%Updens  = ion(NumIon)%Updens  - Clock%Dt * Smeta3
    END SELECT

    DO i = 1, sys%Nx
       !**** Renormalization with modified 
       !**** electron density
       F(i) = F(i) * elec%Ni
       En2 = En2 + F(i)*U(i)**(1.5d0)*sys%Dx
    END DO
    !**** Diagnostic
    diag(9)%EnLoss = diag(9)%EnLoss + (En-En2)
    diag(9)%Tx = diag(9)%Tx + Clock%Dt * Se  
  END SUBROUTINE Diffuz
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Diffuz_Norm (sys,meta,ion,elec,F,U,diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN)    :: sys
    TYPE(Species), INTENT(INOUT) :: elec
    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    TYPE(Species), DIMENSION(:) , INTENT(INOUT) :: ion
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: F
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN) :: U
    !LOCAL
    INTEGER :: i
    REAL(DOUBLE) :: mua, mum, Da, Dm, Damb, En, En2
    REAL(DOUBLE) :: Ng_atm, Coef, Sa, Sm, Se, Smeta1, Smeta2, smeta3
    En = 0.d0 ; En2 = 0.d0

    elec%Ni = 0.d0
    DO i = 1, sys%nx
       elec%Ni = elec%Ni + F(i) * sqrt(U(i)) * sys%Dx
       En = En + F(i) * U(i)**(1.5d0) * sys%Dx
    END DO
    F(:) = F(:) / elec%Ni

    mua = 10.3d0*1d-4 ; mum = 21.d0*1d-4 ! mobility (m2/V/s) (p=1 atm.)
    Ng_atm = 2.45d+19 * 1d+6 ! Gaz density at 1 atm. & 300K
    !**** Ion diffusion coefficients
    mua = mua * Ng_atm / meta(0)%Ni
    mum = mum * Ng_atm / meta(0)%Ni
    Da  = mua * (elec%Tp + meta(0)%Tp)
    Dm  = mum * (elec%Tp + meta(0)%Tp)
    !**** Ambipolar Diffusion Coefficient for 2S3 excited state
    meta(1)%Damb = 8.992d-2 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
    meta(2)%Damb = meta(1)%Damb
    !**** Ambipolar diffusion
    Damb = ( Da*ion(1)%Ni + Dm*ion(2)%Ni ) / (ion(1)%Ni + ion(2)%Ni)

    !**** Sj = Dj.nj / Λ²  (m-3 s-¹)
    Coef = (sys%Ra/2.405d0)**2
    Sa = Damb * ion(1)%Ni / Coef
    Sm = Damb * ion(2)%Ni / Coef
    Se = Damb * elec%Ni   / Coef
    Smeta1 = meta(1)%Damb * meta(1)%Ni / Coef
    Smeta2 = meta(2)%Damb * meta(2)%Ni / Coef
    !**** particle balance
    elec%Ni        = elec%Ni        - Clock%Dt * Se
    ion(1)%Updens  = ion(1)%Updens  - Clock%Dt * Sa
    ion(2)%Updens  = ion(2)%Updens  - Clock%Dt * Sm
    meta(1)%Updens = meta(1)%Updens - Clock%Dt * Smeta1
    meta(2)%Updens = meta(2)%Updens - Clock%Dt * Smeta2
    SELECT CASE (NumIon) 
    CASE (3)   
       ion(NumIon)%Damb  = 7.102d-02 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
       Smeta3 = ion(NumIon)%Damb  * ion(NumIon)%Ni  / Coef
       ion(NumIon)%Updens  = ion(NumIon)%Updens  - Clock%Dt * Smeta3
    END SELECT
    !**** Diagnostic
    ion(1)%mobl = mua ; ion(2)%mobl = mum
    ion(1)%Dfree = Da ; ion(2)%Dfree = Dm

    DO i = 1, sys%Nx
       !**** Renormalization with modified 
       !**** electron density
       F(i) = F(i) * elec%Ni
       En2 = En2 + F(i)*U(i)**(1.5d0)*sys%Dx
    END DO
    !**** Diagnostic
    diag(9)%EnLoss = diag(9)%EnLoss + (En-En2)
    diag(9)%Tx = diag(9)%Tx + Clock%Dt * Se
  END SUBROUTINE Diffuz_Norm

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Diffuz_C(sys,meta,ion,elec,F,U,diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN)    :: sys
    TYPE(Species), INTENT(INOUT) :: elec
    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    TYPE(Species), DIMENSION(:) , INTENT(INOUT) :: ion
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: F
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN) :: U
    !LOCAL
    INTEGER :: i, j
    REAL(DOUBLE) :: mua, mum, Da, Dm, En, En2
    REAL(DOUBLE) :: Coef, Coef2, Sa, Sm, Se, Smeta1, Smeta2, smeta3
    REAL(DOUBLE) :: Mp, M, Mp3, Lambda, a, b, Damb
    REAL(DOUBLE), DIMENSION(sys%nx) :: Dee
    En = 0.d0 ; En2 = 0.d0
    !**** Normalization of EEDF
    !**** Σ f(u)u½ du = 1
    elec%Ni = 0.d0
    DO i = 1, sys%nx
       elec%Ni = elec%Ni + F(i) * sqrt(U(i)) * sys%Dx
       En = En + F(i) * U(i)**(1.5d0) * sys%Dx
    END DO
    !**** Atomic & Molecular Ion
    mua = 2.68d19*1d2  / (2.96d-3 * dsqrt(meta(0)%Tp*qok) + 3.11d-2) / meta(0)%Ni ! cf. Santos
    Da  = mua * meta(0)%Tp
    mum = 2.68d19*1d2 / meta(0)%Ni
    Dm  = mum * meta(0)%Tp

    !Approche2: Formule Delcroix
    !Mobilite des electrons et coefficient de diffusion libre Formule Delcroix
    do i=1,sys%Nx                                                                  
       Dee(i) = gama**2 * U(i)/ (3.d0*meta(0)%Nuel(i) )
    end do

    !*******Calculs des coefficients de diffusion ambipolaire********
    
    !Coefficient de diffusion des ions atomiques Daa Phelps
    ion(1)%Damb = Da*(1.d0 + (elec%Tp/meta(0)%Tp))         

    !Coefficient de diffusion des ions moleculaires  Dma
    ion(2)%Damb = Dm*(1.d0 + (elec%Tp/meta(0)%Tp))      

    !**************Calcul du parametre M = na Daa + nm Dma***********
    Lambda = (sys%Ra/2.405d0)**2
    M      = ( (ion(1)%Ni*ion(1)%Damb) + (ion(2)%Ni*ion(2)%Damb) ) /Lambda

    !*******Calcul du coefficient de diffuision ambipolaire moyen pour les deux especes d ions***********
    Damb=(M*Lambda)/(ion(1)%Ni+ ion(2)%Ni)

    !Coefficient de diffusion pour les metastables 2S3 
    meta(1)%Damb = 8.992d-2 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
    meta(2)%Damb = meta(1)%Damb

    !***Determination de l energie minimale du seuil de perte des
    !***electrons
    Mp  = 0.d0
    Mp3 = 0.d0

    Coef = gama * sys%Dx / (3.d0*meta(0)%Ni)
    do i=sys%nx-1,1,-1
       !***Etape1: Recherche de l intervalle contenant le  point fixe 
       Mp3 = Mp ; coef2 = 0.d0
       if (meta(0)%SecMtm(i) .ne. 0.d0) Coef2 = Coef * U(i) / meta(0)%SecMtm(i)
       Mp= Mp + coef2 * F(i) / Lambda
       if (Mp.gt.M) then
          j=i
          !Interpolation lineaire entre les points M1 et M2                   
          a = (Mp3-Mp)/(2.d0*sys%Dx)    
          b = Mp - U(j)*a
          Vg= (M-b)/a
          exit
       end if
    end do

    !write(*,"(ES15.4)") Vg
    
    !****** *************************Calcul du nombre de particules et de l energie apres diffusion****
    elec%Ni =0.d0
    En2=0.d0
    DO i = 1, sys%nx
       !**** Distribution finale: Approche implicite
       if(i.GE.j)then
          F(i)= F(i)/(1.d0+ (clock%Dt*Dee(i))/Lambda)    
       end if
       elec%Ni = elec%Ni + F(i) * sqrt(U(i)) * sys%Dx    !Densite electronique finale apres diffusion
       En2 = En2 + F(i) * U(i)**(1.5d0) * sys%Dx         !Energie finale apres diffusion 
    END DO

    !**** Sj = Dj.nj / Λ²  (m-3 s-¹)
    Sa = Damb * ion(1)%Ni / Lambda
    Sm = Damb * ion(2)%Ni / Lambda
    Se = Damb * elec%Ni   / Lambda
    Smeta1 = meta(1)%Damb * meta(1)%Ni / Lambda
    Smeta2 = meta(2)%Damb * meta(2)%Ni / Lambda

    !**** particle balance
    elec%Ni        = elec%Ni        - Clock%Dt * Se
    ion(1)%Updens  = ion(1)%Updens  - Clock%Dt * Sa
    ion(2)%Updens  = ion(2)%Updens  - Clock%Dt * Sm
    meta(1)%Updens = meta(1)%Updens - Clock%Dt * Smeta1
    meta(2)%Updens = meta(2)%Updens - Clock%Dt * Smeta2
!    print*, Clock%Dt * Se, Clock%Dt * Sa + Clock%Dt * Sm
!    Stop

    !**** Diagnostic
    ion(1)%mobl = mua ; ion(2)%mobl = mum
    ion(1)%Dfree = Da ; ion(2)%Dfree = Dm
    diag(9)%Tx = diag(9)%Tx + Clock%Dt * Se
    diag(9)%EnLoss = diag(9)%EnLoss + (En-En2)
  END SUBROUTINE Diffuz_C
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
END MODULE MOD_RADIFF
