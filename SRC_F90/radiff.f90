!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 14/07/2015
! Objctv: Radiative and ambipolar diffusion processes in He
! note  : data's and analytic formula in 
!         Luis Alves et al (doi:10.1088/0022-3727/25/12/007)
!         M Santos et al (doi:10.1088/0022-3727/47/26/265201)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_RADIFF
  USE F90_KIND  
  USE MOD_PARAM
  IMPLICIT NONE

CONTAINS

  !***********************************************************************
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

    DO i = 1, NumMeta
       DO j = 0, 3
          IF ((meta(i)%Nl == 1 .and. meta(i)%Ns==1 .and. j==0) .or. &
               (i==8 .and.j==3) .or. ((i==7 .or. i==3) .and. j==1) .or. (i==4 .and. j==2)) THEN
             Eij = meta(i)%En-meta(j)%En
             IF (Eij .GT. 0.d0) THEN
                Kor = 2.876d-10 * 1d2 * Fosc(j,i) * meta(j)%Ni * sys%Ra /&
                     (dsqrt(meta(0)%Tp)*Eij) 
                damp = 0.d0
                IF (Kor .GT. 1.d0) THEN
                   damp = 1.d0 + 3.221d-14* meta(j)%Ni * meta(i)%Deg / (meta(j)%Deg * Eij**3)
                   damp = damp * (6.6379d-2*Fosc(j,i)*Eij*meta(j)%Deg / (meta(i)%Deg*dsqrt(meta(0)%Tp)) )
                   Gdop = 1.6d0 / ( Kor * dsqrt(pi*log(Kor)) )
                   Gcol = (2.0d0 / pi) * dsqrt( dsqrt(pi)*damp / Kor )
                   Gcd  = 2.0d0 * damp / ( pi * dsqrt(log(Kor)) )
                   EscapF = 0.d0
                   IF (Gcd/Gcol .GT. 10.d0) THEN
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
                   diag(3)%DnLoss(i) = diag(3)%DnLoss(i) + Clock%Dt* emitF * meta(i)%Ni
                   IF(j>0) diag(3)%DnProd(j) = diag(3)%DnProd(j) + Clock%Dt* meta(i)%Ni * emitF
                   diag(3)%EnLoss(i) = diag(3)%EnLoss(i) + Clock%Dt* emitF * meta(i)%Ni * Eij
                   !****************
                END IF
             END IF
          END IF
       END DO
    END DO
  END SUBROUTINE Radiat
  !***********************************************************************
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
       if (meta(0)%SecMtm(i) .ne. 0.d0) Coef2 = Coef * U(i) / meta(0)%SecMtm(i)
       De = De + Coef2*F(i)
       if (i < sys%nx-1) mue = mue - Coef2 *(F(i+1)-F(i))
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
    IF (NumIon == 3) THEN
       ion(3)%Damb  = 7.102d-02 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
       Smeta3 = ion(3)%Damb  * ion(3)%Ni  / Coef2
       ion(3)%Updens  = ion(3)%Updens  - Clock%Dt * Smeta3
    END IF
    !**** Diagnostic
    diag(9)%DnLoss(1) = diag(9)%DnLoss(1) + Clock%Dt * Smeta1
    diag(9)%DnLoss(2) = diag(9)%DnLoss(2) + Clock%Dt * Smeta2
    diag(9)%DnLoss(NumMeta+1) = diag(9)%DnLoss(NumMeta+1) + Clock%Dt * Sa
    diag(9)%DnLoss(NumMeta+2) = diag(9)%DnLoss(NumMeta+2) + Clock%Dt * Sm
    diag(9)%DnLoss(3) = diag(9)%DnLoss(3) + Clock%Dt * Se
    diag(9)%Tx = Se
    DO i = 1, sys%Nx
       !**** Renormalization with modified 
       !**** electron density
       F(i) = F(i) * elec%Ni
       En2 = En2 + F(i)*U(i)**(1.5d0)*sys%Dx
    END DO
    !**** Diagnostic
    diag(9)%EnLoss(2) = diag(9)%EnLoss(2) + (En-En2)
  END SUBROUTINE Diffuz
  !***********************************************************************
  SUBROUTINE Alves_Diffuz (sys,meta,ion,elec,F,U,diag)
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
    IF (NumIon == 3) THEN
       ion(3)%Damb  = 7.102d-02 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
       Smeta3 = ion(3)%Damb  * ion(3)%Ni  / Coef
       ion(3)%Updens  = ion(3)%Updens  - Clock%Dt * Smeta3
    END IF
    !**** Diagnostic
    ion(1)%mobl = mua ; ion(2)%mobl = mum
    ion(1)%Dfree = Da ; ion(2)%Dfree = Dm
    diag(9)%DnLoss(1) = diag(9)%DnLoss(1) + Clock%Dt * Smeta1
    diag(9)%DnLoss(2) = diag(9)%DnLoss(2) + Clock%Dt * Smeta2
    diag(9)%DnLoss(NumMeta+1) = diag(9)%DnLoss(NumMeta+1) + Clock%Dt * Sa
    diag(9)%DnLoss(NumMeta+2) = diag(9)%DnLoss(NumMeta+2) + Clock%Dt * Sm
    diag(9)%DnLoss(3) = diag(9)%DnLoss(3) + Clock%Dt * Se
    diag(9)%Tx = Se

    DO i = 1, sys%Nx
       !**** Renormalization with modified 
       !**** electron density
       F(i) = F(i) * elec%Ni
       En2 = En2 + F(i)*U(i)**(1.5d0)*sys%Dx
    END DO
    !**** Diagnostic
    diag(9)%EnLoss(2) = diag(9)%EnLoss(2) + (En-En2)
  END SUBROUTINE Alves_Diffuz
  !***********************************************************************
END MODULE MOD_RADIFF
