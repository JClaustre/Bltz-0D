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
  !**** Radiative transfer ***
  SUBROUTINE Radiat (sys, meta, Fosc, diag)
    !**** INTENT ***
    TYPE(SysVar) , INTENT(IN) :: sys
    Type(Diagnos), DIMENSION(:)    , INTENT(INOUT) :: diag
    TYPE(Species), DIMENSION(0:)   , INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(0:,0:), INTENT(IN)    :: Fosc
    !**** LOCAL ***
    INTEGER :: i, j
    REAL(DOUBLE) :: Eij, damp, EscapF, emitF
    REAL(DOUBLE) :: Kor, Gcol, Gdop, Gcd, Rate
    Rate=0.d0 ; diag(3)%Tx(:)=0.d0
    diag(3)%InM2 =0.d0 ; diag(3)%OutM2 =0.d0 ; diag(3)%InM1 =0.d0
    
    DO i = 3, NumMeta
       DO j = 0, i-1
          IF (meta(i)%Aij(j).NE.0.d0) THEN
             Eij = meta(i)%En-meta(j)%En
             Kor = 2.8764d-10 * 1d-4 * Fosc(j,i) * meta(j)%Ni * sys%Ra /&
                  (dsqrt(meta(0)%Tp*qok)*Eij)
             IF (Kor.GT.1.d0 .and. i.NE.3) THEN ! Add i.NE.3 because of Kor > 1 during some cases!
                damp = (1.d0 + 3.221d-14* meta(j)%Ni *(1d-6) * meta(i)%Deg / (meta(j)%Deg * Eij**3) )&
                     * (6.6379d-2*Fosc(j,i)*Eij*meta(j)%Deg / (meta(i)%Deg*dsqrt(meta(0)%Tp*qok)) )
                Gdop = 1.6d0 / ( Kor * dsqrt(pi*log(Kor)) )
                Gcol = (2.0d0 / pi) * dsqrt( dsqrt(pi)*damp / Kor )
                Gcd  = 2.0d0 * damp / ( pi * dsqrt(log(Kor)) )
                IF (Gcd/Gcol .GT. 8.d0) THEN
                   EscapF = Gcol * erf(Gcd/Gcol)
                   !print*, '>10', i,j,meta(i)%name, meta(j)%Name, meta(i)%Aij(j), Kor, EscapF
                ELSE
                   EscapF = Gdop / exp(Gcd**2/Gcol**2) + Gcol * erf(Gcd/Gcol)
                   !print*, '<10', i,j,meta(i)%name, meta(j)%Name, meta(i)%Aij(j), Kor, EscapF
                END IF
             ELSE
                EscapF = 1.d0
             END IF
             emitF = meta(i)%Aij(j) * EscapF
             IF (i.NE.3.or.j.NE.1) THEN
                meta(i)%Updens = meta(i)%Updens - Clock%Dt* emitF * meta(i)%Ni
                meta(j)%Updens = meta(j)%Updens + Clock%Dt* emitF * meta(i)%Ni
             END IF
             !**** Energy conservation Diagnostic ***
             diag(3)%EnLoss = diag(3)%EnLoss + Clock%Dt* emitF * meta(i)%Ni * Eij
             !**** Rate calcul for adaptative time-step ***
             IF (emitF .GT. MaxR) MaxR = emitF
             
             !***************** Diagnostic for relative importance of reactions (m-3/s)
             IF ((emitF*meta(i)%Ni).GT.Rate) THEN
                Rate = emitF * meta(i)%Ni
                diag(3)%Tx(2) = real(i) ; diag(3)%Tx(3) = real(j)
             END IF
             diag(3)%Tx(1) = diag(3)%Tx(1) + emitF * meta(i)%Ni
             !****************

             !*************** Diagnostic for metastable and 2^3P rates (s-1)
             IF (j.EQ.3) THEN !**** 2P3 <-- N0
                diag(3)%InM2 = diag(3)%InM2 + emitF* meta(i)%Ni
             END IF
             IF (i.EQ.3 .and. j==1) THEN !**** 2P3 --> 2S3 
                diag(15)%OutM2 = diag(15)%OutM2 + emitF
                diag(15)%InM1  = diag(15)%InM1 + emitF
             ELSE IF (i.EQ.3.and.j.NE.1) THEN !**** 2P3 --> N0
                diag(3)%OutM2 = diag(3)%OutM2 + emitF
             END IF
             !*************** Diagnostic for metastable and 2^3P rates (s-1)
             IF (j.EQ.1.and.i.NE.3) THEN !**** No --> 2S3
                diag(3)%InM1 = diag(3)%InM1 + emitF* meta(i)%Ni
             END IF
             !****************
          END IF

       END DO
    END DO

  END SUBROUTINE Radiat

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !**** Santos diffusion ***
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
    DO i = 1, sys%nx-1
       coef2 = 0.d0
       if (meta(0)%SecMtm(i) .ne. 0.d0) Coef2 = Coef * U(i) / meta(0)%SecMtm(i)
       De = De + Coef2*F(i)
       IF (i < sys%nx-1) THEN
          mue = mue - Coef2 * (F(i+1)-F(i)) / sys%Dx
       ELSE IF (i == sys%nx-1) THEN
          !**** linear extrapolation for f(nx) - f(nx-1)
          F(sys%nx) = F(sys%nx-2) + (F(sys%nx-1)-F(sys%nx-2))/(sys%Dx)
          mue = mue - Coef2 * (F(i+1)-F(i)) / sys%Dx
       END IF
    END DO
    F(sys%nx) = 0.d0

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
    elec%Damb = ( ion(1)%Ni * (mua*De + mue*Da) + ion(2)%Ni*(mum*De + mue*Dm) ) / Coef
    !**************************************
    ion(1)%Damb = ( ion(2)%Ni * (mum*Da - mua*Dm) + elec%Ni*(mua*De + mue*Da) ) / Coef
    !**************************************
    ion(2)%Damb = ( ion(1)%Ni * (mua*Dm - mum*Da) + elec%Ni*(mum*De + mue*Dm) ) / Coef
    !**** Diffusion Coefficient for He(2S3) | He(2S1) | He2*
    meta(1)%Damb = 8.922d-02 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
    meta(2)%Damb = meta(1)%Damb
    !**************************************
    !**** Sj = Dj.nj / Λ²  (m-3 s-¹)
    Coef2 = (sys%Ra/2.4048d0)**2
    Sa = ion(1)%Damb * ion(1)%Ni / Coef2
    Sm = ion(2)%Damb * ion(2)%Ni / Coef2
    Se = elec%Damb   * elec%Ni   / Coef2
    Smeta1 = meta(1)%Damb * meta(1)%Ni / Coef2
    Smeta2 = meta(2)%Damb * meta(2)%Ni / Coef2
    !**** particle balance
    elec%Ni        = elec%Ni        - Clock%Dt * Se
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
    diag(9)%SumTx = diag(9)%SumTx + Clock%Dt * Se  
    diag(9)%Tx(1) = Se
  END SUBROUTINE Diffuz

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !**** Alves Diffusion ***
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
    diag(9)%SumTx = diag(9)%SumTx + Clock%Dt * Se
  END SUBROUTINE Diffuz_Norm
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!

  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  !**** Charlotte Diffusion ***
  SUBROUTINE Diffuz_Gaine (sys,meta,ion,elec,F,U,diag)
    !INTENT
    TYPE(SysVar) , INTENT(IN)    :: sys
    TYPE(Species), INTENT(INOUT) :: elec
    Type(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    TYPE(Species), DIMENSION(:) , INTENT(INOUT) :: ion
    REAL(DOUBLE) , DIMENSION(:) , INTENT(INOUT) :: F
    TYPE(Species), DIMENSION(0:), INTENT(INOUT) :: meta
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN) :: U
    !LOCAL
    INTEGER :: i, iVg
    REAL(DOUBLE) :: mua, mum, mue, Da, Dm, Damb, De, En, En2
    REAL(DOUBLE) :: Coef, Coef2, Coef3, ratx
    REAL(DOUBLE) :: Sa, Sm, Se, Smeta1, Smeta2, smeta3
    En = 0.d0 ; En2 = 0.d0 ; mue = 0.d0

    elec%Ni = 0.d0 ; elec%Tp = 0.d0
    DO i = 1, sys%nx
       elec%Ni = elec%Ni + F(i) * sqrt(U(i)) * sys%Dx
       elec%Tp = elec%Tp + ( F(i) * dsqrt(U(i)**3) * sys%Dx )
       En = En + F(i) * U(i)**(1.5d0) * sys%Dx
    END DO
    elec%Tp = elec%Tp * 0.6667d0 / elec%Ni
    F(:) = F(:) / elec%Ni

    !**** Ion diffusion coefficients ***
    mua = 2.68d19*1d2  / (2.96d-3 * dsqrt(meta(0)%Tp*qok) + 3.11d-2) / meta(0)%Ni ! cf. Santos
    mum = 2.68d19*1d2 / meta(0)%Ni
    !**** Ion diffusion coefficients ***
    Da  = mua * (elec%Tp + meta(0)%Tp)
    Dm  = mum * (elec%Tp + meta(0)%Tp)
    !**** Ambipolar Diffusion Coefficient for 2S3 excited state ***
    meta(1)%Damb = 8.992d-2 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
    meta(2)%Damb = meta(1)%Damb
    !**** Ambipolar diffusion ***
    Damb = ( Da*ion(1)%Ni + Dm*ion(2)%Ni ) / (ion(1)%Ni + ion(2)%Ni)

    !**** Sj = Dj.nj / Λ² (m-3 s-¹) ***
    Coef = (sys%Ra/2.405d0)**2
    Sa = Damb * ion(1)%Ni / Coef
    Sm = Damb * ion(2)%Ni / Coef
    Smeta1 = meta(1)%Damb * meta(1)%Ni / Coef
    Smeta2 = meta(2)%Damb * meta(2)%Ni / Coef
    !**** particle balance ***
!    ion(1)%Updens  = ion(1)%Updens  - Clock%Dt * Sa
!    ion(2)%Updens  = ion(2)%Updens  - Clock%Dt * Sm
!    meta(1)%Updens = meta(1)%Updens - Clock%Dt * Smeta1
!    meta(2)%Updens = meta(2)%Updens - Clock%Dt * Smeta2
!    SELECT CASE (NumIon) 
!    CASE (3)   
!       ion(NumIon)%Damb  = 7.102d-02 * 1d-4 * (meta(0)%Tp*qok)**(1.5d0) / meta(0)%Prs
!       Smeta3 = ion(NumIon)%Damb  * ion(NumIon)%Ni  / Coef
!       ion(NumIon)%Updens  = ion(NumIon)%Updens  - Clock%Dt * Smeta3
!    END SELECT
    !**** Electron diffusion ***
    Se = (Sa + Sm)
    !**** Calcul de la diffusion libre moyenne electronique + mobility ***
    Coef = gama * sys%Dx / (3.d0*meta(0)%Ni)
    DO i = 1, sys%nx-1
       coef2 = 0.d0
       if (meta(0)%SecMtm(i) .ne. 0.d0) Coef2 = Coef * U(i) / meta(0)%SecMtm(i)
       De = De + Coef2*F(i)
       IF (i < sys%nx-1) THEN
          mue = mue - Coef2 * (F(i+1)-F(i)) / sys%Dx
       ELSE IF (i == sys%nx-1) THEN
          !**** linear extrapolation for f(nx) - f(nx-1)
          F(sys%nx) = F(sys%nx-2) + (F(sys%nx-1)-F(sys%nx-2))/(sys%Dx)
          mue = mue - Coef2 * (F(i+1)-F(i)) / sys%Dx
       END IF
    END DO
    F(sys%nx) = 0.d0
    elec%Dfree = De ; ion(1)%Dfree = Da ; ion(2)%Dfree = Dm 
    elec%mobl = mue ; ion(1)%mobl = mua ; ion(2)%mobl = mum 

    !**** Calcul du potentiel de gaine dans le cas de la diffusion 
    !**** ambipolaire ***
    Coef2 = 0.d0 ; Coef = (sys%Ra/2.405d0)**2
    DO i = sys%nx, 1, -1
       Coef3 = Coef2
       Coef2 = Coef2 + De * F(i) * elec%Ni * U(i)**0.5d0 * sys%Dx / Coef
       IF (Coef2.GT.Se) THEN
          iVg = i
          Vg = (Se - Coef2) * sys%Dx / (Coef3-Coef2) + sys%Dx*i
          exit
       END IF
    END DO

    DO i = 1, sys%Nx
       F(i) = F(i) * elec%Ni
!       IF (i.EQ.iVg) F(i) = F(i) - Clock%Dt * F(i) * De * ((i+1)*sys%Dx-Vg) / sys%Dx / Coef
!       IF (i.GT.iVg) F(i) = F(i) - Clock%Dt * F(i) * De / Coef
       En2 = En2 + F(i)*U(i)**(1.5d0)*sys%Dx
    END DO
    ratx = De/ Coef ! Rate of change of the EEDF (s-1)
    !**** electron density ***
    elec%Ni = 0.d0
    DO i = 1, sys%Nx
       elec%Ni = elec%Ni + F(i) * sqrt(U(i)) * sys%Dx
    END DO
    IF (ratx .GT. maxR) maxR = ratx

    !**** Energy conservation Diagnostic ***
    diag(9)%EnLoss = diag(9)%EnLoss + (En-En2)
    !***************** Diagnostic for relative importance of reactions (m-3/s)
    diag(9)%Tx(1) = Se
    !*************** Diagnostic for metastable and 2^3S rates (s-1) for MEOP
    diag(9)%OutM1 = meta(1)%Damb / Coef
    !***************
  END SUBROUTINE Diffuz_Gaine
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!

END MODULE MOD_RADIFF
