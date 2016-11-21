!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 15/11/2016
! Objctv: Pumping by laser and sublevels of 2S3 and 2P3 excited states
! note  : 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_PUMP
  USE F90_KIND  
  USE MOD_PARAM
  IMPLICIT NONE

CONTAINS

  Subroutine Sublev_coll (Clock,meta,pop,Tij,lasr)
    !INTENT
    TYPE(TIME)   , INTENT(IN) :: Clock
    TYPE(Laser)  , INTENT(IN) :: lasr
    REAL(DOUBLE) , DIMENSION(Npop2,Npop1,3) :: Tij
    TYPE(Excited), DIMENSION(2), INTENT(INOUT) :: pop
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    !LOCAL
    INTEGER :: i, j, k, switch
    REAL(DOUBLE) :: Dt, nu_ij 
    REAL(DOUBLE) :: gammak, Dop, Omeg, Is, Sec, vm
    REAL(DOUBLE) :: tot1, tot2, N1, N2
    REAL(DOUBLE), DIMENSION(2,18) :: Updens
    !**************************
    switch = 0
    N1 = real(Npop1) ; N2 = real(Npop2)
    Dt = Clock%Dt ; Updens(:,:)  = 0.d0
    pop(1)%T_relax = 3.66667d0 * 100.d0 * (meta(1)%Ni/meta(0)%Ni)
    pop(2)%T_relax = 1d-07
    !**** Laser variables ***
    ! mean velocity of the metastables
    vm   = sqrt(2.d0*kb*meta(0)%Tp*qok/mhe) 
    ! Laser frequency (1083 nm)
    Omeg = ppi*Vcel / lasr%Lwave
    ! Doppler Width due to metastable velocity distribution
    Dop  = Omeg * Vm / (ppi*Vcel)
    ! Laser coefficient 
    gammak = sqrt(pi/2.d0) * fineS * fosc(1,3) * lasr%Is / (me*Omeg*Dop*lasr%Sec)
    !************************
    IF (lasr%OnOff.EQ.1 .and. Clock%SumDt.GT.lasr%Stime) THEN
       switch = 1
    END IF

    !**** Dn due to the radiative transitions (2P3 -> 2S3) ***
    pop(2)%Dn_rad = - pop(1)%Dn_rad
    !**** Dn total after one loop in Boltz collisions ***
    pop(1)%Dn_tot = meta(1)%Ni - pop(1)%Ntot
    pop(2)%Dn_tot = meta(3)%Ni - pop(2)%Ntot
    !**** Dn due to other collisions than radiative transfer ***
    pop(1)%Dn_o   = pop(1)%Dn_tot - pop(1)%Dn_rad
    pop(2)%Dn_o   = pop(2)%Dn_tot - pop(2)%Dn_rad
    !***********************************************************
    !**** Update of sublevels densities due to Boltz collisions ***
    Updens(1,1:Npop1) = Updens(1,1:Npop1) + pop(1)%Dn_o   / real(N1)
    Updens(2,1:Npop2) = Updens(2,1:Npop2) + pop(2)%Dn_tot / real(N2)

    DO j = 1, Npop2
       Do i = 1, Npop1
          !**** Update A_i sublevels due to radiative transitions (2P3 -> 2S3) ***
          Updens(1,i) = Updens(1,i) + (pop(1)%Dn_rad/N2) * Tij(j,i,1)
          !**** Update A_i and B_j sublevels due to laser absorption ***
          IF (switch.EQ.1) THEN ! If laser activated (1)
             DO k = 1, lasr%Ntr
                IF ( int(Tij(j,i,3)).EQ.lasr%Ck(k) ) THEN ! Function of transitions Ck
                   IF ( int(Tij(j,i,2)).EQ.lasr%plz) THEN ! Function of polarization
                      nu_ij = Tij(j,i,1) * gammak
                      Updens(1,i) = Updens(1,i) + Dt * nu_ij * ( pop(2)%Ni(j) - pop(1)%Ni(i) )
                      Updens(2,j) = Updens(2,j) + Dt * nu_ij * ( pop(1)%Ni(i) - pop(2)%Ni(j) )
                   END IF
                END IF
             END DO    
          END IF

       END DO
    END DO

    DO i = 1, Npop1
       !**** Update A_i sublevels due to relaxation from collisions ***
       Updens(1,i) = Updens(1,i) + Dt * ( meta(1)%Ni/N1 - pop(1)%Ni(i) ) / pop(1)%T_relax
    END DO
    DO j = 1, Npop2
       !**** Update B_j sublevels due to relaxation from collisions ***
       Updens(2,j) = Updens(2,j) + Dt * (meta(3)%Ni/N2 - pop(2)%Ni(j) ) / pop(2)%T_relax
    END DO

    !**** Update of the pop densities and the total densities of 2S3 and 2P3 ***
    pop(1)%Ni(1:Npop1) = pop(1)%Ni(1:Npop1) + Updens(1,1:Npop1)
    pop(2)%Ni(1:Npop2) = pop(2)%Ni(1:Npop2) + Updens(2,1:Npop2)

    tot1 = 0.d0 ; tot2 = 0.d0
    do i = 1, Npop2    
       if (i.LE.Npop1) tot1 = tot1 + pop(1)%Ni(i)
       tot2 = tot2 + pop(2)%Ni(i)
    END do
    meta(1)%Ni = tot1 ; meta(3)%Ni = tot2

  END Subroutine Sublev_coll

END MODULE MOD_PUMP
