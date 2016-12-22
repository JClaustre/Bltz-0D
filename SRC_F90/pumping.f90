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

  Subroutine Sublev_coll (Clock,meta,pop,Tij,lasr,iter)
    !INTENT
    INTEGER      , INTENT(IN) :: iter
    TYPE(TIME)   , INTENT(IN) :: Clock
    TYPE(Laser)  , INTENT(IN) :: lasr
    REAL(DOUBLE) , DIMENSION(Npop2,Npop1,3) :: Tij
    TYPE(Excited), DIMENSION(2), INTENT(INOUT) :: pop
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    !LOCAL
    INTEGER :: i, j, k, switch
    REAL(DOUBLE) :: Dt, nu_ij, E_meta, E_2P3, errmax
    REAL(DOUBLE) :: gammak, Dop, Omeg, Is, Sec, vm
    REAL(DOUBLE) :: tot1, tot2, N1, N2, pol
    REAL(DOUBLE), DIMENSION(2,18) :: Updens
    !**************************
    switch = 0 ; errmax = 1.d-10
    N1 = real(Npop1) ; N2 = real(Npop2)
    Dt = Clock%Dt ; Updens(:,:)  = 0.d0
    pop(1)%T_relax = 1.d0/(meta(1)%Damb / (sys%Ra/2.405d0)**2)
    pop(2)%T_relax = 1.d0/(3.2d06*1.33322*meta(0)%Prs)
    pop(1)%tau_e   = 1d-06  ! (s)
    pop(1)%Te      = meta(0)%Ni * pop(1)%tau_e / meta(1)%Ni !(s)
    pop(1)%Tr      = 100.d0 ! (s)

    !**** Laser variables ***
    ! mean velocity of the metastables
    vm   = sqrt(2.d0*kb*meta(0)%Tp*qok/mhe) 
    ! Laser frequency (1083 nm)
    Omeg = ppi*Vcel / lasr%Lwave
    ! Doppler Width due to metastable velocity distribution
    Dop  = Omeg * Vm / (ppi*Vcel)
    ! Laser coefficient 
    gammak = sqrt(pi) * fineS * fosc(1,3) * lasr%Is / (me*Omeg*Dop*lasr%Sec)
    !************************
    IF (lasr%OnOff.EQ.1 .and. Clock%SumDt.GT.lasr%Stime) THEN
       switch = 1
    END IF

    DO j = 1, Npop2
       !**** Update B_j sublevels due to radiative transitions (2P3 -> 2S3) ***
       Updens(2,j) = Updens(2,j) - Dt * meta(3)%Aij(1) * pop(2)%Ni(j)
       Do i = 1, Npop1
          !**** Update A_i sublevels due to radiative transitions (2P3 -> 2S3) ***
          Updens(1,i) = Updens(1,i) + Dt * meta(3)%Aij(1) * ( pop(2)%Ni(j)*Tij(j,i,1) )
          !**** Update A_i and B_j sublevels due to laser absorption ***
          IF (switch.EQ.1) THEN ! If laser activated (1)
             DO k = 1, lasr%Ntr
                IF ( int(Tij(j,i,3)).EQ.lasr%Ck(k) ) THEN ! Function of transitions Ck
                   IF ( int(Tij(j,i,2)).EQ.lasr%plz) THEN ! Function of polarization
                      nu_ij = Tij(j,i,1) * gammak
                      Updens(1,i) = Updens(1,i) + Dt * nu_ij * ( pop(2)%Ni(j) - pop(1)%Ni(i) )
                      Updens(2,j) = Updens(2,j) + Dt * nu_ij * ( pop(1)%Ni(i) - pop(2)%Ni(j) )
                      IF (nu_ij.GT.MaxR) MaxR = nu_ij
                   END IF
                END IF
             END DO    
          END IF
          !**** Update A_i sublevels due to metastability exchange between A_k sublevels ***
          IF (j.LE.Npop1) THEN
             Updens(1,i) = Updens(1,i) + Dt * ( lasr%Eij(i,j)+lasr%Fij(i,j)*pop(1)%polarz ) &
                  * pop(1)%Ni(j) / (18.d0*pop(1)%tau_e)
          END IF
       END DO
    END DO

    !**** Dn total after one loop in Boltz collisions ***
    pop(1)%Dn_o = meta(1)%Ni - pop(1)%Ntot
    pop(2)%Dn_o = meta(3)%Ni - pop(2)%Ntot
    !**** Update of sublevels densities due to Boltz collisions ***
    Updens(1,1:Npop1) = Updens(1,1:Npop1) + pop(1)%Dn_o / N1
    Updens(2,1:Npop2) = Updens(2,1:Npop2) + pop(2)%Dn_o / N2

    pol = 0.d0
    DO i = 1, Npop1
       !**** calculation of the third term in polarization's equation ***
       pol = pol + lasr%lamb(i) * pop(1)%Ni(i) / (3.d0*meta(1)%Ni)
       !**** Update A_i sublevels due to relaxation from collisions ***
       Updens(1,i) = Updens(1,i) + Dt * ( meta(1)%Ni/N1 - pop(1)%Ni(i) ) / pop(1)%T_relax
       !**** Update of the densities and the total densities of 2S3 and 2P3 ***
       pop(1)%Ni(i) = pop(1)%Ni(i) + Updens(1,i)
    END DO
    DO j = 1, Npop2
       !**** Update B_j sublevels due to relaxation from collisions ***
       Updens(2,j) = Updens(2,j) + Dt * (meta(3)%Ni/N2 - pop(2)%Ni(j) ) / pop(2)%T_relax
       !**** Update of the densities and the total densities of 2S3 and 2P3 ***
       pop(2)%Ni(j) = pop(2)%Ni(j) + Updens(2,j)
    END DO

    !**** Calculate the total density of metastable 2S3 and radiative state 2P3 ***
    tot1 = 0.d0 ; tot2 = 0.d0
    do i = 1, Npop2    
       if (i.LE.Npop1) tot1 = tot1 + pop(1)%Ni(i)
       tot2 = tot2 + pop(2)%Ni(i)
    END do
    meta(1)%Ni = tot1 ; meta(3)%Ni = tot2

    !**** Calculate the polarization of the Helium gas ***
    IF (mod(iter,1000) == 0 ) THEN
       E_meta = ABS(1.d0 - meta(1)%NStart / meta(1)%Ni)
       E_2P3  = ABS(1.d0 - meta(3)%NStart / meta(3)%Ni)
       IF (E_meta.LE. errmax .and. E_2P3.LE.errmax) THEN
          pop(1)%polarz = pop(1)%polarz +  (pop(1)%Te/20.d0) * &
               ((-pop(1)%polarz + pol)/pop(1)%Te )!- pop(1)%polarz/pop(1)%Tr)
          write(*,*) "polarization updated"
       END IF
       meta(1)%NStart = meta(1)%Ni
       meta(3)%NStart = meta(3)%Ni
    END IF

  END Subroutine Sublev_coll

END MODULE MOD_PUMP
