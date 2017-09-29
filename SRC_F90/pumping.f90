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

  Subroutine Sublev_coll (Clock,meta,pop,Neut,Tij,lasr,iter)
    !INTENT
    INTEGER      , INTENT(INOUT) :: iter
    TYPE(TIME)   , INTENT(IN) :: Clock
    TYPE(Laser)  , INTENT(IN) :: lasr
    REAL(DOUBLE) , DIMENSION(Npop2,Npop1,3) :: Tij
    TYPE(Excited), DIMENSION(2), INTENT(INOUT) :: pop
    TYPE(Species), DIMENSION(2), INTENT(IN) :: Neut
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    !LOCAL
    INTEGER :: i, j, k, switch, Niter
    REAL(DOUBLE) :: Dt, nu_ij, E_meta, Dt_p, E_max, E_pop
    REAL(DOUBLE) :: gammak, Dop, Omeg, vm
    REAL(DOUBLE) :: tot1, tot2, N1, N2, pol
    REAL(DOUBLE), DIMENSION(2,18) :: Updens
    !**************************
    switch = 0 ! Laser On (1) or Off (0) 
    N1 = real(Npop1) ; N2 = real(Npop2)
    Dt = Clock%Dt ; Updens(:,:)  = 0.d0

    E_max = 2d-4 ; Dt_p = 1d-6 ! E_max == percentage (limit) of variation of DAi/Dt.
    Niter = int(Dt_p / Dt) ! Calcul Niter for 1 micro-sec (== Dt_p).

    !pop(1)%T_relax = 1.d0/ (meta(1)%Damb / (sys%Ra/2.405d0)**2)
    !pop(2)%T_relax = 1.d0/ (3.2d06*1.33322*meta(0)%Prs) !(s avec P(mbar)) (cf. These Batz p. 30)
    pop(1)%tau_e   = 1.d0/ (3.75d6*1.33322*meta(0)%Prs) !(s avec P(mbar)) (cf. These Batz p. 30)
    pop(1)%Te      = meta(0)%Ni * pop(1)%tau_e / meta(1)%Ni !(s) (cf. Nacher 1985)
    pop(1)%Tr      = 100.d0 ! (s)

    !**** Laser variables ***
    ! mean velocity of the metastables
    vm   = sqrt(2.d0*kb*meta(0)%Tp*qok/mhe3) 
    ! Laser frequency (1083 nm)
    Omeg = ppi*Vcel / lasr%Lwave
    ! Doppler Width due to metastable velocity distribution
    Dop  = Omeg * Vm / (ppi*Vcel)
    ! Laser coefficient 
    gammak = sqrt(pi) * fineS * 0.5391d0 * lasr%Is / (me*Omeg*Dop*lasr%Sec)
    !************************
    IF (lasr%OnOff.EQ.1 .and. Clock%SumDt.GT.lasr%Stime) THEN
       switch = 1
    END IF

    DO j = 1, Npop2
       !**** Update B_j sublevels due to radiative transitions (2P3 -> 2S3) ***
       Updens(2,j) = Updens(2,j) - Dt * pop(2)%tau_rad * pop(2)%Ni(j)
       Do i = 1, Npop1
          !**** Update A_i sublevels due to radiative transitions (2P3 -> 2S3) ***
          Updens(1,i) = Updens(1,i) + Dt * pop(2)%tau_rad * ( pop(2)%Ni(j)*Tij(j,i,1) )
          !**** Update A_i and B_j sublevels due to laser absorption ***
          IF (switch.EQ.1) THEN ! If laser activated (1)
             DO k = 1, lasr%Ntr
                IF ( int(Tij(j,i,3)).EQ.lasr%Ck(k) ) THEN ! Function of transitions Ck
                   IF ( int(Tij(j,i,2)).EQ.lasr%plz) THEN ! Function of polarization
                      nu_ij = Tij(j,i,1) * gammak
                      Updens(1,i) = Updens(1,i) + Dt * nu_ij * ( pop(2)%Ni(j) - pop(1)%Ni(i) )
                      Updens(2,j) = Updens(2,j) + Dt * nu_ij * ( pop(1)%Ni(i) - pop(2)%Ni(j) )
                      IF (nu_ij.GT.MaxR) MaxR = nu_ij ! Variation of the time step (max of s-1)
                   END IF
                END IF
             END DO
             !**** Update A_i sublevels due to metastability exchange between A_k sublevels ***
             IF (j.LE.Npop1) THEN
                Updens(1,i) = Updens(1,i) + Dt * ( lasr%Eij(i,j)+lasr%Fij(i,j)*pop(1)%polarz ) &
                     * pop(1)%Ni(j) / (18.d0*pop(1)%tau_e)
             END IF
          END IF
       END DO
    END DO

    !**** Dn total after one loop in Boltz collisions ***
    pop(1)%Dn_o = meta(1)%Ni - pop(1)%Ntot
    pop(2)%Dn_o = meta(3)%Ni - pop(2)%Ntot

    pol = 0.d0
    DO i = 1, Npop1
       !**** Update of sublevels densities due to Boltz collisions *************
       !**** If Delta_n > 0 then mean allocations on sublevels Ai **************
       !**** IF Delta_n < 0 then density prorata depending on each sublevels ***
       IF (pop(1)%Dn_o.GE.0.d0) THEN
          Updens(1,i) = Updens(1,i) + pop(1)%Dn_o / N1
       ELSE
          Updens(1,i) = Updens(1,i) + pop(1)%Dn_o * (pop(1)%Ni(i)/pop(1)%Ntot)
       END IF
       !**** calculation of the third term in polarization's equation ***
       pol = pol + lasr%lamb(i) * pop(1)%Ni(i) / (3.d0*meta(1)%Ni)
       !**** Update of the densities and the total densities of 2S3 and 2P3 ***
       pop(1)%Ni(i) = pop(1)%Ni(i) + Updens(1,i)
    END DO
    DO j = 1, Npop2
       !**** Update of sublevels densities due to Boltz collisions for Bj sublevels ***
       IF (pop(2)%Dn_o.GE.0.d0) THEN
          Updens(2,j) = Updens(2,j) + pop(2)%Dn_o / N2
       ELSE
          Updens(2,j) = Updens(2,j) + pop(2)%Dn_o * (pop(2)%Ni(j)/pop(2)%Ntot)
       END IF
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

    !**** Calcul des populations des Ai avec P fixe.
    IF (lasr%OnOff.EQ.1)THEN
       pop(1)%polarz = 0.002 ! Tagada ...
       !**** Calcul du critere d'arret
       !E_pop = ABS(pop(1)%Ni(4)/meta(1)%NStart - 1d0)
       !IF (E_pop.LT.1d-08.and.iter.GT.5000) THEN
       IF (iter.EQ.Clock%NumIter-1) THEN
          OPEN(UNIT=919,File=TRIM(ADJUSTL(DirFile))//"A_i.dat",ACTION="WRITE",POSITION="APPEND",STATUS="UNKNOWN")
          WRITE(919,"(3ES12.3,6(ES19.10),6(ES19.10))") meta(0)%Prs,lasr%Is, pop(1)%polarz, &
               (pop(1)%Ni(i), i=1,6), (meta(i)%Ni, i=1,6) 
          CLOSE(919)
          !print*, "Critere d'arret = ", E_pop
          !iter = Clock%MaxIter
       END IF
       !**** Pour info sur l'evolution du critere d'arret:
       !IF (mod(iter,Niter*2) == 0 .and. switch==1) THEN
       !   print*, "Critere d'arret = ", E_pop, pop(1)%Ni(4), meta(1)%NStart, pop(1)%Ni(4)/meta(1)%NStart
       !END IF
       !meta(1)%NStart = pop(1)%Ni(4)     
       !    !**** Calculate the polarization of the Helium gas ***
       !    IF (mod(iter,Niter) == 0 .and. switch==1) THEN
       !       E_meta = ABS(pop(1)%Ni(1)/meta(1)%NStart - 1d0)
       !       IF (E_meta.LE. E_max) THEN
       !          !**** [pop(1)%Te/20.d0] == delta T from --> dP/dT
       !          pop(1)%polarz = pop(1)%polarz + (pop(1)%Te/20.d0) * ((-pop(1)%polarz + pol)/pop(1)%Te &
       !               - pop(1)%polarz/pop(1)%Tr) - (pop(1)%polarz * Neut(1)%UpDens / meta(0)%Ni) &
       !               * (pop(1)%Te/20.d0)/(Niter*Dt)
       !       END IF
       !       meta(1)%NStart = pop(1)%Ni(1)
       !       print*, "polariz : DT iter = ", Niter, " P= ", pop(1)%polarz, " err= ", E_meta
       !    END IF
       !
    END IF
  END Subroutine Sublev_coll

END MODULE MOD_PUMP
