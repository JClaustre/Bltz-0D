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
    INTEGER      , INTENT(INOUT) :: iter
    TYPE(TIME)   , INTENT(IN) :: Clock
    TYPE(Laser)  , INTENT(IN) :: lasr
    REAL(DOUBLE) , DIMENSION(Npop2,Npop1,3) :: Tij
    TYPE(Excited), DIMENSION(2), INTENT(INOUT) :: pop
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    !LOCAL
    INTEGER :: i, j, k, switch
    REAL(DOUBLE) :: Dt, nu_ij
    REAL(DOUBLE) :: gammak, Dop, Omeg, vm, Tp
    REAL(DOUBLE) :: tot1, tot2, N1, N2, pol
    REAL(DOUBLE), DIMENSION(2,18) :: Updens
    !****************************************************
    switch = 0 ! Laser On (1) or Off (0) --> don't modify here!
    N1 = real(Npop1) ; N2 = real(Npop2)
    Dt = Clock%Dt ; Updens(:,:)  = 0.d0

    !**** Definition of polarization rate
    pop(1)%tau_e   = 1.d0/ (3.75d6*1.33322*meta(0)%Prs) !(s avec P(mbar)) (cf. These Batz p. 30)
    pop(1)%Te      = meta(0)%Ni * pop(1)%tau_e / meta(1)%Ni !(s) (cf. Nacher 1985)
    pop(1)%Tr      = 1.0/ (Ngpl(1)%UpDens/(meta(0)%Ni*Clock%Dt))
    Tp = 6.876e-30 * meta(0)%Ni*ion(1)%Ni ! Milner et al (Nuclear Instruments and Methods in
                                          ! Physics Research A257 (1987) 286-290) 

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

          END IF
          !**** Update A_i sublevels due to metastability exchange between A_k sublevels ***
          IF (j.LE.Npop1) THEN
             Updens(1,i) = Updens(1,i) + Dt * ( lasr%Eij(i,j)+lasr%Fij(i,j)*pop(1)%polarz ) &
                  * pop(1)%Ni(j) / (18.d0*pop(1)%tau_e)
          END IF
          !*********************************************************************************
       END DO
    END DO

    !**** Dn total after one loop in Boltz collisions ***
    pop(1)%Dn_o = meta(1)%Ni - pop(1)%Ntot
    pop(2)%Dn_o = meta(3)%Ni - pop(2)%Ntot

    pol = 0.d0
    DO i = 1, Npop1
       !**** Update of sublevel densities due to Boltz collisions *************
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
       !**** Update of sublevel densities due to Boltz collisions for Bj sublevels ***
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
       pop(1)%polarz = 0 ! Tagada ...
       !**** Copy into files Ai populations.
       IF (iter.EQ.Clock%NumIter-2) THEN
          OPEN(UNIT=919,File=TRIM(ADJUSTL(DirFile))//"A_i.dat",ACTION="WRITE",POSITION="APPEND",STATUS="UNKNOWN")
          WRITE(919,"(3ES12.3,6(ES19.10),6(ES19.10))") meta(0)%Prs,lasr%Is, pop(1)%polarz, &
               (pop(1)%Ni(i), i=1,6), (meta(i)%Ni, i=1,6) 
          CLOSE(919)
          !iter = Clock%MaxIter
       END IF
       !*************************************
       
       !**** Calcul de la polarisationen fonction du temps.
       !    ELSE IF (lasr%OnOff.EQ.0) THEN
       !       !**** Calcul the (de)polarization of the Helium gas ***
       !       pop(1)%polarz = pop(1)%polarz + Clock%Dt * ((-pop(1)%polarz + pol)/pop(1)%Te &
       !            - pop(1)%polarz/pop(1)%Tr - pop(1)%polarz/Tp)
    END IF

  END Subroutine Sublev_coll

END MODULE MOD_PUMP
