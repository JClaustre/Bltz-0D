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

  Subroutine Sublev_coll (Clock,meta,pop,Tij)
    !INTENT
    TYPE(TIME)   , INTENT(IN) :: Clock
    REAL(DOUBLE) , DIMENSION(18,6,2) :: Tij
    TYPE(Excited), DIMENSION(2), INTENT(INOUT) :: pop
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    !LOCAL
    INTEGER :: i, j
    REAL(DOUBLE) :: Dt, tau_ij
    REAL(DOUBLE) :: tot1, tot2
    !**************************
    Dt = Clock%Dt

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
    pop(1)%Ni(:) = pop(1)%Ni(:) + pop(1)%Dn_o   / 6.d0
    pop(2)%Ni(:) = pop(2)%Ni(:) + pop(2)%Dn_tot / 18.d0


    DO j = 1, 18
       Do i = 1, 6
          nu_ij = 1./ (Tij(j,i,1))
          !**** Update A_i sublevel due to radiative transitions (2P3 -> 2S3) ***
          pop(1)%Ni(i) = pop(1)%Ni(i) + (pop(1)%Dn_rad/18.d0) * Tij(j,i,1)
          !**** Update A_i and B_j sublevel due to laser absorption ***
          pop(1)%Ni(i) = pop(1)%Ni(i) + Dt * nu_ij * ( pop(2)%Ni(j) - pop(1)%Ni(i) )
          pop(2)%Ni(j) = pop(2)%Ni(j) + Dt * nu_ij * ( pop(1)%Ni(i) - pop(2)%Ni(j) )
       END DO
    END DO

    !**** Update of the total density of 2S3 and 2P3 ***
    tot1 = 0.d0 ; tot2 = 0.d0
    do i = 1, 18    
       if (i.LE.6) tot1 = tot1 + pop(1)%Ni(i)
       tot2 = tot2 + pop(2)%Ni(i)
    END do
   
!    write(*,"(A,2ES15.6)") "Verif des densites : ", &
!         (1.d0- (meta(1)%Ni/tot1)), (1.d0 - (meta(3)%Ni/tot2))

  END Subroutine Sublev_coll

END MODULE MOD_PUMP
