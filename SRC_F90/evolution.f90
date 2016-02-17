!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 15/05/2015
! Objctv: Main Loop to performe eedf in time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MOD_EVOL

  USE F90_KIND  
  USE MOD_PARAM
  USE MOD_EXCIT
  USE MOD_XCHANGE
  USE MOD_RECOMB
  USE MOD_PENNASS
  USE MOD_IONIZ
  USE MOD_RADIFF
  USE MOD_CHAUF
  USE MOD_READ
  USE MOD_TPGAZ
  IMPLICIT NONE

  INTEGER :: XcDx = 0 ! 1 == equil | 0 == implic
  INTEGER :: IonX = 0 ! 1 == 50-50 | 0 == 100-0

CONTAINS
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE EVOLUTION ()
    ! LOCAL ******************************************************************!
    INTEGER :: i, j, k, l, Nnull=0                                            !
    INTEGER :: t1, t2, clock_rate                                             !
    REAL(DOUBLE) :: count1, count2, MaxDt                                     !
    REAL(DOUBLE) :: GenPwr                                                    !
    CHARACTER(LEN=250)::fileName                                              !
    count1 = 0.d0 ; count2 = 0.d0                                             !
    !*****************************                                            !
    Clock%NumIter = int( (Clock%SimuTime-Clock%SumDt) /Clock%Dt)              !
    write(*,"(2A,I10)") tabul, "Iterations in Time: ",  Clock%NumIter         !
    l = 0 ; k = 0                                                             !
    IF (Clock%Rstart == 0) THEN                                               !
       OPEN(UNIT=99,File="./datFile/evol.dat",ACTION="WRITE",STATUS="UNKNOWN")!
    ELSE IF (Clock%Rstart == 1) THEN                                          !
       OPEN(UNIT=99,File="./datFile/evol.dat",ACCESS="STREAM",ACTION="WRITE",STATUS="UNKNOWN")
    END IF                                                                    !
    !*************************************************************************!
    MaxDt  = 1.d-07 ! Maximum Time-Step allowed
    sys%IPowr = sys%Powr ! Keep Power init in memory
    GenPwr = 0.5d-6 ! Time constant to start the generator.

    !**** MAIN LOOP ***************************
    DO WHILE (Clock%SumDt .LT. Clock%SimuTime)
       if (l == 200) CALL System_clock (t1, clock_rate)
       l = l + 1
       meta(:)%Updens = 0.d0 ; ion(:)%Updens = 0.d0
       !**** Update Time-step
       IF (1.d0/MaxR .GE. 1.0d-12 .and. l .GT.1) THEN
          Clock%Dt = 1.0d0 / (MaxR*3.)
          IF (Clock%Dt .GT. MaxDt) Clock%Dt = MaxDt ! Maximum Time-Step allowed
          MaxR = 0.d0
       END IF
       !**** Check if there's NaN propagation ... probably due to large Dt (change MaxDt).
       IF (ISnan(MaxR) .or. isNaN(elec%Ni) ) THEN 
           print*,""; print*," NaN ! Pobleme in time step?", elec%Ni, MaxR ; Stop 
       END IF
       !*************************************

       !**** Neutral temperature calculation
       CALL TP_Neutral (sys, elec, meta, Tg_p)
       !**** Increase Power exponantially function of time
       sys%Powr = sys%IPowr * (1.d0 - exp( -real(k*Clock%Dt) / GenPwr) )
       k = k+1
       !**** Heat + Elas + Fk-Pl
       CALL Heating (sys,meta, U, F)
       CALL Elastic      (sys,meta, U, F)
       CALL FP           (sys, elec, F, U)
       !**** Excit + De-excit
       SELECT CASE (XcDx)
       CASE (0) ; CALL Exc_Impli     (sys, meta, U, F, diag)
       CASE (1) ; CALL Exc_Equil     (sys, meta, U, F, diag)
       CASE DEFAULT ; CALL Exc_Begin (sys, meta, U, F, diag)
       END SELECT
       !**** De-excit Dimer molecule (He2*)
       IF (NumIon == 3) CALL Dexc_Dimer (sys, U, ion, F, diag)
       !**** Ioniz He+
       SELECT CASE (IonX)
       CASE (0) ; CALL Ioniz_100    (sys, meta, U, F, diag)
       CASE (1) ; CALL Ioniz_50     (sys, meta, U, F, diag)
       CASE DEFAULT ; CALL Ioniz_100(sys, meta, U, F, diag)
       END SELECT
       !**** Ioniz dimer 
       IF (NumIon == 3) CALL Ioniz_Dimer100 (sys, ion, U, F)
       !**** Disso Recombination
       CALL Recomb       (sys, meta, U, F, Diag)
       !**** 3 Body ionic conversion
       CALL Conv_3Body   (meta, ion)
       !**** Penning + Associative ioniz
       CALL Penn_Assoc   (sys, meta, U, F, Diag)
       !**** Radiative transfert
       CALL Radiat       (sys, meta, Fosc, Diag)
       !**** Diffusion
       CALL Diffuz       (sys, meta, ion,elec,F,U, diag)
       !**** L-Exchange
       CALL l_change     (meta, K_ij)
       !*************************************

       !**** Update densities (Ion + Excited)
       Nnull = 0
       do i = 1, NumMeta
          meta(i)%Ni = meta(i)%Ni + meta(i)%Updens
          if (i .LE. NumIon) ion(i)%Ni = ion(i)%Ni + ion(i)%Updens
          IF (meta(i)%Ni < 0.d0) THEN
             Nnull = Nnull + 1
             meta(i)%Ni = 0.d0
          END IF
       END do
       !**** UpDate Density (electron) + Tpe
       elec%Ni = 0.d0 ; elec%Tp = 0.d0; elec%J = 0.d0
       DO i = 1, sys%nx 
          elec%Ni = elec%Ni + ( F(i) * dsqrt(U(i)) * sys%Dx )
          elec%Tp = elec%Tp + ( F(i) * dsqrt(U(i)**3) * sys%Dx )
          elec%J  = elec%J  + ( F(i) * U(i) * gama*qe * sys%Dx )
       END DO
       elec%Tp = elec%Tp * 0.6667d0 / elec%Ni

       !**** Evaluate Calculation Time
       if (l == 300) CALL System_clock (t2, clock_rate)
       if (l == 300) CALL LoopTime(t1, t2, clock_rate, Clock%NumIter)
       !**** UpDate Simulation Time
       Clock%SumDt = Clock%SumDt + Clock%Dt

       !**** WRITE IN SHELL ******************************************************!
       write(*,"(2A,F8.3,A,F5.1,A,I7,A,ES9.3,A,F5.1,A,I4,A)",advance="no") tabul,&!
       "Time in simulation: ", (Clock%SumDt*1e6), " μs | achieved: ",&            !
            Clock%SumDt/Clock%SimuTime*100.d0, "% [ it = ", l, " | Dt = ",&       !
            Clock%Dt, " Pwr(%): ", (sys%Powr*100./sys%IPowr), "]", Nnull, " \r"   !
                                                                                  !
       IF (modulo(l,int(Clock%SimuTime/Clock%Dt)/10) == 0) then                   !
          write(*,"(2A,F7.2,A,4ES13.4,A,ES10.2,A,F7.1)") tabul, "Time :", &          !
               (Clock%SumDt*1e6), " μs", meta(1)%Ni*1d-06, meta(3)%Ni*1d-06,&     !
               ion(1)%Ni*1d-06, ion(2)%Ni*1d-06, " | E/N(Td)", &
               (sys%E/meta(0)%Ni)/1d-21, " Tg(K)", meta(0)%Tp*qok
       END IF                                                                     !
       !**************************************************************************!

       !**** WRITE IN FILES (Function of TIME) (density in cm^-3) ****************!
       IF (modulo(l,100)==0) Then                                                 !
          SELECT CASE (NumIon)                                                    !
          CASE (3)                                                                !
             write(99,"(48ES15.4)") (Clock%SumDt*1e6), elec%Tp, elec%Ni*1d-06,&   !
                  ion(1)%Ni*1d-06, ion(2)%Ni*1d-06, ion(NumIon)%Ni*1d-06, &       !
                  (meta(i)%Ni*1d-06,i=1,NumMeta)                                  !
          CASE DEFAULT                                                            !
             write(99,"(47ES15.4)") (Clock%SumDt*1e6), elec%Tp, elec%Ni*1d-06,&   !
                  ion(1)%Ni*1d-06, ion(2)%Ni*1d-06, (meta(i)%Ni*1d-06,i=1,NumMeta)!
          END SELECT                                                              !
          !****                                                                   !
          !**** WRITE IN FILES (density in cm^-3) ****                            !
          OPEN(UNIT=98,File="./datFile/density.dat",ACTION="WRITE",STATUS="UNKNOWN")
          DO i = 1, NumMeta                                                       !
             write(98,"(I3,A,F10.4,ES15.4)") i, meta(i)%Name, meta(i)%En, meta(i)%Ni*1d-06
          END DO                                                                  !
          DO i = 1, NumIon                                                        !
             write(98,"(I3,A,F10.4,ES15.4)") i, ion(i)%Name, ion(i)%En, ion(i)%Ni*1d-06 
          END DO                                                                  !
          CLOSE(98)                                                               !
                                                                                  !
       END IF                                                                     !
       !**** WRITE IN FILES (Time Dependent) (EEDF in cm^-3) *********************!
       IF ( modulo(l,int(Clock%TRstart/Clock%Dt)) == 0 ) THEN                     ! 
          write(fileName,"('./datFile/F_evol_',I3.3,'.dat')") j                   !
          OPEN(UNIT=90,File=TRIM(ADJUSTL(fileName)),ACTION="WRITE",STATUS="UNKNOWN")
          DO i = 1, sys%nx                                                        !
             write(90,"(2ES15.6)") real(i)*sys%Dx, F(i)                           !
          END DO                                                                  !
          CLOSE(90)                                                               !
          j = j+1                                                                 !
       END IF                                                                     !
       !**************************************************************************!

       !**** WRITE IN FILES (Time Dependent) (Restart files) *********************!
       IF ( modulo(l,int(Clock%TRstart/Clock%Dt)) == 0 ) THEN                     !
          CALL Rstart_SaveFiles (sys, Clock, ion, elec, meta, F)                  !
       END IF                                                                     !
       !**************************************************************************!

       !**** Max iteration number (Else exit loop)********************************!
       IF (l .GE. Clock%MaxIter) EXIT                                             !
    END DO                                                                        !
    !****End of MAIN LOOP ********************************************************!

    CLOSE(99)

    CALL Consv_Test(sys, U, F, Diag, consv)
    CALL Write_Out1D( F,  "F_final.dat")
    CALL Write_Out1D( Tg_p, "Tg_p.dat")
    write(*,"(2A,F6.2,A)") tabul,"--> Simulation Time : ", real(Clock%SumDt/1.0d-6), " μs"

  END SUBROUTINE EVOLUTION
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE Consv_Test(sys, U, Fi, Diag, consv)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN) :: U
    TYPE(Diagnos), DIMENSION(:) , INTENT(INOUT) :: diag
    REAL(DOUBLE) , DIMENSION(:) , INTENT(IN) :: Fi
    REAL(DOUBLE) , DIMENSION(2) , INTENT(INOUT) :: consv
    !LOCAL
    INTEGER      :: i
    REAL(DOUBLE) :: Coef
    ! **** Particle Conservation : Σ f(i).U(i)^½.ΔU ************************************!
    Coef = 0.d0                                                                         !
    DO i = 1, sys%nx                                                                    !
       Coef = Coef + Fi(i)*U(i)**(0.5d0) * sys%Dx                                       !
    END DO                                                                              !
    elec%Ni = Coef                                                                      !
    write(*,"(3A,ES15.4)") tabul,"Density (cm-3): ", elec%Name  , elec%Ni*1.d-6         !
    write(*,"(3A,ES15.4)") tabul,"Density (cm-3): ", ion(1)%Name, ion(1)%Ni*1.d-6       !
    write(*,"(3A,ES15.4)") tabul,"Density (cm-3): ", ion(2)%Name, ion(2)%Ni*1.d-6       !
    Coef = ABS(1.0d0 - (elec%Ni/(ion(1)%Ni+ion(2)%Ni)))                                 !
    write(*,"(2A,ES15.4)") tabul, "Partcl Error : ", Coef                               !
    write(*,"(2A,ES15.4)") tabul, "Particle Loss due to diffusion: ", diag(9)%Tx/elec%Ni!
    write(*,"(2A,ES15.4)") tabul, "Particle Loss due to recombina: ", diag(8)%Tx/elec%Ni!
    write(*,"(2A,ES15.4)") tabul, "Particle Gain due to AssoIoniz: ", diag(6)%Tx/elec%Ni!
    write(*,"(2A,ES15.4)") tabul, "Particle Gain due to Penning  : ", diag(5)%Tx/(2.d0*elec%Ni)
    write(*,"(2A,ES15.4)") tabul, "Particle Gain due to ionizatio: ", diag(2)%Tx/elec%Ni!
    !***********************************************************************************!

    !**** Energy Conservation : Σf(i).U(i)^3/2.ΔU + N*.Eij *****************************!
    Coef = 0.d0                                                                         !
    DO i = 1, sys%nx                                                                    !
       Coef = Coef + Fi(i)* U(i)**(1.5d0) * sys%Dx                                      !
    END DO                                                                              !
    elec%Tp = Coef*0.66667d0/elec%Ni                                                    !
                                                                                        !
                                                                                        !
    !**** (1)  **** Add energy due to excit/de-excit processes                          !
    Coef = Coef + (Diag(1)%EnLoss-Diag(1)%EnProd)                                       !
    !**** (2)  **** Energy Loss due to ionization processes                             !
    Coef = Coef + Diag(2)%EnLoss                                                        !
    !**** (5)  **** Energy gain due to Penning processes                                !
    Coef = Coef - Diag(5)%EnProd                                                        !
    !**** (6)  **** Energy gain due to Associative processes                            !
    Coef = Coef - Diag(6)%EnProd                                                        !
    !**** (7)  **** Energy loss due to Ionic conversion                                 !
    !Coef = Coef + diag(7)%EnLoss(NumMeta+1)                                            !
    !**** (8)  **** Energy loss due to recombination processes                          !
    Coef = Coef + diag(8)%EnLoss                                                        !
    !**** (9)  **** Energy loss due to diffusion process                                !
    Coef = Coef + Diag(9)%EnLoss                                                        !
    !**** (10-1) **** Energy gain due to heating process                                !
    Coef = Coef - Diag(10)%EnProd                                                       !
    !**** (10-2) **** Energy loss due to Elastic collisions                             !
    Coef = Coef + Diag(11)%EnLoss                                                       !
    !**** (12) **** Energy gain due to De-excit collisions from Dimer molecule          !
    Coef = Coef - Diag(12)%EnProd                                                       !
    !**** (13) Energy gain and loss due to ionization/3-body recombination from Dimer   !
    Coef = Coef + (Diag(13)%EnLoss-Diag(13)%EnProd)                                     !
    !**** (14) Energy gain due to Penning reaction btween Dimer and metastable          !
    Coef = Coef - Diag(14)%EnProd                                                       !
    !**** (3)=radiative trans | (4)=l-xchnge reaction                                   !
    !**** (7)=3 body convert  |                                                         !
    !***********************************************************************************!

    elec%En = Coef
    write(*,"(2A,2ES15.4)") tabul, "Gain Power in Heat : ", Diag(10)%EnProd * qe/(clock%Dt*clock%NumIter),&
         Diag(10)%EnProd * qe * sys%volume/(clock%Dt*clock%NumIter)
    write(*,"(2A,2ES15.4)") tabul, "Loss Power in Elast: ", Diag(11)%EnLoss * qe/(clock%Dt*clock%NumIter),&
         Diag(11)%EnLoss * qe * sys%volume/(clock%Dt*clock%NumIter)

    Coef = ABS(1.0d0 - elec%En/consv(2))
    write(*,"(2A,ES15.4)") tabul, "Energy Error : ", Coef
    write(*,"(A,2(A,F5.2),A)") tabul, "Tpe (eV): | Init ", consv(2)*0.66667d0/consv(1),&
         " | Final ", elec%Tp, " |"

  END SUBROUTINE Consv_Test
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE OutPutMD (sys, meta, ion, elec, diag, consv)
    !INTENT
    TYPE(SysVar) , INTENT(IN) :: sys
    TYPE(Diagnos), DIMENSION(:) , INTENT(IN) :: diag
    TYPE(Species), INTENT(IN) :: elec
    REAL(DOUBLE) , DIMENSION(2) , INTENT(IN) :: consv
    TYPE(Species), DIMENSION(NumIon)   , INTENT(IN) :: ion
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(IN) :: meta
    !LOCAL
    REAL(DOUBLE) :: ne, Dt
    ne = elec%Ni ; Dt = Clock%Dt

    OPEN(UNIT=99,File="./datFile/Output.md",ACTION="WRITE",STATUS="UNKNOWN")
    write(99,"(A)")"     .-.     .-.     .-.     .-.     .-.     .-.     .-.    "
    write(99,"(A)")"    .'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `. "
    write(99,"(A)")"               _         "
    write(99,"(A)")"              [ ] Hope you enjoyed, "
    write(99,"(A)")"             (* *) /  look forward to seeing you again!     "
    write(99,"(A)")"              |>|                  -C0PO-                     "
    write(99,"(A)")"           __/===\__     "
    write(99,"(A)")"          //| o=o |\\    "
    write(99,"(A)")"        <]  | o=o |  [>  "
    write(99,"(A)")"            \=====/      "
    write(99,"(A)")"           / / | \ \     "
    write(99,"(A)")"          <_________>    "
    write(99,"(A)")"     .-.     .-.     .-.     .-.     .-.     .-.     .-.    "
    write(99,"(A)")"    .'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `. "
    write(99,"(A)")""

    write(99,"(A)") ""
    write(99,"(A)") "SYSTEM PARAMETERS"
    write(99,"(A)") "--------------------"
    write(99,"(A,I6)")      "* Number of grid Cells : ", sys%Nx
    write(99,"(A,F8.2)")    "* Energy Max (ev)      : ", sys%Emx
    write(99,"(A,ES10.2)")  "* Energy step (Δu)    : ", sys%Dx
    write(99,"(A,F8.2)")    "* Cylinder radius (cm) : ", sys%Ra*1d2
    write(99,"(A,F8.2)")    "* E/N (Td) : ", (sys%E/meta(0)%Ni) * 1d21
    write(99,"(A,F8.2)")    "* E Field (V/cm) : ", sys%E*1d-2
    write(99,"(A,ES11.3)")  "* heating frequency (Hz) : ", sys%Freq / (2*pi)
    write(99,"(A,2ES11.3)") "* Power (W/cm3) | (W): ", sys%Powr*1d-6, sys%Powr * sys%volume

    write(99,"(A)") ""
    write(99,"(A)") "TIME PARAMETERS"
    write(99,"(A)") "-----------------"
    write(99,"(A,ES11.2)")  "* Time step (Δt) : ", Clock%Dt
    write(99,"(A,I10)")     "* Number of iterations : ", Clock%NumIter
    write(99,"(A,F6.2)")    "* Time Simulation (μs): ", Clock%SimuTime*1.d6
    write(99,"(A,3(I3,A))") "* Elapsed Time in CPU : ", int(Clock%Hours),"H ", &
         int(Clock%Minutes), " Min ", int(Clock%Seconds), " sec" 
    write(99,"(A)") ""
    write(99,"(A)") "NEUTRAL GAS PARAMETERS"
    write(99,"(A)") "--------------------"
    write(99,"(A,ES11.3)")  "* Gas density (cm-3): ", meta(0)%Ni*1d-6
    write(99,"(A,F8.2)")    "* Gas Pressure (Torr): ", meta(0)%Prs
    write(99,"(2(A,F7.2))") "* Gas Temperature (°K | eV): ", meta(0)%Tp*qok, " | ", meta(0)%Tp
    write(99,"(A)") ""
    write(99,"(A)") "ELEC | IONS PARAMETERS"
    write(99,"(A)") "--------------------"
    write(99,"(A,ES15.6)") "* Electron   Density (cm-3) : ", elec%Ni*1d-6
    write(99,"(A,ES15.6)") "* Ion [He+]  Density (cm-3) : ", ion(1)%Ni*1d-6
    write(99,"(A,ES15.6)") "* Ion [He2+] Density (cm-3) : ", ion(2)%Ni*1d-6
    write(99,"(/,A)") "-------------------------------------------------------"
    write(99,"(A,ES15.6)") "* Error partcl : ", ABS(1.d0-elec%Ni/(ion(1)%Ni+ion(2)%Ni))
    write(99,"(A,ES15.6)") "* Error energy : ", ABS(1.d0-elec%En/consv(2))
    write(99,"(2(A,F8.2))") "* Electron Tp° (eV) : init ", consv(2)*0.66667d0/consv(1),&
         " | Final ", elec%Tp
    write(99,"(/,A)") "-------------------------------------------------------"
    write(99,"(A,ES11.3)") "* Electron  mobility (cm2.V-¹.s-¹) : ", elec%mobl
    write(99,"(A,ES11.3)") "* ion[He+]  mobility (cm2.V-¹.s-¹) : ", ion(1)%mobl
    write(99,"(A,ES11.3)") "* ion[He2+] mobility (cm2.V-¹.s-¹) : ", ion(2)%mobl
    write(99,"(A,ES11.3)") "* Electron  Free Diff (cm².s-¹) : ", elec%Dfree
    write(99,"(A,ES11.3)") "* ion[He+]  Free Diff (cm².s-¹) : ", ion(1)%Dfree
    write(99,"(A,ES11.3)") "* ion[He2+] Free Diff (cm².s-¹) : ", ion(2)%Dfree
    write(99,"(A,ES10.2)") "* Ionization degree (Ne/Ng) : ", elec%Ni/meta(0)%Ni
    write(99,"(A)") ""
    write(99,"(A)") "ELECTRON | IONS BALANCE"
    write(99,"(A)") "--------------------"
    write(99,"(A)") "### Power balance :"
    write(99,"(A,ES15.6)") "* Gain Heat :  ", Diag(10)%EnProd * qe/(ne*Dt*clock%NumIter)
    write(99,"(A,ES15.6)") "* Gain Penn :  ", Diag(5)%EnProd  * qe/(ne*Dt*clock%NumIter)
    write(99,"(A,ES15.6)") "* Gain Asso :  ", Diag(6)%EnProd  * qe/(ne*Dt*clock%NumIter)
    write(99,"(A,ES15.6)") "* Gain Exct :  ", Diag(1)%EnProd  * qe/(ne*Dt*clock%NumIter)

    write(99,"(/,A)") "-------------------------------------------------------"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Elas :  ", Diag(11)%EnLoss * qe/(ne*Dt*clock%NumIter)&
         , Diag(11)%EnLoss*100.d0/Diag(10)%EnProd, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Recb :  ", Diag(8)%EnLoss  * qe/(ne*Dt*clock%NumIter)&
         , Diag(8)%EnLoss*100.d0/Diag(10)%EnProd, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Ionz :  ", Diag(2)%EnLoss  * qe/(ne*Dt*clock%NumIter)&
         , Diag(2)%EnLoss*100.d0/Diag(10)%EnProd, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Exct :  ", Diag(1)%EnLoss  * qe/(ne*Dt*clock%NumIter)&
         , Diag(1)%EnLoss*100.d0/Diag(10)%EnProd, " %"
    write(99,"(A,ES15.6,F8.2,A,/)") "* Loss Diff :  ", Diag(9)%EnLoss  * qe/(ne*Dt*clock%NumIter)&
         , Diag(9)%EnLoss*100.d0/Diag(10)%EnProd, " %"

    write(99,"(A)") "### Particle balance :"
    write(99,"(A,ES15.6)") "* Gain ioniz : ", Diag(2)%Tx/ne
    write(99,"(A,ES15.6)") "* Gain Assoc : ", Diag(6)%Tx/ne
    write(99,"(A,ES15.6)") "* Gain Penng : ", Diag(5)%Tx/(2.d0*ne)
    write(99,"(/,A)") "-------------------------------------------------------"
    write(99,"(A,ES15.6)") "* Loss recmb : ", Diag(8)%Tx/ne
    write(99,"(A,ES15.6)") "* Loss diffz : ", Diag(9)%Tx/ne
    write(99,"(/,A)") "-------------------------------------------------------"
    write(99,"(A,2ES15.4)")"* Gain elec total (Pwr | Prtcl) : ", (qe/(ne*Dt*clock%NumIter) ) * &
         (Diag(10)%EnProd+Diag(5)%EnProd+Diag(6)%EnProd+Diag(1)%EnProd), &
         (Diag(2)%Tx+Diag(6)%Tx+Diag(5)%Tx/(2.d0))/ne
    write(99,"(A,2ES15.4)")"* Loss elec total (Pwr | Prtcl) : ", (qe/(ne*Dt*clock%NumIter) ) * &
         (Diag(11)%EnLoss+Diag(8)%EnLoss+Diag(2)%EnLoss+Diag(1)%EnLoss+Diag(9)%EnLoss), &
         (Diag(8)%Tx+Diag(9)%Tx) / ne
    write(99,"(A)") " "
    write(99,"(A)")"![Zozor](http://uploads.siteduzero.com/files/420001_421000/420263.png)"
    Close(99)
  END SUBROUTINE OutPutMD
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
END MODULE MOD_EVOL
