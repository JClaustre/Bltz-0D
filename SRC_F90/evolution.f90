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
  
  REAL(DOUBLE), PRIVATE :: Res
CONTAINS
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
  SUBROUTINE EVOLUTION ()
    ! LOCAL ******************************************************************!
    INTEGER :: i, j, k, l, m                                                  !
    INTEGER :: t1, t2, clock_rate                                             !
    REAL(DOUBLE) :: count1, count2, MxDt                                      !
    REAL(DOUBLE) :: Cgen, Res, Post_D, PwE                                    !
    count1=0.d0 ; count2=0.d0 ; Res=0.d0                                      !
    !*****************************                                            !
    Clock%NumIter = int( (Clock%SimuTime-Clock%SumDt) /Clock%Dt)              !
    write(*,"(2A,I10)") tabul, "Iterations in Time: ",  Clock%NumIter         !
    l = 0 ; k = 0                                                             !
    !*************************************************************************!
    sys%IPowr = sys%Powr ! Keep Power init in memory
    Cgen   = 1d-03 ! Time factor for external source.
    Post_D = 50d-2 ! Time to ignitiate post_discharge (micro-sec)
    MxDt   = 2d-10 ! Maximum time-step
    !**** MAIN LOOP ***************************

    !**** MAIN LOOP ***************************
    DO WHILE (Clock%SumDt .LT. Clock%SimuTime)
       if (l == 200) CALL System_clock (t1, clock_rate)
       l = l + 1
       meta(0:NumMeta)%Updens = 0.d0 ; ion(:)%Updens = 0.d0

       !**** Neutral temperature calculation
       !CALL TP_Neutral (sys, elec, meta, OneD)

       !**** Increase Power exponantially function of time
       CALL POWER_CONTROL (Clock, sys, meta, U, F, Post_D, Cgen,l)

       !**** Heat + Elas + Fk-Pl
       IF (sys%Pcent.NE.0.d0) CALL Heating (sys,meta, U, F)
       CALL Elastic      (sys,meta, U, F)
       CALL FP           (sys, elec, F, U)
       !**** Ioniz He+
       SELECT CASE (IonX)
       CASE (1) ; CALL Ioniz_50     (sys, meta, U, F, diag)
       CASE DEFAULT ; CALL Ioniz_100(sys, meta, U, F, diag)
       END SELECT
       !**** Ioniz dimer 
       IF (NumIon == 3) CALL Ioniz_Dimer100 (sys, ion, U, F, diag)
       !**** Disso Recombination
       CALL Recomb       (sys, meta, U, F, Diag)
       !**** 3 Body ionic conversion
       CALL Conv_3Body   (meta, ion)
       !**** Penning + Associative ioniz
       CALL Penn_Assoc   (sys, meta, U, F, Diag)
       !**** Radiative transfert
       CALL Radiat       (sys, meta, Fosc, Diag)
       !**** Diffusion
       CALL Diffuz_Gaine (sys, meta, ion,elec,F,U, diag)
       !**** Excit + De-excit
       SELECT CASE (XcDx)
       CASE (1) ; CALL Exc_Equil     (sys, meta, U, F, diag)
       CASE (2) ; CALL Exc_Begin (sys, meta, U, F, diag)
       CASE DEFAULT ; CALL Exc_Impli     (sys, meta, U, F, diag)
       END SELECT
       !**** De-excit Dimer molecule (He2*)
       IF (NumIon == 3) CALL Dexc_Dimer (sys, U, ion, F, diag)
       !**** L-Exchange
       CALL l_change     (meta, K_ij)
       !*************************************

       !**** UpDateS
       CALL CHECK_AND_WRITE (Clock, sys, meta, elec, ion, F, diag, l, MxDt)

       !**** Evaluation of Calculation Time
       if (l == 300) CALL System_clock (t2, clock_rate)
       if (l == 300) CALL LoopTime(t1, t2, clock_rate, Clock%NumIter)
       !**** UpDate Simulation Time
       Clock%SumDt = Clock%SumDt + Clock%Dt
       !**** Evaluation of Metastable and 2^3P rates
       !**** Total rate is made in diag(10) and superelas and inelas in diag(15)
       diag(10)%InM1=0.d0 ; diag(10)%InM2=0.d0 ; diag(10)%OutM1=0.d0 ; diag(10)%OutM2=0.d0
       DO i = 1, 14
          IF (i.NE.10) THEN 
             diag(10)%InM1  = diag(10)%InM1  + diag(i)%InM1
             diag(10)%InM2  = diag(10)%InM2  + diag(i)%InM2
             diag(10)%OutM1 = diag(10)%OutM1 + diag(i)%OutM1
             diag(10)%OutM2 = diag(10)%OutM2 + diag(i)%OutM2
          END IF
      END DO
      
       IF ( modulo(l,100) == 0 ) THEN
          !**** ALL rates
          IF (l == 100) THEN
             OPEN(UNIT=92,File=TRIM(ADJUSTL(DirFile))//"rates.dat",ACTION="WRITE",STATUS="UNKNOWN")
             write(92,"(15ES15.6)") Clock%SumDt*1d6, (diag(i)%Tx(1), i=1,14)
             OPEN(UNIT=93,File=TRIM(ADJUSTL(DirFile))//"rates_bis.dat",ACTION="WRITE",STATUS="UNKNOWN")
             write(93,"(ES15.6, 2(14F5.1))") Clock%SumDt*1d6, (diag(i)%Tx(2), diag(i)%Tx(3), i=1,14)
             OPEN(UNIT=94,File=TRIM(ADJUSTL(DirFile))//"MEOP_rates.dat",ACTION="WRITE",STATUS="UNKNOWN")
             write(94,"(ES15.6,4(15ES13.5))") Clock%SumDt*1d6, (diag(i)%InM1, diag(i)%OutM1, diag(i)%InM2, diag(i)%OutM2, i=1,15)
          ELSE
             OPEN(UNIT=92,File=TRIM(ADJUSTL(DirFile))//"rates.dat",ACTION="WRITE",STATUS="UNKNOWN",ACCESS="Append")
             write(92,"(15ES15.6)") Clock%SumDt*1d6, (diag(i)%Tx(1), i=1,14)
             OPEN(UNIT=93,File=TRIM(ADJUSTL(DirFile))//"rates_bis.dat",ACTION="WRITE",STATUS="UNKNOWN",ACCESS="Append")
             write(93,"(ES15.6, 2(14F5.1))") Clock%SumDt*1d6, (diag(i)%Tx(2), diag(i)%Tx(3) , i=1,14)
             OPEN(UNIT=94,File=TRIM(ADJUSTL(DirFile))//"MEOP_rates.dat",ACTION="WRITE",STATUS="UNKNOWN",ACCESS="Append")
             write(94,"(ES15.6,4(15ES13.5))") Clock%SumDt*1d6, (diag(i)%InM1, diag(i)%OutM1, diag(i)%InM2, diag(i)%OutM2, i=1,15)
          END IF
          CLOSE(92)
          CLOSE(93)
          CLOSE(94)
          !**** 23S Metastable and 2^3P rates ONLY
          IF (l == 100) THEN
             OPEN(UNIT=92,File=TRIM(ADJUSTL(DirFile))//"MetaTx.dat",ACTION="WRITE",STATUS="UNKNOWN")
             write(92,"(15ES15.6)") Clock%SumDt*1d6, (diag(i)%TxTmp(1), i=1,14)
             OPEN(UNIT=93,File=TRIM(ADJUSTL(DirFile))//"MetaTx_bis.dat",ACTION="WRITE",STATUS="UNKNOWN")
             write(93,"(ES15.6, 2(14F5.1))") Clock%SumDt*1d6, (diag(i)%TxTmp(2), diag(i)%TxTmp(3), i=1,14)
          ELSE
             OPEN(UNIT=92,File=TRIM(ADJUSTL(DirFile))//"MetaTx.dat",ACTION="WRITE",STATUS="UNKNOWN",ACCESS="Append")
             write(92,"(15ES15.6)") Clock%SumDt*1d6, (diag(i)%TxTmp(1), i=1,14)
             OPEN(UNIT=93,File=TRIM(ADJUSTL(DirFile))//"MetaTx_bis.dat",ACTION="WRITE",STATUS="UNKNOWN",ACCESS="Append")
             write(93,"(ES15.6, 2(14F5.1))") Clock%SumDt*1d6, (diag(i)%TxTmp(2), diag(i)%TxTmp(3) , i=1,14)
          END IF
          CLOSE(92)
          CLOSE(93)
       END IF
       !**************************************************************************!


       !**** Max iteration number (Else exit loop)********************************!
       IF (l .GE. Clock%MaxIter) EXIT                                             !
    END DO                                                                        !
    !****End of MAIN LOOP ********************************************************!

    !****End of MAIN LOOP ********************************************************!
    Clock%NumIter = l
    
    CALL Consv_Test(sys, U, F, Diag, consv)
    CALL Write_Out1D( F, "F_final.dat")
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
    write(*,"(2A,ES15.4)") tabul, "Particle Loss due to diffusion: ", diag(9)%SumTx/elec%Ni!
    write(*,"(2A,ES15.4)") tabul, "Particle Loss due to recombina: ", diag(8)%SumTx/elec%Ni!
    write(*,"(2A,ES15.4)") tabul, "Particle Gain due to AssoIoniz: ", diag(6)%SumTx/elec%Ni!
    write(*,"(2A,ES15.4)") tabul, "Particle Gain due to Penning  : ", diag(5)%SumTx/(2.d0*elec%Ni)
    write(*,"(2A,ES15.4)") tabul, "Particle Gain due to ionizatio: ", diag(2)%SumTx/elec%Ni!
    !***********************************************************************************!
    
    !**** Energy Conservation : Σf(i).U(i)^3/2.ΔU + N*.Eij *****************************!
    Coef = 0.d0                                                                         !
    DO i = 1, sys%nx                                                                    !
       Coef = Coef + Fi(i)* U(i)**(1.5d0) * sys%Dx                                      !
    END DO                                                                              !
    elec%Tp = Coef*0.66667d0/elec%Ni                                                    !
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
    Coef = Coef - Diag(4)%EnProd                                                        !
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
    write (*,"(2A,ES15.4,A,ES15.4)") tabul, "Absorbed Power [W/cm3]: ", &
         sys%Pwmoy*1d-6, ' Emoy(V/cm):', sys%Emoy*1d-2/Clock%NumIter

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
    INTEGER :: i
    REAL(DOUBLE) :: ne, Dt, gainE, gainP, SumMeta
    ne = elec%Ni ; Dt = Clock%Dt
    gainE = (Diag(10)%EnProd+Diag(5)%EnProd+Diag(6)%EnProd + &
         Diag(1)%EnProd+Diag(12)%EnProd+Diag(4)%EnProd)
    gainP = Diag(2)%SumTx + Diag(5)%SumTx + Diag(6)%SumTx + Diag(13)%SumTx + Diag(4)%SumTx 
    Do i = 1, NumMeta
       SumMeta = SumMeta + meta(i)%Ni
    END Do
    
    OPEN(UNIT=99,File=TRIM(ADJUSTL(DirFile))//"Output.md",ACTION="WRITE",STATUS="UNKNOWN")
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
    write(99,"(A,F8.3)")    "* E/N (Td) : ", (sys%E/meta(0)%Ni) * 1d21
    write(99,"(A,F8.3)")    "* E Field (V/cm) : ", sys%E*1d-2
    write(99,"(A,ES11.3)")  "* heating frequency (Hz) : ", sys%Freq / (2*pi)
    write(99,"(A,2ES11.3)") "* Power (W/cm3) | (W): ", sys%Powr*1d-6, sys%Powr * sys%volume
    write(99,"(A,F8.2)")    "* Sheath potential (eV) : ", Vg

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
    write(99,"((A,F7.2))" ) "* Gas Tp at the tube bound (°K): ", meta(0)%Tp*qok
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
    write(99,"(A)") "### Power balance : Electron Energy Saving"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Heat :  ", Diag(10)%EnProd * qe/(ne*Dt*clock%NumIter), &
         Diag(10)%EnProd *100.d0 / gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Penn :  ", Diag(5)%EnProd  * qe/(ne*Dt*clock%NumIter), &
         Diag(5)%EnProd *100.d0 / gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Asso :  ", Diag(6)%EnProd  * qe/(ne*Dt*clock%NumIter), &
         Diag(6)%EnProd *100.d0 / gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Dxct :  ", Diag(1)%EnProd  * qe/(ne*Dt*clock%NumIter), &
         Diag(1)%EnProd *100.d0 / gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Dxim :  ", Diag(12)%EnProd  * qe/(ne*Dt*clock%NumIter), &
         Diag(12)%EnProd *100.d0 / gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Pxcim:  ", Diag(4)%EnProd  * qe/(ne*Dt*clock%NumIter), &
         Diag(4)%EnProd *100.d0 / gainE, " %"

    write(99,"(A)") "### Power balance : Electron Energy Loss"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Elas :  ", Diag(11)%EnLoss * qe/(ne*Dt*clock%NumIter)&
         , Diag(11)%EnLoss*100.d0/gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Recb :  ", Diag(8)%EnLoss  * qe/(ne*Dt*clock%NumIter)&
         , Diag(8)%EnLoss*100.d0/gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Ionz :  ", Diag(2)%EnLoss  * qe/(ne*Dt*clock%NumIter)&
         , Diag(2)%EnLoss*100.d0/gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Exct :  ", Diag(1)%EnLoss  * qe/(ne*Dt*clock%NumIter)&
         , Diag(1)%EnLoss*100.d0/gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Diff :  ", Diag(9)%EnLoss  * qe/(ne*Dt*clock%NumIter)&
         , Diag(9)%EnLoss*100.d0/gainE, " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss Iexc :  ", Diag(13)%EnLoss  * qe/(ne*Dt*clock%NumIter)&
         , Diag(13)%EnLoss*100.d0/gainE, " %"
    
    write(99,"(/,A)") "-------------------------------------------------------"
    write(99,"(A)") "### Particle balance : Electron Saving"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain ioniz : ", Diag(2)%SumTx/(ne*clock%NumIter),  &
         Diag(2)%SumTx*100.d0 / abs(gainP), " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Assoc : ", Diag(6)%SumTx/(ne*clock%NumIter),  &
         Diag(6)%SumTx*100.d0 / abs(gainP), " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Penng : ", Diag(5)%SumTx/(ne*clock%NumIter),  &
         Diag(5)%SumTx*100.d0 / abs(gainP), " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Pexci : ", Diag(4)%SumTx/(ne*clock%NumIter),  &
         Diag(4)%SumTx*100.d0 / abs(gainP), " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Gain Ioexc : ", Diag(13)%SumTx/(ne*clock%NumIter),  &
         Diag(13)%SumTx*100.d0 / abs(gainP), " %"

    write(99,"(A)") "### Particle balance : Electron Loss"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss recmb : ", Diag(8)%SumTx/(ne*clock%NumIter),  &
         Diag(8)%SumTx*100.d0 / abs(gainP), " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss diffz : ", Diag(9)%SumTx/(ne*clock%NumIter),  &
         Diag(9)%SumTx*100.d0 / abs(gainP), " %"
    write(99,"(A,ES15.6,F8.2,A)") "* Loss recex : ", Diag(15)%SumTx/(ne*clock%NumIter),  &
         Diag(15)%SumTx*100.d0 / abs(gainP), " %"
    write(99,"(/,A)") "-------------------------------------------------------"
    write(99,"(A,2ES15.4)")"* Gain elec total (Pwr | Prtcl) : ", (qe/(ne*Dt*clock%NumIter) ) * &
         gainE, gainP
    write(99,"(A,2ES15.4)")"* Loss elec total (Pwr | Prtcl) : ", (qe/(ne*Dt*clock%NumIter) ) * &
         (Diag(11)%EnLoss+Diag(8)%EnLoss+Diag(2)%EnLoss+Diag(1)%EnLoss+Diag(9)%EnLoss), &
         (Diag(8)%SumTx+Diag(9)%SumTx+Diag(15)%SumTx)
    write(99,"(A)") " "
    DO i = 1, NumMeta                                                       !
       write(99,"(I3,A,F10.4,ES15.4,F8.2,A)") i, meta(i)%Name, meta(i)%En, &
            meta(i)%Ni*1d-06, meta(i)%Ni *100.d0/SumMeta, " %"
    END DO                                                                  !
    write(99,"(A)") " "
    write(99,"(A)")"![Zozor](http://uploads.siteduzero.com/files/420001_421000/420263.png)"
    Close(99)
  END SUBROUTINE OutPutMD
  !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!

  SUBROUTINE CHECK_AND_WRITE(Clock, sys, meta, elec, ion, F, diag, iter, MxDt)
    !**** INTENT
    INTEGER      , INTENT(IN)    :: iter
    REAL(DOUBLE) , INTENT(INOUT) :: MxDt
    TYPE(SysVar) , INTENT(IN)    :: sys
    TYPE(TIME)   , INTENT(INOUT) :: Clock
    TYPE(Species), INTENT(INOUT) :: elec
    TYPE(Species), DIMENSION(NumIon)   , INTENT(INOUT) :: ion
    TYPE(Species), DIMENSION(0:NumMeta), INTENT(INOUT) :: meta
    TYPE(Diagnos), DIMENSION(:)        , INTENT(INOUT) :: diag
    REAL(DOUBLE) , DIMENSION(:)        , INTENT(INOUT) :: F
    !**** LOCAL
    INTEGER :: i, nx, Mnul, Switch, mdlus
    REAL(DOUBLE) :: RateSum = 0.d0
    CHARACTER(LEN=250)::fileName
    nx = sys%nx ; Switch = 0 ; mdlus = 500

    !**** CHECK PART *********************************
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !**** CHECK PART *********************************

    !**** Check EEDF Positivty
    DO i = 1, nx
       IF (F(i).LT. 0.d0) THEN
          F(i) = 0.d0 ; Switch = 10
       END IF
    END DO

    !**** UpDate Electron (Density + Temperature)
    elec%Ni = 0.d0 ; elec%Tp = 0.d0
    DO i = 1, sys%nx 
       elec%Ni = elec%Ni + ( F(i) * dsqrt(U(i)) * sys%Dx )
       elec%Tp = elec%Tp + ( F(i) * dsqrt(U(i)**3) * sys%Dx )
    END DO
    elec%Tp = elec%Tp * 0.6667d0 / elec%Ni

    !**** Update densities (Ion + Excited)
    do i = 1, NumMeta
       !**** Update excited densities
       meta(i)%Ni = meta(i)%Ni + meta(i)%Updens
       
       IF (meta(i)%Ni < 0.d0) THEN
          meta(i)%Ni = 0.d0 ; Mnul = Mnul + 1
       END IF

       IF (i .LE. NumIon) THEN
          !**** Update ions (+ He2*) densities
          ion(i)%Ni = ion(i)%Ni + ion(i)%Updens
          IF (ion(i)%Ni < 0.d0) THEN
             ion(i)%Ni = 0.d0 ; Mnul = Mnul + 100
             IF (i == 1) ion(2)%Ni = elec%Ni
             IF (i == 2) ion(1)%Ni = elec%Ni
          END IF
       END IF
    END do

    !**** Report 
    !IF (Switch==10) write(*,"(2A,F8.3)") tabul, "Carefull, EEDF has negative(s) values ! time= ",&
    !     (Clock%SumDt*1e6)

    !**** Update Time-step
    IF ( 1.d0/MaxR.GE.1d-12 .and. iter.GT.1 ) THEN
       Clock%Dt = 1.0d0 / (MaxR*3.d0)
       IF (Clock%Dt .GT. MxDt) Clock%Dt = MxDt ! Maximum Time-Step allowed
    END IF
    !**** Check if there's NaN propagation ... probably due to large Dt (change MaxDt).
    IF (ISnan(MaxR) .or. isNaN(elec%Ni) ) THEN 
       print*,""; print*," NaN ! Pobleme in time-step??", elec%Ni, MaxR ; Stop 
    END IF
    MaxR = 0.d0
    !*************************************

    !**** WRITE PART *********************************
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !**** WRITE PART *********************************

    IF ( mod(iter,mdlus).EQ.0 ) THEN
       RateSum = -(diag(3)%OutM2 + diag(15)%InM1)*meta(3)%Ni + (diag(10)%OutM1)*meta(1)%Ni

       !**** WRITE Frequently IN TERMINAL **************!
!       write(*,"(2A,F8.3,A,F5.1,A,I7,A,ES9.3,A,F5.1,A,F5.1,A,ES10.2,A)",advance="no") &
!            tabul, "RunTime : ", (Clock%SumDt*1e6), " μs | ", Clock%SumDt*100.d0/Clock%SimuTime,&
!            "% [Nloop = ", iter, " | Dt = ", Clock%Dt, " | Pwr(%): ", (sys%Pcent*100./sys%IPowr),&
!            "] Sheath: ", Vg, " Emoy(V/m) ", sys%Emoy/iter, " \r"!
       write(*,"(A,F8.3,A,F5.1,A,ES9.3,A,F5.1,A,ES10.2,A)",advance="no") &
            tabul, Clock%SumDt*1e6, " μs | ", Clock%SumDt*100.d0/Clock%SimuTime,&
            "% [Dt = ", Clock%Dt, " | Pwr(%): ", (sys%Pcent*100./sys%IPowr)," n0/tin=", RateSum," \r"!

       !**** WRITE IN EVOL.DAT *************************!
       IF (Clock%Rstart.EQ.0 .and. iter.EQ.mdlus) THEN
          OPEN(UNIT=99,File=TRIM(ADJUSTL(DirFile))//"evol.dat",ACTION="WRITE",STATUS="UNKNOWN")
       ELSE 
          OPEN(UNIT=99,File=TRIM(ADJUSTL(DirFile))//"evol.dat",ACCESS="APPEND",&
               ACTION="WRITE",STATUS="UNKNOWN")
       END IF
       SELECT CASE (NumIon) 
       CASE (3)
          write(99,"(50ES15.4)") (Clock%SumDt*1e6), elec%Tp, meta(0)%Tp*qok,&
               OneD%Tg(OneD%bnd),elec%Ni*1d-06,ion(1)%Ni*1d-06, ion(2)%Ni*1d-06,&
               ion(NumIon)%Ni*1d-06, (meta(i)%Ni*1d-06,i=1,NumMeta)
       CASE DEFAULT
          write(99,"(48ES15.4)") (Clock%SumDt*1e6), elec%Tp, meta(0)%Tp*qok, elec%Ni*1d-06,&
               ion(1)%Ni*1d-06, ion(2)%Ni*1d-06, (meta(i)%Ni*1d-06,i=1,NumMeta)
       END SELECT
       CLOSE(99)
    ELSE IF ( mod(iter,Clock%NumIter/10) == 0 ) then

       !**** WRITE RARELY IN TERMINAL ******************!
       write(*,"(A,F7.2,A,3ES10.2,A,ES10.2,3(A,F7.1),2(A,2ES9.2))") &
            tabul//"**", (Clock%SumDt*1e6), " μs", meta(1)%Ni*1d-06, meta(3)%Ni*1d-06, elec%Ni," | E/N(Td)",&
            (sys%E/meta(0)%Ni)/1d-21, " | Tg(K)",meta(0)%Tp*qok," | Tg(bnd)",OneD%Tg(OneD%bnd), &
            "\n\t\t Pg(Torr)", meta(0)%Prs, " | M1 In/Out", diag(10)%InM1,diag(10)%OutM1, " | M2 In/Out", &
            diag(10)%InM2,diag(10)%OutM2

       !**** WRITE IN DENSITY.DAT (cm^{-3}) ************!
       OPEN(UNIT=98,File=TRIM(ADJUSTL(DirFile))//"density.dat",ACTION="WRITE",STATUS="UNKNOWN")
       DO i = 1, NumMeta+3
          IF (i.LE.NumMeta) write(98,"(I3,A,2F10.4,ES12.4)") i, meta(i)%Name, meta(i)%En, meta(i)%deg, meta(i)%Ni*1d-06
          IF (i.GT.NumMeta) write(98,"(I3,A,2F10.4,ES12.4)") i-NumMeta, ion(i-NumMeta)%Name, ion(i-NumMeta)%En, &
               ion(i-NumMeta)%deg, ion(i-NumMeta)%Ni*1d-06 
       END DO
       CLOSE(98)
       !************************************************!
    END IF

    IF ( Clock%SumDt.GE.Res ) THEN

       !**** WRITE RESTART FILES ***********************!
       CALL Rstart_SaveFiles (sys, Clock, ion, elec, meta, F)

       !**** WRITE EEDF ********************************!
       write(fileName,"('F_evol_',I3.3,'.dat')") int(Res/Clock%TRstart)
       OPEN(UNIT=90,File=TRIM(ADJUSTL(DirFile))//TRIM(ADJUSTL(fileName)),ACTION="WRITE",STATUS="UNKNOWN")
       DO i = 1, nx
          write(90,"(3ES15.6E3)") real(i)*sys%Dx, F(i), Clock%SumDt*1d06
       END DO 
       CLOSE(90)
       CALL Write_Out1D( OneD%Tg, "Tg.dat")
       Res = Res + Clock%TRstart 
    END IF

  END SUBROUTINE CHECK_AND_WRITE


END MODULE MOD_EVOL
