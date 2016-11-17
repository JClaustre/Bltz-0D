!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Author: Jonathan Claustre
! Date  : 15/05/2015
! Objctv: (Main) Inelastic collisions for global model
! Complt: see makefile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MAIN

  USE F90_KIND
  USE MOD_PARAM
  USE MOD_EVOL
  USE MOD_READ

  IMPLICIT NONE
  INTEGER :: t1, t2, clock_rate

  CALL welcome()
  !**** INIT PARAM ****!
  CALL Init(sys, Clock, ion, elec, meta, lasr)

  CALL SYSTEM_CLOCK (t1, clock_rate)
  CALL Evolution ()
  CALL SYSTEM_CLOCK (t2, clock_rate)
  CALL PrinTime (t1,t2,clock_rate)
  CALL OutPutMD (sys, meta, ion, elec, diag, consv)

  CALL DelocArray()
  write(*,"(2A)") tabul, "***** Goodbye ! *****"


CONTAINS

  SUBROUTINE welcome()
    

    write(*,"(A)")"      .-.     .-.     .-.     .-.     .-.     .-.     .-.    "
    write(*,"(A)")"     .'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `. "
    write(*,"(A)")"                _         "
    write(*,"(A)")"               [ ] Welcome into my program, "
    write(*,"(A)")"              (* *) /   and have fun in Plasma-Land!         "
    write(*,"(A)")"               |>|                  -C0PO-                   "
    write(*,"(A)")"            __/===\__     "
    write(*,"(A)")"           //| o=o |\\    "
    write(*,"(A)")"         <]  | o=o |  [>  "
    write(*,"(A)")"             \=====/      "
    write(*,"(A)")"            / / | \ \     "
    write(*,"(A)")"           <_________>    "
    write(*,"(A)")"      .-.     .-.     .-.     .-.     .-.     .-.     .-.    "
    write(*,"(A)")"     .'   `._.'   `._.'   `._.'   `._.'   `._.'   `._.'   `. "
    write(*,"(A)")""
    
  END SUBROUTINE welcome

END PROGRAM MAIN

