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


  !**** Change the terminal size (used for the "gnome-terminal"):
  call execute_command_line ('resize -s 40 122')
  !**** Ref of My paper! :)
  CALL welcome()
  
  !**** INIT PARAM ****!
  CALL Init(sys, Clock, ion, elec, meta, lasr)
  !**** Enter the main loop
  CALL SYSTEM_CLOCK (t1, clock_rate)
  CALL Evolution ()
  CALL SYSTEM_CLOCK (t2, clock_rate)
  !**** Write results & simulation time in files
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
    write(*,"(A)")"     Ref: J Claustre et al. (doi:10.1088/1361-6595/aa8a16) "
    write(*,"(A)")""
    
  END SUBROUTINE welcome

END PROGRAM MAIN

