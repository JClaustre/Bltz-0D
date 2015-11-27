!  f90_kind.f90
!
!  Module defining the KIND numbers for the NAGWare f90 Compiler.
!
!  Copyright 1991,1993 Numerical Algorithms Group Ltd., Oxford, U.K.
!
!  Malcolm Cohen, Robert Iles, June 1991
!
      module f90_kind
!
      intrinsic kind,selected_int_kind,selected_real_kind  ! we use
      ! these,
      private kind,selected_int_kind,selected_real_kind    ! but do
      ! not force
                                                         ! them on the user.
!
! Indicator that the KIND= is not available for this compiler/host
      integer, parameter :: not_available = -1  
!
! Real and Complex numbers
!   Single precision
      integer, parameter :: single  = kind(0.0)
!   Double precision
      integer, parameter :: double  = kind(0.0d0)
!   Quadruple precision
      integer, parameter :: quad    = selected_real_kind(p=30)
!
! Integers numbers
!   Single byte integer
      integer, parameter :: int8    = selected_int_kind(2)
!   Two byte integer   
      integer, parameter :: int16   = selected_int_kind(4)
!   Four byte integer
      integer, parameter :: int32   = selected_int_kind(9)
!   Eight byte integer
      integer, parameter :: int64   = selected_int_kind(18)
!
! Logical values
!   Single byte logical
      integer, parameter :: byte    = 1
!   Two byte logical
      integer, parameter :: twobyte = 2
!   Four byte logical
      integer, parameter :: word    = kind(.TRUE.)
!
! Character type
!   Normal single byte character (ASCII sequence)
      integer, parameter :: ascii   = kind('x')
      end module f90_kind
