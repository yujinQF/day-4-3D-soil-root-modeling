!> \file mpiutils.f90
!! \brief Module for MPI.

!> Modul: mpiutils
!> \author Horst Hardelauf
!> \date 10.04.2017
Module MPIutils
  Use typedef
  Integer, Save :: myrank
  Logical, Save :: mpi=.false.
  
Contains
  
  !> Stop the program
  !> \param text ... Error text
  Subroutine stop_program(text)
    Implicit None
    Character(*), Intent(in) :: text
    Write(*,'(//,a,a,//)') ' RSWMS Runtime error: ',Trim(text)
    Stop
  End Subroutine stop_program

End Module MPIutils
