!
!    Include file to make available the maximum problem dimensions:
!
!    MXX    : size of the data matrix
!    MAXY   : maximum number of observations
!    MAXVAR : maximum number of variables
!    MAXPAR : maximum number of parameters
!    BLOCKSIZE : maximum number of numbers to be passed in a single MPI
!                routine call (parallel version only)
!
      integer :: mxx, maxy, maxvar, maxpar
      integer, parameter :: blocksize = 50000000
      common /data_space/ mxx, maxy, maxvar, maxpar
