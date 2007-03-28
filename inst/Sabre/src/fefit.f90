      subroutine fefit( beta, cov, y, x, ncov, nsub, nmes, it, nest, nilev, &
                        con, tol, arith, n1sub, ifail, maxcol, mnames, ndum, &
                        gamma, gammse )
!
      use hsl_ma57_double
      use memory
      implicit none
      include 'limitsf90.h'
      logical ifail
      integer ncov, nsub(2), nmes, nest, nilev, n1sub, maxcol, ndum
      double precision beta( MAXPAR*(MAXPAR+5) ), cov( MAXPAR*(MAXPAR+1)/2 ), &
                       y(4*nmes), x(nmes,maxcol), con, tol, gamma(MAXY), &
                       gammse(MAXY)
      integer it( 2, nmes )
      character arith
      character*12 mnames( MAXVAR )

      include 'mpif.h'
      common /mpi_info/ num_processors, this_processor
      integer num_processors, this_processor
      double precision :: dummy, val, &
                          predicted, sigma, error, strip_size
      integer :: row, col, entry, n, ne, i, ierror, &
                 first, last, ncols, icase, &
                 npasses, processor, ipass, tag, status(MPI_STATUS_SIZE), &
                 strip_number, nvars, id, iobs, alloc_stat
      integer, allocatable :: firstcol(:), lastcol(:)
      integer, parameter :: MASTER = 0, MAXMEM = 10000000
      logical :: last_strip
      character*80 outbuf
      type(zd01_type) matrix
      TYPE(MA57_CONTROL) CONTROL
      TYPE(MA57_AINFO) AINFO
      TYPE(MA57_FINFO) FINFO
      TYPE(MA57_SINFO) SINFO
      TYPE(MA57_FACTORS) FACTORS
      double precision, allocatable :: B(:,:), XSOL(:,:), diag(:)
!
      ifail = .false.
!
      open( unit=9, file='fefit.mat' )
      n = nsub(1) - n1sub + nilev
      ne = nilev*(nilev+1)/2 + ( nsub(1) - n1sub )*( nilev + 1 )
      MATRIX%N = N
      MATRIX%NE = NE
!
!    Allocate arrays of appropriate sizes
!
      ALLOCATE(MATRIX%VAL(NE), MATRIX%ROW(NE), MATRIX%COL(NE), stat=alloc_stat )
      if( alloc_stat /= 0 ) then
         ifail = .true.
         call wrtlin( '    *** ERROR *** '// &
                      'Failed to allocate space for matrix in FEFIT' )
         return
      end if
      ALLOCATE( B(N,1), XSOL(N,1), stat=alloc_stat )
      if( alloc_stat /= 0 ) then
         ifail = .true.
         call wrtlin( '    *** ERROR *** '// &
                      'Failed to allocate space for solution vector in FEFIT' )
         return
      end if
      call update_memory( 16*( ne + n ) )
!
!    Work out matrix for the model
!
!    Work out x'x and store in matrix
!
      entry = 0
      do row = 1, nilev
         do col = row, nilev
            entry = entry + 1
            matrix%row(entry) = row
            matrix%col(entry) = col
            val = 0d0
            i = 0
            do icase = 1, nsub(1)
               if( it(1,icase) > 1 ) then
                  do iobs = 1, it(1,icase)
                     i = i + 1               
                     val = val + x(i,row)*x(i,col)
                  end do
               else
                  i = i + 1
               end if
            end do
            matrix%val(entry) = val
            write( 9, * ) entry, row, col, val
         end do
      end do
!
!    Work out x'd and store in matrix
!
      do row = 1, nilev
         i = 0
         col = nilev
         do icase = 1, nsub(1)
            if( it(1,icase) > 1 ) then
               col = col + 1
               entry = entry + 1
               matrix%row(entry) = row
               matrix%col(entry) = col
               val = 0d0
               do iobs = 1, it(1,icase)
                  i = i + 1
                  val = val + x(i,row)
               end do
               matrix%val(entry) = val
               write( 9, * ) entry, row, col, val
            else
               i = i + 1
            end if
         end do
      end do
!
!    Add d'd
!
      row = nilev
      do icase = 1, nsub(1)
!       Don't include cases with just 1 observation
         if( it(1,icase) > 1 ) then
            row = row + 1
            entry = entry + 1
            matrix%row(entry) = row
            matrix%col(entry) = row
            matrix%val(entry) = it( 1, icase )
            write( 9, * ) entry, row, row, matrix%val(entry)
         end if
      end do
!
!    Calculate x'y and put in b
!
      do row = 1, nilev
         val = 0d0
         i = 0
         do icase = 1, nsub(1)
!          Don't include cases with only 1 observation
            if( it(1,icase) > 1 ) then
               do iobs = 1, it(1,icase)
                  i = i + 1
                  val = val + x(i,row)*y(i)
               end do
            else
               i = i + 1
            end if
         end do
         b(row,1) = val
      end do
!
!    Calculate d'y and add to b
!
      i = 0
      entry = nilev
      do row = nilev+1, nilev+nsub(1)
         if( it(1,row-nilev) > 1 ) then
            val = 0d0
            do iobs = 1, it( 1, row-nilev )
               i = i + 1
               val = val + y(i)
            end do
            entry = entry + 1
            b(entry,1) = val
         else
            i = i + 1
         end if
      end do


! Initialize the structures
      CALL MA57_INITIALIZE(FACTORS,CONTROL)

! Analyse
      CALL MA57_ANALYSE(MATRIX,FACTORS,CONTROL,AINFO)
      IF(AINFO%FLAG<0) THEN
         ifail = .true.
         call wrtlin( '    *** HSL ERROR ***' )
         WRITE(outbuf,'(A,I2)') &
            ' Failure of MA57_ANALYSE with AINFO%FLAG=', AINFO%FLAG
         call wrtlin( outbuf )
         return
      END IF

! Factorize
      CALL MA57_FACTORIZE(MATRIX,FACTORS,CONTROL,FINFO)
      IF(FINFO%FLAG<0) THEN
         ifail = .true.
         call wrtlin( '    *** HSL ERROR ***' )
         WRITE(outbuf,'(A,I2)') &
            ' Failure of MA57_FACTORIZE with FINFO%FLAG=', FINFO%FLAG
         call wrtlin( outbuf )
         return
      END IF

! Solve without refinement
     XSOL = B
      CALL MA57_SOLVE(MATRIX,FACTORS,XSOL,CONTROL,SINFO)

! Perform one refinement
      CALL MA57_SOLVE(MATRIX,FACTORS,XSOL,CONTROL,SINFO,B)
!
!    Copy the estimates of the variates into beta and estimates of the
!    dummy variables into gamma
!
      do i = 1, nilev
         beta(i) = xsol(i,1)
      end do
      do i = nilev+1, nilev + nsub(1) - n1sub
         gamma(i-nilev) = xsol(i,1)
      end do
!
!    Calculate sigma
!
      sigma = 0d0
      i = 0
      row = 0
      do icase = 1, nsub(1)
         if( it( 1, icase ) > 1 ) then
            row = row + 1
            do iobs = 1, it( 1, icase )
               i = i + 1
!
!             Calculate the predicted response for this observation
!
               predicted = 0d0
               do col = 1, nilev
                  predicted = predicted + xsol(col,1) * x( i, col )
               end do
!
!             Add in the fixed effect
!
               predicted = predicted + xsol( row+nilev, 1 )
               error = y(i) - predicted
               sigma = sigma + error**2
            end do
         else
            i = i + 1
         end if
      end do
!
!    Divisor is no. of observations - no. of variates - no. of cases - 1
!
!    Note that no. of observations is the total number of observations minus
!    the number of cases where there is just a single observation and no.
!    of cases is the total number minus those with just a single observation.
!
      sigma = sigma/ ( ( nmes-n1sub) - nilev - ( nsub(1)-n1sub ) - 1 )
!
!    Now compute the inverse of x'x
!
!    B becomes the identity matrix, shared across the processors
!
!    Divide the matrix between the processors
!
      strip_size = real(n) / num_processors
      
      npasses = 1
!
!    Further reduce the strip size if it requires too much memory for
!    one processor
!
      if( n * strip_size > MAXMEM ) then
         strip_size = max( 1d0, real( MAXMEM, kind(1d0) ) / n )
!
!       Calculate the number of passes (around all the processors)
!       required to process the whole of the inverse
!
         npasses = ceiling( real( ceiling( (real(n) / strip_size ) ) ) / &
                   num_processors )
      end if
      allocate( diag(n), stat=alloc_stat )
      if( alloc_stat /= 0 ) then
         ifail = .true.
         call wrtlin( '    *** ERROR *** '// &
                      'Failed to allocate space for diag in FEFIT' )
         return
      end if
      call update_memory( n*8 )
      allocate( firstcol( npasses * num_processors ), stat=alloc_stat )
      if( alloc_stat /= 0 ) then
         ifail = .true.
         call wrtlin( '    *** ERROR *** '// &
                      'Failed to allocate space for firstcol in FEFIT' )
         return
      end if
      allocate( lastcol( npasses * num_processors ), stat=alloc_stat )
      if( alloc_stat /= 0 ) then
         ifail = .true.
         call wrtlin( '    *** ERROR *** '// &
                      'Failed to allocate space for lastcol in FEFIT' )
         return
      end if
      call update_memory( npasses*num_processors*8 )
!
      last = 0
      strip_number = 0
      last_strip = .false.
      do ipass = 1, npasses
         do processor = 0, num_processors-1
            strip_number = strip_number + 1
!
!          Update the counters on each processor
!
            firstcol( strip_number ) = last+1
            lastcol( strip_number ) = nint( strip_size * &
             ( ( ipass-1 ) * num_processors + processor+1 ) )
            if( lastcol( strip_number ) >= n ) then
               lastcol( strip_number ) = n
               last_strip = .true.
            end if
            first = firstcol( strip_number )
            last = lastcol( strip_number )
            ncols = last - first + 1
!
!          See if this processor should solve the current strip
!
            if( processor == this_processor ) then
               call update_memory( -( size(b)*8 + size(x)*8 ) )
               deallocate( b, xsol )
               allocate( b(n,ncols), xsol(n,ncols), stat=alloc_stat )
               if( alloc_stat /= 0 ) then
                  ifail = .true.
                  call wrtlin( '    *** ERROR *** '// &
                               'Failed to allocate space for b & xsol in FEFIT')
                  return
               end if
               call update_memory( 16*n*ncols )
               do col = first, last
                  do row = 1, n
                     b( row, col-first+1 ) = 0d0
                  end do
                  b( col, col-first+1 ) = 1d0
               end do
!
!             Solve without refinement
!
               XSOL = B
               CALL MA57_SOLVE(MATRIX,FACTORS,XSOL,CONTROL,SINFO)
!
!             Perform one refinement
!
               CALL MA57_SOLVE(MATRIX,FACTORS,XSOL,CONTROL,SINFO,B)
!
!             Check by pre-multiplying original matrix if this is MASTER
!
!             Store the diagonal elements of X
!
               do col = 1, ncols
                  diag( first+col-1 ) = xsol( first+col-1, col )
               end do
!
!             Slaves send their part of the diagonal for MASTER to collect
!
               if( this_processor /= MASTER ) then
                  call mpi_send( diag(first), ncols, MPI_DOUBLE_PRECISION, &
                   MASTER, 0, MPI_COMM_WORLD, ierror )
               end if
!
!             End of parallel processing part
!
            end if
!
!          MASTER collects (apart from itself)
!
            if( this_processor == MASTER .and. processor /= MASTER ) then
               call mpi_recv( diag(first), ncols, MPI_DOUBLE_PRECISION, &
                processor, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror )
               if( ierror /= 0 ) then
                  write( outbuf, '(a,i7,a)' ) &
                   '    *** ERROR *** Error in receiving diag(', first, &
                   ') in FEFIT'
                  call wrtlin( outbuf )
               end if
            end if
!
!          Finished all the strips
!
            if( last_strip ) go to 10
         end do
      end do
!
  10  continue
!
!    Copy variances of variates into COV and standard errors of dummy variables
!    into gammse
!
      do i = 1, nilev
         cov(i*(i+1)/2) = sigma*diag(i)
      end do
      do i = nilev+1, nilev + nsub(1) - n1sub
         gammse(i-nilev) = sqrt(sigma*diag(i))
      end do
!
      return
      end
