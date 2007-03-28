      subroutine dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                         beta, c, ldc )
!     .. Scalar Arguments ..
      character*1        transa, transb
      integer            m, n, k, lda, ldb, ldc
      double precision   alpha, beta
!     .. Array Arguments ..
      double precision   a( lda, * ), b( ldb, * ), c( ldc, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      logical            lsame
      external           lsame
!     .. External Subroutines ..
      external           xerbla
!     .. Intrinsic Functions ..
      intrinsic          max
!     .. Local Scalars ..
      logical            nota, notb
      integer            i, info, j, l, ncola, nrowa, nrowb
      double precision   temp
!     .. Parameters ..
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      nota  = lsame( transa, 'N' )
      notb  = lsame( transb, 'N' )
      if( nota )then
         nrowa = m
         ncola = k
      else
         nrowa = k
         ncola = m
      end if
      if( notb )then
         nrowb = k
      else
         nrowb = n
      end if
!
!     Test the input parameters.
!
      info = 0
      if(      ( .not.nota                 ).and. &
               ( .not.lsame( transa, 'C' ) ).and. &
               ( .not.lsame( transa, 'T' ) )      )then
         info = 1
      else if( ( .not.notb                 ).and. &
               ( .not.lsame( transb, 'C' ) ).and. &
               ( .not.lsame( transb, 'T' ) )      )then
         info = 2
      else if( m  .lt.0               )then
         info = 3
      else if( n  .lt.0               )then
         info = 4
      else if( k  .lt.0               )then
         info = 5
      else if( lda.lt.max( 1, nrowa ) )then
         info = 8
      else if( ldb.lt.max( 1, nrowb ) )then
         info = 10
      else if( ldc.lt.max( 1, m     ) )then
         info = 13
      end if
      if( info.ne.0 )then
         call xerbla( 'DGEMM ', info )
         return
      end if
!
!     Quick return if possible.
!
      if( ( m.eq.0 ).or.( n.eq.0 ).or. &
          ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) ) &
         return
!
!     And if  alpha.eq.zero.
!
      if( alpha.eq.zero )then
         if( beta.eq.zero )then
            do 20, j = 1, n
               do 10, i = 1, m
                  c( i, j ) = zero
   10          continue
   20       continue
         else
            do 40, j = 1, n
               do 30, i = 1, m
                  c( i, j ) = beta*c( i, j )
   30          continue
   40       continue
         end if
         return
      end if
!
!     Start the operations.
!
      if( notb )then
         if( nota )then
!
!           Form  C := alpha*A*B + beta*C.
!
            do 90, j = 1, n
               if( beta.eq.zero )then
                  do 50, i = 1, m
                     c( i, j ) = zero
   50             continue
               else if( beta.ne.one )then
                  do 60, i = 1, m
                     c( i, j ) = beta*c( i, j )
   60             continue
               end if
               do 80, l = 1, k
                  if( b( l, j ).ne.zero )then
                     temp = alpha*b( l, j )
                     do 70, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
   70                continue
                  end if
   80          continue
   90       continue
         else
!
!           Form  C := alpha*A'*B + beta*C
!
            do 120, j = 1, n
               do 110, i = 1, m
                  temp = zero
                  do 100, l = 1, k
                     temp = temp + a( l, i )*b( l, j )
  100             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  110          continue
  120       continue
         end if
      else
         if( nota )then
!
!           Form  C := alpha*A*B' + beta*C
!
            do 170, j = 1, n
               if( beta.eq.zero )then
                  do 130, i = 1, m
                     c( i, j ) = zero
  130             continue
               else if( beta.ne.one )then
                  do 140, i = 1, m
                     c( i, j ) = beta*c( i, j )
  140             continue
               end if
               do 160, l = 1, k
                  if( b( j, l ).ne.zero )then
                     temp = alpha*b( j, l )
                     do 150, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  150                continue
                  end if
  160          continue
  170       continue
         else
!
!           Form  C := alpha*A'*B' + beta*C
!
            do 200, j = 1, n
               do 190, i = 1, m
                  temp = zero
                  do 180, l = 1, k
                     temp = temp + a( l, i )*b( j, l )
  180             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  190          continue
  200       continue
         end if
      end if
!
      return
!
!     End of DGEMM .
!
      end
      subroutine dtpsv ( uplo, trans, diag, n, ap, x, incx )
!
!***********************************************************************
!
!     .. Scalar Arguments ..
      integer            incx, n
      character*1        diag, trans, uplo
!     .. Array Arguments ..
      double precision   ap( * ), x( * )
!     ..
!
!  Purpose
!  =======
!
!  DTPSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular matrix, supplied in packed form.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  AP     - DOUBLE PRECISION array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular matrix packed sequentially,
!           column by column, so that AP( 1 ) contains a( 1, 1 ),
!           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
!           respectively, and so on.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular matrix packed sequentially,
!           column by column, so that AP( 1 ) contains a( 1, 1 ),
!           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
!           respectively, and so on.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
!     .. Local Scalars ..
      double precision   temp
      integer            i, info, ix, j, jx, k, kk, kx
      logical            nounit
!     .. External Functions ..
      logical            lsame
      external           lsame
!     .. External Subroutines ..
      external           xerbla
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0
      if     ( .not.lsame( uplo , 'U' ).and. &
               .not.lsame( uplo , 'L' )      )then
         info = 1
      else if( .not.lsame( trans, 'N' ).and. &
               .not.lsame( trans, 'T' ).and. &
               .not.lsame( trans, 'C' )      )then
         info = 2
      else if( .not.lsame( diag , 'U' ).and. &
               .not.lsame( diag , 'N' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( incx.eq.0 )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'DTPSV ', info )
         return
      end if
!
!     Quick return if possible.
!
      if( n.eq.0 ) &
         return
!
      nounit = lsame( diag, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
!
!     Start the operations. In this version the elements of AP are
!     accessed sequentially with one pass through AP.
!
      if( lsame( trans, 'N' ) )then
!
!        Form  x := inv( A )*x.
!
         if( lsame( uplo, 'U' ) )then
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit ) &
                        x( j ) = x( j )/ap( kk )
                     temp = x( j )
                     k    = kk     - 1
                     do 10, i = j - 1, 1, -1
                        x( i ) = x( i ) - temp*ap( k )
                        k      = k      - 1
   10                continue
                  end if
                  kk = kk - j
   20          continue
            else
               jx = kx + ( n - 1 )*incx
               do 40, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     if( nounit ) &
                        x( jx ) = x( jx )/ap( kk )
                     temp = x( jx )
                     ix   = jx
                     do 30, k = kk - 1, kk - j + 1, -1
                        ix      = ix      - incx
                        x( ix ) = x( ix ) - temp*ap( k )
   30                continue
                  end if
                  jx = jx - incx
                  kk = kk - j
   40          continue
            end if
         else
            kk = 1
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     if( nounit ) &
                        x( j ) = x( j )/ap( kk )
                     temp = x( j )
                     k    = kk     + 1
                     do 50, i = j + 1, n
                        x( i ) = x( i ) - temp*ap( k )
                        k      = k      + 1
   50                continue
                  end if
                  kk = kk + ( n - j + 1 )
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if( x( jx ).ne.zero )then
                     if( nounit ) &
                        x( jx ) = x( jx )/ap( kk )
                     temp = x( jx )
                     ix   = jx
                     do 70, k = kk + 1, kk + n - j
                        ix      = ix      + incx
                        x( ix ) = x( ix ) - temp*ap( k )
   70                continue
                  end if
                  jx = jx + incx
                  kk = kk + ( n - j + 1 )
   80          continue
            end if
         end if
      else
!
!        Form  x := inv( A' )*x.
!
         if( lsame( uplo, 'U' ) )then
            kk = 1
            if( incx.eq.1 )then
               do 100, j = 1, n
                  temp = x( j )
                  k    = kk
                  do 90, i = 1, j - 1
                     temp = temp - ap( k )*x( i )
                     k    = k    + 1
   90             continue
                  if( nounit ) &
                     temp = temp/ap( kk + j - 1 )
                  x( j ) = temp
                  kk     = kk   + j
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  do 110, k = kk, kk + j - 2
                     temp = temp - ap( k )*x( ix )
                     ix   = ix   + incx
  110             continue
                  if( nounit ) &
                     temp = temp/ap( kk + j - 1 )
                  x( jx ) = temp
                  jx      = jx   + incx
                  kk      = kk   + j
  120          continue
            end if
         else
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 140, j = n, 1, -1
                  temp = x( j )
                  k = kk
                  do 130, i = n, j + 1, -1
                     temp = temp - ap( k )*x( i )
                     k    = k    - 1
  130             continue
                  if( nounit ) &
                     temp = temp/ap( kk - n + j )
                  x( j ) = temp
                  kk     = kk   - ( n - j + 1 )
  140          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 160, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  do 150, k = kk, kk - ( n - ( j + 1 ) ), -1
                     temp = temp - ap( k )*x( ix )
                     ix   = ix   - incx
  150             continue
                  if( nounit ) &
                     temp = temp/ap( kk - n + j )
                  x( jx ) = temp
                  jx      = jx   - incx
                  kk      = kk   - (n - j + 1 )
  160          continue
            end if
         end if
      end if
!
      return
!
!     End of DTPSV .
!
      end
!       Toolpack tool decs employed.
!       Arg dimension set to *.
!
      integer function idamax(n,dx,incx)
!
!     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
!     .. Scalar Arguments ..
      integer incx,n
!     ..
!     .. Array Arguments ..
      double precision dx(*)
!     ..
!     .. Local Scalars ..
      double precision dmax
      integer i,ix
!     ..
!     .. Intrinsic Functions ..
      intrinsic dabs
!     ..
!     .. Executable Statements ..
!
      idamax = 0
      if (n.lt.1) return
      idamax = 1
      if (n.eq.1) return
      if (incx.eq.1) go to 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
        if (dabs(dx(ix)).le.dmax) go to 5
        idamax = i
        dmax = dabs(dx(ix))
    5   ix = ix + incx
   10 continue
      return
!
!        CODE FOR INCREMENT EQUAL TO 1
!
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
        if (dabs(dx(i)).le.dmax) go to 30
        idamax = i
        dmax = dabs(dx(i))
   30 continue
      return

      end
!       A BLAS routine modified for use with HSL
!       Toolpack tool decs employed.
!       Contained comment lines which caused Decs to fail.
!
      subroutine xerbla(srname,info)
!     .. Scalar Arguments ..
      integer info
      character srname*6
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the Level 2 BLAS routines.
!
!  It is called by the Level 2 BLAS routines if an input parameter is
!  invalid.
!
!  Installers should consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Parameters
!  ==========
!
!  SRNAME - CHARACTER*6.
!           On entry, SRNAME specifies the name of the routine which
!           called XERBLA.
!
!  INFO   - INTEGER.
!           On entry, INFO specifies the position of the invalid
!           parameter in the parameter-list of the calling routine.
!
!
!  Auxiliary routine for Level 2 Blas.
!
!  Written on 20-July-1986.
!
!     .. Executable Statements ..
!
      write (*,fmt=99999) srname,info
!
      stop
!
99999 format (' ** On entry to ',a6,' parameter number ',i2, &
             ' had an illegal value')
!
!     End of XERBLA.
!
      end
      logical function lsame ( ca, cb )
!     .. Scalar Arguments ..
      character*1            ca, cb
!     ..
!
!  Purpose
!  =======
!
!  LSAME  tests if CA is the same letter as CB regardless of case.
!
!  N.B. This version of the routine is only correct for ASCII code.
!       Installers must modify the routine for other character-codes.
!
!       For EBCDIC systems the constant IOFF must be changed to -64.
!       For CDC systems using 6-12 bit representations, the system-
!       specific code in comments must be activated.
!
!  Parameters
!  ==========
!
!  CA     - CHARACTER*1
!  CB     - CHARACTER*1
!           On entry, CA and CB specify characters to be compared.
!           Unchanged on exit.
!
!
!  Auxiliary routine for Level 2 Blas.
!
!  -- Written on 11-October-1988.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, Nag Central Office.
!
!     .. Parameters ..
      integer                ioff
      parameter            ( ioff=32 )
!     .. Intrinsic Functions ..
      intrinsic              ichar
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      lsame = ca .eq. cb
!
!     Now test for equivalence
!
      if ( .not.lsame ) then
         lsame = ichar(ca) - ioff .eq. ichar(cb)
      end if
      if ( .not.lsame ) then
         lsame = ichar(ca) .eq. ichar(cb) - ioff
      end if
!
      return
!
!  The following comments contain code for CDC systems using 6-12 bit
!  representations.
!
!     .. Parameters ..
!     INTEGER                ICIRFX
!     PARAMETER            ( ICIRFX=62 )
!     .. Scalar Arguments ..
!     CHARACTER*1            CB
!     .. Array Arguments ..
!     CHARACTER*1            CA(*)
!     .. Local Scalars ..
!     INTEGER                IVAL
!     .. Intrinsic Functions ..
!     INTRINSIC              ICHAR, CHAR
!     .. Executable Statements ..
!
!     See if the first character in string CA equals string CB.
!
!     LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX)
!
!     IF (LSAME) RETURN
!
!     The characters are not identical. Now check them for equivalence.
!     Look for the 'escape' character, circumflex, followed by the
!     letter.
!
!     IVAL = ICHAR(CA(2))
!     IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN
!        LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB
!     END IF
!
!     RETURN
!
!     End of LSAME.
!
      end
      subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx, &
!
!***********************************************************************
!
!     File of the DOUBLE PRECISION  Level-2 BLAS.
!     ===========================================
!
!     SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
!    $                   BETA, Y, INCY )
!
!     SUBROUTINE DGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
!    $                   BETA, Y, INCY )
!
!     SUBROUTINE DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
!    $                   BETA, Y, INCY )
!
!     SUBROUTINE DSBMV ( UPLO, N, K, ALPHA, A, LDA, X, INCX,
!    $                   BETA, Y, INCY )
!
!     SUBROUTINE DSPMV ( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY )
!
!     SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
!
!     SUBROUTINE DTBMV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
!
!     SUBROUTINE DTPMV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
!
!     SUBROUTINE DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
!
!     SUBROUTINE DTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
!
!     SUBROUTINE DTPSV ( UPLO, TRANS, DIAG, N, AP, X, INCX )
!
!     SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!
!     SUBROUTINE DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
!
!     SUBROUTINE DSPR  ( UPLO, N, ALPHA, X, INCX, AP )
!
!     SUBROUTINE DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!
!     SUBROUTINE DSPR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, AP )
!
!     See:
!
!        Dongarra J. J., Du Croz J. J., Hammarling S.  and Hanson R. J..
!        An  extended  set of Fortran  Basic Linear Algebra Subprograms.
!
!        Technical  Memoranda  Nos. 41 (revision 3) and 81,  Mathematics
!        and  Computer Science  Division,  Argonne  National Laboratory,
!        9700 South Cass Avenue, Argonne, Illinois 60439, US.
!
!        Or
!
!        NAG  Technical Reports TR3/87 and TR4/87,  Numerical Algorithms
!        Group  Ltd.,  NAG  Central  Office,  256  Banbury  Road, Oxford
!        OX2 7DE, UK,  and  Numerical Algorithms Group Inc.,  1101  31st
!        Street,  Suite 100,  Downers Grove,  Illinois 60515-1263,  USA.
!
!***********************************************************************
!
                         beta, y, incy )
!     .. Scalar Arguments ..
      double precision   alpha, beta
      integer            incx, incy, lda, m, n
      character*1        trans
!     .. Array Arguments ..
      double precision   a( lda, * ), x( * ), y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
!     .. Local Scalars ..
      double precision   temp
      integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
!     .. External Functions ..
      logical            lsame
      external           lsame
!     .. External Subroutines ..
      external           xerbla
!     .. Intrinsic Functions ..
      intrinsic          max
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0
      if     ( .not.lsame( trans, 'N' ).and. &
               .not.lsame( trans, 'T' ).and. &
               .not.lsame( trans, 'C' )      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( lda.lt.max( 1, m ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'DGEMV ', info )
         return
      end if
!
!     Quick return if possible.
!
      if( ( m.eq.0 ).or.( n.eq.0 ).or. &
          ( ( alpha.eq.zero ).and.( beta.eq.one ) ) ) &
         return
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      if( lsame( trans, 'N' ) )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero ) &
         return
      if( lsame( trans, 'N' ) )then
!
!        Form  y := alpha*A*x + y.
!
         jx = kx
         if( incy.eq.1 )then
            do 60, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  do 50, i = 1, m
                     y( i ) = y( i ) + temp*a( i, j )
   50             continue
               end if
               jx = jx + incx
   60       continue
         else
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy   = ky
                  do 70, i = 1, m
                     y( iy ) = y( iy ) + temp*a( i, j )
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
   80       continue
         end if
      else
!
!        Form  y := alpha*A'*x + y.
!
         jy = ky
         if( incx.eq.1 )then
            do 100, j = 1, n
               temp = zero
               do 90, i = 1, m
                  temp = temp + a( i, j )*x( i )
   90          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  100       continue
         else
            do 120, j = 1, n
               temp = zero
               ix   = kx
               do 110, i = 1, m
                  temp = temp + a( i, j )*x( ix )
                  ix   = ix   + incx
  110          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  120       continue
         end if
      end if
!
      return
!
!     End of DGEMV .
!
      end

