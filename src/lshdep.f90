!
!***********************************************************************
!
      subroutine lshdep(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                        score,hess,ncov,flag,ifail,offlag,offpos,maxcol, &
                        link,iter,cflag,step,scnorm,scores,family,n1lev, &
                        risk,ind,xil,sign,iquad,aquad,mus,taus)
      use accurate_arithmetic
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character link(3),family(3)
      integer nmes,nest,ncov,maxcol,nsub(2),nm(3),it(2,nsub(1)),flag, &
              nilev,n1lev,iter,cflag,offpos(3),risk(nmes),ind,sign,iquad
      double precision y(nmes),x(nmes,maxcol),beta(nest),xll, &
             score(nest),hess(ncov),quad(3,2,256),scnorm,step, &
             scores(nmes,nest),xil,aquad(3,2,256),mus(3,nsub(1)), &
             taus(3,nsub(1))
      logical ifail,offlag(3)
!-----------------------------------------------------------------------
!     function : calculates deviance, score vector and hessian matrix at
!                 the current parameter estimates.
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character (len=320) :: outbuf
      character(len=5) :: str
      double precision theta(maxpar),qprob(256),e2,r,rr,e8(maxpar),b,pi, &
             e10(maxpar*(maxpar+3)/2),meansc(maxpar),e1(2),pil(nmes), &
             sumxb(nmes),e3(2),expme,deri(maxpar),e,phi,listar, &
             d2(maxpar*(maxpar+1)/2),d1(maxpar),qloc(256), &
             e7(maxpar),der1(maxpar),dpieta,e3s,rs,e1s,e1ss,rrs,rrss, &
             e8s(maxpar),qpe2,sigmae
      double precision score_total(maxpar),xll_total, &
             hesinc(maxpar*(maxpar+1)/2)
      integer t,i,irow,l,j,k,m,rowind(nsub(1)),yy,nq,nc,i1,i2
      logical sflag,anyfail,myturn,iflow(2),can
      double precision :: xll_first,score_first(nest),hess_first(ncov)
      integer :: outer,firstcase,lastcase
!-----------------------------------------------------------------------
!     mpi parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
      can = .true.
!
      if (link(1) == 'p' .or. link(1) == 'c') then
          can = .false.
      end if
!
      xll = 0
!
      irow = 1
      sflag = .false.
!
      do 18 k = 1,nest
          theta(k) = beta(k)
          score(k) = 0
   18 end do
!
      do 15 k = 1,ncov
          hess(k) = 0
   15 end do
!
      do 1200 i = 1,nsub(1)
          rowind(i) = irow
          irow = irow + it(1,i)
 1200 end do
!
      nq = 0
      nc = 0
!
      if (ind == 0) then
          i1 = 1
          i2 = nsub(1)
      else
          i1 = ind
          i2 = i1
      end if
!
!---- loop around cases
!
!     Split the 200 loop into the first case (outer=1) and all the rest
!     (outer=2)
!
      do outer = 1,2
!
      if (outer == 1) then
          firstcase = i1
          lastcase = i1
      else
          firstcase = i1+1
          lastcase = i2
      end if
!
      ifail = .false.
!     Values of xll, score & hess after the first case are saved and the
!     variables are reset to zero for the parallel part (outer=2). The
!     initial values are then added after the final case. This is to
!     avoid the values after the first iteration being multiplied by the
!     number of threads in the REDUCTION if OPENMP is being used
!
      if (outer == 2) then
          xll_first = xll
          xll = 0d0
!
          do k = 1,nest
              score_first(k) = score(k)
              score(k) = 0d0
          end do
!
          do k = 1,ncov
              hess_first(k) = hess(k)
              hess(k) = 0d0
          end do
!
      end if
!
!$OMP PARALLEL DO IF(outer == 2), DEFAULT(FIRSTPRIVATE), &
!$OMP SHARED(firstcase,lastcase,i1,i2,ind,this_processor), &
!$OMP SHARED(num_processors,nest,ncov,nm,iquad,quad,mus,taus,can,nilev), &
!$OMP SHARED(rowind,it,risk,family,theta,offlag,offpos,x,link,y,flag), &
!$OMP SHARED(scores,n1lev,z1,z2,sflag,sign,nmes), &
!$OMP SHARED(qprob,qloc,deri,aquad), &
!$OMP REDUCTION(+: xll,score,hess), REDUCTION(.or.: ifail)
!
      do 200 i = firstcase,lastcase
!
!         Any ifail must jump to the end of the loop since jumping out
!         of a parallel section is not allowed in Openmp
!
!         Cycle if a previous iteration has failed
!
          if (ifail) then
              go to 2000
          end if
!
          if (i == i1 .or. ind /= 0) then
              myturn = .true.
          else
              myturn = this_processor == mod(i-2,num_processors)
          end if
!
          if (.not. myturn .and. i /= i2) then
              go to 200
          end if
!
          ifail = .false.
          iflow(1) = .false.
          iflow(2) = .false.
          listar = 0
!
          do 20 k = 1,nest
              d1(k) = 0
   20     end do
!
          do 25 k = 1,ncov
              d2(k) = 0
              hesinc(k) = 0
   25     end do
!
          do 150 j = 1,nm(1)
!
              if (iquad == 0 .and. i == i1) then
                  qloc(j) = quad(1,1,j)
                  qprob(j) = quad(1,2,j)
              else if (iquad /= 0) then
!
                  if (ind == 0  .and. outer == 1) then
                      aquad(1,1,j) = mus(1,i) + taus(1,i)*quad(1,1,j)
                      aquad(1,2,j) = taus(1,i)*quad(1,2,j)* &
                      exp((quad(1,1,j)**2 - aquad(1,1,j)**2)/2)
                  end if
!
                  qloc(j) = aquad(1,1,j)
                  qprob(j) = aquad(1,2,j)
              end if
!
              if (can) then
                  e2 = 0
              else
                  e2 = 1
              end if
!
              e1(1) = 0
              e1(2) = 0
              e1s = 0
              e1ss = 0
              e3(1) = 0
              e3(2) = 0
              e3s = 0
              m = 0
!
              do 40 k = 1,nilev
                  e7(k) = 0
                  e8(k) = 0
                  e8s(k) = 0
!
                  do 30 l = 1,k
                      e10(m+l) = 0
   30             end do
!
                  m = m+k
   40         end do
!
              irow = rowind(i) - 1
!
              do 90 t = 1,it(1,i)
                  irow = irow+1
                  yy = risk(irow)
!
                  if (j == 1) then
!
                      if (.not. myturn .and. t /= it(1,i)) then
                          go to 90
                      else if (.not. myturn) then
                          go to 200
                      end if
!
                      if (family(1) == 'p') then
                          pil(irow) = 0
!
                          do 1425 l = 1,nint(y(irow))
                              pil(irow) = pil(irow) + log(dble(l))
 1425                     end do
!
                      end if
!
                      sumxb(irow) = 0
!
                      do 48 k = 1,nilev
                          sumxb(irow) = sumxb(irow) + theta(k)*x(irow,k)
   48                 end do
!
                      if (offlag(1)) then
                          sumxb(irow) = sumxb(irow) + x(irow,offpos(1))
                      end if
!
                  end if
!
                  if (yy == 1) then
                      b = sumxb(irow) + theta(nest-1)*qloc(j)
                  else
                      b = sumxb(irow) + theta(nest)*qloc(j)
                  end if
!
                  if (family(1) /= 'g' .and. link(1) /= 'p') then
                      e = fexp(b,iflow)
!
                      if (iflow(1) .or. iflow(2)) then
                          ifail = .true.
                          go to 2000
                      end if
!
                  end if
!
                  pi = 3.14159265
!
                  if (family(1) == 'b') then
!
                      if (link(1) == 'c') then
!
                          if (e < 1d-8) then
                              expme = e
                          else
                              expme = 1 - fexp(-e,iflow)
!
                              if (iflow(1) .or. iflow(2)) then
                                  ifail = .true.
                                  go to 2000
                              end if
!
                          end if
!
                          if (y(irow) == 0) then
                              e2 = e2*(1-expme)
                              r = -e
                              rr = e
                          else
                              e2 = e2*expme
                              r = e*(1-expme)/expme
                              rr = r*(e-expme)/expme
                          end if
!
                      else if (link(1) == 'p') then
                          phi = erfc(-b/sqrt(2d0))/2
                          expme = fexp(-b**2/2,iflow)
!
                          if (iflow(1) .or. iflow(2)) then
                              ifail = .true.
                              go to 2000
                          end if
!
                          dpieta = expme/sqrt(2*pi)
!
                          if (y(irow) == 0) then
                              e2 = e2*(1-phi)
!
                              if (phi /= 1.0) then
                                  r = -dpieta/(1-phi)
                              else
                                  r = 1
                              end if
!
                          else
                              e2 = e2*phi
!
                              if (phi /= 0.0) then
                                  r = dpieta/phi
                              else
                                  r = 1
                              end if
!
                          end if
!
                          rr = r*(r+b)
                      else
!
                          if (y(irow) == 0) then
                              e2 = e2 - log(1+e)
                              r = -e/(1+e)
                          else
                              e2 = e2 + log(e/(1+e))
                              r = 1/(1+e)
                          end if
!
                          rr = e/(1+e)**2
                      end if
!
                  else if (family(1) == 'p') then
                      e2 = e2 + y(irow)*b - e - pil(irow)
                      r = y(irow) - e
                      rr = e
                  else if (family(1) == 'g') then
                      sigmae = theta(nilev+1)
                      e2 = e2 - ((y(irow) - b)**2/(sigmae**2) + &
                      log(2*pi*sigmae**2))/2
                      r = (y(irow) - b)/sigmae**2
                      rs = (r*(y(irow) - b) - 1)/sigmae
                      rr = 1/sigmae**2
                      rrs = 2*r/sigmae
                      rrss = rr*(3*r*(y(irow) - b) - 1)
                  end if
!
                  if (flag /= 1) then
                      e3(yy) = e3(yy) + r
!
                      if (family(1) == 'g') then
                          e3s = e3s+rs
                      end if
!
                      if (flag == 3) then
                          e1(yy) = e1(yy) + rr
!
                          if (family(1) == 'g') then
                              e1s = e1s+rrs
                              e1ss = e1ss+rrss
                          end if
!
                      end if
!
                      m = 0
!
                      do 70 k = 1,nilev
                          e7(k) = e7(k) + x(irow,k)*r
                          scores(irow,k) = x(irow,k)*r
!
                          if (flag == 3) then
                              e8(k) = e8(k) + x(irow,k)*rr
!
                              if (family(1) == 'g') then
                                  e8s(k) = e8s(k) + x(irow,k)*rrs
                              end if
!
                              do 60 l = 1,k
                                  e10(m+l) = e10(m+l) + &
                                  x(irow,k)*x(irow,l)*rr
   60                         end do
!
                              m = m+k
                          end if
!
   70                 end do
!
                  end if
!
   90         end do
!
              if (can) then
                  e2 = fexp(e2,iflow)
!
                  if (iflow(1) .or. iflow(2)) then
                      ifail = .true.
                      go to 2000
                  end if
!
              end if
!
              qpe2 = qprob(j)*e2
              listar = listar+qpe2
!
              if (flag /= 1) then
!
                  do 100 k = 1,nilev
                      d1(k) = d1(k) + qpe2*e7(k)
  100             end do
!
                  if (family(1) == 'g') then
                      d1(nilev+1) = d1(nilev+1) + qpe2*e3s
                  end if
!
                  d1(nest-1) = d1(nest-1) + qpe2*qloc(j)*e3(1)
                  d1(nest) = d1(nest) + qpe2*qloc(j)*e3(2)
              end if
!
              if (flag == 3) then
                  m = 0
!
                  do 120 k = 1,nilev
!
                      do 110 l = 1,k
                          d2(m+l) = d2(m+l) + &
                          qpe2*(e7(k)*e7(l) - e10(m+l))
  110                 end do
!
                      m = m+k
  120             end do
!
                  if (family(1) == 'g') then
!
                      do 1120 k = 1,nilev
                          d2(m+k) = d2(m+k) + qpe2*(e7(k)*e3s - e8s(k))
 1120                 end do
!
                      m = m+nilev
                      d2(m+1) = d2(m+1) + qpe2*(e3s**2 - e1ss)
                      m = m+1
                  end if
!
                  do 140 k = 1,n1lev
                      d2(m+k) = d2(m+k) + &
                      qpe2*qloc(j)*(e7(k)*e3(1) - e8(k))
  140             end do
!
                  do 1400 k = n1lev+1,nilev
                      d2(m+k) = d2(m+k) + qpe2*qloc(j)*e7(k)*e3(1)
 1400             end do
!
                  m = m+nilev
!
                  if (family(1) == 'g') then
                      d2(m+1) = d2(m+1) + qpe2*qloc(j)*(e3s*e3(1) - e1s)
                      m = m+1
                  end if
!
                  d2(m+1) = d2(m+1) + qpe2*qloc(j)**2*(e3(1)**2 - e1(1))
!
                  m = m+1
!
                  do 141 k = 1,n1lev
                      d2(m+k) = d2(m+k) + qpe2*qloc(j)*e7(k)*e3(2)
  141             end do
!
                  do 1410 k = n1lev+1,nilev
                      d2(m+k) = d2(m+k) + &
                      qpe2*qloc(j)*(e7(k)*e3(2) - e8(k))
 1410             end do
!
                  m = m+nilev
!
                  if (family(1) == 'g') then
                      d2(m+1) = d2(m+1) + qpe2*qloc(j)*(e3s*e3(2) - e1s)
                      m = m+1
                  end if
!
                  d2(m+1) = d2(m+1) + qpe2*qloc(j)**2*e3(1)*e3(2)
!
                  m = m+1
!
                  d2(m+1) = d2(m+1) + qpe2*qloc(j)**2*(e3(2)**2 - e1(2))
              end if
!
  150     end do
!
          if (abs(listar) > z1 .or. abs(listar) < z2) then
!
              if (abs(listar) < z2) then
                  iflow(1) = .true.
              else
                  iflow(2) = .true.
              end if
!
              ifail = .true.
              go to 2000
          end if
!
          if (flag /= 1) then
!
              do 160 k = 1,nest
                  der1(k) = d1(k)/listar
  160         end do
!
          end if
!
          if (flag == 2) then
!
              if (.not. sflag) then
!
                  do 170 k = 1,nest
                      deri(k) = der1(k)
  170             end do
!
                  sflag = .true.
              end if
!
              m = 0
!
              do 172 k = 1,nest
!
                  do 171 l = 1,k
                       hesinc(m+l) = hesinc(m+l) + &
                       (der1(k) - deri(k))*(der1(l) - deri(l))
  171             end do
!
                  m = m+k
  172         end do
!
          else if (flag == 3) then
              m = 0
!
              do 196 k = 1,nest
!
                  do 195 l = 1,k
                      hesinc(m+l) = hesinc(m+l) + &
                      d2(m+l)/listar - der1(k)*der1(l)
  195             end do
!
                  m = m+k
  196         end do
!
          end if
!
          if (i /= i1 .or. this_processor == boss_processor .or. &
          ind /= 0) then
!
              if (ind == 0) then
                  xll = xll + log(listar)
              else if (listar > 0) then
                  xil = log(listar)
                  sign = 1
              else
                  xil = log(-listar)
                  sign = -1
              end if
!
              if (flag /= 1) then
!
                  do 165 k = 1,nest
                      score(k) = score(k) + der1(k)
  165             end do
!
                  do 199 k = 1,ncov
                      hess(k) = hess(k) + hesinc(k)
  199             end do
!
              end if
!
          end if
!
 2000 continue
!
  200 end do
!
!$OMP END PARALLEL DO
!     end of outer loop
      end do
!
!---- Check to see if numerical failure has occurred
!---- Reduce all the ifails from each process, creating anyfail
!---- (only if there is more than one processor)
  300 if (num_processors > 1 .and. ind == 0) then
!
!---------******************
          call mpi_allreduce(ifail,anyfail,1,mpi_logical,mpi_lor, &
                             mpi_comm_world,ierror)
!---------******************
!
      else
          anyfail = ifail
      end if
!
      if (anyfail) then
!-------- Set the ifail flag for all processors and then return
          ifail = .true.
          return
      end if
!
!     Add in contributions to xll, score & hess from 1st case
!
      xll = xll + xll_first
!
      if (flag /= 1) then
!
          do k = 1,nest
              score(k) = score(k) + score_first(k)
          end do
!
          do k = 1,ncov
              hess(k) = hess(k) + hess_first(k)
          end do
!
      end if
!
      if (num_processors > 1 .and. ind == 0) then
!
!---------******************
          call mpi_allreduce(xll,xll_total,1,mpi_double_precision, &
                             mpi_sum,mpi_comm_world,ierror)
!---------******************
!
          xll = xll_total
!
          if (flag /= 1) then
!
!-------------******************
              call mpi_allreduce(score,score_total,nest, &
                                 mpi_double_precision,mpi_sum, &
                                 mpi_comm_world,ierror)
!-------------******************
!
              do 201 k = 1,nest
                  score(k) = score_total(k)
  201         end do
!
!-------------******************
              call mpi_allreduce(hess,hesinc,ncov,mpi_double_precision, &
                                 mpi_sum,mpi_comm_world,ierror)
!-------------******************
!
              do 202 k = 1,ncov
                  hess(k) = hesinc(k)
  202         end do
!
          end if
!
      end if
!
      if (flag == 2) then
!
          do 210 k = 1,nest
              meansc(k) = score(k)/nsub(1)
  210     end do
!
          m = 0
!
          do 230 k = 1,nest
!
              do 220 l = 1,k
                  hess(m+l) = hess(m+l) - nsub(1)* &
                  (deri(k) - meansc(k))*(deri(l) - meansc(l))
  220         end do
!
              m = m+k
  230    end do
!
      else if (flag == 3) then
!
          do 240 k = 1,ncov
              hess(k) = -hess(k)
  240     end do
!
      end if
!
      if (flag >= 2 .and. (iter == 1 .or. cflag == 2 .or. &
      cflag == 4) .and. tflag) then
!
          if (cflag /= 4 .and. step == 1.0) then
              write (outbuf,'(a,i4)') 'Iteration ',iter
          else if (cflag /= 4) then
              write (outbuf,'(a,i4,a,f7.4,a)') 'Iteration ',iter, &
              ' (step =',step,')'
          else
              write (outbuf,'(a,i4,a)') 'Iteration ',iter,' (converged)'
          end if
!
          call wrtlit(outbuf)
!
          do 9000 k = 1,nest
!
              if (k <= nilev) then
                  write (outbuf,'(a,i3,a,g14.4,a,f14.4,a)') 'beta(', &
                  k,') = ',dble(theta(k)),' [',score(k),']'
              else if (family(1) == 'g' .and. k == nilev+1) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'sigma_e   = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (k == nest-1) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') &
                  'scale1    = ',dble(theta(k)),' [',score(k),']'
              else if (k == nest) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') &
                  'scale2    = ',dble(theta(k)),' [',score(k),']'
              end if
!
              call wrtlit(outbuf)
!
              if (iter == 1 .or. k <= nest .or. theta(k) /= 0) &
              then
                  scnorm = scnorm + score(k)**2
              end if
!
 9000     end do
!
      end if
!
      if (flag >= 2 .and. (iter == 1 .or. cflag == 2 .or. &
      cflag == 4)) then
          scnorm = 0
!
          do 9500 k = 1,nest
              scnorm = scnorm + score(k)**2
 9500     end do
!
          scnorm = sqrt(scnorm)
!
          if (tflag) then
              write (outbuf,'(a,f14.4)') 'gradient norm = ',scnorm
              call wrtlit(outbuf)
              write (outbuf,'(a,f14.4)') 'log likelihood = ',xll
              call wrtlit(outbuf)
              call newlit
          end if
!
      end if
!
      if (xll .lt. -z1) then
          ifail = .true.
      end if
!
      return
!
      end subroutine lshdep
