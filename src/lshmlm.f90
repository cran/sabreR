!
!***********************************************************************
!
      subroutine lshmlm(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                        score,hess,ncov,flag,ifail,offlag,offpos,maxcol, &
                        link,iter,cflag,step,scnorm,scores,family,maxit, &
                        ind,xil,sign,iquad,aquad,mus,taus)
      use accurate_arithmetic
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character link(3),family(3)
      integer nmes,nest,ncov,maxcol,nsub(2),nm(3),it(2,nsub(1)),flag, &
              nilev,iter,cflag,offpos(3),maxit,ind,sign,iquad
      double precision y(nmes),x(nmes,maxcol),beta(nest),xll,xil, &
             score(nest),hess(ncov),quad(3,2,256),scnorm,step, &
             scores(nmes,nest),aquad(3,2,256),mus(3,nsub(1)), &
             taus(3,nsub(1))
      logical ifail,offlag(3)
!-----------------------------------------------------------------------
!     function : calculates deviance, score vector and hessian matrix at
!                 the current parameter estimates, for multilevel model.
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     hesinc is an array to store the terms added to the hessian for one
!     case
      character(len=80) :: outbuf
      character(len=5) :: str
      double precision theta(maxpar),qprob(2,256),e2,r,rr,e8(maxpar),b, &
             pi,e10(maxpar*(maxpar+3)/2),meansc(maxpar),e1,pil(nmes), &
             sumxb(nmes),e3,expme,deri(maxpar),e,phi,listar, &
             d2(maxpar*(maxpar+1)/2),d1(maxpar), &
             qloc(2,256),e7(maxpar),der1(maxpar),dpieta,e3s,rs,e1s, &
             e1ss,rrs,rrss,e8s(maxpar),qpe2,qpxl, &
             li2star,xl1,dd1(maxpar),dd2(maxpar*(maxpar+1)/2), &
             a1(maxpar),a2(maxpar*(maxpar+1)/2),a3(maxpar*(maxpar+1)/2)
      double precision score_total(maxpar),xll_total, &
             hesinc(maxpar*(maxpar+1)/2)
      integer t,irow,l,k,m,rowind(nsub(2),maxit),j1,j2,i1,l1,l2, &
              i2,nit(nmes),nq,nc
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
      do 10 k = 1,nm(1)
          qloc(1,k) = quad(1,1,k)
          qprob(1,k) = quad(1,2,k)
   10 end do
!
      do 11 k = 1,nm(2)
          qloc(2,k) = quad(2,1,k)
          qprob(2,k) = quad(2,2,k)
   11 end do
!
      xll = 0
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
      irow = 1
      nit(1) = 0
!
      do 1200 i2 = 1,nsub(2)
!
          if (i2 /= 1) then
              nit(i2) = nit(i2-1) + it(2,i2-1)
          end if
!
          do 1190 i1 = 1,it(2,i2)
              rowind(i2,i1) = irow
              irow = irow + it(1,nit(i2) + i1)
 1190     end do
!
 1200 end do
!
      nq = 0
      nc = 0
!
      if (ind == 0) then
          l1 = 1
          l2 = nsub(2)
      else
          l1 = ind
          l2 = l1
      end if
!
!---- loop around cases
!
!     Split the 5010 loop into the first case (outer=1) and all the rest
!     (outer=2)
!
      do outer = 1,2
!
      if (outer == 1) then
          firstcase = l1
          lastcase = l1
      else
          firstcase = l1 + 1
          lastcase = l2
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
!$OMP SHARED(firstcase,lastcase,l1,l2,ind,this_processor), &
!$OMP SHARED(num_processors,nest,ncov,nm,iquad,quad,mus,taus,can,nilev), &
!$OMP SHARED(rowind,it,nit,family,y,theta,x,offlag,offpos,link,flag), &
!$OMP SHARED(scores,sflag,xil,sign,deri,aquad), &
!$OMP REDUCTION(+: xll,score,hess), REDUCTION(.or.: ifail)
!
      do 5010 i2 = firstcase,lastcase
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
!          if (i2 == l1 .or. ind /= 0 .or. iquad > 0) then
          if (i2 == l1 .or. ind /= 0) then
              myturn = .true.
          else
              myturn = this_processor == mod(i2-2,num_processors)
          end if
!
          if (.not. myturn .and. i2 /= l2) then
              go to 5010
          end if
!
          ifail = .false.
          iflow(1) = .false.
          iflow(2) = .false.
!
          li2star = 0
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
          do 4010 j2 = 1,nm(2)
!!
              if (iquad == 0 .and. i2 == l1) then
                  qloc(2,j2) = quad(2,1,j2)
                  qprob(2,j2) = quad(2,2,j2)
              else if (iquad /= 0) then
!
                  if (ind == 0  .and. outer == 1) then
                      aquad(2,1,j2) = mus(2,i2) + &
                      taus(2,i2)*quad(2,1,j2)
                      aquad(2,2,j2) = taus(2,i2)*quad(2,2,j2)* &
                      exp((quad(2,1,j2)**2 - aquad(2,1,j2)**2)/2)
                  end if

                  qloc(2,j2) = aquad(2,1,j2)
                  qprob(2,j2) = aquad(2,2,j2)
              end if
!
              xl1 = 1
!
              do 2300 k = 1,nest
                  a1(k) = 0
 2300         end do
!
              do 2310 k = 1,ncov
                  a2(k) = 0
                  a3(k) = 0
 2310         end do
!
              do 200 i1 = 1,it(2,i2)
                  listar = 0
!
                  do 2400 k = 1,nest
                      dd1(k) = 0
 2400             end do
!
                  do 2410 k = 1,ncov
                      dd2(k) = 0
 2410             end do
!
                  do 150 j1 = 1,nm(1)
!                      qloc(1,j1) = quad(1,1,j1)
!                      qprob(1,j1) = quad(1,2,j1)
!
                      if (can) then
                          e2 = 0
                      else
                          e2 = 1
                      end if
!
                      e1 = 0
                      e1s = 0
                      e1ss = 0
                      e3 = 0
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
   30                     end do
!
                          m = m+k
   40                 end do
!
                      irow = rowind(i2,i1) - 1
!
                      do 90 t = 1,it(1,nit(i2) + i1)
                          irow = irow+1
!
                          if (j1 == 1 .and. j2 == 1) then
!
                              if (.not. myturn .and. &
                              t /= it(1,nit(i2) + i1)) then
                                  go to 90
                              else if (.not. myturn) then
                                  go to 5010
                              end if
!
                              if (family(1) == 'p') then
                                  pil(irow) = 0
!
                                  do 1425 l = 1,nint(y(irow))
                                      pil(irow) = pil(irow) + &
                                      log(dble(l))
 1425                             end do
!
                              end if
!
                              sumxb(irow) = 0
!
                              do 48 k = 1,nilev
                                  sumxb(irow) = sumxb(irow) + &
                                  theta(k)*x(irow,k)
   48                         end do
!
                              if (offlag(1)) then
                                  sumxb(irow) = sumxb(irow) + &
                                  x(irow,offpos(1))
                              end if
!
                          end if
!
                          b = sumxb(irow) + &
                          theta(nest-1)*qloc(1,j1) + &
                          theta(nest)*qloc(2,j2)
!
                          if (family(1) /= 'g' .and. link(1) /= 'p') &
                          then
                              e = fexp(b,iflow)
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
                              e2 = e2 - ((y(irow) - b)**2/ &
                              (theta(nilev+1)**2) + &
                              log(2*pi*theta(nilev+1)**2))/2
                              r = (y(irow) - b)/theta(nilev+1)**2
                              rs = (r*(y(irow) - b) - 1)/theta(nilev+1)
                              rr = 1/theta(nilev+1)**2
                              rrs = 2*r/theta(nilev+1)
                              rrss = rr*(3*r*(y(irow) - b) - 1)
                          end if
!
                          if (flag /= 1) then
                              e3 = e3+r
!
                              if (family(1) == 'g') then
                                  e3s = e3s+rs
                              end if
!
                              if (flag == 3) then
                                  e1 = e1+rr
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
                                          e8s(k) = e8s(k) + &
                                          x(irow,k)*rrs
                                      end if
!
                                      do 60 l = 1,k
                                          e10(m+l) = e10(m+l) + &
                                          x(irow,k)*x(irow,l)*rr
   60                                 end do
!
                                      m = m+k
                                  end if
!
   70                         end do
!
                          end if
!
   90                 end do
!
                      if (can) then
                          e2 = fexp(e2,iflow)
                      end if
!
                      qpe2 = qprob(1,j1)*e2
                      listar = listar+qpe2
!
                      if (flag /= 1) then
!
                          do 100 k = 1,nilev
                              dd1(k) = dd1(k) + qpe2*e7(k)
  100                     end do
!
                          if (family(1) == 'g') then
                              dd1(nilev+1) = dd1(nilev+1) + qpe2*e3s
                          end if
!
                          dd1(nest-1) = dd1(nest-1) + qpe2*qloc(1,j1)*e3
!
                          dd1(nest) = dd1(nest) + qpe2*qloc(2,j2)*e3
!
                      end if
!
                      if (flag == 3) then
                          m = 0
!
                          do 120 k = 1,nilev
!
                              do 110 l = 1,k
                                  dd2(m+l) = dd2(m+l) + &
                                  qpe2*(e7(k)*e7(l) - e10(m+l))
  110                         end do
!
                              m = m+k
  120                     end do
!
                          if (family(1) == 'g') then
!
                              do 1125 k = 1,nilev
                                  dd2(m+k) = dd2(m+k) + &
                                  qpe2*(e7(k)*e3s - e8s(k))
 1125                         end do
!
                              m = m+nilev
!
                              dd2(m+1) = dd2(m+1) + qpe2*(e3s**2 - e1ss)
!
                              m = m+1
                          end if
!
                          do 140 k = 1,nilev
                              dd2(m+k) = dd2(m+k) + &
                              qpe2*qloc(1,j1)*(e7(k)*e3 - e8(k))
  140                     end do
!
                          m = m+nilev
!
                          if (family(1) == 'g') then
                              dd2(m+1) = dd2(m+1) + &
                              qpe2*qloc(1,j1)*(e3s*e3 - e1s)
                              m = m+1
                          end if
!
                          dd2(m+1) = dd2(m+1) + &
                          qpe2*qloc(1,j1)**2*(e3**2 - e1)
!
                          m = m+1
!
                          do 145 k = 1,nilev
                              dd2(m+k) = dd2(m+k) + &
                              qpe2*qloc(2,j2)*(e7(k)*e3 - e8(k))
  145                     end do
!
                          m = m+nilev
!
                          if (family(1) == 'g') then
                              dd2(m+1) = dd2(m+1) + &
                              qpe2*qloc(2,j2)*(e3s*e3 - e1s)
                              m = m+1
                          end if
!
                          dd2(m+1) = dd2(m+1) + &
                          qpe2*qloc(1,j1)*qloc(2,j2)*(e3**2 - e1)
!
                          m = m+1
!
                          dd2(m+1) = dd2(m+1) + &
                          qpe2*qloc(2,j2)**2*(e3**2 - e1)
!
                      end if
!
  150             end do
!
                  xl1 = xl1*listar
!
                  if (flag /= 1) then
!
                      do 180 k = 1,nest
                          a1(k) = a1(k) + dd1(k)/listar
  180                 end do
!
                      if (flag == 3) then
!
                          m = 0
!
                          do 192 k = 1,nest
!
                              do 190 l = 1,k
                                  a2(m+l) = a2(m+l) + dd2(m+l)/listar
                                  a3(m+l) = a3(m+l) + &
                                  dd1(k)*dd1(l)/(listar**2)
  190                         end do
!
                              m = m+k
  192                     end do
!
                      end if
!
                  end if
!
  200         end do
!
              qpxl = qprob(2,j2)*xl1
              li2star = li2star+qpxl
!
              if (flag /= 1) then
!
                  do 3000 k = 1,nest
                      d1(k) = d1(k) + qpxl*a1(k)
 3000             end do
!
                  if (flag == 3) then
                      m = 0
!
                      do 3050 k = 1,nest
!
                          do 3040 l = 1,k
                              d2(m+l) = d2(m+l) + &
                              qpxl*(a2(m+l) + a1(k)*a1(l) - a3(m+l))
 3040                     end do
!
                          m = m+k
 3050                 end do
!
                  end if
!
              end if
!
 4010     end do
!
          if (flag /= 1) then
!
              do 160 k = 1,nest
                  der1(k) = d1(k)/li2star
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
                      d2(m+l)/li2star - der1(k)*der1(l)
  195             end do
!
                  m = m+k
  196         end do
!
          end if
!
          if (i2 > 1 .or. this_processor == boss_processor) then
!
              if (ind == 0) then
                  xll = xll + log(li2star)
              else if (li2star > 0) then
                  xil = log(li2star)
                  sign = 1
              else
                  xil = log(-li2star)
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
 5010 end do
!
!$OMP END PARALLEL DO
!     end of outer loop
      end do
!
!---- Check to see if numerical failure has occurred
!---- Reduce all the ifails from each process, creating anyfail
!---- (only if there is more than one processor)
      if (num_processors > 1 .and. iquad == 0) then
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
      if (num_processors > 1 .and. iquad == 0) then
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
              meansc(k) = score(k)/nsub(2)
  210     end do
!
          m = 0
!
          do 230 k = 1,nest
!
              do 220 l = 1,k
                  hess(m+l) = hess(m+l) - nsub(2)* &
                  (deri(k) - meansc(k))*(deri(l) - meansc(l))
  220         end do
!
              m = m+k
  230     end do
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
              else if (family(1) == 'g' .and. k == nilev+1) &
              then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'sigma_e   = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (k == nest-1) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') &
                  'scale2    = ',dble(theta(k)),' [',score(k),']'
              else if (k == nest) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') &
                  'scale3    = ',dble(theta(k)),' [',score(k),']'
              end if
!
              call wrtlit(outbuf)
!
              if (iter == 1 .or. k <= nest .or. &
              theta(k) /= 0) then
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
      end subroutine lshmlm
