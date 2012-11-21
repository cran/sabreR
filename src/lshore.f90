!
!***********************************************************************
!
      subroutine lshore(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                        score,hess,ncov,flag,ifail,offlag,offpos,maxcol, &
                        link,iter,cflag,step,scnorm,scores,ncat,ind,xil, &
                        sign,iquad,aquad,mus,taus,maxcat)
      use accurate_arithmetic
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character link(3)
      integer nmes,nest,ncov,maxcol,nsub(2),it(2,nsub(1)),flag,nilev, &
              iter,cflag,ncat(3),offpos(3),nm(3),ind,sign,iquad,maxcat
      double precision y(nmes),x(nmes,maxcol),beta(nest),xll, &
             score(nest),hess(ncov),scnorm,step,scores(nmes,nest), &
             quad(3,2,256),xil,aquad(3,2,256),mus(3,nsub(1)), &
             taus(3,nsub(1))
      logical ifail,offlag(3)
!-----------------------------------------------------------------------
!     function : calculates deviance, score vector and hessian matrix at
!                 the current parameter estimates, for random effects
!                 ordered response model.
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     hesinc is an array to store the terms added to the hessian for one
!     case
      character(len=80) :: outbuf
      character(len=5) :: str
      double precision theta(maxpar),e2,r,rr,b,pi,bb,meansc(maxpar), &
             sumxb(nmes),expme,deri(maxpar),e,listar, &
             d2(maxpar*(maxpar+1)/2),d1(maxpar),der1(maxpar), &
             phi(0:maxcat),dpieta(0:maxcat),phi2(0:maxcat),qloc(256), &
             qprob(256),qpe2,e1,e3,e4(maxcat),e7(maxpar), &
             e8(maxpar),e10(maxpar*(maxpar+1)/2),d,rd(maxcat), &
             rrd(maxcat),rrdd(maxcat,maxcat),sr1,sr2,s1,s2,srd1(maxcat), &
             srd2(maxcat),e5(maxpar*maxcat),e6(maxcat*(maxcat+1)/2), &
             e9(maxcat)
      double precision hesinc(maxpar*(maxpar+1)/2), &
             score_total(maxpar),xll_total
      integer t,i,irow,l,j,k,m,rowind(nsub(1)),c,cat,ncut,jj,n,nq,nc,i1, &
              i2
      logical sflag,anyfail,myturn,iflow(2)
      double precision :: xll_first,score_first(nest),hess_first(ncov)
      integer :: outer,firstcase,lastcase
!-----------------------------------------------------------------------
!     mpi parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
      xll = 0
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
      ncut = ncat(1) - 1
      c = 1
!
      if (link(1) == 'p') then
          c = -1
      end if
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
!$OMP SHARED(num_processors,nest,ncov,nm,iquad,quad,mus,taus,nilev), &
!$OMP SHARED(ncut,ncat,rowind,it,theta,x,offlag,offpos,y,link,flag), &
!$OMP SHARED(sflag,xil,sign,c,scores), &
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
              e1 = 0
              e2 = 0
              e3 = 0
              m = 0
              n = 0
!
              do 40 k = 1,nilev
                  e7(k) = 0
                  e8(k) = 0
!
                  do 30 l = 1,k
                      e10(m+l) = 0
   30             end do
!
                  m = m+k
!
                  do 35 l = 1,ncut
                      e5(n+l) = 0
   35             end do
!
                  n = n+ncut
   40         end do
!
              m = 0
!
              do 45 k = 1,ncut
                  e4(k) = 0
                  e9(k) = 0
!
                  do 44 l = 1,k
                      e6(m+l) = 0
   44             end do
!
                  m = m+k
   45         end do
!
              irow = rowind(i) - 1
!
              do 90 t = 1,it(1,i)
                  irow = irow+1
!
                  if (j == 1) then
!
                      if (.not. myturn .and. t /= it(1,i)) then
                          go to 90
                      else if (.not. myturn) then
                          go to 200
                      end if
!
                      sumxb(irow) = 0
!
                      do 50 k = 1,nilev
                          sumxb(irow) = sumxb(irow) + theta(k)*x(irow,k)
   50                 end do
!
                      if (offlag(1)) then
                          sumxb(irow) = sumxb(irow) + x(irow,offpos(1))
                      end if
!
                  end if
!
                  b = sumxb(irow) + theta(nilev+ncat(1))*qloc(j)
                  e = fexp(b,iflow)
!
                  if (iflow(1) .or. iflow(2)) then
                      ifail = .true.
                      go to 2000
                  end if
!
                  phi(0) = 0
                  phi(ncat(1)) = 1
                  dpieta(0) = 0
                  dpieta(ncat(1)) = 0
                  phi2(0) = -1
                  phi2(ncat(1)) = 1
                  pi = 3.14159265
!
                  do 52 cat = 1,ncat(1)
!
                      if (y(irow) == cat) then
                          go to 53
                      end if
!
   52             end do
!
   53             do 55 jj = max(cat-1,1),min(cat,ncut)
                      bb = theta(nilev+jj) - b
!
                      if (link(1) == 'p') then
                          phi(jj) = erfc(-bb/sqrt(2d0))/2
                          expme = fexp(-bb**2/2,iflow)
!
                          if (iflow(1) .or. iflow(2)) then
                              ifail = .true.
                              go to 2000
                          end if
!
                          dpieta(jj) = expme/sqrt(2*pi)
                          phi2(jj) = bb
                      else
                          phi(jj) = 1/(1 + fexp(-bb,iflow))
!
                          if (iflow(1) .or. iflow(2)) then
                              ifail = .true.
                              go to 2000
                          end if
!
                          dpieta(jj) = -phi(jj)*(1 - phi(jj))
                          phi2(jj) = 2*phi(jj) - 1
                      end if
!
   55             end do
!
                  d = phi(cat) - phi(cat-1)
!
                  if (d == 0) then
                      d = sqrt(z2)
                  end if
!
                  e2 = e2 + log(d)
                  r = (dpieta(cat) - dpieta(cat-1))/d
                  sr1 = dpieta(cat-1)*phi2(cat-1)
                  sr2 = dpieta(cat)*phi2(cat)
                  rr = r**2 + c*(sr1-sr2)/d
!
                  do 68 k = 1,ncut
                      s1 = 0
                      s2 = 0
                      srd1(k) = 0
                      srd2(k) = 0
!
                      if (k == cat-1) then
                          s1 = dpieta(k)
                          srd1(k) = sr1
                      else if (k == cat) then
                          s2 = dpieta(k)
                          srd2(k) = sr2
                      end if
!
                      rd(k) = (s1-s2)/d
                      rrd(k) = r*rd(k) + c*(srd2(k) - srd1(k))/d
!
                      do 1452 l = 1,k
                          s1 = 0
                          s2 = 0
!
                          if (k == cat-1) then
                              s1 = srd1(l)
                          else if (k == cat) then
                              s2 = srd2(l)
                          end if
!
                          rrdd(k,l) = rd(k)*rd(l) + c*(s1-s2)/d
 1452                 end do
!
   68             end do
!
                  if (flag /= 1) then
                      e3 = e3+r
!
                      if (flag == 3) then
                          e1 = e1+rr
                      end if
!
                      m = 0
!
                      do 100 k = 1,nilev
                          e7(k) = e7(k) + x(irow,k)*r
                          scores(irow,k) = x(irow,k)*r
!
                          if (flag == 3) then
                              e8(k) = e8(k) + x(irow,k)*rr
!
                              do 110 l = 1,k
                                  e10(m+l) = e10(m+l) + &
                                  x(irow,k)*x(irow,l)*rr
  110                         end do
!
                              m = m+k
                          end if
!
  100                 end do
!
                      m = 0
                      n = 0
!
                      do 98 k = 1,ncut
                          e4(k) = e4(k) + rd(k)
!
                          if (flag == 3) then
                              e9(k) = e9(k) + rrd(k)
!
                              do 130 l = 1,nilev
                                  e5(m+l) = e5(m+l) + x(irow,l)*rrd(k)
  130                         end do
!
                              m = m+nilev
!
                              do 145 l = 1,k
                                  e6(n+l) = e6(n+l) + rrdd(k,l)
  145                         end do
!
                              n = n+k
                          end if
!
   98                 end do
!
                  end if
!
   90         end do
!
              e2 = fexp(e2,iflow)
!
              if (iflow(1) .or. iflow(2)) then
                  ifail = .true.
                  go to 2000
              end if
!
              qpe2 = qprob(j)*e2
              listar = listar+qpe2
!
              if (flag /= 1) then
!
                  do 1001 k = 1,nilev
                      d1(k) = d1(k) + c*qpe2*e7(k)
 1001             end do
!
                  do 981 k = 1,ncut
                      d1(nilev+k) = d1(nilev+k) + c*qpe2*e4(k)
  981             end do
!
                  d1(nest) = d1(nest) + c*qpe2*qloc(j)*e3
              end if
!
              if (flag == 3) then
                  m = 0
!
                  do 1201 k = 1,nilev
!
                      do 1101 l = 1,k
                          d2(m+l) = d2(m+l) + &
                          qpe2*(e7(k)*e7(l) - e10(m+l))
 1101                 end do
!
                      m = m+k
 1201             end do
!
                  n = 0
!
                  do 1411 k = 1,ncut
!
                      do 1301 l = 1,nilev
                          d2(m+l) = d2(m+l) + &
                          qpe2*(e7(l)*e4(k) - e5(n+l))
 1301                 end do
!
                      m = m+nilev+k
                      n = n+nilev
 1411             end do
!
                  m = nilev*(nilev+3)/2
                  n = 0
!
                  do 1511 k = 1,ncut
!
                      do 1451 l = 1,k
                          d2(m+l) = d2(m+l) + &
                          qpe2*(e4(k)*e4(l) - e6(n+l))
 1451                 end do
!
                      m = m+nilev+k
                      n = n+k
 1511             end do
!
                  m = (nest-1)*nest/2
!
                  do 1600 k = 1,nilev
                      d2(m+k) = d2(m+k) + &
                      qpe2*qloc(j)*(e7(k)*e3 - e8(k))
 1600             end do
!
                  m = m+nilev
!
                  do 1610 k = 1,ncut
                      d2(m+k) = d2(m+k) + &
                      qpe2*qloc(j)*(e4(k)*e3 - e9(k))
 1610             end do
!
                  m = m+ncut
!
                  d2(m+1) = d2(m+1) + qpe2*qloc(j)**2*(e3**2 - e1)
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
                  write (outbuf,'(a,i3,a,g14.4,a,f14.4,a)') 'beta(',k, &
                  ') = ',dble(theta(k)),' [',score(k),']'
              else if (k /= nest) then
                  write (outbuf,'(a,i3,a,g14.4,a,f14.4,a)') 'cut(', &
                  k-nilev,') = ',dble(theta(k)),' [',score(k),']'
              else
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'scale = ', &
                  dble(theta(k)),' [',score(k),']'
              end if
!
              call wrtlit(outbuf)
              scnorm = scnorm + score(k)**2
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
      end subroutine lshore
