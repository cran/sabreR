!
!***********************************************************************
!
      subroutine lstore_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                            nm,score,hess,ncov,flag,ifail,maxcol,link, &
                            risk,corr,iter,cflag,step,scnorm,scores, &
                            n1lev,n2lev,family,ind,xil,sign,iquad,aquad, &
                            mus,taus,offlag,offpos,ncat,maxcat)
      use accurate_arithmetic
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character link(3),family(3)
      integer nmes,nest,ncov,maxcol,nsub(2),nm(3),it(2,nsub(1)),flag, &
              nilev,risk(nmes),iter,cflag,n1lev,ind,sign,iquad,n2lev, &
              offpos(3),ncat(3),maxcat
      double precision y(nmes),x(nmes,maxcol),beta(nest),xll,xil, &
             score(nest),hess(ncov),quad(3,2,256),scnorm,step, &
             scores(nmes,nest),aquad(3,2,256),mus(3,nsub(1)), &
             taus(3,nsub(1))
      logical ifail,corr,offlag(3)
!-----------------------------------------------------------------------
!     function : calculates deviance, score vector and hessian matrix at
!                 the current parameter estimates, for bivariate model.
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     hesinc is an array to store the terms added to the hessian for one
!     case
      character(len=80) :: outbuf
      character(len=5) :: str
      type (accurate) theta(maxpar),thetanir,qprob(3,256),e2,r,e2s(3), &
             b,pi,pil(nmes),bb,d,s1,s2,meansc(maxpar), &
             rd(maxcat),sumxb(nmes),e3(3),expme,deri(maxpar),e,phi, &
             listar,d1(maxpar),qloc(3,256),e7(3,maxpar),der1(maxpar), &
             dpieta,rs,e3s(3),qpe2,e2s2,phi1(0:maxcat), &
             dpieta1(0:maxcat),e4(3,maxcat)
      double precision hesinc(maxpar*(maxpar+1)/2),score_total(maxpar), &
             xll_total
      integer t,i,irow,l,j1,k,m,j2,rowind(nsub(1)),ir,i1,i2,k1,k2,j3,ic, &
              ncut(3),c(3),cat,jj,n,nic1,nic2,nic3
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
      sflag = .false.
!
      ic = 1
!
      if (corr) then
          ic = -1
      end if
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
!
      do 1200 i = 1,nsub(1)
          rowind(i) = irow
          irow = irow + it(1,i)
 1200 end do
!
      if (ind == 0) then
          i1 = 1
          i2 = nsub(1)
      else
          i1 = ind
          i2 = i1
      end if
!
      ncut(1) = ncat(1) - 1
      ncut(2) = ncat(2) - 1
      ncut(3) = ncat(3) - 1
      nic1 = nilev + ncut(1)
      nic2 = nic1 + ncut(2)
      nic3 = nic2 + ncut(3)
      c(1) = 1
      c(2) = 1
      c(3) = 1
!
      if (link(1) == 'p') then
          c(1) = -1
      end if
!
      if (link(2) == 'p') then
          c(2) = -1
      end if
!
      if (link(3) == 'p') then
          c(3) = -1
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
!     values of xll, score & hess after the first case are saved and the
!     variables are reset to zero for the parallel part (outer=2). the
!     initial values are then added after the final case. this is to
!     avoid the values after the first iteration being multiplied by the
!     number of threads in the reduction if openmp is being used
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
!$OMP SHARED(rowind,it,family,y,theta,x,offlag,offpos,link,flag,scores), &
!$OMP SHARED(ic,z1,z2,sflag,xil,sign,c,risk,n1lev,n2lev,corr,nmes,ncut), &
!$OMP SHARED(ncat,nic1,nic2,nic3), &
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
              hesinc(k) = 0
   25     end do
!
          do 150 j1 = 1,nm(1)
!
              if (iquad == 0 .and. i == i1) then
                  qloc(1,j1) = quad(1,1,j1)
                  qprob(1,j1) = quad(1,2,j1)
              else if (iquad /= 0) then
!
                  if (ind == 0 .and. outer == 1) then
                      aquad(1,1,j1) = mus(1,i) + taus(1,i)*quad(1,1,j1)
                      aquad(1,2,j1) = taus(1,i)*quad(1,2,j1)* &
                      exp((quad(1,1,j1)**2 - aquad(1,1,j1)**2)/2)
                  end if
!
                  qloc(1,j1) = aquad(1,1,j1)
                  qprob(1,j1) = aquad(1,2,j1)
              end if
!
              do 148 j2 = 1,nm(2)
!
                  if (iquad == 0 .and. i == i1) then
                      qloc(2,j2) = quad(2,1,j2)
                      qprob(2,j2) = quad(2,2,j2)
                  else if (iquad /= 0) then
!
                      if (ind == 0  .and. outer == 1) then
                          aquad(2,1,j2) = mus(2,i) + &
                          taus(2,i)*quad(2,1,j2)
                          aquad(2,2,j2) = taus(2,i)*quad(2,2,j2)* &
                          exp((quad(2,1,j2)**2 - aquad(2,1,j2)**2)/2)
                      end if
!
                      qloc(2,j2) = aquad(2,1,j2)
                      qprob(2,j2) = aquad(2,2,j2)
                  end if
!
                  do 145 j3 = 1,nm(3)
!
                      if (iquad == 0 .and. i == i1) then
                          qloc(3,j3) = quad(3,1,j3)
                          qprob(3,j3) = quad(3,2,j3)
                      else if (iquad /= 0) then
!
                          if (ind == 0  .and. outer == 1) then
                              aquad(3,1,j3) = mus(3,i) + &
                              taus(3,i)*quad(3,1,j3)
                              aquad(3,2,j3) = taus(3,i)*quad(3,2,j3)* &
                              exp((quad(3,1,j3)**2 - aquad(3,1,j3)**2)/2)
                          end if
!
                          qloc(3,j3) = aquad(3,1,j3)
                          qprob(3,j3) = aquad(3,2,j3)
                      end if
!
                      e2 = 0
                      e2s(1) = 0
                      e2s(2) = 0
                      e2s(3) = 0
                      e3(1) = 0
                      e3(2) = 0
                      e3(3) = 0
                      e3s(1) = 0
                      e3s(2) = 0
                      e3s(3) = 0
!
                      do 40 k = 1,nilev
                          e7(1,k) = 0
                          e7(2,k) = 0
                          e7(3,k) = 0
   40                 end do
!
                      do 45 k = 1,ncut(1)
                          e4(1,k) = 0
   45                 continue
!
                      do 46 k = 1,ncut(2)
                          e4(2,k) = 0
   46                 end do
!
                      do 47 k = 1,ncut(3)
                          e4(3,k) = 0
   47                 end do
!
                      irow = rowind(i) - 1
!
                      do 90 t = 1,it(1,i)
                          irow = irow+1
                          ir = risk(irow)
!
                          if (j1 == 1 .and. j2 == 1 .and. j3 == 1) then
!
                              if (.not. myturn .and. t /= it(1,i)) then
                                  go to 90
                              else if (.not. myturn) then
                                  go to 200
                              end if
!
                              sumxb(irow) = 0
!
                              if (ir == 1) then
                                  k1 = 1
                                  k2 = n1lev
                              else if (ir == 2) then
                                  k1 = n1lev+1
                                  k2 = n1lev+n2lev
                              else
                                  k1 = n1lev+n2lev+1
                                  k2 = nilev
                              end if
!
                              do 50 k = k1,k2
                                  sumxb(irow) = sumxb(irow) + &
                                  theta(k)*x(irow,k)
   50                         end do
!
                              if (offlag(ir)) then
                                  sumxb(irow) = sumxb(irow) + &
                                  x(irow,offpos(ir))
                              end if
!
                          end if
!
                          if (ir == 1) then
                              b = sumxb(irow) + theta(nic3+1)*qloc(1,j1)
                          else if (ir == 2 .and. .not. corr) then
                              b = sumxb(irow) + theta(nic3+2)*qloc(2,j2)
                          else if (.not. corr) then
                              b = sumxb(irow) + theta(nic3+3)*qloc(3,j3)
                          else if (ir == 2) then
                              b = sumxb(irow) + theta(nic3+2)* &
                              (theta(nic3+4)*qloc(1,j1) + &
                              sqrt(1 - theta(nic3+4)**2)*qloc(2,j2))
                          else
                              b = sumxb(irow) + theta(nic3+3)* &
                              (theta(nic3+5)*qloc(1,j1) + (theta(nic3+6) &
                              - theta(nic3+4)*theta(nic3+5))*qloc(2,j2)/ &
                              sqrt(1 - theta(nic3+4)**2) + &
                              sqrt((1 - theta(nic3+4)**2 - &
                              theta(nic3+5)**2 - theta(nic3+6)**2 + &
                              2*theta(nic3+4)*theta(nic3+5)* &
                              theta(nic3+6))/(1 - theta(nic3+4)**2))* &
                              qloc(3,j3))
                          end if
!
                          e = fexp(b,iflow)
!
                          if (iflow(1) .or. iflow(2)) then
                              ifail = .true.
                              go to 2000
                          end if
!
                          pi = 3.14159265
                          phi1(0) = 0
                          phi1(ncat(ir)) = 1
                          dpieta1(0) = 0
                          dpieta1(ncat(ir)) = 0
!
                          do 52 cat = 1,ncat(ir)
!
                              if (y(irow) == cat) then
                                  go to 53
                              end if
!
   52                     end do
!
   53                     do 55 jj = max(cat-1,1),min(cat,ncut(ir))
!
                              if (ir == 1) then
                                  bb = theta(nilev+jj) - b
                              else if (ir == 2) then
                                  bb = theta(nic1+jj) - b
                              else
                                  bb = theta(nic2+jj) - b
                              end if
!
                              if (link(ir) == 'p') then
                                  phi1(jj) = erfc(-bb/sqrt(2d0))/2
                                  expme = fexp(-bb**2/2,iflow)
!
                                  if (iflow(1) .or. iflow(2)) then
                                      ifail = .true.
                                      go to 2000
                                  end if
!
                                  dpieta1(jj) = expme/sqrt(2*pi)
                              else
                                  phi1(jj) = 1/(1 + fexp(-bb,iflow))
!
                                  if (iflow(1) .or. iflow(2)) then
                                      ifail = .true.
                                      go to 2000
                                  end if
!
                                  dpieta1(jj) = -phi1(jj)*(1 - phi1(jj))
                              end if
!
   55                     end do
!
                          d = phi1(cat) - phi1(cat-1)
!
                          if (d == 0) then
                              d = sqrt(z2)
                          end if
!
                          e2 = e2 + log(d)
                          e2s(ir) = e2s(ir) + log(d)
                          r = (dpieta1(cat) - dpieta1(cat-1))/d
!
                          do 68 k = 1,ncut(ir)
                              s1 = 0
                              s2 = 0
!
                              if (k == cat-1) then
                                  s1 = dpieta1(k)
                              else if (k == cat) then
                                  s2 = dpieta1(k)
                              end if
!
                              rd(k) = (s1-s2)/d
   68                     end do
!
                          if (flag /= 1) then
                              e3(ir) = e3(ir) + r
!
                              do 70 k = 1,nilev
                                  e7(ir,k) = e7(ir,k) + x(irow,k)*r
                                  scores(i,k) = x(irow,k)*r
   70                         end do
!
                              do 98 k = 1,ncut(ir)
                                  e4(ir,k) = e4(ir,k) + rd(k)
   98                         end do
!
                          end if
!
   90                 end do
!
                      e2s2 = e2s(1) + e2s(2) + e2s(3)
                      e2 = fexp(e2s2,iflow)
!
                      if (iflow(1) .or. iflow(2)) then
                          ifail = .true.
                          go to 2000
                      end if
!
                      qpe2 = qprob(1,j1)*qprob(2,j2)*qprob(3,j3)*e2
                      listar = listar+qpe2
!
                      if (flag /= 1) then
!
!------------------------ beta_1
                          do 102 k = 1,n1lev
                              d1(k) = d1(k) + c(1)*qpe2*e7(1,k)
  102                     end do
!
!------------------------ beta_2
                          do 105 k = n1lev+1,n1lev+n2lev
                              d1(k) = d1(k) + c(2)*qpe2*e7(2,k)
  105                     end do
!
!------------------------ beta_3
                          do 108 k = n1lev+n2lev+1,nilev
                              d1(k) = d1(k) + c(3)*qpe2*e7(3,k)
  108                     end do
!
                          do 398 k = 1,ncut(1)
                              d1(nilev+k) = d1(nilev+k) + &
                              c(1)*qpe2*e4(1,k)
  398                     end do
!
                          do 399 k = 1,ncut(2)
                              d1(nic1+k) = d1(nic1+k) + &
                              c(2)*qpe2*e4(2,k)
  399                     end do
!
                          do 401 k = 1,ncut(3)
                              d1(nic2+k) = d1(nic2+k) + &
                              c(3)*qpe2*e4(3,k)
  401                     end do
!
!------------------------ sigma_1 = phi_1
                          d1(nic3+1) = d1(nic3+1) + &
                          c(1)*qpe2*qloc(1,j1)*e3(1)
!
!------------------------ sigma_2 = phi_2
                          d1(nic3+2) = d1(nic3+2) + &
                          c(2)*qpe2*qloc(2,j2)*e3(2)
!
!------------------------ sigma_3 = phi_3
                          d1(nic3+3) = d1(nic3+3) + &
                          ic*c(3)*qpe2*qloc(3,j3)*e3(3)
!
                          if (corr) then
!---------------------------- rho_12 = phi_4
                              d1(nic3+4) = d1(nic3+4) + &
                              c(2)*qpe2*qloc(1,j1)*e3(2)
!
!---------------------------- rho_13 = phi_5
                              d1(nic3+5) = d1(nic3+5) + &
                              c(3)*qpe2*qloc(1,j1)*e3(3)
!
!---------------------------- rho_23 = phi_6
                              d1(nic3+6) = d1(nic3+6) + &
                              c(3)*qpe2*qloc(2,j2)*e3(3)
                          end if
!
                      end if
!
  145             end do
!
  148         end do
!
  150     end do
!
  909     if (abs(listar) > z1 .or. abs(listar) < z2) then
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
          if (flag /= 1) then
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
      if (flag /= 1) then
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
              else if (k <= nilev+ncut(1)) then
                  write (outbuf,'(a,i3,a,g14.4,a,f14.4,a)') 'cut1(', &
                  k-nilev,') = ',dble(theta(k)),' [',score(k),']'
              else if (k <= nilev+ncut(1)+ncut(2)) then
                  write (outbuf,'(a,i3,a,g14.4,a,f14.4,a)') 'cut2(', &
                  k-nilev-ncut(1),') = ',dble(theta(k)),' [',score(k), &
                  ']'
              else if (k == nilev+ncut(1)+ncut(2)+1) then
                  write (outbuf, &
                  '(a,g14.4,a,f14.4,a)') 'scale1           = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (k == nilev+ncut(1)+ncut(2)+2) then
                  write (outbuf, &
                  '(a,g14.4,a,f14.4,a)') 'scale2           = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (corr .and. k == nest) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'rho        = ', &
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
      end subroutine lstore_acc
