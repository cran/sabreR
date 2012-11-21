!
!***********************************************************************
!
      subroutine lshtri(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                        score,hess,ncov,flag,ifail,maxcol,link,risk, &
                        corr,iter,cflag,step,scnorm,scores,n1lev,family, &
                        n2lev,norms,ind,xil,sign,iquad,aquad,mus,taus, &
                        offlag,offpos)
      use accurate_arithmetic
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character link(3),family(3)
      integer nmes,nest,ncov,maxcol,norms,nsub(2),nm(3),it(2,nsub(1)), &
              flag,nilev,risk(nmes),iter,cflag,n1lev,n2lev,ind,sign, &
              iquad,offpos(3)
      double precision y(nmes),x(nmes,maxcol),beta(nest),xll,xil, &
             score(nest),hess(ncov),quad(3,2,256),scnorm,step, &
             scores(nmes,nest),aquad(3,2,256),mus(3,nsub(1)), &
             taus(3,nsub(1))
      logical ifail,corr,offlag(3)
!-----------------------------------------------------------------------
!     Function : Calculates deviance, score vector and Hessian matrix at
!                 the current parameter estimates, for trivariate model.
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     hesinc is an array to store the terms added to the hessian for one
!     case
      character(len=80) :: outbuf
      character(len=5) :: str
      double precision theta(maxpar),thetanir,qprob(3,256),e2,r,rr, &
             e2s(3),e8(3,maxpar),b,pi,pil(nmes), &
             e10(3,maxpar*(maxpar+3)/2),meansc(maxpar),e1(3), &
             sumxb(nmes),e3(3),expme,deri(maxpar),e,phi,listar, &
             d2(maxpar*(maxpar+1)/2),d1(maxpar),qloc(3,256), &
             e7(3,maxpar),der1(maxpar),dpieta,rs,rrs,rrss,e3s(3),e1s(3), &
             e1ss(3),e8s(3,maxpar),qpe2,e2s2
      double precision hesinc(maxpar*(maxpar+1)/2),score_total(maxpar), &
             xll_total
      integer t,i,irow,l,j1,k,m,j2,j3,norm2,rowind(nsub(1)),ir,ic,nn, &
              nq,nc,i1,i2,k1,k2
      logical sflag,anyfail,myturn,iflow(2),can(3)
      double precision :: xll_first,score_first(nest),hess_first(ncov)
      integer :: outer,firstcase,lastcase
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
      can(1) = .true.
!
      if (link(1) == 'p' .or. link(1) == 'c') then
          can(1) = .false.
      end if
!
      can(2) = .true.
!
      if (link(2) == 'p' .or. link(2) == 'c') then
          can(2) = .false.
      end if
!
      can(3) = .true.
!
      if (link(3) == 'p' .or. link(3) == 'c') then
          can(3) = .false.
      end if
!
      xll = 0
      sflag = .false.
      norm2 = norms
!
      if (family(3) == 'g') then
          norm2 = norms-1
      end if
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
!$OMP SHARED(firstcase,lastcase,i1,i2,ind,this_processor,risk,norm2), &
!$OMP SHARED(num_processors,nest,ncov,nm,iquad,quad,mus,taus,can,nilev), &
!$OMP SHARED(rowind,it,nmes,family,y,n1lev,theta,x,offlag,offpos,n2lev), &
!$OMP SHARED(norms,corr,link,flag,scores,ic,z1,z2,sflag,xil,sign), &
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
!
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
          do 150 j1 = 1,nm(1)
!
              if (iquad == 0 .and. i == i1) then
                  qloc(1,j1) = quad(1,1,j1)
                  qprob(1,j1) = quad(1,2,j1)
              else if (iquad /= 0) then
!
                  if (ind == 0  .and. outer == 1) then
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
                      if (can(1)) then
                          e2s(1) = 0
                      else
                          e2s(1) = 1
                      end if
!
                      if (can(2)) then
                          e2s(2) = 0
                      else
                          e2s(2) = 1
                      end if
!
                      if (can(3)) then
                          e2s(3) = 0
                      else
                          e2s(3) = 1
                      end if
!
                      e1(1) = 0
                      e1(2) = 0
                      e1(3) = 0
                      e1s(1) = 0
                      e1s(2) = 0
                      e1s(3) = 0
                      e1ss(1) = 0
                      e1ss(2) = 0
                      e1ss(3) = 0
                      e3(1) = 0
                      e3(2) = 0
                      e3(3) = 0
                      e3s(1) = 0
                      e3s(2) = 0
                      e3s(3) = 0
                      m = 0
!
                      do 40 k = 1,nilev
                          e7(1,k) = 0
                          e7(2,k) = 0
                          e7(3,k) = 0
                          e8(1,k) = 0
                          e8(2,k) = 0
                          e8(3,k) = 0
                          e8s(1,k) = 0
                          e8s(2,k) = 0
                          e8s(3,k) = 0
!
                          do 30 l = 1,k
                              e10(1,m+l) = 0
                              e10(2,m+l) = 0
                              e10(3,m+l) = 0
   30                     end do
!
                          m = m+k
   40                 end do
!
                      irow = rowind(i) - 1
!
                      do 90 t = 1,it(1,i)
                          irow = irow+1
                          ir = risk(irow)
!
                          if (j1 == 1 .and. j2 == 1 .and. j3 == 1) &
                          then
!
                              if (.not. myturn .and. t /= it(1,i)) then
                                  go to 90
                              else if (.not. myturn) then
                                  go to 200
                              end if
!
                              if (family(ir) == 'p') then
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
                          nn = nilev+norms
!
                          if (ir == 1) then
                              b = sumxb(irow) + theta(nn+1)*qloc(1,j1)
                          else if (ir == 2 .and. .not. corr) then
                              b = sumxb(irow) + theta(nn+2)*qloc(2,j2)
                          else if (.not. corr) then
                              b = sumxb(irow) + theta(nn+3)*qloc(3,j3)
                          else if (ir == 2) then
                              b = sumxb(irow) + theta(nn+2)* &
                              (theta(nn+4)*qloc(1,j1) + &
                              sqrt(1 - theta(nn+4)**2)*qloc(2,j2))
                          else
                              b = sumxb(irow) + theta(nn+3)* &
                              (theta(nn+5)*qloc(1,j1) + (theta(nn+6) - &
                              theta(nn+4)*theta(nn+5))*qloc(2,j2)/ &
                              sqrt(1 - theta(nn+4)**2) + &
                              sqrt((1 - theta(nn+4)**2 - &
                              theta(nn+5)**2 - theta(nn+6)**2 + &
                              2*theta(nn+4)*theta(nn+5)*theta(nn+6))/ &
                              (1 - theta(nn+4)**2))*qloc(3,j3))
                          end if
!
                          if (family(ir) /= 'g' .and. &
                          link(ir) /= 'p') then
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
                          if (family(ir) == 'b') then
!
                              if (link(ir) == 'c') then
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
                                      e2s(ir) = e2s(ir)*(1-expme)
                                      r = -e
                                      rr = e
                                  else
                                      e2s(ir) = e2s(ir)*expme
                                      r = e*(1-expme)/expme
                                      rr = r*(e-expme)/expme
                                  end if
!
                              else if (link(ir) == 'p') then
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
                                      e2s(ir) = e2s(ir)*(1-phi)
!
                                      if (phi /= 1.0) then
                                          r = -dpieta/(1-phi)
                                      else
                                          r = 1
                                      end if
!
                                  else
                                      e2s(ir) = e2s(ir)*phi
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
                                      e2s(ir) = e2s(ir) - log(1+e)
                                      r = -e/(1+e)
                                  else
                                      e2s(ir) = &
                                      e2s(ir) + log(e) - log(1+e)
                                      r = 1/(1+e)
                                  end if
!
                                  rr = e/(1+e)**2
                              end if
!
                          else if (family(ir) == 'p') then
                              e2s(ir) = e2s(ir) + &
                              y(irow)*b - e - pil(irow)
                              r = y(irow) - e
                              rr = e
                          else
!
                              if (norms == 1 .or. (norms == 2 .and. &
                              norm2 == 1 .and. ir /= 3)) then
                                  thetanir = theta(nilev+1)
                              else if (norms == 2 .and. norm2 == 1) &
                              then
                                  thetanir = theta(nilev+2)
                              else
                                  thetanir = theta(nilev+ir)
                              end if
!
                              e2s(ir) = e2s(ir) - ((y(irow) - b)**2/ &
                              (thetanir**2) + log(2*pi*thetanir**2))/2
                              r = (y(irow) - b)/thetanir**2
                              rs = (r*(y(irow) - b) - 1)/thetanir
                              rr = 1/thetanir**2
                              rrs = 2*r/thetanir
                              rrss = rr*(3*r*(y(irow) - b) - 1)
                          end if
!
                          if (flag /= 1) then
                              e3(ir) = e3(ir) + r
!
                              if (family(ir) == 'g') then
                                  e3s(ir) = e3s(ir) + rs
                              end if
!
                              if (flag == 3) then
                                  e1(ir) = e1(ir) + rr
!
                                  if (family(ir) == 'g') then
                                      e1s(ir) = e1s(ir) + rrs
                                      e1ss(ir) = e1ss(ir) + rrss
                                  end if
!
                              end if
!
                              m = 0
!
                              do 70 k = 1,nilev
                                  e7(ir,k) = e7(ir,k) + x(irow,k)*r
                                  scores(i,k) = x(irow,k)*r
!
                                  if (flag == 3) then
                                      e8(ir,k) = e8(ir,k) + x(irow,k)*rr
!
                                      if (family(ir) == 'g') then
                                          e8s(ir,k) = &
                                          e8s(ir,k) + x(irow,k)*rrs
                                      end if
!
                                      do 60 l = 1,k
                                          e10(ir,m+l) = e10(ir,m+l) + &
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
                      if (can(1) .and. (.not. can(2) .or. .not. can(3))) &
                      then
                          e2s(1) = fexp(e2s(1),iflow)
!
                          if (iflow(1) .or. iflow(2)) then
                              ifail = .true.
                              go to 2000
                          end if
!
                      end if
!
                      if (can(2) .and. (.not. can(1) .or. .not. can(3))) &
                      then
                          e2s(2) = fexp(e2s(2),iflow)
!
                          if (iflow(1) .or. iflow(2)) then
                              ifail = .true.
                              go to 2000
                          end if
!
                      end if
!
                      if ((.not. can(1) .or. .not. can(2)) .and. can(3)) &
                      then
                          e2s(3) = fexp(e2s(3),iflow)
!
                          if (iflow(1) .or. iflow(2)) then
                              ifail = .true.
                              go to 2000
                          end if
!
                      end if
!
                      if (can(1) .and. can(2) .and. can(3)) then
                          e2s2 = e2s(1) + e2s(2) + e2s(3)
                          e2 = fexp(e2s2,iflow)
!
                          if (iflow(1) .or. iflow(2)) then
                              ifail = .true.
                              go to 2000
                          end if
!
                      else
                          e2 = e2s(1)*e2s(2)*e2s(3)
                      end if
!
                      qpe2 = qprob(1,j1)*qprob(2,j2)*qprob(3,j3)*e2
                      listar = listar+qpe2
!
                      if (flag /= 1) then
!
!------------------------ beta_1
                          do 100 k = 1,n1lev
                              d1(k) = d1(k) + qpe2*e7(1,k)
  100                     end do
!
!------------------------ beta_2
                          do 105 k = n1lev+1,n1lev+n2lev
                              d1(k) = d1(k) + qpe2*e7(2,k)
  105                     end do
!
!------------------------ beta_3
                          do 107 k = n1lev+n2lev+1,nilev
                              d1(k) = d1(k) + qpe2*e7(3,k)
  107                     end do
!
!------------------------ sigma_e1
                          if (family(1) == 'g') then
                              d1(nilev+1) = d1(nilev+1) + qpe2*e3s(1)
                          end if
!
!------------------------ sigma_e2
                          if (family(2) == 'g') then
                              d1(nilev+norm2) = d1(nilev+norm2) + &
                              qpe2*e3s(2)
                          end if
!
!------------------------ sigma_e3
                          if (family(3) == 'g') then
                              d1(nilev+norms) = d1(nilev+norms) + &
                              qpe2*e3s(3)
                          end if
!
!------------------------ sigma_1 = phi_1
                          d1(nilev+norms+1) = d1(nilev+norms+1) + &
                          qpe2*qloc(1,j1)*e3(1)
!
!------------------------ sigma_2 = phi_2
                          d1(nilev+norms+2) = d1(nilev+norms+2) + &
                          qpe2*qloc(2,j2)*e3(2)
!
!------------------------ sigma_3 = phi_3
                          d1(nilev+norms+3) = d1(nilev+norms+3) + &
                          ic*qpe2*qloc(3,j3)*e3(3)
!
                          if (corr) then
!---------------------------- rho_12 = phi_4
                              d1(nilev+norms+4) = d1(nilev+norms+4) + &
                              qpe2*qloc(1,j1)*e3(2)
!
!---------------------------- rho_13 = phi_5
                              d1(nilev+norms+5) = d1(nilev+norms+5) + &
                              qpe2*qloc(1,j1)*e3(3)
!
!---------------------------- rho_23 = phi_6
                              d1(nilev+norms+6) = d1(nilev+norms+6) + &
                              qpe2*qloc(2,j2)*e3(3)
                          end if
!
                      end if
!
                      if (flag == 3) then
                          m = 0
!
!------------------------ beta_1.beta_1
                          do 120 k = 1,n1lev
!
                              do 110 l = 1,k
                                  d2(m+l) = d2(m+l) + &
                                  qpe2*(e7(1,l)*e7(1,k) - e10(1,m+l))
  110                         end do
!
                              m = m+k
  120                     end do
!
!------------------------ beta_1.beta_2
                          do 125 k = n1lev+1,n1lev+n2lev
!
                              do 115 l = 1,n1lev
                                  d2(m+l) = d2(m+l) + &
                                  qpe2*e7(1,l)*e7(2,k)
  115                         end do
!
                              m = m+k
  125                     end do
!
                          m = n1lev*(n1lev+1)/2
!
!------------------------ beta_2.beta_2
                          do 128 k = n1lev+1,n1lev+n2lev
!
                              do 118 l = n1lev+1,k
                                  d2(m+l) = d2(m+l) + &
                                  qpe2*(e7(2,l)*e7(2,k) - e10(2,m+l))
  118                         end do
!
                              m = m+k
  128                     end do
!
!------------------------ beta_1.beta_3
                          do 3125 k = n1lev+n2lev+1,nilev
!
                              do 3115 l = 1,n1lev
                                  d2(m+l) = d2(m+l) + &
                                  qpe2*e7(1,l)*e7(3,k)
 3115                         end do
!
                              m = m+k
 3125                     end do
!
                          m = (n1lev+n2lev)*(n1lev+n2lev+1)/2
!
!------------------------ beta_2.beta_3
                          do 3140 k = n1lev+n2lev+1,nilev
!
                              do 3130 l = n1lev+1,n1lev+n2lev
                                  d2(m+l) = d2(m+l) + &
                                  qpe2*e7(2,l)*e7(3,k)
 3130                         end do
!
                              m = m+k
 3140                     end do
!
                          m = (n1lev+n2lev)*(n1lev+n2lev+1)/2
!
!------------------------ beta_3.beta_3
                          do 3128 k = n1lev+n2lev+1,nilev
!
                              do 3118 l = n1lev+n2lev+1,k
                                  d2(m+l) = d2(m+l) + &
                                  qpe2*(e7(3,l)*e7(3,k) - e10(3,m+l))
 3118                         end do
!
                              m = m+k
 3128                     end do
!
                          if (family(1) == 'g') then
!
!---------------------------- beta_1.sigma_e1
                              do 1141 k = 1,n1lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*(e7(1,k)*e3s(1) - e8s(1,k))
 1141                         end do
!
!---------------------------- beta_2.sigma_e1
                              do 1142 k = n1lev+1,n1lev+n2lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*e7(2,k)*e3s(1)
 1142                         end do
!
!---------------------------- beta_3.sigma_e1
                              do 3142 k = n1lev+n2lev+1,nilev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*e7(3,k)*e3s(1)
 3142                         end do
!
                              m = m+nilev
!
!---------------------------- sigma_e1.sigma_e1
                              d2(m+1) = d2(m+1) + &
                              qpe2*(e3s(1)**2 - e1ss(1))
                              m = m+1
                          end if
!
                          if (family(2) == 'g') then
!
!---------------------------- beta_1.sigma_e2
                              do 1143 k = 1,n1lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*e7(1,k)*e3s(2)
 1143                         end do
!
!---------------------------- beta_2.sigma_e2
                              do 1144 k = n1lev+1,n1lev+n2lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*(e7(2,k)*e3s(2) - e8s(2,k))
 1144                         end do
!
!---------------------------- beta_3.sigma_e2
                              do 3144 k = n1lev+n2lev+1,nilev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*e7(3,k)*e3s(2)
 3144                         end do
!
                              m = m+nilev
!
                              if (family(1) == 'g') then
!-------------------------------- sigma_e1.sigma_e2
                                  d2(m+1) = d2(m+1) + qpe2*e3s(1)*e3s(2)
                                  m = m+1
                              end if
!
!---------------------------- sigma_e2.sigma_e2
                              d2(m+1) = d2(m+1) + &
                              qpe2*(e3s(2)**2 - e1ss(2))
                              m = m+1
                          end if
!
                          if (family(3) == 'g') then
!
!---------------------------- beta_1.sigma_e3
                              do 3143 k = 1,n1lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*e7(1,k)*e3s(3)
 3143                         end do
!
!---------------------------- beta_2.sigma_e3
                              do 3244 k = n1lev+1,n1lev+n2lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*e7(2,k)*e3s(3)
 3244                         end do
!
!---------------------------- beta_3.sigma_e3
                              do 3145 k = n1lev+n2lev+1,nilev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*(e7(3,k)*e3s(3) - e8s(3,k))
 3145                         end do
!
                              m = m+nilev
!
                              if (family(1) == 'g') then
!-------------------------------- sigma_e1.sigma_e3
                                  d2(m+1) = d2(m+1) + qpe2*e3s(1)*e3s(3)
                                  m = m+1
                              end if
!
                              if (family(2) == 'g') then
!-------------------------------- sigma_e2.sigma_e3
                                  d2(m+1) = d2(m+1) + qpe2*e3s(2)*e3s(3)
                                  m = m+1
                              end if
!
!---------------------------- sigma_e3.sigma_e3
                              d2(m+1) = d2(m+1) + &
                              qpe2*(e3s(3)**2 - e1ss(3))
                              m = m+1
                          end if
!
!------------------------ beta_1.sigma_1 = beta_1.phi_1
                          do 141 k = 1,n1lev
                              d2(m+k) = d2(m+k) + qpe2*qloc(1,j1)* &
                              (e7(1,k)*e3(1) - e8(1,k))
  141                     end do
!
!------------------------ beta_2.sigma_1 = beta_2.phi_1
                          do 142 k = n1lev+1,n1lev+n2lev
                              d2(m+k) = d2(m+k) + &
                              qpe2*qloc(1,j1)*e7(2,k)*e3(1)
  142                     end do
!
!------------------------ beta_3.sigma_1 = beta_3.phi_1
                          do 3242 k = n1lev+n2lev+1,nilev
                              d2(m+k) = d2(m+k) + &
                              qpe2*qloc(1,j1)*e7(3,k)*e3(1)
 3242                     end do
!
                          m = m+nilev
!
                          if (family(1) == 'g') then
!---------------------------- sigma_e1.sigma_1 = sigma_e1.phi_1
                              d2(m+1) = d2(m+1) + qpe2*qloc(1,j1)* &
                              (e3s(1)*e3(1) - e1s(1))
                              m = m+1
                          end if
!
                          if (family(2) == 'g') then
!---------------------------- sigma_e2.sigma_1 = sigma_e2.phi_1
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)*e3s(2)*e3(1)
                              m = m+1
                          end if
!
                          if (family(3) == 'g') then
!---------------------------- sigma_e3.sigma_1 = sigma_e3.phi_1
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)*e3s(3)*e3(1)
                              m = m+1
                          end if
!
!------------------------ sigma_1.sigma_1 = phi_1.phi_1
                          d2(m+1) = d2(m+1) + &
                          qpe2*qloc(1,j1)**2*(e3(1)**2 - e1(1))
                          m = m+1
!
!------------------------ beta_1.sigma_2 = beta_1.phi_2
                          do 143 k = 1,n1lev
                              d2(m+k) = d2(m+k) + &
                              qpe2*qloc(2,j2)*e7(1,k)*e3(2)
  143                     end do
!
!------------------------ beta_2.sigma_2 = beta_2.phi_2
                          do 144 k = n1lev+1,n1lev+n2lev
                              d2(m+k) = d2(m+k) + qpe2*qloc(2,j2)* &
                              (e7(2,k)*e3(2) - e8(2,k))
  144                     end do
!
!------------------------ beta_3.sigma_2 = beta_3.phi_2
                          do 3243 k = n1lev+n2lev+1,nilev
                              d2(m+k) = d2(m+k) + &
                              qpe2*qloc(2,j2)*e7(3,k)*e3(2)
 3243                     end do
!
                          m = m+nilev
!
                          if (family(1) == 'g') then
!---------------------------- sigma_e1.sigma_2 = sigma_e1.phi_2
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(2,j2)*e3s(1)*e3(2)
                              m = m+1
                          end if
!
                          if (family(2) == 'g') then
!---------------------------- sigma_e2.sigma_2 = sigma_e2.phi_2
                              d2(m+1) = d2(m+1) + qpe2*qloc(2,j2)* &
                              (e3s(2)*e3(2) - e1s(2))
                              m = m+1
                          end if
!
                          if (family(3) == 'g') then
!---------------------------- sigma_e3.sigma_2 = sigma_e3.phi_2
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(2,j2)*e3s(3)*e3(2)
                              m = m+1
                          end if
!
!------------------------ sigma_1.sigma_2 = phi_1.phi_2
                          d2(m+1) = d2(m+1) + &
                          qpe2*qloc(1,j1)*qloc(2,j2)*e3(1)*e3(2)
                          m = m+1
!
!------------------------ sigma_2.sigma_2 = phi_2.phi_2
                          d2(m+1) = d2(m+1) + &
                          qpe2*qloc(2,j2)**2*(e3(2)**2 - e1(2))
                          m = m+1
!
!------------------------ beta_1.sigma_3 = beta_1.phi_3
                          do 3343 k = 1,n1lev
                              d2(m+k) = d2(m+k) + &
                              ic*qpe2*qloc(3,j3)*e7(1,k)*e3(3)
 3343                     end do
!
!------------------------ beta_2.sigma_3 = beta_2.phi_3
                          do 3443 k = n1lev+1,n1lev+n2lev
                              d2(m+k) = d2(m+k) + &
                              ic*qpe2*qloc(3,j3)*e7(2,k)*e3(3)
 3443                     end do
!
!------------------------ beta_3.sigma_3 = beta_3.phi_3
                          do 3344 k = n1lev+n2lev+1,nilev
                              d2(m+k) = d2(m+k) + ic*qpe2*qloc(3,j3)* &
                              (e7(3,k)*e3(3) - e8(3,k))
 3344                     end do
!
                          m = m+nilev
!
                          if (family(1) == 'g') then
!---------------------------- sigma_e1.sigma_3 = sigma_e1.phi_3
                              d2(m+1) = d2(m+1) + &
                              ic*qpe2*qloc(3,j3)*e3s(1)*e3(3)
                              m = m+1
                          end if
!
                          if (family(2) == 'g') then
!---------------------------- sigma_e2.sigma_3 = sigma_e2.phi_3
                              d2(m+1) = d2(m+1) + &
                              ic*qpe2*qloc(3,j3)*e3s(2)*e3(3)
                              m = m+1
                          end if
!
                          if (family(3) == 'g') then
!---------------------------- sigma_e3.sigma_3 = sigma_e3.phi_3
                              d2(m+1) = d2(m+1) + ic*qpe2*qloc(3,j3)* &
                              (e3s(3)*e3(3) - e1s(3))
                              m = m+1
                          end if
!
!------------------------ sigma_1.sigma_3 = phi_1.phi_3
                          d2(m+1) = d2(m+1) + &
                          ic*qpe2*qloc(1,j1)*qloc(3,j3)*e3(1)*e3(3)
                          m = m+1
!
!------------------------ sigma_2.sigma_3 = phi_2.phi_3
                          d2(m+1) = d2(m+1) + &
                          ic*qpe2*qloc(2,j2)*qloc(3,j3)*e3(2)*e3(3)
                          m = m+1
!
!------------------------ sigma_3.sigma_3 = phi_3.phi_3
                          d2(m+1) = d2(m+1) + &
                          qpe2*qloc(3,j3)**2*(e3(3)**2 - e1(3))
!
                          if (corr) then
                              m = m+1
!
!---------------------------- beta_1.rho_12 = beta_1.phi_4
                              do 3245 k = 1,n1lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*qloc(1,j1)*e7(1,k)*e3(2)
 3245                         end do
!
!---------------------------- beta_2.rho_12 = beta_2.phi_4
                              do 3146 k = n1lev+1,n1lev+n2lev
                                  d2(m+k) = d2(m+k) + qpe2*qloc(1,j1)* &
                                  (e7(2,k)*e3(2) - e8(2,k))
 3146                         end do
!
!---------------------------- beta_3.rho_12 = beta_3.phi_4
                              do 3345 k = n1lev+n2lev+1,nilev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*qloc(1,j1)*e7(3,k)*e3(2)
 3345                         end do
!
                              m = m+nilev
!
                              if (family(1) == 'g') then
!-------------------------------- sigma_e1.rho_12 = sigma_e1.phi_4
                                  d2(m+1) = d2(m+1) + &
                                  qpe2*qloc(1,j1)*e3s(1)*e3(2)
                                  m = m+1
                              end if
!
                              if (family(2) == 'g') then
!-------------------------------- sigma_e2.rho_12 = sigma_e2.phi_4
                                  d2(m+1) = d2(m+1) + qpe2*qloc(1,j1)* &
                                  (e3s(2)*e3(2) - e1s(2))
                                  m = m+1
                              end if
!
                              if (family(3) == 'g') then
!-------------------------------- sigma_e3.rho_12 = sigma_e3.phi_4
                                  d2(m+1) = d2(m+1) + &
                                  qpe2*qloc(1,j1)*e3s(3)*e3(2)
                                  m = m+1
                              end if
!
!---------------------------- sigma_1.rho_12 = phi_1.phi_4
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)**2*e3(1)*e3(2)
                              m = m+1
!
!---------------------------- sigma_2.rho_12 = phi_2.phi_4
                              d2(m+1) = d2(m+1) + qpe2* &
                              qloc(1,j1)*qloc(2,j2)*(e3(2)**2 - e1(2))
                              m = m+1
!
!---------------------------- sigma_3.rho_12 = phi_3.phi_4
                              d2(m+1) = d2(m+1) - &
                              qpe2*qloc(1,j1)*qloc(3,j3)*e3(2)*e3(3)
                              m = m+1
!
!---------------------------- rho_12.rho_12 = phi_4.phi_4
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)**2*(e3(2)**2 - e1(2))
                              m = m+1
!
!---------------------------- beta_1.rho_13 = beta_1.phi_5
                              do 3645 k = 1,n1lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*qloc(1,j1)*e7(1,k)*e3(3)
 3645                         end do
!
!---------------------------- beta_2.rho_13 = beta_2.phi_5
                              do 3646 k = n1lev+1,n1lev+n2lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*qloc(1,j1)*e7(2,k)*e3(3)
 3646                         end do
!
!---------------------------- beta_3.rho_13 = beta_3.phi_5
                              do 3745 k = n1lev+n2lev+1,nilev
                                  d2(m+k) = d2(m+k) + qpe2*qloc(1,j1)* &
                                  (e7(3,k)*e3(3) - e8(3,k))
 3745                         end do
!
                              m = m+nilev
!
                              if (family(1) == 'g') then
!-------------------------------- sigma_e1.rho_13 = sigma_e1.phi_5
                                  d2(m+1) = d2(m+1) + &
                                  qpe2*qloc(1,j1)*e3s(1)*e3(3)
                                  m = m+1
                              end if
!
                              if (family(2) == 'g') then
!-------------------------------- sigma_e2.rho_13 = sigma_e2.phi_5
                                  d2(m+1) = d2(m+1) + &
                                  qpe2*qloc(1,j1)*e3s(2)*e3(3)
                                  m = m+1
                              end if
!
                              if (family(3) == 'g') then
!-------------------------------- sigma_e3.rho_13 = sigma_e3.phi_5
                                  d2(m+1) = d2(m+1) + qpe2*qloc(1,j1)* &
                                  (e3s(3)*e3(3) - e1s(3))
                                  m = m+1
                              end if
!
!---------------------------- sigma_1.rho_13 = phi_1.phi_5
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)**2*e3(1)*e3(3)
                              m = m+1
!
!---------------------------- sigma_2.rho_13 = phi_2.phi_5
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)*qloc(2,j2)*e3(2)*e3(3)
                              m = m+1
!
!---------------------------- sigma_3.rho_13 = phi_3.phi_5
                              d2(m+1) = d2(m+1) - qpe2* &
                              qloc(1,j1)*qloc(3,j3)*(e3(3)**2 - e1(3))
                              m = m+1
!
!---------------------------- rho_12.rho_13 = phi_4.phi_5
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)**2*e3(2)*e3(3)
                              m = m+1
!
!---------------------------- rho_13.rho_13 = phi_5.phi_5
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)**2*(e3(3)**2 - e1(3))
                              m = m+1
!
!---------------------------- beta_1.rho_23 = beta_1.phi_6
                              do 3445 k = 1,n1lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*qloc(2,j2)*e7(1,k)*e3(3)
 3445                         end do
!
!---------------------------- beta_2.rho_23 = beta_2.phi_6
                              do 3246 k = n1lev+1,n1lev+n2lev
                                  d2(m+k) = d2(m+k) + &
                                  qpe2*qloc(2,j2)*e7(2,k)*e3(3)
 3246                         end do
!
!---------------------------- beta_3.rho_23 = beta_3.phi_6
                              do 3545 k = n1lev+n2lev+1,nilev
                                  d2(m+k) = d2(m+k) + qpe2*qloc(2,j2)* &
                                  (e7(3,k)*e3(3) - e8(3,k))
 3545                         end do
!
                              m = m+nilev
!
                              if (family(1) == 'g') then
!-------------------------------- sigma_e1.rho_23 = sigma_e1.phi_6
                                  d2(m+1) = d2(m+1) + &
                                  qpe2*qloc(2,j2)*e3s(1)*e3(3)
                                  m = m+1
                              end if
!
                              if (family(2) == 'g') then
!-------------------------------- sigma_e2.rho_23 = sigma_e2.phi_6
                                  d2(m+1) = d2(m+1) + &
                                  qpe2*qloc(2,j2)*e3s(2)*e3(3)
                                  m = m+1
                              end if
!
                              if (family(3) == 'g') then
!-------------------------------- sigma_e3.rho_23 = sigma_e3.phi_6
                                  d2(m+1) = d2(m+1) + qpe2*qloc(2,j2)* &
                                  (e3s(3)*e3(3) - e1s(3))
                                  m = m+1
                              end if
!
!---------------------------- sigma_1.rho_23 = phi_1.phi_6
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)*qloc(2,j2)*e3(1)*e3(3)
                              m = m+1
!
!---------------------------- sigma_2.rho_23 = phi_2.phi_6
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(2,j2)**2*e3(2)*e3(3)
                              m = m+1
!
!---------------------------- sigma_3.rho_23 = phi_3.phi_6
                              d2(m+1) = d2(m+1) - qpe2* &
                              qloc(2,j2)*qloc(3,j3)*(e3(3)**2 - e1(3))
                              m = m+1
!
!---------------------------- rho_12.rho_23 = phi_4.phi_6
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(1,j1)*qloc(2,j2)*e3(2)*e3(3)
                              m = m+1
!
!---------------------------- rho_13.rho_23 = phi_5.phi_6
                              d2(m+1) = d2(m+1) + qpe2* &
                              qloc(1,j1)*qloc(2,j2)*(e3(3)**2 - e1(3))
                              m = m+1
!
!---------------------------- rho_23.rho_23 = phi_6.phi_6
                              d2(m+1) = d2(m+1) + &
                              qpe2*qloc(2,j2)**2*(e3(3)**2 - e1(3))
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
                  hess(m+l) = hess(m+l) - &
                  nsub(1)*(deri(k) - meansc(k))*(deri(l) - meansc(l))
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
              else if (family(1) == 'g' .and. k == nilev+1) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'sigma1     = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (family(2) == 'g' .and. k == nilev+norm2) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'sigma2     = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (family(3) == 'g' .and. k == nilev+norms) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'sigma3     = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (k == nilev+norms+1) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'scale1     = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (k == nilev+norms+2) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'scale2     = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (k == nilev+norms+3) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'scale3     = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (corr .and. k == nilev+norms+4) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'rho12   = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (corr .and. k == nilev+norms+5) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'rho13   = ', &
                  dble(theta(k)),' [',score(k),']'
              else if (corr .and. k == nilev+norms+6) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') 'rho23   = ', &
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
      end subroutine lshtri
