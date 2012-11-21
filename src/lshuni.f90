
!***********************************************************************
!
      subroutine lshuni(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev, &
                        quad,nm,score,hess,ncov,flag,endind,ifail, &
                        offlag,offpos,maxcol,link,iter,cflag,step, &
                        scnorm,scores,family,inorm,ind,xil,sign,iquad, &
                        aquad,mus,taus)
      use accurate_arithmetic
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character endind,link(3),family(3)
      integer nmes,nest,ncov,maxcol,inorm,nsub(2),nm(3),it(2,nsub(1)), &
              flag,nilev,produc(nsub(1),2),iter,cflag,offpos(3),ind, &
              sign,iquad
      double precision y(nmes),x(nmes,maxcol),beta(nest),xll,xil, &
             score(nest),hess(ncov),quad(3,2,256),scnorm,step, &
             scores(nmes,nest),aquad(3,2,256),mus(3,nsub(1)), &
             taus(3,nsub(1))
      logical ifail,offlag(3)
!-----------------------------------------------------------------------
!     Function : Calculates deviance, score vector and Hessian matrix at
!                the current parameter estimates.
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     hesinc is an array to store the terms added to the hessian for one
!     case
      character(len=320) :: outbuf
      character(len=5) :: str
      double precision theta(maxpar),qprob(256),e2,r,rr,e8(maxpar),b,pi, &
             e10(maxpar*(maxpar+3)/2),meansc(maxpar),e1, &
             sumxb(nmes),e3,expme,deri(maxpar),e,diveta,phi,listar, &
             lefeta,d2(maxpar*(maxpar+1)/2),rigeta,d1(maxpar),qloc(256), &
             e7(maxpar),der1(maxpar),dpieta,e3s,rs,e1s,e1ss,rrs,rrss, &
             e8s(maxpar),qpe2,resid(nmes),sigma,sigmae,pil(nmes)
      double precision score_total(maxpar),xll_total, &
             hesinc(maxpar*(maxpar+1)/2)
      integer t,i,irow,l,j,k,m,nl,nr,rowind(nsub(1)),nq,nc,i1,i2
      logical sflag,anyfail,myturn,iflow(2),can
      double precision :: xll_first,score_first(nest),hess_first(ncov)
      integer :: outer,firstcase,lastcase
!-----------------------------------------------------------------------
!     MPI parameter definitions
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
      nl = 0
      nr = 0
!
!---- left endpoint in model
      if (endind == 'b' .or. endind == 'l') then
          nl = 1
      end if
!
!---- right endpoint in model
      if (endind == 'b' .or. endind == 'r') then
          nr = 1
      end if
!
      lefeta = 0
      rigeta = 0
!
!---- left endpoint, psi_0
      if (nl == 1) then
          lefeta = beta(nest-nr)
      end if
!
!---- right endpoint, psi_1
      if (nr == 1) then
          rigeta = beta(nest)
      end if
!
!---- calculate 1/(1 + psi_0 + psi_1)
      diveta = 1/(1+lefeta+rigeta)
!---- initialise log-likelihood and observation count
      xll = 0
      sflag = .false.
!
!---- set parameter values and initialise score vector
      do 18 k = 1,nest
          theta(k) = beta(k)
          score(k) = 0
   18 end do
!
!---- initialise Hessian matrix
      do 15 k = 1,ncov
          hess(k) = 0
   15 end do
!
      irow = 1
!
!---- loop around cases to store starting values of irow
      do 1200 i = 1,nsub(1)
!-------- Store starting irow for this case
          rowind(i) = irow
!-------- update number of observations
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
      
!
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


!$OMP PARALLEL DO IF(outer == 2), DEFAULT(FIRSTPRIVATE), &
!$OMP SHARED(firstcase,lastcase,i1,i2,ind,this_processor,lefeta,rigeta), &
!$OMP SHARED(num_processors,nest,ncov,nm,iquad,quad,mus,taus,can,nilev), &
!$OMP SHARED(rowind,it,family,y,theta,x,offlag,offpos,inorm,link,flag), &
!$OMP SHARED(scores,z1,z2,cflag,nl,nr,produc,diveta,sflag,xil,sign), &
!$OMP SHARED(qprob,qloc,deri,aquad), &
!$OMP REDUCTION(+: xll,score,hess), REDUCTION(.or.: ifail)
!

      do 200 i = firstcase,lastcase
!
!         Any ifail must jump to the end of the loop since jumping out
!         of a parallel section is not allowed in Openmp
!
!         cycle if a previous iteration has failed
!
          if (ifail) then
              go to 2000
          end if
!
!-------- all processors execute first case. thereafter they take turns
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
!-------- initialise ifail
          ifail = .false.
          iflow(1) = .false.
          iflow(2) = .false.
!
!-------- initialise the marginal likelihood
          listar = 0
!
!-------- initialise first derivatives
          do 20 k = 1,nest
              d1(k) = 0
   20     end do
!
!-------- initialise second derivatives and hessian increment array
          do 25 k = 1,ncov
              d2(k) = 0
              hesinc(k) = 0
   25     end do
!
!-------- loop around mass points


          do 150 j = 1,nm(1)
!
              if (iquad == 0 .and. i == i1) then
!---------------- quadrature locations z_j
                  qloc(j) = quad(1,1,j)
!---------------- quadrature probabilities p_j
                  qprob(j) = quad(1,2,j)
              else if (iquad /= 0) then
!
                  if (ind == 0 .and. outer == 1) then
!-------------------- quadrature locations z_j
                      aquad(1,1,j) = mus(1,i) + taus(1,i)*quad(1,1,j)
!-------------------- quadrature probabilities p_j
                      aquad(1,2,j) = taus(1,i)*quad(1,2,j)* &
                      exp((quad(1,1,j)**2 - aquad(1,1,j)**2)/2)!
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
   30             end do
!
                  m = m+k
   40         end do
!
              irow = rowind(i) - 1
!
!------------ loop around time points
              do 90 t = 1,it(1,i)
!---------------- update number of observations
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
!-------------------- calculate beta'x_it
                      do 48 k = 1,nilev
                          sumxb(irow) = sumxb(irow) + theta(k)*x(irow,k)
   48                 end do
!
!-------------------- calculate beta'x_it = beta'x_it + offset
                      if (offlag(1)) then
                          sumxb(irow) = sumxb(irow) + x(irow,offpos(1))
                      end if
!
                      resid(irow) = y(irow) - sumxb(irow)
                  end if
!
!---------------- do the kernel of the likelihood
!---------------- i.e. the bit involving the linear predictor
!---------------- calculate b = beta'x_it + scale.z_j
                  sigma = theta(nilev+inorm+1)
                  b = sumxb(irow) + sigma*qloc(j)
!
!---------------- calculate e = exp(b)
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
!---------------- binary : calculate y_it - e/(1+e) and e/(1+e)^2
!---------------- binary : update prod{t=1,...,t_i}[e^y_it/(1+e)]
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
!---------------- poisson: calculate y_it - e and e
!---------------- poisson: update sum{t=1,...,t_i}[y_it.ln(e) - e]
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
!-------------------- binary : update sum{t=1,...,t_i}[y_it - e/(1+e)]
!-------------------- poisson: update sum{t=1,...,t_i}[y_it - e]
                      e3 = e3+r
!
                      if (family(1) == 'g') then
                          e3s = e3s+rs
                      end if
!
!-------------------- binary : update sum{t=1,...,t_i}[e/(1+e)^2]
!-------------------- poisson: update sum{t=1,...,t_i}[e]
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
!------------------------ binary : update
!------------------------ sum{t=1,...,t_i}[x_it(r).[y_it - e/(1+e)]]
!------------------------ poisson: update
!------------------------ sum{t=1,...,t_i}[x_it(r).[y_it - e]]
                          e7(k) = e7(k) + x(irow,k)*r
                          scores(irow,k) = x(irow,k)*r
!
                          if (flag == 3) then
!---------------------------- binary : update
!---------------------------- sum{t=1,...,t_i}[x_it(r).e/(1+e)^2]
!---------------------------- poisson: update
!---------------------------- sum{t=1,...,t_i}[x_it(r).e]
                              e8(k) = e8(k) + x(irow,k)*rr
!
                              if (family(1) == 'g') then
                                  e8s(k) = e8s(k) + x(irow,k)*rrs
                              end if
!
!---------------------------- binary : update sum{t=1,...,t_i}
!---------------------------- [x_it(r).x_it(s).e/(1+e)^2]
!---------------------------- poisson: update sum{t=1,...,t_i}
!---------------------------- [x_it(r).x_it(s).e]
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
!------------ let e_0 = exp(beta_0'x_it + alpha_0.omega.z_j)
!------------ let e_1 = exp(beta_1'x_it + alpha_1.omega.z_j)
!------------ let q_0 = e_0^y_it/(1+e_0), q_1 = e_1^y_it/(1+e_1)
!------------ let qq = e^y_it/(1+e)
!------------ let prod_qq = prod{t=2,...,t_i}[q_0^(1-y_it-1).q_1^y_it-1]
!------------ let prod_q = prod{t=1,...,t_i}qq
!
!------------ update marginal likelihood
!------------ original model: l_i = sum{j=1,...,q}[p_j.prod_q]
              qpe2 = qprob(j)*e2
              listar = listar+qpe2
!
              if (flag /= 1) then
!
!---------------- update explanatory variable parameter marginal
!---------------- likelihood derivative
!---------------- original model: d~l_i/d~beta(r)
!---------------- = sum{j=1,...,q}[p_j.prod_q.
!---------------- sum{t=1,...,t_i}[x_it(r)[y_it - e/(1+e)]]]
                  do 100 k = 1,nilev
                      d1(k) = d1(k) + qpe2*e7(k)
  100             end do
!
                  if (family(1) == 'g') then
                      d1(nilev+1) = d1(nilev+1) + qpe2*e3s
                  end if
!
!---------------- update scale parameter marginal likelihood derivative
!---------------- original model: d~l_i/d~omega
!---------------- = sum{j=1,...,q}[p_j.z_j.prod_q.
!---------------- sum{t=1,...,t_i}[y_it - e/(1+e)]]
                  d1(nilev+inorm+1) = d1(nilev+inorm+1) + &
                  qpe2*qloc(j)*e3
              end if
!
              if (flag == 3) then
                  m = 0
!
                  do 120 k = 1,nilev
!
!-------------------- update second derivative of marginal likelihood
!-------------------- with respect to explanatory variable parameters
!-------------------- original model: d2~l_i/d~beta(r).d~beta(s)
!-------------------- = sum{j=1,...,q}[p_j.prod_q.
!-------------------- [sum{t=1,...,t_i}[x_it(r)[y_it - e/(1+e)]].
!-------------------- sum{t=1,...,t_i}[x_it(s)[y_it - e/(1+e)]] -
!-------------------- sum{t=1,...,t_i}[x_it(r).x_it(s)e/(1+e)^2]]]
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
                          d2(m+k) = d2(m+k) + &
                          qpe2*(e7(k)*e3s - e8s(k))
 1120                 end do
!
                      m = m+nilev
                      d2(m+1) = d2(m+1) + qpe2*(e3s**2 - e1ss)
                      m = m+1
                  end if
!
!---------------- update second derivative of marginal likelihood with
!---------------- respect to explanatory variable and scale parameters
!---------------- original model: d2~l_i/d~beta(r).d~omega
!---------------- = sum{j=1,...,q}[p_j.z_j.prod_q.
!---------------- [sum{t=1,...,t_i}[x_it(r)[y_it - e/(1+e)]].
!---------------- sum{t=1,...,t_i}[y_it - e/(1+e)] -
!---------------- sum{t=1,...,t_i}[x_it(r)e/(1+e)^2]]]
                  do 140 k = 1,nilev
                      d2(m+k) = d2(m+k) + &
                      qpe2*qloc(j)*(e7(k)*e3 - e8(k))
  140             end do
!
                  m = m+nilev
!
                  if (family(1) == 'g') then
                      d2(m+1) = d2(m+1) + qpe2*qloc(j)*(e3s*e3 - e1s)
                      m = m+1
                  end if
!
!---------------- update second derivative of marginal likelihood with
!---------------- respect to scale parameter squared
!---------------- original model: d2~l_i/d~omega.d~omega
!---------------- = sum{j=1,...,q}[p_j.z_j^2.prod_q.
!---------------- [[sum{t=1,...,t_i}[y_it - e/(1+e)]]^2 -
!---------------- sum{t=1,...,t_i}[e/(1+e)^2]]]
                  d2(m+1) = d2(m+1) + qpe2*qloc(j)**2*(e3**2 - e1)
!
              end if
!
  150     end do
!
!-------- calculate scaled marginal likelihood l_i**
!-------- = l_i + psi_0.prod{t=1,...,t_i}(max(1,y_it)-y_it)
!-------- + psi_1.prod{t=1,...,t_i}y_it
          listar = listar + produc(i,1)*lefeta + produc(i,2)*rigeta
!
!-------- overflow or underflow detected
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
!-------- let the marginal log-likelihood be
!-------- l_i = ln[l_i**/(1 + psi_0 + psi_1)]
!-------- let the log-likelihood be l = sum{i=1,...,n}l_i
          if (flag /= 1) then
!
!------------ calculate explanatory variable and scale
!------------ parameter marginal log-likelihood derivatives
!------------ d~l_i/d~theta_r = (d~l_i/d~theta_r)/l_i**
              do 160 k = 1,nest-nl-nr
                  der1(k) = d1(k)/listar
  160         end do
!
!------------ calculate left endpoint parameter marginal log-likelihood
!------------ derivative d~l_i/d~psi_0
!------------ = prod{t=1,...,t_i}(max(1,y_it)-y_it)/l_i**
!------------ - 1/(1 + psi_0 + psi_1)
              if (nl == 1) then
                  der1(nest-nr) = produc(i,1)/listar - diveta
              end if
!
!------------ calculate right endpoint parameter marginal
!------------ log-likelihood derivative d~l_i/d~psi_1
!------------ = prod{t=1,...,t_i}y_it/l_i** - 1/(1 + psi_0 + psi_1)
              if (nr == 1) then
                  der1(nest) = produc(i,2)/listar - diveta
              end if
!
!------------ update score vector
!------------ d~l/d~theta_r = sum{i=1,...,n}d~l_i/d~theta_r
!------------ this is now done at the end of the 200 loop
!             (for each case)
          end if
!
          if (flag == 2) then
!
!------------ form the inner product for the scores for each parameter
!------------ for the ith case (correction needed after all cases have
!------------ been processed). first form an estimate of the mean of
!------------ each parameters score from the score for the first case
!------------ first case
              if (.not. sflag) then
!
!---------------- store marginal log-likelihood derivatives with respect
!---------------- to parameters dl_1/d_beta_r, dl_1/d_gamma,
!---------------- dl_1/d_omega, dl_1/d_psi_0, dl_1/d_psi_1
                  do 170 k = 1,nest
                      deri(k) = der1(k)
  170             end do
!
                  sflag = .true.
              end if
!
              m = 0
!
!------------ update sum{i=1,...,n}[(dl_i/d_theta_r - dl_1/d_theta_r)
!------------ .(dl_i/d_theta_s - dl_1/d_theta_s)]
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
!------------ sum the hessian matrix of 2nd derivatives. these are stored
!------------ in a vector representing the lower triangle of nest*nest
!------------ parameters. the betas are first, then gamma and alpha and
!------------ then the two endpoints
              m = 0
!
!------------ update second derivative of log-likelihood with respect to
!------------ non endpoint parameters d2~l/d~theta(r).d~theta(s)
!------------ = sum{i=1,...,n}d2~l_i/d~theta(r).d~theta(s)
!------------ = sum{i=1,...,n}[[l_i**.d2~l_i/d~theta(r).d~theta(s) -
!------------ d~l_i/d~theta(r).d~l_i/d~theta(s)]/l_i**^2]
              do 196 k = 1,nest-nl-nr
!
                  do 195 l = 1,k
                      hesinc(m+l) = hesinc(m+l) + &
                      d2(m+l)/listar - der1(k)*der1(l)
  195             end do
!
                  m = m+k
  196         end do
!
              if (nl == 1) then
!
!---------------- update second derivative of log-likelihood with
!---------------- respect to non endpoint parameter and left endpoint
!---------------- parameter d2~l/d~theta(r).d~psi_0
!---------------- = sum{i=1,...,n}d2~l_i/d~theta(r).d~psi_0
!---------------- = -sum{i=1,...,n}[d~l_i/d~theta(r).
!---------------- prod{t=1,...,t_i}(max(1,y_it)-y_it)/l_i**^2]
                  do 197 k = 1,nest-1-nr
                      hesinc(m+k) = hesinc(m+k) - &
                      der1(k)*produc(i,1)/listar
  197             end do
!
                  m = m+nest-nr
!---------------- update second derivative of log-likelihood with
!---------------- respect to left endpoint parameter squared
!---------------- d~2l/d~psi_0^2
!---------------- = sum{i=1,...,n}d2~l_i/d~psi_0^2
!---------------- = sum{i=1,...,n}[1/(1 + psi_0 + psi_1)^2 -
!---------------- [prod{t=1,...,t_i}(max(1,y_it)-y_it)]^2/l_i**^2]
                  hesinc(m) = hesinc(m) + &
                  diveta**2 - (produc(i,1)/listar)**2
              end if
!
              if (nr == 1) then
!
!---------------- update second derivative of log-likelihood with
!---------------- respect to non endpoint parameter and right endpoint
!---------------- parameter d2~l/d~theta(r).d~psi_1
!---------------- = sum{i=1,...,n}d2~l_i/d~theta(r).d~psi_1
!---------------- = -sum{i=1,...,n}
!---------------- [d~l_i/d~theta(r).prod{t=1,...,t_i}(y_it)/l_i**^2]
                  do 198 k = 1,nest-1-nl
                      hesinc(m+k) = hesinc(m+k) - &
                      der1(k)*produc(i,2)/listar
  198             end do
!
!---------------- update second derivative of log-likelihood with
!---------------- respect to left and right endpoint parameters
!---------------- d2~l/d~psi_0.d~psi_1
!---------------- = sum{i=1,...,n}d2~l_i/d~psi_0.d~psi_1
!---------------- = sum{i=1,...,n}[1/(1 + psi_0 + psi_1)^2]
                  if (nl == 1) then
                      hesinc(ncov-1) = hesinc(ncov-1) + diveta**2
                  end if
!
!---------------- update second derivative of log-likelihood with
!---------------- respect to right endpoint parameter squared
!---------------- d2~l/d~psi_1^2
!---------------- = sum{i=1,...,n}d2~l_i/d~psi_1^2
!---------------- = sum{i=1,...,n}[1/(1 + psi_0 + psi_1)^2 -
!---------------- [prod{t=1,...,t_i}y_it]^2/l_i**^2]
                  hesinc(ncov) = hesinc(ncov) + &
                  diveta**2 - (produc(i,2)/listar)**2
              end if
!
          end if
!
!-------- for the first case accumulation is done by the boss only
!-------- (since all processors perform the first case)
          if (i /= i1 .or. this_processor == boss_processor .or. &
          ind /= 0) then
!
!------------ update the log-likelihood l = sum{i=1,...,n}l_i
!------------ = sum{i=1,...,n}ln[l_i**/(1 + psi_0 + psi_1)]
              if (ind == 0) then
                  xll = xll + log(listar) + log(diveta)
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
!---------------- update score vector
!---------------- d~l/d~theta_r = sum{i=1,...,n}d~l_i/d~theta_r
                  do 165 k = 1,nest
                      score(k) = score(k) + der1(k)
  165             end do
!
!---------------- add the hessian increments to the hessian
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
!---- check to see if numerical failure has occurred
!---- reduce all the ifails from each process, creating anyfail
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
!-------- set the ifail flag for all processors and then return
          ifail = .true.
          return
      end if
!
!     add in contributions to xll, score & hess from 1st case
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
!---- collect each processor's partial sums for xll, score and hess
!---- then copy accumulated xll, score and hess (if reqd.) into correct
!---- places (only if there is more than one processor)
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
!------------ reduce hessian contributions into hesinc to save space
!------------ (instead of creating a new array called, say, hess_total)
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
!-------- calculate mean score d~l/d~theta(r)
!-------- = [sum{i=1,...,n}d~l_i/d~theta(r)]/n
          do 210 k = 1,nest
              meansc(k) = score(k)/nsub(1)
  210     end do
!
          m = 0
!
          do 230 k = 1,nest
!
!------------ calculate approximate hessian matrix values
!------------ = sum{i=1,...,n}[(d~l_i/d~theta(r) - d~l_1/d~theta(r)).
!------------ (d~l_i/d~theta(s) - d~l_1/d~theta(s))] -
!------------ n.(d~l_1/d~theta(r) - d~l/d~theta(r)).
!------------ (d~l_1/d~theta(s) - d~l/d~theta(s))
!------------ = sum{i=1,...,n}[(d~l_i/d~theta(r)).(d~l_i/d~theta(s))] -
!------------ n.(d~l/d~theta(r)).(d~l/d~theta(s))
!------------ = sum{i=1,...,n}[(d~l_i/d~theta(r) - d~l/d~theta(r)).
!------------ (d~l_i/d~theta(s) - d~l/d~theta(s))]
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
!-------- hessian is estimated by minus the matrix of second derivatives
          do 240 k = 1,ncov
!------------ calculate true hessian matrix values
!------------ -d2~l/d~theta(r).d~theta(s)
              hess(k) = -hess(k)
  240     end do
!
      end if
!
      if (flag >= 2 .and. (iter == 1 .or. cflag == 2 .or. &
      cflag == 4) .and. tflag .and. nm(1) > 1) then
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
                  else if (k == nest-nl-nr) then
                  write (outbuf,'(a,g14.4,a,f14.4,a)') &
                  'scale     = ',dble(theta(k)),' [',score(k),']'
              else if (k == nest-nr) then
                  write (outbuf,'(a,g14.4)') 'endpt 0   = ', &
                  dble(theta(k))
              else if (k == nest) then
                  write (outbuf,'(a,g14.4)') 'endpt 1   = ', &
                  dble(theta(k))
              end if
!
              call wrtlit(outbuf)
 9000     end do
!
      end if
!
      if (flag >= 2 .and. (iter == 1 .or. cflag == 2 .or. &
      cflag == 4)) then
          scnorm = 0
!
          do 9500 k = 1,nest-nl-nr
              scnorm = scnorm + score(k)**2
 9500     end do
!
          scnorm = sqrt(scnorm)
!
          if (tflag .and. nm(1) > 1) then
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
      end subroutine lshuni
