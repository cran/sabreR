      module estimatesmod
!-----------------------------------------------------------------------
!     Variables to hold the table of parameter names, estimates and
!     standard errors as constructed by DISEST to pass back to a calling
!     program
!-----------------------------------------------------------------------
      double precision, allocatable :: model_estimates(:), &
             model_errors(:)
      character(len=50), allocatable :: model_names(:)
      integer :: estimate_count = 0
!
      end module estimatesmod
!
!***********************************************************************
!
      subroutine sabrem(x,beta,name,mnames,ilev,ipos,stop)
!     version: 08.03.12, 11:22
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     There are three main storage areas in the program. These are:
!
!     1. The X array/matrix.
!        -------------------
!        Storage area for all variables including the intercept, all
!        variables listed in the DATA command and any new variables or
!        factors formed via the FACTOR and TRANSFORM commands. Also used
!        as workspace for inverting matrices. Note that X is treated as
!        both a 1-dimensional array of length MXX and a 2-dimensional
!        matrix of size NMES*MAXCOL.
!        1               ... NMES      = intercept
!        1+NMES          ... NMES*2    = 1st variable or factor level
!                        ...
!        1+NMES*(NCOL-1) ... NMES*NCOL = (NCOL-1)st "  "      "     "
!        1+NMES*NCOL     ... MXX       = workspace for matrix inversion
!
!     2. The Y array.
!        ------------
!        Storage area for the response variable y_it and related
!        products and the number of observations per case. Note that the
!        y-variate is also still stored in X.
!        1             ... NMES        = y-variate
!        logistic-normal/log-linear-normal mixture model:
!        1+NMES        ... NMES+NSUB   = number of observations T_i in
!                                        ith case
!        1+NMES+NSUB   ... NMES+NSUB*2 = product{t=1,...,T_i}(1-y_it)
!        1+NMES+NSUB*2 ... NMES+NSUB*3 = product{t=1,...,T_i}y_it
!
!     3. The BETA array.
!        ---------------
!        Storage area for the parameter estimates, aliasing indicators,
!        score vector, Hessian matrix and covariances.
!        1                     ... MAXPAR              = parameter
!                                                        estimate array
!        1+MAXPAR              ... MAXPAR*2            = intrinsic
!                                                        aliasing
!                                                        indicator (-1)
!        1+MAXPAR*2            ... MAXPAR*3            = extrinsic
!                                                        aliasing
!                                                        indicator (+1)
!        1+MAXPAR*3            ... MAXPAR*(MAXPAR+7)/2 = variance-
!                                                        covariance
!                                                        matrix
!        logistic-normal/log-linear-normal mixture model:
!        1+MAXPAR*(MAXPAR+7)/2 ... MAXPAR*(MAXPAR+9)/2 = score vector
!        1+MAXPAR*(MAXPAR+9)/2 ... MAXPAR*(MAXPAR+5)   = Hessian matrix
!
!
!     Subroutine parameters.
!     ----------------------
!     ALP       : orthogonality constant
!     ARITH     : arithmetic; f=fast, a=accurate (mantissa-exponent)
!     CMMNDS(.) : list of valid SABRE commands (3 lower case characters)
!     CNAME     : name of the case variable
!     CON       : convergence criterion
!     COV(.)    : lower triangular section of variance-covariance matrix
!     ENDIND    : endpoints; b=both, l=left, r=right, n=none
!     ERROR     : if true, an error has occurred in called subroutine
!     EST0      : initial estimate of the free left endpoint
!     EST1      : initial estimate of the free right endpoint
!     FLAG      : 1 = deviance, 2 = Meilijson Hessian, 3 = true Hessian
!     GAMMA     : estimates of the dummy variables in fixed effects
!                 model
!     GAMMSE    : standard errors of dummy variables in fixed effects
!                 model
!     HESS(.)   : estimated Hessian matrix of second derivatives
!     ICVAR     : if true, the case variable has been set
!     IDATA     : if true, a valid variable list has been defined
!     IEND(.)   : 0 = endpoint fixed at zero, 1 = endpoint free
!     IENDDF    : endpoint degrees of freedom
!     IFAIL     : if true, the fitting algorithm has failed
!     ILEV(.)   : number of levels of the variables indexed by NAME
!     ILFIT     : 0 = mixture, 1 = standard logistic/log-linear
!     ILPREV    : value of ILFIT from the previous model fit
!     IPOS(.)   : X-column positions of the variables indexed by NAME
!     IREAD     : if true, the data has been successfully read in
!     IRET      : if true, a new command line has already been read in
!     IREDF    : scale parameter(s) degrees of freedom
!     IT(.)     : number of observations per case
!     ITOTDF    : residual degrees of freedom
!     IYPOS     : X-column position of the y-variate
!     IYVAR     : if true, the y-variate has been set
!     LINE(.)   : input line converted to lower case
!     LINK      : link function; l=logit, c=complementary log-log,
!                 p=probit
!     MAXCOL    : maximum number of columns in the X array
!     MNAMES(.) : variable names in the model
!     ILEVDF    : explanatory variables degrees of freedom
!     MODE      : 0 = standard fit
!     NAME(.)   : variable names
!     NCOL      : number of columns in the X array
!     NCOV      : size of lower triangular variance-covariance matrix
!     NDUM      : number of dummy variables in fixed effects model
!     NEST      : number of parameter estimates
!     NI        : number of variables in the current model
!     NILEV     : number of X-columns in the current model
!     NITER     : maximum number of iterations
!     NM        : number of mass points used in the Gaussian quadrature
!     NMEIL     : number of approximate (i.e. Meilijson) iterations
!     NMES      : number of observations in dataset
!     NPAR      : number of model parameters prior to intrinsic aliasing
!     NSUB      : number of cases in dataset
!     NUMCOM    : command number as indexed by CMMNDS
!     NVAR      : number of variables in the X array
!     OFFLAG    : if true, the linear predictor should include an offset
!     OFFNME    : name of the offset variable
!     PRODUC(.) : products over observations per case of 1-y_it & y_it
!     QUAD(.)   : mass point locations and probabilities
!     ROWIND(.) : starting indices in the X and Y arrays for each case
!     SCA       : initial estimate of the scale parameter(s)
!     SCORES(.) : estimated score vector of first derivatives
!     TOL       : extrinsic aliasing & positive semi-definite tolerance
!     WORK(.)   : tail of the X array used as workspace
!     XLL       : log-likelihood for current model
!     YNAME     : name of the y-variate
!
!     The relationship between NVAR:NCOL and NI:NILEV is that the latter
!     member of each pair is the sum of the number of levels of each
!     variate/factor in the former member (does anyone understand?).
!     Thus for the second pair, if in a model we have 2 variates, one
!     with 1 level and one with 6 levels (i.e. a factor), NI=2, NILEV=7.
!-----------------------------------------------------------------------
!     Common block ACCMAC (initialised in subroutine BEGIN in macdep.f)
!-----------------------------------------------------------------------
!     Z1  : largest double precision real number allowed on the machine
!     Z2  : smallest double precision real number allowed on the machine
!     ZL1 : natural logarithm of Z1
!     ZL2 : natural logarithm of Z2
!     IZ  : largest integer number allowed on the machine
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
!     Common block CHANS (initialised in subroutine BEGIN in macdep.f)
!-----------------------------------------------------------------------
!     OUTCH  : output channel number
!     INCH   : terminal input channel number
!     CINCH  : current input channel
!     RINCH  : channel used by the read command
!     INPCH  : channel used by the input command
!     OOUTCH : file output channel number
!     tOUTCH : trace file output channel number
!     routch : residuals file output channel number
!     qOUTCH : time file output channel number
!     tFLAG  : trace file output flag
!     rFLAG  : resid file output flag
!     qFLAG  : time file output flag
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!     Include file limitsf90.h defines the following parameters:
!-----------------------------------------------------------------------
!     MXX    : size of the data array
!     MAXY   : maximum number of observations
!     MAXVAR : maximum number of variables
!     MAXPAR : maximum number of parameters prior to intrinsic aliasing
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
!     Data array and workspace. Dimension = MXX
!-----------------------------------------------------------------------
      double precision x(mxx)
!-----------------------------------------------------------------------
!     Response variable and related structures. Dimension = MAXY*4.
!-----------------------------------------------------------------------
      double precision y(maxy*4)
      integer risk(maxy)
!-----------------------------------------------------------------------
!     Ordinal response variable. Dimension = MAXY.
!-----------------------------------------------------------------------
      integer ord(maxy)
!-----------------------------------------------------------------------
!     Parameter related structures. Dimension = MAXPAR*(MAXPAR+5).
!-----------------------------------------------------------------------
      double precision beta(maxpar*(maxpar+5))
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character endind,arith,cmmnds(100)*6,link(3),cquad, &
                family(3),modelt
      character(len=50) :: yname,rname,cname(2),name(maxvar), &
                mnames(maxvar),offnme(3),var_prefix(3),iname(3)
      double precision sca(3),con,xll,alp,tol,est0,est1,xinit(1000), &
             xcut(1000),sig_e(3),sig(3),rho(3),pscale
      integer ipos(maxvar),ilev(maxvar),ilprev,ilfit,maxcol,npar,ndum, &
              numcom,nm(3),inormdf,nmes,nsub(2),nvar,ncol,ni,nilev, &
              iend(2),ilevdf,iredf,ienddf,nmeil,niter,itotdf,icutdf, &
              mode,numcat,ncat(3),n1var,n2var,itemp(maxvar),nargs,nlevel
      logical idata,iread,iyvar,bivar,icvar(2),iret,ifail,order, &
              corr,robust,sigflag(3),trivar,irvar,contin,error, &
              offlag(3),inflag,first,stop,cutflag,univar, &
              scaflag(3),rhoflag(3),eofflag,depend,eqscale,dfirst
      double precision gamma(maxy),gammse(maxy)
!-----------------------------------------------------------------------
!---- For comparison purposes argument prefixes must be in lower case
!---- since user input is converted to lower case.
      data var_prefix /'first','second','third'/
      data first /.true./
      save
!-----------------------------------------------------------------------
      stop = .false.
      pscale = 0.5
!
      if (first) then
!
!-------- define machine specific settings
!---------**********
          call begin
!---------**********
!
!-------- initialise the screen
!---------***********
          call scrset(0)
!---------***********
!
!-------- Initialise line buffer for log file output
!---------***********
          call filini
!---------***********
!
!-------- define valid SABRE commands
!---------***********
          call comlst(cmmnds)
!---------***********
!
!-------- define default settings
!---------***********
          call defalt(alp,nmeil,arith,icvar,nsub,con,idata,endind,est0, &
                      est1,ilfit,mode,nm,niter,iread,sca,tol,iyvar, &
                      bivar,offlag,inflag,link,corr,robust,order, &
                      cutflag,sig,sigflag,n1var,family,trivar,irvar, &
                      n2var,univar,rho,depend,eqscale,scaflag,rhoflag, &
                      dfirst,cquad)
!---------***********
!
          if (.not. terse) then
              call newlns
              call wrtscr('                         '// &
                          'Welcome to SABRE (version 6.0)')
              call newlns
          end if

          first = .false.
          tflag = .false.
          rflag = .false.
          qflag = .false.
      end if
!
!---- process new command line
      iret = .false.
!
   15 contin = .false.
!
!---- get next input line
!-----***********
      call getlin(line,contin,iret,error,eofflag)
!-----***********
!
      if (error) then
!          stop = .true.
          go to 10
      end if
!
      if (eofflag) then
          stop = .true.
          go to 999
      end if
!
!---- get command index
!-----***********
      call comind(numcom,cmmnds,line(1:6))
!-----***********
!
      if (numcom == 0) then
          go to 10
      end if
!
      go to (1010,1020,1030,1040,  10,1050,1060,1070,1080,1090,1100, &
             1110,1120,1130,1135,1140,1160,1170,1180,1190,1202,1210, &
             1215,1220,1230,1240,1250,1255,1260,1165,1300,1400,1235, &
             1212,1500,1600,1700,1900,1930,1955,1960,1970,  10,  10, &
               10,1270,2000,2010,2100,1105,2200,2300,2400), numcom
!
!---- command is ALPHA
!-----***********
 1010 call argrel(line,alp)
!-----***********
!
      go to 10
!
!---- command is APPROXIMATE
!-----*********
 1020 call meil(line,nmeil,niter)
!-----*********
!
      go to 10
!
!---- command is ARITHMETIC
!-----***********
 1030 call arithm(line,arith)
!-----***********
!
      go to 10
!
!---- command is CASE
!-----*********
 1040 call case(line,cname,icvar,iread,name,nvar,var_prefix,2,nlevel)
!-----*********
!
      if (.not. icvar(1)) then
!          stop = .true.
      end if
!
      go to 10
!
!---- command is CONVERGENCE
!-----***********
 1050 call argrel(line,con)
!-----***********
!
      go to 10
!
!---- command is DATA
 1060 if (iread) then
          call wrtlin('    --- new analysis begins')
!
!-------- define default settings as data command implies new analysis
!---------***********
          call defalt(alp,nmeil,arith,icvar,nsub,con,idata,endind,est0, &
                      est1,ilfit,mode,nm,niter,iread,sca,tol,iyvar, &
                      bivar,offlag,inflag,link,corr,robust,order, &
                      cutflag,sig,sigflag,n1var,family,trivar,irvar, &
                      n2var,univar,rho,depend,eqscale,scaflag,rhoflag, &
                      dfirst,cquad)
!---------***********
!
      end if
!
!---- get data variable names
!-----*********
      call data(line,name,nvar,idata,iname)
!-----*********
!
!---- error in DATA command
      if (.not. idata) then
!          stop = .true.
      end if
!
      go to 10
!
!---- command is DELETE
!-----***********
 1070 call delete(line,name,nvar,iread,ilev,ipos,ncol,nmes,maxcol,x, &
                  yname,iyvar,rname,bivar,cname,icvar,irvar,trivar, &
                  nlevel)
!-----***********
!
      go to 10
!
!---- command is DISPLAY
!---- *********
 1080 call disp(line,name,yname,rname,cname,ncol,ni,beta,beta(maxpar+1), &
                beta(3*maxpar+1),ipos,nm,sca,con,alp,ilfit,npar,nmes, &
                nsub,endind,nvar,ilev,iend,iread,xll,ilevdf,iredf, &
                ienddf,mnames,arith,tol,est0,est1,nmeil,niter,itotdf, &
                offlag,offnme,link,corr,bivar,robust,order,ncat, &
                sig,sig_e,inormdf,family,trivar,univar,n1var,n2var,rho, &
                nlevel,ifail,depend,eqscale,cquad,icutdf)
!-----*********
!
      go to 10
!
!---- command is ENDPOINTS
!-----***********
 1090 call endpts(line,endind,est0,est1,family)
!-----***********
!
      go to 10
!
!---- command is FACTOR
!-----***********
 1100 call factor(line,x,maxcol,nmes,ipos,name,ilev,iread,ncol,nvar, &
                  ifail)
!-----***********
!
      if (ifail) then
          stop = .true.
      end if
!
      go to 10
!
!---- command is FEFIT
 1105 continue
      ilprev = ilfit
      ilfit = 3
      univar = .true.
      go to 1115
!
!---- command is FIT
 1110 ilprev = ilfit
      ilfit = 0
!
      if (family(1) == 'p' .and. .not. scaflag(1)) then
          sca(1) = pscale
      end if
!
      if (family(2) == 'p' .and. .not. scaflag(2)) then
          sca(2) = pscale
      end if
!
      if (family(3) == 'p' .and. .not. scaflag(3)) then
          sca(3) = pscale
      end if
!
      if (family(1) == 'p' .and. nlevel == 2 .and. &
      .not. scaflag(2)) then
          sca(2) = pscale
      end if
!
      if (family(1) == 'p' .and. depend .and. .not. scaflag(2)) then
          sca(2) = pscale
      end if 
!
 1115 if (order .and. (.not. univar .or. depend)) then
          dfirst = .true.
      end if
!
!---- fit a logistic-normal/log-linear-normal mixture model
!-----**********
      call model(x,line,name,ncol,mnames,ni,iyvar,yname,bivar,rname, &
                 cname,nmes,ipos,x(nmes*ncol+1),maxcol,nm,beta, &
                 beta(maxpar+1),beta(3*maxpar+1),sca,ilfit,ilprev,con, &
                 alp,npar,nsub,endind,ilev,nvar,nilev,iend,icvar,xll, &
                 ilevdf,iredf,ienddf,tol,est0,est1,nmeil,niter,arith, &
                 itotdf,mode,y,ifail,offlag,offnme,inflag,xinit,link, &
                 risk,corr,robust,order,ncat,cutflag,xcut,sig_e,sig, &
                 sigflag,inormdf,n1var,family,trivar,n2var,univar, &
                 iname,rho,nlevel,irvar,depend,eqscale,dfirst,ndum, &
                 gamma,gammse,cquad,icutdf)
!-----**********
!
      if (ifail) then
!          go to 1240
          go to 10
      end if
!
      go to 10
!
!---- command is HELP or INFO
!-----*********
 1120 call help(numcom)
!-----*********
!
      go to 10
!
!---- command is HISTOGRAM
!-----*********
 1130 call hist(line,x,name,nvar,nmes,maxcol,iread,ilev,ipos,ncol)
!-----*********
!
      go to 10
!
!---- command is INITIAL
!-----**********
 1135 call inits(line,xinit)
!-----**********
!
      inflag = .true.
!
      go to 10
!
!---- command is INPUT
!-----**********
 1140 call input(line)
!-----**********
!
      go to 10
!
!---- command is LFIT
 1160 ilprev = ilfit   
      ilfit = 1
!
      if (order .or. offlag(1) .or. offlag(2) .or. offlag(3)) then
          ilfit = 2
      end if
!
!---- fit a standard logistic/log-linear (GLIM-type) model
      go to 1115
!
!---- command is LINK
!-----**********
 1165 call linkf(line,link(1),var_prefix,3)
!-----**********
!
      go to 10
!
!---- command is LOOK
!-----*********
 1170 call look(line,name,nvar,iread,ilev,ipos,nmes,maxcol,x,yname, &
                rname,cname,nlevel)
!-----*********
!
      go to 10
!
!---- command is MASS
!-----*********
 1180 call mass(line,nm(1),var_prefix,3)
!-----*********
!
      go to 10
!
!---- command is MAXITS
!-----*********
 1190 call iter(line,niter)
!-----*********
!
      go to 10
!
!---- command is OFFSET
!-----***********
 1202 call offset(line,offnme,offlag,iread,name,nvar,var_prefix,3)
!-----***********
!
      go to 10
!
!---- command is OUTPUT
!-----***********
 1210 call output(line)
!-----***********
!
      go to 10
!
!---- command is trace
!-----**********
 1212 call trace(line)
!-----**********
!
      go to 10
!
!---- command is resid
!-----**********
 1270 call resid(line)
!-----**********
!
      go to 10
!
!---- command is time
!-----**********
 2200 call timer(line)
!-----**********
!
      go to 10
!
!---- command is family
!-----************
 1215 call familyf(line,family(1),var_prefix,3)
!-----************
!
      if (family(1) == 'o') then
          order = .true.
      end if
!
      go to 10
!
!---- command is READ
!-----*********
 1220 call read(line,nmes,nvar,ncol,idata,iread,maxcol,ipos,x,ilev,iret, &
                numcat,yname,name,ord)
!-----*********
!
!---- error in READ command
      if (.not. iread) then
!          stop = .true.
          go to 10
      end if
!
      go to 15
!
!---- command is SCALE
!-----****************
 1230 call multi_scale(line,sca(1),var_prefix,3,scaflag)
!-----****************
!
      go to 10
!
!---- command is rho
!-----**************
 1235 call multi_rho(line,rho(1),var_prefix,3,rhoflag)
!-----**************
!
      go to 10
!
!---- command is STOP
!-----***********
 1240 call scrset(1)
!-----***********
!
      stop = .true.
!
      go to 10
!
!---- command is TOLERANCE
!-----***********
 1250 call argrel(line,tol)
!-----***********
!
      go to 10
!
!---- command is TRANSFORM
!-----**********
 1255 call trans(line,x,maxcol,nmes,ipos,name,ilev,iread,ncol,nvar, &
                 yname,rname,cname,nlevel,ifail)
!-----**********
!
      if (ifail) then
          stop = .true.
      end if
!
      go to 10
!
!---- command is YVARIATE
!-----*********
 1260 call yvar(line,yname,name,nvar,idata,iread,iyvar)
!-----*********
!
      if (.not. iyvar) then
!          stop = .true.
      end if
!
      go to 10
!
!---- command is rvariate
!-----*********
 1300 call rvar(line,rname,name,nvar,idata,iread,irvar)
!-----*********
!
      if (.not. irvar) then
!          stop = .true.
      end if
!
      go to 10
!
!---- command is correlated
!-----***********
 1400 call argcar(line,corr)
!-----***********
!
      go to 10
!
!---- command is robust
!-----***********
 1500 call argcar(line,robust)
!-----***********
!
      go to 10
!
!---- command is ordered
!-----***********
 1600 call argcar(line,order)
!-----***********
!
      go to 10
!
!---- command is cutpoints
!-----**********
 1700 call inits(line,xcut)
!-----**********
!
      cutflag = .true.
!
      go to 10
!
!---- command is sigma
!-----****************
 1900 call multi_sigma(line,sig(1),var_prefix,3,sigflag)
!-----****************
!
      go to 10
!
!---- command is nvar
!-----**********
 1930 call numlp(line,itemp,var_prefix,2)
!-----**********
!
      n1var = itemp(1)
      n2var = itemp(2)
!
      go to 10
!
!---- command is model
!-----***********
 1955 call modelf(line,modelt)
!-----***********
!
      univar = .true.
      bivar = .false.
      trivar = .false.
!
      if (modelt == 'b') then
          univar = .false.
          bivar = .true.
      else if (modelt == 't') then
          univar = .false.
          trivar = .true.
      end if
!
      go to 10
!
!---- Command is CONSTANT
!-----**************
 1960 call multi_var(line,iname,var_prefix,3,name,nvar,iread,nargs)
!-----**************
!
      go to 10
!
!---- define default settings
!-----***********
 1970 call defalt(alp,nmeil,arith,icvar,nsub,con,idata,endind,est0,est1, &
                  ilfit,mode,nm,niter,iread,sca,tol,iyvar,bivar, &
                  offlag,inflag,link,corr,robust,order,cutflag,sig, &
                  sigflag,n1var,family,trivar,irvar,n2var,univar,rho, &
                  depend,eqscale,scaflag,rhoflag,dfirst,cquad)
!-----***********
!
      go to 10
!
!---- command is DEPEND
!-----***********
 2000 call argcar(line,depend)
!-----***********
!
      go to 10
!
!---- command is EQSCALE
!-----***********
 2010 call argcar(line,eqscale)
!-----***********
!
      go to 10
!
!---- command is DER1
!-----***********
 2100 call argcar(line,dfirst)
!-----***********
!
      go to 10
!
!     command is QUADRATURE
!-----**********
 2300 call quads(line,cquad)
!-----**********
!
      go to 10
!
!     command is RESTART
!---- define default initial settings
!-----***********
 2400 call defini(sca,sig,sigflag,n1var,n2var,rho,scaflag,rhoflag, &
                  inflag,family,nlevel,pscale,depend)
!-----***********
!
      go to 10
!
   10 continue
!
  999 return
!
      end subroutine sabrem
!
!***********************************************************************
!
      subroutine defalt(alp,nmeil,arith,icvar,nsub,con,idata,endind, &
                        est0,est1,ilfit,mode,nm,niter,iread,sca,tol, &
                        iyvar,bivar,offlag,inflag,link,corr, &
                        robust,order,cutflag,sig,sigflag,n1var,family, &
                        trivar,irvar,n2var,univar,rho,depend,eqscale, &
                        scaflag,rhoflag,dfirst,cquad)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character endind,arith,link(3),family(3),cquad
      double precision con,sca(3),alp,tol,est0,est1,sig(3),rho(3)
      integer nm(3),ilfit,nmeil,niter,nsub(2),mode,n1var,n2var
      logical idata,iread,iyvar,bivar,icvar(2),trivar,corr,order, &
              robust,cutflag,sigflag(3),offlag(3),inflag,irvar, &
              univar,depend,eqscale,scaflag(3),rhoflag(3),dfirst
!-----------------------------------------------------------------------
!     Function : Sets default values for model constants and variables.
!-----------------------------------------------------------------------
!---- default settings
      alp = 1.0d-2
      nmeil = 5
      arith = 'f'
      icvar(1) = .false.
      icvar(2) = .false.
      nsub(1) = 0
      nsub(2) = 0
      con = 5.0d-5
      idata = .false.
      endind = 'n'
      est0 = 0
      est1 = 0
      ilfit = -1
      mode = 0
      nm(1) = 12
      nm(2) = 12
      nm(3) = 12
      niter = 100
      iread = .false.
      sca(1) = 1.0
      sca(2) = 1.0
      sca(3) = 1.0
      rho(1) = 0.0
      rho(2) = 0.0
      rho(3) = 0.0
      tol = 1.0d-6
      iyvar = .false.
      irvar = .false.
      bivar = .false.
      trivar = .false.
      offlag(1) = .false.
      offlag(2) = .false.
      offlag(3) = .false.
      inflag = .false.
      link(1) = 'l'
      link(2) = 'l'
      link(3) = 'l'
      corr = .true.
      robust = .false.
      order = .false.
      cutflag = .false.
      sig(1) = 1
      sig(2) = 1
      sig(3) = 1
      sigflag(1) = .false.
      sigflag(2) = .false.
      sigflag(3) = .false.
      n1var = 0
      n2var = 0
      family(1) = 'b'
      family(2) = 'b'
      family(3) = 'b'
      univar = .true.
      depend = .false.
      eqscale = .false.
      scaflag(1) = .false.
      scaflag(2) = .false.
      scaflag(3) = .false.
      rhoflag(1) = .false.
      rhoflag(2) = .false.
      rhoflag(3) = .false.
      dfirst = .false.
      cquad = 'g'
!
      return
!
      end subroutine defalt
!
!***********************************************************************
!
      subroutine defini(sca,sig,sigflag,n1var,n2var,rho,scaflag,rhoflag, &
                        inflag,family,nlevel,pscale,depend)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character family(3)
      double precision sca(3),sig(3),rho(3),pscale
      integer n1var,n2var,nlevel
      logical sigflag(3),scaflag(3),rhoflag(3),inflag,depend
!-----------------------------------------------------------------------
!     Function : Sets default values for model constants and variables.
!-----------------------------------------------------------------------
!---- default settings
      inflag = .false.
      sig(1) = 1
      sig(2) = 1
      sig(3) = 1
      sigflag(1) = .false.
      sigflag(2) = .false.
      sigflag(3) = .false.
      sca(1) = 1.0
      sca(2) = 1.0
      sca(3) = 1.0
!
      if (family(1) == 'p') then
          sca(1) = pscale
      end if
!
      if (family(2) == 'p') then
          sca(2) = pscale
      end if
!
      if (family(3) == 'p') then
          sca(3) = pscale
      end if
!
      if (family(1) == 'p' .and. nlevel == 2) then
          sca(2) = pscale
      end if
!
      if (family(1) == 'p' .and. depend) then
          sca(2) = pscale
      end if
!
      scaflag(1) = .false.
      scaflag(2) = .false.
      scaflag(3) = .false.
      rho(1) = 0.0
      rho(2) = 0.0
      rho(3) = 0.0
      rhoflag(1) = .false.
      rhoflag(2) = .false.
      rhoflag(3) = .false.
      n1var = 0
      n2var = 0
!
      return
!
      end subroutine defini
!
!***********************************************************************
!
      subroutine comlst(cmmnds)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character(len=6) :: cmmnds(100)
!-----------------------------------------------------------------------
!     Function : Stores a list of valid SABRE commands in array CMMNDS.
!-----------------------------------------------------------------------
!---- list of valid commands abbreviated to the first 6 characters
      cmmnds(1) = 'alpha '
      cmmnds(2) = 'approx'
      cmmnds(3) = 'arithm'
      cmmnds(4) = 'case  '
      cmmnds(5) = 'commen'
      cmmnds(41) = 'consta'
      cmmnds(6) = 'conver'
      cmmnds(32) = 'correl'
      cmmnds(37) = 'cutpoi'
      cmmnds(7) = 'data  '
      cmmnds(42) = 'defaul'
      cmmnds(8) = 'delete'
      cmmnds(47) = 'depend'
      cmmnds(49) = 'der1  '
      cmmnds(9) = 'displa'
      cmmnds(10) = 'endpoi'
      cmmnds(48) = 'eqscal'
      cmmnds(11) = 'factor'
      cmmnds(23) = 'family'
      cmmnds(50) = 'fefit '
      cmmnds(12) = 'fit   '
      cmmnds(13) = 'help  '
      cmmnds(14) = 'histog'
      cmmnds(15) = 'initia'
      cmmnds(16) = 'input '
      cmmnds(17) = 'lfit  '
      cmmnds(30) = 'link  '
      cmmnds(18) = 'look  '
      cmmnds(19) = 'mass  '
      cmmnds(20) = 'maximu'
      cmmnds(40) = 'model '
      cmmnds(39) = 'nvar  '
      cmmnds(21) = 'offset'
      cmmnds(36) = 'ordere'
      cmmnds(22) = 'output'
      cmmnds(52) = 'quadra'
      cmmnds(24) = 'read  '
      cmmnds(53) = 'restar'
      cmmnds(33) = 'rho   '
      cmmnds(31) = 'rvaria'
      cmmnds(25) = 'scale '
      cmmnds(38) = 'sigma '
      cmmnds(26) = 'stop  '
      cmmnds(51) = 'time  '
      cmmnds(27) = 'tolera'
      cmmnds(34) = 'trace '
      cmmnds(28) = 'transf'
      cmmnds(29) = 'yvaria'
!      cmmnds(35) = 'robust'
!      cmmnds(46) = 'resid '
!
      return
!
      end subroutine comlst
!
!***********************************************************************
!
      subroutine comind(numcom,cmmnds,command)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character(len=6) :: cmmnds(100),command
      integer numcom
!-----------------------------------------------------------------------
!     Function : Returns the command number as indexed by CMMNDS.
!-----------------------------------------------------------------------
!     SIX : first 6 characters of the input line
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=6) :: six
      character(len=3) :: three
      integer i
!-----------------------------------------------------------------------
      numcom = 0
!
!---- if no command entered then exit
      if (command == '      ') then
          return
      end if
!
!---- comment issued as 'c ' rather than the full 'com'
      if (command(1:2) == 'c ') then
          command = 'commen'
      end if
!
      six = command(1:6)
      three = command(1:3)
!
      if (command(4:4) == ' ') then
          six(5:6) = '  '
      else if (command(5:5) == ' ') then
          six(6:6) = ' '
      end if
!
      if (command(1:4) == 'con ' .or. command(1:5) == 'cons ' .or. &
      command(1:6) == 'const ') then
          six = 'consta'
      else if (command(1:4) == 'tra ' .or. command(1:5) == 'tran ' .or. &
      command(1:6) == 'trans ') then
          six = 'transf'
      end if
!
!---- look for entered characters in CMMNDS(.)
      do 10 i = 1,100
!
          if (six == cmmnds(i) .or. (three == cmmnds(i)(1:3) .and. &
          three /= 'con' .and. three /= 'tra')) then
!------------ command index
              numcom = i
              return
          end if
!
   10 end do
!
!---- if not found (i.e. NUMCOM=0)
      call wrtlin('    *** ERROR *** NO SUCH COMMAND')
!
      return
!
      end subroutine comind
!
!***********************************************************************
!
      subroutine input(line)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
!-----------------------------------------------------------------------
!     Function : Opens the input file.
!-----------------------------------------------------------------------
!     FNAME : name of the input file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer, parameter :: FILENAME_MAX = 120
      character(len=FILENAME_MAX) fname
      integer istpos,ienpos
      logical arg,error
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** NO FILE NAME GIVEN')
      else if (ienpos-istpos+1 > FILENAME_MAX) then
          call wrtlin('    *** ERROR *** '// &
                      'FILE NAME CANNOT BE MORE THAN 120 CHARACTERS')
      else
!-------- input file name
          fname = line(istpos:ienpos)
!
!-------- open input file
!---------***********
          call opnfil(inpch,fname,'old',error)
!---------***********
!
          if (error) then
              call wrtlin('    *** ERROR *** FILE DOES NOT EXIST')
              return
          else
!------------ set current input channel
              cinch = inpch
          end if
!
      end if
!
      return
!
      end subroutine input
!
!***********************************************************************
!
      subroutine output(line)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
!-----------------------------------------------------------------------
!     Function : Opens the log file.
!-----------------------------------------------------------------------
!     FNAME : name of the log file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer, parameter :: FILENAME_MAX = 120
      character(len=FILENAME_MAX) fname
      integer istpos,ienpos
      logical arg,error
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
!---- no file name given so set the default filename
      if (.not. arg) then
          fname = 'sabre.log'
      else if (ienpos-istpos+1 > FILENAME_MAX) then
          call wrtlin('    *** ERROR *** '// &
                      'FILE NAME CANNOT BE MORE THAN 120 CHARACTERS')
          return
!---- log file name
      else
          fname = line(istpos:ienpos)
      end if
!
!---- open log file
!-----***********
      call opnfil(ooutch,fname,'unknown',error)
!-----***********
!
      return
!
      end subroutine output
!
!***********************************************************************
!
      subroutine trace(line)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
!-----------------------------------------------------------------------
!     Function : Opens the trace file.
!-----------------------------------------------------------------------
!     FNAME : name of the trace file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer, parameter :: FILENAME_MAX = 120
      character(len=FILENAME_MAX) fname
      integer istpos,ienpos
      logical arg,error
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
!---- no file name given so set the default filename
      if (.not. arg) then
          fname = 'sabre.trace'
      else if (ienpos-istpos+1 > FILENAME_MAX) then
          call wrtlin('    *** ERROR *** '// &
                      'FILE NAME CANNOT BE MORE THAN 120 CHARACTERS')
          return
!---- print file name
      else
          fname = line(istpos:ienpos)
      end if
!
!---- open trace file
!-----***********
      call opnfil(toutch,fname,'unknown',error)
!-----***********
!
      tflag = .true.
!
      return
!
      end subroutine trace
!
!***********************************************************************
!
      subroutine resid(line)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
!-----------------------------------------------------------------------
!     Function : Opens the resid file.
!-----------------------------------------------------------------------
!     FNAME : name of the resid file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer, parameter :: FILENAME_MAX = 120
      character(len=FILENAME_MAX) fname
      integer istpos,ienpos
      logical arg,error
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
!---- no file name given so set the default filename
      if (.not. arg) then
          fname = 'sabre.res'
      else if (ienpos-istpos+1 > FILENAME_MAX) then
          call wrtlin('    *** ERROR *** '// &
                      'FILE NAME CANNOT BE MORE THAN 120 CHARACTERS')
          return
!---- print file name
      else
          fname = line(istpos:ienpos)
      end if
!
!---- open resid file
!-----***********
      call opnfil(routch,fname,'unknown',error)
!-----***********
!
      rflag = .true.
!
      return
!
      end subroutine resid
!
!***********************************************************************
!
      subroutine timer(line)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
!-----------------------------------------------------------------------
!     Function : Opens the time file.
!-----------------------------------------------------------------------
!     FNAME : name of the time file
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer, parameter :: FILENAME_MAX = 120
      character(len=FILENAME_MAX) fname
      integer istpos,ienpos
      logical arg,error
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
!---- no file name given so set the default filename
      if (.not. arg) then
          fname = 'sabre.time'
      else if (ienpos-istpos+1 > FILENAME_MAX) then
          call wrtlin('    *** ERROR *** '// &
                      'FILE NAME CANNOT BE MORE THAN 120 CHARACTERS')
          return
!---- print file name
      else
          fname = line(istpos:ienpos)
      end if
!
!---- open time file
!-----***********
      call opnfil(qoutch,fname,'unknown',error)
!-----***********
!
      qflag = .true.
!
      return
!
      end subroutine timer
!
!***********************************************************************
!
      subroutine data(line,name,nvar,idata,iname)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: name(maxvar),iname(3)
      integer nvar
      logical idata
!-----------------------------------------------------------------------
!     Function : Reads in the list of data variable names, declares the
!                default name of the case variable to be the first name
!                in the list and gives the first entry in NAME(.) to the
!                interaction term.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer istpos,ienpos,i
      logical arg,error
!-----------------------------------------------------------------------
      do 10 i = 1,maxvar
          name(i) = ' '
 10   end do
!
!---- first variable name is intercept
      name(1) = 'cons'
      iname(1) = 'cons'
!---- initialise number of variables in dataset
      nvar = 0
      ienpos = 3
!
      do 20 i = 2,maxvar+1
!
!-------- get position of next argument
!---------*********
          call next(istpos,ienpos,ienpos+1,line,arg)
!---------*********
!
          if (.not. arg .and. i == 2) then
              call wrtlin('    *** ERROR *** NO NAMES GIVEN')
              return
          else if (.not. arg) then
!
!------------ check for repeated variable names
!-------------***********
              call repeat(name,nvar+1,error)
!-------------***********
!
              if (.not. error) then
                  idata = .true.
              end if
!
              return
          else if (i <= maxvar) then
!------------ update number of variables in dataset
              nvar = nvar+1
!------------ next variable name
              name(i) = line(istpos:ienpos)
!
!------------ validate variable name
!-------------***********
              call namchk(name(i),error)
!-------------***********
!
              if (error) then
                  return
              end if
!
          else
              call wrtlin('    *** ERROR *** TOO MANY VARIABLES')
              call wrtlin('    Check the parameter MAXVAR in macdep.f')
          end if
!
   20 end do
!
      return
!
      end subroutine data
!
!***********************************************************************
!
      subroutine namchk(vname,error)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character(len=50) :: vname
      logical error
!-----------------------------------------------------------------------
!     Function : Validates the names of all SABRE variables and factors.
!-----------------------------------------------------------------------
!     VNAME : name to be checked
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=10) :: digit
      integer i
!-----------------------------------------------------------------------
      data digit /'0123456789'/
!-----------------------------------------------------------------------
!---- initial character of variable name not allowed to be an integer
      do 10 i = 1,10
!
          if (vname(1:1) == digit(i:i)) then
              error = .true.
              call wrtlin('    *** ERROR *** '// &
                          'NAMES MUST NOT START WITH AN INTEGER')
              return
          end if
!
   10 end do
!
      if (vname == 'cons') then
          error = .true.
          call wrtlin('    *** ERROR *** "CONS" IS A PROHIBITED NAME')
          return
      end if
!
      if (vname == 'exp' .or. vname == 'log') then
          error = .true.
          call wrtlin('    *** ERROR *** '// &
                      '"EXP" AND "LOG" ARE PROHIBITED NAMES')
          return
      end if
!
!---- valid variable name
      error = .false.
!
      return
!
      end subroutine namchk
!
!***********************************************************************
!
      subroutine repeat(name,numv,error)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character(len=50) :: name(maxvar)
      integer numv
      logical error
!-----------------------------------------------------------------------
!     Function : Checks for variable or factor name repetition.
!-----------------------------------------------------------------------
!     NUMV : number of variable names to be checked for repetition
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer i,j
!-----------------------------------------------------------------------
      do 40 i = 1,numv
!
          do 30 j = 1,numv
!
!------------ repeated variable names not allowed
              if (name(j) == name(i) .and. i /= j) then
                  error = .true.
                  call wrtlin('    *** ERROR *** '// &
                              'REPEATED NAMES NOT ALLOWED')
                  call wrtlin('    Please re-enter command with ' // &
                              'valid variable list')
                  return
              end if
!
   30     end do
!
   40 end do
!
!---- no repeated variable names
      error = .false.
!
      return
!
      end subroutine repeat
!
!***********************************************************************
!
      subroutine read(line,nmes,nvar,ncol,idata,iread,maxcol,ipos,x, &
                      ilev,iret,numcat,yname,name,ord)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: yname,name(maxvar)
      double precision x(mxx)
      integer nmes,ncol,ipos(maxvar),maxcol,nvar,ilev(maxvar),numcat, &
              ord(maxy)
      logical idata,iread,iret
!-----------------------------------------------------------------------
!     Function : Reads the data and stores it in the X array.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer, parameter :: FILENAME_MAX = 120
      character(len=FILENAME_MAX) fname
      character(len=80) :: outbuf
      double precision rarg
      integer istpos,ienpos,nobs,n,k,maxobs,i,j,ipis,reinch,iypos, &
              oldnvar,clock1
      logical calpha,arg,error,eof
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer nblocks,iblock,remain,istart,ierror
!-----------------------------------------------------------------------
      if (qflag) then
!
!---------**********
          call clock('READ',clock1,1)
!---------**********
!
      end if
!
      oldnvar = nvar
      if (.not. idata) then
          call wrtlin('    *** ERROR *** '// &
                      'DATA COMMAND MUST BE GIVEN FIRST')
          return
      end if
!
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (ienpos-istpos+1 > FILENAME_MAX) then
          call wrtlin('    *** ERROR *** '// &
                      'FILE NAME CANNOT BE MORE THAN 120 CHARACTERS')
          return
      end if
!
!---- check if file name given; if not then data is read from CINCH -
!---- this must be the secondary input channel.
!---- reading data from an input file
      if (arg) then
          go to 5
      end if
!
!---- reading data from terminal input
      if (cinch == inch) then
          call wrtscr( &
          '   --- please type data items separated by spaces' )
      else
          call wrtscr('   --- data follows')
      end if
!
      reinch = cinch
!
      go to 6
!
!---- filename
    5 fname = line(istpos:ienpos)
!
!---- open data file
!-----***********
      call opnfil(rinch,fname,'old',error)
!-----***********
!
      if (error) then
          call wrtlin('    *** ERROR *** CANNOT OPEN FILE')
          go to 999
      end if
!
      reinch = rinch
!
!---- calculate maximum number of observations (measurements) allowed
!---- = total space divided by (number of variables + constant term)
    6 maxobs = mxx/(nvar+1)
!
!---- read in file such that the first MAXOBS entries are for variable
!---- 1, the next NMES for variable 2,...., for variable NVAR
!---- The first column is filled in below with 1s for the constant term
!---- reading data from input file
      if (reinch == rinch) then
!
!-------- If data read from file, then use free-format read
!-------- Error if characters present on data file
!-------- Only the BOSS reads the data
          if (this_processor == boss_processor) then
!
              do 10 i = 1,maxobs
!
!---------------- read data from file in MAXOBS rows of NVAR columns
!-----------------***********
                  call readob(x,maxobs,i,nvar,reinch,eof,error)
!-----------------***********
!
                  if (eof .or. error) then
                      go to 11
                  end if
!
   10         end do
!
          end if
!
!-------- Broadcast end of file and error flags to all processors
!-------- (omit if there is just one processor)
!
   11     continue
!
          if (num_processors > 1) then
!
!-------------**************
              call mpi_bcast(error,1,mpi_logical,boss_processor, &
                             mpi_comm_world,ierror)
!-------------**************
!
!-------------**************
              call mpi_bcast(eof,1,mpi_logical,boss_processor, &
                             mpi_comm_world,ierror)
!-------------**************
!
          end if
!
          if (error) then
              go to 900
          end if
!
          if (eof) then
              go to 36
          end if
!
!---- reading data from terminal or command file
      else
!
!-------- If data not from file, then check first character of all input
!-------- lines to see if a character is present - evidence of a command
!-------- Then simulate free-format read by using NEXT and CHAREL
          do 25 i = 1,maxobs
              j = 1
!------------ read next input line
   15         continue
!
!-------------***********
              call readln(reinch,line,eof)
!-------------***********
!
              if (eof) then
                  go to 36
              end if
!
!------------ alphabetic character read (function CALPHA is in macdep.f)
              if (calpha(line(1:1))) then
                  go to 35
              end if
!
              call wrtfil(line)
!
              if (line(SABRE_CMDLINE_MAX-1:SABRE_CMDLINE_MAX-1) /= ' ' &
              .or. line(SABRE_CMDLINE_MAX:SABRE_CMDLINE_MAX) /= ' ') &
              then
                  call wrtlin('    *** ERROR *** DATA LINE TOO LONG')
                  call wrtlin('    Use the continuation symbol [&] at ' &
                              //'the end of the line to be continued.')
                  return
              end if
!
!------------ read line from first character
              ipis = 1
!
!------------ get position of next argument
!-------------*********
   20         call next(istpos,ienpos,ipis,line,arg)
!-------------*********
!
              if (.not. arg) then
                  go to 15
              end if
!
!------------ read line from character after current data item
              ipis = ienpos+1
!
!------------ convert character to real
!-------------***********
              call charel(rarg,line(istpos:ienpos),ienpos-istpos+1, &
                          error)
!-------------***********
!
              if (error) then
                  go to 902
              end if
!
!------------ store new data item in ith element of jth block of MAXOBS
              x(j*maxobs + i) = rarg
!------------ update block
              j = j+1
!
!------------ not end of row
              if (j <= nvar) then
                  go to 20
              end if
!
   25     end do
!
      end if
!
!---- If end of file not found after reading MAXOBS measurements, ERROR
      call wrtlin('    *** ERROR *** '// &
                  'TOO MANY OBSERVATIONS IN DATA FILE')
      write (outbuf,'(a,i3,a,i7)') '    With ',nvar, &
      ' variables, the maximum number of observations is ',maxobs
      call wrtlin(outbuf)
!
      go to 999
!
!---- potential command line in data; set indicator to syntax check it
!---- later
   35 iret = .true.
!
!---- close input file
   36 if (reinch == rinch) then
!
!---------***********
          call clsfil(rinch)
!---------***********
!
      end if
!
!---- number of measurements read in for each variable
      nmes = i-1
!
!---- Broadcast nmes to all processors (omit if there is just one
!     processor)
      if (num_processors > 1) then
!
!---------**************
          call mpi_bcast(nmes,1,mpi_integer,boss_processor, &
                         mpi_comm_world,ierror)
!---------**************
!
      end if
!
      call newlin
      write (*,'(a,i7,a)') '    ',nmes,' observations in dataset'
      call wrtlin(outbuf)
      call newlin
!---- update number of variables
      nvar = nvar+1
!---- set number of columns in X-vector to number of variables
      ncol = nvar
!
!---- now store data properly in X-vector in first NCOL*NMES cells
      do 48 j = 1,ncol
!
          do 45 i = 1,nmes
!
!------------ set first column to 1's for the constant term
              if (j == 1) then
                  x(i) = 1
              else
                  x((j-1)*nmes + i) = x((j-1)*maxobs + i)
              end if
!
   45     end do
!
   48 end do
!
!---- Broadcast the data to all processors
!---- (omit if there is just one processor)
      if (num_processors > 1) then
!
!-------- Broadcast the data in blocks of BLOCKSIZE max
          nblocks = ncol*nmes/blocksize
          remain = ncol*nmes - nblocks*blocksize
!
          do 50 iblock = 1,nblocks
              istart = 1 + (iblock-1)*blocksize
!
!-------------**************
              call mpi_bcast(x(istart),blocksize,mpi_double_precision, &
                             boss_processor,mpi_comm_world,ierror)
!-------------**************
!
   50     end do
!
          if (remain > 0) then
              istart = 1 + nblocks * blocksize
!
!-------------**************
              call mpi_bcast(x(istart),remain,mpi_double_precision, &
                             boss_processor,mpi_comm_world,ierror)
!-------------**************
!
          end if
!
      end if
!
!---- initialise positions of variables and number of levels
      do 60 i = 1,nvar
          ipos(i) = i
          ilev(i) = 1
   60 end do
!
!---- indicators to show that data successfully read
      iread = .true.
      idata = .false.
!---- calculate maximum number of columns to fill data space,
!---- remembering that SHUFFL needs an extra column of space
      maxcol = mxx/nmes - 1
!
      if (qflag) then
!
!---------**********
          call clock('',clock1,2)
!---------**********
!
      end if
!
!---- normal exit
      return
!
!---- error exits for READ
  900 continue
!
      call wrtlin('    *** ERROR *** ERROR READING DATA')
      write (outbuf,'(a,i2,a,i7,a)') '    Data element in column ',j, &
      ' of row ',i,' is invalid'
      call wrtlin(outbuf)
!
      go to 999
!
  902 continue
!
      call wrtlin('    *** ERROR *** ERROR READING DATA')
      write (outbuf,'(a,i6)') '    Non-numeric characters found ['// &
      line(istpos:ienpos)//'] while reading measurement ',i
      call wrtlin(outbuf)
!---- tidy up
  999 iread = .false.
      nvar = oldnvar
!
!---- close input file
      if (reinch == rinch) then
!
!---------***********
          call clsfil(rinch)
!---------***********
!
      end if
!
      return
!
      end subroutine read
!
!***********************************************************************
!
      subroutine delete(line,name,nvar,iread,ilev,ipos,ncol,nmes,maxcol, &
                        x,yname,iyvar,rname,bivar,cname,icvar,irvar, &
                        trivar,nlevel)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: name(maxvar),yname,rname,cname(2)
      integer nmes,maxcol,nlevel,nvar,ncol,ilev(maxvar),ipos(maxvar)
      double precision x(nmes,maxcol)
      logical iread,iyvar,bivar,icvar(2),trivar,irvar
!-----------------------------------------------------------------------
!     Function : Deletes a list of variates.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=50) :: dname
      integer istpos,ienpos,i,j,ixpos,ixlev,iarg, length
      logical arg
!-----------------------------------------------------------------------
      if (.not. iread) then
          call wrtlin('    *** ERROR *** NO DATA READ')
          return
      end if
!
!---- initialise start search position
      ienpos = 3
      iarg = 0
    5 iarg = iarg+1
!
!---- get position of argument
!-----*********
      call next(istpos,ienpos,ienpos+1,line,arg)
!-----*********
!
      if (.not. arg .and. iarg == 1) then
          call wrtlin('    *** ERROR *** NO NAMES GIVEN')
          go to 999
      else if (.not. arg) then
          go to 999
      else
!-------- name of variable to be deleted
          dname = line(istpos:ienpos)
      end if
!
!---- check if variate exists
      do 10 i = 1,nvar
!
!-------- variable name OK
          if (name(i) == dname) then
              go to 15
          end if
!
   10 end do
!
      call wrtlin('    *** ERROR *** VARIATE `'//dname(1:length(dname)) &
                  //'` DOES NOT EXIST')
!
      go to 999
!
!---- y-variable no longer included in model
   15 if (dname == yname) then
          iyvar = .false.
      end if
!
      if (dname == rname) then
          irvar = .false.
          bivar = .false.
          trivar = .false.
      end if
!
!---- case-variable no longer included in model
      if (dname == cname(1)) then
          icvar(1) = .false.
      end if
!
      if (nlevel == 2 .and. dname == cname(2)) then
          icvar(2) = .false.
      end if
!
!---- get position of deleted variable in X-vector
!-----***********
      call fndpos(nvar,dname,name,ixpos,ipos)
!-----***********
!
!---- get number of levels of deleted variable
!-----***********
      call fndlev(nvar,dname,name,ixlev,ilev)
!-----***********
!
!---- move deleted variable to end of X-vector and close gaps
!-----***********
      call shuffl(x,ixpos,ncol-ixlev+1,ncol,nmes,ipos,maxcol,ixlev,nvar)
!-----***********
!
!---- shuffle up NAME and corresponding level and position
      do 30 i = 1,nvar
!
          if (name(i) == dname) then
!
!------------ loop for each variable from deleted one onwards
              do 20 j = i,nvar-1
!---------------- update tracking vectors to account for deleted
!                 variable
                  name(j) = name(j+1)
                  ipos(j) = ipos(j+1)
                  ilev(j) = ilev(j+1)
   20         end do
!
              go to 40
!
          end if
!
   30 end do
!
      call wrtlin('    *** ERROR *** VARIATE `'//dname(1:length(dname)) &
      //'` DOES NOT EXIST')
!---- update model information for deleted variable
!---- update number of columns used in X-matrix
   40 ncol = ncol-ixlev
!---- update number of variables in model
      nvar = nvar-1
!
      go to 5
!
  999 return
!
      end subroutine delete
!
!***********************************************************************
!
      subroutine fndpos(nvar,vname,name,ixpos,ipos)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character(len=50) :: name(maxvar),vname
      integer nvar,ixpos,ipos(maxvar)
!-----------------------------------------------------------------------
!     Function : Finds the position of a variable in the X array.
!-----------------------------------------------------------------------
!     VNAME : variable name
!     IXPOS : position of the variable in the X array
!-----------------------------------------------------------------------
      integer i
!-----------------------------------------------------------------------
      do 10 i = 1,nvar
!
          if (vname == name(i)) then
!------------ position of variable in X-vector
              ixpos = ipos(i)
              return
          end if
!
   10 end do
!
      return
!
      end subroutine fndpos
!
!***********************************************************************
!
      subroutine fndlev(nvar,vname,name,ixlev,ilev)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character(len=50) :: name(maxvar),vname
      integer nvar,ixlev,ilev(maxvar)
!-----------------------------------------------------------------------
!     Function : Finds the number of levels in a variable.
!-----------------------------------------------------------------------
!     VNAME : variable name
!     IXLEV : number of levels in the variable
!-----------------------------------------------------------------------
      integer i
!-----------------------------------------------------------------------
      do 10 i = 1,nvar
!
          if (vname == name(i)) then
!------------ number of levels of variable
              ixlev = ilev(i)
              return
          end if
!
   10 end do
!
      return
!
      end subroutine fndlev
!
!***********************************************************************
!
      subroutine shuffl(x,istpos,ienpos,ncol,nmes,ipos,maxcol,nlev,nvar)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer nmes,maxcol,nvar,istpos,ienpos,ncol,ipos(maxvar),nlev
      double precision x(nmes,maxcol+1)
!-----------------------------------------------------------------------
!     Function : Moves the variate at position ISTPOS to IENPOS and then
!                reshuffles the X array to close all the gaps. For a
!                factor, the shuffling is done for NLEV levels.
!-----------------------------------------------------------------------
!     ISTPOS : X-position of variable to be moved, prior to shuffling
!     IENPOS : X-position of moved variable, on completion of shuffle
!     NLEV   : number of levels in variable being shuffled
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer i,j,k,l
!-----------------------------------------------------------------------
!---- no move necessary
      if (istpos == ienpos) then
          return
      end if
!
!---- shuffle error detected if IENPOS which indexes the beginning
!---- column of a factor is too near the end to fit all its levels in.
!---- Similarly for ISTPOS
!---- no room to move factor
      if (ienpos+nlev-1 > ncol .or. istpos+nlev-1 > ncol) then
          call wrtlin('    *** ERROR *** ERROR IN SHUFFLE ROUTINE')
          return
      end if
!
!---- move variable forward in X-vector
      if (ienpos > istpos) then
!
          do 60 i = 1,nmes
!
!------------ repeat the operation NLEV times
              do 50 l = 1,nlev
!
!---------------- shuffle columns forward from (IENPOS+NLEV) to make a
!---------------- space for one of the columns of ISTPOS (the last
!---------------- shuffle puts NCOL into NCOL+1)
!---------------- shuffling variable to end of X-vector
                  if (ienpos+nlev == ncol+1) then
                      go to 30
                  end if
!
!---------------- shuffle columns following end position forward to
!---------------- create space for move
                  do 20 j = ncol,ienpos+nlev,-1
                      x(i,j+1) = x(i,j)
   20             end do
!
!---------------- move column ISTPOS (still at ISTPOS) to IENPOS+NLEV
   30             x(i,ienpos+nlev) = x(i,istpos)
!
!---------------- shuffle columns from start position onwards backward
!---------------- to fill gap created by move
                  do 40 j = istpos,ncol
                      x(i,j) = x(i,j+1)
   40             end do
!
   50         end do
!
   60     end do
!
!-------- reset vector of variable start positions in X-vector
          do 70 k = 1,nvar
!
!------------ update position of intermediate variables in X-vector
              if (ipos(k) <= ienpos+nlev-1 .and. ipos(k) > istpos) &
              then
                  ipos(k) = ipos(k) - nlev
!------------ update position of moved variable in X-vector
              else if (ipos(k) == istpos) then
                  ipos(k) = ienpos
              end if
!
   70     end do
!
!---- move variable backward in X-vector
      else
!
          do 120 i = 1,nmes
!
              do 115 l = 0,nlev-1
!
!---------------- shuffle columns forward to make a space for L+1st
!---------------- level of ISTPOS. The last shuffle puts NCOL into
!---------------- NCOL+1
!---------------- shuffle columns from L after end position forward to
!---------------- create space for move
                  do 110 j = ncol,ienpos+l,-1
                      x(i,j+1) = x(i,j)
  110             end do
!
!---------------- move column ISTPOS (now in ISTPOS+1+L) to IENPOS+L
!---------------- move variable
                  x(i,ienpos+l) = x(i,istpos+1+l)
!
!---------------- shuffle columns from L+2 after start position backward
!---------------- to fill gap created by move
                  do 112 j = istpos+1+l,ncol
                      x(i,j) = x(i,j+1)
  112             end do
!
  115         end do
!
  120     end do
!
!-------- update positions of variables in X-vector
          do 130 k = 1,nvar
!
!------------ update position of intermediate variables in X-vector
              if (ipos(k) >= ienpos .and. ipos(k) < istpos) then
                  ipos(k) = ipos(k) + nlev
!------------ update position of moved variable in X-vector
              else if (ipos(k) == istpos) then
                  ipos(k) = ienpos
              end if
!
  130     end do
!
      end if
!
      return
!
      end subroutine shuffl
!
!***********************************************************************
!
      subroutine look(line,name,nvar,iread,ilev,ipos,nmes,maxcol,x, &
                      yname,rname,cname,nlevel)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: name(maxvar),yname,rname,cname(2)
      integer nmes,maxcol,nlevel,nvar,ilev(maxvar),ipos(maxvar)
      double precision x(nmes,maxcol)
      logical iread
!-----------------------------------------------------------------------
!     Function : Displays the values of the listed variables/factors.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      character eno,cpiece(6)*12,vname(8)*50
      integer istpos,ienpos,i,j,ixpos,ixlev,iarg,level(6),k,posit(6), &
              numarg,line1,line2,scount,index,ivalue,maxlines,length
      logical lrange,arg,error,eof
!-----------------------------------------------------------------------
      parameter (maxlines=50)
!-----------------------------------------------------------------------
      lrange = .false.
!
      if (.not. iread) then
          call wrtlin('    *** ERROR *** NO DATA READ')
          return
      end if
!
!---- initialise start search position
      ienpos = 3
      iarg = 0
!
!---- get position of next argument
!-----*********
   10 call next(istpos,ienpos,ienpos+1,line,arg)
!-----*********
!
!---- If no argument found, go to the printing out stuff. Or, if there
!---- aren't any arguments, return. If argument begins with an integer,
!---- this must be a line marker so look for next one. If just a normal
!---- variate carry on with existence checks, looking for next argument
      if (.not. arg) then
!
          if (iarg == 0) then
              return
          end if
!
          numarg = iarg
          go to 50
      else
!
!-------- check if argument is a line marker
!-------- convert character to integer
!---------***********
          call chaint(line(istpos:ienpos),line1,ienpos-istpos+1,error)
!---------***********
!
          if (.not. error) then
!
              if (iarg == 0) then
                  return
              end if
!
!------------ first line marker greater than number of measurements
              if (line1 > nmes) then
                  call wrtlin('    *** ERROR *** FIRST LINE MARKER IS ' &
                  //'LARGER THAN THE NUMBER OF OBSERVATIONS')
                  return
!------------ first line marker less than 1
              else if (line1 < 1) then
                  call wrtlin('    *** ERROR *** '// &
                  'FIRST LINE MARKER IS LESS THAN 1')
                  return
              end if
!
!------------ get position of next argument
!-------------*********
              call next(istpos,ienpos,ienpos+1,line,arg)
!-------------*********
!
              if (.not. arg) then
                  call wrtlin('    *** ERROR *** '// &
                  'YOU NEED TWO LINE MARKERS')
                  return
              end if
!
!------------ check if the argument really is an integer and store it
!             anyway
!------------ convert character to integer
!-------------***********
              call chaint(line(istpos:ienpos),line2,ienpos-istpos+1, &
                          error)
!-------------***********
!
              if (error) then
                  call wrtlin('    *** ERROR *** '// &
                  'SECOND LINE MARKER IS NOT AN INTEGER')
                  return
              end if
!
!------------ second line marker greater than number of measurements
              if (line2 > nmes) then
                  call wrtlin('    *** ERROR *** SECOND LINE MARKER IS ' &
                  //'LARGER THAN THE NUMBER OF OBSERVATIONS')
                  return
!------------ second line marker less than 1
              else if (line2 < 1) then
                  call wrtlin('    *** ERROR *** '// &
                  'SECOND LINE MARKER IS LESS THAN 1')
                  return
!------------ first line marker is larger than second
              else if (line1 > line2) then
                  call wrtlin('    *** ERROR *** '// &
                  'SECOND LINE MARKER IS SMALLER THAN FIRST')
                  return
              else if (line2-line1 > maxlines-1) then
!------------ too many lines requested
                  call wrtlin('    *** WARNING *** '// &
                  'TOO MANY LINES REQUESTED')
                  write (outbuf,'(a,i3,a)') '    Only ',maxlines, &
                  ' will be displayed'
                  call wrtlin(outbuf)
                  line2 = line1 + maxlines - 1
              end if
!
              lrange = .true.
              numarg = iarg
              go to 50
          else
              iarg = iarg+1
!
!------------ number of variables which may be looked at limited to 6
              if (iarg > 6) then
                  call wrtlin('    *** WARNING *** '// &
                  'YOU ASKED FOR TOO MANY VARIATES')
                  call wrtlin('    Only the first 6 will be displayed')
                  iarg = iarg-1
                  numarg = iarg
                  go to 50
              end if
!
              vname(iarg) = line(istpos:ienpos)
          end if
!
      end if
!
!---- check if variate exists
      do 20 i = 1,nvar
!
!-------- variate name OK
          if (name(i) == vname(iarg)) then
              go to 30
          end if
!
   20 end do
!
      if (length(vname(iarg)) <= 35) then
          call wrtlin('    *** WARNING *** VARIATE `'// &
          vname(iarg)(1:length(vname(iarg)))//'` DOES NOT EXIST')
      else
          call wrtlin('    *** WARNING *** VARIATE ' )
          call wrtlin( vname(iarg)(1:length(vname(iarg))) )
          call wrtlin( 'DOES NOT EXIST')
      end if
      iarg = iarg-1
!
      go to 10
!
!---- get position of variable in X-vector
!-----***********
   30 call fndpos(nvar,vname(iarg),name,ixpos,ipos)
!-----***********
!
!---- get number of levels of variable to be looked at
!-----***********
      call fndlev(nvar,vname(iarg),name,ixlev,ilev)
!-----***********
!
!---- Store the level for later on. This is so that we know whether we
!---- have to do any funny business to convert the 0/1 factor columns to
!---- a column of integer 'levels'. Also, store the position.
      level(iarg) = ixlev
      posit(iarg) = ixpos
!
      go to 10
!
!---- Start printing out the stuff - stop every 20 lines
!---- Make sure that factors are properly dealt with
!---- Write out the column names
   50 continue
!
      call newlin
!    Print the first 12 characters of each name
      write (outbuf,'(7a)') '          ',(vname(j)(1:12),j=1,numarg)
      call wrtlin(outbuf)
      write (outbuf,'(7a)') '        ',('____________',j=1,numarg-1), &
      '___________'
      call wrtlin(outbuf)
!
!---- Set the line dividers if they have not been specified in the
!---- command variable values for all measurements to be printed out
      if (.not. lrange) then
          line1 = 1
          line2 = min(nmes,maxlines)
      end if
!
!---- set the count for scrolling
      scount = 0
!
      do 70 i = line1,line2
          index = 0
!
          do 60 j = 1,numarg
!------------ Store the factor levels and variable values in temporary
!------------ buffers
              index = index+1
!
!------------ factor
              if (level(j) > 1) then
!
!---------------- get factor level, remembering that the factor columns
!---------------- are stored in reverse order
                  do 55 k = 1,level(j)
!
                      if (x(i,posit(j)+k-1) /= 0) then
!------------------------ factor level
                          ivalue = level(j) - k + 1
                          go to 56
                      end if
!
   55             end do
!
                  call wrtlin('    *** ERROR *** '// &
                           'MISTAKE IN FACTOR STORAGE')
                  return
   56             write (cpiece(index),'(i3,a)') ivalue,'    '
!------------ variable
              else
!
!---------------- case variable or y-variate
                  if (vname(j) == cname(1) .or. (nlevel == 2 .and. &
                  vname(j) == cname(2)) .or. vname(j) == yname .or. &
                  vname(j) == rname) then
                      write (cpiece(index),'(i6,a)') &
                      idnint(x(i,posit(j))),' '
!---------------- non-zero value
                  else if (x(i,posit(j)) /= 0) then
                      write (cpiece(index),'(g11.4)') x(i,posit(j))
                  else
                      cpiece(index) = '   0.0'
                  end if
!
              end if
!
   60     end do
!
!-------- print out variable/factor values
          write (outbuf,'(i7,7a)') i,' ',(cpiece(k),k=1,numarg)
          call wrtlin(outbuf)
!-------- Increment line count
          scount = scount+1
!
!-------- Prompt and reply for scrolling removed for compatibility with
!-------- Sabre in R where such interaction is not possible
!-------- scroll ?
!          IF (SCOUNT == 20 .AND. I < LINE2) THEN
!              call newlns
!-------------***********
!              CALL PROMPT('    continue (y/n) ? ')
!-------------***********
!
!-------------***********
!              CALL READLN(inch,eno,eof)
!-------------***********
!
!--------- quit scrolling
!              IF (ENO == 'N' .OR. ENO == 'n') go to 80
!--------- continue scrolling
!              SCOUNT = 0
!              call newlin
!              write (outbuf,'(13a)') '    ',
!     &        ('      ',VNAME(J)(1:12),J=1,NUMARG)
!              call wrtlin(outbuf)
!              write (outbuf,'(7a)') '        ',
!     &        ('____________',J=1,NUMARG-1),'___________'
!              call wrtlin(outbuf)
!          END IF
!
   70 end do
!
   80 continue
!
      call newlin
!
      return
!
      end subroutine look
!
!***********************************************************************
!
      subroutine hist(line,x,name,nvar,nmes,maxcol,iread,ilev,ipos,ncol)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: name(maxvar)
      integer nmes,maxcol,nvar,ilev(maxvar),ncol,ipos(maxvar)
      double precision x(nmes,maxcol)
      logical iread
!-----------------------------------------------------------------------
!     Function : Draws a histogram of a variable or factor.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=50) :: hname
      integer ienpos,istpos,i,j,dummy,ixpos,freq(50),levs,ixlev
      logical groups,arg,error
!-----------------------------------------------------------------------
      if (.not. iread) then
          call wrtlin('    *** ERROR *** NO DATA READ')
          return
      end if
!
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** MISSING ARGUMENT')
          return
      end if
!
!---- reduce name of histogram variable to first 50 characters
      if (ienpos-istpos+1 > 50) then
          hname = line(istpos:istpos+49)
      else
          hname = line(istpos:ienpos)
      end if
!
!---- determine whether HISTOGRAM variate exists
      do 10 i = 1,nvar
!
          if (name(i) == hname) then
              go to 20
          end if
!
   10 end do
!
      call wrtlin('    *** ERROR *** VARIATE DOES NOT EXIST')
      call wrtlin('    Use DISPLAY VARIABLES to check spelling')
!
      return
!
!---- histogram variable name OK
   20 groups = .false.
      dummy = ienpos
!
!---- get position of argument
!-----*********
      call next(istpos,ienpos,dummy,line,arg)
!-----*********
!
      if (arg) then
!
!-------- convert character into integer number of classes
!---------***********
          call chaint(line(istpos:ienpos),levs,ienpos-istpos+1,error)
!---------***********
!
          if (.not. error .and. levs > 1) then
!------------ number of classes specified
              groups = .true.
          else if (.not. error) then
              call wrtlin('    *** ERROR *** '// &
              'THERE MUST BE AT LEAST TWO HISTOGRAM BINS')
              return
          else
              call wrtlin( &
              '    *** ERROR *** ARGUMENT IS NOT AN INTEGER')
          end if
!
      end if
!
!---- get position of histogram variable in X-vector
!-----***********
      call fndpos(nvar,hname,name,ixpos,ipos)
!-----***********
!
!---- get number of levels of histogram variable
!-----***********
      call fndlev(nvar,hname,name,ixlev,ilev)
!-----***********
!
!---- factor - so cannot specify number of classes
      if (ixlev > 1 .and. groups) then
          call wrtlin('    *** WARNING *** YOU '// &
          'CANNOT SPECIFY GROUPS WHEN PLOTTING A FACTOR')
          call wrtlin('    Number of groups set to '// &
          'the number of factor levels')
          groups = .false.
      end if
!
!---- draw the nice histogram. If a factor, calculate the frequencies
!---- and ask for a histogram with IXLEV bars (remembering of course
!---- that the factor columns are stored in reverse order)
!---- variable - so raw data
      if (ixlev == 1) then
!
!-------- set default number of classes
          if (.not. groups) then
              levs = 11
          end if
!
!-------- draw histogram
!---------***********
          call hisfrq(freq,x(1,ixpos),nmes,levs,.true.)
!---------***********
!
!---- factor
      else
!
          do 30 i = 1,ixlev
              freq(i) = 0
   30     end do
!
          do 50 i = 1,nmes
!
              do 40 j = 1,ixlev
!
!---------------- update frequencies for each factor level
                  if (x(i,ixpos-1+j) /= 0) then
                      freq(ixlev+1-j) = freq(ixlev+1-j) + 1
                      go to 50
                  end if
!
   40         end do
!
   50     end do
!
!-------- draw histogram
!---------***********
          call hisfrq(freq,x(1,ncol+1),nmes,ixlev,.false.)
!---------***********
!
      end if
!
      return
!
      end subroutine hist
!
!***********************************************************************
!
      subroutine hisfrq(freq,dat,nmes,m,ind)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer nmes,freq(50),m
      double precision dat(nmes)
      logical ind
!-----------------------------------------------------------------------
!     Function : Calculates the frequencies for each histogram group.
!-----------------------------------------------------------------------
!     ALGORITHM AS45; APPL. STATIST. (1971) VOL. 20, P.332.
!     Published in the book
!     'APPLIED STATISTICS ALGORITHMS' by P. Griffiths and I.D. Hill.
!     PLUS SOME MAJOR AND VITALLY IMPORTANT CHANGES MADE BY JON BARRY IN
!     JAN 1989. INFACT, THE PROGRAM IS HARDLY RECOGNISABLE FROM ITS
!     ORIGINAL FORM. FURTHER ALTERATIONS MADE BY DAVE STOTT IN 1995/96.
!-----------------------------------------------------------------------
!     FREQ   : factor level frequencies
!     DAT    : variable data values (only used if IND is true)
!     M      : number of groups
!     IND    : true if argument is a variable, false if a factor
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      character(len=6) :: cout(11)
      double precision xmin,xmax,r,b,tk,step,c,eta,xm,out(11)
      integer scale,key,k,i,j,kount,imax,iout(11),mm,jmax
!-----------------------------------------------------------------------
!     DEFINE MAGNITUDE OF SMALLEST ACCEPTABLE NUMBER
!-----------------------------------------------------------------------
      data eta /1.0d-6/
!-----------------------------------------------------------------------
      if (ind .and. m > 11) then
          call wrtlin('    *** ERROR *** TOO MANY CLASSES')
          return
      end if
!
      key = 1
!
      if (m > 11) then
          mm = m
          key = (m+11)/11
          k = 1
!
          do 11 i = 1,11
              iout(i) = 0
!
              do 10 j = 1,key
                  iout(i) = iout(i) + freq(k)
                  k = k+1
!
                  if (k > m) then
                      go to 12
                  end if
!
   10         end do
!
   11     end do
!
   12     m = i
!
          do 14 i = 1,m
              freq(i) = iout(i)
   14     end do
!
          call newlin
          write (outbuf,'(a,i1,a)') '    FACTOR LEVELS COMBINED ',key, &
          ' AT A TIME TO FIT WIDTH'
          call wrtlin(outbuf)
      end if
!
      if (.not. ind) then
          go to 120
      end if
!
!---- define a suitable scale
   15 xmin = dat(1)
      xmax = xmin
!
      do 30 i = 2,nmes
!
!-------- update smallest data value
          if (dat(i) < xmin) then
              xmin = dat(i)
          end if
!
!-------- update largest data value
          if (dat(i) > xmax) then
              xmax = dat(i)
          end if
!
   30 end do
!
      if (xmax-xmin < eta*m) then
          call wrtlin('    *** ERROR *** ALL VALUES ARE EQUAL')
          return
      end if
!
      key = 1
      kount = 0
!---- range
   40 r = xmax-xmin
      b = xmin
!
   50 if (r > 1) then
          go to 60
      end if
!
      kount = kount+1
      r = r*10
!
      go to 50
!
   60 if (r <= 10) then
          go to 70
      end if
!
      kount = kount-1
      r = r/10
!
      go to 60
!
!---- range scaled to (1,10]
   70 if (key >= 3) then
          go to 80
      end if
!
      tk = 10d0**kount
      b = b*tk
!
      if (b < 0 .and. b /= dint(b)) then
          b = b-1
      end if
!
      b = dint(b)/tk
      r = (xmax-b)/(m-1)
      kount = 0
      key = key+2
!
      go to 50
!
   80 step = dint(r)
!
      if (step /= r) then
          step = step+1
      end if
!
      if (r < 1.5) then
          step = step-0.5
      end if
!
      step = step/10d0**kount
!
      if (key == 4) then
          go to 90
      end if
!
      if (xmax-xmin > 0.8*m*step) then
          go to 90
      end if
!
      kount = 1
      key = 2
!
      go to 40
!
   90 xmin = b
      c = step*dint(b/step)
!
      if (c < 0 .and. c /= b) then
          c = c-step
      end if
!
      if (c + m*step >= xmax) then
          xmin = c
      end if
!
      do 100 i = 1,m
         freq(i) = 0
  100 end do
!
      jmax = 1
!
!---- calculate frequencies for each interval
      do 110 i = 1,nmes
!-------- class index
          j = idnint((dat(i) - xmin)/step) + 1
!
          if (j == 0) then
              call wrtlin('    *** warning *** unable to produce '// &
              'histogram for specified number of classes')
              return
          else if (j > m) then
              m = m-1
              go to 15
          else if (j > jmax) then
              jmax = j
          end if
!
!-------- update frequencies for each class
          freq(j) = freq(j) + 1
  110 end do
!
      m = jmax
!
!---- print frequency vector
  120 continue
!
      call newlin
      write (outbuf,'(a,11i6)') '    FREQUENCY ',(freq(i),i=1,m)
      call wrtlin(outbuf)
      write (outbuf,'(77a)') '    ',('_',i=1,m*6+10)
      call wrtlin(outbuf)
!---- find largest frequency and scale if necessary
      imax = 0
!
      do 130 i = 1,m
!
!-------- update largest frequency
          if (freq(i) > imax) then
              imax = freq(i)
          end if
!
  130 end do
!
      scale = 1
!
      if (imax > 16) then
          scale = (imax+15)/16
      end if
!
!---- clear output to blanks
      do 150 i = 1,m
          cout(i) = ' '
  150 end do
!
!---- for each line of print, place output characters in their
!---- appropriate positions in the output vector
      imax = idnint(imax/dble(scale))
!
      do 170 i = imax,1,-1
!
          do 160 j = 1,m
!
              if (idnint(freq(j)/dble(scale)) == i) then
                  cout(j) = '[XXXX]'
              end if
!
  160     end do
!
          k = i*scale
!-------- print line of frequencies
          write (outbuf,'(a,i5,12a)') '        ',k,' ',(cout(j),j=1,m)
          call wrtlin(outbuf)
  170 end do
!
      write (outbuf,'(77a)') '    ',('_',i=1,m*6+10)
      call wrtlin(outbuf)
!
!---- factor
      if (.not. ind) then
!
          if (key == 1) then
              write (outbuf,'(a,10(i2,4x),i2)') '    LEVEL         ', &
              (i,i=1,m)
              call wrtlin(outbuf)
          else if (m*key == mm) then
              write (outbuf,'(a,10(I2,''-'',I2,1X),I2,''-'',I2)') &
              '    LEVELS     ',((i-1)*key+1,i*key,i=1,m)
              call wrtlin(outbuf)
          else
              write (outbuf,'(a,10(I2,''-'',I2,1X),I2,''-'',I2)') &
              '    LEVELS     ',((i-1)*key+1,i*key,i=1,m-1),(m-1)*key+1, &
              mm
              call wrtlin(outbuf)
          end if
!
          go to 900
      end if
!
!---- compute interval mid-points and scale if necessary
      j = 0
      xmax = xmin + step*(m-1)
      xm = dmin1(abs(xmin),abs(xmax))
!
      if (xm < eta) then
          xm = xm+step
      end if
!
  180 if (xm >= 1.0d-1) then
          go to 190
      end if
!
      j = j+1
      xm = xm*10
!
      go to 180
!
  190 xm = dmax1(xmax,-xmin)
!
  200 if (xm < 1.0d3) then
          go to 210
      end if
!
      j = j-1
      xm = xm/10
!
      go to 200
!
  210 tk = 10d0**j
      step = step*tk
      out(1) = xmin*tk
!
      do 220 i = 2,m
          out(i) = out(i-1) + step
  220 end do
!
!---- print interval mid-points
      write (outbuf,'(a,5(F7.2,5X),F7.2)') '    INTERVAL ', &
      (out(i),i=1,m,2)
      call wrtlin(outbuf)
      write (outbuf,'(a,5(F7.2,5X))') '    MID-POINT      ', &
      (out(i),i=2,m,2)
      call wrtlin(outbuf)
!
      if (j /= 0) then
          write (outbuf,'(a,i2)') &
          '    THE PRINTED VALUES MUST BE MULTIPLIED BY 10^',-j
          call wrtlin(outbuf)
      end if
!
  900 continue
!
      call newlin
!
      return
!
      end subroutine hisfrq
!
!***********************************************************************
!
      subroutine modelf(line,modelt)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character modelt
!-----------------------------------------------------------------------
!     Function : Reads the model type.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character c
      integer istpos,ienpos
      logical arg
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** no argument')
          return
      else
!-------- argument - model = univariate/bivariate/trivariate
          c = line(istpos:istpos)
      end if
!
!---- invalid argument
      if (c == 'u' .or. c == 'b' .or. c == 't') then
          modelt = c
      else
          call wrtlin('    *** ERROR *** INCORRECT ARGUMENT')
          call wrtlin('    Type U, B or T')
      end if
!
      return
!
      end subroutine modelf
!
!***********************************************************************
!
      subroutine case(line,cname,icvar,iread,name,nvar,prefix,maxargs, &
                      nlevel)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer maxargs,nvar,nlevel
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: cname(2),name(maxvar)
      character*(*) prefix(maxargs)
      logical iread,icvar(2)
!-----------------------------------------------------------------------
!     Function : Reads the name of the case variable.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=50) :: cold(2)
      integer nargs
!-----------------------------------------------------------------------
      cold(1) = cname(1)
      cold(2) = cname(2)
!
!-----**************
      call multi_var(line,cname,prefix,maxargs,name,nvar,iread,nargs)
!-----**************
!
      nlevel = nargs
!
!---- case-variable specified OK
      if (nlevel == 1) then
          icvar(1) = .true.
      else if (nlevel == 2) then
          icvar(1) = .true.
          icvar(2) = .true.
      else
          cname(1) = cold(1)
          cname(2) = cold(2)
      end if
!
      return
!
      end subroutine case
!
!***********************************************************************
!
      subroutine yvar(line,yname,name,nvar,idata,iread,iyvar)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: yname,name(maxvar)
      integer nvar
      logical iread,iyvar,idata
!-----------------------------------------------------------------------
!     Function : Reads the name of the y-variate.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=50) :: yold
      integer ienpos,istpos,i,dummy
      logical arg
!-----------------------------------------------------------------------
      yold = yname
!
      if (.not. idata .and. .not. iread) then
          call wrtlin( &
          '    *** ERROR *** DATA command must be issued first')
          return
      end if
!
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** MISSING ARGUMENT')
          return
      end if
!
!---- reduce to first 50 characters
      if (ienpos-istpos+1 > 50) then
          yname = line(istpos:istpos+49)
      else
          yname = line(istpos:ienpos)
      end if
!
      dummy = 0
!
      if (.not. iread) then
          dummy = 1
      end if
!
!---- determine whether Y variable exists
      do 10 i = 1,nvar+dummy
!
          if (name(i) == yname) then
              go to 20
          end if
!
   10 end do
!
!---- illegal y-variate name
      call wrtlin('    *** ERROR *** Y-VARIABLE DOES NOT EXIST')
      call wrtlin('    Use DISPLAY VARIABLES to check spelling')
!
      yname = yold
      return
!
   20 iyvar = .true.
!
      return
!
      end subroutine yvar
!
!***********************************************************************
!
      subroutine rvar(line,rname,name,nvar,idata,iread,irvar)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: rname,name(maxvar)
      integer nvar
      logical iread,idata,irvar
!-----------------------------------------------------------------------
!     Function : Reads the name of the risk-variate.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=50) :: rold
      integer ienpos,istpos,i,dummy
      logical arg
!-----------------------------------------------------------------------
      rold = rname
!
      if (.not. idata .and. .not. iread) then
          call wrtlin( &
          '    *** ERROR *** DATA command must be issued first')
          return
      end if
!
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** MISSING ARGUMENT')
          return
      end if
!
!---- reduce to first 50 characters
      if (ienpos-istpos+1 > 50) then
          rname = line(istpos:istpos+49)
      else
          rname = line(istpos:ienpos)
      end if
!
      dummy = 0
!
      if (.not. iread) then
          dummy = 1
      end if
!
!---- determine whether risk variable exists
      do 10 i = 1,nvar+dummy
!
          if (name(i) == rname) then
              go to 20
          end if
!
   10 end do
!
!---- illegal risk-variate name
      call wrtlin('    *** ERROR *** RISK-VARIABLE DOES NOT EXIST')
      call wrtlin('    Use DISPLAY VARIABLES to check spelling')
!
      rname = rold
      return
!
   20 irvar = .true.
!
      return
!
      end subroutine rvar
!
!***********************************************************************
!
      subroutine factor(line,x,maxcol,nmes,ipos,name,ilev,iread,ncol, &
                        nvar,ifail)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: name(maxvar)
      integer nmes,maxcol,ilev(maxvar),ipos(maxvar),ncol,nvar, length
      double precision x(nmes,maxcol)
      logical iread,ifail
!-----------------------------------------------------------------------
!     Function : Converts a variable into a categorical variable.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=50) :: vname,fname
      double precision cut(maxpar),flevs(maxpar)
      integer istpos,ienpos,nlev,i,j,ccount,ixpos
      logical arg,error
!-----------------------------------------------------------------------
      ifail = .false.

      if (.not. iread) then
          call wrtlin('    *** ERROR *** NO DATA READ')
          return
      end if
!
!---- store the arguments in temporary variables
!---- get position of first argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** AT LEAST TWO ARGUMENTS NEEDED')
          return
      end if
!
!---- store variable name in VNAME and check that it exists and that it
!---- is a variable, not a factor
!---- variable name
      vname = line(istpos:ienpos)
!
      do 10 i = 1,nvar
!
          if (name(i) == vname) then
!
!------------ variable already a factor
              if (ilev(i) > 1) then
                  if (length(vname) <= 32) then
                      call wrtlin('    *** ERROR *** `'// &
                      vname(1:length(vname))// &
                      '` IS A FACTOR, NOT A VARIABLE')
                  else
                      call wrtlin('    *** ERROR *** ')
                      call wrtlin(vname(1:length(vname)))
                      call wrtlin(' IS A FACTOR, NOT A VARIABLE')
                  end if
!
                  return
              end if
!
!------------ variable specified OK
              go to 15
          end if
!
   10 end do
!
!---- variable non-existent
      call wrtlin('    *** ERROR *** VARIABLE DOES NOT EXIST')
      vname = ' '
!
      return
!
!---- get position of second argument
!-----*********
   15 call next(istpos,ienpos,ienpos+1,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** AT LEAST TWO ARGUMENTS NEEDED')
          return
      end if
!
!---- factor name
      fname = line(istpos:ienpos)
!
!---- validate factor name
!-----***********
      call namchk(fname,error)
!-----***********
!
      if (error) then
          return
      end if
!
!---- check that factor name does not exist as a variable or a factor
!---- (if it exists as a factor, invite the person to delete it and then
!---- reissue the command)
      do 20 i = 1,nvar
!
          if (name(i) == fname) then
!
!------------ factor name already exists
              if (ilev(i) > 1) then
                  call wrtlin('    *** ERROR *** `'//fname(1:12)// &
                  '` IS ALREADY A FACTOR')
                  call wrtlin('    DELETE it first and then use FACTOR ' &
                  //'if you must')
!------------ factor exists as variable
              else
                  call wrtlin('    *** ERROR *** FACTOR `'// &
                   fname(1:12)//'` EXISTS AS A VARIABLE')
              end if
!
              return
          end if
!
   20 end do
!
!---- See if there are any cutpoints. If there are, read them in.
      ccount = 0
!
!---- get position of next argument
!-----*********
   25 call next(istpos,ienpos,ienpos+1,line,arg)
!-----*********
!
      if (.not. arg) then
          go to 30
      else
!-------- update number of cutpoints
          ccount = ccount+1
!
!-------- too many cutpoints
          if (ccount >= maxvar .or. ccount >= maxpar) then
              call wrtlin('    *** ERROR *** '// &
              'TOO MANY CUTPOINTS')
              return
          end if
!
!-------- ERROR is returned false if the argument has been correctly
!-------- read as a real number
!-------- convert character into real cutpoint
!---------***********
          call charel(cut(ccount),line(istpos:ienpos),ienpos-istpos+1, &
                      error)
!---------***********
!
          if (error) then
              call wrtlin('    *** ERROR *** '// &
              'CUTPOINT IS NOT A REAL NUMBER')
              return
          end if
!
          go to 25
      end if
!
!---- get here when there aren't any more cutpoints to read in
!---- if there aren't any cutpoints, don't bother with all the
!---- following error checks
   30 if (ccount == 0) then
          go to 38
      end if
!
      do 35 i = 2,ccount
!
          if (cut(i) > cut(i-1)) then
              go to 35
!-------- cutpoints not in ascending order
          else
              call wrtlin('    *** ERROR *** '// &
              'CUTPOINTS MUST BE IN ASCENDING ORDER')
              return
          end if
!
   35 end do
!
!---- continuation mark for if there aren't any cutpoints
!---- get position of variable in X-vector
!-----***********
   38 call fndpos(nvar,vname,name,ixpos,ipos)
!-----***********
!
!---- no cutpoints specified
      if (ccount == 0) then
!-------- FLEVS(.) will contain a list of the levels found when
!-------- searching through VNAME. FLEVS(.) is only added to when a new
!-------- level turns up. Of course the first entry is always 'a new
!-------- level'. NLEV is a count of the levels found. The loop on 40
!-------- checks the previous levels to see if the current value is a
!-------- new one.
!-------- set first factor level to first measurement of variable
          flevs(1) = x(1,ixpos)
!-------- initialise number of factor levels found
          nlev = 1
!
          do 50 i = 2,nmes
!
              do 40 j = 1,nlev
!
!---------------- current measurement equals existing factor level
                  if (x(i,ixpos) == flevs(j)) then
                      go to 50
                  end if
!
   40         end do
!
!------------ update number of factor levels so far recognised
              nlev = nlev+1
!
              if (nlev > maxpar) then
                  call wrtlin( &
                  '    *** ERROR *** TOO MANY FACTOR LEVELS')
                  call wrtlin('    You need to specify some cutpoints')
                  return
              end if
!
!------------ store new factor level
              flevs(nlev) = x(i,ixpos)
   50     end do
!
!-------- not enough space for new factor
          if (ncol+nlev > maxcol) then
              error = .true.
              go to 100
          end if
!
!-------- sort factor levels into ascending order
!---------*********
          call sort(flevs,nlev)
!---------*********
!
!-------- create the factor columns in reverse order of the levels.
!-------- This is so that the first level gets aliased, not the last.
!-------- 1 if of ith level; 0 otherwise
          do 70 i = 1,nmes
!
              do 60 j = 1,nlev
!
!---------------- variable value is factor level NLEV-J+1, so set
!---------------- indicator to 1
                  if (x(i,ixpos) == flevs(nlev-j+1)) then
                      x(i,ncol+j) = 1
!---------------- variable value is not factor level NLEV-J+1, so set
!---------------- indicator to 0
                  else
                      x(i,ncol+j) = 0
                  end if
!
   60         end do
!
   70     end do
!
!---- this is the other bit of an if statement upwards a long way, where
!---- the factor is created if there are cutpoints. This code is a bit
!---- suspect because rounding errors could conceivably allow a point to
!---- be bunged into two groups.
!---- cutpoints specified
      else
          nlev = ccount+1
!
!-------- not enough space for new factor
          if (ncol+nlev > maxcol) then
              error = .true.
              go to 100
          end if
!
!-------- create the factor columns in reverse order of the levels.
!-------- This is so that the first level gets aliased, not the last.
!-------- 1 if of ith level; 0 otherwise
          do 90 i = 1,nmes
!
!------------ variable value is less than first cutpoint, so set
!------------ indicator to 1
              if (x(i,ixpos) < cut(1)) then
                  x(i,ncol+nlev) = 1
              else
                  x(i,ncol+nlev) = 0
              end if
!
              do 80 j = 2,nlev-1
!
!---------------- variable value is between (j-1)st and jth cutpoints,
!---------------- so set indicator to 1
                  if (x(i,ixpos) < cut(j) .and. &
                  x(i,ixpos) >= cut(j-1)) then
                      x(i,ncol+nlev+1-j) = 1
                  else
                      x(i,ncol+nlev+1-j) = 0
                  end if
!
   80         end do
!
!------------ variable value is greater than or equal to last cutpoint,
!------------ so set indicator to 1
              if (x(i,ixpos) >= cut(nlev-1)) then
                  x(i,ncol+1) = 1
              else
                  x(i,ncol+1) = 0
              end if
!
   90     end do
!
      end if
!
!---- not enough space for new factor
  100 if (error) then
          call wrtlin('    *** ERROR *** NOT ENOUGH SPACE FOR FACTOR')
          call wrtlin('    You had better DELETE some variates or '// &
          'increase the space')
          call wrtlin('    (see DISPLAY LIMITS for space used)')
          return
      end if
!
!---- update model information for new factor
!---- update number of variables in model
      nvar = nvar+1
!
      if (nvar > maxvar) then
          call wrtlin('    *** ERROR *** TOO MANY VARIABLES')
          ifail = .true.
          return
      end if
!
!---- store position of new factor in X-matrix
      ipos(nvar) = ncol+1
!---- store number of levels of new factor
      ilev(nvar) = nlev
!---- store new factor name
      name(nvar) = fname
!---- update number of columns in X-matrix
      ncol = ncol+nlev
!
      return
!
      end subroutine factor
!
!***********************************************************************
!
      subroutine sort(flevs,nlev)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer nlev
      double precision flevs(nlev)
!-----------------------------------------------------------------------
!     Function : Sorts factor levels into ascending order.
!-----------------------------------------------------------------------
!     Shell sort (D.L. Shell). Also called the "Merge-Exchange Sort".
!     With this algorithm the functional form of the running time is not
!     known (and depends on the increment sequence). For the following
!     case, two conjectures are NLEV*((log(NLEV))^2) and NLEV^1.25; the
!     complete sort is proportional to NLEV^2.
!     Programmer: M.A. Bradley, Date: 5-1-89
!-----------------------------------------------------------------------
      double precision t
      integer h,i,j
!-----------------------------------------------------------------------
      h = 1
!
 1000 if (h <= nlev) then
          h = 3*h + 1
          go to 1000
      end if
!
 2000 if (h /= 1) then
          h = h/3
!
!-------- sort vector into ascending order
          do 6000 i = h+1,nlev
              t = flevs(i)
              j = i
!
 4000         if (flevs(j-h) > t) then
                  flevs(j) = flevs(j-h)
                  j = j-h
!
                  if (j <= h) then
                      go to 5000
                  end if
!
                  go to 4000
              end if
!
 5000         flevs(j) = t
 6000     end do
!
          go to 2000
      end if
!
      return
!
      end subroutine sort
!
!***********************************************************************
!
      subroutine trans(line,x,maxcol,nmes,ipos,name,ilev,iread,ncol, &
                       nvar,yname,rname,cname,nlevel,ifail)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: name(maxvar),yname,rname,cname(2)
      integer nmes,maxcol,nlevel,nvar,ipos(maxvar),ncol,ilev(maxvar)
      double precision x(nmes,maxcol)
      logical iread,ifail
!-----------------------------------------------------------------------
!     Function : Performs simple data transformations and also forms
!                interactions between variables and/or factors.
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=3) :: func
      character(len=50) :: vname(3)
      double precision power
      double precision :: cnstnt = 0d0
      integer istpos,ienpos,i,k,dummy,ixpos(3),ixlev(3),marker,nlev,l, &
              ivpos,ifpos,j
      logical arg,error
!-----------------------------------------------------------------------
      ifail = .false.

      if (.not. iread) then
          call wrtlin('    *** ERROR *** NO DATA READ')
          return
      end if
!
      if (ncol+1 > maxcol) then
          call wrtlin('    *** ERROR *** '// &
          'NOT ENOUGH SPACE FOR TRANSFORMED VARIATE')
          return
      end if
!
!---- get position of first argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** NO TRANSFORMED VARIATE GIVEN')
          return
      end if
!
!---- transformed variable name
      vname(1) = line(istpos:ienpos)
!
!---- validate transformed variable name
!-----***********
      call namchk(vname(1),error)
!-----***********
!
      if (error) then
          return
      end if
!
      do 10 i = 1,nvar
!
          if (name(i) == vname(1)) then
              call wrtlin('    *** ERROR *** '// &
              'TRANSFORMED VARIATE ALREADY EXISTS')
              return
          end if
!
   10 end do
!
!---- initialise the number of levels of the transformed variable/factor
      nlev = 1
!---- initialise indicator for command element most recently read
      k = 1
!---- initialise flag for variable or constant arithmetic argument
      marker = 0
!---- initialise the function / arithmetic operator
      func = ' '
      dummy = ienpos+1
!
!---- get position of first transformation argument
!-----*********
      call next(istpos,ienpos,dummy,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** NO TRANSFORMATION GIVEN')
          return
      end if
!
!---- skip this bit if no arguments have yet been read
   15 if (k >= 2) then
!
!-------- store the name of the function argument / variable
          vname(k) = line(istpos:ienpos)
!
!-------- disallow transformations of the y-variate
          if (vname(k) == yname) then
              call wrtlin('    *** ERROR *** TRANSFORMATIONS OF '// &
              'THE Y-VARIATE ARE NOT ALLOWED')
              return
!-------- disallow transformations of the risk-variate
          else if (vname(k) == rname) then
              call wrtlin('    *** ERROR *** TRANSFORMATIONS OF '// &
              'THE risk-VARIATE ARE NOT ALLOWED')
              return
!-------- disallow transformations of the case variate
          else if (vname(k) == cname(1) .or. (nlevel == 2 .and. &
          vname(k) == cname(2))) then
              call wrtlin('    *** ERROR *** TRANSFORMATIONS OF '// &
              'THE CASE VARIATE ARE NOT ALLOWED')
              return
          end if
!
          do 20 i = 1,nvar
!
              if (name(i) == vname(k)) then
                  go to 30
              end if
!
   20     end do
!
          if (func /= ' ' .and. k == 2) then
              call wrtlin('    *** ERROR *** INVALID FUNCTION ARGUMENT')
          else if (k == 2) then
              call wrtlin('    *** ERROR *** INVALID TRANSFORMATION')
!-------- set the flag to indicate non-variable arithmetic argument
          else
              marker = 1
              go to 35
          end if
!
          return
!
!-------- variable name specified OK; get its position in the x-vector
!---------***********
   30     call fndpos(nvar,vname(k),name,ixpos(k),ipos)
!---------***********
!
!-------- get the number of levels of the variable
!---------***********
          call fndlev(nvar,vname(k),name,ixlev(k),ilev)
!---------***********
!
      end if
!
!---- function transformation
   35 if (line(istpos:ienpos) == 'exp' .or. &
      line(istpos:ienpos) == 'log' .or. (func /= ' ' .and. &
      k == 2)) then
!
          if (k == 2 .and. ixlev(2) > 1) then
              call wrtlin('    *** ERROR *** FUNCTIONAL '// &
              'TRANSFORMATIONS OF FACTORS ARE NOT ALLOWED')
              return
          end if
!
          if (func == ' ') then
!------------ store the function name
              func = line(istpos:ienpos)
              dummy = ienpos+1
!
!------------ get the function argument
!-------------*********
              call next(istpos,ienpos,dummy,line,arg)
!-------------*********
!
              if (.not. arg) then
                  call wrtlin('    *** ERROR *** '// &
                  'NO TRANSFORMATION ARGUMENT GIVEN')
                  return
              end if
!
!------------ reset flag for function argument
              k = 2
              go to 15
          end if
!
!-------- exponentiation
          if (func == 'exp') then
!
              do 40 i = 1,nmes
!
!---------------- prevent numerical overflow or underflow
                  if (x(i,ixpos(2)) > zl1 .or. &
                  x(i,ixpos(2)) < zl2) then
                      call wrtlin('    *** ERROR *** VARIATE '// &
                      'LIES OUTSIDE EXPONENTIATION RANGE')
                      return
!---------------- exponential transformation
                  else
                      x(i,ncol+1) = exp(x(i,ixpos(2)))
                  end if
!
   40         end do
!
!-------- natural logarithm
          else
!
              do 50 i = 1,nmes
!
                  if (x(i,ixpos(2)) <= 0) then
                      call wrtlin('    *** ERROR *** '// &
                      'VARIATE IS NON-POSITIVE')
                      return
!---------------- logarithmic transformation
                  else
                      x(i,ncol+1) = log(x(i,ixpos(2)))
                  end if
!
   50         end do
!
          end if
!
!---- arithmetic transformation
      else
!
!-------- reset flag for first arithmetic argument
          if (k == 1) then
              k = 2
              go to 15
!-------- both arithmetic arguments and the operator have been read
          else if (k == 3) then
              go to 60
          end if
!
          dummy = ienpos+1
!
!-------- get position of operator
!---------*********
          call next(istpos,ienpos,dummy,line,arg)
!---------*********
!
          if (.not. arg) then
              call wrtlin('    *** ERROR *** '// &
              'NO ARITHMETIC OPERATOR GIVEN')
              return
          end if
!
!-------- arithmetic operator
          func = line(istpos:ienpos)
!
          if (func /= '^' .and. func /= '*' .and. func /= '/' &
          .and. func /= '+' .and. func /= '-') then
              call wrtlin('    *** ERROR *** '// &
              'INVALID ARITHMETIC OPERATOR')
              return
          end if
!
          dummy = ienpos+1
!
!-------- get position of second variable / constant
!---------*********
          call next(istpos,ienpos,dummy,line,arg)
!---------*********
!
          if (.not. arg) then
              call wrtlin('    *** ERROR *** '// &
              'INCOMPLETE ARITHMETIC OPERATION')
              return
          end if
!
!-------- power transformation
   60     if (func == '^') then
!
              if (ixlev(2) > 1) then
                  call wrtlin('    *** ERROR *** POWER TRANSFORMATIONS ' &
                  //'OF FACTORS ARE NOT ALLOWED')
                  return
              end if
!
!------------ convert the power into a real number
!-------------***********
              call charel(power,line(istpos:ienpos),ienpos-istpos+1, &
                          error)
!-------------***********
!
              if (error) then
                  call wrtlin('    *** ERROR *** '// &
                  'INVALID POWER TRANSFORMATION')
                  return
              end if
!
              do 70 i = 1,nmes
!
!---------------- negative values may only be raised to integer powers
                  if (x(i,ixpos(2)) < 0 .and. &
                  idint(power) /= power) then
                      call wrtlin('    *** ERROR *** CANNOT RAISE '// &
                      'NEGATIVE VALUES TO NON-INTEGER POWERS')
                      return
!---------------- power zero variable value
                  else if (idint(x(i,ixpos(2))) == 0) then
                      x(i,ncol+1) = 0
!---------------- prevent numerical overflow
                  else if (power*log(abs(x(i,ixpos(2)))) > zl1) &
                  then
                      call wrtlin('    *** ERROR *** '// &
                      'POWERED VALUES TOO LARGE')
                      return
!---------------- prevent numerical underflow
                  else if (power*log(abs(x(i,ixpos(2)))) < zl2) &
                  then
                      call wrtlin('    *** ERROR *** '// &
                      'POWERED VALUES TOO SMALL')
                      return
!---------------- power the explanatory variable
                  else
                      x(i,ncol+1) = x(i,ixpos(2))**power
                  end if
!
   70         end do
!
!-------- arithmetic operation
          else
!
!------------ reset flag for second arithmetic argument
              if (k == 2) then
                  k = 3
                  go to 15
              end if
!
              if (((ixlev(2) > 1 .or. (marker == 0 .and. &
              ixlev(3) > 1)) .and. func /= '*') .or. &
              (ixlev(2) > 1 .and. marker == 1 .and. &
              func == '*')) then
                  call wrtlin('    *** ERROR *** NON-INTERACTION '// &
                  'TRANSFORMATIONS OF FACTORS ARE NOT ALLOWED')
                  return
              end if
!
!------------ the second arithmetic argument was a constant
              if (marker == 1) then
!
!---------------- convert the constant into a real number
!-----------------***********
                  call charel(cnstnt,line(istpos:ienpos), &
                              ienpos-istpos+1,error)
!-----------------***********
!
                  if (error) then
                      call wrtlin('    *** ERROR *** '// &
                      'INVALID ARITHMETIC OPERATION')
                      return
                  end if
!
                  if (func == '/' .and. cnstnt == 0) then
                      call wrtlin('    *** ERROR *** GET REAL MATEY'// &
                      ' - DIVISION BY ZERO IS NOT ALLOWED')
                      return
                  end if
!
!---------------- set dummy number of levels for the constant to be one
                  ixlev(3) = 1
!---------------- set dummy X-position for the constant to be one
                  ixpos(3) = 1
              end if
!
!------------ product transformation
              if (func == '*') then
!---------------- number of levels in the interaction
                  nlev = ixlev(2) * ixlev(3)
!
                  if (ncol+nlev > maxcol) then
                      call wrtlin('    *** ERROR *** '// &
                      'NOT ENOUGH SPACE FOR INTERACTION')
                      call wrtlin('    You had better DELETE some '// &
                      'variates or increase the space')
                      call wrtlin( &
                      '    (see DISPLAY LIMITS for space used)')
                      return
                  end if
!
                  if (ixlev(2) == 1 .and. ixlev(3) == 1) then
!
!-------------------- multiply explanatory variable by variate/constant
                      do 130 i = 1,nmes
                          x(i,ncol+1) = x(i,ixpos(2)) * &
                          (marker*cnstnt + (1-marker)*x(i,ixpos(3)))
  130                 end do
!
                  else if (ixlev(2) == 1 .or. ixlev(3) == 1) then
!
!-------------------- variable by factor interaction
                      if (ixlev(2) == 1) then
                          ivpos = ixpos(2)
                          ifpos = ixpos(3)
!-------------------- factor by variable interaction
                      else
                          ivpos = ixpos(3)
                          ifpos = ixpos(2)
                      end if
!
                      do 150 i = 1,nmes
!
!                          IF (X(I,IVPOS) == 0.0) THEN
!                              call wrtlin('    *** ERROR *** VARIABLE '
!     &                        //'MUST BE NON-ZERO EVERYWHERE')
!                              RETURN
!                          END IF
!
                          do 140 j = 1,nlev
                              x(i,ncol+j) = x(i,ivpos)*x(i,ifpos+j-1)
  140                     end do
!
  150                 end do
!
!---------------- REMEMBER THAT THE INTERACTION FACTOR THING NEEDS TO BE
!---------------- STORED IN REVERSE ORDER
!---------------- factor by factor interaction
                  else
!
                      do 180 i = 1,nmes
!
                          do 170 j = 1,ixlev(2)
!
                              do 160 l = 1,ixlev(3)
                                  x(i,ncol+nlev+1-(j-1)*ixlev(3)-l) = &
                                  x(i,ixpos(2)+ixlev(2)-j)* &
                                  x(i,ixpos(3)+ixlev(3)-l)
  160                         end do
!
  170                     end do
!
  180                 end do
!
                  end if
!
!------------ ratio transformation
              else if (func == '/') then
!
                  do 190 i = 1,nmes
!
!-------------------- avoid division by zero
                      if &
                      (marker*cnstnt + (1-marker)*x(i,ixpos(3)) == 0) &
                      then
                          call wrtlin('    *** ERROR *** '// &
                          'DIVISOR VARIATE IS ZERO')
                          return
                      end if
!
!-------------------- divide explanatory variable by variate/constant
                      x(i,ncol+1) = x(i,ixpos(2))/ &
                      (marker*cnstnt + (1-marker)*x(i,ixpos(3)))
  190             end do
!
!------------ sum transformation
              else if (func == '+') then
!
!---------------- add the variate/constant to the explanatory variable
                  do 200 i = 1,nmes
                      x(i,ncol+1) = x(i,ixpos(2)) + &
                      (marker*cnstnt + (1-marker)*x(i,ixpos(3)))
  200             end do
!
!------------ difference transformation
              else
!
!---------------- subtract variate/constant from explanatory variable
                  do 210 i = 1,nmes
                      x(i,ncol+1) = x(i,ixpos(2)) - &
                      (marker*cnstnt + (1-marker)*x(i,ixpos(3)))
  210             end do
!
              end if
!
          end if
!
      end if
!
!---- update model information for transformation
!---- update number of variates
      nvar = nvar+1
!
      if (nvar > maxvar) then
          call wrtlin('    *** ERROR *** TOO MANY VARIABLES')
          ifail = .true.
          return
      end if
!
!---- store position of transformed variate in X-vector
      ipos(nvar) = ncol+1
!---- update number of columns in X-matrix
      ncol = ncol+nlev
!---- store number of levels in transformed variate
      ilev(nvar) = nlev
!----- store the name of the transformed variate
      name(nvar) = vname(1)
!
      return
!
      end subroutine trans
!
!***********************************************************************
!
      subroutine familyf(line,family,prefix,maxargs)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs
      character*(*) family(maxargs),prefix(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads the family/families.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character f
      integer nargs,iarg,arg_index(3)
!-----------------------------------------------------------------------
!-----***************
      call multi_char(line,family,prefix,maxargs,nargs,arg_index)
!-----***************
!
      if (nargs == 0) then
          call wrtlin('    *** ERROR *** no argument')
          return
      else
!
          do 10 iarg = 1,nargs
!
!------------ check for valid argument - family = binomial/gaussian/
!             poisson
              f = family(arg_index(iarg))
!
              if (f /= 'b' .and. f /= 'g' .and. f /= 'p' .and. f /= 'o') &
              then
                  call wrtlin('    *** ERROR *** INCORRECT ARGUMENT')
                  call wrtlin('    Type B, G, P or O')
              end if
!
 10       end do
!
      end if
!
      return
!
      end subroutine familyf
!
!***********************************************************************
!
      subroutine linkf(line,link,prefix,maxargs)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs
      character*(*) link(maxargs),prefix(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads the link function.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character l
      integer nargs,iarg,arg_index(3)
!-----------------------------------------------------------------------
!-----***************
      call multi_char(line,link,prefix,maxargs,nargs,arg_index)
!-----***************
!
      if (nargs == 0) then
          call wrtlin('    *** ERROR *** no argument')
          return
      else
!
          do 10 iarg = 1,nargs
!------------ check for valid argument - logit, cll or probit
              l = link(arg_index(iarg))
!
              if ( l /= 'l' .and. l /= 'c' .and. l /= 'p') then
                  call wrtlin('    *** ERROR *** INCORRECT ARGUMENT')
                  call wrtlin('    Type L, C or P')
              end if
!
 10       end do
!
      end if
!
      return
!
      end subroutine linkf
!
!***********************************************************************
!
      subroutine numlp(line,nnlev,prefix,maxargs)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs,nnlev(maxargs),nargs,iarg
      character*(*) prefix(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads the number of variables in the linear
!                predictors.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer arg_index(3)
!-----------------------------------------------------------------------
!-----**************
      call multi_int(line,nnlev,prefix,maxargs,nargs,arg_index)
!-----**************
!
!---- Check for non-negative values
      do 10 iarg = 1,nargs
!
          if (nnlev(arg_index(iarg)) < 0) then
              call wrtlin('    *** ERROR *** '// &
                          'ARGUMENT IS NOT A NON-NEGATIVE INTEGER')
          end if
!
 10   end do
!
      return
!
      end subroutine numlp
!
!***********************************************************************
!
      subroutine offset(line,offnme,offlag,iread,name,nvar,prefix, &
                        maxargs)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer nvar,maxargs
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: offnme(3),name(maxvar)
      character*(*) prefix(maxargs)
      logical offlag(3),iread
!-----------------------------------------------------------------------
!     Function: To specify an offset included in the linear predictor
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=50) :: oold(3)
      integer istpos,ienpos,i,nargs
      logical arg
!-----------------------------------------------------------------------
      oold(1) = offnme(1)
      oold(2) = offnme(2)
      oold(3) = offnme(3)
!
!----**************
     call multi_var(line,offnme,prefix,maxargs,name,nvar,iread,nargs)
!----**************
!
      if (nargs == 1) then
          offlag(1) = .true.
      else if (nargs == 2) then
          offlag(1) = .true.
          offlag(2) = .true.
      else if (nargs == 3) then
          offlag(1) = .true.
          offlag(2) = .true.
          offlag(3) = .true.
      else
          offnme(1) = oold(1)
          offnme(2) = oold(2)
          offnme(3) = oold(3)
      end if
!
      return
!
      end subroutine offset
!
!***********************************************************************
!
      subroutine offset_old(line,offnme,offlag,iread,name,nvar)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: offnme,name(maxvar)
      integer nvar
      logical offlag,iread
!-----------------------------------------------------------------------
!     Function: To specify an offset included in the linear predictor
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=50) :: oold
      integer istpos,ienpos,i,nargs
      logical arg
!-----------------------------------------------------------------------
      oold = offnme
!
      if (.not. iread) then
          call wrtlin('    *** ERROR *** NO DATA READ')
          return
      end if
!
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          offlag = .false.
          return
      end if
!
      if (ienpos-istpos+1 > 50) then
          offnme = line(istpos:istpos+49)
      else
          offnme = line(istpos:ienpos)
      end if
!
      do 10 i = 1,nvar
!
          if (name(i) == offnme) then
              go to 20
          end if
!
   10 end do
!
      call wrtlin('    *** ERROR *** OFFSET VARIABLE DOES NOT EXIST')
!
      offnme = oold
      return
!
   20 offlag = .true.
!
      return
!
      end subroutine offset_old
!
!***********************************************************************
!
      subroutine arithm(line,arith)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character arith
!-----------------------------------------------------------------------
!     Function : Reads the arithmetic mode.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character c
      integer istpos,ienpos
      logical arg
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** no argument')
          return
      else
!-------- argument - arithmetic = accurate/fast
          c = line(istpos:istpos)
      end if
!
!---- invalid argument
      if (c == 'a' .or. c == 'f') then
          arith = c
      else
          call wrtlin('    *** ERROR *** INCORRECT ARGUMENT')
          call wrtlin('    Type A or F')
      end if
!
      return
!
      end subroutine arithm
!
!***********************************************************************
!
      subroutine quads(line,cquad)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character cquad
!-----------------------------------------------------------------------
!     Function : Reads the quadrature mode.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character c
      integer istpos,ienpos
      logical arg
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** no argument')
          return
      else
!-------- argument - quadrature = Gaussian/adaptive
          c = line(istpos:istpos)
      end if
!
!---- invalid argument
      if (c == 'g' .or. c == 'a' .or. c == 'b' .or. c == 'c' &
      .or. c == 'd') then
          cquad = c
      else
          call wrtlin('    *** ERROR *** INCORRECT ARGUMENT')
          call wrtlin('    Type G or A or B or C or D')
      end if
!
      return
!
      end subroutine quads
!
!***********************************************************************
!
      subroutine iter(line,niter)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer niter
!-----------------------------------------------------------------------
!     Function : Reads the maximum number of iterations.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer itemp
!-----------------------------------------------------------------------
!---- get maximum number of iterations
!-----***********
      call argint(line,niter,itemp)
!-----***********
!
!---- invalid maximum number of iterations
      if (niter <= 0) then
          call wrtlin('    *** ERROR *** '// &
          'YOU MUST HAVE AT LEAST ONE ITERATION')
          niter = itemp
      end if
!
      return
!
      end subroutine iter
!
!***********************************************************************
!
      subroutine meil(line,nmeil,niter)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer nmeil,niter
!-----------------------------------------------------------------------
!     Function : Reads the number of approximate (Meilijson) iterations.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer itemp
!-----------------------------------------------------------------------
!---- get number of approximate iterations
!-----***********
      call argint(line,nmeil,itemp)
!-----***********
!
!---- invalid number of approximate iterations
      if (nmeil < 0) then
          call wrtlin('    *** ERROR *** '// &
          'CANNOT BE NEGATIVE')
          nmeil = itemp
      else if (nmeil > niter) then
          call wrtlin('    *** ERROR *** CANNOT BE GREATER THAN '// &
          'MAXIMUM NUMBER OF ITERATIONS')
          nmeil = itemp
      end if
!
      return
!
      end subroutine meil
!
!***********************************************************************
!
      subroutine mass(line,nmp,prefix,maxargs)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs,nmp(maxargs),nargs,m
      character*(*) prefix(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads the number of mass points for up to maxargs
!                responses.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer iarg,itemp(maxvar),arg_index(3)
!-----------------------------------------------------------------------
!---- Save any mass point values already specified
      do 10 iarg = 1,maxargs
          itemp(iarg) = nmp(iarg)
 10   end do
!
!---- get number of mass points for up to maxargs responses
!-----**************
      call multi_int(line,nmp,prefix,maxargs,nargs,arg_index)
!-----**************
!
!---- check for invalid mass point values
      do 20 iarg = 1,nargs
          m = nmp(arg_index(iarg))
!
          if (m /= 1 .and. m /= 2 .and. m /= 4 .and. m /= 6 &
          .and. m /= 8 .and. m /= 10 .and. m /= 12 .and. m /= 14 &
          .and. m /= 16 .and. m /= 20 .and. m /= 24 .and. &
          m /= 28 .and. m /= 32 .and. m /= 36 .and. m /= 40 &
          .and. m /= 44 .and. m /= 48 .and. m /= 56 .and. &
          m /= 64 .and. m /= 72 .and. m /= 80 .and. m /= 88 &
          .and. m /= 96 .and. m /= 104 .and. m /= 112 .and. &
          m /= 128 .and. m /= 144 .and. m /= 160 .and. m /= 176 &
          .and. m /= 192 .and. m /= 208 .and. m /= 224 .and. &
          m /= 240 .and. m /= 256) then
              call wrtlin('    *** ERROR *** the number of mass '// &
                          'points must be 1,2,4,6,8,10,12,14,16,20,')
              call wrtlin('    24,28,32,36,40,44,48,56,64,72,80,88,'// &
                          '96,104,112,128,144,160,176,192,208,')
              call wrtlin('    224,240, or 256')
              nmp(iarg) = itemp(iarg)
          end if
!
 20   end do
!
      return
!
      end subroutine mass
!
!***********************************************************************
!
      subroutine endpts(line,endind,est0,est1,family)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character endind,family(3)
      double precision est0,est1
!-----------------------------------------------------------------------
!     Function : Reads the endpoint status and their initial estimates.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character c
      double precision temp0,temp1,temp(2)
      integer istpos,ienpos,i,k,dummy
      logical arg,error
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** no argument')
          return
      else
!-------- argument - endpoints = both/left/right/none
          c = line(istpos:ienpos)
      end if
!
      if ((family(1) == 'b' .and. (c == 'b' .or. c == 'l' .or. &
      c == 'r' .or. c == 'n')) .or. (family(1) == 'p' .and. &
      (c == 'l' .or. c == 'n')) .or. family(1) == 'g' .and. &
      c == 'n') then
!-------- valid endpoint argument
          endind = c
      else
          call wrtlin('    *** ERROR *** INCORRECT ARGUMENT')
!
          if (family(1) == 'b') then
              call wrtlin('    Type B, L, R or N')
          else if (family(1) == 'p') then
              call wrtlin('    Type L or N')
          else
              call wrtlin('    Type N')
          end if
!
          go to 200
      end if
!
      if (endind == 'n') then
          return
      else if (endind /= 'b') then
          k = 1
      else
          k = 2
      end if
!
!---- store default settings
      temp0 = est0
      temp1 = est1
!
!---- put the estimates into EST0 and EST1
      do 10 i = 1,k
          dummy = ienpos+1
!
!-------- get position of next argument
!---------*********
          call next(istpos,ienpos,dummy,line,arg)
!---------*********
!
          if (.not. arg) then
!
              if (i == 2) then
                  call wrtlin('    *** ERROR *** '// &
                  'ESTIMATES FOR BOTH ENDPOINTS NEEDED')
              end if
!
              go to 200
          end if
!
!-------- convert character to real initial estimate of endpoint
!---------***********
          call charel(temp(i),line(istpos:ienpos),ienpos-istpos+1,error)
!---------***********
!
          if (error) then
              call wrtlin('    *** ERROR *** VALUE IS NOT REAL')
              go to 200
          end if
!
   10 end do
!
!---- initial estimate of left endpoint
      if (endind /= 'r') then
          est0 = temp(1)
      end if
!
!---- initial estimate of right endpoint
      if (endind /= 'l') then
          est1 = temp(k)
      end if
!
!---- invalid initial endpoint estimates
      if (est0 < 0 .or. est1 < 0) then
          call wrtlin( &
          '    *** ERROR *** ESTIMATES MUST NOT BE NEGATIVE')
          go to 200
      end if
!
      return
!
!---- reset default settings
  200 est0 = temp0
      est1 = temp1
!
      return
!
      end subroutine endpts
!
!***********************************************************************
!
      subroutine inits(line,rarg)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      double precision rarg(maxpar)
!-----------------------------------------------------------------------
!     Function : Reads a single non-negative real argument.
!-----------------------------------------------------------------------
!     RARG : real valued argument
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer istpos,ienpos,i
      logical arg,error
!-----------------------------------------------------------------------
      ienpos = 3
!
      do 10 i = 1,maxpar
!
!-------- get position of argument
!---------*********
          call next(istpos,ienpos,ienpos+1,line,arg)
!---------*********
!
          if (.not. arg) then
              return
          end if
!
!-------- convert character to real
!---------***********
          call charel(rarg(i),line(istpos:ienpos),ienpos-istpos+1,error)
!---------***********
!
!-------- ERROR is returned false if the argument has been correctly
!-------- read as a real
          if (error) then
              call wrtlin('    *** ERROR *** ARGUMENT IS NOT A REAL')
          end if
!
   10 end do
!
      return
!
      end subroutine inits
!
!***********************************************************************
!
      subroutine disp(line,name,yname,rname,cname,ncol,ni,beta,alias, &
                      cov,ipos,nm,sca,con,alp,ilfit,npar,nmes,nsub, &
                      endind,nvar,ilev,iend,iread,xll,ilevdf,iredf, &
                      ienddf,mnames,arith,tol,est0,est1,nmeil,niter, &
                      itotdf,offlag,offnme,link,corr,bivar, &
                      robust,order,ncat,sig,sig_e,inormdf,family,trivar, &
                      univar,n1var,n2var,rho,nlevel,ifail,depend, &
                      eqscale,cquad,icutdf)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character endind,arith,link(3),family(3),cquad
      character(len=50) :: yname,rname,cname(2),name(maxvar), &
                mnames(maxvar),offnme(3)
      double precision beta(maxpar),sca(3),con,xll,alp,tol,est0, &
             cov(maxpar*(maxpar+1)/2),est1,alias(maxpar),sig(3), &
             sig_e(3),rho(3)
      integer ncol,ni,ipos(maxvar),nm(3),npar,ilfit,nvar,icutdf, &
              ilev(maxvar),iend(2),nmes,nsub(2),ilevdf,ncat(3),iredf, &
              ienddf,itotdf,nmeil,niter,inormdf,n1var,n2var,nlevel
      logical iread,trivar,offlag(3),corr,bivar,robust,order, &
              univar,ifail,depend,eqscale
!-----------------------------------------------------------------------
!     Function : Displays information on the status of the analysis.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character c
      integer ienpos,istpos
      logical arg
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** NO ARGUMENT GIVEN')
          return
      end if
!
!---- reduce to initial character
      c = line(istpos:istpos)
!
!---- display variables
      if (c == 'v') then
!
!---------***********
          call disvar(name,nvar,yname,rname,cname,ilev,iread,offlag, &
                      offnme,nlevel)
!---------***********
!
!---- display limits
      else if (c == 'l') then
!
!---------***********
          call dislim(ncol,nmes,npar,nsub,ilfit,family)
!---------***********
!
!---- display settings
      else if (c == 's') then
!
!---------***********
          call disset(nm,sca,con,alp,endind,tol,est0,est1,nmeil,niter, &
                      arith,offlag,link,corr,bivar,robust,order, &
                      sig,family,trivar,univar,n1var,n2var,rho,nlevel, &
                      depend,cquad)
!---------***********
!
!---- display model
      else if (c == 'm' .and. .not. ifail) then
!
!---------***********
          call dismod(xll,ilevdf,iredf,ienddf,ilfit,mnames,ni,cname, &
                      yname,nsub,nmes,itotdf,endind,link,corr, &
                      bivar,robust,order,inormdf,family,trivar,univar, &
                      offlag,offnme,nlevel,depend,icutdf)
!---------***********
!
!---- display estimates
      else if (c == 'e' .and. .not. ifail) then
!
!---------***********
          call disest(beta,cov,alias,ipos,mnames,ni,nvar,ilfit,npar, &
                      ilev,iend,endind,corr,robust,order,ncat,sig_e, &
                      family,trivar,univar,nlevel,depend,eqscale)
!---------***********
!
      else if (c == 'm' .or. c == 'e') then
          call wrtlin('    Current model invalid')
      else
          call wrtlin('    *** ERROR *** INVALID ARGUMENT')
      end if
!
      return
!
      end subroutine disp
!
!***********************************************************************
!
      subroutine dislim(ncol,nmes,npar,nsub,ilfit,family)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character family(3)
      integer nmes,npar,ncol,nsub(2),ilfit
!-----------------------------------------------------------------------
!     Function : Displays the workspace restrictions.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      integer itemp
!-----------------------------------------------------------------------
      call newlin
      call wrtlin('    Data array')
      call wrtlin('    __________')
      write (outbuf,'(a,i9)') &
      '    Space installed:______________________',mxx
      call wrtlin(outbuf)
      write (outbuf,'(a,i9)') &
      '    Space used by data:                   ',nmes*ncol
      call wrtlin(outbuf)
!
      itemp = mxx - nmes*ncol
!
!---- X-array limits for logistic/log-linear model
      if (ilfit == 1) then
          write (outbuf,'(a,i9)') &
          '    Workspace for irls algorithm:         ', &
          nmes*2 + npar*(npar+2)
          call wrtlin(outbuf)
          itemp = itemp - 2*nmes - npar*(npar+2)
      end if
!
      write (outbuf,'(a,i9)') &
      '    Space left:___________________________',itemp
      call wrtlin(outbuf)
      itemp = 4*maxy
      call newlin
      call wrtlin('    Response array')
      call wrtlin('    ______________')
      write (outbuf,'(a,i9)') &
      '    Space installed:______________________',itemp
      call wrtlin(outbuf)
!
      if (ilfit >= 0) then
          write (outbuf,'(a,i9)') &
          '    Space used by y-variate:              ',nmes
          call wrtlin(outbuf)
      end if
!
!---- y-array limits for logistic/log-linear model
      if (ilfit == 1) then
          itemp = itemp-nmes
!---- y-array limits for mixture model
      else if (ilfit == 0) then
          write (outbuf,'(a,i9)') &
          '    Space used by case structure:         ',nsub(1)
          call wrtlin(outbuf)
!
          if (family(1) == 'p') then
              write (outbuf,'(a,i9)') &
              '    Space used by y-variate products:     ',nsub(1)
              call wrtlin(outbuf)
              itemp = itemp - nmes - nsub(1)
          else
              write (outbuf,'(a,i9)') &
              '    Space used by y-variate products:     ',2*nsub(1)
              call wrtlin(outbuf)
              itemp = itemp - nmes - 2*nsub(1)
          end if
!
      end if
!
      write (outbuf,'(a,i9)') &
      '    Space left:___________________________',itemp
      call wrtlin(outbuf)
      itemp = maxpar*(maxpar+5)
      call newlin
      call wrtlin('    Parameter array')
      call wrtlin('    _______________')
      write (outbuf,'(a,i9)') &
      '    Space installed:______________________',itemp
      call wrtlin(outbuf)
!
!---- beta array limits
      if (ilfit >= 0) then
          write (outbuf,'(a,i9)') &
          '    Space used by estimates:              ',npar
          call wrtlin(outbuf)
          write (outbuf,'(a,i9)') &
          '    Space used by aliasing:               ',2*npar
          call wrtlin(outbuf)
          itemp = itemp - 3*npar
      end if
!
!---- beta array limits for logistic/log-linear model
      if (ilfit == 1) then
          write (outbuf,'(a,i9)') &
          '    Space used by covariances:            ',npar*(npar+1)/2
          call wrtlin(outbuf)
          itemp = itemp - npar*(npar+1)/2
!---- beta array limits for mixture model
      else if (ilfit == 0 .or. ilfit == 2) then
          write (outbuf,'(a,i9)') &
          '    Space used by scores:                 ',npar
          call wrtlin(outbuf)
          write (outbuf,'(a,i9)') &
          '    Space used by Hessian:                ',npar*(npar+1)/2
          call wrtlin(outbuf)
          write (outbuf,'(a,i9)') &
          '    Space used by covariances:            ',npar*(npar+1)/2
          call wrtlin(outbuf)
          itemp = itemp - npar*(npar+2)
      end if
!
      write (outbuf,'(a,i9)') &
      '    Space left:___________________________',itemp
      call wrtlin(outbuf)
      call newlin
      write (outbuf,'(a,i6)') '    Maximum number of variables  = ', &
      maxvar
      call wrtlin(outbuf)
      write (outbuf,'(a,i6)') '    Maximum number of parameters = ', &
      maxpar
      call wrtlin(outbuf)
      call newlin
!
      return
!
      end subroutine dislim
!
!***********************************************************************
!
      subroutine disvar(name,nvar,yname,rname,cname,ilev,iread,offlag, &
                        offnme,nlevel)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character(len=50) :: name(maxvar),yname,rname,cname(2),offnme(3)
      integer ilev(maxvar),nvar,nlevel
      logical iread,offlag(3)
!-----------------------------------------------------------------------
!     Function : Displays the list of variables in current worksheet.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      character vartype*12,space*76,underline*76
      integer i,longest,length
!-----------------------------------------------------------------------
      data space /' '/
!-----------------------------------------------------------------------
!
      underline = '________________________________________'// &
                  '____________________________________'
!
!     Find the longest name subject to a minimum of 12
!
      longest = 0
!
      do i = 1,nvar
!
          if (length(name(i)) > longest) then
              longest = length(name(i))
          end if
!
      end do
!
      if (longest < 12) then
          longest = 12
      end if
!
      if (.not. iread) then
          call wrtlin('    *** WARNING *** NO DATA READ')
          return
      end if
!
      call newlin
      call wrtlin('    Name'//space(1:longest)//'Levels    Type')
      call wrtlin('    '//underline(1:20+longest))
!
      do 10 i = 1,nvar
!
!-------- y-variate
          if (name(i) == yname) then
              vartype = 'YVAR'
!-------- risk-variate
          else if (name(i) == rname) then
              vartype = 'RVAR'
!-------- case-variate
          else if (nlevel == 1 .and. name(i) == cname(1)) then
              vartype = 'CASE'
!-------- case-variate
          else if (nlevel == 2 .and. name(i) == cname(1)) then
              vartype = 'CASE_1'
!-------- case-variate
          else if (nlevel == 2 .and. name(i) == cname(2)) then
              vartype = 'CASE_2'
!-------- offset variable
          else if (offlag(1) .and. name(i) == offnme(1)) then
              vartype = 'OFFSET_1'
          else if (offlag(2) .and. name(i) == offnme(2)) then
              vartype = 'OFFSET_2'
          else if (offlag(3) .and. name(i) == offnme(3)) then
              vartype = 'OFFSET_3'
!-------- explanatory variable
          else
              vartype = 'X'
          end if
!
          write (outbuf,'(3a,i2,2a)') '    ',name(i)(1:longest), &
          '      ',ilev(i),'      ',vartype
          call wrtlin(outbuf)
   10 end do
!
      call newlin
!
      return
!
      end subroutine disvar
!
!***********************************************************************
!
      subroutine disset(nm,sca,con,alp,endind,tol,est0,est1,nmeil,niter, &
                        arith,offlag,link,corr,bivar,robust, &
                        order,sig,family,trivar,univar,n1var,n2var,rho, &
                        nlevel,depend,cquad)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character endind,arith,link(3),family(3),cquad
      double precision sca(3),con,alp,tol,est0,est1,sig(3),rho(3)
      integer nm(3),nmeil,niter,n1var,n2var,nlevel
      logical offlag(3),corr,bivar,trivar,robust,univar,depend,order
!-----------------------------------------------------------------------
!     Function : Displays the settings of the model constants.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character ind1*10,ind3*8,ind10*21,ind11*3,ind12*3,ind2*10, &
                ind15*10,ind16*8,ind17*21,ind18*8,ind19*8, &
                ind20*21,ind9*3,ind25*3
      character(len=80) :: outbuf
!-----------------------------------------------------------------------
!---- left endpoint only in model
      if (endind == 'l') then
          ind1 = 'LEFT ONLY'
!---- right endpoint only in model
      else if (endind == 'r') then
          ind1 = 'RIGHT ONLY'
!---- both endpoints in model
      else if (endind == 'b') then
          ind1 = 'BOTH'
!---- no endpoints in model
      else
          ind1 = 'NONE'
      end if
!
!---- accurate arithmetic
      if (arith == 'a') then
          ind3 = 'ACCURATE'
!---- fast arithmetic
      else
          ind3 = 'FAST'
      end if
!
!---- Gaussian quadrature
      if (cquad == 'g') then
          ind2 = 'GAUSSIAN'
!---- adaptive quadrature
      else if (cquad == 'a') then
          ind2 = 'ADAPTIVE'
!---- adaptive-b quadrature
      else if (cquad == 'b') then
          ind2 = 'ADAPTIVE-B'
!---- adaptive-c quadrature
      else if (cquad == 'c') then
          ind2 = 'ADAPTIVE-C'
!---- adaptive-d quadrature
      else if (cquad == 'd') then
          ind2 = 'ADAPTIVE-D'
      end if
!
      if (family(1) == 'b') then
          ind16 = 'BINOMIAL'
      else if (family(1) == 'g') then
          ind16 = 'GAUSSIAN'
      else if (family(1) == 'p') then
          ind16 = 'POISSON'
      else if (family(1) == 'o') then
          ind16 = 'ORDERED'
      end if
!
!---- poisson: log link
      if (family(1) == 'p') then
          ind10 = 'LOG'
!---- normal: identity link
      else if (family(1) == 'g') then
          ind10 = 'IDENTITY'
!---- cll link
      else if (link(1) == 'c') then
          ind10 = 'COMPLEMENTARY LOG-LOG'
!---- probit link
      else if (link(1) == 'p') then
          ind10 = 'PROBIT'
!---- logit link
      else
          ind10 = 'LOGIT'
      end if
!
      if (.not. univar) then
!
          if (family(2) == 'b') then
              ind18 = 'BINOMIAL'
          else if (family(2) == 'g') then
              ind18 = 'GAUSSIAN'
          else if (family(2) == 'p') then
              ind18 = 'POISSON'
          else if (family(2) == 'o') then
              ind18 = 'ORDERED'
          end if
!
!-------- poisson: log link
          if (family(2) == 'p') then
              ind17 = 'LOG'
!-------- normal: identity link
          else if (family(2) == 'g') then
              ind17 = 'IDENTITY'
!-------- cll link
          else if (link(2) == 'c') then
              ind17 = 'COMPLEMENTARY LOG-LOG'
!-------- probit link
          else if (link(2) == 'p') then
              ind17 = 'PROBIT'
!-------- logit link
          else
              ind17 = 'LOGIT'
          end if
!
          if (trivar) then
!
              if (family(3) == 'b') then
                  ind19 = 'BINOMIAL'
              else if (family(3) == 'g') then
                  ind19 = 'GAUSSIAN'
              else if (family(3) == 'p') then
                  ind19 = 'POISSON'
              end if
!
              if (family(3) == 'p') then
                  ind20 = 'LOG'
              else if (family(3) == 'g') then
                  ind20 = 'IDENTITY'
              else if (link(3) == 'c') then
                  ind20 = 'COMPLEMENTARY LOG-LOG'
              else if (link(3) == 'p') then
                  ind20 = 'PROBIT'
              else
                  ind20 = 'LOGIT'
              end if
!
          end if
!
      end if
!
      if (depend) then
          ind25 = 'ON'
      else
          ind25 = 'OFF'
      end if
!
      if (offlag(1) .or. offlag(2) .or. offlag(3)) then
          ind9 = 'ON'
      else
          ind9 = 'OFF'
      end if
!
      if (univar) then
          ind15 = 'UNIVARIATE'
      else if (bivar) then
          ind15 = 'BIVARIATE'
      else
          ind15 = 'TRIVARIATE'
      end if
!
      if (corr) then
          ind11 = 'ON'
      else
          ind11 = 'OFF'
      end if
!
!      if (robust) then
!          ind12 = 'ON'
!      else
!          ind12 = 'OFF'
!      end if
!
!---- display settings
      call newlin
      call wrtlin('    Setting                    Value')
      call wrtlin('    ______________________________________')
      call newlin
      write (outbuf,'(a,g12.5)') '    Orthogonality Constant    ',alp
      call wrtlin(outbuf)
      write (outbuf,'(a,g12.5)') '    Aliasing Tolerance        ',tol
      call wrtlin(outbuf)
      write (outbuf,'(a,g12.5)') '    Convergence Criterion     ',con
      call wrtlin(outbuf)
      write (outbuf,'(a,i5)') '    Maximum Iterations         ',niter
      call wrtlin(outbuf)
      write (outbuf,'(a,i5)') '    Approximate Iterations     ',nmeil
      call wrtlin(outbuf)
      call wrtlin('    Arithmetic                 '//ind3)
      call wrtlin('    Quadrature                 '//ind2)
!
      if (univar) then
          write (outbuf,'(a,i3)') &
          '    Mass Points                  ',nm(1)
          call wrtlin(outbuf)
      else
          write (outbuf,'(a,i3)') &
          '    Mass Points [1]              ',nm(1)
          call wrtlin(outbuf)
          write (outbuf,'(a,i3)') &
          '    Mass Points [2]              ',nm(2)
          call wrtlin(outbuf)
!
          if (trivar) then
              write (outbuf,'(a,i3)') &
              '    Mass Points [3]              ',nm(3)
              call wrtlin(outbuf)
          end if
!
      end if
!
      call newlin
!
!      if (.not. order) then
          call wrtlin('    Model                      '//ind15)
!
          if (univar) then
              call wrtlin('    Family                     '//ind16)
              call wrtlin('    Link                       '//ind10)
          else
              call wrtlin('    Family [1]                 '//ind16)
              call wrtlin('    Link [1]                   '//ind10)
              call wrtlin('    Family [2]                 '//ind18)
              call wrtlin('    Link [2]                   '//ind17)
!
              if (trivar) then
                  call wrtlin('    Family [3]                 '//ind19)
                  call wrtlin('    Link [3]                   '//ind20)
              end if
!
              call wrtlin('    Correlation                '//ind11)
!
!              write (outbuf,'(a,i4)') &
!              '    Number of variables [1]       ',n1var
!              call wrtlin(outbuf)
!
!              if (trivar) then
!                  write (outbuf,'(a,i4)') &
!                  '    Number of variables [2]       ',n2var
!                  call wrtlin(outbuf)
!              end if
!
          end if
!
          call newlin
!      end if
!
      call wrtlin('    Offset                     '//ind9)
      call wrtlin('    Dependent                  '//ind25)
      call wrtlin('    Endpoints                  '//ind1)
!
      if (endind == 'b' .or. endind == 'l') then
          write (outbuf,'(a,g12.5)') &
          '    Initial Left Endpoint     ',est0
          call wrtlin(outbuf)
      end if
!
      if (endind == 'b' .or. endind == 'r') then
          write (outbuf,'(a,g12.5)') &
          '    Initial Right Endpoint    ',est1
          call wrtlin(outbuf)
      end if
!
!      call wrtlin('    Robust                     '//ind12)
      call newlin
!
      if (univar) then
!
          if (family(1) == 'g') then
              write (outbuf,'(a,g12.5)') &
              '    Initial Sigma             ',sig(1)
              call wrtlin(outbuf)
          end if
!
          write (outbuf,'(a,g12.5)') '    Initial Scale             ', &
          sca(1)
          call wrtlin(outbuf)
!
          if (nlevel == 2) then
              write (outbuf,'(a,g12.5)') &
              '    Initial Level 2 Scale     ',sca(2)
              call wrtlin(outbuf)
          end if
!
      else
!
          if (family(1) == 'g') then
              write (outbuf,'(a,g12.5)') &
              '    Initial Sigma [1]         ',sig(1)
              call wrtlin(outbuf)
          end if
!
          if (family(2) == 'g') then
              write (outbuf,'(a,g12.5)') &
              '    Initial Sigma [2]         ',sig(2)
              call wrtlin(outbuf)
          end if
!
          if (trivar .and. family(3) == 'g') then
              write (outbuf,'(a,g12.5)') &
              '    Initial Sigma [3]         ',sig(3)
              call wrtlin(outbuf)
          end if
!
          write (outbuf,'(a,g12.5)') '    Initial Scale [1]         ', &
          sca(1)
          call wrtlin(outbuf)
          write (outbuf,'(a,g12.5)') '    Initial Scale [2]         ', &
          sca(2)
          call wrtlin(outbuf)
!
          if (.not. trivar .and. corr) then
              write (outbuf,'(a,g12.5)') &
              '    Initial Correlation       ',rho(1)
              call wrtlin(outbuf)
          end if
!
          if (trivar) then
              write (outbuf,'(a,g12.5)') &
              '    Initial Scale [3]         ',sca(3)
              call wrtlin(outbuf)
!
              if (corr) then
                  write (outbuf,'(a,g12.5)') &
                  '    Initial Correlation [1,2] ',rho(1)
                  call wrtlin(outbuf)
                  write (outbuf,'(a,g12.5)') &
                  '    Initial Correlation [1,3] ',rho(2)
                  call wrtlin(outbuf)
                  write (outbuf,'(a,g12.5)') &
                  '    Initial Correlation [2,3] ',rho(3)
                  call wrtlin(outbuf)
              end if
!
          end if
!
      end if
!
      call newlin
!
      return
!
      end subroutine disset
!
!***********************************************************************
!
      subroutine dismod(xll,ilevdf,iredf,ienddf,ilfit,mnames,ni,cname, &
                        yname,nsub,nmes,itotdf,endind,link,corr, &
                        bivar,robust,order,inormdf,family,trivar,univar, &
                        offlag,offnme,nlevel,depend,icutdf)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character endind,link(3),family(3)
      character(len=50) :: mnames(maxvar),cname(2),yname,offnme(3)
      double precision xll
      integer ilevdf,iredf,ienddf,ni,nsub(2),nmes,ilfit,icutdf, &
              itotdf,inormdf,nlevel
      logical corr,bivar,robust,trivar,univar,offlag(3),depend,order
!-----------------------------------------------------------------------
!     Function : Displays the model variables, deviance and residual df.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      integer i
!-----------------------------------------------------------------------
      if (ilfit == -1) then
          call wrtlin( &
          '    *** WARNING *** THE CURRENT MODEL IS INVALID')
          return
      end if
!
      if (ilfit == 0 .and. (offlag(1) .or. offlag(2) .or. offlag(3))) &
      then
          call newlin
          call wrtlin( &
          '    X-vars            Y-var             Case-var          '// &
          'Offset')
          call wrtlin('    __________________________________________'// &
          '________________________________')
          call wrtlin('    '//mnames(1)(1:12)//'      '//yname(1:12)// &
          '      '//cname(1)(1:12)//'      '//offnme(1)(1:20))
      else if (ilfit == 0) then
          call newlin
          call wrtlin( &
          '    X-vars            Y-var             Case-var')
          call wrtlin( &
          '    ________________________________________________')
          call wrtlin('    '//mnames(1)(1:12)//'      '//yname(1:12)// &
          '      '//cname(1)(1:12))
      else if (offlag(1) .or. offlag(2) .or. offlag(3)) then
          call newlin
          call wrtlin('    X-vars            Y-var             Offset')
          call wrtlin('    __________________________________________'// &
          '______________')
          call wrtlin('    '//mnames(1)(1:12)//'      '//yname(1:12)// &
          '      '// offnme(1)(1:20))
      else
          call newlin
          call wrtlin('    X-vars            Y-var')
          call wrtlin('    ______________________________')
          call wrtlin('    '//mnames(1)(1:12)//'      '//yname(1:12))
      end if
!
      do 10 i = 2,ni
!
          if (i == 3 .and. ilfit == 0 .and. offlag(3)) then
              call wrtlin('    '//mnames(i)(1:12)// &
              '                       '//' '//'            '//'      '// &
              offnme(3)(1:20))
          else if (i == 3 .and. offlag(3)) then
              call wrtlin('    '//mnames(i)(1:12)// &
              '                       '//' '//offnme(3)(1:20))
          else if (i == 2 .and. ilfit == 0 .and. offlag(2)) then
              call wrtlin('    '//mnames(i)(1:12)// &
              '                       '//' '//'            '//'      '// &
              offnme(2)(1:20))
          else if (i == 2 .and. offlag(2)) then
              call wrtlin('    '//mnames(i)(1:12)// &
              '                       '//' '//offnme(2)(1:20))
          else if (ilfit == 0 .and. nlevel == 2 .and. i == 2) then
              call wrtlin('    '//mnames(i)(1:12)// &
              '                       '//' '//cname(2)(1:12))
          else
              call wrtlin('    '//mnames(i)(1:12))
          end if
!
 10   end do
!
      call newlin
!
      if (univar) then
          call wrtlin('    Univariate model')
      else if (bivar) then
!
          if (corr) then
              call wrtlin('    Correlated bivariate model')
          else
              call wrtlin('    Independent bivariate model')
          end if
!
          call newlin
      else if (trivar) then
!
          if (corr) then
              call wrtlin('    Correlated trivariate model')
          else
              call wrtlin('    Independent trivariate model')
          end if
!
          call newlin
      end if
!
      if (nlevel == 2) then
          call wrtln2('    Multilevel')
      else if (depend) then
          call wrtln2('    Dependent')
      else
          call wrtln2('    Standard')
      end if
!
      if (family(1) == 'b') then
!
          if (link(1) == 'c') then
              call wrtln2(' complementary log-log')
          else if (link(1) == 'p') then
              call wrtln2(' probit')
          else
              call wrtln2(' logit')
          end if
!
      else if (family(1) == 'p') then
          call wrtln2(' Poisson')
      else if (family(1) == 'g') then
          call wrtln2(' linear')
      else if (family(1) == 'o') then
          call wrtln2(' ordered')
!
          if (link(1) == 'p') then
              call wrtln2(' probit')
          else
              call wrtln2(' logit')
          end if
!
      end if
!
      if (.not. univar) then
!
          if (family(2) == 'b') then
!
              if (link(2) == 'c') then
                  call wrtln2('/complementary log-log')
              else if (link(2) == 'p') then
                  call wrtln2('/probit')
              else
                  call wrtln2('/logit')
              end if
!
          else if (family(2) == 'p') then
              call wrtln2('/Poisson')
          else if (family(2) == 'g') then
              call wrtln2('/linear')
          else if (family(2) == 'o') then
              call wrtln2('/ordered')
!
              if (link(1) == 'p') then
                  call wrtln2(' probit')
              else
                  call wrtln2(' logit')
              end if
!
          end if
!
          if (trivar) then
!
              if (family(3) == 'b') then
!
                  if (link(3) == 'c') then
                      call wrtln2('/complementary log-log')
                  else if (link(3) == 'p') then
                      call wrtln2('/probit')
                  else
                      call wrtln2('/logit')
                  end if
!
              else if (family(3) == 'p') then
                  call wrtln2('/Poisson')
              else if (family(3) == 'g') then
                  call wrtln2('/linear')
              end if
!
          end if
!
      end if
!
      if (ilfit == 0) then
          call newlin
          call wrtln2('    Gaussian random effects')
!
          if (endind == 'b') then
              call wrtln2(', with endpoints')
          else if (endind == 'l') then
              call wrtln2(', with left endpoint')
          else if (endind == 'r') then
              call wrtln2(', with right endpoint')
          end if
!
      else if (ilfit == 3) then
          call newlin
          call wrtln2('    Fixed effects')
      end if
!
      call newlin
!
      if (ilfit == 1 .and. robust) then
          call wrtlin('    Robust standard errors')
          call newlin
      end if
!
      call newlin
      write (outbuf,'(a,i7)') &
      '    Number of observations             = ',nmes
      call wrtlin(outbuf)
!
      if (ilfit == 0 .or. ilfit == 3) then
!
          if (nlevel == 1) then
              write (outbuf,'(a,i7)') &
              '    Number of cases                    = ',nsub(1)
              call wrtlin(outbuf)
          else
              write (outbuf,'(a,i7)') &
              '    Number of level 2 cases            = ',nsub(1)
              call wrtlin(outbuf)
              write (outbuf,'(a,i7)') &
              '    Number of level 3 cases            = ',nsub(2)
              call wrtlin(outbuf)
          end if
!
      end if
!
      call newlin
      write (outbuf,'(a,i5)') '    X-var df           = ',ilevdf
      call wrtlin(outbuf)
!
      if (order) then
          write (outbuf,'(a,i5)') '    Cutpoint df        = ',icutdf
          call wrtlin(outbuf)
      end if
!
      if (ilfit == 3) then
          go to 900
      end if
!
      if (ilfit == 1 .and. (family(1) == 'g' .or. family(2) == 'g' &
      .or. family(3) == 'g')) then
          write (outbuf,'(a,i5)') '    Sigma df           = ',inormdf
          call wrtlin(outbuf)
      end if
!
      if (ilfit == 0) then
!
          if (family(1) == 'g' .or. family(2) == 'g' .or. &
          family(3) == 'g') then
              write (outbuf,'(a,i5)') '    Sigma df           = ', &
              inormdf
              call wrtlin(outbuf)
          end if
!
          write (outbuf,'(a,i5)') '    Scale df           = ',iredf
          call wrtlin(outbuf)
!
          if (endind /= 'n') then
              write (outbuf,'(a,i5)') '    Endpoint df        = ',ienddf
              call wrtlin(outbuf)
          end if
!
      end if
!
      call newlin
      write (outbuf,'(a,g18.8,a,i7,a)') '    Log likelihood = ',xll, &
      ' on ',itotdf,' residual degrees of freedom'
      call wrtlin(outbuf)
!
  900 call newlin
!
      return
!
      end subroutine dismod
!
!***********************************************************************
!
      subroutine disest(beta,cov,alias,ipos,mnames,ni,nvar,ilfit,npar, &
                        ilev,iend,endind,corr,robust,order,ncat,sig_e, &
                        family,trivar,univar,nlevel,depend,eqscale)
!-----------------------------------------------------------------------
      use estimatesmod
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character mnames(maxvar)*50,endind,family(3)
      integer npar,nlevel,ipos(maxvar),ni,nvar,ilfit,ilev(maxvar), &
              iend(2),ncat(3),length,longest
      double precision beta(npar),cov(npar*(npar+1)/2),alias(npar), &
             sig_e(3)
      logical corr,robust,trivar,univar,depend,eqscale,order
!-----------------------------------------------------------------------
!     Function : Displays the parameter estimates and standard errors.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      character partxt(12)*50,text(3)*11,endtxt(2)*11, &
                underline*76,space*76
      double precision prob,z
      integer i,j,k,m,index,l,nend,istat
!-----------------------------------------------------------------------
      data partxt /'sigma      ','scale      ','corr       ', &
                   'sigma1     ','sigma2     ','sigma3     ', &
                   'scale1     ','scale2     ','scale3     ', &
                   'corr12     ','corr13     ','corr23     '/
      data endtxt /'endpoint 0 ','endpoint 1 '/
      data text /'ALIASED [I]','ALIASED [E]','FIXED      '/
      data space /' '/
!-----------------------------------------------------------------------
!     Allocate space for the model arrays if they have not already been
!     allocated
!
      if (.not. allocated(model_names)) then
          allocate(model_names(maxpar), stat=istat)
      end if
!
      if (.not. allocated(model_estimates)) then
          allocate(model_estimates(maxpar), stat=istat)
      end if
!
      if (.not. allocated(model_errors)) then
          allocate(model_errors(maxpar), stat=istat)
      end if
!
      estimate_count = 0
!
      underline = '________________________________________'// &
                  '____________________________________'
!
!     Find the longest name subject to a minimum of 12 and maximum of 20
!
      longest = 0
!
      do j = 1,ni
!
          if (length(mnames(j)) > longest) then
              longest = length(mnames(j))
          end if
!
      end do
!
      if (longest < 12) then
          longest = 12
      end if
!
      if (longest > 20) then
          longest = 20
      end if
!
      if (ilfit == -1) then
          call wrtlin( &
          '    *** WARNING *** THE CURRENT MODEL IS INVALID')
          return
      end if
!
      call newlin
!
      if (ilfit == 1 .and. robust) then
         call wrtlin('    Parameter '//space(1:longest)//' Estimate'// &
                     '         Rob. S.E.        Z-score')
      else
         call wrtlin('    Parameter '//space(1:longest)//' Estimate'// &
                     '         Std. Err.        Z-score')
      end if
!
      call wrtlin( '    '//underline(1:longest+56))
!
      index = 0
!
!---- outer loop looks for jth position in model. Inner loop finds the
!---- variable that is in the jth position and prints it out.
      do 120 j = 1,ni
!
          do 110 k = 1,nvar
!
              if (ipos(k) == index+1) then
                  go to 115
              end if
!
  110     end do
!
  115     index = index + ilev(k)
!
          do 118 l = index,index-ilev(k)+1,-1
!
!------------ display parameter estimates and standard errors
              if (idnint(alias(l)) == 0) then
                  z = beta(l)/sqrt(cov((l+1)*l/2))
!
                  if (ilev(k) == 1) then
                      write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
                      mnames(j)(1:longest),beta(l),sqrt(cov((l+1)*l/2)), &
                      z
                      call wrtlin(outbuf)
!                     Update the arrays wherever the estimate etc. is
!                     output
                      estimate_count = estimate_count+1
                      model_names(estimate_count) = mnames(j)
                      model_estimates(estimate_count) = beta(l)
                      model_errors(estimate_count) = &
                      sqrt(cov((l+1)*l/2))
                  else
                      write (outbuf,'(4x,2a,i4,a,1x,g15.5,2(2x,g15.5))') &
                      mnames(j)(1:longest),'(',index-l+1,')',beta(l), &
                      sqrt(cov((l+1)*l/2)),z
                      call wrtlin(outbuf)
                      estimate_count = estimate_count+1
                      model_names(estimate_count) = mnames(j)
                      model_estimates(estimate_count) = beta(l)
                      model_errors(estimate_count) = &
                      sqrt(cov((l+1)*l/2))
                  end if
!
!------------ display estimates for aliased factor levels
              else
                  m = 1
!
                  if (idnint(alias(l)) == 1) then
                      m = 2
                  end if
!
                  if (ilev(k) == 1) then
                      write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
                      mnames(j)(1:longest),beta(l),text(m)
                      call wrtlin(outbuf)
                      estimate_count = estimate_count+1
                      model_names(estimate_count) = mnames(j)
                      model_estimates(estimate_count) = beta(l)
                      model_errors(estimate_count) = 0d0
                  else
                      write (outbuf,'(4x,2a,i4,a,1x,g15.5,6x,a)') &
                      mnames(j)(1:longest),'(',index-l+1,')',beta(l), &
                      text(m)
                      call wrtlin(outbuf)
                      estimate_count = estimate_count+1
                      model_names(estimate_count) = mnames(j)
                      model_estimates(estimate_count) = beta(l)
                      model_errors(estimate_count) = 0d0
                  end if
!
              end if
!
  118     end do
!
  120 end do
!
      if (order) then
!
          do 50 i = 1,ncat(1)-1
!
              if (idnint(alias(index+i)) == 0) then
                  z = beta(index+i)/sqrt(cov((index+i+1)*(index+i)/2))
!
                  if (univar .and. .not. depend) then
                      write (outbuf,'(4x,a,i1,a,3x,g15.5,2(2x,g15.5))') &
                      'cut',i,space(1:longest),beta(index+i), &
                      sqrt(cov((index+i+1)*(index+i)/2)),z
                  else
                      write (outbuf,'(4x,a,i1,a,1x,g15.5,2(2x,g15.5))') &
                      'cut1_',i,space(1:longest),beta(index+i), &
                      sqrt(cov((index+i+1)*(index+i)/2)),z
                  end if
!
              else
!
                  if (univar .and. .not. depend) then
                      write (outbuf,'(4x,a,i1,a,3x,g15.5,2x,g15.5)') &
                      'cut',i,space(1:longest),beta(index+i),text(2)
                  else
                      write (outbuf,'(4x,a,i1,a,1x,g15.5,2x,g15.5)') &
                      'cut1_',i,space(1:longest),beta(index+i),text(2)
                  end if
!
              end if
!
              call wrtlin(outbuf)
              estimate_count = estimate_count+1
!
              if (univar) then
                  model_names(estimate_count) = 'cut'//char(i+48)
              else
                  model_names(estimate_count) = 'cut1_'//char(i+48)
              end if
!
              model_estimates(estimate_count) = beta(index+i)
              model_errors(estimate_count) = &
              sqrt(cov((index+i+1)*(index+i)/2))
   50     end do
!
          index = index + ncat(1) - 1
!
          if (.not. univar .or. depend) then
!
              do 51 i = 1,ncat(2)-1
!
                  if (idnint(alias(index+i)) == 0) then
                      z = beta(index+i)/ &
                      sqrt(cov((index+i+1)*(index+i)/2))
                      write (outbuf,'(4x,a,i1,a,1x,g15.5,2(2x,g15.5))') &
                      'cut2_',i,space(1:longest),beta(index+i), &
                      sqrt(cov((index+i+1)*(index+i)/2)),z
                  else
                      write (outbuf,'(4x,a,i1,a,1x,g15.5,2x,g15.5)') &
                      'cut2_',i,space(1:longest),beta(index+i),text(2)
                  end if
!
                  call wrtlin(outbuf)
                  estimate_count = estimate_count+1
                  model_names(estimate_count) = 'cut2_'//char(i+48)
                  model_estimates(estimate_count) = beta(index+i)
                  model_errors(estimate_count) = &
                  sqrt(cov((index+i+1)*(index+i)/2))
   51         end do
!
              index = index + ncat(2) - 1
!
              if (trivar) then
!
                  do 52 i = 1,ncat(3)-1
!
                      if (idnint(alias(index+i)) == 0) then
                          z = beta(index+i)/ &
                          sqrt(cov((index+i+1)*(index+i)/2))
                          write (outbuf, &
                          '(4x,a,i1,a,1x,g15.5,2(2x,g15.5))') &
                          'cut3_',i,space(1:longest),beta(index+i), &
                          sqrt(cov((index+i+1)*(index+i)/2)),z
                      else
                          write (outbuf,'(4x,a,i1,a,1x,g15.5,2x,g15.5)') &
                          'cut3_',i,space(1:longest),beta(index+i), &
                          text(2)
                      end if
!
                      call wrtlin(outbuf)
                      estimate_count = estimate_count+1
                      model_names(estimate_count) = 'cut3_'//char(i+48)
                      model_estimates(estimate_count) = beta(index+i)
                      model_errors(estimate_count) = &
                      sqrt(cov((index+i+1)*(index+i)/2))
   52             end do
!
                  index = index + ncat(3) - 1
              end if
!
          end if
!
      end if
!
      if (ilfit == 1 .or. ilfit == 3) then
!
          if (family(1) == 'g') then
              m = 1
!
              if (.not. univar) then
                  m = 4
              end if
!
              write (outbuf,'(4x,a,7x,g15.5)') partxt(m)(1:longest), &
              sig_e(1)
              call wrtlin(outbuf)
          end if
!
          if (.not. univar .and. family(2) == 'g') then
              write (outbuf,'(4x,a,7x,g15.5)') partxt(5)(1:longest), &
              sig_e(2)
              call wrtlin(outbuf)
          end if
!
          if (trivar .and. family(3) == 'g') then
              write (outbuf,'(4x,a,7x,g15.5)') partxt(6)(1:longest), &
              sig_e(3)
              call wrtlin(outbuf)
          end if
!
          go to 200
      end if
!
      if (family(1) == 'g') then
          index = index+1
          m = 1
!
          if (.not. univar) then
              m = 4
          end if
!
          if (idnint(alias(index)) == 0) then
              z = beta(index)/sqrt(cov((index+1)*index/2))
              write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
              partxt(m)(1:longest),beta(index), &
              sqrt(cov((index+1)*index/2)),z
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = partxt(m)
              model_estimates(estimate_count) = beta(index)
              model_errors(estimate_count) = &
              sqrt(cov((index+1)*index/2))
          else
              write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
              partxt(m)(1:longest),beta(index),text(2)
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = partxt(m)
              model_estimates(estimate_count) = beta(index)
              model_errors(estimate_count) = 0d0
          end if
!
      end if
!
      if (.not. univar .and. family(2) == 'g') then
          index = index+1
!
          if (idnint(alias(index)) == 0) then
              z = beta(index)/sqrt(cov((index+1)*index/2))
              write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
              partxt(5)(1:longest),beta(index), &
              sqrt(cov((index+1)*index/2)),z
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = partxt(5)
              model_estimates(estimate_count) = beta(index)
              model_errors(estimate_count) = &
              sqrt(cov((index+1)*index/2))
          else
              write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
              partxt(5)(1:longest),beta(index),text(2)
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = partxt(5)
              model_estimates(estimate_count) = beta(index)
              model_errors(estimate_count) = 0d0
          end if
!
      end if
!
      if (trivar .and. family(3) == 'g') then
          index = index+1
!
          if (idnint(alias(index)) == 0) then
              z = beta(index)/sqrt(cov((index+1)*index/2))
              write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
              partxt(6)(1:longest),beta(index), &
              sqrt(cov((index+1)*index/2)),z
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = partxt(6)
              model_estimates(estimate_count) = beta(index)
              model_errors(estimate_count) = &
              sqrt(cov((index+1)*index/2))
          else
              write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
              partxt(6)(1:longest),beta(index),text(2)
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = partxt(6)
              model_estimates(estimate_count) = beta(index)
              model_errors(estimate_count) = 0d0
          end if
!
      end if
!
      m = 2
!
      if ((.not. univar .and. .not. eqscale) .or. depend) then
          m = 7
      end if
!
      if (ilfit == 2) then
          go to 200
      end if
!
      if (.not. eqscale) then
          index = index+1
!
          if (nlevel == 1) then
!
!------------ display scale parameter estimate and standard error
              if (idnint(alias(index)) == 0) then
                  z = beta(index)/sqrt(cov((index+1)*index/2))
                  write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
                  partxt(m)(1:longest),beta(index), &
                  sqrt(cov((index+1)*index/2)),z
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = partxt(m)
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = &
                  sqrt(cov((index+1)*index/2))
!------------ display estimate for aliased scale parameter
              else if (ilfit == 0) then
                  write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
                  partxt(m)(1:longest),beta(index),text(2)
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = partxt(m)
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = 0d0
              end if
!
              if (depend) then
                  index = index+1
!
!---------------- display scale parameter estimate and standard error
                  if (idnint(alias(index)) == 0) then
                      z = beta(index)/sqrt(cov((index+1)*index/2))
                      write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
                      partxt(8)(1:longest),beta(index), &
                      sqrt(cov((index+1)*index/2)),z
                      call wrtlin(outbuf)
                      estimate_count = estimate_count + 1
                      model_names(estimate_count) = partxt(8)
                      model_estimates(estimate_count) = beta(index)
                      model_errors(estimate_count) = &
                      sqrt(cov((index+1)*index/2))
!---------------- display estimate for aliased scale parameter
                  else if (ilfit == 0) then
                      write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
                      partxt(8)(1:longest),beta(index),text(2)
                      call wrtlin(outbuf)
                      estimate_count = estimate_count + 1
                      model_names(estimate_count) = partxt(8)
                      model_estimates(estimate_count) = beta(index)
                      model_errors(estimate_count) = 0d0
                  end if
!
              end if
!
          else
!
!------------ display scale parameter estimate and standard error
              if (idnint(alias(index)) == 0) then
                  z = beta(index)/sqrt(cov((index+1)*index/2))
                  write (outbuf,'(4x,2a,g15.5,2(2x,g15.5))') &
                  'scale2 ',space(1:longest),beta(index), &
                  sqrt(cov((index+1)*index/2)),z
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = 'scale2'
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = &
                  sqrt(cov((index+1)*index/2))
!------------ display estimate for aliased scale parameter
              else if (ilfit == 0) then
                  write (outbuf,'(4x,2a,g15.5,6x,a)') &
                  'scale2 ',space(1:longest),beta(index),text(2)
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = 'scale2'
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = 0d0
              end if
!
!------------ display scale parameter estimate and standard error
              if (idnint(alias(index+1)) == 0) then
                  z = beta(index+1)/sqrt(cov((index+2)*(index+1)/2))
                  write (outbuf,'(4x,2a,g15.5,2(2x,g15.5))') &
                  'scale3 ',space(1:longest),beta(index+1), &
                  sqrt(cov((index+2)*(index+1)/2)),z
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = 'scale3'
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = &
                  sqrt(cov((index+1)*index/2))
!------------ display estimate for aliased scale parameter
              else if (ilfit == 0) then
                  write (outbuf,'(4x,2a,g15.5,6x,a)') &
                  'scale3 ',space(1:longest),beta(index+1),text(2)
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = 'scale3'
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = 0d0
              end if
!
          end if
!
      end if
!
      m = 8
!
      if (eqscale) then
          m = 2
      end if
!
      if (.not. univar) then
!
          index = index+1
!
          if (idnint(alias(index)) == 0) then
              z = beta(index)/sqrt(cov((index+1)*index/2))
              write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
              partxt(m)(1:longest),beta(index), &
              sqrt(cov((index+1)*index/2)),z
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = partxt(m)
              model_estimates(estimate_count) = beta(index)
              model_errors(estimate_count) = &
              sqrt(cov((index+1)*index/2))
!-------- display estimate for aliased scale parameter
          else if (ilfit == 0) then
              write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
              partxt(m)(1:longest),beta(index),text(2)
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = partxt(m)
              model_estimates(estimate_count) = beta(index)
              model_errors(estimate_count) = 0d0
          end if
!
          if (trivar) then
              index = index+1
!
              if (idnint(alias(index)) == 0) then
                  z = beta(index)/sqrt(cov((index+1)*index/2))
                  write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
                  partxt(9)(1:longest),beta(index), &
                  sqrt(cov((index+1)*index/2)),z
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = partxt(9)
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = &
                  sqrt(cov((index+1)*index/2))
              else if (ilfit == 0) then
                  write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
                  partxt(9)(1:longest),beta(index),text(2)
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = partxt(9)
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = 0d0
              end if
!
          end if
!
          if (corr) then
              index = index+1
              m = 3
!
              if (trivar) then
                  m = 10
              end if
!
              if (idnint(alias(index)) == 0) then
                  z = beta(index)/sqrt(cov((index+1)*index/2))
                  write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
                  partxt(m)(1:longest),beta(index), &
                  sqrt(cov((index+1)*index/2)),z
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = partxt(m)
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = &
                  sqrt(cov((index+1)*index/2))
!------------ display estimate for aliased scale parameter
              else if (ilfit == 0) then
                  write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
                  partxt(m)(1:longest),beta(index),text(2)
                  call wrtlin(outbuf)
                  estimate_count = estimate_count + 1
                  model_names(estimate_count) = partxt(m)
                  model_estimates(estimate_count) = beta(index)
                  model_errors(estimate_count) = 0d0
              end if
!
              if (trivar) then
                  m = 10
 1250             m = m+1
                  index = index+1
!
                  if (idnint(alias(index)) == 0) then
                     z = beta(index)/sqrt(cov((index+1)*index/2))
                     write (outbuf,'(4x,a,7x,g15.5,2(2x,g15.5))') &
                     partxt(m)(1:longest),beta(index), &
                     sqrt(cov((index+1)*index/2)),z
                     call wrtlin(outbuf)
                     estimate_count = estimate_count + 1
                     model_names(estimate_count) = partxt(m)
                     model_estimates(estimate_count) = beta(index)
                     model_errors(estimate_count) = &
                     sqrt(cov((index+1)*index/2))
                  else if (ilfit == 0) then
                     write (outbuf,'(4x,a,7x,g15.5,6x,a)') &
                     partxt(m)(1:longest),beta(index),text(2)
                     call wrtlin(outbuf)
                     estimate_count = estimate_count + 1
                     model_names(estimate_count) = partxt(m)
                     model_estimates(estimate_count) = beta(index)
                     model_errors(estimate_count) = 0d0
                  end if
!
                  if (m == 11) then
                      go to 1250
                  end if
!
              end if
!
          end if
!
      end if
!
      index = index+1
!
      if (endind == 'n') then
          go to 200
      end if
!
      nend = 1
!
      if (endind == 'b') then
          nend = 2
      end if
!
      call wrtlin(space(1:61)//'PROBABILITY')
      call wrtlin(space(1:61)//underline(1:11))
!
      do 35 j = index,index-1+nend
!
!-------- end-point probability
          if (endind == 'b') then
              prob = beta(j)/(1 + beta(index) + beta(index+1))
          else
              prob = beta(j)/(1 + beta(j))
          end if
!
          if (endind == 'r') then
              index = index-1
          end if
!
!-------- display estimate for fixed end-point
          if (iend(j-index+1) == 0) then
              prob = 0
              write (outbuf,'(4x,a,8x,f11.5,10x,a,2x,f11.5)') &
              endtxt(j-index+1),beta(j),text(3),prob
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = endtxt(j-index+1)
              model_estimates(estimate_count) = beta(j)
              model_errors(estimate_count) = 0d0
!-------- display estimate for aliased end-point
          else if (idnint(alias(j)) /= 0) then
              prob = 0
              write (outbuf,'(4x,a,8x,f11.5,10x,a,2x,f11.5)') &
              endtxt(j-index+1),beta(j),text(2),prob
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = endtxt(j-index+1)
              model_estimates(estimate_count) = beta(j)
              model_errors(estimate_count) = 0d0
!-------- display estimates and standard errors for free end-points
          else
              write (outbuf,'(4x,a,8x,g15.5,2(2x,g15.5))') &
              endtxt(j-index+1),beta(j),sqrt(cov((j+1)*j/2)),prob
              call wrtlin(outbuf)
              estimate_count = estimate_count + 1
              model_names(estimate_count) = endtxt(j-index+1)
              model_estimates(estimate_count) = beta(j)
              model_errors(estimate_count) = sqrt(cov((j+1)*j/2))
          end if
!
   35 end do
!
  200 continue
!
      call newlin
!
      return
!
      end subroutine disest
!
!***********************************************************************
!
      subroutine robind(scores,nmes,nest,cov,ncov,matrob)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer nmes,nest,ncov
      double precision scores(nmes,nest),cov(ncov),matrob(ncov)
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      double precision colvec(nest,1),rowvec(1,nest),z(1), &
             scomat(nest,nest),scosum(nest,nest),covcov(nest,nest), &
             matrix1(nest,nest),matrix2(nest,nest)
      integer i,j,n,p,m,iz,opt,j1,j2,k,ifail
!-----------------------------------------------------------------------
      do 500 j1 = 1,nest
!
          do 400 j2 = 1,nest
              scosum(j1,j2) = 0
  400     end do
!
  500 end do
!
      n = nest
      p = nest
      m = 1
      iz = 1
      opt = 1
!
      do 4000 i = 1,nmes
!
          do 1000 j = 1,nest
              colvec(j,1) = scores(i,j)
              rowvec(1,j) = scores(i,j)
 1000     end do
!
          ifail = 0
!
!---------***********
!          CALL F01CKF(scomat,colvec,rowvec,n,p,m,z,iz,opt,ifail)
!---------***********
!
          if (ifail /= 0) then
              write (outbuf,'(a,i1)') 'f01ckf[1]: ifail = ',ifail
              call wrtlin(outbuf)
          end if
!
          do 3000 j1 = 1,nest
!
              do 2000 j2 = 1,nest
                  scosum(j1,j2) = &
                  scosum(j1,j2) + scomat(j1,j2)*nmes/(nmes-1)
 2000         end do
!
 3000     end do
!
 4000 end do
!
      k = 0
!
      do 6000 j1 = 1,nest
!
          do 5000 j2 = 1,j1
              k = k+1
              covcov(j1,j2) = cov(k)
!
              if (j1 /= j2) then
                  covcov(j2,j1) = cov(k)
              end if
!
 5000     end do
!
 6000 end do
!
      n = nest
      p = nest
      m = nest
      iz = 1
      opt = 1
      ifail = 0
!
!-----***********
!      CALL F01CKF(matrix1,covcov,scosum,n,p,m,z,iz,opt,ifail)
!-----***********
!
      if (ifail /= 0) then
          write (outbuf,'(a,i1)') 'f01ckf[2]: ifail = ',ifail
          call wrtlin(outbuf)
      end if
!
      n = nest
      p = nest
      m = nest
      iz = 1
      opt = 1
      ifail = 0
!
!-----***********
!      CALL F01CKF(matrix2,matrix1,covcov,n,p,m,z,iz,opt,ifail)
!-----***********
!
      if (ifail /= 0) then
          write (outbuf,'(a,i1)') 'f01ckf[3]: ifail = ',ifail
          call wrtlin(outbuf)
      end if
!
      k = 0
!
      do 8000 j1 = 1,nest
!
          do 7000 j2 = 1,nest
!
              if (j1 >= j2) then
                  k = k+1
                  matrob(k) = matrix2(j1,j2)
              end if
!
 7000     end do
!
 8000 end do
!
      return
!
      end subroutine robind
!
!***********************************************************************
!
      subroutine robfit(scores,nmes,nest,cov,ncov,matrob)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer nmes,nest,ncov
      double precision scores(nmes,nest),cov(ncov),matrob(ncov)
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      double precision colvec(nest,1),rowvec(1,nest),z(1), &
             scomat(nest,nest),scosum(nest,nest),covcov(nest,nest), &
             matrix1(nest,nest),matrix2(nest,nest)
      integer i,j,n,p,m,iz,opt,j1,j2,k,ifail
!-----------------------------------------------------------------------
      do 500 j1 = 1,nest
!
          do 400 j2 = 1,nest
              scosum(j1,j2) = 0
  400     end do
!
  500 end do
!
      n = nest
      p = nest
      m = 1
      iz = 1
      opt = 1
!
      do 4000 i = 1,nmes
!
          do 1000 j = 1,nest
              colvec(j,1) = scores(i,j)
              rowvec(1,j) = scores(i,j)
 1000     end do
!
          ifail = 0
!
!---------***********
!          CALL F01CKF(scomat,colvec,rowvec,n,p,m,z,iz,opt,ifail)
!---------***********
!
          if (ifail /= 0) then
              write (outbuf,'(a,i1)') 'f01ckf[1]: ifail = ',ifail
              call wrtlin(outbuf)
          end if
!
          do 3000 j1 = 1,nest
!
              do 2000 j2 = 1,nest
                  scosum(j1,j2) = &
                  scosum(j1,j2) + scomat(j1,j2)*nmes/(nmes-1)
 2000         end do
!
 3000     end do
!
 4000 end do
!
      k = 0
!
      do 6000 j1 = 1,nest
!
          do 5000 j2 = 1,j1
              k = k+1
              covcov(j1,j2) = cov(k)
!
              if (j1 /= j2) then
                  covcov(j2,j1) = cov(k)
              end if
!
 5000     end do
!
 6000 end do
!
      n = nest
      p = nest
      m = nest
      iz = 1
      opt = 1
      ifail = 0
!
!-----***********
!      CALL F01CKF(matrix1,covcov,scosum,n,p,m,z,iz,opt,ifail)
!-----***********
!
      if (ifail /= 0) then
          write (outbuf,'(a,i1)') 'f01ckf[2]: ifail = ',ifail
          call wrtlin(outbuf)
      end if
!
      n = nest
      p = nest
      m = nest
      iz = 1
      opt = 1
      ifail = 0
!
!-----***********
!      CALL F01CKF(matrix2,matrix1,covcov,n,p,m,z,iz,opt,ifail)
!-----***********
!
      if (ifail /= 0) then
          write (outbuf,'(a,i1)') 'f01ckf[3]: ifail = ',ifail
          call wrtlin(outbuf)
      end if
!
      k = 0
!
      do 8000 j1 = 1,nest
!
          do 7000 j2 = 1,nest
!
              if (j1 >= j2) then
                  k = k+1
                  matrob(k) = matrix2(j1,j2)
              end if
!
 7000     end do
!
 8000 end do
!
      return
!
      end subroutine robfit
!
!***********************************************************************
!
      subroutine clock(string,clock1,index)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) string
      integer clock1,index
!-----------------------------------------------------------------------
!     Function : Calculates CPU time
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      character date*8,time*10,text*7
      integer clock2,hrs,mins,secs,length,timed,rate
!-----------------------------------------------------------------------
      length = len(string)
!
      if (index == 1) then
!
!---------*****************
          call system_clock(clock1,rate)
!---------*****************
!
          text = 'Begin: '
      else
!
!---------*****************
          call system_clock(clock2,rate)
!---------*****************
!
          text = 'End:   '
      end if
!
!-----******************
      call date_and_time(date,time)
!-----******************
!
      if (index == 1) then
          write (outbuf,'(a)') string(1:length)
          call wrtflq(outbuf)
      end if
!
      write (outbuf,'(5(a,a2),a,a4)') text, &
      time(1:2),':',time(3:4),':',time(5:6),' on ', &
      date(7:8),'/',date(5:6),'/',date(1:4)
      call wrtflq(outbuf)
!
      if (index == 2) then
          timed = nint((clock2-clock1)/real(rate))
          hrs = timed/3600
          mins = (timed - hrs*3600)/60
          secs = timed - hrs*3600 - mins*60
          write (outbuf,'(a,i3,a,2(i2,a))') 'Time = ',hrs,'hr ',mins, &
          'min ',secs,'sec'
          call wrtflq(outbuf)
          call newliq
      end if
!
      return
!
      end subroutine clock
!
!***********************************************************************
!
      subroutine multip(man1,exp1,man2,exp2,man3,exp3)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      double precision man1,man2,man3
      integer exp1,exp2,exp3
!-----------------------------------------------------------------------
!     Function : Forms the product of 2 mantissa-exponent numbers.
!-----------------------------------------------------------------------
!     MAN1 : mantissa of the first number in the product
!     EXP1 : exponent of the first number in the product
!     MAN2 : mantissa of the second number in the product
!     EXP2 : exponent of the second number in the product
!     MAN3 : mantissa of the product
!     EXP3 : exponent of the product
!-----------------------------------------------------------------------
!---- mantissa of z = mantissa of x . mantissa of y
      man3 = man1*man2
!---- exponent of z = exponent of x + exponent of y
      exp3 = exp1+exp2
!
!---- standard form of z = x.y
!-----***********
      call manexp(man3,exp3)
!-----***********
!
      return
!
      end subroutine multip
!
!***********************************************************************
!
      subroutine getlin(line,contin,iret,error,eofflag)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=SABRE_CMDLINE_MAX) :: previous_line
      logical contin,iret,error,eofflag
!-----------------------------------------------------------------------
!     Function : Gets the next line from the current input buffer.
!-----------------------------------------------------------------------
!     CONTIN : true if a continuation line
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      logical eof
!-----------------------------------------------------------------------
      error = .false.
      eofflag = .false.
!
!---- save the current line in case it gets overwritten by a
!---- continuation
!
      previous_line = line
!
!---- command line already encountered
      if (iret) then
          go to 15
      end if
!
!---- call prompt if not reading from some other INPUT channel
   10 if (inch == cinch .and. .not. contin) then
!
!-------- print SABRE prompt to terminal (subroutine PROMPT is in
!-------- subsidiary file macdep.f)
!---------***********
          call prompt('<S> ')
!---------***********
!
      else if (inch == cinch) then
!
!-------- print continuation prompt to terminal
!---------***********
          call prompt('<&> ')
!---------***********
!
      end if
!
!---- read input line
!-----***********
      call readln(cinch,line,eof)
!-----***********
!
      if (eof) then
          eofflag = .true.
!          go to 900
          go to 999
      end if
!
!---- print commands from input file to screen
   15 if (cinch /= inch) then
          call wrtscr('<S> '//line(1:SABRE_CMDLINE_MAX-6))
      end if
!
      if (line(1:3) /= 'inp' .and. line(1:3) /= 'out') then
          call wrtfil('<S> '//line(1:SABRE_CMDLINE_MAX-6))
      end if
!
      if (line(SABRE_CMDLINE_MAX-1:SABRE_CMDLINE_MAX-1) /= ' ' .or. &
         line(SABRE_CMDLINE_MAX:SABRE_CMDLINE_MAX) /= ' ') then
          call wrtlin('    *** ERROR *** COMMAND LINE TOO LONG' )
          call wrtlin('    Use the continuation symbol [&] at the end'// &
                      ' of the line to be continued.' )
          error = .true.
          go to 999
      end if
!
      if (line(1:3) /= 'out' .and. line(1:3) /= 'OUT' .and. &
          line(1:5) /= 'trace' .and. line(1:5) /= 'TRACE' .and. &
          line(1:3) /= 'res' .and. line(1:3) /= 'RES' .and. &
          line(1:3) /= 'rea' .and. line(1:3) /= 'REA') then
!
!-------- avoid the continuation of a command with a filename as
!         argument
          if (.not. (contin .and. (previous_line(1:3) == 'out' .or. &
              previous_line(1:5) == 'trace' .or. &
              previous_line(1:3) == 'res' .or. &
              previous_line(1:3) == 'rea'))) then
!
!------------ convert input line to lower case
!------------ (subroutine LOWER is in macdep.f)
!-------------**********
              call lower(line)
!-------------**********
!
          end if
!
      end if
!
      if (.not. contin) then
!
!-------- avoid changing a filename beginning OUT, TRACE or REA if it's
!-------- on a continuation line
          if (line(1:3) == 'OUT') then
              line(1:3) = 'out'
          end if
!
          if (line(1:5) == 'TRACE') then
              line(1:5) = 'trace'
          end if
!
          if (line(1:3) == 'RES') then
              line(1:3) = 'res'
          end if
!
          if (line(1:3) == 'REA') then
              line(1:3) = 'rea'
          end if
      end if
!
      go to 999
!
!---- This is all to do with reading from INPUT channels.
!---- End of file detected on input.
  900 if (contin) then
          call wrtlin('    *** ERROR *** '// &
                      'END OF FILE WHEN LOOKING FOR CONTINUATION LINE')
          contin = .false.
      end if
!
!-----***********
      call clsfil(cinch)
!-----***********
!
      cinch = inch
!
      go to 10
!
  999 return
!
      end subroutine getlin
!
!***********************************************************************
!
      subroutine next(istpos,ienpos,ibeg,line,arg)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer istpos,ienpos,ibeg
      logical arg
!-----------------------------------------------------------------------
!     Function : Returns position of the next argument on command line.
!-----------------------------------------------------------------------
!     ISTPOS : position of the first character of the argument
!     IENPOS :    "     "   "  last      "     "   "     "
!     IBEG   : position for start of search
!     ARG    : true if argument found; false if not
!-----------------------------------------------------------------------
      integer i
      logical contin,error,eofflag
!-----------------------------------------------------------------------
      arg = .false.
      i = ibeg
      istpos = 0
      ienpos = 0
!
!---- If IBEG not first column, skip over all non-blanks until blank
!---- found
!---- reading line from first character
      if (ibeg == 1) then
          go to 20
      end if
!
!---- blank found
   10 if (line(i:i) == ' ') then
          go to 20
      end if
!
!---- no blanks found on line
      if (i == SABRE_CMDLINE_MAX) then
          return
      end if
!
!---- next character
      i = i+1
!
      go to 10
!
!---- Now skip over blanks until next item
!---- non-blank character found
   20 if (line(i:i) /= ' ') then
          go to 25
      end if
!
!---- no non-blank characters on line
      if (i == SABRE_CMDLINE_MAX) then
          return
      end if
!
!---- next character
      i = i+1
!
      go to 20
!
!---- Start of argument found; check for continuation character &
   25 if (line(i:i) /= '&') then
          go to 30
      end if
!
!---- continuation symbol found
      contin = .true.
!
!---- get continuation input line
!-----***********
      call getlin(line,contin,.false.,error,eofflag)
!-----***********
!
      if (.not. contin) then
          return
      end if
!
!---- reading continuation line from first character
      i = 1
!
      go to 20
!
!---- start position
   30 istpos = i
      i = i+1
!
!---- end of argument found or end of line
   32 if (line(i:i) == ' ' .or. i == SABRE_CMDLINE_MAX) then
          go to 40
      end if
!
!---- next character
      i = i+1
!
      go to 32
!
!---- end position
   40 ienpos = i-1
      arg = .true.
!
      return
!
      end subroutine next
!
!***********************************************************************
!
      subroutine charel(value,mchars,nch,error)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character mchars*(*)
      double precision value
      integer nch
      logical error
!-----------------------------------------------------------------------
!     Function : Converts a character string to a floating point number.
!                Format is (blanks) (+/-) (digits) (.) (digits)
!                (D/E and/or +/- digits). Note that embedded and
!                trailing blanks are treated as zeroes.
!-----------------------------------------------------------------------
!     VALUE     : real value returned
!     MCHARS(.) : character string input
!     NCH       : number of characters in the string
!-----------------------------------------------------------------------
      character mchar,mdigit*10
      double precision accum
      integer nchars,ichar,lsign,iscale,index,lsgne,iexpt
!-----------------------------------------------------------------------
      data mdigit /'0123456789'/
!-----------------------------------------------------------------------
!---- number of characters in string
      nchars = nch
!
!---- non-blank character found
  100 if (mchars(nchars:nchars) /= ' ') then
          go to 101
      end if
!
!---- work backwards through string
      nchars = nchars-1
!
!---- stop reading string immediately before initial character
      if (nchars >= 2) then
          go to 100
      end if
!
!---- set up pointer to the first non-blank: if the entire field is
!---- blank, return value = 0.0
  101 value = 0
      error = .false.
      ichar = 0
!---- update character index
  220 ichar = ichar+1
!
      if (ichar > nchars) then
          return
      end if
!
!---- ignore leading blanks
      if (mchars(ichar:ichar) == ' ') then
          go to 220
      end if
!
!---- note the sign, if any
      lsign = 0
!
!---- positive sign detected
      if (mchars(ichar:ichar) == '+') then
          go to 340
      end if
!
!---- negative sign detected
      if (mchars(ichar:ichar) /= '-') then
          go to 360
      end if
!
      lsign = -1
!---- update character index
  340 ichar = ichar+1
!---- initialise accumulator and power
  360 accum = 0
      iscale = 0
!
!---- process the integer part, if any
  420 if (ichar > nchars) then
          go to 800
      end if
!
      mchar = mchars(ichar:ichar)
!
!---- blank character string
      if (mchar == ' ') then
          mchar = mdigit(1:1)
      end if
!
      index = 0
!---- update integer value corresponding to character
  440 index = index+1
!
      if (index > 10) then
          go to 460
      end if
!
      if (mchar /= mdigit(index:index)) then
          go to 440
      end if
!
!---- incorporate next digit into integer part
      accum = 10*accum + index - 1
      ichar = ichar+1
!
      go to 420
!
!---- process the fractional part, if any
  460 if (mchars(ichar:ichar) /= '.') then
          go to 580
      end if
!
      ichar = ichar+1
!
  520 if (ichar > nchars) then
          go to 800
      end if
!
      mchar = mchars(ichar:ichar)
!
      if (mchar == ' ') then
          mchar = mdigit(1:1)
      end if
!
      index = 0
  540 index = index+1
!
      if (index > 10) then
          go to 580
      end if
!
      if (mchar /= mdigit(index:index)) then
          go to 540
      end if
!
!---- incorporate next digit into fractional part
      accum = 10*accum + index - 1
      iscale = iscale-1
      ichar = ichar+1
!
      go to 520
!
!---- process the exponent - if none, error, shouldn't have got here
!---- could be any of E,D,+,-,E+,E-,D+,D-, so skip to the first digit
!---- and note the sign (note: blank(s) allowed between D/E and +/-)
  580 mchar = mchars(ichar:ichar)
!
      if (mchar == 'd') then
          go to 620
      end if
!
      if (mchar /= 'e') then
          go to 640
      end if
!
  620 ichar = ichar+1
!
      if (ichar > nchars) then
          go to 800
      end if
!
      mchar = mchars(ichar:ichar)
!
      if (mchar == ' ') then
          go to 620
      end if
!
  640 lsgne = 0
!
!---- positive exponent sign detected
      if (mchar == '+') then
          go to 660
      end if
!
!---- negative exponent sign detected
      if (mchar /= '-') then
          go to 680
      end if
!
      lsgne = -1
  660 ichar = ichar+1
!---- process the digits and add to 'ISCALE'
  680 iexpt = 0
!
  720 if (ichar > nchars) then
          go to 760
      end if
!
      mchar = mchars(ichar:ichar)
!
      if (mchar == ' ') then
          mchar = mdigit(1:1)
      end if
!
      index = 0
  740 index = index+1
!
      if (index > 10) then
          error = .true.
          return
      end if
!
      if (mchar /= mdigit(index:index)) then
          go to 740
      end if
!
!---- incorporate next digit into exponent
      iexpt = 10*iexpt + index-1
      ichar = ichar+1
!
      go to 720
!
!---- negative exponent
  760 if (lsgne < 0) then
          iexpt = -iexpt
      end if
!
      iscale = iscale+iexpt
!
!---- end of field
!---- combine SIGN, ACCUM, and SCALE
!---- generate real number
  800 if (iscale /= 0) then
          accum = accum*10d0**iscale
      end if
!
!---- negative real number
      if (lsign < 0) then
          accum = -accum
      end if
!
!---- return rounded double precision result - rounding must be done
!---- after signing to allow for both signed-fraction and two's-
!---- complement-fraction machines
!---- double precision
      value = accum
      value = accum + 2*(accum-value)
!
      return
!
      end subroutine charel
!
!***********************************************************************
!
      subroutine chaint(mnum,inum,nchars,error)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character mnum*(*)
      integer inum,nchars
      logical error
!-----------------------------------------------------------------------
!     Function : Converts a character string into an integer.
!-----------------------------------------------------------------------
!     MNUM(.) : character string input
!     INUM    : integer value returned
!     NCHARS  : number of characters in the string
!-----------------------------------------------------------------------
      character mdgts*10,mchar
      integer ichar,kdgt
      logical negtiv
!-----------------------------------------------------------------------
      data mdgts /'0123456789'/
!-----------------------------------------------------------------------
      error = .false.
!---- initialise the output value to zero, and wind past any leading
!---- blanks, jumping to the exit if there are no non-blanks here
      inum = 0
      ichar = 0
  220 ichar = ichar+1
!
      if (ichar > nchars) then
          return
      end if
!
      mchar = mnum(ichar:ichar)
!
      if (mchar == ' ') then
          go to 220
      end if
!
!---- note the sign, if there is one, and advance the pointer past it
      negtiv = .false.
!
      if (ichar == nchars) then
          go to 400
      end if
!
      if (mchar /= '-') then
          go to 340
      end if
!
!---- negative sign detected
      negtiv = .true.
      ichar = ichar+1
!
      go to 400
!
!---- positive sign detected
  340 if (mchar == '+') then
          ichar = ichar+1
      end if
!
!---- convert the digit string or take the error exit if non-digit found
  400 mchar = mnum(ichar:ichar)
!
      if (mchar == ' ') then
          mchar = mdgts(1:1)
      end if
!
      kdgt = 0
  440 kdgt = kdgt+1
!
      if (kdgt > 10) then
          error = .true.
          return
      end if
!
      if (mchar /= mdgts(kdgt:kdgt)) then
          go to 440
      end if
!
!---- incorporate next digit into number
      inum = 10*inum + kdgt - 1
      ichar = ichar+1
!
      if (ichar <= nchars) then
          go to 400
      end if
!
!---- sign, and sign off; negative integer
      if (negtiv) then
          inum = -inum
      end if
!
      return
!
      end subroutine chaint
!
!***********************************************************************
!
      subroutine argrel(line,rarg)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      double precision rarg
!-----------------------------------------------------------------------
!     Function : Reads a single non-negative real argument.
!-----------------------------------------------------------------------
!     RARG : real valued argument
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      double precision rtemp
      integer istpos,ienpos
      logical arg,error
!-----------------------------------------------------------------------
!---- RARG is the integer argument. The previous value of RARG is stored
!---- in RTEMP so that if the argument is not correct, RARG can be
!---- reassigned its previous value and not the mistaken value.
      rtemp = rarg
!
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** MISSING ARGUMENT')
          return
      end if
!
!---- convert character to real
!-----***********
      call charel(rarg,line(istpos:ienpos),ienpos-istpos+1,error)
!-----***********
!
!---- ERROR is returned false if the argument has been correctly read
!---- as a real
      if (error .or. rarg < 0) then
          call wrtlin('    *** ERROR *** '// &
          'ARGUMENT IS NOT A NON-NEGATIVE REAL')
          rarg = rtemp
      end if
!
      return
!
      end subroutine argrel
!
!***********************************************************************
!
      subroutine argint(line,iarg,itemp)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer iarg,itemp
!-----------------------------------------------------------------------
!     Function : Reads a single non-negative integer argument.
!-----------------------------------------------------------------------
!     IARG  : integer valued argument
!     ITEMP : previous value of IARG
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer istpos,ienpos
      logical arg,error
!-----------------------------------------------------------------------
!---- IARG is the integer argument. The previous value of IARG is stored
!---- in ITEMP so that if the argument is not correct, IARG can be
!---- reassigned its previous value and not the mistaken value.
      itemp = iarg
!
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** MISSING ARGUMENT')
          return
      end if
!
!---- convert character into integer
!-----***********
      call chaint(line(istpos:ienpos),iarg,ienpos-istpos+1,error)
!-----***********
!
!---- ERROR is returned false if the argument has been correctly read
!---- as an integer
      if (error .or. iarg < 0) then
          call wrtlin('    *** ERROR *** '// &
          'ARGUMENT IS NOT A NON-NEGATIVE INTEGER')
          iarg = itemp
      end if
!
      return
!
      end subroutine argint
!
!***********************************************************************
!
      subroutine argcar(line,larg)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      logical larg
!-----------------------------------------------------------------------
!     Function : Reads a y/n argument and sets logical to true or false.
!-----------------------------------------------------------------------
!     LARG : true if argument is y, false if n
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character c
      integer istpos,ienpos
      logical arg
!-----------------------------------------------------------------------
!---- get position of argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      if (.not. arg) then
          call wrtlin('    *** ERROR *** no argument')
          return
      else
!-------- argument = y/n
          c = line(istpos:istpos)
      end if
!
      larg = .false.
!
      if (c == 'y') then
          larg = .true.
      else if (c /= 'n') then
          call wrtlin('    *** ERROR *** INCORRECT ARGUMENT')
      end if
!
      return
!
      end subroutine argcar
!
!***********************************************************************
!
      subroutine filini
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'fildat.h'
!-----------------------------------------------------------------------
!     Initialise line buffer for log file output
!-----------------------------------------------------------------------
      chrpos = 0
!
      return
!
      end subroutine filini
!
!***********************************************************************
!
      subroutine wrtlin(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the screen and optionally to the log file
!     as well, moving the cursor to the start of the next line.
!     If terse output required then wrtlin does not output to the screen
!     Any screen output will come from calling wterse
!-----------------------------------------------------------------------
!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
!-----------------------------------------------------------------------
      if (.not. terse) then
          call wrtscr(text(1:length(text)))
      end if

      call wrtfil(text(1:length(text)))
!
      return
!
      end subroutine wrtlin
!
!***********************************************************************
!
      subroutine wterse(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the screen if terse output is required
!-----------------------------------------------------------------------
!     John Pritchard 03-05-07
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
!-----------------------------------------------------------------------
      if (terse) then
          call wrtscr(text(1:length(text)))
      end if
!
      return
!
      end subroutine wterse
!
!***********************************************************************
!
      subroutine wrtln2(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the screen and optionally to the log file
!     as well, leaving the cursor at the end of the line.
!-----------------------------------------------------------------------
!     John Pritchard, 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
!-----***********
      call prompt(text)
!-----***********
!
      call wrtfl2(text)
!
      return
!
      end subroutine wrtln2
!
!***********************************************************************
!
      subroutine wrtlit(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the trace file,
!     moving the cursor to the start of the next line.
!-----------------------------------------------------------------------
!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
!-----------------------------------------------------------------------
      call wrtfit(text(1:length(text)))
!
      return
!
      end subroutine wrtlit
!
!***********************************************************************
!
      subroutine wrtlir(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line of text to the resid file,
!     moving the cursor to the start of the next line.
!-----------------------------------------------------------------------
!     John Pritchard 10-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer length
!-----------------------------------------------------------------------
      call wrtfir(text(1:length(text)))
!
      return
!
      end subroutine wrtlir
!
!***********************************************************************
!
      subroutine wrtfil(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line to the log file, prepending anything in the log file
!     line buffer first
!-----------------------------------------------------------------------
      include 'fildat.h'
!-----------------------------------------------------------------------
      character(len=80) :: temp
      integer length
!-----------------------------------------------------------------------
      if (chrpos > 0) then
          temp = linbuf(1:chrpos)//text
          call wrtfln(temp(1:length(temp)))
          chrpos = 0
      else
          call wrtfln(text(1:length(text)))
      end if
!
      return
!
      end subroutine wrtfil
!
!***********************************************************************
!
      subroutine wrtfl2(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Add text to line buffer for log file output
!-----------------------------------------------------------------------
      include 'fildat.h'
!-----------------------------------------------------------------------
      character(len=80) :: temp
!-----------------------------------------------------------------------
      if (chrpos == 0) then
          linbuf = text
          chrpos = len(text)
      else
          temp = linbuf
          linbuf = temp(1:chrpos)//text
          chrpos = chrpos + len(text)
      end if
!
      return
!
      end subroutine wrtfl2
!
!***********************************************************************
!
      subroutine wrtfit(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line to the trace file, prepending anything in the log
!     file line buffer first
!-----------------------------------------------------------------------
      include 'fildat.h'
!-----------------------------------------------------------------------
      character(len=80) :: temp
      integer length
!-----------------------------------------------------------------------
      if (chrpos > 0) then
          temp = linbuf(1:chrpos)//text
          call wrtflt(temp(1:length(temp)))
          chrpos = 0
      else
          call wrtflt(text(1:length(text)))
      end if
!
      return
!
      end subroutine wrtfit
!
!***********************************************************************
!
      subroutine wrtfir(text)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) text
!-----------------------------------------------------------------------
!     Write a line to the resid file, prepending anything in the log
!     file line buffer first
!-----------------------------------------------------------------------
      include 'fildat.h'
!-----------------------------------------------------------------------
      character(len=80) :: temp
      integer length
!-----------------------------------------------------------------------
      if (chrpos > 0) then
          temp = linbuf(1:chrpos)//text
          call wrtflr(temp(1:length(temp)))
          chrpos = 0
      else
          call wrtflr(text(1:length(text)))
      end if
!
      return
!
      end subroutine wrtfir
!
!***********************************************************************
!
      subroutine newlin
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the screen and optionally to the log file
!     as well.
!
!     Does nothing to screen if terse output is required
!-----------------------------------------------------------------------
!     John Pritchard 17-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      include 'fildat.h'
!-----------------------------------------------------------------------
      if (.not. terse) then
          call newlns
      end if
!
!---- Write anything in the file line buffer first
      if (chrpos > 0) then
          call wrtfln(linbuf(1:chrpos))
          chrpos = 0
      else
          call newlnf
      end if
!
      return
!
      end subroutine newlin
!
!***********************************************************************
!
      subroutine newlit
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the trace file
!-----------------------------------------------------------------------
!     John Pritchard 17-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      include 'fildat.h'
!-----------------------------------------------------------------------
!---- Write anything in the file line buffer first
      if (chrpos > 0) then
          call wrtflt(linbuf(1:chrpos))
          chrpos = 0
      else
          call newlnt
      end if
!
      return
!
      end subroutine newlit
!
!***********************************************************************
!
      subroutine newlir
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the resid file
!-----------------------------------------------------------------------
!     John Pritchard 17-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      include 'fildat.h'
!-----------------------------------------------------------------------
!---- Write anything in the file line buffer first
      if (chrpos > 0) then
          call wrtflr(linbuf(1:chrpos))
          chrpos = 0
      else
          call newlnr
      end if
!
      return
!
      end subroutine newlir
!
!***********************************************************************
!
      subroutine newliq
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!     Write a new line to the time file
!-----------------------------------------------------------------------
!     John Pritchard 17-03-04
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      include 'fildat.h'
!-----------------------------------------------------------------------
!---- Write anything in the file line buffer first
      if (chrpos > 0) then
          call wrtflq(linbuf(1:chrpos))
          chrpos = 0
      else
          call newlnq
      end if
!
      return
!
      end subroutine newliq
!
!***********************************************************************
!
      subroutine multi_int(line,intarg,prefix,maxargs,nargs,arg_index)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs,intarg(maxargs),nargs,arg_index(maxargs)
      character*(*) prefix(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads one or more integer arguments prefixed
!     by keywords in the character array prefix.  In the case of a
!     single integer argument without a prefix, prefix(1) is assumed.
!-----------------------------------------------------------------------
!     INTARG : integer valued argument(s)
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer istpos,ienpos,argpos,start,equals_pos,iarg,itemp
      logical arg,error
!-----------------------------------------------------------------------
      start = 4
      nargs = 0
!
   10 continue
!
!---- get position of next argument
!-----*********
      call next(istpos,ienpos,start,line,arg)
!-----*********
!
      if (.not. arg) then
          go to 100
      end if
!
!---- look for prefix
      argpos = 0
!
      if (ienpos >= istpos+2) then
!
          do 20 iarg = 1,maxargs
!
              if (line(istpos:istpos+2) == prefix(iarg)(1:3)) then
                  argpos = iarg
              end if
!
   20     end do
!
      end if
!
      if (argpos == 0 .and. nargs == 0) then
!-------- No prefix found and no argument read already
!-------- in which case a single integer value is valid
          argpos = 1
!
!---------***********
          call chaint(line(istpos:ienpos),itemp,ienpos-istpos+1,error)
!---------***********
!
          if (error) then
              call wrtlin('    *** ERROR *** '// &
              'ARGUMENT IS NOT AN INTEGER')
              return
          else
!------------ argument is OK, return immediately since only one argument
!------------ is allowed without a prefix
              intarg(argpos) = itemp
              nargs = 1
              arg_index(nargs) = argpos
              return
          end if
!
      else if (argpos == 0 .and. nargs > 0) then
!-------- One or more arguments already read and next argument is not
!-------- a valid prefix
          call wrtlin('    *** ERROR *** Invalid argument')
          return
      end if
!
!---- Get to here if a valid prefix has been read
!---- See if "=" is present in the argument and if so make that the new
!---- starting point.  If "=" is not present then get the next argument
      equals_pos = index(line(istpos:ienpos),'=')
!
      if (equals_pos > 0) then
          istpos = istpos+equals_pos-1
      else
          start = ienpos
!
!---------*********
          call next(istpos,ienpos,start,line,arg)
!---------*********
!
          if (.not. arg) then
              call wrtlin('    *** ERROR *** MISSING ARGUMENT')
              return
          end if
!
      end if
!
!---- Next argument can be "=" by itself, an integer by itself or
!---- an integer prefixed by "="
      if (ienpos == istpos .and. line(istpos:ienpos) == '=') then
!-------- "=" by itself found, get the next argument
          start = ienpos
!
!---------*********
          call next(istpos,ienpos,start,line,arg)
!---------*********
!
          if (.not. arg) then
              call wrtlin('    *** ERROR *** MISSING ARGUMENT')
              return
          else
!
!------------ an argument has been found, now see if it is a valid
!             integer
!-------------***********
              call chaint(line(istpos:ienpos),itemp,ienpos-istpos+1, &
                          error)
!-------------***********
!
              if (error) then
                  call wrtlin('    *** ERROR *** '// &
                  'ARGUMENT IS NOT AN INTEGER')
                  return
              else
                  intarg(argpos) = itemp
              end if
!
!------------ Get to here if a valid integer has been found after "="
              nargs = nargs+1
              arg_index(nargs) = argpos
              start = ienpos
          end if
!
      else if (ienpos > istpos .and. line(istpos:istpos) == '=') &
      then
!
!-------- argument prefixed by "=" found, see if the rest of the
!-------- argument is a valid integer
!---------***********
          call chaint(line(istpos+1:ienpos),itemp,ienpos-istpos,error)
!---------***********
!
          if (error) then
              call wrtlin('    *** ERROR *** '// &
              'ARGUMENT IS NOT AN INTEGER')
              return
          else
              intarg(argpos) = itemp
          end if
!
!-------- Get to here if a valid integer has been found after "=" prefix
          nargs = nargs+1
          arg_index(nargs) = argpos
          start = ienpos
      else
!
!-------- No "=" present, see if the argument is a valid integer
!---------***********
          call chaint(line(istpos:ienpos),itemp,ienpos-istpos+1,error)
!---------***********
!
          if (error) then
              call wrtlin('    *** ERROR *** '// &
              'ARGUMENT IS NOT AN INTEGER')
              return
          else
              intarg(argpos) = itemp
          end if
!
!-------- Get to here if a valid integer has been found after the
!-------- argument prefix
          nargs = nargs+1
          arg_index(nargs) = argpos
          start = ienpos
      end if
!
!---- Get to here if a valid integer argument with prefix has been found
!---- If maximum number of arguments have not been read see if there
!---- is another argument
      if (nargs < maxargs) then
          go to 10
      end if
!
  100 continue
!
      if (nargs == 0) then
          call wrtlin('    *** ERROR *** MISSING ARGUMENT')
      end if
!
      return
!
      end subroutine multi_int
!
!***********************************************************************
!
      subroutine multi_real(line,rarg,prefix,maxargs,nargs,arg_index)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs,nargs,arg_index(maxargs)
      character*(*) prefix(maxargs)
      double precision rarg(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads one or more real arguments prefixed
!     by keywords in the character array prefix.  In the case of a
!     single real argument without a prefix, prefix(1) is assumed.
!-----------------------------------------------------------------------
!     RARG : real valued argument(s)
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      double precision rtemp
      integer istpos,ienpos,argpos,start,equals_pos,iarg
      logical arg,error
!-----------------------------------------------------------------------
      start = 4
      nargs = 0
!
   10 continue
!
!---- get position of next argument
!-----*********
      call next(istpos,ienpos,start,line,arg)
!-----*********
!
      if (.not. arg) then
          go to 100
      end if
!
!---- look for prefix
      argpos = 0
!
      if (ienpos >= istpos+2) then
!
          do 20 iarg = 1,maxargs
!
              if (line(istpos:istpos+2) == prefix(iarg)(1:3)) then
                  argpos = iarg
              end if
!
   20     end do
!
      end if
!
      if (argpos == 0 .and. nargs == 0) then
!-------- No prefix found and no argument read already
!-------- in which case a single real value is valid
          argpos = 1
!
!-------- Read into temporary variable first for checking
!---------***********
          call charel(rtemp,line(istpos:ienpos),ienpos-istpos+1,error)
!---------***********
!
          if (error) then
              call wrtlin('    *** ERROR *** '// &
              'ARGUMENT IS NOT A REAL')
              return
          else
!------------ argument is OK, return immediately since only one argument
!------------ is allowed without a prefix
              rarg(argpos) = rtemp
              nargs = 1
              arg_index(nargs) = argpos
              return
          end if
!
      else if (argpos == 0 .and. nargs > 0) then
!-------- One or more arguments already read and next argument is not
!-------- a valid prefix
          call wrtlin('    *** ERROR *** Invalid argument')
          return
      end if
!
!---- Get to here if a valid prefix has been read
!---- See if "=" is present in the argument and if so make that the new
!---- starting point.  If "=" is not present then get the next argument
      equals_pos = index(line(istpos:ienpos),'=')
!
      if (equals_pos > 0) then
          istpos = istpos+equals_pos-1
      else
          start = ienpos
!
!---------*********
          call next(istpos,ienpos,start,line,arg)
!---------*********
!
          if (.not. arg) then
              call wrtlin('    *** ERROR *** MISSING ARGUMENT')
              return
          end if
!
      end if
!
!---- Next argument can be "=" by itself, a real number by itself or
!---- a real number prefixed by "="
      if (ienpos == istpos .and. line(istpos:ienpos) == '=') then
!-------- "=" by itself found, get the next argument
          start = ienpos
!
!---------*********
          call next(istpos,ienpos,start,line,arg)
!---------*********
!
          if (.not. arg) then
              call wrtlin('    *** ERROR *** MISSING ARGUMENT')
              return
          else
!
!------------ an argument has been found, now see if it is a valid real
!-------------***********
              call charel(rtemp,line(istpos:ienpos),ienpos-istpos+1, &
                          error)
!-------------***********
!
              if (error) then
                  call wrtlin('    *** ERROR *** '// &
                  'ARGUMENT IS NOT A REAL')
                  return
              else
                  rarg(argpos) = rtemp
              end if
!
!------------ Get to here if a valid real has been found after "="
              nargs = nargs+1
              arg_index(nargs) = argpos
              start = ienpos
          end if
!
      else if (ienpos > istpos .and. line(istpos:istpos) == '=') &
      then
!
!-------- argument prefixed by "=" found, see if the rest of the
!-------- argument is a valid real
!---------***********
          call charel(rtemp,line(istpos+1:ienpos),ienpos-istpos,error)
!---------***********
!
          if (error) then
              call wrtlin('    *** ERROR *** '// &
              'ARGUMENT IS NOT A REAL')
              return
          else
              rarg(argpos) = rtemp
          end if
!
!-------- Get to here if a valid real has been found after "=" prefix
          nargs = nargs+1
          arg_index(nargs) = argpos
          start = ienpos
      else
!
!-------- No "=" present, see if the argument is a valid real
!---------***********
          call charel(rtemp,line(istpos:ienpos),ienpos-istpos+1,error)
!---------***********
!
          if (error) then
              call wrtlin('    *** ERROR *** '// &
              'ARGUMENT IS NOT A REAL')
              return
          else
              rarg(argpos) = rtemp
          end if
!
!-------- Get to here if a valid real has been found after the argument
!         prefix
          nargs = nargs+1
          arg_index(nargs) = argpos
          start = ienpos
      end if
!
!---- Get to here if a valid real argument with prefix has been found
!---- If maximum number of arguments have not been read see if there
!---- is another argument
      if (nargs < maxargs) then
          go to 10
      end if
!
  100 continue
!
      return
!
      end subroutine multi_real
!
!***********************************************************************
!
      subroutine multi_char(line,charg,prefix,maxargs,nargs,arg_index)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs,nargs,arg_index(maxargs)
      character*(*) charg(maxargs),prefix(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads one or more string arguments prefixed
!     by keywords in the character array prefix.  In the case of a
!     single string argument without a prefix, prefix(1) is assumed.
!-----------------------------------------------------------------------
!     CHARG : character valued argument(s)
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer istpos,ienpos,argpos,start,equals_pos,iarg
      logical arg
!-----------------------------------------------------------------------
      start = 4
      nargs = 0
!
   10 continue
!
!---- get position of next argument
!-----*********
      call next(istpos,ienpos,start,line,arg)
!-----*********
!
      if (.not. arg) then
          go to 100
      end if
!
!---- look for prefix
      argpos = 0
!
      if (ienpos >= istpos+2) then
!
          do 20 iarg = 1,maxargs
!
              if (line(istpos:istpos+2) == prefix(iarg)(1:3)) then
                  argpos = iarg
              end if
!
   20     end do
!
      end if
!
      if (argpos == 0 .and. nargs == 0) then
!-------- No prefix found and no argument read already
!-------- in which case a single string value is valid
          argpos = 1
          charg(argpos) = line(istpos:ienpos)
          nargs = 1
          arg_index(nargs) = argpos
!-------- return immediately since only one argument
!-------- is allowed without a prefix
          return
      else if (argpos == 0 .and. nargs > 0) then
!-------- One or more arguments already read and next argument is not
!-------- a valid prefix
          call wrtlin('    *** ERROR *** Invalid argument')
          return
      end if
!
!---- Get to here if a valid prefix has been read
!---- See if "=" is present in the argument and if so make that the new
!---- starting point.  If "=" is not present then get the next argument
      equals_pos = index(line(istpos:ienpos),'=')
!
      if (equals_pos > 0) then
          istpos = istpos+equals_pos-1
      else
          start = ienpos
!
!---------*********
          call next(istpos,ienpos,start,line,arg)
!---------*********
!
          if (.not. arg) then
              call wrtlin('    *** ERROR *** MISSING ARGUMENT')
              return
          end if
!
      end if
!
!---- Next argument can be "=" by itself, a string value by itself or
!---- an integer prefixed by "="
      if (ienpos == istpos .and. line(istpos:ienpos) == '=') then
!-------- "=" by itself found, get the next argument
          start = ienpos
!
!---------*********
          call next(istpos,ienpos,start,line,arg)
!---------*********
!
          if (.not. arg) then
              call wrtlin('    *** ERROR *** MISSING ARGUMENT')
              return
          else
!------------ an argument has been found
              charg(argpos) = line(istpos:ienpos)
              nargs = nargs+1
              arg_index(nargs) = argpos
              start = ienpos
          end if
!
      else if (ienpos > istpos .and. line(istpos:istpos) == '=') &
      then
!-------- argument prefixed by "=" found, get the rest of the argument
          charg(argpos) = line(istpos+1:ienpos)
          nargs = nargs+1
          arg_index(nargs) = argpos
          start = ienpos
      else
!-------- No "=" present
          charg(argpos) = line(istpos:ienpos)
          nargs = nargs+1
          arg_index(nargs) = argpos
          start = ienpos
      end if
!
!---- Get to here if a string argument with prefix has been found
!---- If maximum number of arguments have not been read see if there
!---- is another argument
      if (nargs < maxargs) then
          go to 10
      end if
!
  100 continue
!
      if (nargs == 0) then
          call wrtlin('    *** ERROR *** MISSING ARGUMENT')
      end if
!
      return
!
      end subroutine multi_char
!
!***********************************************************************
!
      subroutine multi_var(line,vname,prefix,maxargs,name,nvar,iread, &
                           nargs)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: name(maxvar)
      integer maxargs,nvar,nargs
      character*(*) vname(maxargs),prefix(maxargs)
      logical iread
!-----------------------------------------------------------------------
!     Function : Reads one or more string arguments prefixed
!     by keywords in the character array prefix.  In the case of a
!     single string argument without a prefix, prefix(1) is assumed.
!     These are then checked to see that they are valid variable names.
!-----------------------------------------------------------------------
!     vname : character valued argument(s)
!-----------------------------------------------------------------------
      integer iarg,ivar,length,arg_index(3)
!-----------------------------------------------------------------------
!---- Check to see that variable names have been defined already
      if (iread) then
!
!-------- Names exist so read the arguments
!---------***************
          call multi_char(line,vname,prefix,maxargs,nargs,arg_index)
!---------***************
!
!-------- Check to see that arguments given are valid variable names
          do 30 iarg = 1,nargs
!
              do 10 ivar = 1,nvar
!
                  if (vname(arg_index(iarg))(1:50) == &
                  name(ivar)(1:50)) then
                      go to 20
                  end if
!
   10         end do
!
!------------ Reach here if no match is found for vname(iarg)
              call wrtlin('    *** ERROR *** '// &
                          vname(arg_index(iarg)) &
                          (1:length(vname(iarg)))// &
                          ' IS NOT A VARIABLE NAME')
              return
!
!------------ Reach here if a match is found
   20         continue
!
   30     end do
      else
!-------- No data read yet
          call wrtlin('    *** ERROR *** NO DATA HAS BEEN READ')
      end if
!
      return
!
      end subroutine multi_var
!
!***********************************************************************
!
      subroutine multi_scale(line,scale,prefix,maxargs,scaflag)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs
      character*(*) prefix(maxargs)
      double precision scale(maxargs)
      logical scaflag(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads one or more scale values prefixed
!     by keywords in the character array prefix.  In the case of a
!     single string argument without a prefix, prefix(1) is assumed.
!     These are then checked to see that they are non-negative.
!-----------------------------------------------------------------------
!     scale() : non-negative real valued argument(s)
!-----------------------------------------------------------------------
      integer nargs,iarg,arg_index(3)
!-----------------------------------------------------------------------
!---- Read the arguments
!-----***************
      call multi_real(line,scale,prefix,maxargs,nargs,arg_index)
!-----***************
!
!---- Initialise the flags
      do 1 iarg = 1,maxargs
          scaflag(iarg) = .false.
    1 end do
!
!---- Check to see that arguments given are non-negative
      do 10 iarg = 1,nargs
!
          if (scale(arg_index(iarg)) < 0) then
              call wrtlin('    *** ERROR *** '// &
                          'ARGUMENT IS NOT A NON-NEGATIVE REAL')
          else
              scaflag(arg_index(iarg)) = .true.
          end if
!
   10 end do
!
      return
!
      end subroutine multi_scale
!
!***********************************************************************
!
      subroutine multi_sigma(line,sigma,prefix,maxargs,sigflag)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs
      character*(*) prefix(maxargs)
      double precision sigma(maxargs)
      logical sigflag(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads one or more sigma values prefixed
!     by keywords in the character array prefix.  In the case of a
!     single string argument without a prefix, prefix(1) is assumed.
!     These are then checked to see that they are non-negative.
!-----------------------------------------------------------------------
!     sigma() : non-negative real valued argument(s)
!-----------------------------------------------------------------------
      integer nargs,iarg,arg_index(3)
!-----------------------------------------------------------------------
!---- Read the arguments
!-----***************
      call multi_real(line,sigma,prefix,maxargs,nargs,arg_index)
!-----***************
!
!---- Initialise the flags
      do 1 iarg = 1,maxargs
          sigflag(iarg) = .false.
    1 end do
!
!---- Check to see that arguments given are non-negative and set the
!---- flag for each valid argument given
      do 10 iarg = 1,nargs
!
          if (sigma(arg_index(iarg)) < 0) then
              call wrtlin('    *** ERROR *** '// &
                          'ARGUMENT IS NOT A NON-NEGATIVE REAL')
          else
              sigflag(arg_index(iarg)) = .true.
          end if
!
   10 end do
!
      return
!
      end subroutine multi_sigma
!
!***********************************************************************
!
      subroutine multi_rho(line,rho,prefix,maxargs,rhoflag)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      integer maxargs
      character*(*) prefix(maxargs)
      double precision rho(maxargs)
      logical rhoflag(maxargs)
!-----------------------------------------------------------------------
!     Function : Reads one or more rho values prefixed
!     by keywords in the character array prefix.  In the case of a
!     single string argument without a prefix, prefix(1) is assumed.
!     These are then checked to see that they lie between -1 and 1.
!-----------------------------------------------------------------------
!     rho() : real valued argument(s) between -1 and 1
!-----------------------------------------------------------------------
      integer nargs,iarg,arg_index(3)
!-----------------------------------------------------------------------
!---- Read the arguments
!-----***************
      call multi_real(line,rho,prefix,maxargs,nargs,arg_index)
!-----***************
!
!---- Initialise the flags
      do 1 iarg = 1,maxargs
          rhoflag(iarg) = .false.
    1 end do
!
!---- Check to see that arguments given are between -1 and 1
      do 10 iarg = 1,nargs
!
          if (rho(arg_index(iarg)) < -1d0 .or. &
          rho(arg_index(iarg)) > 1d0) then
              call wrtlin('    *** ERROR *** '// &
                          'ARGUMENT MUST LIE BETWEEN -1 and 1')
          else
              rhoflag(arg_index(iarg)) = .true.
          end if
!
   10 end do
!
      return
!
      end subroutine multi_rho
!
!***********************************************************************
!
      integer function length(string)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character*(*) string
!-----------------------------------------------------------------------
!     Returns the length of the string minus trailing spaces or nulls
!-----------------------------------------------------------------------
      integer i
!-----------------------------------------------------------------------
      length = len(string)
!
      do 10 i = length,1,-1
!
          if (string(i:i) /= ' ') then
              go to 20
          end if
!
          length = length - 1
   10 end do
!
!---- Check for null as end of string
   20 continue
!
      do 30 i = 1,length
!
          if (ichar(string(i:i)) == 0) then
              length = i-1
              return
          end if
!
   30 end do
!
      return
!
      end function length
!
!***********************************************************************
!
      subroutine manexp(mant,exp)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      double precision mant
      integer exp
!-----------------------------------------------------------------------
!     Function : Stores a mantissa-exponent number in standard form.
!-----------------------------------------------------------------------
!     MANT : the mantissa of the number
!     EXP  : the exponent of the number
!-----------------------------------------------------------------------
!---- x = 0
      if (mant == 0) then
          exp = 0
          return
      end if
!
!---- |x| < 1
      if (abs(mant) < 1) then
!-------- update reduction in exponent
   10     exp = exp-1
!-------- update mantissa
          mant = mant*10
!
!-------- absolute value of mantissa of x now lies between 1 and 10;
!-------- i.e. written in standard form x = mant.10^exp,
!-------- where 1 <= |mant| < 10
          if (abs(mant) < 1) then
              go to 10
          end if
!
!---- |x| >= 1
      else if (abs(mant) >= 10) then
!-------- update exponent
   20     exp = exp+1
!-------- update mantissa
          mant = mant/10
!
!-------- absolute value of mantissa of x now lies between 1 and 10;
!-------- i.e. written in standard form x = mant.10^exp,
!-------- where 1 <= |mant| < 10
          if (abs(mant) >= 10) then
              go to 20
          end if
!
      end if
!
      end subroutine manexp
!
!***********************************************************************
!
      subroutine normal(man1,exp1,number)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      double precision man1,number
      integer exp1
!-----------------------------------------------------------------------
!     Function : Changes a mantissa-exponent number back into real form.
!-----------------------------------------------------------------------
!     MAN1   : mantissa of the number to be converted to the normal form
!     EXP1   : exponent of the number to be converted to the normal form
!     NUMBER : real number expressed in normal form
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      double precision num
      integer i
!-----------------------------------------------------------------------
!---- let x = a.10^b
!---- b = 0
      if (exp1 == 0) then
!-------- calculate x = a.10^b = a.10^0 = a
          number = man1
          return
      end if
!
!---- prevent overflow
      exp1 = min(exp1,idint(log10(z1))-1)
!---- prevent underflow
      exp1 = max(exp1,idint(log10(z2)))
      num = 1
!
!---- calculate 10^|b|
      do 10 i = 1,abs(exp1)
          num = num*10
   10 end do
!
!---- b < 0
      if (exp1 <= -1) then
!-------- calculate x = a.10^b = a/10^(-b) = a/10^|b|
          number = man1/num
!---- b > 0
      else
!-------- calculate x = a.10^b = a.10^|b|
          number = man1*num
      end if
!
      end subroutine normal
!
!***********************************************************************
!
      subroutine addnum(man1,exp1,man2,exp2,man3,exp3)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      double precision man1,man2,man3
      integer exp1,exp2,exp3
!-----------------------------------------------------------------------
!     Function : Forms the sum of 2 mantissa-exponent numbers.
!-----------------------------------------------------------------------
!     MAN1 : mantissa of the first number in the summation
!     EXP1 : exponent of the first number in the summation
!     MAN2 : mantissa of the second number in the summation
!     EXP2 : exponent of the second number in the summation
!     MAN3 : mantissa of the sum
!     EXP3 : exponent of the sum
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      double precision tman1,tman2
      integer diff,i
!-----------------------------------------------------------------------
!---- if x = 0, then z = x+y = y
      if (man1 == 0) then
          man3 = man2
          exp3 = exp2
          return
!---- if y = 0, then z = x+y = x
      else if (man2 == 0) then
          man3 = man1
          exp3 = exp1
          return
      end if
!
      tman1 = man1
      tman2 = man2
!
!---- x <= y
      if (exp1 <= exp2) then
!-------- exponent of y - exponent of x; prevent underflow
          diff = min(exp2-exp1,idint(-log10(z2)))
!-------- exponent of z = exponent of y
          exp3 = exp2
!
!-------- scale the mantissa of y so that x and y have equal exponents
          do 10 i = 1,diff
              tman1 = tman1/10
   10     end do
!
!---- x > y
      else
!-------- exponent of x - exponent of y; prevent underflow
          diff = min(exp1-exp2,idint(-log10(z2)))
!-------- exponent of z = exponent of x
          exp3 = exp1
!
!-------- scale the mantissa of y so that x and y have equal exponents
          do 20 i = 1,diff
              tman2 = tman2/10
   20     end do
!
      end if
!
!---- mantissa of z = mantissa of x + mantissa of y
      man3 = tman1+tman2
!
!---- standard form of z = x+y
!-----***********
      call manexp(man3,exp3)
!-----***********
!
      return
!
      end subroutine addnum
!
!***********************************************************************
!
      subroutine emexp(number,man1,exp1)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      double precision number,man1
      integer exp1
!-----------------------------------------------------------------------
!     Function : The mantissa-exponent version of the function FEXP. It
!                takes exponentials in a devious way using the formula
!                exp(z) = 10^(z.log10(e))
!-----------------------------------------------------------------------
!     NUMBER : number to be exponentiated
!     MAN1   : mantissa of the exponent
!     EXP1   : exponent of the exponent (!) - well, you know what I mean
!-----------------------------------------------------------------------
      double precision prod
!-----------------------------------------------------------------------
!---- calculate x.log10(e)
      prod = number*log10(exp(1d0))
!---- calculate [x.log10(e)] = integer_part
      exp1 = idint(prod)
!---- calculate 10^(x.log10(e) - [x.log10(e)]) = 10^remainder
      man1 = 10**(prod-exp1)
!---- standard form is MAN1.10^EXP1 = 10^remainder.10^integer_part
!---- = 10^(remainder+integer_part) = 10^PROD = 10^(x.log10(e)) = exp(x)
!
      return
!
      end subroutine emexp
!
!***********************************************************************
!
      subroutine help(numcom)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer numcom
!-----------------------------------------------------------------------
!     Function : Activates the on-line help facility.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      call wrtscr('    The on-line help is no longer available')
      call wrtscr( &
      '    See the User Manual at http://www.sabre.lancs.ac.uk')
      numcom = 0
!
!---- command is information, so go straight to help facility
      if (numcom == 0) then
          return
      end if
!
!---- (help) info
      call wrtscr ( '   SABRE 2.0 was written in 1989 by Jon Barry, '// &
                       'Richard Davies & Brian Francis' )
      call wrtscr ( '   of the Centre for Applied Statistics, '// &
                       'Lancaster University, U.K.' )
      call newlns
      call wrtscr ( '   Version 3.x was written during '// &
                       '1994-1999 by Dave Stott,' )
      call wrtscr ( '   Brian Francis & Richard Davies of the Centre'// &
                      ' for Applied Statistics,' )
      call wrtscr ( '   Lancaster University, U.K.' )
      call newlns
      call wrtscr ( '   Versions 4-6 were written during 2005-2008 '// &
                       'by Rob Crouchley and the team at the,' )
      call wrtscr ( '   Centre for e-Science, Lancaster University, '// &
                       'U.K.' )
      call newlns
      call wrtscr ( '   SABRE is available for Sparc machines and '// &
                       'PCs running Windows 3.x, 95, NT or XP' )
      call newlns
      call wrtscr ( '   Telephone:    +44 (0) 1524 65201' )
      call wrtscr ( '   E-mail:       statistics@lancaster.ac.uk or '// &
                       'd.stott@lancaster.ac.uk' )
      call wrtscr ( '   Web page:     http://www.sabre.lancs.ac.uk')
      call wrtscr ( '   Mailing list: send an e-mail to '// &
                       'sabre-request@lists.lancs.ac.uk' )
      call wrtscr ( '                 containing the word "subscribe"'// &
                                    ' in the body of the message.' )
      call newlns
      call wrtscr ( '   Its development was initially funded by the '// &
                       'Economic and Social Research' )
      call wrtscr ( '   Council.' )
      call newlns
      call wrtscr ( '   The current release is version 6.0' )
      go to 1300
!
 1300 call newlns
!
      return
!
      end subroutine help
!
!***********************************************************************
!
      subroutine model(x,line,name,ncol,mnames,ni,iyvar,yname,bivar, &
                       rname,cname,nmes,ipos,work,maxcol,nm,beta,alias, &
                       cov,sca,ilfit,ilprev,con,alp,npar,nsub,endind, &
                       ilev,nvar,nilev,iend,icvar,xll,ilevdf,iredf, &
                       ienddf,tol,est0,est1,nmeil,niter,arith,itotdf, &
                       mode,y,ifail,offlag,offnme,inflag,xinit,link, &
                       risk,corr,robust,order,ncat,cutflag,xcut,sig_e, &
                       sig,sigflag,inormdf,n1var,family,trivar, &
                       n2var,univar,iname,rho,nlevel,irvar,depend, &
                       eqscale,dfirst,ndum,gamma,gammse,cquad,icutdf)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character endind,arith,link(3),family(3),cquad
      character(len=50) :: mnames(maxvar),name(maxvar),yname,rname, &
                cname(2),offnme(3),iname(3)
      integer maxcol,nmes,ncol,nlevel,ndum,nm(3),iend(2),ncat(3), &
              n2var,ni,ipos(maxvar),nsub(2),ilev(maxvar),nvar,nilev, &
              niter,ilfit,ilprev,ilevdf,iredf,ienddf,mode,nmeil,itotdf, &
              n1var,risk(nmes),npar,ncov,icutdf,oiend(2),maxcat,inormdf
      double precision x(nmes,maxcol),beta(maxpar*(maxpar+5)),y(4*nmes), &
             alp,con,xll,sca(3),alias(2*maxpar),sig(3),xinit(maxpar), &
             xcut(maxpar),sig_e(3),est0,est1,work(mxx-nmes*ncol), &
             cov(maxpar*(maxpar+1)/2),tol,rho(3),cind(nmes),gamma(maxy), &
             gammse(maxy)
      logical iyvar,bivar,icvar(2),ifail,offlag(3),sigflag(3),trivar, &
              inflag,corr,robust,cutflag,univar,irvar,order, &
              depend,eqscale,dfirst,ocorr
!-----------------------------------------------------------------------
!     Function : To fit the model. Initial parameter estimates are
!                obtained from logistic binary (or log-linear Poisson)
!                regression, ignoring the fact that measurements are
!                repetitions on the same individual. If the routine is
!                called via the FIT (as opposed to LFIT) command, then
!                the full logistic-normal/log-linear-normal mixture
!                model is fitted.
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      double precision quad(3,2,256),qloc(1373),qprb(1373),ymin(3), &
             oquad(3,2,256),ymax(3)
      integer k,i,j,idel(maxvar),ixpos,iypos,irpos,icpos(2),ilevt,l, &
              offpos(3),ninali,ixlev,prod(2*nmes),nest,it(2,nmes),maxit, &
              nnorm,ncorr,nre,ninal1,ninal2,n1lev,n2lev,onm(3),nsca, &
              rmax,r1,n1sub,clock1,icase,iobs,ncut,nend
      logical gmyes,facyes
!-----------------------------------------------------------------------
      data qloc(1:309) / &
!----   1 mass point
       0.00000000, &
!----   2 mass points
       1.00000000, &
!----   4 mass points
       2.33441422,  0.74196378, &
!----   6 mass points
       3.32425743,  1.88917588,  0.61670659, &
!----   8 mass points
       4.14454719,  2.80248586,  1.63651904,  0.53907981, &
!----  10 mass points
       4.85946283,  3.58182348,  2.48432584,  1.46598909,  0.48493571, &
!----  12 mass points
       5.50090170,  4.27182585,  3.22370983,  2.25946445,  1.34037520, &
       0.44440300, &
!----  14 mass points
       6.08740955,  4.89693640,  3.88692458,  2.96303658,  2.08834475, &
       1.24268896,  0.41259046, &
!----  16 mass points
       6.63087820,  5.47222571,  4.49295530,  3.60087362,  2.76024505, &
       1.95198035,  1.16382910,  0.38676060, &
!----  20 mass points
       7.61904854,  6.51059016,  5.57873881,  4.73458133,  3.94396735, &
       3.18901482,  2.45866361,  1.74524732,  1.04294535,  0.34696416, &
!----  24 mass points
       8.50780352,  7.43789067,  6.54167501,  5.73274718,  4.97804137, &
       4.26038361,  3.56930676,  2.89772864,  2.24046785,  1.59348043, &
       0.95342192,  0.31737010, &
!----  28 mass points
       9.32193781,  8.28306954,  7.41512529,  6.63373149,  5.90665633, &
       5.21722367,  4.55534038,  3.91425373,  3.28910697,  2.67620188, &
       2.07258267,  1.47578174,  0.88365256,  0.29425171, &
!----  32 mass points
      10.07742267,  9.06439921,  8.21972877,  7.46075575,  6.75593083, &
       6.08896431,  5.45003327,  4.83260461,  4.23202111,  3.64478125, &
       3.06813517,  2.49984042,  1.93800491,  1.38098020,  0.82728490, &
       0.27554642, &
!----  36 mass points
      10.78525331,  9.79427602,  8.96928653,  8.22911537,  7.54279304, &
       6.89434765,  6.27416833,  5.67588471,  5.09497851,  4.52807699, &
       3.97255734,  3.42630860,  2.88757970,  2.35487772,  1.82689658, &
       1.30246495,  0.78050649,  0.26000793, &
!----  40 mass points
      11.45337784, 10.48156053,  9.67355637,  8.94950454,  8.27894062, &
       7.64616376,  7.04173841,  6.45942338,  5.89480568,  5.34460545, &
       4.80628719,  4.27782616,  3.75755978,  3.24408873,  2.73620834, &
       2.23285922,  1.73309059,  1.23603200,  0.74087073,  0.24683290, &
!----  44 mass points
      12.08776071, 11.13284252, 10.33973351,  9.62974115,  8.97284870, &
       8.35359229,  7.76268639,  7.19399724,  6.64319640,  6.10707582, &
       5.58316483,  5.06950023,  4.56448030,  4.06676779,  3.57522300, &
       3.08885599,  2.60679146,  2.12824212,  1.65248799,  1.17885980, &
       0.70672524,  0.23547712, &
!----  48 mass points
      12.69301232, 11.75317846, 10.97330096, 10.27574158,  9.63088569, &
       9.02348028,  8.44437132,  7.88751732,  7.34866057,  6.82465132, &
       6.31306991,  5.81199999,  5.31988462,  4.83543089,  4.35754411, &
       3.88528112,  3.41781605,  2.95441469,  2.49441474,  2.03721030, &
       1.58223932,  1.12897323,  0.67690814,  0.22555707, &
!----  56 mass points
      13.83002175, 12.91614692, 12.15890778, 11.48250621, 10.85802651, &
      10.27057895,  9.71123280,  9.17410004,  8.65503764,  8.15098741, &
       7.65960702,  7.17904852,  6.70781798,  6.24468248,  5.78860610, &
       5.33870459,  4.89421236,  4.45445784,  4.01884475,  3.58683758, &
       3.15795009,  2.73173605,  2.30778176,  1.88569965,  1.46512304, &
       1.04570145,  0.62709651,  0.20897827, &
!----  64 mass points
      14.88618614, 13.99404991, 13.25564936, 12.59675248, 11.98903661, &
      11.41791806, 10.87465199, 10.35347592,  9.85033846,  9.36225255, &
       8.88693391,  8.42258409,  7.96775298,  7.52124766,  7.08206983, &
       6.64937146,  6.22242253,  5.80058710,  5.38330506,  4.97007811, &
       4.56045873,  4.15404137,  3.75045539,  3.34935918,  2.95043540, &
       2.55338687,  2.15793312,  1.76380738,  1.37075396,  0.97852591, &
       0.58688282,  0.19558891, &
!----  72 mass points
      15.87654063, 15.00303462, 14.28069917, 13.63666257, 13.04311651, &
      12.48574171, 11.95595282, 11.44809396, 10.95819048, 10.48331385, &
      10.02122653,  9.57016895,  9.12872447,  8.69572985,  8.27021384, &
       7.85135365,  7.43844330,  7.03087014,  6.62809706,  6.22964868, &
       5.83510071,  5.44407131,  5.05621428,  4.67121346,  4.28877809, &
       3.90863903,  3.53054548,  3.15426228,  2.77956754,  2.40625057, &
       2.03411013,  1.66295278,  1.29259148,  0.92284423,  0.55353292, &
       0.18448209, &
!----  80 mass points
      16.81197787, 15.95472489, 15.24634858, 14.61517739, 14.03385652, &
      13.48829932, 12.97006014, 12.47357593, 11.99493827, 11.53126851, &
      11.08036846, 10.64051073, 10.21030580,  9.78861405,  9.37448524, &
       8.96711570,  8.56581726,  8.16999410,  7.77912524,  7.39275113, &
       7.01046300,  6.63189458,  6.25671534,  5.88462505,  5.51534927, &
       5.14863571,  4.78425105,  4.42197838,  4.06161491,  3.70297012, &
       3.34586403,  2.99012582,  2.63559250,  2.28210779,  1.92952110, &
       1.57768663,  1.22646246,  0.87570982,  0.52529230,  0.17507520 /
      data qloc(310:725) / &
!----  88 mass points
      17.70068723, 16.85781509, 16.16174900, 15.54188953, 14.97128697, &
      14.43606214, 13.92789371, 13.44130085, 12.97243349, 12.51845662, &
      12.07720611, 11.64698208, 11.22641799, 10.81439384, 10.40997657, &
      10.01237792,  9.62092374,  9.23503121,  8.85419164,  8.47795716, &
       8.10593038,  7.73775611,  7.37311480,  7.01171715,  6.65329976, &
       6.29762150,  5.94446042,  5.59361130,  5.24488339,  4.89809863, &
       4.55309005,  4.20970035,  3.86778076,  3.52718994,  3.18779307, &
       2.84946100,  2.51206950,  2.17549859,  1.83963192,  1.50435619, &
       1.16956062,  0.83513646,  0.50097650,  0.16697463, &
!----  96 mass points
      18.54900896, 17.71900805, 17.03392905, 16.42414008, 15.86305703, &
      15.33698780, 14.83772298, 14.35985560, 13.89958774, 13.45412323, &
      13.02132803, 12.59952643, 12.18737180, 11.78376099, 11.38777557, &
      10.99864015, 10.61569213, 10.23835931,  9.86614278,  9.49860392, &
       9.13535408,  8.77604655,  8.42036997,  8.06804314,  7.71881063, &
       7.37243927,  7.02871513,  6.68744106,  6.34843453,  6.01152588, &
       5.67655674,  5.34337869,  5.01185214,  4.68184527,  4.35323316, &
       4.02589701,  3.69972343,  3.37460381,  3.05043374,  2.72711256, &
       2.40454280,  2.08262983,  1.76128143,  1.44040739,  1.11991923, &
       0.79972980,  0.47975303,  0.15990355, &
!---- 104 mass points
      19.36197013, 18.54359944, 17.86842433, 17.26769351, 16.71515562, &
      16.19728883, 15.70598546, 15.23590550, 14.78329825, 14.34540254, &
      13.92011135, 13.50577050, 13.10105101, 12.70486457, 12.31630545, &
      11.93460937, 11.55912362, 11.18928488, 10.82460244, 10.46464522, &
      10.10903172,  9.75742194,  9.40951102,  9.06502402,  8.72371165, &
       8.38534674,  8.04972134,  7.71664424,  7.38593891,  7.05744168, &
       6.73100030,  6.40647257,  6.08372523,  5.76263299,  5.44307762, &
       5.12494725,  4.80813565,  4.49254164,  4.17806857,  3.86462382, &
       3.55211839,  3.24046646,  2.92958509,  2.61939382,  2.30981444, &
       2.00077064,  1.69218778,  1.38399263,  1.07611314,  0.76847818, &
       0.46101739,  0.15366089, &
!---- 112 mass points
      20.14363610, 19.33585885, 18.66968592, 18.07717317, 17.53237438, &
      17.02192379, 16.53780688, 16.07474494, 15.62903036, 15.19793386, &
      14.77937293, 14.37171285, 13.97364053, 13.58408091, 13.20213954, &
      12.82706193, 12.45820400, 12.09501015, 11.73699664, 11.38373884, &
      11.03486120, 10.69002932, 10.34894366, 10.01133434,  9.67695698, &
       9.34558917,  9.01702764,  8.69108576,  8.36759157,  8.04638596, &
       7.72732122,  7.41025975,  7.09507294,  6.78164020,  6.46984811, &
       6.15958971,  5.85076379,  5.54327435,  5.23703008,  4.93194386, &
       4.62793238,  4.32491577,  4.02281721,  3.72156268,  3.42108064, &
       3.12130177,  2.82215876,  2.52358604,  2.22551962,  1.92789686, &
       1.63065631,  1.33373750,  1.03708081,  0.74062732,  0.44431859, &
       0.14809657, &
!---- 128 mass points
      21.62589891, 20.83679985, 20.18645795, 19.60836376, 19.07711363, &
      18.57961928, 18.10803180, 17.65718123, 17.22343632, 16.80412412, &
      16.39720529, 16.00107925, 15.61446060, 15.23629726, 14.86571422, &
      14.50197368, 14.14444616, 13.79258900, 13.44593006, 13.10405521, &
      12.76659856, 12.43323467, 12.10367235, 11.77764963, 11.45492963, &
      11.13529716, 10.81855588, 10.50452595, 10.19304198,  9.88395138, &
       9.57711285,  9.27239519,  8.96967613,  8.66884147,  8.36978421, &
       8.07240386,  7.77660576,  7.48230058,  7.18940378,  6.89783518, &
       6.60751855,  6.31838129,  6.03035407,  5.74337057,  5.45736718, &
       5.17228283,  4.88805869,  4.60463804,  4.32196605,  4.03998960, &
       3.75865718,  3.47791869,  3.19772532,  2.91802944,  2.63878447, &
       2.35994477,  2.08146552,  1.80330266,  1.52541273,  1.24775285, &
       0.97028058,  0.69295384,  0.41573086,  0.13857005, &
!---- 144 mass points
      23.01657536, 22.24353006, 21.60675497, 21.04098509, 20.52128764, &
      20.03481525, 19.57386266, 19.13335420, 18.70972473, 18.30034994, &
      17.90322765, 17.51678658, 17.13976499, 16.77113030, 16.41002389, &
      16.05572199, 15.70760723, 15.36514759, 15.02788037, 14.69539994, &
      14.36734807, 14.04340637, 13.72329014, 13.40674345, 13.09353508, &
      12.78345517, 12.47631249, 12.17193204, 11.87015313, 11.57082768, &
      11.27381883, 10.97899966, 10.68625214, 10.39546625, 10.10653913, &
       9.81937438,  9.53388147,  9.24997515,  8.96757501,  8.68660502, &
       8.40699315,  8.12867102,  7.85157360,  7.57563894,  7.30080790, &
       7.02702392,  6.75423283,  6.48238266,  6.21142344,  5.94130708, &
       5.67198720,  5.40341901,  5.13555916,  4.86836567,  4.60179779, &
       4.33581591,  4.07038149,  3.80545692,  3.54100549,  3.27699131, &
       3.01337922,  2.75013470,  2.48722388,  2.22461339,  1.96227037, &
       1.70016237,  1.43825733,  1.17652348,  0.91492934,  0.65344364, &
       0.39203528,  0.13067328, &
!---- 160 mass points
      24.33074285, 23.57173605, 22.94679647, 22.39175553, 21.88209561, &
      21.40518174, 20.95343534, 20.52186386, 20.10696035, 19.70614328, &
      19.31744313, 18.93931425, 18.57051556, 18.21003148, 17.85701762, &
      17.51076226, 17.17065847, 16.83618330, 16.50688206, 16.18235626, &
      15.86225408, 15.54626296, 15.23410352, 14.92552472, 14.62029991, &
      14.31822348, 14.01910819, 13.72278282, 13.42909029, 13.13788595, &
      12.84903624, 12.56241742, 12.27791459, 11.99542071, 11.71483589, &
      11.43606661, 11.15902521, 10.88362925, 10.60980113, 10.33746759, &
      10.06655937,  9.79701088,  9.52875988,  9.26174721,  8.99591655, &
       8.73121418,  8.46758882,  8.20499139,  7.94337490,  7.68269426, &
       7.42290615,  7.16396889,  6.90584233,  6.64848773,  6.39186768, &
       6.13594598,  5.88068758,  5.62605847,  5.37202564,  5.11855700, &
       4.86562128,  4.61318802,  4.36122749,  4.10971063,  3.85860900, &
       3.60789475,  3.35754056,  3.10751960,  2.85780547,  2.60837219, &
       2.35919415,  2.11024607,  1.86150297,  1.61294012,  1.36453303, &
       1.11625742,  0.86808915,  0.62000423,  0.37197878,  0.12398897 /
      data qloc(726:1125) / &
!---- 176 mass points
      25.57975196, 24.83319122, 24.21872379, 23.67315934, 23.17235084, &
      22.70385299, 22.26020036, 21.83647458, 21.42922065, 21.03589518, &
      20.65455773, 20.28368548, 19.92205568, 19.56866782, 19.22269005, &
      18.88342132, 18.55026383, 18.22270256, 17.90028981, 17.58263325, &
      17.26938663, 16.96024233, 16.65492552, 16.35318932, 16.05481087, &
      15.75958812, 15.46733713, 15.17788977, 14.89109188, 14.60680159, &
      14.32488794, 14.04522975, 13.76771450, 13.49223751, 13.21870114, &
      12.94701407, 12.67709077, 12.40885092, 12.14221895, 11.87712365, &
      11.61349775, 11.35127766, 11.09040309, 10.83081684, 10.57246456, &
      10.31529448, 10.05925727,  9.80430582,  9.55039513,  9.29748208, &
       9.04552536,  8.79448533,  8.54432390,  8.29500440,  8.04649154, &
       7.79875127,  7.55175071,  7.30545809,  7.05984266,  6.81487463, &
       6.57052512,  6.32676607,  6.08357023,  5.84091106,  5.59876273, &
       5.35710006,  5.11589844,  4.87513386,  4.63478282,  4.39482230, &
       4.15522975,  3.91598305,  3.67706046,  3.43844062,  3.20010250, &
       2.96202537,  2.72418880,  2.48657263,  2.24915691,  2.01192193, &
       1.77484815,  1.53791622,  1.30110694,  1.06440120,  0.82778005, &
       0.59122460,  0.35471602,  0.11823555, &
!---- 192 mass points
      26.77240039, 26.03699800, 25.43190350, 24.89480830, 24.40189961, &
      23.94090324, 23.50445584, 23.08770598, 22.68724550, 22.30056542, &
      21.92575150, 21.56130144, 21.20600900, 20.85888718, 20.51911540, &
      20.18600211, 19.85895767, 19.53747414, 19.22111000, 18.90947839, &
      18.60223791, 18.29908534, 17.99974979, 17.70398797, 17.41158035, &
      17.12232790, 16.83604951, 16.55257968, 16.27176671, 15.99347104, &
      15.71756392, 15.44392622, 15.17244743, 14.90302476, 14.63556237, &
      14.36997073, 14.10616600, 13.84406951, 13.58360729, 13.32470970, &
      13.06731103, 12.81134917, 12.55676534, 12.30350380, 12.05151164, &
      11.80073855, 11.55113661, 11.30266017, 11.05526561, 10.80891126, &
      10.56355722, 10.31916528, 10.07569878,  9.83312250,  9.59140259, &
       9.35050647,  9.11040274,  8.87106113,  8.63245240,  8.39454829, &
       8.15732148,  7.92074550,  7.68479468,  7.44944415,  7.21466972, &
       6.98044791,  6.74675585,  6.51357128,  6.28087251,  6.04863837, &
       5.81684819,  5.58548178,  5.35451938,  5.12394165,  4.89372963, &
       4.66386474,  4.43432873,  4.20510368,  3.97617196,  3.74751623, &
       3.51911940,  3.29096463,  3.06303530,  2.83531500,  2.60778750, &
       2.38043676,  2.15324689,  1.92620214,  1.69928691,  1.47248568, &
       1.24578306,  1.01916372,  0.79261244,  0.56611402,  0.33965333, &
       0.11321526, &
!---- 208 mass points
      27.91567043, 27.19036527, 26.59374009, 26.06428798, 25.57850013, &
      25.12425806, 24.69429128, 24.28380930, 23.88944658, 23.50872545, &
      23.13975547, 22.78105302, 22.43142683, 22.08990215, 21.75566861, &
      21.42804329, 21.10644390, 20.79036887, 20.47938225, 20.17310208, &
      19.87119131, 19.57335063, 19.27931268, 18.98883735, 18.70170802, &
      18.41772835, 18.13671968, 17.85851883, 17.58297620, 17.30995424, &
      17.03932610, 16.77097440, 16.50479032, 16.24067267, 15.97852714, &
      15.71826566, 15.45980579, 15.20307022, 14.94798631, 14.69448567, &
      14.44250383, 14.19197990, 13.94285626, 13.69507835, 13.44859439, &
      13.20335517, 12.95931391, 12.71642602, 12.47464898, 12.23394219, &
      11.99426682, 11.75558573, 11.51786330, 11.28106539, 11.04515921, &
      10.81011324, 10.57589715, 10.34248174, 10.10983885,  9.87794132, &
       9.64676291,  9.41627825,  9.18646282,  8.95729284,  8.72874530, &
       8.50079785,  8.27342882,  8.04661714,  7.82034232,  7.59458442, &
       7.36932405,  7.14454226,  6.92022062,  6.69634109,  6.47288607, &
       6.24983836,  6.02718112,  5.80489786,  5.58297242,  5.36138896, &
       5.14013193,  4.91918607,  4.69853636,  4.47816804,  4.25806659, &
       4.03821769,  3.81860724,  3.59922131,  3.38004618,  3.16106827, &
       2.94227417,  2.72365059,  2.50518439,  2.28686256,  2.06867217, &
       1.85060042,  1.63263457,  1.41476199,  1.19697010,  0.97924639, &
       0.76157838,  0.54395367,  0.32635985,  0.10878456, &
!---- 224 mass points
      29.01521347, 28.29911716, 27.71020644, 27.18770851, 26.70839310, &
      26.26028438, 25.83619702, 25.43139636, 25.04255596, 24.66722680, &
      24.30354027, 23.95002986, 23.60551800, 23.26904117, 22.93979830, &
      22.61711434, 22.30041373, 21.98920070, 21.68304435, 21.38156718, &
      21.08443608, 20.79135528, 20.50206058, 20.21631474, 19.93390377, &
      19.65463371, 19.37832810, 19.10482580, 18.83397909, 18.56565220, &
      18.29971990, 18.03606639, 17.77458430, 17.51517381, 17.25774196, &
      17.00220191, 16.74847242, 16.49647734, 16.24614513, 15.99740847, &
      15.75020393, 15.50447161, 15.26015487, 15.01720010, 14.77555643, &
      14.53517559, 14.29601166, 14.05802093, 13.82116175, 13.58539437, &
      13.35068080, 13.11698471, 12.88427134, 12.65250736, 12.42166078, &
      12.19170091, 11.96259822, 11.73432432, 11.50685188, 11.28015452, &
      11.05420684, 10.82898428, 10.60446313, 10.38062046, 10.15743408, &
       9.93488248,  9.71294482,  9.49160091,  9.27083111,  9.05061637, &
       8.83093816,  8.61177844,  8.39311968,  8.17494477,  7.95723704, &
       7.73998023,  7.52315848,  7.30675626,  7.09075843,  6.87515016, &
       6.65991692,  6.44504451,  6.23051898,  6.01632668,  5.80245419, &
       5.58888834,  5.37561620,  5.16262502,  4.94990230,  4.73743570, &
       4.52521309,  4.31322247,  4.10145206,  3.88989019,  3.67852534, &
       3.46734614,  3.25634133,  3.04549978,  2.83481046,  2.62426244, &
       2.41384490,  2.20354708,  1.99335832,  1.78326802,  1.57326565, &
       1.36334072,  1.15348282,  0.94368157,  0.73392660,  0.52420761, &
       0.31451430,  0.10483639 /
      data qloc(1126:1373) / &
!---- 240 mass points
      30.07567994, 29.36803867, 28.78620308, 28.27007685, 27.79668645, &
      27.35418761, 26.93547353, 26.53586104, 26.15205975, 25.78164710, &
      25.42277461, 25.07399151, 24.73413288, 24.40224553, 24.07753697, &
      23.75933938, 23.44708339, 23.14027853, 22.83849857, 22.54137007, &
      22.24856354, 21.95978642, 21.67477740, 21.39330187, 21.11514820, &
      20.84012459, 20.56805659, 20.29878487, 20.03216342, 19.76805804, &
      19.50634498, 19.24690982, 18.98964648, 18.73445638, 18.48124769, &
      18.22993469, 17.98043718, 17.73268002, 17.48659262, 17.24210858, &
      16.99916535, 16.75770390, 16.51766840, 16.27900605, 16.04166677, &
      15.80560303, 15.57076966, 15.33712369, 15.10462415, 14.87323200, &
      14.64290994, 14.41362232, 14.18533504, 13.95801541, 13.73163213, &
      13.50615512, 13.28155553, 13.05780557, 12.83487856, 12.61274877, &
      12.39139140, 12.17078255, 11.95089912, 11.73171883, 11.51322010, &
      11.29538208, 11.07818456, 10.86160798, 10.64563337, 10.43024230, &
      10.21541690, 10.00113982,  9.78739416,  9.57416349,  9.36143184, &
       9.14918362,  8.93740367,  8.72607717,  8.51518967,  8.30472709, &
       8.09467562,  7.88502181,  7.67575247,  7.46685470,  7.25831587, &
       7.05012361,  6.84226578,  6.63473047,  6.42750600,  6.22058089, &
       6.01394386,  5.80758383,  5.60148989,  5.39565129,  5.19005746, &
       4.98469799,  4.77956259,  4.57464113,  4.36992361,  4.16540013, &
       3.96106094,  3.75689638,  3.55289690,  3.34905305,  3.14535545, &
       2.94179485,  2.73836203,  2.53504787,  2.33184332,  2.12873938, &
       1.92572711,  1.72279764,  1.51994211,  1.31715174,  1.11441776, &
       0.91173145,  0.70908409,  0.50646702,  0.30387157,  0.10128908, &
!---- 256 mass points
      31.10095104, 30.40111778, 29.82580946, 29.31555650, 28.84762330, &
      28.41028757, 27.99651529, 27.60167108, 27.22249796, 26.85659791, &
      26.50214110, 26.15769138, 25.82209555, 25.49440997, 25.17385011, &
      24.85975486, 24.55156053, 24.24878161, 23.95099612, 23.65783441, &
      23.36897031, 23.08411422, 22.80300749, 22.52541791, 22.25113601, &
      21.97997201, 21.71175323, 21.44632202, 21.18353395, 20.92325621, &
      20.66536640, 20.40975133, 20.15630610, 19.90493323, 19.65554193, &
      19.40804744, 19.16237052, 18.91843689, 18.67617682, 18.43552471, &
      18.19641880, 17.95880078, 17.72261557, 17.48781105, 17.25433780, &
      17.02214896, 16.79119999, 16.56144852, 16.33285421, 16.10537858, &
      15.87898493, 15.65363817, 15.42930474, 15.20595252, 14.98355072, &
      14.76206981, 14.54148145, 14.32175840, 14.10287447, 13.88480443, &
      13.66752402, 13.45100982, 13.23523926, 13.02019052, 12.80584255, &
      12.59217499, 12.37916815, 12.16680294, 11.95506089, 11.74392410, &
      11.53337519, 11.32339731, 11.11397408, 10.90508959, 10.69672836, &
      10.48887534, 10.28151587, 10.07463568,  9.86822085,  9.66225782, &
       9.45673336,  9.25163453,  9.04694871,  8.84266357,  8.63876705, &
       8.43524735,  8.23209290,  8.02929242,  7.82683480,  7.62470918, &
       7.42290491,  7.22141152,  7.02021876,  6.81931652,  6.61869489, &
       6.41834413,  6.21825465,  6.01841700,  5.81882188,  5.61946015, &
       5.42032277,  5.22140084,  5.02268557,  4.82416828,  4.62584042, &
       4.42769352,  4.22971920,  4.03190920,  3.83425532,  3.63674945, &
       3.43938356,  3.24214969,  3.04503994,  2.84804649,  2.65116155, &
       2.45437743,  2.25768644,  2.06108097,  1.86455345,  1.66809633, &
       1.47170211,  1.27536332,  1.07907252,  0.88282227,  0.68660519, &
       0.49041388,  0.29424097,  0.09807910/
!-----------------------------------------------------------------------
      data qprb(1:201) / &
!----   1 mass point
      0.10000000d+001, &
!----   2 mass points
      0.50000000d+000, &
!----   4 mass points
      0.45875855d-001, 0.45412415d+000, &
!----   6 mass points
      0.25557844d-002, 0.88615746d-001, 0.40882847d+000, &
!----   8 mass points
      0.11261454d-003, 0.96352201d-002, 0.11723991d+000, &
      0.37301226d+000, &
!----  10 mass points
      0.43106526d-005, 0.75807093d-003, 0.19111581d-001, &
      0.13548370d+000, 0.34464233d+000, &
!----  12 mass points
      0.14999272d-006, 0.48371849d-004, 0.22033807d-002, &
      0.29116688d-001, 0.14696705d+000, 0.32166436d+000, &
!----  14 mass points
      0.48681613d-008, 0.26609913d-005, 0.20033955d-003, &
      0.44289191d-002, 0.38650109d-001, 0.15408334d+000, &
      0.30263463d+000, &
!----  16 mass points
      0.14978147d-009, 0.13094732d-006, 0.15300032d-004, &
      0.52598493d-003, 0.72669376d-002, 0.47284752d-001, &
      0.15833837d+000, 0.28656852d+000, &
!----  20 mass points
      0.12578007d-012, 0.24820624d-009, 0.61274903d-007, &
      0.44021211d-005, 0.12882628d-003, 0.18301031d-002, &
      0.13997837d-001, 0.61506372d-001, 0.16173933d+000, &
      0.26079306d+000, &
!----  24 mass points
      0.93901937d-016, 0.37149742d-012, 0.17186649d-009, &
      0.22674617d-007, 0.12176597d-005, 0.32095006d-004, &
      0.46471872d-003, 0.39766089d-002, 0.21126344d-001, &
      0.72069364d-001, 0.16145951d+000, 0.24087012d+000, &
!----  28 mass points
      0.64325474d-019, 0.46917656d-015, 0.37459010d-012, &
      0.83266098d-010, 0.74793626d-008, 0.33048644d-006, &
      0.80935841d-005, 0.11882854d-003, 0.11043059d-002, &
      0.67524597d-002, 0.27935785d-001, 0.79773366d-001, &
      0.15941819d+000, 0.22488863d+000, &
!----  32 mass points
      0.41246075d-022, 0.52084496d-018, 0.67552902d-015, &
      0.23780649d-012, 0.33475012d-010, 0.23125184d-008, &
      0.88812907d-007, 0.20596221d-005, 0.30559803d-004, &
      0.30255703d-003, 0.20620511d-002, 0.99034617d-002, &
      0.34109848d-001, 0.85344808d-001, 0.15653899d+000, &
      0.21170557d+000, &
!----  36 mass points
      0.25090376d-025, 0.52223662d-021, 0.10502942d-017, &
      0.55871140d-015, 0.11740293d-012, 0.12047446d-010, &
      0.68710834d-009, 0.23738228d-007, 0.52783152d-006, &
      0.78969781d-005, 0.82200256d-004, 0.61075483d-003, &
      0.33041345d-002, 0.13216574d-001, 0.39552370d-001, &
      0.89342498d-001, 0.15329101d+000, 0.20059201d+000, &
!----  40 mass points
      0.14618399d-028, 0.48204679d-024, 0.14486094d-020, &
      0.11222752d-017, 0.33898534d-015, 0.49680885d-013, &
      0.40376386d-011, 0.19891185d-009, 0.63258972d-008, &
      0.13603424d-006, 0.20488974d-005, 0.22211771d-004, &
      0.17707293d-003, 0.10558790d-002, 0.47735449d-002, &
      0.16537844d-001, 0.44274555d-001, 0.92176579d-001, &
      0.14992111d+000, 0.19105901d+000, &
!----  44 mass points
      0.82152509d-032, 0.41536349d-027, 0.18094811d-023, &
      0.19815156d-020, 0.83466770d-018, 0.16934271d-015, &
      0.18985180d-013, 0.12890648d-011, 0.56565652d-010, &
      0.16830458d-008, 0.35220801d-007, 0.53341522d-006, &
      0.59805065d-005, 0.50550243d-004, 0.32690297d-003, &
      0.16368797d-002, 0.64079902d-002, 0.19765681d-001, &
      0.48334227d-001, 0.94144717d-001, 0.14656144d+000, &
      0.18276506d+000, &
!----  48 mass points
      0.44771555d-035, 0.33764561d-030, 0.20790590d-026, &
      0.31394766d-023, 0.17988549d-020, 0.49254651d-018, &
      0.74199936d-016, 0.67566773d-014, 0.39758060d-012, &
      0.15883610d-010, 0.44742872d-009, 0.91540557d-008, &
      0.13927917d-006, 0.16063937d-005, 0.14266092d-004, &
      0.98818049d-004, 0.53958658d-003, 0.23430821d-002, &
      0.81496969d-002, 0.22838212d-001, 0.51805184d-001, &
      0.95463401d-001, 0.14328246d+000, 0.17546354d+000, &
!----  56 mass points
      0.12321604d-041, 0.19323403d-036, 0.22354061d-032, &
      0.60469805d-029, 0.60406163d-026, 0.28362374d-023, &
      0.72528624d-021, 0.11146506d-018, 0.11040774d-016, &
      0.74239232d-015, 0.35264363d-013, 0.12210149d-011, &
      0.31602503d-010, 0.62410691d-009, 0.95654531d-008, &
      0.11540196d-006, 0.11090576d-005, 0.85762918d-005, &
      0.53821424d-004, 0.27609911d-003, 0.11649091d-002, &
      0.40633228d-002, 0.11768190d-001, 0.28400399d-001, &
      0.57276184d-001, 0.96745677d-001, 0.13709202d+000, &
      0.16314957d+000 /
      data qprb(202:453) / &
!----  64 mass points
      0.31231880d-048, 0.94769632d-043, 0.19301704d-038, &
      0.87866357d-035, 0.14384921d-031, 0.10883802d-028, &
      0.44355444d-026, 0.10785651d-023, 0.16829001d-021, &
      0.17784692d-019, 0.13269089d-017, 0.72221536d-016, &
      0.29442931d-014, 0.91869288d-013, 0.22337269d-011, &
      0.42964262d-010, 0.66214234d-009, 0.82660844d-008, &
      0.84376410d-007, 0.70994246d-006, 0.49583797d-005, &
      0.28919958d-004, 0.14160239d-003, 0.58468608d-003, &
      0.20438258d-002, 0.60684460d-002, 0.15347722d-001, &
      0.33140486d-001, 0.61213639d-001, 0.96863364d-001, &
      0.13145323d+000, 0.15310832d+000, &
!----  72 mass points
      0.74135643d-055, 0.41123079d-049, 0.14016288d-044, &
      0.10230108d-040, 0.26173963d-037, 0.30435097d-034, &
      0.18845623d-031, 0.69079701d-029, 0.16160953d-026, &
      0.25517950d-024, 0.28389893d-022, 0.23023525d-020, &
      0.13988496d-018, 0.65122528d-017, 0.23669459d-015, &
      0.68233622d-014, 0.15812368d-012, 0.29798321d-011, &
      0.46122267d-010, 0.59144493d-009, 0.63311995d-008, &
      0.56951075d-007, 0.43300012d-006, 0.27968283d-005, &
      0.15416641d-004, 0.72807840d-004, 0.29562300d-003, &
      0.10351118d-002, 0.31337909d-002, 0.82219159d-002, &
      0.18729950d-001, 0.37107748d-001, 0.64022449d-001, &
      0.96292561d-001, 0.12635187d+000, 0.14471746d+000, &
!----  80 mass points
      0.16676172d-061, 0.16148147d-055, 0.88417337d-051, &
      0.99534648d-047, 0.38322200d-043, 0.65965481d-040, &
      0.59776715d-037, 0.31804454d-034, 0.10736475d-031, &
      0.24360676d-029, 0.38835748d-027, 0.45051594d-025, &
      0.39121582d-023, 0.26028723d-021, 0.13528365d-019, &
      0.55835186d-018, 0.18557109d-016, 0.50266251d-015, &
      0.11213386d-013, 0.20789539d-012, 0.32290148d-011, &
      0.42312639d-010, 0.47070652d-009, 0.44700308d-008, &
      0.36415408d-007, 0.25560810d-006, 0.15519302d-005, &
      0.81787428d-005, 0.37528431d-004, 0.15034405d-003, &
      0.52713231d-003, 0.16210281d-002, 0.43803902d-002, &
      0.10418184d-001, 0.21838989d-001, 0.40396396d-001, &
      0.65999311d-001, 0.95313043d-001, 0.12173798d+000, &
      0.13756965d+000, &
!----  88 mass points
      0.35859827d-068, 0.58341380d-062, 0.49613692d-057, &
      0.83446225d-053, 0.46876562d-049, 0.11586562d-045, &
      0.14906304d-042, 0.11166526d-039, 0.52748649d-037, &
      0.16671728d-034, 0.36899814d-032, 0.59291595d-030, &
      0.71208587d-028, 0.65469945d-026, 0.47011487d-024, &
      0.26813749d-022, 0.12324609d-020, 0.46222920d-019, &
      0.14299102d-017, 0.36833105d-016, 0.79664361d-015, &
      0.14574255d-013, 0.22701173d-012, 0.30282025d-011, &
      0.34774709d-010, 0.34539564d-009, 0.29796900d-008, &
      0.22411391d-007, 0.14746497d-006, 0.85145700d-006, &
      0.43260231d-005, 0.19388514d-004, 0.76823915d-004, &
      0.26965549d-003, 0.83995002d-003, 0.23254747d-002, &
      0.57303903d-002, 0.12583289d-001, 0.24648529d-001, &
      0.43107633d-001, 0.67359030d-001, 0.94094249d-001, &
      0.11755464d+000, 0.13138560d+000, &
!----  96 mass points
      0.74209952d-075, 0.19638545d-068, 0.25211954d-063, &
      0.61687500d-059, 0.49274811d-055, 0.17051468d-051, &
      0.30371581d-048, 0.31239543d-045, 0.20135580d-042, &
      0.86422502d-040, 0.25881642d-037, 0.56118810d-035, &
      0.90770951d-033, 0.11224816d-030, 0.10832275d-028, &
      0.83005916d-027, 0.51262072d-025, 0.25843274d-023, &
      0.10754921d-021, 0.37309582d-020, 0.10882536d-018, &
      0.26893965d-017, 0.56695702d-016, 0.10257921d-014, &
      0.16016030d-013, 0.21685652d-012, 0.25576425d-011, &
      0.26381564d-010, 0.23885495d-009, 0.19044697d-008, &
      0.13412853d-007, 0.83667589d-007, 0.46339963d-006, &
      0.22839710d-005, 0.10037952d-004, 0.39411082d-004, &
      0.13846233d-003, 0.43594693d-003, 0.12317003d-002, &
      0.31265306d-002, 0.71377741d-002, 0.14669121d-001, &
      0.27160034d-001, 0.45334862d-001, 0.68257414d-001, &
      0.92741275d-001, 0.11374796d+000, 0.12596662d+000, &
!---- 104 mass points
      0.14856824d-081, 0.62197723d-075, 0.11765027d-069, &
      0.40935565d-065, 0.45488728d-061, 0.21571308d-057, &
      0.52080233d-054, 0.72017204d-051, 0.62013603d-048, &
      0.35385005d-045, 0.14034696d-042, 0.40184683d-040, &
      0.85637771d-038, 0.13929701d-035, 0.17661104d-033, &
      0.17767201d-031, 0.14399777d-029, 0.95266134d-028, &
      0.52040437d-026, 0.23709709d-024, 0.90896384d-023, &
      0.29554393d-021, 0.82073665d-020, 0.19589628d-018, &
      0.40415564d-017, 0.72442124d-016, 0.11333483d-014, &
      0.15541368d-013, 0.18751142d-012, 0.19975177d-011, &
      0.18847671d-010, 0.15797544d-009, 0.11793402d-008, &
      0.78606959d-008, 0.46883533d-007, 0.25072630d-006, &
      0.12044997d-005, 0.52068790d-005, 0.20285416d-004, &
      0.71324165d-004, 0.22661545d-003, 0.65139195d-003, &
      0.16956826d-002, 0.40012637d-002, 0.85655523d-002, &
      0.16646783d-001, 0.29389412d-001, 0.47159049d-001, &
      0.68808320d-001, 0.91320463d-001, 0.11027014d+000, &
      0.12116700d+000 /
      data qprb(454:725) / &
!---- 112 mass points
      0.28894199d-088, 0.18679162d-081, 0.50974647d-076, &
      0.24731793d-071, 0.37518905d-067, 0.23937124d-063, &
      0.76927779d-060, 0.14045907d-056, 0.15870142d-053, &
      0.11823664d-050, 0.60992883d-048, 0.22643129d-045, &
      0.62414116d-043, 0.13106266d-040, 0.21421628d-038, &
      0.27752360d-036, 0.28945501d-034, 0.24634176d-032, &
      0.17308611d-030, 0.10144333d-028, 0.50046244d-027, &
      0.20951550d-025, 0.74970890d-024, 0.23078851d-022, &
      0.61477339d-021, 0.14245798d-019, 0.28854121d-018, &
      0.51306092d-017, 0.80407057d-016, 0.11147033d-014, &
      0.13715356d-013, 0.15023146d-012, 0.14690477d-011, &
      0.12857289d-010, 0.10095536d-009, 0.71272519d-009, &
      0.45331387d-008, 0.26023301d-007, 0.13506743d-006, &
      0.63480877d-006, 0.27055970d-005, 0.10470886d-004, &
      0.36840685d-004, 0.11797066d-003, 0.34415760d-003, &
      0.91552841d-003, 0.22226572d-002, 0.49280735d-002, &
      0.99854577d-002, 0.18500886d-001, 0.31359216d-001, &
      0.48648334d-001, 0.69095748d-001, 0.89874162d-001, &
      0.10707987d+000, 0.11687712d+000, &
!---- 128 mass points
      0.10150143d-101, 0.14715037d-094, 0.79251134d-089, &
      0.71158382d-084, 0.19189386d-079, 0.21163969d-075, &
      0.11518414d-071, 0.35060118d-068, 0.65229070d-065, &
      0.79238580d-062, 0.66121806d-059, 0.39454189d-056, &
      0.17388760d-053, 0.58138961d-051, 0.15079342d-048, &
      0.30918814d-046, 0.50935262d-044, 0.68367339d-042, &
      0.75685891d-040, 0.69851498d-038, 0.54256465d-036, &
      0.35769132d-034, 0.20166273d-032, 0.97892705d-031, &
      0.41166347d-029, 0.15080357d-027, 0.48367017d-026, &
      0.13644380d-024, 0.33997974d-023, 0.75114388d-022, &
      0.14767423d-020, 0.25918913d-019, 0.40735093d-018, &
      0.57487151d-017, 0.73037059d-016, 0.83738363d-015, &
      0.86831583d-014, 0.81601409d-013, 0.69633096d-012, &
      0.54051148d-011, 0.38228248d-010, 0.24673064d-009, &
      0.14552690d-008, 0.78545950d-008, 0.38842090d-007, &
      0.17618924d-006, 0.73386791d-006, 0.28096215d-005, &
      0.98961594d-005, 0.32095300d-004, 0.95920178d-004, &
      0.26435190d-003, 0.67226789d-003, 0.15785118d-002, &
      0.34239890d-002, 0.68644489d-002, 0.12724907d-001, &
      0.21819442d-001, 0.34619080d-001, 0.50838378d-001, &
      0.69115071d-001, 0.87003921d-001, 0.10142610d+000, &
      0.10950785d+000, &
!---- 144 mass points
      0.32844601d-115, 0.99785274d-108, 0.10003433d-101, &
      0.15746098d-096, 0.71684361d-092, 0.12999701d-087, &
      0.11408011d-083, 0.55151526d-080, 0.16103364d-076, &
      0.30405116d-073, 0.39125048d-070, 0.35765842d-067, &
      0.24019098d-064, 0.12181830d-061, 0.47748650d-059, &
      0.14750040d-056, 0.36515861d-054, 0.73504970d-052, &
      0.12183911d-049, 0.16815732d-047, 0.19515210d-045, &
      0.19211444d-043, 0.16168825d-041, 0.11716257d-039, &
      0.73562435d-038, 0.40252034d-036, 0.19295988d-034, &
      0.81428937d-033, 0.30383180d-031, 0.10064336d-029, &
      0.29706429d-028, 0.78400879d-027, 0.18559932d-025, &
      0.39527162d-024, 0.75939190d-023, 0.13194453d-021, &
      0.20782609d-020, 0.29740662d-019, 0.38746992d-018, &
      0.46046693d-017, 0.50005083d-016, 0.49706964d-015, &
      0.45299502d-014, 0.37904023d-013, 0.29160519d-012, &
      0.20653235d-011, 0.13483217d-010, 0.81228470d-010, &
      0.45206204d-009, 0.23264799d-008, 0.11082072d-007, &
      0.48904152d-007, 0.20009287d-006, 0.75964900d-006, &
      0.26779458d-005, 0.87717960d-005, 0.26714263d-004, &
      0.75686097d-004, 0.19958993d-003, 0.49014426d-003, &
      0.11214180d-002, 0.23913763d-002, 0.47547434d-002, &
      0.88176103d-002, 0.15256312d-001, 0.24634170d-001, &
      0.37129225d-001, 0.52247861d-001, 0.68654088d-001, &
      0.84249201d-001, 0.96562553d-001, 0.10337683d+000, &
!---- 160 mass points
      0.99358171d-129, 0.59860021d-121, 0.10652214d-114, &
      0.28135652d-109, 0.20742508d-104, 0.59413938d-100, &
      0.80834812d-096, 0.59717328d-092, 0.26339002d-088, &
      0.74419069d-085, 0.14219063d-081, 0.19175333d-078, &
      0.18893683d-075, 0.13994321d-072, 0.79795759d-070, &
      0.35739670d-067, 0.12792607d-064, 0.37144031d-062, &
      0.88633464d-060, 0.17581633d-057, 0.29287116d-055, &
      0.41340089d-053, 0.49848856d-051, 0.51723520d-049, &
      0.46486469d-047, 0.36405189d-045, 0.24978140d-043, &
      0.15089538d-041, 0.80628836d-040, 0.38266949d-038, &
      0.16194019d-036, 0.61324489d-035, 0.20849641d-033, &
      0.63838430d-032, 0.17653267d-030, 0.44206148d-029, &
      0.10049237d-027, 0.20786620d-026, 0.39208285d-025, &
      0.67576919d-024, 0.10662870d-022, 0.15430610d-021, &
      0.20514299d-020, 0.25094741d-019, 0.28288429d-018, &
      0.29426935d-017, 0.28285655d-016, 0.25154437d-015, &
      0.20720561d-014, 0.15827398d-013, 0.11222613d-012, &
      0.73940897d-012, 0.45309297d-011, 0.25845478d-010, &
      0.13735282d-009, 0.68059204d-009, 0.31466862d-008, &
      0.13584374d-007, 0.54793627d-007, 0.20662946d-006, &
      0.72891607d-006, 0.24066978d-005, 0.74412344d-005, &
      0.21555297d-004, 0.58524930d-004, 0.14899935d-003, &
      0.35583625d-003, 0.79742774d-003, 0.16774461d-002, &
      0.33132229d-002, 0.61463122d-002, 0.10711396d-001, &
      0.17540396d-001, 0.26994673d-001, 0.39051224d-001, &
      0.53109365d-001, 0.67910800d-001, 0.81654267d-001, &
      0.92325549d-001, 0.98172150d-001 /
      data qprb(726:1013) / &
!---- 176 mass points
      0.28407918d-142, 0.32413023d-134, 0.98428879d-128, &
      0.42072302d-122, 0.48533684d-117, 0.21243221d-112, &
      0.43387881d-108, 0.47455905d-104, 0.30646262d-100, &
      0.12562925d-096, 0.34562876d-093, 0.66686676d-090, &
      0.93500343d-087, 0.98092232d-084, 0.78907812d-081, &
      0.49689036d-078, 0.24932044d-075, 0.10122186d-072, &
      0.33699502d-070, 0.93093074d-068, 0.21561644d-065, &
      0.42261846d-063, 0.70685139d-061, 0.10164154d-058, &
      0.12650769d-056, 0.13713040d-054, 0.13018286d-052, &
      0.10879309d-050, 0.80411259d-049, 0.52793207d-047, &
      0.30911292d-045, 0.16200799d-043, 0.76264297d-042, &
      0.32348313d-040, 0.12399729d-038, 0.43072591d-037, &
      0.13593691d-035, 0.39072280d-034, 0.10251243d-032, &
      0.24602648d-031, 0.54119041d-030, 0.10931894d-028, &
      0.20313491d-027, 0.34780862d-026, 0.54959752d-025, &
      0.80267908d-024, 0.10850273d-022, 0.13593043d-021, &
      0.15802104d-020, 0.17066770d-019, 0.17144103d-018, &
      0.16034986d-017, 0.13978224d-016, 0.11367936d-015, &
      0.86328190d-015, 0.61268884d-014, 0.40672353d-013, &
      0.25273658d-012, 0.14711824d-011, 0.80278299d-011, &
      0.41091402d-010, 0.19742302d-009, 0.89083486d-009, &
      0.37774085d-008, 0.15059816d-007, 0.56479809d-007, &
      0.19935237d-006, 0.66252036d-006, 0.20740078d-005, &
      0.61182560d-005, 0.17014265d-004, 0.44618909d-004, &
      0.11037912d-003, 0.25766142d-003, 0.56771359d-003, &
      0.11809739d-002, 0.23199959d-002, 0.43049288d-002, &
      0.75468199d-002, 0.12501445d-001, 0.19571552d-001, &
      0.28961544d-001, 0.40513891d-001, 0.53581948d-001, &
      0.67004500d-001, 0.79230279d-001, 0.88593773d-001, &
      0.93681831d-001, &
!---- 192 mass points
      0.77409655d-156, 0.16086789d-147, 0.80629006d-141, &
      0.54093564d-135, 0.94860755d-130, 0.61709034d-125, &
      0.18416559d-120, 0.29044729d-116, 0.26756957d-112, &
      0.15509562d-108, 0.59891054d-105, 0.16118296d-101, &
      0.31354502d-098, 0.45428573d-095, 0.50268011d-092, &
      0.43391577d-089, 0.29755374d-086, 0.16466501d-083, &
      0.74554305d-081, 0.27952470d-078, 0.87716976d-076, &
      0.23259299d-073, 0.52560790d-071, 0.10200301d-068, &
      0.17118467d-066, 0.25000890d-064, 0.31958090d-062, &
      0.35943827d-060, 0.35742220d-058, 0.31563226d-056, &
      0.24854461d-054, 0.17518542d-052, 0.11091503d-050, &
      0.63285224d-049, 0.32640513d-047, 0.15261224d-045, &
      0.64856605d-044, 0.25114979d-042, 0.88825877d-041, &
      0.28755901d-039, 0.85386934d-038, 0.23301122d-036, &
      0.58543551d-035, 0.13565924d-033, 0.29040070d-032, &
      0.57516891d-031, 0.10555488d-029, 0.17974166d-028, &
      0.28436612d-027, 0.41851169d-026, 0.57365744d-025, &
      0.73316532d-024, 0.87462172d-023, 0.97487682d-022, &
      0.10162744d-020, 0.99175502d-020, 0.90679935d-019, &
      0.77748937d-018, 0.62560481d-017, 0.47277839d-016, &
      0.33580031d-015, 0.22432090d-014, 0.14102909d-013, &
      0.83496773d-013, 0.46581292d-012, 0.24500828d-011, &
      0.12156600d-010, 0.56928417d-010, 0.25173539d-009, &
      0.10516223d-008, 0.41520999d-008, 0.15500677d-007, &
      0.54737067d-007, 0.18290449d-006, 0.57854152d-006, &
      0.17328415d-005, 0.49162719d-005, 0.13215888d-004, &
      0.33671551d-004, 0.81330137d-004, 0.18628139d-003, &
      0.40468566d-003, 0.83404651d-003, 0.16310742d-002, &
      0.30272498d-002, 0.53331943d-002, 0.89198625d-002, &
      0.14165148d-001, 0.21361440d-001, 0.30593872d-001, &
      0.41617383d-001, 0.53775848d-001, 0.66008778d-001, &
      0.76973597d-001, 0.85275705d-001, 0.89756132d-001, &
!---- 208 mass points
      0.20235148d-169, 0.74059447d-161, 0.59536897d-154, &
      0.61072491d-148, 0.15883652d-142, 0.14996880d-137, &
      0.63912127d-133, 0.14210691d-128, 0.18267116d-124, &
      0.14649079d-120, 0.77702698d-117, 0.28550106d-113, &
      0.75427079d-110, 0.14774910d-106, 0.22016041d-103, &
      0.25503642d-100, 0.23398696d-097, 0.17278061d-094, &
      0.10413802d-091, 0.51867582d-089, 0.21582515d-086, &
      0.75762894d-084, 0.22633311d-081, 0.57994631d-079, &
      0.12836925d-076, 0.24704163d-074, 0.41578644d-072, &
      0.61531774d-070, 0.80464363d-068, 0.93402534d-066, &
      0.96647139d-064, 0.89491807d-062, 0.74423208d-060, &
      0.55772673d-058, 0.37781604d-056, 0.23203477d-054, &
      0.12954644d-052, 0.65918754d-051, 0.30644078d-049, &
      0.13044233d-047, 0.50950413d-046, 0.18298008d-044, &
      0.60535143d-043, 0.18481384d-041, 0.52157648d-040, &
      0.13628642d-038, 0.33021469d-037, 0.74297430d-036, &
      0.15544529d-034, 0.30281051d-033, 0.54990736d-032, &
      0.93205809d-031, 0.14761027d-029, 0.21866167d-028, &
      0.30328575d-027, 0.39425299d-026, 0.48077601d-025, &
      0.55047591d-024, 0.59228031d-023, 0.59932047d-022, &
      0.57077853d-021, 0.51200237d-020, 0.43288970d-019, &
      0.34520352d-018, 0.25980343d-017, 0.18465160d-016, &
      0.12400953d-015, 0.78739909d-015, 0.47294007d-014, &
      0.26885115d-013, 0.14471906d-012, 0.73799202d-012, &
      0.35668483d-011, 0.16345979d-010, 0.71057249d-010, &
      0.29312088d-009, 0.11478574d-008, 0.42686080d-008, &
      0.15079563d-007, 0.50621536d-007, 0.16153214d-006, &
      0.49010140d-006, 0.14142825d-005, 0.38826037d-005, &
      0.10142709d-004, 0.25219137d-004, 0.59696275d-004, &
      0.13455341d-003, 0.28883958d-003, 0.59062503d-003, &
      0.11506237d-002, 0.21359396d-002, 0.37786923d-002, &
      0.63716003d-002, 0.10241474d-001, 0.15693900d-001, &
      0.22929594d-001, 0.31944542d-001, 0.42438953d-001, &
      0.53768573d-001, 0.64969810d-001, 0.74874269d-001, &
      0.82301074d-001, 0.86285859d-001 /
      data qprb(1014:1245) / &
!---- 224 mass points
      0.51008040d-183, 0.31929603d-174, 0.40157722d-167, &
      0.61568714d-161, 0.23245092d-155, 0.31207823d-150, &
      0.18618368d-145, 0.57243634d-141, 0.10074065d-136, &
      0.10969314d-132, 0.78454241d-129, 0.38638386d-125, &
      0.13612701d-121, 0.35400781d-118, 0.69759721d-115, &
      0.10650135d-111, 0.12838532d-108, 0.12422888d-105, &
      0.97882280d-103, 0.63596857d-100, 0.34456055d-097, &
      0.15722190d-094, 0.60960369d-092, 0.20246656d-089, &
      0.58020537d-087, 0.14440946d-084, 0.31405477d-082, &
      0.60006704d-080, 0.10124459d-077, 0.15154483d-075, &
      0.20210289d-073, 0.24109748d-071, 0.25822841d-069, &
      0.24917080d-067, 0.21730092d-065, 0.17178888d-063, &
      0.12345594d-061, 0.80863157d-060, 0.48392715d-058, &
      0.26521916d-056, 0.13340493d-054, 0.61712544d-053, &
      0.26305850d-051, 0.10351546d-049, 0.37669176d-048, &
      0.12697244d-046, 0.39705694d-045, 0.11536110d-043, &
      0.31184624d-042, 0.78537571d-041, 0.18451154d-039, &
      0.40486084d-038, 0.83066623d-037, 0.15953835d-035, &
      0.28713012d-034, 0.48473499d-033, 0.76835046d-032, &
      0.11445711d-030, 0.16037488d-029, 0.21154627d-028, &
      0.26290542d-027, 0.30807309d-026, 0.34063471d-025, &
      0.35564151d-024, 0.35084811d-023, 0.32725900d-022, &
      0.28880244d-021, 0.24127196d-020, 0.19092332d-019, &
      0.14318467d-018, 0.10182373d-017, 0.68696826d-017, &
      0.43991716d-016, 0.26751889d-015, 0.15455455d-014, &
      0.84867140d-014, 0.44310502d-013, 0.22006688d-012, &
      0.10400353d-011, 0.46789068d-011, 0.20044472d-010, &
      0.81798112d-010, 0.31807435d-009, 0.11789164d-008, &
      0.41661293d-008, 0.14041018d-007, 0.45143554d-007, &
      0.13849498d-006, 0.40552535d-006, 0.11335669d-005, &
      0.30256312d-005, 0.77128360d-005, 0.18781331d-004, &
      0.43695130d-004, 0.97142601d-004, 0.20640853d-003, &
      0.41923120d-003, 0.81404811d-003, 0.15113857d-002, &
      0.26833878d-002, 0.45564268d-002, 0.74002113d-002, &
      0.11497000d-001, 0.17087699d-001, 0.24298271d-001, &
      0.33059023d-001, 0.43038178d-001, 0.53615391d-001, &
      0.63916875d-001, 0.72920059d-001, 0.79615038d-001, &
      0.83189266d-001, &
!---- 240 mass points
      0.12451852d-196, 0.12991780d-187, 0.25009722d-180, &
      0.56178121d-174, 0.30217521d-168, 0.56660461d-163, &
      0.46506842d-158, 0.19440620d-153, 0.46069201d-149, &
      0.67009007d-145, 0.63589134d-141, 0.41312954d-137, &
      0.19104674d-133, 0.64929301d-130, 0.16657125d-126, &
      0.32994893d-123, 0.51451459d-120, 0.64229212d-117, &
      0.65133263d-114, 0.54349187d-111, 0.37743984d-108, &
      0.22037980d-105, 0.10917204d-102, 0.46261523d-100, &
      0.16893180d-097, 0.53518709d-095, 0.14800116d-092, &
      0.35927461d-090, 0.76953267d-088, 0.14612588d-085, &
      0.24707546d-083, 0.37350508d-081, 0.50671462d-079, &
      0.61908751d-077, 0.68340957d-075, 0.68371620d-073, &
      0.62169445d-071, 0.51516643d-069, 0.39001376d-067, &
      0.27039686d-065, 0.17206030d-063, 0.10070016d-061, &
      0.54314097d-060, 0.27048531d-058, 0.12459364d-056, &
      0.53174244d-055, 0.21059828d-053, 0.77520450d-052, &
      0.26559092d-050, 0.84809153d-049, 0.25273973d-047, &
      0.70379746d-046, 0.18335063d-044, 0.44737456d-043, &
      0.10234955d-041, 0.21977443d-040, 0.44338003d-039, &
      0.84119266d-038, 0.15022089d-036, 0.25273142d-035, &
      0.40090832d-034, 0.60011673d-033, 0.84833043d-032, &
      0.11333221d-030, 0.14318852d-029, 0.17120835d-028, &
      0.19385961d-027, 0.20800227d-026, 0.21160660d-025, &
      0.20423150d-024, 0.18710737d-023, 0.16280463d-022, &
      0.13460909d-021, 0.10581067d-020, 0.79111460d-020, &
      0.56286436d-019, 0.38125411d-018, 0.24595548d-017, &
      0.15118454d-016, 0.88580504d-016, 0.49489526d-015, &
      0.26374910d-014, 0.13412920d-013, 0.65111328d-013, &
      0.30180832d-012, 0.13362333d-011, 0.56524868d-011, &
      0.22852171d-010, 0.88321319d-010, 0.32641370d-009, &
      0.11538406d-008, 0.39021425d-008, 0.12628187d-007, &
      0.39116122d-007, 0.11599508d-006, 0.32936711d-006, &
      0.89569790d-006, 0.23332613d-005, 0.58231982d-005, &
      0.13926059d-004, 0.31917707d-004, 0.70119137d-004, &
      0.14767395d-003, 0.29818868d-003, 0.57736747d-003, &
      0.10721064d-002, 0.19093920d-002, 0.32618776d-002, &
      0.53455958d-002, 0.84046124d-002, 0.12678456d-001, &
      0.18351561d-001, 0.25489731d-001, 0.33975636d-001, &
      0.43461172d-001, 0.53356223d-001, 0.62868596d-001, &
      0.71098371d-001, 0.77174148d-001, 0.80403773d-001 /
      data qprb(1246:1373) / &
!---- 256 mass points
      0.29540146d-210, 0.50208178d-201, 0.14509335d-193, &
      0.46913583d-187, 0.35359397d-181, 0.91144574d-176, &
      0.10136209d-170, 0.56754753d-166, 0.17847925d-161, &
      0.34184302d-157, 0.42436752d-153, 0.35864591d-149, &
      0.21469376d-145, 0.94051577d-142, 0.30983962d-138, &
      0.78550152d-135, 0.15630453d-131, 0.24832856d-128, &
      0.31972925d-125, 0.33800826d-122, 0.29682323d-119, &
      0.21876645d-116, 0.13658255d-113, 0.72838408d-111, &
      0.33431061d-108, 0.13296528d-105, 0.46114419d-103, &
      0.14025800d-100, 0.37608844d-098, 0.89335232d-096, &
      0.18882727d-093, 0.35662449d-091, 0.60412782d-089, &
      0.92123102d-087, 0.12687510d-084, 0.15830825d-082, &
      0.17947950d-080, 0.18539433d-078, 0.17492932d-076, &
      0.15113381d-074, 0.11983435d-072, 0.87388393d-071, &
      0.58729939d-069, 0.36444489d-067, 0.20919910d-065, &
      0.11127348d-063, 0.54933568d-062, 0.25210019d-060, &
      0.10770598d-058, 0.42899179d-057, 0.15950824d-055, &
      0.55436741d-054, 0.18031128d-052, 0.54949673d-051, &
      0.15707548d-049, 0.42161469d-048, 0.10637272d-046, &
      0.25250874d-045, 0.56449599d-044, 0.11895249d-042, &
      0.23647681d-041, 0.44387936d-040, 0.78731310d-039, &
      0.13205880d-037, 0.20962490d-036, 0.31512301d-035, &
      0.44892461d-034, 0.60646536d-033, 0.77741057d-032, &
      0.94616999d-031, 0.10939939d-029, 0.12023504d-028, &
      0.12567594d-027, 0.12499833d-026, 0.11836003d-025, &
      0.10674923d-024, 0.91745691d-024, 0.75173209d-023, &
      0.58746940d-022, 0.43806011d-021, 0.31180536d-020, &
      0.21193550d-019, 0.13761221d-018, 0.85388712d-018, &
      0.50650722d-017, 0.28731517d-016, 0.15590500d-015, &
      0.80951675d-015, 0.40233406d-014, 0.19145601d-013, &
      0.87255717d-013, 0.38095919d-012, 0.15938055d-011, &
      0.63910596d-011, 0.24569446d-010, 0.90573898d-010, &
      0.32025235d-009, 0.10863148d-008, 0.35357535d-008, &
      0.11044763d-007, 0.33117736d-007, 0.95339500d-007, &
      0.26355296d-006, 0.69970864d-006, 0.17843872d-005, &
      0.43716822d-005, 0.10291000d-004, 0.23279572d-004, &
      0.50612366d-004, 0.10576827d-003, 0.21248156d-003, &
      0.41039425d-003, 0.76214868d-003, 0.13610624d-002, &
      0.23375178d-002, 0.38610584d-002, 0.61343092d-002, &
      0.93748300d-002, 0.13782486d-001, 0.19493216d-001, &
      0.26524913d-001, 0.34726352d-001, 0.43743821d-001, &
      0.53020240d-001, 0.61836747d-001, 0.69397142d-001, &
      0.74943513d-001, 0.77880553d-001/
!-----------------------------------------------------------------------
      if (nmes > maxy) then
          call wrtlin('    *** ERROR *** TOO MANY OBSERVATIONS')
          go to 500
      end if
!
      if ((offlag(1) .or. offlag(2) .or. offlag(3)) .and. ilfit == 1) &
      then
          call newlin
          call wrtlin('    *** WARNING *** To fit a conventional '// &
                      'regression model with offset(s), use')
          call wrtlin('    the FIT command with single mass point(s) '// &
                      'and no endpoints.')
          return
      end if
!
      iend(1) = 0
      iend(2) = 0
!
!---- left endpoint in model
      if (endind == 'b' .or. endind == 'l') then
          iend(1) = 1
      end if
!
!---- right endpoint in model
      if (endind == 'b' .or. endind == 'r') then
          iend(2) = 1
      end if
!
      nend = iend(1) + iend(2)
!
!---- check validity of model specification
!-----***********
      call fitchk(mnames,iyvar,icvar,line,ni,nilev,ilev,name,nvar,yname, &
                  cname,idel,ilfit,offnme,n1var,n2var,n1lev,n2lev,order, &
                  iname,univar)
!-----***********
!
!---- error in model specification
      if (idel(1) == -1) then
          return
      end if
!
!-----***********
      call modchk(endind,offlag,link,corr,order,family,ifail, &
                  univar,cname,yname,rname,offnme,nlevel,irvar,ilfit, &
                  bivar,eqscale,cquad)
!-----***********
!
      if (ifail) then
          return
      end if
!
      if (eqscale .and. .not. corr) then
          call wrtlin('    *** ERROR *** Constrained random effects '// &
                      'models must be correlated')
          ifail = .true.
          return
      end if
!
      if ((bivar .or. depend) .and. n1var == 0) then
          n1var = ni/2
          n1lev = nilev/2
      else if (trivar .and. n1var == 0 .and. n2var == 0) then
          n1var = ni/3
          n2var = n1var
          n1lev = nilev/3
          n2lev = n1lev
      end if
!
!---- get position of y-variate in X-vector
!-----***********
      call fndpos(nvar,yname,name,iypos,ipos)
!-----***********
!
!---- get position of risk-variate in X-vector
!-----***********
      call fndpos(nvar,rname,name,irpos,ipos)
!-----***********
!
      if (.not. univar .and. irpos == 0) then
          call wrtlin( &
          '    *** ERROR *** multivariate model with univariate data')
          ifail = .true.
          return
      end if
!
      rmax = 1
      r1 = 0
      ymin(1) = 1d8
      ymax(1) = 0
      ymin(2) = 1d8
      ymax(2) = 0
      ymin(3) = 1d8
      ymax(3) = 0
!
      if (depend) then
          family(2) = family(1)
          link(2) = link(1)
      end if
!
      do 1 i = 1,nmes
!
!-------- store the risk-variate values in integer form
          if (univar .and. .not. depend) then
              risk(i) = 1
          else
              risk(i) = x(i,irpos)
          end if
!
          if (risk(i) > rmax) then
              rmax = risk(i)
          end if
!
          if (risk(i) == 1) then
              r1 = r1+1
          end if
!
!-------- invalid risk-variate values
          if ((bivar .or. depend) .and. risk(i) /= 1 .and. &
          risk(i) /= 2) then
              call wrtlin( &
              '    *** ERROR *** INVALID RISK VALUES')
              go to 500
          end if
!
!-------- invalid risk-variate values
          if (trivar .and. risk(i) /= 1 .and. risk(i) /= 2 .and. &
          risk(i) /= 3) then
              call wrtlin( &
              '    *** ERROR *** INVALID trivariate risk VALUES')
              go to 500
          end if
!
!-------- response value too large for integer storage
          if (family(risk(i)) == 'p' .and. x(i,iypos) > iz) then
              y(i) = -1
              go to 2
          end if
!
!-------- store the y-variate values
          y(i) = x(i,iypos)
!
!-------- invalid y-variate values
    2     if (ilfit /= 3) then
!
              if ((family(risk(i)) == 'b' .and. .not. order .and. &
              y(i) /= 0 .and. y(i) /= 1) .or. &
              (family(risk(i)) == 'p' .and. (y(i) < 0 .or. &
              y(i) /= int(y(i)))) .or. (order .and. y(i) /= int(y(i)))) &
              then
                  call wrtlin( &
                  '    *** ERROR *** INVALID Y-VARIATE VALUES')
                  go to 500
              end if
!
          end if
!
          if (y(i) < ymin(risk(i))) then
              ymin(risk(i)) = y(i)
          end if
!
          if (y(i) > ymax(risk(i))) then
              ymax(risk(i)) = y(i)
          end if
!
    1 end do
!
      if (bivar .and. rmax <= 1) then
          call wrtlin( &
          '    *** ERROR *** bivariate model with non-bivariate data')
          ifail = .true.
          return
      end if
!
      if (trivar .and. rmax <= 2) then
          call wrtlin( &
          '    *** ERROR *** trivariate model with non-trivariate data')
          ifail = .true.
          return
      end if
!
      nnorm = 0
!
      if (family(1) == 'g') then
          nnorm = 1
      end if
!
      if (.not. univar .and. family(2) == 'g') then
          nnorm = nnorm+1
      end if
!
      if (trivar .and. family(3) == 'g') then
          nnorm = nnorm+1
      end if
!
      ncat(1) = 1
      ncat(2) = 1
      ncat(3) = 1
!
      if (order) then
!
          do 10 i = 1,nmes
!
              if (ymin(risk(i)) == 0 .and. ymax(risk(i)) == 1) then
                  y(i) = y(i) + 1
              end if
!
              if (y(i) > ncat(risk(i))) then
                  ncat(risk(i)) = y(i)
              end if
!
   10     end do
!
          if (depend) then
              ncat(2) = ncat(1)
          end if
!
      end if
!
      maxcat = ncat(1)
!
      if (ncat(2) > maxcat) then
          maxcat = ncat(2)
      end if
!
      if (ncat(3) > maxcat) then
          maxcat = ncat(3)
      end if
!
      ncut = ncat(1) + ncat(2) + ncat(3) - 3
      ncorr = 0
!
      if (ilfit /= 0) then
          nsca = 0
      else if (univar .and. .not. depend .and. nlevel == 1) then
          nsca = 1
      else if (univar) then
          nsca = 2
      else if (bivar .and. .not. corr .and. eqscale) then
          nsca = 1
      else if (bivar .and. .not. corr) then
          nsca = 2
      else if (bivar .and. eqscale) then
          nsca = 1
          ncorr = 1
      else if (bivar) then
          nsca = 2
          ncorr = 1
      else if (trivar .and. .not. corr) then
          nsca = 3
      else if (trivar) then
          nsca = 3
          ncorr = 3
      end if
!
      nre = nsca+ncorr
      npar = nilev + ncut + nnorm + nre + nend
!
      if (order .and. bivar .and. ilfit == 2) then
          npar = npar+3
      else if (order .and. trivar .and. ilfit == 2) then
          npar = npar+6
      end if
!
      if (npar == 0) then
          call wrtlin('    *** ERROR *** '// &
                      'THE NULL MODEL CANNOT BE FITTED IN THIS WAY')
          call wrtlin('    Fit a model with '// &
                      'a single variable which is identically zero')
          go to 500
      else if (npar > maxpar) then
          call wrtlin('    *** ERROR *** TOO MANY PARAMETERS')
          write (outbuf,'(a,i3,a,i10)') '    The maximum is ',maxpar, &
          ' - you want to estimate ',npar
          call wrtlin(outbuf)
          go to 500
      end if
!
      if (qflag) then
!
!---------**********
          call clock('SHUFFLE',clock1,1)
!---------**********
!
      end if
!
!---- move model variates to front of X-vector (in order of LFIT or FIT)
!---- i.e. you have to do it backwards or something
      do 6 j = ni,1,-1
!
          do 5 i = 1,nvar
!
              if (name(i) == mnames(j)) then
!
!---------------- get position of variable in X-vector
!-----------------***********
                  call fndpos(nvar,name(i),name,ixpos,ipos)
!-----------------***********
!
!---------------- get number of levels in variable
!-----------------***********
                  call fndlev(nvar,name(i),name,ixlev,ilev)
!-----------------***********
!
!---------------- move variable to front of X-vector
!-----------------***********
                  call shuffl(x,ixpos,1,ncol,nmes,ipos,maxcol,ixlev, &
                              nvar)
!-----------------***********
!
                  go to 6
              end if
!
    5     end do
!
    6 end do
!
      if (qflag) then
!
!---------**********
          call clock('',clock1,2)
!---------**********
!
      end if
!
      ilevt = 0
!
!---- move intercept term to front of list of model variables
      do 1120 i = 1,ni
!
!-------- get number of levels in model variate
!---------***********
          call fndlev(nvar,mnames(i),name,ixlev,ilev)
!---------***********
!
          ilevt = ilevt+ixlev
!
          if (i >= 2 .and. mnames(i) == iname(1)) then
!
              do 1110 j = i,2,-1
                  mnames(j) = mnames(j-1)
 1110         end do
!
              mnames(1) = iname(1)
!
!------------ get position of intercept in X-vector
!-------------***********
              call fndpos(nvar,iname(1),name,ixpos,ipos)
!-------------***********
!
!------------ move intercept to first position in X-vector
!-------------***********
              call shuffl(x,ixpos,1,ncol,nmes,ipos,maxcol,ixlev,nvar)
!-------------***********
!
          end if
!
          if (.not. univar .or. depend) then
!
              if (i >= n1var+2 .and. mnames(i) == iname(2)) then
!
                  do 2110 j = i,n1var+2,-1
                      mnames(j) = mnames(j-1)
 2110             end do
!
                  mnames(n1var+1) = iname(2)
!
!---------------- get position of intercept in X-vector
!-----------------***********
                  call fndpos(nvar,iname(2),name,ixpos,ipos)
!-----------------***********
!
!---------------- move intercept to first position in X-vector
!-----------------***********
                  call shuffl(x,ixpos,n1lev+1,ncol,nmes,ipos,maxcol, &
                              ixlev,nvar)
!-----------------***********
!
              end if
!
              if (trivar) then
!
                  if (i >= n1var+n2var+2 .and. &
                  mnames(i) == iname(3)) then
!
                      do 3110 j = i,n1var+n2var+2,-1
                          mnames(j) = mnames(j-1)
 3110                 end do
!
                      mnames(n1var+n2var+1) = iname(3)
!
!-------------------- get position of intercept in X-vector
!---------------------***********
                      call fndpos(nvar,iname(3),name,ixpos,ipos)
!---------------------***********
!
!-------------------- move intercept to first position in X-vector
!---------------------***********
                      call shuffl(x,ixpos,n1lev+n2lev+1,ncol,nmes,ipos, &
                                  maxcol,ixlev,nvar)
!---------------------***********
!
                  end if
!
              end if
!
          end if
!
 1120 end do
!
!---- initialise parameter estimates to zero
      do 3 i = 1,npar
          beta(i) = 0
    3 end do
!
!---- initialise the variance-covariance matrix
      do 122 i = 1,npar*(npar+1)/2
          cov(i) = 0
  122 end do
!
!---- initialise ALIAS vector to zero
      do 125 i = 1,2*npar
          alias(i) = 0
  125 end do
!
!---- FACYES detects if a factor has been found
      facyes = .false.
!---- GMYES detects if an intercept has been found
      gmyes = .false.
      k = 1
!
!---- intrinsically alias appropriate model parameters
      do 130 i = 1,ni
!
!-------- get number of levels in variable
!---------***********
          call fndlev(nvar,mnames(i),name,ixlev,ilev)
!---------***********
!
          if (ixlev > 1) then
!
!------------ factor or constant term already found, so intrinsically
!------------ alias first level of new factor
              if (facyes .or. gmyes .or. order) then
                  alias(k+ixlev-1) = -1
              end if
!
              facyes = .true.
          end if
!
          if (mnames(i) == iname(1)) then
              gmyes = .true.
          end if
!
!-------- update number of levels found
          k = k+ixlev
  130 end do
!
!---- move intrinsically aliased columns temporarily to the end of X
!---- start at end of X matrix and work to front so that intrinsic alias
!---- vector points to correct columns throughout the loop
      ninali = 0
      ninal1 = 0
      ninal2 = 0
!
      do 135 i = nilev,1,-1
!
          if (idnint(alias(i)) == -1) then
!------------ update number of intrinsically aliased variables/factor
!             levels
              ninali = ninali+1
!
              if (.not. univar .and. i <= n1lev) then
                  ninal1 = ninal1+1
              else if (trivar .and. i <= n1lev+n2lev) then
                  ninal2 = ninal2+1
              end if
!
!------------ move intrinsically aliased variables/factor levels to end 
!------------ of X-vector
!-------------***********
              call shuffl(x,i,ncol,ncol,nmes,ipos,maxcol,1,nvar)
!-------------***********
!
          end if
!
  135 end do
!
!---- update total number of levels of all parameters in model
      nilev = nilev-ninali
!
      if (.not. univar) then
          n1lev = n1lev-ninal1
      end if
!
      if (trivar) then
          n2lev = n2lev-ninal2
      end if
!
!---- update number of parameters in model
      nest = npar-ninali
!
      if (offlag(1)) then
!
!-------- get position of offset variable in X-vector
!---------***********
          call fndpos(nvar,offnme(1),name,offpos(1),ipos)
!---------***********
!
      end if
!
      if (.not. univar .and. offlag(2)) then
!
!-------- get position of offset variable in X-vector
!---------***********
          call fndpos(nvar,offnme(2),name,offpos(2),ipos)
!---------***********
!
      end if
!
      if (trivar .and. offlag(3)) then
!
!-------- get position of offset variable in X-vector
!---------***********
          call fndpos(nvar,offnme(3),name,offpos(3),ipos)
!---------***********
!
      end if
!
      ifail = .true.
!
      if (nilev*(nilev+2) > nmes*(maxcol-ncol-2)) then
          call wrtlin('    *** ERROR *** NOT ENOUGH WORKSPACE')
          go to 500
      end if
!
!     If Fixed Effects model then workspace within subroutine FEFIT is
!     used and a call to LFIT is not required
!
      if (ilfit /= 3 .and. .not. inflag .and. .not. order .and. &
      .not. offlag(1) .and. .not. offlag(2) .and. .not. offlag(3)) then
!
!-------- fit logistic/log-linear regression model to obtain starting
!-------- values
!---------*********
          call lfit(x,work,nilev,nmes,maxcol,y,beta,alias(maxpar+1),cov, &
                    nilev*(nilev+1)/2,ilfit,xll,con,tol,niter,ifail, &
                    arith,link,bivar,risk,robust,sig_e,n1lev,family, &
                    trivar,n2lev,nnorm,univar,depend)
!---------*********
!
!-------- don't do the rest if logistic/log-linear regression only
!-------- (but calculate degrees of freedom first)
!-------- or if LFIT fails for any reason
          if (ifail) then
              go to 500
          end if
!
      end if
!
!---- logistic/log-linear regression model
      if (ilfit == 1) then
          go to 777
      end if
!
!---- get position of case variate in X-vector
!-----***********
      call fndpos(nvar,cname(1),name,icpos(1),ipos)
!-----***********
!
      if (nlevel == 2) then
!
!-------- get position of case variate in X-vector
!---------***********
          call fndpos(nvar,cname(2),name,icpos(2),ipos)
!---------***********
!
      end if
!
!---- get number of cases and number of measurements per case
      k = 1
      nsub(nlevel) = 0
      n1sub = 0
      maxit = 0
!
      do 144 i = 2,nmes
!
!-------- measurement refers to new case
          if (x(i,icpos(nlevel)) > x(i-1,icpos(nlevel))) then
!------------ update number of cases
              nsub(nlevel) = nsub(nlevel) + 1
!
!------------ too many cases
              if (nmes + nsub(nlevel) > 4*maxy) then
                  call wrtlin('    *** ERROR *** TOO MANY CASES')
                  go to 500
              end if
!
!------------ store number of measurements for current case
              y(nmes + nsub(nlevel)) = k
              it(nlevel,nsub(nlevel)) = k
!
              if (k > maxit) then
                  maxit = k
              end if
!
              if (it(nlevel,nsub(nlevel)) == 1) then
                  n1sub = n1sub + 1
              end if
!
!------------ initialise number of measurements for new case
              k = 1
!-------- measurement refers to same case as before
          else if (nlevel == 1 .and. &
          x(i,icpos(1)) == x(i-1,icpos(1))) then
!------------ update number of measurements for current case
              k = k+1
!-------- measurement refers to same case as before
          else if (nlevel == 2 .and. &
          x(i,icpos(2)) == x(i-1,icpos(2)) .and. &
          x(i,icpos(1)) > x(i-1,icpos(1))) then
!------------ update number of measurements for current case
              k = k+1
          else if (nlevel == 1 .or. (nlevel == 2 .and. &
          (x(i,icpos(2)) < x(i-1,icpos(2)) .or. &
          x(i,icpos(1)) < x(i-1,icpos(1))))) then
              call wrtlin('    *** ERROR *** '// &
                          'INVALID NON-ASCENDING CASE STRUCTURE')
              return
          end if
!
  144 end do
!
!---- update number of cases
      nsub(nlevel) = nsub(nlevel) + 1
!
!---- too many cases
      if (nmes + 3*nsub(nlevel) > 4*maxy) then
          call wrtlin('    *** ERROR *** TOO MANY CASES')
          go to 500
      end if
!
!---- store number of measurements for last case
      y(nmes + nsub(nlevel)) = k
      it(nlevel,nsub(nlevel)) = k
!
      if (k > maxit) then
          maxit = k
      end if
!
      if (it(nlevel,nsub(nlevel)) == 1) then
          n1sub = n1sub + 1
      end if
!
      if (nlevel == 2) then
!-------- get number of cases and number of measurements per case
          k = 1
          nsub(1) = 0
!
          do 1444 i = 2,nmes
!
!------------ measurement refers to same case as before
              if (x(i,icpos(1)) == x(i-1,icpos(1))) then
!---------------- update number of measurements for current case
                  k = k+1
!------------ measurement refers to new case
              else
!
                  do 1325 j = 1,nsub(1)
!
                      if (x(i-1,icpos(1)) == cind(j)) then
                          it(1,j) = it(1,j) + k
                          k = 1
                          go to 1444
                      end if
!
 1325             end do
!
!---------------- update number of cases
                  nsub(1) = nsub(1) + 1
                  cind(nsub(1)) = x(i-1,icpos(1))
!
!---------------- too many cases
                  if (nmes + nsub(1) > 4*maxy) then
                      call wrtlin( &
                      '    *** ERROR *** TOO MANY LEVEL 2 CASES')
                      go to 500
                  end if
!
!---------------- store number of measurements for current case
                  it(1,nsub(1)) = k
!---------------- initialise number of measurements for new case
                  k = 1
              end if
!
 1444     end do
!
!-------- update number of cases
          nsub(1) = nsub(1) + 1
!
!-------- too many cases
          if (nmes + 3*nsub(1) > 4*maxy) then
              call wrtlin('    *** ERROR *** TOO MANY LEVEL 2 CASES')
              go to 500
          end if
!
!-------- store number of measurements for last case
          it(1,nsub(1)) = k
      end if
!
!     If Fixed effect model then go straight to solution routine
!
      if (ilfit == 3) then
          go to 777
      end if
!
!---- Workspace for fitting procedure starts at X(1,NCOL+1) and needs
!---- NEST*(NEST+2) cells for matrix inversion. Are there sufficient
!---- cells available?
      if (nest*(nest+2) > mxx - nmes*ncol) then
          call wrtlin('    *** ERROR *** NOT ENOUGH WORKSPACE')
          go to 500
      end if
!
      if (univar .and. family(1) /= 'p' .and. nsub(1) == nmes .and. &
      ilfit == 0) then
          call wrtlin('    *** ERROR *** One observation per case')
          call wrtlin('    Random effects model cannot be fitted')
          ifail = .true.
          return
      else if (bivar .and. family(1) /= 'p' .and. nsub(1) == r1 .and. &
      .not. eqscale) then
          call wrtlin('    *** ERROR *** Initial conditions')
          ifail = .true.
          return
      end if
!
      do 149 j = 1,3
!
          do 148 i = 1,nm(j)
!
              if (nm(j) == 1) then
!---------------- mass point location
                  quad(j,1,1) = qloc(1)
!---------------- mass point probability
                  quad(j,2,1) = qprb(1)
              else if (nm(j) <= 16 .and. i <= nm(j)/2) then
!---------------- positive mass point locations
                  quad(j,1,i) = qloc(nm(j)*(nm(j) - 2)/8 + 1 + i)
!---------------- associated mass point probabilities
                  quad(j,2,i) = qprb(nm(j)*(nm(j) - 2)/8 + 1 + i)
              else if (nm(j) <= 48 .and. i <= nm(j)/2) then
!---------------- positive mass point locations
                  quad(j,1,i) = qloc(nm(j)*(nm(j) - 4)/16 + 17 + i)
!---------------- associated mass point probabilities
                  quad(j,2,i) = qprb(nm(j)*(nm(j) - 4)/16 + 17 + i)
              else if (nm(j) <= 112 .and. i <= nm(j)/2) then
!---------------- positive mass point locations
                  quad(j,1,i) = qloc(nm(j)*(nm(j) - 8)/32 + 89 + i)
!---------------- associated mass point probabilities
                  quad(j,2,i) = qprb(nm(j)*(nm(j) - 8)/32 + 89 + i)
              else if (nm(j) <= 256 .and. i <= nm(j)/2) then
!---------------- positive mass point locations
                  quad(j,1,i) = qloc(nm(j)*(nm(j) - 16)/64 + 285 + i)
!---------------- associated mass point probabilities
                  quad(j,2,i) = qprb(nm(j)*(nm(j) - 16)/64 + 285 + i)
              else
!---------------- negative mass point locations
                  quad(j,1,i) = -quad(j,1,nm(j) + 1 - i)
!---------------- associated mass point probabilities
                  quad(j,2,i) = quad(j,2,nm(j) + 1 - i)
              end if
!
  148     end do
!
  149 end do
!
!---- calculate matrix containing product (1-Yit) and product (Yit) over
!---- t for each i. It's stored in Y starting from position NMES+NSUB+1
      j = 0
!
      do 160 i = 1,nsub(1)
!
!-------- initialise the products; Poisson model has only left endpoint
          if (family(1) /= 'g') then
              y(nmes + nsub(1) + i) = 1
              prod(i) = 1
!
              if (family(1) /= 'p') then
                  y(nmes + 2*nsub(1) + i) = 1
                  prod(nsub(1) + i) = 1
              end if
!
          end if
!
          do 150 k = 1,it(1,i)
!------------ update number of measurements
              j = j+1
!
              if (family(1) /= 'g') then
!
!---------------- update prod{t=1,...,T_i}(max(1,y_it)-y_it) and store
!---------------- in y-vector
                  if (y(j) /= 0) then
                      y(nmes + nsub(1) + i) = 0
                      prod(i) = 0
                  end if
!
!---------------- update prod{t=1,...,T_i}y_it and store in y-vector
                  if (family(1) /= 'p' .and. y(j) == 0) then
                      y(nmes + 2*nsub(1) + i) = 0
                      prod(nsub(1) + i) = 0
                  end if
!
              end if
!
  150     end do
!
!-------- store dummy value of prod{t=1,...,T_i}y_it for Poisson model
          if (family(1) == 'p') then
              y(nmes + 2*nsub(1) + i) = 0
              prod(nsub(1) + i) = 0
          end if
!
  160 end do
!
!---- initialise beta array
      if (inflag) then
!
          do 200 i = 1,nilev
              beta(i) = xinit(i)
  200     end do
!
      end if
!
      if (order) then
!
          if (cutflag) then
!
              do 202 i = nilev+1,nest-nre
                  beta(i) = xcut(i-nilev)
  202         end do
!
          else
              beta(nilev+1) = -(ncat(1) - 1)/2
!
              do 203 i = nilev+2,nest-nre-ncat(2)-ncat(3)+2
                  beta(i) = beta(i-1) + 1
  203         end do
!
              if (.not. univar .or. depend) then
                  beta(nilev + ncat(1)) = -(ncat(2) - 1)/2
!
                  do 204 i = nilev+ncat(1)+1,nest-nre-ncat(3)+1
                      beta(i) = beta(i-1) + 1
  204             end do
!
                  if (trivar) then
                      beta(nilev + ncat(1) + ncat(2) - 1) = &
                      -(ncat(3) - 1)/2
!
                      do 205 i = nilev+ncat(1)+ncat(2),nest-nre
                          beta(i) = beta(i-1) + 1
  205                 end do
!
                  end if
!
              end if
!
          end if
!
      end if
!
!---- initial estimate of scale parameter, omega
!---- this is done i) if no updating operator is present or
!----             ii) if there is updating from a standard
!----                 logistic/log-linear to a
!----                 logistic-normal/log-linear-normal mixture model or
!----            iii) if there is no updating of estimates from the
!----                 previous model fit
      if (univar .and. ilfit == 0) then
!
          if (family(1) == 'g' .and. sigflag(1)) then
              beta(nest-nsca) = sig(1)
          else if (family(1) == 'g') then
              beta(nest-nsca) = sig_e(1)
          end if
!
          if (nlevel == 1 .and. .not. depend) then
              beta(nest - iend(1) - iend(2)) = sca(1)
          else if (nlevel == 1) then
              beta(nest-1) = sca(1)
              beta(nest) = sca(2)
          else
              beta(nest - 1 - iend(1) - iend(2)) = sca(1)
              beta(nest - iend(1) - iend(2)) = sca(2)
          end if
!
      else if (bivar .and. ilfit == 0) then
!
          if (family(1) == 'g' .and. sigflag(1)) then
              beta(nest-nnorm-nsca-ncorr+1) = sig(1)
          else if (family(1) == 'g') then
              beta(nest-nnorm-nsca-ncorr+1) = sig_e(1)
          end if
!
          if (family(2) == 'g' .and. sigflag(2)) then
              beta(nest-nsca-ncorr) = sig(2)
          else if (family(2) == 'g') then
              beta(nest-nsca-ncorr) = sig_e(2)
          end if
!
          if (.not. eqscale) then
              beta(nest-ncorr-1) = sca(1)
              beta(nest-ncorr) = sca(2)
          else
              beta(nest-ncorr) = sca(1)
          end if
!
          if (corr) then
              beta(nest) = rho(1)
          end if
!
      else if (trivar .and. ilfit == 0) then
!
          if (family(1) == 'g' .and. sigflag(1)) then
              beta(nest-nnorm-ncorr-2) = sig(1)
          else if (family(1) == 'g') then
              beta(nest-nnorm-ncorr-2) = sig_e(1)
          end if
!
          if (family(2) == 'g' .and. family(3) == 'g' .and. &
          sigflag(2)) then
              beta(nest-ncorr-4) = sig(2)
          else if (family(2) == 'g' .and. sigflag(2)) then
              beta(nest-ncorr-3) = sig(2)
          else if (family(2) == 'g' .and. family(3) == 'g') then
              beta(nest-ncorr-4) = sig_e(2)
          else if (family(2) == 'g') then
              beta(nest-ncorr-3) = sig_e(2)
          end if
!
          if (family(3) == 'g' .and. sigflag(3)) then
              beta(nest-ncorr-3) = sig(3)
          else if (family(3) == 'g') then
              beta(nest-ncorr-3) = sig_e(3)
          end if
!
          beta(nest-ncorr-2) = sca(1)
          beta(nest-ncorr-1) = sca(2)
          beta(nest-ncorr) = sca(3)
!
          if (corr) then
              beta(nest-2) = rho(1)
              beta(nest-1) = rho(2)
              beta(nest) = rho(3)
          end if
!
      end if
!
!---- initial estimates of endpoints
      if (iend(1) == 1) then
          beta(nest - iend(2)) = est0
      end if
!
      if (iend(2) == 1) then
          beta(nest) = est1
      end if
!
      if (order .and. ilfit == 0 .and. .not. inflag .and. .not. cutflag) &
      then
          call newlin
          call wrtlin('    Initial Homogeneous Fit:')
!
          if (qflag) then
!
!-------------**********
              call clock('LFIT',clock1,1)
!-------------**********
!
          end if
!
!-------- fit homogeneous model
!---------********
          call fit(beta,alias(maxpar+1),beta(maxpar*(maxpar+7)/2+1), &
                   beta(maxpar*(maxpar+9)/2+1),cov,prod,work,y,x, &
                   nest*(nest+1)/2,nsub,quad,nmes,it,nm,nest,nilev,iend, &
                   con,alp,xll,tol,nmeil,niter,arith,endind,ifail, &
                   offlag,offpos,maxcol,link,risk,corr,bivar,robust, &
                   order,ncat,n1lev,family,trivar,nnorm,ncorr,n2lev,2, &
                   univar,nsca,nlevel,depend,eqscale,dfirst,cquad,maxit, &
                   maxcat)
!---------********
!
          if (qflag) then
!
!-------------**********
              call clock('',clock1,2)
!-------------**********
!
          end if
!
          if (ifail) then
              call wrtlin('    *** FITTING ABANDONED')
              ilfit = ilprev
              return
          end if
!
      else if (((offlag(1) .or. offlag(2) .or. offlag(3)) .and. &
      ilfit == 0 .and. .not. inflag) .or. (.not. univar .and. order &
      .and. ilfit == 0 .and. .not. cutflag)) then
          call newlin
          call wrtlin('    Initial Homogeneous Fit:')
!
          onm(1) = 1
          oquad(1,1,1) = 0
          oquad(1,2,1) = 1
          oiend(1) = 0
          oiend(2) = 0
!
          if (nlevel == 2) then
              onm(2) = 1
              oquad(2,1,1) = 0
              oquad(2,2,1) = 1
          end if
!
          if (.not. univar) then
              onm(2) = 1
              oquad(2,1,1) = 0
              oquad(2,2,1) = 1
              ocorr = .false.
          end if
!
          if (trivar) then
              onm(3) = 1
              oquad(3,1,1) = 0
              oquad(3,2,1) = 1
          end if
!
          if (qflag) then
!
!-------------**********
              call clock('LFIT',clock1,1)
!-------------**********
!
          end if
!
!-------- fit homogeneous model
!---------********
          call fit(beta,alias(maxpar+1),beta(maxpar*(maxpar+7)/2+1), &
                   beta(maxpar*(maxpar+9)/2+1),cov,prod,work,y,x, &
                   nest*(nest+1)/2,nsub,oquad,nmes,it,onm,nest,nilev, &
                   oiend,con,alp,xll,tol,nmeil,niter,arith,'n',ifail, &
                   offlag,offpos,maxcol,link,risk,ocorr,bivar,robust, &
                   order,ncat,n1lev,family,trivar,nnorm,ncorr,n2lev,2, &
                   univar,nsca,nlevel,depend,eqscale,dfirst,cquad,maxit, &
                   maxcat)
!---------********
!
          if (qflag) then
!
!-------------**********
              call clock('',clock1,2)
!-------------**********
!
          end if
!
          if (ifail) then
              call wrtlin('    *** FITTING ABANDONED')
              ilfit = ilprev
              return
          end if
!
      end if
!
      if (((offlag(1) .or. offlag(2) .or. offlag(3)) .and. ilfit == 2) &
      .or. (.not. univar .and. order .and. ilfit == 2)) then
          onm(1) = nm(1)
          oquad(1,1,1) = quad(1,1,1)
          oquad(1,2,1) = quad(1,2,1)
          nm(1) = 1
          quad(1,1,1) = 0
          quad(1,2,1) = 1
!
          if (.not. univar .or. nlevel == 2) then
              onm(2) = nm(2)
              oquad(2,1,1) = quad(2,1,1)
              oquad(2,2,1) = quad(2,2,1)
              nm(2) = 1
              quad(2,1,1) = 0
              quad(2,2,1) = 1
!
              if (.not. univar) then
                  ocorr = corr
                  corr = .false.
              end if
!
              if (trivar) then
                  onm(3) = nm(3)
                  oquad(3,1,1) = quad(3,1,1)
                  oquad(3,2,1) = quad(3,2,1)
                  nm(3) = 1
                  quad(3,1,1) = 0
                  quad(3,2,1) = 1
              end if
!
          end if
!
      end if
!
!---- fit mixture model
!-----********
      call fit(beta,alias(maxpar+1),beta(maxpar*(maxpar+7)/2+1), &
               beta(maxpar*(maxpar+9)/2+1),cov,prod,work,y,x, &
               nest*(nest+1)/2,nsub,quad,nmes,it,nm,nest,nilev,iend,con, &
               alp,xll,tol,nmeil,niter,arith,endind,ifail,offlag,offpos, &
               maxcol,link,risk,corr,bivar,robust,order,ncat,n1lev, &
               family,trivar,nnorm,ncorr,n2lev,ilfit,univar,nsca,nlevel, &
               depend,eqscale,dfirst,cquad,maxit,maxcat)
!-----********
!
      if (((offlag(1) .or. offlag(2) .or. offlag(3)) .and. ilfit == 2) &
      .or. (.not. univar .and. order .and. ilfit == 2)) then
          nm(1) = onm(1)
          quad(1,1,1) = oquad(1,1,1)
          quad(1,2,1) = oquad(1,2,1)
!
          if (nlevel == 2) then
              nm(2) = onm(2)
              quad(2,1,1) = oquad(2,1,1)
              quad(2,2,1) = oquad(2,2,1)
          end if
!
          if (.not. univar) then
              nm(2) = onm(2)
              quad(2,1,1) = oquad(2,1,1)
              quad(2,2,1) = oquad(2,2,1)
              corr = ocorr
          end if
!
          if (trivar) then
              nm(3) = onm(3)
              quad(3,1,1) = oquad(3,1,1)
              quad(3,2,1) = oquad(3,2,1)
          end if
!
      end if
!
      if (.not. ifail) then
          go to 777
      end if
!
  500 continue
!
      call wrtlin('    *** FITTING ABANDONED')
      ilfit = ilprev
      return
!
  777 continue
!
      if (ilfit == 3) then
!
!         Fixed Effects model
!
!---------**********

!          NB  -  removed fefit functionality due to copyright constraints   !!!!!!!
!          call fefit(beta,cov,y,x,ncov,nsub,nmes,it,nest,nilev,con,tol, &
!                     arith,n1sub,ifail,maxcol,ndum,gamma,gammse,sig_e)
!---------**********
!
!         Write the dummy variables & their standard errors to the trace
!         file
!
          call wrtfit( '       FEFIT Dummy variables' )
          call newlit
          call wrtfit( '       Case         Estimate      Std. Err.' )
          call newlit
          i = 0
!
          do icase = 1,nsub(1)
              i = i+1
              write (outbuf,'(3f13.3)') x(i,icpos(1)),gamma(i),gammse(i)
              call wrtfit(outbuf)
!
              do iobs = 2,it(1,icase)
                  i = i+1
                  write (outbuf,'(13x,2f13.3)') gamma(i),gammse(i)
                  call wrtfit(outbuf)
              end do
!
          end do
!
      end if
!
!---- reset total number of levels of all parameters in model
      nilev = nilev+ninali
!---- get covariance matrix, adjusted for intrinsically aliased
!---- parameters
!---- number of entries in covariance matrix plus one
      k = npar*(npar+1)/2 + 1
!---- number of entries in covariance matrix relating to
!---- non-intrinsically aliased parameters plus one
      l = nest*(nest+1)/2 + 1
!
      do 780 i = npar,1,-1
!
          do 779 j = i,1,-1
              k = k-1
!
!------------ covariance of non-intrinsically aliased parameter pair
              if (idnint(alias(i)) == 0 .and. idnint(alias(j)) == 0) &
              then
                  l = l-1
                  cov(k) = cov(l)
!------------ set covariance of parameter pair, at least one of which is
!------------ intrinsically aliased, to zero
              else
                  cov(k) = 0
              end if
!
  779     end do
!
  780 end do
!
!---- get beta and extrinsic alias vectors, adjusted for intrinsically
!---- aliased parameters
!---- number of parameters not intrinsically aliased
      k = nest
!
      do 781 i = npar,1,-1
!
          if (idnint(alias(i)) == 0) then
!------------ estimate of non-intrinsically aliased parameter
              beta(i) = beta(k)
!------------ extrinsic alias of non-intrinsically aliased parameter
              alias(maxpar+i) = alias(maxpar+k)
!------------ update number of remaining non-intrinsically aliased
!------------ parameters
              k = k-1
          else
!------------ set intrinsically aliased parameter to zero
              beta(i) = 0
!------------ set extrinsic alias of intrinsically aliased parameter to
!------------ zero
              alias(maxpar+i) = 0
          end if
!
  781 end do
!
      do 785 i = 1,nilev
!
          if (idnint(alias(i)) == -1) then
!
!------------ move intrinsically aliased variables back from end of
!------------ X-vector to their correct positions
!-------------***********
              call shuffl(x,ncol,i,ncol,nmes,ipos,maxcol,1,nvar)
!-------------***********
!
          end if
!
  785 end do
!
!---- combine intrinsic and extrinsic aliasing indicator vectors
      do 790 i = 1,npar
!
          if (idnint(alias(maxpar+i)) == 1) then
!------------ set extrinsically aliased parameter to zero
              beta(i) = 0
              alias(i) = 1
          end if
!
  790 end do
!
!---- initialise degrees of freedom
      ilevdf = 0
      icutdf = 0
      inormdf = nnorm
      iredf = 0
      ienddf = 0
!
!---- calculate degrees of freedom for explanatory variable parameters
      do 800 j = 1,nilev
          ilevdf = ilevdf + 1 - abs(idnint(alias(j)))
  800 end do
!
      do 805 j = nilev+1,nilev+ncut
!
          if (beta(j) == 0) then
              alias(j) = 1
          end if
!
          icutdf = icutdf + 1 - abs(idnint(alias(j)))
  805 end do
!
      do 810 j = nilev+ncut+nnorm+1,npar-nend
          iredf = iredf + 1 - abs(idnint(alias(j)))
  810 end do
!
!---- left endpoint free and not aliased, so update degrees of freedom
      if (iend(1) == 1 .and. ((endind == 'b' .and. &
      idnint(alias(npar-1)) == 0) .or. (endind == 'l' .and. &
      idnint(alias(npar)) == 0))) then
          ienddf = 1
      end if
!
!---- right endpoint free and not aliased, so update degrees of freedom
      if (iend(2) == 1 .and. idnint(alias(npar)) == 0) then
          ienddf = ienddf+1
      end if
!
      itotdf = nmes-ilevdf-icutdf-inormdf-iredf-ienddf
!
      return
!
      end subroutine model
!
!***********************************************************************
!
      subroutine modchk(endind,offlag,link,corr,order,family, &
                        ifail,univar,cname,yname,rname,offnme,nlevel, &
                        irvar,ilfit,bivar,eqscale,cquad)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character endind,link(3),family(3),cquad
      character(len=50) :: cname(2),yname,rname,offnme(3)
      integer nlevel,ilfit
      logical offlag(3),corr,order,ifail,univar,irvar,bivar,eqscale
!-----------------------------------------------------------------------
!     Function : Checks the settings of the model constants.
!-----------------------------------------------------------------------
      include 'chans.h'
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      ifail = .true.
!
      if (ilfit /= 1 .and. (cname(1) == yname .or. &
      (nlevel == 2 .and. cname(2) == yname))) then
          call wrtlin( &
          '    *** WARNING *** CASE VARIABLE SAME AS Y-VARIATE')
      else if (ilfit /= 1 .and. (cname(1) == rname .or. &
      (nlevel == 2 .and. cname(2) == rname))) then
          call wrtlin( &
          '    *** WARNING *** CASE VARIABLE SAME AS R-VARIATE')
      else if (ilfit /= 1 .and. ((cname(1) == offnme(1) .or. &
      (nlevel == 2 .and. cname(2) == offnme(1))) .or. &
      cname(1) == offnme(2) .or. cname(1) == offnme(3))) then
          call wrtlin( &
          '    *** WARNING *** CASE VARIABLE SAME AS OFFSET')
      else if (nlevel == 2 .and. cname(1) == cname(2)) then
          call wrtlin( &
          '    *** WARNING *** BOTH CASE VARIABLES THE SAME')
      end if
!
      if (yname == rname) then
          call wrtlin( &
          '    *** WARNING *** Y-VARIATE SAME AS R-VARIATE')
      else if (yname == offnme(1) .or. yname == offnme(2) .or. &
      yname == offnme(3)) then
          call wrtlin( &
          '    *** WARNING *** Y-VARIATE SAME AS OFFSET')
      end if
!
      if (irvar .and. (rname == offnme(1) .or. rname == offnme(2) .or. &
      rname == offnme(3))) then
          call wrtlin( &
          '    *** WARNING *** R-VARIATE SAME AS OFFSET')
      end if
!
      if (family(1) /= 'b' .and. family(1) /= 'o' .and. &
      link(1) /= 'l') then
          call wrtlin( &
          '    *** ERROR *** NON-BINOMIAL/PROBIT OR CLOGLOG')
          return
      end if
!
      if (bivar .and. ((family(1) == 'o' .and. family(2) /= 'o') .or. &
      (family(1) /= 'o' .and. family(2) == 'o'))) then
          call wrtlin( &
          '    *** ERROR *** ORDERED/NON-ORDERED NOT ALLOWED')
          return
      end if
!
!      if (bivar .and. family(1) /= 'b' .and. family(2) /= 'b' .and. &
!      order) then
!          call wrtlin('    *** ERROR *** NON-BINOMIAL/ORDERED')
!          return
!      end if
!
      if (family(1) == 'p' .and. endind == 'b') then
          call wrtlin('    *** ERROR *** POISSON/ENDPOINTS')
          return
      else if (family(1) == 'p' .and. endind == 'r') then
          call wrtlin('    *** ERROR *** POISSON/RIGHT ENDPOINT')
          return
      end if
!
      if (family(1) == 'g' .and. endind == 'b') then
          call wrtlin('    *** ERROR *** LINEAR/ENDPOINTS')
          return
      else if (family(1) == 'g' .and. endind == 'l') then
          call wrtlin('    *** ERROR *** LINEAR/LEFT ENDPOINT')
          return
      else if (family(1) == 'g' .and. endind == 'r') then
          call wrtlin('    *** ERROR *** LINEAR/RIGHT ENDPOINT')
          return
      end if
!
      if ((.not. univar .or. order) .and. endind /= 'n') then
          call wrtlin('    *** ERROR *** ENDPOINTS!')
          return
      end if
!
!      if ((.not. univar .or. order) .and. offlag) then
!      if (order .and. offlag) then
!          call wrtlin('    *** ERROR *** OFFSET!')
!          return
!      end if
!
!      if ((univar .or. order) .and. corr) then
!          call wrtlin('    *** ERROR *** CORRELATED!')
!          return
!      end if
!
!      if ((univar .or. order) .and. family(2) /= 'b') then
!          call wrtlin('    *** ERROR *** 2nd FAMILY!')
!          return
!      end if
!
!      if ((univar .or. order) .and. link(2) /= 'l') then
!          call wrtlin('    *** ERROR *** 2nd LINK!')
!          return
!      end if
!
!      if ((univar .or. order) .and. family(3) /= 'b') then
!          call wrtlin('    *** ERROR *** 3rd FAMILY!')
!          return
!      end if
!
!      if ((univar .or. order) .and. link(3) /= 'l') then
!          call wrtlin('    *** ERROR *** 3rd LINK!')
!          return
!      end if
!
      if (eqscale .and. .not. bivar) then
          call wrtlin('    *** ERROR *** Equal scale parameters '// &
                      'available for bivariate models only')
          return
      end if
!
      if (ilfit == 3 .and. nlevel > 1) then
          call wrtlin('    *** ERROR *** Fixed Effects available for '// &
                      'univariate models only')
          return
      end if
!
      if (nlevel == 2 .and. cquad /= 'g') then
!
          if (num_processors > 1) then
              call wrtlin('    *** WARNING *** Parallel processing '// &
                          'not available for AQ 3-level models')
          else
              call wrtlin('    *** WARNING *** Partial adaptive '// &
                          'quadrature, on level-3 only')
          end if
!
      end if
!
      ifail = .false.
!
      return
!
      end subroutine modchk
!
!***********************************************************************
!
      subroutine fitchk(mnames,iyvar,icvar,line,ni,nilev,ilev,name,nvar, &
                        yname,cname,idel,ilfit,offnme,n1var,n2var,n1lev, &
                        n2lev,order,iname,univar)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      integer, parameter :: SABRE_CMDLINE_MAX = 132
      character(len=SABRE_CMDLINE_MAX) :: line
      character(len=50) :: mnames(maxvar),name(maxvar),yname,cname(2), &
                   offnme(3),iname(3)
      integer ni,nvar,idel(maxvar),nilev,ilev(maxvar),ilfit,n1lev,n2lev, &
              n1var,n2var
      logical iyvar,icvar(2),order,univar
!-----------------------------------------------------------------------
!     Function : Checks the syntax of the FIT command.
!-----------------------------------------------------------------------
!     NI       : number of variables in a fit list
!     NILEV    : number of levels in the variables in a fit list
!     IDEL(.)  : model variables to be removed as indexed by MNAMES;
!              : if an error is detected, then IDEL(1) = -1
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=50) :: newnme(maxvar)
      integer istpos,ienpos,i,j,dummy,ixlev
      logical arg,error,null
!-----------------------------------------------------------------------
!---- y-variate not set
      if (.not. iyvar) then
          call wrtlin('    *** ERROR *** Y-VARIATE NOT SET')
          go to 290
      end if
!
!---- case-variate not set for mixture model
      if (ilfit /= 1 .and. .not. icvar(1)) then
          call wrtlin('    *** ERROR *** CASE VARIABLE NOT SET')
          go to 290
      end if
!
!---- get position of first argument
!-----*********
      call next(istpos,ienpos,4,line,arg)
!-----*********
!
      null = .false.
!
      if (.not. order .and. .not. arg) then
          call wrtlin('    *** ERROR *** NO VARIABLE NAMES GIVEN')
          go to 290
      else if (.not. arg) then
          null = .true.
      end if
!
!---- initialise number of levels of X-variables
      nilev = 0
      n1lev = 0
      n2lev = 0
!
!---- initialise model names
      do 10 i = 1,maxvar
          mnames(i) = ' '
   10 end do
!
      do 20 i = 1,maxvar
!-------- initialise error indicator
          idel(i) = 0
!
          if (null) then
              go to 30
          else if (i >= 2) then
              dummy = ienpos+1
!
!------------ get position of next argument
!-------------*********
              call next(istpos,ienpos,dummy,line,arg)
!-------------*********
!
              if (.not. arg) then
                  go to 30
              end if
!
          end if
!
!-------- name of variate to be included in model
          newnme(i) = line(istpos:ienpos)
   20 end do
!
!---- number of variates included in model
   30 ni = i-1
!
!---- check for repeated variable names
!-----***********
      call repeat(newnme,ni,error)
!-----***********
!
      if (error) then
          go to 290
      end if
!
      do 70 i = 1,ni
!
          do 40 j = 1,nvar
!
              if (newnme(i) == name(j)) then
                  go to 45
              end if
!
   40     end do
!
          if (newnme(i)(1:1) == '+') then
              call wrtlin('    *** ERROR *** INVALID "+" SYNTAX. '// &
                          'Use "FIT variable_list"')
          else if (newnme(i)(1:1) == '-') then
              call wrtlin('    *** ERROR *** INVALID "-" SYNTAX. '// &
                          'Use "FIT variable_list"')
          else if (newnme(i)(1:1) == '.') then
              call wrtlin('    *** ERROR *** INVALID "." SYNTAX. '// &
                          'Use "FIT variable_list"')
          else
!------------ non-existent variable name
              call wrtlin('    *** ERROR *** '//'VARIABLE `'// &
                          newnme(i)(1:12)//'` DOES NOT EXIST')
          end if
!
          go to 290
!
!-------- y-variate declared as explanatory variable
   45     if (newnme(i) == yname) then
              call wrtlin('    *** ERROR *** '// &
                          'Y-VARIATE NOT ALLOWED IN MODEL')
              go to 290
!-------- case variable declared as explanatory variable
          else if (newnme(i) == cname(1)) then
              call wrtlin('    *** ERROR *** '// &
                          'CASE VARIABLE NOT ALLOWED IN MODEL')
              go to 290
!-------- constant declared in ordered response model
          ELSE IF (univar .and. NEWNME(I) == iNAME(1) .and. order .and. &
          ni > 1) THEN
              call wrtlin('    *** ERROR *** CONSTANT NOT ALLLOWED '// &
              'IN ORDERED RESP. MODEL WITH COVARIATES')
              go to 290
!-------- constant-only ordered response model
          else if (newnme(i) == iname(1) .and. order .and. ni == 1) &
          then
              call wrtlin('    *** WARNING *** CONSTANT-ONLY ORDERED '// &
              'RESPONSE MODELS FIT ALL THE CUTPOINTS')
!              nilev = 0
!              go to 70
              go to 300
!              go to 290
!-------- offset declared as explanatory variable
          else if (newnme(i) == offnme(1) .or. newnme(i) == offnme(2) &
          .or. newnme(i) == offnme(3)) then
              call wrtlin('    *** ERROR *** '// &
                          'offsets not allowed as covariates in model')
              go to 290
          end if
!
!-------- get number of levels of current variable
!---------***********
          call fndlev(nvar,newnme(i),name,ixlev,ilev)
!---------***********
!
!-------- update number of levels of all variates included in the model
          nilev = nilev+ixlev
!
          if (i <= n1var) then
              n1lev = n1lev+ixlev
          else if (i <= n1var+n2var) then
              n2lev = n2lev+ixlev
          end if
!
!-------- save name of variate
          mnames(i) = newnme(i)
!-------- variate specified OK
   70 end do
!
      return
!
!---- error exit
  290 idel(1) = -1
      ilfit = -1
!
!---- cancel all model settings
      do 291 i = 1,maxvar
          mnames(i) = ' '
  291 end do
!
  300 ni = 0
      nilev = 0
      n1lev = 0
      n2lev = 0
!
      return
!
      end subroutine fitchk
!
!***********************************************************************
!
      subroutine lfit(x,work,nilev,nmes,maxcol,y,beta,alias,cov,ncov, &
                      ilfit,xll,con,tol,niter,ifail,arith,link,bivar, &
                      risk,robust,sig_e,n1lev,family,trivar,n2lev,inorm, &
                      univar,depend)
      use accurate_arithmetic
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character arith,link(3),family(3)
      integer maxcol,nmes,nilev,ncov,ilfit,niter,risk(nmes),n1lev,n2lev, &
              inorm
      double precision y(nmes),x(nmes,maxcol),beta(nilev),tol,xll,con, &
             cov(ncov),alias(nilev),work(2*nmes+nilev*(nilev+2)), &
             sig_e(3)
      logical ifail,robust,bivar,trivar,univar,depend
!-----------------------------------------------------------------------
!     Function : Performs a standard GLIM-type analysis. For the LFIT
!                command this is the fitted model, whereas for the FIT
!                command it simply provides starting values for the main
!                mixture model fitting algorithm.
!-----------------------------------------------------------------------
!     Programmer: Brian Francis
!     History: 27/11/88 Initial version
!              30/12/88 New version with beta array
!              27/01/89 Double precision version
!               5/12/89 Intrinsic aliasing added
!              19/12/89 Vital components added (Jon Barry)
!              03/04/96 Final version for SABRE 3.0 release (Dave Stott)
!              02/11/05 Version 4.0
!                    07 Version 5.0
!                 06/08 Version 6.0
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character(len=80) :: outbuf
      double precision oldlik,expeta,eta,mu,var,diff,loglik,llik,dpieta, &
             pi,scores(nmes,nilev),score(nilev),covrob(ncov),sigma(3), &
             ploglk(3),yfac,loglk(3)
      integer iter,i,j,ir,nn,n1,n2,n3,na,nmes1,nmes2,nmes3,clock1
!-----------------------------------------------------------------------
      if (qflag) then
!
!---------**********
          call clock('LFIT',clock1,1)
!---------**********
!
      end if
!
      oldlik = 1.0d20
      iter = 1
      ploglk(1) = 0
      ploglk(2) = 0
      ploglk(3) = 0
      nmes1 = 0
      nmes2 = 0
      nmes3 = 0
!
      do 20 i = 1,nmes
          ir = risk(i)
!
          if (depend) then
              ir = 1
          end if
!
          if (ir == 1) then
              nmes1 = nmes1+1
          else if (ir == 2) then
              nmes2 = nmes2+1
          else
              nmes3 = nmes3+1
          end if
!
!-------- starting values for iteratively reweighted least squares, mu~
          if (family(ir) == 'b') then
              work(i) = 0
          else if (family(ir) == 'g') then
              work(i) = 1
          else
              work(i) = y(i)
          end if
!
          if (family(ir) /= 'g' .and. y(i) == 0) then
              work(nmes+i) = -(log(3d0) + 1.333333)
          else if (family(ir) == 'b') then
              work(nmes+i) = log(3d0) + 1.333333
          else if (family(ir) == 'g') then
              work(nmes+i) = y(i)
          else
              work(nmes+i) = log(y(i))
          end if
!
          if (family(ir) == 'p') then
              yfac = 0
!
              do 1420 j = 1,nint(y(i))
                  yfac = yfac + log(dble(j))
 1420         end do
!
              ploglk(ir) = ploglk(ir) + yfac
          end if
!
   20 end do
!
!---- mixture model: initial results from logistic/log-linear regression
      if (ilfit == 0) then
          call newlin
          call wrtlin('    Initial Homogeneous Fit:')
      end if
!
      call newlin
      call wrtlin('    Iteration       Log. lik.       Difference')
      call wrtlin('    __________________________________________')
      call wterse('    Iteration       Log. lik.       Difference')
      call wterse('    __________________________________________')
!
!---- end of initialisation; start of main loop - calculate beta vector
!-----********
   30 call wls(x,nilev,nmes,work,work(nmes+1),beta,alias,work(2*nmes+1), &
               cov,ncov,tol,ifail,arith)
!-----********
!
      if (ifail) then
          return
      end if
!
      loglk(1) = 0
      loglk(2) = 0
      loglk(3) = 0
!
      if (family(1) == 'p') then
          loglk(1) = -ploglk(1)
      end if
!
      if (family(2) == 'p') then
          loglk(2) = -ploglk(2)
      end if
!
      if (family(3) == 'p') then
          loglk(3) = -ploglk(3)
      end if
!
      loglik = loglk(1) + loglk(2) + loglk(3)
!
      do 25 j = 1,nilev
          score(j) = 0
   25 end do
!
      do 40 i = 1,nmes
          ir = risk(i)
!
          if (depend) then
              ir = 1
          end if
!
          n1 = 1
          n2 = nilev
!
          if (.not. univar .and. ir == 1) then
              n2 = n1lev
          else if (bivar) then
              n1 = n1lev+1
          else if (trivar .and. ir == 2) then
              n1 = n1lev+1
              n2 = n1lev+n2lev
          else if (trivar) then
              n1 = n1lev+n2lev+1
          end if
!
          eta = 0
!
!-------- calculate beta'x_it
          do 37 j = n1,n2
              eta = eta + beta(j)*x(i,j)
   37     end do
!
!-------- calculate exp(beta'x_it)
          if (eta > zl1) then
              expeta = z1
          else if (eta < zl2) then
              expeta = z2
          else
              expeta = exp(eta)
          end if
!
!-------- pi_i = exp(beta'x_it)/[1 + exp(beta'x_it)], w_i = pi_i(1-pi_i)
!-------- -(log-likelihood) = ln[1 + exp(beta'x_it)] - y_it.beta'x_it
          if (family(ir) == 'b') then
!
              if (link(ir) == 'c') then
                  mu = 1 - exp(-expeta)
              else if (link(ir) == 'p') then
                  mu = erfc(-eta/sqrt(2d0))/2
              else
                  mu = expeta/(1+expeta)
              end if
!
              if (mu == 1) then
                  mu = 1 - 1.0d-16
              end if
!
              var = mu*(1-mu)
              llik = y(i)*log(mu) + (1 - y(i))*log(1-mu)
!-------- pi_i = exp(beta'x_it), w_i = pi_i
!-------- -(log-likelihood) = exp(beta'x_it) - y_it.beta'x_it
          else if (family(ir) == 'p') then
              mu = expeta
              var = mu
              llik = y(i)*log(mu) - mu
          else
              mu = eta
              var = 1
              llik = -(y(i) - mu)**2/2
          end if
!
          pi = 3.14159265
!
          if (link(ir) == 'c') then
              dpieta = -(1-mu)*log(1-mu)
          else if (link(ir) == 'p') then
              dpieta = exp(-eta**2/2)/sqrt(2*pi)
          else
              dpieta = var
          end if
!
!-------- store w_i for next iteration
          work(i) = dpieta**2/var
!-------- store adjusted dependent variate Z = [z_i] for next iteration
          work(nmes+i) = eta + (y(i) - mu)/dpieta
!
!-------- update log-likelihood
          loglik = loglik+llik
          loglk(ir) = loglk(ir) + llik
!
          do 35 j = 1,nilev
              scores(i,j) = (y(i) - mu)*x(i,j)
              score(j) = score(j) + scores(i,j)
   35     end do
!
   40 end do
!
      na = 0
!
      do 45 i = 1,nilev
!
          if (idnint(alias(i)) == 1) then
              na = na+1
          end if
!
   45 end do
!
      n1 = n1lev
      n2 = n2lev
!
      if (univar) then
          nn = 1
          n1 = nilev-na
      else if (bivar) then
          nn = 2
          n2 = nilev-n1lev
      else
          nn = 3
          n3 = nilev-n1lev-n2lev
      end if
!
      if (inorm /= 0) then
!
          if (family(1) == 'g') then
              sig_e(1) = sqrt(-2*loglk(1)/(nmes1-n1))
              sigma(1) = sqrt(-2*loglk(1)/nmes1)
              loglk(1) = -nmes1*(1 + log(2*pi*sigma(1)**2))/2
          end if
!
          if (.not. univar .and. family(2) == 'g') then
              sig_e(2) = sqrt(-2*loglk(2)/(nmes2-n2))
              sigma(2) = sqrt(-2*loglk(2)/nmes2)
              loglk(2) = -nmes2*(1 + log(2*pi*sigma(2)**2))/2
          end if
!
          if (trivar .and. family(3) == 'g') then
              sig_e(3) = sqrt(-2*loglk(3)/(nmes3-n3))
              sigma(3) = sqrt(-2*loglk(3)/nmes3)
              loglk(3) = -nmes3*(1 + log(2*pi*sigma(3)**2))/2
          end if
!
          loglik = loglk(1)
!
          if (.not. univar) then
              loglik = loglik + loglk(2)
!
              if (trivar) then
                  loglik = loglik + loglk(3)
              end if
!
          end if
!
      end if
!
      n1 = n1lev*(n1lev+1)/2
      n2 = (n1lev+n2lev)*(n1lev+n2lev+1)/2
!
      if (univar) then
          n1 = ncov
      else if (bivar) then
          n2 = ncov
      end if
!
      do 50 i = 1,ncov
!
          if (family(1) == 'g' .and. i <= n1) then
              cov(i) = cov(i)*sig_e(1)**2
          end if
!
          if (.not. univar .and. family(2) == 'g' .and. i > n1 &
          .and. i <= n2) then
              cov(i) = cov(i)*sig_e(2)**2
          end if
!
          if (trivar .and. family(3) == 'g' .and. i > n2) then
              cov(i) = cov(i)*sig_e(3)**2
          end if
!
   50 end do
!
!---- calculate change in deviance between current and previous
!---- iteration
      diff = loglik-oldlik
!
      if (iter == 1) then
          write (outbuf,'(a,i5,a,g20.8)') '    ',iter,'    ',loglik
          call wrtlin(outbuf)
          call wterse(outbuf)
      else
          write (outbuf,'(a,i5,a,g20.8,g13.4)') '    ',iter,'    ', &
          loglik,diff
          call wrtlin(outbuf)
          call wterse(outbuf)
      end if
!
!---- not yet converged
!      IF (((univar .and. inorm /= 1) .or. (bivar .and. inorm /= 2) &
!      .or. (trivar .and. inorm /= 3)) .and. abs(DIFF) >= CON .AND. &
!      ITER < NITER) THEN
      if (((univar .and. inorm /= 1) .or. (bivar .and. inorm /= 2) &
      .or. (trivar .and. inorm /= 3)) .and. abs(diff) >= con) then
!-------- update number of iterations
          iter = iter+1
!
!-------- no more iterations allowed
          if (iter > niter) then
              call wrtlin('    *** WARNING *** The maximum number of '// &
                          'iterations has been reached')
              call newlin
!
              if (ilfit == 1) then
                  ifail = .true.
              end if
!
              return
          end if
!
!-------- save current deviance value
          oldlik = loglik
!-------- next weighted least squares iteration to solve
!-------- beta^ = (X'WX)`(X'WZ)
          ifail = .true.
          go to 30
!      ELSE IF (((univar .and. inorm /= 1) .or. &
!      (bivar .and. inorm /= 2) .or. (trivar .and. inorm /= 3)) .and. &
!      abs(DIFF) >= CON .AND. ITER == NITER) THEN
!          call wrtlin('    *** WARNING *** '// &
!                      'THE IRLS ALGORITHM HAS NOT CONVERGED')
      end if
!
      call newlin
!
!---- save deviance
      xll = loglik
!
      if (inorm == 0) then
!
!---------********
          call wls(x,nilev,nmes,work,work(nmes+1),beta,alias, &
                   work(2*nmes+1),cov,ncov,tol,ifail,arith)
!---------********
!
          if (ifail) then
              return
          end if
!
      end if
!
      if (robust) then
!
!---------***********
          call robind(scores,nmes,nilev,cov,ncov,covrob)
!---------***********
!
          do 90 i = 1,ncov
              cov(i) = covrob(i)
   90     end do
!
      end if
!
      if (qflag) then
!
!---------**********
          call clock('',clock1,2)
!---------**********
!
      end if
!
      return
!
      end subroutine lfit
!
!***********************************************************************
!
      subroutine wls(x,nilev,nmes,w,z,beta,alias,wk,cov,ncov,tol,ifail, &
                     arith)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character arith
      integer nilev,nmes,ncov
      double precision x(nmes,nilev),z(nmes),w(nmes),tol,beta(nilev), &
             cov(ncov),wk(nilev*(nilev+2)),alias(nilev)
      logical ifail
!-----------------------------------------------------------------------
!     Function : Performs a single iteration of the IRLS algorithm.
!-----------------------------------------------------------------------
!     W(.)  : weights for iteratively reweighted least squares algorithm
!     Z(.)  : adjusted dependent variates for irls algorithm
!     WK(.) : tail of the X array used as workspace
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      integer i,j,k,l
      logical myturn
!-----------------------------------------------------------------------
!     wk_all is a work array to store the accumulated wk's since each
!     process assigns only a portion of its copy of wk
      double precision wk_all(maxvar*(maxvar+1)/2)
!-----------------------------------------------------------------------
!     MPI parameter definitions
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      integer ierror
!-----------------------------------------------------------------------
!---- calculate X'WX, where X = [x_ij] and W = diag[w_i]
!---- X'WX is a symmetric matrix, so only lower triangular section is
!---- stored as (wk_1,...,wk_n(n+1)/2)
      l = 0
!
      do 12 i = 1,nilev
!
          do 10 j = 1,i
              l = l+1
              wk(l) = 0
!
              myturn = this_processor == mod(l-1,num_processors)
!
              if (.not. myturn) then
                  go to 10
              end if
!
              do 8 k = 1,nmes
                  wk(l) = wk(l) + x(k,j)*w(k)*x(k,i)
    8         end do
!
   10     end do
!
   12 end do
!
      if (num_processors > 1) then
!
!-------- Reduce the partially assigned wk's into wk_all then copy this
!-------- back to wk so that all of wk is then assigned in each copy of
!         the program
!---------******************
          call mpi_allreduce(wk,wk_all,nilev*(nilev+1)/2, &
                             mpi_double_precision,mpi_sum, &
                             mpi_comm_world,ierror)
!---------******************
!
          do 13 l = 1,nilev*(nilev+1)/2
              wk(l) = wk_all(l)
 13       end do
!
      end if
!
!---- invert X'WX and store lower triangular section of this symmetric
!---- matrix as (wk_n(n+1)/2 + 1,...,wk_n(n+1))
!-----***********
      call syminv(wk,nilev,wk(ncov+1),wk(2*ncov+1),alias,ncov,tol,ifail, &
                  arith)
!-----***********
!
      if (ifail) then
          call wrtlin('    *** ERROR *** '// &
                      'HESSIAN MATRIX NOT POSITIVE SEMI-DEFINITE')
          return
      end if
!
!---- store lower triangular section of symmetric covariance matrix
      do 15 i = 1,ncov
          cov(i) = wk(ncov+i)
   15 end do
!
!---- calculate X'WZ, where Z = [z_i], and store this vector as
!---- (wk_1,...,wk_n)
      do 20 i = 1,nilev
          wk(i) = 0
!
          myturn = this_processor == mod(i-1,num_processors)
!
          if (.not. myturn) then
              go to 20
          end if
!
          do 18 j = 1,nmes
              wk(i) = wk(i) + x(j,i)*w(j)*z(j)
   18     end do
!
   20 end do
!
      if (num_processors > 1) then
!
!-------- Reduce the partially assigned wk's into wk_all then copy this
!-------- back to wk so that all of wk is then assigned in each copy of
!         the program
!-------- Note that only nilev wk's are being assigned this time
!---------******************
          call mpi_allreduce(wk,wk_all,nilev,mpi_double_precision, &
                             mpi_sum,mpi_comm_world,ierror)
!---------******************
!
          do 23 l = 1,nilev
              wk(l) = wk_all(l)
 23       end do
!
      end if
!
      do 25 i = 1,nilev
          beta(i) = 0
   25 end do
!
      i = 1
      j = 1
!
!---- calculate BETA = (X'WX)`(X'WZ), where ` signifies matrix inversion
      do 30 k = ncov+1,2*ncov
          beta(i) = beta(i) + wk(j)*wk(k)
!
          if (i == j) then
              i = 1
              j = j+1
          else
              beta(j) = beta(j) + wk(i)*wk(k)
              i = i+1
          end if
!
   30 end do
!
      return
!
      end subroutine wls
!
!***********************************************************************
!
      subroutine syminv(a,nilev,u,wk,alias,ncov,tol,ifail,arith)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character arith
      integer nilev,ncov
      double precision tol,alias(nilev),a(ncov),u(ncov),wk(nilev)
      logical ifail
!-----------------------------------------------------------------------
!     Function : Inverts the lower triangular symmetric A matrix into U.
!-----------------------------------------------------------------------
!     SYMINV and CHOL are actually Applied Statistics algorithms
!     [AS7 & AS6] adapted and corrected by us.
!     Published in the book
!     'APPLIED STATISTICS ALGORITHMS' by P. Griffiths and I.D. Hill.
!-----------------------------------------------------------------------
!     A  : lower triangular part of psd symmetric matrix to be inverted
!     U  : upper triangular part of generalised inverse of A;
!          i.e. AU = I
!     WK : workspace
!-----------------------------------------------------------------------
      double precision v
      integer i,j,k,l,irow,icol,jcol,ndiag,mdiag
!-----------------------------------------------------------------------
!---- use Cholesky decomposition to write positive semi-definite
!---- symmetric matrix A = [a_ij] in form A = U'U, where U = [c_ij] is
!---- upper triangular. Lower triangular section of A stored
!---- (a_11,a_21,a_22,a_31,...,a_n(n+1)/2)
!---- store U as (c_11,c_12,c_22,c_13,c_23,c_33,...,c_n(n+1)/2)
      if (arith == 'f') then
!
!---------***********
          call chofas(a,nilev,u,alias,ncov,tol,ifail)
!---------***********
!
      else
!
!---------***********
!          call choacc(a,nilev,u,alias,ncov,tol,ifail)
          call chofas(a,nilev,u,alias,ncov,tol,ifail)
!---------***********
!
      end if
!
      if (ifail) then
          return
      end if
!
!---- calculate generalised inverse of A. i.e. find symmetric matrix U
!---- such that AU = I. Store lower triangular section of U,
!---- (c_11,c_21,c_22,c_31,...,c_n(n+1)/2)
      irow = nilev
      ndiag = ncov
!---- check if parameter should be aliased
   10 l = ndiag
!
      if (abs(u(ndiag)) == 0) then
          go to 60
      end if
!
      do 20 i = irow,nilev
          wk(i) = u(l)
          l = l+i
   20 end do
!
      icol = nilev
      jcol = ncov
      mdiag = ncov
   30 l = jcol
      v = 0
!
      if (icol == irow) then
          v = 1/wk(irow)
      end if
!
      k = nilev
!
   40 if (k == irow) then
          go to 50
      end if
!
      v = v - wk(k)*u(l)
      k = k-1
      l = l-1
!
      if (l > mdiag) then
          l = l-k+1
      end if
!
      go to 40
!
   50 u(l) = v/wk(irow)
!
      if (icol == irow) then
          go to 80
      end if
!
      mdiag = mdiag-icol
      icol = icol-1
      jcol = jcol-1
!
      go to 30
!
   60 do 70 j = irow,nilev
          u(l) = 0
          l = l+j
   70 end do
!
   80 ndiag = ndiag-irow
      irow = irow-1
!
      if (irow /= 0) then
          go to 10
      end if
!
      return
!
      end subroutine syminv
!
!***********************************************************************
!
      subroutine chofas(a,nilev,u,alias,ncov,tol,ifail)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer nilev,ncov
      double precision tol,alias(nilev),a(ncov),u(ncov)
      logical ifail
!-----------------------------------------------------------------------
!     Function : Performs Cholesky decomposition on the psd matrix A.
!-----------------------------------------------------------------------
!     A : lower triangular part of psd symmetric matrix to be decomposed
!     U : upper triangular part of Cholesky factor matrix; i.e. U'U = A
!-----------------------------------------------------------------------
      double precision w,v
      integer i,j,k,l,m,irow,icol,ii,kk
!-----------------------------------------------------------------------
      j = 1
      k = 0
      ii = 0
      w = 0
!
      do 80 icol = 1,nilev
!-------- initialise extrinsic aliasing indicator
          alias(icol) = 0
          ii = ii+icol
          v = tol**2*a(ii)
          l = 0
          kk = 0
!
          do 40 irow = 1,icol
              kk = kk+irow
              k = k+1
              w = a(k)
              m = j
!
              do 10 i = 1,irow
                  l = l+1
!
                  if (i == irow) then
                      go to 20
                  end if
!
!---------------- update a_ji - sum{k=1,...,i-1}u_ki.u_kj, for j>=i
                  w = w - u(l)*u(m)
                  m = m+1
   10         end do
!
   20         if (irow == icol) then
                  go to 50
              end if
!
!------------ check if parameter should be aliased
              if (abs(u(l)) == 0) then
                  go to 30
              end if
!
!------------ calculate u_ij = (a_ji - sum{k=1,...,i-1}u_ki.u_kj)/u_ii,
!------------ for j>i
              u(k) = w/u(l)
              go to 40
!
!------------ A not positive semi-definite
   30         if (w*w > abs(v*a(kk))) then
                  return
              end if
!
!------------ A singular, so set u_ij = 0 for j>i
              u(k) = 0
   40     end do
!
   50     if (abs(w) <= abs(tol*a(k))) then
              go to 60
          end if
!
!-------- A not positive semi-definite
          if (w < 0) then
              return
          end if
!
!-------- calculate u_ii = sqrt(a_ii - sum{k=1,...,i-1}u_ki^2)
          u(k) = sqrt(w)
          go to 70
!-------- |a_ii - sum{k=1,...,i-1}u_ki^2| < TOL.|a_ii|, so set u_ii = 0
!-------- i.e. A singular
   60     u(k) = 0
!-------- set the extrinsic aliasing indicator
          alias(icol) = 1
   70     j = j+icol
   80 end do
!
      ifail = .false.
!
      return
!
      end subroutine chofas
!
!***********************************************************************
!
      subroutine choacc(a,nilev,u,alias,ncov,tol,ifail)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      integer nilev,ncov
      double precision tol,alias(nilev),a(ncov),u(ncov)
      logical ifail
!-----------------------------------------------------------------------
!     Function : Performs Cholesky decomposition on the psd matrix A.
!-----------------------------------------------------------------------
!     A : lower triangular part of psd symmetric matrix to be decomposed
!     U : upper triangular part of Cholesky factor matrix; i.e. U'U = A
!-----------------------------------------------------------------------
      double precision w,tolsq,dummy,u2,wman,aman(281625),vman,atol, &
             atoman,wsqman,wsq,vaman,va,uman(281625)
      integer i,j,k,l,m,irow,icol,ii,kk,aexp(281625),vexp,dumexp, &
              uexp(281625),u2exp,atoexp,wsqexp,vaexp,wexp
!-----------------------------------------------------------------------
      j = 1
      k = 0
      ii = 0
!
      do 80 icol = 1,nilev
!-------- initialise extrinsic aliasing indicator
          alias(icol) = 0
          ii = ii+icol
          tolsq = tol**2
          aman(ii) = a(ii)
          aexp(ii) = 0
!
!---------***********
          call manexp(aman(ii),aexp(ii))
!---------***********
!
!---------***********
          call multip(tolsq,0,aman(ii),aexp(ii),vman,vexp)
!---------***********
!
          l = 0
          kk = 0
!
          do 40 irow = 1,icol
              kk = kk+irow
              k = k+1
              w = a(k)
              wman = w
              wexp = 0
!
!-------------***********
              call manexp(wman,wexp)
!-------------***********
!
              m = j
!
              do 10 i = 1,irow
                  l = l+1
!
                  if (i == irow) then
                      go to 20
                  end if
!
!---------------- update a_ji - sum{k=1,...,i-1}u_ki.u_kj, for j>=i
                  uman(l) = u(l)
                  uexp(l) = 0
!
!-----------------***********
                  call manexp(uman(l),uexp(l))
!-----------------***********
!
                  uman(m) = u(m)
                  uexp(m) = 0
!
!-----------------***********
                  call manexp(uman(m),uexp(m))
!-----------------***********
!
!-----------------***********
                  call multip(uman(l),uexp(l),uman(m),uexp(m),u2,u2exp)
!-----------------***********
!
                  dummy = wman
                  dumexp = wexp
!
!-----------------***********
                  call addnum(dummy,dumexp,-u2,u2exp,wman,wexp)
!-----------------***********
!
!-----------------***********
                  call normal(wman,wexp,w)
!-----------------***********
!
                  m = m+1
   10         end do
!
   20         if (irow == icol) then
                  go to 50
              end if
!
!------------ check if parameter should be aliased
              if (abs(u(l)) == 0) then
                  go to 30
              end if
!
!------------ calculate u_ij = (a_ji - sum{k=1,...,i-1}u_ki.u_kj)/u_ii,
!------------ for j>i
!
!-------------***********
              call multip(wman,wexp,1d0/uman(l),-uexp(l),uman(k), &
                          uexp(k))
!-------------***********
!
!-------------***********
              call normal(uman(k),uexp(k),u(k))
!-------------***********
!
              go to 40
!
!------------ A not positive semi-definite
!
!-------------***********
   30         call multip(wman,wexp,wman,wexp,wsqman,wsqexp)
!-------------***********
!
!-------------***********
              call normal(wsqman,wsqexp,wsq)
!-------------***********
!
!-------------***********
              call multip(vman,vexp,aman(kk),aexp(kk),vaman,vaexp)
!-------------***********
!
!-------------***********
              call normal(vaman,vaexp,va)
!-------------***********
!
              if (wsq > abs(va)) then
                  return
              end if
!
!------------ A singular, so set u_ij = 0 for j>i
              u(k) = 0
              uman(k) = u(k)
              uexp(k) = 0
   40     end do
!
!---------***********
   50     call multip(tol,0,aman(k),aexp(k),atoman,atoexp)
!---------***********
!
!---------***********
          call normal(atoman,atoexp,atol)
!---------***********
!
          if (abs(w) <= abs(atol)) then
              go to 60
          end if
!
!-------- A not positive semi-definite
          if (w < 0) then
              return
          end if
!
!-------- calculate u_ii = sqrt(a_ii - sum{k=1,...,i-1}u_ki^2)
          u(k) = sqrt(w)
          uman(k) = u(k)
          uexp(k) = 0
!
!---------***********
          call manexp(uman(k),uexp(k))
!---------***********
!
          go to 70
!
!-------- |a_ii - sum{k=1,...,i-1}u_ki^2| < TOL.|a_ii|, so set u_ii = 0
!-------- i.e. A singular
   60     u(k) = 0
          uman(k) = u(k)
          uexp(k) = 0
!-------- set the extrinsic aliasing indicator
          alias(icol) = 1
   70     j = j+icol
   80 end do
!
      ifail = .false.
!
      return
!
      end subroutine choacc
!
!***********************************************************************
!
      subroutine fit(beta,alias,score,hess,cov,produc,work,y,x,ncov, &
                     nsub,quad,nmes,it,nm,nest,nilev,iend,con,alp,xll, &
                     tol,nmeil,niter,arith,endind,ifail,offlag,offpos, &
                     maxcol,link,risk,corr,bivar,robust,order,ncat, &
                     n1lev,family,trivar,inorm,icorr,n2lev,ilfit,univar, &
                     isca,nlevel,depend,eqscale,dfirst,cquad,maxit, &
                     maxcat)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character arith,endind,link(3),family(3),cquad
      integer nest,ncov,nmes,nilev,nlevel,nsub(2),it(2,nsub(1)),nm(3), &
              iend(2),offpos(3),maxcol,nmeil,niter,produc(nsub(1),2), &
              risk(nmes),ncat(3),n1lev,inorm,n2lev,icorr,ilfit,isca, &
              maxit,maxcat
      double precision y(nmes),tol,beta(nest),con,x(nmes,nilev), &
             quad(3,2,256),work(nest*(nest+2)),xll,cov(ncov), &
             score(nest),hess(ncov),alp,alias(nest)
      logical ifail,trivar,offlag(3),corr,bivar,robust,order,appflag, &
              univar,depend,eqscale,dfirst
!-----------------------------------------------------------------------
!     Function : Algorithm for fitting a logistic-normal/log-linear
!                normal mixture model.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      character text(2)*5,text1*4,itstr*3,clkstr*40
      character(len=80) :: outbuf
      double precision step,alpha(maxpar),delta(maxpar),gamma(maxpar), &
             oll,adjust,sing,bottom,covari(maxpar,maxpar),xxll,value, &
             diff,alph,za,ab,temp,oldadj,zb,scomax,scnorm,wll,xil_psi, &
             mu(3,0:100,nsub(1)),tau(3,0:100,nsub(1)),xall, &
             xil,xil_mu(3),xil_tau(3),taui,xal,psi(0:100,nsub(1)), &
             scores(nmes,nest),covrob(ncov),sesig12,sig12,phi(6), &
             varphi(6),betaphi(6),betasig(3),betarho(3),covphi(6,6),v1, &
             v2,v3,v4,v5,v6,sig13,sig23,sesig13,sesig23, &
             aquad(3,2,256),xild,bquad(3,2,256),xill(0:100,nsub(1)), &
             mus(3,nsub(1)),taus(3,nsub(1)),xal_total,dcon, &
             cutalpha(ncat(1))
      integer loc(maxpar),free,nr,nl,iter,i,flag,j,l1,ival,k,nlines,m, &
              iline,clock1,clock2,sign(3),nconflag(0:100),k1,k2,iquad, &
              nosign,ierror,nconflags(0:100)
      logical done,stick,chop,doadd,oldadd,recycl,orthog,wrnflg,dflag, &
              halfstep,conflag(nsub(1)),myturn,anyfail,parflag
!-----------------------------------------------------------------------
      data text /'fixed','free '/
      data xxll /0d0/, nconflag(0) /0/
!-----------------------------------------------------------------------
      include 'sabre_mpi.h'
!-----------------------------------------------------------------------
      if (qflag .and. ilfit == 0) then
          clkstr = 'FIT iteration '
          write (itstr,'(i3)') 1
          clkstr(15:17) = itstr
!
!---------**********
          call clock(clkstr(1:17),clock1,1)
!---------**********
!
          call newliq
      end if
!
      parflag = .false.
!
      if (num_processors > 1 .and. cquad /= 'g' .and. nlevel == 1) then
          parflag = .true.
      end if
!
!---- set iteration counter etc.
      oldadd = .false.
      oldadj = 1
      iter = 1
      step = 1
      alph = alp
      nl = iend(1)
      nr = iend(2)
      free = nest
!
      call newlin
!
      if (nl+nr /= 0) then
          call wrtlin('                      Log           Step'// &
                      '       Endpoints     Orthogonality')
          call wrtlin('    Iteration      likelihood      '// &
                      'length     0         1      criterion')
          call wterse('                      Log           Step'// &
                      '       Endpoints     Orthogonality')
          call wterse('    Iteration      likelihood      '// &
                      'length     0         1      criterion')
      else
          call wrtlin('                      Log           Step'// &
                      '       Gradient      Orthogonality')
          call wrtlin('    Iteration      likelihood      '// &
                      'length        norm          criterion')
          call wterse('                      Log           Step'// &
                      '       Gradient      Orthogonality')
          call wterse('    Iteration      likelihood      '// &
                      'length        norm          criterion')
      end if
!
      call wrtlin('    _________________________________'// &
                  '______________________________________')
      call wterse('    _________________________________'// &
                  '______________________________________')
!
!---- indicator for controlling details of the screen output
      wrnflg = .true.
!
      dcon = con
!
      if (dfirst) then
          dcon = sqrt(con)
      end if
!
      if (cquad == 'g' .or. nm(1) == 1 .or. (order .and. ilfit == 2)) &
      then
          iquad = 0
      else if (cquad == 'a') then
          iquad = 1
      else if (cquad == 'b') then
          iquad = 2
      else if (cquad == 'c') then
          iquad = 3
      else if (cquad == 'd') then
          iquad = 4
      end if
!
!---- Start of main iteration
   10 if (iter <= nmeil .or. dfirst) then
!-------- approximate Hessian
          flag = 2
          appflag = .true.
      else
!-------- true Hessian
          flag = 3
          appflag = .false.
      end if
!
      recycl = .false.
      halfstep = .false.
      done = .false.
!
      if (iquad == 0 .or. (iquad >= 3 .and. iter /= 1)) then
          go to 12
      end if
!
      do 1140 i = 1,nsub(nlevel)
          xill(0,i) = 0
          conflag(i) = .false.
 1140 end do
!
      m = 0
      xall = 0
 1150 m = m+1
      xal = xall
!
      if ((iquad == 1 .or. iquad == 3) .and. m == 1) then
          nconflag(m) = 0
      else if (iquad == 1 .or. iquad == 3) then
          nconflag(m) = nconflag(m-1)
      else if (iquad == 2 .or. iquad == 4) then
          nconflag(m) = nsub(nlevel)
      end if

      nconflags(m) = nconflag(m)
!
      if (tflag .or. qflag) then
          clkstr = 'AQ iteration '
          write (itstr,'(i3)') m
          clkstr(14:16) = itstr
      end if
!
      if (qflag .and. ilfit == 0) then
!
!---------**********
          call clock(clkstr(1:16),clock2,1)
!---------**********
!
      end if
!
      if (nlevel == 2) then
          go to 1502
      end if
!
      do 1500 i = 1,nsub(1)
!
          if (i == 1 .or. .not. parflag) then
              myturn = .true.
          else
              myturn = this_processor == mod(i-2,num_processors)
          end if
!
          if (.not. myturn) then
              go to 1500
          end if
!
          if (conflag(i)) then
              go to 1500
          end if
!
          if (m == 1) then
!
              if (iter == 1) then
                  mu(1,0,i) = 0
                  tau(1,0,i) = 1
              else
                  mu(1,0,i) = mus(1,i)
                  tau(1,0,i) = taus(1,i)
              end if
!
              if (.not. univar) then
!
                  if (iter == 1) then
                      mu(2,0,i) = 0
                      tau(2,0,i) = 1
                  else
                      mu(2,0,i) = mus(2,i)
                      tau(2,0,i) = taus(2,i)
                  end if
!
                  if (trivar) then
!
                      if (iter == 1) then
                          mu(3,0,i) = 0
                          tau(3,0,i) = 1
                      else
                          mu(3,0,i) = mus(3,i)
                          tau(3,0,i) = taus(3,i)
                      end if
!
                  end if
!
              end if
!
          end if
!
          do 1300 k = 1,nm(1)
!------------ adaptive quadrature locations z_k
              aquad(1,1,k) = mu(1,m-1,i) + tau(1,m-1,i)*quad(1,1,k)
!------------ adaptive quadrature weights w_k
              aquad(1,2,k) = tau(1,m-1,i)* &
              exp((quad(1,1,k)**2 - aquad(1,1,k)**2)/2)*quad(1,2,k)
              bquad(1,1,k) = aquad(1,1,k)
              bquad(1,2,k) = aquad(1,2,k)
 1300     end do
!
          if (.not. univar) then
!
              do 1305 k = 1,nm(2)
                  aquad(2,1,k) = mu(2,m-1,i) + tau(2,m-1,i)*quad(2,1,k)
                  aquad(2,2,k) = tau(2,m-1,i)*quad(2,2,k)* &
                  exp((quad(2,1,k)**2 - aquad(2,1,k)**2)/2)
                  bquad(2,1,k) = aquad(2,1,k)
                  bquad(2,2,k) = aquad(2,2,k)
 1305         end do
!
              if (trivar) then
!
                  do 1302 k = 1,nm(3)
                      aquad(3,1,k) = mu(3,m-1,i) + &
                      tau(3,m-1,i)*quad(3,1,k)
                      aquad(3,2,k) = tau(3,m-1,i)*quad(3,2,k)* &
                      exp((quad(3,1,k)**2 - aquad(3,1,k)**2)/2)
                      bquad(3,1,k) = aquad(3,1,k)
                      bquad(3,2,k) = aquad(3,2,k)
 1302             end do
!
              end if
!
          end if
!
          if ((iquad == 2 .or. iquad == 4) .and. m == 2) then
              go to 1500
          end if
!
!---------***********
          call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev,quad, &
                      nm,score,hess,ncov,1,arith,'n',ifail,offlag, &
                      offpos,maxcol,link,risk,corr,bivar,iter,1, &
                      1d0,scnorm,scores,order,ncat,n1lev,family,trivar, &
                      n2lev,ilfit,inorm,nlevel,depend,eqscale,i,xil, &
                      nosign,iquad,aquad,maxit,mus,taus,maxcat)
!---------***********
!
          if (ifail) then
              return
          end if
!
          do 1310 k = 1,nm(1)
              bquad(1,2,k) = aquad(1,1,k)*aquad(1,2,k)
 1310     end do
!
!---------***********
          call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev,quad, &
                      nm,score,hess,ncov,1,arith,'n',ifail,offlag, &
                      offpos,maxcol,link,risk,corr,bivar,iter,2,1d0, &
                      scnorm,scores,order,ncat,n1lev,family,trivar, &
                      n2lev,ilfit,inorm,nlevel,depend,eqscale,i, &
                      xil_mu(1),sign(1),iquad,bquad,maxit,mus,taus, &
                      maxcat)
!---------***********
!
          if (ifail) then
              return
          end if
!
          do 1320 k = 1,nm(1)
              bquad(1,2,k) = aquad(1,1,k)*bquad(1,2,k)
 1320     end do
!
!---------***********
          call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev,quad, &
                      nm,score,hess,ncov,1,arith,'n',ifail,offlag, &
                      offpos,maxcol,link,risk,corr,bivar,iter,1,1d0, &
                      scnorm,scores,order,ncat,n1lev,family,trivar, &
                      n2lev,ilfit,inorm,nlevel,depend,eqscale,i, &
                      xil_tau(1),nosign,iquad,bquad,maxit,mus,taus, &
                      maxcat)
!---------***********
!
          if (ifail) then
              return
          end if
!
          if (.not. univar) then
!
              do 1330 k = 1,nm(1)
                  bquad(1,2,k) = aquad(1,2,k)
 1330         end do
!
              do 1340 k = 1,nm(2)
                  bquad(2,2,k) = aquad(2,1,k)*aquad(2,2,k)
 1340         end do
!
!-------------***********
              call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev, &
                          quad,nm,score,hess,ncov,1,arith,'n',ifail, &
                          offlag,offpos,maxcol,link,risk,corr,bivar, &
                          iter,2,1d0,scnorm,scores,order,ncat,n1lev, &
                          family,trivar,n2lev,ilfit,inorm,nlevel,depend, &
                          eqscale,i,xil_mu(2),sign(2),iquad,bquad,maxit, &
                          mus,taus,maxcat)
!-------------***********
!
              if (ifail) then
                  return
              end if
!
              do 1350 k = 1,nm(2)
                  bquad(2,2,k) = aquad(2,1,k)*bquad(2,2,k)
 1350         end do
!
!-------------***********
              call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev, &
                          quad,nm,score,hess,ncov,1,arith,'n',ifail, &
                          offlag,offpos,maxcol,link,risk,corr,bivar, &
                          iter,1,1d0,scnorm,scores,order,ncat,n1lev, &
                          family,trivar,n2lev,ilfit,inorm,nlevel,depend, &
                          eqscale,i,xil_tau(2),nosign,iquad,bquad,maxit, &
                          mus,taus,maxcat)
!-------------***********
!
              if (ifail) then
                  return
              end if
!
              if (trivar) then
!
                  do 1335 k = 1,nm(2)
                      bquad(2,2,k) = aquad(2,2,k)
 1335             end do
!
                  do 1345 k = 1,nm(3)
                      bquad(3,2,k) = aquad(3,1,k)*aquad(3,2,k)
 1345             end do
!
!-----------------***********
                  call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest, &
                              nilev,quad,nm,score,hess,ncov,1,arith,'n', &
                              ifail,offlag,offpos,maxcol,link,risk,corr, &
                              bivar,iter,2,1d0,scnorm,scores,order,ncat, &
                              n1lev,family,trivar,n2lev,ilfit,inorm, &
                              nlevel,depend,eqscale,i,xil_mu(3),sign(3), &
                              iquad,bquad,maxit,mus,taus,maxcat)
!-----------------***********
!
                  if (ifail) then
                      return
                  end if
!
                  do 1355 k = 1,nm(3)
                      bquad(3,2,k) = aquad(3,1,k)*bquad(3,2,k)
 1355             end do
!
!-----------------***********
                  call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest, &
                              nilev,quad,nm,score,hess,ncov,1,arith,'n', &
                              ifail,offlag,offpos,maxcol,link,risk,corr, &
                              bivar,iter,1,1d0,scnorm,scores,order,ncat, &
                              n1lev,family,trivar,n2lev,ilfit,inorm, &
                              nlevel,depend,eqscale,i,xil_tau(3),nosign, &
                              iquad,bquad,maxit,mus,taus,maxcat)
!-----------------***********
!
                  if (ifail) then
                      return
                  end if
!
              end if
!
          end if
!
          mu(1,m,i) = sign(1)*exp(xil_mu(1) - xil)
          taui = exp(xil_tau(1) - xil) - mu(1,m,i)**2
!
          if (abs(taui) > dcon) then
              tau(1,m,i) = sqrt(abs(taui))
          else
              tau(1,m,i) = tau(1,m-1,i)/2
          end if
!
          mus(1,i) = mu(1,m,i)
          taus(1,i) = tau(1,m,i)
!
          if (.not. univar) then
              mu(2,m,i) = sign(2)*exp(xil_mu(2) - xil)
              taui = exp(xil_tau(2) - xil) - mu(2,m,i)**2
!
              if (abs(taui) > dcon) then
                  tau(2,m,i) = sqrt(abs(taui))
              else
                  tau(2,m,i) = tau(2,m-1,i)/2
              end if
!
              mus(2,i) = mu(2,m,i)
              taus(2,i) = tau(2,m,i)
!
              if (trivar) then
                  mu(3,m,i) = sign(3)*exp(xil_mu(3) - xil)
                  taui = exp(xil_tau(3) - xil) - mu(3,m,i)**2
!
                  if (abs(taui) > dcon) then
                      tau(3,m,i) = sqrt(abs(taui))
                  else
                      tau(3,m,i) = tau(3,m-1,i)/2
                  end if
!
                  mus(3,i) = mu(3,m,i)
                  taus(3,i) = tau(3,m,i)
              end if
!
          end if
!
          if (iquad == 1 .or. iquad == 3) then
              xild = abs(xil - xill(m-1,i))
!
              if (xild < dcon) then
                  conflag(i) = .true.
!
                  if (i /= 1 .or. this_processor == boss_processor) then
                      nconflag(m) = nconflag(m) + 1
                  end if
!
                  xall = xall+xil
              end if
!
              xill(m,i) = xil
          end if
!
          if (i /= 1 .or. this_processor == boss_processor) then
              xal = xal+xil
          end if
!
 1500 end do
!
      if (nlevel == 1) then
          go to 1503
      end if
!
 1502 do 1501 i = 1,nsub(2)
!
!          if (i == 1 .or. iquad > 0) then
          if (i == 1) then
              myturn = .true.
          else
              myturn = this_processor == mod(i-2,num_processors)
          end if
!
          if (.not. myturn) then
              go to 1501
          end if
!
          if (conflag(i)) then
              go to 1501
          end if
!
          if (m == 1) then
!
              if (iter == 1) then
                  mu(2,0,i) = 0
                  tau(2,0,i) = 1
              else
                  mu(2,0,i) = mus(2,i)
                  tau(2,0,i) = taus(2,i)
              end if
!
          end if
!
          do 1301 k = 1,nm(2)
!------------ adaptive quadrature locations z_k
              aquad(2,1,k) = mu(2,m-1,i) + tau(2,m-1,i)*quad(2,1,k)
!------------ adaptive quadrature weights w_k
              aquad(2,2,k) = tau(2,m-1,i)* &
              exp((quad(2,1,k)**2 - aquad(2,1,k)**2)/2)*quad(2,2,k)
              bquad(2,1,k) = aquad(2,1,k)
              bquad(2,2,k) = aquad(2,2,k)
 1301     end do
!
          if ((iquad == 2 .or. iquad == 4) .and. m == 2) then
              go to 1501
          end if
!
!---------***********
          call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev,quad, &
                      nm,score,hess,ncov,1,arith,'n',ifail,offlag, &
                      offpos,maxcol,link,risk,corr,bivar,iter,1,1d0, &
                      scnorm,scores,order,ncat,n1lev,family,trivar, &
                      n2lev,ilfit,inorm,nlevel,depend,eqscale,i,xil, &
                      nosign,iquad,aquad,maxit,mus,taus,maxcat)
!---------***********
!
          if (ifail) then
              return
          end if
!
          do 1311 k = 1,nm(2)
              bquad(2,2,k) = aquad(2,1,k)*aquad(2,2,k)
 1311     end do
!
!---------***********
          call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev,quad, &
                      nm,score,hess,ncov,1,arith,'n',ifail,offlag, &
                      offpos,maxcol,link,risk,corr,bivar,iter,2,1d0, &
                      scnorm,scores,order,ncat,n1lev,family,trivar, &
                      n2lev,ilfit,inorm,nlevel,depend,eqscale,i, &
                      xil_mu(2),sign(2),iquad,bquad,maxit,mus,taus, &
                      maxcat)
!---------***********
!
          if (ifail) then
              return
          end if
!
          do 1321 k = 1,nm(2)
              bquad(2,2,k) = aquad(2,1,k)*bquad(2,2,k)
 1321     end do
!
!---------***********
          call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev,quad, &
                      nm,score,hess,ncov,1,arith,'n',ifail,offlag, &
                      offpos,maxcol,link,risk,corr,bivar,iter,1,1d0, &
                      scnorm,scores,order,ncat,n1lev,family,trivar, &
                      n2lev,ilfit,inorm,nlevel,depend,eqscale,i, &
                      xil_tau(2),nosign,iquad,bquad,maxit,mus,taus, &
                      maxcat)
!---------***********
!
          if (ifail) then
              return
          end if
!
          mu(2,m,i) = sign(2)*exp(xil_mu(2) - xil)
          taui = exp(xil_tau(2) - xil) - mu(2,m,i)**2
!
          if (abs(taui) > dcon) then
              tau(2,m,i) = sqrt(abs(taui))
          else
              tau(2,m,i) = tau(2,m-1,i)/2
          end if
!
          mus(2,i) = mu(2,m,i)
          taus(2,i) = tau(2,m,i)
!
          if (iquad == 1 .or. iquad == 3) then
              xild = abs(xil - xill(m-1,i))
!
              if (xild < dcon) then
                  conflag(i) = .true.
!
                  if (i > 1 .or. this_processor == boss_processor) then
                      nconflag(m) = nconflag(m) + 1
                  end if
!
                  xall = xall+xil
              end if
!
              xill(m,i) = xil
          end if
!
          if (i > 1 .or. this_processor == boss_processor) then
              xal = xal+xil
          end if
!
 1501 end do
!
 1503 if (num_processors > 1 .and. parflag) then

!---------******************
          call mpi_allreduce(ifail,anyfail,1,mpi_logical,mpi_lor, &
                             mpi_comm_world,ierror)
!---------******************

      else
          anyfail = ifail
      end if
!
      if (anyfail) then
          ifail = .true.
          return
      end if
!
      if (num_processors > 1 .and. parflag) then

!---------******************
          call mpi_allreduce(xal,xal_total,1,mpi_double_precision, &
                             mpi_sum,mpi_comm_world,ierror)
!---------******************

          xal = xal_total

!---------******************
          call mpi_allreduce(nconflag(m),nconflags(m),1,mpi_integer, &
                             mpi_sum,mpi_comm_world,ierror)
!---------******************

      else
          nconflags(m) = nconflag(m)
      end if
!
      if (qflag .and. ilfit == 0) then
!
!---------**********
          call clock('',clock2,2)
!---------**********
!
      end if
!
      if (tflag .and. (m == 1 .or. iquad == 1 .or. iquad == 3)) &
      then
          write (outbuf,'(a)') clkstr(1:16)
          call wrtlit(outbuf)
          write (outbuf,'(a,f14.4)') 'log likelihood = ',xal
          call wrtlit(outbuf)
          write (outbuf,'(a,i6)') 'cases converged = ',nconflags(m)
          call wrtlit(outbuf)
          call newlit
      end if
!
      if (m == 1 .or. (m /= 1 .and. nconflags(m) < nsub(nlevel) &
      .and. (nconflags(m) > nconflags(m-1) .or. &
      nconflags(m) <= nsub(nlevel)/2))) then
!
          if (m == 100) then
              call wrtlin( &
              '    *** ERROR *** Adaptive quadrature algorithm failed')
              ifail = .true.
              return
           else
              go to 1150
          end if
!
      end if
!
!---- get log-likelihood, score vector and Hessian matrix at current
!---- estimates
!-----***********
   12 call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev,quad,nm, &
                  score,hess,ncov,flag,arith,endind,ifail,offlag,offpos, &
                  maxcol,link,risk,corr,bivar,iter,1,1d0,scnorm,scores, &
                  order,ncat,n1lev,family,trivar,n2lev,ilfit,inorm, &
                  nlevel,depend,eqscale,0,xil,nosign,iquad,aquad,maxit, &
                  mus,taus,maxcat)
!-----***********
!
      if (ifail) then
          return
      end if
!
      if (iter == 1) then
          wll = xll
      end if
!
!---- STICK is used to prevent the optimisation algorithm 'sticking'
!---- when an endpoint is at or very close to its boundary zero, with a
!---- positive log-likelihood derivative, but a negative N-R increment.
!---- When these conditions occur, STICK is set to TRUE and the
!---- iteration is repeated with the endpoint held at the boundary
      stick = .false.
!---- set OLL to be log-likelihood at start of iteration
      oll = xll
!
!---- this is where you come back to if an endpoint has stuck
!---- i.e. you do the whole thing again with the endpoint held fixed
!---- store endpoint probabilities in default locations indexed by LOC
   20 do 30 i = 1,nest
          loc(i) = i
   30 end do
!
!---- we now sort out how endpoints are to be treated in this iteration
!---- STICK over-rides the rules for releasing endpoints
!---- algorithm stuck at endpoint
      if (stick) then
          go to 35
      end if
!
!---- release fixed endpoint if log-likelihood derivative positive
      if (nl == 1 .and. score(nest-nr) > 0) then
          iend(1) = 1
      end if
!
      if (nr == 1 .and. score(nest) > 0) then
          iend(2) = 1
      end if
!
!---- if the left endpoint is fixed and the right free, scores and
!---- Hessian elements for the right are shuffled forward in the
!---- appropriate vectors with new locations indexed by LOC(.)
   35 if (nl+nr == 2 .and. iend(1) == 0 .and. iend(2) == 1) then
          score(nest-1) = score(nest)
          loc(nest-1) = nest
          l1 = nest*(nest-1)/2
!
          do 40 j = 1,nest-2
              hess(l1+j+1-nest) = hess(l1+j)
   40     end do
!
          hess(l1) = hess(ncov)
!
      end if
!
!---- number of free parameters
      free = nest - nl - nr + iend(1) + iend(2)
!---- default settings for the adjustments to Hessian diagonal elements
      adjust = 1
      doadd = .false.
      value = -1
   45 orthog = .false.
!
!---- get covariance matrix
!-----***********
      call covmat(hess,cov,alias,adjust,free,work,free*(free+1)/2,doadd, &
                  value,tol,ifail,arith,nest)
!-----***********
!
!---- either Hessian matrix not positive semi-definite or orthogonality
   47 if (ifail .or. orthog) then
!
          if (ifail .and. flag == 3) then
              flag = 2
              appflag = .true.
              go to 12
          else if (orthog) then
              go to 50
          end if
!
!         code up to line 50 not used?
          doadd = .false.
!
          if (wrnflg .and. nl+nr /= 0) then
              write (outbuf,'(a,i7,a,g20.8,a,f7.4,5a)') '  ',iter, &
              '    ',wll,' ',step,'     ',text(iend(1)+1),' ', &
              text(iend(2)+1),'      not psd'
              call wrtlin(outbuf)
              call wterse(outbuf)
          else if (wrnflg) then
              write (outbuf,'(a,i7,a,g20.8,a,f7.4,a,g13.6,a)') &
              '  ',iter,'    ',wll,' ',step,'     ',scnorm, &
              '      not psd'
              call wrtlin(outbuf)
              call wterse(outbuf)
          end if
!
          wrnflg = .false.
   50     i = 0
          ival = 0
          value = -1
!
   55     if (i == free*(free+1)/2) then
!
!------------ update scaling factor for diagonal elements of Hessian
              if (.not. doadd) then
                  adjust = 2*adjust
              end if
!
              go to 45
          end if
!
          ival = ival+1
          i = i+ival
!
          if (hess(i) <= 0 .and. -hess(i) > value) then
!------------ increment for diagonal elements of Hessian
              value = 1 - hess(i)
!------------ negative diagonal element found
              doadd = .true.
          end if
!
          go to 55
      end if
!
      k = 0
!
      do 60 i = 1,free
!
          do 58 j = 1,i
              k = k+1
!------------ store the elements of the symmetric covariance matrix
              covari(i,j) = cov(k)
!
              if (i /= j) then
                  covari(j,i) = cov(k)
              end if
!
   58     end do
!
   60 end do
!
      delta(nest-1) = 0
      delta(nest) = 0
!
!---- calculate delta = [d_i] = H`g, where H` = [h_ij] is the inverse
!---- Hessian and g = [g_j] is the score vector
      do 64 i = 1,free
          delta(loc(i)) = 0
!
!-------- calculate sum{j=1,...,N}h_ij.g_j
          do 62 j = 1,free
              delta(loc(i)) = delta(loc(i)) + covari(i,j)*score(j)
   62     end do
!
   64 end do
!
      sing = 0
      bottom = 0
!
      do 65 i = 1,free
!-------- update sum{i=1,...,N}d_i.g_i
          sing = sing + delta(loc(i))*score(i)
!-------- update sum{i=1,...,N}d_i^2
          bottom = bottom + delta(loc(i))*delta(loc(i))
   65 end do
!
!---- calculate orthogonality measure g'delta/delta'delta
      sing = sing/bottom
!
      if (sing < alph) then
          orthog = .true.
!---- reset the orthogonality criterion to its original value
      else
          alph = alp
      end if
!
      if (appflag) then
          text1 = '*   '
      else
          text1 = '    '
      end if
!
      if (wrnflg .and. .not. orthog .and. nl+nr /= 0) then
          write (outbuf,'(a,i7,a,g20.8,a,f7.4,5a,g13.6)') &
          '  ',iter,text1,wll,' ',step,'     ',text(iend(1)+1),' ', &
          text(iend(2)+1),'     ',sing
          call wrtlin(outbuf)
          call wterse(outbuf)
      else if (wrnflg .and. nl+nr /= 0) then
          write (outbuf,'(a,i7,a,g20.8,a,f7.4,5a)') &
          '  ',iter,text1,wll,' ',step,'     ',text(iend(1)+1),' ', &
          text(iend(2)+1),'      orthogonal'
          call wrtlin(outbuf)
          call wterse(outbuf)
      else if (wrnflg .and. .not. orthog) then
          write (outbuf,'(a,i7,a,g20.8,a,f7.4,a,g13.6,a,g13.6)') &
          '  ',iter,text1,wll,' ',step,'     ',scnorm,'   ',sing
          call wrtlin(outbuf)
          call wterse(outbuf)
      else if (wrnflg) then
          write (outbuf,'(a,i7,a,g20.8,a,f7.4,a,g13.6,a)') &
          '  ',iter,text1,wll,' ',step,'     ',scnorm, &
          '    orthogonal'
          call wrtlin(outbuf)
          call wterse(outbuf)
      end if
!
      wrnflg = .false.
!
      if (orthog) then
          go to 47
      end if
!
      do 67 i = 1,free
!
!-------- set delta to zero for extrinsically aliased parameters
          if (idnint(alias(i)) == 1) then
              delta(i) = 0
          end if
!
   67 end do
!
!---- Re-initialise STICK to default value
      stick = .false.
!---- CHOP has default value false and is set to true whenever the lower
!---- bound of zero for the endpoint parameters restricts the step
!---- length at an iteration.  This avoids a spurious convergence due to
!---- severe restrictions on step length. This pragmatic procedure may
!---- benefit from further thought.
      chop = .false.
!---- we now ensure that the initial step length does not carry end
!---- point probabilities beyond lower bound of zero
!---- endpoint(s) fixed
      if (free == nest-nl-nr) then
          step = 2
!---- at least one free endpoint in model
      else
          zb = 1
!
          do 70 i = nest-nl-nr+1,nest
!
!------------ endpoint fixed or aliased
              if (iend(i+1-nest+nr) == 0 .or. idnint(alias(i)) == 1) &
              then
                  go to 70
              end if
!
              za = 1
!------------ update endpoint probability
              ab = beta(i) + delta(i)
!
!------------ decrease step-length to prevent endpoint probability
!------------ becoming negative
              if (ab < 0) then
                  za = -beta(i)/delta(i)
              end if
!
!------------ endpoint probability close to zero with negative
!------------ increment, so fix endpoint and indicate that algorithm
!------------ has stuck
              if (abs(beta(i)) < 1.0d-9 .and. delta(i) < 0) then
                  iend(i+1-nest+nr) = 0
                  stick = .true.
              end if
!
!------------ maximum step-length to prevent endpoint probability
!------------ becoming negative
              if (za < zb) then
                  zb = za
!---------------- set CHOP to TRUE if zero bounds on endpoint
!---------------- parameters have restricted step length
                  chop = .true.
              end if
!
   70     end do
!
          step = 2*zb
      end if
!
!---- algorithm stuck
      if (stick) then
          oldadj = adjust
          oldadd = doadd
          go to 20
      end if
!
      if (qflag .and. ilfit == 0) then
          clkstr = 'FIT iteration '
          write (itstr,'(i3)') iter
          clkstr(15:17) = itstr
          write (outbuf,'(a)') clkstr(1:17)
          call wrtflq(outbuf)
!
!---------**********
          call clock(clkstr(1:17),clock1,2)
!---------**********
!
      end if
!
!---- update iteration number
      iter = iter+1
!
      if (qflag .and. ilfit == 0) then
          clkstr = 'FIT iteration '
          write (itstr,'(i3)') iter
          clkstr(15:17) = itstr
!
!---------**********
          call clock(clkstr(1:17),clock1,1)
!---------**********
!
          call newliq
      end if
!
!---- no more iterations allowed
      if (iter > niter) then
          call wrtlin('    *** WARNING *** The maximum number of '// &
                      'iterations has been reached')
          call newlin
          return
      end if
!
      dflag = .false.
!
!---- now move into crude step length search, halving step length if
!---- necessary until improvement obtained in log-likelihood or
!---- convergence identified
   80 step = step/2
      wrnflg = .true.
!
      if (step < 1.0d-4) then
          call wrtlin('    *** ALGORITHM FAILED, step length too small')
          call wterse('    *** ALGORITHM FAILED, step length too small')
          ifail = .true.
          return
      end if
!
      if (bivar .and. corr) then
!
          if (.not. dflag) then
              betasig(2) = beta(nest-1)
              betarho(1) = beta(nest)
          end if
!
          beta(nest-1) = betasig(2)*sqrt(1 - betarho(1)**2)
          beta(nest) = betasig(2)*betarho(1)
      end if
!
      if (trivar .and. corr) then
          betasig(2) = beta(nest-4)
          betasig(3) = beta(nest-3)
          betarho(1) = beta(nest-2)
          betarho(2) = beta(nest-1)
          betarho(3) = beta(nest)
          beta(nest-4) = betasig(2)*sqrt(1 - betarho(1)**2)
          beta(nest-3) = -betasig(3)*sqrt((1 - betarho(1)**2 - &
          betarho(2)**2 - betarho(3)**2 + &
          2*betarho(1)*betarho(2)*betarho(3))/(1 - betarho(1)**2))
          beta(nest-2) = betasig(2)*betarho(1)
          beta(nest-1) = betasig(3)*betarho(2)
          beta(nest) = betasig(3)*(betarho(3) - betarho(1)*betarho(2))/ &
          sqrt(1 - betarho(1)**2)
      end if
!
      do 90 i = 1,nest
!-------- update parameter estimates, beta' = beta + step.delta
          alpha(i) = beta(i) + step*delta(i)
          dflag = .false.
!
!-------- ensure that scale parameter estimates remain strictly positive
          if (univar) then
!
              if ((i == nest-nl-nr .or. &
              (nlevel == 2 .and. i == nest-1)) .and. &
              ilfit == 0 .and. alpha(i) <= 0) then
                  dflag = .true.
              end if
!
              if (depend .and. i == nest-1 .and. ilfit == 0 .and. &
              alpha(i) <= 0) then
                  dflag = .true.
              end if
!
              if (family(1) == 'g' .and. i == nest-isca .and. &
              alpha(i) <= 0) then
                  dflag = .true.
              end if
!
          else if (bivar) then
!
              if ((corr .and. .not. eqscale) .or. .not. corr) then
!
                  if ((i == nest-icorr-1 .or. i == nest-icorr) .and. &
                  alpha(i) <= 0 .and. ilfit /= 2) then
                      dflag = .true.
                  end if
!
              end if
!
              if (eqscale .and. i == nest-1 .and. alpha(i) <= 0) then
                  dflag = .true.
              end if
!
              if (inorm == 1 .and. i == nest-icorr-isca .and. &
              alpha(i) <= 0) then
                  dflag = .true.
              end if
!
              if (inorm == 2 .and. (i == nest-icorr-isca-1 .or. &
              i == nest-icorr-isca) .and. alpha(i) <= 0) then
                  dflag = .true.
              end if
!
          else if (trivar) then
!
              if (.not. corr) then
!
                  if (i >= nest-2 .and. alpha(i) <= 0 .and. ilfit /= 2) &
                  then
                      dflag = .true.
                  end if
!
              else
!
                  if (((i == nest-5 .or. i == nest-4) .and. &
                  alpha(i) <= 0) .or. (i == nest-3 .and. &
                  alpha(i) >= 0)) then
                      dflag = .true.
                  end if
!
              end if
!
              if (inorm == 1 .and. i == nest-icorr-3 .and. &
              alpha(i) <= 0) then
                  dflag = .true.
              end if
!
              if (inorm == 2 .and. (i == nest-icorr-4 .or. &
              i == nest-icorr-3) .and. alpha(i) <= 0) then
                  dflag = .true.
              end if
!
              if (inorm == 3 .and. i >= nest-icorr-5 .and. &
              i <= nest-icorr-3 .and. alpha(i) <= 0) then
                  dflag = .true.
              end if
!
          end if
!
          if (dflag) then
              go to 80
          end if
!
   90 end do
!
      if (order) then
!
!---------*********
          call sort(alpha(nilev+1),ncat(1)-1)
!---------*********
!
          if (.not. univar) then
!
!-------------*********
              call sort(alpha(nilev+ncat(1)),ncat(2)-1)
!-------------*********
!
              if (trivar) then
!
!-----------------*********
                  call sort(alpha(nilev+ncat(1)+ncat(2)-1),ncat(3)-1)
!-----------------*********
!
              end if
!
          end if
!
      end if
!
      if (bivar .and. corr .and. ilfit /= 2) then
          betaphi(2) = alpha(nest-1)
          betaphi(3) = alpha(nest)
          alpha(nest-1) = sqrt(betaphi(2)**2 + betaphi(3)**2)
          alpha(nest) = betaphi(3)/alpha(nest-1)
!
          if (eqscale .and. (alpha(nest) <= (1d-6)-1 .or. &
          alpha(nest) >= 1-(1d-6))) then
              call wrtlin('    *** ALGORITHM FAILED, rho at limit')
              call wterse('    *** ALGORITHM FAILED, rho at limit')
              ifail = .true.
              return
          end if
!
      end if
!
      if (trivar .and. corr) then
          betaphi(2) = alpha(nest-4)
          betaphi(3) = alpha(nest-3)
          betaphi(4) = alpha(nest-2)
          betaphi(5) = alpha(nest-1)
          betaphi(6) = alpha(nest)
          alpha(nest-4) = sqrt(betaphi(2)**2 + betaphi(4)**2)
          alpha(nest-3) = &
          sqrt(betaphi(3)**2 + betaphi(5)**2 + betaphi(6)**2)
          alpha(nest-2) = betaphi(4)/alpha(nest-4)
          alpha(nest-1) = betaphi(5)/alpha(nest-3)
          alpha(nest) = (betaphi(2)*betaphi(6) + betaphi(4)*betaphi(5))/ &
          (alpha(nest-4)*alpha(nest-3))
      end if
!
!---- get log-likelihood and score vector at updated estimates
!-----***********
      call lshess(x,y,nmes,xll,produc,it,nsub,alpha,nest,nilev,quad,nm, &
                  score,hess,ncov,2,arith,endind,ifail,offlag,offpos, &
                  maxcol,link,risk,corr,bivar,iter,2,step,scnorm,scores, &
                  order,ncat,n1lev,family,trivar,n2lev,ilfit,inorm, &
                  nlevel,depend,eqscale,0,xil,nosign,iquad,aquad,maxit, &
                  mus,taus,maxcat)
!-----***********
!
      if (ifail) then
          return
      end if
!
      wll = xll
!
!---- change in log-likelihood between current and previous iterations
      diff = xll-oll
      scomax = 0
!
      do 95 i = 1,nest-nl-nr
!
          if (abs(score(i)) > scomax) then
              scomax = abs(score(i))
          end if
!
   95 end do
!
!---- convergence maybe?
      if (abs(diff) < dcon .and. .not. chop) then
!      IF ((diff >= 0 .or. halfstep) .and. abs(DIFF) <= CON .AND. &
!      .NOT. CHOP) THEN
          recycl = .true.
!
!-------- halve step-length
          if (.not. halfstep) then
              step = step/2
          end if
!
!-------- update parameter estimates, beta' = beta + step.delta
          do 100 i = 1,nest
              gamma(i) = beta(i) + step*delta(i)
  100     end do
!
          if (bivar .and. corr .and. ilfit /= 2) then
              betaphi(2) = gamma(nest-1)
              betaphi(3) = gamma(nest)
              gamma(nest-1) = sqrt(betaphi(2)**2 + betaphi(3)**2)
              gamma(nest) = betaphi(3)/gamma(nest-1)
          end if
!
          if (trivar .and. corr) then
              betaphi(2) = gamma(nest-4)
              betaphi(3) = gamma(nest-3)
              betaphi(4) = gamma(nest-2)
              betaphi(5) = gamma(nest-1)
              betaphi(6) = gamma(nest)
              gamma(nest-4) = sqrt(betaphi(2)**2 + betaphi(4)**2)
              gamma(nest-3) = &
              sqrt(betaphi(3)**2 + betaphi(5)**2 + betaphi(6)**2)
              gamma(nest-2) = betaphi(4)/gamma(nest-4)
              gamma(nest-1) = betaphi(5)/gamma(nest-3)
              gamma(nest) = (betaphi(2)*betaphi(6) + &
              betaphi(4)*betaphi(5))/(gamma(nest-4)*gamma(nest-3))
          end if
!
!-------- get log-likelihood at updated estimates
!---------***********
          call lshess(x,y,nmes,xxll,produc,it,nsub,gamma,nest,nilev, &
                     quad,nm,score,hess,ncov,1,arith,endind,ifail, &
                     offlag,offpos,maxcol,link,risk,corr,bivar,iter,3, &
                     step,scnorm,scores,order,ncat,n1lev,family,trivar, &
                     n2lev,ilfit,inorm,nlevel,depend,eqscale,0,xil, &
                     nosign,iquad,aquad,maxit,mus,taus,maxcat)
!---------***********
!
          if (ifail) then
              return
          end if
!
!-------- log-likelihood change for 'half-step' from previous iteration
          diff = xxll-oll
!
          if (diff <= dcon .and. scnorm < dcon) then
              done = .true.
          end if
!
!---- decrease in log-likelihood from previous iteration, so 'half-step'
      else if (abs(diff) > dcon .and. diff < dcon) then
!      ELSE IF (DIFF < CON) THEN
          halfstep = .true.
!
          if (bivar .and. corr .and. ilfit /= 2) then
              beta(nest-1) = betasig(2)
              beta(nest) = betarho(1)
          end if
!
          if (trivar .and. corr) then
              beta(nest-4) = betasig(2)
              beta(nest-3) = betasig(3)
              beta(nest-2) = betarho(1)
              beta(nest-1) = betarho(2)
              beta(nest) = betarho(3)
          end if
!
          go to 80
      end if
!
      if (recycl .and. xxll > xll) then
          wll = xxll
!
!-------- set parameter estimates to those from recycling iteration
          do 110 i = 1,nest
              beta(i) = gamma(i)
  110     end do
!
      else
!
!-------- set parameter estimates to those from pre-recycling iteration
          do 120 i = 1,nest
              beta(i) = alpha(i)
  120     end do
!
!-------- reset the steplength to its value prior to recycling
          if (recycl) then
              step = 2*step
          end if
!
      end if
!
!---- Specify fixed/free status of endpoint parameters for next
!---- iteration
      if (nl+nr /= 0) then
!
          do 130 i = nest-nl-nr+1,nest
!
              if (abs(beta(i)) < 1.0d-10) then
                  beta(i) = 0
!
!---------------- fix endpoint(s)
                  if (.not. done) then
                      iend(i+1-nest+nr) = 0
                  end if
!
              end if
!
  130     end do
!
      end if
!
!---- not yet converged
      if (.not. done) then
          go to 10
      end if
!
!---- Convergence achieved. Final checks to make sure Hessian was not
!---- adjusted. If it was, set copy of alpha to be a very small number
!---- and jump back in.
      if (doadd .or. oldadd .or. adjust > 1 .or. oldadj > 1) then
!-------- reset orthogonality criterion (should prevent orthogonality)
          alph = 1.0d-9
!-------- converged maybe? - no way!
          oldadd = .false.
          oldadj = 1
!-------- another iteration needed
          go to 10
      end if
!
      if (appflag) then
          text1 = '*   '
      else
          text1 = '    '
      end if
!
      if (nl+nr /= 0) then
          write (outbuf,'(a,i7,a,g20.8,a,f7.4,4a)') &
          '  ',iter,text1,xll,' ',step,'     ',text(iend(1)+1),' ', &
          text(iend(2)+1)
          call wrtlin(outbuf)
          call wterse(outbuf)
          call newlin
      else
          write (outbuf,'(a,i7,a,g20.8,a,f7.4,a,g13.6)') &
          '  ',iter,text1,xll,' ',step,'     ',scnorm
          call wrtlin(outbuf)
          call wterse(outbuf)
          call newlin
      end if
!
!      IF ((STEP < 1.0D-4 .OR. RECYCL) .AND. SCOMAX > sqrt(CON))
!     &THEN
!          call wrtlin('    *** WARNING *** '//
!     &                'The algorithm has failed to converge')
!          call newlin
!          call wrtlin('    The score vector is:')
!          nlines = ( nest-nl-nr+5 ) / 6
!
!          do 140 iline = 1, nlines - 1
!              write (outbuf,'(4x,6f12.3)') &
!              (score(i),i=6*iline-5,6*iline)
!              call wrtlin(outbuf)
!140       end do
!
!          write (outbuf,'(4x,6f12.3)')
!     &    (score(i),i=6*nlines-5,nest-nl-nr)
!          call wrtlin(outbuf)
!          call newlin
!          write (outbuf,'(a,f12.3)')
!     &    '    The largest element of the score vector is ',scomax
!          call wrtlin(outbuf)
!          write (outbuf,'(a,f8.6)') '    For convergence, all '//
!     &    'elements must be less than or equal to ',sqrt(con)
!          call wrtlin(outbuf)
!          call newlin
!
!          RETURN
!      END IF
!
      if (dfirst) then
          flag = 2
      else
          flag = 3
      end if
!
!---- get log-likelihood, score vector and true Hessian matrix at final
!---- parameter estimates
!-----***********
      call lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev,quad,nm, &
                  score,hess,ncov,flag,arith,endind,ifail,offlag,offpos, &
                  maxcol,link,risk,corr,bivar,iter,4,1d0,scnorm,scores, &
                  order,ncat,n1lev,family,trivar,n2lev,ilfit,inorm, &
                  nlevel,depend,eqscale,0,xil,nosign,iquad,aquad,maxit, &
                  mus,taus,maxcat)
!-----***********
!
      if (ifail) then
          return
      end if
!
!---- if the left endpoint is fixed and the right free, Hessian
!---- elements for the right are shuffled forward in the vector
      if (nl+nr == 2 .and. iend(1) == 0 .and. iend(2) == 1) then
!
          do 216 j = 1,nest-2
              hess(l1+j+1-nest) = hess(l1+j)
  216     end do
!
          hess(l1) = hess(ncov)
      end if
!
!---- get covariance matrix
!-----***********
      call covmat(hess,cov,alias,1d0,free,work,free*(free+1)/2,.false., &
                  value,tol,ifail,arith,nest)
!-----***********
!
      if (ifail) then
          call wrtlin('    *** WARNING *** The final Hessian matrix '// &
                      'is not positive semi-definite')
          call newlin
          ifail = .false.
          return
      end if
!
      if (bivar .and. corr .and. ilfit /= 2) then
          phi(1) = beta(nest-2)
          phi(2) = betaphi(2)
          phi(3) = betaphi(3)
          varphi(2) = cov((nest-1)*nest/2)
          varphi(3) = cov(nest*(nest+1)/2)
          covphi(2,3) = cov(nest*(nest+1)/2 - 1)
!
!-------- var(sigma2)
          cov((nest-1)*nest/2) = abs((phi(2)**2*varphi(2) + &
          phi(3)**2*varphi(3) + 2*phi(2)*phi(3)*covphi(2,3))/ &
          beta(nest-1)**2)
!
!-------- var(rho)
          cov(nest*(nest+1)/2) = abs((phi(3)**2*varphi(2) + &
          phi(2)**2*varphi(3) - 2*phi(2)*phi(3)*covphi(2,3))* &
          (1 - beta(nest)**2)/beta(nest-1)**4)
!
!-------- cov12, se(cov12)
!          sig12 = phi(1)*phi(3)
!          sesig12 = sqrt(phi(1)**2*varphi(3) + phi(3)**2*varphi(1)
!     &    + 2*phi(1)*phi(3)*covphi(1,3))
      end if
!
      if (trivar .and. corr) then
          phi(1) = beta(nest-5)
          phi(2) = betaphi(2)
          phi(3) = betaphi(3)
          phi(4) = betaphi(4)
          phi(5) = betaphi(5)
          phi(6) = betaphi(6)
          varphi(1) = cov((nest-5)*(nest-4)/2)
          varphi(2) = cov((nest-4)*(nest-3)/2)
          varphi(3) = cov((nest-3)*(nest-2)/2)
          varphi(4) = cov((nest-2)*(nest-1)/2)
          varphi(5) = cov((nest-1)*nest/2)
          varphi(6) = cov(nest*(nest+1)/2)
          covphi(1,4) = cov((nest-2)*(nest-1)/2 - 3)
          covphi(1,5) = cov((nest-1)*nest/2 - 4)
          covphi(2,3) = cov((nest-3)*(nest-2)/2 - 1)
          covphi(2,4) = cov((nest-2)*(nest-1)/2 - 2)
          covphi(2,5) = cov((nest-1)*nest/2 - 3)
          covphi(2,6) = cov(nest*(nest+1)/2 - 4)
          covphi(3,4) = cov((nest-2)*(nest-1)/2 - 1)
          covphi(3,5) = cov((nest-1)*nest/2 - 2)
          covphi(3,6) = cov(nest*(nest+1)/2 - 3)
          covphi(4,5) = cov((nest-1)*nest/2 - 1)
          covphi(4,6) = cov(nest*(nest+1)/2 - 2)
          covphi(5,6) = cov(nest*(nest+1)/2 - 1)
          v1 = beta(nest-4)**2
          v2 = beta(nest-3)**2
          v3 = phi(2)*phi(6) + phi(4)*phi(5)
          v4 = phi(2)*phi(5) - phi(4)*phi(6)
          v5 = phi(2)*(phi(3)**2 + phi(5)**2) - phi(4)*phi(5)*phi(6)
          v6 = phi(4)*(phi(3)**2 + phi(6)**2) - phi(2)*phi(5)*phi(6)
!
!-------- var(sigma1)
!
!-------- var(sigma2)
          cov((nest-4)*(nest-3)/2) = abs((phi(2)**2*varphi(2) + &
          phi(4)**2*varphi(4) + 2*phi(2)*phi(4)*covphi(2,4))/v1)
!
!-------- var(sigma3)
          cov((nest-3)*(nest-2)/2) = abs((phi(3)**2*varphi(3) + &
          phi(5)**2*varphi(5) + phi(6)**2*varphi(6) + &
          2*(phi(3)*phi(5)*covphi(3,5) + phi(3)*phi(6)*covphi(3,6) + &
          phi(5)*phi(6)*covphi(5,6)))/v2)
!
!-------- cov12, se(cov12)
          sig12 = phi(1)*phi(4)
          sesig12 = sqrt(phi(4)**2*varphi(1) + phi(1)**2*varphi(4) + &
          2*phi(1)*phi(4)*covphi(1,4))
!
!-------- cov13, se(cov13)
          sig13 = phi(1)*phi(5)
          sesig13 = sqrt(phi(5)**2*varphi(1) + phi(1)**2*varphi(5) + &
          2*phi(1)*phi(5)*covphi(1,5))
!
!-------- cov23, se(cov23)
          sig23 = v3
          sesig23 = sqrt(phi(6)**2*varphi(2) + phi(5)**2*varphi(4) + &
          phi(4)**2*varphi(5) + phi(2)**2*varphi(6) + &
          2*(phi(5)*phi(6)*covphi(2,4) + phi(4)*phi(6)*covphi(2,5) + &
          phi(2)*phi(6)*covphi(2,6) + phi(4)*phi(5)*covphi(4,5) + &
          phi(2)*phi(5)*covphi(4,6) + phi(2)*phi(4)*covphi(5,6)))
!
!-------- var(rho12)
          cov((nest-2)*(nest-1)/2) = abs(phi(2)**2*(phi(4)**2*varphi(2) &
          + phi(2)**2*varphi(4) - 2*phi(2)*phi(4)*covphi(2,4))/(v1**3))
!
!-------- var(rho13)
          cov((nest-1)*nest/2) = abs(((phi(3)*phi(5))**2*varphi(3) + &
          (phi(3)**2 + phi(6)**2)**2*varphi(5) + &
          (phi(5)*phi(6))**2*varphi(6) - &
          2*phi(5)*(phi(3)*(phi(3)**2 + phi(6)**2)*covphi(3,5) - &
          phi(3)*phi(5)*phi(6)*covphi(3,6) + &
          phi(6)*(phi(3)**2 + phi(6)**2)*covphi(5,6)))/(v2**3))
!
!-------- var(rho23)
          cov(nest*(nest+1)/2) = abs(((phi(4)*v4/v1)**2*varphi(2) + &
          (phi(3)*v3/v2)**2*varphi(3) + (phi(2)*v4/v1)**2*varphi(4) + &
          (v6/v2)**2*varphi(5) + (v5/v2)**2*varphi(6) + &
          2*(phi(3)*phi(4)*(v3*v4/(v1*v2))*covphi(2,3) - &
          phi(2)*phi(4)*(v4/v1)**2*covphi(2,4) - &
          phi(4)*(v4*v6/(v1*v2))*covphi(2,5) - &
          phi(4)*(v4*v5/(v1*v2))*covphi(2,6) - &
          phi(2)*phi(3)*(v3*v4/(v1*v2))*covphi(3,4) - &
          phi(3)*(v3*v6/(v2**2))*covphi(3,5) - &
          phi(3)*(v3*v5/(v2**2))*covphi(3,6) + &
          phi(2)*(v4*v6/(v1*v2))*covphi(4,5) + &
          phi(2)*(v4*v5/(v1*v2))*covphi(4,6) + &
          (v5*v6/(v2**2))*covphi(5,6)))/(v1*v2))
      end if
!
      k = 0
!
      do 11160 i = 1,nest
!
          do 11158 j = 1,i
              k = k+1
!------------ store the elements of the symmetric covariance matrix
              covari(i,j) = cov(k)
!
              if (i /= j) then
                  covari(j,i) = cov(k)
              end if
!
11158     end do
!
11160 end do
!
!      do 11170 i = 1,nest
!
!          do 11165 j = 1,nest
!              write (1,*) 'cov(',i,',',j,') = ',covari(i,j)
!11165     end do
!
!11170 end do
!
      j = 0
!
!---- left endpoint fixed, right endpoint free
      if (nl+nr == 2 .and. iend(1) == 0 .and. iend(2) == 1) then
!-------- switch aliasing indicators
          temp = alias(nest-1)
          alias(nest-1) = alias(nest)
          alias(nest) = temp
!
!-------- correctly align entries in covariance matrix so that
!-------- non-negative values refer to right endpoint
          do 220 j = 1,nest-2
              cov(l1+j) = cov(l1+j+1-nest)
              cov(l1+j+1-nest) = 0
  220     end do
!
          cov(ncov-1) = 0
          cov(ncov) = cov(l1)
          cov(l1) = 0
      end if
!
!      if (robust) then
!
!---------***********
!          CALL ROBFIT(scores,nmes,nest,cov,ncov,covrob)
!---------***********
!
!          do 300 i = 1,ncov
!              cov(i) = covrob(i)
!  300     end do
!
!      end if
!
      if (qflag .and. ilfit == 0) then
          clkstr = 'FIT iteration '
          write (itstr,'(i3)') iter
          clkstr(15:17) = itstr
          write (outbuf,'(a)') clkstr(1:17)
          call wrtflq(outbuf)
!
!---------**********
          call clock(clkstr(1:17),clock1,2)
!---------**********
!
      end if
!
      return
!
      end subroutine fit
!
!***********************************************************************
!
      subroutine lshess(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev, &
                        quad,nm,score,hess,ncov,flag,arith,endind,ifail, &
                        offlag,offpos,maxcol,link,risk,corr,bivar,iter, &
                        cflag,step,scnorm,scores,order,ncat,n1lev, &
                        family,trivar,n2lev,ilfit,inorm,nlevel,depend, &
                        eqscale,ind,xil,sign,iquad,aquad,maxit,mus,taus, &
                        maxcat)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      include 'limitsf90.h'
!-----------------------------------------------------------------------
      character arith,endind,link(3),family(3), &
                itstr*3,clkstr*40,derstr*17
      integer nsub(2),it(2,nsub(1)),nm(3),flag,offpos(3),n1lev, &
              produc(nsub(1),2),maxcol,risk(nmes),iter,cflag,ncat(3), &
              n2lev,ilfit,inorm,clock1,nmes,nest,nilev,ncov,nlevel,ind, &
              sign,iquad,maxit,maxcat,onm(3)
      double precision y(nmes),x(nmes,nilev),beta(nest),quad(3,2,256), &
             score(nest),xll,hess(ncov),step,scnorm,scores(nmes,nest), &
             xil,aquad(3,2,256),mus(3,nsub(1)),taus(3,nsub(1)), &
             oquad(3,2,256)
      logical ifail,offlag(3),corr,bivar,order,trivar,depend,eqscale
!-----------------------------------------------------------------------
!     Function : Calls a likelihood routine based upon arithmetic type.
!-----------------------------------------------------------------------
      include 'chans.h'
!-----------------------------------------------------------------------
      if (qflag .and. ind == 0 .and. ilfit == 0) then
          clkstr = 'FIT iteration '
          write (itstr,'(i3)') iter
          clkstr(15:17) = itstr
!
          if (flag == 1) then
              derstr = ', no derivatives'
          else if (flag == 2) then
              derstr = ', 1st derivatives'
          else
              derstr = ', 2nd derivatives'
          end if
!
          clkstr(18:34) = derstr
!
!---------**********
          call clock(clkstr(1:34),clock1,1)
!---------**********
!
      end if
!
      ifail = .false.
!
      if (trivar .and. order) then
!
          if (ilfit == 2) then
              onm(1) = nm(1)
              onm(2) = nm(2)
              onm(3) = nm(3)
              oquad(1,1,1) = quad(1,1,1)
              oquad(1,2,1) = quad(1,2,1)
              oquad(2,1,1) = quad(2,1,1)
              oquad(2,2,1) = quad(2,2,1)
              oquad(3,1,1) = quad(3,1,1)
              oquad(3,2,1) = quad(3,2,1)
              nm(1) = 1
              nm(2) = 1
              nm(3) = 1
              quad(1,1,1) = 0
              quad(1,2,1) = 1
              quad(2,1,1) = 0
              quad(2,2,1) = 1
              quad(3,1,1) = 0
              quad(3,2,1) = 1
          end if
!
          if (arith == 'f') then
!
!-------------***********
              call lstore(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                          score,hess,ncov,flag,ifail,maxcol,link,risk, &
                          corr,iter,cflag,step,scnorm,scores,n1lev, &
                          n2lev,family,ind,xil,sign,iquad,aquad,mus, &
                          taus,offlag,offpos,ncat,maxcat)
!-------------***********
!
          else
!
!-------------***************
              call lstore_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                              nm,score,hess,ncov,flag,ifail,maxcol,link, &
                              risk,corr,iter,cflag,step,scnorm,scores, &
                              n1lev,n2lev,family,ind,xil,sign,iquad, &
                              aquad,mus,taus,offlag,offpos,ncat,maxcat)
!-------------***************
!
          end if
!
          if (ilfit == 2) then
              nm(1) = onm(1)
              nm(2) = onm(2)
              nm(3) = onm(3)
              quad(1,1,1) = oquad(1,1,1)
              quad(1,2,1) = oquad(1,2,1)
              quad(2,1,1) = oquad(2,1,1)
              quad(2,2,1) = oquad(2,2,1)
              quad(3,1,1) = oquad(3,1,1)
              quad(3,2,1) = oquad(3,2,1)
          end if
!
      else if (trivar) then
!
!-------- get log-likelihood (, score vector and Hessian matrix) at
!-------- current estimates
          if (arith == 'f') then
!
!-------------***********
!          write (*,'(a)') 'calling lshtri'
              call lshtri(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                          score,hess,ncov,flag,ifail,maxcol,link,risk, &
                          corr,iter,cflag,step,scnorm,scores,n1lev, &
                          family,n2lev,inorm,ind,xil,sign,iquad,aquad, &
                          mus,taus,offlag,offpos)
!-------------***********
!
          else
!
!-------------***************
!          write (*,'(a)') 'calling lshess_acc'
              call lshtri_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                              nm,score,hess,ncov,flag,ifail,maxcol,link, &
                              risk,corr,iter,cflag,step,scnorm,scores, &
                              n1lev,family,n2lev,inorm,ind,xil,sign, &
                              iquad,aquad,mus,taus,offlag,offpos)
!-------------***************
!
          end if
!
      else if (bivar .and. order) then
!
          if (ilfit == 2) then
              onm(1) = nm(1)
              onm(2) = nm(2)
              oquad(1,1,1) = quad(1,1,1)
              oquad(1,2,1) = quad(1,2,1)
              oquad(2,1,1) = quad(2,1,1)
              oquad(2,2,1) = quad(2,2,1)
              nm(1) = 1
              nm(2) = 1
              quad(1,1,1) = 0
              quad(1,2,1) = 1
              quad(2,1,1) = 0
              quad(2,2,1) = 1
          end if
!
          if (arith == 'f') then
!
!-------------***********
              call lsbore(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                          score,hess,ncov,flag,ifail,maxcol,link,risk, &
                          corr,iter,cflag,step,scnorm,scores,n1lev, &
                          family,ind,xil,sign,iquad,aquad,mus,taus, &
                          offlag,offpos,ncat,maxcat)
!-------------***********
!
          else
!
!-------------***************
              call lsbore_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                              nm,score,hess,ncov,flag,ifail,maxcol,link, &
                              risk,corr,iter,cflag,step,scnorm,scores, &
                              n1lev,family,ind,xil,sign,iquad,aquad,mus, &
                              taus,offlag,offpos,ncat,maxcat)
!-------------***************
!
          end if
!
          if (ilfit == 2) then
              nm(1) = onm(1)
              nm(2) = onm(2)
              quad(1,1,1) = oquad(1,1,1)
              quad(1,2,1) = oquad(1,2,1)
              quad(2,1,1) = oquad(2,1,1)
              quad(2,2,1) = oquad(2,2,1)
          end if
!
      else if (bivar) then
!
!-------- get log-likelihood (, score vector and Hessian matrix) at
!-------- current estimates
          if (arith == 'f' .and. .not. eqscale) then
!
!-------------***********
!          write (*,'(a)') 'calling lshbiv'
              call lshbiv(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                          score,hess,ncov,flag,ifail,maxcol,link,risk, &
                          corr,iter,cflag,step,scnorm,scores,n1lev, &
                          family,inorm,ind,xil,sign,iquad,aquad,mus, &
                          taus,offlag,offpos)
!-------------***********
!
          else if (.not. eqscale) then
!
!-------------***************
!             write (*,'(a)') 'calling lshess_acc'
              call lshbiv_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                              nm,score,hess,ncov,flag,ifail,maxcol,link, &
                              risk,corr,iter,cflag,step,scnorm,scores, &
                              n1lev,family,inorm,ind,xil,sign,iquad, &
                              aquad,mus,taus,offlag,offpos)
!-------------***************
!
          else if (arith == 'f' .and. eqscale) then
!
!-------------***********
!          write (*,'(a)') 'calling lshbeq'
              call lshbeq(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                          score,hess,ncov,flag,ifail,maxcol,link,risk, &
                          iter,cflag,step,scnorm,scores,n1lev,family, &
                          inorm,ind,xil,sign,iquad,aquad,mus,taus, &
                          offlag,offpos)
!-------------***********
!
          else if (eqscale) then
!
!-------------***************
!             write (*,'(a)') 'calling lshbeq_acc'
              call lshbeq_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                              nm,score,hess,ncov,flag,ifail,maxcol,link, &
                              risk,iter,cflag,step,scnorm,scores,n1lev, &
                              family,inorm,ind,xil,sign,iquad,aquad,mus, &
                              taus,offlag,offpos)
!-------------***************
!
          end if
!
      else if (order .and. ilfit == 2) then
!
!-------- get log-likelihood (, score vector and Hessian matrix) at
!-------- current estimates
          if (arith == 'f') then
!
!-------------***********
!          write (*,'(a)') 'calling lshord'
              call lshord(x,y,nmes,xll,it,nsub,beta,nest,nilev,score, &
                          hess,ncov,flag,ifail,offlag,offpos,maxcol, &
                          link,iter,cflag,step,scnorm,ncat,maxcat)
!-------------***********
!
          else
!
!-------------***************
!             write (*,'(a)') 'calling lshord_acc'
              call lshord_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev, &
                              score,hess,ncov,flag,ifail,offlag,offpos, &
                              maxcol,link,iter,cflag,step,scnorm,ncat, &
                              maxcat)
!-------------***************
!
          end if
!
      else if (order .and. depend) then
!
          if (arith == 'f') then
!
!------------ get log-likelihood (, score vector and Hessian matrix) at
!------------ current 
              call lsdore(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                          score,hess,ncov,flag,ifail,offlag,offpos, &
                          maxcol,link,iter,cflag,step,scnorm,scores, &
                          ncat,ind,xil,sign,iquad,aquad,mus,taus,maxcat, &
                          risk)
!-------------***********
!
          else
!
!------------ get log-likelihood (, score vector and Hessian matrix) at
!------------ current estimates
!-------------***************
              call lsdore_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                              nm,score,hess,ncov,flag,ifail,offlag, &
                              offpos,maxcol,link,iter,cflag,step,scnorm, &
                              scores,ncat,ind,xil,sign,iquad,aquad,mus, &
                              taus,maxcat,risk)
!-------------***************
!
          end if
!
      else if (order) then
!
!-------- get log-likelihood (, score vector and Hessian matrix) at
!-------- current estimates
          if (arith == 'f') then
!
!-------------***********
!          write (*,'(a)') 'calling lshore'
              call lshore(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                          score,hess,ncov,flag,ifail,offlag,offpos, &
                          maxcol,link,iter,cflag,step,scnorm,scores, &
                          ncat,ind,xil,sign,iquad,aquad,mus,taus,maxcat)
!-------------***********
!
          else
!
!-------------***************
!             write (*,'(a)') 'calling lshore_acc'
              call lshore_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                              nm,score,hess,ncov,flag,ifail,offlag, &
                              offpos,maxcol,link,iter,cflag,step,scnorm, &
                              scores,ncat,ind,xil,sign,iquad,aquad,mus, &
                              taus,maxcat)
!-------------***************
!
          end if
!
!---- multilevel model
      else if (nlevel == 2) then
!
          if (arith == 'f') then
!
!------------ get log-likelihood (, score vector and Hessian matrix) at
!------------ current estimates
!-------------***********
!          write (*,'(a)') 'calling lshmlm'
              call lshmlm(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                          score,hess,ncov,flag,ifail,offlag,offpos, &
                          maxcol,link,iter,cflag,step,scnorm,scores, &
                          family,maxit,ind,xil,sign,iquad,aquad,mus, &
                          taus)
!-------------***********
!
          else
!
!-------------***************
!          write (*,'(a)') 'calling lshmlm_acc'
              call lshmlm_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                              nm,score,hess,ncov,flag,ifail,offlag, &
                              offpos,maxcol,link,iter,cflag,step,scnorm, &
                              scores,family,maxit,ind,xil,sign,iquad, &
                              aquad,mus,taus)
!-------------***************
!
          end if
!
      else if (depend) then
!
          if (arith == 'f') then
!
!------------ get log-likelihood (, score vector and Hessian matrix) at
!------------ current estimates
!-------------***********
!             write (*,'(a)') 'calling lshdep'
              call lshdep(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad,nm, &
                          score,hess,ncov,flag,ifail,offlag,offpos, &
                          maxcol,link,iter,cflag,step,scnorm,scores, &
                          family,n1lev,risk,ind,xil,sign,iquad,aquad, &
                          mus,taus)
!-------------***********
!
          else
!
!------------ get log-likelihood (, score vector and Hessian matrix) at
!------------ current estimates
!-------------***************
!             write (*,'(a)') 'calling lshdep_acc'
              call lshdep_acc(x,y,nmes,xll,it,nsub,beta,nest,nilev,quad, &
                              nm,score,hess,ncov,flag,ifail,offlag, &
                              offpos,maxcol,link,iter,cflag,step,scnorm, &
                              scores,family,n1lev,risk,ind,xil,sign, &
                              iquad,aquad,mus,taus)
!-------------***************
!
          end if
!
      else
!
          if (arith == 'f') then
!
!------------ get log-likelihood (, score vector and Hessian matrix) at
!------------ current estimates
!-------------***********
!             write (*,'(a)') 'calling lshuni'
              call lshuni(x,y,nmes,xll,produc,it,nsub,beta,nest,nilev, &
                          quad,nm,score,hess,ncov,flag,endind,ifail, &
                          offlag,offpos,maxcol,link,iter,cflag,step, &
                          scnorm,scores,family,inorm,ind,xil,sign,iquad, &
                          aquad,mus,taus)
!-------------***********
!
          else
!
!------------ get log-likelihood (, score vector and Hessian matrix) at
!------------ current estimates
!-------------***************
!             write (*,'(a)') 'calling lshuni_acc'
              call lshuni_acc(x,y,nmes,xll,produc,it,nsub,beta,nest, &
                              nilev,quad,nm,score,hess,ncov,flag,endind, &
                              ifail,offlag,offpos,maxcol,link,iter, &
                              cflag,step,scnorm,scores,family,inorm, &
                              ind,xil,sign,iquad,aquad,mus,taus)
!-------------***************
!
          end if
!
      end if
!
!---- overflow or underflow occurred
      if (arith == 'f' .and. ifail) then
          call newlin
          call wrtlin('    *** ERROR *** Overflow/Underflow. '// &
                      'Use the `ARITHMETIC ACCURATE` command.')
          call wterse('    *** ERROR *** Overflow/Underflow. '// &
                      'Use the `ARITHMETIC ACCURATE` command.')
      else if (ifail) then
          call newlin
          call wrtlin('    *** ALGORITHM FAILED')
          call wterse('    *** ALGORITHM FAILED')
      end if
!
      if (qflag .and. ind == 0 .and. ilfit == 0) then
!
!---------**********
          call clock('',clock1,2)
!---------**********
!
      end if
!
      return
!
      end subroutine lshess
!
!***********************************************************************
!
      subroutine covmat(hess,cov,alias,adjust,free,work,ncov,doadd, &
                        value,tol,ifail,arith,nest)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character arith
      integer ncov,nest,free
      double precision cov(ncov),adjust,tol,hess(ncov),value, &
             alias(nest),work(nest*(nest+2))
      logical doadd,ifail
!-----------------------------------------------------------------------
!     Function : Calculates elements of the variance-covariance matrix.
!-----------------------------------------------------------------------
!     ADJUST : scale for the diagonal elements of the Hessian matrix
!     FREE   : number of free parameters in current model
!     DOADD  : if true, then Hessian diagonals are incremented
!     VALUE  : value by which Hessian diagonals are incremented
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      if (arith == 'f') then
!
!---------***********
          call covfas(hess,cov,alias,adjust,free,work,free*(free+1)/2, &
                      doadd,value,tol,ifail,arith,nest)
!---------***********
!
      else
!
!---------***********
!          call covacc(hess,cov,alias,adjust,free,work,free*(free+1)/2, &
!                      doadd,value,tol,ifail,arith,nest)
          call covfas(hess,cov,alias,adjust,free,work,free*(free+1)/2, &
                      doadd,value,tol,ifail,arith,nest)
!---------***********
!
      end if
!
      return
!
      end subroutine covmat
!
!***********************************************************************
!
      subroutine covfas(hess,cov,alias,adjust,free,work,ncov,doadd, &
                        value,tol,ifail,arith,nest)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character arith
      integer ncov,nest,free
      double precision cov(ncov),adjust,tol,hess(ncov),value, &
             alias(nest),work(nest*(nest+2))
      logical doadd,ifail
!-----------------------------------------------------------------------
!     Function : Calculates elements of the variance-covariance matrix.
!-----------------------------------------------------------------------
!     ADJUST : scale for the diagonal elements of the Hessian matrix
!     FREE   : number of free parameters in current model
!     DOADD  : if true, then Hessian diagonals are incremented
!     VALUE  : value by which Hessian diagonals are incremented
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      integer i,j,l
!-----------------------------------------------------------------------
      l = 0
!
      do 51 i = 1,free
!
          do 50 j = 1,i
              l = l+1
!
              if (i == j) then
!
                  if (doadd) then
!-------------------- increment diagonal elements of Hessian
                      hess(l) = hess(l) + value
                  else
!-------------------- scale diagonal elements of Hessian
                      hess(l) = hess(l) * adjust
                  end if
!
              end if
!
              work(l) = hess(l)
   50     end do
!
   51 end do
!
      ifail = .true.
!
!---- invert Hessian matrix
!-----***********
      call syminv(work,free,work(ncov+1),work(2*ncov+1),alias,ncov,tol, &
                  ifail,arith)
!-----***********
!
      if (ifail) then
          return
      end if
!
!---- store elements of the covariance matrix
      do 80 i = 1,ncov
          cov(i) = work(ncov+i)
   80 end do
!
      return
!
      end subroutine covfas
!
!***********************************************************************
!
      subroutine covacc(hess,cov,alias,adjust,free,work,ncov,doadd, &
                        value,tol,ifail,arith,nest)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      character arith
      integer ncov,nest,free
      double precision cov(ncov),adjust,tol,hess(ncov),value, &
             work(nest*(nest+2)),alias(nest)
      logical doadd,ifail
!-----------------------------------------------------------------------
!     Function : Calculates elements of the variance-covariance matrix.
!-----------------------------------------------------------------------
!     ADJUST : scale for the diagonal elements of the Hessian matrix
!     FREE   : number of free parameters in current model
!     DOADD  : if true, then Hessian diagonals are incremented
!     VALUE  : value by which Hessian diagonals are incremented
!-----------------------------------------------------------------------
      include 'accmac.h'
!-----------------------------------------------------------------------
      double precision dummy
      integer i,j,l,hesexp,dumexp,adjexp
!-----------------------------------------------------------------------
      adjexp = 0
!
!-----***********
      call manexp(adjust,adjexp)
!-----***********
!
      l = 0
!
      do 51 i = 1,free
!
          do 50 j = 1,i
              l = l+1
              hesexp = 0
!
!-------------***********
              call manexp(hess(l),hesexp)
!-------------***********
!
              if (i == j) then
!
                  if (doadd) then
!-------------------- increment diagonal elements of Hessian
                      hess(l) = hess(l) + value
!
                  else
!-------------------- scale diagonal elements of Hessian
                      dummy = hess(l)
                      dumexp = hesexp
!
!---------------------***********
                      call multip(dummy,dumexp,adjust,adjexp, &
                                  hess(l),hesexp)
!---------------------***********
!
                  end if
!
              end if
!
              dummy = hess(l)
!
!-------------***********
              call normal(dummy,hesexp,hess(l))
!-------------***********
!
              work(l) = hess(l)
   50     end do
!
   51 end do
!
      ifail = .true.
!
!---- invert Hessian matrix
!-----***********
      call syminv(work,free,work(ncov+1),work(2*ncov+1),alias,ncov,tol, &
                  ifail,arith)
!-----***********
!
      if (ifail) then
          return
      end if
!
!---- store elements of the covariance matrix
      do 80 i = 1,ncov
          cov(i) = work(ncov+i)
   80 end do
!
      return
!
      end subroutine covacc






















!

