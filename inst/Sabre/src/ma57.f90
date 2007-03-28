! *******************************************************************
! COPYRIGHT (c) 1999 CCLRC Council for the Central Laboratory
!                    of the Research Councils
! All rights reserved.
!
! None of the comments in this Copyright notice between the lines
! of asterisks shall be removed or altered in any way.
!
! This Package is intended for compilation without modification,
! so most of the embedded comments have been removed.
!
! ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
! Licence, see http://hsl.rl.ac.uk/acuk/cou.html
!
! Please note that for a UK ACADEMIC Licence:
!
! 1. The Packages may only be used for academic research or teaching
!    purposes by the Licensee, and must not be copied by the Licensee for
!    use by any other persons. Use of the Packages in any commercial
!    application shall be subject to prior written agreement between
!    Hyprotech UK Limited and the Licensee on suitable terms and
!    conditions, which will include financial conditions.
! 2. All information on the Package is provided to the Licensee on the
!    understanding that the details thereof are confidential.
! 3. All publications issued by the Licensee that include results obtained
!    with the help of one or more of the Packages shall acknowledge the
!    use of the Packages. The Licensee will notify the Numerical Analysis
!    Group at Rutherford Appleton Laboratory of any such publication.
! 4. The Packages may be modified by or on behalf of the Licensee
!    for such use in research applications but at no time shall such
!    Packages or modifications thereof become the property of the
!    Licensee. The Licensee shall make available free of charge to the
!    copyright holder for any purpose all information relating to
!    any modification.
! 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
!    direct or consequential loss or damage whatsoever arising out of
!    the use of Packages by the Licensee.
! *******************************************************************
!
! Original date 13 September 1999
! 01/11/00  Entries in IW initialized to zero in MA57OD to avoid copy
!           of unassigned variables by MA57ED.
!           AINPUT and IINPUT reset in call to MA57ED.
! 06/02/01  Default values for ICNTL(12) and ICNTL(13) changed.
!           Control for direct addressing in solve changed to be
!           on number of rows and number columns in block pivot.
!           Several comments changed as consequence.
!           INFO(31) added to record number of block pivots.
!           Subtroutines MA57XD and MA57YD added for efficiency when only
!           one rhs (equivalent to MA57QD and MA57RD resp).
! 04/07/01  Use of MC41 changed to use of MC71.
! 26/10/01  Printing controls corrected to ensure ICNTL(5) is used and
!           unit number always checked for being positive before printing.
!           Text and comments changed to reflect that D inverse is held in
!           factors and text for solution changed from Right-hand side
!           to solution.
!           Option of choosing two 1 x 1 pivots when 2 x 2 fails
!           removed.
!           MC47B/BD given remaining length in KEEP to avoid compresses
! 20/12/01  INFO(1) initilaized to zero in MA57ED
! 06/12/02  The test for convergence of iterative refinement changed to
!           avoid any problem with comparisons of numbers held in
!           registers.
! 25/03/03  MC50 (AMD with dense row protection) and MA27 (minimum degree)
!           added.  Invoked by ICNTL(6) equal to 2 and 3, respectively.
!           Routines MA57H/HD, MA57V/VD, and MA57Z/ZD have been added to
!           duplicate routines MA27H/HD, MA27G/GD, and MA27U/UD from
!           MA57 and MC50B/BD is another internal routine of MA57.
!           ICNTL(14) has been added to control density of rows regarded
!           as dense by the MC50 and MA27 orderings.
! 24/05/04  Statment functions in MA57U/UD replaced by in-line code.

! 12th July 2004 Version 1.0.0. Version numbering added.

! 20/07/04  Several changes incorporated for HSL 2004 code.
!           Removed unused INT,ABS from MA57U/UD
!           INFO(32), INFO(33), and INFO(34) added
!           INFO(32): number of zeros in the triangle of the factors
!           INFO(33): number of zeros in the rectangle of the factors
!           INFO(34): number of zero columns in rectangle of the factors
!           Static pivoting available (controlled by CNTL(4), CNTL(5))
!           Scaling using symmetrized MC64 (ICNTL(15))
!           Links to METIS_NODEND ordering


! 31st July 2004 Version 2.0.0 established at HSL 2004 release.
! 1st Sept  2004 Version 2.1.0. Default changed to static pivoting off.
! 10th Sept 2004 Version 2.2.0. Defaults for ICNTL(6), ICNTL(9) and
!           CNTL(5) changed. Scaling factors (optionally) printed.
!  4th Nov  2004 Version 2.2.1. Change to assembly of reals in MA57OD
!           leading to more efficient code at suggestion of Stephane
!           Pralet.
! 13th Dec  2004 Version 2.3.0. Several minor changes after field
!           testing.
!           Scale factors (RINFO(16) and RINFO(17) set to 1
!           if scaling not used.
!           Option to handle dense columns invoked for METIS ordering.
!           Value of SCHNAB(1) set to 1. to allow Schnabel-Eskow to
!           work on matrix with a rows of zeros.
!           Some diagnostic printing and STOP statements removed from
!           MC50.
! 2nd March 2005  Version 3.0.0.  A new option has been added for
!           ordering the matrix.  If ICNTL(6) is equal to 5 then the
!           ordering chosen depends on the matrix characteristics.
!           At the moment the choices are MC50 or METIS.
!           INFO(36) is set to ordering used.
!           A minor chnage has been made to the pivot control to reduce
!           the amount of researching on failed pivots (restting of KR).
!           FD05 dependence changed to FD15.
! 15th June 2005  Version 3.0.1.  Setting of ALENB in MA57BD moved
!           before first error exit to avoid undefined variable
!           if error invoked.  INFO(1) initialized to zero in call to
!           MA57CD.

      subroutine ma57id(cntl, icntl)
!****************************************************************
      double precision    cntl(5)
      integer             icntl(20)
      integer i
      double precision zero
      parameter (zero=0.0d0)
!===============================================
!===============================================
      cntl(1)   = 0.01d0
      cntl(2)   = 1.0d-20
      cntl(3)   = 0.5d0
      cntl(4) = zero
      cntl(5) = zero
      icntl(1)  = 6
      icntl(2)  = 6
      icntl(3)  = 6
      icntl(4)  = -1
      icntl(5)  = 2
      icntl(6)  = 5
      icntl(7)  = 1
      icntl(8)  = 0
      icntl(9)  = 10
      icntl(10) = 0
      icntl(11) = 16
      icntl(12) = 16
      icntl(13) = 10
      icntl(14) = 100
      icntl(15) = 1
      do 110 i=16,20
        icntl(i) = 0
  110 continue
      return
      end
!CC
      subroutine ma57ad(n,ne,irn,jcn,lkeep,keep,iwork,icntl,info,rinfo)
      integer n,ne,irn(ne),jcn(ne),iwork(5*n),lkeep,keep(lkeep), &
              icntl(20),info(40)
      double precision rinfo(20)
!**** Still to be updated
      intrinsic min
      external ma57gd,mc47bd,mc50bd,ma57vd,ma57hd,ma57jd,ma57kd, &
               ma57ld,ma57md,ma57nd
      integer i,il,in,ipe,irnprm,count,fils,frere,hold,ifct,invp,ips, &
              iw,iwfr,k,ldiag,lp,lw,lrow,map,expne, &
              mp,ncmpa,nemin,node,nst,nsteps,nv,perm, &
              iw1,iw2,iw3,iw4,iw5,nstk,nd,nelim,nze,alenb, &
              j,jj,j1,j2,size22,oxo
      integer metopt(8),metftn,icntl6
      double precision zero,thresh,avnum,mc50fi
      parameter (zero=0.0d0)
      character*80 outbuf
      lp = icntl(1)
      mp = icntl(3)
      ldiag = icntl(5)
      do 10 i = 1,40
        info(i) = 0
   10 continue
      do 11 i = 1,20
        rinfo(i) = zero
   11 continue
      if (n.lt.1)  go to 20
      if (ne.lt.0) go to 30
      if (lkeep.lt.5*n+ne+max(n,ne)+42) go to 40
      if (icntl(6).eq.1) then
        do 12 i = 1,n
          iwork(i) = 0
   12   continue
        do 14 i=1,n
          k = keep(i)
          if (k.le.0 .or. k.gt.n) go to 80
          if (iwork(k).ne.0) go to 80
          iwork(k) = i
   14   continue
      endif
      if (ldiag.ge.3 .and. mp.ge.0) then
        write(mp,99980) n,ne,(icntl(i),i=1,7),icntl(12),icntl(15)
99980 format (//'Entering analysis phase (MA57AD) with ...'/ &
            'N         Order of matrix                     =',i12/ &
            'NE        Number of entries                   =',i12/ &
            'ICNTL(1)  Stream for errors                   =',i12/ &
            ' --- (2)  Stream for warnings                 =',i12/ &
            ' --- (3)  Stream for monitoring               =',i12/ &
            ' --- (4)  Stream for statistics               =',i12/ &
            ' --- (5)  Level of diagnostic printing        =',i12/ &
            ' --- (6)  Flag for input pivot order          =',i12/ &
            ' --- (7)  Numerical pivoting control (st est) =',i12/ &
            ' --- (12) Node amalgamation parameter         =',i12/ &
            ' --- (15) Scaling control (storage estimate)  =',i12)
        k = min(10,ne)
        if (ldiag.ge.4) k = ne
        write (mp,'(/A/(3(I6,A,2I8,A)))') ' Matrix entries:', &
              (i,': (',irn(i),jcn(i),')',i=1,k)
        if (k.lt.ne) write (mp,'(A)') '     . . .'
        if (icntl(6).eq.1) then
          k = min(10,n)
          if (ldiag.ge.4) k = n
          write (mp,'(A,10I6:/(7X,10I6))') ' KEEP =', (keep(i),i=1,k)
          if (k.lt.n) write (mp,'(7X,A)') '     . . .'
        end if
      end if
      iw1 = 1
      iw2 = iw1 + n
      iw3 = iw2 + n
      iw4 = iw3 + n
      iw5 = iw4 + n
      fils  = iw1
      frere = iw2
      nd    = iw3
      nelim = iw4
      nv    = iw5
      perm = 1
      nsteps = perm + n
      expne  = nsteps + 1
      hold   = expne + 1
      lrow = hold + 40
      node = lrow + n
      nstk = node + n
      map  = nstk + n
      irnprm = map + max(n,ne)
      invp  = node
      iw    = node
      ipe   = lrow
      ifct  = map
      ips   = map
      count = nstk
      keep(hold) = 0
      icntl6 = icntl(6)
      if (icntl(6).gt.5) icntl6 = 5
      if (icntl6.eq.4 .or. icntl6.eq.5) then
        metftn    = 1
        metopt(1) = 0
        keep(ipe)   = 1
        keep(ipe+1) = 2
        keep(ifct)  = 1
        call metis_nodend(1,keep(ipe),keep(ifct),metftn,metopt, &
                          keep(nstk),keep(perm))
        if (keep(perm).eq.-1) then
          if (icntl6 .eq. 4) go to 90
          icntl6 = 2
        endif
      endif
      if (icntl6.ne.1) then
      if (icntl6 .ne. 3) then
        call ma57gd(n,ne,irn,jcn,keep(ifct),keep(ipe),keep(count), &
                    keep(iw),iwfr,icntl,info)
        if (icntl6.eq.5) then
          if (icntl(7).eq.2) then
            avnum = float(iwfr+n-1)/float(n)
            if (n.ge.50000) then
              icntl6 = 4
              go to 97
            endif
            if (n.le.30000) then
              icntl6 = 2
              if (avnum.gt.100.0) icntl6 = 4
              go to 97
            endif
            if (n.gt.30000 .and. n.lt.50000) then
              if (avnum.gt.46.0) then
                icntl6 = 4
              else
                icntl6 = 2
              endif
              go to 97
            endif
          else
            avnum = float(iwfr+n-1)/float(n)
            oxo = 0
            j2 = iwfr - 1
            size22 = 0
            do 100 j = n,1,-1
              j1 = keep(ipe+j-1)
              do  99 jj = j1,j2
                if (keep(ifct+jj-1).gt.j) go to 101
   99         continue
              size22 = size22 + 1
              j2 = j1-1
  100       continue
  101       if (size22 .gt. 0) then
              do 98 i = 1,ne
                if (irn(i) .le. n-size22 &
              .and. jcn(i) .le. n-size22) then
                    avnum = float(iwfr+n-size22-1)/float(n)
                    go to 96
                endif
   98         continue
              oxo = 1
              avnum = float(iwfr-1)/float(n)
            endif
   96       if (n .ge. 100000) then
              if (avnum.gt.5.42) then
                icntl6 = 4
              else
                icntl6 = 2
              endif
              go to 97
            endif
            if (oxo.eq.1) then
              if (float(n-size22)/float(size22) .gt. 1.8d0) then
                icntl6 = 2
              else
                icntl6 = 4
              endif
              go to 97
            endif
            lw = lkeep-ifct+1
            call mc50bd(icntl(14),n,lw,keep(ipe),iwfr,keep(count), &
                        keep(ifct),iwork(nv), &
                        keep(invp),keep(perm),ncmpa,iwork(iw1), &
                        iwork(iw2),iwork(iw3),iwork(iw4))
            info(13) = ncmpa
            icntl6 = 2
            nemin    = icntl(12)
      call ma57ld(n,keep(ipe),iwork(nv),keep(ips),iwork(nelim), &
                  keep(nstk),keep(node),keep(perm), &
                  keep(nsteps),iwork(fils),iwork(frere),iwork(nd), &
                  nemin,keep(irnprm))
            nst = keep(nsteps)
      call ma57md(n,ne,irn,jcn,keep(map),keep(irnprm), &
                  keep(lrow),keep(perm), &
                  iwork(iw2),iwork(iw5))
            keep(expne) = iwork(iw5)
      call ma57nd(n,keep(lrow),keep(nstk),iwork(nelim), &
                  iwork(nd),nst,iwork(iw1),iwork(iw2), &
                  info,rinfo)
            if (float(info(5))/float(ne) .lt. 10.0) then
              go to 93
            else
              mc50fi = float(info(5))/float(ne)
        call ma57gd(n,ne,irn,jcn,keep(ifct),keep(ipe),keep(count), &
                    keep(iw),iwfr,icntl,info)
              keep(ipe+n) = iwfr
              metftn    = 1
              metopt(1) = 0
              if (n.lt.50) go to 92
              do 91 i = 1,n
                if ((keep(ipe+i)-keep(ipe+i-1)) .gt. n/10) then
                  metopt(1) = 1
                  metopt(2) = 3
                  metopt(3) = 1
                  metopt(4) = 2
                  metopt(5) = 0
                  metopt(6) = 1
                  metopt(7) = 200
                  metopt(8) = 1
                  go to 92
                endif
   91         continue
   92     call metis_nodend(n,keep(ipe),keep(ifct),metftn,metopt, &
                            keep(nstk),keep(perm))
        call ma57jd(n,ne,irn,jcn,keep(perm),keep(ifct),keep(ipe), &
                    keep(count),iwork(iw1),iwfr,icntl,info)
              lw = 2*ne
        call ma57kd(n,keep(ipe),keep(ifct),lw,iwfr,keep(perm), &
                    keep(invp),iwork(nv),iwork(iw1),ncmpa)
              info(13) = ncmpa
              nemin = icntl(12)
      call ma57ld(n,keep(ipe),iwork(nv),keep(ips),iwork(nelim), &
                  keep(nstk),keep(node),keep(perm), &
                  keep(nsteps),iwork(fils),iwork(frere),iwork(nd), &
                  nemin,keep(irnprm))
              nst = keep(nsteps)
      call ma57md(n,ne,irn,jcn,keep(map),keep(irnprm), &
                  keep(lrow),keep(perm), &
                  iwork(iw2),iwork(iw5))
              keep(expne) = iwork(iw5)
      call ma57nd(n,keep(lrow),keep(nstk),iwork(nelim), &
                  iwork(nd),nst,iwork(iw1),iwork(iw2), &
                  info,rinfo)
              if (float(info(5))/float(ne).lt.mc50fi) then
                icntl6 = 4
                go to 93
              else
                icntl6=2
        call ma57gd(n,ne,irn,jcn,keep(ifct),keep(ipe),keep(count), &
                    keep(iw),iwfr,icntl,info)
                go to 97
              endif
            endif
          endif
        endif
   97   if (icntl6.eq.4) then
          keep(ipe+n) = iwfr
          metftn    = 1
          metopt(1) = 0
          if (n.lt.50) go to 103
          do 102 i = 1,n
            if ((keep(ipe+i)-keep(ipe+i-1)) .gt. n/10) then
              metopt(1) = 1
              metopt(2) = 3
              metopt(3) = 1
              metopt(4) = 2
              metopt(5) = 0
              metopt(6) = 1
              metopt(7) = 200
              metopt(8) = 1
              go to 103
            endif
  102     continue
  103     call metis_nodend(n,keep(ipe),keep(ifct),metftn,metopt, &
                            keep(nstk),keep(perm))
          go to 111
        endif
      lw = lkeep-ifct+1
        if (icntl6.eq.2) then
          call mc50bd(icntl(14),n,lw,keep(ipe),iwfr,keep(count), &
                      keep(ifct),iwork(nv), &
                      keep(invp),keep(perm),ncmpa,iwork(iw1), &
                      iwork(iw2),iwork(iw3),iwork(iw4))
        else
          call mc47bd(n,lw,keep(ipe),iwfr,keep(count), &
                      keep(ifct),iwork(nv), &
                      keep(invp),keep(perm),ncmpa,iwork(iw1), &
                      iwork(iw2),iwork(iw3),iwork(iw4))
        endif
        info(13) = ncmpa
        else
      lw = lkeep-ifct+1
         call ma57vd(n,ne,irn,jcn,keep(ifct),lw,keep(ipe),iwork(iw1), &
                     iwork(iw2),iwfr,icntl,info)
         thresh = float(icntl(14))/100.0
         call ma57hd(n,keep(ipe),keep(ifct),lw,iwfr,iwork(nv), &
                     iwork(iw1),iwork(iw2),iwork(iw3),iwork(iw4), &
                     2139062143,info(13),thresh)
      do 110 i = 1,n
        if (iwork(nv+i-1).ne.0) go to 110
        in = i
  105   il = in
        in = - keep(ipe+il-1)
        if (iwork(nv+in-1).eq.0) go to 105
        keep(ipe+i-1) = -in
  110 continue
        endif
      endif
  111 if (icntl6.eq.1 .or. icntl6.eq.4) then
        call ma57jd(n,ne,irn,jcn,keep(perm),keep(ifct),keep(ipe), &
                    keep(count),iwork(iw1),iwfr,icntl,info)
        lw = 2*ne
        call ma57kd(n,keep(ipe),keep(ifct),lw,iwfr,keep(perm), &
                    keep(invp),iwork(nv),iwork(iw1),ncmpa)
        info(13) = ncmpa
      end if
      nemin = icntl(12)
      call ma57ld(n,keep(ipe),iwork(nv),keep(ips),iwork(nelim), &
                  keep(nstk),keep(node),keep(perm), &
                  keep(nsteps),iwork(fils),iwork(frere),iwork(nd), &
                  nemin,keep(irnprm))
      nst = keep(nsteps)
      call ma57md(n,ne,irn,jcn,keep(map),keep(irnprm), &
                  keep(lrow),keep(perm), &
                  iwork(iw2),iwork(iw5))
      keep(expne) = iwork(iw5)
      call ma57nd(n,keep(lrow),keep(nstk),iwork(nelim), &
                  iwork(nd),nst,iwork(iw1),iwork(iw2), &
                  info,rinfo)
   93 info(36) = icntl6
      alenb    = 1
      if (icntl(7).eq.4) alenb = alenb + n + 5
      if (icntl(15).eq.1) alenb = alenb + n
      info(9)  = max(info(9)+alenb,alenb+keep(expne)+1)
      info(11) = max(info(11)+alenb,alenb+keep(expne)+1)
      info(10) = max(info(10),keep(expne)+n+5)
      info(12) = max(info(12),keep(expne)+n+5)
      if (icntl(15).eq.1) then
        info(9) = max(info(9),alenb+3*keep(expne)+3*n)
        info(11) = max(info(11),alenb+3*keep(expne)+3*n)
        info(10) = max(info(10),3*keep(expne)+5*n+1)
        info(12) = max(info(12),3*keep(expne)+5*n+1)
      endif
      if (ldiag.ge.3 .and. mp.ge.0) then
        nze = keep(expne)
        write (mp,99999) info(1),nze, &
                        (info(i),i=3,13),info(36),(rinfo(i),i=1,2)
99999 format (/'Leaving analysis phase (MA57AD) with ...'/ &
          'INFO(1)  Error indicator                      =',i12/ &
          'Number of entries in matrix with diagonal     =',i12/ &
          'INFO(3)  Number of out-of-range indices       =',i12/ &
          'INFO(4)  Number of off-diagonal duplicates    =',i12/ &
          'INFO(5)  Forecast real storage for factors    =',i12/ &
          '----(6)  Forecast integer storage for factors =',i12/ &
          '----(7)  Forecast maximum front size          =',i12/ &
          '----(8)  Number of nodes in assembly tree     =',i12/ &
          '----(9)  Size of FACT without compress        =',i12/ &
          '----(10) Size of IFACT without compress       =',i12/ &
          '----(11) Size of FACT with compress           =',i12/ &
          '----(12) Size of IFACT with compress          =',i12/ &
          '----(13) Number of compresses                 =',i12/ &
          '----(36) Ordering strategy used by code       =',i12/ &
          'RINFO(1) Forecast additions for assembly      =',1p,d12.5/ &
          'RINFO(2) Forecast ops for elimination         =',1p,d12.5)
        k = min(10,n)
        if (ldiag.ge.4) k = n
        write (mp,'(/A/(5I12))')  'Permutation array:', &
                        (keep(i),i=1,k)
        if (k.lt.n) write (mp,'(16X,A)') ' . . .'
        write (mp,'(/A/(5I12))') &
              'Number of entries in rows of permuted matrix:', &
              (keep(lrow+i-1),i=1,k)
        if (k.lt.n) write (mp,'(16X,A)') ' . . .'
        k = min(10,nze)
        if (ldiag.ge.4) k = nze
        write (mp,'(/A/(5I12))') &
              'Column indices of permuted matrix:', &
                                 (keep(irnprm+i-1),i=1,k)
        if (k.lt.nze) write (mp,'(16X,A)') '     . . .'
        k = min(10,n)
        if (ldiag.ge.4) k = n
        write (mp,'(/A/(5I12))') &
          'Tree nodes at which variables eliminated:', &
          (keep(node+i-1),i=1,k)
        if (k.lt.n) write (mp,'(16X,A)') ' . . .'
        k = min(10,ne)
        if (ldiag.ge.4) k = ne
        write (mp,'(/A/(5I12))') 'Map array:', &
                                     (keep(i),i=map,map+k-1)
        if (k.lt.ne) write (mp,'(16X,A)') ' . . .'
      end if
      return
   20 info(1) = -1
      info(2) = n
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57AD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') 'N has value ',info(2)
         call wrtlin( outbuf )
      end if
      return
   30 info(1) = -2
      info(2) = ne
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57AD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') 'NE has value',info(2)
         call wrtlin( outbuf )
      end if
      return
   40 info(1) = -15
      info(2) = lkeep
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57AD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') 'LKEEP has value    ',info(2)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') 'Should be at least ',5*n+ne+max(n,ne)+42
         call wrtlin( outbuf )
      end if
      return
   80 info(1) = -9
      info(2) = i
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57AD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         call wrtlin( 'Invalid permutation supplied in KEEP' )
         write (outbuf,'(A,I10,A)') 'Component',info(2),' is faulty'
         call wrtlin( outbuf )
      end if
      return
   90 info(1) = -18
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57AD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         call wrtlin( 'MeTiS ordering requested but MeTiS not linked' )
      end if
      end
!CC
!--------------------------------------------------------------------
!-             Copyright Rutherford Appleton Laboratory
!--------------------------------------------------------------------
      subroutine ma57bd(n, ne, a, fact, lfact, ifact, lifact, &
       lkeep, keep, ppos, icntl, cntl, info, rinfo)
      integer n,ne,lfact,lifact,lkeep
      double precision a(ne),fact(lfact)
      double precision rinfo(20)
      double precision cntl(5)
      integer icntl(20), ifact(lifact)
      integer   info(40), keep(lkeep), ppos(n)
      integer expne,hold,i,irnprm,k,ldiag,llfact,lp,lrow,map,mm1,mm2,mp
      integer j,jj,kk,iscale,num,ne64,idup,imat,ipt,jloop,jnew,nn,ising
      integer nsteps,node,nstk,perm,inew,alenb,biga
      double precision one,zero,rinf,fd15ad,fct,smax,smin,reps
      parameter (one = 1.0d0, zero=0.0d0)
      character*80 outbuf
!?? To identify bug
      intrinsic min
      external ma57od,ma57ud,fd15ad,mc34ad,mc64wd
      rinf = fd15ad('H')
      reps = fd15ad('E')
      lp     = icntl(1)
      mp     = icntl(3)
      ldiag  = icntl(5)
!??
      if (n.le.0)  go to 25
      if (ne.lt.0) go to 30
      if (lkeep.lt.5*n+ne+max(n,ne)+42) go to 40
      if (icntl(7).lt.1 .or. icntl(7).gt.4) go to 35
      nsteps = keep(n+1)
      expne  = keep(n+2)
      perm = 1
      hold = perm + n + 2
      lrow = hold + 40
      node = lrow + n
      nstk = node + n
      map  = nstk + n
      irnprm = map + max(ne,n)
      biga = lfact
      llfact = lfact - 1
      if (icntl(15).eq.1) then
        iscale = llfact - n + 1
        llfact = iscale - 1
      endif
      if (icntl(7).eq.4) then
        llfact = llfact - n - 5
        mm1 = llfact+6
        mm2 = llfact+1
      else
        mm1 = 1
        mm2 = 1
      endif
      alenb = 1
      if (icntl(7).eq.4)  alenb = alenb + n + 5
      if (icntl(15).eq.1) alenb = alenb + n
      if (llfact.lt.expne+1)   go to 85
      if (lifact.lt.expne+n+5)  go to 95
      if (icntl(15).eq.1)  then
        if (lfact .lt. alenb + 3*expne  + 3*n) go to 85
        if (lifact .lt. 3*expne + 5*n + 1) go to 95
      endif
!*****************************
      if (ldiag.ge.3 .and. mp.ge.0) then
        write (mp,99999)
99999 format (//'Entering factorization phase (MA57BD) with ...')
        if (keep(hold).gt.0) write (mp,99998)
99998 format ('Re-entry call after call to MA57ED')
        write (mp,99997) n,ne,expne,(icntl(i),i=1,5),icntl(7),icntl(8), &
               icntl(11),icntl(15),lfact, lifact, nsteps, &
               cntl(1), cntl(2), cntl(4), cntl(5)
99997 format ('N       Order of input matrix               =',i12/ &
              'NE      Entries in input matrix             =',i12/ &
              '        Entries in input matrix (inc diags) =',i12/ &
              'ICNTL(1)  Stream for errors                 =',i12/ &
              ' --- (2)  Stream for warnings               =',i12/ &
              ' --- (3)  Stream for monitoring             =',i12/ &
              ' --- (4)  Stream for statistics             =',i12/ &
              ' --- (5)  Level of diagnostic printing      =',i12/ &
              ' --- (7)  Numerical pivoting control        =',i12/ &
              ' --- (8)  Restart or discard factors        =',i12/ &
              ' --- (11) Block size for Level 3 BLAS       =',i12/ &
              ' --- (15) Scaling control (1 on)            =',i12/ &
              'LFACT   Size of real working space          =',i12/ &
              'LIFACT  Size of integer working space       =',i12/ &
              '        Number nodes in assembly tree       =',i12/ &
              'CNTL(1) Value of threshold parameter        =',d12.5/ &
              'CNTL(2) Threshold for zero pivot            =',d12.5/ &
              'CNTL(4) Control for value of static pivots  =',d12.5/ &
              'CNTL(5) Control for number delayed pivots   =',d12.5)
        k = min(10,ne)
        if (ldiag.ge.4) k = ne
        if (ne.gt.0) then
          write (mp,'(/A/(3(I6,A,1P,D16.8,A)))') 'Matrix entries:', &
           (i,': (',a(i),')',i=1,k)
          if (k.lt.ne) write (mp,'(A)') '     . . .'
        end if
        k = min(10,n)
        if (ldiag.ge.4) k = n
        write (mp,'(/A/(5I12))')  'Permutation array:', &
                          (keep(i),i=1,k)
        if (k.lt.n) write (mp,'(16X,A)') ' . . .'
        write (mp,'(/A/(5I12))') &
                'Number of entries in rows of permuted matrix:', &
                (keep(lrow+i-1),i=1,k)
        if (k.lt.n) write (mp,'(16X,A)') ' . . .'
        write (mp,'(/A/(5I12))') &
          'Tree nodes at which variables eliminated:', &
          (keep(node+i-1),i=1,k)
        if (k.lt.n) write (mp,'(16X,A)') ' . . .'
        k = min(10,nsteps)
        if (ldiag.ge.4) k = nsteps
        if (k.gt.0) write (mp,'(/A/(5I12))') &
           'Number of assemblies at each tree node:', &
           (keep(nstk+i-1),i=1,k)
        if (k.lt.nsteps) write (mp,'(16X,A)') ' . . .'
        k = min(10,ne)
        if (ldiag.ge.4) k = ne
        write (mp,'(/A/(5I12))') 'Map array:', &
                                     (keep(i),i=map,map+k-1)
        if (k.lt.ne) write (mp,'(16X,A)') ' . . .'
        k = min(10,expne)
        if (ldiag.ge.4) k = expne
        write (mp,'(/A/(5I12))') &
                'Column indices of permuted matrix:', &
                                   (keep(irnprm+i-1),i=1,k)
        if (k.lt.expne) write (mp,'(16X,A)') '     . . .'
      endif
      if (keep(hold) .gt. 0) go to 22
!***************************************************
!***************************************************
!?? For the moment to handle missing diagonals
      do 19 k = 1,expne
        fact(llfact-expne+k) = zero
   19 continue
      fact(biga) = zero
      do 20 k = 1,ne
        fact(biga) = max(fact(biga),abs(a(k)))
        fact(keep(map+k-1)+llfact-expne) = a(k)
   20 continue
      rinfo(18) = fact(biga)
      do 21 k = 1,expne
        ifact(lifact-expne+k) = keep(irnprm+k-1)
   21 continue
      do 23 i = 1,n
        ppos(keep(perm+i-1)) = i
   23 continue
      if (icntl(15).eq.1) then
        ipt = 1
        idup = ipt+n+1
        imat = idup+n
        ising = imat + max(ne,expne)
        do 4444 i = 1,n
          ifact(idup+i-1) = 0
 4444   continue
!9999   CONTINUE
        ifact(ipt) = 1
        kk = 1
        k = 1
        do 3333 j = 1,n
          do 2222 jj = 1,keep(lrow+j-1)
            i = keep(perm+ifact(lifact-expne+k)-1)
            if (ifact(idup+i-1).ge.ifact(ipt+j-1)) then
              fact(ifact(idup+i-1)) = &
                fact(ifact(idup+i-1)) + fact(llfact-expne+k)
            else
              if (fact(llfact-expne+k).ne.zero) then
                ifact(idup+i-1) = kk
                fact(kk) = fact(llfact-expne+k)
                ifact(imat-1+kk) = i
                kk = kk+1
              endif
            endif
            k = k + 1
 2222     continue
          ifact(ipt+j) = kk
 3333   continue
        call mc34ad(n,ifact(imat),ifact(ipt),.true.,fact,keep(perm))
        ne64 = ifact(ipt+n)-1
        do 75 j = 1,n
          fct = zero
          do 60 k = ifact(ipt+j-1),ifact(ipt+j)-1
            fact(k) = abs(fact(k))
            if (fact(k).gt.fct) fct = fact(k)
   60     continue
          fact(ne64+2*n+j) = fct
          if (fct.ne.zero) then
            fct = log(fct)
          else
            fct = rinf/n
          endif
          do 70 k = ifact(ipt+j-1),ifact(ipt+j)-1
            if (fact(k).ne.zero) then
              fact(k) = fct - log(fact(k))
            else
              fact(k) = rinf/n
            endif
   70     continue
   75   continue
        call mc64wd(n,ne64,ifact(ipt),ifact(imat),fact,keep(perm),num, &
            ifact(idup),ifact(imat+ne64),ifact(imat+ne64+n), &
            ifact(imat+ne64+2*n),ifact(imat+ne64+3*n), &
            fact(ne64+1),fact(ne64+n+1))
        if (num.eq.n) then
          do 80 j = 1,n
            if (fact(ne64+2*n+j).ne.zero) then
              fact(ne64+n+j) = fact(ne64+n+j) - log(fact(ne64+2*n+j))
            else
              fact(ne64+n+j) = zero
            endif
   80     continue
          do 5555 i=1,n
            fact(iscale+ppos(i)-1) = &
              sqrt(exp(fact(ne64+i)+fact(ne64+n+i)))
 5555     continue
        else
        k = 0
        do 3501 i = 1,n
          if (keep(perm+i-1).lt.0) then
            ppos(i) = -ppos(i)
            ifact(ising+i-1) = 0
          else
            k = k + 1
            ifact(ising+i-1) = k
          endif
 3501   continue
        do 3502 i = 1,n
          keep(perm+abs(ppos(i))-1) = i
 3502   continue
        do 3503 i = 1,n
          ifact(idup+i-1) = 0
 3503   continue
        ifact(ipt) = 1
        kk = 1
        k = 1
        jnew = 0
        nn = n
        do 3505 j = 1,n
          if (ppos(j).lt.0) then
            nn = nn - 1
            k = k + keep(lrow+j-1)
            go to 3505
          endif
          jnew = jnew + 1
          do 3504 jj = 1,keep(lrow+j-1)
            i = keep(perm+ifact(lifact-expne+k)-1)
            if (ppos(i).gt.0) then
              if (ifact(idup+i-1).ge.ifact(ipt+j-1)) then
                fact(ifact(idup+i-1)) = &
                  fact(ifact(idup+i-1)) + fact(llfact-expne+k)
              else
                if (fact(llfact-expne+k).ne.zero) then
                  ifact(idup+i-1) = kk
                  fact(kk) = fact(llfact-expne+k)
                  ifact(imat-1+kk) = ifact(ising+i-1)
                  kk = kk+1
                endif
              endif
            endif
            k = k + 1
 3504     continue
          ifact(ipt+jnew) = kk
 3505   continue
      ne64 = ifact(ipt+nn)-1
        call mc34ad(nn,ifact(imat),ifact(ipt),.true.,fact,keep(perm))
        ne64 = ifact(ipt+nn)-1
        do 3508 j = 1,nn
          fct = zero
          do 3506 k = ifact(ipt+j-1),ifact(ipt+j)-1
            fact(k) = abs(fact(k))
            if (fact(k).gt.fct) fct = fact(k)
 3506     continue
          fact(ne64+2*n+j) = fct
          if (fct.ne.zero) then
            fct = log(fct)
          else
            fct = rinf/nn
          endif
          do 3507 k = ifact(ipt+j-1),ifact(ipt+j)-1
            if (fact(k).ne.zero) then
              fact(k) = fct - log(fact(k))
            else
              fact(k) = rinf/nn
            endif
 3507     continue
 3508   continue
        call mc64wd(nn,ne64,ifact(ipt),ifact(imat),fact,keep(perm),num, &
            ifact(idup),ifact(imat+ne64),ifact(imat+ne64+n), &
            ifact(imat+ne64+2*n),ifact(imat+ne64+3*n), &
            fact(ne64+1),fact(ne64+n+1))
        do 3509 j = 1,nn
            if (fact(ne64+2*n+j).ne.zero) then
              fact(ne64+n+j) = fact(ne64+n+j) - log(fact(ne64+2*n+j))
            else
              fact(ne64+n+j) = zero
            endif
 3509     continue
          k=0
          do 3510 i=1,n
            if (ppos(i).lt.0) then
              k = k + 1
              fact(iscale-ppos(i)-1) = zero
            else
              fact(iscale+ppos(i)-1) = &
                sqrt(exp(fact(ne64+i-k)+fact(ne64+n+i-k)))
            endif
 3510     continue
          do 3516 i = 1,n
            keep(perm+abs(ppos(i))-1) = i
 3516     continue
          k = 1
          do 3514 jj = 1,n
            j = ppos(jj)
            if (j.gt.0) then
              do 3511 jloop = 1,keep(lrow+jj-1)
                i = ifact(lifact-expne+k)
                inew = keep(perm+i-1)
                if (ppos(inew).lt.0) &
                  fact(iscale+i-1) = max(fact(iscale+i-1), &
                       abs(fact(llfact-expne+k))*fact(iscale+j-1))
                k = k + 1
 3511         continue
            else
              do 3512 jloop = 1,keep(lrow+jj-1)
                i = ifact(lifact-expne+k)
                inew = keep(perm+i-1)
                if (i .ne. -j)  then
                fact(iscale-j-1) = &
                    max(fact(iscale-j-1), &
                    abs(fact(llfact-expne+k))*fact(iscale+i-1))
                endif
                k = k + 1
 3512         continue
            endif
 3514     continue
          do 3513 i = 1,n
            inew = keep(perm+i-1)
            if (ppos(inew) .lt. 0) then
              ppos(inew) = - ppos(inew)
              if (fact(iscale+i-1) .eq. zero) then
                fact(iscale+i-1) = one
              else
                fact(iscale+i-1) = one/fact(iscale+i-1)
              endif
            endif
 3513     continue
        endif
!8888     CONTINUE
          smax = fact(iscale)
          smin = fact(iscale)
          do 5566 i = 1,n
            smax = max(smax,fact(iscale+i-1))
            smin = min(smin,fact(iscale+i-1))
 5566     continue
          rinfo(16) = smin
          rinfo(17) = smax
          k = 1
          fact(biga) = zero
          do 6666 jj = 1,n
            j = ppos(jj)
            do 7777 jloop = 1,keep(lrow+jj-1)
              i = ifact(lifact-expne+k)
              fact(llfact-expne+k) = &
                fact(iscale+i-1)*fact(llfact-expne+k)*fact(iscale+j-1)
              fact(biga) = max(fact(biga), abs(fact(llfact-expne+k)))
              k = k + 1
 7777       continue
 6666     continue
!PRINT
!6661     CONTINUE
!6663     CONTINUE
      else
        rinfo(16) = one
        rinfo(17) = one
      endif
!**********************************
!**********************************
   22 call ma57od(n, expne, fact, llfact, ifact, lifact, keep(lrow), &
                  ppos, &
                  nsteps, keep(nstk), keep(node), fact(mm1), &
                  fact(mm2), &
                  keep(perm), &
                  cntl, icntl, &
                  info, rinfo, keep(hold), fact(biga))
      if (info(1).eq.10 .or. info(1).eq.11) then
        if (ldiag.gt.2 .and. mp.ge.0)  then
          if (info(1).eq.10) write (mp,99982) info(1)
99982 format (/'Leaving factorization phase (MA57BD) with ...'/ &
        'Factorization suspended because of lack of real space'/ &
        'INFO (1) = ',i3)
          if (info(1).eq.11) write (mp,99983) info(1)
99983 format (/'Leaving factorization phase (MA57BD) with ...'/ &
        'Factorization suspended because of lack of integer space'/ &
        'INFO (1) = ',i3)
        endif
        return
      endif
      do 24 i = 1,n
        keep(perm+ppos(i)-1) = i
   24 continue
        info(17) = alenb + info(17)
        info(19) = alenb + info(19)
      if (icntl(15).eq.1) then
        info(17) = max(info(17),alenb + 3*expne+3*n)
        info(19) = max(info(19),alenb + 3*expne+3*n)
        info(18) = max(info(18),3*expne+5*n+1)
        info(20) = max(info(18),3*expne+5*n+1)
      endif
      if (info(1).lt.0) return
      go to 100
!************************
!************************
   25 info(1) = -1
      info(2) =  n
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57BD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') 'N has value ',info(2)
         call wrtlin( outbuf )
      end if
      return
   30 info(1) = -2
      info(2) = ne
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57BD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') 'NE has value',info(2)
         call wrtlin( outbuf )
      end if
      return
   40 info(1) = -15
      info(2) = lkeep
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57BD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') 'LKEEP has value    ',info(2)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') 'Should be at least ',5*n+ne+max(n,ne)+42
         call wrtlin( outbuf )
      end if
      return
   35 info(1) = -10
      info(2) = icntl(7)
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57BD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') 'ICNTL(7) has value',icntl(7)
         call wrtlin( outbuf )
      end if
      return
   85 info(1) = -3
      info(2) = lfact
      info(17) = alenb + expne + 1
      if (icntl(15).eq.1) info(17) = alenb + 3*expne + 3*n
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write (outbuf,'(A,I3)') &
          '**** Error return from MA57BD ****  INFO(1) =',info(1)
         call wrtlin( outbuf )
         write (outbuf,'(A,I10)') &
          'Insufficient real space in FACT, LFACT = ',info(2)
         call wrtlin( outbuf )
      end if
      return
   95 info(1) = -4
      info(2) = lifact
      info(18) = expne+n+5
      if (icntl(15).eq.1) info(18) = 3*expne + 5*n + 1
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
	 write (outbuf,'(A,I3)') &
	  '**** Error return from MA57BD ****  INFO(1) =',info(1)
	 call wrtlin( outbuf )
	 write (outbuf,'(A,I10)') &
	  'Insufficient integer space in IFACT, LIFACT = ',info(2)
	 call wrtlin( outbuf )
      end if
      return
!****************
!****************
 100  if (ldiag.le.2 .or. mp.lt.0) return
      write (mp,99980) info(1), info(2), &
	  (info(i),i=14,25),info(28),info(29)
      write (mp,99984) (info(i),i=31,35),rinfo(3), rinfo(4), &
		       rinfo(5), rinfo(18)
99980 format (/'Leaving factorization phase (MA57BD) with ...'/ &
	'INFO (1)				       =',i12/ &
	' --- (2)				       =',i12/ &
	' --- (14) Number of entries in factors        =',i12/ &
	' --- (15) Real storage for factors	       =',i12/ &
	' --- (16) Integer storage for factors	       =',i12/ &
	' --- (17) Min LFACT with compresses	       =',i12/ &
	' --- (18) Min LIFACT with compresses	       =',i12/ &
	' --- (19) Min LFACT without compresses        =',i12/ &
	' --- (20) Min LIFACT without compresses       =',i12/ &
	' --- (21) Order of largest frontal matrix     =',i12/ &
	' --- (22) Number of 2x2 pivots 	       =',i12/ &
	' --- (23) Number of delayed pivots	       =',i12/ &
	' --- (24) Number of negative eigenvalues      =',i12/ &
	' --- (25) Rank of factorization	       =',i12/ &
	' --- (28) Number compresses on real data      =',i12/ &
	' --- (29) Number compresses on integer data   =',i12)
      if (icntl(15).eq.1) write (mp,99985) rinfo(16),rinfo(17)
99985 format ( &
	'RINFO(16) Minimum value of scaling factor     =  ',1pd10.3/ &
	'-----(17) Maximum value of scaling factor     =  ',1pd10.3)
99984 format ( &
	' --- (31) Number of block pivots in factors   =',i12/ &
	' --- (32) Number of zeros factors triangle    =',i12/ &
	' --- (33) Number of zeros factors rectangle   =',i12/ &
	' --- (34) Number of zero cols factors rect    =',i12/ &
	' --- (35) Number of static pivots	       =',i12/ &
	'RINFO(3)  Operations during node assembly     =  ',1pd10.3/ &
	'-----(4)  Operations during node elimination  =  ',1pd10.3/ &
	'-----(5)  Extra operations because of BLAS    =  ',1pd10.3/ &
	'-----(18) Largest modulus of entry in matrix  =  ',1pd10.3)
      if (info(27).gt.0) write (mp,99981) info(27),rinfo(14),rinfo(15)
99981 format (/'Matrix modification performed'/ &
	'INFO (27) Step at which matrix first modified =',i12/ &
	'RINFO(14) Maximum value added to diagonal     =  ',1pd10.3/ &
	'RINFO(15) Smallest pivot in modified matrix   =  ',1pd10.3)
      call ma57ud(fact,lfact,ifact,lifact,icntl)
      if (icntl(15).ne.1) return
      k = min(10,n)
      if (ldiag.ge.4) k = n
      write (mp,'(/A/(5D12.5))')  'Scaling factors:', &
			  (fact(iscale+i-1),i=1,k)
      if (k.lt.n) write (mp,'(16X,A)') ' . . .'
      end
      subroutine ma57cd(job,n,fact,lfact,ifact,lifact,nrhs,rhs,lrhs,w, &
			lw,iw1,icntl,info)
      integer job,n,lfact
      double precision fact(lfact)
      integer lifact,ifact(lifact),nrhs,lrhs,lw
      double precision w(lw),rhs(lrhs,nrhs)
      integer iw1(n),icntl(20),info(40)
      intrinsic min
      external ma57qd,ma57rd,ma57sd,ma57td,ma57ud,ma57xd,ma57yd
      double precision scale,one
      parameter (one = 1.0d0)
      integer i,j,k,ldiag,llw,lp,mp,iscale
      character*80 outbuf
      lp = icntl(1)
      mp = icntl(3)
      ldiag = icntl(5)
      info(1) = 0
      if (n.le.0) then
	info(1) = -1
	info(2) = n
	if (ldiag.gt.0 .and. lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           write (outbuf,'(A,I3)') &
            '**** Error return from MA57CD ****  INFO(1) =',info(1)
           call wrtlin( outbuf )
           write (outbuf,'(A,I10)') 'N has value ',n
           call wrtlin( outbuf )
        end if
        goto 500
      endif
      if (nrhs.lt.1) then
        info(1) = -16
        info(2) = nrhs
        if (ldiag.gt.0 .and. lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           write (outbuf,'(A,I4)') &
            '**** Error return from MA57CD ****  INFO(1) =',info(1)
           call wrtlin( outbuf )
           write (outbuf,'(A,I10,A)') 'value of NRHS =',nrhs,' is less than 1'
           call wrtlin( outbuf )
         end if
        goto 500
      endif
      if (lrhs.lt.n) then
        info(1) = -11
        info(2) = lrhs
        if (ldiag.gt.0 .and. lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           write (outbuf,'(A,I4)') &
            '**** Error return from MA57CD ****  INFO(1) =',info(1)
           call wrtlin( outbuf )
           write (outbuf,'(A,I10,A,I10)') &
          'value of LRHS =',lrhs,' is less than N=',n
           call wrtlin( outbuf )
        end if
        goto 500
      endif
      if (lw.lt.n*nrhs) then
        info(1) = -17
        info(2) = n*nrhs
        if (ldiag.gt.0 .and. lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           write (outbuf,'(A,I4)') &
            '**** Error return from MA57CD ****  INFO(1) =',info(1)
           call wrtlin( outbuf )
           write (outbuf,'(A,I10,A,I10)') &
          'value of LW =',lw,' is less than', n*nrhs
           call wrtlin( outbuf )
        end if
        goto 500
      endif
      if (ldiag.ge.3 .and. mp.ge.0) then
        write (mp,99999) job,n,(icntl(i),i=1,5),lfact,lifact,nrhs, &
               lrhs,lw,icntl(13)
99999 format(/'Entering solution phase (MA57CD) with ...'/ &
          'JOB       Control on coefficient matrix       =',i12/ &
          'N         Order of matrix                     =',i12/ &
          'ICNTL(1)  Stream for errors                   =',i12/ &
          ' --- (2)  Stream for warnings                 =',i12/ &
          ' --- (3)  Stream for monitoring               =',i12/ &
          ' --- (4)  Stream for statistics               =',i12/ &
          ' --- (5)  Level of diagnostic printing        =',i12/ &
          'LFACT     Length of array FACT                =',i12/ &
          'LIFACT    Length of array IFACT               =',i12/ &
          'NRHS      Number of right-hand sides          =',i12/ &
          'LRHS      Leading dimension of RHS array      =',i12/ &
          'LW        Leading dimension of work array     =',i12/ &
          'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',i12)
        call ma57ud(fact,lfact,ifact,lifact,icntl)
        if (icntl(15).eq.1) then
          iscale = lfact-n
          k = min(10,n)
          if (ldiag.ge.4) k = n
          write (mp,'(/A/(5D12.5))')  'Scaling factors:', &
                              (fact(iscale+i-1),i=1,k)
          if (k.lt.n) write (mp,'(16X,A)') ' . . .'
        endif
        k = min(10,n)
        if (ldiag.ge.4) k = n
        do 10 j = 1,nrhs
          write(mp,'(/A,I10)') 'Right-hand side',j
          write (mp,'((1P,5D13.3))') (rhs(i,j),i=1,k)
          if (k.lt.n) write (mp,'(A)') '     . . .'
   10   continue
      end if
      llw = lw/nrhs
      if (icntl(15).eq.1) then
        iscale = lfact-n
        do 5555 i = 1, n
          scale = fact(iscale+i-1)
          if (job.ge.4) scale = one/fact(iscale+i-1)
          do 4444 j = 1, nrhs
            rhs(i,j) = scale*rhs(i,j)
 4444     continue
 5555   continue
      endif
      if (job.le.2) then
        if (nrhs.eq.1) then
          call ma57xd(n,fact,lfact,ifact,lifact,rhs,lrhs, &
                      w,llw,iw1,icntl)
        else
          call ma57qd(n,fact,lfact,ifact,lifact,nrhs,rhs,lrhs, &
                      w,llw,iw1,icntl)
        endif
        if (job.eq.2) go to 15
        if (nrhs.eq.1) then
          call ma57yd(n,fact,lfact,ifact,lifact,rhs,lrhs, &
                      w,llw,iw1,icntl)
        else
          call ma57rd(n,fact,lfact,ifact,lifact,nrhs,rhs,lrhs, &
                      w,llw,iw1,icntl)
        endif
      endif
      if (job.eq.3) &
        call ma57sd(fact,lfact,ifact,lifact,nrhs,rhs,lrhs, &
                    w,llw,icntl)
      if (job.ge.4) &
        call ma57td(n,fact,lfact,ifact,lifact,nrhs,rhs,lrhs, &
                    w,llw,iw1,icntl)
   15 if (icntl(15).eq.1) then
        iscale = lfact-n
        do 6666 i = 1, n
          scale = fact(iscale+i-1)
          if (job.eq.2) scale = one/fact(iscale+i-1)
          do 7777 j = 1, nrhs
            rhs(i,j) = scale*rhs(i,j)
 7777     continue
 6666   continue
      endif
      if (ldiag.ge.3 .and. mp.ge.0) then
        write (mp,'(//A)') &
             'Leaving solution phase (MA57CD) with ...'
        do 20 j = 1,nrhs
          write(mp,'(/A,I10)') 'Solution       ',j
          write (mp,'(1P,5D13.3)') (rhs(i,j),i=1,k)
          if (k.lt.n) write (mp,'(A)') '     . . .'
   20   continue
      endif
  500 return
      end
      subroutine ma57qd(n,fact,lfact,ifact,lifact,nrhs,rhs,lrhs, &
                        w,lw,iw1,icntl)
      integer n,lfact
      double precision fact(lfact)
      integer lifact,ifact(lifact),nrhs,lrhs,lw
      double precision w(lw,nrhs),rhs(lrhs,nrhs)
      integer iw1(n),icntl(20)
      intrinsic abs
      external dgemm,dtpsv
      double precision one
      parameter (one=1.0d0)
      integer apos,i,iblk,ii,ipiv,irhs,iwpos,j,j1,j2,k, &
              ncols,nrows
      double precision w1
      apos = 1
      iwpos = 4
      do 270 iblk = 1,ifact(3)
        iw1(iblk) = iwpos
        ncols = ifact(iwpos)
        nrows = ifact(iwpos+1)
        iwpos = iwpos + 2
        if (nrows.gt.4 .and. ncols.gt.icntl(13)) then
          do 10 i = 1,ncols
            ii = abs(ifact(iwpos+i-1))
            do 11 j = 1,nrhs
              w(i,j) = rhs(ii,j)
   11       continue
   10     continue
          do 12 j = 1,nrhs
            call dtpsv('L','N','U',nrows,fact(apos),w(1,j),1)
   12     continue
          apos = apos + (nrows* (nrows+1))/2
          if (ncols.gt.nrows) call dgemm('N','N',ncols-nrows,nrhs,nrows, &
                                        one,fact(apos),ncols-nrows, &
                                        w,lw,one,w(nrows+1,1),lw)
          apos = apos + nrows* (ncols-nrows)
          do 35 i = 1,ncols
            ii = abs(ifact(iwpos+i-1))
            do 36 j = 1,nrhs
              rhs(ii,j) = w(i,j)
   36       continue
   35     continue
        else
        j1 = iwpos
        j2 = iwpos + nrows - 1
        do 130 ipiv = 1,nrows
          apos = apos + 1
          do 101 ii = 1,nrhs
            w1 = rhs(abs(ifact(j1)),ii)
            k = apos
            do 100 j = j1+1,j2
              irhs = abs(ifact(j))
              rhs(irhs,ii) = rhs(irhs,ii) - fact(k)*w1
              k = k + 1
  100       continue
  101     continue
          apos = k
          j1 = j1 + 1
  130   continue
        j2 = iwpos + ncols - 1
        do 136 ipiv = 1,nrows
          do 135 ii = 1,nrhs
            k = apos
            w1 = rhs(abs(ifact(iwpos+ipiv-1)),ii)
            do 133 j = j1,j2
              irhs = abs(ifact(j))
              rhs(irhs,ii) = rhs(irhs,ii) + w1*fact(k)
              k = k + 1
  133       continue
  135     continue
          apos = k
  136   continue
      end if
      iwpos = iwpos + ncols
  270 continue
      end
      subroutine ma57rd(n,fact,lfact,ifact,lifact,nrhs,rhs,lrhs, &
                        w,lw,iw1,icntl)
      integer n,lfact
      double precision fact(lfact)
      integer lifact,ifact(lifact),nrhs,lrhs,lw
      double precision w(lw,nrhs),rhs(lrhs,nrhs)
      integer iw1(n),icntl(20)
      intrinsic abs
      external dgemm,dtpsv
      double precision one
      parameter (one=1.0d0)
      integer apos,apos2,i,iblk,ii,ipiv,irhs,irhs1, &
              irhs2,iwpos,j,jpiv,j1,j2,k,kk,lrow,ncols,nrows
      double precision w1
      apos = ifact(1)
      apos2 = ifact(2)
      do 380 iblk = ifact(3),1,-1
        iwpos = iw1(iblk)
        ncols = abs(ifact(iwpos))
        nrows = abs(ifact(iwpos+1))
        apos = apos - nrows* (ncols-nrows)
        iwpos = iwpos + 2
        if (nrows.gt.4 .and. ncols.gt.icntl(13)) then
          do 5 i = nrows + 1,ncols
            ii = abs(ifact(iwpos+i-1))
            do 3 j = 1,nrhs
              w(i,j) = rhs(ii,j)
    3       continue
    5     continue
          do 10 ipiv = nrows,1,-1
            irhs = abs(ifact(iwpos+ipiv-1))
            apos = apos - (nrows+1-ipiv)
            do 9 j = 1,nrhs
              w(ipiv,j) = rhs(irhs,j)*fact(apos)
    9       continue
   10     continue
          jpiv = -1
          do 20 ipiv = nrows,1,-1
            irhs = ifact(iwpos+ipiv-1)
            if (irhs.lt.0) then
              irhs1 = -ifact(iwpos+ipiv-1+jpiv)
              do 19 j = 1,nrhs
                w(ipiv,j) = rhs(irhs1,j)*fact(apos2) + w(ipiv,j)
   19         continue
              if (jpiv.eq.1) apos2 = apos2 - 1
              jpiv = -jpiv
            end if
   20     continue
          k = ncols - nrows
          if (k.gt.0) call dgemm('T','N',nrows,nrhs,k,one, &
                                 fact(apos+(nrows*(nrows+1))/2),k, &
                                 w(nrows+1,1),lw,one,w,lw)
          do 22 j = 1,nrhs
            call dtpsv('L','T','U',nrows,fact(apos),w(1,j),1)
   22     continue
          do 60 i = 1,nrows
            ii = abs(ifact(iwpos+i-1))
            do 59 j = 1,nrhs
              rhs(ii,j) = w(i,j)
   59       continue
   60     continue
        else
          j1 = iwpos
          j2 = iwpos + ncols - 1
          jpiv = -1
          do 210 ipiv = nrows,1,-1
            irhs = ifact(iwpos+ipiv-1)
            lrow = nrows + 1 - ipiv
            if (irhs.gt.0) then
              apos = apos - lrow
              do 65 j = 1,nrhs
                rhs(irhs,j) = rhs(irhs,j)*fact(apos)
   65         continue
            else
              if (jpiv.eq.-1) then
                irhs1 = -ifact(iwpos+ipiv-2)
                irhs2 = -irhs
                apos = apos - lrow - lrow - 1
                do 68 j = 1,nrhs
                  w1 = rhs(irhs1,j)*fact(apos) + &
                       rhs(irhs2,j)*fact(apos2)
                  rhs(irhs2,j) = rhs(irhs1,j)*fact(apos2) + &
                                 rhs(irhs2,j)*fact(apos+lrow+1)
                  rhs(irhs1,j) = w1
   68           continue
                apos2 = apos2 - 1
              end if
              jpiv = -jpiv
            end if
  210     continue
          apos = apos + (nrows* (nrows+1))/2
          kk = apos
          j1 = iwpos + nrows
          do 220 ipiv = 1,nrows
            irhs = abs(ifact(iwpos+ipiv-1))
            do 218 ii = 1,nrhs
              w1 = rhs(irhs,ii)
              k = kk
              do 215 j = j1,j2
                w1 = w1 + fact(k)*rhs(abs(ifact(j)),ii)
                k = k + 1
  215         continue
              rhs(irhs,ii) = w1
  218       continue
            kk = k
  220     continue
          j2 = iwpos + nrows - 1
          do 260 ipiv = 1,nrows
            irhs = abs(ifact(j1-1))
            apos = apos - ipiv
            do 240 ii = 1,nrhs
              w1 = rhs(irhs,ii)
              k = apos + 1
              do 230 j = j1,j2
                w1 = w1 - fact(k)*rhs(abs(ifact(j)),ii)
                k = k + 1
  230         continue
              rhs(irhs,ii) = w1
  240       continue
            j1 = j1 - 1
  260     continue
        end if
  380 continue
      end
      subroutine ma57ud(fact,lfact,ifact,lifact,icntl)
      integer lfact
      double precision fact(lfact)
      integer lifact,ifact(lifact),icntl(20)
      intrinsic min,sign
      character*72 line
      integer apos,apos2,iblk,iline,irow,iwpos,j,jpiv,j1,j2,k, &
              ldiag,len,mp,nblk,ncols,nrows
      character*1 pm(-2:2)
      data pm/'*','-','.','+','.'/
      double precision zero,tiny,fd15ad
      parameter (zero=0.0d0)
      external fd15ad
      mp = icntl(3)
      ldiag = icntl(5)
      tiny = fd15ad('T')
      apos2 = ifact(1)
      nblk = ifact(3)
      if (ldiag.eq.3) nblk = min(1,nblk)
      len = 12
      if (ldiag.eq.5) len = 1
      if (len.eq.12) then
        if (nblk.eq.ifact(3)) then
          write (mp,'(/A)') &
            'For each block, the following information is provided:'
        else
          write (mp,'(/A,A)') 'For the first block only,', &
            ' the following information is provided:'
        end if
      end if
      if (len.eq.12) write (mp,'(A)') &
          '   1. Block number, number of rows, number of columns', &
          '   2. List of indices for the pivot, each negated if part of' &
          ,'      a 2x2 pivot', &
          '   3. The factorized block pivot', &
          '      It has the form', &
          '            -1  T', &
          '        L  D   L ', &
          '                         -1    T', &
          '      and is printed as D and L  packed together.', &
          '   4. List of indices for the non-pivot columns', &
          '   5. The non-pivot part as rectangular block by rows'
      iwpos = 4
      apos = 1
      do 300 iblk = 1,nblk
        ncols = ifact(iwpos)
        nrows = ifact(iwpos+1)
        iwpos = iwpos + 2
        write (mp,'(/4(A,I6))') 'Block pivot',iblk,' with',nrows, &
              ' rows and', ncols,' columns'
        if (len.eq.12) write (mp,'(6I12)') &
                             (ifact(k),k=iwpos,iwpos+nrows-1)
        if (len.eq.1) write (mp,'(72A1)') (pm(sign(1,ifact(k))), &
            k=iwpos,iwpos+nrows-1)
        jpiv = 0
        do 30 irow = 1,nrows
          if (jpiv.eq.1) then
            jpiv = 0
          else
            if (ifact(iwpos+irow-1).lt.0) jpiv = 1
          end if
          iline = 1
          do 10 j = 1,irow - 1
            write (line(iline:iline+len-1),'(A)') ' '
            iline = iline + len
            if (iline.gt.72) then
              write (mp,'(A)') line
              iline = 1
            end if
   10     continue
          do 20 j = irow,nrows
            if (len.eq.12) write (line(iline:iline+11), &
                '(1P,D12.4)') fact(apos)
            if (len.eq.1) then
               if (fact(apos).eq.zero) then
                  write (line(iline:iline),'(A)') '.'
               else
                  write (line(iline:iline),'(A)') '*'
               end if
            end if
            apos = apos + 1
            if (j.eq.irow+1) then
              if (jpiv.eq.1) then
                if (len.eq.12) write (line(iline:iline+11), &
                    '(1P,D12.4)') fact(apos2)
                if (len.eq.1) then
                    if (fact(apos2).eq.zero) then
                       write (line(iline:iline),'(A)') '.'
                    else
                       write (line(iline:iline),'(A)') '*'
                    end if
                end if
                apos2 = apos2 + 1
              end if
            end if
            iline = iline + len
            if (iline.gt.72) then
              write (mp,'(A)') line
              iline = 1
            end if
   20     continue
          if (iline.gt.1) then
            line(iline:) = ' '
            write (mp,'(A)') line
          end if
   30   continue
        iwpos = iwpos + nrows
        if (len.eq.12) write (mp,'(6I12)') (ifact(k),k=iwpos, &
            iwpos+ncols-nrows-1)
        if (len.eq.1) write (mp,'(72A1)') (pm(sign(1,ifact(k))), &
            k=iwpos,iwpos+ncols-nrows-1)
        iwpos = iwpos + ncols - nrows
        do 280 irow = 1,nrows
          j1 = nrows
          j2 = ncols
          iline = 1
          do 110 j = j1 + 1,j2
            if (len.eq.12) write (line(iline:iline+11), &
                '(1P,D12.4)') fact(apos)
            if (len.eq.1) then
               if (fact(apos).eq.zero) then
                  write (line(iline:iline),'(A)') '.'
               else
                  write (line(iline:iline),'(A)') '*'
               end if
            end if
            apos = apos + 1
            iline = iline + len
            if (iline.gt.72) then
              write (mp,'(A)') line
              iline = 1
            end if
  110     continue
          if (iline.gt.1) then
            line(iline:) = ' '
            write (mp,'(A)') line
          end if
  280   continue
  300 continue
      end
      subroutine ma57sd(fact,lfact,ifact,lifact,nrhs,rhs,lrhs, &
                        w,lw,icntl)
      integer lfact
      double precision fact(lfact)
      integer lifact,ifact(lifact),nrhs,lrhs,lw
      double precision w(lw,nrhs),rhs(lrhs,nrhs)
      integer icntl(20)
      intrinsic abs
      external dgemm,dtpsv
      integer apos,apos2,i,iblk,ii,ipiv,irhs,irhs1, &
              irhs2,iwpos,j,jpiv,ncols,nrows
      double precision w1
      apos = 1
      apos2 = ifact(1)
      iwpos = 4
      do 380 iblk = 1,ifact(3)
        ncols = ifact(iwpos)
        nrows = ifact(iwpos+1)
        iwpos = iwpos + 2
        if (nrows.gt.4 .and. ncols.gt.icntl(13)) then
          do 10 ipiv = 1,nrows
            irhs = abs(ifact(iwpos+ipiv-1))
            do 9 j = 1,nrhs
              w(ipiv,j) = rhs(irhs,j)*fact(apos)
    9       continue
            apos = apos + (nrows+1-ipiv)
   10     continue
          jpiv = 1
          do 20 ipiv = 1,nrows
            irhs = ifact(iwpos+ipiv-1)
            if (irhs.lt.0) then
              irhs1 = -ifact(iwpos+ipiv-1+jpiv)
              do 19 j = 1,nrhs
                w(ipiv,j) = rhs(irhs1,j)*fact(apos2) + w(ipiv,j)
   19         continue
              if (jpiv.eq.-1) apos2 = apos2 + 1
              jpiv = -jpiv
            end if
   20     continue
          do 60 i = 1,nrows
            ii = abs(ifact(iwpos+i-1))
            do 59 j = 1,nrhs
              rhs(ii,j) = w(i,j)
   59       continue
   60     continue
        else
          jpiv = 1
          do 210 ipiv = 1,nrows
            irhs = ifact(iwpos+ipiv-1)
            if (irhs.gt.0) then
              do 65 j = 1,nrhs
                rhs(irhs,j) = rhs(irhs,j)*fact(apos)
   65         continue
              apos = apos + nrows - ipiv + 1
            else
              if (jpiv.eq.1) then
                irhs1 = -irhs
                irhs2 = -ifact(iwpos+ipiv)
                do 68 j = 1,nrhs
                  w1 = rhs(irhs1,j)*fact(apos) + &
                       rhs(irhs2,j)*fact(apos2)
                  rhs(irhs2,j) = rhs(irhs1,j)*fact(apos2) + &
                                 rhs(irhs2,j)*fact(apos+nrows-ipiv+1)
                  rhs(irhs1,j) = w1
   68           continue
                apos2 = apos2 + 1
              end if
              jpiv = -jpiv
              apos = apos + nrows - ipiv + 1
            end if
  210     continue
        end if
        iwpos = iwpos + ncols
        apos = apos + nrows*(ncols-nrows)
  380 continue
      end
      subroutine ma57td(n,fact,lfact,ifact,lifact,nrhs,rhs,lrhs, &
                        w,lw,iw1,icntl)
      integer n,lfact
      double precision fact(lfact)
      integer lifact,ifact(lifact),nrhs,lrhs,lw
      double precision w(lw,nrhs),rhs(lrhs,nrhs)
      integer iw1(n),icntl(20)
      intrinsic abs
      external dgemm,dtpsv
      double precision one
      parameter (one=1.0d0)
      integer apos,i,iblk,ii,ipiv,irhs, &
              iwpos,j,j1,j2,k,kk,ncols,nrows
      double precision w1
      apos = ifact(1)
      iwpos = 4
      do 10 i = 1,ifact(3)-1
        iw1(i) = iwpos
        iwpos = iwpos + abs(ifact(iwpos))+2
   10 continue
      iw1(ifact(3)) = iwpos
      do 380 iblk = ifact(3),1,-1
        iwpos = iw1(iblk)
        ncols = abs(ifact(iwpos))
        nrows = abs(ifact(iwpos+1))
        apos = apos - nrows* (ncols-nrows)
        iwpos = iwpos + 2
        if (nrows.gt.4 .and. ncols.gt.icntl(13)) then
          do 5 i = 1,ncols
            ii = abs(ifact(iwpos+i-1))
            do 3 j = 1,nrhs
              w(i,j) = rhs(ii,j)
    3       continue
    5     continue
          k = ncols - nrows
          if (k.gt.0) call dgemm('T','N',nrows,nrhs,k,one, &
                                 fact(apos),k, &
                                 w(nrows+1,1),lw,one,w,lw)
          apos = apos-(nrows*(nrows+1))/2
          do 22 j = 1,nrhs
            call dtpsv('L','T','U',nrows,fact(apos),w(1,j),1)
   22     continue
          do 60 i = 1,nrows
            ii = abs(ifact(iwpos+i-1))
            do 59 j = 1,nrhs
              rhs(ii,j) = w(i,j)
   59       continue
   60     continue
        else
          j1 = iwpos
          j2 = iwpos + ncols - 1
          kk = apos
          j1 = iwpos + nrows
          do 220 ipiv = 1,nrows
            irhs = abs(ifact(iwpos+ipiv-1))
            do 218 ii = 1,nrhs
              w1 = rhs(irhs,ii)
              k = kk
              do 215 j = j1,j2
                w1 = w1 + fact(k)*rhs(abs(ifact(j)),ii)
                k = k + 1
  215         continue
              rhs(irhs,ii) = w1
  218       continue
            kk = k
  220     continue
          j2 = iwpos + nrows - 1
          do 260 ipiv = 1,nrows
            irhs = abs(ifact(j1-1))
            apos = apos - ipiv
            do 240 ii = 1,nrhs
              w1 = rhs(irhs,ii)
              k = apos + 1
              do 230 j = j1,j2
                w1 = w1 - fact(k)*rhs(abs(ifact(j)),ii)
                k = k + 1
  230         continue
              rhs(irhs,ii) = w1
  240       continue
            j1 = j1 - 1
  260     continue
        end if
  380 continue
      end
      subroutine ma57dd(job,n,ne,a,irn,jcn,fact,lfact,ifact,lifact, &
                        rhs,x,resid,w,iw,icntl,cntl,info,rinfo)
      integer job,n,ne
      double precision a(ne)
      integer irn(ne),jcn(ne),lfact
      double precision fact(lfact)
      integer lifact,ifact(lifact)
      double precision rhs(n),x(n),resid(n),w(n,*)
      integer iw(n),icntl(20)
      double precision cntl(5)
      integer info(40)
      double precision rinfo(20)
      double precision zero,one
      parameter (zero=0.d0,one=1.0d0)
      double precision cond(2),ctau,dxmax,error,oldomg(2),oldom2, &
                       omega(2),om2,tau
      integer i,icntlc(20),iter,j,k,kase,kk,ldiag,lp,mp,keep71(5)
      logical lcond(2)
      intrinsic min
      external ma57cd,ma57ud,fd15ad,mc71ad
      double precision eps,fd15ad
      character*80 outbuf
      lp = icntl(1)
      mp = icntl(3)
      ldiag = icntl(5)
      info(1) = 0
      if (n.le.0) then
        info(1) = -1
        if (ldiag.gt.0 .and. lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           write (outbuf,'(A,I3)') &
            '**** Error return from MA57DD ****  INFO(1) =',info(1)
           call wrtlin( outbuf )
           write (outbuf,'(A,I12)') 'N has value ',n
           call wrtlin( outbuf )
        end if
        info(2) = n
        goto 500
      endif
      if (ne.lt.0) then
        info(1) = -2
        if (ldiag.gt.0 .and. lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           write (outbuf,'(A,I3)') &
            '**** Error return from MA57DD ****  INFO(1) =',info(1)
           call wrtlin( outbuf )
           write (outbuf,'(A,I12)') 'NE has value ',ne
           call wrtlin( outbuf )
        end if
        info(2) = ne
        goto 500
      endif
      if (icntl(9).lt.1) then
        info(1) = -13
        if (ldiag.gt.0 .and. lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           write (outbuf,'(A,I4)') &
            '**** Error return from MA57DD ****  INFO(1) =',info(1)
           call wrtlin( outbuf )
           write (outbuf,'(A,I12)') 'ICNTL(9) has value ',icntl(9)
           call wrtlin( outbuf )
        end if
        info(2) = icntl(9)
        goto 500
      endif
      if (job.lt.0 .or. job.gt.4 .or. (icntl(9).gt.1 .and. &
          (job.ne.0 .and. job.ne.2)))  then
        info(1) = -12
        info(2) = job
        if (ldiag.gt.0 .and. lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           write (outbuf,'(A,I4)') &
            '**** Error return from MA57DD ****  INFO(1) =',info(1)
           call wrtlin( outbuf )
           write (outbuf,'(A,I12)') 'JOB has value ',job
           call wrtlin( outbuf )
        end if
        if (icntl(9).gt.1 .and. ldiag.gt.0 .and. lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           write (outbuf,'(A,I3)') 'and ICNTL(9) =',icntl(9)
           call wrtlin( outbuf )
        end if
        goto 500
      endif
      if (ne.eq.0) then
        if (job.ne.3) then
          do 8 i = 1,n
            resid(i) = zero
  8       continue
        endif
        do 9 i = 1,n
          x(i) = zero
  9     continue
        info(30)=0
        do 10 i = 6,13
          rinfo(i) = zero
 10     continue
        go to 500
      endif
      if (ldiag.ge.3 .and. mp.ge.0) then
        write (mp,99999) job,n,ne,(icntl(i),i=1,5),lfact,lifact, &
         icntl(9),icntl(10),icntl(13),cntl(3)
99999 format(/'Entering iterative refinement solution phase ', &
        '(MA57DD) with ...'/ &
        'JOB       Control for coefficient matrix      =',i12/ &
        'N         Order of matrix                     =',i12/ &
        'NE        Number of entries in matrix         =',i12/ &
        'ICNTL(1)  Stream for errors                   =',i12/ &
        ' --- (2)  Stream for warnings                 =',i12/ &
        ' --- (3)  Stream for monitoring               =',i12/ &
        ' --- (4)  Stream for statistics               =',i12/ &
        ' --- (5)  Level of diagnostic printing        =',i12/ &
        'LFACT     Length of array FACT                =',i12/ &
        'LIFACT    Length of array IFACT               =',i12/ &
        'ICNTL(9)  Number steps iterative refinement   =',i12/ &
        'ICNTL(10) Control for error analysis          =',i12/ &
        'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',i12/ &
        'CNTL(3)   Convergence test for IR             =',1p,d12.4)
        call ma57ud(fact,lfact,ifact,lifact,icntl)
        k = min(10,n)
        if (ldiag.ge.4) k = n
        write(mp,'(/A)') 'Right-hand side'
        write (mp,'((4X, 1P,5D13.3))') (rhs(i),i=1,k)
        if (k.lt.n) write (mp,'(A)') '     . . .'
      end if
      do 15 i=1,5
        icntlc(i) = icntl(i)
   15 continue
      icntlc(13) = icntl(13)
      icntlc(15) = icntl(15)
      icntlc(3) = -1
      if (job.le.2) then
        if (job .le. 1) then
          do 14 i = 1,n
            x(i) = rhs(i)
            resid(i) = rhs(i)
   14     continue
          call ma57cd(1,n,fact,lfact,ifact,lifact,1,x,n,w,n,iw, &
                      icntlc,info)
        else
          do 13 i = 1,n
            resid(i) = rhs(i)
   13     continue
        endif
        if (icntl(9).eq.1) then
          do 16 kk = 1,ne
            i = irn(kk)
            j = jcn(kk)
            if (min(i,j).lt.1 .or. max(i,j).gt.n) go to 16
            resid(j) = resid(j) - a(kk)*x(i)
            if (i.ne.j) resid(i) = resid(i) - a(kk)*x(j)
   16     continue
          if (job.eq.0) go to 340
        else
          do 18 i = 1,n
            w(i,1) = zero
            w(i,3) = zero
   18     continue
          do 17 kk = 1,ne
            i = irn(kk)
            j = jcn(kk)
            if (min(i,j).lt.1 .or. max(i,j).gt.n) go to 17
            resid(j) = resid(j) - a(kk)*x(i)
            w(j,1) = w(j,1) + abs(a(kk)*x(i))
            w(j,3) = w(j,3) + abs(a(kk))
            if (i.ne.j) then
              resid(i) = resid(i) - a(kk)*x(j)
              w(i,1) = w(i,1) + abs(a(kk)*x(j))
              w(i,3) = w(i,3) + abs(a(kk))
            endif
   17     continue
        dxmax = zero
        do 221 i = 1,n
          dxmax = max(dxmax,abs(x(i)))
  221   continue
      eps = fd15ad('E')
        ctau = 1000.*eps
          omega(1) = zero
          omega(2) = zero
          do 231 i = 1,n
            tau = (w(i,3)*dxmax+abs(rhs(i)))*n*ctau
            if ((w(i,1)+abs(rhs(i))).gt.tau) then
              omega(1) = max(omega(1),abs(resid(i))/ &
                         (w(i,1)+abs(rhs(i))))
              iw(i) = 1
            else
              if (tau.gt.zero) then
                omega(2) = max(omega(2),abs(resid(i))/ &
                           (w(i,1)+w(i,3)*dxmax))
              end if
              iw(i) = 2
            end if
  231     continue
          om2 = omega(1) + omega(2)
          iter = 0
          if (om2.le.eps) then
            go to 270
          endif
          do 251 i = 1,n
            w(i,2) = x(i)
  251     continue
          oldomg(1) = omega(1)
          oldomg(2) = omega(2)
          oldom2 = om2
        endif
      endif
      do 260 iter = 1,icntl(9)
        call ma57cd(1,n,fact,lfact,ifact,lifact,1,resid,n,w,n,iw, &
                    icntlc,info)
        do 141 i = 1,n
          x(i) = x(i) + resid(i)
  141   continue
        if (job.lt.4 .and. icntl(9).eq.1) go to 340
        if (icntl(9).eq.1) then
          do 151 i = 1,n
            resid(i) = rhs(i)
  151     continue
          do 181 kk = 1,ne
            i = irn(kk)
            j = jcn(kk)
            if (min(i,j).lt.1 .or. max(i,j).gt.n) go to 181
            resid(j) = resid(j) - a(kk)*x(i)
            if (i.ne.j) resid(i) = resid(i) - a(kk)*x(j)
  181     continue
          go to 340
        else
          do 153 i = 1,n
            resid(i) = rhs(i)
            w(i,1) = zero
  153     continue
          do 183 kk = 1,ne
            i = irn(kk)
            j = jcn(kk)
            if (min(i,j).lt.1 .or. max(i,j).gt.n) go to 183
            resid(j) = resid(j) - a(kk)*x(i)
            w(j,1) = w(j,1) + abs(a(kk)*x(i))
            if (i.ne.j) then
              resid(i) = resid(i) - a(kk)*x(j)
              w(i,1) = w(i,1) + abs(a(kk)*x(j))
            endif
  183     continue
        endif
        dxmax = zero
        do 220 i = 1,n
          dxmax = max(dxmax,abs(x(i)))
  220   continue
        omega(1) = zero
        omega(2) = zero
        do 230 i = 1,n
          tau = (w(i,3)*dxmax+abs(rhs(i)))*n*ctau
          if ((w(i,1)+abs(rhs(i))).gt.tau) then
            omega(1) = max(omega(1),abs(resid(i))/ &
                       (w(i,1)+abs(rhs(i))))
            iw(i) = 1
          else
            if (tau.gt.zero) then
              omega(2) = max(omega(2),abs(resid(i))/ &
                         (w(i,1)+w(i,3)*dxmax))
            end if
            iw(i) = 2
          end if
  230   continue
        om2 = omega(1) + omega(2)
        if ((om2+one).le.one) then
          go to 270
        endif
        if (om2.gt.oldom2*cntl(3)) then
          if (om2.gt.oldom2) then
            omega(1) = oldomg(1)
            omega(2) = oldomg(2)
            do 240 i = 1,n
              x(i) = w(i,2)
  240       continue
          end if
          go to 270
        else
          do 250 i = 1,n
            w(i,2) = x(i)
  250     continue
          oldomg(1) = omega(1)
          oldomg(2) = omega(2)
          oldom2 = om2
        end if
  260 continue
      info(1) = -8
      if (lp.ge.0 .and. ldiag.ge.1) then
         call wrtlin( '    *** HSL ERROR ***' )
         call wrtlin( 'Error return from MA57D/DD because of nonconvergence'// &
                      ' of iterative refinement' )
         write( outbuf, '(a,i2,a,i10)' ) 'Error INFO(1) = ',info(1), &
                                         '  with ICNTL(9) = ',icntl(9)
         call wrtlin( outbuf )
      end if
  270 rinfo(6)  = omega(1)
      rinfo(7)  = omega(2)
      rinfo(8) = zero
      do 271 i=1,n
        rinfo(8) = max(rinfo(8),w(i,3))
  271 continue
      rinfo(9) = dxmax
      rinfo(10) = zero
      do 272 i=1,n
        rinfo(10) = max(rinfo(10),abs(resid(i)))
  272 continue
      if (rinfo(8)*rinfo(9).ne.zero) &
      rinfo(10) = rinfo(10)/(rinfo(8)*rinfo(9))
      info(30) = iter
      if (info(1).lt.0) go to 340
      if (icntl(10).le.0) go to 340
      lcond(1) = .false.
      lcond(2) = .false.
      error    = zero
      do 280 i = 1,n
        if (iw(i).eq.1) then
          w(i,1) = w(i,1) + abs(rhs(i))
          w(i,2) = zero
          lcond(1) = .true.
        else
          w(i,2) = w(i,1) + w(i,3)*dxmax
          w(i,1) = zero
          lcond(2) = .true.
        end if
  280 continue
      do 330 k = 1,2
        if (lcond(k)) then
          kase = 0
          do 310 kk = 1,40
            call mc71ad(n,kase,w(1,3),cond(k),w(1,4),iw,keep71)
            if (kase.eq.0) go to 320
            if (kase.eq.1) then
              call ma57cd(1,n,fact,lfact,ifact,lifact,1,w(1,3), &
                          n,w(1,4),n,iw,icntlc,info)
              do 290 i = 1,n
                w(i,3) = w(i,k)*w(i,3)
  290         continue
            end if
            if (kase.eq.2) then
              do 300 i = 1,n
                w(i,3) = w(i,k)*w(i,3)
  300         continue
              call ma57cd(1,n,fact,lfact,ifact,lifact,1,w(1,3),n, &
                          w(1,4),n,iw,icntlc,info)
            end if
  310     continue
          info(1) = -14
          if (lp.ge.0 .and. ldiag.ge.1) then
             call wrtlin( '    *** HSL ERROR ***' )
             call wrtlin( 'Error return from MA57D/DD because of '// &
                          'error in MC71A/AD' )
             call wrtlin( 'Error not calculated' )
          end if
  320     if (dxmax.gt.zero) cond(k) = cond(k)/dxmax
          error = error + omega(k)*cond(k)
        else
          cond(k) = zero
        endif
  330 continue
      rinfo(11)  = cond(1)
      rinfo(12)  = cond(2)
      rinfo(13)  = error
 340  if (ldiag.ge.3 .and. mp.ge.0) then
        write(mp,99980) info(1)
99980 format (/'Leaving iterative refinement solution phase ', &
        '(MA57DD) with ...'/ &
            'INFO (1)                                      =',i12/)
        if (info(1).lt.0) go to 500
        if (icntl(9).gt.1) then
          write(mp,99981) info(30),(rinfo(i),i=6,10)
99981     format( &
           'INFO(30)  Number steps iterative ref   =',i10/ &
           'RINFO(6)  Backward errors  (OMEGA(1))  =',1pd10.3/ &
           '-----(7)  Backward errors  (OMEGA(2))  =',1pd10.3/ &
           '-----(8)  Infinity norm of matrix      =',1pd10.3/ &
           '-----(9)  Infinity norm of solution    =',1pd10.3/ &
           '-----(10) Norm of scaled residuals     =',1pd10.3)
          if (icntl(10).gt.0) write(mp,99979) (rinfo(i),i=11,13)
99979       format ( &
             'RINFO(11) Condition number (COND(1))   =',1pd10.3/ &
             'RINFO(12) Condition number (COND(2))   =',1pd10.3/ &
             'RINFO(13) Error in solution            =',1pd10.3)
          write(mp,'(/A,I10)') 'Residual'
          k=min(n,10)
          if (ldiag.ge.4) k = n
          write (mp,'(1P,5D13.3)') (resid(i),i=1,k)
          if (k.lt.n) write (mp,'(A)') '     . . .'
        else
          if (job.ge.1 .and. job.le.3) then
            write(mp,'(/A,I10)') 'Correction to solution'
          else
            write(mp,'(/A,I10)') 'Residual'
          endif
          k=min(n,10)
          if (ldiag.ge.4) k = n
          write (mp,'(1P,5D13.3)') (resid(i),i=1,k)
          if (k.lt.n) write (mp,'(A)') '     . . .'
        end if
        k=min(n,10)
        if (ldiag.ge.4) k = n
        write(mp,'(/A,I10)') 'Solution'
        write (mp,'(1P,5D13.3)') (x(i),i=1,k)
        if (k.lt.n) write (mp,'(A)') '     . . .'
      end if
 500  return
      end
      subroutine ma57ed(n,ic,keep,fact,lfact,newfac,lnew, &
                        ifact,lifact,newifc,linew,info)
      integer n,ic,keep(*),lfact,lnew,lifact,linew,info(40)
      double precision fact(lfact),newfac(lnew)
      integer ifact(lifact),newifc(linew)
      integer aposbb,astk,hold,i,istk,iwpos,move,nfront
      hold = n + 3
      info(1) = 0
      info(2) = 0
      if (ic.ge.1) then
        if (linew.le.lifact) then
          info(1) = -7
          info(2) = linew
          return
        endif
        iwpos = keep(hold+7)
        istk  = keep(hold+14)
        nfront = keep(hold+23)
        do 10 i = 1,iwpos+nfront-1
          newifc(i) = ifact(i)
   10   continue
        move = linew - lifact
        do 20 i = istk+1,lifact
          newifc(i+move) = ifact(i)
   20   continue
          keep(hold+13) = keep(hold+13) + move
          keep(hold+14) = istk + move
          keep(hold+18) = keep(hold+18) + move
      endif
      if (ic.ne.1) then
        if (lnew.le.lfact) then
          info(1) = -7
          info(2) = lnew
          return
        endif
        aposbb = keep(hold+9)
        astk   = keep(hold+15)
        do 60 i = 1, aposbb-1
          newfac(i) = fact(i)
   60   continue
        move = lnew - lfact
        do 70 i = astk+1,lfact
          newfac(i+move) = fact(i)
   70   continue
        keep(hold+12) = keep(hold+12) + move
        keep(hold+15) = astk + move
        keep(hold+19) = keep(hold+19) + move
      endif
      return
      end
      subroutine ma57gd(n,ne,irn,jcn,iw,ipe,count,flag,iwfr,icntl,info)
      integer n,ne,irn(ne),jcn(ne),iw(ne*2+n),ipe(n),count(n), &
              flag(n),iwfr,icntl(20),info(40)
      intrinsic max,min
      integer i,j,k,l,ldiag,wp
      character*80 outbuf
      wp = icntl(2)
      ldiag = icntl(5)
      if (wp.lt.0) ldiag = 0
      info(1) = 0
      info(3) = 0
      do 10 i = 1,n
        flag(i) = 0
        count(i) = 0
   10 continue
      do 20 k = 1,ne
        i = irn(k)
        j = jcn(k)
        if (min(i,j).lt.1 .or. max(i,j).gt.n) then
          info(3) = info(3) + 1
          info(1) = 1
          if (info(3).eq.1 .and. ldiag.gt.1) then
             call wrtlin( '    *** HSL WARNING ***' )
             write (outbuf,'(A,I2)') &
              '*** Warning message from subroutine MA57AD *** INFO(1) =',info(1)
             call wrtlin( outbuf )
          end if
          if (info(3).le.10 .and. ldiag.gt.1) then
             call wrtlin( '    *** HSL WARNING ***' )
             write (outbuf,'(3(I10,A))') &
              k, 'th entry (in row', i, ' and column', j, ') ignored'
             call wrtlin( outbuf )
          end if
        else if (i.ne.j) then
          count(i) = count(i) + 1
          count(j) = count(j) + 1
        end if
   20 continue
      ipe(1) = count(1)+1
      do 30 i = 2,n
        ipe(i) = ipe(i-1) + count(i)
   30 continue
      do 40 k = 1,ne
        i = irn(k)
        j = jcn(k)
        if (min(i,j).lt.1 .or. max(i,j).gt.n .or. i.eq.j) go to 40
        ipe(i) = ipe(i) - 1
        iw(ipe(i)) = j
        ipe(j) = ipe(j) - 1
        iw(ipe(j)) = i
   40 continue
      info(4) = 0
      iwfr = 1
      do 60 i = 1,n
        l = ipe(i)
        ipe(i) = iwfr
        do 50 k = l,l+count(i)-1
          j = iw(k)
          if (flag(j).ne.i) then
            flag(j) = i
            iw(iwfr) = j
            iwfr = iwfr + 1
          else
            if (i.lt.j) info(4) = info(4) + 1
          end if
   50   continue
        count(i) = iwfr - ipe(i)
   60 continue
      if (info(4).gt.0) then
        info(1) = info(1) + 2
        if (ldiag.gt.1 .and. wp.ge.0) then
           call wrtlin( '    *** HSL WARNING ***' )
           call wrtlin( '*** Warning message from subroutine MA57AD ***' )
           write( outbuf, '(i10,a)' ) info(4), &
                                      ' off-diagonal duplicate entries found'
           call wrtlin( outbuf )
        end if
      end if
      end
      subroutine ma57jd(n,ne,irn,jcn,perm,iw,ipe,count,flag,iwfr, &
                        icntl,info)
      integer n,ne,irn(ne),jcn(ne),iw(ne+n),ipe(n),count(n), &
              perm(n),flag(n),iwfr,icntl(20),info(40)
      intrinsic max,min
      integer i,j,k,l,ldiag,wp
      character*80 outbuf
      wp = icntl(2)
      ldiag = icntl(5)
      if (wp.lt.0) ldiag = 0
      info(1) = 0
      do 10 i = 1,n
        flag(i) = 0
        count(i) = 0
   10 continue
      info(3) = 0
      do 30 k = 1,ne
        i = irn(k)
        j = jcn(k)
        if (min(i,j).lt.1 .or. max(i,j).gt.n) then
          irn(k) = 0
          jcn(k) = 0
          info(3) = info(3) + 1
          info(1) = 1
          if (info(3).eq.1 .and. ldiag.gt.1) then
             call wrtlin( '    *** HSL WARNING ***' )
             write (outbuf,'(2A,I2)') &
              '*** Warning message from subroutine MA57AD ***', &
              ' INFO(1) =',info(1)
             call wrtlin( outbuf )
          end if
          if (info(3).le.10 .and. ldiag.gt.1) then
             call wrtlin( '    *** HSL WARNING ***' )
             write (outbuf,'(3(I10,A))') &
              k,'th entry (in row',i,' and column',j,') ignored'
             call wrtlin( outbuf )
          end if
        else if (perm(i).le.perm(j)) then
          count(i) = count(i) + 1
        else
          count(j) = count(j) + 1
        end if
   30 continue
      ipe(1) = count(1) + 1
      do 40 i = 2,n
        ipe(i) = ipe(i-1) + count(i) + 1
   40 continue
      do 50 k = 1,ne
        i = irn(k)
        j = jcn(k)
        if (min(i,j).lt.1 .or. max(i,j).gt.n) go to 50
        if (perm(i).le.perm(j)) then
          iw(ipe(i)) = j
          ipe(i) = ipe(i) - 1
        else
          iw(ipe(j)) = i
          ipe(j) = ipe(j) - 1
        end if
   50 continue
      iwfr = 1
      info(4) = 0
      do 70 i = 1,n
        l = ipe(i)
        ipe(i) = iwfr
        do 60 k = l + 1,l + count(i)
          j = iw(k)
          if (flag(j).ne.i) then
            flag(j) = i
            iwfr = iwfr + 1
            iw(iwfr) = j
          else
            if (i.lt.j) info(4) = info(4) + 1
          end if
   60   continue
        if (iwfr.gt.ipe(i)) then
          iw(ipe(i)) = iwfr - ipe(i)
          iwfr = iwfr + 1
        else
          ipe(i) = 0
        end if
   70 continue
      if (info(4).gt.0) then
        info(1) = info(1) + 2
        if (ldiag.gt.1 .and. wp.ge.0) then
           call wrtlin( '    *** HSL WARNING ***' )
           call wrtlin( '*** Warning message from subroutine MA57AD ***' )
           write ( outbuf,'(I10,A)') &
            info(4),' off-diagonal duplicate entries found'
           call wrtlin( outbuf )
        end if
      end if
      end
      subroutine ma57kd(n, ipe, iw, lw, iwfr, perm, ips, nv, flag, &
                        ncmpa)
      integer n,lw,iwfr,ncmpa
      integer ipe(n)
      integer iw(lw), perm(n), ips(n), nv(n), flag(n)
      integer i,j,ml,ms,me,ip,minjs,ie,kdummy,jp
      integer ln,jp1,js,lwfr,jp2,je
      external ma57fd
      do 10 i=1,n
        flag(i) = 0
        nv(i) = 0
        j = perm(i)
        ips(j) = i
   10 continue
      ncmpa = 0
      do 100 ml=1,n
        ms = ips(ml)
        me = ms
        flag(ms) = me
        ip = iwfr
        minjs = n
        ie = me
        do 70 kdummy=1,n
          jp = ipe(ie)
          ln = 0
          if (jp.le.0) go to 60
          ln = iw(jp)
          do 50 jp1=1,ln
            jp = jp + 1
            js = iw(jp)
            if (flag(js).eq.me) go to 50
            flag(js) = me
            if (iwfr.lt.lw) go to 40
            ipe(ie) = jp
            iw(jp) = ln - jp1
            call ma57fd(n, ipe, iw, ip-1, lwfr, ncmpa)
            jp2 = iwfr - 1
            iwfr = lwfr
            if (ip.gt.jp2) go to 30
            do 20 jp=ip,jp2
              iw(iwfr) = iw(jp)
              iwfr = iwfr + 1
   20       continue
   30       ip = lwfr
            jp = ipe(ie)
   40       iw(iwfr) = js
            minjs = min0(minjs,perm(js)+0)
            iwfr = iwfr + 1
   50     continue
   60     ipe(ie) = -me
          je = nv(ie)
          nv(ie) = ln + 1
          ie = je
          if (ie.eq.0) go to 80
   70   continue
   80   if (iwfr.gt.ip) go to 90
        ipe(me) = 0
        nv(me) = 1
        go to 100
   90   minjs = ips(minjs)
        nv(me) = nv(minjs)
        nv(minjs) = me
        iw(iwfr) = iw(ip)
        iw(ip) = iwfr - ip
        ipe(me) = ip
        iwfr = iwfr + 1
  100 continue
      return
      end
!** end of MA57KD**
      subroutine ma57fd(n, ipe, iw, lw, iwfr, ncmpa)
      integer n,lw,iwfr,ncmpa
      integer ipe(n)
      integer   iw(lw)
      integer i,k1,lwfr,ir,k,k2
      ncmpa = ncmpa + 1
      do 10 i=1,n
        k1 = ipe(i)
        if (k1.le.0) go to 10
        ipe(i) = iw(k1)
        iw(k1) = -i
   10 continue
      iwfr = 1
      lwfr = iwfr
      do 60 ir=1,n
        if (lwfr.gt.lw) go to 70
        do 20 k=lwfr,lw
          if (iw(k).lt.0) go to 30
   20   continue
        go to 70
   30   i = -iw(k)
        iw(iwfr) = ipe(i)
        ipe(i) = iwfr
        k1 = k + 1
        k2 = k + iw(iwfr)
        iwfr = iwfr + 1
        if (k1.gt.k2) go to 50
        do 40 k=k1,k2
          iw(iwfr) = iw(k)
          iwfr = iwfr + 1
   40   continue
   50   lwfr = k2 + 1
   60 continue
   70 return
      end
!--------------------------------------------------------------------
!-             Copyright CCLRC Rutherford Appleton Laboratory
!--------------------------------------------------------------------
      subroutine ma57ld(n, ipe, nv, ips, ne, na, node, perm, nsteps, &
                        fils, frere, nd, nemin, subord)
      integer n, nsteps
      integer nd(n)
      integer ipe(n), fils(n), frere(n), subord(n)
      integer nv(n), ips(n), ne(n), na(n), node(n), perm(n)
      integer nemin
      integer i,if,is,nr,nr1,ins,inl,inb,inf,infs,insw
      integer k,l,ison,in,ifson,ino
      integer inos,ib,il,int
      integer iperm
      do 10 i=1,n
        ips(i) = 0
        ne(i) = 0
        node(i) = 0
        subord(i) = 0
   10 continue
      nr = n + 1
      do 50 i=1,n
        if = -ipe(i)
        if (nv(i).eq.0) then
          if (subord(if).ne.0) subord(i) = subord(if)
          subord(if) = i
          node(if) = node(if)+1
        else
          if (if.ne.0) then
            is = -ips(if)
            if (is.gt.0) ipe(i) = is
            ips(if) = -i
          else
            nr = nr - 1
            ne(nr) = i
          endif
        endif
   50 continue
      do 999 i=1,n
       fils(i) = ips(i)
 999  continue
      nr1 = nr
      ins = 0
 1000 if (nr1.gt.n) go to 1151
      ins = ne(nr1)
      nr1 = nr1 + 1
 1070 inl = fils(ins)
      if (inl.lt.0) then
       ins = -inl
       go to 1070
      endif
 1080 if (ipe(ins).lt.0) then
       ins       = -ipe(ins)
       fils(ins) = 0
       go to 1080
      endif
      if (ipe(ins).eq.0) then
       ins = 0
       go to 1000
      endif
      inb = ipe(ins)
!?? I think this test is the wrong way round
      if (nv(inb).ge.nv(ins)) then
!?? So reversed
       ins = inb
       go to 1070
      endif
      inf = inb
 1090 inf = ipe(inf)
      if (inf.gt.0) go to 1090
      inf  = -inf
      infs = -fils(inf)
      if (infs.eq.ins) then
        fils(inf) = -inb
        ips(inf)  = -inb
        ipe(ins)  = ipe(inb)
        ipe(inb)  = ins
      else
        insw = infs
 1100   infs = ipe(insw)
        if (infs.ne.ins) then
          insw = infs
          go to 1100
        endif
        ipe(ins) = ipe(inb)
        ipe(inb) = ins
        ipe(insw)= inb
      endif
        ins      = inb
        go to 1070
 1151 do 51 i=1,n
       frere(i) = ipe(i)
       fils(i) = ips(i)
 51   continue
      is = 1
      i = 0
      iperm = 1
      do 160 k=1,n
        if (i.gt.0) go to 60
        if (nr.gt.n) go to 161
        i = ne(nr)
        ne(nr) = 0
        nr = nr + 1
        il = n
        na(n) = 0
   60   continue
        do 70 l=1,n
          if (ips(i).ge.0) go to 80
          ison = -ips(i)
          ips(i) = 0
          i = ison
          il = il - 1
          na(il) = 0
   70   continue
   80   continue
!?? Do we want to expand for subordinate variables
        ips(i) = k
        ne(is) = ne(is) + node(i) + 1
        if (il.lt.n) na(il+1) = na(il+1) + 1
        na(is) = na(il)
        nd(is) = nv(i)
        node(i) = is
        perm(i) = iperm
        iperm = iperm + 1
        in = i
  777   if (subord(in).eq.0) go to 778
          in = subord(in)
          node(in) = is
          perm(in) = iperm
          iperm = iperm + 1
          go to 777
  778   if (na(is).ne.1) go to 90
        if (nd(is-1)-ne(is-1).eq.nd(is)) go to 100
   90   if (ne(is).ge.nemin) go to 110
        if (na(is).eq.0) go to 110
        if (ne(is-1).ge.nemin) go to 110
  100   na(is-1) = na(is-1) + na(is) - 1
        nd(is-1) = nd(is) + ne(is-1)
        ne(is-1) = ne(is) + ne(is-1)
        ne(is) = 0
        node(i) = is-1
        ifson = -fils(i)
        in = ifson
 102    ino = in
        in =  frere(in)
        if (in.gt.0) go to 102
        nv(ino) = 0
        in = i
  888   if (subord(in).eq.0) go to 889
        in = subord(in)
        node(in) = is-1
        go to 888
  889   subord(in) = ino
        in = ino
        if (subord(in).eq.0) go to 887
        in = subord(in)
        ipe(in) = -i
  887   continue
      inos = -fils(ino)
      if (ifson.eq.ino) go to 107
      in = ifson
 105  ins = in
      in =  frere(in)
      if (in.ne.ino) go to 105
        if (inos.eq.0) then
          frere(ins) = -i
          go to 120
        else
          frere(ins) =  inos
        endif
 107    in = inos
        if (in.eq.0) go to 120
 108    int = in
        in =  frere(in)
        if (in.gt.0) go to 108
        frere(int) = -i
        go to 120
  110   is = is + 1
  120   ib = ipe(i)
        if (ib.ge.0) then
          if (ib.gt.0) na(il) = 0
          i = ib
          go to 160
        else
          i = -ib
          il = il + 1
        endif
  160 continue
  161 nsteps = is - 1
      return
      end
      subroutine ma57md(n,ne,irn,jcn,map,irnprm, &
                        lrow,perm,count,idiag)
      integer n,ne
      integer irn(ne),jcn(ne),map(ne),irnprm(n+ne),lrow(n),perm(n), &
              count(n), &
              idiag(n)
      integer expne,i,j,k
      do 10 i = 1,n
        count(i) = 1
        idiag(i) = 0
   10 continue
      expne = ne + n
      do 20 k = 1,ne
        i = irn(k)
        j = jcn(k)
        if (max(i,j).gt.n .or. min(i,j).lt.1) then
          expne = expne - 1
          go to 20
        endif
        if (i.eq.j) then
          i = perm(i)
          if (idiag(i).ge.1) then
            count(i) = count(i) + 1
            idiag(i) = idiag(i) + 1
          else
            idiag(i) = 1
            expne = expne - 1
          endif
          go to 20
        endif
        if (perm(i).lt.perm(j)) then
          i = perm(i)
          count(i) = count(i) + 1
        else
          j = perm(j)
          count(j) = count(j) + 1
        end if
   20 continue
      lrow(1) = count(1)
      idiag(1) = max(idiag(1),1)
      do 30 i = 2,n
        lrow(i) = count(i)
        count(i) = count(i-1) + lrow(i)
        idiag(i) = count(i-1) + max(idiag(i),1)
   30 continue
      do 35 i = 1,n
        k = perm(i)
        irnprm(idiag(k)) = i
   35 continue
      do 40 k = 1,ne
        i = irn(k)
        j = jcn(k)
        if (min(i,j).lt.1 .or. max(i,j).gt.n) then
          map(k) = 0
          go to 40
        endif
        i = perm(irn(k))
        j = perm(jcn(k))
        if (i.eq.j) then
          map(k) = idiag(i)
          irnprm(idiag(i)) = irn(k)
          idiag(i) = idiag(i) - 1
        else
          if (i.gt.j) then
            map(k) = count(j)
            irnprm(count(j)) = irn(k)
            count(j) = count(j) - 1
          else
            map(k) = count(i)
            irnprm(count(i)) = jcn(k)
            count(i) = count(i) - 1
          endif
        endif
   40 continue
      idiag(1) = expne
      return
      end
      subroutine ma57nd(n,lenr,na,ne,nd,nsteps,lstki,lstkr, &
                        info,rinfo)
      integer n,nsteps
      integer lenr(n),lstki(n),lstkr(n),na(nsteps), &
              nd(nsteps),ne(nsteps),info(40)
      double precision rinfo(20)
      integer i,iorg,istki,istkr,itop,itree,jorg,k, &
              lstk,nassr,nelim,nfr,nstk,ntotpv,nz1,nz2
      double precision delim
      intrinsic max
      double precision ops,opsass
      integer niradu,nirnec,nirtot,nrladu,nrlnec,nrltot,maxfrt
      nz1 = 0
      do 40 i = 1,n
        nz1 = nz1 + lenr(i)
   40 continue
      nz2 = nz1
      istki = 0
      istkr = 0
      ops = 0.0d0
      opsass = 0.0d0
      nrladu = 0
      niradu = 3
      nirtot = nz1+n+5
      nrltot = nz1
      nirnec = nz2+n+5
      nrlnec = nz2
      ntotpv = 0
      itop = 0
      maxfrt = 0
      do 100 itree = 1,nsteps
        nelim = ne(itree)
        delim = nelim
        nfr = nd(itree)
        maxfrt = max(maxfrt,nfr)
        nstk = na(itree)
        nassr = nelim*(nelim+1)/2 + nfr*nfr
        nrltot = max(nrltot,nrladu+nassr+istkr+nz1)
        nrlnec = max(nrlnec,nrladu+nassr+istkr+nz2)
        do 70 iorg = 1,nelim
          jorg = ntotpv + iorg
          opsass = opsass + lenr(jorg)
          nz2 = nz2 - lenr(jorg)
   70   continue
        ntotpv = ntotpv + nelim
        do 80 k = 1,nstk
          lstk = lstkr(itop)
          istkr = istkr - lstk
          opsass = opsass + lstk
          lstk = lstki(itop)
          istki = istki - lstk
          itop = itop - 1
   80   continue
        nrladu = nrladu + (nelim*(nelim+1))/2 + (nfr-nelim)*nelim
        niradu = niradu + 2 + nfr
        ops = ops + (delim* (12*nfr+6*nfr*nfr - (delim+1)* &
              (6*nfr+6-(2*delim+1))))/6 + delim
        if (nfr.gt.nelim) then
          itop = itop + 1
          lstkr(itop) = ((nfr-nelim)*(nfr-nelim+1))/2
          lstki(itop) = nfr - nelim + 1
          istki = istki + lstki(itop)
          istkr = istkr + lstkr(itop)
        endif
        if (itree.eq.nsteps) then
          nirtot = max(nirtot,niradu+istki+nz1)
          nirnec = max(nirnec,niradu+istki+nz2)
        else
          nirtot = max(nirtot,niradu+(n-ntotpv+2)+istki+nz1)
          nirnec = max(nirnec,niradu+(n-ntotpv+2)+istki+nz2)
        endif
  100 continue
      info(5)   = nrladu
      info(6)   = niradu
      info(7)   = maxfrt
      info(8)   = nsteps
      info(9)   = nrltot
      info(10)  = nirtot
      info(11)  = nrlnec
      info(12)  = nirnec
      rinfo(1)  = opsass
      rinfo(2)  = ops
      return
      end
      subroutine ma57od(n,ne,a,la,iw,liw,lrow,perm,nsteps,nstk,node, &
                        diag,schnab,ppos,cntl,icntl,info,rinfo,hold, &
                        biga)
      integer n,ne,la
      double precision a(la),diag(n),schnab(*),cntl(5),rinfo(20),biga
      integer liw,iw(liw),lrow(n),perm(n),nsteps,nstk(nsteps), &
              node(n),ppos(n),icntl(20),info(40),hold(40)
      integer zcol,rpos
      double precision zero,half,one
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0)
      integer ainput
      double precision amax,amult1,amult2
      integer apos,aposa,aposb,aposbb,aposbk,aposc,aposi,aposj,aposm, &
              apos1,apos2,apos3,apos4,astk,atrash,blk
      double precision delta,detpiv
      integer elt
      double precision flopsa,flopsb,flopsx
      integer i,i1,iass,ibeg,iell,iend,iexch,iinput,intspa, &
              iorg,ipiv,ipos,irow,isnpiv,istk,iswop,iwnfs,iwpos, &
              j,jay,ja1,jcol,jj,jjj,jmax,j1,j2,k, &
              kb,kblk,kct,kr,krow,k1,k2,l,laspiv,ldiag,liell, &
              lp,lpiv, nbstatic
      logical lastbk,ltwo
      integer maxfrt
      double precision maxpiv
      integer nass,nblk,nbloc,ncmpbi,ncmpbr,neig,nell,nfront,nirbdu
      double precision normj
      integer ntwo
      logical schur,lstat
      integer mpiv,npiv,npotpv,nrlbdu,nsc1,nst, &
              nstack(2),nstkac(2),ntotpv, &
              numorg,offdag,phase,pivblk
      double precision pivot
      integer pivsiz,poselt,pospv1,pospv2,ptra,ptrirn,rlspa, &
              sizblk,sizc,sizf,trlspa,tinspa,totsta(2),wp
      double precision rmax,swop,tmax,tol,uu,uloc,utarg,stctol
      double precision fd15ad
      character*80 outbuf
!?? To identify bug
      intrinsic min,max,abs
      external dgemm,fd15ad,ma57pd,ma57wd
      nbloc = icntl(11)
      tol = cntl(2)
      lp = icntl(1)
      wp = icntl(2)
      ldiag = icntl(5)
      uu = min(cntl(1),half)
      uu = max(uu,zero)
      lstat = .false.
      if (cntl(4).gt.zero) then
        if (cntl(5).eq.zero) lstat = .true.
        utarg = sqrt(uu/cntl(4))*cntl(4)
        stctol = biga*cntl(4)
      endif
      if (hold(1).gt.0) then
        info(1) = 0
        nblk = hold(2)
        ntwo = hold(3)
        info(23) = hold(4)
        ncmpbr = 0
        ncmpbi = 0
        neig   = hold(6)
        maxfrt = hold(7)
        iwpos  = hold(8)
        apos   = hold(9)
        aposbb = hold(10)
        nstkac(1) = hold(11)
        nstkac(2) = hold(12)
        ainput  = hold(13)
        iinput  = hold(14)
        istk    = hold(15)
        astk    = hold(16)
        intspa  = hold(17)
        rlspa   = hold(18)
        ptrirn  = hold(19)
        ptra    = hold(20)
        ntotpv  = hold(21)
        npotpv  = hold(22)
        numorg  = hold(23)
        nfront  = hold(24)
        nass    = hold(25)
        if (hold(1).eq.1) nell    = hold(27)
        if (hold(1).eq.2) npiv    = hold(27)
        iass    = hold(28)
        tinspa  = hold(29)
        trlspa  = hold(30)
        totsta(1) = hold(31)
        totsta(2) = hold(32)
        nstack(1) = hold(33)
        nstack(2) = hold(34)
        info(32)  = hold(37)
        info(33)  = hold(38)
        info(34)  = hold(39)
        nbstatic  = hold(40)
        if (icntl(7).eq.2 .or.icntl(7).eq.3) isnpiv = hold(35)
        if (icntl(7).eq.4) phase = hold(36)
        if (hold(1).eq.2) nsc1    = nfront-npiv
        flopsa = rinfo(3)
        flopsb = rinfo(4)
        flopsx = rinfo(5)
        if (hold(1).eq.1) then
          hold(1) = 0
          go to 333
        else
          hold(1) = 0
          go to 444
        endif
      endif
      nbstatic = 0
      nblk = 0
      ntwo = 0
      ncmpbr = 0
      ncmpbi = 0
      flopsa = zero
      flopsb = zero
      flopsx = zero
      neig = 0
      maxfrt  = 0
      info(1) = 0
      info(2) = 0
      info(17) = 0
      info(40) = 0
      info(26) = 0
      info(27) = 0
      info(32) = 0
      info(33) = 0
      info(34) = 0
      info(23) = 0
      rinfo(3) = zero
      rinfo(4) = zero
      rinfo(5) = zero
      rinfo(15) = zero
      do 10 i = 1,n
        ppos(i) = n + 1
   10 continue
      iwpos = 6
      iw(1) = 0
      iw(2) = 0
      iw(3) = 0
      iw(4) = 0
      iw(5) = 0
      aposbb = 1
      nstack(1) = 0
      nstack(2) = 0
      nstkac(1) = ne
      nstkac(2) = ne
      totsta(1) = ne
      totsta(2) = ne
      intspa = ne+5+n
      rlspa = ne
      tinspa = ne+5+n
      trlspa = ne
      ptrirn = liw - ne + 1
      ptra = la - ne + 1
      istk = ptrirn - 1
      astk = ptra - 1
      ainput = ptra
      iinput = ptrirn
      ntotpv = 0
      npotpv = 0
      if (icntl(7).eq.2 .or. icntl(7).eq.3) isnpiv = 0
      if (icntl(7).eq.4) then
        phase = 1
        do 19 i = 1,n
          diag(i) = zero
   19   continue
        apos1 = ptra-1
        j1 = ptrirn
        do 20 i = 1,n
          j2 = j1 + lrow(i) - 1
          do 25 jj = j1,j2
            j = iw(jj)
            apos1 = apos1 + 1
            if (j.eq.perm(i)) diag(j) = diag(j) + a(apos1)
   25     continue
          j1 = j2 + 1
   20   continue
        schnab(1) = one
        schnab(5) = zero
        do 21 i = 1,n
          schnab(1) = max(schnab(1),abs(diag(i)))
          schnab(5) = min(schnab(5),diag(i))
   21   continue
        schnab(4) = schnab(1)
        schnab(2) = fd15ad('E')**(1.0/3.0)
        schnab(3) = 0.1
        rinfo(15) = fd15ad('H')
        delta     = zero
      endif
      iass = 1
 2160 continue
        numorg = 0
        do 30 i = npotpv + 1,n
          j = perm(i)
          if (abs(node(j)).gt.iass) go to 40
          iw(iwpos+numorg) = j
          numorg = numorg + 1
          ppos(j) = numorg
   30   continue
   40   nass = numorg
        nell = nstk(iass)
        iell = istk + 1
        do 70 elt = 1,nell
          do 50 jj = iell + 1,iell + iw(iell)
            j = iw(jj)
            if (node(j).gt.iass) go to 50
            if (ppos(j).le.n) go to 50
            iw(iwpos+nass) = j
            nass = nass + 1
            ppos(j) = nass
   50     continue
          iell = iell + iw(iell) + 1
   70   continue
        iwnfs = iwpos + nass
        j1 = ptrirn
        do 90 iorg = 1,numorg
          j2 = j1 + lrow(npotpv+iorg) - 1
          do 80 jj = j1,j2
            j = iw(jj)
            if (ppos(j).le.n) go to 80
            iw(iwnfs) = j
            iwnfs = iwnfs + 1
            ppos(j) = iwnfs - iwpos
   80     continue
          j1 = j2 + 1
   90   continue
        iell = istk + 1
        do 170 elt = 1,nell
          j1 = iell+1
          j2 = iell+iw(iell)
          do 150 jj = j1,j2
            j = iw(jj)
            if (ppos(j).le.n) go to 150
            iw(iwnfs) = j
            iwnfs = iwnfs + 1
            ppos(j) = iwnfs - iwpos
  150     continue
          iell = j2 + 1
  170   continue
        nfront = iwnfs - iwpos
        maxfrt = max(maxfrt,nfront)
        if (info(1).ne.-3) then
          apos = aposbb + (nass*(nass+1))/2
        else
          apos = 1
        end if
        rlspa  = max(rlspa,info(40)+apos+nfront*nfront-1+nstkac(1))
        trlspa = max(trlspa,info(40)+apos+nfront*nfront-1+totsta(1))
  333   if (apos+nfront*nfront-1.gt.astk) then
          call ma57pd(a,iw,astk,ainput,ptra,.true.)
          ncmpbr = ncmpbr + 1
          if (apos+nfront*nfront-1.gt.astk) then
            if (icntl(8).ne.0) then
              hold(1) = 1
              hold(2) = nblk
              hold(3) = ntwo
              hold(4) = info(23)
              hold(5) = ncmpbi
              hold(6) = neig
              hold(7) = maxfrt
              hold(8) = iwpos
              hold(9) = apos
              hold(10) = aposbb
              hold(11) = nstkac(1)
              hold(12) = nstkac(2)
              hold(13) = ainput
              hold(14) = iinput
              hold(15) = istk
              hold(16) = astk
              hold(17) = intspa
              hold(18) = rlspa
              hold(19) = ptrirn
              hold(20) = ptra
              hold(21) = ntotpv
              hold(22) = npotpv
              hold(23) = numorg
              hold(24) = nfront
              hold(25) = nass
              hold(27) = nell
              hold(28) = iass
              hold(29) = tinspa
              hold(30) = trlspa
              hold(31) = totsta(1)
              hold(32) = totsta(2)
              hold(33) = nstack(1)
              hold(34) = nstack(2)
              if (icntl(7).eq.2 .or.icntl(7).eq.3) hold(35) = isnpiv
              if (icntl(7).eq.4) hold(36) = phase
              hold(37) = info(32)
              hold(38) = info(33)
              hold(39) = info(34)
              rinfo(3) =flopsa
              rinfo(4) =flopsb
              rinfo(5) =flopsx
              hold(40) = nbstatic
              info(35) = hold(40)
              info(1) = 10
              return
            else
              info(40) = info(40) + apos - 1
              apos = 1
              aposbb = 1
              info(1) = -3
              if (nfront*nfront.gt.astk) then
                info(17) = max(info(17),rlspa)
                if (icntl(7).eq.4) info(17) = max(info(17),rlspa + n)
                info(2) = la
                return
              endif
            endif
          endif
        end if
        atrash = apos + nfront*nfront - 1
        do 210 jj = apos,atrash
          a(jj) = zero
  210   continue
        j1 = ptrirn
        do 230 iorg = 1,numorg
          j = perm(npotpv+iorg)
          aposi = apos + (ppos(j)-1)*nfront - 1
          j2 = j1 + lrow(npotpv+iorg) - 1
          flopsa = flopsa + j2 - j1 + 1
          do 220 jj = j1,j2
            jay = iw(jj)
            if (ppos(jay).ge.ppos(j)) then
              apos2 = aposi + ppos(jay)
            else
              apos2 = apos + (ppos(jay)-1)*nfront + ppos(j) - 1
            endif
            a(apos2) = a(apos2) + a(ptra)
            ptra = ptra + 1
  220     continue
          nstkac(1) = nstkac(1) - j2 + j1 - 1
          j1 = j2 + 1
  230   continue
        nstkac(2) = nstkac(2) - j1 + ptrirn
        ptrirn = j1
        npotpv = npotpv + numorg
!???
!?? Depends if we need lower triangle and whether all entries are already
        do 380 elt = 1,nell
          poselt = astk + 1
          liell = iw(istk+1)
          j1 = istk + 2
          j2 = istk+1 + liell
          flopsa = flopsa + (liell*(liell+1))/2
          do 250 jj = j1,j2
            j = iw(jj)
            apos2 = apos + (ppos(j)-1)*nfront
            apos1 = poselt
            do 240 jjj=jj,j2
              jay = iw(jjj)
!???          APOS3 = APOS2 + PPOS(JAY) - 1
!???          A(APOS3) = A(APOS3) + A(APOS1)
!???          APOS5 = APOS+(PPOS(JAY)-1)*NFRONT+PPOS(J)-1
!???          IF (APOS3.NE.APOS5) A(APOS5) = A(APOS5) + A(APOS1)
              if (ppos(jay) .ge. ppos(j)) then
                apos3 = apos2 + ppos(jay) - 1
              else
                apos3 = apos+(ppos(jay)-1)*nfront+ppos(j)-1
              endif
              a(apos3) = a(apos3) + a(apos1)
              apos1 = apos1 + 1
  240       continue
            poselt = poselt + liell - (jj-j1)
  250     continue
          nstkac(2) = nstkac(2) - (j2-istk)
          nstack(2) = nstack(2) - (j2-istk)
          totsta(2) = totsta(2) - (j2-istk)
          istk = j2
          astk = astk + (liell*(liell+1))/2
          nstkac(1) = nstkac(1) - (liell*(liell+1))/2
          nstack(1) = nstack(1) - (liell*(liell+1))/2
          totsta(1) = totsta(1) - (liell*(liell+1))/2
  380   continue
!1122     CONTINUE
        pivblk = min(nbloc,nass)
        aposbk = apos
        npiv = 0
        uloc = uu
        do 918 blk = 1,nass
        if (npiv+pivblk .ge. nass) then
          lastbk = .true.
          sizblk = nass - npiv
        else
          lastbk = .false.
          sizblk = pivblk
        endif
        laspiv = npiv
        mpiv = 0
!CCCCCCC
        kr = 0
!CC Set to following to force 2 by 2 pivots in Nocedal examples
        kct = sizblk + 1
  920   continue
          kr = kr + 1
          kct = kct - 1
          if (kct.eq.0) go to 930
          if (kr.gt.sizblk) kr = mpiv + 1
          ipiv = laspiv + kr
            aposi = apos + (ipiv-1)*nfront
            pospv1 = aposi + ipiv - 1
            pivot = a(pospv1)
   29       if (icntl(7).eq.4 .and. phase.eq.2) then
              if (info(27).eq.0) info(27) = ntotpv + 1
              normj = zero
              do 28 i = pospv1+1,pospv1+nfront-npiv-1
                normj = normj + abs(a(i))
   28         continue
              delta = max(zero, &
                          - a(pospv1) + max(normj,schnab(2)*schnab(1)))
              a(pospv1) = a(pospv1) + delta
              if (a(pospv1).eq.zero) go to 970
              rinfo(15) = min(rinfo(15),a(pospv1))
              diag(perm(ntotpv+1)) = delta
              pivsiz = 1
              go to 811
            endif
            if (icntl(7).gt.1) then
              if (abs(pivot).le.cntl(2)) then
                if (icntl(7).lt.4) go to 970
                phase = 2
                go to 29
              endif
              if (ntotpv.eq.0) then
                if (pivot.gt.zero) isnpiv = 1
                if (pivot.lt.zero) isnpiv = -1
              else
                if (icntl(7).eq.2 .and. isnpiv*pivot.lt.zero) go to 980
                if (icntl(7).eq.3 .and. isnpiv*pivot.lt.zero) then
                    info(26) = info(26) + 1
                    isnpiv = -isnpiv
                endif
              endif
              if (icntl(7).eq.4) then
                if (pivot.ge.schnab(1)*schnab(2) .and. &
                    schnab(5).ge.-schnab(3)*schnab(4)) then
                  schnab(5) = zero
                  schnab(4) = zero
                  do 22 i = pospv1+1,pospv1+nfront-npiv-1
                    j = iw(iwpos+npiv+i-pospv1)
                    diag(j) = diag(j) - a(i)*a(i)/pivot
                    schnab(5) = min(diag(j),schnab(5))
                    schnab(4) = max(diag(j),schnab(4))
                    if (diag(j).lt.-schnab(3)*schnab(1)) then
                      phase = 2
                      go to 29
                    endif
   22             continue
                  diag(perm(ntotpv+1)) = zero
                  rinfo(15) = min(rinfo(15),pivot)
                else
                  phase = 2
                  go to 29
                endif
              endif
              pivsiz = 1
              go to 811
            endif
            amax = zero
            jmax = 0
            do 110 k = 1, ipiv - npiv - 1
              if (abs(a(pospv1-k*nfront)).gt.amax) then
                amax = abs(a(pospv1-k*nfront))
                jmax = ipiv - k
              endif
  110       continue
            do 111 k =  1, min(nass,laspiv+pivblk) - ipiv
              if (abs(a(pospv1+k)).gt.amax) then
                amax = abs(a(pospv1+k))
                jmax = ipiv + k
              endif
  111       continue
            rmax = zero
            do 112 k = min(nass,laspiv+pivblk)-ipiv+1,nfront-ipiv
               rmax = max(rmax,abs(a(pospv1+k)))
 112        continue
            if (max(amax,rmax,abs(pivot)).le.tol) then
              go to 920
            end if
            if (max(amax,abs(pivot)).le.tol) go to 920
            pivsiz = 0
            if (abs(pivot).gt.uloc*max(rmax,amax)) then
              pivsiz = 1
              a(pospv1) = pivot
              go to 810
            end if
            if (npiv+1.eq.nass) then
              a(pospv1) = pivot
              go to 920
            end if
            if (amax.le.tol) go to 920
            if (rmax.lt.amax) then
              rmax = zero
              do 113 k = 1, ipiv - npiv - 1
                if (ipiv-k.eq.jmax) go to 113
                rmax=max(rmax,abs(a(pospv1-k*nfront)))
  113         continue
              do 114 k =  1, nfront - ipiv
                if (ipiv+k.eq.jmax) go to 114
                rmax = max(rmax,abs(a(pospv1+k)))
  114         continue
            endif
            aposj = apos + (jmax-1)*nfront
            pospv2 = aposj + jmax - 1
            if (ipiv.gt.jmax) then
              offdag = aposj + ipiv - 1
            else
              offdag = aposi + jmax - 1
            end if
            tmax = zero
            do 115 k = 1, jmax - npiv - 1
              if (jmax-k.eq.ipiv) go to 115
              tmax=max(tmax,abs(a(pospv2-k*nfront)))
  115       continue
            do 116 k =  1, nfront - jmax
              if (jmax+k.eq.ipiv) go to 116
              tmax = max(tmax,abs(a(pospv2+k)))
  116       continue
            detpiv = a(pospv1)*a(pospv2) - amax*amax
            maxpiv = max(abs(a(pospv1)),abs(a(pospv2)))
            if (maxpiv.eq.zero) maxpiv = one
            if (abs(detpiv)/maxpiv.le.tol) go to 920
            pivsiz = 2
            if ((abs(a(pospv2))*rmax+amax*tmax)*uloc.gt. &
                abs(detpiv)) go to 920
            if ((abs(a(pospv1))*tmax+amax*rmax)*uloc.gt. &
                abs(detpiv)) go to 920
  810       lpiv = ipiv
            if (pivsiz.eq.2) lpiv = min(ipiv,jmax)
!CC         KR = MAX(KR,NPIV+PIVSIZ)
            kr = max(kr,mpiv+pivsiz)
            kct = sizblk - mpiv - pivsiz + 1
            do 860 krow = npiv,npiv + pivsiz - 1
              if (lpiv.eq.krow+1) go to 850
              ja1 = apos + (lpiv-1)
              j1 = apos + krow
              do 820 jj = 1,krow
                swop = a(ja1)
                a(ja1) = a(j1)
                a(j1) = swop
                ja1 = ja1 + nfront
                j1 = j1 + nfront
  820         continue
              ja1 = ja1 + nfront
              j1 = j1 + 1
              do 830 jj = 1,lpiv - krow - 2
                swop = a(ja1)
                a(ja1) = a(j1)
                a(j1) = swop
                ja1 = ja1 + nfront
                j1 = j1 + 1
  830         continue
              swop = a(apos+krow* (nfront+1))
              a(apos+krow* (nfront+1)) = a(ja1)
              a(ja1) = swop
              do 840 jj = 1,nfront - lpiv
                ja1 = ja1 + 1
                j1 = j1 + 1
                swop = a(ja1)
                a(ja1) = a(j1)
                a(j1) = swop
  840         continue
              ipos = iwpos + krow
              iexch = iwpos + lpiv - 1
              iswop = iw(ipos)
              iw(ipos) = iw(iexch)
              iw(iexch) = iswop
  850         lpiv = max(ipiv,jmax)
  860       continue
  811       pospv1 = apos + npiv* (nfront+1)
            pospv2 = pospv1 + nfront + 1
            if (pivsiz.eq.1) then
              flopsb = flopsb + one
              a(pospv1) = one/a(pospv1)
              if (a(pospv1).lt.zero) neig = neig + 1
              j1 = pospv1 + 1
              j2 = pospv1 + nass - (npiv+1)
              ibeg = pospv1 + nfront + 1
              iend = apos + (npiv+1)*nfront + nfront - 1
              do 880 jj = j1,j2
                amult1 = -a(jj)*a(pospv1)
                if (.not.lastbk) a(pospv1+(jj-j1+1)*nfront) = a(jj)
                jcol = jj
                flopsb = flopsb + (iend-ibeg+1)*2 + 1
                if (mpiv+jj-j1+2.gt.pivblk) go to 871
!DIR$            IVDEP
                do 870 irow = ibeg,iend
                  a(irow) = a(irow) + amult1*a(jcol)
                  jcol = jcol + 1
  870           continue
  871           a(jj) = amult1
                ibeg = ibeg + nfront + 1
                iend = iend + nfront
  880         continue
              npiv = npiv + 1
              mpiv = mpiv + 1
              ntotpv = ntotpv + 1
              if (mpiv.eq.sizblk) go to 930
            else
              offdag = pospv1 + 1
              flopsb = flopsb + 6.0
              swop = a(pospv2)
              if (detpiv.lt.zero) then
                neig = neig + 1
              else
                if (swop.lt.zero) neig = neig + 2
              end if
              a(pospv2) = a(pospv1)/detpiv
              a(pospv1) = swop/detpiv
              a(offdag) = -a(offdag)/detpiv
              j1 = pospv1 + 2
              j2 = pospv1 + nass - (npiv+1)
              ibeg = pospv2 + nfront + 1
              iend = apos + (npiv+2)*nfront + nfront - 1
              do 900 jj = j1,j2
                k1 = jj
                k2 = jj + nfront
                amult1 = - (a(pospv1)*a(k1)+a(pospv1+1)*a(k2))
                amult2 = - (a(pospv1+1)*a(k1)+a(pospv2)*a(k2))
                if (.not.lastbk) then
                  a(pospv1 + (jj-j1+2)*nfront) = a(k1)
                  a(pospv1 + (jj-j1+2)*nfront + 1) = a(k2)
                endif
                flopsb = flopsb + (iend-ibeg+1)*4 + 6
                if (mpiv+jj-j1+3.gt.pivblk) go to 891
!DIR$            IVDEP
                do 890 irow = ibeg,iend
                  a(irow) = a(irow) + amult1*a(k1) + amult2*a(k2)
                  k1 = k1 + 1
                  k2 = k2 + 1
  890           continue
  891           a(jj) = amult1
                a(jj+nfront) = amult2
                ibeg = ibeg + nfront + 1
                iend = iend + nfront
  900         continue
              ipos = iwpos + npiv
              iw(ipos) = -iw(ipos)
              iw(ipos+1) = -iw(ipos+1)
              npiv = npiv + 2
              mpiv = mpiv + 2
              ntotpv = ntotpv + 2
              ntwo = ntwo + 1
              if (mpiv.eq.sizblk) go to 930
            end if
        go to 920
 930    if (lastbk) then
          if (npiv.eq.nass) go to 935
          if (.not. lstat)  go to 935
          uloc = uloc/10.0d0
          if (uloc.lt.utarg) then
            uloc = uloc * 10.0d0
            go to 9919
          endif
          kct = sizblk + 1 - mpiv
          go to 920
        endif
        if (mpiv.eq.0) then
          pivblk = 2*pivblk
          go to 918
        endif
        kblk = (nass-(laspiv+pivblk))/pivblk
        l = nass - (laspiv+pivblk)
        apos4 = apos+(laspiv+pivblk)*(nfront+1)
        do 931 kb = 1,kblk
          flopsx = flopsx + pivblk*(pivblk-1)*mpiv
          call dgemm('N','N',l-(kb-1)*pivblk,pivblk,mpiv,one, &
                     a(aposbk+pivblk*kb),nfront, &
                     a(aposbk+pivblk*kb*nfront),nfront,one, &
                     a(apos4+pivblk*(kb-1)*(nfront+1)),nfront)
          if (nfront.gt.nass) &
          call dgemm('N','T',nfront-nass,pivblk,mpiv,one, &
                     a(aposbk+nass-laspiv),nfront, &
                     a(aposbk+pivblk*kb),nfront,one, &
                     a(aposbk+kb*nfront*pivblk+nass-laspiv),nfront)
  931   continue
       sizc = nass - (kblk+1)*pivblk - laspiv
       sizf = nfront - (kblk+1)*pivblk - laspiv
       aposa = aposbk + (kblk+1)*pivblk
       do 934 k = 1,mpiv
         aposb = aposbk + nfront*pivblk*(kblk+1) + k - 1
         aposm = aposbk + pivblk*(kblk+1) + (k-1)*nfront
         aposc = aposbk + pivblk*(kblk+1)*(nfront+1)
         do 933 jj = 1,sizc
            do 932 j = jj,sizc
              a(aposc+j-1) = a(aposc+j-1) + a(aposa+j-1)*a(aposb)
  932       continue
            do 936 j = sizc+1,sizf
              a(aposc+j-1) = a(aposc+j-1) + a(aposa+j-1)*a(aposm)
  936       continue
            aposc = aposc + nfront
            aposb = aposb + nfront
            aposm = aposm + 1
  933     continue
          aposa = aposa + nfront
  934   continue
        aposbk = aposbk + mpiv*(nfront+1)
        laspiv = npiv
  918   continue
        if (lp.ge.0) then
           call wrtlin( '    *** HSL ERROR ***' )
           call wrtlin( '****** BE WORRIED LOOP 918' )
        end if
 9919      ipiv = laspiv+mpiv
 9920      ipiv = ipiv + 1
!ADD Probably not needed .. use only IPIV
           aposi = apos + (ipiv-1)*nfront
           pospv1 = aposi + ipiv - 1
           pivot = a(pospv1)
!ADD
           pivsiz = 1
           lpiv = ipiv
           amax = zero
           do 9876 k = 1, ipiv - npiv - 1
             amax = max(amax,abs(a(pospv1-k*nfront)))
 9876      continue
           do 9878 k =  1, nfront - ipiv
             amax = max(amax,abs(a(pospv1+k)))
 9878      continue
           if (abs(a(pospv1)).lt.stctol) then
               pivot = stctol
              if (a(pospv1) .lt. zero) then
                 a(pospv1) = -pivot
                 pivot     = -pivot
              else
                 a(pospv1) = pivot
              endif
              nbstatic = nbstatic + 1
           endif
           flopsb = flopsb + one
           a(pospv1) = one/a(pospv1)
           if (a(pospv1).lt.zero) neig = neig + 1
           j1 = pospv1 + 1
           j2 = pospv1 + nass - (npiv+1)
           ibeg = pospv1 + nfront + 1
!BUG
           iend = aposi + 2*nfront - 1
           do 9880 jj = j1,j2
              amult1 = -a(jj)*a(pospv1)
!ADD  Not needed since LASTBK always true
              jcol = jj
              flopsb = flopsb + (iend-ibeg+1)*2 + 1
!ADD Not necessary because control is on setting of J2
!DIR$            IVDEP
              do 9870 irow = ibeg,iend
                 a(irow) = a(irow) + amult1*a(jcol)
                 jcol = jcol + 1
 9870         continue
 9871         a(jj) = amult1
              ibeg = ibeg + nfront + 1
              iend = iend + nfront
 9880      continue
           npiv = npiv + 1
           mpiv = mpiv + 1
           ntotpv = ntotpv + 1
           if (mpiv.lt.sizblk) go to 9920
!ADD
  935   schur = (nbloc.lt.(nfront-nass) .and. npiv.ge.nbloc)
        nsc1 = nfront - npiv
        if (iass.ne.nsteps) info(23) = info(23) + nass - npiv
        if (cntl(4).gt.zero .and. info(23).gt.cntl(5)*n) lstat = .true.
        if (nsc1.eq.0) go to 1830
        if (.not.schur) then
          rlspa = max(rlspa,info(40)+apos+nfront*nfront-1+ &
                            nstkac(1))
          trlspa = max(trlspa,info(40)+apos+nfront*nfront-1+ &
                            totsta(1))
          nstkac(1) = nstkac(1) + ((nsc1+1)*nsc1)/2
          nstack(1) = nstack(1) + ((nsc1+1)*nsc1)/2
          totsta(1) = totsta(1) + ((nsc1+1)*nsc1)/2
          aposi = apos + nfront*nfront - 1
          do 1370 jj = 1,nfront-npiv
            j = aposi
            do 1360 jjj = 1,jj
                a(astk) = a(j)
                astk = astk - 1
                j = j - 1
 1360       continue
            aposi = aposi - nfront
 1370     continue
          apos4 = astk + 1
          j1 = iwpos
          ltwo = .false.
          pospv1 = apos
          do 1450 i1 = 1,npiv
            if (ltwo) go to 1440
            aposi = apos + (i1-1)*nfront + nass
            j2 = apos + nfront* (i1-1) + nfront - 1
!CC What happens here ??
            aposc = apos4 + ((nass-npiv)*(2*nfront-npiv-nass+1))/2
            if (iw(j1).gt.0) then
              flopsb = flopsb + (nfront-nass) + &
                                (nfront-nass)* (nfront-nass+1)
              do 1410 jj = aposi,j2
                amult1 = -a(jj)*a(pospv1)
                do 1400 jjj = jj,j2
                  a(aposc) = a(aposc) + amult1*a(jjj)
                  aposc = aposc + 1
 1400           continue
                a(jj) = amult1
 1410         continue
              j1 = j1 + 1
            else
              pospv2 = pospv1 + nfront + 1
              offdag = pospv1 + 1
              flopsb = flopsb + 6* (nfront-nass) + &
                       2* (nfront-nass)* (nfront-nass+1)
              do 1430 jj = aposi,j2
                amult1 = - (a(pospv1)*a(jj)+a(offdag)*a(jj+nfront))
                amult2 = -a(pospv2)*a(jj+nfront) - a(offdag)*a(jj)
                do 1420 jjj = jj,j2
                  a(aposc) = a(aposc) + amult1*a(jjj) + &
                             amult2*a(jjj+nfront)
                  aposc = aposc + 1
 1420           continue
                a(jj) = amult1
                a(jj+nfront) = amult2
 1430         continue
              j1 = j1 + 2
              pospv1 = pospv2
              ltwo = .true.
              go to 1450
            end if
 1440       ltwo = .false.
            pospv1 = pospv1 + nfront + 1
 1450     continue
        else
          apos4 = apos+nass*(nfront+1)
        apos3 = apos+nass*nfront
        j1 = iwpos
        ltwo = .false.
        pospv1 = apos
          do 1490 i = 1,npiv
            if (ltwo) go to 1480
            aposi = apos + (i-1)*nfront + nass
            poselt = apos3 + i - 1
            if (iw(j1).gt.0) then
              flopsb = flopsb + (nfront-nass)
              do 1460 jj = aposi,apos + nfront*i - 1
                a(poselt) = a(jj)
                a(jj) = -a(jj)*a(pospv1)
                poselt = poselt + nfront
 1460         continue
              j1 = j1 + 1
            else
              pospv2 = pospv1 + nfront + 1
              offdag = pospv1 + 1
              flopsb = flopsb + 6* (nfront-nass)
              do 1470 jj = aposi,apos + nfront*i - 1
                a(poselt) = a(jj)
                a(poselt+1) = a(jj+nfront)
                a(jj) = - (a(pospv1)*a(jj)+a(offdag)*a(jj+nfront))
                a(jj+nfront) = -a(pospv2)*a(jj+nfront) - &
                               a(offdag)*a(poselt)
                poselt = poselt + nfront
 1470         continue
              j1 = j1 + 2
              pospv1 = pospv2
              ltwo = .true.
              go to 1490
            end if
 1480       ltwo = .false.
            pospv1 = pospv1 + nfront + 1
 1490     continue
          flopsb = flopsb + npiv* (nfront-nass)**2 + &
                            npiv* (nfront-nass)
          kblk = ( nfront-nass)/nbloc
          l =  nfront - nass
          do 1500 kb = 1,kblk
            flopsx = flopsx + nbloc* (nbloc-1)* (npiv)
            call dgemm('N','N',l-(kb-1)*nbloc,nbloc,npiv,one, &
                       a(apos+nass+nbloc*(kb-1)),nfront, &
                       a(apos3+nbloc*(kb-1)*nfront),nfront,one, &
                       a(apos4+nbloc*(nfront+1)*(kb-1)),nfront)
 1500     continue
          do 1550 i = 1 + kblk*nbloc,l
            aposa = apos + nass
            aposb = apos3 +(i-1)*nfront
            aposc = apos4 + (i-1)*nfront - 1
            do 1540 k = 1,npiv
              do 1530 j = i,l
                a(aposc+j) = a(aposc+j) + a(aposa+j-1)*a(aposb)
 1530         continue
              aposa = aposa + nfront
              aposb = aposb + 1
 1540       continue
 1550     continue
          ja1 = apos+nfront*nfront-1
          nstkac(1) = nstkac(1) + ((nsc1+1)* (nsc1))/2
          nstack(1) = nstack(1) + ((nsc1+1)* (nsc1))/2
          totsta(1) = totsta(1) + ((nsc1+1)* (nsc1))/2
          do 1710 i = nsc1,1,-1
            do 1700 jj = ja1,ja1-(nsc1-i),-1
              a(astk) = a(jj)
              astk = astk - 1
 1700       continue
            ja1 = ja1 - nfront
 1710     continue
        end if
        nstkac(2) = nstkac(2) + nsc1 + 1
        nstack(2) = nstack(2) + nsc1 + 1
        totsta(2) = totsta(2) + nsc1 + 1
 1830   if (iass.eq.nsteps) then
          intspa = max(intspa,iwpos+nfront-1+nstkac(2))
          tinspa = max(tinspa,iwpos+nfront-1+totsta(2))
          go to 2158
        else
          intspa = max(intspa,iwpos+nfront-1+(n-ntotpv+2)+nstkac(2))
          tinspa = max(tinspa,iwpos+nfront-1+(n-ntotpv+2)+totsta(2))
        endif
  444   nst = 0
        if (nsc1.gt.0) nst = nsc1 + 1
        if (iwpos+nfront-1+(n-ntotpv+2)+nst.gt.istk) then
          call ma57pd(a,iw,istk,iinput,ptrirn,.false.)
          ncmpbi = ncmpbi + 1
          if (iwpos+nfront-1+(n-ntotpv+2)+nst.gt.istk) then
            if (icntl(8).ne.0) then
              hold(1) = 2
              hold(2) = nblk
              hold(3) = ntwo
              hold(4) = info(23)
              hold(5) = ncmpbi
              hold(6) = neig
              hold(7) = maxfrt
              hold(8) = iwpos
              hold(9) = apos
              hold(10) = aposbb
              hold(11) = nstkac(1)
              hold(12) = nstkac(2)
              hold(13) = ainput
              hold(14) = iinput
              hold(15) = istk
              hold(16) = astk
              hold(17) = intspa
              hold(18) = rlspa
              hold(19) = ptrirn
              hold(20) = ptra
              hold(21) = ntotpv
              hold(22) = npotpv
              hold(23) = numorg
              hold(24) = nfront
              hold(25) = nass
              hold(27) = npiv
              hold(28) = iass
              hold(29) = tinspa
              hold(30) = trlspa
              hold(31) = totsta(1)
              hold(32) = totsta(2)
              hold(33) = nstack(1)
              hold(34) = nstack(2)
              if (icntl(7).eq.2 .or.icntl(7).eq.3) hold(35) = isnpiv
              if (icntl(7).eq.4) hold(36) = phase
              hold(37) = info(32)
              hold(38) = info(33)
              hold(39) = info(34)
              nsc1    = nfront-npiv
              rinfo(3) =flopsa
              rinfo(4) =flopsb
              rinfo(5) =flopsx
              info(1) = 11
              hold(40) = nbstatic
              info(35) = hold(40)
            else
              info(1)  = -4
              info(2)  = liw
              info(18) = intspa
            endif
            return
          end if
        end if
        if (nsc1.gt.0) then
          do 1720 i = 1,nsc1
            iw(istk) = iw(iwpos+nfront-i)
            istk = istk - 1
 1720     continue
          iw(istk) = nsc1
          istk = istk - 1
        endif
        do 1840 jj = iwpos + npiv,iwpos + nfront - 1
          j = abs(iw(jj))
          ppos(j) = n + 1
 1840   continue
!********************************
!********************************
 2158   if (npiv.eq.0) go to 2159
        nblk = nblk + 1
        iw(iwpos-2) = nfront
        iw(iwpos-1) = npiv
        iwpos = iwpos + nfront + 2
        if (info(1).eq.-3) then
          info(40) = info(40) + (npiv * (2*nfront-npiv+1))/2
          go to 2159
        end if
        apos2 = aposbb
        do 2130 i = 1,npiv
          ja1 = apos + (i-1)* (nfront+1)
          do 2120 j = i,npiv
            a(apos2) = a(ja1)
            if (a(apos2).eq.zero) info(32) = info(32) + 1
            apos2 = apos2 + 1
            ja1 = ja1 + 1
 2120     continue
 2130   continue
        rpos = apos2
        do 2150 i = 1,npiv
          ja1 = apos + (i-1)*nfront + npiv
          do 2140 j = 1,nfront - npiv
            a(apos2) = a(ja1)
            apos2 = apos2 + 1
            ja1 = ja1 + 1
 2140     continue
 2150   continue
        aposbb = apos2
        do 2152 j = 1,nfront-npiv
        apos2 = rpos+j-1
        zcol = 1
          do 2151 i = 1,npiv
            if (a(apos2).eq.zero) info(33) = info(33)+1
            if (a(apos2).ne.zero) zcol = 0
            apos2 = apos2 + nfront - npiv
 2151     continue
        if (zcol.eq.1) info(34) = info(34)+1
 2152   continue
 2159   iass = iass + 1
      if (iass.le.nsteps) then
        iw(iwpos-2) = 0
        iw(iwpos-1) = 0
        go to 2160
      endif
!2160 CONTINUE
      info(35) = nbstatic
      if (info(1).eq.-3) then
        info(2)  = la
        info(17) = max(info(17),rlspa)
        if (icntl(7).eq.4) info(17) = max(info(17),rlspa + n)
        return
      end if
      go to 1000
 970  info(1) = -5
      info(2) = ntotpv + 1
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write( outbuf, '(2a,i3)' ) '*** Error message from routine MA57BD **',&
                                    '   INFO(1) = ', info(1)
         call wrtlin( outbuf )
         write( outbuf, '(a,d16.8,a,d16.8)' ) 'Pivot has value ', pivot, &
          ' when CNTL(2) has value ', cntl(2)
         call wrtlin( outbuf )
         write( outbuf, '(a,i11,a,i3)' ) 'at stage', info(2), &
          '  when ICNTL(7) =', icntl(7)
         call wrtlin( outbuf )
      end if
      return
 980  info(1) = -6
      info(2) = ntotpv + 1
      if (ldiag.gt.0 .and. lp.ge.0) then
         call wrtlin( '    *** HSL ERROR ***' )
         write( outbuf, '(2a,i3)' ) '*** Error message from routine MA57BD **',&
          '   INFO(1) = ', info(1)
         call wrtlin( outbuf )
         write( outbuf, '(a,i10,a,i3)' )  'Change in sign of pivot at stage', &
          info(2), '  when ICNTL(7) = ', icntl(7)
         call wrtlin( outbuf )
      end if
      return
 1000 nrlbdu = aposbb - 1
      nirbdu = iwpos - 3
      if (ntotpv.ne.n) then
        info(1) = 4
        if (ldiag.gt.0 .and. wp.ge.0) then
           call wrtlin( '    *** HSL WARNING ***' )
           write( outbuf, '(2a,i2)' ) &
            '*** Warning message from routine MA57BD **', &
            '   INFO(1) =', info(1)
           call wrtlin( outbuf )
           write( outbuf, '(a,i5)' ) '     Matrix is singular, rank =', ntotpv
           call wrtlin( outbuf )
        end if
        do 3331 i = 1,n
          ppos(i) = 0
 3331   continue
        iwpos = 4
        do 3332 i = 1,nblk
          nfront = iw(iwpos)
          npiv = iw(iwpos+1)
          do 3330 j = iwpos+2,iwpos+npiv+1
            ppos(abs(iw(j))) = 1
 3330     continue
          iwpos = iwpos + nfront + 2
 3332   continue
        k= 0
        do 3333 i=1,n
          if (ppos(i).eq.0) then
            k=k+1
            nblk = nblk + 1
            nrlbdu = nrlbdu+1
            a(nrlbdu) = one
            iw(nirbdu+1) = 1
            iw(nirbdu+2) = 1
            iw(nirbdu+3) = i
            nirbdu = nirbdu+3
          endif
 3333   continue
      endif
      info(14) = nrlbdu
      iw(1) = nrlbdu + 1
      iw(2) = nrlbdu + ntwo
      info(15) = iw(2)
      iw(3) = nblk
      info(31) = nblk
      call ma57wd(a,la,iw,liw,nrlbdu)
      info(16) = nirbdu
      info(18) = intspa
      info(20) = tinspa
      info(17) = rlspa
      info(19) = trlspa
      info(21) = maxfrt
      info(22) = ntwo
      info(24) = neig
      info(25) = ntotpv
      info(28) = ncmpbr
      info(29) = ncmpbi
      rinfo(3) = flopsa
      rinfo(4) = flopsb
      rinfo(5) = flopsx
      if (info(27).gt.0) then
        rinfo(14) = zero
        do 332 i = 1,n
          rinfo(14) = max(rinfo(14),diag(i))
 332    continue
      endif
      return
      end
      subroutine ma57pd(a,iw,j1,j2,itop,real)
      integer itop,j1,j2
      logical real
      double precision a(*)
      integer iw(*)
      integer ipos,jj
      if (j2.eq.itop) go to 50
      ipos = itop - 1
      if (real) then
        do 10 jj = j2-1,j1+1,-1
          a(ipos) = a(jj)
          ipos = ipos - 1
   10   continue
      else
        do 20 jj = j2-1,j1+1,-1
          iw(ipos) = iw(jj)
          ipos = ipos - 1
   20   continue
      endif
      j2 = itop
      j1 = ipos
   50 return
      end
      subroutine ma57wd(a,la,iw,liw,nrlbdu)
      integer la,liw
      double precision a(la)
      integer iw(liw)
      integer nrlbdu
      double precision zero
      parameter (zero=0.0d0)
      integer apos,iblk,irow,iwpos,j,jpiv,ncols,nrows
      apos = 1
      iwpos = 6
      do 40 iblk = 1,iw(3)
        ncols = iw(iwpos-2)
        nrows = iw(iwpos-1)
        jpiv = 1
        do 30 irow = 1,nrows
          jpiv = jpiv - 1
          if (jpiv.eq.1) go to 10
          if (iw(iwpos+irow-1).lt.0) then
            jpiv = 2
            nrlbdu = nrlbdu + 1
            a(nrlbdu) = a(apos+1)
            a(apos+1) = zero
          end if
   10     do 20 j = apos + 1,apos + nrows - irow
            a(j) = -a(j)
   20     continue
          apos = apos + nrows - irow + 1
   30   continue
        apos = apos + nrows* (ncols-nrows)
        iwpos = iwpos + ncols + 2
   40 continue
      end
      subroutine ma57xd(n,fact,lfact,ifact,lifact,rhs,lrhs, &
                        w,lw,iw1,icntl)
      integer n,lfact
      double precision fact(lfact)
      integer lifact,ifact(lifact),lrhs,lw
      double precision w(lw),rhs(lrhs)
      integer iw1(n),icntl(20)
      intrinsic abs
      external dgemv,dtpsv
      double precision one
      parameter (one=1.0d0)
      integer apos,i,iblk,ii,ipiv,irhs,iwpos,j,j1,j2,k,k1,k2, &
              ncols,nrows
      double precision w1,w2
      apos = 1
      iwpos = 4
      do 270 iblk = 1,ifact(3)
        iw1(iblk) = iwpos
        ncols = ifact(iwpos)
        nrows = ifact(iwpos+1)
        iwpos = iwpos + 2
        if (nrows.gt.4 .and. ncols.gt.icntl(13)) then
          do 10 i = 1,ncols
            ii = abs(ifact(iwpos+i-1))
            w(i) = rhs(ii)
   10     continue
          call dtpsv('L','N','U',nrows,fact(apos),w,1)
          apos = apos + (nrows* (nrows+1))/2
          if (ncols.gt.nrows) call dgemv('N',ncols-nrows,nrows, &
                                        one,fact(apos),ncols-nrows, &
                                        w,1,one,w(nrows+1),1)
          apos = apos + nrows* (ncols-nrows)
          do 35 i = 1,ncols
            ii = abs(ifact(iwpos+i-1))
            rhs(ii) = w(i)
   35     continue
        else
        j1 = iwpos
        j2 = iwpos + nrows - 1
        do 130 ipiv = 1,nrows
          apos = apos + 1
          w1 = rhs(abs(ifact(j1)))
          k = apos
          do 100 j = j1+1,j2
            irhs = abs(ifact(j))
            rhs(irhs) = rhs(irhs) - fact(k)*w1
            k = k + 1
  100     continue
          apos = k
          j1 = j1 + 1
  130   continue
        j2 = iwpos + ncols - 1
        do 136 ipiv = 1,nrows-1,2
          k1 = apos
          k2 = apos+ncols-nrows
          w1 = rhs(abs(ifact(iwpos+ipiv-1)))
          w2 = rhs(abs(ifact(iwpos+ipiv)))
          do 133 j = j1,j2
            irhs = abs(ifact(j))
            rhs(irhs) = rhs(irhs) + w1*fact(k1) + w2*fact(k2)
            k1 = k1 + 1
            k2 = k2 + 1
  133     continue
          apos = k2
  136   continue
        if (mod(nrows,2).eq.1) then
          k = apos
          w1 = rhs(abs(ifact(iwpos+ipiv-1)))
          do 137 j = j1,j2
            irhs = abs(ifact(j))
            rhs(irhs) = rhs(irhs) + w1*fact(k)
            k = k + 1
  137     continue
          apos = k
        endif
      end if
      iwpos = iwpos + ncols
  270 continue
      end
      subroutine ma57yd(n,fact,lfact,ifact,lifact,rhs,lrhs, &
                        w,lw,iw1,icntl)
      integer n,lfact
      double precision fact(lfact)
      integer lifact,ifact(lifact),lrhs,lw
      double precision w(lw),rhs(lrhs)
      integer iw1(n),icntl(20)
      intrinsic abs
      external dgemv,dtpsv
      double precision one
      parameter (one=1.0d0)
      integer apos,apos2,i,iblk,ii,ipiv,irhs,irhs1, &
              irhs2,iwpos,j,jpiv,j1,j2,k,k2,lrow,ncols,nrows
      double precision w1,w2
      apos = ifact(1)
      apos2 = ifact(2)
      do 380 iblk = ifact(3),1,-1
        iwpos = iw1(iblk)
        ncols = abs(ifact(iwpos))
        nrows = abs(ifact(iwpos+1))
        apos = apos - nrows* (ncols-nrows)
        iwpos = iwpos + 2
        if (nrows.gt.4 .and. ncols.gt.icntl(13)) then
          do 5 i = nrows + 1,ncols
            ii = abs(ifact(iwpos+i-1))
            w(i) = rhs(ii)
    5     continue
          do 10 ipiv = nrows,1,-1
            irhs = abs(ifact(iwpos+ipiv-1))
            apos = apos - (nrows+1-ipiv)
            w(ipiv) = rhs(irhs)*fact(apos)
   10     continue
          jpiv = -1
          do 20 ipiv = nrows,1,-1
            irhs = ifact(iwpos+ipiv-1)
            if (irhs.lt.0) then
              irhs1 = -ifact(iwpos+ipiv-1+jpiv)
              w(ipiv) = rhs(irhs1)*fact(apos2) + w(ipiv)
              if (jpiv.eq.1) apos2 = apos2 - 1
              jpiv = -jpiv
            end if
   20     continue
          k = ncols - nrows
          if (k.gt.0) call dgemv('T',k,nrows,one, &
                                 fact(apos+(nrows*(nrows+1))/2),k, &
                                 w(nrows+1),1,one,w,1)
          call dtpsv('L','T','U',nrows,fact(apos),w,1)
          do 60 i = 1,nrows
            ii = abs(ifact(iwpos+i-1))
            rhs(ii) = w(i)
   60     continue
        else
          j1 = iwpos
          j2 = iwpos + ncols - 1
          jpiv = -1
          do 210 ipiv = nrows,1,-1
            irhs = ifact(iwpos+ipiv-1)
            lrow = nrows + 1 - ipiv
            if (irhs.gt.0) then
              apos = apos - lrow
              rhs(irhs) = rhs(irhs)*fact(apos)
            else
              if (jpiv.eq.-1) then
                irhs1 = -ifact(iwpos+ipiv-2)
                irhs2 = -irhs
                apos = apos - lrow - lrow - 1
                w1 = rhs(irhs1)*fact(apos) + &
                     rhs(irhs2)*fact(apos2)
                rhs(irhs2) = rhs(irhs1)*fact(apos2) + &
                               rhs(irhs2)*fact(apos+lrow+1)
                rhs(irhs1) = w1
                apos2 = apos2 - 1
              end if
              jpiv = -jpiv
            end if
  210     continue
          apos = apos + (nrows* (nrows+1))/2
          k = apos
          j1 = iwpos + nrows
          do 220 ipiv = 1,nrows-1,2
            irhs = abs(ifact(iwpos+ipiv-1))
            w1 = rhs(irhs)
            irhs1 = abs(ifact(iwpos+ipiv))
            w2 = rhs(irhs1)
            k2 = k+(ncols-nrows)
            do 215 j = j1,j2
              ii = abs(ifact(j))
              w1 = w1 + fact(k)*rhs(ii)
              w2 = w2 + fact(k2)*rhs(ii)
              k = k + 1
              k2 = k2 + 1
  215       continue
            rhs(irhs) = w1
            rhs(irhs1) = w2
            k = k2
  220     continue
          if (mod(nrows,2).eq.1) then
            irhs = abs(ifact(iwpos+ipiv-1))
            w1 = rhs(irhs)
            do 216 j = j1,j2
              w1 = w1 + fact(k)*rhs(abs(ifact(j)))
              k = k + 1
  216       continue
            rhs(irhs) = w1
          endif
          j2 = iwpos + nrows - 1
          do 260 ipiv = 1,nrows
            irhs = abs(ifact(j1-1))
            apos = apos - ipiv
            w1 = rhs(irhs)
            k = apos + 1
            do 230 j = j1,j2
              w1 = w1 - fact(k)*rhs(abs(ifact(j)))
              k = k + 1
  230       continue
            rhs(irhs) = w1
            j1 = j1 - 1
  260     continue
        end if
  380 continue
      end
      subroutine ma57vd(n,nz,irn,icn,iw,lw,ipe,iq,flag,iwfr, &
                       icntl,info)
      integer iwfr,lw,n,nz
      integer flag(n),icn(*),ipe(n),iq(n),irn(*),iw(lw)
      integer icntl(*),info(*)
      integer i,id,j,jn,k,k1,k2,l,last,lr,n1,ndup
      character*80 outbuf
      info(2) = 0
      do 10 i = 1,n
        ipe(i) = 0
   10 continue
      lr = nz
      if (nz.eq.0) go to 120
      do 110 k = 1,nz
        i = irn(k)
        j = icn(k)
        if (i.lt.j) then
          if (i.ge.1 .and. j.le.n) go to 90
        else if (i.gt.j) then
          if (j.ge.1 .and. i.le.n) go to 90
        else
          if (i.ge.1 .and. i.le.n) go to 80
        end if
        info(2) = info(2) + 1
        info(1) = 1
        if (info(2).le.1 .and. icntl(2).gt.0) then
           call wrtlin( '    *** HSL WARNING ***' )
           write( outbuf, '(2a,i2)' ) &
            ' *** WARNING MESSAGE FROM SUBROUTINE MA57AD', &
            '  *** INFO(1) =', info(1)
           call wrtlin( outbuf )
        end if
        if (info(2).le.10 .and. icntl(2).gt.0) then
           call wrtlin( '    *** HSL WARNING ***' )
           write( outbuf, '(3(i6,a))' ) k, 'TH NON-ZERO (IN ROW', i, &
            ' AND COLUMN', j, ') IGNORED'
           call wrtlin( outbuf )
        end if
   80   i = 0
        j = 0
        go to 100
   90   ipe(i) = ipe(i) + 1
        ipe(j) = ipe(j) + 1
  100   iw(k) = j
        lr = lr + 1
        iw(lr) = i
  110 continue
  120 iq(1) = 1
      n1 = n - 1
      if (n1.le.0) go to 140
      do 130 i = 1,n1
        flag(i) = 0
        if (ipe(i).eq.0) ipe(i) = -1
        iq(i+1) = ipe(i) + iq(i) + 1
        ipe(i) = iq(i)
  130 continue
  140 last = ipe(n) + iq(n)
      flag(n) = 0
      if (lr.ge.last) go to 160
      k1 = lr + 1
      do 150 k = k1,last
        iw(k) = 0
  150 continue
  160 ipe(n) = iq(n)
      iwfr = last + 1
      if (nz.eq.0) go to 230
      do 220 k = 1,nz
        j = iw(k)
        if (j.le.0) go to 220
        l = k
        iw(k) = 0
        do 210 id = 1,nz
          if (l.gt.nz) go to 170
          l = l + nz
          go to 180
  170     l = l - nz
  180     i = iw(l)
          iw(l) = 0
          if (i.lt.j) go to 190
          l = iq(j) + 1
          iq(j) = l
          jn = iw(l)
          iw(l) = -i
          go to 200
  190     l = iq(i) + 1
          iq(i) = l
          jn = iw(l)
          iw(l) = -j
  200     j = jn
          if (j.le.0) go to 220
  210   continue
  220 continue
  230 ndup = 0
      do 280 i = 1,n
        k1 = ipe(i) + 1
        k2 = iq(i)
        if (k1.le.k2) go to 240
        ipe(i) = 0
        iq(i) = 0
        go to 280
  240   do 260 k = k1,k2
          j = -iw(k)
          if (j.le.0) go to 270
          l = iq(j) + 1
          iq(j) = l
          iw(l) = i
          iw(k) = j
          if (flag(j).ne.i) go to 250
          ndup = ndup + 1
          iw(l) = 0
          iw(k) = 0
  250     flag(j) = i
  260   continue
  270   iq(i) = iq(i) - ipe(i)
        if (ndup.eq.0) iw(k1-1) = iq(i)
  280 continue
      if (ndup.eq.0) go to 310
      iwfr = 1
      do 300 i = 1,n
        k1 = ipe(i) + 1
        if (k1.eq.1) go to 300
        k2 = iq(i) + ipe(i)
        l = iwfr
        ipe(i) = iwfr
        iwfr = iwfr + 1
        do 290 k = k1,k2
          if (iw(k).eq.0) go to 290
          iw(iwfr) = iw(k)
          iwfr = iwfr + 1
  290   continue
        iw(l) = iwfr - l - 1
  300 continue
  310 return
      end
      subroutine ma57hd(n,ipe,iw,lw,iwfr,nv,nxt,lst,ipd,flag,iovflo, &
                       ncmpa,fratio)
      double precision fratio
      integer iwfr,lw,n,iovflo,ncmpa
      integer flag(n),ipd(n),ipe(n),iw(lw),lst(n),nv(n),nxt(n)
      integer i,id,idl,idn,ie,ip,is,jp,jp1,jp2,js,k,k1,k2,ke,kp,kp0,kp1, &
              kp2,ks,l,len,limit,ln,ls,lwfr,md,me,ml,ms,nel,nflg,np, &
              np0,ns,nvpiv,nvroot,root
      external ma57zd
      intrinsic abs,min
      do 10 i = 1,n
        ipd(i) = 0
        nv(i) = 1
        flag(i) = iovflo
   10 continue
      md = 1
      ncmpa = 0
      nflg = iovflo
      nel = 0
      root = n+1
      nvroot = 0
      do 30 is = 1,n
        k = ipe(is)
        if (k.le.0) go to 20
        id = iw(k) + 1
        ns = ipd(id)
        if (ns.gt.0) lst(ns) = is
        nxt(is) = ns
        ipd(id) = is
        lst(is) = -id
        go to 30
   20   nel = nel + 1
        flag(is) = -1
        nxt(is) = 0
        lst(is) = 0
   30 continue
      do 340 ml = 1,n
        if (nel+nvroot+1.ge.n) go to 350
        do 40 id = md,n
          ms = ipd(id)
          if (ms.gt.0) go to 50
   40   continue
   50   md = id
        nvpiv = nv(ms)
        ns = nxt(ms)
        nxt(ms) = 0
        lst(ms) = 0
        if (ns.gt.0) lst(ns) = -id
        ipd(id) = ns
        me = ms
        nel = nel + nvpiv
        idn = 0
        kp = ipe(me)
        flag(ms) = -1
        ip = iwfr
        len = iw(kp)
        do 140 kp1 = 1,len
          kp = kp + 1
          ke = iw(kp)
          if (flag(ke).le.-2) go to 60
          if (flag(ke).le.0) then
             if (ipe(ke).ne.-root) go to 140
             ke = root
             if (flag(ke).le.0) go to 140
          end if
          jp = kp - 1
          ln = len - kp1 + 1
          ie = ms
          go to 70
   60     ie = ke
          jp = ipe(ie)
          ln = iw(jp)
   70     do 130 jp1 = 1,ln
            jp = jp + 1
            is = iw(jp)
            if (flag(is).le.0) then
               if (ipe(is).eq.-root) then
                  is = root
                  iw(jp) = root
                  if (flag(is).le.0) go to 130
               else
                  go to 130
               end if
            end if
            flag(is) = 0
            if (iwfr.lt.lw) go to 100
            ipe(ms) = kp
            iw(kp) = len - kp1
            ipe(ie) = jp
            iw(jp) = ln - jp1
            call ma57zd(n,ipe,iw,ip-1,lwfr,ncmpa)
            jp2 = iwfr - 1
            iwfr = lwfr
            if (ip.gt.jp2) go to 90
            do 80 jp = ip,jp2
              iw(iwfr) = iw(jp)
              iwfr = iwfr + 1
   80       continue
   90       ip = lwfr
            jp = ipe(ie)
            kp = ipe(me)
  100       iw(iwfr) = is
            idn = idn + nv(is)
            iwfr = iwfr + 1
            ls = lst(is)
            lst(is) = 0
            ns = nxt(is)
            nxt(is) = 0
            if (ns.gt.0) lst(ns) = ls
            if (ls.lt.0) then
              ls = -ls
              ipd(ls) = ns
            else if (ls.gt.0) then
              nxt(ls) = ns
            end if
  130     continue
          if (ie.eq.ms) go to 150
          ipe(ie) = -me
          flag(ie) = -1
  140   continue
  150   nv(ms) = idn + nvpiv
        if (iwfr.eq.ip) go to 330
        k1 = ip
        k2 = iwfr - 1
        limit = nint(fratio*(n-nel))
        do 310 k = k1,k2
          is = iw(k)
          if (is.eq.root) go to 310
          if (nflg.gt.2) go to 170
          do 160 i = 1,n
            if (flag(i).gt.0) flag(i) = iovflo
            if (flag(i).le.-2) flag(i) = -iovflo
  160     continue
          nflg = iovflo
  170     nflg = nflg - 1
          id = idn
          kp1 = ipe(is) + 1
          np = kp1
          kp2 = iw(kp1-1) + kp1 - 1
          do 220 kp = kp1,kp2
            ke = iw(kp)
          if (flag(ke).eq.-1) then
             if (ipe(ke).ne.-root) go to 220
             ke = root
             iw(kp) = root
             if (flag(ke).eq.-1) go to 220
          end if
          if (flag(ke).ge.0) go to 230
            jp1 = ipe(ke) + 1
            jp2 = iw(jp1-1) + jp1 - 1
            idl = id
            do 190 jp = jp1,jp2
              js = iw(jp)
              if (flag(js).le.nflg) go to 190
              id = id + nv(js)
              flag(js) = nflg
  190       continue
            if (id.gt.idl) go to 210
            do 200 jp = jp1,jp2
              js = iw(jp)
              if (flag(js).ne.0) go to 210
  200       continue
            ipe(ke) = -me
            flag(ke) = -1
            go to 220
  210       iw(np) = ke
            flag(ke) = -nflg
            np = np + 1
  220     continue
          np0 = np
          go to 250
  230     kp0 = kp
          np0 = np
          do 240 kp = kp0,kp2
            ks = iw(kp)
            if (flag(ks).le.nflg) then
               if (ipe(ks).eq.-root) then
                  ks = root
                  iw(kp) = root
                  if (flag(ks).le.nflg) go to 240
               else
                  go to 240
               end if
            end if
            id = id + nv(ks)
            flag(ks) = nflg
            iw(np) = ks
            np = np + 1
  240     continue
  250     if (id.ge.limit) go to 295
          iw(np) = iw(np0)
          iw(np0) = iw(kp1)
          iw(kp1) = me
          iw(kp1-1) = np - kp1 + 1
          js = ipd(id)
          do 280 l = 1,n
            if (js.le.0) go to 300
            kp1 = ipe(js) + 1
            if (iw(kp1).ne.me) go to 300
            kp2 = kp1 - 1 + iw(kp1-1)
            do 260 kp = kp1,kp2
              ie = iw(kp)
              if (abs(flag(ie)+0).gt.nflg) go to 270
  260       continue
            go to 290
  270       js = nxt(js)
  280     continue
  290     ipe(js) = -is
          nv(is) = nv(is) + nv(js)
          nv(js) = 0
          flag(js) = -1
          ns = nxt(js)
          ls = lst(js)
          if (ns.gt.0) lst(ns) = is
          if (ls.gt.0) nxt(ls) = is
          lst(is) = ls
          nxt(is) = ns
          lst(js) = 0
          nxt(js) = 0
          if (ipd(id).eq.js) ipd(id) = is
          go to 310
  295     if (nvroot.eq.0) then
            root = is
            ipe(is) = 0
          else
            iw(k) = root
            ipe(is) = -root
            nv(root) = nv(root) + nv(is)
            nv(is) = 0
            flag(is) = -1
          end if
          nvroot = nv(root)
          go to 310
  300     ns = ipd(id)
          if (ns.gt.0) lst(ns) = is
          nxt(is) = ns
          ipd(id) = is
          lst(is) = -id
          md = min(md,id)
  310   continue
        do 320 k = k1,k2
          is = iw(k)
          if (nv(is).eq.0) go to 320
          flag(is) = nflg
          iw(ip) = is
          ip = ip + 1
  320   continue
        iwfr = k1
        flag(me) = -nflg
        iw(ip) = iw(k1)
        iw(k1) = ip - k1
        ipe(me) = k1
        iwfr = ip + 1
        go to 335
  330   ipe(me) = 0
  335   continue
  340 continue
  350 do 360 is = 1,n
        if(nxt(is).ne.0 .or. lst(is).ne.0) then
          if (nvroot.eq.0) then
            root = is
            ipe(is) = 0
          else
            ipe(is) = -root
          end if
          nvroot = nvroot + nv(is)
          nv(is) = 0
         end if
  360 continue
      do 370 ie = 1,n
        if (ipe(ie).gt.0) ipe(ie) = -root
  370 continue
      if(nvroot.gt.0)nv(root)=nvroot
      end
      subroutine ma57zd(n,ipe,iw,lw,iwfr,ncmpa)
      integer iwfr,lw,n,ncmpa
      integer ipe(n),iw(lw)
      integer i,ir,k,k1,k2,lwfr
      ncmpa = ncmpa + 1
      do 10 i = 1,n
        k1 = ipe(i)
        if (k1.le.0) go to 10
        ipe(i) = iw(k1)
        iw(k1) = -i
   10 continue
      iwfr = 1
      lwfr = iwfr
      do 60 ir = 1,n
        if (lwfr.gt.lw) go to 70
        do 20 k = lwfr,lw
          if (iw(k).lt.0) go to 30
   20   continue
        go to 70
   30   i = -iw(k)
        iw(iwfr) = ipe(i)
        ipe(i) = iwfr
        k1 = k + 1
        k2 = k + iw(iwfr)
        iwfr = iwfr + 1
        if (k1.gt.k2) go to 50
        do 40 k = k1,k2
          iw(iwfr) = iw(k)
          iwfr = iwfr + 1
   40   continue
   50   lwfr = k2 + 1
   60 continue
   70 return
      end
!dense    Version 1
!AMD      SUBROUTINE AMDD (N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
!AMD     $                   LAST, NCMPA, DEGREE, HEAD, NEXT, W)
      subroutine mc50bd (thresh, n, iwlen, pe, pfree, len, iw, nv, &
                         elen, last, ncmpa, degree, head, next, w)
!dense
      integer n, iwlen, pe(n), pfree, len(n), iw(iwlen), nv(n), &
              elen(n), last(n), ncmpa, degree(n), head(n), next(n), &
              w(n),ndense(n)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!dense
      integer thresh
      integer thresm, minden, maxden, ndme
      integer nbd,nbed, nbdm, lastd, nelme
      logical idense
      double precision relden
!dense
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!dense
!dense
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      integer deg, degme, dext, dmax, e, elenme, eln, hash, hmod, i, &
              ilast, inext, j, jlast, jnext, k, knt1, knt2, knt3, &
              lenj, ln, maxmem, me, mem, mindeg, nel, newmem, &
              nleft, nvi, nvj, nvpiv, slenme, we, wflg, wnvi, x
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      integer p, p1, p2, p3, pdst, pend, pj, pme, pme1, pme2, pn, psrc
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      intrinsic max, min, mod
!=======================================================================
!=======================================================================
!dense
      if (thresh.gt.0) then
         thresm  = 0
!
         relden = 0.0
         do i=1,n
             relden = relden + float(len(i))/float(n)
             thresm = max(thresm, len(i))
          enddo
         thresm =  int(relden)*10 + (thresm-int(relden))/10 + 1
      else
         thresm = thresh
      endif
      if (thresm.ge.0) then
       if ((thresm.ge.n).or.(thresm.lt.2)) then
          thresm = n
       endif
      endif
      lastd = 0
      nbd   = 0
      nbed  = 0
      nbdm  = 0
!dense
      wflg = 2
      mindeg = 1
      ncmpa = 0
      nel = 0
      hmod = max (1, n-1)
      dmax = 0
      mem = pfree - 1
      maxmem = mem
      do 10 i = 1, n
!dense
        ndense(i)= 0
!dense
        last (i) = 0
        head (i) = 0
        nv (i) = 1
        w (i) = 1
        elen (i) = 0
        degree (i) = len (i)
   10 continue
      do 20 i = 1, n
        deg = degree (i)
        if (deg .gt. 0) then
!dense
          if ( (thresm.ge.0) .and. &
               (deg+1.ge.thresm) ) then
            nbd = nbd+1
            if (deg+1.ne.n-nel) then
             degree(i) = degree(i)+n+1
             deg = n
             inext = head (deg)
             if (inext .ne. 0) last (inext) = i
             next (i) = inext
             head (deg) = i
             last(i)  = 0
             if (lastd.eq.0) lastd=i
            else
             nbed = nbed+1
             degree(i) = n+1
             deg = n
             if (lastd.eq.0) then
               lastd     = i
               head(deg) = i
               next(i)   = 0
               last(i)   = 0
             else
               next(lastd) = i
               last(i)     = lastd
               lastd       = i
               next(i)     = 0
             endif
            endif
          else
!dense
            inext = head (deg)
            if (inext .ne. 0) last (inext) = i
            next (i) = inext
            head (deg) = i
          endif
        else
          nel = nel + 1
          elen (i) = -nel
          pe (i) = 0
          w (i) = 0
        endif
   20 continue
          if (nbd.eq.0) thresm = n
!=======================================================================
!=======================================================================
   30 if (nel .lt. n) then
!=======================================================================
!=======================================================================
        do 40 deg = mindeg, n
          me = head (deg)
          if (me .gt. 0) go to 50
   40   continue
   50   mindeg = deg
!dense
        if (deg.lt.n)  then
!dense
          inext = next (me)
          if (inext .ne. 0) last (inext) = 0
          head (deg) = inext
!dense
        else
          nbdm = max(nbdm,nbd)
          if (degree(me).gt.n+1) then
            minden = nbd
            maxden = 0
            if (wflg+nbd+1 .le. wflg) then
             do  52 x = 1, n
              if (w (x) .ne. 0) w (x) = 1
  52         continue
             wflg = 2
            endif
            wflg = wflg + 1
  51        continue
            inext = next (me)
            if (inext .ne. 0) then
               last (inext) = 0
            else
               lastd = 0
            endif
            ndense(me) = 0
            w(me)      = wflg
            p1 = pe(me)
            p2 = p1 + len(me) -1
            ln       = p1
            eln      = p1
            do 55 p=p1,p2
              e= iw(p)
              if (w(e).eq.wflg) goto 55
              w(e) = wflg
              if (pe(e).lt.0) then
                x = e
  53            x = -pe(x)
                if (w(x) .eq.wflg) goto 55
                w(x) = wflg
                if ( pe(x) .lt. 0 ) goto 53
                e = x
              endif
              if (elen(e).lt.0) then
               ndense(e) = ndense(e) - nv(me)
               iw(ln) = iw(eln)
               iw(eln) = e
               ln  = ln+1
               eln = eln + 1
               pme1 = pe(e)
               do 54 pme = pme1, pme1+len(e)-1
                x = iw(pme)
                if ((elen(x).ge.0).and.(w(x).ne.wflg)) then
                 ndense(me) = ndense(me) + nv(x)
                 w(x) = wflg
                endif
 54            continue
              else
               ndense(me) = ndense(me) + nv(e)
               iw(ln)=e
               ln = ln+1
              endif
  55        continue
            wflg     = wflg + 1
            len(me)  = ln-p1
            elen(me) = eln- p1
            ndme = ndense(me)+nv(me)
            minden = min (minden, ndme)
            maxden = max (maxden, ndme)
            if (ndense(me).eq.0) ndense(me) =1
            degree(me) = ndense(me)
            deg = degree(me)
            mindeg = min(deg,mindeg)
            jnext = head(deg)
            if (jnext.ne. 0) last (jnext) = me
            next(me) = jnext
            head(deg) = me
            me    = inext
            if (me.ne.0) then
              if (degree(me).gt.(n+1) ) goto 51
            endif
            head (n) = me
            thresm = max(thresm*2, minden+(maxden-minden)/2)
            thresm = min(thresm,nbd)
            if (thresm.ge.nbd) thresm=n
            nbd    = nbed
!
            goto 30
          endif
          if (degree(me).eq.n+1) then
!CC Diagnostic printing and stop removed
!CC        IF (NBD.NE.NBED) THEN
!CC         write(6,*) ' Error -1 quasi dense rows remains'
!CC         stop
!CC        ENDIF
           nelme    = -(nel+1)
           do 59 x=1,n
            if ((pe(x).gt.0) .and. (elen(x).lt.0)) then
             pe(x) = -me
            elseif (degree(x).eq.n+1) then
             nel   = nel + nv(x)
             pe(x) = -me
             elen(x) = 0
             nv(x) = 0
            endif
   59      continue
           elen(me) = nelme
           nv(me)   = nbd
           pe(me)   = 0
!CC Diagnostic printing and stop removed
!CC        IF (NEL.NE.N) THEN
!CC         write(6,*) 'ERROR 2 detected in AMDD'
!CC         write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
!CC         stop
!CC        ENDIF
           goto 265
          endif
        endif
!dense
!dense traces
        elenme = elen (me)
        elen (me) = - (nel + 1)
        nvpiv = nv (me)
        nel = nel + nvpiv
!dense
        ndense(me) = 0
!dense
!=======================================================================
!=======================================================================
        nv (me) = -nvpiv
        degme = 0
        if (elenme .eq. 0) then
          pme1 = pe (me)
          pme2 = pme1 - 1
          do 60 p = pme1, pme1 + len (me) - 1
            i = iw (p)
            nvi = nv (i)
            if (nvi .gt. 0) then
              degme = degme + nvi
              nv (i) = -nvi
              pme2 = pme2 + 1
              iw (pme2) = i
!dense
              if (degree(i).le.n) then
!dense
              ilast = last (i)
              inext = next (i)
              if (inext .ne. 0) last (inext) = ilast
              if (ilast .ne. 0) then
                next (ilast) = inext
              else
                head (degree (i)) = inext
              endif
!dense
              else
               ndense(me) = ndense(me) + nvi
              endif
!dense
            endif
   60     continue
          newmem = 0
        else
          p = pe (me)
          pme1 = pfree
          slenme = len (me) - elenme
          do 120 knt1 = 1, elenme + 1
            if (knt1 .gt. elenme) then
              e = me
              pj = p
              ln = slenme
            else
              e = iw (p)
              p = p + 1
              pj = pe (e)
              ln = len (e)
            endif
            do 110 knt2 = 1, ln
              i = iw (pj)
              pj = pj + 1
              nvi = nv (i)
              if (nvi .gt. 0) then
                if (pfree .gt. iwlen) then
                  pe (me) = p
                  len (me) = len (me) - knt1
                  if (len (me) .eq. 0) pe (me) = 0
                  pe (e) = pj
                  len (e) = ln - knt2
                  if (len (e) .eq. 0) pe (e) = 0
                  ncmpa = ncmpa + 1
                  do 70 j = 1, n
                    pn = pe (j)
                    if (pn .gt. 0) then
                      pe (j) = iw (pn)
                      iw (pn) = -j
                    endif
   70             continue
                  pdst = 1
                  psrc = 1
                  pend = pme1 - 1
   80             continue
                  if (psrc .le. pend) then
                    j = -iw (psrc)
                    psrc = psrc + 1
                    if (j .gt. 0) then
                      iw (pdst) = pe (j)
                      pe (j) = pdst
                      pdst = pdst + 1
                      lenj = len (j)
                      do 90 knt3 = 0, lenj - 2
                        iw (pdst + knt3) = iw (psrc + knt3)
   90                 continue
                      pdst = pdst + lenj - 1
                      psrc = psrc + lenj - 1
                    endif
                    go to 80
                  endif
                  p1 = pdst
                  do 100 psrc = pme1, pfree - 1
                    iw (pdst) = iw (psrc)
                    pdst = pdst + 1
  100             continue
                  pme1 = p1
                  pfree = pdst
                  pj = pe (e)
                  p = pe (me)
                endif
                degme = degme + nvi
                nv (i) = -nvi
                iw (pfree) = i
                pfree = pfree + 1
!dense
                if (degree(i).le.n) then
!dense
                ilast = last (i)
                inext = next (i)
                if (inext .ne. 0) last (inext) = ilast
                if (ilast .ne. 0) then
                  next (ilast) = inext
                else
                  head (degree (i)) = inext
                endif
!dense
                else
                 ndense(me) = ndense(me) + nvi
                endif
!dense
              endif
  110       continue
            if (e .ne. me) then
              pe (e) = -me
              w (e) = 0
            endif
  120     continue
          pme2 = pfree - 1
          newmem = pfree - pme1
          mem = mem + newmem
          maxmem = max (maxmem, mem)
        endif
        degree (me) = degme
        pe (me) = pme1
        len (me) = pme2 - pme1 + 1
        if (wflg+n .le. wflg) then
          do 130 x = 1, n
            if (w (x) .ne. 0) w (x) = 1
  130     continue
          wflg = 2
        endif
!=======================================================================
!dense
!dense
!=======================================================================
!dense
!dense
        do 150 pme = pme1, pme2
          i = iw (pme)
!dense
          if (degree(i).gt.n) goto 150
!dense
          eln = elen (i)
          if (eln .gt. 0) then
            nvi = -nv (i)
            wnvi = wflg - nvi
            do 140 p = pe (i), pe (i) + eln - 1
              e = iw (p)
              we = w (e)
              if (we .ge. wflg) then
                we = we - nvi
              else if (we .ne. 0) then
!dense
                we = degree (e) + wnvi - ndense(e)
!dense
              endif
              w (e) = we
  140       continue
          endif
  150   continue
!=======================================================================
!=======================================================================
!dense Lme should be read Lme(G')
        do 180 pme = pme1, pme2
          i = iw (pme)
!dense
          if (degree(i).gt.n) goto 180
!dense
          p1 = pe (i)
          p2 = p1 + elen (i) - 1
          pn = p1
          hash = 0
          deg = 0
          do 160 p = p1, p2
            e = iw (p)
            dext = w (e) - wflg
            if (dext .gt. 0) then
              deg = deg + dext
              iw (pn) = e
              pn = pn + 1
              hash = hash + e
!AMD            ELSE IF (DEXT .EQ. 0) THEN
!AMDC             aggressive absorption: e is not adjacent to me, but
!AMDC             the |Le \ Lme| is 0, so absorb it into me
!AMD              PE (E) = -ME
!AMD              W (E) = 0
!dense
              else if ((dext .eq. 0) .and. &
                      (ndense(me).eq.nbd)) then
                pe (e) = -me
                w (e)  = 0
              else if (dext.eq.0) then
                  iw(pn) = e
                  pn     = pn+1
                  hash   = hash + e
!dense
            endif
  160     continue
          elen (i) = pn - p1 + 1
          p3 = pn
          do 170 p = p2 + 1, p1 + len (i) - 1
            j = iw (p)
            nvj = nv (j)
            if (nvj .gt. 0) then
!AMD              DEG = DEG + NVJ
!dense
              if (degree(j).le.n) deg=deg+nvj
!dense
              iw (pn) = j
              pn = pn + 1
              hash = hash + j
            endif
  170     continue
!AMD          IF (DEG .EQ. 0) THEN
!dense
          if ((deg .eq. 0).and.(ndense(me).eq.nbd)) then
!dense
            pe (i) = -me
            nvi = -nv (i)
            degme = degme - nvi
            nvpiv = nvpiv + nvi
            nel = nel + nvi
            nv (i) = 0
            elen (i) = 0
          else
!AMD            DEGREE (I) = MIN (DEGREE (I), DEG)
!AMD        modified test moved to loop 260
!dense
            degree(i) = min (deg+nbd-ndense(me), &
                             degree(i))
!dense
            iw (pn) = iw (p3)
            iw (p3) = iw (p1)
            iw (p1) = me
            len (i) = pn - p1 + 1
            hash = mod (hash, hmod) + 1
            j = head (hash)
            if (j .le. 0) then
              next (i) = -j
              head (hash) = -i
            else
              next (i) = last (j)
              last (j) = i
            endif
            last (i) = hash
          endif
  180   continue
        degree (me) = degme
        dmax = max (dmax, degme)
        wflg = wflg + dmax
        if (wflg+n .le. wflg) then
          do 190 x = 1, n
            if (w (x) .ne. 0) w (x) = 1
  190     continue
          wflg = 2
        endif
!=======================================================================
!=======================================================================
        do 250 pme = pme1, pme2
          i = iw (pme)
!AMD          IF (NV (I) .LT. 0) THEN
!dense
          if ( (nv(i).lt.0) .and. (degree(i).le.n) ) then
!dense
            hash = last (i)
            j = head (hash)
            if (j .eq. 0) go to 250
            if (j .lt. 0) then
              i = -j
              head (hash) = 0
            else
              i = last (j)
              last (j) = 0
            endif
            if (i .eq. 0) go to 250
  200       continue
            if (next (i) .ne. 0) then
              ln = len (i)
              eln = elen (i)
              do 210 p = pe (i) + 1, pe (i) + ln - 1
                w (iw (p)) = wflg
  210         continue
              jlast = i
              j = next (i)
  220         continue
              if (j .ne. 0) then
                if (len (j) .ne. ln) go to 240
                if (elen (j) .ne. eln) go to 240
                do 230 p = pe (j) + 1, pe (j) + ln - 1
                  if (w (iw (p)) .ne. wflg) go to 240
  230           continue
                pe (j) = -i
                nv (i) = nv (i) + nv (j)
                nv (j) = 0
                elen (j) = 0
                j = next (j)
                next (jlast) = j
                go to 220
  240           continue
                jlast = j
                j = next (j)
              go to 220
              endif
              wflg = wflg + 1
              i = next (i)
              if (i .ne. 0) go to 200
            endif
          endif
  250   continue
!=======================================================================
!=======================================================================
        p = pme1
        nleft = n - nel
        do 260 pme = pme1, pme2
          i = iw (pme)
          nvi = -nv (i)
          if (nvi .gt. 0) then
            nv (i) = nvi
!dense
            if (degree(i).le.n) then
!dense
            deg = min (degree (i)+ degme - nvi, nleft - nvi)
            degree (i) = deg
            idense = .false.
            if (thresm.ge.0) then
             if (deg+nvi .ge. thresm) then
              if (thresm.eq.n) then
               if ((elen(i).le.2) .and. ((deg+nvi).eq.nleft) ) then
                degree(i) = n+1
                idense = .true.
               endif
              else
               idense = .true.
               if ((elen(i).le.2).and.((deg+nvi).eq.nleft) ) then
                 degree(i) = n+1
               else
                 degree(i) = n+1+degree(i)
               endif
              endif
             endif
             if (idense) then
               p1 = pe(i)
               p2 = p1 + elen(i) - 1
               if (p2.ge.p1) then
               do 264 pj=p1,p2
                 e= iw(pj)
                 ndense (e) = ndense(e) + nvi
 264           continue
               endif
               nbd = nbd+nvi
               deg = n
               if (degree(i).eq.n+1) then
                nbed = nbed +nvi
                if (lastd.eq.0) then
                 lastd     = i
                 head(deg) = i
                 next(i)   = 0
                 last(i)   = 0
                else
                 next(lastd) = i
                 last(i)     = lastd
                 lastd       = i
                 next(i)     = 0
                endif
               else
                inext = head(deg)
                if (inext .ne. 0) last (inext) = i
                next (i) = inext
                head (deg) = i
                last(i)    = 0
                if (lastd.eq.0) lastd=i
               endif
             endif
            endif
            if (.not.idense) then
!dense
            inext = head (deg)
            if (inext .ne. 0) last (inext) = i
            next (i) = inext
            last (i) = 0
            head (deg) = i
!dense
            endif
!dense
            mindeg = min (mindeg, deg)
!dense
            endif
!dense
            iw (p) = i
            p = p + 1
          endif
  260   continue
!=======================================================================
!=======================================================================
        nv (me) = nvpiv + degme
        len (me) = p - pme1
        if (len (me) .eq. 0) then
          pe (me) = 0
          w (me) = 0
        endif
        if (newmem .ne. 0) then
          pfree = p
          mem = mem - newmem + len (me)
        endif
!dense
!=======================================================================
      go to 30
      endif
!=======================================================================
!dense
  265 continue
!dense
!=======================================================================
!=======================================================================
      do 290 i = 1, n
        if (elen (i) .eq. 0) then
          j = -pe (i)
  270     continue
            if (elen (j) .ge. 0) then
              j = -pe (j)
              go to 270
            endif
            e = j
            k = -elen (e)
            j = i
  280       continue
            if (elen (j) .ge. 0) then
              jnext = -pe (j)
              pe (j) = -e
              if (elen (j) .eq. 0) then
                elen (j) = k
                k = k + 1
              endif
              j = jnext
            go to 280
            endif
          elen (e) = -k
        endif
  290 continue
      do 300 i = 1, n
        k = abs (elen (i))
        last (k) = i
        elen (i) = k
  300 continue
!=======================================================================
!=======================================================================
      pfree = maxmem
      return
      end

