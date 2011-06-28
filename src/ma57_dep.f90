! *******************************************************************
! COPYRIGHT (c) 1995 Timothy A. Davis, Patrick Amestoy and
!             Council for the Central Laboratory of the Research Councils
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
! Original date 30 November 1995
!  April 2001: call to MC49 changed to MC59 to make routine threadsafe
! 20/2/02 Cosmetic changes applied to reduce single/double differences

! 12th July 2004 Version 1.0.0. Version numbering added.

      subroutine mc47ad(n, ne, pe, iw, iwlen, mp, info)
      integer n, ne, pe(n+1), iwlen, iw(iwlen), mp, info(7)
      integer degree
      double precision dummy(1)
      integer elen,head,i,iflag,ii,i1,i2,j,last,len,leniw,ncmpa, &
              next,nv,pfree,w
      integer ict59(10),info59(10),iout,jout,idup,nzout
      external mc59ad,mc34ad,mc47bd
      info(1) = 0
      if (n.lt.1) then
        info(1) = -1
        go to 1000
      endif
      if (pe(1).lt.1) then
        if (2*ne+n.gt.iwlen) then
          info(1) = -2
          go to 1000
        endif
      else
        if (ne+n.gt.iwlen) then
          info(1) = -2
          go to 1000
        endif
      endif
      if (mp.gt.0) then
        write(mp,'(/A)') 'Entry to MC47A/AD'
        write(mp,'(A,I10,A,I10,A)') 'Matrix of order',n,' with',ne, &
                                  ' entries'
        if (pe(1).lt.0)  then
          write(mp,'(A)') 'Matrix input in coordinate form'
          write(mp,'(A/(4(I8,I8,4X)))') 'Row and column indices', &
                (iw(i),iw(ne+i),i=1,ne)
        else
          write(mp,'(A)') 'Matrix input by columns'
          do 10 j=1,n
            write(mp,'(A,I4/(10I8))') 'Column',j, &
                                      (iw(i),i=pe(j),pe(j+1)-1)
   10     continue
        endif
      endif
      last   = iwlen  - n + 1
      elen   = last   - n
      nv     = elen   - n
      w      = nv     - n
      degree = w      - n
      head   = degree - n
      next   = head   - n
      len    = next   - n
      leniw = len-1
      info(6) = 0
      info(7) = 0
      if (pe(1).lt.0) then
        do 20 i=1,ne
          if (iw(i).le.iw(ne+i)) then
            if (iw(i).eq.iw(ne+i) .and. iw(i).ne.0) then
              info(7) = info(7) + 1
            else
              if (iw(i).gt.0) info(6) = info(6) + 1
            endif
            iw(i)=0
          endif
   20   continue
        ict59(1) = 0
        ict59(2) = 1
        ict59(3) = 1
        ict59(4) = mp
        ict59(5) = -1
        ict59(6) = 0
        call mc59ad(ict59,n,n,ne,iw,ne,iw(ne+1),1,dummy, &
                    n+1,pe,n+1,iw(2*ne+1),info59)
        iflag = info59(1)
        idup  = info59(3)
        iout  = info59(4)
        jout  = info59(5)
        nzout = info59(6)
      else
        idup = 0
        iout = 0
        jout = 0
        do 30 i = 1,n
          iw(ne+i) = 0
   30   continue
        do 50 j=1,n
          i1 = pe(j)
          pe(j) = i1-(iout+idup)
          i2 = pe(j+1)-1
          if (i2.lt.i1-1) then
            info(1) = -3
            go to 1000
          endif
          do 40 ii = i1,i2
            i = iw(ii)
            if (i.le.j .or. i.gt.n) then
              if (i.eq.j) info(7) = info(7) + 1
              if (i.gt.0 .and. i.lt.j) info(6) = info(6) + 1
              iout = iout + 1
            else
              if (iw(ne+i).eq.j) then
                idup = idup + 1
              else
                iw(ne+i)=j
                iw(ii-(iout+idup)) = i
              endif
            endif
   40     continue
   50   continue
        pe(n+1) = ne - (iout+idup) + 1
      endif
      if (idup.gt.0) then
        info(1) = 1
        info(4) = idup
      else
        info(4) = 0
      endif
      if (iout.gt.0 .or. jout.gt.0) then
        info(1) = 1
        info(5) = iout + jout - info(7)
      else
        info(5) = 0
      endif
      if (info(6).gt.0 .or. info(7).gt.0) info(1) = 1
      if (ne-(iout+idup).eq.0) then
        info(1) = -4
        go to 1000
      endif
      if (leniw.lt.2*(pe(n+1)-1)) then
        info(1) = -2
        go to 1000
      endif
      call mc34ad(n,iw,pe,.false.,dummy,iw(w))
      pfree = pe(n+1)
      do 60 i=1,n
        iw(len+i-1) = pe(i+1) - pe(i)
   60 continue
      call mc47bd(n,leniw,pe,pfree,iw(len),iw,iw(nv),iw(elen), &
                  iw(last),ncmpa,iw(degree),iw(head),iw(next),iw(w))
      info(2) = ncmpa
      info(3) = pfree+8*n
      if (mp.gt.0) then
        write(mp,'(/A)') 'Exit from MC47A/AD'
        write(mp,'(A/7I10)') 'INFO(1-7):',(info(i),i=1,7)
        write(mp,'(A/(8I10))') 'Parent array',(pe(i),i=1,n)
        write(mp,'(A/(8I10))') 'Permutation',(iw(elen+i-1),i=1,n)
        write(mp,'(A/(8I10))') 'Inverse permutation', &
                               (iw(last+i-1),i=1,n)
        write(mp,'(A/(8I10))') 'Degree array',(iw(nv+i-1),i=1,n)
      endif
 1000 return
      end
      subroutine mc47bd (n, iwlen, pe, pfree, len, iw, nv, elen, &
                         last, ncmpa, degree, head, next, w)
      integer n, iwlen, pe(n), pfree, len(n), iw(iwlen), nv(n), &
              elen(n), last(n), ncmpa, degree(n), head(n), next(n), &
              w(n)
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
!-----------------------------------------------------------------------
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
      wflg = 2
      mindeg = 1
      ncmpa = 0
      nel = 0
      hmod = max (1, n-1)
      dmax = 0
      mem = pfree - 1
      maxmem = mem
      do 10 i = 1, n
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
          inext = head (deg)
          if (inext .ne. 0) last (inext) = i
          next (i) = inext
          head (deg) = i
        else
          nel = nel + 1
          elen (i) = -nel
          pe (i) = 0
          w (i) = 0
        endif
   20 continue
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
        inext = next (me)
        if (inext .ne. 0) last (inext) = 0
        head (deg) = inext
        elenme = elen (me)
        elen (me) = - (nel + 1)
        nvpiv = nv (me)
        nel = nel + nvpiv
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
              ilast = last (i)
              inext = next (i)
              if (inext .ne. 0) last (inext) = ilast
              if (ilast .ne. 0) then
                next (ilast) = inext
              else
                head (degree (i)) = inext
              endif
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
                ilast = last (i)
                inext = next (i)
                if (inext .ne. 0) last (inext) = ilast
                if (ilast .ne. 0) then
                  next (ilast) = inext
                else
                  head (degree (i)) = inext
                endif
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
!=======================================================================
        do 150 pme = pme1, pme2
          i = iw (pme)
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
                we = degree (e) + wnvi
              endif
              w (e) = we
  140       continue
          endif
  150   continue
!=======================================================================
!=======================================================================
        do 180 pme = pme1, pme2
          i = iw (pme)
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
            else if (dext .eq. 0) then
              pe (e) = -me
              w (e) = 0
            endif
  160     continue
          elen (i) = pn - p1 + 1
          p3 = pn
          do 170 p = p2 + 1, p1 + len (i) - 1
            j = iw (p)
            nvj = nv (j)
            if (nvj .gt. 0) then
              deg = deg + nvj
              iw (pn) = j
              pn = pn + 1
              hash = hash + j
            endif
  170     continue
          if (deg .eq. 0) then
            pe (i) = -me
            nvi = -nv (i)
            degme = degme - nvi
            nvpiv = nvpiv + nvi
            nel = nel + nvi
            nv (i) = 0
            elen (i) = 0
          else
            degree (i) = min (degree (i), deg)
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
          if (nv (i) .lt. 0) then
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
            deg = min (degree (i) + degme - nvi, nleft - nvi)
            inext = head (deg)
            if (inext .ne. 0) last (inext) = i
            next (i) = inext
            last (i) = 0
            head (deg) = i
            mindeg = min (mindeg, deg)
            degree (i) = deg
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
!=======================================================================
      go to 30
      endif
!=======================================================================
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
! *******************************************************************
! COPYRIGHT (c) 1988 Hyprotech UK and
! Council for the Central Laboratory of the Research Councils
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
! Original date 14 June 2001
!  June 2001: threadsafe version of MC41
! 20/2/02 Cosmetic changes applied to reduce single/double differences

! 12th July 2004 Version 1.0.0. Version numbering added.

      subroutine mc71ad(n,kase,x,est,w,iw,keep)
      integer itmax
      parameter (itmax=5)
      double precision zero,one
      parameter (zero=0.0d0,one=1.0d0)
      double precision est
      integer kase,n
      double precision w(*),x(*)
      integer iw(*),keep(5)
      double precision altsgn,temp
      integer i,iter,j,jlast,jump
      integer idamax
      external idamax
      intrinsic abs,sign,nint,dble
      if (n.le.0) then
        kase = -1
        return
      end if
      if (kase.eq.0) then
        do 10 i = 1,n
          x(i) = one/dble(n)
   10   continue
        kase = 1
        jump = 1
        keep(1) = jump
        keep(2) = 0
        keep(3) = 0
        keep(4) = 0
        return
      end if
      jump  = keep(1)
      iter  = keep(2)
      j     = keep(3)
      jlast = keep(4)
      go to (100,200,300,400,500) jump
  100 continue
      if (n.eq.1) then
        w(1) = x(1)
        est = abs(w(1))
        go to 510
      end if
      do 110 i = 1,n
        x(i) = sign(one,x(i))
        iw(i) = nint(x(i))
  110 continue
      kase = 2
      jump = 2
      go to 1010
  200 continue
      j = idamax(n,x,1)
      iter = 2
  220 continue
      do 230 i = 1,n
        x(i) = zero
  230 continue
      x(j) = one
      kase = 1
      jump = 3
      go to 1010
  300 continue
      do 310 i = 1,n
        w(i) = x(i)
  310 continue
      do 320 i = 1,n
        if (nint(sign(one,x(i))).ne.iw(i)) go to 330
  320 continue
      go to 410
  330 continue
      do 340 i = 1,n
        x(i) = sign(one,x(i))
        iw(i) = nint(x(i))
  340 continue
      kase = 2
      jump = 4
      go to 1010
  400 continue
      jlast = j
      j = idamax(n,x,1)
      if ((abs(x(jlast)).ne.abs(x(j))) .and. (iter.lt.itmax)) then
        iter = iter + 1
        go to 220
      end if
  410 continue
      est = zero
      do 420 i = 1,n
        est = est + abs(w(i))
  420 continue
      altsgn = one
      do 430 i = 1,n
        x(i) = altsgn* (one+dble(i-1)/dble(n-1))
        altsgn = -altsgn
  430 continue
      kase = 1
      jump = 5
      go to 1010
  500 continue
      temp = zero
      do 520 i = 1,n
        temp = temp + abs(x(i))
  520 continue
      temp = 2.0*temp/dble(3*n)
      if (temp.gt.est) then
        do 530 i = 1,n
          w(i) = x(i)
  530   continue
        est = temp
      end if
  510 kase = 0
 1010 continue
      keep(1) = jump
      keep(2) = iter
      keep(3) = j
      keep(4) = jlast
      return
      end
! *******************************************************************
! COPYRIGHT (c) 1987 Hyprotech UK
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
! Original date 10 Feb 1993
!       Toolpack tool decs employed.
! 20/2/02 Cosmetic changes applied to reduce single/double differences

! 12th July 2004 Version 1.0.0. Version numbering added.

      subroutine mc34ad(n,irn,jcolst,yesa,a,iw)
      integer n
      logical yesa
      double precision a(*)
      integer irn(*),iw(*),jcolst(*)
      integer ckp1,i,i1,i2,ii,ipkp1,ipos,j,jstart,lenk,ndiag,newtau, &
              oldtau
      oldtau = jcolst(n+1) - 1
      do 5 i = 1,n
        iw(i) = 0
    5 continue
      ndiag = 0
      do 20 j = 1,n
        i1 = jcolst(j)
        i2 = jcolst(j+1) - 1
        iw(j) = iw(j) + i2 - i1 + 1
        do 10 ii = i1,i2
          i = irn(ii)
          if (i.ne.j) then
            iw(i) = iw(i) + 1
          else
            ndiag = ndiag + 1
          end if
   10   continue
   20 continue
      newtau = 2*oldtau - ndiag
      ipkp1 = oldtau + 1
      ckp1 = newtau + 1
      do 40 j = n,1,-1
        i1 = jcolst(j)
        i2 = ipkp1
        lenk = i2 - i1
        jstart = ckp1
        ipkp1 = i1
        i2 = i2 - 1
        do 30 ii = i2,i1,-1
          jstart = jstart - 1
          if (yesa) a(jstart) = a(ii)
          irn(jstart) = irn(ii)
   30   continue
        jcolst(j) = jstart
        ckp1 = ckp1 - iw(j)
        iw(j) = lenk
   40 continue
      do 80 j = n,1,-1
        i1 = jcolst(j)
        i2 = jcolst(j) + iw(j) - 1
        do 60 ii = i1,i2
          i = irn(ii)
          if (i.eq.j) go to 60
          jcolst(i) = jcolst(i) - 1
          ipos = jcolst(i)
          if (yesa) a(ipos) = a(ii)
          irn(ipos) = j
   60   continue
   80 continue
      jcolst(n+1) = newtau + 1
      return
      end
! *******************************************************************
! COPYRIGHT (c) 1993 Council for the Central Laboratory
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
! Original date 29 Jan 2001
! 29 January 2001. Modified from MC49 to be threadsafe.

! 12th July 2004 Version 1.0.0. Version numbering added.

      subroutine mc59ad(icntl,nc,nr,ne,irn,ljcn,jcn,la,a,lip,ip, &
                        liw,iw,info)
      integer la,lip,liw,ljcn,nc,ne,nr
      double precision a(la)
      integer icntl(10),ip(lip),info(10),irn(ne),iw(liw),jcn(ljcn)
      integer i,icntl1,icntl2,icntl3,icntl6,laa
      integer idup,iout,iup,jout,lp,mp,kne,part
      logical lcheck
      external mc59bd,mc59cd,mc59dd,mc59ed,mc59fd
      intrinsic max
      do 10 i = 1,10
         info(i) = 0
   10 continue
      icntl1 = icntl(1)
      icntl2 = icntl(2)
      icntl3 = icntl(3)
      icntl6 = icntl(6)
      lcheck = (icntl1.eq.0)
      lp = icntl(4)
      mp = icntl(5)
      if (icntl2.gt.2 .or. icntl2.lt.0) then
         info(1) = -1
         info(2) = icntl2
         if (lp.gt.0) then
            write (lp,fmt=9000) info(1)
            write (lp,fmt=9010) icntl2
         end if
         go to 70
      end if
      if (icntl6.gt.2 .or. icntl6.lt.-2) then
         info(1) = -11
         info(2) = icntl6
         if (lp.gt.0) then
            write (lp,fmt=9000) info(1)
            write (lp,fmt=9150) icntl6
         end if
         go to 70
      end if
      if (nc.lt.1) then
        info(1) = -2
        info(2) = nc
        if (lp.gt.0) then
          write (lp,fmt=9000) info(1)
          write (lp,fmt=9020) nc
        end if
        go to 70
      end if
      if (nr.lt.1) then
        info(1) = -3
        info(2) = nr
        if (lp.gt.0) then
          write (lp,fmt=9000) info(1)
          write (lp,fmt=9030) nr
        end if
        go to 70
      end if
      if (icntl6.ne.0 .and. nr.ne.nc) then
        info(1) = -3
        info(2) = nr
        if (lp.gt.0) then
          write (lp,fmt=9000) info(1)
          write (lp,fmt=9035) nc,nr
        end if
        go to 70
      end if
      if (ne.lt.1) then
        info(1) = -4
        info(2) = ne
        if (lp.gt.0) then
          write (lp,fmt=9000) info(1)
          write (lp,fmt=9040) ne
        end if
        go to 70
      end if
      if (icntl2.eq.0 .or. icntl2.eq.1) then
        if (ljcn.lt.ne) then
          info(1) = -5
          info(2) = ne
        end if
      else
        if (ljcn.lt.1) then
          info(1) = -5
          info(2) = 1
        end if
      end if
      if (info(1).eq.-5) then
         if (lp.gt.0) then
            write (lp,fmt=9000) info(1)
            write (lp,fmt=9050) ljcn,info(2)
         end if
         go to 70
      end if
      if (icntl3.eq.0) then
        if (la.lt.ne) then
          info(1) = -6
          info(2) = ne
        end if
      else
        if (la.lt.1) then
          info(1) = -6
          info(2) = 1
        end if
      end if
      if (info(1).eq.-6) then
         if (lp.gt.0) then
            write (lp,fmt=9000) info(1)
            write (lp,fmt=9060) la,info(2)
         end if
         go to 70
      end if
      if (icntl2.eq.0 .or. icntl2.eq.2) then
        if (lip.lt.nc+1) then
          info(1) = -7
          info(2) = nc+1
        end if
      else if (lip.lt.max(nr,nc)+1) then
        info(1) = -7
        info(2) = max(nr,nc)+1
      end if
      if (info(1).eq.-7) then
        if (lp.gt.0) then
          write (lp,fmt=9000) info(1)
          write (lp,fmt=9065) lip,info(2)
        end if
        go to 70
      end if
      if (liw.lt.max(nr,nc)+1) then
        info(1) = -8
        info(2) = max(nr,nc)+1
        if (lp.gt.0) then
          write (lp,fmt=9000) info(1)
          write (lp,fmt=9070) liw,info(2)
        end if
        go to 70
      end if
      laa = ne
      if (icntl3.ne.0) laa = 1
      iout = 0
      jout = 0
      idup = 0
      iup = 0
      part = 0
      if (icntl6.ne.0) part = 1
      if (icntl2.eq.0) then
        call mc59bd(lcheck,part,nc,nr,ne,irn,jcn,laa,a,ip,iw, &
                    iout,jout,kne)
        if (kne.eq.0) go to 50
        if (lcheck) call mc59ed(nc,nr,ne,irn,lip,ip,laa,a,iw,idup, &
                                kne,icntl6)
      else if (icntl2.eq.1) then
        if (icntl6.ne.0) part = -1
        call mc59bd(lcheck,part,nr,nc,ne,jcn,irn,laa,a,iw,ip, &
                    jout,iout,kne)
        if (kne.eq.0) go to 50
        if (lcheck) call mc59ed(nr,nc,ne,jcn,nr+1,iw,laa,a,ip, &
                                idup,kne,icntl6)
        call mc59cd(nc,nr,kne,irn,jcn,laa,a,ip,iw)
      else if (icntl2.eq.2) then
        if (lcheck) then
          call mc59fd(nc,nr,ne,irn,nc+1,ip,laa,a,liw,iw,idup, &
                      iout,iup,kne,icntl6,info)
          if (info(1).eq.-9) go to 40
          if (kne.eq.0) go to 50
        else
           kne = ne
        end if
        call mc59dd(nc,kne,irn,ip,laa,a)
      end if
      info(3) = idup
      info(4) = iout
      info(5) = jout
      info(6) = kne
      info(7) = iup
      if (idup.gt.0) info(1) = info(1) + 1
      if (iout.gt.0) info(1) = info(1) + 2
      if (jout.gt.0) info(1) = info(1) + 4
      if (info(1).gt.0 .and. mp.gt.0) then
        write (mp,fmt=9080) info(1)
        if (iout.gt.0) write (mp,fmt=9090) iout
        if (jout.gt.0) write (mp,fmt=9110) jout
        if (idup.gt.0) write (mp,fmt=9100) idup
        if (iup.gt.0)  write (mp,fmt=9130) iup
      end if
      go to 70
   40 info(3) = idup
      info(4) = iout
      info(7) = iup
      if (lp.gt.0) then
        write (lp,fmt=9000) info(1)
        write (lp,fmt=9140)
      end if
      go to 70
   50 info(1) = -10
      info(4) = iout
      info(5) = jout
      info(2) = iout + jout
      if (lp.gt.0) then
        write (lp,fmt=9000) info(1)
        write (lp,fmt=9120)
      end if
   70 return
 9000 format (/,' *** Error return from MC59AD *** INFO(1) = ',i3)
 9010 format (1x,'ICNTL(2) = ',i2,' is out of range')
 9020 format (1x,'NC = ',i6,' is out of range')
 9030 format (1x,'NR = ',i6,' is out of range')
 9035 format (1x,'Symmetric case. NC = ',i6,' but NR = ',i6)
 9040 format (1x,'NE = ',i10,' is out of range')
 9050 format (1x,'Increase LJCN from ',i10,' to at least ',i10)
 9060 format (1x,'Increase LA from ',i10,' to at least ',i10)
 9065 format (1x,'Increase LIP from ',i8,' to at least ',i10)
 9070 format (1x,'Increase LIW from ',i8,' to at least ',i10)
 9080 format (/,' *** Warning message from MC59AD *** INFO(1) = ',i3)
 9090 format (1x,i8,' entries in IRN supplied by the user were ', &
             /,'       out of range and were ignored by the routine')
 9100 format (1x,i8,' duplicate entries were supplied by the user')
 9110 format (1x,i8,' entries in JCN supplied by the user were ', &
             /,'       out of range and were ignored by the routine')
 9120 format (1x,'All entries out of range')
 9130 format (1x,i8,' of these entries were in the upper triangular ', &
             /,'       part of matrix')
 9140 format (1x,'Entries in IP are not monotonic increasing')
 9150 format (1x,'ICNTL(6) = ',i2,' is out of range')
      end
!***********************************************************************
      subroutine mc59bd(lcheck,part,nc,nr,ne,irn,jcn,la,a,ip,iw,iout, &
                        jout,kne)
      integer la,nc,ne,nr,iout,jout,kne,part
      logical lcheck
      double precision a(la)
      integer ip(nc+1),irn(ne),iw(nc+1),jcn(ne)
      double precision ace,acep
      integer i,ice,icep,j,jce,jcep,k,l,loc
      do 10 j = 1,nc + 1
        iw(j) = 0
   10 continue
      kne = 0
      iout = 0
      jout = 0
      if (lcheck) then
        if (la.gt.1) then
          if (part.eq.0) then
            do 20 k = 1,ne
              i = irn(k)
              j = jcn(k)
              if (i.gt.nr .or. i.lt.1) then
                iout = iout + 1
                if (j.gt.nc .or. j.lt.1)  jout = jout + 1
              else if (j.gt.nc .or. j.lt.1) then
                jout = jout + 1
              else
                kne = kne + 1
                irn(kne) = i
                jcn(kne) = j
                a(kne) = a(k)
                iw(j) = iw(j) + 1
              end if
   20       continue
          else if (part.eq.1) then
            do 21 k = 1,ne
              i = irn(k)
              j = jcn(k)
              if (i.gt.nr .or. i.lt.1) then
                iout = iout + 1
                if (j.gt.nc .or. j.lt.1)  jout = jout + 1
              else if (j.gt.nc .or. j.lt.1) then
                jout = jout + 1
              else
                kne = kne + 1
                if (i.lt.j) then
                  irn(kne) = j
                  jcn(kne) = i
                  iw(i) = iw(i) + 1
                else
                  irn(kne) = i
                  jcn(kne) = j
                  iw(j) = iw(j) + 1
                end if
                a(kne) = a(k)
              end if
   21       continue
          else if (part.eq.-1) then
            do 22 k = 1,ne
              i = irn(k)
              j = jcn(k)
              if (i.gt.nr .or. i.lt.1) then
                iout = iout + 1
                if (j.gt.nc .or. j.lt.1)  jout = jout + 1
              else if (j.gt.nc .or. j.lt.1) then
                jout = jout + 1
              else
                kne = kne + 1
                if (i.gt.j) then
                  irn(kne) = j
                  jcn(kne) = i
                  iw(i) = iw(i) + 1
                else
                  irn(kne) = i
                  jcn(kne) = j
                  iw(j) = iw(j) + 1
                end if
                a(kne) = a(k)
              end if
   22       continue
          end if
        else
          if (part.eq.0) then
            do 25 k = 1,ne
              i = irn(k)
              j = jcn(k)
              if (i.gt.nr .or. i.lt.1) then
                iout = iout + 1
                if (j.gt.nc .or. j.lt.1)  jout = jout + 1
              else if (j.gt.nc .or. j.lt.1) then
                jout = jout + 1
              else
                kne = kne + 1
                irn(kne) = i
                jcn(kne) = j
                iw(j) = iw(j) + 1
              end if
   25       continue
          else if (part.eq.1) then
            do 26 k = 1,ne
              i = irn(k)
              j = jcn(k)
              if (i.gt.nr .or. i.lt.1) then
                iout = iout + 1
                if (j.gt.nc .or. j.lt.1)  jout = jout + 1
              else if (j.gt.nc .or. j.lt.1) then
                jout = jout + 1
              else
                kne = kne + 1
                if (i.lt.j) then
                  irn(kne) = j
                  jcn(kne) = i
                  iw(i) = iw(i) + 1
                else
                  irn(kne) = i
                  jcn(kne) = j
                  iw(j) = iw(j) + 1
                end if
              end if
   26       continue
          else if (part.eq.-1) then
            do 27 k = 1,ne
              i = irn(k)
              j = jcn(k)
              if (i.gt.nr .or. i.lt.1) then
                iout = iout + 1
                if (j.gt.nc .or. j.lt.1)  jout = jout + 1
              else if (j.gt.nc .or. j.lt.1) then
                jout = jout + 1
              else
                kne = kne + 1
                if (i.gt.j) then
                  irn(kne) = j
                  jcn(kne) = i
                  iw(i) = iw(i) + 1
                else
                  irn(kne) = i
                  jcn(kne) = j
                  iw(j) = iw(j) + 1
                end if
              end if
   27       continue
          end if
        end if
        if (kne.eq.0) go to 130
      else
        kne = ne
        if (part.eq.0) then
          do 30 k = 1,ne
            j = jcn(k)
            iw(j) = iw(j) + 1
   30     continue
        else if (part.eq.1) then
          do 35 k = 1,ne
            i = irn(k)
            j = jcn(k)
            if (i.lt.j) then
               irn(k) = j
               jcn(k) = i
               iw(i) = iw(i) + 1
            else
              iw(j) = iw(j) + 1
            end if
   35     continue
        else if (part.eq.-1) then
          do 36 k = 1,ne
            i = irn(k)
            j = jcn(k)
            if (i.gt.j) then
               irn(k) = j
               jcn(k) = i
               iw(i) = iw(i) + 1
            else
              iw(j) = iw(j) + 1
            end if
   36     continue
        end if
      end if
      ip(1) = 1
      do 37 j = 2,nc + 1
        ip(j) = iw(j-1) + ip(j-1)
        iw(j-1) = ip(j-1)
   37 continue
      if (la.eq.1) then
        do 70 l = 1,nc
          do 60 k = iw(l),ip(l+1) - 1
            ice = irn(k)
            jce = jcn(k)
            do 40 j = 1,ne
              if (jce.eq.l) go to 50
              loc = iw(jce)
              jcep = jcn(loc)
              icep = irn(loc)
              iw(jce) = loc + 1
              jcn(loc) = jce
              irn(loc) = ice
              jce = jcep
              ice = icep
   40       continue
   50       jcn(k) = jce
            irn(k) = ice
   60     continue
   70   continue
      else
        do 120 l = 1,nc
          do 110 k = iw(l),ip(l+1) - 1
            ice = irn(k)
            jce = jcn(k)
            ace = a(k)
            do 90 j = 1,ne
              if (jce.eq.l) go to 100
              loc = iw(jce)
              jcep = jcn(loc)
              icep = irn(loc)
              iw(jce) = loc + 1
              jcn(loc) = jce
              irn(loc) = ice
              jce = jcep
              ice = icep
              acep = a(loc)
              a(loc) = ace
              ace = acep
   90       continue
  100       jcn(k) = jce
            irn(k) = ice
            a(k) = ace
  110     continue
  120   continue
      end if
  130 continue
      return
      end
!**********************************************************
      subroutine mc59cd(nc,nr,ne,irn,jcn,la,a,ip,iw)
      integer la,nc,ne,nr
      double precision a(la)
      integer ip(nc+1),irn(ne),iw(nr+1),jcn(ne)
      double precision ace,acep
      integer i,ice,icep,j,j1,j2,k,l,loc,locp
      do 10 j = 1,nc
        ip(j) = 0
   10 continue
      if (la.gt.1) then
        do 20 k = 1,ne
          i = jcn(k)
          ip(i) = ip(i) + 1
          irn(k) = jcn(k)
   20   continue
        ip(nc+1) = ne + 1
        ip(1) = ip(1) + 1
        do 30 j = 2,nc
          ip(j) = ip(j) + ip(j-1)
   30   continue
        do 50 i = nr,1,-1
          j1 = iw(i)
          j2 = iw(i+1) - 1
          do 40 j = j1,j2
            k = irn(j)
            l = ip(k) - 1
            jcn(j) = l
            irn(j) = i
            ip(k) = l
   40     continue
   50   continue
        ip(nc+1) = ne + 1
        do 70 j = 1,ne
          loc = jcn(j)
          if (loc.eq.0) go to 70
          ice = irn(j)
          ace = a(j)
          jcn(j) = 0
          do 60 k = 1,ne
            locp = jcn(loc)
            icep = irn(loc)
            acep = a(loc)
            jcn(loc) = 0
            irn(loc) = ice
            a(loc) = ace
            if (locp.eq.0) go to 70
            ice = icep
            ace = acep
            loc = locp
   60     continue
   70   continue
      else
        do 90 k = 1,ne
          i = jcn(k)
          ip(i) = ip(i) + 1
   90   continue
        ip(nc+1) = ne + 1
        ip(1) = ip(1) + 1
        do 100 j = 2,nc
          ip(j) = ip(j) + ip(j-1)
  100   continue
        do 120 i = nr,1,-1
          j1 = iw(i)
          j2 = iw(i+1) - 1
          do 110 j = j1,j2
            k = jcn(j)
            l = ip(k) - 1
            irn(l) = i
            ip(k) = l
  110     continue
  120   continue
      end if
      return
      end
!**********************************************************
      subroutine mc59dd(nc,ne,irn,ip,la,a)
      integer la,nc,ne
      double precision a(la)
      integer irn(ne),ip(nc)
      double precision ace
      integer ice,ik,j,jj,k,kdummy,klo,kmax,kor
      intrinsic abs
      if (la.gt.1) then
        kmax = ne
        do 50 jj = 1,nc
          j = nc + 1 - jj
          klo = ip(j) + 1
          if (klo.gt.kmax) go to 40
          kor = kmax
          do 30 kdummy = klo,kmax
            ace = a(kor-1)
            ice = irn(kor-1)
            do 10 k = kor,kmax
              ik = irn(k)
              if (abs(ice).le.abs(ik)) go to 20
              irn(k-1) = ik
              a(k-1) = a(k)
   10       continue
            k = kmax + 1
   20       irn(k-1) = ice
            a(k-1) = ace
            kor = kor - 1
   30     continue
   40     kmax = klo - 2
   50   continue
      else
        kmax = ne
        do 150 jj = 1,nc
          j = nc + 1 - jj
          klo = ip(j) + 1
          if (klo.gt.kmax) go to 140
          kor = kmax
          do 130 kdummy = klo,kmax
            ice = irn(kor-1)
            do 110 k = kor,kmax
              ik = irn(k)
              if (abs(ice).le.abs(ik)) go to 120
              irn(k-1) = ik
  110       continue
            k = kmax + 1
  120       irn(k-1) = ice
            kor = kor - 1
  130     continue
  140     kmax = klo - 2
  150   continue
      end if
      end
!***********************************************************************
      subroutine mc59ed(nc,nr,ne,irn,lip,ip,la,a,iw,idup,kne,icntl6)
      integer icntl6,idup,kne,lip,la,nc,nr,ne
      double precision a(la)
      integer irn(ne),ip(lip),iw(nr)
      integer i,j,k,kstart,kstop,nzj
      idup = 0
      kne = 0
      do 10 i = 1,nr
        iw(i) = 0
   10 continue
      kstart = ip(1)
      if (la.gt.1) then
        nzj = 0
        do 30 j = 1,nc
          kstop = ip(j+1)
          ip(j+1) = ip(j)
          do 20 k = kstart,kstop - 1
            i = irn(k)
            if (iw(i).le.nzj) then
              kne = kne + 1
              irn(kne) = i
              a(kne) = a(k)
              ip(j+1) = ip(j+1) + 1
              iw(i) = kne
            else
              idup = idup + 1
              if (icntl6.ge.0) a(iw(i)) = a(iw(i)) + a(k)
            end if
   20     continue
          kstart = kstop
          nzj = kne
   30   continue
      else
        do 50 j = 1,nc
          kstop = ip(j+1)
          ip(j+1) = ip(j)
          do 40 k = kstart,kstop - 1
            i = irn(k)
            if (iw(i).lt.j) then
              kne = kne + 1
              irn(kne) = i
              ip(j+1) = ip(j+1) + 1
              iw(i) = j
            else
              idup = idup + 1
            end if
   40     continue
          kstart = kstop
   50   continue
      end if
      return
      end
!***********************************************************************
      subroutine mc59fd(nc,nr,ne,irn,lip,ip,la,a,liw,iw,idup,iout, &
                        iup,kne,icntl6,info)
      integer icntl6,idup,iout,iup,kne,la,lip,liw,nc,nr,ne
      double precision a(la)
      integer irn(ne),ip(lip),iw(liw),info(2)
      integer i,j,k,kstart,kstop,nzj,lower
      idup = 0
      iout = 0
      iup = 0
      kne = 0
      do 10 i = 1,nr
        iw(i) = 0
   10 continue
      kstart = ip(1)
      lower = 1
      if (la.gt.1) then
        nzj = 0
        do 30 j = 1,nc
          if (icntl6.ne.0) lower = j
          kstop = ip(j+1)
          if (kstart.gt.kstop) then
            info(1) = -9
            info(2) = j
            return
          end if
          ip(j+1) = ip(j)
          do 20 k = kstart,kstop - 1
            i = irn(k)
            if (i.gt.nr .or. i.lt.lower) then
              iout = iout + 1
              if (icntl6.ne.0 .and. i.lt.j) iup = iup + 1
            else if (iw(i).le.nzj) then
              kne = kne + 1
              irn(kne) = i
              a(kne) = a(k)
              ip(j+1) = ip(j+1) + 1
              iw(i) = kne
            else
              idup = idup + 1
              if (icntl6.ge.0) a(iw(i)) = a(iw(i)) + a(k)
            end if
   20     continue
          kstart = kstop
          nzj = kne
   30   continue
      else
        do 50 j = 1,nc
          if (icntl6.ne.0) lower = j
          kstop = ip(j+1)
          if (kstart.gt.kstop) then
            info(1) = -9
            info(2) = j
            return
          end if
          ip(j+1) = ip(j)
          do  40 k = kstart,kstop - 1
            i = irn(k)
            if (i.gt.nr .or. i.lt.lower) then
              iout = iout + 1
              if (icntl6.ne.0 .and. i.gt.1) iup = iup + 1
            else if (iw(i).lt.j) then
              kne = kne + 1
              irn(kne) = i
              ip(j+1) = ip(j+1) + 1
              iw(i) = j
            else
              idup = idup + 1
            end if
   40     continue
          kstart = kstop
   50   continue
      end if
      return
      end
! COPYRIGHT (c) 1982 AEA Technology
!######DATE 20 September 2001
!  September 2001: threadsafe version of MA27
!  19/3/03. Array ICNTL in MA27GD made assumed size.

      subroutine ma27id(icntl,cntl)
      integer icntl(30)
      double precision cntl(5)

      integer ifrlvl
      parameter ( ifrlvl=5 )

! Stream number for error messages
      icntl(1) = 6
! Stream number for diagnostic messages
      icntl(2) = 6
! Control the level of diagnostic printing.
!   0 no printing
!   1 printing of scalar parameters and first parts of arrays.
!   2 printing of scalar parameters and whole of arrays.
      icntl(3) = 0
! The largest integer such that all integers I in the range
! -ICNTL(4).LE.I.LE.ICNTL(4) can be handled by the shortest integer
! type in use.
      icntl(4) = 2139062143
! Minimum number of eliminations in a step that is automatically
! accepted. if two adjacent steps can be combined and each has less
! eliminations then they are combined.
      icntl(5) = 1
! Control whether direct or indirect access is used by MA27C/CD.
! Indirect access is employed in forward and back substitution
! respectively if the size of a block is less than
! ICNTL(IFRLVL+MIN(10,NPIV)) and ICNTL(IFRLVL+10+MIN(10,NPIV))
! respectively, where NPIV is the number of pivots in the block.
      icntl(ifrlvl+1)  = 32639
      icntl(ifrlvl+2)  = 32639
      icntl(ifrlvl+3)  = 32639
      icntl(ifrlvl+4)  = 32639
      icntl(ifrlvl+5)  = 14
      icntl(ifrlvl+6)  = 9
      icntl(ifrlvl+7)  = 8
      icntl(ifrlvl+8)  = 8
      icntl(ifrlvl+9)  = 9
      icntl(ifrlvl+10) = 10
      icntl(ifrlvl+11) = 32639
      icntl(ifrlvl+12) = 32639
      icntl(ifrlvl+13) = 32639
      icntl(ifrlvl+14) = 32689
      icntl(ifrlvl+15) = 24
      icntl(ifrlvl+16) = 11
      icntl(ifrlvl+17) = 9
      icntl(ifrlvl+18) = 8
      icntl(ifrlvl+19) = 9
      icntl(ifrlvl+20) = 10
! Not used
      icntl(26) = 0
      icntl(27) = 0
      icntl(28) = 0
      icntl(29) = 0
      icntl(30) = 0

! Control threshold pivoting.
      cntl(1) = 0.1d0
! If a column of the reduced matrix has relative density greater than
! CNTL(2), it is forced into the root. All such columns are taken to
! have sparsity pattern equal to their merged patterns, so the fill
! and operation counts may be overestimated.
      cntl(2) = 1.0d0
! An entry with absolute value less than CNTL(3) is never accepted as
! a 1x1 pivot or as the off-diagonal of a 2x2 pivot.
      cntl(3) = 0.0d0
! Not used
      cntl(4) = 0.0
      cntl(5) = 0.0

      return
      end

      subroutine ma27ad(n,nz,irn,icn,iw,liw,ikeep,iw1,nsteps,iflag, &
                       icntl,cntl,info,ops)
! THIS SUBROUTINE COMPUTES A MINIMUM DEGREE ORDERING OR ACCEPTS A GIVEN
!     ORDERING. IT COMPUTES ASSOCIATED ASSEMBLY AND ELIMINATION
!     INFORMATION FOR MA27B/BD.
! N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
! NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
!     ALTERED.
! IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW NUMBERS OF THE
!     NON-ZEROS. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
!     TO IW (SEE DESCRIPTION OF IW).
! ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN NUMBERS OF THE
!     NON-ZEROS. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
!     TO IW (SEE DESCRIPTION OF IW).
! IW NEED NOT BE SET ON INPUT. IT IS USED FOR WORKSPACE.
!     IRN(1) MAY BE EQUIVALENCED TO IW(1) AND ICN(1) MAY BE
!     EQUIVALENCED TO IW(K), WHERE K.GT.NZ.
! LIW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST 2*NZ+3*N
!     FOR THE IFLAG=0 ENTRY AND AT LEAST NZ+3*N FOR THE IFLAG=1
!     ENTRY. IT IS NOT ALTERED.
! IKEEP NEED NOT BE SET UNLESS AN ORDERING IS GIVEN, IN WHICH CASE
!     IKEEP(I,1) MUST BE SET TO THE POSITION OF VARIABLE I IN THE
!     ORDER. ON OUTPUT IKEEP CONTAINS INFORMATION NEEDED BY MA27B/BD.
!     IKEEP(I,1) HOLDS THE POSITION OF VARIABLE I IN THE PIVOT ORDER.
!     IKEEP(I,2), IKEEP(I,3) HOLD THE NUMBER OF ELIMINATIONS, ASSEMBLIES
!     AT MAJOR STEP I, I=1,2,...,NSTEPS. NOTE THAT WHEN AN ORDER IS
!     GIVEN IT MAY BE REPLACED BY ANOTHER ORDER THAT GIVES IDENTICAL
!     NUMERICAL RESULTS.
! IW1 IS USED FOR WORKSPACE.
! NSTEPS NEED NOT BE SET. ON OUTPUT IT CONTAINS THE NUMBER OF MAJOR
!     STEPS NEEDED FOR A DEFINITE MATRIX AND MUST BE PASSED UNCHANGED
!     TO MA27B/BD.
! IFLAG MUST SET TO ZERO IF THE USER WANTS THE PIVOT ORDER CHOSEN
!     AUTOMATICALLY AND TO ONE IF HE WANTS TO SPECIFY IT IN IKEEP.
! ICNTL is an INTEGER array of length 30 containing user options
!     with integer values (defaults set in MA27I/ID)
!   ICNTL(1) (LP) MUST BE SET TO THE STREAM NUMBER FOR ERROR MESSAGES.
!     ERROR MESSAGES ARE SUPPRESSED IF ICNTL(1) IS NOT POSITIVE.
!     IT IS NOT ALTERED.
!   ICNTL(2) (MP) MUST BE SET TO THE STREAM NUMBER FOR DIAGNOSTIC
!     MESSAGES.  DIAGNOSTIC MESSAGES ARE SUPPRESSED IF ICNTL(2) IS NOT
!     POSITIVE.  IT IS NOT ALTERED.
!   ICNTL(3) (LDIAG) CONTROLS THE LEVEL OF DIAGNOSTIC PRINTING.
!     0 NO PRINTING
!     1 PRINTING OF SCALAR PARAMETERS AND FIRST PARTS OF ARRAYS.
!     2 PRINTING OF SCALAR PARAMETERS AND WHOLE OF ARRAYS.
!   ICNTL(4) (IOVFLO) IS THE LARGEST INTEGER SUCH THAT ALL INTEGERS
!     I IN THE RANGE -IOVFLO.LE.I.LE.IOVFLO CAN BE HANDLED BY THE
!     SHORTEST INTEGER TYPE IN USE.
!   ICNT(5) (NEMIN) MUST BE SET TO THE MINIMUM NUMBER OF ELIMINATIONS
!     IN A STEP THAT IS AUTOMATICALLY ACCEPTED. IF TWO ADJACENT STEPS
!     CAN BE COMBINED AND EACH HAS LESS ELIMINATIONS THEN THEY ARE
!     COMBINED.
!   ICNTL(IFRLVL+I) I=1,20, (IFRLVL) MUST BE SET TO CONTROL WHETHER
!     DIRECT OR INDIRECT ACCESS IS USED BY MA27C/CD.  INDIRECT ACCESS
!     IS EMPLOYED IN FORWARD AND BACK SUBSTITUTION RESPECTIVELY IF THE
!     SIZE OF A BLOCK IS LESS THAN ICNTL(IFRLVL+(MIN(10,NPIV)) AND
!     ICNTL(IFRLVL+10+MIN(10,NPIV)) RESPECTIVELY, WHERE NPIV IS THE
!     NUMBER OF PIVOTS IN THE BLOCK.
!   ICNTL(I) I=26,30 are not used.
! CNTL is an DOUBLE PRECISION array of length 5 containing user options
!     with real values (defaults set in MA27I/ID)
!   CNTL(1) (U) IS USED TO CONTROL THRESHOLD PIVOTING. IT IS NOT
!     ALTERED.
!   CNTL(2) (FRATIO) has default value 1.0.  If a column of the
!      reduced matrix has relative density greater than CNTL(2), it
!      is forced into the root. All such columns are taken to have
!      sparsity pattern equal to their merged patterns, so the fill
!      and operation counts may be overestimated.
!   CNTL(3) (PIVTOL) has default value 0.0. An entry with absolute
!      value less than CNTL(3) is never accepted as a 1x1 pivot or
!      as the off-diagonal of a 2x2 pivot.
!   CNTL(I) I=4,5 are not used.
! INFO is an INTEGER array of length 20 which is used to return
!     information to the user.
!   INFO(1) (IFLAG) is an error return code, zero for success, greater
!      than zero for a warning and less than zero for errors, see
!      INFO(2).
!   INFO(2) (IERROR) HOLDS ADDITIONAL INFORMATION IN THE EVENT OF ERRORS.
!     IF INFO(1)=-3 INFO(2) HOLDS A LENGTH THAT MAY SUFFICE FOR IW.
!     IF INFO(1)=-4 INFO(2) HOLDS A LENGTH THAT MAY SUFFICE FOR A.
!     IF INFO(1)=-5 INFO(2) IS SET TO THE PIVOT STEP AT WHICH SINGULARITY
!                 WAS DETECTED.
!     IF INFO(1)=-6 INFO(2) IS SET TO THE PIVOT STEP AT WHICH A CHANGE OF
!                 PIVOT SIGN WAS FOUND.
!     IF INFO(1)= 1 INFO(2) HOLDS THE NUMBER OF FAULTY ENTRIES.
!     IF INFO(1)= 2 INFO(2) IS SET TO THE NUMBER OF SIGNS CHANGES IN
!                 THE PIVOTS.
!     IF INFO(1)=3 INFO(2) IS SET TO THE RANK OF THE MATRIX.
!   INFO(3) and INFO(4) (NRLTOT and NIRTOT) REAL AND INTEGER STRORAGE
!     RESPECTIVELY REQUIRED FOR THE FACTORIZATION IF NO COMPRESSES ARE
!     ALLOWED.
!   INFO(5) and INFO(6) (NRLNEC and NIRNEC) REAL AND INTEGER STORAGE
!     RESPECTIVELY REQUIRED FOR THE FACTORIZATION IF COMPRESSES ARE
!     ALLOWED AND THE MATRIX IS DEFINITE.
!   INFO(7) and INFO(8) (NRLADU and NIRADU) REAL AND INTEGER STORAGE
!     RESPECTIVELY FOR THE MATRIX FACTORS AS CALCULATED BY MA27A/AD
!     FOR THE DEFINITE CASE.
!   INFO(9) and INFO(10) (NRLBDU and NIRBDU) REAL AND INTEGER STORAGE
!     RESPECTIVELY FOR THE MATRIX FACTORS AS FOUND  BY MA27B/BD.
!   INFO(11) (NCMPA) ACCUMULATES THE NUMBER OF TIMES THE ARRAY IW IS
!     COMPRESSED BY MA27A/AD.
!   INFO(12) and INFO(13) (NCMPBR and NCMPBI) ACCUMULATE THE NUMBER
!     OF COMPRESSES OF THE REALS AND INTEGERS PERFORMED BY MA27B/BD.
!   INFO(14) (NTWO) IS USED BY MA27B/BD TO RECORD THE NUMBER OF 2*2
!     PIVOTS USED.
!   INFO(15) (NEIG) IS USED BY ME27B/BD TO RECORD THE NUMBER OF
!     NEGATIVE EIGENVALUES OF A.
!   INFO(16) to INFO(20) are not used.
! OPS ACCUMULATES THE NO. OF MULTIPLY/ADD PAIRS NEEDED TO CREATE THE
!     TRIANGULAR FACTORIZATION, IN THE DEFINITE CASE.
!
!     .. Scalar Arguments ..
      integer iflag,liw,n,nsteps,nz
!     ..
!     .. Array Arguments ..
      double precision cntl(5),ops
      integer icntl(30),info(20)
      integer icn(*),ikeep(n,3),irn(*),iw(liw),iw1(n,2)
!     ..
!     .. Local Scalars ..
      integer i,iwfr,k,l1,l2,lliw
!     ..
!     .. External Subroutines ..
      external ma27gd,ma27hd,ma27jd,ma27kd,ma27ld,ma27md,ma27ud
!     ..
!     .. Intrinsic Functions ..
      intrinsic min
!     ..
!     .. Executable Statements ..
      do 5 i = 1,15
        info(i) = 0
    5 continue

      if (icntl(3).le.0 .or. icntl(2).le.0) go to 40
! PRINT INPUT VARIABLES.
      write (icntl(2),fmt=10) n,nz,liw,iflag

   10 format(/,/,' ENTERING MA27AD WITH      N     NZ      LIW  IFLAG', &
             /,21x,i7,i7,i9,i7)

      nsteps = 0
      k = min(8,nz)
      if (icntl(3).gt.1) k = nz
      if (k.gt.0) write (icntl(2),fmt=20) (irn(i),icn(i),i=1,k)

   20 format (' MATRIX NON-ZEROS',/,4 (i9,i6),/, &
             (i9,i6,i9,i6,i9,i6,i9,i6))

      k = min(10,n)
      if (icntl(3).gt.1) k = n
      if (iflag.eq.1 .and. k.gt.0) then
        write (icntl(2),fmt=30) (ikeep(i,1),i=1,k)
      end if

   30 format (' IKEEP(.,1)=',10i6,/, (12x,10i6))

   40 if (n.lt.1 .or. n.gt.icntl(4)) go to 70
      if (nz.lt.0) go to 100
      lliw = liw - 2*n
      l1 = lliw + 1
      l2 = l1 + n
      if (iflag.eq.1) go to 50
      if (liw.lt.2*nz+3*n+1) go to 130
! SORT
      call ma27gd(n,nz,irn,icn,iw,lliw,iw1,iw1(1,2),iw(l1),iwfr, &
                 icntl,info)
! ANALYZE USING MINIMUM DEGREE ORDERING
      call ma27hd(n,iw1,iw,lliw,iwfr,iw(l1),iw(l2),ikeep(1,2), &
                 ikeep(1,3),ikeep,icntl(4),info(11),cntl(2))
      go to 60
! SORT USING GIVEN ORDER
   50 if (liw.lt.nz+3*n+1) go to 120
      call ma27jd(n,nz,irn,icn,ikeep,iw,lliw,iw1,iw1(1,2),iw(l1),iwfr, &
                 icntl,info)
! ANALYZE USING GIVEN ORDER
      call ma27kd(n,iw1,iw,lliw,iwfr,ikeep,ikeep(1,2),iw(l1),iw(l2), &
                 info(11))
! PERFORM DEPTH-FIRST SEARCH OF ASSEMBLY TREE
   60 call ma27ld(n,iw1,iw(l1),ikeep,ikeep(1,2),ikeep(1,3),iw(l2), &
                 nsteps,icntl(5))
! EVALUATE STORAGE AND OPERATION COUNT REQUIRED BY MA27B/BD IN THE
!     DEFINITE CASE.
! SET IW(1) SO THAT ARRAYS IW AND IRN CAN BE TESTED FOR EQUIVALENCE
!     IN MA27M/MD.
      if(nz.ge.1) iw(1) = irn(1) + 1
      call ma27md(n,nz,irn,icn,ikeep,ikeep(1,3),ikeep(1,2),iw(l2), &
                 nsteps,iw1,iw1(1,2),iw,info,ops)
      go to 160

   70 info(1) = -1
      if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)

   80 format (' **** ERROR RETURN FROM MA27AD **** INFO(1)=',i3)

      if (icntl(1).gt.0) write (icntl(1),fmt=90) n

   90 format (' VALUE OF N OUT OF RANGE ... =',i10)

      go to 160

  100 info(1) = -2
      if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)
      if (icntl(1).gt.0) write (icntl(1),fmt=110) nz

  110 format (' VALUE OF NZ OUT OF RANGE .. =',i10)

      go to 160

  120 info(2) = nz + 3*n + 1
      go to 140

  130 info(2) = 2*nz + 3*n + 1
  140 info(1) = -3
      if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)
      if (icntl(1).gt.0) write (icntl(1),fmt=150) liw,info(2)

  150 format (' LIW TOO SMALL, MUST BE INCREASED FROM',i10, &
             ' TO AT LEAST',i10)

  160 if (icntl(3).le.0 .or. icntl(2).le.0 .or. info(1).lt.0) go to 200
! PRINT PARAMETER VALUES ON EXIT.
      write (icntl(2),fmt=170) nsteps,info(1),ops,info(2),info(3), &
        info(4),info(5),info(6),info(7),info(8),info(11)

  170 format (/,' LEAVING MA27AD WITH NSTEPS  INFO(1)    OPS IERROR', &
                ' NRLTOT NIRTOT', &
              /,20x,2i7,f7.0,3i7, &
              /,20x,' NRLNEC NIRNEC NRLADU NIRADU  NCMPA', &
              /,20x,6i7)

      k = min(9,n)
      if (icntl(3).gt.1) k = n
      if (k.gt.0) write (icntl(2),fmt=30) (ikeep(i,1),i=1,k)
      k = min(k,nsteps)
      if (k.gt.0) write (icntl(2),fmt=180) (ikeep(i,2),i=1,k)

  180 format (' IKEEP(.,2)=',10i6,/, (12x,10i6))

      if (k.gt.0) write (icntl(2),fmt=190) (ikeep(i,3),i=1,k)

  190 format (' IKEEP(.,3)=',10i6,/, (12x,10i6))

  200 continue

      return
      end

      subroutine ma27bd(n,nz,irn,icn,a,la,iw,liw,ikeep,nsteps,maxfrt, &
                       iw1,icntl,cntl,info)
! THIS SUBROUTINE COMPUTES THE FACTORISATION OF THE MATRIX INPUT IN
!     A,IRN,ICN USING INFORMATION (IN IKEEP) FROM MA27A/AD.
! N MUST BE SET TO THE ORDER OF THE MATRIX.  IT IS NOT ALTERED.
! NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT.  IT IS NOT
!     ALTERED.
! IRN,ICN,A.  ENTRY K (K=1,NZ) OF IRN,ICN MUST BE SET TO THE ROW
!     AND COLUMN INDEX RESPECTIVELY OF THE NON-ZERO IN A.
!     IRN AND ICN ARE UNALTERED BY MA27B/BD.
!     ON EXIT, ENTRIES 1 TO NRLBDU OF A HOLD REAL INFORMATION
!     ON THE FACTORS AND SHOULD BE PASSED UNCHANGED TO MA27C/CD.
! LA LENGTH OF ARRAY A.  AN INDICATION OF A SUITABLE VALUE,
!     SUFFICIENT FOR FACTORIZATION OF A DEFINITE MATRIX, WILL
!     HAVE BEEN PROVIDED IN NRLNEC AND NRLTOT BY MA27A/AD.
!     IT IS NOT ALTERED BY MA27B/BD.
! IW NEED NOT BE SET ON ENTRY.  USED AS A WORK ARRAY BY MA27B/BD.
!     ON EXIT, ENTRIES 1 TO NIRBDU HOLD INTEGER INFORMATION ON THE
!     FACTORS AND SHOULD BE PASSED UNCHANGED TO MA27C/CD.
! LIW LENGTH OF ARRAY IW.  AN INDICATION OF A SUITABLE VALUE WILL
!     HAVE BEEN PROVIDED IN NIRNEC AND NIRTOT BY MA27A/AD.
!     IT IS NOT ALTERED BY MA27B/BD.
! IKEEP MUST BE UNCHANGED SINCE THE CALL TO MA27A/AD.
!     IT IS NOT ALTERED BY MA27B/BD.
! NSTEPS MUST BE UNCHANGED SINCE THE CALL TO MA27A/AD.
!     IT IS NOT ALTERED BY MA27B/BD.
! MAXFRT NEED NOT BE SET AND MUST BE PASSED UNCHANGED TO MA27C/CD.
!     IT IS THE MAXIMUM SIZE OF THE FRONTAL MATRIX GENERATED DURING
!     THE DECOMPOSITION.
! IW1 USED AS WORKSPACE BY MA27B/BD.
! ICNTL is an INTEGER array of length 30, see MA27A/AD.
! CNTL is a DOUBLE PRECISION array of length 5, see MA27A/AD.
! INFO is an INTEGER array of length 20, see MA27A/AD.
!
!     .. Scalar Arguments ..
      integer la,liw,maxfrt,n,nsteps,nz
!     ..
!     .. Array Arguments ..
      double precision a(la),cntl(5)
      integer icn(*),ikeep(n,3),irn(*),iw(liw),iw1(n)
      integer icntl(30),info(20)
!     ..
!     .. Local Scalars ..
      integer i,iapos,iblk,ipos,irows,j1,j2,jj,k,kblk,kz,len,ncols, &
              nrows,nz1
!     ..
!     .. External Subroutines ..
      external ma27nd,ma27od
!     ..
!     .. Intrinsic Functions ..
      intrinsic abs,min
!     ..
!     .. Executable Statements ..
      info(1) = 0

      if (icntl(3).le.0 .or. icntl(2).le.0) go to 60
! PRINT INPUT PARAMETERS.
      write (icntl(2),fmt=10) n,nz,la,liw,nsteps,cntl(1)

   10 format (/,/, &
       ' ENTERING MA27BD WITH      N     NZ       LA      LIW', &
             ' NSTEPS      U',/,21x,i7,i7,i9,i9,i7,1pd10.2)

      kz = min(6,nz)
      if (icntl(3).gt.1) kz = nz
      if (nz.gt.0) write (icntl(2),fmt=20) (a(k),irn(k),icn(k),k=1,kz)

   20 format (' MATRIX NON-ZEROS',/,1x,2 (1p,d16.3,2i6),/, &
             (1x,1p,d16.3,2i6,1p,d16.3,2i6))

      k = min(9,n)
      if (icntl(3).gt.1) k = n
      if (k.gt.0) write (icntl(2),fmt=30) (ikeep(i,1),i=1,k)

   30 format (' IKEEP(.,1)=',10i6,/, (12x,10i6))

      k = min(k,nsteps)
      if (k.gt.0) write (icntl(2),fmt=40) (ikeep(i,2),i=1,k)

   40 format (' IKEEP(.,2)=',10i6,/, (12x,10i6))

      if (k.gt.0) write (icntl(2),fmt=50) (ikeep(i,3),i=1,k)

   50 format (' IKEEP(.,3)=',10i6,/, (12x,10i6))

   60 if (n.lt.1 .or. n.gt.icntl(4)) go to 70
      if (nz.lt.0) go to 100
      if (liw.lt.nz) go to 120
      if (la.lt.nz+n) go to 150
      if (nsteps.lt.1 .or. nsteps.gt.n) go to 175
! SORT
      call ma27nd(n,nz,nz1,a,la,irn,icn,iw,liw,ikeep,iw1,icntl,info)
      if (info(1).eq.-3) go to 130
      if (info(1).eq.-4) go to 160
! FACTORIZE
      call ma27od(n,nz1,a,la,iw,liw,ikeep,ikeep(1,3),nsteps,maxfrt, &
                 ikeep(1,2),iw1,icntl,cntl,info)
      if (info(1).eq.-3) go to 130
      if (info(1).eq.-4) go to 160
      if (info(1).eq.-5) go to 180
      if (info(1).eq.-6) go to 200
! **** WARNING MESSAGE ****
      if (info(1).eq.3 .and. icntl(2).gt.0) then
        write (icntl(2),fmt=65) info(1),info(2)
      end if

   65 format (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD', &
              '  *** INFO(1) =',i2, &
              /,5x,'MATRIX IS SINGULAR. RANK=',i5)

      go to 220
! **** ERROR RETURNS ****
   70 info(1) = -1
      if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)

   80 format (' **** ERROR RETURN FROM MA27BD **** INFO(1)=',i3)

      if (icntl(1).gt.0) write (icntl(1),fmt=90) n

   90 format (' VALUE OF N OUT OF RANGE ... =',i10)

      go to 220

  100 info(1) = -2
      if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)
      if (icntl(1).gt.0) write (icntl(1),fmt=110) nz

  110 format (' VALUE OF NZ OUT OF RANGE .. =',i10)

      go to 220

  120 info(1) = -3
      info(2) = nz
  130 if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)
      if (icntl(1).gt.0) write (icntl(1),fmt=140) liw,info(2)

  140 format (' LIW TOO SMALL, MUST BE INCREASED FROM',i10,' TO', &
             ' AT LEAST',i10)

      go to 220

  150 info(1) = -4
      info(2) = nz + n
  160 if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)
      if (icntl(1).gt.0) write (icntl(1),fmt=170) la,info(2)

  170 format (' LA TOO SMALL, MUST BE INCREASED FROM ',i10,' TO', &
             ' AT LEAST',i10)

      go to 220

  175 info(1) = -7
      if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)
      if (icntl(1).gt.0) then
        write (icntl(1),fmt='(A)') ' NSTEPS is out of range'
      end if
      go to 220

  180 if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)
      if (icntl(1).gt.0) write (icntl(1),fmt=190) info(2)

  190 format (' ZERO PIVOT AT STAGE',i10, &
              ' WHEN INPUT MATRIX DECLARED DEFINITE')

      go to 220

  200 if (icntl(1).gt.0) write (icntl(1),fmt=80) info(1)
      if (icntl(1).gt.0) write (icntl(1),fmt=210)

  210 format (' CHANGE IN SIGN OF PIVOT ENCOUNTERED', &
              ' WHEN FACTORING ALLEGEDLY DEFINITE MATRIX')

  220 if (icntl(3).le.0 .or. icntl(2).le.0 .or. info(1).lt.0) go to 310
! PRINT OUTPUT PARAMETERS.
      write (icntl(2),fmt=230) maxfrt,info(1),info(9),info(10),info(12), &
        info(13),info(14),info(2)

  230 format (/,' LEAVING MA27BD WITH', &
              /,10x,'  MAXFRT  INFO(1) NRLBDU NIRBDU NCMPBR', &
               ' NCMPBI   NTWO IERROR', &
              /,11x,8i7)

      if (info(1).lt.0) go to 310
! PRINT OUT MATRIX FACTORS FROM MA27B/BD.
      kblk = abs(iw(1)+0)
      if (kblk.eq.0) go to 310
      if (icntl(3).eq.1) kblk = 1
      ipos = 2
      iapos = 1
      do 300 iblk = 1,kblk
        ncols = iw(ipos)
        nrows = iw(ipos+1)
        j1 = ipos + 2
        if (ncols.gt.0) go to 240
        ncols = -ncols
        nrows = 1
        j1 = j1 - 1
  240   write (icntl(2),fmt=250) iblk,nrows,ncols

  250   format (' BLOCK PIVOT =',i8,' NROWS =',i8,' NCOLS =',i8)

        j2 = j1 + ncols - 1
        ipos = j2 + 1
        write (icntl(2),fmt=260) (iw(jj),jj=j1,j2)

  260   format (' COLUMN INDICES =',10i6,/, (17x,10i6))

        write (icntl(2),fmt=270)

  270   format (' REAL ENTRIES .. EACH ROW STARTS ON A NEW LINE')

        len = ncols
        do 290 irows = 1,nrows
          j1 = iapos
          j2 = iapos + len - 1
          write (icntl(2),fmt=280) (a(jj),jj=j1,j2)

  280     format (1p,5d13.3)

          len = len - 1
          iapos = j2 + 1
  290   continue
  300 continue
  310 return
      end

      subroutine ma27cd(n,a,la,iw,liw,w,maxfrt,rhs,iw1,nsteps, &
       icntl,info)
! THIS SUBROUTINE USES THE FACTORISATION OF THE MATRIX IN A,IW TO
!     SOLVE A SYSTEM OF EQUATIONS.
! N MUST BE SET TO THE ORDER OF THE MATRIX.  IT IS NOT ALTERED.
! A,IW HOLD INFORMATION ON THE FACTORS AND MUST BE UNCHANGED SINCE
!     THE CALL TO MA27B/BD. THEY ARE NOT ALTERED BY MA27C/CDD.
! LA,LIW MUST BE SET TO THE LENGTHS OF A,IW RESPECTIVELY.  THEY
!     ARE NOT ALTERED.
! W USED AS A WORK ARRAY.
! MAXFRT IS THE LENGTH OF W AND MUST BE PASSED UNCHANGED FROM THE
!     CALL TO MA27B/BD.  IT IS NOT ALTERED.
! RHS MUST BE SET TO THE RIGHT HAND SIDE FOR THE EQUATIONS BEING
!     SOLVED.  ON EXIT, THIS ARRAY WILL HOLD THE SOLUTION.
! IW1 USED AS A WORK ARRAY.
! NSTEPS IS THE LENGTH OF IW1 WHICH MUST BE AT LEAST THE ABSOLUTE
!     VALUE OF IW(1).  IT IS NOT ALTERED.
! ICNTL is an INTEGER array of length 30, see MA27A/AD.
! INFO is an INTEGER array of length 20, see MA27A/AD.
!
!     .. Scalar Arguments ..
      integer la,liw,maxfrt,n,nsteps
!     ..
!     .. Array Arguments ..
      double precision a(la),rhs(n),w(maxfrt)
      integer iw(liw),iw1(nsteps),icntl(30),info(20)
!     ..
!     .. Local Scalars ..
      integer i,iapos,iblk,ipos,irows,j1,j2,jj,k,kblk,latop,len,nblk, &
              ncols,nrows
!     ..
!     .. External Subroutines ..
      external ma27qd,ma27rd
!     ..
!     .. Intrinsic Functions ..
      intrinsic abs,min
!     ..
!     .. Executable Statements ..
      info(1) = 0

      if (icntl(3).le.0 .or. icntl(2).le.0) go to 110
! PRINT INPUT PARAMETERS.
      write (icntl(2),fmt=10) n,la,liw,maxfrt,nsteps

   10 format (/,/,' ENTERING MA27CD WITH      N     LA    LIW MAXFRT', &
             '  NSTEPS',/,21x,5i7)
! PRINT OUT MATRIX FACTORS FROM MA27B/BD.
      kblk = abs(iw(1)+0)
      if (kblk.eq.0) go to 90
      if (icntl(3).eq.1) kblk = 1
      ipos = 2
      iapos = 1
      do 80 iblk = 1,kblk
        ncols = iw(ipos)
        nrows = iw(ipos+1)
        j1 = ipos + 2
        if (ncols.gt.0) go to 20
        ncols = -ncols
        nrows = 1
        j1 = j1 - 1
   20   write (icntl(2),fmt=30) iblk,nrows,ncols

   30   format (' BLOCK PIVOT =',i8,' NROWS =',i8,' NCOLS =',i8)

        j2 = j1 + ncols - 1
        ipos = j2 + 1
        write (icntl(2),fmt=40) (iw(jj),jj=j1,j2)

   40   format (' COLUMN INDICES =',10i6,/, (17x,10i6))

        write (icntl(2),fmt=50)

   50   format (' REAL ENTRIES .. EACH ROW STARTS ON A NEW LINE')

        len = ncols
        do 70 irows = 1,nrows
          j1 = iapos
          j2 = iapos + len - 1
          write (icntl(2),fmt=60) (a(jj),jj=j1,j2)

   60     format (1p,5d13.3)

          len = len - 1
          iapos = j2 + 1
   70   continue
   80 continue
   90 k = min(10,n)
      if (icntl(3).gt.1) k = n
      if (n.gt.0) write (icntl(2),fmt=100) (rhs(i),i=1,k)

  100 format (' RHS',1p,5d13.3,/, (4x,1p,5d13.3))

  110 if (iw(1).lt.0) go to 130
      nblk = iw(1)
      if (nblk.gt.0) go to 140
! WE HAVE A ZERO MATRIX
      do 120 i = 1,n
        rhs(i) = 0.0d0
  120 continue
      go to 150

  130 nblk = -iw(1)
! FORWARD SUBSTITUTION
  140 call ma27qd(n,a,la,iw(2),liw-1,w,maxfrt,rhs,iw1,nblk,latop,icntl)
! BACK SUBSTITUTION.
      call ma27rd(n,a,la,iw(2),liw-1,w,maxfrt,rhs,iw1,nblk,latop,icntl)
  150 if (icntl(3).le.0 .or. icntl(2).le.0) go to 170
! PRINT OUTPUT PARAMETERS.
      write (icntl(2),fmt=160)

  160 format (/,/,' LEAVING MA27CD WITH')

      if (n.gt.0) write (icntl(2),fmt=100) (rhs(i),i=1,k)
  170 continue

      return
      end
      subroutine ma27gd(n,nz,irn,icn,iw,lw,ipe,iq,flag,iwfr, &
                       icntl,info)
!
! SORT PRIOR TO CALLING ANALYSIS ROUTINE MA27H/HD.
!
! GIVEN THE POSITIONS OF THE OFF-DIAGONAL NON-ZEROS OF A SYMMETRIC
!     MATRIX, CONSTRUCT THE SPARSITY PATTERN OF THE OFF-DIAGONAL
!     PART OF THE WHOLE MATRIX (UPPER AND LOWER TRIANGULAR PARTS).
! EITHER ONE OF A PAIR (I,J),(J,I) MAY BE USED TO REPRESENT
!     THE PAIR. DIAGONAL ELEMENTS AND DUPLICATES ARE IGNORED.
!
! N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
! NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
!     ALTERED.
! IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW NUMBERS OF THE
!     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
!     TO IW (SEE DESCRIPTION OF IW).
! ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN NUMBERS OF THE
!     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
!     TO IW (SEE DESCRIPTION OF IW).
! IW NEED NOT BE SET ON INPUT. ON OUTPUT IT CONTAINS LISTS OF
!     COLUMN INDICES, EACH LIST BEING HEADED BY ITS LENGTH.
!     IRN(1) MAY BE EQUIVALENCED TO IW(1) AND ICN(1) MAY BE
!     EQUIVALENCED TO IW(K), WHERE K.GT.NZ.
! LW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST 2*NZ+N.
!     IT IS NOT ALTERED.
! IPE NEED NOT BE SET ON INPUT. ON OUTPUT IPE(I) POINTS TO THE START OF
!     THE ENTRY IN IW FOR ROW I, OR IS ZERO IF THERE IS NO ENTRY.
! IQ NEED NOT BE SET.  ON OUTPUT IQ(I),I=1,N CONTAINS THE NUMBER OF
!     OFF-DIAGONAL N0N-ZEROS IN ROW I INCLUDING DUPLICATES.
! FLAG IS USED FOR WORKSPACE TO HOLD FLAGS TO PERMIT DUPLICATE ENTRIES
!     TO BE IDENTIFIED QUICKLY.
! IWFR NEED NOT BE SET ON INPUT. ON OUTPUT IT POINTS TO THE FIRST
!     UNUSED LOCATION IN IW.
! ICNTL is an INTEGER array of length 30, see MA27A/AD.
! INFO is an INTEGER array of length 20, see MA27A/AD.
!
!     .. Scalar Arguments ..
      integer iwfr,lw,n,nz
!     ..
!     .. Array Arguments ..
      integer flag(n),icn(*),ipe(n),iq(n),irn(*),iw(lw)
      integer icntl(*),info(20)
!     ..
!     .. Local Scalars ..
      integer i,id,j,jn,k,k1,k2,l,last,lr,n1,ndup
!     ..
!     .. Executable Statements ..
!
! INITIALIZE INFO(2) AND COUNT IN IPE THE
!     NUMBERS OF NON-ZEROS IN THE ROWS AND MOVE ROW AND COLUMN
!     NUMBERS INTO IW.
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
          write (icntl(2),fmt=60) info(1)
        end if

   60   format (' *** WARNING MESSAGE FROM SUBROUTINE MA27AD', &
                '  *** INFO(1) =',i2)

        if (info(2).le.10 .and. icntl(2).gt.0) then
          write (icntl(2),fmt=70) k,i,j
        end if

   70   format (i6,'TH NON-ZERO (IN ROW',i6,' AND COLUMN',i6, &
               ') IGNORED')

   80   i = 0
        j = 0
        go to 100

   90   ipe(i) = ipe(i) + 1
        ipe(j) = ipe(j) + 1
  100   iw(k) = j
        lr = lr + 1
        iw(lr) = i
  110 continue
!
! ACCUMULATE ROW COUNTS TO GET POINTERS TO ROW STARTS IN BOTH IPE AND IQ
!     AND INITIALIZE FLAG
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
!
! RUN THROUGH PUTTING THE MATRIX ELEMENTS IN THE RIGHT PLACE
!     BUT WITH SIGNS INVERTED. IQ IS USED FOR HOLDING RUNNING POINTERS
!     AND IS LEFT HOLDING POINTERS TO ROW ENDS.
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
!
! RUN THROUGH RESTORING SIGNS, REMOVING DUPLICATES AND SETTING THE
!     MATE OF EACH NON-ZERO.
! NDUP COUNTS THE NUMBER OF DUPLICATE ELEMENTS.
  230 ndup = 0
      do 280 i = 1,n
        k1 = ipe(i) + 1
        k2 = iq(i)
        if (k1.le.k2) go to 240
! ROW IS EMPTY. SET POINTER TO ZERO.
        ipe(i) = 0
        iq(i) = 0
        go to 280
! ON ENTRY TO THIS LOOP FLAG(J).LT.I FOR J=1,2,...,N. DURING THE LOOP
!     FLAG(J) IS SET TO I IF A NON-ZERO IN COLUMN J IS FOUND. THIS
!     PERMITS DUPLICATES TO BE RECOGNIZED QUICKLY.
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
!
! COMPRESS IW TO REMOVE DUMMY ENTRIES CAUSED BY DUPLICATES.
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
      subroutine ma27hd(n,ipe,iw,lw,iwfr,nv,nxt,lst,ipd,flag,iovflo, &
                       ncmpa,fratio)
!
! ANALYSIS SUBROUTINE
!
! GIVEN REPRESENTATION OF THE WHOLE MATRIX (EXCLUDING DIAGONAL)
!     PERFORM MINIMUM DEGREE ORDERING, CONSTRUCTING TREE POINTERS.
!     IT WORKS WITH SUPERVARIABLES WHICH ARE COLLECTIONS OF ONE OR MORE
!     VARIABLES, STARTING WITH SUPERVARIABLE I CONTAINING VARIABLE I FOR
!     I = 1,2,...,N. ALL VARIABLES IN A SUPERVARIABLE ARE ELIMINATED
!     TOGETHER. EACH SUPERVARIABLE HAS AS NUMERICAL NAME THAT OF ONE
!     OF ITS VARIABLES (ITS PRINCIPAL VARIABLE).
!
! N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
! IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
!     START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
!     DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS. IF
!     SUPERVARIABLE I IS ABSORBED INTO SUPERVARIABLE J THEN IPE(I)=-J.
!     IF SUPERVARIABLE I IS ELIMINATED THEN IPE(I) EITHER POINTS TO THE
!     LIST OF SUPERVARIABLES FOR CREATED ELEMENT I OR IS ZERO IF
!     THE CREATED ELEMENT IS NULL. IF ELEMENT I
!     IS ABSORBED INTO ELEMENT J THEN IPE(I)=-J.
! IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
!     ROWS, EACH LIST BEING HEADED BY ITS LENGTH.
!     DURING EXECUTION THESE LISTS ARE REVISED AND HOLD
!     LISTS OF ELEMENTS AND SUPERVARIABLES. THE ELEMENTS
!     ALWAYS HEAD THE LISTS. WHEN A SUPERVARIABLE
!     IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF SUPERVARIABLES
!     IN THE NEW ELEMENT.
! LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
! IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
!     IT IS REVISED DURING EXECUTION AND CONTINUES TO HAVE THIS MEANING.
! NV(JS) NEED NOT BE SET. DURING EXECUTION IT IS ZERO IF
!     JS IS NOT A PRINCIPAL VARIABLE AND IF IT IS IT HOLDS
!     THE NUMBER OF VARIABLES IN SUPERVARIABLE JS. FOR ELIMINATED
!     VARIABLES IT IS SET TO THE DEGREE AT THE TIME OF ELIMINATION.
! NXT(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE NEXT
!     SUPERVARIABLE HAVING THE SAME DEGREE AS JS, OR ZERO
!     IF IT IS LAST IN ITS LIST.
! LST(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE
!     LAST SUPERVARIABLE HAVING THE SAME DEGREE AS JS OR
!     -(ITS DEGREE) IF IT IS FIRST IN ITS LIST.
! IPD(ID) NEED NOT BE SET. DURING EXECUTION IT
!     IS THE FIRST SUPERVARIABLE WITH DEGREE ID OR ZERO
!     IF THERE ARE NONE.
! FLAG IS USED AS WORKSPACE FOR ELEMENT AND SUPERVARIABLE FLAGS.
!     WHILE THE CODE IS FINDING THE DEGREE OF SUPERVARIABLE IS
!     FLAG HAS THE FOLLOWING VALUES.
!     A) FOR THE CURRENT PIVOT/NEW ELEMENT ME
!           FLAG(ME)=-1
!     B) FOR VARIABLES JS
!           FLAG(JS)=-1 IF JS IS NOT A PRINCIPAL VARIABLE
!           FLAG(JS)=0 IF JS IS A SUPERVARIABLE IN THE NEW ELEMENT
!           FLAG(JS)=NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
!                 ELEMENT THAT HAS BEEN COUNTED IN THE DEGREE
!                 CALCULATION
!           FLAG(JS).GT.NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
!                 ELEMENT THAT HAS NOT BEEN COUNTED IN THE DEGREE
!                 CALCULATION
!     C) FOR ELEMENTS IE
!           FLAG(IE)=-1 IF ELEMENT IE HAS BEEN MERGED INTO ANOTHER
!           FLAG(IE)=-NFLG IF ELEMENT IE HAS BEEN USED IN THE DEGREE
!                 CALCULATION FOR IS.
!           FLAG(IE).LT.-NFLG IF ELEMENT IE HAS NOT BEEN USED IN THE
!                 DEGREE CALCULATION FOR IS
! IOVFLO see ICNTL(4) in MA27A/AD.
! NCMPA see INFO(11) in MA27A/AD.
! FRATIO see CNTL(2) in MA27A/AD.
!
!     .. Scalar Arguments ..
      double precision fratio
      integer iwfr,lw,n,iovflo,ncmpa
!     ..
!     .. Array Arguments ..
      integer flag(n),ipd(n),ipe(n),iw(lw),lst(n),nv(n),nxt(n)
!     ..
!     .. Local Scalars ..
      integer i,id,idl,idn,ie,ip,is,jp,jp1,jp2,js,k,k1,k2,ke,kp,kp0,kp1, &
              kp2,ks,l,len,limit,ln,ls,lwfr,md,me,ml,ms,nel,nflg,np, &
              np0,ns,nvpiv,nvroot,root
! LIMIT  Limit on number of variables for putting node in root.
! NVROOT Number of variables in the root node
! ROOT   Index of the root node (N+1 if none chosen yet).
!     ..
!     .. External Subroutines ..
      external ma27ud
!     ..
!     .. Intrinsic Functions ..
      intrinsic abs,min
!     ..
! If a column of the reduced matrix has relative density greater than
! CNTL(2), it is forced into the root. All such columns are taken to
! have sparsity pattern equal to their merged patterns, so the fill
! and operation counts may be overestimated.
!
! IS,JS,KS,LS,MS,NS ARE USED TO REFER TO SUPERVARIABLES.
! IE,JE,KE ARE USED TO REFER TO ELEMENTS.
! IP,JP,KP,K,NP ARE USED TO POINT TO LISTS OF ELEMENTS.
!     OR SUPERVARIABLES.
! ID IS USED FOR THE DEGREE OF A SUPERVARIABLE.
! MD IS USED FOR THE CURRENT MINIMUM DEGREE.
! IDN IS USED FOR THE NO. OF VARIABLES IN A NEWLY CREATED ELEMENT
! NEL IS USED TO HOLD THE NO. OF VARIABLES THAT HAVE BEEN
!     ELIMINATED.
! ME=MS IS THE NAME OF THE SUPERVARIABLE ELIMINATED AND
!     OF THE ELEMENT CREATED IN THE MAIN LOOP.
! NFLG IS USED FOR THE CURRENT FLAG VALUE IN ARRAY FLAG. IT STARTS
!     WITH THE VALUE IOVFLO AND IS REDUCED BY 1 EACH TIME IT IS USED
!     UNTIL IT HAS THE VALUE 2 WHEN IT IS RESET TO THE VALUE IOVFLO.
!
!     .. Executable Statements ..
! INITIALIZATIONS
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
!
! LINK TOGETHER VARIABLES HAVING SAME DEGREE
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
! WE HAVE A VARIABLE THAT CAN BE ELIMINATED AT ONCE BECAUSE THERE IS
!     NO OFF-DIAGONAL NON-ZERO IN ITS ROW.
   20   nel = nel + 1
        flag(is) = -1
        nxt(is) = 0
        lst(is) = 0
   30 continue
!
! START OF MAIN LOOP
!
      do 340 ml = 1,n
! LEAVE LOOP IF ALL VARIABLES HAVE BEEN ELIMINATED.
        if (nel+nvroot+1.ge.n) go to 350
!
! FIND NEXT SUPERVARIABLE FOR ELIMINATION.
        do 40 id = md,n
          ms = ipd(id)
          if (ms.gt.0) go to 50
   40   continue
   50   md = id
! NVPIV HOLDS THE NUMBER OF VARIABLES IN THE PIVOT.
        nvpiv = nv(ms)
!
! REMOVE CHOSEN VARIABLE FROM LINKED LIST
        ns = nxt(ms)
        nxt(ms) = 0
        lst(ms) = 0
        if (ns.gt.0) lst(ns) = -id
        ipd(id) = ns
        me = ms
        nel = nel + nvpiv
! IDN HOLDS THE DEGREE OF THE NEW ELEMENT.
        idn = 0
!
! RUN THROUGH THE LIST OF THE PIVOTAL SUPERVARIABLE, SETTING TREE
!     POINTERS AND CONSTRUCTING NEW LIST OF SUPERVARIABLES.
! KP IS A POINTER TO THE CURRENT POSITION IN THE OLD LIST.
        kp = ipe(me)
        flag(ms) = -1
! IP POINTS TO THE START OF THE NEW LIST.
        ip = iwfr
! LEN HOLDS THE LENGTH OF THE LIST ASSOCIATED WITH THE PIVOT.
        len = iw(kp)
        do 140 kp1 = 1,len
          kp = kp + 1
          ke = iw(kp)
! JUMP IF KE IS AN ELEMENT THAT HAS NOT BEEN MERGED INTO ANOTHER.
          if (flag(ke).le.-2) go to 60
! JUMP IF KE IS AN ELEMENT THAT HAS BEEN MERGED INTO ANOTHER OR IS
!     A SUPERVARIABLE THAT HAS BEEN ELIMINATED.
          if (flag(ke).le.0) then
             if (ipe(ke).ne.-root) go to 140
! KE has been merged into the root
             ke = root
             if (flag(ke).le.0) go to 140
          end if
! WE HAVE A SUPERVARIABLE. PREPARE TO SEARCH REST OF LIST.
          jp = kp - 1
          ln = len - kp1 + 1
          ie = ms
          go to 70
! SEARCH VARIABLE LIST OF ELEMENT KE, USING JP AS A POINTER TO IT.
   60     ie = ke
          jp = ipe(ie)
          ln = iw(jp)
!
! SEARCH FOR DIFFERENT SUPERVARIABLES AND ADD THEM TO THE NEW LIST,
!     COMPRESSING WHEN NECESSARY. THIS LOOP IS EXECUTED ONCE FOR
!     EACH ELEMENT IN THE LIST AND ONCE FOR ALL THE SUPERVARIABLES
!     IN THE LIST.
   70     do 130 jp1 = 1,ln
            jp = jp + 1
            is = iw(jp)
! JUMP IF IS IS NOT A PRINCIPAL VARIABLE OR HAS ALREADY BEEN COUNTED.
            if (flag(is).le.0) then
               if (ipe(is).eq.-root) then
! IS has been merged into the root
                  is = root
                  iw(jp) = root
                  if (flag(is).le.0) go to 130
               else
                  go to 130
               end if
            end if
            flag(is) = 0
            if (iwfr.lt.lw) go to 100
! PREPARE FOR COMPRESSING IW BY ADJUSTING POINTERS AND
!     LENGTHS SO THAT THE LISTS BEING SEARCHED IN THE INNER AND OUTER
!     LOOPS CONTAIN ONLY THE REMAINING ENTRIES.
            ipe(ms) = kp
            iw(kp) = len - kp1
            ipe(ie) = jp
            iw(jp) = ln - jp1
! COMPRESS IW
            call ma27ud(n,ipe,iw,ip-1,lwfr,ncmpa)
! COPY NEW LIST FORWARD
            jp2 = iwfr - 1
            iwfr = lwfr
            if (ip.gt.jp2) go to 90
            do 80 jp = ip,jp2
              iw(iwfr) = iw(jp)
              iwfr = iwfr + 1
   80       continue
! ADJUST POINTERS FOR THE NEW LIST AND THE LISTS BEING SEARCHED.
   90       ip = lwfr
            jp = ipe(ie)
            kp = ipe(me)
! STORE IS IN NEW LIST.
  100       iw(iwfr) = is
            idn = idn + nv(is)
            iwfr = iwfr + 1
! REMOVE IS FROM DEGREE LINKED LIST
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
! JUMP IF WE HAVE JUST BEEN SEARCHING THE VARIABLES AT THE END OF
!     THE LIST OF THE PIVOT.
          if (ie.eq.ms) go to 150
! SET TREE POINTER AND FLAG TO INDICATE ELEMENT IE IS ABSORBED INTO
!     NEW ELEMENT ME.
          ipe(ie) = -me
          flag(ie) = -1
  140   continue

! STORE THE DEGREE OF THE PIVOT.
  150   nv(ms) = idn + nvpiv
! JUMP IF NEW ELEMENT IS NULL.
        if (iwfr.eq.ip) go to 330
        k1 = ip
        k2 = iwfr - 1
!
! RUN THROUGH NEW LIST OF SUPERVARIABLES REVISING EACH ASSOCIATED LIST,
!     RECALCULATING DEGREES AND REMOVING DUPLICATES.
        limit = nint(fratio*(n-nel))
        do 310 k = k1,k2
          is = iw(k)
          if (is.eq.root) go to 310
          if (nflg.gt.2) go to 170
! RESET FLAG VALUES TO +/-IOVFLO.
          do 160 i = 1,n
            if (flag(i).gt.0) flag(i) = iovflo
            if (flag(i).le.-2) flag(i) = -iovflo
  160     continue
          nflg = iovflo
! REDUCE NFLG BY ONE TO CATER FOR THIS SUPERVARIABLE.
  170     nflg = nflg - 1
! BEGIN WITH THE DEGREE OF THE NEW ELEMENT. ITS VARIABLES MUST ALWAYS
!     BE COUNTED DURING THE DEGREE CALCULATION AND THEY ARE ALREADY
!     FLAGGED WITH THE VALUE 0.
          id = idn
! RUN THROUGH THE LIST ASSOCIATED WITH SUPERVARIABLE IS
          kp1 = ipe(is) + 1
! NP POINTS TO THE NEXT ENTRY IN THE REVISED LIST.
          np = kp1
          kp2 = iw(kp1-1) + kp1 - 1
          do 220 kp = kp1,kp2
            ke = iw(kp)
! TEST WHETHER KE IS AN ELEMENT, A REDUNDANT ENTRY OR A SUPERVARIABLE.
          if (flag(ke).eq.-1) then
             if (ipe(ke).ne.-root) go to 220
! KE has been merged into the root
             ke = root
             iw(kp) = root
             if (flag(ke).eq.-1) go to 220
          end if
          if (flag(ke).ge.0) go to 230
! SEARCH LIST OF ELEMENT KE, REVISING THE DEGREE WHEN NEW VARIABLES
!     FOUND.
            jp1 = ipe(ke) + 1
            jp2 = iw(jp1-1) + jp1 - 1
            idl = id
            do 190 jp = jp1,jp2
              js = iw(jp)
! JUMP IF JS HAS ALREADY BEEN COUNTED.
              if (flag(js).le.nflg) go to 190
              id = id + nv(js)
              flag(js) = nflg
  190       continue
! JUMP IF ONE OR MORE NEW SUPERVARIABLES WERE FOUND.
            if (id.gt.idl) go to 210
! CHECK WHETHER EVERY VARIABLE OF ELEMENT KE IS IN NEW ELEMENT ME.
            do 200 jp = jp1,jp2
              js = iw(jp)
              if (flag(js).ne.0) go to 210
  200       continue
! SET TREE POINTER AND FLAG TO INDICATE THAT ELEMENT KE IS ABSORBED
!     INTO NEW ELEMENT ME.
            ipe(ke) = -me
            flag(ke) = -1
            go to 220
! STORE ELEMENT KE IN THE REVISED LIST FOR SUPERVARIABLE IS AND FLAG IT.
  210       iw(np) = ke
            flag(ke) = -nflg
            np = np + 1
  220     continue
          np0 = np
          go to 250
! TREAT THE REST OF THE LIST ASSOCIATED WITH SUPERVARIABLE IS. IT
!     CONSISTS ENTIRELY OF SUPERVARIABLES.
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
! ADD TO DEGREE, FLAG SUPERVARIABLE KS AND ADD IT TO NEW LIST.
            id = id + nv(ks)
            flag(ks) = nflg
            iw(np) = ks
            np = np + 1
  240     continue
! MOVE FIRST SUPERVARIABLE TO END OF LIST, MOVE FIRST ELEMENT TO END
!     OF ELEMENT PART OF LIST AND ADD NEW ELEMENT TO FRONT OF LIST.
  250     if (id.ge.limit) go to 295
          iw(np) = iw(np0)
          iw(np0) = iw(kp1)
          iw(kp1) = me
! STORE THE NEW LENGTH OF THE LIST.
          iw(kp1-1) = np - kp1 + 1
!
! CHECK WHETHER ROW IS IS IDENTICAL TO ANOTHER BY LOOKING IN LINKED
!     LIST OF SUPERVARIABLES WITH DEGREE ID AT THOSE WHOSE LISTS HAVE
!     FIRST ENTRY ME. NOTE THAT THOSE CONTAINING ME COME FIRST SO THE
!     SEARCH CAN BE TERMINATED WHEN A LIST NOT STARTING WITH ME IS
!     FOUND.
          js = ipd(id)
          do 280 l = 1,n
            if (js.le.0) go to 300
            kp1 = ipe(js) + 1
            if (iw(kp1).ne.me) go to 300
! JS HAS SAME DEGREE AND IS ACTIVE. CHECK IF IDENTICAL TO IS.
            kp2 = kp1 - 1 + iw(kp1-1)
            do 260 kp = kp1,kp2
              ie = iw(kp)
! JUMP IF IE IS A SUPERVARIABLE OR AN ELEMENT NOT IN THE LIST OF IS.
              if (abs(flag(ie)+0).gt.nflg) go to 270
  260       continue
            go to 290

  270       js = nxt(js)
  280     continue
! SUPERVARIABLE AMALGAMATION. ROW IS IS IDENTICAL TO ROW JS.
! REGARD ALL VARIABLES IN THE TWO SUPERVARIABLES AS BEING IN IS. SET
!     TREE POINTER, FLAG AND NV ENTRIES.
  290     ipe(js) = -is
          nv(is) = nv(is) + nv(js)
          nv(js) = 0
          flag(js) = -1
! REPLACE JS BY IS IN LINKED LIST.
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
! Treat IS as full. Merge it into the root node.
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
! INSERT IS INTO LINKED LIST OF SUPERVARIABLES OF SAME DEGREE.
  300     ns = ipd(id)
          if (ns.gt.0) lst(ns) = is
          nxt(is) = ns
          ipd(id) = is
          lst(is) = -id
          md = min(md,id)
  310   continue
!
! RESET FLAGS FOR SUPERVARIABLES IN NEWLY CREATED ELEMENT AND
!     REMOVE THOSE ABSORBED INTO OTHERS.
        do 320 k = k1,k2
          is = iw(k)
          if (nv(is).eq.0) go to 320
          flag(is) = nflg
          iw(ip) = is
          ip = ip + 1
  320   continue
        iwfr = k1
        flag(me) = -nflg
! MOVE FIRST ENTRY TO END TO MAKE ROOM FOR LENGTH.
        iw(ip) = iw(k1)
        iw(k1) = ip - k1
! SET POINTER FOR NEW ELEMENT AND RESET IWFR.
        ipe(me) = k1
        iwfr = ip + 1
        go to 335

  330   ipe(me) = 0
!
  335   continue
  340 continue
!

! Absorb any remaining variables into the root
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
! Link any remaining elements to the root
      do 370 ie = 1,n
        if (ipe(ie).gt.0) ipe(ie) = -root
  370 continue
      if(nvroot.gt.0)nv(root)=nvroot
      end

      subroutine ma27ud(n,ipe,iw,lw,iwfr,ncmpa)
! COMPRESS LISTS HELD BY MA27H/HD AND MA27K/KD IN IW AND ADJUST POINTERS
!     IN IPE TO CORRESPOND.
! N IS THE MATRIX ORDER. IT IS NOT ALTERED.
! IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
!     ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
! IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
!     LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
! LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
! IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
!     LOCATION IN IW.
!     ON RETURN IT IS SET TO THE FIRST FREE LOCATION IN IW.
! NCMPA see INFO(11) in MA27A/AD.
!
!     .. Scalar Arguments ..
      integer iwfr,lw,n,ncmpa
!     ..
!     .. Array Arguments ..
      integer ipe(n),iw(lw)
!     ..
!     .. Local Scalars ..
      integer i,ir,k,k1,k2,lwfr
!     ..
!     .. Executable Statements ..
      ncmpa = ncmpa + 1
! PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
!     LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
!     -(LIST NUMBER).
      do 10 i = 1,n
        k1 = ipe(i)
        if (k1.le.0) go to 10
        ipe(i) = iw(k1)
        iw(k1) = -i
   10 continue
!
! COMPRESS
! IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
! LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.
      iwfr = 1
      lwfr = iwfr
      do 60 ir = 1,n
        if (lwfr.gt.lw) go to 70
! SEARCH FOR THE NEXT NEGATIVE ENTRY.
        do 20 k = lwfr,lw
          if (iw(k).lt.0) go to 30
   20   continue
        go to 70
! PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW POINTER
!     AND PREPARE TO COPY LIST.
   30   i = -iw(k)
        iw(iwfr) = ipe(i)
        ipe(i) = iwfr
        k1 = k + 1
        k2 = k + iw(iwfr)
        iwfr = iwfr + 1
        if (k1.gt.k2) go to 50
! COPY LIST TO NEW POSITION.
        do 40 k = k1,k2
          iw(iwfr) = iw(k)
          iwfr = iwfr + 1
   40   continue
   50   lwfr = k2 + 1
   60 continue
   70 return

      end
      subroutine ma27jd(n,nz,irn,icn,perm,iw,lw,ipe,iq,flag,iwfr, &
                       icntl,info)
!
! SORT PRIOR TO CALLING ANALYSIS ROUTINE MA27K/KD.
!
! GIVEN THE POSITIONS OF THE OFF-DIAGONAL NON-ZEROS OF A SYMMETRIC
!     MATRIX AND A PERMUTATION, CONSTRUCT THE SPARSITY PATTERN
!     OF THE STRICTLY UPPER TRIANGULAR PART OF THE PERMUTED MATRIX.
!     EITHER ONE OF A PAIR (I,J),(J,I) MAY BE USED TO REPRESENT
!     THE PAIR. DIAGONAL ELEMENTS ARE IGNORED. NO CHECK IS MADE
!     FOR DUPLICATE ELEMENTS UNLESS ANY ROW HAS MORE THAN ICNTL(4)
!     NON-ZEROS, IN WHICH CASE DUPLICATES ARE REMOVED.
!
! N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
! NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
!     ALTERED.
! IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW INDICES OF THE
!     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS EQUIVALENCED WITH IW.
!     IRN(1) MAY BE EQUIVALENCED WITH IW(1).
! ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN INDICES OF THE
!     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS EQUIVALENCED
!     WITH IW.ICN(1) MAY BE EQUIVELENCED WITH IW(K),K.GT.NZ.
! PERM(I) MUST BE SET TO HOLD THE POSITION OF VARIABLE I IN THE
!     PERMUTED ORDER. IT IS NOT ALTERED.
! IW NEED NOT BE SET ON INPUT. ON OUTPUT IT CONTAINS LISTS OF
!     COLUMN NUMBERS, EACH LIST BEING HEADED BY ITS LENGTH.
! LW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST
!     MAX(NZ,N+(NO. OF OFF-DIAGONAL NON-ZEROS)). IT IS NOT ALTERED.
! IPE NEED NOT BE SET ON INPUT. ON OUTPUT IPE(I) POINTS TO THE START OF
!     THE ENTRY IN IW FOR ROW I, OR IS ZERO IF THERE IS NO ENTRY.
! IQ NEED NOT BE SET. ON OUTPUT IQ(I) CONTAINS THE NUMBER OF
!     OFF-DIAGONAL NON-ZEROS IN ROW I, INCLUDING DUPLICATES.
! FLAG IS USED FOR WORKSPACE TO HOLD FLAGS TO PERMIT DUPLICATE
!     ENTRIES TO BE IDENTIFIED QUICKLY.
! IWFR NEED NOT BE SET ON INPUT. ON OUTPUT IT POINTS TO THE FIRST
!     UNUSED LOCATION IN IW.
! ICNTL is an INTEGER array of length 30, see MA27A/AD.
! INFO is an INTEGER array of length 20, see MA27A/AD.
!
!     .. Scalar Arguments ..
      integer iwfr,lw,n,nz
!     ..
!     .. Array Arguments ..
      integer flag(n),icn(*),ipe(n),iq(n),irn(*),iw(lw),perm(n)
      integer icntl(30),info(20)
!     ..
!     .. Local Scalars ..
      integer i,id,in,j,jdummy,k,k1,k2,l,lbig,len
!     ..
!     .. Intrinsic Functions ..
      intrinsic max
!     ..
!     .. Executable Statements ..
!
! INITIALIZE INFO(1), INFO(2) AND IQ
      info(1) = 0
      info(2) = 0
      do 10 i = 1,n
        iq(i) = 0
   10 continue
!
! COUNT THE NUMBERS OF NON-ZEROS IN THE ROWS, PRINT WARNINGS ABOUT
!     OUT-OF-RANGE INDICES AND TRANSFER GENUINE ROW NUMBERS
!     (NEGATED) INTO IW.
      if (nz.eq.0) go to 110
      do 100 k = 1,nz
        i = irn(k)
        j = icn(k)
        iw(k) = -i
        if(i.lt.j) then
          if (i.ge.1 .and. j.le.n) go to 80
        else if(i.gt.j) then
          if (j.ge.1 .and. i.le.n) go to 80
        else
          iw(k) = 0
          if (i.ge.1 .and. i.le.n) go to 100
        end if
        info(2) = info(2) + 1
        info(1) = 1
        iw(k) = 0
        if (info(2).le.1 .and. icntl(2).gt.0) then
          write (icntl(2),fmt=60) info(1)
        end if

   60   format (' *** WARNING MESSAGE FROM SUBROUTINE MA27AD', &
                '  *** INFO(1) =',i2)

        if (info(2).le.10 .and. icntl(2).gt.0) then
          write (icntl(2),fmt=70) k,i,j
        end if

   70   format (i6,'TH NON-ZERO (IN ROW',i6,' AND COLUMN',i6, &
               ') IGNORED')

        go to 100

   80   if (perm(j).gt.perm(i)) go to 90
        iq(j) = iq(j) + 1
        go to 100

   90   iq(i) = iq(i) + 1
  100 continue
!
! ACCUMULATE ROW COUNTS TO GET POINTERS TO ROW ENDS
!     IN IPE.
  110 iwfr = 1
      lbig = 0
      do 120 i = 1,n
        l = iq(i)
        lbig = max(l,lbig)
        iwfr = iwfr + l
        ipe(i) = iwfr - 1
  120 continue
!
! PERFORM IN-PLACE SORT
      if (nz.eq.0) go to 250
      do 160 k = 1,nz
        i = -iw(k)
        if (i.le.0) go to 160
        l = k
        iw(k) = 0
        do 150 id = 1,nz
          j = icn(l)
          if (perm(i).lt.perm(j)) go to 130
          l = ipe(j)
          ipe(j) = l - 1
          in = iw(l)
          iw(l) = i
          go to 140

  130     l = ipe(i)
          ipe(i) = l - 1
          in = iw(l)
          iw(l) = j
  140     i = -in
          if (i.le.0) go to 160
  150   continue
  160 continue
!
! MAKE ROOM IN IW FOR ROW LENGTHS AND INITIALIZE FLAG.
      k = iwfr - 1
      l = k + n
      iwfr = l + 1
      do 190 i = 1,n
        flag(i) = 0
        j = n + 1 - i
        len = iq(j)
        if (len.le.0) go to 180
        do 170 jdummy = 1,len
          iw(l) = iw(k)
          k = k - 1
          l = l - 1
  170   continue
  180   ipe(j) = l
        l = l - 1
  190 continue
      if (lbig.ge.icntl(4)) go to 210
!
! PLACE ROW LENGTHS IN IW
      do 200 i = 1,n
        k = ipe(i)
        iw(k) = iq(i)
        if (iq(i).eq.0) ipe(i) = 0
  200 continue
      go to 250
!
!
! REMOVE DUPLICATE ENTRIES
  210 iwfr = 1
      do 240 i = 1,n
        k1 = ipe(i) + 1
        k2 = ipe(i) + iq(i)
        if (k1.le.k2) go to 220
        ipe(i) = 0
        go to 240

  220   ipe(i) = iwfr
        iwfr = iwfr + 1
        do 230 k = k1,k2
          j = iw(k)
          if (flag(j).eq.i) go to 230
          iw(iwfr) = j
          iwfr = iwfr + 1
          flag(j) = i
  230   continue
        k = ipe(i)
        iw(k) = iwfr - k - 1
  240 continue
  250 return

      end
      subroutine ma27kd(n,ipe,iw,lw,iwfr,ips,ipv,nv,flag,ncmpa)
!
! USING A GIVEN PIVOTAL SEQUENCE AND A REPRESENTATION OF THE MATRIX THAT
!     INCLUDES ONLY NON-ZEROS OF THE STRICTLY UPPER-TRIANGULAR PART
!     OF THE PERMUTED MATRIX, CONSTRUCT TREE POINTERS.
!
! N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
! IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
!     START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
!     DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS.
!     IF VARIABLE I IS ELIMINATED THEN IPE(I) POINTS TO THE LIST
!     OF VARIABLES FOR CREATED ELEMENT I. IF ELEMENT I IS
!     ABSORBED INTO NEWLY CREATED ELEMENT J THEN IPE(I)=-J.
! IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
!     ROWS, EACH LIST BEING HEADED BY ITS LENGTH. WHEN A VARIABLE
!     IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF VARIABLES
!     IN THE NEW ELEMENT.
! LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
! IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
!     IT IS REVISED DURING EXECUTION, CONTINUING TO HAVE THIS MEANING.
! IPS(I) MUST BE SET TO THE POSITION OF VARIABLE I IN THE REQUIRED
!     ORDERING. IT IS NOT ALTERED.
! IPV NEED NOT BE SET. IPV(K) IS SET TO HOLD THE K TH VARIABLE
!     IN PIVOT ORDER.
! NV NEED NOT BE SET. IF VARIABLE J HAS NOT BEEN ELIMINATED THEN
!     THE LAST ELEMENT WHOSE LEADING VARIABLE (VARIABLE EARLIEST
!     IN THE PIVOT SEQUENCE) IS J IS ELEMENT NV(J). IF ELEMENT J
!     EXISTS THEN THE LAST ELEMENT HAVING THE SAME LEADING
!     VARIABLE IS NV(J). IN BOTH CASES NV(J)=0 IF THERE IS NO SUCH
!     ELEMENT. IF ELEMENT J HAS BEEN MERGED INTO A LATER ELEMENT
!     THEN NV(J) IS THE DEGREE AT THE TIME OF ELIMINATION.
! FLAG IS USED AS WORKSPACE FOR VARIABLE FLAGS.
!     FLAG(JS)=ME IF JS HAS BEEN INCLUDED IN THE LIST FOR ME.
! NCMPA see INFO(11) in MA27A/AD.
!
!     .. Scalar Arguments ..
      integer iwfr,lw,n,ncmpa
!     ..
!     .. Array Arguments ..
      integer flag(n),ipe(n),ips(n),ipv(n),iw(lw),nv(n)
!     ..
!     .. Local Scalars ..
      integer i,ie,ip,j,je,jp,jp1,jp2,js,kdummy,ln,lwfr,me,minjs,ml,ms
!     ..
!     .. External Subroutines ..
      external ma27ud
!     ..
!     .. Intrinsic Functions ..
      intrinsic min
!     ..
!     .. Executable Statements ..
!
! INITIALIZATIONS
      do 10 i = 1,n
        flag(i) = 0
        nv(i) = 0
        j = ips(i)
        ipv(j) = i
   10 continue
      ncmpa = 0
!
! START OF MAIN LOOP
!
      do 100 ml = 1,n
! ME=MS IS THE NAME OF THE VARIABLE ELIMINATED AND
!     OF THE ELEMENT CREATED IN THE MAIN LOOP.
        ms = ipv(ml)
        me = ms
        flag(ms) = me
!
! MERGE ROW MS WITH ALL THE ELEMENTS HAVING MS AS LEADING VARIABLE.
! IP POINTS TO THE START OF THE NEW LIST.
        ip = iwfr
! MINJS IS SET TO THE POSITION IN THE ORDER OF THE LEADING VARIABLE
!     IN THE NEW LIST.
        minjs = n
        ie = me
        do 70 kdummy = 1,n
! SEARCH VARIABLE LIST OF ELEMENT IE.
! JP POINTS TO THE CURRENT POSITION IN THE LIST BEING SEARCHED.
          jp = ipe(ie)
! LN IS THE LENGTH OF THE LIST BEING SEARCHED.
          ln = 0
          if (jp.le.0) go to 60
          ln = iw(jp)
!
! SEARCH FOR DIFFERENT VARIABLES AND ADD THEM TO LIST,
!     COMPRESSING WHEN NECESSARY
          do 50 jp1 = 1,ln
            jp = jp + 1
! PLACE NEXT VARIABLE IN JS.
            js = iw(jp)
! JUMP IF VARIABLE HAS ALREADY BEEN INCLUDED.
            if (flag(js).eq.me) go to 50
            flag(js) = me
            if (iwfr.lt.lw) go to 40
! PREPARE FOR COMPRESSING IW BY ADJUSTING POINTER TO AND LENGTH OF
!     THE LIST FOR IE TO REFER TO THE REMAINING ENTRIES.
            ipe(ie) = jp
            iw(jp) = ln - jp1
! COMPRESS IW.
            call ma27ud(n,ipe,iw,ip-1,lwfr,ncmpa)
! COPY NEW LIST FORWARD
            jp2 = iwfr - 1
            iwfr = lwfr
            if (ip.gt.jp2) go to 30
            do 20 jp = ip,jp2
              iw(iwfr) = iw(jp)
              iwfr = iwfr + 1
   20       continue
   30       ip = lwfr
            jp = ipe(ie)
! ADD VARIABLE JS TO NEW LIST.
   40       iw(iwfr) = js
            minjs = min(minjs,ips(js)+0)
            iwfr = iwfr + 1
   50     continue
! RECORD ABSORPTION OF ELEMENT IE INTO NEW ELEMENT.
   60     ipe(ie) = -me
! PICK UP NEXT ELEMENT WITH LEADING VARIABLE MS.
          je = nv(ie)
! STORE DEGREE OF IE.
          nv(ie) = ln + 1
          ie = je
! LEAVE LOOP IF THERE ARE NO MORE ELEMENTS.
          if (ie.eq.0) go to 80
   70   continue
   80   if (iwfr.gt.ip) go to 90
! DEAL WITH NULL NEW ELEMENT.
        ipe(me) = 0
        nv(me) = 1
        go to 100
! LINK NEW ELEMENT WITH OTHERS HAVING SAME LEADING VARIABLE.
   90   minjs = ipv(minjs)
        nv(me) = nv(minjs)
        nv(minjs) = me
! MOVE FIRST ENTRY IN NEW LIST TO END TO ALLOW ROOM FOR LENGTH AT
!     FRONT. SET POINTER TO FRONT.
        iw(iwfr) = iw(ip)
        iw(ip) = iwfr - ip
        ipe(me) = ip
        iwfr = iwfr + 1
  100 continue
      return

      end
      subroutine ma27ld(n,ipe,nv,ips,ne,na,nd,nsteps,nemin)
!
! TREE SEARCH
!
! GIVEN SON TO FATHER TREE POINTERS, PERFORM DEPTH-FIRST
!     SEARCH TO FIND PIVOT ORDER AND NUMBER OF ELIMINATIONS
!     AND ASSEMBLIES AT EACH STAGE.
! N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
! IPE(I) MUST BE SET EQUAL TO -(FATHER OF NODE I) OR ZERO IF
!      NODE IS A ROOT. IT IS ALTERED TO POINT TO ITS NEXT
!      YOUNGER BROTHER IF IT HAS ONE, BUT OTHERWISE IS NOT
!      CHANGED.
! NV(I) MUST BE SET TO ZERO IF NO VARIABLES ARE ELIMINATED AT NODE
!      I AND TO THE DEGREE OTHERWISE. ONLY LEAF NODES CAN HAVE
!      ZERO VALUES OF NV(I). NV IS NOT ALTERED.
! IPS(I) NEED NOT BE SET. IT IS USED TEMPORARILY TO HOLD
!      -(ELDEST SON OF NODE I) IF IT HAS ONE AND 0 OTHERWISE. IT IS
!      EVENTUALLY SET TO HOLD THE POSITION OF NODE I IN THE ORDER.
! NE(IS) NEED NOT BE SET. IT IS SET TO THE NUMBER OF VARIABLES
!      ELIMINATED AT STAGE IS OF THE ELIMINATION.
! NA(IS) NEED NOT BE SET. IT IS SET TO THE NUMBER OF ELEMENTS
!      ASSEMBLED AT STAGE IS OF THE ELIMINATION.
! ND(IS) NEED NOT BE SET. IT IS SET TO THE DEGREE AT STAGE IS OF
!     THE ELIMINATION.
! NSTEPS NEED NOT BE SET. IT IS SET TO  THE NUMBER OF ELIMINATION
!      STEPS.
! NEMIN see ICNTL(5) in MA27A/AD.
!
!     .. Scalar Arguments ..
      integer n,nsteps,nemin
!     ..
!     .. Array Arguments ..
      integer ipe(n),ips(n),na(n),nd(n),ne(n),nv(n)
!     ..
!     .. Local Scalars ..
      integer i,ib,if,il,is,ison,k,l,nr
!     ..
!     .. Executable Statements ..
! INITIALIZE IPS AND NE.
      do 10 i = 1,n
        ips(i) = 0
        ne(i) = 0
   10 continue
!
! SET IPS(I) TO -(ELDEST SON OF NODE I) AND IPE(I) TO NEXT YOUNGER
!     BROTHER OF NODE I IF IT HAS ONE.
! FIRST PASS IS FOR NODES WITHOUT ELIMINATIONS.
      do 20 i = 1,n
        if (nv(i).gt.0) go to 20
        if = -ipe(i)
        is = -ips(if)
        if (is.gt.0) ipe(i) = is
        ips(if) = -i
   20 continue
! NR IS DECREMENTED FOR EACH ROOT NODE. THESE ARE STORED IN
!     NE(I),I=NR,N.
      nr = n + 1
! SECOND PASS TO ADD NODES WITH ELIMINATIONS.
      do 50 i = 1,n
        if (nv(i).le.0) go to 50
! NODE IF IS THE FATHER OF NODE I.
        if = -ipe(i)
        if (if.eq.0) go to 40
        is = -ips(if)
! JUMP IF NODE IF HAS NO SONS YET.
        if (is.le.0) go to 30
! SET POINTER TO NEXT BROTHER
        ipe(i) = is
! NODE I IS ELDEST SON OF NODE IF.
   30   ips(if) = -i
        go to 50
! WE HAVE A ROOT
   40   nr = nr - 1
        ne(nr) = i
   50 continue
!
! DEPTH-FIRST SEARCH.
! IL HOLDS THE CURRENT TREE LEVEL. ROOTS ARE AT LEVEL N, THEIR SONS
!     ARE AT LEVEL N-1, ETC.
! IS HOLDS THE CURRENT ELIMINATION STAGE. WE ACCUMULATE THE NUMBER
!     OF ELIMINATIONS AT STAGE IS DIRECTLY IN NE(IS). THE NUMBER OF
!     ASSEMBLIES IS ACCUMULATED TEMPORARILY IN NA(IL), FOR TREE
!     LEVEL IL, AND IS TRANSFERED TO NA(IS) WHEN WE REACH THE
!     APPROPRIATE STAGE IS.
      is = 1
! I IS THE CURRENT NODE.
      i = 0
      do 160 k = 1,n
        if (i.gt.0) go to 60
! PICK UP NEXT ROOT.
        i = ne(nr)
        ne(nr) = 0
        nr = nr + 1
        il = n
        na(n) = 0
! GO TO SON FOR AS LONG AS POSSIBLE, CLEARING FATHER-SON POINTERS
!     IN IPS AS EACH IS USED AND SETTING NA(IL)=0 FOR ALL LEVELS
!     REACHED.
   60   do 70 l = 1,n
          if (ips(i).ge.0) go to 80
          ison = -ips(i)
          ips(i) = 0
          i = ison
          il = il - 1
          na(il) = 0
   70   continue
! RECORD POSITION OF NODE I IN THE ORDER.
   80   ips(i) = k
        ne(is) = ne(is) + 1
! JUMP IF NODE HAS NO ELIMINATIONS.
        if (nv(i).le.0) go to 120
        if (il.lt.n) na(il+1) = na(il+1) + 1
        na(is) = na(il)
        nd(is) = nv(i)
! CHECK FOR STATIC CONDENSATION
        if (na(is).ne.1) go to 90
        if (nd(is-1)-ne(is-1).eq.nd(is)) go to 100
! CHECK FOR SMALL NUMBERS OF ELIMINATIONS IN BOTH LAST TWO STEPS.
   90   if (ne(is).ge.nemin) go to 110
        if (na(is).eq.0) go to 110
        if (ne(is-1).ge.nemin) go to 110
! COMBINE THE LAST TWO STEPS
  100   na(is-1) = na(is-1) + na(is) - 1
        nd(is-1) = nd(is) + ne(is-1)
        ne(is-1) = ne(is) + ne(is-1)
        ne(is) = 0
        go to 120

  110   is = is + 1
  120   ib = ipe(i)
        if (ib.ge.0) then
! NODE I HAS A BROTHER OR IS A ROOT
          if (ib.gt.0) na(il) = 0
          i = ib
        else
! GO TO FATHER OF NODE I
          i = -ib
          il = il + 1
        end if
  160 continue
      nsteps = is - 1
      return

      end
      subroutine ma27md(n,nz,irn,icn,perm,na,ne,nd,nsteps,lstki,lstkr, &
                       iw,info,ops)
!
! STORAGE AND OPERATION COUNT EVALUATION.
!
! EVALUATE NUMBER OF OPERATIONS AND SPACE REQUIRED BY FACTORIZATION
!     USING MA27B/BD.  THE VALUES GIVEN ARE EXACT ONLY IF NO NUMERICAL
!     PIVOTING IS PERFORMED AND THEN ONLY IF IRN(1) WAS NOT
!     EQUIVALENCED TO IW(1) BY THE USER BEFORE CALLING MA27A/AD.  IF
!     THE EQUIVALENCE HAS BEEN MADE ONLY AN UPPER BOUND FOR NIRNEC
!     AND NRLNEC CAN BE CALCULATED ALTHOUGH THE OTHER COUNTS WILL
!     STILL BE EXACT.
!
! N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
! NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT.  IT IS NOT ALTERED.
! IRN,ICN.  UNLESS IRN(1) HAS BEEN EQUIVALENCED TO IW(1)
!     IRN,ICN MUST BE SET TO THE ROW AND COLUMN INDICES OF THE
!     NON-ZEROS ON INPUT.  THEY ARE NOT ALTERED BY MA27M/MD.
! PERM MUST BE SET TO THE POSITION IN THE PIVOT ORDER OF EACH ROW.
!     IT IS NOT ALTERED.
! NA,NE,ND MUST BE SET TO HOLD, FOR EACH TREE NODE, THE NUMBER OF STACK
!     ELEMENTS ASSEMBLED, THE NUMBER OF ELIMINATIONS AND THE SIZE OF
!     THE ASSEMBLED FRONT MATRIX RESPECTIVELY.  THEY ARE NOT ALTERED.
! NSTEPS MUST BE SET TO HOLD THE NUMBER OF TREE NODES. IT IS NOT
!     ALTERED.
! LSTKI IS USED AS A WORK ARRAY BY MA27M/MD.
! LSTKR.  IF IRN(1) IS EQUIVALENCED TO IW(1)  THEN LSTKR(I)
!     MUST HOLD THE TOTAL NUMBER OF OFF-DIAGONAL ENTRIES (INCLUDING
!     DUPLICATES) IN ROW I (I=1,..,N) OF THE ORIGINAL MATRIX.  IT
!     IS USED AS WORKSPACE BY MA27M/MD.
! IW IS A WORKSPACE ARRAY USED BY OTHER SUBROUTINES AND PASSED TO THIS
!     SUBROUTINE ONLY SO THAT A TEST FOR EQUIVALENCE WITH IRN CAN BE
!     MADE.
!
! COUNTS FOR OPERATIONS AND STORAGE ARE ACCUMULATED IN VARIABLES
!     OPS,NRLTOT,NIRTOT,NRLNEC,NIRNEC,NRLADU,NRLNEC,NIRNEC.
! OPS NUMBER OF MULTIPLICATIONS AND ADDITIONS DURING FACTORIZATION.
! NRLADU,NIRADU REAL AND INTEGER STORAGE RESPECTIVELY FOR THE
!     MATRIX FACTORS.
! NRLTOT,NIRTOT REAL AND INTEGER STRORAGE RESPECTIVELY REQUIRED
!     FOR THE FACTORIZATION IF NO COMPRESSES ARE ALLOWED.
! NRLNEC,NIRNEC REAL AND INTEGER STORAGE RESPECTIVELY REQUIRED FOR
!     THE FACTORIZATION IF COMPRESSES ARE ALLOWED.
! INFO is an INTEGER array of length 20, see MA27A/AD.
! OPS ACCUMULATES THE NO. OF MULTIPLY/ADD PAIRS NEEDED TO CREATE THE
!     TRIANGULAR FACTORIZATION, IN THE DEFINITE CASE.
!
!     .. Scalar Arguments ..
      double precision ops
      integer n,nsteps,nz
!     ..
!     .. Array Arguments ..
      integer icn(*),irn(*),iw(*),lstki(n),lstkr(n),na(nsteps), &
              nd(nsteps),ne(nsteps),perm(n),info(20)
!     ..
!     .. Local Scalars ..
      integer i,inew,iold,iorg,irow,istki,istkr,itop,itree,jold,jorg,k, &
              lstk,nassr,nelim,nfr,nstk,numorg,nz1,nz2
      double precision delim
      integer nrladu,niradu,nirtot,nrltot,nirnec,nrlnec
!     ..
!     .. Intrinsic Functions ..
      intrinsic max,min
!     ..
!     .. Executable Statements ..
!
      if (nz.eq.0) go to 20
! JUMP IF IW AND IRN HAVE NOT BEEN EQUIVALENCED.
      if (irn(1).ne.iw(1)) go to 20
! RESET IRN(1).
      irn(1) = iw(1) - 1
! THE TOTAL NUMBER OF OFF-DIAGONAL ENTRIES IS ACCUMULATED IN NZ2.
! LSTKI IS SET TO HOLD THE TOTAL NUMBER OF ENTRIES (INCUDING
!     THE DIAGONAL) IN EACH ROW IN PERMUTED ORDER.
      nz2 = 0
      do 10 iold = 1,n
        inew = perm(iold)
        lstki(inew) = lstkr(iold) + 1
        nz2 = nz2 + lstkr(iold)
   10 continue
! NZ1 IS THE NUMBER OF ENTRIES IN ONE TRIANGLE INCLUDING THE DIAGONAL.
! NZ2 IS THE TOTAL NUMBER OF ENTRIES INCLUDING THE DIAGONAL.
      nz1 = nz2/2 + n
      nz2 = nz2 + n
      go to 60
! COUNT (IN LSTKI) NON-ZEROS IN ORIGINAL MATRIX BY PERMUTED ROW (UPPER
!     TRIANGLE ONLY). INITIALIZE COUNTS.
   20 do 30 i = 1,n
        lstki(i) = 1
   30 continue
! ACCUMULATE NUMBER OF NON-ZEROS WITH INDICES IN RANGE IN NZ1
!     DUPLICATES ON THE DIAGONAL ARE IGNORED BUT NZ1 INCLUDES ANY
!     DIAGONALS NOT PRESENT ON INPUT.
! ACCUMULATE ROW COUNTS IN LSTKI.
      nz1 = n
      if (nz.eq.0) go to 50
      do 40 i = 1,nz
        iold = irn(i)
        jold = icn(i)
! JUMP IF INDEX IS OUT OF RANGE.
        if (iold.lt.1 .or. iold.gt.n) go to 40
        if (jold.lt.1 .or. jold.gt.n) go to 40
        if (iold.eq.jold) go to 40
        nz1 = nz1 + 1
        irow = min(perm(iold)+0,perm(jold)+0)
        lstki(irow) = lstki(irow) + 1
   40 continue
   50 nz2 = nz1
! ISTKR,ISTKI CURRENT NUMBER OF STACK ENTRIES IN
!     REAL AND INTEGER STORAGE RESPECTIVELY.
! OPS,NRLADU,NIRADU,NIRTOT,NRLTOT,NIRNEC,NRLNEC,NZ2 ARE DEFINED ABOVE.
! NZ2 CURRENT NUMBER OF ORIGINAL MATRIX ENTRIES NOT YET PROCESSED.
! NUMORG CURRENT TOTAL NUMBER OF ROWS ELIMINATED.
! ITOP CURRENT NUMBER OF ELEMENTS ON THE STACK.
   60 istki = 0
      istkr = 0
      ops = 0.0d0
      nrladu = 0
!     ONE LOCATION IS NEEDED TO RECORD THE NUMBER OF BLOCKS
!     ACTUALLY USED.
      niradu = 1
      nirtot = nz1
      nrltot = nz1
      nirnec = nz2
      nrlnec = nz2
      numorg = 0
      itop = 0
!
! EACH PASS THROUGH THIS LOOP PROCESSES A NODE OF THE TREE.
      do 100 itree = 1,nsteps
        nelim = ne(itree)
        delim = nelim
        nfr = nd(itree)
        nstk = na(itree)
! ADJUST STORAGE COUNTS ON ASSEMBLY OF CURRENT FRONTAL MATRIX.
        nassr = nfr* (nfr+1)/2
        if (nstk.ne.0) nassr = nassr - lstkr(itop) + 1
        nrltot = max(nrltot,nrladu+nassr+istkr+nz1)
        nirtot = max(nirtot,niradu+nfr+2+istki+nz1)
        nrlnec = max(nrlnec,nrladu+nassr+istkr+nz2)
        nirnec = max(nirnec,niradu+nfr+2+istki+nz2)
! DECREASE NZ2 BY THE NUMBER OF ENTRIES IN ROWS BEING ELIMINATED AT
!     THIS STAGE.
        do 70 iorg = 1,nelim
          jorg = numorg + iorg
          nz2 = nz2 - lstki(jorg)
   70   continue
        numorg = numorg + nelim
! JUMP IF THERE ARE NO STACK ASSEMBLIES AT THIS NODE.
        if (nstk.le.0) go to 90
! REMOVE ELEMENTS FROM THE STACK.  THERE ARE ITOP ELEMENTS ON THE
!     STACK WITH THE APPROPRIATE ENTRIES IN LSTKR,LSTKI GIVING
!     THE REAL AND INTEGER STORAGE RESPECTIVELY FOR EACH STACK
!     ELEMENT.
        do 80 k = 1,nstk
          lstk = lstkr(itop)
          istkr = istkr - lstk
          lstk = lstki(itop)
          istki = istki - lstk
          itop = itop - 1
   80   continue
! ACCUMULATE NON-ZEROS IN FACTORS AND NUMBER OF OPERATIONS.
   90   nrladu = nrladu + (nelim* (2*nfr-nelim+1))/2
        niradu = niradu + 2 + nfr
        if (nelim.eq.1) niradu = niradu - 1
        ops = ops + ((nfr*delim*(nfr+1)-(2*nfr+1)*delim*(delim+1)/2+ &
              delim* (delim+1)* (2*delim+1)/6)/2)
        if (itree.eq.nsteps) go to 100
! JUMP IF ALL OF FRONTAL MATRIX HAS BEEN ELIMINATED.
        if (nfr.eq.nelim) go to 100
! STACK REMAINDER OF ELEMENT.
        itop = itop + 1
        lstkr(itop) = (nfr-nelim)* (nfr-nelim+1)/2
        lstki(itop) = nfr - nelim + 1
        istki = istki + lstki(itop)
        istkr = istkr + lstkr(itop)
! WE DO NOT NEED TO ADJUST THE COUNTS FOR THE REAL STORAGE BECAUSE
!     THE REMAINDER OF THE FRONTAL MATRIX IS SIMPLY MOVED IN THE
!     STORAGE FROM FACTORS TO STACK AND NO EXTRA STORAGE IS REQUIRED.
        nirtot = max(nirtot,niradu+istki+nz1)
        nirnec = max(nirnec,niradu+istki+nz2)
  100 continue
!
! ADJUST THE STORAGE COUNTS TO ALLOW FOR THE USE OF THE REAL AND
!     INTEGER STORAGE FOR PURPOSES OTHER THAN PURELY THE
!     FACTORIZATION ITSELF.
! THE SECOND TWO TERMS ARE THE MINUMUM AMOUNT REQUIRED BY MA27N/ND.
      nrlnec = max(nrlnec,n+max(nz,nz1))
      nrltot = max(nrltot,n+max(nz,nz1))
      nrlnec = min(nrlnec,nrltot)
      nirnec = max(nz,nirnec)
      nirtot = max(nz,nirtot)
      nirnec = min(nirnec,nirtot)

      info(3) = nrltot
      info(4) = nirtot
      info(5) = nrlnec
      info(6) = nirnec
      info(7) = nrladu
      info(8) = niradu
      return

      end
      subroutine ma27nd(n,nz,nz1,a,la,irn,icn,iw,liw,perm,iw2,icntl, &
                       info)
!
! SORT PRIOR TO FACTORIZATION USING MA27O/OD.
!
! THIS SUBROUTINE REORDERS THE USER'S INPUT SO THAT THE UPPER TRIANGLE
!     OF THE PERMUTED MATRIX, INCLUDING THE DIAGONAL, IS
!     HELD ORDERED BY ROWS AT THE END OF THE STORAGE FOR A AND IW.
!     IT IGNORES ENTRIES WITH ONE OR BOTH INDICES OUT OF RANGE AND
!     ACCUMULATES DIAGONAL ENTRIES.
!     IT ADDS EXPLICIT ZEROS ON THE DIAGONAL WHERE NECESSARY.
! N      - MUST BE SET TO THE ORDER OF THE MATRIX.
!          IT IS NOT ALTERED BY MA27N/ND.
! NZ     - ON ENTRY NZ MUST BE SET TO THE NUMBER
!          OF NON-ZEROS INPUT.  NOT ALTERED BY MA27N/ND.
! NZ1    - ON EXIT NZ1 WILL BE EQUAL TO THE NUMBER OF ENTRIES IN THE
!          SORTED MATRIX.
! A      - ON ENTRY A(I) MUST
!          HOLD THE VALUE OF THE ORIGINAL MATRIX ELEMENT IN POSITION
!          (IRN(I),ICN(I)),I=1,NZ.  ON EXIT A(LA-NZ1+I),I=1,NZ1 HOLDS
!          THE UPPER TRIANGLE OF THE PERMUTED MATRIX BY ROWS WITH
!          THE DIAGONAL ENTRY FIRST ALTHOUGH THERE IS NO FURTHER
!          ORDERING WITHIN THE ROWS THEMSELVES.
! LA     - LENGTH OF ARRAY A. MUST BE AT LEAST N+MAX(NZ,NZ1).
!          IT IS NOT ALTERED BY MA27N/ND.
! IRN    - IRN(I) MUST BE SET TO
!          HOLD THE ROW INDEX OF ENTRY A(I),I=1,NZ.  IRN WILL BE
!          UNALTERED BY MA27N/ND, UNLESS IT IS EQUIVALENCED WITH IW.
! ICN    - ICN(I) MUST BE SET TO
!          HOLD THE COLUMN INDEX OF ENTRY A(I),I=1,NZ.  ICN WILL BE
!          UNALTERED BY MA27N/ND, UNLESS IT IS EQUIVALENCED WITH IW.
! IW     - USED AS WORKSPACE AND ON
!          EXIT, ENTRIES IW(LIW-NZ1+I),I=1,NZ1 HOLD THE COLUMN INDICES
!          (THE ORIGINAL UNPERMUTED INDICES) OF THE CORRESPONDING ENTRY
!          OF A WITH THE FIRST ENTRY FOR EACH ROW FLAGGED NEGATIVE.
!          IRN(1) MAY BE EQUIVALENCED WITH IW(1) AND ICN(1) MAY BE
!          EQUIVALENCED WITH IW(K) WHERE K.GT.NZ.
! LIW    - LENGTH OF ARRAY IW. MUST BE AT LEAST AS
!          GREAT AS THE MAXIMUM OF NZ AND NZ1.
!          NOT ALTERED BY MA27N/ND.
! PERM   - PERM(I) HOLDS THE
!          POSITION IN THE TENTATIVE PIVOT ORDER OF ROW I IN THE
!          ORIGINAL MATRIX,I=1,N. NOT ALTERED BY MA27N/ND.
! IW2    - USED AS WORKSPACE.
!          SEE COMMENTS IN CODE IMMEDIATELY PRIOR TO
!          EACH USE.
! ICNTL is an INTEGER array of length 30, see MA27A/AD.
! INFO is an INTEGER array of length 20, see MA27A/AD.
!   INFO(1)  - ON EXIT FROM MA27N/ND, A ZERO VALUE OF
!          INFO(1) INDICATES THAT NO ERROR HAS BEEN DETECTED.
!          POSSIBLE NON-ZERO VALUES ARE ..
!          +1  WARNING.  INDICES OUT OF RANGE.  THESE ARE IGNORED,
!              THEIR NUMBER IS RECORDED IN INFO(2) OF MA27E/ED AND
!              MESSAGES IDENTIFYING THE FIRST TEN ARE OUTPUT ON UNIT
!              ICNTL(2).
!          -3  INTEGER ARRAY IW IS TOO SMALL.
!          -4  DOUBLE PRECISION ARRAY A IS TOO SMALL.
!
!     .. Parameters ..
      double precision zero
      parameter (zero=0.0d0)
!     ..
!     .. Scalar Arguments ..
      integer la,liw,n,nz,nz1
!     ..
!     .. Array Arguments ..
      double precision a(la)
      integer icn(*),irn(*),iw(liw),iw2(n),perm(n),icntl(30),info(20)
!     ..
!     .. Local Scalars ..
      double precision anext,anow
      integer i,ia,ich,ii,iiw,inew,iold,ipos,j1,j2,jj,jnew,jold,jpos,k
!     ..
!     .. Intrinsic Functions ..
      intrinsic min
!     ..
!     .. Executable Statements ..
      info(1) = 0
! INITIALIZE WORK ARRAY (IW2) IN PREPARATION FOR
!     COUNTING NUMBERS OF NON-ZEROS IN THE ROWS AND INITIALIZE
!     LAST N ENTRIES IN A WHICH WILL HOLD THE DIAGONAL ENTRIES
      ia = la
      do 10 iold = 1,n
        iw2(iold) = 1
        a(ia) = zero
        ia = ia - 1
   10 continue
! SCAN INPUT COPYING ROW INDICES FROM IRN TO THE FIRST NZ POSITIONS
!     IN IW.  THE NEGATIVE OF THE INDEX IS HELD TO FLAG ENTRIES FOR
!     THE IN-PLACE SORT.  ENTRIES IN IW CORRESPONDING TO DIAGONALS AND
!     ENTRIES WITH OUT-OF-RANGE INDICES ARE SET TO ZERO.
!     FOR DIAGONAL ENTRIES, REALS ARE ACCUMULATED IN THE LAST N
!     LOCATIONS OF A.
!     THE NUMBER OF ENTRIES IN EACH ROW OF THE PERMUTED MATRIX IS
!     ACCUMULATED IN IW2.
! INDICES OUT OF RANGE ARE IGNORED  AFTER BEING COUNTED AND
!     AFTER APPROPRIATE MESSAGES HAVE BEEN PRINTED.
      info(2) = 0
! NZ1 IS THE NUMBER OF NON-ZEROS HELD AFTER INDICES OUT OF RANGE HAVE
!     BEEN IGNORED AND DIAGONAL ENTRIES ACCUMULATED.
      nz1 = n
      if (nz.eq.0) go to 80
      do 70 k = 1,nz
        iold = irn(k)
        if (iold.gt.n .or. iold.le.0) go to 30
        jold = icn(k)
        if (jold.gt.n .or. jold.le.0) go to 30
        inew = perm(iold)
        jnew = perm(jold)
        if (inew.ne.jnew) go to 20
        ia = la - n + iold
        a(ia) = a(ia) + a(k)
        go to 60

   20   inew = min(inew,jnew)
! INCREMENT NUMBER OF ENTRIES IN ROW INEW.
        iw2(inew) = iw2(inew) + 1
        iw(k) = -iold
        nz1 = nz1 + 1
        go to 70
! ENTRY OUT OF RANGE.  IT WILL BE IGNORED AND A FLAG SET.
   30   info(1) = 1
        info(2) = info(2) + 1
        if (info(2).le.1 .and. icntl(2).gt.0) then
          write (icntl(2),fmt=40) info(1)
        endif

   40   format (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD', &
                '  *** INFO(1) =',i2)

        if (info(2).le.10 .and. icntl(2).gt.0) then
          write (icntl(2),fmt=50) k,irn(k),icn(k)
        end if

   50   format (i6,'TH NON-ZERO (IN ROW',i6,' AND COLUMN',i6, &
               ') IGNORED')

   60   iw(k) = 0
   70 continue
! CALCULATE POINTERS (IN IW2) TO THE POSITION IMMEDIATELY AFTER THE END
!     OF EACH ROW.
   80 if (nz.lt.nz1 .and. nz1.ne.n) go to 100
! ROOM IS INCLUDED FOR THE DIAGONALS.
      k = 1
      do 90 i = 1,n
        k = k + iw2(i)
        iw2(i) = k
   90 continue
      go to 120
! ROOM IS NOT INCLUDED FOR THE DIAGONALS.
  100 k = 1
      do 110 i = 1,n
        k = k + iw2(i) - 1
        iw2(i) = k
  110 continue
! FAIL IF INSUFFICIENT SPACE IN ARRAYS A OR IW.
  120 if (nz1.gt.liw) go to 210
      if (nz1+n.gt.la) go to 220
! NOW RUN THROUGH NON-ZEROS IN ORDER PLACING THEM IN THEIR NEW
! POSITION AND DECREMENTING APPROPRIATE IW2 ENTRY.  IF WE ARE
! ABOUT TO OVERWRITE AN ENTRY NOT YET MOVED, WE MUST DEAL WITH
! THIS AT THIS TIME.
      if (nz1.eq.n) go to 180
      do 140 k = 1,nz
        iold = -iw(k)
        if (iold.le.0) go to 140
        jold = icn(k)
        anow = a(k)
        iw(k) = 0
        do 130 ich = 1,nz
          inew = perm(iold)
          jnew = perm(jold)
          inew = min(inew,jnew)
          if (inew.eq.perm(jold)) jold = iold
          jpos = iw2(inew) - 1
          iold = -iw(jpos)
          anext = a(jpos)
          a(jpos) = anow
          iw(jpos) = jold
          iw2(inew) = jpos
          if (iold.eq.0) go to 140
          anow = anext
          jold = icn(jpos)
  130   continue
  140 continue
      if (nz.ge.nz1) go to 180
! MOVE UP ENTRIES TO ALLOW FOR DIAGONALS.
      ipos = nz1
      jpos = nz1 - n
      do 170 ii = 1,n
        i = n - ii + 1
        j1 = iw2(i)
        j2 = jpos
        if (j1.gt.jpos) go to 160
        do 150 jj = j1,j2
          iw(ipos) = iw(jpos)
          a(ipos) = a(jpos)
          ipos = ipos - 1
          jpos = jpos - 1
  150   continue
  160   iw2(i) = ipos + 1
        ipos = ipos - 1
  170 continue
! RUN THROUGH ROWS INSERTING DIAGONAL ENTRIES AND FLAGGING BEGINNING
!     OF EACH ROW BY NEGATING FIRST COLUMN INDEX.
  180 do 190 iold = 1,n
        inew = perm(iold)
        jpos = iw2(inew) - 1
        ia = la - n + iold
        a(jpos) = a(ia)
        iw(jpos) = -iold
  190 continue
! MOVE SORTED MATRIX TO THE END OF THE ARRAYS.
      ipos = nz1
      ia = la
      iiw = liw
      do 200 i = 1,nz1
        a(ia) = a(ipos)
        iw(iiw) = iw(ipos)
        ipos = ipos - 1
        ia = ia - 1
        iiw = iiw - 1
  200 continue
      go to 230
! **** ERROR RETURN ****
  210 info(1) = -3
      info(2) = nz1
      go to 230

  220 info(1) = -4
      info(2) = nz1 + n
!
  230 return

      end
      subroutine ma27od(n,nz,a,la,iw,liw,perm,nstk,nsteps,maxfrt,nelim, &
                       iw2,icntl,cntl,info)
!
! FACTORIZATION SUBROUTINE
!
! THIS SUBROUTINE OPERATES ON THE INPUT MATRIX ORDERED BY MA27N/ND AND
!     PRODUCES THE FACTORS OF U AND D ('A'=UTRANSPOSE*D*U) FOR USE IN
!     SUBSEQUENT SOLUTIONS.  GAUSSIAN ELIMINATION IS USED WITH PIVOTS
!     CHOSEN FROM THE DIAGONAL.  TO ENSURE STABILITY, BLOCK PIVOTS OF
!     ORDER 2 WILL BE USED IF THE DIAGONAL ENTRY IS NOT LARGE ENOUGH.
!
! N      - MUST BE SET TO THE ORDER OF THE MATRIX. IT IS NOT ALTERED.
! NZ     - MUST BE SET TO THE NUMBER OF NON-ZEROS IN UPPER TRIANGLE OF
!          PERMUTED MATRIX.  NOT ALTERED BY MA27O/OD.
! A      - MUST BE SET ON INPUT TO MATRIX HELD BY ROWS REORDERED BY
!          PERMUTATION FROM MA27A/AD IN A(LA-NZ+I),I=1,NZ.   ON
!          EXIT FROM MA27O/OD, THE FACTORS OF U AND D ARE HELD IN
!          POSITIONS 1 TO POSFAC-1.
! LA     - LENGTH OF ARRAY A.  A VALUE FOR LA
!          SUFFICIENT FOR DEFINITE SYSTEMS
!          WILL HAVE BEEN PROVIDED BY MA27A/AD. NOT ALTERED BY MA27O/OD.
! IW     - MUST BE SET SO THAT,ON INPUT, IW(LIW-NZ+I),I=1,NZ
!          HOLDS THE COLUMN INDEX OF THE ENTRY IN A(LA-NZ+I).  ON EXIT,
!          IW HOLDS INTEGER INDEXING INFORMATION ON THE FACTORS.
!          THE ABSOLUTE VALUE OF THE FIRST ENTRY IN IW WILL BE SET TO
!          THE NUMBER OF BLOCK PIVOTS ACTUALLY USED.  THIS MAY BE
!          DIFFERENT FROM NSTEPS SINCE NUMERICAL CONSIDERATIONS
!          MAY PREVENT US CHOOSING A PIVOT AT EACH STAGE.  IF THIS ENTRY
!          IN IW IS NEGATIVE, THEN AT LEAST ONE TWO BY TWO
!          PIVOT HAS BEEN USED DURING THE DECOMPOSITION.
!          INTEGER INFORMATION ON EACH BLOCK PIVOT ROW FOLLOWS.  FOR
!          EACH BLOCK PIVOT ROW THE COLUMN INDICES ARE PRECEDED BY A
!          COUNT OF THE NUMBER OF ROWS AND COLUMNS IN THE BLOCK PIVOT
!          WHERE, IF ONLY ONE ROW IS PRESENT, ONLY THE NUMBER OF
!          COLUMNS TOGETHER WITH A NEGATIVE FLAG IS HELD.  THE FIRST
!          COLUMN INDEX FOR A TWO BY TWO PIVOT IS FLAGGED NEGATIVE.
! LIW    - LENGTH OF ARRAY IW.  A VALUE FOR LIW SUFFICIENT FOR
!          DEFINITE SYSTEMS
!          WILL HAVE BEEN PROVIDED BY MA27A/AD.  NOT ALTERED BY MA27O/OD
! PERM   - PERM(I) MUST BE SET TO POSITION OF ROW/COLUMN I IN THE
!          TENTATIVE PIVOT ORDER GENERATED BY MA27A/AD.
!          IT IS NOT ALTERED BY MA27O/OD.
! NSTK   - MUST BE LEFT UNCHANGED SINCE OUTPUT FROM MA27A/AD. NSTK(I)
!          GIVES THE NUMBER OF GENERATED STACK ELEMENTS ASSEMBLED AT
!          STAGE I.  IT IS NOT ALTERED BY MA27O/OD.
! NSTEPS - LENGTH OF ARRAYS NSTK AND NELIM. VALUE IS GIVEN ON OUTPUT
!          FROM MA27A/AD (WILL NEVER EXCEED N). IT IS NOT ALTERED BY
!          MA27O/OD.
! MAXFRT - NEED NOT BE SET ON INPUT.  ON OUTPUT
!          MAXFRT WILL BE SET TO THE MAXIMUM FRONT SIZE ENCOUNTERED
!          DURING THE DECOMPOSITION.
! NELIM  - MUST BE UNCHANGED SINCE OUTPUT FROM MA27A/AD. NELIM(I)
!          GIVES THE NUMBER OF ORIGINAL ROWS ASSEMBLED AT STAGE I.
!          IT IS NOT ALTERED BY MA27O/OD.
! IW2    - INTEGER ARRAY OF LENGTH N. USED AS WORKSPACE BY MA27O/OD.
!          ALTHOUGH WE COULD HAVE USED A SHORT WORD INTEGER IN THE IBM
!          VERSION, WE HAVE NOT DONE SO BECAUSE THERE IS A SPARE
!          FULL INTEGER ARRAY (USED IN THE SORT BEFORE MA27O/OD)
!          AVAILABLE WHEN MA27O/OD IS CALLED FROM MA27B/BD.
! ICNTL is an INTEGER array of length 30, see MA27A/AD.
! CNTL is a DOUBLE PRECISION array of length 5, see MA27A/AD.
! INFO is an INTEGER array of length 20, see MA27A/AD.
!   INFO(1)  - INTEGER VARIABLE.  DIAGNOSTIC FLAG.  A ZERO VALUE ON EXIT
!          INDICATES SUCCESS.  POSSIBLE NEGATIVE VALUES ARE ...
!          -3  INSUFFICIENT STORAGE FOR IW.
!          -4  INSUFFICIENT STORAGE FOR A.
!          -5  ZERO PIVOT FOUND IN FACTORIZATION OF DEFINITE MATRIX.
!
!     .. Parameters ..
      double precision zero,half,one
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0)
!     ..
!     .. Scalar Arguments ..
      integer la,liw,maxfrt,n,nsteps,nz
!     ..
!     .. Array Arguments ..
      double precision a(la),cntl(5)
      integer iw(liw),iw2(n),nelim(nsteps),nstk(nsteps),perm(n)
      integer icntl(30),info(20)
!     ..
!     .. Local Scalars ..
      double precision amax,amult,amult1,amult2,detpiv,rmax,swop, &
              thresh,tmax,uu
      integer ainput,apos,apos1,apos2,apos3,astk,astk2,azero,i,iass, &
              ibeg,idummy,iell,iend,iexch,ifr,iinput,ioldps,iorg,ipiv, &
              ipmnp,ipos,irow,isnpiv,istk,istk2,iswop,iwpos,ix,iy,j,j1, &
              j2,jcol,jdummy,jfirst,jj,jj1,jjj,jlast,jmax,jmxmip,jnew, &
              jnext,jpiv,jpos,k,k1,k2,kdummy,kk,kmax,krow,laell,lapos2, &
              liell,lnass,lnpiv,lt,ltopst,nass,nblk,newel,nfront,npiv, &
              npivp1,ntotpv,numass,numorg,numstk,pivsiz,posfac,pospv1, &
              pospv2
      integer ntwo,neig,ncmpbi,ncmpbr,nrlbdu,nirbdu
!     ..
!     .. External Subroutines ..
      external ma27pd
!     ..
!     .. Intrinsic Functions ..
      intrinsic abs,max,min
!     ..
!     .. Statement Functions ..
      integer idiag
!     ..
!     .. Statement Function definitions ..
! THE FOLLOWING ARITHMETIC FUNCTION GIVES THE DISPLACEMENT FROM
!     THE START OF THE ASSEMBLED MATRIX(OF ORDER IX) OF THE DIAGONAL
!     ENTRY IN ITS ROW IY.
      idiag(ix,iy) = ((iy-1)* (2*ix-iy+2))/2
!     ..
!     .. Executable Statements ..
! INITIALIZATION.
! NBLK IS THE NUMBER OF BLOCK PIVOTS USED.
      nblk = 0
      ntwo = 0
      neig = 0
      ncmpbi = 0
      ncmpbr = 0
      maxfrt = 0
      nrlbdu = 0
      nirbdu = 0
! A PRIVATE VARIABLE UU IS SET TO CNTL(1), SO THAT CNTL(1) WILL REMAIN
! UNALTERED.
      uu = min(cntl(1),half)
      uu = max(uu,-half)
      do 10 i = 1,n
        iw2(i) = 0
   10 continue
! IWPOS IS POINTER TO FIRST FREE POSITION FOR FACTORS IN IW.
! POSFAC IS POINTER FOR FACTORS IN A. AT EACH PASS THROUGH THE
!     MAJOR LOOP POSFAC INITIALLY POINTS TO THE FIRST FREE LOCATION
!     IN A AND THEN IS SET TO THE POSITION OF THE CURRENT PIVOT IN A.
! ISTK IS POINTER TO TOP OF STACK IN IW.
! ISTK2 IS POINTER TO BOTTOM OF STACK IN IW (NEEDED BY COMPRESS).
! ASTK IS POINTER TO TOP OF STACK IN A.
! ASTK2 IS POINTER TO BOTTOM OF STACK IN A (NEEDED BY COMPRESS).
! IINPUT IS POINTER TO CURRENT POSITION IN ORIGINAL ROWS IN IW.
! AINPUT IS POINTER TO CURRENT POSITION IN ORIGINAL ROWS IN A.
! AZERO IS POINTER TO LAST POSITION ZEROED IN A.
! NTOTPV IS THE TOTAL NUMBER OF PIVOTS SELECTED. THIS IS USED
!     TO DETERMINE WHETHER THE MATRIX IS SINGULAR.
      iwpos = 2
      posfac = 1
      istk = liw - nz + 1
      istk2 = istk - 1
      astk = la - nz + 1
      astk2 = astk - 1
      iinput = istk
      ainput = astk
      azero = 0
      ntotpv = 0
! NUMASS IS THE ACCUMULATED NUMBER OF ROWS ASSEMBLED SO FAR.
      numass = 0
!
! EACH PASS THROUGH THIS MAIN LOOP PERFORMS ALL THE OPERATIONS
!     ASSOCIATED WITH ONE SET OF ASSEMBLY/ELIMINATIONS.
      do 760 iass = 1,nsteps
! NASS WILL BE SET TO THE NUMBER OF FULLY ASSEMBLED VARIABLES IN
!     CURRENT NEWLY CREATED ELEMENT.
        nass = nelim(iass)
! NEWEL IS A POINTER INTO IW TO CONTROL OUTPUT OF INTEGER INFORMATION
!     FOR NEWLY CREATED ELEMENT.
        newel = iwpos + 1
! SYMBOLICALLY ASSEMBLE INCOMING ROWS AND GENERATED STACK ELEMENTS
!     ORDERING THE RESULTANT ELEMENT ACCORDING TO PERMUTATION PERM.  WE
!     ASSEMBLE THE STACK ELEMENTS FIRST BECAUSE THESE WILL ALREADY BE
!     ORDERED.
! SET HEADER POINTER FOR MERGE OF INDEX LISTS.
        jfirst = n + 1
! INITIALIZE NUMBER OF VARIABLES IN CURRENT FRONT.
        nfront = 0
        numstk = nstk(iass)
        ltopst = 1
        lnass = 0
! JUMP IF NO STACK ELEMENTS ARE BEING ASSEMBLED AT THIS STAGE.
        if (numstk.eq.0) go to 80
        j2 = istk - 1
        lnass = nass
        ltopst = ((iw(istk)+1)*iw(istk))/2
        do 70 iell = 1,numstk
! ASSEMBLE ELEMENT IELL PLACING
!     THE INDICES INTO A LINKED LIST IN IW2 ORDERED
!     ACCORDING TO PERM.
          jnext = jfirst
          jlast = n + 1
          j1 = j2 + 2
          j2 = j1 - 1 + iw(j1-1)
! RUN THROUGH INDEX LIST OF STACK ELEMENT IELL.
          do 60 jj = j1,j2
            j = iw(jj)
! JUMP IF ALREADY ASSEMBLED
            if (iw2(j).gt.0) go to 60
            jnew = perm(j)
! IF VARIABLE WAS PREVIOUSLY FULLY SUMMED BUT WAS NOT PIVOTED ON
!     EARLIER BECAUSE OF NUMERICAL TEST, INCREMENT NUMBER OF FULLY
!     SUMMED ROWS/COLUMNS IN FRONT.
            if (jnew.le.numass) nass = nass + 1
! FIND POSITION IN LINKED LIST FOR NEW VARIABLE.  NOTE THAT WE START
!     FROM WHERE WE LEFT OFF AFTER ASSEMBLY OF PREVIOUS VARIABLE.
            do 20 idummy = 1,n
              if (jnext.eq.n+1) go to 30
              if (perm(jnext).gt.jnew) go to 30
              jlast = jnext
              jnext = iw2(jlast)
   20       continue
   30       if (jlast.ne.n+1) go to 40
            jfirst = j
            go to 50

   40       iw2(jlast) = j
   50       iw2(j) = jnext
            jlast = j
! INCREMENT NUMBER OF VARIABLES IN THE FRONT.
            nfront = nfront + 1
   60     continue
   70   continue
        lnass = nass - lnass
! NOW INCORPORATE ORIGINAL ROWS.  NOTE THAT THE COLUMNS IN THESE ROWS
!     NEED NOT BE IN ORDER. WE ALSO PERFORM
!     A SWOP SO THAT THE DIAGONAL ENTRY IS THE FIRST IN ITS
!     ROW.  THIS ALLOWS US TO AVOID STORING THE INVERSE OF ARRAY PERM.
   80   numorg = nelim(iass)
        j1 = iinput
        do 150 iorg = 1,numorg
          j = -iw(j1)
          do 140 idummy = 1,liw
            jnew = perm(j)
! JUMP IF VARIABLE ALREADY INCLUDED.
            if (iw2(j).gt.0) go to 130
! HERE WE MUST ALWAYS START OUR SEARCH AT THE BEGINNING.
            jlast = n + 1
            jnext = jfirst
            do 90 jdummy = 1,n
              if (jnext.eq.n+1) go to 100
              if (perm(jnext).gt.jnew) go to 100
              jlast = jnext
              jnext = iw2(jlast)
   90       continue
  100       if (jlast.ne.n+1) go to 110
            jfirst = j
            go to 120

  110       iw2(jlast) = j
  120       iw2(j) = jnext
! INCREMENT NUMBER OF VARIABLES IN FRONT.
            nfront = nfront + 1
  130       j1 = j1 + 1
            if (j1.gt.liw) go to 150
            j = iw(j1)
            if (j.lt.0) go to 150
  140     continue
  150   continue
! NOW RUN THROUGH LINKED LIST IW2 PUTTING INDICES OF VARIABLES IN NEW
!     ELEMENT INTO IW AND SETTING IW2 ENTRY TO POINT TO THE RELATIVE
!     POSITION OF THE VARIABLE IN THE NEW ELEMENT.
        if (newel+nfront.lt.istk) go to 160
! COMPRESS IW.
        call ma27pd(a,iw,istk,istk2,iinput,2,ncmpbr,ncmpbi)
        if (newel+nfront.lt.istk) go to 160
        info(2) = liw + 1 + newel + nfront - istk
        go to 770

  160   j = jfirst
        do 170 ifr = 1,nfront
          newel = newel + 1
          iw(newel) = j
          jnext = iw2(j)
          iw2(j) = newel - (iwpos+1)
          j = jnext
  170   continue
!
! ASSEMBLE REALS INTO FRONTAL MATRIX.
        maxfrt = max(maxfrt,nfront)
        iw(iwpos) = nfront
! FIRST ZERO OUT FRONTAL MATRIX AS APPROPRIATE FIRST CHECKING TO SEE
!     IF THERE IS SUFFICIENT SPACE.
        laell = ((nfront+1)*nfront)/2
        apos2 = posfac + laell - 1
        if (numstk.ne.0) lnass = lnass* (2*nfront-lnass+1)/2
        if (posfac+lnass-1.ge.astk) go to 180
        if (apos2.lt.astk+ltopst-1) go to 190
! COMPRESS A.
  180   call ma27pd(a,iw,astk,astk2,ainput,1,ncmpbr,ncmpbi)
        if (posfac+lnass-1.ge.astk) go to 780
        if (apos2.ge.astk+ltopst-1) go to 780
  190   if (apos2.le.azero) go to 220
        apos = azero + 1
        lapos2 = min(apos2,astk-1)
        if (lapos2.lt.apos) go to 210
        do 200 k = apos,lapos2
          a(k) = zero
  200   continue
  210   azero = apos2
! JUMP IF THERE ARE NO STACK ELEMENTS TO ASSEMBLE.
  220   if (numstk.eq.0) go to 260
! PLACE REALS CORRESPONDING TO STACK ELEMENTS IN CORRECT POSITIONS IN A.
        do 250 iell = 1,numstk
          j1 = istk + 1
          j2 = istk + iw(istk)
          do 240 jj = j1,j2
            irow = iw(jj)
            irow = iw2(irow)
            apos = posfac + idiag(nfront,irow)
            do 230 jjj = jj,j2
              j = iw(jjj)
              apos2 = apos + iw2(j) - irow
              a(apos2) = a(apos2) + a(astk)
              a(astk) = zero
              astk = astk + 1
  230       continue
  240     continue
! INCREMENT STACK POINTER.
          istk = j2 + 1
  250   continue
! INCORPORATE REALS FROM ORIGINAL ROWS.
  260   do 280 iorg = 1,numorg
          j = -iw(iinput)
! WE CAN DO THIS BECAUSE THE DIAGONAL IS NOW THE FIRST ENTRY.
          irow = iw2(j)
          apos = posfac + idiag(nfront,irow)
! THE FOLLOWING LOOP GOES FROM 1 TO NZ BECAUSE THERE MAY BE DUPLICATES.
          do 270 idummy = 1,nz
            apos2 = apos + iw2(j) - irow
            a(apos2) = a(apos2) + a(ainput)
            ainput = ainput + 1
            iinput = iinput + 1
            if (iinput.gt.liw) go to 280
            j = iw(iinput)
            if (j.lt.0) go to 280
  270     continue
  280   continue
! RESET IW2 AND NUMASS.
        numass = numass + numorg
        j1 = iwpos + 2
        j2 = iwpos + nfront + 1
        do 290 k = j1,j2
          j = iw(k)
          iw2(j) = 0
  290   continue
! PERFORM PIVOTING ON ASSEMBLED ELEMENT.
! NPIV IS THE NUMBER OF PIVOTS SO FAR SELECTED.
! LNPIV IS THE NUMBER OF PIVOTS SELECTED AFTER THE LAST PASS THROUGH
!     THE THE FOLLOWING LOOP.
        lnpiv = -1
        npiv = 0
        do 650 kdummy = 1,nass
          if (npiv.eq.nass) go to 660
          if (npiv.eq.lnpiv) go to 660
          lnpiv = npiv
          npivp1 = npiv + 1
! JPIV IS USED AS A FLAG TO INDICATE WHEN 2 BY 2 PIVOTING HAS OCCURRED
!     SO THAT IPIV IS INCREMENTED CORRECTLY.
          jpiv = 1
! NASS IS MAXIMUM POSSIBLE NUMBER OF PIVOTS.
! WE EITHER TAKE THE DIAGONAL ENTRY OR THE 2 BY 2 PIVOT WITH THE
!     LARGEST OFF-DIAGONAL AT EACH STAGE.
! EACH PASS THROUGH THIS LOOP TRIES TO CHOOSE ONE PIVOT.
          do 640 ipiv = npivp1,nass
            jpiv = jpiv - 1
! JUMP IF WE HAVE JUST PROCESSED A 2 BY 2 PIVOT.
            if (jpiv.eq.1) go to 640
            apos = posfac + idiag(nfront-npiv,ipiv-npiv)
! IF THE USER HAS INDICATED THAT THE MATRIX IS DEFINITE, WE
!     DO NOT NEED TO TEST FOR STABILITY BUT WE DO CHECK TO SEE IF THE
!     PIVOT IS NON-ZERO OR HAS CHANGED SIGN.
!     IF IT IS ZERO, WE EXIT WITH AN ERROR. IF IT HAS CHANGED SIGN
!     AND U WAS SET NEGATIVE, THEN WE AGAIN EXIT IMMEDIATELY.  IF THE
!     PIVOT CHANGES SIGN AND U WAS ZERO, WE CONTINUE WITH THE
!     FACTORIZATION BUT PRINT A WARNING MESSAGE ON UNIT ICNTL(2).
! ISNPIV HOLDS A FLAG FOR THE SIGN OF THE PIVOTS TO DATE SO THAT
!     A SIGN CHANGE WHEN DECOMPOSING AN ALLEGEDLY DEFINITE MATRIX CAN
!     BE DETECTED.
            if (uu.gt.zero) go to 320
            if (abs(a(apos)).le.cntl(3)) go to 790
! JUMP IF THIS IS NOT THE FIRST PIVOT TO BE SELECTED.
            if (ntotpv.gt.0) go to 300
! SET ISNPIV.
            if (a(apos).gt.zero) isnpiv = 1
            if (a(apos).lt.zero) isnpiv = -1
  300       if (a(apos).gt.zero .and. isnpiv.eq.1) go to 560
            if (a(apos).lt.zero .and. isnpiv.eq.-1) go to 560
            if (info(1).ne.2) info(2) = 0
            info(2) = info(2) + 1
            info(1) = 2
            i = ntotpv + 1
            if (icntl(2).gt.0 .and. info(2).le.10) then
              write (icntl(2),fmt=310) info(1),i
            end if

  310       format (' *** WARNING MESSAGE FROM SUBROUTINE MA27BD', &
                    '  *** INFO(1) =',i2,/,' PIVOT',i6, &
                   ' HAS DIFFERENT SIGN FROM THE PREVIOUS ONE')

            isnpiv = -isnpiv
            if (uu.eq.zero) go to 560
            go to 800

  320       amax = zero
            tmax = amax
! FIND LARGEST ENTRY TO RIGHT OF DIAGONAL IN ROW OF PROSPECTIVE PIVOT
!     IN THE FULLY-SUMMED PART.  ALSO RECORD COLUMN OF THIS LARGEST
!     ENTRY.
            j1 = apos + 1
            j2 = apos + nass - ipiv
            if (j2.lt.j1) go to 340
            do 330 jj = j1,j2
              if (abs(a(jj)).le.amax) go to 330
              jmax = ipiv + jj - j1 + 1
              amax = abs(a(jj))
  330       continue
! DO SAME AS ABOVE FOR NON-FULLY-SUMMED PART ONLY HERE WE DO NOT NEED
!     TO RECORD COLUMN SO LOOP IS SIMPLER.
  340       j1 = j2 + 1
            j2 = apos + nfront - ipiv
            if (j2.lt.j1) go to 360
            do 350 jj = j1,j2
              tmax = max(abs(a(jj)),tmax)
  350       continue
! NOW CALCULATE LARGEST ENTRY IN OTHER PART OF ROW.
  360       rmax = max(tmax,amax)
            apos1 = apos
            kk = nfront - ipiv
            lt = ipiv - (npiv+1)
            if (lt.eq.0) go to 380
            do 370 k = 1,lt
              kk = kk + 1
              apos1 = apos1 - kk
              rmax = max(rmax,abs(a(apos1)))
  370       continue
! JUMP IF STABILITY TEST SATISFIED.
  380       if (abs(a(apos)).gt.max(cntl(3),uu*rmax)) go to 450
! CHECK BLOCK PIVOT OF ORDER 2 FOR STABILITY.
            if (abs(amax).le.cntl(3)) go to 640
            apos2 = posfac + idiag(nfront-npiv,jmax-npiv)
            detpiv = a(apos)*a(apos2) - amax*amax
            thresh = abs(detpiv)
! SET THRESH TO U TIMES THE RECIPROCAL OF THE MAX-NORM OF THE INVERSE
!     OF THE PROSPECTIVE BLOCK.
            thresh = thresh/ (uu*max(abs(a(apos))+amax, &
                     abs(a(apos2))+amax))
! CHECK 2 BY 2 PIVOT FOR STABILITY.
! FIRST CHECK AGAINST ROW IPIV.
            if (thresh.le.rmax) go to 640
! FIND LARGEST ENTRY IN ROW JMAX.
! FIND MAXIMUM TO THE RIGHT OF THE DIAGONAL.
            rmax = zero
            j1 = apos2 + 1
            j2 = apos2 + nfront - jmax
            if (j2.lt.j1) go to 400
            do 390 jj = j1,j2
              rmax = max(rmax,abs(a(jj)))
  390       continue
! NOW CHECK TO THE LEFT OF THE DIAGONAL.
! WE USE TWO LOOPS TO AVOID TESTING FOR ROW IPIV INSIDE THE LOOP.
  400       kk = nfront - jmax + 1
            apos3 = apos2
            jmxmip = jmax - ipiv - 1
            if (jmxmip.eq.0) go to 420
            do 410 k = 1,jmxmip
              apos2 = apos2 - kk
              kk = kk + 1
              rmax = max(rmax,abs(a(apos2)))
  410       continue
  420       ipmnp = ipiv - npiv - 1
            if (ipmnp.eq.0) go to 440
            apos2 = apos2 - kk
            kk = kk + 1
            do 430 k = 1,ipmnp
              apos2 = apos2 - kk
              kk = kk + 1
              rmax = max(rmax,abs(a(apos2)))
  430       continue
  440       if (thresh.le.rmax) go to 640
            pivsiz = 2
            go to 460

  450       pivsiz = 1
  460       irow = ipiv - npiv
!
! PIVOT HAS BEEN CHOSEN.  IF BLOCK PIVOT OF ORDER 2, PIVSIZ IS EQUAL TO
!     TWO OTHERWISE PIVSIZ EQUALS ONE..
! THE FOLLOWING LOOP MOVES THE PIVOT BLOCK TO THE TOP LEFT HAND CORNER
!     OF THE FRONTAL MATRIX.
            do 550 krow = 1,pivsiz
! WE JUMP IF SWOP IS NOT NECESSARY.
              if (irow.eq.1) go to 530
              j1 = posfac + irow
              j2 = posfac + nfront - (npiv+1)
              if (j2.lt.j1) go to 480
              apos2 = apos + 1
! SWOP PORTION OF ROWS WHOSE COLUMN INDICES ARE GREATER THAN LATER ROW.
              do 470 jj = j1,j2
                swop = a(apos2)
                a(apos2) = a(jj)
                a(jj) = swop
                apos2 = apos2 + 1
  470         continue
  480         j1 = posfac + 1
              j2 = posfac + irow - 2
              apos2 = apos
              kk = nfront - (irow+npiv)
              if (j2.lt.j1) go to 500
! SWOP PORTION OF ROWS/COLUMNS WHOSE INDICES LIE BETWEEN THE TWO ROWS.
              do 490 jjj = j1,j2
                jj = j2 - jjj + j1
                kk = kk + 1
                apos2 = apos2 - kk
                swop = a(apos2)
                a(apos2) = a(jj)
                a(jj) = swop
  490         continue
  500         if (npiv.eq.0) go to 520
              apos1 = posfac
              kk = kk + 1
              apos2 = apos2 - kk
! SWOP PORTION OF COLUMNS WHOSE INDICES ARE LESS THAN EARLIER ROW.
              do 510 jj = 1,npiv
                kk = kk + 1
                apos1 = apos1 - kk
                apos2 = apos2 - kk
                swop = a(apos2)
                a(apos2) = a(apos1)
                a(apos1) = swop
  510         continue
! SWOP DIAGONALS AND INTEGER INDEXING INFORMATION
  520         swop = a(apos)
              a(apos) = a(posfac)
              a(posfac) = swop
              ipos = iwpos + npiv + 2
              iexch = iwpos + irow + npiv + 1
              iswop = iw(ipos)
              iw(ipos) = iw(iexch)
              iw(iexch) = iswop
  530         if (pivsiz.eq.1) go to 550
! SET VARIABLES FOR THE SWOP OF SECOND ROW OF BLOCK PIVOT.
              if (krow.eq.2) go to 540
              irow = jmax - (npiv+1)
              jpos = posfac
              posfac = posfac + nfront - npiv
              npiv = npiv + 1
              apos = apos3
              go to 550
! RESET VARIABLES PREVIOUSLY SET FOR SECOND PASS.
  540         npiv = npiv - 1
              posfac = jpos
  550       continue
!
            if (pivsiz.eq.2) go to 600
! PERFORM THE ELIMINATION USING ENTRY (IPIV,IPIV) AS PIVOT.
! WE STORE U AND DINVERSE.
  560       a(posfac) = one/a(posfac)
            if (a(posfac).lt.zero) neig = neig + 1
            j1 = posfac + 1
            j2 = posfac + nfront - (npiv+1)
            if (j2.lt.j1) go to 590
            ibeg = j2 + 1
            do 580 jj = j1,j2
              amult = -a(jj)*a(posfac)
              iend = ibeg + nfront - (npiv+jj-j1+2)
! THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
!DIR$ IVDEP
              do 570 irow = ibeg,iend
                jcol = jj + irow - ibeg
                a(irow) = a(irow) + amult*a(jcol)
  570         continue
              ibeg = iend + 1
              a(jj) = amult
  580       continue
  590       npiv = npiv + 1
            ntotpv = ntotpv + 1
            jpiv = 1
            posfac = posfac + nfront - npiv + 1
            go to 640
! PERFORM ELIMINATION USING BLOCK PIVOT OF ORDER TWO.
! REPLACE BLOCK PIVOT BY ITS INVERSE.
! SET FLAG TO INDICATE USE OF 2 BY 2 PIVOT IN IW.
  600       ipos = iwpos + npiv + 2
            ntwo = ntwo + 1
            iw(ipos) = -iw(ipos)
            pospv1 = posfac
            pospv2 = posfac + nfront - npiv
            swop = a(pospv2)
            if (detpiv.lt.zero) neig = neig + 1
            if (detpiv.gt.zero .and. swop.lt.zero) neig = neig + 2
            a(pospv2) = a(pospv1)/detpiv
            a(pospv1) = swop/detpiv
            a(pospv1+1) = -a(pospv1+1)/detpiv
            j1 = pospv1 + 2
            j2 = pospv1 + nfront - (npiv+1)
            if (j2.lt.j1) go to 630
            jj1 = pospv2
            ibeg = pospv2 + nfront - (npiv+1)
            do 620 jj = j1,j2
              jj1 = jj1 + 1
              amult1 = - (a(pospv1)*a(jj)+a(pospv1+1)*a(jj1))
              amult2 = - (a(pospv1+1)*a(jj)+a(pospv2)*a(jj1))
              iend = ibeg + nfront - (npiv+jj-j1+3)
! THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
!DIR$ IVDEP
              do 610 irow = ibeg,iend
                k1 = jj + irow - ibeg
                k2 = jj1 + irow - ibeg
                a(irow) = a(irow) + amult1*a(k1) + amult2*a(k2)
  610         continue
              ibeg = iend + 1
              a(jj) = amult1
              a(jj1) = amult2
  620       continue
  630       npiv = npiv + 2
            ntotpv = ntotpv + 2
            jpiv = 2
            posfac = pospv2 + nfront - npiv + 1
  640     continue
  650   continue
! END OF MAIN ELIMINATION LOOP.
!
  660   if (npiv.ne.0) nblk = nblk + 1
        ioldps = iwpos
        iwpos = iwpos + nfront + 2
        if (npiv.eq.0) go to 690
        if (npiv.gt.1) go to 680
        iw(ioldps) = -iw(ioldps)
        do 670 k = 1,nfront
          j1 = ioldps + k
          iw(j1) = iw(j1+1)
  670   continue
        iwpos = iwpos - 1
        go to 690

  680   iw(ioldps+1) = npiv
! COPY REMAINDER OF ELEMENT TO TOP OF STACK
  690   liell = nfront - npiv
        if (liell.eq.0 .or. iass.eq.nsteps) go to 750
        if (iwpos+liell.lt.istk) go to 700
        call ma27pd(a,iw,istk,istk2,iinput,2,ncmpbr,ncmpbi)
  700   istk = istk - liell - 1
        iw(istk) = liell
        j1 = istk
        kk = iwpos - liell - 1
! THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
!DIR$ IVDEP
        do 710 k = 1,liell
          j1 = j1 + 1
          kk = kk + 1
          iw(j1) = iw(kk)
  710   continue
! WE COPY IN REVERSE DIRECTION TO AVOID OVERWRITE PROBLEMS.
        laell = ((liell+1)*liell)/2
        kk = posfac + laell
        if (kk.ne.astk) go to 720
        astk = astk - laell
        go to 740
! THE MOVE AND ZEROING OF ARRAY A IS PERFORMED WITH TWO LOOPS SO
! THAT THE CRAY-1 WILL VECTORIZE THEM SAFELY.
  720   kmax = kk - 1
! THE FOLLOWING SPECIAL COMMENT FORCES VECTORIZATION ON THE CRAY-1.
!DIR$ IVDEP
        do 730 k = 1,laell
          kk = kk - 1
          astk = astk - 1
          a(astk) = a(kk)
  730   continue
        kmax = min(kmax,astk-1)
        do 735 k = kk,kmax
          a(k) = zero
  735   continue
  740   azero = min(azero,astk-1)
  750   if (npiv.eq.0) iwpos = ioldps
  760 continue
!
! END OF LOOP ON TREE NODES.
!
      iw(1) = nblk
      if (ntwo.gt.0) iw(1) = -nblk
      nrlbdu = posfac - 1
      nirbdu = iwpos - 1
      if (ntotpv.eq.n) go to 810
      info(1) = 3
      info(2) = ntotpv
      go to 810
! **** ERROR RETURNS ****
  770 info(1) = -3
      go to 810

  780 info(1) = -4
      info(2) = la + max(posfac+lnass,apos2-ltopst+2) - astk
      go to 810

  790 info(1) = -5
      info(2) = ntotpv + 1
      go to 810

  800 info(1) = -6
      info(2) = ntotpv + 1
  810 continue
      info(9) = nrlbdu
      info(10) = nirbdu
      info(12) = ncmpbr
      info(13) = ncmpbi
      info(14) = ntwo
      info(15) = neig

      return
      end
      subroutine ma27pd(a,iw,j1,j2,itop,ireal,ncmpbr,ncmpbi)
! THIS SUBROUTINE PERFORMS A VERY SIMPLE COMPRESS (BLOCK MOVE).
!     ENTRIES J1 TO J2 (INCL.) IN A OR IW AS APPROPRIATE ARE MOVED TO
!     OCCUPY THE POSITIONS IMMEDIATELY PRIOR TO POSITION ITOP.
! A/IW HOLD THE ARRAY BEING COMPRESSED.
! J1/J2 DEFINE THE ENTRIES BEING MOVED.
! ITOP DEFINES THE POSITION IMMEDIATELY AFTER THE POSITIONS TO WHICH
!     J1 TO J2 ARE MOVED.
! IREAL MUST BE SET BY THE USER TO 2 IF THE MOVE IS ON ARRAY IW,
!     ANY OTHER VALUE WILL PERFORM THE MOVE ON A.
! NCMPBR and NCMPBI, see INFO(12) and INFO(13) in MA27A/AD (ACCUMULATE
!     THE NUMBER OF COMPRESSES OF THE REALS AND INTEGERS PERFORMED BY
!     MA27B/BD.
!
!     .. Scalar Arguments ..
      integer ireal,itop,j1,j2,ncmpbr,ncmpbi
!     ..
!     .. Array Arguments ..
      double precision a(*)
      integer iw(*)
!     ..
!     .. Local Scalars ..
      integer ipos,jj,jjj
!     ..
!     .. Executable Statements ..
      ipos = itop - 1
      if (j2.eq.ipos) go to 50
      if (ireal.eq.2) go to 20
      ncmpbr = ncmpbr + 1
      if (j1.gt.j2) go to 40
      do 10 jjj = j1,j2
        jj = j2 - jjj + j1
        a(ipos) = a(jj)
        ipos = ipos - 1
   10 continue
      go to 40

   20 ncmpbi = ncmpbi + 1
      if (j1.gt.j2) go to 40
      do 30 jjj = j1,j2
        jj = j2 - jjj + j1
        iw(ipos) = iw(jj)
        ipos = ipos - 1
   30 continue
   40 j2 = itop - 1
      j1 = ipos + 1
   50 return

      end
      subroutine ma27qd(n,a,la,iw,liw,w,maxfnt,rhs,iw2,nblk,latop,icntl)
! THIS SUBROUTINE PERFORMS FORWARD ELIMINATION
!     USING THE FACTOR U TRANSPOSE STORED IN A/IA AFTER MA27B/BD.
!
! N      - MUST BE SET TO THE ORDER OF THE MATRIX. NOT ALTERED
!          BY MA27Q/QD.
! A      - MUST BE SET TO HOLD THE REAL VALUES
!          CORRESPONDING TO THE FACTORS OF DINVERSE AND U.  THIS MUST BE
!          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
!          BY MA27Q/QD.
! LA     - LENGTH OF ARRAY A.  NOT ALTERED BY MA27Q/QD.
! IW     - HOLDS THE INTEGER INDEXING
!          INFORMATION FOR THE MATRIX FACTORS IN A.  THIS MUST BE
!          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
!          BY MA27Q/QD.
! LIW    - LENGTH OF ARRAY IW.  NOT ALTERED BY MA27Q/QD.
! W      - USED
!          AS WORKSPACE BY MA27Q/QD TO HOLD THE COMPONENTS OF THE RIGHT
!          HAND SIDES CORRESPONDING TO CURRENT BLOCK PIVOTAL ROWS.
! MAXFNT - MUST BE SET TO THE LARGEST NUMBER OF
!          VARIABLES IN ANY BLOCK PIVOT ROW.  THIS VALUE WILL HAVE
!          BEEN OUTPUT BY MA27B/BD.  NOT ALTERED BY MA27Q/QD.
! RHS    - ON INPUT,
!          MUST BE SET TO HOLD THE RIGHT HAND SIDES FOR THE EQUATIONS
!          WHICH THE USER DESIRES TO SOLVE.  ON OUTPUT, RHS WILL HOLD
!          THE MODIFIED VECTORS CORRESPONDING TO PERFORMING FORWARD
!          ELIMINATION ON THE RIGHT HAND SIDES.
! IW2    - NEED NOT BE SET ON ENTRY. ON EXIT IW2(I) (I = 1,NBLK)
!          WILL HOLD POINTERS TO THE
!          BEGINNING OF EACH BLOCK PIVOT IN ARRAY IW.
! NBLK   - NUMBER OF BLOCK PIVOT ROWS. NOT ALTERED BY MA27Q/QD.
! LATOP  - NEED NOT BE SET ON ENTRY. ON EXIT, IT IS THE POSITION IN
!          A OF THE LAST ENTRY IN THE  FACTORS. IT MUST BE PASSED
!          UNCHANGED TO MA27R/RD.
! ICNTL is an INTEGER array of length 30, see MA27A/AD.
!   ICNTL(IFRLVL+I) I=1,20 IS USED TO CONTROL WHETHER DIRECT OR INDIRECT
!     ACCESS IS USED BY MA27C/CD.  INDIRECT ACCESS IS EMPLOYED
!     IN FORWARD AND BACK SUBSTITUTION RESPECTIVELY IF THE SIZE OF
!     A BLOCK IS LESS THAN ICNTL(IFRLVL+MIN(10,NPIV)) AND
!     ICNTL(IFRLVL+10+MIN(10,NPIV)) RESPECTIVELY, WHERE NPIV IS THE
!     NUMBER OF PIVOTS IN THE BLOCK.
!
      integer ifrlvl
      parameter ( ifrlvl=5 )
!     .. Scalar Arguments ..
      integer la,latop,liw,maxfnt,n,nblk
!     ..
!     .. Array Arguments ..
      double precision a(la),rhs(n),w(maxfnt)
      integer iw(liw),iw2(nblk),icntl(30)
!     ..
!     .. Local Scalars ..
      double precision w1,w2
      integer apos,iblk,ifr,ilvl,ipiv,ipos,irhs,irow,ist,j,j1,j2,j3,jj, &
              jpiv,k,k1,k2,k3,liell,npiv
!     ..
!     .. Intrinsic Functions ..
      intrinsic abs,min
!     ..
!     .. Executable Statements ..
! APOS. RUNNING POINTER TO CURRENT PIVOT POSITION IN ARRAY A.
! IPOS. RUNNING POINTER TO BEGINNING OF BLOCK PIVOT ROW IN IW.
      apos = 1
      ipos = 1
      j2 = 0
      iblk = 0
      npiv = 0
      do 140 irow = 1,n
        if (npiv.gt.0) go to 90
        iblk = iblk + 1
        if (iblk.gt.nblk) go to 150
        ipos = j2 + 1
! SET UP POINTER IN PREPARATION FOR BACK SUBSTITUTION.
        iw2(iblk) = ipos
! ABS(LIELL) IS NUMBER OF VARIABLES (COLUMNS) IN BLOCK PIVOT ROW.
        liell = -iw(ipos)
! NPIV IS NUMBER OF PIVOTS (ROWS) IN BLOCK PIVOT.
        npiv = 1
        if (liell.gt.0) go to 10
        liell = -liell
        ipos = ipos + 1
        npiv = iw(ipos)
   10   j1 = ipos + 1
        j2 = ipos + liell
        ilvl = min(npiv,10)
        if (liell.lt.icntl(ifrlvl+ilvl)) go to 90
!
! PERFORM OPERATIONS USING DIRECT ADDRESSING.
!
! LOAD APPROPRIATE COMPONENTS OF RIGHT HAND SIDES INTO ARRAY W.
        ifr = 0
        do 20 jj = j1,j2
          j = abs(iw(jj)+0)
          ifr = ifr + 1
          w(ifr) = rhs(j)
   20   continue
! JPIV IS USED AS A FLAG SO THAT IPIV IS INCREMENTED CORRECTLY AFTER
! THE USE OF A 2 BY 2 PIVOT.
        jpiv = 1
        j3 = j1
! PERFORM OPERATIONS.
        do 70 ipiv = 1,npiv
          jpiv = jpiv - 1
          if (jpiv.eq.1) go to 70
! JUMP IF WE HAVE A 2 BY 2 PIVOT.
          if (iw(j3).lt.0) go to 40
! PERFORM FORWARD SUBSTITUTION USING 1 BY 1 PIVOT.
          jpiv = 1
          j3 = j3 + 1
          apos = apos + 1
          ist = ipiv + 1
          if (liell.lt.ist) go to 70
          w1 = w(ipiv)
          k = apos
          do 30 j = ist,liell
            w(j) = w(j) + a(k)*w1
            k = k + 1
   30     continue
          apos = apos + liell - ist + 1
          go to 70
! PERFORM OPERATIONS WITH 2 BY 2 PIVOT.
   40     jpiv = 2
          j3 = j3 + 2
          apos = apos + 2
          ist = ipiv + 2
          if (liell.lt.ist) go to 60
          w1 = w(ipiv)
          w2 = w(ipiv+1)
          k1 = apos
          k2 = apos + liell - ipiv
          do 50 j = ist,liell
            w(j) = w(j) + w1*a(k1) + w2*a(k2)
            k1 = k1 + 1
            k2 = k2 + 1
   50     continue
   60     apos = apos + 2* (liell-ist+1) + 1
   70   continue
! RELOAD W BACK INTO RHS.
        ifr = 0
        do 80 jj = j1,j2
          j = abs(iw(jj)+0)
          ifr = ifr + 1
          rhs(j) = w(ifr)
   80   continue
        npiv = 0
        go to 140
!
! PERFORM OPERATIONS USING INDIRECT ADDRESSING.
!
! JUMP IF WE HAVE A 2 BY 2 PIVOT.
   90   if (iw(j1).lt.0) go to 110
! PERFORM FORWARD SUBSTITUTION USING 1 BY 1 PIVOT.
        npiv = npiv - 1
        apos = apos + 1
        j1 = j1 + 1
        if (j1.gt.j2) go to 140
        irhs = iw(j1-1)
        w1 = rhs(irhs)
        k = apos
        do 100 j = j1,j2
          irhs = abs(iw(j)+0)
          rhs(irhs) = rhs(irhs) + a(k)*w1
          k = k + 1
  100   continue
        apos = apos + j2 - j1 + 1
        go to 140
! PERFORM OPERATIONS WITH 2 BY 2 PIVOT
  110   npiv = npiv - 2
        j1 = j1 + 2
        apos = apos + 2
        if (j1.gt.j2) go to 130
        irhs = -iw(j1-2)
        w1 = rhs(irhs)
        irhs = iw(j1-1)
        w2 = rhs(irhs)
        k1 = apos
        k3 = apos + j2 - j1 + 2
        do 120 j = j1,j2
          irhs = abs(iw(j)+0)
          rhs(irhs) = rhs(irhs) + w1*a(k1) + w2*a(k3)
          k1 = k1 + 1
          k3 = k3 + 1
  120   continue
  130   apos = apos + 2* (j2-j1+1) + 1
  140 continue
  150 latop = apos - 1
      return

      end
      subroutine ma27rd(n,a,la,iw,liw,w,maxfnt,rhs,iw2,nblk,latop,icntl)
! THIS SUBROUTINE PERFORMS BACKWARD ELIMINATION OPERATIONS
!     USING THE FACTORS DINVERSE AND U
!     STORED IN A/IW AFTER MA27B/BD.
!
! N      - MUST BE SET TO THE ORDER OF THE MATRIX. NOT ALTERED
!          BY MA27R/RD.
! A      - MUST BE SET TO HOLD THE REAL VALUES CORRESPONDING
!          TO THE FACTORS OF DINVERSE AND U.  THIS MUST BE
!          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
!          BY MA27R/RD.
! LA     - LENGTH OF ARRAY A. NOT ALTERED BY MA27R/RD.
! IW     - HOLDS THE INTEGER INDEXING
!          INFORMATION FOR THE MATRIX FACTORS IN A.  THIS MUST BE
!          UNCHANGED SINCE THE PRECEDING CALL TO MA27B/BD.  NOT ALTERED
!          BY MA27R/RD.
! LIW    - LENGTH OF ARRAY IW.  NOT ALTERED BY MA27R/RD.
! W      - USED
!          AS WORKSPACE BY MA27R/RD TO HOLD THE COMPONENTS OF THE RIGHT
!          HAND SIDES CORRESPONDING TO CURRENT BLOCK PIVOTAL ROWS.
! MAXFNT - INTEGER VARIABLE.  MUST BE SET TO THE LARGEST NUMBER OF
!          VARIABLES IN ANY BLOCK PIVOT ROW.  THIS VALUE WAS GIVEN
!          ON OUTPUT FROM MA27B/BD.  NOT ALTERED BY MA27R/RD.
! RHS    - ON INPUT,
!          MUST BE SET TO HOLD THE RIGHT HAND SIDE MODIFIED BY THE
!          FORWARD SUBSTITUTION OPERATIONS.  ON OUTPUT, RHS WILL HOLD
!          THE SOLUTION VECTOR.
! IW2    - ON ENTRY IW2(I) (I = 1,NBLK)
!          MUST HOLD POINTERS TO THE
!          BEGINNING OF EACH BLOCK PIVOT IN ARRAY IW, AS SET BY
!          MA27Q/QD.
! NBLK   - NUMBER OF BLOCK PIVOT ROWS. NOT ALTERED BY MA27R/RD.
! LATOP  - IT IS THE POSITION IN
!          A OF THE LAST ENTRY IN THE  FACTORS. IT MUST BE UNCHANGED
!          SINCE THE CALL TO MA27Q/QD.  IT IS NOT ALTERED BY MA27R/RD.
! ICNTL is an INTEGER array of length 30, see MA27A/AD.
!   ICNTL(IFRLVL+I) I=1,20 IS USED TO CONTROL WHETHER DIRECT OR INDIRECT
!     ACCESS IS USED BY MA27C/CD.  INDIRECT ACCESS IS EMPLOYED
!     IN FORWARD AND BACK SUBSTITUTION RESPECTIVELY IF THE SIZE OF
!     A BLOCK IS LESS THAN ICNTL(IFRLVL+MIN(10,NPIV)) AND
!     ICNTL(IFRLVL+10+MIN(10,NPIV)) RESPECTIVELY, WHERE NPIV IS THE
!     NUMBER OF PIVOTS IN THE BLOCK.
!
      integer ifrlvl
      parameter ( ifrlvl=5 )
!
!     .. Scalar Arguments ..
      integer la,latop,liw,maxfnt,n,nblk
!     ..
!     .. Array Arguments ..
      double precision a(la),rhs(n),w(maxfnt)
      integer iw(liw),iw2(nblk),icntl(30)
!     ..
!     .. Local Scalars ..
      double precision w1,w2
      integer apos,apos2,i1rhs,i2rhs,iblk,ifr,iipiv,iirhs,ilvl,ipiv, &
              ipos,irhs,ist,j,j1,j2,jj,jj1,jj2,jpiv,jpos,k,liell,loop, &
              npiv
!     ..
!     .. Intrinsic Functions ..
      intrinsic abs,min
!     ..
!     .. Executable Statements ..
! APOS. RUNNING POINTER TO CURRENT PIVOT POSITION IN ARRAY A.
! IPOS. RUNNING POINTER TO BEGINNING OF CURRENT BLOCK PIVOT ROW.
      apos = latop + 1
      npiv = 0
      iblk = nblk + 1
! RUN THROUGH BLOCK PIVOT ROWS IN THE REVERSE ORDER.
      do 180 loop = 1,n
        if (npiv.gt.0) go to 110
        iblk = iblk - 1
        if (iblk.lt.1) go to 190
        ipos = iw2(iblk)
! ABS(LIELL) IS NUMBER OF VARIABLES (COLUMNS) IN BLOCK PIVOT ROW.
        liell = -iw(ipos)
! NPIV IS NUMBER OF PIVOTS (ROWS) IN BLOCK PIVOT.
        npiv = 1
        if (liell.gt.0) go to 10
        liell = -liell
        ipos = ipos + 1
        npiv = iw(ipos)
   10   jpos = ipos + npiv
        j2 = ipos + liell
        ilvl = min(10,npiv) + 10
        if (liell.lt.icntl(ifrlvl+ilvl)) go to 110
!
! PERFORM OPERATIONS USING DIRECT ADDRESSING.
!
        j1 = ipos + 1
! LOAD APPROPRIATE COMPONENTS OF RHS INTO W.
        ifr = 0
        do 20 jj = j1,j2
          j = abs(iw(jj)+0)
          ifr = ifr + 1
          w(ifr) = rhs(j)
   20   continue
! JPIV IS USED AS A FLAG SO THAT IPIV IS INCREMENTED CORRECTLY AFTER
!     THE USE OF A 2 BY 2 PIVOT.
        jpiv = 1
! PERFORM ELIMINATIONS.
        do 90 iipiv = 1,npiv
          jpiv = jpiv - 1
          if (jpiv.eq.1) go to 90
          ipiv = npiv - iipiv + 1
          if (ipiv.eq.1) go to 30
! JUMP IF WE HAVE A 2 BY 2 PIVOT.
          if (iw(jpos-1).lt.0) go to 60
! PERFORM BACK-SUBSTITUTION USING 1 BY 1 PIVOT.
   30     jpiv = 1
          apos = apos - (liell+1-ipiv)
          ist = ipiv + 1
          w1 = w(ipiv)*a(apos)
          if (liell.lt.ist) go to 50
          jj1 = apos + 1
          do 40 j = ist,liell
            w1 = w1 + a(jj1)*w(j)
            jj1 = jj1 + 1
   40     continue
   50     w(ipiv) = w1
          jpos = jpos - 1
          go to 90
! PERFORM BACK-SUBSTITUTION OPERATIONS WITH 2 BY 2 PIVOT
   60     jpiv = 2
          apos2 = apos - (liell+1-ipiv)
          apos = apos2 - (liell+2-ipiv)
          ist = ipiv + 1
          w1 = w(ipiv-1)*a(apos) + w(ipiv)*a(apos+1)
          w2 = w(ipiv-1)*a(apos+1) + w(ipiv)*a(apos2)
          if (liell.lt.ist) go to 80
          jj1 = apos + 2
          jj2 = apos2 + 1
          do 70 j = ist,liell
            w1 = w1 + w(j)*a(jj1)
            w2 = w2 + w(j)*a(jj2)
            jj1 = jj1 + 1
            jj2 = jj2 + 1
   70     continue
   80     w(ipiv-1) = w1
          w(ipiv) = w2
          jpos = jpos - 2
   90   continue
! RELOAD WORKING VECTOR INTO SOLUTION VECTOR.
        ifr = 0
        do 100 jj = j1,j2
          j = abs(iw(jj)+0)
          ifr = ifr + 1
          rhs(j) = w(ifr)
  100   continue
        npiv = 0
        go to 180
!
! PERFORM OPERATIONS USING INDIRECT ADDRESSING.
!
  110   if (npiv.eq.1) go to 120
! JUMP IF WE HAVE A 2 BY 2 PIVOT.
        if (iw(jpos-1).lt.0) go to 150
! PERFORM BACK-SUBSTITUTION USING 1 BY 1 PIVOT.
  120   npiv = npiv - 1
        apos = apos - (j2-jpos+1)
        iirhs = iw(jpos)
        w1 = rhs(iirhs)*a(apos)
        j1 = jpos + 1
        if (j1.gt.j2) go to 140
        k = apos + 1
        do 130 j = j1,j2
          irhs = abs(iw(j)+0)
          w1 = w1 + a(k)*rhs(irhs)
          k = k + 1
  130   continue
  140   rhs(iirhs) = w1
        jpos = jpos - 1
        go to 180
! PERFORM OPERATIONS WITH 2 BY 2 PIVOT
  150   npiv = npiv - 2
        apos2 = apos - (j2-jpos+1)
        apos = apos2 - (j2-jpos+2)
        i1rhs = -iw(jpos-1)
        i2rhs = iw(jpos)
        w1 = rhs(i1rhs)*a(apos) + rhs(i2rhs)*a(apos+1)
        w2 = rhs(i1rhs)*a(apos+1) + rhs(i2rhs)*a(apos2)
        j1 = jpos + 1
        if (j1.gt.j2) go to 170
        jj1 = apos + 2
        jj2 = apos2 + 1
        do 160 j = j1,j2
          irhs = abs(iw(j)+0)
          w1 = w1 + rhs(irhs)*a(jj1)
          w2 = w2 + rhs(irhs)*a(jj2)
          jj1 = jj1 + 1
          jj2 = jj2 + 1
  160   continue
  170   rhs(i1rhs) = w1
        rhs(i2rhs) = w2
        jpos = jpos - 2
  180 continue
  190 return

      end
! *******************************************************************
! COPYRIGHT (c) 1999 Council for the Central Laboratory
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
! Original date July 1999
!CCCC PACKAGE MC64A/AD
!CCCC AUTHORS Iain Duff (i.duff@rl.ac.uk) and Jacko Koster (jak@ii.uib.no)

! 12th July 2004 Version 1.0.0. Version numbering added.

! 30/07/04  Version 1.1.0. Permutation array flagged negative to indicate
!           dependent columns in singular case.  Calls to MC64F changed
!           to avoid unsafe reference to array L.
! 21st February 2005 Version 1.2.0. FD05 dependence changed to FD15.

      subroutine mc64id(icntl)
      implicit none
      integer icntl(10)
      integer i
      icntl(1) = 6
      icntl(2) = 6
      icntl(3) = -1
      do 10 i = 4,10
        icntl(i) = 0
   10 continue
      return
      end
!**********************************************************************
      subroutine mc64ad(job,n,ne,ip,irn,a,num,cperm,liw,iw,ldw,dw, &
                 icntl,info)
      implicit none
      integer job,n,ne,num,liw,ldw
      integer ip(n+1),irn(ne),cperm(n),iw(liw),icntl(10),info(10)
      double precision a(ne),dw(ldw)
      integer i,j,k
      double precision fact,zero,rinf
      parameter (zero=0.0d+00)
      external fd15ad,mc21ad,mc64bd,mc64rd,mc64sd,mc64wd
      double precision fd15ad
      intrinsic abs,log
      rinf = fd15ad('H')
      if (job.lt.1 .or. job.gt.5) then
        info(1) = -1
        info(2) = job
        if (icntl(1).ge.0) write(icntl(1),9001) info(1),'JOB',job
        go to 99
      endif
      if (n.lt.1) then
        info(1) = -2
        info(2) = n
        if (icntl(1).ge.0) write(icntl(1),9001) info(1),'N',n
        go to 99
      endif
      if (ne.lt.1) then
        info(1) = -3
        info(2) = ne
        if (icntl(1).ge.0) write(icntl(1),9001) info(1),'NE',ne
        go to 99
      endif
      if (job.eq.1) k = 5*n
      if (job.eq.2) k = 4*n
      if (job.eq.3) k = 10*n + ne
      if (job.eq.4) k = 5*n
      if (job.eq.5) k = 5*n
      if (liw.lt.k) then
        info(1) = -4
        info(2) = k
        if (icntl(1).ge.0) write(icntl(1),9004) info(1),k
        go to 99
      endif
      if (job.gt.1) then
        if (job.eq.2) k = n
        if (job.eq.3) k = ne
        if (job.eq.4) k = 2*n + ne
        if (job.eq.5) k = 3*n + ne
        if (ldw.lt.k) then
          info(1) = -5
          info(2) = k
          if (icntl(1).ge.0) write(icntl(1),9005) info(1),k
          go to 99
        endif
      endif
      if (icntl(4).eq.0) then
        do 3 i = 1,n
          iw(i) = 0
    3   continue
        do 6 j = 1,n
          do 4 k = ip(j),ip(j+1)-1
            i = irn(k)
            if (i.lt.1 .or. i.gt.n) then
              info(1) = -6
              info(2) = j
              if (icntl(1).ge.0) write(icntl(1),9006) info(1),j,i
              go to 99
            endif
            if (iw(i).eq.j) then
              info(1) = -7
              info(2) = j
              if (icntl(1).ge.0) write(icntl(1),9007) info(1),j,i
              go to 99
            else
              iw(i) = j
            endif
    4     continue
    6   continue
      endif
      if (icntl(3).ge.0) then
        write(icntl(3),9020) job,n,ne
        write(icntl(3),9021) (ip(j),j=1,n+1)
        write(icntl(3),9022) (irn(j),j=1,ne)
        if (job.gt.1) write(icntl(3),9023) (a(j),j=1,ne)
      endif
      do 8 i=1,10
        info(i) = 0
    8 continue
      if (job.eq.1) then
        do 10 j = 1,n
          iw(j) = ip(j+1) - ip(j)
   10   continue
        call mc21ad(n,irn,ne,ip,iw(1),cperm,num,iw(n+1))
        go to 90
      endif
      if (job.eq.2) then
        call mc64bd(n,ne,ip,irn,a,cperm,num, &
           iw(1),iw(n+1),iw(2*n+1),iw(3*n+1),dw)
        go to 90
      endif
      if (job.eq.3) then
        do 20 k = 1,ne
          iw(k) = irn(k)
          dw(k) = abs(a(k))
   20   continue
        call mc64rd(n,ne,ip,iw,dw)
        call mc64sd(n,ne,ip,iw(1),dw,cperm,num,iw(ne+1), &
           iw(ne+n+1),iw(ne+2*n+1),iw(ne+3*n+1),iw(ne+4*n+1), &
           iw(ne+5*n+1),iw(ne+6*n+1))
        go to 90
      endif
      if (job.eq.4) then
        do 50 j = 1,n
          fact = zero
          do 30 k = ip(j),ip(j+1)-1
            if (abs(a(k)).gt.fact) fact = abs(a(k))
   30     continue
          do 40 k = ip(j),ip(j+1)-1
            dw(2*n+k) = fact - abs(a(k))
   40     continue
   50   continue
        call mc64wd(n,ne,ip,irn,dw(2*n+1),cperm,num, &
           iw(1),iw(n+1),iw(2*n+1),iw(3*n+1),iw(4*n+1), &
           dw(1),dw(n+1))
        go to 90
      endif
      if (job.eq.5) then
        do 75 j = 1,n
          fact = zero
          do 60 k = ip(j),ip(j+1)-1
            dw(3*n+k) = abs(a(k))
            if (dw(3*n+k).gt.fact) fact = dw(3*n+k)
   60     continue
          dw(2*n+j) = fact
          if (fact.ne.zero) then
            fact = log(fact)
          else
            fact = rinf/n
          endif
          do 70 k = ip(j),ip(j+1)-1
            if (dw(3*n+k).ne.zero) then
              dw(3*n+k) = fact - log(dw(3*n+k))
            else
              dw(3*n+k) = rinf/n
            endif
   70     continue
   75   continue
        call mc64wd(n,ne,ip,irn,dw(3*n+1),cperm,num, &
           iw(1),iw(n+1),iw(2*n+1),iw(3*n+1),iw(4*n+1), &
           dw(1),dw(n+1))
        if (num.eq.n) then
          do 80 j = 1,n
            if (dw(2*n+j).ne.zero) then
              dw(n+j) = dw(n+j) - log(dw(2*n+j))
            else
              dw(n+j) = zero
            endif
   80     continue
        endif
        fact = 0.5*log(rinf)
        do 86 j = 1,n
          if (dw(j).lt.fact .and. dw(n+j).lt.fact) go to 86
          info(1) = 2
          go to 90
   86   continue
      endif
   90 if (info(1).eq.0 .and. num.lt.n) then
        info(1) = 1
        if (icntl(2).ge.0) write(icntl(2),9011) info(1)
      endif
      if (info(1).eq.2) then
        if (icntl(2).ge.0) write(icntl(2),9012) info(1)
      endif
      if (icntl(3).ge.0) then
        write(icntl(3),9030) (info(j),j=1,2)
        write(icntl(3),9031) num
        write(icntl(3),9032) (cperm(j),j=1,n)
        if (job.eq.5) then
          write(icntl(3),9033) (dw(j),j=1,n)
          write(icntl(3),9034) (dw(n+j),j=1,n)
        endif
      endif
   99 return
 9001 format (' ****** Error in MC64A/AD. INFO(1) = ',i2, &
              ' because ',(a),' = ',i10)
 9004 format (' ****** Error in MC64A/AD. INFO(1) = ',i2/ &
              '        LIW too small, must be at least ',i8)
 9005 format (' ****** Error in MC64A/AD. INFO(1) = ',i2/ &
              '        LDW too small, must be at least ',i8)
 9006 format (' ****** Error in MC64A/AD. INFO(1) = ',i2/ &
              '        Column ',i8, &
              ' contains an entry with invalid row index ',i8)
 9007 format (' ****** Error in MC64A/AD. INFO(1) = ',i2/ &
              '        Column ',i8, &
              ' contains two or more entries with row index ',i8)
 9011 format (' ****** Warning from MC64A/AD. INFO(1) = ',i2/ &
              '        The matrix is structurally singular.')
 9012 format (' ****** Warning from MC64A/AD. INFO(1) = ',i2/ &
              '        Some scaling factors may be too large.')
 9020 format (' ****** Input parameters for MC64A/AD:'/ &
              ' JOB = ',i8/' N   = ',i8/' NE  = ',i8)
 9021 format (' IP(1:N+1)  = ',8i8/(14x,8i8))
 9022 format (' IRN(1:NE)  = ',8i8/(14x,8i8))
 9023 format (' A(1:NE)    = ',4(1pd14.4)/(14x,4(1pd14.4)))
 9030 format (' ****** Output parameters for MC64A/AD:'/ &
              ' INFO(1:2)  = ',2i8)
 9031 format (' NUM        = ',i8)
 9032 format (' CPERM(1:N) = ',8i8/(14x,8i8))
 9033 format (' DW(1:N)    = ',5(f11.3)/(14x,5(f11.3)))
 9034 format (' DW(N+1:2N) = ',5(f11.3)/(14x,5(f11.3)))
      end
!**********************************************************************
      subroutine mc64bd(n,ne,ip,irn,a,iperm,num,jperm,pr,q,l,d)
      implicit none
      integer n,ne,num
      integer ip(n+1),irn(ne),iperm(n),jperm(n),pr(n),q(n),l(n)
      double precision a(ne),d(n)
      integer i,ii,j,jj,jord,q0,qlen,idum,jdum,isp,jsp, &
              k,kk,kk1,kk2,i0,up,low,lpos
      double precision csp,di,dnew,dq0,ai,a0,bv
      double precision rinf,zero,minone
      parameter (zero=0.0d+0,minone=-1.0d+0)
      intrinsic abs,min
      external fd15ad,mc64dd,mc64ed,mc64fd
      double precision fd15ad
      rinf = fd15ad('H')
      num = 0
      bv = rinf
      do 10 k = 1,n
        iperm(k) = 0
        jperm(k) = 0
        pr(k) = ip(k)
        d(k) = zero
   10 continue
      do 20 j = 1,n
        a0 = minone
        do 30 k = ip(j),ip(j+1)-1
          i = irn(k)
          ai = abs(a(k))
          if (ai.gt.d(i)) d(i) = ai
          if (jperm(j).ne.0) go to 30
          if (ai.ge.bv) then
            a0 = bv
            if (iperm(i).ne.0) go to 30
            jperm(j) = i
            iperm(i) = j
            num = num + 1
          else
            if (ai.le.a0) go to 30
            a0 = ai
            i0 = i
          endif
   30   continue
        if (a0.ne.minone .and. a0.lt.bv) then
          bv = a0
          if (iperm(i0).ne.0) go to 20
          iperm(i0) = j
          jperm(j) = i0
          num = num + 1
        endif
   20 continue
      do 25 i = 1,n
        bv = min(bv,d(i))
   25 continue
      if (num.eq.n) go to 1000
      do 95 j = 1,n
        if (jperm(j).ne.0) go to 95
        do 50 k = ip(j),ip(j+1)-1
          i = irn(k)
          ai = abs(a(k))
          if (ai.lt.bv) go to 50
          if (iperm(i).eq.0) go to 90
          jj = iperm(i)
          kk1 = pr(jj)
          kk2 = ip(jj+1) - 1
          if (kk1.gt.kk2) go to 50
          do 70 kk = kk1,kk2
            ii = irn(kk)
            if (iperm(ii).ne.0) go to 70
            if (abs(a(kk)).ge.bv) go to 80
   70     continue
          pr(jj) = kk2 + 1
   50   continue
        go to 95
   80   jperm(jj) = ii
        iperm(ii) = jj
        pr(jj) = kk + 1
   90   num = num + 1
        jperm(j) = i
        iperm(i) = j
        pr(j) = k + 1
   95 continue
      if (num.eq.n) go to 1000
      do 99 i = 1,n
        d(i) = minone
        l(i) = 0
   99 continue
      do 100 jord = 1,n
        if (jperm(jord).ne.0) go to 100
        qlen = 0
        low = n + 1
        up = n + 1
        csp = minone
        j = jord
        pr(j) = -1
        do 115 k = ip(j),ip(j+1)-1
          i = irn(k)
          dnew = abs(a(k))
          if (csp.ge.dnew) go to 115
          if (iperm(i).eq.0) then
            csp = dnew
            isp = i
            jsp = j
            if (csp.ge.bv) go to 160
          else
            d(i) = dnew
            if (dnew.ge.bv) then
              low = low - 1
              q(low) = i
            else
              qlen = qlen + 1
              l(i) = qlen
              call mc64dd(i,n,q,d,l,1)
            endif
            jj = iperm(i)
            pr(jj) = j
          endif
  115   continue
        do 150 jdum = 1,num
          if (low.eq.up) then
            if (qlen.eq.0) go to 160
            i = q(1)
            if (csp.ge.d(i)) go to 160
            bv = d(i)
            do 152 idum = 1,n
              call mc64ed(qlen,n,q,d,l,1)
              l(i) = 0
              low = low - 1
              q(low) = i
              if (qlen.eq.0) go to 153
              i = q(1)
              if (d(i).ne.bv) go to 153
  152       continue
          endif
  153     up = up - 1
          q0 = q(up)
          dq0 = d(q0)
          l(q0) = up
          j = iperm(q0)
          do 155 k = ip(j),ip(j+1)-1
            i = irn(k)
            if (l(i).ge.up) go to 155
            dnew = min(dq0,abs(a(k)))
            if (csp.ge.dnew) go to 155
            if (iperm(i).eq.0) then
              csp = dnew
              isp = i
              jsp = j
              if (csp.ge.bv) go to 160
            else
              di = d(i)
              if (di.ge.bv .or. di.ge.dnew) go to 155
              d(i) = dnew
              if (dnew.ge.bv) then
                if (di.ne.minone) then
                  lpos = l(i)
                  call mc64fd(lpos,qlen,n,q,d,l,1)
                endif
                l(i) = 0
                low = low - 1
                q(low) = i
              else
                if (di.eq.minone) then
                  qlen = qlen + 1
                  l(i) = qlen
                endif
                call mc64dd(i,n,q,d,l,1)
              endif
              jj = iperm(i)
              pr(jj) = j
            endif
  155     continue
  150   continue
  160   if (csp.eq.minone) go to 190
        bv = min(bv,csp)
        num = num + 1
        i = isp
        j = jsp
        do 170 jdum = 1,num+1
          i0 = jperm(j)
          jperm(j) = i
          iperm(i) = j
          j = pr(j)
          if (j.eq.-1) go to 190
          i = i0
  170   continue
  190   do 191 kk = up,n
          i = q(kk)
          d(i) = minone
          l(i) = 0
  191   continue
        do 192 kk = low,up-1
          i = q(kk)
          d(i) = minone
  192   continue
        do 193 kk = 1,qlen
          i = q(kk)
          d(i) = minone
          l(i) = 0
  193   continue
  100 continue
      if (num.eq.n) go to 1000
      do 300 j = 1,n
        jperm(j) = 0
  300 continue
      k = 0
      do 310 i = 1,n
        if (iperm(i).eq.0) then
          k = k + 1
          pr(k) = i
        else
          j = iperm(i)
          jperm(j) = i
        endif
  310 continue
      k = 0
      do 320 i = 1,n
        if (jperm(i).ne.0) go to 320
        k = k + 1
        jdum = pr(k)
        iperm(jdum) = i
  320 continue
 1000 return
      end
!**********************************************************************
      subroutine mc64dd(i,n,q,d,l,iway)
      implicit none
      integer i,n,iway
      integer q(n),l(n)
      double precision d(n)
      integer idum,k,pos,posk,qk
      parameter (k=2)
      double precision di
      di = d(i)
      pos = l(i)
      if (iway.eq.1) then
        do 10 idum = 1,n
          if (pos.le.1) go to 20
          posk = pos/k
          qk = q(posk)
          if (di.le.d(qk)) go to 20
          q(pos) = qk
          l(qk) = pos
          pos = posk
   10   continue
      else
        do 15 idum = 1,n
          if (pos.le.1) go to 20
          posk = pos/k
          qk = q(posk)
          if (di.ge.d(qk)) go to 20
          q(pos) = qk
          l(qk) = pos
          pos = posk
   15   continue
      endif
   20 q(pos) = i
      l(i) = pos
      return
      end
!**********************************************************************
      subroutine mc64ed(qlen,n,q,d,l,iway)
      implicit none
      integer qlen,n,iway
      integer q(n),l(n)
      double precision d(n)
      integer i,idum,k,pos,posk
      parameter (k=2)
      double precision dk,dr,di
      i = q(qlen)
      di = d(i)
      qlen = qlen - 1
      pos = 1
      if (iway.eq.1) then
        do 10 idum = 1,n
          posk = k*pos
          if (posk.gt.qlen) go to 20
          dk = d(q(posk))
          if (posk.lt.qlen) then
            dr = d(q(posk+1))
            if (dk.lt.dr) then
              posk = posk + 1
              dk = dr
            endif
          endif
          if (di.ge.dk) go to 20
          q(pos) = q(posk)
          l(q(pos)) = pos
          pos = posk
   10   continue
      else
        do 15 idum = 1,n
          posk = k*pos
          if (posk.gt.qlen) go to 20
          dk = d(q(posk))
          if (posk.lt.qlen) then
            dr = d(q(posk+1))
            if (dk.gt.dr) then
              posk = posk + 1
              dk = dr
            endif
          endif
          if (di.le.dk) go to 20
          q(pos) = q(posk)
          l(q(pos)) = pos
          pos = posk
   15   continue
      endif
   20 q(pos) = i
      l(i) = pos
      return
      end
!**********************************************************************
      subroutine mc64fd(pos0,qlen,n,q,d,l,iway)
      implicit none
      integer pos0,qlen,n,iway
      integer q(n),l(n)
      double precision d(n)
      integer i,idum,k,pos,posk,qk
      parameter (k=2)
      double precision dk,dr,di
      if (qlen.eq.pos0) then
        qlen = qlen - 1
        return
      endif
      i = q(qlen)
      di = d(i)
      qlen = qlen - 1
      pos = pos0
      if (iway.eq.1) then
        do 10 idum = 1,n
          if (pos.le.1) go to 20
          posk = pos/k
          qk = q(posk)
          if (di.le.d(qk)) go to 20
          q(pos) = qk
          l(qk) = pos
          pos = posk
   10   continue
   20   q(pos) = i
        l(i) = pos
        do 30 idum = 1,n
          posk = k*pos
          if (posk.gt.qlen) go to 40
          dk = d(q(posk))
          if (posk.lt.qlen) then
            dr = d(q(posk+1))
            if (dk.lt.dr) then
              posk = posk + 1
              dk = dr
            endif
          endif
          if (di.ge.dk) go to 40
          qk = q(posk)
          q(pos) = qk
          l(qk) = pos
          pos = posk
   30   continue
      else
        do 32 idum = 1,n
          if (pos.le.1) go to 34
          posk = pos/k
          qk = q(posk)
          if (di.ge.d(qk)) go to 34
          q(pos) = qk
          l(qk) = pos
          pos = posk
   32   continue
   34   q(pos) = i
        l(i) = pos
        do 36 idum = 1,n
          posk = k*pos
          if (posk.gt.qlen) go to 40
          dk = d(q(posk))
          if (posk.lt.qlen) then
            dr = d(q(posk+1))
            if (dk.gt.dr) then
              posk = posk + 1
              dk = dr
            endif
          endif
          if (di.le.dk) go to 40
          qk = q(posk)
          q(pos) = qk
          l(qk) = pos
          pos = posk
   36   continue
      endif
   40 q(pos) = i
      l(i) = pos
      return
      end
!**********************************************************************
      subroutine mc64rd(n,ne,ip,irn,a)
      implicit none
      integer n,ne
      integer ip(n+1),irn(ne)
      double precision a(ne)
      integer thresh,tdlen
      parameter (thresh=15,tdlen=50)
      integer j,ipj,k,len,r,s,hi,first,mid,last,td
      double precision ha,key
      integer todo(tdlen)
      do 100 j = 1,n
        len = ip(j+1) - ip(j)
        if (len.le.1) go to 100
        ipj = ip(j)
        if (len.lt.thresh) go to 400
        todo(1) = ipj
        todo(2) = ipj + len
        td = 2
  500   continue
        first = todo(td-1)
        last = todo(td)
        key = a((first+last)/2)
        do 475 k = first,last-1
          ha = a(k)
          if (ha.eq.key) go to 475
          if (ha.gt.key) go to 470
          key = ha
          go to 470
  475   continue
        td = td - 2
        go to 425
  470   mid = first
        do 450 k = first,last-1
          if (a(k).le.key) go to 450
          ha = a(mid)
          a(mid) = a(k)
          a(k) = ha
          hi = irn(mid)
          irn(mid) = irn(k)
          irn(k) = hi
          mid = mid + 1
  450   continue
        if (mid-first.ge.last-mid) then
          todo(td+2) = last
          todo(td+1) = mid
          todo(td) = mid
        else
          todo(td+2) = mid
          todo(td+1) = first
          todo(td) = last
          todo(td-1) = mid
        endif
        td = td + 2
  425   continue
        if (td.eq.0) go to 400
        if (todo(td)-todo(td-1).ge.thresh) go to 500
        td = td - 2
        go to 425
  400   do 200 r = ipj+1,ipj+len-1
          if (a(r-1) .lt. a(r)) then
            ha = a(r)
            hi = irn(r)
            a(r) = a(r-1)
            irn(r) = irn(r-1)
            do 300 s = r-1,ipj+1,-1
              if (a(s-1) .lt. ha) then
                a(s) = a(s-1)
                irn(s) = irn(s-1)
              else
                a(s) = ha
                irn(s) = hi
                go to 200
              end if
  300       continue
            a(ipj) = ha
            irn(ipj) = hi
          end if
  200   continue
  100 continue
      return
      end
!**********************************************************************
      subroutine mc64sd(n,ne,ip,irn,a,iperm,numx, &
                 w,len,lenl,lenh,fc,iw,iw4)
      implicit none
      integer n,ne,numx
      integer ip(n+1),irn(ne),iperm(n), &
              w(n),len(n),lenl(n),lenh(n),fc(n),iw(n),iw4(4*n)
      double precision a(ne)
      integer num,nval,wlen,ii,i,j,k,l,cnt,mod,idum1,idum2,idum3
      double precision bval,bmin,bmax,rinf
      external fd15ad,mc64qd,mc64ud
      double precision fd15ad
      rinf = fd15ad('H')
      do 20 j = 1,n
        fc(j) = j
        iw(j) = 0
        len(j) = ip(j+1) - ip(j)
   20 continue
      cnt = 1
      mod = 1
      numx = 0
      call mc64ud(cnt,mod,n,irn,ne,ip,len,fc,iw,numx,n, &
                  iw4(1),iw4(n+1),iw4(2*n+1),iw4(3*n+1))
      num = numx
      if (num.ne.n) then
        bmax = rinf
      else
        bmax = rinf
        do 30 j = 1,n
          bval = 0.0
          do 25 k = ip(j),ip(j+1)-1
            if (a(k).gt.bval) bval = a(k)
   25     continue
          if (bval.lt.bmax) bmax = bval
   30   continue
        bmax = 1.001 * bmax
      endif
      bval = 0.0
      bmin = 0.0
      wlen = 0
      do 48 j = 1,n
        l = ip(j+1) - ip(j)
        lenh(j) = l
        len(j) = l
        do 45 k = ip(j),ip(j+1)-1
          if (a(k).lt.bmax) go to 46
   45   continue
        k = ip(j+1)
   46   lenl(j) = k - ip(j)
        if (lenl(j).eq.l) go to 48
        wlen = wlen + 1
        w(wlen) = j
   48 continue
      do 90 idum1 = 1,ne
        if (num.eq.numx) then
          do 50 i = 1,n
            iperm(i) = iw(i)
   50     continue
          do 80 idum2 = 1,ne
            bmin = bval
            if (bmax .eq. bmin) go to 99
            call mc64qd(ip,lenl,len,w,wlen,a,nval,bval)
            if (nval.le.1) go to 99
            k = 1
            do 70 idum3 = 1,n
              if (k.gt.wlen) go to 71
              j = w(k)
              do 55 ii = ip(j)+len(j)-1,ip(j)+lenl(j),-1
                if (a(ii).ge.bval) go to 60
                i = irn(ii)
                if (iw(i).ne.j) go to 55
                iw(i) = 0
                num = num - 1
                fc(n-num) = j
   55         continue
   60         lenh(j) = len(j)
              len(j) = ii - ip(j) + 1
              if (lenl(j).eq.lenh(j)) then
                w(k) = w(wlen)
                wlen = wlen - 1
              else
                k = k + 1
              endif
   70       continue
   71       if (num.lt.numx) go to 81
   80     continue
   81     mod = 1
        else
          bmax = bval
          call mc64qd(ip,len,lenh,w,wlen,a,nval,bval)
          if (nval.eq.0 .or. bval.eq.bmin) go to 99
          k = 1
          do 87 idum3 = 1,n
            if (k.gt.wlen) go to 88
            j = w(k)
            do 85 ii = ip(j)+len(j),ip(j)+lenh(j)-1
              if (a(ii).lt.bval) go to 86
   85       continue
   86       lenl(j) = len(j)
            len(j) = ii - ip(j)
            if (lenl(j).eq.lenh(j)) then
              w(k) = w(wlen)
              wlen = wlen - 1
            else
              k = k + 1
            endif
   87     continue
   88     mod = 0
        endif
        cnt = cnt + 1
        call mc64ud(cnt,mod,n,irn,ne,ip,len,fc,iw,num,numx, &
                  iw4(1),iw4(n+1),iw4(2*n+1),iw4(3*n+1))
   90 continue
   99 if (numx.eq.n) go to 1000
      do 300 j = 1,n
        w(j) = 0
  300 continue
      k = 0
      do 310 i = 1,n
        if (iperm(i).eq.0) then
          k = k + 1
          iw(k) = i
        else
          j = iperm(i)
          w(j) = i
        endif
  310 continue
      k = 0
      do 320 j = 1,n
        if (w(j).ne.0) go to 320
        k = k + 1
        idum1 = iw(k)
        iperm(idum1) = j
  320 continue
 1000 return
      end
!**********************************************************************
      subroutine mc64qd(ip,lenl,lenh,w,wlen,a,nval,val)
      implicit none
      integer wlen,nval
      integer ip(*),lenl(*),lenh(*),w(*)
      double precision a(*),val
      integer xx,j,k,ii,s,pos
      parameter (xx=10)
      double precision split(xx),ha
      nval = 0
      do 10 k = 1,wlen
        j = w(k)
        do 15 ii = ip(j)+lenl(j),ip(j)+lenh(j)-1
          ha = a(ii)
          if (nval.eq.0) then
            split(1) = ha
            nval = 1
          else
            do 20 s = nval,1,-1
              if (split(s).eq.ha) go to 15
              if (split(s).gt.ha) then
                pos = s + 1
                go to 21
              endif
  20        continue
            pos = 1
  21        do 22 s = nval,pos,-1
              split(s+1) = split(s)
  22        continue
            split(pos) = ha
            nval = nval + 1
          endif
          if (nval.eq.xx) go to 11
  15    continue
  10  continue
  11  if (nval.gt.0) val = split((nval+1)/2)
      return
      end
!**********************************************************************
      subroutine mc64ud(id,mod,n,irn,lirn,ip,lenc,fc,iperm,num,numx, &
                 pr,arp,cv,out)
      implicit none
      integer id,mod,n,lirn,num,numx
      integer arp(n),cv(n),irn(lirn),ip(n), &
              fc(n),iperm(n),lenc(n),out(n),pr(n)
      integer i,ii,in1,in2,j,j1,jord,k,kk,last,nfc, &
              num0,num1,num2,id0,id1
      if (id.eq.1) then
        do 5 i = 1,n
          cv(i) = 0
          arp(i) = 0
    5   continue
        num1 = n
        num2 = n
      else
        if (mod.eq.1) then
          do 8 i = 1,n
            arp(i) = 0
    8     continue
        endif
        num1 = numx
        num2 = n - numx
      endif
      num0 = num
      nfc = 0
      id0 = (id-1)*n
      do 100 jord = num0+1,n
        id1 = id0 + jord
        j = fc(jord-num0)
        pr(j) = -1
        do 70 k = 1,jord
          if (arp(j).ge.lenc(j)) go to 30
          in1 = ip(j) + arp(j)
          in2 = ip(j) + lenc(j) - 1
          do 20 ii = in1,in2
            i = irn(ii)
            if (iperm(i).eq.0) go to 80
   20     continue
          arp(j) = lenc(j)
   30     out(j) = lenc(j) - 1
          do 60 kk = 1,jord
            in1 = out(j)
            if (in1.lt.0) go to 50
            in2 = ip(j) + lenc(j) - 1
            in1 = in2 - in1
            do 40 ii = in1,in2
              i = irn(ii)
              if (cv(i).eq.id1) go to 40
              j1 = j
              j = iperm(i)
              cv(i) = id1
              pr(j) = j1
              out(j1) = in2 - ii - 1
              go to 70
   40       continue
   50       j1 = pr(j)
            if (j1.eq.-1) then
              nfc = nfc + 1
              fc(nfc) = j
              if (nfc.gt.num2) then
                last = jord
                go to 101
              endif
              go to 100
            endif
            j = j1
   60     continue
   70   continue
   80   iperm(i) = j
        arp(j) = ii - ip(j) + 1
        num = num + 1
        do 90 k = 1,jord
          j = pr(j)
          if (j.eq.-1) go to 95
          ii = ip(j) + lenc(j) - out(j) - 2
          i = irn(ii)
          iperm(i) = j
   90   continue
   95   if (num.eq.num1) then
          last = jord
          go to 101
        endif
  100 continue
      last = n
  101 do 110 jord = last+1,n
        nfc = nfc + 1
        fc(nfc) = fc(jord-num0)
  110 continue
      return
      end
!**********************************************************************
      subroutine mc64wd(n,ne,ip,irn,a,iperm,num, &
                 jperm,out,pr,q,l,u,d)
      implicit none
      integer n,ne,num
      integer ip(n+1),irn(ne),iperm(n), &
              jperm(n),out(n),pr(n),q(n),l(n)
      double precision a(ne),u(n),d(n)
      integer i,i0,ii,j,jj,jord,q0,qlen,jdum,isp,jsp, &
              k,k0,k1,k2,kk,kk1,kk2,up,low,lpos
      double precision csp,di,dmin,dnew,dq0,vj
      double precision rinf,zero
      parameter (zero=0.0d+0)
      external fd15ad,mc64dd,mc64ed,mc64fd
      double precision fd15ad
      rinf = fd15ad('H')
      num = 0
      do 10 k = 1,n
        u(k) = rinf
        d(k) = zero
        iperm(k) = 0
        jperm(k) = 0
        pr(k) = ip(k)
        l(k) = 0
   10 continue
      do 30 j = 1,n
        do 20 k = ip(j),ip(j+1)-1
          i = irn(k)
          if (a(k).gt.u(i)) go to 20
          u(i) = a(k)
          iperm(i) = j
          l(i) = k
   20   continue
   30 continue
      do 40 i = 1,n
        j = iperm(i)
        if (j.eq.0) go to 40
        iperm(i) = 0
        if (jperm(j).ne.0) go to 40
        num = num + 1
        iperm(i) = j
        jperm(j) = l(i)
   40 continue
      if (num.eq.n) go to 1000
      do 95 j = 1,n
        if (jperm(j).ne.0) go to 95
        k1 = ip(j)
        k2 = ip(j+1) - 1
        if (k1.gt.k2) go to 95
        vj = rinf
        do 50 k = k1,k2
          i = irn(k)
          di = a(k) - u(i)
          if (di.gt.vj) go to 50
          if (di.lt.vj .or. di.eq.rinf) go to 55
          if (iperm(i).ne.0 .or. iperm(i0).eq.0) go to 50
   55     vj = di
          i0 = i
          k0 = k
   50   continue
        d(j) = vj
        k = k0
        i = i0
        if (iperm(i).eq.0) go to 90
        do 60 k = k0,k2
          i = irn(k)
          if (a(k)-u(i).gt.vj) go to 60
          jj = iperm(i)
          kk1 = pr(jj)
          kk2 = ip(jj+1) - 1
          if (kk1.gt.kk2) go to 60
          do 70 kk = kk1,kk2
            ii = irn(kk)
            if (iperm(ii).gt.0) go to 70
            if (a(kk)-u(ii).le.d(jj)) go to 80
   70     continue
          pr(jj) = kk2 + 1
   60   continue
        go to 95
   80   jperm(jj) = kk
        iperm(ii) = jj
        pr(jj) = kk + 1
   90   num = num + 1
        jperm(j) = k
        iperm(i) = j
        pr(j) = k + 1
   95 continue
      if (num.eq.n) go to 1000
      do 99 i = 1,n
        d(i) = rinf
        l(i) = 0
   99 continue
      do 100 jord = 1,n
        if (jperm(jord).ne.0) go to 100
        dmin = rinf
        qlen = 0
        low = n + 1
        up = n + 1
        csp = rinf
        j = jord
        pr(j) = -1
        do 115 k = ip(j),ip(j+1)-1
          i = irn(k)
          dnew = a(k) - u(i)
          if (dnew.ge.csp) go to 115
          if (iperm(i).eq.0) then
            csp = dnew
            isp = k
            jsp = j
          else
            if (dnew.lt.dmin) dmin = dnew
            d(i) = dnew
            qlen = qlen + 1
            q(qlen) = k
          endif
  115   continue
        q0 = qlen
        qlen = 0
        do 120 kk = 1,q0
          k = q(kk)
          i = irn(k)
          if (csp.le.d(i)) then
            d(i) = rinf
            go to 120
          endif
          if (d(i).le.dmin) then
            low = low - 1
            q(low) = i
            l(i) = low
          else
            qlen = qlen + 1
            l(i) = qlen
            call mc64dd(i,n,q,d,l,2)
          endif
          jj = iperm(i)
          out(jj) = k
          pr(jj) = j
  120   continue
        do 150 jdum = 1,num
          if (low.eq.up) then
            if (qlen.eq.0) go to 160
            i = q(1)
            if (d(i).ge.csp) go to 160
            dmin = d(i)
  152       call mc64ed(qlen,n,q,d,l,2)
            low = low - 1
            q(low) = i
            l(i) = low
            if (qlen.eq.0) go to 153
            i = q(1)
            if (d(i).gt.dmin) go to 153
            go to 152
          endif
  153     q0 = q(up-1)
          dq0 = d(q0)
          if (dq0.ge.csp) go to 160
          up = up - 1
          j = iperm(q0)
          vj = dq0 - a(jperm(j)) + u(q0)
          do 155 k = ip(j),ip(j+1)-1
            i = irn(k)
            if (l(i).ge.up) go to 155
            dnew = vj + a(k)-u(i)
            if (dnew.ge.csp) go to 155
            if (iperm(i).eq.0) then
              csp = dnew
              isp = k
              jsp = j
            else
              di = d(i)
              if (di.le.dnew) go to 155
              if (l(i).ge.low) go to 155
              d(i) = dnew
              if (dnew.le.dmin) then
                lpos = l(i)
                if (lpos.ne.0) &
                  call mc64fd(lpos,qlen,n,q,d,l,2)
                low = low - 1
                q(low) = i
                l(i) = low
              else
                if (l(i).eq.0) then
                  qlen = qlen + 1
                  l(i) = qlen
                endif
                call mc64dd(i,n,q,d,l,2)
              endif
              jj = iperm(i)
              out(jj) = k
              pr(jj) = j
            endif
  155     continue
  150   continue
  160   if (csp.eq.rinf) go to 190
        num = num + 1
        i = irn(isp)
        iperm(i) = jsp
        jperm(jsp) = isp
        j = jsp
        do 170 jdum = 1,num
          jj = pr(j)
          if (jj.eq.-1) go to 180
          k = out(j)
          i = irn(k)
          iperm(i) = jj
          jperm(jj) = k
          j = jj
  170   continue
  180   do 185 kk = up,n
          i = q(kk)
          u(i) = u(i) + d(i) - csp
  185   continue
  190   do 191 kk = low,n
          i = q(kk)
          d(i) = rinf
          l(i) = 0
  191   continue
        do 193 kk = 1,qlen
          i = q(kk)
          d(i) = rinf
          l(i) = 0
  193   continue
  100 continue
 1000 do 200 j = 1,n
        k = jperm(j)
        if (k.ne.0) then
          d(j) = a(k) - u(irn(k))
        else
          d(j) = zero
        endif
        if (iperm(j).eq.0) u(j) = zero
  200 continue
      if (num.eq.n) go to 1100
      do 300 j = 1,n
        jperm(j) = 0
  300 continue
      k = 0
      do 310 i = 1,n
        if (iperm(i).eq.0) then
          k = k + 1
          out(k) = i
        else
          j = iperm(i)
          jperm(j) = i
        endif
  310 continue
      k = 0
      do 320 j = 1,n
        if (jperm(j).ne.0) go to 320
        k = k + 1
        jdum = out(k)
        iperm(jdum) = - j
  320 continue
 1100 return
      end
! *******************************************************************
! COPYRIGHT (c) 1977 Hyprotech UK
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
! Original date 8 Oct 1992
!######8/10/92 Toolpack tool decs employed.
!######8/10/92 D version created by name change only.
! 13/3/02 Cosmetic changes applied to reduce single/double differences
!
! 12th July 2004 Version 1.0.0. Version numbering added.

      subroutine mc21ad(n,icn,licn,ip,lenr,iperm,numnz,iw)
      integer licn,n,numnz
      integer icn(licn),ip(n),iperm(n),iw(n,4),lenr(n)
      external mc21bd
      call mc21bd(n,icn,licn,ip,lenr,iperm,numnz,iw(1,1),iw(1,2), &
                  iw(1,3),iw(1,4))
      return
      end
      subroutine mc21bd(n,icn,licn,ip,lenr,iperm,numnz,pr,arp,cv,out)
      integer licn,n,numnz
      integer arp(n),cv(n),icn(licn),ip(n),iperm(n),lenr(n),out(n),pr(n)
      integer i,ii,in1,in2,ioutk,j,j1,jord,k,kk
      do 10 i = 1,n
        arp(i) = lenr(i) - 1
        cv(i) = 0
        iperm(i) = 0
   10 continue
      numnz = 0
      do 100 jord = 1,n
        j = jord
        pr(j) = -1
        do 70 k = 1,jord
          in1 = arp(j)
          if (in1.lt.0) go to 30
          in2 = ip(j) + lenr(j) - 1
          in1 = in2 - in1
          do 20 ii = in1,in2
            i = icn(ii)
            if (iperm(i).eq.0) go to 80
   20     continue
          arp(j) = -1
   30     continue
          out(j) = lenr(j) - 1
          do 60 kk = 1,jord
            in1 = out(j)
            if (in1.lt.0) go to 50
            in2 = ip(j) + lenr(j) - 1
            in1 = in2 - in1
            do 40 ii = in1,in2
              i = icn(ii)
              if (cv(i).eq.jord) go to 40
              j1 = j
              j = iperm(i)
              cv(i) = jord
              pr(j) = j1
              out(j1) = in2 - ii - 1
              go to 70
   40       continue
   50       continue
            j = pr(j)
            if (j.eq.-1) go to 100
   60     continue
   70   continue
   80   continue
        iperm(i) = j
        arp(j) = in2 - ii - 1
        numnz = numnz + 1
        do 90 k = 1,jord
          j = pr(j)
          if (j.eq.-1) go to 100
          ii = ip(j) + lenr(j) - out(j) - 2
          i = icn(ii)
          iperm(i) = j
   90   continue
  100 continue
      if (numnz.eq.n) return
      do 110 i = 1,n
        arp(i) = 0
  110 continue
      k = 0
      do 130 i = 1,n
        if (iperm(i).ne.0) go to 120
        k = k + 1
        out(k) = i
        go to 130
  120   continue
        j = iperm(i)
        arp(j) = i
  130 continue
      k = 0
      do 140 i = 1,n
        if (arp(i).ne.0) go to 140
        k = k + 1
        ioutk = out(k)
        iperm(ioutk) = i
  140 continue
      return
      end
! *******************************************************************
! COPYRIGHT (c) 1988 Hyprotech UK
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
! Original date 17 Feb 2005

! 17th February 2005 Version 1.0.0. Replacement for FD05.

      double precision function fd15ad(t)
!----------------------------------------------------------------
!----------------------------------------------------------------
      character t
      if ( t.eq.'E' ) then
         fd15ad = epsilon(1.0d0)
      else if ( t.eq.'T' ) then
         fd15ad = tiny(1.0d0)
      else if ( t.eq.'H' ) then
         fd15ad = huge(1.0d0)
      else if ( t.eq.'R' ) then
         fd15ad = dble(radix(1.0d0))
      else
         fd15ad = 0.0d0
      endif
      return
      end

