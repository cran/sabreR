      subroutine metis_nodend(n,iptr,irn,metftn,metopt,invprm,perm)
! Dummy routine that is called if MeTiS is not linked.
      integer n
      integer iptr(n+1),irn(*),metftn,metopt(8),invprm(n),perm(n)
      perm(1) = -1
      return
      end

