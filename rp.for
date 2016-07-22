*
      subroutine rp(t,u,v,w,ut,vt,wt,ox,or,ot,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      t=t
      call bc_om(u,v,w,ox,or,ot,Imax,Jmax)
      call visc(ox,or,ot,ut,vt,wt,Imax,Jmax)
      call add_nl(u,v,w,ox,or,ot,ut,vt,wt,Imax,Jmax)
      return
      end
