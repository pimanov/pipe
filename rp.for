*
      subroutine rp(t,u,v,w,ut,vt,wt,ox,or,ot,buf,Imax,Jmax)
*
      t=t
      call bc_om(u,v,w,ox,or,ot,buf,Imax,Jmax)
      call visc(ox,or,ot,ut,vt,wt,Imax,Jmax)
      call add_nl(u,v,w,ox,or,ot,ut,vt,wt,Imax,Jmax)
*
      return
      end
