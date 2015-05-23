*
      subroutine servis(t,u,v,w,ox,or,ot,p,ii,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
      common
*     >/dim/Xmax,epsr
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/servst/iserv
*
      a=t
      a=u(0,0,0)
      a=v(0,0,0)
      a=w(0,0,0)
      a=ox(0,0,0)
      a=or(0,0,0)
      a=ot(0,0,0)
      a=p(0,0,0)
      i=ii
      return
      end
