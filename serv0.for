*
      subroutine servis(t,u,v,w,ox,or,ot,p,ii,Jmax)
      implicit real*8 (a-h,o-z)
      complex*16 u,v,w,ox,or,ot,p
      dimension
     > u(0:Jmax,0:*)
     >,v(0:Jmax,0:*)
     >,w(0:Jmax,0:*)
     >,ox(0:Jmax,0:*)
     >,or(0:Jmax,0:*)
     >,ot(0:Jmax,0:*)
     >,p(0:Jmax,0:*)
      common
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/servst/iserv
     >/alpha/alpha
*
      t=t
      u(0,0)=u(0,0)
      v(0,0)=v(0,0)
      w(0,0)=w(0,0)
      ox(0,0)=ox(0,0)
      or(0,0)=or(0,0)
      ot(0,0)=ot(0,0)
      p(0,0)=p(0,0)
      ii=ii
      return
      end
