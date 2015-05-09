*
      function rrt(x,i)
      implicit real*8 (a-h,o-z)
      common/set/y0,y01,aset,bset,iset
      yp(a,b)=0.01d0*sqrt(0.3d0*(b+1))*
     > (a*(58.d0*b-21.d0-21.d0*b*b)+100.d0)
      y1p(a,b)=1.d0+a*(b-0.45d0*(b+1.d0)**2)
      dyda(a,b)=0.01d0*sqrt(0.3d0*(b+1))*(58.d0*b-21.d0-21.d0*b*b)
      dydb(a,b)=0.01d0*sqrt(0.3d0*(b+1))*
     > (0.5d0*(a*(58.d0*b-21.d0-21.d0*b*b)+100.d0)/(b+1)+a*(58.d0-42*b))
      dy1da(a,b)=b-0.45d0*(b+1.d0)**2
      dy1db(a,b)=a*(1.d0-0.9d0*(b+1.d0))
      if(iset.ne.0) goto 2
      a=0
      b=10.d0/3.d0*y0*y0-1.d0
      nit=0
1     continue
      c1=dyda(a,b)
      c2=dydb(a,b)
      c3=y0-yp(a,b)
      c4=dy1da(a,b)
      c5=dy1db(a,b)
      c6=y01-y1p(a,b)
      d=c1*c5-c2*c4
      da=(c3*c5-c2*c6)/d
      db=(c1*c6-c3*c4)/d
      nit=nit+1
      a=a+da
      b=b+db
      if(abs(da).gt.1.e-5.or.abs(db).gt.1.e-5) goto 1
      aset=a
      bset=b
      iset=1
2     continue
      x2=x*x
      res=-1.d0
      if(i.eq.0)res=x*(1.d0+aset*(1.d0-x2)*(bset-x2))
      if(i.eq.1)res=1.d0+aset*(bset-3.d0*x2*(1.d0+bset)+5.d0*x2*x2)
      rrt=res
      return
      end
