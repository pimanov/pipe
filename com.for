*
      subroutine com
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      common
     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt1(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/rlx/rlx(2048)
     >/rlt/rlt(256)
     >/pry/apy(128),bpy(128),cpy(128)
     >/prr/apr(128),bpr(128),cpr(128)
     >/prv/apv(128),bpv(128),cpv(128)
     >/prw/apw(128),bpw(128),cpw(128)
     >/servst/iserv
     >/pi/pi
     >/set/set0,set1,aset,bset,iset
     >/proc/Np,Npm
*
      iserv=0
      one=1.d0
      pi=4.d0*atan(one)
*
* dimt
      Km=2**lt
      ht=2.4/(dsym*Km)
* dimx
      Imm=2**lx
      Im=Imm/Npm
      hx=Xmax/Imm
* dimr
      iset=0
      set0=1.d0
      set1=epsr
      hr=1.d0/Jm
      do j=0,Jm+1
        ro=j*hr
        rt1(j)=rrt(ro,1)*hr
        ro=(j-0.5d0)*hr
        yt1(j)=rrt(ro,1)*hr
      end do
      do j=1,Jm
        c=1.d0/yt1(j)
        apy(j)=c/rt1(j-1)
        cpy(j)=c/rt1(j)
        bpy(j)=-apy(j)-cpy(j)
      end do
      do j=1,Jm-1
        c=1.d0/rt1(j)
        apr(j)=c/yt1(j)
        cpr(j)=c/yt1(j+1)
        bpr(j)=-apr(j)-cpr(j)
      end do
      do j=1,Jm-1
        c=1.d0/rt1(j)
        apv(j)=c/yt1(j)
        cpv(j)=c/yt1(j+1)
        bpv(j)=-apv(j)-cpv(j)
      end do
      j=1
        c=1.d0/yt1(j)
        apw(j)=0.d0
        cpw(j)=c/rt1(j)
        bpw(j)=-apw(j)-cpw(j)
      do j=2,Jm
        c=1.d0/yt1(j)
        apw(j)=c/rt1(j-1)
        cpw(j)=c/rt1(j)
        bpw(j)=-apw(j)-cpw(j)
      end do
* rl
      do i=0,Imm/2
        rlx(i+1)=(2.d0/hx*sin(i*pi/Imm))**2
      end do
      do i=1,Imm/2-1
        i1=Imm+1-i
        rlx(i1)=rlx(i+1)
      end do
      do k=0,Km-1
        rlt(k+1)=(2.d0/ht*sin(0.5d0*k*pi/Km))**2
      end do
      return
      end
