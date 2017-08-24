*
      subroutine pres(u,v,w,p,Ub,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,a1(2048),a2(2048)
     >,b1(2048),b2(2048)
     >,bp(2048),dp(2048),ep(2048)
      common
     >/dim/Xmax,epsr,dsym
     >/dimx/hx,Im,lx
     >/dimr/rt1(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
     >/rlx/rlx(2048)
     >/rlt/rlt(256)
     >/pry/apy(128),bpy(128),cpy(128)
*
      c0=0.d0
      Im2=Im/2
      cik=4.d0/(Im2*Km)
*
* Boundary conditions
      do k=1,Km
        do j=1,Jm
          u(0,j,k)=u(Im,j,k)
        end do
      end do
      do j=1,Jm
        do i=1,Im
          w(i,j,0)=c0
          w(i,j,Km)=c0
        end do
      end do
      do k=1,Km
        do i=1,Im
          v(i,Jm,k)=c0
          v(i,0,k)=c0
        end do
      end do
*
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            p(i,j,k)=d
          end do
        end do
      end do
*
*   FFT in tt-direction
      do j=1,Jm
        do i=1,Im 
          do k=1,Km
            a1(k)=p(i,j,k)
          end do
          call ftc05d(a1,b1,lt)
          do k=1,Km
            p(i,j,k)=b1(k)
          end do
        end do
      end do
*
*  FFT in x-direction
      do j=1,Jm
        do k=1,Km
          do i=1,Im
            b1(i)=p(i,j,k)
          end do 
          do i=1,Im2
            i1=Im+1-i
            a1(i)=0.5d0*(b1(i)+b1(i1))
            a2(i)=0.5d0*(b1(i)-b1(i1))
          end do
          call ftc05d(a1,b1,lx-1)
          call fts05d(a2,b2,lx-1)
          do i=1,Im2
            i1=Im+1-i
            a1(i)=b1(i)*cik
            a1(i1)=b2(i)*cik
          end do
          do i=1,Im
            p(i,j,k)=a1(i)
          end do           
        end do
      end do
*
*   Solution in wall-normal coordinate
      do k=1,Km
        do i=1,Im
          do j=1,Jm
            bp(j)=bpy(j)-rlx(i)-rlt(k)
            dp(j)=p(i,j,k)
          end do
          bp(1)=bp(1)+apy(1)
          bp(Jm)=bp(Jm)+cpy(Jm)
          Jm1=Jm
          ep(Jm)=0.d0
          if(rlx(i).eq.0.d0.and.rlt(k).eq.0.d0)Jm1=Jm-1
          call prog3(apy,bp,cpy,dp,ep,Jm1)
          do j=1,Jm
            p(i,j,k)=ep(j)
          end do
        end do 
      end do
*
*  Inverse FFT in x-direction
      do j=1,Jm
        do k=1,Km
          do i=1,Im
            b1(i)=p(i,j,k)
          end do 
          do i=1,Im2
            i1=Im+1-i
            a1(i)=b1(i)
            a2(i)=b1(i1)
          end do
          call ftc05b(a1,b1,lx-1)
          call fts05b(a2,b2,lx-1)
          do i=1,Im2
            i1=Im+1-i
            a1(i)=b1(i)+b2(i)
            a1(i1)=b1(i)-b2(i)
          end do
          do i=1,Im
            p(i,j,k)=a1(i)
          end do           
        end do
      end do
*
*   Inverse FFT in tt-direction
      do j=1,Jm
        do i=1,Im 
          do k=1,Km
            a1(k)=p(i,j,k)
          end do
          call ftc05b(a1,b1,lt)
          do k=1,Km
            p(i,j,k)=b1(k)
          end do
        end do
      end do
*
* Boundary conditions
      do k=1,Km
        do j=1,Jm
          p(Im+1,j,k)=p(1,j,k)
        end do
      end do
*
*  Mean pressure gradient
      ss=0.d0
      su=0.d0
      do j=1,Jm
        ss=ss+yt1(j)
        ssu=0.d0
        do k=1,Km
          do i=1,Im
            ssu=ssu+u(i,j,k)
          end do
        end do
        su=su+ssu*yt1(j)
      end do
      Dpp=Ub-su/(Im*Km*ss)
      p(0,0,0)=Dpp
*
      call gradp(u,v,w,p,Imax,Jmax)
      return
      end
*
      subroutine gradp(u,v,w,p,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,lx
     >/dimr/rt1(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
      Dp=p(0,0,0)
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k)=u(i,j,k)-(p(i+1,j,k)-p(i,j,k))/hx+Dp
          end do
        end do
      end do
*
      do k=1,Km
        do j=1,Jm-1
          do i=1,Im
            v(i,j,k)=v(i,j,k)-(p(i,j+1,k)-p(i,j,k))/rt1(j)
          end do
        end do
      end do
*
      do k=1,Km-1
        do j=1,Jm
          do i=1,Im
            w(i,j,k)=w(i,j,k)-(p(i,j,k+1)-p(i,j,k))/ht
          end do
        end do
      end do
      return
      end
*
      subroutine div(i,j,k,u,v,w,d,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,lx
     >/dimr/rt1(0:129),yt1(0:129),hr,Jm
     >/dimt/ht,Km,lt
      d=(u(i,j,k)-u(i-1,j,k))/hx
     > +(v(i,j,k)-v(i,j-1,k))/yt1(j)
     > +(w(i,j,k)-w(i,j,k-1))/ht
      return
      end
