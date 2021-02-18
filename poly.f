        subroutine polyse(px,py,n,xx,yy,m,flag)
        implicit real*8 (a-h, o-z)
        double precision px(1:n+1),py(1:n+1),x,y,temp,esp
        double precision xx(1:m),yy(1:m)
        integer te,n,ip,sgn,is,m,flag(1:m)
        esp=1e-8
        
        do k=1,m
        x=xx(k)
        y=yy(k)
        ip=0
        te=-1
        do i=1,n
           k2 = i+1
           if(k2.gt.n) k2=1
c          write(*,*)py(i),y,py(i+1)
c          write(*,*)'sgn',sgn(y-py(i))*sgn(py(i+1)-y)
          is=sgn(y-py(i))*sgn(py(k2)-y) 
         if(is.ge.0) then
           if(is.gt.0) then
              temp=px(k2)*abs(y-py(i))+px(i)*abs(py(k2)-y)
     &              -abs(py(k2)-py(i))*x
c                write(*,*) 'temp=',temp
                if(temp.gt.0) then 
                    ip=ip+1
c                    write(*,*)i,'ip=',ip
                  elseif(temp.eq.0) then
                    te=0
                    goto 20
                  else
                endif
            elseif(y.eq.py(i).and.py(k2).gt.y)  then
               if (x-px(i).lt.0 ) then 
                     ip=ip+1
c                    write(*,*)i,'ip=',ip
               elseif (x-px(i).eq.0) then
                     te=0
                     goto 20
               endif
            elseif((y.eq.py(k2)).and.py(i).gt.y)  then
               if (x-px(k2).lt.0 ) then 
                     ip=ip+1
cc                     write(*,*)i,'ip=',ip
               elseif (x-px(k2).eq.0) then
                     te=0
                     goto 20
               endif
            elseif(y.eq.py(i).and.py(k2).lt.y)  then
               if (x-px(i).eq.0) then
                     te=0
                     goto 20
               endif
            elseif((y.eq.py(k2)).and.py(i).lt.y)  then
               if (x-px(k2).eq.0) then
                     te=0
                     goto 20
               endif
            elseif(y.eq.py(k2).and.py(i).eq.y)  then
                if(sgn(x-px(i))*sgn(px(k2)-x).gt.0) then
                     te=0
                     goto 20
                endif
             else
          endif
        endif
        enddo
cc        write(*,*)i,'ip=',ip
        if (ip-ip/2*2.eq.1) then 
           te=1
         else 
           te=-1
        endif
 20     continue
        flag(k)=te
        enddo
        return
        end

        integer function sgn(x)
        double precision x
          if (x .gt. 0.0d0) then 
              sgn=1
            elseif (x .lt. 0.0d0) then
              sgn=-1
            else 
              sgn=0
           endif
          return
         end

        function ycentroid(tx,ty, np)
        implicit real*8 (a-h, o-z)

        real *8 tx(np+1), ty(np+1)
        
        area=0
        ymoment=0

        do i=1, np
          k2 = i+1
          if(k2.gt.np)k2=1
          temp=tx(i)*ty(k2)-tx(k2)*ty(i)
          area=area+temp
          ymoment=ymoment+(ty(i)+ty(k2))*temp
        enddo
     
        area=area/2
C        write(*,*) area, ymoment
        ycentroid=ymoment/6d0/area
        return
        end

        subroutine calcularea(tx,ty, np,area)
        implicit real*8 (a-h, o-z)

        real *8 tx(np+1), ty(np+1)

        area=0d0

        do i=1, np
          k2 = i+1
          if(k2.gt.np)k2=1
          temp=tx(i)*ty(k2)-tx(k2)*ty(i)
          area=area+temp
        enddo

        area=abs(area)/2
C        write(*,*) area, ymoment
        return
        end

        
