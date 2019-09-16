      subroutine virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,
     &           rho,u,p,itype, nwall,mother) 

c----------------------------------------------------------------------
c   Subroutine to determine the information of virtual particles
c   Here only the Monaghan type virtual particles for the 2D shear
c   cavity driven problem are generated.
c     itimestep : Current time step                                 [in]
c     ntotal : Number of particles                                  [in]
c     nvirt  : Number of virtual particles                         [out]
c     hsml   : Smoothing Length                                 [in|out]
c     mass   : Particle masses                                  [in|out]
c     x      : Coordinates of all particles                     [in|out]
c     vx     : Velocities of all particles                      [in|out]
c     rho    : Density                                          [in|out]
c     u      : internal energy                                  [in|out]
c     itype   : type of particles                               [in|out]
c     mother  : indicating which particle it was generated from  [out]
c     ogn     : number of over generation                        [out]

      implicit none
      include 'param.inc'

      integer itimestep, ntotal, itype(maxn), nvirt, mother(maxn)
      double precision hsml(maxn),mass(maxn),x(dim,maxn),vx(dim,maxn),
     &                 rho(maxn), u(maxn), p(maxn)
      integer i, j, d, im,mp,np, qp,scale_k, nnp, nwall, flag
      double precision xl, dx, v_gate, tiny, b, a, y1, y2, dps
c      common nvirt
      real corner(2,5)

      if (skf.eq.1) then 
        scale_k = 2 
      else if (skf.eq.2) then 
        scale_k = 3 
      else if (skf.eq.3) then 
        scale_k = 3 
      endif 
     
      if (vp_input) then          
                        
        open(1,file="../data/xv_vp.dat")
        open(2,file="../data/state_vp.dat")
        open(3,file="../data/other_vp.dat")            
        read(1,*) nvirt
        do j = 1, nvirt   
          i = ntotal + j      
          read(1,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)                     
          read(2,*)im, mass(i), rho(i), p(i), u(i)        
          read(3,*)im, itype(i), hsml(i), mother(i)                            
        enddo  
        close(1)
        close(2) 
        close(3) 
      
	else 
c   in this case, the computing domain is defaultly set as square
       
	nvirt = 0
      if(geometry) then
       a = tan(pi/3)
c       print *,a
       mp = 64 
       np = 32
       nnp = 10
       qp = int(nnp/a)+1
      else        
        np = mmp
      endif

	xl = x_maxgeom-x_mingeom
	dx = 2.e-2
c  speed of the gate(dambreak)/ speed of the top(cavityflow)        
      v_gate = 1.5
c      h = hsml(1)
c   coordinates of the corners
      corner(1,1:4)=(/x_mingeom,x_maxgeom,xl/2-nnp*dx/a,xl/2+nnp*dx/a/)
      corner(:,5)=(/xl/2,nnp*dx/)
c      corner(:,2)=(/xl/2-nnp*dx/a,y_mingeom/)
c      print *,corner
      if(.not.dynamic.and..not.dummy.and..not.mirror) then
c     repulsive boundary   --fixed solid particle
c     Monaghan type virtual particle on the Upper side
      
c        do i = 1, mp+1
c   	  nvirt = nvirt + 1
c	  x(1, ntotal + nvirt) = (i-1)*dx
c          x(2, ntotal + nvirt) = xl
c          vx(1, ntotal + nvirt) = 0.
c	  vx(2, ntotal + nvirt) = 0.
c        enddo
     
c     Monaghan type virtual particle on the Lower side

        do i = 1, 2*mp+1
          if (i.gt.mp-qp*2+1.and.i.le.mp) then 
           nvirt = nvirt + 1
           x(1,ntotal+nvirt)=x_mingeom+(i-1)*dx/2
           x(2,ntotal+nvirt)=a*x(1,ntotal+nvirt)-a*xl/2+nnp*dx
          else if(i.gt.mp.and.i.le.mp+qp*2) then
            nvirt = nvirt + 1
            x(1,ntotal+nvirt)=x_mingeom+(i-1)*dx/2
            x(2,ntotal+nvirt)=-a*x(1,ntotal+nvirt)+a*xl/2+nnp*dx
          else        
   	   nvirt = nvirt + 1
           x(1, ntotal + nvirt) = x_mingeom+(i-1)*dx/2
           x(2, ntotal + nvirt) = y_mingeom  
          endif
        enddo

c     Monaghan type virtual particle on the Left side

        do i = 1, np*2
   	  nvirt = nvirt + 1
 	  x(1, ntotal + nvirt) = x_mingeom 
         x(2, ntotal + nvirt) = y_mingeom+i*dx/2
        enddo

c     Monaghan type virtual particle on the Right side

       do i = 1, np*2
    	  nvirt = nvirt + 1
	  x(1, ntotal + nvirt) = x_maxgeom 
          x(2, ntotal + nvirt) = y_mingeom+i*dx/2  
       enddo
c    Monaghan type virtual particle as obsatacle 
c       do i = 1, 20
c          nvirt = nvirt + 1
c           x(1,ntotal+nvirt) = x_mingeom+(np+i)*dx
c          x(2,ntotal+nvirt) = y_mingeom+i*dx
c        enddo   
        nwall = nvirt
	do i = 1, nvirt
         vx(1, ntotal + i) = 0.
	  vx(2, ntotal + i) = 0.
	  if(itimestep.eq.1)rho (ntotal + i) = 1000.
	   mass(ntotal + i) = rho (ntotal + i) * dx * dx
          if(itimestep.eq.1) p(ntotal + i) = 0. 
	   u(ntotal + i) = 357.1
	   itype(ntotal + i) = 0
	    hsml(ntotal + i) = 1.3*dx/2
        enddo
      endif

      if(dynamic.or.dummy)then
c--- staggered grid on the boundary, left-down-right
c    upside        
c       do i = 1, mp
c           nvirt = nvirt + 1
c	      x(1, ntotal + nvirt) = (i-1)*dx
c            x(2, ntotal + nvirt) = xl
c        enddo
      if (static)then
         do i = 1,np+3
            do j = 1,3
            nvirt = nvirt + 1
            x(1, ntotal+nvirt) = x_maxgeom+(j-1)*dx+dx/2
            x(2, ntotal+nvirt) = y_mingeom+(i-3)*dx-dx/2
            enddo
         enddo

         do i = 1,mp+3
            do j = 1,3
            nvirt = nvirt + 1
           x(1, ntotal+nvirt) = x_maxgeom-i*dx+dx/2
            x(2, ntotal+nvirt) = y_mingeom-(j-1)*dx-dx/2
            enddo
         enddo

         do i = 1,np
            do j = 1,3
            nvirt = nvirt + 1
           x(1, ntotal+nvirt) = x_mingeom-(j-1)*dx-dx/2
            x(2, ntotal+nvirt) = y_mingeom+i*dx-dx/2
            enddo
         enddo
c    add obstacle
         if(geometry)then
c  on the left wall
c          do i = 1,qp
c            do j =1,i     
c             nvirt = nvirt+1
c              x(1,ntotal+nvirt) = x_mingeom+(qp-i)*dx+dx/2
c              x(2,ntotal+nvirt) = y_mingeom+(nnp-j)*dx+dx/2
c            enddo
c          enddo

c at the bottom
         do i = np+1-qp,np+1+qp
           do j = 1, nnp
            y1 = a*(i*dx-dx/2)+nnp*dx-a*xl/2
            y2 = a*(dx/2-i*dx)+nnp*dx+a*xl/2
             if(j*dx-dx/2.lt.y1.and.j*dx-dx/2.lt.y2) then
                nvirt = nvirt + 1
               x(1,ntotal+nvirt) = x_mingeom+i*dx-dx/2
               x(2,ntotal+nvirt) = y_mingeom+j*dx-dx/2
             endif 
            enddo
          enddo
        endif

        do i =1,nvirt
          vx(1, ntotal + i) = 0.
          vx(2, ntotal +i) = 0.
          rho(ntotal + i) = 1000.
          p(ntotal + i) = 0.
          mass(ntotal + i) = rho (ntotal + i) * dx * dx
c  if(itimestep.eq.1)  p(ntotal + i)= 1000*9.8*(xl-x(2,ntotal+i))
          u(ntotal + i) = 357.1
          itype(ntotal + i) = -2
          hsml(ntotal + i) = 1.3*dx
        enddo

      endif

      if(dambreak)then     
       do i = 1, nnp
          do j = 1,2
   	   nvirt = nvirt + 1
	   x(1, ntotal + nvirt) = x_maxgeom+(j-1)*dx/2
           x(2, ntotal + nvirt) = y_maxgeom-(i-1)*dx-(j-1)*dx/2
           p(ntotal+nvirt) = 0
           if(x(2,ntotal+nvirt).le.bedheight) then
             p(ntotal+nvirt)=9.8*1000*(bedheight-x(2,ntotal+nvirt))
           endif
          enddo
       enddo
     
       do i = 1, mp
          do j = 1,2
   	    nvirt = nvirt + 1
	    x(1, ntotal + nvirt) = x_maxgeom-(i-1)*dx-(j-1)*dx/2
            x(2, ntotal + nvirt) = y_mingeom-(j-1)*dx/2
            if (x(1,ntotal+nvirt).le.damlength)then
             p(ntotal+nvirt) = ((j-1)*dx/2+damheight)*9.8*1000
            else
             p(ntotal+nvirt) = ((j-1)*dx/2+bedheight)*9.8*1000
            endif
           enddo 
       enddo
      
      do i = 1, nnp
         do j =1,2
   	    nvirt = nvirt + 1
	    x(1, ntotal + nvirt) = x_mingeom-(j-1)*dx/2
           x(2, ntotal + nvirt) = y_mingeom+(i-1)*dx+(j-1)*dx/2
           p(ntotal+nvirt) = 0
           if(x(2,ntotal+nvirt).le.damheight) then
             p(ntotal+nvirt)=9.8*1000*(damheight-x(2,ntotal+nvirt))
           endif
         enddo
      enddo
      if (gate) then
          do i = 1,2*np
            nvirt = nvirt + 1
            x(1,ntotal+nvirt) = damlength
            x(2,ntotal+nvirt) = y_mingeom+i*dx/2
            p(ntotal+nvirt) = 0
             if(x(2,ntotal+nvirt).le.damheight) then
             p(ntotal+nvirt)=9.8*1000*(damheight-x(2,ntotal+nvirt))
            endif
         enddo
      endif

      do i = 1, nvirt
        vx(1, ntotal + i) = 0.
	  vx(2, ntotal +i) = 0.
         if (i.gt.nvirt-np*2)  vx(2, ntotal +i) = v_gate  
c        vx(1,ntotal + nvirt - 2) = v_inf
	  rho(ntotal + i) = 1000.
	  mass(ntotal + i) = rho (ntotal + i) * dx * dx
c       if(itimestep.eq.1)  p(ntotal + i)= 1000*9.8*(xl-x(2,ntotal+i))
	  u(ntotal + i) = 357.1
	  itype(ntotal + i) = -2
	  hsml(ntotal + i) = 1.3*dx
        enddo
       endif

      endif
c     fixed ghost particles around the cavity, including particles in corner
      if (mirror) then
        do i = 1, ntotal
c    upside
c         if((x(2,i).gt.y_maxgeom-scale_k*hsml(i)).and.
c    &  (x(2,i).lt.y_maxgeom))  then
c           nvirt=nvirt+1
c           x(1, ntotal + nvirt)=x(1,i)
c          x(2, ntotal + nvirt)=2*xl-x(2,i)
c          vx(1, ntotal + nvirt)=vx(1,i)
c          vx(2, ntotal + nvirt)=-vx(2,i)
c          rho(ntotal + nvirt)=rho(i)
c           mass(ntotal + nvirt)=mass(i)
c          p(ntotal + nvirt)= p(i)-9.8*1000*(1e-3-x(2,i))
c          u(ntotal + nvirt)=u(i)
c          endif
c    rightside
      if ((x(1,i).gt.x_maxgeom-scale_k*hsml(i)).and.
     &  (x(1,i).lt.x_maxgeom)) then
           nvirt=nvirt+1 
           x(1, ntotal + nvirt)=2*x_maxgeom-x(1,i)
           x(2, ntotal + nvirt)=x(2,i)
           vx(1, ntotal + nvirt)=-vx(1,i)
            vx(2, ntotal + nvirt)=vx(2,i)
            rho(ntotal + nvirt)=rho(i)
            mass(ntotal + nvirt)=mass(i)
            p(ntotal + nvirt)=p(i)
            u(ntotal + nvirt)=u(i)
            mother(ntotal + nvirt)=i
           endif
c    downside  except the wedge
        if ((x(2,i).gt.y_mingeom).and.
     &    (x(2,i).lt.y_mingeom+scale_k*hsml(i)))then

          if ((x(2,i).lt.a*x(1,i)-a*xl/2-nnp*dx).or.
     &    (x(2,i).lt.-a*x(1,i)+a*xl/2-nnp*dx))  then
           nvirt=nvirt+1
           x(1, ntotal + nvirt) = x(1,i)
           x(2, ntotal + nvirt) = 2*y_mingeom-x(2,i)
           vx(1, ntotal + nvirt) = vx(1,i)
           vx(2, ntotal + nvirt) = -vx(2,i)
           mass(ntotal + nvirt) = mass(i)
             p(ntotal + nvirt) = p(i)
            rho(ntotal + nvirt) = rho(i)
           mother(ntotal + nvirt)=i
           u(ntotal + nvirt) = u(i)
           endif
       endif
c   leftside
          if ((x(1,i).gt.x_mingeom).and.
     &    (x(1,i).lt.x_mingeom+scale_k*hsml(i))) then
           nvirt=nvirt+1
           x(1, ntotal + nvirt) = 2*x_mingeom-x(1,i)
           x(2, ntotal + nvirt) = x(2,i)
           vx(1, ntotal + nvirt) = -vx(1,i)
           vx(2, ntotal + nvirt) = vx(2,i)
           rho(ntotal + nvirt) = rho(i)
           mass(ntotal + nvirt) = mass(i)
           p(ntotal + nvirt) = p(i)
           u(ntotal + nvirt) = u(i)
           mother(ntotal + nvirt) = i
           endif
        if(geometry)then
c wedge side
         if  ((x(2,i).gt.-a*x(1,i)+a*xl/2-nnp*dx).and.
c     &     (x(2,i).gt.-a*x(1,i)+a*xl/2-nnp*dx).and.
     &     (x(2,i).lt.nnp*dx).and.(x(1,i).le.xl/2)) then
               if (x(2,i).lt.a*x(1,i)-a*xl/2+nnp*dx+4*hsml(i))then 
                 nvirt = nvirt + 1
                 x(1,ntotal + nvirt)=((1-a**2)*x(1,i)+2*a*
     &          x(2,i)- 2*a*(-a*xl/2+nnp*dx))/(a**2+1)
                 x(2,ntotal + nvirt)=(2*a*x(1,i)+(a**2-1)*x(2,i)+
     &          2*(-a*xl/2+nnp*dx))/(a**2+1)

                 vx(1,ntotal + nvirt)=((sin(pi/3)+cos(pi/3))*
     &          vx(1,i)**2+ (-cos(pi/3)+sin(pi/3))*vx(1,i)*vx(2,i))/
     &           sqrt(vx(1,i)**2+vx(2,i)**2)            
                 vx(2,ntotal + nvirt)=((sin(pi/3)+cos(pi/3))*
     &          vx(1,i)*vx(2,i)+ (-cos(pi/3)+sin(pi/3))*vx(2,i)**2)/
     &           sqrt(vx(1,i)**2+vx(2,i)**2)
                 if(itimestep.eq.1)then
                   do d=1,dim
                    vx(d,ntotal+nvirt)=0
                   enddo
                 endif 
              rho(ntotal + nvirt) = rho(i)
c   over creation, in this case, twice
               if(x(2,i).gt.nnp*dx-scale_k*hsml(i))then
                  mass(ntotal+nvirt)=mass(i)/2
               else
                   mass(ntotal + nvirt) = mass(i)
               endif
               p(ntotal + nvirt) = p(i)
               u(ntotal + nvirt) = u(i)
               mother(ntotal + nvirt) = i
c               print *,ntotal+nvirt,mother(ntotal+nvirt)
              endif
           endif        
c             if  ((x(2,i).gt.a*x(1,i)-a*xl/2-nnp*dx).and.
              if  ((x(2,i).gt.a*x(1,i)-a*xl/2-nnp*dx).and.
     &     (x(2,i).lt.nnp*dx).and.(x(1,i).ge.xl/2)) then
               if (x(2,i).lt.-a*x(1,i)+a*xl/2+nnp*dx+4*hsml(i))then
                 nvirt = nvirt + 1
                 x(1,ntotal + nvirt)=((1-a**2)*x(1,i)-2*a*
     &          x(2,i)+2*a*(a*xl/2+nnp*dx))/(a**2+1)
                 x(2,ntotal + nvirt)=(-2*a*x(1,i)+(a**2-1)*x(2,i)+
     &          2*(a*xl/2+nnp*dx))/(a**2+1)
                 vx(1,ntotal + nvirt)=((-sin(pi/3)-cos(pi/3))*
     &          vx(1,i)**2+ (-cos(pi/3)+sin(pi/3))*vx(1,i)*vx(2,i))/
     &           sqrt(vx(1,i)**2+vx(2,i)**2)
                 vx(2,ntotal + nvirt)=((-sin(pi/3)-cos(pi/3))*
     &          vx(1,i)*vx(2,i)+ (-cos(pi/3)+sin(pi/3))*vx(2,i)**2)/
     &           sqrt(vx(1,i)**2+vx(2,i)**2)
                  if(itimestep.eq.1)then
                   do d=1,dim
                    vx(d,ntotal+nvirt)=0
                   enddo
                 endif
              rho(ntotal + nvirt) = rho(i)
c   over creation, in this case, twice
               if(x(2,i).gt.nnp*dx-scale_k*hsml(i))then
                  mass(ntotal+nvirt)=mass(i)/2
               else
                   mass(ntotal + nvirt) = mass(i)
               endif
               p(ntotal + nvirt) = p(i)
               u(ntotal + nvirt) = u(i)
               mother(ntotal + nvirt) = i
c               print *,ntotal+nvirt,mother(ntotal+nvirt)
              endif
           endif
c   sharp corner (central immetry)
        dps = sqrt((x(1,i)-corner(1,5))**2+(x(2,i)-corner(2,5))**2)
        if ((x(2,i).gt.a*x(1,i)-a*xl/2+nnp*dx).and.
     &    (x(2,i).gt.-a*x(1,i)+a*xl/2+nnp*dx).and.
     &     (dps.lt.scale_k*hsml(i)))then
            nvirt = nvirt + 1
            do d = 1,dim
             x(d,ntotal+nvirt)=2*corner(d,5)-x(d,i)
             vx(d,ntotal+nvirt)=-vx(d,i)
            enddo
            mass(ntotal + nvirt)= mass(i)
            p(ntotal + nvirt)= p(i)
            rho(ntotal + nvirt) = rho(i)
            u(ntotal + nvirt) = u(i)
            mother(ntotal + nvirt) = i
          endif
        endif
c  mirror particles in corner including corners of the wedge

          do j=1,4
           dps = sqrt((x(1,i)-corner(1,j))**2+(x(2,i)-corner(2,j))**2)             
           if (dps.lt.scale_k*hsml(i)) then
                nvirt=nvirt+1
                do d=1,dim
                   x(d,ntotal+nvirt)=2*corner(d,j)-x(d,i)
                   vx(d,ntotal+nvirt)=-vx(d,i)
                enddo
              mass(ntotal + nvirt)=mass(i)
              p(ntotal + nvirt)= p(i)
              rho(ntotal + nvirt) =rho(i)
               u(ntotal + nvirt)=u(i)
               mother(ntotal + nvirt) = i
               endif
           enddo
         enddo

         do i=ntotal+1,ntotal+nvirt
            itype(i) = -2
            hsml(i)= 1.3*dx
           dps = sqrt((x(1,i)-corner(1,5))**2+(x(2,i)-corner(2,5))**2) 
            if (dps.lt.scale_k*hsml(i))then
               mass(i) = mass(mother(i))/3
            endif
         enddo
       endif
      endif   

c      if(itimestep.eq.1) then

      if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
c        print *,' >> Statistics: Virtual boundary particles:'
         print *,'   Number of virtual particles:',nvirt
        endif     
      endif

      end

      
