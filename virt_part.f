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
c     integer, intent(out):: nvirt

      double precision hsml(maxn),mass(maxn),x(dim,maxn),vx(dim,maxn),
     &                 rho(maxn), u(maxn), p(maxn)
      integer i, j, d, im,mp,np, qp,scale_k, nnp, nwall, ogn(maxn)
      double precision xl, dx, v_gate, tiny, b, a, y1, y2
c      common nvirt
      real corner(2,4)

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
       mp = 65
       
       np = 32
       nnp = 10
       qp = int(nnp/a)+1
      else        
        np = mmp
      endif

	xl = x_maxgeom-x_mingeom
	dx = xl / mmp
c  speed of the gate(dambreak)/ speed of the top(cavityflow)        
      v_gate = 1.5
c      h = hsml(1)
c   coordinates of the corners
      corner(:,1)=(/x_mingeom,y_mingeom/)
      corner(:,2)=(/x_maxgeom,y_mingeom/)
c      print *,corner
      if(.not.dynamic.and..not.dummy) then

c     Monaghan type virtual particle on the Upper side
      
c        do i = 1, mp+1
c   	  nvirt = nvirt + 1
c	  x(1, ntotal + nvirt) = (i-1)*dx
c          x(2, ntotal + nvirt) = xl
c          vx(1, ntotal + nvirt) = 0.
c	  vx(2, ntotal + nvirt) = 0.
c        enddo
     

c     Monaghan type virtual particle on the Lower side

c        do i = 1, mp+1
c   	  nvirt = nvirt + 1
cc	   x(1, ntotal + nvirt) = x_mingeom+(i-1)*dx
c          x(2, ntotal + nvirt) = y_mingeom  
c        enddo

c     Monaghan type virtual particle on the Left side

c        do i = 1, np
c   	  nvirt = nvirt + 1
c 	  x(1, ntotal + nvirt) = x_mingeom 
c         x(2, ntotal + nvirt) = y_mingeom+i*dx
c        enddo

c     Monaghan type virtual particle on the Right side

c       do i = 1, mp-1
c    	  nvirt = nvirt + 1
c	  x(1, ntotal + nvirt) = xl 
c         x(2, ntotal + nvirt) = i*dx  
c        enddo
c    Monaghan type virtual particle as obsatacle 
c       do i = 1, 20
c          nvirt = nvirt + 1
c           x(1,ntotal+nvirt) = x_mingeom+(np+i)*dx
c          x(2,ntotal+nvirt) = y_mingeom+i*dx
c        enddo   
        nwall = nvirt
c	do i = 1, nvirt
c         vx(1, ntotal + i) = 0.
c	  vx(2, ntotal + i) = 0.
c	  if(itimestep.eq.1)rho (ntotal + i) = 1000.
c	     mass(ntotal + i) = rho (ntotal + i) * dx * dx
c          if(itimestep.eq.1) p(ntotal + i) = 0. 
c	   u(ntotal + i) = 357.1
c	   itype(ntotal + i) = 0
c	    hsml(ntotal + i) = 1.3*dx
c        enddo
      else
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
            x(2, ntotal+nvirt) = y_maxgeom-(i-1)*dx-dx/2
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
      do i = 1, ntotal+nwall
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
c    downside
           if ((x(2,i).gt.y_mingeom).and.
     &    (x(2,i).lt.y_mingeom+scale_k*hsml(i)))  then
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
c  mirror particles in corner

           do j=1,2
c              print *,corner(:,j)
              if(abs(x(1,i)-corner(1,j)).lt.scale_k*hsml(i).and.
c     &         abs(x(1,i)-corner(1,j)).ne.0.and.
     &         x(2,i)-corner(2,j).lt.scale_k*hsml(i)) then
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

         do i=1,nvirt
            itype = -2
            hsml(i+ntotal)=1.3*dx
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

c      print *,"inside virt_part"
c      print *,x(2,ntotal+1)

      end

      
