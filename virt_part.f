      subroutine virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,
     &           rho,u,p,itype, nwall) 

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

      implicit none
      include 'param.inc'
    
      
      integer itimestep, ntotal, itype(maxn),nvirt
c     integer, intent(out):: nvirt

      double precision hsml(maxn),mass(maxn),x(dim,maxn),vx(dim,maxn),
     &                 rho(maxn), u(maxn), p(maxn)
      integer i, j, d, im, mp, scale_k, line1(2), line2(2), nwall
      double precision xl, dx, v_inf, tiny, beta
c      common nvirt
      real corner(2,4)
      corner=0.
      line1=(/2,4/)
      line2=(/3,4/)

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
          read(3,*)im, itype(i), hsml(i)                            
        enddo  
        close(1)
        close(2) 
        close(3) 
      
	else 
       
	nvirt = 0
      mp = 39
      beta = 315.e0
	xl = 3.9e-1
	dx = xl / mp
	v_inf = 1.e-3
      tiny = 0
c   coordinate of the corners
      corner(1,line1)=xl
      corner(2,line2)=xl
     
      if(.not.dynamic) then

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
c	  x(1, ntotal + nvirt) = (i-1)*dx
c          x(2, ntotal + nvirt) = 0.  
c        enddo

c     Monaghan type virtual particle on the Left side

c        do i = 1, mp-1
c   	  nvirt = nvirt + 1
c	  x(1, ntotal + nvirt) = 0. 
c         x(2, ntotal + nvirt) = i*dx
c        enddo

c     Monaghan type virtual particle on the Right side

c       do i = 1, mp-1
c    	  nvirt = nvirt + 1
c	  x(1, ntotal + nvirt) = xl 
c         x(2, ntotal + nvirt) = i*dx  
c        enddo

        nwall = nvirt
c	do i = 1, nvirt
c         vx(1, ntotal + i) = 0.
c	  vx(2, ntotal + i) = 0.
c	  if(itimestep.eq.1)rho (ntotal + i) = 1000.
c	  mass(ntotal + i) = rho (ntotal + i) * dx * dx
c        if(itimestep.eq.1) p(ntotal + i) = 9.8*1000*(1e-3-x(2,ntotal+i)) 
c	  u(ntotal + i) = 357.1
c	   itype(ntotal + i) = 0
c	    hsml(ntotal + i) = 1.3*dx
c        enddo
      else
c--- staggered grid on the boundary, up-right-down-left
        
       do i = 1, mp
           nvirt = nvirt + 1
	      x(1, ntotal + nvirt) = (i-1)*dx
            x(2, ntotal + nvirt) = xl
        enddo
     
       do i = 1, mp
   	     nvirt = nvirt + 1
	     x(1, ntotal + nvirt) = xl
           x(2, ntotal + nvirt) = xl-(i-1)*dx
         enddo
     
       do i = 1, mp
   	    nvirt = nvirt + 1
	    x(1, ntotal + nvirt) = xl-(i-1)*dx
          x(2, ntotal + nvirt) = 0.
        enddo
      
       do i = 1, mp
   	    nvirt = nvirt + 1
	    x(1, ntotal + nvirt) = 0.
          x(2, ntotal + nvirt) = (i-1)*dx
         enddo
     
c     add four particles in the corner
c      nvirt = nvirt +4
c      x(1,ntotal + nvirt - 3) = -dx
c      x(2,ntotal + nvirt - 3) = xl+dx/2
c       x(1,ntotal + nvirt - 2) = xl+dx/2
c      x(2,ntotal + nvirt - 2) = xl+dx
c       x(1,ntotal + nvirt - 1) = xl+dx
c      x(2,ntotal + nvirt - 1) = -dx/2
c       x(1,ntotal + nvirt) = -dx/2
c      x(2,ntotal + nvirt) = -dx

      do i = 1, nvirt
        vx(1, ntotal + i) = 0.
	  vx(2, ntotal +i) = 0.
c        if (i.le.mp+1)  vx(1, ntotal +i) = v_inf   
c        vx(1,ntotal + nvirt - 2) = v_inf
	 if(itimestep.eq.1)  rho(ntotal + i) = 1000.
	  mass(ntotal + i) = rho (ntotal + i) * dx * dx
       if(itimestep.eq.1)  p(ntotal + i)= 0.
	  u(ntotal + i) = 357.1
	  itype(ntotal + i) = -2
	  hsml(ntotal + i) = 1.3*dx
        enddo
      endif
c     fixed ghost particles around the cavity, including particles in corner
      if (mirror) then
      do i = 1, ntotal+nwall
c    upside
c         if ((x(2,i).ge.xl-scale_k*hsml(i)).and.(x(2,i).lt.xl))  then
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
          if ((x(1,i).ge.xl-scale_k*hsml(i)).and.(x(1,i).lt.xl)) then
           nvirt=nvirt+1 
           x(1, ntotal + nvirt)=2*xl-x(1,i)
           x(2, ntotal + nvirt)=x(2,i)
           vx(1, ntotal + nvirt)=-vx(1,i)
            vx(2, ntotal + nvirt)=vx(2,i)
            rho(ntotal + nvirt)=rho(i)
            mass(ntotal + nvirt)=mass(i)
            p(ntotal + nvirt)= p(i)
            u(ntotal + nvirt)=u(i)
           endif
c    downside
           if ((x(2,i).gt.0).and.(x(2,i).le.scale_k*hsml(i)))  then
           nvirt=nvirt+1
           x(1, ntotal + nvirt) = x(1,i)
           x(2, ntotal + nvirt) = -x(2,i)
           vx(1, ntotal + nvirt) = vx(1,i)
           vx(2, ntotal + nvirt) = -vx(2,i)
           
           mass(ntotal + nvirt) = mass(i)
           p(ntotal + nvirt) = p(i)+2*9.8*1000*x(2,i)
           rho(ntotal + nvirt) =1000*(p(ntotal+nvirt)/beta+1)**(1/7)
           u(ntotal + nvirt) = u(i)
           endif
c   leftside
          if ((x(1,i).gt.0).and.(x(1,i).le.scale_k*hsml(i))) then
           nvirt=nvirt+1
           
           x(1, ntotal + nvirt) = -x(1,i)
           x(2, ntotal + nvirt) = x(2,i)
           vx(1, ntotal + nvirt) = -vx(1,i)
           vx(2, ntotal + nvirt) = vx(2,i)
           rho(ntotal + nvirt) = rho(i)
           mass(ntotal + nvirt) = mass(i)
           p(ntotal + nvirt) = p(i)
           u(ntotal + nvirt) = u(i)
           endif


           do j=1,2
              if(abs(x(1,i)-corner(1,j)).le.scale_k*hsml(i).and.
     &         abs(x(1,i)-corner(1,j)).ne.0.and.
     &         abs(x(2,i)-corner(2,j)).le.scale_k*hsml(i)) then
                nvirt=nvirt+1
                do d=1,dim
                   x(d,ntotal+nvirt)=2*corner(d,j)-x(d,i)
                   vx(d,ntotal+nvirt)=-vx(d,i)
                enddo
           
               mass(ntotal + nvirt)=mass(i)
               p(ntotal + nvirt)= p(i)+2*9.8*1000*x(2,i)
               rho(ntotal + nvirt) =1000*(p(ntotal+nvirt)/beta+1)**(1/7)
               u(ntotal + nvirt)=u(i)
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
        
  
      if (mod(itimestep,save_step).eq.0) then
        open(4,file="../data/xv_vp.dat")
        open(5,file="../data/state_vp.dat")
        open(6,file="../data/other_vp.dat")            
        write(4,*) nvirt
        do i = ntotal + 1, ntotal + nvirt         
          write(4,1004) i, (x(d, i), d=1,dim), (vx(d, i), d = 1, dim)              
          write(5,1005) i, mass(i), rho(i), p(i), u(i)
          write(6,1006) i, itype(i), hsml(i)                               
        enddo       
1004    format(1x, I6, 4(2x, e14.8))
1005    format(1x, I6, 4(2x, e14.8)) 
1006    format(1x, I6, 2x, I4, 2x, e14.8)
        close(4)
        close(5) 
        close(6) 
      endif 

      if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
c        print *,' >> Statistics: Virtual boundary particles:'
         print *,'          Number of virtual particles:',Nvirt
        endif     
      endif

      print *,"inside virt_part"
      print *,x(2,ntotal+1)

      end

      