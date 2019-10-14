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
     &                 rho(maxn), u(maxn), p(maxn),slope(maxn)
      integer i, j, d, im,mp,np, qp,scale_k, nnp, nwall, flag,ocn(maxn)
      double precision xl, dx, v_gate, tiny, b, a, y1, y2, dps,period
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
        nwall = 0
        do i = 1,ntotal
          ocn(i) =1
        enddo
c      if(geometry) then
       a = tan(pi/3)
c       print *,a
       mp = 32 
       np = 16
       nnp = 5
       qp = int(nnp/a)+1
       if(.not.geometry) qp = 0     
c        np = mmp
c      endif
      

	xl = x_maxgeom-x_mingeom
	dx = xl/mmp
        period = sqrt(3.)*pi/(2*nnp*dx)

c  speed of the gate(dambreak)/ speed of the top(cavityflow)        
      v_gate = 1.5
c      h = hsml(1)
c   coordinates of the corners
      corner(1,1:4)=(/x_mingeom,x_maxgeom,xl/2-qp*dx,xl/2+qp*dx/)
      corner(:,5)=(/xl/2,a*qp*dx/)
c      corner(:,2)=(/xl/2-nnp*dx/a,y_mingeom/)
c      print *,corner
      if(.not.dynamic.and..not.dummy) then
c     repulsive boundary   --fixed solid particle
c     Monaghan type virtual particle on the Upper side
      
c        do i = 1, mp+1
c   	  nvirt = nvirt + 1
c	  x(1, ntotal + nvirt) = (i-1)*dx
c          x(2, ntotal + nvirt) = xl
c          vx(1, ntotal + nvirt) = 0.
c	  vx(2, ntotal + nvirt) = 0.
c        enddo
       if(geom.eq.1.and.ex_force) then   
c     Monaghan type virtual particle on the Lower side

        do i = 1, 2*mp+1
          if ((i-1)*dx/2.le.xl/2-nnp*dx/a.or.(i-1)*dx/2.ge.xl/2+
     &   nnp*dx/a)then
         nvirt = nvirt + 1
           x(1, ntotal + nvirt) = x_mingeom+(i-1)*dx/2
           x(2, ntotal + nvirt) = y_mingeom
c          endif
c        enddo
c        do i = 1,2*mp+1
c         if((i-1)*dx/2.gt.xl/2-nnp*dx/a.and.(i-1)*dx/2.lt.xl/2+
c     &   nnp*dx/a) then
c           nvirt = nvirt + 1
c           x(1,ntotal+nvirt)=x_mingeom+(i-1)*dx/2
c           x(2,ntotal+nvirt)=nnp*dx*sin(period*(x(1,ntotal+nvirt)-
c     &     xl/2+nnp*dx/a))
c           slope(ntotal+nvirt)=nnp*dx*period*cos(period*
c     &     (x(1,ntotal+nvirt)-xl/2+nnp*dx/a))
c           dps = sqrt((x(1,ntotal+nvirt)-x(1,ntotal+nvirt-1))**2+
c     &      (x(2,ntotal+nvirt)-x(2,ntotal+nvirt-1))**2)
c           if(dps.gt.dx/2)then
c             do j =1,int(2*dps/dx)
c               ntotal = ntotal+1
c               x(1,ntotal+nvirt)=x(1,ntotal+nvirt-1)-
c    &          dps/(int(2*dps/dx)+1)
c              x(2,ntotal+nvirt)=nnp*dx*sin(period*(x(1,ntotal+nvirt)-
c     &     xl/2+nnp*dx/a))
c             slope(ntotal+nvirt)=nnp*dx*period*cos(period*
c     &     (x(1,ntotal+nvirt)-xl/2+nnp*dx/a))
c             enddo
c           endif
c          else
c           nvirt = nvirt + 1
c           x(1, ntotal + nvirt) = x_mingeom+(i-1)*dx/2
c           x(2, ntotal + nvirt) = y_mingeom

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
       do i = 1, qp*2
          nvirt = nvirt + 2
          x(1,ntotal+nvirt-1) = x_mingeom+(np-qp)*dx+i*dx/2
          x(1,ntotal+nvirt) = xl-x(1,ntotal+nvirt-1)
          x(2,ntotal+nvirt-1) = y_mingeom+a*i*dx/2
          x(2,ntotal+nvirt) = x(2,ntotal+nvirt-1)
       enddo

        nwall = nvirt
	do i = 1, nvirt
         vx(1, ntotal + i) = 0.
	  vx(2, ntotal + i) = 0.
	  if(itimestep.eq.1)rho (ntotal + i) = 1000.
	   mass(ntotal + i) = rho(ntotal+i)*dx*dx
          if(itimestep.eq.1) p(ntotal + i) = 0. 
	   u(ntotal + i) = 357.1
	   itype(ntotal + i) = 0
	    hsml(ntotal + i) = 1.3*dx
        enddo
       endif
      endif

      if(dynamic.or.dummy)then
c--- staggered grid on the boundary, left-down-right
c    upside        
c       do i = 1, mp
c           nvirt = nvirt + 1
c	      x(1, ntotal + nvirt) = (i-1)*dx
c            x(2, ntotal + nvirt) = xl
c        enddo
       call fix_particle(xl,dx,nnp,ntotal,nvirt,v_gate,x,vx,p,rho,mass,
     &     u,hsml,itype)
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
            mother(ntotal + nvirt)=i
           endif
c    downside  except the wedge
        if ((x(2,i).gt.y_mingeom).and.
     &    (x(2,i).lt.y_mingeom+scale_k*hsml(i)))then

c          if ((x(2,i).lt.a*x(1,i)-a*xl/2-a*qp*dx).or.
c     &    (x(2,i).lt.-a*x(1,i)+a*xl/2-a*qp*dx))  then
           nvirt=nvirt+1
           x(1, ntotal + nvirt) = x(1,i)
           x(2, ntotal + nvirt) = 2*y_mingeom-x(2,i)
           vx(1, ntotal + nvirt) = vx(1,i)
           vx(2, ntotal + nvirt) = -vx(2,i)
           mother(ntotal + nvirt)=i
c           endif
       endif
c   leftside
          if ((x(1,i).gt.x_mingeom).and.
     &    (x(1,i).lt.x_mingeom+scale_k*hsml(i))) then
           nvirt=nvirt+1
           x(1, ntotal + nvirt) = 2*x_mingeom-x(1,i)
           x(2, ntotal + nvirt) = x(2,i)
           vx(1, ntotal + nvirt) = -vx(1,i)
           vx(2, ntotal + nvirt) = vx(2,i)
           mother(ntotal + nvirt) = i
           endif
        do j = 1,2
           dps = sqrt((x(1,i)-corner(1,j))**2+(x(2,i)-corner(2,j))**2)
           if(dps.lt.scale_k*hsml(i).and.dps.gt.1e-9)then
               nvirt = nvirt+1
               do d =1,dim
                 x(d,ntotal+nvirt)=2*corner(d,j)-x(d,i)
                 vx(d,ntotal+nvirt)=-vx(d,i)
                 enddo
                mother(ntotal+nvirt)=i
            endif
           enddo
      enddo

        if(geom.ne.0)then
          call geom_generate(itimestep,scale_k,corner,itype, hsml,slope,
     &   qp,ntotal,nvirt,nwall,x,vx,mother,ocn)
        endif
       do i=ntotal+nwall+1,ntotal+nvirt
            itype(i) = -2
            hsml(i)= 1.3*dx
            p(i)=p(mother(i))
            rho(i)=rho(mother(i))
            u(i)=u(mother(i))
           mass(i)=mass(mother(i))/ocn(mother(i))
          if(geom.eq.1)then
           dps = sqrt((x(1,i)-corner(1,5))**2+(x(2,i)-corner(2,5))**2) 
            if (x(2,i).gt.a*x(1,i)-a*xl/2+a*qp*dx-4*hsml(i).and.
     &         x(2,i).gt.-a*x(1,i)+a*xl/2+a*qp*dx-4*hsml(i))then
               mass(i)=mass(mother(i))/2
              if (dps.lt.scale_k*hsml(i))then
                 mass(i) = mass(mother(i))/3
              endif
            endif
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

      
