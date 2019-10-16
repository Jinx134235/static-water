      subroutine input(x, vx, mass, rho, p, u, itype, hsml,
     &    ntotal)
    
c----------------------------------------------------------------------
c     Subroutine for loading or generating initial particle information

c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c     hsml-- smoothing lengths of particles                        [out]
c     ntotal-- total particle number                               [out]

      implicit none     
      include 'param.inc'

      integer itype(maxn), ntotal, np
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn), 
     &                 p(maxn), u(maxn), hsml(maxn), rho(maxn), dx
      integer i, d, im       

c     load initial particle information from external disk file

      if(config_input) then    
                        
        open(10,file="../data/f_xv.dat")
        open(20,file="../data/f_state.dat")
        open(30,file="../data/f_other.dat")        
      
        write(*,*)'  **************************************************'
        write(*,*)'      Loading initial particle configuration...   '       
	  read (10,*) ntotal 
        write(*,*)'      Total number of particles   ', ntotal    	
        write(*,*)'  **************************************************'	
        do i = 1, ntotal         
          read(10,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)                     
          read(20,*)im, mass(i), rho(i), p(i), u(i)        
          read(30,*)im, itype(i), hsml(i)                                        
        enddo          
          
      else 
          
        open(11,file="../data/ini_xv.dat")
       open(22,file="../data/ini_state.dat")
       open(33,file="../data/ini_other.dat") 
       
      
      if(dambreak) call dam_break(x, vx, mass, rho, p, u, 
     &                      itype, hsml, ntotal)
      if(static.and..not.geometry) call static_water(x, vx, mass, rho,
     &                     p,u, itype, hsml, ntotal)
      if(geometry) call two_phase(x, vx, mass, rho, p, u,
     &                      itype, hsml, ntotal)
      
      do i = 1, ntotal 
          write(11,1001) i, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim) 
          write(22,1002) i, mass(i), rho(i), p(i), u(i)         
          write(33,1003) i, itype(i), hsml(i)    
        enddo   
1001    format(1x, I5, 4(2x, e14.8)) 
1002    format(1x, I5, 4(2x, e14.8)) 
1003    format(1x, I5, 2x, I2, 2x, e14.8) 
        write(*,*)'  **************************************************'
        write(*,*)'      Initial particle configuration generated   '       
        write(*,*)'      Total number of particles(fluid)   ', ntotal    	
        write(*,*)'  **************************************************' 

      endif

      close(11)
      close(22) 
      close(33) 

      end              
       
      
      subroutine dam_break(x, vx, mass, rho, p, u,itype, hsml, ntotal)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     2-d  dambreak benchmark
c     x-- coordinates of particles                                 [out]
c     vx-- velocities of particles                                 [out]
c     mass-- mass of particles                                     [out]
c     rho-- dnesities of particles                                 [out]
c     p-- pressure  of particles                                   [out]
c     u-- internal energy of particles                             [out]
c     itype-- types of particles                                   [out]
c          =2   water
c     h-- smoothing lengths of particles                           [out]
c     ntotal-- total particle number                               [out]
c     dx-- initial interval among particles                        [out]
c     np-- total particle number in one column                     [out]
      implicit none     
      include 'param.inc'
      
      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn), hsml(maxn)
      integer i, j, d, m, n, mp, np,mnp, nnp,k
      double precision xl, yl, dx, dy, gate_x


      m = 1872
      n = 104
  
      xl = x_maxgeom-x_mingeom
      yl = y_maxgeom-y_mingeom
      dx = xl/m
      dy = yl/n

      mp = int(damlength/dx)
      np = int(damheight/dy)
      nnp = int(bedheight/dy)
      mnp = (x_maxgeom-damlength)/dx
      ntotal = mp*np+(mp-1)*(np-1)+mnp*nnp

c      print *,mp,np
c   particles of dam part 
      do i = 1, mp
	  do j = 1, np
	      k = j + (i-1)*np
	      x(1, k) = (i-1)*dx+dx/2+x_mingeom
	      x(2, k) = (j-1)*dy+dy/2+y_mingeom
         enddo
      enddo
c  staggered tencil
      do i = 1,mp-1
         do j = 1,np-1
           k = mp*np + j +(i-1)*(np-1)
           x(1,k) = i*dx+x_mingeom
           x(2,k) = j*dy+y_mingeom
         enddo
      enddo 
    
c  particles of  bed part
      do i = 1,mnp
        do j = 1,nnp
         k = mp*np+(mp-1)*(np-1)+j+(i-1)*nnp
         x(1,k) = damlength+dx/2+(i-1)*dx
         x(2,k) = y_mingeom+dy/2+(j-1)*dy
         enddo
      enddo

      do i = 1, ntotal
       	vx(1, i) = 0.
	  vx(2, i) = 0.     
c--- original density,pressure & mass of the particles    
c--- zero pressure
c        p(i) = 0
       if (x(1,k).lt.gate_x) then
          p(i)=9.8*1000*(np*dy-x(2,i))
       else
         p(i)=9.8*1000*(nnp*dy-x(2,i))
       endif   
c        rho(i)= 1000*(p(i)/20+1)**(1/7)
        rho(i) = 1000   
        mass(i) = dx*dy*rho(i)  
        u(i)=357.1
        itype(i) = 2
        hsml(i) = 1.3*dx
      enddo  

      end

      subroutine static_water(x, vx, mass, rho, p, u,
     &           itype, hsml, ntotal)
c---------------------------------------------------------------
c   input for static water benchmark
c   parameter list is the same as above, except for the digitals

              
      implicit none
      include 'param.inc'

      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn), hsml(maxn)
      integer i, j, d, np,k
      double precision xl, yl, dx, dy

c      m = 104
c      n = 104
      np = mmp

      xl = x_maxgeom-x_mingeom
      yl = y_maxgeom-y_mingeom

      dx = xl/mmp
      dy = yl/np
c      mp = m
c      np = n
      ntotal = mmp*np

      do i = 1, mmp
         do j= 1, np
           k = j+(i-1)*np
           x(1,k) = x_mingeom+(i-1)*dx+dx/2
           x(2,k) = y_mingeom+(j-1)*dy+dy/2
         enddo
      enddo

      do i = 1, ntotal
        vx(1, i) = 0.
        vx(2, i) = 0.
c--- original density,pressure & mass of the particles    
c--- either zero pressure or hydrostatic works
c        p(i) = 0
        p(i) = 9.8*1000*(yl-x(2,i))
        rho(i) = 1000
        mass(i) = dx*dy*rho(i)
        u(i)=357.1
        itype(i) = 2
        hsml(i) = 1.3*dx
      enddo

      end


      subroutine two_phase(x, vx, mass, rho, p, u, itype, hsml, ntotal)
c---------------------------------------------------------------
c   input for static water benchmark with two phases, and an obstacle
c   on the wall, which can be in any geomatrical shape.In this case,we
c   choose it as isoceles right triangle
c   parameter list is the same as above, except for the digitals

c   this subroutine contains two set of testing geometry, which are like:
c   |                         I
c   |                         |
c   L__________               |
c             /               |             /\
c            /                |            /  \
c           /                 |           /    \
c          /                  |     _____/      \_____
c         /                   |
c        /                    |
c       /                     I
c   on the left wall                  in the bottom  

c   mp:rows of particle above the obstacle
c   np:rows of particle below the obstacle
c   qp:colums of particle right below the obstacle
c   hp:rows of praticle in high-density domain
c   theta:angle of the obstacle, default:pi/3 
c   a: slope of the line(obstacle)

      implicit none
      include 'param.inc'

      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn), hsml(maxn)
      integer i, j, d, m, n, mp, np,qp,hp,k
      double precision xl, yl, dx, dy, theta, a, y1, y2, y3, dis,x1


      m = 35
      n = 60
      theta = pi/3
      a = tan((pi-theta)/2)

      xl = x_maxgeom-x_mingeom
      dx = xl/mmp
c      dy = (y_maxgeom-y_mingeom)/n
      
      mp = 40
      np = int(n/a)
      qp = int((n-mp)/a)
c      hp = 20
      
c      ntotal = m*(mp+np)-qp*(qp+1)/2
c    geometry 1      
c      do i = 1,qp
c         do j = 1,np-qp+i-1
c           k = j+(np-qp)*(i-1)+(i-1)*(i-2)/2
c          x(1,k) = (i-1)*dx+dx/2+x_mingeom
c          x(2,k) = (j-1)*dy+dy/2+y_mingeom
c          enddo          
c       enddo

c      do i = 1,m-qp
c         do j = 1,np
c           k = qp*((np-qp)*2+qp-1)/2+(i-1)*np+j
c           x(1,k) = (qp+i)*dx-dx/2+x_mingeom 
c           x(2,k) = (j-1)*dy+dy/2+y_mingeom
c         enddo
c       enddo

c       do i = 1,m
c          do j = 1,mp
c            k = m*np-qp*(qp+1)/2+(i-1)*mp+j
c            x(1,k) = (i-1)*dx+dx/2+x_mingeom
c            x(2,k) = (np+j)*dy-dy/2+y_mingeom
c           enddo
c       enddo

c  geometry 2
      if(indis.eq.0) then
        do i = 1,m
          do j = 1,n
c         y1 = a*(i*dx-dx/2)+(n-mp)*dx-a*xl/2
c          y2 = a*(dx/2-i*dx)+(n-mp)*dx+a*xl/2
c          if(j*dx-dx/2.gt.y1.or.j*dx-dx/2.gt.y2)then
            ntotal = ntotal + 1
            x(1,ntotal) = x_mingeom + i*dx-dx/2
            x(2,ntotal) = y_mingeom + j*dx-dx/2
c          endif
         enddo
        enddo

c   adjust the y-coordinate a little for those near to wedge
c   to make the paricles distribute more evenly while not adding or
c   removing any particle, this may be useful in minor grid     
c       do i = 1,ntotal
c          if (x(1,i).lt.xl/2)then
c             dis = abs((a*x(1,i)-x(2,i)+(n-mp)*dx-a*xl/2))/sqrt(a**2+1)
c             if(dis.lt.dx/2.and.x(2,i)+dx-dis*2.le.x(2,i+1)-dx/2)then
c               x(2,i) = x(2,i)+dx-dis*2
c             endif
c          else
c            dis = abs((-a*x(1,i)-x(2,i)+(n-mp)*dx+a*xl/2))/sqrt(a**2+1)
c             if(dis.lt.dx/2.and.x(2,i)+dx-dis*2.le.x(2,i+1)-dx/2) then
c               x(2,i) = x(2,i)+dx-dis*2
c             endif
c          endif
c       enddo
      else if(indis.eq.1) then     
c   distribution 2:
       do i = 1,m
        do j = 1,mp 
c          y1 = a*(i*dx-dx/2)+(n-mp)*dx-a*xl/2
c          y2 = a*(dx/2-i*dx)+(n-mp)*dx+a*xl/2
c          y3 = sqrt(((a+.5)*dx)**2-(i*dx-dx/2-xl/2)**2)+(n-mp+1)*dx 
c         if (j*dx-dx/2.gt.y1+5.5*dx.or.j*dx-dx/2.gt.y2+5.5*dx.or.
c          if( j*dx-(a-1)*dx.gt.(n-mp)*dx) then
          ntotal = ntotal + 1
          x(1,ntotal) = x_mingeom + i*dx-dx/2
          x(2,ntotal) = y_mingeom + (j-1)*dx+(n-mp)*dx+dx/2
c         endif 
        enddo
       enddo
c   distribute particles along wall via transformation of coordinates 
c       do j = 1,2*qp+2
c        do i =1,m/2-qp+1
c            x1 = cos(theta)*((j-1)*dx+dx/4-(i-1)*dx/a)-
c     &    sin(theta)*((i-1)*dx+a*dx/4)+xl/2-(n-mp)*dx/a
c            y1 = cos(theta)*((i-1)*dx+a*dx/4)+sin(theta)*
c     &    ((j-1)*dx+dx/4-(i-1)*dx/a)
c            if (x1.gt.x_mingeom) then  
c            ntotal = ntotal+ 2
c            x(1,ntotal-1) = x1
c             x(1,ntotal) = xl-x1
c            x(2,ntotal-1) = y1
c            x(2,ntotal) = x(2,ntotal-1)
c           endif
c        enddo
c       enddo
        do i = 1,n-mp
          do j = 1,m/2
            y1 = (i-1)*dx+dx/2
            x1 = (y1-(n-mp)*dx+a*xl/2)/a-(j-1)*dx-dx/2
            if(x1.gt.x_mingeom) then
               ntotal = ntotal +2 
               x(1,ntotal-1) = x1
               x(2,ntotal-1) = y1
               x(1,ntotal) = xl-x1
               x(2,ntotal) = y1
             endif
          enddo
        enddo

c    some particle around the corner(coordinates change)            
    
c          hp = int(pi*(a+1)/2)
c          print *,3./4
c          do j =1,hp-1
c           ntotal = ntotal +1
c           print *,real(j)/hp
c   !!! caution: j is an int, when it is devided by a bigger int,
c   convert it into real type first           
c            x(1,ntotal) = xl/2+(a+1)*dx/2*cos((real(j)/hp)*pi)
c            x(2,ntotal) = (n-mp)*dx+(a+1)*dx/2*sin((real(j)/hp)*pi)
c          enddo
c          hp = int(pi*(a+1./2)) 
c          print *,hp
c          do j =1,hp-1
c           ntotal = ntotal +1
c            x(1,ntotal) = xl/2+(2*a+1)*dx/2*cos((real(j)/hp)*pi)
c            x(2,ntotal) = (n-mp)*dx+(2*a+1)*dx/2*sin((real(j)/hp)*pi)
c          enddo

      else if(indis.eq.2) then
       do i = 1,m
        do j = 1,np
         y1 = a*(i*dx-dx/2)+(n-mp)*dx-a*xl/2
          y2 = a*(dx/2-i*dx)+(n-mp)*dx+a*xl/2
c   distribution 3:another set of grid
c    . . . .
c     . . .
c    . . . .
          if (a*(j-.5)*dx.gt.y1.or.a*(j-.5)*dx.gt.y2) then
            ntotal = ntotal + 1
             x(1,ntotal) = x_mingeom+i*dx-dx/2
             x(2,ntotal) = y_mingeom+a*(j-1)*dx+a*dx/2
          endif
         enddo
       enddo

       do i = 1,m+1
        do j = 1,np+1
          y1 = a*(i-1)*dx+(n-mp)*dx-a*xl/2
           y2 = -a*(i-1)*dx+(n-mp)*dx+a*xl/2
          if (a*(j-1)*dx.gt.y1.or.a*(j-1)*dx.gt.y2) then
             ntotal = ntotal + 1
             x(1,ntotal) = x_mingeom+(i-1)*dx
             x(2,ntotal) = y_mingeom+a*(j-1)*dx
          endif
         enddo
        enddo
      endif
      
      do i = 1,ntotal   
        vx(1, i) = 0.
        vx(2, i) = 0.
c--- original density,pressure & mass of the particles    
c--- zero pressure
        p(i) = 0
c--- hydrostatic pressure        
c        p(i) = 9.8*1000*(yl-x(2,i))
        rho(i) = 1000
        itype(i)  = 2
c     low-density phase        
c       if(x(1,i).gt.dx*m/2.or.x(2,i).gt.hp*dy)then
c           rho(i)= 250
c           itype(i) = 3
c        endif
        mass(i) = dx*dx*rho(i)
        u(i) = 357.1
        hsml(i) = 1.3*dx
      enddo

      end



