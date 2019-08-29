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
                        
        open(1,file="../data/f_xv.dat")
        open(2,file="../data/f_state.dat")
        open(3,file="../data/f_other.dat")        
      
        write(*,*)'  **************************************************'
        write(*,*)'      Loading initial particle configuration...   '       
	  read (1,*) ntotal 
        write(*,*)'      Total number of particles   ', ntotal    	
        write(*,*)'  **************************************************'	
        do i = 1, ntotal         
          read(1,*)im, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)                     
          read(2,*)im, mass(i), rho(i), p(i), u(i)        
          read(3,*)im, itype(i), hsml(i)                                        
        enddo          
          
      else 
          
        open(1,file="../data/ini_xv.dat")
       open(2,file="../data/ini_state.dat")
       open(3,file="../data/ini_other.dat") 
       
      
      if(dambreak) call dam_break(x, vx, mass, rho, p, u, 
     &                      itype, hsml, ntotal)
      if(static) call static_water(x, vx, mass, rho, p, u,
     &                      itype, hsml, ntotal)

        do i = 1, ntotal 
          write(1,1001) i, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim) 
          write(2,1002) i, mass(i), rho(i), p(i), u(i)         
          write(3,1003) i, itype(i), hsml(i)    
        enddo   
1001    format(1x, I5, 4(2x, e14.8)) 
1002    format(1x, I5, 4(2x, e14.8)) 
1003    format(1x, I5, 2x, I2, 2x, e14.8) 
        write(*,*)'  **************************************************'
        write(*,*)'      Initial particle configuration generated   '       
        write(*,*)'      Total number of particles(fluid)   ', ntotal    	
        write(*,*)'  **************************************************' 

      endif

      close(1)
      close(2) 
      close(3) 

      end              
       
      
      subroutine dam_break(x, vx, mass, rho, p, u, 
     &                        itype, hsml, ntotal)

c----------------------------------------------------------------------     
c     This subroutine is used to generate initial data for the 
c     2-d static water benchmark as well as dambreak with Re = 1
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
      integer i, j, d, m, n, mp, np,mmp, nnp,k
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
      mmp = (x_maxgeom-damlength)/dx
      ntotal = mp*np+(mp-1)*(np-1)+mmp*nnp

      print *,mp,np
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
      do i = 1,mmp
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
     &                        itype, hsml, ntotal)
c---------------------------------------------------------------
c   input for static water benchmark
c   parameter list is the same as above, except for the digitals

              
      implicit none
      include 'param.inc'

      integer itype(maxn), ntotal
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn),
     &     rho(maxn), p(maxn), u(maxn), hsml(maxn)
      integer i, j, d, m, n, mp, np,k
      double precision xl, yl, dx, dy

      m = 39
      n = 39

      xl = x_maxgeom-x_mingeom
      yl = y_maxgeom-y_mingeom

      dx = xl/m
      dy = yl/n
      mp = m
      np = n
      ntotal = mp*np

      do i = 1 , mp
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
c--- zero pressure
        p(i) = 0
c        p(i) = 9.8*1000*(yl-x(2,i))
        rho(i) = 1000
        mass(i) = dx*dy*rho(i)
        u(i)=357.1
        itype(i) = 2
        hsml(i) = 1.3*dx
      enddo

      end
