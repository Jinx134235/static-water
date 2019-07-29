      subroutine sum_density(ntotal,hsml,mass,niac,pair_i,pair_j,w,
     &           itype,rho)

C----------------------------------------------------------------------
C   Subroutine to calculate the density with SPH summation algorithm.
c   See Equ.(4.35)

C     ntotal : Number of particles                                  [in]
C     hsml   : Smoothing Length                                     [in]
C     mass   : Particle masses                                      [in]
C     niac   : Number of interaction pairs                          [in]
C     pair_i : List of first partner of interaction pair            [in]
C     pair_j : List of second partner of interaction pair           [in]
C     w      : Kernel for all interaction pairs                     [in]
c     itype   : type of particles                                   [in]
c     x       : Coordinates of all particles                        [in]
c     rho    : Density                                             [out]
    
      implicit none
      include 'param.inc'
      
      integer ntotal, niac, pair_i(max_interaction),
     &        pair_j(max_interaction), itype(maxn), itimestep  
      double precision hsml(maxn),mass(maxn), w(max_interaction),
     &       rho(maxn) 
      integer i, j, k, d      
      double precision selfdens, hv(dim), r, wi(maxn), w_tol     

     
      
c     wi(maxn)---integration of the kernel itself
        
      do d=1,dim
        hv(d) = 0.e0
      enddo

c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.
     
c     Firstly calculate the integration of the kernel over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        wi(i)=selfdens*mass(i)/rho(i)
      enddo

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        if (mirror) then
           if (itype(i).ne.0.and.itype(j).ne.0)then
               wi(i) = wi(i) + mass(j)/rho(j)*w(k)
               wi(j) = wi(j) + mass(i)/rho(i)*w(k)
           endif
        else
         wi(i) = wi(i) + mass(j)/rho(j)*w(k)
         wi(j) = wi(j) + mass(i)/rho(i)*w(k)
        endif
      enddo

      open(1,file="../data/kernel.dat")
c     print *,int(0.5)
      do i=1,ntotal
          write (1,1001) i,wi(i)
c       print *,i,wi(i)
c         if (i.le.1600.and.abs(1-wi(i)).gt.1.e-2) then
c             print *,
c     &          ' >>> Error <<< : normalization condition unsatisfied',
c     &          i,wi(i)
c              stop
c          endif
      enddo
1001  format(1x,I6,2x,e14.8)
      close(1)
c     Secondly calculate the rho integration over the space

      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        rho(i) = selfdens*mass(i)
      enddo

c     Calculate SPH sum for rho:
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        rho(i) = rho(i) + mass(j)*w(k)
        rho(j) = rho(j) + mass(i)*w(k)
      enddo

c     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
     
      if (nor_density) then 
        do i=1, ntotal
          rho(i)=rho(i)/wi(i)
        enddo
      endif 
      
      end
      
      subroutine con_density(ntotal,mass,niac,pair_i,pair_j,
     &          hsml,w, dwdx,vx,itype,x,rho,wi,drhodt)

c----------------------------------------------------------------------
c     Subroutine to calculate the density with SPH continuity approach.
c     See Equ.(4.34)

c     ntotal : Number of particles                                  [in]
c     mass   : Particle masses                                      [in]
c     niac   : Number of interaction pairs                          [in]
c     c      : Sound speed of particles                             
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     dwdx   : derivation of Kernel for all interaction pairs       [in]
c     vx     : Velocities of all particles                          [in]
c     itype   : type of particles                                   [in]
c     x      : Coordinates of all particles                         [in]
c     rho    : Density                                              [in]
c     drhodt : Density change rate of each particle                [out]
C     w      : Kernel for all interaction pairs                     [in]
c     norrho : normalized density of all particles                 [out]   

      implicit none
      include 'param.inc'
     

      integer ntotal,niac,pair_i(max_interaction),
     &        pair_j(max_interaction), itype(maxn), itimestep    
      double precision mass(maxn), dwdx(3,max_interaction),
     &       vx(dim,maxn), x(dim,maxn), rho(maxn), drhodt(maxn),
     &       hsml(maxn), w(max_interaction)
      integer i,j,k,d    
      double precision   vcc, dvx(dim),delta, r,c, dx(dim),hv(dim),
     &       selfdens,wi(maxn), psi(dim), xcc, b, rhoh
c      double precision, intent,output:: wi(maxn) 
c      real wi(maxn)

      do d=1,dim
        hv(d) = 0.e0
      enddo

c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.
      do i =1,ntotal
         wi(i) = 0.
      enddo

c     Firstly calculate the integration of the kernel over the space
c     checking the normalization condition 
      do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        wi(i)=selfdens*mass(i)/rho(i)
        
c        print *,i,selfdens,rho(i),wi(i)    
      enddo

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
c        if (pair_i(k).eq.1) print *, pair_j(k)
         wi(i) = wi(i) + mass(j)/rho(j)*w(k)
         wi(j) = wi(j) + mass(i)/rho(i)*w(k)
       
      enddo
      
      open(10, file="../data/kernel.dat")
    
      do i = 1,ntotal
         write (10,1001) i,wi(i)
       enddo
 1001     format(2x, I4, 2x, e14.8)
      close(10)

c  check normalization condition
c      do i=1,ntotal
c       if (i.eq.moni_particle) print *,'wi(1600):',wi(moni_particle)
c         if (i.le.120) print *,i,wi(i)
c         if (abs(1-wi(i)).gt.1.e-2) then
c             print *,
c     &          ' >>> Error <<< : normalization condition unsatisfied'
c              stop
c         endif
c      enddo

c     Secondly calculate the rho integration over the space

   
      do i = 1, ntotal
        drhodt(i) = 0.
c     density correction(Shepard filter)
c        if (nor_density.and.mod(itimestep,30).eq.0) then 
c          rho(i) = rho(i)/wi(i)
c       endif
      enddo
     
      delta  = 0.1
      c = 29.32
      b = c**2*1000/7
      do k=1,niac      
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,i) - vx(d,j) 
          dx(d) = x(d,i) - x(d,j)
          r = r+dx(d)**2
        enddo
        vcc = dvx(1)*dwdx(1,k) 
         do d=2,dim
          vcc = vcc + dvx(d)*dwdx(d,k)
        enddo

        drhodt(i) = drhodt(i) + mass(j)*vcc
        drhodt(j) = drhodt(j) + mass(i)*vcc
c  add filter to the continuity equation(Molteni,2009)
        if (filt.eq.1) then
         do d=1,dim
            psi(d) = 2*(rho(j)-rho(i))*dx(d)/sqrt(r)
          enddo
         xcc = psi(1)*dwdx(1,k)
         do d=2,dim
           xcc = xcc+psi(d)*dwdx(d,k)
         enddo 
          drhodt(i) = drhodt(i) + delta*hsml(i)*c*xcc*mass(j)/rho(j)
c           if (i.eq.1) print *,'after filter',drhodt(i) 
          drhodt(j) = drhodt(j) + delta*hsml(j)*c*xcc*mass(i)/rho(i)
        elseif (filt.eq.2) then
           rhoh = 1000*(1000*9.8*dx(dim)/b+1)**(1/7)
c           if(k.le.10)print *,rhoh  
           do d=1,dim
            psi(d) = 2*(rho(j)-rho(i)-rhoh)*dx(d)/sqrt(r)
          enddo
c          if(k.eq.1) print*,psi
           xcc = psi(1)*dwdx(1,k)
           do d=2,dim
           xcc = xcc+psi(d)*dwdx(d,k)
           enddo
c           if (i.eq.1) print *,'before filter',drhodt(i) 
          drhodt(i) = drhodt(i) + delta*hsml(i)*c*xcc*mass(j)/rho(j)
c           if (i.eq.1) print *,'after filter',drhodt(i) 
          drhodt(j) = drhodt(j) + delta*hsml(j)*c*xcc*mass(i)/rho(i)
c        elseif (filt.eq.3) then
           
         endif

       enddo    
    
      end
