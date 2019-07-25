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
     
      if (nor_density.and.mod(itimestep,30).eq.0) then 
        do i=1, ntotal
          rho(i)=rho(i)/wi(i)
        enddo
      endif 
      
      end
      
      subroutine con_density(ntotal,mass,niac,pair_i,pair_j,hsml,w,
     &           dwdx,vx,itype,x,rho,drhodt)

c----------------------------------------------------------------------
c     Subroutine to calculate the density with SPH continuity approach.
c     See Equ.(4.34)

c     ntotal : Number of particles                                  [in]
c     mass   : Particle masses                                      [in]
c     niac   : Number of interaction pairs                          [in]
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
     &        pair_j(max_interaction), itype(maxn)    
      double precision mass(maxn), dwdx(3,max_interaction),
     &       vx(dim,maxn), x(dim,maxn), rho(maxn), drhodt(maxn),
     &       hsml(maxn), w(max_interaction)
      integer i,j,k,d    
      double precision   vcc, dvx(dim),delta, r,c, dx(dim),hv(dim),
     &       selfdens, wi(maxn),psi(dim), xcc
c      real wi(maxn)

      do d=1,dim
        hv(d) = 0.e0
      enddo

c     Self density of each particle: Wii (Kernel for distance 0)
c     and take contribution of particle itself:

      r=0.
      
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
      
      open(1, file="../data/kernel.dat")
    
      do i = 1,ntotal
         write (1,1001) i,wi(i)
       enddo
 1001     format(2x, I4, 2x, e14.8)
      close(1)

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
      enddo
     
      delta  = 0.
      c = 29.32
      do k=1,niac      
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,i) - vx(d,j) 
          dx(d) = x(d,i) - x(d,j)
          r = r+dx(d)**2
        enddo
        do d=1,dim
          psi(d) = 2*(rho(j)-rho(i))*dx(d)/sqrt(r)
        enddo  
        vcc = dvx(1)*dwdx(1,k) 
        xcc = psi(1)*dwdx(1,k)
        do d=2,dim
          vcc = vcc + dvx(d)*dwdx(d,k)
          xcc = xcc + psi(d)*dwdx(d,k)
        enddo    
        drhodt(i) = drhodt(i) + mass(j)*vcc
        drhodt(j) = drhodt(j) + mass(i)*vcc
c  add filter to the continuity equation(Molteni,2009)
        if (filt_density) then
c              if (i.eq.123)print *,"before filter:",drhodt(i)      
          drhodt(i) = drhodt(i) + delta*hsml(i)*c*xcc*mass(j)/rho(j)
c              if (i.eq.123)print *,"after filtered:",drhodt(i)
          drhodt(j) = drhodt(j) + delta*hsml(j)*c*xcc*mass(i)/rho(i)
        endif
       enddo    
    
c       print *,"after filtered" ,drhodt(1)
      end
