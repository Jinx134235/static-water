      subroutine ext_force(ntotal,nvirt,nwall,mass,x,vx,
     &       itype,hsml,maxvel,dvxdt)

c--------------------------------------------------------------------------
c     Subroutine to calculate the external forces, e.g. gravitational forces.      
c     The forces from the interactions with boundary virtual particles 
c     are also calculated here as external forces.

c     here as the external force. 
c     ntotal  : Number of particles                                 [in]
c     mass    : Particle masses                                     [in]
c     x       : Coordinates of all particles                        [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     itype   : type of particles                                   [in]
c     hsml   : Smoothing Length                                     [in]
c     maxvel : maximum of velocity                                  [in]
c     dvxdt   : Acceleration with respect to x, y and z            [out] 

      implicit none
      include 'param.inc'
      
      integer ntotal, itype(maxn), niac,nvirt,nwall,
     &        pair_i(max_interaction), pair_j(max_interaction)
      double precision mass(maxn), x(dim,maxn), hsml(maxn),          
     &       dvxdt(dim,maxn),vx(dim,maxn),maxvel
      integer i, j, k, d
      double precision dx(dim), rr, f, rr0, dd, p1, p2     
           
      do i = 1, ntotal
        do d = 1, dim
          dvxdt(d, i) = 0.
	enddo
      enddo
        
c     Boundary particle force and penalty anti-penetration force. 
      rr0 = 5.e-3
      dd = c0**2
c      print *,dd
      p1 = 12
      p2 = 4
      
c      do  k=1,niac
c        i = pair_i(k)
c        j = pair_j(k)
c     only for the wall particles
      do i = 1,ntotal
        do j = ntotal+nvirt+1,ntotal+nvirt+nwall     
c        if(itype(ntotal+j).eq.0) then  
          rr = 0.              
          do d=1,dim
            dx(d) =  x(d,i) -  x(d,j)
c           rd(d) = vx(d,i) - vx(d,j)
            rr = rr + dx(d)*dx(d)
          enddo  
          rr = sqrt(rr)
          if(rr.lt.rr0) then
            f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
            do d = 1, dim
              dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f
            enddo
          endif
c        endif 
        enddo       
      enddo   
       
      end         
