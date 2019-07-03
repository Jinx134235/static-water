      subroutine single_step(itimestep, dt, ntotal,nvirt, hsml, mass, x,
     &           vx,u, s,rho, p, t, tdsdt,dx,dvx, du, ds, drho,itype,av,  
     &           niac, pair_i, pair_j) 

c----------------------------------------------------------------------
c   Subroutine to determine the right hand side of a differential 
c   equation in a single step for performing time integration 

c   In this routine and its subroutines the SPH algorithms are performed.
c     itimestep: Current timestep number                            [in]
c     dt       : Timestep                                           [in]
c     ntotal   :  Number of particles                               [in]
c     hsml     :  Smoothing Length                                  [in]
c     mass     :  Particle masses                                   [in]
c     x        :  Particle position                                 [in]
c     vx       :  Particle velocity                                 [in]
c     u        :  Particle internal energy                          [in]
c     s        :  Particle entropy (not used here)                  [in]
c     rho      :  Density                                       [in/out]
c     p        :  Pressure                                         [out]
c     t        :  Temperature                                   [in/out]
c     tdsdt    :  Production of viscous entropy t*ds/dt            [out]
c     dx       :  dx = vx = dx/dt                                  [out]
c     dvx      :  dvx = dvx/dt, force per unit mass                [out]
c     du       :  du  = du/dt                                      [out]
c     ds       :  ds  = ds/dt                                      [out]     
c     drho     :  drho =  drho/dt                                  [out]
c     itype    :  Type of particle                                 [in]
c     av       :  Monaghan average velocity                        [out]

      implicit none
      include 'param.inc'

      integer itimestep, ntotal, itype(maxn)
      double precision dt, hsml(maxn), mass(maxn),  u(maxn), s(maxn), 
     &        rho(maxn), p(maxn),  t(maxn), tdsdt(maxn),  du(maxn),
     &        ds(maxn), drho(maxn)          
      integer i, d,j,k, nvirt, niac, pair_i(max_interaction),
     &        pair_j(max_interaction), ns(maxn), nwall, maxi
      double precision w(max_interaction), dwdx(3,max_interaction),  
     &       indvxdt(dim,maxn),exdvxdt(dim,maxn),ardvxdt(dim,maxn),  
     &       avdudt(maxn), ahdudt(maxn), c(maxn), eta(maxn),dis_x, 
     &       dis_y 
      double precision x(dim,maxn), vx(dim,maxn),dx(dim,maxn), 
     &       dvx(dim,maxn), av(dim,maxn), maxvel, beta                          

      do  i=1,ntotal
        avdudt(i) = 0.
        ahdudt(i) = 0.
        do  d=1,dim
          indvxdt(d,i) = 0.
          ardvxdt(d,i) = 0.
          exdvxdt(d,i) = 0.
        enddo
      enddo  
 
c---  Positions of virtual (boundary) particles: 

     
      if (virtual_part) then 
       
         call virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,
     &       rho,u,p,itype, nwall)
       
      endif 
     
c---  Interaction parameters, calculating neighboring particles
c     and optimzing smoothing length

      if (nnps.eq.1) then 
        call direct_find(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
      else if (nnps.eq.2) then
        call link_list(itimestep, ntotal+nvirt,hsml(1),x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
c---   else if (nnps.eq.3) then 
c       call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
c     &       pair_j,w,dwdx,ns)
      endif         
                        
c---  Density approximation or change rate
c---  summation_density:(4.26)
c---  con_density: calculting density through continuity equation (4.31)/(4.34)      
      if (summation_density) then      
        call sum_density(ntotal+nvirt,hsml,mass,niac,pair_i,pair_j,w,
     &       itype,rho)          
        call p_art_water(rho(i), x(2,i),p(i),c(i))     
c      else   
        
      endif

      if(dynamic) then
c---  Dynamic viscosity:

      if (visc) call viscosity(ntotal+nvirt,itype,x,rho,eta)
       
c---  Internal forces:(4.42)/(4.43)  4.58/4.59
 
      call int_force(itimestep,dt,ntotal+nvirt,hsml,mass,vx,niac,rho,
     &     eta, pair_i,pair_j,dwdx,u,itype,x,t,c,p,indvxdt,tdsdt,du) 
                  
c---  Artificial viscosity:(4.66)

      if (visc_artificial) call art_visc(ntotal+nvirt,hsml,
     &      mass,x,vx,niac,rho,c,pair_i,pair_j,w,dwdx,ardvxdt,avdudt)
   
      else
         if (visc) call viscosity(ntotal,itype,x,rho,eta)
       
c---  Internal forces:
 
         call int_force(itimestep,dt,ntotal+nwall,hsml,mass,vx,niac,rho,
     &     eta, pair_i,pair_j,dwdx,u,itype,x,t,c,p,indvxdt,tdsdt,du) 
                  
c---  Artificial viscosity:

         if (visc_artificial) call art_visc(ntotal+nwall,hsml,
     &      mass,x,vx,niac,rho,c,pair_i,pair_j,w,dwdx,ardvxdt,avdudt)
      
c---  External forces:(4.93)

         if (ex_force) call ext_force(ntotal+nvirt,mass,x,vx,niac,
     &                  pair_i,pair_j,itype, hsml, exdvxdt)
      endif

c     Calculating the neighboring particles and undating HSML (4.80)/(4.81)
      
         if (sle.ne.0) call h_upgrade(dt,ntotal, mass, vx, rho, niac, 
     &                   pair_i, pair_j, dwdx, hsml)
c     Calculating artificial heat attached to the energy equation (4.74)
         if (heat_artificial) call art_heat(ntotal+nvirt,hsml,
     &         mass,x,vx,niac,rho,u, c,pair_i,pair_j,w,dwdx,ahdudt)
     
c     Calculating average velocity of each partile for avoiding penetration (4.92)

         if (average_velocity) call av_vel(ntotal,mass,niac,pair_i,
     &                           pair_j, w, vx, rho, av) 
c---  Convert velocity, force, and energy to f and dfdt  

      maxvel = 0.e0
     
      do i=1,ntotal+nwall

        do d=1,dim
c          dvx(1,i)=0
          dvx(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i)          
        enddo
        if (self_gravity) then
          dvx(dim, i) = dvx(dim,i)-9.8
         endif

          du(i) = du(i) + avdudt(i) + ahdudt(i)
           u(i) = u(i) + dt*du(i)
            if(u(i).lt.0)  u(i) = 0.         
            do d = 1, dim                   
              vx(d, i) = vx(d, i) + dt * dvx(d, i) + av(d, i)
              x(d, i) = x(d, i) + dt * vx(d, i)       
            enddo
            if (abs(vx(2,i)).gt.maxvel)then
                maxvel = vx(2,i)
                maxi = i
            endif
c   update density & pressure
       enddo

c       print *,p(1)

      call con_density(ntotal+nvirt,mass,niac,pair_i,pair_j,hsml,w,
     &       dwdx,vx,itype,x,rho,drho)
      
      do i = 1,ntotal                   
           rho(i) = rho(i) + dt*drho(i)
	     call p_art_water(rho(i),x(2,i),c(i),p(i))
c           p(i)=p(i)+9.8*1000*(1.e-3-x(2,i))
c          if (p(i).lt.0) p(i)=0
      enddo

c   correct pressure for mirror particle
      beta = 315.e0
      do i = 1,nvirt
          if (itype(ntotal+i).ne.0.and.x(2,ntotal+i).lt.0)then
           p(ntotal + i) = p(ntotal+i)-2*9.8*1000*x(2,ntotal+i)
            rho(ntotal + i)= 1000*(p(ntotal+i)/beta+1)**(1/7)
           endif
      enddo

      if (mod(itimestep,print_step).eq.0) then      
          write(*,*)
          write(*,*) '**** particle moving fastest ****', 
     &        		maxi         
          write(*,101)'velocity(y)','internal(y)','total(y)'   
          write(*,100)vx(2,maxi),indvxdt(2,maxi),dvx(2,maxi)          
      endif
101   format(1x,3(2x,a12))      
100   format(1x,3(2x,e12.6))      

      end
