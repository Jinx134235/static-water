      subroutine single_step(itimestep,nstart, dt, ntotal, nvirt,hsml, 
     &           mass, x, vx,u, s,rho, p, t, tdsdt,dx,dvx, du, ds, drho, 
     &           itype, av, niac, pair_i, pair_j) 

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
c     p_record :  Record the pressure of center                    [out]
      implicit none
      include 'param.inc'

      integer itimestep, ntotal, itype(maxn), maxtimestep, scale_k,
     &        nstart
      double precision dt, hsml(maxn), mass(maxn), u(maxn), s(maxn), 
     &        rho(maxn), p(maxn),  t(maxn), tdsdt(maxn), du(maxn),
     &        ds(maxn), drho(maxn)          
      integer i, d,j,k,nvirt, niac, pair_i(max_interaction),mini,
     &        pair_j(max_interaction), ns(maxn), nwall, maxi,ntotalvirt,
     &         mother(maxn)
c      common nvirt
      double precision w(max_interaction), dwdx(3,max_interaction),  
     &       indvxdt(dim,maxn),exdvxdt(dim,maxn),ardvxdt(dim,maxn),  
     &       avdudt(maxn), ahdudt(maxn), c(maxn), eta(maxn),dis_x, 
     &       dis_y, wi(maxn) 
      double precision x(dim,maxn), vx(dim,maxn),dx(dim,maxn), 
     &       dvx(dim,maxn), av(dim,maxn), maxvel, b, minvel, vel                          
     

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
c      print *,"begining timestep"
c      print *,itimestep
c      print *,x(2,ntotal+1)
c      nvirt=0
      if (itimestep.eq.1) then 
       
         call virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,
     &       rho,u,p,itype, nwall,mother)
c        print *,"after virt_part"
c        print *,nvirt 
       open(13,file="../data/ini_virt.dat")
        do i = ntotal+1, ntotal+nvirt 
          write(13,1001) i, (x(d, i),d = 1, dim), p(i)
        enddo   
1001    format(1x, I5, 5(2x, e14.8)) 
        close(13)
      endif
       

      
c---  Interaction parameters, calculating neighboring particles
c     and optimzing smoothing length
      if (itimestep.eq.nstart+1.and.nstart.ne.0) then
         
         open(15,file="../data/xv_vp.dat")
         read(15,*) nvirt
         close(15)  
      endif
     
      ntotalvirt = ntotal + nvirt
c      print *,ntotal,nvirt
      if (nnps.eq.1) then 
        call direct_find(itimestep, ntotal,nvirt, hsml,x,niac,pair_i,
     &       pair_j,w,dwdx,ns)

            
      else if (nnps.eq.2) then
        call link_list(itimestep, ntotalvirt, hsml(1),x,niac,pair_i,
     &       pair_j,w,dwdx,ns)
c---   else if (nnps.eq.3) then 
c       call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,
c     &       pair_j,w,dwdx,ns)
      endif         
c      do k=1,50                 
c         print *,k,dwdx(dim,k),pair_i(k),pair_j(k)
c      enddo
c---  Density approximation or change rate
c---  summation_density:(4.26)
c---  con_density: calculting density through continuity equation (4.31)/(4.34)      
 
  
      call con_density(ntotal+nvirt,mass,niac,pair_i,pair_j,
     &       hsml,w,dwdx,vx,itype,x,rho,wi,drho)
   
      do i = 1,ntotal              
c          if(mod(i,39).le.10.and.i.lt.ntotal) print *,i,drho(i)`
             rho(i) = rho(i) + dt*drho(i)	 
             call p_art_water(rho(i),x(2,i),c(i),p(i))   
      enddo
c      print *,hsml(1)

      if (nor_density.and.mod(itimestep,30).eq.0) then      
        call sum_density(ntotal,hsml,mass,niac,pair_i,pair_j,w,
     &       itype,rho)         
      endif
  
        b = 1.2281e5
       do i = 1,nvirt
           p(ntotal + i) = p(mother(ntotal + i))
           rho(ntotal + i) = rho(mother(ntotal + i))
          if (itype(ntotal+i).ne.0.and.x(2,ntotal+i).lt.y_mingeom) then
c             print *,p(ntotal+i)
           p(ntotal + i) = p(mother(ntotal + i))+2*9.8*1000*
     &     (y_mingeom-x(2,ntotal+i))
           rho(ntotal + i)= 1000*(p(ntotal+i)/b+1)**(1/7)

           endif
      enddo


      if(dynamic) then
c---  Dynamic viscosity:

      if (visc) call viscosity(ntotalvirt,itype,x,rho,eta)
       
c---  Internal forces:(4.42)/(4.43)  4.58/4.59
 
      call int_force(itimestep,dt,ntotal,nvirt,hsml,mass,vx,niac,rho,
     &     eta, pair_i,pair_j,dwdx,u,itype,x,t,c,p,indvxdt,tdsdt,du) 
                  
c---  Artificial viscosity:(4.66)

      if (visc_artificial) call art_visc(ntotalvirt,hsml,
     &      mass,x,vx,niac,rho,c,pair_i,pair_j,w,dwdx,ardvxdt,avdudt)
   
      else
         if (visc) call viscosity(ntotal,itype,x,rho,eta)
       
c---  Internal forces:
 
         call int_force(itimestep,dt,ntotal,nvirt,hsml,mass,vx,niac,rho,
     &     eta, pair_i,pair_j,dwdx,u,itype,x,t,c,p,indvxdt,tdsdt,du) 
                  
c---  Artificial viscosity:

         if (visc_artificial) call art_visc(ntotal,hsml,
     &      mass,x,vx,niac,rho,c,pair_i,pair_j,w,dwdx,ardvxdt,avdudt)
      
c---  External forces:(4.93)

         if (ex_force) call ext_force(ntotalvirt,mass,x,vx,niac,
     &                  pair_i,pair_j,itype, hsml, exdvxdt)
      endif

c     Calculating the neighboring particles and undating HSML (4.80)/(4.81)
      
         if (sle.ne.0) call h_upgrade(dt,ntotal, mass, vx, rho, niac, 
     &                   pair_i, pair_j, dwdx, hsml)
c     Calculating artificial heat attached to the energy equation (4.74)
         if (heat_artificial) call art_heat(ntotalvirt,hsml,
     &         mass,x,vx,niac,rho,u, c,pair_i,pair_j,w,dwdx,ahdudt)
     
c     Calculating average velocity of each partile for avoiding penetration (4.92)

         if (average_velocity) call av_vel(ntotal,mass,niac,pair_i,
     &                           pair_j, w, vx, rho, av) 
c---  Convert velocity, force, and energy to f and dfdt  

    
     
      maxvel = 0.e0
      minvel = 1.e1

      do i=1,ntotal
c          print *,indvxdt(dim,i)
        do d=1,dim
c          dvx(1,i)=0
          dvx(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i)          
        enddo
c        if(mod(i,39).eq.1) print*,indvxdt(1,i),p(i),rho(i)
c         if(i.eq.1) print*,indvxdt(2,i),ardvxdt(2,i)
c        gravity
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
c            if(mod(i,39).le.10) print *,i,vx(1,i),vx(2,i)
c               print *,itimestep
c               stop
c             endif
c            if (x(2,i).lt.x_mingeom+scale_k*hsml(i))then
                vel = sqrt(vx(1,i)**2+vx(2,i)**2)
                if (vel.gt.maxvel)then
                maxvel = vel
                maxi = i
                endif
                if (vel.lt.minvel)then
                minvel = vel
                mini = i
                endif
c            endif
       enddo

c   update velocity of the mirror particles
       call virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,
     &       rho,u,p,itype, nwall,mother)


      if (mod(itimestep,save_step).eq.0) then
c        open(30,file="../data/trace_p.dat")
        open(40,file="../data/xv_vp.dat")
        open(50,file="../data/state_vp.dat")
        open(60,file="../data/other_vp.dat")            
        write(40,*) nvirt
        do i = ntotal + 1, ntotal + nvirt         
          write(40,1004) i, (x(d, i), d=1,dim), (vx(d, i), d = 1, dim)              
          write(50,1005) i, mass(i), rho(i), p(i), u(i)
          write(60,1006) i, itype(i), hsml(i), mother(i)                               
        enddo       
1004    format(1x, I6, 4(2x, e14.8))
1005    format(1x, I6, 4(2x, e14.8)) 
1006    format(1x, I6, 2x, I4, 2x, e14.8, 2x, I6)
        close(40)
        close(50) 
        close(60) 
      endif 
     


      if (mod(itimestep,print_step).eq.0) then      
          write(*,*) 'pressure of center',p(int(ntotal/2))
c          p_record(1) = itimestep
c          p_record(2) = p(int(ntotal/2))
c          write(*,*) '**** particle moving fastest ****', maxi         
c          write(*,101)'velocity(y)','internal(y)','total(y)'   
c          write(*,100) x(1,maxi),x(2,maxi),vx(1,maxi),vx(2,maxi)
c          write(*,101) dvx(1,maxi),dvx(2,maxi),p(maxi)
c           write(*,*) '**** particle moving slowest ****', mini         
c         write(*,102)'velocity(y)','internal(y)','total(y)'   
c          write(*,103)  x(1,mini),x(2,mini),vx(1,mini),vx(2,mini) 
      endif
      
c100   format(1x,4(2x,e12.6))     
c101   format(1x,4(2x,e12.6))          
c103   format(1x,4(2x,e12.6))      

      end
