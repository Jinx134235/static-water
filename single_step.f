      subroutine single_step(itimestep,nstart, dt, ntotal, nvirt,nwall, 
     &           hsml,mass, x, vx,u, s,rho, p, t, tdsdt, du,ds,c,
     &           itype, av, niac, pair_i, pair_j, sumvel) 

c----------------------------------------------------------------------
c   Subroutine to determine the right hand side of a differential 
c   equation in a single step for performing time integration 

c   In this routine and its subroutines the SPH algorithms are performed.
c     itimestep: Current timestep number                            [in]
c     dt       : Timestep                                           [in]
c     ntotal   :  Number of particles                               [in]
c     hsml     :  Smoothing Length                                  [in]
c     mass     :  Particle masses                                   [in]
c     x        :  Particle position                             [in/out]
c     vx       :  Particle velocity                             [in/out]
c     u        :  Particle internal energy                          [in]
c     s        :  Particle entropy (not used here)                  [in]
c     rho      :  Density                                       [in/out]
c     p        :  Pressure                                      [in/out]
c     t        :  Temperature                                   [in/out]
c     tdsdt    :  Production of viscous entropy t*ds/dt            [out]
c     du       :  du  = du/dt                                      [out]
c     ds       :  ds  = ds/dt                                      [out]     
c     drho     :  drho =  drho/dt                                  [out]
c     itype    :  Type of particle                                 [in]
c     av       :  Monaghan average velocity                        [out]
c     sumvel   :  Summation of velocity                            [out]
      implicit none
      include 'param.inc'

      integer itimestep, ntotal, itype(maxn), maxtimestep, scale_k,
     &        nstart, np,nnp
      double precision dt, hsml(maxn), mass(maxn), u(maxn), s(maxn), 
     &        rho(maxn), p(maxn),  t(maxn), tdsdt(maxn), du(maxn),
     &        ds(maxn), drho(maxn), pp(maxn),sumw(maxn)
      integer i, d,j,k,nvirt, niac, pair_i(max_interaction),mini,
     &        pair_j(max_interaction), ns(maxn), nwall, maxi,ntotalvirt,
     &         mother(maxn)
c      common nvirt
      double precision w(max_interaction), dwdx(3,max_interaction),  
     &       indvxdt(dim,maxn),exdvxdt(dim,maxn),ardvxdt(dim,maxn),  
     &       avdudt(maxn), ahdudt(maxn), c(maxn), eta(maxn),dis_x, 
     &       dis_y, wi(maxn), nvx(dim,maxn),grap(dim,maxn),egrd(maxn)
      double precision x(dim,maxn),vx(dim,maxn),ddx(dim),dvx(dim,maxn),
     &       av(dim,maxn), maxvel, b, minvel, vel,sumvel, norp,
     &       kai,v_inf,cita,a,xl,dx, delta_r(dim,maxn)                          
     

      do i=1,maxn
        avdudt(i) = 0.
        ahdudt(i) = 0.
        sumw(i) = 0.
        pp(i) = 0.
        egrd(i) = 0.
        do  d=1,dim
          indvxdt(d,i) = 0.
          ardvxdt(d,i) = 0.
          exdvxdt(d,i) = 0.
          nvx(d,i) = 0.
          grap(d,i) = 0.
        enddo
      enddo  

      np = 31
      nnp = 10
      xl = x_maxgeom-x_mingeom
      dx = xl/mmp
      a  = tan(pi/3)   
      b = c0**2*1000/7
      v_inf = 0.
      kai = 0.
      cita = 0.
       maxvel = 0.e0
      minvel = 1.e1
      sumvel = 0.e0
c  velocity statistic
      do i =1,ntotal
        vel = sqrt(vx(1,i)**2+vx(2,i)**2)
         sumvel = sumvel + vel
         if (vel.gt.maxvel)then
            maxvel = vel
            maxi = i
         endif
         if (vel.lt.minvel)then
           minvel = vel
            mini = i
        endif
      enddo
c  generate virtual particles
      if(dynamic.or.dummy.or..not.mirror)then

         if(itimestep.eq.1) call virt_part(itimestep, ntotal,nvirt,hsml,
     &       mass,x,vx,rho,u,p,itype, nwall,mother)
          
       else
         call virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,
     &       rho,u,p,itype, nwall,mother)
       endif   
c      nwall = nwall   
      if(itimestep.eq.1)then
       open(13,file="../data/ini_virt.dat")
        do i = ntotal+1, ntotal+nvirt+nwall 
          write(13,1001) i, (x(d, i),d = 1, dim), p(i)
        enddo   
1001    format(1x, I5, 5(2x, e14.8)) 
        close(13)
      endif

      
c---  Interaction parameters, calculating neighboring particles
c     and optimzing smoothing length
      if (itimestep.eq.nstart+1.and.nstart.ne.0) then
         open(15,file="../data/xv_vp.dat")
         read(15,*) nvirt, nwall
         close(15)  
      endif
     
      ntotalvirt = ntotal + nvirt + nwall 
c      print *,ntotal,nvirt,nwall
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
c---  Density approximation or change rate
c---  summation_density:(4.26)
c---  con_density: calculting density through continuity equation (4.31)/(4.34)      
      if(dummy)then 
         do i = 1,nvirt
           vx(1,ntotal+i) = v_inf
           vx(2,ntotal+i) = v_inf
         enddo
      endif

      call con_density(ntotalvirt,mass,niac,pair_i,pair_j,
     &       hsml,w,dwdx,vx,itype,x,rho,wi,drho)
   
      if (dynamic) then
         do i = 1,ntotal+nvirt
           rho(i)=rho(i)+dt*drho(i)
           call p_art_water(rho(i),x(2,i),c(i),p(i))
         enddo
      else 
         do i = 1,ntotal
            rho(i)=rho(i)+dt*drho(i)
c          if (i.eq.75)  print *,drho(i)
            call p_art_water(rho(i),x(2,i),c(i),p(i))
         enddo   
      endif

      if (dummy) then
        do k = 1,niac
            i = pair_i(k)
            j = pair_j(k)
            if (itype(j).lt.0.and.itype(i).gt.0)then
              do d= 1,dim
                nvx(d,j) = nvx(d,j)+vx(d,i)*w(k)
              enddo
              sumw(j) = sumw(j)+w(k)
            endif
         enddo
         do i = ntotal+1,ntotal+nvirt
c           if(i.eq.ntotal+1) print *,nvx(1,i),sumw(i)
           do d = 1,dim
             if (sumw(i).ne.0)then
             vx(d,i) = 2*v_inf-nvx(d,i)/sumw(i)
             endif
           enddo
         enddo
      endif

      if(mirror) then
          if (nor_density.and.mod(itimestep,30).eq.0) then
            call sum_density(ntotal,hsml,mass,niac,pair_i,pair_j,w,
     &         itype,rho)
          endif
c   pressure correction as well as density     
c        b = c0**2*1000/7
        do i = ntotal+1,ntotal+nvirt
           p(i) = p(mother(i))
           rho(i) = rho(mother(i))
           if((x(2,i).lt.y_mingeom).or.(x(2,i).lt.a*x(1,i)-a*xl/2+nnp*dx
     &    .and.x(2,i).lt.-a*x(1,i)+a*xl/2+nnp*dx)) then
           p(i) = p(mother(i))+9.8*1000*(x(2,mother(i))-x(2,i))
           rho(i)= 1000*(p(i)/b+1)**(1/7)
           endif
         enddo
      endif
c  Shepard filter
       if (nor_density.and.mod(itimestep,30).eq.0) then
         call sum_density(ntotalvirt,hsml,mass,niac,pair_i,pair_j,w,
     &        itype,rho)
       endif

c---  Dynamic viscosity:

      if (visc) call viscosity(ntotalvirt,itype,x,rho,eta)
       
c---  Internal forces:(4.42)/(4.43)  4.58/4.59
 
      call int_force(itimestep,dt,ntotal,nvirt,hsml,mass,vx,niac,rho,
     &     eta, pair_i,pair_j,dwdx,u,itype,x,t,c,p,grap,indvxdt,
     &     tdsdt,du)

c          enddo
c       print *,grap(1,1),grap(2,1) 

c---  Artificial viscosity:(4.66)
      if (visc_artificial) call art_visc(ntotal,hsml,mass,x,vx,
     &      niac,rho,c,pair_i,pair_j,itype,w,dwdx,ardvxdt,avdudt)
   
      
c---  External forces:(4.93)
         if (ex_force) call ext_force(ntotal,nvirt,nwall,mass,x,vx,niac,
     &                  itype,pair_i,pair_j, hsml,maxvel, exdvxdt)

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
c---  Correction for dummy particles(pressure & density)
       if(dummy) then
        do i = 1,ntotal
          do d = 1,dim
          grap(d,i)=-grap(d,i)
c          if(d.eq.dim) grap(d,i)= grap(d,i)+9.8
          enddo
       enddo

        do k = 1, niac
           i = pair_i(k)
           j = pair_j(k)
           if (itype(j).lt.0.and.itype(i).gt.0) then

             pp(j) = pp(j) + p(i)*w(k)
             do d= 1,dim
               ddx(d) = x(d,j)-x(d,i)
             egrd(j) = egrd(j)+rho(i)*ddx(d)*grap(d,i)*w(k)
c  print the process of summation for debug
c             if(itimestep.eq.100.and.j.eq.ntotal+1)then 
c                print *,rho(i),dx(d),grap(d,i),w(k)    
              enddo
            endif
        enddo
c       if(mod(itimestep,print_step).eq.0)then
c         print *,pp(ntotal+1), egrd(ntotal+1),sumw(ntotal+1)
c         print *,

        do i = ntotal+1,ntotal+nvirt
          if(sumw(i).ne.0)then
                  
           p(i) = (pp(i)+egrd(i))/sumw(i)
c      background pressure   
c           kai = 1000*9.8*(y_maxgeom-x(2,i))
           rho(i) = 1000*((p(i)-kai)/b+1)**(1/7)
          endif
         enddo

c      print *,indvxdt(2,7555),indvxdt(2,7556)
c      print *,ardvxdt(2,7555),ardvxdt(2,7556)
       endif

      do i=1,ntotal 
       if(itype(i).gt.0) then
        do d=1,dim
          dvx(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i)   
c          if(i.eq.52.and.d.eq.dim) print *,indvxdt(d,i),exdvxdt(d,i)       
        enddo
c     gravity, damping technique(Adami,2012)
        if (self_gravity) then
           if(damp_t.ne.0.and.itimestep*dt.le.damp_t)then
               cita = 0.5*(sin((-0.5+itimestep*dt/damp_t)*pi)+1)
c     different samping function
c                cita = (itimestep*dt/damp_t)**3
c                dvx(dim,i) = dvx(dim,i)-9.8*cita
c            else if(itimestep*dt.gt.damp_t/2.and.
c     & itimestep*dt.le.damp_t)then
c                cita = 4*(itimestep*dt/damp_t-1)**3+1
                dvx(dim,i) = dvx(dim,i)-9.8*cita
            else
               dvx(dim,i) = dvx(dim,i)-9.8
            endif
         endif
c         if(abs(dvx(2,int(ntotal/2))).le.1e-7) print *,itimestep
          du(i) = du(i) + avdudt(i) + ahdudt(i)
           u(i) = u(i) + dt*du(i)
        if(u(i).lt.0)  u(i) = 0.         
        do d = 1, dim                   
          vx(d, i) = vx(d, i) + dt * dvx(d, i) + av(d, i)
c        if(i.eq.40) print *,dvx(d,i),av(d,i)
          x(d, i) = x(d, i) + dt * vx(d, i)       
        enddo
       endif
      enddo
c      if(abs(dvx(2,int(ntotal/2))).le.1e-6) print *,itimestep
      if (shifting) then
        call shift_position(ntotal,nvirt,ns,maxvel,pair_i,pair_j,niac,
     &   dt,x,delta_r)

        do i = 1,ntotal
            do d = 1,dim
               x(d,i) = x(d,i)+delta_r(d,i)
            enddo
            
            do d = 1,dim
               vx(d,i) = vx(d,i)+dvx(d,i)*delta_r(d,i)
               p(i) = p(i)+grap(d,i)*delta_r(d,i)
            enddo
            rho(i) = 1000*((p(i)-kai)/b+1)**(1/7)
         enddo

       endif
c      print *,p(75),wi(75)

c     keep the gate moving to a certain height
       if(gate.and.itimestep.le.2000)then
         do i = ntotal+nvirt-np*2+1,ntotal+nvirt
             x(2,i) = x(2,i) + dt*vx(2,i)
         enddo
      endif

c     output data of virtual particles
      if (mod(itimestep,save_step).eq.0) then
c       open(30,file="../data/trace_p.dat")
        open(40,file="../data/xv_vp.dat")
        open(50,file="../data/state_vp.dat")
        open(60,file="../data/other_vp.dat")            
        write(40,*) nvirt, nwall
        do i = ntotal + 1, ntotal + nvirt + nwall         
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
          write(*,*) '**** particle moving fastest ****', maxi, maxvel         
c          write(*,*) dvx(2,int(ntotal/2))
c          write(*,101)'velocity(y)','internal(y)','total(y)'   
c          write(*,100) x(1,maxi),x(2,maxi),vx(1,maxi),vx(2,maxi)
          write(*,*) '**** average velocity:', real(sumvel/ntotal)
c           write(*,*) '**** particle moving slowest ****', mini         
c         write(*,102)'velocity(y)','internal(y)','total(y)'   
c          write(*,103)  x(1,mini),x(2,mini),vx(1,mini),vx(2,mini) 
      endif
      
c100   format(1x,4(2x,e12.6))     
c101   format(1x,4(2x,e12.6))          
c103   format(1x,4(2x,e12.6))      

      end
