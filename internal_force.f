      subroutine int_force(itimestep,dt,ntotal,nvirt,hsml,mass,vx,niac,
     &      rho,eta,pair_i,pair_j,dwdx,u,itype,x,t,c,p,grap,dvxdt,tdsdt,
     &      dedt)

c----------------------------------------------------------------------
c   Subroutine to calculate the internal forces on the right hand side 
c   of the Navier-Stokes equations, i.e. the pressure gradient and the
c   gradient of the viscous stress tensor, used by the time integration. 
c   Moreover the entropy production due to viscous dissipation, tds/dt, 
c   and the change of internal energy per mass, de/dt, are calculated. 
 
c     itimestep: Current timestep number                            [in]
c     dt     : Time step                                            [in]
c     ntotal : Number of particles                                  [in]
c     hsml   : Smoothing Length                                     [in]
c     mass   : Particle masses                                      [in]
c     vx     : Velocities of all particles                          [in]
c     niac   : Number of interaction pairs                          [in]
c     rho    : Density                                              [in]
c     eta    : Dynamic viscosity                                    [in]
c     pair_i : List of first partner of interaction pair            [in]
c     pair_j : List of second partner of interaction pair           [in]
c     dwdx   : Derivative of kernel with respect to x, y and z      [in]
c     itype  : Type of particle (material types)                    [in]
c     u      : Particle internal energy                             [in]
c     x      : Particle coordinates                                 [in]
c     itype  : Particle type                                        [in]
c     t      : Particle temperature                             [in/out]
c     c      : Particle sound speed                                [out]
c     p      : Particle pressure                                   [out]
c     dvxdt  : Acceleration with respect to x, y and z             [out] 
c     tdsdt  : Production of viscous entropy                       [out]
c     dedt   : Change of specific internal energy                  [out]
c     grap   : Pressure gradient of a single particle              [out]
      implicit none
      include 'param.inc'
      
     
      integer itimestep, ntotal,niac,pair_i(max_interaction),
     &        pair_j(max_interaction), itype(maxn),nvirt 
      double precision dt, hsml(maxn), mass(maxn), vx(dim,maxn),      
     &       rho(maxn), eta(maxn), dwdx(3,max_interaction), u(maxn),
     &       x(dim,maxn), t(maxn), c(maxn), p(maxn), dvxdt(dim,maxn),          
     &       tdsdt(maxn),dedt(maxn),wi(maxn),hv(dim),selfdens,
     &       w(max_interaction), grap(dim,maxn)
      integer i, j, k, d
      double precision  dvx(3), txx(maxn), tyy(maxn), tzz(maxn),
     &       txy(maxn), txz(maxn), tyz(maxn), vcc(maxn),
     &       hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij, r,
     &       dwdr,mp,meta,mvol

c     Initialization of shear tensor, velocity divergence, 
c     viscous energy, internal energy, acceleration 
c      print *,ntotal

      do i=1,ntotal      
        txx(i) = 0.e0
        tyy(i) = 0.e0
        tzz(i) = 0.e0
        txy(i) = 0.e0
        txz(i) = 0.e0
        tyz(i) = 0.e0
        vcc(i) = 0.e0
        tdsdt(i) = 0.e0
        dedt(i) = 0.e0
        do d=1,dim
          dvxdt(d,i) = 0.e0
          grap(d,i) = 0.e0
        enddo 
      enddo

c     Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c

      if (visc) then 
        do k=1,niac
          i = pair_i(k)
          j = pair_j(k)
          
          do d=1,dim
            dvx(d) = vx(d,j) - vx(d,i)
          enddo
c          if (i.eq.52) print *,dvx(1),dvx(2), dwdx(1,k),mass(j)
          if (dim.eq.1) then 
            hxx = 2.e0*dvx(1)*dwdx(1,k)        
          else if (dim.eq.2) then           
            hxx = 2.e0*dvx(1)*dwdx(1,k) -  dvx(2)*dwdx(2,k) 
            hxy = dvx(1)*dwdx(2,k) + dvx(2)*dwdx(1,k)
            hyy = 2.e0*dvx(2)*dwdx(2,k) - dvx(1)*dwdx(1,k)
          else if (dim.eq.3) then
            hxx = 2.e0*dvx(1)*dwdx(1,k) - dvx(2)*dwdx(2,k) 
     &                                  - dvx(3)*dwdx(3,k) 
            hxy = dvx(1)*dwdx(2,k) + dvx(2)*dwdx(1,k)
            hxz = dvx(1)*dwdx(3,k) + dvx(3)*dwdx(1,k)          
            hyy = 2.e0*dvx(2)*dwdx(2,k) - dvx(1)*dwdx(1,k)        
     &                                  - dvx(3)*dwdx(3,k)
            hyz = dvx(2)*dwdx(3,k) + dvx(3)*dwdx(2,k)
            hzz = 2.e0*dvx(3)*dwdx(3,k) - dvx(1)*dwdx(1,k)
     &                                  - dvx(2)*dwdx(2,k)
          endif                              
          hxx = 2.e0/3.e0*hxx
          hyy = 2.e0/3.e0*hyy
          hzz = 2.e0/3.e0*hzz
          if (dim.eq.1) then 
            txx(i) = txx(i) + mass(j)*hxx/rho(j)
            txx(j) = txx(j) + mass(i)*hxx/rho(i)                 
          else if (dim.eq.2) then           
            txx(i) = txx(i) + mass(j)*hxx/rho(j)
            txx(j) = txx(j) + mass(i)*hxx/rho(i)   
            txy(i) = txy(i) + mass(j)*hxy/rho(j)
            txy(j) = txy(j) + mass(i)*hxy/rho(i)            
            tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
            tyy(j) = tyy(j) + mass(i)*hyy/rho(i)          
          else if (dim.eq.3) then
            txx(i) = txx(i) + mass(j)*hxx/rho(j)
            txx(j) = txx(j) + mass(i)*hxx/rho(i)   
            txy(i) = txy(i) + mass(j)*hxy/rho(j)
            txy(j) = txy(j) + mass(i)*hxy/rho(i) 
            txz(i) = txz(i) + mass(j)*hxz/rho(j)
            txz(j) = txz(j) + mass(i)*hxz/rho(i)                     
            tyy(i) = tyy(i) + mass(j)*hyy/rho(j)
            tyy(j) = tyy(j) + mass(i)*hyy/rho(i)
            tyz(i) = tyz(i) + mass(j)*hyz/rho(j)
            tyz(j) = tyz(j) + mass(i)*hyz/rho(i)   
            tzz(i) = tzz(i) + mass(j)*hzz/rho(j)
            tzz(j) = tzz(j) + mass(i)*hzz/rho(i)                 
          endif                              

c     Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz:

         hvcc = 0.
         do d=1,dim
           hvcc = hvcc + dvx(d)*dwdx(d,k)
         enddo
         vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
         vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
        enddo
      endif   

c      print *,txx(52),txy(52),tyy(52)
      do i=1,ntotal
      
c     Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab

        if (visc) then 
          if (dim.eq.1) then 
            tdsdt(i) = txx(i)*txx(i)                             
          else if (dim.eq.2) then           
            tdsdt(i) = txx(i)*txx(i) + 2.e0*txy(i)*txy(i)  
     &                               + tyy(i)*tyy(i) 
          else if (dim.eq.3) then
            tdsdt(i) = txx(i)*txx(i) + 2.e0*txy(i)*txy(i)  
     &                               + 2.e0*txz(i)*txz(i)
     &               + tyy(i)*tyy(i) + 2.e0*tyz(i)*tyz(i) 
     &                               + tzz(i)*tzz(i)
          endif   
          tdsdt(i) = 0.5e0*eta(i)/rho(i)*tdsdt(i)
        endif  

c        call p_art_water(rho(i),x(2,i),p(i),c(i))
      enddo
   


c      Calculate SPH sum for pressure force -p,a/rho
c      and viscous force (eta Tab),b/rho
c      and the internal energy change de/dt due to -p/rho vc,c
c      do d=1,dim
cc        hv(d) = 0.e0
c      enddo
       r=0. 
   
c       do i=1,ntotal
c              call kernel(r,hv,hsml(i),selfdens,hv)
c                wi(i)=selfdens*mass(i)/rho(i)
c       enddo

     
c  change the momentum equation when the particle has approached the boundary
c  applying the formulae put forward by A.J.Crespo(2007)               
           
c           if (itype(i).lt.0.or.itype(j).lt.0) then
               
c              he = he + (vx(d,j) - vx(d,i))*h
c              do d=1,dim
c              dvxdt(d,i)=dvxdt(d,i)-2*c(i)**2*w(k)/((w(k)+wi(i))**2)
c     &*dwdx(d,k)
c              dvxdt(d,j)=dvxdt(d,j)+2*c(j)**2*w(k)/((w(k)+wi(j))**2)
c     &*dwdx(d,k)
c                if(d.eq.dim.and.self_gravity)then
c                    dvxdt(d,i)=dvxdt(d,i)-9.8
c                    dvxdt(d,j)=dvxdt(d,j)-9.8
c                endif
c              enddo         
c              he = he*rhoij
c             dedt(i) = dedt(i) + mass(j)*he
c             dedt(j) = dedt(j) + mass(i)*he                  
c            else
           
c      if(pa_sph.eq.1) then
c       call dividition(rho, dwdx,pair_i,pair_j, eta, mass,p,vx, txx,  
c     &   txy, tyz, tyy, txz, tzz,dvxdt, dedt)
c        else
      do k=1,niac
           i = pair_i(k)
           j = pair_j(k)
           he = 0.e0
         if(pa_sph.eq.1)then
          do d=1,dim                
            h = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(d,k) 
c     Pressure gradient term
            grap(d,i) = grap(d,i) + mass(j)*h
            grap(d,j) = grap(d,j) + mass(i)*h
c     Viscous force

          if (visc) then             
             if (d.eq.1) then
                       
c     x-coordinate of acceleration

               h = h + (eta(i)*txx(i)/rho(i)**2 +
     &                  eta(j)*txx(j)/rho(j)**2)*dwdx(1,k)
               if (dim.ge.2) then
                 h = h + (eta(i)*txy(i)/rho(i)**2 + 
     &                    eta(j)*txy(j)/rho(j)**2)*dwdx(2,k)
                 if (dim.eq.3) then
                   h = h + (eta(i)*txz(i)/rho(i)**2 + 
     &                      eta(j)*txz(j)/rho(j)**2)*dwdx(3,k)
                 endif
               endif            
             elseif (d.eq.2) then
            
c     y-coordinate of acceleration

               h = h + (eta(i)*txy(i)/rho(i)**2  
     &               +  eta(j)*txy(j)/rho(j)**2)*dwdx(1,k)
     &               + (eta(i)*tyy(i)/rho(i)**2  
     &               +  eta(j)*tyy(j)/rho(j)**2)*dwdx(2,k)
               if (dim.eq.3) then
                 h = h + (eta(i)*tyz(i)/rho(i)**2  
     &                 +  eta(j)*tyz(j)/rho(j)**2)*dwdx(3,k)
               endif              
             elseif (d.eq.3) then
            
c     z-coordinate of acceleration

               h = h + (eta(i)*txz(i)/rho(i)**2 + 
     &                  eta(j)*txz(j)/rho(j)**2)*dwdx(1,k)
     &               + (eta(i)*tyz(i)/rho(i)**2 + 
     &                  eta(j)*tyz(j)/rho(j)**2)*dwdx(2,k)
     &               + (eta(i)*tzz(i)/rho(i)**2 + 
     &                  eta(j)*tzz(j)/rho(j)**2)*dwdx(3,k)            
             endif            
           endif              
           dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
           dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
          enddo
c        elseif(pa_sph.eq.2) then
c         call dividition(rho, dwdx,pair_i, pair_j, eta, mass,p, vx,
c     &    txx, txy, tyz, tyy, txz, tzz,dvxdt, dedt)

        elseif(pa_sph.eq.3) then
          mp = (rho(i)*p(j)+rho(j)*p(i))/(rho(i)+rho(j))
          meta = 2*eta(i)*eta(j)/(eta(i)+eta(j))
          mvol = (mass(i)/rho(i))**2+(mass(j)/rho(j))**2
          do d = 1,dim
            r = r + (x(d,i)-x(d,j))**2
            dwdr = dwdr + dwdx(d,k)*(x(d,j)-x(d,i))
          enddo
         
          do d = 1,dim
            h = -mvol*mp*dwdx(d,k)
            grap(d,i) = grap(d,i)+h/mass(i)
            grap(d,j) = grap(d,j)-h/mass(j)
          if(visc)then
            dvx(d) = vx(d,j)-vx(d,i)
            h = h+meta*mvol*dvx(d)*dwdr/r
           endif
           dvxdt(d,i) = dvxdt(d,i) + h/mass(i)
           dvxdt(d,j) = dvxdt(d,j) - h/mass(j)
           enddo
         endif
c    Energy equation         
         do d = 1,dim
          he = he + (vx(d,j) - vx(d,i))*h
         enddo
          dedt(i) = dedt(i) + mass(j)*he
          dedt(j) = dedt(j) + mass(i)*he
      enddo  

       if(pa_sph.eq.2) then
         call dividition(rho, dwdx,pair_i, pair_j, eta, mass,p, vx,
     &    txx, txy, tyz, tyy, txz, tzz,dvxdt, dedt) 
       endif

c     Change of specific internal energy de/dt = T ds/dt - p/rho vc,c:

      do i=1,ntotal
c         grap(i) = dvxdt(dim,i)*rho(i)
         dedt(i) = tdsdt(i) + 0.5e0*dedt(i)
c       if(i.eq.40) print *,dvxdt(2,i), dedt(i)
      enddo

      end

c----------------------------------------------------------------------
      subroutine dividition(rho, dwdx,pair_i, pair_j, eta, mass,p, vx,
     &    txx, txy, tyz, tyy, txz, tzz,dvxdt, dedt)

      implicit none
      include 'param.inc'
      double precision rho(maxn),dwdx(3,max_interaction),vx(dim,maxn),
     &               mass(maxn),txx(maxn),txy(maxn),tyy(maxn),tyz(maxn),
     &               txz(maxn),tzz(maxn),dvxdt(dim,maxn),dedt(maxn),
     &               h, he, rhoij,eta(maxn),p(maxn)
      integer i,j,k,d,niac,pair_i(max_interaction),
     &         pair_j(max_interaction)


      do k = 1,niac
         i=pair_i(k)
         j=pair_j(k)
         he = 0.e0
         rhoij = 1.e0/(rho(i)*rho(j))        
             
         do d=1,dim
            h = -(p(i) + p(j))*dwdx(d,k)
            he = he + (vx(d,i)-vx(d,j))*h
            if (visc) then             
              if (d.eq.1) then
                  h = h + (eta(i)*txx(i) + eta(j)*txx(j))*dwdx(1,k)
                if (dim.ge.2) then
                    h = h + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(2,k)
                 if (dim.eq.3) then
                    h = h + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(3,k)
                 endif
                endif            
               elseif (d.eq.2) then
           
                 h = h + (eta(i)*txy(i) + eta(j)*txy(j))*dwdx(1,k)
     &           + (eta(i)*tyy(i) + eta(j)*tyy(j))*dwdx(2,k)
                 if (dim.eq.3) then
                   h = h + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(3,k)
                 endif             
               elseif (d.eq.3) then
 
                 h = h + (eta(i)*txz(i) + eta(j)*txz(j))*dwdx(1,k)
     &                 + (eta(i)*tyz(i) + eta(j)*tyz(j))*dwdx(2,k)
     &                 + (eta(i)*tzz(i) + eta(j)*tzz(j))*dwdx(3,k)            
               endif
              endif             
             h = h*rhoij
             dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
             dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
          enddo
          he = he*rhoij
          dedt(i) = dedt(i) + mass(j)*he
          dedt(j) = dedt(j) + mass(i)*he                  
      enddo


      end

c----------------------------------------------------------------------
      subroutine plusition(rho, dwdx,pair_i, pair_j, eta, mass,p, vx,
     &    txx, txy, tyz, tyy, txz, tzz,dvxdt, dedt)

      implicit none
      include 'param.inc'
      double precision rho(maxn),dwdx(3,max_interaction),vx(dim,maxn),
     &               mass(maxn),txx(maxn),txy(maxn),tyy(maxn),tyz(maxn),
     &               txz(maxn),tzz(maxn),dvxdt(dim,maxn),dedt(maxn),
     &               h, he, rhoij,eta(maxn),p(maxn)
      integer i,j,k,d,niac,ntotal,pair_i(max_interaction),
     &         pair_j(max_interaction)

      do k=1,niac
           i = pair_i(k)
           j = pair_j(k)
c          print *,i,j
           he = 0.e0
          do d=1,dim                
            h = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(d,k) 
c            print *,h
            he = he + (vx(d,j) - vx(d,i))*h

c     Viscous force

          if (visc) then             
             if (d.eq.1) then
                       
c     x-coordinate of acceleration

               h = h + (eta(i)*txx(i)/rho(i)**2 +
     &                  eta(j)*txx(j)/rho(j)**2)*dwdx(1,k)
               if (dim.ge.2) then
                 h = h + (eta(i)*txy(i)/rho(i)**2 + 
     &                    eta(j)*txy(j)/rho(j)**2)*dwdx(2,k)
                 if (dim.eq.3) then
                   h = h + (eta(i)*txz(i)/rho(i)**2 + 
     &                      eta(j)*txz(j)/rho(j)**2)*dwdx(3,k)
                 endif
               endif            
             elseif (d.eq.2) then
            
c     y-coordinate of acceleration

               h = h + (eta(i)*txy(i)/rho(i)**2  
     &               +  eta(j)*txy(j)/rho(j)**2)*dwdx(1,k)
     &               + (eta(i)*tyy(i)/rho(i)**2  
     &               +  eta(j)*tyy(j)/rho(j)**2)*dwdx(2,k)
               if (dim.eq.3) then
                 h = h + (eta(i)*tyz(i)/rho(i)**2  
     &                 +  eta(j)*tyz(j)/rho(j)**2)*dwdx(3,k)
               endif              
             elseif (d.eq.3) then
            
c     z-coordinate of acceleration

               h = h + (eta(i)*txz(i)/rho(i)**2 + 
     &                  eta(j)*txz(j)/rho(j)**2)*dwdx(1,k)
     &               + (eta(i)*tyz(i)/rho(i)**2 + 
     &                  eta(j)*tyz(j)/rho(j)**2)*dwdx(2,k)
     &               + (eta(i)*tzz(i)/rho(i)**2 + 
     &                  eta(j)*tzz(j)/rho(j)**2)*dwdx(3,k)            
             endif            
           endif              
           dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
           dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
          enddo
          dedt(i) = dedt(i) + mass(j)*he
          dedt(j) = dedt(j) + mass(i)*he       
        enddo

c       do i=1,ntotal
         
c        enddo

        end
