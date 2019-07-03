      subroutine time_integration(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt)
     
c----------------------------------------------------------------------
c      x-- coordinates of particles                       [input/output]
c      vx-- velocities of particles                       [input/output]
c      mass-- mass of particles                                  [input]
c      rho-- dnesities of particles                       [input/output]
c      p-- pressure  of particles                         [input/output]
c      u-- internal energy of particles                   [input/output]
c      c-- sound velocity of particles                          [output]
c      s-- entropy of particles, not used here                  [output]
c      e-- total energy of particles                            [output]
c      itype-- types of particles                               [input]
c           =1   ideal gas
c           =2   water
c           =3   tnt
c      hsml-- smoothing lengths of particles              [input/output]
c      ntotal-- total particle number                            [input]  
c      maxtimestep-- maximum timesteps                           [input]
c      dt-- timestep                                             [input]
   

      implicit none     
      include 'param.inc'
      
      integer itype(maxn), ntotal, maxtimestep, niac
      double precision x(3, maxn), vx(3, maxn), mass(maxn), 
     &       rho(maxn), p(maxn), u(maxn), c(maxn), s(maxn), e(maxn), 
     &       hsml(maxn), dt
      integer i, j, k, im,itimestep, d, current_ts, nstart, nvirt,
     &      pair_i(max_interaction),   pair_j(max_interaction)       
      double precision  x_min(dim, maxn), v_min(dim, maxn), u_min(maxn),
     &       rho_min(maxn), dx(3,maxn), dvx(3, maxn), du(maxn),  
     &       drho(maxn),  av(3, maxn), ds(maxn),
     &       t(maxn), tdsdt(maxn), temp_u, temp_rho 
      double precision  time
               
      do i = 1, ntotal
        do d = 1, dim
          av(d, i) = 0.
        enddo
      enddo  
     
      nvirt = 0

       

      do itimestep = 1, maxtimestep   
	   
        current_ts=current_ts+1
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',
     &           itimestep,'     current time=', real(time+dt)
         write(*,*)'______________________________________________'
        endif      
       

c---  Definition of variables out of the function vector:    
      
        call single_step(itimestep, dt, ntotal,nvirt, hsml, mass, x, vx,
     &        u, s, rho, p, t, tdsdt, dx,dvx, du, ds, drho, itype, av,
     &        niac, pair_i, pair_j)  
                  
        
         if(dynamic) then
              do i=1,nvirt
                 rho(ntotal+i) = rho(ntotal+i) + dt*drho(ntotal+i)
               call p_art_water(rho(ntotal+i),x(2,ntotal+i),c(ntotal+i),
     &         p(ntotal+i))
         enddo

         endif
        
     
        time = time + dt

	if (mod(itimestep,save_step).eq.0) then
          call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)
c         print x
	endif 

        if (mod(itimestep,print_step).eq.0) then
          write(*,*)
          write(*,101)'x','velocity', 'dvx'
c      
         write(*,100)x(2,ntotal), vx(2,ntotal), dvx(2,ntotal)
c         enddo    
        endif
        
101     format(1x,3(2x,a12))	 
100     format(1x,3(2x,e12.6))
	 
      enddo

c--      nstart=current_ts

      end
