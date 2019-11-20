      subroutine time_integration(x,vx, mass, rho, p, u, c, s, e, itype, 
     &           hsml, ntotal, maxtimestep, dt, current_ts)
     
c----------------------------------------------------------------------
c    !!! this subroutine has been modified a lot
c    !!! original integration scheme is leap-frog, removed
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
c      ntotal-- fluid particle number                            [input]  
c      nvirt-- 
c      nwall-- wall particle number                       [input/output]
c      maxtimestep-- maximum timesteps                           [input]
c      dt-- timestep                                             [input]
   

      implicit none     
      include 'param.inc'
      
      integer itype(maxn), ntotal, maxtimestep, niac,np
      double precision x(dim, maxn), vx(dim, maxn), mass(maxn), 
     &       rho(maxn), p(maxn), u(maxn), c(maxn), s(maxn), e(maxn), 
     &       hsml(maxn), dt
      integer i, j, k, im,itimestep, d, current_ts, nstart, nvirt,nwall,
     &      pair_i(max_interaction), pair_j(max_interaction),scale_k       
      double precision  x_min(dim, maxn), v_min(dim, maxn), u_min(maxn),
     &       rho_min(maxn), dx(dim,maxn), dvx(dim, maxn), du(maxn),  
     &       drho(maxn),  av(dim, maxn), ds(maxn),dxiac(dim),hv(dim),
     &       t(maxn), tdsdt(maxn), temp_u, temp_rho, sumvel, maxvel 
      double precision  time, water_h, dis, xl, wij
      real :: bound(2,2)
    
      real, allocatable :: p_record(:),v_record(:,:),x_record(:,:),
     &        pro(:,:)
      integer, allocatable :: step(:)
      allocate(p_record(200))
      allocate(v_record(2,200))
      allocate(x_record(2,200))
      allocate(step(4))
c      common nvirt
      np = 60
      allocate(pro(4,100))
      water_h = 0.3      
      xl = x_maxgeom-x_mingeom
c      dx = xl/mmp
      bound(1,:) = (/x_mingeom,x_maxgeom/)
      bound(2,:) = (/y_mingeom,y_maxgeom/)  
c    timepoint for output and observation in dambreak     
      step(:) = (/815,1595,2410,3673/) 
      do i = 1, ntotal
        do d = 1, dim
          av(d, i) = 0.
        enddo
      enddo  
   
      do j = 1,np
        pro(1,j) = 0.
        pro(2,j) = 0.
        pro(3,j) = 0.
        pro(4,j) = 0.
      enddo
      
      if(skf.eq.1) then
          scale_k = 2
      elseif(skf.eq.2) then
          scale_k = 3
      elseif(skf.eq.3) then
          scale_k = 3
      endif

      nstart = current_ts
      do itimestep = nstart+1, nstart+maxtimestep   
	   
        current_ts=current_ts+1
c        time=current_ts*dt
        if (mod(itimestep,print_step).eq.0) then
         write(*,*)'______________________________________________'
         write(*,*)'  current number of time step =',
     &           itimestep,'     current time=', real(time+dt)
         write(*,*)'_____________________________________________'
        endif      
       
c---  Definition of variables out of the function vector:    

        call single_step(itimestep,nstart, dt, ntotal,nvirt,nwall, hsml,
     &        mass,x, vx,u, s, rho, p, t, tdsdt, du, ds,c, itype, av,
     &        niac, pair_i, pair_j, sumvel,maxvel)  
     
c      deal with those particles out of domain, maybe not work
c        do i =1,ntotal
c          do d  =1,dim
c           if(x(d,i).gt.bound(d,2).or.x(d,i).lt.bound(d,1)) then
c            itype(i) = 0
c             mass(i) = 0
c           endif
c          enddo
c        enddo
c       do i =1,4 
        if (mod(itimestep,print_step).eq.0) then
c              p_record(1) = itimestep
c         if(itimestep.eq.step(i)) then
             i = itimestep/print_step 
             p_record(i) = p(1)
             v_record(1,i) = sumvel/ntotal
             v_record(2,i) = maxvel
c          wavefront of dambreak
             x_record(1,i) = maxval(x(1,1:ntotal))/water_h
c          waterheight of dambreak          
             x_record(2,i) = maxval(x(2,1:ntotal))/water_h
         endif
c       enddo
c        print *,vx(2,1)
c        if (vx(2,1).gt.0) then
c            print *,itimestep
c        endif

c         if(dynamic) then
c           do i=1,nvirt
c                 rho(ntotal+i) = rho(ntotal+i) + dt*drho(ntotal+i)
c               call p_art_water(rho(ntotal+i),x(2,ntotal+i),c(ntotal+i),
c     &         p(ntotal+i))
c            enddo
c         endif
      
c---  Judge if the system achieves the balancing point      
c      if  (abs(p(1)-9.8e3*(y_maxgeom-x(2,1))).le.10) then
c           print *,p(1)
c           print *,current_ts    
c           stop
c      endif
        time = time + dt

        if(itimestep.eq.nstart+maxtimestep) then
            open(21,file="../data/record.dat")
            do i =1,(nstart+maxtimestep)/print_step
c                step = i*print_step
c            do i = 1,4
              write(21,1002) i,p_record(i),v_record(1,i),v_record(2,i),
     &    x_record(1,i),x_record(2,i)
            enddo
1002      format(1x, I6, 5(2x, e14.8))
          close(21)
                
c  extract pressure profile in wedge case, along the wedge
c                    |
c                    |
cc                   |
c                    |
cc                   |
c                    |
c                     \  
c                      \
c                       \
c                        \ 
c       
       if(geometry)then
        do j =1,np 
         if (j.le.20) then
           pro(1,j) =xl/2+.2/sqrt(.3)-j*xl/mmp/sqrt(.3)
         else
           pro(1,j) = xl/2
         endif
         pro(2,j) = j*xl/mmp
         enddo
         do i = 1,ntotal
           do j = 1,np
             dis = (x(1,i)-pro(1,j))**2+(x(2,i)-pro(2,j))**2
             dis = sqrt(dis)
             do d = 1,dim
                dxiac(d) = x(d,i)-pro(d,j)
             enddo
             if (dis.lt.scale_k*hsml(i))then
                call kernel(dis,dxiac,hsml(i),wij,hv)
                 pro(3,j) = pro(3,j)+p(i)*wij
                  pro(4,j) = pro(4,j)+wij
             endif
           enddo
          enddo

          do i = 1,np
            pro(3,i) = pro(3,i)/pro(4,i)
          enddo

          open(32,file="../data/pre_pro.dat")
            do i = 1,np
              write(32,1003) i,pro(1,i),pro(2,i),pro(3,i)
            enddo
1003      format(1x, I6, 3(2x, e14.8))
          close(32)
         endif
c  extract velocity profile in cavity case, central line in
c  both  x-&y-directions
         if(cavity) then
c  reassign profile points            
           do i = 1,mmp
              pro(1,i) = x_mingeom+(i-1)*xl/mmp
              pro(2,i) = xl/2
           enddo
           do i = mmp+1,2*mmp
              pro(1,i) = xl/2
              pro(2,i) = y_mingeom+(i-mmp)*xl/mmp
           enddo

           do i = 1,ntotal
             do j = 1,2*mmp
               dis = sqrt((x(1,i)-pro(1,j))**2+(x(2,i)-pro(2,j))**2)
               if(dis.lt.scale_k*hsml(i)) then               
                do d = 1,dim
                 dxiac(d) = x(d,i) - pro(d,j)
                enddo
                call kernel(dis,dxiac,hsml(i),wij,hv)
                  if(j.le.mmp) then
                     pro(3,j) = pro(3,j)+vx(2,i)*wij
                     pro(4,j) = pro(4,j)+wij
                   else
                     pro(3,j) = pro(3,j)+vx(1,i)*wij
                     pro(4,j) = pro(4,j)+wij
                   endif
                endif
             enddo
            enddo

            do i = 1,2*mmp
               pro(3,i) = pro(3,i)/pro(4,i)
            enddo
          open(43,file="../data/vel_pro.dat") 
            do i = 1,2*mmp
              write(43,1004) i,pro(1,i),pro(2,i),pro(3,i)
            enddo
1004      format(1x, I6, 3(2x, e14.8))
          close(43)
         endif
        endif

	if (mod(itimestep,save_step).eq.0) then
          call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)
c         print x
	endif 

c        if (mod(itimestep,print_step).eq.0) then
c          write(*,*)
c          write(*,101)'x','velocity', 'dvx'
c         write(*,100)x(2,ntotal), vx(2,ntotal), dvx(2,ntotal)
c        endif
        
c101     format(1x,3(2x,a12))	 
c100     format(1x,3(2x,e12.6))
	 
      enddo
c      print *,current_ts
      deallocate(pro)
      deallocate(step)
      deallocate(x_record)
      deallocate(v_record)
      deallocate(p_record)
      end
