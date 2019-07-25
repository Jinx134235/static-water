      program SPH

c----------------------------------------------------------------------
c     This is a three dimensional SPH code. the followings are the 
c     basic parameters needed in this code or calculated by this code

c     mass-- mass of particles                                      [in]
c     ntotal-- total particle number ues                            [in]
c     dt--- Time step used in the time integration                  [in]
c     itype-- types of particles                                    [in]
c     x-- coordinates of particles                              [in/out]
c     vx-- velocities of particles                              [in/out]
c     rho-- densities of particles                              [in/out]
c     p-- pressure  of particles                                [in/out]
c     u-- internal energy of particles                          [in/out]
c     hsml-- smoothing lengths of particles                     [in/out]
c     c-- sound velocity of particles                              [out]
c     s-- entropy of particles                                     [out]
c     e-- total energy of particles                                [out]

      implicit none     
      include 'param.inc'

      integer ntotal, itype(maxn), maxtimestep, d, m, i, yesorno, 
     & current_ts, nstart,np    
      double precision x(3,maxn), vx(3, maxn), mass(maxn),rho(maxn),
     &     p(maxn), u(maxn), c(maxn), s(maxn), e(maxn), hsml(maxn), dt
      double precision s1, s2, c0, dx

      current_ts=0
      call time_print

c  this timestep is set according to dt<=c*h/c0, where c=0.2 is the crount
c  number, and c0 is the speed of sound      
     
      
      call input(x, vx, mass, rho, p, u, itype, hsml, ntotal,np,dx)     
c the speed of sound under water is calculated with c0=15sqrt(gh), which
c is also implied in the equation of state
      c0 = 15*sqrt(9.8*np*dx)
      dt = 0.2*hsml(1)/c0 
  1   write(*,*)'  ***************************************************' 
      write(*,*)'          Please input the maximal time steps        '
      write(*,*)'  ***************************************************'
      read(*,*) maxtimestep      
c      
      call time_integration(x, vx, mass, rho, p, u, c, s, e, itype, 
     &     hsml, ntotal, maxtimestep, dt, current_ts )
      call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)      
c      write(*,*)'  ***************************************************' 
c      write(*,*)'          Run more time steps?(1=yes, 0=no)        '
c      write(*,*)'  ***************************************************'
c      read(*,*) yesorno
c      if(yesorno.eq.1) go to 1
      call time_print
c      call time_elapsed(s2)      
c      write (*,*)'        Elapsed CPU time = ', s2-s1
                           
      end
