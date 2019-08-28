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
      double precision s1, s2, dx

      current_ts=0
      call time_print
      
      call input(x, vx, mass, rho, p, u, itype, hsml, ntotal)     
c     timestep according to CFL condition, crount number =0.2
      dt = 7.8e-5
      write(*,*)'timestep(s):', real(dt)
  1   write(*,*)'  ***************************************************' 
      write(*,*)'          Please input the maximal time steps        '
      write(*,*)'  ***************************************************'
      read(*,*) maxtimestep      
c      
      call time_integration(x, vx, mass, rho, p, u, c, s, e, itype, 
     &     hsml, ntotal, maxtimestep, dt, current_ts )
      call output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal)      
      call time_print
      write(*,*)'  ***************************************************' 
      write(*,*)'          Run more time steps?(1=yes, 0=no)        '
      write(*,*)'  ***************************************************'
      read(*,*) yesorno
      if(yesorno.eq.1) go to 1
                           
      end
