c---------------------------------------------------------
c     Including file for parameters and constants used 
c     in the entire SPH software packages.
c---------------------------------------------------------

c     dim : Dimension of the problem (1, 2 or 3)
      integer dim
      parameter ( dim = 2 )

c     maxn    : Maximum number of particles
c     max_interation : Maximum number of interaction pairs
      integer maxn,max_interaction
      parameter ( maxn    = 30000 ,
     &            max_interaction = 100 * maxn )

c     Parameters for the computational geometry,  
c     x_maxgeom : Upper limit of allowed x-regime 
c     x_mingeom : Lower limit of allowed x-regime 
c     y_maxgeom : Upper limit of allowed y-regime 
c     y_mingeom : Lower limit of allowed y-regime 
c     z_maxgeom : Upper limit of allowed z-regime 
c     z_mingeom : Lower limit of allowed z-regime 
      double precision x_maxgeom,x_mingeom,y_maxgeom,
     &                 y_mingeom,z_maxgeom,z_mingeom,
     &                 damheight,damlength,bedheight
      parameter ( x_maxgeom =  .39e0     ,
     &            x_mingeom =  0.e0     ,
     &            y_maxgeom =  .39e0     ,
     &            y_mingeom =  0.e0     ,
     &            z_maxgeom =  10.e0     ,
     &            z_mingeom = -10.e0     ,
     &            damheight =  .15e0    ,
     &            damlength =  .38e0    ,
     &            bedheight =  0.e0    )
    
c     SPH algorithm for particle approximation (pa_sph)
c     pa_sph = 1 : (e.g. (p(i)+p(j))/(rho(i)*rho(j))
c              2 : (e.g. (p(i)/rho(i)**2+p(j)/rho(j)**2)
c              3 : (e.g. (p(i)*rho(j)+p(j)*rho(i))/(rho(i)+rho(j))
      integer pa_sph 
      parameter(pa_sph = 1)

c     Nearest neighbor particle searching (nnps) method
c     nnps = 1 : Simplest and direct searching
c            2 : Sorting grid linked list
c            3 : Tree algorithm
      integer nnps 
      parameter(nnps = 2)

c     Smoothing length evolution (sle) algorithm
c     sle = 0 : Keep unchanged,
c           1 : h = fac * (m/rho)^(1/dim)
c           2 : dh/dt = (-1/dim)*(h/rho)*(drho/dt)
c           3 : Other approaches (e.g. h = h_0 * (rho_0/rho)**(1/dim) ) 

      integer sle 
      parameter(sle = 0)

c     Smoothing kernal
c     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
c         = 2, Gauss kernel   (Gingold and Monaghan 1981) 
c         = 3, Quintic kernel (Morris 1997)
      integer skf 
      parameter(skf = 1)

c     Density filter
c     filt = 0, no filter
c          = 1, Moltine & Clagrossi,2009
c          = 2, Georgios,2019
c          = 3, Antuono.etc,2010
      integer filt
      parameter(filt = 0)

c     Switches for different senarios

c     summation_density = .TRUE. : Use density summation model in the code, 
c                        .FALSE.: Use continuiity equation
c     average_velocity = .TRUE. : Monaghan treatment on average velocity,
c                       .FALSE.: No average treatment.
c     config_input = .TRUE. : Load initial configuration data,
c                   .FALSE.: Generate initial configuration.
c     virtual_part = .TRUE. : Use vritual particle,
c                   .FALSE.: No use of vritual particle.
c     vp_input = .TRUE. : Load virtual particle information,
c               .FALSE.: Generate virtual particle information.
c     visc = .true. : Consider viscosity,
c           .false.: No viscosity.
c     ex_force =.true. : Consider external force,
c               .false.: No external force.
c     visc_artificial = .true. : Consider artificial viscosity,
c                      .false.: No considering of artificial viscosity.
c     heat_artificial = .true. : Consider artificial heating,
c                      .false.: No considering of artificial heating.
c     self_gravity = .true. : Considering self_gravity,
c                    .false.: No considering of self_gravity
c     nor_density =  .true. : Density normalization by using CSPM,
c                    .false.: No normalization.
c     dynamic= .true. : apply dynamic boundary into simulation,Crespo ги2005
c     mirror = .true. : mirror boundary                      
c     gate = .true. : add moving gate in the dambreak simulation,only 
c                     for dynamic
c     dummy = .true. : dummy boundary

      logical summation_density, average_velocity, config_input,
     &        virtual_part, vp_input, visc, ex_force, heat_artificial,
     &        visc_artificial, self_gravity, nor_density, dynamic,
     &        mirror, gate, dummy

      parameter ( summation_density  = .false. )
      parameter ( average_velocity  = .true. )
      parameter ( config_input  = .false. )
      parameter ( virtual_part  = .true. )
      parameter ( vp_input  = .false.)
      parameter ( visc  = .false. )
      parameter ( ex_force  = .false.)
      parameter ( visc_artificial  = .true. )
      parameter ( heat_artificial  = .false. )
      parameter ( self_gravity  = .true. )      
      parameter ( nor_density  = .true. )      
      parameter ( dynamic = .false. )
      parameter ( gate = .false.  )
      parameter ( mirror = .true. )
      parameter ( dummy = .false. )

c     Symmetry of the problem
c     nsym = 0 : no symmetry,
c          = 1 : axis symmetry,
c          = 2 : center symmetry.     
      integer    nsym
      parameter ( nsym = 0)

c     Control parameters for output 
c     int_stat = .true. : Print statistics about SPH particle interactions.
c                        including virtual particle information.
c     print_step: Print Timestep (On Screen)
c     save_step : Save Timestep    (To Disk File)
c     moni_particle: The particle number for information monitoring.
      logical int_stat
      parameter ( int_stat = .true. )
      integer print_step, save_step, moni_particle
      parameter ( print_step = 100,
     &            save_step = 1 ) 

c    Some constant
c    c0: sound speed according to a*sqrt(g*h),a:10~20
c    damp_t: damping time in damping technique     
      double precision pi,c0,damp_t
      parameter ( pi = 3.14159265358979323846,
     &            c0 = 2.632425e1,   
     &            damp_t = 0.2 )

c    Switch for different benchmarks
      logical dambreak, static, cavity
      parameter ( dambreak = .false.)
      parameter ( static = .true.)
      parameter ( cavity = .false.)
