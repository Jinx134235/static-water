      subroutine p_art_water(rho,h, c, p)
      
c----------------------------------------------------------------------
c   Artificial equation of state for the artificial compressibility 

c     rho    : Density                                              [in]
c     u      : Internal energy                                      [in]
c     p      : Pressure                                            [out]
c     c      : sound velocity                                      [out]
c     Equation of state for artificial compressibility   

      implicit none
      include 'param.inc'
      double precision rho, u, p, c, b, h
      double precision gamma, rho0

c     Artificial EOS, Form 1 (Monaghan, 1994) 
c     See Equ.(4.88)
       gamma=7.
       rho0=1000.       
       c = c0
       b = 100*9.8*rho0*damheight/gamma 
       p = b*((rho/rho0)**gamma-1)      
c       print *,p
c       c = 0.03
       
c     coeff in c fluctuates from 15 to 20

c     static pressure underwater 
c      p = rho*9.8*(1.e-3-h)

c     Artificial EOS, Form 2 (Morris, 1997)
c     See Equ.(4.89)
c      c = 0.01
c      p = c**2 * rho      
      
      end 
