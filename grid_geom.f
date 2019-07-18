      subroutine grid_geom(i,x,hsml,scale_k,maxgridx,mingridx,xgcell,
     &             itimestep)

c----------------------------------------------------------------------
c   Subroutine to calculate the coordinates (xgcell) of the cell of 
c   the sorting  grid, in which the particle with coordinates (x) lies.

c     x        : Coordinates of particle                            [in]    
c     ngridx   : Number of sorting grid cells in x, y, z-direction  [in]
c     maxgridx : Maximum x-, y- and z-coordinate of grid range      [in]
c     mingridx : Minimum x-, y- and z-coordinate of grid range      [in]
c     dgeomx   : x-, y- and z-expansion of grid range               [in]
c     xgcell   : x-, y- and z-coordinte of sorting grid cell       [out]

      implicit none
      include 'param.inc'

      integer i,xgcell(3), scale_k, itimestep
      double precision x(dim), maxgridx(dim), mingridx(dim),
     &    cll, hsml
      integer d

      do d=1,3
        xgcell(d) = 1
      enddo

      cll = scale_k*hsml
c     print *,cll
      do d=1,dim
        if ((x(d).gt.maxgridx(d)).or.(x(d).lt.mingridx(d))) then
          print *,' >>> ERROR <<< : Particle out of range'
          print *,'    Particle position: x(',i,d,') = ',x(d)
          print *,'    Range: [xmin,xmax](',D,') = 
     &         [',mingridx(d),',',maxgridx(d),']'
          print *,'current timestep:', itimestep
          stop
        else
          xgcell(d) = int((x(d)-mingridx(d))/cll + 1.e0)
        endif
      enddo

c      print *,xgcell(2)

      end