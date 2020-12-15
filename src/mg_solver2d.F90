! c
! c     file tmud2cr.f
! c
! c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! c  .                                                             .
! c  .                  copyright (c) 1999 by UCAR                 .
! c  .                                                             .
! c  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
! c  .                                                             .
! c  .                      all rights reserved                    .
! c  .                                                             .
! c  .                                                             .
! c  .                      MUDPACK version 5.0                    .
! c  .                                                             .
! c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! c
! c
! c ... author and specialist
! c
! c          John C. Adams (National Center for Atmospheric Research)
! c          email: johnad@ucar.edu, phone: 303-497-1213
! 
! c ... For MUDPACK 5.0 information, visit the website:
! c     (http://www.scd.ucar.edu/css/software/mudpack)
! c
! c
! c ... purpose
! c
! c     test program for the mudpack solver mud2cr
! c
! c ... required MUDPACK files
! c
! c     mud2cr.f, mudcom.f
! c
! c *********************************************************
! c *********************************************************
! c
! c     sample program/test driver for mud2cr
! c
! c **********************************************************
! c **********************************************************
! c
! c
! c     a sample program/test driver for mud2cr is listed below.  it
! c     can be executed as an initial test.  the output is listed
! c     for the test case described.
! c
! c     test mud2cr below by solving the nonseparable elliptic pde
! c     with cross derivative term
! c
! c          (1.+y**2)*pxx + (1.+x**2)*pyy + 2.*x*y*pxy +
! c
! c          y*px + x*py - (x*y)*pe = r(x,y)
! c
! c     on a grid as close to 50 by 64 as the mudpack size constraints
! c     allow.  the solution region is the unit square.  assume a
! c     mixed derivative boundary condition at y=1 of the form
! c
! c          -x * dp/dx + (1+x) * dp/dy - x * pe = gbdyd(x).
! c
! c     and specified (Dirchlet) boundary conditions elsewhere.  the
! c     exact solution
! c
! c          p(x,y) = (x*y)**5
! c
! c     is used to set the right hand side, boundary conditions, and
! c     compute the error.
! c
! c     red/black gauss-seidel point relaxation is used along with the
! c     the default multigrid options.  one full multigrid cycle reaches
! c     discretization level error for this problem.
! c
! c
! c
  subroutine MULTIGRID_MUDPACK(phi)

  use ntypes
  use Domain, only: iixp, jjyq, iiex, jjey, nnx, nny, llwork, ny, nz
  
      implicit none
! c
! c     set grid size params
! c
! c      integer iixp,jjyq,iiex,jjey,nnx,nny,llwork
! c      parameter (iixp = 8 , jjyq = 4, iiex = 8, jjey = 7 )
! c      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
! c
! c     estimate work space for point relaxation (see mud2cr.d)
! c
! c      parameter (llwork=(7*(nnx+2)*(nny+2)+44*nnx*nny)/3 )

      real*8 phi(nnx,nny),rhs(nnx,nny),work(llwork)
! c
! c     put integer and floating point argument names in contiguous
! c     storeage for labelling in vectors iprm,fprm
! c
      integer iprm(16),mgopt(4)
      real*8 fprm(6)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny_t, & 
                   iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny_t, & 
                   iguess,maxcy,method,nwork,lwrkqd,itero
      real*8 xa,xb,yc,yd,tolmax,relmax
      common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer ierror
      real*8 dlx,dly
! c
! c     declare coefficient and boundary condition input subroutines external
! c
      external cofcr,bndcr
      
! !       call  JACOBI_TRANS
! c
! c     set input integer arguments
! c
      intl = 0
! c
! c     set boundary condition flags
! c
!        nxa = 0
!        nxb = 0
!        nyc = 1
!        nyd = 1
! c
! c     set grid sizes from parameter statements
! c
      ixp = iixp
      jyq = jjyq
      iex = iiex
      jey = jjey
      nx = nnx
      ny_t = nny

! c   checking the dimensions     

      IF ( ((nz+2) .NE. nx) .OR. ((ny+2) .NE. ny_t) ) THEN
	STOP 'Mudpack Multigrid Parameters are not Matching Grid Dimensions.'
      ENDIF
          
! c
! c     set multigrid arguments (w(2,1) cycling with fully weighted
! c     residual restriction and cubic prolongation)
! c
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
! c
! c     set for one cycle
! c
      maxcy = 6 
! c
! c     set no initial guess forcing full multigrid cycling
! c
      iguess = 0
! c
! c     set work space length approximation from parameter statement
! c
      nwork = llwork
! c
! c     set point relaxation
! c
      method = 0
! c
! c     set end points of solution rectangle in (x,y) space
! c
!      if(nxa .eq. 0 .and. nxb .eq. 0) then
!       xa = 0.0; xb = nnx-2
!      else
       xa = 0.0; xb = nnx-1
!      endif
      yc = 0.0
      yd = nny-1
! c
! c     set mesh increments
! c
!      if(nxa .eq. 0 .and. nxb .eq. 0) then
!       dlx = (xb-xa)/float(nx-2)
!      else 
       dlx = (xb-xa)/float(nx-1)
!     endif
      dly = (yd-yc)/float(ny_t-1)
! c
! c     set for no error control flag
! c
      tolmax = 10.0**(-6)
  
      phi=0.0d0
      call phi_bc(rhs,nxa,nxb,nyc,nyd) !later boundary mj
      


      call mud2cr(iprm,fprm,work,cofcr,bndcr,rhs,phi,mgopt,ierror)
   
       intl = 1

       call mud2cr(iprm,fprm,work,cofcr,bndcr,rhs,phi,mgopt,ierror)

      end

      
      
      
       subroutine bndcr(kbdy,xory,alfa,beta,gama,gbdy)
! c
! c     input mixed "oblique" derivative b.c. to mud2cr
! c     at upper y boundary
! c
      implicit none
      integer kbdy
      real*8 xory,alfa,beta,gama,gbdy
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny_t, & 
                   iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny_t, & 
                   iguess,maxcy,method,nwork,lwrkqd,itero
      real*8 xa,xb,yc,yd,tolmax,relmax
      common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
      real*8 x,y
      
      if (kbdy.eq.3) then
! c
! c     y=yd boundary (nyd must equal 2 if this code is to be executed).
! c     b.c. has the form alfyd(x)*px+betyd(x)*py+gamyd(x)*pe = gbdyd(x)
! c     where x = yorx.   alfa,beta,gama,gbdy corresponding to alfyd(x),
! c     betyd(x),gamyd(x),gbdyd(y) must be output.
! c
      y = yd
      x = xory
      alfa = 0
      beta = 1
      gama = 0
!      call exacr(x,y,pxx,pxy,pyy,px,py,pe)
      gbdy = 0.0d0 !alfa*px + beta*py + gama*pe
      return
      
      
      else if (kbdy.eq.4) then
! c
! c     y=yd boundary (nyd must equal 2 if this code is to be executed).
! c     b.c. has the form alfyd(x)*px+betyd(x)*py+gamyd(x)*pe = gbdyd(x)
! c     where x = yorx.   alfa,beta,gama,gbdy corresponding to alfyd(x),
! c     betyd(x),gamyd(x),gbdyd(y) must be output.
! c
      y = yd
      x = xory
      alfa = 0
      beta = 1
      gama = 0
!      call exacr(x,y,pxx,pxy,pyy,px,py,pe)
      gbdy = 0.0d0 !alfa*px + beta*py + gama*pe
      return
      end if
      end
