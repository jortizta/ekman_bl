      subroutine mgd9v(levels, nxc, nyc, nxf, nyf, nm,
     +                 iout, istart, iprep, maxit, tol,
     +                 rhs, rhsc, a, ac, v, vc, vb, vbc, work,
     +                 ldu, lduc, wa, wac, wb, wbc, resno, ndid, ifail)
      integer levels, nxc, nyc, nxf, nyf, nm, iout(6), istart, iprep,
     +        maxit, ndid, ifail
      double precision
     +     tol, rhs(nxf*nyf), rhsc(nm), a(nxf*nyf*9), ac(nm*9),
     +     v(nxf*nyf), vc(nm), vb(nxf*nyf), vbc(nm), work(nxf*12),
     +     ldu(nxf*nyf*3), lduc(nm*3), wa(nxf*nyf), wac(nm),
     +     wb(nxf*nyf), wbc(nm), resno
c-----------------------------------------------------------------------
c123456789012345678901234567890123456789012345678901234567890123456789
c
c     ******************************   ******************************
c    ******************************** ********************************
c   *******************************************************************
c   *******                         ***                         *******
c   ******* multigrid program mgd9v *** multigrid program mgd9v *******
c   *******  version june 6, 1994   *** version june 6, 1994    *******
c   *******                         ***                         *******
c   *******************************************************************
c    ******************************** ********************************
c     ******************************   ******************************
c
c     purpose
c     -------
c
c     this program solves a user provided linear system of 9-point
c     difference equations on a rectangular (curvulinear) grid.
c
c     mathematical method
c     -------------------
c
c     sawtooth multigrid cycling
c     (i.e. one smoothing-sweep after each coarse grid correction)
c
c     with
c
c          smoothing by incomplete line lu-decomposition,
c          weighted 9-point prolongation and restriction,
c          galerkin approximation of coarse grid matrices.
c
c     comment
c     -------
c
c -   this version is written in ansi fortran 77, however only a subset
c     of the full fortran 77 language has been used.
c -   for an efficient performance it is advisable to take the longest
c     side of a rectangular grid along the x-axis.
c -   double precision type is very easily changed into double precision
c     ( but check the parameter statements )
c -   compared with previous versions of mgd9v, one notes that now
c     data corresponding to the finest grid (e.g. the matrix)
c     are stored in arrays different from the ones that contain data
c     of the coarser grids. this enhances flexibility of use of mgd9v
c     as a module in enveloping codes
c
c     author and programmer
c     ---------------------
c
c         p.m. de zeeuw
c         department of numerical mathematics  (nw)
c         the centre for mathematics and computer science  (cwi)
c         kruislaan 413
c         1098 sj  amsterdam
c         the netherlands
c
c     references
c     ----------
c
c     1.  p.m. de zeeuw
c         matrix-dependent prolongations and restrictions in a
c         blackbox multigrid solver
c         journal of computational and applied mathematics
c         vol. 33 (1990) pp.1-27.
c     2.  p.m. de zeeuw
c         multigrid and advection
c         c.b. vreugdenhil, b. koren, (eds.)
c         numerical methods for advection-diffusion problems
c         notes on numerical fluid mechanics, vol. 45,
c         vieweg verlag, braunschweig, germany (1993)
c         issn 0179-9614, isbn 3-528-07645-3.
c
c     Copyright (C) Stichting Mathematisch Centrum, Amsterdam
c***********************************************************************
c
c                  ****   parameters   ****
c
c***********************************************************************
c ---
c     levels       integer.
c                  on entry, number of levels in multigrid method,
c                  should be .ge.1 and .le.12
c     nxc,nyc      integer.
c                  on entry, number of vertical, horizontal grid-lines
c                  on coarsest grid, nxc should be .ge.3
c                                and nyc should be .ge.3
c     nxf,nyf      integer.
c                  on entry, number of vertical, horizontal grid-lines
c                  on finest grid.
c                  the following relations should hold between
c                   levels, nxc, nxf  and between  levels, nyc, nyf:
c                   nxf = 2**(levels-1)*(nxc-1) + 1,
c                   nyf = 2**(levels-1)*(nyc-1) + 1.
c     nm           integer.
c                  on entry,
c                   if levels equals 1 then 1,
c                   else the number of all grid-points in the x,y-plane
c                   on all coarser grids together.
c
c                  examples
c                  --------
c
c                  levels =   1    2    3     4     5      6      7
c                  nxc    =   5    5    5     5     5      5      5
c                  nyc    =   5    5    5     5     5      5      5
c                  nxf    =   5    9   17    33    65    129    257
c                  nyf    =   5    9   17    33    65    129    257
c                  nm     =   1   25  106   395  1484   5709  22350
c
c
c                  levels =   1    2    3     4     5      6      7
c                  nxc    =   5    5    5     5     5      5      5
c                  nyc    =   3    3    3     3     3      3      3
c                  nxf    =   5    9   17    33    65    129    257
c                  nyf    =   3    5    9    17    33     65    129
c                  nm     =   1   15   60   213   774   2919  11304
c
c
c                  levels =   1    2    3     4     5      6      7
c                  nxc    =   6    6    6     6     6      6      6
c                  nyc    =   6    6    6     6     6      6      6
c                  nxf    =   6   11   21    41    81    161    321
c                  nyf    =   6   11   21    41    81    161    321
c                  nm     =   1   36  157   598  2279   8840  34761
c
c               levels,nxc,nyc,nxf,nyf,nm remain unchanged on exit.
c ---
c     istart    integer.
c               on entry,
c               =2 if the user provides an initial estimate
c                  of the solution in v and the residual of v in vb,
c               =1 if the user provides an initial estimate
c                  of the solution in v,
c               =0 if no initial estimate is provided.
c               unchanged on exit.
c ---
c     iprep     integer.
c               on entry,
c               =1 if  mgd9v has been called before and
c                  neither a nor ldu have been overwritten.
c                 ( this is a useful option if the user desires to solve
c                  a sequence of problems with the same matrices but
c                  with different righthandsides, it prevents a new
c                  set up of coarsegridmatrices and decompositions which
c                  are computed before in a previous call of mgd9v )
c               =0 if mgd9v has not been called before or
c                  either a, ac, ldu, lduc, wa, wac, wb, wbc
c                  have been overwritten.
c               unchanged on exit.
c ---
c     maxit     integer.
c               on entry,
c               maximum number of multigrid iterations
c               a value between 5 and 20 might be reasonable.
c               unchanged on exit.
c ---
c     tol       double precision.
c               on entry,
c               this is the tolerance desired by the user, tol is a
c               bound for the l2-norm of the residual.
c               clearly tol should not be less than
c               nxf * nyf * machineaccuracy * average size of entries
c               in the matrix a.
c               remark  if either maxit iterations or the tolerance
c               ------  have been achieved multigrid cycling is stopped.
c               unchanged on exit.
c ---
c     iout      integer array of dimension 6.
c               on entry,
c               it governs the amount of output desired by the user.
c               smaller iout-values mean less output,
c               possible values are ,
c               iout(1)=1 confirmation of input data
c                       0 none
c               iout(2)=2 matrices and right-hand sides on all levels
c                       1 matrix and right-hand side on highest level
c                       0 none
c               iout(3)=2 matrix-decompositions on all levels
c                       1 matrix-decomposition on highest level
c                       0 none
c               iout(4)=4 norms of residuals, reduction factors,
c                         final residual, final solution
c                       3 norms of residuals, reduction factors,
c                         final residual
c                       2 norms of residuals, reduction factors
c                       1 final norm of residual, number of iterations.
c                       0 none
c               iout(5)=1 the time spent in various subroutines
c                       0 none
c                         remark  clock routines are not standard
c                         ------  fortran. to obtain timings the user
c                                 must adapt the subroutine timing,
c                         it should deliver the cpu-time elapsed.
c               iout(6)=2 prolongation-weights in wb() and wa()
c                         useful only for inquisitive user with detailed
c                         knowledge of the program
c                       1 prolongation-weights in wa()
c                         useful only for inquisitive user with detailed
c                         knowledge of the program
c                       0 none
c               unchanged on exit.
c ---
c     a         double precision array,
c               dimensioned as a(nxf*nyf*9)
c               before entry,
c               the user has to initialize a(1),...,a(nxf*nyf*9) as
c               follows, first declare
c                        double precision a(nxf, nyf, 9)
c               a is the matrix corresponding to the finest grid.
c               the ordering of the points in the grid is as follows
c               the subscript ( i, j ) corresponds to the point
c
c               (x,y) = ( i*h , j*h )
c                            x     y
c                                     i=1,...,nxf  j=1,...,nyf
c               the 9-point difference molecule at the point with
c               subscript ( i, j ) is positioned in the x,y-plane
c               as follows
c
c
c                       y,j
c                        +
c                        +
c                        +   a(i,j,7)   a(i,j,8)   a(i,j,9)
c                        +      .        .
c                        +   a(i,j,4)   a(i,j,5)   a(i,j,6)
c                        +      .        .        .
c                        +   a(i,j,1)   a(i,j,2)   a(i,j,3)
c                        +      .        .        .
c                        +
c                        o+ + + + + + + + + + + + + + + x, i
c
c     important the user has to provide the matrix a only on the finest
c     --------- grid.
c     important the user has to take care that parts of the molecules
c     --------- outside the domain are initialized to zero, otherwise
c                                                     ----
c               wrong results are produced.
c        remark the matrix a remains unchanged on exit.
c        ------
c ---
c     ac        double precision array,
c               dimensioned as ac(nm*9)
c               on exit,
c               ac contains all matrices on the coarser grids
c ---
c     ldu       double precision array,
c               dimensioned as ldu(nxf*nyf*3)
c               on exit,                                              -
c               ldu contains decompositions of all tridiagonal blocks d
c                                                                      j
c ---
c     lduc      double precision array,
c               dimensioned as lduc(nm*3)
c               on exit,                                              -
c               ldu contains decompositions of all tridiagonal blocks d
c                                                                      j
c               on the coarser grids
c ---
c     wa        double precision array,
c               dimensioned as wa(nxf*nyf).
c               on exit,
c               wa contains (part of) the prolongation-weights.
c ---
c     wac       double precision array,
c               dimensioned as wac(nm).
c               on exit,
c               wa contains (part of) the prolongation-weights
c               on the coarser grids
c ---
c     wb        double precision array,
c               dimensioned as wb(nxf*nyf).
c               on exit,
c               wb contains (part of) the prolongation-weights.
c ---
c     wbc       double precision array,
c               dimensioned as wbc(nm).
c               on exit,
c               wa contains (part of) the prolongation-weights
c               on the coarser grids
c ---
c     rhs       double precision array,
c               dimensioned as rhs(nxf*nyf).
c               before entry, the user has to initialize
c               rhs(1),...,rhs(nxf*nyf) with
c               the right-hand side of the equation as follows,
c               first declare
c                             double precision rhs(nxf,nyf)
c               then rhs(i,j) corresponds to the righthandside of the
c               equation at the ( i, j )-th gridpoint.
c     important the user has to provide the right-hand side of the
c     --------- discretized equation only on the finest grid.
c        remark rhs remains unchanged on exit.
c        ------
c ---
c     rhsc      double precision array,
c               dimensioned as rhsc(nm)
c               used as workspace
c ---
c     v         double precision array,
c               dimensioned as v(nxf*nyf).
c               before entry,
c               if istart.ge.1 then v(1),..,v(nxf*nyf) should contain an
c               initial estimate of the solution provided by the user.
c               if istart=0 then no initializing is needed because then
c               v is initialized to zero within mgd9v.
c               after a call of mgd9v, v(1),..,v(nxf*nyf) contains
c               the (approximate) numerical solution.
c ---
c     vc        double precision array,
c               dimensioned as vc(nxf*nyf).
c               used as workspace.
c ---
c     vb        double precision array,
c               dimensioned as vb(nxf*nyf).
c               before entry,
c               if istart=2 then vb(1),...,vb(nxf*nyf) should contain
c               the residual ( rhs - a * v )
c               if istart.le.1 then no initializing is required.
c               after a call of mgd9v, vb contains the residual of the
c               numerical solution v.
c ---
c     vbc       double precision array,
c               dimensioned as vbc(nm).
c               used as workspace.
c ---
c     work      double precision array,
c               dimensioned as work(nxf*12).
c               it is used as a (small) scratch array.
c ---
c     resno     double precision.
c               on exit,
c               this variable contains the l2-norm of the residual at
c               the end of execution of mgd9v.
c ---
c     ndid      integer.
c               on exit,
c               this variable contains the number of multigrid
c               iterations performed
c ---
c     ifail     integer.
c               on entry,
c               the value assigned to ifail must have the
c               decimal expansion  cba  where each of the decimal digits
c               c, b and a has the value 0 or 1.
c                         a = 0 specifies hard failure,
c                               ( i.e. execution of the program will
c                                 terminate if mgd9v detects an error)
c                         a = 1 specifies soft failure,
c                               ( i.e. if mgd9v detects an error, ifail
c                               is reset to the associated error number
c                               and control returns to the calling
c                               program )
c                         b = 0 suppresses error messages,
c                         b = 1 specifies that error messages are to be
c                               output,
c                         c = 0 suppresses advisory messages,
c                         c = 1 specifies that advisory messages are to
c                               be output.
c               ifail = 110 is the normal recommended value for inexper-
c               ienced users, it gives a hard failure with error and
c               advisory output.
c               unless the routine detects an error, ifail contains 0 on
c               exit.
c               errors detected by the routine
c                ( on exit ) ifail = 1
c                                     on entry levels is out of bounds.
c                            ifail = 2
c                                     on entry nxc or nyc is less than 3
c                            ifail = 3
c                                     on entry no consistency among
c                                     levels nxc nxf.
c                            ifail = 4
c                                     on entry no consistency among
c                                     levels nyc nyf.
c                            ifail = 5
c                                     on entry nm is wrong.
c                            ifail = 6
c                                     divergence.
c                            ifail = 7
c                                     poor convergence,but no divergence
c                                     maybe due to user discretization
c                            ifail = 8
c                                     maxit iterations performed without
c                                     reaching tol.
c
c-----------------------------------------------------------------------
      integer i, ifailb, j, l, lev, nadv, nerr, nmcomp
      integer p01cwi
      external p01cwi
      double precision cpb, cpe
      character*8 srname
      intrinsic mod
      integer ngp, ngx, ngy, maxlev
      common /poi/ ngp(12),ngx(12),ngy(12),maxlev
      double precision cp
      common /cpu/ cp(5)
      data srname/'mgd9v   '/
c
      call timing(cpb)
      call x04aaf(0,nerr)
      call x04abf(0,nadv)
c
c     initializing of cp times
c
      do 10 j=1,5
        cp(j)= 0.0d0
   10 continue
      if(iout(1).gt.0) then
         write(nadv,11) levels,nxc,nyc,nxf,nyf,nm,maxit,tol,iout,
     +               istart,iprep,ifail
   11    format(/' multigrid program mgd9v, version june 6 1994'//
     +           ' p.m. de zeeuw '/
     +           ' the centre for mathematics and computerscience   '/
     +           ' department of numerical mathematics              '/
     +           ' kruislaan 413, 1098 sj amsterdam, the netherlands'//
     +      ' levels    nxc    nyc     nxf     nyf        nm'/
     +       3x,i4,1x,i6,1x,i6,1x,i7,1x,i7,1x,i9/
     +      '  maxit                 tol'/i7,7x,d13.6/
     +      18x,'iout  istart  iprep  ifail'/4x,6i3,5x,i3,4x,i3,3x,i4)
      end if
c
c     ngp, ngx, ngy are pointer-arrays that determine how
c     the matrix and vectors on the various grids are stored
c     computation of pointer arrays
c     lev=1, coarsest grid, lev=levels, finest grid
c     ngx(lev), ngy(lev), number of vertical, horizontal
c     grid-lines on grid(lev)
c     (ngp(lev)+1) points to first grid-point on grid(lev)
c     points are counted along horizontal lines in direction
c     of increasing x and y
c
c     ifailc=ifail/100
      ifailb=mod(ifail/10,10)
c     ifaila=ifail-100*ifailc-10*ifailb
c
      ndid = 0
c
      if(levels.ge.1.and.levels.le.12) goto 15
        if(ifailb.gt.0) write(nerr,14)
   14   format(/' levels out of bounds')
        ifail=p01cwi(ifail,1,srname)
        return
   15 if(nxc.ge.3.and.nyc.ge.3) goto 19
        if(ifailb.gt.0) write(nerr,18)
   18   format(/' nxc or nyc too small')
        ifail=p01cwi(ifail,2,srname)
        return
   19 ngx(1)=nxc
      ngy(1)=nyc
      do 20 lev=2,levels
        ngx(lev)= 2*ngx(lev-1)-1
        ngy(lev)= 2*ngy(lev-1)-1
   20 continue
      maxlev=levels
c     verification of nxc nxf nyc nyf
      if (nxf.eq.ngx(levels)) goto 22
         if (ifailb.gt.0) write(nerr,21)
   21    format(/' no consistency among  levels nxc nxf')
         ifail=p01cwi(ifail,3,srname)
         return
   22 if (nyf.eq.ngy(levels)) goto 26
         if (ifailb.gt.0) write(nerr,23)
   23    format(/' no consistency among  levels nyc nyf')
         ifail=p01cwi(ifail,4,srname)
         return
c     computation of pointer array ngp()
   26 ngp(levels)=0
      if (levels.gt.1) then
        ngp(levels-1)=0
      end if
      do 30 l = (levels-2), 1, -1
        ngp(l)=ngp(l+1)+ngx(l+1)*ngy(l+1)
   30 continue
c     verification of nm
      if (levels.gt.1) then
        nmcomp=ngp(1)+nxc*nyc
      else
        nmcomp=1
      end if
      if(nm.eq.nmcomp) goto 32
        if(ifailb.gt.0) write(nerr,31) nmcomp
   31   format(/' nm is wrong should be ',i9)
        ifail=p01cwi(ifail,5,srname)
        return
c     input data are consistent
   32 call prepar(levels, nxf, nyf, nm, iout, istart, iprep,
     +            v, rhs, a, ac, vb,
     +            ldu, lduc, wa, wac, wb, wbc, resno, work)
      call cycles(levels, nxf, nyf, nm, iout, maxit, tol,
     +            rhs, rhsc, a, ac, ldu, lduc, wa, wac, wb, wbc,
     +            v, vc, vb, vbc, resno, work, ndid, ifail)
      call timing(cpe)
      cp(5)=(cpe-cpb)
      cp(2)=cp(2)-cp(1)
      if(iout(5).gt.0) write(nadv,47) (cp(i),i=1,5)
   47 format(//
     +       ' weights  ',f12.3//
     +       ' galerkin ',f12.3//
     +       ' decompose',f12.3//
     +       ' mg-cycles',f12.3//
     +       ' total    ',f12.3)
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine prepar(levels, nxf, nyf, nm, iout, istart, iprep,
     +            v, rhs, a, ac, vb,
     +            ldu, lduc, wa, wac, wb, wbc, resno, work)
      integer levels, nxf, nyf, nm, iout(6), istart, iprep
      double precision
     +     v(nxf*nyf), rhs(nxf*nyf),
     +     a(nxf*nyf*9), ac(nm*9),
     +     vb(nxf*nyf), ldu(nxf*nyf*3), lduc(nm*3),
     +     wa(nxf*nyf), wac(nm), wb(nxf*nyf), wbc(nm),
     +     resno, work(nxf*12)
c-----------------------------------------------------------------------
c     this is the preparational phase, coarse grid matrices and
c     decompositions on all levels are computed.
c     if istart=0 the first approximation of the solution is initialized
c                 to zero in array v.
c     if istart=1 prepar expects a user estimate of the solution in v.
c     if istart=2 as =1 and moreover the residual (rhs - a * v) in vb.
c     after completion vb contains the residual.
c     resno becomes the l2-norm of the residual.
c
c     by  p.m. de zeeuw,  version 930913
c-----------------------------------------------------------------------
      integer lev, nadv, nf, np1, nx, ny
      double precision vl2nor
      external vl2nor
      integer ngp, ngx, ngy, maxlev
      common /poi/ ngp(12),ngx(12),ngy(12),maxlev
c
      call x04abf(0,nadv)
      if ((iout(2).gt.0).and.(iprep.lt.1)) then
         write(nadv,23) levels
   23    format(/' matrix by user, level = ', i3)
         call outmat(nadv,a, nxf, nyf, 9)
      end if
      if (iout(2).gt.0) then
         write(nadv,24)
   24    format(/' right-hand side (rhs) by user')
         call outvec(nadv,rhs,nxf,nyf)
      end if
      if( iprep.ge.1 ) goto 55
c     computation of coarse grid matrices
c
      call rap( levels, nm, nxf, nyf,
     +          wa, wac, wb, wbc, a, ac)
c
      if (iout(6).ge.1) then
          write(nadv,26) levels
          call outvec( nadv, wa, nxf, nyf)
          write(nadv,27) levels
          call outvec( nadv, wb, nxf, nyf)
         do 30 lev = levels-1, 2,-1
          np1=ngp(lev)+1
          nx=ngx(lev)
          ny=ngy(lev)
          write(nadv,26) lev
   26     format('1'/' weights wa, level = ', i3)
          call outvec( nadv, wac(np1), nx, ny)
          if( iout(6).ge.2) then
             write(nadv,27) lev
   27        format('1'/' weights wb, level = ', i3)
             call outvec( nadv, wbc(np1), nx, ny)
          end if
   30    continue
      end if
      if (iout(2).gt.1) then
         do 40 lev = (levels-1), 1,-1
          write(nadv,36) lev
   36     format(/' matrix by galerkin (rap), level = ', i3)
          call outmat(nadv, ac(ngp(lev)*9+1), ngx(lev), ngy(lev), 9)
   40    continue
      end if
c     computation of line l u decompositions
c
      call decomp(levels, nm, nxf, nyf, a, ac, ldu, lduc,
     +            (nxf*(7+5)), work)
c
      if (iout(3).gt.0) then
         write(nadv,43) levels
   43    format(/' line l u decomposition, level = ',i3)
         call outmat(nadv,ldu, nxf, nyf, 3)
      end if
      if (iout(3).gt.1) then
         do 50 lev = (levels-1), 1,-1
          write(nadv,47) lev
   47     format(/' line l u decomposition, level = ',i3)
          call outmat(nadv, lduc(ngp(lev)*3+1), ngx(lev), ngy(lev), 3)
   50    continue
      end if
   55 continue
      nf=nxf*nyf
      if(istart.le.0) then
          call zeros(nf, v)
          call copy(nf,rhs,vb)
      else if(istart.eq.1) then
          call residu(nxf,nf,rhs,a,v,vb)
      end if
      resno=vl2nor(vb,nf)
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine cycles(levels, nxf, nyf, nm, iout, maxit,tol,
     +            rhs, rhsc, a, ac, ldu, lduc, wa, wac, wb, wbc,
     +            v, vc, vb, vbc,
     +            resno, work, ndid, ifail)
      integer levels, nxf, nyf, nm, iout(6), maxit, ndid, ifail
      double precision
     +     tol, rhs(nxf*nyf), rhsc(nm), a(nxf*nyf*9), ac(nm*9),
     +     ldu(nxf*nyf*3), lduc(nm*3),
     +     wa(nxf*nyf), wac(nm), wb(nxf*nyf), wbc(nm),
     +     v(nxf*nyf), vc(nm), vb(nxf*nyf), vbc(nm), resno, work(nxf*2)
c-----------------------------------------------------------------------
c     governs the process of sawtooth multigrid cycling.
c
c     cycles expects the 9-point difference matrix in a and the right-
c     hand side in rhs.
c     further cycles expects an approximate solution in v, the residual
c     of v in vb and the l2-norm of vb in resno.
c
c     after a call of mgd9v or a previous call of cycles, cycles can be
c     used to refine the solution obtained.
c     the parameters have the same meaning as in mgd9v.
c
c     by  p.m. de zeeuw,  version 930913
c-----------------------------------------------------------------------
      integer ifailb, ifailc, k, lev, mgit, nadv, nerr, nf,
     +     npl, nplm, npl1, nxl, nxlm, nxl1, nyl, nylm, nyl1
      double precision
     +     cpb, cpe, resnav, resn1, resn2, resn3, resn32, eps
      parameter( eps = 1.0d-19 )
      integer p01cwi
      external p01cwi
      double precision vl2nor
      external vl2nor
      logical pr
      character*8 srname
      intrinsic exp, log, mod
      integer ngp, ngx, ngy, maxlev
      common /poi/ ngp(12),ngx(12),ngy(12),maxlev
      double precision cp
      common /cpu/ cp(5)
      data srname/'mgd9v   '/
c
      call timing(cpb)
      call x04aaf(0,nerr)
      call x04abf(0,nadv)
c
      ndid = 0
c
      ifailc=ifail/100
      ifailb=mod(ifail/10,10)
c     ifaila=ifail-100*ifailc-10*ifailb
c
      if(maxit.le.0) return
      resn1=resno
      resn2=resno
      if(iout(4).ge.1) write(nadv,9) resn1
    9 format(/' l2-norm of initial residual= ',d11.3)
      nf=nxf*nyf
c     start of iterations
c
      do 100 mgit=1, maxit
c
        if(levels.le.1) goto 42
        npl=ngp(levels-1)+1
        nxl=ngx(levels-1)
        nyl=ngy(levels-1)
        call restri( nxf, nyf, wa, wb, vb,
     +               nxl, nyl, rhsc)
        do 10 lev = (levels-2), 1,-1
          npl1=ngp(lev+1)+1
          nxl1=ngx(lev+1)
          nyl1=ngy(lev+1)
          npl=ngp(lev)+1
          nxl=ngx(lev)
          nyl=ngy(lev)
          call restri(nxl1, nyl1, wac(npl1), wbc(npl1), rhsc(npl1),
     +                nxl , nyl , rhsc(npl) )
   10   continue
        if((iout(2).ge.2).and.(mgit.eq.1)) then
           do 20 lev = (levels-1), 1,-1
             write(nadv,17)  lev
   17        format(/' restriction of residual, level = ',i3)
             call outvec(nadv,rhsc(ngp(lev)+1),ngx(lev),ngy(lev))
   20      continue
        end if
        if( levels.eq.1 ) then
            call smooth(1,v,vb,a,rhs,ldu,work,.false.)
        else
            call smooth(1,vc,vbc,ac,rhsc,lduc,work,.false.)
        end if
c       luctor et emergo
        do 30 lev = 2, (levels-1)
          nplm=ngp(lev-1)+1
          nxlm=ngx(lev-1)
          nylm=ngy(lev-1)
          npl=ngp(lev)+1
          nxl=ngx(lev)
          nyl=ngy(lev)
          call prolon(nxlm, nylm, vc(nplm),
     +                nxl , nyl , wac(npl), wbc(npl), vc(npl))
          call smooth(lev,vc,vbc,ac,rhsc,lduc,work,.false.)
   30   continue
        nplm=ngp(levels-1)+1
        nxlm=ngx(levels-1)
        nylm=ngy(levels-1)
        call prolon(nxlm, nylm,  vc(nplm),
     +              nxf , nyf , wa, wb, vb)
        do 37 k = 1, nf
          v(k)=v(k)+vb(k)
   37   continue
   42   continue
        if( levels.eq.1 ) then
           call smooth(levels,v,vb,a,rhs,ldu,work,.true. )
        else
           call smooth(levels,v,vb,a,rhs,ldu,work,.false.)
        end if
        call residu(nxf,nf,rhs,a,v,vb)
        ndid = ndid + 1
        resn3=vl2nor(vb,nf)
c
        if( resn3.ge.resn2 ) then
          resn32=(resn3+eps)/(resn2+eps)
        else
          resn32= resn3/resn2
        end if
c
        resn2=resn3
        if( resn3.ge.resn1 ) then
          resnav=exp(log((resn3+eps)/(resn1+eps))/mgit)
        else
c         resnav=exp(log(resn3/resn1)/mgit)
c Dit was nog niet waterdicht, denk aan log van 0
          resnav=resn3/resn1
          if( resnav .lt. eps ) then
            resnav = eps
          end if
          resnav=exp(log(resnav)/mgit)
        end if
c
        pr= (iout(4).ge.2)
        if(.not.pr) then
           pr= ((iout(4).ge.1).and.((mgit.eq.maxit).or.(resn3.lt.tol)))
        end if
        if(pr) write(nadv,45) mgit,resn3,resn32,resnav
   45   format(/' mgd9v, iteration number =  ',i6/
     +          ' l2-norm of residual     = ',d11.3/
     +          ' current reduction factor= ',d11.3/
     +          ' average reduction factor= ',d11.3)
        if(resn3.le.tol) goto 111
c
  100 continue
c
  111 resno=resn3
      if (iout(4).gt.2) then
         write(nadv,333)
  333    format(/' residual')
         call outvec(nadv,vb,nxf,nyf)
      end if
      if (iout(4).gt.3) then
         write(nadv,444)
  444    format(/' numerical solution')
         call outvec(nadv,v ,nxf,nyf)
      end if
      call timing(cpe)
      cp(4)=cp(4)+(cpe-cpb)
      if( (resnav.ge.1.0d0).and.(resno.gt.tol) ) then
         if(ifailb.gt.0) write(nerr,555)
  555    format(/' divergence.')
         ifail=p01cwi(ifail,6,srname)
         return
      end if
      if( (resnav.ge.0.85d0).and.(resno.gt.tol) ) then
         if(ifailb.gt.0) write(nerr,666)
  666    format(/' poor convergence.')
         ifail=p01cwi(ifail,7,srname)
         return
      end if
      if(resn3.gt.tol) then
         if(ifailb.gt.0) write(nerr,888)
  888    format(/' maxit iterations performed without reaching tol.')
         if((ifailc.gt.0).and.(resnav.lt.0.55d0)) write(nadv,999)
  999    format(/' good convergencerate, so v and vb are valuable,'
     +          /' at this point restart mgd9v with larger maxit and',
     +           ' with istart=2 and iprep=1')
         ifail=p01cwi(ifail,8,srname)
         return
      end if
      ifail=0
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine rap(levels,nm,nxf,nyf,
     +               wa,wac,wb,wbc,a,ac)
      integer levels, nm, nxf, nyf
      double precision
     +     wa(nxf*nyf), wac(nm), wb(nxf*nyf), wbc(nm),
     +     a(nxf*nyf*9), ac(nm*9)
c-----------------------------------------------------------------------
c     constructs coarse grid matrices on all lower levels by means of
c     galerkin approximation.
c-----------------------------------------------------------------------
      integer levc, levf, npac, npaf, nw, nxc, nx, nyc, ny
      double precision cpb, cpe
      integer ngp, ngx, ngy, maxlev
      common /poi/ ngp(12),ngx(12),ngy(12),maxlev
      double precision cp
      common /cpu/ cp(5)
      call timing(cpb)
      if(levels.le.1) return
      nx=ngx(levels)
      ny=ngy(levels)
      nyc=ngy(levels-1)
      nxc=ngx(levels-1)
      call zeros( 9*nm, ac)
      call weight( nx, ny, a, wa, wb)
      call glrka9( nx, ny, a, wa, wb,
     +               nxc, nyc, ac)
      do 100 levc = (levels-2), 1,-1
        levf=levc+1
        nx=ngx(levf)
        ny=ngy(levf)
        npaf=9*ngp(levf)+1
        nw=ngp(levf)+1
        nyc=ngy(levc)
        nxc=ngx(levc)
        npac=9*ngp(levc)+1
        call weight( nx, ny, ac(npaf), wac(nw), wbc(nw))
        call glrka9( nx, ny, ac(npaf), wac(nw), wbc(nw),
     +               nxc, nyc, ac(npac))
  100 continue
      call timing(cpe)
      cp(2)=cp(2)+(cpe-cpb)
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine weight(nx,ny,a,wa,wb)
      integer nx, ny
      double precision
     +     a(nx,ny,9), wa(nx,ny), wb(nx,ny)
      double precision cpb, cpe
      double precision cp
      common /cpu/ cp(5)
      call timing(cpb)
      call whoriz(nx,ny,a,wa,wb)
      call wverti(nx,ny,a,wa,wb)
      call wdiag(nx,ny,a,wa,wb)
      call timing(cpe)
      cp(1)=cp(1)+(cpe-cpb)
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine whoriz(nx,ny,a,wa,wb)
      integer nx, ny
      double precision
     +     a(nx,ny,9), wa(nx,ny), wb(nx,ny)
c-----------------------------------------------------------------------
c     computes weights on horizontal lines.
c-----------------------------------------------------------------------
      double precision zero, one, half, eps, twoeps
      parameter( zero = 0.0d0  , one    = 1.0d0, half = 0.5d0,
     +           eps  = 1.0d-14, twoeps = 2.0d-14 )
      integer i, j
      double precision
     +     as1, as3, as7, as9, a1, a3, a4, a6, a7, a9, c1,
     +     de, dn, ds, dw, sig, s1, s2, s3, s4, s6, s7, s8, s9, ww
      logical notlow, notupp
      intrinsic abs, max
c
      do 20 j = 1,  ny   , 2
        notlow=(j.gt.1 )
        notupp=(j.lt.ny)
        do 10 i = 2, (nx-1), 2
          if(notlow) then
                     s1=a(i,j,1)+a(i-1,j-1,9)
                     a1=a(i,j,1)-a(i-1,j-1,9)
                     s2=a(i,j,2)+a(i  ,j-1,8)
                     s3=a(i,j,3)+a(i+1,j-1,7)
                     a3=a(i,j,3)-a(i+1,j-1,7)
          else
                     s1=zero
                     a1=zero
                     s2=zero
                     s3=zero
                     a3=zero
          end if
          as1=abs(s1)
          as3=abs(s3)
           s4=a(i,j,4)+a(i-1,j  ,6)
           a4=a(i,j,4)-a(i-1,j  ,6)
           s6=a(i,j,6)+a(i+1,j,4)
           a6=a(i,j,6)-a(i+1,j,4)
          if(notupp) then
                     s7=a(i,j,7)+a(i-1,j+1,3)
                     a7=a(i,j,7)-a(i-1,j+1,3)
                     s8=a(i,j,8)+a(i  ,j+1,2)
                     s9=a(i,j,9)+a(i+1,j+1,1)
                     a9=a(i,j,9)-a(i+1,j+1,1)
          else
                     s7=zero
                     a7=zero
                     s8=zero
                     s9=zero
                     a9=zero
          end if
          as7=abs(s7)
          as9=abs(s9)
c
          dw=max( abs(s1+s4+s7), as1, as7)
c
          de=max( abs(s3+s6+s9), as3, as9)
c
          dn=max( abs(s7+s8+s9), as7, as9)
c
          ds=max( abs(s1+s2+s3), as1, as3)
c
          c1=(a3+a6+a9)-(a1+a4+a7)
c
          sig=abs( (a(i,j,7)+a(i,j,8)+a(i,j,9)+
     +              a(i,j,4)         +a(i,j,6)+
     +              a(i,j,1)+a(i,j,2)+a(i,j,3))/a(i,j,5) )
          if(sig.gt.one) sig=one
          ww= sig*half*( one+(dw-de)/(dw+de+eps)+
     +                            c1/(dw+de+dn+ds+twoeps) )
          if(ww.lt.zero) ww=zero
          if(ww.gt.sig ) ww=sig
          wb(i,j)= ww
          wa(i,j)= sig-ww
   10   continue
   20 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine wverti(nx,ny,a,wa,wb)
      integer nx, ny
      double precision
     +     a(nx,ny,9), wa(nx,ny), wb(nx,ny)
c-----------------------------------------------------------------------
c     computes weights on vertical lines.
c-----------------------------------------------------------------------
      double precision zero, one, half, eps, twoeps
      parameter( zero = 0.0d0  , one    = 1.0d0, half = 0.5d0,
     +           eps  = 1.0d-14, twoeps = 2.0d-14 )
      integer i, j
      double precision
     +     as1, as3, as7, as9, a1, a2, a3, a7, a8, a9, c2,
     +     de, dn, ds, dw, sig, s1, s2, s3, s4, s6, s7, s8, s9, ws
      logical ntleft, ntrght
      intrinsic abs, max
c
      do 20 j = 2, (ny-1), 2
        do 10 i = 1,  nx   , 2
          ntleft=(i.gt.1 )
          ntrght=(i.lt.nx)
c
          if(ntleft) then
                     s1=a(i,j,1)+a(i-1,j-1,9)
                     a1=a(i,j,1)-a(i-1,j-1,9)
                     s4=a(i,j,4)+a(i-1,j  ,6)
                     s7=a(i,j,7)+a(i-1,j+1,3)
                     a7=a(i,j,7)-a(i-1,j+1,3)
          else
                     s1=zero
                     a1=zero
                     s4=zero
                     s7=zero
                     a7=zero
          end if
          as1=abs(s1)
          as7=abs(s7)
c
          if(ntrght) then
                     s3=a(i,j,3)+a(i+1,j-1,7)
                     a3=a(i,j,3)-a(i+1,j-1,7)
                     s6=a(i,j,6)+a(i+1,j  ,4)
                     s9=a(i,j,9)+a(i+1,j+1,1)
                     a9=a(i,j,9)-a(i+1,j+1,1)
          else
                     s3=zero
                     a3=zero
                     s6=zero
                     s9=zero
                     a9=zero
          end if
          as3=abs(s3)
          as9=abs(s9)
c
          s2=a(i,j,2)+a(i  ,j-1,8)
          a2=a(i,j,2)-a(i  ,j-1,8)
          s8=a(i,j,8)+a(i  ,j+1,2)
          a8=a(i,j,8)-a(i  ,j+1,2)
c
          dw=max( abs(s1+s4+s7), as1, as7)
c
          de=max( abs(s3+s6+s9), as3, as9)
c
          dn=max( abs(s7+s8+s9), as7, as9)
c
          ds=max( abs(s1+s2+s3), as1, as3)
c
          c2=(a7+a8+a9)-(a1+a2+a3)
c
          sig=abs( (a(i,j,7)+a(i,j,8)+a(i,j,9)+
     +              a(i,j,4)         +a(i,j,6)+
     +              a(i,j,1)+a(i,j,2)+a(i,j,3))/a(i,j,5) )
          if(sig.gt.one) sig=one
          ws= sig*half*(one+(ds-dn)/(ds+dn+eps)+
     +                           c2/(dw+de+dn+ds+twoeps) )
          if(ws.lt.zero) ws=zero
          if(ws.gt.sig ) ws=sig
          wa(i,j)= ws
          wb(i,j)= sig-ws
   10   continue
   20 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine wdiag(nx,ny,a,wa,wb)
      integer nx, ny
      double precision
     +     a(nx,ny,9), wa(nx,ny), wb(nx,ny)
c-----------------------------------------------------------------------
c     wdiag assumes that weights on horizontal and vertical lines are
c     already known.
c-----------------------------------------------------------------------
      double precision one
      parameter( one = 1.0d0 )
      integer i, j
      double precision oa5
c
      do 20 j = 2, (ny-1), 2
        do 10 i = 2, (nx-1), 2
          oa5= -one/a(i,j,5)
          wa(i,j)= (a(i,j,2)*wa(i,j-1)+a(i,j,3)+a(i,j,6)*wa(i+1,j)
     +                     )*oa5
          wb(i,j)= (a(i,j,4)*wb(i-1,j)+a(i,j,7)+a(i,j,8)*wb(i,j+1)
     +                     )*oa5
c         wd(i,j)
          wb(i-1,j-1)= (a(i,j,1)+a(i,j,2)*wb(i,j-1)+a(i,j,4)*wa(i-1,j)
     +                     )*oa5
c         wc(i,j)
          wa(i-1,j-1)= (a(i,j,6)*wb(i+1,j)+a(i,j,8)*wa(i,j+1)+a(i,j,9)
     +                     )*oa5
   10   continue
   20 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine glrka9( nxf, nyf, f, wa, wb, nxc, nyc, c)
      integer nxf, nyf, nxc, nyc
      double precision
     +     f(nxf,nyf,9), wa(nxf,nyf), wb(nxf,nyf), c(nxc,nyc,9)
c-----------------------------------------------------------------------
c     derives coarse grid matrix from fine grid matrix in a
c     finite-element-method-like way
c     by p.m. de zeeuw
c-----------------------------------------------------------------------
      integer ic, jc, jf, jff, jfp, nxcm, nycm
      double precision
     +     a01, a10, a11, a12, a21, b01, b10, b11, b12, b21, c11, d11
c
c     ora et labora
      jff=1
      nxcm=nxc-1
      nycm=nyc-1
      do 200 jc=1,nycm
c  ---------------------------------------------------------------------
c     left hand boundary at ( ic+1=1, jc )
c  ---------------------------------------------------------------------
      jf=jff
      jfp=jf+1
      jff=jf+2
c     iff=1
      a21=wa(   1,jfp )
      b21=wb(   1,jfp )
c
c
c      (2,0) to (2,0)
      c(   1,jc  , 5)= c(   1,jc  , 5)                    +
     +                             f(  1,jf ,5)           +
     +                   a21 * (   f(  1,jf ,8)         +
     +                             f(  1,jfp,2)         +
     +                             f(  1,jfp,5) * a21   )
c
c      (2,0) to (2,2)
      c(   1,jc  , 8)= c(   1,jc  , 8)                    +
     +                             f(  1,jf ,8) * b21     +
     +                   a21 * (   f(  1,jfp,5) * b21   +
     +                             f(  1,jfp,8)         )
c
c      (2,2) to (2,0)
      c(   1,jc+1, 2)= c(   1,jc+1, 2)                    +
     +                   b21 * (   f(  1,jfp,2)         +
     +                             f(  1,jfp,5) * a21   ) +
     +                             f(  1,jff,2) * a21
c
c      (2,2) to (2,2)
      c(   1,jc+1, 5)= c(   1,jc+1, 5)                    +
     +                   b21 * (   f(  1,jfp,5) * b21   +
     +                             f(  1,jfp,8)         +
     +                             f(  1,jff,2)         )
c  ---------------------------------------------------------------------
c     inner area
c  ---------------------------------------------------------------------
      do  10 ic=1,nxcm
      b10=wb(2*ic  ,jf  )
      a01=wa(2*ic-1,jfp )
      d11=wb(2*ic-1,jf  )
c
c     (0,0) to (0,0)
      c(ic  ,jc  , 5)= c(ic  ,jc  , 5)                       +
     +                   b10 * (   f(2*ic-1,jf ,6)         +
     +                             f(2*ic  ,jf ,4)         +
     +                             f(2*ic  ,jf ,5) * b10   +
     +                             f(2*ic  ,jf ,8) * d11   +
     +                           ( f(2*ic  ,jf ,7) +
     +                             f(2*ic-1,jfp,3) ) * a01 ) +
     +                   d11 * (   f(2*ic-1,jf ,9)         +
     +                             f(2*ic-1,jfp,6) * a01   )
   10 continue
      do  20 ic=1,nxcm
      a10=wa(2*ic  ,jf  )
      b10=wb(2*ic  ,jf  )
      a01=wa(2*ic-1,jfp )
      a11=wa(2*ic  ,jfp )
      a21=wa(2*ic+1,jfp )
c
c     (0,0) to (2,0)
      c(ic  ,jc  , 6)= c(ic  ,jc  , 6)                       +
     +                   b10 * (   f(2*ic  ,jf ,5) * a10   +
     +                             f(2*ic  ,jf ,6)         +
     +                             f(2*ic  ,jf ,8) * a11   +
     +                             f(2*ic  ,jf ,9) * a21   ) +
     +                             f(2*ic-1,jf ,6) * a10     +
     +                             f(2*ic-1,jf ,9) * a11     +
     +                   a01 * (   f(2*ic-1,jfp,3) * a10   +
     +                             f(2*ic-1,jfp,6) * a11   )
   20 continue
      do  30 ic=1,nxcm
      b10=wb(2*ic  ,jf  )
      a01=wa(2*ic-1,jfp )
      b01=wb(2*ic-1,jfp )
      b11=wb(2*ic  ,jfp )
      b12=wb(2*ic  ,jff )
c
c     (0,0) to (0,2)
      c(ic  ,jc  , 8)= c(ic  ,jc  , 8)                       +
     +                   b11 * (   f(2*ic-1,jf ,9)         +
     +                             f(2*ic-1,jfp,6) * a01   +
     +                             f(2*ic  ,jf ,8) * b10   ) +
     +                       b10 * f(2*ic  ,jf ,7) * b01     +
     +                       a01 * f(2*ic-1,jfp,9) * b12
   30 continue
      do  40 ic=1,nxcm
      b10=wb(2*ic  ,jf  )
      a01=wa(2*ic-1,jfp )
      c11=wa(2*ic-1,jf  )
      b21=wb(2*ic+1,jfp )
      a12=wa(2*ic  ,jff )
c
c     (0,0) to (2,2)
      c(ic  ,jc  , 9)= c(ic  ,jc  , 9)                       +
     +                   c11 * (   f(2*ic-1,jf ,9)         +
     +                             f(2*ic  ,jf ,8) * b10   +
     +                             f(2*ic-1,jfp,6) * a01   ) +
     +                       b10 * f(2*ic  ,jf ,9) * b21     +
     +                       a01 * f(2*ic-1,jfp,9) * a12
   40 continue
      do  50 ic=1,nxcm
      a10=wa(2*ic  ,jf  )
      b10=wb(2*ic  ,jf  )
      a01=wa(2*ic-1,jfp )
      d11=wb(2*ic-1,jf  )
      a21=wa(2*ic+1,jfp )
c
c     (2,0) to (0,0)
      c(ic+1,jc  , 4)= c(ic+1,jc  , 4)                       +
     +                   a10 * (   f(2*ic  ,jf ,4)         +
     +                             f(2*ic  ,jf ,5) * b10   +
     +                             f(2*ic  ,jf ,7) * a01   +
     +                             f(2*ic  ,jf ,8) * d11   ) +
     +                             f(2*ic+1,jf ,4) * b10     +
     +                             f(2*ic+1,jf ,7) * d11     +
     +                   a21 * (   f(2*ic+1,jfp,1) * b10   +
     +                             f(2*ic+1,jfp,4) * d11   )
   50 continue
      do  60 ic=1,nxcm
      a10=wa(2*ic  ,jf  )
      a11=wa(2*ic  ,jfp )
      a21=wa(2*ic+1,jfp )
c
c     (2,0) to (2,0)
      c(ic+1,jc  , 5)= c(ic+1,jc  , 5)                       +
     +                   a10 * (   f(2*ic  ,jf ,5) * a10   +
     +                             f(2*ic  ,jf ,6)         +
     +                             f(2*ic  ,jf ,8) * a11   +
     +                             f(2*ic+1,jf ,4)         +
     +                           ( f(2*ic  ,jf ,9) +
     +                             f(2*ic+1,jfp,1) ) * a21 ) +
     +                             f(2*ic+1,jf ,7) * a11     +
     +                             f(2*ic+1,jf ,5)           +
     +                   a21 * (   f(2*ic+1,jf ,8)         +
     +                             f(2*ic+1,jfp,2)         +
     +                             f(2*ic+1,jfp,4) * a11   +
     +                             f(2*ic+1,jfp,5) * a21   )
   60 continue
      do  70 ic=1,nxcm
      a10=wa(2*ic  ,jf  )
      b01=wb(2*ic-1,jfp )
      b11=wb(2*ic  ,jfp )
      a21=wa(2*ic+1,jfp )
      b12=wb(2*ic  ,jff )
c
c     (2,0) to (0,2)
      c(ic+1,jc  , 7)= c(ic+1,jc  , 7)                       +
     +                   b11 * (   f(2*ic  ,jf ,8) * a10   +
     +                             f(2*ic+1,jf ,7)         +
     +                             f(2*ic+1,jfp,4) * a21   ) +
     +                       a10 * f(2*ic  ,jf ,7) * b01     +
     +                       a21 * f(2*ic+1,jfp,7) * b12
   70 continue
      do  80 ic=1,nxcm
      a10=wa(2*ic  ,jf  )
      c11=wa(2*ic-1,jf  )
      a21=wa(2*ic+1,jfp )
      b21=wb(2*ic+1,jfp )
      a12=wa(2*ic  ,jff )
c
c     (2,0) to (2,2)
      c(ic+1,jc  , 8)= c(ic+1,jc  , 8)                       +
     +                   a10 * (   f(2*ic  ,jf ,8) * c11   +
     +                             f(2*ic  ,jf ,9) * b21   ) +
     +                             f(2*ic+1,jf ,7) * c11     +
     +                             f(2*ic+1,jf ,8) * b21     +
     +                   a21 * (   f(2*ic+1,jfp,4) * c11   +
     +                             f(2*ic+1,jfp,5) * b21   +
     +                             f(2*ic+1,jfp,7) * a12   +
     +                             f(2*ic+1,jfp,8)         )
   80 continue
      do  90 ic=1,nxcm
      b10=wb(2*ic  ,jf  )
      a01=wa(2*ic-1,jfp )
      b01=wb(2*ic-1,jfp )
      d11=wb(2*ic-1,jf  )
      b12=wb(2*ic  ,jff )
c
c     (0,2) to (0,0)
      c(ic  ,jc+1, 2)= c(ic  ,jc+1, 2)                       +
     +                       b01 * f(2*ic-1,jfp,3) * b10     +
     +                   d11 * (   f(2*ic-1,jfp,6) * b01   +
     +                             f(2*ic-1,jff,3)         +
     +                             f(2*ic  ,jff,2) * b12   ) +
     +                       b12 * f(2*ic  ,jff,1) * a01
   90 continue
      do 100 ic=1,nxcm
      a10=wa(2*ic  ,jf  )
      b01=wb(2*ic-1,jfp )
      a11=wa(2*ic  ,jfp )
      a21=wa(2*ic+1,jfp )
      b12=wb(2*ic  ,jff )
c
c     (0,2) to (2,0)
      c(ic  ,jc+1, 3)= c(ic  ,jc+1, 3)                       +
     +                       b01 * f(2*ic-1,jfp,3) * a10     +
     +                   a11 * (   f(2*ic-1,jfp,6) * b01   +
     +                             f(2*ic-1,jff,3)         +
     +                             f(2*ic  ,jff,2) * b12   ) +
     +                       b12 * f(2*ic  ,jff,3) * a21
  100 continue
      do 110 ic=1,nxcm
      b01=wb(2*ic-1,jfp )
      b11=wb(2*ic  ,jfp )
      b12=wb(2*ic  ,jff )
c
c     (0,2) to (0,2)
      c(ic  ,jc+1, 5)= c(ic  ,jc+1, 5)                       +
     +                   b11 * (   f(2*ic-1,jfp,6) * b01   +
     +                             f(2*ic-1,jff,3)         ) +
     +                   b12 * ( ( f(2*ic-1,jfp,9) +
     +                             f(2*ic  ,jff,1) ) * b01 +
     +                             f(2*ic  ,jff,2) * b11   )
  110 continue
      do 120 ic=1,nxcm
      b01=wb(2*ic-1,jfp )
      c11=wa(2*ic-1,jf  )
      b21=wb(2*ic+1,jfp )
      a12=wa(2*ic  ,jff )
      b12=wb(2*ic  ,jff )
c
c     (0,2) to (2,2)
      c(ic  ,jc+1, 6)= c(ic  ,jc+1, 6)                       +
     +                       b01 * f(2*ic-1,jfp,9) * a12     +
     +                   c11 * (   f(2*ic-1,jfp,6) * b01 +
     +                             f(2*ic-1,jff,3)       +
     +                             f(2*ic  ,jff,2) * b12   ) +
     +                       b12 * f(2*ic  ,jff,3) * b21
  120 continue
      do 130 ic=1,nxcm
      b10=wb(2*ic  ,jf  )
      a01=wa(2*ic-1,jfp )
      d11=wb(2*ic-1,jf  )
      b21=wb(2*ic+1,jfp )
      a12=wa(2*ic  ,jff )
c
c     (2,2) to (0,0)
      c(ic+1,jc+1, 1)= c(ic+1,jc+1, 1)                       +
     +                       b21 * f(2*ic+1,jfp,1) * b10     +
     +                       a12 * f(2*ic  ,jff,1) * a01     +
     +                   d11 * (   f(2*ic+1,jfp,4) * b21   +
     +                             f(2*ic  ,jff,2) * a12   +
     +                             f(2*ic+1,jff,1)         )
  130 continue
      do 140 ic=1,nxcm
      a10=wa(2*ic  ,jf  )
      a11=wa(2*ic  ,jfp )
      a21=wa(2*ic+1,jfp )
      b21=wb(2*ic+1,jfp )
      a12=wa(2*ic  ,jff )
c     (2,2) to (2,0)
      c(ic+1,jc+1, 2)= c(ic+1,jc+1, 2)                       +
     +                   b21 * (   f(2*ic+1,jfp,1) * a10   +
     +                             f(2*ic+1,jfp,2)         +
     +                             f(2*ic+1,jfp,4) * a11   +
     +                             f(2*ic+1,jfp,5) * a21   ) +
     +                   a12 * (   f(2*ic  ,jff,2) * a11   +
     +                             f(2*ic  ,jff,3) * a21   ) +
     +                             f(2*ic+1,jff,1) * a11     +
     +                             f(2*ic+1,jff,2) * a21
  140 continue
      do 150 ic=1,nxcm
      b01=wb(2*ic-1,jfp )
      b11=wb(2*ic  ,jfp )
      b21=wb(2*ic+1,jfp )
      a12=wa(2*ic  ,jff )
      b12=wb(2*ic  ,jff )
c
c     (2,2) to (0,2)
      c(ic+1,jc+1, 4)= c(ic+1,jc+1, 4)                       +
     +                       a12 * f(2*ic  ,jff,1) * b01     +
     +                       b21 * f(2*ic+1,jfp,7) * b12     +
     +                    b11 * (  f(2*ic+1,jfp,4) * b21   +
     +                             f(2*ic  ,jff,2) * a12   +
     +                             f(2*ic+1,jff,1)         )
  150 continue
      do 160 ic=1,nxcm
      c11=wa(2*ic-1,jf  )
      b21=wb(2*ic+1,jfp )
      a12=wa(2*ic  ,jff )
c
c     (2,2) to (2,2)
      c(ic+1,jc+1, 5)= c(ic+1,jc+1, 5)                       +
     +                   b21 * (   f(2*ic+1,jfp,4) * c11   +
     +                             f(2*ic+1,jfp,5) * b21   +
     +                             f(2*ic+1,jfp,8)         +
     +                             f(2*ic+1,jff,2)         +
     +                           ( f(2*ic+1,jfp,7) +
     +                             f(2*ic  ,jff,3) ) * a12 ) +
     +                   c11 * (   f(2*ic  ,jff,2) * a12   +
     +                             f(2*ic+1,jff,1)         )
  160 continue
  200 continue
c
c  ---------------------------------------------------------------------
c     upper left hand corner point at ( ic+1=1, jc= nycm+1=nyc )
c  ---------------------------------------------------------------------
      c(   1,nyc , 5)= c(   1,nyc ,5) + f(   1,nyf ,5)
c
c  ---------------------------------------------------------------------
c     upper boundary at ( ic, jc=nyc ) and jf=nyf
c  ---------------------------------------------------------------------
      do 210 ic=1,nxcm
      b10=wb(2*ic  ,nyf )
c
c      (0,0) to (0,0)
      c(ic  ,nyc , 5)= c(ic  ,nyc , 5)                       +
     +                   b10 * (   f(2*ic-1,nyf,6)         +
     +                             f(2*ic  ,nyf,4)         +
     +                             f(2*ic  ,nyf,5) * b10   )
  210 continue
      do 220 ic=1,nxcm
      a10=wa(2*ic  ,nyf )
      b10=wb(2*ic  ,nyf )
c
c      (0,0) to (2,0)
      c(ic  ,nyc , 6)= c(ic  ,nyc , 6)                       +
     +                             f(2*ic-1,nyf,6) * a10     +
     +                   b10 * (   f(2*ic  ,nyf,5) * a10   +
     +                             f(2*ic  ,nyf,6)         )
  220 continue
      do 230 ic=1,nxcm
      a10=wa(2*ic  ,nyf )
      b10=wb(2*ic  ,nyf )
c
c      (2,0) to (0,0)
      c(ic+1,nyc , 4)= c(ic+1,nyc , 4)                       +
     +                   a10 * (   f(2*ic  ,nyf,4)         +
     +                             f(2*ic  ,nyf,5) * b10   ) +
     +                             f(2*ic+1,nyf,4) * b10
  230 continue
      do 240 ic=1,nxcm
      a10=wa(2*ic  ,nyf )
c
c      (2,0) to (2,0)
      c(ic+1,nyc , 5)= c(ic+1,nyc , 5)                       +
     +                   a10 * (   f(2*ic  ,nyf,5) * a10   +
     +                             f(2*ic  ,nyf,6)         +
     +                             f(2*ic+1,nyf,4)         ) +
     +                             f(2*ic+1,nyf,5)
  240 continue
c  ---------------------------------------------------------------------
c     upper right hand corner point at ( ic=nxc, jc= nycm+1=nyc )
c     already completed
c  ---------------------------------------------------------------------
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine decomp(levels, nm, nxf, nyf, a, ac, ldu, lduc,
     +                  dimwrk, work)
      integer levels, nm, nxf, nyf, dimwrk
      double precision
     +     a(nxf*nyf*9), ac(nm*9), ldu(nxf*nyf*3), lduc(nm*3),
     +     work(dimwrk)
c-----------------------------------------------------------------------
c     incomplete line l u decomposition of the matrices on all levels.
c
c     by  p.m. de zeeuw, version 911010
c-----------------------------------------------------------------------
      integer lev, npa, npl, nx, nxy, ny
      double precision cpb, cpe
      integer ngp, ngx, ngy, maxlev
      common /poi/ ngp(12),ngx(12),ngy(12),maxlev
      double precision cp
      common /cpu/ cp(5)
c
      call timing(cpb)
c
      nxy=nxf*nyf
      call illudc( nxf, nyf, a,
     +            ldu, ldu(1+nxy), ldu(1+2*nxy), work)
      do 200 lev=(levels-1), 1, -1
        nx=ngx(lev)
        ny=ngy(lev)
        nxy=nx*ny
        npa=9*ngp(lev)+1
        npl=3*ngp(lev)+1
        call illudc( nx, ny, ac(npa),
     +              lduc(npl), lduc(npl+nxy), lduc(npl+2*nxy), work)
  200 continue
      call timing(cpe)
      cp(3)=cp(3)+(cpe-cpb)
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine illudc( nx, ny, a, l, d, u, work)
      integer nx, ny
      double precision
     +     a(nx,ny,9), l(nx,ny), d(nx,ny), u(nx,ny), work(nx,12)
c-----------------------------------------------------------------------
c     performs illu-decomposition on a.
c     a remains intact, l,d,u are initialized with the decompositions
c            -
c     of the d
c             j
c
c     by  p.m. de zeeuw,  version 911010
c-----------------------------------------------------------------------
      integer j
c
      call tridec(nx,a(1,1,4),a(1,1,5),a(1,1,6),l,d,u)
      do 100 j=2,ny
        call blocks( nx*ny, a(1, j-1, 1), a(1, j  , 1),
     +                  nx, l(1, j-1)   , d(1, j-1)   , u(1, j-1),
     +                      l(1, j  )   , d(1, j  )   , u(1, j  ),
     +              work(1,1),work(1,2),work(1,3),work(1,4),work(1,5),
     +              work(1,6),work(1,7),work(1,8))
  100 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine tridec(nx,dm,dz,dp,lj,dj,uj)
      integer nx
      double precision
     +     dm(nx), dz(nx), dp(nx), lj(nx), dj(nx), uj(nx)
c-----------------------------------------------------------------------
c                                                 -
c     performs decomposition of tridiagonal block d .
c                                                  j
c
c     by  p.m. de zeeuw,  version 911010
c-----------------------------------------------------------------------
      double precision one, zero
      parameter( one = 1.0d0, zero = 0.0d0 )
      integer i, nx1
      double precision djim1
c
      dj(1)=one/dz(1)
      djim1=dj(1)
c     --------------------
c     loop 10 is recursive
c     --------------------
      lj(1)=zero
      do 10 i=2,nx
        lj(i)=-dm(i)*djim1
        djim1=one/(dz(i)+lj(i)*dp(i-1))
        dj(i)=djim1
   10 continue
      nx1=nx-1
      do 20 i=1,nx1
        uj(i)=-dp(i)*dj(i)
   20 continue
      uj(nx)=zero
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine blocks( dima, ajm1, aj  ,
     +                     nx,
     +                   ljm1, djm1, ujm1, lj  , dj  , uj  ,
     +                    qm3,  qm2,  qm1,  qze,  qp1,  qp2, qp3, ld  )
      integer dima, nx
      double precision
     +     ajm1(dima,9), aj(dima,9),
     +     ljm1(nx), djm1(nx), ujm1(nx),
     +     lj  (nx), dj  (nx), uj  (nx),
     +     qm3(nx),qm2(nx),qm1(nx), qze(nx), qp1(nx),qp2(nx),qp3(nx),
     +     ld(nx,5)
c-----------------------------------------------------------------------
c     performs illu-decomposition of j-th row of blocks of a.
c
c     by  p.m. de zeeuw,  version 911010
c-----------------------------------------------------------------------
      double precision zero
      parameter( zero = 0.0d0 )
      integer i, nxm, nx2
      double precision qzeip
      intrinsic mod
c
c-----------------------------------------------------------------------
c                                         - -1
c     first step - computation of 7-diag( d    ),
c                                          j-1
c              resulting diagonals are qm3, qm2, qm1, qze, qp1, qp2, qp3
c-----------------------------------------------------------------------
      nxm=nx-1
      do 10 i=1,nxm
        qze(i)=ujm1(i)*ljm1(i+1)
   10 continue
c     --------------------
c     loop 20 is recursive
c     --------------------
      qze(nx)=djm1(nx)
      qzeip=qze(nx)
      do 20 i=nxm,2,(-2)
        qzeip =  djm1(i)+qze(i)*qzeip
        qze(i)=  qzeip
        qzeip =  djm1(i-1)+qze(i-1)*qzeip
        qze(i-1)=qzeip
   20 continue
      if(mod(nx,2).eq.0) then
        qze(1)=  djm1(1)+qze(1)*qze(2)
      end if
c
      qm1(1)= zero
      qp1(1)=ujm1(1)*qze(2)
      do 30 i=2,nxm
        qm1(i)=ljm1(i)*qze(i)
        qp1(i)=ujm1(i)*qze(i+1)
   30 continue
      qm1(nx)=ljm1(nx)*qze(nx)
      qp1(nx)= zero
c
      qm2(1)= zero
      qp2(1)=ujm1(1)*qp1(2)
      do 40 i=2,nxm
        qm2(i)=ljm1(i-1)*qm1(i)
        qp2(i)=ujm1(i)*qp1(i+1)
   40 continue
      qm2(nx)=ljm1(nxm)*qm1(nx)
      qp2(nx)= zero
c
      qp3(1)=ujm1(1)*qp2(2)
      qp3(2)=ujm1(2)*qp2(3)
      qm3(1)= zero
      qm3(2)= zero
      do 50 i = 3, nxm
        qm3(i)=ljm1(i-2)*qm2(i  )
        qp3(i)=ujm1(i  )*qp2(i+1)
   50 continue
      qp3(nx)= zero
      qm3(nx)=ljm1(nx-2)*qm2(nx)
c-----------------------------------------------------------------------
c                                                     - -1
c     second step - computation of  5 diagonals of  l d
c                                                    j j-1
c-----------------------------------------------------------------------
      ld( 1,1)= zero
      do 51 i=2,nxm
        ld( i,1)=aj( i,1)*qm1(i-1)+aj( i,2)*qm2( i)+aj( i,3)*qm3(i+1)
   51 continue
      ld(nx,1)=  aj(nx,1)*qm1(nxm)+aj(nx,2)*qm2(nx)
      ld( 1,2)= zero
      do 52 i=2,nxm
        ld( i,2)=aj( i,1)*qze(i-1)+aj( i,2)*qm1( i)+aj( i,3)*qm2(i+1)
   52 continue
      ld(nx,2)=  aj(nx,1)*qze(nxm)+aj(nx,2)*qm1(nx)
      ld(1,3)=                     aj( 1,2)*qze( 1)+aj( 1,3)*qm1(  2)
      do 53 i=2,nxm
        ld( i,3)=aj( i,1)*qp1(i-1)+aj( i,2)*qze( i)+aj( i,3)*qm1(i+1)
   53 continue
      ld(nx,3)=  aj(nx,1)*qp1(nxm)+aj(nx,2)*qze(nx)
      ld( 1,4)=                    aj( 1,2)*qp1( 1)+aj( 1,3)*qze(  2)
      do 54 i=2,nxm
        ld( i,4)=aj( i,1)*qp2(i-1)+aj( i,2)*qp1( i)+aj( i,3)*qze(i+1)
   54 continue
      ld(nx,4)= zero
      ld( 1,5)=                    aj( 1,2)*qp2( 1)+aj( 1,3)*qp1(  2)
      do 55 i=2,nxm
        ld( i,5)=aj( i,1)*qp3(i-1)+aj( i,2)*qp2( i)+aj( i,3)*qp1(i+1)
   55 continue
      ld(nx,5)= zero
c-----------------------------------------------------------------------
c                                            -                 - -1
c     third and fourth step - computation of d = d - 3-diag( l d    u  )
c                                             j   j           j j-1  j-1
c     -
c     d  is represented by qm1, qze, qp1
c      j
c-----------------------------------------------------------------------
      qm1(1)= zero
      qm1(2)= aj(2,4)          -ld( 2,2)*ajm1(  1,8)-ld(2,3)*ajm1(  2,7)
      do 60 i=3,nx
        qm1( i)=aj( i,4)
     +     -ld(i,1)*ajm1(i-2,9)-ld( i,2)*ajm1(i-1,8)-ld(i,3)*ajm1(i  ,7)
   60 continue
      qze( 1)=aj( 1,5)         -ld( 1,3)*ajm1(  1,8)-ld(1,4)*ajm1(  2,7)
      do 70 i=2,nxm
        qze( i)=aj( i,5)
     +     -ld(i,2)*ajm1(i-1,9)-ld( i,3)*ajm1(i  ,8)-ld(i,4)*ajm1(i+1,7)
   70 continue
      qze(nx)=aj(nx,5)
     +    -ld(nx,2)*ajm1(nxm,9)-ld(nx,3)*ajm1(nx ,8)
      nx2=nx-2
      do 80 i=1,nx2
        qp1( i)=aj( i,6)
     +     -ld(i,3)*ajm1(i  ,9)-ld( i,4)*ajm1(i+1,8)-ld(i,5)*ajm1(i+2,7)
   80 continue
      qp1(nxm)=aj(nxm,6)
     +   -ld(nxm,3)*ajm1(nxm,9)-ld(nxm,4)*ajm1( nx,8)
      qp1(nx)= zero
c-----------------------------------------------------------------------
c                                                              -
c     fifth step -  computation of decomposition l ,d ,u   of  d
c                                                 j  j  j       j
c-----------------------------------------------------------------------
      call tridec(nx,qm1,qze,qp1,lj,dj,uj)
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine residu(nx, nxy, f, a, u, r)
      integer nx, nxy
      double precision
     +     f(nxy), a(nxy,9), u(nxy), r(nxy)
c-----------------------------------------------------------------------
c     computes r = f - a * u
c
c     ( this routine is the vectorized version of residu )
c
c     by  p.m. de zeeuw,  version 911010
c-----------------------------------------------------------------------
      integer i, nxm1, nxp1, nxp2, nxynx, nxynxm, nxy1
c
      nxy1=nxy-1
      nxp1=nx+1
      nxp2=nx+2
      nxm1=nx-1
      nxynx=nxy-nx
      nxynxm=nxy-nx-1
c
      r(1)=f(1)-a(1,5)*u(1)-a(1,6)*u(2)-a(1,8)*u(nxp1)-a(1,9)*u(nxp2)
      do 10 i=2,nxy1
        r(i)=f(i)-a(i,4)*u(i-1)-a(i,5)*u(i)-a(i,6)*u(i+1)
   10 continue
      r(nxp1)=r(nxp1)-a(nxp1,2)*u(1)-a(nxp1,3)*u(2)
      do 20 i=nxp2,nxy1
        r(i)=r(i)-a(i,1)*u(i-nxp1)-a(i,2)*u(i-nx)-a(i,3)*u(i-nxm1)
   20 continue
      do 30 i=2,nxynxm
        r(i)=r(i)-a(i,7)*u(i+nxm1)-a(i,8)*u(i+nx)-a(i,9)*u(i+nxp1)
   30 continue
      r(nxynx)=r(nxynx)-a(nxynx,7)*u(nxy1)-a(nxynx,8)*u(nxy)
      r(nxy)=f(nxy)-a(nxy,4)*u(nxy1)  -a(nxy,5)*u(nxy)
     +             -a(nxy,1)*u(nxynxm)-a(nxy,2)*u(nxynx)
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine smooth(level,v,vs,a,rhs,ldu,w,resvs)
      integer level
      double precision
     +     v(*), vs(*), a(*), rhs(*), ldu(*), w(*)
      logical resvs
c-----------------------------------------------------------------------
c     if logical resvs equals .true. then vs is expected to contain
c     the residual rhs-a*v as input, else vs is just scratch.
c     the purpose of this subroutine is to overwrite v with a better
c     (smoothed) approximation of   a * v     =  rhs
c                                        exact
c     by performing an incomplete line l u relaxation sweep.
c-----------------------------------------------------------------------
      integer ir, k, na, ne, nl, np, nx, nxy, ny
      integer ngp, ngx, ngy, maxlev
      common /poi/ ngp(12),ngx(12),ngy(12),maxlev
c
      nx=ngx(level)
      ny=ngy(level)
      nxy=nx*ny
      np=ngp(level)+1
      ne=  (np-1)+nxy
      na=9*(np-1)+1
      nl=3*(np-1)+1
      if((level.eq.1).and.(maxlev.gt.1)) then
         call copy( nxy, rhs(np), v(np))
         call solve( a(na),ldu(nl),ldu(nl+nxy),ldu(nl+2*nxy),
     +               v(np),w,w(nx+1),nx,ny)
         do 7 ir = 1, 8
           call residu(nx,nxy,rhs(np),a(na),v(np),vs(np))
           call solve( a(na),ldu(nl),ldu(nl+nxy),ldu(nl+2*nxy),
     +                vs(np),w,w(nx+1),nx,ny)
           do 5 k=np,ne
             v(k)=v(k)+vs(k)
    5      continue
    7    continue
      else
           if(.not.resvs) then
              call residu(nx,nxy,rhs(np),a(na),v(np),vs(np))
           end if
           call solve( a(na),ldu(nl),ldu(nl+nxy),ldu(nl+2*nxy),
     +                vs(np),w,w(nx+1),nx,ny)
         do 10 k=np,ne
           v(k)=v(k)+vs(k)
   10    continue
      end if
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine solve(a,l,d,u,r,t,tt,nx,ny)
      integer nx, ny
      double precision
     +     a(nx*ny,9), l(nx*ny), d(nx*ny), u(nx*ny),
     +     r(nx*ny), t(nx), tt(nx)
c-----------------------------------------------------------------------
c                  -   - -1  -
c     solves ( l + d ) d   ( d + u ) r      = r
c                                     (new)    (old)
c
c     ( this routine is the vectorized version of solve )
c
c     by  p.m. de zeeuw,  version 911010
c-----------------------------------------------------------------------
      integer i, ib, ibm1, ie, j, nxm, nxp
      double precision
     +     ri, rii, rip, ri1, ti, tip, ti1, ti1p, tti1
      logical even
      intrinsic mod
c
c-----------------------------------------------------------------------
c              - -1
c     r      = d    r
c      1        1    1
c       (new)         (old)
c-----------------------------------------------------------------------
      even=(mod(nx,2).eq.0)
      rip=r(1)
      nxm=nx-1
      nxp=nx+1
      do 10 i=2,nxm,2
        ri=r(i)+l(i)*rip
        r(i)=ri
        rip=r(i+1)+l(i+1)*ri
        r(i+1)=rip
   10 continue
      if(even) then
        r(nx)=r(nx)+l(nx)*r(nxm)
      end if
      do 20 i=1,nx
        r(i)=r(i)*d(i)
   20 continue
      rii=r(nx)
      do 30 i=nxm,2,(-2)
        rii=r(i)+u(i)*rii
        r(i)=rii
        rii=r(i-1)+u(i-1)*rii
        r(i-1)=rii
   30 continue
      if(even) then
        r(1)=r(1)+u(1)*r(2)
      end if
c-----------------------------------------------------------------------
c              - -1
c     r      = d   ( r      - l  r        )         j=2(1)ny
c      j        j     j        j  j-1
c       (new)          (old)         (new)
c-----------------------------------------------------------------------
      ib=1
      ie=nx
      do 1000 j=2,ny
        ib=ib+nx
        ie=ie+nx
        ibm1=ib-1
        t(1)=     r(ib)             -a(ib,2)*r(ib-nx)-a(ib,3)*r(ib-nxm)
        do 100 i=(ib+1),ie
          t(i-ibm1)=r(i)-a(i,1)*r(i-nxp)-a(i,2)*r(i-nx)-a(i,3)*r(i-nxm)
  100   continue
        tip=t(1)
        do 200 i=2,(nx-1),2
          ti=t(i)+l(i+ibm1)*tip
          t(i)=ti
          tip=t(i+1)+l(i+ibm1+1)*ti
          t(i+1)=tip
  200   continue
        if(even) then
           t(nx)=t(nx)+l(nx+ibm1)*t(nx-1)
        end if
        r(ie)=t(nx)*d(ie)
        ri1=r(ie)
        do 400 i=(ie-1),ib,(-1)
          ri1=t(i-ibm1)*d(i)+u(i)*ri1
          r(i)=ri1
  400   continue
 1000 continue
c-----------------------------------------------------------------------
c                       - -1
c     r      = r      - d   u   r                j=ny-1(-1)1
c      j        j        j   j   j+1
c       (new)    (old)              (new)
c-----------------------------------------------------------------------
      do 2000 j=2,ny
        ib=ib-nx
        ie=ie-nx
        ibm1=ib-1
        do 1100 i=ib,(ie-1)
          t(i-ibm1)=a(i,7)*r(i+nxm)+a(i,8)*r(i+nx)+a(i,9)*r(i+nxp)
 1100   continue
        t(nx)=  a(ie,7)*r(ie+nxm)+a(ie,8)*r(ie+nx)
        ti1p=t(1)
        do 1200 i=2,(nx-1),2
          ti1=t(i)+l(i+ibm1)*ti1p
          t(i)=ti1
          ti1p=t(i+1)+l(i+ibm1+1)*ti1
          t(i+1)=ti1p
 1200   continue
        if(even) then
          t(nx)=t(nx)+l(nx+ibm1)*t(nx-1)
        end if
        tti1=t(nx)*d(ie)
        tt(nx)=tti1
        do 1400 i=(nx-1),1,(-1)
          tti1=t(i)*d(i+ibm1)+u(i+ibm1)*tti1
          tt(i)=tti1
 1400   continue
        do 1900 i=ib,ie
          r(i)=r(i)-tt(i-ibm1)
 1900   continue
 2000 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      double precision function vl2nor(u,np)
      integer np
      double precision u(np)
c-----------------------------------------------------------------------
c     computes l2-norm of vector u.
c-----------------------------------------------------------------------
      double precision zero
      parameter( zero = 0.0d0 )
      integer k
      double precision v
      intrinsic sqrt
c
      v= zero
      do 10 k=1,np
        v=v+u(k)*u(k)
   10 continue
      vl2nor= sqrt(v)
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine prolon(nxc,nyc,c,nxf,nyf,wa,wb,f)
      integer nxc, nyc, nxf, nyf
      double precision
     +     c(nxc,nyc),  wa(nxf,nyf), wb(nxf,nyf), f(nxf,nyf)
c-----------------------------------------------------------------------
c
c     prolongates coarse-grid-function c onto fine-grid-function f
c
c     by p.m. de zeeuw, version 911010
c
c-----------------------------------------------------------------------
      integer i, j, jf
c
      do 20 j = 1, nyc
        jf=2*j-1
        do 10 i = 1, nxc
          f(2*i-1,jf)= c(i,j)
   10   continue
   20 continue
c
      do 40 j = 1, nyc
        jf=2*j-1
        do 30 i = 1, (nxc-1)
          f(2*i  ,jf)=
     +                 wb(2*i  ,jf)*c(i  ,j  )+wa(2*i  ,jf)*c(i+1,j  )
   30   continue
   40 continue
c
      do 60 j = 1, (nyc-1)
        jf=2*j
        do 50 i = 1, nxc
          f(2*i-1,jf)=
     +                 wb(2*i-1,jf)*c(i  ,j+1)+wa(2*i-1,jf)*c(i  ,j  )
   50 continue
   60 continue
c
      do 80 j = 1, (nyc-1)
        jf=2*j
        do 70 i = 1, (nxc-1)
          f(2*i,jf)=
     +             wb(2*i  ,jf  )*c(i  ,j+1)+wa(2*i-1,jf-1)*c(i+1,j+1)+
     +             wb(2*i-1,jf-1)*c(i  ,j  )+wa(2*i  ,jf  )*c(i+1,j  )
   70   continue
   80 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine restri(nxf, nyf, wa, wb, f, nxc, nyc, c)
      integer nxf, nyf, nxc, nyc
      double precision
     +     wa(nxf,nyf), wb(nxf,nyf), f(nxf,nyf), c(nxc,nyc)
c-----------------------------------------------------------------------
c     computes c = restriction of f
c     according to the molecule
c
c                                y
c                                +  wa   wa    wd
c                                +  wa    1    wb
c                                +  wc   wb    wb
c                                +
c                                + + + + + + + + + x
c
c     wc has been stored in wa and wd has been stored in wb as follows
c     wc(i,j) in wa(i-1,j-1) and wd(i,j) in wb(i-1,j-1).
c
c     by p.m. de zeeuw, version 911010
c-----------------------------------------------------------------------
      integer ic, jc, jf, jfm, jfp, nxc1, nyc1
c
c     boundary below
      c(1, 1)= f(1, 1)+wb(2,1)*f(2,1)+wa(1,2)*f(1,2)+wb(1,1)*f(2,2)
      nxc1=nxc-1
      do 10 ic = 2, nxc1
         c(ic, 1)=  f(2*ic-2, 1) * wa(2*ic-2, 1) +
     +              f(2*ic-1, 1)                 +
     +              f(2*ic  , 1) * wb(2*ic  , 1) +
     +              f(2*ic-2, 2) * wa(2*ic-2, 2) +
     +              f(2*ic-1, 2) * wa(2*ic-1, 2) +
     +              f(2*ic  , 2) * wb(2*ic-1, 1)
   10 continue
      c(nxc, 1)=  f(nxf-1, 1) * wa(nxf-1, 1) +
     +            f(nxf  , 1)                +
     +            f(nxf-1, 2) * wa(nxf-1, 2) +
     +            f(nxf  , 2) * wa(nxf  , 2)
c     inner area
      nyc1=nyc-1
      do 30 jc = 2, nyc1
         jf =2*jc-1
         jfm=jf-1
         jfp=jf+1
         c( 1, jc)=  f(     1, jfm) * wb(     1, jfm) +
     +               f(     2, jfm) * wb(     2, jfm) +
     +               f(     1, jf )                   +
     +               f(     2, jf ) * wb(     2, jf ) +
     +               f(     1, jfp) * wa(     1, jfp) +
     +               f(     2, jfp) * wb(     1, jf )
         do 20 ic = 2, nxc1
          c(ic, jc)=  f(2*ic-2, jfm) * wa(2*ic-3, jfm-1) +
     +                f(2*ic-1, jfm) * wb(2*ic-1, jfm) +
     +                f(2*ic  , jfm) * wb(2*ic  , jfm) +
     +                f(2*ic-2, jf ) * wa(2*ic-2, jf ) +
     +                f(2*ic-1, jf )                   +
     +                f(2*ic  , jf ) * wb(2*ic  , jf ) +
     +                f(2*ic-2, jfp) * wa(2*ic-2, jfp) +
     +                f(2*ic-1, jfp) * wa(2*ic-1, jfp) +
     +                f(2*ic  , jfp) * wb(2*ic-1, jf )
   20    continue
         c(nxc,jc)=  f(nxf-1, jfm) * wa(nxf-2, jfm-1) +
     +               f(nxf  , jfm) * wb(nxf  , jfm) +
     +               f(nxf-1, jf ) * wa(nxf-1, jf ) +
     +               f(nxf  , jf )                  +
     +               f(nxf-1, jfp) * wa(nxf-1, jfp) +
     +               f(nxf  , jfp) * wa(nxf  , jfp)
   30 continue
c     upper boundary
      c(  1, nyc)=  f(  1, nyf-1) * wb(  1, nyf-1) +
     +              f(  2, nyf-1) * wb(  2, nyf-1) +
     +              f(  1, nyf  )                  +
     +              f(  2, nyf  ) * wb(  2, nyf  )
      do 40 ic = 2, nxc1
         c(ic,nyc)= f(2*ic-2,nyf-1) * wa(2*ic-3,nyf-2) +
     +              f(2*ic-1,nyf-1) * wb(2*ic-1,nyf-1) +
     +              f(2*ic  ,nyf-1) * wb(2*ic  ,nyf-1) +
     +              f(2*ic-2,nyf  ) * wa(2*ic-2,nyf  ) +
     +              f(2*ic-1,nyf  )                    +
     +              f(2*ic  ,nyf  ) * wb(2*ic  ,nyf  )
   40 continue
      c(nxc,nyc)=   f(nxf-1, nyf-1) * wa(nxf-2, nyf-2) +
     +              f(nxf  , nyf-1) * wb(nxf  , nyf-1) +
     +              f(nxf-1, nyf  ) * wa(nxf-1, nyf  ) +
     +              f(nxf  , nyf  )
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine zeros( n, a)
      integer n
      double precision a(n)
      double precision zero
      parameter( zero = 0.0d0 )
      integer nn
c
      do 10 nn = 1, n
        a(nn)= zero
   10 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine copy( n, a, ca)
      integer n
      double precision a(n), ca(n)
      integer nn
c
      do 10 nn = 1, n
        ca(nn)=a(nn)
   10 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine outmat(nout,a,nx,ny,nd)
      integer nout, nx, ny, nd
      double precision a(nx,ny,nd)
c-----------------------------------------------------------------------
c     prints matrix a with nd diagonals on file tapenout
c-----------------------------------------------------------------------
      integer i, j, k
c
      do 10 j = 1, ny
        write(nout,2) j
    2   format(' y-index=',i4)
        do 10 i = 1, nx
          write(nout,3) i,(a(i,j,k),k=1,nd)
    3     format(' x-index=',i4,1x,9d11.3)
   10 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine outvec(nout,u,nx,ny)
      integer nout, nx, ny
      double precision u(nx,ny)
c-----------------------------------------------------------------------
c     prints vector u on file tapenout.
c-----------------------------------------------------------------------
      integer i, j
c
      do 10 j = 1, ny
        write(nout,2)j,(u(i,j),i=1,nx)
    2   format(' y-index=',i4,1x,10d11.3/(12x,10d11.3))
   10 continue
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine timing(cp)
      double precision cp
c-----------------------------------------------------------------------
c     important  user provided time measuring  (machine dependent)
c     ---------
c
c-----------------------------------------------------------------------
      cp= 0.0d0
c     if cyber .and. fortran 77 then
c       cp=second()
c     end if
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      integer function p01cwi(ifail,ierror,srname)
      integer ifail, ierror
      character*8 srname
      integer ia, ib, nerr
      intrinsic mod
c     to err is human
c
      p01cwi=ierror
      ia=mod(ifail,10)
      ib=mod(ifail/10,10)
c     ia is the last digit of ifail
      call x04aaf(0,nerr)
      if (ib.eq.1) then
        write(nerr,1) srname,ierror
    1   format(/' error detected by subroutine ',a8,
     +        ' - ifail = ',i2)
      end if
      if (ia.eq.1) then
        return
      else
        stop
      end if
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine x04aaf(iflag,nerr)
      integer iflag, nerr
      integer nsav
      save nsav
c
      data nsav /6/
      if (iflag.eq.0) then
         nerr=nsav
      else
         nsav=nerr
      end if
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine x04abf(iflag,nadv)
      integer iflag, nadv
      integer nsav
      save nsav
c
      data nsav /6/
      if (iflag.eq.0) then
         nadv=nsav
      else
         nsav=nadv
      end if
      return
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
