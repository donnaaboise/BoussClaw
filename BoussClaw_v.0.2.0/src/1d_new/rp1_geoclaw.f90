!!======================================================================
SUBROUTINE rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,&
               qr,auxl,auxr,fwave,s,amdq,apdq)

!!======================================================================
!!
!! Solves normal Riemann problems for the 2D SHALLOW WATER equations
!!     with topography:
!!     #        h_t + (hu)_x + (hv)_y = 0                           #
!!     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
!!     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #
!!
!! On input, ql contains the state vector at the left edge of each cell
!!     qr contains the state vector at the right edge of each cell
!!
!! This data is along a slice in the x-direction if ixy=1
!!     or the y-direction if ixy=2.
!!
!!  Note that the i'th Riemann problem has left state qr(i-1,:)
!!     and right state ql(i,:)
!!  From the basic clawpack routines, this routine is called with
!!     ql = qr
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water equations.            !
!                                                                           !
!       It allows the user to easily select a Riemann solver in             !
!       riemannsolvers_geo.f. this routine initializes all the variables    !
!       for the shallow water equations, accounting for wet dry boundary    !
!       dry cells, wave speeds etc.                                         !
!                                                                           !
!           David George, Vancouver WA, Feb. 2009                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!  use geoclaw_module, only: g => grav, drytol => dry_tolerance, rho
!!  use geoclaw_module, only: earth_radius, deg2rad
!!  use amr_module, only: mcapa
!!
!!  use storm_module, only: pressure_forcing, pressure_index

  USE geoclaw_module, only: g => grav, drytol => dry_tolerance, rho  
  USE geoclaw_module, only: earth_radius, deg2rad
  USE amr_module, only: mcapa

  IMPLICIT NONE

  !input
  INTEGER maxm,meqn,maux,mwaves,mbc,mx

  DOUBLE PRECISION  fwave(meqn, mwaves, 1-mbc:maxm+mbc)
  DOUBLE PRECISION  s(mwaves, 1-mbc:maxm+mbc)
  DOUBLE PRECISION  ql(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION  qr(meqn, 1-mbc:maxm+mbc)
  DOUBLE PRECISION  apdq(meqn,1-mbc:maxm+mbc)
  DOUBLE PRECISION  amdq(meqn,1-mbc:maxm+mbc)
  DOUBLE PRECISION  auxl(maux,1-mbc:maxm+mbc)
  DOUBLE PRECISION  auxr(maux,1-mbc:maxm+mbc)


  !local only
  INTEGER m,i,mw,maxiter,mu,nv
  DOUBLE PRECISION wall(2)
  DOUBLE PRECISION fw(2,2)
  DOUBLE PRECISION sw(2)

  DOUBLE PRECISION hR,hL,huR,huL,uR,uL,phiR,phiL,pL,pR
  DOUBLE PRECISION bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
  double precision hvl, hvr, vl, vr
  DOUBLE PRECISION s1m,s2m
  DOUBLE PRECISION hstar,hstartest,hstarHLL,sLtest,sRtest
  DOUBLE PRECISION tw,dxdc

  LOGICAL rare1,rare2

  pL = 0
  pR = 0


  !loop through Riemann problems at each grid cell
  DO i = 2-mbc,mx+mbc

      !! -----------------------Initializing-----------------------------------
      !! inform of a bad riemann problem from the start
     IF((qr(1,i-1).LT.0.d0).OR.(ql(1,i) .LT. 0.d0)) THEN
        WRITE(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
     ENDIF

     !!Initialize Riemann problem for grid interface
     DO mw=1,mwaves
        s(mw,i)=0.d0
        fwave(1,mw,i)=0.d0
        fwave(2,mw,i)=0.d0
     ENDDO

     !!zero (small) negative values if they exist
     IF (qr(1,i-1).LT.0.d0) THEN
        qr(1,i-1)=0.d0
        qr(2,i-1)=0.d0
     ENDIF

     IF (ql(1,i).LT.0.d0) THEN
        ql(1,i)=0.d0
        ql(2,i)=0.d0
     ENDIF

     !!skip problem if in a completely dry area
     IF (qr(1,i-1) <= drytol .AND. ql(1,i) <= drytol) THEN
        go to 30
     ENDIF

     !! Riemann problem variables
     hL = qr(1,i-1)
     hR = ql(1,i)
     huL = qr(2,i-1)
     huR = ql(2,i)
     bL = auxr(1,i-1)
     bR = auxl(1,i)

     !!check for wet/dry boundary
     IF (hR .GT. drytol) THEN
        uR = huR/hR
        phiR = 0.5d0*g*hR**2 + huR**2/hR
     ELSE
        hR = 0.d0
        huR = 0.d0
        uR = 0.d0
        phiR = 0.d0
     ENDIF

     IF (hL .gt. drytol) THEN
        uL = huL/hL
        phiL = 0.5d0*g*hL**2 + huL**2/hL
     ELSE
        hL=0.d0
        huL=0.d0
        uL=0.d0
        phiL = 0.d0
     ENDIF

     wall(1) = 1.d0
     wall(2) = 1.d0
     IF (hR .LE. drytol) THEN
        CALL riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m, &
             rare1,rare2,1,drytol,g)

        hstartest = MAX(hL,hstar)
        IF (hstartest+bL .LT. bR) THEN
           !!right state should become ghost values that mirror left for wall problem
           !! bR=hstartest+bL
           wall(2) = 0.d0
           hR = hL
           huR = -huL
           bR = bL
           phiR = phiL
           uR = -uL
        ELSEIF (hL+bL.LT.bR) THEN
           bR=hL+bL
        ENDIF
     ELSEIF (hL .LE. drytol) THEN ! right surface is lower than left topo
        CALL riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m, &
             rare1,rare2,1,drytol,g)
        hstartest=MAX(hR,hstar)
        IF (hstartest+bR.LT.bL) THEN
           !!left state should become ghost values that mirror right
           !! bL=hstartest+bR
           wall(1) = 0.d0
           wall(2) = 0.d0
           hL = hR
           huL = -huR
           bL = bR
           phiL = phiR
           uL = -uR
        ELSEIF (hR + bR .LT. bL) THEN
           bL = hR + bR
        ENDIF
     ENDIF

     !!determine wave speeds
     sL = uL - SQRT(g*hL) ! 1 wave speed of left state
     sR = uR + SQRT(g*hR) ! 2 wave speed of right state

     uhat = (SQRT(g*hL)*uL + SQRT(g*hR)*uR)/(SQRT(g*hR)+SQRT(g*hL)) ! Roe average
     chat = SQRT(g*0.5d0*(hR+hL)) ! Roe average
     sRoe1 = uhat-chat ! Roe wave speed 1 wave
     sRoe2 = uhat+chat ! Roe wave speed 2 wave

     sE1 = MIN(sL,sRoe1) ! Eindfeldt speed 1 wave
     sE2 = MAX(sR,sRoe2) ! Eindfeldt speed 2 wave

     !!--------------------end initializing...finally----------
     !!solve Riemann problem.

     maxiter = 1
     hvL = 0
     vL = 0
     vR = 0

     CALL riemann_aug_JCP(maxiter,3,3,hL,hR,huL, &
          huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2, &
          drytol,g,rho,sw,fw)

     !!CALL riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR, &
     !!     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g, &
     !!     rho,sw,fw)

     !! CALL riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR, &
     !!     bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2,drytol,g,rho,sw, &
     !!     fw)

     !! eliminate ghost fluxes for wall
     DO mw=1,mwaves
        sw(mw) = sw(mw)*wall(mw)
        fw(1,mw) = fw(1,mw)*wall(mw)
        fw(2,mw) = fw(2,mw)*wall(mw)
     ENDDO

     DO mw = 1,mwaves
        s(mw,i) = sw(mw)
        fwave(1,mw,i) = fw(1,mw)
        fwave(2,mw,i) = fw(2,mw)
        !!            write(51,515) sw(mw),fw(1,mw),fw(2,mw),fw(3,mw)
        !!515         format("++sw",4e25.16)
     ENDDO

30   CONTINUE
  ENDDO


  !!===============================================================================


  !!============= compute fluctuations=============================================
  amdq(1:2,:) = 0.d0
  apdq(1:2,:) = 0.d0
  DO i = 2-mbc,mx+mbc
     DO  mw = 1,mwaves
        IF (s(mw,i) < 0.d0) THEN
           amdq(:,i) = amdq(:,i) + fwave(:,mw,i)
        ELSE IF (s(mw,i) > 0.d0) THEN
           apdq(:,i)  = apdq(:,i) + fwave(:,mw,i)
        ELSE
           amdq(:,i) = amdq(:,i) + 0.5d0 * fwave(:,mw,i)
           apdq(:,i) = apdq(:,i) + 0.5d0 * fwave(:,mw,i)
        ENDIF
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE rp1
