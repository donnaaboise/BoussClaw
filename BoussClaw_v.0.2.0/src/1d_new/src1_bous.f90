
  subroutine src1_bous(meqn,mbc,mx, &
                  xlower,dx,q,maux,aux,t,dt)

    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient,sea_level
    use bous_module

    implicit none
            
    integer(kind=4) i,mx,my,k,ii,nstep,kk
    integer(kind=4), parameter :: ndim = 1
      
    real(kind=8)   q(meqn,1-mbc:mx+mbc)
    real(kind=8)  q0(meqn,1-mbc:mx+mbc)
    real(kind=8) aux(maux,1-mbc:mx+mbc)
    real(kind=8) psi(mx+2)
     
    real(kind=8)   xx0(1:mx,4)
            
    real(kind=8) :: delt

    real(kind=8) tol,t,dt,x,y
    real(kind=8) dx,dy,xlower,ylower
    integer(kind=4) INFO,maux,mbc,meqn
     
    INTEGER            LDB
    INTEGER            IPIV( mx+2 )
    DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1), DU2(1:mx)             
      
    logical verbose

    verbose=.False.
           
    !! -------------------------------------------------
    !!    Boussinesq type
    !! -------------------------------------------------
    
    nstep=ceiling(dt/dx*2*sqrt(2.))

    delt=dt/nstep

    LDB = mx+2
    D  =0.d0
    DU =0.d0
    DL =0.d0
    DU2=0.d0
    
    do ii = 1,nstep

        call read_diag(mx,meqn,mbc,dx,q,maux,aux, DL, D, DU)

        call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )
    
        psi = 0.d0
        xx0 = 0.d0
        q0  = q

        !! RK4   

        !! First Stage
      
        call read_psi(mx,meqn,mbc,dx,q0,maux,aux,psi,g)

        call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

        XX0(1:mx,1) = psi(2:mx+1)
        
        q0(2,1:mx)=q(2,1:mx)-delt/2.d0*xx0(1:mx,1)
        
        !! Second Stage

        call read_psi(mx,meqn,mbc,dx,q0,maux,aux,psi,g )

        call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                        , IPIV, psi, LDB,INFO )

        XX0(1:mx,2)=psi(2:mx+1)
        
        q0(2,1:mx)=q(2,1:mx)-delt/2.d0*xx0(1:mx,2)
        
        !! Third Stage
        
        call read_psi(mx,meqn,mbc,dx,q0,maux,aux,psi,g )

        call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )
                
        XX0(1:mx,3)=psi(2:mx+1)
        
        q0(2,1:mx)=q(2,1:mx)-delt*xx0(1:mx,3)
        
        !! Fourth Stage
        
        call read_psi(mx,meqn,mbc,dx,q0,maux,aux,psi,g )

        call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2 &
                      , IPIV, psi, LDB,INFO )

        XX0(1:mx,4)=psi(2:mx+1)

        !!=======================================================================

        q(2,1:mx) = q(2,1:mx)- delt/6.d0*(xx0(1:mx,1) + 2.d0*xx0(1:mx,2) &
                                + 2.d0*xx0(1:mx,3) + xx0(1:mx,4))

    enddo
      
      
  end subroutine src1_bous
  
!! =========================================================
subroutine read_diag(mx,meqn,mbc,dx,q,maux,aux,DL,D,DU)
!! =========================================================
    use bous_module
    use geoclaw_module, only: sea_level

    implicit none
     
    integer(kind=4) mx,meqn,mbc,maux,mxy
    integer(kind=4) i,k
    real(kind=8) dx

    real(kind=8)   q(meqn,1-mbc:mx+mbc)
    real(kind=8) aux(maux,1-mbc:mx+mbc)

    DOUBLE PRECISION   D(mx+2), DL(mx+1), DU(mx+1)                
      
    real(kind=8)  HH(1-mbc:mx+mbc)
    real(kind=8) HHH(1-mbc:mx+mbc)
    real(kind=8) slope(1-mbc:mx+mbc)
      
    slope = aux(1,:)-sea_level  
     
    do i = 1-mbc,mx+mbc
        HH(i) = (max(0.,-slope(i)))**2
        HHH(i) = (max(0.,-slope(i)))**3
    enddo
      
    D = 1.d0
    DU= 0.d0
    DL= 0.d0

    do i = 1,mx
        
        if (q(1,i).lt.1.d-1 .or. slope(i)>-1.d-1) then

            !! do nothing
        
        elseif (maxval(slope(i-2:i+2))  >-1d-1) then

            !! do nothing

        elseif (maxval(q(1,i-2:i+2)) < 1.d-1) then

            !! do nothing

        elseif (maxval(slope(i-1:i+1)) < -1d-1 ) then              

            !!! D1 part     

            D(i+1) = 1.d0 + 2.d0*(B_param+.5d0)*HH(i)/dx**2 &
                   -2.d0/6.d0*HHH(i)/(-slope(i))/dx**2
     
            DU(i+1)=-(B_param+.5d0)*HH(i)/dx**2 &
              +1.d0/6.d0*HHH(i)/(-slope(i+1))/dx**2
     
            DL(i)=-(B_param+.5d0)*HH(i)/dx**2 &
              +1.d0/6.d0*HHH(i)/(-slope(i-1))/dx**2
     
        endif
            
    enddo
             
    return
end subroutine read_diag
   
!======================================================================
subroutine read_psi(mx,meqn,mbc,dx,q,maux,aux,psi,g )

    use geoclaw_module, only: sea_level
    use bous_module
     
    implicit none

    integer(kind=4) mx,meqn,mbc,maux,i,k,kk,iL,iR
    real   (kind=8) dx,dy,g
     
    real(kind=8)   q(meqn,1-mbc:mx+mbc)
    real(kind=8) aux(maux,1-mbc:mx+mbc)
    real(kind=8) slope(1-mbc:mx+mbc)
    real(kind=8) psi(mx+2)

    real(kind=8)  hh(1-mbc:mx+mbc)
    real(kind=8)  hhh(1-mbc:mx+mbc)
    real(kind=8)  hetax(-1:mx+2)
    real(kind=8)  hetay(-1:mx+2)
    real(kind=8)  eta(1-mbc:mx+mbc)
    real(kind=8)  hu2(1-mbc:mx+mbc)
    real(kind=8)  detax(-1:mx+2)
    real(kind=8)  s1(0:mx+1)
    real(kind=8)  s2(0:mx+1)
    real(kind=8)  s1_H(0:mx+1)
    real(kind=8)  s2_H(0:mx+1)
    real(kind=8)  tol,topo
    
    slope = aux(1,:)-sea_level
     
    do i=1-mbc,mx+mbc
        if (q(1,i) .gt. 1d-4) then
            hu2(i)= q(2,i)**2/q(1,i)
        else
            hu2(i)= 0.d0
        endif
        eta(i) = q(1,i) + slope(i) 
        HH(i) = (max(0.d0,-slope(i)))**2
        HHH(i) = (max(0.d0,-slope(i)))**3
    enddo
     
    hetax = 0.d0
    detax = 0.d0

    do i= 1,mx
        if (slope(i+1) < 0.d0.and.slope(i-1) < 0.d0.and.slope(i) < 0.d0) then
            if (i == 1) then
                hetax(i) = q(1,i)*(eta(i+1)-eta(i))/dx
                detax(i) = -slope(i)*(eta(i+1)-eta(i))/dx
            elseif (i==mx) then
                hetax(i) = q(1,i)*(eta(i)-eta(i-1))/dx
                detax(i) = -slope(i)*(eta(i)-eta(i-1))/dx
            else
                hetax(i)= q(1,i)*(eta(i+1)-eta(i-1))/dx/2.
                detax(i)=-slope(i)*(eta(i+1)-eta(i-1))/dx/2.
            endif
        endif
    enddo
     
    s1=0.d0
    s2=0.d0
    s1_H=0.d0
    s2_H=0.d0

    do i=1,mx
        if (i==1) then
            s1(i)= (hu2(i+1)-hu2(i))/dx
        elseif (i==mx) then
            s1(i)= (hu2(i)-hu2(i-1))/dx
        else
            s1(i)= (hu2(i+1)-hu2(i-1))/2.d0/dx
        endif
                
        if (slope(i) < -1d-1) then
            s1_H(i) = s1(i)/(-slope(i))
        else
            s1_H(i) = 0.d0
        endif
    enddo
     
    tol = 1d-8

    psi=0.d0

    do i=1,mx

        k = i+1

        topo =-slope(i)

        if (q(1,i) .lt. 1.d-1) then
            psi(k) = 0.d0

        elseif (slope(i)>-1.d-1) then
            psi(k) = 0.d0
         
        elseif (use_bous_sw_thresh.and.abs(q(1,i)-topo)>bous_sw_thresh*topo) then
            psi(k) = 0.d0
              !go to 994

        !elseif (-eta(i,j)>.1*q(1,i,j)) then
        !   psi(k) = 0.d0

        elseif (maxval(slope(i-2:i+2))> sw_depth ) then
            
            psi(k)=0.d0

        elseif (maxval(q(1,i-2:i+2)) < 1.d-1) then

            psi(k) = 0.d0

        else
            if (i == 1) then
                psi(k) = (B_param +.5d0)*hh(i)*(s1(i+2)-2.d0*s1(i+1)+s1(i))/dx**2 &
                              -B_param*g*HH(i)*(detax(i+2)-2.d0*detax(i+1)+detax(i))/dx**2 &
                              -HHH(i)/6.d0*(s1_H(i+2)-2.d0*s1_H(i+1)+s1_H(i))/dx**2
              
            elseif (i == mx) then

                psi(k)=(B_param +.5d0)*hh(i)*(s1(i)-2.d0*s1(i-1)+s1(i-2))/dx**2. &
                       -B_param*g*HH(i)*(detax(i)-2.d0*detax(i-1)+detax(i-2))/dx**2 &
                         -HHH(i)/6.d0*(s1_H(i)-2.d0*s1_H(i-1) + s1_H(i-2))/dx**2                          
            else
               
                psi(k)=(B_param+.5d0)*hh(i)*(s1(i+1)-2.d0*s1(i)+s1(i-1))/dx**2. &
                      -B_param*g*HH(i)*(detax(i+1)-2.d0*detax(i)+detax(i-1))/dx**2 &
                        -HHH(i)/6.d0*(s1_H(i+1)-2.d0*s1_H(i) + s1_H(i-1))/dx**2                
            endif
        endif
                
    enddo

    !! if a wave packet intercept the slope, then use the SWE

    do i = 2,mx-1
        if (q(1,i) > 0.01d0 .and. &
               (eta(i)-eta(i-1))*(eta(i+1)-eta(i)) <  0.d0 ) then
                ! local min or max

            iL = i - abs(eta(i))*1./dx
            iL = max(1,iL)
            iR = i + abs(eta(i))*1./dx
            iR = min(mx,iR)

            if (maxval(slope(iL:iR))>0.d0) then
                psi(iL:iR) = 0.d0
            endif

        endif

    end do

        !! Limit the dispersion for large troughs
    do i=2,mx-1

        if (-eta(i)>.5*q(1,i)) then

            iL = i - abs(eta(i))*2/dx
            iL = max(1,iL)
            iR = i + abs(eta(i))*2/dx
            iR = min(mx,iR)

            !psi(iL:iR) = 0.d0

        endif

    end do

    if (.not.use_bous_sw_thresh) go to 994

    do i = 2,mx-1
        k= i+1

        topo = -slope(i)

        if (q(1,i)>1.d-1 .and. topo>1.d-1.and. &
            q(1,i)-topo>bous_sw_thresh*topo) then
 
            psi(k) = 0.d0

            xaxis1: do kk=i,1,-1
                if (eta(kk)-eta(kk-1) >  0.d0 ) then
                    !if ((q(1,kk,j)+aux(1,kk,j))> sea_level ) then
                    psi(kk) = 0.d0
                else
                    exit xaxis1
                endif
            end do xaxis1

            xaxis2: do kk=i,mx
                if (eta(kk+1)-eta(kk)< 0.d0 ) then
                    !if ((q(1,kk,j)+aux(1,kk,j))> sea_level ) then
                    psi(kk) = 0.d0
                else
                    exit xaxis2
                endif
            end do xaxis2

        elseif (q(1,i)>1.d-1 .and. topo > 1.d-1.and. &
            q(1,i)-topo < -bous_sw_thresh*topo) then
 
            psi(k) = 0.d0

            xaxis11: do kk=i,1,-1
                if (eta(kk)-eta(kk-1)< 0.d0 ) then
                    psi(kk) = 0.d0
                else
                    exit xaxis11
                endif
            end do xaxis11

            xaxis12: do kk=i,mx
                if (eta(kk+1)-eta(kk)> 0.d0 ) then
                    psi(kk) = 0.d0
                else
                    exit xaxis12
                endif
            end do xaxis12
        endif
    enddo

994 continue

    return
end subroutine read_psi
