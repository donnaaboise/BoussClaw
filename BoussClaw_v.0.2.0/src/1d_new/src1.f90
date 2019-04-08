subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
    use geoclaw_module, only: rad2deg
    use bous_module

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,maux
    double precision, intent(in) :: xlower,dx,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc)

    ! Locals
    integer :: i, nman
    real(kind=8) :: h, hu, gamma, dgamma, fdt, a(2,2), coeff

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    logical :: verbose
    real(kind=8) :: wave_max_slope,eta(1:mx),eta_max

    ! Friction source term
    if (friction_forcing) then
        do i=1,mx
            ! Extract appropriate momentum
            if (q(1,i) < depth_tolerance) then
                q(2,i) = 0.d0
            else
                !! Apply friction source term only if in shallower water
                if (q(1,i) <= friction_depth) then
                    do nman = num_manning, 1, -1
                        if (aux(1,i) .lt. manning_break(nman)) then
                            coeff = manning_coefficient(nman)
                        endif
                    enddo
                    !! Calculate source term
                    gamma = sqrt(q(2,i)**2) * g     &   
                       * coeff**2 / (q(1,i)**(7.d0/3.d0))
                    dgamma = 1.d0 + dt * gamma
                    q(2, i) = q(2, i) / dgamma
                endif
            endif
        enddo
    endif
    !! End of friction source term
          
    !! Computation of Dispersive Terms
     
    if (bouss) then
         call src1_bous(meqn,mbc,mx, &
                  xlower,dx,q,maux,aux,t,dt)
    endif

    verbose = .True.

    if (verbose) then
        wave_max_slope = 0.d0
        eta = 0.d0
        do i = 1,mx
            if (q(1,i) > 0.d0) then
                wave_max_slope=max(wave_max_slope, &
                        abs( (q(1,i) + aux(1,i) - (q(1,i-1) + aux(1,i-1)) )/dx/1.d0))
               eta(i) = q(1,i)+aux(1,i)
            endif
         end do
         eta_max=maxval(eta)/dx
         !wave_max_slope = wave_max_slope*rad2deg
         !print*,'maximum wave slope =', wave_max_slope
         !print*,'max possible slope =', eta_max/dx
         !print*,'ratio              =', wave_max_slope/eta_max,eta_max*dx,dx
    endif

    return
end subroutine src1
      
