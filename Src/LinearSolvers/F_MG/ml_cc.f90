module ml_cc_module

  use bl_constants_module
  use mg_module
  use ml_layout_module
  use bndry_reg_module
! use mg_hypre_module

  implicit none

  private :: ml_fill_fluxes, ml_fill_n_fluxes
  private :: scale_residual_1d, scale_residual_2d, scale_residual_3d

contains

  subroutine ml_cc(mla, mgt, rh, full_soln, fine_mask, ref_ratio, &
                   do_diagnostics, rel_eps, abs_eps_in, need_grad_phi_in, final_resnorm)

    use bl_prof_module
    use ml_util_module        , only : ml_norm_inf
    use ml_restriction_module , only : ml_restriction
    use ml_prolongation_module, only : ml_prolongation, ml_interp_bcs

    type(ml_layout), intent(in   ) :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: rh(:)
    type( multifab), intent(inout) :: full_soln(:)
    type(lmultifab), intent(in   ) :: fine_mask(:)
    integer        , intent(in   ) :: ref_ratio(:,:)
    integer        , intent(in   ) :: do_diagnostics
    real(dp_t)     , intent(in   ) :: rel_eps

    real(dp_t)     , intent(in   ), optional :: abs_eps_in
    logical        , intent(in   ), optional :: need_grad_phi_in
    real(dp_t)     , intent(  out), optional :: final_resnorm

    integer :: nlevs
    type(multifab), allocatable  ::      soln(:)
    type(multifab), allocatable  ::        uu(:)
    type(multifab), allocatable  ::   uu_hold(:)
    type(multifab), allocatable  ::       res(:)
    type(multifab), allocatable  ::  temp_res(:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box) :: pd, pdc
    type(layout) :: la, lac
    integer :: i, n, dm
    integer :: mglev, mglev_crse, iter
    logical :: fine_converged,need_grad_phi,lcross

    real(dp_t) :: Anorm, bnorm, abs_eps, ni_res
    real(dp_t) :: tres, tres0, max_norm
    real(dp_t) :: rho, sum

    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_cc")

    dm = get_dim(rh(1))
    nlevs = mla%nlevel

    if ( present(abs_eps_in) ) then
       abs_eps = abs_eps_in 
    else
       abs_eps = mgt(nlevs)%abs_eps
    end if

    if ( present(need_grad_phi_in) ) then
       need_grad_phi = need_grad_phi_in 
    else
       need_grad_phi = .false.
    end if

    allocate(soln(nlevs), uu(nlevs), res(nlevs), temp_res(nlevs))
    allocate(uu_hold(2:nlevs-1))
    allocate(brs_flx(2:nlevs))
    allocate(brs_bcs(2:nlevs))

    do n = 2,nlevs-1
       la = mla%la(n)
       call build(uu_hold(n),la,1,1)
       call setval( uu_hold(n), ZERO,all=.true.)
    end do
    
    do n = nlevs, 1, -1

       la = mla%la(n)
       call build(    soln(n), la, 1, nghost(full_soln(1)))
       call build(      uu(n), la, 1, nghost(full_soln(1)))
       call build(     res(n), la, 1, 0)
       call build(temp_res(n), la, 1, 0)
       call setval(    soln(n), ZERO,all=.true.)
       call setval(      uu(n), ZERO,all=.true.)
       call setval(     res(n), ZERO,all=.true.)
       call setval(temp_res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, &
            width = 0)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, &
            width = 2, other = .false.)

    end do

    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels
       call ml_restriction(rh(n-1), rh(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
    end do
    bnorm = ml_norm_inf(rh,fine_mask)

    lcross = ((ncomp(mgt(nlevs)%ss(mgt(nlevs)%nlevels)) == 5) .or. &
         (ncomp(mgt(nlevs)%ss(mgt(nlevs)%nlevels)) == 7))

    Anorm = stencil_norm(mgt(nlevs)%ss(mgt(nlevs)%nlevels))
    do n = 1, nlevs-1
       Anorm = max(stencil_norm(mgt(n)%ss(mgt(n)%nlevels),fine_mask(n)),Anorm)
    end do

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in 
    !  its ghost cells before we evaluate the initial residual  
    do n = 2,nlevs
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, &
                             ref_ratio(n-1,:), -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, &
                             ref_ratio(n-1,:), +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do

    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), &
                      mgt(n)%mm(mglev))
    end do

    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
                                  ref_ratio(n-1,:))

       call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
    enddo

    do n = 1,nlevs
       call multifab_copy(rh(n),res(n),ng = nghost(rh(n)))
    end do

    tres0 = ml_norm_inf(rh,fine_mask)
    if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) then
       write(unit=*, &
             fmt='("F90mg: Initial rhs                  = ",g15.8)') bnorm
       write(unit=*, &
             fmt='("F90mg: Initial residual (resid0)    = ",g15.8)') tres0
    end if

    ! ************************************************************************
    !  Define norm to be used for convergence testing that is the maximum
    !    of bnorm (norm of rhs) and tres0 (norm of resid0)
    ! ************************************************************************
   
    max_norm = max(bnorm,tres0)

    fine_converged = .false.

    if ( ml_converged(res, soln, fine_mask, bnorm, Anorm, rel_eps, abs_eps, ni_res, mgt(nlevs)%verbose) ) then
       if ( parallel_IOProcessor() .and. mgt(nlevs)%verbose > 0 ) &
            write(unit=*, fmt='("F90mg: No iterations needed ")')

!   else if (mgt%use_hypre .eq. 1) then

!      if (nlevs .gt. 1) then
!         call bl_error("ml_cc: can't use HYPRE with nlevs > 1")
!      else
!         call cc_hypre_solve(mgt, res, soln, eps)
!      end if

    else 

     do iter = 1, mgt(nlevs)%max_iter

       if ( fine_converged ) then
          if ( ml_converged(res, soln, fine_mask, max_norm, Anorm, rel_eps, abs_eps, ni_res, mgt(nlevs)%verbose) ) exit
       end if

       ! Set: uu = 0
       do n = 1,nlevs
          call setval(uu(n), ZERO, all=.true.)
       end do

       ! Set: uu_hold = 0
       do n = 2,nlevs-1
          call setval(uu_hold(n), ZERO, all=.true.)
       end do

       !   Down the V-cycle
       do n = nlevs,1,-1

          mglev = mgt(n)%nlevels

          ! Enforce solvability if appropriate
          if (nlevs .eq. 1 .and. mgt(1)%bottom_singular) then

             sum = multifab_sum(res(1))  / boxarray_dvolume(get_boxarray(res(1)))

             ! Set this to all one for use in saxpy 
             call setval( uu(1),  ONE, all=.true.)

             ! Subtract "sum" from res(1) in order to make this solvable
             call  saxpy(res(1), -sum, uu(1))
  
             ! Return this to zero
             call setval( uu(1), ZERO, all=.true.)

             if ( parallel_IOProcessor() .and. (do_diagnostics == 1) ) then
                write(unit=*, fmt='("F90mg: Subtracting from rhs ",g15.8)') sum
             end if

          end if

          if ( do_diagnostics == 1 ) then
             tres = norm_inf(res(n))
             if ( parallel_ioprocessor() ) then
                print *,'DWN: RES BEFORE GSRB AT LEVEL ',n, tres
             end if
          end if

          ! Relax ...
          if (iter < mgt(nlevs)%max_iter) then
             if (n > 1) then
                call mini_cycle(mgt(n), mglev, &
                                mgt(n)%ss(mglev), uu(n), res(n), &
                                mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2)
             else 
               call mg_tower_cycle(mgt(n), mgt(n)%cycle_type, mglev, &
                                   mgt(n)%ss(mglev), uu(n), res(n), &
                                   mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
                                   mgt(n)%gamma)
             end if
          end if

          ! Add: Soln += uu
          call saxpy(soln(n), ONE, uu(n))

          if (n > 1) then

             mglev_crse = mgt(n-1)%nlevels

             ! Compute COARSE Res = Rh - Lap(Soln)
             call mg_defect(mgt(n-1)%ss(mglev_crse),res(n-1), &
                  rh(n-1),soln(n-1),mgt(n-1)%mm(mglev_crse))

             ! Compute FINE Res = Res - Lap(uu)
             call mg_defect(mgt(n)%ss(mglev),temp_res(n), &
                  res(n),uu(n),mgt(n)%mm(mglev))
             call multifab_copy(res(n), temp_res(n), ng = nghost(res(n)))

             if (do_diagnostics == 1 ) then
                tres = norm_inf(res(n))
                if ( parallel_ioprocessor()) then
                   print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, tres
                end if
             end if

             ! Compute CRSE-FINE Res = Res - Crse Flux(soln) + Fine Flux(soln)
             pdc = layout_get_pd(mla%la(n-1))
             call crse_fine_residual_cc(n,mgt,soln,res(n-1),brs_flx(n), &
                                        pdc,ref_ratio(n-1,:))

             ! Restrict FINE Res to COARSE Res (important to do this last
             !     so we overwrite anything extra which may have been defined
             !     above near fine-fine interfaces)
             call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev), &
                  mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))

             ! Copy u_hold = uu
             if (n < nlevs) call multifab_copy(uu_hold(n), uu(n), ng = nghost(uu(n)))

             ! Set: uu = 0
             call setval(uu(n), ZERO, all=.true.)

          else

             if (do_diagnostics == 1 ) then
                call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), &
                               mgt(n)%mm(mglev))
                tres = norm_inf(temp_res(n))
                if ( parallel_ioprocessor()) then
                   print *,'DWN: RES AFTER  GSRB AT LEVEL ',n, tres
                end if
             else
                ! Seem to need this in periodic case to get right answer...
                call multifab_fill_boundary(uu(n), cross = lcross)
             end if

          end if

       end do

       !   Back up the V-cycle
       do n = 2, nlevs

          pd = layout_get_pd(mla%la(n))
          mglev = mgt(n)%nlevels

          ! Interpolate uu from coarser level
          call ml_prolongation(uu(n), uu(n-1), ref_ratio(n-1,:))

          ! Add: soln += uu
          call saxpy(soln(n), ONE, uu(n), .true.)

          ! Add: uu_hold += uu
          if (n < nlevs) call saxpy(uu_hold(n), ONE, uu(n), .true.)

          ! Interpolate uu to supply boundary conditions for new 
          ! residual calculation
          call bndry_reg_copy(brs_bcs(n), uu(n-1))
          do i = 1, dm
             call ml_interp_bcs(uu(n), brs_bcs(n)%bmf(i,0), pd, &
                                ref_ratio(n-1,:), -i)
             call ml_interp_bcs(uu(n), brs_bcs(n)%bmf(i,1), pd, &
                                ref_ratio(n-1,:), +i)
          end do
          call multifab_fill_boundary(uu(n))

          ! Compute Res = Res - Lap(uu)
          call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), &
                         mgt(n)%mm(mglev))
          call multifab_copy(res(n), temp_res(n), ng = nghost(res(n)))

          if (do_diagnostics == 1 ) then
             tres = norm_inf(temp_res(n))
             if ( parallel_ioprocessor() ) then
                print *,'UP : RES BEFORE GSRB AT LEVEL ',n, tres
             end if
          end if

          ! Set: uu = 0
          call setval(uu(n), ZERO, all=.true.)

          ! Relax ...
          call mini_cycle(mgt(n), mglev, mgt(n)%ss(mglev), &
               uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2)

          ! Compute Res = Res - Lap(uu)

          if (do_diagnostics == 1 ) then
             call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), &
                            mgt(n)%mm(mglev))
             call multifab_copy(res(n), temp_res(n), ng = nghost(res(n)))
             tres = norm_inf(res(n))
             if ( parallel_ioprocessor() ) then
                print *,'UP : RES AFTER  GSRB AT LEVEL ',n, tres
                if (n == nlevs) print *,' '
             end if
          end if

          ! Add: soln += uu
          call saxpy(soln(n), ONE, uu(n), .true.)

          ! Add: uu += uu_hold so that it will be interpolated too
          if (n < nlevs) call saxpy(uu(n), ONE, uu_hold(n), .true.)

          ! Only do this as long as tangential interp looks under fine grids
          mglev_crse = mgt(n-1)%nlevels
          call ml_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
               mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))

       end do

       !    Average the solution to the coarser grids.
       do n = nlevs,2,-1
          mglev      = mgt(n)%nlevels
          mglev_crse = mgt(n-1)%nlevels
          call ml_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
               mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
       end do

       do n = 1,nlevs
          call multifab_fill_boundary(soln(n), cross = lcross)
       end do

       ! Interpolate soln to supply boundary conditions 
       do n = 2,nlevs
          pd = layout_get_pd(mla%la(n))
          call bndry_reg_copy(brs_bcs(n), soln(n-1))
          do i = 1, dm
             call ml_interp_bcs(soln(n), brs_bcs(n)%bmf(i,0), pd, &
                                ref_ratio(n-1,:), -i)
             call ml_interp_bcs(soln(n), brs_bcs(n)%bmf(i,1), pd, &
                                ref_ratio(n-1,:), +i)
          end do
          call multifab_fill_boundary(soln(n))
       end do

       !    Optimization so don't have to do multilevel convergence test 
       !    each time

       !    Compute the residual on just the finest level
       n = nlevs
       mglev = mgt(n)%nlevels
       call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n),mgt(n)%mm(mglev))

       if ( ml_fine_converged(res, soln, max_norm, Anorm, rel_eps, abs_eps) ) then

          fine_converged = .true.

          !      Compute the residual on every level
          do n = 1,nlevs-1
             mglev = mgt(n)%nlevels
             call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),soln(n), &
                            mgt(n)%mm(mglev))
          end do

          !      Compute the coarse-fine residual 
          do n = nlevs,2,-1
             pdc = layout_get_pd(mla%la(n-1))
             call crse_fine_residual_cc(n,mgt,soln,res(n-1),brs_flx(n),pdc, &
                                        ref_ratio(n-1,:))
          end do

          !      Average the fine residual onto the coarser level
          do n = nlevs,2,-1
             mglev      = mgt(n  )%nlevels
             mglev_crse = mgt(n-1)%nlevels
             call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
                  mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
          end do

          if ( mgt(nlevs)%verbose > 1 ) then
             do n = 1,nlevs
                tres = norm_inf(res(n))
                if ( parallel_ioprocessor() ) then
!                  write(unit=*, fmt='(i3,": Level ",i2,"  : SL_Ninf(defect) = ",g15.8)') &
!                       iter,n,tres
                   write(unit=*, fmt='("F90mg: Iteration   ",i3," Lev ",i1," resid/resid0 = ",g15.8)') &
                        iter,n,tres/tres0
                end if
             end do
!            tres = ml_norm_inf(res,fine_mask)
!            if ( parallel_ioprocessor() ) then
!               write(unit=*, fmt='(i3,": All Levels: ML_Ninf(defect) = ",g15.8)') iter, tres
!            end if
          end if

       else 

          fine_converged = .false.
          if ( mgt(nlevs)%verbose > 1 ) then 
             tres = norm_inf(res(nlevs))
             if ( parallel_IOProcessor() ) then
!               write(unit=*, fmt='(i3,": FINE_Ninf(defect) = ",g15.8)') iter, tres
                write(unit=*, fmt='("F90mg: Iteration   ",i3," Fine  resid/resid0 = ",g15.8)') iter,tres/tres0
             end if
          end if

       end if


     end do

    ! ****************************************************************************

     iter = iter-1
     if (iter < mgt(nlevs)%max_iter) then
        if ( mgt(nlevs)%verbose > 0 ) then
          tres = ml_norm_inf(res,fine_mask)
          if ( parallel_IOProcessor() ) then
             if (tres0 .gt. 0.0_dp_t) then
               write(unit=*, fmt='("F90mg: Final Iter. ",i3," resid/resid0 = ",g15.8)') iter,tres/tres0
               write(unit=*, fmt='("")') 
             else
               write(unit=*, fmt='("F90mg: Final Iter. ",i3," resid/resid0 = ",g15.8)') iter,0.0_dp_t
               write(unit=*, fmt='("")') 
             end if
          end if
        end if
     else
        call bl_error("Multigrid Solve: failed to converge in max_iter iterations")
     end if

    end if

    ! Add: soln += full_soln
    do n = 1,nlevs
       call saxpy(full_soln(n),ONE,soln(n))
       call multifab_fill_boundary(full_soln(n))
    end do

    do n = 2,nlevs-1
       call multifab_destroy(uu_hold(n))
    end do

    if (need_grad_phi) then

       !   Interpolate boundary conditions of soln in order to get correct grad(phi) at
       !   crse-fine boundaries
       do n = 2,nlevs
          pd = layout_get_pd(mla%la(n))
          call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
          do i = 1, dm
             call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, ref_ratio(n-1,:), -i)
             call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, ref_ratio(n-1,:), +i)
          end do
          call multifab_fill_boundary(full_soln(n))
       end do

    end if

    !   Make sure all periodic and internal boundaries are filled
    do n = 1,nlevs
       call multifab_fill_boundary(full_soln(n))
    end do

    do n = nlevs, 1, -1
       call multifab_destroy(    soln(n))
       call multifab_destroy(      uu(n))
       call multifab_destroy(     res(n))
       call multifab_destroy(temp_res(n))
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_flx(n))
       call bndry_reg_destroy(brs_bcs(n))
    end do

    call destroy(bpt)

    if ( present(final_resnorm) ) &
       final_resnorm = ni_res

  contains

    function ml_fine_converged(res, sol, bnorm, Anorm, rel_eps, abs_eps) result(r)
      logical :: r
      type(multifab), intent(in) :: res(:), sol(:)
      real(dp_t), intent(in) :: Anorm, rel_eps, abs_eps, bnorm
      real(dp_t) :: ni_res, ni_sol
      integer    :: nlevs
      nlevs = size(res)
      ni_res = norm_inf(res(nlevs))
      ni_sol = norm_inf(sol(nlevs))
!      r =  ni_res <= rel_eps*(Anorm*ni_sol + bnorm) .or. &
!           ni_res <= abs_eps .or. &
!           ni_res <= epsilon(Anorm)*Anorm
      r =  ni_res <= rel_eps*(bnorm) .or. &
           ni_res <= abs_eps
    end function ml_fine_converged

    function ml_converged(res, sol, mask, bnorm, Anorm, rel_eps, abs_eps, ni_res, verbose) result(r)

      use ml_util_module, only : ml_norm_inf

      logical :: r
      integer :: verbose
      type(multifab), intent(in) :: res(:), sol(:)
      type(lmultifab), intent(in) :: mask(:)
      real(dp_t), intent(in   ) :: Anorm, rel_eps, abs_eps, bnorm
      real(dp_t), intent(  out) :: ni_res
      real(dp_t) :: ni_sol

      ni_res = ml_norm_inf(res, mask)
      ni_sol = ml_norm_inf(sol, mask)
!      r =  ni_res <= rel_eps*(Anorm*ni_sol + bnorm) .or. &
!           ni_res <= abs_eps .or. &
!           ni_res <= epsilon(Anorm)*Anorm
      r =  ni_res <= rel_eps*(bnorm) .or. &
           ni_res <= abs_eps 
      if ( r .and. parallel_IOProcessor() .and. verbose > 1) then
         if (ni_res <= rel_eps*Anorm*ni_sol) then
            print *,'Converged res < rel_eps*Anorm*sol'
         else if (ni_res <= abs_eps) then
            print *,'Converged res < abs_eps '
         else if (ni_res <= rel_eps*bnorm) then
            print *,'Converged res < rel_eps*bnorm '
         else 
            print *,'Converged res < epsilon(Anorm)*Anorm'
         end if
      end if
    end function ml_converged

  end subroutine ml_cc

!
! ******************************************************************************************
!

  subroutine crse_fine_residual_cc(n, mgt, uu, crse_res, brs_flx, pdc, ref_ratio)

      use ml_interface_stencil_module, only : ml_interface

      integer        , intent(in   ) :: n
      type(mg_tower) , intent(inout) :: mgt(:)
      type(bndry_reg), intent(inout) :: brs_flx
      type(multifab) , intent(inout) :: uu(:)
      type(multifab) , intent(inout) :: crse_res
      type(box)      , intent(in   ) :: pdc
      integer        , intent(in   ) :: ref_ratio(:)

      integer :: i, dm, mglev

      dm = brs_flx%dim
      mglev = mgt(n)%nlevels

      call multifab_fill_boundary(uu(n))

      do i = 1, dm
         call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,0), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), -1, i)
         call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,1), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), 1, i)
      end do
      call bndry_reg_copy_to_other(brs_flx)
      do i = 1, dm
         call ml_interface(crse_res, brs_flx%obmf(i,0), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, -1, i, ONE)
         call ml_interface(crse_res, brs_flx%obmf(i,1), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, +1, i, ONE)
      end do

  end subroutine crse_fine_residual_cc

 subroutine crse_fine_residual_n_cc(n, mgt, uu, crse_res, brs_flx, pdc, ref_ratio)

      use ml_interface_stencil_module, only : ml_interface

      integer        , intent(in   ) :: n
      type(mg_tower) , intent(inout) :: mgt(:)
      type(bndry_reg), intent(inout) :: brs_flx
      type(multifab) , intent(inout) :: uu(:)
      type(multifab) , intent(inout) :: crse_res
      type(box)      , intent(in   ) :: pdc
      integer        , intent(in   ) :: ref_ratio(:)

      integer :: i, dm, mglev

      dm = brs_flx%dim
      mglev = mgt(n)%nlevels

      call multifab_fill_boundary(uu(n))

      do i = 1, dm
         call ml_fill_n_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,0), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), -1, i)
         call ml_fill_n_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,1), &
              uu(n), mgt(n)%mm(mglev), ref_ratio(i), 1, i)
      end do
      call bndry_reg_copy_to_other(brs_flx)
      do i = 1, dm
         call ml_interface(crse_res, brs_flx%obmf(i,0), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, -1, i, ONE)
         call ml_interface(crse_res, brs_flx%obmf(i,1), uu(n-1), &
              mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, +1, i, ONE)
      end do

  end subroutine crse_fine_residual_n_cc

!
! ******************************************************************************************
!

  subroutine ml_resid(mla, mgt, rh, res, full_soln, ref_ratio)

    use ml_restriction_module  , only : ml_restriction
    use ml_prolongation_module , only : ml_interp_bcs

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: rh(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box)    :: pd, pdc
    type(layout) :: la, lac
    integer      :: i, n, dm, nlevs, mglev, mglev_crse

    dm = get_dim(rh(1))

    nlevs = mla%nlevel

    allocate(brs_bcs(2:nlevs))
    allocate(brs_flx(2:nlevs))

    do n = nlevs, 1, -1

       la = mla%la(n)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, width = 0)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2, other = .false.)

    end do

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in its ghost cells 
    !   before we evaluate the initial residual  
    do n = 2,nlevs
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, ref_ratio(n-1,:), -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, ref_ratio(n-1,:), +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do

    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n),mgt(n)%mm(mglev))
    end do

    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc,ref_ratio(n-1,:))

       call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
    enddo

    do n = nlevs, 1, -1
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_flx(n))
       call bndry_reg_destroy(brs_bcs(n))
    end do

  end subroutine ml_resid

!
! ******************************************************************************************
!

  subroutine ml_cc_applyop(mla, mgt, res, full_soln, ref_ratio)

    use bl_prof_module
    use ml_restriction_module , only : ml_restriction
    use ml_prolongation_module, only : ml_prolongation, ml_interp_bcs

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    integer :: nlevs
    type(multifab), allocatable  ::      soln(:)
    type(multifab), allocatable  ::        uu(:)
    type(multifab), allocatable  ::   uu_hold(:)
    type(multifab), allocatable  ::        rh(:) ! this will be set to zero
    type(multifab), allocatable  ::  temp_res(:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box) :: pd, pdc
    type(layout) :: la, lac
    integer :: i, n, dm
    integer :: mglev, mglev_crse

    type(bl_prof_timer), save :: bpt
    integer                   :: lo(get_dim(res(1))),hi(get_dim(res(1))),ng
    real(kind=dp_t),  pointer :: resp(:,:,:,:)

    call build(bpt, "ml_cc_applyop")

    nlevs = mla%nlevel

    allocate(soln(nlevs), uu(nlevs), rh(nlevs), temp_res(nlevs))
    allocate(uu_hold(2:nlevs-1))

    allocate(brs_flx(2:nlevs))
    allocate(brs_bcs(2:nlevs))

    do n = 2,nlevs-1
       la = mla%la(n)
       call build(uu_hold(n),la,1,1)
       call setval( uu_hold(n), ZERO,all=.true.)
    end do

    do n = nlevs, 1, -1

       la = mla%la(n)
       call build(    soln(n), la, 1, 1)
       call build(      uu(n), la, 1, 1)
       call build(      rh(n), la, 1, 0)
       call build(temp_res(n), la, 1, 0)
       call setval(    soln(n), ZERO,all=.true.)
       call setval(      uu(n), ZERO,all=.true.)
       call setval(      rh(n), ZERO,all=.true.)
       call setval(temp_res(n), ZERO,all=.true.)

       ! zero residual just to be safe
       call setval(     res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2, other = .false.)
       call bndry_reg_rr_build(brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, width = 0)

    end do

    dm = get_dim(rh(1))

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in 
    !  its ghost cells before we evaluate the initial residual  
    do n = 2,nlevs
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, &
                             ref_ratio(n-1,:), -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, &
                             ref_ratio(n-1,:), +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do


    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), &
                      mgt(n)%mm(mglev))
    end do

    ! Make sure to correct the coarse cells immediately next to fine grids
    !   using the averaged fine grid fluxes
    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
                                  ref_ratio(n-1,:))

       call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
    enddo


    ! still need to multiply residual by -1 to get (alpha - del dot beta grad)
    do n=1,nlevs
       ng = nghost(res(n))
       
       do i=1, nboxes(res(n))
          if (multifab_remote(res(n),i)) cycle
          resp  => dataptr(res(n),i)
          lo =  lwb(get_box(res(n), i))
          hi =  upb(get_box(res(n), i))
          select case (dm)
          case (1)
             call scale_residual_1d(lo,hi,ng,resp(:,1,1,1))
          case (2)
             call scale_residual_2d(lo,hi,ng,resp(:,:,1,1))
          case (3)
             call scale_residual_3d(lo,hi,ng,resp(:,:,:,1))
          end select
       end do
    enddo

    do n = 2,nlevs-1
       call destroy(uu_hold(n))
    end do

    do n = nlevs, 1, -1
       call destroy(soln(n))
       call destroy(uu(n))
       call destroy(rh(n))
       call destroy(temp_res(n))
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_bcs(n))
       call bndry_reg_destroy(brs_flx(n))
    end do

    call destroy(bpt)

  end subroutine ml_cc_applyop

 subroutine ml_cc_n_applyop(mla, mgt, res, full_soln, ref_ratio)

    use bl_prof_module
    use ml_restriction_module , only : ml_restriction
    use ml_prolongation_module, only : ml_prolongation, ml_interp_bcs

    type(ml_layout), intent(in)    :: mla
    type(mg_tower) , intent(inout) :: mgt(:)
    type( multifab), intent(inout) :: res(:)
    type( multifab), intent(inout) :: full_soln(:)
    integer        , intent(in   ) :: ref_ratio(:,:)

    integer :: nlevs
    type(multifab), allocatable  ::      soln(:)
    type(multifab), allocatable  ::        uu(:)
    type(multifab), allocatable  ::   uu_hold(:)
    type(multifab), allocatable  ::        rh(:) ! this will be set to zero
    type(multifab), allocatable  ::  temp_res(:)

    type(bndry_reg), allocatable :: brs_flx(:)
    type(bndry_reg), allocatable :: brs_bcs(:)

    type(box) :: pd, pdc
    type(layout) :: la, lac
    integer :: i, n, dm, nComp
    integer :: mglev, mglev_crse

    type(bl_prof_timer), save :: bpt
    integer                   :: lo(get_dim(res(1))),hi(get_dim(res(1))),ng
    real(kind=dp_t),  pointer :: resp(:,:,:,:)

    call build(bpt, "ml_cc_applyop")

    nlevs = mla%nlevel

    allocate(soln(nlevs), uu(nlevs), rh(nlevs), temp_res(nlevs))
    allocate(uu_hold(2:nlevs-1))

    allocate(brs_flx(2:nlevs))
    allocate(brs_bcs(2:nlevs))

    do n = 2,nlevs-1
       la = mla%la(n)
       call build(uu_hold(n),la,1,1)
       call setval( uu_hold(n), ZERO,all=.true.)
    end do

    dm    = 2
    nComp = 2

    do n = nlevs, 1, -1

       la = mla%la(n)
       call build(    soln(n), la, 1, 1)
       call build(      uu(n), la, 1, 1)
       call build(      rh(n), la, 1, 0)
       call build(temp_res(n), la, 1, 0)
       call setval(    soln(n), ZERO,all=.true.)
       call setval(      uu(n), ZERO,all=.true.)
       call setval(      rh(n), ZERO,all=.true.)
       call setval(temp_res(n), ZERO,all=.true.)

       ! zero residual just to be safe
       call setval(     res(n), ZERO,all=.true.)

       if ( n == 1 ) exit

       ! Build the (coarse resolution) flux registers to be used in computing
       !  the residual at a non-finest AMR level.

       pdc = layout_get_pd(mla%la(n-1))
       lac = mla%la(n-1)
       call bndry_reg_rr_build(brs_bcs(n), la, lac, ref_ratio(n-1,:), pdc, width = 2, other = .false.)
       call bndry_reg_rr_build(brs_flx(n), la, lac, ref_ratio(n-1,:), pdc, nc = nComp, width = 0)

    end do

    dm = get_dim(rh(1))

    !  Make sure full_soln at fine grid has the correct coarse grid bc's in 
    !  its ghost cells before we evaluate the initial residual  
    do n = 2,nlevs
       pd = layout_get_pd(mla%la(n))
       call bndry_reg_copy(brs_bcs(n), full_soln(n-1))
       do i = 1, dm
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,0), pd, &
                             ref_ratio(n-1,:), -i)
          call ml_interp_bcs(full_soln(n), brs_bcs(n)%bmf(i,1), pd, &
                             ref_ratio(n-1,:), +i)
       end do
       call multifab_fill_boundary(full_soln(n))
    end do

    !   Make sure all periodic and internal boundaries are filled as well
    do n = 1,nlevs   
       call multifab_fill_boundary(full_soln(n))
    end do


    do n = 1,nlevs,1
       mglev = mgt(n)%nlevels
       call mg_defect(mgt(n)%ss(mglev),res(n),rh(n),full_soln(n), &
                      mgt(n)%mm(mglev))
    end do

    ! Make sure to correct the coarse cells immediately next to fine grids
    !   using the averaged fine grid fluxes
    do n = nlevs,2,-1
       mglev      = mgt(n  )%nlevels
       mglev_crse = mgt(n-1)%nlevels

       pdc = layout_get_pd(mla%la(n-1))
       call crse_fine_residual_n_cc(n,mgt,full_soln,res(n-1),brs_flx(n),pdc, &
                                    ref_ratio(n-1,:))
       call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),&
            mgt(n-1)%mm(mglev_crse), ref_ratio(n-1,:))
    enddo


    ! still need to multiply residual by -1 to get (alpha - del dot beta grad)
    do n=1,nlevs
       ng = nghost(res(n))
       
       do i=1, nboxes(res(n))
          if (multifab_remote(res(n),i)) cycle
          resp  => dataptr(res(n),i)
          lo =  lwb(get_box(res(n), i))
          hi =  upb(get_box(res(n), i))
          select case (dm)
          case (1)
             call scale_residual_1d(lo,hi,ng,resp(:,1,1,1))
          case (2)
             call scale_residual_2d(lo,hi,ng,resp(:,:,1,1))
          case (3)
             call scale_residual_3d(lo,hi,ng,resp(:,:,:,1))
          end select
       end do
    enddo

    do n = 2,nlevs-1
       call destroy(uu_hold(n))
    end do

    do n = nlevs, 1, -1
       call destroy(soln(n))
       call destroy(uu(n))
       call destroy(rh(n))
       call destroy(temp_res(n))
       if ( n == 1 ) exit
       call bndry_reg_destroy(brs_bcs(n))
       call bndry_reg_destroy(brs_flx(n))
    end do

    call destroy(bpt)

  end subroutine ml_cc_n_applyop

!
! ******************************************************************************************
!
  subroutine ml_fill_fluxes(ss, flux, uu, mm, ratio, face, dim)
    use bl_prof_module
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    integer :: ratio
    integer :: face, dim
    integer :: i, n, dm
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: ng
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fill_fluxes")

    ng = nghost(uu)

    if ( ncomp(uu) /= ncomp(flux) ) then
       call bl_error("ML_FILL_FLUXES: uu%nc /= flux%nc")
    end if

    dm = get_dim(ss)

    do i = 1, nboxes(flux)
       if ( remote(flux, i) ) cycle
       fp => dataptr(flux, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, ncomp(uu)
          select case(dm)
          case (1)
             call stencil_flux_1d(sp(:,1,1,:), fp(:,1,1,n), up(:,1,1,n), &
                  mp(:,1,1,1), ng, ratio, face, dim)
          case (2)
             call stencil_flux_2d(sp(:,:,1,:), fp(:,:,1,n), up(:,:,1,n), &
                  mp(:,:,1,1), ng, ratio, face, dim)
          case (3)
             call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
                  mp(:,:,:,1), ng, ratio, face, dim)
          end select
       end do
    end do
    call destroy(bpt)
  end subroutine ml_fill_fluxes

!
! ******************************************************************************************
!

  subroutine ml_fill_n_fluxes(ss, flux, uu, mm, ratio, face, dim)
    use bl_prof_module
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    integer :: ratio
    integer :: face, dim
    integer :: i, n
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: ng
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fill_fluxes")

    ng = nghost(uu)

    do i = 1, nboxes(flux)
       if ( remote(flux, i) ) cycle
       fp => dataptr(flux, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, ncomp(uu)
          select case(get_dim(ss))
          case (1)
             call stencil_flux_1d(sp(:,1,1,:), fp(:,1,1,n), up(:,1,1,n), &
                  mp(:,1,1,1), ng, ratio, face, dim)
          case (2)
             if ( ncomp(flux) > 1 ) then
                call stencil_flux_n_2d(sp(:,:,1,:), fp(:,:,1,:), up(:,:,1,n), &
                     mp(:,:,1,1), ng, ratio, face, dim)
             else
                call stencil_flux_2d(sp(:,:,1,:), fp(:,:,1,n), up(:,:,1,n), &
                     mp(:,:,1,1), ng, ratio, face, dim)
             end if
          case (3)
             call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
                  mp(:,:,:,1), ng, ratio, face, dim)
          end select
       end do
    end do
    call destroy(bpt)
  end subroutine ml_fill_n_fluxes

!
! ******************************************************************************************
!

! Multiply residual by -1 in 1d
  subroutine scale_residual_1d(lo,hi,ng,res)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(inout) :: res(lo(1)-ng:)

! Local
  integer :: i

  do i=lo(1),hi(1)
     res(i) = -res(i)
  enddo
  
  end subroutine scale_residual_1d

!
! ******************************************************************************************
!

! Multiply residual by -1 in 2d
  subroutine scale_residual_2d(lo,hi,ng,res)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(inout) :: res(lo(1)-ng:,lo(2)-ng:)

! Local
  integer :: i,j

  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        res(i,j) = -res(i,j)
     enddo
  enddo
  
  end subroutine scale_residual_2d

!
! ******************************************************************************************
!

! Multiply residual by -1 in 3d
  subroutine scale_residual_3d(lo,hi,ng,res)

  integer        , intent(in   ) :: lo(:),hi(:),ng
  real(kind=dp_t), intent(inout) :: res(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)

! Local
  integer :: i,j,k

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           res(i,j,k) = -res(i,j,k)
        enddo
     enddo
  enddo

  end subroutine scale_residual_3d

end module ml_cc_module
