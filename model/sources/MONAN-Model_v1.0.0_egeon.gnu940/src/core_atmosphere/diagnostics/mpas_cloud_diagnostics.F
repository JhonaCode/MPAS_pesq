! Copyright (c) 2022, University Corporation for Atmospheric Research (UCAR).
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html
!
module mpas_cloud_diagnostics

    use mpas_dmpar
    use mpas_kind_types, only : RKIND
    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
    use mpas_pool_routines
    use mpas_constants
    use mpas_log, only : mpas_log_write
    use VerticalInterpolation

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_pool_type), pointer :: diag_physics
    type (MPAS_pool_type), pointer :: tend_physics
    type (MPAS_clock_type), pointer :: clock

    public :: cloud_diagnostics_setup, &
              cloud_diagnostics_update, &
              cloud_diagnostics_compute,&
              cloud_diagnostics_reset, & 
              cloud_diagnostics_cleanup

    private

    logical :: need_cldfrac_isobaric
    logical :: need_cldfrac_low_UPP
    logical :: need_cldfrac_mid_UPP
    logical :: need_cldfrac_high_UPP
    logical :: need_cldfrac_tot_UPP

    contains


    !-----------------------------------------------------------------------
    !  routine cloud_diagnostics_setup
    !
    !> \brief Initialize the cloud diagnostic module
    !> \author G. Dylan Dickerson
    !> \date   23 August 2022
    !> \details
    !>  Initialize the diagnostic and save pointers to subpools for
    !>  reuse in this module
    !
    !-----------------------------------------------------------------------
    subroutine cloud_diagnostics_setup(all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock
        !local pointers:
        real (kind=RKIND), dimension(:  ), pointer :: t_iso_levels
        real (kind=RKIND), dimension(:,:), pointer :: cldfrac_isobaric
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_low_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_mid_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_high_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_tot_UPP
        real (kind=RKIND), dimension(:)  , pointer :: nsteps_accum_cloud

        clock => simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
        call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)
        call mpas_pool_get_subpool(all_pools, 'tend_physics', tend_physics)
        call mpas_pool_get_array(diag, 't_iso_levels', t_iso_levels)

        call mpas_pool_get_array(diag, 'cldfrac_isobaric', cldfrac_isobaric)
        call mpas_pool_get_array(diag, 'cldfrac_low_UPP',  cldfrac_low_UPP)
        call mpas_pool_get_array(diag, 'cldfrac_mid_UPP',  cldfrac_mid_UPP)
        call mpas_pool_get_array(diag, 'cldfrac_high_UPP', cldfrac_high_UPP)
        call mpas_pool_get_array(diag, 'cldfrac_tot_UPP',  cldfrac_tot_UPP)
        call mpas_pool_get_array(diag, 'nsteps_accum_cloud',  nsteps_accum_cloud)

              t_iso_levels(:) = (/100000.0, &
                                  97500.0, &
                                  95000.0, &
                                  92500.0, &
                                  90000.0, &
                                  87500.0, &
                                  85000.0, &
                                  82500.0, &
                                  80000.0, &
                                  77500.0, &
                                  75000.0, &
                                  70000.0, &
                                  65000.0, &
                                  60000.0, &
                                  55000.0, &
                                  50000.0, &
                                  45000.0, &
                                  40000.0, &
                                  35000.0, &
                                  30000.0, &
                                  25000.0, &
                                  22500.0, &
                                  20000.0, &
                                  17500.0, &
                                  15000.0, &
                                  12500.0, &
                                  10000.0/)        
      cldfrac_isobaric=0.0
      cldfrac_low_UPP=0.0
      cldfrac_mid_UPP=0.0
      cldfrac_high_UPP=0.0
      cldfrac_tot_UPP=0.0
      nsteps_accum_cloud=0.0
    end subroutine cloud_diagnostics_setup

    !-----------------------------------------------------------------------
    !  routine radiation_diagnostic_update
    !
    !> \brief Handle diagnostic calculation needed after each timestep
    !> \author 
    !> \date   
    !> \details
    !>  Handle diagnostic calculation needed after each timestep
    !
    !-----------------------------------------------------------------------
    subroutine cloud_diagnostics_update()
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none

        ! Called at the end of every timestep
        ! Update extrema, accumulations, etc.
       
        !inout arguments:
        integer                          , pointer :: nCells
        integer                          , pointer :: nCellsSolve
        integer                          , pointer :: nVertLevels
        integer                          , pointer :: nVertices
        integer                          , pointer :: vertexDegree
        integer                          , pointer :: nIsoLevelsT
       
        integer, dimension(:,:)          , pointer :: cellsOnVertex
        real (kind=RKIND), dimension(:)  , pointer :: areaTriangle
        real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
        real (kind=RKIND), dimension(:)  , pointer :: latCell
        real (kind=RKIND), dimension(:)  , pointer :: lonCell
       
        real (kind=RKIND), dimension(:,:), pointer :: height
        real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p 
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_low_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_mid_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_high_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_tot_UPP
        real (kind=RKIND), dimension(:)  , pointer :: nsteps_accum_cloud
	
        real (kind=RKIND), dimension(:  ), pointer :: t_iso_levels
        real (kind=RKIND), dimension(:,:), pointer :: cldfrac
        real (kind=RKIND), dimension(:,:), pointer :: cldfrac_isobaric 

        real (kind=RKIND), dimension(:)  , allocatable :: cld_low_UPP
        real (kind=RKIND), dimension(:)  , allocatable :: cld_mid_UPP
        real (kind=RKIND), dimension(:)  , allocatable :: cld_high_UPP
        real (kind=RKIND), dimension(:)  , allocatable :: cld_tot_UPP

        real (kind=RKIND), dimension(:,:), pointer :: pressure_base        
        !-----------------------------------------------------------! 
        real (kind=RKIND), dimension(:,:), allocatable :: pressure, pressure2
       
        !local interpolated fields:
        real (kind=RKIND) :: w1,w2,z0,z1,z2
        real (kind=RKIND), dimension(:,:), allocatable :: field_interp
        integer :: nVertLevelsP1,iCell,k

       !--------------------------------------------------------------------------------------------------

        ! levels for low/mid/high cloud fraction - UPP method
        real (kind=RKIND), parameter :: ptop_low = 64200.0, ptop_mid = 35000.0, ptop_high = 15000.0

        logical :: need_cldfrac_UPP
        logical :: need_any_diags

        need_any_diags = .false.
        need_cldfrac_low_UPP= MPAS_field_will_be_written('cldfrac_low_UPP')
        need_any_diags = need_any_diags .or. need_cldfrac_low_UPP
        need_cldfrac_mid_UPP  = MPAS_field_will_be_written('cldfrac_mid_UPP') 
        need_any_diags = need_any_diags .or. need_cldfrac_mid_UPP
        need_cldfrac_high_UPP = MPAS_field_will_be_written('cldfrac_high_UPP') 
        need_any_diags = need_any_diags .or. need_cldfrac_high_UPP
        need_cldfrac_tot_UPP  = MPAS_field_will_be_written('cldfrac_tot_UPP')  
        need_any_diags = need_any_diags .or. need_cldfrac_tot_UPP

      !  if (need_any_diags) then
           call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
           call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

           call mpas_pool_get_array(diag, 'cldfrac_low_UPP',  cldfrac_low_UPP)
           call mpas_pool_get_array(diag, 'cldfrac_mid_UPP',  cldfrac_mid_UPP)
           call mpas_pool_get_array(diag, 'cldfrac_high_UPP', cldfrac_high_UPP)
           call mpas_pool_get_array(diag, 'cldfrac_tot_UPP',  cldfrac_tot_UPP)
           call mpas_pool_get_array(diag, 'nsteps_accum_cloud',  nsteps_accum_cloud)

           call mpas_pool_get_array(diag_physics, 'cldfrac', cldfrac)

           call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
           call mpas_pool_get_array(diag, 'pressure_p'   , pressure_p)

           if(.not.allocated(pressure)    ) allocate(pressure(nVertLevels,nCellsSolve))
           if(.not.allocated(cld_low_UPP) ) allocate(cld_low_UPP(nCellsSolve)) 
           if(.not.allocated(cld_mid_UPP) ) allocate(cld_mid_UPP(nCellsSolve)) 
           if(.not.allocated(cld_high_UPP)) allocate(cld_high_UPP(nCellsSolve))
           if(.not.allocated(cld_tot_UPP) ) allocate(cld_tot_UPP(nCellsSolve)) 
           do iCell = 1, nCellsSolve
               cld_low_UPP (iCell) = 0.0
               cld_mid_UPP (iCell) = 0.0
               cld_high_UPP(iCell) = 0.0
               cld_tot_UPP (iCell) = 0.0
               do k = 1, nVertLevels
                  pressure(k,iCell) = pressure_p(k,iCell) + pressure_base(k,iCell)

                  if ( pressure(k,iCell) >= ptop_low ) then
                      cld_low_UPP(iCell)  = max(cld_low_UPP(iCell), cldfrac(k,iCell))
                  else if ( pressure(k,iCell) < ptop_low .and. pressure(k,iCell) >= ptop_mid  ) then
                      cld_mid_UPP(iCell)  = max(cld_mid_UPP(iCell), cldfrac(k,iCell))
                  else if ( pressure(k,iCell) < ptop_mid .and. pressure(k,iCell) >= ptop_high ) then
                      cld_high_UPP(iCell) = max(cld_high_UPP(iCell), cldfrac(k,iCell))
                  end if
                  cld_tot_UPP(iCell) = max(cld_tot_UPP(iCell), cldfrac(k,iCell))
               end do
           end do           
           do iCell = 1, nCellsSolve
              cldfrac_low_UPP (iCell) = cldfrac_low_UPP (iCell) + cld_low_UPP (iCell)
              cldfrac_mid_UPP (iCell) = cldfrac_mid_UPP (iCell) + cld_mid_UPP (iCell)
              cldfrac_high_UPP(iCell) = cldfrac_high_UPP(iCell) + cld_high_UPP(iCell)
              cldfrac_tot_UPP (iCell) = cldfrac_tot_UPP (iCell) + cld_tot_UPP (iCell)
              nsteps_accum_cloud(iCell) = nsteps_accum_cloud (iCell) + 1.0
           end do   
           if(allocated(pressure    )) deallocate(pressure)
           if(allocated(cld_low_UPP )) deallocate(cld_low_UPP) 
           if(allocated(cld_mid_UPP )) deallocate(cld_mid_UPP) 
           if(allocated(cld_high_UPP)) deallocate(cld_high_UPP)
           if(allocated(cld_tot_UPP )) deallocate(cld_tot_UPP) 

      !  end if

        need_any_diags = .false.
        need_cldfrac_isobaric = MPAS_field_will_be_written('cldfrac_isobaric')
        need_any_diags = need_any_diags .or. need_cldfrac_isobaric

       ! if (need_any_diags) then
           ! call interp_diagnostics(mesh, state, 1, diag,diag_physics, tend_physics)

       
           call mpas_pool_get_dimension(mesh, 'nCells', nCells)
           call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
           call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
           call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
           call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
           call mpas_pool_get_dimension(mesh, 'nIsoLevelsT', nIsoLevelsT)
       
           nVertLevelsP1 = nVertLevels + 1
           call mpas_pool_get_array(mesh,'latCell',latCell)
           call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
           call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
           call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
       
           call mpas_pool_get_array(mesh, 'zgrid', height)
              
           call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
           call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
       
           call mpas_pool_get_array(diag_physics, 'cldfrac', cldfrac)

           call mpas_pool_get_array(diag, 't_iso_levels', t_iso_levels)
           call mpas_pool_get_array(diag, 'cldfrac_isobaric', cldfrac_isobaric)

           if(.not.allocated(pressure)    ) allocate(pressure(nVertLevels,nCells)      )
           if(.not.allocated(pressure2)   ) allocate(pressure2(nVertLevelsP1,nCells)   )
           !calculation of total pressure at cell centers (at mass points):
           do iCell = 1, nCells
              do k = 1, nVertLevels
                 pressure(k,iCell)    = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
              enddo
           enddo
 
           !calculation of total pressure at cell centers (at vertical velocity points):
           k = nVertLevelsP1
           do iCell = 1, nCells
              z0 = height(k,iCell)
              z1 = 0.5*(height(k,iCell)+height(k-1,iCell)) 
              z2 = 0.5*(height(k-1,iCell)+height(k-2,iCell))
              w1 = (z0-z2)/(z1-z2)
              w2 = 1.-w1
              !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
              pressure2(k,iCell) = exp(w1*log(pressure(k-1,iCell))+w2*log(pressure(k-2,iCell)))
           enddo
           do k = 2, nVertLevels
              do iCell = 1, nCells
                 w1 = (height(k,iCell)-height(k-1,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
                 w2 = (height(k+1,iCell)-height(k,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
                 ! pressure2(k,iCell) = w1*pressure(k,iCell) + w2*pressure(k-1,iCell)
                 !
                 ! switch to use ln(pressure) for more accurate vertical interpolation, WCS 20230407
                 pressure2(k,iCell) = exp(w1*log(pressure(k,iCell))+w2*log(pressure(k-1,iCell)))
              enddo
           enddo
           k = 1
           do iCell = 1, nCells
              z0 = height(k,iCell)
              z1 = 0.5*(height(k,iCell)+height(k+1,iCell)) 
              z2 = 0.5*(height(k+1,iCell)+height(k+2,iCell))
              w1 = (z0-z2)/(z1-z2)
              w2 = 1.-w1
              ! pressure2(k,iCell) = w1*pressure(k,iCell)+w2*pressure(k+1,iCell)
              !
              ! switch to use ln(pressure) for more accurate vertical interpolation, WCS 20230407
              pressure2(k,iCell) = exp(w1*log(pressure(k,iCell))+w2*log(pressure(k+1,iCell)))
           enddo
          !!!!!!!!!!! Additional temperature levels for vortex tracking !!!!!!!!!!!
                     
                 
        if(allocated(field_interp)) deallocate(field_interp)
        if(.not.allocated(field_interp))allocate(field_interp(nIsoLevelsT,nCells))
                    
        CALL VerticalInterpolation_Single(nCells                                                ,&
                                          nVertLevels                                           ,&
                                          nVertLevels                                           ,&
                                          nIsoLevelsT                                           ,&
                                          nIsoLevelsT                                           ,&
                                          cldfrac            (1:nVertLevels,1:nCells)           ,&
                                          field_interp       (1:nIsoLevelsT,1:nCells)           ,&
                                          pressure           (1:nVertLevels,1:nCells) * 100.0   ,&
                                          t_iso_levels       (1:nIsoLevelsT)                     )

           do iCell = 1, nCells
              do k=1, nIsoLevelsT
                 cldfrac_isobaric (k,iCell) = cldfrac_isobaric (k,iCell) +  field_interp(k,iCell)
              end do
           end do   

         if(allocated(field_interp ))deallocate(field_interp)

    end subroutine cloud_diagnostics_update
    !-----------------------------------------------------------------------
    !  routine cloud_diagnostics_compute
    !
    !> \brief Compute diagnostic before model output is written
    !> \author G. Dylan Dickerson
    !> \date   23 August 2022
    !> \details
    !>  Compute diagnostic before model output is written
    !>  The following fields are computed by this routine:
    !>     cldfrac_low_UPP
    !>     cldfrac_mid_UPP
    !>     cldfrac_high_UPP
    !>     cldfrac_tot_UPP
    !
    !-----------------------------------------------------------------------
    subroutine cloud_diagnostics_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

        implicit none

        logical :: need_any_diags

        need_any_diags = .false.
        need_cldfrac_low_UPP= MPAS_field_will_be_written('cldfrac_low_UPP')
        need_any_diags = need_any_diags .or. need_cldfrac_low_UPP
        need_cldfrac_mid_UPP  = MPAS_field_will_be_written('cldfrac_mid_UPP') 
        need_any_diags = need_any_diags .or. need_cldfrac_mid_UPP
        need_cldfrac_high_UPP = MPAS_field_will_be_written('cldfrac_high_UPP') 
        need_any_diags = need_any_diags .or. need_cldfrac_high_UPP
        need_cldfrac_tot_UPP  = MPAS_field_will_be_written('cldfrac_tot_UPP')  
        need_any_diags = need_any_diags .or. need_cldfrac_tot_UPP
        need_cldfrac_isobaric = MPAS_field_will_be_written('cldfrac_isobaric')
        need_any_diags = need_any_diags .or. need_cldfrac_isobaric

        if (need_any_diags) then
            call interp_diagnostics(mesh, state, diag,diag_physics, tend_physics)
        end if

    end subroutine cloud_diagnostics_compute





   !==================================================================================================
    subroutine interp_diagnostics(mesh, state, diag, diag_physics, tend_physics)
   !==================================================================================================

       !input arguments:
        type (mpas_pool_type), intent(in)  :: mesh
        type (mpas_pool_type), intent(in) :: state
        type (mpas_pool_type), intent(in) :: diag_physics
        type (mpas_pool_type), intent(in) :: tend_physics
       
       !inout arguments:
        type (mpas_pool_type), intent(inout) :: diag
       
       !local variables:
        integer :: iCell,iVert,iVertD,k,kk
        integer, pointer :: nCells, nCellsSolve, nVertLevels, nVertices, vertexDegree, nIsoLevelsT
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_low_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_mid_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_high_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_tot_UPP
        real (kind=RKIND), dimension(:)  , pointer :: nsteps_accum_cloud
        real (kind=RKIND), dimension(:,:), pointer :: cldfrac_isobaric 
                                                      
                      
       ! call mpas_log_write('')
       ! call mpas_log_write('--- enter subroutine interp_diagnostics:')
       
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
        call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
        call mpas_pool_get_dimension(mesh, 'nIsoLevelsT', nIsoLevelsT)

        call mpas_pool_get_array(diag, 'cldfrac_isobaric', cldfrac_isobaric)
        call mpas_pool_get_array(diag, 'cldfrac_low_UPP' , cldfrac_low_UPP )
        call mpas_pool_get_array(diag, 'cldfrac_mid_UPP' , cldfrac_mid_UPP )
        call mpas_pool_get_array(diag, 'cldfrac_high_UPP', cldfrac_high_UPP)
        call mpas_pool_get_array(diag, 'cldfrac_tot_UPP' , cldfrac_tot_UPP )
        call mpas_pool_get_array(diag, 'nsteps_accum_cloud',  nsteps_accum_cloud)
       
       
        do iCell=1,nCells          
           if(nsteps_accum_cloud(iCell) /= 0._RKIND)then
              do iVert=1,nIsoLevelsT
                  cldfrac_isobaric(iVert,iCell) = cldfrac_isobaric(iVert,iCell)/nsteps_accum_cloud(iCell)
              end do
           endif
        end do   
        do iCell=1,nCells
          if(nsteps_accum_cloud(iCell) /= 0._RKIND)then
            !short-wave radiation:
            cldfrac_low_UPP (iCell) = cldfrac_low_UPP (iCell)/nsteps_accum_cloud(iCell) 
            cldfrac_mid_UPP (iCell) = cldfrac_mid_UPP (iCell)/nsteps_accum_cloud(iCell) 
            cldfrac_high_UPP(iCell) = cldfrac_high_UPP(iCell)/nsteps_accum_cloud(iCell) 
            cldfrac_tot_UPP (iCell) = cldfrac_tot_UPP (iCell)/nsteps_accum_cloud(iCell) 
           endif
        end do

    end subroutine interp_diagnostics


    !-----------------------------------------------------------------------
    !  routine cloud_diagnostic_reset
    !
    !> \brief Reset diagnostic after it has been written
    !> \author 
    !> \date   
    !> \details
    !>  Reset diagnostic after it has been written
    !
    !-----------------------------------------------------------------------
    subroutine cloud_diagnostics_reset()
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        ! Called immediately after diagnostics have been written
        ! Reset counters, accumulations, etc.

        implicit none
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_low_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_mid_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_high_UPP
        real (kind=RKIND), dimension(:)  , pointer :: cldfrac_tot_UPP
        real (kind=RKIND), dimension(:)  , pointer :: nsteps_accum_cloud
        real (kind=RKIND), dimension(:,:), pointer :: cldfrac_isobaric 
        logical :: need_any_diags

        need_any_diags = .false.
        need_cldfrac_low_UPP= MPAS_field_will_be_written('cldfrac_low_UPP')
        need_any_diags = need_any_diags .or. need_cldfrac_low_UPP
        need_cldfrac_mid_UPP  = MPAS_field_will_be_written('cldfrac_mid_UPP') 
        need_any_diags = need_any_diags .or. need_cldfrac_mid_UPP
        need_cldfrac_high_UPP = MPAS_field_will_be_written('cldfrac_high_UPP') 
        need_any_diags = need_any_diags .or. need_cldfrac_high_UPP
        need_cldfrac_tot_UPP  = MPAS_field_will_be_written('cldfrac_tot_UPP')  
        need_any_diags = need_any_diags .or. need_cldfrac_tot_UPP
        need_cldfrac_isobaric = MPAS_field_will_be_written('cldfrac_isobaric')
        need_any_diags = need_any_diags .or. need_cldfrac_isobaric
        if(need_any_diags)then
           call mpas_pool_get_array(diag, 'cldfrac_isobaric', cldfrac_isobaric)
           call mpas_pool_get_array(diag, 'cldfrac_low_UPP' , cldfrac_low_UPP )
           call mpas_pool_get_array(diag, 'cldfrac_mid_UPP' , cldfrac_mid_UPP )
           call mpas_pool_get_array(diag, 'cldfrac_high_UPP', cldfrac_high_UPP)
           call mpas_pool_get_array(diag, 'cldfrac_tot_UPP' , cldfrac_tot_UPP )
           call mpas_pool_get_array(diag, 'nsteps_accum_cloud',  nsteps_accum_cloud)
           cldfrac_isobaric=0.0
           cldfrac_low_UPP =0.0
           cldfrac_mid_UPP =0.0
           cldfrac_high_UPP=0.0
           cldfrac_tot_UPP =0.0
           nsteps_accum_cloud=0.0
        end if 

       call reset_diagnostics(diag_physics )

    end subroutine cloud_diagnostics_reset

   !==================================================================================================
    subroutine reset_diagnostics( diag_physics)
   !==================================================================================================

       !input arguments:
       type (mpas_pool_type), intent(in) :: diag_physics
       ! Called immediately after diagnostics have been written
       ! Reset counters, accumulations, etc.


    end subroutine reset_diagnostics

    !-----------------------------------------------------------------------
    !  routine radiation_diagnostic_cleanup
    !
    !> \brief Finalizes diagnostic
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Finalizes diagnostic
    !
    !-----------------------------------------------------------------------
    subroutine cloud_diagnostics_cleanup()

        implicit none

        ! Deallocate scratch arrays, etc.
   
    end subroutine cloud_diagnostics_cleanup




   !==================================================================================================
    subroutine interp_tofixed_pressure(ncol,nlev_in,nlev_out,pres_in,field_in,pres_out,field_out)
   !==================================================================================================
   
   !input arguments:
    integer,intent(in):: ncol,nlev_in,nlev_out
   
    real(kind=RKIND),intent(in),dimension(ncol,nlev_in) :: pres_in,field_in
    real(kind=RKIND),intent(in),dimension(ncol,nlev_out):: pres_out
   
   !output arguments:
    real(kind=RKIND),intent(out),dimension(ncol,nlev_out):: field_out
   
   !local variables:
   ! integer:: i1,i2,icol,k,kk
    integer:: icol,k,kk
    integer:: kkstart,kount
    integer,dimension(ncol):: kupper
   
    real(kind=RKIND):: dpl,dpu
   
   !--------------------------------------------------------------------------------------------------
   
   !call mpas_log_write('')
   !call mpas_log_write('--- enter subroutine interp_tofixed_pressure:')
   !call mpas_log_write('... ncol     = $i',intArgs=(/ncol/))
   !call mpas_log_write('... nlev_in  = $i',intArgs=(/nlev_in/))
   !call mpas_log_write('... nlev_out = $i',intArgs=(/nlev_out/))
   !i1=1 ; i2=ncol
   !do k = 1, nlev_in
   !   call mpas_log_write('$i $r $r $r $r', intArgs=(/k/), realArgs=(/pres_in(i1,k),field_in(i1,k),pres_in(i2,k),field_in(i2,k)/))
   !enddo
   !call mpas_log_write('')
   
    do icol = 1, ncol
       kupper(icol) = 1
    enddo
   
    do k = 1, nlev_out 
   
       kkstart = nlev_in
       do icol = 1, ncol
          kkstart = min0(kkstart,kupper(icol))
       enddo
       kount = 0
   
       do kk = kkstart, nlev_in-1
          do icol = 1, ncol
             if(pres_out(icol,k).gt.pres_in(icol,kk).and.pres_out(icol,k).le.pres_in(icol,kk+1)) then
                kupper(icol) = kk
                kount = kount + 1
   !            call mpas_log_write('$i $i $r $r $r', intArgs=(/k,kupper(icol)/), realArgs=(/pres_out(icol,k),pres_in(icol,kk),pres_in(icol,kk+1)/))
             endif
          enddo
   
          if(kount.eq.ncol) then
             do icol = 1, ncol
                dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
                dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
                field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                                  + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
             end do
             goto 35
           end if
       enddo
   
       do icol = 1, ncol
          if(pres_out(icol,k) .lt. pres_in(icol,1)) then
             field_out(icol,k) = field_in(icol,1)*pres_out(icol,k)/pres_in(icol,1)
          elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in)) then
             field_out(icol,k) = field_in(icol,nlev_in)
          else
             dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
             dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
             field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                               + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
          endif
       enddo
   
    35 continue
   !   call mpas_log_write('$i $r $r $r $r $r $r', intArgs=(/kupper(i1)/), &
   !                       realArgs=(/pres_out(i1,k),pres_in(i1,kupper(i1)),pres_in(i1,kupper(i1)+1),field_out(i1,k),field_in(i1,kupper(i1)),field_in(i1,kupper(i1)+1)/))
   !   call mpas_log_write('$i $r $r $r $r $r $r', intArgs=(/kupper(i2)/), &
   !                       realArgs=(/pres_out(i2,k),pres_in(i2,kupper(i2)),pres_in(i2,kupper(i2)+1),field_out(i2,k),field_in(i2,kupper(i2)),field_in(i2,kupper(i2)+1)/))
   
    enddo
   
    end subroutine interp_tofixed_pressure
   

end module mpas_cloud_diagnostics
