! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_radiation_diagnostics

    use mpas_dmpar
    use mpas_kind_types
    use mpas_derived_types
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

    public :: radiation_diagnostics_setup, &
              radiation_diagnostics_update, &
              radiation_diagnostics_compute,&
              radiation_diagnostics_reset,&
              radiation_diagnostics_cleanup

    private

    logical :: need_rthratenlw_isobaric   , &
               need_rthratensw_isobaric   , &
               need_acswdnb               , &  
               need_acswdnbc              , &  
               need_acswdnt               , &  
               need_acswdntc              , &  
               need_acswupb               , &  
               need_acswupbc              , &  
               need_acswupt               , &  
               need_acswuptc              , &  
               need_aclwdnb               , &  
               need_aclwdnbc              , &  
               need_aclwdnt               , &  
               need_aclwdntc              , &  
               need_aclwupb               , &  
               need_aclwupbc              , &  
               need_aclwupt               , &  
               need_aclwuptc
    real (kind=RKIND),pointer:: config_bucket_radt
    character(len=StrKIND),pointer:: config_bucket_update

    contains


    !-----------------------------------------------------------------------
    !  routine radiation_diagnostics_setup
    !
    !> \brief Set up the isobaric diagnostics module
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  This routine sets up the isobaric diagnostics module, principally by
    !>  saving pointers to pools that are used in the computation of diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine radiation_diagnostics_setup(all_pools, simulation_clock,domain)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock
        type(domain_type),intent(in):: domain
        !local pointers:
        real (kind=RKIND), dimension(:  ), pointer :: t_iso_levels
        real (kind=RKIND), dimension(:,:), pointer :: rthratenlw_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rthratensw_isobaric
        real (kind=RKIND), dimension(:)  , pointer :: acswdnb,acswdnbc,acswdnt,acswdntc,  &
                                                      acswupb,acswupbc,acswupt,acswuptc,  &
                                                      aclwdnb,aclwdnbc,aclwdnt,aclwdntc,  &
                                                      aclwupb,aclwupbc,aclwupt,aclwuptc
        real (kind=RKIND), dimension(:)  , pointer :: nsteps_accum_rad


        clock => simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
	call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)
        call mpas_pool_get_subpool(all_pools, 'tend_physics', tend_physics)
        call mpas_pool_get_config(domain%configs,'config_bucket_radt'      ,config_bucket_radt      )
        call mpas_pool_get_config(domain%configs,'config_bucket_update'    ,config_bucket_update    )

        call mpas_pool_get_array(diag, 'rthratenlw_isobaric', rthratenlw_isobaric)
        call mpas_pool_get_array(diag, 'rthratensw_isobaric', rthratensw_isobaric)
        !short-wave radiation:
        call mpas_pool_get_array(diag_physics,'acswdnb'   , acswdnb   )
        call mpas_pool_get_array(diag_physics,'acswdnbc'  , acswdnbc  )
        call mpas_pool_get_array(diag_physics,'acswdnt'   , acswdnt   )
        call mpas_pool_get_array(diag_physics,'acswdntc'  , acswdntc  )
        call mpas_pool_get_array(diag_physics,'acswupb'   , acswupb   )
        call mpas_pool_get_array(diag_physics,'acswupbc'  , acswupbc  )
        call mpas_pool_get_array(diag_physics,'acswupt'   , acswupt   )
        call mpas_pool_get_array(diag_physics,'acswuptc'  , acswuptc  )

        !long-wave radiation:
        call mpas_pool_get_array(diag_physics,'aclwdnb'   , aclwdnb   )
        call mpas_pool_get_array(diag_physics,'aclwdnbc'  , aclwdnbc  )
        call mpas_pool_get_array(diag_physics,'aclwdnt'   , aclwdnt   )
        call mpas_pool_get_array(diag_physics,'aclwdntc'  , aclwdntc  )
        call mpas_pool_get_array(diag_physics,'aclwupb'   , aclwupb   )
        call mpas_pool_get_array(diag_physics,'aclwupbc'  , aclwupbc  )
        call mpas_pool_get_array(diag_physics,'aclwupt'   , aclwupt   )
        call mpas_pool_get_array(diag_physics,'aclwuptc'  , aclwuptc  )

        call mpas_pool_get_array(diag_physics,'nsteps_accum_rad',nsteps_accum_rad)

        call mpas_pool_get_array(diag, 't_iso_levels', t_iso_levels)

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
      rthratenlw_isobaric=0.0
      rthratensw_isobaric=0.0
      acswdnb  =0.0 
      acswdnbc =0.0 
      acswdnt  =0.0 
      acswdntc =0.0 
      acswupb  =0.0 
      acswupbc =0.0 
      acswupt  =0.0 
      acswuptc =0.0 

      aclwdnb  =0.0 
      aclwdnbc =0.0 
      aclwdnt  =0.0 
      aclwdntc =0.0 
      aclwupb  =0.0 
      aclwupbc =0.0 
      aclwupt  =0.0 
      aclwuptc =0.0 
      nsteps_accum_rad=0.0
    end subroutine radiation_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine radiation_diagnostics_update
    !
    !> \brief Handle diagnostic calculation needed after each timestep
    !> \author 
    !> \date   
    !> \details
    !>  Handle diagnostic calculation needed after each timestep
    !
    !-----------------------------------------------------------------------
    subroutine radiation_diagnostics_update()
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
        real (kind=RKIND), dimension(:,:), pointer :: rthratenlw, rthratensw, rthblten,rqvblten
 
       !local variables:
     
        real (kind=RKIND), dimension(:  ), pointer :: t_iso_levels
        real (kind=RKIND), dimension(:,:), pointer :: rthratenlw_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rthratensw_isobaric
	
        real (kind=RKIND), dimension(:)  , pointer :: acswdnb,acswdnbc,acswdnt,acswdntc,  &
                                                      acswupb,acswupbc,acswupt,acswuptc,  &
                                                      aclwdnb,aclwdnbc,aclwdnt,aclwdntc,  &
                                                      aclwupb,aclwupbc,aclwupt,aclwuptc
        real(kind=RKIND),dimension(:),pointer:: swdnb,swdnbc,swdnt,swdntc,  &
                                                swupb,swupbc,swupt,swuptc,  &
                                                lwdnb,lwdnbc,lwdnt,lwdntc,  &
                                                lwupb,lwupbc,lwupt,lwuptc

	real (kind=RKIND), dimension(:)  , pointer :: nsteps_accum_rad
                                                      
               
        !-----------------------------------------------------------! 
        
        !local interpolated fields:
        real (kind=RKIND), dimension(:,:), allocatable :: pressure, pressure2

        real (kind=RKIND)                              :: w1,w2,z0,z1,z2
        real (kind=RKIND), dimension(:,:), allocatable :: field_in,press_in
        real (kind=RKIND), dimension(:,:), allocatable :: field_interp,press_interp

        integer               :: iCell,iVert,iVertD,k,kk
        integer               :: nVertLevelsP1
              
         
       !--------------------------------------------------------------------------------------------------
       
       ! call mpas_log_write('')
       ! call mpas_log_write('--- enter subroutine interp_diagnostics:')
       
        call mpas_pool_get_dimension(mesh, 'nCells'      , nCells)
        call mpas_pool_get_dimension(mesh, 'nCellsSolve' , nCellsSolve)
        call mpas_pool_get_dimension(mesh, 'nVertLevels' , nVertLevels)
        call mpas_pool_get_dimension(mesh, 'nVertices'   , nVertices)
        call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
        call mpas_pool_get_dimension(mesh, 'nIsoLevelsT' , nIsoLevelsT)
       
        nVertLevelsP1 = nVertLevels + 1

        call mpas_pool_get_array(mesh, 'cellsOnVertex'    , cellsOnVertex)
        call mpas_pool_get_array(mesh, 'areaTriangle'     , areaTriangle)
        call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
        call mpas_pool_get_array(mesh, 'latCell'          , latCell)
	call mpas_pool_get_array(mesh, 'lonCell'          , lonCell)  
 
        call mpas_pool_get_array(mesh, 'zgrid'        , height)             
        call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
        call mpas_pool_get_array(diag, 'pressure_p'   , pressure_p)

       
        call mpas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
        call mpas_pool_get_array(tend_physics, 'rthratensw', rthratensw)
        call mpas_pool_get_array(tend_physics, 'rthblten'  , rthblten)
        call mpas_pool_get_array(tend_physics, 'rqvblten'  , rqvblten)

        call mpas_pool_get_array(diag, 't_iso_levels', t_iso_levels)
        call mpas_pool_get_array(diag, 'rthratenlw_isobaric', rthratenlw_isobaric)
        call mpas_pool_get_array(diag, 'rthratensw_isobaric', rthratensw_isobaric)

        !short-wave radiation:
        call mpas_pool_get_array(diag_physics,'acswdnb'   , acswdnb   )
        call mpas_pool_get_array(diag_physics,'acswdnbc'  , acswdnbc  )
        call mpas_pool_get_array(diag_physics,'acswdnt'   , acswdnt   )
        call mpas_pool_get_array(diag_physics,'acswdntc'  , acswdntc  )
        call mpas_pool_get_array(diag_physics,'acswupb'   , acswupb   )
        call mpas_pool_get_array(diag_physics,'acswupbc'  , acswupbc  )
        call mpas_pool_get_array(diag_physics,'acswupt'   , acswupt   )
        call mpas_pool_get_array(diag_physics,'acswuptc'  , acswuptc  )

        !long-wave radiation:
        call mpas_pool_get_array(diag_physics,'aclwdnb'   , aclwdnb   )
        call mpas_pool_get_array(diag_physics,'aclwdnbc'  , aclwdnbc  )
        call mpas_pool_get_array(diag_physics,'aclwdnt'   , aclwdnt   )
        call mpas_pool_get_array(diag_physics,'aclwdntc'  , aclwdntc  )
        call mpas_pool_get_array(diag_physics,'aclwupb'   , aclwupb   )
        call mpas_pool_get_array(diag_physics,'aclwupbc'  , aclwupbc  )
        call mpas_pool_get_array(diag_physics,'aclwupt'   , aclwupt   )
        call mpas_pool_get_array(diag_physics,'aclwuptc'  , aclwuptc  )

        !short-wave radiation:

        call mpas_pool_get_array(diag_physics,'swdnb'     , swdnb     )
        call mpas_pool_get_array(diag_physics,'swdnbc'    , swdnbc    )
        call mpas_pool_get_array(diag_physics,'swdnt'     , swdnt     )
        call mpas_pool_get_array(diag_physics,'swdntc'    , swdntc    )
        call mpas_pool_get_array(diag_physics,'swupb'     , swupb     )
        call mpas_pool_get_array(diag_physics,'swupbc'    , swupbc    )
        call mpas_pool_get_array(diag_physics,'swupt'     , swupt     )
        call mpas_pool_get_array(diag_physics,'swuptc'    , swuptc    )
        !long-wave radiation:

        call mpas_pool_get_array(diag_physics,'lwdnb'     , lwdnb     )
        call mpas_pool_get_array(diag_physics,'lwdnbc'    , lwdnbc    )
        call mpas_pool_get_array(diag_physics,'lwdnt'     , lwdnt     )
        call mpas_pool_get_array(diag_physics,'lwdntc'    , lwdntc    )
        call mpas_pool_get_array(diag_physics,'lwupb'     , lwupb     )
        call mpas_pool_get_array(diag_physics,'lwupbc'    , lwupbc    )
        call mpas_pool_get_array(diag_physics,'lwupt'     , lwupt     )
        call mpas_pool_get_array(diag_physics,'lwuptc'    , lwuptc    )

        call mpas_pool_get_array(diag_physics,'nsteps_accum_rad',nsteps_accum_rad)

        
       
        if(.not.allocated(pressure)    ) allocate(pressure (nVertLevels  ,nCells) )
        if(.not.allocated(pressure2)   ) allocate(pressure2(nVertLevelsP1,nCells) )
       
            !Changed by RRB to add more isobaric levels
             !t_iso_levels(:) = (/10000.0, 12500.0, 15000.0, 17500.0, &
             !                    20000.0, 22500.0, 25000.0, &
             !                    30000.0, 35000.0, &
             !                    40000.0, 45000.0, &
             !                    50000.0, 55000.0, &
             !                    60000.0, 65000.0, &
             !                    70000.0, 75000.0, 77500.0, &
             !                    80000.0, 82500.0, 85000.0, 87500.0, & 
             !                    90000.0, 92500.0, 95000.0, 97500.0, 100000.0/)
             !t_iso_levels(:) = (/100000.0, &
             !                    97500.0, &
             !                    95000.0, &
             !                    92500.0, &
             !                    90000.0, &
             !                    87500.0, &
             !                    85000.0, &
             !                    82500.0, &
             !                    80000.0, &
             !                    77500.0, &
             !                    75000.0, &
             !                    70000.0, &
             !                    65000.0, &
             !                    60000.0, &
             !                    55000.0, &
             !                    50000.0, &
             !                    45000.0, &
             !                    40000.0, &
             !                    35000.0, &
             !                    30000.0, &
             !                    25000.0, &
             !                    22500.0, &
             !                    20000.0, &
             !                    17500.0, &
             !                    15000.0, &
             !                    12500.0, &
             !                    10000.0/)
!
            
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
           z1 = 0.5*(height(k  ,iCell)+height(k-1,iCell)) 
           z2 = 0.5*(height(k-1,iCell)+height(k-2,iCell))
           w1 = (z0-z2)/(z1-z2)
           w2 = 1.0-w1
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
           z1 = 0.5*(height(k  ,iCell)+height(k+1,iCell)) 
           z2 = 0.5*(height(k+1,iCell)+height(k+2,iCell))
           w1 = (z0-z2)/(z1-z2)
           w2 = 1.0-w1
           ! pressure2(k,iCell) = w1*pressure(k,iCell)+w2*pressure(k+1,iCell)
           !
           ! switch to use ln(pressure) for more accurate vertical interpolation, WCS 20230407
           pressure2(k,iCell) = exp(w1*log(pressure(k,iCell))+w2*log(pressure(k+1,iCell)))
        enddo
                      
          
        if(allocated(field_in)) deallocate(field_in)
        if(allocated(press_in)) deallocate(press_in)
       
        if(allocated(field_interp)) deallocate(field_interp)
        if(allocated(press_interp)) deallocate(press_interp)
                      
        !!!!!!!!!!! Additional temperature levels for vortex tracking !!!!!!!!!!!
        if(.not.allocated(field_interp) ) allocate(field_interp(nIsoLevelsT,nCells)   )
          
 
 

        CALL VerticalInterpolation_Single(nCells                                                ,&
                                          nVertLevels                                           ,&
                                          nVertLevels                                           ,&
                                          nIsoLevelsT                                           ,&
                                          nIsoLevelsT                                           ,&
                                          rthratenlw         (1:nVertLevels,1:nCells)           ,&
                                          field_interp       (1:nIsoLevelsT,1:nCells)           ,&
                                          pressure           (1:nVertLevels,1:nCells) * 100.0   ,&
                                          t_iso_levels       (1:nIsoLevelsT)                     )
        do iCell = 1,nCells
           do k = 1,nIsoLevelsT
              rthratenlw_isobaric(k,iCell) = rthratenlw_isobaric(k,iCell) + field_interp(k,iCell) 
           end do
        end do


        CALL VerticalInterpolation_Single(nCells                                                ,&
                                          nVertLevels                                           ,&
                                          nVertLevels                                           ,&
                                          nIsoLevelsT                                           ,&
                                          nIsoLevelsT                                           ,&
                                          rthratensw         (1:nVertLevels,1:nCells)           ,&
                                          field_interp       (1:nIsoLevelsT,1:nCells)           ,&
                                          pressure           (1:nVertLevels,1:nCells) * 100.0   ,&
                                          t_iso_levels       (1:nIsoLevelsT)                     )

        do iCell = 1,nCells
           do k = 1,nIsoLevelsT
              rthratensw_isobaric(k,iCell) = rthratensw_isobaric(k,iCell) + field_interp(k,iCell) 
           end do
        end do
 
 
 
 
        do iCell = 1, nCells
           !short-wave radiation:
           acswdnb (iCell) = acswdnb (iCell) + swdnb (iCell)
           acswdnbc(iCell) = acswdnbc(iCell) + swdnbc(iCell)
           acswdnt (iCell) = acswdnt (iCell) + swdnt (iCell)
           acswdntc(iCell) = acswdntc(iCell) + swdntc(iCell)
           acswupb (iCell) = acswupb (iCell) + swupb (iCell)
           acswupbc(iCell) = acswupbc(iCell) + swupbc(iCell)
           acswupt (iCell) = acswupt (iCell) + swupt (iCell)
           acswuptc(iCell) = acswuptc(iCell) + swuptc(iCell)
           !long-wave radiation:
           aclwdnb (iCell) = aclwdnb (iCell) + lwdnb (iCell)
           aclwdnbc(iCell) = aclwdnbc(iCell) + lwdnbc(iCell)
           aclwdnt (iCell) = aclwdnt (iCell) + lwdnt (iCell)
           aclwdntc(iCell) = aclwdntc(iCell) + lwdntc(iCell)
           aclwupb (iCell) = aclwupb (iCell) + lwupb (iCell)
           aclwupbc(iCell) = aclwupbc(iCell) + lwupbc(iCell)
           aclwupt (iCell) = aclwupt (iCell) + lwupt (iCell)
           aclwuptc(iCell) = aclwuptc(iCell) + lwuptc(iCell)

           nsteps_accum_rad(iCell) = nsteps_accum_rad(iCell) + 1.0

        enddo
 
        if(allocated(pressure) ) deallocate(pressure )
        if(allocated(pressure2)) deallocate(pressure2)
        if(allocated(field_in) ) deallocate(field_in)
        if(allocated(press_in) ) deallocate(press_in)
        if(allocated(field_interp)) deallocate(field_interp)
        if(allocated(press_interp)) deallocate(press_interp)

   
    end subroutine radiation_diagnostics_update

    !-----------------------------------------------------------------------
    !  routine radiation_diagnostics_compute
    !
    !> \brief Compute isobaric diagnostic before model output is written
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  Compute isobaric diagnostic before model output is written. Code called
    !>  from here was previously in mpas_atm_interp_diagnostics.F.
    !
    !-----------------------------------------------------------------------
    subroutine radiation_diagnostics_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none

        logical :: need_any_diags

        need_any_diags = .false.

        need_any_diags = need_any_diags
        need_rthratenlw_isobaric = MPAS_field_will_be_written('rthratenlw_isobaric')
        need_any_diags = need_any_diags .or. need_rthratenlw_isobaric
        need_rthratensw_isobaric = MPAS_field_will_be_written('rthratensw_isobaric')
        need_any_diags = need_any_diags .or. need_rthratensw_isobaric
        need_acswdnb   = MPAS_field_will_be_written('acswdnb') 
        need_any_diags = need_any_diags .or. need_acswdnb
        need_acswdnbc  = MPAS_field_will_be_written('acswdnbc')
        need_any_diags = need_any_diags .or. need_acswdnbc 
        need_acswdnt   = MPAS_field_will_be_written('acswdnt') 
        need_any_diags = need_any_diags .or.need_acswdnt  
        need_acswdntc  = MPAS_field_will_be_written('acswdntc')
        need_any_diags = need_any_diags .or.need_acswdntc  
        need_acswupb   = MPAS_field_will_be_written('acswupb') 
        need_any_diags = need_any_diags .or.need_acswupb  
        need_acswupbc  = MPAS_field_will_be_written('acswupbc')
        need_any_diags = need_any_diags .or.need_acswupbc
        need_acswupt   = MPAS_field_will_be_written('acswupt') 
        need_any_diags = need_any_diags .or. need_acswupt
        need_acswuptc  = MPAS_field_will_be_written('acswuptc')
        need_any_diags = need_any_diags .or. need_acswuptc 
        need_aclwdnb   = MPAS_field_will_be_written('aclwdnb') 
        need_any_diags = need_any_diags .or. need_aclwdnb 
        need_aclwdnbc  = MPAS_field_will_be_written('aclwdnbc')
        need_any_diags = need_any_diags .or. need_aclwdnbc 
        need_aclwdnt   = MPAS_field_will_be_written('aclwdnt') 
        need_any_diags = need_any_diags .or. need_aclwdnt 
        need_aclwdntc  = MPAS_field_will_be_written('aclwdntc')
        need_any_diags = need_any_diags .or. need_aclwdntc 
        need_aclwupb   = MPAS_field_will_be_written('aclwupb') 
        need_any_diags = need_any_diags .or. need_aclwupb 
        need_aclwupbc  = MPAS_field_will_be_written('aclwupbc')
        need_any_diags = need_any_diags .or. need_aclwupbc 
        need_aclwupt   = MPAS_field_will_be_written('aclwupt')
        need_any_diags = need_any_diags .or. need_aclwupt 
        need_aclwuptc  = MPAS_field_will_be_written('aclwuptc')
        need_any_diags = need_any_diags .or. need_aclwuptc 
        if (need_any_diags) then
            call interp_diagnostics(mesh, state, diag,diag_physics, tend_physics)
        end if
   
    end subroutine radiation_diagnostics_compute
    !-----------------------------------------------------------------------
    !  routine radiation_diagnostics_reset
    !
    !> \brief Reset diagnostic after it has been written
    !> \author 
    !> \date   
    !> \details
    !>  Reset diagnostic after it has been written
    !
    !-----------------------------------------------------------------------
    subroutine radiation_diagnostics_reset()
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        ! Called immediately after diagnostics have been written
        ! Reset counters, accumulations, etc.

        implicit none
        real (kind=RKIND), dimension(:,:), pointer :: rthratenlw_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rthratensw_isobaric
	
        real (kind=RKIND), dimension(:)  , pointer :: acswdnb,acswdnbc,acswdnt,acswdntc,  &
                                                      acswupb,acswupbc,acswupt,acswuptc,  &
                                                      aclwdnb,aclwdnbc,aclwdnt,aclwdntc,  &
                                                      aclwupb,aclwupbc,aclwupt,aclwuptc

	real (kind=RKIND), dimension(:)  , pointer :: nsteps_accum_rad
        logical :: need_any_diags

        need_any_diags = .false.
        need_rthratenlw_isobaric = MPAS_field_will_be_written('rthratenlw_isobaric')
        need_rthratensw_isobaric = MPAS_field_will_be_written('rthratensw_isobaric')
        need_acswdnb = MPAS_field_will_be_written('acswdnb') 
        need_acswdnbc= MPAS_field_will_be_written('acswdnbc')
        need_acswdnt = MPAS_field_will_be_written('acswdnt') 
        need_acswdntc= MPAS_field_will_be_written('acswdntc')
        need_acswupb = MPAS_field_will_be_written('acswupb') 
        need_acswupbc= MPAS_field_will_be_written('acswupbc')
        need_acswupt = MPAS_field_will_be_written('acswupt') 
        need_acswuptc= MPAS_field_will_be_written('acswuptc')
        need_aclwdnb = MPAS_field_will_be_written('aclwdnb') 
        need_aclwdnbc= MPAS_field_will_be_written('aclwdnbc')
        need_aclwdnt = MPAS_field_will_be_written('aclwdnt') 
        need_aclwdntc= MPAS_field_will_be_written('aclwdntc')
        need_aclwupb = MPAS_field_will_be_written('aclwupb') 
        need_aclwupbc= MPAS_field_will_be_written('aclwupbc')
        need_aclwupt = MPAS_field_will_be_written('aclwupt') 
        need_aclwuptc= MPAS_field_will_be_written('aclwuptc')
                                                   
        need_any_diags = need_any_diags .or. need_rthratenlw_isobaric
        need_any_diags = need_any_diags .or. need_rthratensw_isobaric
        need_any_diags = need_any_diags .or. need_acswdnb 
        need_any_diags = need_any_diags .or. need_acswdnbc
        need_any_diags = need_any_diags .or. need_acswdnt 
        need_any_diags = need_any_diags .or. need_acswdntc
        need_any_diags = need_any_diags .or. need_acswupb 
        need_any_diags = need_any_diags .or. need_acswupbc
        need_any_diags = need_any_diags .or. need_acswupt 
        need_any_diags = need_any_diags .or. need_acswuptc
        need_any_diags = need_any_diags .or. need_aclwdnb 
        need_any_diags = need_any_diags .or. need_aclwdnbc
        need_any_diags = need_any_diags .or. need_aclwdnt 
        need_any_diags = need_any_diags .or. need_aclwdntc
        need_any_diags = need_any_diags .or. need_aclwupb 
        need_any_diags = need_any_diags .or. need_aclwupbc
        need_any_diags = need_any_diags .or. need_aclwupt 
        need_any_diags = need_any_diags .or. need_aclwuptc

        if(need_any_diags)then
           call mpas_pool_get_array(diag, 'rthratenlw_isobaric', rthratenlw_isobaric)
           call mpas_pool_get_array(diag, 'rthratensw_isobaric', rthratensw_isobaric)

        !short-wave radiation:
           call mpas_pool_get_array(diag_physics,'acswdnb'   , acswdnb   )
           call mpas_pool_get_array(diag_physics,'acswdnbc'  , acswdnbc  )
           call mpas_pool_get_array(diag_physics,'acswdnt'   , acswdnt   )
           call mpas_pool_get_array(diag_physics,'acswdntc'  , acswdntc  )
           call mpas_pool_get_array(diag_physics,'acswupb'   , acswupb   )
           call mpas_pool_get_array(diag_physics,'acswupbc'  , acswupbc  )
           call mpas_pool_get_array(diag_physics,'acswupt'   , acswupt   )
           call mpas_pool_get_array(diag_physics,'acswuptc'  , acswuptc  )

           !long-wave radiation:
           call mpas_pool_get_array(diag_physics,'aclwdnb'   , aclwdnb   )
           call mpas_pool_get_array(diag_physics,'aclwdnbc'  , aclwdnbc  )
           call mpas_pool_get_array(diag_physics,'aclwdnt'   , aclwdnt   )
           call mpas_pool_get_array(diag_physics,'aclwdntc'  , aclwdntc  )
           call mpas_pool_get_array(diag_physics,'aclwupb'   , aclwupb   )
           call mpas_pool_get_array(diag_physics,'aclwupbc'  , aclwupbc  )
           call mpas_pool_get_array(diag_physics,'aclwupt'   , aclwupt   )
           call mpas_pool_get_array(diag_physics,'aclwuptc'  , aclwuptc  )
           call mpas_pool_get_array(diag_physics,'nsteps_accum_rad',nsteps_accum_rad)

           rthratenlw_isobaric=0.0
           rthratensw_isobaric=0.0
           acswdnb =0.0
           acswdnbc=0.0
           acswdnt =0.0
           acswdntc=0.0
           acswupb =0.0
           acswupbc=0.0
           acswupt =0.0
           acswuptc=0.0


           aclwdnb =0.0
           aclwdnbc=0.0
           aclwdnt =0.0
           aclwdntc=0.0
           aclwupb =0.0
           aclwupbc=0.0
           aclwupt =0.0
           aclwuptc=0.0
   
           nsteps_accum_rad=0.0
        end if 

       call reset_diagnostics(diag_physics )

    end subroutine radiation_diagnostics_reset

   !==================================================================================================
    subroutine reset_diagnostics( diag_physics)
   !==================================================================================================

       !input arguments:
       type (mpas_pool_type), intent(in) :: diag_physics
       ! Called immediately after diagnostics have been written
       ! Reset counters, accumulations, etc.


    end subroutine reset_diagnostics


    !-----------------------------------------------------------------------
    !  routine radiation_diagnostics_cleanup
    !
    !> \brief Finalizes diagnostic
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Finalizes diagnostic
    !
    !-----------------------------------------------------------------------
    subroutine radiation_diagnostics_cleanup()

        implicit none

        ! Deallocate scratch arrays, etc.
   
    end subroutine radiation_diagnostics_cleanup





   !==================================================================================================
    subroutine interp_diagnostics(mesh, state, diag, diag_physics, tend_physics)
   !==================================================================================================

       !input arguments:
        type (mpas_pool_type), intent(in) :: mesh
        type (mpas_pool_type), intent(in) :: state
        type (mpas_pool_type), intent(in) :: diag_physics
        type (mpas_pool_type), intent(in) :: tend_physics
       !inout arguments:
        type (mpas_pool_type), intent(inout) :: diag
       
       !local variables:
        integer          :: iCell,iVert,iVertD,k
        integer, pointer :: nCells, nVertLevels, nVertices, nIsoLevelsT

        real (kind=RKIND), dimension(:,:), pointer :: rthratenlw_isobaric
        real (kind=RKIND), dimension(:,:), pointer :: rthratensw_isobaric
	
        real (kind=RKIND), dimension(:)  ,pointer :: acswdnb,acswdnbc,acswdnt,acswdntc,  &
                                                     acswupb,acswupbc,acswupt,acswuptc,  &
                                                     aclwdnb,aclwdnbc,aclwdnt,aclwdntc,  &
                                                     aclwupb,aclwupbc,aclwupt,aclwuptc

	real (kind=RKIND), dimension(:)  , pointer:: nsteps_accum_rad
                                                      
               

       !--------------------------------------------------------------------------------------------------
       
       ! call mpas_log_write('')
       ! call mpas_log_write('--- enter subroutine interp_diagnostics:')
       
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
        call mpas_pool_get_dimension(mesh, 'nIsoLevelsT', nIsoLevelsT)
       
        call mpas_pool_get_array(diag, 'rthratenlw_isobaric', rthratenlw_isobaric)
        call mpas_pool_get_array(diag, 'rthratensw_isobaric', rthratensw_isobaric)
      
        !short-wave radiation:
        call mpas_pool_get_array(diag_physics,'acswdnb'   , acswdnb   )
        call mpas_pool_get_array(diag_physics,'acswdnbc'  , acswdnbc  )
        call mpas_pool_get_array(diag_physics,'acswdnt'   , acswdnt   )
        call mpas_pool_get_array(diag_physics,'acswdntc'  , acswdntc  )
        call mpas_pool_get_array(diag_physics,'acswupb'   , acswupb   )
        call mpas_pool_get_array(diag_physics,'acswupbc'  , acswupbc  )
        call mpas_pool_get_array(diag_physics,'acswupt'   , acswupt   )
        call mpas_pool_get_array(diag_physics,'acswuptc'  , acswuptc  )

        !long-wave radiation:
        call mpas_pool_get_array(diag_physics,'aclwdnb'   , aclwdnb   )
        call mpas_pool_get_array(diag_physics,'aclwdnbc'  , aclwdnbc  )
        call mpas_pool_get_array(diag_physics,'aclwdnt'   , aclwdnt   )
        call mpas_pool_get_array(diag_physics,'aclwdntc'  , aclwdntc  )
        call mpas_pool_get_array(diag_physics,'aclwupb'   , aclwupb   )
        call mpas_pool_get_array(diag_physics,'aclwupbc'  , aclwupbc  )
        call mpas_pool_get_array(diag_physics,'aclwupt'   , aclwupt   )
        call mpas_pool_get_array(diag_physics,'aclwuptc'  , aclwuptc  )

        call mpas_pool_get_array(diag_physics,'nsteps_accum_rad',nsteps_accum_rad)
       
                      
        !!!!!!!!!!! Additional temperature levels for vortex tracking !!!!!!!!!!!
        do iCell=1,nCells          
           if(nsteps_accum_rad(iCell) /= 0._RKIND)then
              do iVert=1,nIsoLevelsT
                  rthratenlw_isobaric(iVert,iCell) = rthratenlw_isobaric(iVert,iCell)/nsteps_accum_rad(iCell)
                  rthratensw_isobaric(iVert,iCell) = rthratensw_isobaric(iVert,iCell)/nsteps_accum_rad(iCell) 
              end do
           endif
        end do   
        do iCell=1,nCells
          if(nsteps_accum_rad(iCell) /= 0._RKIND)then
            !short-wave radiation:
            acswdnb(iCell)  = acswdnb (iCell)/nsteps_accum_rad(iCell) 
            acswdnbc(iCell) = acswdnbc(iCell)/nsteps_accum_rad(iCell) 
            acswdnt(iCell)  = acswdnt (iCell)/nsteps_accum_rad(iCell) 
            acswdntc(iCell) = acswdntc(iCell)/nsteps_accum_rad(iCell) 
            acswupb(iCell)  = acswupb (iCell)/nsteps_accum_rad(iCell) 
            acswupbc(iCell) = acswupbc(iCell)/nsteps_accum_rad(iCell) 
            acswupt(iCell)  = acswupt (iCell)/nsteps_accum_rad(iCell) 
            acswuptc(iCell) = acswuptc(iCell)/nsteps_accum_rad(iCell)
            !long-wave radiation:
            aclwdnb(iCell)  = aclwdnb (iCell)/nsteps_accum_rad(iCell) 
            aclwdnbc(iCell) = aclwdnbc(iCell)/nsteps_accum_rad(iCell) 
            aclwdnt(iCell)  = aclwdnt (iCell)/nsteps_accum_rad(iCell) 
            aclwdntc(iCell) = aclwdntc(iCell)/nsteps_accum_rad(iCell) 
            aclwupb(iCell)  = aclwupb (iCell)/nsteps_accum_rad(iCell)
            aclwupbc(iCell) = aclwupbc(iCell)/nsteps_accum_rad(iCell) 
            aclwupt(iCell)  = aclwupt (iCell)/nsteps_accum_rad(iCell)
            aclwuptc(iCell) = aclwuptc(iCell)/nsteps_accum_rad(iCell) 
           endif
        end do

       
    end subroutine interp_diagnostics


 

   
end module mpas_radiation_diagnostics
