! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_sfclayer_diagnostics

    use mpas_dmpar
    use mpas_kind_types
    use mpas_derived_types
    use mpas_pool_routines
    use mpas_constants
    use mpas_log, only : mpas_log_write


    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_pool_type), pointer :: diag_physics
    type (MPAS_pool_type), pointer :: sfc_input    
    type (MPAS_pool_type), pointer :: tend_physics
    type (MPAS_clock_type), pointer :: clock

    public :: sfclayer_diagnostics_setup, &
              sfclayer_diagnostics_update, &
              sfclayer_diagnostics_compute,&
              sfclayer_diagnostics_reset,&
              sfclayer_diagnostics_cleanup

    private
    real(kind=RKIND),pointer:: bucket_sfcpbldt 
    logical :: need_alh      , &
               need_ahfx     , &   
               need_sfc_runoff    , &
               need_udr_runoff    , &
               need_sms_total     , &
               need_smois_lys01   , &
               need_smois_lys02   , &
               need_smois_lys03   , &
               need_smois_lys04   , &
               need_tslb_lys01    , &
               need_tslb_lys02    , &
               need_tslb_lys03    , &
               need_tslb_lys04    , &
               need_t02mt         , &
               need_q02mt         , &
               need_v10mt         , &
               need_u10mt

    real (kind=RKIND),pointer:: config_bucket_sfcpbldt
    character(len=StrKIND),pointer:: config_bucket_update

    contains


    !-----------------------------------------------------------------------
    !  routine sfclayer_diagnostics_setup
    !
    !> \brief Set up the isobaric diagnostics module
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  This routine sets up the isobaric diagnostics module, principally by
    !>  saving pointers to pools that are used in the computation of diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine sfclayer_diagnostics_setup(all_pools, simulation_clock,domain)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock
        type(domain_type),intent(in):: domain
        !local pointers:
        integer,pointer:: nCellsSolve 
        real(kind=RKIND),dimension(:),pointer:: alh,ahfx, sfc_runoff    , &
                                                          udr_runoff    , &
                                                          sms_total     , &
                                                          smois_lys01   , &
                                                          smois_lys02   , &
                                                          smois_lys03   , &
                                                          smois_lys04   , &
                                                          tslb_lys01    , &
                                                          tslb_lys02    , &
                                                          tslb_lys03    , &
                                                          tslb_lys04    , &
                                                          t02mt         , &
                                                          q02mt         , &
                                                          v10mt         , &
                                                          u10mt
        real(kind=RKIND),dimension(:),pointer:: nsteps_accum_sfcpbl

 
        !local variables and arrays:
        integer:: iCell

        clock => simulation_clock
        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
        call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)
        call mpas_pool_get_subpool(all_pools, 'sfc_input',sfc_input)
        call mpas_pool_get_subpool(all_pools, 'tend_physics', tend_physics)
        call mpas_pool_get_config(domain%configs,'config_bucket_sfcpbldt'      ,config_bucket_sfcpbldt      )
        call mpas_pool_get_config(domain%configs,'config_bucket_update'    ,config_bucket_update    )

        call mpas_pool_get_array(diag,'alh'           ,alh   )
        call mpas_pool_get_array(diag,'ahfx'          ,ahfx   )
        call mpas_pool_get_array(diag,'sfc_runoff'    ,sfc_runoff  ) 
        call mpas_pool_get_array(diag,'udr_runoff'    ,udr_runoff  ) 
        call mpas_pool_get_array(diag,'sms_total'     ,sms_total   ) 
        call mpas_pool_get_array(diag,'smois_lys01'   ,smois_lys01  )
        call mpas_pool_get_array(diag,'smois_lys02'   ,smois_lys02  )
        call mpas_pool_get_array(diag,'smois_lys03'   ,smois_lys03  )
        call mpas_pool_get_array(diag,'smois_lys04'   ,smois_lys04  )
        call mpas_pool_get_array(diag,'tslb_lys01'    ,tslb_lys01  ) 
        call mpas_pool_get_array(diag,'tslb_lys02'    ,tslb_lys02  ) 
        call mpas_pool_get_array(diag,'tslb_lys03'    ,tslb_lys03  ) 
        call mpas_pool_get_array(diag,'tslb_lys04'    ,tslb_lys04  )
        call mpas_pool_get_array(diag,'t02mt'         ,t02mt )
        call mpas_pool_get_array(diag,'q02mt'         ,q02mt )
        call mpas_pool_get_array(diag,'v10mt'         ,v10mt )
        call mpas_pool_get_array(diag,'u10mt'         ,u10mt )
        call mpas_pool_get_array(diag,'nsteps_accum_sfcpbl',nsteps_accum_sfcpbl)

       !surface layer boundary layer:
       alh(:)  = 0.0
       ahfx(:) = 0.0
       sfc_runoff  (:) = 0.0 
       udr_runoff  (:) = 0.0  
       sms_total   (:) = 0.0  
       smois_lys01 (:) = 0.0 
       smois_lys02 (:) = 0.0 
       smois_lys03 (:) = 0.0 
       smois_lys04 (:) = 0.0 
       tslb_lys01  (:) = 0.0 
       tslb_lys02  (:) = 0.0 
       tslb_lys03  (:) = 0.0 
       tslb_lys04  (:) = 0.0
       t02mt(:) = 0.0
       q02mt(:) = 0.0
       v10mt(:) = 0.0
       u10mt(:) = 0.0
       nsteps_accum_sfcpbl(:) =0.0 
    end subroutine sfclayer_diagnostics_setup





    !-----------------------------------------------------------------------
    !  routine sfclayer_diagnostics_update
    !
    !> \brief Handle diagnostic calculation needed after each timestep
    !> \author 
    !> \date   
    !> \details
    !>  Handle diagnostic calculation needed after each timestep
    !
    !-----------------------------------------------------------------------
    subroutine sfclayer_diagnostics_update()
    
        implicit none

        real(kind=RKIND),dimension(:),pointer:: lh,hfx,sfcrunoff       ,& 
                                                       smstot          ,&
                                                       udrunoff
        real(kind=RKIND),dimension(:,:),pointer:: smois
        real(kind=RKIND),dimension(:,:),pointer:: tslb  
        real(kind=RKIND),dimension(:)  ,pointer:: q2,t2m,u10,v10

        real(kind=RKIND),dimension(:),pointer:: alh,ahfx,sfc_runoff    , &
                                                         udr_runoff    , &
                                                         sms_total     , &
                                                         smois_lys01   , &
                                                         smois_lys02   , &
                                                         smois_lys03   , &
                                                         smois_lys04   , &
                                                         tslb_lys01    , &
                                                         tslb_lys02    , &
                                                         tslb_lys03    , &
                                                         tslb_lys04    , &
                                                         t02mt         , &
                                                         q02mt         , &
                                                         v10mt         , &
                                                         u10mt
        real(kind=RKIND),dimension(:),pointer:: nsteps_accum_sfcpbl

       if(config_bucket_update /= 'none' .and. config_bucket_sfcpbldt .gt. 0._RKIND) then
           call mpas_pool_get_array(diag_physics,'lh'   , lh   )
           call mpas_pool_get_array(diag_physics,'hfx'  , hfx   )
           call mpas_pool_get_array(diag_physics,'sfcrunoff' ,sfcrunoff )
           call mpas_pool_get_array(diag_physics,'smstot'    ,smstot    )
           call mpas_pool_get_array(diag_physics,'udrunoff'  ,udrunoff  )
           call mpas_pool_get_array(sfc_input,'smois'     ,smois     )
           call mpas_pool_get_array(sfc_input,'tslb'      ,tslb      )
           call mpas_pool_get_array(diag_physics,'q2'    ,q2    )
           call mpas_pool_get_array(diag_physics,'t2m'   ,t2m   )
           call mpas_pool_get_array(diag_physics,'u10'   ,u10   )
           call mpas_pool_get_array(diag_physics,'v10'   ,v10   )

           !surface layer pbl :
           call mpas_pool_get_array(diag,'alh'   , alh   )
           call mpas_pool_get_array(diag,'ahfx'  , ahfx   )
           call mpas_pool_get_array(diag,'sfc_runoff'    ,sfc_runoff  ) 
           call mpas_pool_get_array(diag,'udr_runoff'    ,udr_runoff  ) 
           call mpas_pool_get_array(diag,'sms_total'     ,sms_total   ) 
           call mpas_pool_get_array(diag,'smois_lys01'   ,smois_lys01  )
           call mpas_pool_get_array(diag,'smois_lys02'   ,smois_lys02  )
           call mpas_pool_get_array(diag,'smois_lys03'   ,smois_lys03  )
           call mpas_pool_get_array(diag,'smois_lys04'   ,smois_lys04  )
           call mpas_pool_get_array(diag,'tslb_lys01'    ,tslb_lys01  ) 
           call mpas_pool_get_array(diag,'tslb_lys02'    ,tslb_lys02  ) 
           call mpas_pool_get_array(diag,'tslb_lys03'    ,tslb_lys03  ) 
           call mpas_pool_get_array(diag,'tslb_lys04'    ,tslb_lys04  )
           call mpas_pool_get_array(diag,'t02mt'         ,t02mt  )
           call mpas_pool_get_array(diag,'q02mt'         ,q02mt  )
           call mpas_pool_get_array(diag,'v10mt'         ,v10mt  )
           call mpas_pool_get_array(diag,'u10mt'         ,u10mt  )
           call mpas_pool_get_array(diag,'nsteps_accum_sfcpbl',nsteps_accum_sfcpbl)

           ! Called at the end of every timestep
           ! Update extrema, accumulations, etc.
           !surface layer boundary layer:
           alh(:)  = alh (:) + lh (:)
           ahfx(:) = ahfx(:) + hfx(:)
           sfc_runoff(:) = sfc_runoff(:)  + sfcrunoff(:)
           udr_runoff(:) = udr_runoff(:)  + udrunoff(:)
           sms_total(:)  = sms_total (:)  + smstot (:)
   
           smois_lys01(:)= smois_lys01(:) + smois(1,:)
           smois_lys02(:)= smois_lys02(:) + smois(2,:)
           smois_lys03(:)= smois_lys03(:) + smois(3,:)
           smois_lys04(:)= smois_lys04(:) + smois(4,:)
   
           tslb_lys01(:) = tslb_lys01(:)  + tslb(1,:)
           tslb_lys02(:) = tslb_lys02(:)  + tslb(2,:)
           tslb_lys03(:) = tslb_lys03(:)  + tslb(3,:)
           tslb_lys04(:) = tslb_lys04(:)  + tslb(4,:)


           t02mt(:) = t02mt(:) + t2m(:) 
           q02mt(:) = q02mt(:) + q2(:)
           v10mt(:) = v10mt(:) + v10(:)
           u10mt(:) = u10mt(:) + u10(:)
           nsteps_accum_sfcpbl(:) = nsteps_accum_sfcpbl(:) + 1.0
        end if 

    end subroutine sfclayer_diagnostics_update


    !-----------------------------------------------------------------------
    !  routine sfclayer_diagnostics_compute
    !
    !> \brief Compute isobaric diagnostic before model output is written
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  Compute isobaric diagnostic before model output is written. Code called
    !>  from here was previously in mpas_atm_interp_diagnostics.F.
    !
    !-----------------------------------------------------------------------
    subroutine sfclayer_diagnostics_compute()
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none
        logical :: need_any_diags


        need_any_diags = .false.
        need_any_diags = need_any_diags
        need_alh = MPAS_field_will_be_written('alh')
        need_any_diags = need_any_diags .or. need_alh
        need_ahfx = MPAS_field_will_be_written('ahfx')
        need_any_diags = need_any_diags .or. need_ahfx
        need_sfc_runoff  = MPAS_field_will_be_written('sfc_runoff')
        need_any_diags = need_any_diags .or. need_sfc_runoff
        need_udr_runoff  = MPAS_field_will_be_written('udr_runoff')
        need_any_diags = need_any_diags .or. need_udr_runoff 
        need_sms_total   = MPAS_field_will_be_written('sms_total')
        need_any_diags = need_any_diags .or. need_sms_total 
        need_smois_lys01 = MPAS_field_will_be_written('smois_lys01')
        need_any_diags = need_any_diags .or. need_smois_lys01
        need_smois_lys02 = MPAS_field_will_be_written('smois_lys02')
        need_any_diags = need_any_diags .or. need_smois_lys02 
        need_smois_lys03 = MPAS_field_will_be_written('smois_lys03')
        need_any_diags = need_any_diags .or. need_smois_lys03
        need_smois_lys04 = MPAS_field_will_be_written('smois_lys04')
        need_any_diags = need_any_diags .or. need_smois_lys04 
        need_tslb_lys01  = MPAS_field_will_be_written('tslb_lys01')
        need_any_diags = need_any_diags .or.  need_tslb_lys01
        need_tslb_lys02  = MPAS_field_will_be_written('tslb_lys02')
        need_any_diags = need_any_diags .or. need_tslb_lys02
        need_tslb_lys03  = MPAS_field_will_be_written('tslb_lys03')
        need_any_diags = need_any_diags .or. need_tslb_lys03 
        need_tslb_lys04  = MPAS_field_will_be_written('tslb_lys04')
        need_any_diags = need_any_diags .or. need_tslb_lys04
        need_t02mt = MPAS_field_will_be_written('t02mt')
        need_any_diags = need_any_diags .or. need_t02mt
        need_q02mt = MPAS_field_will_be_written('q02mt')
        need_any_diags = need_any_diags .or. need_q02mt
        need_v10mt = MPAS_field_will_be_written('v10mt')
        need_any_diags = need_any_diags .or. need_v10mt
        need_u10mt = MPAS_field_will_be_written('u10mt')
        need_any_diags = need_any_diags .or. need_u10mt

        if (need_any_diags) then
            call interp_diagnostics(mesh, state, 1, diag,diag_physics, tend_physics)
        end if
   
    end subroutine sfclayer_diagnostics_compute
    !-----------------------------------------------------------------------
    !  routine sfclayer_diagnostics_reset
    !
    !> \brief Reset diagnostic after it has been written
    !> \author 
    !> \date   
    !> \details
    !>  Reset diagnostic after it has been written
    !
    !-----------------------------------------------------------------------
    subroutine sfclayer_diagnostics_reset()
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none

        real(kind=RKIND),dimension(:),pointer:: alh,ahfx,sfc_runoff    , &
                                                         udr_runoff    , &
                                                         sms_total     , &
                                                         smois_lys01   , &
                                                         smois_lys02   , &
                                                         smois_lys03   , &
                                                         smois_lys04   , &
                                                         tslb_lys01    , &
                                                         tslb_lys02    , &
                                                         tslb_lys03    , &
                                                         tslb_lys04    , &
                                                         t02mt    , &
                                                         q02mt    , &
                                                         v10mt    , &
                                                         u10mt
        real(kind=RKIND),dimension(:),pointer:: nsteps_accum_sfcpbl
        logical :: need_any_diags

        need_any_diags = .false.

        need_any_diags = need_any_diags
        need_alh = MPAS_field_will_be_written('alh')
        need_any_diags = need_any_diags .or. need_alh
        need_ahfx = MPAS_field_will_be_written('ahfx')
        need_any_diags = need_any_diags .or. need_ahfx
        need_sfc_runoff  = MPAS_field_will_be_written('sfc_runoff')
        need_any_diags = need_any_diags .or. need_sfc_runoff
        need_udr_runoff  = MPAS_field_will_be_written('udr_runoff')
        need_any_diags = need_any_diags .or. need_udr_runoff 
        need_sms_total   = MPAS_field_will_be_written('sms_total')
        need_any_diags = need_any_diags .or. need_sms_total 
        need_smois_lys01 = MPAS_field_will_be_written('smois_lys01')
        need_any_diags = need_any_diags .or. need_smois_lys01
        need_smois_lys02 = MPAS_field_will_be_written('smois_lys02')
        need_any_diags = need_any_diags .or. need_smois_lys02 
        need_smois_lys03 = MPAS_field_will_be_written('smois_lys03')
        need_any_diags = need_any_diags .or. need_smois_lys03
        need_smois_lys04 = MPAS_field_will_be_written('smois_lys04')
        need_any_diags = need_any_diags .or. need_smois_lys04 
        need_tslb_lys01  = MPAS_field_will_be_written('tslb_lys01')
        need_any_diags = need_any_diags .or.  need_tslb_lys01
        need_tslb_lys02  = MPAS_field_will_be_written('tslb_lys02')
        need_any_diags = need_any_diags .or. need_tslb_lys02
        need_tslb_lys03  = MPAS_field_will_be_written('tslb_lys03')
        need_any_diags = need_any_diags .or. need_tslb_lys03 
        need_tslb_lys04  = MPAS_field_will_be_written('tslb_lys04')
        need_any_diags = need_any_diags .or. need_tslb_lys04
        need_t02mt = MPAS_field_will_be_written('t02mt')
        need_any_diags = need_any_diags .or. need_t02mt
        need_q02mt = MPAS_field_will_be_written('q02mt')
        need_any_diags = need_any_diags .or. need_q02mt
        need_v10mt = MPAS_field_will_be_written('v10mt')
        need_any_diags = need_any_diags .or. need_v10mt
        need_u10mt = MPAS_field_will_be_written('u10mt')
        need_any_diags = need_any_diags .or. need_u10mt
 
        if(config_bucket_update /= 'none' .and. config_bucket_sfcpbldt .gt. 0._RKIND) then
           !surface layer pbl :
           if(need_any_diags)then
              call mpas_pool_get_array(diag,'alh'   , alh   )
              call mpas_pool_get_array(diag,'ahfx'  , ahfx   )
              call mpas_pool_get_array(diag,'sfc_runoff'    ,sfc_runoff  ) 
              call mpas_pool_get_array(diag,'udr_runoff'    ,udr_runoff  ) 
              call mpas_pool_get_array(diag,'sms_total'     ,sms_total   ) 
              call mpas_pool_get_array(diag,'smois_lys01'   ,smois_lys01  )
              call mpas_pool_get_array(diag,'smois_lys02'   ,smois_lys02  )
              call mpas_pool_get_array(diag,'smois_lys03'   ,smois_lys03  )
              call mpas_pool_get_array(diag,'smois_lys04'   ,smois_lys04  )
              call mpas_pool_get_array(diag,'tslb_lys01'    ,tslb_lys01  ) 
              call mpas_pool_get_array(diag,'tslb_lys02'    ,tslb_lys02  ) 
              call mpas_pool_get_array(diag,'tslb_lys03'    ,tslb_lys03  ) 
              call mpas_pool_get_array(diag,'tslb_lys04'    ,tslb_lys04  )
              call mpas_pool_get_array(diag,'t02mt'         ,t02mt )
              call mpas_pool_get_array(diag,'q02mt'         ,q02mt )
              call mpas_pool_get_array(diag,'v10mt'         ,v10mt )
              call mpas_pool_get_array(diag,'u10mt'         ,u10mt )
              call mpas_pool_get_array(diag,'nsteps_accum_sfcpbl',nsteps_accum_sfcpbl)
              ! Called immediately after diagnostics have been written
              ! Reset counters, accumulations, etc.
              alh=0.0
              ahfx=0.0
              sfc_runoff=0.0
              udr_runoff=0.0
              sms_total=0.0
              smois_lys01=0.0
              smois_lys02=0.0
              smois_lys03=0.0
              smois_lys04=0.0
              tslb_lys01=0.0
              tslb_lys02=0.0
              tslb_lys03=0.0
              tslb_lys04=0.0
              t02mt=0.0
              q02mt=0.0
              v10mt=0.0
              u10mt=0.0
              nsteps_accum_sfcpbl=0.0
           end if 
        end if 

       !call reset_diagnostics(diag_physics )

    end subroutine sfclayer_diagnostics_reset

   !==================================================================================================
    subroutine reset_diagnostics( diag_physics)
   !==================================================================================================

       !input arguments:
	type (mpas_pool_type), intent(in) :: diag_physics

    end subroutine reset_diagnostics


    !-----------------------------------------------------------------------
    !  routine sfclayer_diagnostics_cleanup
    !
    !> \brief Finalizes diagnostic
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Finalizes diagnostic
    !
    !-----------------------------------------------------------------------
    subroutine sfclayer_diagnostics_cleanup()

        implicit none

        ! Deallocate scratch arrays, etc.
   
    end subroutine sfclayer_diagnostics_cleanup





   !==================================================================================================
    subroutine interp_diagnostics(mesh, state, time_lev, diag, diag_physics, tend_physics)
   !==================================================================================================

       !input arguments:
        type (mpas_pool_type), intent(in)  :: mesh
        type (mpas_pool_type), intent(in) :: state
	type (mpas_pool_type), intent(in) :: diag_physics
        type (mpas_pool_type), intent(in) :: tend_physics

        integer, intent(in) :: time_lev              ! which time level to use from state
       
       !inout arguments:
        type (mpas_pool_type), intent(inout) :: diag
       
       !local variables:
        integer :: iCell
        integer, pointer :: nCells


        real(kind=RKIND),dimension(:),pointer:: alh,ahfx, sfc_runoff    , &
                                                          udr_runoff    , &
                                                          sms_total     , &
                                                          smois_lys01   , &
                                                          smois_lys02   , &
                                                          smois_lys03   , &
                                                          smois_lys04   , &
                                                          tslb_lys01    , &
                                                          tslb_lys02    , &
                                                          tslb_lys03    , &
                                                          tslb_lys04    , &
                                                          t02mt    , &
                                                          q02mt    , &
                                                          v10mt    , &
                                                          u10mt
        real (kind=RKIND), dimension(:)  , pointer:: nsteps_accum_sfcpbl
                                                      
               
        !--------------------------------------------------------------------------------------------------
       
       ! call mpas_log_write('')
       ! call mpas_log_write('--- enter subroutine interp_diagnostics:')
       
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)

   
        if(config_bucket_update /= 'none' .and. config_bucket_sfcpbldt .gt. 0._RKIND) then

           call mpas_pool_get_array(diag,'alh'   , alh   )
           call mpas_pool_get_array(diag,'ahfx'  , ahfx  )
           call mpas_pool_get_array(diag,'sfc_runoff'    ,sfc_runoff  ) 
           call mpas_pool_get_array(diag,'udr_runoff'    ,udr_runoff  ) 
           call mpas_pool_get_array(diag,'sms_total'     ,sms_total   ) 
           call mpas_pool_get_array(diag,'smois_lys01'   ,smois_lys01  )
           call mpas_pool_get_array(diag,'smois_lys02'   ,smois_lys02  )
           call mpas_pool_get_array(diag,'smois_lys03'   ,smois_lys03  )
           call mpas_pool_get_array(diag,'smois_lys04'   ,smois_lys04  )
           call mpas_pool_get_array(diag,'tslb_lys01'    ,tslb_lys01  ) 
           call mpas_pool_get_array(diag,'tslb_lys02'    ,tslb_lys02  ) 
           call mpas_pool_get_array(diag,'tslb_lys03'    ,tslb_lys03  ) 
           call mpas_pool_get_array(diag,'tslb_lys04'    ,tslb_lys04  )
           call mpas_pool_get_array(diag,'t02mt'         ,t02mt  )
           call mpas_pool_get_array(diag,'q02mt'         ,q02mt  )
           call mpas_pool_get_array(diag,'v10mt'         ,v10mt  )
           call mpas_pool_get_array(diag,'u10mt'         ,u10mt  )
           call mpas_pool_get_array(diag,'nsteps_accum_sfcpbl',nsteps_accum_sfcpbl)
        endif       
       
        if(config_bucket_update /= 'none' .and. config_bucket_sfcpbldt .gt. 0._RKIND) then
           do iCell=1,nCells
              if(nsteps_accum_sfcpbl(iCell) /= 0._RKIND)then
                !short-wave radiation:
                alh(iCell)         = alh (iCell)/nsteps_accum_sfcpbl(iCell) 
                ahfx(iCell)        = ahfx(iCell)/nsteps_accum_sfcpbl(iCell) 
                sfc_runoff(iCell)  = sfc_runoff(iCell)/nsteps_accum_sfcpbl(iCell)  
                udr_runoff(iCell)  = udr_runoff(iCell)/nsteps_accum_sfcpbl(iCell)  
                sms_total (iCell)   = sms_total (iCell)/nsteps_accum_sfcpbl(iCell)  
                smois_lys01(iCell) = smois_lys01(iCell)/nsteps_accum_sfcpbl(iCell) 
                smois_lys02(iCell) = smois_lys02(iCell)/nsteps_accum_sfcpbl(iCell) 
                smois_lys03(iCell) = smois_lys03(iCell)/nsteps_accum_sfcpbl(iCell) 
                smois_lys04(iCell) = smois_lys04(iCell)/nsteps_accum_sfcpbl(iCell) 
                tslb_lys01(iCell)  = tslb_lys01(iCell)/nsteps_accum_sfcpbl(iCell)  
                tslb_lys02(iCell)  = tslb_lys02(iCell)/nsteps_accum_sfcpbl(iCell)  
                tslb_lys03(iCell)  = tslb_lys03(iCell)/nsteps_accum_sfcpbl(iCell)  
                tslb_lys04(iCell)  = tslb_lys04(iCell)/nsteps_accum_sfcpbl(iCell)  
                t02mt(iCell)  = t02mt(iCell)/nsteps_accum_sfcpbl(iCell)  
                q02mt(iCell)  = q02mt(iCell)/nsteps_accum_sfcpbl(iCell)  
                v10mt(iCell)  = v10mt(iCell)/nsteps_accum_sfcpbl(iCell)  
                u10mt(iCell)  = u10mt(iCell)/nsteps_accum_sfcpbl(iCell)  
               endif
            end do
        endif
       
    end subroutine interp_diagnostics



   
end module mpas_sfclayer_diagnostics
