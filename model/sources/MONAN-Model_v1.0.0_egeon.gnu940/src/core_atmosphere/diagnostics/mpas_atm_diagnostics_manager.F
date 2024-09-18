! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_atm_diagnostics_manager

    private

    public :: mpas_atm_diag_setup, &
              mpas_atm_diag_update, &
              mpas_atm_diag_compute, &
              mpas_atm_diag_reset, &
              mpas_atm_diag_cleanup

    contains


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_setup
    !
    !> \brief Initialize the diagnostics manager and all diagnostics
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Initialize the diagnostics manager and all diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_setup(stream_mgr, configs, structs, clock, dminfo,domain)

        use mpas_atm_diagnostics_utils, only : mpas_atm_diag_utils_init
        use mpas_derived_types, only : MPAS_streamManager_type, MPAS_pool_type, MPAS_clock_type, dm_info,domain_type
        use mpas_diagnostic_template, only : diagnostic_template_setup
        use mpas_isobaric_diagnostics, only : isobaric_diagnostics_setup
        use mpas_cloud_diagnostics, only : cloud_diagnostics_setup
        use mpas_radiation_diagnostics, only : radiation_diagnostics_setup
        use mpas_sfclayer_diagnostics, only : sfclayer_diagnostics_setup
        use mpas_pbl_diagnostics, only : pbl_diagnostics_setup
        use mpas_convective_diagnostics, only : convective_diagnostics_setup
        use mpas_pv_diagnostics, only : pv_diagnostics_setup
        use mpas_soundings, only : soundings_setup

        implicit none

        type (MPAS_streamManager_type), target, intent(inout) :: stream_mgr
        type (MPAS_pool_type), pointer :: configs
        type (MPAS_pool_type), pointer :: structs
        type (MPAS_clock_type), pointer :: clock
        type (dm_info), intent(in) :: dminfo
        type (domain_type),intent(in):: domain

        !
        ! Prepare the diagnostics utilities module for later use by diagnostics
        !
        call mpas_atm_diag_utils_init(stream_mgr)

        call diagnostic_template_setup(configs, structs, clock)
        call isobaric_diagnostics_setup(structs, clock ,domain)
        call cloud_diagnostics_setup(structs, clock)
        call radiation_diagnostics_setup(structs, clock,domain)
        call sfclayer_diagnostics_setup(structs, clock,domain)
        call pbl_diagnostics_setup(configs,structs, clock)
        call convective_diagnostics_setup(configs,structs, clock)
        call pv_diagnostics_setup(structs, clock)
        call soundings_setup(configs, structs, clock, dminfo)

    end subroutine mpas_atm_diag_setup


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_update
    !
    !> \brief Handle diagnostics accumulation at the end of each timestep
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  MPAS_atm_diag_update.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_update()

        use mpas_diagnostic_template, only : diagnostic_template_update
        use mpas_convective_diagnostics, only : convective_diagnostics_update
        use mpas_radiation_diagnostics, only : radiation_diagnostics_update
        use mpas_cloud_diagnostics, only : cloud_diagnostics_update
        use mpas_sfclayer_diagnostics, only : sfclayer_diagnostics_update
        use mpas_pbl_diagnostics, only : pbl_diagnostics_update
        use mpas_isobaric_diagnostics, only : isobaric_diagnostic_update
        implicit none


        call diagnostic_template_update()
	call isobaric_diagnostic_update()
        call convective_diagnostics_update()
        call radiation_diagnostics_update()
        call cloud_diagnostics_update()
        call sfclayer_diagnostics_update()
        call pbl_diagnostics_update()

    end subroutine mpas_atm_diag_update


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_compute
    !
    !> \brief Compute diagnostics before they are written to output streams
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  MPAS_atm_diag_compute.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_compute()

        use mpas_diagnostic_template, only : diagnostic_template_compute
        use mpas_isobaric_diagnostics, only : isobaric_diagnostics_compute
        use mpas_cloud_diagnostics, only : cloud_diagnostics_compute
        use mpas_radiation_diagnostics, only : radiation_diagnostics_compute
        use mpas_sfclayer_diagnostics, only : sfclayer_diagnostics_compute
        use mpas_pbl_diagnostics, only : pbl_diagnostics_compute
        use mpas_convective_diagnostics, only : convective_diagnostics_compute
        use mpas_pv_diagnostics, only : pv_diagnostics_compute
        use mpas_soundings, only : soundings_compute

        implicit none


        call diagnostic_template_compute()
        call isobaric_diagnostics_compute()
        call cloud_diagnostics_compute()
        call radiation_diagnostics_compute()
        call sfclayer_diagnostics_compute()
        call pbl_diagnostics_compute()
        call convective_diagnostics_compute()
        call pv_diagnostics_compute()
        call soundings_compute()

    end subroutine mpas_atm_diag_compute


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_reset
    !
    !> \brief Resets a diagnostic after it has been computed and written
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  MPAS_atm_diag_reset.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_reset()

        use mpas_diagnostic_template, only : diagnostic_template_reset
        use mpas_convective_diagnostics, only : convective_diagnostics_reset
        use mpas_isobaric_diagnostics, only : isobaric_diagnostic_reset
        use mpas_radiation_diagnostics, only : radiation_diagnostics_reset
        use mpas_cloud_diagnostics, only : cloud_diagnostics_reset
        use mpas_sfclayer_diagnostics, only : sfclayer_diagnostics_reset
        use mpas_pbl_diagnostics, only : pbl_diagnostics_reset

        implicit none


        call diagnostic_template_reset()
        call convective_diagnostics_reset()
        call isobaric_diagnostic_reset()
        call radiation_diagnostics_reset()
        call cloud_diagnostics_reset()
        call sfclayer_diagnostics_reset()
        call pbl_diagnostics_reset()

    end subroutine mpas_atm_diag_reset


    !-----------------------------------------------------------------------
    !  routine MPAS_atm_diag_cleanup
    !
    !> \brief Finalizes diagnostics manager and all diagnostics
    !> \author Michael Duda
    !> \date   6 September 2016
    !> \details
    !>  Finalizes diagnostics manager and all diagnostics.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_diag_cleanup()

        use mpas_atm_diagnostics_utils, only : mpas_atm_diag_utils_finalize
        use mpas_diagnostic_template, only : diagnostic_template_cleanup
        use mpas_soundings, only : soundings_cleanup
        use mpas_radiation_diagnostics, only : radiation_diagnostics_cleanup
        use mpas_cloud_diagnostics, only : cloud_diagnostics_cleanup
        use mpas_sfclayer_diagnostics, only : sfclayer_diagnostics_cleanup
        use mpas_pbl_diagnostics, only : pbl_diagnostics_cleanup
        use mpas_convective_diagnostics, only : convective_diagnostic_cleanup
        use mpas_isobaric_diagnostics, only : isobaric_diagnostic_cleanup

        implicit none


        call diagnostic_template_cleanup()
        call soundings_cleanup()
	call radiation_diagnostics_cleanup()	
	call cloud_diagnostics_cleanup()
	call sfclayer_diagnostics_cleanup()
	call pbl_diagnostics_cleanup()
	call convective_diagnostic_cleanup()
	call isobaric_diagnostic_cleanup()
        !
        ! Take care of any needed cleanup in the diagnostics utility module
        !
        call mpas_atm_diag_utils_finalize()

    end subroutine mpas_atm_diag_cleanup

end module mpas_atm_diagnostics_manager
