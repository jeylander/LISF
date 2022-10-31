!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: ISCCPHXG_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! ISCCP IR Land Surface Temperature observations
!
!
!   
! !REVISION HISTORY: 
!  21 Aug 2016: Sujay Kumar, Initial Specification
!  29 Dec 2021: John Eylander, ERDC CHL - used this templated to create
!                  module to handle ISCCP HXG land surface temperature
!                  data
!
module ISCCP_HXG_obsMod
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ISCCPHXG_stobsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ISCCPHXGstobs
!EOP
  type, public :: ISCCPHXGstobsdec

     character*100          :: odir
     character*20           :: data_designation
     character*3             :: release_number
     real                   :: search_radius
     integer                :: mo
     real,    allocatable   :: LST_obs(:,:)
     integer                :: nc, nr
     type(proj_info)        :: proj
     integer, allocatable   :: n11(:)
  end type ISCCPHXGstobsdec

  type(ISCCPHXGstobsdec), allocatable:: ISCCPHXGstobs(:)

contains
  
!BOP
! 
! !ROUTINE: ISCCPHXG_stobsinit
! \label{ISCCPHXG_stobsinit}
! 
! !INTERFACE: 
  subroutine ISCCPHXG_stobsinit()
! !USES: 
    use ESMF
    use LDT_coreMod,    only : LDT_rc, LDT_config
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod, only : LDT_clock, LDT_calendar
    use LDT_logMod,     only : LDT_verify, LDT_logunit

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading NASASMAP soil moisture data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(ISCCPHXGstobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
        "ISCCP HXG skin temperature data directory:", rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, ISCCPHXGstobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'ISCCP HXG skin temperature data directory: not defined')
    enddo


    call ESMF_ConfigFindLabel(LDT_config, &
        "ISCCP skin temperature data designation:", rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, &
            ISCCPHXGstobs(n)%data_designation, &
            rc=status)
       call LDT_verify(status, &
            'ISCCP skin temperature data designation: not defined')
    enddo

    

    do n=1,LDT_rc%nnest

       allocate(ISCCPHXGstobs(n)%LST_obs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       ISCCPHXGstobs(n)%LST_obs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%lst_obs, &
            "K",1,1)

       LDT_DAobsData(n)%lst_obs%selectStats = 1
    
       if(ISCCPHXGstobs(n)%data_designation.eq."HXG") then
          ISCCPHXGstobs(n)%nc = 3600
          ISCCPHXGstobs(n)%nr = 1800

          gridDesci = 0
          gridDesci(1) = 0     ! 0 -> means lat/lon grid
          gridDesci(2) = 3600   !  x resolution
          gridDesci(3) = 1800   !  y resolution
          gridDesci(4) = -89.95   ! Lower Left Lat
          !gridDesci(5) = 0.05   ! Lower Left Lon
          gridDesci(5) = -179.95   ! Lower Left Lon
          gridDesci(6) = 128
          gridDesci(7) = 89.95   ! Upper Right Lat
          !gridDesci(8) = 359.95   ! Upper Right Lon
          gridDesci(8) = 179.95   ! Upper Right Lon
          gridDesci(9) = .10      ! dx
          gridDesci(10) = .10    ! dy
          gridDesci(20) = 64
          gridDesci(11) = 1 !for the global switch

     
          allocate(ISCCPHXGstobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input (n, gridDesci,&
               ISCCPHXGstobs(n)%n11)


    !!!!   A FEW NOTES
    !!!!   FOUND A BUG IN THE NEIGHBOR_INTERP ROUTINES
    !!!!   SO in the reader, after data is read in I had to shift the ISCCP data to be from
    !!!!   -180 to 180 instead of 0-360 longitude.
    !!!!   SO...the commented sections in gridDesci above are for the native data
    !!!!   and the replacements are for the grid after shifted.
    !!!!   hopefully sometime in the future the bug is fixed or the regrid routines are
    !!!!   replaced

       endif

    enddo
  end subroutine ISCCPHXG_stobsinit
     
end module ISCCP_HXG_obsMod
