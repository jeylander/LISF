!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: readNASASMAPsmObs
! \label{readNASASMAPsmObs}
!
! !REVISION HISTORY:
!  21 July 2010: Sujay Kumar, Initial Specification
!  12 Feb 2018: Mahdi Navari, openwater proximity detection was added
!                         edited to read New version of the SPL3SMP_R14 (file structure
!                          differs from the previous versions)
!  31 Aug 2018: Mahdi Navari, Edited to read SPL3SMP.005 & SPL3SMP_E.002
!  04 Jun 2019: Sujay Kumar, Updated to support SMAP L2 retrievals
!  15 Aug 2019: Mahdi Navari, There are several version of SMAP sm data available in each directory
!                  with different Release number and different CRID Version Number. The reader was
!                  modified to read the latest version of data (the reader no longer reads the symbolic
!                  link to the SMAP sm data)
!
! !INTERFACE:
subroutine readISCCPHXGSTObs(n)
! !USES:
   use ESMF
   use LDT_coreMod
   use LDT_logMod
   use LDT_timeMgrMod
   use LDT_DAobsDataMod
   use ISCCP_HXG_obsMod
   use map_utils

   implicit none
! !ARGUMENTS:
   integer, intent(in) :: n
!
! !DESCRIPTION:
!
! This subroutine provides the data reader for the ESACCI
! soil moisture retrieval product.
!
!EOP

   real*8            :: timenow
   logical           :: alarmCheck
   logical           :: file_exists
   integer           :: c, r, i, j
   character*100     :: fname
   integer           :: mn_ind
   integer           :: yr, mo, da, hr, mn, ss
   integer           :: doy
   integer           :: ftn
   integer           :: ierr
   real              :: gmt
   character*8       :: yyyymmdd
   character*4       :: yyyy
   character*2       :: mm, dd, hh
   character*100     :: ISCCP_HXG_filename(10)
   
   real              :: LST_obs(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   
   character(len=3) :: CRID

!-----------------------------------------------------------------------
! Read in the ISCCP observations.  Main routine below.
!-----------------------------------------------------------------------
   ISCCPHXGstobs(n)%LST_obs = LDT_rc%udef
   LST_obs = LDT_rc%udef

   if (ISCCPHXGstobs(n)%data_designation .eq. "HXG") then

      if (LDT_rc%ts .ne. 10800) then
         write (LDT_logunit, *) '[ERR] Please set the LDT timestep to 3hr'
         write (LDT_logunit, *) '[ERR] This is required for ISCCP data processing'
         call LDT_endrun()
      endif

      write (yyyymmdd, '(i4.4,2i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
      write (yyyy, '(i4.4)') LDT_rc%yr
      write (mm, '(i2.2)') LDT_rc%mo
      write (dd, '(i2.2)') LDT_rc%da
      write (hh, '(i2.2)') LDT_rc%hr

     call create_ISCCP_HXG_filename(ISCCPHXGstobs(n)%odir, LDT_rc%yr, LDT_rc%mo, &
        LDT_rc%da, LDT_rc%hr, fname)

   
      write(LDT_logunit,*) '[ALERT] ISCCP file info: ',fname, ' ',ISCCPHXGstobs(n)%odir, &
        'LDT Year ', LDT_rc%yr, 'LDT Month ', LDT_rc%mo, &
        'LDT Day ', LDT_rc%da, &
        'LDT Hour ', LDT_rc%hr, &
        'ISCCP Filename: ', fname

      inquire(file=fname,exist=file_exists)

      if(file_exists) then
        write(LDT_logunit,*) '[INFO] Reading ',trim(fname)
        call read_ISCCP_HXG_data(n,fname,LST_obs)
      else
        write(LDT_logunit,*) '[WARN] Missing ISCCP file: ',trim(fname)
      endif

      ISCCPHXGstobs(n)%LST_obs  = LDT_rc%udef

!-------------------------------------------------------------------------
!   Grab the soil temperature obs and put them into the ISCCP data
!   structure
!-------------------------------------------------------------------------
      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            if(LST_obs(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then
                ISCCPHXGstobs(n)%LST_obs(c,r) = &
                        LST_obs(c+(r-1)*LDT_rc%lnc(n))
            endif
         enddo
      enddo

      call LDT_logSingleDAobs(n, LDT_DAobsData(n)%lst_obs, &
            ISCCPHXGstobs(n)%LST_obs, vlevel=1)

  endif


end subroutine readISCCPHXGSTObs


!BOP
!
! !ROUTINE: read_ISCCP_HXG_data
! \label{read_ISCCP_HXG_data}
!
! !INTERFACE:
subroutine read_ISCCP_HXG_data(n, fname, tskin_obs_ip)
!
! !USES:
  use LDT_coreMod
  use LDT_logMod
  use map_utils
  use ISCCP_HXG_obsMod, only : ISCCPHXGstobs
  use LDT_timeMgrMod
  
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS:
!
  integer                       :: n
  character (len=*)             :: fname
  real, intent(inout)           :: tskin_obs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1                     :: tskin_array_logical_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))


! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION:
!  This subroutine reads the NOAA/NESDIS ISCCP NETCDF file
!       and applies the data  quality flags to filter the
!       data. !\normalsize
!
!
!  The arguments are:
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the ISCCP file
!  \item[stobs_obs_ip]    ISCCP skin temperature data processed to the LDT domain
! \end{description}
!
!
!EOP

#if (defined USE_NETCDF4)

  integer             :: i, j, istat, nid, ios
  integer             :: ncid, nrid, nc, nr, ncount
  integer             :: tskinid, cloudsid, tmptabid, latid, lonid

  integer,          allocatable  :: dims(:)
  integer,          dimension(2) :: dimst
  integer,          dimension(2) :: count_file
  integer,          dimension(2) :: count_mem
  integer                        :: memspace
  integer                        :: dataspace
  integer                        :: memrank = 2 ! scaler--> rank = 0 ; 2D array--> rank = 2
  integer,          dimension(2) :: offset_mem = (/0,0/)
  integer,          dimension(2) :: offset_file = (/0,0/)
  integer                        :: file_id
  integer                        :: rec_length

  integer                        :: c,r,t
  integer                        :: bit_temp
  integer                        :: status
  integer                        :: temp_tab_size

  real                           :: assigned_temp

  integer*2, allocatable            :: ts_temp(:,:)
  real*8, allocatable               :: tmp_lat(:)
  real*8, allocatable               :: tmp_lon(:)
  integer*2, allocatable            :: clouds(:,:)
  integer, allocatable              :: temp_array(:,:)
  integer, allocatable              :: new_temp_array(:,:)
  integer, allocatable              :: tmp_tab(:)

  logical*1    :: tskin_array_logical(ISCCPHXGstobs(n)%nc*ISCCPHXGstobs(n)%nr)
  real         :: tskin_data(ISCCPHXGstobs(n)%nc*ISCCPHXGstobs(n)%nr)

  dimst      = (/ISCCPHXGstobs(n)%nc, ISCCPHXGstobs(n)%nr/)
  count_file = (/ISCCPHXGstobs(n)%nc, ISCCPHXGstobs(n)%nr/)
  count_mem  = (/ISCCPHXGstobs(n)%nc, ISCCPHXGstobs(n)%nr/)
  
  allocate(dims(2))
  
  allocate(tmp_lat(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
  allocate(tmp_lon(LDT_rc%lnc(n)*LDT_rc%lnr(n)))


  dims(1) = ISCCPHXGstobs(n)%nc
  dims(2) = ISCCPHXGstobs(n)%nr

!-------------------------------------------------------------------------
!   Open the NETCDF file
!-------------------------------------------------------------------------

  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
  call LDT_VERIFY(ios, 'Error opening file'//trim(fname))

!-------------------------------------------------------------------------
!   Obtain the NETCDF variable record ID, then get the variable for TSKIN
!   and cloud top temperature.  We will use the cloud mask as a filter to
!   make sure we are only using the TSKIN data
!-------------------------------------------------------------------------

  ios=nf90_inq_varid(nid, "itmp", tskinid)
  call LDT_VERIFY(ios,'[ERR] Error nf90_inq_varid: ISCCP skin temperature')

  ios=nf90_inq_dimid(nid,"lon",ncId)
  call LDT_verify(ios,'[ERR] Error in nf90_inq_dimid in ISCCP Skin Temp')

  ios=nf90_inquire_dimension(nid,ncId, len=nc)
  call LDT_verify(ios, '[ERR] Error in nf90_inquire_dimension in ISCCP Skin Temp')

  ios=nf90_inq_dimid(nid,"lat",nrId)
  call LDT_verify(ios, '[ERR] Error in nf90_inq_dimid in ISCCP Skin Temp')

  ios=nf90_inquire_dimension(nid,nrId, len=nr)
  call LDT_verify(ios, '[ERR] Error in nf90_inquire_dimension in ISCCP Skin Temp')

  ISCCPHXGstobs(n)%nc = nc
  ISCCPHXGstobs(n)%nr = nr

  allocate(ts_temp(ISCCPHXGstobs(n)%nc, ISCCPHXGstobs(n)%nr))
  allocate(clouds(ISCCPHXGstobs(n)%nc, ISCCPHXGstobs(n)%nr))
  allocate(temp_array(ISCCPHXGstobs(n)%nc, ISCCPHXGstobs(n)%nr))
  allocate(new_temp_array(ISCCPHXGstobs(n)%nc, ISCCPHXGstobs(n)%nr))

  ios=nf90_get_var(nid, tskinid, ts_temp)
  call LDT_VERIFY(ios, '[ERR] Error nf90_get_var: ISCCP Surface Skin Temperature')

!-------------------------------------------------------------------------
!   read in the cloud mask as a filter
!-------------------------------------------------------------------------

  ios=nf90_inq_varid(nid, "cloud", cloudsid)
  call LDT_VERIFY(ios, '[ERR] Error nf90_inq_varid: ISCCP Cloud Mask')

  ios=nf90_get_var(nid, cloudsid, clouds)
  call LDT_VERIFY(ios, '[ERR] Error nf90_get_var: ISCCP Cloud mask')

!-------------------------------------------------------------------------
!   read in the conversion table to convert the byte data from the ISCCP
! temperature record into a floating point actual temperature value
!-------------------------------------------------------------------------

  ios=nf90_inq_dimid(nid, "count", ncount)
  call LDT_verify(ios, &
    '[ERR] Error in nf90_inq_dimid obtaining temp counts in ISCCP Skin Temp')

  ios=nf90_inquire_dimension(nid, ncount, len=temp_tab_size)
  call LDT_verify(ios,'[ERR] Error in nf90_inquire_dimension in ISCCP Skin Temp')

  allocate(tmp_tab(temp_tab_size))
  
  ios=nf90_inq_varid(nid,"tmptab", tmptabid)
  call LDT_VERIFY(ios, '[ERR] Error nf90_inq_varid: ISCCP temp conversion table')

  ios=nf90_get_var(nid, tmptabid, tmp_tab)
  call LDT_VERIFY(ios, '[ERR] Error nf90_get_var: ISCCP Surface Skin Temperature')

  write(LDT_logunit,*) '[ALERT] INSIDE ISCCP reader: ', &
     'COLUMNS/ROWS ', nc, nr

  write(LDT_logunit,*) '[ALERT] INSIDE ISCCP reader: ', &
   'clouds values 2513/1209 ', clouds(2513,1209)

  write(LDT_logunit,*) '[ALERT] INSIDE ISCCP reader: ', &
    'ts_temp values 2513/1209 ', ts_temp(2513,1209)

! =============================================================================
! JBE the following section added
!   1 - to convert the byte values to temperature values
!   2 - to filter the cloudy grid points

  
  write(LDT_logunit,*) '[ALERT] INSIDE ISCCP reader: ', &
    'temp_table= ', tmp_tab, LDT_rc%lnc(n), LDT_rc%lnr(n), ISCCPHXGstobs(n)%nc, &
         ISCCPHXGstobs(n)%nr


  temp_array(:,:) = LDT_rc%udef
  do r=1,ISCCPHXGstobs(n)%nr
     do c=1,ISCCPHXGstobs(n)%nc
        if(clouds(c,r) .eq. 0 ) then !JBE filter out value that are cloudy
           bit_temp=ts_temp(c,r)
           assigned_temp=tmp_tab(bit_temp)
           temp_array(c,r) = assigned_temp
        endif
     enddo
  enddo

  !!!  THIS NEXT SECTION SHOULD EVENTUALLY BE REMOVED
  !!!  There is a feature in the Neighbor_Interp routines
  !!!  that can't currently handle 0-360 longitude domains data
  !!!  so I use these next two lines of code to swap the domain into a
  !!!  -180 to 180 domain.  Hopefully can be removed at some point.  If so
  !!!  then the MOD will need to be updated -- griddesci section.

  new_temp_array(1801:3600,:) = temp_array(1:1800,:)
  new_temp_array(1:1800,:) = temp_array(1801:3600,:)

  tskin_array_logical = .false.
  t = 1

  do r=1,ISCCPHXGstobs(n)%nr
     do c=1,ISCCPHXGstobs(n)%nc
        tskin_data(t) = new_temp_array(c,r)
        if(tskin_data(t).ne.-9999.0) then
  
           tskin_array_logical(t) = .true.
        endif
        t = t+1
     enddo
  enddo


! JBE end added section
! =============================================================================


!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!--------------------------------------------------------------------------
  call neighbor_interp(LDT_rc%gridDesc,&
       tskin_array_logical, tskin_data, tskin_array_logical_ip,  &
       tskin_obs_ip, ISCCPHXGstobs(n)%nc*ISCCPHXGstobs(n)%nr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       ISCCPHXGstobs(n)%n11, LDT_rc%udef, ios)

  deallocate(clouds)
  deallocate(tmp_tab)
  deallocate(dims)
  deallocate(temp_array)
  deallocate(ts_temp)
  deallocate(tmp_lat)
  deallocate(tmp_lon)
  deallocate(new_temp_array)

#endif

end subroutine read_ISCCP_HXG_data



!BOP
! !ROUTINE: create_ISCCP_HXG_filename
! \label{create_ISCCP_HXG_filename}
!
! !INTERFACE:
subroutine create_ISCCP_HXG_filename(ndir, yr, mo,da, hr, filename)
! !USES:
  use LDT_logMod
  implicit none
! !ARGUMENTS:
  character(len=*)  :: filename

  integer           :: yr, mo, da, hr
  character (len=*) :: ndir
!
! !DESCRIPTION:
!  This subroutine creates the ISCCP filename (from NESDIS)
!  based on the time and date
!
!  The arguments are:
!  \begin{description}
!  \item[ndir] name of the ISCCP Skin Temperature data directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated ISCCP filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr

  write(LDT_logunit,*) '[ALERT] INSIDE ISCCP file generator: ', &
        'Input Directory ',ndir, &
        'LDT Year ', yr, 'LDT Month ', mo, &
        'LDT Day ', da, &
        'LDT HR ', hr

  
  filename = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//&
          '/ISCCP.HXG.v01r00.GLOBAL.'//trim(fyr)//'.'//trim(fmo)//&
          '.'//trim(fda)//'.'//trim(fhr)//'00.GPC.10KM.EQ0.10.nc'

   write(LDT_logunit,*) '[ALERT] INSIDE ISCCP file generator: ', &
        'filename= ', filename

end subroutine create_ISCCP_HXG_filename

