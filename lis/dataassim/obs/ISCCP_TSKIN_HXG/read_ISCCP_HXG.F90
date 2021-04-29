!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_ISCCP_HXG
! \label{read_ISCCP_HXG}
!
! !REVISION HISTORY:
!  17 Jun 2010: Sujay Kumar; Updated for use with LPRM AMSRE Version 5. 
!  6 Apr 2021:  John Eylander; Updated for use with NOAA ISCCP HXG data
!
! !INTERFACE: 
subroutine read_ISCCP_HXG(n, k, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use map_utils
  use LIS_pluginIndices
  use ISCCP_HXG_Mod, only : ISCCP_HXG_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the NOAA ISCCP land surface temperature observations
!  from NETCDF files and applies the spatial masking for dense
!  vegetation, rain and clouds. The data is then rescaled
!  to the land surface model's climatology using rescaling
!  algorithms. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  real, parameter        :: minssdev = 1.0
  real, parameter        :: maxssdev = 5.0
  real,  parameter       :: MAX_STO_VALUE=330.0, MIN_STO_VALUE=230.0
  integer                :: status
  integer                :: grid_index
  character*100          :: tskindir
  character*100          :: fname
  logical                :: alarmCheck, file_exists
  integer                :: t,c,r,i,j,p,jj
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: tskinfield, pertField
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  real                   :: obs_unsc(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  real                   :: stobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                   :: st_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                   :: dt
  real                   :: lon
  real                   :: lhour
  real                   :: gmt
  integer                :: zone
  integer                :: fnd
  real, allocatable      :: ssdev(:)
  integer                :: lis_julss
  real                   :: stvalue
  real                   :: model_delta(LIS_rc%obs_ngrid(k))
  real                   :: obs_delta(LIS_rc%obs_ngrid(k))
 
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       tskindir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  data_upd = .false. 
  obs_unsc = LIS_rc%udef

  alarmCheck = LIS_isAlarmRinging(LIS_rc, "ISCCP Tskin read alarm")

  stobs = LIS_rc%udef

  if(alarmCheck.or.ISCCP_HXG_struc(n)%startMode) then
     ISCCP_HXG_struc(n)%startMode = .false.
     
        call create_ISCCP_HXG_filename(tskindir, LIS_rc%yr, LIS_rc%mo, &
             LIS_rc%da, fname)
        
        inquire(file=fname,exist=file_exists)
        
        if(file_exists) then 
           write(LIS_logunit,*) '[INFO] Reading ',trim(fname)
           call read_ISCCP_HXG_data(n,k,fname,stobs)
        else
           write(LIS_logunit,*) '[WARN] Missing ISCCP file: ',trim(fname)
        endif

        ISCCP_HXG_struc(n)%STO_obs  = LIS_rc%udef
        ISCCP_HXG_struc(n)%STO_time = -1


!------------------------------------------------------------------------- 
!   Grab the soil temperature obs and put them into the ISCCP data
!   structure
!-------------------------------------------------------------------------
        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              grid_index = LIS_obs_domain(n,k)%gindex(c,r)
              if(grid_index.ne.-1) then 
                 
                 if(stobs(c+(r-1)*LIS_rc%obs_lnc(k)).ne.-9999.0) then
                    ISCCP_HXG_struc(n)%STO_obs(c,r) = &
                         stobs(c+(r-1)*LIS_rc%obs_lnc(k))
                    lon = LIS_obs_domain(n,k)%lon(c+(r-1)*LIS_rc%obs_lnc(k))
                    
                 endif

              endif
           enddo
        enddo


  endif
  
  
  call ESMF_StateGet(OBS_State,"Observation01",tskinfield,&
       rc=status)
  call LIS_verify(status, 'Error: StateGet Observation01')
  
  call ESMF_FieldGet(tskinfield,localDE=0,farrayPtr=obsl,rc=status)
  call LIS_verify(status, 'Error: FieldGet')
  
  fnd = 0 
  st_current = LIS_rc%udef
 
  ! dt is not defined as absolute value of the time difference to avoid
  ! double counting of the data in assimilation. 

  do r=1,LIS_rc%obs_lnr(k)
     do c=1,LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           grid_index = c+(r-1)*LIS_rc%obs_lnc(k)

           dt = (LIS_rc%gmt - ISCCP_HXG_struc(n)%STO_time(c,r))*3600.0
           if(dt.ge.0.and.dt.lt.LIS_rc%ts) then 
              st_current(c,r) = &
              ISCCP_HXG_struc(n)%STO_obs(c,r)
              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                 obs_unsc(LIS_obs_domain(n,k)%gindex(c,r)) = &
                      st_current(c,r)
              endif
              if(st_current(c,r).ne.LIS_rc%udef) then
                 fnd = 1
              endif
           endif
        endif
     enddo
  enddo
 
!-------------------------------------------------------------------------
!  Transform data to the LSM climatology using a CDF-scaling approach
!-------------------------------------------------------------------------     

  ! Read monthly CDF (only for the current month)
  if (ISCCP_HXG_struc(n)%ntimes.gt.1.and.ISCCP_HXG_struc(n)%cdf_read_opt.eq.1) then
     if (.not. ISCCP_HXG_struc(n)%cdf_read_mon.or.LIS_rc%da.eq.1.and.LIS_rc%hr.eq.0.and.LIS_rc%mn.eq.0.and.LIS_rc%ss.eq.0) then
        call LIS_readMeanSigmaData(n,k,&
             ISCCP_HXG_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             ISCCP_HXG_struc(n)%modelcdffile, &
             "SoilTemp",&
             ISCCP_HXG_struc(n)%model_mu,&
             ISCCP_HXG_struc(n)%model_sigma,&
             LIS_rc%mo)

        call LIS_readMeanSigmaData(n,k,&
             ISCCP_HXG_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             ISCCP_HXG_struc(n)%obscdffile, &
             "SoilTemp",&
             ISCCP_HXG_struc(n)%obs_mu,&
             ISCCP_HXG_struc(n)%obs_sigma,&
             LIS_rc%mo)

        call LIS_readCDFdata(n,k,&
             ISCCP_HXG_struc(n)%nbins,&
             ISCCP_HXG_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             ISCCP_HXG_struc(n)%modelcdffile, &
             "SoilTemp",&
             ISCCP_HXG_struc(n)%model_xrange,&
             ISCCP_HXG_struc(n)%model_cdf,&
             LIS_rc%mo)

        call LIS_readCDFdata(n,k,&
             ISCCP_HXG_struc(n)%nbins,&
             ISCCP_HXG_struc(n)%ntimes,&
             LIS_rc%obs_ngrid(k), &
             ISCCP_HXG_struc(n)%obscdffile, &
             "SoilTemp",&
             ISCCP_HXG_struc(n)%obs_xrange,&
             ISCCP_HXG_struc(n)%obs_cdf,&
             LIS_rc%mo)

         ISCCP_HXG_struc(n)%cdf_read_mon = .true.
     endif
  endif

  if(LIS_rc%dascaloption(k).eq."CDF matching".and.fnd.ne.0) then
     if (ISCCP_HXG_struc(n)%ntimes.gt.1.and.ISCCP_HXG_struc(n)%cdf_read_opt.eq.1) then
        call LIS_rescale_with_CDF_matching(     &
             n,k,                               &
             ISCCP_HXG_struc(n)%nbins,         &
             1,                                 &
             MAX_ST_VALUE,                      &
             MIN_ST_VALUE,                      &
             ISCCP_HXG_struc(n)%model_xrange,  &
             ISCCP_HXG_struc(n)%obs_xrange,    &
             ISCCP_HXG_struc(n)%model_cdf,     &
             ISCCP_HXG_struc(n)%obs_cdf,       &
             st_current)
     else
        call LIS_rescale_with_CDF_matching(     &
             n,k,                               & 
             ISCCP_HXG_struc(n)%nbins,         &
             ISCCP_HXG_struc(n)%ntimes,        &
             MAX_ST_VALUE,                      &
             MIN_ST_VALUE,                      &
             ISCCP_HXG_struc(n)%model_xrange,  &
             ISCCP_HXG_struc(n)%obs_xrange,    &
             ISCCP_HXG_struc(n)%model_cdf,     &
             ISCCP_HXG_struc(n)%obs_cdf,       &
             st_current)
     endif  
  endif
  
  obsl = LIS_rc%udef 
  do r=1, LIS_rc%obs_lnr(k)
     do c=1, LIS_rc%obs_lnc(k)
        if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
           obsl(LIS_obs_domain(n,k)%gindex(c,r))=st_current(c,r)
        endif
     enddo
  enddo
  !-------------------------------------------------------------------------
  !  Apply LSM based QC and screening of observations
  !-------------------------------------------------------------------------  
  call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
       //trim(LIS_ISCCP_HXGobsId)//char(0),n,k,OBS_state)

  call LIS_checkForValidObs(n,k,obsl,fnd,st_current)

  if(fnd.eq.0) then 
     data_upd_flag_local = .false. 
  else
     data_upd_flag_local = .true. 
  endif
        
#if (defined SPMD)
  call MPI_ALLGATHER(data_upd_flag_local,1, &
       MPI_LOGICAL, data_upd_flag(:),&
       1, MPI_LOGICAL, LIS_mpi_comm, status)
#endif
  data_upd = .false.
  do p=1,LIS_npes
     data_upd = data_upd.or.data_upd_flag(p)
  enddo
  
  if(data_upd) then 

     do t=1,LIS_rc%obs_ngrid(k)
        gid(t) = t
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo
  
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true. , rc=status)
     call LIS_verify(status)

     if(LIS_rc%obs_ngrid(k).gt.0) then 
        call ESMF_AttributeSet(tskinField,"Grid Number",&
             gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)
        
        call ESMF_AttributeSet(tskinField,"Assimilation Flag",&
             assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
        call LIS_verify(status)

        call ESMF_AttributeSet(tskinfield, "Unscaled Obs",&
             obs_unsc, itemCount=LIS_rc%obs_ngrid(k), rc=status)
        call LIS_verify(status, 'Error in setting Unscaled Obs attribute')

     endif

     if(LIS_rc%dascaloption(k).eq."CDF matching") then 
        if(ISCCP_HXG_struc(n)%useSsdevScal.eq.1) then
           call ESMF_StateGet(OBS_Pert_State,"Observation01",pertfield,&
                rc=status)
           call LIS_verify(status, 'Error: StateGet Observation01')
           
           allocate(ssdev(LIS_rc%obs_ngrid(k)))
           ssdev = ISCCP_HXG_struc(n)%ssdev_inp
           
           if(ISCCP_HXG_struc(n)%ntimes.eq.1) then
              jj = 1
           else
              jj = LIS_rc%mo
           endif
           do t=1,LIS_rc%obs_ngrid(k)
              if(ISCCP_HXG_struc(n)%obs_sigma(t,jj).gt.0) then
                 ssdev(t) = ssdev(t)*ISCCP_HXG_struc(n)%model_sigma(t,jj)/&
                    ISCCP_HXG_struc(n)%obs_sigma(t,jj)
                 if(ssdev(t).lt.minssdev) then 
                    ssdev(t) = minssdev
                 endif
              endif
           enddo
           
           if(LIS_rc%obs_ngrid(k).gt.0) then 
              call ESMF_AttributeSet(pertField,"Standard Deviation",&
                   ssdev,itemCount=LIS_rc%obs_ngrid(k),rc=status)
              call LIS_verify(status)
           endif
           deallocate(ssdev)
        endif
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif

end subroutine read_ISCCP_HXG

!BOP
! 
! !ROUTINE: read_ISCCP_HXG_data
! \label{read_ISCCP_HXG_data}
!
! !INTERFACE:
subroutine read_ISCCP_HXG_data(n, k, fname, tskin_obs_ip)
! 
! !USES:   

  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod
  use LIS_timeMgrMod
  use ISCCP_HXG_Mod, only : ISCCP_HXG_struc
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: n 
  integer                       :: k
  character (len=*)             :: fname
  real                          :: tskin_obs_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))


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
!  \item[stobs_obs_ip]    ISCCP skin temperature data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_NETCDF4)

  integer             :: i, j, istat, nid, ios
  integer             :: tskinid, cloudsid, tmptabid, latid, lonid

  integer(hsize_t), allocatable  :: dims(:)
  integer(hsize_t), dimension(2) :: dimst
  integer(hsize_t), dimension(2) :: count_file
  integer(hsize_t), dimension(2) :: count_mem
  integer(hid_t)                 :: memspace
  integer(hid_t)                 :: dataspace
  integer                        :: memrank = 2 ! scaler--> rank = 0 ; 2D array--> rank = 2
  integer(hsize_t), dimension(2) :: offset_mem = (/0,0/)
  integer(hsize_t), dimension(2) :: offset_file = (/0,0/)
  integer(hid_t)                 :: file_id


  integer                        :: c,r,t
  integer                        :: bit_temp
  integer                        :: status

  real                           :: assigned_temp

  dimst      = (/ISCCP_HXG_struc(n)%nc, ISCCP_HXG_struc(n)%nr/)
  count_file = (/ISCCP_HXG_struc(n)%nc, ISCCP_HXG_struc(n)%nr/)
  count_mem  = (/ISCCP_HXG_struc(n)%nc, ISCCP_HXG_struc(n)%nr/)
  
  allocate(dims(2))

  
  integer, allocatable            :: ts_temp(:,:)
  real, allocatable               :: tmp_lat(:,:)
  real, allocatable               :: tmp_lon(:,:)
  real, allocatable               :: clouds(:,:)
  real, allocatable               :: temp_array(:,:)
  real, allocatable               :: tmp_tab

  logical*1    :: tskin_array_logical(ISCCP_HXG_struc(n)%nc*ISCCP_HXG_struc(n)%nr)
  real         :: tskin_array(ISCCP_HXG_struc(n)%nc*ISCCP_HXG_struc(n)%nr)

  allocate(tskin_array(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  
  allocate(tmp_lat(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
  allocate(tmp_lon(LIS_rc%lnc(n)*LIS_rc%lnr(n)))


  dims(1) = ISCCP_HXG_struc(n)%nc
  dims(2) = ISCCP_HXG_struc(n)%nr

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

!-------------------------------------------------------------------------
!   Open the NETCDF file
!-------------------------------------------------------------------------

     ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
     call LIS_VERIFY(ios, 'Error opening file'//trim(fname)

!-------------------------------------------------------------------------
!   Obtain the NETCDF variable record ID, then get the variable for TSKIN
!   and cloud top temperature.  We will use the cloud mask as a filter to
!   make sure we are only using the TSKIN data
!-------------------------------------------------------------------------

     ios = nf90_inq_varid(nid, 'IR retrieved cloud top temperature or surface skin temperature', tskinid)
     call LIS_VERIFY(ios, 'Error nf90_inq_varid: ISCCP skin temperature')

     ios = nf90_inq_dimid(nid,"lon",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in ISCCP Skin Temp')

     ios = nf90_inq_dimid(nid,"lat",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in ISCCP Skin Temp')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in ISCCP Skin Temp')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in ISCCP Skin Temp')

     ISCCP_HXG_struc(n)%nc = nc
     ISCCP_HXG_struc(n)%nr = nr

     allocate(ts_temp(ISCCP_HXG_struc(n)%nc, ISCCP_HXG_struc(n)%nr))
     allocate(clouds(ISCCP_HXG_struc(n)%nc, ISCCP_HXG_struc(n)%nr))
     allocatable(temp_array(ISCCP_HXG_struc(n)%nc, ISCCP_HXG_struc(n)%nr))

     ios = nf90_get_var(nid, tskinid, ts_temp)
     call LIS_VERIFY(ios, 'Error nf90_get_var: ISCCP Surface Skin Temperature')

!-------------------------------------------------------------------------
!   read in the cloud mask as a filter
!-------------------------------------------------------------------------

     ios = nf90_inq_varid(nid, 'Cloud mask', cloudsid)
     call LIS_VERIFY(ios, 'Error nf90_inq_varid: ISCCP Cloud Mask')

     ios = nf90_get_var(nid, cloudsid, clouds)
     call LIS_VERIFY(ios, 'Error nf90_get_var: ISCCP Cloud mask')

!-------------------------------------------------------------------------
!   read in the conversion table to convert the byte data from the ISCCP
! temperature record into a floating point actual temperature value
!-------------------------------------------------------------------------

     ios = nf90_inq_dimid(nid,"count",ncount)
     call LIS_verify(ios,'Error in nf90_inq_dimid obtaining temp counts in ISCCP Skin Temp')

     ios = nf90_inquire_dimension(nid,ncount, len=temp_tab_size)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in ISCCP Skin Temp')

     allocate(tmptab(temp_tab_size))
     
     ios = nf90_inq_varid(nid, 'Count to temperature conversion table', tmptabid)
     call LIS_VERIFY(ios, 'Error nf90_inq_varid: ISCCP temp conversion table')

     ios = nf90_get_var(nid, tmptabid, tmptab)
     call LIS_VERIFY(ios, 'Error nf90_get_var: ISCCP Surface Skin Temperature')


#if 0 
! =============================================================================
! JBE the following section added
!   1 - to convert the byte values to temperature values
!   2 - to filter the cloudy grid points

     
  temp_array(:,:) = LIS_rc%udef
  do r=1,ISCCP_HXG_struc(n)%nr
     do c=1,ISCCP_HXG_struc(n)%nc
        if(clouds(c,r) .eq. 1 ) then !JBE filter out value that are cloudy
           bit_temp=ts_temp(c,r)
           assigned_temp=tmptab(bit_temp)
           temp_array(c,r) = assigned_temp
        endif
     enddo
  enddo

  tskin_array_logical = .false.
  t = 1

  do r=1,ISCCP_HXG_struc(n)%nr
     do c=1,ISCCP_HXG_struc(n)%nc
        tskin_data(t) = temp_array(c,r)
        if(tskin_data(t).ne.-9999.0) then
  
           tskin_array_logical(t) = .true.
        endif
        t = t+1
     enddo
  enddo


! JBE end added section
! =============================================================================
#endif

!--------------------------------------------------------------------------
! Interpolate to the LIS running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LIS_rc%obs_gridDesc(k,:),&
       tskin_array_logical, tskin_data, tskin_array_logical_ip,  &
       tskin_obs_ip, ISCCP_HXG_struc(n)%nc*ISCCP_HXG_struc(n)%nr, &
       LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
       ISCCP_HXG_struc(n)%rlat, ISCCP_HXG_struc(n)%rlon,&
       ISCCP_HXG_struc(n)%n11, LIS_rc%udef, ios)

  deallocate(clouds)
  deallocate(tmptab)
  deallocate(dims)
  deallocate(temp_array)
  deallocate(ts_temp)
  deallocate(tmp_lat)
  deallocate(tmp_lon)

#endif

end subroutine read_ISCCP_HXG_data


!BOP
! !ROUTINE: create_ISCCP_filename
! \label{create_ISCCP_filename}
! 
! !INTERFACE: 
subroutine create_ISCCP_filename(ndir, yr, mo,da, hr, filename)
! !USES:   

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
!  \item[filename] Generated SMAP filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr

  
  filename = trim(ndir)//'/Y'//trim(fyr)//'/M'//trim(fmo)//&
          '/ISCCP.HXG.v01r00.GLOBAL.'//trim(fyr)//'.'//trim(fmo)//&
          '.'//trim(fda)//'.'//trim(fhr)//'00.GPC.10KM.EQ0.10.nc'

  endif

end subroutine create_ISCCP_filename





