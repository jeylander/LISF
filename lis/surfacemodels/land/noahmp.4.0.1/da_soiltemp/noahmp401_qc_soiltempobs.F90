!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_qc_soiltempobs
! \label{NoahMP401_qc_soiltempobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 15 Dec 2018: Mahdi Navari; Modified for NoahMP401
! 21 Mar 2021: John Eylander (ERDC); Modified for Soil Temperature
!
! !INTERFACE:
subroutine NoahMP401_qc_soiltempobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use LIS_DAobservationsMod
  use NoahMP401_lsmMod


  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the soil temperature observations
!  are flagged when LSM indicates that (1) rain is falling
!  (2) ground is fully or partially covere with snow and
!  (3) ground is covered with vegatation (more than 50%) and/or
!  (4) vegetation water content is higher than xx.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_tslb_field

  real, pointer            :: tslb_obs(:)
  integer                  :: t
  integer                  :: gid
  integer                  :: status
  real                     :: lat,lon

! mn
  real                     :: tsk(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: tslb1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: tslb2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: tslb3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: tslb4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  
  real                     :: stc1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: stc4(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                     :: vegt(LIS_rc%npatch(n,LIS_rc%lsm_index))
  
  real                     :: rainf_obs(LIS_rc%obs_ngrid(k))
  real                     :: sneqv_obs(LIS_rc%obs_ngrid(k))
  real                     :: sca_obs(LIS_rc%obs_ngrid(k))
  real                     :: shdfac_obs(LIS_rc%obs_ngrid(k))
  real                     :: t1_obs(LIS_rc%obs_ngrid(k))
  real                     :: tsk_obs(LIS_rc%obs_ngrid(k))
  real                     :: tslb1_obs(LIS_rc%obs_ngrid(k))
  real                     :: tslb2_obs(LIS_rc%obs_ngrid(k))
  real                     :: tslb3_obs(LIS_rc%obs_ngrid(k))
  real                     :: tslb4_obs(LIS_rc%obs_ngrid(k))
  
  real                     :: stc1_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc2_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc3_obs(LIS_rc%obs_ngrid(k))
  real                     :: stc4_obs(LIS_rc%obs_ngrid(k))
  real                     :: vegt_obs(LIS_rc%obs_ngrid(k))


  call ESMF_StateGet(OBS_State,"Observation01",obs_tslb_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in NoahMP401_qc_soiltempobs")
  call ESMF_FieldGet(obs_tslb_field,localDE=0,farrayPtr=tslb_obs,rc=status)
  call LIS_verify(status,& 
       "ESMF_FieldGet failed in NoahMP401_qc_soiltempobs")
  
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     tsk(t) = noahmp401_struc(n)%noahmp401(t)%tsk
     tslb1(t) = noahmp401_struc(n)%noahmp401(t)%tslb(1)
     tslb2(t) = noahmp401_struc(n)%noahmp401(t)%tslb(2)
     tslb3(t) = noahmp401_struc(n)%noahmp401(t)%tslb(3)
     tslb4(t) = noahmp401_struc(n)%noahmp401(t)%tslb(4)

     vegt(t) = noahmp401_struc(n)%noahmp401(t)%vegetype
     
  enddo

  call LIS_convertPatchSpaceToObsSpace(n,k,&       
       LIS_rc%lsm_index, &
       noahmp401_struc(n)%noahmp401(:)%prcp,&
       rainf_obs)  ! MN prcp is total precip 
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noahmp401_struc(n)%noahmp401(:)%sneqv,&
       sneqv_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noahmp401_struc(n)%noahmp401(:)%snowc,&   ! MP36 fsno
       sca_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noahmp401_struc(n)%noahmp401(:)%fveg,&
       shdfac_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       noahmp401_struc(n)%noahmp401(:)%tg,&
       t1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       tslb1,&
       tslb1_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       tslb2,&
       tslb2_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       tslb3,&
       tslb3_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       tslb4,&
       tslb4_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       tsk,&
       tsk_obs)
  call LIS_convertPatchSpaceToObsSpace(n,k,&
       LIS_rc%lsm_index, &
       vegt,&
       vegt_obs)

  do t = 1,LIS_rc%obs_ngrid(k)

     if(tslb_obs(t).ne.LIS_rc%udef) then
! MN: check for rain
        if(rainf_obs(t).gt.3E-6) then   ! Var name Noah36 --> rainf 
           tslb_obs(t) = LIS_rc%udef
!           print*, 'rainf ',gid,t,noahmp401_struc(n)%noahmp401(t)%prcp
! MN: check for frozen soil
        elseif(tsk_obs(t).le.LIS_CONST_TKFRZ) then
           tsk_obs(t) = LIS_rc%udef
        elseif(tslb1_obs(t).le.LIS_CONST_TKFRZ) then
           tslb_obs(t) = LIS_rc%udef
        elseif(tslb2_obs(t).le.LIS_CONST_TKFRZ) then
           tslb_obs(t) = LIS_rc%udef
        elseif(tslb3_obs(t).le.LIS_CONST_TKFRZ) then
           tslb_obs(t) = LIS_rc%udef
        elseif(tslb4_obs(t).le.LIS_CONST_TKFRZ) then
           tslb_obs(t) = LIS_rc%udef
        elseif(t1_obs(t).le.LIS_CONST_TKFRZ) then ! Var name Noah36 --> t1
           tslb_obs(t) = LIS_rc%udef
        elseif(vegt_obs(t).le.4) then !forest types ! Var name Noah36 --> vegt
           tslb_obs(t) = LIS_rc%udef
 ! MN: check for snow  
       elseif(sneqv_obs(t).gt.0.001) then 
           tslb_obs(t) = LIS_rc%udef
        elseif(sca_obs(t).gt.0.0001) then  ! Var name sca 
           tslb_obs(t) = LIS_rc%udef
 ! MN: check for green vegetation fraction NOTE: threshold incerased from 0.5 to 0.7 
       elseif(shdfac_obs(t).gt.0.7) then  ! var name Noah36 shdfac 12-month green veg. frac.  
           tslb_obs(t) = LIS_rc%udef

        endif
     endif
  enddo

end subroutine NoahMP401_qc_soiltempobs

