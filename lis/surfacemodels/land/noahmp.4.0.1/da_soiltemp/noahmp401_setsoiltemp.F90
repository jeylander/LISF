!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_setsoiltemp
!  \label{NoahMP401_setsoiltemp}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari: Modified for NoahMP401
! 21 Mar 2021: John Eylander (ERDC); Modified for soil temperature assimilation
! 
! Apply the update if it met the update conditions
! Update conditions: 
!                  1- Prior TSLB + increment > MIN_THRESHOLD
!                  2- Prior TSLB + increment < MAX_THRESHOLD
! There are 3 cases 
! 1- If all the ensemble members met the update conditions --> apply the update
! 2- If more than 50% of the ensemble members met the update condition --> 
!    apply the update for that members and set the other member to the mean 
!    value of the ensemble (i.e. mean of the members that met the conditions)
! 3- If less then 50% of the ensemble members met the update conditions --> 
!    adjust the states    


! !INTERFACE:
subroutine NoahMP401_setsoiltemp(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use NoahMP401_lsmMod



  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil temperature prognostic variables to noah's
!  model space. 
! 
!EOP
  real, parameter        :: MIN_THRESHOLD = 268.0
  real, parameter        :: MAX_THRESHOLD = 323.0
  type(ESMF_Field)       :: tslb1Field
  type(ESMF_Field)       :: tslb2Field
  type(ESMF_Field)       :: tslb3Field
  type(ESMF_Field)       :: tslb4Field
  type(ESMF_Field)       :: tskinField
  real, pointer          :: tskin(:)
  real, pointer          :: soilt1(:)
  real, pointer          :: soilt2(:)
  real, pointer          :: soilt3(:)
  real, pointer          :: soilt4(:)
  integer                :: t, j,i, gid, m, t_unpert
  integer                :: status
  real                   :: delta(4)
  real                   :: delta1,delta2,delta3,delta4
  real                   :: tmpval
  logical                :: bounds_violation
  integer                :: nIter
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: ens_flag(LIS_rc%nensem(n))
! mn
  real                   :: tmp(LIS_rc%nensem(n)), tmp0(LIS_rc%nensem(n))
  real                   :: tmp1(LIS_rc%nensem(n)),tmp2(LIS_rc%nensem(n))
  real                   :: tmp3(LIS_rc%nensem(n)),tmp4(LIS_rc%nensem(n))
  real                   :: tmp5(LIS_rc%nensem(n))
  logical                :: update_flag_tile(LIS_rc%npatch(n,LIS_rc%lsm_index))
  logical                :: flag_ens(LIS_rc%ngrid(n))
  logical                :: flag_tmp(LIS_rc%nensem(n))
  logical                :: update_flag_ens(LIS_rc%ngrid(n))
  logical                :: update_flag_new(LIS_rc%ngrid(n))
  integer                :: RESULT, pcount, icount
  real                   :: MaxEnsST1, MaxEnsST2, MaxEnsST3, MaxEnsST4, MaxEnsST5
  real                   :: MinEnsST1, MinEnsST2, MinEnsST3, MinEnsST4, MinEnsST5
  real                   :: stc_rnd, stc_tmp
  INTEGER, DIMENSION (1) :: seed 

  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 1",tslb1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Temperature Layer 1 failed in NoahMP401_setsoiltemp")
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 2",tslb2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Temperature Layer 2 failed in NoahMP401_setsoiltemp")
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 3",tslb3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Temperature Layer 3 failed in NoahMP401_setsoiltemp")
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 4",tslb4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Temperature Layer 4 failed in NoahMP401_setsoiltemp")
  call ESMF_StateGet(LSM_State,"Skin Temperature",tskinField,rc=status)
  call LIS_verify(status,&
     "ESMF_StateSet: Skin Temperature failed in NoahMP401_setsoiltemp")

  call ESMF_FieldGet(tslb1Field,localDE=0,farrayPtr=soilt1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 1 failed in NoahMP401_setsoiltemp")
  call ESMF_FieldGet(tslb2Field,localDE=0,farrayPtr=soilt2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 2 failed in NoahMP401_setsoiltemp")
  call ESMF_FieldGet(tslb3Field,localDE=0,farrayPtr=soilt3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 3 failed in NoahMP401_setsoiltemp")
  call ESMF_FieldGet(tslb4Field,localDE=0,farrayPtr=soilt4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 4 failed in NoahMP401_setsoiltemp")
  call ESMF_FieldGet(tskinField,localDE=0,farrayPtr=tskin,rc=status)
  call LIS_verify(status,&
     "ESMF_FieldGet: Skin Temperature failed in NoahMP401_setsoiltemp")

  update_flag = .true. 
  update_flag_tile= .true. 

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row) 
     
     !MN: delta = X(+) - X(-)
     !NOTE: "noahmp401_updatesoilm.F90" updates the soilm_(t)
     deltats = tskin(t)-noahmp401_struc(n)%noahmp401(t)%tsk
     delta1 = soilt1(t)-noahmp401_struc(n)%noahmp401(t)%tslb(1)
     delta2 = soilt2(t)-noahmp401_struc(n)%noahmp401(t)%tslb(2)
     delta3 = soilt3(t)-noahmp401_struc(n)%noahmp401(t)%tslb(3)
     delta4 = soilt4(t)-noahmp401_struc(n)%noahmp401(t)%tslb(4)

     ! MN: check    MIN_THRESHOLD < Skin Temperature < threshold
     if(noahmp401_struc(n)%noahmp401(t)%tsk+deltats.gt.MIN_THRESHOLD .and.&
          noahmp401_struc(n)%noahmp401(t)%tsk+deltats.lt.&
            MAX_THRESHOLD) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        ! MN save the flag for each tile (col*row*ens)   (64*44)*20
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     ! MN: check    MIN_THRESHOLD < Soil Temperature Layer 1 < threshold
     if(noahmp401_struc(n)%noahmp401(t)%tslb(1)+delta1.gt.MIN_THRESHOLD .and.&
          noahmp401_struc(n)%noahmp401(t)%tslb(1)+delta1.lt.&
            MAX_THRESHOLD) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        ! MN save the flag for each tile (col*row*ens)   (64*44)*20
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     !  Layer 2
     if(noahmp401_struc(n)%noahmp401(t)%tslb(2)+delta2.gt.MIN_THRESHOLD .and.&
          noahmp401_struc(n)%noahmp401(t)%tslb(2)+delta2.lt.&
            MAX_THRESHOLD) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     !Layer 3
     if(noahmp401_struc(n)%noahmp401(t)%tslb(3)+delta3.gt.MIN_THRESHOLD .and.&
          noahmp401_struc(n)%noahmp401(t)%tslb(3)+delta3.lt.&
            MAX_THRESHOLD) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
     ! Layer 4
     if(noahmp401_struc(n)%noahmp401(t)%tslb(4)+delta4.gt.MIN_THRESHOLD .and.&
          noahmp401_struc(n)%noahmp401(t)%tslb(4)+delta4.lt.&
            MAX_THRESHOLD) then
        update_flag(gid) = update_flag(gid).and.(.true.)
        update_flag_tile(t) = update_flag_tile(t).and.(.true.)
     else
        update_flag(gid) = update_flag(gid).and.(.false.)
        update_flag_tile(t) = update_flag_tile(t).and.(.false.)
     endif
  enddo

!-----------------------------------------------------------------------------------------
! MN create new flag: if update falg for 50% of the ensemble members is true 
! then update the state variabels 
!-----------------------------------------------------------------------------------------
  update_flag_ens = .True.
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row) 
      flag_tmp=update_flag_tile(i:i+LIS_rc%nensem(n)-1)
      pcount = COUNT(flag_tmp) ! Counts the number of .TRUE. elements
      if (pcount.lt.LIS_rc%nensem(n)*0.5) then   ! 50%
         update_flag_ens(gid)= .False.
      endif
      update_flag_new(gid)= update_flag(gid).or.update_flag_ens(gid)  ! new flag
  enddo 
  
#if 0   
if(i.eq.66) then !i=66  ! --> domain's center  1376
  if(LIS_rc%hr.eq.12) then
     write(2001,'(I4, 2x, 3(I2,x), 2x, 23(L1,2x))'),&
          i, LIS_rc%mo, LIS_rc%da, LIS_rc%hr,update_flag_tile&
          ((i-1)*LIS_rc%nensem(n)+1:(i)*LIS_rc%nensem(n)),&
          update_flag_ens(i), update_flag_new(i), update_flag(i) 
  endif !mn
  endif
#endif 
  
  ! update step
  ! loop over grid points 
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row) 
     
     !if(update_flag(gid)) then
     if(update_flag_new(gid)) then 
!-----------------------------------------------------------------------------------------
       ! 1- update the states
       ! 1-1- if update flag for tile is TRUE --> apply the DA update    
       ! 1-2- if update flag for tile is FALSE --> resample the states  
       ! 2- adjust the states
!-----------------------------------------------------------------------------------------
        ! store update value for  cases that update_flag_tile & update_flag_new are TRUE
        ! update_flag_tile = TRUE --> means met the both min and max threshold 
      
        tmp1 = LIS_rc%udef
        tmp2 = LIS_rc%udef
        tmp3 = LIS_rc%udef
        tmp4 = LIS_rc%udef
        tmp5 = LIS_rc%udef
	!icount = 1
        do m=1,LIS_rc%nensem(n)
           t = i+m-1
           !t = (i-1)*LIS_rc%nensem(n)+m
 
	  if(update_flag_tile(t)) then
          
           tmp1(m) = soilt1(t)
           tmp2(m) = soilt2(t)
           tmp3(m) = soilt3(t)
           tmp4(m) = soilt4(t)
           tmp5(m) = tskin(t)
           !icount = icount + 1 
          endif
        enddo
        
        MaxEnsST1 = -10000
        MaxEnsST2 = -10000
        MaxEnsST3 = -10000
        MaxEnsST4 = -10000
        MaxEnsST5 = -10000

        MinEnsST1 = 10000
        MinEnsST2 = 10000
        MinEnsST3 = 10000
        MinEnsST4 = 10000
        MinEnsST5 = 10000

        do m=1,LIS_rc%nensem(n)
           if(tmp1(m).ne.LIS_rc%udef) then 
              MaxEnsST1 = max(MaxEnsST1, tmp1(m))
              MaxEnsST2 = max(MaxEnsST2, tmp2(m))
              MaxEnsST3 = max(MaxEnsST3, tmp3(m))
              MaxEnsST4 = max(MaxEnsST4, tmp4(m))
              MaxEnsST5 = max(MaxEnsST5, tmp5(m))

              MinEnsST1 = min(MinEnsST1, tmp1(m))
              MinEnsST2 = min(MinEnsST2, tmp2(m))
              MinEnsST3 = min(MinEnsST3, tmp3(m))
              MinEnsST4 = min(MinEnsST4, tmp4(m))
              MinEnsST5 = min(MinEnsST5, tmp5(m))
              
           endif
        enddo

        
        ! loop over tile       
        do m=1,LIS_rc%nensem(n)
           t = i+m-1           
           !t = (i-1)*LIS_rc%nensem(n)+m
           
           ! MN check update status for each tile  
           if(update_flag_tile(t)) then
              
              delta1 = soilt1(t)-noahmp401_struc(n)%noahmp401(t)%tslb(1)
              delta2 = soilt2(t)-noahmp401_struc(n)%noahmp401(t)%tslb(2)
              delta3 = soilt3(t)-noahmp401_struc(n)%noahmp401(t)%tslb(3)
              delta4 = soilt4(t)-noahmp401_struc(n)%noahmp401(t)%tslb(4)
              deltats = tskin(t)-noahmp401_struc(n)%noahmp401(t)%tsk

              noahmp401_struc(n)%noahmp401(t)%tslb(1) = noahmp401_struc(n)%noahmp401(t)%tslb(1)+&
                   delta1

              if(soilt1(t).lt.MIN_THRESHOLD) then
                 print*, 'setsoilt1 ',t,soilt1(t)
                 stop
              endif
              if(noahmp401_struc(n)%noahmp401(t)%tslb(2)+delta2.gt.MIN_THRESHOLD .and.&
                   noahmp401_struc(n)%noahmp401(t)%tslb(2)+delta2.lt.MAX_THRESHOLD) then
                 noahmp401_struc(n)%noahmp401(t)%tslb(2) = noahmp401_struc(n)%noahmp401(t)%tslb(2)+&
                      delta2

                 if(soilt2(t).lt.MIN_THRESHOLD) then
                    print*, 'setsoilt2 ',t,soilt2(t)
                    stop
                 endif
              endif
  
              if(noahmp401_struc(n)%noahmp401(t)%tslb(3)+delta3.gt.MIN_THRESHOLD .and.&
                   noahmp401_struc(n)%noahmp401(t)%tslb(3)+delta3.lt.MAX_THRESHOLD) then
                 noahmp401_struc(n)%noahmp401(t)%tslb(3) = noahmp401_struc(n)%noahmp401(t)%tslb(3)+&
                      delta3

                 if(soilt3(t).lt.MIN_THRESHOLD) then
                    print*, 'setsoilt3 ',t,soilt3(t)
                    stop
                 endif
              endif

              if(noahmp401_struc(n)%noahmp401(t)%tslb(4)+delta4.gt.MIN_THRESHOLD .and.&
                   noahmp401_struc(n)%noahmp401(t)%tslb(4)+delta4.lt.MAX_THRESHOLD) then
                 noahmp401_struc(n)%noahmp401(t)%tslb(4) = noahmp401_struc(n)%noahmp401(t)%tslb(4)+&
                      delta4

                 if(soilt4(t).lt.MIN_THRESHOLD) then
                    print*, 'setsoilt4 ',t,soilt4(t)
                    stop
                 endif
              endif

              if(noahmp401_struc(n)%noahmp401(t)%tsk+deltats.gt.MIN_THRESHOLD .and.&
                  noahmp401_struc(n)%noahmp401(t)%tsk+deltats.lt.MAX_THRESHOLD) then
              noahmp401_struc(n)%noahmp401(t)%tsk = noahmp401_struc(n)%noahmp401(t)%tsk+&
                      deltats
              
              if(tskin(t).lt.MIN_THRESHOLD) then
                  print*, 'settskin ',t,tskin(t)
                  stop
              endif
              endif
              
              
!-----------------------------------------------------------------------------------------              
              ! randomly resample the soil temperature from [MIN_THRESHOLD,
!               Max value from DA @ that time step]
!-----------------------------------------------------------------------------------------
           else 
          
!-----------------------------------------------------------------------------------------  
! set the soil moisture to the ensemble mean  
!-----------------------------------------------------------------------------------------
              
              ! use mean value
              ! Assume sh2o = smc (i.e. ice content=0) 
              stc_tmp = (MaxEnsST1 - MinEnsST1)/2 + MinEnsST1
              noahmp401_struc(n)%noahmp401(t)%tslb(1) = stc_tmp
                          
              stc_tmp = (MaxEnsST2 - MinEnsST2)/2 + MinEnsST2
              noahmp401_struc(n)%noahmp401(t)%tslb(2) = stc_tmp
               
              stc_tmp = (MaxEnsST3 - MinEnsST3)/2 + MinEnsST3
              noahmp401_struc(n)%noahmp401(t)%tslb(3) = stc_tmp
              
              stc_tmp = (MaxEnsST4 - MinEnsST4)/2 + MinEnsST4
              noahmp401_struc(n)%noahmp401(t)%tslb(4) = stc_tmp
    
              stc_tmp = (MaxEnsST5 - MinEnsST5)/2 + MinEnsST5
              noahmp401_struc(n)%noahmp401(t)%tsk = stc_tmp
              
 
           endif ! flag for each tile

        enddo ! loop over tile
       
     else ! if update_flag_new(gid) is FALSE   
        if(LIS_rc%pert_bias_corr.eq.1) then           
           !--------------------------------------------------------------------------
           ! if no update is made, then we need to readjust the ensemble if pert bias
           ! correction is turned on because the forcing perturbations may cause 
           ! biases to persist. 
           !--------------------------------------------------------------------------
           bounds_violation = .true. 
           nIter = 0
           ens_flag = .true. 
           
           do while(bounds_violation) 
              niter = niter + 1
              !t_unpert = i*LIS_rc%nensem(n)
	      t_unpert = i+LIS_rc%nensem(n)-1
              do j=1,4
                 delta(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                     t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    
                    if(m.ne.LIS_rc%nensem(n)) then 
                       delta(j) = delta(j) + &
                            (noahmp401_struc(n)%noahmp401(t)%tslb(j) - &
                            noahmp401_struc(n)%noahmp401(t_unpert)%tslb(j))
                    endif
                    
                 enddo
              enddo
              delta(j) = 0.0
              do m=1,LIS_rc%nensem(n)-1
                t = i+m-1
   
                if(m.ne.LIS_rc%nensem(n)) then
                    delta(j) = delta(j) + &
                    (noahmp401_struc(n)%noahmp401(t)%tsk - &
                        noahmp401_struc(n)%noahmp401(t_unpert)%tsk)
                endif
   
              enddo
              
              do j=1,4
                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                     t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    
                    tmpval = noahmp401_struc(n)%noahmp401(t)%tslb(j) - &
                         delta(j)
                    if(tmpval.le.MAX_THRESHOLD) then
                       noahmp401_struc(n)%noahmp401(t)%tslb(j) = &
                            max(noahmp401_struc(n)%noahmp401(t_unpert)%tslb(j), &
                            MIN_THRESHOLD)
                       ens_flag(m) = .false. 
                    elseif(tmpval.ge.MIN_THRESHOLD) then
                       noahmp401_struc(n)%noahmp401(t)%tslb(j) = &
                            min(noahmp401_struc(n)%noahmp401(t_unpert)%tslb(j), &
                            MIN_THRESHOLD)
                       ens_flag(m) = .false. 
                    endif
                 enddo

                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                     t = i+m-1
                     tmpval = noahmp401_struc(n)%noahmp401(t)%tsk - delta(j)
                     if(tmpval.le.MAX_THRESHOLD) then
                        noahmp401_struc(n)%noahmp401(t)%tsk = &
                            max(noahmp401_struc(n)%noahmp401(t_unpert)%tsk, &
                            MIN_THRESHOLD)
                        ens_flag(m) = .false.
                     elseif(tmpval.ge.MIN_THRESHOLD) then
                        noahmp401_struc(n)%noahmp401(t)%tsk = &
                            min(noahmp401_struc(n)%noahmp401(t_unpert)%tsk,&
                            MIN_THRESHOLD)
                        ens_flag(m) = .false.
                     endif
                 enddo

              enddo
              
              !--------------------------------------------------------------------------
              ! Recalculate the deltas and adjust the ensemble
              !--------------------------------------------------------------------------
              do j=1,4
                 delta(j) = 0.0
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    if(m.ne.LIS_rc%nensem(n)) then 
                       delta(j) = delta(j) + &
                            (noahmp401_struc(n)%noahmp401(t)%tslb(j) - &
                            noahmp401_struc(n)%noahmp401(t_unpert)%tslb(j))
                    endif
                 enddo
              enddo
              delta(j) = 0.0
              do m=1,LIS_rc%nensem(n)-1
                 t = i+m-1
                 !t = (i-1)*LIS_rc%nensem(n)+m
                 if(m.ne.LIS_rc%nensem(n)) then
                     delta(j) = delta(j) + &
                         (noahmp401_struc(n)%noahmp401(t)%tsk - &
                         noahmp401_struc(n)%noahmp401(t_unpert)%tsk)
                 endif
              enddo
              
              do j=1,4
                 delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
                 do m=1,LIS_rc%nensem(n)-1
                    t = i+m-1
                    
                    if(ens_flag(m)) then 
                       tmpval = noahmp401_struc(n)%noahmp401(t)%tslb(j) - &
                            delta(j)
                       
                       if(.not.(tmpval.le.MIN_THRESHOLD .or.&
                            tmpval.gt.(MAX_THRESHOLD))) then 
                          
                          noahmp401_struc(n)%noahmp401(t)%tslb(j) = &
                               noahmp401_struc(n)%noahmp401(t)%tslb(j) - delta(j)
                         
                          bounds_violation = .false.
                       endif
                    endif
                    
                    tmpval = noahmp401_struc(n)%noahmp401(t)%tslb(j)
                    
                    if(tmpval.le.MIN_THRESHOLD .or.&
                         tmpval.gt.(MAX_THRESHOLD)) then 
                       bounds_violation = .true. 
                    else
                       bounds_violation = .false.
                    endif
                 enddo
              enddo

              delta(j) =delta(j)/(LIS_rc%nensem(n)-1)
              do m=1,LIS_rc%nensem(n)-1
                 t = i+m-1
                 
                 if(ens_flag(m)) then
                     tmpval = noahmp401_struc(n)%noahmp401(t)%tsk - delta(j)
                     if(.not.(tmpval.le.MIN_THRESHOLD .or.&
                         tmpval.gt.(MAX_THRESHOLD))) then
                         
                         noahmp401_struc(n)%noahmp401(t)%tsk = &
                             noahmp401_struc(n)%noahmp401(t)%tsk - delta(j)
                         
                         bounds_violation = .false.
                     endif
                 endif
                 
                 tmpval = noahmp401_struc(n)%noahmp401(t)%tsk
                 
                 if(tmpval.le.MIN_THRESHOLD .or.&
                         tmpval.gt.(MAX_THRESHOLD)) then
                     bounds_violation = .true.
                 else
                     bounds_violation = .false.
                 endif
              enddo
              
              if(nIter.gt.10.and.bounds_violation) then 
                 !--------------------------------------------------------------------------
                 ! All else fails, set to the bounds
                 !--------------------------------------------------------------------------
                 
                 write(LIS_logunit,*) '[ERR] Ensemble structure violates physical bounds '
                 write(LIS_logunit,*) '[ERR] Please adjust the perturbation settings ..'

                 do j=1,4
                    do m=1,LIS_rc%nensem(n)
                       t = i+m-1
                       !t = (i-1)*LIS_rc%nensem(n)+m
                       
                       if(noahmp401_struc(n)%noahmp401(t)%tslb(j).gt.MAX_THRESHOLD) then
                          noahmp401_struc(n)%noahmp401(t)%tslb(j) = MAX_THRESHOLD
                       endif
                       
                       if(noahmp401_struc(n)%noahmp401(t)%tslb(j).lt.MIN_THRESHOLD) then
                          noahmp401_struc(n)%noahmp401(t)%tslb(j) = MIN_THRESHOLD
                       endif

                    enddo
                 enddo
                 do m=1,LIS_rc%nensem(n)
                    t = i+m-1
                    !t = (i-1)*LIS_rc%nensem(n)+m
                    
                    if(noahmp401_struc(n)%noahmp401(t)%tsk.gt.MAX_THRESHOLD) then
                        noahmp401_struc(n)%noahmp401(t)%tsk = MAX_THRESHOLD
                    endif
                    
                    if(noahmp401_struc(n)%noahmp401(t)%tsk.lt.MIN_THRESHOLD) then
                        noahmp401_struc(n)%noahmp401(t)%tsk = MIN_THRESHOLD
                    endif
                 
                 enddo
              endif          
           end do
        endif        
     endif
  enddo


end subroutine NoahMP401_setsoiltemp

