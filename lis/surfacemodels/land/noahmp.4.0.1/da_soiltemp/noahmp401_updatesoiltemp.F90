!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_updatesoiltemp
!  \label{NoahMP401_updatesoiltemp}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari; Modified for NoahMP401
! 21 Mar 2021: John Eylander (ERDC); Modified for soil temperature assimilation
!
! !INTERFACE:
subroutine NoahMP401_updatesoiltemp(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use NoahMP401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP

  type(ESMF_Field)       :: tslb1Field
  type(ESMF_Field)       :: tslb2Field
  type(ESMF_Field)       :: tslb3Field
  type(ESMF_Field)       :: tslb4Field
  type(ESMF_Field)       :: tskinField
  type(ESMF_Field)       :: tslb1IncrField
  type(ESMF_Field)       :: tslb2IncrField
  type(ESMF_Field)       :: tslb3IncrField
  type(ESMF_Field)       :: tslb4IncrField
  type(ESMF_Field)       :: tskinIncrField

  real, pointer          :: tskin(:)
  real, pointer          :: soilt1(:)
  real, pointer          :: soilt2(:)
  real, pointer          :: soilt3(:)
  real, pointer          :: soilt4(:)
  real, pointer          :: tskin_Incr(:)
  real, pointer          :: soilT_Incr1(:)
  real, pointer          :: soilT_Incr2(:)
  real, pointer          :: soilT_Incr3(:)
  real, pointer          :: soilT_Incr4(:)
  integer                :: t,i,m
  integer                :: status

  call ESMF_StateGet(LSM_State,"Skin Temperature",tskinField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for tskin in NoahMP401_updatesoiltemp')

  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 1",tslb1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Temperature Layer 1 failed in NoahMP401_updatesoiltemp")
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 2",tslb2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Temperature Layer 2 failed in NoahMP401_updatesoiltemp")
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 3",tslb3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Temperature Layer 3 failed in NoahMP401_updatesoiltemp")
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 4",tslb4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Temperature Layer 4 failed in NoahMP401_updatesoiltemp")

  call ESMF_FieldGet(tskinField,localDE=0,farrayPtr=tskin,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for tskin in NoahMP401_updatesoiltemp')

  call ESMF_FieldGet(tslb1Field,localDE=0,farrayPtr=soilt1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 1 failed in NoahMP401_updatesoiltemp")
  call ESMF_FieldGet(tslb2Field,localDE=0,farrayPtr=soilt2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 2 failed in NoahMP401_updatesoiltemp")
  call ESMF_FieldGet(tslb3Field,localDE=0,farrayPtr=soilt3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 3 failed in NoahMP401_updatesoiltemp")
  call ESMF_FieldGet(tslb4Field,localDE=0,farrayPtr=soilt4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 4 failed in NoahMP401_updatesoiltemp")

  call ESMF_StateGet(LSM_Incr_State,"Skin Temperature",tskinIncrField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for tskin in NoahMP401_updatesoiltemp')

  call ESMF_StateGet(LSM_Incr_State,"Soil Temperature Layer 1",tslb1IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Temperature Layer 1 failed in NoahMP401_updatesoiltemp")
  call ESMF_StateGet(LSM_Incr_State,"Soil Temperature Layer 2",tslb2IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Temperature Layer 2 failed in NoahMP401_updatesoiltemp")
  call ESMF_StateGet(LSM_Incr_State,"Soil Temperature Layer 3",tslb3IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Temperature Layer 3 failed in NoahMP401_updatesoiltemp")
  call ESMF_StateGet(LSM_Incr_State,"Soil Temperature Layer 4",tslb4IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Temperature Layer 4 failed in NoahMP401_updatesoiltemp")

  call ESMF_FieldGet(tskinIncrField,localDE=0,farrayPtr=tskin_Incr,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for tskin in NoahMP401_updatesoiltemp')

  call ESMF_FieldGet(tslb1IncrField,localDE=0,farrayPtr=soilT_Incr1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 1 failed in NoahMP401_updatesoiltemp")
  call ESMF_FieldGet(tslb2IncrField,localDE=0,farrayPtr=soilT_Incr2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 2 failed in NoahMP401_updatesoiltemp")
  call ESMF_FieldGet(tslb3IncrField,localDE=0,farrayPtr=soilT_Incr3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 3 failed in NoahMP401_updatesoiltemp")
  call ESMF_FieldGet(tslb4IncrField,localDE=0,farrayPtr=soilT_Incr4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Temperature Layer 4 failed in NoahMP401_updatesoiltemp")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tskin(t) = tskin(t) + tskin_Incr(t)
     soilt1(t) = soilt1(t) + soilT_Incr1(t)
     soilt2(t) = soilt2(t) + soilT_Incr2(t)
     soilt3(t) = soilt3(t) + soilT_Incr3(t)
     soilt4(t) = soilt4(t) + soilT_Incr4(t)
  enddo
end subroutine NoahMP401_updatesoiltemp

