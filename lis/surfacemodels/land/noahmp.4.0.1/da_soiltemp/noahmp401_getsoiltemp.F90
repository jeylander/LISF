!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_getsoiltemp
! \label{NoahMP401_getsoiltemp}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari; Modified for NoahMP401
! 21 Mar 2021: John Eylander (ERDC); Modified for soil temperature assimilation

! !INTERFACE:
subroutine NoahMP401_getsoiltemp(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use NoahMP401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soil/skin temperature related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: tskinField
  type(ESMF_Field)       :: st1Field
  type(ESMF_Field)       :: st2Field
  type(ESMF_Field)       :: st3Field
  type(ESMF_Field)       :: st4Field
  integer                :: t
  integer                :: status
  real, pointer          :: tskin(:)
  real, pointer          :: soilt1(:)
  real, pointer          :: soilt2(:)
  real, pointer          :: soilt3(:)
  real, pointer          :: soilt4(:)
  character*100          :: lsm_state_objs(4)

  call ESMF_StateGet(LSM_State,"Skin Temperature",tskinField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for tskin in NoahMP401_getsoiltemp')
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 1",st1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for st1 in NoahMP401_getsoiltemp')
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 2",st2Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for st2 in NoahMP401_getsoiltemp')
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 3",st3Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for st3 in NoahMP401_getsoiltemp')
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 4",st4Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for st4 in NoahMP401_getsoiltemp')

  call ESMF_FieldGet(tskinField,localDE=0,farrayPtr=tskin,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for tskin in NoahMP401_getsoiltemp')
  call ESMF_FieldGet(st1Field,localDE=0,farrayPtr=soilt1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for st1 in NoahMP401_getsoiltemp')
  call ESMF_FieldGet(st2Field,localDE=0,farrayPtr=soilt2,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for st2 in NoahMP401_getsoiltemp')
  call ESMF_FieldGet(st3Field,localDE=0,farrayPtr=soilt3,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for st3 in NoahMP401_getsoiltemp')
  call ESMF_FieldGet(st4Field,localDE=0,farrayPtr=soilt4,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for st4 in NoahMP401_getsoiltemp')


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tskin(t) = noahmp401_struc(n)%noahmp401(t)%tsk
     soilt1(t) = noahmp401_struc(n)%noahmp401(t)%tslb(1)
     soilt2(t) = noahmp401_struc(n)%noahmp401(t)%tslb(2)
     soilt3(t) = noahmp401_struc(n)%noahmp401(t)%tslb(3)
     soilt4(t) = noahmp401_struc(n)%noahmp401(t)%tslb(4)
  enddo

end subroutine NoahMP401_getsoiltemp

