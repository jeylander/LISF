!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_qcsoiltemp
! \label{NoahMP401_qcsoiltemp}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari; Modified for NoahMP401
! 21 Mar 2021: John Eylander (ERDC); Modified for soil temperature assimilation
!
! !INTERFACE:
subroutine NoahMP401_qcsoiltemp(n, LSM_State)

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
  type(ESMF_Field)       :: tslb1Field
!  type(ESMF_Field)       :: sm2Field
!  type(ESMF_Field)       :: sm3Field
!  type(ESMF_Field)       :: sm4Field
  integer                :: t
  integer                :: status
  real, pointer          :: soilt1(:)
!  real, pointer          :: soilm2(:)
!  real, pointer          :: soilm3(:)
!  real, pointer          :: soilm4(:)
  real                   :: stmax1!,smmax2,smmax3,smmax4
  real                   :: stmin1!,smmin2,smmin3,smmin4
 
  call ESMF_StateGet(LSM_State,"Soil Temperature Layer 1",tslb1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 1 failed in NoahMP401_qcsoiltemp")
 
  call ESMF_FieldGet(tslb1Field,localDE=0,farrayPtr=soilT1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Temperature Layer 1 failed in NoahMP401_qcsoiltemp")

  call ESMF_AttributeGet(tslb1Field,"Max Value",stmax1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Max Value failed in NoahMP401_qcsoiltemp")

  call ESMF_AttributeGet(tslb1Field,"Min Value",stmin1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Min Value failed in NoahMP401_qcsoiltemp")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(soilt1(t).gt.stmax1) soilt1(t) = stmax1
     if(soilt1(t).lt.stmin1) soilt1(t) = stmin1
  enddo

end subroutine NoahMP401_qcsoiltemp

