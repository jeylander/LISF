!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP401_write_soiltemp
! \label{NoahMP401_write_soiltemp}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 15 Dec 2018: Mahdi Navari; Modified for NoahMP401
! 21 Mar 2021:   John Eylander (ERDC); ISCCP DA soil/skin temp added
!
! !INTERFACE:
subroutine NoahMP401_write_soiltemp(ftn,n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use NoahMP401_lsmMod
  use LIS_historyMod, only : LIS_writevar_restart
  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: ftn
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
  integer :: t
  real, allocatable    :: tmp(:)

  allocate(tmp(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
   tmp(t) = noahmp401_struc(n)%noahmp401(t)%tsk
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = noahmp401_struc(n)%noahmp401(t)%tslb(1)
  enddo  
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = noahmp401_struc(n)%noahmp401(t)%tslb(2)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = noahmp401_struc(n)%noahmp401(t)%tslb(3)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tmp(t) = noahmp401_struc(n)%noahmp401(t)%tslb(4)
  enddo
  call LIS_writevar_restart(ftn,n,1,tmp)
  deallocate(tmp)

end subroutine NoahMP401_write_soiltemp

