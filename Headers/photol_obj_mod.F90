!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: photol_obj_mod.F90
!
! !DESCRIPTION: Module PHOTOL\_OBJ\_MOD contains the derived type used to
!  define the Photolysis object for GEOS-Chem that is included as a member
!  of the State_Chm object
!\\
!\\
! !INTERFACE:
!
MODULE Photol_Obj_Mod
!
! USES:
!
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_Photol_Obj
  PUBLIC :: Cleanup_Photol_Obj
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Grid State
  !=========================================================================
  TYPE, PUBLIC :: PhotolState

     !----------------------------------------
     ! User-defined photolysis fields
     !----------------------------------------
     !CHARACTER(LEN=255) :: tmp       ! desc
     !REAL(fp)           :: tmp       ! desc
     INTEGER             :: W_        ! desc
     !LOGICAL            :: tmp       ! desc

     !----------------------------------------
     ! Photolysis fields computed in ___
     !----------------------------------------
     !INTEGER            :: tmp       ! desc

     ! Arrays
     !REAL(fp),  POINTER :: tmp (:,:) ! desc

  END TYPE PhotolState
!
! !REMARKS:
!
! !REVISION HISTORY:
!  28 Nov 2022 - E. Lundgren- Initial version, based on state_grid_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Photol_Obj
!
! !DESCRIPTION: Subroutine INIT\_PHOTOL\_OBJ allocates and initializes
! the Photol object
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Photol_Obj( Input_Opt, Photol, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)  :: Input_Opt ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(PhotolState), POINTER     :: Photol    ! Photolysis state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  28 Nov 2022 - E. Lundgren- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    !INTEGER :: AS

    !======================================================================
    ! Allocate and initialize module variables
    !======================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! NOTE: State_Grid%GlobalXMid and State_Grid%GlobalYMid are allocated
    ! in gc_grid_mod.F90 after computing State_Grid%GlobalNX and
    ! State_Grid%GlobalNY

    !ALLOCATE( Photol, STAT=RC )
    !CALL GC_CheckVar( 'Photol', 0, RC )
    !IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize
    Photol%W_   = 18


  END SUBROUTINE Init_Photol_Obj
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Photol_Obj
!
! !DESCRIPTION: Subroutine CLEANUP\_STATE\_PHOTOL deallocates all fields
!  of the photolysis state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Photol_Obj( Photol, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(PhotolState), POINTER     :: Photol ! Obj for photolysis state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC            ! Return code
!
! !REVISION HISTORY:
!  28 Nov 2022 - E. Lundgren- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC      = GC_SUCCESS

    !=======================================================================
    ! Deallocate arrays
    !=======================================================================
    ! Will need to change this to just do arrays etc
    IF ( ASSOCIATED( Photol ) ) THEN
       !DEALLOCATE( Photol, STAT=RC )
       !CALL GC_CheckVar( 'Photol', 2, RC )
       !IF ( RC /= GC_SUCCESS ) RETURN
       Photol => NULL()
    ENDIF

  END SUBROUTINE Cleanup_Photol_Obj
!EOC
END MODULE Photol_Obj_Mod
