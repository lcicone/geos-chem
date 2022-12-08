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
  ! Parameters used for allocation

  !=========================================================================
  ! Derived type for Grid State
  !=========================================================================
  TYPE, PUBLIC :: PhotolState

     ! NOTE: we might want to make some of the parameters input
     ! values rather than hard-code them during initialization.
     ! Several of the them are physical constants that could be
     ! moved to physconstants. And some of these might be able to
     ! be removed if only used in Cloud-J routines and they are defined
     ! in cloud-j.

     !---------------------------------------------------
     ! Fields that are also defined in Cloud-J
     !---------------------------------------------------

     ! Scalars used to allocate array sizes
     INTEGER  :: JVN_      ! Max # J-values
     INTEGER  :: AN_       ! Number of separate aerosols per layer
     INTEGER  :: WX_       ! # WLs in input file
     INTEGER  :: W_        ! # WL bins
     INTEGER  :: M_        ! # Gaussian points used (must be 4 in Fast-JX)
     INTEGER  :: A_        ! # input aerosol/cloud Mie data sets

     ! Other scalars
     INTEGER  :: L_        ! Number of CTM layers
     INTEGER  :: L1_       ! Number of CTM layer edges
     INTEGER  :: L2_       ! Number of levels with both edges and mid-pt
     INTEGER  :: JVL_      ! Vertical J-value levels
     INTEGER  :: JXL_      ! # levels in J-values array
     INTEGER  :: JXL1_     ! ??
     INTEGER  :: JXL2_     ! Max # levels in basic Fast-JX/Cloud-J grid
     INTEGER  :: X_        ! Max # input X-section data sets
     INTEGER  :: N_        ! # levels in Mie scattering arrays
     INTEGER  :: M2_       ! ??
     INTEGER  :: NJX       ! ??
     INTEGER  :: NW1       ! ??
     INTEGER  :: NW2       ! ??
     INTEGER  :: NAA       ! # categories for scattering phase functions

     ! Physical constants
     ! (Candidates for physconst or use from Cloud-J)
     REAL(fp) :: ZZHT      ! Scale height [cm]
     REAL(fp) :: RAD       ! Radius of Earth [cm]
     REAL(fp) :: ATAU      ! Heating rate (factor increase between levels)
     REAL(fp) :: ATAU0     ! Minimum heating rate

     ! Scalar set outside of this module
     INTEGER  :: NRATJ     ! # photolysis rxns in chemistry (.LE. JVN_)

     ! Arrays
     CHARACTER(LEN=6), ALLOCATABLE :: TITLEJX(:)   ! Input cross sections title
     CHARACTER(LEN=1), ALLOCATABLE :: SQQ    (:)   ! Input cross sections flag
     INTEGER,  ALLOCATABLE :: JIND   (:)   ! Index mapping J-values onto rates
     REAL(fp), ALLOCATABLE :: LQQ    (:)   ! Categorical interpolation options
     REAL(fp), ALLOCATABLE :: WBIN   (:)   ! Boundaries of WL bins
     REAL(fp), ALLOCATABLE :: WL     (:)   ! Centers of WL bins (eff WL)
     REAL(fp), ALLOCATABLE :: FL     (:)   ! Solar flx incident on TOA [cm-2s-1]
     REAL(fp), ALLOCATABLE :: QRAYL  (:)   ! Rayleigh params (eff Xs) [cm2]
     REAL(fp), ALLOCATABLE :: EMU    (:)   ! 4 Gauss pts ??
     REAL(fp), ALLOCATABLE :: WT     (:)   ! 4 Gauss pts ??
     REAL(fp), ALLOCATABLE :: JFACTA (:)   ! Multiplication factor for J-values
     REAL(fp), ALLOCATABLE :: JLABEL (:)   ! J-value label in main chem model
     REAL(fp), ALLOCATABLE :: SAA    (:,:) ! Single scattering albedo
     REAL(fp), ALLOCATABLE :: QAA    (:,:) ! Aerosol scatting phase fnctns
     REAL(fp), ALLOCATABLE :: WAA    (:,:) ! WLa for supplied phase functions 
     REAL(fp), ALLOCATABLE :: RAA    (:,:) ! Aerosol type effective radius
     REAL(fp), ALLOCATABLE :: QO2    (:,:) ! O2 X-sections
     REAL(fp), ALLOCATABLE :: QO3    (:,:) ! O3 X-sections
     REAL(fp), ALLOCATABLE :: Q1D    (:,:) ! O3 => O(1D) quantum yield
     REAL(fp), ALLOCATABLE :: PAA    (:,:,:) ! Phase fnctn (first 8 terms)

     !--------------------------------------------------
     ! Fields that are not in Cloud-J
     !--------------------------------------------------

     ! Scalars
     INTEGER  :: IND999      ! Index in RAA & QAA of 999 nm
     INTEGER  :: JTAUMX      ! max # divisions

     ! Photo-reaction flags for reactions adjusted in PhotRate_Adj
     INTEGER  :: RXN_O2      ! O2  + jv --> O   + O
     INTEGER  :: RXN_O3_1    ! O3  + hv --> O2  + O
     INTEGER  :: RXN_O3_2    ! O3  + hv --> O2  + O(1D)
     INTEGER  :: RXN_H2SO4   ! SO4 + hv --> SO2 + 2OH
     INTEGER  :: RXN_NO2     ! NO2 + hv --> NO  + O
     INTEGER  :: RXN_JHNO3   ! HNO3 + hv --> OH + NO2
     INTEGER  :: RXN_JNITSa  ! NITs  + hv --> HNO2
     INTEGER  :: RXN_JNITSb  ! NITs  + hv --> NO2
     INTEGER  :: RXN_JNITa   ! NIT + hv --> HNO2
     INTEGER  :: RXN_JNITb   ! NIT + hv --> NO2
     INTEGER  :: RXN_NO      ! For ucx_mod
     INTEGER  :: RXN_NO3     ! For ucx_mod
     INTEGER  :: RXN_N2O     ! For ucx_mod
     INTEGER  :: RXN_BrO     ! For Hg chem
     INTEGER  :: RXN_ClO     ! For Hg chem

     ! Arrays
     CHARACTER(LEN=10), ALLOCATABLE :: RNAMES(:)  ! Photolysis spc names
     CHARACTER(LEN=80), ALLOCATABLE :: TITLEAA(:) ! Scattering data title 

     INTEGER,  ALLOCATABLE :: BRANCH     (:) ! Photolysis spc branches
     INTEGER,  ALLOCATABLE :: RINDEX     (:) ! GC to UCI spc name index mapping
     INTEGER,  ALLOCATABLE :: GC_Photo_Id(:) ! GC id per photolysis species
     INTEGER,  ALLOCATABLE :: MIEDX      (:) ! Interface indices for GC/FJX spc

     REAL(fp), ALLOCATABLE :: UVXFACTOR(:) ! Photons/cm2s -> W/m2 conv factors
     REAL(fp), ALLOCATABLE :: QAA_AOD  (:) ! Single scattering albedo        
     REAL(fp), ALLOCATABLE :: WAA_AOD  (:) ! Aerosol scattering phase fnctns 
     REAL(fp), ALLOCATABLE :: PAA_AOD  (:) ! WLs for supplied phase functions
     REAL(fp), ALLOCATABLE :: RAA_AOD  (:) ! Phase fnctn (first 8 terms)     
     REAL(fp), ALLOCATABLE :: SAA_AOD  (:) ! Aerosol type effective radius

     REAL(fp), ALLOCATABLE :: TQQ      (:,:)       ! Temp for X-sections
     REAL(fp), ALLOCATABLE :: QQQ      (:,:,:)     ! Xs in each WL bin [cm2]
     REAL(fp), ALLOCATABLE :: TREF     (:,:,:)     ! Temp reference profile
     REAL(fp), ALLOCATABLE :: OREF     (:,:,:)     ! Ozone reference profile
     REAL(fp), ALLOCATABLE :: ZPJ      (:,:,:,:)   ! J-values

     !--------------------------------------------------
     ! Fields for RRTMG and optical depth diagnostics
     !--------------------------------------------------

     ! RRTMG scalars

     ! Scalars used for allocating arrays
     INTEGER :: NWVAA        ! # LUT wavelengths
     INTEGER :: NSPAA        ! # LUT species
     INTEGER :: NRAA         ! # LUT aerosol sizes    

     ! Other scalars
     INTEGER :: NWVAA0       ! # non-RRTMG wavelengths
     INTEGER :: NWVAART      ! # RRTMG wavelengths
     INTEGER :: NALBD        ! ??                     
     INTEGER :: NEMISS       ! ??                     
     INTEGER :: NASPECRAD    ! # RRTMG aerosol species
     INTEGER :: NSPECRAD     ! # RRTMG aerosol+gas species
     INTEGER :: NWVREQUIRED  ! # WLs needed for interpolation
     INTEGER :: NRTWVREQUIRED! # WLs needed for RT interpolation
     INTEGER :: IWV1000      ! WL index for 1000 nm

     ! RRTMG allocatable arrays
     INTEGER, ALLOCATABLE :: SPECMASK     (:)     ! binary switches for spc flux
     INTEGER, ALLOCATABLE :: IWVREQUIRED  (:)     ! WL indexes for interpolation
     INTEGER, ALLOCATABLE :: IRTWVREQUIRED(:)     ! WL indexes for RT interp
     INTEGER, ALLOCATABLE :: IWVSELECT    (:,:)   ! Indexes of requested WLs
     INTEGER, ALLOCATABLE :: IRTWVSELECT  (:,:)   ! Indexes of requested RT WLs
     INTEGER, ALLOCATABLE :: IRHARR       (:,:,:) ! Relative humidity indices

     REAL*8,  ALLOCATABLE :: ACOEF_WV  (:)   ! Coeffs for WL interpolation
     REAL*8,  ALLOCATABLE :: BCOEF_WV  (:)   ! Coeffs for WL interpolation
     REAL*8,  ALLOCATABLE :: CCOEF_WV  (:)   ! Coeffs for WL interpolation
     REAL*8,  ALLOCATABLE :: ACOEF_RTWV(:)   ! Coeffs for RT WL interpolation
     REAL*8,  ALLOCATABLE :: BCOEF_RTWV(:)   ! Coeffs for RT WL interpolation
     REAL*8,  ALLOCATABLE :: CCOEF_RTWV(:)   ! Coeffs for RT WL interpolation
     REAL*8,  ALLOCATABLE :: WVAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: RHAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: RDAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: RWAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: SGAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: REAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: NRLAA     (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: NCMAA     (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: QQAA      (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: ALPHAA    (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: SSAA      (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: ASYMAA    (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: PHAA      (:,:,:,:) ! ??

     ! For optical depth diagnostics
     REAL(fp), ALLOCATABLE :: ISOPOD   (:,:,:,:)   ! Isoprene optical depth
     REAL(fp), ALLOCATABLE :: ODMDUST  (:,:,:,:,:) ! Dust optical depth
     REAL(fp), ALLOCATABLE :: ODAER    (:,:,:,:,:) ! Aerosol optical depth

#ifdef RRTMG
     REAL*8,  ALLOCATABLE :: RTODAER   (:,:,:,:,:) ! Optical dust
     REAL*8,  ALLOCATABLE :: RTSSAER   (:,:,:,:,:) ! ??
     REAL*8,  ALLOCATABLE :: RTASYMAER (:,:,:,:,:) ! ??
#endif

  END TYPE PhotolState
!
! !REMARKS:
! Acronyms used in this file (may appear as upper or lower-case):
!   CLDJ      : Cloud-J
!   CONV      : conversion
!   EFF       : effective
!   FJX       : Fast-JX
!   FLX       : flux
!   FNCTN     : function
!   GC        : GEOS-Chem
!   LUT       : look-up table
!   RT        : radiative transfer
!   RXN       : reaction
!   SPC       : species
!   TEMP      : temperature
!   TOA       : top-of-atmosphere
!   UCI       : University of California, Irvine
!   WL        : wavelegnth
!   Xs        : cross-section
!   X-section : cross-section
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
  SUBROUTINE Init_Photol_Obj( Input_Opt, State_Grid, Photol, RC )
!
! !USES:
!
    USE CMN_Size_Mod,   ONLY : NDUST, NAER
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(GrdState),    INTENT(IN)  :: State_Grid ! Grid object
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
    CHARACTER(LEN=255)   :: errMsg, thisLoc

    !======================================================================
    ! Allocate and initialize module variables
    !======================================================================

    ! Assume success
    RC        = GC_SUCCESS
    thisLoc   = ' -> at Init_Photol_Obj (in module Headers/photol_obj_mod.F90)'

    !---------------------------------------------------
    ! Photolysis fields that are also in Cloud-J
    !---------------------------------------------------

    ! Integer scalars

    ! Used in array dimensions
    Photol%JVN_ = 18  ! Max # J-values                               
    Photol%AN_  = 37  ! Number of separate aerosols per layer        
    Photol%WX_  = 18  ! # WLs in input file                          
    Photol%W_   = 18  ! # WL bins                                    
    Photol%M_   = 4   ! # Gaussian points used (must be 4 in Fast-JX)
    Photol%A_   = 56  ! # input aerosol/cloud Mie data sets          

    ! Other scalars
    Photol%L_   = State_Grid%NZ
    Photol%L1_  = State_Grid%NZ + 1
    Photol%L2_  = 2 * ( State_Grid%NZ + 1 )
    Photol%JVL_ = State_Grid%NZ
    Photol%JXL_ = State_Grid%NZ
    Photol%JXL1_= State_Grid%NZ + 1
    Photol%JXL2_= 2 * ( State_Grid%NZ + 1 ) ! ewl: check this

    Photol%X_   = 123
!ewl: consider making this a config input??
    Photol%N_ = 601
    Photol%M2_ = 2*Photol%M_

    ! These are set somewhere else? (ewl)
    Photol%NJX = 0
    Photol%NW1 = 0
    Photol%NW2 = 0
    Photol%NAA = 0

    ! Real(fp) scalars

!ewl: consider putting these in physconst if not already there
    Photol%ZZHT  = 5.e+5_fp    ! Scale height [cm]
    Photol%RAD   = 6375.e+5_fp ! Radius of Earth [cm]
    Photol%ATAU  = 1.120e+0_fp ! Heating rate (factor increase between layers)
    Photol%ATAU0 = 0.010e+0_fp ! Minimum heating rate

    ! Allocate arrays
    IF ( .not. Input_Opt%DryRun ) THEN

       ! Character arrays

       ! Photol%TITLEJX (:)
       ALLOCATE( Photol%TITLEJX(Photol%X_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array TITLEJX!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%TITLEJX = ''

       ! Photol%SQQ     (:)
       ALLOCATE( Photol%SQQ(Photol%X_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SQQ!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%SQQ = ''

       ! Integer arrays

       ! Photol%JIND(:)
       ALLOCATE( Photol%JIND(Photol%JVN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array JIND!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%JIND = 0

       ! Real(fp) arrays

       ! Photol%LQQ     (:)
       ALLOCATE( Photol%LQQ(Photol%X_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
           errMsg = 'Error allocating array LQQ!'
           CALL GC_Error( errMsg, RC, thisLoc )
           RETURN
        ENDIF
       Photol%LQQ = 0e+0_fp

       ! Photol%WBIN    (:)
       ALLOCATE( Photol%WBIN(Photol%WX_+1), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array WBIN!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%WBIN = 0e+0_fp

       ! Photol%WL      (:)
       ALLOCATE( Photol%WL(Photol%WX_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array WL!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%WL = 0e+0_fp

       ! Photol%FL      (:)
       ALLOCATE( Photol%FL(Photol%WX_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array FL!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%FL = 0e+0_fp

       ! Photol%QRAYL   (:)
       ALLOCATE( Photol%QRAYL(Photol%WX_+1), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array QRAYL!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%QRAYL = 0e+0_fp

       ! Photol%EMU     (:)
       ALLOCATE( Photol%EMU(Photol%M_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array EMU!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%EMU = 0e+0_fp

       ! Photol%WT      (:)
       ALLOCATE( Photol%WT(Photol%M_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array WT!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%WT = 0e+0_fp

       ! Photol%JFACTA   (:)
       ALLOCATE( Photol%JFACTA(Photol%JVN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array JFACTA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%JFACTA = 0e+0_fp

       ! Photol%JLABEL   (:)
       ALLOCATE( Photol%JLABEL(Photol%JVN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array JLABEL!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%JLABEL = 0e+0_fp

       ! Photol%SAA     (:,:)
       ALLOCATE( Photol%SAA(5,Photol%A_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%SAA = 0e+0_fp

       ! Photol%QAA     (:,:)
       ALLOCATE( Photol%QAA(5,Photol%A_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array QAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%QAA = 0e+0_fp

       ! Photol%WAA     (:,:)
       ALLOCATE( Photol%WAA(5,Photol%A_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array WAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%WAA = 0e+0_fp

       ! Photol%RAA     (:,:)
       ALLOCATE( Photol%RAA(5,Photol%A_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RAA = 0e+0_fp

       ! Photol%QO2     (:,:)
       ALLOCATE( Photol%QO2(Photol%WX_,3), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array QO2!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%QO2 = 0e+0_fp

       ! Photol%QO3     (:,:)
       ALLOCATE( Photol%QO3(Photol%WX_,3), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array QO3!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%QO3 = 0e+0_fp

       ! Photol%Q1D     (:,:)
       ALLOCATE( Photol%Q1D(Photol%WX_,3), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array Q1D!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%Q1D = 0e+0_fp

       ! Photol%PAA     (:,:,:)
       ALLOCATE( Photol%PAA(8,5,Photol%A_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array PAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%PAA = 0e+0_fp

    ENDIF

    !--------------------------------------------------
    ! General photolysis fields that are not in Cloud-J
    !--------------------------------------------------

    ! Integer scalars

    Photol%IND999 = 5
    Photol%JTAUMX = (Photol%N_ - 4*Photol%JXL_) / 2

    Photol%RXN_O2     = -1
    Photol%RXN_O3_1   = -1
    Photol%RXN_O3_2   = -1
    Photol%RXN_H2SO4  = -1
    Photol%RXN_NO2    = -1
    Photol%RXN_JHNO3  = -1
    Photol%RXN_JNITSa = -1
    Photol%RXN_JNITSb = -1
    Photol%RXN_JNITa  = -1
    Photol%RXN_JNITb  = -1
    Photol%RXN_NO     = -1
    Photol%RXN_NO3    = -1
    Photol%RXN_N2O    = -1
    Photol%RXN_BrO    = -1
    Photol%RXN_ClO    = -1

    ! Allocate arrays
    IF ( .not. Input_Opt%DryRun ) THEN

       ! Character arrays

       ! Photol%RNAMES      (:)
       ALLOCATE( Photol%RNAMES(Photol%JVN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RNAMES!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RNAMES = ''

       ! Photol%TITLEAA     (:)
       ALLOCATE( Photol%TITLEAA(Photol%A_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array TITLEAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%TITLEAA = ''

       ! Integer arrays

       ! Photol%BRANCH      (:)
       ALLOCATE( Photol%BRANCH(Photol%JVN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array BRANCH!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%BRANCH = 0

       ! Photol%RINDEX      (:)
       ALLOCATE( Photol%RINDEX(Photol%JVN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RINDEX!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RINDEX = 0

       ! Photol%GC_Photo_Id (:)
       ALLOCATE( Photol%GC_Photo_Id(Photol%JVN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array GC_Photo_Id!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%GC_Photo_Id = 0

       ! Photol%MIEDX       (:)
       ALLOCATE( Photol%MIEDX(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array MIEDX!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%MIEDX = 0

       ! Real(fp) arrays

       ! Photol%UVXFACTOR(:)
       ALLOCATE( Photol%UVXFACTOR(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array UVXFACTOR!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%UVXFACTOR = 0e+0_fp

       ! Photol%QAA_AOD  (:)
       ALLOCATE( Photol%QAA_AOD(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array QAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%QAA_AOD = 0e+0_fp

       ! Photol%WAA_AOD  (:)
       ALLOCATE( Photol%WAA_AOD(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array WAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%WAA_AOD = 0e+0_fp

       ! Photol%PAA_AOD  (:)
       ALLOCATE( Photol%PAA_AOD(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array PAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%PAA_AOD = 0e+0_fp

       ! Photol%RAA_AOD  (:)
       ALLOCATE( Photol%RAA_AOD(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RAA_AOD = 0e+0_fp

       ! Photol%SAA_AOD  (:)
       ALLOCATE( Photol%SAA_AOD(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%SAA_AOD = 0e+0_fp

       ! Photol%TQQ      (:,:)
       ALLOCATE( Photol%TQQ(3,Photol%X_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array TQQ!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%TQQ = 0e+0_fp

       ! Photol%QQQ      (:,:,:)
       ALLOCATE( Photol%QQQ(Photol%WX_,3,Photol%X_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array QQQ!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%QQQ = 0e+0_fp

       ! Photol%TREF     (:,:,:)
       ALLOCATE( Photol%TREF(51,18,12), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array TREF!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%TREF = 0e+0_fp

       ! Photol%OREF     (:,:,:)
       ALLOCATE( Photol%OREF(51,18,12), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array OREF!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%OREF = 0e+0_fp

       ! Photol%ZPJ      (:,:,:,:)
       ALLOCATE( Photol%ZPJ( State_Grid%NZ, Photol%JVN_, State_Grid%NX, &
                 State_Grid%NY), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ZPJ!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%ZPJ = 0e+0_fp

    ENDIF

    !--------------------------------------------------
    ! Fields for RRTMG and optical depth diags (none in Cloud-J)
    !--------------------------------------------------
    
    ! RRTMG integer scalars

    ! Scalars used for allocating arrays
    Photol%NWVAA     = 41           ! # LUT wavelengths      
    Photol%NSPAA     = 8            ! # LUT species
    Photol%NRAA      = 7            ! # LUT aerosol sizes    

    ! Other scalars
    Photol%NWVAA0    = 11           ! # non-RRTMG wavelengths
    Photol%NWVAART   = Photol%NWVAA-Photol%NWVAA0 ! # RRTMG wavelengths    
    Photol%NALBD     = 2            ! ??                     
    Photol%NEMISS    = 16           ! ??                     
    Photol%NASPECRAD = 16           ! # RT aerosol species
    Photol%NSPECRAD  = 18           ! # RT aerosol+gas species

    ! These are set elsewhere? (ewl)
    Photol%NWVREQUIRED = 0
    Photol%NRTWVREQUIRED = 0
    Photol%IWV1000 = 0

    ! Allocate arrays
    IF ( .not. Input_Opt%DryRun ) THEN
    
       ! RRTMG integer arrays

       ! Photol%SPECMASK     (:)
       ALLOCATE( Photol%SPECMASK(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SPECMASK!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%SPECMASK = 0
    
       ! Photol%IWVREQUIRED  (:)
       ALLOCATE( Photol%IWVREQUIRED(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IWVREQUIRED!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%IWVREQUIRED = 0
    
       ! Photol%IRTWVREQUIRED(:)
       ALLOCATE( Photol%IRTWVREQUIRED(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IRTWVREQUIRED!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF 
       Photol%IRTWVREQUIRED = 0
    
       ! Photol%IWVSELECT    (:,:)
       ALLOCATE( Photol%IWVSELECT(2,3), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IWVSELECT!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%IWVSELECT = 0
    
       ! Photol%IRTWVSELECT  (:,:)
       ALLOCATE( Photol%IRTWVSELECT(2,3), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IRTWVSELECT!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%IRTWVSELECT = 0
    
       ! Photol%IRHARR    (:,:,:)
       ALLOCATE( Photol%IRHARR( State_Grid%NX, State_Grid%NY, &
                                State_Grid%NZ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IRHARR!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%IRHARR = 0d0
    
       ! RRTMG real*8 arrays
    
       ! Photol%ACOEF_WV  (:)
       ALLOCATE( Photol%ACOEF_WV(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ACOEF_WV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%ACOEF_WV = 0d0
    
       ! Photol%BCOEF_WV  (:)
       ALLOCATE( Photol%BCOEF_WV(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array BCOEF_WV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%BCOEF_WV = 0d0
    
       ! Photol%CCOEF_WV  (:)
       ALLOCATE( Photol%CCOEF_WV(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array CCOEF_WV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%CCOEF_WV = 0d0
    
       ! Photol%ACOEF_RTWV(:)
       ALLOCATE( Photol%ACOEF_RTWV(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ACOEF_RTWV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%ACOEF_RTWV = 0d0
    
       ! Photol%BCOEF_RTWV(:)
       ALLOCATE( Photol%BCOEF_RTWV(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array BCOEF_RTWV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%BCOEF_RTWV = 0d0
    
       ! Photol%CCOEF_RTWV(:)
       ALLOCATE( Photol%CCOEF_RTWV(Photol%AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array CCOEF_RTWV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%CCOEF_RTWV = 0d0
    
       ! Photol%WVAA      (:,:)
       ALLOCATE( Photol%WVAA(Photol%NWVAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array WVAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%WVAA = 0d0
    
       ! Photol%RHAA      (:,:)
       ALLOCATE( Photol%RHAA(Photol%NRAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RHAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RHAA = 0d0
    
       ! Photol%RDAA      (:,:)!
       ALLOCATE( Photol%RDAA(Photol%NRAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RDAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RDAA = 0d0
    
       ! Photol%RWAA      (:,:)
       ALLOCATE( Photol%RWAA(Photol%NRAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RWAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RWAA = 0d0
    
       ! Photol%SGAA      (:,:)
       ALLOCATE( Photol%SGAA(Photol%NRAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SGAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%SGAA = 0d0
    
       ! Photol%REAA      (:,:)
       ALLOCATE( Photol%REAA(Photol%NRAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array REAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%REAA = 0d0
    
       ! Photol%NRLAA     (:,:,:)
       ALLOCATE( Photol%NRLAA(Photol%NWVAA,Photol%NRAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array NRLAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%NRLAA = 0d0
    
       ! Photol%NCMAA     (:,:,:)
       ALLOCATE( Photol%NCMAA(Photol%NWVAA,Photol%NRAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array NCMAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%NCMAA = 0d0
    
       ! Photol%QQAA      (:,:,:)
       ALLOCATE( Photol%QQAA(Photol%NWVAA,Photol%NRAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array QQAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%QQAA = 0d0
    
       ! Photol%ALPHAA    (:,:,:)
       ALLOCATE( Photol%ALPHAA(Photol%NWVAA,Photol%NRAA,Photol%NSPAA), &
                 STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ALPHAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%ALPHAA = 0d0
    
       ! Photol%SSAA      (:,:,:)
       ALLOCATE( Photol%SSAA(Photol%NWVAA,Photol%NRAA,Photol%NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SSAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%SSAA = 0d0
    
       ! Photol%ASYMAA    (:,:,:)
       ALLOCATE( Photol%ASYMAA(Photol%NWVAA,Photol%NRAA,Photol%NSPAA), &
                 STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ASYMAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%ASYMAA = 0d0
    
       ! Photol%PHAA      (:,:,:,:)
       ALLOCATE( Photol%PHAA(Photol%NWVAA,Photol%NRAA,Photol%NSPAA,8), &
                 STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array PHAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%PHAA = 0d0

       ! Photol%ISOPOD   (:,:,:,:)
       ALLOCATE( Photol%ISOPOD( State_Grid%NX, State_Grid%NY, &
                                State_Grid%NZ, Photol%NWVAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ISOPOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%ISOPOD = 0e+0_fp
       
       ! Photol%ODMDUST  (:,:,:,:,:)
       ALLOCATE( Photol%ODMDUST( State_Grid%NX, State_Grid%NY, &
                 State_Grid%NZ, Photol%NWVAA, NDUST), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ODMDUST!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%ODMDUST = 0e+0_fp
       
       ! Photol%ODAER    (:,:,:,:,:)
       ALLOCATE( Photol%ODAER( State_Grid%NX, State_Grid%NY, &
                 State_Grid%NZ, Photol%NWVAA, NAER), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ODAER!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%ODAER = 0e+0_fp

#ifdef RRTMG
       ! Photol%RTODAER   (:,:,:,:,:)
       ! +2 to split SNA into SU, NI and AM
       ALLOCATE( Photol%RTODAER( State_Grid%NX, State_Grid%NY, &
                 State_Grid%NZ, Photol%NWVAA,NAER+2+NDUST), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RTODAER!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RTODAER = 0d0

       ! Photol%RTSSAER   (:,:,:,:,:)
       ALLOCATE( Photol%RTSSAER( State_Grid%NX, State_Grid%NY, &
                 State_Grid%NZ, Photol%NWVAA, NAER+2+NDUST ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RTSSAER!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RTSSAER = 0d0
       
       ! Photol%RTASYMAER (:,:,:,:,:)
       ALLOCATE( Photol%RTASYMAER( State_Grid%NX, State_Grid%NY, &
                 State_Grid%NZ, Photol%NWVAA, NAER+2+NDUST ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RTASYMAER!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Photol%RTASYMAER = 0d0
#endif
     
    ENDIF

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
       IF (ALLOCATED(Photol%JIND          )) DEALLOCATE(Photol%JIND          )
       IF (ALLOCATED(Photol%LQQ           )) DEALLOCATE(Photol%LQQ           )
       IF (ALLOCATED(Photol%TITLEJX       )) DEALLOCATE(Photol%TITLEJX       )
       IF (ALLOCATED(Photol%SQQ           )) DEALLOCATE(Photol%SQQ           )
       IF (ALLOCATED(Photol%WBIN          )) DEALLOCATE(Photol%WBIN          )
       IF (ALLOCATED(Photol%WL            )) DEALLOCATE(Photol%WL            )
       IF (ALLOCATED(Photol%FL            )) DEALLOCATE(Photol%FL            )
       IF (ALLOCATED(Photol%EMU           )) DEALLOCATE(Photol%EMU           )
       IF (ALLOCATED(Photol%WT            )) DEALLOCATE(Photol%WT            )
       IF (ALLOCATED(Photol%JFACTA        )) DEALLOCATE(Photol%JFACTA        )
       IF (ALLOCATED(Photol%SAA           )) DEALLOCATE(Photol%SAA           )
       IF (ALLOCATED(Photol%QAA           )) DEALLOCATE(Photol%QAA           )
       IF (ALLOCATED(Photol%WAA           )) DEALLOCATE(Photol%WAA           )
       IF (ALLOCATED(Photol%RAA           )) DEALLOCATE(Photol%RAA           )
       IF (ALLOCATED(Photol%QO2           )) DEALLOCATE(Photol%QO2           )
       IF (ALLOCATED(Photol%QO3           )) DEALLOCATE(Photol%QO3           )
       IF (ALLOCATED(Photol%Q1D           )) DEALLOCATE(Photol%Q1D           )
       IF (ALLOCATED(Photol%PAA           )) DEALLOCATE(Photol%PAA           )
       IF (ALLOCATED(Photol%QRAYL         )) DEALLOCATE(Photol%QRAYL         )
       IF (ALLOCATED(Photol%RNAMES        )) DEALLOCATE(Photol%RNAMES        )
       IF (ALLOCATED(Photol%TITLEAA       )) DEALLOCATE(Photol%TITLEAA       )
       IF (ALLOCATED(Photol%BRANCH        )) DEALLOCATE(Photol%BRANCH        )
       IF (ALLOCATED(Photol%RINDEX        )) DEALLOCATE(Photol%RINDEX        )
       IF (ALLOCATED(Photol%GC_Photo_Id   )) DEALLOCATE(Photol%GC_Photo_Id   )
       IF (ALLOCATED(Photol%MIEDX         )) DEALLOCATE(Photol%MIEDX         )
       IF (ALLOCATED(Photol%UVXFACTOR     )) DEALLOCATE(Photol%UVXFACTOR     )
       IF (ALLOCATED(Photol%QAA_AOD       )) DEALLOCATE(Photol%QAA_AOD       )
       IF (ALLOCATED(Photol%WAA_AOD       )) DEALLOCATE(Photol%WAA_AOD       )
       IF (ALLOCATED(Photol%PAA_AOD       )) DEALLOCATE(Photol%PAA_AOD       )
       IF (ALLOCATED(Photol%RAA_AOD       )) DEALLOCATE(Photol%RAA_AOD       )
       IF (ALLOCATED(Photol%SAA_AOD       )) DEALLOCATE(Photol%SAA_AOD       )
       IF (ALLOCATED(Photol%TQQ           )) DEALLOCATE(Photol%TQQ           )
       IF (ALLOCATED(Photol%QQQ           )) DEALLOCATE(Photol%QQQ           )
       IF (ALLOCATED(Photol%TREF          )) DEALLOCATE(Photol%TREF          )
       IF (ALLOCATED(Photol%OREF          )) DEALLOCATE(Photol%OREF          )
       IF (ALLOCATED(Photol%ISOPOD        )) DEALLOCATE(Photol%ISOPOD        )
       IF (ALLOCATED(Photol%ZPJ           )) DEALLOCATE(Photol%ZPJ           )
       IF (ALLOCATED(Photol%ODMDUST       )) DEALLOCATE(Photol%ODMDUST       )
       IF (ALLOCATED(Photol%ODAER         )) DEALLOCATE(Photol%ODAER         )
       IF (ALLOCATED(Photol%SPECMASK      )) DEALLOCATE(Photol%SPECMASK      )
       IF (ALLOCATED(Photol%IWVREQUIRED   )) DEALLOCATE(Photol%IWVREQUIRED   )
       IF (ALLOCATED(Photol%IRTWVREQUIRED )) DEALLOCATE(Photol%IRTWVREQUIRED )
       IF (ALLOCATED(Photol%IWVSELECT     )) DEALLOCATE(Photol%IWVSELECT     )
       IF (ALLOCATED(Photol%IRTWVSELECT   )) DEALLOCATE(Photol%IRTWVSELECT   )
       IF (ALLOCATED(Photol%IRHARR        )) DEALLOCATE(Photol%IRHARR        )
       IF (ALLOCATED(Photol%ACOEF_WV      )) DEALLOCATE(Photol%ACOEF_WV      )
       IF (ALLOCATED(Photol%BCOEF_WV      )) DEALLOCATE(Photol%BCOEF_WV      )
       IF (ALLOCATED(Photol%CCOEF_WV      )) DEALLOCATE(Photol%CCOEF_WV      )
       IF (ALLOCATED(Photol%ACOEF_RTWV    )) DEALLOCATE(Photol%ACOEF_RTWV    )
       IF (ALLOCATED(Photol%BCOEF_RTWV    )) DEALLOCATE(Photol%BCOEF_RTWV    )
       IF (ALLOCATED(Photol%CCOEF_RTWV    )) DEALLOCATE(Photol%CCOEF_RTWV    )
       IF (ALLOCATED(Photol%WVAA          )) DEALLOCATE(Photol%WVAA          )
       IF (ALLOCATED(Photol%RHAA          )) DEALLOCATE(Photol%RHAA          )
       IF (ALLOCATED(Photol%RDAA          )) DEALLOCATE(Photol%RDAA          )
       IF (ALLOCATED(Photol%RWAA          )) DEALLOCATE(Photol%RWAA          )
       IF (ALLOCATED(Photol%SGAA          )) DEALLOCATE(Photol%SGAA          )
       IF (ALLOCATED(Photol%REAA          )) DEALLOCATE(Photol%REAA          )
       IF (ALLOCATED(Photol%NRLAA         )) DEALLOCATE(Photol%NRLAA         )
       IF (ALLOCATED(Photol%NCMAA         )) DEALLOCATE(Photol%NCMAA         )
       IF (ALLOCATED(Photol%QQAA          )) DEALLOCATE(Photol%QQAA          )
       IF (ALLOCATED(Photol%ALPHAA        )) DEALLOCATE(Photol%ALPHAA        )
       IF (ALLOCATED(Photol%SSAA          )) DEALLOCATE(Photol%SSAA          )
       IF (ALLOCATED(Photol%ASYMAA        )) DEALLOCATE(Photol%ASYMAA        )
       IF (ALLOCATED(Photol%PHAA          )) DEALLOCATE(Photol%PHAA          )
#ifdef RRTMG 
       IF (ALLOCATED(Photol%RTODAER       )) DEALLOCATE(Photol%RTODAER   )
       IF (ALLOCATED(Photol%RTSSAER       )) DEALLOCATE(Photol%RTSSAER   )
       IF (ALLOCATED(Photol%RTASYMAER     )) DEALLOCATE(Photol%RTASYMAER )
#endif

       Photol => NULL()
    ENDIF

  END SUBROUTINE Cleanup_Photol_Obj
!EOC
END MODULE Photol_Obj_Mod
