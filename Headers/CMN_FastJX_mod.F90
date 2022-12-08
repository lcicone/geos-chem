!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_FastJX_mod.F90
!
! !DESCRIPTION: Module CMN\_FastJX\_MOD contains parameters and global variables
!  used to interface between Harvard chemistry and UC-Irvine photolysis
!  programs (Fast-J/Fast-JX).
!\\
!\\
! !INTERFACE:
!
MODULE CMN_FastJX_MOD
!
! !USES:
!
  USE PRECISION_MOD      ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PUBLIC
!
! !DEFINED PARAMETERS:
!
  ! Required size of aerosol arrays
  INTEGER            :: L_             ! Number of CTM layers

  INTEGER            :: L1_            ! Number of CTM layer edges

  INTEGER            :: L2_            ! Number of levels in FJX grid that
                                       ! inc. both edges and mid-points

  INTEGER            :: JVL_           ! Vertical levels for J-values

  INTEGER, PARAMETER :: JVN_ = 166     ! Max number of J-values

  INTEGER            :: AN_            ! # of separate aerosols per layer
                                       ! Now set in Init_CMN_FastJX below


  !-----------------------------------------------------------------------
  ! variables used to map fast-JX J's onto CTM J's
  !-----------------------------------------------------------------------

  ! Multiplication factor for fast-JX calculated J
  REAL(fp)           :: JFACTA(JVN_)

  ! Index arrays that map Jvalue(j) onto rates
  INTEGER            :: JIND(JVN_)

  ! Mumber of Photolysis reactions in CTM chemistry, derived here NRATJ
  ! must be .le. JVN_
  INTEGER            :: NRATJ

  ! Label of J-value used in the main chem model
  CHARACTER*50       :: JLABEL(JVN_)

  ! JXL_: vertical(levels) dim for J-values computed within fast-JX
  INTEGER            :: JXL_
  INTEGER            :: JXL1_

  ! JXL2_: 2*JXL_ + 2 = mx no. levels in the basic Fast-JX grid (mid-level)
  INTEGER            :: JXL2_

  ! WX_  = dim = no. of wavelengths in input file
  INTEGER, PARAMETER :: WX_   = 18

  ! X_   = dim = max no. of X-section data sets (input data)
  ! SDE 2016-11-04: Increased from 75 to 123 for new halogen
  ! chemistry (iodine cross-sections)
  INTEGER, PARAMETER :: X_    = 123

  ! A_   = dim = no. of Aerosol/cloud Mie sets (input data)
  INTEGER, PARAMETER :: A_    = 56

  ! W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only
  ! Now set below in Init_CMN_FastJX (mps, 1/3/18)
  INTEGER            :: W_

  ! N_  = no. of levels in Mie scattering arrays
  !     = 2*NC+1 = 4*(L_+1) + 1 + 2*sum(JADDLV)
#ifdef MODEL_GEOS
!!!INTEGER, PARAMETER :: N_    = 601
!!!INTEGER, PARAMETER :: N_    = 1201
  INTEGER            :: N_
#else
  INTEGER, PARAMETER :: N_    = 601
#endif

  ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
  INTEGER, PARAMETER :: M_    = 4

  ! M2_ = 2*M_ = 8, replaces MFIT
  INTEGER, PARAMETER :: M2_   = 2*M_

  !-----------------------------------------------------------------------
  ! 4 Gauss pts = 8-stream
  REAL(fp), DIMENSION(M_), PARAMETER  :: &
       EMU = [.06943184420297e+0_fp, .33000947820757e+0_fp, &
              .66999052179243e+0_fp, .93056815579703e+0_fp]
  REAL(fp), DIMENSION(M_), PARAMETER  :: &
       WT  = [.17392742256873e+0_fp, .32607257743127e+0_fp, &
              .32607257743127e+0_fp, .17392742256873e+0_fp]
  !-----------------------------------------------------------------------

  ! ZZHT: scale height (cm)
  REAL(fp), PARAMETER  :: ZZHT   = 5.e+5_fp

  ! RAD: Radius of Earth (cm)
  REAL(fp), PARAMETER  :: RAD    = 6375.e+5_fp

  ! ATAU: heating rate (factor increase from one layer to the next)
  REAL(fp), PARAMETER  :: ATAU   = 1.120e+0_fp
#ifdef MODEL_GEOS
  !REAL(fp), PARAMETER  :: ATAU   = 1.180e+0_fp
#endif

  ! ATAU0: minimum heating rate
  REAL(fp), PARAMETER  :: ATAU0  = 0.010e+0_fp

  ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
  INTEGER              :: JTAUMX

  !-----------------------------------------------------------------------
  ! Variables in file 'FJX_spec.dat' (RD_XXX)
  !-----------------------------------------------------------------------

  ! WBIN: Boundaries of wavelength bins
  REAL(fp)             :: WBIN(WX_+1)

  ! WL: Centres of wavelength bins - 'effective wavelength'
  REAL(fp)             :: WL(WX_)

  ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
  REAL(fp)             :: FL(WX_)

  REAL(fp)             :: QO2(WX_,3)   ! QO2: O2 cross-sections
  REAL(fp)             :: QO3(WX_,3)   ! QO3: O3 cross-sections
  REAL(fp)             :: Q1D(WX_,3)   ! Q1D: O3 => O(1D) quantum yield

  ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
  REAL(fp)             :: QRAYL(WX_+1)

  ! LQQQ = 1, 2, or 3 to determine interpolation with T or P
  INTEGER              :: LQQ(X_)

  ! TITLEJX: Title for supplied cross sections, from 'FJX_spec.dat'
  CHARACTER*6          :: TITLEJX(X_)

  ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
  CHARACTER*1          :: SQQ(X_)

  ! TQQ: Temperature for supplied cross sections
  REAL(fp)             :: TQQ(3,X_)

  !-----------------------------------------------------------------------
  ! Variables in file 'jv_spec_mie.dat' (RD_MIE)
  !-----------------------------------------------------------------------

  ! QAA: Aerosol scattering phase functions
  REAL(fp)             :: QAA(5,A_)

  ! WAA: 5 Wavelengths for the supplied phase functions
  REAL(fp)             :: WAA(5,A_)

  ! PAA: Phase function: first 8 terms of expansion
  REAL(fp)             :: PAA(8,5,A_)

! RAA is in Cloud-J but with different dimensions. Use this in fast-jx for now.
  ! RAA: Effective radius associated with aerosol type
  REAL(fp)             :: RAA(5,A_)

  ! SAA: Single scattering albedo
  REAL(fp)             :: SAA(5,A_)

  ! NAA: Number of categories for scattering phase functions
  INTEGER              :: NAA

  !-----------------------------------------------------------------------
  ! Variables in file 'atmos_std.dat' (RD_PROF)
  !-----------------------------------------------------------------------

  ! TITLAA: Title for scattering data
  CHARACTER*80, DIMENSION(A_) :: TITLAA

  INTEGER NJX,NW1,NW2

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Cmn_FastJX
!
! !DESCRIPTION: Routine INIT\_CMN\_FastJX initializes quantities based on
!  the grid-independent size parameters.
!\\
!\\
! !INTERFACE:

  SUBROUTINE Init_CMN_FastJX( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! INIT_CMN_FastJX begins here!
    !=================================================================

    L_     = State_Grid%NZ ! Number of CTM layers
    L1_    = L_+1          ! Number of CTM layer edges
    L2_    = L1_*2         ! Number of levels in FJX grid that
                           ! inc. both edges and mid-points
    JVL_   = State_Grid%NZ ! Vertical levs for J-values

    JXL_   = State_Grid%NZ ! Vertical levs for J-values computed w/in Fast-JX
    JXL1_  = JXL_+1        ! Vertical levs edges for J-values
    JXL2_  = 2*JXL_+2      ! Max # levs in the basic Fast-JX grid (mid-level)

    JTAUMX = ( N_ - 4*JXL_ ) / 2  ! Maximum number of divisions ( i.e., may
                                  ! not get to ATAUMN)

#ifdef MODEL_GEOS
    ! N_  = no. of levels in Mie scattering arrays
    IF ( Input_Opt%LLFASTJX > 0 ) THEN
       N_ = Input_Opt%LLFASTJX
    ELSE
       N_ = 601
    ENDIF
#endif

    AN_       = 37  ! # of separate aerosols per layer; Including PSCs
    W_        = 18  ! # of wavelength bins

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE Init_CMN_FastJX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Cmn_FastJX
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_FastJX deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_CMN_FastJX( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Feb 2014 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !IF ( ALLOCATED( ZPJ       ) ) DEALLOCATE( ZPJ       )

    ! Return successfully
    RC = GC_SUCCESS

  END SUBROUTINE Cleanup_CMN_FastJX
!EOC
END MODULE CMN_FastJX_MOD
