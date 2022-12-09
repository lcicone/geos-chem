!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fjx_interface_mod.F90
!
! !DESCRIPTION: Module FJX\_INTERFACE\_MOD contains routines and variables
!  for interfacing with the Fast-JX scheme (Prather et al) that calculates
!  photolysis rates. Current implementation is version 7.0a.
!\\
!\\
! !INTERFACE:
!
MODULE FJX_INTERFACE_MOD
!
! !USES:
!
  USE FJX_Mod          ! ewl new
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
#if defined( MODEL_CESM ) && defined( SPMD )
  USE MPISHORTHAND
  USE SPMD_UTILS
#endif

  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: INIT_FJX
  PUBLIC  :: FAST_JX
  PUBLIC  :: PHOTRATE_ADJ
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: GC_EXITC ! ewl new
  PRIVATE :: RD_PROF_NC
  PRIVATE :: RD_AOD
  PRIVATE :: CALC_AOD
  PRIVATE :: SET_PROF
  PRIVATE :: SET_AER
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: int_fjx
!
! !DESCRIPTION: Subroutine INIT\_FJX initializes Fast-JX variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_FJX( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE Charpak_Mod,    ONLY : CSTRIP
    USE CMN_FastJX_Mod, ONLY : JVN_, NJX, NRATJ, W_, WL
    USE CMN_FastJX_Mod, ONLY : TITLEJX, JLABEL, RNAMES, JFACTA
    USE CMN_Phot_Mod,   ONLY : GC_Photo_ID, UVXFACTOR
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE inquireMod,     ONLY : findFreeLUN
    USE PhysConstants,  ONLY : UVXPlanck, UVXCConst
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
#if defined( MODEL_CESM )
    USE UNITS,          ONLY : freeUnit
#endif

!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: notDryRun
    INTEGER            :: JXUNIT, J, NJXX, PhotoId
    REAL(fp)           :: ND64MULT

    ! Strings
    CHARACTER(LEN=6)   :: TITLEJXX(JVN_)
    CHARACTER(LEN=255) :: DATA_DIR
    CHARACTER(LEN=255) :: FILENAME
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! ewl: Added due to code taken out of rd_fs_fx
    INTEGER            :: K
    CHARACTER(LEN=50 ) :: TEXT

    !=================================================================
    ! INIT_FJX begins here!
    !=================================================================

    ! Initialize
    RC          = GC_SUCCESS
    notDryRun   = ( .not. Input_Opt%DryRun )
    ErrMsg      = ''
    ThisLoc     = ' -> at Init_FJX (in module GeosCore/fjx_interface_mod.F90)'

#ifdef CLOUDJ
    ! For now, if using cloud-j we are still using this init function
    ! Overwrite certain values from Cloud-J for consistency with GEOS-Chem
    ZZHT  = 5.e+5_fp
    RAD   = 6375.e+5_fp
    ATAU  = 1.120e+0_fp
    ATAU0 = 0.010e+0_fp
#endif

    ! Skip these opterations when running in dry-run mode
    IF ( notDryRun ) THEN

       ! Print info
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) ' Initializing Fast-JX v7.0 standalone CTM code.'

          if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
             ErrMsg =  ' INIT_FJX: invalid no. wavelengths'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          endif
       ENDIF

#if defined( MODEL_CESM )
       IF ( Input_Opt%amIRoot ) THEN
          JXUNIT = findFreeLUN()
       ELSE
          JXUNIT = 0
       ENDIF
#else
       ! Get a free LUN
       JXUNIT = findFreeLUN()
#endif

    ENDIF

    ! Define data directory for FAST-JX input
    DATA_DIR = TRIM( Input_Opt%FAST_JX_DIR )

    !=====================================================================
    ! Read in fast-J X-sections (spectral data)
    !=====================================================================
    FILENAME = TRIM( DATA_DIR ) // 'FJX_spec.dat'

    ! Read file, or just print filename if we are in dry-run mode
    CALL RD_XXX( Input_Opt%amIRoot, Input_Opt%DryRun, JXUNIT, &
                 TRIM( FILENAME ), RC)

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_XXX"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Compute factors for UV flux diagnostics
    IF ( notDryRun ) THEN
       IF ( State_Diag%Archive_UVFluxNet      .or. &
            State_Diag%Archive_UVFluxDirect   .or. &
            State_Diag%Archive_UVFluxDiffuse ) THEN
          UVXFACTOR = 0e+0_fp
          ND64MULT  = UVXPLANCK*UVXCCONST*1.0e+13_fp
          DO J = 1, W_
             UVXFACTOR(J) = ND64MULT/WL(J)
          ENDDO
       ENDIF
    ENDIF

    !=====================================================================
    ! Read in 5-wavelength scattering data
    ! (or just print file name if in dry-run mode)
    !=====================================================================
    FILENAME = TRIM( DATA_DIR ) // 'jv_spec_mie.dat'

    ! Read data
    CALL RD_MIE( Input_Opt%amIRoot, Input_Opt%DryRun, Input_Opt%LBRC, &
                 JXUNIT, TRIM( FILENAME ), RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_MIE"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=====================================================================
    ! Read in AOD data
    ! (or just print file name if in dry-run mode)
    !=====================================================================
    CALL RD_AOD( JXUNIT, Input_Opt, RC )

!!ewl: took this out of RD_AOD
    ! Only do the following if we are not running in dry-run mode
    IF ( .not. Input_Opt%DryRun ) THEN

       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, * ) 'Optics read for all wavelengths successfully'
       ENDIF

       ! Now calculate the required wavelengths in the LUT to calculate
       ! the requested AOD
       CALL CALC_AOD( Input_Opt )
    ENDIF
!!ewl

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_AOD"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Set up MIEDX array to interpret between GC and FJX aerosol indexing
    IF ( notDryRun ) THEN
       CALL SET_AER( Input_Opt )
    ENDIF

    !=====================================================================
    ! Read in T & O3 climatology used to fill e.g. upper layers
    ! or if O3 not calc.
    !=====================================================================
    CALL RD_PROF_NC( Input_Opt, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Rd_Prof_Nc"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Skip if not in dry-run mode
    IF ( notDryRun ) THEN
       NJXX = NJX
       do J = 1,NJX
          TITLEJXX(J) = TITLEJX(J)
       enddo
    ENDIF

    !=====================================================================
    ! Read in photolysis rates used in chemistry code and mapping onto
    ! FJX J's CTM call:  read in J-values names and link to fast-JX names
    !=====================================================================
    FILENAME = TRIM( DATA_DIR ) // 'FJX_j2j.dat'

    ! Read mapping information
    CALL RD_JS_JX( Input_Opt%amIRoot, Input_Opt%DryRun, JXUNIT, &
                   TRIM( FILENAME ), TITLEJXX, NJXX, RC )

    ! Store # of photolysis reactions in state_chm
    State_Chm%Photol%NRatJ = NRatJ

! ewl: bring out of RD_JS_JX
    !========================================================================
    ! Flag special reactions that will be later adjusted by
    ! routine PHOTRATE_ADJ (called from FlexChem)
    !========================================================================

    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       ! Loop over all photolysis reactions
       DO K = 1, NRATJ

          ! Strip all blanks from the reactants and products list
          TEXT = JLABEL(K)
          CALL CSTRIP( TEXT )

          !IF ( Input_Opt%amIRoot ) THEN
          !   WRITE(*,*) K, TRIM( TEXT )
          !ENDIF

          ! Look for certain reactions
          SELECT CASE( TRIM( TEXT ) )
             CASE( 'O2PHOTONOO' )
                State_Chm%Photol%RXN_O2 = K     ! O2 + hv -> O + O
             CASE( 'O3PHOTONO2O' )
                State_Chm%Photol%RXN_O3_1 = K   ! O3 + hv -> O2 + O
             CASE( 'O3PHOTONO2O(1D)' )
                State_Chm%Photol%RXN_O3_2 = K   ! O3 + hv -> O2 + O(1D)
             CASE( 'SO4PHOTONSO2OHOH' )
                State_Chm%Photol%RXN_H2SO4 = K  ! SO4 + hv -> SO2 + OH + OH
             CASE( 'NO2PHOTONNOO' )
                State_Chm%Photol%RXN_NO2 = K    ! NO2 + hv -> NO + O
             CASE( 'NOPHOTONNO' )
                State_Chm%Photol%RXN_NO = K     ! NO + hv -> N + O
             CASE( 'NO3PHOTONNO2O' )
                State_Chm%Photol%RXN_NO3 = K    ! NO3 + hv -> NO2 + O
             CASE( 'N2OPHOTONN2O' )
                State_Chm%Photol%RXN_N2O = K    ! N2O + hv -> N2 + O
             CASE( 'NITsPHOTONHNO2' )
                State_Chm%Photol%RXN_JNITSa = K ! NITs + hv -> HNO2
             CASE( 'NITsPHOTONNO2' )
                State_Chm%Photol%RXN_JNITSb = K ! NITs + hv -> NO2
             CASE( 'NITPHOTONHNO2' )
                State_Chm%Photol%RXN_JNITa = K  ! NIT + hv -> HNO2
             CASE( 'NITPHOTONNO2' )
                State_Chm%Photol%RXN_JNITb = K  ! NIT + hv -> NO2
             CASE( 'HNO3PHOTONNO2OH' )
                State_Chm%Photol%RXN_JHNO3 = K  ! HNO3 + hv = OH + NO2
             CASE DEFAULT
                ! Nothing
          END SELECT
       ENDDO

       !---------------------------------------------------------------------
       ! Error check the various rxn flags
       !---------------------------------------------------------------------
       IF ( State_Chm%Photol%RXN_O2 < 0 ) THEN
          ErrMsg = 'Could not find rxn O2 + hv -> O + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_O3_1 < 0 ) THEN
          ErrMsg = 'Could not find rxn O3 + hv -> O2 + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_O3_2 < 0 ) THEN
          ErrMsg = 'Could not find rxn O3 + hv -> O2 + O(1D)'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
       ENDIF

       IF ( State_Chm%Photol%RXN_NO2 < 0 ) THEN
          ErrMsg = 'Could not find rxn NO2 + hv -> NO + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_NO2 < 0 ) THEN
          ErrMsg = 'Could not find rxn NO2 + hv -> NO + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_JNITSa < 0 ) THEN
          ErrMsg = 'Could not find rxn NITS + hv -> HNO2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_JNITSb < 0 ) THEN
          ErrMsg = 'Could not find rxn NITS + hv -> NO2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_JNITa < 0 ) THEN
          ErrMsg = 'Could not find rxn NIT + hv -> HNO2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_JNITb < 0 ) THEN
          ErrMsg = 'Could not find rxn NIT + hv -> NO2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_H2SO4  < 0 ) THEN
          ErrMsg = 'Could not find rxn SO4 + hv -> SO2 + OH + OH!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_NO3 < 0 ) THEN
          ErrMsg = 'Could not find rxn NO3 + hv -> NO2 + O'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_NO < 0 ) THEN
          ErrMsg = 'Could not find rxn NO + hv -> O + N'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( State_Chm%Photol%RXN_N2O < 0 ) THEN
          ErrMsg = 'Could not find rxn N2O + hv -> N2 + O(1D)'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! Print out saved rxn flags for fullchem simulations
       !---------------------------------------------------------------------
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 100 ) REPEAT( '=', 79 )
          WRITE( 6, 110 )
          WRITE( 6, 120 ) State_Chm%Photol%RXN_O2
          WRITE( 6, 130 ) State_Chm%Photol%RXN_O3_1
          WRITE( 6, 140 ) State_Chm%Photol%RXN_O3_2
          WRITE( 6, 180 ) State_Chm%Photol%RXN_JNITSa
          WRITE( 6, 190 ) State_Chm%Photol%RXN_JNITSb
          WRITE( 6, 200 ) State_Chm%Photol%RXN_JNITa
          WRITE( 6, 210 ) State_Chm%Photol%RXN_JNITb
          WRITE( 6, 160 ) State_Chm%Photol%RXN_H2SO4
          WRITE( 6, 170 ) State_Chm%Photol%RXN_NO2
          WRITE( 6, 100 ) REPEAT( '=', 79 )
       ENDIF
    ENDIF

    !========================================================================
    ! Flag reactions for diagnostics (only in Hg chem)
    !========================================================================
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
        ! Loop over all photolysis reactions
        DO K = 1, NRATJ

           ! Strip all blanks from the reactants and products list
           TEXT = JLABEL(K)
           CALL CSTRIP( TEXT )

           ! Look for certain reactions
           SELECT CASE( TRIM( TEXT ) )
              CASE( 'O3PHOTONO2O' )
                 State_Chm%Photol%RXN_O3_1 = K ! O3 + hv -> O2 + O
              CASE( 'O3PHOTONO2O(1D)' )
                 State_Chm%Photol%RXN_O3_2 = K ! O3 + hv -> O2 + O(1D)
              CASE( 'NO2PHOTONNOO' )
                 State_Chm%Photol%RXN_NO2 = K  ! NO2 + hv -> NO + O
              CASE( 'BrOPHOTONBrO' )
                 State_Chm%Photol%RXN_BrO = K  ! BrO + hv -> Br + O
              CASE( 'ClOPHOTONClO' )
                 State_Chm%Photol%RXN_ClO = K  ! ClO + hv -> Cl + O
              CASE DEFAULT
                 ! Nothing
           END SELECT
        ENDDO

        !--------------------------------------------------------------------
        ! Error check the various rxn flags
        !--------------------------------------------------------------------
        IF ( State_Chm%Photol%RXN_O3_1 < 0 ) THEN
           ErrMsg = 'Could not find rxn O3 + hv -> O2 + O'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

        IF ( State_Chm%Photol%RXN_O3_2 < 0 ) THEN
           ErrMsg = 'Could not find rxn O3 + hv -> O2 + O(1D) #1'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

        IF ( State_Chm%Photol%RXN_NO2 < 0 ) THEN
           ErrMsg = 'Could not find rxn NO2 + hv -> NO + O'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

        IF ( State_Chm%Photol%RXN_BrO < 0 ) THEN
           ErrMsg = 'Could not find rxn BrO + hv -> Br + O'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

        IF ( State_Chm%Photol%RXN_ClO < 0 ) THEN
           ErrMsg = 'Could not find rxn ClO + hv -> Cl + O'
           CALL GC_Error( ErrMsg, RC, ThisLoc )
           RETURN
        ENDIF

       !---------------------------------------------------------------------
       ! Print out saved rxn flags for Hg simulation
       !---------------------------------------------------------------------
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 100 ) REPEAT( '=', 79 )
          WRITE( 6, 110 )
          WRITE( 6, 130 ) State_Chm%Photol%RXN_O3_1
          WRITE( 6, 140 ) State_Chm%Photol%RXN_O3_2
          WRITE( 6, 170 ) State_Chm%Photol%RXN_NO2
          WRITE( 6, 220 ) State_Chm%Photol%RXN_BrO
          WRITE( 6, 230 ) State_Chm%Photol%RXN_ClO
          WRITE( 6, 100 ) REPEAT( '=', 79 )
       ENDIF
    ENDIF

    ! FORMAT statements
100 FORMAT( a                                                 )
110 FORMAT( 'Photo rxn flags saved for use in PHOTRATE_ADJ:', / )
120 FORMAT( 'RXN_O2     [ O2   + hv -> O + O         ]  =  ', i5 )
130 FORMAT( 'RXN_O3_1   [ O3   + hv -> O2 + O        ]  =  ', i5 )
140 FORMAT( 'RXN_O3_2a  [ O3   + hv -> O2 + O(1D) #1 ]  =  ', i5 )
150 FORMAT( 'RXN_O3_2b  [ O3   + hv -> O2 + O(1D) #2 ]  =  ', i5 )
160 FORMAT( 'RXN_H2SO4  [ SO4  + hv -> SO2 + OH + OH ]  =  ', i5 )
170 FORMAT( 'RXN_NO2    [ NO2  + hv -> NO + O        ]  =  ', i5 )
180 FORMAT( 'RXN_JNITSa [ NITS + hv -> HNO2          ]  =  ', i5 )
190 FORMAT( 'RXN_JNITSb [ NITS + hv -> NO2           ]  =  ', i5 )
200 FORMAT( 'RXN_JNITa  [ NIT  + hv -> HNO2          ]  =  ', i5 )
210 FORMAT( 'RXN_JNITb  [ NIT  + hv -> NO2           ]  =  ', i5 )
220 FORMAT( 'RXN_BrO    [ BrO  + hv -> Br + O        ]  =  ', i5 )
230 FORMAT( 'RXN_ClO    [ ClO  + hv -> Cl + O        ]  =  ', i5 )

!ewl end

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Rd_Js_Jx"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Skip further processing if we are in dry-run mode
    IF ( notDryRun ) THEN

       ! Get the GEOS-Chem photolysis index for each of the 1...JVN_ entries
       ! in the FJX_j2j.dat file.  We'll use this for the diagnostics.
       DO J = 1, JVN_

          IF ( J == State_Chm%Photol%Rxn_O3_2 ) THEN

             !------------------------------------------------------------
             ! O3 + hv = O + O(1D)
             !
             ! Save this as JO3_O1D in the nPhotol+1 slot
             !------------------------------------------------------------
             GC_Photo_Id(J) = State_Chm%nPhotol + 1

          ELSE IF ( J == State_Chm%Photol%Rxn_O3_1 ) THEN

             !------------------------------------------------------------
             ! O3 + hv -> O + O
             !
             ! Save this as JO3_O3P in the nPhotol+2 slot
             !-------------------------------------------------------------
             GC_Photo_Id(J) = State_Chm%nPhotol + 2

          ELSE

             !------------------------------------------------------------
             ! Everything else
             !
             ! Find the matching GEOS-Chem photolysis species number
             !------------------------------------------------------------
             GC_Photo_Id(J) = Ind_( RNAMES(J), 'P' )

          ENDIF

          ! Print the mapping
          IF ( Input_Opt%amIRoot ) THEN
             IF ( GC_Photo_Id(J) > 0 ) THEN
! ewl: changed format from 200 to 240
                WRITE(6, 240) RNAMES(J), J, GC_Photo_Id(J), JFACTA(J)
240             FORMAT( a10, ':', i7, 2x, i7, 2x, f7.4 )
             ENDIF
          ENDIF
       ENDDO

#if defined( MODEL_CESM )
       IF ( Input_Opt%amIRoot ) THEN
         CALL freeUnit(JXUnit)
       ENDIF
#endif

    ENDIF

  END SUBROUTINE INIT_FJX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fast_jx
!
! !DESCRIPTION: Subroutine FAST\_JX loops over longitude and latitude, and
!  calls PHOTO\_JX to compute J-Values for each column at every chemistry
!  time-step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FAST_JX( WLAOD, Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
!
! !USES:
!
    USE CMN_FastJX_Mod,     ONLY : A_, L_, L1_, W_, JVN_, JXL_, JXL1_
    USE CMN_FastJX_Mod,     ONLY : NRATJ, JIND, JFACTA, FL
    USE CMN_Phot_Mod,       ONLY : ZPJ, IRHARR, UVXFACTOR, IWV1000
    USE CMN_Phot_Mod,       ONLY : ODAER, ODMDUST
    USE CMN_SIZE_MOD,       ONLY : NDUST, NRH, NAER
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP, ALLOC_ERR
    USE ERROR_MOD,          ONLY : DEBUG_MSG
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR
    USE TIME_MOD,           ONLY : GET_TAU,   GET_YEAR
    USE TOMS_MOD,           ONLY : GET_OVERHEAD_O3

    IMPLICIT NONE

!==============================================================================
! Uncomment the appropriate #define statement to denote which of the
! available cloud overlap options that you wish to use.

!! Linear overlap
!#define USE_LINEAR_OVERLAP 1

! Approximate random overlap (balance between accuracy & speed)
#define USE_APPROX_RANDOM_OVERLAP 1

!! Maximum random cloud overlap (most computationally intensive)
!#define USE_MAXIMUM_RANDOM_OVERLAP 1
!==============================================================================
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: WLAOD       ! AOD calculated how?
                                                 ! (1: 550 nm, 0: 999 nm)
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
!  Parameter to choose cloud overlap algorithm:
!  ============================================================================
!  (1 ) OVERLAP (INTEGER) : 1 - Linear Approximation (used up to v7-04-12)
!                           2 - Approximate Random Overlap (default)
!                           3 - Maximum Random Overlap (computation intensive)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, SAVE :: LASTMONTH = -1
    INTEGER       :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR, L, N, J
    INTEGER       :: IOPT, LCHEM
    REAL(fp)      :: CSZA, PRES, SFCA, YLAT,  O3_TOMS, SZA, SOLF
    REAL(fp)      :: O3_CTM(State_Grid%NZ+1)
    REAL(fp)      :: T_CTM(State_Grid%NZ+1), OPTD(State_Grid%NZ)
    REAL(fp)      :: OPTDUST(State_Grid%NZ,NDUST)
    REAL(fp)      :: OPTAER(State_Grid%NZ,A_)

    ! Local variables for cloud overlap (hyl, phs)
    INTEGER       :: NUMB, KK, I
    INTEGER       :: INDIC(State_Grid%NZ+1)
    INTEGER       :: INDGEN(State_Grid%NZ+1)! = (/ (i,i=1,State_Grid%NZ+1) /)
    INTEGER       :: KBOT(State_Grid%NZ)
    INTEGER       :: KTOP(State_Grid%NZ)
    INTEGER       :: INDICATOR(State_Grid%NZ+2)
    REAL(fp)      :: FMAX(State_Grid%NZ)    ! maximum cloud fraction
                                              !  in a block, size can be to
                                              !  FIX(State_Grid%NZ)+1
    REAL(fp)      :: CLDF1D(State_Grid%NZ)
    REAL(fp)      :: ODNEW(State_Grid%NZ)
    REAL(fp)      :: P_CTM(State_Grid%NZ+2)
    REAL(fp)      :: AERX_COL(A_,L1_)
    REAL(fp)      :: T_CLIM(L1_)
    REAL(fp)      :: O3_CLIM(L1_)
    REAL(fp)      :: Z_CLIM(L1_+1)
    REAL(fp)      :: AIR_CLIM(L1_)
    REAL(fp)      :: VALJXX(L_,JVN_)

    LOGICAL       :: AOD999
    LOGICAL, SAVE :: FIRST = .true.
    LOGICAL       :: prtDebug

    ! Species ID flags
    INTEGER, SAVE :: id_O3

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! ewl: input to set_prof, but I question whether we need it
    REAL(fp) :: ODCLOUD_COL(L_)

    ! ewl: output from photo_jx for use in diagnostics
    REAL(fp), DIMENSION(W_)         :: FJBOT,FSBOT
    REAL(fp), DIMENSION(JXL1_,W_)   :: FLXD
    REAL(fp), DIMENSION(JXL_, W_)   :: FJFLX

    ! UVFlux* diagnostics
    REAL(fp) :: FDIRECT (JXL1_)
    REAL(fp) :: FDIFFUSE(JXL1_)
    REAL(fp) :: UVX_CONST
    INTEGER  :: S, K

    !=================================================================
    ! FAST_JX begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Fast_JX (in module GeosCore/fjx_interface_mod.F90)'
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot)

    ! Get day of year (0-365 or 0-366)
    DAY_OF_YR = GET_DAY_OF_YEAR()

    ! Get current month
    MONTH     = GET_MONTH()

    ! Get day of month
    DAY       = GET_DAY()

    ! Was AOD calculated at 999 nm or reference?
    AOD999    = ( WLAOD == 0 )

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF ( State_Diag%Archive_UVFluxDiffuse ) THEN
       State_Diag%UVFluxDiffuse = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_UVFluxDirect ) THEN
       State_Diag%UVFluxDirect = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_UVFluxNet ) THEN
       State_Diag%UVFluxNet = 0.0_f4
    ENDIF

    !-----------------------------------------------------------------
    ! Special handling for first-time setup
    !-----------------------------------------------------------------
    IF ( FIRST ) THEN

       ! Get the species ID for O3
       id_O3 = Ind_('O3')
       IF ( id_O3 < 0 ) THEN
          ErrMsg = 'O3 is not a defined species!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

#ifdef USE_MAXIMUM_RANDOM_OVERLAP
       ! Special setup only for max random overlap
       DO i = 1,State_Grid%NZ+1
          INDGEN(i) = i       !(/(i,i=1,State_Grid%NZ+1)/)
       ENDDO
#endif

       ! Reset first-time flag
       FIRST = .FALSE.

    ENDIF

    !=================================================================
    ! For each (NLON,NLAT) location, call subroutine PHOTO_JX (in a
    ! parallel loop to compute J-values for the entire column.
    ! J-values will be stored in the common-block variable ZPJ, and
    ! will be later accessed via function FJXFUNC.
    !=================================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( NLAT,    NLON,   YLAT,      CSZA,    L       ) &
    !$OMP PRIVATE( P_CTM ,  T_CTM,  SFCA,      O3_TOMS, O3_CTM  ) &
    !$OMP PRIVATE( LCHEM,   OPTAER, N,         IOPT,    J       ) &
    !$OMP PRIVATE( OPTDUST, OPTD,   CLDF1D                      ) &
#ifdef USE_MAXIMUM_RANDOM_OVERLAP
    !$OMP PRIVATE( FMAX,    KK,     NUMB,      KBOT             ) &
    !$OMP PRIVATE( KTOP     ODNEW,  INDICATOR, INDIC            ) &
#endif
    !$OMP PRIVATE( SZA, SOLF, ODCLOUD_COL                       ) &
    !$OMP PRIVATE( AERX_COL,  T_CLIM, O3_CLIM, Z_CLIM, AIR_CLIM ) &
    !$OMP PRIVATE( VALJXX,    FSBOT,  FJBOT,   FLXD,   FJFLX    ) &
    !$OMP SCHEDULE( DYNAMIC )

    ! Loop over latitudes and longitudes
    DO NLAT = 1, State_Grid%NY
    DO NLON = 1, State_Grid%NX

       ! Grid box latitude [degrees]
       YLAT = State_Grid%YMid( NLON, NLAT )

       ! Cosine of solar zenith angle [unitless] at (NLON,NLAT)
       CSZA = State_Met%SUNCOSmid(NLON,NLAT)

       ! Define the P array here
       DO L = 1, State_Grid%NZ+1
          P_CTM(L) = State_Met%PEDGE( NLON, NLAT, L )
       ENDDO

       ! Top edge of P_CTM is top of atmosphere (bmy, 2/13/07)
       P_CTM(State_Grid%NZ+2) = 0e+0_fp

       ! Temperature profile [K] at (NLON,NLAT)
       T_CTM(1:State_Grid%NZ) = State_Met%T( NLON, NLAT, 1:State_Grid%NZ)

       ! Top of atmosphere
       T_CTM(State_Grid%NZ+1) = T_CTM(State_Grid%NZ)

       ! Surface albedo [unitless] at (NLON,NLAT)
       SFCA = State_Met%UVALBEDO(NLON,NLAT)

       ! Overhead ozone column [DU] at (NLON, NLAT)
       ! These values are either from the met fields or TOMS/SBUV,
       ! depending on the settings in geoschem_config.yml
       O3_TOMS = GET_OVERHEAD_O3( State_Chm, NLON, NLAT )

       ! CTM ozone densities (molec/cm3) at (NLON, NLAT)
       O3_CTM = 0e+0_fp
       LCHEM  = State_Met%ChemGridLev(NLON,NLAT)
       DO L = 1, LCHEM
          O3_CTM(L) = State_Chm%Species(id_O3)%Conc(NLON,NLAT,L)
       ENDDO

       ! Aerosol OD profile [unitless] at (NLON,NLAT)
       ! and at 1000nm, IWV1000 (DAR)
       !OPTAER wants NAER*NRH values but ODAER is now NAER
       !use IRHARR to map to correct OPTAER bin (DAR 08/13)
       OPTAER = 0.0e+0_fp
       DO N = 1, NAER
       DO L = 1, State_Grid%NZ
          IOPT = ( (N-1) * NRH ) + IRHARR(NLON,NLAT,L)
          OPTAER(L,IOPT) = ODAER(NLON,NLAT,L,IWV1000,N)
       ENDDO
       ENDDO
       DO N = 1, NDUST
       DO L = 1, State_Grid%NZ
          OPTDUST(L,N) = ODMDUST(NLON,NLAT,L,IWV1000,N)
       ENDDO
       ENDDO

       ! Mineral dust OD profile [unitless] at (NLON,NLAT)
       ! and at 1000nm, IWV1000 (DAR)
       !OPTDUST = ODMDUST(NLON,NLAT,:,IWV1000,:)

       ! Cloud OD profile [unitless] at (NLON,NLAT)
       OPTD = State_Met%OPTD(NLON,NLAT,1:State_Grid%NZ)

       !-----------------------------------------------------------
       !### If you want to exclude aerosol OD, mineral dust OD,
       !### or cloud OD, then uncomment the following lines:
       !OPTAER  = 0d0
       !OPTDUST = 0d0
       !OPTD(:)    = 0d0
       !-----------------------------------------------------------

#if defined( MODEL_GEOS )
       ! Initialize diagnostics arrays
       IF ( State_Diag%Archive_EXTRALNLEVS ) THEN
          State_Diag%EXTRALNLEVS(ILON,ILAT) = 0.0
       ENDIF
       IF ( State_Diag%Archive_EXTRALNITER ) THEN
          State_Diag%EXTRALNITER(ILON,ILAT) = 0.0
       ENDIF
#endif

       if (State_Grid%NZ+1 .gt. JXL1_) then
          call GC_EXITC(' PHOTO_JX: not enough levels in JX')
       endif

       ! Input conversion (SDE 03/29/13)
       ! Calculate solar zenith angle (degrees)
       CALL SOLAR_JX(DAY_OF_YR,CSZA,SZA,SOLF)
       
       ! check for dark conditions SZA > 98.0 deg => tan ht = 63 km or
       !                                 99.0                 80 km
       if (SZA .gt. 98.e+0_fp) cycle

! ewl: better to put these options in geoschem_config.yml?
#if defined( USE_LINEAR_OVERLAP )
       !===========================================================
       ! %%%% CLOUD OVERLAP: LINEAR ASSUMPTION %%%%
       !
       ! Directly use OPTDEPTH = TAUCLD * CLDTOT
       !===========================================================

       ! Column cloud fraction (not less than zero)
       CLDF1D = State_Met%CLDF(NLON,NLAT,1:State_Grid%NZ)
       WHERE ( CLDF1D < 0e+0_fp ) CLDF1D = 0e+0_fp

       ! NOTE: For GEOS-FP and MERRA-2 met fields, the optical
       ! depth is the in-cloud optical depth.  At this point it has
       ! NOT been multiplied by cloud fraction yet.  Therefore, we can
       ! just apply the ! we can just apply the linear overlap formula
       ! as written above (i.e. multiply by cloud fraction).
       OPTD = OPTD * CLDF1D

#elif defined( USE_APPROX_RANDOM_OVERLAP )
       !===========================================================
       ! %%%% CLOUD OVERLAP: APPROX RANDOM OVERLAP ASSUMPTION %%%%
       !
       ! Use OPTDEPTH = TAUCLD * CLDTOT**1.5
       !===========================================================

       ! Column cloud fraction (not less than zero)
       CLDF1D = State_Met%CLDF(NLON,NLAT,1:State_Grid%NZ)
       WHERE ( CLDF1D < 0e+0_fp ) CLDF1D = 0e+0_fp

       ! NOTE: For GEOS-FP and MERRA-2 met fields, the optical
       ! depth is the in-cloud optical depth.  At this point it has
       ! NOT been multiplied by cloud fraction yet.  Therefore, we can
       ! just apply the approximate random overlap formula as written
       ! above (i.e. multiply by cloud fraction^1.5).
       OPTD = OPTD * ( CLDF1D )**1.5e+0_fp

#elif defined( USE_MAXIMUM_RANDOM_OVERLAP )
       ! See commented out code in github history
       CALL ERROR_STOP('MMRAN_16 not yet FJX compatible.', &
                       'fjx_interface_mod.F90')


#endif

       ! Copy cloud OD data to a variable array (ewl: why???)
       DO L=1,L_
          ODCLOUD_COL(L) = OPTD(L)
       ENDDO
       
       ! Use GEOS-Chem methodology to set vertical profiles of:
       ! Pressure      (PPJ)    [hPa]
       ! Temperature   (T_CLIm) [K]
       ! Path density  (DDJ)    [# molec/cm2]
       ! New methodology for:
       ! Ozone density (OOJ)    [# O3 molec/cm2]
       CALL SET_PROF (YLAT,        MONTH,     DAY,         &
                      T_CTM,       P_CTM,     OPTD,        &
                      OPTDUST,     OPTAER,    O3_CTM,      &
                      O3_TOMS,     AERX_COL,  T_CLIM,      &
                      O3_CLIM,     Z_CLIM,    AIR_CLIM,    &
                      Input_Opt,   State_Grid )

       ! Call FAST-JX routines to compute J-values
       CALL PHOTO_JX( Input_Opt%amIRoot, Input_Opt%DryRun,          &
                      CSZA,      SFCA,       SZA,       SOLF,       &
                      P_CTM,     T_CTM,      AOD999,    NLON,       &
                      NLAT,      AERX_COL,   T_CLIM,    O3_CLIM,    &
                      Z_CLIM,    AIR_CLIM,   State_Grid%maxChemLev, &
                      VALJXX,     FSBOT,     FJBOT,     FLXD,       &
                      FJFLX                                        )

       ! Fill out common-block array of J-rates using PHOTO_JX output
       DO L=1,State_Grid%MaxChemLev
          DO J=1,NRATJ
             IF (JIND(J).gt.0) THEN
                ZPJ(L,J,NLON,NLAT) = VALJXX(L,JIND(J))*JFACTA(J)
             ELSE
                ZPJ(L,J,NLON,NLAT) = 0.e+0_fp
             ENDIF
          ENDDO
       ENDDO
       
       ! Set J-rates outside the chemgrid to zero
       IF (State_Grid%MaxChemLev.lt.L_) THEN
          DO L=State_Grid%MaxChemLev+1,L_
             DO J=1,NRATJ
                ZPJ(L,J,NLON,NLAT) = 0.e+0_fp
             ENDDO
          ENDDO
       ENDIF

       !=================================================================
       ! UV radiative fluxes (direct, diffuse, net) [W/m2]
       !
       ! Updated for netCDF from nd64 (JMM 2019-09-11)
       ! Use it to calculate fluxes for output if necessary
       !
       ! Get net direct and net diffuse fluxes separately
       ! Order:
       !    1 - Net flux
       !    2 - Direct flux
       !    3 - Diffuse flux
       ! Convention: negative is downwards
       !=================================================================
       IF ( State_Diag%Archive_UVFluxDiffuse .or. &
            State_Diag%Archive_UVFluxDirect .or. &
            State_Diag%Archive_UVFluxNet ) THEN
       
          ! Loop over wavelength bins
          DO K = 1, W_
       
             ! Initialize
             FDIRECT  = 0.0_fp
             FDIFFUSE = 0.0_fp
       
             ! Direct & diffuse fluxes at each level
             FDIRECT(1)  = FSBOT(K)                    ! surface
             FDIFFUSE(1) = FJBOT(K)                    ! surface
             DO L = 2, State_Grid%NZ
                FDIRECT(L) = FDIRECT(L-1) + FLXD(L-1,K)
                FDIFFUSE(L) = FJFLX(L-1,K)
             ENDDO
       
             ! Constant to multiply UV fluxes at each wavelength bin
             UVX_CONST = SOLF * FL(K) * UVXFACTOR(K)
       
             ! Archive into diagnostic arrays
             DO L = 1, State_Grid%NZ
       
                IF ( State_Diag%Archive_UVFluxNet ) THEN
                   S = State_Diag%Map_UvFluxNet%id2slot(K)
                   IF ( S > 0 ) THEN
                      State_Diag%UVFluxNet(NLON,NLAT,L,S) =  &
                      State_Diag%UVFluxNet(NLON,NLAT,L,S) +  &
                           ( ( FDIRECT(L) + FDIFFUSE(L) ) * UVX_CONST )
                   ENDIF
                ENDIF
       
                IF ( State_Diag%Archive_UVFluxDirect ) THEN
                   S = State_Diag%Map_UvFluxDirect%id2slot(K)
                   IF ( S > 0 ) THEN
                      State_Diag%UVFluxDirect(NLON,NLAT,L,S) =  &
                      State_Diag%UVFluxDirect(NLON,NLAT,L,S) +  &
                           ( FDIRECT(L) * UVX_CONST )
                   ENDIF
                ENDIF
       
                IF ( State_Diag%Archive_UVFluxDiffuse ) THEN
                   S = State_Diag%Map_UvFluxDiffuse%id2slot(K)
                   IF ( S > 0 ) THEN
                      State_Diag%UVFluxDiffuse(NLON,NLAT,L,S) =  &
                      State_Diag%UVFluxDiffuse(NLON,NLAT,L,S) +  &
                           ( FDIFFUSE(L) * UVX_CONST )
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Reset first-time flag
    FIRST=.FALSE.

  END SUBROUTINE FAST_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photrate_adj
!
! !DESCRIPTION: Subroutine PHOTRATE\_ADJ adjusts certain photolysis rates
!  for chemistry.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PHOTRATE_ADJ( Input_Opt, State_Chm,  State_Diag, State_Met,     &
                           I,         J,          L,                         &
                           FRAC,      RC                                    )
!
! !USES:
!
    USE CMN_Phot_Mod,   ONLY : ZPJ
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input_Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm  ! Chemistry State object
    TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
    INTEGER,        INTENT(IN)    :: I, J, L    ! Lon, lat, lev indices
    REAL(fp),       INTENT(IN)    :: FRAC       ! Result of SO4_PHOTFRAC,
                                                !  called from DO_FLEXCHEM
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
!
! !REMARKS:
!  NOTE: The netCDF diagnostics are attached in DO_FLEXCHEM so that we have
!  access to the adjusted rates.  Only the bpch diagnostics are updated
!  here.
!    -- Bob Yantosca, 19 Dec 2017
!
!  %%%% NOTE: WE SHOULD UPDATE THE COMMENTS TO MAKE SURE THAT WE DO      %%%%
!  %%%% NOT KEEP ANY CONFLICTING OR INCORRECT INFORMATION (bmy, 3/28/16) %%%%
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: RXN_JNITSa, RXN_JNITSb, RXN_JNITa, RXN_JNITb
    INTEGER  :: RXN_JHNO3,  RXN_H2SO4,  RXN_O3_1,  RXN_O3_2
    REAL(fp) :: JscaleNITs, JscaleNIT,  JNITChanA, JNITChanB
    REAL(fp) :: C_O2,     C_N2, C_H2,   ITEMPK, RO1DplH2O
    REAL(fp) :: RO1DplH2, RO1D, NUMDEN, TEMP,   C_H2O


    !=================================================================
    ! PHOTRATE_ADJ begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    TEMP    = State_Met%T(I,J,L)                                 ! K
    NUMDEN  = State_Met%AIRNUMDEN(I,J,L)                         ! molec/cm3
    C_H2O   = State_Met%AVGW(I,J,L) * State_Met%AIRNUMDEN(I,J,L) ! molec/cm3
    RXN_JNITSa = State_Chm%Photol%RXN_JNITSa
    RXN_JNITSb = State_Chm%Photol%RXN_JNITSb
    RXN_JNITa  = State_Chm%Photol%RXN_JNITa
    RXN_JNITb  = State_Chm%Photol%RXN_JNITb
    RXN_JHNO3  = State_Chm%Photol%RXN_JHNO3
    RXN_H2SO4  = State_Chm%Photol%RXN_H2SO4
    RXN_O3_1   = State_Chm%Photol%RXN_O3_1
    RXN_O3_2   = State_Chm%Photol%RXN_O3_2

    ! For all mechanisms. Set the photolysis rate of NITs and NIT to a
    ! scaled value of JHNO3. NOTE: this is set in geoschem_config.yml
    IF ( Input_Opt%hvAerNIT ) THEN

       ! Get the photolysis scalars read in from geoschem_config.yml
       JscaleNITs = Input_Opt%hvAerNIT_JNITs
       JscaleNIT  = Input_Opt%hvAerNIT_JNIT
       ! convert reaction channel % to a fraction
       JNITChanA  = Input_Opt%JNITChanA
       JNITChanB  = Input_Opt%JNITChanB
       JNITChanA  = JNITChanA / 100.0_fp
       JNITChanB  = JNITChanB / 100.0_fp
       ! Set the photolysis rate of NITs
       ZPJ(L,RXN_JNITSa,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNITs
       ZPJ(L,RXN_JNITSb,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNITs
       ! Set the photolysis rate of NIT
       ZPJ(L,RXN_JNITa,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNIT
       ZPJ(L,RXN_JNITb,I,J) = ZPJ(L,RXN_JHNO3,I,J) * JscaleNIT
       ! Adjust to scaling for channels set in geoschem_config.yml
       ! NOTE: channel scaling is 1 in FJX_j2j.dat, then updated here
       ZPJ(L,RXN_JNITSa,I,J) = ZPJ(L,RXN_JNITSa,I,J) * JNITChanA
       ZPJ(L,RXN_JNITa,I,J) = ZPJ(L,RXN_JNITa,I,J) * JNITChanA
       ZPJ(L,RXN_JNITSb,I,J) = ZPJ(L,RXN_JNITSb,I,J) * JNITChanB
       ZPJ(L,RXN_JNITb,I,J) = ZPJ(L,RXN_JNITb,I,J) * JNITChanB

    ! Gotcha to set JNIT and JNITs to zero if hvAerNIT switch is off
    ELSE

       ! Set the photolysis rate of NITs to zero
       ZPJ(L,RXN_JNITSa,I,J) = 0.0_fp
       ZPJ(L,RXN_JNITSb,I,J) = 0.0_fp
       ! Set the photolysis rate of NIT to zero
       ZPJ(L,RXN_JNITa,I,J) = 0.0_fp
       ZPJ(L,RXN_JNITb,I,J) = 0.0_fp

    ENDIF

    !==============================================================
    ! SPECIAL TREATMENT FOR H2SO4+hv -> SO2 + 2OH
    !
    ! Only allow photolysis of H2SO4 when gaseous (SDE 04/11/13)
    !==============================================================

    ! Calculate if H2SO4 expected to be gaseous or aqueous
    ! Only allow photolysis above 6 hPa
    ! RXN_H2SO4 specifies SO4 + hv -> SO2 + OH + OH
    ZPJ(L,RXN_H2SO4,I,J) = ZPJ(L,RXN_H2SO4,I,J) * FRAC

    !==============================================================
    ! SPECIAL TREATMENT FOR O3+hv -> O+O2
    !
    ! [O1D]ss=J[O3]/(k[H2O]+k[N2]+k[O2])
    ! SO, THE EFFECTIVE J-VALUE IS J*k[H2O]/(k[H2O]+k[N2]+k[O2])
    !
    ! We don't want to do this if strat-chem is in use, as all
    ! the intermediate reactions are included - this would be
    ! double-counting (SDE 04/01/13)
    !==============================================================

    ! Need to subtract O3->O1D from rate
    ! RXN_O3_1 specifies: O3 + hv -> O2 + O
    ! RXN_O3_2 specifies: O3 + hv -> O2 + O(1D)
    ZPJ(L,RXN_O3_1,I,J) = ZPJ(L,RXN_O3_1,I,J) &
                        - ZPJ(L,RXN_O3_2,I,J)

  END SUBROUTINE PHOTRATE_ADJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_aod
!
! !DESCRIPTION: Subroutine RD\_AOD reads aerosol phase functions that are
!  used to scale diagnostic output to an arbitrary wavelengh.  This
!  facilitates comparing with satellite observations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_AOD( NJ1, Input_Opt, RC )
!
! !USES:
!
    USE CMN_Phot_Mod,  ONLY : RDAA, RWAA, WVAA, NRAA, RHAA, SGAA, QQAA, REAA
    USE CMN_Phot_Mod,  ONLY : SSAA, ASYMAA, PHAA
    USE CMN_Phot_Mod,  ONLY : IWV1000, NSPAA, NWVAA, NRLAA, NCMAA, ALPHAA
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: NJ1         ! Unit # of file to open
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  The .dat files for each species contain the optical properties
!  at multiple wavelengths to be used in the online calculation of the aerosol
!  optical depth diagnostics.
!  These properties have been calculated using the same size and optical
!  properties as the FJX_spec.dat file used for the FAST-J photolysis
!  calculations (which is now redundant for aerosols, the values in the .dat
!  files here are now used). The file currently contains 11 wavelengths
!  for Fast-J and other commonly used wavelengths for satellite and
!  AERONET retrievals. 30 wavelengths follow that map onto RRTMG
!  wavebands for radiaitive flux calculations (not used if RRTMG is off).
!  A complete set of optical properties from 250-2000 nm for aerosols is
!  available at:
!  ftp://ftp.as.harvard.edu/geos-chem/data/aerosol_optics/hi_spectral_res
!                                                                             .
!     -- Colette L. Heald, 05/10/10)
!     -- David A. Ridley, 05/10/13 (update for new optics files)
!
! !REVISION HISTORY:
!  10 May 2010 - C. Heald      - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Scalars
    INTEGER            :: I, J, K, N
    INTEGER            :: IOS
    LOGICAL            :: LBRC, FileExists

    ! Strings
    CHARACTER(LEN=78 ) :: TITLE0
    CHARACTER(LEN=255) :: DATA_DIR
    CHARACTER(LEN=255) :: THISFILE
    CHARACTER(LEN=255) :: FileMsg
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! String arrays
    CHARACTER(LEN=30)  :: SPECFIL(8)

    !================================================================
    ! RD_AOD begins here!
    !================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at RD_AOD (in module GeosCore/fast_jx_mod.F90)'
    LBRC     = Input_Opt%LBRC
    DATA_DIR = TRIM( Input_Opt%FAST_JX_DIR )

    ! IMPORTANT: aerosol_mod.F and dust_mod.F expect aerosols in this order
    !
    ! Treating strat sulfate with GADS data but modified to match
    ! the old Fast-J values size (r=0.09um, sg=0.6) - I think there's
    ! evidence that this is too smale and narrow e.g. Deshler et al. 2003
    ! NAT should really be associated with something like cirrus cloud
    ! but for now we are just treating the NAT like the sulfate... limited
    ! info but ref index is similar e.g. Scarchilli et al. (2005)
    !(DAR 05/2015)
    DATA SPECFIL /"so4.dat","soot.dat","org.dat", &
                  "ssa.dat","ssc.dat",            &
                  "h2so4.dat","h2so4.dat",        &
                  "dust.dat"/

    ! Loop over the array of filenames
    DO k = 1, NSPAA

       ! Choose different set of input files for standard (trop+strat chenm)
       ! and tropchem (trop-only chem) simulations
       THISFILE = TRIM( DATA_DIR ) // TRIM( SPECFIL(k) )

       !--------------------------------------------------------------
       ! In dry-run mode, print file path to dryrun log and cycle.
       ! Otherwise, print file path to stdout and continue.
       !--------------------------------------------------------------

       ! Test if the file exists
       INQUIRE( FILE=TRIM( ThisFile ), EXIST=FileExists )

       ! Test if the file exists and define an output string
       IF ( FileExists ) THEN
          FileMsg = 'FAST-JX (RD_AOD): Opening'
       ELSE
          FileMsg = 'FAST-JX (RD_AOD): REQUIRED FILE NOT FOUND'
       ENDIF

       ! Write to stdout for both regular and dry-run simulations
       IF ( Input_Opt%amIRoot ) THEN
          WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
300       FORMAT( a, ' ', a )
       ENDIF

       ! For dry-run simulations, cycle to next file.
       ! For regular simulations, throw an error if we can't find the file.
       IF ( Input_Opt%DryRun ) THEN
          CYCLE
       ELSE
          IF ( .not. FileExists ) THEN
             WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       !--------------------------------------------------------------
       ! If not a dry-run, read data from each species file
       !--------------------------------------------------------------

#if defined( MODEL_CESM )
       ! Only read file on root thread if using CESM
       IF ( Input_Opt%amIRoot ) THEN
#endif

       ! Open file
       OPEN( NJ1, FILE=TRIM( THISFILE ), STATUS='OLD', IOSTAT=RC )

       ! Error check
       IF ( RC /= 0 ) THEN
          ErrMsg = 'Error opening file: ' // TRIM( ThisFile )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Read header lines
       READ(  NJ1, '(A)' ) TITLE0
       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0

       ! Second header line added for more info
       READ(  NJ1, '(A)' ) TITLE0
       IF ( Input_Opt%amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0

       READ(  NJ1, '(A)' ) TITLE0
110    FORMAT( 3x, a20 )

       DO i = 1, NRAA
       DO j = 1, NWVAA

          READ(NJ1,*) WVAA(j,k),RHAA(i,k),NRLAA(j,i,k),NCMAA(j,i,k), &
                      RDAA(i,k),RWAA(i,k),SGAA(i,k),QQAA(j,i,k),   &
                      ALPHAA(j,i,k),REAA(i,k),SSAA(j,i,k),         &
                      ASYMAA(j,i,k),(PHAA(j,i,k,n),n=1,8)

          ! make note of where 1000nm is for FAST-J calcs
          IF (WVAA(j,k).EQ.1000.0) IWV1000=J

       ENDDO
       ENDDO

       ! Close file
       CLOSE( NJ1 )

#if defined( MODEL_CESM )
       ENDIF
#endif

    ENDDO

#if defined( MODEL_CESM ) && defined( SPMD )
    CALL MPIBCAST( WVAA,      Size(WVAA),     MPIR8,   0, MPICOM )
    CALL MPIBCAST( RHAA,      Size(RHAA),     MPIR8,   0, MPICOM )
    CALL MPIBCAST( NRLAA,     Size(NRLAA),    MPIR8,   0, MPICOM )
    CALL MPIBCAST( NCMAA,     Size(NCMAA),    MPIR8,   0, MPICOM )
    CALL MPIBCAST( RDAA,      Size(RDAA),     MPIR8,   0, MPICOM )
    CALL MPIBCAST( RWAA,      Size(RWAA),     MPIR8,   0, MPICOM )
    CALL MPIBCAST( SGAA,      Size(SGAA),     MPIR8,   0, MPICOM )
    CALL MPIBCAST( QQAA,      Size(QQAA),     MPIR8,   0, MPICOM )
    CALL MPIBCAST( ALPHAA,    Size(ALPHAA),   MPIR8,   0, MPICOM )
    CALL MPIBCAST( REAA,      Size(REAA),     MPIR8,   0, MPICOM )
    CALL MPIBCAST( SSAA,      Size(SSAA),     MPIR8,   0, MPICOM )
    CALL MPIBCAST( ASYMAA,    Size(ASYMAA),   MPIR8,   0, MPICOM )
    CALL MPIBCAST( PHAA,      Size(PHAA),     MPIR8,   0, MPICOM )
    CALL MPIBCAST( IWV1000,   1,              MPIINT,  0, MPICOM )
#endif

! ewl: this and subroutine calc_aod moved to fast_jx_interface_mod.F90
!    !=================================================================
!    ! Only do the following if we are not running in dry-run mode
!    !=================================================================
!    IF ( .not. Input_Opt%DryRun ) THEN
!
!       IF ( Input_Opt%amIRoot ) THEN
!          WRITE( 6, * ) 'Optics read for all wavelengths successfully'
!       ENDIF
!
!       ! Now calculate the required wavelengths in the LUT to calculate
!       ! the requested AOD
!       CALL CALC_AOD( Input_Opt )
!    ENDIF

  END SUBROUTINE RD_AOD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_aod
!
! !DESCRIPTION: Subroutine CALC\_AOD works out the closest tie points
! in the optics LUT wavelengths and the coefficients required to
! calculate the angstrom exponent for interpolating optics to the requested
! wavelength. If the wavelength requested matches a standard wavelength
! in the LUT then we skip the interpolation (DAR 09/2013)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_AOD( Input_Opt )
!
! !USES:
!
    USE CMN_Phot_Mod,  ONLY : NWVAA, NWVAA0, WVAA
    USE CMN_Phot_Mod,  ONLY : IWVSELECT
    USE CMN_Phot_Mod,  ONLY : IRTWVSELECT
    USE CMN_Phot_Mod,  ONLY : ACOEF_WV, BCOEF_WV, CCOEF_WV
    USE CMN_Phot_Mod,  ONLY : ACOEF_RTWV, BCOEF_RTWV, CCOEF_RTWV
    USE CMN_Phot_Mod,  ONLY : NWVREQUIRED, IWVREQUIRED
    USE CMN_Phot_Mod,  ONLY : NRTWVREQUIRED, IRTWVREQUIRED
    USE Input_Opt_Mod, ONLY : OptInput
#ifdef RRTMG
    USE PARRRTM,     ONLY : NBNDLW
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN) :: Input_Opt
!
! !REMARKS:
!  Now the user is able to select any 3 wavelengths for optics
!  output in the geoschem_config.yml file we need to be able to interpolate
!  to those wavelengths based on what is available in the optics
!  look-up table.
!                                                                             .
!  The standard lookup table currently has values for
!  11 common wavelengths followed by 30 that are required by RRTMG.
!  Only those required to interpolate to user requested
!  wavelengths are selected from the standard wavelengths. RRTMG
!  wavelengths are not used in the interpolation for AOD output
!  (DAR 10/2013)
!                                                                             .
!   UPDATE: because the RT optics output doesnt have access to the
!   standard wavelengths we now calculate two sets of values: one
!   for the ND21 and diag3 outputs that use the standard wavelengths
!   and one for RRTMG diagnostics that interpolate the optics from RRTMG
!   wavelengths. Perhaps a switch needs adding to switch off the RT
!   optics output (and interpolation) if this ends up costing too
!   much and is not used, but it is ideal to have an optics output
!   that matches exactly what RRTMG uses to calculate the fluxes
!
! !REVISION HISTORY:
!  18 Jun 2013 - D. Ridley   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER             :: MINWV, MAXWV, N, N0, N1, W, NSTEP
    REAL(fp)            :: WVDIF

    !================================================================
    ! CALC_AOD begins here!
    !================================================================

    !cycle over standard wavelengths
    N0=1
    N1=NWVAA0
    NSTEP=1
    NWVREQUIRED=0
    DO W=1,Input_Opt%NWVSELECT
       MINWV     = -999
       MAXWV     =  999
       DO N=N0,N1,NSTEP
          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
             MINWV = WVDIF
             IWVSELECT(1,W)=N
          ENDIF
          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
             MAXWV = WVDIF
             IWVSELECT(2,W)=N
          ENDIF
       ENDDO
       IF (IWVSELECT(2,W).EQ.IWVSELECT(1,W)) THEN
          !we have a match!
          MINWV=0
          MAXWV=0
          !add this wavelength to those for output
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ELSE
          !we are going to have to interpolate to the requested wavelength
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
          NWVREQUIRED=NWVREQUIRED+1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(2,W)
       ENDIF

       !Error check - ensure we have a match or requested wavelength
       !falls within two LUT tie points
       IF (MINWV.EQ.-999) THEN
          ! requested wavelength is shorter than min wv in LUT
          ! set to min
          write(6,*) 'ERROR requested wavelength is too short!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(1,1)
          IWVSELECT(1,W)=1
          IWVSELECT(2,W)=1 !300nm
          NWVREQUIRED=NWVREQUIRED-1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ENDIF
       IF (MAXWV.EQ.999) THEN
          ! requested wavelength is longer than min wv in LUT
          ! set to max
          write(6,*) 'ERROR requested wavelength is too long!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0,1)
          IWVSELECT(1,W)=NWVAA0
          IWVSELECT(2,W)=NWVAA0 !1020nm
          NWVREQUIRED=NWVREQUIRED-1
          IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
       ENDIF

       !now calcualte the angstrom exponent coefs for interpolation -
       !this is done here to save time and repetition in aerosol_mod.F
       IF (IWVSELECT(1,W).NE.IWVSELECT(2,W)) THEN
          ACOEF_WV(W) = WVAA(IWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
          BCOEF_WV(W) =1.0d0/(LOG(WVAA(IWVSELECT(2,W),1)/ &
                                  WVAA(IWVSELECT(1,W),1)))
          !relative location of selected wavelength between tie points
          !for interpolating SSA and ASYM for output in aerosol_mod.F and
          !dust_mod.F
          CCOEF_WV(W) =(Input_Opt%WVSELECT(W)-WVAA(IWVSELECT(1,W),1))/ &
                      (WVAA(IWVSELECT(2,W),1)-WVAA(IWVSELECT(1,W),1))
       ENDIF
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) 'N WAVELENGTHS: ',Input_Opt%NWVSELECT
          write(6,*) 'WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
          write(6,*) 'WAVELENGTH REQUIRED:', NWVREQUIRED
          !write(6,*) IWVSELECT(1,W),WVAA(IWVSELECT(1,W),1)
          !write(6,*) IWVSELECT(2,W),WVAA(IWVSELECT(2,W),1)
          !write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
          write(6,*) '*********************************'
       ENDIF
    ENDDO !Input_Opt%NWVSELECT
#ifdef RRTMG
    !repeat for RRTMG wavelengths to get the closest wavelength
    !indices and the interpolation coefficients
    !Indices are relative to all wavelengths in the LUT i.e. the RRTMG
    !wavelengths start at NWVAA0+1
    N0=NWVAA0+1
    N1=NWVAA
    NSTEP=1
    NRTWVREQUIRED=0
    DO W=1,Input_Opt%NWVSELECT
       MINWV     = -999
       MAXWV     =  999
       DO N=N0,N1,NSTEP
          WVDIF = WVAA(N,1)-Input_Opt%WVSELECT(W)
          IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
             MINWV = WVDIF
             IRTWVSELECT(1,W)=N
          ENDIF
          IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
             MAXWV = WVDIF
             IRTWVSELECT(2,W)=N
          ENDIF
       ENDDO
       IF (IRTWVSELECT(2,W).EQ.IRTWVSELECT(1,W)) THEN
          !we have a match!
          MINWV=0
          MAXWV=0
          !add this wavelength to those for output
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ELSE
          !we are going to have to interpolate to the requested
          !wavelength
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
          NRTWVREQUIRED=NRTWVREQUIRED+1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(2,W)
       ENDIF

       !Error check - ensure we have a match or requested wavelength
       !falls within two LUT tie points
       IF (MINWV.EQ.-999) THEN
          ! requested wavelength is shorter than min wv in LUT
          ! set to min
          write(6,*) 'ERROR requested wavelength is too short!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA-1,1)
          IRTWVSELECT(1,W)=NWVAA-1
          IRTWVSELECT(2,W)=NWVAA-1
          NRTWVREQUIRED=NRTWVREQUIRED-1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ENDIF
       IF (MAXWV.EQ.999) THEN
          ! requested wavelength is longer than min wv in LUT
          ! set to max
          write(6,*) 'ERROR requested wavelength is too long!!'
          write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA0+1,1)
          IRTWVSELECT(1,W)=NWVAA0+1
          IRTWVSELECT(2,W)=NWVAA0+1
          NRTWVREQUIRED=NRTWVREQUIRED-1
          IRTWVREQUIRED(NRTWVREQUIRED)=IRTWVSELECT(1,W)
       ENDIF

       !now calcualte the angstrom exponent coefs for interpolation -
       !this is done here to save time and repetition in aerosol_mod.F
       IF (IRTWVSELECT(1,W).NE.IRTWVSELECT(2,W)) THEN
          ACOEF_RTWV(W) = WVAA(IRTWVSELECT(2,W),1)/Input_Opt%WVSELECT(W)
          BCOEF_RTWV(W) =1.0d0/(LOG(WVAA(IRTWVSELECT(2,W),1)/ &
                                    WVAA(IRTWVSELECT(1,W),1)))
          !relative location of selected wavelength between tie points
          !for interpolating SSA and ASYM for output in aerosol_mod.F and
          !dust_mod.F
          CCOEF_RTWV(W) =(Input_Opt%WVSELECT(W)-WVAA(IRTWVSELECT(1,W),1))/ &
                      (WVAA(IRTWVSELECT(2,W),1)-WVAA(IRTWVSELECT(1,W),1))
       ENDIF
       !convert wavelength index to that required by rrtmg_rad_transfer
       !i.e. without the standard and LW wavelengths
       IRTWVSELECT(1,W) = IRTWVSELECT(1,W) - NWVAA0 - NBNDLW
       IRTWVSELECT(2,W) = IRTWVSELECT(2,W) - NWVAA0 - NBNDLW
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) 'N RT WAVELENGTHS: ',Input_Opt%NWVSELECT
          write(6,*) 'RT WAVELENGTH REQUESTED:',Input_Opt%WVSELECT(W)
          write(6,*) 'RT WAVELENGTH REQUIRED:', NRTWVREQUIRED
          write(6,*) IRTWVSELECT(1,W),WVAA(IRTWVSELECT(1,W)+NWVAA0+NBNDLW,1)
          write(6,*) IRTWVSELECT(2,W),WVAA(IRTWVSELECT(2,W)+NWVAA0+NBNDLW,1)
          write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
          write(6,*) '*********************************'
       ENDIF
    ENDDO !Input_Opt%NWVSELECT
#endif
  END SUBROUTINE CALC_AOD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_prof_nc
!
! !DESCRIPTION: Subroutine RAD\_PROF\_NC reads in the reference climatology
!  from a NetCDF file rather than an ASCII .dat.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_PROF_NC( Input_Opt, RC )
!
! !USES:
!
    USE CMN_Phot_Mod,  ONLY : OREF, TREF
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput

#if defined( MODEL_CESM )
    USE CAM_PIO_UTILS,     ONLY : CAM_PIO_OPENFILE
    USE IOFILEMOD,         ONLY : GETFIL
    USE PIO,               ONLY : PIO_CLOSEFILE
    USE PIO,               ONLY : PIO_INQ_DIMID
    USE PIO,               ONLY : PIO_INQ_DIMLEN
    USE PIO,               ONLY : PIO_INQ_VARID
    USE PIO,               ONLY : PIO_GET_VAR
    USE PIO,               ONLY : PIO_NOERR
    USE PIO,               ONLY : PIO_NOWRITE
    USE PIO,               ONLY : FILE_DESC_T
#else
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This file was automatically generated by the Perl scripts in the
!  NcdfUtilities package (which ships w/ GEOS-Chem) and was subsequently
!  hand-edited.
!
! !REVISION HISTORY:
!  19 Apr 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists          ! Does input file exist?
    INTEGER            :: fId                 ! netCDF file ID

    ! Strings
    CHARACTER(LEN=255) :: nc_dir              ! netCDF directory name
    CHARACTER(LEN=255) :: nc_file             ! netCDF file name
    CHARACTER(LEN=255) :: nc_path             ! netCDF path name
    CHARACTER(LEN=255) :: v_name              ! netCDF variable name
    CHARACTER(LEN=255) :: a_name              ! netCDF attribute name
    CHARACTER(LEN=255) :: a_val               ! netCDF attribute value
    CHARACTER(LEN=255) :: FileMsg
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! Arrays
    INTEGER            :: st3d(3), ct3d(3)    ! For 3D arrays

#if defined( MODEL_CESM )
    type(FILE_DESC_T)  :: ncid
    INTEGER            :: vId, iret
#endif

    !=================================================================
    ! RD_PROF_NC begins here!
    !=================================================================

    ! Initialize
    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at RD_PROF_NC (in module GeosCore/fjx_interface_mod.F90)'

    ! Directory and file names
    nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // 'FastJ_201204/'
    nc_file = 'fastj.jv_atms_dat.nc'
    nc_path = TRIM( nc_dir ) // TRIM( nc_file )

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Test if the file exists
    INQUIRE( FILE=TRIM( nc_path ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'FAST-JX (RD_PROF_NC): Opening'
    ELSE
       FileMsg = 'FAST-JX (RD_PROF_NC): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( nc_path )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( Input_Opt%DryRun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( nc_path )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=========================================================================
    ! Open and read data from the netCDF file
    !=========================================================================

    ! Open netCDF file
#if defined( MODEL_CESM )
    CALL CAM_PIO_OPENFILE( ncid, TRIM(nc_path), PIO_NOWRITE )
#else
    CALL Ncop_Rd( fId, TRIM(nc_path) )
#endif

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) REPEAT( '%', 79 )
       WRITE( 6, 110 ) TRIM(nc_file)
       WRITE( 6, 120 ) TRIM(nc_dir)
    ENDIF

    !----------------------------------------
    ! VARIABLE: T
    !----------------------------------------

    ! Variable name
    v_name = "T"

    ! Read T from file
    st3d   = (/  1,  1,  1 /)
    ct3d   = (/ 51, 18, 12 /)
#if defined( MODEL_CESM )
    iret = PIO_INQ_VARID( ncid, trim(v_name), vid  )
    iret = PIO_GET_VAR(   ncid, vid, TREF          )
#else
    CALL NcRd( TREF, fId, TRIM(v_name), st3d, ct3d )

    ! Read the T:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF
#endif

    !----------------------------------------
    ! VARIABLE: O3
    !----------------------------------------

    ! Variable name
    v_name = "O3"

    ! Read O3 from file
    st3d   = (/  1,  1,  1 /)
    ct3d   = (/ 51, 18, 12 /)
#if defined( MODEL_CESM )
    iret = PIO_INQ_VARID( ncid, trim(v_name), vid  )
    iret = PIO_GET_VAR(   ncid, vid, OREF          )
#else
    CALL NcRd( OREF, fId, TRIM(v_name), st3d, ct3d )

    ! Read the O3:units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
    ENDIF
#endif

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Close netCDF file
#if defined( MODEL_CESM )
    CALL PIO_CLOSEFILE( ncid )
#else
    CALL NcCl( fId )
#endif

    ! Echo info to stdout
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 140 )
       WRITE( 6, 100 ) REPEAT( '%', 79 )
    ENDIF

    ! FORMAT statements
100 FORMAT( a                                              )
110 FORMAT( '%% Opening file  : ',         a               )
120 FORMAT( '%%  in directory : ',         a, / , '%%'     )
130 FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
140 FORMAT( '%% Successfully closed file!'                 )

  END SUBROUTINE RD_PROF_NC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_prof
!
! !DESCRIPTION: Subroutine SET\_PROF sets vertical profiles for a given
!  latitude and longitude.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_PROF( YLAT,      MONTH,  DAY,     T_CTM,  P_CTM,    &
                       CLDOD,     DSTOD,  AEROD,   O3_CTM, O3_TOMS,  &
                       AERCOL,    T_CLIM, O3_CLIM, Z_CLIM, AIR_CLIM, &
                       Input_Opt, State_Grid )
!
! !USES:
!
    USE CMN_FastJX_Mod,     ONLY : L_, L1_, A_, ZZHT
    USE CMN_Phot_Mod,       ONLY : OREF, TREF
    USE CMN_SIZE_Mod,       ONLY : NAER, NRH, NDUST
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AIRMW, AVO, g0, BOLTZ
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)       :: YLAT              ! Latitude (degrees)
    INTEGER,  INTENT(IN)       :: MONTH             ! Month
    INTEGER,  INTENT(IN)       :: DAY               ! Day *of month*
    REAL(fp), INTENT(IN)       :: T_CTM(L1_)        ! CTM temperatures (K)
    REAL(fp), INTENT(IN)       :: O3_TOMS           ! O3 column (DU)
    REAL(fp), INTENT(IN)       :: P_CTM(L1_+1)      ! CTM edge pressures (hPa)
    REAL(fp), INTENT(INOUT)    :: CLDOD(L_)         ! Cloud optical depth
    REAL(fp), INTENT(IN)       :: DSTOD(L_,NDUST)   ! Mineral dust OD
    REAL(fp), INTENT(IN)       :: AEROD(L_,A_)      ! Aerosol OD
    REAL(fp), INTENT(IN)       :: O3_CTM(L1_)       ! CTM ozone (molec/cm3)
    TYPE(OptInput), INTENT(IN) :: Input_Opt         ! Input options
    TYPE(GrdState), INTENT(IN) :: State_Grid        ! Grid State object
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT)      :: AERCOL(A_,L1_)    ! Aerosol column
    REAL(fp), INTENT(OUT)      :: T_CLIM(L1_)       ! Clim. temperatures (K)
    REAL(fp), INTENT(OUT)      :: Z_CLIM(L1_+1)     ! Edge altitudes (cm)
    REAL(fp), INTENT(OUT)      :: O3_CLIM(L1_)      ! O3 column depth (#/cm2)
    REAL(fp), INTENT(OUT)      :: AIR_CLIM(L1_)     ! O3 column depth (#/cm2)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  30 Mar 2013 - S. D. Eastham - Adapted from J. Mao code
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                  :: I, K, L, M, N, LCTM
    REAL(fp)                 :: DLOGP,F0,T0,B0,PB,PC,XC,MASFAC,SCALEH
    REAL(fp)                 :: PSTD(52),OREF2(51),TREF2(51)
    REAL(fp)                 :: PROFCOL, ODSUM
    REAL(fp), PARAMETER      :: ODMAX = 200.0e+0_fp

    ! Local variables for quantities from Input_Opt
    LOGICAL :: USE_ONLINE_O3

    !=================================================================
    ! SET_PROF begins here!
    !=================================================================

    ! Copy fields from INPUT_OPT
    USE_ONLINE_O3   = Input_Opt%USE_ONLINE_O3

    ! Zero aerosol column
    DO K=1,A_
       DO I=1,L1_
          AERCOL(K,I) = 0.e+0_fp
       ENDDO
    ENDDO

    ! Scale optical depths to stay within limits
    ODSUM = 0.e+0_fp
    DO I=1,L_
       CLDOD(I) = DBLE(CLDOD(I))
       ODSUM = ODSUM + CLDOD(I)
    ENDDO
    IF (ODSUM.gt.ODMAX) THEN
       ODSUM = ODMAX/ODSUM ! Temporary
       DO I=1,L_
          CLDOD(I) = CLDOD(I)*ODSUM
       ENDDO
       ODSUM = ODMAX
    ENDIF

    !=================================================================
    ! Set up pressure levels for O3/T climatology - assume that value
    ! given for each 2 km z* level applies from 1 km below to 1 km
    ! above, so select pressures at these boundaries. Surface level
    ! values at 1000 mb are assumed to extend down to the actual
    ! surface pressure for this lat/lon.
    !=================================================================
    PSTD(1)  = MAX(P_CTM(1),1000.e+0_fp)
    PSTD(2)  = 1000.e+0_fp * 10.e+0_fp ** (-1.e+0_fp/16.e+0_fp)
    DLOGP    = 10.e+0_fp**(-2.e+0_fp/16.e+0_fp)
    DO I=3,51
       PSTD(I) = PSTD(I-1) * DLOGP
    ENDDO
    PSTD(52) = 0.e+0_fp

    ! Mass factor - delta-Pressure [hPa] to delta-Column [molec/cm2]
    MASFAC = 100.e+0_fp * AVO / ( AIRMW * g0 * 10.e+0_fp )

    ! Select appropriate monthly and latitudinal profiles
    ! Now use YLAT instead of Oliver's YDGRD(NSLAT) (bmy, 9/13/99)
    M = MAX( 1, MIN( 12, MONTH                   ) )
    L = MAX( 1, MIN( 18, ( INT(YLAT) + 99 ) / 10 ) )

    ! Temporary arrays for climatology data
    DO I = 1, 51
       OREF2(I) = OREF(I,L,M)
       TREF2(I) = TREF(I,L,M)
    ENDDO

    ! Apportion O3 and T on supplied climatology z* levels onto CTM levels
    ! with mass (pressure) weighting, assuming constant mixing ratio and
    ! temperature half a layer on either side of the point supplied.
    DO I = 1, L1_
       F0 = 0.e+0_fp
       T0 = 0.e+0_fp
       DO K = 1, 51
          PC = MIN( P_CTM(I),   PSTD(K)   )
          PB = MAX( P_CTM(I+1), PSTD(K+1) )
          IF ( PC .GT. PB ) THEN
             XC = ( PC - PB ) / ( P_CTM(I) - P_CTM(I+1) )
             F0 = F0 + OREF2(K)*XC
             T0 = T0 + TREF2(K)*XC
          ENDIF
       ENDDO
       T_CLIM(I)  = T0
       O3_CLIM(I) = F0 * 1.e-6_fp
    ENDDO

    !=================================================================
    ! Calculate effective altitudes using scale height at each level
    !=================================================================
    Z_CLIM(1) = 0.e+0_fp
    DO I = 1, L_
       SCALEH = BOLTZ * 1.e+4_fp * MASFAC * T_CLIM(I)

       Z_CLIM(I+1) = Z_CLIM(I) - ( LOG( P_CTM(I+1) / P_CTM(I) ) * SCALEH )
    ENDDO
    Z_CLIM(L1_+1)=Z_CLIM(L1_) + ZZHT

    !=================================================================
    ! Add Aerosol Column - include aerosol types here. Currently use
    ! soot water and ice; assume black carbon x-section of 10 m2/g,
    ! independent of wavelength; assume limiting temperature for
    ! ice of -40 deg C.
    !=================================================================
    DO I = 1, L_
       ! Turn off uniform black carbon profile (rvm, bmy, 2/27/02)
       AERCOL(1,I) = 0e+0_fp

       IF ( T_CTM(I) .GT. 233.e+0_fp ) THEN
          AERCOL(2,I) = CLDOD(I)
          AERCOL(3,I) = 0.e+0_fp
       ELSE
          AERCOL(2,I) = 0.e+0_fp
          AERCOL(3,I) = CLDOD(I)
       ENDIF

       ! Also add in aerosol optical depth columns (rvm, bmy, 9/30/00)
       DO N = 1, NDUST
          AERCOL(3+N,I) = DSTOD(I,N)
       ENDDO

       ! Also add in other aerosol optical depth columns (rvm, bmy, 2/27/02)
       DO N = 1, NAER*NRH
          AERCOL(3+N+NDUST,I) = AEROD(I,N)
       ENDDO

    ENDDO

    DO K = 1,(3+NDUST+(NAER))
       AERCOL(K,L1_    ) = 0.e+0_fp
    ENDDO

    !=================================================================
    ! Calculate column quantities for FAST-JX
    !=================================================================
    PROFCOL = 0e+0_fp

    DO I = 1, L1_

       ! Monthly mean air Column [molec/cm2]
       AIR_CLIM(I)  = ( P_CTM(I) - P_CTM(I+1) ) * MASFAC

       ! Monthly mean O3 column [molec/cm2]
       O3_CLIM(I) = O3_CLIM(I) * AIR_CLIM(I)

       ! Monthly mean O3 column [DU]
       PROFCOL = PROFCOL + ( O3_CLIM(I) / 2.69e+16_fp )
    ENDDO

    !! Top values are special (do not exist in CTM data)
    !AIR_CLIM(L1_)     = P_CTM(L1_) * MASFAC
    !O3_CLIM(L1_) = O3_CLIM(L1_) * AIR_CLIM(L1_)

    !=================================================================
    ! Now weight the O3 column by the observed monthly mean TOMS.
    ! Missing data is denoted by the flag -999. (mje, bmy, 7/15/03)
    !
    ! TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
    ! Resolution:  5 x 10 deg.
    !
    ! Methodology (bmy, 2/12/07)
    ! ----------------------------------------------------------------
    ! FAST-J comes with its own default O3 column climatology (from
    ! McPeters 1992 & Nagatani 1991), which is stored in the input
    ! file "jv_atms.dat".  These "FAST-J default" O3 columns are used
    ! in the computation of the actinic flux and other optical
    ! quantities for the FAST-J photolysis.
    !
    ! The TOMS/SBUV O3 columns and 1/2-monthly O3 trends (contained
    ! in the TOMS_200701 directory) are read into GEOS-Chem by routine
    ! READ_TOMS in "toms_mod.f".  Missing values (i.e. locations where
    ! there are no data) in the TOMS/SBUV O3 columns are defined by
    ! the flag -999.
    !
    ! After being read from disk in routine READ_TOMS, the TOMS/SBUV
    ! O3 data are then passed to the FAST-J routine "set_prof.f".  In
    ! "set_prof.f", a test is done to make sure that the TOMS/SBUV O3
    ! columns and 1/2-monthly trends do not have any missing values
    ! for (lat,lon) location for the given month.  If so, then the
    ! TOMS/SBUV O3 column data is interpolated to the current day and
    ! is used to weight the "FAST-J default" O3 column.  This
    ! essentially "forces" the "FAST-J default" O3 column values to
    ! better match the observations, as defined by TOMS/SBUV.
    !
    ! If there are no TOMS/SBUV O3 columns (and 1/2-monthly trends)
    ! at a (lat,lon) location for given month, then FAST-J will revert
    ! to its own "default" climatology for that location and month.
    ! Therefore, the TOMS O3 can be thought of as an  "overlay" data
    ! -- it is only used if it exists.
    !
    ! Note that there are no TOMS/SBUV O3 columns at the higher
    ! latitudes.  At these latitudes, the code will revert to using
    ! the "FAST-J default" O3 columns.
    !
    ! As of February 2007, we have TOMS/SBUV data for 1979 thru 2005.
    ! 2006 TOMS/SBUV data is incomplete as of this writing.  For years
    ! 2006 and onward, we use 2005 TOMS O3 columns.
    !
    ! This methodology was originally adopted by Mat Evans.  Symeon
    ! Koumoutsaris was responsible for creating the downloading and
    ! processing the TOMS O3 data files from 1979 thru 2005 in the
    ! TOMS_200701 directory.
    !=================================================================

    ! Since we now have stratospheric ozone calculated online, use
    ! this instead of archived profiles for all chemistry-grid cells
    ! The variable O3_CTM is obtained from State_Met%Species, and will be 0
    ! outside the chemgrid (in which case we use climatology)

    ! Scale monthly O3 profile to the daily O3 profile (if available)
    DO I = 1, L1_

       ! Use online O3 values in the chemistry grid if selected
       IF ( (USE_ONLINE_O3) .and. &
            (I <= State_Grid%MaxChemLev) .and. &
            (O3_CTM(I) > 0e+0_fp) ) THEN

          ! Convert from molec/cm3 to molec/cm2
          O3_CLIM(I) = O3_CTM(I) * (Z_CLIM(I+1)-Z_CLIM(I))

       ! Otherwise, use O3 values from the met fields or TOMS/SBUV
       ELSEIF (O3_TOMS > 0e+0_fp) THEN

          O3_CLIM(I) = O3_CLIM(I) * ( O3_TOMS / PROFCOL )

       ENDIF

    ENDDO

  END SUBROUTINE SET_PROF
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_aer
!
! !DESCRIPTION: Subroutine SET\_AER fills out the array MIEDX.
!  Each entry connects a GEOS-Chem aerosol to its Fast-JX counterpart:
!  MIEDX(Fast-JX index) = (GC index)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_AER( Input_Opt )
!
! !USES:
!
    USE CMN_FastJX_Mod, ONLY : AN_, NAA, TITLAA
    USE CMN_Phot_Mod,   ONLY : MIEDX
    USE CMN_SIZE_Mod,   ONLY : NRHAER, NRH
    USE Input_Opt_Mod,  ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt ! Input options
!
! !REVISION HISTORY:
!  31 Mar 2013 - S. D. Eastham - Adapted from J. Mao FJX v6.2 implementation
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I, J, K
    INTEGER               :: IND(NRHAER)

    !=================================================================
    ! SER_AER begins here!
    !=================================================================

    ! Taken from aerosol_mod.F
    IND = (/22,29,36,43,50/)

    DO I=1,AN_
       MIEDX(I) = 0
    ENDDO

    ! Select Aerosol/Cloud types to be used - define types here
    ! Each of these types must be listed in the order used by OPMIE.F

    ! Clouds
    MIEDX(1)  =  3   !  Black carbon absorber
    MIEDX(2)  = 10   !  Water Cloud (Deirmenjian 8 micron)
    MIEDX(3)  = 14   !  Irregular Ice Cloud (Mishchenko)

    ! Dust
    MIEDX(4)  = 15   !  Mineral Dust  .15 micron    (rvm, 9/30/00)
    MIEDX(5)  = 16   !  Mineral Dust  .25 micron    (rvm, 9/30/00)
    MIEDX(6)  = 17   !  Mineral Dust  .4  micron    (rvm, 9/30/00)
    MIEDX(7)  = 18   !  Mineral Dust  .8  micron    (rvm, 9/30/00)
    MIEDX(8)  = 19   !  Mineral Dust 1.5  micron    (rvm, 9/30/00)
    MIEDX(9)  = 20   !  Mineral Dust 2.5  micron    (rvm, 9/30/00)
    MIEDX(10) = 21   !  Mineral Dust 4.0  micron    (rvm, 9/30/00)

    ! Aerosols
    DO I=1,NRHAER
       DO J=1,NRH
          MIEDX(10+((I-1)*NRH)+J)=IND(I)+J-1
       ENDDO
    ENDDO

    ! Stratospheric aerosols - SSA/STS and solid PSCs
    MIEDX(10+(NRHAER*NRH)+1) = 4  ! SSA/LBS/STS
    MIEDX(10+(NRHAER*NRH)+2) = 14 ! NAT/ice PSCs

    ! Ensure all 'AN_' types are valid selections
    do i=1,AN_
       IF (Input_Opt%amIRoot) write(6,1000) MIEDX(i),TITLAA(MIEDX(i))
       if (MIEDX(i).gt.NAA.or.MIEDX(i).le.0) then
          if (Input_Opt%amIRoot) then
             write(6,1200) MIEDX(i),NAA
          endif
          CALL GC_EXITC('Bad MIEDX value.')
       endif
    enddo

1000 format('Using Aerosol type: ',i3,1x,a)
1200 format('Aerosol type ',i3,' unsuitable; supplied values must be ', &
            'between 1 and ',i3)

  END SUBROUTINE SET_AER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_exitc
!
! !DESCRIPTION: Subroutine GC_EXITC forces an error in GEOS-Chem and quits.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_EXITC (T_EXIT)
!
! !USES:
!
    USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) ::  T_EXIT
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CALL ERROR_STOP( T_EXIT, 'fjx_interface_mod.F90' )

  END SUBROUTINE GC_EXITC
!EOC
END MODULE FJX_INTERFACE_MOD
