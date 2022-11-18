!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fjx_interface_gc_mod.F90
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
MODULE FJX_INTERFACE_GC_MOD
!
! !USES:
!
  USE CMN_FJX_MOD
  USE FAST_JX_MOD
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
  PUBLIC  :: GC_EXITC ! not used... move exitc back???
  PUBLIC  :: CALC_AOD
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
  ! Flags for certain photo-reactions that will be adjusted by
  ! subroutine PHOTRATE_ADJ, which is called by FlexChem (bmy 3/29/16)
  INTEGER, PUBLIC :: RXN_O2    = -1   ! O2  + jv --> O   + O
  INTEGER, PUBLIC :: RXN_O3_1  = -1   ! O3  + hv --> O2  + O
  INTEGER, PUBLIC :: RXN_O3_2  = -1   ! O3  + hv --> O2  + O(1D)
  INTEGER, PUBLIC :: RXN_H2SO4 = -1   ! SO4 + hv --> SO2 + 2OH
  INTEGER, PUBLIC :: RXN_NO2   = -1   ! NO2 + hv --> NO  + O

  INTEGER, PUBLIC :: RXN_JHNO3  = -1   ! HNO3 + hv --> OH + NO2
  INTEGER, PUBLIC :: RXN_JNITSa = -1   ! NITs  + hv --> HNO2
  INTEGER, PUBLIC :: RXN_JNITSb = -1   ! NITs  + hv --> NO2
  INTEGER, PUBLIC :: RXN_JNITa  = -1   ! NIT + hv --> HNO2
  INTEGER, PUBLIC :: RXN_JNITb  = -1   ! NIT + hv --> NO2

  ! Needed for UCX_MOD
  INTEGER, PUBLIC :: RXN_NO    = -1
  INTEGER, PUBLIC :: RXN_NO3   = -1
  INTEGER, PUBLIC :: RXN_N2O   = -1

  ! For Hg chem
  INTEGER, PUBLIC :: RXN_BrO   = -1
  INTEGER, PUBLIC :: RXN_ClO   = -1

  ! Species ID flags
  INTEGER :: id_CH2IBr, id_IBr,  id_CH2ICl, id_ICl,   id_I2
  INTEGER :: id_HOI,    id_IO,   id_OIO,    id_INO,   id_IONO
  INTEGER :: id_IONO2,  id_I2O2, id_CH3I,   id_CH2I2, id_I2O4
  INTEGER :: id_I2O3

  ! Needed for scaling JNIT/JNITs photolysis to JHNO3
  REAL(fp)      :: JscaleNITs, JscaleNIT, JNITChanA, JNITChanB

CONTAINS

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
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE inquireMod,     ONLY : findFreeLUN
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
    LOGICAL            :: amIRoot
    LOGICAL            :: notDryRun, dryRun
    LOGICAL            :: is_fullchem, is_mercury
    INTEGER            :: JXUNIT, J, NJXX, PhotoId
    REAL(fp)           :: ND64MULT

    ! Strings
    CHARACTER(LEN=6)   :: TITLEJXX(JVN_)
    CHARACTER(LEN=255) :: FILENAME
    CHARACTER(LEN=255) :: fjx_data_dir
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_FJX begins here!
    !=================================================================

    ! Initialize
    RC          = GC_SUCCESS
    notDryRun   = ( .not. Input_Opt%DryRun )
    ErrMsg      = ''
    ThisLoc     = ' -> at Init_FJX (in module GeosCore/fast_jx_mod.F90)'

    ! Store Input_Opt options
    amIRoot          = Input_Opt%amIRoot
    dryRun           = Input_Opt%DryRun
    fjx_data_dir     = Input_Opt%FASTJX_DATA_DIR
    is_fullchem      = Input_Opt%ITS_A_FULLCHEM_SIM
    is_mercury       = Input_Opt%ITS_A_MERCURY_SIM

#ifdef CLOUDJ
    ! For now, if using cloud-j we are still using this init function
    ! Overwrite certain values from Cloud-J for consistency with GEOS-Chem
    ZZHT  = 5.e+5_fp
    RAD   = 6375.e+5_fp
    ATAU  = 1.120e+0_fp
    ATAU0 = 0.010e+0_fp
#endif

    ! Skip these operations when running in dry-run mode
    IF ( notDryRun ) THEN

       ! Define species IDs
       id_CH2IBr   = IND_('CH2IBr'  )
       id_IBr      = IND_('IBr'     )
       id_CH2ICl   = IND_('CH2ICl'  )
       id_ICl      = IND_('ICl'     )
       id_I2       = IND_('I2'      )
       id_HOI      = IND_('HOI'     )
       id_IO       = IND_('IO'      )
       id_OIO      = IND_('OIO'     )
       id_INO      = IND_('INO'     )
       id_IONO     = IND_('IONO'    )
       id_IONO2    = IND_('IONO2'   )
       id_I2O2     = IND_('I2O2'    )
       id_CH3I     = IND_('CH3i'    )
       id_CH2I2    = IND_('CH2I2'   )
       id_I2O4     = IND_('I2O4'    )
       id_I2O3     = IND_('I2O3'    )

       ! Print info
       IF ( Input_Opt%amIRoot ) THEN
          write(6,*) ' Initializing Fast-JX v7.0 standalone CTM code.'
!ewl: is this needed?
!          if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
!             ErrMsg =  ' INIT_FJX: invalid no. wavelengths'
!             CALL GC_Error( ErrMsg, RC, ThisLoc )
!             RETURN
!          endif
       ENDIF

#if defined( MODEL_CESM )
       JXUNIT = 0
       IF ( Input_Opt%amIRoot ) JXUNIT = findFreeLUN()
#else
       JXUNIT = findFreeLUN()
#endif

    ENDIF

    !=====================================================================
    ! Read in fast-J X-sections (spectral data)
    !=====================================================================
    FILENAME = TRIM(fjx_data_dir) // 'FJX_spec.dat'

    ! Read file, or just print filename if we are in dry-run mode
    CALL RD_XXX( JXUNIT, TRIM(FILENAME), amIRoot, dryRun, RC)

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
    FILENAME = TRIM(fjx_data_dir) // 'jv_spec_mie.dat'

    ! Read data
    CALL RD_MIE( JXUNIT, TRIM(FILENAME), amIRoot, dryRun, &
                 Input_Opt%LBRC, RC )

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
    CALL RD_AOD( JXUNIT, amIRoot, dryRun, fjx_data_dir, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in FAST-JX routine "RD_AOD"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Now calculate the required wavelengths in the LUT to calculate
    ! the requested AOD, and set up MIEDX array to interpret between GC
    ! and FJX aerosol indexing
    IF ( notDryRun ) THEN
       CALL CALC_AOD( Input_Opt )
       CALL SET_AER( amIRoot )
    ENDIF

    !=====================================================================
    ! Read in T & O3 climatology used to fill e.g. upper layers
    ! or if O3 not calc.
    !=====================================================================
    nc_dir  = TRIM( Input_Opt%CHEM_INPUTS_DIR ) // 'FastJ_201204/'
    FILENAME = TRIM( nc_dir ) // 'fastj.jv_atms_dat.nc'

    CALL RD_PROF_NC( amIRoot, dryRun, TRIM(FILENAME), RC )

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
    FILENAME = TRIM( Input_Opt%FAST_JX_DIR ) // 'FJX_j2j.dat'

    ! Read mapping information
    CALL RD_JS_JX( JXUNIT, TRIM(FILENAME), TITLEJXX, NJXX, amIRoot, &
                   dryRun, is_fullchem, is_mercury, RC )

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

          IF ( J == Rxn_O3_2 ) THEN

             !------------------------------------------------------------
             ! O3 + hv = O + O(1D)
             !
             ! Save this as JO3_O1D in the nPhotol+1 slot
             !------------------------------------------------------------
             GC_Photo_Id(J) = State_Chm%nPhotol + 1

          ELSE IF ( J == Rxn_O3_1 ) THEN

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
                WRITE(6, 200) RNAMES(J), J, GC_Photo_Id(J), JFACTA(J)
200             FORMAT( a10, ':', i7, 2x, i7, 2x, f7.4 )
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
  SUBROUTINE FAST_JX( Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD,       ONLY : NDUST, NRH
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

!==============================================================================
!
! !INPUT PARAMETERS:
!
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
    INTEGER       :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR, L, N
    INTEGER       :: IOPT, LCHEM
    REAL(fp)      :: CSZA, PRES, SFCA, YLAT,  O3_TOMS
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

    LOGICAL       :: AOD999    ! AOD calculated how? (1: 550 nm, 0: 999 nm)

    LOGICAL, SAVE :: FIRST = .true.
    LOGICAL       :: prtDebug

    ! Species ID flags
    INTEGER, SAVE :: id_O3

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! FAST_JX begins here!
    !=================================================================

    ! Initialize
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at Fast_JX (in module GeosCore/fast_jx_mod.F)'
    prtDebug  = ( Input_Opt%LPRT .and. Input_Opt%amIRoot)

    ! Get day of year (0-365 or 0-366)
    DAY_OF_YR = GET_DAY_OF_YEAR()

    ! Get current month
    MONTH     = GET_MONTH()

    ! Get day of month
    DAY       = GET_DAY()

    ! Assume AOD calculated at 999 nm
    AOD999    = .TRUE.

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
    !$OMP PRIVATE( LCHEM,   OPTAER, N,         IOPT             ) &
    !$OMP PRIVATE( OPTDUST, OPTD,   CLDF1D                      ) &
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

       ! Column cloud fraction (not less than zero)
       CLDF1D = State_Met%CLDF(NLON,NLAT,1:State_Grid%NZ)
       WHERE ( CLDF1D < 0e+0_fp ) CLDF1D = 0e+0_fp

       ! NOTE: For GEOS-FP and MERRA-2 met fields, the optical
       ! depth is the in-cloud optical depth.  At this point it has
       ! NOT been multiplied by cloud fraction yet.  Therefore, we can
       ! just apply the linear overlap formula as written above 
       ! (i.e. multiply by cloud fraction), or the approximate random 
       ! overlap formula as written (i.e. multiply by cloud fraction^1.5). 
#if defined( USE_LINEAR_OVERLAP )
       ! %%%% CLOUD OVERLAP: LINEAR ASSUMPTION %%%%
       ! Directly use OPTDEPTH = TAUCLD * CLDTOT
       OPTD = OPTD * CLDF1D
#elif defined( USE_APPROX_RANDOM_OVERLAP )
       ! %%%% CLOUD OVERLAP: APPROX RANDOM OVERLAP ASSUMPTION %%%%
       ! Use OPTDEPTH = TAUCLD * CLDTOT**1.5
       OPTD = OPTD * ( CLDF1D )**1.5e+0_fp
#endif

       ! Call FAST-JX routines to compute J-values
       CALL PHOTO_JX( CSZA,       SFCA,      P_CTM,  T_CTM,  &
                      O3_CTM,     O3_TOMS,   AOD999, OPTAER, &
                      OPTDUST,    OPTD,      NLON,   NLAT,   &
                      YLAT,       DAY_OF_YR, MONTH,  DAY,    &
                      Input_Opt, State_Diag, State_Grid, State_Met )

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
  SUBROUTINE PHOTRATE_ADJ( Input_Opt, State_Diag, State_Met,                 &
                           I,         J,          L,                         &
                           FRAC,      RC                                    )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input_Options object
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
! !IROUTINE: gc_exitc
!
! !DESCRIPTION: Subroutine GC_EXITC forces an error in GEOS-Chem and quits.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXITC_GC (T_EXIT)
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
    CALL ERROR_STOP( T_EXIT, 'fast_jx_mod.F90' )

  END SUBROUTINE GC_EXITC
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
    USE CMN_FJX_MOD, ONLY : NWVAA, NWVAA0, WVAA
    USE CMN_FJX_MOD, ONLY : IWVSELECT
    USE CMN_FJX_MOD, ONLY : IRTWVSELECT
    USE CMN_FJX_MOD, ONLY : ACOEF_WV, BCOEF_WV, CCOEF_WV
    USE CMN_FJX_MOD, ONLY : ACOEF_RTWV, BCOEF_RTWV, CCOEF_RTWV
    USE CMN_FJX_MOD, ONLY : NWVREQUIRED, IWVREQUIRED
    USE CMN_FJX_MOD, ONLY : NRTWVREQUIRED, IRTWVREQUIRED
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

END MODULE FJX_INTERFACE_GC_MOD


