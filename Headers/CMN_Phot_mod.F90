
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_Phot_mod.F90
!
! !DESCRIPTION: Module CMN\_Phot\_Mod contains parameters and global variables
!  related to photolysis and radiative transfer but not used to interface
!  between Harvard chemistry and UC-Irvine photolysis programs (Fast-J/Fast-JX)
!\\
!\\
! !INTERFACE:
!
MODULE CMN_Phot_Mod
!
! !USES:
!
  USE CMN_FastJX_Mod, ONLY: A_, AN_, JVN_, WX_
  USE CMN_SIZE_MOD,   ONLY : NDUST, NAER
  USE PRECISION_MOD      ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PUBLIC
!
! !DEFINED PARAMETERS:
!
  ! New (SDE 03/28/13)
  ! Index in RAA & QAA of 999 nm wavelength
  INTEGER, PARAMETER :: IND999  = 5

  ! Mapping array from Harvard species names to UCI species names
  INTEGER            :: RINDEX(JVN_)

  ! GEOS-Chem "ModelId" corresponding to each photolysis species
  INTEGER            :: GC_Photo_Id(JVN_)

  ! Output J values
  REAL(fp), ALLOCATABLE :: ZPJ(:,:,:,:)

  !-----------------------------------------------------------------------
  ! variables used to map fast-JX J's onto CTM J's
  !-----------------------------------------------------------------------

  ! Conversion factors from photons/cm2s to W/m2
  REAL(fp), DIMENSION(WX_) :: UVXFACTOR

  !-----------------------------------------------------------------------
  ! Variables in file 'jv_spec_aod.dat' (RD_AOD)
  !-----------------------------------------------------------------------

  ! QAA_AOD: Aerosol scattering phase functions
  REAL(fp)             :: QAA_AOD(A_)

  ! WAA: 5 Wavelengths for the supplied phase functions
  REAL(fp)             :: WAA_AOD(A_)

  ! PAA: Phase function: first 8 terms of expansion
  REAL(fp)             :: PAA_AOD(8,A_)

  ! RAA: Effective radius associated with aerosol type
  REAL(fp)             :: RAA_AOD(A_)

  ! SAA: Single scattering albedo
  REAL(fp)             :: SAA_AOD(A_)

  !-----------------------------------------------------------------------
  ! Variables in file 'atmos_std.dat' (RD_PROF)
  !-----------------------------------------------------------------------

  ! T and O3 reference profiles
  REAL(fp), DIMENSION(51,18,12) :: TREF, OREF

  ! Interfacing indices for GC and FJX aerosols
  INTEGER, ALLOCATABLE  :: MIEDX(:)

  ! Dust and aerosol optical depths
  REAL(fp), ALLOCATABLE :: ODMDUST(:,:,:,:,:)
  REAL(fp), ALLOCATABLE :: ODAER(:,:,:,:,:)
  REAL(fp), ALLOCATABLE :: ISOPOD(:,:,:,:)   ! eam, 2014

  !-----------------------------------------------------------------------
  !  Variables added for RRTMG (dar, mps, 12/5/14)
  !-----------------------------------------------------------------------

  INTEGER, PARAMETER :: NWVAA   = 41     !number of wavelengths in LUT
  INTEGER, PARAMETER :: NWVAA0  = 11     !number of non-RRTMG wavelengths
  INTEGER, PARAMETER :: NWVAART = NWVAA-NWVAA0 !number of RRTMG wvs
  INTEGER, PARAMETER :: NRAA    = 7      !number of aer sizes in LUT
  INTEGER, PARAMETER :: NALBD   = 2
  INTEGER, PARAMETER :: NEMISS  = 16

  ! Now set the following in Init_CMN_Phot below (mps, 1/3/18)
  INTEGER            :: NSPAA            !number of species in LUT
  INTEGER            :: NASPECRAD        !aerosol species in RT
  INTEGER            :: NSPECRAD         !aerosol+gas species in RT

  ! New optical arrays
  REAL*8, ALLOCATABLE :: WVAA(:,:)
  REAL*8, ALLOCATABLE :: RHAA(:,:)
  REAL*8, ALLOCATABLE :: NRLAA(:,:,:)
  REAL*8, ALLOCATABLE :: NCMAA(:,:,:)
  REAL*8, ALLOCATABLE :: RDAA(:,:)
  REAL*8, ALLOCATABLE :: RWAA(:,:)
  REAL*8, ALLOCATABLE :: SGAA(:,:)
  REAL*8, ALLOCATABLE :: QQAA(:,:,:)
  REAL*8, ALLOCATABLE :: ALPHAA(:,:,:)
  REAL*8, ALLOCATABLE :: REAA(:,:)
  REAL*8, ALLOCATABLE :: SSAA(:,:,:)
  REAL*8, ALLOCATABLE :: ASYMAA(:,:,:)
  REAL*8, ALLOCATABLE :: PHAA(:,:,:,:)
  INTEGER :: IWVSELECT(2,3) !index of requested wavelengths
  INTEGER :: IRTWVSELECT(2,3) !index of requested RT wavelengths

  ! max of 3 but need 2 per wavelength if interpolating
  INTEGER :: NWVREQUIRED !number of wvs required for interpolation
  INTEGER :: IWVREQUIRED(6) !index of wavelengths for interpo.
  INTEGER :: NRTWVREQUIRED !number of wvs required for RT interpolation
  INTEGER :: IRTWVREQUIRED(6) !index of wavelengths for RT interpo.
  ! list of required wavelengths, up to max of 3 x 2

  INTEGER :: IWV1000 !Store the wavelength index for 1000nm for Fast-J

  !coefficients for interpolation of wavelength (and for RT too)
  REAL*8      :: ACOEF_WV(3),BCOEF_WV(3),CCOEF_WV(3)
  REAL*8      :: ACOEF_RTWV(3),BCOEF_RTWV(3),CCOEF_RTWV(3)
  INTEGER, ALLOCATABLE :: SPECMASK(:) !list of binary switches for different
                                      !species flux output

  ! RH indices
  INTEGER, ALLOCATABLE :: IRHARR(:,:,:)

#ifdef RRTMG
  !to pass to RT code
  !one for each hydrophilic/hydrophobic aerosol and optical dust bin
  !and also sulfate, nitrate and ammonia are separate too
  REAL*8, ALLOCATABLE :: RTODAER(:,:,:,:,:)
  REAL*8, ALLOCATABLE :: RTSSAER(:,:,:,:,:)
  REAL*8, ALLOCATABLE :: RTASYMAER(:,:,:,:,:)
#endif

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Cmn_Phot
!
! !DESCRIPTION: Routine INIT\_CMN\_Phot initializes quantities based on
!  the grid-independent size parameters.
!\\
!\\
! !INTERFACE:

  SUBROUTINE Init_CMN_Phot( Input_Opt, State_Grid, RC )
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
    ! INIT_CMN_Phot begins here!
    !=================================================================

    ! For RRTMG:
    NSPAA     = 8   ! number of species in LUT
    NASPECRAD = 16  ! aerosol species in RT
    NSPECRAD  = 18  ! aerosol+gas species in RT

    IF ( .not. Input_Opt%DryRun ) THEN
       !-----------------------------------------------------------------------
       !  Allocate arrays
       !-----------------------------------------------------------------------

       ALLOCATE( ZPJ( State_Grid%NZ, JVN_, State_Grid%NX, State_Grid%NY), &
            STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       ZPJ = 0e+0_fp

       ALLOCATE( ODMDUST( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
            NWVAA, NDUST), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       ODMDUST = 0e+0_fp

       ALLOCATE( ODAER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, NWVAA, NAER),&
            STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       ODAER = 0e+0_fp

       ALLOCATE( MIEDX(AN_), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       MIEDX = 0

       ! Allocate array for isoprene SOA AOD (eam, 2014):
       ALLOCATE( ISOPOD( State_Grid%NX, State_Grid%NY, State_Grid%NZ, NWVAA), &
            STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       ISOPOD = 0e+0_fp

       !-----------------------------------------------------------------------
       !  Variables added for RRTMG (dar, mps, 12/5/14)
       !-----------------------------------------------------------------------

       ALLOCATE( IRHARR( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
            STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       IRHARR = 0d0

       ALLOCATE( WVAA(NWVAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       WVAA = 0d0

       ALLOCATE( RHAA(NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       RHAA = 0d0

       ALLOCATE( NRLAA(NWVAA,NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       NRLAA = 0d0

       ALLOCATE( NCMAA(NWVAA,NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       NCMAA = 0d0

       ALLOCATE( RDAA(NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       RDAA = 0d0

       ALLOCATE( RWAA(NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       RWAA = 0d0

       ALLOCATE( SGAA(NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SGAA = 0d0

       ALLOCATE( QQAA(NWVAA,NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       QQAA = 0d0

       ALLOCATE( ALPHAA(NWVAA,NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       ALPHAA = 0d0

       ALLOCATE( REAA(NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       REAA = 0d0

       ALLOCATE( SSAA(NWVAA,NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SSAA = 0d0

       ALLOCATE( ASYMAA(NWVAA,NRAA,NSPAA), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       ASYMAA = 0d0

       ALLOCATE( PHAA(NWVAA,NRAA,NSPAA,8), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       PHAA = 0d0

       ALLOCATE( SPECMASK(NSPECRAD), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SPECMASK = 0

#ifdef RRTMG
       ! +2 to split SNA into SU, NI and AM
       ALLOCATE( RTODAER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
            NWVAA,NAER+2+NDUST), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       RTODAER = 0d0

       ALLOCATE( RTSSAER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
            NWVAA, NAER+2+NDUST ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       RTSSAER = 0d0

       ALLOCATE( RTASYMAER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
            NWVAA, NAER+2+NDUST ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       RTASYMAER = 0d0
#endif
    ENDIF

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE Init_CMN_Phot
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Cmn_Phot
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_Phot deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_CMN_Phot( RC )
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
    IF ( ALLOCATED( ZPJ       ) ) DEALLOCATE( ZPJ       )
    IF ( ALLOCATED( ODMDUST   ) ) DEALLOCATE( ODMDUST   )
    IF ( ALLOCATED( ODAER     ) ) DEALLOCATE( ODAER     )
    IF ( ALLOCATED( MIEDX     ) ) DEALLOCATE( MIEDX     )
    IF ( ALLOCATED( ISOPOD    ) ) DEALLOCATE( ISOPOD    )
    IF ( ALLOCATED( IRHARR    ) ) DEALLOCATE( IRHARR    )
    IF ( ALLOCATED( WVAA      ) ) DEALLOCATE( WVAA      )
    IF ( ALLOCATED( RHAA      ) ) DEALLOCATE( RHAA      )
    IF ( ALLOCATED( NRLAA     ) ) DEALLOCATE( NRLAA     )
    IF ( ALLOCATED( NCMAA     ) ) DEALLOCATE( NCMAA     )
    IF ( ALLOCATED( RDAA      ) ) DEALLOCATE( RDAA      )
    IF ( ALLOCATED( RWAA      ) ) DEALLOCATE( RWAA      )
    IF ( ALLOCATED( SGAA      ) ) DEALLOCATE( SGAA      )
    IF ( ALLOCATED( QQAA      ) ) DEALLOCATE( QQAA      )
    IF ( ALLOCATED( ALPHAA    ) ) DEALLOCATE( ALPHAA    )
    IF ( ALLOCATED( REAA      ) ) DEALLOCATE( REAA      )
    IF ( ALLOCATED( SSAA      ) ) DEALLOCATE( SSAA      )
    IF ( ALLOCATED( ASYMAA    ) ) DEALLOCATE( ASYMAA    )
    IF ( ALLOCATED( PHAA      ) ) DEALLOCATE( PHAA      )
    IF ( ALLOCATED( SPECMASK  ) ) DEALLOCATE( SPECMASK  )
#ifdef RRTMG
    IF ( ALLOCATED( RTODAER   ) ) DEALLOCATE( RTODAER   )
    IF ( ALLOCATED( RTSSAER   ) ) DEALLOCATE( RTSSAER   )
    IF ( ALLOCATED( RTASYMAER ) ) DEALLOCATE( RTASYMAER )
#endif

    ! Return successfully
    RC = GC_SUCCESS

  END SUBROUTINE Cleanup_CMN_Phot
!EOC
END MODULE CMN_Phot_Mod
