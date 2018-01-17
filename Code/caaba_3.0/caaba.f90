! Time-stamp: <2011-02-14 16:34:55 sander>

! CAABA = Chemistry As A Boxmodel Application

! CAABA is a box model that uses MECCA chemistry, plus simplified
! calculations for emission, deposition, and photolysis

! Authors:
! Rolf Sander, MPICH, Mainz, 2003-2010
! Hella Riede, MPICH, Mainz, 2007

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

!*****************************************************************************

MODULE caaba_module

  USE messy_main_constants_mem,   ONLY: DP, HLINE1, HLINE2

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: modstr = 'caaba'  ! name of module

  INTEGER :: ncid_messy

CONTAINS

  !***************************************************************************

  SUBROUTINE caaba_read_nml(status, iou)

    ! read coupling namelist (based on dradon_read_nml_cpl by P. Joeckel)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close, &
                                strcrack, ucase, &
                                init_convect_tables, & ! mz_hr_20100704
                                cair_trad, cair_wmo ! mz_hr_20100705
    USE messy_main_constants_mem, ONLY: R_gas, N_A, OneDay, &
                                        STRLEN_SHORT, STRLEN_MEDIUM
    USE caaba_mem, ONLY: &
      cair, runtime, timesteplen, l_input_jval, l_runtime_str, l_spechum, &
      ! namelist:
      USE_JVAL, USE_MECCA, USE_READJ,                & ! MESSy submodels
      USE_SAPPHO, USE_SEMIDEP, USE_TRAJECT,          & ! MESSy submodels
#ifdef E4CHEM
      USE_E4CHEM,                                    & ! MESSy submodels
#endif
      init_scenario, photo_scenario,                 & ! scenarios
      emission_scenario, drydep_scenario,            & ! ...
      temp, press, relhum, zmbl,                     & ! CAABA (meteorology)
      l_ignore_relhum, l_psat_emac, l_relhum_wmo,    & ! CAABA (meteorology)
      degree_lat, degree_lon, l_ff,                  & ! CAABA (location)
      model_start_day, runtime_str, timesteplen_str, & ! CAABA (time)
      Ca_precip, init_spec,                          & ! MECCA-specific
      photrat_channel, l_skipkpp,                    & ! MECCA-specific
      l_steady_state_stop,                           & ! MECCA-specific
      init_j, init_j_index,                          & ! READJ-specific
      l_injectNOx, t_NOxon, t_NOxoff,                & ! SEMIDEP-specific
      runlast, input_physc, input_jval,              & ! TRAJECT-specific
      vlon, vlat, vpress, vtemp,                     & ! TRAJECT-specific
      vrelhum, vspechum                                ! TRAJECT-specific

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, TRIM

    ! I/O:
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    ! local:
    CHARACTER(LEN=*), PARAMETER :: substr = 'caaba_read_nml'
    INTEGER, PARAMETER :: MAX_SCENARIOS = 12
    CHARACTER(LEN=12), PARAMETER, DIMENSION(MAX_SCENARIOS) :: &
      list_of_scenarios = (/ &
      '            ', 'FF_ANTARCTIC', 'FF_ARCTIC   ', 'FREE_TROP   ', &
      'HOOVER      ', 'LAB         ', 'LAB_C15     ', 'MBL         ', &
      'MIM2        ', 'OOMPH       ', 'STRATO      ', 'AUS_FOR     ' /)
    LOGICAL :: l_init_scenario_ok     = .FALSE.
    LOGICAL :: l_photo_scenario_ok    = .FALSE.
    LOGICAL :: l_emission_scenario_ok = .FALSE.
    LOGICAL :: l_drydep_scenario_ok   = .FALSE.
    LOGICAL :: lex   ! file exists?
    LOGICAL :: l_ok
    INTEGER :: fstat ! file status
    INTEGER :: i
    INTEGER                        :: nosub
    CHARACTER(LEN=STRLEN_SHORT)    :: tsunit    = ''! time step length unit
    CHARACTER(LEN=STRLEN_SHORT)    :: rtunit    = ''! runtime_str unit
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: field => NULL()

    NAMELIST /CAABA/ &
      USE_JVAL, USE_MECCA, USE_READJ,                & ! MESSy submodels
      USE_SAPPHO, USE_SEMIDEP, USE_TRAJECT,          & ! MESSy submodels
#ifdef E4CHEM
      USE_E4CHEM,                                    & ! MESSy submodels
#endif
      init_scenario, photo_scenario,                 & ! scenarios
      emission_scenario, drydep_scenario,            & ! ...
      temp, press, relhum, zmbl,                     & ! CAABA (meteorology)
      l_ignore_relhum, l_psat_emac, l_relhum_wmo,    & ! CAABA (meteorology)
      degree_lat, degree_lon, l_ff,                  & ! CAABA (location)
      model_start_day, runtime_str, timesteplen_str, & ! CAABA (time)
      Ca_precip, init_spec,                          & ! MECCA-specific
      photrat_channel, l_skipkpp,                    & ! MECCA-specific
      l_steady_state_stop,                           & ! MECCA-specific
      init_j, init_j_index,                          & ! READJ-specific
      l_injectNOx, t_NOxon, t_NOxoff,                & ! SEMIDEP-specific
      runlast, input_physc, input_jval,              & ! TRAJECT-specific
      vlon, vlat, vpress, vtemp,                     & ! TRAJECT-specific
      vrelhum, vspechum                                ! TRAJECT-specific

    status = 1
    CALL read_nml_open(lex, substr, iou, 'CAABA', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist
    READ(iou, NML=CAABA, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CAABA', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    CALL read_nml_close(substr, iou, modstr)
    status = 0  ! no error

    WRITE(*,*)
    WRITE(*,*) HLINE2
    WRITE(*,*) 'Selected MESSy submodels:'
    IF (USE_JVAL)    WRITE(*,*) '  JVAL'
    IF (USE_MECCA)   WRITE(*,*) '  MECCA'
    IF (USE_READJ)   WRITE(*,*) '  READJ'
    IF (USE_SAPPHO)  WRITE(*,*) '  SAPPHO'
    IF (USE_SEMIDEP) WRITE(*,*) '  SEMIDEP'
    IF (USE_TRAJECT) WRITE(*,*) '  TRAJECT'
#ifdef E4CHEM
    IF (USE_E4CHEM)  WRITE(*,*) '  E4CHEM'
#endif
    WRITE(*,*) HLINE2

    ! scenarios:
    DO i=1, MAX_SCENARIOS
      IF (TRIM(list_of_scenarios(i))==TRIM(init_scenario)) &
        l_init_scenario_ok = .TRUE.
      IF (TRIM(list_of_scenarios(i))==TRIM(photo_scenario)) &
        l_photo_scenario_ok = .TRUE.
      IF (TRIM(list_of_scenarios(i))==TRIM(emission_scenario)) &
        l_emission_scenario_ok = .TRUE.
      IF (TRIM(list_of_scenarios(i))==TRIM(drydep_scenario)) &
        l_drydep_scenario_ok = .TRUE.
    ENDDO
    WRITE(*,*) 'Selected scenarios:'
    IF (l_init_scenario_ok) THEN
      WRITE(*,*) '  Init:       ', TRIM(init_scenario)
    ELSE
      WRITE(*,*) 'ERROR, unknown init scenario ', TRIM(init_scenario)
      STOP
    ENDIF
    IF (l_photo_scenario_ok) THEN
      WRITE(*,*) '  Photo:      ', TRIM(photo_scenario)
    ELSE
      WRITE(*,*) 'ERROR, unknown photo scenario ', TRIM(photo_scenario)
      STOP
    ENDIF
    IF (l_emission_scenario_ok) THEN
      WRITE(*,*) '  Emission:   ', TRIM(emission_scenario)
    ELSE
      WRITE(*,*) 'ERROR, unknown emission scenario ', TRIM(emission_scenario)
      STOP
    ENDIF
    IF (l_drydep_scenario_ok) THEN
      WRITE(*,*) '  Deposition: ', TRIM(drydep_scenario)
    ELSE
      WRITE(*,*) 'ERROR, unknown deposition scenario ', TRIM(drydep_scenario)
      STOP
    ENDIF
    WRITE(*,*) HLINE2

    WRITE(*,*) 'SATURATION PARTIAL PRESSURE OF WATER VAPOR'
    IF (l_psat_emac) THEN
      WRITE(*,*) '  using EMAC init_convect look-up tables'
      WRITE(*,*) '  alternative/default: psat_mk, '// &
                 'do not use l_psat_emac in CAABA nml'
      CALL init_convect_tables
    ELSE
      WRITE(*,*) '  calculated in psat_mk'
      WRITE(*,*) '  alternative: EMAC init_convect look-up tables, '// &
                 'use l_psat_emac in CAABA nml'
    ENDIF
    IF (l_ignore_relhum) THEN
      WRITE(*,*) 'l_ignore_relhum=TRUE, relative humidity is ignored.'
    ELSE
      WRITE(*,*) 'RELATIVE HUMIDITY DEFINITION'
      ! used for mixing ratio of water vapor<->concentration of water vapor
      ! and for conversion spechum<->relhum
      IF (l_relhum_wmo) THEN
        WRITE(*,*) '  WMO definition:'
        WRITE(*,*) '  relhum_WMO = mass mixing ratio H2O / saturation'// &
          ' mass mixing ratio H2O (in dry air)'
        WRITE(*,*) '  alternative: standard definition, do not use'// &
          ' l_relhum_wmo in CAABA nml'
      ELSE
        WRITE(*,*) '  standard definition: '
        WRITE(*,*) '  relhum = p(H2O) / psat(H2O) (in humid air)'
        WRITE(*,*) '  alternative: WMO definition, use l_relhum_wmo in'// &
          ' CAABA nml'
      ENDIF
    ENDIF

    ! cair = c(air) in [mcl/cc]
    IF (l_ignore_relhum) THEN
      ! usually, l_ignore_relhum is set to TRUE when c(H2O) is
      ! initialized directly (in SUBROUTINE x0) and not calculated from
      ! the relative humidity. In this case, the humidity is not known
      ! yet at this point, and for the calculation of cair it is assumed
      ! that the air is dry:
      cair = (N_A/1.E6) * press / (R_gas*temp)
    ELSE
      ! first, check that relhum is okay:
      IF (relhum < 0._dp) THEN
        PRINT *, 'ERROR: relhum = ', relhum, '< 0'
        STOP
      ENDIF
      IF ((relhum > 1.1_dp)) THEN
        PRINT *, 'ERROR: relhum = ', relhum, '> 1.1'
        STOP
      ENDIF
      IF (l_relhum_wmo) THEN
        cair = cair_wmo(status, relhum, temp, press, l_psat_emac)
      ELSE
        cair = cair_trad(status, relhum, temp, press, l_psat_emac)
      ENDIF
    ENDIF

    ! plausibility checks:
    IF (USE_TRAJECT.AND.(TRIM(input_physc)=='')) THEN
      WRITE(*,*) HLINE1
      PRINT *, 'ERROR: Submodel TRAJECT is still under construction.'
      PRINT *, 'If you want to use it, please contact Hella Riede.'
      PRINT *, 'Submodel TRAJECT requires input file in variable input_physc'
      WRITE(*,*) HLINE1
      STOP
    ENDIF
    IF (USE_TRAJECT.AND.USE_SEMIDEP) THEN
      PRINT *, 'Submodel TRAJECT does not work with SEMIDEP.'
      STOP
    ENDIF

    ! set runtime
    !mz_hr_20070508+
    IF (TRIM(runtime_str) /= '') THEN
      ! crack string into value and unit
      CALL strcrack(TRIM(runtime_str), " ", field, nosub)

      rtunit        = TRIM(ADJUSTL(field(2)(1:7)))
      runtime_str   = TRIM(field(1))
      READ(runtime_str, *) runtime

      CALL ucase(rtunit)
      SELECT CASE (rtunit)
        CASE ('SECONDS', 'SECOND')
          runtime = runtime/OneDay
        CASE ('MINUTES', 'MINUTE')
          runtime = runtime/1440._DP
        CASE ('HOURS', 'HOUR')
          runtime = runtime/24._DP
        CASE ('DAYS', 'DAY')
          ! already in days...
        CASE DEFAULT
          WRITE(*,*) 'Error in caaba_read_nml: unknown unit for runtime:'
          WRITE(*,*) rtunit
          STOP
      END SELECT
      l_runtime_str = .TRUE.
    ELSE
      runtime = 8._DP ! days
    ENDIF
    !mz_hr_20070508-

    ! set time step if specified in namelist
    IF (TRIM(timesteplen_str) /= '') THEN
      ! crack string into value and unit
      CALL strcrack(TRIM(timesteplen_str), " ", field, nosub)

      tsunit          = TRIM(ADJUSTL(field(2)(1:7)))
      timesteplen_str = TRIM(field(1))
      READ(timesteplen_str, *) timesteplen

      !mz_hr_20100611+
      CALL ucase(tsunit)
      SELECT CASE (tsunit)
        CASE ('HOURS', 'HOUR')
          timesteplen = timesteplen * 3600.0_DP
        CASE ('MINUTES', 'MINUTE')
          timesteplen = timesteplen * 60.0_DP
        CASE ('SECONDS', 'SECOND')
          ! already in seconds...
        CASE DEFAULT
          WRITE(*,*) 'ERROR in caaba_read_nml: unknown unit for timesteplen:'
          WRITE(*,*) tsunit
          STOP
      END SELECT

    ELSE
      timesteplen = 20._DP * 60._DP    ! seconds
    ENDIF

    IF (USE_TRAJECT) THEN
      WRITE(*,*) HLINE2
      WRITE(*,*) 'NETCDF INPUT FOR PHYSICAL DATA:'
      WRITE(*,*) TRIM(input_physc)
      IF (TRIM(input_jval)/='') THEN
        l_input_jval = .TRUE.
        WRITE(*,*) 'NETCDF INPUT FOR PHOTOLYSIS RATE COEFFICIENTS:'
        WRITE(*,*) '(AT THE MOMENT ONLY J_NO2)'
        WRITE(*,*) TRIM(input_jval)
      ENDIF
      IF (TRIM(vrelhum)/='' .AND. TRIM(vspechum)/='') THEN
        WRITE(*,*) 'Submodel TRAJECT allows only ONE active external'// &
         ' humidity variable'
        STOP
      ELSEIF (TRIM(vspechum) /= '') THEN
        l_spechum = .TRUE.
        WRITE(*,*) 'WATER VAPOUR CONTENT GIVEN AS SPECIFIC HUMIDITY'
      ELSEIF (TRIM(vrelhum) /= '') THEN
        WRITE(*,*) 'WATER VAPOUR CONTENT GIVEN AS RELATIVE HUMIDITY'
      ELSE ! default: relhum with name 'RELHUM'
        WRITE(*,*) "ASSUMED THAT WATER VAPOUR CONTENT GIVEN AS RELATIVE"// &
                   " HUMIDITY ('RELHUM')"
        vrelhum = "RELHUM"
      ENDIF
    ELSE
      IF (TRIM(input_jval)/='') THEN
        WRITE(*,*) 'ERROR: CURRENTLY, INPUT_JVAL CAN ONLY BE USED FOR'
        WRITE(*,*) 'TRAJECTORY MODEL RUNS, I.E. IF USE_TRAJECT=T'
        STOP
      ENDIF
      WRITE(*,'(A,F10.1)')     ' model_start_day = ', model_start_day
      WRITE(*,'(A,F10.1,A)')   ' runtime         = ', runtime,        ' days'
      WRITE(*,'(A,F10.1,A)')   ' timesteplen     = ', timesteplen,    ' s'
      WRITE(*,'(A,F10.1,A)')   ' latitude        = ', degree_lat,     ' degree'
      WRITE(*,'(A,F10.1,A)')   ' longitude       = ', degree_lon,     ' degree'
      WRITE(*,'(A,F10.1,A)')   ' T               = ', temp,           ' K'
      WRITE(*,'(A,F10.1,A)')   ' p               = ', press,          ' Pa'
      WRITE(*,'(A,F10.1,A)')   ' relhum          = ', 100.*relhum,    ' %'
      WRITE(*,'(A,F10.1,A)')   ' zmbl            = ', zmbl,           ' m'
      WRITE(*,'(A,1PE10.3,A)') ' c(air)          = ', cair,           ' mcl/cm3'
    ENDIF
    WRITE(*,*) HLINE2

    WRITE(*,*) 'l_skipkpp           = ', l_skipkpp
    WRITE(*,*) 'l_steady_state_stop = ', l_steady_state_stop
    WRITE(*,*) 'l_ff                = ', l_ff
    WRITE(*,*) 'Ca_precip           = ', Ca_precip, '(CaCO3 precipitation)'
    WRITE(*,*) HLINE2

    IF (TRIM(init_spec)/='') THEN
      WRITE(*,*) 'NETCDF INPUT FOR CHEMICAL TRACER INITIALIZATION:'
      WRITE(*,*) TRIM(init_spec)
      WRITE(*,*) 'NOTE THAT ONLY SPECIES BELONGING TO THE CHOSEN MECHANISM'
      WRITE(*,*) 'ARE INITIALIZED! SEE mecca.spc,'
      WRITE(*,*) 'messy_mecca_kpp_parameters.f90, OR YOUR MECCANISM FILE'
      WRITE(*,*) 'OTHER SPECIES (NOT initialized) ARE INITIALIZED WITH'
      WRITE(*,*) 'THE DEFAULT DEFINED IN messy_mecca_box.f90/x0'
    ELSE
      WRITE(*,*) 'init_spec = empty'
      WRITE(*,*) 'NO EXTERNAL INPUT FOR CHEMICAL TRACER INITIALIZATION'
    ENDIF
    WRITE(*,*) HLINE2

    IF (USE_READJ) THEN
      IF (TRIM(init_j)/='') THEN
        WRITE(*,*) 'netcdf input for READJ:'
        WRITE(*,*) TRIM(init_j)
        WRITE(*,*) 'index =', init_j_index
      ELSE
        PRINT *, 'ERROR: ', 'init_j is empty'
        STOP
      ENDIF
      WRITE(*,*) HLINE2
    ENDIF

    IF ((TRIM(photrat_channel)=='jval') .OR. &
      (TRIM(photrat_channel)=='readj') .OR. &
      (TRIM(photrat_channel)=='sappho')) THEN
      WRITE(*,*) 'PHOTOLYSIS RATE COEFFICIENTS ARE TAKEN FROM'
      WRITE(*,*) TRIM(photrat_channel)
    ELSE
      PRINT *, 'ERROR: ', TRIM(photrat_channel), &
        'is not a valid photolysis submodel'
      STOP
    ENDIF
    IF ((TRIM(photrat_channel)=='jval').AND.(.NOT.USE_JVAL)) THEN
      PRINT *, 'ERROR: photrat_channel=jval but USE_JVAL=F'
      STOP
    ENDIF
    IF ((TRIM(photrat_channel)=='readj').AND.(.NOT.USE_READJ)) THEN
      PRINT *, 'ERROR: photrat_channel=readj but USE_READJ=F'
      STOP
    ENDIF
    IF ((TRIM(photrat_channel)=='sappho').AND.(.NOT.USE_SAPPHO)) THEN
      PRINT *, 'ERROR: photrat_channel=sappho but USE_SAPPHO=F'
      STOP
    ENDIF
    WRITE(*,*) HLINE2

    IF (ASSOCIATED(field)) THEN
      DEALLOCATE(field)
      NULLIFY(field)
    ENDIF

  END SUBROUTINE caaba_read_nml

  !***************************************************************************

  SUBROUTINE calc_sza

    USE messy_main_constants_mem, ONLY: OneDay, PI
    USE messy_main_timer,         ONLY: utc2lt, julian2gregor
    USE caaba_mem,                ONLY: cossza, model_start, model_time, &
                                        degree_lon, degree_lat, localtime, &
                                        firstjan_jul, time0_jul, lyear, &
                                        lmonth, lday, lhour, lmin, lsec, &
                                        l_steady_state_stop

    IMPLICIT NONE
    INTRINSIC :: SIN, COS

    REAL(DP) :: localtime_start ! localtime at model start
    REAL(DP), PARAMETER :: InitSpr = 80._DP ! first day of spring in NH [d]
    REAL(DP), PARAMETER :: cancer  = 23.441_DP*(PI/180._DP) ! Earth's obliquity
    REAL(DP) :: rad_lat, dayreal, sodecli, localtime_jul
    INTEGER  :: status
    LOGICAL, PARAMETER :: l_groundhogday = .FALSE.

    ! local time now:
    localtime = utc2lt(status, model_time, degree_lon)
    IF (status /= 0) THEN
      WRITE(*,*) 'ERROR in utc2lt'
      STOP
    ENDIF
    ! local time at model start:
    localtime_start = utc2lt(status, model_start, degree_lon)
    IF (l_steady_state_stop) THEN
      ! for steady-state calculations, freeze the sun in the sky at the
      ! position where it was at the start of the model run:
      localtime = localtime_start
    ENDIF
    ! localtime as Julian date
    localtime_jul = time0_jul + localtime/OneDay

    ! year, month, ... of localtime_jul for jval, localtime_jul modified to
    CALL julian2gregor(localtime_jul, lyear, lmonth, lday, lhour, lmin, lsec)
    !print *, 'calc: localtime_jul = ', localtime_jul
    !print *, 'lyear  = ', lyear
    !print *, 'lmonth = ', lmonth
    !print *, 'lday   = ', lday
    !print *, 'lhour  = ', lhour
    !print *, 'lmin   = ', lmin
    !print *, 'lsec   = ', lsec

    ! day as a real value (at start of spring dayreal = 0):
    ! = (julian date)-(1st January of current year)-(start of spring)
    dayreal = localtime_jul - firstjan_jul - InitSpr

    IF (l_groundhogday) THEN
      ! repetition of the first day:
      dayreal = AINT(localtime_start/OneDay + time0_jul - firstjan_jul) + &
        MODULO(localtime/OneDay,1._DP) - InitSpr
      ! print "normal" dayreal and "groundhog" dayreal:
      ! WRITE(*,'(2F16.6)') localtime_jul-firstjan_jul-InitSpr, dayreal 
    ENDIF

    ! seasonal cycle:
    sodecli = cancer * SIN (2._DP * PI * dayreal / 365.25_DP)
    ! diurnal cycle of psi, the solar elevation angle
    rad_lat  = degree_lat * (PI/180.)
    cossza = SIN(rad_lat) * SIN(sodecli) &
      - COS(rad_lat) * COS(sodecli) * COS(2._DP * PI * dayreal)

  END SUBROUTINE calc_sza

  !***************************************************************************

  SUBROUTINE caaba_result

    USE caaba_mem,                ONLY: model_time, model_start, &
                                        model_end, USE_TRAJECT, &
                                        temp, press, cossza, &
                                        percent_done, degree_sza
    USE caaba_io,                 ONLY: write_output_file
    USE messy_main_constants_mem, ONLY: OneDay, PI, FLAGGED_BAD

    degree_sza   = ACOS(cossza)*180./PI
    percent_done = 100. * (model_time-model_start) / (model_end-model_start)
    IF (percent_done>0.) THEN
      WRITE(*,'(A,F9.4,A,F6.2,A,F6.2,A)') &
      ' day = ', model_time/OneDay, &
      '    sza = ', degree_sza, &
      '     (', percent_done, '% done)'
      IF (.NOT. USE_TRAJECT) THEN
        CALL write_output_file(ncid_messy, model_time, &
          (/ press, temp, degree_sza /) )
      ENDIF
    ELSE ! 1st output, sza not calc yet
      ! sza at model_start is meaningless
      WRITE(*,'(A,F9.4,A,F6.2,A)') &
      ' day = ', model_time/OneDay, &
      '                     (', percent_done, '% done)'
      IF (.NOT. USE_TRAJECT) THEN
        CALL write_output_file(ncid_messy, model_time, &
          (/ press, temp, FLAGGED_BAD /) )
      ENDIF
    ENDIF

  END SUBROUTINE caaba_result

  !***************************************************************************

  SUBROUTINE caaba_init

#ifndef E4CHEM
    USE messy_mecca_kpp,            ONLY: NSPEC
#else
    USE messy_mecca_kpp,            ONLY: NSPEC_mecca => NSPEC
    USE messy_e4chem,               ONLY: NSPEC_fchem => NSPEC
#endif
    USE messy_mecca,                ONLY: mecca_version => modver
    USE caaba_mem,                  ONLY: model_time, model_start,     &
                                          model_end, USE_TRAJECT,      &
                                          temp, press, caaba_version,  &
#ifdef E4CHEM
                                          USE_MECCA, USE_E4CHEM,       &
#endif
                                          c, runtime, model_start_day, &
                                          caaba_version
    USE messy_main_constants_mem,   ONLY: OneDay, TINY_DP
    USE caaba_io,                   ONLY: open_output_file

#ifdef E4CHEM
    INTEGER             :: NSPEC
#endif
    INTEGER, PARAMETER  :: iou = 999   ! I/O unit
    INTEGER             :: status ! error status

    ! set the version number of CAABA to that of MECCA:
    caaba_version = mecca_version

    PRINT *, HLINE1
    PRINT *, '*** START OF CAABA BOX MODEL RUN (VERSION '// &
      TRIM(caaba_version)//')'
    PRINT *, HLINE1
    PRINT *

    CALL caaba_read_nml(status, iou)
    IF (status /= 0) STOP

    model_start = model_start_day * OneDay
    model_time  = model_start
    model_end   = model_time + runtime * OneDay ! time end

#ifdef E4CHEM
    IF (USE_MECCA)    NSPEC = NSPEC_mecca
    IF (USE_E4CHEM)   NSPEC = NSPEC_fchem
#endif

    ALLOCATE(c(NSPEC))

    ! open output file mecca_messy.nc:
    IF (.NOT. USE_TRAJECT) THEN
      CALL open_output_file(ncid_messy, 'caaba_messy', &
        (/ 'press', 'temp ', 'sza  ' /), &
        (/ 'Pa ', 'K  ', 'deg' /) )
    ENDIF

  END SUBROUTINE caaba_init

  !***************************************************************************

  SUBROUTINE caaba_physc

    ! USE caaba_mem,                ONLY: model_time, temp, press, relhum, cair
    ! USE messy_main_constants_mem, ONLY: R_gas, N_A
    ! USE messy_main_tools,         ONLY: psat
    ! REAL(dp) :: temp_old, press_old

    CALL calc_sza

    ! If you want to modify p,T,rh during the model run, do it here:
    !
    ! temp_old = temp
    ! press_old = press
    ! temp   = myfunction(model_time)
    ! press  = myfunction(model_time)
    ! relhum = myfunction(model_time)
    !
    ! After changing temp and/or press, the concentration of 'air'
    ! molecules [mcl/cc] must be updated:
    !
    ! cair = ...
    !
    ! If you want to keep the mixing ratios [mol/mol] of chemical
    ! species constant, it is necessary to calculate their new
    ! concentrations [mcl/cc]:
    !
    ! c(:) = c(:) * (press*temp_old) / (press_old*temp)
    !
    ! Do not use this correction for water, instead calculate it directly:
    !
    ! c(ind_H2O) = cair * relhum * psat(temp) / press

  END SUBROUTINE caaba_physc

  !***************************************************************************

  SUBROUTINE caaba_finish

    USE caaba_mem, ONLY: c, USE_TRAJECT, caaba_version
    USE caaba_io,  ONLY: close_file
    USE messy_mecca,  ONLY: modver

    DEALLOCATE(c)

    IF (.NOT. USE_TRAJECT) THEN
      CALL close_file(ncid_messy)
    ENDIF

    PRINT *
    PRINT *, HLINE1
    PRINT *, '*** END OF CAABA BOX MODEL RUN (VERSION '// &
      TRIM(caaba_version)//')'
    PRINT *, HLINE1

  END SUBROUTINE caaba_finish

  !***************************************************************************

END MODULE caaba_module

!*****************************************************************************

PROGRAM caaba

  ! This is the base model layer (BML), i.e. the main box model program

  USE caaba_module,           ONLY: caaba_init, caaba_result, caaba_physc, &
                                    caaba_finish
  USE caaba_mem,              ONLY: model_time, model_end, timesteplen
  USE messy_main_control_cb,  ONLY: messy_init, messy_result, messy_physc, &
                                    messy_finish

  IMPLICIT NONE
  CHARACTER(LEN=*), PARAMETER :: statusfile = 'status.log'

  ! write exit status 1 to statusfile:
  OPEN(10, FILE=statusfile, status='UNKNOWN')
  WRITE(10,'(A)') "1"
  CLOSE(10)

  ! initialization:
  CALL caaba_init
  CALL messy_init
  CALL caaba_result
  CALL messy_result

  ! time loop:
  DO WHILE (model_time < model_end)
    CALL caaba_physc
    CALL messy_physc
    model_time = model_time + timesteplen
    CALL caaba_result
    CALL messy_result
  END DO

  ! final clean up:
  CALL messy_finish
  CALL caaba_finish

  ! program has finished successfully, write exit status 0 to statusfile:
  OPEN(10, FILE=statusfile, status='UNKNOWN')
  WRITE(10,'(A)') "0"
  CLOSE(10)

END PROGRAM caaba

!*****************************************************************************
