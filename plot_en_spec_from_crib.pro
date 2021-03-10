;+
;PROCEDURE: plot_en_spec_from_crib
;PURPOSE: 
;  It is a crib sheet for plotting cluster cis/codif energy spectra using
;  the TPLOT package
;
;INPUT
; PARAMETERS:   sat:        Satellite number
;               inst:       Instrument (0: CODIF, 1: HIA)
;               time:       Date/time in tplot format
;               timespan:   Time span
;               units_name: Units for energy spectra
;                           'Counts', 'NCOUNTS', 'RATE', 'NRATE', 
;                            'DIFF FLUX', 'EFLUX'
;               angle:      Angle range to sum over [thetarange, phirange]
;
;                           THETA_88: [-78.75, -56.25, -33.75, -11.25,
;                                      11.25, 33.75, 56.25, 78.75] 
; 
;                           PHI_88:   [11.25, 22.50, 33.75, 45.00, 56.25,
;                                      67.50, 78.75, 101.25, 112.50, 123.75,
;                                      135.00, 146.25, 157.50, 168.75, 191.25,
;                                      202.50, 213.75, 225.00, 236.25, 247.50,
;                                      258.75, 281.25, 292.50, 303.75, 315.00,
;                                      326.25, 337.50, 348.75]
; 
;                           THETA_24: [-78.75, -56.25, -33.75, -11.25,
;                                      11.25, 33.75, 56.25, 78.75] 
;
;                           PHI_24:   [45.0, 90.0, 135.0, 225.0, 270.0, 315.0]
;
; KEYWORDS: OLD_EFF -> this keyword allows use of the old efficiency
;                      file format
;           HR -> produce energy spectra at phi angle resolution
;
; CREATED BY: C Mouikis
;
; LAST MODIFICATION: 10/26/01
;
; MODIFICATION HISTORY: 
;   09/26/01: If pre-processed files exist, they are read otherwise
;             the energy spectra are calculated from the L1 files
;   10/26/01: The RECALC keyword is introduced allowing the
;             possibility to recalculate the energy spectra even if a
;             pre-processed file does exist
;   10/26/01: If tplot save file is compressed it uncompresses it
;             in the directory declared by the CCAT_TEMP
;             environment variable
;   08/05/02: BKG keyword introduced
;   12/31/02: HR  keyword introduced
;   03/05/03: SPILL keyword introduced
;   09/11/03: INCRATES keyword introduced
;   06/05/06: Added phi range in variable name
;   06/22/06: Phi ongle in variable name included
;   11/22/06: Combine spectra variable name corrected
;-

PRO plot_en_spec_from_crib, $
                            sat,                 $
                            specie,              $
                            inst,                $
                            units_name,          $
                            angle,               $
                            eff_table,           $
                            recalc = recalc,     $
                            OLD_EFF=OLD_EFF,     $
                            BKG=BKG,             $
                            INCRATES=INCRATES,   $
                            HR=HR,               $
                            SPILL=SPILL, $
                            data_input = data_input
  
  COMMON get_error, get_err_no, get_err_msg, default_verbose
  
  get_err_no = 0 & get_err_msg=''

  ; Check if sat and specie arrays have the same number of elements
  IF n_elements(sat) NE n_elements(specie) THEN BEGIN
    print, 'The sat array and the specie array MUST have'
    print, 'the same number of elements.'
    stop
  ENDIF
  
  ; Check if input parameters are within the allowed limits
  FOR ic = 0, n_elements(sat)-1 DO BEGIN
    IF sat(ic) LT 1 OR sat(ic) GT 4 THEN BEGIN
      print, 'Allowed values for sat are: 1,2,3 or 4'
      STOP
    ENDIF
    IF specie(ic) LT 0 OR specie(ic) GT 3 THEN BEGIN
      print, 'Allowed values for specie are: 0,1,2 or 3'
      STOP
    ENDIF
  ENDFOR
  
  
  IF eff_table NE 0 AND eff_table NE 1 THEN BEGIN
    print, 'Allowed values for eff_table are: 0 or 1'
    STOP
  ENDIF

  specie_str = ['H!U+!N','He!U++!N','He!U+!N','O!U+!N']
  
  ; Loop over all sat/specie combinations
  FOR ii=0,N_ELEMENTS(specie)-1 DO BEGIN

    get_err_no = 0 ; reset error indicator

    ; Read pre-processed file
    IF NOT KEYWORD_SET(RECALC) THEN BEGIN
      IF eff_table EQ 0 THEN BEGIN
        IF units_name EQ 'DIFF FLUX' OR units_name EQ 'EFLUX' THEN BEGIN
          r_en_spec_onefile, $
            sat(ii), specie(ii), inst, units_name, angle, eff_table
        ENDIF ELSE get_err_no = 1
      ENDIF ELSE get_err_no = 1

      IF get_err_no EQ 0 THEN BEGIN
        
      ;----------------------------------------------------------------
      ; Limit time interval
      ;----------------------------------------------------------------
;      get_timespan,time_interval
      
;      t_s=gettime(time_interval(0)) ; start time in tplot-time
;      t_e=gettime(time_interval(1)) ; end time in tplot-time  
      
      ;get_data, 
      
;      tind=where(dtime GE t_s AND dtime LT t_e, tc) ; check if there are data
;      IF tc LT 2 THEN BEGIN     ; changed to include    ; (more than one packet)
;        get_err_no = 1          ; one packet            ; for time interval requested
;        get_err_msg = 'No data in time interval'
;        print, get_err_msg
;        RETURN
;      ENDIF
      
;      ind_s=tind(0)             ; index that corresponds to start time
;      ind_e=tind(tc-1)          ; index that corresponds to end time
      ;----------------------------------------------------------------
      
      
        
        GOTO, read
      ENDIF
    ENDIF
    
    RECALC = 1
    ; Products that correspond to a specie  

;Jing: Add the HIA prod number if inst is set to be 1,
;  and add the the "return", if inst is not set correctly.
    IF inst EQ 0 OR inst EQ 2 THEN BEGIN 
        CASE specie(ii) OF
            0: prod = [12, 13]
            1: prod = [15, 16]
            2: prod = [17, 18, 46, 48] ; [17,18,46,48]
            3: prod = [17, 18, 47, 49] ; [17,18,47,49]
            ELSE: BEGIN
                print, 'Specie numbers must be 0 for H+, 1 for He++,'
                print, '2 for He+ or 3 for O+'
                stop
            END
        ENDCASE
    ENDIF ELSE BEGIN 
        IF inst EQ 1 THEN prod = [6, 15, 23] $ ;17
        ELSE BEGIN  
            print, 'Instrument number must be 0 (CODIF) or 1 (HIA) or 2 (RPA)'
            RETURN
        ENDELSE 
    ENDELSE 
    ;Loop over all products that correspond to a specie
    FOR jj=0,n_elements(prod)-1 DO BEGIN

      name = 'ENSPEC' + $
        '_SC' + strcompress(sat(ii),/remove_all) + $
        '_IN' + STRCOMPRESS(inst,/REMOVE_ALL) + $
        '_PHI' + STRCOMPRESS(FIX(angle(0,1)), /REMOVE_ALL) + $
        '_' + STRCOMPRESS(FIX(angle(1,1)), /REMOVE_ALL) + $
        '_UN' + STRUPCASE(strcompress(units_name,/remove_all)) + $
        '_PR' + strcompress(string(prod(jj)),/remove_all) + $
        '_SP' + strcompress(specie(ii),/remove_all) + $
        '_ET' + strcompress(eff_table,/remove_all) ;name of TPLOT structure
      
      IF NOT KEYWORD_SET(HR) THEN BEGIN
        get_spec_cis, $
          sat(ii), $
          inst, $
          prod(jj), $
          specie=specie(ii), $
          eff_table, $
          angle=angle,      $
          units=units_name, $
          name=name, $
          OLD_EFF=OLD_EFF, $
          BKG=BKG, $
          INCRATES=INCRATES, $
          SPILL=SPILL
      ENDIF ELSE BEGIN
        get_spec_cis_hr, $
          sat(ii), $
          inst, $
          prod(jj), $
          specie=specie(ii), $
          eff_table, $
          angle=angle,      $
          units=units_name, $
          name=name, $
          OLD_EFF=OLD_EFF, $
          BKG=BKG, $
          INCRATES=INCRATES, $
          HR=HR
      ENDELSE

      IF get_err_no GT 0 THEN GOTO, next   
      
      next:
    ENDFOR

    ; Combine the different product of the same specie
    combine_en_spec, sat(ii), inst, specie(ii), $
      strcompress(STRUPCASE(units_name), /remove_all), name_all, eff_table, $
      angle=angle
    
    ; Set plot attributes
    CASE STRUPCASE(units_name) OF
      'COUNTS': uname = 'COUNTS'
      'BCOUNTS': uname = '1/tof_bin'
      'NCOUNTS': uname = '1/bin'
      'RATE': uname = '1/s'
      'NRATE': uname = '1/s-bin'
      'EFLUX': uname = 'eV/cm!E2!N-s-sr-eV'
      'DIFF FLUX': uname = '1/cm!E2!N-s-sr-(eV/e)'
      'DIST FUNC': uname = 's!E3!N/cm!E3!N-km!E3!N'
    ENDCASE
      
    options, name_all, 'spec',1
    options, name_all, 'x_no_interp',1
    options, name_all, 'y_no_interp',1
    options, name_all, 'ytitle', $
      'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
      specie_str(specie(ii))+' (eV)'
    options, name_all, 'ztitle', uname
    ylim,    name_all,  20., 4.5e4, 1
    
    CASE specie(ii) OF
      0: zlim, name_all,  1.e+0, 1.e3, 1
      1: zlim, name_all,  1.e-1, 1.e2, 1
      2: zlim, name_all,  1.e-2, 1.e1, 1
      3: zlim, name_all,  1.e-1, 1.e2, 1
    ENDCASE
    
    read:
    
  ENDFOR
  tplot_names, names=names_old
  IF N_ELEMENTS(names_old) GT 0 AND KEYWORD_SET(RECALC) THEN BEGIN
      FOR nn = 0, N_ELEMENTS(names_old)-1 DO BEGIN
          IF (STRMID(names_old(nn), 0, 6) EQ 'ENSPEC') THEN BEGIN
;Jing: change "0" to string(inst,format='(i1.1)' ), so it will also
;work for hia
              IF (STRMID(names_old(nn), STRPOS(names_old(nn), 'IN')+2, 1) EQ string(inst, format = '(i1.1)')) AND $
                (STRPOS(names_old(nn), '_PR') NE -1) AND $  ; JING: ADD THIS LINE TO NOT DELETE OTHER STRING OTHER THAN PRODUCTS SPECTRA
                (STRMID(names_old(nn), STRLEN(names_old(nn))-4, 4) NE '_All') THEN $
                store_data, names_old(nn), /DELETE
          ENDIF
      ENDFOR
  ENDIF  
END 

;+
; PROCEDURE r_en_spec_onefile
;
; PURPOSE: To restore pre-processed energy spectra. Files are saved
;          and restored using the tplot_save/tplot_restore functions
;
; INPUT: sat        -> Satellite number
;        specie     -> 0: H+, 1: He++, 2: He+, 3: O+
;        inst       -> Instrument (0: CODIF, 1: HIA)
;        units_name -> Units for energy spectra
;        angle      -> Angle range to sum over [thetarange, phirange]
;        eff_table  -> 0: GROUND, 1: ONBOARD
;
; OUTPUT: It creates a tplot variable
;
; CREATED BY: C. Mouikis (09/14/01)
;
; LAST MODIFIED: 09/14/01
;
; MODIFICATION HISTORY:
;
;-
PRO r_en_spec_onefile, sat, specie, inst, units_name, angle, eff_table
  
  COMMON get_error, get_err_no, get_err_msg, default_verbose

  ;----------------------------------------------------------------------
  ; Read pre-processed data
  ;----------------------------------------------------------------------
  sp_name = ['H1','He2','He1','O1']

  path = getenv('L4DATA')
  path = path + '/en_spec'
  
  IF angle(1,1) EQ 270.0 THEN BEGIN
    path = getenv('L4DATA')
    path = path + '../L4_en_spec_90_270/en_spec'
  ENDIF
  
  IF angle(1,1) EQ 90.0 THEN BEGIN
    path = getenv('L4DATA')
    path = path + '../L4_en_spec_270_90/en_spec'
  ENDIF
  
  get_timespan, time_interval
  t_s=gettime(time_interval(0)) ; start time in tplot-time
  t_e=gettime(time_interval(1)) ; end time in tplot-time  
  
  t_s_str = time_struct(t_s)    ; start_time tplot time structure
  t_e_str = time_struct(t_e)    ; end_time tplot time structure
  
  mjd_s = julday(t_s_str.month, t_s_str.date, t_s_str.year) ;start julian day
  mjd_e = julday(t_e_str.month, t_e_str.date, t_e_str.year) ; end julian day
  
  no_of_files = (mjd_e - mjd_s) + 1 ; number of days to be loaded
  
  ;Last day is not included if hour=min=sec=0
  IF t_e_str.hour EQ 0 AND t_e_str.min EQ 0 AND t_e_str.sec EQ 0 THEN $
    no_of_files = no_of_files - 1
  
  ;--------------------------------------------------------------------
  ; Read all 1 day files that correspond to requested time interval
  ;--------------------------------------------------------------------
  ffc = 0                        ; Files-found counter
  FOR nd = 0 , no_of_files-1 DO BEGIN ; Loop trough all days
    
    caldat, mjd_s + nd, month, day, year ; find caledar date
    

    tplot_var_dflux_name = $ ; units: DIFF FLUX (units files are saved in)
      'ENSPEC' + $
      '_SC' + strcompress(sat,/remove_all) + $
      '_IN' + STRCOMPRESS(inst,/REMOVE_ALL) + '*' + $
      '_UN' + STRUPCASE(strcompress('DIFF FLUX',/remove_all)) + $
      '_SP' + strcompress(specie,/remove_all) + $
      '_ET' + strcompress(eff_table,/remove_all) + $
      '_All'
    
    tplot_var_name = $ ; units selected
      'ENSPEC' + $
      '_SC' + strcompress(sat,/remove_all) + $
      '_IN' + STRCOMPRESS(inst,/REMOVE_ALL) + '*' + $
      '_UN' + STRUPCASE(strcompress(units_name,/remove_all)) + $
      '_SP' + strcompress(specie,/remove_all) + $
      '_ET' + strcompress(eff_table,/remove_all) + $
      '_All'
    
    
    IF (angle(0,1) EQ 90 AND angle(1,1) EQ 270) OR $
      (angle(0,1) EQ 270 AND angle(1,1) EQ 90) THEN BEGIN
      
      tplot_var_dflux_name = $ ; units: DIFF FLUX (units files are saved in)
        'ENSPEC' + $
        '_SC'  + strcompress(sat,/remove_all) + $
        '_IN'  + STRCOMPRESS(inst,/REMOVE_ALL) + $
        '_PHI' + STRCOMPRESS(FIX(angle(0,1)), /REMOVE_ALL) + $
        '_' + STRCOMPRESS(FIX(angle(1,1)), /REMOVE_ALL) + $
        '_UN'  + STRUPCASE(strcompress('DIFF FLUX',/remove_all)) + $
        '_SP'  + strcompress(specie,/remove_all) + $
        '_ET'  + strcompress(eff_table,/remove_all) + $
        '_All'
      
      tplot_var_name = $        ; units selected
        'ENSPEC' + $
        '_SC'  + strcompress(sat,/remove_all) + $
        '_IN'  + STRCOMPRESS(inst,/REMOVE_ALL) + $
        '_PHI' + STRCOMPRESS(FIX(angle(0,1)), /REMOVE_ALL) + $
        '_' + STRCOMPRESS(FIX(angle(1,1)), /REMOVE_ALL) + $
        '_UN'  + STRUPCASE(strcompress(units_name,/remove_all)) + $
        '_SP'  + strcompress(specie,/remove_all) + $
        '_ET'  + strcompress(eff_table,/remove_all) + $
        '_All'
      
    ENDIF
    
    filename = 'COD_EN_SPEC_SC' + string(sat,format='(i1.1)') + $
      '_' + sp_name(specie) + '_' + $
      string(year, month, day, format='(i4.4,i2.2,i2.2)') + $
      '.tplot*'
    
    ; Set the directory where the pre-processed data exist
    IF sat EQ 1 THEN sc_dir = '/CLUSTER1/' + sp_name(specie) + '/'
    IF sat EQ 2 THEN sc_dir = '/CLUSTER2/' + sp_name(specie) + '/'
    IF sat EQ 3 THEN sc_dir = '/CLUSTER3/' + sp_name(specie) + '/'
    IF sat EQ 4 THEN sc_dir = '/CLUSTER4/' + sp_name(specie) + '/'
    
    ff = findfile(path + sc_dir + filename, COUNT=fc)

    IF fc GT 0 THEN BEGIN
      
      file_path = ff(0)
      zipflag=0
      IF STRMID(ff(0),STRLEN(ff(0))-2,2) EQ 'gz' THEN BEGIN
        fln_from_path, ff(0), fln
        unziped_filename = strmid(fln, 0, strlen(fln)-3)
        unziped_filename = getenv('CCAT_TEMP') + '/' + unziped_filename
        PRINT, 'gzip -dc ' + ff(0) + ' > ' + unziped_filename
        SPAWN, 'gzip -dc ' + ff(0) + ' > ' + unziped_filename
        file_path = unziped_filename
        zipflag=1
      ENDIF

      ffc = ffc + 1
      IF ffc GT 1 THEN BEGIN
        ;get old data
        tplot_names, tplot_var_name, names=names
        tplot_var_name = names(0)
        get_data, tplot_var_name, $
          data=data_old, dlim=dlim, lim=lim
        store_data, tplot_var_name, /delete
        
        ;get new data
        tplot_restore, filename = file_path
        IF units_name EQ 'DIFF FLUX' THEN BEGIN ; if units DIFF FLUX
          tplot_names, tplot_var_name, names=names
          tplot_var_name = names(0)
          get_data, tplot_var_name, data=data_new
          store_data, tplot_var_name, /delete
        ENDIF ELSE BEGIN ; if units EFLUX
          tplot_names, tplot_var_dflux_name, names=names
          tplot_var_dflux_name = names(0)
          get_data, tplot_var_dflux_name, data=data_new
          ; transform DIFF FLUX in EFLUX and save as tplot variable
          FOR jj = 0, N_ELEMENTS(data_new.v(0,*))-1 DO BEGIN
            data_new.y(*,jj) = data_new.y(*,jj) * data_new.v(0,jj)
          ENDFOR
          store_data, tplot_var_dflux_name, /delete
        ENDELSE

        ;combine the two
        en_spec_combine_days, $
          tplot_var_name, units_name, data_old, data_new, dlim, lim
      ENDIF ELSE BEGIN ; restore the first file
        
        tplot_restore, filename = file_path
        IF units_name EQ 'EFLUX' THEN BEGIN ; if units EFLUX
          tplot_names, tplot_var_dflux_name, names=names
          tplot_var_dflux_name = names(0)
          tplot_var_name = $
            STRMID(tplot_var_dflux_name, 0, STRPOS(tplot_var_dflux_name, 'UN')+2) + $
            'EFLUX' + $
            STRMID(tplot_var_dflux_name, STRPOS(tplot_var_dflux_name, '_SP'), $
                   STRLEN(tplot_var_dflux_name))
          get_data, tplot_var_dflux_name, data=data_new, dlim=dlim, lim=lim
          ; transform DIFF FLUX in EFLUX and save as tplot variable
          FOR jj = 0, N_ELEMENTS(data_new.v(0,*))-1 DO BEGIN
            data_new.y(*,jj) = data_new.y(*,jj) * data_new.v(0,jj)
          ENDFOR
          store_data, tplot_var_name, data=data_new, dlim=dlim, lim=lim
          store_data, tplot_var_dflux_name, /delete
        ENDIF

      ENDELSE
      
      
      IF zipflag EQ 1 THEN BEGIN
        PRINT, '/bin/rm -f ' + file_path
        spawn, '/bin/rm -f ' + file_path
      ENDIF
      
    ENDIF
  ENDFOR
  
  IF ffc eq 0 THEN BEGIN
    get_err_no = 1
    PRINT, 'No Processed files found'
  ENDIF

END


;+
;-
PRO en_spec_combine_days, tplot_var_name, $
                          units_name, $
                          data_old, $
                          data_new, $
                          dlim, lim
  
  ;Initialize the data arrays with the data_old values
  xall = data_old.x
  yall = data_old.y
  vall = data_old.v
  
  nold = n_elements(data_old.y(0,*))  ; number of energy bins for old
  nnew = n_elements(data_new.y(0,*)) ; number of energy bins for new
  
  IF nold EQ nnew THEN BEGIN    ; if both have the same # of en. bins
    xall = [xall, data_new.x]
    yall = [yall, data_new.y]
    vall = [vall, data_new.v]
  ENDIF ELSE BEGIN
    IF nold EQ 31 AND  nnew EQ 16  THEN BEGIN ; old:31 & new:16
      
      yall_new = dblarr(n_elements(data_new.x),31)
      vall_new = vall
      FOR i2 = 0, 14 DO BEGIN
        yall_new(*,i2*2)   = data_new.y(*,i2)
        yall_new(*,i2*2+1) = data_new.y(*,i2)
      ENDFOR
      yall_new(*,30) = data_new.y(*,15)
      
      xall = [xall, data_new.x]
      yall = [yall, yall_new]
      vall = [vall, vall_new]
      
    ENDIF
    
    IF nold EQ 16 AND  nnew EQ 31  THEN BEGIN ; old:16 & new:31
      
      yall_old = yall
      vall_old = vall
      
      energy_new = data_new.v(0,*)
      
      yall = dblarr(n_elements(xall),31)
      vall = dblarr(n_elements(xall),31)
      FOR i2 = 0, 14 DO BEGIN
        yall(*,i2*2)   = yall_old(*,i2)
        yall(*,i2*2+1) = yall_old(*,i2)
        vall(*,i2*2)   = energy_new(i2*2)
        vall(*,i2*2+1) = energy_new(i2*2+1)
      ENDFOR
      yall(*,30) = yall_old(*,i2)
      vall(*,30) = energy_new(30)
      
      xall = [xall, data_new.x]
      yall = [yall, data_new.y]
      vall = [vall, data_new.v]
      
    ENDIF
  ENDELSE
  
  sort_ind = SORT(xall)
    
  xall_new = xall(sort_ind)
  yall_new = yall(sort_ind,*)
  vall_new = vall(sort_ind,*)
  
  datastr = {x:xall_new,y:yall_new,v:vall_new}
  
  store_data, tplot_var_name, data=datastr, dlim=dlim, lim=lim
END
