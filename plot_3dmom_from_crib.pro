;+
;PROCEDURE: plot_3dmom_from_crib
;PURPOSE: 
;  It is a crib sheet for plotting cluster cis 3D moments using
;  the TPLOT package
;
;INPUT
;PARAMETERS:   sat:        Satellite number
;              inst:       Instrument (0: CODIF, 1: HIA)
;              prod:       Product number
;              time:       Date/time in tplot format
;              angle:      Angle range to sum over [thetarange, phirange]
;
;                          THETA_88: [-78.75, -56.25, -33.75, -11.25,
;                                     11.25, 33.75, 56.25, 78.75] 
;
;                          PHI_88:   [11.25, 22.50, 33.75, 45.00, 56.25,
;                                     67.50, 78.75, 101.25, 112.50, 123.75,
;                                     135.00, 146.25, 157.50, 168.75, 191.25,
;                                     202.50, 213.75, 225.00, 236.25, 247.50,
;                                     258.75, 281.25, 292.50, 303.75, 315.00,
;                                     326.25, 337.50, 348.75]
;
;                          THETA_24: [-78.75, -56.25, -33.75, -11.25,
;                                     11.25, 33.75, 56.25, 78.75] 
;
;                          PHI_24:   [45.0, 90.0, 135.0, 225.0, 270.0, 315.0]
;
;KEYWORDS: INST_COORD -> allows moments to be calculated in the
;                        instrument coordinate system
;          RECALC     -> forces the recalculation of the moments
;                        even if a pre-processed file does exist
;          NEW_NAME       -> sets a variable new_name different than the default
;
;CREATED BY: C. Mouikis
;
;LAST MODIFICATION: 06/05/12
;
;MODIFICATION HISTORY:
;        05/22/01: Temperature and Pressure are included
;        07/17/01: The tplot variable names retain the efficiency
;                  table information
;        09/19/01: The keyword INST_COORD is added.
;        10/26/01: If pre-processed files exist, they are read otherwise
;                  the 3d moments are calculated from the L1 files
;        10/26/01: The RECALC keyword is introduced allowing the
;                  possibility to recalculate the moments even if a
;                  pre-processed file does exist
;        10/26/01: If tplot save file is compressed it uncompresses it
;                  in the directory declared by the CCAT_TEMP
;                  environment variable
;        11/16/01: Keyword new_name is added.
;        11/16/01: variable moment changed to moments because moment
;                  is an idl function
;        12/05/01: default units for EFLUX changed to ev/cm^2-s
;        12/05/01: If pre-processed files are read the eflux units are
;                  converted to ev/cm^2-s
;        02/21/03: Added: keyword INCRATES
;        06/15/03: added: keyword SPILL
;        05/06/05: added: keyword PREPROCESSED
;        07/21/10: added: keyword HE1_CLEAN
;        03/07/11: Added: keyword ERROR_BAR
;        10/01/11: Added: keyword FRS
;        06/05/12: Added: Keyword RAPID. The tplot variable names
;                         created had to be extended to describe the 
;                         longer energy range
;       2019/11/19 Spacecraft potential keyward is passed through
;-
PRO plot_3dmom_from_crib, sat,                   $
                          specie,                $
                          inst,                  $
                          moments,               $
                          angle,                 $
                          energy,                $
                          eff_table,             $
                          scp = scp, $
                          frs=frs,               $
                          inst_coord=inst_coord, $
                          recalc=recalc,         $
                          error_bar=error_bar,   $
                          spill=spill,           $
                          bkg=bkg,               $
                          rapid=rapid,           $
                          he1_clean=he1_clean,   $
                          no_phi_cor=no_phi_cor, $
                          new_name=new_name,     $
                          old_eff=old_eff,       $
                          incrates=incrates,     $
                          interp_rates=interp_rates, $
                          preprocessed=preprocessed, $
                          bins = bins_input, $                  ;jing
                          use_bins = use_bins, $ ;jing
                          e_r = e_r,  $ ; jing
                          a_r = a_r, $ ;jing
                          sum_up = sum_up,$ ;jing
                          diffflux_threshold=diffflux_threshold ;jing

  COMMON get_error, get_err_no, get_err_msg, default_verbose

  specie_str = ['H!U+!N','He!U++!N','He!U+!N','O!U+!N']

  FOR ii=0, n_elements(moments)-1 DO BEGIN

    CASE moments(ii) OF
      'D': moments(ii)='DENSITY'
      'V': moments(ii)='VELOCITY'
      'P': moments(ii)='PRESSURE'
      'T': moments(ii)='TEMPERATURE'
      'J': moments(ii)='JFLUX'
      'E': moments(ii)='EFLUX' ; energy flux, in eV/(cm^2-s)
      'A': moments(ii)='ALL'  
    ENDCASE
  ENDFOR
  
  ; Check if sat, specie and moments arrays have the same number of elements
  IF n_elements(sat) NE n_elements(specie) OR $
    n_elements(sat) NE n_elements(moments)  OR $
    n_elements(specie) NE n_elements(moments)  THEN BEGIN
    print, 'The sat array, the specie array and the moments'
    print, 'array MUST have the same number of elements.'
    stop
  ENDIF
  
  ; Loop over all sat/specie/moments combinations
  FOR ii=0,n_elements(specie)-1 DO BEGIN
    
    get_err_no = 0 ; reset error indicator

    ; Read pre-processed file
    IF NOT KEYWORD_SET(RECALC) THEN BEGIN
      IF eff_table EQ 0 AND $  ; this list has to be completed
        inst EQ 0       THEN BEGIN
        r_3dmom_onefile, sat(ii), specie(ii)
      ENDIF

      IF get_err_no EQ 0 THEN BEGIN
;        vn = 'TDMOM_EN00040_40000_SC' + $
;          STRING(sat(ii),FORMAT='(i1.1)') + $
;          '_MTEFLUX_SP' + $
;          STRING(specie(ii),FORMAT='(i1.1)') + $
;          '_ET0_All'
;        get_data, vn, data=data, dlim=dlim, lim=lim
;        data.y=data.y/1.0622e-12
;        lim.ytitle = $
;          'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
;          specie_str(specie(ii))+' JE (eV cm!U-2!N s!U-1!N)' + '!C!C' + $
;          STRCOMPRESS(STRING(energy(0), format='(i5)'),/REMOVE_ALL) + $
;          '-' + $
;          STRCOMPRESS(STRING(energy(1), format='(i5)'),/REMOVE_ALL) + $
;          ' (eV)'
;
;        store_data, vn, data=data, dlim=dlim, lim=lim
        GOTO, read
      ENDIF ELSE BEGIN
        IF KEYWORD_SET(PREPROCESSED) THEN BEGIN
          get_err_no = 2 ; No pre-processed file found
          RETURN
        ENDIF
      ENDELSE
    ENDIF

    ; Products that correspond to a specie
    CASE specie(ii) OF
      0: prod = [12, 13]
      1: prod = [15, 16]
      2: prod = [17, 18, 46, 48]
      3: prod = [17, 18, 47, 49]
      ELSE: BEGIN
        print, 'Specie numbers must be 0 for H+, 1 for He++,'
        print, '2 for He+ or 3 for O+'
        stop
      END
   ENDCASE

;Jing add prod for hia
    IF inst EQ 1 THEN  prod = [6, 15, 23] ;17

    ;Loop over all products that correspond to a specie
    FOR jj=0,n_elements(prod)-1 DO BEGIN
      
      IF NOT KEYWORD_SET(NEW_NAME) THEN BEGIN
        IF KEYWORD_SET(INCRATES) THEN BEGIN
          name= 'TDMOM_INCRATES_' + $
            '_EN' + STRCOMPRESS(STRING(energy(0), $
                                       format='(i7.7)'),/REMOVE_ALL) + $
            '_'   + STRCOMPRESS(STRING(energy(1), $
                                       format='(i7.7)'),/REMOVE_ALL) + $
            '_SC' + strcompress(sat(ii),/remove_all) + $
            '_MT' + strmid(moments(ii),0,1) + $
            '_PR' + strcompress(string(prod(jj)),/remove_all) + $
            '_SP' + strcompress(specie(ii),/remove_all) + $
            '_ET' + strcompress(eff_table,/remove_all) ;name of TPLOT structure
        ENDIF ELSE BEGIN
          name= 'TDMOM' + $
            '_EN' + STRCOMPRESS(STRING(energy(0), $
                                       format='(i7.7)'),/REMOVE_ALL) + $
            '_'   + STRCOMPRESS(STRING(energy(1), $
                                       format='(i7.7)'),/REMOVE_ALL) + $
            '_SC' + strcompress(sat(ii),/remove_all) + $
            '_MT' + strmid(moments(ii),0,1) + $
            '_PR' + strcompress(string(prod(jj)),/remove_all) + $
            '_SP' + strcompress(specie(ii),/remove_all) + $
            '_ET' + strcompress(eff_table,/remove_all) ;name of TPLOT structure
        ENDELSE
      ENDIF ELSE BEGIN
        name = new_name
      ENDELSE
;Jing: If the keyword use_bins is set, then set the angle to 0, bins
;to bins_input. If not, then set the bin to 0. 
   IF KEYWORD_SET(use_bins) THEN BEGIN 
            angle = 0
            bins = bins_input
        ENDIF  ELSE bin = 0     

      ; Calculate moments
      codif_moments, $
        sat(ii), $
        inst, $
        prod(jj), $
        moments(ii), $
        scp = scp, $
        frs=frs, $
        specie=specie(ii), $
        eff_table, $
        name=name, $
        angle=angle, $
        energy=energy, $
        inst_coord=inst_coord, $
        no_phi_cor=no_phi_cor, $
        old_eff=old_eff, $
        incrates=incrates, $
        interp_rates=interp_rates, $
        spill=spill, $
        bkg=bkg, $
        rapid=rapid, $
        he1_clean=he1_clean, $
        bins = bins, $         ;jing
        e_r = e_r, $          ;Jing
        a_r = a_r, $          ;Jing
        sum_up = sum_up,$       ;jing
        diffflux_threshold=diffflux_threshold

      IF get_err_no GT 0 THEN GOTO, next

      enrange=STRCOMPRESS(STRING(energy(0), format='(i7.7)'),/REMOVE_ALL) + $
        '_'   + STRCOMPRESS(STRING(energy(1), format='(i7.7)'),/REMOVE_ALL)
      
      IF moments(ii) EQ 'ALL' THEN BEGIN
        
        tplot_names, 'TDMOM*', names=names_all
        FOR im = 0, N_ELEMENTS(names_all)-1 DO BEGIN
          mom = STRMID(names_all(im), STRPOS(names_all(im), '_MT')+3, 1)
          CASE mom OF
            'D': mom='DENSITY'
            'V': mom='VELOCITY'
            'P': mom='PRESSURE'
            'T': mom='TEMPERATURE'
            'J': mom='JFLUX'
            'E': mom='EFLUX'
          ENDCASE
          combine_3dmoments, sat(ii), specie(ii), $
            mom, name_all, enrange, $
            strcompress(string(eff_table), /remove_all), $
            incrates=incrates, error_bar=error_bar
        
          ; Set plot attributes
          CASE mom OF
            'DENSITY' : BEGIN
              yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
                specie_str(specie(ii))+' n (cm!U-3!N)'+ '!C!C' + $
                STRCOMPRESS(STRING(energy(0), format='(i7)'),/REMOVE_ALL) + '-' + $
                STRCOMPRESS(STRING(energy(1), format='(i7)'),/REMOVE_ALL) + ' (eV)'
              ylim,    name_all,  0, 0, 1
            END
            'VELOCITY': BEGIN
              yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
                specie_str(specie(ii))+' V (km s!U-1!N)'
              ylim,    name_all,  0, 0, 0
            END
            'TEMPERATURE': BEGIN
              yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
                specie_str(specie(ii))+' T (eV)'
              ylim,    name_all,  0, 0, 1
            END
            'PRESSURE': BEGIN
              yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
                specie_str(specie(ii))+' P (nPa)'
              ylim,    name_all,  0, 0, 1
            END
            'JFLUX': BEGIN
              yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
                specie_str(specie(ii))+' J (cm!U-2!N s!U-1!N)'
              ylim,    name_all,  0, 0, 0
            END
            'EFLUX': BEGIN
;          yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
;            specie_str(specie(ii))+' JE (mW m!U-2!N)' + '!C!C' + $
;            STRCOMPRESS(STRING(energy(0), format='(i7)'),/REMOVE_ALL) $
;            + '-' + $
;            STRCOMPRESS(STRING(energy(1), format='(i7)'),/REMOVE_ALL) $
;            + ' (eV)'
          
              yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
                specie_str(specie(ii))+' JE (eV cm!U-2!N s!U-1!N)' + '!C!C' + $
                STRCOMPRESS(STRING(energy(0), format='(i7)'),/REMOVE_ALL) + $
                '-' + $
                STRCOMPRESS(STRING(energy(1), format='(i7)'),/REMOVE_ALL) + $
                ' (eV)'
              
              ylim,    name_all,  0, 0, 0
            END
          ENDCASE       
          options, name_all, 'ytitle', yt
        ENDFOR
        
      ENDIF ELSE BEGIN
        
        combine_3dmoments, sat(ii), specie(ii), $
          moments(ii), name_all, enrange, $
          strcompress(string(eff_table), /remove_all), $
          incrates=incrates, error_bar=error_bar
        
        ; Set plot attributes
        CASE moments(ii) OF
          'DENSITY' : BEGIN
            yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
              specie_str(specie(ii))+' n (cm!U-3!N)'+ '!C!C' + $
              STRCOMPRESS(STRING(energy(0), format='(i7)'),/REMOVE_ALL) + '-' + $
              STRCOMPRESS(STRING(energy(1), format='(i7)'),/REMOVE_ALL) + ' (eV)'
            ylim,    name_all,  0, 0, 1
          END
          'VELOCITY': BEGIN
            yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
              specie_str(specie(ii))+' V (km s!U-1!N)'
            ylim,    name_all,  0, 0, 0
          END
          'TEMPERATURE': BEGIN
            yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
              specie_str(specie(ii))+' T (eV)'
            ylim,    name_all,  0, 0, 1
          END
          'PRESSURE': BEGIN
            yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
              specie_str(specie(ii))+' P (nPa)'
            ylim,    name_all,  0, 0, 1
          END
          'JFLUX': BEGIN
            yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
              specie_str(specie(ii))+' J (cm!U-2!N s!U-1!N)'
            ylim,    name_all,  0, 0, 0
          END
          'EFLUX': BEGIN
;          yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
;            specie_str(specie(ii))+' JE (mW m!U-2!N)' + '!C!C' + $
;            STRCOMPRESS(STRING(energy(0), format='(i7)'),/REMOVE_ALL) $
;            + '-' + $
;            STRCOMPRESS(STRING(energy(1), format='(i7)'),/REMOVE_ALL) $
;            + ' (eV)'
          
            yt = 'SC' + string(sat(ii), format='(i1.1)') + '!C!C' + $
              specie_str(specie(ii))+' JE (eV cm!U-2!N s!U-1!N)' + '!C!C' + $
              STRCOMPRESS(STRING(energy(0), format='(i7)'),/REMOVE_ALL) + $
              '-' + $
              STRCOMPRESS(STRING(energy(1), format='(i7)'),/REMOVE_ALL) + $
              ' (eV)'
            
            ylim,    name_all,  0, 0, 0
          END
        ENDCASE      
        options, name_all, 'ytitle', yt
        
      ENDELSE
      
      options, name_all, 'ytitle', yt
      
      next:
    ENDFOR
    
    read:

  ENDFOR    

  ; Delete all product related structures
  tplot_names, names=names_old
  IF n_elements(names_old) GT 0 THEN BEGIN
    FOR nn=0, n_elements(names_old)-1 DO BEGIN
      IF (strmid(names_old(nn),strlen(names_old(nn))-3,2) EQ 'ET') THEN $
          store_data, names_old(nn), /DELETE
    ENDFOR
  ENDIF

  ; If keyword ERROR_BAR is not set then
  ; delete the error variable from the
  ; data structure

  IF NOT KEYWORD_SET(ERROR_BAR) THEN BEGIN
      tplot_names, names=names_old
      
      IF n_elements(names_old) GT 0 THEN BEGIN
          FOR nn=0, n_elements(names_old)-1 DO BEGIN
              IF (strmid(names_old(nn),0,5) EQ 'TDMOM') AND $
                (strpos(names_old(nn),'DENSITY') NE -1) THEN BEGIN
                  
                  get_data, names_old(nn), data=data, lim=lim, dlim=dlim
                  store_data, names_old(nn), data={x:data.x, y:data.y}, dlim=dlim, lim=lim

              ENDIF
          ENDFOR
          
      ENDIF

  ENDIF


  ; Store the vector components in different structures
  tplot_names, names=names_old

  IF n_elements(names_old) GT 0 THEN BEGIN
    FOR nn=0, n_elements(names_old)-1 DO BEGIN
      IF (strmid(names_old(nn),0,5) EQ 'TDMOM') AND $
        (strmid(names_old(nn),strlen(names_old(nn))-2,2) NE '_X') AND $
        (strmid(names_old(nn),strlen(names_old(nn))-2,2) NE '_Y') AND $
        (strmid(names_old(nn),strlen(names_old(nn))-2,2) NE '_Z') AND $
        (strmid(names_old(nn),strlen(names_old(nn))-2,2) NE '_T') THEN BEGIN
        
        IF (strpos(names_old(nn),'VELOCITY') NE -1) THEN BEGIN
          sat = strmid(names_old(nn),strpos(names_old(nn),'_SC')+3,1)*1
          specie = strmid(names_old(nn),strpos(names_old(nn),'_SP')+3,1)*1
          plot_vector, names_old(nn), 'V', sat, /polar, $
            specie=specie_str(specie), units='km s!U-1!N', /erange 
        ENDIF
        
        IF (strpos(names_old(nn),'PRESSURE') NE -1) THEN BEGIN
          sat = strmid(names_old(nn),strpos(names_old(nn),'_SC')+3,1)*1
          specie = strmid(names_old(nn),strpos(names_old(nn),'_SP')+3,1)*1
          plot_vector, names_old(nn), 'P', sat, $
            specie=specie_str(specie), units='nPa', /ylog, /erange
        ENDIF
        
        IF (strpos(names_old(nn),'TEMPERATURE') NE -1) THEN BEGIN
          sat = strmid(names_old(nn),strpos(names_old(nn),'_SC')+3,1)*1
          specie = strmid(names_old(nn),strpos(names_old(nn),'_SP')+3,1)*1
          plot_vector, names_old(nn), 'T', sat, $
            specie=specie_str(specie), units='eV', /ylog, /erange
        ENDIF
        
        IF (strpos(names_old(nn),'JFLUX') NE -1) THEN BEGIN
          sat = strmid(names_old(nn),strpos(names_old(nn),'_SC')+3,1)*1
          specie = strmid(names_old(nn),strpos(names_old(nn),'_SP')+3,1)*1
          plot_vector, names_old(nn), 'J', sat, $
            specie=specie_str(specie), units='cm!U-2!N s!U-1!N', /erange
        ENDIF
        
        IF (strpos(names_old(nn),'EFLUX') NE -1) THEN BEGIN
          sat = strmid(names_old(nn),strpos(names_old(nn),'_SC')+3,1)*1
          specie = strmid(names_old(nn),strpos(names_old(nn),'_SP')+3,1)*1
          plot_vector, names_old(nn), 'JE', sat, $
            specie=specie_str(specie), units='eV cm!U-2!N s!U-1!N', /erange
;          plot_vector, names_old(nn), 'JE', sat, $
;            specie=specie_str(specie), units='mW m!U-2!N', /erange
        ENDIF

      ENDIF
    ENDFOR
  ENDIF
  
END

PRO r_3dmom_onefile, sat, specie
  
  
  COMMON get_error, get_err_no, get_err_msg, default_verbose
  
  ;----------------------------------------------------------------------
  ; Read pre-processed data
  ;----------------------------------------------------------------------
  sp_name = ['H1','He2','He1','O1']
  
  path = getenv('L4DATA')
  path = path + '/3d_mom'
  
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
    
    filename = 'COD_3DMOM_SC' + string(sat,format='(i1.1)') + $
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
        SPAWN, 'gzip -dc ' + ff(0) + ' > ' + unziped_filename
        file_path = unziped_filename
        zipflag=1
      ENDIF
      
;      tplot_restore, filename = file_path, name=tvar_name, /enqname
      
;      findvar = WHERE(STRPOS(tvar_name, moments) GE 0)

      ffc = ffc + 1
      IF ffc GT 1 THEN BEGIN
        
        tplot_restore, filename = file_path, /APPEND
        
      ENDIF ELSE BEGIN          ; restore the first file
        
        tplot_restore, filename = file_path

      ENDELSE
      
      
      IF zipflag EQ 1 THEN BEGIN
        spawn, '/bin/rm -f ' + file_path
      ENDIF
      
      
    ENDIF
  ENDFOR

  IF ffc eq 0 THEN BEGIN
    get_err_no = 1
    PRINT, 'No Processed files found'
  ENDIF



END
