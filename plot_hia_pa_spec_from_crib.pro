;+
;PROCEDURE: plot_pa_spec_from_crib
;PURPOSE: 
;  It is a crib sheet for plotting cluster cis energy spectra using
;  the TPLOT package
;
;INPUT
;PARAMETERS:   sat:        Satellite number
;              inst:       Instrument (0: CODIF, 1: HIA)
;              time:       Date/time in tplot format
;              timespan:   Time span
;              units_name: Units for energy spectra
;                          'Counts', 'NCOUNTS', 'RATE', 'NRATE', 
;                           'DIFF FLUX', 'EFLUX'
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
;CREATED BY: C Mouikis
;
;LAST MODIFIED: 05/22/01
;
;MODIFICATION HISTORY:
;   05/22/01 - PABIN keyword included
;   05/22/01 - EDITA keyword included
;   05/22/01 - IC keyword included
;   08/01/01 - RECALC keyword included
;   02/05/02 - The keyword ALL_ENERGY_BINS is introduced
;   03/13/02 - OLD_EFF keyword added
;   08/08/02 - BKG keyword added
;-

PRO plot_hia_pa_spec_from_crib, sat, specie, inst, units_name, $
                            energy, eff_table, $
                            PABIN=PABIN, $
                            IC=IC, $
                            CNES=CNES, $
                            FORNACON=FORNACON, $
                            EDITA=EDITA, $
                            RECALC=RECALC, $
                            ALL_ENERGY_BINS=ALL_ENERGY_BINS, $
                            OLD_EFF=OLD_EFF, $
                            BKG=BKG, $
                            COMBINE=COMBINE
  
  COMMON get_error, get_err_no, get_err_msg, default_verbose
  
  specie_str = ['H!U+!N','He!U++!N','He!U+!N','O!U+!N']
  
  ; Check if sat and specie arrays have the same number of elements
  IF n_elements(sat) NE n_elements(specie) THEN BEGIN
    print, 'The sat array and the specie array MUST have'
    print, 'the same number of elements.'
    stop
  ENDIF
  
  ; Loop over all sat/specie combinations
  FOR ii=0,n_elements(specie)-1 DO BEGIN
    
    get_err_no = 0 ; reset error indicator
    prod = [6, 15, 23] 
    ;Loop over all products that correspond to a specie
    FOR jj=0,n_elements(prod)-1 DO BEGIN

      ; Read pre-processed file
      IF NOT KEYWORD_SET(RECALC) THEN BEGIN
        IF eff_table EQ 0 THEN BEGIN
          IF units_name EQ 'DIFF FLUX' THEN BEGIN
            r_pa_spec_onefile, $
              sat(ii), specie(ii), inst, units_name, $
              energy, prod(jj), eff_table
          ENDIF ELSE get_err_no = 1
        ENDIF ELSE get_err_no = 1 
     
        get_err_no = 0
        IF get_err_no EQ 0 THEN GOTO, read
      ENDIF
      
      name = 'PASPEC' + $
        '_EN' + STRCOMPRESS(STRING(energy(0), format='(i5.5)'),/REMOVE_ALL) + $
        '_'   + STRCOMPRESS(STRING(energy(1), format='(i5.5)'),/REMOVE_ALL) + $
        '_SC' + STRCOMPRESS(sat(ii),/REMOVE_ALL) + $
        '_UN' + STRUPCASE(strcompress(units_name,/REMOVE_ALL)) + $
        '_PR' + STRCOMPRESS(string(prod(jj)),/REMOVE_ALL) + $
        '_SP' + STRCOMPRESS(specie(ii),/REMOVE_ALL) + $
        '_ET' + STRCOMPRESS(eff_table,/REMOVE_ALL) ;name of TPLOT structure
      
      ;------------------------------------------------------------------
      ; dat structure is needed for the 
      ; GSE -> CODIF coordinate transformation
      ;------------------------------------------------------------------
      dat = call_function('get_cis_hia_data',prod(jj), sat)
      
      IF get_err_no EQ 0 THEN BEGIN ;check if data were found for time interval
        found = 1
        
        get_theta_phi, sat, dat, mag_theta, mag_phi, inst, $
          FORNACON=FORNACON, $
          EDITA=EDITA, $
          IC=IC, CNES=CNES, /ALL            ; get magnetic field data
        ;all_energy_bins = 1
        get_hia_pa_spec,   $
          sat(ii),           $
          inst,              $
          prod(jj),          $
          mag_theta,         $
          mag_phi,           $
          specie=specie(ii), $
          eff_table,         $
          angle=angle,       $
          energy=energy,     $
          units=units_name,  $
          name=name,         $
          pabin=pabin,       $
          all_energy_bins=all_energy_bins, $
          old_eff=old_eff, $
          bkg=bkg, $
          combine=combine
      ENDIF
;stop
      read:

      IF get_err_no GT 0 THEN GOTO, next
      
      IF NOT KEYWORD_SET(ALL_ENERGY_BINS) THEN BEGIN
      
        ; Combine the different product of the same specie
        enrange=STRCOMPRESS(STRING(energy(0), format='(i5.5)'),/REMOVE_ALL) + $
          '_'   + STRCOMPRESS(STRING(energy(1), format='(i5.5)'),/REMOVE_ALL)
        combine_pa_spec, sat(ii), specie(ii), $
          strcompress(STRUPCASE(units_name), /remove_all), name_all, $
          enrange, eff_table
        
        ; Set plot attributes
        CASE STRUPCASE(units_name) OF
          'COUNTS': uname = 'COUNTS'
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
          specie_str(specie(ii)) + '!C!C' + $
          STRCOMPRESS(STRING(energy(0), format='(i5)'),/REMOVE_ALL) + '-' + $
          STRCOMPRESS(STRING(energy(1), format='(i5)'),/REMOVE_ALL) + ' (eV)'
        options, name_all, 'ztitle', uname
        IF KEYWORD_SET(COMBINE) THEN BEGIN
          ylim,    name_all,  0., 180., 0
        ENDIF ELSE BEGIN
          ylim, name_all, 0., 360., 0
        ENDELSE
        
        next:
      ENDIF
    ENDFOR
  ENDFOR

  IF NOT KEYWORD_SET(ALL_ENERGY_BINS) THEN BEGIN
    tplot_names, names=names_old
    IF n_elements(names_old) GT 0 THEN BEGIN
      FOR nn=0, n_elements(names_old)-1 DO BEGIN
        IF ((strmid(names_old(nn),0,6)) EQ 'PASPEC') AND $
          (strmid(names_old(nn),strlen(names_old(nn))-3,2) EQ 'ET') THEN $
          store_data, names_old(nn), /DELETE
      ENDFOR
    ENDIF
  ENDIF

END

;+
; PROCEDURE r_pa_spec_onefile
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
PRO r_pa_spec_onefile, sat, specie, inst, units_name, energy, prod, eff_table

  ;----------------------------------------------------------------------
  ; Read pre-processed data
  ;----------------------------------------------------------------------
  sp_name = ['H1','He2','He1','O1']
  
  path = getenv('L1DATA')
  path = path + '../L4/pa_spec'
  
  get_timespan, time_interval
  t_s=gettime(time_interval(0)) ; start time in tplot-time
  t_e=gettime(time_interval(1)) ; end time in tplot-time  
  
  t_s_str = time_struct(t_s)    ; start_time tplot time structure
  t_e_str = time_struct(t_e)    ; end_time tplot time structure
  
  mjd_s = julday(t_s_str.month, t_s_str.date, t_s_str.year) ;start julian day
  mjd_e = julday(t_e_str.month, t_e_str.date, t_e_str.year) ; end julian day
  
  no_of_files = (mjd_e - mjd_s) + 1 ; number of days to be loaded

  ;--------------------------------------------------------------------
  ; Read all 1 day files that correspond to requested time interval
  ;--------------------------------------------------------------------
  ffc = 0                        ; Files-found counter
  FOR nd = 0 , no_of_files-1 DO BEGIN ; Loop trough all days
    
    caldat, mjd_s + nd, month, day, year ; find caledar date

    tplot_var_name = $
      'PASPEC' + $
      '_SC' + strcompress(sat,/remove_all) + $
      '_UN' + STRUPCASE(strcompress(units_name,/remove_all)) + $
      '_PR' + STRCOMPRESS(string(prod),/REMOVE_ALL) + $
      '_SP' + strcompress(specie,/remove_all) + $
      '_ET' + STRCOMPRESS(eff_table,/REMOVE_ALL) ;name of TPLOT structure

    filename = 'COD_PA_SPEC_SC' + string(sat,format='(i1.1)') + $
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

    ENDIF
    
    store_name = $
      'PASPEC' + $
      '_EN' + STRCOMPRESS(STRING(energy(0), format='(i5.5)'),/REMOVE_ALL) + $
      '_'   + STRCOMPRESS(STRING(energy(1), format='(i5.5)'),/REMOVE_ALL) + $
      '_SC' + strcompress(sat,/remove_all) + $
      '_UN' + STRUPCASE(strcompress(units_name,/remove_all)) + $
      '_PR' + STRCOMPRESS(string(prod),/REMOVE_ALL) + $
      '_SP' + strcompress(specie,/remove_all) + $
      '_ET' + STRCOMPRESS(eff_table,/REMOVE_ALL) ;name of TPLOT structure

    IF fc GT 0 THEN BEGIN
      ffc = ffc + 1
      IF ffc EQ 1 THEN BEGIN ; first recovered file
        
        tplot_restore, filename = file_path
        get_data, tplot_var_name, data = data, index = index
        IF index EQ 0 THEN BEGIN
          ffc = ffc - 1
        ENDIF ELSE BEGIN
          
            ier = WHERE(data.e GE energy(0) AND data.e LE energy(1), ierc)
            IF ierc GT 0 THEN BEGIN
              
              strdum = {x: data.x, $
                        y: TOTAL(data.y(ier,*,*),1) / ierc, $
                        v: data.v}
              
              store_data, store_name, data=strdum
              
            ENDIF
            
          
        ENDELSE
        
      ENDIF ELSE BEGIN ; subsequently recovered files
        
        ; get_old_data
        get_data, store_name, data=data_old, dlim=dlim, lim=lim
        
        ;get new data
        tplot_restore, filename = file_path
        get_data, tplot_var_name, data = data, index = index
        IF index EQ 0 THEN BEGIN
          ffc = ffc - 1
        ENDIF ELSE BEGIN
          
            ier = WHERE(data.e GE energy(0) AND data.e LE energy(1), ierc)
            IF ierc GT 0 THEN BEGIN
              
              strdum = {x: data.x, $
                        y: TOTAL(data.y(ier,*,*),1) / ierc, $
                        v: data.v}
              
              data_new = strdum
              
              ;combine the two
              pa_spec_combine_days, $
                store_name, $
                data_old, $
                data_new, $
                dlim, $
                lim
              
            ENDIF
            
          
        ENDELSE
        
      ENDELSE
    ENDIF
  ENDFOR
  
END


;+
;-
PRO pa_spec_combine_days, tplot_var_name, data_old, data_new, dlim, lim

  ;Initialize the data arrays with the data_old values
  xall = data_old.x
  yall = data_old.y
  vall = data_old.v
  
  xall = [xall, data_new.x]
  yall = [yall, data_new.y]
  vall = [vall, data_new.v]
      
  sort_ind = SORT(xall)
  
  xall_new = xall(sort_ind)
  yall_new = yall(sort_ind,*)
  vall_new = vall(sort_ind,*)
  
  datastr = {x:xall_new,y:yall_new,v:vall_new}
  
  store_data, tplot_var_name, data=datastr, dlim=dlim, lim=lim
  
END

PRO bin_size, var
  
  tplot_names, name=name
  get_data, name(var-1), data=data, dlim=dlim, lim=lim
  
  FOR jj = 0, N_ELEMENTS(data.x)-1 DO BEGIN
    FOR ii=0, 15 DO BEGIN
      data.y(jj,ii) = MEAN(data.y(jj,ii*2:ii*2+1),/NaN)
      data.v(jj,ii) = data.v(jj,ii)*2.
    ENDFOR
  ENDFOR
  
  store_data, name(var-1), data={x:data.x, y:data.y(*,0:15), v:data.v(*,0:15)},$
    dlim=dlim, lim=lim

END
