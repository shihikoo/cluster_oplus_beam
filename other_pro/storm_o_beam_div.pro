PRO storm_o_beam_div

  sc = 4
  sc_str = STRING(sc, FORMAT = '(i1.1)')
;  beam_energy = 100 ;eV
;  beam_v = - sqrt(2*beam_energy*1.6e-19/(16*1.66e-27))/1000 ; km/s
;---------------------------------------------------------------------
; Read storm minimum Dst list from file: storm_min_dst.dat
; (this list is extracted from storm list spreadsheet)
; Store these times in variable min_dst (in tplot time format)
;---------------------------------------------------------------------
  OPENR, unit, 'storm_min_dst.dat', /GET_LUN
  
  min_dst = DBLARR(300)         ; assume no more than 300 events
  dummy = ''
  jj = 0
  WHILE NOT EOF(unit) DO BEGIN
      
      READF, unit, dummy
      min_dst(jj) = time_double(dummy)
      jj = jj + 1
      
  ENDWHILE

  min_dst = min_dst(0:jj-1)
  
  CLOSE, unit

;---------------------------------------------------------------------
; Read CLUSTER perigee times list from file: sc4_perigee_times.dat
; (different list for each S/C)
; Store these times in variable petime (in tplot time format)
;---------------------------------------------------------------------
  OPENR, unit, 'sc' + sc_str + '_perigee_times.dat', /GET_LUN

  petime = DBLARR(3000)         ; assumes no more than 3000 perigee passes
  dummy = ''
  jj = 0l
  WHILE NOT EOF(unit) DO BEGIN
    
      READF, unit, dummy
      petime(jj) = time_double(dummy)
      jj = jj + 1
    
  ENDWHILE

  petime = petime(0:jj-1)

  CLOSE, unit

  OPENU, unit, 'log_plotted_d.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[START]'
  FREE_LUN, unit         
 
  OPENU, unit, 'log_errors_d.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[START]'
  FREE_LUN, unit         

;---------------------------------------------------------------------
; For each minimum Dst time choose the time intarval of the plots.
; The time interval is selected is of three orbits duration
;---------------------------------------------------------------------
  FOR ii = 76, N_ELEMENTS(min_dst)-1   DO BEGIN 
   ;-----------------------------------------------------------------
    ; find the time of the perigee pass closest to the minimum Dst
    ;-----------------------------------------------------------------
      before_index = where ((petime - min_dst(ii)) < 0, counts)
      after_index = where ((petime - min_dst(ii)) > 0)          
      before_pe = petime(before_index(counts-1))
      after_pe = petime(after_index(0))
    ;-----------------------------------------------------------------
    ; time interval start time, end time and dt
    ;-----------------------------------------------------------------
      time_start = petime(before_index(counts-2))
      time_end = petime(after_index(1))
      dt = time_end - time_start      
      displaytime =  6. * 3600. 
      idt = CEIL(dt / displaytime)
      
      FOR kk = 0, idt-1  DO BEGIN   

          tplot_names, names = names
          store_data, DELETE = names

          PRINT, STRING(ii) + '   --' + STRING(kk)
      ;----------------------------------------------------------------
      ; Set time interval in tplot
      ;-----------------------------------------------------------------
          t_s = time_start + kk * displaytime
          t_e = t_s + displaytime
          t_dt = t_e - t_s
          timespan, t_s, t_dt, /SECONDS  
;stop
   ;--------------------------Loading data----------------------------------
    
     ;-----------------------------------------------------------------
     ; Load Dst
     ;-----------------------------------------------------------------
;      read_omni
           
     ;------------------------------------------------------------
     ; Load energy spectra - tailward
     ;---------------------------------------------
          sat = [sc, sc]
          specie = [0, 3]
          angle = [[-90, 90], [90, 270]]
          inst = 0 & units_name = 'DIFF FLUX' & eff_table = 0
          plot_en_spec_from_crib, sat, specie, inst, $
            units_name, angle, eff_table, recalc = 0
  
     ;------------------------------------------------------------
     ; Load energy spectra - earthward
     ;----------------------------------------------------------
          sat = [sc, sc]
          specie = [0, 3]
          angle = [[-90, 90], [270, 90]]
          inst = 0 & units_name = 'DIFF FLUX' & eff_table = 0
          plot_en_spec_from_crib, sat, specie, inst, $
            units_name, angle, eff_table, recalc = 0
 
      ;-----------------------------------------------------------------
      ; Load CLUSTER Magnetic field
      ;----------------------------------------------------------------
          sat = sc
          plot_mag_from_crib, $
            sat, $
            POLAR = 1, $
            GSM = 1
      ;-----------------------------------------------------------------
      ; Load CLUSTER moments (n, P)
      ;----------------------------------------------------------------   
          sat = [sc, sc]
          specie = [0, 3]
          moments = ['D', 'D']
          angle = [[-90.0, 90.0], [0., 360.]]
          energy = [40.0, 40000.0]
          inst = 0 & eff_table = 0
          plot_3dmom_from_crib, sat, specie, inst, $
            moments, angle, energy, eff_table, recalc = 0
   ;------------------------------------------------------------------
    ;  Checking weather all data required are loaded
    ;----------------------------------------------------------------  
          s_e = STRARR(4)
          nerror = 0

          tplot_names, 'ENSPEC_SC' + sc_str + $
                       '_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All', names = names
          IF names(0) EQ '' THEN BEGIN 
              s_e(nerror) = ' energy,'
              nerror = nerror+1
          ENDIF 
          
          tplot_names, 'TDMOM_EN00040_40000_SC'+sc_str + $ 
                       '_MTPRESSURE_SP0_ET0_All', names = names
          IF names(0) EQ '' THEN BEGIN 
              s_e(nerror) = ' H_pressure,'
              nerror = nerror+1
          ENDIF 

          tplot_names, 'TDMOM_EN00040_40000_SC'+sc_str + $
                       '_MTPRESSURE_SP3_ET0_All', names = names
          IF names(0) EQ '' THEN BEGIN 
              s_e(nerror) = ' O_pressure,'
              nerror = nerror+1
          ENDIF 
        
          tplot_names,  'MAG_SC'+ sc_str +'_B_xyz_gse_MAG_PR', names = names
          IF names(0) EQ '' THEN BEGIN 
              s_e(nerror)  = ' Mag,'  
              nerror = nerror+1
          ENDIF 

          ;nerror = 0 means all data have been loaded
          IF nerror EQ 0 THEN BEGIN
      ;-------------------------------------------------------------------
      ; Calculate total pressure & beta
      ;--------------------------------------------------------------
              h1_press = 'TDMOM_EN00040_40000_SC' + sc_str + '_MTPRESSURE_SP0_ET0_All'
              o1_press = 'TDMOM_EN00040_40000_SC' + sc_str + '_MTPRESSURE_SP3_ET0_All'
              mag_press = 'MAG_SC' + sc_str +'_B_xyz_gse_MAG_PR'

              plasma_beta, h1_press, mag_press, O1_PRESSURE = o1_press
    
      ;-----------------------------------------------------------------
      ; Load CLUSTER ephemeris
      ;-----------------------------------------------------------------    
              sat = [sc] 
              get_cluster_ephemeris, sat, /GSE_X, /GSE_Y, /GSE_Z, /DIST, /MLT 
     
           
    ;--------------------------------------------------------------------
    ; Using pitch angle to find the O+ beam 
    ;-------------------------------------------------------------------
       ;------------------------------------------------------------- 
       ;Finding energy peak with tailward and earthward diff spectra
       ;Plot pitch angle around energy peak in given range
       ;------------------------------------------------------------

        ;---Read energy spctra data of tailward and earthward
              o_t = 'ENSPEC_SC' + sc_str + $
                    '_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All'
          
              get_data,  o_t,  data = ot, dlim = dlim, lim = lim
              
              ot_start = sort(ABS(ot.x - t_s))
              ot_end = sort(ABS(ot.x - t_e))

              time_ot = ot.x(ot_start(0):ot_end(0))
              flux_ot = ot.y(ot_start(0):ot_end(0), *)
              energy_ot = ot.v(ot_start(0):ot_end(0), *)
              ntime_ot = N_ELEMENTS(time_ot)

              str = {x:time_ot, y:flux_ot, v:energy_ot}
              store_data, o_t, data = str, dlim = dlim, lim = lim
              
                    
              o_e = 'ENSPEC_SC' $ 
                    + STRCOMPRESS(sc,  /REMOVE_ALL)$
                    + '_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All'

              get_data, o_e, data = oe, dlim = dlim, lim = lim
              
              oe_start = sort(ABS(oe.x - t_s))
              oe_end = sort(ABS(oe.x - t_e))

              time_oe = oe.x(oe_start(0):oe_end(0))
              flux_oe = oe.y(oe_start(0):oe_end(0), *)
              energy_oe = oe.v(oe_start(0):oe_end(0), *)
              ntime_oe = N_ELEMENTS(time_oe)

              str = {x:time_oe, y:flux_oe, v:energy_oe}
              store_data, o_e, data = str, dlim = dlim, lim = lim
   
         ;---average them over average_time
              average_time = 5 * 60 ;seconds
              
              average_time = fix((time_ot(ntime_ot-1) - time_ot(0))/2 <  $
                                 (time_oe(ntime_oe-1) - time_oe(0))/2 < $
                                 average_time)
                    
              IF average_time GT 0 THEN BEGIN 
                  IF average_time GE (time_ot(ntime_ot-1)-time_ot(0))/ntime_ot AND $
                     average_time GE (time_oe(ntime_oe-1)-time_oe(0))/ntime_oe THEN BEGIN
                  
                  at_str = STRCOMPRESS(average_time,  /REMOVE_ALL) 
   
                  average_tplot_variable, o_t, at_str, /new
                  average_tplot_variable, o_e, at_str, /new
    
        ;---------------------------------------------------------
         ;Calculate tailward diff from the average data 
 ;         time_diff = time_ota
  ;        energy_diff = energy_ota
   ;       flux_diff = (flux_ota-flux_oea) > 0
         
    ;      diff_name = 'ENSPEC_SC'+ STRCOMPRESS(sc,  /REMOVE_ALL)$
     ;                + '_IN0_PHI90_270_DIFF_UNDIFFFLUX_SP3_ET0_All_AVG'+at_str
        
      ;    str = {x: time_diff, y: flux_diff, v: energy_diff}
       ;   store_data, diff_name, data = str, dlim = dlim, lim = lim

        ; Calculate earthward diff from the average data 
    ;      time_diff = time_oea
     ;     energy_diff = energy_oea
      ;    flux_diff = (flux_oea-flux_ota) > 0
          
       ;   diff_name = 'ENSPEC_SC'+ STRCOMPRESS(sc,  /REMOVE_ALL)$
        ;             + '_IN0_PHI270_90_DIFF_UNDIFFFLUX_SP3_ET0_All_AVG'+at_str
        
         ; str = {x: time_diff, y: flux_diff, v: energy_diff}
        ;  store_data, diff_name, data = str, dlim = dlim, lim = lim

         ;---------------------------------------------------------------

       ;      earth_en_name = diff_name
        ;  tail_en_name = diff_name          

                  en_names = ['', '_h', '_m', '_l']
                  n_en_names = N_ELEMENTS(en_names)

                  tail_en_name = o_t + '_AVG'+at_str              
                  energy_div, tail_en_name             
                                   
                  FOR i_en = 0, n_en_names-1 DO BEGIN 
                      find_energy_peak, tail_en_name+en_names(i_en), tail_epcut
                  
                      sat  = [sc] 
                      specie = [3]        
                      units_name = 'DIFF FLUX' & inst = 0 & eff_table = 0
        
                      plot_pa_spec_around_energy_peak, sat, specie, inst, $
                                                   units_name, eff_table, $
                                                   invar = tail_epcut, $
                                                   outvar = tail_pa, $
                                                   average_time = average_time, $
                                                   recalc = 1, $
                                                   n_range = 1, $
                                                   fixtherange = 1
                      timespan, t_s, t_dt, /SECONDS

                      find_pa_peak, tail_pa, tail_pap
                      
                      beam_filter, tail_pap, tail_epcut, tail_pap_beam, $
                                   tail_epcut_beam, direction = 't'
                  ENDFOR 

             ;     combine_hml_epcut, tail_epcut_beam+en_names(0), $
              ;                       tail_epcut_beam+en_names(1), $
               ;                      tail_epcut_beam+en_names(2), $
                ;                     tail_epcut_beam+en_names(3), $
                 ;                    tail_epcut_beam+'_c'

           ;----------------------------------------
                  earth_en_name = o_e + '_AVG'+at_str
                  energy_div, earth_en_name 

                  FOR i_en = 0, n_en_names-1 DO BEGIN 
                      find_energy_peak, earth_en_name+en_names(i_en), earth_epcut
                  
                      sat  = [sc] 
                      specie = [3]        
                      units_name = 'DIFF FLUX' & inst = 0 & eff_table = 0
        
                      plot_pa_spec_around_energy_peak, sat, specie, inst, $
                                                   units_name, eff_table, $
                                                   invar = earth_epcut, $
                                                   outvar = earth_pa, $
                                                   average_time = average_time, $
                                                   recalc = 1, $
                                                   n_range = 1, $
                                                   fixtherange = 1
                      timespan, t_s, t_dt, /SECONDS

                      find_pa_peak, earth_pa, earth_pap
                      
                      beam_filter, earth_pap, earth_epcut, earth_pap_beam, $
                                   earth_epcut_beam, direction = 'e'
                   
                  ENDFOR 
       
           ;       combine_hml_epcut, earth_epcut_en_names(0), $
            ;                         earth_epcut_beam+en_names(1), $
             ;                        earth_epcut_beam+en_names(2), $
              ;                       earth_epcut_beam+en_names(3), $
               ;                      earth_epcut_beam+'_c'
                 
       ;--------------------------------------------------------------
       ;Overview plots
       ;--------------------------------------------------------------
   
               ;  p01 = 'TDMOM_EN00040_40000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_P_total'
                  p02 = 'TDMOM_EN00040_40000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_beta'
              
               ;  p03 = 'ENSPEC_SC'+ sc_str +'_IN0_PHI90_270_UNDIFFFLUX_SP0_ET0_All'
                  p04 = 'ENSPEC_SC'+ sc_str +'_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All'
               ;  p05 = 'ENSPEC_SC'+ sc_str +'_IN0_PHI90_270_UNDIFFFLUX_SP0_ET0_All'
                  p06 = 'ENSPEC_SC'+ sc_str +'_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All'
            
               ;  p07 = 'Dst_Index'               
               ;  p08 = 'MAG_SC4_B_xyz_gse_GSM_X'

                  p09 = tail_en_name ;avg en
                  p10 = earth_en_name ;avg en

                  pos = STREGEX(p09, '_AVG')

                  p11 = 'PA'+ STRMID(P09, 2, pos+5)  ;avg pa
                  p12 = 'PA'+ STRMID(P10, 2, pos+5)  ;avg pa
 
                  options, '*', 'panel_size', 1
                  options, [p06], 'ztitle', ''
             ;    ylim, p01, 0.01, 3, 1
                  ylim, p02, 0.01, 10
               ;  zlim, p03, 0.1, 100
                  zlim, p04, 0.1, 100
                ; zlim, p05, 0.1, 100 
                  zlim, p06, 0.1, 100

 ;                options, p03, 'ytitle', 'SC' + sc_str + '  H!U+!N (eV)' + '!C!C' + 'Tailward'
                  options, p04, 'ytitle', 'SC' + sc_str + '   O!U+!N (eV)' + '!C!C' + 'Tailward'
   ;              options, p05, 'ytitle', 'SC' + sc_str + '   H!U+!N (eV)' + '!C!C' + 'Earthward'
                  options, p06, 'ytitle', 'SC' + sc_str + '   O!U+!N (eV)' + '!C!C' + 'Earthward'
     ;            options, p08, 'ytitle', 'SC' + sc_str + $
           ;                   '!C!CBx (nT)'

                  FOR i_en = 0, n_en_names-1 DO BEGIN 

                      options, [p09+en_names(i_en), p10+en_names(i_en), $
                                p11+en_names(i_en), p12+en_names(i_en), $
                                p11+en_names(i_en)+'_pap', p12+en_names(i_en)+'_pap', $
                                p11+en_names(i_en)+'_pap_beam', $
                                p12+en_names(i_en)+'_pap_beam'], 'ztitle', ''

                      zlim, p09+en_names(i_en), 0.1, 100
                      zlim, p10+en_names(i_en), 0.1, 100
                      zlim, p11+en_names(i_en), 0.1, 100
                      zlim, p12+en_names(i_en), 0.1, 100
                      zlim, p11+en_names(i_en)+'_pap', 0.1, 100
                      zlim, p12+en_names(i_en)+'_pap', 0.1, 100
                      zlim, p11+en_names(i_en)+'_pap_beam', 0.1, 100
                      zlim, p12+en_names(i_en)+'_pap_beam', 0.1, 100  
    
                      options, p09+en_names(i_en), 'ytitle', 'SC' + sc_str + $
                               '  O!U+!N (eV)' + '!C!C' + 'T'+' (avg '+at_str+'s)'
                      options, p10+en_names(i_en), 'ytitle', 'SC' + sc_str + $
                               '  O!U+!N (eV)' + '!C!C' + 'E'+' (avg '+at_str+'s)'
                      options, p11+en_names(i_en), 'ytitle', $
                               'Pitch Angle)' + '!C!C' + ' (avg '+at_str+'s)'
                      options, p12+en_names(i_en), 'ytitle', $
                               'Pitch Angle' + '!C!C' + ' (avg '+at_str+'s)'

                      options, p11+en_names(i_en)+'_pap', 'ytitle', $
                               'PA  Peak' + '!C!C' + 'E----------T '
                      options, p12+en_names(i_en)+'_pap', 'ytitle', $
                               'PA  Peak' + '!C!C' + 'E----------T '

                      options, p11+en_names(i_en)+'_pap_beam', 'ytitle', $
                               'beam(T'+en_names(i_en)+')' + '!C!C' + 'E----------T '
                      options, p12+en_names(i_en)+'_pap_beam', 'ytitle', $
                               'beam(E'+en_names(i_en)+')' + '!C!C' + 'E----------T '
     
                  ENDFOR 
                                    
                  var_label = 'EPH_SC' + sc_str + '_'
                  var_label = var_label + ['MLT', 'GSE_X', 'GSE_Y', 'GSE_Z', 'DIST']

      ;-----------------------------------------------------------------
      ; Plot in PS file
      ;-----------------------------------------------------------------
                  dd = time_string(t_s)  
                  FOR i_en = 0, n_en_names-1 DO BEGIN 
                      
                      fln = 'plots/tem_d/storm_o_beam_' + $
                            STRMID(dd, 0, 4) + STRMID(dd, 5, 2) + STRMID(dd, 8, 2) + '_' +$
                            STRMID(dd, 11, 2) + STRMID(dd, 14, 2) + STRMID(dd, 17, 2) + $ 
                            'tail_'+en_names(i_en) +'.ps'
      
                      popen, fln
         
                      tplot, [p02, p09+en_names(i_en), p11+en_names(i_en), p11+en_names(i_en)+'_pap', $
                              p04, p11+en_names(i_en)+'_pap_beam'], var_label = var_label
                                       
                      tplot_panel, v = p09+en_names(i_en), $
                                   o = p09+en_names(i_en) + '_epcut', psym = -2
                      tplot_panel, v = p04, o = p09+en_names(i_en)+'_epcut_beam', psym = 0
 
                      pclose

                      ;----earth----
                      fln = 'plots/tem_d/storm_o_beam_' + $
                            STRMID(dd, 0, 4) + STRMID(dd, 5, 2) + STRMID(dd, 8, 2) + '_' +$
                            STRMID(dd, 11, 2) + STRMID(dd, 14, 2) + STRMID(dd, 17, 2) + $ 
                            'earth_'+en_names(i_en) +'.ps'
                   
                      popen, fln
                                         
                      tplot, [p02, p10+en_names(i_en), p12+en_names(i_en), p12+en_names(i_en)+'_pap',$
                              p06, p12+en_names(i_en)+'_pap_beam'], var_label = var_label

                      tplot_panel, v = p10+en_names(i_en), $
                                   o = p10+en_names(i_en) + '_epcut', psym = -2
                      tplot_panel, v = p06, o = p10+en_names(i_en)+'_epcut_beam', psym = 0
                      
                      pclose
                  ENDFOR 

        ;---------------------------------------------------
        ;write energy spectra errors into the log
        ;------------------------------------------------------
                  OPENU, unit, 'log_plotted_d.txt', /GET_LUN, /APPEND
                  PRINTF, unit,  TIME_STRING(min_dst(ii)) + ' [' + $ 
                          STRING(ii, FORMAT = '(i2.2)') + ',' $
                          + STRING(kk, FORMAT = '(i2.2)')+'] '+ $
                          TIME_STRING(t_s)+'   --------------OK'
                  FREE_LUN, unit
              ENDIF  ELSE BEGIN 
                  OPENU, unit, 'log_errors_d.txt', /GET_LUN, /APPEND
                  PRINTF, unit, TIME_STRING(min_dst(ii)) + ' [' + $ 
                          STRING(ii, FORMAT = '(i2.2)') + ',' $
                          + STRING(kk, FORMAT = '(i2.2)')+'] '+ $
                          TIME_STRING(t_s)+' Energy data doubt '
                  FREE_LUN, unit 
              ENDELSE 
          ENDIF ELSE BEGIN 
              OPENU, unit, 'log_errors_d.txt', /GET_LUN, /APPEND
              PRINTF, unit, TIME_STRING(min_dst(ii)) + ' [' + $ 
                      STRING(ii, FORMAT = '(i2.2)') + ',' $
                      + STRING(kk, FORMAT = '(i2.2)')+'] '+ $
                      TIME_STRING(t_s)+' Energy data incorrect '
              FREE_LUN, unit 
          ENDELSE 
          ENDIF ELSE BEGIN 
        ;---------------------------------------------------------
        ;write errors info into log
        ;-----------------------------------------------------------
              OPENU, unit, 'log_errors_d.txt', /GET_LUN, /APPEND
              PRINTF, unit, TIME_STRING(min_dst(ii)) + ' [' $ 
                      + STRING(ii, FORMAT = '(i2.2)') + ',' $
                      + STRING(kk, FORMAT = '(i2.2)')+'] ' $
                      + TIME_STRING(t_s) + ' Missing data: ' $
                      , s_e
              FREE_LUN, unit         
          ENDELSE 
      ENDFOR 
  ENDFOR  

  OPENU, unit, 'log_plotted_d.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[END]'
  FREE_LUN, unit         
 
  OPENU, unit, 'log_errors_d.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[END]'
  FREE_LUN, unit        

END
