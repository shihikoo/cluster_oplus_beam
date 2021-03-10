PRO storm_o_beam_c

  sc = 4
  sc_str = STRING(sc, FORMAT = '(i1.1)')
  ps = 1 ; program will plot ps files only when ps EQ 1
  dumpdata = 0; program will dump data out only when dumpdata EQ 1 
   ;seconds. 
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

  OPENU, unit, 'log_plotted_c.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[START]'
  FREE_LUN, unit         
 
  OPENU, unit, 'log_errors_c.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[START]'
  FREE_LUN, unit         

;---------------------------------------------------------------------
; For each minimum Dst time choose the time intarval of the plots.
; The time interval is selected is of three orbits duration
;---------------------------------------------------------------------
  FOR ii = 0, N_ELEMENTS(min_dst)-1   DO BEGIN 
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
      
      FOR kk = 3, idt-1  DO BEGIN   

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

   ;--------------------------Loading data----------------------------------
    
     ;-----------------------------------------------------------------
     ; Load Dst
     ;-----------------------------------------------------------------
;      read_omni
           
     ;------------------------------------------------------------
     ; Load energy spectra - tailward
     ;---------------------------------------------
          sat = [sc]
          specie = [3]
          angle = [[-90, 90], [90, 270]]
          inst = 0 & units_name = 'DIFF FLUX' & eff_table = 0
          plot_en_spec_from_crib, sat, specie, inst, $
            units_name, angle, eff_table, recalc = 0
  
     ;------------------------------------------------------------
     ; Load energy spectra - earthward
     ;----------------------------------------------------------
          sat = [sc]
          specie = [3]
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
       ;Finding energy peaks with tailward and earthward spectra in 
       ;different three different energy range: 
       ; low: 30-200eV, mid: 200-6000eV, high: 6000-40000eV           
       ;Plot pitch angle around all energy peaks in given range
       ;------------------------------------------------------------

           ;---Read energy spctra data of tailward and earthward 
           ;   keep only data during displaytime 
           ;   and restore them in original files but only-----------

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
              average_time = 5 * 60
              
              average_time = fix((time_ot(ntime_ot-1) - time_ot(0))/2 <  $
                                 (time_oe(ntime_oe-1) - time_oe(0))/2 < $
                                 average_time)
              at_str = STRCOMPRESS(average_time,  /REMOVE_ALL) 
      
              IF average_time GT 0 THEN BEGIN    ; if energy data has more than one
                                                 ; data during the display time
                  IF average_time GE (time_ot(ntime_ot-1)-time_ot(0))/ntime_ot AND $
                     average_time GE (time_oe(ntime_oe-1)-time_oe(0))/ntime_oe THEN BEGIN

               ; and if the avergae time is longer than
               ; the time difference between two time data points 
               ; avergae the energy plot over average time

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
       ;      earth_en_name = diff_name
        ;  tail_en_name = diff_name     
         ;---------------------------------------------------------------   

           ; Divide tailward energy spectra into 3 energy range 
                  en_names = ['', '_h', '_m', '_l']  
                  n_en_names = N_ELEMENTS(en_names)
                  tail_en_name = o_t + '_AVG'+at_str              

                  energy_div, tail_en_name        
     
           ; run a set of routines to get the beam plot for 
           ; all the energy spectra divisions          
                  tail_epcut_beam_names = STRARR(4)                 
                  
                  FOR i_en = 0, n_en_names-1 DO BEGIN 
                      
                      find_energy_peak_from_enspec, tail_en_name+en_names(i_en), tail_epcut, $
                        drop_bad_bin = 0
                  
                      sat  = [sc] 
                      specie = [3]        
                      units_name = 'DIFF FLUX' & inst = 0 & eff_table = 0
        
                      plot_pa_spec_around_energy_peak, sat, specie, inst, $
                                                   units_name, eff_table, $
                                                   invar = tail_epcut, $
                                                   outvar = tail_pa, $
                                                   average_time = average_time,$
                                                   recalc = 1, $
                                                   n_range = 1, $
                                                   fixtherange = 1, $
                                                   PaBin = 22.5
                      tplot_names, 'PASPEC_EN*SC4_UNDIFFFLUX_SP3_All', names = names
                      store_data, delete = names
                      
                      timespan, t_s, t_dt, /SECONDS
;stop
                      find_pa_peak, tail_pa, tail_pap
                      
                      beam_filter, tail_pap, tail_epcut, tail_pap_beam, $
                                   tail_epcut_beam, direction = 't'

                      tail_epcut_beam_names(i_en) = tail_epcut_beam 
                  ENDFOR 

                  ; combine all the beam cuts
                  combine_hml_epcut, tail_epcut_beam_names(0), $
                                     tail_epcut_beam_names(1), $
                                     tail_epcut_beam_names(2), $
                                     tail_epcut_beam_names(3), $
                                     tail_epcut_beam_names(0)+'_c'

           ;----------------------------------------
             ;Repeat everything for earthward
                
                  earth_en_name = o_e + '_AVG'+at_str
                  energy_div, earth_en_name 

                  earth_epcut_beam_names = STRARR(4)

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
                                                   fixtherange = 1, $
                                                   PaBin = 22.5
                      timespan, t_s, t_dt, /SECONDS

                      find_pa_peak, earth_pa, earth_pap
                      
                      beam_filter, earth_pap, earth_epcut, earth_pap_beam, $
                                   earth_epcut_beam, direction = 'e'
                      
                      earth_epcut_beam_names(i_en) = earth_epcut_beam
                   
                  ENDFOR 
       
                  combine_hml_epcut, earth_epcut_beam_names(0), $
                                     earth_epcut_beam_names(1), $
                                     earth_epcut_beam_names(2), $
                                     earth_epcut_beam_names(3), $
                                     earth_epcut_beam_names(0)+'_c'
                 
       ;--------------------------------------------------------------
       ;Overview plots
       ;--------------------------------------------------------------
          ;----------------------
          ;reset plots options
          ;----------------------

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
                  
                  p13 = tail_epcut_beam_names(0)+'_c'
                  p14 = earth_epcut_beam_names(0)+'_c'

                  options, '*', 'panel_size', 1
                  options, [p06], 'ztitle', ''
             ;    ylim, p01, 0.01, 3, 1
                  ylim, p02, 0.01, 10
               ;  zlim, p03, 0.1, 100
                  zlim, p04, 0.1, 100
                ; zlim, p05, 0.1, 100 
                  zlim, p06, 0.1, 100

 ;                options, p03, 'ytitle', 'SC' + sc_str + '  H!U+!N (eV)' + '!C!C' + 'Tailward'
                  options, p04, 'ytitle', 'SC' + sc_str + '  O!U+!N (eV)' + '!C!C' + 'Tailward'
   ;              options, p05, 'ytitle', 'SC' + sc_str + '  H!U+!N (eV)' + '!C!C' + 'Earthward'
                  options, p06, 'ytitle', 'SC' + sc_str + '  O!U+!N (eV)' + '!C!C' + 'Earthward'
          ;       options, p08, 'ytitle', 'SC' + sc_str + '!C!CBx (nT)'
                                  
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
                               '  O!U+!N (eV)' + '!C!C' + 'T'+en_names(i_en)+' (avg'+at_str+'s)'
                      options, p10+en_names(i_en), 'ytitle', 'SC' + sc_str + $
                               '  O!U+!N (eV)' + '!C!C' + 'E'+en_names(i_en)+' (avg'+at_str+'s)'
                      options, p11+en_names(i_en), 'ytitle', 'SC' + sc_str + $
                               '  O!U+!N (PA)' + '!C!C' + 'T'+en_names(i_en)+' (avg'+at_str+'s)'
                      options, p12+en_names(i_en), 'ytitle', 'SC' + sc_str + $
                               '  O!U+!N (PA)' + '!C!C' + 'E'+en_names(i_en)+' (avg'+at_str+'s)'

                      options, p11+en_names(i_en)+'_pap', 'ytitle', $
                               'O!U+!N PA  Peak' + '!C!C' + 'E-------T'
                      options, p12+en_names(i_en)+'_pap', 'ytitle', $
                               'O!U+!N PA  Peak' + '!C!C' + 'E-------T'

                      options, p11+en_names(i_en)+'_pap_beam', 'ytitle', $
                               'T'+en_names(i_en)+ '!C!C!C' + 'E-------T'
                      options, p12+en_names(i_en)+'_pap_beam', 'ytitle', $
                               'E'+en_names(i_en)+ '!C!C!C' + 'E-------T'
stop
                  ENDFOR 
                                    
                  var_label = 'EPH_SC' + sc_str + '_'
                  var_label = var_label + ['MLT', 'GSE_X', 'GSE_Y', 'GSE_Z', 'DIST']
                  
              ;-------------------------------
              ; Plots the graph in idl window
              ;-------------------------------    
                  tplot, [p02, p04, p11+'_h_pap_beam', p11+'_m_pap_beam', p11+'_l_pap_beam', $
                          p06, p12+'_h_pap_beam', p12+'_m_pap_beam', p12+'_l_pap_beam'], $
                         var_label = var_label
                  tplot_panel, v = p04, o = p13, psym = 0
                  tplot_panel, v = p06, o = p14, psym = 0

              ;---------------------------------------
              ; Plot the graph in PS file if ps is set to 1 
              ;---------------------------------------
                  ts = time_string(t_s)  
                  te = time_string(t_e)
                  date_s = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
                  time_s = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
                  date_e = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2)
                  time_e = STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2)
                  
                  IF ps EQ 1 THEN BEGIN  
                     fln = 'plots/tem_c/storm_o_beam_' + date_s + '_' + time_s + '.ps'
      
                     popen, fln
                     
                     tplot, [p02, p04, p11+'_h_pap_beam', p11+'_m_pap_beam', $
                             p11+'_l_pap_beam', $
                             p06, p12+'_h_pap_beam', p12+'_m_pap_beam', p12+'_l_pap_beam'], $
                            var_label = var_label
                     tplot_panel, v = p04, o = p13, psym = 0
                     tplot_panel, v = p06, o = p14, psym = 0
                     
                     pclose

                 ENDIF         
               
               ;--------------------------------------
               ;dump the data out if dumpdata is set to 1  
               ;--------------------------------------
  
                 IF dumpdata EQ 1 THEN BEGIN 

                     get_data, tail_epcut_beam_names(0), data = data
                     time_dd = data.x
                     
                     n_time = N_ELEMENTS(time_dd) 
                     
                     title_dd = STRARR(n_time, 19 )
                     data_dd = DBLARR(n_time, 19 )
                     
                     FOR i_time = 0, n_time-1 DO BEGIN 
                         title_dd(i_time, *) = ['        flag      ', $
                                                '    en_t_high  ', $
                                                '    en_t_mid   ', $
                                                '    en_t_low   ', $
                                                '    en_e_hige  ',$
                                                '    en_e_mid   ',$
                                                '    en_e_low   ', $
                                                '    pa_t_hige  ',$
                                                '    pa_t_mid   ',$
                                                '    pa_t_low   ', $
                                                '    pa_e_high  ',$
                                                '    pa_e_mid   ',$
                                                '    pa_e_low   ', $
                                                '   flux_t_high ',$
                                                '   flux_t_mid  ',$
                                                '   flux_t_low  ', $
                                                '   flux_e_high ',$
                                                '   flux_e_mid  ',$
                                                '   flux_e_low  ']  
                     ENDFOR 
                     
                     FOR i_en = 1, n_en_names-1 DO BEGIN 
                       ;tail energy  
                         get_data, tail_epcut_beam_names(i_en), data = data
                        
                         index = where (~FINITE(data.y), ct)
                         IF ct GT 0 THEN BEGIN 
                             data.y(index) = 0
                         ENDIF 
                         
                         data_dd(*, 0) = data_dd(*, 0) + (data.y GT 0)
                         data_dd(*, i_en) = data.y
    
                       ; earth energy
                         get_data, earth_epcut_beam_names(i_en), data = data
                         
                         index = where ( ~FINITE(data.y), ct)
                         IF ct GT 0 THEN BEGIN 
                             data.y(index) = 0
                         ENDIF 

                         data_dd(*, 0) = data_dd(*, 0) + (data.y GT 0)
                         data_dd(*, 3+i_en) = data.y
  
                       ; tail pap
                         get_data, p11+en_names(i_en)+'_pap_beam', data = data

                         index = where( ~FINITE(data.y), ct)
                         IF ct GT 0 THEN BEGIN 
                             data.y(index) = 0
                             data.v(index) = 0
                         ENDIF 

                         data_dd(*, 6+i_en ) = total(data.v(*, *), 2)
                         data_dd(*, 12+i_en) = total(data.y(*, *), 2)
            
                        ; earth pap
                         get_data, p12+en_names(i_en)+'_pap_beam', data = data
                     
                         index = where( ~FINITE(data.y), ct)
                         IF ct GT 0 THEN BEGIN 
                             data.y(index) = 0
                             data.v(index) = 0
                         ENDIF 
                       
                         data_dd(*, 9+i_en) = total(data.v(*, *), 2)
                         data_dd(*, 15+i_en) = total(data.y(*, *), 2)  
                
                     ENDFOR 
                 
                     IF date_s EQ date_e THEN BEGIN 
                         str = {x:time_dd, y:data_dd, v:title_dd}
                         store_data, 'dump_data', data = str
                         dump_data, 'dump_data', $
                                    file_out = 'data/storm_o_beam_'+date_s+'.dat'

                     ENDIF ELSE BEGIN 
                         
                         midnight = time_double(STRMID(te, 0, 10)+'/00:00:00')
                         fday = where(time_dd LT midnight)
                         sday = where(time_dd GE midnight)
                         
                         str = {x:time_dd(fday), y:data_dd(fday, *), v:title_dd(fday, *)}
                         store_data, 'dump_data_f', data = str
                         dump_data, 'dump_data_f', $
                                    file_out = 'data/storm_o_beam_'+date_s+'.dat'
                         
                         str = {x:time_dd(sday), y:data_dd(sday, *), v:title_dd(sday, *)}
                         store_data, 'dump_data_s', data = str
                         dump_data, 'dump_data_s', $
                                    file_out = 'data/storm_o_beam_'+date_e+'.dat'
                 ENDELSE 
;stop
             ENDIF 
;stop
        ;---------------------------------------------------
        ;write plotted info into the plotted log or
        ;energy spectra errors into the error log 
        ;------------------------------------------------------
                  OPENU, unit, 'log_plotted_c.txt', /GET_LUN, /APPEND
                  PRINTF, unit,  TIME_STRING(min_dst(ii)) + ' [' + $ 
                          STRING(ii, FORMAT = '(i2.2)') + ',' $
                          + STRING(kk, FORMAT = '(i2.2)')+'] '+ $
                          TIME_STRING(t_s)+'   --------------OK'
                  FREE_LUN, unit
              ENDIF  ELSE BEGIN 
                  OPENU, unit, 'log_errors_c.txt', /GET_LUN, /APPEND
                  PRINTF, unit, TIME_STRING(min_dst(ii)) + ' [' + $ 
                          STRING(ii, FORMAT = '(i2.2)') + ',' $
                          + STRING(kk, FORMAT = '(i2.2)')+'] '+ $
                          TIME_STRING(t_s)+' Energy data doubt '
                  FREE_LUN, unit 
              ENDELSE 
          ENDIF ELSE BEGIN 
              OPENU, unit, 'log_errors_c.txt', /GET_LUN, /APPEND
              PRINTF, unit, TIME_STRING(min_dst(ii)) + ' [' + $ 
                      STRING(ii, FORMAT = '(i2.2)') + ',' $
                      + STRING(kk, FORMAT = '(i2.2)')+'] '+ $
                      TIME_STRING(t_s)+' Energy data incorrect '
              FREE_LUN, unit 
          ENDELSE 
          ENDIF ELSE BEGIN 
        ;---------------------------------------------------------
        ;write missing data errors info into log
        ;-----------------------------------------------------------
              OPENU, unit, 'log_errors_c.txt', /GET_LUN, /APPEND
              PRINTF, unit, TIME_STRING(min_dst(ii)) + ' [' $ 
                      + STRING(ii, FORMAT = '(i2.2)') + ',' $
                      + STRING(kk, FORMAT = '(i2.2)')+'] ' $
                      + TIME_STRING(t_s) + ' Missing data: ' $
                      , s_e
              FREE_LUN, unit         
          ENDELSE 
      ENDFOR 
  ENDFOR  

  OPENU, unit, 'log_plotted_c.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[END]'
  FREE_LUN, unit         
 
  OPENU, unit, 'log_errors_c.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[END]'
  FREE_LUN, unit        

END
