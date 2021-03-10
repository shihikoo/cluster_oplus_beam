PRO storm_o_beam

  sc = 4
  sc_str = STRING(sc, FORMAT = '(i1.1)')
  ps = 1 ; program will plot ps files only when ps EQ 1

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

  OPENU, unit, 'log_plotted.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[START]'
  FREE_LUN, unit         
 
  OPENU, unit, 'log_errors.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[START]'
  FREE_LUN, unit       
;---------------------------------------------------------------------
; For each minimum Dst time choose the time intarval of the plots.
; The time interval is selected is of three orbits duration
;---------------------------------------------------------------------
  FOR ii = 30, N_ELEMENTS(min_dst)-1   DO BEGIN 
   ;-----------------------------------------------------------------
    ; find the time of the perigee pass closest to the minimum Dst
    ;------------------------------- ----------------------------------
      before_index = where ((petime - min_dst(ii)) < 0, counts)
      after_index = where ((petime - min_dst(ii)) > 0)          
      before_pe = petime(before_index(counts-1))
      after_pe = petime(after_index(0))
    ;-----------------------------------------------------------------
    ; time interval start time, end time and dt
    ;-----------------------------------------------------------------
  ;    time_start = petime(before_index(counts-2))
  ;    time_end = petime(after_index(1))
      
      time_start = petime(before_index(counts-2))
      time_end = petime(after_index(1))

      dt = time_end - time_start      
      displaytime =  6. * 3600. 
      idt = CEIL(dt / displaytime)
      
      FOR kk = 10, 19  DO BEGIN

          tplot_names, names = names
          store_data, DELETE = names

          PRINT, STRING(ii) + '   --' + STRING(kk)
      ;----------------------------------------------------------------
      ; Set time interval in tplot
      ;-----------------------------------------------------------------
          t_s = time_start + kk * displaytime
          t_e = t_s + 6.*3600.
      ;    t_s = time_double('2002-09-30/12:49:30')
      ;    t_s = time_double('2002-10-01/10:00:00')
       ;   t_e = t_s + 1*3600  ;  /second
     
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
          moments = ['D', 'P']
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
              average_time = 5 * 60 ;seconds
              
              average_time = fix((time_ot(ntime_ot-1) - time_ot(0))/2 <  $
                                 (time_oe(ntime_oe-1) - time_oe(0))/2 < $
                                 average_time)
              at_str = STRCOMPRESS(average_time,  /REMOVE_ALL) 
      
              IF average_time GT 0 THEN BEGIN    ; if energy data has more than one
                                                 ; data during the display time
                  IF average_time GE (time_ot(ntime_ot-1)-time_ot(0))/ntime_ot AND $
                     average_time GE (time_oe(ntime_oe-1)-time_oe(0))/ntime_oe THEN BEGIN
          
                      average_tplot_variable, o_t, at_str, /new
                      average_tplot_variable, o_e, at_str, /new
                      
                      tail_en_name = o_t + '_AVG'+at_str  
                      earth_en_name = o_e+'_AVG'+at_str
      
                      get_data, tail_en_name, data = data
                      time_avg = data.x
                      energy_avg = data.v

                      nenergy = N_ELEMENTS(energy_avg(0, *))-1 ; the last one is broken        
   ;-----------------------------------------------------------------   
                      PaBin = 22.5
                      PaBin_str = STRING(PaBin*10, format = '(i3.3)')
                      npa = FIX(180/PaBin)

                      highen = 0 
                      lowen = nenergy-1
                      nen = lowen-highen+1

                      en_pe = REFORM(energy_avg(0,highen:lowen))
                      pa_pe = FINDGEN(npa)*PaBin+PaBin/2.
                      time_pe = time_avg ; pe = pitch angle vs energy

                      flux_pe = DBLARR(N_ELEMENTS(time_pe), nen, npa)

                      FOR i_en = highen, lowen  DO BEGIN                        
                          i_en_str = STRING(i_en, format = '(i2.2)')
                          
                          str = {x: time_avg, y:energy_avg(*, i_en)}
                          store_data,'ENERGY_AVG',data = str

                          sat  = [sc] 
                          specie = [3]        
                          units_name = 'DIFF FLUX' & inst = 0 & eff_table = 0
        
                          plot_pa_spec_around_energy_peak, sat, specie, inst, $
                                                   units_name, eff_table, $
                                                   invar = 'ENERGY_AVG', $
                                                   outvar = pa_name, $
                                                   average_time = average_time, $
                                                   recalc = 1, $
                                                   n_range = 0, $
                                                   fixtherange = 1, $
                                                   PaBin = PaBin
     
                          timespan, t_s, t_dt, /SECONDS   
                   ;       store_data, DELETE = 'ENERGY_AVG'         
                      
                          get_data, pa_name, data = data
                          flux_pe(*, i_en, *) = data.y
                      ENDFOR    

                      str = {x:en_pe, y:flux_pe, v:pa_pe}
                      store_data, 'pa_en_all', data = str
                      
                      FOR i_t = 0, N_ELEMENTS(time_pe)-1  DO BEGIN 
                          
                          flux_pe_local = REFORM(flux_pe(i_t, *, *))
             ;           flux_pe_local(where(flux_pe_local LE 0)) = !VALUES.F_NAN

                          data = {x: en_pe, y:flux_pe_local, v:pa_pe}
             ;    stop
                ;          window, 1
                 ;         specplot, data = data, NO_INTERP = 1, $
                  ;                  lim = $
                   ;                 {xlog:1, xrange:[40, 40000],$
                    ;                 ystyle:1, yrange:[0, 180], $
                     ;                zlog:1, zrange:[0.1, 100], $
                      ;               title:time_string(time_pe(i_t))+' - ' +$
                       ;              time_string(time_pe(i_t)+average_time), $
                        ;             xtitle:'Energy ( eV )', ytitle:'Pitch  Angle'}
                    
                          tpe = time_string(time_pe(i_t))  
                          dt_pe = STRMID(tpe, 0, 4) + STRMID(tpe, 5, 2) + STRMID(tpe, 8, 2)
                          tm_pe = STRMID(tpe, 11, 2) + STRMID(tpe, 14, 2) + STRMID(tpe, 17, 2)
                       
                          fln = 'plots/pa_en/'+strcompress(ii, /remove_all)$
                                +'/pa_en_'+dt_pe+'_'+tm_pe+'.ps'
                          popen, fln
                         
                          specplot, data = data, NO_INTERP = 1, $
                                    lim = {xlog:1, xrange:[40, 40000], $
                                     ystyle:1, yrange:[0, 180], $
                                     zlog:1, zrange:[0.1, 100], $
                                     title:time_string(time_pe(i_t))+' - ' +$
                                     time_string(time_pe(i_t)+average_time), $
                                     xtitle:'Energy ( eV )', $
                                           ytitle:'T         Pitch  Angle         E', $
                                          ztitle:'0.1                 100'}
;stop
                          pclose
                      ENDFOR 
      
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
                  p08 = 'MAG_SC4_B_xyz_gse_GSM'

                  p09 = tail_en_name ;avg en
                  p10 = earth_en_name ;avg en
                  
                  options, '*', 'panel_size', 1
                  options, [p06, p09, p10], 'ztitle', ''
             ;    ylim, p01, 0.01, 3, 1
                  ylim, p02, 0.01, 10
               ;  zlim, p03, 0.1, 100
                  zlim, p04, 0.1, 100
                ; zlim, p05, 0.1, 100 
                  zlim, p06, 0.1, 100

 ;                options, p03, 'ytitle', 'SC' + sc_str + '  H!U+!N (eV)' + '!C!C' + 'Tailward'
                  options, p04, 'ytitle', '!CT'
   ;              options, p05, 'ytitle', 'SC' + sc_str + '  H!U+!N (eV)' + '!C!C' + 'Earthward'
                  options, p06, 'ytitle', '!CE'
          ;       options, p08, 'ytitle', 'SC' + sc_str + '!C!CBx (nT)'
                  
                  options, 'DR*', 'ztitle', ''
                  zlim, 'PA*', 0.1, 100
                  zlim, 'DR*', 0.1, 100              
                 
                  var_label = 'EPH_SC' + sc_str + '_'
                  var_label = var_label + ['MLT', 'GSE_X', 'GSE_Y', 'GSE_Z', 'DIST']
              ;-------------------------------
              ; Plots the graph in idl window
              ;-------------------------------    
                 tplot, [p08, p04, p06, 'DR*'], var_label = var_label
; stop
              ;---------------------------------------
              ; Plot the graph in PS file if ps is set to 1 
              ;---------------------------------------
                  ts = time_string(t_s)  
                  date_s = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
                  time_s = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
                      
                  IF ps EQ 1 THEN BEGIN  
                      fln = 'plots/pa_en/'+strcompress(ii, /remove_all)$
                            +'/paspec_' + date_s + '_' + time_s $
                            + '_' + PaBin_str $
                            + '_d' + STRCOMPRESS(fix(t_dt/3600),  /REMOVE_ALL)+'hour'$
                            + '_a' + STRCOMPRESS(fix(average_time/60),  /REMOVE_ALL) +'min_'$
                            + STRCOMPRESS(nen, /REMOVE_ALL)+'en.ps'
;    +'_'+STRCOMPRESS(nen,  /REMOVE_ALL) +'en.ps'               
                      popen, fln
                      tplot, [p08, p04, p06, 'DR*'], var_label = var_label
                     
                      pclose
                  ENDIF      
;stop
        ;---------------------------------------------------
        ;write plotted info into the plotted log or
        ;energy spectra errors into the error log 
        ;------------------------------------------------------
                      OPENU, unit, 'log_plotted.txt', /GET_LUN, /APPEND
                      PRINTF, unit,  TIME_STRING(min_dst(ii)) + ' [' + $ 
                              STRING(ii, FORMAT = '(i2.2)') + ',' $
                              + STRING(kk, FORMAT = '(i2.2)')+'] '+ $
                              TIME_STRING(t_s)+'   --------------OK'
                      FREE_LUN, unit
                  ENDIF  ELSE BEGIN  
                      OPENU, unit, 'log_errors.txt', /GET_LUN, /APPEND
                      PRINTF, unit, TIME_STRING(min_dst(ii)) + ' [' + $ 
                              STRING(ii, FORMAT = '(i2.2)') + ',' $
                              + STRING(kk, FORMAT = '(i2.2)')+'] '+ $
                              TIME_STRING(t_s)+' Energy data doubt '
                      FREE_LUN, unit 
                  ENDELSE 
              ENDIF  ELSE BEGIN 
                  OPENU, unit, 'log_errors.txt', /GET_LUN, /APPEND
                  PRINTF, unit, TIME_STRING(min_dst(ii)) + ' [' + $ 
                          STRING(ii, FORMAT = '(i2.2)') + ',' $
                          + STRING(kk, FORMAT = '(i2.2)')+'] '+ $
                          TIME_STRING(t_s)+' Energy data incorrect '
                  FREE_LUN, unit 
              ENDELSE 
          ENDIF  ELSE BEGIN 
        ;---------------------------------------------------------
        ;write missing data errors info into log
        ;-----------------------------------------------------------
              OPENU, unit, 'log_errors.txt', /GET_LUN, /APPEND
              PRINTF, unit, TIME_STRING(min_dst(ii)) + ' [' $ 
                      + STRING(ii, FORMAT = '(i2.2)') + ',' $
                      + STRING(kk, FORMAT = '(i2.2)')+'] ' $
                      + TIME_STRING(t_s) + ' Missing data: ' $
                      , s_e
              FREE_LUN, unit         
          ENDELSE 
      ENDFOR  
  ENDFOR  

  OPENU, unit, 'log_plotted.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[END]'
  FREE_LUN, unit         
 
  OPENU, unit, 'log_errors.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[END]'
  FREE_LUN, unit        

END
