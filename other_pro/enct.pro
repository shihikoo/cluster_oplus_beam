PRO enct
  sc = 4 
  sc_str = STRING(sc, FORMAT = '(i1.1)')

  ps = 0 ; program will plot ps files only when ps EQ 1
  average_time_input = 5 * 60 ;5 min average time
  low_counts_line = 27
  lcl_str = string(low_counts_line,  format = '(i2.2)')
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
;--------------------------------
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
;---------------------------------------------------------------------
  average_time = average_time_input
  at_str = STRCOMPRESS(average_time,  /REMOVE_ALL)

  FOR ii = 26, 43 DO BEGIN ; N_ELEMENTS(min_dst)-1   DO BEGIN 
      before_index = where ((petime - min_dst(ii)) < 0, counts)
      after_index = where ((petime - min_dst(ii)) > 0)          
      before_pe = petime(before_index(counts-1))
      after_pe = petime(after_index(0))

      time_start = petime(before_index(counts-2))
      time_end = petime(after_index(1))
      dt = time_end - time_start      
      displaytime =  6. * 3600. 
      idt = CEIL(dt / displaytime)

      FOR kk = 23, idt-1  DO BEGIN   

          tplot_names, names = names
          store_data, DELETE = names

          PRINT, STRING(ii) + '   --' + STRING(kk)
     
          t_s = time_start + kk * displaytime
          t_e = t_s + displaytime
          t_dt = t_e - t_s       

          timespan, t_s, t_dt, /SECONDS  

          ts = time_string(t_s)  
          te = time_string(t_e)
          date_s = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
          time_s = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
          date_e = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2)
          time_e = STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2)
;----------------------------------------------------------------------------
          sat = [sc]
          specie = [3]
          angle = [[-90, 90], [90, 270]]
          inst = 0 & units_name = 'DIFF FLUX' & eff_table = 0
          plot_en_spec_from_crib, sat, specie, inst, $
            units_name, angle, eff_table, recalc = 1
          ; limit the time of preposseced file
          flux_name = 'ENSPEC_SC' + sc_str + $
                      '_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All'
          
          tplot_names, flux_name, names = names
          IF names(0) NE '' THEN BEGIN 
              get_data, flux_name,  data = data, dlim = dlimf, lim = limf
              index = where(data.x GE (t_s+4) AND data.x LE t_e,  ct)
              IF ct GT 0 THEN BEGIN 
                  time_flux = data.x(index)
                  flux_flux = data.y(index, *)
                  energy_flux = data.v(index, *)
                  store_data, flux_name, $
                              data = {x:time_flux, y:flux_flux, v:energy_flux},  $
                              dlim = dlimf, lim = limf
 ;-------
                  sat = [sc]
                  specie = [3]
                  angle = [[-90, 90], [90, 270]]
                  inst = 0 & units_name = 'Counts'  & eff_table = 0
                  plot_en_spec_from_crib, sat, specie, inst, $
                    units_name, angle, eff_table, recalc = 1
                  
                  counts_name = 'ENSPEC_SC' + sc_str + $
                                '_IN0_PHI90_270_UNCOUNTS_SP3_ET0_All'

                  tplot_names, counts_name, names = names
                  IF names(0) NE '' THEN BEGIN 
                      get_data, counts_name,  data = data, dlim = dlimc, lim = limc
                      index = where(data.x GE (t_s+4) AND data.x LE t_e,  ct)
                      IF ct GT 0 THEN BEGIN 
                          time_counts = data.x(index)
                          counts_counts = data.y(index, *)
                          energy_counts = data.v(index, *)
;stop
 ;-----------------------------------------------------------------------
                          average_tplot_variable, flux_name, at_str, /new
                          
                          tplot_names, flux_name+'_AVG'+at_str, names = names
                          IF names(0) NE '' THEN BEGIN 
                              
                              get_data, flux_name+'_AVG'+at_str, data = data
                              time_flux_a = data.x
                              flux_flux_a = data.y
                              energy_flux_a = data.v
                              ntime = N_ELEMENTS(time_flux_a)
                              nenergy = N_ELEMENTS(energy_counts(0, *))
;stop
                              counts_counts_new = DBLARR(ntime, nenergy)
                              FOR i = 0, ntime-1  DO BEGIN 
                                  index = where(time_counts-time_counts(0) GE i*average_time AND $
                                                time_counts-time_counts(0) LT (i+1)*average_time, ct)
                                  IF ct GT 0  THEN BEGIN 
                                      counts_counts_new(i, *) = total(counts_counts(index, * ), 1)
                                  ENDIF 
                              ENDFOR 
                              
                              store_data, counts_name+'_AVG'+at_str, $
                                          data = {x:time_flux_a, y:counts_counts_new, v:energy_counts},  $
                                          dlim = dlimc, lim = limc
                              
                              counts_filter = counts_counts_new GE low_counts_line
                              
                              store_data, 'ocunts_filter',  $
                                          data = {x:time_flux_a, y:counts_filter, v:energy_flux_a},  $
                                          dlim = dlimc, lim = limc
                              
                              index = where(counts_filter EQ 0, ct)
                              IF ct GT 0 THEN  flux_flux_a(index) = 0
                              flux_flux_a(*, n_elements(flux_flux_a(0, *))-1) = 0

                              store_data, 'final_en',  $
                                          data = {x:time_flux_a, y:flux_flux_a, v:energy_flux_a},  $
                                          dlim = dlimf, lim = limf
                              
                              
                              zlim, 2, 0, 1, 0
                              zlim, 4, 0, low_counts_line, 0
                              zlim, 5, 0, 1, 0
                             ; zlim, 6, 1, 100, 1
                              tplot, [1, 3, 2, 4, 5, 6]

                              stop

                              IF ps EQ 1 THEN BEGIN 
                                  popen, 'plots/enct/'+lcl_str +'/enct_' + date_s + '_' + time_s + '.ps'
                                  
                                  tplot, [1, 3, 2, 4, 5, 6]
                                  
                                  pclose
                              ENDIF 
                          ENDIF 
                      ENDIF 
                  ENDIF 
              ENDIF 
;stop
          ENDIF 
          
      ENDFOR  
  ENDFOR  
  spawn = 'mogrify -format png *.ps'
;stop
END
