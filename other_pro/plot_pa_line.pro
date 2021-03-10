PRO plot_pa_line
  sc = 4
  sc_str = STRING(sc, FORMAT = '(i1.1)')

  start_time = '2002-09-30/07:56:30'
  t_s = time_double(start_time) 
  t_dt =  6 * 3600               ;seconds
  
  timespan, t_s, t_dt, /SECONDS  
;--------------------------------------------------------------------
;Load pitch angle spectra
;-------------------------------------------------------------------
  sat   = [sc]
  specie = [3]        
  units_name = 'DIFF FLUX' & inst = 0 & eff_table = 0 
 
  energy = [40., 40000.]
 

;energy = [500., 3000.]
;energy = [1000., 10000.]
;energy = [100., 1000.]
;energy = [40., 100.]
 
  plot_pa_spec_from_crib, sat, specie, inst, units_name, $
                          energy, eff_table, $
                          PaBin = 22.5, $
                          energy = [40., 100.]  $
                          recalc = 1, $
                          BKG = 0, $
                          COMBINE = 1
;--------------------------------------------------------------------
;decide type of the pitch angle distributions
;-------------------------------------------------------------------
  e_min = STRCOMPRESS(STRING(energy(0), FORMAT = '(i5.5)'), /REMOVE_ALL)
  e_max = STRCOMPRESS(STRING(energy(1), FORMAT = '(i5.5)'), /REMOVE_ALL)
  pa_name = 'PASPEC_EN'+ e_min +'_'+ e_max +'_SC4_UNDIFFFLUX_SP3_All'
  
  get_data, pa_name, data = pa
  
  time_pa = pa.x
  flux_pa = pa.y
  data_pa = pa.v
 
  zlim, pa_name, 0.1, 100.
  tplot, pa_name
  
  cut_time = 5 * 60            ;second 
  icut = CEIL(t_dt / cut_time)-1
          
  FOR jj = 0, icut-1 DO BEGIN 
      
      cut_s = t_s + jj * cut_time
      cut_e = cut_s + cut_time * 2
      cut_dt = cut_e - cut_s
        
      index = where(time_pa GE cut_s AND time_pa LE cut_e, tc)
      
      IF index(0) NE -1 THEN BEGIN 

        ;  timebar, (cut_s + cut_e)/2, color = 1

          y_data = TOTAL(flux_pa(index(0):index(tc-1), *), 1)/tc
          x_data = REFORM(data_pa(0, *))
                                                             
          pa_name_new = pa_name +'_' + time_string(cut_s)
                  
          store_data, pa_name_new, data = {x:x_data, y:y_data}
                
          dd = time_string(t_s)
          cc = time_string(cut_s)

          fln = 'PA_plots_'+ STRMID(dd, 0, 4) + STRMID(dd, 5, 2) $
                + STRMID(dd, 8, 2) + '_' + STRMID(dd, 11, 2) $
                + STRMID(dd, 14, 2) + STRMID(dd, 17, 2) + '/PA_' $
                + STRMID(cc, 0, 4) + STRMID(cc, 5, 2) + STRMID(cc, 8, 2) $
                + '_' + STRMID(cc, 11, 2) + STRMID(cc, 14, 2) $
                + STRMID(cc, 17, 2) +'.ps'
  
          popen, fln

          plot, x_data, y_data, xlog = 0, ylog = 1, $ 
                xrange = [0, 180 ], yrange = [1e0, 1e4]
             
          pclose

      ENDIF ELSE BEGIN 
          print, 'no pitch angle data between' + time_string(cut_s) $
                 + 'and' + time_string(cut_e)                  
      ENDELSE
  ENDFOR
END
