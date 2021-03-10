PRO plot_pa_spec
  sc = 4
  sc_str = STRING(sc, FORMAT = '(i1.1)')
;-----------------------------------------------------------------
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
;-----------------------------------------------------------------
  
  FOR ii = 0, N_ELEMENTS(min_dst)-1 DO BEGIN
      
      before_index = where ((petime - min_dst(ii)) < 0, counts)
      after_index = where ((petime - min_dst(ii)) > 0)    
     
      before_pe = petime(before_index(counts-1))
      after_pe = petime(after_index(0))

    ;-----------------------------------------------------------------
      time_start = petime(before_index(counts-2))
      time_end = petime(after_index(1))
      dt = time_end - time_start

      displaytime = 6. * 60.* 60.;second
      
      idt = CEIL(dt / displaytime)
      
      FOR kk = 0, idt-1 DO BEGIN
  
          t_s = time_start + kk * displaytime
          t_e = t_s + displaytime
          t_dt = t_e - t_s
          timespan, t_s, t_dt, /SECONDS  

;--------------------------------------------------------------------
;Load pitch angle spectra
;-------------------------------------------------------------------
      sat   = [sc]
      specie = [3]        
      units_name = 'DIFF FLUX' 
      inst = 0 
      eff_table = 0 

      all_energy = [[30000., 40000.], [15000., 20000.], [10000., 15000.], $
                    [5000., 10000.],  [4000., 5000.], [2000., 4000.], $
                    [1500., 2000.], [1000., 1500.], [500., 1000.], $
                    [400., 500.], [200., 400.], [100., 200.], $
                    [60., 100.], [40., 60.], [30., 40.], [20., 30.]]


      all_energy = [[31000., 32000.], [19000., 20000.], [11000., 12000.], $
                    [7300., 7400.],  [4500., 4600.], [2800., 2900.], $
                    [1700., 1800.], [1000., 1200.], [650., 660.], $
                    [400., 410.], [250., 260.], [150., 160.], $
                    [95., 96.], [58., 59.], [36., 37.], [25., 26.]]
            
      FOR i = 12, N_ELEMENTS(all_energy)/2 - 1 DO BEGIN 
         
          energy = all_energy[*, i]

          plot_pa_spec_from_crib, sat, specie, inst, units_name, $
                                  energy, eff_table, $
                                  PaBin = 22.5, $
                                  recalc = 1, $
                                  BKG = 0, $
                                  COMBINE = 1

          average_tplot_variable, 1, '600'
          zlim, 1, 0.1, 100
       ;   e_min = STRCOMPRESS(STRING(energy(0), $
        ;                             FORMAT = '(i5.5)'), /REMOVE_ALL)
         ; e_max = STRCOMPRESS(STRING(energy(1),$
         ;                            FORMAT = '(i5.5)'), /REMOVE_ALL)
       ;   pa_name = 'PASPEC_EN'+ e_min +'_'+ e_max +'_SC4_UNDIFFFLUX_SP3_All'
stop    
      ENDFOR 
      
      zlim, '*', 0.1, 100.

      window, j
      tplot, 'PASPEC*'
      


print, 'start time: '+time_string(t_s)+'---end time: '+time_string(t_e)

stop
;--------------------------------------------------------------------
;decide type of the pitch angle distributions
;-------------------------------------------------------------------
e_min = STRCOMPRESS(STRING(energy(0, i), FORMAT = '(i5.5)'), /REMOVE_ALL)
e_max = STRCOMPRESS(STRING(energy(1, i), FORMAT = '(i5.5)'), /REMOVE_ALL)
pa_name = 'PASPEC_EN'+ e_min +'_'+ e_max +'_SC4_UNDIFFFLUX_SP3_All'
 
  time_pa = pa.x
  flux_pa = pa.y
  data_pa = pa.v
 
  zlim, pa_name, 0.1, 100.

;  tplot, pa_name
                  
  dd = time_string(min_dst(ii))

  fln = 'plots/storm_pa_' + $
        STRMID(dd, 0, 4) + STRMID(dd, 5, 2) + STRMID(dd, 8, 2) + '_' +$
        STRMID(dd, 11, 2) + STRMID(dd, 14, 2) + STRMID(dd, 17, 2) + '.ps'
  
  popen, fln, /land

  tplot, pa_name
             
  pclose

 ENDFOR 
ENDFOR
END
