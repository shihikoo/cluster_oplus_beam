
PRO plot_pa_spec_around_energy_peak, sat, specie, inst, units_name, $
                                     eff_table, $
                                     invar = invar, $
                                     outvar = outvar, $
                                     average_time = average_time, $
                                     n_range = n_range, $
                                     fixtherange = fixtherange, $
                                     recalc = recalc, $
                                     PaBin = PaBin

;----------------------------------------------------------
; check keywords and set the energybins range to nrange           
;---------------------------------------------------------
 IF N_ELEMENTS(invar) THEN BEGIN 
     get_data, invar, data = data

     time_ep = data.x
     energy_peak = data.y
    
     ncut = N_ELEMENTS(time_ep)
     
     IF N_ELEMENTS(n_range) THEN BEGIN 
         nrange = ABS(n_range)
     ENDIF ELSE BEGIN 
         nrange = 0             ;default energybins range
     ENDELSE 
;stop
     ; adjust the range if it is required and nrange is not 0 
     IF nrange GT 0 THEN BEGIN 
         IF N_ELEMENTS(fixtherange) EQ 0 THEN BEGIN 
             fixtherange = 0    ; undefined fixtherange equals fixtherange=0
         ENDIF 
         IF fixtherange EQ 0 THEN BEGIN ;if don't wanna fix the range


             nrange = ROUND((nenergybins_good/15.)*nrange)
         ENDIF ELSE BEGIN 
             nrange = ROUND(nrange)
         ENDELSE 
     ENDIF 
 ENDIF ELSE BEGIN 
     energy = [40, 40000]
     plot_pa_spec_from_crib, sat, specie, inst, units_name, $
                             energy, eff_table, $
                             PaBin = PaBin, $
                             recalc = recalc, $
                             BKG = 0, $
                             COMBINE = 1     
     PRINT, 'NO CHOSEN ENERGY PEAK, PLOT PA AT ALL ENERGY'
     Return
 ENDELSE 


;-----------------------------------------------------------------
;plot pa in the energy range and store pa specta data
;-----------------------------------------------------------------
 npa = CEIL(180/PaBin)
 energy_range = DBLARR(2, ncut)
 time_pa = time_ep
 flux_pa = DBLARR(ncut, npa)
 pa_pa = DBLARR(ncut, npa)

 ;store the energy range into a 2D array according to energybins range(nrange)
 IF nrange EQ 0 THEN BEGIN 
     energy_range(0, *) = energy_peak
     energy_range(1, *) = energy_peak
 ENDIF ELSE BEGIN 
     FOR jjj = 0, ncut-1 DO BEGIN 
         energybins = data.energybins         
         nenergybins = N_ELEMENTS(energybins)
         nenergybins_good = nenergybins ; nenergybins_good is those ebins which are good
         ebin = where(energybins EQ energy_peak(jjj))
         IF ebin NE -1 THEN BEGIN 
             energy_range(0, jjj) = energybins((ebin+nrange) < (nenergybins_good-1)) 
             energy_range(1, jjj) = energybins((ebin-nrange) > 0)
         ENDIF ELSE BEGIN 
             energy_range(*, jjj) = !VALUES.F_NAN
         ENDELSE 
     ENDFOR 
 ENDELSE 

 ; plot pa with energy_range
 IF N_ELEMENTS(average_time) EQ 0 THEN BEGIN 
     average_time = round((time_ep(ncut-1)-time_ep(0))/(ncut-1)) ;second
 ENDIF 

 n_mag_error = 0
 FOR jjj = 0, ncut-1 DO BEGIN 
     
     IF TOTAL(energy_range(*, jjj)) GT 0 THEN BEGIN 
        
         energy = energy_range(*, jjj)

         time = time_ep(jjj) - average_time/2
         timespan, time, average_time, /SECONDS

       ;------------------check mag data in this interval-------------------  
          get_cluster_mag_gse, sat, Btime, Bxyz

          IF Btime(0) NE -9999.9 THEN BEGIN 
       ;---------------------------------------------------------

              plot_pa_spec_from_crib, sat, specie, inst, units_name, $
                                      energy, eff_table, $
                                      PaBin = PaBin, $
                                      recalc = recalc, $
                                      BKG = 0, $
                                      COMBINE = 1     
           ;       stop
  
              e_min = STRCOMPRESS(STRING(energy(0), $
                                         FORMAT = '(i5.5)'), /REMOVE_ALL)
              e_max = STRCOMPRESS(STRING(energy(1), $
                                         FORMAT = '(i5.5)'), /REMOVE_ALL)
              sc_str = STRING(sat, FORMAT = '(i1.1)')

              s_pa_name = 'PASPEC_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                          + '_UNDIFFFLUX_SP3_All'
         
              get_data, s_pa_name, data = s_pa, dlim = dlim, lim = lim
              
              s_time_pa = s_pa.x
              s_flux_pa = s_pa.y
              s_pa_pa = s_pa.v
                   
       ;average flux data over average_time and save them into arrays
        
              y_data = TOTAL(s_flux_pa(*, 0:npa-1), 1)/N_ELEMENTS(s_flux_pa(*, 0))
              x_data = REFORM(s_pa_pa(0, 0:npa-1))

              flux_pa(jjj, *) = y_data
              pa_pa(jjj, * ) = x_data
          
      ;            PLOT, x_data, y_data, xlog = 0, ylog = 1,$ 
       ;                 xrange = [0, 180], yrange = [1e0, 1e4], psym = -2, $
        ;                xtitle = 'PITCH ANGLE', ytitle = 'Diff Flux', $
         ;               title = 'ENERGY RANGE:'+ e_min+'_'+ e_max +'!Ctimespan' $
          ;              + time_string(time)+' to '+time_string(time+average_time)
          ENDIF ELSE BEGIN 
              n_mag_error = n_mag_error + 1
              flux_pa(jjj, *) = !VALUES.F_NAN
              pa_pa(jjj, * ) = !VALUES.F_NAN
          ENDELSE 
      ENDIF ELSE BEGIN 
          flux_pa(jjj, *) = !VALUES.F_NAN
          pa_pa(jjj, * ) = !VALUES.F_NAN
      ENDELSE 
  ENDFOR 
  
 ;pos = STREGEX(invar, '_epcut')
 ;outvar = 'PA'+ STRMID(invar, 2, pos-2)
  pa_name = 'PASPEC_EN' + e_min + '_' + e_max + '_SC' + sc_str $
           + '_UNDIFFFLUX_SP3_All'
  str = {x: time_pa, y: flux_pa, v: pa_pa}
  store_data, pa_name, data = str, dlim = dlim, lim = lim

;--------------------------------------------------------
  ; use bx to decide tail-earth direction
  bx_name = 'MAG_SC4_B_xyz_gse_GSM_X'
  get_data, bx_name, data = data
  time_bx = data.x
  data_bx = data.y

  data_bx = INTERPOL(data_bx, time_bx, time_pa)
  time_bx = time_pa

  str = {x:time_bx, y:data_bx}
  store_data, bx_name+'_AVG', data = str, dlim = {psym:-7}
  
  time_et = time_pa
  pa_et = pa_pa
  flux_et = flux_pa
;stop
  FOR i = 0, ncut-1 DO BEGIN 
      IF data_bx(i) LE 0 THEN BEGIN 
          FOR ipa = 0, npa-1 DO BEGIN 
              flux_et(i, ipa) = flux_pa(i, npa-1-ipa)
          ENDFOR      
      ENDIF
  ENDFOR 
  index = where(flux_et EQ 0, ct)
  IF ct GT 0 THEN BEGIN 
      flux_et(index) = !VALUES.F_NAN ;consider 0 data as nan data for graphing
  ENDIF 
  et_name = 'DRSPEC_EN' + e_min + '_' + e_max + '_SC' + sc_str $
           + '_UNDIFFFLUX_SP3_All'
  str = {x:time_et, y:flux_et, v:pa_et}
  store_data, et_name ,data = str, dlim = dlim, lim = lim
  options, et_name, 'ytitle', ' '+STRING(energy(0), FORMAT = '(i5.2)')+'!C!CE---T'

  outvar = et_name
END
