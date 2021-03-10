PRO plot_pa_spec_around_energy_peak, sat, specie, inst, units_name, $
  eff_table, $
  invar = invar, $
  outvar = outvar, $
  average_time = average_time, $
  start_time = start_time, $
  END_time = END_time, $
  n_range = n_range, $
  PaBin = PaBin

;----------------------------------------------------------
; check keywords and set the energybins range to nrange           
;---------------------------------------------------------
CASE units_name OF 
    'DIFF FLUX': units_str = 'DIFFFLUX' 
    'Counts':units_str = 'COUNTS'  
    'EFLUX': units_str = 'EFLUX'
ENDCASE 

IF N_ELEMENTS(invar) THEN BEGIN 
    get_data, invar, data = data

    time_ep = data.x
    energy_peak = data.y
    
    ncut = N_ELEMENTS(energy_peak)
    ntime = floor((end_time-start_time)/average_time)

    energybins = data.energybins         
    nenergybins = N_ELEMENTS(energybins)
    nenergybins_good = nenergybins - 1 ;the last energybin is borken

    IF KEYWORD_SET(n_range) THEN nrange = ROUND((nenergybins_good/15.)*ABS(n_range)) $
    ELSE nrange = 0             ;default energybins range    
ENDIF ELSE BEGIN 
    energy = [40., 40000.]
    plot_pa_spec_from_crib, sat, specie, inst, units_name, $
      energy, eff_table, $
      PaBin = PaBin, $
      recalc = recalc, $
      BKG = 0, $
      COMBINE = 1     
    PRINT, 'NO CHOSEN ENERGY PEAK, PLOT PA AT ALL ENERGY'
    Return
ENDELSE 

get_timespan, interval
;-----------------------------------------------------------------
;plot pa in the energy range and store pa specta data
;-----------------------------------------------------------------
energy_range = DBLARR(2, ncut)
time_pa = time_ep
flux_pa = DBLARR(ncut, 8)
pa_pa = DBLARR(ncut, 8)

;store the energy range into a 2D array according to energybins range(nrange)
IF nrange EQ 0 THEN BEGIN 
    energy_range(0, *) = energy_peak
    energy_range(1, *) = energy_peak
ENDIF ELSE BEGIN 
    FOR jjj = 0, ncut-1 DO BEGIN 
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
IF N_ELEMENTS(average_time) EQ 0 THEN average_time = 300 ;default average time: 300s
n_mag_error = 0

FOR jjj = 0, ntime-1 DO BEGIN 

    time = start_time+jjj*average_time
    index = where(time_ep GT time AND time_ep LE time+average_time, ct)
    loc = index(0)

    IF ct EQ 1 THEN BEGIN 
        IF TOTAL(energy_range(*, loc)) GT 0 THEN BEGIN 
            energy = energy_range(*, loc)
            timespan, time, average_time, /SECONDS

;------------------check mag data in this interval-------------------  
            get_cluster_mag_gse, sat, Btime, Bxyz
            IF N_ELEMENTS(Btime) GE 2  THEN BEGIN 
;---------------------------------------------------------
                IF inst EQ 0 OR inst EQ 2 THEN BEGIN 
                    plot_pa_spec_from_crib, sat, specie, inst, units_name, $
                      energy, eff_table, $
                      PaBin = PaBin, $
                      recalc = 1, $
                      BKG = 0, $
                      COMBINE = 1     
                ENDIF ELSE BEGIN 
                    plot_hia_pa_spec_from_crib, sat, specie, inst, units_name, $
                      energy, eff_table, $
                      PaBin = PaBin, $
                      recalc = 1, $
                      BKG = 0, $
                      COMBINE = 1     
                ENDELSE 

                e_min = STRCOMPRESS(STRING(energy(0), $
                                           FORMAT = '(i5.5)'), /REMOVE_ALL)
                e_max = STRCOMPRESS(STRING(energy(1), $
                                           FORMAT = '(i5.5)'), /REMOVE_ALL)
                sc_str = STRING(sat, FORMAT = '(i1.1)')
                
                s_pa_name = 'PASPEC_EN' + e_min + '_' + e_max $
                            + '_SC' + sc_str + '_UN'+units_str+'_SP3_All'

                tplot_names, s_pa_name, names = names
                IF names(0) NE '' THEN BEGIN 
                    get_data, s_pa_name, data = s_pa, dlim = dlim, lim = lim
                    store_data, delete = s_pa_name

                    s_time_pa = s_pa.x
                    s_flux_pa = s_pa.y
                    s_pa_pa = s_pa.v                   
;average flux data over average_time and save them into arrays 
                    IF units_name EQ 'DIFF FLUX' OR units_name EQ 'EFLUX' THEN  y_data = TOTAL(s_flux_pa(*, 0:7), 1)/N_ELEMENTS(s_flux_pa(*, 0))
                    IF units_name EQ 'Counts' THEN   y_data = TOTAL(s_flux_pa(*, 0:7), 1)

                    x_data = REFORM(s_pa_pa(0, 0:7))
                    flux_pa(loc, *) = y_data
                    pa_pa(loc, * ) = x_data

                                ;            PLOT, x_data, y_data, xlog = 0, ylog = 1,$ 
                                ;                 xrange = [0, 180], yrange = [1e0, 1e4], psym = -2, $
                                ;                xtitle = 'PITCH ANGLE', ytitle = 'Diff Flux', $
                                ;               title = 'ENERGY RANGE:'+ e_min+'_'+ e_max +'!Ctimespan' $
                                ;              + time_string(time)+' to '+time_string(time+average_time)
                ENDIF ELSE BEGIN
                    flux_pa(loc, *) = !VALUES.F_NAN
                    pa_pa(loc, * ) = !VALUES.F_NAN
                ENDELSE 
            ENDIF  ELSE BEGIN 
                n_mag_error = n_mag_error + 1
                flux_pa(loc, *) = !VALUES.F_NAN
                pa_pa(loc, * ) = !VALUES.F_NAN
            ENDELSE 
        ENDIF     ELSE BEGIN 
            flux_pa(loc, *) = !VALUES.F_NAN
            pa_pa(loc, * ) = !VALUES.F_NAN
        ENDELSE 
    ENDIF 
    IF ct GE  2 THEN stop
ENDFOR  

pos = STREGEX(invar, '_epcut') 
outvar = 'PA'+ STRMID(invar, 2, 25)+units_str+STRMID(invar, pos-19, 19)
str = {x: time_pa, y: flux_pa, v: pa_pa}
store_data, outvar, data = str, dlim = dlim, lim = lim
options,  outvar, 'ytitle', 'Pitch Angle!C!CVarious En'
zlim, outvar, 0.1, 100

;set the timespan as before 
timespan, interval(0), interval(1)-interval(0), /SECONDS
END
