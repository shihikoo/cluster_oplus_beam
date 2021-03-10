
; Input: invar = spectrum name
; Output: stored energy peak data with the original name + 'epcut"
;         array all_energy with all the energy bins
;add drop_bad_bin keyword to decide if to drop the bad bin (bin 15)or not

PRO find_energy_peak_from_enspec, invar, outvar, drop_bad_bin = drop_bad_bin

 get_data, invar, data = data, dlim = dlim, lim = lim

 time_diff = data.x
 flux_diff = data.y
 energy_diff = data.v
 
 n_energybins = N_ELEMENTS(flux_diff(0, *))
 n_cut = N_ELEMENTS(flux_diff(*, 0))

 energybins = reform(energy_diff(0, 0:n_energybins-1)) 

 IF KEYWORD_SET(drop_bad_bin) THEN BEGIN 
     bad = 1                    ;the last bin was broken
 ENDIF ELSE BEGIN 
     bad = 0
 ENDELSE 
  
 n_energybins_good = n_energybins-bad

 energy_peak = FLTARR(n_cut)

 FOR iii = 0, n_cut - 1 DO BEGIN 
     index = WHERE(flux_diff(iii, 0:n_energybins_good-1) EQ $
                   MAX(flux_diff(iii, 0:n_energybins_good-1)) $
                   AND MAX(flux_diff(iii, 0:n_energybins_good-1)) GT 0, ct)
     IF ct NE 0 THEN BEGIN 
         energy_peak(iii) = energy_diff(iii, index(0))
     ENDIF ELSE BEGIN 
         energy_peak(iii) = !VALUES.F_NAN
     ENDELSE 
 ENDFOR

 outvar = invar+'_epcut'
 str = {x: time_diff, y: energy_peak, energybins: energybins }
 store_data, outvar, data = str, dlim = {psym: -2}

END
