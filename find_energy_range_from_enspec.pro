; Input: enspec_name : energy spectrum name
; Output: epcut_name:stored energy peak data with the original name + '_epcut'
;         erange_name: stored energy range data with the original name
;                      +'_erange'
;         array all_energy with all the energy bins
;
; by Jing Liao 11/01/2007

PRO find_energy_range_from_enspec, enspec_name, epcut_name, erange_name

 get_data, enspec_name, data = data, dlim = dlim, lim = lim

 time = data.x
 flux = data.y
 energy = data.v
 
 n_energybins = N_ELEMENTS(flux(0, *))
 n_cut = N_ELEMENTS(flux(*, 0))

 energybins = reform(energy(0, 0:n_energybins-1)) 

 n_energybins_good = n_elements(energy(0, where(energy(0, *) GT 35))) 
;make sure energy is higer than 35eV to avoid the bad energy bin

 energy_peak = FLTARR(n_cut)
 energy_range = FLTARR(n_cut, 2)

 FOR iii = 0, n_cut - 1 DO BEGIN 
     index = WHERE(flux(iii, 0:n_energybins_good-1) EQ $
                   MAX(flux(iii, 0:n_energybins_good-1)) $
                   AND MAX(flux(iii, 0:n_energybins_good-1)) GT 0, ct)
     IF ct NE 0 THEN BEGIN 
         energy_peak(iii) = energy(iii, index(0))
         energy_range(iii, *) = [energy(iii, (index(0)+1) < (n_energybins_good-1)),$
                                 energy(iii, (index(0)-1) > 0)]

         i_f = index(0)-round(n_energybins/16.) > 0
         WHILE i_f GT 0 DO BEGIN
             IF flux(iii, index(0))/flux(iii, i_f) LE 10. THEN BEGIN  
                 energy_range(iii, 1) = energy(iii, i_f)
                 i_f = i_f-1
             ENDIF  ELSE BEGIN  
                 i_f = -1
             ENDELSE  
         ENDWHILE  
    
         i_f = index(0)+round(n_energybins/16.) < (n_energybins_good-1)
         WHILE i_f LT n_energybins_good-1  DO BEGIN
             IF flux(iii, index(0))/flux(iii, i_f) LE 10. THEN BEGIN  
                 energy_range(iii, 0) = energy(iii, i_f)
                 i_f = i_f+1
             ENDIF ELSE BEGIN 
                 i_f = 100
             ENDELSE 
         ENDWHILE 
     ENDIF  ELSE BEGIN 
         energy_peak(iii) = !VALUES.F_NAN
         energy_range(iii, *) = energy_peak(iii)
     ENDELSE 
 ENDFOR 
 epcut_name = enspec_name+'_epcut'
 str = {x: time, y: energy_peak, energybins: energybins }
 store_data, epcut_name, data = str, dlim = {psym: -3}
     
 erange_name = enspec_name+'_erange'
 str = {x: time, y: energy_range, energybins: energybins }
 store_data, erange_name, data = str, dlim = {psym: -3}
END
