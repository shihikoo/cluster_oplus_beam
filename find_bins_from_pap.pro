PRO find_bins_from_pap, sat, specie, inst, units_name, eff_table, $
  average_time, pap_name, bins_name, start_time = start_time, END_time = END_time, $
  beam_angle_range = beam_angle_range

IF NOT keyword_set(beam_angle_range) THEN beam_angle_range = 11.25
IF NOT keyword_set(bins_name) THEN BEGIN 
    pos = STREGEX(pap_name, 'AVG') 
    bins_name = 'GLOBE'+STRMID(pap_name, 6, pos)+'_bin_range'
ENDIF 
;get info from pitch angle peak file, so the time can be set
get_data, pap_name, data = data
time_pap = data.x
pap = data.v
pap(where( ~FINITE(data.y))) = !VALUES.F_NAN
ntime = floor((end_time-start_time)/average_time) ;number of cuts during all tiem
n_avg = N_ELEMENTS(time_pap)    ;number of averaged data

th_ph_index = INTARR(n_avg, 88)

;get the original time span  
get_timespan, interval

FOR it = 0, ntime-1 DO BEGIN 
;load mag_theta, mag_phi data in codif coordinates 
    time = start_time + it*average_time
    index = where(time_pap GT time AND time_pap LE time+average_time, ct)
    pos = index(0)
    IF ct EQ 1 THEN BEGIN 
; calculate the magnetic field data
        index_pap = where(pap(pos, *) GE 0, ct_pap)
        IF ct_pap GT 0 THEN BEGIN 
            timespan, time, average_time, /SECONDS
            COMMON get_error, get_err_no, get_err_msg, default_verbose
            CASE inst OF 
               0: prod = [17, 18, 47, 49]
               ELSE: BEGIN  
                   print, '0 for CODIF'
                   stop
               ENDELSE  
           ENDCASE

           jj = 0
           found = 0
           WHILE found EQ 0 AND jj LT n_elements(prod) DO BEGIN
              
; dat structure is needed for use 'get_theta_phi' program 
               IF inst EQ 0 THEN dat = call_function('get_cis_cod_data', prod(jj), specie = specie, sat)  
; GSE -> CODIF coordinate transformation

               IF get_err_no EQ 0 THEN BEGIN ; check if data were found for time interval
                   found = 1
                   mag_theta_s = 0.
                   mag_phi_s = 180.
                   inst = inst
                   get_theta_phi, sat, dat, mag_theta_s, mag_phi_s, inst ;get_theta_phi gives mag_theta and mag_phi in instrument coodinates
               ENDIF
               jj = jj + 1
           ENDWHILE

           theta = reform(dat.theta(0, *))
           phi = reform(dat.phi(0, *))

           Bx_norm = cos(!DTOR*mag_theta_s)*sin(!DTOR*mag_phi_s)
           By_norm = cos(!DTOR*mag_theta_s)*cos(!DTOR*mag_phi_s)
           Bz_norm = sin(!DTOR*mag_theta_s)
           pa_cal = fltarr(88)

           Pa_cal = acos(Bx_norm*cos(!DTOR*theta)*sin(!DTOR*phi)$
                         +By_norm*cos(!DTOR*theta)*cos(!DTOR*phi) $
                         +Bz_norm*sin(!DTOR*theta))/!DTOR
        
           FOR i_pap = 0, ct_pap-1 DO BEGIN 
               IF pap(pos, index_pap(i_pap)) NE 90 THEN BEGIN  
                   pn = (90-pap(pos, index_pap(i_pap)))/ABS(90-pap(pos, index_pap(i_pap))) ;=1for pap<90, =-1for pap>90
                   index =  where ((pa_cal - pap(pos, index_pap(i_pap)))*pn LE beam_angle_range) 
                   IF index(0) NE -1 THEN th_ph_index(pos, index) = 1
               ENDIF ELSE BEGIN 
                   th_ph_index(pos, *) = 1
               ENDELSE 
           ENDFOR
       ENDIF 
   ENDIF  ELSE BEGIN 
       IF ct GT  1 THEN stop
   ENDELSE  
ENDFOR     

store_data, bins_name, data = {x:time_pap, y:th_ph_index}

tplot_names, 'B_xyz_*', names = names
store_data, delete = names

timespan, interval(0), interval(1)-interval(0), /seconds
END 
