PRO find_pa_peak, pa_name_counts, pa_name, pap_name, def_pap = def, pa_counts_line = pa_counts_line,flux_threshold=flux_threshold

if not keyword_set (flux_threshold) or n_elements(flux_threshold) ne 3 then flux_threshold=[10,15,18]
IF NOT keyword_set(pa_counts_line) THEN pa_counts_line = 9/88.
get_data, pa_name_counts, data = data
counts_pa = data.y(*, 0:7)

get_data, pa_name, data = data, dlim = dlim, lim = lim
time_pa = data.x
flux_pa = data.y(*, 0:7)
pa_pa = data.v(*, 0:7)

ntime = N_ELEMENTS(time_pa)
npa = N_ELEMENTS(pa_pa(0, *))

; calculate pitch angle peak cutting line 'def_pap' if keywords def is not set
average_flux_pa = TOTAL(flux_pa(*, *), 2)/npa
def_pap = DBLARR(ntime)
IF N_ELEMENTS(def) EQ 0 THEN BEGIN 
    tplot_names, '*beta', names = pb_name ; pb--plasma beta
    pb_name = pb_name(0)
    get_data, pb_name, data = pb
    time_pb = pb.x
    data_pb = pb.y
    data_pb = INTERPOL(data_pb, time_pb, time_pa)
    time_pb = time_pa
;    store_data, pb_name+'_INTERPOL', $
;                data = {x:time_pb, y:data_pb}, dlim = {psym:-7}
    FOR i = 0, ntime-1 DO BEGIN 
        FOR j = 0, npa-1 DO BEGIN               

            IF data_pb(i) LE 0.05 THEN BEGIN 
                def_pap(i) = total(flux_pa(i, *))*3/npa > flux_threshold(0)
            ENDIF 
            
            IF data_pb(i) GT 0.05 AND data_pb(i) LE 1 THEN BEGIN 
                def_pap(i) = total(flux_pa(i, *))*2/npa > flux_threshold(1)
            ENDIF 
                
            IF data_pb(i) GT 1 THEN BEGIN 
                def_pap(i) = total(flux_pa(i, *))*1.1/npa > flux_threshold(2)
            ENDIF 
      
        ENDFOR 
    ENDFOR 
ENDIF ELSE BEGIN    
    def_pap(*) = ABS(def)
ENDELSE 

; pitch angle peak have to be
; 1. local peak for unit flux
; 2. counts greater than a counts threshold  set before

flux_peak_pa = DBLARR(ntime, npa)
flux_peak_pa(*, *) = !VALUES.F_NAN

FOR i = 0, ntime-1 DO BEGIN 
    FOR j = 0, npa-1 DO BEGIN 
        IF flux_pa(i, j) GE flux_pa(i, j-1 > 0) AND $
          flux_pa(i, j) GE flux_pa(i, j+1 < (npa-1)) AND $
          flux_pa(i, j) GT def_pap(i) AND $
          counts_pa(i, j) GT pa_counts_line $
          THEN BEGIN
            flux_peak_pa(i, j) = flux_pa(i, j)
        ENDIF
    ENDFOR 
ENDFOR 

str = {x:time_pa, y:flux_peak_pa, v:pa_pa}
pap_name = pa_name+'_PAP'
store_data, pap_name, data = str, dlim = dlim, lim = lim
options, pap_name, 'ytitle', 'PAP'
zlim, pap_name, 0.1, 100
END
