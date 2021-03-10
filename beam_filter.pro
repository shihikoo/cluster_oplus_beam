PRO beam_filter, sc, pap_name, epcut_name, pap_beam_et_name, epcut_beam_name

sc_str = STRING(sc, FORMAT = '(i1.1)')
;-----------------------------------
;load input data and basic settings
;-------------------------------------
get_data, pap_name, data = data, dlim = dlim, lim = lim
time_avg= data.x
flux_pap = data.y
pa_pap = data.v

get_data, epcut_name, data = data
energy_peak = data.y
nenergy = N_ELEMENTS(data.energybins)
energybins = data.energybins(0:nenergy-1)

ntime = N_ELEMENTS(time_avg)
npa = N_ELEMENTS(pa_pap(0, *)) 
n_e = 2                         ;accepted energy bin difference 
;-------------------------------------------------------------------
;Use x gse divide data into two regions: near earth and tail region
;Change pitch angle direction data to earth-tail direction data with 
;Bx (tail region) and z_gsm (near earth region)
;------------------------------------------------------------------
;load data 
bx_name = 'MAG_SC'+sc_str+'_B_xyz_gse_X'
get_data,bx_name, data = data
time_bx = data.x
data_bx = data.y
data_bx = INTERPOL(data_bx, time_bx, time_avg)
time_bx = time_avg

x_gse_name = 'EPH_SC'+sc_str+'_GSE_X'
get_data, x_gse_name, data = data
time_x_gse = data.x
data_x_gse = data.y
data_x_gse = INTERPOL(data_x_gse, time_x_gse, time_avg)
time_x_gse = time_avg

z_gsm_name = 'EPH_SC'+sc_str+'_GSM_Z'
get_data, z_gsm_name, data = data
time_z_gsm = data.x
data_z_gsm = data.y
data_z_gsm = INTERPOL(data_z_gsm, time_z_gsm, time_avg)
time_z_gsm = time_avg

; For Tail region (x gse < -1)
; bx > 0 --> 0 --> earthward,  180--> tailward
; bx < 0 --> 0 --> tailward, 180--> earthward
; For near earth region (x gse >= -1)
; z_gsm > 0(north) 0 --> earthward, 180--> outward
; z_gsm < 0(south) 0 --> outward, 180--> earthward

flux_pap_et = flux_pap
FOR i = 0, ntime-1 DO  BEGIN 
    IF data_x_gse(i) LT -1 THEN BEGIN 
        IF data_bx(i) LT 0 then flux_pap_et(i, *) = REVERSE(flux_pap(i, *), 2)
    ENDIF ELSE BEGIN 
        IF data_z_gsm(i) LT 0 then flux_pap_et(i, *) = REVERSE(flux_pap(i, *), 2)
    ENDELSE 
ENDFOR 

str = {x:time_avg, y:flux_pap_et, v:pa_pap}
pos = STREGEX(pap_name, 'PAP')
pap_et_name = STRMID(pap_name, 0, pos+3)+'_ET'
store_data, pap_et_name, data = str, dlim = dlim, lim = lim
zlim, pap_et_name, 0.1, 100
options, pap_et_name, 'ytitle', 'PAP!C!C E------T'
;----------------------------
;check for the beam
;----------------------------
flux_beam = flux_pap_et
epcut_beam = energy_peak

;check energy (longer than 2 times of average_time)
index_ep = INTARR(ntime)
FOR k = 0, ntime-1 DO  index_ep(k) = where(energybins EQ epcut_beam(k)) 
ind=where(index_ep eq -1,ct)
;if ct gt 0 then index_ep(ind)=-1 ; -100 as a flag for no energy peak or no beam
ep_beam=index_ep
ep_beam(*)=!values.f_nan
IF ntime GT 2 THEN BEGIN
    FOR i = 1, ntime-2 DO BEGIN
        if index_ep(i) ne -1 then begin 
            l=ABS(index_ep(i)-index_ep(i-1)) le n_e and index_ep(i-1) ne -1
            r=ABS(index_ep(i)-index_ep(i+1)) le n_e and index_ep(i+1) ne -1
            ep_beam(i)=l+r
        endif else ep_beam(i)=0
    ENDFOR
; add the edge point back
    FOR i = 1, ntime-2 DO BEGIN 
        IF ep_beam(i) eq 2 THEN BEGIN 
            IF ep_beam(i-1) eq 1 THEN ep_beam(i-1)=ep_beam(i-1)+1
            IF ep_beam(i+1) eq 1 THEN ep_beam(i+1)=ep_beam(i+1)+1
        ENDIF
    ENDFOR
endif  else ep_beam(*)=0

; only the points with both sides having beam remain
gap = where(ep_beam lt 2, ct)
;put data with no beam into Nan
IF ct GT 0 THEN BEGIN
    index_ep(gap)=-1
    flux_beam(gap, *) = !VALUES.F_NAN 
ENDIF  

;check flow direction  (longer than 3 times of average_time)
nbeam = intarr(ntime, npa)
FOR i = 1, ntime-2 DO BEGIN   
    FOR j = 0, npa-1 DO BEGIN 
        IF flux_beam(i, j) GT 0 THEN BEGIN 

            l1 = flux_beam(i-1, j-1 > 0)GT 0
            l2 = flux_beam(i-1, j)GT 0
            l3 = flux_beam(i-1, j+1 < (npa-1))GT 0

            r1 = flux_beam(i+1, j-1 > 0)GT 0
            r2 = flux_beam(i+1, j)GT 0
            r3 = flux_beam(i+1, j+1 < (npa-1))GT 0
                                ; for the case pitch angle = 0 / 180
            IF j EQ 0 THEN BEGIN 
                l1 = 0 & r1 = 0 
            ENDIF 
            IF j EQ npa-1 THEN BEGIN 
                l3 = 0 & r3 = 0 
            ENDIF  
            nbeam(i, j) = ((l1+l2+l3) GT 0) + ((r1+r2+r3) GT 0)
        ENDIF 
    ENDFOR 
ENDFOR

; add the edge point back as 2
FOR i = 1, ntime-2 DO BEGIN 
    FOR j = 0, npa-1 DO BEGIN 
        IF nbeam(i, j) GE 2 THEN BEGIN 
            IF total(nbeam(i-1, (j-1 > 0):(j+1 < (npa-1))), 2) EQ 1 THEN BEGIN 
                nbeam(i-1, (j-1 > 0):(j+1 < (npa-1))) =  $
                  nbeam(i-1, (j-1 > 0):(j+1 < (npa-1))) + 1
            ENDIF
            IF total(nbeam(i+1, (j-1 > 0):(j+1 < (npa-1))), 2) EQ 1 THEN BEGIN 
                nbeam(i+1, (j-1 > 0):(j+1 < (npa-1))) = $
                  nbeam(i+1, (j-1 > 0):(j+1 < (npa-1)))+ 1
            ENDIF
        ENDIF 
    ENDFOR 
ENDFOR 
; so anything here marked as 0 or 1 is not beam
ind=where(nbeam le 1, ct)
if ct gt 0 then begin 
    nbeam(ind) = 0
    flux_beam(ind) = !VALUES.F_NAN
endif 

;record all non-beam filtered with direction also into energy peak string
loc = where(total(nbeam(*, *), 2,/nan) EQ 0, ct)
IF ct GT 0 THEN BEGIN 
    index_ep(loc) = -1
    epcut_beam(where(index_ep EQ -1)) = !VALUES.F_NAN
ENDIF 

;store data into input_name +'_ET_beam'
str = {x:time_avg, y:flux_beam, v:pa_pap}
pos = STREGEX(pap_name, 'PAP')
pap_beam_et_name = STRMID(pap_name, 0, pos+3)+'_ET_beam'
store_data,pap_beam_et_name , data = str, dlim = dlim, lim = lim
options, pap_beam_et_name, 'ytitle', 'BEAM!C!C E------T'

str = {x:time_avg, y:epcut_beam, energybins:energybins}
epcut_beam_name = epcut_name+'_beam'
store_data, epcut_beam_name, data = str, dlim = {psym:0}
; calculate data back into pa direction 
;flux_beam_pa = flux_beam
;FOR i = 0, ntime-1 DO  BEGIN 
;    IF data_x_gse(i) LT -1 THEN BEGIN 
;        IF data_bx(i) LT 0 then flux_beam_pa(i, *) = REVERSE(flux_beam(i, *), 2)
;    ENDIF ELSE BEGIN 
;        IF data_z_gsm(i) LT 0 then flux_beam_pa(i, *) = REVERSE(flux_beam(i, *), 2)
;    ENDELSE 
;ENDFOR 
END
