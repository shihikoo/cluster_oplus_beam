PRO find_pa_peak, invar, outvar, def_pap = def

get_data, invar, data = data, dlim = dlim, lim = lim
time_pa = data.x
flux_pa = data.y(*, 0:7)
pa_pa = data.v(*, 0:7)

ntime = N_ELEMENTS(time_pa)
npa = N_ELEMENTS(pa_pa(0, *));have to be even number

bx_name = 'MAG_SC4_B_xyz_gse_GSM_X'
get_data,bx_name, data = data
time_bx = data.x
data_bx = data.y
;stop
data_bx = INTERPOL(data_bx, time_bx, time_pa)
time_bx = time_pa

str = {x:time_bx, y:data_bx}
store_data, bx_name+'_AVG', data = str, dlim = {psym:-7}

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
    store_data, pb_name+'_AVG', data = {x:time_pb, y:data_pb}, dlim = {psym:-7}
    FOR i = 0, ntime-1 DO BEGIN 
        FOR j = 0, npa-1 DO BEGIN               

            IF data_pb(i) LE 0.05 THEN BEGIN 
                def_pap(i) = total(flux_pa(i, *))*3/npa > 10
            ENDIF 
            
            IF data_pb(i) GT 0.05 AND data_pb(i) LE 1 THEN BEGIN 
                def_pap(i) = total(flux_pa(i, *))*2/npa > 15
            ENDIF 
                
            IF data_pb(i) GT 1 THEN BEGIN 
                def_pap(i) = total(flux_pa(i, *))*1.1/npa > 18
            ENDIF 
              
        ENDFOR 
    ENDFOR 
ENDIF ELSE BEGIN    
    def_pap(*) = ABS(def)
ENDELSE 

;time_beam = DBLARR(300)
;judge_beam = INTARR(300)

time_peak = time_pa
flux_peak = DBLARR(ntime, npa)
flux_peak(*, *) = !VALUES.F_NAN
pa_peak = pa_pa


; bx > 0 --> 0 --> earthward(-1)
;            180--> tailward(1)
; bx < 0 --> 0 --> tailward(1)
;            180--> earthward(-1)


;npeak = 0
FOR i = 0, ntime-1 DO BEGIN 
  ;  np = 0

    IF flux_pa(i, 0) GE flux_pa(i, 1) $ 
      AND flux_pa(i, 0) GT def_pap(i) THEN BEGIN 
      ;  time_beam(npeak) = time_pa(i)
        IF data_bx(i) GE 0 THEN BEGIN  
            flux_peak(i, 0) = flux_pa(i, 0)
          ;  judge_beam(npeak) = pa_pa(i, 0)
        ENDIF ELSE BEGIN 
            flux_peak(i, npa-1) = flux_pa(i, 0)
          ;  judge_beam(npeak) = pa_pa(i, npa-1)
        ENDELSE 
     ;   npeak = npeak+1        
     ;   np = np+1

    ENDIF 
  
    FOR j = 1, npa/2-1 DO BEGIN 
        IF flux_pa(i, j) GE flux_pa(i, j-1) AND $
          flux_pa(i, j) GE flux_pa(i, j+1) AND $
          flux_pa(i, j) GT def_pap(i) THEN BEGIN

         ;   time_beam(npeak) = time_pa(i)
            IF data_bx(i) GE 0 THEN BEGIN  
                flux_peak(i, j) = flux_pa(i, j)
              ;  judge_beam(npeak) = pa_pa(i, j)
            ENDIF ELSE BEGIN 
                flux_peak(i, npa-1-j) = flux_pa(i, j)
              ;  judge_beam(npeak) = pa_pa(i, npa-1-j)
            ENDELSE 
        ;    npeak = npeak+1 
         ;   np = np+1

        ENDIF
    ENDFOR 

    FOR j = npa/2, npa-2 DO BEGIN 
        IF flux_pa(i, j) GE flux_pa(i, j-1) AND $
          flux_pa(i, j) GE flux_pa(i, j+1) AND $
          flux_pa(i, j) GT def_pap(i) THEN BEGIN
            
        ;    time_beam(npeak) = time_pa(i)
            IF data_bx(i) GE 0 THEN BEGIN  
                flux_peak(i, j) = flux_pa(i, j)
           ;     judge_beam(npeak) = pa_pa(i, j)
            ENDIF ELSE BEGIN 
                flux_peak(i, npa-1-j) = flux_pa(i, j)
            ;    judge_beam(npeak) = pa_pa(i, npa-1-j)
            ENDELSE 
          
        ;    npeak = npeak+1 
        ;    np = np+1
        ENDIF
    ENDFOR 

    IF flux_pa(i, npa-1) GE flux_pa(i, npa-2) $ 
      AND flux_pa(i, npa-1) GT def_pap(i) THEN BEGIN
       
     ;   time_beam(npeak) = time_pa(i)
        
        IF data_bx(i) GE 0 THEN BEGIN  
            flux_peak(i, npa-1) = flux_pa(i, npa-1)
         ;   judge_beam(npeak) = pa_pa(i, npa-1)
        ENDIF ELSE BEGIN 
            flux_peak(i, 0) = flux_pa(i, npa-1)
         ;   judge_beam(npeak) = pa_pa(i, 0)

        ENDELSE 

      ;  npeak = npeak+1 
      ;  np = np+1
    ENDIF 
ENDFOR 

;IF npeak GT 0 THEN BEGIN 
;judge_beam = judge_beam(0:npeak-1)
;time_beam = time_beam(0:npeak-1)
;ENDIF

str = {x:time_peak, y:flux_peak, v:pa_peak}
pap_name = invar+'_pap'
store_data,pap_name, data = str, dlim = dlim, lim = lim
zlim, pap_name, 0.1, 100

outvar = pap_name

;str = {x:time_beam, y:judge_beam}
;judge_name = invar+'_BEAMJUDGEMENT' 
;store_data, judge_name, data = str, dlim = {psym:2}

;stop
END
