;+
; PROCEDURE:  combine_pa_spec
;
; PURPOSE: For pith angle spectra to combine different products of the
;          same specie in one
;
; INPUT:
;        sat -> spacecraft number
;        specie -> specie number (0: H+, 1: He++, 2: He+, 3: O+)
;        units -> units
;        enrange -> Energy range
;        eff_table -> efficiency table
;
; KEYWORDS: type 'PASPEC' or 'PAFLUX'
;
; CREATED BY: C. Mouikis
;
; LAST MODIFICATION: 09/13/01
;
; MODIFICATION HISTORY:
;   09/10/01 - Adjusted for pa spec
;   09/13/01 - Plot type check introduced. This way there are less
;              variables that it tries to combine
;   20/11/01 - type keyword MF
;-
PRO combine_pa_spec, sat, specie, units, name_all, enrange, eff_table, $
                     type = type

IF NOT KEYWORD_SET(type) THEN type = 'PASPEC'

sat       = string(sat,    format = '(i1.1)')
specie    = string(specie, format = '(i1.1)')
eff_table = string(eff_table, format = '(i1.1)')

                                ; Read names of saves tplot structures
tplot_names, names = names

icomb = 0
FOR i1 = 0, n_elements(names)-1 DO BEGIN ; Loop through all structures
    
    plot_type   = strmid(names(i1), 0, 6)
    IF plot_type EQ type THEN BEGIN
        sc_numb     = strmid(names(i1), strpos(names(i1), '_SC')+3, 1)
        units_type  = strmid(names(i1), strpos(names(i1), '_UN')+3, $
                             strpos(names(i1), '_PR')-strpos(names(i1), 'UN') -2)
        specie_numb = strmid(names(i1), strpos(names(i1), '_SP')+3, 1)
        en_range = strmid(names(i1), strpos(names(i1), '_EN')+3, $
                          strpos(names(i1), '_SC')-strpos(names(i1), 'EN') -2)
        eff_t       = strmid(names(i1), strpos(names(i1), '_ET')+3, 1)
        
        IF units EQ units_type $ ; for units, sat, species do ...
          AND sat EQ sc_numb $
          AND specie EQ specie_numb $
          AND enrange EQ en_range $
          AND eff_table EQ eff_t THEN BEGIN
            
            icomb = icomb + 1
            
            get_data, names(i1), data = data, dlim = dlim ; get data
            
            IF icomb EQ 1 THEN BEGIN ; first time just set arrays
                xall = data.x
                yall = data.y
                vall = data.v
            ENDIF ELSE   BEGIN 
                nold = n_elements(yall(0, *)) ; number of energy bins for prev.
                nnew = n_elements(data.y(0, *)) ; number of energy bins for curr.
                IF nold EQ nnew THEN BEGIN ; if both have the same # of en. bins
                    xall = [xall, data.x]
                    yall = [yall, data.y]
                    vall = [vall, data.v]
                ENDIF ELSE BEGIN
                    IF nold EQ 16 AND  nnew EQ 8  THEN BEGIN ; old:16 & new:8
                        yall_new = dblarr(n_elements(data.x), 16)
                        vall_new = vall
                        FOR i2 = 0, 7 DO BEGIN
                            yall_new(*, i2*2)   = data.y(*, i2)
                            yall_new(*, i2*2+1) = data.y(*, i2)
                        ENDFOR 
                        xall = [xall, data.x]
                        yall = [yall, yall_new]
                        vall = [vall, vall_new]
                    ENDIF 
                    IF nold EQ 8 AND  nnew EQ 16  THEN BEGIN ; old:8 & new:16
                        yall_old = yall
                        vall_old = vall
                        pitchangle_new = vall_old(0, *)
                        yall = dblarr(n_elements(xall), 16)
                        vall = dblarr(n_elements(xall), 16)
                        FOR i2 = 0, 7 DO BEGIN
                            yall(*, i2*2)   = yall_old(*, i2)
                            yall(*, i2*2+1) = yall_old(*, i2)
                            vall(*, i2*2)   = pitchangle_new(i2)
                            vall(*, i2*2+1) = pitchangle_new(i2+1)
                        ENDFOR
                        xall = [xall, data.x]
                        yall = [yall, data.y]
                        vall = [vall, data.v]
                    ENDIF 
                ENDELSE   
            ENDELSE
        ENDIF
    ENDIF  
ENDFOR    

IF icomb GT 0 THEN BEGIN
    sort_ind = SORT(xall)
    
    xall_new = xall(sort_ind)
    yall_new = yall(sort_ind, *)
    vall_new = vall(sort_ind, *)
    
    datastr = {x:xall_new, y:yall_new, v:vall_new}
    labels = ''
    name_all = type + $
               '_EN' + en_range + $
               '_SC' + sat + $
               '_UN' + units + $
               '_SP' + specie + $
               '_All'
    store_data, name_all, data = datastr, dlim = dlim
ENDIF

END

