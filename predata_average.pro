FUNCTION predata_average, imf, time
pretime = 3600.

ntime = N_ELEMENTS(time)
imf_pre = FLTARR(ntime)
FOR itime = 0l, ntime-1 DO BEGIN 
    index = where(time GT (time(itime)-pretime) AND time LT time(itime), ct)
    IF ct GE 6 THEN BEGIN 
        IF  TOTAL(ABS(imf(index)) GE 0) GE 6 THEN BEGIN 
            imf_pre(itime) = TOTAL(imf(index), /NAN)/ct
        ENDIF ELSE imf_pre(itime) = !VALUES.F_NAN
    ENDIF ELSE imf_pre(itime) = !VALUES.F_NAN
ENDFOR 
return, imf_pre
END 
