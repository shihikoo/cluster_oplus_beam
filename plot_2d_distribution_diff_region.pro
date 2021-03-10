FUNCTION plot_2d_distribution_diff_region, sta_name, xrange = xrange, para_name = para_name, region_name = region_name, yrange = yrange, path = path, ps = ps

IF NOT keyword_set(yrange) THEN yrange = [0.01, 1]
n_sta = N_ELEMENTS(sta_name) 
IF n_sta LT 2 THEN BEGIN 
    print,  'Statistical data from more than 2 regions are required'
    return, 0
ENDIF 
IF NOT keyword_set(para_name)THEN BEGIN 
    pos = STREGEX(sta_name(0), '_x_')
    para_name = STRMID(sta_name(0), 0, pos)
ENDIF 

datax_all = FLTARR(3000, n_sta)
n_data = INTARR(n_sta)
sta_all = FLTARR(3000, 6, n_sta)
ratio_storm_all = FLTARR(3000, n_sta)
ratio_nonstorm_all = FLTARR(3000, n_sta)
error_storm_all = FLTARR(3000, n_sta)
error_nonstorm_all = FLTARR(3000, n_sta)
region_name = STRARR(n_sta)
phase_name_ = ['all_data',  'beam_for_all_data', 'nonstorm', 'beam_for_nonstorm_data', $
              'storm', 'beam_for_storm_data']

FOR i = 0, n_sta-1 DO BEGIN 
    get_data, sta_name(i), data = data
    n_data(i) = N_ELEMENTS(data.x)
    IF n_data(i) GT 3000 THEN BEGIN 
        print,  'Too many data, need to change the code'
        return, 0
    ENDIF 
    datax_all(0:n_data(i)-1, i) = data.x
    sta_all(0:n_data(i)-1, *, i) = REFORM(data.sta)
    pos1 = STREGEX(sta_name(i), '_x_')
    pos2 = STREGEX(sta_name(i), 'data_distribution') 
    region_name(i) = STRMID(sta_name(i), pos1+3, pos2-pos1-4)
ENDFOR 

;---------Ratio
FOR i = 0, n_sta-1 DO BEGIN 
    sta = sta_all(0:n_data(i)-1, *, i)
    datax = datax_all(0:n_data(i)-1, i)

    error = sqrt(sta)
    error_nonstorm = sqrt(error(*, 3)^2/sta(*, 2)^2+ $
                          sta(*, 3)^2*error(*, 2)^2/sta(*, 2)^4)
    error_storm = sqrt(error(*, 5)^2/sta(*, 4)^2+ $
                       sta(*, 5)^2*error(*, 4)^2/sta(*, 4)^4)
    index_nonstorm = where(sta(*, 2) EQ 0, ct)
    ratio_nonstorm = sta(*, 3)/sta(*, 2)
    IF ct GT 0 THEN  BEGIN 
        ratio_nonstorm(index_nonstorm) = !VALUES.F_NAN
        error_nonstorm(index_nonstorm) = !VALUES.F_NAN
    ENDIF 
    index_storm = where(sta(*, 4) EQ 0, ct)
    ratio_storm = sta(*, 5)/sta(*, 4)
    IF ct GT 0 THEN BEGIN 
        ratio_storm(index_storm) = !VALUES.F_NAN
        error_storm(index_storm) = !VALUES.F_NAN
    ENDIF 
    index = where(sta(*, 2) EQ 0 AND sta(*, 4) EQ 0, ct)
    IF ct GT 0 THEN  datax(index) = !VALUES.F_NAN

    datax_all(0:n_data(i)-1, i) = datax
    ratio_storm_all(0:n_data(i)-1, i) = ratio_storm
    ratio_nonstorm_all(0:n_data(i)-1, i) = ratio_nonstorm
    error_storm_all(0:n_data(i)-1, i) = error_storm
    error_nonstorm_all(0:n_data(i)-1, i) = error_nonstorm
ENDFOR 

;---write it into ps file
;IF keyword_set(ratio_region) THEN BEGIN 
IF NOT keyword_set(xrange) THEN xrnage = [min(datax_all)-1, max(datax_all)+1]

IF keyword_set(ps) THEN popen, path+para_name+'_ratio_nonstorm.ps', /land ELSE window, /free
plot, datax_all(*, 0), ratio_nonstorm(*, 0), xrange = xrange, yrange = yrange $
      , xstyle = 1, title = 'Nonstorm Beam Ratio Distribution for ' + para_name, $
      xtitle = para_name, ytitle = 'Beam Ratio', ylog = 0, $
      position = [0.15, 0.15, 0.95, 0.85], charsize = 2, charthick = 3, /nodata
FOR i = 0, n_sta-1 DO BEGIN 
    oplot, datax_all(0:n_data(i)-1, i), ratio_nonstorm_all(0:n_data(i)-1, i), psym = -2, col = i+1, thick = 4
    errplot, datax_all(0:n_data(i)-1, i), $
             ratio_nonstorm_all(0:n_data(i)-1, i)-error_nonstorm_all(0:n_data(i)-1, i), $
             ratio_nonstorm_all(0:n_data(i)-1, i)+error_nonstorm_all(0:n_data(i)-1, i), col = i+1, thick = 4
    xyouts, 1, 0.7+0.1*i, region_name(i), col = i+1, charsize = 3, charthick = 3
ENDFOR 
pclose
;stop
; window,/free
IF keyword_set(ps) THEN popen, path+para_name+'_ratio_storm.ps', /land ELSE window, /free
plot, datax_all(*, 0), ratio_storm(*, 0), xrange = xrange, yrange = yrange $
      , xstyle = 1, title = 'Storm Beam Ratio Distribution for ' + para_name $
      , xtitle = para_name, ytitle = 'Beam Ratio', ylog = 0, $
      position = [0.15, 0.15, 0.95, 0.85], charsize = 2, charthick = 3, /nodata
FOR i = 0, n_sta-1 DO BEGIN 
    oplot, datax_all(0:n_data(i)-1, i), ratio_storm_all(0:n_data(i)-1, i), psym = -2, col = i+1, thick = 4
    errplot, datax_all(0:n_data(i)-1, i), $
             ratio_storm_all(0:n_data(i)-1, i)-error_storm_all(0:n_data(i)-1, i), $
             ratio_storm_all(0:n_data(i)-1, i)+error_storm_all(0:n_data(i)-1, i), col = i+1, thick = 4
    xyouts, 1+i*2, 0.8, region_name(i), col = i+1, charsize = 3, charthick = 34
ENDFOR 
pclose
;ENDIF 

;  ratio difference
IF n_sta EQ 2 THEN BEGIN 
    IF min(datax_all(*, 0), /nan) GE min(datax_all(*, 1), /nan) THEN BEGIN
        ii = 0  &    ct = 0
        WHILE ct EQ 0  DO BEGIN 
        IF datax_all(ii, 0) GT -999 THEN index = where(datax_all(*, 1) EQ datax_all(ii, 0), ct)
        ii = ii+1
    ENDWHILE 
    min0 = ii-1 & min1 = index(0)

    ii = n_data(0)-1  &    ct = 0
    WHILE ct EQ 0  DO BEGIN 
        IF datax_all(ii, 0) GT -999 THEN index = where(datax_all(*, 1) EQ datax_all(ii, 0), ct)
        ii = ii-1
    ENDWHILE 
    max0 = ii-1  &  max1 = index(0)
ENDIF ELSE BEGIN 
    ii = 0 &  ct = 0
    WHILE ct EQ 0  DO BEGIN 
        IF datax_all(ii, 1) GT -999 THEN index = where(datax_all(*, 0) EQ datax_all(ii, 1), ct)
        ii = ii+1
    ENDWHILE 
    min1 = ii-1  &  min0 = index(0)
    ii = n_data(1)-1 &  ct = 0
    WHILE ct EQ 0  DO BEGIN 
        IF datax_all(ii, 1) GT -999 THEN index = where(datax_all(*, 0) EQ datax_all(ii, 1), ct)
        ii = ii-1
    ENDWHILE 
    max1 = ii-1  &  max0 = index(0)
ENDELSE  

new_datax = datax_all(min0:max0, 0)
ratio_storm_diff = ratio_storm_all(min1:max1, 1)/ratio_storm_all(min0:max0, 0)
ratio_nonstorm_diff = ratio_nonstorm_all(min1:max1, 1)/ratio_nonstorm_all(min0:max0, 0)
error_storm_diff = sqrt(error_storm_all(min1:max1, 1)^2/ratio_storm_all(min0:max0, 0)^2+ $
                        ratio_storm_all(min1:max1, 1)^2*error_storm_all(min0:max0, 0)^2/ratio_storm_all(min0:max0, 0)^4)
error_nonstorm_diff = sqrt(error_nonstorm_all(min1:max1, 1)^2/ratio_nonstorm_all(min0:max0, 0)^2+ $
                           ratio_nonstorm_all(min1:max1, 1)^2*error_nonstorm_all(min0:max0, 0)^2/ratio_nonstorm_all(min0:max0, 0)^4)

IF keyword_set(ps) THEN popen, path+para_name+'_ratio_div.ps', /land ELSE window, /free
plot, new_datax, ratio_storm_diff,  $
      yrange = [0.1, 2], xrange = xrange, xstyle = 1, ylog = 0, ystyle = 1, $ 
      title = 'Storm Beam Ratio Difference Distribution for ' + para_name, $
      xtitle = para_name, ytitle = 'Beam Ratio Difference', $
      charsize = 2, charthick = 3, /nodata

oplot, new_datax, ratio_storm_diff, psym = -2, col = 1, thick = 4
errplot, new_datax, ratio_storm_diff-error_storm_diff, $
         ratio_storm_diff+error_storm_diff, col = 1, $
         thick = 4
xyouts, 0, 1.75, 'storm', col = 1,  charsize = 3, charthick = 3

oplot, new_datax, ratio_nonstorm_diff, psym = -2, col = 2, thick = 4
errplot, new_datax, ratio_nonstorm_diff-error_nonstorm_diff, $
         ratio_nonstorm_diff+error_nonstorm_diff, col = 2, $
         thick = 4
xyouts, 0, 1.8, 'nonstorm', col = 2, charsize = 3, charthick = 3

oplot, [-999, 999], [1, 1]

pclose
ENDIF 
return, para_name
END 
