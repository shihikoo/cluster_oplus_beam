; use to calculated the counts distribution of input data and can also
; give the ratio distribution with flag information
; data_y input will give a distribution in 2d, which I don't use often
FUNCTION  data_distribution, x, flag=flag, storm_phase=storm_phase, data_y = y, xgrid = grid_x, ygrid = grid_y, $
                             para_name = para_name, xrange = xrange, yrange = yrange, xlog=xlog,ylog=ylog,region_name = region_name, $
                             write_data = write_data, ratio_distribution = ratio_distribution, events_distribution = events_distribution, $
                             path = path, ps = ps, plot_single_phase_events_distribution = plot_single_phase_events_distribution, $
                             barcolor_input = barcolor_input, bar_yrange_input = bar_yrange_input, single_phase = single_phase, $
                             charsize_input = charsize_input, ratio_plot_range=ratio_plot_range,rotate=rotate,quartile_range=quartile_range

IF NOT KEYWORD_SET(barcolor_input) THEN barcolor_input = 1
IF NOT KEYWORD_SET(charsize_input) THEN charsize_input = 1.5
IF NOT KEYWORD_SET(y) THEN  BEGIN 
    y = 0 & grid_y = 1
ENDIF 
IF keyword_set(y_range) THEN BEGIN 
    y_max = max(yrange) & y_min = min(yrange)
ENDIF ELSE BEGIN
    y_max = CEIL(MAX(y, /nan)) & y_min = FLOOR(MIN(y, /NAN) < 0)
    y_range = [0.1, 1]
ENDELSE 
IF keyword_set(xrange) THEN BEGIN 
    x_max = max(xrange) & x_min = min(xrange) 
ENDIF ELSE BEGIN
    x_max = CEIL(MAX(x, /NAN)) & x_min = FLOOR(MIN(x, /NAN) < 0)
ENDELSE
if not keyword_set(xlog) then xlog=0
if not keyword_set(ylog) then ylog=0
IF NOT keyword_set(region_name) THEN region_name = 'All_Region'

IF NOT KEYWORD_SET(grid_x) THEN grid_x = (x_max-x_min)/50.
IF NOT KEYWORD_SET(grid_y) THEN grid_y = (y_max-y_min)/50.

IF NOT KEYWORD_SET(flag) THEN flag=replicate(1,n_elements(x))
if not keyword_set(storm_phase) then storm_phase=replicate(0,n_elements(x))
if not keyword_set(para_name) then para_name=''
if not keyword_set(region_name) then region_name=''
if not keyword_set(path) then path=''

nx = CEIL((x_max - x_min)/grid_x)+1
ny = CEIL((y_max - y_min)/grid_y)+1

sta = FLTARR(nx, ny, 6)

FOR  i = 0l, nx-1 DO BEGIN 
    FOR j = 0l, ny-1 DO BEGIN
        index = where(x GE (i*grid_x+x_min) AND x LT ((i+1)*grid_x+x_min) AND $
                      y GE (j*grid_y+y_min) AND y LT ((j+1)*grid_y+y_min) AND $
                      abs(flag) GE 0, ct) 
        sta(i, j, 0) = ct
        index = where(x GE (i*grid_x+x_min) AND x LT ((i+1)*grid_x+x_min) AND $
                      y GE (j*grid_y+y_min) AND y LT ((j+1)*grid_y+y_min) AND $
                      abs(flag) GT 0, ct)
        sta(i, j, 1) = ct
        index = where(x GE (i*grid_x+x_min) AND x LT ((i+1)*grid_x+x_min) AND $
                      y GE (j*grid_y+y_min) AND y LT ((j+1)*grid_y+y_min) AND $
                      abs(flag) GE 0 AND  (storm_phase EQ 0 or storm_phase eq 5), ct) 
        sta(i, j, 2) = ct
        index = where(x GE (i*grid_x+x_min) AND x LT ((i+1)*grid_x+x_min) AND $
                      y GE (j*grid_y+y_min) AND y LT ((j+1)*grid_y+y_min) AND $
                      abs(flag) GE 1 AND  (storm_phase EQ 0 or storm_phase eq 5), ct) 
        sta(i, j, 3) = ct
        index = where(x GE (i*grid_x+x_min) AND x LT ((i+1)*grid_x+x_min) AND $
                      y GE (j*grid_y+y_min) AND y LT ((j+1)*grid_y+y_min) AND $
                      abs(flag) GE 0 AND  (storm_phase GE 1 AND storm_phase LE 3), ct) 
        sta(i, j, 4) = ct
        index = where(x GE (i*grid_x+x_min) AND x LT ((i+1)*grid_x+x_min) AND $
                      y GE (j*grid_y+y_min) AND y LT ((j+1)*grid_y+y_min) AND $
                      abs(flag) GE 1 AND  (storm_phase GE 1 AND storm_phase LE 3), ct) 
        sta(i, j, 5) = ct
    ENDFOR
ENDFOR

;sta(*, *, 4) = sta(*, *, 0)-sta(*, *, 2)
;sta(*, *, 5) = sta(*, *, 1)-sta(*, *, 3)
datax = indgen(nx)*grid_x+x_min+ grid_x*0.5
datay = indgen(ny)*grid_y+y_min
name = ['all_data',  'beam_for_all_data', 'nonstorm', 'beam_for_nonstorm_data', $
        'storm', 'beam_for_storm_data']

sta_str =  para_name+'_'+region_name+'_data_distribution'
;stop
IF N_ELEMENTS(datay) GT 1 THEN $
  store_data, sta_str, data = {x:datax, y:datay, sta:sta, name:name} $
ELSE store_data, sta_str, data = {x:datax, sta:sta, name:name}

IF NOT keyword_set(y) THEN BEGIN 
;---------data
    IF keyword_set(write_data) THEN BEGIN 
        output = FLTARR(nx, 7)
        output(*, 0) = datax
        output(*, 1:6) = reform(sta)
        str = {x:datax, y:output, v:para_name+'  '+region_name}
        store_data, 'dump_data_s', data = str
        fln_dump = 'sta/'+para_name+'_sta.dat'
        dump_data, 'dump_data_s', file_out = fln_dump
    ENDIF 
;---------Ratio
    IF keyword_set(ratio_distribution) THEN BEGIN 
        error = sqrt(sta)
        error_nonstorm = sqrt(error(*, 0, 3)^2/sta(*, 0, 2)^2+ $
                              sta(*, 0, 3)^2*error(*, *, 2)^2/sta(*, *, 2)^4)
        error_storm = sqrt(error(*, 0, 5)^2/sta(*, 0, 4)^2+ $
                           sta(*, 0, 5)^2*error(*, *, 4)^2/sta(*, *, 4)^4)
        index_nonstorm = where(sta(*, 0, 2) EQ 0 OR sta(*, 0, 3) EQ 0, ct)

        ratio_nonstorm = sta(*, 0, 3)/sta(*, 0, 2)
        IF ct GT 0 THEN BEGIN 
            ratio_nonstorm(index_nonstorm) = !VALUES.F_NAN
            error_nonstorm(index_nonstorm) = !VALUES.F_NAN
        ENDIF 
        index_storm = where(sta(*, 0, 4) EQ 0, ct)
        ratio_storm = sta(*, 0, 5)/sta(*, 0, 4)
        IF ct GT 0 THEN BEGIN 
            ratio_storm(index_storm) = !VALUES.F_NAN
            error_storm(index_storm) = !VALUES.F_NAN
        ENDIF 
        index = where(sta(*, 0, 2) EQ 0 AND sta(*, 0, 4)EQ 0, ct)
        IF ct GT 0 THEN datax(index) = !VALUES.F_NAN
;---write it into ps file
        IF NOT keyword_set(xrange) THEN xrnage = [min(datax)-1, max(datax)+1]
        if not keyword_set (ratio_plot_range) then begin 
            if keyword_set(ylog) then ratio_plot_range = [0.001,1] else ratio_plot_range=[0,1]  
        endif 
        IF keyword_set(ps) THEN $
          popen, path+para_name+'_'+region_name+'_ratio.ps', /land ;ELSE window, /free
        plot, datax, ratio_nonstorm, yrange = ratio_plot_range, xrange = xrange, ylog = ylog $
          , xstyle = 1, charsize = 2, charthick = 3, xlog=xlog,$
          title = 'Beam Ratio Distribution for '+para_name+'!C'+region_name, $
          xtitle = para_name, ytitle = 'Beam Ratio',$
          position = [0.15, 0.15, 0.95, 0.85], /nodata
        oplot, datax, ratio_nonstorm, psym = -7, col = 2, thick = 8
        errplot, datax, ratio_nonstorm-error_nonstorm,$
          ratio_nonstorm+error_nonstorm, col = 2, thick = 4
        oplot, datax, ratio_storm, psym = -7, col = 1, thick = 8
        errplot, datax, ratio_storm-error_storm, ratio_storm+error_storm, col = 1, thick = 4
        xyouts, (xrange(0)+xrange(1))/2., 0.005, 'non storm', col = 2, charsize = 3, charthick = 3
        xyouts, (xrange(0)+xrange(1))/2., 0.01, 'storm', col = 1, charsize = 3, charthick = 3
        IF keyword_set(ps) THEN pclose 
        store_data, sta_str+'_ratio_storm', data = {x:datax, y:ratio_storm, dy:error_storm}
        store_data, sta_str+'_ratio_nonstorm', $
          data = {x:datax, y:ratio_nonstorm, dy:error_nonstorm}   
    ENDIF  
;----------Events
    IF keyword_set(events_distribution) THEN BEGIN 
        loc = sort(abs(datax))  
        IF x_min LT 0 THEN  index_bar = loc(0)+INDGEN(51)-24 ELSE index_bar = loc(0)+INDGEN(51)
        
        barcolor = sta(index_bar, 0, 2)
        index = where(datax(index_bar) GT 0, ct)
        IF ct GT 0 THEN  barcolor(index) = 1
        index = where(datax(index_bar) LT  0, ct)
        IF ct GT 0 THEN barcolor(index) = 2
                                ; stop
        IF keyword_set(ps) THEN popen, path+para_name+'_'+region_name+'_events.ps', /land ; ELSE window, /free
        bar_plot, sta(index_bar, 0, 2), baserange = 0.25, $
          xtitle = '  all events, nonstorm        beam event, nonstorm       all events, storm        beam events, storm ', $
          ytitle = 'Events', colors = barcolor, title = para_name+'   '+region_name+'  Data Distribution'
        bar_plot, sta(index_bar, 0, 3), baserange = 0.25, $
          baroffset = 60, colors = barcolor+1, /overplot
        bar_plot, sta(index_bar, 0, 4), baserange = 0.25, $
          baroffset = 120, colors = barcolor+2, /overplot
        bar_plot, sta(index_bar, 0, 5), baserange = 0.25, $
          baroffset = 180, colors = barcolor+3, /overplot
        
        max_num = fltarr(4)
        FOR im = 0, 3 DO BEGIN 
            max_num(im) =  max(sta(index_bar, 0, im+2))
            IF keyword_set(neg)THEN offsets = 6 ELSE offsets = 1
            xyouts, offsets+11*im, max_num(im), $
              string(datax(where(sta(*, 0, im+2)EQ max_num(im))), format = '(f5.1)')
        ENDFOR 
        IF keyword_set(ps) THEN pclose  
    ENDIF  
;plot events bar plot into single page
    IF keyword_set(plot_single_phase_events_distribution) THEN BEGIN 
        loc = sort(datax)   
; loc = sort(datax) may be right for property plot
;    IF x_min LT 0 THEN index_bar = loc(0)+INDGEN(nx)-(nx/2)+1 else$
        index_bar = loc(0)+INDGEN(nx)
        
        barcolor = sta(index_bar, 0, 2)
        index = where(datax(index_bar) GT 0, ct)
        IF ct GT 0 THEN  barcolor(index) = 1
        index = where(datax(index_bar) LT  0, ct)
        IF ct GT 0 THEN barcolor(index) = 2
        
        IF NOT  keyword_set(single_phase) THEN single_phase = 'alltime'
        IF single_phase EQ 'alltime_alldata' THEN BEGIN 
            is = 0 & ind_sp = where(storm_phase GE 0) 
        ENDIF 
        IF single_phase EQ 'alltime' THEN BEGIN 
            is = 1 & ind_sp = where(storm_phase GE 0 and ABS(flag) gt 0)
        ENDIF 
        IF single_phase EQ 'nonstorm_alldata' THEN BEGIN 
            is = 2 & ind_sp = where(storm_phase EQ 0 or storm_phase eq 5)
        ENDIF 
        IF single_phase EQ 'nonstorm' THEN BEGIN 
            is = 3 & ind_sp = where((storm_phase EQ 0 or storm_phase eq 5) and ABS(flag) gt 0)
        ENDIF 
        IF single_phase EQ 'storm_alldata' THEN BEGIN 
            is = 4 & ind_sp = where(storm_phase Ge 1 and storm_phase le 4)
        ENDIF
        IF single_phase EQ 'storm' THEN BEGIN 
            is = 5 & ind_sp = where(storm_phase Ge 1 and storm_phase le 4 and ABS(flag) gt 0)
        ENDIF 

        barnames = strarr(n_elements(index_bar))
        max_num =  max(sta(index_bar, 0, is), /nan)
        max_num_loc = datax(where(sta(*, 0, is)EQ max_num))
        
        index=where(finite(x(ind_sp)),ct)
        x_new=x(ind_sp)
        if ct gt 0 then begin 
            median_value = median(x_new)
            mean_value=mean(x_new,/nan)
            index=where(x_new le median_value,ct)
            if ct gt 0 then low_quartile=median(x_new(index)) else low_quartile = !values.f_nan
            index=where(x_new ge median_value,ct)
            if ct gt 0 then high_quartile=median(x_new(index)) else high_quartile=!values.f_nan
            error = sqrt(total((x_new-mean_value)^2,/nan)/ct)
        endif else begin 
            median_value=!values.f_nan
            mean_value=!values.f_nan
            error=!values.f_nan
            low_quartile=!values.f_nan
            high_quartile=!values.f_nan
        endelse   

        store_data, delete= 'mean_median_error'
        store_data, 'mean_median_error', data = {mean:mean_value, median:median_value, error:error}
        
        IF NOT keyword_set(bar_yrange_input) THEN bar_yrange = [0,max_num] ELSE bar_yrange = bar_yrange_input
        if keyword_set(quartile_range) then xrange=[low_quartile,high_quartile]

        if keyword_set(rotate) then begin
            plot, [x_min, x_min+(nx)*grid_x], [0., 0], xrange = bar_yrange, yrange=xrange,$
              charsize = charsize_input, ticklen = 0, xstyle = 1, $
              ytitle=para_name,xtitle = 'Events', title = region_name+' '+single_phase +' Data Distribution' 
                                ;      position=[0.15,0.15,0.95,0.85]
            bar_plot, sta(index_bar, 0, is), baserange = 1, colors = barcolor+barcolor_input, over = 1,rotate=rotate
            oplot,  [ 0, max_num], [median_value, median_value], color = 1, thick = 10 
            oplot,  [ 0, max_num/2.], [low_quartile, low_quartile], color = 1, thick = 10 
            oplot,  [ 0, max_num/2.], [high_quartile, high_quartile], color = 1, thick = 10 

            xyouts, max(bar_yrange)*0.5, mean(xrange)*1.1, $
              'meidan:'+ strcompress(string(median_value,format='(f10.3)'),/remove_all) $
              +'!Cmean:'+ strcompress(string(mean_value,format='(f10.3)'),/remove_all) $
              +'!Cerror:'+ strcompress(string(error,format='(f10.3)'),/remove_all) $
              , charsize = 2    ;3-charsize_input
        endif else begin 
            plot, [x_min, x_min+(nx)*grid_x], [0., 0], xrange = xrange, yrange=bar_yrange, $
              charsize = charsize_input, ticklen = 0, xstyle = 1, $
              ytitle = 'Events', title = para_name+'  '+region_name+' '+single_phase +' Data Distribution'
            bar_plot, sta(index_bar, 0, is), baserange = 1, colors = barcolor+barcolor_input, over = 1,rotate=rotate
            oplot,  [median_value, median_value], [ 0, max_num], color = 1, thick = 10 
            oplot,  [low_quartile,low_quartile], [ 0, max_num/2.], color = 1, thick = 10 
            oplot,  [high_quartile,high_quartile], [ 0, max_num/2.], color = 1, thick = 10 
            xyouts, mean(xrange)*1.1, max(bar_yrange)*0.9, $
              'meidan:'+ strcompress(string(median_value,format='(f10.3)'),/remove_all) $
              +'!Cmean:'+ strcompress(string(mean_value,format='(f10.3)'),/remove_all) $
              +'!Cerror:'+ strcompress(string(error,format='(f10.3)'),/remove_all) $
                                ;    + '!Cbar grid:'+strcompress(grid_x, /remove_all) $
            , charsize = 1      ;3-charsize_input
        endelse 
    ENDIF 
ENDIF  
;stop
RETURN, sta_str
END
