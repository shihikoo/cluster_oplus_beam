PRO PLOT_DIFF_YEAR
path = 'output/o_beam/test_year_dependence_with_HIA/'
tplot_restore, filename = path+'year_efficience.tplot'
aa = 1
ps = 0
IF  aa EQ  0 THEN  BEGIN  
    tplot_names, '*ratio_nonstorm', names = names
    n = n_elements(names)
    FOR  i = 0, n-1 DO  BEGIN  
        name = names(i)
        get_data, name, data = data
        store_data, name, data = {x: time_double(string(data.x)), y:data.y, dy:data.dy}
    ENDFOR  
ENDIF  

IF  aa EQ  1 THEN  BEGIN  
    FOR iphase = 0, 1 DO BEGIN 
        IF iphase EQ 0 THEN phase = 'storm' ELSE phase = 'nonstorm'
        IF keyword_set(ps) THEN popen, path+'ratio_year_'+phase+'.ps', /land 

        get_data, 'CODIF_north_lobe_data_distribution_ratio_'+phase, data = data
        year = data.x
        codif_north = data.y
        codif_north_error = data.dy
        get_data, 'CODIF_south_lobe_data_distribution_ratio_'+phase, data = data
        codif_south = data.y
        codif_south_error = data.dy
        get_data, 'HIA_north_lobe_data_distribution_ratio_'+phase, data = data
        hia_north = data.y
        hia_north_error = data.dy
        get_data, 'HIA_south_lobe_data_distribution_ratio_'+phase, data = data
        hia_south = data.y
        hia_south_error = data.dy

        plot, [2001, 2006], [0.001, 1], ylog = 1, $
              title = 'CODIF/HIA Beam Fraction for different years at north/south lobe during '+phase, /nodata

        oplot, year, codif_north, color = 1, psym = -1, linestyle = 2, thick = 4
        xyouts, 2002, 0.04, 'CODIF NORTH - -', color = 1
        errplot, year, codif_north-codif_north_error, $
                 codif_north+codif_north_error, color = 1, thick = 4

        oplot, year, codif_south, color = 1, psym = -1, thick = 4
        xyouts, 2002, 0.03, 'CODIF SOUTH __', color = 1
        errplot, year, codif_south-codif_south_error, $
                 codif_south+codif_south_error, color = 1, thick = 4

        oplot, year, hia_north, color = 4, psym = -1, linestyle = 2, thick = 4
        xyouts, 2002, 0.02, 'HIA NORTH - -', color = 4
        errplot, year, hia_north-hia_north_error, $
                 hia_north+hia_north_error, color = 4, thick = 4

        oplot, year, hia_south, color = 4, psym = -1, thick = 4
        xyouts, 2002, 0.015, 'HIA SOUTH __', color = 4
        errplot, year, hia_south-hia_south_error, $
                 hia_south+hia_south_error, color = 4, thick = 4
        IF keyword_set(ps) THEN pclose
    ENDFOR 
ENDIF   

spawn, 'mogrify -format png '+path+'*.ps '
spawn, 'mogrify -rotate -90 '+path+'*.png &'

stop
END
