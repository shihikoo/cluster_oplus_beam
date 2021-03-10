PRO plot_3d_distribution, sta_name, xtitle = xtitle, ytitle = ytitle, region = region

get_data, sta_name, data = data

IF NOT KEYWORD_SET(xtitle) THEN xtitle = ''
IF NOT KEYWORD_SET(ytitle) THEN ytitle = ''

    FOR i = 0, 5 DO BEGIN 
        popen, region+'/events_'+data.name(i)+'.ps', /port
        specplot, data.x, data.y,  data.sta(*, *, i), $
                  no_interp = 1, $
                  lim = { zlog:1, xlog:0, ylog:0, $
                          zrange:[1., 1000.], $
                          title: data.name(i), $
                          xtitle: xtitle, $
                          ytitle: ytitle, $
                          xrange: [0, 0.13], $
                          yrange: [0, 10], $
                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                          position: [0.13, 0.08, 0.88, 0.96]} 
        pclose
    ENDFOR 

    FOR i = 0, 2 DO BEGIN 
        popen, region+'/ratio_'+data.name(i*2+1)+'.ps', /port
        specplot, data.x, data.y,  data.sta(*, *, i*2+1)/data.sta(*, *, i*2), $
                  no_interp = 1, $
                  lim = { zlog:0, xlog:0, ylog:0, $
                          zrange:[0, 1.], $
                          title: data.name(i*2+1)+'  %', $
                          xtitle: xtitle, $
                          ytitle: ytitle, $
                          xrange: [0, 0.13], $ ;  1.3  0.13
                          yrange: [0, 10], $   ; 67,  10
                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                          position: [0.13, 0.08, 0.88, 0.96]} 
        pclose
    ENDFOR 


END 
