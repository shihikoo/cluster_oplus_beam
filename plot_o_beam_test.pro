;Purpose: Find O+ beam over time interval (input) with displaytime 
;         for every displaytime (input)
;
;Input: sc           : Cluster no. if not set the default is 4
;       calc_time    : the calculate time for each loop
;       time_start   : in idl time format
;       time_end     : in idl time format 
;       displaytime  : the displaytime for single plot
;
;Keywords:
;       average_time : time for averaging data, if not set the default is 5 min
;       idl_plot     : plot the result plot in idl_window
;       ps           : plot the result plot in ps files,
;       dumpdata     : output data file
;       globe_plot   : plot a set of globe plot to show the selected
;                      range for plotting mom
;       store_data   : store_data into .tplot      default: 1 
;       plot_mom     : plot mom of the identifed O+ beam
;
;Output: Depends on Keywords settings 
;        There will also be two .log files
;
;Written by Jing Liao  02/20/2008
;
PRO plot_o_beam_test, time_start = time_start, time_end = time_end
;20020611,12,13
sc = 4                          ;set the satallite number 
sc_str = STRING(sc, FORMAT = '(i1.1)')
s = '2'
e = '4'
IF NOT keyword_set(time_start) THEN  time_start = '2002-10-01/0'+s+':00:00'
IF NOT keyword_set(time_end) THEN time_end = '2002-10-01/0'+e+':00:00'
calc_time = 2.*60. *60.         ;in seconds
;-----------------------------------------
;Set the keywords used in find_o_beam.pro

;------------------------------------------
beam_recalc = 0
mom_recalc = 1
find_phase = 0
add_imf = 0

plot_imf = 0
plot_mom = 0
idl_plot = 0
ps = 0
dumpdata = 0
store_data = 1

globe_plot = 0
plot_lowcount_filter = 0

show_fit = 0
fit_result = 0
dfit_temperature = 1

dont_plot_sw_sheath = 1
display_time = 2.*60*60
path_set = ['original_setting', 'limit_no', 'limit_en', 'ar_2250', $
            'ar_3375', 'avg10min', 'avg15min']
FOR ipath = 0, 6 DO  BEGIN
    tplot_names, names = names
    store_data, delete = names

    path = 'test_mom/'+ path_set(ipath)+'/' 

    spawn,  'mkdir '+path
    spawn, 'mkdir '+path+'tplot_restore/'
    spawn, 'mkdir '+path+'data/'
    spawn, 'mkdir '+path+'plots/'

    use_angle_range = 1
    use_energy_range = 1
    average_time = 5 * 60       ;in seconds  
    beam_angle_range = 11.25

    IF ipath EQ 1 THEN BEGIN 
        use_angle_range = 0
        use_energy_range = 0
    ENDIF 

    IF ipath EQ 2 THEN use_angle_range = 0
    IF ipath EQ 3 THEN beam_angle_range = 22.5
    IF ipath EQ 4 THEN beam_angle_range = 33.75
    IF ipath EQ 5 THEN average_time = 10* 60   
    IF ipath EQ 6 THEN average_time = 15* 60   

;-----------------------------------------------
;Write [START] in log files
;-----------------------------------------------
    OPENU, unit, path+'log_plotted.txt', /GET_LUN, /APPEND
    PRINTF, unit, SYSTIME(), '[START]'
    FREE_LUN, unit         

    OPENU, unit, path+'log_errors.txt', /GET_LUN, /APPEND
    PRINTF, unit, SYSTIME(), '[START]'
    FREE_LUN, unit     
;---------------------------------------------------
; Set the loop as requested
;---------------------------------------------------
    ts = time_double(time_start)
    te = time_double(time_end)
    ntime = CEIL((te - ts)/calc_time) 
;------------------------------------------------------------ 
    FOR i = 0l, ntime-1 DO BEGIN  
; Timespan over each displaytime
        timespan, ts + i*calc_time, calc_time, /seconds
;some days cannot be plotted because of the orbit or else
        time_useless = [['2002-06-11/00:00:00', '2002-06-14/00:00:00'], $
                        ['2005-10-10/00:00:00', '2005-10-11/00:00:00'], $
                        ['2006-03-23/00:00:00', '2006-03-24/00:00:00']]           
        nouse = 0
        FOR j = 0, N_ELEMENTS(time_useless(0, *))-1 DO BEGIN 
            nouse = nouse+(ts+i*calc_time  GE  time_double(time_useless(0, j)) AND $
                           ts+(i+1)*calc_time LE  time_double(time_useless(1, j))) 
        ENDFOR 
        
        IF nouse EQ 0 THEN BEGIN 
; identify O+ beam plot
            find_o_beam, sc = sc, $
                         average_time = average_time, $ 
                         idl_plot = idl_plot, $
                         ps = ps, $
                         dumpdata = dumpdata, $
                         globe_plot = globe_plot, $
                         store_data = store_data, $
                         plot_mom = plot_mom, $
                         path = path, $
                         plot_lowcount_filter =  plot_lowcount_filter, $         
                         beam_recalc = beam_recalc, $
                         mom_recalc = mom_recalc,  $
                         find_phase = find_phase, $
                         displaytime = display_time, $
                         add_imf = add_imf, $
                         dont_plot_sw_sheath = dont_plot_sw_sheath, $
                         use_angle_range = use_angle_range, $
                         use_energy_range = use_energy_range, $
                         plot_imf = plot_imf, $
                         beam_angle_range =  beam_angle_range, $
                         dfit_temperature = dfit_temperature, $
                         show_fit = show_fit
        ENDIF 

        at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL) 
        t_name_tail = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_PHI90_270_MTTEMPERATURE_SP3_ET0_All_T_AVG'+at_str
        t_name_tail_para = 'TDMOM_ENVARIOUS_SC'+sc_str $
                           +'_PHI90_270_MTTEMPERATURE_SP3_ET0_All_X_AVG'+at_str
        t_name_tail_perp1 = 'TDMOM_ENVARIOUS_SC' +sc_str $
                            +'_PHI90_270_MTTEMPERATURE_SP3_ET0_All_Y_AVG'+at_str
        t_name_tail_perp2 = 'TDMOM_ENVARIOUS_SC'+sc_str $
                            +'_PHI90_270_MTTEMPERATURE_SP3_ET0_All_Z_AVG'+at_str
        t_fit_name_tail_para = 'TDMOM_ENVARIOUS_SC' + sc_str $
          +'_PHI90_270_MTTEMPERATURE_SP3_ET0_All_para_dfit_AVG'+at_str
        t_fit_name_tail_perp = 'TDMOM_ENVARIOUS_SC'+sc_str $
          +'_PHI90_270_MTTEMPERATURE_SP3_ET0_All_perp_dfit_AVG'+at_str
        t_fit_name_tail = 'TDMOM_ENVARIOUS_SC'+sc_str $
                          +'_PHI90_270_MTTEMPERATURE_SP3_ET0_All_T_dfit_AVG'+at_str
        
        t_name_earth = 'TDMOM_ENVARIOUS_SC'+sc_str $
                       +'_PHI270_90_MTTEMPERATURE_SP3_ET0_All_T_AVG'+at_str
        t_name_earth_para = 'TDMOM_ENVARIOUS_SC'+sc_str $
                            +'_PHI270_90_MTTEMPERATURE_SP3_ET0_All_X_AVG'+at_str
        t_name_earth_perp1 = 'TDMOM_ENVARIOUS_SC'+sc_str $
                             +'_PHI270_90_MTTEMPERATURE_SP3_ET0_All_Y_AVG'+at_str
        t_name_earth_perp2 = 'TDMOM_ENVARIOUS_SC'+sc_str $
                             +'_PHI270_90_MTTEMPERATURE_SP3_ET0_All_Z_AVG'+at_str
        t_fit_name_earth_para = 'TDMOM_ENVARIOUS_SC'+sc_str $
          +'_PHI270_90_MTTEMPERATURE_SP3_ET0_All_para_dfit_AVG'+at_str
        t_fit_name_earth_perp = 'TDMOM_ENVARIOUS_SC'+sc_str $
          +'_PHI270_90_MTTEMPERATURE_SP3_ET0_All_perp_dfit_AVG'+at_str
        t_fit_name_earth = 'TDMOM_ENVARIOUS_SC'+sc_str $
                           +'_PHI270_90_MTTEMPERATURE_SP3_ET0_All_T_dfit_AVG'$
                           +at_str

        d_name_tail = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_PHI90_270_MTDENSITY_SP3_ET0_All_AVG'+at_str
        
        d_fit_name_tail = 'TDMOM_ENVARIOUS_SC'+sc_str $
                          +'_PHI90_270_MTDENSITY_SP3_ET0_All_dfit_AVG'+at_str
        d_name_earth = 'TDMOM_ENVARIOUS_SC'+sc_str $
                       +'_PHI270_90_MTDENSITY_SP3_ET0_All_AVG'+at_str
        d_fit_name_earth = 'TDMOM_ENVARIOUS_SC'+sc_str $
                           +'_PHI270_90_MTDENSITY_SP3_ET0_All_dfit_AVG'+at_str

        ylim, [t_fit_name_tail, t_fit_name_tail_para, $
               t_fit_name_tail_perp, $
               t_fit_name_earth, t_fit_name_earth_para, $
               t_fit_name_earth_perp], 0.1, 10000
        options, [t_name_tail, t_name_tail_para,$
                  t_name_tail_perp1, t_name_tail_perp2, $
                  t_name_earth, t_name_earth_para, $
                  t_name_earth_perp1, t_name_earth_perp2], $
                 'color', 6
;tail        
        IF NOT keyword_set(fit_result) THEN $
          popen, path+'20021001_'+s+'_'+e+'temperature_compare_tail.ps', /land ELSE window, /free
        
        tplot,  [t_fit_name_tail, t_fit_name_tail_para, $
                 t_fit_name_tail_perp], title = path
        
        tplot_panel, v = t_fit_name_tail, o = t_name_tail, psym = 7
        tplot_panel, v = t_fit_name_tail_para, o = t_name_tail_para, psym = 7
        tplot_panel, v = t_fit_name_tail_perp, o = t_name_tail_perp1, psym = 7
        tplot_panel, v = t_fit_name_tail_perp, o = t_name_tail_perp2, psym = 7

        xyouts, 7000, 50, 'T (fit)', size = 2
        xyouts, 7000, 1000, 'T (3dmom)', color = 6, size = 2
        IF NOT keyword_set(show_fit) THEN   pclose
;earth
        IF NOT keyword_set(fit_result) THEN $
          popen, path+'20021001_'+s+'_'+e+'temperature_compare_earth.ps', /land ELSE window, /free
        
        tplot,  [t_fit_name_earth, t_fit_name_earth_para, $
                 t_fit_name_earth_perp], title = path
        
        tplot_panel, v = t_fit_name_earth, o = t_name_earth, psym = 7
        tplot_panel, v = t_fit_name_earth_para, o = t_name_earth_para, psym = 7
        tplot_panel, v = t_fit_name_earth_perp, o = t_name_earth_perp1, psym = 7
        tplot_panel, v = t_fit_name_earth_perp, o = t_name_earth_perp2, psym = 7

        xyouts, 7000, 50, 'T (fit)', size = 2
        xyouts, 7000, 1000, 'T (3dmom)', color = 6, size = 2
        IF NOT keyword_set(show_fit) THEN   pclose
;density
        IF NOT keyword_set(fit_result) THEN $
          popen, path+'20021001_'+s+'_'+e+'density_compare.ps', /land ELSE window, /free
        
        tplot,  [d_fit_name_tail, d_fit_name_earth], title = path
        
        tplot_panel, v = d_fit_name_tail, o = d_name_tail, psym = 7
        tplot_panel, v = d_fit_name_earth, o = d_name_earth, psym = 7

        xyouts, 7000, 50, 'D (fit)', size = 2
        xyouts, 7000, 1000, 'D (3dmom)', color = 6, size = 2
        IF NOT keyword_set(show_fit) THEN   pclose
        
        spawn, 'mogrify -format png '+path+'*.ps'
        spawn, 'mogrify -rotate -90 '+path+'*.png'
    ENDFOR 

ENDFOR   
;------------------------------------------------------------
;----write [END] in log files
OPENU, unit, path+'log_plotted.txt', /GET_LUN, /APPEND
PRINTF, unit, SYSTIME(), '[END]'
FREE_LUN, unit         

OPENU, unit, path+'log_errors.txt', /GET_LUN, /APPEND
PRINTF, unit, SYSTIME(), '[END]'
FREE_LUN, unit  

stop
END 
