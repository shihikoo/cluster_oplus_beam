PRO map_o_beam, time, data, $
                sort_flag = sort_flag, sort_title = sort_title, $
                sc = sc, $
                ps_plot = ps_plot, idl_plot = idl_plot, $
                plot_path = plot_path, $
                coor_set = coor_set, $
                grid_set = grid_set, slice_grid_set = slice_grid_set, $
                diff_beta = diff_beta, $
                direction_set = direction_set, storm_phase_set = storm_phase_set, $
                plot_2d = plot_2d, slice_plot = slice_plot, waterdrop_plot = waterdrop_plot, $
                point_plot = point_plot, $
                events_map = events_map,  $
                property_map_set = property_map_set, $
                symmetry_template = symmetry_template, $
                mogrify_ps = mogrify_ps
;-- default keywords settings if keywords are not set-
IF NOT KEYWORD_SET(sc) THEN SC = 4
IF NOT KEYWORD_SET(coor_set) THEN $
  coor_set = [['X_GSE', 'Y_GSE'], ['X_GSE', 'Z_GSE'], ['X_GSM', 'Y_GSM'], ['X_GSM', 'Z_GSM']]
IF NOT KEYWORD_SET( grid_set)THEN   grid_set = [1., 2.]
IF NOT KEYWORD_SET(slice_grid_set)THEN  slice_grid_set = [20., 8.]  
IF NOT KEYWORD_SET( direction_set)THEN   direction_set =  ['all_directions', 'tailward', 'earthward', 'both']
IF NOT KEYWORD_SET( storm_phase_set)THEN  $
  storm_phase_set = ['all_phases', 'nonstorm_time', 'storm_time', 'prestorm', 'main_phase', 'recovery'] 
IF NOT KEYWORD_SET( PROPERTY_MAP_SET) THEN $
  PROPERTY_MAP_SET = ['']       ;['energy', 'flux','pitch_angle','density','temperature','velocity']
IF NOT KEYWORD_SET(plot_path) THEN BEGIN 
    plot_path = 'unknow_sort/'
    spawn, 'mkdir '+ plot_path
ENDIF 
IF NOT KEYWORD_SET(sort_title) THEN sort_title = ''

; -- other basic settings --
sc_str = STRING(sc, format = '(i1.1)')

X_GSE_RANGE = [20., -20.] & Y_GSE_RANGE = [20., -20.] & Z_GSE_RANGE = [-20., 20.]
X_GSM_RANGE = [20., -20.] & Y_GSM_RANGE = [20., -20.] & Z_GSM_RANGE = [-20., 20.]

EVENTS_V_LOG = 0 & EVENTS_V_RANGE = [0, 100.]
RATIO_V_LOG = 0  & RATIO_V_RANGE = [0, 1.]
ENERGY_V_LOG = 1 & ENERGY_V_RANGE = [40, 1000.]
FLUX_V_LOG = 0 & FLUX_V_RANGE = [0, 1000.]
Density_v_log = 1 &  density_v_range = [0.001, 0.1]
velocity_v_log = 0 &  velocity_v_range = [0, 100]
pitch_angle_v_log = 0 & pitch_angle_v_range = [0, 180]
nvpara_over_b_v_log = 1 &  nvpara_over_b_v_range = [0.001, 0.01]
temperature_v_log = 1 &  temperature_v_range = [1, 100]
anodes_v_log = 0 &  anodes_v_range = [0, 8]

la_x = -10 & la_y = -20         ; legend position
plot_path =  plot_path       &   spawn, 'mkdir ' + plot_path 

;-- time settings-- 
get_timespan, interval
ts_str = time_struct(interval(0)) ; start_time tplot time structure
te_str = time_struct(interval(1)) ; end_time tplot time structure
ts_date = string(ts_str.year, format = '(i4.4)')+ $
          string(ts_str.month, format = '(i2.2)')+ $
          string(ts_str.date, format = '(i2.2)')
te_date = string(te_str.year, format = '(i4.4)')+ $
          string(te_str.month, format = '(i2.2)')+ $
          string(te_str.date, format = '(i2.2)')

;-- save input data into different arraies --
flag_o = data(*, 0)
beta = ABS(data(*, 1))
x_gse = data(*, 2)
y_gse = data(*, 3)
storm_phase = data(*, 42)    
h_density = data(*, 29)
h_velocity = sqrt((data(*, 30))^2+(data(*, 31))^2 + (data(*, 32))^2)
energy_tail = data(*, 5)
energy_earth = data(*, 14)
mass_o = 16*1.6e-27*(1e3)^2/(1.6e-19) ; unit: ev/(km/s) 
energy_v_tail = sqrt(2*energy_tail/mass_o)
energy_v_earth = sqrt(2*energy_earth/mass_o)
flux_tail = data(*, 7)
flux_earth = data(*, 16)
density_tail = data(*, 8)
density_earth = data(*, 17)
velocity_tail = data(*, 9)
velocity_earth = data(*, 18)
pitch_angle_tail = data(*, 6)
pitch_angle_earth = data(*, 15)
v_para_tail =  data(*, 10)
v_para_earth = data(*, 19) 
B =  sqrt(data(*, 26)^2+data(*, 27)^2+data(*, 28)^2)
temperature_tail = data(*, 12)
temperature_earth = data(*, 21)
anodes_tail = (data(*, 50)+11.25)/22.5
anodes_earth = (data(*, 51)+11.25)/22.5

ntime = N_ELEMENTS(time)
; reset with sort_flag
flag_o = flag_o*sort_flag

; run for different storms
FOR i_storm_phase = 0, N_ELEMENTS(storm_phase_set)-1 DO BEGIN 
;reset for different storm phases
    phase = storm_phase_set(i_storm_phase)
    flag_phase = flag_o
; combine flag and storm_phase 
    IF phase EQ 'prestorm' THEN phase_flag = 1
    IF phase EQ 'main_phase' THEN phase_flag = 2
    IF phase EQ 'recovery' THEN phase_flag = 3
    IF phase EQ 'storm_time' THEN phase_flag = -1
    IF phase EQ 'nonstorm_time' THEN phase_flag = 0
    IF phase EQ 'all_phases' THEN phase_flag = -2
; the flags of data out of certain storm phase  are set to be 0  
    IF phase_flag GE 0 THEN flag_phase = flag_phase/(storm_phase EQ (phase_flag)) $
    ELSE BEGIN 
        IF phase_flag EQ -1 THEN flag_phase = flag_phase/(storm_phase GT 0)
    ENDELSE 

; set all infinite flag value to nan 
    index = WHERE( ~FINITE(flag_phase), ct)
    IF ct GT 0 THEN flag_phase(index) = !VALUES.F_NAN

;------------------------------------------------------------------------
;find the third axis and read data position info from "data" into new
;arrar data_pos and set the range, x,y are plotting axes and z is
;slice axis.  
;All useful info will be put into new array 
;------------------------------------------------------------------------
;if keyword plot_all set , plot all coordinates group set before in ALL
    nplot = N_ELEMENTS(coor_set)/2
    FOR iplot = 0, nplot-1 DO BEGIN 
        PLOT_AXIS = coor_set(*, iplot)
        IF  ARRAY_EQUAL(PLOT_AXIS, ['X_GSE', 'Y_GSE']) OR $
          ARRAY_EQUAL(PLOT_AXIS, ['Y_GSE', 'X_GSE']) THEN  plot_axis_3 = 'Z_GSE'
        IF  ARRAY_EQUAL(PLOT_AXIS, ['X_GSE', 'Z_GSE']) OR $
          ARRAY_EQUAL(PLOT_AXIS, ['Z_GSE', 'X_GSE']) THEN  plot_axis_3 = 'Y_GSE'
        IF  ARRAY_EQUAL(PLOT_AXIS, ['Y_GSE', 'Z_GSE']) OR $ 
          ARRAY_EQUAL(PLOT_AXIS, ['Z_GSE', 'Y_GSE']) THEN  plot_axis_3 = 'X_GSE'
        IF  ARRAY_EQUAL(PLOT_AXIS, ['X_GSM', 'Y_GSM']) OR $
          ARRAY_EQUAL(PLOT_AXIS, ['Y_GSM', 'X_GSM']) THEN  plot_axis_3 = 'Z_GSM'
        IF  ARRAY_EQUAL(PLOT_AXIS, ['X_GSM', 'Z_GSM']) OR $
          ARRAY_EQUAL(PLOT_AXIS, ['Z_GSM', 'X_GSM']) THEN  plot_axis_3 = 'Y_GSM'
        IF  ARRAY_EQUAL(PLOT_AXIS, ['Y_GSM', 'Z_GSM']) OR $ 
          ARRAY_EQUAL(PLOT_AXIS, ['Z_GSM', 'Y_GSM']) THEN  plot_axis_3 = 'X_GSM'     
        IF NOT KEYWORD_SET(plot_axis_3) THEN BEGIN 
            PRINT, 'AXES CAN ONLY BE X,Y,Z IN GSE OR GSM COORDINATE SYSTEM'+ $
                   '!C!CCONTINUE WITH DEFAUL AXES!C!CX_GSE VS Y_GSE !C!C INPUT ".CONT"'
            STOP
            PLOT_AXIS = ['X_GSE', 'Y_GSE']
            plot_axis_3 = 'Z_GSE'      
        ENDIF 
        plot_axis = [PLOT_AXIS, plot_axis_3]   
        
        data_pos = DBLARR(ntime, 3) 
        range = FLTARR(2, 3)
        FOR i = 0, 2 DO BEGIN 
            IF PLOT_AXIS(i) EQ 'X_GSE' THEN BEGIN 
                data_pos(*, i) = data(*, 2) 
                range(*, i) = X_GSE_RANGE
            ENDIF 
            IF PLOT_AXIS(i) EQ 'Y_GSE' THEN BEGIN 
                data_pos(*, i) = data(*, 3) 
                range(*, i) = Y_GSE_RANGE
            ENDIF 
            IF PLOT_AXIS(i) EQ 'Z_GSE' THEN BEGIN 
                data_pos(*, i) = data(*, 4) 
                range(*, i) = Z_GSE_RANGE
            ENDIF 
            IF PLOT_AXIS(i) EQ 'X_GSM' THEN BEGIN 
                data_pos(*, i) = data(*, 23) 
                range(*, i) = X_GSM_RANGE
            ENDIF 
            IF PLOT_AXIS(i) EQ 'Y_GSM' THEN BEGIN 
                data_pos(*, i) = data(*, 24)
                range(*, i) = Y_GSM_RANGE
            ENDIF 
            IF PLOT_AXIS(i) EQ 'Z_GSM' THEN BEGIN 
                data_pos(*, i) = data(*, 25) 
                range(*, i) = Z_GSM_RANGE
            ENDIF 
            PRINT, '    '+PLOT_AXIS(i)+'     '
        ENDFOR   
        X_RANGE = range(*, 0)
        Y_RANGE = range(*, 1)
        Z_RANGE = range(*, 2)

;--------------------------------------------------------------------------------------------
;                                      Calculation and plot                                 ;
;--------------------------------------------------------------------------------------------

; run for different directions as set in direction_set
        FOR idirection = 0, N_ELEMENTS(direction_set)-1 DO BEGIN 
            flag = flag_phase
            direction = direction_set(idirection)
; the flags of data out of direction are set as no event            
            IF direction EQ 'tailward' THEN flag(where(flag EQ -1)) = !VALUES.F_NAN
            IF direction EQ 'earthward' THEN flag(where(flag EQ 1)) = !VALUES.F_NAN
            IF direction EQ 'both' THEN flag(where(ABS(flag) LE 1 )) = !VALUES.F_NAN             
            spawn, 'mkdir '+ plot_path+direction+'/'            
; run for different grids as set in grid_set
            FOR i_grid = 0, N_ELEMENTS(grid_set)-1 DO BEGIN 
                grid = grid_set(i_grid) & GRID_STR = STRING(GRID, FORMAT = '(i1)')
                
                path = plot_path+direction+'/'+'grid_'+grid_str+'/' & spawn, 'mkdir '+ path
                path = path+storm_phase_set(i_storm_phase)+'/' &  spawn, 'mkdir '+ path
                
                IF KEYWORD_SET(slice_plot) OR KEYWORD_SET(waterdrop_plot) THEN $
                  nslice =  N_ELEMENTS(slice_grid_set) ELSE nslice = 1

;------------------------- POINT PLOT -------------------------
                IF KEYWORD_SET(point_plot) THEN BEGIN
                    path_ev = path+'events/' &  spawn, 'mkdir '+ path_ev
;--- Calculation ---
;classify data into different catergary
                    no_events = DBLARR(ntime, 3) 
                    lobe_t = DBLARR(ntime, 3) 
                    bl_t = DBLARR(ntime, 3) 
                    ps_t = DBLARR(ntime, 3)
                    lobe_e = DBLARR(ntime, 3) 
                    bl_e = DBLARR(ntime, 3) 
                    ps_e = DBLARR(ntime, 3)
                    lobe_2 = DBLARR(ntime, 3) 
                    bl_2 = DBLARR(ntime, 3) 
                    ps_2 = DBLARR(ntime, 3)

                    FOR i = 0l, ntime-1 DO BEGIN                  
                        no_events(i, *) =  data_pos(i, *) * (flag(i) EQ 0)
                        
                        lobe_t(i, *) = data_pos(i, *) *(flag(i) EQ 1)*(beta(i) LE 0.05)
                        bl_t(i, *) = data_pos(i, *) * (flag(i) EQ 1) *  $
                                     (beta(i) LE 1) * (beta(i) GT 0.05)
                        ps_t(i, *) = data_pos(i, *) * (flag(i) EQ 1)*(beta(i) GT 1)
                        
                        lobe_e(i, *) = data_pos(i, *) *(flag(i) EQ -1)*(beta(i) LE 0.05)
                        bl_e(i, *) = data_pos(i, *) * (flag(i) EQ -1) * $
                                     (beta(i) LE 1) * (beta(i) GT 0.05)
                        ps_e(i, *) = data_pos(i, *) * (flag(i) EQ -1)*(beta(i) GT 1)
                        
                        lobe_2(i, *) = data_pos(i, *) * (flag(i) EQ 2) * (beta(i) LE 0.05)
                        bl_2(i, *) = data_pos(i, *) * (flag(i) EQ 2) * $
                                     (beta(i) LE 1) * (beta(i) GT 0.05)
                        ps_2(i, *) = data_pos(i, *) * (flag(i) EQ 2) * (beta(i) GT 1)
                    ENDFOR     
;  IDL
                    IF KEYWORD_SET(IDL_PLOT)THEN BEGIN
                        window, /free
                        PLOT, [0, 0, -100, 100], [-100, 100, 0, 0], $
                              title = sort_title+'         '+PHASE $
                              +'       O!U+!N BEAM EVENS!Cfrom ' + ts_date+ $
                              '  to  ' +te_date,  $
                              xtitle = PLOT_AXIS(0), $
                              ytitle = PLOT_AXIS(1),  $
                              xrange = X_RANGE, yrange = Y_RANGE, $
                              XSTYLE = 1, ystyle = 1, charsize = 1.2, $
                              position = [0.1, 0.1, 0.9, 0.9]
                                ; data with no beam observed                                
                        oplot, no_events(*, 0), no_events(*, 1), color = 5, psym = 3
                                ; tailward events
                        oplot, lobe_t(*, 0), lobe_t(*, 1), color = 3, psym = 3           
                        oplot, bl_t(*, 0), bl_t(*, 1), color = 2, psym = 3
                        oplot, ps_t(*, 0), ps_t(*, 1), color = 1, psym = 3
                                ; earthward events
                        oplot, lobe_e(*, 0), lobe_e(*, 1), color = 3, psym = 3
                        oplot, bl_e(*, 0), bl_e(*, 1), color = 2, psym = 3
                        oplot, ps_e(*, 0), ps_e(*, 1), color = 1, psym = 3
                                ; events with beam in both directions
                        oplot, lobe_2(*, 0), lobe_2(*, 1), color = 3, psym = 3
                        oplot, bl_2(*, 0), bl_2(*, 1), color = 2, psym = 3
                        oplot, ps_2(*, 0), ps_2(*, 1), color = 1, psym = 3
                                ;legend
                        xyouts, la_x, la_y-1.5*grid, 'no events', color = 5
                        xyouts, la_x, la_y, 'lobe', color = 3
                        xyouts, la_x, la_y+1.5*grid, 'boundary layer', color = 2
                        xyouts, la_x, la_y+3*grid, 'plasma sheet', color = 1
                                ;xyouts,la_x, la_y+4.5*grid, 'double beams in lobe', color = 4    
                                ;xyouts,la_x,la_y+6*grid,'double beams in boundary layer',color=5
                                ;xyouts,la_x,la_y+7.5*grid,'double beams in plasma sheet',color=6 
                    ENDIF 
; ps 
                    IF KEYWORD_SET(ps_plot) THEN BEGIN 
                        popen,  path_ev+'events_' $
                                + ts_date+'_to_'+ te_date+'_' $
                                + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                +'.ps', /land
                        PLOT, [0, 0, -100, 100], [-100, 100, 0, 0], $
                              title = sort_title+'         '+PHASE+ '       O!U+!N  BEAM  EVENS    ' $
                              +'!Cfrom  ' +ts_date+'  to  ' +te_date,  $
                              xtitle = PLOT_AXIS(0), ytitle = PLOT_AXIS(1),  $
                              xrange = X_RANGE, yrange = Y_RANGE, $
                              XSTYLE = 1, ystyle = 1, charsize = 2, $
                              position = [0.15, 0.15, 0.85, 0.85]
                                ; no events
                                ;       oplot, no_events(*, 0), no_events(*, 1), color = 2, psym = 3
                                ;               ; tailward events
                        oplot, lobe_t(*, 0), lobe_t(*, 1), color = 2, psym = 1
                        oplot, bl_t(*, 0), bl_t(*, 1), color = 2, psym = 1
                        oplot, ps_t(*, 0), ps_t(*, 1), color = 2, psym = 1     
                                ; earthward events
                        oplot, lobe_e(*, 0), lobe_e(*, 1), color = 2, psym = 1
                        oplot, bl_e(*, 0), bl_e(*, 1), color = 2, psym = 1
                        oplot, ps_e(*, 0), ps_e(*, 1), color = 2, psym = 1    
                                ; events with beam in both directions
                        oplot, lobe_2(*, 0), lobe_2(*, 1), color = 2, psym = 1
                        oplot, bl_2(*, 0), bl_2(*, 1), color = 2, psym = 1
                        oplot, ps_2(*, 0), ps_2(*, 1), color = 2, psym = 1
                                ; legend         
                                ;           xyouts, la_x, la_y, 'lobe', color = 3
                                ;          xyouts, la_x, la_y+1.5*grid, 'boundary layer', color = 2
                                ;         xyouts, la_x, la_y+3*grid, 'plasma sheet', color = 1
                                ;xyouts,la_x,la_y+4.5*grid,'double beams in lobe', color = 4
                                ;xyouts,la_x,la_y+6*grid,'double beams in boundary layer',color=5
                                ;xyouts,la_x,la_y+7.5*grid,'double beams in plasma sheet',color=6
                                ;        xyouts, la_x, la_y-1.5*grid, 'no events', color = 5
                                ; draw the grid lines
                                ;                       FOR i = min(x_range), max(x_range) DO BEGIN 
                                ;                          oplot, [i, i], [-100, 100]
                                ;                     ENDFOR 
                                ;                    FOR i = min(y_range), max(y_range) DO BEGIN 
                                ;                       oplot, [-100, 100], [i, i]
                                ;                  ENDFOR 
                        OPLOT, [0, 0, -100, 100], [-100, 100, 0, 0], col = 0, thick = 8
                        pclose
                        IF keyword_set(mogrify_ps) THEN BEGIN  
                            spawn, 'mogrify -format png '+ path_ev +'*.ps'
                            spawn, 'mogrify -rotate -90 '+ path_ev +'*.png'
                        ENDIF 
                    ENDIF    
                ENDIF   

; run for different slice cut
                FOR i_slice_grid = 0, nslice-1 DO BEGIN 
                    slice_grid = slice_grid_set(i_slice_grid)
                    slice_grid_str = STRING(slice_grid, format = '(i2.2)')

                                ;     stop          
;--------------------------------- EVENTS MAP ---------------------
                    IF KEYWORD_SET(EVENTS_MAP) THEN BEGIN 
                        path_ev = path+'events/' &  spawn, 'mkdir '+ path_ev
;-----  Calculation ------  
                        nx = CEIL(ABS(x_range(1) - x_range(0))/grid)
                        ny = CEIL(ABS(y_range(1) - y_range(0))/grid)
                        nz = CEIL(ABS(z_range(1)-z_range(0))/SLICE_GRID) 
                        x_axis = INDGEN(nx)*grid+(x_range(0) < x_range(1))+grid*0.5
                        y_axis = INDGEN(ny)*grid+(y_range(0) < y_range(1))+grid*0.5
                        z_axis = INDGEN(nz)*slice_grid+(z_range(0) < z_range(1))+ $
                                 slice_grid*0.5
                        
                        event_counts = FLTARR(nx, ny, nz)
                        total_counts = FLTARR(nx, ny, nz)
                        lobe_counts =  FLTARR(nx, ny, nz) 
                        lobe_total =  FLTARR(nx, ny, nz) 
                        bl_counts =  FLTARR(nx, ny, nz) 
                        bl_total =  FLTARR(nx, ny, nz) 
                        ps_counts = FLTARR(nx, ny, nz)       
                        ps_total =  FLTARR(nx, ny, nz) 
                        
                        FOR i = 0l, ntime-1 DO BEGIN 
                            index = SORT(ABS(x_axis - data_pos(i, 0))) 
                            index_x = index(0)
                            index = SORT(ABS(y_axis - data_pos(i, 1)))
                            index_y = index(0)
                            index = SORT(ABS(z_axis - data_pos(i, 2)))
                            index_z = index(0)
                                ;all 
                            total_counts(index_x, index_y, index_z) =  $
                              total_counts(index_x, index_y, index_z)+(ABS(flag(i)) GE 0)
                            event_counts(index_x, index_y, index_z) =  $
                              event_counts(index_x, index_y, index_z )+ (ABS(flag(i)) GE 1) 
                                ;lobe
                            lobe_total(index_x, index_y, index_z) = $
                              lobe_total(index_x, index_y, index_z)+ $
                              (beta(i) LE 0.05 AND ABS(flag(i))GE 0) 
                            lobe_counts(index_x, index_y, index_z) = $
                              lobe_counts(index_x, index_y, index_z)+ $
                              (beta(i) LE 0.05 AND ABS(flag(i))GE 1)                 
                                ; boundary layer
                            bl_total(index_x, index_y, index_z) =  $
                              bl_total(index_x, index_y, index_z)+ $
                              (beta(i) GT 0.05 AND beta(i) LE 1 AND ABS(flag(i))GE 0)
                            bl_counts(index_x, index_y, index_z) = $
                              bl_counts(index_x, index_y, index_z)+ $
                              (beta(i) GT 0.05 AND beta(i) LE 1 AND ABS(flag(i))GE 1)     
                                ; plasma sheet
                            ps_total(index_x, index_y, index_z) = $
                              ps_total(index_x, index_y, index_z) $
                              + (beta(i) GT 1 AND ABS(flag(i))GE 0)
                            ps_counts(index_x, index_y, index_z) = $
                              ps_counts(index_x, index_y, index_z) $
                              +(beta(i) GT 1 AND ABS(flag(i))GE 1) 
                        ENDFOR 
                        event_ratio = event_counts/total_counts
                        lobe_ratio = lobe_counts/lobe_total
                        bl_ratio = bl_counts/bl_total
                        ps_ratio = ps_counts/ps_total                         
;------ 2D plot --------
                        IF KEYWORD_SET(PLOT_2D) THEN BEGIN     
                            path_2d = path_ev+'2d/' 
                            spawn, 'mkdir '+path_2d
                            
                            event_counts_2d =  TOTAL(event_counts, 3, /nan)
                            lobe_counts_2d = TOTAL(lobe_counts, 3, /nan)
                            bl_counts_2d =  TOTAL(bl_counts, 3, /nan)
                            ps_counts_2d =  TOTAL(ps_counts, 3, /nan)
                            
                            event_counts_2d(where(event_counts_2d EQ 0)) =  !VALUES.F_NAN
                            lobe_counts_2d(where(lobe_counts_2d EQ 0)) = !VALUES.F_NAN
                            bl_counts_2d(where(bl_counts_2d EQ 0)) = !VALUES.F_NAN
                            ps_counts_2d(where(ps_counts_2d EQ 0)) = !VALUES.F_NAN

                            event_ratio_2d = TOTAL(event_counts, 3, /nan)/TOTAL(total_counts, 3, /nan)
                            lobe_ratio_2d = TOTAL(lobe_counts, 3, /nan)/TOTAL(lobe_total, 3, /nan)
                            bl_ratio_2d = TOTAL(bl_counts, 3, /nan)/TOTAL(bl_total, 3, /nan)
                            ps_ratio_2d = TOTAL(ps_counts, 3, /nan)/TOTAL(ps_total, 3, /nan)
                            
                            IF keyword_set(symmetry_template) THEN BEGIN
                                total_counts_2d = TOTAL(total_counts, 3, /nan)
                                positive_half = total_counts_2d(*, where(y_axis GT 0)) GT  0
                                positive_half = [[reverse(positive_half, 2)], [positive_half]]
                                negative_half = total_counts_2d(*, where(y_axis LT 0)) GT  0
                                negative_half = [[negative_half], [reverse(negative_half, 2)]]
                                template_2d = FLOAT(positive_half*negative_half)

                                template_2d(where(template_2d EQ 0 )) =  !VALUES.F_NAN

                                event_counts_2d = event_counts_2d*template_2d
                                lobe_counts_2d = lobe_counts_2d*template_2d
                                bl_counts_2d = bl_counts_2d*template_2d
                                ps_counts_2d = ps_counts_2d*template_2d

                                event_ratio_2d = event_ratio_2d*template_2d
                                lobe_ratio_2d = lobe_ratio_2d*template_2d
                                bl_ratio_2d = bl_ratio_2d*template_2d
                                ps_ratio_2d = ps_ratio_2d*template_2d
                            ENDIF                     
;     idl
                            IF KEYWORD_SET(IDL_PLOT) THEN BEGIN 
                                window, /free ; all events
                                specplot, x_axis, y_axis, event_counts_2d, $
                                          no_interp = 1, $
                                          lim = { zlog:EVENTS_V_LOG, $
                                                  zrange: EVENTS_V_RANGE, $
                                                  title: sort_title+'         '+PHASE+$
                                                  '       O!U+!N BEAM EVENTS !CFROM  ' $
                                                  + ts_date+'  TO  ' +te_date,  $
                                                  xtitle: PLOT_AXIS(0), $
                                                  ytitle: PLOT_AXIS(1), $
                                                  xrange: X_RANGE, yrange: Y_RANGE, $
                                                  XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                  position: [0.1, 0.1, 0.9, 0.9]}   
                                oplot, [0, 0, -100, 100], [-100, 100, 0, 0]

                                IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                    window, /free ;lobe events
                                    specplot, x_axis, y_axis,  lobe_counts_2d, $
                                              no_interp = 1, $
                                              lim = { zlog:EVENTS_V_LOG, $
                                                      zrange: EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM EVENTS in LOBE!CFROM  '$
                                                      + ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    window, /free ; bl events
                                    specplot, x_axis, y_axis,  bl_counts_2d, $
                                              no_interp = 1, $
                                              lim = { zlog:EVENTS_V_LOG, $
                                                      zrange: EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM EVENTS  in BL!CFROM  ' $
                                                      + ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    window, /free ;ps events
                                    specplot, x_axis, y_axis,  ps_counts_2d, $
                                              no_interp = 1, $
                                              lim = { zlog:EVENTS_V_LOG, $
                                                      zrange: EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM EVENTS  in PS!CFROM  ' $
                                                      + ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                ENDIF  
                            ENDIF 
; ps
                            IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                popen,  path_2d+'all_events_' $ ; all events
                                        +ts_date+'_to_' + te_date+'_' $
                                        +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                        +'.ps', /land 
                                specplot, x_axis, y_axis,  event_counts_2d, $
                                          no_interp = 1, $
                                          lim = { zlog:EVENTS_V_LOG, zrange: EVENTS_V_RANGE, $
                                                  title: sort_title+'         '+PHASE+$
                                                  '       O!U+!N BEAM EVENTS !Cfrom '$
                                                  + ts_date+'  to  ' +te_date,  $
                                                  xtitle: PLOT_AXIS(0),  $
                                                  ytitle: PLOT_AXIS(1), $
                                                  xrange: X_RANGE, yrange: Y_RANGE, $
                                                  XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                  position: [0.1, 0.1, 0.9, 0.9]}     
                                oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                pclose 
                                IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                    popen,  path_2d+'lobe_events_'+ $ ;lobe events
                                            ts_date+'_to_' + te_date+'_'+ $
                                            PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                            +'.ps', /land 
                                    specplot, x_axis, y_axis, lobe_counts_2d, $
                                              no_interp = 1, $
                                              lim = { zlog: EVENTS_V_LOG, $
                                                      zrange: EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM EVENTS in LOBE!Cfrom  '$
                                                      + ts_date+'  to  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    pclose 
                                    popen,  path_2d+'bl_events_'+ $ ; bl events
                                            ts_date+'_to_' + te_date+'_'+ $
                                            PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                            +'.ps', /land 
                                    specplot, x_axis, y_axis,  bl_counts_2d, $
                                              no_interp = 1, $
                                              lim = { zlog: EVENTS_V_LOG, $
                                                      zrange: EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM EVENTS  in BL!CFROM  '$
                                                      + ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    pclose
                                    popen,  path_2d+'ps_events_'+ $ ;ps events
                                            ts_date+'_to_' + te_date+'_'+ $
                                            PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                            +'.ps', /land 
                                    specplot, x_axis, y_axis, ps_counts_2d, $
                                              no_interp = 1, $
                                              lim = { zlog:EVENTS_V_LOG, $
                                                      zrange: EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+ $
                                                      '       O!U+!N BEAM EVENTS  in PS!CFROM  ' $
                                                      + ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    pclose
                                ENDIF 
                                IF keyword_set(mogrify_ps) THEN BEGIN  
                                    spawn, 'mogrify -format png '+ path_2d +'*.ps'
                                    spawn, 'mogrify -rotate -90 '+ path_2d +'*.png'
                                ENDIF 
                            ENDIF                                
; idl
                            IF KEYWORD_SET(IDL_PLOT) THEN BEGIN 
                                window, /free ; all ratio
                                specplot, x_axis, y_axis, event_ratio_2d, $
                                          no_interp = 1, $
                                          lim = { zlog:RATIO_V_LOG, $
                                                  zrange:RATIO_V_RANGE, $
                                                  title: sort_title+'         '+PHASE+$
                                                  '       O!U+!N BEAM RATIO!CFROM  ' + $
                                                  ts_date+'  TO  ' +te_date,  $
                                                  xtitle: PLOT_AXIS(0), $
                                                  ytitle: PLOT_AXIS(1), $
                                                  xrange: X_RANGE, yrange: Y_RANGE, $
                                                  XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                  position: [0.1, 0.1, 0.9, 0.9]}   
                                oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                    window, /free ;lobe ratio
                                    specplot, x_axis, y_axis, lobe_ratio_2d, $
                                              no_interp = 1, $
                                              lim = { zlog: RATIO_V_LOG, $
                                                      zrange: RATIO_V_RANGE, $
                                                      title:sort_title+'         '+PHASE+ $
                                                      '       O!U+!N BEAM RATIO  in LOBE!CFROM '+$
                                                      ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    window, /free ; bl ratio
                                    specplot, x_axis, y_axis, bl_ratio_2d, $
                                              no_interp = 1, $
                                              lim = { zlog: RATIO_V_LOG, $
                                                      zrange: RATIO_V_RANGE, $
                                                      title:sort_title+'         '+PHASE+ $
                                                      '       O!U+!N BEAM RATIO  in BL!CFROM ' + $
                                                      ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    window, /free ; ps ratio
                                    specplot, x_axis, y_axis, ps_ratio_2d, $
                                              no_interp = 1, $
                                              lim = { zlog: RATIO_V_LOG, $
                                                      zrange:RATIO_V_RANGE, $
                                                      title:sort_title+'         '+PHASE+$
                                                      '      O!U+!N BEAM RATIO  in PS!CFROM ' + $
                                                      ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                ENDIF
                            ENDIF 
; ps
                            IF KEYWORD_SET(PS_PLOT) THEN BEGIN 
                                popen, path_2d+'all_ratio_' $ ; all ratio
                                       +ts_date+'_to_' + te_date+'_'+ $
                                       PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                       +'.ps', /land 
                                specplot, x_axis, y_axis, event_ratio_2d, $
                                          no_interp = 1, $
                                          lim = { zlog:RATIO_V_LOG, $
                                                  zrange: RATIO_V_RANGE, $
                                                  title: sort_title+'         '+PHASE+ $
                                                  '       O!U+!N BEAM RATIO!CFROM  ' + $
                                                  ts_date+'  TO  ' +te_date,  $
                                                  xtitle: PLOT_AXIS(0), $
                                                  ytitle: PLOT_AXIS(1), $
                                                  xrange: X_RANGE, yrange: Y_RANGE, $
                                                  XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                  position: [0.1, 0.1, 0.9, 0.9]}   
                                oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                pclose   
                                IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                    popen,  path_2d+'lobe_ratio_' $ ; lobe ratio
                                            +ts_date+'_to_' +te_date+'_' $
                                            + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                            +'.ps', /land 
                                    specplot, x_axis, y_axis, lobe_ratio_2d, $
                                              no_interp = 1, $
                                              lim = { zlog:RATIO_V_LOG, $
                                                      zrange:RATIO_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+ $
                                                      '       O!U+!N BEAM RATIO  in LOBE!CFROM '+$
                                                      ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    pclose
                                    popen,  path_2d+'bl_ratio_' $ ;bl ratio
                                            + ts_date+'_to_' + te_date+'_' $
                                            +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                            +'.ps', /land 
                                    specplot, x_axis, y_axis, bl_ratio_2d, $
                                              no_interp = 1, $
                                              lim = { zlog: RATIO_V_LOG, $
                                                      zrange: RATIO_V_RANGE, $
                                                      title:sort_title+'         '+PHASE+ $
                                                      '       O!U+!N BEAM RATIO  in BL!CFROM ' + $
                                                      ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    pclose
                                    popen,  path_2d+'ps_ratio_' $ ;ps ratio
                                            +ts_date+'_to_' + te_date+'_' $
                                            +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                            +'.ps', /land 
                                    specplot, x_axis, y_axis, ps_ratio_2d, $
                                              no_interp = 1, $
                                              lim = { zlog:RATIO_V_LOG, $
                                                      zrange:RATIO_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+ $
                                                      '       O!U+!N BEAM RATIO  in PS!CFROM ' + $
                                                      ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    pclose
                                ENDIF 
                                IF keyword_set(mogrify_ps) THEN BEGIN  
                                    spawn, 'mogrify -format png '+ path_2d +'*.ps'
                                    spawn, 'mogrify -rotate -90 '+ path_2d +'*.png'
                                ENDIF 
                            ENDIF     
                        ENDIF                                 
;--------- slice-plots -----------
; z is the direction been sliced
                        IF KEYWORD_SET(SLICE_PLOT) THEN BEGIN 
                            path_slice_main = path_ev+'slice/'
                            spawn, 'mkdir ' +path_slice_main
;Calulation
                            FOR iz = 0, nz-1 DO BEGIN 
                                slice_block = [STRING(z_axis(iz)-SLICE_GRID*0.5, $
                                                      format = '(f5.1)'), $
                                               STRING(z_axis(iz)+SLICE_GRID*0.5, $
                                                      format = '(f5.1)')]

                                slice_event_counts = event_counts(*, *, iz)
                                slice_lobe_counts = lobe_counts(*, *, iz)
                                slice_bl_counts = bl_counts(*, *, iz)
                                slice_ps_counts = ps_counts(*, *, iz)
                                
                                slice_event_ratio = event_ratio(*, *, iz)
                                slice_lobe_ratio = lobe_ratio(*, *, iz)
                                slice_bl_ratio = bl_ratio(*, *, iz)
                                slice_ps_ratio = ps_ratio(*, *, iz)
                                
                                slice_event_counts(where(slice_event_counts EQ 0)) $
                                  = !VALUES.F_NAN
                                slice_lobe_counts(where(slice_lobe_counts EQ 0)) $
                                  = !VALUES.F_NAN
                                slice_bl_counts(where(slice_bl_counts EQ 0)) = !VALUES.F_NAN
                                slice_ps_counts(where(slice_ps_counts EQ 0)) = !VALUES.F_NAN

                                IF keyword_set(symmetry_template) THEN BEGIN
                                    slice_total_counts = total_counts(*, *, iz)
                                    positive_half = slice_total_counts(*, where(y_axis GT 0)) GT 0
                                    positive_half = [[reverse(positive_half, 2)], [positive_half]]
                                    negative_half = slice_total_counts(*, where(y_axis LT 0)) GT 0
                                    negative_half = [[negative_half], [reverse(negative_half, 2)]]
                                    template_slice = FLOAT(positive_half*negative_half)
                                    
                                ;                             window, /free ; all events
                                ;                            specplot, x_axis, y_axis, template_slice, $
                                ;                                     no_interp = 1, $
                                ;                                    lim = { zlog:0, $
                                ;                                           zrange: [0, 1], $
                                ;                                          title: sort_title+'         '+PHASE+$
                                ;                                         '       O!U+!N BEAM EVENTS !CFROM  ' $
                                ;                                        + ts_date+'  TO  ' +te_date,  $
                                ;                                       xtitle: PLOT_AXIS(0), $
                                ;                                      ytitle: PLOT_AXIS(1), $
                                ;                                     xrange: x_range, yrange: Y_RANGE, $
                                ;                                    XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                ;                                   position: [0.1, 0.1, 0.9, 0.9]}   
                                    template_slice(where(template_slice EQ 0 )) =  !VALUES.F_NAN
                                    slice_event_counts = slice_event_counts*template_slice
                                    slice_lobe_counts = slice_lobe_counts*template_slice
                                    slice_bl_counts = slice_bl_counts*template_slice
                                    slice_ps_counts = slice_ps_counts*template_slice

                                    slice_event_ratio = slice_event_ratio*template_slice
                                    slice_lobe_ratio = slice_lobe_ratio*template_slice
                                    slice_bl_ratio = slice_bl_ratio*template_slice
                                    slice_ps_ratio = slice_ps_ratio*template_slice                                 
                                ENDIF       
;   idl
                                IF KEYWORD_SET(idl_plot) THEN BEGIN 
                                    window, /free ; all events
                                    specplot, x_axis, y_axis,  slice_event_counts, $
                                              no_interp = 1, $
                                              lim = { zlog:EVENTS_V_LOG, $
                                                      zrange: EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM EVENTS at  ' $
                                                      + plot_axis(2) $
                                                      + ':( ' +slice_block(0)+',' $
                                                      + slice_block(1)+' ) ' $
                                                      + '!CFROM  ' + ts_date+'  TO  ' $
                                                      +te_date, $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    
                                    window, /free ;all ratio
                                    specplot, x_axis, y_axis,  slice_event_ratio, $
                                              no_interp = 1, $
                                              lim = { zlog:RATIO_V_LOG, $
                                                      zrange: RATIO_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM RATIO at  ' $
                                                      + plot_axis(2) $
                                                      + ':( ' +slice_block(0)+',' $
                                                      +slice_block(1)+' ) ' $
                                                      + '!CFROM  ' + ts_date+'  TO  ' $
                                                      +te_date, $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]     
                                    
                                    IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                        window, /free ; lobe event
                                        specplot, x_axis, y_axis,  slice_lobe_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog: EVENTS_V_LOG, $
                                                          zrange: EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+',' $
                                                          +slice_block(1)+' ) ' $
                                                          + 'in LOBE'$
                                                          + '!CFROM  ' + ts_date+'  TO  ' $
                                                          +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        window, /free ;lobe ratio
                                        specplot, x_axis, y_axis,  slice_lobe_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog: RATIO_V_LOG, $
                                                          zrange:RATIO_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+',' $
                                                          +slice_block(1)+' ) ' $
                                                          + 'in LOBE'$
                                                          + '!CFROM  ' + ts_date+'  TO  ' $
                                                          +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]     
                                        window, /free ;bl events
                                        specplot, x_axis, y_axis,  slice_bl_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog:EVENTS_V_LOG, $
                                                          zrange:EVENTS_V_RANGE,  $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+','$
                                                          +slice_block(1)+' ) ' $
                                                          + 'in BL'$
                                                          + '!CFROM  ' + ts_date+'  TO  ' $
                                                          +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        window, /free ;bl ratio
                                        specplot, x_axis, y_axis,  slice_bl_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog:RATIO_V_LOG, $
                                                          zrange:RATIO_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+','$
                                                          +slice_block(1)+' ) ' $
                                                          + 'in BL'$
                                                          + '!CFROM  ' + ts_date+'  TO  ' $
                                                          +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]     
                                        window, /free ;plasma sheet events
                                        specplot, x_axis, y_axis,  slice_ps_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog:EVENTS_V_LOG, $
                                                          zrange:EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+','$
                                                          +slice_block(1)+' ) ' $
                                                          + 'in PS'$
                                                          + '!CFROM  ' + ts_date+'  TO  ' $
                                                          +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        window, /free ;plasma sheet ratio
                                        specplot, x_axis, y_axis,  slice_ps_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog: RATIO_V_LOG, $
                                                          zrange:RATIO_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+','$
                                                          +slice_block(1)+' ) ' $
                                                          + 'in PS'$
                                                          + '!CFROM  ' + ts_date+'  TO  ' $
                                                          +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]     
                                    ENDIF  
                                ENDIF 
;  ps
                                IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                    path_slice = path_slice_main $
                                                 +'slice_zgrid_'+slice_grid_str+'/'
                                    spawn, 'mkdir '+ path_slice
                                    popen, path_slice+'all_events_' $ ; all events
                                           + ts_date+'_to_' + te_date+'_' $
                                           + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $ 
                                           +'_at_' + plot_axis(2)+'_' $
                                           +slice_block(0)+'_'+slice_block(1) $
                                           +'.ps', /land 
                                    specplot, x_axis, y_axis,  slice_event_counts $
                                              , no_interp = 1, $
                                              lim = { zlog:EVENTS_V_LOG, $
                                                      zrange:EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM EVENTS at  ' $
                                                      + plot_axis(2) $
                                                      + ':( ' +slice_block(0)+','$
                                                      +slice_block(1)+' ) ' $
                                                      + '!CFROM  ' + ts_date+'  TO  ' $
                                                      +te_date, $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    pclose
                                    popen,  path_slice+'all_ratio_' $ ; all ratio      
                                            +ts_date+'_to_' $
                                            + te_date+'_'+PLOT_AXIS(0)+'_vs_'$
                                            +PLOT_AXIS(1)+'_at_' $
                                            + plot_axis(2)+'_'+$
                                            slice_block(0)+'_'+slice_block(1) $
                                            +'.ps', /land 
                                    specplot, x_axis, y_axis,  slice_event_ratio, $
                                              no_interp = 1, $
                                              lim = { zlog: RATIO_V_LOG, $
                                                      zrange: RATIO_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM RATIO at  ' $
                                                      + plot_axis(2) $
                                                      + ':( ' +slice_block(0)+',' $
                                                      +slice_block(1)+' ) ' $
                                                      + '!CFROM  ' + ts_date $
                                                      +'  TO  ' +te_date, $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]    
                                    pclose          
                                    IF KEYWORD_SET(diff_beta) THEN BEGIN     
                                        popen,  path_slice+'lobe_events_' $ ;lobe events
                                                +ts_date+'_to_' $
                                                + te_date+'_'$
                                                +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1)+'_at_' $
                                                + plot_axis(2)+'_'$
                                                +slice_block(0)+'_'+slice_block(1) +'.ps', /land 
                                        specplot, x_axis, y_axis,  slice_lobe_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog:EVENTS_V_LOG, $
                                                          zrange:EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+','$
                                                          +slice_block(1)+' ) ' $
                                                          + 'in LOBE'$
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        popen,  path_slice+'lobe_ratio_' $ ; lobe ratio        
                                                +ts_date+'_to_' $
                                                + te_date+'_'+PLOT_AXIS(0)+'_vs_'$
                                                +PLOT_AXIS(1)+'_at_' $
                                                + plot_axis(2)+'_' $
                                                +slice_block(0)+'_'+slice_block(1) +'.ps', /land 
                                        specplot, x_axis, y_axis,  slice_lobe_ratio $
                                                  , no_interp = 1, $
                                                  lim = { zlog:RATIO_V_LOG, $
                                                          zrange:RATIO_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+ $
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          + 'in LOBE'$
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]  
                                        pclose
                                        popen,   path_slice+'bl_events_' $ ;boundary layer events
                                                 +ts_date+'_to_' $
                                                 + te_date+'_'+PLOT_AXIS(0)$
                                                 +'_vs_'+PLOT_AXIS(1)+'_at_' $
                                                 + plot_axis(2)+'_'+slice_block(0)$
                                                 +'_'+slice_block(1) +'.ps', /land 
                                        specplot, x_axis, y_axis,  slice_bl_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog:EVENTS_V_LOG, $
                                                          zrange:EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          + 'in BL'$
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        popen,   path_slice+'bl_ratio_' $ ; boundary layer ratio 
                                                 +ts_date+'_to_' $
                                                 + te_date+'_'+PLOT_AXIS(0)+'_vs_' $
                                                 +PLOT_AXIS(1)+'_at_' $
                                                 + plot_axis(2)+'_' $
                                                 +slice_block(0)+'_'+slice_block(1) +'.ps', /land 
                                        specplot, x_axis, y_axis,  slice_bl_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog:RATIO_V_LOG, $
                                                          zrange: RATIO_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          + 'in BL'$
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]   
                                        pclose
                                        popen,  path_slice+'ps_events_' $ ;plasma sheet events
                                                +ts_date+'_to_' $
                                                + te_date+'_'+PLOT_AXIS(0)+'_vs_' $
                                                +PLOT_AXIS(1)+'_at_' $
                                                + plot_axis(2)+'_' $
                                                +slice_block(0)+'_'+slice_block(1) +'.ps', /land 
                                        specplot, x_axis, y_axis,  slice_ps_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog:EVENTS_V_LOG, $
                                                          zrange:EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+',' $
                                                          +slice_block(1)+' ) ' $
                                                          + 'in PS'$
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        popen, path_slice+'ps_ratio_' $ ;    plasma sheet ratio  
                                               +ts_date+'_to_' $
                                               + te_date+'_'+PLOT_AXIS(0)+'_vs_' $
                                               + PLOT_AXIS(1)+'_at_' $
                                               + plot_axis(2)+'_' $
                                               +slice_block(0)+'_'+slice_block(1) +'.ps', /land 
                                        specplot, x_axis, y_axis,  slice_ps_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog:RATIO_V_LOG, $
                                                          zrange:RATIO_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          + 'in PS'$
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0] 
                                        pclose
                                    ENDIF 
                                    IF keyword_set(mogrify_ps) THEN BEGIN 
                                        spawn, 'mogrify -format png '+ path_slice +'*.ps'
                                        spawn, 'mogrify -rotate -90 '+ path_slice +'*.png'
                                    ENDIF 
                                ENDIF    
                            ENDFOR    
                        ENDIF      
; ------- WATERDROP PLOT ------------
;SLICE_GRID WILL BE ALSO USED HERE
                        IF KEYWORD_SET(WATERDROP_PLOT) AND (nz NE 2) THEN BEGIN 
                            path_wd_main = path_ev+'waterdrop/'
                            spawn, 'mkdir ' + path_wd_main
;run for all slices
                            FOR iz = 0, nz-1 DO BEGIN 
; Calculation
                                IF iz LT CEIL(nz/2.) THEN index = INDGEN(CEIL(nz/2.)-iz)+iz  $
                                ELSE index = INDGEN(iz+1-FLOOR(nz/2.))+FLOOR(nz/2.)

                                IF N_ELEMENTS(index) GT 1 THEN BEGIN 
                                    wd_event_counts = TOTAL(event_counts(*, *, index), 3, /nan)
                                    wd_lobe_counts = TOTAL(lobe_counts(*, *, index), 3, /nan)
                                    wd_bl_counts = TOTAL(bl_counts(*, *, index), 3, /nan)
                                    wd_ps_counts = TOTAL(ps_counts(*, *, index), 3, /nan)

                                    wd_event_ratio = TOTAL(event_counts(*, *, index), 3, /nan)/ $
                                                     TOTAL(total_counts(*, *, index), 3, /nan)
                                    wd_lobe_ratio = TOTAL(lobe_counts(*, *, index), 3, /nan)/ $
                                                    TOTAL(lobe_total(*, *, index), 3, /nan)
                                    wd_bl_ratio = TOTAL(bl_counts(*, *, index), 3, /nan)/ $
                                                  TOTAL(bl_total(*, *, index), 3, /nan)
                                    wd_ps_ratio = TOTAL(ps_counts(*, *, index), 3, /nan)/ $
                                                  TOTAL(ps_total(*, *, index), 3, /nan)
                                ENDIF ELSE BEGIN
                                    wd_event_counts = event_counts(*, *, index)
                                    wd_lobe_counts = lobe_counts(*, *, index)
                                    wd_bl_counts = bl_counts(*, *, index)
                                    wd_ps_counts = ps_counts(*, *, index)
                                    
                                    wd_event_ratio = event_ratio(*, *, index)
                                    wd_lobe_ratio = lobe_ratio(*, *, index)
                                    wd_bl_ratio = bl_ratio(*, *, index)
                                    wd_ps_ratio = ps_ratio(*, *, index) 
                                ENDELSE     
                                wd_event_counts(where(wd_event_counts EQ 0)) =  !VALUES.F_NAN
                                wd_lobe_counts(where(wd_lobe_counts EQ 0)) = !VALUES.F_NAN
                                wd_bl_counts(where(wd_bl_counts EQ 0)) = !VALUES.F_NAN
                                wd_ps_counts(where(wd_ps_counts EQ 0)) = !VALUES.F_NAN
                                
                                slice_block = [ STRING(z_axis(index(0))-SLICE_GRID*0.5, $
                                                       format = '(f5.1)'), $
                                                STRING(z_axis(index(N_ELEMENTS(index)-1)) $
                                                       +SLICE_GRID*0.5,  $
                                                       format = '(f5.1)')]
;idl plot
                                IF KEYWORD_SET(idl_plot) THEN BEGIN 
                                    window, /free ; all EVENTS
                                    specplot, x_axis, y_axis, wd_event_counts, $
                                              no_interp = 1, $
                                              lim = { zlog:EVENTS_V_LOG, $
                                                      zrange: EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM EVENTS at  ' $
                                                      + plot_axis(2) $
                                                      + ':( ' +slice_block(0)+',' $
                                                      +slice_block(1)+' ) ' $
                                                      +'in LOBE' $
                                                      + '!CFROM  ' + ts_date+'  TO  ' $
                                                      +te_date, $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    window, /free ; all RATIO
                                    specplot, x_axis, y_axis, wd_event_ratio, $
                                              no_interp = 1, $
                                              lim = { zlog: RATIO_V_LOG, $
                                                      zrange: RATIO_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM RATIO at  ' $
                                                      + plot_axis(2) $
                                                      + ':( ' +slice_block(0) $
                                                      +','+slice_block(1)+' ) ' $
                                                      +'in LOBE' $
                                                      + '!CFROM  ' + ts_date+'  TO  ' $
                                                      +te_date, $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]     
                                    IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                        window, /free ;lobe EVENTS
                                        specplot, x_axis, y_axis, wd_lobe_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog:EVENTS_V_LOG, $
                                                          zrange: EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  '$
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+',' $
                                                          +slice_block(1)+' ) ' $
                                                          +'in LOBE' $
                                                          + '!CFROM  ' + ts_date+'  TO  ' $
                                                          +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        window, /free ; LOBE RATIO
                                        specplot, x_axis, y_axis, wd_lobe_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog:RATIO_V_LOG, $
                                                          zrange:RATIO_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          +'in LOBE' $
                                                          + '!CFROM  ' + ts_date+'  TO  ' $
                                                          +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]     
                                        window, /free ;boundary layer EVENTS
                                        specplot, x_axis, y_axis, wd_bl_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog:EVENTS_V_LOG, $
                                                          zrange: EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          +'in BL' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        window, /free ; boundary layer RATIO
                                        specplot, x_axis, y_axis, wd_bl_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog:RATIO_V_LOG, $
                                                          zrange:RATIO_V_RANGE, $
                                                          title:sort_title+'         '+PHASE+ $
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +',' +slice_block(1)+' ) ' $
                                                          +'in BL' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]     
                                        window, /free ;plasma sheet events
                                        specplot, x_axis, y_axis, wd_ps_counts, $ 
                                                  no_interp = 1, $
                                                  lim = { zlog:EVENTS_V_LOG, $
                                                          zrange: EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          +'in PS' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        window, /free ;plasma sheet ratio
                                        specplot, x_axis, y_axis, wd_ps_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog:RATIO_V_LOG, $
                                                          zrange: RATIO_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          +'in PS' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]     
                                    ENDIF 
                                ENDIF 
;ps plot
                                IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                    path_wd = path_wd_main+'wd_zgrid_' + slice_grid_str +'/'
                                    spawn, 'mkdir ' + path_wd 
                                    
                                    popen,  path_wd+'all_events_'+ $ ; all events
                                            ts_date+'_to_' + te_date+'_' $
                                            +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                            +'_at_'+ plot_axis(2) $
                                            +'_'+slice_block(0)+'_'+slice_block(1) $
                                            +'.ps', /land 
                                    specplot, x_axis, y_axis, wd_event_counts, $
                                              no_interp = 1, $
                                              lim = { zlog: EVENTS_V_LOG, $
                                                      zrange: EVENTS_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM EVENTS at  ' $
                                                      + plot_axis(2) $
                                                      + ':( ' +slice_block(0) $
                                                      +','+slice_block(1)+' ) ' $
                                                      + '!CFROM  ' + ts_date $
                                                      +'  TO  ' +te_date, $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    pclose
                                    popen,   path_wd+'all_ratio_'+ $ ; all ratio
                                             ts_date+'_to_' + te_date+'_' $
                                             +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                             +'_at_'+ plot_axis(2) $
                                             +'_'+slice_block(0)+'_'+slice_block(1) $
                                             +'.ps', /land 
                                    specplot, x_axis, y_axis, wd_event_ratio, $
                                              no_interp = 1, $
                                              lim = { zlog:RATIO_V_LOG, $
                                                      zrange:RATIO_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N BEAM RATIO at  ' $
                                                      + plot_axis(2) $
                                                      + ':( ' +slice_block(0) $
                                                      +','+slice_block(1)+' ) ' $
                                                      + '!CFROM  ' + ts_date $
                                                      +'  TO  ' +te_date, $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]    
                                    pclose              
                                    IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                        popen,   path_wd+'lobe_events_'+ $ ;lobe events
                                                 ts_date+'_to_' + te_date+'_' $
                                                 +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                 +'_at_'+ plot_axis(2) $
                                                 +'_'+slice_block(0)+'_'+slice_block(1) $
                                                 +'.ps', /land  
                                        specplot, x_axis, y_axis, wd_lobe_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog: EVENTS_V_LOG, $
                                                          zrange:  EVENTS_V_RANGE, $
                                                          title:sort_title+'         '+PHASE+ $
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          +'in LOBE' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        popen,   path_wd+'lobe_ratio_'+ $ ;lobe ratio
                                                 ts_date+'_to_' + te_date+'_' $
                                                 +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                 +'_at_'+ plot_axis(2)$
                                                 +'_'+slice_block(0)+'_'+slice_block(1) $
                                                 +'.ps', /land  
                                        specplot, x_axis, y_axis, wd_lobe_ratio, no_interp = 1, $
                                                  lim = { zlog: RATIO_V_LOG, $
                                                          zrange:RATIO_V_RANGE, $
                                                          title:sort_title+'         '+PHASE+ $
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          +'in LOBE' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]  
                                        pclose
                                        popen,  path_wd+'bl_events_'+ $ ;bl events
                                                ts_date+'_to_' + te_date+'_' $
                                                +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                +'_at_'+ plot_axis(2)+'_' $
                                                +slice_block(0)+'_'+slice_block(1) $
                                                +'.ps', /land  
                                        specplot, x_axis, y_axis, wd_bl_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog: EVENTS_V_LOG, $
                                                          zrange:  EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          +'in BL' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        popen,  path_wd+'bl_ratio_'+ $ ;bl ratio
                                                ts_date+'_to_' + te_date+'_' $
                                                +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                +'_at_'+ plot_axis(2)+'_' $
                                                +slice_block(0)+'_'+slice_block(1) $
                                                +'.ps', /land  
                                        specplot, x_axis, y_axis, wd_bl_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog:RATIO_V_LOG, $
                                                          zrange:RATIO_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          +'in BL' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]   
                                        pclose
                                        popen,  path_wd+'ps_events_'+ $ ;plasma sheet events
                                                ts_date+'_to_' + te_date+'_' $
                                                +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                +'_at_'+ plot_axis(2)+'_' $
                                                +slice_block(0)+'_'+slice_block(1) $
                                                +'.ps', /land  
                                        specplot, x_axis, y_axis, wd_ps_counts, $
                                                  no_interp = 1, $
                                                  lim = { zlog: EVENTS_V_LOG, $
                                                          zrange: EVENTS_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N BEAM EVENTS at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+',' $
                                                          +slice_block(1)+' ) ' $
                                                          +'in PS' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        popen,  path_wd+'ps_ratio_'+ $ ;plasma sheet ratio
                                                ts_date+'_to_' + te_date+'_' $
                                                +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                +'_at_'+ plot_axis(2)+'_' $
                                                +slice_block(0)+'_'+slice_block(1) $
                                                +'.ps', /land  
                                        specplot, x_axis, y_axis, wd_ps_ratio, $
                                                  no_interp = 1, $
                                                  lim = { zlog: RATIO_V_LOG, $
                                                          zrange: RATIO_V_RANGE, $
                                                          title:sort_title+'         '+PHASE+ $
                                                          '       O!U+!N BEAM RATIO at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+','$
                                                          +slice_block(1)+' ) ' $
                                                          +'in PS' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0] 
                                        pclose
                                    ENDIF 
                                    IF keyword_set(mogrify_ps) THEN BEGIN 
                                        spawn, 'mogrify -format png '+ path_wd +'*.ps'
                                        spawn, 'mogrify -rotate -90 '+ path_wd +'*.png'
                                    ENDIF 
                                ENDIF    
                            ENDFOR           
                        ENDIF 
                    ENDIF                                                                

;--------------------------------- PROPERTY_MAP ---------------------
                    IF KEYWORD_SET(PROPERTY_MAP_SET) THEN BEGIN 
;run for different properties as set in property_map_set
                        FOR ipp = 0, N_ELEMENTS(PROPERTY_MAP_SET)-1 DO  BEGIN 
; setup the plot saving directory
                            path_pp_main = path+property_map_set(ipp)+'/'
                            spawn, 'mkdir ' + path_pp_main
; input different data for different properties
                            property = property_map_set(ipp)
                            IF property EQ 'energy' THEN BEGIN 
                                property_tail = energy_tail
                                property_earth = energy_earth
                                PROPERTY_V_LOG = ENERGY_V_LOG
                                PROPERTY_V_RANGE = ENERGY_V_RANGE 
                            ENDIF 
                            IF property EQ 'energy_v' THEN BEGIN 
                                property_tail = energy_v_tail
                                property_earth = energy_v_earth
                                PROPERTY_V_LOG = velocity_V_LOG
                                PROPERTY_V_RANGE = velocity_V_RANGE 
                            ENDIF 

                            IF property EQ 'flux' THEN BEGIN 
                                property_tail = flux_tail
                                property_earth = flux_earth
                                PROPERTY_V_LOG = FLUX_V_LOG
                                PROPERTY_V_RANGE = FLUX_V_RANGE 
                            ENDIF 
                            IF property EQ 'density' THEN BEGIN 
                                property_tail = density_tail
                                property_earth = density_earth
                                PROPERTY_V_LOG = density_V_LOG
                                PROPERTY_V_RANGE = density_V_RANGE 
                            ENDIF 
                            IF property EQ 'velocity' THEN BEGIN 
                                property_tail = velocity_tail
                                property_earth = velocity_earth
                                PROPERTY_V_LOG = velocity_V_LOG
                                PROPERTY_V_RANGE = velocity_V_RANGE 
                            ENDIF 
                            IF property EQ 'pitch_angle' THEN BEGIN 
                                property_tail = pitch_angle_tail
                                property_earth = pitch_angle_earth
                                PROPERTY_V_LOG = pitch_angle_V_LOG
                                PROPERTY_V_RANGE = pitch_angle_V_RANGE 
                            ENDIF 

                            IF property EQ 'nVpara_over_B' THEN BEGIN 
                                property_tail = (density_tail*ABS(v_para_tail))/B
                                property_earth = (density_earth*ABS(v_para_earth))/B
                                PROPERTY_V_LOG = nvpara_over_b_v_log
                                PROPERTY_V_RANGE = nvpara_over_b_v_range
                            ENDIF 
                            IF property EQ 'temperature' THEN BEGIN 
                                property_tail = temperature_tail
                                property_earth = temperature_earth
                                PROPERTY_V_LOG = temperature_v_log
                                PROPERTY_V_RANGE = temperature_v_range
                            ENDIF
                            IF property EQ 'anodes' THEN BEGIN 
                                property_tail = anodes_tail
                                property_earth = anodes_earth
                                PROPERTY_V_LOG = anodes_v_log
                                PROPERTY_V_RANGE = anodes_v_range
                            ENDIF 
                            
                            IF NOT keyword_set(property_tail) THEN BEGIN 
                                print, 'no ' +property+' avaliabe for mapping'
                                stop
                                GOTO, next_property
                            ENDIF 
;----- Calculation -------
                            nx = CEIL(ABS(x_range(1) - x_range(0))/grid)
                            ny = CEIL(ABS(y_range(1) - y_range(0))/grid)
                            nz = CEIL(ABS(z_range(1)-z_range(0))/SLICE_GRID) 
                            x_axis = INDGEN(nx)*grid+(x_range(0) < x_range(1))+grid*0.5
                            y_axis = INDGEN(ny)*grid+(y_range(0) < y_range(1))+grid*0.5
                            z_axis = INDGEN(nz)*slice_grid+(z_range(0) < z_range(1))$
                                     +slice_grid*0.5
                            
                            property_total = REPLICATE(!VALUES.F_NAN, ntime, nx, ny, nz)
                            event_counts = REPLICATE(!VALUES.F_NAN, ntime, nx, ny, nz)
                            lobe_property = REPLICATE(!VALUES.F_NAN, ntime, nx, ny, nz)
                            lobe_counts =  REPLICATE(!VALUES.F_NAN, ntime, nx, ny, nz) 
                            bl_property = REPLICATE(!VALUES.F_NAN, ntime, nx, ny, nz)
                            bl_counts =  REPLICATE(!VALUES.F_NAN, ntime, nx, ny, nz)
                            ps_property = REPLICATE(!VALUES.F_NAN, ntime, nx, ny, nz)
                            ps_counts = REPLICATE(!VALUES.F_NAN, ntime, nx, ny, nz)       
                            FOR i = 0l, ntime-1 DO BEGIN 
                                index = SORT(ABS(x_axis - data_pos(i, 0))) 
                                index_x = index(0)
                                index = SORT(ABS(y_axis - data_pos(i, 1)))
                                index_y = index(0)
                                index = SORT(ABS(z_axis - data_pos(i, 2)))
                                index_z = index(0)
                                ;all
                                ;change abs(flag(i)) ge 1 to abs(flag(i))>0 so bidirection
                                ;beam will be recognized as two beams
                                property_total(i, index_x, index_y, index_z) = ((property_tail(i) > 0)+ (property_earth(i) > 0))* (ABS(flag(i)) GE 1)
                                ;lobe
                                lobe_property(i, index_x, index_y, index_z) = ((property_tail(i) > 0) +(property_earth(i) > 0)) * $
                                  (beta(i) LE 0.05 AND ABS(flag(i))GE 1)
                                ; boundary layer
                                bl_property(i, index_x, index_y, index_z) = ((property_tail(i) > 0)+(property_earth(i) > 0)) * $
                                  (beta(i) GT 0.05 AND beta(i) LE 1 AND ABS(flag(i))GE 1)      
                                ; plasma sheet
                                ps_property(I, index_x, index_y, index_z) = ((property_tail(i) > 0)+(property_earth(i) > 0))* $
                                  (beta(i) GT 1 AND  ABS(flag(i))GE 1) 
                            ENDFOR 

                                IF keyword_set(medium_value) THEN BEGIN 
                                    

                                    property_total(index_x, index_y, index_z)

                                ENDIF 
                                IF keyword_set(mean_value)THEN BEGIN 
                                ;all 
                                property_total(index_x, index_y, index_z) =  $         
                                  property_total(index_x, index_y, index_z) $
                                  +((property_tail(i) > 0)+ (property_earth(i) > 0))*(ABS(flag(i)) GE 1)
                                event_counts(index_x, index_y, index_z) =  $
                                  event_counts(index_x, index_y, index_z ) $
                                  + ABS(flag(i) > 0) 

                                ;change abs(flag(i)) ge 1 to abs(flag(i))>0 so bidirection
                                ;beam will be recognized as two beams

                                ;lobe
                                lobe_property(index_x, index_y, index_z) = $
                                  lobe_property(index_x, index_y, index_z)+ $
                                  ((property_tail(i) > 0) +(property_earth(i) > 0)) * $
                                  (beta(i) LE 0.05 AND ABS(flag(i))GE 1) 
                                lobe_counts(index_x, index_y, index_z) = $
                                  lobe_counts(index_x, index_y, index_z)+ $
                                  (beta(i) LE 0.05 AND (ABS(flag(i))) > 0)                 
                                ; boundary layer
                                bl_property(index_x, index_y, index_z) = $
                                  bl_property(index_x, index_y, index_z)+ $
                                  ((property_tail(i) > 0)+(property_earth(i) > 0))* $
                                  (beta(i) GT 0.05 AND beta(i) LE 1 $
                                   AND ABS(flag(i))GE 1)   
                                bl_counts(index_x, index_y, index_z) = $
                                  bl_counts(index_x, index_y, index_z)+ $
                                  (beta(i) GT 0.05 AND beta(i) LE 1 AND (ABS(flag(i)) > 0))     
                                ; plasma sheet
                                ps_property(index_x, index_y, index_z) = $
                                  ps_property(index_x, index_y, index_z)+ $
                                  ((property_tail(i) > 0)+(property_earth(i) > 0))* $
                                  (beta(i) GT 1 AND  ABS(flag(i))GE 1) 
                                ps_counts(index_x, index_y, index_z) = $
                                  ps_counts(index_x, index_y, index_z) + $
                                  (beta(i) GT 1 AND (ABS(flag(i)) > 0)) 
                            ENDIF 
                            ENDFOR  

;-------- 2d plot --------
                            IF KEYWORD_SET(PLOT_2D) THEN BEGIN 
;claculation
                                IF keyword_set(mean_value) THEN BEGIN 
                                    property_map_2d = TOTAL(property_total, 3, /nan)/ $
                                                      TOTAL(event_counts, 3, /nan)
                                    property_lobe_map_2d = TOTAL(lobe_property, 3, /nan)/ $
                                      TOTAL(lobe_counts, 3, /nan)
                                    property_bl_map_2d = TOTAL(bl_property, 3, /nan)/TOTAL(bl_counts, 3, /nan)
                                    property_ps_map_2d = TOTAL(ps_property, 3, /nan)/TOTAL(ps_counts, 3, /nan)
                                ENDIF 
;idl      
                                IF KEYWORD_SET(IDL_PLOT) THEN BEGIN 
                                    window, /free ;all
                                    specplot, x_axis, y_axis, property_map_2d, $
                                              no_interp = 1, $
                                              lim = { zlog: PROPERTY_V_LOG, $
                                                      zrange:PROPERTY_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N  beam   ' $
                                                      + property $
                                                      + 'distribution !CFROM  ' $
                                                      + ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]

; try to overplot the circle of magnetosphere in order to define solawind
;                                    x = indgen(91)*0.1-1
                                ;                                   x = [x, !values.f_nan]
                                ;                                  y = sqrt((1-x^2/8^2)*(15^2))+1
                                ;                                 x = [x, x]
                                ;                                y = [-y, y]
                                ;                               oplot, x, y
                                    IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                        window, /free ;lobe
                                        specplot, x_axis, y_axis, property_lobe_map_2d, $
                                                  no_interp = 1, $
                                                  lim = { zlog:PROPERTY_V_LOG, $
                                                          zrange:PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   ' $
                                                          +PROPERTY $
                                                          +' distribution IN LOBE!C' $
                                                          +'FROM  ' + ts_date $
                                                          +'  TO  ' +te_date,  $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        window, /free ; boundary layer
                                        specplot, x_axis, y_axis, property_bl_map_2d, $
                                                  no_interp = 1, $
                                                  lim = { zlog:PROPERTY_V_LOG, $
                                                          zrange:PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   ' $
                                                          +PROPERTY $
                                                          +' distribution IN BL!CFROM  ' $
                                                          + ts_date+'  TO  ' +te_date,  $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        window, /free ;plasma sheet
                                        specplot, x_axis, y_axis, property_ps_map_2d, $
                                                  no_interp = 1, $
                                                  lim = { zlog:PROPERTY_V_LOG, $
                                                          zrange:PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   '$
                                                          +PROPERTY $
                                                          +' distribution IN PS!CFROM  ' $
                                                          + ts_date+'  TO  ' +te_date,  $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    ENDIF              
                                ENDIF 
; ps          
                                IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                    path_pp_2d = path_pp_main+'2d/'
                                    spawn, 'mkdir ' + path_pp_2d 
                                    
                                    popen,  path_pp_2d+property+'_all_'+ $ ; all 
                                            ts_date+'_to_' + te_date+'_' $
                                            +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                            +'.ps', /land 
                                    specplot, x_axis, y_axis, property_map_2d, $
                                              no_interp = 1, $
                                              lim = { zlog:PROPERTY_V_LOG, $
                                                      zrange:PROPERTY_V_RANGE, $
                                                      title: sort_title+'         '+PHASE+$
                                                      '       O!U+!N  beam   ' + PROPERTY $
                                                      +' distribution!CFROM  ' $
                                                      + ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), $
                                                      ytitle: PLOT_AXIS(1), $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.9, 0.9]}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    pclose
                                    IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                        popen,  path_pp_2d+property+'_lobe_'+ $ ; lobe       
                                                ts_date+'_to_' + te_date+'_' $
                                                +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                +'.ps', /land 
                                        specplot, x_axis, y_axis, property_lobe_map_2d, $
                                                  no_interp = 1, $
                                                  lim = { zlog:PROPERTY_V_LOG, $
                                                          zrange: PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   '+PROPERTY $
                                                          +' distribution IN LOBE!CFROM  ' $
                                                          + ts_date+'  TO  ' +te_date,  $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        popen,  path_pp_2d+property+'_bl_'+ $ ;boundary layer
                                                ts_date+'_to_' + te_date+'_' $
                                                +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                +'.ps', /land 
                                        specplot, x_axis, y_axis, property_bl_map_2d, $
                                                  no_interp = 1, $
                                                  lim = { zlog:PROPERTY_V_LOG, $
                                                          zrange:PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   '+PROPERTY $
                                                          +' distribution IN BL!CFROM  ' $
                                                          + ts_date+'  TO  ' +te_date,  $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        popen,  path_pp_2d+property+'_ps_'+ $ ;plasma sheet
                                                ts_date+'_to_' + te_date+'_' $
                                                +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                +'.ps', /land 
                                        specplot, x_axis, y_axis, property_ps_map_2d, $
                                                  no_interp = 1, $
                                                  lim = { zlog: PROPERTY_V_LOG, $
                                                          zrange: PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   '+PROPERTY $
                                                          +' distribution IN PS!CFROM  ' $
                                                          + ts_date+'  TO  ' +te_date,  $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                    ENDIF 
                                    IF keyword_set(mogrify_ps) THEN BEGIN 
                                        spawn, 'mogrify -format png '+ path_pp_2d +'*.ps'
                                        spawn, 'mogrify -rotate -90 '+ path_pp_2d +'*.png'
                                    ENDIF                             
                                ENDIF  
                            ENDIF 
; --------------- slice plots --------------------
; z is the direction been sliced
                            path_pp_slice_main = path_pp_main+'slice/'
                            spawn, 'mkdir ' + path_pp_slice_main
                            IF KEYWORD_SET(SLICE_PLOT) THEN BEGIN 
                                FOR iz = 0, nz-1 DO BEGIN 
                                    slice_block = [STRING(z_axis(iz)-SLICE_GRID*0.5, $
                                                          format = '(f5.1)'), $
                                                   STRING(z_axis(iz)+SLICE_GRID*0.5, $
                                                          format = '(f5.1)')]
                                    
                                    slice_property = property_total(*, *, iz) / $
                                                     event_counts(*, *, iz)
                                    slice_lobe_property = lobe_property(*, *, iz) / $
                                      lobe_counts(*, *, iz)
                                    slice_bl_property = bl_property(*, *, iz) / $
                                                        bl_counts(*, *, iz)
                                    slice_ps_property = ps_property(*, *, iz) / $
                                                        ps_counts(*, *, iz)
;idl
                                    IF KEYWORD_SET(idl_plot) THEN BEGIN 
                                        window, /free ;all
                                        specplot, x_axis, y_axis,  slice_property, $
                                                  no_interp = 1, $
                                                  lim = { zlog: PROPERTY_V_LOG, $
                                                          zrange: PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   '+PROPERTY $
                                                          +' distribution at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+',' $
                                                          +slice_block(1)+' ) ' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, $
                                                          charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                            window, /free ;lobe
                                            specplot, x_axis, y_axis,  slice_lobe_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange:PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              +':( ' +slice_block(0)+',' $
                                                              +slice_block(1)+' ) ' $
                                                              + 'in LOBE'$
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                            window, /free ;boundary layer
                                            specplot, x_axis, y_axis,  slice_bl_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog: PROPERTY_V_LOG, $
                                                              zrange:PROPERTY_V_RANGE,  $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0)+',' $
                                                              +slice_block(1)+' ) ' $
                                                              + 'in BL'$
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                            window, /free ;plasma sheet
                                            specplot, x_axis, y_axis,  slice_ps_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange:PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0)+',' $
                                                              +slice_block(1)+' ) ' $
                                                              + 'in PS'$
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        ENDIF 
                                    ENDIF  
;ps
                                    IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                        path_pp_slice = path_pp_slice_main $
                                                        +'slice_zgrid_'+slice_grid_str+'/'
                                        spawn, 'mkdir '+ path_pp_slice
                                        
                                        popen, path_pp_slice+'all_'+property+'_' $ ;all
                                               + ts_date+'_to_' + te_date+'_' $
                                               + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $ 
                                               +'_at_' + plot_axis(2)+'_' $
                                               +slice_block(0)+'_'+slice_block(1) $
                                               +'.ps', /land 
                                        specplot, x_axis, y_axis,  slice_property, $
                                                  no_interp = 1, $
                                                  lim = { zlog: PROPERTY_V_LOG, $
                                                          zrange:PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   '+PROPERTY $
                                                          +' distribution at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+',' $
                                                          +slice_block(1)+' ) ' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, $
                                                          charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                            popen,  path_pp_slice+'lobe_'+property+'_' $ ;lobe
                                                    +ts_date+'_to_' + te_date $
                                                    +'_' +PLOT_AXIS(0)+'_vs_' +PLOT_AXIS(1) $
                                                    +'_at_' + plot_axis(2) $
                                                    +'_' +slice_block(0)+'_'+slice_block(1) $
                                                    +'.ps', /land 
                                            specplot, x_axis, y_axis,  slice_lobe_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange:PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0)+',' $
                                                              +slice_block(1)+' ) ' $
                                                              + 'in LOBE'$
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                            pclose
                                            popen,  path_pp_slice+'bl_' $ ;boundary layer
                                                    +property+'_' $ 
                                                    +ts_date+'_to_' + te_date $
                                                    +'_'+PLOT_AXIS(0) + '_vs_'+PLOT_AXIS(1) $
                                                    +'_at_' + plot_axis(2) $
                                                    +'_' +slice_block(0)+'_'+slice_block(1) $
                                                    +'.ps', /land 
                                            specplot, x_axis, y_axis,  slice_bl_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog: PROPERTY_V_LOG, $
                                                              zrange: PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0)+',' $
                                                              +slice_block(1)+' ) ' $
                                                              + 'in BL'$
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                            pclose                                                
                                            popen,  path_pp_slice+'ps_'+property+'_'$ ;plasma sheet
                                                    +ts_date+'_to_' + te_date+'_' $
                                                    +PLOT_AXIS(0) +'_vs_'+PLOT_AXIS(1) $
                                                    +'_at_' + plot_axis(2) $
                                                    +'_'+slice_block(0)+'_'+slice_block(1) $
                                                    +'.ps', /land 
                                            specplot, x_axis, y_axis,  slice_ps_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange:PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0) $
                                                              +','+slice_block(1)+' ) ' $
                                                              + 'in PS'$
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                            pclose
                                        ENDIF 
                                        IF keyword_set(mogrify_ps) THEN BEGIN 
                                            spawn, 'mogrify -format png '+ path_pp_slice +'*.ps'
                                            spawn, 'mogrify -rotate -90 '+ path_pp_slice +'*.png'
                                        ENDIF 
                                    ENDIF   
                                ENDFOR      
                            ENDIF         
; ------ WATER_DROP PLOT ------
;SLICE_GRID WILL BE ALSO USED HERE
                            IF KEYWORD_SET(WATERDROP_PLOT) THEN BEGIN 
                                path_pp_wd_main = path_pp_main +'waterdrop/'
                                spawn, 'mkdir ' + path_pp_wd_main
;run for all slice
                                FOR iz = 0, nz-1 DO BEGIN 
                                    IF iz LT CEIL(nz/2.) THEN  $
                                      index = INDGEN(CEIL(nz/2.)-iz)+iz  $
                                    ELSE index = INDGEN(iz+1-FLOOR(nz/2.))+FLOOR(nz/2.)   

                                    IF N_ELEMENTS(index) GT 1 THEN BEGIN 
                                        wd_property = TOTAL(property_total(*, *, index), 3, /nan)/ $
                                                      TOTAL(event_counts(*, *, index), 3, /nan)
                                        wd_lobe_property = $
                                          TOTAL(lobe_property(*, *, index), 3, /nan)/ $
                                          TOTAL(lobe_counts(*, *, index), 3, /nan)
                                        wd_bl_property = TOTAL(bl_property(*, *, index), 3, /nan)/$
                                                         TOTAL(bl_counts(*, *, index), 3, /nan)
                                        wd_ps_property = TOTAL(ps_property(*, *, index), 3, /nan)/$
                                                         TOTAL(ps_counts(*, *, index), 3, /nan)
                                    ENDIF ELSE BEGIN
                                        wd_property = property_total(*, *, index)/$
                                                      event_counts(*, *, index)
                                        wd_lobe_property = lobe_property(*, *, index)/$
                                                           (lobe_counts(*, *, index))
                                        wd_bl_property = bl_property(*, *, index)/$
                                                         bl_counts(*, *, index)
                                        wd_ps_property = ps_property(*, *, index)/$
                                                         ps_counts(*, *, index)
                                    ENDELSE   
                                    slice_block = [ STRING(z_axis(index(0))-SLICE_GRID*0.5, $
                                                           format = '(f5.1)'), $
                                                    STRING(z_axis(index(N_ELEMENTS(index)-1))$
                                                           + SLICE_GRID*0.5, $
                                                           format = '(f5.1)')]
;idl
                                    IF KEYWORD_SET(idl_plot) THEN BEGIN 
                                        window, /free ;all
                                        specplot, x_axis, y_axis, wd_property, $
                                                  no_interp = 1, $
                                                  lim = { zlog:PROPERTY_V_LOG, $
                                                          zrange: PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   '+PROPERTY $
                                                          +' distribution at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0)+',' $
                                                          +slice_block(1)+' ) ' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, $
                                                          charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                            window, /free ;lobe
                                            specplot, x_axis, y_axis, wd_lobe_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange: PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0) $
                                                              +','+slice_block(1)+' ) ' $
                                                              +'in LOBE' $
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]   
                                            window, /free ;boundary layer
                                            specplot, x_axis, y_axis, wd_bl_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange: PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0) $
                                                              +','+slice_block(1)+' ) ' $
                                                              +'in BL' $
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                            window, /free ;plasma sheet
                                            specplot, x_axis, y_axis, wd_ps_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange: PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0) $
                                                              +','+slice_block(1)+' ) ' $
                                                              +'in PS' $
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]   
                                        ENDIF 
                                    ENDIF 
; ps
                                    IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                        path_pp_wd = path_pp_wd_main $
                                                     +'wd_zgrid_' + slice_grid_str +'/'
                                        spawn, 'mkdir ' + path_pp_wd 
                                        
                                        popen,  path_pp_wd+'all_'+property+'_'+ $ ; all 
                                                ts_date+'_to_' + te_date+'_' $
                                                +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                +'_at_'+ plot_axis(2) $
                                                +'_'+slice_block(0)+'_'+slice_block(1) $
                                                +'.ps', /land 
                                        specplot, x_axis, y_axis, wd_property, $
                                                  no_interp = 1, $
                                                  lim = { zlog:PROPERTY_V_LOG, $
                                                          zrange: PROPERTY_V_RANGE, $
                                                          title: sort_title+'         '+PHASE+$
                                                          '       O!U+!N  beam   '+PROPERTY $
                                                          + ' distribution at  ' $
                                                          + plot_axis(2) $
                                                          + ':( ' +slice_block(0) $
                                                          +','+slice_block(1)+' ) ' $
                                                          + '!CFROM  ' + ts_date $
                                                          +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, $
                                                          charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.9, 0.9]}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        pclose
                                        IF KEYWORD_SET(diff_beta) THEN BEGIN 
                                            popen,   path_pp_wd+'lobe_'+property+'_'+ $ ;lobe
                                                     ts_date+'_to_' + te_date+'_' $
                                                     +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                     +'_at_'+ plot_axis(2) $
                                                     +'_'+slice_block(0)+'_'+slice_block(1) $
                                                     +'.ps', /land  
                                            specplot, x_axis, y_axis, wd_lobe_property,  $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange: PROPERTY_V_RANGE, $
                                                              title:sort_title+'         '+PHASE+ $
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0) $
                                                              +','+slice_block(1)+' ) ' $
                                                              +'in LOBE' $
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                            pclose                                 
                                            popen,  path_pp_wd+'bl_'+property+'_'+ $ ;bl
                                                    ts_date+'_to_' + te_date+'_' $
                                                    +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                    +'_at_'+ plot_axis(2)+'_'  $
                                                    +slice_block(0)+'_'+slice_block(1) $
                                                    +'.ps', /land  
                                            specplot, x_axis, y_axis, wd_bl_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange: PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0) $
                                                              +','+slice_block(1)+' ) ' $
                                                              +'in BL' $
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                            pclose
                                            popen,  path_pp_wd+'ps_'+property+'_'+ $ ;plasma sheet
                                                    ts_date+'_to_' + te_date+'_' $
                                                    +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                    +'_at_'+ plot_axis(2) $
                                                    +'_'+slice_block(0)+'_'+slice_block(1) $
                                                    +'.ps', /land  
                                            specplot, x_axis, y_axis, wd_ps_property, $
                                                      no_interp = 1, $
                                                      lim = { zlog:PROPERTY_V_LOG, $
                                                              zrange: PROPERTY_V_RANGE, $
                                                              title: sort_title+'         '+PHASE+$
                                                              '       O!U+!N  beam   '+PROPERTY $
                                                              +' distribution at  ' $
                                                              + plot_axis(2) $
                                                              + ':( ' +slice_block(0) $
                                                              +','+slice_block(1)+' ) ' $
                                                              +'in PS' $
                                                              + '!CFROM  ' + ts_date $
                                                              +'  TO  ' +te_date, $
                                                              xtitle: PLOT_AXIS(0), $
                                                              ytitle: PLOT_AXIS(1), $
                                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                                              XSTYLE:1, ystyle: 1, $
                                                              charsize: 1.2, $
                                                              position: [0.1, 0.1, 0.9, 0.9]}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                            pclose
                                        ENDIF 
                                        IF keyword_set(mogrify_ps) THEN BEGIN 
                                            spawn, 'mogrify -format png '+ path_pp_wd +'*.ps'
                                            spawn, 'mogrify -rotate -90 '+ path_pp_wd +'*.png'
                                        ENDIF 
                                    ENDIF     
                                ENDFOR               
                            ENDIF
                            next_property:    
                        ENDFOR        
                    ENDIF     
                ENDFOR   
            ENDFOR  
        ENDFOR         
    ENDFOR      
ENDFOR               
;stop
END
