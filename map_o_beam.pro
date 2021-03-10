FUNCTION cal_property_map_value, property_value, property_map_type
IF size(property_value, /n_dim) NE 4 AND size(property_value, /n_dim) NE 5 THEN stop
IF size(property_value, /n_dim) EQ 5 THEN BEGIN 
    nx = n_elements(property_value(0, 0, *, 0, 0))
    ny = n_elements(property_value(0, 0, 0, *, 0))
ENDIF ELSE BEGIN 
    nx = n_elements(property_value(0, 0, *, 0))
    ny = n_elements(property_value(0, 0, 0, *))
ENDELSE
property_result = REPLICATE(!VALUES.F_NAN, nx, ny)
;-------- mean
IF property_map_type EQ 'mean' THEN $
  IF size(property_value, /n_dim) EQ 5 $
    THEN property_result = TOTAL(TOTAL(TOTAL(property_value, 2, /NAN), 1, /NAN), 3, /NAN)/ $
      TOTAL(TOTAL(TOTAL(property_value GE 0, 2, /NAN), 1, /NAN), 3, /NAN) $
  ELSE  property_result = TOTAL(TOTAL(property_value, 2, /NAN), 1, /NAN)/ $
    TOTAL(TOTAL(property_value GE 0, 2, /NAN), 1, /NAN)
;-------- median
IF property_map_type EQ 'median' THEN $
  FOR ix = 0, nx-1 DO $
  FOR iy = 0, ny-1 DO BEGIN 
    IF size(property_value, /n_dim) EQ 5 THEN dummy = property_value(*, *, ix, iy, *) $
    ELSE  dummy = property_value(*, *, ix, iy)
    index = where(finite(dummy), ct)
    IF ct GT 0 THEN dummy = dummy(index)
    property_result(ix, iy) = median(dummy)
ENDFOR
;-------- peak ------
;IF property_map_type EQ 'peak' THEN $
 ; FOR ix = 0, nx-1 DO $
  ;FOR iy = 0, ny-1 DO BEGIN 
;    p_max = CEIL(MAX(property_value))
 ;   p_min = FLOOR(min(property_value)) < 0 ; min property gt 0 set as 0
  ;  IF size(property_value, /n_dim) EQ 5 THEN dummy = property_value(*, *, ix, iy, *) $
 ;   ELSE  dummy = property_value(*, *, ix, iy)
  ;  index = where(finite(dummy), ct)
   ; IF ct GT 0 THEN dummy = dummy(index)
    
 ;   np = 50
  ;  grid_p = (p_max-p_min)/50.
   ; p = indgen(np)*grid_p+p_min
 ;   sta = FLTARR(np)
  ;  FOR  i = 0l, np-1 DO sta(i) = Total(dummy GE (i*grid_x+x_min) AND dummy LT ((i+1)*grid_x+x_min)) 
    
 ;   index = where(sta EQ max(sta), ct)
  ;  property_result(ix, iy) = mean(p(index))
;ENDFOR
;-------- varition ------
IF property_map_type EQ 'varition' THEN $
  FOR ix = 0, nx-1 DO $
  FOR iy = 0, ny-1 DO BEGIN

  ;  p_max = CEIL(MAX(property_value))
   ; p_min = FLOOR(min(property_value)) < 0 ; min property gt 0 set as 0
;    IF size(property_value, /n_dim) EQ 5 THEN dummy = property_value(*, *, ix, iy, *) $
 ;   ELSE  dummy = property_value(*, *, ix, iy)
  ;  index = where(finite(dummy), ct)
   ; IF ct GT 0 THEN dummy = dummy(index)
        
;    np = 50
 ;   grid_p = (p_max-p_min)/50.
  ;  p = indgen(np)*grid_p+p_min
;    sta = FLTARR(np)
 ;   FOR  i = 0l, np-1 DO sta(i) = Total(dummy GE (i*grid_x+x_min) AND dummy LT ((i+1)*grid_x+x_min)) 
        
;    index = where(sta EQ max(sta), ct)
 ;   property_result(ix, iy) = mean(p(index))
ENDFOR 
return, property_result
END 

PRO map_o_beam, time, data, $
                sort_flag = sort_flag, sort_title = sort_title, $
                sc = sc, $
                ps_plot = ps_plot, $
                plot_path = plot_path, $
                coor_set = coor_set, $
                grid_set = grid_set, slice_grid_set = slice_grid_set, $
                diff_beta = diff_beta, $
                direction_set = direction_set, storm_phase_set = storm_phase_set, $
                plot_2d= plot_2d, slice_plot = slice_plot, waterdrop_plot = waterdrop_plot, $
                point_plot = point_plot, $
                events_map = events_map, $
                property_map_set = property_map_set, $
                symmetry_template = symmetry_template, $
                property_map_type_set = property_map_type_set,$
                make_table=make_table
;-- default keywords settings if keywords are not set-
IF NOT KEYWORD_SET(sc) THEN SC = 4
IF NOT KEYWORD_SET(coor_set) THEN $
  coor_set = [['X_GSE', 'Y_GSE'], ['X_GSE', 'Z_GSE'], ['X_GSM', 'Y_GSM'], ['X_GSM', 'Z_GSM']]
IF NOT KEYWORD_SET( grid_set)THEN   grid_set = [1., 2.]
IF NOT KEYWORD_SET(slice_grid_set)THEN  slice_grid_set = [20., 8.]  
IF NOT KEYWORD_SET( direction_set)THEN   direction_set =  ['all_directions', 'tailward', 'earthward', 'both']
IF NOT KEYWORD_SET( storm_phase_set)THEN  $
  storm_phase_set = ['prestorm', 'nonstorm_time', 'storm_time', 'initial_phase', 'main_phase', 'recovery'] 
IF NOT KEYWORD_SET( PROPERTY_MAP_SET) THEN $
  PROPERTY_MAP_SET = [''] ;['energy', 'flux','pitch_angle','density','temperature','velocity']
IF NOT KEYWORD_SET(plot_path) THEN BEGIN 
    plot_path = 'unknow_sort/'
    spawn, 'mkdir '+ plot_path
ENDIF
IF NOT KEYWORD_SET(sort_title) THEN sort_title = ''
property_map_type_set_all = ['mean', 'median', 'peak']
IF NOT KEYWORD_SET(property_map_type_set) THEN property_map_type_set = property_map_type_set_all

; 1: mean value; 2: median value 3: peak value
region_name_set = ['ALL', 'LOBE', 'BL', 'PS','le1','le05','le01','gt005le01','gt01','lt002']
; -- other basic settings --
sc_str = STRING(sc, format = '(i1.1)')

X_GSE_RANGE = [20., -20.] & Y_GSE_RANGE = [20., -20.] & Z_GSE_RANGE = [-20., 20.] 
X_GSM_RANGE = [20, -20.] & Y_GSM_RANGE = [-20., 20.] & Z_GSM_RANGE = [-20., 20.]
X_GSM_mod_RANGE = [20., -20.] & Y_GSM_mod_RANGE = [20., -20.] & Z_GSM_mod_RANGE = [-20., 20.]
MLT_RANGE = [0,24] & ILAT_RANGE =[60,90] 
EVENTS_V_LOG = 0 & EVENTS_V_RANGE = [1, 100.] & events_units = '# of events'
sample_units='# of samples'
RATIO_V_LOG = 0  & RATIO_V_RANGE = [0, 1.] & ratio_units='Occurance Frequency'
ENERGY_V_LOG = 1 & ENERGY_V_RANGE = [40, 400.] & energy_units = 'eV'
; for v_energy: lobe [40,400], psbl[200,2000],ps[500,5000]
FLUX_V_LOG = 1 & FLUX_V_RANGE = [10., 1000.] & flux_units = '1/cm!U-3!N-s-sr-(eV/e)'
EFLUX_V_LOG = 1 & EFLUX_V_RANGE = [1000.,10000.] & eflux_units = 'eV/cm!U-3!N-s-sr(eV/e)'
Density_v_log = 1 &  density_v_range = [0.0001, 0.01] & density_units = 'cm!U-3'
;0.001,0.1
velocity_v_log = 1 &  velocity_v_range = [20, 100] & velocity_units = 'km/s'
velocity_vpara_log = 1 &  velocity_vpara_range = [20, 100]
velocity_vperp_log = 1 &  velocity_vperp_range = [3, 30]
distfunc_v_log = 1 &  distfunc_v_range = [1e-10, 1e-6] & distfunc_units='(s!E3!N/cm!E3!N-km!E3!N)'
normal_distfunc_v_log=1 & normal_distfunc_v_range=[1e-10,1e-8] & normal_distfunc_units='(s!E3!N/cm!E3!N-km!E3!N)'
;1e-11,1e-9
nv_v_log = 1 & nv_v_range = [0.1, 10]*1e9 & nv_units = 'm!U-2!Ns!U-1'
pitch_angle_v_log = 0 & pitch_angle_v_range = [0, 90] & pitch_angle_units = 'deg'
nvpara_over_b_v_log = 1 &  nvpara_over_b_v_range = [0.001, 0.01] & nvpara_over_b_units = 'cm!U-3!N-km-s!U-1!N/gauss'
temperature_v_log = 1 &  temperature_v_range = [1000, 10000] & temperature_units = 'eV'
; or: lobe[1,100],psbl[10,1000], ps[100,10000]
anodes_v_log = 0 &  anodes_v_range = [1, 4] & anodes_units = 'anodes #'
By_v_log = 0 & By_v_range = [-50., 50.] & By_units='nT'
Beta_v_log=1 & Beta_v_range=[1e-4,1.] & Beta_units=''

la_x = 20 & la_y = -20          ; legend position
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
storm_phase = data(*, 42)    
energy_tail = data(*, 5)
energy_earth = data(*, 14)
mass_o = 16*1.6e-27*(1e3)^2/(1.6e-19) ; unit: ev/(km/s) 
energy_v_tail = sqrt(2*energy_tail/mass_o)
energy_v_earth = sqrt(2*energy_earth/mass_o)
flux_tail = data(*, 7)
flux_earth = data(*, 16)
eflux_tail = data(*, 48)
eflux_earth = data(*, 49)
density_tail = data(*, 8)
density_earth = data(*, 17)
velocity_tail = data(*, 9)
velocity_earth = data(*, 18)
v_energy_tail=mass_o*velocity_tail^2/2
v_energy_earth=mass_o*velocity_earth^2/2
pitch_angle_tail = data(*, 6)
pitch_angle_earth = data(*, 15)
mag = sqrt(data(*, 26)^2+data(*, 27)^2+data(*, 28)^2)
temperature_tail = data(*, 12)
temperature_earth = data(*, 21)
anodes_tail = (data(*, 50)+11.25)/22.5
anodes_earth = (data(*, 51)+11.25)/22.5
v_para_tail = data(*, 10)
v_perp_tail=data(*,11)
v_para_earth = data(*, 19)
v_perp_earth=data(*,20)
By=data(*,27)
distfunc_tail = data(*,59)
distfunc_earth= data(*,60)
ntime = N_ELEMENTS(time)
; reset with sort_flag
flag_o = flag_o*sort_flag
norm_factor_mlt=1/24.*360/180*!PI
; run for different storms

averaged_with_ratio=0 ; keyword setting for making table and the plots with it, if set then average from ratio in different grid rather than O+ counts/all counts
selected_region_name_x = ['polar','tail']
selected_region_x = [[10,-5],[-5,-20]]
; here z and y are opposite for the
; start settings are set in this way
selected_region_name_z = ['all'] ;['dawn','dusk','center']
selected_region_z = [-20,20]    ;[[-4,-12],[4,12],[-4,4]]
selected_region_name_y = ['all'] ;['south','north']

nnx=n_elements(selected_region_name_x)
nny=n_elements(selected_region_name_y)
nnz=n_elements(selected_region_name_z)
if keyword_set(averaged_with_ratio) then begin 
    selected_region_ratio_mean= FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
    selected_region_ratio_mean_error= FLTARR(n_elements(storm_phase_set),nnx,nny,nnz) 
    selected_region_ratio_median = FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
    selected_region_ratio_median_error = FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
endif else begin 
    a=FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
    b=FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
    selected_region_ratio = FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
    selected_region_ratio_error = FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
endelse 
FOR i_storm_phase = 0, N_ELEMENTS(storm_phase_set)-1 DO BEGIN 
;reset for different storm phases
    phase = storm_phase_set(i_storm_phase)
    flag_phase = flag_o
; combine flag and storm_phase 
    IF phase EQ 'initial_phase' THEN phase_flag = 1
    IF phase EQ 'main_phase' THEN phase_flag = 2
    IF phase EQ 'recovery' THEN phase_flag = 3
    IF phase EQ 'storm_time' THEN phase_flag = -1
    IF phase EQ 'nonstorm_time' THEN phase_flag = 0
    IF phase EQ 'all_time' THEN phase_flag = -2
    if phase eq 'prestorm' then phase_flag = 5
; the flags of data out of certain storm phase  are set to be 0  
    IF phase_flag GE 1 THEN flag_phase = flag_phase/(storm_phase EQ phase_flag) $
    ELSE BEGIN 
        if phase_flag eq 0 then flag_phase = flag_phase/(storm_phase eq 0 or storm_phase eq 5)
        IF phase_flag EQ -1 THEN flag_phase = flag_phase/(storm_phase ge 1 and storm_phase le 3)
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
        IF ARRAY_EQUAL(PLOT_AXIS, ['MLT', 'ILAT'])  THEN plot_axis_3 = 'Z_GSM'
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
        IF  ARRAY_EQUAL(PLOT_AXIS, ['X_GSM_mod', 'Y_GSM_mod']) OR $
          ARRAY_EQUAL(PLOT_AXIS, ['Y_GSM_mod', 'X_GSM_mod']) THEN  plot_axis_3 = 'Z_GSM_mod'
        IF  ARRAY_EQUAL(PLOT_AXIS, ['X_GSM_mod', 'Z_GSM_mod']) OR $
          ARRAY_EQUAL(PLOT_AXIS, ['Z_GSM_mod', 'X_GSM_mod']) THEN  plot_axis_3 = 'Y_GSM_mod'
        IF  ARRAY_EQUAL(PLOT_AXIS, ['Y_GSM_mod', 'Z_GSM_mod']) OR $ 
          ARRAY_EQUAL(PLOT_AXIS, ['Z_GSM_mod', 'Y_GSM_mod']) THEN  plot_axis_3 = 'X_GSM_mod'   
        
        IF NOT KEYWORD_SET(plot_axis_3) THEN BEGIN 
            PRINT, 'AXES CAN ONLY BE X,Y,Z IN GSE OR GSM COORDINATE SYSTEM OR MLT/ILAT'+ $
              '!C!CCONTINUE WITH DEFAUL AXES!C!CX_GSM VS Y_GSM !C!C INPUT ".CONT"'
            STOP
            PLOT_AXIS = ['X_GSM', 'Y_GSM']
            plot_axis_3 = 'Z_GSM'      
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
            IF PLOT_AXIS(i) EQ 'MLT' THEN BEGIN 
                data_pos(*, i) = data(*, 52)
                range(*, i) = MLT_RANGE
            ENDIF 
            IF PLOT_AXIS(i) EQ 'ILAT' THEN BEGIN 
                data_pos(*, i) = data(*,53 )
                range(*, i) = ILAT_RANGE
            ENDIF 
            IF PLOT_AXIS(i) EQ 'X_GSM_mod' THEN BEGIN 
                r=sqrt(data(*,23)^2+data(*,24)^2)
                sita=atan(data(*,24)/data(*,23))
                
                index=where(sita le 0 and data(*,24) gt 0,ct)
                if ct gt 0 then sita(index)=sita(index)+!PI
                index=where(sita ge 0 and data(*,24) lt 0,ct)
                if ct gt 0 then sita(index)=sita(index)-!PI
                data_pos(*, i) = r*cos(sita+4*!PI/180)
                range(*, i) = X_GSM_mod_RANGE
            ENDIF 
            IF PLOT_AXIS(i) EQ 'Y_GSM_mod' THEN BEGIN 
                r=sqrt(data(*,23)^2+data(*,24)^2)
                sita=atan(data(*,24)/data(*,23))
                
                index=where(sita le 0 and data(*,24) gt 0,ct)
                if ct gt 0 then sita(index)=sita(index)+!PI
                index=where(sita ge 0 and data(*,24) lt 0,ct)
                if ct gt 0 then sita(index)=sita(index)-!PI
                data_pos(*, i) = r*sin(sita+4*!PI/180)
                range(*, i) = X_GSM_mod_RANGE
                range(*, i) = Y_GSM_mod_RANGE
            ENDIF 
            IF PLOT_AXIS(i) EQ 'Z_GSM_mod' THEN BEGIN 
                data_pos(*, i) = data(*, 25) 
                range(*, i) = Z_GSM_mod_RANGE
            ENDIF 
            PRINT, '    '+PLOT_AXIS(i)+'     '
        ENDFOR   
        X_RANGE = range(*, 0)
        Y_RANGE = range(*, 1)
        Z_RANGE = range(*, 2)
        r_range=ABS(ILAT_RANGE(0)-ILAT_RANGE(1))
;--------------------------------------------------------------------------------------------
;                                      Calculation and plot                                 ;
;--------------------------------------------------------------------------------------------
;For MLT/ILAT plot, since magnetic field used for MLT and ILAT is
;dipole field, only data at polar region
        flag_phase_axis=flag_phase
        IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
          THEN BEGIN
            index = where(sqrt(data(*, 2)^2+data(*,3)^2+data(*,4)^2) gt 5, ct) ; only keep data close to the Earth for this map
                                ;        index=where(  	data(*,23) le -2 and data(*,23) ge -5 $
                                ;                           and data(*,25) ge -2 and data(*,25) le 3 $
                                ;                          and data(*,24) ge -4 and data(*,24) le 4 ,ct)
                                ;       a= flag_phase
                                ;      a(*)=!VALUES.F_NAN
                                ;     a(index)=flag_phase(index)
                                ;    flag_phase=a

            IF ct GT 0 THEN flag_phase_axis(index) = !VALUES.F_NAN
            index=where(data(*,53) eq 0,ct)
            if ct gt 0 then flag_phase_axis(index)=!VALUES.F_NAN
        endif
; run for different directions as set in direction_set
        FOR idirection = 0, N_ELEMENTS(direction_set)-1 DO BEGIN 
            flag = flag_phase_axis
            direction = direction_set(idirection)
; the flags of data out of direction are set as no event            
            IF direction EQ 'tailward' THEN begin 
                index= where(flag EQ -1,ct) 
                if ct gt 0 then flag(index)= !VALUES.F_NAN
            endif 
            IF direction EQ 'earthward' THEN begin 
                index= where(flag EQ 1,ct) 
                if ct gt 0 then flag(index)= !VALUES.F_NAN
            endif 
            IF direction EQ 'both' THEN begin 
                index= where(flag eq  1 or flag eq -1,ct) 
                if ct gt 0 then flag(index)= !VALUES.F_NAN
            endif 
            
            spawn, 'mkdir '+ plot_path+direction+'/'            
; run for different grids as set in grid_set
            FOR i_grid = 0, N_ELEMENTS(grid_set)-1 DO BEGIN 
                grid = grid_set(i_grid) & GRID_STR = STRING(GRID, FORMAT = '(i1)')
                
                path = plot_path+direction+'/'+'grid_'+grid_str+'/' & spawn, 'mkdir '+ path
                path = path+storm_phase_set(i_storm_phase)+'/' &  spawn, 'mkdir '+ path
                path_ev = path+'events/' &  spawn, 'mkdir '+ path_ev
                spawn, 'mkdir '+ path_ev +'png/'
                IF KEYWORD_SET(slice_plot) OR KEYWORD_SET(waterdrop_plot) THEN $
                  nslice =  N_ELEMENTS(slice_grid_set) ELSE nslice = 1

;------------------------- POINT PLOT -------------------------
                IF KEYWORD_SET(point_plot) THEN BEGIN
;--- Calculation ---
;classify data into different catergary
                    index_valid=where(ABS(flag) ge 0,nvalid)
                    if nvalid gt 0 then begin 
                        point_map_data = REPLICATE(!VALUES.F_NAN, nvalid, 3, 3)
                        FOR i = 0l, nvalid-1 DO BEGIN   
                            point_map_data(i, *, 0) = data_pos(index_valid(i), *)*(beta(index_valid(i)) LE 0.05)
                            point_map_data(i, *, 1) = data_pos(index_valid(i), *)*(beta(index_valid(i)) LE 1)*(beta(index_valid(i)) GT 0.05)
                            point_map_data(i, *, 2) = data_pos(index_valid(i), *)*(beta(index_valid(i)) GT 1)
                        ENDFOR     
                        
;--- draw the plot --- 
                        IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
                          THEN BEGIN
; if axis input is ilat/mlt, change it back to mlt/ilat since the
; later part of point map code are run for mlt/ilat
                            if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin
                                temp=point_map_data(*,*,0)
                                point_map_data(*,*,0)=point_map_data(*,*,1)
                                point_map_data(*,*,1)=point_map_data(*,*,0)
                            endif 
                            
                            IF KEYWORD_SET(ps_plot) THEN popen,  $
                              path_ev+'events_' + ts_date+'_to_'+ te_date+'_' $
                              + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'.ps', /land 
                            PLOT, [0, 0, -100, 100], [-100, 100, 0, 0], $
                              title = sort_title+'         '+PHASE $
                              + '       O!U+!N  BEAM  EVENS    ' $
                              +'!Cfrom  ' +ts_date+'  to  ' +te_date,  $
                              xtitle = PLOT_AXIS(0), ytitle = PLOT_AXIS(1),  $
                              xrange = [-r_range,r_range], yrange = [-r_range,r_range], $
                              XSTYLE = 5, ystyle = 5, charsize = 1.5, $
                              position = [0.15, 0.15, 0.85, 0.85]

; plot data at different region: lobe, ps, bl                         
                            if keyword_set(ps_plot) then psym_point_map=1 else psym_point_map=4
                            FOR iregion = 0, 2 DO begin 
; plot data that has no beam
                                index=where(flag(index_valid) eq 0,ct)
                                if ct gt 0 then oplot, 90-point_map_data(index,1,iregion), $
                                  point_map_data(index,0,iregion)*norm_factor_mlt,color=5,psym=psym_point_map,/polar
                                index=where(ABS(flag(index_valid)) gt 0,ct)
                                if ct gt 0 then oplot, 90-point_map_data(index,1,iregion), $
                                  point_map_data(index,0,iregion)*norm_factor_mlt,psym = psym_point_map, color = 3-iregion,/polar
                            endfor 
; legend         
                            xyouts, la_x+7, la_y+4.5*grid, 'no events', color = 5
                            xyouts, la_x+7, la_y, 'lobe', color = 3
                            xyouts, la_x+7, la_y+1.5*grid, 'boundary layer', color = 2
                            xyouts, la_x+7, la_y+3*grid, 'plasma sheet', color = 1
                            xyouts, r_range, 0, '0'
                            xyouts, 0, r_range, '6'
                            xyouts, -r_range*1.05,0, '12'
                            xyouts, 0 ,-r_range*1.05, '18'
                            for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                            for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                            IF KEYWORD_SET(ps_plot) THEN pclose ELSE stop
                        ENDIF  ELSE BEGIN 
                            IF KEYWORD_SET(ps_plot) THEN popen,  $
                              path_ev+'events_' + ts_date+'_to_'+ te_date+'_' $
                              + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'.ps', /land 
                            PLOT, [0, 0, -100, 100], [-100, 100, 0, 0], $
                              title = sort_title+'         '+PHASE $
                              + '       O!U+!N  BEAM  EVENS    ' $
                              +'!Cfrom  ' +ts_date+'  to  ' +te_date,  $
                              xtitle = PLOT_AXIS(0), ytitle = PLOT_AXIS(1),  $
                              xrange = X_RANGE, yrange = Y_RANGE, $
                              charsize = 1.5, $
                              position = [0.15, 0.15, 0.85, 0.85]
; plot data with or withnot beam at different region: lobe, ps, bl   ;1
                            if keyword_set(ps_plot) then psym_point_map=1 else psym_point_map=3
                            FOR iregion = 0, 2 DO begin 
                                oplot, point_map_data(*, 0, iregion)*(flag(index_valid) EQ 0), $
                                  point_map_data(*, 1, iregion)*(flag(index_valid) EQ 0), $
                                  symsize=0.5,$
                                  psym = psym_point_map, color = 2;5

                            endfor 
                            for iregion= 0 , 2 do begin 
                                oplot, point_map_data(*, 0, iregion)*(ABS(flag(index_valid)) GE 1), $
                                  point_map_data(*, 1, iregion)*(ABS(flag(index_valid)) GE 1), $
                                  symsize=0.5,$
                                  psym = psym_point_map, color = 2;3-iregion
                            endfor 
; legend         
                                ;              xyouts, la_x, la_y+4.5*grid, 'no events', color = 5
                                ;             xyouts, la_x, la_y, 'lobe', color = 3
                                ;            xyouts, la_x, la_y+1.5*grid, 'boundary layer', color = 2
                                ;           xyouts, la_x, la_y+3*grid, 'plasma sheet', color = 1
                            
                            OPLOT, [0, 0, -100, 100], [-100, 100, 0, 0], col = 0, thick = 8
                            IF KEYWORD_SET(ps_plot) THEN pclose ELSE stop 
                            spawn, 'mogrify -format png '+ path_ev +'*.ps'
                            spawn, 'mogrify -rotate -90 '+ path_ev +'*.png'
                            spawn, 'mv -f '+ path_ev +'*.png '+path_ev+'png/'
                            spawn, 'gzip -9f ' + path_ev+'*.ps'
                        ENDELSE  
                    endif 
                ENDIF           
; run for different slice cut
                FOR i_slice_grid = 0, nslice-1 DO BEGIN 
                    slice_grid = slice_grid_set(i_slice_grid)
                    slice_grid_str = STRING(slice_grid, format = '(i2.2)')        
;--------------------------------- EVENTS MAP ---------------------
                    IF KEYWORD_SET(EVENTS_MAP) THEN BEGIN 
                        path_ev = path+'events/' &  spawn, 'mkdir '+ path_ev
;-----  Calculation ------  
                        IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                          or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) THEN BEGIN
                            if ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT'])  then begin 
                                grid_x= 1*grid & grid_y=2.5*grid
                            endif else begin 
                                grid_x=1*grid & grid_y=2.5*grid
                            endelse 
                        endif  else begin 
                            grid_x=grid & grid_y=grid
                        endelse 
                        nx = CEIL(ABS(x_range(1) - x_range(0))/grid_x)
                        ny = CEIL(ABS(y_range(1) - y_range(0))/grid_y)
                        nz = CEIL(ABS(z_range(1)-z_range(0))/SLICE_GRID) 
                        x_axis = INDGEN(nx)*grid_x+(x_range(0) < x_range(1))+grid_x*0.5
                        y_axis = INDGEN(ny)*grid_y+(y_range(0) < y_range(1))+grid_y*0.5
                        z_axis = INDGEN(nz)*slice_grid+(z_range(0) < z_range(1))+ slice_grid*0.5
                        
                        IF NOT KEYWORD_SET(diff_beta) THEN  diff_beta = [0]
                        ns = n_elements(diff_beta)
                        FOR iis = 0, ns-1 DO BEGIN
                            is = diff_beta(iis)
                            IF is EQ 0 THEN ind_beam = where(ABS(flag) GE 0, ct_beam)
                            IF is EQ 1 THEN ind_beam = where(ABS(flag) GE 0 AND beta LE 0.05, ct_beam)
                            IF is EQ 2 THEN ind_beam = where(ABS(flag) GE 0 AND beta LE 1 AND beta GT 0.05, ct_beam)
                            IF is EQ 3 THEN ind_beam = where(ABS(flag) GE 0 AND beta GT 1, ct_beam)

                            IF is EQ 4 THEN ind_beam = where(ABS(flag) GE 0 AND beta LE 1, ct_beam)
                            IF is EQ 5 THEN ind_beam = where(ABS(flag) GE 0 AND beta LE 0.5, ct_beam)
                            IF is EQ 6 THEN ind_beam = where(ABS(flag) GE 0 AND beta LE 0.1, ct_beam)

                            IF is EQ 7 THEN ind_beam = where(ABS(flag) GE 0 AND beta gt 0.05 and beta LE 0.1, ct_beam)
                            IF is EQ 8 THEN ind_beam = where(ABS(flag) GE 0 AND beta gt 0.1, ct_beam)
                            IF is EQ 9 THEN ind_beam = where(ABS(flag) GE 0 AND beta LE 0.02, ct_beam)
                            IF ct_beam eq 0 then begin 
                                OPENU, unit, path+'log_nodata.txt', /GET_LUN, /APPEND
                                PRINTF, unit, 'no data for'+ region_name_set(is)
                                FREE_LUN, unit     
                            endif else begin 
                                event_counts = fltarr(nx, ny, nz)
                                total_counts = fltarr(nx, ny, nz)
                                FOR i = 0l, ct_beam-1 DO BEGIN 
                                    if data_pos(ind_beam(i),0) ge min(x_range) and data_pos(ind_beam(i),0) le max (x_range) $
                                     and  data_pos(ind_beam(i),1) ge min(y_range) and data_pos(ind_beam(i),1) le max(y_range) $
                                      and  data_pos(ind_beam(i),2) ge min(z_range) and data_pos(ind_beam(i),2) le max(z_range) $
                                      then begin 
                                        index = SORT(ABS(x_axis - data_pos(ind_beam(i), 0))) &  index_x = index(0)
                                        index = SORT(ABS(y_axis - data_pos(ind_beam(i), 1))) &  index_y = index(0)
                                        index = SORT(ABS(z_axis - data_pos(ind_beam(i), 2))) &  index_z = index(0)
                                        
                                        total_counts(index_x, index_y, index_z) = total_counts(index_x, index_y, index_z) + 1
                                        event_counts(index_x, index_y, index_z) = event_counts(index_x, index_y, index_z )+ $
                                          (ABS(flag(ind_beam(i))) GE 1) 
                                    endif 
                                ENDFOR 
                     ;           index = where(total_counts le 4,ct)
                      ;          if ct gt 0 then begin 
                       ;             total_counts(index) = 0
                        ;            event_counts(index) = 0
                         ;       endif 
                                event_ratio = event_counts/total_counts                                
;------ make table ----
                                if KEYWORD_SET(make_table) then begin
                                    path_table = path+'table/' &  spawn, 'mkdir '+path_table
                                    spawn, 'mkdir '+ path_table +'png/'
                                    
                                    for i_selected_region_x = 0, n_elements(selected_region_name_x)-1 do begin
                                        for i_selected_region_y = 0, n_elements(selected_region_name_y)-1 do begin
                                            if selected_region_name_x(i_selected_region_x) eq 'tail' then selected_region_y= [-20,20] 
                                            if selected_region_name_x(i_selected_region_x) eq 'polar' then  selected_region_y= [[-4,-20],[4,20]]
                                            for i_selected_region_z = 0,n_elements(selected_region_name_z)-1 do begin
                                                selected_region_name = selected_region_name_x(i_selected_region_x)+'_'+selected_region_name_y(i_selected_region_y)+'_'+selected_region_name_z(i_selected_region_z)
                                                index_x = where(x_axis lt min(selected_region_x(*,i_selected_region_x)) or x_axis gt max(selected_region_x(*,i_selected_region_x)),ct_x)
                                                if n_elements(selected_region_y) eq n_elements(selected_region_name_y)*2 then index_y = where(y_axis lt min(selected_region_y(*,i_selected_region_y)) or y_axis gt max(selected_region_y(*,i_selected_region_y)),ct_y)
                                                
                                                if n_elements(selected_region_y) eq n_elements(selected_region_name_y)*4 then index_y = where (y_axis lt 4 and y_axis gt -4 ,ct_y)

                                                index_z = where(z_axis lt min(selected_region_z(*,i_selected_region_z)) or z_axis gt max(selected_region_z(*,i_selected_region_z)),ct_z)
                                                if keyword_set(averaged_with_ratio) then begin 
                                                    selected_region_event_ratio = event_ratio
                                                    if ct_x gt 0 then selected_region_event_ratio(index_x,*,*) = !VALUES.F_NAN
                                                    if ct_y gt 0 then selected_region_event_ratio(*,index_y,*) = !VALUES.F_NAN
                                                    if ct_z gt 0 then selected_region_event_ratio(*,*,index_z) = !VALUES.F_NAN
                                                    
                                                    if keyword_set(ps_plot) then popen, path_table+'table_'+selected_region_name+'.ps', /land                                       
                                                    surface, DIST(5), /nodata, /save, $
                                                      xrange = [7, -20], yrange = [20,-20], zrange = [-20, 20], $
                                                      xstyle = 1, ystyle = 1, zstyle = 1, charsize = 1.2, $
                                                      position = [0.1, 0.1, 0.9, 0.9, 0.1, 0.9], $
                                                      xticklen = 1, yticklen = 1, zticklen = 1, $
                                                      xgridstyle = 1, ygridstyle = 1, zgridstyle = 1, $
                                                      az = 30, ax = 30, $333
                                                      xtitle = plot_axis(0)+' (R!DE!N)', $
                                                      ytitle = plot_axis(1)+' (R!DE!N)', $
                                                      ztitle =plot_axis(2)+' (R!DE!N)', $
                                                      title = '  '+phase+', '+selected_region_name+'!C!Cmean value: ' $
                                                      + string(100*mean(selected_region_event_ratio,/nan),format='(f6.2)')+'%, median value: ' $
                                                      + string(100*median(selected_region_event_ratio),format='(f6.2)')+'%'                                                
                                                    plots,[-20,7],[0,0],[0,0],color = -1,psym = -3,/t3d
                                                    plots,[0,0],[-20,20],[0,0],color = -1,psym = -3,/t3d
                                                    plots,[0,0],[0,0],[-15,15],color = -1,psym = -3,/t3d
                                                    for ix=0,n_elements(x_axis)-1 do for iy=0,n_elements(y_axis)-1 do for iz=0,n_elements(z_axis)-1 do  if ABS(selected_region_event_ratio(ix,iy,iz)) ge 0  and total_counts(ix,iy,iz) ge 0 then plots, (ix-(n_elements(x_axis)-1)/2.)*slice_grid,(iy-(n_elements(y_axis)-1)/2)*slice_grid,(iz-(n_elements(z_axis)-1)/2.)*slice_grid, psym = 1, color = selected_region_event_ratio(ix,iy,iz)*(254.-7.)+7., symsize = 3,/t3d
                                                    if keyword_set(ps_plot) then pclose else stop
                                                    
                                                    selected_region_ratio_mean(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                                                      =mean(selected_region_event_ratio,/nan)
                                                    selected_region_ratio_median(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                                                      =median(selected_region_event_ratio)
                                                    
                                                    OPENU, unit, plot_path+'table_occurrence_frequency.txt', /GET_LUN, /APPEND
                                                    PRINTF, unit, +phase+', '+selected_region_name $
                                                      +', mean value: '+string(100*mean(selected_region_event_ratio,/nan),format='(f6.2)') $
                                                      +'%, median value: '+string(100*median(selected_region_event_ratio),format='(f6.2)')+'%' 
                                                    FREE_LUN, unit
                                                endif else begin
                                                    selected_region_event_counts=event_counts
                                                    selected_region_total_counts=total_counts

                                                    if ct_x gt 0 then begin 
                                                        selected_region_event_counts(index_x,*,*) = !VALUES.F_NAN
                                                        selected_region_total_counts(index_x,*,*) = !values.f_nan
                                                    endif 
                                                    if ct_y gt 0 then begin 
                                                        selected_region_event_counts(*,index_y,*) = !VALUES.F_NAN
                                                        selected_region_total_counts(*,index_y,*) = !values.f_nan
                                                    endif
                                                    if ct_z gt 0 then begin 
                                                        selected_region_event_counts(*,*,index_z) = !VALUES.F_NAN
                                                        selected_region_total_counts(*,*,index_z) = !values.f_nan
                                                    endif 
                                                    
                                                    a(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                                                      = total(selected_region_event_counts,/nan)
                                                    b(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                                                      = total(selected_region_total_counts,/nan)

                                                    selected_region_ratio(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                                                      = total(selected_region_event_counts,/nan)/total(selected_region_total_counts,/nan)
                                                    
                                                    selected_region_ratio_error(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                                                      = selected_region_ratio(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                                                      *sqrt(1/total(selected_region_event_counts,/nan)+1/total(selected_region_total_counts,/nan))
                                                    
                                                    OPENU, unit, plot_path+'table_occurrence_frequency.txt', /GET_LUN, /APPEND
                                                    PRINTF, unit, phase+', '+selected_region_name+', ratio: '$
                                                      +string(100*selected_region_ratio(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z),format='(f6.2)')+'%'
                                                    FREE_LUN, unit
                                                endelse
                                            endfor 
                                        endfor 
                                    endfor 
                                    if keyword_set(ps_plot) then begin 
                                        spawn, 'mogrify -format png '+ path_table +'*.ps'
                                        spawn, 'mogrify -rotate -90 '+ path_table +'*.png'
                                        spawn, 'mv -f '+ path_table +'*.png '+ path_table +'png/'
                                        spawn, 'gzip -9f ' + path_table +'*.ps'        
                                    endif  
                                endif       
;------ 2D plot --------
                                IF KEYWORD_SET(PLOT_2D) THEN BEGIN     
                                    path_2d = path_ev+'2d/' &  spawn, 'mkdir '+path_2d
                                    spawn, 'mkdir '+ path_2d +'png/'

                                    total_counts_2d = TOTAL(total_counts, 3, /nan)
                                    index=where(total_counts_2d eq 0,ct)
                                    if ct gt 0 then  total_counts_2d(index) =  !VALUES.F_NAN 

                                    event_counts_2d =  TOTAL(event_counts, 3, /nan)
                                    index=where(event_counts_2d eq 0,ct)
                                    if ct gt 0 then  event_counts_2d(index) =  !VALUES.F_NAN 
                                    event_ratio_2d = TOTAL(event_counts, 3, /nan)/TOTAL(total_counts, 3, /nan)
                                ; use a symmetry template if required
                                    IF keyword_set(symmetry_template) THEN BEGIN
                                      
                                        positive_half = total_counts_2d(*, where(y_axis GT 0)) GT  0
                                        positive_half = [[reverse(positive_half, 2)], [positive_half]]
                                        negative_half = total_counts_2d(*, where(y_axis LT 0)) GT  0
                                        negative_half = [[negative_half], [reverse(negative_half, 2)]]
                                        template_2d = FLOAT(positive_half*negative_half)
                                        template_2d(where(template_2d EQ 0 )) =  !VALUES.F_NAN

                                        total_counts_2d=total_counts_2d*template_2d
                                        event_counts_2d = event_counts_2d*template_2d
                                        event_ratio_2d = event_ratio_2d*template_2d
                                    ENDIF
;samples
                                    IF KEYWORD_SET(ps_plot) THEN     popen,  path_2d+region_name_set(is) +'_samples_' $
                                      +ts_date+'_to_' + te_date+'_'  +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'.ps', /land 
                                    specplot, x_axis, y_axis,  total_counts_2d, $
                                      no_interp = 1, $
                                      lim = {  $
                                              zlog:EVENTS_V_LOG, zrange: EVENTS_V_RANGE, $
                                              title: sort_title+'  '+PHASE+$
                                              ' O!U+!N BEAM SAMPLES in '+region_name_set(is) +' Region!Cfrom ' $
                                              + ts_date+'  to  ' +te_date,  $
                                              xtitle: PLOT_AXIS(0), ytitle: PLOT_AXIS(1), $
                                              ztitle: sample_units, zticklen: -1,$
                                              xrange: x_range, yrange: y_range, $
                                              XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                              position: [0.1, 0.1, 0.85, 0.85], $
                                              zticks: 4}     
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                    IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT'])  or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) THEN BEGIN
                                        IF KEYWORD_SET(ps_plot) THEN     popen,  path_2d+region_name_set(is) +'_samples_' $
                                          +ts_date+'_to_' + te_date+'_'  +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'_polar.ps', /land 
                                        
                                        if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                            polar_spec,(90-x_axis),y_axis*norm_factor_mlt,event_counts_2d, $
                                              r_range=90-x_range,zrange=events_v_range,charsize=1.2,$
                                              ztitle=events_units,zlog=events_v_log,zticklen= 1,zticks=4,$
                                              title=sort_title+'  '+PHASE+$
                                              ' O!U+!N BEAM SAMPLESS in '+region_name_set(is) +' Region!Cfrom ' $
                                              + ts_date+'  to  ' +te_date
                                        endif  else begin 
                                            polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(total_counts_2d), $
                                              r_range=90-y_range,zrange=events_v_range, charsize=1.2,$
                                              ztitle= events_units,zlog=events_v_log,zticklen= -1,zticks=4,$
                                              title=sort_title+'  '+PHASE+$
                                              ' O!U+!N BEAM EVENTS in '+region_name_set(is) +' Region!Cfrom ' $
                                              + ts_date+'  to  ' +te_date
                                        endelse
                                        for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                        for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                                        IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                    endif                       
; events                                              
                                    IF KEYWORD_SET(ps_plot) THEN     popen,  path_2d+region_name_set(is) +'_events_' $
                                      +ts_date+'_to_' + te_date+'_'  +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'.ps', /land 
                                    specplot, x_axis, y_axis,  event_counts_2d, $
                                      no_interp = 1, $
                                      lim = {  $
                                              zlog:EVENTS_V_LOG, zrange: EVENTS_V_RANGE, $
                                              title: sort_title+'  '+PHASE+$
                                              ' O!U+!N BEAM EVENTS in '+region_name_set(is) +' Region!Cfrom ' $
                                              + ts_date+'  to  ' +te_date,  $
                                              xtitle: PLOT_AXIS(0), ytitle: PLOT_AXIS(1), $
                                              ztitle: events_units, zticklen: -1,$
                                              xrange: x_range, yrange: y_range, $
                                              XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                              position: [0.1, 0.1, 0.85, 0.85], $
                                              zticks: 4}     
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                    IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT'])  or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) THEN BEGIN
                                        IF KEYWORD_SET(ps_plot) THEN     popen,  path_2d+region_name_set(is) +'_events_' $
                                          +ts_date+'_to_' + te_date+'_'  +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'_polar.ps', /land 
                                        
                                        if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                            polar_spec,(90-x_axis),y_axis*norm_factor_mlt,event_counts_2d, $
                                              r_range=90-x_range,zrange=events_v_range,charsize=1.2,$
                                              ztitle=events_units,zlog=events_v_log,zticklen= 1,zticks=4,$
                                              title=sort_title+'  '+PHASE+$
                                              ' O!U+!N BEAM EVENTS in '+region_name_set(is) +' Region!Cfrom ' $
                                              + ts_date+'  to  ' +te_date
                                        endif  else begin 
                                            polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(event_counts_2d), $
                                              r_range=90-y_range,zrange=events_v_range, charsize=1.2,$
                                              ztitle= events_units,zlog=events_v_log,zticklen= -1,zticks=4,$
                                              title=sort_title+'  '+PHASE+$
                                              ' O!U+!N BEAM EVENTS in '+region_name_set(is) +' Region!Cfrom ' $
                                              + ts_date+'  to  ' +te_date
                                        endelse
                                        for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                        for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                                        IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                    endif  
; ratio  
                                    IF KEYWORD_SET(PS_PLOT) THEN popen, path_2d+region_name_set(is)+'_ratio_' $
                                      +ts_date+'_to_' + te_date+'_'+ PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1)  +'.ps', /land 
                                    specplot, x_axis, y_axis, event_ratio_2d, $
                                      no_interp = 1, $
                                      lim = { zlog: RATIO_V_LOG, $
                                              zrange: RATIO_V_RANGE, $
                                              title: sort_title+'  '+PHASE+ $
                                              ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cfrom ' $
                                              +ts_date+'  TO  ' +te_date,  $
                                              xtitle: PLOT_AXIS(0), ytitle: PLOT_AXIS(1), $
                                              ztitle: ratio_units, zticklen: -1,$
                                              xrange: X_RANGE, yrange: Y_RANGE, $
                                              XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                              position: [0.1, 0.1, 0.85, 0.85],$
                                              zticks:4}   
                                    oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                    IF KEYWORD_SET(ps_plot) THEN    pclose   ELSE stop
                                    IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT'])  or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) THEN BEGIN
                                        IF KEYWORD_SET(ps_plot) THEN     popen,  path_2d+region_name_set(is) +'_ratio_' $
                                          +ts_date+'_to_' + te_date+'_'  +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'_polar.ps', /land 
                                        
                                        if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                            polar_spec,(90-x_axis),y_axis*norm_factor_mlt,ratio_counts_2d, $
                                              r_range= 90-x_range,zrange=RATIO_V_RANGE,zlog=ratio_v_log,zticklen= -1,zticks=4, $
                                              ztitle= ratio_units ,charsize=1.2, title= sort_title+'  '+PHASE+ $
                                              ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cfrom ' $
                                              +ts_date+'  TO  ' +te_date
                                        endif  else begin 
                                            polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(event_ratio_2d), $
                                              r_range=90-y_range,zrange=ratio_v_range,zlog=ratio_v_log, $
                                              charsize=1.2,zticklen= -1, zticks=4,$
                                              ztitle= ratio_units, title= sort_title+'  '+PHASE+ $
                                              ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cfrom ' $
                                              +ts_date+'  TO  ' +te_date
                                        endelse 
                                        for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                        for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                                        IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                    endif 
                                    spawn, 'mogrify -format png '+ path_2d +'*.ps'
                               
                                    spawn, 'mogrify -rotate -90 '+ path_2d +'*.png'
                                    spawn, 'mv -f '+ path_2d +'*.png '+ path_2d +'png/'
                                    spawn, 'gzip -9f ' + path_2d +'*.ps'
                                  
                                ENDIF
;--------- slice-plots -----------
; z is the direction been sliced
                                IF KEYWORD_SET(SLICE_PLOT) THEN BEGIN 
                                    path_slice_main = path_ev+'slice/'
                                    spawn, 'mkdir ' +path_slice_main
;Calulation
                                    FOR iz = 0, nz-1 DO BEGIN 
                                        slice_block = [strcompress(STRING(z_axis(iz)-SLICE_GRID*0.5, format = '(f5.1)'),/remove_all), $
                                                       strcompress(STRING(z_axis(iz)+SLICE_GRID*0.5, format = '(f5.1)'),/remove_all)]
                                        slice_total_counts=total_counts(*,*,iz)
                                        slice_event_counts = event_counts(*, *, iz)
                                        slice_event_ratio = event_ratio(*, *, iz)
                                        index= where(slice_event_counts eq 0,ct)
                                        if ct gt 0 then slice_event_counts(index) = !VALUES.F_NAN
                                        index= where(slice_total_counts eq 0,ct)
                                        if ct gt 0 then slice_total_counts(index) = !VALUES.F_NAN

                                        IF keyword_set(symmetry_template) THEN BEGIN
                                            positive_half = slice_total_counts(*, where(y_axis GT 0)) GT 0
                                            positive_half = [[reverse(positive_half, 2)], [positive_half]]
                                            negative_half = slice_total_counts(*, where(y_axis LT 0)) GT 0
                                            negative_half = [[negative_half], [reverse(negative_half, 2)]]
                                            template_slice = FLOAT(positive_half*negative_half)
;draw the template
                                            draw_the_template = 0
                                            IF keyword_set(draw_the_template) THEN BEGIN 
                                                window, /free 
                                                specplot, x_axis, y_axis, template_slice, $
                                                  no_interp = 1, $
                                                  lim = { zlog:0, $
                                                          zrange: [0, 1], zticklen: -1,$
                                                          title: sort_title+'  '+PHASE+$
                                                          ' O!U+!N BEAM EVENTS in '+region_name_set(is) $
                                                          +' Region!Cfrom ' $
                                                          + ts_date+'  TO  ' +te_date,  $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          xrange: x_range, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.85, 0.85],$
                                                          zticks:4}   
                                            ENDIF
                                            template_slice(where(template_slice EQ 0 )) =  !VALUES.F_NAN
                                            slice_total_counts= slice_total_counts*template_slice
                                            slice_event_counts = slice_event_counts*template_slice
                                            slice_event_ratio = slice_event_ratio*template_slice
                                        ENDIF       
;slice
                                        path_slice = path_slice_main+'slice_zgrid_'+slice_grid_str+'/'
                                        spawn, 'mkdir '+ path_slice
                                        spawn, 'mkdir '+ path_slice +'png/'
;samples
                                        IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                            popen, path_slice+region_name_set(is) + '_samples_' $ 
                                              + ts_date+'_to_' + te_date+'_' $
                                              + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $ 
                                              +'_at_' + plot_axis(2)+'_' $
                                              +slice_block(0)+'_'+slice_block(1) $
                                              +'.ps', /land 
                                        ENDIF
                                        specplot, x_axis, y_axis,  slice_total_counts $
                                          , no_interp = 1, $
                                          lim = { zlog:EVENTS_V_LOG, zrange:EVENTS_V_RANGE, $
                                                  title: sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM SAMPLES in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ': [' +slice_block(0)+',' +slice_block(1)+']' $
                                                  + '!CFROM  ' + ts_date+'  TO  '  +te_date, $
                                                  xtitle:PLOT_AXIS(0),ytitle:PLOT_AXIS(1),$
                                                  ztitle:events_units, zticklen: -1,$
                                                  xrange: X_RANGE, yrange: Y_RANGE, $
                                                  XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                  position: [0.1, 0.1, 0.85, 0.85],$
                                                  zticks:4}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        IF KEYWORD_SET(ps_plot) THEN pclose ELSE stop
                                        IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                                          or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
                                          THEN BEGIN
                                            IF KEYWORD_SET(ps_plot) THEN    popen, path_slice+region_name_set(is) + '_samples_' $ 
                                              + ts_date+'_to_' + te_date+'_' $
                                              + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $ 
                                              +'_at_' + plot_axis(2)+'_' $
                                              +slice_block(0)+'_'+slice_block(1) $
                                              +'_polar.ps', /land
                                            if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                                polar_spec,(90-x_axis),y_axis*norm_factor_mlt,sliec_total_counts, $
                                                  r_range=90-x_range,zrange=events_v_range,charsize=1.2,$
                                                  ztitle=events_units,zlog=events_v_log,zticklen= -1,zticks=4,$
                                                  title= sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM SAMPLES in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ': [' +slice_block(0)+','+slice_block(1)+']' $
                                                  + '!CFROM  ' + ts_date+'  TO  ' +te_date
                                            endif  else begin 
                                                polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(slice_total_counts), $
                                                  r_range=90-y_range,zrange=events_v_range, charsize=1.2,$
                                                  ztitle=events_units,zlog=events_v_log,zticklen= -1,zticks=4,$
                                                  title= sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM SAMPLES in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ': [' +slice_block(0)+',' +slice_block(1)+']' $
                                                  + '!CFROM  ' + ts_date+'  TO  ' $
                                                  +te_date
                                            endelse
                                            for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                            for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                                            IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                        endif 

;events
                                        IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                            popen, path_slice+region_name_set(is) + '_events_' $ 
                                              + ts_date+'_to_' + te_date+'_' $
                                              + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $ 
                                              +'_at_' + plot_axis(2)+'_' $
                                              +slice_block(0)+'_'+slice_block(1) $
                                              +'.ps', /land 
                                        ENDIF
                                        specplot, x_axis, y_axis,  slice_event_counts $
                                          , no_interp = 1, $
                                          lim = { zlog:EVENTS_V_LOG, zrange:EVENTS_V_RANGE, $
                                                  title: sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM EVENTS in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ': [' +slice_block(0)+',' +slice_block(1)+']' $
                                                  + '!CFROM  ' + ts_date+'  TO  '  +te_date, $
                                                  xtitle:PLOT_AXIS(0),ytitle:PLOT_AXIS(1),$
                                                  ztitle:events_units, zticklen: -1,$
                                                  xrange: X_RANGE, yrange: Y_RANGE, $
                                                  XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                  position: [0.1, 0.1, 0.85, 0.85],$
                                                  zticks:4}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        IF KEYWORD_SET(ps_plot) THEN pclose ELSE stop
                                        IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                                          or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
                                          THEN BEGIN
                                            IF KEYWORD_SET(ps_plot) THEN    popen, path_slice+region_name_set(is) + '_events_' $ 
                                              + ts_date+'_to_' + te_date+'_' $
                                              + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $ 
                                              +'_at_' + plot_axis(2)+'_' $
                                              +slice_block(0)+'_'+slice_block(1) $
                                              +'_polar.ps', /land
                                            if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                                polar_spec,(90-x_axis),y_axis*norm_factor_mlt,sliec_event_counts, $
                                                  r_range=90-x_range,zrange=events_v_range,charsize=1.2,$
                                                  ztitle=events_units,zlog=events_v_log,zticklen= -1,zticks=4,$
                                                  title= sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM EVENTS in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ': [' +slice_block(0)+','+slice_block(1)+']' $
                                                  + '!CFROM  ' + ts_date+'  TO  ' +te_date
                                            endif  else begin 
                                                polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(slice_event_counts), $
                                                  r_range=90-y_range,zrange=events_v_range, charsize=1.2,$
                                                  ztitle=events_units,zlog=events_v_log,zticklen= -1,zticks=4,$
                                                  title= sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM EVENTS in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ': [' +slice_block(0)+',' +slice_block(1)+']' $
                                                  + '!CFROM  ' + ts_date+'  TO  ' $
                                                  +te_date
                                            endelse
                                            for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                            for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                                            IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                        endif 
;ratio                                                       
                                        IF KEYWORD_SET(ps_plot) THEN popen, path_slice+region_name_set(is)+'_ratio_' $    
                                          + ts_date+'_to_' + te_date+'_'$ 
                                          + PLOT_AXIS(0)+'_vs_' +PLOT_AXIS(1)+'_at_' + plot_axis(2)+'_' $
                                          + slice_block(0)+'_'+slice_block(1) $
                                          +'.ps', /land 
                                        specplot, x_axis, y_axis,  slice_event_ratio, no_interp = 1, $
                                          lim = { zlog: RATIO_V_LOG,zrange: RATIO_V_RANGE, $
                                                  title: sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ': [' +slice_block(0)+',' +slice_block(1)+']' $
                                                  + '!CFROM  ' + ts_date +'  TO  ' +te_date, $
                                                  xtitle:PLOT_AXIS(0),ytitle:PLOT_AXIS(1), $
                                                  ztitle: ratio_units,zticklen: -1, $
                                                  xrange: X_RANGE, yrange: Y_RANGE, $
                                                  XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                  position: [0.1, 0.1, 0.85, 0.85],$
                                                  zticks:4}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]    
                                        IF KEYWORD_SET(ps_plot) THEN pclose ELSE stop
                                        IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                                          or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
                                          THEN BEGIN
                                            IF KEYWORD_SET(ps_plot) THEN  popen, path_slice+region_name_set(is)+'_ratio_' $    
                                              + ts_date+'_to_' + te_date+'_'$ 
                                              + PLOT_AXIS(0)+'_vs_' +PLOT_AXIS(1)+'_at_' + plot_axis(2)+'_' $
                                              + slice_block(0)+'_'+slice_block(1) $
                                              +'_polar.ps', /land  
                                            
                                            if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                                polar_spec,(90-x_axis),y_axis*norm_factor_mlt,slice_event_ratio, $
                                                  r_range= 90-x_range,zrange= RATIO_V_RANGE,zlog=ratio_v_log, $
                                                  ztitle=ratio_units ,charsize=1.2, zticklen= -1, zticks=4,$
                                                  title=sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ': [' +slice_block(0)+',' +slice_block(1)+']' $
                                                  + '!CFROM  ' + ts_date +'  TO  ' +te_date
                                            endif  else begin 
                                                polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(slice_event_ratio), $
                                                  r_range=90-y_range,zrange= RATIO_V_RANGE,zlog=ratio_v_log, $
                                                  ztitle=ratio_units,charsize=1.2,zticklen= -1, zticks=4,$
                                                  title= sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ': [' +slice_block(0)+','  +slice_block(1)+']' $
                                                  + '!CFROM  ' + ts_date  +'  TO  ' +te_date
                                            endelse
                                            for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                            for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                                            IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                        endif 
                                        spawn, 'mogrify -format png '+ path_slice +'*.ps'
                                     
                                        spawn, 'mogrify -rotate -90 '+ path_slice +'*.png'
                                        spawn, 'mv -f '+ path_slice +'*.png '+path_slice+'png/'
                                        spawn, 'gzip -9f ' + path_slice+'*.ps'  
                                      
                                    ENDFOR       
                                ENDIF         
; ------- WATERDROP PLOT ------------  not useable right now
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
                                            wd_event_ratio = TOTAL(total_counts(*, *, index), 3, /nan)/ $
                                              TOTAL(event_counts(*, *, index), 3, /nan)
                                        ENDIF ELSE BEGIN
                                            wd_event_counts = event_counts(*, *, index)
                                            wd_event_ratio = event_ratio(*, *, index)
                                        ENDELSE     
                                        wd_event_counts(where(wd_event_counts EQ 0)) =  !VALUES.F_NAN
                                        
                                        slice_block = [ STRCOMPRESS(STRING(z_axis(index(0))-SLICE_GRID*0.5, format = '(f5.1)'),/REMOVE_ALL), $
                                                        STRCOMPRESS(STRING(z_axis(index(N_ELEMENTS(index)-1)) +SLICE_GRID*0.5,  format = '(f5.1)'),/REMOVE_ALL)]
;waterdrop plot
                                        IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                            path_wd = path_wd_main+'wd_zgrid_' + slice_grid_str +'/'
                                            spawn, 'mkdir ' + path_wd 
                                            spawn, 'mkdir ' + path_wd+'png/'
                                            popen,  path_wd+region_name_set(is)+'_events_'+ $
                                              ts_date+'_to_' + te_date+'_' $
                                              +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'_at_'+ plot_axis(2) $
                                              +'_'+slice_block(0)+'_'+slice_block(1) $
                                              +'.ps', /land 
                                        ENDIF
                                        specplot, x_axis, y_axis, wd_event_counts, $
                                          no_interp = 1, $
                                          lim = { zlog: EVENTS_V_LOG, zrange: EVENTS_V_RANGE, $
                                                  title: sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ':[ ' +slice_block(0) +','+slice_block(1)+' ] ' $
                                                  + '!CFROM  ' + ts_date +'  TO  ' +te_date, $
                                                  xtitle:PLOT_AXIS(0),ytitle:PLOT_AXIS(1),ztitle:events_units, $
                                                  xrange: X_RANGE, yrange: Y_RANGE, zticklen: -1,$
                                                  XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                  position: [0.1, 0.1, 0.85, 0.85],$
                                                  zticks:4}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                        IF KEYWORD_SET(ps_plot) THEN   pclose ELSE stop
                                        IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                                          or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
                                          THEN BEGIN
                                            IF KEYWORD_SET(ps_plot) THEN  popen,  path_wd+region_name_set(is)+'_events_'+ $
                                              ts_date+'_to_' + te_date+'_' $
                                              +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1)+'_at_'+ plot_axis(2) $
                                              +'_'+slice_block(0)+'_'+slice_block(1) $
                                              +'_polar.ps', /land 
                                            if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                                polar_spec,(90-x_axis),y_axis*norm_factor_mlt,wd_event_counts, $
                                                  r_range=90-x_range,zrange=events_v_range,charsize=1.2,$
                                                  ztitle=events_units,zlog=events_v_log,zticklen= -1,zticks=4,$
                                                  title= sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2)+ ':[ ' +slice_block(0) +','+slice_block(1)+' ] ' $
                                                  + '!CFROM  ' + ts_date  +'  TO  ' +te_date
                                            endif  else begin 
                                                polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(wd_event_counts), $
                                                  r_range=90-y_range,zrange=events_v_range,charsize=1.2, $
                                                  ztitle=events_units,zlog=events_v_log,zticklen= -1,zticks=4,$
                                                  title= sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ':[ ' +slice_block(0)  +','+slice_block(1)+' ] ' $
                                                  + '!CFROM  ' + ts_date +'  TO  ' +te_date
                                            endelse
                                            for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                            for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                                            IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                        endif 
;ratio
                                        IF KEYWORD_SET(ps_plot) THEN   popen,   path_wd+region_name_set(is)+'_ratio_'+ $ 
                                          ts_date+'_to_' + te_date+'_' +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'_at_'+ plot_axis(2) $
                                          +'_'+slice_block(0)+'_'+slice_block(1)  +'.ps', /land 
                                        specplot, x_axis, y_axis, wd_event_ratio, no_interp = 1, $
                                          lim = { zlog:RATIO_V_LOG, zrange:RATIO_V_RANGE,zticklen: -1, $
                                                  title: sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2)  + ':[ ' +slice_block(0) +','+slice_block(1)+' ] ' $
                                                  + '!CFROM  ' + ts_date  +'  TO  ' +te_date, $
                                                  xtitle:PLOT_AXIS(0),ytitle:PLOT_AXIS(1), ztitle:ratio_units, $
                                                  xrange: X_RANGE, yrange: Y_RANGE, XSTYLE:1, ystyle: 1, $
                                                  charsize: 1.2, position: [0.1, 0.1, 0.85, 0.85],$
                                                  zticks:4}   
                                        oplot, [0, 0, -100, 100], [-100, 100, 0, 0]    
                                        IF KEYWORD_SET(ps_plot) THEN pclose ELSE stop  
                                        IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                                          or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
                                          THEN BEGIN
                                            IF KEYWORD_SET(ps_plot) THEN popen,  path_wd+region_name_set(is)+'_ratio_'+ $ 
                                              ts_date+'_to_' + te_date+'_' +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1)+'_at_'+plot_axis(2) $
                                              +'_'+slice_block(0)+'_'+slice_block(1)  +'_polar.ps', /land 
                                            if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                                polar_spec,(90-x_axis),y_axis*norm_factor_mlt,sliced_event_ratio, $
                                                  r_range= 90-x_range,zrange= RATIO_V_RANGE,zlog=ratio_v_log, $
                                                  ztitle=ratio_units ,charsize=1.2, zticklen= -1,zticks=4,$
                                                  title= sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM RATIO in '+region_name_set(is) +' Region!Cat ' $
                                                  + plot_axis(2) + ':[ ' +slice_block(0)  +','+slice_block(1)+' ] ' $
                                                  + '!CFROM  ' + ts_date   +'  TO  ' +te_date
                                            endif  else begin 
                                                polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(sliced_event_ratio), $
                                                  r_range=90-y_range,zrange= RATIO_V_RANGE,zlog=ratio_v_log, $
                                                  ztitle=ratio_units ,charsize=1.2, zticklen= -1,zticks=4,$
                                                  title= sort_title+'  '+PHASE+$
                                                  ' O!U+!N BEAM RATIO in '+region_name_set(is) +'s Region!Cat ' $
                                                  + plot_axis(2) + ':[ ' +slice_block(0)  +','+slice_block(1)+' ] ' $
                                                  + '!CFROM  ' + ts_date  +'  TO  ' +te_date
                                            endelse
                                            for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                            for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                                            IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                        endif   
                                        spawn, 'mogrify -format png '+ path_wd +'*.ps'
                                       
                                        spawn, 'mogrify -rotate -90 '+ path_wd +'*.png'
                                        spawn, 'mv -f '+ path_wd +'*.png '+path_wd+'png/'           
                                        spawn, 'gzip -9f ' + path_slice+'*.ps'  
                                       
                                    ENDFOR             
                                ENDIF   
                            endelse     
                        ENDFOR 
                    ENDIF                                                                  
;--------------------------------- PROPERTY_MAP ---------------------
                    IF KEYWORD_SET(PROPERTY_MAP_SET)  THEN BEGIN 
;run for different property map type
                        FOR ipmt = 0, n_elements(property_map_type_set)-1 DO BEGIN 
                            property_map_type = property_map_type_set(ipmt)
                            IF TOTAL(property_map_type_set_all EQ  property_map_type) NE 1  THEN BEGIN 
                                stop
                                print, 'no property_map_type: '+property_map_type $
                                  +'exists!Cinput .cont to continue to next property_map_type'
                                GOTO,  next_property_map_type
                            ENDIF 
;run for different properties as set in property_map_set
                            FOR ipp = 0, N_ELEMENTS(PROPERTY_MAP_SET)-1 DO  BEGIN 
; setup the plot saving directory
                                path_pp_main = path+property_map_set(ipp)+'/' &  spawn, 'mkdir ' + path_pp_main
; input different data for different properties
                                property = property_map_set(ipp)
                                IF property EQ 'distfunc' THEN BEGIN 
                                    property_tail = distfunc_tail & property_earth = distfunc_earth
                                    PROPERTY_V_LOG = distfunc_V_LOG & PROPERTY_V_RANGE = distfunc_V_RANGE 
                                    property_units = distfunc_units
                                ENDIF
                                IF property EQ 'normal_distfunc' THEN BEGIN 
                                    property_tail = distfunc_tail/mag & property_earth = distfunc_earth/mag
                                    PROPERTY_V_LOG = normal_distfunc_V_LOG & PROPERTY_V_RANGE = normal_distfunc_V_RANGE 
                                    property_units = distfunc_units
                                ENDIF
                                IF property EQ 'energy' THEN BEGIN 
                                    property_tail = energy_tail & property_earth = energy_earth
                                    PROPERTY_V_LOG = ENERGY_V_LOG & PROPERTY_V_RANGE = ENERGY_V_RANGE 
                                    property_units = energy_units
                                ENDIF 
                                IF property EQ 'energy_v' THEN BEGIN 
                                    property_tail = energy_v_tail & property_earth = energy_v_earth
                                    PROPERTY_V_LOG = velocity_V_LOG & PROPERTY_V_RANGE = velocity_V_RANGE 
                                    property_units = velocity_units
                                ENDIF 
                                 IF property EQ 'v_energy' THEN BEGIN 
                                    property_tail = v_energy_tail & property_earth = v_energy_earth
                                    PROPERTY_V_LOG = energy_V_LOG & PROPERTY_V_RANGE = energy_V_RANGE 
                                    property_units = energy_units
                                ENDIF 
                                IF property EQ 'flux' THEN BEGIN 
                                    property_tail = flux_tail  &  property_earth = flux_earth
                                    PROPERTY_V_LOG = FLUX_V_LOG &  PROPERTY_V_RANGE = FLUX_V_RANGE 
                                    property_units = FLUX_units
                                ENDIF 
                                IF property EQ 'eflux' THEN BEGIN 
                                    property_tail = eflux_tail  &  property_earth = eflux_earth
                                    PROPERTY_V_LOG = EFLUX_V_LOG &  PROPERTY_V_RANGE = EFLUX_V_RANGE 
                                    property_units = EFLUX_units
                                ENDIF 
                                IF property EQ 'density' THEN BEGIN 
                                    property_tail = density_tail &   property_earth = density_earth
                                    PROPERTY_V_LOG = density_V_LOG &   PROPERTY_V_RANGE = density_V_RANGE 
                                    property_units = density_units
                                ENDIF 
                                IF property EQ 'velocity' THEN BEGIN 
                                    property_tail = velocity_tail &  property_earth = velocity_earth
                                    PROPERTY_V_LOG = velocity_V_LOG &   PROPERTY_V_RANGE = velocity_V_RANGE 
                                    property_units = velocity_units
                                ENDIF 
                                IF property EQ 'Vpara' THEN BEGIN 
                                    property_tail = ABS(v_para_tail) & property_earth = ABS(v_para_earth)
                                    PROPERTY_V_LOG = velocity_Vpara_LOG &   PROPERTY_V_RANGE = velocity_Vpara_RANGE 
                                    property_units = velocity_units
                                ENDIF 
                                IF property EQ 'Vperp' THEN BEGIN 
                                    property_tail = v_perp_tail &  property_earth = v_perp_earth
                                    PROPERTY_V_LOG = velocity_Vperp_LOG &   PROPERTY_V_RANGE = velocity_Vperp_RANGE 
                                    property_units = velocity_units
                                ENDIF 
                                IF property EQ 'nV' THEN BEGIN 
                                    property_tail = density_tail*velocity_tail*1e9 &  property_earth = density_earth*velocity_earth*1e9
                                    PROPERTY_V_LOG = nv_V_LOG &  PROPERTY_V_RANGE = nv_V_RANGE 
                                    property_units = nv_units
                                ENDIF 
                                IF property EQ 'pitch_angle' THEN BEGIN 
                                    property_tail = pitch_angle_tail & property_earth = pitch_angle_earth
                                    PROPERTY_V_LOG = pitch_angle_V_LOG & PROPERTY_V_RANGE = pitch_angle_V_RANGE 
                                    property_units = pitch_angle_units
                                ENDIF 
                                IF property EQ 'nVpara_over_B' THEN BEGIN 
                                    property_tail = (density_tail*ABS(v_para_tail))/mag & property_earth = (density_earth*ABS(v_para_earth))/mag
                                    PROPERTY_V_LOG = nvpara_over_b_v_log & PROPERTY_V_RANGE = nvpara_over_b_v_range
                                    property_units = nvpara_over_b_units
                                ENDIF 
                                IF property EQ 'temperature' THEN BEGIN 
                                    property_tail = ABS(temperature_tail) & property_earth = ABS(temperature_earth)
                                    PROPERTY_V_LOG = temperature_v_log & PROPERTY_V_RANGE = temperature_v_range
                                    property_units = temperature_units
                                ENDIF
                                IF property EQ 'anodes' THEN BEGIN 
                                    property_tail = 4.5-ABS(4.5-anodes_tail) & property_earth = 4.5-ABS(4.5-anodes_earth)
                                    PROPERTY_V_LOG = anodes_v_log &  PROPERTY_V_RANGE = anodes_v_range
                                    property_units = anodes_units
                                ENDIF 
                                IF property EQ 'By' THEN BEGIN 
                                    property_tail = By & property_earth = replicate(!values.f_nan,n_elements(by))
                                    PROPERTY_V_LOG = By_v_log &  PROPERTY_V_RANGE = By_v_range
                                    property_units = By_units
                                ENDIF 
                                IF property EQ 'Beta' THEN BEGIN 
                                    property_tail = Beta & property_earth = replicate(!values.f_nan,n_elements(Beta))
                                    PROPERTY_V_LOG = Beta_v_log &  PROPERTY_V_RANGE = Beta_v_range
                                    property_units = Beta_units
                                ENDIF 
                                
                                IF NOT keyword_set(property_tail) THEN BEGIN 
                                    print, 'no ' +property+' avaliabe for mapping'
                                    stop
                                    GOTO, next_property
                                ENDIF 
;----- Calculation -------     
                                IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                                  or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) THEN BEGIN
                                    if ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT'])  then begin 
                                        grid_x= 1*grid & grid_y=2.5*grid
                                    endif else begin 
                                        grid_x=1*grid & grid_y=2.5*grid
                                    endelse 
                                endif  else begin 
                                    grid_x=grid & grid_y=grid
                                endelse 
                                nx = CEIL(ABS(x_range(1) - x_range(0))/grid_x)
                                ny = CEIL(ABS(y_range(1) - y_range(0))/grid_y)
                                nz = CEIL(ABS(z_range(1)-z_range(0))/SLICE_GRID) 
                                x_axis = INDGEN(nx)*grid_x+(x_range(0) < x_range(1))+grid_x*0.5
                                y_axis = INDGEN(ny)*grid_y+(y_range(0) < y_range(1))+grid_y*0.5
                                z_axis = INDGEN(nz)*slice_grid+(z_range(0) < z_range(1))$
                                  +slice_grid*0.5

                                IF NOT KEYWORD_SET(diff_beta) THEN  diff_beta = [0]
                                ns = n_elements(diff_beta)
                                FOR iis = 0, ns-1 DO BEGIN 
                                    is = diff_beta(iis)
                                    IF is EQ 0 THEN ind_beam = where(ABS(flag) GT 0, ct_beam)
                                    IF is EQ 1 THEN ind_beam = where(ABS(flag) GT 0 AND beta LE 0.05, ct_beam)
                                    IF is EQ 2 THEN ind_beam = where(ABS(flag) GT 0 AND beta LE 1 AND beta GT 0.05, ct_beam)
                                    IF is EQ 3 THEN ind_beam = where(ABS(flag) GT 0 AND beta GT 1, ct_beam)
                                    IF is EQ 9 THEN ind_beam = where(ABS(flag) GE 0 AND beta LE 0.02, ct_beam)
                                        
                                    IF ct_beam eq 0 then begin 
                                        OPENU, unit, path+'log_nodata.txt', /GET_LUN, /APPEND
                                        PRINTF, unit, 'no data for'+ region_name_set(is)
                                        FREE_LUN, unit     
                                    endif else begin 
                                        property_value = REPLICATE(!VALUES.F_NAN, ct_beam, 2, nx, ny, nz)
                                        FOR i = 0l, ct_beam-1 DO BEGIN 
                                            index = SORT(ABS(x_axis - data_pos(ind_beam(i), 0))) &  index_x = index(0)
                                            index = SORT(ABS(y_axis - data_pos(ind_beam(i), 1))) &  index_y = index(0)
                                            index = SORT(ABS(z_axis - data_pos(ind_beam(i), 2))) &  index_z = index(0)
                                            
                                            property_value(i, *, index_x, index_y, index_z) = $
                                              [property_tail(ind_beam(i)), property_earth(ind_beam(i))]
                                        ENDFOR  
;-------- 2d plot --------
                                        IF KEYWORD_SET(PLOT_2D) THEN BEGIN   
                                            property_map_2d = cal_property_map_value( property_value, property_map_type)
                                            path_pp_2d = path_pp_main+'2d/'
                                            IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                                spawn, 'mkdir ' + path_pp_2d 
                                                spawn, 'mkdir '+ path_pp_2d +'png/'
                                                popen,  path_pp_2d+property_map_type+'_' $
                                                  +property+'_'+region_name_set(is)+'_'+ $ ; all 
                                                  ts_date+'_to_' + te_date+'_' $
                                                  +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                  +'.ps', /land 
                                            ENDIF 
                                            specplot, x_axis, y_axis, property_map_2d, $
                                              no_interp = 1, $
                                              lim = { zlog:PROPERTY_V_LOG, $
                                                      zrange:PROPERTY_V_RANGE, $
                                                      title: sort_title+'  '+PHASE+$
                                                      ' O!U+!N  beam!C' +property_map_type +' ' +PROPERTY $
                                                      +' distribution in '+region_name_set(is)+' Region !Cfrom  ' $
                                                      + ts_date+'  TO  ' +te_date,  $
                                                      xtitle: PLOT_AXIS(0), ytitle: PLOT_AXIS(1), $
                                                      ztitle: property_units,zticklen: -1, $
                                                      xrange: X_RANGE, yrange: Y_RANGE, $
                                                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                                      position: [0.1, 0.1, 0.85, 0.85],$
                                                      zticks:4}   
                                            oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
; try to overplot the circle of magnetosphere in order to define solawind
                                ;                                    x = indgen(91)*0.1-1
                                ;                                   x = [x, !values.f_nan]
                                ;                                  y = sqrt((1-x^2/8^2)*(15^2))+1
                                ;                                 x = [x, x]
                                ;                                y = [-y, y]
                                ;                               oplot,
;                               x, y
                                            IF keyword_set(ps_plot) THEN pclose ELSE stop

                                            IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                                              or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
                                              THEN BEGIN
                                                IF KEYWORD_SET(ps_plot) THEN   popen,  path_pp_2d+property_map_type+'_' $
                                                  +property+'_'+region_name_set(is)+'_'+ $ ; all 
                                                  ts_date+'_to_' + te_date+'_' $
                                                  +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                  +'_polar.ps', /land 
                                                
                                                if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                                    polar_spec,(90-x_axis),y_axis*norm_factor_mlt,property_map_2d, $
                                                      r_range=90-x_range,zrange=property_v_range,charsize=1.2,$
                                                      ztitle=property_units, zlog=property_v_log,zticklen= -1,zticks=4,$
                                                      title=sort_title+'  '+PHASE+$
                                                      ' O!U+!N  beam!C' +property_map_type +' ' +PROPERTY $
                                                      +' distribution in '+region_name_set(is)+' Region !Cfrom  ' $
                                                      + ts_date+'  TO  ' +te_date
                                                endif  else begin 
                                                    polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(property_map_2d), $
                                                      r_range=90-y_range,zrange=property_v_range, charsize=1.2,$
                                                      ztitle=property_units,zlog=property_v_log,zticklen= -1,zticks=4,$
                                                      title=sort_title+'  '+PHASE+$
                                                      ' O!U+!N  beam!C' +property_map_type +' ' +PROPERTY $
                                                      +' distribution in '+region_name_set(is)+' Region !Cfrom  ' $
                                                      + ts_date+'  TO  ' +te_date
                                                endelse
                                                for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                                for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar

                                                IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                            endif 

                                            spawn, 'mogrify -format png '+ path_pp_2d +'*.ps'
                                        
                                            spawn, 'mogrify -rotate -90 '+ path_pp_2d +'*.png'
                                          
                                            spawn, 'mv -f '+ path_pp_2d +'*.png '+path_pp_2d+'png/'
                                            spawn, 'gzip -9f ' + path_pp_2d+'*.ps'                     
                                          
                                        ENDIF
; --------------- slice plots --------------------
; z is the direction been sliced
                                        IF KEYWORD_SET(SLICE_PLOT) THEN BEGIN                                            
                                            path_pp_slice_main = path_pp_main+'slice/'
                                            spawn, 'mkdir ' + path_pp_slice_main
                                            FOR iz = 0, nz-1 DO BEGIN 
                                                slice_block = [STRCOMPRESS(STRING(z_axis(iz)-SLICE_GRID*0.5, $
                                                                                  format = '(f5.1)'),/REMOVE_ALL), $
                                                               STRCOMPRESS(STRING(z_axis(iz)+SLICE_GRID*0.5, $
                                                                                  format = '(f5.1)'),/REMOVE_ALL)]
                                                property_map_slice = $
                                                  cal_property_map_value(property_value(*, *, *, *, iz), property_map_type)
; draw the plot
                                                path_pp_slice = path_pp_slice_main +'slice_zgrid_'+slice_grid_str+'/'
                                                spawn, 'mkdir '+ path_pp_slice
                                                spawn, 'mkdir '+ path_pp_slice +'png/'
                                                IF KEYWORD_SET(ps_plot) THEN  popen, path_pp_slice+property_map_type+'_' $
                                                  +property+'_'+region_name_set(is)+'_' $
                                                  + ts_date+'_to_' + te_date+'_' $
                                                  + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $ 
                                                  +'_at_' + plot_axis(2)+'_' $
                                                  +slice_block(0)+'_'+slice_block(1) $
                                                  +'.ps', /land 

                                                specplot, x_axis, y_axis,  property_map_slice, $
                                                  no_interp = 1, $
                                                  lim = { zlog: PROPERTY_V_LOG, $
                                                          zrange:PROPERTY_V_RANGE, $
                                                          title: sort_title+'  '+PHASE+$
                                                          ' O!U+!N  beam!C' +property_map_type +' '+ PROPERTY $
                                                          +' distribution in ' +region_name_set(is) $
                                                          +' at  '+ plot_axis(2) $
                                                          + ':[ ' +slice_block(0)+','  +slice_block(1)+' ] ' $
                                                          + '!CFROM  ' + ts_date  +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          ztitle: property_units, zticklen: -1,$
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, $
                                                          charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.85, 0.85],$
                                                          zticks:4}   
                                                oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                                IF KEYWORD_SET(ps_plot) THEN pclose ELSE stop
                                                IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                                                  or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
                                                  THEN BEGIN
                                                    IF KEYWORD_SET(ps_plot) THEN     popen, path_pp_slice+property_map_type+'_' $
                                                      +property+'_'+region_name_set(is)+'_' $
                                                      + ts_date+'_to_' + te_date+'_' $
                                                      + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $ 
                                                      +'_at_' + plot_axis(2)+'_' $
                                                      +slice_block(0)+'_'+slice_block(1) $
                                                      +'_polar.ps', /land 
                                                    if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                                        polar_spec,(90-x_axis),y_axis*norm_factor_mlt,property_map_slice, $
                                                          r_range=90-x_range,zrange=property_v_range,charsize=1.2,$
                                                          ztitle=property_units, zlog=property_v_log,zticklen= -1,zticks=4,$
                                                          title= sort_title+'  '+PHASE+$
                                                          ' O!U+!N  beam!C' +property_map_type +' '+ PROPERTY $
                                                          +' distribution in ' +region_name_set(is) $
                                                          +' at  '+ plot_axis(2) $
                                                          + ':[ ' +slice_block(0)+','  +slice_block(1)+' ] ' $
                                                          + '!CFROM  ' + ts_date  +'  TO  ' +te_date
                                                    endif  else begin 
                                                        polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(property_map_slice), $
                                                          r_range=90-y_range,zrange=property_v_range,charsize=1.2, $
                                                          ztitle=property_units,zlog=property_v_log,zticklen= -1,zticks=4,$
                                                          title= sort_title+'  '+PHASE+$
                                                          ' O!U+!N  beam!C' +property_map_type +' '+ PROPERTY $
                                                          +' distribution in ' +region_name_set(is) $
                                                          +' at  '+ plot_axis(2) $
                                                          + ':[ ' +slice_block(0)+','  +slice_block(1)+' ] ' $
                                                          + '!CFROM  ' + ts_date  +'  TO  ' +te_date
                                                    endelse
                                                    for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
                                                    for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar
                                                    IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                                endif 
                                                spawn, 'mogrify -format png '+ path_pp_slice +'*.ps'
                                               
                                                spawn, 'mogrify -rotate -90 '+ path_pp_slice +'*.png'
                                               
                                                spawn, 'mv -f '+ path_pp_slice +'*.png '+path_pp_slice+'png/'
                                                spawn, 'gzip -9f ' + path_pp_slice+'*.ps'     
                                            ENDFOR          
                                        ENDIF    
; ------ WATER_DROP PLOT ------
;SLICE_GRID WILL BE ALSO USED HERE
                                        IF KEYWORD_SET(WATERDROP_PLOT) THEN BEGIN 
                                            path_pp_wd_main = path_pp_main +'waterdrop/'
                                            spawn, 'mkdir ' + path_pp_wd_main
;run for all slices
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
                                                slice_block = [STRCOMPRESS(STRING(z_axis(index(0))-SLICE_GRID*0.5, $
                                                                                  format = '(f5.1)'),/REMOVE_ALL), $
                                                               STRCOMPRESS(STRING(z_axis(index(N_ELEMENTS(index)-1))$
                                                                                  + SLICE_GRID*0.5, $
                                                                                  format = '(f5.1)'),/REMOVE_ALL)]
                                                IF KEYWORD_SET(ps_plot) THEN BEGIN 
                                                    path_pp_wd = path_pp_wd_main $
                                                      +'wd_zgrid_' + slice_grid_str +'/'
                                                    spawn, 'mkdir ' + path_pp_wd 
                                                    spawn, 'mkdir '+ path_pp_wd +'png/'
                                                    popen,  path_pp_slice+property_map_type+'_' $
                                                      +property+'_'+region_name_set(is)+'_' + $
                                                      ts_date+'_to_' + te_date+'_' $
                                                      +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                      +'_at_'+ plot_axis(2) $
                                                      +'_'+slice_block(0)+'_'+slice_block(1) $
                                                      +'.ps', /land 
                                                ENDIF
                                                specplot, x_axis, y_axis, wd_property, $
                                                  no_interp = 1, $
                                                  lim = { zlog:PROPERTY_V_LOG, $
                                                          zrange: PROPERTY_V_RANGE, $
                                                          title: sort_title+'  '+PHASE+$
                                                          ' O!U+!N  beam!C' +property_map_type +' '+ PROPERTY $
                                                          +' distribution in ' +region_name_set(is) $
                                                          +' at  '+ plot_axis(2) $
                                                          + ':[ ' +slice_block(0)+','  +slice_block(1)+' ] ' $
                                                          + '!CFROM  ' + ts_date  +'  TO  ' +te_date, $
                                                          xtitle: PLOT_AXIS(0), $
                                                          ytitle: PLOT_AXIS(1), $
                                                          ztitle: property_units,zticklen: -1, $
                                                          xrange: X_RANGE, yrange: Y_RANGE, $
                                                          XSTYLE:1, ystyle: 1, $
                                                          charsize: 1.2, $
                                                          position: [0.1, 0.1, 0.85, 0.85],$
                                                          zticks:4}   
                                                oplot, [0, 0, -100, 100], [-100, 100, 0, 0]
                                                IF KEYWORD_SET(ps_plot) THEN pclose ELSE stop
                                                IF ARRAY_EQUAL(PLOT_AXIS(0:1), ['MLT', 'ILAT']) $
                                                  or  ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) $
                                                  THEN BEGIN
                                                    IF KEYWORD_SET(ps_plot) THEN   popen,  path_pp_slice+property_map_type+'_' $
                                                      +property+'_'+region_name_set(is)+'_' + $
                                                      ts_date+'_to_' + te_date+'_' $
                                                      +PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) $
                                                      +'_at_'+ plot_axis(2) $
                                                      +'_'+slice_block(0)+'_'+slice_block(1) $
                                                      +'_polar.ps', /land 
                                                    if ARRAY_EQUAL(PLOT_AXIS(0:1), ['ILAT', 'MLT']) then begin 
                                                        polar_spec,(90-x_axis),y_axis*norm_factor_mlt,property_map_slice, $
                                                          r_range=90-x_range,zrange=property_v_range,charsize=1.2,$
                                                          ztitle=property_units, zlog=property_v_log,zticklen= -1,zticks=4,$
                                                          title=sort_title+'  '+PHASE+$
                                                          ' O!U+!N  beam!C' +property_map_type +' '+ PROPERTY $
                                                          +' distribution in ' +region_name_set(is) $
                                                          +' at  '+ plot_axis(2) $
                                                          + ':[ ' +slice_block(0)+','  +slice_block(1)+' ] ' $
                                                          + '!CFROM  ' + ts_date  +'  TO  ' +te_date
                                                    endif  else begin 
                                                        polar_spec,(90-y_axis),x_axis*norm_factor_mlt,transpose(property_map_slice), $
                                                          r_range=90-y_range,zrange=property_v_range, charsize=1.2,$
                                                          ztitle=property_units,zlog=property_v_log,zticklen= -1,zticks=4,$
                                                          title= sort_title+'  '+PHASE+$
                                                          ' O!U+!N  beam!C' +property_map_type +' '+ PROPERTY $
                                                          +' distribution in ' +region_name_set(is) $
                                                          +' at  '+ plot_axis(2) $
                                                          + ':[ ' +slice_block(0)+','  +slice_block(1)+' ] ' $
                                                          + '!CFROM  ' + ts_date  +'  TO  ' +te_date
                                                    endelse
                                                    IF KEYWORD_SET(ps_plot) THEN    pclose  ELSE stop
                                                endif 
                                                spawn, 'mogrify -format png '+ path_pp_wd +'*.ps'
                                                spawn, 'mogrify -rotate -90 '+ path_pp_wd +'*.png'
                                                spawn, 'mv -f '+ path_pp_wd +'*.png '+path_pp_wd+'png/'
                                                spawn, 'gzip -9f ' + path_pp_wd+'*.ps'     
                                            ENDFOR                
                                        ENDIF
                                    endelse    
                                ENDFOR         
                                next_property:    
                            ENDFOR
                            next_property_map_type:
                        ENDFOR       
                    ENDIF ;------------finish property map--            
                ENDFOR 
            endfor     
        ENDFOR    
    ENDFOR          
ENDFOR
if keyword_set(make_table) then begin
    !p.multi=[0,1,2]
    if keyword_set(averaged_with_ratio) then begin 
        if keyword_set(ps_plot) then popen,plot_path+'mean_ygrid'+slice_grid_str+'.ps',/port else  window,/free
        for ix=0, n_elements(selected_region_name_x)-1 do begin      
            plot,[0,0],[0,0],/nodata,xrange=[-1,N_ELEMENTS(storm_phase_set)],yrange=[0,1],xtitle='storm phases',ytitle='Occurrence Frequenct',title=selected_region_name_x(ix)
            for iy=0, n_elements(selected_region_name_y)-1 do begin
                for iz=0, n_elements(selected_region_name_z)-1 do begin
                    oplot,selected_region_ratio_mean(*,ix,iy,iz),psym=-(iy+1),color=2*(iz+1),symsize=2,thick=6 
                    xyouts, -0.8, 0.95-iz*0.05,selected_region_name_z(iz),color=2*(iz+1),charsize=1.2
                    xyouts, -0.88, 0.8-iy*0.05,selected_region_name_y(iy),charsize=1.2
                    oplot,[-0.28,-0.28],[0.81-iy*0.05,0.81-iy*0.05],psym=-(iy+1),symsize=1
                endfor 
            endfor 
        endfor 

        if keyword_set(ps_plot) then pclose else stop
        if keyword_set(ps_plot) then popen, plot_path+'median_ygrid'+slice_grid_str+'.ps',/port else window,/free
        for ix=0,1 do begin
            plot,[0,0],[0,0],/nodata,xrange=[-1,3],yrange=[0,1],ytitle='Occurrence Frequenct',title=selected_region_name_x(ix),xtickname=[' ','non-storm time','storm main phase','recovery phase',' ']
            for iy=0,1 do begin
                for iz=0,2 do begin             
                    oplot,selected_region_ratio_median(*,ix,iy,iz),psym=-(iy+1),color=2*(iz+1),symsize=2,thick=6
                    xyouts, -0.8, 0.95-iz*0.05,selected_region_name_z(iz),color=2*(iz+1),charsize=1.2
                    xyouts, -0.88, 0.8-iy*0.05,selected_region_name_y(iy),charsize=1.2
                    oplot,[-0.28,-0.28],[0.81-iy*0.05,0.81-iy*0.05],psym=-(iy+1),symsize=1
                endfor 
            endfor 
        endfor
        if keyword_set(ps_plot) then pclose else stop
    endif  else begin 
        if  not (selected_region_name_y eq 'all' and selected_region_name_z eq 'all') then begin
            for iy=0, n_elements(selected_region_name_y)-1 do begin 
                if keyword_set(ps_plot) then popen, plot_path+selected_region_name_y(iy)+'_table_plot.ps',/port else window,/free
                for ix=0, n_elements(selected_region_name_x)-1 do begin                    
                    plot,[0,0],[0,0],/nodata,xrange=[-1,3],yrange=[0.1,1],ytitle='Occurrence Frequenct',title=selected_region_name_y(iy)+' '+selected_region_name_x(ix),xtickname=[' ','non-storm','main phase','recovery phase',' '],charsize=1.2,ylog=1
                    for iz=0, n_elements(selected_region_name_z)-1 do begin             
                        oplot,[0,1,2],selected_region_ratio(*,ix,iy,iz),color=2*(iz+1),symsize=2,thick=6
                        xyouts, -0.8, 0.95-iz*0.08,selected_region_name_z(iz),color=2*(iz+1),charsize=1.5
                        errplot, [0,1,2],$
                          selected_region_ratio(*,ix,iy,iz)-selected_region_ratio_error(*,ix,iy,iz),$
                          selected_region_ratio(*,ix,iy,iz)+selected_region_ratio_error(*,ix,iy,iz), col = 2*(iz+1), thick = 4
                    endfor
                endfor 
                if keyword_set(ps_plot) then pclose else stop
            endfor 
        endif else begin 
            !p.multi=[0,1,1] 
            for iy=0, n_elements(selected_region_name_y)-1 do begin 
                if keyword_set(ps_plot) then popen, plot_path+'table_plot.ps',/land
                plot,[0,0],[0,0],/nodata,xrange=[-0.5,n_elements(storm_phase_set)-0.5],yrange=[0,1],ytitle='Occurrence Frequency',xtickname=storm_phase_set,charsize=1,ylog=0,position=[0.15,0.1,0.80,0.95],xstyle=1,ystyle=1
                
                x_axis=indgen(n_elements(storm_phase_set))
                for ix=0, n_elements(selected_region_name_x)-1 do begin
                    for iz=0, n_elements(selected_region_name_z)-1 do begin           
                        oplot,x_axis,selected_region_ratio(*,ix,iy,iz),color=2*(ix+1),symsize=2,thick=10
                        xyouts, 0, 0.7-ix*0.2,selected_region_name_x(ix),color=2*(ix+1),charsize=4
                        errplot, x_axis,$
                          selected_region_ratio(*,ix,iy,iz)-selected_region_ratio_error(*,ix,iy,iz),$
                          selected_region_ratio(*,ix,iy,iz)+selected_region_ratio_error(*,ix,iy,iz), col = 2*(ix+1), thick = 8
                    endfor
                endfor 
                if keyword_set(ps_plot) then pclose else stop
            endfor
   !p.multi=[0,1,2] 
        endelse  
    endelse 
    !p.multi=[0,1,1] 
    if keyword_set(ps_plot) then begin 
        spawn, 'mogrify -format png '+plot_path+'*.ps'
    endif  
endif  
END
