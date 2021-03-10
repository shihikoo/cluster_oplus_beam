PRO sort_o_beam_map, stop = stop,ps_plot=ps_plot,short_time_plot=short_time_plot
sc = 4  & sc_str = STRING(sc, format = '(i1.1)')

time_start = '2001-01-01/00:00:00' 
;time_start='2005-09-02/10'
time_end = '2002-12-31/12:59:59'
;time_end='2005-09-03'
;inst_name = 'HIA'
inst_name = 'CODIF'
if inst_name eq 'HIA' then sc =1
IF inst_name EQ 'CODIF' and sc eq 4 THEN path = 'output/o_beam/new_beam_filter/'
;IF inst_name EQ 'CODIF' and sc eq 4 THEN path = 'output/o_beam/sc4_5min_ar225/'
IF inst_name EQ 'CODIF' and sc eq 1 then path = 'output/o_beam/test_year_dependence_with_HIA/sc1_hiabeta_lobe/'
IF inst_name EQ 'HIA' THEN path = 'output/o_beam/test_year_dependence_with_HIA/hia_hiabeta_lobe/'
average_time = 300
;--- keywords settings ---
data_read_from_dat = 0 & save_data = 0
sort_recalc = 0 & save_sort = 0

non_sort = 0 & sort_season = 0
sort_kp = 0 & sort_ae = 0 & sort_f107 = 0
sort_imf_by = 0 & sort_imf_bz = 0 & sort_IMF_By_Bz = 0 & sort_IMF_By_smallBx = 0
sort_sw_p = 0 & sort_IMF_B = 0 & sort_sw_v = 0
sort_season_imf_by = 0 & sort_long_imf_bz = 0 & sort_long_imf_by_bz = 0
sort_substorm = 0

make_table = 0 & sort_imf_by_with_eflux_filter = 0 & eflux_filter = 0 & energy_filter = 0
en_vs_distfunc = 2

sort_by = 0 & sort_bx = 0 &  en_vs_beta = 0

map = 0 & sort_distribution = 0 & year_distribution = 0 & property_distribution = 1 & sort_anodes = 0

; 'all_direcitons','all_directions', 'tailward', 'earthward', 'both'
direction_set = ['all_directions'] ;['all_directions'] ;'tailward','earthward','both']

storm_phase_set = ['all_time'] ;,'storm_time','nonstorm_time'] ;,'main_phase','recovery','iniital_phase'] ;, 'main_phase', 'recovery','storm_time','all_time']

;storm_phase_set = ['nonstorm_time', 'prestorm', 'initial_phase', 'main_phase', 'recovery', 'storm_time','all_time']

plot_2d = 0 & slice_plot = 0 & waterdrop_plot = 0

grid_set = [1]                  ;, 2.]
slice_grid_set = [20.]           ;,20.]

point_plot = 1 & events_map = 0
property_map_set = ['density','velocity','nV','temperature'];];['normal_distfunc'];,'distfunc']
;property_map_set =['energy', 'velocity', 'flux', 'density', 'nVpara_over_B', 'temperature','anodes','eflux','nV','pitch_angle','Vperp','Vpara','IMF_By','sw_p','sw_v','VxBz','normal_distfunc','distfunc','v_energy','energy_v','E']
property_map_type_set = ['median'] ; ['mean','median','peak','varition']
diff_beta =[0] ; ['ALL', 'LOBE', 'BL', 'PS','le1','le05','le01','gt005le01','gt01','lt002']  ; for map only

coor_set = [['X_GSM','Y_GSM'] $;,['X_GSM','Y_GSM'] , ['Y_GSM', 'Z_GSM'] $
;  , ['MLT', 'ILAT'] $ 
;       , ['X_GSM_mod', 'Z_GSM_mod'], ['X_GSM_mod', 'Y_GSM_mod'], ['Y_GSM_mod', 'Z_GSM_mod'] $
                                ;           ['X_GSE','Y_GSE']  ,['X_GSE','Z_GSE'],['Y_GSE','Z_GSE'] $
]
if not keyword_set(ps_plot) then ps_plot = 1 else begin 
    if ps_plot lt 0 then ps_plot =0 else ps_plot=1
endelse
;---- basic settings ---
no_magnetosheath = 1
no_solarwind = 1

delete_uncertain_time = 1
delete_HIA_wrong_time = 1
symmetry_template = 0
delay_to_cluster=0
pre1h = 1

mass_o = 16*1.6e-27*(1e3)^2/(1.6e-19) ; unit: ev/(km/s)^2

mag_normal = 31200.*(6370./(6370+1000))^3*sqrt(1+3*sin(80*3.1415926/180)^2) ; =39832.1 ; dipole field at 1000km and 80 invariant latitude ;33695.9 ;nT at 60 invariant latitude degree from Seki 1998
;--- set the time info ---
ts = time_double(time_start)
te = time_double(time_end)
dt = te-ts

timespan, time_start, dt,  /SECONDS

ts_str = time_struct(ts)        ; start_time tplot time structure
te_str = time_struct(te)        ; end_time tplot time structure
jd_s = julday(ts_str.month, ts_str.date, ts_str.year) ;start julian day
jd_e = julday(te_str.month, te_str.date, te_str.year) ;end julian day

nyear = te_str.year-ts_str.year

time_str = strcompress(ts_str.year, /remove_all)$
  +string(ts_str.month, format = '(i2.2)') $
  +string(ts_str.date, format = '(i2.2)')+'_to_' $
  + strcompress(te_str.year, /remove_all) $
  +strcompress(te_str.month, /remove_all) $
  +strcompress(te_str.date, /remove_all)

ndays = (jd_e - jd_s) + 1       ; number of days to be loaded
;Last day is not included if hour=min=sec=0
IF te_str.hour EQ 0 AND te_str.min EQ 0 AND te_str.sec EQ 0 THEN ndays = ndays - 1
;eflux threshold setting
eflux_threshold = 0. & eflux_threshold_str = STRING(eflux_threshold, format = '(i4.4)')
;eflux_threshold = 2000. & eflux_threshold_str = STRING(eflux_threshold, format = '(i4.4)')
energy_threshold =[0,100.]  & energy_threshold_str = strcompress(STRING(energy_threshold, format = '(i5.5)'),/remove_all)
;------------------------------------------
;check prestore data 
;------------------------------------
flndata = path+	'tplot_restore/data_'+time_str
print, FINDFILE(flndata+'.tplot', COUNT = ct_data)

IF KEYWORD_SET(data_read_from_dat) THEN ct_data = 0
IF ct_data GT 0 THEN BEGIN 
    tplot_restore, filenames = flndata+'.tplot' 
    get_data, 'data', data = data
    time = data.x
    title = data.title
    data = data.y
    ntime = n_elements(time)
ENDIF ELSE BEGIN 
;-----------------------
;if there is no tplot store data then read from .dat files
;--------------------------
;--- known titles info in data files ---
    title =  [ '         flag  ', $ ;0
               '         Beta  ', $ ;1
               '      GSE_X(Re)', $ ;2
               '      GSE_Y(Re)', $ ;3
               '      GSE_Z(Re)', $ ;4
               '       en_tail ', $ ;5
               '       pa_tail ', $ ;6
               '      flux_tail', $ ;7
               '   Density_tail', $ ;8
               '   V_total_tail', $ ;9
               '     V_par_tail', $ ;10
               '    V_perp_tail', $ ;11
               '   T_total_tail', $ ;12
               '   P_total_tail', $ ;13
               '       en_earth', $ ;14
               '       pa_earth', $ ;15
               '     flux_earth', $ ;16
               '  Density_earth', $ ;17
               '  V_total_earth', $ ;18
               '    V_par_earth', $ ;19
               '   V_perp_earth', $ ;20
               '  T_total_earth', $ ;21
               '  P_total_earth', $ ;22
               '      GSM_X(Re)', $ ;23
               '      GSM_Y(Re)', $ ;24
               '      GSM_Z(Re)', $ ;25
               '     MAG_X(GSE)', $ ;26
               '     MAG_Y(GSE)', $ ;27
               '     MAG_Z(GSE)', $ ;28
               '      H_DENSITY', $ ;29
               '       H_V_X   ', $ ;30
               '       H_V_Y   ', $ ;31
               '       H_V_Z   ', $ ;32
               '       H_T_X   ', $ ;33
               '       H_T_Y   ', $ ;34
               '       H_T_Z   ', $ ;35
               '      T_x_tail ', $ ;36
               '      T_y_tail ', $ ;37
               '      T_z_tail ', $ ;38
               '      T_x_earth', $ ;39
               '      T_y_earth', $ ;40
               '      T_z_earth', $ ;41
               '    Storm_Phase', $ ;42
               '      IMF_Bx   ', $ ;43
               '      IMF_By   ', $ ;44
               '      IMF_Bz   ', $ ;45
               '       SW_V    ', $ ;46
               '       SW_P    ', $ ;47
               '    eflux_tail ', $ ;48
               '    eflux_earth', $ ;49
               '    theta_tail ', $ ;50
               '    theta_earth', $ ;51
               '         MLT   ', $ ;52
               '       ILAT_D  ', $ ;53
               '    IMF_Bx_DC  ', $ ;54
               '    IMF_By_DC  ', $ ;55
               '    IMF_Bz_DC  ', $ ;56
               '      SW_V_DC  ', $ ;57
               '      SW_P_DC  ', $ ;58
               ' DistFunc_tail ', $ ;59
               ' DistFunc_earth'] ;, $ ;60
;               '  Vgse_tail_x  ', $ ;61
;               '  Vgse_tail_y  ', $ ;62
;               '  Vgse_tail_z  ', $ ;63
;               '  Vgse_earth_x ', $ ;64
;               '  Vgse_earth_y ', $ ;65
;               '  Vgse_earth_z '] ;66
; Setup data arrays 
    title_lenth = STRLEN(title)
    nterms =  N_ELEMENTS(title)
    ntime = ndays*24*12
    time = DBLARR(ntime)
    data = DBLARR(ntime, nterms)
;---------------------------------------------------------------
; Read data into time and data arraies
;---------------------------------------------------------------
; Read all 1 day files that correspond to requested time interval
    jj = 0l                  
    FOR iday = 0l, ndays-1 DO BEGIN ; Loop trough all days   
        caldat, jd_s + iday, month, day, year ; find caledar date
        month_str = string(month, format = '(i2.2)')
        day_str = string(day, format = '(i2.2)')
        year_str = string(year, format = '(i4.4)')
        fln = path+'data/'+'storm_o_beam_'+year_str+month_str+day_str  +'.dat'
        
        names = FINDFILE(fln)
        IF names(0) NE '' THEN BEGIN 
            OPENR, unit, names(0), /GET_LUN
            dummy = ''
            a = DBLARR(nterms)
            WHILE NOT EOF(unit) DO BEGIN
                READF, unit, dummy
                IF STRMID(dummy, 0, 3) EQ  '200' THEN BEGIN
                    time(jj) = time_double(STRMID(dummy, 0, 10) + $
                                           '/'+STRMID(dummy, 11, 12))
                    IF time(jj) GE ts AND time(jj) LE te THEN BEGIN  
                        READS, STRMID(dummy, 31), a  
                        data(jj, *) = a
                        jj = jj+1
                    ENDIF 
                ENDIF   
            ENDWHILE
            CLOSE, unit, /all
        ENDIF   
    ENDFOR 
    ntime = jj
    IF ntime EQ 0 THEN BEGIN 
        print, 'no files found' & stop
    ENDIF 
    time = time(0:ntime-1)
    data = data(0:ntime-1, *)
    
    IF KEYWORD_SET(save_data) THEN BEGIN 
        store_data, 'data', data = {x:time, y:data, title:title}
        tplot_save, 'data', filename = flndata
    ENDIF 
ENDELSE 
; when data stored, flux and eflux data were divided by a factor so
; most of the information can be written. Here we put the number back
data(*, 7) = data(*, 7)*10
data(*, 16) = data(*, 16)*10
data(*, 48) = data(*, 48)*1.e3
data(*, 49) = data(*, 49)*1.e3
data(*,6) = 90-ABS(data(*,6)-90)
data(*,15) = 90-ABS(data(*,15)-90)
data(*,59) = data(*,59)/1e6
data(*,60) = data(*,60)/1e6
;--------------------------------------------------------
;Elimilate data during in the magnetosheath, solar wind and uncertain time intervals
;--------------------------------------------------------
h_density = data(*, 29)
h_velocity = sqrt((data(*, 30))^2+(data(*, 31))^2 + (data(*, 32))^2)
x_gse = data(*, 2)
y_gse = data(*, 3)
z_gse = data(*,4)
z_gsm = data(*,25)
plasma_beta = abs(data(*,1))
;- elimilate magnetosheath data by set the flag:data(*,0) into INF 
IF KEYWORD_SET(no_magnetosheath) THEN begin 
    locate_sheath = float((h_density GT 3 and h_velocity GT 65) or (x_gse gt 1 and ABS(z_gsm) gt 5 and plasma_beta gt 0.05) or (x_gse le 1 and ABS(z_gsm) gt 10 and plasma_beta gt 1))
    data(*, 0) = data(*, 0)/abs(1-locate_sheath)
endif 

;- elimilate solarwind data: define solarwind region as 
;x_gse > = -1 and outside an ellips deceided by observations 
IF KEYWORD_SET(no_solarwind) THEN  data(*, 0) = data(*, 0)/((x_gse LE -1.) OR (((x_gse+1.)^2/8.^2+(y_gse-1.)^2/14.^2) LE 1))

;- elimilate the data durig uncertain time intervals
;  read the uncertain time interval file
IF KEYWORD_SET(delete_uncertain_time) THEN BEGIN
    OPENR, unit, 'uncertain_time_interval.dat', /GET_LUN
    utime_start = DBLARR(3000)
    utime_end = DBLARR(3000)
    jj = 0l
    dummy = ''
    WHILE NOT EOF(unit) DO BEGIN
        READF, unit, dummy
        IF STRMID(dummy, 0, 5) NE 'Start' THEN BEGIN  
            utime_start(jj) = time_double(STRMID(dummy, 0, 20))
            utime_end(jj) = time_double(STRMID(dummy, 20, 20))
        ENDIF 
        jj = jj+1      
    ENDWHILE
    CLOSE, unit, /all 
    utime_start = utime_start(0:jj-1)
    utime_end = utime_end(0:jj-1)
;  clean up the uncertain time data
    FOR iut = 0, jj-1 DO data(*, 0) = $ 
      data(*, 0)/(time LT utime_start(iut) OR time GT utime_end(iut))
ENDIF 

IF KEYWORD_SET(delete_HIA_wrong_time) AND inst_name EQ 'HIA' THEN BEGIN 
    OPENR, unit, 'HIA_wrong_time.dat', /GET_LUN
    utime_start = DBLARR(3000)
    utime_end = DBLARR(3000) 
    jj = 0l
    dummy = ''
    WHILE NOT EOF(unit) DO BEGIN
        READF, unit, dummy
        IF STRMID(dummy, 0, 5) NE 'Start' THEN BEGIN  
            utime_start(jj) = time_double(STRMID(dummy, 0, 20))
            utime_end(jj) = time_double(STRMID(dummy, 20, 20))  
        ENDIF 
        jj = jj+1      
    ENDWHILE
    CLOSE, unit, /all 
    utime_start = utime_start(0:jj-1)
    utime_end = utime_end(0:jj-1)
;  clean up the uncertain time data
    FOR iut = 0, jj-1 DO data(*, 0) = $ 
      data(*, 0)/(time LT utime_start(iut) OR time GT utime_end(iut))
ENDIF 

;- set all infinite flag value to nan 
index = WHERE( ~FINITE(data(*, 0)), ct)
IF ct GT 0 THEN data(index, 0) = !VALUES.F_NAN

;---------------------------------------------
; save data into array with meanful names
;--------------------------------------------
flag = data(*, 0)
beta = ABS(data(*, 1))
storm_phase = data(*, 42) 
imf_bx = data(*, 43)
imf_by = data(*, 44)
imf_bz = data(*, 45)
imf_b = sqrt(data(*, 43)^2+data(*, 44)^2+data(*, 45)^2)
sw_v =  data(*, 46)
sw_p =  data(*, 47)
bx = data(*, 26)
By = data(*, 27)
Bz = data(*, 28)
b = sqrt(data(*, 26)^2+data(*, 27)^2+data(*, 28)^2)
y_gsm= data(*,24)
z_gsm = data(*, 25)
anodes_tail = (data(*, 50)+101.25)/22.5
anodes_earth = (data(*, 51)+101.25)/22.5
energy_tail = data(*, 5)
energy_earth = data(*, 14)
pa_tail = data(*, 6)
pa_earth = data(*, 15)
flux_tail = data(*, 7)
flux_earth = data(*, 16)
density_tail = data(*, 8)
density_earth = data(*, 17)
velocity_tail = data(*, 9)
velocity_earth = data(*, 18)
v_para_tail = data(*, 10)
v_perp_tail=data(*,11)
v_para_earth = data(*, 19)
v_perp_earth=data(*,20)
temperature_tail = data(*, 12)
temperature_earth = data(*, 21)
pressure_tail=data(*,13)
pressure_earth=data(*,22)
eflux_tail = data(*, 48)
eflux_earth = data(*, 49)
v_energy_tail = velocity_tail^2*mass_o/2
v_energy_earth = velocity_earth^2*mass_o/2
proton_n=data(*,29)
proton_Vx=data(*,30)
proton_Vy=data(*,31)
proton_Vz=data(*,32)
proton_v=sqrt(proton_vx^2+proton_vy^2+proton_vz^2)
proton_Tx=data(*,33)
proton_Ty=data(*,34)
proton_Tz=data(*,35)
proton_T=sqrt(proton_Tx+proton_Ty+proton_Tz)

VxBz=sw_v*(imf_bz<0)*(1e3*1e-9*1e3) ;calcualted SW Ey. sw velocity are considered to be as sw Vx here and only keep southward Bz since dayside reconnection brings the Ey in. unit is mV/m


;Vx_tail=data(*,61)  & Vy_tail=data(*,62)  & Vz_tail=data(*,63)
;Vx_earth=data(*,64) & Vy_earth=data(*,65) & Vz_earth=data(*,66)
;trying to adding oxygen V into bulk V 
;prepare data for calculation
;index=where(~finite(Vx_tail),ct)
;if ct gt 0 then begin 
;    Vx_tail(index)=0 & Vy_tail(index)=0 & Vz_tail(index)=0
;endif 
;index=where(~finite(Vx_earth),ct)
;if ct gt 0 then begin 
;    Vx_earth(index)=0 & Vy_earth(index)=0 & Vz_earth(index)=0
;endif 
;Vbulk_x=(proton_Vx*proton_n>0)+(Vx_tail*density_tail>0)+(Vx_earth*density_earth>0)/((proton_n>0)+(density_tail>0)+(density_earth>0))
;Ex = -(Vbulk_y * Bz - Vbulk_z * By) * 1e3 * 1e-9 * 1e3 ;mV/m  E=VXB
;Ey = -(Vbulk_z * Bx - Vbulk_x * Bz) * 1e3 * 1e-9 * 1e3
;Ez = -(Vbluk_x * By - Vbulk_y * Bx) * 1e3 * 1e-9 * 1e3
if not keyword_set(use_proton) then begin
    E_t=(V_perp_tail * B) * 1e3 * 1e-9 * 1e3  ;mV/m  E=-VXB
    EXB_t=E_t*B/B^2 * 1e-3*1e9*1e-3 ;km/s V=EXB
endif else begin
    Ex = -(proton_Vy * Bz - Proton_Vz * By) * 1e3 * 1e-9 * 1e3 ;mV/m  E=VXB
    Ey = -(proton_Vz * Bx - Proton_Vx * Bz) * 1e3 * 1e-9 * 1e3
    Ez = -(proton_Vx * By - Proton_Vy * Bx) * 1e3 * 1e-9 * 1e3 

    E = [[Ex],[Ey],[Ez]]
    e_t=sqrt(ex^2+ey^2+ez^2)
    EXB_x= -(Ey * Bz - Ez * By)/b^2*(1e3*1e9*1e-3*1e-3) ;km/s drift velocity EXB
    EXB_y= -(Ez * Bx - Ex * Bz)/b^2*(1e3*1e9*1e-3*1e-3)
    EXB_z= -(Ex * By - Ey * Bx)/b^2*(1e3*1e9*1e-3*1e-3)
    EXB=[[EXB_X],[EXB_Y],[EXB_Z]]
    exb=sqrt(exb_x^2+exb_y^2+exb_z^2)
endelse 
;reload V vector for further calculation
;Vx_tail=data(*,61)  & Vy_tail=data(*,62)  & Vz_tail=data(*,63)
;Vx_earth=data(*,64) & Vy_earth=data(*,65) & Vz_earth=data(*,66)
;V_paraE_tail=(Vx_tail*Ex+Vy_tail*Ey+Vz_tail*Ez)/sqrt(Ex^2+Ey^2+Ez^2)
;V_paraE_earth=(Vx_earth*Ex+Vy_earth*Ey+Vz_earth*Ez)/sqrt(Ex^2+Ey^2+Ez^2)
;V_paraE_proton=(proton_Vx*Ex+proton_Vy*Ey+proton_Vz*Ez)/sqrt(Ex^2+Ey^2+Ez^2)
dist=data(*,23)^2+data(*,24)^2+data(*,25)^2
;index=where(sqrt(x_gse^2+y_gse^2+z_gse^2) gt 5,ct)
;if ct gt 0 then data(index,0) = !values.F_nan

IF keyword_set(stop) THEN stop

If keyword_set(short_time_plot) then begin 
    store_data,'imf_bz',data={x:time,y:imf_bz}
    store_data,'storm_phase',data={x:time,y:storm_phase}
endif 

;en vs plasma beta plot
IF keyword_set(en_vs_beta) THEN BEGIN 
    FOR ip = 0, 1 DO BEGIN 
        IF ip EQ 0 THEN BEGIN 
            ene = data(*, 14)*(x_gse LE -15)*(storm_phase EQ 0 or storm_phase eq 5)
            ent = data(*, 5)*(x_gse LE -15)*(storm_phase EQ 0 or storm_phase eq 5)
            phase_beta = 'nonstorm'
        ENDIF ELSE BEGIN 
            ene = data(*, 14)*(x_gse LE -15)*(storm_phase GT  0 and storm_phase lt 4)
            ent = data(*, 5)*(x_gse LE -15)*(storm_phase GT  0 and storm_phase lt 4)
            phase_beta = 'storm'
        ENDELSE 
        beta_scale_plot = [1/10^2.5, 1/10^1.5, 1/10^0.5, 10^0.5, 10^1.5]
        beta_scale = [0.01, 0.1, 1, 10]
        en_scale = [ 40000., 31444., 19398.,  11966., 7382., $
                     4554.,  2809.,  1733., 1069., $
                     659, 406., 251., 154.,  95.,  58., 36.]
        en_scale_plot = ((en_scale(0:14)+en_scale(1:15)))/2
        
        counts = dblarr(5, 15)
        FOR ien = 0, 14 DO BEGIN
            FOR ibt = 0, 2 DO BEGIN 
                counts(ibt+1, ien) = total(beta ge beta_scale(ibt) AND beta lT beta_scale(ibt+1)AND $
                                           ent lt en_scale(ien) AND ent ge en_scale(ien+1))+ $
                  total(beta gE beta_scale(ibt) AND beta LT beta_scale(ibt+1) AND $
                        ene lt en_scale(ien) AND ene ge en_scale(ien+1))
            ENDFOR
            counts(0, ien) = total(beta LE beta_scale(0) AND $
                                   ent lt en_scale(ien) AND ent ge en_scale(ien+1))+ $
              total(beta le beta_scale(0) AND $
                    ene lt en_scale(ien) AND ene ge en_scale(ien+1))
            counts(4, ien) =  total(beta gt beta_scale(3) AND $
                                    ent lt en_scale(ien) AND ent ge en_scale(ien+1))+ $
              total(beta Gt beta_scale(3) AND $
                    ene lt en_scale(ien) AND ene ge en_scale(ien+1))
        ENDFOR 
        
        store_data, 'plot', data = {x:beta_scale, y:en_scale, v:counts}
        t_counts = dblarr(5, 15)
        
        FOR i = 0, 14 DO t_counts(*, i) = total(counts, 2)
        n_counts = counts/(t_counts)
        popen, path+'plots/beta_energy/energy_vs_beta_'+phase_beta+'.ps', /land
        specplot, beta_scale_plot, en_scale_plot, counts, $
          no_interp = 1, $
          lim = { title:'O+ Beam Energy vs beta '+phase_beta+' (x<-15Re)', $
                  xtitle: 'beta', $
                  ytitle: 'Energy', $
                  ztitle: 'Normalized Beam Numbers (log)', $
                  xrange: [0.001, 100], xlog:1, $
                  yrange: [40, 40000], ylog:1, $
                  zrange: [1, 1000], zlog: 1, $
                  XSTYLE:1, ystyle: 1, zstyle:1, $
                  position:[0.2, 0.2, 0.8, 0.8], charsize: 1.5}
        pclose
        popen, path+'plots/beta_energy/energy_vs_beta_normalized_'+phase_beta+'.ps', /land
        specplot, beta_scale_plot, en_scale_plot, n_counts, $
          no_interp = 1, $
          lim = { title:'O+ Beam Energy vs beta (Normalized)'+phase_beta+' (x<-15Re)', $
                  xtitle: 'beta', $
                  ytitle: 'Energy', $
                  ztitle: 'Normalized Beam Numbers (log)', $
                  xrange: [0.001, 100], xlog:1, $
                  yrange: [40, 40000], ylog:1, $
                  zrange: [0.01, 1], zlog: 1, $
                  XSTYLE:1, ystyle: 1, zstyle:1, $
                  position:[0.2, 0.2, 0.8, 0.8], charsize: 1.5}
        pclose
    ENDFOR 
ENDIF 
;--------------------------------------
;set the path for plots
;----------------------------------------
main_path = path+'plots/'  &  spawn, 'mkdir ' + main_path
main_path = main_path+'map/' &  spawn, 'mkdir ' + main_path

IF keyword_set(symmetry_template) THEN  main_path = main_path+'symmetry/' &  spawn, 'mkdir ' + main_path
;------------------------------
; non_sort plots 
;----------------------------------
IF KEYWORD_SET(non_sort) THEN BEGIN
    IF keyword_set(map) THEN BEGIN 
        plot_path = main_path +'non_sort/' & spawn, 'mkdir ' + plot_path
        
        sort_flag = 1
                                ;      sort_flag= 1./(data(*,5) gt 0. and data(*,5) lt 40. and beta lt 0.02 and data(*,2) lt -5.)
        
        sort_title = ''
        map_o_beam, time, data, $
          sort_flag = sort_flag, sort_title = sort_title, $
          sc = sc, $
          ps_plot = ps_plot, $
          plot_path = plot_path, $
          coor_set = coor_set, $
          grid_set = grid_set, slice_grid_set = slice_grid_set, $
          direction_set = direction_set, $
          storm_phase_set = storm_phase_set, $
          plot_2d = plot_2d, slice_plot = slice_plot, $
          waterdrop_plot = waterdrop_plot, $
          point_plot = point_plot, $
          events_map = events_map,  $
          property_map_set = property_map_set, $
          diff_beta = diff_beta, $
          symmetry_template = symmetry_template, $
          property_map_type_set =  property_map_type_set
    ENDIF 
ENDIF
;------------------------------
;make table of ratio in big regions
;----------------------------------
IF KEYWORD_SET(make_table) THEN BEGIN
    plot_path = main_path +'table/' & spawn, 'mkdir ' + plot_path
    sort_flag = 1
    sort_title = ''
    map_o_beam, time, data, $
      sort_flag = sort_flag, sort_title = sort_title, $
      sc = sc, $
      ps_plot = ps_plot, $
      plot_path = plot_path, $
      coor_set = ['X_GSM', 'Z_GSM'] , $
      grid_set = grid_set, slice_grid_set = slice_grid_set, $
      direction_set = direction_set, $
      storm_phase_set = ['nonstorm_time','prestorm','initial_phase','main_phase','recovery'], $
      plot_2d = 0, slice_plot = 0, $
      waterdrop_plot = 0, $
      point_plot = 0, $
      events_map = 1, $
      property_map_set = property_map_set, $
      diff_beta = [1], $
      symmetry_template = symmetry_template, $
      property_map_type_set =  property_map_type_set,$
      make_table=make_table
    
ENDIF 

;------------------------------
; bx divide the plot instead of z gsm (useful for dividing plamsasheet)
;----------------------------------
IF KEYWORD_SET(sort_bx) THEN BEGIN
    IF keyword_set(map) THEN BEGIN 
        plot_path = main_path +'sort_bx/' & spawn, 'mkdir ' + plot_path
        
        FOR ibx = 0, 1 DO BEGIN 
            IF ibx EQ 0 THEN  BEGIN 
                sort_path = 'bx_ge_0/'
                sort_flag = 1./(bx GE 0)
                sort_title = 'bx > 0'
            ENDIF ELSE BEGIN 
                sort_path = 'bx_lt_0/'
                sort_flag = 1./(bx LT 0)
                sort_title = 'bx < 0'
            ENDELSE 
            plot_path_bx = plot_path+sort_path & spawn, 'mkdir ' + plot_path_bx
            
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = plot_path_bx, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              symmetry_template = symmetry_template, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF  
ENDIF
;-----------------------------------------
;sort by different conditions 
;------------------------------------------
; restore data of index and solar wind & IMF pre 1 hour data
;------------------------------------------------------------
flndata_sort = path+'tplot_restore/sort_condition_data_'+time_str
print, FINDFILE(flndata_sort+'.tplot', COUNT = ct_sort_store)
IF ct_sort_store GT 0 THEN tplot_restore, filenames = flndata_sort+'.tplot' 

;-------------------------------------
; plot within energy range threshold
;-------------------------------------
IF KEYWORD_SET(energy_filter) THEN BEGIN
    if keyword_set(eflux_filter) then begin 
        eflux_threshold_str = STRING(eflux_threshold, format = '(i4.4)')
        eflux_tail=data(*,48)
        eflux_earth=data(*,49)
    endif 
    en_set = [31444.7, 19398.3, 11966.9, 7382.39, 4554.22, 2809.51, 1733.19, 1069.21, $
              659.599, 406.909, 251.023, 154.857, 95.5315, 58.9337, 36.3563]
    den_set = [14900.6, 9192.19,5670.69,3498.26,2158.09,1331.33,821.301,506.663,$
               312.562, 192.820,118.951,73.3813,45.2691,27.9267,17.2280,5.63791]
;    energy_range_set = [[en_set-den_set/2],[en_set+den_set/2]]
    energy_range_set=[[20,100,300,2000],[100,300,2000.,40000]]

    for ien=0,n_elements(energy_range_set)/2-1 do begin 
        energy_range=energy_range_set(ien,*)
        energy_range_str = STRING(min(energy_range,/nan), format = '(i5.5)')+'_'+STRING(max(energy_range,/nan),format='(i5.5)')
        energy_tail=data(*,5)
        energy_earth=data(*,14)
; MAP under certain condition IF NEEDED
        IF KEYWORD_SET(MAP) THEN BEGIN 
            if keyword_set(eflux_filter) then begin 
                plot_path = main_path+'energy_filter_eflux_threshold/'  & spawn, 'mkdir ' + plot_path
                sort_path = plot_path+'eflux_ge_'+ eflux_threshold_str+'_'+'energy_'+ energy_range_str+'/' & spawn,'mkdir '+ sort_path
                sort_flag = 1./(((eflux_tail GE eflux_threshold or eflux_earth ge eflux_threshold) or flag eq 0) and (energy_tail GE min(energy_range,/nan) and energy_tail le max(energy_range,/nan) or (energy_earth ge min(energy_range,/nan) and energy_earth le max(energy_range,/nan)) or flag eq 0))
                sort_title = 'eflux > '+ eflux_threshold_str +', '+'energy range: '+ energy_range_str
            endif else begin 
                plot_path = main_path+'energy_filter/'  & spawn, 'mkdir ' + plot_path
                sort_path = plot_path+'energy_'+ energy_range_str+'/' & spawn,'mkdir '+ sort_path
                sort_flag = 1./((energy_tail GE min(energy_range,/nan) and energy_tail le max(energy_range,/nan)) or (energy_earth ge min(energy_range,/nan) and energy_earth le max(energy_range,/nan)) or flag eq 0)
                sort_title = 'energy range: '+ energy_range_str
            endelse 
            

            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = sort_path, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              symmetry_template = symmetry_template, $
              property_map_type_set =  property_map_type_set
        endif             
    endfor  
ENDIF    

;-------------------------------------
; filter with eflux threshold
;-------------------------------------
IF KEYWORD_SET(eflux_filter) THEN BEGIN
    eflux_threshold_str = STRING(eflux_threshold, format = '(i4.4)')
    eflux_tail=data(*,48)
    eflux_earth=data(*,49)
; MAP under certain condition IF NEEDED
    IF KEYWORD_SET(MAP) THEN BEGIN 
        plot_path = main_path+'eflux_threshold/'  & spawn, 'mkdir ' + plot_path
        
        sort_path = plot_path+'eflux_ge_'+ eflux_threshold_str+'/' & spawn,'mkdir '+ sort_path
        
        sort_flag = 1./((eflux_tail GE eflux_threshold or eflux_earth ge eflux_threshold) or flag eq 0)
        sort_title = 'eflux > '+ eflux_threshold_str
        
        map_o_beam, time, data, $
          sort_flag = sort_flag, sort_title = sort_title, $
          sc = sc, $
          ps_plot = ps_plot,  $
          plot_path = sort_path, $
          coor_set = coor_set, $
          grid_set = grid_set, slice_grid_set = slice_grid_set, $
          direction_set = direction_set, $
          storm_phase_set = storm_phase_set, $
          plot_2d = plot_2d, slice_plot = slice_plot, $
          waterdrop_plot = waterdrop_plot, $
          point_plot = point_plot, $
          events_map = events_map,  $
          property_map_set = property_map_set, $
          diff_beta = diff_beta, $
          symmetry_template = symmetry_template, $
          property_map_type_set =  property_map_type_set
    endif           
ENDIF     
;------------------------------------------
; sort by season
;------------------------------------------
IF KEYWORD_SET(sort_season) THEN BEGIN
    if keyword_set(eflux_filter) then begin 
        eflux_threshold_str = STRING(eflux_threshold, format = '(i4.4)')
        eflux_tail=data(*,48)
        eflux_earth=data(*,49)
    endif
    hemi_set=['north_lobe','south_lobe'] ;'all'  ;plot for different hemispehres
    if keyword_set(hemi_set) then nhemi=n_elements(hemi_set) else nhemi=1
    for ihemi=0,nhemi-1 do begin
        hemi_str=hemi_set(ihemi)
; MAP under certain condition IF NEEDED 
    plot_path = main_path+'sort_season/'  & spawn, 'mkdir ' + plot_path
    if keyword_set(eflux_filter) then plot_path=plot_path+'eflux_ge_'+ eflux_threshold_str+'/' & spawn,'mkdir '+ plot_path
    if keyword_set(hemi_set)  then plot_path=plot_path+hemi_str+'/' & spawn,'mkdir '+ plot_path

    IF KEYWORD_SET(MAP) THEN BEGIN
        season=['spring','summer','fall','winter','Engwall','July_Oct','Nov','July_Aug','Sep_Oct']
        for iseason = 5, 5 do begin ;n_elements(season)-1 do begin 
            case iseason of 
                0: sort_flag = 1./((time ge time_double('2001-03-20') and time lt time_double('2001-06-20')) or $
                                   (time ge time_double('2002-03-20') and time lt time_double('2002-06-20')) or $
                                   (time ge time_double('2003-03-20') and time lt time_double('2003-06-20')) or $
                                   (time ge time_double('2004-03-20') and time lt time_double('2004-06-20')) or $
                                   (time ge time_double('2005-03-20') and time lt time_double('2005-06-20')) or $
                                   (time ge time_double('2006-03-20') and time lt time_double('2006-06-20')) or $
                                   (time ge time_double('2007-03-20') and time lt time_double('2007-06-20')) or $
                                   (time ge time_double('2008-03-20') and time lt time_double('2008-06-20')) or $
                                   (time ge time_double('2009-03-20') and time lt time_double('2009-06-20')))

                1: sort_flag = 1./((time ge time_double('2001-06-20') and time lt time_double('2001-09-22')) or $
                                   (time ge time_double('2002-06-20') and time lt time_double('2002-09-22')) or $
                                   (time ge time_double('2003-06-20') and time lt time_double('2003-09-22')) or $
                                   (time ge time_double('2004-06-20') and time lt time_double('2004-09-22')) or $
                                   (time ge time_double('2005-06-20') and time lt time_double('2005-09-22')) or $
                                   (time ge time_double('2006-06-20') and time lt time_double('2006-09-22')) or $
                                   (time ge time_double('2007-06-20') and time lt time_double('2007-09-22')) or $
                                   (time ge time_double('2008-06-20') and time lt time_double('2008-09-22')) or $
                                   (time ge time_double('2009-06-20') and time lt time_double('2009-09-22')))

                2: sort_flag = 1./((time ge time_double('2001-09-22') and time lt time_double('2001-12-21')) or $
                                   (time ge time_double('2002-09-22') and time lt time_double('2002-12-21')) or $
                                   (time ge time_double('2003-09-22') and time lt time_double('2003-12-21')) or $
                                   (time ge time_double('2004-09-22') and time lt time_double('2004-12-21')) or $
                                   (time ge time_double('2005-09-22') and time lt time_double('2005-12-21')) or $
                                   (time ge time_double('2006-09-22') and time lt time_double('2006-12-21')) or $
                                   (time ge time_double('2007-09-22') and time lt time_double('2007-12-21')) or $
                                   (time ge time_double('2008-09-22') and time lt time_double('2008-12-21')) or $
                                   (time ge time_double('2009-09-22') and time lt time_double('2009-12-21')))

                3: sort_flag = 1./((time ge time_double('2001-12-21') and time lt time_double('2002-03-20')) or $
                                   (time ge time_double('2002-12-21') and time lt time_double('2003-03-20')) or $
                                   (time ge time_double('2003-12-21') and time lt time_double('2004-03-20')) or $
                                   (time ge time_double('2004-12-21') and time lt time_double('2005-03-20')) or $
                                   (time ge time_double('2005-12-21') and time lt time_double('2006-03-20')) or $
                                   (time ge time_double('2006-12-21') and time lt time_double('2007-03-20')) or $
                                   (time ge time_double('2007-12-21') and time lt time_double('2008-03-20')) or $
                                   (time ge time_double('2008-12-21') and time lt time_double('2009-03-20')) or $
                                   (time ge time_double('2009-12-21') and time lt time_double('2010-03-20')))

                4: sort_flag = 1./((time ge time_double('2001-07-03') and time lt time_double('2001-11-03')) or $
                                   (time ge time_double('2002-07-03') and time lt time_double('2002-11-03')))
                5: sort_flag = 1./((time ge time_double('2001-07-01') and time lt time_double('2001-11-01')) or $
                                   (time ge time_double('2002-07-01') and time lt time_double('2002-11-01')) or $
                                   (time ge time_double('2003-07-01') and time lt time_double('2003-11-01')) or $
                                   (time ge time_double('2004-07-01') and time lt time_double('2004-11-01')) or $
                                   (time ge time_double('2005-07-01') and time lt time_double('2005-11-01')) or $
                                   (time ge time_double('2006-07-01') and time lt time_double('2006-11-01')) or $
                                   (time ge time_double('2007-07-01') and time lt time_double('2007-11-01')) or $
                                   (time ge time_double('2008-07-01') and time lt time_double('2008-11-01')) or $
                                   (time ge time_double('2009-07-01') and time lt time_double('2009-11-01')))
                6: sort_flag = 1./((time ge time_double('2001-11-01') and time lt time_double('2001-12-01')) or $
                                   (time ge time_double('2002-11-01') and time lt time_double('2002-12-01')) or $
                                   (time ge time_double('2003-11-01') and time lt time_double('2003-12-01')) or $
                                   (time ge time_double('2004-11-01') and time lt time_double('2004-12-01')) or $
                                   (time ge time_double('2005-11-01') and time lt time_double('2005-12-01')) or $
                                   (time ge time_double('2006-11-01') and time lt time_double('2006-12-01')) or $
                                   (time ge time_double('2007-11-01') and time lt time_double('2007-12-01')) or $
                                   (time ge time_double('2008-11-01') and time lt time_double('2008-12-01')) or $
                                   (time ge time_double('2009-11-01') and time lt time_double('2009-12-01')))
                7: sort_flag = 1./((time ge time_double('2001-07-01') and time lt time_double('2001-08-01')) or $
                                   (time ge time_double('2002-07-01') and time lt time_double('2002-08-01')) or $
                                   (time ge time_double('2003-07-01') and time lt time_double('2003-08-01')) or $
                                   (time ge time_double('2004-07-01') and time lt time_double('2004-08-01')) or $
                                   (time ge time_double('2005-07-01') and time lt time_double('2005-08-01')) or $
                                   (time ge time_double('2006-07-01') and time lt time_double('2006-08-01')) or $
                                   (time ge time_double('2007-07-01') and time lt time_double('2007-08-01')) or $
                                   (time ge time_double('2008-07-01') and time lt time_double('2008-08-01')) or $
                                   (time ge time_double('2009-07-01') and time lt time_double('2009-08-01')))
                8: sort_flag = 1./((time ge time_double('2001-09-01') and time lt time_double('2001-11-01')) or $
                                   (time ge time_double('2002-09-01') and time lt time_double('2002-11-01')) or $
                                   (time ge time_double('2003-09-01') and time lt time_double('2003-11-01')) or $
                                   (time ge time_double('2004-09-01') and time lt time_double('2004-11-01')) or $
                                   (time ge time_double('2005-09-01') and time lt time_double('2005-11-01')) or $
                                   (time ge time_double('2006-09-01') and time lt time_double('2006-11-01')) or $
                                   (time ge time_double('2007-09-01') and time lt time_double('2007-11-01')) or $
                                   (time ge time_double('2008-09-01') and time lt time_double('2008-11-01')) or $
                                   (time ge time_double('2009-09-01') and time lt time_double('2009-11-01')))
            endcase

            if keyword_set(eflux_filter) then sort_flag= 1./((1./((eflux_tail GE eflux_threshold or eflux_earth ge eflux_threshold) or flag eq 0)) and sort_flag eq 1)
            if keyword_set(hemi_set) then begin
                if hemi_str eq 'north_lobe' then sort_flag=1./((z_gsm gt 0 and beta le 0.05) and sort_flag eq 1)
                if hemi_str eq 'south_lobe' then sort_flag=1./((z_gsm lt 0 and beta le 0.05) and sort_flag eq 1)
            endif 
            
            sort_title = season(iseason)
            sort_path = sort_title+'/'
            
            if keyword_set(eflux_filter) then sort_title = sort_title+' '+'eflux > '+ eflux_threshold_str
            if keyword_set(hemi_set) then sort_title = sort_title + '  '+hemi_str
          
            plot_path_season = plot_path+sort_path & spawn, 'mkdir  ' + plot_path_season

            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot, $
              plot_path = plot_path_season, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map, $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF 
endfor  
ENDIF

;------------------------------------
; F10.7 index (for quiet time/non) 
;------------------------------------
IF KEYWORD_SET(sort_F107) THEN BEGIN
    sort_f107_str = STRING(sort_f107, format = '(i3.1)')
; Read F10.7 Index
;sort_recalc=1
    IF ct_sort_store EQ 0 OR keyword_set(sort_recalc) THEN BEGIN 
        read_omni, all = 1
        get_data, 'F10_7_Index', data = sort_data
        sort_data_y = FLOAT(sort_data.y)
        index = where(sort_data.y GE 900, ct)
        IF ct GT 0 THEN sort_data_y(index) = !VALUES.F_NAN  
        f107_index = INTERPOL(sort_data_y, sort_data.x, time)
        store_data, 'F10_7_Index', data = {x:time, y:f107_index}
    ENDIF ELSE BEGIN 
        get_data, 'F10_7_Index', data = sort_data
        f107_index = sort_data.y
    ENDELSE 

    plot_path = main_path+'f107_sort/'  & spawn, 'mkdir ' + plot_path
; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(MAP)THEN BEGIN 
        FOR if107 = 0, 2 DO BEGIN
            IF if107 EQ 0 THEN  BEGIN 
                sort_path = 'f107_ge_'+sort_f107_str+'/' 
                sort_flag = 1./(f107_index GE sort_f107)
                sort_title = 'f107 > '+sort_f107_str
            ENDIF ELSE BEGIN 
                sort_path = 'f107_lt_'+sort_f107_str+'/' 
                sort_flag = 1./(f107_index LT sort_f107)
                sort_title = 'f107 < '+sort_f107_str
            ENDELSE 

            plot_path_f107 = plot_path+sort_path & spawn, 'mkdir ' + plot_path_f107
            
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = plot_path_f107, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF 
ENDIF 

;--------------------------------------
; KP index (for quiet time/non) 
;------------------------------------------
IF KEYWORD_SET(sort_kp) THEN BEGIN
    sort_kp_str = STRING(sort_kp, format = '(i1.1)')
    storm_phase_set_kp = ['all_time', 'nonstorm_time']

; Read KP Index
    IF ct_sort_store EQ 0 OR keyword_set(sort_recalc) THEN BEGIN 
        read_omni
        get_data, 'Kp_Index', data = sort_data
        sort_data_y = FLOAT(sort_data.y)
        index = where(sort_data.y EQ 99, ct)
        IF ct GT 0 THEN sort_data_y(index) = !VALUES.F_NAN  
        kp_index = (INTERPOL(sort_data_y, sort_data.x, time))/10.
        store_data, 'KP_Index', data = {x:time, y:kp_index}
    ENDIF ELSE BEGIN 
        get_data, 'KP_Index', data = sort_data
        kp_index = sort_data.y
    ENDELSE  
    plot_path = main_path+'kp_sort/'  & spawn, 'mkdir ' + plot_path
; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(MAP)THEN BEGIN 
        FOR ikp = 0, 1 DO BEGIN 
            IF ikp EQ 0 THEN  BEGIN 
                sort_path = 'kp_ge_'+sort_kp_str+'/' 
                sort_flag = 1./(kp_index GE sort_kp)
                sort_title = 'kp > '+sort_kp_str
            ENDIF ELSE BEGIN 
                sort_path = 'kp_lt_'+sort_kp_str+'/' 
                sort_flag = 1./(kp_index LT sort_kp)
                sort_title = 'kp < '+sort_kp_str
            ENDELSE 
            plot_path_kp = plot_path+sort_path & spawn, 'mkdir ' + plot_path_kp
            
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = plot_path_kp, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set_kp, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF 
ENDIF 

if keyword_set(sort_substorm) then begin
    tplot_names,'substorm_flag',names=names
    if names eq '' or keyword_set(sort_recalc) then begin 
        expan_time=DBLARR(1000,2)
        jj = 0l                  
        fln = 'substorm_list.dat'
        names = FINDFILE(fln)
        IF names(0) NE '' THEN BEGIN 
            OPENR, unit, names(0), /GET_LUN
            dummy = ''
            a = STRARR(2)
            WHILE NOT EOF(unit) DO BEGIN
                READF, unit, dummy
                IF STRMID(dummy, 0, 3) EQ  '200' THEN BEGIN
                    expan_time(jj,0) = time_double(strmid(dummy,0,20))
                    expan_time(jj,1) = time_double(strmid(dummy,20,20))
                    jj = jj+1
                ENDIF   
            ENDWHILE
            CLOSE, unit, /all
        ENDIF  else stop 
        expan_time = expan_time(0:jj-1,*)
        substorm_flag=DBLARR(n_elements(time))
        for it=0l,n_elements(time)-1 do begin
            iexpan=0
            while iexpan le n_elements(expan_time(*,0))-1 and substorm_flag(it) eq 0 do begin              
                if time(it) gt expan_time(iexpan,0) and time(it) lt expan_time(iexpan,1) then substorm_flag(it)=1
                iexpan=iexpan+1
            endwhile  
        endfor 
        store_data,'substorm_flag',data={x:time,y:substorm_flag,expan_time:expan_time}
    endif

    get_data,'substorm_flag',data=dd
    substorm_flag = dd.y
    region_str=['lobe']
    plot_path = main_path+'substorm_sort/'  & spawn, 'mkdir ' + plot_path

    for istorm=0,n_elements(storm_phase_set)-1 do begin 
        if storm_phase_set(istorm) eq 'nonstorm_time' then storm_flag = storm_phase eq 0 or storm_phase eq 5
        if storm_phase_set(istorm) eq 'storm_time' then storm_flag = storm_phase ge 1 and storm_phase le 3
        for ire=0,n_elements(region_str)-1 do begin
            if region_str(ire) eq 'lobe' then begin
                beta_flag = beta lt 0.05 & color_input=2
            endif 
            if region_str(ire) eq 'bl' then begin 
                beta_flag = beta gt 0.05 and beta lt 1 & color_input=4
            endif 
            if region_str(ire) eq 'ps' then beta_flag = beta gt 1 

            if region_str(ire) eq 'tail_lobe' then begin
                beta_flag = beta lt 0.05 and x_gse lt -5 & color_input=2
            endif 
            if region_str(ire) eq 'polar_cap' then begin 
                beta_flag = beta lt 0.05 and x_gse gt -5 & color_input=3
            endif 
            if region_str(ire) eq 'tail_bl' then begin 
                beta_flag = beta gt 0.05 and beta lt 1 and x_gse lt -5 & color_input=1
            endif 
            if region_str(ire) eq 'tail_ps' then begin 
                beta_flag = beta gt 1 and x_gse lt -5 & color_input=4
            endif 
            ind_tail = where((flag ge 1 or flag eq 2) and beta_flag eq 1 and storm_flag eq 1,ct_tail)
            ind_earth = where((flag eq -1 or flag eq 2) and beta_flag eq 1 and storm_flag eq 1,ct_earth)
            ratio_tail=density_tail/h_density
            ratio_earth=density_earth/h_density
            if ct_tail then ratio_tail(ind_tail) = !values.f_nan
            if ct_earth then ratio_earth(ind_earth)=!values.f_nan
            store_data,'ratio_tail_'+storm_phase_set(istorm), data={x:time,y:ratio_tail},lim={ylog:1,yrange:[1e-4,1e3]},dlim={psym:3,ytitle:storm_phase_set(istorm)}
            store_data,'ratio_earth_'+storm_phase_set(istorm),data={x:time,y:ratio_earth},lim={ylog:1,yrange:[1e-4,1e3]},dlim={psym:3,ytitle:storm_phase_set(istorm)}
        endfor
    endfor 
    tplot,['ratio_tail_nonstorm_time','ratio_tail_storm_time']
    tplot_panel,v='ratio_tail_nonstorm_time',o='ratio_earth_nonstorm_time'
    tplot_panel,v='ratio_tail_storm_time',o='ratio_earth_storm_time'
    stop
; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(MAP) THEN BEGIN 
        sort_flag = substorm_flag
        sort_title = 'substorm expantion phase'
        
        plot_path_substorm = plot_path & spawn, 'mkdir ' + plot_path_substorm
        
        map_o_beam, time, data, $
          sort_flag = sort_flag, sort_title = sort_title, $
          sc = sc, $
          ps_plot = ps_plot,  $
          plot_path = plot_path_substorm, $
          coor_set = coor_set, $
          grid_set = grid_set, slice_grid_set = slice_grid_set, $
          direction_set = direction_set, $
          storm_phase_set = storm_phase_set, $
          plot_2d = plot_2d, slice_plot = slice_plot, $
          waterdrop_plot = waterdrop_plot, $
          point_plot = point_plot, $
          events_map = events_map,  $
          property_map_set = property_map_set, $
          diff_beta = diff_beta, $
          property_map_type_set =  property_map_type_set
    ENDIF 
endif 
;------------------------------------------
; AE index (for substorms/non) 
;------------------------------------------
IF KEYWORD_SET(sort_ae) THEN BEGIN
    sort_ae_str = STRING(sort_ae, format = '(i3.3)')
    storm_phase_set_ae = ['nonstorm_time', 'all_time']
;Load AE Index if not restored
    IF ct_sort_store EQ 0  OR keyword_set(sort_recalc)THEN BEGIN 
        plot_ae, AE = 1, AL = 0, AO = 0, AU = 0
        get_data, 'AUX_ae', data = sort_data

        sort_data_y = FLOAT(sort_data.y)
        index = where(sort_data.y GE 8000, ct)
        IF ct GT 0 THEN sort_data_y(index) = !VALUES.F_NAN  
;average AE Index for 5 minutes around the time data 
        ae_index = dblARR(ntime)
        FOR itime = 0l, ntime-1 DO BEGIN 
            index = where(sort_data.x GE time(itime)-average_time/2. AND $
                          sort_data.x LT time(itime)+average_time/2., ct)
            IF ct GT 0 THEN BEGIN 
                IF TOTAL(ABS(sort_data_y(index)) GE 0) GT 0  THEN $
                  ae_index(itime) = total(sort_data_y(index), /NAN)/ct $
                ELSE ae_index(itime) = !VALUES.F_NAN
            ENDIF ELSE  ae_index(itime) = !VALUES.F_NAN
            IF (itime/1000.-itime/1000) EQ 0 THEN print,  itime, ntime
        ENDFOR 
        store_data, 'AE_Index', data = {x:time, y:ae_index}
    ENDIF ELSE BEGIN 
        get_data,  'AE_Index', data = sort_data
        ae_index = sort_data.y
    ENDELSE 
    plot_path = main_path+'ae_sort/'  & spawn, 'mkdir ' + plot_path
; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(MAP) THEN BEGIN 
        plot_path = main_path+'ae_sort/'  & spawn, 'mkdir ' + plot_path
        FOR iae = 0, 1 DO BEGIN 
            IF iae EQ 0 THEN  BEGIN 
                sort_path = 'ae_ge_'+sort_ae_str+'/' 
                sort_flag = 1./(ae_index GE sort_ae)
                sort_title = 'ae > '+sort_ae_str
            ENDIF ELSE BEGIN 
                sort_path = 'ae_lt_'+sort_ae_str+'/' 
                sort_flag = 1./(ae_index LT sort_ae)
                sort_title = 'ae < '+sort_ae_str
            ENDELSE 
            plot_path_ae = plot_path+sort_path & spawn, 'mkdir ' + plot_path_ae
            
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = plot_path_ae, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set_ae, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF 
ENDIF 

;-------------------------------------
; sort by IMF By
;-------------------------------------
IF KEYWORD_SET(sort_IMF_By) THEN BEGIN
    sort_imf_by_set=sort_imf_by
    for i_sort_imf_by=0,n_elements(sort_imf_by_set)-1 do begin 
        sort_imf_by= sort_imf_by_set(i_sort_imf_by)
        
        if sort_IMF_By lt 0 then sort_IMF_By =0
        sort_IMF_By_str = STRING(sort_IMF_By, format = '(i1.1)')
;average solar wind By over 1 hour before each time point
        if keyword_set(delay_to_cluster) then  begin 
            imf_by_pre=data(*,55)  
            delay_location='Delay_to_CLUSTER'
        endif else begin 
            IF keyword_set(pre1h) THEN BEGIN 
                IF ct_sort_store EQ 0  OR keyword_set(sort_recalc) THEN  imf_By_pre = predata_average(imf_by, time) $
                ELSE BEGIN 
                    get_data, 'By', data = sort_data
                    imf_by_pre = sort_data.y
                ENDELSE 
                delay_location=''
            ENDIF ELSE begin 
                imf_by_pre = data(*,44)
                delay_location='nodelay'
            endelse 
        endelse 

; MAP UNDER DIFFERENT CONDITION IF NEEDED
        IF KEYWORD_SET(MAP) THEN BEGIN 
            plot_path = main_path+'IMF_By_'+delay_location+'/'  & spawn, 'mkdir ' + plot_path
            FOR iby = 0, 1 DO BEGIN 
                IF iby EQ 0 THEN  BEGIN 
                    sort_path = 'IMF_By_ge_'+ sort_IMF_By_str+'/' 
                    sort_flag = 1./(imf_by_pre GE sort_IMF_By)
                    sort_title = 'IMF By > '+ sort_IMF_By_str+'  '+delay_location
                    
                ENDIF ELSE BEGIN 
                    sort_path = 'IMF_By_lt_neg'+ sort_IMF_By_str+'/' 
                    sort_flag = 1./(imf_by_pre LT -sort_IMF_By)
                    sort_title = 'IMF By < -'+ sort_IMF_By_str+'  '+delay_location
                ENDELSE 

                plot_path_imf_by = plot_path+sort_path & spawn, 'mkdir ' + plot_path_imf_by
                map_o_beam, time, data, $
                  sort_flag = sort_flag, sort_title = sort_title, $
                  sc = sc, $
                  ps_plot = ps_plot,  $
                  plot_path = plot_path_imf_by, $
                  coor_set = coor_set, $
                  grid_set = grid_set, slice_grid_set = slice_grid_set, $
                  direction_set = direction_set, $
                  storm_phase_set = storm_phase_set, $
                  plot_2d = plot_2d, slice_plot = slice_plot, $
                  waterdrop_plot = waterdrop_plot, $
                  point_plot = point_plot, $
                  events_map = events_map,  $
                  property_map_set = property_map_set, $
                  diff_beta = diff_beta, $
                  symmetry_template = symmetry_template, $
                  property_map_type_set =  property_map_type_set
            ENDFOR  
        ENDIF  
    endfor    
ENDIF   

;-------------------------------------
; sort by IMF By with eflux filter
;-------------------------------------
IF KEYWORD_SET(sort_IMF_By_with_eflux_filter) THEN BEGIN
    eflux_tail=data(*,48)
    eflux_earth=data(*,49)
    
    sort_imf_by_set=sort_imf_by_with_eflux_filter 
    for i_sort_imf_by=0,n_elements(sort_imf_by_set)-1 do begin 
        sort_imf_by= sort_imf_by_set(i_sort_imf_by)
        
        if sort_IMF_By lt 0 then sort_IMF_By =0
        sort_IMF_By_str = STRING(sort_IMF_By, format = '(i1.1)')
;average solar wind By over 1 hour before each time point
        if keyword_set(delay_to_cluster) then  begin 
            imf_by_pre=data(*,55)  
            delay_location='Delay_to_CLUSTER'
        endif else begin 
            IF keyword_set(pre1h) THEN BEGIN 
                IF ct_sort_store EQ 0  OR keyword_set(sort_recalc) THEN  imf_By_pre = predata_average(imf_by, time) $
                ELSE BEGIN 
                    get_data, 'By', data = sort_data
                    imf_by_pre = sort_data.y
                ENDELSE 
                delay_location=''
            ENDIF ELSE begin 
                imf_by_pre = data(*,44)
                delay_location='nodelay'
            endelse 
        endelse 

; MAP UNDER DIFFERENT CONDITION IF NEEDED
        IF KEYWORD_SET(MAP) THEN BEGIN 
            plot_path = main_path+'IMF_By_'+delay_location+'_with_eflux_threshold/'  & spawn, 'mkdir ' + plot_path
            plot_path = plot_path+'eflux_gt_'+eflux_threshold_str+'/'  & spawn, 'mkdir ' + plot_path
            FOR iby = 0, 1 DO BEGIN 
                IF iby EQ 0 THEN  BEGIN 
                    sort_path = 'IMF_By_ge_'+ sort_IMF_By_str+'/' 
                    sort_flag = 1./(imf_by_pre GE sort_IMF_By $
                                    and ((eflux_tail GE eflux_threshold or eflux_earth ge eflux_threshold) or flag eq 0) $
                                   )
                    sort_title ='eflux > '+ eflux_threshold_str+ ', IMF By > '+ sort_IMF_By_str+'  '+delay_location
                    
                ENDIF ELSE BEGIN 
                    sort_path = 'IMF_By_lt_neg'+ sort_IMF_By_str+'/' 
                    sort_flag = 1./(imf_by_pre LT -sort_IMF_By $
                                    and ((eflux_tail GE eflux_threshold or eflux_earth ge eflux_threshold) or flag eq 0) $
                                   )
                    sort_title = 'eflux > '+ eflux_threshold_str+', IMF By < -'+ sort_IMF_By_str+'  '+delay_location
                ENDELSE 

                plot_path_imf_by = plot_path+sort_path & spawn, 'mkdir ' + plot_path_imf_by
                map_o_beam, time, data, $
                  sort_flag = sort_flag, sort_title = sort_title, $
                  sc = sc, $
                  ps_plot = ps_plot,  $
                  plot_path = plot_path_imf_by, $
                  coor_set = coor_set, $
                  grid_set = grid_set, slice_grid_set = slice_grid_set, $
                  direction_set = direction_set, $
                  storm_phase_set = storm_phase_set, $
                  plot_2d = plot_2d, slice_plot = slice_plot, $
                  waterdrop_plot = waterdrop_plot, $
                  point_plot = point_plot, $
                  events_map = events_map,  $
                  property_map_set = property_map_set, $
                  diff_beta = diff_beta, $
                  symmetry_template = symmetry_template, $
                  property_map_type_set =  property_map_type_set
            ENDFOR  
        ENDIF  
    endfor     
ENDIF   

;------------------------------------
; sort by IMF Bz
;----------------------------------
IF KEYWORD_SET(sort_IMF_Bz) THEN BEGIN 
    if sort_imf_bz lt 0 then sort_imf_bz=0
    sort_IMF_Bz_str = STRING(sort_IMF_Bz, format = '(i1.1)')
;pre the solar wind data
    if keyword_set(pre1h) then begin 
        IF ct_sort_store EQ 0 OR keyword_set(sort_recalc) THEN $
          imf_Bz_pre = predata_average(imf_bz, time) $
        ELSE BEGIN 
            get_data, 'Bz', data = sort_data
            imf_bz_pre = sort_data.y
        ENDELSE  
    endif else imf_bz_pre=data(*,45)

; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(MAP) THEN BEGIN
        plot_path = main_path+'IMF_Bz/'  & spawn, 'mkdir ' + plot_path
        FOR ikp = 0, 1 DO BEGIN
            IF ikp EQ 0 THEN  BEGIN
                sort_path = 'IMF_Bz_ge_'+sort_IMF_Bz_str+'/' 
                sort_flag = 1./(imf_bz_pre GE sort_IMF_Bz)
                sort_title = 'IMF Bz > '+sort_IMF_Bz_str+'  '
            ENDIF ELSE BEGIN 
                sort_path = 'IMF_Bz_lt_neg'+sort_IMF_Bz_str+'/' 
                sort_flag = 1./(imf_bz_pre LT -sort_IMF_Bz)
                sort_title = 'IMF Bz < -'+sort_IMF_Bz_str+'  '
            ENDELSE 
            plot_path_imf_bz = plot_path+sort_path & spawn, 'mkdir ' + plot_path_imf_bz
            
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = plot_path_imf_bz, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              symmetry_template = symmetry_template, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF
    if sort_imf_bz eq 0 then sort_imf_bz =-1    
ENDIF  

;------------------------------------
; sort by long time IMF Bz
;----------------------------------
IF KEYWORD_SET(sort_long_IMF_Bz) THEN BEGIN 
    if sort_long_imf_bz lt 0 then sort_long_imf_bz=0
    sort_long_IMF_Bz_str = STRING(sort_long_IMF_Bz, format = '(i1.1)')
;pre the solar wind data
    if keyword_set(pre1h) then begin 
        IF ct_sort_store EQ 0 OR keyword_set(sort_recalc) THEN $
          imf_Bz_pre = predata_average(imf_bz, time) $
        ELSE BEGIN 
            get_data, 'Bz', data = sort_data
            imf_bz_pre = sort_data.y
        ENDELSE
    endif else imf_bz_pre=data(*,45)

    pretime = 3600.
    qualify_time=3600.

    imf=data(*,45)
    ntime = N_ELEMENTS(time)
    imf_pre = DBLARR(ntime)
    imf_pre(*) = !VALUES.F_NAN
    qualify_imf = DBLARR(ntime)
    qualify_imf(*) = !VALUES.F_NAN
    qualify_num = DBLARR(ntime)
    qualify_num(*) = !VALUES.F_NAN
    FOR itime = 0l, ntime-1 DO BEGIN
        index = where(time GT (time(itime)-pretime) AND time LT time(itime), ct)
        IF ct GE 6 THEN BEGIN
            IF  TOTAL(ABS(imf(index)) GE 0) GE 6 THEN BEGIN
                imf_pre(itime) = TOTAL(imf(index), /NAN)/ct
                if imf_pre(itime) gt 0 and (1.*total(imf(index) gt 0)/ct) ge qualify_time/pretime then begin 
                    qualify_imf(itime) = 1 
                    qualify_num(itime) = ct 
                endif else begin 
                    if imf_pre(itime) lt 0 and (1.*total(imf(index) lt 0)/ct) ge qualify_time/pretime then begin 
                        qualify_imf(itime) = -1 
                        qualify_num(itime) = ct
                    endif else qualify_imf(itime) =0
                endelse
            ENDIF
        Endif
    ENDFOR
; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(MAP) THEN BEGIN
        plot_path = main_path+'IMF_long_Bz/'  & spawn, 'mkdir ' + plot_path
        FOR ikp = 0, 1 DO BEGIN
            IF ikp EQ 0 THEN  BEGIN 
                sort_path = 'IMF_Bz_ge_'+sort_long_IMF_Bz_str+'/' 
                sort_flag = 1./(imf_bz_pre GE sort_long_IMF_Bz and qualify_imf gt 0 and qualify_num ge 10)
                sort_title = '(long) IMF Bz > '+sort_long_IMF_Bz_str+'  '
            ENDIF ELSE BEGIN
                sort_path = 'IMF_Bz_lt_neg'+sort_long_IMF_Bz_str+'/' 
                sort_flag = 1./(imf_bz_pre LT -sort_long_IMF_Bz and qualify_imf lt 0 and qualify_num ge 10)
                sort_title = '(long) IMF Bz < -'+sort_long_IMF_Bz_str+'  '
            ENDELSE 
            plot_path_long_imf_bz = plot_path+sort_path & spawn, 'mkdir ' + plot_path_long_imf_bz
            
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = plot_path_long_imf_bz, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              symmetry_template = symmetry_template, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF
ENDIF   

;------------------------------------
; sort by sw velocity
;----------------------------------
IF KEYWORD_SET(sort_sw_v) THEN BEGIN 
    sort_sw_v_str = STRING(sort_sw_v, format = '(i3.3)')
;pre the solar wind data
;sort_recalc=1
    if keyword_set(pre1h) then begin 
        IF ct_sort_store EQ 0 OR keyword_set(sort_recalc) THEN $
          sw_v_pre = predata_average(sw_v, time) $
        ELSE BEGIN 
            get_data, 'V', data = sort_data
            sw_v_pre = sort_data.y
        ENDELSE  
    endif else sw_v_pre=data(*,46)

; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(MAP) THEN BEGIN
        plot_path = main_path+'sw_v/'  & spawn, 'mkdir ' + plot_path
        FOR ikp = 0, 1 DO BEGIN
            IF ikp EQ 0 THEN  BEGIN
                sort_path = 'sw_v_ge_'+sort_sw_v_str+'/' 
                sort_flag = 1./(sw_v_pre GE sort_sw_V)
                sort_title = 'sw V > '+sort_sw_V_str+'  '
            ENDIF ELSE BEGIN 
                sort_path = 'sw_v_lt_'+sort_sw_V_str+'/' 
                sort_flag = 1./(sw_v_pre LT sort_sw_V)
                sort_title = 'sw V < '+sort_sw_V_str+'  '
            ENDELSE 
            plot_path_sw_v = plot_path+sort_path & spawn, 'mkdir ' + plot_path_sw_v
            
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = plot_path_sw_v, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              symmetry_template = symmetry_template, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF  
ENDIF  

;---------------------------------
; sort by Solar Wind Pressure
;---------------------------------
IF KEYWORD_SET(sort_sw_p) THEN BEGIN 
    sort_swp_str = strcompress(sort_sw_p, /remove_all)
;pre the solar wind data
    IF ct_sort_store EQ 0  OR keyword_set(sort_recalc) THEN $
      sw_p_pre = predata_average(sw_p, time) $
    ELSE BEGIN 
        get_data, 'P', data = sort_data
        sw_p_pre = sort_data.y
    ENDELSE 
; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(map) THEN BEGIN 
        plot_path = main_path+'SW_P/'  & spawn, 'mkdir ' + plot_path
        FOR iswp = 0, 1 DO BEGIN 
            IF iswp EQ 0 THEN  BEGIN 
                sort_path = 'SW_P_ge_'+sort_swp_str+'/' 
                sort_flag =  1./(sw_p_pre GE sort_sw_p) 
                sort_title = 'SW_P > '+sort_swp_str
            ENDIF ELSE BEGIN 
                sort_path = 'SW_P_ge_'+sort_swp_str+'/' 
                sort_flag =  1./(sw_p_pre LT sort_sw_p)
                sort_title = 'SW_P < '+sort_swp_str
            ENDELSE 
            plot_path_swp = plot_path+sort_path & spawn, 'mkdir ' + plot_path_swp         
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = plot_path_swp, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF 
ENDIF   
;-----------------------------------------
; sort by B magnitud
;----------------------------------------------
IF KEYWORD_SET(sort_IMF_B) THEN BEGIN 
    sort_IMF_mag_str = STRING(sort_IMF_b, format = '(i1.1)')
;pre the olar wind data
    IF ct_sort_store EQ 0  OR keyword_set(sort_recalc) THEN $
      imf_b_pre = predata_average(imf_b, time) $
    ELSE BEGIN 
        get_data, 'B', data = sort_data
        imf_b_pre = sort_data.y
    ENDELSE 
; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(map) THEN BEGIN 
        plot_path = main_path+'IMF_MAG/'  & spawn, 'mkdir ' + plot_path
        FOR imag = 0, 1 DO BEGIN 
            IF imag EQ 0 THEN  BEGIN 
                sort_path = 'IMF_MAG_ge_'+sort_IMF_mag_str+'/' 
                sort_flag =  1./(imf_b_pre GE sort_IMF_B) 
                sort_title = 'IMF MAG > '+sort_IMF_mag_str
            ENDIF ELSE BEGIN 
                sort_path = 'IMF_MAG_ge_'+sort_IMF_mag_str+'/' 
                sort_flag =  1./(imf_b_pre LT sort_IMF_B)
                sort_title = 'IMF MAG < '+sort_IMF_mag_str
            ENDELSE        
            plot_path_imf_mag = plot_path+sort_path & spawn, 'mkdir ' + plot_path_imf_mag        
            
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot,  $
              plot_path = plot_path_imf_mag, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF  
ENDIF 

; sort by IMF By and Bz
IF KEYWORD_SET(sort_IMF_By_Bz) THEN BEGIN
    sort_imf_by_bzy=sort_imf_by_bz(0) > 0 ;negative input is considered to be 0 here
    sort_IMF_By_Bzy_str = STRING(sort_IMF_By_Bzy, format = '(i1.1)')
    sort_IMF_By_Bzz= sort_imf_by_bz(1) ; set the Bz sort point different as IMF By
    sort_IMF_By_Bzz_str =  STRING(sort_IMF_By_Bzz, format = '(i1.1)')
;average solar wind By over 1 hour before each time point
    IF ct_sort_store EQ 0  OR keyword_set(sort_recalc) THEN  BEGIN 
        imf_By_pre = predata_average(imf_by, time) 
        imf_Bz_pre = predata_average(imf_bz, time)
    ENDIF  ELSE BEGIN 
        get_data, 'By', data = sort_data
        imf_by_pre = sort_data.y
        get_data, 'Bz', data = sort_data
        imf_bz_pre = sort_data.y
    ENDELSE    
; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(MAP) THEN BEGIN 
        plot_path = main_path+'IMF_By_Bz/'  & spawn, 'mkdir ' + plot_path
        byz_path = ['IMF_negBy_'+sort_IMF_By_Bzy_str+'_negBz_'+ sort_IMF_By_Bzz_str+'/',  $
                    'IMF_negBy_'+sort_IMF_By_Bzy_str+'_posBz_'+ sort_IMF_By_Bzz_str+'/', $
                    'IMF_posBy_'+sort_IMF_By_Bzy_str+'_posBz_'+ sort_IMF_By_Bzz_str+'/', $
                    'IMF_posBy_'+sort_IMF_By_Bzy_str+'_negBz_'+ sort_IMF_By_Bzz_str+'/']
        byz_flag = [[1./(imf_by_pre LT -sort_IMF_By_Bzy AND imf_bz_pre LT -sort_IMF_By_Bzz)], $
                    [1./(imf_by_pre LT -sort_IMF_By_Bzy AND imf_bz_pre GE sort_IMF_By_Bzz)], $
                    [1./(imf_by_pre GE sort_IMF_By_Bzy AND imf_bz_pre GE sort_IMF_By_Bzz)], $
                    [1./(imf_by_pre GE sort_IMF_By_Bzy AND imf_bz_pre LT -sort_IMF_By_Bzz)]]
        byz_title = ['IMF By < -'+ sort_IMF_By_Bzy_str+' and IMF Bz < -'+ sort_IMF_By_Bzz_str, $
                     'IMF By < -'+ sort_IMF_By_Bzy_str +' and IMF Bz > '+ sort_IMF_By_Bzz_str, $
                     'IMF By > '+ sort_IMF_By_Bzy_str+' and IMF Bz > '+ sort_IMF_By_Bzz_str, $
                     'IMF By > '+ sort_IMF_By_Bzy_str +' and IMF Bz < -'+ sort_IMF_By_Bzz_str]
        FOR ibyz = 0,  3 DO BEGIN 
            sort_path = byz_path(ibyz) 
            sort_flag = byz_flag(*, ibyz)
            sort_title = byz_title(ibyz)           
            plot_path_imf_byz = plot_path+sort_path & spawn, 'mkdir ' + plot_path_imf_byz
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot, $
              plot_path = plot_path_imf_byz, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              symmetry_template = symmetry_template, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF
ENDIF    

; sort by IMF By and long Bz
IF KEYWORD_SET(sort_long_IMF_By_Bz) THEN BEGIN
    sort_long_imf_by_bzy = sort_long_imf_by_bz(0) > 0 ;negative input is considered to be 0 here
    sort_long_IMF_By_Bzy_str = STRING(sort_long_IMF_By_Bzy, format = '(i1.1)')
    sort_long_IMF_By_Bzz = sort_long_imf_by_bz(1) > 0 ; set the Bz sort point different as IMF By
    sort_long_IMF_By_Bzz_str =  STRING(sort_long_IMF_By_Bzz, format = '(i1.1)')
;average solar wind By over 1 hour before each time point
    IF ct_sort_store EQ 0  OR keyword_set(sort_recalc) THEN  BEGIN 
        imf_By_pre = predata_average(imf_by, time)
        imf_Bz_pre = predata_average(imf_bz, time)
    ENDIF  ELSE BEGIN
        get_data, 'By', data = sort_data
        imf_by_pre = sort_data.y
        get_data, 'Bz', data = sort_data
        imf_bz_pre = sort_data.y
    ENDELSE    

    pretime = 3600.
    qualify_time=3600.

    imf=data(*,45)
    ntime = N_ELEMENTS(time)
    imf_pre = DBLARR(ntime)
    imf_pre(*) = !VALUES.F_NAN
    qualify_imf = DBLARR(ntime)
    qualify_imf(*) = !VALUES.F_NAN
    FOR itime = 0l, ntime-1 DO BEGIN
        index = where(time GT (time(itime)-pretime) AND time LT time(itime), ct)
        IF ct GE 6 THEN BEGIN
            IF  TOTAL(ABS(imf(index)) GE 0) GE 6 THEN BEGIN
                imf_pre(itime) = TOTAL(imf(index), /NAN)/ct
                if imf_pre(itime) gt 0 and (1.*total(imf(index) gt 0)/ct) ge qualify_time/pretime then begin 
                    qualify_imf(itime) = 1
                    qualify_num(itime) = ct 
                endif else begin 
                    if imf_pre(itime) lt 0 and (1.*total(imf(index) lt 0)/ct) ge qualify_time/pretime then begin 
                        qualify_imf(itime) = -1 
                        qualify_num(itime) = ct 
                    endif  else qualify_imf(itime) =0
                endelse
            ENDIF
        Endif
    ENDFOR
; MAP UNDER DIFFERENT CONDITION IF NEEDED
    IF KEYWORD_SET(MAP) THEN BEGIN 
        plot_path = main_path+'IMF_By_long_Bz/'  & spawn, 'mkdir ' + plot_path
        byz_path = ['IMF_negBy_'+sort_long_IMF_By_Bzy_str+'_negBz_'+ sort_long_IMF_By_Bzz_str+'/',  $
                    'IMF_negBy_'+sort_long_IMF_By_Bzy_str+'_posBz_'+ sort_long_IMF_By_Bzz_str+'/', $
                    'IMF_posBy_'+sort_long_IMF_By_Bzy_str+'_posBz_'+ sort_long_IMF_By_Bzz_str+'/', $
                    'IMF_posBy_'+sort_long_IMF_By_Bzy_str+'_negBz_'+ sort_long_IMF_By_Bzz_str+'/']
        byz_flag = [[1./(imf_by_pre LT -sort_long_IMF_By_Bzy AND imf_bz_pre LT -sort_long_IMF_By_Bzz $
                         and (qualify_imf) LT 0 and qualify_num ge 10) ], $
                    [1./(imf_by_pre LT -sort_long_IMF_By_Bzy AND imf_bz_pre GE sort_long_IMF_By_Bzz $
                         and (qualify_imf) GT 0 and qualify_num ge 10)  ], $
                    [1./(imf_by_pre GE sort_long_IMF_By_Bzy AND imf_bz_pre GE sort_long_IMF_By_Bzz $
                         and (qualify_imf) GT 0 and qualify_num ge 10)  ], $
                    [1./(imf_by_pre GE sort_long_IMF_By_Bzy AND imf_bz_pre LT -sort_long_IMF_By_Bzz $
                         and (qualify_imf) LT 0 and qualify_num ge 10)  ]]
        byz_title = ['IMF By < -'+ sort_long_IMF_By_Bzy_str+' and long IMF Bz < -'+ sort_long_IMF_By_Bzz_str, $
                     'IMF By < -'+ sort_long_IMF_By_Bzy_str +' and long IMF Bz > '+ sort_long_IMF_By_Bzz_str, $
                     'IMF By > '+ sort_long_IMF_By_Bzy_str+' and long IMF Bz > '+ sort_long_IMF_By_Bzz_str, $
                     'IMF By > '+ sort_long_IMF_By_Bzy_str +' and long IMF Bz < -'+ sort_long_IMF_By_Bzz_str]
        FOR ibyz = 0,  3 DO BEGIN 
            sort_path = byz_path(ibyz) 
            sort_flag = byz_flag(*, ibyz)
            sort_title = byz_title(ibyz)           
            plot_path_long_imf_byz = plot_path+sort_path & spawn, 'mkdir ' + plot_path_long_imf_byz
            map_o_beam, time, data, $
              sort_flag = sort_flag, sort_title = sort_title, $
              sc = sc, $
              ps_plot = ps_plot, $
              plot_path = plot_path_long_imf_byz, $
              coor_set = coor_set, $
              grid_set = grid_set, slice_grid_set = slice_grid_set, $
              direction_set = direction_set, $
              storm_phase_set = storm_phase_set, $
              plot_2d = plot_2d, slice_plot = slice_plot, $
              waterdrop_plot = waterdrop_plot, $
              point_plot = point_plot, $
              events_map = events_map,  $
              property_map_set = property_map_set, $
              diff_beta = diff_beta, $
              symmetry_template = symmetry_template, $
              property_map_type_set =  property_map_type_set
        ENDFOR 
    ENDIF
ENDIF    

; Because of the typical Parker spiral orientation of IMF, we use
; small IMF Bx here and then plot IMF By map to tell the real
; influence of IMF By

IF KEYWORD_SET(sort_IMF_By_smallBx) THEN BEGIN 
    small_imf_bx=3
    small_imf_bx_str= STRING(small_imf_bx, format = '(i1.1)')

    sort_imf_by_smallbx_set=sort_imf_by_smallbx
    for i_sort_imf_by_smallbx=0,n_elements(sort_imf_by_smallbx_set)-1 do begin 
        sort_imf_by_smallbx= sort_imf_by_smallbx_set(i_sort_imf_by_smallbx)
        if sort_IMF_By_smallbx lt 0 then sort_IMF_By_smallbx =0
        sort_IMF_By_smallbx_str = STRING(sort_IMF_by_smallbx, format = '(i1.1)')
;average solar wind By over 1 hour before each time point
        if keyword_set(delay_to_cluster) then  begin 
            imf_by_pre=data(*,55)  
            imf_bx_pre=data(*,54)
            delay_location='Delay_to_CLUSTER'
        endif else begin 
            IF keyword_set(pre1h) THEN BEGIN 
                IF ct_sort_store EQ 0  OR keyword_set(sort_recalc) THEN  begin 
                    imf_By_pre = predata_average(imf_by, time) 
                    imf_Bx_pre = predata_average(imf_bx, time)
                endif ELSE BEGIN 
                    get_data, 'By', data = sort_data
                    imf_by_pre = sort_data.y
                    get_data, 'Bx', data = sort_data
                    imf_bx_pre = sort_data.y
                ENDELSE 
                delay_location=''
            ENDIF ELSE begin 
                imf_by_pre = data(*,44)
                imf_bx_pre = data(*,43)
                delay_location='nodelay'
            endelse 
        endelse 
; MAP UNDER DIFFERENT CONDITION IF NEEDED
        IF KEYWORD_SET(MAP) THEN BEGIN 
            plot_path = main_path+'SmallBx_IMF_By_'+delay_location+'/'  & spawn, 'mkdir ' + plot_path
            FOR iby = 0, 1 DO BEGIN 
                IF iby EQ 0 THEN  BEGIN 
                    sort_path = 'IMF_Bxmag_lt_'+small_imf_bx_str+'_By_ge_'+ sort_IMF_By_smallBx_str+'/' 
                    sort_flag = 1./(imf_by_pre GE sort_IMF_By_smallBx and ABS(imf_bx_pre) lt small_imf_bx)
                    sort_title = '|Bx| < '+small_imf_bx_str +'IMF By > '+ sort_IMF_By_smallBx_str+'  '+delay_location
                ENDIF ELSE BEGIN 
                    sort_path = 'IMF_Bxmag_lt_'+small_imf_bx_str+'_By_lt_neg'+ sort_IMF_By_smallBx_str+'/' 
                    sort_flag = 1./(imf_by_pre LT -sort_IMF_By_smallBx and ABS(imf_bx_pre) lt small_imf_bx)
                    sort_title = '|Bx| < '+small_imf_bx_str+ 'IMF By < -'+ sort_IMF_By_smallBx_str+'  '+delay_location
                ENDELSE 
                plot_path_imf_by_smallBx = plot_path+sort_path & spawn, 'mkdir ' + plot_path_imf_by_smallBx
                map_o_beam, time, data, $
                  sort_flag = sort_flag, sort_title = sort_title, $
                  sc = sc, $
                  ps_plot = ps_plot,  $
                  plot_path = plot_path_imf_by_smallBx, $
                  coor_set = coor_set, $
                  grid_set = grid_set, slice_grid_set = slice_grid_set, $
                  direction_set = direction_set, $
                  storm_phase_set = storm_phase_set, $
                  plot_2d = plot_2d, slice_plot = slice_plot, $
                  waterdrop_plot = waterdrop_plot, $
                  point_plot = point_plot, $
                  events_map = events_map,  $
                  property_map_set = property_map_set, $
                  diff_beta = diff_beta, $
                  symmetry_template = symmetry_template, $
                  property_map_type_set =  property_map_type_set
            ENDFOR  
        ENDIF  
    endfor     
ENDIF

IF KEYWORD_SET(sort_season_imf_by) THEN BEGIN
    plot_path = main_path+'sort_season_imf_by/'  & spawn, 'mkdir ' + plot_path
    season=['spring','summer','fall','winter','Engwall']
    sort_season_imf_by_set= [0,3]
    for iseason = 4, 4 do begin 
        case iseason of 
            0: season_flag = 1./((time ge time_double('2001-05-01') and time lt time_double('2001-08-01')) or $
                                 (time ge time_double('2002-05-01') and time lt time_double('2002-08-01')))
            1: season_flag = 1./((time ge time_double('2001-08-01') and time lt time_double('2001-11-01')) or $
                                 (time ge time_double('2002-08-01') and time lt time_double('2002-11-01')))
            2: season_flag = 1./((time ge time_double('2001-11-01') and time lt time_double('2002-02-01')) or $
                                 (time ge time_double('2002-11-01') and time lt time_double('2003-02-01')))
            3: season_flag = 1./((time ge time_double('2001-02-01') and time lt time_double('2001-05-01')) or $
                                 (time ge time_double('2002-02-01') and time lt time_double('2002-05-01')))
            4: season_flag = 1./((time ge time_double('2001-07-03') and time lt time_double('2001-11-04')) or $
                                 (time ge time_double('2002-07-03') and time lt time_double('2002-11-04')))
        endcase
        
        for i_sort_season_imf_by=0,n_elements(sort_season_imf_by_set)-1 do begin 
            sort_season_imf_by=sort_season_imf_by_set(i_sort_season_imf_by)
            sort_season_IMF_By_str = STRING(sort_season_imf_by, format = '(i1.1)')
;average solar wind By over 1 hour before each time point 
            if keyword_set(delay_to_cluster) then  begin 
                imf_by_pre=data(*,55)  
                delay_location='Delay_to_CLUSTER'
            endif else begin 
                IF keyword_set(pre1h) THEN BEGIN 
                    IF ct_sort_store EQ 0  OR keyword_set(sort_recalc) THEN  imf_By_pre = predata_average(imf_by, time) $
                    ELSE BEGIN 
                        get_data, 'By', data = sort_data
                        imf_by_pre = sort_data.y
                    ENDELSE 
                    delay_location=''
                ENDIF ELSE begin 
                    imf_by_pre = data(*,44)
                    delay_location='nodelay'
                endelse 
            endelse 
; MAP UNDER DIFFERENT CONDITION IF NEEDED
            IF KEYWORD_SET(MAP) THEN BEGIN 
                FOR iby = 0, 1 DO BEGIN 
                    IF iby EQ 0 THEN  BEGIN 
                        sort_path = season(iseason)+'_IMF_By_ge_'+ sort_season_IMF_By_str+'/' 
                        sort_flag = 1./(imf_by_pre GE sort_season_imf_by and season_flag eq 1)
                        sort_title = season(iseason)+'IMF By > '+ sort_season_imf_by_str+'  '+delay_location
                    ENDIF ELSE BEGIN 
                        sort_path = season(iseason)+'_IMF_By_lt_neg'+ sort_season_IMF_By_str+'/' 
                        sort_flag = 1./(imf_by_pre LT -sort_season_imf_by and season_flag eq 1)
                        sort_title = season(iseason)+'IMF By < -'+ sort_season_imf_by_str+'  '+delay_location
                    ENDELSE 
                    plot_path_season_imf_by = plot_path+sort_path & spawn, 'mkdir ' + plot_path_season_imf_by

                    map_o_beam, time, data, $
                      sort_flag = sort_flag, sort_title = sort_title, $
                      sc = sc, $
                      ps_plot = ps_plot, $
                      plot_path = plot_path_season_imf_by, $
                      coor_set = coor_set, $
                      grid_set = grid_set, slice_grid_set = slice_grid_set, $
                      direction_set = direction_set, $
                      storm_phase_set = storm_phase_set, $
                      plot_2d = plot_2d, slice_plot = slice_plot, $
                      waterdrop_plot = waterdrop_plot, $
                      point_plot = point_plot, $
                      events_map = events_map, $
                      property_map_set = property_map_set, $
                      diff_beta = diff_beta, $
                      property_map_type_set =  property_map_type_set
                    
                ENDFOR                 
            ENDIF 
        endfor
    endfor   
ENDIF   

;----------------------------------------------
; save pre imf data
;----------------------------------------------
IF keyword_set(save_sort) THEN BEGIN 
    if n_elements(kp_index) gt 0 then store_data, 'KP_Index', data = {x:time, y:kp_index} ;this is not pre 1 hour data
    if n_elements(ae_index) gt 0 then store_data, 'AE_Index', data = {x:time, y:ae_index} ;this is not pre 1 hour data 
    if n_elements(imf_by_pre) gt 0 then store_data, 'By', data = {x:time, y:imf_by_pre}
    if n_elements(imf_bz_pre) gt 0 then store_data, 'Bz', data = {x:time, y:imf_bz_pre}
    if n_elements(imf_bx_pre) gt 0 then store_data, 'Bx', data = {x:time, y:imf_bx_pre}
    if n_elements(imf_b_pre) gt 0 then store_data, 'B', data = {x:time, y:imf_b_pre}
    if n_elements(sw_p_pre) gt 0 then store_data, 'P', data = {x:time, y:sw_p_pre}
    if n_elements(sw_v_pre) gt 0 then store_data,'V', data ={x:time,y:sw_v_pre}
    if n_elements(storm_phase) gt 0 then store_data, 'storm_phase', data = {x:time, y:storm_phase} 
    if n_elements(flag) gt 0 then store_data, 'flag', data = {x:time, y:flag} 
    if n_elements(x_gse) gt 0 then store_data, 'x_gse', data = {x:time, y:x_gse}
    if n_elements(f107_index) gt 0 then store_data,'F10_7_Index',data={x:time,y:f107_index}
    tplot_save, ['AE_Index', 'KP_Index', 'x_gse', 'Bx','By', 'Bz', 'B', 'P', 'V', 'storm_phase', 'flag','F10_7_Index','Substorm_flag'], $
      filename = flndata_sort
ENDIF 
;---------------------------------------
;calc distributions for different years and plot them
;---------------------------------------
IF keyword_set(property_distribution) THEN BEGIN  
  fluxdrop_effect=0
    if keyword_set(fluxdrop_effect) then eflux_filter=1
;IMF and solar wind parameters distribution can also be made here.
    index=where(vxbz eq 0,ct)
    if ct gt 0 then vxbz(index) = !values.f_nan
    mean_median_error = 1
    nsort = 1
    sort_condi_set = ''
    sort_path = 'non_sort/'
    ps=ps_plot
    if keyword_set(eflux_filter) then begin 
        sort_path = 'eflux_filter/'
        sort_condi_set = 'eflux_ge_'+eflux_threshold_str
    ENDIF
    IF keyword_set(sort_anodes) THEN BEGIN 
        sort_path = 'anodes_sort/'
        sort_condi_set = 'anodes_eq_'+string(indgen(4)+1, format = '(i1.1)')+'_or_'+string(8-indgen(4), format = '(i1.1)')
    ENDIF
    if keyword_set(eflux_filter) and  keyword_set(sort_anodes) then begin 
        sort_path = 'anodes_sort_'+'eflux_filter/'
        sort_condi_set = 'anodes_eq_'+string(indgen(4)+1, format = '(i1.1)')+'_or_'+string(8-indgen(4), format = '(i1.1)')+'_eflux_ge_'+eflux_threshold_str
    endif 
    if keyword_set(sort_imf_by) then begin 
        if sort_imf_by lt 0 then sort_imf_by=0
        sort_path = 'IMFBy_sort/'
        sort_condi_set = ['IMFBy_gt_'+string(sort_imf_by, format = '(i1.1)'),'IMFBy_lt_neg'+string(sort_imf_by,format='(i1.1)')]
        if sort_imf_by eq 0 then sort_imf_by =-1
    endif 

    if keyword_set(sort_imf_bz) then begin 
        if sort_imf_bz lt 0 then sort_imf_bz = 0
        sort_path = 'IMFBz_sort/'
        sort_condi_set = ['IMFBz_gt_'+string(sort_imf_bz, format = '(i1.1)'),'IMFBz_lt_neg'+string(sort_imf_bz,format='(i1.1)')]
        if sort_imf_bz eq 0 then sort_imf_bz =-1
    endif 

    if keyword_set(energy_filter) then begin 
        sort_path = 'energy_filter/'
        sort_condi_set = 'energy_gt'+ energy_threshold_str(0)+'_lt'+energy_threshold_str(1)
    endif
    if  keyword_set(eflux_filter) and keyword_set(energy_filter) then begin 
        sort_path = 'energy_filter_'+'eflux_filter/' 
        sort_condi_set = 'energy_gt'+ energy_threshold_str(0)+'_lt'+energy_threshold_str(1) +'_eflux_ge_'+eflux_threshold_str
    endif 
    FOR isort = 0, n_elements(sort_condi_set)-1 DO BEGIN
        sort_condi = sort_condi_set(isort)
        sort_ind_tail = REPLICATE(1, ntime)
        sort_ind_earth = REPLICATE(1, ntime)
        if keyword_set(eflux_filter) then begin 
            sort_ind_tail = (eflux_tail ge eflux_threshold or flag eq 0)
            sort_ind_earth = (eflux_earth ge eflux_threshold or flag eq 0)
        ENDIF   
        IF keyword_set(sort_anodes) THEN BEGIN
            sort_ind_tail = (anodes_tail EQ (isort+1)) OR (anodes_tail EQ (8-isort))
            sort_ind_earth = (anodes_earth EQ (isort+1)) OR (anodes_earth EQ (8-isort))
        ENDIF
        if keyword_set(eflux_filter) and  keyword_set(sort_anodes) then begin 
            sort_ind_tail = (anodes_tail EQ (isort+1)) OR (anodes_tail EQ (8-isort)) and  (eflux_tail ge eflux_threshold or flag eq 0)
            sort_ind_earth = (anodes_earth EQ (isort+1)) OR (anodes_earth EQ (8-isort)) and  (eflux_earth ge eflux_threshold or flag eq 0)
        endif 
        if keyword_set(energy_filter) then begin 
            sort_ind_tail = (v_energy_tail ge energy_threshold(0) and v_energy_tail le energy_threshold(1)) or flag eq 0
            sort_ind_earth = (v_energy_earth ge energy_threshold(0) and v_energy_tail le energy_threshold(1)) or flag eq 0
        ENDIF
        if keyword_set(energy_filter) and keyword_set(eflux_filter) then begin 
            sort_ind_tail = ((energy_tail ge energy_threshold(0)  and v_energy_tail le energy_threshold(1) and eflux_tail ge eflux_threshold)) or flag eq 0
            sort_ind_earth = ((energy_earth ge energy_threshold(0) and v_energy_tail le energy_threshold(1) and eflux_earth ge eflux_threshold))  or flag eq 0
        ENDIF  
        if keyword_set(sort_imf_by) then begin 
            if isort eq 0 then begin 
                sort_ind_tail = (imf_by ge sort_imf_by)
                sort_ind_earth = (imf_by ge sort_imf_by)
            endif else begin 
                sort_ind_tail = (imf_by le -sort_imf_by)
                sort_ind_earth = (imf_by le -sort_imf_by)
            endelse 
        ENDIF
        if keyword_set(sort_imf_bz) then begin 
            if sort_imf_bz lt 0 then sort_imf_bz =0
            if isort eq 0 then begin 
                sort_ind_tail = (imf_bz ge sort_imf_bz)
                sort_ind_earth = (imf_bz ge sort_imf_bz)
            endif else begin 
                sort_ind_tail = (imf_bz le -sort_imf_bz)
                sort_ind_earth = (imf_bz le -sort_imf_bz)
            endelse 
            if sort_imf_bz eq 0 then sort_imf_bz =-1
        ENDIF
        
        plot_path = path+'plots/property_distribution/'
        spawn, 'mkdir '+plot_path
        plot_path = plot_path+time_str+'/'
        spawn, 'mkdir '+plot_path
        plot_path = plot_path +sort_path
        spawn, 'mkdir '+plot_path
        range_scale = ['small']
        plot_path = plot_path+range_scale + '_scale/'
        spawn, 'mkdir '+plot_path
        spawn, 'mkdir '+plot_path+'diff_years/'
        spawn, 'mkdir '+plot_path+'median_line_plots/'
        
        property_names =  property_map_set

        year = INTARR(ntime)
;   single_phase_set= ['alltime_alldata','alltime','storm_alldata','nonstorm_alldata','nonstorm', 'storm']
        single_phase_set= ['nonstorm','storm'] ;'nonstorm','storm'];'nonstorm','storm']
        region_name_set=['lobe'];'lobe','bl','ps']
;        region_name_set = ['polar_cap'];,'tail_lobe'];,'tail_lobe'];,'bl','All_Region'];['All_Region','outflow'];south_lobe','north_lobe','south_polar','north_polar'] ;, 'Xgse_le_neg5', 'Xgse_gt_neg5']
        if (te_str.year-ts_str.year+1) eq 7 and keyword_set(mean_median_error) then mme = dblarr(te_str.year-ts_str.year+1,n_elements(single_phase_set), n_elements(region_name_set),n_elements(property_names),3)
        
        FOR is = 0, n_elements(single_phase_set)-1 DO BEGIN
            single_phase = single_phase_set(is)
            FOR iy = ts_str.year, te_str.year DO BEGIN 
                year_str = strcompress(iy, /remove_all)
                index = where(time GE  time_double(year_str+'-01-01') $
                              AND time LT time_double(year_str+'-12-31/13:59:59'), ct)
                IF ct GT 0 THEN year(index) = iy
            ENDFOR
; plot_en part can be used to look
; into data with certain condition for different years
            plot_en = 0
            IF plot_en EQ 1 THEN BEGIN
                in2001 = where(energy_tail*flux_tail lt 1400 and year eq 2001 and abs(flag) ge 1) ; AND beta LT 0.05 AND x_gse LE -5)
                in2002 = where(energy_tail*flux_tail lt 1400 and year eq 2002 and abs(flag) ge 1) ; AND beta LT 0.05 AND x_gse LE -5)
                in2003 = where(energy_tail*flux_tail lt 1400 and year eq 2003 and abs(flag) ge 1) ; AND beta LT 0.05 AND x_gse LE -5)
                in2004 = where(energy_tail*flux_tail lt 1400 and year eq 2004 and abs(flag) ge 1) ; AND beta LT 0.05 AND x_gse LE -5)
                in2005 = where(energy_tail*flux_tail lt 1400 and year eq 2005 and abs(flag) ge 1) ; AND beta LT 0.05 AND x_gse LE -5)
                store_data, '2001', data = {x:time(in2001), y:energy_tail(in2001)}
                store_data, '2002', data = {x:time(in2002), y:energy_tail(in2002)}
                store_data, '2003', data = {x:time(in2003), y:energy_tail(in2003)}
                store_data, '2004', data = {x:time(in2004), y:energy_tail(in2004)}
                store_data, '2005', data = {x:time(in2005), y:energy_tail(in2005)}
                
                timespan, '2001/01/01', 5*365+1, /days
                
                options, '2001', 'psym', 7
                options, '2001', 'title', 'Energy of beams with flux Lt 12 at the tail lobe'
                ylim, '2001', 40, 40000, 1
                options, '2002', 'color', 1
                options, '2003', 'color', 2
                options, '2004', 'color', 3
                options, '2005', 'color', 4
                popen, 'fluxlt12.ps', /land
                tplot, '2001'
                tplot_panel, v = '2001', o = '2002', psym = 7
                tplot_panel, v = '2001', o = '2003', psym = 7
                tplot_panel, v = '2001', o = '2004', psym = 7
                tplot_panel, v = '2001', o = '2005', psym = 7
                pclose
            ENDIF 
                                ;      !p.multi = [0, 1, fix(nyear)+1]
            FOR ireg = 0, n_elements(region_name_set)-1 DO  BEGIN 
                region_name = region_name_set(ireg)

                FOR ip = 0, n_elements(property_names)-1 DO BEGIN 
                    IF property_names(ip) EQ 'By' THEN BEGIN 
                        property_tail = By
                        property_earth = !values.f_nan
                        CASE range_scale OF 
                            'small': property_range = [-100, 100] 
                            'large': property_range = [-200, 200]
                            ELSE: stop
                        ENDCASE 
                        property_grid = (max(property_range, /NAN)-min(property_range,/nan))/50.
                                ;     barcolor_input = 6
                                ;     bar_yrange_input = 750
                    ENDIF 
                    IF property_names(ip) EQ 'IMF_By' THEN BEGIN 
                        property_tail = imf_By
                        property_earth = !values.f_nan
                        CASE range_scale OF 
                            'small': property_range = [-10, 10] 
                            'large': property_range = [-40, 40]
                            ELSE: stop
                        ENDCASE 
                        property_grid = (max(property_range, /NAN)-min(property_range,/nan))/50.
                                ;     barcolor_input = 6
                                ;     bar_yrange_input = 750
                    ENDIF 
                    IF property_names(ip) EQ 'sw_p' THEN BEGIN 
                        property_tail = sw_p
                        property_earth = !values.f_nan
                        CASE range_scale OF 
                            'small': property_range = [1, 5] 
                            'large': property_range = [1, 40]
                            ELSE: stop
                        ENDCASE 
                        property_grid = (max(property_range, /NAN)-min(property_range,/nan))/50.
                                ;     barcolor_input = 6
                                ;     bar_yrange_input = 750
                    ENDIF 
                    IF property_names(ip) EQ 'sw_v' THEN BEGIN 
                        property_tail = sw_v
                        property_earth = !values.f_nan
                        CASE range_scale OF 
                            'small': property_range = [200, 1000] 
                            'large': property_range = [200, 1000]
                            ELSE: stop
                        ENDCASE 
                        property_grid = (max(property_range, /NAN)-min(property_range,/nan))/50.
                                ;     barcolor_input = 6
                                ;     bar_yrange_input = 750
                    ENDIF 
                    IF property_names(ip) EQ 'VxBz' THEN BEGIN 
                        property_tail = VxBz
                        property_earth = !values.f_nan
                        CASE range_scale OF 
                            'small': property_range = [-3, 0] 
                            'large': property_range = [0, 100]
                            ELSE: stop
                        ENDCASE 
                        property_grid = (max(property_range, /NAN)-min(property_range,/nan))/20.
                                ;     barcolor_input = 6
                                ;     bar_yrange_input = 750
                    ENDIF 
                    IF property_names(ip) EQ 'E' THEN BEGIN 
                        property_tail = E
                        property_earth = !values.f_nan
                        CASE range_scale OF 
                            'small': property_range = [-1, 1] 
                            'large': property_range = [0, 100]
                            ELSE: stop
                        ENDCASE 
                        property_grid = (max(property_range, /NAN)-min(property_range,/nan))/20.
                                ;     barcolor_input = 6
                                ;     bar_yrange_input = 750
                    ENDIF 
                    IF property_names(ip) EQ 'energy' THEN BEGIN 
                        property_tail = energy_tail
                        property_earth = energy_earth
                        CASE range_scale OF 
                            'small': property_range = [0, 2000] 
                            'large': property_range = [2000, 40000]
                            ELSE: stop
                        ENDCASE 
                        property_grid = max(property_range, /NAN)/50.
                                ;     barcolor_input = 6
                                ;     bar_yrange_input = 750
                    ENDIF 
                    IF property_names(ip) EQ 'pitch_angle' THEN BEGIN 
                        property_tail = pa_tail
                        property_earth = pa_earth
                        property_range = [0, 90] 
                        property_grid = 22.5
                                ;     barcolor_input = 7
                                ;     bar_yrange_input = 2500
                    ENDIF 
                    IF property_names(ip) EQ 'flux' THEN BEGIN 
                        property_tail = flux_tail
                        property_earth = flux_earth  
                        CASE range_scale OF 
                            'small': property_range = [0, 2600] 
                            'large': property_range = [0, 8000]
                            ELSE: stop
                        ENDCASE  
                        property_grid = max(property_range, /NAN)/50.
                                ;     barcolor_input = 3
                                ;     bar_yrange_input = 70
                    ENDIF   
                    IF property_names(ip) EQ 'density' THEN BEGIN 
                        property_tail = density_tail
                        property_earth = density_earth
                        CASE range_scale OF
                            'small': property_range = [0, 0.08] 
                            'large': property_range = [0, 0.2]
                            ELSE: stop
                        ENDCASE  
                        property_grid = max(property_range, /NAN)/20.
                                ;    bar_yrange_input = 180
                                ;    barcolor_input = 2
                    ENDIF 
                    IF property_names(ip) EQ 'velocity' THEN BEGIN 
                        property_tail = velocity_tail
                        property_earth = velocity_earth
                        CASE range_scale OF 
                            'small': property_range = [0, 80] 
                            'large': property_range = [0, 400]
                            ELSE: stop
                        ENDCASE  
                        property_grid = max(property_range, /NAN)/20.                
                                ;barcolor_input = 5                     
                                ;    bar_yrange_input = 250       
                    ENDIF 
                    IF property_names(ip) EQ 'nV' THEN BEGIN 
                        property_tail = Alog10(density_tail*velocity_tail*1e5) ;cm-2s-1
                        property_earth = Alog10(density_earth*velocity_earth*1e5) ;cm-2s-1
                        CASE range_scale OF 
                            'small': property_range = [2, 7] 
                            'large': property_range = [1, 8]
                            ELSE: stop
                        ENDCASE  
                        property_grid = max(property_range, /NAN)/30.                
                                ;    barcolor_input = 1
                                ;    bar_yrange_input = 120     
                    ENDIF                      
                    IF property_names(ip) EQ 'temperature' THEN BEGIN 
                        property_tail = temperature_tail
                        property_earth = temperature_earth
                        CASE range_scale OF 
                            'small': property_range = [0, 30] 
                            'large': property_range = [0, 100]
                            ELSE: stop
                        ENDCASE  
                        property_grid = max(property_range, /NAN)/15.                
                                ;     barcolor_input = 5
                                ;    bar_yrange_input = 250       
                    ENDIF 
; for v parallel, the sign only means tailward/earthward, so take the
; ABS value mix them up to get the distribution of the parallel value
                    IF property_names(ip) EQ 'Vpara' THEN BEGIN 
                        property_tail = abs(v_para_tail) 
                        property_earth = abs(v_para_earth)
                        CASE range_scale OF 
                            'small': property_range = [0, 100.] 
                            'large': property_range = [0., 400.]
                            ELSE: stop
                        ENDCASE  
                        property_grid = (max(property_range, /NAN)-min(property_range,/nan))/50.   
                    ENDIF 
                    IF property_names(ip) EQ 'Vperp' THEN BEGIN 
                        property_tail = v_perp_tail
                        property_earth = v_perp_earth
                        CASE range_scale OF 
                            'small': property_range = [0, 30.] 
                            'large': property_range = [0, 400.]
                            ELSE: stop
                        ENDCASE  
                        property_grid = max(property_range, /NAN)/50.   
                    ENDIF 
                    IF property_names(ip) EQ 'Vperp_over_Vpara' THEN BEGIN 
                        property_tail = v_perp_tail/v_para_tail
                        property_earth = v_perp_earth/v_para_earth
                        CASE range_scale OF 
                            'small': property_range = [-1., 1.] 
                            'large': property_range = [-40., 40.]
                            ELSE: stop
                        ENDCASE  
                        property_grid = (max(property_range, /NAN)-min(property_range,/nan))/50.   
                    ENDIF 
                    IF property_names(ip) EQ 'nVpara_over_B' THEN BEGIN 
                        property_tail = (density_tail*ABS(v_para_tail))/B
                        property_earth = (density_earth*ABS(v_para_earth))/B
                        CASE range_scale OF 
                            'small': property_range = [0, 0.025] 
                            'large': property_range = [0, 0.4]
                            ELSE: stop
                        ENDCASE  
                        property_grid = max(property_range, /NAN)/50. 
                                ;      barcolor_input = 1
                                ;bar_yrange_input = 250     
                    ENDIF 
                    IF property_names(ip) EQ 'eflux' THEN BEGIN 
                        property_tail = eflux_tail ; energy_tail*flux_tail 
                        property_earth = eflux_earth ; energy_earth*energy_earth
                        CASE range_scale OF 
                            'small': property_range = [0, 1e5] 
                            'large': property_range = [0, 1e5]
                            ELSE: stop
                        ENDCASE   
                        property_grid = max(property_range, /NAN)/50.  
                                ;       barcolor_input = 5
                                ;      bar_yrange_input = 80
                    ENDIF 
                    IF property_names(ip) EQ 'anodes' THEN BEGIN 
                        property_tail = anodes_tail
                        property_earth = anodes_earth
                        property_range = [1, 8] 
                        property_grid = max(property_range, /NAN)/20.  
                                ;       barcolor_input = 5
                                ;      bar_yrange_input = 80
                    ENDIF 
                    IF property_names(ip) EQ 'v_energy' THEN BEGIN 
                        property_tail =  velocity_tail^2*mass_o/2
                        property_earth = velocity_earth^2*mass_o/2
                        CASE range_scale OF 
                            'small': property_range = [0, 1600.] 
                            'large': property_range = [0, 4e4]
                            ELSE: stop
                        ENDCASE  
                        property_grid = max(property_range, /NAN)/50.         
                    ENDIF 
                    IF property_names(ip) EQ 'energy_v' THEN BEGIN 
                        property_tail =  sqrt(2*energy_tail/mass_o)
                        property_earth =  sqrt(2*energy_earth/mass_o)
                        CASE range_scale OF 
                            'small': property_range = [0, 100] 
                            'large': property_range = [0, 500]
                            ELSE: stop
                        ENDCASE  
                        property_grid = max(property_range, /NAN)/50.  
                    ENDIF 
                    a=0
                    if a eq 3 then begin 
                        IF keyword_set(ps_plot) THEN popen, plot_path+'diff_years/diff_years_'+sort_condi+'_' $
                          +property_names(ip)+'_'+region_name+'_'+single_phase+'.ps', /port ; else window,/free                       
                        old_pmulti = !p.multi
                        !p.multi = [0, 1, fix(nyear)+1]
                        
                        FOR iy = ts_str.year-2001,te_str.year-2001  DO BEGIN 
                            IF region_name EQ 'dusklobe_mainphase' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy and x_gse lT -5 and x_gse gt -15 $
                                                   AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail AND $
                                                   z_gsm lt -9 and z_gsm gt -11 and y_gsm gt 4 and y_gsm lt 12 $
                                                   and storm_phase eq 2, ctt)
                                index_earth = where(year EQ 2001+iy and x_gse lT -5 AND x_gse gt -15 $
                                                    and beta lT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth AND $
                                                    z_gsm lt -9 and z_gsm gt -11 and y_gsm gt 4 and y_gsm lt 12 $ 
                                                    and storm_phase eq 2, cte)
                            ENDIF
                            IF region_name EQ 'dusklobe_nonstorm' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy and x_gse lT -5 and x_gse gt -15 $
                                                   AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail AND $
                                                   z_gsm lt -9 and z_gsm gt -11 and y_gsm gt 4 and y_gsm lt 12 $
                                                   and (storm_phase eq 0 or storm_phase eq 5), ctt)
                                index_earth = where(year EQ 2001+iy and x_gse lT -5 AND x_gse gt -15 $
                                                    and beta lT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth AND $
                                                    z_gsm lt -9 and z_gsm gt -11 and y_gsm gt 4 and y_gsm lt 12 $ 
                                                    and (storm_phase eq 0 or storm_phase eq 5), cte)
                            ENDIF

                            IF region_name EQ 'polar_cap' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND x_gse GT -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail AND ABS(z_gsm) ge 4, ctt)
                                index_earth = where(year EQ 2001+iy AND x_gse GT -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth AND ABS(z_gsm) ge 4, cte)
                            ENDIF 
                            IF region_name EQ 'tail_lobe' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND x_gse LE -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                                index_earth = where(year EQ 2001+iy AND x_gse LE -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                            ENDIF
                            IF region_name EQ 'south_lobe' THEN BEGIN
                                index_tail = where(year EQ 2001+iy AND x_gse LE -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail and z_gsm lt 0, ctt)
                                index_earth = where(year EQ 2001+iy AND x_gse LE -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth and z_gsm lt 0, cte)
                            ENDIF
                            IF region_name EQ 'north_lobe' THEN BEGIN
                                index_tail = where(year EQ 2001+iy AND x_gse LE -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail and z_gsm gt 0, ctt)
                                index_earth = where(year EQ 2001+iy AND x_gse LE -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth and z_gsm gt 0, cte)
                            ENDIF                        
                            IF region_name EQ 'Xgse_le_neg5' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND x_gse LE -5 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                                index_earth = where(year EQ 2001+iy AND x_gse LE -5 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                            ENDIF 
                            IF region_name EQ 'Xgse_gt_neg5' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND x_gse GT -5 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                                index_earth = where(year EQ 2001+iy AND x_gse GT -5 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                            ENDIF 
                            IF region_name EQ 'All_Region'THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND ABS(flag) Ge 0 AND sort_ind_tail EQ 1, ctt)
                                index_earth = where(year EQ 2001+iy AND ABS(flag) Ge 0 AND sort_ind_earth EQ 1, cte)
                            ENDIF 
                            IF region_name EQ 'south_polar' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND x_gse GT -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail AND z_gsm lt -4, ctt)
                                index_earth = where(year EQ 2001+iy AND x_gse GT -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth AND z_gsm lt -4, cte)
                            ENDIF 
                            IF region_name EQ 'north_polar' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND x_gse GT -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail AND z_gsm ge 4, ctt)
                                index_earth = where(year EQ 2001+iy AND x_gse GT -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth AND z_gsm gt 4, cte)
                            ENDIF
                            IF region_name EQ 'lobe' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                                index_earth = where(year EQ 2001+iy AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                            ENDIF
                            IF region_name EQ 'bl' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND beta GT 0.05 AND beta LT 1 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                                index_earth = where(year EQ 2001+iy AND beta GT 0.05 AND beta LT 1 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                            ENDIF 
                            IF region_name EQ 'ps' THEN BEGIN 
                                index_tail = where(year EQ 2001+iy AND beta GT 1 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                                index_earth = where(year EQ 2001+iy AND beta GT 1 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                            ENDIF
                            IF (ctt+cte) GT 0 THEN BEGIN
                                index=where(finite(property_tail),ctpt)
                                index=where(finite(property_earth),ctpe)
                                IF ctt GT 0 AND cte GT 0 and ctpt gt 0 and ctpe gt 0 THEN BEGIN 
                                    property = [property_tail(index_tail), property_earth(index_earth)]
                                    flag_1 =  [flag(index_tail)>0, flag(index_earth)<0]
                                    storm_phase_1 =  [storm_phase(index_tail), storm_phase(index_earth)]
                                ENDIF 
                                IF (ctt GT 0 and ctpt gt 0) AND (cte EQ  0 or ctpe eq 0) THEN BEGIN 
                                    property = [property_tail(index_tail)]
                                    flag_1 =  [flag(index_tail)>0]
                                    storm_phase_1 =  [storm_phase(index_tail)]
                                ENDIF 
                                IF (ctt EQ 0 or ctpt eq 0) AND (cte GT 0 and ctpe eq 0) THEN BEGIN 
                                    property = [property_earth(index_earth)]
                                    flag_1 =  [flag(index_earth)<0]
                                    storm_phase_1 = [storm_phase(index_earth)]
                                ENDIF 
                                index=where(ABS(property) gt 0 and ABS(flag gt 0), ct)
                                if ct gt 0 then result = data_distribution(property, flag=flag_1, storm_phase=storm_phase_1, xgrid = property_grid, region_name = region_name,para_name = strcompress(2001+iy, /remove_all)+sort_condi +'_'+property_names(ip), xrange = property_range, write_data = 0, events_distribution = 0, ratio_distribution = 0, ps = 0,  plot_single_phase_events_distribution = 1, barcolor_input = barcolor_input, bar_yrange_input = bar_yrange_input,single_phase = single_phase, charsize_input = 2) ;1+nyear*0.2)
                                if (te_str.year-ts_str.year+1) eq 7 and keyword_set(mean_median_error) then begin 
                                    tplot_names,'mean_median_error',names=names
                                    if names(0) ne '' then begin 
                                        get_data,'mean_median_error',data=data
                                        mme(iy,is,ireg,ip,0)=data.mean
                                        mme(iy,is,ireg,ip,1)=data.median
                                        mme(iy,is,ireg,ip,2)=data.error
                                    endif else begin
                                        mme(iy,is,ireg,ip,0)=!values.f_nan
                                        mme(iy,is,ireg,ip,1)=!values.f_nan
                                        mme(iy,is,ireg,ip,2)=!values.f_nan
                                    endelse 
                                endif 
                            ENDIF  
                        ENDFOR  
                        
                        IF keyword_set(ps_plot) THEN pclose  else stop 
                        !p.multi = old_pmulti
                    endif    
; now draw the the property distribution for all years 
                    IF keyword_set(ps_plot) THEN popen, plot_path+time_str+'_'+sort_condi+'_'+property_names(ip)+'_'+region_name+'_'+single_phase+'.ps', /land
                    
                    IF region_name EQ 'dusklobe_mainphase' THEN BEGIN 
                        index_tail = where( x_gse lT -13 and x_gse gt -15 $
                                            and beta LT 0.05 AND flag ge 1 AND sort_ind_tail $
                                            and  z_gsm lt -9 and z_gsm gt -11  and y_gsm gt 4 and y_gsm lt 12 $
                                            and storm_phase eq 2, ctt)
                        index_earth = where( x_gse lT -13 AND x_gse gt -15 $
                                             and  beta lT 0.05 AND (flag eq -1 or flag eq 2) AND sort_ind_earth $
                                             and  z_gsm lt -9 and z_gsm gt -11 and y_gsm gt 4 and y_gsm lt 12 $ 
                                             and storm_phase eq 2, cte)
                    ENDIF
                    IF region_name EQ 'dusklobe_nonstorm' THEN BEGIN 
                        index_tail = where( x_gse lT -13 and x_gse gt -15 $
                                            and beta LT 0.05 AND flag Ge 1 AND sort_ind_tail $
                                            and  z_gsm lt -9 and z_gsm gt -11 and y_gsm gt 4 and y_gsm lt 12 $
                                            and (storm_phase eq 0 or storm_phase eq 5), ctt)
                        index_earth = where(x_gse lT -13 AND x_gse gt -15 $
                                            and  beta lT 0.05 AND (flag eq -1 or flag eq 2) AND sort_ind_earth $
                                            and  z_gsm lt -9 and z_gsm gt -11 and y_gsm gt 4 and y_gsm lt 12 $ 
                                            and (storm_phase eq 0 or storm_phase eq 5), cte)
                    ENDIF
                    IF region_name EQ 'bl' THEN BEGIN 
                        index_tail = where(beta gt 0.05 AND beta LT 1 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                        index_earth = where(beta GT 0.05 AND beta LT 1 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                    ENDIF 
                    IF region_name EQ 'polar_cap' THEN BEGIN 
                        index_tail = where(x_gse GT -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail AND ABS(z_gsm) ge 4, ctt)
                        index_earth = where(x_gse GT -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth AND ABS(z_gsm) ge 4, cte)
                    ENDIF 
                    IF region_name EQ 'tail_lobe' THEN BEGIN 
                        index_tail = where(x_gse LE -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                        index_earth =  where(x_gse LE -5 AND beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                    ENDIF 
                    IF region_name EQ 'Xgse_le_neg5' THEN BEGIN 
                        index_tail = where(x_gse LE -5 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                        index_earth = where(x_gse LE -5 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                    ENDIF 
                    IF region_name EQ 'Xgse_gt_neg5' THEN BEGIN 
                        index_tail = where(x_gse GT -5 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                        index_earth = where(x_gse GT -5 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                    ENDIF 
                    IF region_name EQ 'All_Region'THEN BEGIN 
                        index_tail = where(ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                        index_earth = where(ABS(flag) Ge  0 AND sort_ind_earth, cte)
                    ENDIF 
                    IF region_name EQ 'lobe' THEN BEGIN 
                        index_tail = where( beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                        index_earth = where(beta LT 0.05 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                    ENDIF
                    
                    IF region_name EQ 'ps' THEN BEGIN 
                        index_tail = where(beta GT 1 AND ABS(flag) Ge 0 AND sort_ind_tail, ctt)
                        index_earth = where(beta GT 1 AND ABS(flag) Ge 0 AND sort_ind_earth, cte)
                    ENDIF

                    IF (ctt+cte) GT 0 THEN BEGIN    
                        IF ctt GT 0 AND cte GT 0 THEN BEGIN 
                            property = [property_tail(index_tail), property_earth(index_earth)]
                            flag_1 =  [flag(index_tail)>0, flag(index_earth)<0]
                            storm_phase_1 =  [storm_phase(index_tail), storm_phase(index_earth)]
                        ENDIF 
                        IF ctt GT 0 AND cte EQ  0 THEN BEGIN 
                            property = [property_tail(index_tail)]
                            flag_1 =  [flag(index_tail)>0]
                            storm_phase_1 =  [storm_phase(index_tail)]
                        ENDIF 
                        IF ctt EQ  0 AND cte GT 0 THEN BEGIN 
                            property = [property_earth(index_earth)]
                            flag_1 =  [flag(index_earth)<0]
                            storm_phase_1 = [storm_phase(index_earth)]
                        ENDIF                            
;plot,x_gse(index_tail),property,psym=1,xrange=[-5,-15],xstyle=1,xtitle='x_gse(Re)',ytitle='velocity',title=region_name,yrange=[0,100]
                        c=3
                        if c eq 3 then begin 
                            result = data_distribution(property, flag=flag_1,storm_phase= storm_phase_1, $
                                                       xgrid = property_grid, region_name = region_name, $
                                                       para_name = time_str +sort_condi $
                                                       +'_'+property_names(ip), xrange = property_range, $
                                                       write_data = 0, events_distribution = 0, $
                                                       ratio_distribution = 0,  ps = 0, $
                                                       plot_single_phase_events_distribution = 1, $
                                                       barcolor_input = barcolor_input, $
                                                       bar_yrange_input = bar_yrange_input, $
                                                       single_phase = single_phase, charsize_input = 1.2,quartile_range=0)
                           
                            if keyword_set(fluxdrop_effect) then begin
                                get_data,result,data=stadata
                                index1=where(stadata.x ge eflux_threshold,ct) 
                                
                                if ct gt 0 then begin 
                                    if keyword_set(ps_plot)then popen,'thesis_figure/solarcycle_fluxdrop_effect_nonstorm.ps'
                                ;100,1e5
                                    plot,stadata.x(index1),stadata.sta(index1,0,3),psym=-2,xlog=1,xrange=[100,1e5],xstyle=1,xtitle=property_names,ytitle='# of Events',title=region_name+time_str+sort_condi+' nonstorm'
                                  ;  oplot,stadata.x(index1)/8.,stadata.sta(index1,0,3),psym=-5,color=2
                                  ;  xyouts,350,650,'1/8',color=2,charsize=2
                                ;  oplot,stadata.x(index1)/16.,stadata.sta(index1,0,3),psym=-4,color=4
                                ;  xyouts,150,650,'1/16',color=4,charsize=2

                                    oplot,stadata.x(index1)*0.33,stadata.sta(index1,0,3),psym=-5,color=2
                                    xyouts,350,650,'33%',color=2,charsize=2
                                    oplot,stadata.x(index1)*0.24,stadata.sta(index1,0,3),psym=-4,color=4
                                    xyouts,150,650,'24%',color=4,charsize=2

                                    oplot,[eflux_threshold,eflux_threshold],[0,5000],colo=1
                                    num0=total(stadata.sta(index1,0,3))
                                 ;   index2=where(stadata.x/8 ge eflux_threshold,ct)
                                    index2=where(stadata.x*0.33 ge eflux_threshold,ct)
                                    if ct gt 0 then begin
                                        num2=total(stadata.sta(index2,0,3))
                                        ratio2=num2/num0
                                        print,num2,num0,ratio2
                                    endif 
                                 ;   index3=where(stadata.x/16. ge eflux_threshold,ct)
                                    index3=where(stadata.x*0.24 ge eflux_threshold,ct)
                                    if ct gt 0 then begin
                                        num3=total(stadata.sta(index3,0,3))
                                        ratio3=num3/num0
                                        print,num3,num0,ratio3
                                    endif 

                                    if keyword_set(ps_plot) then pclose else stop
                                    if keyword_set(ps_plot) then popen, 'thesis_figure/solarcycle_fluxdrop_effect_storm.ps'
                                    plot,stadata.x(index1),stadata.sta(index1,0,5),psym=-2,xlog=1,xrange=[100,1e5],xstyle=1,xtitle=property_names,ytitle='# of Events',title=region_name+time_str+sort_condi+' storm'
                                ;          oplot,stadata.x(index1)/2,stadata.sta(index1,0,5),psym=-5,color=3
                                 
                                    ;oplot,stadata.x(index1)/4.,stadata.sta(index1,0,5),psym=-4,color=2
                              ;      xyouts,400,130,'1/4',color=2,charsize=2
                               ;     oplot,stadata.x(index1)/13.,stadata.sta(index1,0,5),psym=-4,color=4
                                ;    xyouts,200,130,'1/13',color=4,charsize=2
                                 ;   oplot,[eflux_threshold,eflux_threshold],[0,5000],colo=1 

                                    oplot,stadata.x(index1)*0.19,stadata.sta(index1,0,5),psym=-4,color=2
                                    xyouts,400,130,'19%',color=2,charsize=2
                                    oplot,stadata.x(index1)*0.09,stadata.sta(index1,0,5),psym=-4,color=4
                                    xyouts,200,130,'9%',color=4,charsize=2
                                    oplot,[eflux_threshold,eflux_threshold],[0,5000],colo=1 

                                    num0=total(stadata.sta(index1,0,5)) 
                                 ;   index2=where(stadata.x/4 ge eflux_threshold,ct)
                                    index2=where(stadata.x*0.19 ge eflux_threshold,ct)
                                    if ct gt 0 then begin
                                        num2=total(stadata.sta(index2,0,5))
                                        ratio2=num2/num0
                                        print,num2,num0,ratio2
                                    endif 
                                 ;   index3=where(stadata.x/13. ge eflux_threshold,ct)
                                    index3=where(stadata.x*0.09 ge eflux_threshold,ct)
                                    if ct gt 0 then begin
                                        num3=total(stadata.sta(index3,0,5))
                                        ratio3=num3/num0
                                        print,num3,num0,ratio3
                                    endif 
                                    if keyword_set(ps_plot) then pclose else stop
                                endif      
                            endif   
                        endif  
                    ENDIF  
                    IF keyword_set(ps_plot) THEN pclose else stop
                ENDFOR               
            ENDFOR    
        ENDFOR 
        stop
        if (te_str.year-ts_str.year+1) eq 7 and keyword_set (mean_median_error) then begin 
            store_data,'mme_'+sort_condi,data={title:'iy,is,ireg,ip,ivalue',mme:mme}
            FOR ip = 0, n_elements(property_names)-1 DO BEGIN 
                FOR is = 0, n_elements(single_phase_set)-1 DO BEGIN
                    FOR ireg = 0, n_elements(region_name_set)-1 DO  BEGIN 
; have to be 7 years here
                        data = {x:[time_double('2001-07-01'),time_double('2002-07-01'), $
                                   time_double('2003-07-01'),time_double('2004-07-01'), $
                                   time_double('2005-07-01'), time_double('2006-07-01'),$
                                   time_double('2007-07-01') ],y:mme(*,is,ireg,ip,0),dy:mme(*,is,ireg,ip,2)}
                        store_data, single_phase_set(is)+'_'+region_name_set(ireg)+'_'+property_names(ip)+'_mean', $
                          data=data,lim={yrange:[0,max(mme(*,is,ireg,ip,0)+mme(*,is,ireg,ip,2),/nan)*1.1]}, $
                          dlim={ytitle: 'mean     '+property_names(ip)+'!C'+region_name_set(ireg)+'     '+ $x
                                single_phase_set(is)+'!C'+'energy: '+energy_threshold_str(0)+' - '+energy_threshold_str(1), $
                                psym:-7,thick:3}
                        
                        data = {x:[time_double('2001-07-01'),time_double('2002-07-01'),$
                                   time_double('2003-07-01'),time_double('2004-07-01'), $
                                   time_double('2005-07-01'), time_double('2006-07-01'),$
                                   time_double('2007-07-01')],y:mme(*,is,ireg,ip,1)}
                        store_data, single_phase_set(is)+'_'+region_name_set(ireg)+'_'+property_names(ip)+'_median', $
                          data=data ,lim={yrange:property_range}, $
                          dlim={ytitle: 'median     '+property_names(ip)+'!C'+region_name_set(ireg)+'     ' $
                                +single_phase_set(is)+'!C'+'energy: '+energy_threshold_str(0)+' - '+energy_threshold_str(1), $
                                psym:-7,thick:3}
                        
                        OPENU, unit, plot_path+'median_line_plots/data.txt', /GET_LUN, /APPEND
                        PRINTF, unit, sort_condi+'_'+property_names(ip)+'_'+region_name_set(ireg)+'_' $
                          +single_phase_set(is)+' mean  median  error'
                        for iy =0,4 do  PRINTF, unit, mme(iy,is,ireg,ip,0),mme(iy,is,ireg,ip,1),mme(iy,is,ireg,ip,2)
                        FREE_LUN, unit
                    endfor
                    
                    CASE property_names(ip) OF 
                        'flux': property_range = [0, 2600.] 
                        'density': property_range = [0, 0.07]
                        'nV': property_range = [0, 3]
                        'v_energy': property_range = [0, 1600.]
                        'nVpara_over_B': property_range=[0.0005,0.05]
                        ELSE: stop
                    ENDCASE  
                    
                    ps_plot=1
                    if keyword_set(ps_plot) then popen,plot_path+'median_line_plots/'+sort_condi+'_'+property_names(ip) +'_'+single_phase_set(is)+'_median.ps' ,/land 
                    plot,[2000.5,2007.5],[0.1,1],/nodata,xstyle=1,title = single_phase_set(is)+',  Energy: '+energy_threshold_str(0)+'-'+energy_threshold_str(1), xtitle='year', ytitle='median  '+property_names(ip),yrange = [10,5e3],charsize=1.5, position=[0.15,0.15,0.95,0.85],ylog=1,ystyle=1

                    oplot, [2001,2002,2003,2004,2005,2006,2007], mme(*,is,0,ip,1), psym=-1, color=6, thick=10
                                ;       oplot, [2001,2002,2003,2004,2005,2006,2007], mme(*,is,1,ip,1), psym=-1, color=2, thick=10
                    xyouts, 2004, property_range(1)/1.5, region_name_set(0), color=6, charsize=3, charthick = 3
                                ;        xyouts, 2004, property_range(1)/2., region_name_set(1), color=2, charsize=3, charthick = 3
                    
                    if keyword_set(ps_plot) then pclose else stop
                    
                    if keyword_set(ps_plot)  then popen, plot_path+'median_line_plots/'+sort_condi+'_'+property_names(ip)$
                      +'_'+single_phase_set(is)+'_mean.ps' ,/land else window,2
                    tplot, [single_phase_set(is)+'_'+region_name_set(0)+'_'+property_names(ip)+'_mean'] ;, single_phase_set(is) $
                                ;    +'_'+region_name_set(1)+'_'+property_names(ip)+'_mean']
                    if keyword_set(ps_plot) then pclose else stop
                endfor 
            endfor
        endif     
    ENDFOR          
    IF keyword_set(ps_plot) THEN   begin 
        spawn, 'mogrify -format png '+plot_path+'*.ps'  
        spawn, 'mogrify -rotate -90 '+plot_path+'*.png' 
        spawn, 'mkdir '+plot_path+'png/'
        spawn, 'mv -f '+plot_path+'*.png '+plot_path+'png/'
        spawn, 'gzip -9f '+plot_path+'*.ps'  
        spawn, 'mogrify -format png '+plot_path+'diff_years/*.ps'
        spawn, 'mkdir '+plot_path+'diff_years/png/'
        spawn, 'mv -f '+plot_path+'diff_years/*.png '+plot_path+'diff_years/png/'
        spawn, 'gzip -9f '+plot_path+'diff_years/*.ps'  
        spawn, 'mogrify -format png '+plot_path+'median_line_plots/*.ps'
        spawn, 'mogrify -rotate -90 '+plot_path+'median_line_plots/*.png' 
        spawn, 'mkdir '+plot_path+'median_line_plots/png/'
        spawn, 'mv -f '+plot_path+'median_line_plots/*.png '+plot_path+'median_line_plots/png/'
        spawn, 'gzip -9f '+plot_path+'median_line_plots/*.ps'  
    endif 
ENDIF        
;---------------------------------------
;calc distributions for different years and plot them
;---------------------------------------
IF keyword_set(year_distribution) THEN BEGIN
month_settings=[0] ; 0:all month; 1: tail_season; 2: tail_season_plus; 3: July_Aug; 4: Sep_oct;
path_year_main = 'output/o_beam/test_year_dependence_with_newfilter/'
eflux_filter=1
if keyword_set(eflux_filter) then  path_year_main =path_year_main+'eflux_filter_result/'+'eflux_ge_'+eflux_threshold_str+'/' else path_year_main=path_year_main+'result/'
for im=0,n_elements(month_settings)-1 do begin
    month_setting=month_settings(im)
    case month_setting of
        0:path_year=path_year_main
        1:path_year=path_year_main+'tail_season/'
        2:path_year=path_year_main+'tail_season_plus/'
        3:path_year=path_year_main+'July_Aug/'
        4:path_year=path_year_main+'Sep_Oct/'
    endcase 
    print,month_settings(im)
    spawn, 'mkdir '+path_year
    fln_year = FINDFILE(path_year+'year_efficience.tplot', count = ct)
    IF ct GT 0 THEN  BEGIN 
        tplot_restore, filename = fln_year  
        tplot_names, inst_name+'year', names = names
        IF names(0) NE '' THEN BEGIN 
            get_data, inst_name+'year', data = data
            year = data.y
        ENDIF ELSE rc = 1
    ENDIF ELSE rc = 1 
    rc=1
    IF rc EQ 1 THEN BEGIN
        year = INTARR(ntime)
        case month_setting of
            0: begin
                year(where(time GE  time_double('2001-01-01') AND time LT time_double('2001-12-31'))) = 2001
                year(where(time GE  time_double('2002-01-01') AND time LT time_double('2002-12-31'))) = 2002
                year(where(time GE  time_double('2003-01-01') AND time LT time_double('2003-12-31'))) = 2003
                year(where(time GE  time_double('2004-01-01') AND time LT time_double('2004-12-31'))) = 2004
                year(where(time GE  time_double('2005-01-01') AND time LT time_double('2005-12-31'))) = 2005
                year(where(time GE  time_double('2006-01-01') AND time LT time_double('2006-12-31'))) = 2006
                year(where(time GE  time_double('2007-01-01') AND time LT time_double('2007-12-31'))) = 2007
                year(where(time GE  time_double('2008-01-01') AND time LT time_double('2008-12-31'))) = 2008
                year(where(time GE  time_double('2009-01-01') AND time LT time_double('2009-12-31'))) = 2009
            end 
            1: begin
                year(where(time GE  time_double('2001-07-01') AND time LT time_double('2001-11-01'))) = 2001
                year(where(time GE  time_double('2002-07-01') AND time LT time_double('2002-11-01'))) = 2002
                year(where(time GE  time_double('2003-07-01') AND time LT time_double('2003-11-01'))) = 2003
                year(where(time GE  time_double('2004-07-01') AND time LT time_double('2004-11-01'))) = 2004
                year(where(time GE  time_double('2005-07-01') AND time LT time_double('2005-11-01'))) = 2005
                year(where(time GE  time_double('2006-07-01') AND time LT time_double('2006-11-01'))) = 2006
                year(where(time GE  time_double('2007-07-01') AND time LT time_double('2007-11-01'))) = 2007
                year(where(time GE  time_double('2008-07-01') AND time LT time_double('2008-11-01'))) = 2008
                year(where(time GE  time_double('2009-07-01') AND time LT time_double('2009-11-01'))) = 2009
            end 
            2: begin
                year(where(time GE  time_double('2001-07-01') AND time LT time_double('2001-12-01'))) = 2001
                year(where(time GE  time_double('2002-07-01') AND time LT time_double('2002-12-01'))) = 2002
                year(where(time GE  time_double('2003-07-01') AND time LT time_double('2003-12-01'))) = 2003
                year(where(time GE  time_double('2004-07-01') AND time LT time_double('2004-12-01'))) = 2004
                year(where(time GE  time_double('2005-07-01') AND time LT time_double('2005-12-01'))) = 2005
                year(where(time GE  time_double('2006-07-01') AND time LT time_double('2006-12-01'))) = 2006
                year(where(time GE  time_double('2007-07-01') AND time LT time_double('2007-12-01'))) = 2007
                year(where(time GE  time_double('2008-07-01') AND time LT time_double('2008-12-01'))) = 2008
                year(where(time GE  time_double('2009-07-01') AND time LT time_double('2009-12-01'))) = 2009
            end 
            3: begin
                year(where(time GE  time_double('2001-07-01') AND time LT time_double('2001-09-01'))) = 2001
                year(where(time GE  time_double('2002-07-01') AND time LT time_double('2002-09-01'))) = 2002
                year(where(time GE  time_double('2003-07-01') AND time LT time_double('2003-09-01'))) = 2003
                year(where(time GE  time_double('2004-07-01') AND time LT time_double('2004-09-01'))) = 2004
                year(where(time GE  time_double('2005-07-01') AND time LT time_double('2005-09-01'))) = 2005
                year(where(time GE  time_double('2006-07-01') AND time LT time_double('2006-09-01'))) = 2006
                year(where(time GE  time_double('2007-07-01') AND time LT time_double('2007-09-01'))) = 2007
                year(where(time GE  time_double('2008-07-01') AND time LT time_double('2008-09-01'))) = 2008
                year(where(time GE  time_double('2009-07-01') AND time LT time_double('2009-09-01'))) = 2009
            end 
            4: begin
                year(where(time GE  time_double('2001-09-01') AND time LT time_double('2001-11-01'))) = 2001
                year(where(time GE  time_double('2002-09-01') AND time LT time_double('2002-11-01'))) = 2002
                year(where(time GE  time_double('2003-09-01') AND time LT time_double('2003-11-01'))) = 2003
                year(where(time GE  time_double('2004-09-01') AND time LT time_double('2004-11-01'))) = 2004
                year(where(time GE  time_double('2005-09-01') AND time LT time_double('2005-11-01'))) = 2005
                year(where(time GE  time_double('2006-09-01') AND time LT time_double('2006-11-01'))) = 2006
                year(where(time GE  time_double('2007-09-01') AND time LT time_double('2007-11-01'))) = 2007
                year(where(time GE  time_double('2008-09-01') AND time LT time_double('2008-11-01'))) = 2008
                year(where(time GE  time_double('2009-09-01') AND time LT time_double('2009-11-01'))) = 2009
            end 
        endcase   
    ENDIF  
    
    region_name_all = ['north_lobe', 'south_lobe','north_polar','south_polar']
    n_re = n_elements(region_name_all)   &   ny = n_elements(year_all)
    FOR i_re = 0, n_re-1 DO BEGIN   
        IF region_name_all(i_re) EQ 'north_lobe' THEN  region_flag=x_gse LT -5 AND z_gsm GT 0 AND beta LT 0.05
        IF region_name_all(i_re) EQ 'south_lobe' THEN  region_flag=x_gse LT -5 AND z_gsm LT 0 AND beta LT 0.05 
        IF region_name_all(i_re) EQ 'north_polar' THEN  region_flag=x_gse GT -5 AND z_gsm GT 0 and beta lt 0.05
        IF region_name_all(i_re) EQ 'south_polar' THEN  region_flag=x_gse GT -5 AND z_gsm LT 0 and beta lt 0.05
        
        if keyword_set(eflux_filter) then filter_flag= (eflux_tail ge eflux_threshold or eflux_earth ge eflux_threshold or flag eq 0) else filter_flag = fltarr(n_elements(flag))+1
        index=where(region_flag and filter_flag and ABS(flag) ge 0,ct)
        
        if keyword_set(eflux_filter) then region_name= region_name_all(i_re)+'_and_eflux_ge_'+eflux_threshold_str else  region_name= region_name_all(i_re) 
        IF ct GT 0 THEN region_name_all(i_re) = data_distribution(year(index), flag=flag(index), storm_phase=storm_phase(index), xgrid = 1, para_name = inst_name, region_name = region_name, xrange = [2001, 2010], yrange = [0.01, 1],ylog=1,write_data = 0, events_distribution = 0, ratio_distribution = 1, path = path_year, ps = ps_plot,ratio_plot_range=[0.0001,1])
;        get_data,region_name_all(i_re)+'_ratio_nonstorm',data=data
;        store_data, name2(i),data={x:time_double(['2001-6','2002-6','2003-6','2004-6','2005-6','2006-6','2007-6','2008-6','2009-6']),y:data.y,dy:data.dy},lim={ylog:1,yrange:[0.001,1]},dlim={ytitle:region_name+'!C!C'+inst_name}
;        get_data,region_name_all(i_re)+'_ratio_storm',data=data
;        store_data, name2(i),data={x:time_double(['2001-6','2002-6','2003-6','2004-6','2005-6','2006-6','2007-6','2008-6','2009-6']),y:data.y,dy:data.dy},lim={ylog:1,yrange:[0.001,1]},dlim={ytitle:region_name+'!C!C'+inst_name}
        stop
    ENDFOR
    store_data, inst_name+'_year', data = {x:time, y:year}
                                ;   tplot_names, 'HIA*', names = name1
    tplot_names, 'CODIF*', names = name2
;    tplot_save, [name1], filename = path_year+'year_efficience'
    tplot_save, [name2], filename = path_year+'year_efficience'

    spawn, 'mogrify -format png '+path_year+'/*.ps'
    spawn, 'mogrify -rotate -90 '+path_year+'/*.png'    
    spawn, 'mv *.png '+path_year+'png/'
endfor 
ENDIF     

;---------------------------------------
;calc other distributions and plot them
;---------------------------------------
IF keyword_set(sort_distribution) THEN BEGIN 
    sta_path = path+'plots/sta/'
    spawn, 'mkdir '+sta_path
    sta_path = sta_path+time_str+'/'
    spawn, 'mkdir '+sta_path

    get_data, 'P', data = pre_data & time = pre_data.x & sw_p_pre = pre_data.y
    get_data, 'V', data = pre_data & sw_v_pre = pre_data.y
    get_data, 'B', data = pre_data & imf_b_pre = pre_data.y
    get_data, 'storm_phase', data = pre_data & storm_phase = pre_data.y
    get_data, 'flag', data = pre_data & flag = pre_data.y
    get_data, 'By', data = pre_data & imf_by_pre = pre_data.y
    get_data, 'Bz', data = pre_data & imf_bz_pre = pre_data.y
    get_data, 'x_gse', data = pre_data & x_gse = pre_data.y
    get_data, 'KP_Index', data = pre_data & kp_index = pre_data.y
    get_data, 'AE_Index', data = pre_data & ae_index = pre_data.y

    mu = 1.2566e-6
    p_mag = ((imf_b_pre*1e-9)^2/2/mu)*1e9 & p_kine = sw_p_pre ;nPa
    beta_imf = p_kine/p_mag  &  p_total = p_kine+p_mag    
    clock_angle = ATAN(imf_By_pre, imf_Bz_pre)*180/3.1415926
    
    index = where(clock_angle LT 0, ct)
    IF ct GT 0 THEN clock_angle(index) = clock_angle(index)+360

    region_name_all = ['outflow','tail_lobe'] ;,'all_region','x_le_neg5']
    n_re = N_ELEMENTS(region_name_all)
    imf_by_name = STRARR(n_re) & imf_bz_name = STRARR(n_re) & sw_p_name = STRARR(n_re) & sw_v_name = STRARR(n_re) 
    imf_b_name = STRARR(n_re) & clockangle_name = STRARR(n_re) & pb_name = STRARR(n_re) & Ey_name=STRARR(n_re)
    kp_name = STRARR(n_re) & ae_name = STRARR(n_re)
    FOR i_re = 0, n_re-1 DO BEGIN
        IF region_name_all(i_re) EQ 'all_region' THEN  index = where(x_gse GE -20, ct)
        IF region_name_all(i_re) EQ 'outflow' THEN  index = where(x_gse GT -5 AND beta Le 0.05 AND ABS(z_gsm) ge 4, ct)
        IF region_name_all(i_re) EQ 'x_gt_neg5' THEN  index = where(x_gse GT -5, ct)
        IF region_name_all(i_re) EQ 'tail_lobe' THEN  index = where(x_gse LE -5 AND beta Le 0.05, ct)
        IF region_name_all(i_re) EQ 'x_le_neg5' THEN  index = where(x_gse LE -5, ct)

        IF ct GT 0 THEN BEGIN 
            region_name = region_name_all(i_re) 
;            ae_name(i_re) = data_distribution(ae_index(index), flag=flag(index), storm_phase=storm_phase(index), xgrid = 100, para_name = 'AE_Index', region_name = region_name, xrange = [0, 1000],yrange = [0, 1], write_data = 0, events_distribution = 1, ratio_distribution = 1, path = sta_path, ps = ps_plot)
;            kp_name(i_re) = data_distribution(kp_index(index), flag=flag(index), storm_phase=storm_phase(index), xgrid = 1, para_name = 'KP_Index', region_name = region_name, xrange = [0, 10], yrange = [0, 1],write_data = 0, events_distribution = 1, ratio_distribution = 1, path = sta_path, ps = ps_plot)
;            imf_by_name(i_re) = data_distribution(imf_by_pre(index), flag=flag(index), storm_phase=storm_phase(index), xgrid = 1, para_name = 'IMF_By', region_name = region_name, xrange = [-10, 10],yrange = [0, 1], write_data = 0, events_distribution = 0, ratio_distribution = 1, path = sta_path, ps = ps_plot)
;           imf_bz_name(i_re) = data_distribution(imf_bz_pre(index), flag=flag(index), storm_phase=storm_phase(index), xgrid = 1, para_name = 'IMF_Bz', region_name = region_name, xrange = [-10, 10], yrange = [0, 1],write_data = 0, events_distribution = 0, ratio_distribution = 1, path = sta_path, ps = ps_plot)
;            sw_p_name(i_re) = data_distribution(sw_p_pre(index), flag=flag(index), storm_phase=storm_phase(index), xgrid = 1, para_name = 'SW_P', region_name = region_name, xrange = [0, 9],yrange = [0, 1], write_data = 0, events_distribution = 1, ratio_distribution = 1, path = sta_path, ps = ps_plot)
;            sw_v_name(i_re) = data_distribution(sw_v_pre(index), flag=flag(index), storm_phase=storm_phase(index), xgrid = 30, para_name = 'SW_V', region_name = region_name, xrange = [300, 600],yrange = [0, 1], write_data = 0, events_distribution = 1, ratio_distribution = 1, path = sta_path, ps = ps_plot)
;            imf_b_name(i_re) = data_distribution(imf_b_pre(index), flag=flag(index), storm_phase=storm_phase(index), xgrid = 1, para_name = 'IMF_B', region_name = region_name, xrange = [0, 14], yrange = [0, 1],write_data = 0, events_distribution = 0, ratio_distribution = 1, path = sta_path, ps = ps_plot)
            clockangle_name(i_re) = data_distribution(clock_angle(index), flag=flag(index), storm_phase=storm_phase(index), xgrid = 30, para_name = 'clock_angle', region_name = region_name, write_data = 0, events_distribution = 0, ratio_distribution = 1, xrange = [0, 360], yrange = [0, 1],path = sta_path, ps = ps_plot) 
            
;            Ey_name(i_re) = data_distribution(-sw_v_pre(index)*(imf_bz_pre(index) <0)*(1e3*1e-9*1e3), flag=flag(index), storm_phase=storm_phase(index), xgrid = 0.5, para_name = 'Ey', region_name = region_name, write_data = 0, events_distribution = 0, ratio_distribution = 1, xrange = [0, 6], yrange = [0, 1],path = sta_path, ps = ps_plot)         
;            v_tail=data(*,)
;            flux_tail=data(*,)
;            v_en_tail=mass_o*(v_tail),
;            sta_name = data_distribution(v_en_tail(index), flag(index), storm_phase(index), data_y = flux_tail(index), xgrid = 0.01, ygrid = 1)
;            plot_3d_distribution, sta_name, dim = 3, xtitle = 'flux', ytitle = 'en', region = region_name
            IF keyword_set(ps_plot) THEN   begin 
                spawn, 'mogrify -format png '+sta_path+'*.ps'  
                spawn, 'mogrify -rotate -90 '+sta_path+'*.png' 
                spawn, 'mkdir '+sta_path+'png/'
                spawn, 'mkdir '+sta_path+'events/'
                spawn, 'mv -f '+sta_path+'*events*.p*'+sta_path+'events/'
                spawn, 'mv -f '+sta_path+'*.png '+sta_path+'png/'
            endif 
        ENDIF  
    ENDFOR    
    plot_2d_distribution=0
    if keyword_set(plot_2d_distribution) then begin 
        print, plot_2d_distribution_diff_region(imf_by_name, xrange = [-10, 10], yrange = [0, 1], path = sta_path, ps = 1)
        print, plot_2d_distribution_diff_region(imf_bz_name, xrange = [-10, 10], yrange = [0, 1], path = sta_path, ps = 1)
        print, plot_2d_distribution_diff_region(sw_p_name, xrange = [0, 12], yrange = [0, 1], path = sta_path, ps = 1)
        print, plot_2d_distribution_diff_region(imf_b_name, xrange = [0, 15], yrange = [0, 1], path = sta_path, ps = 1)
;  print, plot_2d_distribution_diff_region(pb_name, xrange = [0, 0.3], yrange = [0.01, 1], path = path, ps = 1)
        print, plot_2d_distribution_diff_region(clockangle_name, xrange = [0, 360], yrange = [0, 1], path = sta_path, ps = 1)
    endif   
;    spawn, 'mogrify -format png '+sta_path+'*.ps'
                                ;   spawn, 'mogrify -rotate -90 '+sta_path+'*.png'
ENDIF 

if keyword_set(en_vs_distfunc) then begin  
    mean_value=0
    barplot=0
    histo_distfunc=0
    globe_plot=0
    shift=0
    no_outflows=0
    no_beam=0
    en_type=1
    vxb_study=1
    Vpar_study=1
    Vpar_study2=1
    Vpar_study3=1
  
    if en_type eq 1 then cusp_no_exb=1 else cusp_no_exb=0
    if en_vs_distfunc eq 1 then unit_name='DIFF FLUX'
    if en_vs_distfunc eq 2 then unit_name='DIST FUNC'
    if en_vs_distfunc eq 3 then unit_name='FLUX'
    if keyword_set(ps_plot) then thickness = 6 else thickness=3
    en_set = DOUBLE([31444.7, 19398.3, 11966.9, 7382.39, 4554.22, 2809.51, 1733.19, 1069.21, 659.599, 406.909, 251.023, 154.857, 95.5315, 58.9337, 36.3563])
    den_set = ([14900.6, 9192.19,5670.69,3498.26,2158.09,1331.33,821.301,506.663,312.562,192.820,118.951,73.3813,45.2691,27.9267,17.2280,5.63791])
    erange_set=[[en_set-den_set/2],[en_set+den_set/2]]

    distfunc_tail = data(*,59)
    distfunc_earth = data(*,60)
    diffflux_tail=data(*,7)
    diffflux_earth=data(*,16)
    flux_tail=data(*,8)*ABS(data(*,10))*1e5 ;cm-2s-1
    flux_earth=data(*,17)*ABS(data(*,19))*1e5 ;cm-2s-1
    mag = sqrt(data(*,26)^2+data(*,27)^2+data(*,28)^2)
    Vpar_tail=ABS(data(*,10))
    Vpar_earth=ABS(data(*,19))
    V_perp_tail=data(*,11)
    V_perp_earth=data(*,20)
    proton_Vx=data(*,30)
    proton_Vy=data(*,31)
    proton_Vz=data(*,32)
    v_tail=data(*,9)
    v_earth=data(*,18)
    bz=data(*,45)
    if not keyword_set(en_type) then en_type=0 
    case en_type of 
        0: begin
            energy_tail = data(*,5)
            energy_earth = data(*,14)
            en_type_str='all'
        end
        1: begin
            energy_tail=mass_o*(Vpar_tail)^2./2
            energy_earth=mass_o*(Vpar_earth)^2./2
            en_type_str='par'
        end
        2: begin
            energy_tail=mass_o*(V_perp_tail)^2./2.
            energy_earth=mass_o*(V_perp_earth)^2./2.
            en_type_str='perp'
        end
    endcase
    if keyword_set(eflux_filter) then begin
        index=where(eflux_tail lt eflux_threshold or eflux_earth lt eflux_threshold or flag eq 0,ct)
        if ct gt 0 then flag(index)=!values.f_nan
    endif 
    if unit_name eq 'DIST FUNC' THEN begin 
        cusp_outflow_max = [2.19618e-14,  6.38528e-14 , 1.83410e-13,  5.19320e-13 , 1.59421e-12 , 3.88443e-12,1.69021e-11,  6.10463e-11 , 1.54715e-10,  3.96661e-10 , 1.08688e-09 , 7.87800e-09,2.19171e-08,  3.07271e-08 , 1.16339e-08]*mag_normal
        cusp_outflow_median = [6.89262e-15,  3.46449e-14,  1.15268e-13,  3.13622e-13,  1.94826e-13 , 1.20396e-12,7.66728e-12 , 1.37056e-11,  3.87153e-11,  1.01449e-10 , 1.80336e-10,  3.34962e-10,3.40584e-10,  3.71697e-10,  1.67845e-10]*mag_normal
        lobe_beam_median=[4.65407e-16,  1.48719e-15 , 2.03593e-15 , 1.55242e-14 , 4.32470e-14,1.20995e-13 , 3.07697e-13 , 5.27818e-13,  6.66293e-12,  4.83978e-11, 1.74588e-10 , 2.70604e-10 , 4.34590e-10,  9.21287e-10,  1.14907e-09]*mag_normal
    endif 
    if unit_name eq 'DIFF FLUX' then begin 
        cusp_outflow_max=[166.830, 299.229,530.227, 926.172,1753.96, 2636.43, 7076.96, 15768.2,24653.1,38992.1,65910.4,294717.,505812.,437466.,102180.]*mag_normal
        cusp_outflow_median=[52.3590,162.354, 333.232, 559.324,214.348, 817.147,3210.32,3540.13,6169.11,9972.46,10935.9,12531.0,7860.14,5291.90,1474.17]*mag_normal
        lobe_beam_median=[!values.f_nan, !values.f_nan,!values.f_nan,!values.f_nan,!values.f_nan,!values.f_nan,!values.f_nan ,7354.56, 18207.0,27461.1,33464.7, 41585.2, 69578.1,102847,37763.6]*mag_normal
;    j0=[0.00000,0.00000,0.00000,0.00000 ,0.00000,0.00000 , 0.00000, 0.00000  , 0.00000 , 0.00000,0.0195301,21788.2 , 9.58386e+06 , 4.41944e+08 , 2.41601e+09]
    endif
    if unit_name eq 'FLUX' then begin 
                                ;     cusp_outflow_max=[!values.f_nan]
                                ;    cusp_outflow_median=[!values.f_nan]
                                ;   lobe_beam_median=[!values.f_nan, !values.f_nan,!values.f_nan,!values.f_nan,!values.f_nan,!values.f_nan,!values.f_nan]
                                ;   j0=[0.00000,0.00000,0.00000,0.00000 ,0.00000,0.00000 , 0.00000, 0.00000  , 0.00000 , 0.00000,0.0195301,21788.2 , 9.58386e+06 , 4.41944e+08 , 2.41601e+09]
    endif 
    if unit_name eq 'DIST FUNC' then begin 
        distfunc_tail_normal = distfunc_tail*mag_normal/mag
        distfunc_earth_normal = distfunc_earth*mag_normal/mag
                                ;       distfunc_tail_normal = (diffflux_tail*(0.16^2)/2e5/energy_tail)*mag_normal/mag/8
                                ;      distfunc_earth_normal = (diffflux_earth*(0.16^2)/2e5/energy_earth)*mag_normal/mag/8
                                ; !!!
        ytitle_input='Distribution Function (s!E3!N/cm!E3!N-km!E3!N)'
        yrange_input=[1e-9,1e-2] ;   [1e-10,1e-4]*mag_normal/100.
    endif
    if unit_name eq 'DIFF FLUX' then begin 
        distfunc_tail_normal = diffflux_tail ;*mag_normal/mag
        distfunc_earth_normal = diffflux_earth ;*mag_normal/mag
        ytitle_input='Differential Flux (1/cm!U2!N-s-sr-(eV/e))'
        yrange_input=[1e-1,1e5] ;*mag_normal/100.
    endif 
    if unit_name eq 'FLUX' then begin
        distfunc_tail_normal = flux_tail*mag_normal/mag
        distfunc_earth_normal = flux_earth*mag_normal/mag
        ytitle_input='Flux (cm!U-2!Ns!U-1!N)'
        yrange_input=[1e-5,1e10]*mag_normal/100.
    endif 
    if keyword_set(barplot) then begin 
        yrange_input=[0.,4100.] & ylog_input=0 & ytitle_input='# of O!U+!N beams'
    endif else ylog_input=1
    if keyword_set(shift)  then begin
        yrange_input=[1e-9,1e-2] & ylog_input=1
    endif 
    xlog_input=1
    xrange_input=[30,4e4]

    region_str = ['tail'];,'tail_lobe','polar_cap'] ;['tail_bl'];['polar_cap','tail_lobe','tail_bl'];,'ps','bl']
    storm_phase_set=['nonstorm','storm']
;set up strong_beam_ratio to record O+ beams that have higher distfunc
;than cusp outflows. 3 dimentions for beams that are higher, beams
;that are lower and the ratio of stronger beams over all beams and the error
    if keyword_set(number_ratio) then strong_beam_ratio=dblarr(15,4,n_elements(region_str),n_elements(storm_phase_set))
    for istorm=0,n_elements(storm_phase_set)-1 do begin 
        if storm_phase_set(istorm) eq 'nonstorm' then storm_flag = storm_phase eq 0 or storm_phase eq 5
        if storm_phase_set(istorm) eq 'storm' then storm_flag = storm_phase ge 1 and storm_phase le 3

        for ire=0,n_elements(region_str)-1 do begin
            if region_str(ire) eq 'lobe' then begin
                beta_flag = beta lt 0.02 & color_input=2
            endif 
            if region_str(ire) eq 'bl' then begin 
                beta_flag = beta gt 0.05 and beta lt 0.5 & color_input=1 ;PSBL has a different definition here with beta
            endif 
            if region_str(ire) eq 'ps' then beta_flag = beta gt 1
            if region_str(ire) eq 'tail_lobe' then begin
                beta_flag = beta lt 0.02 and x_gse lt -5 & color_input=2
            endif 
            if region_str(ire) eq 'polar_cap' then begin 
                beta_flag = beta lt 0.02 and x_gse gt -5 & color_input=3
            endif 
            if region_str(ire) eq 'tail_bl' then begin 
                beta_flag = beta gt 0.05 and beta lt 0.5 and x_gse lt -10 & color_input=1 ;PSBL has a different definition here with beta
            endif 
            if region_str(ire) eq 'tail_ps' then begin 
                beta_flag = beta gt 1 and x_gse lt -10 & color_input=4
            endif 
            if region_str(ire) eq 'tail' then begin
                beta_flag= (x_gse lt -5 and beta lt 0.05) or (x_gse lt -10 and beta gt 0.05) & color_input=6
            endif 

            ind_tail = where((flag eq 1 or flag eq 2) and beta_flag eq 1 and storm_flag eq 1,ct_tail) ; and (ABS(time-time_double('2002-09-11/9')) lt 1.*3600),ct_tail)
            ind_earth = where((flag eq -1 or flag eq 2) and beta_flag eq 1 and storm_flag eq 1,ct_earth) ; and (ABS(time-time_double('2002-09-11/9')) lt 1.*3600) ,ct_earth)
            plot_type=''
            if keyword_set(barplot) then plot_type='barplot' 
            if keyword_set(mean_value) then plot_type='mean'
            if keyword_set(histo_distfunc) then plot_type='histo_distfunc'
            if keyword_set(globe_plot) then plot_type='globe_plot'
            if keyword_set(shift) then plot_type='shift'
            
            if keyword_set(ps_plot) then begin
                plot_path = path+'plots/en_distfunc/'
                spawn, 'mkdir ' + plot_path
                if keyword_set(eflux_filter) then begin 
                    plot_path = plot_path+'eflux_gt_'+eflux_threshold_str+'/' 
                    spawn, 'mkdir ' + plot_path
                endif 
                plot_path = plot_path + time_str + '/'
                spawn, 'mkdir ' + plot_path
                plot_path = plot_path + strcompress(unit_name,/remove_all)+'/'
                spawn, 'mkdir ' + plot_path
                if not keyword_set(histo_distfunc) and not keyword_set(globe_plot) and not keyword_set(vxb_study) and not keyword_set(Vpar_study) then $
                  popen, plot_path+'en'+en_type_str+'_distfunc_'+region_str(ire)+'_'+storm_phase_set(istorm)+'_'+plot_type+'.ps',/land
            endif              
;            ylog_input=1
;            yrange_input= [1e-9,1e-2] ; [1e-10,1e-3]
                                ;     xrange_input=[500,4e4]
            if not keyword_set(histo_distfunc) and not keyword_set(globe) and not keyword_set(vxb_study) and not keyword_set(Vpar_study) then $
              plot,[0.,0.],[0.,0.],xlog=xlog_input,ylog=ylog_input,xrange=xrange_input,yrange=yrange_input,xstyle=1,ystyle=1,charsize=1.5,yminor=10,$
              title=region_str(ire)+'    '+storm_phase_set(istorm),xtitle='energy(eV)',ytitle=ytitle_input,/nodata
            
            dist_en,J0,unit=en_vs_distfunc,other_inputs=0,storm_phase=storm_phase_cusp,vperp=vperp
;apply storm_phase, so as to only plot
            temp=j0
            for i=0, 14 do temp(i,*)=storm_phase_cusp
            storm_phase_cusp=temp
            if storm_phase_set(istorm) eq 'nonstorm' then j0=j0*(storm_phase_cusp eq 0 or storm_phase_cusp eq 5)
            if storm_phase_set(istorm) eq 'storm' then j0=j0*(storm_phase_cusp ge 1 and storm_phase_cusp le 3)
            if storm_phase_set(istorm) eq 'main_phase' then j0=j0*(storm_phase_cusp eq 2)
            if storm_phase_set(istorm) eq 'recovery_phase' then j0=j0*(storm_phase_cusp eq 3)
            j0(where(j0 eq 0)) = !values.D_nan
            if keyword_set( cusp_no_exb) then begin
                en_set_new=DBLARR(n_elements(j0(*,0)),n_elements(j0(0,*)))
                for ievent=0, n_elements(j0(0,*))-1 do en_set_new(*,ievent)=en_set-mass_o*vperp(ievent)^2/2
                j0_median=DBLARR(n_elements(j0(*,0)))
                j0_up_median=DBLARR(n_elements(j0(*,0)))
                j0_down_median=DBLARR(n_elements(j0(*,0)))
                j0_max = DBLARR(n_elements(j0(*,0)))
                j0_error= DBLARR(n_elements(j0(*,0)))
                j0_mean= DBLARR(n_elements(j0(*,0)))
                for ien=0,n_elements(j0(*,0))-1 do begin
                    index=where(en_set_new gt erange_set(ien,0) and en_set_new lt erange_set(ien,1),ct)
                    if ct gt 0 then begin
                        j0_median(ien)=median(j0(index))
                        j0_max(ien)=max(j0(index),/nan)
                        j0_mean(ien)=mean(j0(index),/nan)
                        j0_error=sdom(j0(index))
                    endif 
                    index=where(en_set_new gt erange_set(ien,0) and en_set_new lt erange_set(ien,1) and j0 gt j0_median(ien),ct)
                    if ct gt 0 then j0_up_median(ien)=median(j0(index))
                    
                    index=where(en_set_new gt erange_set(ien,0) and en_set_new lt erange_set(ien,1) and j0 lt j0_median(ien),ct)
                    if ct gt 0 then j0_down_median(ien)=median(j0(index))
                    
                endfor 
            endif else begin 
                j0_median = median(j0,dim=2,/even)
                j0_up_median=REPLICATE(!VALUES.D_NAN,n_elements(j0_median))
                j0_down_median=REPLICATE(!VALUES.D_NAN,n_elements(j0_median))
                for ien=0,n_elements(j0_median)-1 do begin
                    index=where(j0(ien,*) gt j0_median(ien),ct)
                    if ct gt 0 then j0_up_median(ien)=median(j0(ien,index),/even)
                    index=where(j0(ien,*) le j0_median(ien),ct)
                    if ct gt 0 then j0_down_median(ien)=median(j0(ien,index),/even)
                endfor 
                j0_max = max(j0,dim=2,/nan)
                j0_error = DBLARR(n_elements(j0(*,0)))
                j0_mean = DBLARR(n_elements(j0(*,0)))
                for ie=0, n_elements(j0_max)-1 do begin 
                    j0_mean(ie)=mean(j0(ie,*),/nan)
                    j0_error(ie) = sdom(j0(ie,*))
                endfor  
            endelse 
            if plot_type eq '' or plot_type eq 'mean' or plot_type eq 'shift' then begin 
                if keyword_set(other_inputs) then begin 
                    oplot,en_set,j0(*,0),color=1,thick=4
                    xyouts, 3e2,1e-2,'outflow, calc, Seki',color=1,charsize=2
;   oplot,en_set,j0(*,1),color=3,thick=4
;     xyouts, 3e2,0.3e-2,'outflow, calc, Su 5000km',color=3,charsize=2
                    oplot,en_set,j0(*,2),color=4,thick=4
                    xyouts, 3e2,0.3e-2,'outflow calc, Su 8Re',color=4,charsize=2
                    oplot,en_set,j0(*,3),color=5,thick=4
                    xyouts, 3e2,1e-3,'outflow calc, Abe 6000km-9000km',color=5,charsize=2
                    oplot,en_set,j0(*,4),color=3,thick=4
                    xyouts, 3e2,0.3e-2,'outflow calc, Bouhram, 3.5-5Re',color=3,charsize=2
                endif  else begin 
                    if keyword_set(mean_value) then begin
                        oplot,en_set,j0_mean,thick=thickness,color=4
                        errplot,en_set,j0_mean+j0_error,j0_mean-j0_error,thick=thickness,color=4
                        xyouts, 1e3,1e-4,'cusp outflow mean value',charsize=2,color=245
                    endif else begin
                        if not keyword_set(no_outflows) then begin 
                            if en_vs_distfunc eq 2 and not keyword_set(shift) and not keyword_set(vxb_study) and not keyword_set(Vpar_study) then begin
                                if keyword_set(cusp_no_EXB) then begin
                                    for input=0,n_elements(j0(0,*))-1 do oplot,en_set_new(0:8,input),j0(0:8,input),color=230
                                endif else begin
                                    for input=0,n_elements(j0(0,*))-1 do oplot,en_set(0:8),j0(0:8,input),color=230
                                endelse
                            endif
                            oplot,en_set(0:8),j0_median(0:8),color=245,thick=thickness,psym=-2
                            xyouts, 1e3,1e-4,'cusp outflow median value',color=245,charsize=2
                        endif  
                        if keyword_set(shift) then begin
                ;            errplot,en_set,j0_down_median,j0_up_median,color=245,thick=thickness                            
;    oplot,en_set,j0_up_median,color=1,thick=thickness,psym=-2
;   oplot,en_set,j0_down_median,color=1,thick=thickness,psym=-2
                        endif
                    endelse  
                endelse  
                if not keyword_set(mean_value) and not keyword_set(shift) and not keyword_set(no_beam) then begin 
                    if not keyword_set(vxb_study) and not keyword_set(vxb_study) then begin 
                        if ct_tail gt 0 then oplot,energy_tail(ind_tail),distfunc_tail_normal(ind_tail),psym=1,color=color_input,symsize=0.5,thick=thickness
                        if ct_earth gt 0 then oplot,energy_earth(ind_earth),distfunc_earth_normal(ind_earth),psym=1,color=color_input,symsize=0.5,thick=thickness
                    endif 
;enpar_tail=mass_o*(vpar_tail)^2/2
;plot,energy_tail(ind_tail),enpar_tail(ind_tail),psym=3,xlog=1,ylog=1
;oplot,[1,1e6],[1,1e6],color=1
;----------- plot VXB study -----
                    if keyword_set(VXB_study) then begin
                        old=!p.multi
                        !p.multi=[0,2,2]
                        Xrange_input=[1e-5,1e3]
                        if keyword_set(ps_plot) then popen, path+'plots/beta_dependence_'+region_str(ire)+'.ps',/land
                        plot,beta(ind_tail),B(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Plasma Beta',ytitle='B (nT)',xrange=xrange_input,xstyle=1,ystyle=1,yrange=[1,1000],charsize=1.2,title=region_str(ire)
                        oplot,[0.02,0.02],[0.001,1000]
                        oplot,[0.05,0.05],[0.001,1000]
                        oplot,[0.5,0.5],[0.001,1000]
                        oplot,[1,1],[0.001,1000]
                        median_line,beta(ind_tail),B(ind_tail),xrange_input,result,xlog=1
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;---------------------
                        plot,beta(ind_tail),v_perp_tail(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Plasma Beta',ytitle='Vperp (km/s)',xrange=xrange_input,yrange=[1,1e3],xstyle=1,ystyle=1,charsize=1.2
                        oplot,[0.02,0.02],[0.001,1000]
                        oplot,[0.5,0.5],[0.001,1000]
                        oplot,[0.05,0.05],[0.001,1000]
                        oplot,[1,1],[0.001,1000]
                        median_line,beta(ind_tail),v_perp_tail(ind_tail),xrange_input,result,xlog=1
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;---------------------------
                        plot,beta(ind_tail),e_t(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Plasma Beta',ytitle='E (mV/m)',xrange=xrange_input,xstyle=1,ystyle=1,yrange=[0.01,10],charsize=1.2
                        oplot,[0.02,0.02],[0.001,100]
                        oplot,[0.5,0.5],[0.001,100]
                        oplot,[0.05,0.05],[0.001,1000]
                        oplot,[1,1],[0.001,1000]
                        median_line,beta(ind_tail),e_t(ind_tail),xrange_input,result,xlog=1
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;                        plot,beta(ind_tail),exb_t(ind_tail),xlog=1,ylog=1,psym=3,xtitle='beta',ytitle='EXB',xrange=xrange_input,yrange=[0.1,1000],xstyle=1,ystyle=1
;                        oplot,[0.02,0.02],[0.001,1000]
;                        oplot,[0.5,0.5],[0.001,1000]
;                        oplot,[0.2,0.2],[0.001,1000]
;                        oplot,[1,1],[0.001,1000]
                        
;                        plot,beta(ind_tail),ABS(vpar_tail(ind_tail)),xlog=1,ylog=1,psym=3,xtitle='beta',ytitle='Vpar',xrange=xrange_input,yrange=[1,1000],xstyle=1,ystyle=1,charsize=1.2
                        plot,beta(ind_tail),e_t(ind_tail)*sqrt(mag_normal/b(ind_tail)),xlog=1,ylog=1,psym=3,xtitle='Plasma Beta',ytitle='Normalized E (mV/m)',xrange=xrange_input,yrange=[1,1000],xstyle=1,ystyle=1,charsize=1.2
;                        plot,beta(ind_tail),exb_t(ind_tail),xlog=1,ylog=1,psym=3,xtitle='beta',ytitle='E(normal)',xrange=xrange_input,yrange=[1,1000],xstyle=1,ystyle=1,charsize=1.2
                        oplot,[0.02,0.02],[0.001,1000]
                        oplot,[0.5,0.5],[0.001,1000]
                        oplot,[0.05,0.05],[0.001,1000]
                        oplot,[1,1],[0.001,1000]
                        median_line,beta(ind_tail),e_t(ind_tail)*sqrt(mag_normal/b(ind_tail)),xrange_input,result,xlog=1
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
                        if keyword_set(ps_plot) then pclose else stop
                    endif
                endif 
;----------- plot Vpar study -----
                    if keyword_set(Vpar_study) then begin
                        old=!p.multi
                        !p.multi=[0,2,2]
                        Xrange_input=[10,500]
                        xgrid_input=0.1
                        if keyword_set(ps_plot) then popen, path+'plots/Vpar_dependence_'+region_str(ire)+'.ps',/land
                        plot,Vpar_tail(ind_tail),beta(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='O+ beta',xrange=xrange_input,xstyle=1,ystyle=1,yrange=[1e-4,10],charsize=1.2,title=region_str(ire)
;                        oplot,[0.02,0.02],[0.001,1000]
 ;                       oplot,[0.05,0.05],[0.001,1000]
  ;                      oplot,[0.5,0.5],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),beta(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;---------------------
mu=1.26e-6
pressure_B=B^2/2/mu*1e-9
                        plot,Vpar_tail(ind_tail),pressure_B(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='B pressure',xrange=xrange_input,yrange=[0.001,100],xstyle=1,ystyle=1,charsize=1.2
;                        oplot,[0.02,0.02],[0.001,1000]
 ;                       oplot,[0.5,0.5],[0.001,1000]
  ;                      oplot,[0.05,0.05],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),pressure_B(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;---------------------------
                        pressure_total=beta*pressure_B
                        plot,Vpar_tail(ind_tail),pressure_total(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='Plasma Pressure',xrange=xrange_input,xstyle=1,ystyle=1,yrange=[1e-5,1],charsize=1.2
;                        oplot,[0.02,0.02],[0.001,100]
 ;                       oplot,[0.5,0.5],[0.001,100]
  ;                      oplot,[0.05,0.05],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),pressure_total(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;                        plot,beta(ind_tail),exb_t(ind_tail),xlog=1,ylog=1,psym=3,xtitle='beta',ytitle='EXB',xrange=xrange_input,yrange=[0.1,1000],xstyle=1,ystyle=1
;                        oplot,[0.02,0.02],[0.001,1000]
;                        oplot,[0.5,0.5],[0.001,1000]
;                        oplot,[0.2,0.2],[0.001,1000]
;                        oplot,[1,1],[0.001,1000]
;-----------------------
                        plot,Vpar_tail(ind_tail),pressure_tail(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='O+ Pressure',xrange=xrange_input,yrange=[1e-6,1e-1],xstyle=1,ystyle=1,charsize=1.2
;                        oplot,[0.02,0.02],[0.001,1000]
 ;                       oplot,[0.5,0.5],[0.001,1000]
  ;                      oplot,[0.05,0.05],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),pressure_tail(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
                        if keyword_set(ps_plot) then pclose else stop
                    endif 
;----------- plot Vpar study -----
                    if keyword_set(Vpar_study2) then begin
                        old=!p.multi
                        !p.multi=[0,2,2]
                        Xrange_input=[10,500]
                        xgrid_input=0.1
                        if keyword_set(ps_plot) then popen, path+'plots/Vpar_dependence2_'+region_str(ire)+'.ps',/land
                        pressure_h=pressure_total-pressure_tail
                        plot,Vpar_tail(ind_tail),pressure_h(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='P H+',xrange=xrange_input,xstyle=1,ystyle=1,yrange=[1e-5,1],charsize=1.2,title=region_str(ire)
;                        oplot,[0.02,0.02],[0.001,1000]
 ;                       oplot,[0.05,0.05],[0.001,1000]
  ;                      oplot,[0.5,0.5],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),pressure_h(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;---------------------
                        plot,Vpar_tail(ind_tail),proton_n(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='n H+',xrange=xrange_input,yrange=[0.0001,10],xstyle=1,ystyle=1,charsize=1.2
;                        oplot,[0.02,0.02],[0.001,1000]
 ;                       oplot,[0.5,0.5],[0.001,1000]
  ;                      oplot,[0.05,0.05],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),proton_n(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;---------------------------
kb=1.38e-23
;1.16e4 K/eV is the temperature for a 1eV particle 
temperature_h=pressure_h/proton_n/kb*1e-8*1e-9/1.16e4
                        plot,Vpar_tail(ind_tail),temperature_h(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='T H+',xrange=xrange_input,xstyle=1,ystyle=1,yrange=[0.1,1e4],charsize=1.2
;                        oplot,[0.02,0.02],[0.001,100]
 ;                       oplot,[0.5,0.5],[0.001,100]
  ;                      oplot,[0.05,0.05],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),temperature_h(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;                        plot,beta(ind_tail),exb_t(ind_tail),xlog=1,ylog=1,psym=3,xtitle='beta',ytitle='EXB',xrange=xrange_input,yrange=[0.1,1000],xstyle=1,ystyle=1
;                        oplot,[0.02,0.02],[0.001,1000]
;                        oplot,[0.5,0.5],[0.001,1000]
;                        oplot,[0.2,0.2],[0.001,1000]
;                        oplot,[1,1],[0.001,1000]
;-----------------------
                        plot,Vpar_tail(ind_tail),temperature_tail(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='T O+',xrange=xrange_input,yrange=[0.1,1e4],xstyle=1,ystyle=1,charsize=1.2
;                        oplot,[0.02,0.02],[0.001,1000]
 ;                       oplot,[0.5,0.5],[0.001,1000]
  ;                      oplot,[0.05,0.05],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),temperature_tail(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
                        if keyword_set(ps_plot) then pclose else stop
                    endif 
 ;----------- plot Vpar study -----

                    if keyword_set(Vpar_study3) then begin
                        old=!p.multi
                        !p.multi=[0,2,2]
                        Xrange_input=[10,500]
                        xgrid_input=0.1
                        if keyword_set(ps_plot) then popen, path+'plots/Vpar_dependence3_'+region_str(ire)+'.ps',/land
                        pressure_h=pressure_total-pressure_tail
                        plot,Vpar_tail(ind_tail),density_tail(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='n O+',xrange=xrange_input,xstyle=1,ystyle=1,yrange=[1e-5,1],charsize=1.2,title=region_str(ire)
;                        oplot,[0.02,0.02],[0.001,1000]
 ;                       oplot,[0.05,0.05],[0.001,1000]
  ;                      oplot,[0.5,0.5],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),density_tail(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;---------------------
                        plot,proton_v(ind_tail),proton_n(ind_tail),xlog=1,ylog=1,psym=3,xtitle='V H+',ytitle='n H+',xrange=xrange_input,yrange=[0.0001,10],xstyle=1,ystyle=1,charsize=1.2
;                        oplot,[0.02,0.02],[0.001,1000]
 ;                       oplot,[0.5,0.5],[0.001,1000]
  ;                      oplot,[0.05,0.05],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,proton_v(ind_tail),proton_n(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;---------------------------
kb=1.38e-23
;1.16e4 K/eV is the temperature for a 1eV particle 
temperature_h=pressure_h/proton_n/kb*1e-8*1e-9/1.16e4
                        plot,Vpar_tail(ind_tail),temperature_h(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='T H+',xrange=xrange_input,xstyle=1,ystyle=1,yrange=[0.1,1e4],charsize=1.2
;                        oplot,[0.02,0.02],[0.001,100]
 ;                       oplot,[0.5,0.5],[0.001,100]
  ;                      oplot,[0.05,0.05],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),temperature_h(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
;                        plot,beta(ind_tail),exb_t(ind_tail),xlog=1,ylog=1,psym=3,xtitle='beta',ytitle='EXB',xrange=xrange_input,yrange=[0.1,1000],xstyle=1,ystyle=1
;                        oplot,[0.02,0.02],[0.001,1000]
;                        oplot,[0.5,0.5],[0.001,1000]
;                        oplot,[0.2,0.2],[0.001,1000]
;                        oplot,[1,1],[0.001,1000]
;-----------------------
                        plot,Vpar_tail(ind_tail),temperature_tail(ind_tail),xlog=1,ylog=1,psym=3,xtitle='Vpar',ytitle='T O+',xrange=xrange_input,yrange=[0.1,1e4],xstyle=1,ystyle=1,charsize=1.2
;                        oplot,[0.02,0.02],[0.001,1000]
 ;                       oplot,[0.5,0.5],[0.001,1000]
  ;                      oplot,[0.05,0.05],[0.001,1000]
   ;                     oplot,[1,1],[0.001,1000]
                        median_line,Vpar_tail(ind_tail),temperature_tail(ind_tail),xrange_input,result,xlog=1,xgrid=xgrid_input
                        oplot,result(*,0),result(*,1),psym=-2,color=1,thick=3
                        if keyword_set(ps_plot) then pclose else stop
                    endif 
                    if keyword_set(flux_threshold) then begin 
                    oplot, en_set,10*(0.16^2)/2e5/en_set/36.26*mag_normal/8,psym=-1,color=color_input,thick=4
                    xyouts,0.2e4,2e-6,'flux threshold',color=color_input,charsize=2 
                endif
            endif 
; calculate how many beams are stronger than the upper envelope of
; cusp outflows
            analysis_of_accelerated_beams=0
            if keyword_set(analysis_of_accelerated_beams) then begin
                for ien=0,8 do begin 
                    ind_tail=where(energy_tail gt erange_set(ien,0) and energy_tail lt erange_set(ien,1) $
                                   and (flag eq 1 or flag eq 2) and beta_flag eq 1 and storm_flag eq 1 $
                                   and distfunc_tail_normal gt j0_max(ien) $
                                   , ct_tail)
                    ind_earth=where(energy_earth gt erange_set(ien,0) and energy_earth lt erange_set(ien,1) $
                                    and (flag eq -1 or flag eq 2) and beta_flag eq 1 and storm_flag eq 1 $
                                    and distfunc_earth_normal gt j0_max(ien)$
                                    ,ct_earth)
                    print,en_set(ien),ct_tail,ct_earth
                    if ct_tail gt 0 then begin 
                        if keyword_set(ind_tail_acce) then ind_tail_acce=[ind_tail_acce,ind_tail] else ind_tail_acce=ind_tail
                    endif 
                    if ct_earth gt 0 then begin 
                        if keyword_set(ind_earth_acce) then ind_earth_acce=[ind_earth_acce,ind_earth] else ind_earth_acce=ind_earth 
                    endif 
                endfor 
                                ;  plot,proton_n(ind_tail_acce),density_tail(ind_tail_acce),psym=1,xrange=[0,1],yrange=[0,1]
                                ;         plot,proton_Vy(ind_tail_acce)
                 
            endif       
; calculate number, median,mean,max,min,sdom value of distfunc and save in distfunc_avg           
            distfunc_avg=DBLARR(15,8)
            for ien=0, 14 do begin
                ind_tail=where(energy_tail gt erange_set(ien,0) and energy_tail lt erange_set(ien,1) $
                               and (flag eq 1 or flag eq 2) and beta_flag eq 1 and storm_flag eq 1 $
                               , ct_tail)
                ind_earth=where(energy_earth gt erange_set(ien,0) and energy_earth lt erange_set(ien,1) $
                                and (flag eq -1 or flag eq 2) and beta_flag eq 1 and storm_flag eq 1 $
                                , ct_earth)
                
                if keyword_set(globe_plot) and (ct_tail+ct_earth) gt 0 then begin
                    if ct_tail gt 0 then ind_beam_tail=where(diffflux_tail(ind_tail) lt 10,nbeam_tail) else nbeam_tail=0
                    if ct_earth gt 0 then ind_beam_earth=where(diffflux_earth(ind_earth) lt 10,nbeam_earth) else nbeam_earth=0
                    nbeam=nbeam_tail+nbeam_earth
                    if nbeam gt 0 then begin 
                        for ibeam=0,nbeam-1 do begin 
                            if ibeam lt nbeam_tail then begin 
                                index=ind_tail(ind_beam_tail(ibeam))
                            endif else begin 
                                index=ind_earth(ind_beam_earth(ibeam-nbeam_tail))
                            endelse 

                            sat = sc
                            specie= 3
                            timespan, time(index)-150, 5, /min
                            units_name='DIFF FLUX' ;'Counts', 'NCOUNTS', 'RATE', 'NRATE', 'DIFF FLUX', 'EFLUX'
                            inst=0 &  eff_table=0 
                            
                            plot_globe_from_crib, sat, specie, inst,units_name, eff_table
                            name = 'GLOBE_SC'+string(sat,format='(i1.1)')+'_'+strcompress(units_name, /remove_all)  +'*'+'SP'+string(specie,format='(i1.1)')
                            tplot_names, name, names=gname
                            get_data, gname(0), data=data
                            plot3d_options, log = 1
                            if keyword_set(ps_plot) then begin
                                plot_path=path+'plots/en_distfunc/globe/'
                                spawn, 'mkdir '+plot_path
                                popen, plot_path+'en_distfunc_'+region_str(ire)+'_'+storm_phase_set(istorm)+'_'+plot_type+'_'+strcompress(fix(en_set(ien)),/remove_all)+'_'+strcompress(ibeam,/remove_all)+'.ps',/land
                            endif
                            print, en_set(ien),nbeam
                            if ien eq 0 then plot3d_codif, data, zrange=[0.1,10] , ebins = [0,1,2]
                            if ien eq 1 then plot3d_codif, data, zrange=[0.1,10] , ebins = [0,1,2,3]
                            if ien eq 14 then plot3d_codif, data, zrange=[0.1,10] , ebins = [12,13,14]
                            if ien eq 13 then plot3d_codif, data, zrange=[0.1,10] , ebins = [11,12,13,14]
                            if ien ge 2 and ien le 12 then plot3d_codif, data, zrange=[0.1,10] , ebins = [ien-2,ien-1,ien,ien+1,ien+2]
                            if keyword_set(ps_plot) then pclose else stop
                        endfor   
                    endif  
                endif   

                if keyword_set(histo_distfunc) and (ct_tail+ct_earth gt 0) then begin 
                    if keyword_set(ps_plot) then begin
                        plot_path=path+'plots/en_distfunc/histo_distfunc/'
                        spawn, 'mkdir '+plot_path
                        popen, plot_path+'en_distfunc_'+region_str(ire)+'_'+storm_phase_set(istorm)+'_'+plot_type+'_'+strcompress(fix(en_set(ien)),/remove_all)+'.ps',/land
                    endif
                    if ct_tail gt 0 and ct_earth gt 0 then $
                      a=data_distribution(alog10([distfunc_tail_normal(ind_tail),distfunc_earth_normal(ind_earth)]),flag=[flag(ind_tail),flag(ind_earth)],storm_phase=[storm_phase(ind_tail),storm_phase(ind_earth)],events_distribution=0,plot_single_phase_events_distribution=1,single_phase=storm_phase_set(istorm),xgrid=0.2,xrange=[-8,-2],region_name=strcompress(fix(en_set(ien)),/remove_all)+'eV_'+region_str(ire),para_name='Dist_Func',path=plot_path,rotate=90) ;,bar_yrange_input=[0,200])
                    if ct_tail gt 0 and ct_earth eq 0 then $
                      a=data_distribution(alog10([distfunc_tail_normal(ind_tail)]),flag=[flag(ind_tail)],storm_phase=[storm_phase(ind_tail)],events_distribution=0,plot_single_phase_events_distribution=1,single_phase=storm_phase_set(istorm),xgrid=0.2,xrange=[-8,-2],region_name=strcompress(fix(en_set(ien)),/remove_all)+'eV_'+region_str(ire),para_name='Dist_Func',path=plot_path,rotate=90) ;,bar_yrange_input=[0,200])
                    if ct_tail eq 0 and ct_earth gt 0 then $
                      a=data_distribution(alog10([distfunc_tail_normal(ind_earth)]),flag=[flag(ind_earth)],storm_phase=[storm_phase(ind_earth)],events_distribution=0,plot_single_phase_events_distribution=1,single_phase=storm_phase_set(istorm),xgrid=0.2,xrange=[-8,-2],region_name=strcompress(fix(en_set(ien)),/remove_all)+'eV_'+region_str(ire),para_name='Dist_Func',path=plot_path,rotate=90) ;,bar_yrange_input=[0,200])
                    if keyword_set(ps_plot) then pclose else stop
                endif  
                if ct_tail gt 0  or ct_earth gt 0 then begin
                    if ct_tail gt 0 and ct_earth gt 0 then Var = [distfunc_tail_normal(ind_tail),distfunc_earth_normal(ind_earth)]
                    if ct_tail gt 0 and ct_earth eq 0 then Var = distfunc_tail_normal(ind_tail )
                    if ct_tail eq 0 and ct_earth gt 0 then Var = distfunc_earth_normal(ind_earth )
                    if total(var,/nan) ne 0 then begin 
                        distfunc_avg(ien,0) = median(Var,/even)
                        distfunc_avg(ien,1) = mean(Var,/nan)
                        distfunc_avg(ien,2) = max(Var,/nan)
                        distfunc_avg(ien,3) = min(Var,/nan)
                        distfunc_avg(ien,4) = sdom(Var)
                        distfunc_avg(ien,5)= ct_tail+ct_earth
                        distfunc_avg(ien,6) = median(Var(where(Var le distfunc_avg(ien,0))),/even)
                        distfunc_avg(ien,7) = median(Var(where(Var ge distfunc_avg(ien,0))),/even)
                    endif else distfunc_avg(ien,*) = !VALUES.F_NAN
                endif else distfunc_avg(ien,*) = !VALUES.F_NAN
            endfor
            
            if  keyword_set(barplot) then begin
                bar_plot,n_elements(j0(0,*)),outline=1,/overplot,baroffset=0.56,color=4
                bar_plot,reverse(distfunc_avg(*,5)),outline=1,/overplot,color=reverse(distfunc_avg(*,5) gt 0)*color_input
            endif else begin 
                if en_vs_distfunc eq 2 then begin
                    if keyword_set(mean_value) then begin
                        oplot,en_set,distfunc_avg(*,1),psym=0,color=color_input,thick=thickness
                        xyouts, 1e3,1e-3,'O!U+!N beam mean value',color=color_input,charsize=2
                        errplot,en_set,distfunc_avg(*,1)+distfunc_avg(*,4),distfunc_avg(*,1)-distfunc_avg(*,4),color=color_input,thick=thickness
                    endif else begin  
                        oplot,en_set,distfunc_avg(*,0),psym=-2,color=6,thick=thickness
                                ;       calculated_density=(sqrt(9./10)*sqrt(1./20)*sqrt(1./20))*(1/mass_o)^1.5*Total(distfunc_avg(*,0)*(1/en_set^1.5)*(den_set^3),/nan) ;cm-3
                                ;      calculated_flux=(sqrt(9./10)*sqrt(1./20)*sqrt(1./20))*(sqrt(2)/mass_o^2)*Total(distfunc_avg(*,0)*(1/en_set)*(den_set^3),/nan)*double(1e5) ;cm-2s-1
                        calculated_density=(0.3*3.1415926/(2*mass_o)^1.5)*Total(distfunc_avg(*,0)*en_set^1.5,/nan)*(sqrt(2.)/2)
                        calculated_flux=(0.3*3.1415926/(2*mass_o)^2)*Total(distfunc_avg(*,0)*en_set^2,/nan)*double(1.e5)*(sqrt(2.)/2)
                        print, calculated_density, calculated_flux
                        xyouts, 1e3,1e-3,'O!U+!N beam median value',color=6,charsize=2
                  ;      xyouts, 1e3,5e-4,strcompress(calculated_density,/remove_all)+', '+strcompress(calculated_flux,/remove_all),color=6,charsize=2
                        if keyword_set(shift) then begin
                            errplot,en_set,distfunc_avg(*,6),distfunc_avg(*,7),thick=thickness,color=6
;    oplot,en_set,distfunc_avg(*,6),psym=-2,color=6,thick=thickness
;   oplot,en_set,distfunc_avg(*,7),psym=-2,color=6,thick=thickness
                        endif
                    endelse   
; max and min         
;   oplot,en_set,distfunc_avg(*,2),psym=0,color=4
; oplot,en_set,distfunc_avg(*,3),psym=0,color=5 
                                ; case value
                                ;        oplot, en_set,cusp_outflow_max,color=3,linestyle=5,thick=6
                                ;         xyouts, 5e2,4e-3,'outflow envelop from cases',color=3,charsize=2
                                ;       oplot, en_set,cusp_outflow_median,color=2,linestyle=5,thick=6
                                ;      xyouts, 5e2,9e-3,'outflow median from cases',color=2,charsize=2
                                ;    oplot, en_set,lobe_beam_median,color=3
                                ;     xyouts, 5e2,2.5e6,'lobe beam median from cases',color=3,charsize=2
                endif  
            endelse                
            
            if keyword_set(shift) then begin
                energy_shift=REPLICATE(!values.d_nan,n_elements(distfunc_avg(*,0)),3)
                compare_var=j0_median
                
                if storm_phase_set(istorm) eq 'nonstorm' then compare_var(9:14)= [3.3834524e-06,8.7104033e-06,1.9936867e-05,4.7596766e-05,0.00010271290,7.6125051e-05]
                if storm_phase_set(istorm) eq 'storm' then compare_var(9:14)=[6.3186379e-06,3.8432709e-05,0.00011813651,0.00033697372,0.00035912757,0.00013837380]
                x=sqrt(2*en_set(0:8)/mass_o) & y = j0_median(0:8)
                if storm_phase_set(istorm) eq 'storm' then estimates=[8e-5,50,90]
                if storm_phase_set(istorm) eq 'nonstorm' then estimates = [1e-7,-50,90]

                measure_error=(j0_up_median(0:8)-j0_down_median(0:8))/2
                yfit_para = GAUSSFIT(x, y, co, NTERMS = 3,estimates=estimates,yerror = yerror, sigma = sigma, chisq = chisq,measure_error=measure_error)
                
                                ;plot,en_set,y,xlog=1,ylog=0,xstyle=1,ystyle=1,xrange=[40,4e4],yrange=[1e-11,1e-5]
                oplot,en_set,co(0)*exp(-((sqrt(2*en_set/mass_o)-co(1))/co(2))^2/2),color=2
                for jj=0,2 do begin
                    if jj eq 0 then var=distfunc_avg(*,0)
                    if jj eq 1 then var=distfunc_avg(*,6)
                    if jj eq 2 then var=distfunc_avg(*,7)
                    for ien=0,n_elements(var)-1 do begin
                        if en_set(ien) gt 0 and compare_var(ien) lt var(ien) then begin 
                            ind1=where(compare_var le var(ien),ct1)
                            ind2=where(compare_var gt var(ien),ct2)
                            print,en_set(ien),var(ien)
                            if ct1 ne 0 and ct2 ne 0 then begin
                                x0=en_set(ien)  & y0=var(ien)
                                x1=en_set(ind1(ct1-1))  &  y1=compare_var(ind1(ct1-1))
                                x2=en_set(ind2(0))  &  y2=compare_var(ind2(0))
                                k=(y2-y1)/(x2-x1)  &  b=y1-k*x1
                                x=(y0-b)/k
                                energy_shift(ien,jj)=x0-x
                            endif  
                        endif else energy_shift(ien,jj)=0
                    endfor   
                endfor
                if keyword_set(ps_plot) then pclose else stop
                if keyword_set(ps_plot) then  popen, plot_path+'en_distfunc_'+region_str(ire)+'_'+storm_phase_set(istorm)+'_'+plot_type+'_dE.ps',/land

;                plot,[0.,0.],[0.,0.],xlog=1,ylog=0,xrange=[500,4e4],yrange=[0,1],xstyle=1,ystyle=1,charsize=1.5,$
 ;                 title=region_str(ire)+'    '+storm_phase_set(istorm),xtitle='E (eV)',ytitle='dE (eV)',/nodata
  ;              oplot,en_set,energy_shift(*,0)/en_set,psym=-1,thick=thickness
   ;             errplot,en_set,energy_shift(*,1)/en_set,energy_shift(*,2)/en_set,thick=thickness

                                   plot,[0.,0.],[0.,0.],xlog=1,ylog=1,xrange=[40,4e4],yrange=[1,3e4],xstyle=1,ystyle=1,charsize=1.5,$
                                     title=region_str(ire)+'    '+storm_phase_set(istorm),xtitle='E(eV)',ytitle='dE (eV)',/nodata
                                   oplot,en_set,energy_shift(*,0),psym=-7,thick=thickness
                                   errplot,en_set,energy_shift(*,1),energy_shift(*,2),thick=thickness
            endif  
            
;            x=[!values.f_NaN ,           !values.f_NaN  ,           !values.f_NaN       ,      !values.f_nan , 1.3585331e-08,   3.0600020e-08,$
;            1.2461401e-07 ,  2.9932620e-07  , 1.4845807e-06,   3.3834524e-06 ,  8.7104033e-06,   1.9936867e-05,$
;            4.7596766e-05 ,  0.00010271290 ,  7.6125051e-05]
                                ;NaN             NaN             NaN             NaN      0.21516875       1.0152053
                                ;     0.18268627      0.14990889      0.28781352      0.53547180      0.22664038      0.16876128
                                ;     0.14124770      0.28600672      0.55014061

;r=   [ NaN             NaN             NaN             NaN       4.6475151      0.98502248
;       5.4738651       6.6707186       3.4744720       1.8675120       4.4122765       5.9255299
;       7.0797611       3.4964214       1.8177171]
            if keyword_set(ps_plot) then pclose else stop
        endfor
    endfor
    
    if keyword_set(ps_plot) then begin 
        spawn,'mogrify -format png '+plot_path+'*.ps'
        spawn,'mogrify -rotate -90 '+plot_path+'*.png'
    endif 
endif
print,'Program End'     
stop   
END 
