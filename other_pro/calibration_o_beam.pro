PRO calibration_o_beam

sc = 4  & sc_str = STRING(sc, format = '(i1.1)')

time_start = '2006-01-01/00:00:00' 
time_end = '2008-01-01/00:00:00'

;keywords settings
data_read_from_dat = 1
save_data = 1

;---- basic settings ---
average_time = 300

;--- set the time info ---
ts = time_double(time_start)
te = time_double(time_end)
dt = te-ts

timespan, time_start, dt,  /SECONDS

ts_str = time_struct(ts)        ; start_time tplot time structure
te_str = time_struct(te)        ; end_time tplot time structure
jd_s = julday(ts_str.month, ts_str.date, ts_str.year) ;start julian day
jd_e = julday(te_str.month, te_str.date, te_str.year) ;end julian day

ndays = (jd_e - jd_s) + 1       ; number of days to be loaded
;Last day is not included if hour=min=sec=0
IF te_str.hour EQ 0 AND te_str.min EQ 0 AND te_str.sec EQ 0 THEN ndays = ndays - 1

;------------------------------------------
;check prestore data 
;------------------------------------
flndata = 'tplot_restore/2006_2007data'
print, FINDFILE(flndata+'.tplot', COUNT = ct_data)
IF KEYWORD_SET(data_read_from_dat) THEN ct_data = 0
IF ct_data GT 0 THEN BEGIN 
    tplot_restore, filenames = flndata+'.tplot' 
    get_data, '2006_2007data', data = data
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
               '       SW_P    '] ;47

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
        fln = 'data/'+'storm_o_beam_'+year_str+month_str+day_str  +'.dat'
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
        store_data, '2006_2007data', data = {x:time, y:data, title:title}
        tplot_save, filename = flndata
    ENDIF
ENDELSE  

;- all data = -32768. is default invalid data 
    index = where(data EQ -32768., ct)
    IF ct GT 0 THEN data(index) = !VALUES.F_NAN  

;- elimilate magnetosheath data by set the flag:data(*,0) into INF 
    IF KEYWORD_SET(no_magnetosheath) THEN BEGIN 
        h_density = data(*, 29)
        h_velocity = sqrt((data(*, 30))^2+(data(*, 31))^2 + (data(*, 32))^2)
        data(*, 0) = data(*, 0)/((h_density LE 3) OR (h_velocity LE 65)) 
    ENDIF     

;- elimilate solarwind data: define solarwind region as 
;x_gse > = -1 and outside an ellips deceided by observations 
    IF KEYWORD_SET(no_solarwind) THEN BEGIN 
        x_gse = data(*, 2)
        y_gse = data(*, 3)
        data(*, 0) = data(*, 0)/((x_gse LE -1) OR ((x_gse^2/8^2+(y_gse-1)^2/15^2) LE 1))
    ENDIF 

; save data into more meanful names
flag = data(*, 0)
distance = SQRT(data(*, 2)^2+data(*, 3)^2+data(*, 4)^2)
orbit_time = 57*3600.

petime = time_double('2005-12-29/23:22:30')
pedistance = -1.
npe = 0l
WHILE (petime(npe)+orbit_time) LT time_double('2008-01-01/00:00:00') DO BEGIN   
    index = WHERE(ABS(time - (petime(npe)+orbit_time)) LE 1.5*60.*60., ct)
    IF ct GT 0 THEN BEGIN 
        distance_temp = distance(index)
        time_temp = time(index)
        loc = SORT(distance_temp)
        
        petime = [petime, time_temp(loc(0))]
        pedistance = [pedistance, distance_temp(loc(0))]
        npe = npe+1
    ENDIF ELSE BEGIN 
        petime = [petime, petime(npe)+orbit_time]
        pedistance = [pedistance, -1]
        npe = npe+1
    ENDELSE 
ENDWHILE 

index = where(pedistance NE  -1, ct)
IF ct GT 0 THEN BEGIN 
    petime = petime(index)
    pedistance = pedistance(index)
    npe = ct
ENDIF 
;window, /free
;plot, petime, pedistance, psym = 3

o_beam_around_perigee = FLTARR(3000)
io = 0
FOR ipe = 0, npe-1 DO BEGIN 
    index = where(ABS(time - petime(ipe)+30*60.) GE 1.*60.*60. AND $
                  ABS(time - petime(ipe)+30*60.) LE 3.5*60.*60., ct)
    time_temp = time(index)
    a = flag(index) NE 0

    IF TOTAL(a, /NAN) GE 8 THEN BEGIN 
        o_beam_around_perigee(io) = petime(ipe)
        io = io+1
    ENDIF 
ENDFOR 
o_beam_around_perigee = o_beam_around_perigee(0:io-1)

print, time_string(o_beam_around_perigee)

;OPENU, unit, '2006_2007calibration_o_beam_time.txt', /GET_LUN, /APPEND
;PRINTF, unit, time_string(o_beam_around_perigee)+'                                              '
;FREE_LUN, unit      

stop
END 
