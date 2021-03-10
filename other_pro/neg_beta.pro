PRO neg_beta
sc = 4 
sc_str = STRING(sc, FORMAT = '(i1.1)')
average_time = 5 * 60 ;in seconds
at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL) 

OPENR, unit, 'neg_beta_day.dat', /GET_LUN
neg_beta_days = DBLARR(300) 
jj = 0l
dummy = ''
WHILE NOT EOF(unit) DO BEGIN
    READF, unit, dummy
    neg_beta_days(jj) = time_double(STRMID(dummy, 0, 20))
    jj = jj+1      
ENDWHILE
CLOSE, unit, /all
ndays = jj
neg_beta_days = neg_beta_days(0:ndays-1)

FOR i = 0, ndays-1 DO BEGIN 

    tplot_names, names = names
    store_data, DELETE = names

    timespan, neg_beta_days(i), 1, /days
    ts = time_string(neg_beta_days(i))
    te = time_string(neg_beta_days(i)+ 24.*3600.)
    date_s = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
    time_s = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
    date_e = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2)
    time_e = STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2)

;-- Load CLUSTER ephemeris--
    sat = [sc] 
    get_cluster_ephemeris, sat, /GSE_X, /GSE_Y, /GSE_Z, /DIST, /MLT $
      , /GSM_X, /GSM_Y, /GSM_Z 
    
;-- Load CLUSTER Magnetic field--
    sat = sc
    plot_mag_from_crib, sat, POLAR = 1, GSM = 1

;-- Load CLUSTER H+ and O+ moments--
    sat = [sc, sc]
    specie = [0, 3]
    moments = ['A', 'A']
    angle = [[-90, 90], [0, 360]]
    energy = [40., 40000.]
    inst = 0 & eff_table = 0
    
    plot_3dmom_from_crib, sat, specie, inst, $
      moments, angle, energy, eff_table, recalc = 0

;---------Check weather data required to calculate beta are loaded 
;Including : H pressure, O pressure and mag pressure 
;nerror = 0 = > no problem to calculate beta      ------------
    h_press =  'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP0_ET0_All'
    o_press = 'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP3_ET0_All'
    mag_press = 'MAG_SC' + sc_str +'_B_xyz_gse_MAG_PR'      

    s_e = STRARR(10)
    nerror = 0
    
    tplot_names, h_press, names = names
    IF names(0) EQ '' THEN BEGIN 
        s_e(nerror) = ' H_pressure,'
        nerror = nerror+1
    ENDIF  ELSE BEGIN 
        get_data, names(0), data = data
        IF N_ELEMENTS(data.x) LT 2 THEN BEGIN 
            s_e(nerror) = 'H_pressure not enough data'
            nerror = nerror+1
        ENDIF 
    ENDELSE 
    
    tplot_names,o_press, names = names
    IF names(0) EQ '' THEN BEGIN 
        s_e(nerror) = ' O_pressure,'
        nerror = nerror+1
    ENDIF ELSE BEGIN 
        get_data, names(0), data = data
        IF N_ELEMENTS(data.x) LT 2 THEN BEGIN
            s_e(nerror) = 'O_pressure not enough data'
            nerror = nerror+1
        ENDIF 
    ENDELSE 
    
    tplot_names, mag_press, names = names
    IF names(0) EQ '' THEN BEGIN 
        s_e(nerror)  = ' Mag,'  
        nerror = nerror+1
    ENDIF ELSE BEGIN 
        get_data, names(0), data = data
        IF N_ELEMENTS(data.x) LT 2 THEN BEGIN 
            s_e(nerror) = 'Mag not enough data'
            nerror = nerror+1
        ENDIF 
    ENDELSE 
;---------------------------------------------------------------------
    IF nerror EQ 0  THEN BEGIN
        h1_press = h_press & o1_press = o_press & mag_press = mag_press

        plasma_beta, h1_press, mag_press, O1_PRESSURE = o1_press   

        beta_name = 'TDMOM_EN00040_40000_SC'+sc_str+'_MTPRESSURE_SP0_ET0_All_O1_beta'

        get_data, beta_name, data = data, lim = lim, dlim = dlim
        store_data, beta_name+'_neg', data = data, lim = lim, dlim = dlim
        ylim, beta_name+'_neg', min(data.y), 0, 0
        options, beta_name+'_neg', 'psym', 7
        options, beta_name+'_neg', 'ytitle', 'neg beta'

        get_data, h_press, data = data, lim = lim, dlim = dlim
        store_data, h_press+'_neg', data = data, lim = lim, dlim = dlim
        ylim, h_press+'_neg', min(data.y) < (-1e-20), 0, 0
        options, h_press+'_neg', 'psym', 7
        options, h_press+'_neg', 'ytitle', 'neg h press'
        
        get_data, o_press, data = data, lim = lim, dlim = dlim 
        store_data, o_press+'_neg', data = data, lim = lim, dlim = dlim
        ylim, o_press+'_neg', min(data.y) < (-1e-20), 0, 0
        options, o_press+'_neg', 'psym', 7
        options, o_press+'_neg', 'ytitle', 'neg o press'
        
        get_data, mag_press, data = data, lim = lim, dlim = dlim
        store_data, mag_press+'_neg', data = data, lim = lim, dlim = dlim
        ylim, mag_press+'_neg', min(data.y) < (-1e-20), 0, 0
        options, mag_press+'_neg', 'psym', 7
        options, mag_press+'_neg', 'ytitle', 'neg mag press'     
  
        popen, 'plots/neg_beta/'+date_s+'.ps'
        tplot, [h_press, h_press+'_neg', o_press, o_press+'_neg', $
                mag_press, mag_press+'_neg', beta_name, beta_name+'_neg']   
        pclose
    ENDIF
    
ENDFOR 
spawn, 'mogrify -format png plots/neg_beta/*.ps'
END 

