PRO replot
sc = 4
sc_str = STRING(sc, FORMAT = '(i1.1)')

plot_mom = 1
mom_recalc = 1
idl_plot = 0
ps = 1
dumpdata = 0
store_new_data = 1
find_phase = 1

path = 'output_storm/'
;--------------------------------------------
; Read storm minimum Dst list from file: storm_min_dst.dat
; (this list is extracted from storm list spreadsheet)
; Store these times in variable min_dst (in tplot time format)
;---------------------------------------------------------------------
OPENR, unit, 'storm_min_dst.dat', /GET_LUN
min_dst = DBLARR(300)           ; assume no more than 300 events
dummy = ''
jj = 0
WHILE NOT EOF(unit) DO BEGIN   
    READF, unit, dummy
    min_dst(jj) = time_double(dummy)
    jj = jj + 1
ENDWHILE
min_dst = min_dst(0:jj-1)
CLOSE, unit
;---------------------------------------------------------------------
; Read CLUSTER perigee times list from file: sc4_perigee_times.dat
; (different list for each S/C)
; Store these times in variable petime (in tplot time format)
;---------------------------------------------------------------------
OPENR, unit, 'sc' + sc_str + '_perigee_times.dat', /GET_LUN
petime = DBLARR(3000)           ; assumes no more than 3000 perigee passes
dummy = ''
jj = 0l
WHILE NOT EOF(unit) DO BEGIN  
    READF, unit, dummy
    petime(jj) = time_double(dummy)
    jj = jj + 1
ENDWHILE
petime = petime(0:jj-1)
CLOSE, unit
;-------------------------------------------------------------------
;Read storm phase(prestorm, storm time and recovery time) form 
;file 'storm_phase.dat' and store them into arrays
;-------------------------------------------------------------------
OPENR, unit, 'storm_phase.dat', /GET_LUN
prestorm_start = DBLARR(300)
storm_onset = DBLARR(300) 
min_dst_new = DBLARR(300)
recovery_end = DBLARR(300) 
jj = 0l
dummy = ''
WHILE NOT EOF(unit) DO BEGIN

    READF, unit, dummy
    IF STRMID(dummy, 0, 5) NE 'Start' THEN BEGIN 
        prestorm_start(jj) = time_double(STRMID(dummy, 0, 20))
        storm_onset(jj) = time_double(STRMID(dummy, 20, 20))
        min_dst_new(jj) =  time_double(STRMID(dummy, 40, 24))
        recovery_end(jj) =  time_double(STRMID(dummy, 60, 24))

        jj = jj+1      
    ENDIF  
ENDWHILE
CLOSE, unit, /all
nstorm = jj

prestorm_start = prestorm_start(0:nstorm-1)
storm_onset = storm_onset(0:nstorm-1)
min_dst_new = min_dst_new(0:nstorm-1)
recovery_end = recovery_end (0:nstorm-1)
;stop
;---------------------------------------------------------------------
; For each minimum Dst time choose the time intarval of the plots.
; The time interval is selected is of three orbits duration
;---------------------------------------------------------------------

FOR ii = 17, N_ELEMENTS(min_dst)-1   DO BEGIN 

                                ;-----------------------------------------------------------------
                                ; find the time of the perigee pass closest to the minimum Dst
                                ;-----------------------------------------------------------------
    before_index = where ((petime - min_dst(ii)) < 0, counts)
    after_index = where ((petime - min_dst(ii)) > 0)          
    before_pe = petime(before_index(counts-1))
    after_pe = petime(after_index(0))
                                ;-----------------------------------------------------------------
                                ; time interval start time, end time and dt
                                ;-----------------------------------------------------------------
    time_start = petime(before_index(counts-2))
    time_end = petime(after_index(1))
    dt = time_end - time_start      
    displaytime =  6. * 3600. 
    idt = CEIL(dt / displaytime)
;stop
    FOR kk = 0,  idt-1  DO BEGIN   
        
        tplot_names, names = names
        store_data, DELETE = names

        PRINT, STRING(ii) + '   --' + STRING(kk)
                                ;-------------------------------------,63---------------------------
                                ; Set time interval in tplot
                                ;-----------------------------------------------------------------
        t_s = time_start + kk * displaytime
        t_e = t_s + displaytime
        t_dt = t_e - t_s
        timespan, t_s, t_dt, /SECONDS  

                                ;--------------------------Loading data----------------------------------
                                ; load other o_beam tplot
        ts = time_string(t_s)  
        te = time_string(t_e)
        date_s = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
        time_s = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
        date_e = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2)
        time_e = STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2) 

        flndata = path+'tplot_restore/o_beam_'+date_s+'_'+time_s
        print, FINDFILE(flndata+'.tplot', COUNT = ct)

        IF ct GT 0 THEN BEGIN 
            tplot_restore, filenames = flndata+'.tplot' 
            
                                ;find average time
            beam_name = 'PASPEC_SC'+sc_str+ $
                        '_IN0_PHICOMBINED_UNDIFFFLUX_SP3_ET0_All_AVG*_PAP_ET_beam'
            tplot_names, beam_name, names = names
            IF names(0) EQ '' THEN stop ELSE BEGIN 
                beam_name = names(0)
                get_data, names(0), data = data
                average_time = data.average_time
                start_time = data.start_time
                time = data.x
                ntime = N_ELEMENTS(time)
                at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL) 
                                ;-----------------------------------------------------------
                tb = 'PASPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP_ET_beam'
                eb = 'PASPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP_ET_beam' 
                cb = 'PASPEC_SC4_IN0_PHICOMBINED_UNDIFFFL X_SP3_ET0_All_AVG300_PAP_ET_beam'
                combine_et_pap, tb, eb, cb, $
                                start_time = start_time, average_time = average_time
                
                                ;--------------------------------------------------------------
                                ; plot mom with energy range and pa reange if keyword mom_recalc is set 
                
                IF KEYWORD_SET(mom_recalc) THEN BEGIN      

                    direction = ['t', 'e']
                    
                    FOR id = 0, 1  DO BEGIN 
                        IF direction(id) EQ 't' THEN phi_str = 'PHI90_270'
                        IF direction(id) EQ 'e' THEN phi_str = 'PHI270_90'
                        
                        en_range_name = 'ENSPEC_SC' + sc_str + '_IN0_'+ phi_str $
                                        +'_UNDIFFFLUX_SP3_ET0_All_AVG'+ at_str $
                                        +'_erange'
                        
                        get_data, en_range_name, data = data
                        energy_range_time = data.x
                        energy_range = data.y
                        energy_range_bins = data.energybins


                        pa_beam_name = 'PASPEC_SC'+sc_str+'_IN0_'+phi_str+ $
                                       '_UNDIFFFLUX_SP3_ET0_All_AVG'+at_str+'_PAP_PA_beam'

                        sat = sc
                        specie = 3
                        inst = 0 & eff_table = 0
                        units_name = 'DIFF FLUX'
                        
                        find_bins_from_pap, sat, specie, inst, units_name, eff_table, $
                                            average_time, pa_beam_name, bins_name, $
                                            start_time = start_time


                        get_data, bins_name, data = data
                        time_tp = data.x
                        angle_bins = data.y
                        
                        units_name = 'DIFF FLUX' ; 'Counts', 'NCOUNTS','RATE', 
                                ; 'NRATE', 'DIFF FLUX', 'EFLUX'
                        inst = 0 & eff_table = 0
                        angle = [[-90.0, 90.0], [0., 360.]] 
                        energy = [40., 40000.] ; !!!       
                        ntime = N_ELEMENTS(energy_range_time) 
                        
                        time_Tt = FLTARR(ntime) &  data_Tt = FLTARR(ntime)
                        time_D = FLTARR(ntime) & data_D = FLTARR(ntime)
                        time_Vt = FLTARR(ntime) &  data_Vt = FLTARR(ntime)
                        time_Pt = FLTARR(ntime) &  data_Pt = FLTARR(ntime)  
                        time_V = FLTARR(ntime) &  data_V = FLTARR(ntime, 3)
                        time_Tx = FLTARR(ntime) &  data_Tx = FLTARR(ntime)
                        time_Ty = FLTARR(ntime) &  data_Ty = FLTARR(ntime)
                        time_Tz = FLTARR(ntime) &  data_Tz = FLTARR(ntime)

                        FOR i_m = 0, ntime-1 DO BEGIN
                            moments = ['A']

                                ;                e_r = [FLOOR(energy_range(i_m, 0)),  $
                                ;                      CEIL(energy_range(i_m, 1))] 
                                ;               a_r =  REFORM(angle_bins(i_m, *))
                            
                            energy = [FLOOR(energy_range(i_m, 0)),  $
                                      CEIL(energy_range(i_m, 1))] 
                            use_bins = 1                      
                            
                            timespan, start_time + average_time*i_m, average_time, /SECONDS
                            bins = REFORM(angle_bins(i_m, *))                  
                            
                            index = where(bins GT 0, ct)
                            
                            IF ct GT 0 THEN BEGIN 
                                sat = [sc]
                                specie = [3]
                                
                                plot_3dmom_from_crib, sat, specie, inst, $
                                  moments, angle, energy, eff_table, recalc = 1, $
                                  bins = bins, use_bins = use_bins, $
                                  e_r = e_r, a_r = a_r
                                
                                e_min = STRCOMPRESS(STRING(energy(0), $
                                                           FORMAT = '(i5.5)'), /REMOVE_ALL)
                                e_max = STRCOMPRESS(STRING(energy(1), $
                                                           FORMAT = '(i5.5)'), /REMOVE_ALL)
                                
                                Tt_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                                          +'_MTTEMPERATURE_SP3_ET0_All_T'
                                Vt_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                                          + '_MTVELOCITY_SP3_ET0_All_T'
                                D_name = 'TDMOM_EN'+ e_min + '_' + e_max + '_SC' + sc_str  $
                                         +'_MTDENSITY_SP3_ET0_All'
                                Pt_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                                          + '_MTPRESSURE_SP3_ET0_All_T'
                                V_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                                         + '_MTVELOCITY_SP3_ET0_All'
                                Tx_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                                          +'_MTTEMPERATURE_SP3_ET0_All_X'
                                Ty_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                                          +'_MTTEMPERATURE_SP3_ET0_All_Y'
                                Tz_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                                          +'_MTTEMPERATURE_SP3_ET0_All_Z'
                                
                                ;---------check weather data has been loaded or not
                                nerror = 0
                                
                                tplot_names, Tt_name, names = names
                                IF names(0) EQ '' THEN nerror = nerror+1
                                tplot_names, Vt_name, names = names
                                IF names(0) EQ '' THEN nerror = nerror+1
                                tplot_names, D_name, names = names
                                IF names(0) EQ '' THEN nerror = nerror+1
                                tplot_names, Pt_name, names = names
                                IF names(0) EQ '' THEN nerror = nerror+1
                                
                                IF nerror EQ 0 THEN BEGIN 
                                    get_data, Tt_name, data = data, dlim = dlim_Tt, lim = lim_Tt
                                    time_Tt(i_m) = data.x(0)
                                    data_Tt(i_m) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                    
                                    get_data, D_name, data = data, dlim = dlim_D, lim = lim_D
                                    time_D(i_m) = data.x(0)
                                    data_D(i_m) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                    
                                    get_data, Vt_name, data = data, dlim = dlim_Vt, lim = lim_Vt
                                    time_Vt(i_m) = data.x(0)
                                    data_Vt(i_m) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                    
                                    get_data, Pt_name, data = data, dlim = dlim_Pt, lim = lim_Pt
                                    time_Pt(i_m) = data.x(0)
                                    data_Pt(i_m) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                    
                                    get_data, V_name, data = data, dlim = dlim_V, lim = lim_V
                                    time_V(i_m) = data.x(0)
                                    data_V(i_m, *) = TOTAL(data.y(*, *), 1)/N_ELEMENTS(data.y(*, 0))
                                    
                                    get_data, Tx_name, data = data, dlim = dlim_Tx, lim = lim_Tx
                                    time_Tx(i_m) = data.x(0)
                                    data_Tx(i_m) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                    
                                    get_data, Ty_name, data = data, dlim = dlim_Ty, lim = lim_Ty
                                    time_Ty(i_m) = data.x(0)
                                    data_Ty(i_m) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                    
                                    get_data, Tz_name, data = data, dlim = dlim_Tz, lim = lim_Tz
                                    time_Tz(i_m) = data.x(0)
                                    data_Tz(i_m) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                    
                                    tplot_names, 'TDMOM_EN'+e_min+'_'+e_max+'*', names = names
                                    store_data, delete = names
                                ENDIF  ELSE BEGIN 
                                    energy_range(i_m, *) = !VALUES.F_NAN
                                ENDELSE 
                            ENDIF ELSE BEGIN 
                                energy_range(i_m, *) = !VALUES.F_NAN
                            ENDELSE 
                            store_data, en_range_name, $
                                        data = {x:energy_range_time, y:energy_range, $
                                                energybins:energy_range_bins}
                        ENDFOR                   
                        
                        timespan, t_s, t_dt, /SECONDS  
                        
                        Tt_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                                  +phi_str+'_MTTEMPERATURE_SP3_ET0_All_T'+'_AVG'+at_str
                        Vt_name = 'TDMOM_ENVARIOUS'+'_SC' + sc_str+'_' $
                                  +phi_str + '_MTVELOCITY_SP3_ET0_All_T'+'_AVG'+at_str
                        D_name =  'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_'  $
                                  +phi_str +'_MTDENSITY_SP3_ET0_All'+'_AVG'+at_str
                        Pt_name = 'TDMOM_ENVARIOUS'+'_SC' + sc_str+'_' $
                                  +phi_str +'_MTPRESSURE_SP3_ET0_All_T'+'_AVG'+at_str
                        V_name = 'TDMOM_ENVARIOUS'+'_SC' + sc_str+'_' $
                                 +phi_str + '_MTVELOCITY_SP3_ET0_All'+'_AVG'+at_str
                        Tx_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                                  +phi_str+'_MTTEMPERATURE_SP3_ET0_All_X'+'_AVG'+at_str
                        Ty_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                                  +phi_str+'_MTTEMPERATURE_SP3_ET0_All_Y'+'_AVG'+at_str
                        Tz_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                                  +phi_str+'_MTTEMPERATURE_SP3_ET0_All_Z'+'_AVG'+at_str
                        
                        store_data,  Tt_name, data = {x:time_Tt,  y:data_Tt}, $
                                     dlim = dlim_Tt, lim = lim_Tt
                        store_data,  Vt_name, data =  {x:time_Vt,  y:data_Vt}, $
                                     dlim = dlim_Vt, lim = lim_Vt
                        store_data,  D_name, data =  {x:time_D,  y:data_D}, $
                                     dlim = dlim_D, lim = lim_D
                        store_data,  Pt_name, data =  {x:time_Pt,  y:data_Pt}, $
                                     dlim = dlim_Pt, lim = lim_Pt
                        store_data,  V_name, data =  {x:time_V,  y:data_V}, $
                                     dlim = dlim_V, lim = lim_V
                        store_data,  Tx_name, data = {x:time_Tx,  y:data_Tx}, $
                                     dlim = dlim_Tx, lim = lim_Tx
                        store_data,  Ty_name, data = {x:time_Ty,  y:data_Ty}, $
                                     dlim = dlim_Ty, lim = lim_Ty
                        store_data,  Tz_name, data = {x:time_Tz,  y:data_Tz}, $
                                     dlim = dlim_Tz, lim = lim_Tz
                        
                        options, Tt_name, 'ytitle', 'SC'+sc_str+' O!U+!N!C!CT!Davg!N (eV)'
                        options, Vt_name, 'ytitle', 'SC'+sc_str+' O!U+!N!C!CV!DT!N (kms!U-1!N)!C'
                        options, D_name, 'ytitle', 'SC'+sc_str+'!C!CO!U+!N!C!Cn (cm!U-3!N)' 
                        options, Pt_name, 'ytitle', 'SC'+sc_str+' O!U+!N!C!cP!Davg!N (nPa)'  
                        options, V_name, 'ytitle', 'SC'+sc_str+' O!U+!N!C!CV!D!N (kms!U-1!N)!C'

                        options, [D_name, Vt_name, Tt_name, Pt_name, V_name, Tx_name, $
                                  Ty_name, Tz_name], 'psym', 1
                                ; if negative T and P exist, plot and record in log 
                        index = (where (data_tt LT 0))
                        IF index(0) GT 0 THEN BEGIN 
                            store_data,  Tt_name+'_neg', data = {x:time_Tt,  y:data_Tt}, $
                                         dlim = dlim_Tt, lim = lim_Tt       
                            store_data,  Pt_name+'_neg', data =  {x:time_Pt,  y:data_Pt}, $
                                         dlim = dlim_Pt, lim = lim_Pt
                            
                            options, Tt_name+'_neg', 'ytitle',  $
                                     'SC' + sc_str +'!C!CO!U+!N T!Davg!N (eV)'  
                            options, Pt_name+'_neg', 'ytitle', $ 
                                     'SC'+sc_str+'!C!CO!U+!N P!Davg!N (nPa)' 
                            
                            options, [Tt_name+'_neg', Pt_name+'_neg'], 'psym', 7
                            
                            ylim, Tt_name+'_neg', min(data_tt(index))*10, 0, 0 
                            ylim, Pt_name+'_neg', min(data_pt (index))*10, 0, 0
                            
                            OPENU, unit, path+'log_neg_T.txt', /GET_LUN, /APPEND
                            PRINTF, unit,  ii, kk, '-----negative Tt and Pt-----'+direction(id)
                            FREE_LUN, unit         
                        ENDIF
;stop
                                ;----- calculate V_perp and V_par
                        v_perp, V_name

                    ENDFOR   
                ENDIF     

                IF KEYWORD_SET(find_phase)THEN BEGIN 
                    storm_phase = INTARR(ntime)
                    FOR  itime = 0, ntime-1 DO BEGIN
                        FOR istorm = 0, nstorm-1 DO BEGIN 
                            belong = 0
                            IF time(itime) GE prestorm_start(istorm) AND $
                              time(itime) LT storm_onset(istorm) THEN BEGIN 
                                storm_phase(itime) = 1
                                belong = belong+1
                                ;                stop
                            ENDIF 
                            IF time(itime) GE storm_onset(istorm) AND $
                              time(itime)  LT min_dst_new(istorm) THEN BEGIN 
                                storm_phase(itime) = 2
                                belong = belong+1
                            ENDIF 
                            IF time(itime) GE min_dst_new(istorm) AND $
                              time(itime) LT recovery_end(istorm) THEN BEGIN 
                                storm_phase(itime) = 3
                                belong = belong+1
                            ENDIF

                            IF belong GT 1 THEN stop
                        ENDFOR
                    ENDFOR 
                    store_data, 'storm_phase', data = {x:time, y:storm_phase}
                    ylim, 'storm_phase', -1, 4
                ENDIF  

                                ;--------------------------Plot-----------------------------------------
                
                                ; set options
                p01 = 'TDMOM_EN00040_40000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_P_total'
                p02 = 'TDMOM_EN00040_40000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_beta'
                
                                ;    p03 = 'ENSPEC_SC'+ sc_str +'_IN0_PHI90_270_UNDIFFFLUX_SP0_ET0_All'
                p04 = 'ENSPEC_SC'+ sc_str +'_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All'
                                ;   p05 = 'ENSPEC_SC'+ sc_str +'_IN0_PHI270_90_UNDIFFFLUX_SP0_ET0_All'
                p06 = 'ENSPEC_SC'+ sc_str +'_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All'
                
                                ;     p07 = 'Dst_Index'               
                p08 = 'MAG_SC4_B_xyz_gse_X'
                
                p09 = 'ENSPEC_SC'+sc_str+'_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG'+at_str
                p10 = 'ENSPEC_SC'+sc_str+'_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG'+at_str
                
                p11 = 'PASPEC_SC'+sc_str+'_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG'+at_str
                p12 = 'PASPEC_SC'+sc_str+'_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG'+at_str
                
                p13 = p11+'_PAP'
                p14 = P12+'_PAP'
                p15 = p13 +'_ET'
                p16 = p14 +'_ET'
                p17 = p15+'_beam'
                p18 = p16+'_beam'
                p19 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI90_270_MTDENSITY_SP3_ET0_All' $
                      + '_AVG'+at_str
                p20 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI270_90_MTDENSITY_SP3_ET0_All' $
                      +'_AVG'+at_str
                p21 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI90_270_MTVELOCITY_SP3_ET0_All_T' $
                      +'_AVG'+at_str
                p22 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI270_90_MTVELOCITY_SP3_ET0_All_T' $
                      +'_AVG'+at_str
                p23 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI90_270_MTTEMPERATURE_SP3_ET0_All_T' $
                      +'_AVG'+at_str
                p24 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI270_90_MTTEMPERATURE_SP3_ET0_All_T' $
                      +'_AVG'+at_str
                p25 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI90_270_MTTEMPERATURE_SP3_ET0_All_T' $
                      +'_AVG'+at_str+'_neg'
                p26 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI270_90_MTTEMPERATURE_SP3_ET0_All_T' $
                      +'_AVG'+at_str+'_neg'
                p27 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI90_270_MTPRESSURE_SP3_ET0_All_T' $
                      +'_AVG'+at_str
                p28 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI270_90_MTPRESSURE_SP3_ET0_All_T' $
                      +'_AVG'+at_str
                p29 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI90_270_MTPRESSURE_SP3_ET0_All_T' $
                      +'_AVG'+at_str+'_neg'
                p30 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI270_90_MTPRESSURE_SP3_ET0_All_T' $
                      +'_AVG'+at_str+'_neg'
                p31 = 'EPH_SC'+ sc_str+'_GSE_X'
                p32 = 'EPH_SC'+ sc_str+'_GSE_Y'
                p33 = 'EPH_SC'+ sc_str+'_GSE_Z'

                p34 = 'PASPEC_SC4_IN0_PHICOMBINED_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP_ET_beam'

                p35 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI90_270_MTVELOCITY_SP3_ET0_All' $
                      +'_AVG'+at_str+'_V_PAR_T'
                p36 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI270_90_MTVELOCITY_SP3_ET0_All' $
                      +'_AVG'+at_str+'_V_PAR_T'
                p37 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI90_270_MTVELOCITY_SP3_ET0_All' $
                      +'_AVG'+at_str+'_V_PERP_T'
                p38 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_PHI270_90_MTVELOCITY_SP3_ET0_All' $
                      +'_AVG'+at_str+'_V_PERP_T'
                p39 = 'EPH_SC'+ sc_str+'_GSM_X'
                p40 = 'EPH_SC'+ sc_str+'_GSM_Y'
                p41 = 'EPH_SC'+ sc_str+'_GSM_Z'
                p42 = 'MAG_SC'+sc_str+'_B_xyz_gse'            
                p43 = 'TDMOM_EN00040_40000_SC'+sc_str+'_MTDENSITY_SP0_ET0_All'
                p44 = 'TDMOM_EN00040_40000_SC'+sc_str+'_MTVELOCITY_SP0_ET0_All'
                p45 = 'TDMOM_EN00040_40000_SC'+sc_str+'_MTTEMPERATURE_SP0_ET0_All'
                
                p46 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +'PHI90_270'+'_MTTEMPERATURE_SP3_ET0_All_X'+'_AVG'+at_str
                p47 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +'PHI90_270'+'_MTTEMPERATURE_SP3_ET0_All_Y'+'_AVG'+at_str
                p48 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +'PHI90_270'+'_MTTEMPERATURE_SP3_ET0_All_Z'+'_AVG'+at_str
                
                p49 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +'PHI270_90'+'_MTTEMPERATURE_SP3_ET0_All_X'+'_AVG'+at_str
                p50 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +'PHI270_90'+'_MTTEMPERATURE_SP3_ET0_All_Y'+'_AVG'+at_str
                p51 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +'PHI270_90'+'_MTTEMPERATURE_SP3_ET0_All_Z'+'_AVG'+at_str
                p52 = 'storm_phase'                   
                
                options, '*', 'panel_size', 1
                options, [p06, p09, p10, p11, p12, p13, p14, p17, p18, p34], 'ztitle', ''
                
                ylim, p01, 0.01, 3, 1
                ylim, p02, 0.01, 10
                zlim, [p04, p06, p11, p12], 0.1, 100
                
                                ;        options, p03, 'ytitle', 'SC' + sc_str + ' H!U+!N!C!C(eV)' + '!C!C' + 'Tailward'
                options, p04, 'ytitle', 'SC' + sc_str + ' O!U+!N!C!C(eV)' + '!C!C' + 'Tailward'
                                ;       options, p05, 'ytitle', 'SC' + sc_str + ' H!U+!N!C!C(eV)' + '!C!C' + 'Earthward'
                options, p06, 'ytitle', 'SC' + sc_str + ' O!U+!N!C!C(eV)' + '!C!C' + 'Earthward'
                options, p08, 'ytitle', 'SC' + sc_str + '!C!CBx (nT)'
                options, p09, 'ytitle', 'SC' + sc_str + ' O!U+!N (eV)!C!CTailward!C!CAVG-'+at_str
                options, p10, 'ytitle', 'SC' + sc_str + ' O!U+!N (eV)!C!CEarthward!C!CAVG-'+at_str
                
                options, p34, 'ytitle', 'SC'+sc_str+' O!U+!C!CBEAM!C!CE-----T'
                options, [p09+'_erange', p10+'_erange'], 'color', 2                         
                options,  p11, 'ytitle', 'Pitch Angle!C!CVarious EN'
                options,  p12, 'ytitle', 'Pitch Angle!C!CVarious EN'                  
                
                var_label = 'EPH_SC' + sc_str + '_'
                var_label = var_label + ['MLT', 'GSE_X', 'GSE_Y', 'GSE_Z', 'DIST']
                
                IF plot_mom EQ 1 THEN BEGIN 
                    ylim, p19, 0.001, 1, 1
                    ylim, p20, 0.001, 1, 1
                    ylim, p21, 0, 300, 0
                    ylim, p22, 0, 300, 0
                    ylim, p23, 0.1, 10000, 1
                    ylim, p24, 0.1, 10000, 1
                    ylim, p27, 1e-12, 1e-2, 1
                    ylim, p28, 1e-12, 1e-2, 1
                    ylim, [p35, p36], -100, 100
                    ylim, [p37, p38], 0, 100
                    options, [p20, p22, p24, p26, p28, p30, p36, p38, P49, P50, P51], 'color', 2
                    options, [p35, p36, p37, p38], 'labels', ''

                    options, [p19, p20], 'ytitle', 'SC'+sc_str+' O!U+!N!C!Cn (cm!U-3!N)'
                    options, [p21, p22], 'ytitle', 'SC'+sc_str+' O!U+!N!C!CV!LT!N (kms!U-1!N)!C' 
                    options, [p35, p36], 'ytitle', 'SC'+sc_str+' O!N+!N!C!CV!L//!N (kms!U-1!N)!C'
                    options, [p37, p38], 'ytitle', 'SC'+sc_str+' O!N+!N!C!CV!DL!N (kms!U-1!N)!C'
                    options, [p23, p24, p25, p26], 'ytitle', 'SC' + sc_str +' O!U+!N!C!CT!Davg!N (eV)'
                    options, [p27, p28, p29, p30], 'ytitle', 'SC'+sc_str+' O!U+!N!C!cP!Davg!N (nPa)'   
                    options, [p46, p49], 'ytitle', 'SC' + sc_str +' O!U+!N!C!CT!D//!N (eV)'
                    options, [p47, p48], 'ytitle', 'SC' + sc_str +' O!U+!N!C!CT!DL1!N (eV)'
                    options, [p50, p51], 'ytitle', 'SC' + sc_str +' O!U+!N!C!CT!DL2!N (eV)'
                ENDIF 

                                ;plot in idl windows
                IF KEYWORD_SET(idl_plot) THEN BEGIN 
                    tplot, [p02, p04, p06, p34, p19, p21, p35, p37, p23, $
                            p25, p27, p29, p46, p47, p48, p52], $
                           var_label = var_label
                    tplot_panel, v = p04, o = p09+'_epcut_beam', psym = 0
                    tplot_panel, v = p06, o = p10+'_epcut_beam', psym = 0
                    tplot_panel, v = p04, o = p09+'_erange'
                    tplot_panel, v = p06, o = p10+'_erange'
                    tplot_panel, v = p19, o = p20, psym = 1
                    tplot_panel, v = p21, o = p22, psym = 1
                    tplot_panel, v = p23, o = p24, psym = 1
                    tplot_panel, v = p27, o = p28, psym = 1
                    tplot_panel, v = p35, o = p36, psym = 1
                    tplot_panel, v = p37, o = p38, psym = 1
                    tplot_panel, v = p25, o = p26, psym = 1
                    tplot_panel, v = p29, o = p30, psym = 1
                    tplot_panel, v = p46, o = p49, psym = 1
                    tplot_panel, v = p47, o = p50, psym = 1
                    tplot_panel, v = p48, o = p51, psym = 1
                    
                    yline, p02, offset = 0.05, col = 1
                    yline, p02, offset = 1, col = 1
                ENDIF 
;stop
                                ;-----------------------------------------------------------------------
                                ; Plot the graph in PS file if ps is set to 1 
                                ;--------------------------------------------------------------------------
                
                IF ps EQ 1 THEN BEGIN  
                    IF plot_mom EQ 1 THEN BEGIN 
                        fln = path+'plots/mom_with_neg/storm_o_beam_mom' + $
                              date_s + '_' + time_s + '.ps' 
                        popen, fln, /port
                        
                        tplot, [p09, p10, p34, p19, p21, p35, p37, p23, p25, p27, p29, p46, p47, p48], $
                               var_label = var_label
                        tplot_panel, v = p09, o = p09+'_epcut_beam', psym = 0
                        tplot_panel, v = p10, o = p10+'_epcut_beam', psym = 0
                        tplot_panel, v = p09, o = p09+'_erange', psym = 0
                        tplot_panel, v = p10, o = p10+'_erange', psym = 0
                        tplot_panel, v = p19, o = p20, psym = 1
                        tplot_panel, v = p21, o = p22, psym = 1
                        tplot_panel, v = p23, o = p24, psym = 1 
                        tplot_panel, v = p25, o = p26, psym = 1 
                        tplot_panel, v = p27, o = p28, psym = 1
                        tplot_panel, v = p29, o = p30, psym = 1 
                        tplot_panel, v = p35, o = p36, psym = 1
                        tplot_panel, v = p37, o = p38, psym = 1
                        tplot_panel, v = p46, o = p49, psym = 1
                        tplot_panel, v = p47, o = p50, psym = 1
                        tplot_panel, v = p48, o = p51, psym = 1

                        yline, p02, offset = 0.05, col = 1
                        yline, p02, offset = 1, col = 1                    
                        
                        pclose
                        
                        fln = path+'plots/mom_2pages/storm_o_beam_mom' +$
                              date_s + '_' + time_s + '_page1.ps' 
                        popen, fln, /port
                        
                        tplot, [p52, p02, p04, p06, p34, p19, p21, p35, p37, p23, p27, P46, P47, P48], $
                               var_label = var_label
                        tplot_panel, v = p04, o = p09+'_epcut_beam', psym = 0
                        tplot_panel, v = p06, o = p10+'_epcut_beam', psym = 0
                        tplot_panel, v = p04, o = p09+'_erange'
                        tplot_panel, v = p06, o = p10+'_erange'
                        tplot_panel, v = p19, o = p20, psym = 1
                        tplot_panel, v = p21, o = p22, psym = 1
                        tplot_panel, v = p23, o = p24, psym = 1 
                        tplot_panel, v = p27, o = p28, psym = 1
                        tplot_panel, v = p35, o = p36, psym = 1
                        tplot_panel, v = p37, o = p38, psym = 1
                        tplot_panel, v = p46, o = p49, psym = 1
                        tplot_panel, v = p47, o = p50, psym = 1
                        tplot_panel, v = p48, o = p51, psym = 1

                        yline, p02, offset = 0.05, col = 1
                        yline, p02, offset = 1, col = 1       
                        
                        pclose
                        
                        fln = path+'plots/mom_2pages/storm_o_beam_mom' + $
                              date_s + '_' + time_s + '_page2.ps' 
                        popen, fln, /port
                        
                        tplot, [p34, p09, p11, p13, p17, p10, p12, p14, p18, p08], $
                               var_label = var_label
                        tplot_panel, v = p09, o = p09+'_epcut', psym = -7
                        tplot_panel, v = p10, o = p10+'_epcut', psym = -7
                        yline, p08, col = 3
                        
                        pclose
                    ENDIF 
                ENDIF                       
;stop
                                ;------------------------------------------------------------------------
                                ;dump the data out if dumpdata is set to 1  
                                ;---------------------------------------------------------------------
                IF dumpdata EQ 1 THEN BEGIN 
                    title_set =  ['         flag  ', $
                                  '         Beta  ', $
                                  '      GSE_X(Re)', $
                                  '      GSE_Y(Re)', $
                                  '      GSE_Z(Re)', $
                                  '       en_tail ', $
                                  '       pa_tail ', $
                                  '      flux_tail', $
                                  '   Density_tail', $
                                  '   V_total_tail', $
                                  '     V_par_tail', $
                                  '    V_perp_tail', $
                                  '   T_total_tail', $
                                  '   P_total_tail', $
                                  '       en_earth', $
                                  '       pa_earth', $
                                  '     flux_earth', $
                                  '  Density_earth', $
                                  '  V_total_earth', $
                                  '    V_par_earth', $
                                  '   V_perp_earth', $
                                  '  T_total_earth', $
                                  '  P_total_earth', $
                                  '      GSM_X(Re)', $
                                  '      GSM_Y(Re)', $
                                  '      GSM_Z(Re)', $
                                  '     MAG_X(GSE)', $
                                  '     MAG_Y(GSE)', $
                                  '     MAG_Z(GSE)', $
                                  '      H_DENSITY', $
                                  '       H_V_X   ', $
                                  '       H_V_Y   ', $
                                  '       H_V_Z   ', $
                                  '       H_T_X   ', $
                                  '       H_T_Y   ', $
                                  '       H_T_Z   ', $
                                  '      T_x_tail ', $
                                  '      T_y_tail ', $
                                  '      T_z_tail ', $
                                  '      T_x_earth', $
                                  '      T_y_earth', $
                                  '      T_z_earth', $
                                  '    Storm_Phase']

                    nterm = N_ELEMENTS(title_set)
                    time_dd = time
                    index_valid = where(time_dd GE t_s AND time_dd  LE t_e, ct)
                    IF index_valid(0) GE 0 THEN BEGIN 
                        time_dd = time_dd(index_valid)
                        
                        n_time = N_ELEMENTS(time_dd) 
                        title_dd = STRARR(n_time, nterm )
                        data_dd = DBLARR(n_time, nterm)
                        FOR i_time = 0, n_time-1 DO BEGIN 
                            title_dd(i_time, *) = title_set
                        ENDFOR  
                        
                                ;Beta
                        get_data, p02, data = data
                        data_y = INTERPOL(data.y, data.x, time_dd)
                        index = where (~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 1) = data_y

                                ;GSE X
                        get_data, p31, data = data  
                        data_y = INTERPOL(data.y, data.x, time_dd)
                        index = where (~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF 
                        data_dd(*, 2) = data_y
                        
                                ;GSE Y
                        get_data, p32, data = data  
                        data_y = INTERPOL(data.y, data.x, time_dd)
                        index = where (~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF 
                        data_dd(*, 3) = data_y
                        
                                ;GSE Z
                        get_data, p33, data = data    
                        data_y = INTERPOL(data.y, data.x, time_dd)
                        index = where (~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF 
                        data_dd(*, 4) = data_y

                                ;GSM X
                        get_data, p39, data = data    
                        data_y = INTERPOL(data.y, data.x, time_dd)
                        index = where (~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF 
                        data_dd(*, 23 ) = data_y/6370.
                        
                                ;GSM Y
                        get_data, p40, data = data    
                        data_y = INTERPOL(data.y, data.x, time_dd)
                        data_x = time_dd
                        index = where (~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF 
                        data_dd(*, 24 ) = data_y/6370.
                        
                                ;GSM Z
                        get_data, p41, data = data    
                        data_y = INTERPOL(data.y, data.x, time_dd)
                        index = where (~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF 
                        data_dd(*, 25 ) = data_y/6370.
                        
                                ;tail energy
                        get_data, p09+'_epcut_beam', data = data
                        data_y = data.y(index_valid)
                        index = where (~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF 
                        data_dd(*, 0) = data_dd(*, 0) + (data_y GT 0)
                        data_dd(*, 5) = data_y                           
                        
                                ; tail pap
                        get_data, p17, data = data
                        data_y = data.y(index_valid, *)
                        data_v = data.v(index_valid, *)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                            data_v(index) = 0
                        ENDIF 
                        data_dd(*, 6) = total(data_v(*, *), 2)
                        data_dd(*, 7) = total(data_y(*, *), 2)
                        IF plot_mom EQ 1 THEN BEGIN           
                                ;tail Density
                            get_data, p19, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 8) = data_y                
                            
                                ;tail V_total
                            get_data, p21, data = data
                            data_y =  data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 9) = data_y
                            
                                ;tail V_par
                            get_data, p35, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 10) = data_y
                            
                                ;tail V_perp
                            get_data, p37, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 11) = data_y
                            
                                ;tail T_total
                            get_data, p23, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 12) = data_y
                            
                                ;tail P_total
                            get_data, p27, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 13) = data_y
                        ENDIF      
                                ;   earth energy
                        get_data, p10+'_epcut_beam', data = data
                        data_y = data.y(index_valid)
                        index = where ( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF 
                        data_dd(*, 0) = data_dd(*, 0) - (data_y GT 0)*10
                        data_dd(*, 14) = data_y
                        
                                ; dealing with flag 
                        index = where(data_dd(*, 0) EQ -10)
                        IF index(0) GE 0 THEN data_dd(index, 0) = -1
                        index = where(data_dd(*, 0) EQ -9)
                        IF index(0) GE 0 THEN data_dd(index, 0) = 2 
                        
                                ; earth pap
                        get_data, p18, data = data
                        data_y = data.y(index_valid, *)
                        data_v = data.v(index_valid, *)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                            data_v(index) = 0
                        ENDIF 
                        data_dd(*, 15) = total(data_v(*, *), 2)
                        data_dd(*, 16) = total(data_y(*, *), 2)  
                        IF plot_mom EQ 1 THEN BEGIN       
                                ;earth Density
                            get_data, p20, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 17) = data_y                
                            
                                ;earth V_total
                            get_data, p22, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data.y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data.y(index) = 0
                            ENDIF       
                            data_dd(*, 18) = data_y
                            
                                ;earth V_par
                            get_data, p36, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 19) = data_y
                            
                                ;earth V_perp
                            get_data, p38, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 20) = data_y
                            
                                ;earth T_total
                            get_data, p24, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 21) = data_y
                            
                                ;earth P_total
                            get_data, p28, data = data
                            data_y = data.y(index_valid)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 22) = data_y
                        ENDIF       
                                ; mag X_gse
                        get_data, p42, data = data
                        data_y = INTERPOL(data.y(*, 0), data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 26 ) = data_y
                        
                                ; mag Y_gse
                        data_y = INTERPOL(data.y(*, 1), data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 27 ) = data_y
                        
                                ; mag Z_gse
                        data_y = INTERPOL(data.y(*, 2), data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 28 ) = data_y
                        
                                ; H+ Density
                        get_data, p43, data = data
                        data_y = INTERPOL(data.y, data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 29 ) = data_y
                        
                                ;  H+ Vx
                        get_data, p44, data = data
                        data_y = INTERPOL(data.y(*, 0), data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 30 ) = data_y
                        
                                ;  H+ Vy
                        data_y = INTERPOL(data.y(*, 1), data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 31 ) = data_y

                                ;  H+ Vz
                        data_y = INTERPOL(data.y(*, 2), data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 32 ) = data_y

                                ;  H+ Tx
                        get_data, p45, data = data
                        data_y = INTERPOL(data.y(*, 0), data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 33 ) = data_y
                        
                                ;  H+ Ty
                        data_y = INTERPOL(data.y(*, 1), data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 34 ) = data_y
                        
                                ;  H+ Tz
                        data_y = INTERPOL(data.y(*, 2), data.x, time_dd)
                        index = where( ~FINITE(data_y), ct)
                        IF ct GT 0 THEN BEGIN 
                            data_y(index) = 0
                        ENDIF       
                        data_dd(*, 35 ) = data_y

                        IF plot_mom EQ 1 THEN BEGIN  
                                ;  T_x_tail
                            get_data, p46, data = data
                            data_y = INTERPOL(data.y(*, 0), data.x, time_dd)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 36 ) = data_y

                                ;  T_y_tail
                            data_y = INTERPOL(data.y(*, 1), data.x, time_dd)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 37 ) = data_y

                                ;  T_z_tail
                            data_y = INTERPOL(data.y(*, 2), data.x, time_dd)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 38 ) = data_y

                                ;  T_x_earth
                            get_data, p47, data = data
                            data_y = INTERPOL(data.y(*, 0), data.x, time_dd)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 39 ) = data_y

                                ;  T_y_earth
                            data_y = INTERPOL(data.y(*, 1), data.x, time_dd)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 40 ) = data_y

                                ;  T_z_earth
                            data_y = INTERPOL(data.y(*, 2), data.x, time_dd)
                            index = where( ~FINITE(data_y), ct)
                            IF ct GT 0 THEN BEGIN 
                                data_y(index) = 0
                            ENDIF       
                            data_dd(*, 41 ) = data_y
                        ENDIF 
                                ; Storm_phase
                        data_dd(*, 42) = storm_phase(index_valid)                    
                        
                        IF date_s EQ date_e THEN BEGIN 
                            str = {x:time_dd, y:data_dd, v:title_dd}
                            store_data, 'dump_data', data = str   
                            fln_dump = path+'data/'+ '/storm_o_beam_'+date_s+'.dat'
                            dump_data, 'dump_data', file_out = fln_dump
                        ENDIF ELSE BEGIN 
                            
                            midnight = time_double(STRMID(te, 0, 10)+'/00:00:00')
                            fday = where(time_dd LT midnight)
                            sday = where(time_dd GE midnight)
                            
                            IF fday(0) GE 0 THEN BEGIN 
                                str = {x:time_dd(fday), y:data_dd(fday, *), v:title_dd(fday, *)}
                                store_data, 'dump_data_f', data = str
                                fln_dump = path+'data/'+ '/storm_o_beam_'+date_s+'.dat'                                                                           
                                dump_data,  'dump_data_f', file_out = fln_dump
                            ENDIF 
                            
                            IF sday(0) GE 0 THEN BEGIN 
                                str = {x:time_dd(sday), y:data_dd(sday, *), v:title_dd(sday, *)}
                                store_data, 'dump_data_s', data = str
                                fln_dump = path+'data/'+ '/storm_o_beam_'+date_e+'.dat'                                                                           
                                dump_data, 'dump_data_s', file_out = fln_dump
                            ENDIF 
                        ENDELSE   
                        
                        tplot_names, 'dump_data*', names = names
                        store_data, delete = names
                        
                    ENDIF       
                    
                ENDIF            
                
                                ; restore the data into the original name if needed
                IF KEYWORD_SET(store_new_data) THEN BEGIN 
                    
                    tplot_names, 'TDMOM*X', names = names
                    store_data, delete = names
                    tplot_names, 'TDMOM*Y', names = names
                    store_data, delete = names
                    tplot_names, 'TDMOM*Z', names = names
                    store_data, delete = names
                    tplot_names, 'TDMOM*SP0*T', names = names
                    store_data, delete = names
                    tplot_names, 'TDMOM*H', names = names
                    store_data, delete = names
                    tplot_names, '*VcrossB*', names = names
                    store_data, delete = names
                    tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+ '*MTVELOCITY_SP3_ET0_All_AVG' $
                                 +at_str+'_V_PERP1', names = names
                    store_data, delete = names
                    tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP3_ET0_All_AVG' $
                                 +at_str+'_V_PERP2', names = names
                    store_data, delete = names
                    tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+ '*MTVELOCITY_SP3_ET0_All_AVG' $
                                 +at_str, names = names
                    store_data, delete = names
                    tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP3_ET0_All_AVG' $
                                 +at_str+'_T', names = names
                    store_data, delete = names 
                    
                    tplot_save, filename = flndata
                ENDIF  
                
                OPENU, unit, path+'log_replot.txt', /GET_LUN, /APPEND
                PRINTF, unit, ii, kk, 'file loaded'
                FREE_LUN, unit       
            ENDELSE 
        ENDIF  ELSE  BEGIN 
            PRINT, 'no restore file found'
            OPENU, unit, path+'log_replot.txt', /GET_LUN, /APPEND
            PRINTF, unit, ii, kk, 'no file found'
            FREE_LUN, unit       
        ENDELSE     
    ENDFOR      
ENDFOR     
stop
END
