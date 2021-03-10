;Purpose: Identify O+ beam using from energy spec, pitch angle spec
;         and then make corresponding mom plot in page1 the whole procedure
;         plot in page2
;
;Input: sc           : Cluster no. if not set the default is 4
;       average_time : in seconds , if not set the default is 5 min
;       
;Keywords: idl_plot  : plot the result plot in idl_window
;          ps        : plot the result plot in dumpdata,
;          dumpdata  : output data
;          globe_plot: plot a set of globe plot to show the selected
;                      range for plotting mom
;          store_data: store_data into .tplot  default: 1 
;
;Output: Depends on Keywords settings 
;        There will also be two .log files
;
;Written by Jing Liao  02/01/2008
;
PRO find_o_beam, sc = sc, $
                 sp=sp,$
                 average_time = average_time, $
                 idl_plot = idl_plot, $
                 ps = ps, $
                 dumpdata = dumpdata, $
                 globe_plot = globe_plot, $ 
                 store_data = store_data, $
                 beam_recalc = beam_recalc, $
                 mom_recalc = mom_recalc,  $
                 path = path, $
                 find_phase = find_phase, $
                 plot_lowcount_filter =  plot_lowcount_filter, $
                 displaytime = displaytime, $
                 plot_mom = plot_mom, $
                 add_imf = add_imf, $
                 plot_sw_sheath = plot_sw_sheath, $
                 use_angle_range = use_angle_range, $
                 use_energy_range = use_energy_range, $
                 plot_imf = plot_imf, $
                 beam_angle_range =  beam_angle_range, $
                 dfit_temperature = dfit_temperature, $
                 show_fit = show_fit, $
                 inst_input = inst_input, $
                 add_eflux = add_eflux, $
                 add_anodes = add_anodes,$
                 use_hiabeta=use_hiabeta,$
                 only_in_lobe=only_in_lobe,$
                 plot_add_eflux_procedure=plot_add_eflux_procedure,$
                 add_distfunc = add_distfunc,$
                 plot_add_distfunc_procedure=plot_add_distfunc_procedure,$
                 flux_threshold=flux_threshold,$
                 add_Vgse=add_Vgse

;-----------------------------------------------------
;check keywords  
;---------------------------------------------------
IF NOT keyword_set(sc) THEN sc = 4 
sc_str = STRING(sc, FORMAT = '(i1.1)')
if not keyword_set(sp) then sp=3
;sp=0
sp_str= STRING(sp, FORMAT = '(i1.1)')

IF NOT KEYWORD_SET(AVERAGE_TIME) THEN average_time = 5 * 60 ;in seconds
at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL) 
average_time = FLOAT(average_time)

IF NOT KEYWORD_SET(inst_input) THEN inst_input = 0
inst_str = STRING(inst_input, FORMAT = '(i1.1)')
; some basic settings
n_delay=1
low_counts_line = 9
pa_counts_line = low_counts_line/88.
IF  inst_input EQ  0 THEN  BEGIN 
    phi_str_set = ['PHI90_270', 'PHI270_90'] 
    angle_set = [[[-90, 90], [90, 270]], [[-90, 90], [270, 90]]] 
ENDIF  ELSE  BEGIN 
    IF inst_input EQ 1 THEN BEGIN 
        phi_str_set = ['PHI270_90', 'PHI90_270'] 
        angle_set = [[[-90, 90], [270, 90]], [[-90, 90], [90, 270]]]
    ENDIF ELSE stop
ENDELSE 

if keyword_set(use_hiabeta) THEN beta_name = 'HIA_L2_MOM_SC'+sc_str+'_MTPRESSURE_SP0_PR202_beta' $
else  beta_name =  'TDMOM_EN00040_40000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_beta' 

If Not keyword_set(use_angle_range) then use_angle_range=1
If Not keyword_set(beam_angle_range) then beam_angle_range=22.5
If Not keyword_set(use_energy_range) then use_energy_range=1
;-------------------------------------------------------------
;Delete all the string stored data in
;order to make sure the program can run correctly
;-----------------------------------------------------------
tplot_names, names = names
store_data, DELETE = names
;----------------------------------------------------------------
;Get the time interval from timespan
;-----------------------------------------------------------------
get_timespan, interval

t_s = interval(0)  
t_e = interval(1)
t_dt = t_e-t_s
ts = time_string(t_s)  
te = time_string(t_e)
date_s = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
time_s = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
date_e = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2)
time_e = STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2)
flndata = path+'tplot_restore/o_beam_'+date_s+'_'+time_s

;----------------------------------------------------------------
;Loading data
;----------------------------------------------------------------
IF NOT KEYWORD_SET(beam_recalc) THEN BEGIN
    print, FINDFILE(flndata+'.tplot.gz', COUNT = ct_beam)
    if ct_beam eq 0 then print, FINDFILE(flndata+'.tplot', COUNT = ct_beam)
    IF ct_beam GT 0 THEN BEGIN 
        spawn,'gzip -d '+flndata+'.tplot.gz'
        tplot_restore, filenames = flndata+'.tplot' 
    ;    spawn,'gzip -9 '+flndata+'.tplot'
;find average time
        beam_name = 'PASPEC_SC'+sc_str+ $
          '_IN'+inst_str+'_PHICOMBINED_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG*_PAP_ET_beam'
        tplot_names, beam_name, names = names
        IF names(0) NE  '' THEN BEGIN 
            beam_name = names(0)
            get_data, beam_name, data = data
            average_time = data.average_time
            start_time = data.start_time
            END_time = data.end_time
            ntime = floor((END_time-start_time)/average_time)
            time_avg = data.x
            n_avg = N_ELEMENTS(time)
            at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL)

        ENDIF ELSE  ct_beam = -1
    ENDIF ELSE ct_beam = -1
ENDIF ELSE ct_beam = 0 

IF ct_beam EQ 0 THEN BEGIN      
;-- Load CLUSTER ephemeris--
    sat = [sc] 
    get_cluster_ephemeris, sat, /GSE_X, /GSE_Y, /GSE_Z, /DIST, /MLT $
      , /GSM_X, /GSM_Y, /GSM_Z, /ILAT_D
    
;-- Load energy spectra - tailward (or earthward for HIA)--
    sat = [sc]
    specie = [sp]
    angle = [[-90, 90], [90, 270]]
    inst = inst_input & units_name = 'DIFF FLUX' 
    eff_table = 0
    
    plot_en_spec_from_crib, sat, specie, inst,  $
      units_name, angle, eff_table, recalc = 1
    
    sat = [sc]
    specie = [sp]
    angle = [[-90, 90], [90, 270]]
    inst = inst_input & units_name = 'Counts' & eff_table = 0

    plot_en_spec_from_crib, sat, specie, inst, $
      units_name, angle, eff_table, recalc = 1
    
;-- Load energy spectra - earthward (or earthward for HIA)--
    sat = [sc]
    specie = [sp]
    angle = [[-90, 90], [270, 90]]
    inst = inst_input & units_name = 'DIFF FLUX' & eff_table = 0

    plot_en_spec_from_crib, sat, specie, inst, $
      units_name, angle, eff_table, recalc = 1
    
    sat = [sc]
    specie = [sp]
    angle = [[-90, 90], [270, 90]]
    inst = inst_input & units_name = 'Counts' & eff_table = 0

    plot_en_spec_from_crib, sat, specie, inst, $
      units_name, angle, eff_table, recalc = 1

;-- Load CLUSTER Magnetic field--
    sat = sc
    plot_mag_from_crib, sat, POLAR = 1, GSM = 1

;-- Load CLUSTER H+ and O+ moments--
    IF sc EQ 4 THEN BEGIN 
        sat = [sc, sc]
        specie = [0, 3]
        moments = ['A', 'A']
        angle = [[-90, 90], [0, 360]]
        energy = [40., 40000.]
        inst = 0 &  eff_table = 0
        
        plot_3dmom_from_crib, sat, specie, inst, $
          moments, angle, energy, eff_table, recalc = 0
    ENDIF ELSE BEGIN 
        IF sc EQ 1 OR sc EQ 3 THEN BEGIN 
            sat = [sc]
            specie = [0]
            moments = ['P']
            inst = 1

            plot_l2mom_from_crib, sat, specie, inst, moments, units = units

            sat = [sc]
            specie = [0]
            moments = ['A']
            angle = [[-90, 90], [0, 360]]
            energy = [40., 40000.]
            inst = 0 &  eff_table = 0
            
            plot_3dmom_from_crib, sat, specie, inst, $
              moments, angle, energy, eff_table, recalc = 0
        ENDIF ELSE stop 
    ENDELSE 
;------------------------------------------------------------------
;  Checking the errors
;----------------------------------------------------------------  
;---------Check weather data required to calculate beta are loaded 
;Including : H pressure, O pressure and mag pressure 
;nerror = 0 = > no problem to calculate beta      ------------
    h_press =  'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP0_ET0_All'
    o_press = 'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP3_ET0_All'
    mag_press =  'MAG_SC' + sc_str +'_B_xyz_gse_MAG_PR'
    hia_press =  'HIA_L2_MOM_SC'+sc_str+'_MTPRESSURE_SP0_PR202'

    s_e = STRARR(10)
    nerror = 0
    IF sc EQ 4 THEN BEGIN 
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
        
        tplot_names, o_press, names = names
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
    ENDIF ELSE BEGIN  
        IF sc EQ 1 OR sc EQ 3 THEN BEGIN 
            tplot_names, hia_press, names = names
            IF names(0) EQ '' THEN BEGIN 
                s_e(nerror) = ' HIA_pressure,'
                nerror = nerror+1
            ENDIF ELSE BEGIN 
                get_data, names(0), data = data
                IF N_ELEMENTS(data.x) LT 2 THEN BEGIN 
                    s_e(nerror) = 'HIA_pressure not enough data'
                    nerror = nerror+1
                ENDIF 
            ENDELSE 
        ENDIF ELSE stop
    ENDELSE 

    tplot_names, mag_press, names = names
    IF names(0) EQ '' THEN BEGIN 
        s_e(nerror)  = ' Mag,'  
        nerror = nerror+1
    ENDIF ELSE BEGIN 
        get_data, names(0), data = data
        IF N_ELEMENTS(data.x) LT 2 THEN BEGIN 
            s_e(nerror) = 'Mag not enought data'
            nerror = nerror+1
        ENDIF 
    ENDELSE 

;------ check the energy spec data during display interval for all 
;keep only data during displaytime 
;and resave them in original names  -------   
    flux_tail_name = 'ENSPEC_SC' + sc_str + '_IN'+ $
      inst_str+'_'+PHI_STR_SET(0)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All'
    counts_tail_name = 'ENSPEC_SC' + sc_str + '_IN'+ $
      inst_str+'_'+phi_str_set(0)+'_UNCOUNTS_SP'+sp_str+'_ET0_All'
    flux_earth_name = 'ENSPEC_SC'  + sc_str + '_IN'+ $
      inst_str+'_'+phi_str_set(1)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All'
    counts_earth_name = 'ENSPEC_SC' + sc_str + '_IN'+ $
      inst_str+'_'+phi_str_set(1)+'_UNCOUNTS_SP'+sp_str+'_ET0_All'
    
    energy_name = [flux_tail_name, counts_tail_name, $
                   flux_earth_name, counts_earth_name]
    
    avg_error = FLTARR(4)
    avg_error_str = ['no energy data loaded', 'data less than 2', $
                     'data too short to be averaged', 'empty data']
; if there is one error in one energy data, there will be the same
; error in the other three, but 
    FOR i_en = 0, N_ELEMENTS(energy_name)-1 DO BEGIN        
; wether energy data loaded or not  
        tplot_names, energy_name(i_en), names = names   
        IF names(0) EQ '' THEN BEGIN 
            s_e(nerror) = 'no energy data loaded'+'('+string(i_en)+')'
            nerror = nerror+1
        ENDIF ELSE BEGIN 
; data during display time need to be more than 2 
            get_data, energy_name(i_en),  data = data, dlim = dlim, lim = lim
            index = where(data.x GT (t_s+1) AND data.x LE t_e, ct)
            IF ct LT 2 THEN BEGIN 
                s_e(nerror) = 'energy data less than 2'
                nerror = nerror+1
            ENDIF ELSE BEGIN   
; original data interval need to be shorter than average_time/3
                n_time = N_ELEMENTS(index)
                stime = data.x(index(n_time-1))
                etime = data.x(index(0))
                IF FLOOR((stime-etime)/average_time) LE 3 THEN BEGIN  
                    s_e(nerror) = 'data too short to be averaged'
                    nerror = nerror+1
                ENDIF ELSE BEGIN 
; data cannot be averaged if data.y are all 0
                    IF  TOTAL(data.y(index, *)) EQ 0 THEN BEGIN 
                        s_e(nerror) = 'empty energy data'
                        nerror = nerror+1 
                    ENDIF ELSE BEGIN 
                        str = {x:data.x(index), y:data.y(index, *), v:data.v(index, *)}
                        store_data, energy_name(i_en), data = str, dlim = dlim, lim = lim
                    ENDELSE 
                ENDELSE 
            ENDELSE 
        ENDELSE   
    ENDFOR  
;-----------------------------------------------------------------
;If there is no error then begin the calculatation
;---------------------------------------------------------------
    IF nerror EQ 0  THEN BEGIN
;-------------------------------------------------------------------
;Calculate total pressure & beta
;--------------------------------------------------------------
        IF sc EQ 4 THEN BEGIN 
            h1_press = h_press & o1_press = o_press
            
            plasma_beta, h1_press, mag_press, O1_PRESSURE = o1_press 
            
       ;     tplot_names, 'TDMOM*SP'+sp_str+'*', names = names
       ;     store_data, delete = names 
        ENDIF ELSE BEGIN 
            IF sc EQ 1 OR sc EQ 3 THEN BEGIN 
                get_data, hia_press, data = data
                store_data, hia_press, data = {x:data.x, y: 2*[[data.y(*, 1)], [data.y(*, 0)], [data.y(*, 0)]]}
                plasma_beta, hia_press, mag_press    
            ENDIF ELSE stop
        ENDELSE 

;------------------------------------------------------------------------
;Since it is possible prepossed data have different energy bins with   
;recalculated data. Here, we give a check. If the energybin in 
;flux data and counts data are different, recalculate flux. 
;----------------------------------------------------------------------
;Check runs for both tailward and earthward
;------------------------------------------------------------------------   
        FOR id = 0, 1  DO BEGIN 
            phi_str = phi_str_set(id)
            angle = angle_set(*, *, id)

            flux_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+phi_str+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All'
            counts_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+phi_str+'_UNCOUNTS_SP'+sp_str+'_ET0_All'
            
            get_data, flux_name,  data = data_flux, dlim = dlimf, lim = limf
            get_data, counts_name,  data = data_counts, dlim = dlimc, lim = limc

;if prepossed file have different energy channels then 
;recalc energy in flux and resave data in same time
;range as before          
            IF N_ELEMENTS(data_flux.v(0, *)) NE N_ELEMENTS(data_counts.v(0, *)) OR $
              N_ELEMENTS(data_flux.x) NE  N_ELEMENTS(data_counts.x) THEN BEGIN 
                sat = [sc]
                specie = [sp]
                
                inst = inst_input & units_name = 'DIFF FLUX' & eff_table = 0
                plot_en_spec_from_crib, sat, specie, inst,  $
                  units_name, angle, eff_table, recalc = 1
                
                get_data, flux_name,  data = data_flux, dlim = dlimf, lim = limf

                index = where(data_flux.x GT (t_s+1) AND data_flux.x LE t_e)
                str = {x:data_flux.x(index), y:data_flux.y(index, *), $
                       v:data_flux.v(index, *)}
                store_data, flux_name, data = str, dlim = dlimf, lim = limf
            ENDIF
        ENDFOR  

;------------------------------------------------------------------------
;Now, the main part of the program:
;= > clean up the low counts data from averaged energy spectra
;= > find Energy peak from filted energy spectra 
;= > plot pitch angle around the them  
;= > find the pitch angle peak
;= > filter the beam out by cleaning up the uncontineous pitch angle 
;-----------------------------------------------------------------------
;Do the Calculation seperately for tailward and earthward
;------------------------------------------------------------------------   
        FOR id = 0, 1 DO BEGIN 
            phi_str = phi_str_set(id)
            angle = angle_set(*, *, id)

;-----------------------------------------------------------------------------
;Clean up the low counts data in flux data, for energy data in both dirctions

; get data from loaded data        
            flux_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+phi_str+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All'
            counts_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+phi_str+'_UNCOUNTS_SP'+sp_str+'_ET0_All'

            get_data, flux_name,  data = data_flux, dlim = dlimf, lim = limf
            get_data, counts_name,  data = data_counts, dlim = dlimc, lim = limc

; record the original 4 sec time and energy info 
            index = where(data_flux.x GT(t_s+1) AND data_flux.x LE t_e)
            time_original = data_flux.x(index)
            start_time = time_original(0)
            end_time = time_original(N_ELEMENTS(time_original)-1)
            energybins = data_flux.v(0, *)
            nenergy = N_ELEMENTS(energybins)    
            ntime = floor((end_time - start_time)/average_time)

; average the flux data into new name and sumup the countsx data into new name
            average_tplot_variable, flux_name, at_str, /new
            average_tplot_variable, counts_name, at_str, /sumup, /new

; save the data into arries
            tplot_names, flux_name+'_AVG'+at_str, names = names
            IF names(0) EQ  '' THEN stop            
            get_data, flux_name+'_AVG'+at_str, data = data
            time_avg = data.x
            flux_flux = data.y
            energy_avg = data.v
            n_avg = N_ELEMENTS(time_avg) 

            tplot_names, counts_name+'_AVG'+at_str, names = names
            IF names(0) EQ  '' THEN stop             
            get_data, counts_name+'_AVG'+at_str, data = data
            counts_counts = data.y        
            
; save the original flux data into another string
            str = {x:time_avg, y:flux_flux, v:energy_avg}
            store_data, flux_name+'_AVG'+at_str+'_OLD', data = str, dlim = dlimf, lim = limf 
            
; filtered all the data which have counts lower than low_counts_line
; also clean the flux in the last energybin since the compression problem (or 2
; bins if energy bin # is 31)
; and store them in the original name

            index = where(counts_counts LT low_counts_line, ct)
            IF ct GT 0 THEN  flux_flux(index) = 0
            flux_flux(*, N_ELEMENTS(flux_flux(0, *))-(nenergy/15)) = 0
            
            str = {x:time_avg, y:flux_flux, v:energy_avg}
            store_data, flux_name+'_AVG'+at_str, data = str, dlim = dimf, lim = limf
            
;draw all those counts and flux plots into ps file 
            IF KEYWORD_SET(plot_lowcount_filter) THEN BEGIN 
                filepath = path+'plots/enct/'
                spwan = 'mkdir '+ filepath
                popen, filepath+'enct_' + date_s + '_' + time_s + '.ps'
                tplot, [flux_name, flux_name+'_AVG'+at_str+'_OLD', $
                        counts_name, counts_name+'_AVG'+at_str, $
                        flux_name+'_AVG'+at_str ]     
                pclose
            ENDIF 
;------------------------------------------------------------------- 
;mark magnetosheath and solar wind data 
;with h_density & h_velocity and x_gse & y_gse
;-----------------------------------------------------------
; get H+ Density and H+ V from saved string 
            tplot_names, 'location', names = names
            IF names(0) EQ ''  THEN BEGIN 
                h_density = fltarr(n_avg)
                h_velocity = fltarr(n_avg)
                get_data, 'TDMOM_EN00040_40000_SC'+sc_str+'_MTDENSITY_SP0_ET0_All', data = data_n
                get_data, 'TDMOM_EN00040_40000_SC'+sc_str+'_MTVELOCITY_SP0_ET0_All', data = data_v     
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data_n.x GE time_avg(itime)-average_time/2 $
                                  AND data_n.x LT time_avg(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data_n.y(index)) GE 0) GT 0  THEN $
                          h_density(itime) = total(data_n.y(index), /NAN)/ct $
                        ELSE h_density(itime) = !VALUES.F_NAN
                    ENDIF ELSE  h_density(itime) = !VALUES.F_NAN
                    
                    index = where(data_v.x GE time_avg(itime)-average_time/2 $
                                  AND data_v.x LT time_avg(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data_v.y(index, 0)) GE 0) GT 0 THEN $
                          h_velocity(itime) = sqrt((total(data_v.y(index, 0), /NAN))^2 $
                                                   + (total(data_v.y(index, 1), /NAN))^2 $
                                                   + (total(data_v.y(index, 2), /NAN))^2)/ct $
                        ELSE h_velocity(itime) = !VALUES.F_NAN
                    ENDIF ELSE h_velocity(itime) = !VALUES.F_NAN
                ENDFOR 
                store_data, 'TDMOM_EN00040_40000_SC'+sc_str+'_MTDENSITY_SP0_ET0_All_AVG'+at_str, $
                  data = {x:  time_avg, y:h_density}
                store_data, 'TDMOM_EN00040_40000_SC'+sc_str+'_MTVELOCITY_SP0_ET0_All_AVG'+at_str, $
                  data = {x:  time_avg, y:h_velocity}
; get plasma beta
                tplot_names, beta_name+'_AVG'+at_str, names = names
                IF names(0) EQ '' THEN  BEGIN 
                    get_data, beta_name, data = data_b
                    beta = FLTARR(n_avg)
                    FOR itime = 0l, n_avg-1 DO BEGIN
                        index = where(data_b.x GE time_avg(itime)-average_time/2 $
                                      AND data_b.x LT time_avg(itime)+average_time/2, ct)
                        IF ct GT 0 THEN BEGIN 
                            IF TOTAL(ABS(data_b.y(index)) GE 0) GT 0  THEN $
                              beta(itime) = total(data_b.y(index), /NAN)/ct $
                            ELSE beta(itime) = !VALUES.F_NAN
                        ENDIF ELSE beta(itime) = !VALUES.F_NAN
                    ENDFOR
                    store_data, beta_name+'_AVG'+at_str, data = {x:time_avg, y:beta}
                    beta=abs(beta)
                ENDIF 

; get x,y gse data from saved string
                get_data, 'EPH_SC'+ sc_str+'_GSE_X', data = data_gse_x 
                get_data, 'EPH_SC'+ sc_str+'_GSE_Y', data = data_gse_y
                get_data, 'EPH_SC'+ sc_str+'_GSM_Z', data = data_gsm_z
                x_gse = FLTARR(n_avg)
                y_gse = FLTARR(n_avg)
                z_gsm = FLTARR(n_avg)
                FOR itime = 0l, n_avg-1 DO BEGIN
                    index = where(data_gse_x.x GE time_avg(itime)-average_time/2 $
                                  AND data_gse_x.x LT time_avg(itime)+average_time/2, ct)
                    IF ct GT 0 AND TOTAL(ABS(data_gse_x.y(index)) GE 0) GT 0  THEN  $
                      x_gse(itime) = total(data_gse_x.y(index), /NAN)/ct $
                    ELSE x_gse(itime) = !VALUES.F_NAN
                    
                    index = where(data_gse_y.x GE time_avg(itime)-average_time/2 $
                                  AND data_gse_y.x LT time_avg(itime)+average_time/2, ct)
                    IF ct GT 0 AND TOTAL(ABS(data_gse_y.y(index)) GE 0) GT 0  THEN $
                      y_gse(itime) = total(data_gse_y.y(index), /NAN)/ct $
                    ELSE y_gse(itime) = !VALUES.F_NAN

                    index = where(data_gsm_z.x GE time_avg(itime)-average_time/2 $
                                  AND data_gsm_z.x LT time_avg(itime)+average_time/2, ct)
                    IF ct GT 0 AND TOTAL(ABS(data_gsm_z.y(index)) GE 0) GT 0  THEN $
                      z_gsm(itime) = total(data_gsm_z.y(index), /NAN)/ct $
                    ELSE z_gsm(itime) = !VALUES.F_NAN
                    
                ENDFOR  
                store_data, 'EPH_SC'+ sc_str+'_GSE_X_AVG'+at_str, data = {x: time_avg, y:x_gse}
                store_data, 'EPH_SC'+ sc_str+'_GSE_Y_AVG'+at_str, data = {x: time_avg, y:y_gse}   
                
;sheath is defined as H+ density gt 3 and velocity gt 65
;or for dayside polar region (x_gse gt 1 Re and ABS(z_gsm) gt 5) defined as plasma beta gt 0.05
;or for hig latitude region (x_gse le 1 and ABS(z_gsm) gt 10) defined as plasma beta gt 1
                locate_sheath = float((h_density GT 3 and h_velocity GT 65) $
                                      or (x_gse gt 1 and ABS(z_gsm) gt 5 and beta gt 0.05) $
                                      or (x_gse le 1 and ABS(z_gsm) gt 10 and beta gt 1))
                locate_no_sheath = ABS(1-(locate_sheath))
                
; sw is definded with location as found by observing the map
; which is: x_gse > = -1 and outside an ellipse deceided by observations 
                locate_no_sw = float(x_gse LE -1) OR  (((x_gse+1)^2/8.^2+(y_gse-1)^2/14.^2) LE 1)
                locate_magnetosphere = float((locate_no_sw + locate_no_sheath) EQ 2)
                
                magnetosphere_data_index = where(locate_magnetosphere EQ 1, ct_magnetosphere)
                index = where(locate_magnetosphere EQ 0, ct)
                IF ct GT 0 THEN locate_magnetosphere(index) = !VALUES.F_NAN
                index = where(locate_no_sheath EQ 1, ct)
                IF ct GT 0 THEN locate_no_sheath(index) = !VALUES.F_NAN
                index = where(locate_no_sw EQ 1, ct)
                IF ct GT 0 THEN locate_no_sw(index) = !VALUES.F_NAN
                
                store_data, 'location', $
                  data = {x:time_avg, $
                          y: [[locate_no_sw], [locate_no_sheath], [locate_magnetosphere]]}, $
                  lim = {yrange:[-1, 2], psym:10}
            ENDIF       
;use beta to deceide regions in magnetosphere   
            get_data, 'location', data = data
            locate_magnetosphere = data.y(*, 2)

            IF inst_input EQ 1 or keyword_set(only_in_lobe) THEN BEGIN 
                get_data, beta_name+'_AVG'+at_str, data = data
                beta = data.y
                locate_mag = (beta LT 0.05)*locate_magnetosphere 
            ENDIF ELSE locate_mag = locate_magnetosphere
            index = where(locate_mag EQ 0, ct)
            IF ct GT 0 THEN locate_mag(index) = !VALUES.F_NAN
            if not keyword_set(plot_sw_sheath) then begin 
                FOR itime = 0, n_avg-1 DO flux_flux(itime, *) = flux_flux(itime, *)*locate_mag(itime)
            endif 

            str = {x:time_avg, y:flux_flux, v:energy_avg}
            store_data, flux_name+'_AVG'+at_str, data = str, dlim = dimf, lim = limf

;-----------------------------------------------------------------
;find the evergy range and energy peak from average energy spectra
            en_name = flux_name
            en_name_avg = en_name + '_AVG' + at_str
            
            find_energy_range_from_enspec, en_name_avg, epcut, en_range
;------------------------------------------------------
; plot pitch angle with routine plot_pa_spec_around_energy_peak

            sat  = [sc]
            specie = [sp]
            units_name = 'DIFF FLUX' & inst = inst_input & eff_table = 0
            
            plot_pa_spec_around_energy_peak, sat, specie, inst, $
              units_name, eff_table, $
              invar = epcut, $
              outvar = pa_name, $
              average_time = average_time, $
              start_time = start_time, $
              END_time = END_time, $
              n_range = 1, $
              PaBin = 22.5
            
            sat  = [sc]
            specie = [sp]
            units_name = 'Counts' & inst = inst_input & eff_table = 0
            
            plot_pa_spec_around_energy_peak, sat, specie, inst, $
              units_name, eff_table, $
              invar = epcut, $
              outvar = pa_name_counts, $
              average_time = average_time, $
              start_time = start_time, $
              END_time = END_time, $
              n_range = 1, $
              PaBin = 22.5
            
            find_pa_peak, pa_name_counts, pa_name, pap_name, pa_counts_line = pa_counts_line,flux_threshold=flux_threshold
            
            beam_filter, sc, pap_name, epcut, et_beam, epcut_beam
        ENDFOR   
;---------------------------------------------------
; combine tail and earth beam results
        tail_beam_et_name = 'PASPEC_SC'+sc_str+ $
          '_IN'+inst_str+'_'+phi_str_set(0)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG' $ 
          +at_str +'_PAP_ET_beam'
        earth_beam_et_name = 'PASPEC_SC'+sc_str $
          +'_IN'+inst_str+'_'+phi_str_set(1)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG' $
          +at_str +'_PAP_ET_beam'
        tail_epcut_beam_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+phi_str_set(0)+''  $
          +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str $
          +'_epcut_beam'
        earth_epcut_beam_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+phi_str_set(1)+''  $
          +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str $
          +'_epcut_beam'
        tail_erange_beam_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+phi_str_set(0)+''  $
          +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str $
          +'_erange'
        earth_erange_beam_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+phi_str_set(1)+''  $
          +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str $
          +'_erange'

        combine_et_pap, sc, $
          tail_beam_et_name, earth_beam_et_name, $
          pap_beam_combine_et, pap_beam_combine_pa, $
          tail_epcut_beam_name, earth_epcut_beam_name, $
          tail_erange_beam_name, earth_erange_beam_name, $
          start_time = start_time, END_time = END_time, $
          average_time = average_time
        
        OPENU, unit, path+'log_plotted.txt', /GET_LUN, /APPEND
        PRINTF, unit, TIME_STRING(t_s)+' TO '+ TIME_STRING(t_e)+ $
          '   --------Find O+ Beam------'
        FREE_LUN, unit
    ENDIF  ELSE BEGIN 
        OPENU, unit, path+'log_errors.txt', /GET_LUN, /APPEND       
        PRINTF, unit, TIME_STRING(t_s) + ' TO '+ TIME_STRING(t_e), s_e
        FREE_LUN, unit         
    ENDELSE        
ENDIF         

;-----------------------------------------------------
; PLOT Moments 
;------------------------------------------------
;check whether beam has already been 
;-------------------------------------------------
IF ct_beam GE 0  THEN BEGIN 
    beam_name = 'PASPEC_SC'+sc_str+ $
      '_IN'+inst_str+'_PHICOMBINED_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG*_PAP_ET_beam'
    tplot_names, beam_name, names = names

    IF names(0) NE  '' THEN BEGIN     
        get_data, names(0), data = data
        average_time = data.average_time
        start_time = data.start_time
        END_time = data.end_time
        ntime = floor((END_time-start_time)/average_time)
        time_avg = data.x
        n_avg = N_ELEMENTS(time_avg)

;------------------------------------------------------------
;plot mom for both direction
;-------------------------------------------------------------------
        IF KEYWORD_SET(mom_recalc)  THEN BEGIN     
            FOR id = 0, 1 DO BEGIN 
                phi_str = phi_str_set(id)

                pa_beam_name = 'PASPEC_SC'+sc_str+'_IN'+inst_str+'_'+phi_str+ $
                  '_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str+'_PAP_PA_beam'
                
                en_range_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str $
                  +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str $
                  +'_erange'
                sat = sc
                specie = sp
                inst = inst_input & eff_table = 0
                units_name = 'DIFF FLUX'
                
                find_bins_from_pap, sat, specie, inst, units_name, eff_table, $
                  average_time, pa_beam_name, bins_name, $
                  start_time = start_time, END_time = END_time, $
                  beam_angle_range = beam_angle_range

                get_data, en_range_name, data = data
                energy_range_time = data.x
                energy_range = data.y
                energy_range_bins = data.energybins
                
                get_data, bins_name, data = data
                time_bins = data.x
                angle_bins = data.y

                ntime = floor((end_time-start_time)/average_time)
                n_avg = N_ELEMENTS(energy_range_time) 
                
                time_Tt = FLTARR(n_avg) &  data_Tt = FLTARR(n_avg)
                time_D = FLTARR(n_avg) & data_D = FLTARR(n_avg)
                time_Vt = FLTARR(n_avg) &  data_Vt = FLTARR(n_avg)
                time_Pt = FLTARR(n_avg) &  data_Pt = FLTARR(n_avg)  
                time_V = FLTARR(n_avg) &  data_V = FLTARR(n_avg, 3)
                time_Tx = FLTARR(n_avg) &  data_Tx = FLTARR(n_avg)
                time_Ty = FLTARR(n_avg) &  data_Ty = FLTARR(n_avg)
                time_Tz = FLTARR(n_avg) &  data_Tz = FLTARR(n_avg)

                t_dfit =  FLTARR(n_avg, 2)   &  t_error =  FLTARR(n_avg, 2)
                d_dfit = FLTARR(n_avg)    &  d_error = FLTARR(n_avg)
                FOR i_m = 0, ntime-1 DO BEGIN
                    
                    time = start_time+i_m*average_time
                    index = where(time_bins GT time AND time_bins LE time+average_time, ct)
                    loc = index(0)
                    IF ct EQ 1 THEN BEGIN 
                        timespan, time, average_time, /SECONDS
                        moments = ['A']
                                ;                e_r = [FLOOR(energy_range(loc, 0)),  $
                                ;                      CEIL(energy_range(loc, 1))] 
                                ;               a_r =  REFORM(angle_bins(loc, *))
                        use_bins = use_angle_range
                        IF keyword_set(use_angle_range) THEN $
                          bins = REFORM(angle_bins(loc, *)) $
                        ELSE bins = REPLICATE(1, 88)
                        
                        IF KEYWORD_SET(use_energy_range) THEN $
                          energy = [FLOOR(energy_range(loc, 0)), CEIL(energy_range(loc, 1))] $
                        ELSE energy = [40., 40000.]
                        
                        index = where(bins GT 0, ct1)
                        index = where(energy_range(loc, *)GE 0, ct2)

                        IF ct1*ct2 EQ 0 AND ct1+ct2 GT 0 THEN stop

                        IF ct1 GT 0 AND ct2 GT 0 THEN BEGIN 
                            IF keyword_set(globe_plot)THEN BEGIN 
                                sat = sc
                                specie = sp
                                bins_input = bins
                                en_range = energy
                                globe_plot, sat = sat, specie = specie, $
                                  units_name = units_name, $
                                  bins = bins_input, energy = en_range, $
                                  path = path, ps = 1
                            ENDIF 
                            IF keyword_set(dfit_temperature) THEN BEGIN 
                                sat = sc
                                specie = sp
                                bins_input = bins
                                en_range = energy
                                IF NOT keyword_set(show_fit) THEN path_input = path
                                t_calc_from_disf, sat = sat, specie = specie, $
                                  energy = en_range, $
                                  bins = bins_input, $
                                  t_dfit_name = t_dfit_name, path = path_input

                                get_data, t_dfit_name, data = data_fit
                                t_dfit(loc(0), *) = [data_fit.t(0), data_fit.t(1)]
                                t_error(loc(0), *) = [data_fit.t_error(0), data_fit.t_error(1)]
                                d_dfit(loc(0)) = [data_fit.n]
                                d_error(loc(0)) = [data_fit.n_error]
                                store_data, delete = t_dfit_name
                            ENDIF 

                            sat = [sc]
                            specie = [sp]
                                ; a_r,e_r here are just for verify that
                                ; using bins and energy is correct or not
                            units_name = 'DIFF FLUX' 
                            inst = inst_input & eff_table = 0
                            angle = [[-90.0, 90.0], [0., 360.]] 

                            plot_3dmom_from_crib, sat, specie, inst, $
                              moments, angle, energy, eff_table, recalc = 1, $
                              bins = bins, use_bins = use_bins, sum_up = 1 ;, $
                                ;                      e_r = e_r, a_r = a_r
                            
                            e_min = STRCOMPRESS(STRING(energy(0), $
                                                       FORMAT = '(i5.5)'), /REMOVE_ALL)
                            e_max = STRCOMPRESS(STRING(energy(1), $
                                                       FORMAT = '(i5.5)'), /REMOVE_ALL)
                            
                            Tt_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                              +'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T'
                            Vt_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                              + '_MTVELOCITY_SP'+sp_str+'_ET0_All_T'
                            D_name = 'TDMOM_EN'+ e_min + '_' + e_max + '_SC' + sc_str  $
                              +'_MTDENSITY_SP'+sp_str+'_ET0_All'
                            Pt_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                              + '_MTPRESSURE_SP'+sp_str+'_ET0_All_T'
                            V_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                              + '_MTVELOCITY_SP'+sp_str+'_ET0_All'
                            Tx_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                              +'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_X'
                            Ty_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                              +'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Y'
                            Tz_name = 'TDMOM_EN' + e_min + '_' + e_max + '_SC' + sc_str $
                              +'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Z'
                            
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
                                time_Tt(loc) = data.x(0)
                                data_Tt(loc) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                
                                get_data, D_name, data = data, dlim = dlim_D, lim = lim_D
                                time_D(loc) = data.x(0)
                                data_D(loc) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                
                                get_data, Vt_name, data = data, dlim = dlim_Vt, lim = lim_Vt
                                time_Vt(loc) = data.x(0)
                                data_Vt(loc) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                
                                get_data, Pt_name, data = data, dlim = dlim_Pt, lim = lim_Pt
                                time_Pt(loc) = data.x(0)
                                data_Pt(loc) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                
                                get_data, V_name, data = data, dlim = dlim_V, lim = lim_V
                                time_V(loc) = data.x(0)
                                data_V(loc, *) = TOTAL(data.y(*, *), 1)/N_ELEMENTS(data.y(*, 0))
                                
                                get_data, Tx_name, data = data, dlim = dlim_Tx, lim = lim_Tx
                                time_Tx(loc) = data.x(0)
                                data_Tx(loc) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                
                                get_data, Ty_name, data = data, dlim = dlim_Ty, lim = lim_Ty
                                time_Ty(loc) = data.x(0)
                                data_Ty(loc) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                
                                get_data, Tz_name, data = data, dlim = dlim_Tz, lim = lim_Tz
                                time_Tz(loc) = data.x(0)
                                data_Tz(loc) = TOTAL(data.y)/N_ELEMENTS(data.y)
                                
                                tplot_names, 'TDMOM_EN'+e_min+'_'+e_max+'*SP'+sp_str+'*', names = names
                                store_data, delete = names
                            ENDIF  ELSE    invalid_mom_data = 1
                        ENDIF  ELSE    invalid_mom_data = 1

                        IF keyword_set(invalid_mom_data) THEN BEGIN     
                            energy_range(loc, *) = !VALUES.F_NAN
                            time_tt(loc) = time_bins(loc)
                            data_Tt(loc) = !VALUES.F_NAN
                            time_d(loc) = time_bins(loc)
                            data_d(loc) = !VALUES.F_NAN
                            time_vt(loc) = time_bins(loc)
                            data_vt(loc) = !VALUES.F_NAN
                            time_pt(loc) = time_bins(loc)
                            data_pt(loc) = !VALUES.F_NAN
                            time_v(loc) = time_bins(loc)
                            data_v(loc,*) = !VALUES.F_NAN
                            time_tx(loc) = time_bins(loc)
                            data_tx(loc) = !VALUES.F_NAN
                            time_ty(loc) = time_bins(loc)
                            data_ty(loc) = !VALUES.F_NAN
                            time_tz(loc) = time_bins(loc)
                            data_tz(loc) = !VALUES.F_NAN
                            invalid_mom_data = 0
                            t_dfit(loc) = !VALUES.F_NAN
                            t_error(loc) =  !VALUES.F_NAN
                            d_dfit(loc) = !VALUES.F_NAN
                            d_error(loc) =  !VALUES.F_NAN
                        ENDIF     
                    ENDIF 
                    tplot_names, 'TDMOM_EN0*SP'+sp_str+'*', names = names
                    store_data, delete = names        
                    store_data, delete = 'B_xyz*'
                ENDFOR                    
                timespan, t_s, t_dt, /SECONDS  
                
                Tt_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                  +phi_str+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T'+'_AVG'+at_str
                Vt_name = 'TDMOM_ENVARIOUS'+'_SC' + sc_str+'_' $
                  +phi_str + '_MTVELOCITY_SP'+sp_str+'_ET0_All_T'+'_AVG'+at_str
                D_name =  'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_'  $
                  +phi_str +'_MTDENSITY_SP'+sp_str+'_ET0_All'+'_AVG'+at_str
                Pt_name = 'TDMOM_ENVARIOUS'+'_SC' + sc_str+'_' $
                  +phi_str +'_MTPRESSURE_SP'+sp_str+'_ET0_All_T'+'_AVG'+at_str
                V_name = 'TDMOM_ENVARIOUS'+'_SC' + sc_str+'_' $
                  +phi_str + '_MTVELOCITY_SP'+sp_str+'_ET0_All'+'_AVG'+at_str
                Tx_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                  +phi_str+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_X'+'_AVG'+at_str
                Ty_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                  +phi_str+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Y'+'_AVG'+at_str
                Tz_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                  +phi_str+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Z'+'_AVG'+at_str
                
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
                index = where (data_tt LT 0, ct)
                IF ct GT 0 THEN BEGIN                   
                    OPENU, unit, path+'log_neg_T.txt', /GET_LUN, /APPEND
                    PRINTF, unit,  ts+' - '+te, '-----negative Tt and Pt-----'+phi_str
                    FREE_LUN, unit  
                ENDIF 
                
                IF keyword_set(dfit_temperature) THEN BEGIN 
                    Tt_fit_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +phi_str+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_para_dfit'+'_AVG'+at_str
                    store_data,  tt_fit_name, data = {x:time_tt, $
                                                      y:t_dfit(*, 0), dy:t_error(*, 0)}, $
                      dlim = dlim_Tt, lim = lim_Tt
                    options, Tt_fit_name, 'ytitle', 'SC'+sc_str+' O!U+!N!C!CT!D//!N (eV)'

                    Tt_fit_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +phi_str+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_perp_dfit'+'_AVG'+at_str
                    store_data,  tt_fit_name, data = {x:time_tt, $
                                                      y:t_dfit(*, 1), dy:t_error(*, 1)}, $
                      dlim = dlim_Tt, lim = lim_Tt
                    options, Tt_fit_name, 'ytitle', 'SC'+sc_str+' O!U+!N!C!CT!DL!N (eV)'

                    Tt_fit_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +phi_str+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T_dfit'+'_AVG'+at_str
                    store_data,  tt_fit_name, $
                      data = {x:time_tt, $
                              y:(t_dfit(*, 0)+2*t_dfit(*, 1))/3, dy:t_error}, $
                      dlim = dlim_Tt, lim = lim_Tt
                    options, Tt_fit_name, 'ytitle', 'SC'+sc_str+' O!U+!N!C!CT!Davg!N (eV)'

                    d_fit_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                      +phi_str+'_MTDENSITY_SP'+sp_str+'_ET0_All_dfit'+'_AVG'+at_str
                    store_data,  d_fit_name, data = {x:time_d, $
                                                     y:d_dfit, dy:d_error}, $
                      dlim = dlim_d, lim = lim_d
                    options, d_fit_name, 'ytitle',  'SC'+sc_str+'!C!CO!U+!N!C!Cn (cm!U-3!N)' 
                ENDIF 
;----- calculate V_perp and V_par
                v_perp, V_name               
            ENDFOR  
;---------------------------------------------------
;write plotted info into the plotted log or
;write data errors info into log
;-----------------------------------------------------------  
            OPENU, unit, path+'log_plotted.txt', /GET_LUN, /APPEND
            PRINTF, unit, ts+' TO '+te+ $
              '   -------mom plotted---------'
            FREE_LUN, unit
        ENDIF     
;-------------------------------------------------------------------
;Adding Storm phase information
;Read storm phase(prestorm, storm time and recovery time) form 
;file 'storm_phase.dat' and store them into arrays 
;Deceide the storm phase for each data
;-------------------------------------------------------------------
        IF KEYWORD_SET(find_phase) THEN BEGIN 
            OPENR, unit, 'storm_phase_long.dat', /GET_LUN
            prestorm_start = DBLARR(300)
            storm_onset = DBLARR(300) 
            min_dst_new = DBLARR(300)
            recovery_fast_end = DBLARR(300) 
            recovery_early_end = DBLARR(300) 
            recovery_long_end = DBLARR(300)

            jj = 0l
            dummy = ''
            WHILE NOT EOF(unit) DO BEGIN
                
                READF, unit, dummy
                IF STRMID(dummy, 0, 5) NE 'Start' THEN BEGIN 
                    prestorm_start(jj) = time_double(STRMID(dummy, 0, 20))
                    storm_onset(jj) = time_double(STRMID(dummy, 20, 20))
                    min_dst_new(jj) =  time_double(STRMID(dummy, 40, 24))
                    recovery_fast_end(jj) =  time_double(STRMID(dummy, 60, 24))
                    recovery_early_end(jj) =  time_double(STRMID(dummy, 80, 24))
                    recovery_long_end(jj) = time_double(STRMID(dummy, 100, 24))

                    jj = jj+1      
                ENDIF  
            ENDWHILE
            CLOSE, unit, /all
            nstorm = jj
            
            prestorm_start = prestorm_start(0:nstorm-1)
            storm_onset = storm_onset(0:nstorm-1)
            min_dst_new = min_dst_new(0:nstorm-1)
            recovery_fast_end = recovery_fast_end (0:nstorm-1)
            recovery_early_end =  recovery_early_end (0:nstorm-1)
            recovery_long_end =  recovery_long_end (0:nstorm-1)

            storm_phase = INTARR(n_avg)
            FOR  itime = 0, n_avg-1 DO BEGIN
                FOR istorm = 0, nstorm-1 DO BEGIN 
                    belong = 0
                    IF time_avg(itime) GE prestorm_start(istorm) AND $
                      time_avg(itime) LT storm_onset(istorm) THEN BEGIN 
                        storm_phase(itime) = 1 ;initial phase
                        belong = belong+1
                    ENDIF 
                    IF time_avg(itime) GE storm_onset(istorm) AND $
                      time_avg(itime)  LT min_dst_new(istorm) THEN BEGIN 
                        storm_phase(itime) = 2 ; mian phase
                        belong = belong+1
                    ENDIF 
                    IF time_avg(itime) GE min_dst_new(istorm) AND $
                      time_avg(itime) LT recovery_early_end(istorm) THEN BEGIN 
                        storm_phase(itime) = 3 ; recovery phase
                        belong = belong+1
                    ENDIF
                    IF time_avg(itime) GE recovery_early_end(istorm) AND $
                      time_avg(itime) LT recovery_long_end(istorm) THEN BEGIN 
                        storm_phase(itime) = 4 ; later recovery phase (not used for any map)
                        belong = belong+1
                    ENDIF

                    IF belong GT 1 THEN stop
                    if belong eq 0 then begin 
                        IF time_avg(itime) GE (prestorm_start(istorm)-60.*60.*4) AND $
                          time_avg(itime) LT prestorm_start(istorm) THEN BEGIN 
                            storm_phase(itime) = 5 ; pre storm
                            belong = belong+1
                        endif 
                    endif  
                ENDFOR 
            ENDFOR  
            store_data, 'storm_phase', data = {x:time_avg, y:storm_phase}
            ylim, 'storm_phase', -1, 6
        ENDIF   

;---------------------------------
;Adding IMF and solar wind information 
;------------------------
        IF KEYWORD_SET(Add_IMF) THEN BEGIN 
;read IMF from ace 
            plot_ace_swepam, /ALL ;Load ACE SWEPAM
            plot_ace_mag, /ALL  ; Load ACE MAG

;Calculate SW delay to subsolar point (or cluster position if required)
            for i_delay=0, n_delay-1 do begin
                ACE_x_time = r_data('ACE_SWEPAM_X_GSM', /X)
                ACE_x_data = r_data('ACE_SWEPAM_X_GSM', /Y)
                SW_velocity_time = r_data('ACE_SWEPAM_VELOCITY', /X)
                SW_velocity_data = r_data('ACE_SWEPAM_VELOCITY', /Y)
                ACE_x_data = INTERPOL(ACE_x_data, ACE_x_time, SW_velocity_time)
                sw_delay_time = SW_velocity_time
                
                if i_delay eq 0 then begin
                                ;delay to CLUSTER position if required
                    CLUSTER_x_time = r_data('EPH_SC'+ sc_str+'_GSE_X',/X)
                    CLUSTER_x_data = r_data('EPH_SC'+ sc_str+'_GSE_X',/Y)
                    CLUSTER_x_data = INTERPOL(CLUSTER_x_data,CLUSTER_x_time,SW_velocity_time)
                    sw_delay = (ACE_x_data - CLUSTER_x_data) / SW_velocity_data
                    delay_location='_DC'
                endif else begin
                    subsolar_point = 10. * 6370.0 ; km  ;set subsolar point at 10 Re            
                    sw_delay = (ACE_x_data - subsolar_point) / SW_velocity_data
                    delay_location=''
                endelse 
                                ; Apply delay to data
                ace_swepam_velocity_time = r_data('ACE_SWEPAM_VELOCITY', /X)
                ace_swepam_velocity_data = r_data('ACE_SWEPAM_VELOCITY', /Y)
                sw_delay_n = INTERPOL(sw_delay, sw_delay_time, ace_swepam_velocity_time)
                ace_swepam_velocity_time = ace_swepam_velocity_time + SW_delay_n
                sort_SW_velocity_time = SORT(ace_swepam_velocity_time)
                ace_swepam_velocity_time = ace_swepam_velocity_time(sort_SW_velocity_time)
                ace_swepam_velocity_data = ace_swepam_velocity_data(sort_SW_velocity_time)
                store_data, 'ACE_SWEPAM_VELOCITY'+delay_location, $
                  data = {x:ace_swepam_velocity_time, y:ace_swepam_velocity_data}
                
                ace_swepam_pressure_time = r_data('ACE_SWEPAM_PRESSURE', /X)
                ace_swepam_pressure_data = r_data('ACE_SWEPAM_PRESSURE', /Y)
                sw_delay_n = INTERPOL(sw_delay, sw_delay_time, ace_swepam_pressure_time)
                ace_swepam_pressure_time = ace_swepam_pressure_time + SW_delay_n
                sort_SW_pressure_time = SORT(ace_swepam_pressure_time)
                ace_swepam_pressure_time = ace_swepam_pressure_time(sort_SW_pressure_time)
                ace_swepam_pressure_data = ace_swepam_pressure_data(sort_SW_pressure_time)
                store_data, 'ACE_SWEPAM_PRESSURE'+delay_location, $
                  data = {x:ace_swepam_pressure_time, y:ace_swepam_pressure_data}
                
                ace_mag_gsm_z_time = r_data('ACE_MAG_GSM_Z', /X)
                ace_mag_gsm_z_data = r_data('ACE_MAG_GSM_Z', /Y)
                sw_delay_n = INTERPOL(sw_delay, sw_delay_time, ace_mag_gsm_z_time)
                ace_mag_gsm_z_time = ace_mag_gsm_z_time + SW_delay_n
                sort_SW_gsm_z_time = SORT(ace_mag_gsm_z_time)
                ace_mag_gsm_z_time = ace_mag_gsm_z_time(sort_SW_gsm_z_time)
                ace_mag_gsm_z_data = ace_mag_gsm_z_data(sort_SW_gsm_z_time)
                store_data, 'ACE_MAG_GSM_Z'+delay_location, $
                  data = {x:ace_mag_gsm_z_time, y:ace_mag_gsm_z_data}
                
                ace_mag_gsm_y_time = r_data('ACE_MAG_GSM_Y', /X)
                ace_mag_gsm_y_data = r_data('ACE_MAG_GSM_Y', /Y)
                sw_delay_n = INTERPOL(sw_delay, sw_delay_time, ace_mag_gsm_y_time)
                ace_mag_gsm_y_time = ace_mag_gsm_y_time + SW_delay_n
                sort_SW_gsm_y_time = SORT(ace_mag_gsm_y_time)
                ace_mag_gsm_y_time = ace_mag_gsm_y_time(sort_SW_gsm_y_time)
                ace_mag_gsm_y_data = ace_mag_gsm_y_data(sort_SW_gsm_y_time)
                store_data, 'ACE_MAG_GSM_Y'+delay_location, $
                  data = {x:ace_mag_gsm_y_time, y:ace_mag_gsm_y_data}
                
                ace_mag_gsm_y_time = r_data('ACE_MAG_GSM_X', /X)
                ace_mag_gsm_y_data = r_data('ACE_MAG_GSM_X', /Y)
                sw_delay_n = INTERPOL(sw_delay, sw_delay_time, ace_mag_gsm_y_time)
                ace_mag_gsm_y_time = ace_mag_gsm_y_time + SW_delay_n
                sort_SW_gsm_y_time = SORT(ace_mag_gsm_y_time)
                ace_mag_gsm_y_time = ace_mag_gsm_y_time(sort_SW_gsm_y_time)
                ace_mag_gsm_y_data = ace_mag_gsm_y_data(sort_SW_gsm_y_time)
                store_data, 'ACE_MAG_GSM_X'+delay_location, $
                  data = {x:ace_mag_gsm_y_time, y:ace_mag_gsm_y_data}
;reduce the ace data into the timespan interval
                get_data, 'ACE_SWEPAM_VELOCITY'+delay_location, data = data, dlim = dlim, lim = lim
                index = where(data.x GE ts AND data.y LT te, ct)
                IF ct GT 0 THEN BEGIN 
                    index_s = sort(ABS(data.x-t_s))
                    index_e = sort(ABS(data.x-t_e))
                    time = data.x(index_s(0):index_e(0))
                    data = data.y(index_s(0):index_e(0))
                ENDIF ELSE BEGIN 
                    time = (ts+te)/2
                    data = !VALUES.F_NAN
                ENDELSE 

                store_data, 'ACE_SWEPAM_VELOCITY'+delay_location,  data = {x:time, y:data}, $
                  dlim = dlim, lim = lim
                
                get_data, 'ACE_SWEPAM_PRESSURE'+delay_location, data = data, dlim = dlim, lim = lim
                index = where(data.x GE ts AND data.y LT te, ct)
                IF ct GT 0 THEN BEGIN 
                    index_s = sort(ABS(data.x-t_s))
                    index_e = sort(ABS(data.x-t_e))
                    time = data.x(index_s(0):index_e(0))
                    data = data.y(index_s(0):index_e(0))
                ENDIF ELSE BEGIN 
                    time = (ts+te)/2
                    data = !VALUES.F_NAN
                ENDELSE 
                store_data, 'ACE_SWEPAM_PRESSURE'+delay_location,  data = {x:time, y: data}, $
                  dlim = dlim, lim = lim
                
                get_data, 'ACE_MAG_GSM_X'+delay_location, data = data, dlim = dlim, lim = lim
                index = where(data.x GE ts AND data.y LT te, ct)
                IF ct GT 0 THEN BEGIN 
                    index_s = sort(ABS(data.x-t_s))
                    index_e = sort(ABS(data.x-t_e))
                    time = data.x(index_s(0):index_e(0))
                    data = data.y(index_s(0):index_e(0))
                ENDIF ELSE BEGIN 
                    time = (ts+te)/2
                    data = !VALUES.F_NAN
                ENDELSE 
                store_data, 'ACE_MAG_GSM_X'+delay_location,  data = {x:time, y: data}, $
                  dlim = dlim, lim = lim
                
                get_data, 'ACE_MAG_GSM_Y'+delay_location, data = data, dlim = dlim, lim = lim
                index = where(data.x GE ts AND data.y LT te, ct)
                IF ct GT 0 THEN BEGIN 
                    index_s = sort(ABS(data.x-t_s))
                    index_e = sort(ABS(data.x-t_e))
                    time = data.x(index_s(0):index_e(0))
                    data = data.y(index_s(0):index_e(0))
                ENDIF ELSE BEGIN 
                    time = (ts+te)/2
                    data = !VALUES.F_NAN
                ENDELSE 
                store_data, 'ACE_MAG_GSM_Y'+delay_location,  data = {x:time, y: data}, $
                  dlim = dlim, lim = lim
                
                get_data, 'ACE_MAG_GSM_Z'+delay_location, data = data, dlim = dlim, lim = lim
                index = where(data.x GE ts AND data.y LT te, ct)
                IF ct GT 0 THEN BEGIN 
                    index_s = sort(ABS(data.x-t_s))
                    index_e = sort(ABS(data.x-t_e))
                    time = data.x(index_s(0):index_e(0))
                    data = data.y(index_s(0):index_e(0))
                ENDIF ELSE BEGIN 
                    time = (ts+te)/2
                    data = !VALUES.F_NAN
                ENDELSE 
                store_data, 'ACE_MAG_GSM_Z'+delay_location,  data = {x:time, y: data}, $
                  dlim = dlim, lim = lim
            endfor 
        ENDIF  
 
;--------------------------
;Adding eflux information
;1. plot energy spectra in unit eflux 
;2. extract the eflux correcponding to the energy peak
;-------------------------- 
        IF KEYWORD_SET(add_eflux) THEN BEGIN 
            tplot_names,'*EFLUX*',names=names
            store_data, delete=names
            tplot_names, 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_*_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str +'_epcut_beam_OLD', names = names                
            IF names(0) NE  '' then begin 
                FOR id = 0, 1 DO BEGIN
                    phi_str = phi_str_set(id)
                    en_peak_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str +'_epcut_beam'
                    options,en_peak_name+'','ytitle', 'en peak!C'+phi_str
                    enspec_eflux_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str  +'_UNEFLUX_SP'+sp_str+'_ET0_All' 
                    eflux_peak_name =  'EFLUX_SC' + sc_str + '_IN'+inst_str+'_' + phi_str +'_UNEFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'
                    
                    sat = [sc]
                    specie = [sp]
                    angle = angle_set(*, *, id)
                    inst = inst_input & units_name = 'EFLUX' 
                    eff_table = 0
                    
                    plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 1
                    zlim, enspec_eflux_name, 1e3, 1e7
                    average_tplot_variable, enspec_eflux_name, average_time, /new
                    
                    get_data, enspec_eflux_name+'_AVG'+at_str, data = data, dlim = dlim, lim = lim
                    en_eflux_time = data.x
                    en_eflux_eflux = data.y
                    en_eflux_en = data.v

                    get_data, en_peak_name+'_OLD', data = data
                    en_peak_old_time = data.x
                    en_peak_old = data.y
                    
                    IF total(en_eflux_time NE en_peak_old_time) EQ 0 THEN BEGIN 
                        FOR iavg = 0, n_avg-1 DO  en_eflux_eflux(iavg, *) = en_eflux_eflux(iavg, *)*(en_eflux_en(iavg, *) EQ en_peak_old(iavg))                         
                        en_eflux_cut = total(en_eflux_eflux, 2,/nan)
                        index= where(en_eflux_cut eq 0,ct) & if ct gt 0 then en_eflux_cut(index) = !VALUES.F_NAN
                        store_data, enspec_eflux_name+'_AVG'+at_str+'_peak_OLD', data = {x:en_eflux_time, y:en_eflux_eflux, v:en_eflux_en}, dlim = dlim, lim = lim
                        store_data, eflux_peak_name+'_OLD', data = {x:en_eflux_time, y:en_eflux_cut}, dlim = {psym:-1,ytitle:'eflux old!C'+phi_str},lim={ylog:1}
                        
                        options,enspec_eflux_name+'_AVG'+at_str+'_peak_OLD','ytitle','energy(eV) old!C'+phi_str
                    ENDIF ELSE stop
                    
                ENDFOR
                
                en_peak_tail_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str_set(0) +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str +'_epcut_beam'
                en_peak_earth_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str_set(1) +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str +'_epcut_beam'              
                eflux_peak_tail_name =  'EFLUX_SC' + sc_str + '_IN'+inst_str+'_' + phi_str_set(0) +'_UNEFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'
                eflux_peak_earth_name =  'EFLUX_SC' + sc_str + '_IN'+inst_str+'_' + phi_str_set(1) +'_UNEFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'

                get_data,eflux_peak_tail_name+'_OLD',data=data
                eflux_time=data.x
                eflux_old_tail = data.y
     
                get_data,eflux_peak_earth_name+'_OLD',data=data
                eflux_old_earth = data.y
               
                get_data,en_peak_tail_name,data=data
                en_peak_tail=data.y
                get_data,en_peak_earth_name,data=data
                en_peak_earth=data.y
                
                eflux_tail=eflux_old_tail*(en_peak_tail/en_peak_tail)
                index=where( ~finite(eflux_old_tail) and en_peak_tail gt 0,ct )
                if ct gt 0 then eflux_tail(index)=eflux_old_earth(index)
                store_data, eflux_peak_tail_name, data = {x:eflux_time, y:eflux_tail}, dlim = {psym:-1,ytitle:'eflux!Ctail'},lim={ylog:1}
                eflux_earth=eflux_old_earth*(en_peak_earth/en_peak_earth)
                index=where( ~finite(eflux_old_tail) and en_peak_tail gt 0,ct )
                if ct gt 0 then eflux_tail(index)=eflux_old_earth(index)
                store_data, eflux_peak_earth_name, data = {x:eflux_time, y:eflux_earth}, dlim = {psym:-1,ytitle:'eflux!Cearth'},lim={ylog:1}
                if keyword_set(plot_add_eflux_procedure) then begin 
                    options,'*','panel_size',1
                    spawn, 'mkdir '+path+'plots/'+'eflux_plots'
                    fln = path+'plots/eflux_plots/o_beam_eflux_' +  date_s + '.ps' 
                    popen, fln, /port 
                    tplot,['EN*beam','EN*EFLUX*OLD','EFLUX*']
                    pclose 
                endif 
            ENDIF  
        endif  
 
;--------------------------
;Adding distribution function information
;1. plot energy spectra in unit dist func 
;2. extract the dist func correcponding to the energy peak
;-------------------------- 
        IF KEYWORD_SET(add_distfunc) THEN BEGIN
            tplot_names,'*DISTFUNC*',names=names
            store_data, delete=names
            tplot_names, 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_*_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str +'_epcut_beam_OLD', names = names
            IF names(0) NE  '' then begin
                FOR id = 0, 1 DO BEGIN
                    phi_str = phi_str_set(id)
                    en_peak_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str +'_epcut_beam'
                    options,en_peak_name+'','ytitle', 'en peak!C'+phi_str
                    enspec_distfunc_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str  +'_UNDISTFUNC_SP'+sp_str+'_ET0_All' 
                    distfunc_peak_name =  'DISTFUNC_SC' + sc_str + '_IN'+inst_str+'_' + phi_str +'_UNDISTFUNC_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'
                    
                    sat = [sc]
                    specie = [sp]
                    angle = angle_set(*, *, id)
                    inst = inst_input & units_name = 'DIST FUNC' 
                    eff_table = 0
                    
                    plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 1
                    zlim, enspec_distfunc_name, 1e-11, 1e-4,1
                    average_tplot_variable, enspec_distfunc_name, average_time, /new
                    
                    get_data, enspec_distfunc_name+'_AVG'+at_str, data = data, dlim = dlim, lim = lim
                    en_distfunc_time = data.x
                    en_distfunc_distfunc = data.y
                    en_distfunc_en = data.v

                    get_data, en_peak_name+'_OLD', data = data
                    en_peak_old_time = data.x
                    en_peak_old = data.y
                    
                    IF total(en_distfunc_time NE en_peak_old_time) EQ 0 THEN BEGIN 
                        FOR iavg = 0, n_avg-1 DO  en_distfunc_distfunc(iavg, *) = en_distfunc_distfunc(iavg, *)*(en_distfunc_en(iavg, *) EQ en_peak_old(iavg)) 
                        en_distfunc_cut = total(en_distfunc_distfunc, 2,/nan)
                        index= where(en_distfunc_cut eq 0,ct) & if ct gt 0 then en_distfunc_cut(index) = !VALUES.F_NAN
                        store_data, enspec_distfunc_name+'_AVG'+at_str+'_peak_OLD', data = {x:en_distfunc_time, y:en_distfunc_distfunc, v:en_distfunc_en}, dlim = dlim, lim = lim
                        store_data, distfunc_peak_name+'_OLD', data = {x:en_distfunc_time, y:en_distfunc_cut}, dlim = {psym:-1,ytitle:'distfunc old!C'+phi_str},lim={ylog:1}
                        
                        options,enspec_distfunc_name+'_AVG'+at_str+'_peak_OLD','ytitle','energy(eV) old!C'+phi_str
                    ENDIF ELSE stop
                    
                ENDFOR
                
                en_peak_tail_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str_set(0) +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str +'_epcut_beam'
                en_peak_earth_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str_set(1) +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str +'_epcut_beam'              
                distfunc_peak_tail_name =  'DISTFUNC_SC' + sc_str + '_IN'+inst_str+'_' + phi_str_set(0) +'_UNDISTFUNC_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'
                distfunc_peak_earth_name =  'DISTFUNC_SC' + sc_str + '_IN'+inst_str+'_' + phi_str_set(1) +'_UNDISTFUNC_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'

                get_data,distfunc_peak_tail_name+'_OLD',data=data
                distfunc_time=data.x
                distfunc_old_tail = data.y
     
                get_data,distfunc_peak_earth_name+'_OLD',data=data
                distfunc_old_earth = data.y
               
                get_data,en_peak_tail_name,data=data
                en_peak_tail=data.y
                get_data,en_peak_earth_name,data=data
                en_peak_earth=data.y
                
                distfunc_tail=distfunc_old_tail*(en_peak_tail/en_peak_tail)
                index=where( ~finite(distfunc_old_tail) and en_peak_tail gt 0,ct )
                if ct gt 0 then distfunc_tail(index)=distfunc_old_earth(index)
                store_data, distfunc_peak_tail_name, data = {x:distfunc_time, y:distfunc_tail}, dlim = {psym:-1,ytitle:'distfunc!Ctail'},lim={ylog:1}
                distfunc_earth=distfunc_old_earth*(en_peak_earth/en_peak_earth)
                index=where( ~finite(distfunc_old_tail) and en_peak_tail gt 0,ct )
                if ct gt 0 then distfunc_tail(index)=distfunc_old_earth(index)
                store_data, distfunc_peak_earth_name, data = {x:distfunc_time, y:distfunc_earth}, dlim = {psym:-1,ytitle:'distfunc!Cearth'},lim={ylog:1}
                if keyword_set(plot_add_distfunc_procedure) then begin 
                    options,'*','panel_size',1
                    spawn, 'mkdir '+path+'plots/'+'distfunc_plots'
                    fln = path+'plots/distfunc_plots/o_beam_distfunc_' +  date_s + '.ps' 
                    zlim, 'EN*beam', 1e-11,1e-4,1
                    ylim,'EN*beam',40,40000.,1
                    ylim,'EN*DISTFUNC*OLD',40,40000.,1
                    ylim,'DISTFUNC*',1e-11,1e-4,1
                    popen, fln, /port 
                    tplot,['EN*beam','EN*DISTFUNC*OLD','DISTFUNC*']
                    pclose 
                endif 
            ENDIF   
        endif    
      
;-----------------
;Adding anodes information
;1. load the energy peak and pitch angle peak from file
;2. plot globe plots for all time
;-----------------
        IF keyword_set(add_anodes) THEN BEGIN 
                                ;          FOR id = 0, 1 DO BEGIN 
                                ;             theta_avg_name =      'THETA_ENVARIOUS'+ '_SC' + sc_str+'_' +phi_str_set(id)+'_SP'+sp_str+'_ET0_All_'+'_AVG'+at_str
                                ;            theta_avg_name_new =  'THETA_ENVARIOUS'+ '_SC' + sc_str+'_' +phi_str_set(id)+'_SP'+sp_str+'_ET0_All_'+'AVG'+at_str
                                ;           tplot_names, theta_avg_name, names = names
                                ;          IF names NE '' THEN BEGIN 
                                ;             get_data, theta_avg_name, data = data
                                ;            store_data,  theta_avg_name_new , data = data
                                ;           store_data, delete = theta_avg_name
                                ;      ENDIF 
                                ; ENDFOR 

            FOR id = 0, 1 DO BEGIN 
                phi_str = phi_str_set(id)
                
                pa_beam_name = 'PASPEC_SC'+sc_str+'_IN'+inst_str+'_'+phi_str+ $
                  '_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str+'_PAP_PA_beam'
                
                en_peak_name = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_'+ phi_str $
                  +'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+ at_str $
                  +'_epcut_beam'
                pos = STREGEX(pa_beam_name, 'AVG') 
                bins_name = 'GLOBE'+STRMID(pa_beam_name, 6, pos)+'_bin_peak'
                
                get_data, en_peak_name, data = data
                energy_peak_time = data.x
                energy_peak = data.y
                
                get_data, pa_beam_name, data = data
                pa_peak = REPLICATE(!VALUES.F_NAN, n_avg)
                FOR i = 0, n_avg-1 DO  BEGIN 
                    IF total(data.y(i, *), /nan) GT 0 THEN $
                      pa_peak(i) = data.v(i, where(data.y(i, *)EQ max(data.y(i, *), /nan)))
                ENDFOR
                
                ntime = floor((end_time-start_time)/average_time)
                n_avg = N_ELEMENTS(energy_peak_time)                 
                theta_avg = FLTARR(n_avg)
                theta_avg(*) = !VALUES.F_NAN
                FOR i_m = 0, ntime-1 DO BEGIN
                    time = start_time+i_m*average_time
                    timespan, time, average_time, /SECONDS
                    index = where(energy_peak_time GT  time AND energy_peak_time LE time+average_time, ct)
                    loc = index(0)
                    IF ct EQ 1 THEN BEGIN 
                        IF keyword_set(use_angle_range) THEN $
                          pitch_angle_range = [pa_peak(loc(0))-22.5, pa_peak(loc(0))+22.5] $
                        ELSE pitch_angle_range = [0, 180]
                        IF KEYWORD_SET(use_energy_range) THEN $ 
                          energy = [FLOOR(energy_peak(loc)), CEIL(energy_peak(loc))] $
                        ELSE energy = [40., 40000.]
                        
                        index = where(pitch_angle_range GT 0, ct1)
                        index = where(energy_peak(loc) GE 0, ct2)
                        IF ct1*ct2 EQ 0 AND ct1+ct2 GT 0 THEN stop
                        IF  ct1 GT 0 AND ct2 GT 0 THEN BEGIN 
                            sat = sc
                            specie = sp
                            units_name = 'DIFF FLUX' & inst = inst_input &  eff_table = 0
                            prod = [17, 18, 47, 49]
                            COMMON get_error, get_err_no, get_err_msg, default_verbose
                            get_err_no = 0 ; reset error indicator
                            found = 0 ;reset prod found indicator
                            gname = strarr(4)
                            
                            FOR jj = 0, n_elements(prod)-1 DO BEGIN
                                dat = call_function('get_cis_cod_data', prod(jj), specie = specie, sat)
                                IF get_err_no EQ 0 THEN BEGIN ;check if data were found for time interval
                                    mag_theta = 0. &  mag_phi = 180.
                                    get_theta_phi, sat, dat, mag_theta, mag_phi, inst

                                    gname(found) = 'GLOBE_SC'+strcompress(sat, /remove_all)+'_'+$
                                      strcompress(units_name, /remove_all)+'_PR'+$
                                      strcompress(string(prod(jj)), /remove_all)+'_SP' +$
                                      strcompress(specie, /remove_all)
                                    codif_globe, sat, inst, prod(jj), eff_table, $
                                      SPECIE = SPECIE, UNITS = units_name, $
                                      NAME = gname(found), $
                                      MAG_THETA = mag_theta, MAG_PHI = mag_phi
                                    found = found+1
                                ENDIF
                            ENDFOR
                            
                            IF found GT 0 THEN BEGIN 
                                max_flux = 0
                                FOR ip = 0, found-1 DO BEGIN
                                ;load the data from globe plot
                                    get_data, gname(ip), data = dat
                                ;set up ebins
                                    ebins = replicate(0, dat.nenergy)
                                    er = [energy_to_ebin(dat, energy)]
                                    IF er(0) GT er(1) THEN er = reverse(er)
                                    n_erange = ROUND(dat.nenergy/15.)
                                    ebins(((er(0)-n_erange) > 0):((er(1)+n_erange) < (dat.nenergy-2))) = 1
                                    ebins = ebins#replicate(1, 88)
                                ;set up bins
                                    pitch_angle =  get_pitch_angle(dat, dat.mag_theta, dat.mag_phi, combine = 1)
                                    pa_ind = pitch_angle GE  min(pitch_angle_range) AND $
                                      pitch_angle LT  max(pitch_angle_range)                            

                                    str_element, dat, 'data', dat.data*(ebins*pa_ind), add_replace = 1
                                    
                                    IF max(dat.data) GT max_flux THEN BEGIN 
                                        max_flux = max(dat.data)
                                        index = where(dat.data EQ max(dat.data), ct)
                                        IF ct EQ 1 THEN theta_avg(loc) = dat.theta(index(0)) ELSE BEGIN 
                                            theta_avg(loc) = dat.theta(index(0))
                                            FOR i = 1, ct-1 DO IF theta_avg(loc) NE dat.theta(index(i)) THEN stop
                                        ENDELSE 
                                    ENDIF 
                                ENDFOR  
                                
                                FOR ip = 0, found-1 DO BEGIN 
                                    tplot_names, gname(ip)+'*', names = names
                                    store_data, delete = names
                                    tplot_names, 'B_xyz*', names = names
                                    store_data, delete_names
                                ENDFOR 
                            ENDIF ELSE stop      
                        ENDIF 
                    ENDIF  
                ENDFOR    
                timespan, t_s, t_dt, /SECONDS  
                theta_avg_name = 'THETA_ENVARIOUS'+ '_SC' + sc_str+'_' $
                  +phi_str+'_SP'+sp_str+'_ET0_All_'+'AVG'+at_str
                store_data,  theta_avg_name, data = {x:time_avg,  y:theta_avg}
            ENDFOR    
        ENDIF       
;gsm_change=1
                                ;  if keyword_set(gsm_change) then begin 
  ;      sat=sc
 ;       get_cluster_ephemeris, sat, /GSM_X, /GSM_Y, /GSM_Z
                                ;      for igsm=0, 2 do begin 
                                ;         if igsm eq 0 then  gsm_co= 'EPH_SC'+ sc_str+'_GSM_X'
                                ;        if igsm eq 1 then   gsm_co  = 'EPH_SC'+ sc_str+'_GSM_Y'
                                ;       if igsm eq 2 then  gsm_co   = 'EPH_SC'+ sc_str+'_GSM_Z'
        
                                ;      get_data,gsm_co,data=data,dlim=dlim
;stop
                                ;     store_data,gsm_co,data={x:data.x,y:data.y/6370.},dlim=dlim
                                ;     endfor 
                                ;endif 
;stop
;--------------------------------------------------------------
;Overview plots
;--------------------------------------------------------------
;reset plots options
        p01 = 'TDMOM_EN00040_40000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_P_total'
        p02 = beta_name
        p03 = 'EFLUX_SC' + sc_str + '_IN'+inst_str+'_' + phi_str_set(0) +'_UNEFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'
        p03DF = 'DISTFUNC_SC' + sc_str + '_IN'+inst_str+'_' + phi_str_set(0) +'_UNDISTFUNC_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'
        p04 = 'ENSPEC_SC' + sc_str + '_IN'+inst_str+'_' + phi_str_set(0)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All'
        p05 = 'EFLUX_SC' + sc_str + '_IN'+inst_str+'_' + phi_str_set(1) +'_UNEFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'
        p05DF = 'DISTFUNC_SC' + sc_str + '_IN'+inst_str+'_' + phi_str_set(1) +'_UNDISTFUNC_SP'+sp_str+'_ET0_All_AVG'+at_str+'_peak'
        p06 = 'ENSPEC_SC'+ sc_str +'_IN'+inst_str+'_'+phi_str_set(1)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All'
;     p07 = 'Dst_Index'               
        p08 = 'MAG_SC'+sc_str+'_B_xyz_gse_X'
        p09 = 'ENSPEC_SC'+sc_str+'_IN'+inst_str+'_'+phi_str_set(0)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str
        p10 = 'ENSPEC_SC'+sc_str+'_IN'+inst_str+'_'+phi_str_set(1)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str
        p11 = 'PASPEC_SC'+sc_str+'_IN'+inst_str+'_'+phi_str_set(0)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str
        p12 = 'PASPEC_SC'+sc_str+'_IN'+inst_str+'_'+phi_str_set(1)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str
        
        P111 = 'PASPEC_SC'+sc_str+'_IN'+inst_str+'_'+phi_str_set(0)+'_UNCOUNTS_SP'+sp_str+'_ET0_All_AVG'+at_str
        P121 = 'PASPEC_SC'+sc_str+'_IN'+inst_str+'_'+phi_str_set(1)+'_UNCOUNTS_SP'+sp_str+'_ET0_All_AVG'+at_str

        p13 = p11+'_PAP'
        p14 = P12+'_PAP'
        p15 = p13 +'_ET'
        p16 = p14 +'_ET'
        p17 = p15+'_beam'
        p18 = p16+'_beam'
        p19 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(0)+'_MTDENSITY_SP'+sp_str+'_ET0_All' $
          + '_AVG'+at_str
        p20 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(1)+'_MTDENSITY_SP'+sp_str+'_ET0_All' $
          +'_AVG'+at_str
        p21 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(0)+'_MTVELOCITY_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str
        p22 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(1)+'_MTVELOCITY_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str
        p23 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(0)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str
        p24 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(1)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str
        p25 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(0)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str+'_neg'
        p26 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(1)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str+'_neg'
        p27 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(0)+'_MTPRESSURE_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str
        p28 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(1)+'_MTPRESSURE_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str
        p29 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(0)+'_MTPRESSURE_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str+'_neg'
        p30 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(1)+'_MTPRESSURE_SP'+sp_str+'_ET0_All_T' $
          +'_AVG'+at_str+'_neg'
        p31 = 'EPH_SC'+ sc_str+'_GSE_X'
        p32 = 'EPH_SC'+ sc_str+'_GSE_Y'
        p33 = 'EPH_SC'+ sc_str+'_GSE_Z'
        p34 = 'PASPEC_SC'+sc_str+'_IN'+inst_str+'_PHICOMBINED_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str+'_PAP_ET_beam'
        p35 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(0)+'_MTVELOCITY_SP'+sp_str+'_ET0_All' $
          +'_AVG'+at_str+'_V_PAR_T'
        p36 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(1)+'_MTVELOCITY_SP'+sp_str+'_ET0_All' $
          +'_AVG'+at_str+'_V_PAR_T'
        p37 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(0)+'_MTVELOCITY_SP'+sp_str+'_ET0_All' $
          +'_AVG'+at_str+'_V_PERP_T'
        p38 = 'TDMOM_ENVARIOUS_SC'+ sc_str +'_'+phi_str_set(1)+'_MTVELOCITY_SP'+sp_str+'_ET0_All' $
          +'_AVG'+at_str+'_V_PERP_T'
        p39 = 'EPH_SC'+ sc_str+'_GSM_X'
        p40 = 'EPH_SC'+ sc_str+'_GSM_Y'
        p41 = 'EPH_SC'+ sc_str+'_GSM_Z'
        p42 = 'MAG_SC'+sc_str+'_B_xyz_gse'            
        p43 = 'TDMOM_EN00040_40000_SC'+sc_str+'_MTDENSITY_SP0_ET0_All'
        p44 = 'TDMOM_EN00040_40000_SC'+sc_str+'_MTVELOCITY_SP0_ET0_All'
        p45 = 'TDMOM_EN00040_40000_SC'+sc_str+'_MTTEMPERATURE_SP0_ET0_All'
        p46 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
          +''+phi_str_set(0)+''+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_X'+'_AVG'+at_str
        p47 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
          +''+phi_str_set(0)+''+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Y'+'_AVG'+at_str
        p48 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
          +''+phi_str_set(0)+''+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Z'+'_AVG'+at_str
        p49 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
          +''+phi_str_set(1)+''+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_X'+'_AVG'+at_str
        p50 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
          +''+phi_str_set(1)+''+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Y'+'_AVG'+at_str
        p51 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
          +''+phi_str_set(1)+''+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Z'+'_AVG'+at_str
        p52 = 'storm_phase'     
        p53 = 'ACE_MAG_GSM_X'
        p54 = 'ACE_MAG_GSM_Y'
        p55 = 'ACE_MAG_GSM_Z'
        p56 = 'ACE_SWEPAM_VELOCITY'
        p57 = 'ACE_SWEPAM_PRESSURE'
        p58 = 'THETA_ENVARIOUS'+ '_SC' + sc_str+'_' +phi_str_set(0)+'_SP'+sp_str+'_ET0_All_'+'AVG'+at_str
        p59 = 'THETA_ENVARIOUS'+ '_SC' + sc_str+'_' +phi_str_set(1)+'_SP'+sp_str+'_ET0_All_'+'AVG'+at_str
        p60 = 'EPH_SC'+sc_str+'_MLT'
        p61 = 'EPH_SC'+sc_str+'_ILAT_D'
        p62 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' +''+phi_str_set(0)+''+'_MTVELOCITY_SP'+sp_str+'_ET0_All'+'_AVG'+at_str
        p63 = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' +''+phi_str_set(1)+''+'_MTVELOCITY_SP'+sp_str+'_ET0_All'+'_AVG'+at_str

        options, '*', 'panel_size', 1
        options, '*', 'zticks', 3
        options, [p06, p09, p10, p11, p12, p13, p14, p17, p18, p34], 'ztitle', ''

        ylim, p01, 0.01, 3, 1
        ylim, p02, 0.01, 10, 1
        zlim, [p04, p06], 0.1, 100, 1
        zlim, [p11, p12], 0, 100, 0

        IF inst_input EQ 0 THEN BEGIN 
            options, p04, 'ytitle', 'SC' + sc_str + ' O!U+!N!C!C(eV)' + '!C!C' + 'Tailward'
            options, p06, 'ytitle', 'SC' + sc_str + ' O!U+!N!C!C(eV)' + '!C!C' + 'Earthward'
            options, p09, 'ytitle', 'SC' + sc_str + ' O!U+!N (eV)!C!CTailward!C!CAVG-'+at_str
            options, p10, 'ytitle', 'SC' + sc_str + ' O!U+!N (eV)!C!CEarthward!C!CAVG-'+at_str
        ENDIF ELSE BEGIN 
            IF inst_input EQ 1 THEN BEGIN 
                options, p04, 'ytitle', 'SC' + sc_str + ' HIA!C!C(eV)' + '!C!C' + 'Tailward'
                options, p06, 'ytitle', 'SC' + sc_str + ' HIA!C!C(eV)' + '!C!C' + 'Earthward'
                options, p09, 'ytitle', 'SC' + sc_str + ' HIA (eV)!C!CTailward!C!CAVG-'+at_str
                options, p10, 'ytitle', 'SC' + sc_str + ' HIA (eV)!C!CEarthward!C!CAVG-'+at_str
            ENDIF 
        ENDELSE 
        options, p34, 'ytitle', 'SC'+sc_str+' O!U+!C!CBEAM!C!CE-----T'
        options, [p09+'_erange', p10+'_erange'], 'color', 2            
        options, p08, 'ytitle', 'SC' + sc_str + '!C!CBx (nT)'             
        options,  p11, 'ytitle', 'Pitch Angle!C!CVarious EN'
        options,  p12, 'ytitle', 'Pitch Angle!C!CVarious EN'                  
        options, [p13, p14], 'ytitle', 'Pitch Angle!C!CLocal Peak'
        var_label = 'EPH_SC' + sc_str + '_'
        var_label = var_label + ['MLT', 'GSM_X', 'GSM_Y', 'GSM_Z', 'DIST']
; set mom data options if plot_mom is required
        IF KEYWORD_SET(plot_mom)  THEN BEGIN 
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
            options, [p37, p38], 'ytitle', 'SC'+sc_str+' O!N+!N!C!CV!LL!N (kms!U-1!N)!C'
            options, [p23, p24, p25, p26], 'ytitle', 'SC' + sc_str +' O!U+!N!C!CT!Davg!N (eV)'
            options, [p27, p28, p29, p30], 'ytitle', 'SC'+sc_str+' O!U+!N!C!cP!Davg!N (nPa)'   
            options, [p19, p20], 'ytitle', 'SC'+sc_str+' O!U+!N!C!Cn (cm!U-3!N)'
            options, [p21, p22], 'ytitle', 'SC'+sc_str+' O!U+!N!C!CV!LT!N (kms!U-1!N)!C' 
            options, [p35, p36], 'ytitle', 'SC'+sc_str+' O!N+!N!C!CV!L//!N (kms!U-1!N)!C'
            options, [p37, p38], 'ytitle', 'SC'+sc_str+' O!N+!N!C!CV!LL!N (kms!U-1!N)!C'
            options, [p23, p24, p25, p26], 'ytitle', 'SC' + sc_str +' O!U+!N!C!CT!Davg!N (eV)'
            options, [p27, p28, p29, p30], 'ytitle', 'SC'+sc_str+' O!U+!N!C!cP!Davg!N (nPa)'  
        ENDIF 

        IF NOT KEYWORD_SET(displaytime) THEN displaytime = t_dt
        get_data, 'location', data = data_m
        index = where(data_m.y(*, 2) EQ 1, ct_magnetosphere)
        IF ct_magnetosphere ne 0 or KEYWORD_SET(plot_sw_sheath) THEN BEGIN  
            FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
                ts_plot = time_string(t_s+idisplay*displaytime)
                te_plot = time_string(t_s+(idisplay+1)*displaytime)
                date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
                time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
                date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
                time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
                year = STRMID(ts_plot, 0, 4)
                timespan, t_s+idisplay*displaytime, displaytime, /SECONDS
;plot in idl windows
                IF KEYWORD_SET(idl_plot) THEN BEGIN 
                  ;  window, idisplay
                    tplot, [p02, p04, p06, p34, p11, p12,  p13, p14], var_label = var_label
                    tplot_panel, v = p04, o = p09+'_epcut_beam', psym = 0;,thick=2
                    tplot_panel, v = p06, o = p10+'_epcut_beam', psym = 0;,thick=2
                    tplot_panel, v = p04, o = p09+'_erange', psym = 0
                    tplot_panel, v = p06, o = p10+'_erange', psym = 0
                    
                    yline, p02, offset = 0.05, col = 1
                    yline, p02, offset = 1, col = 1
                    stop
                ENDIF 
;---------------------------------------
; Plot the graph in PS file if ps is set to 1 
;---------------------------------------
                IF KEYWORD_SET(ps) THEN BEGIN  
                    spawn, 'mkdir '+path+'plots'
                    spawn, 'mkdir '+path+'plots/'+'obeam_3pages'
                    spawn, 'mkdir '+path+'plots/'+'obeam_3pages/'+year
                    spawn, 'mkdir '+path+'plots/'+'obeam_3pages/'+year+'/page1_general_result'
                    spawn, 'mkdir '+path+'plots/'+'obeam_3pages/'+year+'/page2_mom'
                    spawn, 'mkdir '+path+'plots/'+'obeam_3pages/'+year+'/page3_identification_process'
                    
                    IF keyword_set(plot_mom) THEN BEGIN 
                        get_data, 'location', data = data
                        magnetosphere_data_index = where(data.y NE 0)

                        get_data, p23, data = data
                        index = where(data.y(magnetosphere_data_index) LT 0, ct1)
                        get_data, p24, data = data
                        index = where(data.y(magnetosphere_data_index) LT 0, ct2)
                        get_data, p27, data = data
                        index = where(data.y(magnetosphere_data_index) LT 0, ct3)
                        get_data, p28, data = data
                        index = where(data.y(magnetosphere_data_index) LT 0, ct4)
                        IF (ct1+ct2+ct3+ct4) GT 0 THEN BEGIN 
                            spawn, 'mkdir '+path+'/plots/'+'mom_with_neg'
                            spawn, 'mkdir '+path+'/plots/'+'mom_with_neg/'+year

                            get_data, p23, data = data, dlim = dlim, lim = lim
                            store_data, p25, data = data, dlim = dlim, lim = lim
                            low_limit_t = min(data.y)*10
                            get_data, p24, data = data, dlim = dlim, lim = lim
                            store_data, p26, data = data, dlim = dlim, lim = lim
                            low_limit_t = low_limit_t < min(data.y)*10
                            get_data, p27, data = data, dlim = dlim, lim = lim
                            store_data, p29, data = data, dlim = dlim, lim = lim
                            low_limit_p = min(data.y)*10
                            get_data, p28, data = data, dlim = dlim, lim = lim
                            store_data, p30, data = data, dlim = dlim, lim = lim
                            low_limit_p = low_limit_p < min(data.y)*10
                            
                            options, [p25, p26], 'ytitle',  $
                              'SC' + sc_str +'!C!CO!U+!N T!Davg!N (eV)'  
                            options, [p29, p30], 'ytitle', $ 
                              'SC'+sc_str+'!C!CO!U+!N P!Davg!N (nPa)' 
                            options, [p25, p26, p29, p30], 'psym', 7
                            
                            ylim, [p25, p26], low_limit_t, 0, 0
                            ylim, [p29, p30], low_limit_p, 0, 0
                            
                            spawn, 'mkdir '+path+'plots/mom_with_neg/'
                            fln = path+'plots/mom_with_neg/' + year $
                              +'/storm_o_beam_mom_neg_' + $
                              date_s_plot + '_' + time_s_plot + '_to_'+ $
                              date_e_plot + '_' + time_e_plot + 'sc'+ $
                              sc_str+'.ps' 

                            popen, fln, /port               
                            tplot, ['location', p09, p10, p34, p19, $
                                    p21, p23, p25, p27, p29], $
                              var_label = var_label
                            tplot_panel, v = p09, o = p09+'_epcut_beam', psym = 0
                            tplot_panel, v = p10, o = p10+'_epcut_beam', psym = 0
                            tplot_panel, v = p09, o = p09+'_erange', psym = 0
                            tplot_panel, v = p10, o = p10+'_erange', psym = 0
                            tplot_panel, v = p19, o = p20, psym = 1
                            tplot_panel, v = p21, o = p22, psym = 1
                            tplot_panel, v = p23, o = p24, psym = 1
                            tplot_panel, v = p25, o = p26, psym = 7 
                            tplot_panel, v = p27, o = p28, psym = 1
                            tplot_panel, v = p29, o = p30, psym = 7
                            yline, p02, offset = 0.05, col = 1
                            yline, p02, offset = 1, col = 1                    
                            pclose
                            
                            IF idisplay EQ 0 THEN BEGIN 
                                OPENU, unit, path+$
                                  'log_neg_T.txt', /GET_LUN, /APPEND
                                PRINTF, unit,  ts+' - '+te, $
                                  '-----negative Tt and Pt-----sc', sc_str
                                FREE_LUN, unit 
                            ENDIF     
                        ENDIF 

                        fln = path+'plots/obeam_3pages/'+year+'/page1_general_result/storm_o_beam_mom' +$
                          date_s_plot + '_' + time_s_plot + '_to_'+ $
                          date_e_plot + '_' + time_e_plot + '_page1.ps' 
                        popen, fln, /port     
                        tplot, [p52, 'location', p02, p09, p10, p34, p19, p21, p23], $
                          var_label = var_label
                        tplot_panel, v = p09, o = p09+'_epcut_beam', psym = 0
                        tplot_panel, v = p10, o = p10+'_epcut_beam', psym = 0
                        tplot_panel, v = p09, o = p09+'_erange', psym = 0
                        tplot_panel, v = p10, o = p10+'_erange', psym = 0
                        tplot_panel, v = p19, o = p20, psym = 1
                        tplot_panel, v = p21, o = p22, psym = 1
                        tplot_panel, v = p23, o = p24, psym = 1 
                        yline, p02, offset = 0.05, col = 1
                        yline, p02, offset = 1, col = 1       
                        pclose

                        fln = path+'plots/obeam_3pages/'+year+'/page2_mom/storm_o_beam_mom' +$
                          date_s_plot + '_' + time_s_plot + '_to_'+ $
                          date_e_plot + '_' + time_e_plot + '_page2.ps' 
                        popen, fln, /port
                        tplot, ['location', p19, p21, p35, p37, p23, p27, P46, P47, P48], $
                          var_label = var_label
                        tplot_panel, v = p19, o = p20, psym = 1
                        tplot_panel, v = p21, o = p22, psym = 1
                        tplot_panel, v = p23, o = p24, psym = 1 
                        tplot_panel, v = p27, o = p28, psym = 1
                        tplot_panel, v = p35, o = p36, psym = 1
                        tplot_panel, v = p37, o = p38, psym = 1
                        tplot_panel, v = p46, o = p49, psym = 1
                        tplot_panel, v = p47, o = p50, psym = 1
                        tplot_panel, v = p48, o = p51, psym = 1
                        pclose
                    ENDIF          
                    fln = path+'plots/obeam_3pages/'+year $
                      +'/page3_identification_process/storm_o_beam_mom'+ $
                      date_s_plot + '_' + time_s_plot + '_to_'+ $
                      date_e_plot + '_' + time_e_plot + '_page3.ps' 
                    popen, fln, /port
                    
                    tplot, [p02, p04, p06, p34, p09, p11, p17+'_OLD', $
                            p10, p12, p18+'_OLD'], $
                      var_label = var_label
                    tplot_panel, v = p09, o = p09+'_epcut', psym = -7
                    tplot_panel, v = p10, o = p10+'_epcut', psym = -7
                    tplot_panel, v = p04, o = p09+'_epcut_beam'
                    tplot_panel, v = p06, o = p10+'_epcut_beam'
                    yline, p08, col = 3
                    pclose
                ENDIF  

; fit t
                IF keyword_set(dfit_temperature) THEN BEGIN 
                    t_name_tail = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(0)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T_AVG'+at_str
                    t_name_tail_para = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(0)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_X_AVG'+at_str
                    t_name_tail_perp1 = 'TDMOM_ENVARIOUS_SC' +sc_str $
                      +'_'+phi_str_set(0)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Y_AVG'+at_str
                    t_name_tail_perp2 = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(0)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Z_AVG'+at_str
                    t_fit_name_tail_para = 'TDMOM_ENVARIOUS_SC' + sc_str $
                      +'_'+phi_str_set(0)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_para_dfit_AVG'+at_str
                    t_fit_name_tail_perp = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(0)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_perp_dfit_AVG'+at_str
                    t_fit_name_tail = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(0)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T_dfit_AVG'+at_str
                    
                    t_name_earth = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(1)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T_AVG'+at_str
                    t_name_earth_para = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(1)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_X_AVG'+at_str
                    t_name_earth_perp1 = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(1)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Y_AVG'+at_str
                    t_name_earth_perp2 = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(1)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_Z_AVG'+at_str
                    t_fit_name_earth_para = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(1)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_para_dfit_AVG'+at_str
                    t_fit_name_earth_perp = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(1)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_perp_dfit_AVG'+at_str
                    t_fit_name_earth = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(1)+'_MTTEMPERATURE_SP'+sp_str+'_ET0_All_T_dfit_AVG'+at_str
                    d_name_tail = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(0)+'_MTDENSITY_SP'+sp_str+'_ET0_All_AVG'+at_str
                    
                    d_fit_name_tail = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(0)+'_MTDENSITY_SP'+sp_str+'_ET0_All_dfit_AVG'+at_str
                    d_name_earth = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(1)+'_MTDENSITY_SP'+sp_str+'_ET0_All_AVG'+at_str
                    d_fit_name_earth = 'TDMOM_ENVARIOUS_SC'+sc_str $
                      +'_'+phi_str_set(1)+'_MTDENSITY_SP'+sp_str+'_ET0_All_dfit_AVG'+at_str
                    
                    options, d_fit_name_earth, 'color'
                    ylim, [d_fit_name_tail, d_fit_name_earth], 0.001, 1
                    options, [ d_name_tail, d_name_earth, $ 
                               t_name_tail, t_name_tail_para, $
                               t_name_tail_perp1, t_name_tail_perp2, $
                               t_name_earth, t_name_earth_para, $
                               t_name_earth_perp1, t_name_earth_perp2], $
                      'color', 6
                    ylim, [t_fit_name_tail, t_fit_name_tail_para, $
                           t_fit_name_tail_perp, $
                           t_fit_name_earth, t_fit_name_earth_para, $
                           t_fit_name_earth_perp], 0.1, 10000
                    
                    IF NOT keyword_set(show_fit) THEN $
                      popen, path+'temperature_compare'+ date_s_plot + '_' + time_s_plot + '_to_'+ $
                      date_e_plot + '_' + time_e_plot+'_tail.ps', /land ELSE window, /free
                    tplot, [p09, d_fit_name_tail, t_fit_name_tail, $
                            t_fit_name_tail_para, t_fit_name_tail_perp], title = path

                    tplot_panel, v = d_fit_name_tail, o = d_name_tail, psym = 7
                    tplot_panel, v = t_fit_name_tail, o = t_name_tail, psym = 7
                    tplot_panel, v = t_fit_name_tail_para, o = t_name_tail_para, psym = 7
                    tplot_panel, v = t_fit_name_tail_perp, o = t_name_tail_perp1, psym = 7
                    tplot_panel, v = t_fit_name_tail_perp, o = t_name_tail_perp2, psym = 7
                    IF NOT keyword_set(show_fit) THEN  pclose

                    IF NOT keyword_set(show_fit) THEN $
                      popen, path+'temperature_compare'+ date_s_plot + '_' + time_s_plot + '_to_'+ $
                      date_e_plot + '_' + time_e_plot+'_earth.ps', /land ELSE window, /free
                    tplot, [p10, d_fit_name_earth, t_fit_name_earth, $
                            t_fit_name_earth_para, t_fit_name_earth_perp], title = path
                    
                    tplot_panel, v = t_fit_name_earth, o = t_name_earth, psym = 7
                    tplot_panel, v = t_fit_name_earth_para, o = t_name_earth_para, psym = 7
                    tplot_panel, v = t_fit_name_earth_perp, o = t_name_earth_perp1, psym = 7
                    tplot_panel, v = t_fit_name_earth_perp, o = t_name_earth_perp2, psym = 7

                    tplot_panel, v = d_fit_name_earth, o = d_name_earth, psym = 7
                    
                    xyouts, 7000, 50, 'T (fit)', size = 2
                    xyouts, 7000, 1000, 'T (3dmom)', color = 6, size = 2
                    IF NOT keyword_set(show_fit) THEN  pclose
                ENDIF  
            ENDFOR       
            timespan, t_s, t_dt, /SECONDS
        ENDIF                     
;------------- --------------
;restore useful tplot
;---------------------------
        IF KEYWORD_SET(store_data)  THEN  BEGIN 
            tplot_names, 'TDMOM_EN00040_40000*SP'+sp_str+'*', names = names
            store_data, delete = names
            tplot_names, 'TDMOM*SP0*X', names = names
            store_data, delete = names
            tplot_names, 'TDMOM*SP0*Y', names = names
            store_data, delete = names
            tplot_names, 'TDMOM*SP0*Z', names = names
            store_data, delete = names
            tplot_names, 'TDMOM*SP0*H', names = names
            store_data, delete = names
            tplot_names, 'TDMOM*SP0*T', names = names
            store_data, delete = names
            tplot_names, 'TDMOM*FLUX*', names = names
            store_data, delete = names
            tplot_names, '*gse_GXM*', names = names
            store_data, delete = names
            store_data, delete = 'gse_gsm__CL_SP_AUX'
            tplot_names, '*VcrossB*', names = names
            store_data, delete = names      
;            tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP'+sp_str+'_ET0_All_AVG' +at_str, names = names
;            store_data, delete = names
            tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP'+sp_str+'_ET0_All_AVG' +at_str+'_T', names = names
            store_data, delete = names
            tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP'+sp_str+'_ET0_All_AVG' +at_str+'*H', names = names
            store_data, delete = names
            tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP'+sp_str+'_ET0_All_AVG'+at_str+'_V_PERP1', names = names
            store_data, delete = names
            tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP'+sp_str+'_ET0_All_AVG' +at_str+'_V_PERP2', names = names
            store_data, delete = names
            tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP'+sp_str+'_ET0_All_AVG' +at_str+'*X', names = names
            store_data, delete = names
            tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP'+sp_str+'_ET0_All_AVG' +at_str+'*Y', names = names
            store_data, delete = names
            tplot_names, 'TDMOM_ENVARIOUS_SC'+sc_str+'*MTVELOCITY_SP'+sp_str+'_ET0_All_AVG' +at_str+'*Z', names = names
            store_data, delete = names
            tplot_names, 'ACE_SWEPAM_V*GSE', names = names
            store_data, delete = names
            store_data, delete = ['ACE_MAG_GSM_T', 'ACE_MAG_GSM_PR', $
                                  'ACE_SWEPAM_DENSITY', 'ACE_SWEPAM_TEMPERATURE', $
                                  'ACE_SWEPAM_THERM_PRESS', 'ACE_SWEPAM_THERM_PRESS', $
                                  'ACE_SWEPAM_X_GSM']
            flndata = path+'tplot_restore/o_beam_'+date_s+'_'+time_s
            tplot_save, filename = flndata
     ;       spawn,'gzip -9 '+flndata+'.tplot'
        ENDIF      
;--------------------------------------
;dump the data out if required
;--------------------------------------
      
        IF keyword_set(dumpdata) THEN BEGIN 
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
                          '    Storm_Phase', $
                          '      IMF_Bx   ', $
                          '      IMF_By   ', $
                          '      IMF_Bz   ', $
                          '       SW_V    ', $
                          '       SW_P    ', $
                          '  eflux_tail   ', $ 
                          '  eflux_earth  ', $ 
                          '  theta_tail   ', $ 
                          '  theta_earth  ', $
                          '       MLT     ', $
                          '     ILAT_D    ', $
                          '    IMF_Bx_DC  ', $
                          '    IMF_By_DC  ', $
                          '    IMF_Bz_DC  ', $
                          '      SW_V_DC  ', $
                          '      SW_P_DC  ', $
                          ' DistFunc_tail ', $
                          ' DistFunc_earth', $
                          '  Vgse_tail_x  ', $
                          '  Vgse_tail_y  ', $
                          '  Vgse_tail_z  ', $
                          '  Vgse_earth_x ', $
                          '  Vgse_earth_y ', $
                          '  Vgse_earth_z ']

            nterm = N_ELEMENTS(title_set)
            tplot_names, '*AVG*', names = names
            get_data, names(0), data = data
            time_dd = data.x

            n_avg = N_ELEMENTS(time_dd)
            title_dd = STRARR(n_avg, nterm)
            data_dd = DBLARR(n_avg, nterm)
            FOR i_time = 0, n_avg-1 DO title_dd(i_time, *) = title_set
            index_valid = where(time_dd GE t_s AND time_dd  LE t_e, ct)
            IF ct GT 0 THEN BEGIN 
                                ;Beta
                get_data, p02, data = data
                data_y = fltarr(n_avg)
                FOR itime = 0l, n_avg-1 DO BEGIN
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
; if there is more than one valid data in the average time interval
; then average them
                    IF ct GT 0 THEN BEGIN 
                        IF  TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                          data_y(itime) = total(data.y(index), /NAN)/ct $
                        ELSE data_y(itime) =  !VALUES.F_NAN
                    ENDIF ELSE  data_y(itime) =  !VALUES.F_NAN
                ENDFOR  
                data_dd(*, 1) = data_y
                                ;GSE X
                get_data, p31, data = data  
                data_y = fltarr(n_avg)
                FOR itime = 0l, n_avg-1 DO BEGIN
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN  $
                      data_y(itime) = total(data.y(index), /NAN)/ct $
                    ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR  
                data_dd(*, 2) = data_y               
                                ;GSE Y
                get_data, p32, data = data  
                data_y = fltarr(n_avg)
                FOR itime = 0l, n_avg-1 DO BEGIN
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                      data_y(itime) = total(data.y(index), /NAN)/ct $
                    ELSE data_y(itime) = !VALUES.F_NAN
                ENDFOR  
                data_dd(*, 3) = data_y             
                                ;GSE Z
                get_data, p33, data = data
                data_y = fltarr(n_avg)
                FOR itime = 0l, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                      data_y(itime) = total(data.y(index), /NAN)/ct $
                    ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR 
                data_dd(*, 4) = data_y             
                                ;GSM X
                get_data, p39, data = data  
                data_y = fltarr(n_avg)          
                FOR itime = 0, n_avg-1 DO BEGIN
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                      data_y(itime) = (total(data.y(index), /NAN)/ct) $
                    ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR  
                data_dd(*, 23 ) = data_y              
                                ;GSM Y
                get_data, p40, data = data 
                data_y = fltarr(n_avg) 
                FOR itime = 0, n_avg-1 DO BEGIN
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                      data_y(itime) = (total(data.y(index), /NAN)/ct) $
                    ELSE data_y(itime) =   !VALUES.F_NAN
                ENDFOR  
                data_dd(*, 24 ) = data_y   
                                ;GSM Z
                get_data, p41, data = data  
                data_y = fltarr(n_avg) 
                FOR itime = 0, n_avg-1 DO BEGIN
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct) 
                    IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                      data_y(itime) = (total(data.y(index), /NAN)/ct) $
                    ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR     
                data_dd(*, 25 ) = data_y              
                                ;MLT
                get_data, p60 , data = data  
                data_y = fltarr(n_avg) 
                FOR itime = 0, n_avg-1 DO BEGIN
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct) 
                    IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                      data_y(itime) = (total(data.y(index), /NAN)/ct) $
                    ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR     
                data_dd(*, 52 ) = data_y
                                ;ILAT_D
                get_data, p61, data = data  
                data_y = fltarr(n_avg) 
                FOR itime = 0, n_avg-1 DO BEGIN
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct) 
                    IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                      data_y(itime) = (total(data.y(index), /NAN)/ct) $
                    ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR     
                data_dd(*, 53 ) = data_y               
                                ;tail energy and flag
                get_data, p09+'_epcut_beam', data = data
                data_y = data.y(index_valid)
                data_dd(*, 0) = data_dd(*, 0) + (data_y GT 0)
                data_dd(*, 5) = data_y                                           
                                ; tail pap and flux
                get_data, p17, data = data
                data_y = data.y(index_valid, *)
                data_v = data.v(index_valid, *)
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data_y(itime, *) EQ max(data_y(itime, *), /nan))
                    IF index(0) EQ -1 THEN begin 
                        data_dd(itime, 6) = !VALUES.F_NAN
                        data_dd(itime, 7) = !VALUES.F_NAN
                    endif else begin 
                        data_dd(itime, 6) = data_v(itime, index(0))
                        data_dd(itime, 7) = (data_y(itime, index(0))/10.)
                    endelse 
                ENDFOR 
                IF keyword_set(plot_mom) THEN BEGIN 
                                ;tail Density
                    get_data, p19, data = data
                    data_dd(*, 8) = data.y(index_valid) 
                                ;tail V_total
                    get_data, p21, data = data
                    data_dd(*, 9) = data.y(index_valid) 
                                ;tail V_par
                    get_data, p35, data = data
                    data_dd(*, 10) = data.y(index_valid) 
                                ;tail V_perp
                    get_data, p37, data = data
                    data_dd(*, 11) = data.y(index_valid) 
                                ;tail T_total
                    get_data, p23, data = data
                    data_y = data.y(index_valid)                    
                    data_dd(*, 12) = data_y                   
                                ;tail P_total
                    get_data, p27, data = data
                    data_y = data.y(index_valid)
                    data_dd(*, 13) = data_y
                                ; tail Vgse
                    get_data, p62, data = data
                    data_y = data.y(index_valid,0)
                    data_dd(*, 61) = data_y
                    data_y = data.y(index_valid,1)
                    data_dd(*,62 ) = data_y
                    data_y = data.y(index_valid,2)
                    data_dd(*,63 ) = data_y
                ENDIF ELSE data_dd(*, 8:13) =  !VALUES.F_NAN
                                ; earth energy and flag
                get_data, p10+'_epcut_beam', data = data
                data_y = data.y(index_valid)
                data_dd(*, 0) = data_dd(*, 0) - (data_y GT 0)*10
                data_dd(*, 14) = data_y
                                ; dealing with flag 
                index = where(data_dd(*, 0) EQ -10)
                IF index(0) GE 0 THEN data_dd(index, 0) = -1
                index = where(data_dd(*, 0) EQ -9)
                IF index(0) GE 0 THEN data_dd(index, 0) = 2                 
                                ; earth pitch angle peak and flux
                get_data, p18, data = data
                data_y = data.y(index_valid, *)
                data_v = data.v(index_valid, *)
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data_y(itime, *) EQ max(data_y(itime, *), /nan))
                    IF index(0) EQ -1 THEN begin 
                        data_dd(itime, 15) =  !VALUES.F_NAN
                        data_dd(itime, 16) =  !VALUES.F_NAN
                    endif else begin 
                        data_dd(itime, 15) = data_v(itime, index(0))
                        data_dd(itime, 16) = (data_y(itime, index(0))/10.)
                    endelse 
                ENDFOR 
                IF KEYWORD_SET(plot_mom) THEN BEGIN 
                                ;earth Density
                    get_data, p20, data = data
                    data_y = data.y(index_valid)
                    data_dd(*, 17) = data_y                                  
                                ;earth V_total
                    get_data, p22, data = data
                    data_y = data.y(index_valid)
                    data_dd(*, 18) = data_y                    
                                ;earth V_par
                    get_data, p36, data = data
                    data_y = data.y(index_valid)
                    data_dd(*, 19) = data_y                 
                                ;earth V_perp
                    get_data, p38, data = data
                    data_y = data.y(index_valid)
                    data_dd(*, 20) = data_y            
                                ;earth T_total
                    get_data, p24, data = data
                    data_y = data.y(index_valid)
                    data_dd(*, 21) = data_y
                                ;earth P_total
                    get_data, p28, data = data
                    data_y = data.y(index_valid)
                    data_dd(*, 22) = data_y
                                ; earth Vgse
                    get_data, p63, data = data
                    data_y = data.y(index_valid,0)
                    data_dd(*, 64) = data_y
                    data_y = data.y(index_valid,1)
                    data_dd(*, 65) = data_y
                    data_y = data.y(index_valid,2)
                    data_dd(*, 66) = data_y
                ENDIF ELSE data_dd(*, 17:22) =  !VALUES.F_NAN
                                ; mag X_gse
                get_data, p42, data = data
                data_y = fltarr(n_avg)
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data.y(index, 0)) GE 0) GT 0  THEN $
                          data_y(itime) = total(data.y(index, 0), /NAN)/ct $
                        ELSE data_y(itime) =  !VALUES.F_NAN
                    ENDIF  ELSE  data_y(itime) =  !VALUES.F_NAN
                ENDFOR 
                index = where(data_y LE -99999, ct)
                IF ct GT 0 THEN data_y(index) =  !VALUES.F_NAN
                data_dd(*, 26 ) = data_y               
                                ; mag Y_gse
                data_y = fltarr(n_avg)
                FOR itime = 0, n_avg-1 DO BEGIN  
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF  TOTAL(ABS(data.y(index, 1)) GE 0) GT 0  THEN $
                          data_y(itime) = total(data.y(index, 1), /NAN)/ct $
                        ELSE data_y(itime) =  !VALUES.F_NAN
                    ENDIF ELSE data_y(itime) = !VALUES.F_NAN
                ENDFOR 
                index = where(data_y LE -99999, ct)
                IF ct GT 0 THEN data_y(index) =  !VALUES.F_NAN
                data_dd(*, 27 ) = data_y          
                                ; mag Z_gse
                data_y = fltarr(n_avg)
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data.y(index, 2)) GE 0) GT 0  THEN $
                          data_y(itime) = total(data.y(index, 2), /NAN)/ct $
                        ELSE data_y(itime) =  !VALUES.F_NAN
                    ENDIF ELSE  data_y(itime) =  !VALUES.F_NAN
                ENDFOR 
                index = where(data_y LE -99999, ct)
                IF ct GT 0 THEN data_y(index) =  !VALUES.F_NAN
                data_dd(*, 28 ) = data_y      
                                ; H+ Density
                data_y = fltarr(n_avg)
                get_data, p43, data = data
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                          data_y(itime) = total(data.y(index), /NAN)/ct $
                        ELSE data_y(itime) =  !VALUES.F_NAN
                    ENDIF ELSE  data_y(itime) =  !VALUES.F_NAN
                ENDFOR 
                data_dd(*, 29 ) = data_y           
                                ;  H+ Vx
                data_y = fltarr(n_avg)
                get_data, p44, data = data
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data.y(index, 0)) GE 0) GT 0  THEN $
                          data_y(itime) = total(data.y(index, 0), /NAN)/ct $
                        ELSE data_y(itime) =  !VALUES.F_NAN
                    ENDIF  ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR 
                data_dd(*, 30 ) = data_y         
                                ;  H+ Vy
                data_y = fltarr(n_avg)
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data.y(index, 1)) GE 0) GT 0  THEN $
                          data_y(itime) = total(data.y(index, 1), /NAN)/ct $
                        ELSE data_y(itime) =  !VALUES.F_NAN
                    ENDIF  ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR 
                data_dd(*, 31 ) = data_y
                                ;  H+ Vz
                data_y = fltarr(n_avg)
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF   TOTAL(ABS(data.y(index, 2)) GE 0) GT 0 THEN $
                          data_y(itime) = total(data.y(index, 2), /NAN)/ct $
                        ELSE data_y(itime) = !VALUES.F_NAN
                    ENDIF  ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR 
                data_dd(*, 32 ) = data_y               
                                ;  H+ Tx
                get_data, p45, data = data
                data_y = fltarr(n_avg)
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data.y(index, 0)) GE 0) GT 0  THEN $
                          data_y(itime) = total(data.y(index, 0), /NAN)/ct $
                        ELSE data_y(itime) =  !VALUES.F_NAN
                    ENDIF  ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR 
                data_dd(*, 33 ) = data_y              
                                ;  H+ Ty
                data_y = fltarr(n_avg)
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data.y(index, 1)) GE 0) GT 0 THEN $
                          data_y(itime) = total(data.y(index, 1), /NAN)/ct $
                        ELSE data_y(itime) =  !VALUES.F_NAN
                    ENDIF  ELSE data_y(itime) =  !VALUES.F_NAN
                ENDFOR 
                data_dd(*, 34 ) = data_y             
                                ;  H+ Tz
                data_y = fltarr(n_avg)
                FOR itime = 0, n_avg-1 DO BEGIN 
                    index = where(data.x GE time_dd(itime)-average_time/2 $
                                  AND data.x LT time_dd(itime)+average_time/2, ct)
                    IF ct GT 0 THEN BEGIN 
                        IF TOTAL(ABS(data.y(index, 2)) GE 0) GT 0  THEN $
                          data_y(itime) = total(data.y(index, 2), /NAN)/ct $
                        ELSE data_y(itime) = !VALUES.F_NAN
                    ENDIF  ELSE data_y(itime) = !VALUES.F_NAN
                ENDFOR 
                data_dd(*, 35 ) = data_y
                IF keyword_set(plot_mom) THEN BEGIN  
                                ;  T_x_tail
                    get_data, p46, data = data
                    data_y = data.y(index_valid)
                    data_dd(*, 36 ) = data_y
                                ;  T_y_tail
                    get_data, p47, data = data
                    data_y = data.y(index_valid)                    
                    data_dd(*, 37 ) = data_y
                                ;  T_z_tail
                    get_data, p48, data = data
                    data_y = data.y(index_valid)                  
                    data_dd(*, 38 ) = data_y
                                ;  T_x_earth
                    get_data, p49, data = data
                    data_y = data.y(index_valid)                   
                    data_dd(*, 39 ) = data_y
                                ;  T_y_earth
                    get_data, p50, data = data
                    data_y = data.y(index_valid)                 
                    data_dd(*, 40 ) = data_y
                                ;  T_z_earth
                    get_data, p51, data = data
                    data_dd(*, 41 ) = data.y(index_valid)
                ENDIF ELSE data_dd(*, 36:41) = !VALUES.F_NAN
                                ; Storm_phase
                tplot_names,p52, names=names
                if names(0) ne '' then begin 
                    get_data, p52, data = data
                    data_y = data.y(index_valid)
                    data_dd(*, 42 ) = data_y
                endif else data_y=make_array(n_avg,value=!VALUES.F_NAN)
                data_dd(*, 48) = data_y/1.e3
; eflux
; note: since eflux is very large, we devide it by 1e3 for the
; coviniency of damping data
                                ;tail eflux
                tplot_names,p03,names=names
                if names(0) ne '' then begin 
                    get_data, p03, data = data
                    data_y = data.y(index_valid, *)
                endif else data_y=make_array(n_avg,value=!VALUES.F_NAN)
                data_dd(*, 48) = data_y/1.e3 
                
                                ; earth eflux
                tplot_names,p05,names=names
                if names(0) ne '' then begin 
                    get_data, p05, data = data
                    data_y = data.y(index_valid, *)
                endif else data_y=make_array(n_avg,value=!values.f_nan)                
                data_dd(*, 49) = data_y/1.e3       
;  distribution func
; note: since dist func is very small, we multiply it by 1e6 for the
; coviniency of dumping data
                                ;tail dist func
                tplot_names,p03DF,names=names
                if names(0) ne '' then begin 
                    get_data, p03DF, data = data
                    data_y = data.y(index_valid, *)
                endif else data_y=make_array(n_avg,value=!VALUES.F_NAN)
                data_dd(*, 59) = data_y*1.e6
                
                                ; earth dist func
                tplot_names,p05DF,names=names
                if names(0) ne '' then begin 
                    get_data, p05DF, data = data
                    data_y = data.y(index_valid, *)
                endif else data_y=make_array(n_avg,value=!values.f_nan)                
                data_dd(*, 60) = data_y*1.e6  
         
                                ; tail theta
                tplot_names,p58,names=names
                if names(0) ne '' then begin 
                    get_data, p58, data = data
                    data_y = data.y(index_valid, *)
                endif else  data_y=make_array(n_avg,value=!values.f_nan)
                data_dd(*, 50) = data_y
                
                                ; earth theta
                tplot_names,p59,names=names
                if names(0) ne '' then begin 
                    get_data, p59, data = data
                    data_y = data.y(index_valid, *)
                endif else data_y = make_array(n_avg,value=!values.f_nan)
                data_dd(*, 51) = data_y 
                
            ENDIF   
            IF keyword_set(plot_imf) THEN BEGIN 
                for i_delay=0, 1 do begin 
                    if i_delay eq 0 then begin 
                        delay_location='_DC'
                        ipos=11
                    endif else begin delay_location=''
                        ipos=0
                    endelse
                                ;IMF Bx
                    get_data, p53+delay_location, data = data
                    data_y_old = data.y
                    data_y = fltarr(n_avg)
                    index = where(data_y_old EQ -999.9, ct)
                    IF ct GT 0 THEN data_y_old(index) = !VALUES.F_NAN
                    FOR itime = 0, n_avg-1 DO BEGIN
                        index = where(data.x GE time_dd(itime)-average_time/2 $
                                      AND data.x LT time_dd(itime)+average_time/2, ct)
                        IF ct GT 0 THEN BEGIN 
                            IF TOTAL(ABS(data_y_old(index)) GE 0) GT 0 THEN $
                              data_y(itime) = (total(data_y_old(index), /NAN)/ct) $
                            ELSE data_y(itime) =  !VALUES.F_NAN
                        ENDIF ELSE  data_y(itime) =  !VALUES.F_NAN
                    ENDFOR  
; There are times when ACE mag data is not good     
                    index = where(time_dd GT time_double('2001-12-24/01:08:16') $
                                  AND time_dd LT time_double('2002-03-15/01:02:30'), ct)
                    IF ct GT 0 THEN data_y(index) =  !VALUES.F_NAN
                    index = where( time_dd GT time_double('2001-04-25/00:54:00') $
                                   AND time_dd LT time_double('2001-05-22/01:15:30'), ct)
                    IF ct GT 0 THEN data_y(index) =  !VALUES.F_NAN
                    index = where( time_dd GT time_double('2002-12-10/00:50:00') $
                                   AND time_dd LT time_double('2003-01-01/00:59:00'), ct)
                    IF ct GT 0 THEN data_y(index) =  !VALUES.F_NAN
                    data_dd(*, 43 +ipos) = data_y
                                ;IMF By
                    get_data, p54+delay_location, data = data
                    data_y = fltarr(n_avg)
                    data_y_old = data.y
                    index = where(data_y_old EQ -999.9, ct)
                    IF ct GT 0 THEN data_y_old(index) = !VALUES.F_NAN
                    FOR itime = 0, n_avg-1 DO BEGIN
                        index = where(data.x GE time_dd(itime)-average_time/2 $
                                      AND data.x LT time_dd(itime)+average_time/2, ct)
                        IF ct GT 0 THEN BEGIN 
                            IF TOTAL(ABS(data_y_old(index)) GE 0) GT 0 THEN $
                              data_y(itime) = (total(data_y_old(index), /NAN)/ct) $
                            ELSE data_y(itime) = !VALUES.F_NAN
                        ENDIF ELSE data_y(itime) = !VALUES.F_NAN
                    ENDFOR    
                    index = where(time_dd GT time_double('2001-12-24/01:08:16') $
                                  AND time_dd LT time_double('2002-03-15/01:02:30'), ct)
                    IF ct GT 0 THEN data_y(index) = !VALUES.F_NAN
                    index = where( time_dd GT time_double('2001-04-25/00:54:00') $
                                   AND time_dd LT time_double('2001-05-22/01:15:30'), ct)
                    IF ct GT 0 THEN data_y(index) = !VALUES.F_NAN
                    index = where( time_dd GT time_double('2002-12-10/00:50:00') $
                                   AND time_dd LT time_double('2003-01-01/00:59:00'), ct)
                    IF ct GT 0 THEN data_y(index) = !VALUES.F_NAN
                    data_dd(*, 44 +ipos) = data_y  
                                ;IMF Bz
                    get_data, p55+delay_location, data = data
                    data_y = fltarr(n_avg)
                    data_y_old = data.y
                    index = where(data_y_old EQ -999.9, ct)
                    IF ct GT 0 THEN data_y_old(index) = !VALUES.F_NAN
                    FOR itime = 0, n_avg-1 DO BEGIN
                        index = where(data.x GE time_dd(itime)-average_time/2 $
                                      AND data.x LT time_dd(itime)+average_time/2, ct)
                        IF ct GT 0 THEN BEGIN 
                            IF   TOTAL(ABS(data_y_old(index)) GE 0) GT 0 THEN $
                              data_y(itime) = (total(data_y_old(index), /NAN)/ct) $
                            ELSE data_y(itime) = !VALUES.F_NAN
                        ENDIF ELSE data_y(itime) = !VALUES.F_NAN
                    ENDFOR  
                    index = where(time_dd GT time_double('2001-12-24/01:08:16') $
                                  AND time_dd LT time_double('2002-03-15/01:02:30'), ct)
                    IF ct GT 0 THEN data_y(index) = !VALUES.F_NAN
                    index = where( time_dd GT time_double('2001-04-25/00:54:00') $
                                   AND time_dd LT time_double('2001-05-22/01:15:30'), ct)
                    IF ct GT 0 THEN data_y(index) = !VALUES.F_NAN
                    index = where( time_dd GT time_double('2002-12-10/00:50:00') $
                                   AND time_dd LT time_double('2003-01-01/00:59:00'), ct)
                    IF ct GT 0 THEN data_y(index) = !VALUES.F_NAN
                    data_dd(*, 45 +ipos) = data_y   
                                ;SW V
                    get_data, p56+delay_location, data = data
                    data_y = fltarr(n_avg)
                    data_y_old = data.y
                    index = where(data_y_old EQ -999.9, ct)
                    IF ct GT 0 THEN data_y_old(index) = !VALUES.F_NAN
                    FOR itime = 0, n_avg-1 DO BEGIN
                        index = where(data.x GE time_dd(itime)-average_time/2 $
                                      AND data.x LT time_dd(itime)+average_time/2, ct)
                        IF ct GT 0 THEN BEGIN 
                            IF  TOTAL(ABS(data_y_old(index)) GE 0) GT 0 THEN $
                              data_y(itime) = (total(data_y_old(index), /NAN)/ct) $
                            ELSE data_y(itime) = !VALUES.F_NAN
                        ENDIF ELSE  data_y(itime) = !VALUES.F_NAN
                    ENDFOR  
                    data_dd(*, 46 +ipos) = data_y   
                                ;SW P
                    get_data, p57+delay_location, data = data
                    data_y = fltarr(n_avg)
                    data_y_old = data.y
                    index = where(data_y_old EQ -999.9, ct)
                    IF ct GT 0 THEN data_y_old(index) = !VALUES.F_NAN
                    FOR itime = 0, n_avg-1 DO BEGIN
                        index = where(data.x GE time_dd(itime)-average_time/2 $
                                      AND data.x LT time_dd(itime)+average_time/2, ct)
                        IF ct GT 0 THEN BEGIN 
                            IF TOTAL(ABS(data_y_old(index)) GE 0) GT 0 THEN $
                              data_y(itime) = (total(data_y_old(index), /NAN)/ct) $
                            ELSE data_y(itime) = !VALUES.F_NAN
                        ENDIF ELSE data_y(itime) = !VALUES.F_NAN
                    ENDFOR  
                    data_dd(*, 47 +ipos) = data_y 
                endfor 
                                 
            ENDIF  
;dump data using routine dump_data
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
ENDIF        ELSE BEGIN  
    OPENU, unit, path+'log_errors.txt', /GET_LUN, /APPEND
    PRINTF, unit, ts + ' TO '+ te, '-----BEAM NOT FOUND------'
    FREE_LUN, unit         
ENDELSE 
close, /all
END 
