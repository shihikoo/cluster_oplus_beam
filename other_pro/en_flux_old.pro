PRO en_flux, recalc = recalc, all_recalc = all_recalc, ps = ps, single_plot = single_plot, $
             sc = sc, en_perp_plot = en_perp_plot, flux_plot = flux_plot, $
             cut_spec_plot = cut_spec_plot, t_perp_plot = t_perp_plot, pa_plot = pa_plot,$
             output=output,normalize_flux=normalize_flux, ntime=ntime,add_sw_dst=add_sw_dst,$
             bl=bl

storm_time_group=['2002-09-10/13:00:00','2002-09-30/10:00:00','2002-10-17/10:00:00','2001-09-23/10:00:00','2001-10-11/10:00:00','2001-10-27/14:00:00',$
                  '2004-10-13/12:00:00','2005-08-25/18:00:00']
storm_time_duration_group=[24,36,28,24,24,36,36,36]

nonstorm_time_group=['2001-07-04/12','2001-07-14','2001-07-25/22','2001-08-14','2001-09-21/1','2001-09-22/04:00:00']
nonstorm_time_duration_group=[12,18,18,20,27,30]

time_group_Bogdanova = ['2001-07-02','2001-07-30','2001-08-04','2001-08-23','2001-09-28','2001-07-16','2001-09-28','2001-09-30','2001-10-19','2001-11-14']
time_duration_group_Bogdanova = [36,36,36,36,36,36,36,36,36,36]

time_group = nonstorm_time_group
time_duration_group = nonstorm_time_duration_group
if not keyword_set(ntime) then ntime=0
if not keyword_set(add_sw_dst) then add_sw_dst=1
if not keyword_set(normalize_flux)then normalize_flux=1
if not keyword_set(add_beta) then add_beta = 1 
if not keyword_set(add_pa_for_high_low_energy) then add_pa_for_high_low_energy=1
if not keyword_set(units_name_input) then units_name_input = 'DIST FUNC' ;'DIFF FLUX'

for itime = ntime,n_elements(time_group)-1 do begin
    time_duration=time_duration_group(itime)
    IF NOT keyword_set(sc) THEN sc = 4
    if total(time_group eq storm_time_group) eq n_elements(time_group) then begin 
        if itime eq 0 then begin
            if sc eq 4 then cut_time = ['2002-09-10/14:09:10', '2002-09-10/14:31:00', '2002-09-10/15:00:00', '2002-09-10/16:25:43', $
                                        '2002-09-10/18:37:45', '2002-09-10/20:27:03', '2002-09-10/23:26:27', '2002-09-11/01:03:45', $
                                        '2002-09-11/01:41:43', '2002-09-11/03:50:00', '2002-09-11/06:32:00', '2002-09-11/07:54:30', $
                                        '2002-09-11/09:07:00']
            if sc eq 4 and keyword_set(bl) then  cut_time = ['2002-09-10/14:09:10', '2002-09-11/07:53:05', '2002-09-11/07:57:00']
            if sc eq 1 then cut_time = ['2002-09-10/13:57:12', '2002-09-10/14:29:52', '2002-09-10/16:29:55', '2002-09-10/18:02:25', $
                                        '2002-09-10/20:16:40', '2002-09-10/23:11:10', '2002-09-11/01:03:08', '2002-09-11/01:59:52', $
                                        '2002-09-11/06:22:05', '2002-09-11/06:32:20', '2002-09-11/07:52:40', '2002-09-11/09:09:52']
        endif  
        if itime eq 1 then cut_time =  ['2002-10-01/21:32:30', $
                                        '2002-09-30/11:14:00', '2002-09-30/11:19:45', '2002-09-30/12:19:47', '2002-09-30/12:34:25', $
                                        '2002-09-30/13:27:00', '2002-09-30/14:10:00', '2002-09-30/16:23:00', '2002-09-30/16:34:00', $
                                        '2002-09-30/17:33:30', '2002-09-30/18:05:30', '2002-09-30/22:19:30', '2002-10-01/04:35:00', $
                                        '2002-10-01/14:30:00', '2002-10-01/20:44:30']
        if itime eq 2 then cut_time = ['2002-10-18/13:40:00', $
                                       '2002-10-17/12:54:15', '2002-10-17/16:38:20', '2002-10-17/17:32:20', '2002-10-17/18:09:00', $
                                       '2002-10-17/22:06:20', '2002-10-18/00:24:20', '2002-10-18/02:35:00', '2002-10-18/06:29:00', $
                                       '2002-10-18/08:22:10', '2002-10-18/08:30:10', '2002-10-18/11:05:40', '2002-10-18/11:55:00']   
        if itime eq 3 then cut_time = ['2001-09-23/09:40:00', $
                                       '2001-09-17/12:54:15', '2001-09-17/16:38:20', '2001-09-17/17:32:20', '2001-09-17/18:09:00', $
                                       '2001-09-17/22:06:20', '2001-09-18/00:24:20', '2001-09-18/02:35:00', '2001-09-18/06:29:00', $
                                       '2001-09-18/08:22:10', '2001-09-18/08:30:10', '2001-09-18/11:05:40', '2001-09-18/11:55:00']
        if itime eq 4 then cut_time = ['2001-10-12/10:0:00', $
                                       '2001-10-17/12:54:15', '2001-09-17/16:38:20', '2001-09-17/17:32:20', '2001-09-17/18:09:00', $
                                       '2001-10-17/22:06:20', '2001-09-18/00:24:20', '2001-09-18/02:35:00', '2001-09-18/06:29:00', $
                                       '2001-10-18/08:22:10', '2001-09-18/08:30:10', '2001-09-18/11:05:40', '2001-09-18/11:55:00']
        if itime eq 5 then cut_time = ['2001-10-29/01:00:00', $
                                       '2001-10-17/12:54:15', '2001-09-17/16:38:20', '2001-09-17/17:32:20', '2001-09-17/18:09:00', $
                                       '2001-10-17/22:06:20', '2001-09-18/00:24:20', '2001-09-18/02:35:00', '2001-09-18/06:29:00', $
                                       '2001-10-18/08:22:10', '2001-09-18/08:30:10', '2001-09-18/11:05:40', '2001-09-18/11:55:00']
        if itime eq 6 then cut_time = ['2004-10-14/21:00:00', $
                                       '2001-10-17/12:54:15', '2001-09-17/16:38:20', '2001-09-17/17:32:20', '2001-09-17/18:09:00', $
                                       '2001-10-17/22:06:20', '2001-09-18/00:24:20', '2001-09-18/02:35:00', '2001-09-18/06:29:00', $
                                       '2001-10-18/08:22:10', '2001-09-18/08:30:10', '2001-09-18/11:05:40', '2001-09-18/11:55:00']
        if itime eq 7 then cut_time = ['2005-08-25/20:00:00', $
                                       '2001-10-17/12:54:15', '2001-09-17/16:38:20', '2001-09-17/17:32:20', '2001-09-17/18:09:00', $
                                       '2001-10-17/22:06:20', '2001-09-18/00:24:20', '2001-09-18/02:35:00', '2001-09-18/06:29:00', $
                                       '2001-10-18/08:22:10', '2001-09-18/08:30:10', '2001-09-18/11:05:40', '2001-09-18/11:55:00']
    endif
    if total(time_group eq time_group_Bogdanova) eq n_elements(time_group) then begin 
        if itime eq 4 then cut_time = ['2001-09-28/05:50:16', '2001-09-28/06:02:40', '2001-09-28/07:44:14', '2001-09-28/08:42:08s', $
                                       '2001-09-29/09:13:00', '2001-09-29/09:23:50', '2001-09-29/09:30:30', '2001-09-29/09:41:50', $
                                       '2001-09-29/09:59:30', '2001-09-29/10:03:50', '2001-09-29/10:25:26', '2001-09-29/10:23:26', $
                                       '2001-09-29/10:35:36']
    endif 
    if total(time_group eq nonstorm_time_group) eq n_elements(time_group) then begin 
        if itime eq 0 then cut_time = ['2001-07-04/13:40:34', $
                                       '2001-07-04/14:09:40', '2001-07-04/15:01:40', '2001-07-04/15:38:55', $
                                       '2001-07-04/15:44:30', '2001-07-04/16:52:10', '2001-07-04/17:43:45', '2001-07-04/17:54:22', $
                                       '2001-07-04/17:59:06', '2001-07-04/18:42:26', '2001-07-04/18:49:02', '2001-07-04/19:31:10', $
                                       '2001-07-04/19:53:40']
        if itime eq 0 and keyword_set(bl) then cut_time = ['2001-07-04/13:40:34', $
                                                           '2001-07-04/17:48:00', '2001-07-04/17:45:40', '2001-07-04/17:57:04',$
                                                           '2001-07-04/20:15:30','2001-07-04/20:34:10','2001-07-04/21:34:10']
        if itime eq 1 then cut_time = ['2001-07-14/01:31:45', $
                                       '2001-07-14/01:52:30', '2001-07-14/03:33:00', '2001-07-14/04:56:00', $
                                       '2001-07-14/06:22:00', '2001-07-14/07:13:40', '2001-07-14/07:43:00', '2001-07-14/08:19:20', $
                                       '2001-07-14/08:43:10', '2001-07-14/09:00:30', '2001-07-14/09:25:10', '2001-07-14/15:35:20', $
                                       '2001-07-14/15:58:30']
        if itime eq 1 and keyword_set(bl) then cut_time = ['2001-07-14/01:31:45', $
                                                           '2001-07-14/07:22:40']
        if itime eq 2 then cut_time = ['2001-07-25/23:10:04', $
                                       '2001-07-25/23:19:10', '2001-07-25/23:28:20', '2001-07-25/23:37:10', $
                                       '2001-07-26/00:06:10', '2001-07-26/00:27:20', '2001-07-26/01:23:00', '2001-07-26/03:23:40', $
                                       '2001-07-26/04:51:30', '2001-07-26/06:10:40', '2001-07-26/06:17:30', '2001-07-26/06:36:00', $
                                       '2001-07-26/07:05:30']
        if itime eq 2 and keyword_set(bl) then cut_time = ['2001-07-25/23:10:04', $
                                                           '2001-07-26/07:15:00','2001-07-26/07:31:20']
        if itime eq 3 then cut_time = ['2001-08-14/00:28:24', $
                                       '2001-08-14/00:55:40', '2001-08-14/01:10:10', '2001-08-14/01:29:14', $
                                       '2001-08-14/02:00:30', '2001-08-14/02:29:30', '2001-08-14/03:28:32', '2001-08-14/05:05:00', $
                                       '2001-08-14/06:07:00', '2001-08-14/07:19:00', '2001-08-14/07:29:30', '2001-08-14/13:35:00', $
                                       '2001-08-14/16:14:05']
        if itime eq 3 and keyword_set(bl) then cut_time = ['2001-08-14/00:28:24', $
                                                           '2001-08-14/15:23:25','2001-08-14/18:01:00']
        if itime eq 4 then cut_time = ['2001-09-21/02:51:04s', $
                                       '2001-09-21/03:07:30', '2001-09-21/03:50:30', '2001-09-21/17:12:50', $
                                       '2001-09-21/18:10:50', $
                                       '2001-09-21/18:56:10', '2001-09-21/19:47:00', '2001-09-21/20:43:20', $
                                       '2001-09-21/21:54:20', '2001-09-21/22:57:00']
        if itime eq 4 and keyword_set(bl) then cut_time = ['2001-09-21/02:51:04',$
                                                           '2001-09-21/18:27:40', '2001-09-21/18:35:50', $
                                                           '2001-09-21/18:56:10','2001-09-21/19:49:38']
        if itime eq 5 then cut_time = ['2001-09-23/08:05:44', $
                                       '2001-09-23/07:20:10', '2001-09-23/02:36:16', '2001-09-23/01:44:00', '2001-09-22/18:37:40', $
                                       '2001-09-22/14:20:30', '2001-09-22/11:58:55', '2001-09-22/11:05:00','2001-09-22/09:41:47', $
                                       '2001-09-22/08:45:46', '2001-09-22/08:24:28', '2001-09-22/08:11:45', '2001-09-22/07:51:04']
        if itime eq 5 and keyword_set(bl) then cut_time = ['2001-09-23/08:05:44', $
                                                            '2001-09-22/08:41:44', '2001-09-22/08:22:28', '2001-09-22/08:09:07']
        
    endif 
    avg_time = 40
    at_str = strcompress(avg_time, /remove_all)
    
    sc_str = strcompress(sc, /remove_all)
    mass_o = 16*1.67262158e-27
    tplot_names, '*SC'+sc_str+'*AVG*', names = names

    time = time_group(itime)
    timespan, time, time_duration, /hours ; SECONDS, MINUTES, HOURS, DAYS (DEFAULT)    
    ts = time
    date_str = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
    plot_path = 'output/o_beam/en_flux/'+date_str+'/'
    spawn, 'mkdir '+plot_path
    
    beta_name='TDMOM_EN00040_40000_SC'+sc_str+'_MTPRESSURE_SP0_ET0_All_O1_beta'
    p_total='TDMOM_EN00040_40000_SC'+sc_str+'_MTPRESSURE_SP0_ET0_All_O1_P_total'
    proton_en_name='ENSPEC_SC'+sc_str+'_IN0_PHI0_360_UN'+STRUPCASE(strcompress(units_name_input,/remove_all))+'_SP0_ET0_All'
    enspec_name = 'ENSPEC_SC'+sc_str+'_IN0_PHI0_360_UN'+STRUPCASE(strcompress(units_name_input,/remove_all))+'_SP3_ET0_All'
    energy_low = [35., 1000.] &   energy_high = [1000., 40000.]
    energy_low_str = [STRING(min(energy_low), format = '(i5.5)'), STRING(max(energy_low), format = '(i5.5)') ]
    energy_high_str = [STRING(min(energy_high), format = '(i5.5)'), STRING(max(energy_high), format = '(i5.5)') ]
    paspec_name_low = 'PASPEC_EN'+ energy_low_str(0)+'_'+energy_low_str(1)+'_SC'+sc_str+'_UN'+STRUPCASE(strcompress(units_name_input,/remove_all))+'_SP3_All'
    paspec_name_high = 'PASPEC_EN'+ energy_high_str(0)+'_'+energy_high_str(1)+'_SC'+sc_str+'_UN'+STRUPCASE(strcompress(units_name_input,/remove_all))+'_SP3_All'

    v_perp_name_edi = 'Vvec_t'+'_AVG'+at_str 
    v_name_tot = 'TDMOM_EN00040_40000_SC' + sc_str+'_MTVELOCITY_SP3_ET0_All'
    v_perp_name = v_name_tot+'_AVG'+at_str+'_V_PERP_T'
    v_para_name = v_name_tot+'_AVG'+at_str+'_V_PAR_T'
    t_name_tot = 'TDMOM_EN00040_40000_SC'+sc_str+'_MTTEMPERATURE_SP3_ET0_All'
;en_set = [31444.7, 19398.3, 11966.9, 7382.39, 4554.22, 2809.51, 1733.19, 1069.21, $s
;          659.599, 406.909, 251.023, 154.857, 95.5315, 58.9337, 36.3563] 
    erange_set = [[30000., 40000.], [15000., 30000.], [10000., 15000.], [8000., 10000.], $
                  [4000., 8000.], [2000., 3000.], [1500., 2000.], [1000., 1500.], [500., 1000.], $
                  [400., 500.], [200., 400.], [100., 200.], [60., 100.], [50., 60.], [30., 50.]]

    erange_set_str = STRING(long(erange_set), format = '(i5.5)')
    fln = plot_path+'save_sc'+sc_str
    var_label = 'EPH_SC'+STRING(sc, FORMAT = '(i1.1)')+'_' $
      + ['MLT', 'GSM_X', 'GSM_Y', 'GSM_Z', 'DIST']
    mag_name = 'MAG_SC'+sc_str+'_B_xyz_gse_T'

;load locations
    get_cluster_ephemeris, sc, /GSE_X, /GSE_Y, /GSE_z, /GSM_X, /GSM_Y, /GSM_Z, /DIST, /MLT
    
    IF NOT keyword_set(names) OR keyword_set(recalc) THEN BEGIN
        print, FINDFILE(fln+'.tplot', COUNT = count)
        tplot_names, names = names
        store_data, delete = names
        IF count GT 0 AND NOT keyword_set(all_recalc) THEN tplot_restore, file = fln+'.tplot' else begin 
;load energy spectra
            sat = [sc,sc]
            units_name = units_name_input ; 'Counts', 'NCOUNTS', 'RATE', 'NRATE', 'DIFF FLUX', 'EFLUX'
            inst = 0            ; 0: CODIF, 1: HIA
            angle = [[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
            eff_table = 0       ; 0: GROUND, 1: ONBOARD
            specie=[0,3]
            
            plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 1
;load locations 
            sat=sc
            get_cluster_ephemeris, sat, /GSE_X, /GSE_Y, /GSE_z, /GSM_X, /GSM_Y, /GSM_Z, /DIST, /MLT
;average energy spectra 
            average_tplot_variable, enspec_name, avg_time, /new
            options, enspec_name+'_AVG'+at_str, 'ytitle', 'SC'+sc_str $
              +'!C!CO!U+!N(eV)!C!CAVG-'+strcompress(avg_time, /remove_all)+'s'

            ylim, 'ENSPEC*', 40, 4e4, 1
            zlim, 'ENSPEC*', 0.1, 100, 1
; plot pitch angle for high energy range and low energy range and
; average them
            if keyword_set(add_pa_for_high_low_energy) then begin 
                sat=sc
                specie=3
                inst=0
                units_name= units_name_input
                plot_pa_spec_from_crib, sat, specie, inst, units_name, $
                  energy_low, eff_table,  recalc = 1, COMBINE = 1  
                tplot_names, paspec_name_low, names = names
                paspec_name_low = names(0)
                
                plot_pa_spec_from_crib, sat, specie, inst, units_name, $
                  energy_high, eff_table,  recalc = 1, COMBINE = 1
                tplot_names, paspec_name_high, names = names
                paspec_name_high = names(0)

                average_tplot_variable, paspec_name_low, avg_time, /new
                average_tplot_variable, paspec_name_high, avg_time, /new

                zlim, 'PASPEC*', 1, 1e4, 1
                ylim, 'PASPEC*', 0, 180, 0
                options, paspec_name_low+'_AVG'+at_str, 'ytitle', 'SC'+sc_str $
                  +'!C!CO!U+!N(deg)!C!C'+energy_low_str(0)+' - ' $
                  +energy_low_str(1) +'!C!CAVG-'+strcompress(avg_time, /remove_all)+'s'
                options, paspec_name_high+'_AVG'+at_str, 'ytitle', 'SC'+sc_str $
                  +'!C!CO!U+!N(deg)!C!C'+energy_high_str(0)+' - ' $
                  +energy_high_str(1) +'!C!CAVG-'+strcompress(avg_time, /remove_all)+'s'
                sat = sc
                plot_mag_from_crib, sat
                average_tplot_variable, mag_name, avg_time, /new 
            endif  
            
; add proton energy spectrum and plasma beta information
            if keyword_set(add_beta)then begin
               
                average_tplot_variable,proton_en_name,avg_time,/new

                sat = [sc, sc]
                specie = [0, 3]
                moments = ['P', 'P']
                angle = [[-90, 90], [0, 360]]
                energy = [40., 40000.]
                inst = 0 &  eff_table = 0
                plot_3dmom_from_crib, sat, specie, inst, moments, angle, energy, eff_table, recalc = 0             
                
                h_press =  'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP0_ET0_All'
                o_press = 'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP3_ET0_All'
                mag_press =  'MAG_SC' + sc_str +'_B_xyz_gse_MAG_PR'
                hia_press =  'HIA_L2_MOM_SC'+sc_str+'_MTPRESSURE_SP0_PR202'
                
        ;        tplot_names,names=mag_press
         ;       if names(0) eq '' then begin 
          ;          sat=sc
           ;         plot_mag_from_crib,sat
            ;    endif                    

                h1_press = h_press & o1_press = o_press
                
                plasma_beta, h1_press, mag_press, O1_PRESSURE = o1_press 
                
                                ;      tplot_names, 'TDMOM*SP3*', names = names
                                ;     store_data, delete = names
                beta_name='TDMOM_EN00040_40000_SC4_MTPRESSURE_SP0_ET0_All_O1_beta'
                average_tplot_variable, beta_name,avg_time,/new
                p_total='TDMOM_EN00040_40000_SC4_MTPRESSURE_SP0_ET0_All_O1_P_total'
                average_tplot_variable,p_total,avg_time,/new
            endif   

; load pitch angle plots for all different energy bins
            if keyword_set(add_pa_for_diff_energy) then begin 
                FOR ie = 0, n_elements(erange_set(0, *))-1 DO BEGIN 
                    sat = sc &  specie = 3
                    inst = 0 &  eff_table = 0
                    units_name = units_name_input
                    energy = erange_set(*, ie)
                    
                    plot_pa_spec_from_crib, sat, specie, inst, units_name, $
                      energy, eff_table,  recalc = 1, COMBINE = 1  
                    pa_name = 'PASPEC_EN'+ erange_set_str(2*ie)+'_'+erange_set_str(2*ie+1)+'_SC'+sc_str+'_UN'+STRUPCASE(strcompress(units_name_input,/remove_all))+'_SP3_All'
                    average_tplot_variable, pa_name, avg_time, /new
                ENDFOR
            endif
; load velocity and temperature for all different energy bins
            if keyword_set(add_V_T_for_diff_energy) then begin 
                FOR ie = 0, n_elements(erange_set(0, *))-1 DO BEGIN 
                    sat = [sc, sc] &  specie = [3, 3]
                    inst = 0 &  eff_table = 0
                    moments = ['V', 'T']
                    angle = [[-90, 90], [0, 360]] 
                    energy = [erange_set(0, (ie+1) < 14), erange_set(1, (ie-1) > 0)]

                    plot_3dmom_from_crib, sat, specie, inst, moments, angle, energy, eff_table, recalc = 1
                    V_name = 'TDMOM_EN'+erange_set_str(2*(ie+1 < 14))+'_'+erange_set_str(2*(ie-1 > 0)+1)$
                      +'_SC' + sc_str+'_MTVELOCITY_SP3_ET0_All'
                    T_name = 'TDMOM_EN'+erange_set_str(2*(ie+1 < 14))+'_'+erange_set_str(2*(ie-1 > 0)+1)$
                      +'_SC' + sc_str+'_MTTEMPERATURE_SP3_ET0_All'
                    average_tplot_variable, V_name, avg_time, /new
                    v_perp, V_name+'_AVG'+at_str
                    average_tplot_variable, T_name, avg_time, /new
                ENDFOR 
            endif  
; load velocity and temperature for all energy 
            if keyword_set(add_V_T) then begin
                sat = [sc, sc]  &  specie = [3, 3]
                inst = 0 &  eff_table = 0
                moments = ['V', 'T']
                angle = [[-90, 90], [0, 360]]
                energy = [40, 40000]
                plot_3dmom_from_crib, sat, specie, inst, moments, angle, energy, eff_table, recalc = 1

                average_tplot_variable, V_name_tot, avg_time, /new
                v_perp, V_name_tot+'_AVG'+at_str    

                average_tplot_variable, T_name_tot, avg_time, /new
            endif 
; load edi data if using satallite one data
            IF sc EQ 1 THEN BEGIN 
                scnum_in = sc_str
                path_in = 'edidata'
                print, extract_edi_cdf_data( '20020910', scnum_in, $ ; INPUT
                                             epoch, time, vvec, evec, $ ; OUTPUT
                                             path = path_in, quality = quality)
                time_combine = time_double('2002-09-10:/00')+time
                vvec_combine = transpose(vvec)
                evec_combine = transpose(evec)
                quality_combine = quality
                date_in = '20020910'
                
                print, extract_edi_cdf_data( '20020911', scnum_in, $ ; INPUT
                                             epoch, time, vvec, evec, $ ; OUTPUT
                                             path = path_in, quality = quality)
                time_combine = [time_combine, time_double('2002-09-11/00')+time]
                vvec_combine = [vvec_combine, transpose(vvec)]
                evec_combine = [evec_combine, transpose(evec)]
                quality_combine = [quality_combine, quality]

                                ;   vvec_combine(where(quality_combine EQ 1), *) = !VALUES.F_NAN
                                ;  evec_combine(where(quality_combine EQ 1), *) = !VALUES.F_NAN

                store_data, 'quality_edi', data = {x:time_combine, y:quality_combine}
                store_data, 'Vvec', data = {x:time_combine, y:vvec_combine}
                store_data, 'Vvec_x', data = {x:time_combine, y:vvec_combine(*, 0)}
                store_data, 'Vvec_y', data = {x:time_combine, y:vvec_combine(*, 1)}
                store_data, 'Vvec_z', data = {x:time_combine, y:vvec_combine(*, 2)}
                store_data, 'Vvec_t', data = {x:time_combine, y:(sqrt(vvec_combine(*, 0)^2+vvec_combine(*, 1)^2+vvec_combine(*, 2)^2))}
                store_data, 'Evec', data = {x:time_combine, y:evec_combine}
                store_data, 'Evec_x', data = {x:time_combine, y:evec_combine(*, 0)}
                store_data, 'Evec_y', data = {x:time_combine, y:evec_combine(*, 1)}
                store_data, 'Evec_z', data = {x:time_combine, y:evec_combine(*, 2)}
                store_data, 'Evec_t', data = {x:time_combine, y:(sqrt(evec_combine(*, 0)^2+evec_combine(*, 1)^2+evec_combine(*, 2)^2))}
                average_tplot_variable, 'Vvec_t', avg_time, /new
                average_tplot_variable, 'quality_edi', avg_time, /new

                options, v_para_name, 'ytitle', 'SC4!C!CO!U+!N V!D//!N (kms!U-1!N)'
                options, v_perp_name, 'ytitle', 'SC4!C!CO!U+!N V!DL!N (kms!U-1!N)'
                options, v_perp_name_edi, 'color', 2
                options, [v_perp_name, v_perp_name_edi, 'Vvec_t', 'quality_edi'], 'psym'
                ylim, v_perp_name, 0, 100, 0
                ylim, v_perp_name_edi, 0, 100, 0
                ylim, 'Vvec_t', 0, 100, 0
                ylim, 'quality_edi', 0, 3
                average_tplot_variable, 'EPH_SC'+sc_str+'_GSE_X', avg_time, /new
            ENDIF 
; plot temperature for all energy range 
;        sat = sc  &  specie = 3
                                ;       inst = 0 &  eff_table = 0
                                ;      moments = ['T']
                                ;     angle = [[-90, 90], [0, 360]]
                                ;    energy = [40, 40000]
                                ;   plot_3dmom_from_crib, sat, specie, inst, moments, angle, energy, eff_table, recalc = 1
                                ;  average_tplot_variable, T_name_tot, avg_time, /new
            
; plot temperature for different energy range
;    FOR ie = 0, n_elements(erange_set(0, *))-1 DO BEGIN 
                                ;       sat = sc &  specie = 3
                                ;      inst = 0 &  eff_table = 0
                                ;     units_name = units_name_input
                                ;    energy = [erange_set(0, (ie+1) < 14), erange_set(1, (ie-1) > 0)]
                                ;   
                                ;  angle = [[-90, 90], [0, 360]]     
                                ; moments = ['T']
                                ;       plot_3dmom_from_crib, sat, specie, inst, moments, angle, energy, eff_table, recalc = 1
                                ;      
                                ;     T_name = 'TDMOM_EN'+erange_set_str(2*(ie+1 < 14))+'_'+erange_set_str(2*(ie-1 > 0)+1)$
                                ;             +'_SC' + sc_str+'_MTTEMPERATURE_SP3_ET0_All'
            
                                ;   average_tplot_variable, T_name, avg_time, /new
                                ; ENDFOR
 
     

;        sat = sc
 ;       units_name = units_name_input ; 'Counts', 'NCOUNTS', 'RATE', 'NRATE', 'DIFF FLUX', 'EFLUX'
  ;      inst = 0                ; 0: CODIF, 1: HIA
   ;     angle = [[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
    ;    eff_table = 0           ; 0: GROUND, 1: ONBOARD
     ;   plot_en_spec_from_crib, sat, 3, inst, units_name, angle, eff_table, recalc = 1

                                ;      IF keyword_set(ps) THEN popen, plot_path+'SC'+sc_str+'_all_time_avg.ps', /land ELSE window, 1
        
                                ;     tplot, [enspec_name+'_AVG'+at_str, paspec_name_low+'_AVG'+at_str, paspec_name_high+'_AVG'+at_str $
                                ;           ],var_label=var_label
                                ;          ,mag_name+'_AVG'+at_str, v_perp_name], var_label = var_label
;
                                ;       IF sc EQ 1 THEN tplot_panel, v = v_perp_name, o = v_perp_name_edi
                                ;      IF keyword_set(ps) THEN pclose
                                ;      ELSE stop
        if keyword_set(add_sw_dst) then begin 
            read_omni,hr=1,all=1
            read_omni
        endif
        tplot_save, file = fln   
    endelse    
ENDIF          

    IF keyword_set(single_plot) THEN BEGIN  
        FOR icut = 0, n_elements(cut_time)-1 DO BEGIN 
            cal_time = time_double(cut_time(icut))
            timespan, cal_time-150, 300, /seconds
            zlim, [enspec_name, enspec_name+'_AVG'+at_str], 0.1, 100, 1
            get_data, enspec_name+'_AVG'+at_str, data = dd_avg
            ind_avg = sort(ABS(dd_avg.x-cal_time))
            get_data, enspec_name, data = dd
            ind = sort(ABS(dd.x - cal_time))

            te = time_string(cal_time)
            time_str = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2)+'_'+ STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2)
            
            IF keyword_set(ps) THEN popen, plot_path+time_str+'_spectra'+region_str+'.ps', /port ELSE window, 2
            options, '*', 'panel_size', 1
            options, [v_perp_name, v_perp_name_edi, 'Vvec_t', 'quality_edi'], 'psym', -1
            options, v_perp_name_edi, 'color', 2
            tplot, [enspec_name, enspec_name+'_AVG'+at_str, v_perp_name, 'Vvec_t', 'quality_edi'], var_label = var_label
            tplot_panel, v = v_perp_name, o = v_perp_name_edi, psym = -7
            timebar, dd_avg.x(ind_avg(0))-avg_time/2
            timebar, dd_avg.x(ind_avg(0))+avg_time/2
            IF keyword_set(ps) THEN pclose ELSE stop

            num = 10
            time =  dd.x((ind(0)-num/2):(num/2+ind(0))-1)
            flux = dd.y((ind(0)-num/2):(num/2+ind(0)-1), 0:14)
            en = dd.v(ind(0)-num/2:num/2+ind(0), 0:14)

            IF keyword_set(ps) THEN  popen, plot_path+time_str+'_en_flux.ps', /land  ELSE window, 3
            plot, [0, 0], [0, 0], xlog = 1, xstyle = 1, ylog = 1, ystyle = 1, xrange = [40, 4e4], yrange = [1e-2, 3e3], $
              title = time_string(time_string(cal_time)), $
              /nodata

            FOR i = 0, num-1 DO BEGIN 
                oplot, en(i, *), flux(i, *), color = 1+i/2, psym = -1
                xyouts, 1000, 2/(10^(i/3.)), time_string(time(i)), color = 1+i/2, charsize = 1.2
            ENDFOR 
            oplot, dd_avg.v(ind_avg(0), 0:14), dd_avg.y(ind_avg(0), 0:14), thick = 3
            IF keyword_set(ps) THEN pclose ELSE stop
        ENDFOR 
    ENDIF  
    timespan, time, time_duration, /hours 

    get_data, enspec_name+'_AVG'+at_str, data = dd_avg
    get_data,  'MAG_SC'+sc_str+'_B_xyz_gse_T_AVG'+at_str, data = data
    mag_time = data.x
    mag_avg = data.y

    na = n_elements(cut_time)

    IF keyword_set(flux_plot) THEN BEGIN 
        IF keyword_set(en_over_b) THEN BEGIN 
            x_name= 'en_over_B' & x_title='Energy (eV) / B' & xrange = [2e-3,2e2] & xyout_x = 0.1
            y_name='flux' & y_title = 'Differential flux  (1/cm!U2!N-s-sr-(eV/e))' & yrange =  [6e-3, 8e3]
        ENDIF else begin 
            x_name = 'en' & x_title = 'Energy (eV)' & xrange = [30, 4e4] & xyout_x = 36
            y_name='flux' & y_title = 'Differential flux  (1/cm!U2!N-s-sr-(eV/e))' & yrange =  [6e-3, 8e3]
        ENDELSE 
        if units_name_input eq 'DIST FUNC' then begin 
            y_name='flux' & y_title = 'Distribution Function (s!E3!N/cm!E3!N-km!E3!N)' & yrange =  [1e-11, 1e-5]
        endif 

        if keyword_set(bl) then region_str='_Boundary_Layer' else region_str=''
        IF keyword_set(ps) THEN begin 
            if keyword_set(normalize_flux) then normalize_str='_normalized' else normalize_str=''
            if keyword_set(output) then popen,output,/land $
            else popen, plot_path+'sc'+sc_str+'_'+x_name+'_flux'+normalize_str +region_str+'.ps', /land 
        endif ELSE window, 4
       
        if keyword_set(normalize_flux) then title='Normalized Flux '+region_str else title=''
        
        plot, [0, 0], [0, 0], xrange = xrange, yrange = yrange, xlog = 1, ylog = 1, charsize = 1.3, $ 
          xstyle = 1, ystyle = 1, xtitle = x_title, ytitle = y_title, title=title, position=[0.15,0.15,0.95,0.9],/nodata
        FOR ia = 0, na-1 DO BEGIN
            index = sort(ABS(dd_avg.x - time_double(cut_time(ia)))) 
            index_mag = sort(ABS(mag_time - time_double(cut_time(ia))))
            x = dd_avg.v(index(0), 0:14)
            y = dd_avg.y(index(0),0:14)
            IF keyword_set(en_over_b) THEN BEGIN
                n_avg = n_elements(dd_avg.x)
                v_para = fltarr(n_avg, 15)
                v_perp = fltarr(n_avg, 15)
                for ie = 0, n_elements(erange_set(0,*))-1 DO BEGIN 
                    V_name = 'TDMOM_EN'+erange_set_str(2*ie)+'_'+erange_set_str(2*ie+1) +'_SC' + sc_str+'_MTVELOCITY_SP3_ET0_All'
                    get_data, v_name+'_AVG'+at_str+'_V_PAR_T', data = data
                    v_para(*, ie) = data.y
                    get_data, v_name+'_AVG'+at_str+'_V_PERP_T', data = data
                    v_perp(*, ie) = data.y
                    tplot_names, '*'+at_str+'_V*'
                ENDFOR
                x = fltarr(15)
                for ie = 0, n_elements(erange_set(0,*))-1 DO  $
                  x(ie) = mass_o*(v_perp(index(0), ie)*1e3)^2/2/1.6e-19/mag_avg(index_mag(0))
            ENDIF

            if keyword_set(normalize_flux) then begin
                if ia eq 0 then mag_normal=mag_avg(index_mag(0)) else y=y/(mag_avg(index_mag(0))/mag_normal)
            endif 
            if ia eq 0 then thickness= 10 else thickness=4
            if ia eq 0 and keyword_set(ps) eq 0 then color_input = -1 else color_input = (ia+1)/2

            oplot, x, y, color = color_input, psym = -(ia-2*fix(ia/2)+1), $
              linestyle =  ABS((ia-2*fix((ia)/2))*2), thick = thickness
                                ;     average_tplot_variable, T_name, avg_time, /new
            IF  ia-2*fix(ia/2) EQ 0 THEN line_style = ' __' ELSE line_style = ' ---'
            xyout_y = 0.3/(10^(ia/5.7))
            xyouts, xyout_x, xyout_y, time_string(dd_avg.x(index(0)))+line_style, color = color_input, charsize = 1.2
        ENDFOR  
        IF keyword_set(ps) THEN pclose ELSE stop        
    ENDIF
    
    IF keyword_set(cut_spec_plot) THEN BEGIN 
        IF keyword_set(ps) THEN popen, plot_path+'sc'+sc_str+'_spectra'+region_str+'.ps', /land ELSE window, 5 
        zlim,proton_en_name+'_AVG'+at_str,1,1000.,1
        ylim,proton_en_name+'_AVG'+at_str,40,40000.,1
        ylim, [paspec_name_low+'_AVG'+at_str, paspec_name_high+'_AVG'+at_str],0,180,0
        zlim,[paspec_name_low+'_AVG'+at_str, paspec_name_high+'_AVG'+at_str],0.1,100,1
        ylim,'OMNI_HR_flow_speed',200,400,0
        ylim,'OMNI_HR_flow_pressure',1,20,1
        ylim,beta_name+'_AVG'+at_str,0.01,10,1
        ylim,p_total+'_AVG'+at_str,0.1,10,1
        options,'*','panel_size',1
        options,'*SPEC*','panel_size',2
        options,'Dst_Index','ytitle','Dst'
        options,'OMNI_HR_flow_pressure','ytitle','sw P'
        options,'OMNI_HR_flow_speed','ytitle','sw V'
        options,'Kp_Index','ytitle','Kp'
        options, ['PASPEC*'],'ztitle',''
        options,proton_en_name+'_AVG'+at_str,'ztitle',''

        if units_name_input eq 'DIST FUNC' then begin 
            zlim,'*SPEC*',1e-11,1e-5,1
        endif 

        timespan,time,time_duration,/hours
        tplot, [proton_en_name+'_AVG'+at_str,enspec_name+'_AVG'+at_str, paspec_name_low+'_AVG'+at_str, paspec_name_high+'_AVG'+at_str, $
                'OMNI_HR_flow_pressure','OMNI_HR_flow_speed','Dst_Index', 'Kp_Index',beta_name+'_AVG'+at_str,p_total+'_AVG'+at_str $
               ], var_label = var_label 
        tplot_panel, v = v_perp_name, o = v_perp_name_edi
        FOR ia = 0, na-1 DO BEGIN 
            index = sort(ABS(dd_avg.x - time_double(cut_time(ia)))) 
            if keyword_set(ps) then thickness=10 else thickness=1
            timebar,  dd_avg.x(index(0)),thick = thickness,$
              color = (ia+1)/2
        ENDFOR 
        
        IF keyword_set(ps) THEN pclose ELSE stop
    ENDIF 
    IF keyword_set(en_perp_plot) THEN BEGIN
        get_data, enspec_name+'_AVG'+at_str, data = dd_avg
        n_avg = n_elements(dd_avg.x)

        get_data, v_perp_name, data = data
        v_perp_tot_time = data.x
        v_perp_tot = data.y

        get_data, v_perp_name_edi, data = data
        v_perp_time_edi = data.x
        v_perp_edi = data.y

        v_perp_edi = INTERPOL(v_perp_edi, v_perp_time_edi, v_perp_tot_time )

        get_data, 'EPH_SC'+sc_str+'_GSE_X'+'_AVG'+at_str, data = data
        x_gse_time = data.x 
        x_gse = data.y

        v_perp_tot_cut = fltarr(na)
        v_perp_edi_cut = fltarr(na)

        IF keyword_set(ps) THEN popen, plot_path+'sc'+sc_str+'_ENperp_over_B.ps', /land  ELSE window, 6
        plot, [0, 0], [0, 0], xrange = [5, -20], yrange = [1e-3, 2e2], xlog = 0, ylog = 1, $
          charsize = 1.3, xstyle = 1, ystyle = 1, $
          position = [0.15, 0.1, 0.95, 0.92], $
          xtitle = 'X GSE', ytitle = 'Perp Energy /B', $
          title = 'CODIF sc'+sc_str+'   '+ts+'    (flux > 3) ', /nodata
        x = fltarr(na)
        FOR ia = 0, na-1 DO BEGIN 
            index = sort(ABS(dd_avg.x - time_double(cut_time(ia)))) 
            index_mag = sort(ABS(mag_time - time_double(cut_time(ia)))) 
            index_edi = sort(ABS(v_perp_time_edi - time_double(cut_time(ia)))) 
            index_x_gse = sort(ABS(x_gse_time -time_double(cut_time(ia))))
                                ;   v_para = fltarr(n_avg, 15)
            v_perp = fltarr(n_avg, 15)
            for ie = 0, n_elements(erange_set(0,*))-1 DO BEGIN 
                V_name = 'TDMOM_EN'+erange_set_str(2*ie)+'_'+erange_set_str(2*ie+1)$
                  +'_SC' + sc_str+'_MTVELOCITY_SP3_ET0_All'
                                ;       get_data, v_name+'_AVG'+at_str+'_V_PAR_T', data = data
                                ;        v_para(*, ie) = data.y
                get_data, v_name+'_AVG'+at_str+'_V_PERP_T', data = data
                v_perp(*, ie) = data.y
;            tplot_names, '*'+at_str+'_V*'
            ENDFOR 
            y = fltarr(15)

            for ie = 0, n_elements(erange_set(0,*))-1 DO BEGIN 
                x(ia) =  x_gse(index_x_gse(0)) ;replicate(ia+1, 15)
                
                IF dd_avg.y(index(0), ie) GT  3 THEN BEGIN 
                    y(ie) = mass_o*(v_perp(index(0), ie)*1e3)^2/2/1.6e-19/mag_avg(index_mag(0))
                ENDIF ELSE y(ie) = !values.f_nan
            ENDFOR         
            
            oplot, replicate(x(ia), 15), y, color = 1+ia/2, psym = (ia-2*fix(ia/2)+1), thick = 2
            v_perp_tot_cut(ia) = mass_o*(v_perp_tot(index(0))*1e3)^2/2/1.6e-19/mag_avg(index_mag(0))
            v_perp_edi_cut(ia) = mass_o*(v_perp_edi(index_edi(0))*1e3)^2/2/1.6e-19/mag_avg(index_mag(0))
            xyouts, 0.5, 100/(10^(ia/5.5)), time_string(dd_avg.x(index(0))), color = 1+ia/2, charsize = 1.2
;      xyouts, 8, 100/(10^(ia/5.5)), mag_avg(index_mag(0)), color = 1+ia/2, charsize = 1.2
            
        ENDFOR
        oplot, x, v_perp_tot_cut
        oplot, x, v_perp_edi_cut, color = 2
        IF keyword_set(ps) THEN pclose ELSE stop
    ENDIF 

    IF keyword_set(t_perp_plot) THEN BEGIN
        get_data, enspec_name+'_AVG'+at_str, data = dd_avg
        n_avg = n_elements(dd_avg.x)

        get_data, t_name_tot+'_AVG'+at_str, data = data
        t_perp_tot_time = data.x
        t_perp_tot = (data.y(*, 1)+data.y(*, 2))/2.

        get_data, 'EPH_SC'+sc_str+'_GSE_X'+'_AVG'+at_str, data = data
        x_gse_time = data.x 
        x_gse = data.y

        t_perp_tot_cut = fltarr(na)
        
        IF keyword_set(ps) THEN popen, plot_path+'sc'+sc_str+'_Tperp.ps', /land  ELSE window, 6
        plot, [0, 0], [0, 0], xrange = [30, 40000], yrange = [0.01, 20000.], xlog = 1, ylog = 1, $
          charsize = 1.3, xstyle = 1, ystyle = 1, $
          position = [0.15, 0.1, 0.95, 0.92], $
          xtitle = 'TOTAL ENERGY (eV)', ytitle = 'T perp (eV)', $
          title = 'CODIF sc'+sc_str+'   '+ts+'    (flux > 3) ', /nodata
        x = fltarr(na, 15)
        x_tip = fltarr(na) 
        FOR ia = 0, na-1 DO BEGIN 
            index = sort(ABS(dd_avg.x - time_double(cut_time(ia)))) 
            index_mag = sort(ABS(mag_time - time_double(cut_time(ia)))) 
            index_x_gse = sort(ABS(x_gse_time -time_double(cut_time(ia))))
            
            t_perp = fltarr(n_avg, 15)
            for ie = 0, n_elements(erange_set(0,*))-1 DO BEGIN 
                T_name = 'TDMOM_EN'+erange_set_str(2*(ie+1 < 14))+'_'+erange_set_str(2*(ie-1 > 0)+1)$
                  +'_SC' + sc_str+'_MTTEMPERATURE_SP3_ET0_All'
                get_data, t_name+'_AVG'+at_str, data = data
                t_perp(*, ie) = (data.y(*, 1)+data.y(*, 2))/2
            ENDFOR 
            y = fltarr(15)

            for ie = 0, n_elements(erange_set(0,*))-1 DO BEGIN 
                x(ia, ie) = dd_avg.v(index(0), ie)
                
                IF dd_avg.y(index(0), ie) EQ max(dd_avg.y(index(0), *)) OR ia EQ 0 THEN BEGIN 
                    y(ie) = t_perp(index(0), ie)
                    x_tip(ia) = x(ia, ie)
                ENDIF ELSE y(ie) = !values.f_nan
            ENDFOR
            
            oplot, x(ia, *) $   ;replicate(x(ia), 15)$
              , y, color = 1+ia/2, psym = -(ia-2*fix(ia/2)+1), thick = 2
            t_perp_tot_cut(ia) = t_perp_tot(index(0))
            xyouts, 4e3, 100/(10^(ia/5.5)), time_string(dd_avg.x(index(0))), color = 1+ia/2, charsize = 1.2
        ENDFOR

        oplot, x_tip(1:11), t_perp_tot_cut(1:11), psym = 6
        
        IF keyword_set(ps) THEN pclose ELSE stop
    ENDIF  

    IF keyword_set(pa_plot) THEN BEGIN 
        pa_names = strarr(n_elements(erange_set(0, *)))
        FOR  ie = 5, n_elements(erange_set(0, *))-1 DO BEGIN 
            pa_names(ie) = $
              'PASPEC_EN'+ erange_set_str(2*ie)+'_'+erange_set_str(2*ie+1)$
              +'_SC'+sc_str+'_UN'+STRUPCASE(strcompress(units_name_input,/remove_all))+'_SP3_All'+'_AVG'+at_str
            options, pa_names(ie), 'ytitle', 'SC'+sc_str+' O+!C'+erange_set_str(2*ie)+'!C'+erange_set_str(2*ie+1)+'!CPA (deg)'
        ENDFOR 
        options, '*', 'panel_size', 1
        zlim, pa_names, 0.1, 1000, 1
        zlim, [enspec_name+'_AVG'+at_str], 0.1, 100, 1

        ylim, pa_names, 180, 0
                                ;   tplot, [enspec_name, pa_names]
        FOR icut = 0, n_elements(cut_time)-1 DO BEGIN 
            cal_time = time_double(cut_time(icut))
            timespan, cal_time-150, 300, /seconds
            te = time_string(cal_time)
            time_str = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2) $
              +'_'+ STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2)
            
            get_data, enspec_name+'_AVG'+at_str, data = dd_avg
            ind_avg = sort(ABS(dd_avg.x-cal_time))
            
            IF keyword_set(ps) THEN popen, plot_path+time_str+'_pa_spectra.ps', /port 
            
            tplot, [enspec_name+'_AVG'+at_str, pa_names], var_label = var_label
            timebar, dd_avg.x(ind_avg(0))-avg_time/2
            timebar, dd_avg.x(ind_avg(0))+avg_time/2
            IF keyword_set(ps) THEN pclose ELSE stop
        ENDFOR  
        timespan, ts, time_duration, /hours 
    ENDIF  

    tplot_names,names=names
    store_data,delete = names
endfor   
END 
