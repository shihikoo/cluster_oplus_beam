PRO  sim_data

show_fit = 0
cal_mom = 0
dfit_temperature = 0
globe_plot = 0
single_plot = 0
cal_eflux = 1

plot_ps = 1
sc = 4
sc_str = STRING(sc, FORMAT = '(i1.1)')

inst_input = 0
path = 'output/o_beam/sim/'
spawn, 'mkdir '+path
average_time = 300
at_str = strcompress(average_time, /remove_all)
units_name = 'DIFF FLUX'
eflux_all = fltarr(5, 12,16)
FOR iy = 2001, 2005 DO BEGIN 
    year = string(iy, format = '(i4.0)')
 ;   tplot_restore, file = 'sample.tplot'
  ;  get_data, 'sample', data = dat
    time = year+'-10-01/06:00:00'
    IF year EQ '2005' THEN time = year+'-09-20/06:00:00'
    IF year EQ '2001' THEN time = year+'-08-17/06:00:00'
    timespan, time, average_time, /SEC   ; SECONDS, MINUTES, HOURS, DAYS (DEFAULT)
    units_name = 'Counts'

    sat   = 4
    specie = 3
    inst = 0                    ; 0: CODIF, 1: HIA (this is not supported for the moment)
    eff_table = 0               ; 0: GROUND, 1: ONBOARD
    plot_globe_from_crib, sat, specie, inst, units_name, BKG = 0,  eff_table, OLD_EFF = 0, CNES = 0, IC = 0
    name = 'GLOBE_SC'+string(sat, format = '(i1.1)')+$
           '_'+strcompress(units_name, /remove_all)  +$
           '*'+'SP'+string(specie, format = '(i1.1)')
    tplot_names, name, names = gname
    get_data, gname(0), data = dat

FOR ie = 0, 0 DO BEGIN 
    FOR ia = 0, 0 DO BEGIN 
        ie_str = strcompress(ie+1, /remove_all)
        ia_str = strcompress(ia+1, /remove_all)
        IF keyword_set(single_plot) THEN BEGIN 
            plot_path = paht+ie_str+'ebin_'+ia_str+'bins/'
            spawn, 'mkdir '+plot_path
        ENDIF 
;------sim data ------------------------
    
        eflux = fltarr(16, 88)

        pa = fltarr(88)
        abin = fltarr(88)
        tem_ea = fltarr(16, 88)
        dt = fltarr(16, 88)
        tem_ea(*) =  !VALUES.F_NAN
        abin(*) = !VALUES.F_NAN
        pa(*) = !VALUES.F_NAN
        dt(*) = !VALUES.F_NAN
        FOR ii = 0, 15-ie DO BEGIN 
            FOR jj = 0, 87-ia DO BEGIN                 
                ntime = 1
                
                data_t = FLTARR(ntime, 3) &  time_t = FLTARR(ntime)
                data = fltarr(16, 88)
                input_data = fltarr(ntime)
                t_dfit =  FLTARR(ntime, 2)   &  t_error =  FLTARR(ntime, 2)
                d_dfit = FLTARR(ntime)    &  d_error = FLTARR(ntime)

                FOR i_m = 0, ntime-1 DO BEGIN
                    data(ii:(ii+ie), jj:(jj+ia)) = 9
                    str_element, dat, 'data', data, add_replace = 1   
                    
                    pa(jj) = 90-Total(dat.theta(ii, jj:(jj+ia)))/(ia+1)
                    abin(jj) = jj+ia/2.
                    
                    input_data(i_m) = i_m*10+1
;------- globe plot --------------
                    IF keyword_set(globe_plot)  THEN BEGIN 
                        
                        IF keyword_set(plot_ps) $
                          THEN popen, path+'counts/data'+strcompress(ii, /remove_all) $
                                      +'_'+strcompress(jj, /remove_all) $
                                      +'_eq_counts_globe.ps', /land
                        plot3d_options, log = 1
                        plot3d_codif, dat, zrange = [1, 100]
                        IF  keyword_set(plot_ps)  THEN pclose 
                    ENDIF
;------- dif temerature ---------------
                    IF keyword_set(dfit_temperature) THEN BEGIN 
                        sat = sc
                        specie = 3
                        en_range = [40., 40000.]
                        IF NOT keyword_set(show_fit) THEN path_input = path
                        t_calc_from_disf, sat = sat, specie = specie, $
                                          energy = en_range, $
                                          t_dfit_name = t_dfit_name, path = path_input
                        
                        get_data, t_dfit_name, data = data_fit
                        t_dfit(i_m, *) = [data_fit.t(0), data_fit.t(1)]
                        t_error(i_m, *) = [data_fit.t_error(0), data_fit.t_error(1)]
                        d_dfit(i_m) = [data_fit.n]
                        d_error(i_m) = [data_fit.n_error]
                        store_data, delete = t_dfit_name
                    ENDIF 
;-------- cal mom ----------------------  
                    IF keyword_set(cal_mom) THEN BEGIN 
                        sim_name = 'sim_tem'
                        eff_table = 0
                        energy = [40., 40000.]
                        sat = sc
                        dat_ef = convert_codif_units(dat, 'EFLUX', 'codif_ts_eff_corr', eff_table, $
                                                     packets = n_elements(dat.time), sat = sat)
                        angle = [[-90.0, 90.0], [0., 360.]] 
                        temperature = compute_temperature(dat_ef, sat, NAME = sim_name, $
                                                          ENERGY = energy, ANGLE = angle) 
                        get_data, sim_name, data = tem
                        data_t(i_m, *) = tem.y
                        time_t(i_m) = tem.x
                    ENDIF
;-------- cal Eflux -----------
                    IF keyword_set(cal_eflux) THEN BEGIN 
                        sat = sc  & specie = 3
                        eff_table = 0 &   inst = 0
                        units_name = 'EFLUX'
                        angle = [[-90, 90], [90, 270]]
  
                        dat_ef = convert_codif_units(dat, 'EFLUX', 'codif_ts_eff_corr', eff_table, $
                                                     packets = n_elements(dat.time), sat = sat)
                        ; divided by 88 for normalization for 'eflux'
                        eflux(ii, jj) = dat_ef.data(ii:(ii+ie), jj:(jj+ia))/88.
                        
                        IF keyword_set(globe_plot) THEN  BEGIN 
                            IF keyword_set(plot_ps) THEN BEGIN 
                                spawn, 'mkdir '+path+'eflux/'
                                popen, path+'eflux/'+year+'data_'+strcompress(ii, /remove_all) $
                                       +'_'+strcompress(jj, /remove_all) $
                                       +'_eflux_globe.ps', /land
                            ENDIF 
                            plot3d_options, log = 1
                            plot3d_codif, dat_ef, zrange = [1, 1e3], $
                                          title = year+'EFLUX: '+strcompress(eflux(ii, jj), /remove_all)
                            IF  keyword_set(plot_ps)  THEN pclose 
                        ENDIF 
                    ENDIF   
                ENDFOR 
                IF keyword_set(cal_mom) THEN BEGIN 
                    temp = (total(data_t, 2, /nan))/3.
                    tem_ea(ii, jj) = Total(temp)/N_elements(temp)
                    dt(ii, jj) = sqrt(total((temp-tem_ea(ii, jj))^2)/(n_elements(temP)-1))
                    
                    IF keyword_set(single_plot) THEN BEGIN 
                        IF keyword_set(plot_ps)  THEN $
                          popen, plot_path+year+'data' +strcompress(ii, /remove_all) +'_' $
                                 +strcompress(jj, /remove_all)+'_eq_counts_tem.ps', /land $
                        ELSE window, 1
                        store_data, ''
                        plot, input_data, temp, psym = 1, ystyle = 1, $
                              xtitle = 'Counts', ytitle = 'T (eV)'
                        IF keyword_set(plot_ps) THEN pclose
                    ENDIF 
                ENDIF
            ENDFOR 
        ENDFOR
        IF keyword_set(cal_eflux) THEN BEGIN 
                                ;get the eflux data 
            store_data, 'dat_ef'+year, data = {x:dat_ef.time, y:eflux}
            eflux_avg = total(eflux, 2)/n_elements(eflux(0, *))
                                ; minimum eflux
            eflux_min = fltarr(16)
            FOR ief = 0, 15 DO eflux_min(ief) = min(eflux(ief, *))
                                ; all 8 different eflux
            eflux_8 = fltarr(8, 16)
            eflux_8(0, *) = total(eflux(*, 0:3), 2)/4
            eflux_8(1, *) = total(eflux(*, 4:11), 2)/8
            eflux_8(2, *) = total(eflux(*, 12:27), 2)/16
            eflux_8(3, *) = total(eflux(*, 28:43), 2)/16
            eflux_8(4, *) = total(eflux(*, 44:59), 2)/16
            eflux_8(5, *) = total(eflux(*, 60:75), 2)/16
            eflux_8(6, *) = total(eflux(*, 76:83), 2)/8
            eflux_8(7, *) = total(eflux(*, 84:87), 2)/4
                                ; find the err of average eflux
            err_eflux_avg = fltarr(16)
            FOR iebin = 0, 15 DO err_eflux_avg(iebin) = sqrt(total((eflux(iebin, *)-eflux_avg(iebin))^2)/n_elements(x))
            
            store_data, 'eflux'+year, data = {x:dat.energy(*, 0), y:eflux_avg, dy:err_eflux_avg}
            eflux_all(iy-2001, 0, *) = dat.energy(*, 0)
            eflux_all(iy-2001, 1, *) = eflux_avg
            eflux_all(iy-2001, 2, *) = err_eflux_avg
            eflux_all(iy-2001, 3, *) = eflux_min
            eflux_all(iy-2001, 4:11, *) = eflux_8
        ENDIF 

        IF keyword_set(plot_mom) THEN BEGIN 
            IF keyword_set(plot_ps) THEN    popen, 'sim/'+ie_str+'a'+ia_str+'_T.ps', /land $
            ELSE window, 1
            specplot, dat.energy(*, 0), abin, $
                      (tem_ea), $
                      no_interp = 1, $
                      lim = { zlog:0, xlog:1, $
                              title: 'T/dE for '+ie_str+'ebins'+ia_str+'bins', $
                              xtitle: 'Energy (eV)', $
                              ytitle:'Angle bins', $
                              xrange:[30, 40000], yrange: [0, 87], $
                              XSTYLE:1, ystyle: 1, charsize: 1.2, $
                              position: [0.1, 0.1, 0.82, 0.9]}  
            IF keyword_set(plot_ps)  THEN    pclose 
            IF keyword_set(plot_ps)  THEN   popen, 'sim/'+ie_str+'a'+ia_str+'dT.ps', /land $ 
            ELSE window, 2

            specplot, dat.energy(*, 0), abin, $
                      dT, $
                      no_interp = 1, $
                      lim = { zlog:0, xlog:1, $
                              title: 'dT (eV) for '+ie_str+'ebins'+ia_str+'bins', $
                              xtitle: 'Energy (eV)', $
                              ytitle:'Angle bins', $
                              xrange:[30, 40000], yrange: [0, 87], $
                              XSTYLE:1, ystyle: 1, charsize: 1.2, $
                              position: [0.1, 0.1, 0.82, 0.9]}  
            IF keyword_set(plot_ps)  THEN   pclose
            
; if negative T and P exist, plot and record in log 
            index = where (data_t LT 0, ct)
            IF ct GT 0 THEN BEGIN                   
                OPENU, unit, 'sim/log_neg_T.txt', /GET_LUN, /APPEND
                PRINTF, unit,  '---negative Tt and Pt---'+'ebin:' $
                        +ie_str+'abin:'+ia_str
                FREE_LUN, unit  
            ENDIF 
        ENDIF 
        
        IF keyword_set(dfit_temperature) THEN BEGIN 
            Tt_fit_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                          +phi_str+'_MTTEMPERATURE_SP3_ET0_All_para_dfit'+'_AVG'+at_str
            store_data,  tt_fit_name, data = {x:time_tt, $
                                              y:t_dfit(*, 0), dy:t_error(*, 0)}, $
                         dlim = dlim_Tt, lim = lim_Tt
            options, Tt_fit_name, 'ytitle', 'SC'+sc_str+' O!U+!N!C!CT!D//!N (eV)'
            
            d_fit_name = 'TDMOM_ENVARIOUS'+ '_SC' + sc_str+'_' $
                         +phi_str+'_MTDENSITY_SP3_ET0_All_dfit'+'_AVG'+at_str
            store_data,  d_fit_name, data = {x:time_d, $
                                             y:d_dfit, dy:d_error}, $
                         dlim = dlim_d, lim = lim_d
            options, d_fit_name, 'ytitle',  'SC'+sc_str+'!C!CO!U+!N!C!Cn (cm!U-3!N)' 
        ENDIF 
    ENDFOR 
ENDFOR 
ENDFOR

store_data, 'eflux_all', data = eflux_all
ps = 1

FOR iy = 0, 4 DO BEGIN 
    IF keyword_set(ps) THEN popen, $
      path+strcompress(iy+2001, /remove_all) $
      +'_9counts_eflux_vs_energy_avg'+at_str+'.ps', /land
    
    plot, [40, 40000], [6e5, 6e5], xstyle = 1, yrange = [20, 20000], $
          ystyle = 1, xlog = 1, ylog = 1, xrange = [40, 40000], /nodata, $
          title = strcompress(iy+2001, /remove_all) $
          +',   9 counts, average time: '+ at_str, $
          xtitle = 'Energy  (eV)', ytitle = 'EFLUX    (eV/cm!E2!N-s-sr-eV)'
    FOR ia = 4, 11 DO oplot, eflux_all(iy, 0, *), $
      eflux_all(iy, ia, *), color = abs(ia-7.5)+0.5, thick = 2
    FOR ia = 4, 11 DO xyouts, ia*400-600, 5000, ia-3, color = abs(ia-7.5)+0.5, charsize = 2
    IF keyword_set(ps) THEN pclose
ENDFOR 

;IF keyword_set(ps) THEN popen, path+'9counts_eflux_vs_energy_avg'+at_str+'.ps', /land
;plot, [40, 40000], [6e5, 6e5], xstyle = 1, yrange = [10, 1e4], ystyle = 1, xlog = 1, ylog = 1, xrange = [40, 40000], /nodata, title = '9 counts, average time: '+ at_str, xtitle = 'Energy  (eV)', ytitle = 'EFLUX    (eV/cm!E2!N-s-sr-eV)'
;average
;FOR iy = 0, 4 DO BEGIN
 ;   oplot, eflux_all(iy, 0, *), eflux_all(iy, 1, *), color = iy+2, thick = 2
  ;  xyouts, 1000, 10*(iy*0.5+1.2), iy+2001, color = iy+2
;ENDFOR 
;min
;FOR iy = 0, 4 DO  oplot, eflux_all(iy, 0, *), eflux_all(iy, 3, *), color = iy+2, linestyle = 2, thick = 2
;IF keyword_set(ps) THEN pclose

stop
END 
