; This program is used to identify cusp pass of cluster with O+ perp
; energy and others. 

pro cusp_identification_original, cusp_storm=cusp_storm,ps_plot=ps_plot,phase_set=phase_set
;-----------------------------------------
;Set the keywords
;-----------------------------------------
; it can plot all passes close the the perigees for storm or nonstorm
; time or it can also plot only those identifiecd cusp passes for
; storm time (I only did the storm time cusps identification).
if not keyword_set(phase_set) then phase_set=['storm'];'nonstorm','all']
if keyword_set(ps_plot) and not keyword_set(phase_set) and not keyword_set(cusp_storm) then phase_set=['all']
if keyword_set(cusp_storm) and not keyword_set(phase_set) then phase_set=['storm']

sc_set = [1,3,4]                ;set the satallite number 
for isc=0,2 do begin
    sc= sc_set(isc)
    sc_str = STRING(sc, FORMAT = '(i1.1)')
    ion=3                       ;0 is proton, 3 is O+
    ion_str=strcompress(ion,/remove_all)
    plot_time =  7 * 3600.      ;in sec
    average_time = 1 * 60       ;in seconds  
    var_label = 'EPH_SC' + sc_str + '_'
    var_label = var_label + ['MLT', 'GSM_X', 'GSM_Y', 'GSM_Z', 'DIST']
    diffflux_threshold=0
    plot_path='output/cusp/'
    spawn,'mkdir '+plot_path
    spawn,'mkdir '+plot_path+'nonstorm/'
    spawn,'mkdir '+plot_path+'storm/'
    if not keyword_set(cusp_storm) then begin 
;---------------------------------------------------------------------
; Read CLUSTER perigee times list from file: sc4_perigee_times.dat
; (different list for each S/C)
; Store these times in variable petime (in tplot time format)
;-------------------------------------------------------------------
        OPENR, unit, 'sc' + sc_str + '_perigee_times.dat', /GET_LUN
        petime = DBLARR(3000) ; assumes no more than 3000 perigee passes
        dummy = ''
        ii = 0l
        WHILE NOT EOF(unit) DO BEGIN
            READF, unit, dummy
        petime(ii) = time_double(dummy)
        ii = ii + 1
    ENDWHILE                            ;2002-05-12
    index=where(petime gt time_double('2002-09-10') and petime lt time_double('2002-12-31'),ct)
    if ct gt 0 then petime = petime(index) ;all perigee passes time
    n_petime=n_elements(petime) ; number of perigee passes
    CLOSE, unit
    ntimes=n_petime
    nj=1
endif else begin   
;---------------------------------------------------------------------
; Read cusp from list
;-------------------------------------------------------------------
    OPENR, unit, 'cusp_list_sc'+sc_str+'.dat', /GET_LUN
    cusp_data = DBLARR(1000,6)   ; assumes no more than 1000 perigee passes
    dummy = ''
    ii = 0l & jj=0l
    WHILE NOT EOF(unit) DO BEGIN
        READF, unit, dummy
        IF STRMID(dummy, 0, 3) EQ  '200' THEN BEGIN
            cusp_data(jj,0) = time_double(STRMID(dummy, 0, 19))
            cusp_data(jj,1) = time_double(STRMID(dummy, 20, 19))
            READS,STRMID(dummy,39,2),a
            if a ne sc then begin
                print,'sc is not correct' & stop
            endif 
            cusp_data(jj,2) = a
            READS,STRMID(dummy,41,2),a
            cusp_data(jj,3) = a
            cusp_data(jj,4) = time_double(STRMID(dummy,44,19))
            cusp_data(jj,5) = time_double(STRMID(dummy,64,19))
            jj=jj+1
        ENDIF   
        ii = ii + 1
    ENDWHILE
    n_cusp = jj-1 ; number of perigee passes
    CLOSE, unit 
    ntimes=n_cusp 
    nj=0
endelse
;----------------------------------------------------
;Write [START] in log files
;-----------------------------------------------------
;OPENU, unit, path+'cusp_plotted.txt', /GET_LUN, /APPEND
;PRINTF, unit, SYSTIME(), '[START]'
;FREE_LUN, unit
    FOR ii = 0, ntimes-1  DO BEGIN ; run over all 
        for jj =0,nj do begin    ; for north and south cusps passes
;-------------------------------------------------------------
;Delete all the string stored data in
;order to make sure the program can run correctly
;-----------------------------------------------------------
            tplot_names, names = names  & store_data, DELETE = names
            PRINT, STRING(ii) + '   --' + STRING(jj)
; Timespan over each displaytime
            if keyword_set(cusp_storm) then begin 
                t_s=cusp_data(ii,4) &  t_e=cusp_data(ii,5) 
            endif else begin 
                t_s = petime(ii) - (plot_time)*(1-jj)  & t_e = t_s + plot_time
            endelse 
            timespan, t_s, t_e-t_s, /SECONDS 
                                ;      timespan, '2001-09-28/5:40',50,/min
            ts = time_string(t_s)  
            te = time_string(t_e)
            date_s = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
            time_s = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
            date_e = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2)
            time_e = STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2)
;-- Load CLUSTER ephemeris--
            sat = [sc] 
            get_cluster_ephemeris, sat, /GSE_X, /GSE_Y, /GSE_Z, /DIST, /MLT, /GSM_X, /GSM_Y, /GSM_Z, /ILAT_D
;-- Load CLUSTER O+ pitch angle spectrum for energy > 1keV
;            sat = [sc] &  specie=[ion] & energy=[1000,40000.]
;            inst = 0 & eff_table = 0  &  units_name = 'DIFF FLUX'
;            plot_pa_spec_from_crib, sat, specie, inst, units_name,  energy, eff_table,  recalc = 1, COMBINE = 1
;-- Load CLUSTER O+ pitch angle spectrum for energy < 1keV
;            sat = [sc] &  specie=[ion] & energy=[40,1000.]
;            inst = 0 & eff_table = 0  &  units_name = 'DIFF FLUX'
;            plot_pa_spec_from_crib, sat, specie, inst, units_name,  energy, eff_table,  recalc = 1, COMBINE = 1
;-- Load energy spectra -
            sat = [sc,sc]
            specie = [0,ion]
            angle = [[-90, 90], [0, 360]]
            inst = 0 & units_name = 'DIFF FLUX'
            eff_table = 0
            plot_en_spec_from_crib, sat, specie, inst,  units_name, angle, eff_table, recalc = 1
            tplot_names,'EN*DIFFFLUX*',names=names
            get_data,names(0),data=dd,dlim=dlim,lim=lim
            if names(0) ne '' then begin 
                time=dd.x  & ntime=n_elements(time)
                ddy=dd.y
                index=where(ddy lt diffflux_threshold,ct)
                if ct gt 0 then  ddy(index)=!values.f_nan
                store_data,names(0),data={x:dd.x,y:ddy,v:dd.v},dlim=dlim,lim=lim
;-- Load CLUSTER Magnetic field--
                sat = sc
                plot_mag_from_crib, sat, POLAR = 1, GSM = 1
;-- Load CLUSTER H+ and O+ moments--
                sat = [sc]
                specie = [ion]
                moments = ['A']
                angle = [[-90, 90], [0, 360]]
                energy = [40., 40000.]
                inst = 0 &  eff_table = 0      
                plot_3dmom_from_crib, sat, specie, inst, moments, angle, energy, eff_table, recalc = 1,diffflux_threshold=diffflux_threshold
;-----------------------------------------------------
; read the storm phase and flag the cusp into different storm phases
;-------------------------------------------------------
                OPENR, unit_temp, 'storm_phase_long.dat', /GET_LUN
                prestorm_start = DBLARR(300)
                storm_onset = DBLARR(300) 
                min_dst_new = DBLARR(300)
                recovery_fast_end = DBLARR(300) 
                recovery_early_end = DBLARR(300) 
                recovery_long_end = DBLARR(300)
                jj_storm = 0l
                dummy = ''
                WHILE NOT EOF(unit_temp) DO BEGIN
                    READF, unit_temp, dummy
                    IF STRMID(dummy, 0, 5) NE 'Start' THEN BEGIN 
                        prestorm_start(jj_storm) = time_double(STRMID(dummy, 0, 20))
                        storm_onset(jj_storm) = time_double(STRMID(dummy, 20, 20))
                        min_dst_new(jj_storm) =  time_double(STRMID(dummy, 40, 24))
                        recovery_fast_end(jj_storm) =  time_double(STRMID(dummy, 60, 24))
                        recovery_early_end(jj_storm) =  time_double(STRMID(dummy, 80, 24))
                        recovery_long_end(jj_storm) = time_double(STRMID(dummy, 100, 24))        
                        jj_storm = jj_storm+1      
                    ENDIF  
                ENDWHILE
                CLOSE, unit_temp, /all
                nstorm = jj_storm
                prestorm_start = prestorm_start(0:nstorm-1)
                storm_onset = storm_onset(0:nstorm-1)
                min_dst_new = min_dst_new(0:nstorm-1)
                recovery_fast_end = recovery_fast_end (0:nstorm-1)
                recovery_early_end =  recovery_early_end (0:nstorm-1)
                recovery_long_end =  recovery_long_end (0:nstorm-1)
                storm_phase = INTARR(ntime)
                FOR  itime = 0, ntime-1 DO BEGIN
                    FOR istorm = 0, nstorm-1 DO BEGIN 
                        belong = 0
                        IF time(itime) GE prestorm_start(istorm) AND $
                          time(itime) LT storm_onset(istorm) THEN BEGIN 
                            storm_phase(itime) = 1 ;initial phase
                            belong = belong+1
                        ENDIF 
                        IF time(itime) GE storm_onset(istorm) AND $
                          time(itime)  LT min_dst_new(istorm) THEN BEGIN 
                            storm_phase(itime) = 2 ; mian phase
                            belong = belong+1
                        ENDIF 
                        IF time(itime) GE min_dst_new(istorm) AND $
                          time(itime) LT recovery_early_end(istorm) THEN BEGIN 
                            storm_phase(itime) = 3 ; recovery phase
                            belong = belong+1
                        ENDIF
                        IF time(itime) GE recovery_early_end(istorm) AND $
                          time(itime) LT recovery_long_end(istorm) THEN BEGIN 
                            storm_phase(itime) = 4 ; later recovery phase (not used for any map)
                            belong = belong+1
                        ENDIF
                        IF belong GT 1 THEN stop
                        if belong eq 0 then begin 
                            IF time(itime) GE (prestorm_start(istorm)-60.*60.*4) AND $
                              time(itime) LT prestorm_start(istorm) THEN BEGIN 
                                storm_phase(itime) = 5 ; pre storm
                                belong = belong+1
                            endif 
                        endif  
                    ENDFOR 
                ENDFOR   
                store_data, 'storm_phase', data = {x:time, y:storm_phase}
                ylim, 'storm_phase', -1, 6 

                V_name = 'TDMOM_EN'+string(min(energy),format='(i5.5)')+'_'+string(max(energy),format='(i5.5)')+ '_SC' + sc_str + '_MTVELOCITY_SP'+ion_str+'_ET0_All'
                tplot_names,v_name,names=names
                if names ne '' then v_perp,v_name

                options,'*','panel_size',1
                options,'PA*','ztitle',''
                ylim,'*TEM*SP'+ion_str+'*T',10,1e4
                
                zlim,'EN*DIFFFLUX*SP0*',1,1e3,1
                zlim,'EN*EFLUX*SP3*',1e4,1e6
                zlim,'EN*DIFFFLUX*SP3*',1,1e2
                zlim,'PA*',0.1,100,1
                zlim,'*COUNTS*',1,100,1
                ylim,'*VEL*SP'+ion_str+'*All_T',0,200
                tplot_names,'*VEL*SP'+ion_str+'*PERP_T',names=names
                options,names(0),'ytitle','Vperp'
                get_data,names(0), data=vperp
                tplot_names,'*VEL*SP'+ion_str+'*PAR_T',names=names
                ylim,names(0),-200,200
                options,names(0),'ytitle','Vpar'
                get_data,names(0),data=vpar
                store_data,names(0),data={x:vpar.x,y:-vpar.y}
                tplot_names,'*TEM*SP'+ion_str+'*X',names=names
                get_data,names(0),data=tpar
                tplot_names,'*TEM*SP'+ion_str+'*Y',names=names
                get_data,names(0),data=tperp1
                tplot_names,'*TEM*SP'+ion_str+'*Z',names=names
                get_data,names(0),data=tperp2
                mass = 16*1.6e-27*(1e3)^2/(1.6e-19) ; unit: ev/(km/s)^2
                Epar=mass*(vpar.y)^2/2+tpar.y
                store_data,'Epar',data={x:vpar.x,y:Epar}
                Eperp=mass*(Vperp.y)^2/2+sqrt(tperp1.y^2+tperp2.y^2)
                store_data,'Eperp',data={x:Vperp.x,y:Eperp}
                
                store_data,'apex',data={x:vpar.x,y:atan(ABS(Eperp/Epar))*180./3.14}
                ylim,'apex',0,90,0
                ylim,'Ep*',1,1e4,1
                index=where(storm_phase gt 1 and storm_phase le 3,ct)
                if ct gt 0 then storm_flag=1 else storm_flag=0
                if phase_set ne 'storm' or storm_flag ge 0 then begin
                    if keyword_set(ps_plot) then begin 
                        if keyword_set(cusp_storm) then popen, plot_path+'cusp_storm/cusppass_'+date_s+'_'+time_s+'_sc'+sc_str+'.ps',/land else begin 
                            if  storm_flag eq 1 then popen, plot_path+'perigee_storm/perigee_pass_'+date_s+'_'+time_s+'_sc'+sc_str+'.ps',/land else popen,plot_path+'perigee_nonstorm/perigee_pass_'+date_s+'_'+time_s+'_sc'+sc_str+'.ps',/land
                        endelse 
                    endif
                    tplot, ['storm_phase','EN*','*DEN*SP'+ion_str+'*','*PAR_T','Eperp','Epar','*TEM*SP'+ion_str+'*T','apex'], var_label = var_label
                    
                    if keyword_set(cusp_storm) then begin 
                        timebar,cusp_data(*,0)
                        timebar,cusp_data(*,1)
                    endif 
                    if keyword_set(ps_plot) then pclose else stop
                endif   
            endif
        endfor     
    ENDFOR  
endfor  
stop 
end
