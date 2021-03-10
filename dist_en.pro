pro dist_en,output,func_plot=func_plot,unit=unit,ps=ps,spec_plot=spec_plot,other_inputs=other_inputs,direct_cusp=direct_cusp,storm_phase=storm_phase,vperp=vperp
inst_input=0
direct_cusp=1
if not keyword_set(unit) then unit = 2
if unit eq 1 then unit_name='DIFF FLUX'
if unit eq 2 then unit_name='DIST FUNC'
if unit eq 3 then unit_name='FLUX'
pi=3.1415926
mass = 16*1.6e-27*(1e3)^2/(1.6e-19) ; unit: ev/(km/s)^2
en_set = [31444.7, 19398.3, 11966.9, 7382.39, 4554.22, 2809.51, 1733.19, 1069.21, $
          659.599, 406.909, 251.023, 154.857, 95.5315, 58.9337, 36.3563]
den_set = [14900.6, 9192.19,5670.69,3498.26,2158.09,1331.33,821.301,506.663,$
           312.562, 192.820,118.951,73.3813,45.2691,27.9267,17.2280,5.63791]
en_range = [[en_set-den_set/2],[en_set+den_set/2]]
Bt= 31200.*(6370./(6370+1000))^3*sqrt(1+3*sin(80*3.1415926/180)^2) ; =39832.1 ; dipole field at 1000km and 80 invariant latitude ;33695.9 ;nT at 60 invariant latitude degree from Seki 1998; 
;if keyword_set(func_plot) then !p.multi=[0,3,1]

if not keyword_set(other_inputs) then begin 
;--- known titles info in data files ---
    title_input =  [ 'DIST ', $ ;0
                     'DENS ', $ ;1
                     'Vpar ', $ ;2
                     'Tpar ', $ ;3
                     'Tperp ', $ ;4
                     ' T ', $   ;5
                     'Emean ', $ ;6
                     ' Apex ', $ ;7
                     '  sc  ' $ ;8
                   ]
; Setup data arrays 
    title_lenth = STRLEN(title_input)
    nterms =  N_ELEMENTS(title_input)
    ntime = 300
    time = DBLARR(ntime)
    data = DBLARR(ntime, nterms)
;---------------------------------------------------------------
; Read data into time and data arraies
;---------------------------------------------------------------
; Read all 1 day files that correspond to requested time interval
    jj = 0l                  
    fln = 'cusp_list_bouhram.dat'
    names = FINDFILE(fln)
    IF names(0) NE '' THEN BEGIN 
        OPENR, unit_temp, names(0), /GET_LUN
        dummy = ''
        a = DBLARR(nterms)
        WHILE NOT EOF(unit_temp) DO BEGIN
            READF, unit_temp, dummy
            IF STRMID(dummy, 0, 3) EQ  '200' THEN BEGIN
                year = STRMID(dummy, 0, 4)
                doy = STRMID(dummy, 5, 4)
                caldat, julday(01,00,year) + doy, month, day
                HH=STRMID(dummy,10,2)
                MM=STRMID(dummy,13,2)
                SS=STRMID(dummy,16,2)
                time(jj) = time_double(string(year, format = '(i4.4)')+'-'+string(month, format = '(i2.2)')+'-'+string(day, format = '(i2.2)') $
                                       +'/'+string(HH, format = '(i2.2)')+':'+string(MM, format = '(i2.2)')+':'+string(SS, format = '(i2.2)'))
                READS, STRMID(dummy,21), a 
                data(jj, *) = a
                jj = jj+1
            ENDIF   
        ENDWHILE
        CLOSE, unit_temp, /all
    ENDIF     
    ntime = jj
    IF ntime EQ 0 THEN BEGIN 
        print, 'no files found' & stop
    ENDIF 
    time = time(0:ntime-1)
    data = data(0:ntime-1, *)
ENDIF
density=data(*,1)
Vpara_input=data(*,2)
Tem=data(*,5)
Emean=data(*,6)
sc=data(*,8)

ievent = 0
event_flag = DBLARR(ntime)
for itime=1,ntime-1 do begin 
    if (time(itime)-time(itime-1)) ge 0 and (time(itime)-time(itime-1)) le 24*3600. $
      and sc(itime) eq sc(itime-1) then event_flag(itime)=ievent else begin 
        ievent=ievent+1
        event_flag(itime)=ievent
    endelse
endfor
nevent=ievent+1

; read the storm phase and flag the cusp into different storm phases
OPENR, unit_temp, 'storm_phase_long.dat', /GET_LUN
prestorm_start = DBLARR(300)
storm_onset = DBLARR(300) 
min_dst_new = DBLARR(300)
recovery_fast_end = DBLARR(300) 
recovery_early_end = DBLARR(300) 
recovery_long_end = DBLARR(300)
jj = 0l
dummy = ''
WHILE NOT EOF(unit_temp) DO BEGIN
    READF, unit_temp, dummy
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
CLOSE, unit_temp, /all
nstorm = jj
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

mag=DBLARR(ntime)
diffflux=DBLARR(ntime,15)
distfunc=DBLARR(ntime,15)
vperp=DBLARR(ntime)
for ievent=0,nevent-1 do begin
    index = where(event_flag eq ievent,ct)
    if ct le 0 then stop
    dt=time(index(ct-1))-time(index(0))
    at = ROUND(dt/(ct-1))
    at_str = STRCOMPRESS(ROUND(at),  /REMOVE_ALL)
    timespan,time(index(0))-at*15/2.,dt+at*15,/sec
    
    sat = fix([sc(index(0))])
    sc_str=strcompress(sat,/remove_all)
    specie = [3]
    angle = [[-90, 90], [0, 360]]
    inst = inst_input & units_name = 'DIFF FLUX' & eff_table = 0
    plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 1
    en_name='ENSPEC_SC'+sc_str+'_IN0_PHI0_360_UNDIFFFLUX_SP3_ET0_All'
    ylim,en_name,30.,4.e4,1
    average_tplot_variable,en_name,at,/new
    tplot_names,en_name+'_AVG'+at_str,names=names
    get_data,names(0), data=dd
    for i=0,ct-1 do begin
        ind=sort(ABS(dd.x-time(index(i))))
        diffflux(index(i),*)=dd.y(ind(0),0:14)
    endfor 
;--add distribution function
    sat = fix([sc(index(0))])
    sc_str=strcompress(sat,/remove_all)
    specie = [3]
    angle = [[-90, 90], [0, 360]]
    inst = inst_input & units_name = 'DIST FUNC' & eff_table = 0
    plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 1
    distfunc_name='ENSPEC_SC'+sc_str+'_IN0_PHI0_360_UNDISTFUNC_SP3_ET0_All'
    ylim,distfunc_name,30.,4.e4,1
    zlim,distfunc_name,1e-10,1e-5,1
    average_tplot_variable,distfunc_name,at,/new
    tplot_names,distfunc_name+'_AVG'+at_str,names=names
    get_data,names(0),data=dd
    for i=0,ct-1 do begin
        ind=sort(ABS(dd.x-time(index(i))))
        distfunc(index(i),*)=dd.y(ind(0),0:14)
    endfor
;--add B field
    sat=fix(sc(index(0)))
    sc_str=strcompress(sat,/remove_all)
    plot_mag_from_crib,sat
    mag_total_name='MAG_SC'+sc_str+'_B_xyz_gse_T'
    average_tplot_variable,mag_total_name,at,/new
    get_data,mag_total_name+'_AVG'+at_str,data=dd
    for i=0, ct-1 do begin
        ind=sort(ABS(dd.x-time(index(i))))
        mag(index(i))=dd.y(ind(0))
    endfor
;--add vperp
    sat = fix(sc(index(0)))
    sc_str=strcompress(sat,/remove_all)
    specie = [3] & moments = ['V']
    inst=0  & eff_table=0 
    angle=[[-90.0, 90.0], [0., 360.]] & energy=[30., 40000.]
    plot_3dmom_from_crib, sat, specie, inst, moments, angle,energy, eff_table, recalc = 1,frs = 'MAG'
    v_name='TDMOM_EN0000030_0040000_SC'+sc_str+'_MTVELOCITY_SP3_ET0_All'
    v_perp,v_name
    get_data,v_name+'_V_PERP_T',data=dd
    for i=0,ct-1 do begin
        ind=sort(ABS(dd.x-time(index(i))))
        Vperp(index(i))=dd.y(ind(0))
    endfor
;--------
    if keyword_set(spec_plot)then begin
        if keyword_set(ps) then popen,'output/o_beam/cusp/cusp'+strcompress(ievent,/remove_all)+'_sc'+sc_str+'.ps',/land
        tplot,[en_name+'_AVG'+at_str,mag_total_name+'_AVG'+at_str]
        for i=0,ct-1 do timebar,time(index(i)),thick=6
        if keyword_set(ps) then pclose else stop
    endif
    tplot_names,names=names
    store_data,delete=names
endfor
;plot,en_set,distfunc(0,*),xrange=[30.,4.e4],yrange=[1e-11,1e-4],xlog=1,xstyle=1,ylog=1,ystyle=1,/nodata
;for itime=0,ntime-1 do begin 
;    oplot,en_set,distfunc(itime,*)
;endfor 
;stop
if not keyword_set(direct_cusp) then begin 
    if keyword_set(other_inputs) then n_input = 4 else n_input = ntime
    output=DBLARR(15,n_input)
    for input= 0,n_input-1 do begin     
;   if keyword_set(func_plot) then window,input
        if keyword_set(other_inputs) then begin
            if input eq 0 then begin
                n0=10. & T0=260 & V0=100 & B0 = 33695.9 & title='Seki' 
            endif
            if input eq 1  then begin
                n0=7.7 & T0=1 & V0=0.5 & B0 = 10000. & title='Su 5000km' 
            endif 
            if input eq 2 then begin ;cm-3, eV, km/s, nT
                n0=0.05  & T0 = 10  & V0=16.8 & B0 = 100. & title='Su 8Re'
            endif 
            if input eq 3 then begin 
                n0=7.7 & T0=2.6 & V0=10 & B0 = 10000. & title='Abe 6000km-9000km' 
            endif 
            if input eq 4 then begin 
                n0=0.172 & T0=672.4 & V0=82.8 & B0 = 392.8 & title='Bouhram 3.5-5Re mean' 
            endif 
        endif else begin 
            n0=density(input)
            T0=Tem(input)
            V0=Vpara_input(input)
            B0=mag(input)
            title='Bouhram 3.5-5Re, Normalized'
        endelse

        Vth0=sqrt(2*T0/mass)    ;km/s
        Vperp = 0      ; assume that O+ is already cooled down at cusp

        n = fltarr(15)             
        J0 = fltarr(15)             
        F0 = fltarr(15)             
        for ien=0,14 do begin 
            Vpara = sqrt(2*En_range(ien,*)/mass) ;km/s
            Vpara_center = sqrt(2*En_set(ien)/mass) ;km/s
            x=(Vpara-V0)/Vth0       &      y=-((Vpara-V0)/Vth0)^2
            I1=(n0/2)*(erf(x(1))-erf(x(0)))
            I2=(-n0*Vth0/2/pi^0.5)*(exp(y(1))-exp(y(0)))+(n0*V0/2)*(erf(x(1))-erf(x(0))) ;cm-3*(km/s)  ;s3/cm3km3
            n(ien)=I1           ; cm-3
            J0(ien)=I2*1e5      ;cm-2s-1
            f0(ien)=(n0/(sqrt(pi)*Vth0)^3)*exp(-(Vperp^2+(Vpara_center-V0)^2)/Vth0^2) 
        endfor 
        index=where(n eq 0,ct) & if ct gt 0 then n(index) = !values.f_nan
        index=where(j0 eq 0,ct) & if ct gt 0 then j0(index) = !values.f_nan
        index=where(f0 eq 0,ct) & if ct gt 0 then f0(index) = !values.f_nan
        if keyword_set(func_plot) then begin 
;      print,'n',n                 
;     plot,en_set,n,ytitle='Density cm!U-3!N',xtitle='Energy(eV)',xrange=[30.,4.e4],xstyle=1,xlog=1,psym=-1,ylog=1,/nodata,yrange=[1e-5,1e2],title=title,charsize=2,/noerase
;    oplot,en_set,n,color=1,psym=-1,thick=4
;     xyouts,40,1,'Total n0 = '+string(total(n,/nan))
;      print,'J0',J0
;     plot,en_set,J0,ytitle='Flux cm!U-2!Ns!U-1!N',xtitle='Energy(eV)',xrange=[30.,4.e4],xstyle=1,xlog=1,psym=-1,ylog=1,/nodata,yrange=[1e-1,1e10],title=title,charsize=2,/noerase
;    oplot,en_set,J0,color=2,psym=-1,thick=4
;     xyouts,40,1e-1,'Total J0 = '+string(total(J0,/nan))
;      print,'F0',F0
            plot,en_set,F0*ABS(Bt/B0),ytitle='Distribution Function (s!E3!N/cm!E3!N-km!E3!N)',xtitle='Energy(eV)',xrange=[30.,4.e4],xstyle=1,xlog=1,psym=-1,ylog=1,/nodata,yrange=[1e-14,1e-4],title=title,charsize=2,/noerase
            oplot,en_set,F0,color=2,psym=-1 ;,thick=1
;     xyouts,40,1e-1,'Total F0 = '+string(total(F0,/nan))
        endif 
        if unit eq 3 then output(*,input) = J0*ABS(Bt/B0) ; normalized to Bt level 
        if unit eq 2 then output(*,input) = F0*ABS(Bt/B0) ; normalized to Bt level
    endfor 
endif else begin 
    mag_all=DBLARR(15,ntime)
    for ien=0,14 do mag_all(ien,*)=mag
    if unit eq 1 then output=transpose(diffflux)*ABS(Bt/mag_all)
    if unit eq 2 then output=transpose(distfunc)*ABS(Bt/mag_all)
    if unit eq 3 then output=transpose(flux)*ABS(Bt/mag_all)
 ;   plot,en_set,distfunc(*,0),xrange=[30.,4.e4],yrange=[1e-8,1e-2],xlog=1,xstyle=1,ylog=1,ystyle=1,/nodata
   ; for itime=0,ntime-1 do begin 
  ;      oplot,en_set,output(*,itime)
 ;   endfor 
  ;  stop
endelse 
end
