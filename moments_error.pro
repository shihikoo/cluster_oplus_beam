FUNCTION test,data,ps,input_str,path
IF keyword_set(ps) THEN begin 
    popen, path +'globe/test_globe_'+input_str, /land 
    plot3d_options, log = 1
    plot3d_codif, data, zrange = [1, 100], title = input_str
    pclose
endif  
sat=4 & specie=3 & inst=0 & eff_table=0
angle=[[-90.0, 90.0], [0.0, 360.0]]  & energy=[40,40000]
moments=['V']

plot_3dmom_from_crib, sat, specie, inst, moments, angle, energy, eff_table,NEW_NAME='v_cod', INST_COORD = 1,RECALC = 1
IF keyword_set(ps) THEN popen, path +'psd/test_psd_'+input_str
slice2d_mpe, data,$
  units = 'DIST FUNC', $
  thebdata = 'B_xyz_codif', $
  vel = 'v_cod', $
  xrange=[-200,200], $
  nosun=1, $
  range=[1e-15, 1e-7], $
onecnt=0, $
  plotenergy=0, $
  nocross=1, $
  nosmooth=0, $
  nosubtract=1, $
  noolines=0, $
  plotlabel = 1, $
  cut_perp = 0, $
  cut_par = 0, $
  cut_bulk_vel = 1
IF keyword_set(ps) THEN pclose 
print,input_str
if not keyword_set(ps) then  stop

if total(data.data) ne 0 then flag =1 else flag=-1
return ,flag
end  

PRO moments_error, cal_sim = cal_sim, t_plot=t_plot, v_plot=v_plot,nt = nt, ps = ps, test=test

;if keyword_set(test) then !p.multi=[2,2,2]
IF NOT keyword_set(nt) THEN nt = 100
n_condi_set = ['n002','n20'];['n20', 'n002','n005','n01']
v_condi_set = ['typical_v']     ;,'iso_v']
t_condi_set = ['iso_t'];,'typical1_t','typical2_t']

path = 'output/o_beam/test_mom_pnoise_new/'
;path = 'output/o_beam/test_mom_new/'+strcompress(nt, /remove_all)+'/'
spawn,  'mkdir '+path
path = path+strcompress(nt, /remove_all)+'/'
spawn,  'mkdir '+path

ntt = 1 ; # of different scales(0->100 and 0-> 1000) of T are needed to loop
sat   = 4  & specie = 3
inst = 0   &  eff_table = 0       

year_set = ['2001','2002','2003','2004','2005','2006','2007','2008','2009']
;year_set = ['2009','2008','2007','2006']

fln = path+'slice2d_sim_tem'

energy_set = [31444.7d, 19398.3d, 11966.9d, 7382.39d, 4554.22d, 2809.51d, 1733.19d, 1069.21d, 659.599d, 406.909d, 251.023d, 154.857d, 95.5315d, 58.9337d, 36.3563d]
mass_o = 16d*1.67262158e-27
;----------------------------------------------------------------------
IF keyword_set(cal_sim) THEN BEGIN
    FOR iy = 0, n_elements(year_set)-1 DO BEGIN
        next_year:
        year = year_set(iy)
        time = year+'-10-01/06:00:00'
        IF year EQ '2005' THEN time = year+'-09-20/06:00:00'
        IF year EQ '2001' THEN time = year+'-08-17/06:00:00'
        timespan, time, 300, /SEC 
        units_name = 'Counts'
        plot_globe_from_crib, sat,specie,inst,units_name, eff_table     
        name = 'GLOBE_SC'+string(sat, format = '(i1.1)')+  '_'+strcompress(units_name, /remove_all)  +'*'+'SP'+string(specie, format = '(i1.1)')
        tplot_names, name, names = gname
        FOR in = 0, n_elements(n_condi_set)-1 DO BEGIN 
            next_n:
            n_condi = n_condi_set(in)
            for iv=0,n_elements(v_condi_set)-1 do begin
                next_v:
                v_condi=v_condi_set(iv)
                for item=0,n_elements(t_condi_set)-1 do begin
                    next_t:
                    t_condi=t_condi_set(item) 
                    
               ;     plot_path = path+year+'/' & spawn,  'mkdir '+plot_path
               ;     plot_path = path+year+'/'+n_condi+'/' & spawn,  'mkdir '+plot_path
               ;     plot_path = plot_path+v_condi+'_'+t_condi+'/'& spawn,  'mkdir '+plot_path
                    if keyword_set(test) then begin
                        spawn,  'mkdir '+path+'test/'
                        spawn,  'mkdir '+path+'test/'+year+'/'
                        spawn,  'mkdir '+path+'test/'+year+'/'+n_condi+'/'
                        spawn,  'mkdir '+path+'test/'+year+'/'+n_condi+'/globe/'
                        spawn,  'mkdir '+path+'test/'+year+'/'+n_condi+'/psd/'
                    endif                     
                    n_en_sort = n_elements(energy_set)       
                    t_sim = fltarr(nt*ntt, 3, n_elements(energy_set))
                    t_in = fltarr(nt*ntt, 3, n_elements(energy_set))
                    t_tot_in = fltarr(nt*ntt, n_elements(energy_set))
                    t_tot_sim = fltarr(nt*ntt, n_elements(energy_set))
                    v_sim = fltarr(nt*ntt, 2, n_elements(energy_set))
                    v_in = fltarr(nt*ntt, 2, n_elements(energy_set))
                    v_tot_in = fltarr(nt*ntt, n_elements(energy_set))
                    v_tot_sim = fltarr(nt*ntt, n_elements(energy_set))
                    IF n_condi EQ 'n002' THEN n_input = 0.002d ; # / cm^-3
                    IF n_condi EQ 'n20' THEN n_input = 20d ; # / cm^-3
                    IF n_condi EQ 'n005' THEN n_input = 0.005d ; # / cm^-3
                    IF n_condi EQ 'n01' THEN n_input = 0.01d ; # / cm^-3
                    
               ;     if keyword_set(test) then n_input=0.005d
                    FOR ien_sort = 0, n_en_sort-1 DO BEGIN
                        energy_input = energy_set(ien_sort)
                        if keyword_set(test) then v_tot=69 else v_tot = sqrt(2*energy_input*1.6e-19/mass_o)/1e3
                        FOR itt = 0, ntt-1 DO BEGIN 
                            FOR it = 0, nt-1 DO BEGIN 
;loop for different temperature input with it                        
                                if keyword_set(test) then t_input=45 else t_input = 40000.^(1.*it/nt)
                                get_data, gname(0), data = data
                                n0 = n_input & nn = n0 & n = n0
                                next:
                                IF v_condi EQ 'typical_v' THEN BEGIN
                                ;   v01 = v_tot*0.7 &  v02 = v_tot*0.15 & v03 = v_tot*0.15 ; GSE km/s
                                ; It will be the best if we can set up
                                ; the parallel and perpendicular velocity
                                    v01=v_tot*sqrt(9/10) & v02=v_tot*sqrt(1/20) & v03=v_tot*sqrt(1/20) ;GSE km/s
                                ENDIF 
                                IF v_condi EQ 'iso_v' THEN BEGIN
                                    v01 = v_tot*sqrt(1/3)  & v02 = v_tot*sqrt(1/3) & v03 = v_tot*sqrt(1/3) ;GSE km/s
                                ENDIF
                                IF t_condi EQ 'iso_t'THEN BEGIN
                                    T01 = t_input*sqrt(1/3) & T02 = t_input*sqrt(1/3) & T03 = t_input*sqrt(1/3) ; eV
                                ENDIF
                                IF t_condi EQ 'typical1_t'THEN BEGIN 
                                    T01 = t_input*0.6 & T02 = t_input*0.2 & T03 = t_input*0.2 ; eV
                                ENDIF
                                IF t_condi EQ 'typical2_t'THEN BEGIN 
                                    T01 = t_input*0.14 & T02 = t_input*0.43 & T03 = t_input*0.43 ; eV
                                ENDIF

                                v_true=[v01,v02,v03] ;record the true value of v and t
                                t_true=[t01,t02,t03]
                                units_name = 'EFLUX'
                                full = 0d   &     pnoise = 1d &  noise = 0d
                                mp  = 1.67e-27 ; proton mass kg
                                pi  = !pi  &  K   = 1.38e-23 ; J / Kelvin
                                v01 = v01 * 1000. & v02 = v02 * 1000. & v03 = v03 * 1000. ; m/s                                                        
; Calculate velocity in instrument coordinates
                                datastr = {x:data.time, y:[[v01], [v02], [v03]]}
                                store_data, 'dd', data = datastr,  dlim = {inst_num:0, sens:data.sensitivity, sat:data.sat,  phase_instr:data.phase_instr}
                                cis_coord_trans, DATA_IN = 'dd', TRANS = 'GSE->CODIF', DATA_OUT = 'input_vel_instr'
                                get_data, 'input_vel_instr', data = in_vel
                                v01 = in_vel.y(0) & v02 = in_vel.y(1) & v03 = in_vel.y(2)
                                IF specie EQ 0 THEN m = mp
                                IF specie Eq 3 THEN m = 16 * mp
                                n = n * 1e6 ; # / m^3
                                n = n / 8640. ; CORRECTION FACTOR!!!!
                                T01 = T01 * 11600.0 &   T02 = T02 * 11600.0 &  T03 = T03 * 11600.0 ; Kelvin
                                nbins   = data.nbins & nenergy = data.nenergy
                                corfac = [      4,      4,      4,      4, $
                                                2,  2,  2,  2,  2,  2,  2,  2, $
                                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, $
                                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, $
                                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, $
                                                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, $
                                                2,  2,  2,  2,  2,  2,  2,  2, $
                                                4,      4,      4,      4]
; Calculate the f(v)
                                J = FLTARR(nenergy, nbins)
                                FOR ien = 0, nenergy-1 DO BEGIN
                                    FOR ibin = 0, nbins-1 DO BEGIN
                                        theta  = data.theta(ien, ibin)  *!dtor
                                        phi    = data.phi(ien, ibin)    *!dtor
                                        energ  = data.energy(ien, ibin)
                                        u = SQRT(2 * energ * 1.6e-19 / m)
                                        ux = u*COS(theta)*COS(phi) - v01
                                        uy = u*COS(theta)*SIN(phi) - v02
                                        uz = u*SIN(theta)          - v03
                                        f = (n/SQRT(T01*T02*T03))*(m/(2*pi*K))^(3./2)*EXP(-m*ux^2/(2*K*T01)-m*uy^2/(2*K*T02)-m*uz^2/(2*K*T03)) ;*corfac(ibin)
                                        J(ien, ibin) = 0.5 * u^4 * f ; Differential Energy Flux
                                    ENDFOR
                                ENDFOR
                                
; Convert EFLUX in Counts
                                data.data = (J)
                                data = reverse_convert_codif_units(data, units_name, 'codif_ts_eff_corr',  eff_table,  specie = specie,  packets = 1, sat = sat)
                                data.units_name = 'Counts'
                                IF keyword_set(test)  THEN valid=test(data,ps,'addefficiency_n'+strcompress(n_input,/remove_all)+'_v'+strcompress(v_tot,/remove_all)+'_t'+strcompress(t_input,/remove_all)+'year'+year+'_'+v_condi+'_'+t_condi,path+'test/'+year+'/'+n_condi+'/')

; Randomize and integerize
                                IF NOT KEYWORD_SET(full) THEN BEGIN
                                    data.data = data.data + pnoise * (SQRT(data.data)*(RANDOMU(seed, 31, 88))) + noise * (2*RANDOMU(seed, 31, 88)-1)
                                    IF keyword_set(test)  THEN valid=test(data,ps,'noise_n'+strcompress(n_input,/remove_all)+'_v'+strcompress(v_tot,/remove_all)+'_t'+strcompress(t_input,/remove_all)+'year'+year+'_'+v_condi+'_'+t_condi,path+'test/'+year+'/'+n_condi+'/')
                                    ibad = WHERE(data.data LT 1.0, cibad)
                                    IF cibad GT 0 THEN data.data(ibad) = 0.0
                                    IF keyword_set(test)  THEN valid=test(data,ps,'deletelowdata_n'+strcompress(n_input,/remove_all)+'_v'+strcompress(v_tot,/remove_all)+'_t'+strcompress(t_input,/remove_all)+'_year'+year+'_'+v_condi+'_'+t_condi,path+'test/'+year+'/'+n_condi+'/')
                                    data.data = ROUND(data.data)
                                ENDIF
                                IF keyword_set(test)  THEN valid=test(data,ps,'digitalize_n'+strcompress(n_input,/remove_all)+'_v'+strcompress(v_tot,/remove_all)+'_t'+strcompress(t_input,/remove_all)+'year'+year+'_'+v_condi+'_'+t_condi,path+'test/'+year+'/'+n_condi+'/')
; calculate the new density, if density is lower than input loop the
; calculation with a higher density
                                datadummy = data
                                datef = convert_codif_units(datadummy, 'EFLUX', 'codif_ts_eff_corr',  eff_table, packets = 1, sat = sat)
                                angle = [[-90.0, 90.0], [0.0, 360.0]] & energy = [30.0, 40000.0]
                                density_sim = n_3d_cis(datef, ENERGY = energy, ANGLE = angle)
                                
                                if keyword_set(test) then begin
                                    if item ne n_elements(t_condi_set)-1 then begin
                                        item=item+1 
                                        goto, next_t 
                                    endif  else begin
                                   
                                        if iv ne n_elements(v_condi_set)-1 then begin
                                            iv=iv+1
           
                                            goto, next_v 
                                        endif  else begin
                                     
                                            if in ne n_elements(n_condi_set)-1 then begin
                                                in=in+1
                                                goto, next_n 
                                            endif  else begin
                                                if iy ne n_elements(year_set) -1 then begin
                                                    iy=iy+1
                                                    goto, next_year 
                                                endif  else stop
                                            endelse 
                                        endelse 
                                    endelse 
                                endif   
                                IF density_sim LT n0  THEN BEGIN
                                    nn = 1.1 * nn   &   n = nn ; & print, nn, density_sim
                                    GOTO, next
                                ENDIF
;--Jing: add temperature calculation part 
                                angle = [[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
                                energy = [30.0, 40000.0]
                                tt = compute_temperature(datef, sat, NAME = 'sim_tem',  ENERGY = energy, ANGLE = angle)
                                t_sim(itt*nt+it, *, ien_sort) = tt
                                t_in(itt*nt+it, *, ien_sort) = t_true
;--Jing add the velocity calculation part
                                angle = [[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
                                energy = [30.0, 40000.0]
                                vv = compute_velocity(datef, sat, inst, NAME = 'sim_vel', ENERGY = energy, ANGLE = angle)
                                plot_mag_from_crib, sat ; Load CLUSTER Magnetic field for calculate v_perp 
                                v_perp,'sim_vel',mag='MAG_SC4_B_xyz_gse'
                                get_data,'sim_vel_V_PAR_T',data=dd   & v_para=dd.y
                                get_data,'sim_vel_V_PERP_T',data=dd  & v_perp=dd.y
                                v_sim(itt*nt+it, *, ien_sort) = [v_para,v_perp]
                                store_data,'input_vel',data={x:dd.x,y:reform(v_true,1,3)}
                                v_perp,'input_vel',mag='MAG_SC4_B_xyz_gse'
                                get_data,'input_vel_V_PAR_T',data=dd & v_para=dd.y
                                get_data,'input_vel_V_PERP_T',data=dd  & v_perp=dd.y
                                v_in(itt*nt+it, *, ien_sort) = [v_para,v_perp]
                            ENDFOR    
                        ENDFOR   
                        index = sort(t_in(*, 0))
                        t_sim(*, 0, ien_sort) = t_sim(index, 0, ien_sort)
                        t_sim(*, 1, ien_sort) = t_sim(index, 1, ien_sort)
                        t_sim(*, 2, ien_sort) = t_sim(index, 2, ien_sort)
                        t_in(*, 0, ien_sort) = t_in(index, 0, ien_sort)
                        t_in(*, 1, ien_sort) = t_in(index, 1, ien_sort)
                        t_in(*, 2, ien_sort) = t_in(index, 2, ien_sort)
                        t_tot_in(*, ien_sort) = (t_in(*, 0, ien_sort)+t_in(*, 1, ien_sort)+t_in(*, 2, ien_sort))/3
                        t_tot_sim(*, ien_sort) = (t_sim(*, 0, ien_sort)+t_sim(*, 1, ien_sort)+t_sim(*, 2, ien_sort))/3
                        
                        index = sort(v_in(*, 0))
                        v_sim(*, 0, ien_sort) = v_sim(index, 0, ien_sort)
                        v_sim(*, 1, ien_sort) = v_sim(index, 1, ien_sort)
                        v_in(*, 0, ien_sort) = v_in(index, 0, ien_sort)
                        v_in(*, 1, ien_sort) = v_in(index, 1, ien_sort)

                        v_tot_in(*, ien_sort) = sqrt(v_in(index,0,ien_sort)^2+v_in(index,1,ien_sort)^2)
                        v_tot_sim(*, ien_sort) = sqrt(v_sim(index,0,ien_sort)^2+v_sim(index,1,ien_sort)^2)
                    ENDFOR        
                    store_data, year+'_'+n_condi+'_T_tot', data = {x:t_tot_in, y:t_tot_sim/t_tot_in}
                    store_data, year+'_'+n_condi+'_T_para', data = {x:t_in(*, 0, *), y:t_sim(*, 0, *)/t_in(*, 0, *)}
                    store_data, year+'_'+n_condi+'_T_perp1', data = {x:t_in(*, 1, *), y:t_sim(*, 1, *)/t_in(*, 1, *)}
                    store_data, year+'_'+n_condi+'_T_perp2', data = {x:t_in(*, 2, *), y:t_sim(*, 2, *)/t_in(*, 2, *)}
                    
                    store_data, year+'_'+n_condi+'_V_tot', data = {x:v_tot_in, y:ABS(v_tot_sim-v_tot_in)/v_tot_in}
                    store_data, year+'_'+n_condi+'_V_para', data = {x:v_in(*, 0, *), y:ABS(v_sim(*, 0, *)-v_in(*,0,*))/v_in(*, 0, *)}
                    store_data, year+'_'+n_condi+'_V_perp', data = {x:v_in(*, 1, *), y:ABS(v_sim(*, 1, *)-v_in(*,0,*))/v_in(*, 1, *)}

                    tplot_names, '*'+n_condi+'_T*', names = names1
                    tplot_names, '*'+n_condi+'_V*',names=names2
                    tplot_save, [names1,names2], filename = fln

                ENDFOR     
            ENDFOR 
        endfor
    endfor     
ENDIF   

IF keyword_set(v_plot) THEN BEGIN 
    n_en_sort = n_elements(energy_set)
    n_year=n_elements(year_set)
    t_tot_in_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    v_tot_in_set = FLTARR(n_year, 2, nt*ntt, n_en_sort) 
    v_ratio_tot_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    v_para_in_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    v_ratio_para_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    v_perp_in_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    v_ratio_perp1_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    tplot_restore, filename = fln+'.tplot'
    FOR iy = 0, n_year-1 DO BEGIN 
        year = year_set(iy)
        FOR  in = 0, n_elements(n_condi_set)-1 DO BEGIN
            n_condi = n_condi_set(in)
            get_data, year+'_'+n_condi+'_V_tot', data = data
            v_tot_in_set(iy,in,*,*)=data.x
            v_ratio_tot_set(iy, in, *, *)=data.y

            get_data, year+'_'+n_condi+'_T_tot', data = data
            t_tot_in_set(iy,in,*,*)=data.x
        ENDFOR 
    ENDFOR 
    FOR iy = 0, n_year-1 DO BEGIN 
        FOR in = 0, n_elements(n_condi_set)-1 DO BEGIN 
            spawn, 'mkdir ' + path+n_condi_set(in)+'/'
            year = year_set(iy)
            v_input = reform(t_tot_in_set(iy, in, *, *))
            error = ABS(reform((v_ratio_tot_set(iy, in, *, *))))
            x_axis = v_input(*, 0)
            y_axis = energy_set
            v = error
            
            IF keyword_set(ps) THEN popen, path+n_condi_set(in)+'/'+Year_set(iy)+'_'+n_condi_set(in)+'_3d_Verror.ps', /land
            specplot, x_axis, y_axis, v, no_interp = 1, $
              lim = { zlog:0, zrange: [0, 1], xlog:1, ylog:1,  xrange: [1, 40000.], yrange: [30, 40000.], $
                      title: year_set(iy)+'    '+n_condi_set(in), xtitle: 'T (eV)', ytitle: 'Energy (eV)', ztitle: 'V error', $
                      XSTYLE:1, ystyle: 1, charsize: 1.2, position: [0.1, 0.1, 0.9, 0.9]}   
            IF keyword_set(ps) THEN pclose else stop
        ENDFOR 
    ENDFOR   
endif 
IF keyword_set(t_plot) THEN BEGIN 
    n_en_sort = n_elements(energy_set)
    n_year=n_elements(year_set)
    t_tot_in_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    t_ratio_tot_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    t_para_in_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    t_ratio_para_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    t_perp1_in_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    t_ratio_perp1_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    t_perp2_in_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    t_ratio_perp2_set = FLTARR(n_year, 2, nt*ntt, n_en_sort)
    tplot_restore, filename = fln+'.tplot'
    FOR iy = 0, n_year-1 DO BEGIN 
        year = year_set(iy)
        FOR  in = 0, n_elements(n_condi_set)-1 DO BEGIN
            n_condi = n_condi_set(in)
            get_data, year+'_'+n_condi+'_T_tot', data = data
            t_tot_in_set(iy, in, *, *) = data.x
            t_ratio_tot_set(iy, in, *, *) = data.y
            get_data, year+'_'+n_condi+'_T_para', data = data
            t_para_in_set(iy, in, *, *) = data.x
            t_ratio_para_set(iy, in, *, *) = data.y
            get_data, year+'_'+n_condi+'_T_perp1', data = data
            t_perp1_in_set(iy, in, *, *) = data.x
            t_ratio_perp1_set(iy, in, *, *) = data.y
            get_data, year+'_'+n_condi+'_T_perp1', data = data
            t_perp2_in_set(iy, in, *, *) = data.x
            t_ratio_perp2_set(iy, in, *, *) = data.y
        ENDFOR 
    ENDFOR 
    FOR iy = 0, n_year-1 DO BEGIN 
        FOR in = 0, n_elements(n_condi_set)-1 DO BEGIN 
            year = year_set(iy)
            t_input = reform(t_tot_in_set(iy, in, *, *))
            error = ABS(reform((t_ratio_tot_set(iy, in, *, *)))-1)
            x_axis = t_input(*, 0)
            y_axis = energy_set
            v = error
            IF keyword_set(ps) THEN popen, path+n_condi_set(in)+'/'+Year_set(iy)+'_'+n_condi_set(in)+'_3d_Terror.ps', /land
            specplot, x_axis, y_axis, v,no_interp = 1, $
              lim = { zlog:0, zrange: [0, 1], $
                      title: year_set(iy)+'    '+n_condi_set(in),  $
                      xtitle: 'T (eV)', $
                      ytitle: 'Energy (eV)', $
                      ztitle: 'error', $
                      xrange: [1, 40000.], yrange: [30, 40000.], $
                      xlog:1, ylog:1, $
                      XSTYLE:1, ystyle: 1, charsize: 1.2, $
                      position: [0.1, 0.1, 0.9, 0.9]}   
            iF keyword_set(ps) THEN pclose
            IF NOT keyword_set(ps) THEN stop
        ENDFOR 
    ENDFOR  
ENDIF     
stop
END 
