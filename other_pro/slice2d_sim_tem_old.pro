PRO slice2d_sim_tem,cal_sim_tem = cal_sim_tem, plot_sim_tem = plot_sim_tem, ps = ps, interp_terr_plot = interp_terr_plot

; cal_sim_tem = 0
; plot_sim_tem = 0
Terr_plot_3d = 1 ; plot the Temperature error plot in 3d with energy and temperature as parameters
sim = 1         ; set as 1 if you wanna sim data
raw_counts = 0  ; used for plot the globe counts data before apply efficients

plot_globe = 0
plot_dt = 0
plot_log = 0
v_condi_set = ['sort_en', 'typical_v', 'neg_v', 'iso_v', 'iso_neg_v']
v_condi = v_condi_set(0)
t_condi_set = ['typical_t', 'para_t', 'perp_t']
t_condi = t_condi_set(0)
path = 'output/o_beam/test_mom/'
spawn,  'mkdir '+path
path = path+v_condi+'_'+t_condi+'/'
spawn,  'mkdir '+path

ntt = 1 ; # of different scales(0->100 and 0-> 1000) of T are needed to loop
sat   = 4
specie = 3
inst = 0                        ; 0: CODIF, 1: HIA (this is not supported for the moment)
eff_table = 0                   ; 0: GROUND, 1: ONBOARD

year_set = ['2001', '2002', '2003', '2004', '2005']
n_condi_set = ['low_n', 'high_n']
fln = path+'slice2d_sim_tem'
spawn,  'mkdir '+path

; energy bins 
energy_set = [31444.7d, 19398.3d, 11966.9d, 7382.39d, 4554.22d, 2809.51d, 1733.19d, 1069.21d, $
              659.599d, 406.909d, 251.023d, 154.857d, 95.5315d, 58.9337d, 36.3563d] 
mass_o = 16d*1.67262158e-27
nt = 100
;----------------------------------------------------------------------
IF keyword_set(cal_sim_tem)THEN BEGIN
    FOR iy = 0, 4 DO BEGIN 
        newyear:
        year = year_set(iy)
        spawn,  'mkdir '+path+year+'/'
        time = year+'-10-01/06:00:00'
        IF year EQ '2005' THEN time = year+'-09-20/06:00:00'
        IF year EQ '2001' THEN time = year+'-08-17/06:00:00'
        timespan, time, 300, /SEC ; SECONDS, MINUTES, HOURS, DAYS (DEFAULT)
        units_name = 'Counts'
;----------------------------------------------------------------------
; Load/calculate the globe plot data
; Keywords: CNES -> use CNES magnetic field data
;           IC   -> use Imperial College magnetic field data
; The PP (cdf format) magnetic field data is the default setting
;----------------------------------------------------------------------
        plot_globe_from_crib, sat, $
          specie, $
          inst, $
          units_name, $
          BKG = 0, $
          eff_table, $
          OLD_EFF = 0, $
          CNES = 0, $
          IC = 0
        name = 'GLOBE_SC'+string(sat, format = '(i1.1)')+$
               '_'+strcompress(units_name, /remove_all)  +$
               '*'+'SP'+string(specie, format = '(i1.1)')
        tplot_names, name, names = gname

        FOR in = 0, 1 DO BEGIN 
            n_condi = n_condi_set(in)
            plot_path = path+year+'/'
            spawn,  'mkdir '+plot_path
            plot_path = plot_path+n_condi+'/'
            spawn,  'mkdir '+plot_path
            spawn,  'mkdir '+plot_path+'globe/'

            IF v_condi EQ 'sort_en' THEN  n_en_sort = n_elements(energy_set) ELSE  n_en_sort = 1             

            t_sim = fltarr(nt*ntt, 3, n_elements(energy_set))
            t_in = fltarr(nt*ntt, 3, n_elements(energy_set))
            t_tot_in = fltarr(nt*ntt, n_elements(energy_set))
            t_tot_sim = fltarr(nt*ntt, n_elements(energy_set))

    ;        nV = 0.036d         ;0.013466d *71.078383d   ;0.012d   
       
            IF n_condi EQ 'low_n' THEN n_input = 0.002d   ; # / cm^3
            IF n_condi EQ 'high_n' THEN n_input = 20d ; # / cm^3
           
            FOR ien_sort = 0, n_en_sort-1 DO BEGIN
;default enregy set to 440 eV
                IF v_condi EQ 'sort_en' THEN energy_input = energy_set(ien_sort) ELSE energy_input = 453.
                v_tot = sqrt(2*energy_input*1.6e-19/mass_o)/1e3
                FOR itt = 0, ntt-1 DO BEGIN 
                FOR it = 0, nt-1 DO BEGIN 
                    IF keyword_set(raw_counts) THEN t_input = 15 ELSE t_input = energy_input^(it/nt)
; t_input = (it+1)*(energy_input/nt)*(0.25^itt)
               
                    get_data, gname(0), data = data
;----------------------------------------------------------------------
; Fill the globes with sim values
;----------------------------------------------------------------------
                    IF KEYWORD_SET(sim) THEN BEGIN
;-----------------------------------------------------------------
; Input
;-----------------------------------------------------------------
;Density_tail     V_total_tail    V_par_tail      V_perp_tail     T_total_tail
;0.013466        72.908554        71.078383        16.233350        12.227400    
;T_x_tail         T_y_tail         T_z_tail          
;35.773937         0.454132        11.774950
                        n0 = n_input
                        nn = n0
                        n = n0

                        next:
; in the condition sort_en, we run over different energies to get a
; seris of tempreture distribution for all  those energies 
                        IF v_condi EQ 'sort_en' THEN BEGIN 
                            v01 = v_tot*0.7 ; GSE km/s
                            v02 = v_tot*0.15 ; GSE km/s
                            v03 = v_tot*0.15 ; GSE km/s
                        ENDIF 
                        IF v_condi EQ 'typical_v' THEN BEGIN 
                            v01 = 71.078383d ; GSE km/s
                            v02 = 11.4784d ; GSE km/s
                            v03 = 11.4784d ; GSE km/s
                        ENDIF 
                        IF v_condi EQ 'neg_v' THEN BEGIN 
                            v01 = -71.078383d ; GSE km/s
                            v02 = 11.4784d ; GSE km/s
                            v03 = 11.4784d ; GSE km/s
                        ENDIF 
                        IF v_condi EQ 'iso_v' THEN BEGIN 
                            v01 = 42.093773d ; GSE km/s
                            v02 = 42.093773d ; GSE km/s
                            v03 = 42.093773d ; GSE km/s
                        ENDIF 
                        IF v_condi EQ 'iso_neg_v' THEN BEGIN 
                            v01 = -42.093773d ; GSE km/s
                            v02 = 42.093773d ; GSE km/s
                            v03 = 42.093773d ; GSE km/s
                        ENDIF 

                        IF t_condi EQ 'typical_t'THEN BEGIN 
                            T01 = t_input *0.001 ; keV
                            T02 = t_input *0.001 ; keV
                            T03 = t_input *0.001 ; keV
                        ENDIF 
                        IF t_condi EQ 'para_t'THEN BEGIN 
                            T01 = t_input *0.001 ; keV
                            T02 = 25.2*0.001 ; keV
                            T03 = 25.2*0.001 ; keV
                        ENDIF 
                        IF t_condi EQ 'perp_t'THEN BEGIN 
                            T01 = 13.5 *0.001 ; keV
                            T02 = t_input *0.001 ; keV
                            T03 = t_input *0.001 ; keV
                        ENDIF 
                        units_name = 'EFLUX'
                        full = 0d
                        pnoise = 0d
                        noise = 0d
;-----------------------------------------------------------------
                        mp  = 1.67e-27 ; proton mass kg
                        pi  = !pi
                        K   = 1.38e-23 ; J / Kelvin

                        v01 = v01 * 1000. ; m/s
                        v02 = v02 * 1000. ; m/s
                        v03 = v03 * 1000. ; m/s
                        
                                ; Calculate velocity in instrument coordinates
                        datastr = {x:data.time, y:[[v01], [v02], [v03]]}
                        store_data, 'dd', data = datastr, $
                                    dlim = {inst_num:0, sens:data.sensitivity, sat:data.sat, $
                                            phase_instr:data.phase_instr}
                        cis_coord_trans, DATA_IN = 'dd', TRANS = 'GSE->CODIF', $
                                         DATA_OUT = 'input_vel_instr'

                        get_data, 'input_vel_instr', data = in_vel
                        v01 = in_vel.y(0)
                        v02 = in_vel.y(1)
                        v03 = in_vel.y(2)

                        IF specie EQ 0 THEN m = mp
                        IF specie Eq 3 THEN m = 16 * mp

                        n = n * 1e6 ; # / m^3
                        n = n / 8640. ; CORRECTION FACTOR!!!!
                        T01 = 1e3 * T01 * 11600.0 ; Kelvin
                        T02 = 1e3 * T02 * 11600.0 ; Kelvin
                        T03 = 1e3 * T03 * 11600.0 ; Kelvin

                        nbins   = data.nbins
                        nenergy = data.nenergy

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

                                f = (n/SQRT((T01*T02*T03)))*(m/(2*pi*K))^(3./2)*$
                                    EXP(-m*ux^2/(2*K*T01)-m*uy^2/(2*K*T02)-m*uz^2/(2*K*T03)) ;*corfac(ibin)

                                J(ien, ibin) = 0.5 * u^4 * f ; Differential Energy Flux

                            ENDFOR
                        ENDFOR

                                ; Convert EFLUX in Counts
                        data.data = (J)
                        data = reverse_convert_codif_units(data, units_name, 'codif_ts_eff_corr', $
                                                           eff_table, $
                                                           specie = specie, $
                                                           packets = 1, $
                                                           old_eff = old_eff, $
                                                           incrates = incrates, $
                                                           sat = sat)
                        
                        data.units_name = 'Counts'

                                ; Randomize and integerize
                        IF NOT KEYWORD_SET(full) THEN BEGIN
                            data.data = data.data + $
                                        pnoise * (SQRT(data.data)*(RANDOMU(seed, 31, 88))) + $
                                        noise * (2*RANDOMU(seed, 31, 88)-1)
                            ibad = WHERE(data.data LT 1.0, cibad)
                            IF cibad GT 0 THEN BEGIN
                                data.data(ibad) = 0.0
                            ENDIF
                            data.data = ROUND(data.data)
                        ENDIF

;-------------------------------------------------------------------
; Calculate new density and velocity
;-------------------------------------------------------------------
                                ; Density
                        datadummy = data
                        datef = convert_codif_units(datadummy, 'EFLUX', 'codif_ts_eff_corr', $
                                                    eff_table, packets = 1, $
                                                    old_eff = old_eff, sat = sat)
                        angle = [[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
                        energy = [30.0, 40000.0]
                        density_sim = n_3d_cis(datef, ENERGY = energy, ANGLE = angle)

                                ; Velocity In instrument co-ordinates
                        angle = [[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
                        energy = [30.0, 40000.0]
                        vel = compute_velocity(datef, $
                                               sat, $
                                               inst, $
                                               NAME = 'v_cod', $
                                               ENERGY = energy, $
                                               ANGLE = angle, $
                                               INST_COORD = 1)
                        
                                ; Velocity In GSE co-ordinates
                        vel = compute_velocity(datef, $
                                               sat, $
                                               inst, $
                                               NAME = 'v_cod_gse', $
                                               ENERGY = energy, $
                                               ANGLE = angle, $
                                               INST_COORD = 0)
                        
                        get_data, 'v_cod_gse', data = ss
                        velocity_sim = ss.y
                        
                        IF keyword_set(raw_counts) AND iy LE 4 THEN BEGIN 
                            IF keyword_set(ps) THEN popen, path +'counts'+year, /land ELSE window, 1
                            plot3d_options, log = 1
                            plot3d_codif, data, zrange = [1, 100], title = year+', n_min:'+string(n_input) +', nV:'+string(nV)             
                            IF keyword_set(ps) THEN pclose
                            iy = iy+1
                            IF iy GT 4 THEN BEGIN 
                                IF keyword_set(ps) THEN BEGIN 
                                    spawn, 'mogrify -format png '+path+'counts*.ps '
                                    spawn, 'mogrify -rotate -90 '+path+'counts*.png &'
                                ENDIF 
                                stop
                            ENDIF 
                            GOTO, newyear
                        ENDIF 

                        IF density_sim LT n0 THEN BEGIN
                            nn = 1.1 * nn
                            n = nn
                            print, nn, density_sim
                            GOTO, next
                        ENDIF
                    ENDIF   

;Jing: add temperature calculation part 
                    angle = [[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
                    energy = [30.0, 40000.0]
                    tt = compute_temperature(datef, sat, NAME = 'sim_tem', $
                                             ENERGY = energy, ANGLE = angle)

                    IF keyword_set(plot_globe) THEN BEGIN 
                        IF keyword_set(ps) THEN popen, plot_path +'globe/tem_simu'+STRING(t_input, format = '(i4.4)')+'.ps', /land ELSE window, /free
                        plot3d_options, log = 1
                        plot3d_codif, data, zrange = [1, 100], $
                                      title = 'INPUT: T!Lpara!N = '+string(t01/11600., format = '(f6.2)') $
                                      +'eV   T!Lperp!N = '+ string(t02/11600., format = '(f6.2)')+'eV, ' $
                                      + string(t03/11600., format = '(f6.2)')+'eV!C'+ $    
                                      ' SIM  : T!Lpara!N = '+string(tt(0), format = '(f6.2)') $
                                      +'eV   T!Lperp!N = '+ string(tt(1), format = '(f6.2)')+'eV, ' $
                                      + string(tt(2), format = '(f6.2)')+'eV'     
                        IF keyword_set(ps) THEN pclose
                    ENDIF 
                    t_sim(itt*nt+it, *, ien_sort) = tt
                    t_in(itt*nt+it, *, ien_sort) = [t01, t02, t03]/11600.
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
            ENDFOR        
        
          !p.multi = [0, 1, 2]
;draw for dt
            IF keyword_set(plot_dt) THEN BEGIN 
                dt = abs(t_sim-t_in)/t_in
                dt_tot = abs(t_tot_sim-t_tot_in)/t_tot_in
                IF keyword_set(ps) THEN popen, plot_path+'dT_test_tot.ps', /port ELSE window, 2
                plot, t_tot_in, dt_tot, $
                      xrange = [0, energy_input/4], yrange = [0, 2], xstyle = 1, $
                      xtitle = 'T input ( eV )', ytitle = 'dT / Tinput'
                xyouts, 80, 1.3, 'T total'
                xyouts, 80, 1.5, year
                oplot, [0, 100], [0, 0], color = 2
                plot, t_tot_in, dt_tot, $
                      xrange = [0, energy_input], yrange = [0, 2], xstyle = 1, $
                      xtitle = 'T input ( eV )', ytitle = 'dT / T_input'
                xyouts, 800, 1.3, 'T total'
                xyouts, 800, 1.5, year
                oplot, [0, 450], [0, 0], color = 2
                IF keyword_set(ps) THEN pclose

                FOR it = 0, 2 DO BEGIN 
                    IF keyword_set(ps) THEN $
                      popen, plot_path+'dT_test_t'+string(it, format = '(i1)')+'.ps', /port $
                    ELSE window, it+3
                    plot, t_in(*, it), dt(*, it), $
                          xrange = [0, energy_input/4], yrange = [0, 2], xstyle = 1, $
                          xtitle = 'T input ( eV )', ytitle = 'dT / T_input'
                    xyouts, 80, 1.3, 'T'+string(it, format = '(i1)')
                    xyouts, 80, 1.5, year
                    oplot, [0, 100], [0, 0], color = 2
                    plot, t_in(*, it), dt(*, it), $
                          xrange = [0, energy_input], yrange = [0, 2], xstyle = 1, $
                          xtitle = 'T input ( eV )', ytitle = 'dT / T_input'
                    xyouts, 800, 1.3, 'T'+string(it, format = '(i1)')
                    xyouts, 800, 1.5, year
                    oplot, [0, 450], [0, 0], color = 2 
                    IF keyword_set(ps) THEN pclose
                ENDFOR  
            
            ENDIF 
;for T
            FOR j_en_sort = 0, n_en_sort-1 DO BEGIN 
                energy_input = energy_set(j_en_sort)
                energy_input_str = string(energy_input, format = '(i5.5)')
                IF keyword_set(ps) THEN popen, plot_path+'T_test_tot_'$
                  + energy_input_str+'.ps', /port ELSE window, 2
                plot, t_tot_in(*, j_en_sort), t_tot_sim(*, j_en_sort)/t_tot_in(*, j_en_sort), $
                      xrange = [0, energy_input/4], yrange = [0, 2], xstyle = 1, $
                      xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
                xyouts, 80, 1.3, 'T total'
                xyouts, 80, 1.5, year
                oplot, [0, 40000], [1, 1], color = 2
                plot, t_tot_in(*, j_en_sort), t_tot_sim(*, j_en_sort)/t_tot_in(*, j_en_sort), $
                      xrange = [0, energy_input], yrange = [0, 2], xstyle = 1, $
                      xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
                xyouts, 800, 1.3, 'T total'
                xyouts, 800, 1.5, year
                oplot, [0, 40000], [1, 1], color = 2
                IF keyword_set(ps) THEN pclose
             
                FOR it = 0, 2 DO BEGIN 
                    IF keyword_set(ps) THEN $
                      popen, plot_path+'T_test_t'+string(it, format = '(i1)')+'.ps', /port $
                    ELSE window, it+3
                    plot, t_in(*, it, j_en_sort), T_sim(*, it, j_en_sort)/t_in(*, it,j_en_sort), $
                          xrange = [0, energy_input/4], yrange = [0, 2], xstyle = 1, $
                          xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
                    xyouts, 80, 1.3, 'T'+string(it, format = '(i1)')
                    xyouts, 80, 1.5, year
                    oplot, [0, 40000], [1, 1], color = 2
                    plot, t_in(*, it, j_en_sort), T_sim(*, it,j_en_sort)/t_in(*, it,j_en_sort), $
                          xrange = [0, energy_input], yrange = [0, 2],  xstyle = 1, $
                          xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
                    xyouts, 800, 1.3, 'T'+string(it, format = '(i1)')
                    xyouts, 800, 1.5, year
                    oplot, [0, 40000], [1, 1], color = 2 
                    IF keyword_set(ps) THEN pclose
                ENDFOR
            ENDFOR  
;for T in log

            IF keyword_set(plot_log) THEN BEGIN 
                IF keyword_set(ps) THEN popen, plot_path+'T_test_tot_log.ps', /port ELSE window, 2
                plot, t_tot_in, t_tot_sim/t_tot_in, $
                      xrange = [0, energy_input/4], yrange = [0.1, 10], ylog = 1, ystyle = 1,  xstyle = 1, $
                      xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
                xyouts, 80, 1.3, 'T total'
                xyouts, 80, 1.5, year
                oplot, [0, 100], [1, 1], color = 2
                plot, t_tot_in, t_tot_sim/t_tot_in, $
                      xrange = [0, energy_input], yrange = [0.1, 10], ylog = 1, ystyle = 1, xstyle = 1, $
                      xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
                xyouts, 800, 1.3, 'T total'
                xyouts, 800, 1.5, year
                oplot, [0, 450], [1, 1], color = 2
                IF keyword_set(ps) THEN pclose

                FOR it = 0, 2 DO BEGIN 
                    IF keyword_set(ps) THEN $
                      popen, plot_path+'T_test_t'+string(it, format = '(i1)')+'log.ps', /port $
                    ELSE window, it+3
                    plot, t_in(*, it), T_sim(*, it)/t_in(*, it), $
                          xrange = [0, energy_input/4], yrange = [0.1, 10], ylog = 1, ystyle = 1, xstyle = 1, $
                          xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
                    xyouts, 80, 1.3, 'T'+string(it, format = '(i1)')
                    xyouts, 80, 1.5, year
                    oplot, [0, 100], [1, 1], color = 2
                    plot, t_in(*, it), T_sim(*, it)/t_in(*, it), ylog = 1, ystyle = 1, xstyle = 1, $
                          xrange = [0, energy_input], yrange = [0.1, 10], $
                          xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
                    xyouts, 800, 1.3, 'T'+string(it, format = '(i1)')
                    xyouts, 800, 1.5, year
                    oplot, [0, 450], [1, 1], color = 2 
                    IF keyword_set(ps) THEN pclose
                ENDFOR
            ENDIF  
            spawn, 'mogrify -format png '+plot_path+'*.ps &'
            IF keyword_set(plot_globe) THEN BEGIN 
                spawn, 'mogrify -format png '+plot_path+'globe/*.ps '
                spawn, 'mogrify -rotate -90 '+plot_path+'globe/*.png &'
            ENDIF 
            !p.multi = 0
            store_data, year+'_'+n_condi+'_T_tot', data = {x:t_tot_in, y:t_tot_sim/t_tot_in}
            store_data, year+'_'+n_condi+'_T_para', data = {x:t_in(*, 0, *), y:t_sim(*, 0, *)/t_in(*, 0, *)}
            store_data, year+'_'+n_condi+'_T_perp1', data = {x:t_in(*, 1, *), y:t_sim(*, 1, *)/t_in(*, 1, *)}
            store_data, year+'_'+n_condi+'_T_perp2', data = {x:t_in(*, 2, *), y:t_sim(*, 2, *)/t_in(*, 2, *)}

          ;  print, FINDFILE(fln+'.tplot', COUNT = ct_save)
          ;  IF ct_save GT 0 THEN  tplot_restore, filename = fln+'.tplot'
            tplot_names, '*n_T*', names = names
            tplot_save, names, filename = fln
        ENDFOR   
    ENDFOR   
ENDIF  

IF keyword_set(plot_sim_tem) THEN BEGIN 
    IF v_condi EQ 'sort_en' THEN  n_en_sort = n_elements(energy_set) ELSE  n_en_sort = 1             
    t_tot_in_set = FLTARR(5, 2, nt*ntt, n_en_sort)
    t_ratio_tot_set = FLTARR(5, 2, nt*ntt, n_en_sort)
    t_para_in_set = FLTARR(5, 2, nt*ntt, n_en_sort)
    t_ratio_para_set = FLTARR(5, 2, nt*ntt, n_en_sort)
    t_perp1_in_set = FLTARR(5, 2, nt*ntt, n_en_sort)
    t_ratio_perp1_set = FLTARR(5, 2, nt*ntt, n_en_sort)
    t_perp2_in_set = FLTARR(5, 2, nt*ntt, n_en_sort)
    t_ratio_perp2_set = FLTARR(5, 2, nt*ntt, n_en_sort)

    tplot_restore, filename = fln+'.tplot'

    FOR iy = 0, 4 DO BEGIN 
        year = year_set(iy)
        FOR  in = 0, 1 DO BEGIN
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

    IF keyword_set(Terr_plot_3d) THEN BEGIN 
        FOR iy = 0, 4 DO BEGIN 
            FOR in = 0, 1 DO BEGIN 
                t_input = reform(t_tot_in_set(iy, in, *, *))
                error = ABS(reform((t_ratio_tot_set(iy, in, *, *)))-1)
                v = dblarr(40000, 15)
                y_axis = energy_set
                IF keyword_set(interp_terr_plot) THEN BEGIN 
                    x_axis = [indgen(30000)*1.+1, indgen(10000)*1.+30001.] 
                    FOR ien_sort = 0, n_en_sort-1 DO  BEGIN 
                        v(*, ien_sort) = INTERPOL(error(*, ien_sort), t_input(*, ien_sort), x_axis)
                        index = where(x_axis GT energy_set(ien_sort), ct)                        
                        IF ct GT 0 THEN v(index, ien_sort) = !VALUES.F_NAN
                    ENDFOR 
                ENDIF ELSE BEGIN 
                    x_axis = fltarr(ntt*nt*n_en_sort)
                    v = REPLICATE(!VALUES.F_NAN, ntt*nt*n_en_sort, n_en_sort)
                    FOR ien_sort = 0, n_en_sort-1 DO BEGIN 
                        energy_input = energy_set(ien_sort)
                        x_axis((ien_sort*ntt*nt):(ien_sort*ntt*nt+nt*ntt-1)) = t_input(*, ien_sort)
                        v((ien_sort*ntt*nt):(ien_sort*ntt*nt+nt*ntt-1), ien_sort) = error( *, ien_sort)
                    ENDFOR
                    index = sort(x_axis)
                    x_axis = x_axis(index)
                    FOR ien_sort = 0, n_en_sort-1 DO v(*, ien_sort) = v(index, ien_sort)
                ENDELSE 
                
               IF keyword_set(ps) THEN popen, path+year_set(iy)+'_'+n_condi_set(in)+'_3d_Terror.ps', /land
               specplot, x_axis, y_axis, v*100, $
                          no_interp = 1, $
                          lim = { zlog:0, zrange: [0, 100], $
                                  title: year_set(iy)+'    '+n_condi_set(in),  $
                                  xtitle: 'T (eV)', $
                                  ytitle: 'Energy (eV)', $
                                  ztitle: 'T error (%)', $
                                  xrange: [3, 40000.], yrange: [30, 40000.], $
                                  xlog:1, ylog:1, $
                                  XSTYLE:1, ystyle: 1, charsize: 1.2, $
                                  position: [0.12, 0.12, 0.88, 0.9]}   
                
                oplot, [3, 4000.], [30, 40000.]
            ;    oplot, [1, 40000], [453, 453]
                IF keyword_set(ps) THEN pclose
                IF keyword_set(os) THEN stop
            ENDFOR 
        ENDFOR  
ENDIF  
;plot the line plot of Temperature for different year on one plot
    IF keyword_set(line_plot) THEN BEGIN 
    !p.multi = [0, 1, 2]
    FOR in = 0, 1 DO BEGIN 
        n_condi = n_condi_set(in)
        IF keyword_set(ps) THEN popen, path+'different_year_'+n_condi+'.ps', /port ELSE window, in
        plot, [0, 100], [1, 1], yrange = [0.5, 1.5], ystyle = 1, title = n_condi, $
              xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
        oplot, [10, 10], [0.5, 1.5], col = 6
        oplot, [30, 30], [0.5, 1.5], col = 6
        FOR iy = 0, 4 DO BEGIN 
            year = year_set(iy)
            oplot, t_tot_in_set(iy, in, *), t_ratio_tot_set(iy, in, *), color = iy
            xyouts, 80, 1.25+iy*0.05, year, color = iy
        ENDFOR
        plot, [0, 450], [1, 1], yrange = [0.5, 1.5], xstyle = 1, title = n_condi, ystyle = 1, xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
        
        FOR iy = 0, 4 DO BEGIN 
            year = year_set(iy)
            oplot, t_tot_in_set(iy, in, *), t_ratio_tot_set(iy, in, *), color = iy
            xyouts, 300, 1.25+iy*0.05, year, color = iy
        ENDFOR
        IF keyword_set(ps) THEN  pclose
    ENDFOR     

    FOR iy = 0, 4 DO BEGIN 
        year = year_set(iy)
        IF keyword_set(ps) THEN popen, path+'different_density_'+year+'.ps', /port ELSE window, iy
        plot, [0, 100], [1, 1], yrange = [0.5, 1.5], ystyle = 1, title = year, xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
        oplot, [10, 10], [0.5, 1.5], col = 6
        oplot, [30, 30], [0.5, 1.5], col = 6
        FOR in = 0, 1 DO BEGIN 
            n_condi = n_condi_set(in)
            oplot, t_tot_in_set(iy, in, *), t_ratio_tot_set(iy, in, *), color = in+1
            xyouts, 80, 1.25+in*0.05, n_condi, color = in+1
        ENDFOR 
        plot, [0, 450], [1, 1], yrange = [0.5, 1.5],  xstyle = 1, title = year, ystyle = 1, xtitle = 'T input ( eV )', ytitle = 'T_output / T_input'
        
        FOR in = 0, 1 DO BEGIN 
            n_condi = n_condi_set(in)
            oplot, t_tot_in_set(iy, in, *), t_ratio_tot_set(iy, in, *), color = in+1
            xyouts, 300, 1.25+in*0.05, n_condi, color = in+1
        ENDFOR 
        IF keyword_set(ps) THEN  pclose
    ENDFOR

    IF keyword_set(ps) THEN spawn, 'mogrify -format png '+path+'*.ps &'
    !p.multi = 0
ENDIF 
ENDIF     
stop
;----------------------------------------------------------------------
; computing codif moments needed for slice2d_mpe
; Keywords: NAME -> specify other that the default tplot variable name
;           INST_COORD -> Calculate velocities in instrument coordinates
;            RECALC -> Force recalculation instead of reading
;                      pre-processed data
;----------------------------------------------------------------------
IF NOT KEYWORD_SET(sim) THEN BEGIN
    angle = [[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
    energy = [30.0, 40000.0]
    moments = ['V']

    plot_3dmom_from_crib, sat, specie, inst, moments, angle, $
      energy, eff_table, $
      NEW_NAME = 'v_cod',  $
      INST_COORD = 1,    $
      RECALC = 1,        $
      OLD_EFF = 0

                                ; Density
    datadummy = data
    datef = convert_codif_units(datadummy, 'EFLUX', 'codif_ts_eff_corr', $
                                eff_table, packets = 1, $
                                old_eff = old_eff, sat = sat)
    angle = [[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
    energy = [30.0, 40000.0]
    density_real = n_3d_cis(datef, ENERGY = energy, ANGLE = angle)

ENDIF

;----------------------------------------------------------------------
; to run slice2d_mpe, the 3D "data" must be in counts
; Keyword: CUT_BULK_VEL -> get 1-D cut in the bulk velocity frame
;          CUT_PERP     -> the value of vperp to make the 1d cut of vpara
;          CUT_PARA     -> the value of vpara to make the 1d cut of vperp
;----------------------------------------------------------------------
window, /free, ysize = 900
IF specie EQ 3 THEN BEGIN
    range  = [1e-14, 1e-7]
    xrange = [-700, 700]
ENDIF ELSE BEGIN
    range  = [1e-16, 1e-12]
    xrange = [-2700, 2700]
ENDELSE
slice2d_mpe, data, $
             units = 'DIST FUNC', $ ; DIST FUNC, EFLUX, DIFF FLUX
             thebdata = 'B_xyz_codif', $
             vel = 'v_cod', $
                                ;       xrange=xrange, $
             xrange = [-150, 150], $
             nosun = 1, $
             range = range, $
                                ;         onecnt=1, $
             circ = 1, $
                                ;          gsexy = 1, $
                                ;         gsexz = 0, $
                                ;        gseyz = 0, $
             plotenergy = 0, $
             nocross = 0, $
             nosmooth = 0, $
             nosubtract = 1, $
             plotlabel = 0, $
             cut_perp = 38, $
             cut_par = -67, $
             cut_bulk_vel = 0, $
             angle = 20, $
             erange = [40, 40000]

stop

END
