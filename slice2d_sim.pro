PRO slice2d_sim, sim=sim
sim = 1
sat   = 4
specie= 3

time='2002-08-16/11:36:41'
timespan, time, 4, /SEC    ; SECONDS, MINUTES, HOURS, DAYS (DEFAULT)

units_name='Counts'

inst=0       ; 0: CODIF, 1: HIA (this is not supported for the moment)

eff_table=0                     ; 0: GROUND, 1: ONBOARD
;----------------------------------------------------------------------

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
  BKG=0, $
  eff_table, $
  OLD_EFF=0, $
  CNES=0, $
  IC=0

name = 'GLOBE_SC'+string(sat,format='(i1.1)')+$
  '_'+strcompress(units_name, /remove_all)  +$
  '*'+'SP'+string(specie,format='(i1.1)')

tplot_names,name, names=gname

get_data, gname(0), data=data

;----------------------------------------------------------------------
; Fill the globes with sim values
;----------------------------------------------------------------------
IF KEYWORD_SET(sim) THEN BEGIN
;-----------------------------------------------------------------
; Input
;-----------------------------------------------------------------
                ; # / cm^3
;    n0 = 0.01154d               ;+0.0026   ;39
    n0= 0.0096d;+0.0024       ;41
    nn = n0
    n = n0
    loop_times=0
    next:
    
    units_name = 'EFLUX'

;the real velocity
;    v01 = -36.61d; -61.83d  ;6.6519243d           ; GSE km/s
;    v02 = 60.74d; -50.75d  ; -55.826774d          ; GSE km/s
;    v03 = -105.65d; -56.651d ;-33.665597d         ; GSE km/s
; later 4s
;    v01=-43d
;    v02=52d
;    v03=-111d
;;simulated velocity
;    v01 =  -61.83d              ;6.6519243d           ; GSE km/s
;    v02 =  -50.75d               ; -55.826774d          ; GSE km/s
;    v03 =  -56.651d             ;-33.665597d         ; GSE km/s
;41
    V01= -79.2+21.73d
    V02= 8.576-92.8d
    v03= -30.9-90.88d

;use magnetic field
;    T01 = 0.7598d               ;1.6417186d             ; keV
;    T02 = 0.55878d              ;0.007888d             ; keV
;    T03 = 0.55878d              ;0.22775902d             ; keV
;41
    T01=0.062d
    T02=0.1589
    T03=0.1589
;diagnized
;    T01 = 1.6417186d             ; keV
;    T02 = 0.007888d             ; keV
;   T03 =  0.22775902d

    full = 0
    pnoise = 0
    noise = 0.0
 ;-----------------------------------------------------------------
    mp  = 1.67e-27              ; proton mass kg
    pi  = !pi
    K   = 1.38e-23              ; J / Kelvin

    v01 = v01 * 1000.           ; m/s
    v02 = v02 * 1000.           ; m/s
    v03 = v03 * 1000.           ; m/s
                                ; Calculate velocity in instrument coordinates
    datastr = {x:data.time, y:[[v01],[v02],[v03]]}
    store_data, 'dd', data=datastr, $
      dlim={inst_num:0, sens:data.sensitivity, sat:data.sat, $
            phase_instr:data.phase_instr}
    cis_coord_trans, DATA_IN = 'dd', TRANS='GSE->CODIF', $
      DATA_OUT = 'input_vel_instr'

    get_data, 'input_vel_instr', data=in_vel
    v01 = in_vel.y(0)
    v02 = in_vel.y(1)
    v03 = in_vel.y(2)

    IF specie EQ 0 THEN m = mp
    IF specie Eq 3 THEN m = 16 * mp

    n = n * 1e6                 ; # / m^3
    n = n / 8640.               ; CORRECTION FACTOR!!!!
    T01 = 1e3 * T01 * 11600.0   ; Kelvin
    T02 = 1e3 * T02 * 11600.0   ; Kelvin
    T03 = 1e3 * T03 * 11600.0   ; Kelvin

    nbins   = data.nbins
    nenergy = data.nenergy

    corfac = [      4,      4,      4,      4,$
                    2,  2,  2,  2,  2,  2,  2,  2,$
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,$
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,$
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,$
                    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,$
                    2,  2,  2,  2,  2,  2,  2,  2,$
                    4,      4,      4,      4]

                                ; Calculate the f(v)
    J = FLTARR(nenergy,nbins)
    FOR ien = 0, nenergy-1 DO BEGIN
        FOR ibin = 0, nbins-1 DO BEGIN

            theta  = data.theta(ien,ibin)  *!dtor
            phi    = data.phi(ien,ibin)    *!dtor
            energ  = data.energy(ien,ibin)

            u = SQRT(2 * energ * 1.6e-19 / m)

            ux = u*COS(theta)*COS(phi) - v01
            uy = u*COS(theta)*SIN(phi) - v02
            uz = u*SIN(theta)          - v03

            f = (n/SQRT((T01*T02*T03)))*(m/(2*pi*K))^(3./2)*$
              EXP(-m*ux^2/(2*K*T01)-m*uy^2/(2*K*T02)-m*uz^2/(2*K*T03)) ; *corfac(ibin)
            J(ien,ibin) = 0.5 * u^4 * f ; Differential Energy Flux

        ENDFOR
    ENDFOR

                                ; Convert EFLUX in Counts
    data.data = (J)
    data = reverse_convert_codif_units(data,units_name,'codif_ts_eff_corr', $
                                       eff_table, $
                                       specie=specie, $
                                       packets=1, $
                                       old_eff=old_eff, $
                                       incrates=incrates, $
                                       sat=sat)

    data.units_name = 'Counts'

                                ; Randomize and integerize
    IF NOT KEYWORD_SET(full) THEN BEGIN
        data.data = data.data + $
          pnoise * (SQRT(data.data)*(RANDOMU(seed,31,88))) + $
          noise * (2*RANDOMU(seed,31,88)-1)
        ibad = WHERE(data.data LT 1.0, cibad)
        IF cibad GT 0 THEN BEGIN
            data.data(ibad) = 0.0
        ENDIF
        data.data = ROUND(data.data)
    ENDIF

;---------
;Jing: plot the counts data 
;---------

   
    plot3d_options, log = 0
    plot3d_codif, data, zrange = [0, 10];,ebins=[3,4,5,6,7,8,9,10,11,12,13,14]
    stop               
                   


                                ;-------------------------------------------------------------------
                                ; Calculate new density and velocity
                                ;-------------------------------------------------------------------

                                ; Density
    datadummy = data
    datef = convert_codif_units(datadummy,'EFLUX','codif_ts_eff_corr', $
                                eff_table, packets=1, $
                                old_eff=old_eff, sat=sat)
    angle=[[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
    energy=[30.0, 40000.0]
    density_sim = n_3d_cis(datef, ENERGY=energy, ANGLE=angle)

                                ; Velocity In instrument co-ordinates
    angle=[[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
    energy=[30.0, 40000.0]
    vel = compute_velocity(datef, $
                           sat, $
                           inst, $
                           NAME='v_cod', $
                           ENERGY=energy, $
                           ANGLE=angle, $
                           INST_COORD=1)

                                ; Velocity In GSE co-ordinates
    vel = compute_velocity(datef, $
                           sat, $
                           inst, $
                           NAME='v_cod_gse', $
                           ENERGY=energy, $
                           ANGLE=angle, $
                           INST_COORD=0)
    
    get_data, 'v_cod_gse', data=ss
    velocity_sim = ss.y
;   popen,'2005.ps'
;   plot3d_options, log = 1
;   plot3d_codif, data, zrange = [1, 100],title='counts'
;   pclose

    IF density_sim LT n0 THEN BEGIN
        
        nn = 1.1 * nn
        n = nn
        print, nn, density_sim
        loop_times=loop_times+1
        GOTO, next
    ENDIF
    print,loop_times 
      
    stop
ENDIF
;----------------------------------------------------------------------
; computing codif moments needed for slice2d_mpe
; Keywords: NAME -> specify other that the default tplot variable name
;           INST_COORD -> Calculate velocities in instrument coordinates
;            RECALC -> Force recalculation instead of reading
;                      pre-processed data
;----------------------------------------------------------------------
IF NOT KEYWORD_SET(sim) THEN BEGIN
    angle=[[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
    energy=[30.0, 40000.0]
    moments=['V']

    plot_3dmom_from_crib, sat, specie, inst, moments, angle, $
      energy, eff_table, $
      NEW_NAME='v_cod',  $
      INST_COORD = 1,    $
      RECALC = 1,        $
      OLD_EFF = 0

                                ; Density
    datadummy = data
    datef = convert_codif_units(datadummy,'EFLUX','codif_ts_eff_corr', $
                                eff_table, packets=1, $
                                old_eff=old_eff, sat=sat)
    angle=[[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
    energy=[30.0, 40000.0]
    density_real = n_3d_cis(datef, ENERGY=energy, ANGLE=angle)

ENDIF

tt = compute_temperature(data, sat, NAME = 'sim_tem', $
                         ENERGY = energy, ANGLE = angle)

plot3d_options, log = 1
plot3d_codif, data, zrange = [1, 100], $
  title = 'INPUT: T!Lpara!N = '+string(t01/11600, format = '(f6.2)') $
  +'eV   T!Lperp!N = '+ string(t02/11600, format = '(f6.2)')+'eV, ' $
  + string(t03/11600, format = '(f6.2)')+'eV!C'+ $    
  ' SIM  : T!Lpara!N = '+string(tt(0), format = '(f6.2)') $
  +'eV   T!Lperp!N = '+ string(tt(1), format = '(f6.2)')+'eV, ' $
  + string(tt(2), format = '(f6.2)')+'eV'       

stop
;----------------------------------------------------------------------
; to run slice2d_mpe, the 3D "data" must be in counts
; Keyword: CUT_BULK_VEL -> get 1-D cut in the bulk velocity frame
;          CUT_PERP     -> the value of vperp to make the 1d cut of vpara
;          CUT_PARA     -> the value of vpara to make the 1d cut of vperp
;----------------------------------------------------------------------
window, /free, ysize=900
IF specie EQ 3 THEN BEGIN
    range  = [1e-14, 1e-7]
    xrange = [-700,700]
ENDIF ELSE BEGIN
    range  = [1e-16, 1e-12]
    xrange = [-2700,2700]
ENDELSE
slice2d_mpe, data,$
  units = 'DIST FUNC', $        ; DIST FUNC, EFLUX, DIFF FLUX
  thebdata = 'B_xyz_codif', $
  vel = 'v_cod', $
  xrange=xrange, $
  nosun=1, $
  range=range, $
  onecnt=1, $
  circ = 1, $
  gsexy = 1, $
  gsexz = 0, $
  gseyz = 0, $
  plotenergy=0, $
  nocross=0, $
  nosmooth=0, $
  nosubtract=1, $
  plotlabel = 0, $
  cut_perp = 0, $
  cut_par = 0, $
  cut_bulk_vel = 0, $
  angle = 20, $
  erange = [40, 40000]

stop

END
