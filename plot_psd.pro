pro plot_psd,ps=ps

sat   = 4
specie= 3
erange=[40,40000]
;erange=[251, 660]
units_name='Counts'
inst=0       ; 0: CODIF, 1: HIA (this is not supported for the moment)
eff_table=0                     ; 0: GROUND, 1: ONBOARD

time_start='2002-09-11/09:00:00'
time_end='2002-09-11/13:00:00'
nloop=(time_double(time_end)-time_double(time_start))/4.
plot_path='thesis_figure/plots/psd/'
for ii=14, nloop-1 do begin 

    time=time_double(time_start)+4.*ii
    timespan, time, 8, /sec  ; SECONDS, MINUTES, HOURS, DAYS (DEFAULT)
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
    angle=[[-90.0, 90.0], [0.0, 360.0]] ; bin range to sum over
    energy=erange
    moments=['V']

    plot_3dmom_from_crib, sat, specie, inst, moments, angle, $
      energy, eff_table, $
      NEW_NAME='v_cod',  $
      INST_COORD = 1,    $
      RECALC = 1,        $
      OLD_EFF = 0 
    if total(data.data,/nan) ne 0 then begin 
    if keyword_set(ps) then popen, plot_path+strcompress(ii,/remove_all)+'.ps' else window, ysize=900

    slice2d_mpe, data,$
      units = 'DIST FUNC', $
      thebdata = 'B_xyz_codif', $
      vel = 'v_cod', $
      xrange=[-600,600], $
      nosun=1, $
                                ;      range=[max(data.data)/1000., max(data.data)], $
    onecnt=0, $
      ang = 0, $
      gsexy = 0, $
      gsexz = 0, $
      gseyz = 0, $
      plotenergy=0, $
      nocross=0, $
      nosmooth=1, $
      nosubtract=1, $
      noolines=1, $
      finished=0, $
      plotlabel = 0, $
      cut_perp = 0, $
      cut_par = 0, $
      cut_bulk_vel = 1, $
      erange = erange, $
      double=0,$
      showdata=0,$	
      nozlog=1
    if keyword_set(ps) then pclose else stop 
endif 
endfor 
end 
