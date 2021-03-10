PRO globe_plot, sat = sat, specie = specie, units_name = units_name, $
                bins = bins, energy = energy, path = path, ps = ps

IF NOT KEYWORD_SET(sat) THEN sat = 4
IF NOT KEYWORD_SET(specie) THEN specie = 3
IF NOT KEYWORD_SET(units_name) THEN units_name='Counts'
IF NOT KEYWORD_SET(bins) THEN bins = REPLICATE(1, 16, 88)
IF NOT KEYWORD_SET(energy) THEN energy = [40., 40000.]
inst=0 ; 0: CODIF, 1: HIA (this is not supported for the moment)
eff_table=0 ; 0: GROUND, 1: ONBOARD
;----------------------------------------------------------------------
; Load/calculate the globe plot data
plot_globe_from_crib, sat, specie, inst, units_name, eff_table

name = 'GLOBE_SC' + string(sat, format = '(i1.1)') + $
       '_' + strcompress(units_name, /remove_all) + $
       '*' + 'SP'+string(specie,format='(i1.1)')
tplot_names, name, names=names
get_data, names(0), data=dat

plot3d_options, log = 1
get_timespan, in
ts = time_string(in(0))
date_str = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
time_str = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)

fln = path+'plots/globe_plot/'
spawn, 'mkdir '+fln
IF keyword_set(ps)THEN popen, fln+'globe_plot_'+date_str+'_'+time_str+'_original.ps', /land ELSE window, /free
plot3d_codif, dat, zrange=[1,100]
IF keyword_set(ps) THEN pclose

ebins = replicate(0, 16)
er = [energy_to_ebin(dat, energy)]
IF er(0) GT er(1) THEN er = reverse(er)
ebins(er(0):er(1)) = 1

str_element, dat, 'data', dat.data*(ebins#bins), add_replace = 1
IF keyword_set(ps) THEN popen, fln+'globe_plot_'+date_str+'_'+time_str+'_limited.ps', /land $
ELSE window, /free
plot3d_codif, dat, zrange=[1,100]
IF keyword_set(ps) THEN pclose

str_element, dat, 'data', ebins#bins, add_replace = 1
IF keyword_set(ps) THEN popen, fln+'globe_plot_'+date_str+'_'+time_str+'_limits.ps', /land $
ELSE window, /free
plot3d_codif, dat, zrange=[0.1,1]
IF keyword_set(ps) THEN pclose
;stop
END 
