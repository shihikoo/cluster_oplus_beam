; Purpose: Test the effect when change units in 
; espetially, diff flux and dist func
; also check how they works under plot_en routine and plot_pa routine
; written by Jing Liao 12/02/2010
pro units_test
timespan,'2002-07-29/06:48:09',12,/sec
sc=4
inst_input=0 ;'CODIF'
;-- Load energy spectra - tailward (or earthward for HIA)--
sat = [sc]
specie = [3]
angle = [[-90, 90], [90, 270]]
inst = inst_input & units_name = 'DIFF FLUX' 
eff_table = 0

plot_en_spec_from_crib, sat, specie, inst,  units_name, angle, eff_table, recalc = 1
zlim,'*DIFFFLUX*',1,10,1
;sat = [sc]
;specie = [3]
;angle = [[-90, 90], [90, 270]]
;inst = inst_input & units_name = 'Counts' & eff_table = 0

;plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 1
;zlim,'*COUNTS*',1,10,1

sat = [sc]
specie = [3]
angle = [[-90, 90], [90, 270]]
inst = inst_input & units_name = 'DIST FUNC' & eff_table = 0

plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 1
zlim,'*DISTFUNC*',1e-11,1e-7

sat = sc
specie = 3
inst = inst_input & units_name = 'DIST FUNC' & eff_table = 0
plot_globe_from_crib,sat,specie,inst,units_name,eff_table

get_data,1, data=dd1 ;diff flux
get_data,2,data=dd2  ;dist func
get_data,'GLOBE_SC4_DISTFUNC_PR47_SP3',data=dd3

print,dd1.y(2,10)*(dd3.mass)^2/2e5/dd1.v(2,10)  ;flux*mass^2/2.e5/energy
print,dd2.y(2,10)
;tplot,'*'



stop
end
