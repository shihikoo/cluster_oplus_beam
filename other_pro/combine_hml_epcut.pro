PRO combine_hml_epcut, cut_all, cut_h, cut_m, cut_l, cut_c

get_data, cut_h, data = data, dlim = dlim, lim = lim
time_h = data.x
energy_h = data.y

get_data, cut_m, data = data
time_m = data.x
energy_m = data.y

get_data, cut_l, data = data
time_l = data.x
energy_l = data.y

time_nan = !VALUES.F_NAN 
energy_nan = !VALUES.F_NAN 

time_c = [time_h,time_nan,time_m,time_nan, time_l]
energy_c = [energy_h, energy_nan, energy_m,energy_nan, energy_l]

str = {x:time_c , y:energy_c }
cut_c = cut_all+'_c'
store_data, cut_c, data = str, dlim = dlim, lim = lim

;tplot, cut_c
;stop
END 
