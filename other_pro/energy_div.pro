PRO energy_div, en_name

get_data, en_name, data = data, dlim = dlim, lim = lim
time = data.x
flux = data.y
energy = data.v

flux = data.y(*, 0:14)
energy = data.v(*, 0:14)
str = {x:time, y: flux, v:energy}
store_data, en_name, data = str, dlim = dlim, lim = lim

flux_h = data.y(*, 0:3)
energy_h = data.v(*, 0:3)
str = {x:time, y:flux_h, v:energy_h}
store_data, en_name + '_h', data = str, dlim = dlim, lim = lim
         
flux_m = data.y(*, 4:10)
energy_m = data.v(*, 4:10)
str = {x:time, y:flux_m, v:energy_m}
store_data, en_name + '_m', data = str, dlim = dlim, lim = lim

flux_l = data.y(*, 11:14)
energy_l = data.v(*, 11:14)
str = {x:time, y:flux_l, v:energy_l}
store_data, en_name + '_l', data = str, dlim = dlim, lim = lim

END 
