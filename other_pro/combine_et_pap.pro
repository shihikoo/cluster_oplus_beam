PRO combine_pap, tail_pap, earth_pap, combine_pap

get_data, tail_pap, data = data, dlim = dlim, lim = lim
time_tpap = data.x
flux_tpap = data.y
pap_tpap = data.v

get_data, earth_pap, data = data
time_epap = data.x
flux_epap = data.y
pap_epap = data.v

time_cpap = time_tpap
flux_cpap = flux_tpap
flux_cpap(*, 0:3) = flux_epap(*, 0:3)
pap_cpap = pap_tpap

pos_1 = STREGEX(tail_pap, '_IN0')
pos_2 = STREGEX(tail_pap, '_AVG')

combine_pap_name = STRMID(tail_pap, 0, pos_1) + $
              STRMID(tail_pap, pos_1+14, pos_2 + 7 - pos_1 - 14 )
str = {x:time_cpap, y:flux_cpap, v:pap_cpap}
store_data, combine_pap_name, data = str, dlim = dlim, lim = lim

END 
