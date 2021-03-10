PRO combine_beam_judge, names_0, names_1

get_data, names_0, data = data, dlim = dlim, lim = lim
time_one = data.x
data_one = data.y

get_data, names_1, data = data
time_two = data.x
data_two = data.y

time_combine = [time_one, time_two]
data_combine = [data_one, data_one]

data_combine = data_combine(sort(time_combine))
time_combine = time_combine(sort(time_combine))

combine_name = 'PASPEC_SC4_DIFF_UNDIFFFLUX_SP3_ET0_All_AVG'+'600'$
              +'_BEAMJUDGEMENT'
str = {x:time_combine, y:data_combine}
store_data, combine_name, data = str, dlim = {psym:1}

END 
