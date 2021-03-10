function sdom, array, sd=sd
; procedure to give the standard deviation, and the standard deviation
; of the mean, of an array.
index=where(finite(array))
array=array(index)
num=n_elements(array)
mean=mean(array,/nan)
delta=array - mean
sd=sqrt(total(delta*delta,/nan)/num) 
sdm=sd/sqrt(num-1.0)
return, sdm
end
