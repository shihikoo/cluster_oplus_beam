pro fit
mass_o = 16*1.6e-27*(1e3)^2/(1.6e-19) ; unit: ev/(km/s)^2
en_set = [31444.7,19398.3, 11966.9, 7382.39, 4554.22, 2809.51, 1733.19, $
          1069.21, 659.599, 406.909, 251.023, 154.857, 95.5315, 58.9337, $
          36.3563]
x=sqrt(2*en_set(0:8)/mass_o)
x_l=sqrt(2*en_set/mass_o)

y_storm=[2.7365799e-10,9.6366303e-10,3.2452863e-09 , 9.5359473e-09,$
         2.9868257e-08,1.0044489e-07,4.9461761e-07,1.5644051e-06,2.9906197e-06]
y_nonstorm=[7.5049506e-11  , 4.3972948e-10,   1.1537785e-09  , 2.5381126e-09 ,  $
            6.8111613e-09,   3.8969664e-08  , 1.7054779e-07,   5.5797452e-07  ,$
            1.4932526e-06]
measure_error_storm= [3.0771076e-10,1.3456066e-09,3.5141782e-09,7.9588580e-09,$
                      2.5614816e-08,1.1052919e-07,3.6387542e-07,7.2192550e-07,$
                      1.0845574e-06]
measure_error_nonstorm=[1.0934163e-10,   4.2826504e-10  , 1.5870872e-09 , $
                        3.5025926e-09,   7.0291425e-09,   3.7047853e-08  ,$
                        1.4889533e-07  , 5.5000534e-07  , 1.4862890e-06]

estimates_storm=[8e-5,50,90]
estimates_nonstorm=[8e-5,-50,90]
storm=1
if storm eq 1 then begin 
    estimates=estimates_nonstorm
    y=y_nonstorm
    measure_error=measure_error_nonstorm
endif else begin 
    estimates=estimates_storm
    y=y_storm
    measure_error=measure_error_storm
endelse 
yfit_para = GAUSSFIT(x, y, co, NTERMS = 3,estimates=estimates,yerror = yerror, sigma = sigma, chisq = chisq,measure_error=measure_error)

plot,en_set,y,xlog=1,ylog=1,xstyle=1,ystyle=1,xrange=[40,4e4],yrange=[1e-11,1e-5]
oplot,en_set,co(0)*exp(-((x_l-co(1))/co(2))^2/2),color=1
print,co,chisq

stop
end 
