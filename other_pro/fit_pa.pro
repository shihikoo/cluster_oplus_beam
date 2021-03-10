              ; fit the plot with sin function and overplot on the original plot
                       
                  fit_coef = [MAX(y_data), 1] ; y = coef(0) * sin( coef(1) * x)
                 weight =1./(TOTAL((y_data(0:7))^2)/N_ELEMENTS(y_data(0:7)) $ 
                               -(TOTAL(y_data(0:7))/N_ELEMENTS(y_data(0:7)))^2)
                  
                  fit = CURVEFIT(x_data, y_data, weight, fit_coef, $
                                 FUNCTION_NAME = 'sinfunc', /NODERIVATIVE)

           ;       PRINT, 'fit func: y=' + STRCOMPRESS(fit_coef(0), /REMOVE_ALL) $ 
           ;              + 'sin(' + STRCOMPRESS(fit_coef(1), /REMOVE_ALL) + ' * x )'

            ;      y_fit = FLTARR(N_ELEMENTS(y_data))
             ;     test_coef = fit_coef
              ;    test_pa = pa
               ;   sinfunc, x_data, fit_coef, y_fit
                ;  OPLOT, x_data, y_fit

;stop
                 ;plot them on file
 ;                 dd = time_string(time)
  ;                fln = 'plots/pa_' + STRMID(dd, 0, 4) + STRMID(dd, 5, 2) $
   ;                     + STRMID(dd, 8, 2) + '_' + STRMID(dd, 11, 2) $
    ;                    + STRMID(dd, 14, 2) + STRMID(dd, 17, 2) + '.ps'
                  
     ;             popen, fln
                  
      ;            PLOT, x_data, y_data, xlog = 0, ylog = 0, xrange = [0, 180],$
       ;                 yrange = [0, 1600], psym = -2, xtitle = 'PITCH ANGLE', $ 
        ;                ytitle = 'Diff Flux', $
         ;               title = 'ENERGY RANGE:'+ e_min+'_'+ e_max +'!Ctimespan' $
          ;              + time_string(time)+' to '+time_string(time+average_time)
                        
          ;        tplot, pa_name                               
               
           ;       pclose  
