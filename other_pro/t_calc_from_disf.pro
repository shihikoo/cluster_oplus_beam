PRO t_calc_from_disf, sat = sat, specie = specie, $
  energy = energy, bins = bins, $
  path = path, t_dfit_name = t_dfit_name

IF NOT keyword_set(sat) THEN sat = 4
IF NOT keyword_set(specie) THEN specie = 3
IF NOT keyword_set(energy) THEN energy = [40, 40000]
IF NOT keyword_set(bins) THEN bins = REPLICATE(1, 88)
get_timespan, interv
ts = time_string(interv(0))
date_str = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
time_str = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
;-------------------------------------------
;get the raw data from plot_globe_from_crib
;-------------------------------------------
units_name = 'Counts'
inst = 0    ; 0: CODIF, 1: HIA (this is not supported for the moment)
eff_table = 0       ; 0: GROUND, 1: ONBOARD
plot_globe_from_crib, sat,  specie,   inst,   units_name,   eff_table 
    
name = 'GLOBE_SC'+string(sat, format = '(i1.1)')+$
       '_'+strcompress(units_name, /remove_all)  +$
       '*'+'SP'+string(specie, format = '(i1.1)')
tplot_names, name, names = gname
get_data, gname(0), data = data
;----------------------------------------------------------------------
;find bulk velocity
;----------------------------------------------
angle = [[-90.0, 90.0], [0.0, 360.0]] ;
moments = ['V']
bins_input = bins
erange = energy
   ; stop
plot_3dmom_from_crib, sat, specie, inst, moments, angle, erange, eff_table, $
  NEW_NAME = 'v_cod',   INST_COORD = 1,   RECALC = 1,        $
  bins = bins_input,   use_bins = 1

;----------------------------------------------------------------------
;change bins structrue data into input bins to limit the angular bins range
;-----------------------------------------------------------------
str_element, data, 'bins', REPLICATE(1, 16)#bins, add_replace = 1

IF size(path, /type) EQ 7 THEN BEGIN 
    plot_path = path+'plots/fit_calc/'
    spawn, 'mkdir '+ plot_path
    popen, plot_path+'disf_'+date_str+'_'+time_str+'.ps', /port
ENDIF ELSE window, 1, ysize=900
;stop
erange = energy
slice2d_mpe_fit, data, $
                 thebdata = 'B_xyz_codif', $
                 vel = 'v_cod', $
                 xrange = [-150, 150], $
                 nosun = 1, $        
                 nosubtract = 1, $          
                 erange = erange, $       
                 showdata = 1, $
                 xout1 = xout1, $
                 xout2 = xout2, $
                 xout3 = xout3, $
                 xout4 = xout4, $
                 yout1 = yout1, $
                 yout2 = yout2, $
                 yout3 = yout3, $
                 yout4 = yout4, $
                 resolution = 101

IF size(path, /type) EQ 7 THEN pclose
;-----------------------------------------------------------
;fit the distribution_function
;----------------------------------------
;change to plot screen settigs
oldplot = !p.multi
!p.multi = [0, 2, 1]

; get the data from xout,yout before
x_para = [xout1, xout2]
y_para = [yout1, yout2]
y_para = y_para(sort(x_para))
x_para = x_para(sort(x_para))
x_perp = [xout3, xout4]
y_perp = [yout3, yout4]
y_perp = y_perp(sort(x_perp))
x_perp = x_perp(sort(x_perp))

max_vel_mag = max([max(abs(x_perp)), max(abs(x_para))])

;-----------------calculation ---------------------
;para
max_loc = where(y_para EQ max(y_para))
top_loc = where(y_para GE  max(y_para)/3 $
                AND ABS(x_para-x_para(max_loc(0))) LE 30 $
                , ct_para)
nloc = n_elements(top_loc)
IF nloc GT  3 THEN BEGIN 
    y_para_new = y_para(top_loc)
    x_para_new = x_para(top_loc)
ENDIF   ELSE BEGIN 
    top_loc = where(y_para EQ max(y_para))
    nloc = n_elements(top_loc)
    y_para_new = $
      y_para(((top_loc(0)-2) > 0):((top_loc(nloc-1)+2) < (n_elements(y_para)-1)))
    x_para_new = $
      x_para(((top_loc(0)-2) > 0):((top_loc(nloc-1)+2) < (n_elements(x_para)-1)))
ENDELSE 
yfit_para = GAUSSFIT(x_para_new, y_para_new, co_para, NTERMS = 3, $
                     yerror = yerror_para, sigma = sigma_para, chisq = chisq_para)
;perp
max_loc = where(y_perp EQ max(y_perp))
top_loc = where(y_perp GE  max(y_perp)/3d  $
                AND ABS(x_perp-x_perp(max_loc(0))) LE 30, ct_perp)
nloc = n_elements(top_loc)
IF nloc GT 3 THEN BEGIN 
    y_perp_new = y_perp(top_loc)
    x_perp_new = x_perp(top_loc)
ENDIF ELSE BEGIN 
    top_loc = where(y_perp EQ max(y_perp))
    nloc = n_elements(top_loc)
    y_perp_new = y_perp(((top_loc(0)-2) > 0): ((top_loc(nloc-1)+2) < (n_elements(y_perp)-1)))
    x_perp_new = x_perp(((top_loc(0)-2) > 0):((top_loc(nloc-1)+2) < (n_elements(x_perp)-1)))
ENDELSE 
yfit_perp = GAUSSFIT(x_perp_new, y_perp_new, co_perp, NTERMS = 3, $
                     yerror = yerror_perp, sigma = sigma_perp, chisq = chisq_perp)

; only if you wanna fit with hand input
IF keyword_set(hand_fit) THEN BEGIN 
    yfit_para = GAUSSFIT(x_para, y_para, co_para, NTERMS = 3 $
                         ,  estimate = [1.63e-9, 0, 14] $
                        )
    plot, x_para, y_para, ylog = 0, xrange = [-80, 80]
    oplot, x_para, yfit_para, col = 3
    xyouts, 40, 2e-10, 'Tperp = '+string(0.163*co_para(2)^2)+'eV', col = 2
    yfit_para = GAUSSFIT(x_para, y_para, co_para, NTERMS = 3 $
                         ,  estimate = [1.63e-9, 0, 14] $
                        )
    plot, x_para, y_para, ylog = 0, xrange = [-80, 80]
    oplot, x_para, yfit_para, col = 3
    xyouts, 40, 2e-10, 'Tperp = '+string(0.163*co_perp(2)^2) +'eV', col = 2
ENDIF  

;save the data into string
bad_fit = DBLARR(2, 2)
IF sigma_para(0)/co_para(0) GT 0.1 THEN bad_fit(0, 0) = 1
IF sigma_para(2)/co_para(2) GT 0.1 THEN bad_fit(1, 0) = 1
IF sigma_perp(0)/co_perp(0) GT 0.1 THEN bad_fit(0, 0) = 1
IF sigma_perp(2)/co_perp(2) GT 0.1 THEN bad_fit(1, 0) = 1

co_para(0) = co_para(0)*1e7
co_perp(0) = co_perp(0)*1e7
sigma_para(0) = sigma_para(0)*1e7
sigma_perp(0) = sigma_perp(0)*1e7

t = [0.163*(co_para(2)^2), 0.163*(co_perp(2)^2)]

n = (co_para(0)*co_para(2)*(co_perp(0)^2)*(co_perp(2)^2))*((3.1415926*2)^1.5)

t_error = 0.163*[2*co_para(2)*sigma_para(2), 2*co_perp(2)*sigma_perp(2)]

n_error = sqrt(((co_para(2)*sigma_para(0)*(co_perp(0)^2)*(co_perp(2)^2))^2 $
                +(co_para(0)*sigma_para(2)*(co_perp(0)^2)*(co_perp(2)^2))^2 $
                +(2*co_para(0)*co_para(2)*co_perp(0)*(co_perp(2)^2)*sigma_perp(0))^2 $
                +(2*co_para(0)*co_para(2)*co_perp(2)*(co_perp(0)^2)*sigma_perp(2))^2)) $
          *((3.1415926*2)^1.5)
;stop
chisq = [chisq_para, chisq_perp]
yerror =[yerror_para, yerror_perp]
sigma = [[sigma_para(2)], [sigma_perp(2)]]
fod = [n_elements(x_para_new), n_elements(x_perp_new)]  ; freedome of degree

IF NOT  keyword_set(t_dfit_name) THEN t_dfit_name = 'Dfit_temeperature'
store_data, t_dfit_name, data = {n:n, n_error:n_error, $
                                 t:t, t_error:t_error, $
                                 chisq:chisq, $
                                yerror:yerror, sigma:sigma, $
                                fod:fod}

;----------------------- draw fitting plot ----------------
IF size(path, /type) EQ 7 THEN $
  popen,  plot_path+'dfit_'+date_str+'_'+time_str+'.ps', /land $
ELSE  window, 2, xsize = 900
plot, x_para, y_para, ylog = 0, psym = -1, $
      xrange = [-100, 100], xstyle = 1, $
;     xrange = [-max_vel_mag, max_vel_mag], $
      title = '                                                 '+ $
      'Distribution Function Fit for           ' $
      +time_string(interv(0))+'   to   '+ time_string(interv(1))
oplot, x_para_new, y_para_new, psym = 1, col = 6
oplot, x_para_new, yfit_para, col = 2, psym = -1
oplot, [-100, 100], [max(y_para)/2, max(y_para)/2]
xyouts, -70, max(y_para), 'Tpara = '+STRCOMPRESS(t(0), /remove_all)+ $
        ' +- '+STRCOMPRESS(t_error(0), /remove_all) +' eV'
xyouts, -70, max(y_para)*1.1, 'n = '+STRCOMPRESS(n, /remove_all)+ $
        ' +- '+STRCOMPRESS(n_error, /remove_all)+' cm-3'
;--------------------------- draw perp----------------------
plot, x_perp, y_perp, ylog = 0, psym = -1, $
; xrange = [-max_vel_mag, max_vel_mag] , $
      xrange = [-100, 100], xstyle = 1
oplot, x_perp_new, y_perp_new, psym = 1, col = 6
oplot, x_perp_new, yfit_perp, col = 2, psym = -1
oplot, [-100, 100], [max(y_perp)/2, max(y_perp)/2]
xyouts, -70, max(y_perp), 'Tperp = '+string(t(1))+ '  +-'+string(t_error(1)) +'eV'
IF size(path, /type)  EQ 7 THEN pclose

;-back to the normal screen set
!p.multi = oldplot
;stop
END 
