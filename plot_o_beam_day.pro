;Purpose O+ beam over time interval (input) with displaytime 
;         for every displaytime (input)
;
;Input: sc           : Cluster no. if not set the default is 4
;       calc_time    : the calculate time for each loop
;       time_start   : in idl time format
;       time_end     : in idl time format 
;       displaytime  : the displaytime for single plot
;
;Keywords:1
;       average_time : time for averaging data, if not set the default is 5 min
;       idl_plot     : plot the result plot in idl_window
;       ps           : plot the result plot in ps files,
;       dumpdata     : output data file
;       globe_plot   : plot a set of globe plot to show the selected
;                      range for plotting mom
;       store_data   : store_data into .tplot      default: 1 
;       plot_mom     : plot mom of the identifed O+ beam
;
;Output: Depends on Keywords settings 
;        There will also be two .log files
;
;Written by Jing Liao  02/20/2008
;
PRO plot_o_beam_day, time_start = time_start, time_end= time_end,stop=stop,sc=sc,sp=sp
  if not keyword_set(sc) then sc = 4 ;set the satallite number 
  sc_str = STRING(sc, FORMAT = '(i1.1)')
  
  if not keyword_set(sp) then sp=3
  
  inst_input = 0                ;set the instrument number  
  IF NOT keyword_set(time_start) THEN  time_start = '2002-09-10/00:00:00'
  IF NOT keyword_set(time_end) THEN time_end = '2002-09-10/23:59:59'
  dt= time_double(time_end)-time_double(time_start)
  calc_time = 24.*60. *60. < dt  ;in seconds
;-----------------------------------------
;Set the keywords used in find_o_beam.pro
;------------------------------------------
  beam_recalc = 1
  mom_recalc = 1
  find_phase = 1
  add_imf = 1                   ;0
  add_eflux = 1  &  plot_add_eflux_procedure = 0
  add_anodes = 1 
  add_distfunc = 1 & plot_add_distfunc_procedure = 0
  
  idl_plot = 1
  ps = 1
  
  dumpdata = 1
  plot_imf = 1
  plot_mom = 1
  
  store_data = 1
  
  dfit_temperature = 0
  
  use_angle_range = 1
  beam_angle_range = 22.5
  use_energy_range = 1 
  
  globe_plot = 0
  plot_lowcount_filter = 0
  plot_sw_sheath = 0
  
  flux_threshold=[0,0,0]
  
  use_hiabeta = 0
  only_in_lobe = 0
  if sc eq 1 then begin 
     path = 'output/o_beam/test_year_dependence_with_HIA/sc1_hiabeta_lobe/'
     use_hiabeta = 1
     only_in_lobe = 1
     plot_imf=0 
     plot_mom=0
  endif 
  if sc eq 4 then path = 'output/'
;test_flux_threshold/'+strcompress(flux_threshold(0),/remove_all)+'/' 
  spawn, 'mkdir '+path
  spawn, 'mkdir '+path+'tplot_restore/'
  spawn, 'mkdir '+path+'data/'
  spawn, 'mkdir '+path+'plots/'
  spawn, 'mkdir '+path+'plots/mom_with_neg'
  spawn, 'mkdir '+path+'log/'
  
  display_time = 6.*60*60 < dt
  average_time = 5 * 60         ;in seconds 
;-----------------------------------------------
;Write [START] in log files
;-----------------------------------------------
  OPENU, unit, path+'log_plotted.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[START]'
  FREE_LUN, unit         
  
  OPENU, unit, path+'log_errors.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[START]'
  FREE_LUN, unit     
;---------------------------------------------------
; Set the loop as requested
;---------------------------------------------------
  ts = time_double(time_start)
  te = time_double(time_end)
  ntime = CEIL((te - ts)/calc_time) 
;------------------------------------------------------------ 
  FOR i = 0l, ntime-1 DO BEGIN  
; Timespan over each displaytime
     timespan, ts + i*calc_time, calc_time, /seconds
;some days cannot be plotted because of the orbit or else
     time_useless = [['2002-06-11/00:00:00', '2002-06-14/00:00:00'], $
                     ['2005-10-10/00:00:00', '2005-10-11/00:00:00'], $
                     ['2006-03-23/00:00:00', '2006-03-24/00:00:00'], $
                     ['2009-09-22/00:00:00', '2009-09-23/00:00:00'],$
                     ['2009-09-26/00:00:00', '2009-09-27//00:00:00']]   
     nouse = 0
     FOR j = 0, N_ELEMENTS(time_useless(0, *))-1 DO BEGIN 
        nouse = nouse+(ts+i*calc_time  GE  time_double(time_useless(0, j)) AND  $
                       ts+(i+1)*calc_time LE  time_double(time_useless(1, j))) 
     ENDFOR 
     
     IF nouse EQ 0 THEN BEGIN 
; identify O+ beam plot
        find_o_beam, sc = sc, $
                     sp=sp,$
                     average_time = average_time, $ 
                     idl_plot = idl_plot, $
                     ps = ps, $
                     dumpdata = dumpdata, $
                     globe_plot = globe_plot, $
                     store_data = store_data, $
                     plot_mom = plot_mom, $
                     path = path, $
                     plot_lowcount_filter =  plot_lowcount_filter, $      
                     beam_recalc = beam_recalc, $
                     mom_recalc = mom_recalc,  $
                     find_phase = find_phase, $
                     displaytime = display_time, $
                     add_imf = add_imf, $
                     plot_sw_sheath = plot_sw_sheath, $
                     use_angle_range = use_angle_range, $
                     use_energy_range = use_energy_range, $
                     plot_imf = plot_imf, $
                     dfit_temperature = dfit_temperature, $
                     inst_input = inst_input, $
                     beam_angle_range =  beam_angle_range, $
                     add_eflux = add_eflux, $
                     add_anodes = add_anodes, $
                     use_hiabeta = use_hiabeta,$
                     only_in_lobe = only_in_lobe,$
                     add_distfunc = add_distfunc,$
                     plot_add_distfunc_procedure = plot_add_distfunc_procedure,$
                     flux_threshold=flux_threshold
     ENDIF 
     if keyword_set(stop) then stop
  ENDFOR   
;------------------------------------------------------------
;----write [END] in log files
  OPENU, unit, path+'log_plotted.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[END]'
  FREE_LUN, unit         
  
  OPENU, unit, path+'log_errors.txt', /GET_LUN, /APPEND
  PRINTF, unit, SYSTIME(), '[END]'
  FREE_LUN, unit  
  
  print, time_start, time_end
  stop
END 
