;Purpose: Find O+ beam over time interval (input) with displaytime 
;         for every displaytime (input)
;
;Input: sc           : Cluster no. if not set the default is 4
;       calc_time    : the calculate time for each loop
;       time_start   : in idl time format
;       time_end     : in idl time format 
;       displaytime  : the displaytime for single plot
;
;Keywords:
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
PRO plot_o_beam_test, time_start = time_start, time_end= time_end
;20020611,12,13
sc = 4                          ;set the satallite number 
sc_str = STRING(sc, FORMAT = '(i1.1)')

IF NOT keyword_set(time_start) THEN  time_start = '2002-10-01/06:00:00'
IF NOT keyword_set(time_end) THEN time_end = '2002-10-01/23:59:59'
calc_time = 2.*60. *60.        ;in seconds
;-----------------------------------------
;Set the keywords used in find_o_beam.pro
;------------------------------------------
beam_recalc = 0
mom_recalc = 1
find_phase = 0
add_imf = 0

plot_imf = 0
plot_mom = 0
idl_plot = 0
ps = 0
dumpdata = 0
store_data = 1

globe_plot = 0
plot_lowcount_filter = 0
dfit_temperature = 1

use_angle_range = 0
use_energy_range = 0

path = 'test_mom/limit_no/' 
;spawn,  'mkdir '+path
;spawn, 'mkdir '+path+'tplot_restore/'
;spawn, 'mkdir '+path+'data/'
;spawn, 'mkdir '+path+'plots/'
;spawn, 'mkdir '+path+'plots/mom_with_neg/'

display_time = 2.*60*60
average_time = 5 * 60           ;in seconds  
beam_angle_range = 11.25

dont_plot_sw_sheath = 1
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
                    ['2006-03-23/00:00:00', '2006-03-24/00:00:00']]              
    nouse = 0
    FOR j = 0, N_ELEMENTS(time_useless(0, *))-1 DO BEGIN 
        nouse = nouse+(ts+i*calc_time  GE  time_double(time_useless(0, j)) AND  $
                       ts+(i+1)*calc_time LE  time_double(time_useless(1, j))) 
    ENDFOR 
    
    IF nouse EQ 0 THEN BEGIN 
; identify O+ beam plot
        find_o_beam, sc = sc, $
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
                     dont_plot_sw_sheath = dont_plot_sw_sheath, $
                     use_angle_range = use_angle_range, $
                     use_energy_range = use_energy_range, $
                     plot_imf = plot_imf, $
                     beam_angle_range =  beam_angle_range, $
                     dfit_temperature = dfit_temperature
    ENDIF 

t_name='TDMOM_ENVARIOUS_SC4_PHI90_270_MTTEMPERATURE_SP3_ET0_All_T_AVG300'
fit_name='TDMOM_ENVARIOUS_SC4_PHI90_270_MTTEMPERATURE_SP3_ET0_All_T_dfit_AVG300'
ylim, fit_name, 1, 100
options,t_name,'color', 3

;popen, path+'temperature_compare.ps', /land
;tplot, fit_name, title = 'angular angle: 33.75'
tplot_panel,v=fit_name,o=t_name, psym = 1
xyouts, 3600, 50, 'T (fit)', size = 2
xyouts, 3600, 70, 'T (3dmom)', color = 3, size = 2
;pclose

stop
ENDFOR  
;------------------------------------------------------------
;----write [END] in log files
OPENU, unit, path+'log_plotted.txt', /GET_LUN, /APPEND
PRINTF, unit, SYSTIME(), '[END]'
FREE_LUN, unit         

OPENU, unit, path+'log_errors.txt', /GET_LUN, /APPEND
PRINTF, unit, SYSTIME(), '[END]'
FREE_LUN, unit  

;stop
END 
