;Purpose: Find O+ beam during storm time(three perigee interval around min Dst)
;
;Input: sc           : Cluster no. if not set the default is 4
;       calc_time    : the calculate time for each loop
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

PRO plot_o_beam_storm

sc = 4                          ;set the satallite number 
sc_str = STRING(sc, FORMAT = '(i1.1)')

calc_time =  6. * 3600.         ;in sec
;-----------------------------------------
;Set the keywords used in find_o_beam.pro
;------------------------------------------
beam_recalc = 1
mom_recalc = 0
find_phase = 1
index = 1

plot_mom = 0
idl_plot = 0
ps = 1
dumpdata = 0
store_data = 1

globe_plot = 0
plot_lowcount_flter = 0       

display_time = 6.*60*60
average_time = 5 * 60           ;in seconds  
path = 'output_storm/' 
;-------------------------------------------------------------------
;Read storm phase(prestorm, storm time and recovery time) form 
;file 'storm_phase.dat' and store them into arrays
;-------------------------------------------------------------------
OPENR, unit, 'storm_phase.dat', /GET_LUN
prestorm_start = DBLARR(300)
storm_onset = DBLARR(300) 
min_dst = DBLARR(300)
recovery_end = DBLARR(300) 
jj = 0l
dummy = ''
WHILE NOT EOF(unit) DO BEGIN
    
    READF, unit, dummy
    IF STRMID(dummy, 0, 5) NE 'Start' THEN BEGIN 
        prestorm_start(jj) = time_double(STRMID(dummy, 0, 20))
        storm_onset(jj) = time_double(STRMID(dummy, 20, 20))
        min_dst(jj) =  time_double(STRMID(dummy, 40, 24))
        recovery_end(jj) =  time_double(STRMID(dummy, 60, 24))
        
        jj = jj+1      
    ENDIF  
ENDWHILE
CLOSE, unit, /all
nstorm = jj

prestorm_start = prestorm_start(0:nstorm-1)
storm_onset = storm_onset(0:nstorm-1)
min_dst = min_dst(0:nstorm-1)
recovery_end = recovery_end (0:nstorm-1)
;---------------------------------------------------------------------
; Read CLUSTER perigee times list from file: sc4_perigee_times.dat
; (different list for each S/C)
; Store these times in variable petime (in tplot time format)
;-------------------------------------------------------------------
OPENR, unit, 'sc' + sc_str + '_perigee_times.dat', /GET_LUN
petime = DBLARR(3000)           ; assumes no more than 3000 perigee passes
dummy = ''
jj = 0l
WHILE NOT EOF(unit) DO BEGIN
    READF, unit, dummy
    petime(jj) = time_double(dummy)
    jj = jj + 1
ENDWHILE
petime = petime(0:jj-1)
CLOSE, unit

first_pe = DBLARR(jj)
second_pe = DBLARR(jj)
third_pe = DBLARR(jj) 
forth_pe = DBLARR(jj)

;----------------------------------------------------
;Write [START] in log files
;-----------------------------------------------------
OPENU, unit, path+'log_plotted.txt', /GET_LUN, /APPEND
PRINTF, unit, SYSTIME(), '[START]'
FREE_LUN, unit         

OPENU, unit, path+'log_errors.txt', /GET_LUN, /APPEND
PRINTF, unit, SYSTIME(), '[START]'
FREE_LUN, unit      

;---------------------------------------------------------------------
; For each minimum Dst time choose the time intarval of the plots.
; The time interval is selected is of three orbits duration
FOR  ii = 6, N_ELEMENTS(min_dst)-1  DO BEGIN 
    
; find the time of the four perigee passes closest (2 before and 
; 2 after) to the minimum Dst
    before_index = where ((petime - min_dst(ii)) < 0, counts)
    after_index = where ((petime - min_dst(ii)) > 0)          
    before_pe = petime(before_index(counts-1))
    after_pe = petime(after_index(0))
    
; time interval start time, end time and dt
    first_pe(ii) =  petime(before_index(counts-2))
    second_pe(ii) =  petime(before_index(counts-1))
    third_pe(ii) =  petime(after_index(0))
    forth_pe(ii) =  petime(after_index(1))

    time_start = petime(before_index(counts-2))
    time_end = petime(after_index(1))
    dt = time_end - time_start      
    idt = CEIL(dt / calc_time)
    
;------------------------------------------------------
;Run for every single plot displaytime over the whole time interval 
    FOR kk = 0, idt-1  DO BEGIN   
        PRINT, STRING(ii) + '   --' + STRING(kk)

; Timespan over each displaytime
        t_s = time_start + kk * calc_time
        t_e = t_s + calc_time
        t_dt = t_e - t_s
        timespan, t_s, t_dt, /SECONDS  
;identify O+ beam plot 
        find_o_beam, sc = sc, $
                     average_time = average_time, $ 
                     idl_plot = idl_plot, $
                     ps = ps, $
                     dumpdata = dumpdata, $
                     globe_plot = globe_plot, $
                     store_data = store_data, $
                     plot_mom = plot_mom, $
                     path = path, $
                     plot_lowcount_flter =  plot_lowcount_filter, $      
                     beam_recalc = beam_recalc, $
                     mom_recalc = mom_recalc,  $
                     find_phase = find_phase, $
                     displaytime = display_time, $
                     index = index

        OPENU, unit, path+'log.txt', /GET_LUN, /APPEND
        PRINTF, unit, SYSTIME()+' ['+STRING(ii, FORMAT = '(i2.2)')+ $
                ',' +STRING(kk, FORMAT = '(i2.2)') +'] - '+ $
                time_string(t_s) +' - '+time_string(t_e)
        FREE_LUN, unit     
    ENDFOR 
;-------------------------------------------------------
ENDFOR  
;------------------------------------------------------------------

;----write [END] in log files
OPENU, unit, path+'log_plotted.txt', /GET_LUN, /APPEND
PRINTF, unit, SYSTIME(), '[END]'
FREE_LUN, unit         

OPENU, unit, path+'log_errors.txt', /GET_LUN, /APPEND
PRINTF, unit, SYSTIME(), '[END]'
FREE_LUN, unit     
END 
