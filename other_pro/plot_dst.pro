PRO plot_dst

fln = 'storm_phase.dat'

prestorm_start = DBLARR(300)
storm_onset = DBLARR(300) 
min_dst = DBLARR(300)
recovery_end = DBLARR(300) 

OPENR, unit, fln, /GET_LUN
jj = 0
dummy = ''
WHILE NOT EOF(unit) DO BEGIN
    READF, unit, dummy
    IF STRMID(dummy, 0, 5) EQ 'Start' THEN BEGIN 
        title = dummy 
    ENDIF ELSE BEGIN 
        prestorm_start(jj) = time_double(STRMID(dummy, 0, 20))
        storm_onset(jj) = time_double(STRMID(dummy, 20, 20))
        min_dst(jj) =  time_double(STRMID(dummy, 40, 24))
        recovery_end(jj) =  time_double(STRMID(dummy, 60, 24))

        jj = jj+1      
    ENDELSE 
ENDWHILE
CLOSE, unit, /all
nstorm = jj

prestorm_start = prestorm_start(0:nstorm-1)
storm_onset = storm_onset(0:nstorm-1)
min_dst = min_dst(0:nstorm-1)
recovery_end = recovery_end (0:nstorm-1)

FOR iy = 6, 6 DO BEGIN  
    FOR im = 1, 5 DO BEGIN 
        year = string(2000+iy, format = '(i4.4)')
        month = string(im, format = '(i2.2)')
        next_month = string(im+1,  format = '(i2.2)')
        st = time_double(year+'-'+month+'-'+'01/00:00:00')
        et = time_double(year+'-'+next_month+'-'+'01/00:00:00')
        dt = et-st              ;in seconds
        
        index = where(min_dst GE st AND min_dst LE  et, ct)
        
        timespan, st-3*24.*3600., dt+2*24.*3600., /SECONDS
        read_omni          
        
        ylim, 'Bz_gse', -30, 30, 0
        ylim, 'flow_pressure', 0.1, 100, 1 
        ylim, 'Dst_Index', -300, 100, 0 
        
        tplot, ['Bz_gse', 'flow_pressure', 'Dst_Index']
        yline, 'Bz_gse', col = 3
        yline, 'Dst_Index', col = 3 &  yline, 'Dst_Index', offset = -50, col = 3
        IF ct GT 0 THEN BEGIN 
            timebar, prestorm_start(index), color = 4
            timebar, storm_onset(index), color = 6
            timebar, min_dst(index), color = 1
            timebar, recovery_end(index), color = 2
        ENDIF 
stop
        popen, 'output_storm/plots/storm_dst/Dst_'+year+'_'+month+'.ps', /land       
        tplot, ['Bz_gse', 'flow_pressure', 'Dst_Index']
        yline, 'Bz_gse'
        yline, 'Dst_Index', col = 3
        yline, 'Dst_Index', col = 3 & yline,'Dst_Index', offset = -50, col = 3
        IF ct GT 0 THEN BEGIN
            timebar, prestorm_start(index), color = 4
            timebar, storm_onset(index), color = 6
            timebar, min_dst(index), color = 1
            timebar, recovery_end(index), color = 2
        ENDIF                 
        pclose
        stop    
    ENDFOR 
ENDFOR 
stop
END 
