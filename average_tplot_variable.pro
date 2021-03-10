;+
; PROCEDURE: average_tplot_variable
;
; PURPOSE: to time average any tplot variable
;
; INPUT:
;        var -> variable name - Array of strings or numbers. The '*'
;                               is also accepted.
;        avg -> averaging time in seconds
;
; KEYWORD:
;         NEW_NAME -> the averaged variable is stored with a new name that
;                     includes the avereging time
;
; CREATED by: C. Mouikis
;
; LAST MODIFIED: 09/25/01
;
; MODIFICATION HISTORY:
;   09/25/01 - The time step loop index is declared as a long.
;   06/11/08 - Keyword total added            
;-
PRO average_tplot_variable, var, avg, NEW_NAME = NEW_NAME, sumup = sumup
  
COMMON get_error, get_err_no, get_err_msg, default_verbose
  
get_err_no = 0 & get_err_msg = ''

  IF avg LE 0 THEN BEGIN
    PRINT, 'the avg time must be > 0'
    RETURN
  ENDIF
  
  ;--------------------------------------------------------------------
  ; Find the names of the variables to be averaged
  ;--------------------------------------------------------------------
  IF size(var, /TYPE) EQ 7 THEN BEGIN ; var -> string
    IF var EQ '*' THEN BEGIN
      tplot_names, NAMES=var_names
    ENDIF ELSE BEGIN
      var_names=var
    ENDELSE
  ENDIF ELSE BEGIN                    ; var -> integer
    tplot_names, var, NAMES=var_names
  ENDELSE
  
  ;--------------------------------------------------------------------
  ; Loop through all variables
  ;--------------------------------------------------------------------
  FOR iv = 0, n_elements(var_names)-1 DO BEGIN
    
    ;------------------------------------------------------------------
    ; Already averaged variables are excluded
    ;------------------------------------------------------------------
    IF strpos(var_names(iv), '_AVG') NE -1 THEN GOTO, next
    
    ;------------------------------------------------------------------
    ; Extract specie and s/c information from the variable name
    ;------------------------------------------------------------------
    IF strpos(var_names(iv), '_SP') NE -1 THEN BEGIN
      specie = strmid(var_names(iv), strpos(var_names(iv), '_SP') + 3, 1)
    ENDIF
    sat    = strmid(var_names(iv), strpos(var_names(iv), '_SC') + 3, 1)

    ;------------------------------------------------------------------
    ; Extract data information from tplot variable
    ;------------------------------------------------------------------
    get_data, var_names(iv), data=data, dlim=dlim, lim=lim
    
    str_element, data, 'V', value, SUCCESS=vfound
    
    time = data.x
    yray = data.y
    
    ;------------------------------------------------------------------
    ; Set variables/arrays
    ;------------------------------------------------------------------
    specie_str = ['H!U+!N','He!U++!N','He!U+!N','O!U+!N']

    dt = time(n_elements(time)-1) - time(0) ; time interval in seconds
    n_avg = floor(dt/avg) 
    
    IF n_avg LE 3 THEN GOTO, next

    yray_size = size(yray)
    dim_numb = yray_size(0)
    
    time_avg = dblarr(n_avg)
    IF dim_numb EQ 1 THEN BEGIN
      yray_avg = dblarr(n_avg)
    ENDIF ELSE BEGIN
      IF dim_numb EQ 2 THEN BEGIN
        yray_avg = dblarr(n_avg,yray_size(2))
      ENDIF ELSE BEGIN
        PRINT, 'Data array has more than two dimensions'
      ENDELSE
    ENDELSE
    
    ;------------------------------------------------------------------
    ; Time average
    ;------------------------------------------------------------------
    jj = 0l
    FOR ii = 0l, n_avg - 1 DO BEGIN 
      
      av_ind = where((time - time(0)) GE (ii*avg) AND $
                     (time - time(0)) LT ((ii+1)*avg), av_cnt)
      IF av_cnt GT 0 THEN BEGIN ; if there are measurements within 
                                ; an AVG sec interval
        
          time_avg(jj)   = TOTAL(time(av_ind))   / av_cnt
          IF dim_numb EQ 1 THEN BEGIN
              yray_avg(jj, *) = MEAN(yray(av_ind),  /NaN)
          ENDIF ELSE BEGIN
              IF KEYWORD_SET(sumup) THEN BEGIN 
                  yray_avg(jj, *) = TOTAL(yray(av_ind, *), 1)
              ENDIF ELSE BEGIN 
                  yray_avg(jj, *) = TOTAL(yray(av_ind, *), 1) /av_cnt
              ENDELSE 
          ENDELSE       
          jj = jj + 1
      ENDIF
      
    ENDFOR

;   IF jj LT 3 THEN GOTO, next  ; disabled by Jing
;   Adding for the condition jj eq 1 or 2
    IF jj GE 2 THEN BEGIN 
        time = time_avg(0:jj-2)
        IF dim_numb EQ 1 THEN BEGIN
            yray = yray_avg(0:jj-2)
        ENDIF ELSE BEGIN
            yray = yray_avg(0:jj-2, *)
            IF vfound EQ 1 THEN vray = data.v(0:jj-2, *)
        ENDELSE
    ENDIF ELSE BEGIN 
        time = time_avg(0)
        IF dim_numb EQ 1 THEN BEGIN
            yray = yray_avg(0)
        ENDIF ELSE BEGIN
            yray = yray_avg(0, *)
            IF vfound EQ 1 THEN vray = data.v(0, *)
        ENDELSE
    ENDELSE 
    ;------------------------------------------------------------------
    ; If keyword NEW_NAME is set the time averaged variable is stored
    ; with a new name that includes the averaging time
    ;------------------------------------------------------------------
    IF KEYWORD_SET(NEW_NAME) THEN BEGIN
      new_var_name = var_names(iv) + '_AVG' + $
        strcompress(string(avg, format='(i5)'), /REMOVE_ALL)
    ENDIF ELSE BEGIN
      new_var_name = var_names(iv)
    ENDELSE
    
    IF dim_numb EQ 1 THEN BEGIN
      datastr = {x:time, y:yray}
    ENDIF ELSE BEGIN
      IF vfound EQ 1 THEN BEGIN
        datastr = {x:time, y:yray, v:vray}
      ENDIF ELSE BEGIN
        datastr = {x:time, y:yray}
      ENDELSE
    ENDELSE
    
    store_data, new_var_name, data=datastr, dlim=dlim, lim=lim
    
    ;------------------------------------------------------------------
    ; Set plot attributes for new variable names
    ;------------------------------------------------------------------
    IF KEYWORD_SET(NEW_NAME) THEN BEGIN
      
      IF strpos(var_names(iv), '_UN') NE -1 THEN BEGIN
        dumstr = strmid(var_names(iv), $
                        strpos(var_names(iv), '_UN')+3,$
                        strlen(var_names(iv))-1)
        units_name = strmid(dumstr, 0, strpos(dumstr, '_'))
        CASE STRUPCASE(units_name) OF
          'COUNTS': uname = 'COUNTS'
          'BCOUNTS': uname = '1/tof_bin'
          'NCOUNTS': uname = '1/bin'
          'RATE': uname = '1/s'
          'NRATE': uname = '1/s-bin'
          'EFLUX': uname = 'eV/cm!E2!N-s-sr-eV'
          'DIFFFLUX': uname = '1/cm!E2!N-s-sr-(eV/e)'
          'DISTFUNC': uname = 's!E3!N/cm!E3!N-km!E3!N'
        ENDCASE
        
      ENDIF
      
      IF STRMID(var_names(iv),0,6) EQ 'ENSPEC' THEN BEGIN
        options, new_var_name, 'spec',1
        options, new_var_name, 'x_no_interp',1
        options, new_var_name, 'y_no_interp',1
        options, new_var_name, 'ytitle', $
          'SC' + string(sat, format='(i1.1)') + '!C!C' + $
          specie_str(specie)+' (eV)'
        options, new_var_name, 'ztitle', uname
        ylim,    new_var_name,  20., 4.5e4, 1
        
;        zlim, '*SP0*', 1.e0,  1.e3, 1
;        zlim, '*SP1*', 1.e-1, 1.e2, 1
;        zlim, '*SP2*', 1.e-2, 1.e1, 1
;        zlim, '*SP3*', 1.e-1, 1.e1, 1
      ENDIF
      
      IF STRMID(var_names(iv),0,6) EQ 'PASPEC' THEN BEGIN
        options, new_var_name, 'spec',1
        options, new_var_name, 'x_no_interp',1
        options, new_var_name, 'y_no_interp',1
        options, new_var_name, 'ytitle', $
          'SC' + string(sat, format='(i1.1)') + '!C!C' + $
          specie_str(specie)+' (deg)'
        options, new_var_name, 'ztitle', uname
        ylim,    new_var_name,  0., 360., 0
        
;        zlim, '*SP0*', 1.e0,  1.e3, 1
;        zlim, '*SP1*', 1.e-1, 1.e2, 1
;        zlim, '*SP2*', 1.e-2, 1.e1, 1
;        zlim, '*SP3*', 1.e-1, 1.e1, 1
      ENDIF
      
      IF STRMID(var_names(iv),0,3) EQ 'MAG' THEN BEGIN
        ylim,    new_var_name,  0., 0., 0
      ENDIF
      
    ENDIF
    
    next:
  ENDFOR
END
  
