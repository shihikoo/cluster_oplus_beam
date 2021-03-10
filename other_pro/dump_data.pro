;+
; PROCEDURE  dump_data
;
; PURPOSE: to dump in a file data stored as tplot variables. The
;          first column is the time in yyyy-mm-dd/hh:mm:ss format
;
; INPUT:
;        tplot_name: name or number of tplot variable as it appears 
;                    when the tplot_names command is used. Only one 
;                    variable is allowed.
;
;
; OUTPUT:
;         data are saved in the specified file
;
; KEYWORDS: FILE_OUT -> specify filename to dump the data
;           APPEND   -> Append variable in an already existing dump
;                       file. This allows the possibility to dump more
;                       than one tplot variables
;           SECONDS  -> Include a column with seconds from the
;                       begining of day (does not work with the APPEND keyword)
;           NODATA   -> Create columns with time information only
;           INTERP   -> Interpolate data on times that it reads in
;                       file FILE_OUT (works only with APPEND kweyword)
;                       At the moment it only handles 1D data and the
;                       purpose was to allow the possibility of appending
;                       ephemeris data in a dump file
;
; CREATED BY: C. Mouikis
;
; LAST MODIFICATION: 11/07/2001
;
; MODIFICATION HISTORY:
;     07/20/01: the tplot variable name is written at the top of the
;               file
;     09/17/01: dump list can be appended using the APPEND
;               keyword. This assumes that the number of time steps is
;               exactly the same
;     10/29/01: In case of 2D data the V information (i.e energy) is included
;     10/29/01: Keyword SECONDS introduced.
;     11/07/01: Keyword NODATA introduced.
;     11/07/01: Keyword INTERP introduced. At the moment it only
;               handles 1D data and the purpose was to allow the
;               possibility of appending ephemeris data in a dump file
;-
PRO dump_data, tplot_name, $
               FILE_OUT = FILE_OUT, $
               APPEND = APPEND, $
               SECONDS = SECONDS, $
               NODATA = NODATA, $
               INTERP = INTERP, $
               TITLE = TITLE, $
               HELP=HELP

  ;--------------------------------------------------------------------
  ; Help text
  ;--------------------------------------------------------------------
  IF KEYWORD_SET(HELP) THEN BEGIN
    
    PRINT, '- tplot_name'
    PRINT, '- FILE_OUT = FILE_OUT'
    PRINT, '- APPEND = APPEND'
    PRINT, '- SECONDS = SECONDS'
    PRINT, '- NODATA = NODATA'
    PRINT, '- INTERP = INTERP'
    PRINT, '- TITLE = TITLE'
    PRINT, '- HELP = HELP'

    RETURN
  ENDIF
  
  ;--------------------------------------------------------------------
  ; Find the name of the tplot variable to be dumped
  ;--------------------------------------------------------------------
  IF size5(tplot_name, /TYPE) EQ 7 THEN BEGIN ; var -> string
    IF tplot_name EQ '*' THEN BEGIN
      tplot_names, NAMES=var_names
    ENDIF ELSE BEGIN
      var_names=tplot_name
    ENDELSE
  ENDIF ELSE BEGIN                           ; var -> integer
    tplot_names, tplot_name, NAMES=var_names
  ENDELSE
  
  IF N_ELEMENTS(var_names) NE 1 THEN BEGIN
    MESSAGE, 'Only one variable is allowed. To output more variables the APPEND keyword should be used', /INF
    RETURN
  ENDIF
  
  ;--------------------------------------------------------------------
  ; Set output filename
  ;--------------------------------------------------------------------
  IF NOT KEYWORD_SET(FILE_OUT) THEN file_out = 'dump.dat'
  
  ;--------------------------------------------------------------------
  ; If keyword APPEND is set, check if output file exists. If file_out
  ; does not exist then disregard keyword APPEND
  ;--------------------------------------------------------------------
  IF KEYWORD_SET(APPEND) THEN BEGIN
    find = FINDFILE(file_out, COUNT=ff)
    IF ff EQ 0 THEN BEGIN
      append = 0
    ENDIF
  ENDIF
  
  ;--------------------------------------------------------------------
  ; Get the data
  ;--------------------------------------------------------------------
  get_data, var_names(0), data=data, dlim=dlim
  time = time_string_unh(data.x, /SQL)
  IF KEYWORD_SET(SECONDS) AND NOT KEYWORD_SET(APPEND) THEN BEGIN
    secs = time_tplot_secs_of_day(data.x)
  ENDIF
  
  ;--------------------------------------------------------------------
  ; Adjust TITLE LENGTH. Set dashes if no title is set.
  ;--------------------------------------------------------------------
  IF KEYWORD_SET(TITLE) THEN BEGIN
    title = STRING(title, FORMAT='(a15)')
  ENDIF ELSE BEGIN
    title = '---------------'
  ENDELSE
  
  ;--------------------------------------------------------------------
  ; Write the data into file
  ;--------------------------------------------------------------------
  IF NOT KEYWORD_SET(APPEND) THEN BEGIN ; If new file
    OPENW, unitw, file_out, /GET_LUN, /APPEND
    array_dim = SIZE(data.y, /N_DIMENSIONS)
    
    ;------------------------------------------------------------------
    ; Write tplot variable name as header. Nothing if NODATA keyword
    ; is set.
    ;------------------------------------------------------------------
  ;  IF NOT KEYWORD_SET(NODATA) THEN BEGIN
  ;    lstr = var_names(0)
  ;    PRINTF, unitw, lstr
  ;  ENDIF ELSE BEGIN
  ;    PRINTF, unitw, ''
  ;  ENDELSE
    
    ;------------------------------------------------------------------
    ; Write column titles. If variable is a spectrum the values of
    ; variable V (of tplot structure) are used as column titles
    ;------------------------------------------------------------------
    lstr = 'Date       ' + $
      STRING('Time', FORMAT='(a12)')
    IF KEYWORD_SET(SECONDS) THEN BEGIN
      lstr = lstr + STRING('Seconds', FORMAT='(a14)')
    ENDIF
    
    IF NOT KEYWORD_SET(NODATA) THEN BEGIN
      IF array_dim EQ 1 THEN BEGIN ; scalar data just write title
        lstr = lstr + '  ' + title
      ENDIF ELSE BEGIN
        IF array_dim EQ 2 THEN BEGIN
          str_element, data, 'v', v, success=vs
          IF vs EQ 1 THEN BEGIN ; Spectra write V value
            FOR jj = 0, N_ELEMENTS(data.v(0,*))-1 DO BEGIN
                lstr = lstr + '  '+ data.v(0, jj) ;Jing: since I use data.v as
            ENDFOR                                ;title string. I remove the 
                                                  ;coversion to string
          ENDIF ELSE BEGIN ; Vector data write title in each col.
            FOR jj = 0, N_ELEMENTS(data.y(0,*))-1 DO BEGIN
              lstr = lstr + '  ' + title
            ENDFOR
          ENDELSE
        ENDIF ELSE BEGIN ; More than 2D are not handled
          MESSAGE, 'Dimensions more than 2 are not supported', /INF
          STOP
        ENDELSE
      ENDELSE
    ENDIF
    
    PRINTF, unitw, lstr
    
    ;------------------------------------------------------------------
    ; Write data. Format f15.6 is assumed for all data.
    ;------------------------------------------------------------------
    FOR ii=0l, N_ELEMENTS(time)-1 DO BEGIN
      
      IF KEYWORD_SET(SECONDS) THEN BEGIN
        lstr = time(ii) + '  ' + STRING(secs(ii), FORMAT='(f12.3)')
      ENDIF ELSE BEGIN
        lstr = time(ii)
      ENDELSE
      
      IF NOT KEYWORD_SET(NODATA) THEN BEGIN
        IF array_dim EQ 1 THEN BEGIN ; 1D data arrays
          lstr = lstr + '  ' + STRING(data.y(ii), FORMAT='(f15.6)')
          PRINTF, unitw, lstr
        ENDIF ELSE BEGIN
          IF array_dim EQ 2 THEN BEGIN ; 2D data arrays
            FOR jj = 0, N_ELEMENTS(data.y(0,*))-1 DO BEGIN
              lstr = lstr + '  ' + STRING(data.y(ii,jj), FORMAT='(f15.6)')
            ENDFOR
            PRINTF, unitw, lstr
          ENDIF ELSE RETURN     ; Greater than 2D data arrays are not supported
        ENDELSE
      ENDIF ELSE BEGIN
        PRINTF, unitw, lstr
      ENDELSE
    ENDFOR
    FREE_LUN, unitw, /FORCE
        
  ENDIF ELSE BEGIN ; Append data in existing file
    IF KEYWORD_SET(INTERP) THEN BEGIN ; Interpolate data
      
      ;----------------------------------------------------------------
      ; open file to be appended/read #lines
      ;----------------------------------------------------------------
      lstr = ''
      OPENR, unitr, file_out, /GET_LUN
      
      l1 = 0
      WHILE NOT(EOF(unitr)) DO BEGIN
        READF, unitr, lstr
        l1 = l1 + 1
      ENDWHILE
      l1 = l1-2 
      
      FREE_LUN, unitr, /FORCE
      ;----------------------------------------------------------------
      ; open file to be appended/read times
      ;----------------------------------------------------------------
      time_str = STRARR(l1)
      OPENR, unitr, file_out, /GET_LUN
      
      READF, unitr, lstr        ; read header
      READF, unitr, lstr        ; read header
      FOR ii = 0, l1-1 DO BEGIN
        READF, unitr, lstr
        time_str(ii) = STRMID(lstr, 0, 23)
      ENDFOR
      time_old = time_double(time_str)
      FREE_LUN, unitr, /FORCE
      ;----------------------------------------------------------------
      ; Interpolate data
      ;----------------------------------------------------------------
      data_new = INTERPOL(data.y, data.x, time_old)
      
      ;----------------------------------------------------------------
      ; Append file with interpolated data
      ;----------------------------------------------------------------
      OPENR, unitr, file_out, /GET_LUN ; open file to be appended
      OPENW, unitw, file_out+'_tmp', /GET_LUN ; open tmp file
    
      READF, unitr, lstr        ; read header
      PRINTF, unitw, lstr + '/' + var_names(0) ; write header with new variable
      
      array_dim = SIZE(data.y, /N_DIMENSIONS)
      
      READF, unitr, lstr        ; read header
      
      IF array_dim EQ 1 THEN BEGIN
        lstr = lstr + '  ' + title
      ENDIF ELSE BEGIN
        IF array_dim EQ 2 THEN BEGIN
          str_element, data, 'v', v, success=vs
          IF vs EQ 1 THEN BEGIN
            FOR jj = 0, N_ELEMENTS(data.v(0,*))-1 DO BEGIN
              lstr = lstr + ' '+ data.v(0,jj) ; Jing: also here
            ENDFOR
          ENDIF ELSE BEGIN
            FOR jj = 0, N_ELEMENTS(data.y(0,*))-1 DO BEGIN
              lstr = lstr + '  ' + title
            ENDFOR
          ENDELSE
        ENDIF ELSE BEGIN
          MESSAGE, 'Dimensions more than 2 are not supported', /INF
          STOP
        ENDELSE
      ENDELSE
      
      PRINTF, unitw, lstr
      
      FOR ii=0, N_ELEMENTS(time_old)-1 DO BEGIN
        READF, unitr, lstr
        
        IF array_dim EQ 1 THEN BEGIN
          lstr = lstr + '  ' + STRING(data_new(ii), FORMAT='(f15.6)')
          PRINTF, unitw, lstr
        ENDIF ELSE BEGIN
          IF array_dim EQ 2 THEN BEGIN ; 2D data arrays
            FOR jj = 0, N_ELEMENTS(data.y(0,*))-1 DO BEGIN
              lstr = lstr + '  ' + STRING(data.y(ii,jj), FORMAT='(f15.6)')
            ENDFOR
            PRINTF, unitw, lstr
          ENDIF ELSE RETURN     ; Greater than 2D data arrays are not supported
        ENDELSE
      ENDFOR      
      
      
      FREE_LUN, unitr, /FORCE
      FREE_LUN, unitw, /FORCE
      
      SPAWN, 'mv ' + file_out+'_tmp ' + file_out
      
    ENDIF ELSE BEGIN ; Append data to a file (assumes same # of time steps)
      OPENR, unitr, file_out, /GET_LUN ; open file to be appended
      OPENW, unitw, file_out+'_tmp', /GET_LUN ; open tmp file
    
      lstr = ''
      READF, unitr, lstr        ; read header
      PRINTF, unitw, lstr + '/' + var_names(0) ; write header with new variable
      
      array_dim = SIZE(data.y, /N_DIMENSIONS)
      
      READF, unitr, lstr        ; read header
      
      IF array_dim EQ 1 THEN BEGIN
        lstr = lstr + '  ' + title
      ENDIF ELSE BEGIN
        IF array_dim EQ 2 THEN BEGIN
          str_element, data, 'v', v, success=vs
          IF vs EQ 1 THEN BEGIN
            FOR jj = 0, N_ELEMENTS(data.v(0,*))-1 DO BEGIN
                lstr = lstr + '  ' + data.v(0, jj) ; Jing: I also removed the string
            ENDFOR                                     ; string conversion here
          ENDIF ELSE BEGIN
            FOR jj = 0, N_ELEMENTS(data.y(0,*))-1 DO BEGIN
              lstr = lstr + '  ' + title
            ENDFOR
          ENDELSE
        ENDIF ELSE BEGIN
          MESSAGE, 'Dimensions more than 2 are not supported', /INF
          STOP
        ENDELSE
      ENDELSE
      
      PRINTF, unitw, lstr
      
      FOR ii=0l, N_ELEMENTS(time)-1 DO BEGIN
        READF, unitr, lstr
        
        IF array_dim EQ 1 THEN BEGIN
          lstr = lstr + '  ' + STRING(data.y(ii), FORMAT='(f15.6)')
          PRINTF, unitw, lstr
        ENDIF ELSE BEGIN
          IF array_dim EQ 2 THEN BEGIN ; 2D data arrays
            FOR jj = 0, N_ELEMENTS(data.y(0,*))-1 DO BEGIN
              lstr = lstr + '  ' + STRING(data.y(ii,jj), FORMAT='(f15.6)')
            ENDFOR
            PRINTF, unitw, lstr
          ENDIF ELSE RETURN     ; Greater than 2D data arrays are not supported
        ENDELSE
      ENDFOR
      FREE_LUN, unitr, /FORCE
      FREE_LUN, unitw, /FORCE
      SPAWN, 'mv -f ' + file_out+'_tmp ' + file_out
    ENDELSE
  ENDELSE
   
END

