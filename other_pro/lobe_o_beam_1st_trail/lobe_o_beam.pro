PRO lobe_o_beam

  time = '2001-09-01/00:00:00'
  timespan, time, 1, /DAYS
  lobe_def = 1 * 60 * 60  ;seconds

  sc = 4
  sc_str = STRING(sc, FORMAT = '(i1.1)')

;---------------------------------------------------------------------
; Read storm minimum Dst list from file: storm_min_dst.dat
; (this list is extracted from storm list spreadsheet)
; Store these times in variable min_dst (in tplot time format)
;---------------------------------------------------------------------

  OPENR, unit, 'storm_min_dst.dat', /GET_LUN
  
  min_dst = DBLARR(300)         ; assume no more than 300 events
  dummy = ''
  jj = 0
  WHILE NOT EOF(unit) DO BEGIN
      
      READF, unit, dummy
      min_dst(jj) = time_double(dummy)
      jj = jj + 1
      
  ENDWHILE

  min_dst = min_dst(0:jj-1)
  
  CLOSE, unit

;---------------------------------------------------------------------
; Read CLUSTER perigee times list from file: sc4_perigee_times.dat
; (different list for each S/C)
; Store these times in variable petime (in tplot time format)
;---------------------------------------------------------------------
  OPENR, unit, 'sc' + sc_str + '_perigee_times.dat', /GET_LUN

  petime = DBLARR(3000)         ; assumes no more than 3000 perigee passes
  dummy = ''
  jj = 0l
  WHILE NOT EOF(unit) DO BEGIN
    
      READF, unit, dummy
      petime(jj) = time_double(dummy)
      jj = jj + 1
    
  ENDWHILE

  petime = petime(0:jj-1)

  CLOSE, unit

;---------------------------------------------------------------------
; For each minimum Dst time choose the time intarval of the plots.
; The time interval is selected is of three orbits duration
;---------------------------------------------------------------------
  FOR ii = 0, N_ELEMENTS(min_dst)-1 DO BEGIN

    ;-----------------------------------------------------------------
    ; find the time of the perigee pass closest to the minimum Dst
    ;-----------------------------------------------------------------
      tt = SORT(ABS(petime - min_dst(ii)))
      closest_per = petime(tt(0))

    ;-----------------------------------------------------------------
    ; time interval start time, end time and dt
    ;-----------------------------------------------------------------
      time_start = petime(tt(0)-1)
      time_end = petime(tt(0)+2)
      dt = time_end - time_start
      time = time_start
      timespan, time_start, dt, /SECONDS

     ;-----------------------------------------------------------------
     ; Load Dst
     ;-----------------------------------------------------------------
;      read_omni
    
     ;------------------------------------------------------------
     ; Load energy spectra - tailward
     ;---------------------------------------------
      sat = [sc, sc]
      specie = [0, 3]
      angle = [[-90, 90], [90, 270]]
      inst = 0 & units_name = 'DIFF FLUX' & eff_table = 0
      plot_en_spec_from_crib, sat, specie, inst, $
        units_name, angle, eff_table, recalc = 0
  
     ;------------------------------------------------------------
     ; Load energy spectra - earthward
     ;----------------------------------------------------------
      sat = [sc, sc]
      specie = [0, 3]
      angle = [[-90, 90], [270, 90]]
      inst = 0 & units_name = 'DIFF FLUX' & eff_table = 0
      plot_en_spec_from_crib, sat, specie, inst, $
        units_name, angle, eff_table, recalc = 0
   
      ;-----------------------------------------------------------------
      ; Load CLUSTER Magnetic field
      ;----------------------------------------------------------------
      sat = sc
      plot_mag_from_crib, $
        sat, $
        POLAR = 1, $
        GSM = 1
      
      ;-----------------------------------------------------------------
      ; Load CLUSTER moments (n, P)
      ;----------------------------------------------------------------   
      sat = [sc, sc]
      specie = [0, 3]
      moments = ['D', 'D']
      angle = [[-90.0, 90.0], [0., 360.]]
      energy = [40.0, 40000.0]
      inst = 0 & eff_table = 0
      plot_3dmom_from_crib, sat, specie, inst, $
        moments, angle, energy, eff_table
    
      ;-------------------------------------------------------------------
      ; Calculate total pressure & beta
      ;--------------------------------------------------------------
      h1_press = 'TDMOM_EN00040_40000_SC' + $
                 STRCOMPRESS(sc, /REMOVE_ALL)$
                 + '_MTPRESSURE_SP0_ET0_All'
      o1_press = 'TDMOM_EN00040_40000_SC' + $
                 STRCOMPRESS(sc, /REMOVE_ALL)$
                 + '_MTPRESSURE_SP3_ET0_All'
      mag_press = 'MAG_SC'+ $
                  STRCOMPRESS(sc, /REMOVE_ALL)$
                  + '_B_xyz_gse_MAG_PR'
      plasma_beta, h1_press, mag_press, O1_PRESSURE = o1_press
    
      ;-----------------------------------------------------------------
      ; Load CLUSTER ephemeris
      ;-----------------------------------------------------------------    
      sat = [sc]
      get_cluster_ephemeris, sat, /GSE_X, /GSE_Y, /GSE_Z, /DIST, /MLT

      
       ;-----------------------------------------------------------------
       ;select the lobe range
       ;---------------------------------------------------------------
       
       ;Read plasama_beta data
        p_beta = 'TDMOM_EN00040_40000_SC' + $
                  STRCOMPRESS(sc, /REMOVE_ALL)$
                  + '_MTPRESSURE_SP0_ET0_All_O1_beta'
     
        IF NOT(size(p_beta, /TYPE)) EQ 7 THEN BEGIN ; var -> integer
            tplot_names, p_beta, NAMES = var_name
            p_beta = var_name(0)
        ENDIF
        
        get_data, p_beta, data = pb
        
        time_pb = pb.x
        data_pb = pb.y
        

        ;find lobes
      
        lobe_start = DBLARR(3000)
        lobe_end = DBLARR(3000)
        lobe_dt = lobe_start - lobe_end
        kk_lobe = 0
        kk_pb = 0L
        WHILE kk_pb LT N_ELEMENTS(data_pb) - 1 DO BEGIN
            IF data_pb(kk_pb) GT  0.1 THEN BEGIN
                kk_pb = kk_pb + 1
            ENDIF ELSE BEGIN 
                lobe_start (kk_lobe) = time_pb(kk_pb)
                WHILE data_pb(kk_pb)LE 0.1 AND kk_pb LT N_ELEMENTS(data_pb)-1 DO BEGIN
                    kk_pb = kk_pb + 1
                ENDWHILE
                IF data_pb(kk_pb) GE 0.1 THEN BEGIN
                    lobe_end(kk_lobe) = time_pb(kk_pb)
                    lobe_dt(kk_lobe) = lobe_end(kk_lobe)-lobe_start(kk_lobe)
       
                  ;chose lobe time longer than lobe_def        
                    IF lobe_dt(kk_lobe) GE lobe_def THEN BEGIN  
                        kk_lobe = kk_lobe + 1
                    ENDIF
                ENDIF
                kk_pb = kk_pb + 1
            ENDELSE
        ENDWHILE
      
        lobe_start = lobe_start(0:kk_lobe-1)
        lobe_end = lobe_end(0:kk_lobe-1)
        lobe_dt = lobe_dt(0:kk_lobe-1)
  
       ;----------------------------------------------------------------
       ;judge whether there is a O+ beam or not by the velocity
       ;------------------------------------------------------------ ---

       ;Read velocity data 
        o_vx = 'TDMOM_EN00040_40000_SC' $
               + STRCOMPRESS(sc, /REMOVE_ALL)$
               + '_MTVELOCITY_SP3_ET0_All_X'

        IF NOT(size(o_vx, /TYPE)) EQ 7 THEN BEGIN ; var -> integer
            tplot_names, o_vx, NAMES = var_name
            o_vx = var_name(0)
        ENDIF
        
        get_data, o_vx, data = ovx
        
        time_ovx = ovx.x
        data_ovx = ovx.y

        ; chose the data sets which has more than required number (beam_density_def) of 
        ; fast data with tailward (negative) velocity smaller than beam_v_def
        
        beam_density_def = 5
        beam_v_def = -200 ;km/s
        kk_o = 0
        lobe_o_beam_start = DBLARR(300)
        lobe_o_beam_end = DBLARR(300)
        lobe_o_beam_dt = DBLARR(300)

        FOR jjj = 0, N_ELEMENTS(lobe_dt)-1 DO BEGIN
            lloc = where(time_ovx GE lobe_start(jjj) AND time_ovx LE lobe_end(jjj))
            IF total (lloc) NE -1 THEN begin
                o_n = total( data_ovx(lloc) LE beam_v_def) 
                IF o_n GE beam_density_def THEN BEGIN
                    lobe_o_beam_start(kk_o) = lobe_start(jjj)
                    lobe_o_beam_end(kk_o) = lobe_end(jjj)
                    lobe_o_beam_dt(kk_o) = lobe_o_beam_end(kk_o) - lobe_o_beam_start(kk_o)
                    kk_o = kk_o + 1
                ENDIF
            ENDIF
        ENDFOR

        lobe_o_beam_start = lobe_o_beam_start(0:kk_o-1)
        lobe_o_beam_end = lobe_o_beam_end(0:kk_o-1)
        lobe_o_beam_dt = lobe_o_beam_dt(0:kk_o-1)

       ;---------------------------------------------------------------
       ;Overview plots
       ;--------------------------------------------------------------
   


        FOR jjj = 0, N_ElEMENTS(lobe_o_beam_dt)-1 DO BEGIN
            time = lobe_o_beam_start(jjj)
            timespan, time, lobe_o_beam_dt(jjj), /SECONDS
            
            p01 = 'TDMOM_EN00040_40000_SC' $
                  + STRCOMPRESS(sc, /REMOVE_ALL)$
                  + '_MTPRESSURE_SP0_ET0_All_O1_P_total'
            
            p02 = 'TDMOM_EN00040_40000_SC' $ 
                  + STRCOMPRESS(sc, /REMOVE_ALL)$ 
                  + '_MTPRESSURE_SP0_ET0_All_O1_beta'
  
   ;    p03 = 'ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP0_ET0_All'
      
            p04 = 'ENSPEC_SC' $ 
                  + STRCOMPRESS(sc,  /REMOVE_ALL)$
                  +'_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All'

   ;    p05 = 'ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP0_ET0_All'
       
            p06 = 'ENSPEC_SC' $ 
                  + STRCOMPRESS(sc,  /REMOVE_ALL)$
                  + '_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All'
      
            p07 = 'TDMOM_EN00040_40000_SC' $
                  + STRCOMPRESS(sc, /REMOVE_ALL) $
                  + '_MTDENSITY_SP0_ET0_All'
        
            p08 = 'TDMOM_EN00040_40000_SC' $
                  +  STRCOMPRESS(sc, /REMOVE_ALL)$
                  + '_MTPRESSURE_SP0_ET0_All_T' 
              
  ;     p09 = 'Dst_Index'
              
            p10 = 'TDMOM_EN00040_40000_SC' $
                  + STRCOMPRESS(sc, /REMOVE_ALL)$
                  + '_MTDENSITY_SP3_ET0_All'
            
            p11 = 'TDMOM_EN00040_40000_SC' $
                  + STRCOMPRESS(sc, /REMOVE_ALL)$
                  + '_MTPRESSURE_SP3_ET0_All_T' 

            p12 = 'TDMOM_EN00040_40000_SC' $
                  + STRCOMPRESS(sc, /REMOVE_ALL)$
                  + '_MTVELOCITY_SP3_ET0_All_X' 
        
            options, '*', 'panel_size', 1
            options, [p10, p11], 'color', 2
        
            options, [p06], 'ztitle', ''
      
            ylim, p01, 0.01, 3, 1
            ylim, p02, 0.001, 1
            ylim, p07, 0.01, 10
            ylim, p08, 0.001, 10

   ;    options, p03, 'ytitle', 'SC' + STRCOMPRESS(sc, /REMOVE_ALL) + $
   ;              '!C!CH!U+!N (eV)' + '!C!C' + 'Tailward'
            options, p04, 'ytitle', 'SC' + STRCOMPRESS(sc, /REMOVE_ALL) + $
                     '!C!CO!U+!N (eV)' + '!C!C' + 'Tailward'
   ;    options, p05, 'ytitle', 'SC' + STRCOMPRESS(sc, /REMOVE_ALL) + $
   ;              '!C!CH!U+!N (eV)' + '!C!C' + 'Earthward'
            options, p06, 'ytitle', 'SC' + STRCOMPRESS(sc, /REMOVE_ALL) + $
                     '!C!CO!U+!N (eV)' + '!C!C' + 'Earthward'

            var_label = 'EPH_SC' + sc_str + '_'
            var_label = var_label + ['MLT', 'GSE_X', 'GSE_Y', 'GSE_Z', 'DIST']
 
       ; IF KEYWORD_SET(PS) THEN BEGIN
        ;    output_name = 'storm_page1_' + STRMID(time, 0, 10)
         ;   popen, output_name
       ; ENDIF ELSE WINDOW, /FREE
        
     
            tplot, [p01, p02, p04, p06, p07, p08, p12], var_label = var_label
            
            tplot_panel, v = p07, o = p10
            tplot_panel, v = p08, o = p11
        
            yline, p02, offset = 0.1
            yline, p12, offset = -200

     ;        yline, p09
   
     
      
      ;IF KEYWORD_SET(PS) THEN pclose
     

      ;-----------------------------------------------------------------
      ; Plot in PS file
      ;-----------------------------------------------------------------
        dds = time_string(lobe_o_beam_start(jjj))
        dde = time_string(lobe_o_beam_end(jjj))
        fln = 'plots/storm_lobe_o_beam' + $
              STRMID(dds, 0, 4) + STRMID(dds, 5, 2) + STRMID(dds, 8, 2) + '_' $
              + STRMID(dds, 11, 2) + STRMID(dds, 14, 2) + STRMID(dds, 17, 2) + 'to' $
              + STRMID(dde, 0, 4) + STRMID(dde, 5, 2) + STRMID(dde, 8, 2) + '_' $
              + STRMID(dde, 11, 2) + STRMID(dde, 14, 2) + STRMID(dde, 17, 2)  $
              + '.ps'
        popen, fln
        
        tplot, [p01, p02, p04, p06, p07, p08, p12], var_label = var_labe
        tplot_panel, v = p07, o = p10
        tplot_panel, v = p08, o = p11
        
        yline, p02, offset = 0.1
        yline, p12, offset = -200
 ;       yline, p09
        
        pclose
    ENDFOR
    
ENDFOR

END
