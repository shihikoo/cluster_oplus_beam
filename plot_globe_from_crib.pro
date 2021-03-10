;+
; PROCEDURE: plot_globe_from_crib
;
; PURPOSE:
;         to interface the globe routine with the crib sheets
;
; Created by: C. Mouikis
;
; LAST MODIFIED: 03/06/02
; 
; MODIFICATION HISTORY:
;        03/06/02 -> Keyword OLD_EFF introduced
;        06/20/02 -> get_theta_phi routine is called with a new input
;                    parameter: inst
;        09/12/03 -> Keyword INCRATES introduced
;
;-
PRO plot_globe_from_crib, sat, specie, inst, units_name, eff_table, $
                          INCRATES=INCRATES, $
                          FORNACON=FORNACON, $
                          IC=IC, $
                          CNES=CNES, $
                          MAG_THETA = MAG_THETA, $
                          MAG_PHI = MAG_PHI, $
                          OLD_EFF = OLD_EFF, $
                          BKG=BKG

  COMMON get_error, get_err_no, get_err_msg, default_verbose
  
  ;--------------------------------------------------------------------
  ; Products that correspond to a specie
  ;--------------------------------------------------------------------
  CASE specie OF
    0: prod = [12, 13]
    1: prod = [15, 16]
    2: prod = [17, 18, 46, 48]
    3: prod = [17, 18, 47, 49]
    ELSE: BEGIN
      print, 'Specie numbers must be 0 for H+, 1 for He++,'
      print, '2 for He+ or 3 for O+'
      stop
    END
  ENDCASE

  ;--------------------------------------------------------------------
  ;Loop over all products that correspond to a specie
  ;--------------------------------------------------------------------
  jj = 0
  found = 0
  WHILE found EQ 0 AND jj LT n_elements(prod) DO BEGIN
    
    name = 'GLOBE_SC'+strcompress(sat,/remove_all)+'_'+$
      strcompress(units_name,/remove_all)+'_PR'+$
      strcompress(string(prod(jj)),/remove_all)+'_SP' +$
      strcompress(specie,/remove_all) ;name of TPLOT structure
    
    ;------------------------------------------------------------------
    ; dat structure is needed for the
    ; GSE -> CODIF coordinate transformation
    ;------------------------------------------------------------------
    dat = call_function('get_cis_cod_data', prod(jj), specie = specie, sat)

; jing: some strange errors happens
; when only zero data in product47 loaded, since they don't load 49 then
    IF get_err_no EQ 0 THEN BEGIN $
      IF total(where(dat.data)) EQ -1 THEN get_err_no = 1
    ENDIF 

    IF get_err_no EQ 0 THEN BEGIN ; check if data were found for time interval
      found = 1
      
      IF NOT(KEYWORD_SET(MAG_THETA)) OR NOT(KEYWORD_SET(MAG_PHI)) THEN BEGIN
        mag_theta = 0.
        mag_phi = 180.
        inst = 0
        get_theta_phi, sat, dat, mag_theta, mag_phi, inst, $
          IC=IC, CNES=CNES, $
          FORNACON=FORNACON     ;get magnetic field data
      ENDIF

      codif_globe, sat, inst, prod(jj), eff_table, $
        INCRATES=INCRATES, $
        SPECIE=SPECIE, $
        UNITS=units_name, $
        NAME = name, $
        MAG_THETA = mag_theta, $
        MAG_PHI = mag_phi, $
        OLD_EFF = OLD_EFF, $
        BKG=BKG

    ENDIF
    jj = jj + 1
  ENDWHILE
  
END
