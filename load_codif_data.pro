;+
;PROCEDURE: plot_pa_spec_from_crib
;PURPOSE: 
;  It is a crib sheet for plotting cluster cis energy spectra using
;  the TPLOT package
;
;INPUT
;PARAMETERS:   sat:        Satellite number
;              inst:       Instrument (0: CODIF, 1: HIA)
;              time:       Date/time in tplot format
;              timespan:   Time span
;              units_name: Units for energy spectra
;                          'Counts', 'NCOUNTS', 'RATE', 'NRATE', 
;                           'DIFF FLUX', 'EFLUX'
;              angle:      Angle range to sum over [thetarange, phirange]
;
;                          THETA_88: [-78.75, -56.25, -33.75, -11.25,
;                                     11.25, 33.75, 56.25, 78.75] 
;
;                          PHI_88:   [11.25, 22.50, 33.75, 45.00, 56.25,
;                                     67.50, 78.75, 101.25, 112.50, 123.75,
;                                     135.00, 146.25, 157.50, 168.75, 191.25,
;                                     202.50, 213.75, 225.00, 236.25, 247.50,
;                                     258.75, 281.25, 292.50, 303.75, 315.00,
;                                     326.25, 337.50, 348.75]
;
;                          THETA_24: [-78.75, -56.25, -33.75, -11.25,
;                                     11.25, 33.75, 56.25, 78.75] 
;
;                          PHI_24:   [45.0, 90.0, 135.0, 225.0, 270.0, 315.0]

PRO plot_pa_spec_from_crib, sat, specie, inst, units_name, $
                            energy, eff_table
  
  COMMON get_error, get_err_no, get_err_msg, default_verbose
  
  specie_str = ['H!U+!N','He!U++!N','He!U+!N','O!U+!N']
  
  ; Check if sat and specie arrays have the same number of elements
  IF n_elements(sat) NE n_elements(specie) THEN BEGIN
    print, 'The sat array and the specie array MUST have'
    print, 'the same number of elements.'
    stop
  ENDIF
  
  ; Loop over all sat/specie combinations
  FOR ii=0,n_elements(specie)-1 DO BEGIN
    
    get_err_no = 0 ; reset error indicator
    ; Products that correspond to a specie
    CASE specie(ii) OF
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
    
    ;Loop over all products that correspond to a specie
    FOR jj = 0, n_elements(prod)-1 DO BEGIN
      
      dat = call_function('get_cis_cod_data',prod(jj),specie=specie(ii), sat)
      
      IF get_err_no EQ 0 THEN BEGIN ;check if data were found for time interval
        inst = 0
        get_theta_phi, sat, dat, mag_theta, mag_phi, inst
      ENDIF
    
    ENDFOR
  ENDFOR

END
