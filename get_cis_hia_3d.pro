;+
;FUNCTION:      get_cis_hia_3d, prod_num, sat, time=time
;INPUT:
;       prod_num: Product number
;	sat: Satellite number (1-4)
;	time: If set, return time array of data values (index ignored)
;
;PURPOSE:   Retrieve CIS/HIA 3D data.
;
;CREATED BY:    Peter Schroeder - Modified by C. mouikis
;
;LAST MODIFICATION:  06/05/03
;  09/12/01 - The correction of phi angle as a function of energy is
;             included 
;  01/29/02 - Bug fix. In the phi correction routine and for 31 energy
;             products the code should loop over 32 bins and not 31
;             (although the 32nd is not used)
;  01/29/02 - Keyword NO_PHI_COR introduced.
;  06/05/03 - Debugging
;-
FUNCTION get_cis_hia_3d, prod_num, sat, specie=specie, time=time, $
                         frencheff=frencheff
  
  COMMON get_error, get_err_no, get_err_msg, default_verbose
  
  
  get_err_no = 0
  get_err_msg= 'Data found in time interval'
  
  ;--------------------------------------------------------------------
  ; Number of energies and angles for each 3D product
  ;--------------------------------------------------------------------
  CASE prod_num OF
    23: enang = [31,8,16]
  ENDCASE
  
  nenergy = enang[0]
  nangles = enang[1]*enang[2]
  ntheta  = enang[1]
  nphi    = enang[2]
  
  ;--------------------------------------------------------------------
  ; Get the data
  ;--------------------------------------------------------------------
  options = [long(sat),1l,long(prod_num),0l]

  IF keyword_set(time) THEN RETURN,get_cis_data_times(options)

  basic_data = get_cis_all_data(options)
  IF NOT KEYWORD_SET(basic_data) THEN BEGIN ; Test if the product exists
    get_err_no = 1
    get_err_msg = 'Product '+ $
      string(prod_num,format='(i2.2)') + $
      ' does not exist'
    PRINT, get_err_msg
    RETURN, 0
  ENDIF
  packets = N_ELEMENTS(basic_data.data)

  ;--------------------------------------------------------------------
  ; Prepare data arrays
  ;--------------------------------------------------------------------
  dtime = cis_tplot_time(basic_data.header.time_in_ms) ; time -> tplot-time
  data = float(basic_data.data.data)
  data = REFORM(data, nenergy, nangles, packets)
  
  ;--------------------------------------------------------------------
  ; Limit time interval
  ;--------------------------------------------------------------------
  get_timespan,time_interval
  
  t_s=gettime(time_interval(0)) ; start time in tplot-time
  t_e=gettime(time_interval(1))   ; end time in tplot-time
  
  tind=where(dtime GE t_s AND dtime LT t_e, tc) ; check if there are data
  IF tc LT 1 THEN BEGIN ; changed to include    ; (more than one packet)
    get_err_no = 1      ; one packet            ; for time interval requested
    get_err_msg = 'No data in time interval'
    print, get_err_msg
    RETURN, 0
  ENDIF
  
  ind_s=tind(0)       ; index that corresponds to start time
  ind_e=tind(tc-1)    ; index that corresponds to end time
  
  ;--------------------------------------------------------------------

  dtime = dtime(ind_s:ind_e) ; limit dtime to time interval
  data = data(*,*,ind_s:ind_e) ; limit data to time interval
  ndatapoints = n_elements(dtime)
  
  op_mode = basic_data.header(ind_s:ind_e).op_mode ; limit op_mode to tint
  sensitivity = basic_data.header(ind_s:ind_e).sensitivity ; limit sens to tint

  ;--------------------------------------------------------------------
  ; get theta,phi,dtheta,dphi for a nangles product
  ;--------------------------------------------------------------------
  get_cis_hia_angles,ntheta,nphi,theta,phi,dtheta,dphi

  theta_arr = fltarr(nangles)
  phi_arr = fltarr(nangles)
  dtheta_arr = fltarr(nangles)
  dphi_arr = fltarr(nangles)

  ian = 0
  FOR i = 0,nphi-1 DO BEGIN
    FOR j = 0,ntheta-1 DO BEGIN
      theta_arr[ian] = theta[j]
      phi_arr[ian] = phi[i]
      dtheta_arr[ian] = dtheta[j]
      dphi_arr[ian] = dphi[i]
      ian = ian + 1
    ENDFOR
  ENDFOR
  theta=theta_arr
  phi=phi_arr
  dtheta=dtheta_arr
  dphi=dphi_arr
  ;--------------------------------------------------------------------
  ; Replicate angle definitions for all energy bins creating 2D arrays  
  ;--------------------------------------------------------------------
  th     = REPLICATE(1, 1, nangles) * theta
  theta  = REBIN(th, nenergy, nangles)
  ph     = REPLICATE(1, 1, nangles) * phi
  phi    = REBIN(ph, nenergy, nangles)
  dth    = REPLICATE(1, 1, nangles) * dtheta
  dtheta = REBIN(dth, nenergy, nangles)
  dph    = REPLICATE(1, 1, nangles) * dphi
  dphi   = REBIN(dph, nenergy, nangles)

  energy_arr = fltarr(nenergy,nangles,ndatapoints)
  denergy_arr = energy_arr
  geom_factor = dblarr(ndatapoints)
  delta_t = dblarr(ndatapoints)
  tlm_rate = basic_data.header(ind_s:ind_e).tlm_rate
  pspin = basic_data.header(ind_s:ind_e).pspin
  k1 = basic_data.data(ind_s:ind_e).k1
  k2 = basic_data.data(ind_s:ind_e).k2

  ;--------------------------------------------------------------------
  ; For all records get the geometric factor, energies and number of 
  ; spins
  ;--------------------------------------------------------------------
  FOR j = 0,ndatapoints-1 DO BEGIN
    IF KEYWORD_SET(frencheff) THEN BEGIN
      geom_factor[j] = get_cis_hia_geom_factor(dtime[j], sat, sensitivity[j])
    ENDIF ELSE BEGIN
      geom_factor[j] = jim_geom_factor(sat, sensitivity[j])
    ENDELSE
    get_cis_hia_energies, $
      dtime[j], $
      sat, $
      nenergy, $
      op_mode[j], $
      basic_data.data(j).hvtbl, $
      energy, $
      denergy, $
      sensitivity[j]
    
    ;------------------------------------------------------------------
    ; Check compressed/uncompressed mode and product numbers
    ;------------------------------------------------------------------
    new_prod_num = prod_num
    CASE prod_num OF
      
      23: IF k1(j) NE 255 THEN new_prod_num = 23
      
    ENDCASE

    ;------------------------------------------------------------------------
    ; get number of spins accumulated from calibration files
    num_spins = get_cis_hia_acc_spin(dtime[j], sat, tlm_rate[j], $
                                     op_mode[j], $
                                     new_prod_num, $
                                     basic_data)
;-------------
;Jing: due to function get_cis_hia_acc_spin, when there is no spin to
;calc, it will automatically read from the basic_data, which will
;result in num_spins is not the data at the time point "j". Instead,
;it's the num_spins of all time points, which means num_spins become
;an array of same dimention as delta_t. In consideration of this
;case, I added the following one line.
    IF n_elements(num_spins) EQ n_elements(delta_t) THEN num_spins = num_spins(j)
;--------
    delta_t[j] = num_spins*pspin[j]/1d3 ; delta time
    FOR i = 0,nangles-1 DO BEGIN
      energy_arr[*,i,j] = energy
      denergy_arr[*,i,j] = denergy
    ENDFOR
  ENDFOR

  ; now corect geometry factor for number of anodes summed and number of energy levels summed???
  ; prod 8 and 24 are corrected onboard
  IF prod_num EQ 23 THEN geom_factor=geom_factor*(16./ntheta);*(62./nenergy)

  ;--------------------------------------------------------------------
  ; Get efficiency
  ;--------------------------------------------------------------------
  gf = cis_hia_efficiency(dtime, $
                          sat, $
                          energy_arr, $
                          theta, $
                          anode_128_map(), $
                          frencheff=frencheff)

;----------------------------------------------------------------------
; Construct the data structure
;----------------------------------------------------------------------
  retdata = {                                                                 $
              project_name:        'CLUSTER CIS/HIA'    ,                         $
              data_name:           'Product '+strcompress(prod_num,/rem),     $
              data_product:        prod_num,                                  $
              units_name:          'Counts',                                  $
              units_procedure:     'convert_hia_units',                       $
              valid:               1,                                         $
              time:                dtime,                                     $
              end_time:            dtime+delta_t,                             $
              delta_t:             delta_t,                                   $
              integ_t:             delta_t,                                   $
              geom_factor:         geom_factor,                               $
              nenergy:             nenergy,                                   $
              nbins:               nangles,                                   $
              bins:                replicate(1b, nenergy, nangles),           $
              energy:              energy_arr,                                $
              denergy:             denergy_arr,                               $
              theta:               theta,                                     $
              phi:                 phi,                                       $
              dtheta:              dtheta,                                    $
              dphi:                dphi,                                      $
              k1:                  k1,                                        $
              k2:                  k2,                                        $
              data:                data,                                      $
              scale:               fltarr(nenergy,nangles),                   $
              phase_inst:         basic_data.header(ind_s:ind_e).phase_instr, $
              sensitivity:        sensitivity,                                $
              op_mode:            op_mode,                                    $
              phase_instr:        basic_data.header(ind_s:ind_e).phase_instr,  $
              gf: gf, $
              mass: 0.010438871 $
            }

  basic_data=0
  RETURN, retdata
END
