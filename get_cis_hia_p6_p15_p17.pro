;+
;FUNCTION:      get_cis_hia_p6, prod_num, sat, time=time
;INPUT:
;       prod_num: Product number
;	sat: Satellite number (1-4)
;	time: If set, return time array of data values (index ignored)
;
;PURPOSE:   Retrieve CIS/HIA 3D data.
;
;CREATED BY:    Peter Schroeder - Modified by C. mouikis
;
;LAST MODIFICATION:  09/12/01
;  09/12/01 - The correction of phi angle as a function of energy is
;             included
;  01/29/02 - Bug fix. In the phi correction routine and for 31 energy
;             products the code should loop over 32 bins and not 31
;             (although the 32nd is not used)
;  01/29/02 - Keyword NO_PHI_COR introduced.
;  04/15/02 - addpted from p6 MF
;  08/07/02 - give calibration time MF
;  14/01/03 - angle calc changed MF
;  15/01/03 - sw_stop MF
;  06/03/03 - gf MF
;-
FUNCTION get_cis_hia_p6_p15_p17, prod_num, sat,  time=time, $
                         frencheff=frencheff

  COMMON get_error, get_err_no, get_err_msg, default_verbose


  get_err_no = 0
  get_err_msg= 'Data found in time interval'

  ;--------------------------------------------------------------------
  ; Number of energies and angles for each 3D product
  ;--------------------------------------------------------------------

  CASE prod_num OF
    6: enang = [31,88]
    15: enang = [16,88]
    17: enang = [62,88]
  ENDCASE

  nenergy = enang[0]
  nangles = enang[1]
  anode_map=cis_hia_anode_map(88)

  ;--------------------------------------------------------------------
  ; Get the data
  ;--------------------------------------------------------------------
  options = [long(sat),1l,long(prod_num),0l]

  IF keyword_set(time) THEN RETURN,get_cis_data_times(options)

  basic_data = get_cis_all_data(options)
;  packets = N_ELEMENTS(basic_data.data)
  IF NOT KEYWORD_SET(basic_data) THEN BEGIN ; Test if the product exists
    get_err_no = 1
    get_err_msg = 'Product '+ $
      string(prod_num,format='(i2.2)') + $
      ' does not exist'
    PRINT, get_err_msg
    RETURN, 0
  ENDIF

  ;--------------------------------------------------------------------
  ; Prepare data arrays
  ;--------------------------------------------------------------------
  dtime = cis_tplot_time(basic_data.header.time_in_ms) ; time -> tplot-time
  data = float(basic_data.data.data)

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
  theta_arr = fltarr(nangles)
  phi_arr = fltarr(nangles)
  dphi_arr = fltarr(nangles)
  
  ntheta = 8
  nphi = 4
  get_cis_hia_angles,ntheta,nphi,theta,phi,dtheta,dphi,an=anode_map,calt=dtime(0),sat=sat
  FOR i = 0,nphi-1 DO BEGIN
    theta_arr[i] = theta[0]
    phi_arr[i] = phi[i]
    dphi_arr[i] = dphi[i]
  ENDFOR
  
  ntheta = 8
  nphi = 8
  get_cis_hia_angles,ntheta,nphi,theta,phi,dtheta,dphi,an=anode_map,calt=dtime(0),sat=sat
  FOR i = 0,nphi-1 DO BEGIN
    theta_arr[4+i] = theta[1]
    phi_arr[4+i] = phi[i]
    dphi_arr[4+i] = dphi[i]
  ENDFOR
  
  ntheta = 8
  nphi = 16
  get_cis_hia_angles,ntheta,nphi,theta,phi,dtheta,dphi,an=anode_map,calt=dtime(0),sat=sat
  FOR j = 0,3 DO BEGIN
    FOR i = 0,nphi-1 DO BEGIN
      theta_arr[12 + i + j*nphi] = theta[j+2]
      phi_arr[12 + i + j*nphi] = phi[i]
      dphi_arr[12 + i + j*nphi] = dphi[i]
    ENDFOR
  ENDFOR
  
  ntheta = 8
  nphi = 8
  get_cis_hia_angles,ntheta,nphi,theta,phi,dtheta,dphi,an=anode_map,calt=dtime(0),sat=sat
  FOR i = 0,nphi-1 DO BEGIN
    theta_arr[76+i] = theta[6]
    phi_arr[76+i] = phi[i]
    dphi_arr[76+i] = dphi[i]
  ENDFOR
  
  ntheta = 8
  nphi = 4
  get_cis_hia_angles,ntheta,nphi,theta,phi,dtheta,dphi,an=anode_map,calt=dtime(0),sat=sat
  FOR i = 0,nphi-1 DO BEGIN
    theta_arr[84+i] = theta[7]
    phi_arr[84+i] = phi[i]
    dphi_arr[84+i] = dphi[i]
  ENDFOR

  theta=theta_arr
  phi=phi_arr
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

  ;--------------------------------------------------------------------
  ; For all records get the geometric factor, energies and number of 
  ; spins
  ;--------------------------------------------------------------------
  sysin=SYSTIME(1)
  ; the following line masks the solar wind phi sectors for use with the
  ; sw_sweep_stop, phi limits are currently only experimental, but results
  ; agree with CL software.
  sw_phi_index=WHERE(phi(0,*) GT 190 AND phi(0,*) LT 220)

  hvtblset=TOTAL(TAG_NAMES(basic_data.data) EQ 'HVTBL') EQ 1
  FOR j = 0,ndatapoints-1 DO BEGIN
    IF KEYWORD_SET(frencheff) THEN BEGIN
      geom_factor[j] = get_cis_hia_geom_factor(dtime[j], sat, sensitivity[j])
    ENDIF ELSE BEGIN
      geom_factor[j] = jim_geom_factor(sat, sensitivity[j])
    ENDELSE
    IF hvtblset THEN hvtbl=basic_data.data(j).hvtbl ELSE hvtbl=0

    IF hvtbl EQ -32768 THEN hvtbl = -1 ; Jing: since hvtbl is not recognized as the "bad hvtbl" for program "get_cis_hia_energies".

    get_cis_hia_energies, $
      dtime[j], $
      sat, $
      nenergy, $
      op_mode[j], $
      hvtbl, $
      energy, $
      denergy, $
      sensitivity[j], $
      sw_stop_index
    
    ;------------------------------------------------------------------
    ; Check compressed/uncompressed mode and product numbers
    ;------------------------------------------------------------------
    new_prod_num = prod_num

    ;------------------------------------------------------------------------
    ; get number of spins accumulated from calibration files
    IF tlm_rate(j) NE 0 THEN BEGIN ;jing: add this situation for the case tlm_rate(j) eq 0
        num_spins = get_cis_hia_acc_spin(dtime[j], sat, tlm_rate[j], $
                                         op_mode[j], $
                                         new_prod_num, $
                                         basic_data)
    ENDIF ELSE str_element, basic_data.data, 'real_acc_spins', num_spins ; jing
        
;----------------
;Jing: due to function get_cis_hia_acc_spin, when there is no spin to
;calc, it will automatically read from the basic_data, which will
;result in num_spins is not the data at the time point "j". Instead,
;it's the num_spins of all time points, which means num_spins become
;an array of same dimention as delta_t. In consideration of this
;case, I added the following one line.
    IF n_elements(num_spins) EQ n_elements(delta_t) THEN num_spins = num_spins(j) 
;-------------------
    delta_t[j] = num_spins*pspin[j]/1d3 ; delta time 

    FOR i = 0,nangles-1 DO BEGIN
      energy_arr[*,i,j] = energy
      IF sw_stop_index GT 0 AND TOTAL(sw_phi_index EQ i) GT 0 THEN energy_arr[sw_stop_index+1:*,i,j]=0.0
     denergy_arr[*,i,j] = denergy
    ENDFOR
  ENDFOR
  ; now corect geometry factor for number of anodes summed
  ; number of energy levels corrected onbaord
  geom_factor=geom_factor*(16./ntheta);*(62./nenergy)

  ;--------------------------------------------------------------------
  ; Get efficiency
  ;--------------------------------------------------------------------
  gf = cis_hia_efficiency(dtime, $
                          sat, $
                          energy_arr, $
                          theta, $
                          anode_map, $
                          frencheff=frencheff)
  MESSAGE,'Calibration time:'+STRTRIM(systime(1)-sysin,2),/inf

;----------------------------------------------------------------------
; Construct the data structure
;----------------------------------------------------------------------
  retdata = {                                                                 $
              project_name:        'CLUSTER CIS/HIA',                         $
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
              k1:                  INTARR(ndatapoints),                       $
              k2:                  INTARR(ndatapoints),                       $
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




