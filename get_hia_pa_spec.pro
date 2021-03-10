;+
; PROCEDURE:	get_codif_pa_spec
;
; PURPOSE:	
;	Generates pitch angle-time spectrogram data structures for tplot
;
; INPUT:		
;
; Created by: Dimple P. Patel (for EQS/ESIC)
;             Space Science Group
;             Unviversity of New Hampshire
;             dpatel@teams.sr.unh.edu
; Date: 6/30/99
; Version: 1.43
; Last Modification: 05/22/01
; Modification History:
;     07/14/99 - removed while loop that computed time averaged 
;                pitch angles 
;     07/22/99 - add keyword 'gap_time' for DATA GAP time check
;     08/19/99 - add keyword 'EFF_FILE'
;     09/23/99 - changed check for keyword FILTER => 0=Yes, 1=No
;     09/30/99 - added keyword 'PaBin'
;     10/01/99 - mulitply data by aweight just after theta's and phi's
;                are selected out.  only these are used in computing
;                aweight.
;     10/07/99 - fixed miscalculation of flux*aweight
;     11/09/99 - added option for backgroud subtraction for all non
;                proton species.
;     04/09/01 - Modified for the needs of CLUSTER/CODIF data
;     05/22/01 - C.M. the keyword PABIN is utilised
;     02/05/02 - C.M. The keyword ALL_ENERGY_BINS is introduced
;     03/13/02 - C.M. OLD_EFF keyword added
;     08/08/02 - C.M. BKG keyword added
;-

PRO get_hia_pa_spec, sat, inst, prod,           $
                       mag_theta, mag_phi,        $
                       specie=specie,             $
                       eff_table,                 $
                       data_str,  	          $
                       ENERGY=energy, 	          $
                       ANGLE=an, 		  $
                       ARANGE=ar, 	          $
                       BINS=bins, 	          $
                       units = units,  	          $
                       name  = name, 	          $
                       missing = missing,         $
                       s_t=s_t,		          $
                       e_t=e_t,		          $
                       PRODUCT = product,	  $
                       EFF_ROUTINE = eff_routine, $
                       FILTER = FILTER,           $
                       DELTAT = deltat,           $
                       EFF_FILE = eff_file,       $
                       PABIN = PaBin,             $
                       BACKGROUND = background,   $
                       ALL_ENERGY_BINS=ALL_ENERGY_BINS, $
                       OLD_EFF=OLD_EFF, $
                       BKG=BKG, $
                       COMBINE=COMBINE
  
  COMMON get_error, get_err_no, get_err_msg, default_verbose
  
  ex_start = systime(1)         ; strat timing execution time
  
;----------------------------------------------------------------------------
; Choose reading routine (CODIF or HIA) - Read data. Variable INST
;----------------------------------------------------------------------------
  CASE inst OF
      0:	routine = 'get_cis_cod_data'
      1:	routine = 'get_cis_hia_data'
      2:        routine = 'get_cis_cod_data'
      ELSE: BEGIN
          print, 'Instrument number must be 0 (CODIF) or 1 (HIA)'
          RETURN
      END
  ENDCASE
  
  IF inst EQ 1 THEN   dat = call_function(routine,prod,sat) $ ; Jing: get_cis_hia_data does not have kspecie
  ELSE  dat = call_function(routine, prod, specie = specie, sat) 
                                ; parameters: PROD, SPECIE, SAT
  
  IF get_err_no GT 0 THEN RETURN
  
;----------------------------------------------------------------------------
; Substract background. Keyword BKG
;----------------------------------------------------------------------------
  IF KEYWORD_SET(bkg) AND   inst NE 1 THEN   dat = sub_bkg_3d(dat,sat,specie)
;----------------------------------------------------------------------------
  
  packets=n_elements(dat.time) ; number of packets
  nenergy = dat.nenergy ; number of energy bins
  nbins = dat.nbins ; number of angle bins
  
;----------------------------------------------------------------------------
; Convert units. Keyword UNITS
;----------------------------------------------------------------------------  
  IF NOT KEYWORD_SET(units) THEN units = 'Counts'
  IF units NE 'Counts' AND inst EQ 0 THEN $
    dat = convert_codif_units(dat,units,'codif_ts_eff_corr', $
                              eff_table, $
                              specie=specie, $
                              packets=packets, old_eff=old_eff)
  
  IF units NE 'Counts' AND inst EQ 2 THEN $
    dat = convert_rpa_units(dat,units,'codif_ts_eff_corr', $
                            eff_table, $
                            specie=specie, $
                            packets=packets, $
                            INCRATES = incrates)
;jing: add the convert units part for hia.
; convert_hia_units is a pro not a function
  IF units NE 'Counts' AND inst EQ 1 THEN $
    convert_hia_units, dat,units

;  stop
;-----------------------------------------------------------------------
; Set domega
;-----------------------------------------------------------------------
  theta = dat.theta/!radeg
  phi = dat.phi/!radeg
  dtheta = dat.dtheta/!radeg
  dphi = dat.dphi/!radeg
  
  str_element,dat,"domega",value=domega,index=ind
  IF ind GE 0 THEN BEGIN
    IF ndimen(domega) EQ 1 THEN domega=replicate(1.,na)#domega
  ENDIF ELSE BEGIN
    IF ndimen(dtheta) EQ 1 THEN dtheta=replicate(1.,na)#dtheta
    IF ndimen(dphi) EQ 1 THEN dphi=replicate(1.,na)#dphi
    domega=2.*dphi*cos(theta)*sin(.5*dtheta)
  ENDELSE
;-----------------------------------------------------------------------
  
  IF NOT(KEYWORD_SET(paBin)) THEN PaBin = 22.50
  missing = !values.f_nan
  packets = n_elements(dat.time)
  t = (dat.time + dat.end_time)/2.
  
  pa = get_pitch_angle(dat, mag_theta, mag_phi,combine=combine)
  
  fldat = REFORM(dat.data)
  padat = REFORM(pa)
  
  IF KEYWORD_SET(paBin) THEN BEGIN ; set PA bins according to PaBin
    IF dat.nbins NE 88 AND dat.nbins NE 128 THEN paRange = 180 ELSE paRange = 360 ;Jing add 128 for pro23
    BinStartStop = 0.
    FOR i = 1, paRange/paBin DO BinStartStop = [BinStartStop, i*paBin]
    BinCenter = 0.
    FOR i = 1, ((PaRange/paBin)*2)-1,2 DO BinCenter = [BinCenter, i*(PaBin/2.)]
    BinCenter = BinCenter(1:*)
    nmax = paRange / PaBin
  ENDIF
  
  time = dblarr(packets)
  nvar = nmax
  var = fltarr(packets,nvar)
  data = dblarr(packets,nvar)
  
  IF KEYWORD_SET(ALL_ENERGY_BINS) THEN BEGIN
    FOR iebin = 0, dat.nenergy-1 DO BEGIN
      
      data = dblarr(dat.nenergy,packets,nvar)
      
      fldat = REFORM(dat.data)
      padat = REFORM(pa)
      
      fldat = fldat(iebin,*,*)
      fldat = REFORM(fldat)
      padat = padat(iebin,*,*)
      padat = REFORM(padat)
      
      FOR n = 0, packets-1 DO BEGIN ; Loop over all time steps
        
        datSort0 = padat(*,n)
        domegaSort = domega(0,*)
        fldatSort = fldat(*,n)
        
        var(n,0:nvar-1) = BinCenter(*)
        
        FOR k = 0, nmax-1 DO BEGIN ; Loop over all pa bins
          
          exi = WHERE(datSort0 LT BinStartStop(k+1) AND $
                      datSort0 GE BinStartStop(k), ct)
          IF exi(0) NE -1 THEN BEGIN
            IF STRUPCASE(units) EQ 'COUNTS' THEN $
              aweight = REPLICATE(1,N_ELEMENTS(exi)) $
            ELSE $
              aweight = domegaSort(exi)
            IF ct(0) GT 1 THEN $
              data(iebin,n,k) = TOTAL((fldatSort(exi)*aweight),/NaN) / $
              TOTAL(aweight, /NaN)$
            ELSE $
              data(iebin,n,k) = (fldatSort(exi)*aweight)/aweight
          ENDIF ELSE $
            data(iebin,n,k) = missing
          
        ENDFOR
        
      ENDFOR
      
    ENDFOR
    
    datastr = {x: (dat.time + dat.end_time)/2., $
               y: data, v:var, e:dat.energy(*,0,0)}
    labels =''
    name = STRMID(name,0,STRPOS(name,'EN')-1) + $
      STRMID(name, STRPOS(name, '_SC'), STRLEN(name))
    store_data, name, $
      data=datastr,dlim={ylog:1,labels:labels,panel_size:2.}
    
  ENDIF

  IF NOT KEYWORD_SET(ALL_ENERGY_BINS) THEN BEGIN
    
    IF KEYWORD_SET(ENERGY) THEN BEGIN
      ;----------------------------------------------------------------
      ; min,max energy range for integration. Keyword ENERGY
      ;----------------------------------------------------------------
      er2=[energy_to_ebin(dat,energy)] ; find min & max energy bin
      IF er2(0) GT er2(1) THEN er2=reverse(er2)
      n_en = er2(1)-er2(0)+1
      
      fldat = TOTAL(fldat(er2(0):er2(1),*,*),1)/n_en
      padat = TOTAL(padat(er2(0):er2(1),*,*),1)/n_en
    ENDIF ELSE BEGIN
      fldat = TOTAL(fldat,1)/dat.nenergy
      padat = TOTAL(padat,1)/dat.nenergy
    ENDELSE
    
    IF NOT(KEYWORD_SET(paBin)) THEN PaBin = 22.50
    missing = !values.f_nan
    packets = n_elements(dat.time)
    t = (dat.time + dat.end_time)/2.
    
    IF KEYWORD_SET(paBin) THEN BEGIN ; set PA bins according to PaBin
      IF dat.nbins NE 88 AND  dat.nbins NE 128 THEN paRange = 180 ELSE paRange = 360
      BinStartStop = 0.
      FOR i = 1, paRange/paBin DO BinStartStop = [BinStartStop, i*paBin]
      BinCenter = 0.
      FOR i = 1, ((PaRange/paBin)*2)-1,2 DO BinCenter = [BinCenter, i*(PaBin/2.)]
      BinCenter = BinCenter(1:*)
      nmax = paRange / PaBin
    ENDIF
    
    time = dblarr(packets)
    nvar = nmax
    var = fltarr(packets,nvar)
    data = dblarr(packets,nvar)
    
    
    FOR n = 0, packets-1 DO BEGIN ; Loop over all time steps
      
      datSort0 = padat(*,n)
      domegaSort = domega(0,*)
      fldatSort = fldat(*,n)
 
      var(n,0:nvar-1) = BinCenter(*)
      FOR k = 0, nmax-1 DO BEGIN ; Loop over all pa bins
        
        exi = WHERE(datSort0 LT BinStartStop(k+1) AND $
                    datSort0 GE BinStartStop(k), ct)
        IF exi(0) NE -1 THEN BEGIN
          IF STRUPCASE(units) EQ 'COUNTS' THEN $
            aweight = REPLICATE(1,N_ELEMENTS(exi)) $
          ELSE $
            aweight = domegaSort(exi)
          IF ct(0) GT 1 THEN $
            data(n,k) = TOTAL((fldatSort(exi)*aweight),/NaN) / $
            TOTAL(aweight, /NaN)$
          ELSE $
            data(n,k) = (fldatSort(exi)*aweight)/aweight
        ENDIF ELSE $
          data(n,k) = missing
      ENDFOR
      
  ENDFOR  
    datastr = {x: (dat.time + dat.end_time)/2., $
               y: data, v:var}
    labels =''
    store_data, name, data=datastr,dlim={ylog:1,labels:labels,panel_size:2.}
    
  ENDIF
END
