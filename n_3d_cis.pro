;+
;FUNCTION:	n_3d_cis(dat,ENERGY=en,ERANGE=er,EBINS=ebins,
;                        ANGLE=an,ARANGE=ar,BINS=bins)
;INPUT:	
;	dat:	structure,	2d data structure filled 
;                               by get_eesa_surv, get_eesa_burst, etc.
;KEYWORDS
;	ENERGY:	fltarr(2),	optional, min,max energy range for integration
;	ERANGE:	fltarr(2),	optional, min,max energy bin numbers
;	                                  for integration
;	EBINS:	bytarr(na),	optional, energy bins array for integration
;					0,1=exclude,include,  
;					na = dat.nenergy
;	ANGLE:	fltarr(2,2),	optional, angle range for integration
;				theta min,max (0,0),(1,0) -90<theta<90 
;				phi   min,max (0,1),(1,1)   0<phi<360 
;	ARANGE:	fltarr(2),	optional, min,max angle bin numbers
;	                                  for integration
;	BINS:	bytarr(nb),	optional, angle bins array for integration
;					0,1=exclude,include,  
;					nb = dat.ntheta
;	BINS:	bytarr(na,nb),	optional, energy/angle bins array for 
;                                         integration
;					0,1=exclude,include
;PURPOSE:
;	Returns the density, n, 1/cm^3
;NOTES:	
;	Function normally called by "get_3dt" or "get_2dt" to 
;	generate time series data for "tplot.pro".
;
;CREATED BY:
;	J.McFadden	95-7-27	
;LAST MODIFICATION:
;	96-7-6		J.McFadden	added more keywords
;       20/02/02     MF  energy array bug removed
;       28/01/03     MF comments on units added
;       06/03/03     enfac moved to loop MF
;       09/24/2003      C. Mouikis  Introduced normalization
;                                   factor when a limit angle bin range is usde
;       02/18/2004   bug related to the calculation of omega_tot corrected
;-
FUNCTION n_3d_cis,dat2, $
                  ENERGY=en,ERANGE=er,EBINS=ebins,$
                  ANGLE=an,ARANGE=ar,BINS=bins
  
;-----------------------------------------------------------------------  
; Check for valid data  
;-----------------------------------------------------------------------  
  density = 0.
  IF dat2.valid EQ 0 THEN BEGIN
    print,'Invalid Data'
    RETURN, density
  ENDIF
;-----------------------------------------------------------------------  

  dat    = dat2
  na     = dat.nenergy ; number of energy bins
  nb     = dat.nbins ; number of angle bins
  ebins2 = replicate(1b,na) ; used for integration over certain energy bins
  bins2  = replicate(1b,nb) ; used for integration of certain angle bins
  
;-----------------------------------------------------------------------
; min,max energy range for integration. Keyword ENERGY
;-----------------------------------------------------------------------  
  IF keyword_set(en) THEN BEGIN
    ebins2(*)=0
    er2=[energy_to_ebin(dat,en)]
    IF er2(0) GT er2(1) THEN er2=reverse(er2)
    ebins2(er2(0):er2(1))=1
  ENDIF

;-----------------------------------------------------------------------
; min,max energy bin numbers for integration. Keyword ERANGE
;-----------------------------------------------------------------------
  IF keyword_set(er) THEN IF er(0) GE 0 AND er(1) GE 0 THEN BEGIN
    ebins2(*)=0
    er2=er
    IF er2(0) GT er2(1) THEN er2=reverse(er2)
    ebins2(er2(0):er2(1))=1
  ENDIF

;-----------------------------------------------------------------------
; energy bins array for integration. Keyword EBINS
;-----------------------------------------------------------------------
  IF keyword_set(ebins) THEN ebins2=ebins 
  
;-----------------------------------------------------------------------
; angle range for integration. Keyword ANGLE
;-----------------------------------------------------------------------
  IF keyword_set(an) THEN BEGIN
    IF ndimen(an) NE 2 THEN BEGIN
      print,'Error - angle keyword must be (2,2)'
    ENDIF ELSE BEGIN
      bins2=angle_to_bins(dat,an)

;      IF sat EQ 3 THEN BEGIN
;        sc3_half_instrument_time = time_double('2003-02-24/00:00:00')
;        IF dat.time(0) GT  sc3_half_instrument_time THEN bins2(44:87) = 0
;      ENDIF
    ENDELSE
  ENDIF
  
;-----------------------------------------------------------------------
; angle bins numbers for integration. Keyword ARANGE
;---------------------------------------------------help, bins2--------------------
  IF keyword_set(ar) THEN BEGIN ; min,max angle bin numbers for integration
    bins2(*)=0
    IF ar(0) GT ar(1) THEN BEGIN
      bins2(ar(0):nb-1)=1
      bins2(0:ar(1))=1
    ENDIF ELSE BEGIN
      bins2(ar(0):ar(1))=1
    ENDELSE
  ENDIF
  
;-----------------------------------------------------------------------  
; angle bins array for integration
; energy/angle bins array for integration. Keyword BINS
;-----------------------------------------------------------------------
  IF keyword_set(bins) THEN bins2=bins ; 
; stop
  IF ndimen(bins2) NE 2 THEN bins2=ebins2#bins2
  packets = n_elements(dat.time)
  data = DBLARR(dat.nenergy, dat.nbins, packets)
  FOR i = 0, packets-1 DO $
    data(*,*,i) = dat.data(*,*,i) * bins2
;-----------------------------------------------------------------------

  energy = dat.energy(*,*,0)
  denergy = dat.denergy(*,*,0)
  theta = dat.theta/!radeg
  phi = dat.phi/!radeg
  dtheta = dat.dtheta/!radeg
  dphi = dat.dphi/!radeg
  mass = dat.mass    ; units eV/(km/s)^2
  Const = (mass/2)^(.5)*1e-5 ; units sqrt(eV)s/cm

;-----------------------------------------------------------------------  
; Set domega  
;-----------------------------------------------------------------------  
  str_element,dat,"domega",value=domega,index=ind
  IF ind GE 0 THEN BEGIN
    IF ndimen(domega) EQ 1 THEN domega=replicate(1.,na)#domega
  ENDIF ELSE BEGIN
    IF ndimen(dtheta) EQ 1 THEN dtheta=replicate(1.,na)#dtheta
    IF ndimen(dphi) EQ 1 THEN dphi=replicate(1.,na)#dphi
    ; this formula comes from domega=dphi*h where h is the height
    ; h = sin(theta+dtheta/2)-sin(theta-dtheta/2) = 2cos(theta)sin(dtheta/2)
    domega=2.*dphi*cos(theta)*sin(0.5*dtheta)
    ; Penou uses with negligible difference:
    ; domega=dphi*cos(theta)*dtheta
  ENDELSE

  omega_tot = TOTAL(domega(er2(0),*)*bins2(er2(0),*)) ; if a limit number of angle bins
                                            ; is used a normalization is introduced

;-----------------------------------------------------------------------  
; Calculate density
;-----------------------------------------------------------------------  
  sumdata = DBLARR(na, packets)
  density = DBLARR(packets)
  ; Jim McFadden uses dE/E :
  ;enfac=denergy(*,0)/energy(*,0)/sqrt(energy(*,0)) ; units 1/sqrt(ev)
  ; Penou/Pallochia use 2dV/V instead, which would be identical for continuous V(E)
  ; one can emulate that by:
  ; enfac=2.*(sqrt(energy(*,0)+denergy(*,0)/2)-sqrt(energy(*,0)-denergy(*,0)/2)) $
  ;        /energy(*,0)
  ; but transferring the actual limits which agree with Penou gives closer result for HIA:
  ;
  ; experimental: activate emin,emax return in get_cis_hia_3d for this
  ; enfac=2.*(sqrt(dat.emax(*,0,0))-sqrt(dat.emin(*,0,0)))/energy(*,0) $
  ; Regard, that HIA energies change from time to time in SW modes!

  FOR i = 0, packets - 1 DO BEGIN
    sumdata(*,i) = total(data(*,*,i)*domega,2, /NaN); / (omega_tot / (4 * !pi)) ; Jing
;stop
    sumdata(na-1,i)= sumdata(na-1,i)*(dat.k1(i) EQ 255)
    enfac=dat.denergy(*,0,i)/dat.energy(*,0,i)/sqrt(dat.energy(*,0,i)) ; units 1/sqrt(ev)
    density(i) = Const*total(enfac*sumdata(*,i), /NaN)
  ENDFOR       ; Density units are 1/cm^3

  RETURN, density
END
