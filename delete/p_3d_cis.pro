;+
;FUNCTION:	p_3d_cis(dat,ENERGY=en,ERANGE=er,EBINS=ebins,ANGLE=an,ARANGE=ar,BINS=bins)
;INPUT:	
;	dat:	structure,	2d data structure filled by get_eesa_surv, get_eesa_burst, etc.
;	pack:   packet number
;KEYWORDS
;	ENERGY:	fltarr(2),	optional, min,max energy range for integration
;	ERANGE:	fltarr(2),	optional, min,max energy bin numbers for integration
;	EBINS:	bytarr(na),	optional, energy bins array for integration
;					0,1=exclude,include,  
;					na = dat.nenergy
;	ANGLE:	fltarr(2,2),	optional, angle range for integration
;				theta min,max (0,0),(1,0) -90<theta<90 
;				phi   min,max (0,1),(1,1)   0<phi<360 
;	ARANGE:	fltarr(2),	optional, min,max angle bin numbers for integration
;	BINS:	bytarr(nb),	optional, angle bins array for integration
;					0,1=exclude,include,  
;					nb = dat.ntheta
;	BINS:	bytarr(na,nb),	optional, energy/angle bins array for integration
;					0,1=exclude,include
;PURPOSE:
;	Returns the pressure tensor, [Pxx,Pyy,Pzz,Pxy,Pxz,Pyz], eV/cm^3 
;NOTES:	
;	Function normally called by "get_3dt" or "get_2dt" to
;	generate time series data for "tplot.pro".
;
;	This version has been taken from 'p_3d.pro' to work with Equator-S data. 7/21/98
;
;
;ORIGINAL CREATED BY:
;	J.McFadden	95-7-27
;LAST MODIFICATION: 05/22/01
;
;MODIFICATION HISTORY:
;	07/21/98: D.P.P.  added parameter pack which stroess the
;	          current packet number
;		  The routine still acts on one packet of data at a time
;		  but the structure
;		  contains the full data set (data = [E,A,packets])
;       04/14/99: added /NaN to TOTAL (checks for missing data))
;       05/22/01: C.M. Adapted from the EQS software for the codif
;                 needs
;       02/01/03: C.M. Routine is vectorized
;       09/09/03: C.M. Bugs related to density/velocity calculation
;                 are corrected 
;-

FUNCTION p_3d_cis,dat,         $
                  ENERGY=en,   $
                  ERANGE=er,   $
                  EBINS=ebins, $
                  ANGLE=an,    $
                  ARANGE=ar,   $
                  BINS=bins

  p3dxx = 0.
  p3dyy = 0.
  p3dzz = 0.
  p3dxy = 0.
  p3dxz = 0.
  p3dyz = 0.
  
  IF dat.valid EQ 0 THEN BEGIN
    print,'Invalid Data'
    RETURN, [p3dxx,p3dyy,p3dzz,p3dxy,p3dxz,p3dyz]
  ENDIF
  
;FOR pack = 0, dat2.packets-1 DO BEGIN
;	dat = convert_codif_units(dat2,"EFLUX")	; Use Energy Flux
;ENDFOR
  packets = n_elements(dat.time)
  na = dat.nenergy
  nb = dat.nbins
	
  ebins2=replicate(1b,na)
  IF keyword_set(en) THEN BEGIN
    ebins2(*)=0
    er2=[energy_to_ebin(dat,en)]
    IF er2(0) GT er2(1) THEN er2=reverse(er2)
    ebins2(er2(0):er2(1))=1
  ENDIF
  IF keyword_set(er) THEN BEGIN
    ebins2(*)=0
    er2=er
    IF er2(0) GT er2(1) THEN er2=reverse(er2)
    ebins2(er2(0):er2(1))=1
  ENDIF
  IF keyword_set(ebins) THEN ebins2=ebins
  
  bins2=replicate(1b,nb)
  IF keyword_set(an) THEN BEGIN
    IF ndimen(an) NE 2 THEN BEGIN
      print,'Error - angle keyword must be (2,2)'
    ENDIF ELSE BEGIN
      bins2=angle_to_bins(dat,an)
;      IF sat EQ 3 THEN BEGIN
;        sc3_half_instrument_time = time_double('2003-02-23/22:00:00')
;        IF dat.time(0) GT  sc3_half_instrument_time THEN bins(44:87) = 0
;      ENDIF
    ENDELSE
  ENDIF
  IF keyword_set(ar) THEN BEGIN
    bins2(*)=0
    IF ar(0) GT ar(1) THEN BEGIN
      bins2(ar(0):nb-1)=1
      bins2(0:ar(1))=1
    ENDIF ELSE BEGIN
      bins2(ar(0):ar(1))=1
    ENDELSE
  ENDIF
  IF keyword_set(bins) THEN bins2=bins
 ; stop
  IF ndimen(bins2) NE 2 THEN bins2=ebins2#bins2
 ; stop
;data = dat.data(*,*,pack)*bins2 	;7/20/98 DPP

  data = DBLARR(dat.nenergy, dat.nbins, packets)
  
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(bins2(*,0)), $
                               N_ELEMENTS(bins2(0,*))) * bins2)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))))
  
  vdummy = TRANSPOSE(vdummy, [1,2,0])

;  store_data, 'bins2', data = {ebins:ebins2, bins:bins, vdummy:vdummy}

  data = dat.data * vdummy  ; apply the bins2 on the dat.data and put into array data with dimention of [dat.nenergy,dat.nbins,packets)
 ; stop
;  FOR i = 0, packets-1 DO $
;    data(*,*,i) = dat.data(*,*,i) * bins2

  energy = dat.energy
  denergy = dat.denergy
  theta = dat.theta/!radeg
  phi = dat.phi/!radeg
  dtheta = dat.dtheta/!radeg
  dphi = dat.dphi/!radeg
  mass = dat.mass * 1.6e-22   ;1u=1.66e-27
  Const = (mass/(2.*1.6e-12))^(-.5)  
;  domega = dat.domega
;-----------------------------------------------------------------------
; Set domega  
;-----------------------------------------------------------------------
  str_element,dat,"domega",value=domega,index=ind
  IF ind GE 0 THEN BEGIN
    IF ndimen(domega) EQ 1 THEN domega=replicate(1.,na)#domega
  ENDIF ELSE BEGIN
    IF ndimen(dtheta) EQ 1 THEN dtheta=replicate(1.,na)#dtheta
    IF ndimen(dphi) EQ 1 THEN dphi=replicate(1.,na)#dphi
    domega=2.*dphi*cos(theta)*sin(.5*dtheta)
  ENDELSE
  
  cth = cos(theta)
  sth = sin(theta)
  cph = cos(phi)
  sph = sin(phi)
  cth2 = cth^2
  cthsth = cth*sth
  
  dnrg = Const*denergy*(energy^(-.5))
  flux=[0.,0.,0.]
  flux = FLTARR(3,packets)
  vel  = FLTARR(3,packets)
  
;sumxx = total(data*cph*cph*domega*cth2,2)
;sumyy = total(data*sph*sph*domega*cth2,2)
;sumzz = total(data*domega*sth*sth,2)
;sumxy = total(data*cph*sph*domega*cth2,2)
;sumxz = total(data*cph*domega*cthsth,2)
;sumyz = total(data*sph*domega*cthsth,2)

;dnrg = Const*denergy*(energy^(-.5))
;p3dxx = total(dnrg*sumxx)
;p3dyy = total(dnrg*sumyy)
;p3dzz = total(dnrg*sumzz)
;p3dxy = total(dnrg*sumxy)
;p3dxz = total(dnrg*sumxz)
;p3dyz = total(dnrg*sumyz)

;flux=[0.,0.,0.]
;;flux = j_3d_cis(dat2,pack,ENERGY=en,ERANGE=er,EBINS=ebins,ANGLE=an,ARANGE=ar,BINS=bins)
;sumdatax = total(data*cph*domega*cth,2)
;sumdatay = total(data*sph*domega*cth,2)
;sumdataz = total(data*domega*sth,2)
;dnrg=denergy*(energy^(-1))
;flux(0) = total(dnrg*sumdatax)
;flux(1) = total(dnrg*sumdatay)
;flux(2) = total(dnrg*sumdataz)

;;density = n_3d_cis(dat2,pack,ENERGY=en,ERANGE=er,EBINS=ebins,ANGLE=an,ARANGE=ar,BINS=bins)
;sumdata = total(data*domega,2)
;density = ((mass/(2.*1.6e-12))^(.5))*total(denergy*(energy^(-1.5))*sumdata)

;if density eq 0. then begin
;	vel=[0.,0.,0.]
;endif else begin
;	vel = flux/density
;endelse
;p3dxx = mass*(p3dxx-vel(0)*flux(0))/1.6e-12
;p3dyy = mass*(p3dyy-vel(1)*flux(1))/1.6e-12
;p3dzz = mass*(p3dzz-vel(2)*flux(2))/1.6e-12
;p3dxy = mass*(p3dxy-vel(0)*flux(1))/1.6e-12
;p3dxz = mass*(p3dxz-vel(0)*flux(2))/1.6e-12
;p3dyz = mass*(p3dyz-vel(1)*flux(2))/1.6e-12

;	Pressure is in units of eV/cm**3

;if keyword_set(diag) then begin
;print, " This section not tested yet!!!!!"
;if diag eq "diag" then begin

;	p = [[p3dxx,p3dxy,p3dxz],[p3dxy,p3dyy,p3dyz],[p3dxz,p3dyz,p3dzz]]
;	nr_tred2,p,d,e
;	nr_tqli,d,e,p
;	print,"d =",d
;	print,"p =",p(0,0),p(1,1),p(2,2),p(0,1),p(0,2),p(1,2)

;endif
;endif

;return, [p3dxx,p3dyy,p3dzz,p3dxy,p3dxz,p3dyz]

  p3dxx = dblarr(packets)
  p3dyy = dblarr(packets)
  p3dzz = dblarr(packets)
  p3dxy = dblarr(packets)
  p3dxz = dblarr(packets)
  p3dyz = dblarr(packets)
  
  totalp = dblarr(packets, 6)
   
  cp2omct2 = cph * cph * domega * cth2
  sp2omct2 = sph * sph * domega * cth2
  st2om = domega * sth * sth
  cpspomct2 = cph * sph * domega * cth2
  cpomctst = cph * domega * cthsth
  spomctst = sph * domega * cthsth
  
  cpomct = cph * domega * cth
  spomct = sph * domega * cth 
  omst = domega * sth
  
  ; cp2omct2
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(cp2omct2(*,0)), $
                               N_ELEMENTS(cp2omct2(0,*))) * cp2omct2)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )

  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
 
  sumxx = TOTAL(data * vdummy, 2, /NaN)

  ; sp2omct2
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(sp2omct2(*,0)), $
                               N_ELEMENTS(sp2omct2(0,*))) * sp2omct2)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )
  
  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
  
  sumyy = TOTAL(data * vdummy, 2, /NaN)
  
  ; st2om
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(st2om(*,0)), $
                               N_ELEMENTS(st2om(0,*))) * st2om)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )
  
  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
  
  sumzz = TOTAL(data * vdummy, 2, /NaN)
  
  
  ; cpspomct2
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(cpspomct2(*,0)), $
                               N_ELEMENTS(cpspomct2(0,*))) * cpspomct2)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )
  
  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
  
  sumxy = TOTAL(data * vdummy, 2, /NaN)
  
  
  ; cpomctst
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(cpomctst(*,0)), $
                               N_ELEMENTS(cpomctst(0,*))) * cpomctst)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )
  
  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
  
  sumxz = TOTAL(data * vdummy, 2, /NaN)
  
  ; spomctst
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(spomctst(*,0)), $
                               N_ELEMENTS(spomctst(0,*))) * spomctst)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )
  
  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
  
  sumyz = TOTAL(data * vdummy, 2, /NaN)
  
  
  ; cpomct
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(cpomct(*,0)), $
                               N_ELEMENTS(cpomct(0,*))) * cpomct)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )
  
  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
  
  sumdatax = TOTAL(data * vdummy, 2, /NaN)
  
  ; spomct
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(spomct(*,0)), $
                               N_ELEMENTS(spomct(0,*))) * spomct)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )
  
  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
  
  sumdatay = TOTAL(data * vdummy, 2, /NaN)
  
  ; omst
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(omst(*,0)), $
                               N_ELEMENTS(omst(0,*))) * omst)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )
  
  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
  
  sumdataz = TOTAL(data * vdummy,   2, /NaN)
  
  ; domega
  vdummy = TEMPORARY(REPLICATE(1, 1, $
                               N_ELEMENTS(domega(*,0)), $
                               N_ELEMENTS(domega(0,*))) * domega)
  vdummy = TEMPORARY( $
                      REBIN(vdummy, packets, $
                            N_ELEMENTS(vdummy(0,*,0)), $
                            N_ELEMENTS(vdummy(0,0,*))) )
  
  vdummy = TEMPORARY(TRANSPOSE(vdummy, [1,2,0]))
  
  sumdata = TOTAL(data * vdummy, 2, /NaN)
 
  dnrg = Const*REFORM(denergy(*,0,*))*(REFORM(energy(*,0,*))^(-0.5))

  p3dxx = TOTAL(dnrg * sumxx, 1, /NaN)
  p3dyy = TOTAL(dnrg * sumyy, 1, /NaN)
  p3dzz = TOTAL(dnrg * sumzz, 1, /NaN)
  p3dxy = TOTAL(dnrg * sumxy, 1, /NaN)
  p3dxz = TOTAL(dnrg * sumxz, 1, /NaN)
  p3dyz = TOTAL(dnrg * sumyz, 1, /NaN)

  dnrg=REFORM(denergy(*,0,*)) * REFORM((energy(*,0,*)^(-1)))
    
  flux(0,*) = total(dnrg * sumdatax, 1, /NaN)
  flux(1,*) = total(dnrg * sumdatay, 1, /NaN)
  flux(2,*) = total(dnrg * sumdataz, 1, /NaN)
  
  density = FLTARR(packets)
  Const = (mass/(2.*1.6e-12))^(.5)

  FOR ii=0, packets-1 DO BEGIN
    density(ii) = Const * TOTAL(REFORM(denergy(*,*,ii)) * $
                                (REFORM(energy(*,*,ii))^(-1.5)) * $
                                REFORM(sumdata(*,ii)), /NaN)
  ENDFOR
  
  zden = WHERE(density GT 0., czden)
  IF czden GT 0 THEN BEGIN
    vel(0,zden) = flux(0,zden) / density(zden)
    vel(1,zden) = flux(1,zden) / density(zden)
    vel(2,zden) = flux(2,zden) / density(zden)
  ENDIF

  p3dxx = mass*(p3dxx - vel(0,*) * flux(0,*))/1.6e-12
  p3dyy = mass*(p3dyy - vel(1,*) * flux(1,*))/1.6e-12
  p3dzz = mass*(p3dzz - vel(2,*) * flux(2,*))/1.6e-12
  p3dxy = mass*(p3dxy - vel(0,*) * flux(1,*))/1.6e-12
  p3dxz = mass*(p3dxz - vel(0,*) * flux(2,*))/1.6e-12
  p3dyz = mass*(p3dyz - vel(1,*) * flux(2,*))/1.6e-12
  
;	Pressure is in units of eV/cm**3

  IF KEYWORD_SET(diag) THEN BEGIN
    PRINT, " This section not tested yet!!!!!"
    IF diag EQ "diag" THEN BEGIN
      p = [[p3dxx,p3dxy,p3dxz],[p3dxy,p3dyy,p3dyz],[p3dxz,p3dyz,p3dzz]]
      nr_tred2,p,d,e
      nr_tqli,d,e,p
      print,"d =",d
      print,"p =",p(0,0),p(1,1),p(2,2),p(0,1),p(0,2),p(1,2)
      
    ENDIF
  ENDIF
  
  totalp(*, 0) = p3dxx(*)
  totalp(*, 1) = p3dyy(*)
  totalp(*, 2) = p3dzz(*)
  totalp(*, 3) = p3dxy(*)
  totalp(*, 4) = p3dxz(*)
  totalp(*, 5) = p3dyz(*)
;stop
  RETURN, totalp
  
END

