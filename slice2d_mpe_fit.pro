;+
;NAME:			slice2d
;PURPOSE:		creates a velocity or energy spectrogram
;			with v perp and v para as x and y, and the
;			specified units as z (color axis).
;			DOESN'T fold data points; generates true 2D slice
;			through the B-V plane
;CALL:			ex: slice2d,get_el('20:31'),[keywords]
;KEYWORDS:		ANG: the angle of the wedge
;			FINISHED: makes the output publication quality
;			when using ps.
;			ZHEIGHT: the zrange to use.  Def is 100 km/s
;			XRANGE: vector specifying the xrange
;			RANGE: vector specifying the color range
;			UNITS: specifies the units ('eflux','df',etc.)
;                              (Def. is 'df')
;			NOZLOG: specifies a linear Z axis
;			THEBDATA: specifies b data to use (def is B3_gse)
;			VAR_LA: vector of tplot variables to show on plot
;			POSITION: positions the plot using a 4-vector
;			ERANGE: specifies the energy range to be used
;			NOFILL: doesn't fill the contour plot with colors
;			NLINES: says how many lines to use if using
;			        NOFILL (def 60, max 60)
;			NOOLINES: suppresses the black contour lines
;			NUMOLINES: how many black contour lines 
;                                  (def. 20, max 60)
;			SHOWDATA: plots all the data points over the contour
;			PLOTENERGY: plots using energy instead of velocity
;			VEL: tplot variable containing the velocity data
;			     (default is calculated with v_3d)
;			NOGRID: forces no triangulation
;			NOCROSS: suppresses cross section line plots
;			RESOLUTION: resolution of the mesh (default is 51)
;			NOSMOOTH: suppresses smoothing
;			NOSUN: suppresses the sun direction line
;			RMBINS: removes the sun noise by cutting out 
;                               certain bins
;			THETA: specifies the theta range for RMBINS (def 20)
;			PHI: specifies the phi range for RMBINS (def 40)
;			NR: removes background noise from ph using noise_remove
;			NOISELEVEL: background level in eflux
;			BOTTOM: level to set as min eflux for
;			        background. def. is 0. 
;			SR, RS, RM2: removes the sun noise using subtraction
;				REQUIRES write_ph.doc to run
;				Note: automatically sets /nosmooth
;				for smoothing, set /smooth
;			NLOW: used with rm2.  Sets bottom of eflux noise level
;				def. 1e4
;			M: marks the tplot at the current time
;			NOVELLINE: suppresses the velocity line
;			LOGLINES: makes the cross sections plot with log
;			CUT_PARA: the value of vpara to make the 1d
;			          cut of vperp 
;			CUT_PERP: the value of vperp to make the 1d
;			          cut of vpara 
;			CUT_BULK_VEL: get 1-D cut in the bulk velocity frame
;			V2: takes a 3-vector velocity and puts it on the plot
;LAST MODIFIED:		2001-11-20 by Tai Phan
;CREATED BY:		Arjun Raj
;EXAMPLES:
;			slice2d_mpe,data,thebdata='B_xyz_codif',
;			vel='v_cod',xrange=[-100,100] 
;			slice2d_mpe,data,thebdata='B_xyz_codif',
;			vel='v_cod',xrange=[-100,100],/cut_bulk_vel ;
;			1-D cut is done in the bulk frame 
;			slice2d_mpe,data,thebdata='B_xyz_codif',
;			vel='v_cod',xrange=[-100,100],cut_para=50,cut_perp=50 
;                       1-D cut is done at specified vpara and vperp

; modified by Jing:
; 2. auto set find the top point of 2d distribution as the bulk vel 
; 3. auto set the zrange as max of f and max/100 as the min

PRO slice2d_mpe_fit, $
  thedata2, $
  nofix = nofix, $
  nsteps = nsteps, $
  xrange = xrange, $
  range = range, $
  units = units, $
  nozlog = nozlog, $
  zlog = zlog, $
  thebdata = thebdata, $
  b3 = b3, $
  position = position, $
  erange = erange, $
  nofill = nofill, $
  var_label = var_label, $
  nlines = nlines, $
  showdata = showdata, $
  plotenergy = plotenergy, $
  vel = vel, $
  nosubtract = nosubtract, $
  nosmooth = nosmooth, $
  smooth = smooth, $
  nocross = nocross, $
  cross = cross, $
  nogrid = nogrid, $
  grid = grid, $
  resolution = resolution, $
  pos2 = pos2, $
  nosun = nosun, $
  sundir = sundir, $
  olines = olines, $
  noolines = noolines, $
  numolines = numolines, $
  setpos = setpos, $
  leavezero = leavezero, $
  rmbins = rmbins, $
  nr = nr, $
  noiselevel = noiselevel, $
  bottom = bottom, $
  zheight = zheight, $
  theta = theta, $
  phi = phi, $
  m = m, $
  rm2 = rm2, $
  nlow = nlow, $
  angle = angle, $
  phb = phb, $
  filename = filename, $
  novelline = novelline, $
  subtract = subtract, $
  double = double, $
  outfile = outfile, $
  cut_perp = cut_perp, $
  cut_para = cut_para, $
  rs = rs, $
  outcuts = outcuts, $
  sr = sr, $
  loglines = loglines, $
  oldlog = oldlog, $
  logplot = logplot, $
  finished = finished, $
  plotlabel = plotlabel, $
  v2 = v2, $
  cut_bulk_vel = cut_bulk_vel, $
  onecnt = onecnt, $
  circ = circ, $
  _EXTRA = e, $
  gsexy = gsexy, $
  gsexz = gsexz, $
  gseyz = gseyz, $
  xout1 = xout1, $
  xout2 = xout2, $
  xout3 = xout3, $
  xout4 = xout4, $
  yout1 = yout1, $
  yout2 = yout2, $
  yout3 = yout3, $
  yout4 = yout4

;-------------------------------------------------------------------
; Check if data are in 'Counts' units
;-------------------------------------------------------------------
IF thedata2.units_name NE 'Counts' THEN BEGIN
    PRINT, 'the 3D data must be in Counts'
    RETURN
ENDIF ELSE BEGIN
    PRINT, 'the 3D data is in Counts'
ENDELSE

;-------------------------------------------------------------------
; IF codif data and velocity not supplied -> return
;-------------------------------------------------------------------
IF thedata2.project_name EQ 'CLUSTER CODIF' THEN BEGIN
    IF NOT KEYWORD_SET(vel) THEN BEGIN
        PRINT, 'you must supply bulk velocity via VEL keyword'
        RETURN
    ENDIF
ENDIF

;-------------------------------------------------------------------
; Set character size
;-------------------------------------------------------------------
!p.charsize = 0.9

;-------------------------------------------------------------------
; filename (is not used)
;-------------------------------------------------------------------
IF KEYWORD_SET(phb) THEN filename = 'write_phb.doc'

;-------------------------------------------------------------------
; Protect the original data by renaming the data structure
;-------------------------------------------------------------------
thedata = thedata2

;-------------------------------------------------------------------
; SR, RS, RM2: removes the sun noise using subtraction
;-------------------------------------------------------------------
IF KEYWORD_SET(sr) THEN rm2 = 1
IF KEYWORD_SET(rs) THEN rm2 = 1
IF KEYWORD_SET(nofill) THEN noolines = 1


;-------------------------------------------------------------------
; temporary (leavezero seems to make a big difference at the origin.
; The best option seems to be leavezero=1 to not remove zeros)
;-------------------------------------------------------------------
leavezero = 1

;-------------------------------------------------------------------
; M: marks the tplot at the current time
;-------------------------------------------------------------------
IF KEYWORD_SET(m) THEN $
  new_time, 'cut2d', thedata.time

;-------------------------------------------------------------------
; Related to outstring ?
;-------------------------------------------------------------------
numperrow = 4

;-------------------------------------------------------------------
; MODIFICATIONS TO MAKE COMMAND LINE SMALLER
;-------------------------------------------------------------------
IF NOT KEYWORD_SET(zheight) AND NOT KEYWORD_SET(angle) THEN angle = 1

IF NOT KEYWORD_SET(nozlog) THEN zlog = 1
IF NOT KEYWORD_SET(nogrid) THEN grid = 1
IF NOT KEYWORD_SET(nocross) THEN cross = 1
IF NOT KEYWORD_SET(nosmooth) THEN smooth = 1
IF NOT KEYWORD_SET(noolines) THEN BEGIN
    IF KEYWORD_SET(numolines) THEN olines = numolines ELSE olines = 20
ENDIF
IF NOT KEYWORD_SET(thebdata) THEN thebdata = 'B3_gse'

IF KEYWORD_SET(b3) THEN thebdata = 'B3_gse'
PRINT, 'B field used is '+thebdata

;IF NOT KEYWORD_SET(subtract) THEN nosubtract = 1 ; C.M.
IF NOT KEYWORD_SET(nosun) THEN sundir = 0

IF KEYWORD_SET(zlog) THEN PRINT, 'zl'
IF KEYWORD_SET(grid) THEN PRINT, 'grid'
IF KEYWORD_SET(cross) THEN PRINT, 'cross'
IF KEYWORD_SET(smooth) THEN PRINT, 'smooth'
IF KEYWORD_SET(olines) THEN PRINT, 'olines'

IF NOT KEYWORD_SET(nsteps) THEN nsteps = 18
IF NOT KEYWORD_SET(units) THEN units = 'DIST FUNC'
IF NOT KEYWORD_SET(nlines) THEN nlines = 60

;-------------------------------------------------------------------
; 
;-------------------------------------------------------------------
perpsym = byte(94)
perpsymbol = STRING(perpsym)

parasym = byte(47)
parasymbol = STRING(parasym) + STRING(parasym)
;-------------------------------------------------------------------
; changing units to s^-3 m^-6
;-------------------------------------------------------------------
;IF KEYWORD_SET(finished) AND units eq 'df' THEN $

IF thedata.project_name NE 'CLUSTER CIS/HIA' THEN BEGIN
    thedata.data = thedata.data / 1000.
ENDIF
;-------------------------------------------------------------------
; 
;-------------------------------------------------------------------
IF KEYWORD_SET(finished) AND KEYWORD_SET(plotlabel) AND $
  !d.name eq 'PS' THEN BEGIN
;	device,/bold
;	xyouts, 0.0,.95,plotlabel+'!N!7',/normal
ENDIF
;-------------------------------------------------------------------
; Load color table for ps output
;-------------------------------------------------------------------
;IF !d.name eq 'PS' THEN loadct,39
;-------------------------------------------------------------------
; Set resolution
 ;-------------------------------------------------------------------
IF NOT KEYWORD_SET(resolution) THEN resolution = 51
IF resolution MOD 2 EQ 0 THEN resolution = resolution + 1
;------------------------------------------------------------------
; Adjust output and labels
;-------------------------------------------------------------------
oldplot = !p.multi

IF KEYWORD_SET(cross) THEN BEGIN ;and  !d.name ne 'PS' THEN BEGIN
    !p.multi = [0, 2, 1]
    grid = 1
ENDIF

IF NOT KEYWORD_SET(position) THEN BEGIN
    x_size = !d.x_size & y_size = !d.y_size
    xsize = .77
    yoffset = 0.
    d = 1.
    IF KEYWORD_SET(cross) THEN BEGIN
        yoffset = yoffset + .5
        xsize = xsize/2.+.13/1.5
        y_size = y_size/2.
        x_size = x_size/2.
        d = .5
        IF y_size LE x_size THEN $
          pos2 = [.13*d+.05, $
                  .03+.13*d, $
                  .05+.13*d + xsize * y_size/x_size, $
                  .13*d + xsize+.03] $
        ELSE $
          pos2 = [.13*d+.05, $
                  .03+.13*d, $
                  .05+.13*d + xsize, $
                  .13*d + xsize *x_size/y_size+.03]
        
    ENDIF
    IF y_size LE x_size THEN $
      position = [.13*d+.05, $
                  .13*d+yoffset, $
                  .05+.13*d + xsize * y_size/x_size, $
                  .13*d + xsize + yoffset] $
    ELSE $
      position = [.13*d+.05, $
                  .13*d+yoffset, $
                  .05+.13*d + xsize, $
                  .13*d + xsize *x_size/y_size + yoffset]
ENDIF ELSE BEGIN
    IF NOT KEYWORD_SET(pos2) THEN BEGIN
        pos2 = position
        pos2(0) = position(0)
        pos2(2) = position(2)
        pos2(3) = position(1)-.08
        pos2(1) = .1
    ENDIF
ENDELSE

IF KEYWORD_SET(var_label) $
  AND NOT KEYWORD_SET(setpos) THEN BEGIN ; AND !d.name ne 'PS' THEN BEGIN
    IF KEYWORD_SET(cross) THEN BEGIN
        pos2(1) = pos2(1) + .04
        pos2(3) = pos2(3) + .04
    ENDIF
    position(1) = position(1) + .02
    position(3) = position(3) + .02
ENDIF

;-------------------------------------------------------------------
; Set up the structure for the one count level
;-------------------------------------------------------------------
IF KEYWORD_SET(ONECNT) THEN BEGIN
    theonecnt = thedata
    FOR i = 0, theonecnt.nenergy-1 DO BEGIN
        theonecnt.data(i, *) = 1./88.0
    ENDFOR
    theonecnt = conv_units(theonecnt, units)
    IF theonecnt.units_name EQ 'Counts' THEN theonecnt.data(*, *) = 1.
ENDIF
;-------------------------------------------------------------------
; Convert Counts -> DF
;-------------------------------------------------------------------

thedata    =  conv_units(thedata, units)

;**********************************************
;bad_bins=where((thedata.dphi eq 0) or (thedata.dtheta eq 0) or $
;	((thedata.data(0,*) eq 0.) AND (thedata.theta(0,*) eq 0.) AND $
;	(thedata.phi(0,*) eq 180.)),n_bad)
;good_bins=where(((thedata.dphi ne 0) AND (thedata.dtheta ne 0)) AND NOT $
;	((thedata.data(0,*) eq 0.) AND (thedata.theta(0,*) eq 0.) AND $
;	(thedata.phi(0,*) eq 180.)),n_good)

;if n_bad ne 0 THEN PRINT,'There are bad bins'


IF thedata.valid ne 1 THEN BEGIN
    PRINT, 'Not valid data'
    RETURN
ENDIF

;bad120 = where(good_bins eq 120,count)
;if count eq 1 AND thedata.data_name eq 'Pesa High' THEN BEGIN
;	PRINT, 'Fixing bad 120 bin'
;	if n_bad eq 0 THEN bad_bins = [120] ELSE bad_bins = [bad_bins,120]
;	good_bins = good_bins(where(good_bins ne 120))
;	n_bad = n_bad + 1
;	n_good = n_good -1
;ENDIF



;**********************************************

;-------------------------------------------------------------------
; get the magnetic field into a variable
;-------------------------------------------------------------------
get_data, thebdata, data = mgf


;************EXPERIMENTAL INTERPOLATION FIX************
;get_data,thebdata,data = bdata
;index = where(bdata.x LE thedata.time + 600 AND bdata.x ge thedata.time - 600)
;store_data,thebdata+'cut',data={x:bdata.x(index),y:bdata.y(index,*)}
;********

;-------------------------------------------------------------------
; get the time information into a variable
;-------------------------------------------------------------------
store_data, 'time', data = {x:thedata.time+thedata.integ_t*.5}

;PRINT, thedata.integ_t, ' Thedata.integ_t'
;interpolate,'time',thebdata+'cut','Bfield'
;get_data,'Bfield',data = mgf
;bfield = FLTARR(3)
;bfield[0] = mgf.y(0,0)
;bfield[1] = mgf.y(0,1)
;bfield[2] = mgf.y(0,2)

;-------------------------------------------------------------------
; get the average value of magnetic field data for interval
;-------------------------------------------------------------------
get_data, thebdata, data = bf_data
IF N_ELEMENTS(bf_data.y) EQ 3 THEN $
  bfield = bf_data.y $
ELSE $
  bfield = dat_avg(thebdata, thedata.time, thedata.end_time)

;-------------------------------------------------------------------
; Start DF calculation.
                                ; In order to find out how many particles there are at all the
                                ; different locations, we must transform the data into cartesian
                                ; coordinates.
                                ;-------------------------------------------------------------------
totalx = FLTARR(1) & totaly = FLTARR(1) & totalz = FLTARR(1)
ncounts = FLTARR(1)

                                ;-------------------------------------------------------------------
                                ; Keyword: ERANGE
                                ;-------------------------------------------------------------------
IF NOT KEYWORD_SET(erange) THEN BEGIN
    erange = [thedata.energy(thedata.nenergy-1, 0), thedata.energy(0, 0)]
    erange = [min(thedata.energy), max(thedata.energy)]
    eindex = indgen(thedata.nenergy)
ENDIF ELSE BEGIN
    eindex = where(thedata.energy(*, 0) ge erange(0) AND $
                   thedata.energy(*, 0) LE erange(1), ct)
    IF ct LE 0 THEN stop
    erange = [min(thedata.energy(eindex, 0)), max(thedata.energy(eindex, 0))]
ENDELSE

mass = thedata.mass / 6.2508206e24

;-------------------------------------------------------------------
; For each energy and angle bin calculate the corresponding velocity
; and DF value
;-------------------------------------------------------------------
FOR i = 0, thedata.nenergy-1 do begin
    currbins = WHERE(thedata.bins(i, *) NE 0 AND $
                     thedata.energy(i, *) LE erange(1) AND $
                     thedata.energy(i, *) GE erange(0) AND $
                     FINITE(thedata.data(i, *)) EQ 1 $
                     , nbins)  
    IF nbins ne 0 THEN BEGIN
        x = FLTARR(nbins) & y = FLTARR(nbins) & z = FLTARR(nbins)
        
        sphere_to_cart, 1, $    ; convert from spherical to cartesian
                        REFORM(thedata.theta(i, currbins)), $
                        REFORM(thedata.phi(i, currbins)), x, y, z
        
        IF NOT KEYWORD_SET(plotenergy) THEN BEGIN
            totalx = [totalx, x * $
                      REFORM(sqrt(2*1.6e-19*thedata.energy(i, currbins)/mass))]
            totaly = [totaly, y * $
                      REFORM(sqrt(2*1.6e-19*thedata.energy(i, currbins)/mass))]
            totalz = [totalz, z * $
                      REFORM(sqrt(2*1.6e-19*thedata.energy(i, currbins)/mass))]
        ENDIF ELSE BEGIN
            totalx = [totalx, x * REFORM(thedata.energy(i, currbins))]
            totaly = [totaly, y * REFORM(thedata.energy(i, currbins))]
            totalz = [totalz, z * REFORM(thedata.energy(i, currbins))]
        ENDELSE
        
        ncounts = [ncounts, REFORM(thedata.data(i, currbins))]
    ENDIF
ENDFOR

totalx = totalx(1:*)
totaly = totaly(1:*)
totalz = totalz(1:*)
ncounts = ncounts(1:*)

;*****HERES SOMETHING NEW
newdata = {dir:FLTARR(N_ELEMENTS(totalx), 3), n:FLTARR(N_ELEMENTS(totalx))}



IF KEYWORD_SET(GSEXY) OR $
  KEYWORD_SET(GSEXZ) OR $
  KEYWORD_SET(GSEYZ) THEN BEGIN
    
    ft = DBLARR(N_ELEMENTS(ncounts)) + thedata.time
    fd = [[totalx], [totaly], [totalz]]
    fs = INTARR(N_ELEMENTS(ncounts)) + thedata.sensitivity
    fp = INTARR(N_ELEMENTS(ncounts)) + thedata.phase_inst
    
    datastr = {x:ft, y:fd}
    store_data, 'dd', dat = datastr, $
                dlim = {inst_num:0, sens:fs, sat:thedata.sat, $
                        phase_instr:fp}

    cis_coord_trans, DATA_IN = 'dd', TRANS = 'CODIF->GSE', $
                     DATA_OUT = 'dd_new'
    get_data, 'dd_new', data = d, dlim = dlim
    store_data, 'dd_new', /DELETE
    store_data, 'dd', /DELETE
    
    totalx = REFORM(d.y(*, 0))
    totaly = REFORM(d.y(*, 1))
    totalz = REFORM(d.y(*, 2))
    
    get_data, vel, data = bv    ; transform the bulk velocity in gse
    datastr = {x:bv.x, y:bv.y}
    store_data, 'dd', dat = datastr, $
                dlim = {inst_num:0, sens:fs, sat:thedata.sat, $
                        phase_instr:fp}
    cis_coord_trans, DATA_IN = 'dd', TRANS = 'CODIF->GSE', $
                     DATA_OUT = vel
    
ENDIF

newdata.dir(*, 0) = totalx
newdata.dir(*, 1) = totaly
newdata.dir(*, 2) = totalz
newdata.n = ncounts

IF KEYWORD_SET(nosubtract) THEN PRINT, 'No velocity transform' ELSE BEGIN
    IF KEYWORD_SET(vel) THEN $
      PRINT, 'Velocity used for subtraction is '+vel $
    ELSE $
      PRINT, 'Velocity used for subtraction is V_3D'
ENDELSE


IF KEYWORD_SET(vel) THEN BEGIN
    PRINT, 'Using '+vel+' for velocity vector'
    
;	get_data,vel,data = dummy, index = theindex
;	IF theindex eq 0 THEN BEGIN
;		PRINT, 'Loading velocity data....'
;		get_3dt,'v_3d','ph',/nr,/rm2
;	ENDIF

;	interpolate,'time',vel,'value'
;	get_data,'value',data = thevalue
;	thevel = 1000.* REFORM(thevalue.y)

    thevel = 1000. * dat_avg(vel, thedata.time, thedata.end_time)
    
    IF KEYWORD_SET(plotenergy) THEN $
      factor = sqrt(total(thevel(*)^2))*mass/2./1.6e-19 ELSE factor = 1.
    
ENDIF ELSE BEGIN
    PRINT, 'Calculating V with v_3d...'
    
    original_thedata2 = thedata2
    
;	thedata_eflux= $
;    convert_codif_units(original_thedata2,$
;                        "EFLUX",'codif_ts_eff_corr',0,packet=1)
    
    thedata_eflux = conv_units(original_thedata2, "EFLUX")
    
    
;	thevel = 1000. * v_3d_cis(thedata_eflux)
    thevel = 1000. * v_3d(thedata_eflux)
    
    
    IF KEYWORD_SET(plotenergy) THEN $
      factor = sqrt(total(thevel(*)^2))*mass/2./1.6e-19 ELSE factor = 1.
ENDELSE


IF NOT KEYWORD_SET(nosubtract) THEN BEGIN
    newdata.dir(*, 0) = newdata.dir(*, 0) - thevel(0)*factor
    newdata.dir(*, 1) = newdata.dir(*, 1) - thevel(1)*factor
    newdata.dir(*, 2) = newdata.dir(*, 2) - thevel(2)*factor
ENDIF ELSE BEGIN
    newdata.dir(*, 0) = newdata.dir(*, 0)
    newdata.dir(*, 1) = newdata.dir(*, 1)
    newdata.dir(*, 2) = newdata.dir(*, 2)
ENDELSE

;**************NOW CONVERT TO THE DATA SET REQUIRED*****************

                                ;-------------------------------------------------------------------
                                ; Calculate the rotation matrix in which B parallel is on the X axis
                                ; and the bulk velocity vector is in the X-Y plane
                                ;-------------------------------------------------------------------

IF NOT KEYWORD_SET(GSEXY) AND $
  NOT KEYWORD_SET(GSEXZ) AND $
  NOT KEYWORD_SET(GSEYZ) THEN BEGIN
    rot = cal_rot(bfield, thevel)
ENDIF

                                ;-------------------------------------------------------------------
                                ; Apply the rotation to the velocity vectors of each angle bin
                                ; and calculate Vpar, Vperp, DF arrays
                                ;-------------------------------------------------------------------
IF NOT KEYWORD_SET(GSEXY) AND $
  NOT KEYWORD_SET(GSEXZ) AND $
  NOT KEYWORD_SET(GSEYZ)THEN BEGIN
    newdata.dir = newdata.dir#rot
ENDIF

IF KEYWORD_SET(plotenergy) THEN factor = 1. ELSE factor = 1000.

IF KEYWORD_SET(GSEXY) THEN BEGIN ; GSE X-Y
                                ;vperp = $
                                ;  (newdata.dir(*,1)^2 + newdata.dir(*,2)^2)^.5*$
                                ;  newdata.dir(*,1)/abs(newdata.dir(*,1))/factor
    vperp = newdata.dir(*, 1) / factor

    vpara = newdata.dir(*, 0)/factor
    zdata = newdata.n
ENDIF ELSE BEGIN
    IF KEYWORD_SET(GSEXZ) THEN BEGIN ; GSE X-Z
                                ;vperp = $
                                ;  (newdata.dir(*,1)^2 + newdata.dir(*,2)^2)^.5*$
                                ;  newdata.dir(*,2)/abs(newdata.dir(*,2))/factor
        vperp = newdata.dir(*, 2) / factor
        
        vpara = newdata.dir(*, 0)/factor
        zdata = newdata.n
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(GSEYZ) THEN BEGIN ; GSE Y-Z
                                ;vperp = $
                                ;  (newdata.dir(*,0)^2 + newdata.dir(*,2)^2)^.5*$
                                ;  newdata.dir(*,2)/abs(newdata.dir(*,2))/factor
            vperp = newdata.dir(*, 2) / factor

            vpara = newdata.dir(*, 1)/factor
            zdata = newdata.n
        ENDIF ELSE BEGIN        ; Vpara - Vperp
                                ;vperp = $
                                ;  (newdata.dir(*,1)^2 + newdata.dir(*,2)^2)^.5*$
                                ;  newdata.dir(*,1)/abs(newdata.dir(*,1))/factor
            vperp = newdata.dir(*, 1) / factor
            vpara = newdata.dir(*, 0) / factor
            zdata = newdata.n
        ENDELSE
    ENDELSE
ENDELSE

                                ;-------------------------------------------------------------------
                                ; Keyword: ANGLE
                                ;-------------------------------------------------------------------
IF NOT KEYWORD_SET(angle) THEN BEGIN
    
    zmag = newdata.dir(*, 2)/factor
    index = where(abs(zmag) le zheight, count)
    IF count ne 0 THEN BEGIN
        vperp = vperp(index)
        vpara = vpara(index)
        zdata = zdata(index)
    ENDIF ELSE BEGIN
        message, 'NO DATA POINTS AT THAT ZHEIGHT!'
        RETURN
    ENDELSE
    PRINT, 'zheight = ', zheight
    
ENDIF ELSE BEGIN

    IF angle eq 1 THEN angle = 20.
    
    IF KEYWORD_SET(GSEXY)THEN BEGIN
        zmag = newdata.dir(*, 2)/factor
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(GSEXZ) THEN BEGIN
            zmag = newdata.dir(*, 1)/factor
        ENDIF ELSE BEGIN
            IF KEYWORD_SET(GSEYZ) THEN BEGIN
                zmag = newdata.dir(*, 0)/factor
            ENDIF ELSE BEGIN
                zmag = newdata.dir(*, 2)/factor
            ENDELSE
        ENDELSE
    ENDELSE

    
    r = sqrt(vperp(*)^2 + vpara(*)^2)
    
    eachangle = atan(zmag/r)

    index = where(abs(eachangle)/!dtor le angle, count)
    IF count ne 0 THEN BEGIN
        vperp = vperp(index)
        vpara = vpara(index)
        zdata = zdata(index)
    ENDIF ELSE BEGIN
        PRINT, 'NO DATA POINTS AT THAT ANGLE!'
        RETURN
    ENDELSE
    PRINT, 'angle = ', angle
ENDELSE

                                ;-------------------------------------------------------------------
                                ; Keyword: SUNDIR
                                ;-------------------------------------------------------------------
IF KEYWORD_SET(sundir) THEN BEGIN
    sund = [1, 0, 0]
    sund = sund#rot
    vperpsun = (sund(1)^2 + sund(2)^2)^.5*sund(1)/abs(sund(1))
    vparasun = sund(0)
ENDIF

IF NOT KEYWORD_SET(vel) THEN BEGIN
    
    original_thedata2 = thedata2
    
;thedata_eflux= $
;convert_codif_units(original_thedata2,"EFLUX",'codif_ts_eff_corr',0,packet=1)
    
    thedata_eflux = conv_units(original_thedata2, "EFLUX")
    
    PRINT, 'computing veldir using v_3d'
    
;veldir = v_3d_cis(thedata_eflux) 
    veldir = v_3d(thedata_eflux) 
    
ENDIF ELSE BEGIN
    veldir = thevel/1000.
ENDELSE


PRINT, 'here 2******'

IF NOT KEYWORD_SET(GSEXY) AND $
  NOT KEYWORD_SET(GSEXZ) AND $
  NOT KEYWORD_SET(GSEYZ)THEN BEGIN
    veldir = veldir#rot
ENDIF

                                ;-------------------------------------------------------------------
                                ; EXPERIMENTAL GET RID OF 0 THING
                                ;-------------------------------------------------------------------
IF NOT KEYWORD_SET(leavezero) THEN BEGIN
    index = where(zdata ne 0)
    vperp = vperp(index)
    vpara = vpara(index)
    zdata = zdata(index)
ENDIF ELSE PRINT, 'Zeros left in plot'

                                ;-------------------------------------------------------------------
                                ; MAKE SURE THERE ARE NO NEGATIVE VALUES!!
                                ;-------------------------------------------------------------------
index2 = where(zdata lt 0., count)
IF count ne 0 THEN PRINT, 'THERE ARE NEGATIVE DATA VALUES'

index = where(zdata ge 0, count)
IF count ne 0 THEN BEGIN
    vperp = vperp(index)
    vpara = vpara(index)
    zdata = zdata(index)
ENDIF

;******************NOW TO PLOT THE DATA********************

                                ;-------------------------------------------------------------------
                                ; Get the X-axis range
                                ;-------------------------------------------------------------------
IF NOT KEYWORD_SET(xrange) THEN BEGIN
    themax = max(abs([vperp, vpara]))
    xrange = [-1*themax, themax]
ENDIF ELSE themax = max(abs(xrange))

                                ;-------------------------------------------------------------------
                                ; Get the Z-axis range
                                ;-------------------------------------------------------------------
IF NOT KEYWORD_SET(range) THEN BEGIN
    IF NOT KEYWORD_SET(xrange) THEN BEGIN	
        maximum = max(zdata)
        minimum = min(zdata(where(zdata ne 0)))
    ENDIF ELSE BEGIN
        maximum = $
          max(zdata(where(abs(vperp) le themax AND abs(vpara) le themax)))
        minimum = maximum/100    ; changed by Jing from: 
;        min(zdata(where(zdata ne 0 AND abs(vperp) le themax AND $
;                        abs(vpara) le themax)))
    ENDELSE
ENDIF ELSE BEGIN
    maximum = range(1)
    minimum = range(0)
ENDELSE 

IF KEYWORD_SET(zlog) THEN $
  thelevels = 10.^(indgen(nlines)/float(nlines)* $
                   (alog10(maximum) - alog10(minimum)) + alog10(minimum)) $
ELSE $
  thelevels = (indgen(nlines)/float(nlines)*(maximum-minimum)+minimum)

                                ;-------------------------------------------------------------------
                                ; EXTRA STUFF FOR THE CONTOUR LINE OVERPLOTS
                                ;-------------------------------------------------------------------
IF KEYWORD_SET(olines) THEN BEGIN
    IF KEYWORD_SET(zlog) THEN $
      thelevels2 = 10.^(indgen(olines)/float(olines)* $
                        (alog10(maximum) - alog10(minimum)) + $
                        alog10(minimum)) $
    ELSE $
      thelevels2 = (indgen(olines)/float(olines)*(maximum-minimum)+minimum)
    
ENDIF

                                ;-------------------------------------------------------------------
                                ; Set number of colors
                                ;-------------------------------------------------------------------
thecolors = round((indgen(nlines)+1)*(!d.table_size-9)/nlines)+7

IF NOT KEYWORD_SET(nofill) THEN fill = 1 ELSE fill = 0

                                ;-------------------------------------------------------------------
                                ; Set labels
                                ;-------------------------------------------------------------------
IF NOT KEYWORD_SET(finished) THEN BEGIN
    IF NOT KEYWORD_SET(plotenergy) THEN BEGIN
        xtitle = 'V Para (km/sec)'
        ytitle = 'V Perp (km/sec)'
        IF KEYWORD_SET(GSEXY) THEN BEGIN
            xtitle = 'Vx GSE (km/sec)'
            ytitle = 'Vy GSE (km/sec)'
        ENDIF
        IF KEYWORD_SET(GSEXZ) THEN BEGIN
            xtitle = 'Vx GSE (km/sec)'
            ytitle = 'Vz GSE (km/sec)'
        ENDIF
        IF KEYWORD_SET(GSEYZ) THEN BEGIN
            xtitle = 'Vy GSE (km/sec)'
            ytitle = 'Vz GSE (km/sec)'
        ENDIF
    ENDIF ELSE BEGIN
        xtitle = 'E Para (eV)'
        ytitle = 'E Perp (eV)'
        IF KEYWORD_SET(GSEXY) THEN BEGIN
            xtitle = 'Ex GSE (eV)'
            ytitle = 'Ey GSE (eV)'
        ENDIF
        IF KEYWORD_SET(GSEXZ) THEN BEGIN
            xtitle = 'Ex GSE (eV)'
            ytitle = 'Ey GSE (eV)'
        ENDIF
        IF KEYWORD_SET(GSEYZ) THEN BEGIN
            xtitle = 'Ey GSE (eV)'
            ytitle = 'Ez GSE (eV)'
        ENDIF
    ENDELSE
ENDIF ELSE BEGIN
    IF NOT KEYWORD_SET(plotenergy) THEN BEGIN
        xtitle = 'V!19!D'+parasymbol+'!N!7 (km/sec)'
        ytitle = 'V!19!D'+perpsymbol+'!N!7 (km/sec)'
    ENDIF ELSE BEGIN
        xtitle = 'E!19!D'+parasymbol+'!N!7 (eV)'
        ytitle = 'E!19!D'+perpsymbol+'!N!7 (eV)'
    ENDELSE
ENDELSE


                                ;-------------------------------------------------------------------
                                ; Keyword: GRID
                                ;-------------------------------------------------------------------
IF KEYWORD_SET(grid) THEN BEGIN

    x = findgen(resolution)/(resolution-1)*(xrange(1)-xrange(0)) + xrange(0)
    spacing = (xrange(1)-xrange(0))/(resolution-1)
    triangulate, vpara, vperp, tr, b
    thesurf = trigrid(vpara, vperp, zdata, tr, [spacing, spacing], $
                      [xrange(0), xrange(0), xrange(1), xrange(1)], $
                      xgrid = xg, ygrid = yg )
    IF KEYWORD_SET(smooth) THEN thesurf = SMOOTH(thesurf, 3)
    IF N_ELEMENTS(xg) mod 2 ne 1 THEN $
      PRINT, 'The line plots are invalid', N_ELEMENTS(xg)
    
                                ;-----------------------------------------------------------------
                                ; EXPERIMENTAL THINGS HERE. Keyword: LOGPLOT
                                ;-----------------------------------------------------------------
    IF KEYWORD_SET(logplot) THEN BEGIN
        
        vpara2 = vpara
        vperp2 = vperp
        
        magnitude = .5 * alog(vpara^2 + vperp^2)
        
        vpara2 = vpara / sqrt(vpara^2 + vperp^2)
        vperp2 = vperp / sqrt(vpara^2 + vperp^2)
        
        vpara2 = vpara2 * magnitude
        vperp2 = vperp2 * magnitude
        
;CONVERT BACK AS A CHECK
;magnitude = exp(sqrt(vpara2^2 + vperp2^2) )
;vpara3 = vpara2 / sqrt(vpara2^2 + vperp2^2)
;vperp3 = vperp2 / sqrt(vpara2^2 + vperp2^2)
;vpara4 = vpara3 * magnitude
;vperp4 = vperp3 * magnitude

        xrangeold = xrange
        
        xrange(0) = -alog(abs(xrange(0)))
        xrange(1) = alog(xrange(1))
        
        x = findgen(resolution)/(resolution-1)*(xrange(1)-xrange(0)) + xrange(0)
        spacing = (xrange(1)-xrange(0))/(resolution-1)
        triangulate, vpara2, vperp2, tr, b
        thesurf = trigrid(vpara2, vperp2, zdata, tr, [spacing, spacing], $
                          [xrange(0), xrange(0), xrange(1), xrange(1)], $
                          xgrid = xg, ygrid = yg )
        IF KEYWORD_SET(smooth) THEN thesurf = smooth(thesurf, 3)
        IF N_ELEMENTS(xg) mod 2 ne 1 THEN $
          PRINT, 'The line plots are invalid', N_ELEMENTS(xg)
                                ;PRINT,N_ELEMENTS(xg)

        xrange = xrangeold
        
        indexminus = where(xg lt 0.)
        indexplus = where(xg gt 0.)
        
        xg(indexminus) = -exp(abs(xg(indexminus)))
        xg(indexplus) = exp(xg(indexplus))
        yg(indexminus) = -exp(abs(yg(indexminus)))
        yg(indexplus) = exp(yg(indexplus))
        
    ENDIF

;********************************************************
;********************************************************
                                ;-----------------------------------------------------------------
                                ; Set specie & s/c labels
                                ;-----------------------------------------------------------------
    IF thedata.project_name eq 'CLUSTER CIS/HIA' THEN BEGIN
        label_species = 'all ions'
        spacecraft = thedata.sat
    ENDIF ELSE BEGIN
        spacecraft = thedata.sat
        IF thedata.species eq 0 THEN label_species = 'H+'
        IF thedata.species eq 1 THEN label_species = 'He++'
        IF thedata.species eq 2 THEN label_species = 'He+'
        IF thedata.species eq 3 THEN label_species = 'O+'
    ENDELSE
    
    timetitle = 'Sat '+strtrim(spacecraft, 2)+' ' + $
                thedata.project_name+' '+label_species+' ('+ $
                thedata.data_name+')'+'!c'+time_string(thedata.time) + $
                '->' + STRMID(time_string(thedata.end_time), 11, 8)
    
;IF KEYWORD_SET(finished) AND KEYWORD_SET(plotlabel) THEN timetitle = '!B
;IF KEYWORD_SET(finished) AND KEYWORD_SET(plotlabel) THEN timetitle = '!B'
    
                                ;-----------------------------------------------------------------
                                ; Plot contour for GRID
                                ;-----------------------------------------------------------------

    CONTOUR, thesurf, xg, yg, $
             /closed, levels = thelevels, c_color = thecolors, fill = fill, $
             title = timetitle, $
             ystyle = 1, $
             ticklen = -0.01, $
             xstyle = 1, $
             xrange = xrange, $
             yrange = xrange, $
             xtitle = xtitle, $
             ytitle = ytitle, position = position
    IF KEYWORD_SET(olines) THEN BEGIN
        IF !d.name eq 'PS' THEN somecol = !p.color ELSE somecol = 0
        CONTOUR, thesurf, xg, yg, /closed, levels = thelevels2, ystyle = 1+4, $
                 xstyle = 1+4, xrange = xrange, yrange = xrange, $
                 ticklen = 0, /noerase, position = position, col = somecol
    ENDIF

ENDIF  ELSE BEGIN               ; NO GRID
                                ;-----------------------------------------------------------------
                                ; Plot contour for NO GRID
                                ;-----------------------------------------------------------------

    CONTOUR, zdata, vpara, vperp, /irregular, $
             /closed, levels = thelevels, c_color = thecolors, fill = fill, $
             title = timetitle, $
             ystyle = 1, $
             ticklen = -0.01, $
             xstyle = 1, $
             xrange = xrange, $
             yrange = xrange, $
             xtitle = xtitle, $
             ytitle = ytitle, position = position
    IF KEYWORD_SET(olines) THEN BEGIN
        IF !d.name eq 'PS' THEN somecol = !p.color ELSE somecol = 0
        CONTOUR, zdata, vpara, vperp, /irregular, /closed, levels = thelevels2, $
                 ystyle = 1+4, xstyle = 1+4, ticklen = 0, xrange = xrange, $
                 yrange = xrange, position = position, /noerase, col = somecol
    ENDIF
ENDELSE                         ; END GRID

IF NOT KEYWORD_SET(cut_para) THEN cut_para = 0.
IF NOT KEYWORD_SET(cut_perp) THEN cut_perp = 0.
IF KEYWORD_SET(cut_bulk_vel) THEN BEGIN
    cut_para = veldir(0)
    cut_perp = veldir(1)
ENDIF
;------------
; put cut at the top of distribution func,  by Jing
;---------------
index = where(thesurf EQ max(thesurf))
loc_top = array_indices(thesurf, index)
cut_para = xg(loc_top(0))
cut_perp = yg(loc_top(1))
;-----------
oplot, [cut_para, cut_para], xrange, linestyle = 2, thick = 2
oplot, xrange, [cut_perp, cut_perp], linestyle = 2, thick = 2
;oplot,[0,0],xrange,linestyle = 1
;oplot,xrange,[0,0],linestyle = 1

IF KEYWORD_SET(sundir) THEN $
  oplot, [0, vparasun*max(xrange)], [0, vperpsun*max(xrange)], linestyle = 2


                                ;-------------------------------------------------------------------
                                ; Not used
                                ;-------------------------------------------------------------------
IF KEYWORD_SET(vel2) THEN BEGIN
    
    vel2 = vel2#rot
    vperpvel2 = (vel2(1)^2 + vel2(2)^2)^.5*vel2(1)/abs(vel2(1))
    vparavel2 = vel2(0)
    
    bbbb = findgen(36)*(!pi*2/32.)
    usersym, 1.5*cos(bbbb), 1.5*sin(bbbb), /fill
    
                                ;oplot,[vparavel2],[vperpvel2],psym = 8,col= !d.table_size - 10,symsize =1
    oplot, [vparavel2], [vperpvel2], psym = 8, col = 2, symsize = 1
ENDIF


                                ;-------------------------------------------------------------------
                                ; NOVELLINE: suppresses the velocity line
                                ;-------------------------------------------------------------------
IF NOT KEYWORD_SET(novelline) THEN $
  oplot, [0, veldir(0)], [0, veldir(1)] ;,col= !d.table_size-9

                                ;-------------------------------------------------------------------
                                ; Overplot min and max velocity (energy) circles
                                ;-------------------------------------------------------------------
IF KEYWORD_SET(CIRC) THEN BEGIN
    IF NOT KEYWORD_SET(plotenergy) THEN BEGIN
        circy = sin(findgen(360)*!dtor)*sqrt(2.*1.6e-19*erange(0)/mass)/1000.
        circx = cos(findgen(360)*!dtor)*sqrt(2.*1.6e-19*erange(0)/mass)/1000.
        oplot, circx, circy, thick = 2
        
        circy = sin(findgen(360)*!dtor)*sqrt(2.*1.6e-19*erange(1)/mass)/1000.
        circx = cos(findgen(360)*!dtor) * $
                sqrt(2.*1.6e-19*erange(1)/mass)/1000.
        oplot, circx, circy, thick = 2
    ENDIF ELSE BEGIN
        circy = sin(findgen(360)*!dtor)*erange(0)
        circx = cos(findgen(360)*!dtor)*erange(0)
        oplot, circx, circy, thick = 2
        
        circy = sin(findgen(360)*!dtor)*erange(1)
        circx = cos(findgen(360)*!dtor)*erange(1)
        oplot, circx, circy, thick = 2
    ENDELSE
ENDIF
;-------------------------------------------------------------------
; Set Z axis units
;-------------------------------------------------------------------
IF NOT KEYWORD_SET(finished) THEN $
  thetitle = units_string(thedata.units_name) $
ELSE thetitle = units_finalstring(thedata.units_name)

;-------------------------------------------------------------------
; Keyword: PLOTLABEL
;-------------------------------------------------------------------
IF KEYWORD_SET(plotlabel) THEN $
  xyouts, 0.05, .95, plotlabel+'!N!7', /normal, charsize = 1.5


;IF KEYWORD_SET(zlog) THEN thetitle = thetitle + ' (log)'

;-------------------------------------------------------------------
; Draw colorscale
;-------------------------------------------------------------------
draw_color_scale, range = [minimum, maximum], log = zlog, yticks = 10, title = thetitle
;-------------------------------------------------------------------
; Keyword: SHOWDATA
;-------------------------------------------------------------------
IF KEYWORD_SET(showdata) THEN oplot, vpara, vperp, psym = 1


IF KEYWORD_SET(var_label) THEN BEGIN
    get_data, 'time', data = theparticulartime
                                ;PRINT, 'Vars interp. to '+time_string(theparticulartime.x)

    IF KEYWORD_SET(vel) THEN velname = vel ELSE velname = 'V_3D'
    
    IF KEYWORD_SET(nosubtract) THEN $
      outstring = velname + ' used (no trans)' $
    ELSE $
      outstring = vel + 'used (trans)'
    FOR i = 0, N_ELEMENTS(var_label)-1 do begin
        PRINT, var_label(i)
        
;get_data,var_label(i),data = bdata
;index = $
;  where(bdata.x le thedata.time + 600 AND bdata.x ge thedata.time - 600,count)
;IF count ne 0 THEN BEGIN
;store_data,var_label(i)+'cut',data={x:bdata.x(index),y:bdata.y(index,*)}
;interpolate,'time',var_label(i)+'cut','value'
;get_data,'value',data = thevalue
;outstring = $
;   outstring + '  ' + var_label(i) +'= '+ $
;   strtrim(STRING(format = '(G11.4)',thevalue.y(0)),2)
;ENDIF

        thevalue = dat_avg(var_label(i), thedata.time, thedata.end_time)
        
        outstring = outstring + '  ' + var_label(i) +'= '+ $
                    strtrim(STRING(format = '(G11.4)', thevalue), 2)
        
        IF i mod numperrow eq 1 THEN outstring = outstring + '!c' 
    ENDFOR
    xyouts, .13, .06, outstring, /normal
ENDIF

;-------------------------------------------------------------------
; Plot the cross section
;-------------------------------------------------------------------
IF KEYWORD_SET(cross) THEN BEGIN
    n_elem = N_ELEMENTS(thesurf(*, 0))   
    
    IF NOT KEYWORD_SET(cut_perp) THEN perpval = n_elem/2 ELSE BEGIN
        ind = where(xg ge cut_perp)
        IF (xg(ind(0)) - cut_perp) le (cut_perp - xg(ind(0)-1) ) THEN $
          perpval = ind(0) ELSE perpval = ind(0)-1	
    ENDELSE
    
    IF NOT KEYWORD_SET(cut_para) THEN paraval = n_elem/2 ELSE BEGIN
        ind = where(xg ge cut_para)
        IF (xg(ind(0)) - cut_para) le (cut_para - xg(ind(0)-1) ) THEN $
          paraval = ind(0) ELSE paraval = ind(0)-1	
    ENDELSE
    
;	IF KEYWORD_SET(zlog) THEN thetitle = thetitle + ' (log)'

    IF KEYWORD_SET(plotenergy) THEN BEGIN
        xtitle = 'Energy (eV)'
        vore = 'E'
    ENDIF ELSE BEGIN
        xtitle = 'Velocity (km/sec)'
        vore = 'V'
    ENDELSE
    

;HERE COMES SOME NEW COLOR STUFF
    IF !d.name eq 'PS' THEN $
      thecolors = round((indgen(4)+1)*(!d.table_size-9)/4)+7 ELSE BEGIN
        thecolors = indgen(4)
        thecolors = thecolors + 3
    ENDELSE

    IF KEYWORD_SET(double) THEN BEGIN
                                ;first plot vpara on the + side
        plot, xg, [reverse(thesurf(n_elem/2:*, perpval)), $
                   thesurf(n_elem/2+1:*, perpval)], $
              xstyle = 1, ystyle = 1, $
              xrange = xrange, yrange = [minimum, maximum], ylog = zlog, $
              title = 'Cross Sections', xtitle = xtitle, ytitle = thetitle, $
              position = pos2
                                ;overplot vpara on the minus side
        oplot, xg, [thesurf(0:n_elem/2, perpval), $
                    reverse(thesurf(0:n_elem/2-1, perpval))], color = thecolors(0)
                                ;now vperp on the + side
        oplot, xg, [reverse(REFORM(thesurf(paraval, n_elem/2+1:*))), $
                    REFORM(thesurf(paraval, n_elem/2:*))], color = thecolors(1)
                                ;and now vper on the - side
        oplot, xg, [REFORM(thesurf(paraval, 0:n_elem/2)), $
                    reverse(REFORM(thesurf(paraval, 0:n_elem/2-1)))], $
               color = thecolors(2)
        
    ENDIF ELSE BEGIN

        IF NOT KEYWORD_SET(loglines) THEN BEGIN
                                ;first plot vpara on the + side
            plot, xg(n_elem/2:*)-cut_para, thesurf(n_elem/2:*, perpval), $
                  xstyle = 1, ystyle = 1, $
                  xrange = xrange, yrange = [minimum, maximum], ylog = 0, $
                  title = 'Cross Sections', xtitle = xtitle, ytitle = thetitle, $
                  position = pos2
                                ;overplot vpara on the minus side
            oplot, xg(0:n_elem/2)-cut_para, thesurf(0:n_elem/2, perpval)
                                ;,color = thecolors(0)
                                ;now vperp on the + side
            oplot, xg(n_elem/2:*)-cut_perp, REFORM(thesurf(paraval, n_elem/2:*)), $
                   color = thecolors(1)
                                ;and now vperp on the - side
            oplot, xg(0:n_elem/2)-cut_perp, REFORM(thesurf(paraval, 0:n_elem/2)), $
                   color = thecolors(1)
            
            xout1 = xg(n_elem/2:*)-cut_para
            xout2 = xg(0:n_elem/2)-cut_para
            xout3 = xg(n_elem/2:*)-cut_perp
            xout4 = xg(0:n_elem/2)-cut_perp

            yout1 = thesurf(n_elem/2:*, perpval)
            yout2 = thesurf(0:n_elem/2, perpval)
            yout3 = REFORM(thesurf(paraval, n_elem/2:*))
            yout4 = REFORM(thesurf(paraval, 0:n_elem/2))
            
        ENDIF  ELSE BEGIN

;*****PUT IN LOGLINES STUFF HERE********
;********EXPERIMENTAL THINGS HERE************

            IF NOT KEYWORD_SET(oldlog) THEN BEGIN       
                vpara2 = vpara
                vperp2 = vperp
                
                magnitude = .5 * alog(vpara^2 + vperp^2)
                
                vpara2 = vpara / sqrt(vpara^2 + vperp^2)
                vperp2 = vperp / sqrt(vpara^2 + vperp^2)
                
                vpara2 = vpara2 * magnitude
                vperp2 = vperp2 * magnitude         
;CONVERT BACK AS A CHECK
;magnitude = exp(sqrt(vpara2^2 + vperp2^2) )
;vpara3 = vpara2 / sqrt(vpara2^2 + vperp2^2)
;vperp3 = vperp2 / sqrt(vpara2^2 + vperp2^2)
;vpara4 = vpara3 * magnitude
;vperp4 = vperp3 * magnitude
                xrangeold = xrange
                
                xrange(0) = -alog(abs(xrange(0)))
                xrange(1) = alog(xrange(1))        
;stop
                x = findgen(resolution)/(resolution-1)*(xrange(1)-xrange(0)) + $
                    xrange(0)
                spacing = (xrange(1)-xrange(0))/(resolution-1)
                triangulate, vpara2, vperp2, tr, b
                thesurf = trigrid(vpara2, vperp2, zdata, tr, [spacing, spacing], $
                                  [xrange(0), xrange(0), xrange(1), xrange(1)], $
                                  xgrid = xg, ygrid = yg )
                IF KEYWORD_SET(smooth) THEN thesurf = smooth(thesurf, 3)
                IF N_ELEMENTS(xg) mod 2 ne 1 THEN $
                  PRINT, 'The line plots are invalid', N_ELEMENTS(xg)
                                ;PRINT,N_ELEMENTS(xg)

                xrange = xrangeold
                
                indexminus = where(xg lt 0.)
                indexplus = where(xg gt 0.)
                
                xg(indexminus) = -exp(abs(xg(indexminus)))
                xg(indexplus) = exp(xg(indexplus))
                
            ENDIF

;********************************************************

;plot,xg,[reverse(thesurf(n_elem/2:*,perpval)),thesurf(n_elem/2+1:*,perpval)],$
            plot, xg(n_elem/2:*), thesurf(n_elem/2:*, perpval), $
                  xstyle = 1, ystyle = 1, $
                  xrange = [min(thedata.energy(*, 0)), xrange[1]], $
                  yrange = [minimum, maximum], ylog = zlog, /xlog, $
                  title = 'Cross Sections', xtitle = xtitle, ytitle = thetitle, $
                  position = pos2
                                ;vpara on the minus side
            oplot, xg(n_elem/2:*), $
                   reverse(thesurf(0:n_elem/2, perpval)) ;, color = thecolors(0)
                                ;vperp on the + side
            oplot, xg(n_elem/2:*), $
                   REFORM(thesurf(paraval, n_elem/2:*)), color = thecolors(1)
                                ;vperp on the - side
            oplot, xg(n_elem/2:*), $
                   reverse(REFORM(thesurf(paraval, 0:n_elem/2))), color = thecolors(1)
        ENDELSE 
    ENDELSE
;stop
                                ;put a dotted line
    oplot, [0, 0], [minimum, maximum], linestyle = 1
    IF NOT KEYWORD_SET(plotenergy) THEN BEGIN
;oplot,[sqrt(2.*1.6e-19*erange(0)/mass)/1000., $
;sqrt(2.*1.6e-19*erange(0)/mass)/1000.],[minimum,maximum],linestyle = 5
;oplot,-[sqrt(2.*1.6e-19*erange(0)/mass)/1000., $
;sqrt(2.*1.6e-19*erange(0)/mass)/1000.],[minimum,maximum],linestyle = 5
;oplot,[sqrt(2.*1.6e-19*erange(1)/mass)/1000., $
;sqrt(2.*1.6e-19*erange(1)/mass)/1000.],[minimum,maximum],linestyle = 5
;oplot,-[sqrt(2.*1.6e-19*erange(1)/mass)/1000., $
;sqrt(2.*1.6e-19*erange(1)/mass)/1000.],[minimum,maximum],linestyle = 5
        IF KEYWORD_SET(onecnt) THEN BEGIN
            oplot, sqrt(2.*1.6e-19*theonecnt.energy(*, 0)/mass)/1000., $
                   theonecnt.data(*, 0), color = thecolors(3), linestyle = 3
            oplot, -sqrt(2.*1.6e-19*theonecnt.energy(*, 0)/mass)/1000., $
                   theonecnt.data(*, 0), color = thecolors(3), linestyle = 3
        ENDIF
    ENDIF ELSE BEGIN
;		oplot,[erange(0),erange(0)],[minimum,maximum],linestyle = 5
;		oplot,-[erange(0),erange(0)],[minimum,maximum],linestyle = 5
;		oplot,[erange(1),erange(1)],[minimum,maximum],linestyle = 5
;		oplot,-[erange(1),erange(1)],[minimum,maximum],linestyle = 5
        IF KEYWORD_SET(onecnt) THEN BEGIN
            oplot, theonecnt.energy(*, 0), theonecnt.data(*, 0), $
                   color = thecolors(3), linestyle = 3
            oplot, -theonecnt.energy(*, 0), theonecnt.data(*, 0), $
                   color = thecolors(3), linestyle = 3
        ENDIF
    ENDELSE
    
                                ;now put the titles on the side of the graph
    positions = -findgen(5)*(pos2(3)-pos2(1))/5 + pos2(3)-.03
    
    IF KEYWORD_SET(GSEXY) THEN BEGIN
        xyouts, pos2(2) + .03, positions(1), vore+'x GSE', /norm, charsize = 1.01 ;.5
        xyouts, pos2(2) + .03, positions(2), vore+'y GSE', /norm, $
                color = thecolors(1), charsize = 1.01 ;.5
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(GSEXZ) THEN BEGIN
            xyouts, pos2(2) + .03, positions(1), vore+'x GSE', /norm, charsize = 1.01 ;.5
            xyouts, pos2(2) + .03, positions(2), vore+'z GSE', /norm, $
                    color = thecolors(1), charsize = 1.01 ;.5
        ENDIF ELSE BEGIN
            IF KEYWORD_SET(GSEYZ) THEN BEGIN
                xyouts, pos2(2) + .03, positions(1), vore+'y GSE', /norm, charsize = 1.01 ;.5
                xyouts, pos2(2) + .03, positions(2), vore+'z GSE', /norm, $
                        color = thecolors(1), charsize = 1.01 ;.5
            ENDIF ELSE BEGIN
                xyouts, pos2(2) + .03, positions(1), vore+' para', /norm, charsize = 1.01 ;.5
                xyouts, pos2(2) + .03, positions(2), vore+' perp', /norm, $
                        color = thecolors(1), charsize = 1.01 ;.5
            ENDELSE
        ENDELSE
    ENDELSE
    
    IF KEYWORD_SET(onecnt) THEN $
      xyouts, pos2(2) + .03, positions(4), 'One count', /norm, $
              color = thecolors(3), charsize = .5
ENDIF

;stop

IF KEYWORD_SET(outfile) THEN BEGIN
    openw, thefile, outfile, /get_lun
    PRINTF, thefile, time_string(thedata.time)
    
    IF NOT KEYWORD_SET(cut_perp) THEN perpval = n_elem/2 ELSE BEGIN
        ind = where(xg ge cut_perp)
        IF (xg(ind(0)) - cut_perp) le (cut_perp - xg(ind(0)-1) ) THEN $
          perpval = ind(0) ELSE perpval = ind(0)-1	
    ENDELSE
    
    IF NOT KEYWORD_SET(cut_para) THEN paraval = n_elem/2 ELSE BEGIN
        ind = where(xg ge cut_para)
        IF (xg(ind(0)) - cut_para) le (cut_para - xg(ind(0)-1) ) THEN $
          paraval = ind(0) ELSE paraval = ind(0)-1	
    ENDELSE
    
    
    filedata = FLTARR(3, N_ELEMENTS(xg))
    
    filedata(0, *) = xg
    filedata(1, *) = thesurf(*, perpval)
    filedata(2, *) = thesurf(paraval, *)
    
    PRINTF, thefile, filedata
    close, /all
    
ENDIF

;********EXTRA PART*********

;IF NOT KEYWORD_SET(cut_perp) THEN perpval = n_elem/2 ELSE BEGIN
;ind = where(xg ge cut_perp)
;IF (xg(ind(0)) - cut_perp) le (cut_perp - xg(ind(0)-1) ) THEN $
;perpval = ind(0) ELSE perpval = ind(0)-1	
;ENDELSE

;IF NOT KEYWORD_SET(cut_para) THEN paraval = n_elem/2 ELSE BEGIN
;ind = where(xg ge cut_para)
;IF (xg(ind(0)) - cut_para) le (cut_para - xg(ind(0)-1) ) THEN $
;paraval = ind(0) ELSE paraval = ind(0)-1	
;ENDELSE


;filedata = FLTARRx(3,N_ELEMENTS(xg))

;filedata(0,*) = xg
;filedata(1,*) = thesurf(*,perpval)
;filedata(2,*) = thesurf(paraval,*)
;outcuts = filedata

;********END EXTRA PART*******
IF !d.name NE 'PS' THEN !p.multi = oldplot

END
