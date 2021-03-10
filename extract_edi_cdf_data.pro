;+
; NAME: extract_edi_cdf_data.pro
;
; PURPOSE: Extract Cluster EDI data from a variety of CDF products
;
; CATEGORY: Data extraction
;
; CALLING SEQUENCE: 
;    istat = extract_edi_cdf_data('20010213','3', $
;                                 epoch,time,vvec,evec, $
;                                 path='./',ftype='mp',quality=quality)
;
; INPUTS:
;    '20010213'          8-character date string
;    '3'                 1-character spacecraft number string           
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;    path='your_own_path_to_file'     String, path to the EDI CDF file
;                                     Default:  See code below for default paths
;
;    filename='your_own_filename.cdf'  String, filename of CDF file
;                                     NOTE:  Overrides date, scnum, vers and ftype
;
;    vers='V01'    String, version # of EDI CDF file
;                     Default: Highest version is returned
;
;    ftype =       'pp'   4-sec, prime parameter, no "BAD" data
;                  '1p'   4-sec, all qualities
;                  '2p'   2-sec, all qualities
;                  '4p'   1-sec, all qualities
;                  'mp'   Variable resolution, all qualities (Mark Chutter product)
;                  'ppp'  4-sec, all qualities, has "winner" and "loser" data both
;                         (NOTE:  only "winner" data is returned at this time)
;
;    Quality filters:  "gc" and "qflag" keywords can be used
;                      (One or the other; not at same time)
;
;       gc =    not set    All non-fill data returned regardless of quality
;       gc =    0          Same as not set
;       gc =    1          Only Good and Caution data returned
;                             (i.e. sbytes(0,*) = 1 or 2)
;
;       qflag = not set:   All non-fill data returned regardless of quality
;               0:         Same as not set
;               1:         Only Good and Caution data returned
;                             (i.e. sbytes(0,*) = 1 or 2)
;               2:         Only Good data returned
;                             (i.e. sbytes(0,*) = 2)
;
;    trange=[t0,t1]      All data within the time range [t0,t1] is
;                        returned; t1 and t0 are in decimal seconds since
;                        midnight (ssm)
;
; OUTPUTS:
;    istat         0=Failure, 1=Success
;    epoch         Time in Epoch
;    time          Time in decimal seconds since midnight
;    vvec          EDI drift velocity vector, dimensions [3,ntime], GSE, Inertial frame, km/s
;    evec          EDI electric field vector, dimensions [3,ntime], GSE, Intertial frame, mV/m
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;    fnout=fnout   String, the name of the file that was ultimately opened (no path)
;
;    msg=msg_out   String, message pertaining to failure (for istat=0)
;
;    rchi2=rchi2   Float, dim=[ntime]; the reduced chi-squared of the triangulation analysis
;
;    sbytes=sbytes Integer array of dimension [7,ntime], the status
;                  bytes for every data point:
;                      sbytes(0,*) = data quality (0=bad,1=caution,2=good)
;                      sbytes(1,*) = Percentage of 1keV beams used in
;                                    entire spin (100=all 1kev beams;  0=all
;                                    0.5keV beams;  Else, we're in energy
;                                    switching mode)
;                      sbytes(2,*) = Percentage of Class-A beams used
;                                    in entire spin
;                      sbytes(3,*) = Method Papertrail and Ambiguity
;                                    Flag (see EDI_PISO_OUTPUT_DEFS.txt
;                                    from the EDI_PISO library)
;                      sbytes(4,*) = Percentage of Triangulation
;                                    outliers (only reported for
;                                    triangulation results, not TOF results) 
;                      sbytes(5,*) = Fractional drift step magnitude
;                                    error (percent) in SC frame (0-254)
;                      sbytes(6,*) = Drift step azimuthal error in degrees
;                                    in SC frame (0-254)
;
;    quality=quality  Integer array;  dim=[ntime], data quality (0=bad,1=caution,2=good)
;
;    nbeam=nbeam      Integer array;  dim=[ntime], FOR ftype='ppp' ONLY:  # of beams used in
;                                                  analysis (not
;                                                  defined for other
;                                                  ftypes because it
;                                                  exists only in the
;                                                  ppp files)
;    meth_usd         Integer array;  dim=[ntime], Designates the
;                                                  method finally used
;                                                  for the result:
;                                                  0 = Triangulation
;                                                  1 = Poorman's ToF
;                                                  2 = Simultan ToF
;                                                  3 = Richman's ToF

; MODIFICATION HISTORY: Written Jan. 28th, 2002 by
;                       Pamela A. Puhl-Quinn, ppq@mpe.mpg.de
;                       UPDATED:  March 28th, 2008, pamela.puhlquinn@unh.edu
;                                 April 2nd, 2008, PPQ
;
;-


function extract_edi_cdf_data, date_in, scnum_in, $ ; INPUT
                               epoch, time, vvec, evec, $ ; OUTPUT
                               $ ; OPTIONAL INPUT
                               path=path_in, filename=filename_in, vers=vers_in, ftype=ftype_in, $
                               gc=gc, qflag=qflag, trange=trange_in, verbose=verbose, $
                               $ ; OPTIONAL OUTPUT
                               sbytes=sbytes,$
                               fnout=fnout, $
                               msg=msg_out, $
                               rchi2=rchi2, quality=quality, $
                               nbeam=nbeam, meth_usd=meth_usd

if (n_elements(nbeam) ne 0) then dum = temporary(nbeam) ; This shouldn't be defined yet
verbose = keyword_set(verbose)

; Set time range, in seconds since midnight
if (n_elements(trange_in) eq 0) then trange=[0.,24.]*3600. else trange=trange_in

; Quality flag keywords
gc = keyword_set(gc)
if (gc and n_elements(qflag) ne 0) then begin
    msg_out = 'Use either /gc or qflag=N but not both'
    if (verbose) then message, msg_out, /cont
    return, 0
endif

; Trim date and scnum strings (why?  not sure)
date = strtrim(date_in,2)
scnum = strtrim(scnum_in,2)

; File type, ftype (lowercase)
if (n_elements(ftype_in) eq 0) then ftype='pp' else ftype=strlowcase(ftype_in)
msg_out = ''

; Define the path
; Default path is ftype dependent
case ftype of
    'pp':def_path='/nfs/cluster4/CDF/edi/c'+scnum+'/'
    '1p':def_path='/nfs/cluster4/CDF_all/latest/c'+scnum+'/'
    '2p':def_path='/nfs/cluster4/chunk_2_CDF/c'+scnum+'/'
    '4p':def_path='/nfs/cluster4/chunk_4_CDF/c'+scnum+'/'
    'mp':def_path='/nfs/cluster4/unh_CDF/latest/c'+scnum+'/'
    'ppp':def_path='/nfs/cluster4/ppplus/'
    else:message, 'ftype not recognized: '+ftype
endcase
if (n_elements(path_in) eq 0) then path = def_path else path = path_in

if (n_elements(vers_in) eq 0) then vers = '*' else vers=vers_in

;==========================================================

; Define string variables for the CDF variable names, and for the filename
htype = strupcase(ftype)
hstr = 'C'+scnum+'_'+htype+'_EDI'
fstr = 'c'+scnum+'_'+ftype+'_edi'
if (ftype eq 'ppp') then fstr = strupcase(fstr)

; Define the filename
if (n_elements(filename_in) ne 0) then filename=filename_in else $
  filename=fstr+'_'+date+'*'+vers+'.cdf'

; Try to find the file
if (verbose) then message, 'Searching for: '+path+'/'+filename, /cont
f = findfile(path+'/'+filename)
if (f(0) ne '') then begin
    fname = f(n_elements(f)-1)  ; Take highest version
    goto, foundone
endif
msg_out = 'File not found: '+path+'/'+filename
if (verbose) then message, msg_out, /cont
return, 0
foundone:

if (verbose) then message, 'Opening '+fname, /cont
id = cdf_open(fname)

a = strsplit(fname,'/',/extract)
fnout = a(n_elements(a)-1)

cdf_control, id, var='Epoch__'+hstr, get_var_info=r
nrec = r.maxrec+1

if (nrec le 1) then begin
    msg_out = 'No data in '+fname
    if (verbose) then message, msg_out, /cont
    return, 0
endif

cdf_varget, id, 'Epoch__'+hstr, epoch, rec_start=0, rec_count=nrec

; Time in seconds since midnight
yr = long(strmid(date,0,4))
mo = long(strmid(date,4,2))
da = long(strmid(date,6,2))
cdf_epoch, epoch0, yr, mo, da, 0, 0, 0, 0, /compute_epoch
time = reform((epoch - epoch0)/1000.d0) ; seconds since midnight

cdf_varget, id, 'E_xyz_gse__'+hstr, evec, rec_start=0, rec_count=nrec
cdf_varget, id, 'V_ed_xyz_gse__'+hstr, vvec, rec_start=0, rec_count=nrec
cdf_varget, id, 'Status__'+hstr, sbytes, rec_start=0, rec_count=nrec
quality = reform(sbytes(0,*)); Data quality (0=bad, 1=caution, 2=good)
cdf_varget, id, 'Reduced_chi_sq__'+hstr, rchi2, rec_start=0, rec_count=nrec
rchi2=reform(rchi2)

if (ftype eq 'ppp') then begin
    cdf_varget, id, 'Nbeam__'+hstr, nbeam, rec_start=0, rec_count=nrec
    nbeam = reform(nbeam)
endif

cdf_close, id

; Apply data quality filter if desired
if (gc or n_elements(qflag) ne 0) then begin
    
    if (gc) then gd = where(quality ne 255 and quality ne 0)
    
    if (n_elements(qflag) ne 0) then begin
        case qflag of
            0: gd = where(quality ne 255) ; All valid data
            1: gd = where(quality eq 1 or quality eq 2) ; Caution and Good
            2: gd = where(quality eq 2) ; Good data
            else:begin
                msg_out = 'qflag='+string(qflag)+' not handled'
                if (verbose) then message, msg_out, /cont
                return, 0
            end
        endcase
    endif
    
    if (gd(0) ne -1) then begin
        epoch = epoch(gd)
        time = time(gd)
        evec = evec(0:2,gd)
        vvec = vvec(0:2,gd)
        sbytes = sbytes(*,gd)
        rchi2 = rchi2(gd)
        quality = quality(gd)
        if (ftype eq 'ppp') then nbeam = nbeam(gd)
    endif else begin
        msg_out = 'No data in '+fname+' after quality filter applied'
        if (verbose) then message, msg_out, /cont
        return, 0
    endelse
    
endif else begin
    gd = where(quality ne 255)
    if (gd(0) ne -1) then begin
        epoch = epoch(gd)
        time = time(gd)
        evec = evec(0:2,gd)
        vvec = vvec(0:2,gd)
        sbytes = sbytes(*,gd)
        rchi2 = rchi2(gd)
        quality = quality(gd)
        if (ftype eq 'ppp') then nbeam = nbeam(gd)
    endif else begin
        msg_out = 'No non-fill data in '+fname
        if (verbose) then message, msg_out, /cont
        return, 0
    endelse
endelse

; Apply time range filter
gd = where(time ge trange(0) and time le trange(1))

if (gd(0) ne -1) then begin
    epoch = epoch(gd)
    time = time(gd)
    evec = evec(0:2,gd)
    vvec = vvec(0:2,gd)
    sbytes = sbytes(*,gd)
    rchi2 = rchi2(gd)
    quality = quality(gd)
    if (ftype eq 'ppp') then nbeam = nbeam(gd)
endif else begin
    msg_out = 'No data within specified time range'
    if (verbose) then message, msg_out, /cont
    return, 0
endelse

; Figure out which method was used

;===============================================================
; FROM EDI_PISO:================================================
;===============================================================
; S-byte 3
;
; If Bits 0-4 aren't set, then pp_method = 0 -- Forced TRI
; Bit 0 set if pp_method = 1 -- TRI/PMT/SMT examined - Logic Chain 1
; Bit 1 set if pp_method = 2 -- Forced TOF
; Bit 2 set if pp_method = 3 -- TRI/PMT/SMT examined - Logic Chain 2
; Bit 3 set if pp_method = 4 -- TRI/PMT/SMT examined - Logic Chain 3
; Bit 4 set if pp_method = 5 -- Forced SMT
; Bit 0 and Bit 1 set if pp_method = 6 -- Forced PMT/SMT (SMT result
;                                         preferred)
; Bit 0 and Bit 2 set if pp_method = 7 -- TRI/RMT examined
; Bit 0 and Bit 3 set if pp_method = 8 -- Forced RMT
; Bit 0 and Bit 4 set if pp_method = 9 -- Both methods forced
;
; If Bits 5-6 aren't set, then method used in the end was Triangulation (TRI)
; Bit 5 set if method used in the end was Poorman's ToF (PMT)
; Bit 6 set if method used in the end was Simultan ToF (SMT)
; Bit 5 and 6 set if method used in the end was Richman's ToF (RMT)
;
; Bit 7:  Not set = No 180-degree ambiguity in drift step
;         Set     = 180-degree ambiguity exists
;    
;    sum = 0
;    if (pp_method eq 1) then sum = sum + 1
;    if (pp_method eq 2) then sum = sum + 2
;    if (pp_method eq 3) then sum = sum + 4
;    if (pp_method eq 4) then sum = sum + 8
;    if (pp_method eq 5) then sum = sum + 16
;    if (pp_method eq 6) then sum = sum + 1 + 2
;    if (pp_method eq 7) then sum = sum + 1 + 4
;    if (pp_method eq 8) then sum = sum + 1 + 8
;    if (pp_method eq 9) then sum = sum + 1 + 16
;
;    if (ep_out.method eq 1) then sum = sum + 32 ; PMT method used
;    if (ep_out.method eq 2) then sum = sum + 64 ; SMT method used
;    if (ep_out.method eq 3) then sum = sum + 32 + 64 ; RMT method used
;   
;    if (ep_out.ambig_180) then sum = sum + 128 ; 180-degree ambiguity
;
;    sout.sbyte3 = sum
;===============================================================
;===============================================================
; byte:  8-bit, unsigned ranging from 0 to 255

s3 = reform(sbytes(3,*))
nn = n_elements(s3)

; (This is klunky [but accurate] because I'm not a bit person)
; Assumption is made here that pp_method=9 (Bit 0 and 4 set), so check
; this first
w1 = where(s3 and 1)
w16 = where(s3 and 16)
if (n_elements(w1) ne nn or $
    n_elements(w16) ne nn) then message, 'S3 not expected'

; TRI can be 17 or 145
; PMT can be 49 or 177
; SMT can be 81 or 209
; RMT can be 113 or 241

itri = where(s3 eq 17 or s3 eq 145)
ipmt = where(s3 eq 49 or s3 eq 177)
ismt = where(s3 eq 81 or s3 eq 209)
irmt = where(s3 eq 113 or s3 eq 241)

;  Method Used:
;            0 = Triangulation
;            1 = Poorman's ToF
;            2 = Simultan ToF
;            3 = Richman's ToF
meth_usd = intarr(nn)+255       ; Set all values to 255
if (itri(0) ne -1) then meth_usd(itri) = 0
if (ipmt(0) ne -1) then meth_usd(ipmt) = 1
if (ismt(0) ne -1) then meth_usd(ismt) = 2
if (irmt(0) ne -1) then meth_usd(irmt) = 3

icheck = where(meth_usd eq 255)
if (icheck(0) ne -1) then message, 'Problem with S3'

return, 1
end
