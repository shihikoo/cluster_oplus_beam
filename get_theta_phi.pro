;+ 
; PROCEDURE:  get_theta_phi
;
; PURPOSE: get the average or all magnetic theta and phi data in 
;          codif coordinates
;
; INPUT:
;        sat -> spacesraft number
;        dat -> instrument data structure
;        inst -> instrument number 0: CODIF, 1:HIA
;        
; OUTPUT: 
;         mag_theta -> magnetic field theta in instrument coordinates
;         mag_phi   -> magnetic field phi in instrument coordinates
;
; CREATED BY: C. Mouikis
;
; LAST MODIFICATION: 05/22/01
;
; MODIFICATION HISTORY:
;     04/11/01 - New keyword ALL to allows to output all the magnetic
;                field data and not the average value only.
;     04/30/01 - The keyword SPIN_AVG is taken out of the call to the
;                CNES data
;     05/22/01 - CNES magnetic field data are the default
;     06/20/02 - input variable inst introduced
;-
PRO get_theta_phi, sat, dat, mag_theta, mag_phi, inst, $
                   IC=IC, $
                   CNES=CNES, $
                   FORNACON=FORNACON, $
                   EDITA=EDITA, $
                   ALL = ALL
  
  IF KEYWORD_SET(FORNACON) THEN BEGIN
    get_cluster_mag_fornacon, sat, Btime, Bxyz ; get magnetic field data
  ENDIF ELSE BEGIN
    IF KEYWORD_SET(EDITA) THEN BEGIN
      get_cluster_mag_edita, sat, Btime, Bxyz 
    ENDIF ELSE BEGIN
      get_cluster_mag_gse, sat, Btime, Bxyz, $
        IC=IC, CNES=CNES
    ENDELSE
  ENDELSE
;jing: some Bxyz eq -1e31, which seems to be fause data, so I set them all to nan data. Those data influence dramtically on the average magnetic field but cause less problem to non-average data. 
if total(where(Bxyz eq -1e31)) gt 0 then Bxyz(where(Bxyz eq -1e31))=!values.f_nan  

  IF NOT(KEYWORD_SET(ALL)) THEN BEGIN
    IF Btime(0) NE -9999.9 THEN BEGIN
      ; calculate average value
      Btavg = (Btime(0)+Btime(n_elements(Btime)-1))/2.
      Bxavg = mean(Bxyz(*,0), /NaN)
      Byavg = mean(Bxyz(*,1), /NaN)
      Bzavg = mean(Bxyz(*,2), /NaN)
      
      ; store data in tplot structure
      Bdata_gse = {x:Btavg, y:[[Bxavg],[Byavg],[Bzavg]]}
    
      store_data, 'B_xyz_gse', dat = Bdata_gse, $
        dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
              phase_instr:dat.phase_instr}
    
      ; transform data from gse to codif coordinates
      cis_coord_trans, DATA_IN = 'B_xyz_gse', $
        TRANS = 'GSE->CODIF', $
        OUT_FILE = 'No', DATA_OUT = 'B_xyz_codif'
      
      ; calculate mag_theta & mag_phi
      get_data, 'B_xyz_codif', data = Bdata_cod
      
      ravg = SQRT(Bdata_cod.y(0,0)^2+Bdata_cod.y(0,1)^2)
      mag_theta = (ATAN(Bdata_cod.y(0,2)/ravg))*!RADEG 
      mag_phi = (ATAN(Bdata_cod.y(0,1),Bdata_cod.y(0,0)))*!RADEG
      mag_phi = mag_phi - 360*FLOOR(mag_phi/360.)
    ENDIF
  
    IF Btime(0) EQ -9999.9 THEN BEGIN
      print, 'NO CORRESPONDING MAGNETIC FIELD FILES FOUND.'
      print, 'MAGNETIC FIELD PHI & THETA VALUES CAN BE PROVIDED '
      print, 'MANUALLY USING THE CORRESPONDING KEYWORDS'
    ENDIF
    
  ENDIF ELSE BEGIN

    ; interpolate B values on codif data times
    Bx = INTERPOL(Bxyz(*,0), Btime, (dat.time + dat.end_time)/2.)
    By = INTERPOL(Bxyz(*,1), Btime, (dat.time + dat.end_time)/2.)
    Bz = INTERPOL(Bxyz(*,2), Btime, (dat.time + dat.end_time)/2.)
    
    ; store data in tplot structure
    Bdata_gse = {x:(dat.time + dat.end_time)/2., y:[[Bx], [By], [Bz]]}
    
    store_data, 'B_xyz_gse', dat = Bdata_gse, $
      dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
            phase_instr:dat.phase_instr}

    ; transform data from gse to codif coordinates
    cis_coord_trans, DATA_IN = 'B_xyz_gse', $
      TRANS = 'GSE->CODIF', $
      OUT_FILE = 'No', DATA_OUT = 'B_xyz_codif'
    
    ; calculate mag_theta & mag_phi
    get_data, 'B_xyz_codif', data = Bdata_cod
    
    ravg = SQRT(Bdata_cod.y(*,0)^2+Bdata_cod.y(*,1)^2)
    mag_theta = (ATAN(Bdata_cod.y(*,2)/ravg))*!RADEG 
    mag_phi = (ATAN(Bdata_cod.y(*,1),Bdata_cod.y(*,0)))*!RADEG
    mag_phi = mag_phi - 360*FLOOR(mag_phi/360.)
    
    store_data, 'B_xyz_gse', /DELETE
    store_data, 'B_xyz_codif', /DELETE
  ENDELSE
    
END
