; To try to read all everything in data and then process it.
;------------------ INPUT PARAMETERS ----------------------------------
;sat:        Satellite number
;            A single integer (1, 2, 3 or 4) indicating the 4 s/c
;
;specie:     0: H+, 1: He++, 2: He+, 3: O+
;
;inst:       Instrument (0: CODIF, 1: HIA)
;
;time:       Start date/time string in tplot format ('yyyy-mm-dd/hh:mm:ss')
;timespan:   Set the time span
;            (keywords: SECONDS, MINUTES, HOURS, DAYS (DEFAULT))
;
;units_name: The units should be 'Counts'
;
;eff_table:  0: Ground table, 1: OnBoard Table
;----------------------------------------------------------------------
FUNCTION sign, data
if (data eq 0) then return, 0  else return, abs(data)/data
end

;----MAIN PROGAME -------------------
;Read the data, turn the data into distribution function and then plot them as an image
; Plot in V_// and V_perp
sat   = 4
specie= 3       ;specie:     0: H+, 1: He++, 2: He+, 3: O+
time='2001-03-31/17:17:10'

timespan, time, 200, /SECONDS ; SECONDS, MINUTES, HOURS, DAYS (DEFAULT)
units_name='Counts'
inst=0; 0: CODIF, 1: HIA (this is not supported for the moment)
eff_table=0 ; 0: GROUND, 1: ONBOARD 

eph_sc = sat(n_elements(sat)-1) ; (1,2,3,4)

case specie of 
    0: prod = [12, 13]
    1: prod = [15, 16]
    2: prod = [17, 18, 46, 48]
    3: prod = [17, 18, 47, 49]
endcase
;-------------------------------------------------------------------------------
;Get the magnetic field 
;-------------------------------------------------------------------------------
plot_mag_from_crib, sat  0
tplot_names, 1, names = names
get_data, names(0), data = mgf
;----------------------------------------------------------------------------------
;Get the dat
;----------------------------------------------------------------------------------
dat = call_function('get_cis_cod_data',specie=specie,prod,sat,no_phi_cor=no_phi_cor)
packets=n_elements(dat.time)
Tpara=dblarr(packets) & Tperp= dblarr(packets)&Vtpara=dblarr(packets) & Vtperp= dblarr(packets)
index=where (dat.data eq 1)
dat.data(index)=0 
units='DIST FUNC'
dat = convert_codif_units(dat,units,'codif_ts_eff_corr', $
                              eff_table, packets=packets, $
                              old_eff=old_eff, sat=sat, $
                              incrates=incrates)
;--------------------------------------------------------------------------------
;Energy and other unchanged parameters like mass and velocity and n_elements
;------------------------------------------------------------------------------------
  na     = dat.nenergy(0) ; number of energy bins
  nb     = dat.nbins(0) ; number of angle bins
  energy = dat.energy(*,*,0)
  denergy = dat.denergy(*,*,0)
  theta = dat.theta/!radeg
  phi = dat.phi/!radeg
  dtheta = dat.dtheta/!radeg
  dphi = dat.dphi/!radeg
  mass = dat.mass
  Const = (mass/(2.*1.6e-12))^(.5)
  domega=2.*dphi*cos(theta)*sin(.5*dtheta)
;-----------------------------------------------------------
;All the velocity in the spacecraft frame they are also unchanged for the points we consider
;---------------------------------------------------------------------------- 
   vel_x=dblarr(packets) 
   vel_y=dblarr(packets)
   vel_z=dblarr(packets) 
   vx=cos(theta)*cos(phi)*sqrt(2*energy/mass)
   vy=cos(theta)*sin(phi)*sqrt(2*energy/mass)
   vz=sin(theta)*sqrt(2*energy/mass)
;-----------------------------------------
;Cordinate transform for inst_coord to GSE
;-------------------------------------------
  vel=dblarr(3,1) 
  vxg=dblarr(na, nb)       ; vx in GSE coodinate
  vyg=dblarr(na, nb)       ; vy in GSE coodinate
  vzg=dblarr(na, nb)       ; vz in GSE coodinate
  mid_times = dat.time+(dat.end_time-dat.time)/2d
  inst=0
 
  FOR j = 0, na-1 DO BEGIN
      FOR k=0, nb-1 do begin     
        vel=[[vx(j,k)],[vy(j,k)],[vz(j,k)]]
        datastr = {x:mid_times(0),y:vel}
        store_data, 'dd', dat=datastr, $
        dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
            phase_instr:dat.phase_instr}
        IF NOT KEYWORD_SET(INST_COORD) THEN BEGIN
             cis_coord_trans, DATA_IN = 'dd', TRANS='CODIF->GSE', $
          DATA_OUT = 'dd_new'
          get_data, 'dd_new', data=d, dlim=dlim
          vxg(j,k)=d.y(0)
         vyg(j,k)=d.y(1)
          vzg(j,k)=d.y(2)
          store_data, 'dd_new', /DELETE
        ENDIF
        store_data, 'dd', /DELETE
      ENDFOR 
   ENDFOR
;-------------------------------------------------------------------------------
;Do a contour plot on vxg, vyg, vzg and Disf
;-------------------------------------------------------------------------------
  vxo=dblarr(na*nb)       ; vx in gse coodinate in ONE dimension
  vyo=dblarr(na*nb)       ; vy in gse coodinate in ONE dimension
  vzo=dblarr(na*nb)       ; vz in gse coodinate in ONE dimension 
  vxo=vxg(*) & vyo=vyg(*) & vzo=vzg(*) 
;--------------------------------
; Here we us the HIA moments of velocity 
 vxo=vxo+200 & xyo=vyo-0 & vzo=vzo-100

disped=create_struct('time',0.,'v',dblarr(51),'Fx',dblarr(51), 'Fy',dblarr(51), 'Surf', dblarr(51,51))
disped_s=replicate(disped, packets)
;stop
FOR ii=0, packets-1  DO BEGIN
  mfield=mgf.y(ii*2,*)
 bx=mfield(0, 0) & by=mfield(0,1) & bz=mfield(0,2) 
  e0=dblarr(1,3) &e1=dblarr(1,3)& e2=dblarr(1,3)                     ;Direction units
  e0=[bx,by,bz]/sqrt(bx^2+by^2+bz^2)
  e1=crossp(e0, [1,0,0])
  e2=crossp(e0, e1)
  e1=e1/sqrt(e1(0)^2+e1(1)^2+e1(2)^2)
  e2=e2/sqrt(e2(0)^2+e2(1)^2+e2(2)^2)
  vb0=dblarr(na*nb)& vb1=dblarr(na*nb)& vb2=dblarr(na*nb)& vbp=dblarr(na*nb) 
    ; vbp Perpendicular speed
  ;-----------------------------------------------------------------
  ;The direction of perpedicular and parallet to the magnetic field 
  ;---------------------------------------------------------------
  for i=0, na*nb-1 do begin
     vb0(i)=total([vxo(i),vyo(i),vzo(i)]*e0, /NaN)
     vb1(i)=total([vxo(i),vyo(i),vzo(i)]*e1, /NaN)
     vb2(i)=total([vxo(i),vyo(i),vzo(i)]*e2, /NaN)
     vbp(i)=sqrt(vb1(i)^2+vb2(i)^2)*sign(vb1(i))
   end 
  triangulate, vb0, vbp, tri
 ;------------------------------------------------------------------------------------
 ;Set the coutour parameters to make the plot more beautiful
 ;---------------------------------------------------------------------------
 range=[1.e-15, 1.e-6]
 nlines=20
 maximum=range(1)
 minimum=range(0)
 resolution=51 & xrange=[-2000, 2000]& spacing = (xrange(1)-xrange(0))/(resolution-1)
 data=dblarr(na*nb)
 data1 = dat.data(*,*,ii) & data=data1(*)  ; Read the data form the orginal dat and
thesurf = trigrid(vb0,vbp,data,tri,[spacing,spacing], [xrange(0),xrange(0),xrange(1),xrange(1)],xgrid =xg,ygrid = yg)
  thesurf=smooth(thesurf, 3)
  F_vx=dblarr(resolution) & F_vy=dblarr(resolution)
  FOR i=0, resolution-1 do  F_vx(i)=total(thesurf(i, *)*abs(yg))*80
  FOR i=0, resolution-1 do  F_vy(i)=total(thesurf(*, i))*80
  disped.v=xg &  disped.Fx=F_vx &  disped.Fy=F_vy & disped.surf=thesurf
  disped.time=dat.time(ii)-11412.*3600*24-17.*3600-17*60
  disped_s(ii)=disped
   Vtpara(ii)=sqrt(total(F_vx*xg^2)/total(F_vx))
   Vtperp(ii)=sqrt(total(F_vy*abs(yg)^3)/total(F_vx)/2)
 ; PLot the distribution function
 if  (ii eq 0) then plot, yg, F_vy/max(F_vy), /ylog,xtitle='VPerp', ytitle='Normalized Distribution function'
 ;if  (ii eq 1) then oplot, yg, F_vy/max(F_vy), linestyle=1
 if  (ii eq 1) then oplot, xg, F_vx/max(F_vx), linestyle=1
 if  (ii eq 4) then oplot, xg, F_vy/max(F_vy), linestyle=4
  ;if  (ii eq 1) then plot, xg, F_vx/max(F_vx), /ylog, xtitle='VPara', ytitle='Normalized Distribution function'
; if  (ii eq 1) then oplot, xg, F_vy/max(F_vy), linestyle=1
 ;if  (ii eq 2) then oplot, xg, F_vx/max(F_vx), linestyle=2
 ;if  (ii eq 7) then oplot, xg, F_vx/max(F_vx), linestyle=7
 end 
 anis=(Vtperp/Vtpara)^2
 Tperp=Vtperp^2*mass
 Tpara=Vtpara^2*mass
 CASE specie of 
    0: begin
       store_data, 'Anis_H', data={x:dat.time, y:anis} 
       store_data, 'Tpara_H', data={x:dat.time, y:Tpara}
       store_data, 'Tperp_H', data={x:dat.time, y:Tperp}
       store_data, 'Total_H', data={x:dat.time, y:Tperp*2+Tpara}
       end
     1: begin
       store_data, 'Anis_He', data={x:dat.time, y:anis} 
       store_data, 'Tpara_He', data={x:dat.time, y:Tpara}
       store_data, 'Tperp_He', data={x:dat.time, y:Tperp}
       store_data, 'Total_He', data={x:dat.time, y:Tperp*2+Tpara}
       end
     endcase 
;openw, 1, '/home/yliu/CCAT/ccat_user/T_and_anis.txt'
;printf, 1, 'Time      ', 'Tperp          ', 'Tpara     ', 'Anis '
;FOR ii=0, packets-1  DO  printf, 1, time_string(dat.time(ii)), Tperp(ii), Tpara(ii), anis(ii)
;close, 1
END





