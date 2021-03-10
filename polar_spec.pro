; r is the raduis, ang is the angle 
; units_angle : 0 - rad, 1 - degree, 2 - clock
pro polar_spec,r,angle, z,$
               r_range=r_range,angle_range=angle_range, zrange=zrange, $
               zlog=zlog, zticks=zticks,zticklen=zticklen,units_angle=units_angle,$
               ztitle=ztitle,title=title,xtitle=xtitle,ytitle=ytitle,$
               charsize=charsize
IF n_elements(SIZE(z,/dim)) ne 2 then begin 
    print,'z has to be an array of 2'
    stop
    goto, out
endif else begin 
    if n_elements(r) ne n_elements(z(*,0)) or n_elements(angle) ne n_elements(z(0,*)) then begin 
        print, 'size of the array z has to be the same of # of r by #  of angle'
        stop
        goto, out
    endif 
endelse 

if not keyword_set(r_range) then r_range=[min(r),max(r)]*1.1
if not keyword_set(angle_range) then angle_range=[min(angle),max(angle)]*1.1
if not keyword_set(zrange) then zrange=[min(z),max(z)]
if not keyword_set(zrange) then zticks=0
if not keyword_set(zticklen) then zticklen=0

if keyword_set(units_angle) then begin 
    if units_angle eq 1 then begin 
        angle=angle/180.*!pi
        grid_angle= grid_angle/180.*!pi
        angle_range= angle_range/180.*!pi
    endif 
    if untis_angle eq 2 then begin 
        angle=angle/24.*360./180.*!pi
        grid_angle= grid_angle/24.*360./180.*!pi
        angle_range= angle_range/24.*360./180.*!pi
    endif 
endif else print,'angle is in unit: rad'
nr = n_elements(r)
nangle=n_elements(angle)
index=sort(r)
r=r(index)
for iangle=0,nangle-1 do  z(*,iangle)=z(index,iangle)
index=sort(angle)
angle=angle(index)
for ir=0,nr-1 do z(ir,*)=z(ir,index)

plot,[0,0],[-30,30], /nodata, $
  xrange=max(ABS(r_range))*[-1,1],yrange=max(ABS(r_range))*[-1,1], $
  xstyle=5,ystyle=5, charsize=charsize,$
  position=[0.1,0.1,0.85,0.85],$
  title=title, xtitle=xtitle,ytitle=ytitle,ztitle=ztitle
for ir = 0, nr-1 do begin 
    if ir eq 0 then r0 = r(0)-(r(1)-r(0))/2.  else r0=(r(ir)+r(ir-1))/2.
    if ir eq nr-1 then r1 = r(ir)+(r(ir)-r(ir-1))/2. else    r1=(r(ir)+r(ir+1))/2.
    for iangle = 0, nangle-1 do begin 
        if iangle eq 0 then  angle0=angle(0)-(angle(1)-angle(0))/2. else angle0=(angle(iangle)+angle(iangle-1))/2. 
        if iangle eq nangle-1 then angle1= angle(iangle)+(angle(iangle)-angle(iangle-1))/2. $
        else  angle1=(angle(iangle)+angle(iangle+1))/2.

        pixpos=fltarr(2,4)
        pixpos(*,0)= [r0*cos(angle0),r0*sin(angle0)]
        pixpos(*,1)= [r1*cos(angle0),r1*sin(angle0)]
        pixpos(*,2)= [r1*cos(angle1),r1*sin(angle1)]
        pixpos(*,3)= [r0*cos(angle1),r0*sin(angle1)]
        
        x=fltarr(100,4) & y=fltarr(100,4)
        x(*,0)=(indgen(100)*0.01)*(pixpos(0,1)-pixpos(0,0))+pixpos(0,0)
        y(*,0)=(indgen(100)*0.01)*(pixpos(1,1)-pixpos(1,0))+pixpos(1,0)
        x(*,1)=r1*cos((indgen(100)*0.01)*(angle1-angle0)+angle0)
        y(*,1)=r1*sin((indgen(100)*0.01)*(angle1-angle0)+angle0)
        x(*,2)=(indgen(100)*0.01)*(pixpos(0,3)-pixpos(0,2))+pixpos(0,2)
        y(*,2)=(indgen(100)*0.01)*(pixpos(1,3)-pixpos(1,2))+pixpos(1,2)
        x(*,3)=r0*cos((indgen(100)*0.01)*(angle1-angle0)+angle0)
        y(*,3)=r0*sin((indgen(100)*0.01)*(angle1-angle0)+angle0)
    ;    oplot,x,y,psym=3
     ;   oplot,pixpos(0,*),pixpos(1,*),color=1,psym=1

        z1 = z
        zrange1=zrange
        if keyword_set(zlog) then begin
            bad = where(z1 le 0,ct)
            if ct ne 0 then z1(bad) = !values.f_nan
            z1 = alog10(z1)
            zrange1 = alog10(zrange1)
        endif
        color_z= ((1.*z1(ir,iangle)-min(zrange1))/(max(zrange1)-min(zrange1)))*(254-7)+7 
        if color_z lt 7 then color_z = 7
        if color_z gt 254 then color_z =254
        if finite(color_z) then polyfill,x,y,color= color_z
    endfor 
endfor 
xyouts,max(r_range),0,'0'
xyouts,-(max(r_range)*1.05),0,'12'
xyouts,0,(max(r_range)),'6'
xyouts,0,-(max(r_range*1.05)),'18'
if keyword_set(zlog) then zrange1=10.^zrange1
draw_color_scale,range=zrange1,log=zlog,title=ztitle,charsize=charsize,yticks=zticks,ticklen=zticklen
out:
end 
