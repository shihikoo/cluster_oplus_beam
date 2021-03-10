pro draw_b_field

if keyword_set(ps) then popen,'b_field.ps',/land
plot,[1,1],[2,2],xrange=[-2,2],yrange=[-2,2],xstyle=16,ystyle=16,/nodata
t=[10-indgen(1000)*0.01,indgen(1000)*0.01]

x=t
y=sqrt(1-(x-1)^2)
oplot,x,y
oplot,-x,-y
oplot,x,-y
oplot,-x,y



stop
end
