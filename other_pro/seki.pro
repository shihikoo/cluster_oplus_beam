pro seki
pi=3.1415926
m=16*1.6e-27/(1.6e-19) ; unit: keV/(m/s)^2
;v = 18.e3  ; - 154eV
;print, m*v^2/2

nt=(1.e-3)*1e6 ;m-3
B0= 33695.9 ;nT

Bt=43.;nT
Vparat = 44.e3 ;m/s 
dvt = 9.e3 ;m/s
dvu = 9.e3 ;m/s
Tparat = 1.*(1.6e-19)  ;eV

v0 = (indgen(2e3)+1)*1e-3*40.e3 ;m/s
vth0 = (indgen(1e2)+1)*1e-2*260 ;m/s ;0.08e3

a=(vparat-dvt-V0)#(1/Vth0)
b=(vparat+dvu-V0)#(1/Vth0)

W1=erf(b)-erf(a)
;W2=-2*Vth0#(a*exp(-a^2)-b*(exp(-b^2)))/((4/m)*Tparat+2*(Vparat-V0)^2-Vth0^2)
!p.multi=[0,2,2]
specplot,v0*1e-3,(vth0^2*m)/2*1e3,erf(a),no_interp=1,lim={xlog:1,ylog:1,zlog:1,ztitle:'a',no_erase:1,xstyle:1,ystyle:1}

specplot,v0*1e-3,(vth0^2*m)/2*1e3,erf(b),no_interp=1,lim={xlog:1,ylog:1,zlog:1,ztitle:'b',no_erase:1,xstyle:1,ystyle:1}
specplot,v0*1e-3,(vth0^2*m)/2*1e3,erf(b),no_interp=1,lim={xlog:1,ylog:1,zlog:0,ztitle:'W1',no_erase:1,xstyle:1,ystyle:1}
;specplot,v0,vth0,w2,no_interp=1,lim={xlog:1,ylog:1,zlog:0,ztitle:'W2',no_erase:1}
!p.multi=[0,0,0]
;print,(Vth0(exp(-b^2)-exp(-a^2))/sqrt(pi)/(vparat-V0))
;, (exp(-b^2)-exp(-a^2))/(Vparat-V0)
;print,-2*Vth0*(a*exp(-a^2)-b(exp(-b^2)))/((4/m)*Tparat+2*(Vparat-V0)^2-Vth0^2)

x=V0/Vth0
D=exp(-x^2)+sqrt(pi)*x*(1+erf(x))

S0=Vth0^2*nt*D/((sqrt(pi)*(Bt/B0)*Vth0*(erf(b)-erf(a))))
stop
end
