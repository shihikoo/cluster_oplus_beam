PRO pastudy
alfa10 = 180-145 ;degree
alfa20 = 180-172 ;degree
energy1 = 530 ;eV
energy2 = 330 ;eV
pi = 3.1415926

B0 = 41 ;nT
x0 = 17.5
y0 = 7.6
z0 = 0.9
r0 = 19.1
;----------------------------------

namda0 = ATAN(SQRT(x0^2+y0^2)/z0)

cr = r0/COS(namda0)^2
c=COS(namda0)^6/(1+3*SIN(namda0)^2)^0.5

namda = (85+findgen(100000)*0.00005)*3.14/180

r = cr*COS(namda)^2

s = sqrt(1+3*sin(namda)^2)/cos(namda)^6

B = B0*c*s
;-----------------------------
;!p.multi = [0, 1, 2, 0, 0]
;window, 1
f1 = sin(alfa10*pi/180)^2*c*s
index1 = where((f1-1) < 0, ct1)

PLOT, [0, 90], [90, 90],$
      title = 'Pitch Angle  & E!D//!N    vs    Latitude!C', $
      xtitle = 'Latitude ( degree )', $
      ytitle = 'Pitch Angle ( degee ) -- Blue Lines',  $
      xrange = [89, 86], yrange = [0, 90], $
      XSTYLE = 1, ystyle = 1,charsize = 1.2,$
      position = [0.1, 0.1, 0.9, 0.9], /nodata

;AXIS, XAXIS=1, XTITLE='Distance    ( R!D E !N )',xstyle = 1,$
 ;     charsize = 1,xrange = cr*cos(!X.CRANGE*pi/180)^2

alfa1 = ASIN(SQRT(f1(0:index1(ct1-1))))
oplot, namda(0:index1(ct1-1))*180/pi, alfa1*180/pi, color = 3

namda_max1 = namda(index1(ct1-1))*180/pi
oplot,[namda_max1,namda_max1], [0, 90] , color = 3
;-------
f2 = sin(alfa20*pi/180)^2*c*s
index2 = where((f2-1) < 0, ct2)

alfa2 = ASIN(SQRT(f2(0:index2(ct2-1))))
oplot, namda(0:index2(ct2-1))*180/pi, alfa2*180/pi, color = 2

namda_max2 = namda(index2(ct2-1))*180/pi
oplot, [namda_max2, namda_max2],[0, 90], color = 2

oplot, [namda0*180/pi,namda0*180/pi], [0, 90], color = 4

xyouts,87.25,10,'Beam2 ( E!DTotal!N:530eV, PA:35!Uo!N )',charsize=1.2
xyouts,87.25,35,'Beam1 ( E!DTotal!N:330eV, PA:8!Uo!N )',charsize=1.2
xyouts, 87.5, 85, 'Green Line: Cluster Data point!C!C B = 41 nT, r = 19.1!C!C Latitude = 87.3!Uo!N',$
        charsize = 1.2
xyouts, 88.7, 80, $
        'r ='+STRING(cr*cos(namda_max2*pi/180)^2, format = '(i3.1)')+'!C!C'+$
        'B!DR!N ='+ $
        STRING(B(where(namda EQ namda_max2*pi/180)), format = '(i4.1)'), $
                    charsize = 1.2
xyouts, 87.9, 38,$
        'r ='+STRING(cr*cos(namda_max1*pi/180)^2,format = '(i3.1)')+'!C!C'+$
        'B!DR!N ='+ $
        STRING(B(where(namda EQ namda_max1*pi/180)), format = '(i4.1)'), $
                    charsize = 1.2
;--------------------------------
ell10 = energy1*cos(alfa10*pi/180)^2
ell20 = energy2*cos(alfa20*pi/180)^2 
print,  ell10, ell20

ell1 = energy1*cos(alfa1)^2
ell2 = energy2*cos(alfa2)^2

plot,  [1, 1], [2, 2], ystyle =1+4,xstyle = 4, $
       xrange = [89, 86], yrange = [0, 600],$
       position = [0.1, 0.1, 0.9, 0.9], /nodata, /noerase
oplot,  namda(0:index2(ct2-1))*180/pi, ell2, color = 1
oplot,  namda(0:index1(ct1-1))*180/pi, ell1, color = 6

axis, yaxis = 1, yrange = [0, 600], yticks = 6, $
      ytitle = 'E!D//!N (eV) -- Red Lines', charsize = 1.2

xyouts, 87.1, 400, 'E!D//!N of Beam1: 355.635 eV !C!C!C!CE!D//!Nof Beam2: 323.608 eV '


stop

END 
