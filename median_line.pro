pro median_line,x,y,xrange,result,xlog=xlog,xgrid=xgrid
if not keyword_set(xgrid) then xgrid =1
xlog=1
if keyword_set(xlog) then begin 
    n=alog10(CEIL(max(xrange)/min(xrange)))/xgrid
    ymedian=DBLARR(n)
    xmiddle=DBLARR(n)
    for i = 0, n-1 do begin
        index = where(x gt 10^(alog10(min(xrange))+i*xgrid) and x lt 10^(alog10(min(xrange))+(i+1)*xgrid),ct)
        if ct gt 0 then ymedian(i) = median(y(index))
;        if ct gt 0 then ymedian(i) = mean(y(index),/nan) 
        xmiddle(i)= 10^(alog10(min(xrange))+(i+0.5)*xgrid)
    endfor
    result=[[xmiddle],[ymedian]]
endif 
end
