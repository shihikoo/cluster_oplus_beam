pro plot_flux_year,plot_ps=plot_ps
path='median/'
spawn,'mkdir '+path
spawn,'mkdir '+path+'log'
spawn,'mkdir '+path+'linear'

tplot_restore,filenames='median.tplot'
tplot_names
;property_name_set=['nVpara_over_B','flux','density','nV']
property_name_set=['eflux'];,'nV']
region_name_set=['outflow','lobe']
storm_phase_set=['nonstorm','storm']
energy_threshold_set=[[0,100],[100,1000],[1000,40000]]
amme=fltarr(n_elements(energy_threshold_set(0,*)),7, n_elements(storm_phase_set), n_elements(region_name_set), n_elements(property_name_set), 3)
energy_threshold_str = strarr(2,3)
energy_threshold_strr = strarr(2,3)
for ien =0, n_elements(energy_threshold_set(0,*))-1 do begin 
    energy_threshold_str(*,ien)=strcompress(STRING(energy_threshold_set(*,ien), format = '(i5.1)'),/remove_all)
 energy_threshold_strr(*,ien)=strcompress(STRING(energy_threshold_set(*,ien), format = '(i5.5)'),/remove_all)
    var='mme_'+'energy_gt'+energy_threshold_strr(0,ien)+'_lt'+energy_threshold_strr(1,ien)+'_eflux_ge_1100'
    get_data, var ,data=data
    amme(ien,*,*,*,*,*)=data.mme
endfor 


for ip=0,n_elements(property_name_set)-1 do begin 
    for is=0,n_elements(storm_phase_set)-1 do begin
        for ir=0,n_elements(region_name_set)-1 do begin 
            CASE property_name_set(ip) OF 
                'flux': property_range = [10, 2600.] 
                'density': property_range = [0.001, 0.07]
                'nV': property_range = [0.1, 3]
                'nVpara_over_B': property_range=[0.0005,0.05]
                'eflux': property_range = [1e3, 1e5] 
                ELSE: stop
            ENDCASE  

            if keyword_set(plot_ps) then popen,path+'log/'+'lineplots_'+property_name_set(ip)+'_'+storm_phase_set(is)+'_'+region_name_set(ir)+'.ps' ,/land else window,ir
    
            plot,[2000.5,2007.5],[0,1],/nodata,xstyle=1,$
              title = storm_phase_set(is)+',  '+region_name_set(ir),xtitle='year', ytitle='median  ' $
              +property_name_set(ip),yrange=property_range,charsize=1.5,   position=[0.15,0.15,0.95,0.85],ystyle=1,ylog=1
            oplot, [2001,2002,2003,2004,2005,2006,2007], amme(0,*,is,ir,ip,1), psym=-1, color=2, thick=10
            oplot, [2001,2002,2003,2004,2005,2006,2007], amme(1,*,is,ir,ip,1), psym=-1, color=1, thick=10
            oplot, [2001,2002,2003,2004,2005,2006,2007], amme(2,*,is,ir,ip,1), psym=-1, color=6, thick=10
            xyouts, 2000.6, property_range(1)/(10^1.1), energy_threshold_str(0,0)+'eV to '+energy_threshold_str(1,0)+'eV', color=2, charsize=3, charthick = 3
            xyouts, 2000.6, property_range(1)/(10^1.2), energy_threshold_str(0,1)+'eV to '+energy_threshold_str(1,1)+'eV', color=1, charsize=3, charthick = 3
            xyouts, 2000.6, property_range(1)/(10^1.3), energy_threshold_str(0,2)+'eV to '+energy_threshold_str(1,2)+'eV', color=6, charsize=3, charthick = 3
            if keyword_set(plot_ps) then pclose else stop
        endfor 
    endfor 
endfor 

spawn,'mogrify -format png median/log/*.ps'
spawn,'mogrify -rotate -90 median/log/*.png'
;    if keyword_set(plot_ps)  then popen, plot_path+'median_line_plots/'+sort_condi+'_'+property_names(ip)$
;      +'_'+single_phase_set(is)+'_mean.ps' ,/land else window,2
;    tplot, [single_phase_set(is)+'_'+region_name_set(0)+'_'+property_names(ip)+'_mean', single_phase_set(is) $
;            +'_'+region_name_set(1)+'_'+property_names(ip)+'_mean']
;    if keyword_set(plot_ps) then pclose else stop
end
