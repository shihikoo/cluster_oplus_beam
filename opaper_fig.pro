PRO opaper_fig, figure=figure,ps_plot=ps_plot

if keyword_set(figure)  then n_fig=n_elements(figure) else figure=1
output_path = '~/o_beam/opaper/'
spawn,  'mkdir '+output_path
spawn, 'mkdir '+output_path+'tplot_restore/'
spawn, 'mkdir '+output_path+'plots/'
;------------------------------------------------------------
IF figure eq 1 THEN BEGIN 
    save_data = 0
    en_recalc = 0
    iden_recalc=0

    sc = 4  &  sc_str = STRING(sc, FORMAT = '(i1.1)')
    average_time = 5* 60 &    at_str = STRCOMPRESS(average_time, /REMOVE_ALL)   
    time_start = time_double('2002-09-10/13:00:00')
    time_end =time_double('2002-09-11/13:00:00')
    dt= time_end-time_start
    timespan, time_start, dt, /SECONDS
; check the restored data           
    flndata = output_path+'/tplot_restore/spectra'
    print, FINDFILE(flndata+'.tplot', COUNT = ct)   
    IF ct GT 0 AND (iden_recalc EQ 0 or en_recalc EQ 0) THEN tplot_restore, filenames = flndata+'.tplot' 

    IF KEYWORD_SET(en_recalc) OR ct EQ 0 THEN BEGIN 
; Load CLUSTER H+ and O+ energy spectra  
        sat = [sc, sc] & specie = [0, 3] & angle = [[-90, 90], [0, 360]]
        inst = 0 & units_name = 'DIFF FLUX' & eff_table = 0     
        plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 1
; Load CLUSTER O+ pitch angle spectra for low energy
        sat = [sc] &  specie=[3] & energy=[40,40000.]
        inst = 0 & eff_table = 0  &  units_name = 'DIFF FLUX'
        plot_pa_spec_from_crib, sat, specie, inst, units_name,  energy, eff_table,  recalc = 1, COMBINE = 1
;stop
;load Dst Index
        read_omni
    endif

    IF KEYWORD_SET(iden_recalc) OR ct EQ 0 THEN BEGIN 
        beam_recalc = 1
        find_phase = 1
        add_imf = 1
        display_time = 24.*60*60
        find_o_beam, sc = sc,  average_time = average_time,   path = output_path,   beam_recalc = beam_recalc,  find_phase = find_phase,   displaytime = display_time,   add_imf = add_imf
    endif 

    p1 = 'ACE_SWEPAM_VELOCITY'
    p2 = 'ACE_SWEPAM_PRESSURE'
    p3 = 'ACE_MAG_GSM_Z'       
    p7 = 'Dst_Index'
  ;  stop
    p4 = 'ENSPEC_SC4_IN0_PHI0_360_UNDIFFFLUX_SP0_ET0_All'
    p5 = 'ENSPEC_SC4_IN0_PHI0_360_UNDIFFFLUX_SP3_ET0_All'
    p6 = 'PASPEC_EN00040_40000_SC4_UNDIFFFLUX_SP3_All'

    p10='ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300'
;    p20='ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300'
    p11='PASPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300'
 ;   p21='PASPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300'
    p12='PASPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP'
  ;  p22='PASPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP'
    p13='PASPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP_ET_beam'
 ;   p23='PASPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP_ET_beam'
    p14='ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut'
;    p24='ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut'
    p15='ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut_beam'
;    p25='ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut_beam'

    options,p14,'thick',6
    ylim,p4,30,40000.,1
    ylim, p1, 350, 450, 0
    ylim, p2, 0.1, 1e1, 1
    ylim, p3, -20, 20, 0
    ylim,p5,40.,40000.
    ylim,p10,40.,40000.
;    ylim,p7,-80,-20
    zlim,[p6,p13],0,100,0
    zlim, p5, 0.1, 1e2, 1           
    options, '*', 'panel_size', 1
    options,[p4,p5,p6,p10,p11,p13],'xticklen',0.15
    options,[p4,p5,p6,p10,p11,p13],'yticklen',0.02
    options,[p4,p5,p6,p10,p11,p13],'zticklen',0
    options,[p6,p11,p13],'zticks',2
    options,[p4,p5,p10],'zticks',3
    
    var_label = 'EPH_SC' + sc_str + '_'
    var_label = var_label +    ['MLT', 'GSE_X', 'GSE_Y', 'GSE_Z', 'DIST']

    options, p5, 'ztitle', ''   
    options,p6,'ztitle',''
    options, p4, 'ytitle', 'SC'+sc_str+'!C!CH!U+!N (eV)'
    options, p5, 'ytitle', 'SC'+sc_str+'!C!CO!U+!N (eV)'
    options,p6,'ytitle','SC'+sc_str+' O!U+!N!C!CPitch Angle!C!C40-40k (eV)'
    options,p10,'ytitle','Taiward!C!CO!U+!N (eV)!C!CAvg-300s'
    options, p11,'ytitle','O!U+!N!C!CPitch Angle!C!CAt Energy Peak'
    options,p13,'ytitle','O!U+!N Beam!C!CIdentification!C!CE------T'

    IF KEYWORD_SET(ps_plot) THEN   popen,  output_path+'plots/spectra', /land
    time_stamp,/off
    tplot,[p4,p5,p6,p10,p11,p13], var_label = var_label

    tplot_panel,v=p10,o=p14
;    tplot_panel,v=p20,o=p24
;    tplot_panel,v=p10,o=p15,psym=7
;    tplot_panel,v=p20,o=p25,psym=7
;   yline, [p3, p10], col = 1
;    yline, p10, offset = -50, col = 1
    if keyword_set(ps_plot) then pclose 
    IF KEYWORD_SET(save_data) THEN  tplot_save, filename = flndata   
stop
endif  

IF figure eq 100 THEN BEGIN 
en_flux,sc=4,ps=plot_ps,recalc=0,flux_plot=1,output=output_path+'plots/en_flux'
stop
endif 

END
