PRO plot_plasma_beta
sc = 4
sc_str = STRING(sc, format = '(i1.1)')

found = 0

FOR it = 0, 433 DO BEGIN 
    timespan, time_double('2004-10-25/00:00:00')+it*(24.*60.*60.), 1, /days
    tplot_names, names = names
    store_data, delete = names  
    IF it LE  67 $ 
      OR (it GE 219 AND it LE 349) $
      OR (it GE 351 AND it LE 401) THEN BEGIN 

        sat = sc

        sat = [sc] 
        get_cluster_ephemeris, sat, /GSE_X, /GSE_Y, /GSE_Z

        plot_mag_from_crib, sat

;-- Load CLUSTER H+ and O+ moments--
        sat = [sc, sc]
        specie = [0, 3]
        moments = ['A', 'A']
        angle = [[-90, 90], [0, 360]]
        energy = [40., 40000.]
        inst = 0 &  eff_table = 0

        plot_3dmom_from_crib, sat, specie, inst, $
          moments, angle, energy, eff_table, recalc = 0
        
        h_press =  'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP0_ET0_All'
        o_press = 'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP3_ET0_All'
        mag_press =  'MAG_SC' + sc_str +'_B_xyz_gse_MAG_PR'
        nerror = 0
        tplot_names, h_press, names = names
        IF names(0) EQ '' THEN  nerror = nerror+1 $
        ELSE BEGIN 
            get_data, names(0), data = data
            IF N_ELEMENTS(data.x) LT 2 THEN  nerror = nerror+1
        ENDELSE 
        
        tplot_names, o_press, names = names
        IF names(0) EQ '' THEN BEGIN 

            nerror = nerror+1
        ENDIF ELSE BEGIN 
            get_data, names(0), data = data
            IF N_ELEMENTS(data.x) LT 2 THEN BEGIN 

                nerror = nerror+1
            ENDIF 
        ENDELSE 
        
        tplot_names, mag_press, names = names
        IF names(0) EQ '' THEN BEGIN 
            nerror = nerror+1
        ENDIF ELSE BEGIN 
            get_data, names(0), data = data
            IF N_ELEMENTS(data.x) LT 2 THEN BEGIN 
                nerror = nerror+1
            ENDIF 
        ENDELSE 
        IF nerror EQ 0 THEN BEGIN 
            h1_press = h_press & mag_press = mag_press & o1_press = o_press
            plasma_beta, h1_press, mag_press, O1_PRESSURE = o1_press 
            get_data, 'TDMOM_EN00040_40000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_beta', data = data_b

            get_data, 'EPH_SC4_GSE_X', data = data_x
            get_data, 'EPH_SC4_GSE_Y', data = data_y
            get_data, 'EPH_SC4_GSE_Z', data = data_z

            data_xx = INTERPOL(data_x.y, data_x.x, data_b.x)
            data_yy = INTERPOL(data_y.y, data_y.x, data_b.x)
            data_zz = INTERPOL(data_z.y, data_z.x, data_b.x)

            IF found EQ  0 THEN BEGIN 
                time = data_b.x
                beta = data_b.y
                x_gse = data_xx
                y_gse = data_yy
                z_gse = data_zz
                found = 1
            ENDIF ELSE BEGIN 
                time = [time, data_b.x]
                beta = [beta, data_b.y]
                x_gse = [x_gse, data_xx]
                y_gse = [y_gse, data_yy]
                z_gse = [z_gse, data_zz]

            ENDELSE      
        ENDIF 
    ENDIF  
ENDFOR 
stop
store_data, 'beta', data = {x:time, y:beta, x_gse:x_gse, y_gse:y_gse, z_gse:z_gse}
tplot_save, 'beta', filename = 'output/o_beam/hia/beta'

stop


END
