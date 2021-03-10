PRO find_dayom_from_doy, year, doy, mon, day
  
    IF FIX(year/4.) EQ year/4. THEN BEGIN
  
    day = [$
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
    
    mon = [$
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, $
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, $
            2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, $
            2,2,2,2,2,2,2,2,2,2,2,2,2,2, $
            3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, $
            3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, $
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, $
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, $
            5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, $
            5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, $
            6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, $
            6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, $
            7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, $
            7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, $
            8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, $
            8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, $
            9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, $
            9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, $
            10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, $
            10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, $
            11,11,11,11,11,11,11,11,11,11,11,11,11,11,11, $
            11,11,11,11,11,11,11,11,11,11,11,11,11,11,11, $
            12,12,12,12,12,12,12,12,12,12,12,12,12,12,12, $
            12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12 ]
    
  ENDIF ELSE BEGIN
    day = [$
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30, $
            1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, $
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
    
    mon = [$
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, $
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, $
            2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, $
            2,2,2,2,2,2,2,2,2,2,2,2,2, $
            3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, $
            3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, $
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, $
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, $
            5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, $
            5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, $
            6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, $
            6,6,6,6,6,6,6,6,6,6,6,6,6,6,6, $
            7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, $
            7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, $
            8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, $
            8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8, $
            9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, $
            9,9,9,9,9,9,9,9,9,9,9,9,9,9,9, $
            10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, $
            10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10, $
            11,11,11,11,11,11,11,11,11,11,11,11,11,11,11, $
            11,11,11,11,11,11,11,11,11,11,11,11,11,11,11, $
            12,12,12,12,12,12,12,12,12,12,12,12,12,12,12, $
            12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12 ]
    
  ENDELSE
  
  mon = mon(doy-1)
  day = day(doy-1)
  
END

PRO read_omni_hr_idl, filename, omni_time, omni_data

restore, f = filename+'.idl'

END


PRO read_omni_hr_ascii, filename, omni_time, omni_data

omni_hr_struct = {omni_hr, $
                  Year:0, $
                  Day:0, $
                  Hour:0, $
                  Minute:0, $
                  ID_IMF_spacecraft:0, $
                  ID_SW_Plasma_spacecraft:0, $
                  Numb_points_IMF_averages:0, $
                  Numb_points_Plasma_averages:0, $
                  Percent_OF_Interpolation:0, $
                  Timeshift:0, $
                  RMS_Timeshift:0, $
                  RMS_Min_var:0.0, $
                  Time_btwn_observations:0, $
                  Field_magnitude_average:0.0, $
                  BX_GSE:0.0, $
                  BY_GSE:0.0, $
                  BZ_GSE:0.0, $
                  BY_GSM:0.0, $
                  BZ_GSM:0.0, $
                  RMS_SD_B:0.0, $
                  RMS_SD_field:0.0, $
                  Speed:0.0, $
                  Vx:0.0, $
                  Vy:0.0, $
                  Vz:0.0, $
                  Proton_Density:0.0, $
                  Temperature:0.0, $
                  Flow_pressure:0.0, $
                  Electric_field:0.0, $
                  Plasma_beta:0.0, $
                  Alfven_mach_number:0.0, $
                  SC_Xgse:0.0, $
                  SC_Ygse:0.0, $
                  SC_Zgse:0.0, $
                  BSN_location_Xgse:0.0, $
                  BSN_location_Ygse:0.0, $
                  BSN_location_Zgse:0.0, $
                  AE_index:0, $
                  AL_index:0, $
                  AU_index:0, $
                  SYM_D:0, $
                  SYM_H:0, $
                  ASY_D:0, $
                  ASY_H:0, $
                  PCN_index:0.0 $
                 }

omni_data = replicate(omni_hr_struct, 9000000)
omni_time = DBLARR(9000000)

dummy = ''

OPENR, unit, filename+'.lst', /GET_LUN
ii = 0l
WHILE NOT EOF(unit) DO BEGIN
    
    READF, unit, dummy
    
    READS, dummy, omni_hr_struct
    
    find_dayom_from_doy, omni_hr_struct.year, omni_hr_struct.day, month, day
    
    omni_time(ii) = time_double($
                    STRING(omni_hr_struct.year, FORMAT = '(i4.4)')+'-' + $
                    STRING(month, FORMAT = '(i2.2)')+'-' + $
                    STRING(day, FORMAT = '(i2.2)')+'/' + $
                    STRING(omni_hr_struct.hour, FORMAT = '(i2.2)') + ':' + $
                    STRING(omni_hr_struct.minute, FORMAT = '(i2.2)') + ':00')
    omni_data(ii) = omni_hr_struct
    ii = ii + 1
    
ENDWHILE

CLOSE, unit, /FORCE
FREE_LUN, unit, /FORCE

omni_time = omni_time(0:ii-1)
omni_data = omni_data(0:ii-1)

save, omni_time,  omni_data, f = filename+'.idl'

END


;+
;  PROCEDURE: read_omni
;
;  PURPOSE: To read ascci data as downloaded from OMNIweb site
;
;  MODIFICATIONS:
;   2006-11-27  New 'month' format is read
;
;-
PRO read_omni, HR = HR, $
               ALL = ALL, $
               BX_GSE = BX_GSE, $
               BY_GSE = BY_GSE, $
               BZ_GSE = BZ_GSE, $
               BY_GSM = BY_GSM, $
               BZ_GSM = BZ_GSM, $
               BETA = BETA, $
               MACH = MACH, $
               PRESSURE = PRESSURE, $
               DENSITY = DENSITY, $
               TEMPERATURE = TEMPERATURE, $
               VELOCITY = VELOCITY, $
               AE = AE, $
               Ey = Ey

path = GETENV('OMNI_DATA')
  
;    1  Epoch Time                      DD-MMM-YYYY_hr:mm   A17
;    2  Bartels Rotation Number                             I4
;    3  ID for IMF spacecraft                               I2
;    4  ID for SW Plasma spacecraft                         I2
;    5  # fine time scale IMF PTS                           I4
;    6  # fine time scale plasma PTS                        I4
;    7  Field Magnitude Avg.            nT                  F5.1
;    8  Magnitude of avg. field vector  nT                  F5.1
;    9  Lat. Angle of AV.               Deg                 F5.1
;   10  Long. Angle of AV.              Deg                 F5.1
;   11  Bx, GSE                         nT                  F5.1
;   12  By, GSE                         nT                  F5.1
;   13  Bz, GSE                         nT                  F5.1
;   14  By, GSM                         nT                  F5.1
;   15  Bz, GSM                         nT                  F5.1
;   16  sigma-|B|                       nT                  F5.1
;   17  sigma-B                         nT                  F5.1
;   18  sigma-Bx                        nT                  F5.1
;   19  sigma-By                        nT                  F5.1
;   20  sigma-Bz                        nT                  F5.1
;   21  Proton temperature              Deg K               F8.0
;   22  Proton density                  Per cc              F5.1
;   23  Flow speed                      Km/s                F5.0
;   24  Flow longitude                  Deg                 F5.1
;   25  Flow latitude                   Deg                 F5.1
;   26  Alpha/prot. ratio                                   F5.3
;   27  Flow pressure                   nPa                 F5.2
;   28  sigma-T                         Deg K               F8.0
;   29  sigma-n                         Per cc              F5.1
;   30  sigma-V                         Km/s                F5.0
;   31  sigma-phi-V                     Deg                 F5.1
;   32  sigma-theta-V                   Deg                 F5.1
;   33  sigma-ratio                                         F5.3
;   34  Electric Field                  mV/m                F6.2
;   35  Plasma beta                                         F6.2
;   36  Alfen mach number                                   F5.1
;   37  Kp*10                                               I2
;   38  Sunspot number                                      I3
;   39  DST Index                       nT                  I5
;   40  AE-index                        nT                  I4
;   41  PROT Flux > 1 MEV               1/(SQcm-ster-s)     F9.2
;   42  PROT Flux > 2 MEV               1/(SQcm-ster-s)     F8.2
;   43  PROT Flux > 4 MEV               1/(SQcm-ster-s)     F8.2
;   44  PROT Flux >10 MEV               1/(SQcm-ster-s)     F8.2
;   45  PROT Flux >30 MEV               1/(SQcm-ster-s)     F8.2
;   46  PROT Flux >60 MEV               1/(SQcm-ster-s)     F8.2
;   47  M'SPH Flux Flag                                     I2
;
;NOTE: Between each word is a single blank separating the values
;
;DATA FORMAT=ASCII                   RECORD SIZE=302
  
IF NOT KEYWORD_SET(HR) THEN BEGIN
    
    omni_struct = {omni, $
                    year:0, $
                    doy:0, $
                    hour:0, $
                    bartels_rot_num:0, $
                    id_imf_sc:0, $
                    id_sw_sc:0, $
                    fin_tim_imf:0, $
                    fin_tim_plasma:0, $
                    field_mag_avg:0.0, $
                    mag_avg_field_vec:0.0, $
                    lat_ang:0.0, $
                    long_ang:0.0, $
                    bx_gse:0.0, $
                    by_gse:0.0, $
                    bz_gse:0.0, $
                    by_gsm:0.0, $
                    bz_gsm:0.0, $
                    sigma_totb:0.0, $
                    sigma_b:0.0, $
                    sigma_bx:0.0, $
                    sigma_by:0.0, $
                    sigma_bz:0.0, $
                    proton_temp:0.0, $
                    proton_density:0.0, $
                    flow_speed:0.0, $
                    flow_long:0.0, $
                    flow_lat:0.0, $
                    alpha_prot_ratio:0.0, $
                    flow_pressure: 0.0, $
                    sigma_t:0.0, $
                    sigma_n:0.0, $
                    sigma_v:0.0, $
                    sigma_phi_v:0.0, $
                    sigma_theta_v:0.0, $
                    sigma_ratio:0.0, $
                    electric_field:0.0, $
                    plasma_beta:0.0, $
                    alfven_mach:0.0, $
                    kp_index:0, $
                    sunspot_num:0, $
                    dst_index:0l, $
                    ae_index:0l, $
                    pf1:0.0, $
                    pf2:0.0, $
                    pf4:0.0, $
                    pf10:0.0, $
                    pf30:0.0, $
                    pf60:0.0, $
                    msph_flux_flag:0, $
                    ap_index:0, $
                    f10_7_index:0.0, $
                    pc_index:0.0, $
                    AL_index:0, $
                    AU_index:0 $
                   }
    
    omni_data = replicate(omni_struct, 9000000)
    omni_time = DBLARR(9000000)
    
    dummy = ''
    
    month_str = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', $
                 'Sep', 'Oct', 'Nov', 'Dec']
    
    get_timespan, tt
    t_s = tt(0)
    t_e = tt(1)
    
    t1 = time_struct(t_s)
    t2 = time_struct(t_e)
    year_s = t1.year
    year_e = t2.year
    
    no_years = (year_e - year_s) + 1
    
    ii = 0l
    FOR iyear = 0, no_years-1 DO BEGIN
        
        year_str = STRING(year_s + iyear, FORMAT = '(i4.4)')
        OPENR, unit, path + '/ascii/new/new/omni2_' + year_str + '.lst', /GET_LUN
        
        WHILE NOT EOF(unit) DO BEGIN
            
            READF, unit, dummy
            
            
            READS, dummy, omni_struct
            
            find_dayom_from_doy, omni_struct.year, omni_struct.doy, month, day
            
            omni_time(ii) = time_double($
                            STRING(omni_struct.year, FORMAT = '(i4.4)')+'-' + $
                            STRING(month, FORMAT = '(i2.2)')+'-' + $
                            STRING(day, FORMAT = '(i2.2)')+'/' + $
                            STRING(omni_struct.hour, FORMAT = '(i2.2)')+':00:00')
            omni_data(ii) = omni_struct
            
            ii = ii + 1
            
        ENDWHILE
        CLOSE, unit, /FORCE
        FREE_LUN, unit, /FORCE
        
    ENDFOR
    
    omni_data = omni_data(0:ii-1)
    omni_time = omni_time(0:ii-1)
    
    store_data, 'Bx_gse', $
                data = {x:omni_time, y:omni_data.bx_gse}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'By_gse', $
                data = {x:omni_time, y:omni_data.by_gse}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'Bz_gse', $
                data = {x:omni_time, y:omni_data.bz_gse}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'By_gsm', $
                data = {x:omni_time, y:omni_data.by_gsm}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'Bz_gsm', $
                data = {x:omni_time, y:omni_data.bz_gsm}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'plasma_beta', $
                data = {x:omni_time, y:omni_data.plasma_beta}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'alfven_mach', $
                data = {x:omni_time, y:omni_data.alfven_mach}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'flow_pressure', $
                data = {x:omni_time, y:omni_data.flow_pressure}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'proton_density', $
                data = {x:omni_time, y:omni_data.proton_density}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'proton_temp', $
                data = {x:omni_time, y:omni_data.proton_temp}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'flow_speed', $
                data = {x:omni_time, y:omni_data.flow_speed}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'Sunspot_number', $
                data = {x:omni_time, y:omni_data.sunspot_num}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'Dst_Index', $
                data = {x:omni_time, y:omni_data.dst_index}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'AE_Index', $
                data = {x:omni_time, y:omni_data.ae_index}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'Kp_Index', $
                data = {x:omni_time, y:omni_data.kp_index}, $
                dlim = {panel_size:2, ylog:0}
    
    store_data, 'F10_7_Index', $
                data = {x:omni_time, y:omni_data.f10_7_index}, $
                dlim = {panel_size:2, ylog:0}

    store_data, 'Ey', $
                 data = {x:omni_time, y:omni_data.electric_field}, $
	         dlim = {panel_size:2, ylog:0}

ENDIF ELSE BEGIN

    
    get_timespan, tt
    t_s = tt(0)
    t_e = tt(1)
    
    t1 = time_struct(t_s)
    t2 = time_struct(t_e)
    year_s = t1.year
    year_e = t2.year
    
    no_years = (year_e - year_s) + 1
    
    FOR iyear = 0, no_years-1 DO BEGIN
        
        year_str = STRING(year_s + iyear, FORMAT = '(i4.4)')

        filename =  path + '/ascii/new/new/hr/omni_min_' + year_str
        
        ff = ''
        ff = file_search(filename+'.idl')
        IF ff NE '' THEN BEGIN

            read_omni_hr_idl, filename, omni_hr_time_one_year, omni_hr_data_one_year

        ENDIF ELSE BEGIN
            
            read_omni_hr_ascii, filename, omni_hr_time_one_year, omni_hr_data_one_year
            
        ENDELSE
            
        IF iyear EQ 0 THEN BEGIN
            
            omni_hr_time = omni_hr_time_one_year
            omni_hr_data = omni_hr_data_one_year

        ENDIF ELSE BEGIN

            omni_hr_time = [omni_hr_time, omni_hr_time_one_year]
            omni_hr_data = [omni_hr_data, omni_hr_data_one_year]

        ENDELSE

    ENDFOR
    
    varstr = 'OMNI_HR_'

    IF KEYWORD_SET(BX_GSE) OR KEYWORD_SET(ALL) THEN BEGIN

        f_fill = WHERE(omni_hr_data.bx_gse NE 9999.99)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).bx_gse
        store_data, varstr+'Bx_gse', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(BY_GSE) OR KEYWORD_SET(ALL) THEN BEGIN
        
        f_fill = WHERE(omni_hr_data.by_gse NE 9999.99)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).by_gse
        store_data, varstr+'By_gse', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(BZ_GSE) OR KEYWORD_SET(ALL) THEN BEGIN

        f_fill = WHERE(omni_hr_data.bz_gse NE 9999.99)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).bz_gse
        store_data, varstr+'Bz_gse', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(BY_GSM) OR KEYWORD_SET(ALL) THEN BEGIN
        
        f_fill = WHERE(omni_hr_data.by_gsm NE 9999.99)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).by_gsm
        store_data, varstr+'By_gsm', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(BZ_GSM) OR KEYWORD_SET(ALL) THEN BEGIN

        f_fill = WHERE(omni_hr_data.bz_gsm NE 9999.99)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).bz_gsm
        store_data, varstr+'Bz_gsm', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(BETA) OR KEYWORD_SET(ALL) THEN BEGIN

        f_fill = WHERE(omni_hr_data.plasma_beta NE 999.990)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).plasma_beta
        store_data, varstr+'plasma_beta', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(MACH) OR KEYWORD_SET(ALL) THEN BEGIN

        f_fill = WHERE(omni_hr_data.alfven_mach_number NE 999.900)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).alfven_mach_number
        store_data, varstr+'alfven_mack_number', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(PRESSURE) OR KEYWORD_SET(ALL) THEN BEGIN

        f_fill = WHERE(omni_hr_data.flow_pressure NE 99.9900)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).flow_pressure
        store_data, varstr+'flow_pressure', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(DENSITY) OR KEYWORD_SET(ALL) THEN BEGIN
        
        f_fill = WHERE(omni_hr_data.proton_density NE 999.990)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).proton_density
        store_data, varstr+'proton_density', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(TEMPERATURE) OR KEYWORD_SET(ALL) THEN BEGIN

        f_fill = WHERE(omni_hr_data.temperature NE 9999999.0)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).temperature
        store_data, varstr+'temperature', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(VELOCITY) OR KEYWORD_SET(ALL) THEN BEGIN
        
        f_fill = WHERE(omni_hr_data.speed NE 99999.9)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).speed
        store_data, varstr+'flow_speed', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(AE) OR KEYWORD_SET(ALL) THEN BEGIN
        store_data, varstr+'AE_Index', $
                    data = {x:omni_hr_time, y:omni_hr_data.ae_index}, $
                    dlim = {panel_size:2, ylog:0}
    ENDIF
    
    IF KEYWORD_SET(Ey) OR KEYWORD_SET(ALL) THEN BEGIN

        f_fill = WHERE(omni_hr_data.Electric_field NE 999.990)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).Electric_field
        store_data, varstr+'Ey', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}

    ENDIF

    IF KEYWORD_SET(ID_SW) OR KEYWORD_SET(ALL) THEN BEGIN

        f_fill = WHERE(omni_hr_data.ID_SW_Plasma_spacecraft NE 999.990)
        xarr = omni_hr_time(f_fill)
        yarr = omni_hr_data(f_fill).ID_SW_Plasma_spacecraft
        store_data, varstr+'SC_ID_SW', $
                    data = {x:xarr, y:yarr}, $
                    dlim = {panel_size:2, ylog:0}

    ENDIF



    
ENDELSE
  
END
