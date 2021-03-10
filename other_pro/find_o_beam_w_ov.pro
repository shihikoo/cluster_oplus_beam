PRO find_o_beam_w_ov

       ;--------------------------------------------
       ;Read plasama_beta and O+ velocity into array 
       ;---------------------------------------------
sc_str = '4'
beam_v = -30
;---Read plasama_beta data
 p_beta = 'TDMOM_EN00040_40000_SC'+ sc_str $
          +'_MTPRESSURE_SP0_ET0_All_O1_beta'
      
 get_data, p_beta, data = pb
 
 time_pb = pb.x
 data_pb = pb.y
        
;---Read velocity data 
 o_vx = 'TDMOM_EN00040_40000_SC'+ sc_str +'_MTVELOCITY_SP3_ET0_All_X'
 
 average_tplot_variable, o_vx, '120', /new ;average velocity of O+
 get_data, o_vx+'_AVG120', data = ovx
          
 time_ovx = ovx.x               ; get the data from the average data ovx
 data_ovx = ovx.y
                
       ;------------------------------------------------------------ 
       ;If there is no O+ beam:        data_beam = 0
       ;If O+ beam in lobe(<0.05):            data_beam = 1
       ;If O+ beam in plasma boundary(0.05~1): data_beam = 2
       ;If O+ beam in plasama sheet(>1):   data_beam = 3
       ;-------------------------------------------------------------
          
 data_ovx = INTERPOL(data_ovx, time_ovx, time_pb)
 time_ovx = time_pb
          
 data_beam = DBLARR(N_ELEMENTS(time_pb))
 time_beam = time_pb
            
 FOR iii = 0L, N_ELEMENTS(time_beam)-1 DO BEGIN
             
     v_judge = data_ovx(iii)
     b_judge = ABS(data_pb(iii))
                
     IF v_judge LE beam_v AND b_judge GE 1 THEN BEGIN 
         data_beam(iii) = 3 
     ENDIF
              
     IF v_judge LE beam_v AND b_judge GE 0.05 AND b_judge LT 1 THEN BEGIN 
         data_beam(iii) = 2 
     ENDIF  
              
     IF v_judge LE beam_v AND b_judge GE 0 AND b_judge LE 0.05 THEN BEGIN 
         data_beam(iii) = 1 
     ENDIF  
          
     IF v_judge GT beam_v THEN BEGIN 
         data_beam(iii) = 0
     ENDIF 
             
 ENDFOR 
          
 o_beam = 'TDMOM_EN00040_40000_SC'+ sc_str +'_BEAMJUDGEMENT_SP3_ET0_All'
 str = {x: time_beam, y: data_beam}
 store_data, o_beam, data = str, dlim = {psym:3}
END  
