;+
;FUNCTION:      get_cis_hia_data, prod_num, sat, t, frencheff=frencheff, time=time
;INPUT:
;       prod_num: Product number
;       sat: Satellite number (1-4)
;       time: If set, return time array of data values (index ignored)
;
;PURPOSE:   Retrieve CIS/HIA data.
;
;CREATED BY:    Peter Schroeder
;LAST MODIFICATION:  September 2000
;  04/15/02 - more prducts added MF
;  14/01/03 - more prducts added MF
;-
FUNCTION get_cis_hia_data, prod_num, sat,  time = time, frencheff = frencheff

IF NOT keyword_set(t) THEN t = 0

CASE prod_num OF
;    0: RETURN, get_cis_hia_p0(sat, t, frencheff=frencheff, time=time)
;    1: RETURN, get_cis_hia_p1(sat, t, frencheff=frencheff, time=time)
    2: RETURN, get_cis_hia_p2p4(sat, t, prod_num, time = time)
    4: RETURN, get_cis_hia_p2p4(sat, t, prod_num, time = time)
;    5: RETURN, get_cis_hia_p5(sat, t, frencheff=frencheff, time=time)
    6: RETURN, get_cis_hia_p6_p15_p17(prod_num, sat, time = time, frencheff = frencheff)
;    7: RETURN, get_cis_hia_p7(sat, t, frencheff=frencheff, time=time)
    8: RETURN, get_cis_hia_3d(prod_num, sat, time = time, frencheff = frencheff)
;    9: RETURN, get_cis_hia_p9_p18(prod_num, sat, time=time, frencheff=frencheff)
;    10: RETURN, get_cis_hia_p10_11_13_14(prod_num, sat,  frencheff=frencheff, time=time)
;    11: RETURN, get_cis_hia_p10_11_13_14(prod_num, sat,  frencheff=frencheff, time=time)
;    12: RETURN, get_cis_hia_p12(sat, t, frencheff=frencheff, time=time)
;    13: RETURN, get_cis_hia_p10_11_13_14(prod_num, sat,  frencheff=frencheff, time=time)
;    14: RETURN, get_cis_hia_p10_11_13_14(prod_num, sat,  frencheff=frencheff, time=time)
    15: RETURN, get_cis_hia_p6_p15_p17(prod_num, sat, time = time, frencheff = frencheff )
;    16: RETURN, get_cis_hia_p16(sat, t, frencheff=frencheff, time=time)
    17: RETURN, get_cis_hia_p6_p15_p17(prod_num, sat, time = time, frencheff = frencheff )
;    18: RETURN, get_cis_hia_p9_p18(prod_num, sat, time = time, frencheff = frencheff)
;    19: RETURN, get_cis_hia_p19(sat, t, frencheff=frencheff, time=time)
;    20: RETURN, get_cis_hia_p20(sat, t, frencheff=frencheff, time=time)
    21: RETURN, get_cis_hia_3d(prod_num, sat, time = time, frencheff = frencheff)
    22: RETURN, get_cis_hia_3d(prod_num, sat, time = time, frencheff = frencheff)
    23: RETURN, get_cis_hia_3d(prod_num, sat, time = time, frencheff = frencheff)
    24: RETURN, get_cis_hia_3d(prod_num, sat, time = time, frencheff = frencheff)
;    61: RETURN, get_cis_hia_p61(sat, t, frencheff=frencheff, time=time)
;    62: RETURN, get_cis_hia_p62(sat, t, frencheff=frencheff, time=time)
  ENDCASE
  
  print, 'No Product Matching Product Number ', prod_num
RETURN,0
END



