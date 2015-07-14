im_akari9  = read_fit("im_akari9.fits",hdr_akari9)
im_akari18 = read_fits("im_akari18.fits",hdr_akari18)
im_spire250 = read_fits("im_spire250.fits".hdr_spire250)



img_ref = im_akari9
hdr_ref = hdr_akari9

img_old = im_spire250
hdr_old = hdr_spire250

  ;; We use the IDLastro procedure HASTROM to reproject.

 HASTROM, img_old, hdr_old, img_new, hdr_new, hdr_ref, $
           MISSING=!VALUES.D_NaN
  img_stack = DBLARR[*,*,2]
  img_stack[*,*,0] = im_akari9[*,*]
  img_stack[*,*,1] = im_akari18
  img_stack[*,*,2]
