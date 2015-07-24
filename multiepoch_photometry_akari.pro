pro multiepoch_photometry_akari, inputlist, maplist=maplist, radius=radius, galactic=galactic, decimal=decimal, rinner=rinner, router=router

; By default this code takes as input a list of source coordinates 
; which are stored in inputlist with the following form:
; sname, RAHR, RAMIN, RASEC, DECDEG, DECMIN, DECSEC
;  
; However, if the keyword decimal is set, then the assumed format is:
; sname, RA, DEC - in decimal degrees. 
;
; Furthermore if the keyword galactic is set, then the assumed format is:
; sname, GLON, GLAT - in decimal degrees.
;
;
; INPUTS
; inputlist - CSV (or other tabular format) file containing the list of target sources along wiht coordinates (can be RA and DEC or GLON GLAT)
; maplist - file containing names of input HEALPix maps
;     freq [GHz] is taken from the FREQ header keyword in the extension
;     units string is taken from the TUNIT1 header keyword in the
;     extension
; Note: actual frequency and FWHM parameters are from the RIMO for DX9
; radius - radius of the source aperture in arcmin, if 0 then set to
;          Planck FWHM for this band.
; galactic - see above.
; rinner - inner radius of background annulus in units of <radius>,
;          defaults to 2.0.
; router - outer radius of background annulus in units of <radius>,
;          defaults to 3.0.

; OUTPUTS
; This procedure creates an output file named <inputlist>.photo which
; contains the list of the all-sky maps used and the aperture photometry
; results for each source.
; 
; The flux density is given in Jy.  It is measured in an aperture of i
; radius = FWHM if radius is not specified.
;
; The data columns presented in the output file are:
;    Source_Name  Map_number  GLON   GLAT   Flux (Jy)  Flux_RMS (Jy)   Median_Background_Flux (Jy)
;
; HISTORY
;
; 10-Apr-2013  P. McGehee     Corrected ten() -> tenv() error, added /decimal keyword.
; 11-Feb-2013  P. McGehee     Change of input arguments.
; 20-Sep-2012  P. McGehee     Ingested into IPAC SVN, formatting changes
;------------------------------------------------------------
    freqlist = ['30','44','70','100','143','217','353','545','857','1874','2141','2998','3331','4612','4997','11992','24983']
    freqval = [28.405889,44.072241,70.421396,100,143,217,353,545,857.,1874.,2141.,2998.,3331.,4612.,4997.,11992.,24983.]
    fwhmlist = [33.1587,28.0852,13.0812,9.88,7.18,4.87,4.65,4.72,4.39,4.3,4.0,3.8,0.62,0.65,3.8,0.97,1.02] ; fwhm in arcminutes

    if (not keyword_set(rinner)) then rinner = 2.0
    if (not keyword_set(router)) then router = 3.0

    k0 = 1.0 & k1 = rinner & k2 = router 
    apcor = ((1 - (0.5)^(4*k0^2)) - $ 
             ((0.5)^(4*k1^2) - (0.5)^(4*k2^2)))^(-1)
  
; 'galactic' overrules 'decimal' 
    if (keyword_set(galactic)) then begin 
        readcol,inputlist,sname,glon,glat,format='A,D,D'
        euler, glon, glat, ra, dec, 2
    endif else if (keyword_set(decimal)) then begin
        readcol,inputlist,sname,ra,dec,format='A,D,D'
        euler,ra,dec,glon,glat,1
    endif else begin
        readcol,inputlist,sname,rah,ram,ras,decd,decm,decs,format='A,I,I,D,I,I,D'
        ra = 15*tenv(rah,ram,ras)
        dec = tenv(decd,decm,decs)
        euler,ra,dec,glon,glat,1
    endelse

    ns = n_elements(glat)

    fd3 = -1
    fd_err3 = -1

    readcol, maplist, format='(A)', fn
    nmaps = n_elements(fn)
        PhotoResult = DBLARR[ns,nmaps*2]
    openw,1,file_basename(inputlist+'.photo'),width=200

    if (not keyword_set(radius)) then begin
        printf,1,'; A multiplicative aperture correction factor of ',apcor,' has been applied to the flux densities'
        printf,1, '; assuming that the source is a point source. Flux densities are in Jy.'
    endif else $
        printf,1,';No aperture correction factor has been applied. Flux densities are in Jy.'

    printf,1, ';'
    printf,1, ';Output format:'
    printf,1, '; Source_Name  Map_number  GLON   GLAT   Flux (Jy) Flux_RMS (Jy) Median_Background_Flux (Jy)'
    printf, 1, ';'
    printf, 1, '; Map List:'
    for i = 0, nmaps-1 do begin
        printf, 1, '; Map #',i, fn[i]
    endfor
    printf, 1, ';'

    for ct2 = 0,nmaps-1 do begin
        xtmp = mrdfits(fn[ct2], 1, hdr1)
        freq = strtrim(sxpar(hdr1, 'FREQ'),2 )
        units = strtrim(sxpar(hdr1, 'TUNIT1'),2 )
      
        idx = where(freqlist eq freq, cnt)
        if (cnt gt 0) then begin
            currfreq = freqval[idx[0]]
            if (not keyword_set(radius)) then $
                radval = fwhmlist[idx[0]] $
            else radval = radius
        endif else begin
            print, 'Invalid frequency ', freq, ' in ', fn[ct2]
            exit
        endelse

        for ct=0L,ns-1 do begin
            haperflux, fn[ct2], currfreq, fwhm, glon[ct], glat[ct], $
                1.*radval, rinner*radval, router*radval, units, $
                fd, fd_err, fd_bg, /nested,/noise_mod

            if (finite(fd_err) eq 0) then begin
                fd = -1
                fd_err = -1
            endif else begin
                if not keyword_set(radius) then begin
                    fd = fd*apcor
                    fd_err = fd_err*apcor
                endif
            endelse

  ;;;;;;;;;Let's change this next line so that it puts the data into a structure I can use for the blackbody fitting,
  ;;;;;;;;;;;rather than a huge list of photometry results
            printf,1,sname[ct],ct2,glon[ct],glat[ct],fd,fd_err,fd_bg,format='(A18,X,I2,X,F12.7,X,F12.7,X,E11.3,X,E11.3,X,E11.3)'
            PhotoResult[ct,ct2+0] = fd
            PhotoResult[ct,ct2+1] = fd_err
        endfor
    endfor

    close,1
end
