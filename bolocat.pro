pro bolocat, filein, props = props, zero2nan = zero2nan, obj = obj, $
             watershed = watershed, delta = delta, $
             all_neighbors = all_neighbors, thresh = thresh, $
             expand = expand, minpix = minpix, clip = clip, grs_run = grs, $
             absthresh = absthresh, absdelta = absdelta, $
             absexpand = absexpand, grs_dir = grs_dir, errgen = errgen, $
             error = error, minbeam = minbeam, $
             residual = residual, round = round, corect = corect, $
             sp_minpix = sp_minpix, id_minpix = id_minpix, $
             idmap = idmap, smoothmap = smoothmap, beamuc = beamuc, $
             labelmask_in=labelmask_in, apertures=apertures, $
             ct_threshold=ct_threshold, friends = friends
;+
; NAME:
;    BOLOCAT
; PURPOSE:
;    Generate a catalog from a BOLOCAM image.
; CALLING SEQUENCE:
;    
;    BOLOCAT, filename, props = props [, /zero2nan, obj = obj, 
;             /watershed, /clip, delta = delta, /all_neighbors
;             thresh = thresh, expand = expand, minpix
;             = minpix, labelmask_in=labelmask_in
;
; INPUTS:
;    FILENAME -- Name of the file to process.
;
; KEYWORD PARAMETERS:
;    /ZERO2NAN -- Replace values of zero with NaNs assuming that a
;                 zero represents bad data.  This keyword should
;                 (MUST) be set if zero is the blanking value.
;    /WATERSHED -- Use a seeded watershed decomposition on regions to
;                  determine objects
;    /CLIP -- Identify objects by clipping with no decomposition.
;    DELTA -- Saddle point criterion parameter: Set this parameter to
;             the difference (IN UNITS OF THE LOCAL RMS) a local
;             maximum must be above a saddle-point to represent a
;             unique object.  Default: 2
;    /ALL_NEIGHBORS -- Swith that controls the number of neighbors a
;                      pixel has.  Default is 4 neighbors; set the
;                      switch if a pixel should have 8 neighbors.
;    THRESH -- Initial thresholding value for determining significant
;              emission in units of the local RMS.  Default: 3
;    MINPIX -- Significant regions with fewer pixels than MINPIX are
;              rejected.  Default: 10.  
;    EXPAND -- After rejecting small regions, the remaining regions
;              are expanded to include all connected emission down to
;              this threshold (in units of the local RMS).  Default: 2.
;    ABSDELTA, ABSEXPAND, ABSTHRESH -- Function as the DELTA, EXPAND,
;                                      and THRESH keywords excep
;    GRS_RUN -- Calculate the GRS spectrum along the line of sight of
;               each core.
;    GRS_DIR -- Location of the GRS data cubes.
;    RESIDUAL -- An array the same size as the data containing an
;                estimate of signal-free noise.
;    ERROR -- An array the same size as the data image containing an
;             estimate of the standard deviation at every point.
;    ROUND -- Size of rounding element (radius).  Set to zero for no rounding.
;    BEAMUC -- FRACTIONAL beam size uncertainty.
;    labelmask_in -- An input labelmask from which source properties 
;                will be derived.  This prevents objectid from being called
;    apertures -- A 3-element list of aperture diameters in arcseconds to use
;           *WARNING* As of 3/3/2013, these will be stored in flux_40, flux_80,
;           flux_120, etc.  no matter what the true diameters are.
;    ct_threshold -- Default to 25.  Minimum number of pixels in a source to count
;           it as a match.  Should be ~ ppbeam?  Passed to object_photometry
;
; OUTPUTS:
;   PROPS -- An array of structures with each element corresponding to
;            a signficant object in the BOLOCAM image.
;
;
; OPTIONAL OUTPUTS:
;   OBJ -- An object mask of the same dimensions as the input image
;          with the pixels corresponding to the kth element of props
;          labeled with the value k.
;
; MODIFICATION HISTORY:
;
;       Tue May 29 12:22:12 2007, Erik <eros@yggdrasil.local>
;
;		Documented and tidied up.
;
;-

  corect = 0
  if not file_test(filein, /regular) then begin
    message, 'Error: File not found!', /con
    return
  endif
  

; Set default to do Seeded watershed  
  if n_elements(watershed) eq 0 and n_elements(clip) eq 0 then $
     watershed = 1b

; Assume no uncertainty in the beam size
  if n_elements(beamuc) eq 0 then beamuc = 0.0

  data = mrdfits(filein, 0, hd)

; If blanked data are set to 0 move them to NaNs
  if keyword_set(zero2nan) then begin 
    badind = where(data eq 0, ct)
    if ct gt 0 then data[badind] = !values.f_nan
  endif

  if n_elements(apertures) eq 0 then begin
      apertures = [40.,80.,120.]
  endif else if n_elements(apertures) ne 3 then begin
      message,"Error: Need 3 aperture sizes."
  endif
  
  
  
  if n_elements(error) ne n_elements(data) then begin  
    if strcompress(string(sxpar(hd, 'XTEN1')), /rem) eq 'MapError' and (not keyword_set(errgen)) then begin
      err = mrdfits(filein, 1, ehd)
      if n_elements(err) gt 1 then begin 
        message, 'Using Error Extension for Noise Variance Estimate!', /con
        err = median(err, 5)
      endif else begin
        err = bolocam_emap(datA)
        err = median(err, 11)
      endelse
    endif
    if n_elements(residual) eq n_elements(data) then err = bolocam_emap2(residual, box = 2) 
    if n_elements(err) eq n_elements(data) then error = err else begin
      error = bolocam_emap(data)
      error = median(error, 11)
    endelse
  endif
; Swap in SFL for GLS in headers
  if stregex(sxpar(hd, 'CTYPE1'), 'GLS', /bool) then begin
        ct1 = sxpar(hd, 'CTYPE1')
        ct2 = sxpar(hd, 'CTYPE2')
        pos = stregex(sxpar(hd, 'CTYPE1'), 'GLS')
        sxaddpar, hd, 'CTYPE2', strmid(ct2, 0, pos)+'SFL'
        sxaddpar, hd, 'CTYPE1', strmid(ct1, 0, pos)+'SFL'
  endif
    
; if the pixels per beam is not specified in the header, 
; add a header keyword that contains the pixels per beam 
; (if beam FWHM not in header, use PLW beam size)
    if sxpar(hd, 'PPBEAM') eq 0 then begin
      ;the major and minor FWHM of the beam in degrees (BMAJ and BMIN)

      if sxpar(hd,'BMAJ') eq 0 then sxaddpar, hd, 'BMAJ', 35.2/3600.
      if sxpar(hd,'BMIN') eq 0 then sxaddpar, hd, 'BMIN', 35.2/3600.

      bmaj = sxpar(hd, 'BMAJ')
      bmin = sxpar(hd, 'BMIN')

      ;beam major axis position angle
      sxaddpar, hd, 'BPA', 0.0

      ;getrot derives "the counterclockwise rotation angle, and the X and Y scale
;     factors of an image, from a FITS image header."... rot will contain
      ;"Scalar giving the counterclockwise rotation of NORTH in DEGREES 
;               from the +Y axis of the image." and cdv will contain
      ;"2 element vector giving the scale factors in DEGREES/PIXEL in 
;               the X and Y directions.   CDELT[1] is always positive, whereas
;               CDELT[0] is negative for a normal left-handed coordinate system,
;               and positive for a right-handed system. "
      ;basically since the images are not rotated at all, cdv[0] and cdv[1] are
      ;equivalent to cdelt1 and cdelt2 in the header, in that they contain the size
      ;of a pixel on either side in degrees
      getrot, hd, rot, cdv, /silent
   
      ppbeam = abs((bmaj*bmin)/(cdv[0]*cdv[1])*$
                   2*!pi/(8*alog(2)))

      sxaddpar, hd, 'PPBEAM', ppbeam

    endif

;penguins
;  getrot, hd, rot, cd, /silent
;  if strtrim(string(sxpar(hd,'BUNIT'))) eq 'MJy/sr' then data*=1d6*abs(cd[0]*cd[1])*(!pi/180.0)^2*ppbeam
;  sxaddpar,hd,'BUNIT','Jy/beam'
;  if strtrim(string(sxpar(hd,'BUNIT'))) eq 'MJy/sr' then print, 'YAHOO'
;endpenguins

  if n_elements(minpix) eq 0 then begin
    if n_elements(minbeam) gt 0 then begin
      minpix = minbeam*sxpar(hd, 'PPBEAM')
    endif else begin
      minpix = 2*sxpar(hd, 'PPBEAM')
    endelse
  endif 

; Call object identification /segementation routine
  if n_elements(idmap) eq n_elements(data) then map = idmap else map = data

  ; Allow the label mask to be specified instead of derived
  if n_elements(labelmask_in) gt 0 then obj = labelmask_in else begin
      obj = objectid(map, error = error, watershed = watershed, delta = delta, $
                     all_neighbors = all_neighbors, thresh = thresh, $
                     expand = expand, minpix = minpix, absdelta = absdelta, $
                     absthresh = absthresh, absexpand = absexpand, $
                     round = round, sp_minpix = sp_minpix, $
                     id_minpix = id_minpix, original = data, friends = friends)
  endelse

  if n_elements(obj) ne n_elements(data) then begin
    if n_elements(labelmask_in) gt 0 then message,"Labelmask_in was specified.",/info
    message,"Size of labelmask ("+string(n_elements(obj))+") and input data ("+string(n_elements(data))+") are different.",/con
    return
  endif
  if max(obj) gt 0 then begin
; Call catalog routine
    props = propgen(data, hd, obj, error, smoothmap = smoothmap)
    corect = n_elements(props) 
; Do photometry in annuli (apertures are in diameters)
    props.flux_40 = object_photometry(data, hd, error, props, apertures[0], $
                                      fluxerr = fe40, ct_threshold=ct_threshold, backgrounds=backgrounds)
    props.eflux_40 = sqrt((fe40)^2+4*beamuc^2*(props.flux_40)^2)
    props.bg_40 = backgrounds
    
    props.flux_40_nobg = object_photometry(data, hd, error, props, apertures[0], $
                                           fluxerr = fe40, /nobg, ct_threshold=ct_threshold)
    props.eflux_40_nobg = sqrt((fe40)^2+4*beamuc^2*(props.flux_40_nobg)^2)
    
    props.flux_80 = object_photometry(data, hd, error, props, apertures[1], $
                                      fluxerr = fe80, ct_threshold=ct_threshold, backgrounds=backgrounds)
    props.eflux_80 = sqrt((fe80)^2+4*beamuc^2*(props.flux_80)^2)
    props.bg_80 = backgrounds
    props.flux_120 = object_photometry(data, hd, error, props, apertures[2], $
                                       fluxerr = fe120, ct_threshold=ct_threshold, backgrounds=backgrounds)
    props.eflux_120 = sqrt((fe120)^2+4*beamuc^2*(props.flux_120)^2)
    props.bg_120 = backgrounds
    props.flux_80_nobg = object_photometry(data, hd, error, props, apertures[1], $
                                      fluxerr = fe80, /nobg, ct_threshold=ct_threshold)
    props.eflux_80_nobg = sqrt((fe80)^2+4*beamuc^2*(props.flux_80_nobg)^2)
    props.flux_120_nobg = object_photometry(data, hd, error, props, apertures[2], $
                                       fluxerr = fe120, /nobg, ct_threshold=ct_threshold)
    props.eflux_120_nobg = sqrt((fe120)^2+4*beamuc^2*(props.flux_120_nobg)^2)
    
    props.flux_obj = object_photometry(data, hd, error, props, props.rad_as_nodc*2, fluxerr = feobj, ct_threshold=ct_threshold, backgrounds=backgrounds)
    props.flux_obj_err = feobj
; Fill in basic properties that were used in the analysis
; Array compatibility?
    props.exp_thresh = expand[0]
    props.delta = delta[0]
    props.threshold = thresh[0]
    props.minpix = minpix[0]
    props.decomp_alg = (keyword_set(watershed)) ? 'WATERSHED' : 'CLIP'
    props.filename = filein
    
  endif else begin
      corect = 0
      message,"No objects found in bolocat (obj has no IDs)",/info
  endelse

  if keyword_set(grs) then begin
; Tag each core with a GRS spectrum if available and desired
    grslookup, props, hd, obj = obj, data = data, grs_dir = grs_dir
  endif

;  help,props,/struct
  return
end
