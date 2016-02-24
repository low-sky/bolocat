function object_photometry, data, hd, error, props, diam_as, fluxerr = innerflux_err, nobg = nobg, $
    ct_threshold=ct_threshold, backgrounds=backgrounds

  bootiter = 100
  
  if n_elements(diam_as) eq 0 then diam_as = 40.0
  if ~keyword_set(ct_threshold) then ct_threshold = 25
  
  getrot, hd, rot, cd, /silent
  rdhd, hd, s = h
  if h.ppbeam eq 0 then ppbeam = 1.0 else ppbeam = h.ppbeam
  if string(sxpar(hd, 'BUNIT')) eq 'JY/PIX' then ppbeam = 1.0
  if sxpar(hd, 'PPBEAM') gt 0 then ppbeam = sxpar(hd, 'PPBEAM')
  ;penguins
  ;if strtrim(string(sxpar(hd,'BUNIT'))) eq 'MJy/sr' then data*=1d6*abs(cd[0]*cd[1])*(!pi/180.0)^2*ppbeam
  ;endpenguins
  rad_pix = diam_as/abs(cd[1]*3.6d3)/2.0 ; Convert from diam -> radius
  if n_elements(rad_pix) eq 1 then rad_pix = rebin([rad_pix], n_elements(props))
  sz = size(data)
  x = findgen(sz[1])#replicatE(1, sz[2])
  y = replicate(1, sz[1])#findgen(sz[2])
  innerflux = fltarr(n_elements(props))+!values.f_nan
  innerflux_err = fltarr(n_elements(props))+!values.f_nan
  backgrounds = fltarr(n_elements(props))

  ;xokrange = [min(where(total(data,2,/nan) ne 0))-max(rad_pix)*4,max(where(total(data,2,/nan) ne 0))+max(rad_pix)*4]
  ;yokrange = [min(where(total(data,1,/nan) ne 0))-max(rad_pix)*4,max(where(total(data,1,/nan) ne 0))+max(rad_pix)*4]
  ;print,"X OK range: ",xokrange,"Y OK range: ",yokrange

  matches = 0
  for k = 0, n_elements(props)-1 do begin
; Do gcirc distance here?   
    if (props[k].maxxpix eq 0 and props[k].maxypix eq 0) then begin
        counter, k+1, n_elements(props), 'Object photometry. x=0,y=0     '+' Matches:'+string(matches,format='(6I)')+' Source: '
        continue
    endif
    ;else if (props[k].maxxpix gt xokrange[1] or props[k].maxxpix < xokrange[0] or $
    ;        props[k].maxypix gt yokrange[1] or props[k].maxypix lt yokrange[0]) then begin
    ;    counter, k+1, n_elements(props), 'Object photometry. x,y outside '+' Matches:'+string(matches,format='(6I)')+' Source: '
    ;    continue
    ;endif
    dist = sqrt((x-props[k].maxxpix)^2+(y-props[k].maxypix)^2)

    ; how big should the outer apertures be for background subtraction?
    ; 20" = 2.77776 pix
    if rad_pix[k] lt 3 then begin
        bgscale = 2
    endif else begin
        bgscale = 2 ; this is probably more effective with bgscale=1 but that's not how it was done.
    endelse
    ;print,bgscale,rad_pix[k]


    bgind = where(dist ge rad_pix[k]*bgscale and dist le rad_pix[k]*bgscale*2 and data eq data, ct)

    if ct lt ct_threshold then begin
        counter, k+1, n_elements(props), 'Object photometry. f=no match  '+' Matches:'+string(matches,format='(6I)')+' Source: '
        continue
    endif
; Set BGs equal to zero!
    if keyword_set(nobg) then begin
        background = 0.0
    endif else begin
        mmm, data[bgind], background
        backgrounds[k] = background
        ; I really want to use this.... background = median(data[bgind])
    endelse
;    background = background < 0
;    background = 0
    wtmask = dist le (rad_pix[k]-1)
    border = (dist le (rad_pix[k]+1))-wtmask
    ind = where(border)
    xborder = ind mod sz[1]
    yborder = ind / sz[1]
    border_wt = pixwt(props[k].maxxpix, props[k].maxypix, $
                      rad_pix[k], xborder, yborder)
    wtmask = float(wtmask)
    wtmask[ind] = border_wt
    ind = where(wtmask gt 0, inner_ct)

    innerflux[k] = total(wtmask[ind]*data[ind], /nan)-$
                   background*total(wtmask[ind])

    ;print,"NPIX = ",total(wtmask[ind])," aperture= ",diam_as," background= ",background," medbg= ",median(data[bgind])," flux= ",innerflux[k]," totalbackground= ",background*total(wtmask[ind])

               ; shouldn't this be inner_ct, and not threshold*2?  Swapped with the other?
    if ct lt ct_threshold*2 then begin
        counter, k+1, n_elements(props), 'Object photometry. f=no match  '+' Matches:'+string(matches,format='(6I)')+' Source: '
        continue
    endif


;;     bgrun = fltarr(bootiter)
;;     for i = 0, bootiter-1 do begin
;;       subsample = floor(randomu(seed, ct)*ct)
;;       mmm, data[bgind[subsample]], bg
;;       bgrun[i] = bg
;;     endfor
;;     bgerr = mad(bgrun)
    bgerr = 0.0
    if keyword_set(nobg) then bgerr = 0.0

; Calculate the error as the weighted error over the object * sqrt(nbeams)    
    innerflux_err[k] =  sqrt((((total(wtmask[ind]*error[ind])/$
                                total(wtmask[ind])))^2+bgerr^2)*$
                             (total(wtmask[ind])/ppbeam))
    matches += 1
    counter, k+1, n_elements(props), 'Object photometry. f='+string(innerflux(k),format='(G10.3)')+' Matches:'+string(matches,format='(6I)')+' Source: '
  endfor 
  print
  print,"Total ",matches," matches"
  if matches eq 0 then message,'ERROR: No matches.'
  innerflux = innerflux/ppbeam

  return, innerflux
end
