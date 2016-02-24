function cull, bgps, file = file, verbose=verbose

; Culls a catalog based on a given set of boundary files.
  if n_elements(file) eq 0 then file = 'bounds.txt'
  readcol, file, filename, lmin, lmax, bmin, bmax, format = 'A,F,F,F,F'

  if size(bgps,/type) eq 7 then restore,bgps
  
; Trim to the root filename

  slashpos = strpos(filename, '/', /reverse_search)
  root = strmid(filename, slashpos[0]+1, 40)
  
  catroots = stregex(bgps.filename,'v2.0_ds2.*.fits',/extract)
  ;catroots = (bgps.filename)
  ;slashpos = strpos(catroots, '/', /reverse_search)
  ;catroots = strmid(catroots, slashpos[0]+1)
  catroots = catroots[uniq(catroots, sort(catroots))]

  print,"There are ",n_elements(catroots)," unique fields"

  for i = 0, n_elements(catroots)-1 do begin
    ind = where(root eq catroots[i], ct)

    if ct eq 0 then begin
        print,"No matches for ",catroots[i]
        continue
    endif
    for j = 0, ct-1 do begin
      if lmin[ind[j]] eq -1 then begin
          print,"Skipped ",catroots[i]," because of invalid lmin"
          continue
      endif
      whmatch = where(strpos(bgps.filename,catroots[i]) ge 0, nmatch)
      if lmin[ind[j]] gt lmax[ind[j]] or lmax[ind[j]]-lmin[ind[j]] gt 180 then begin
          nkeep = 0
          lower = max([lmin[ind[j]], lmax[ind[j]]])
          upper = min([lmin[ind[j]], lmax[ind[j]]])
          keep = [where(strpos(bgps.filename, catroots[i]) ge 0 and $
                       bgps.glon_max ge lower and $
                       bgps.glon_max lt 360.0 and $
                       bgps.glat_max ge bmin[ind[j]] and $
                       bgps.glat_max lt bmax[ind[j]], nkeep1),$
                 where(strpos(bgps.filename, catroots[i]) ge 0 and $
                       bgps.glon_max ge 0.0 and $
                       bgps.glon_max lt upper and $
                       bgps.glat_max ge bmin[ind[j]] and $
                       bgps.glat_max lt bmax[ind[j]], nkeep2)]
          nkeep = nkeep1+nkeep2
      endif else begin
          keep = where(strpos(bgps.filename, catroots[i]) ge 0 and $
                       bgps.glon_max ge lmin[ind[j]] and $
                       bgps.glon_max lt lmax[ind[j]] and $
                       bgps.glat_max ge bmin[ind[j]] and $
                       bgps.glat_max lt bmax[ind[j]], nkeep)
       endelse
      if keyword_set(verbose) then print,"Keeping ",strtrim(nkeep)," of ",strtrim(nmatch)," for ",catroots[i],".  Catout has ",strtrim(n_elements(catout))," elements"
;    filk = where(strpos(bgps.filename, catroots[i]) ge 0)
;    print, catroots[i], root[ind], lmin[ind[0]], lmax[ind[0]]
;    stop
      if nkeep eq 0 then stop
      if nkeep gt 0 then catout = (n_elements(catout) eq 0) ? bgps[keep] : [catout, bgps[keep]]
    endfor
  endfor

  return,catout
end
