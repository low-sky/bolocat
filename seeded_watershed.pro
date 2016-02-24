function seeded_watershed, data, kernels, levels = levels, $
                           nlevels = nlevels, $
                           all_neighbors = all_neighbors

  if n_elements(levels) eq 0 then $
    levels = contour_values(data, nlevels = nlevels)

  sz = size(data)
  object = ulonarr(sz[1], sz[2])
  object[kernels] = indgen(n_elements(kernels))+1

  for k = 0, n_elements(levels)-1 do begin
    mask = data ge levels[k] ; GT vs GE?
    object[kernels] = indgen(n_elements(kernels))+1
    object = dilator(object, [-1, kernels], constraint = mask, /loop)
  endfor 
  l = label_region(data eq data, all_neighbors = all_neighbors)
  ctr = max(object)+1

  for k = 1,  max(l)-1 do begin
    ind = where(l eq k)
    assigns = object[ind]
    if max(assigns) eq 0 then begin
      object[ind] = ctr
      ctr = ctr+1
    endif
  endfor
  
;  stop
; Reject unconnected pixels
  ind = where(object gt n_elements(kernels), ct)
  if ct gt 0 then object[ind] = 0

  return, object
end
