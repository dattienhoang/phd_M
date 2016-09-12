PRO MSIM_3D

print, 'BEGIN M-MODULATION SIMULATION (in 3-DIMENSIONS)'
print, '...set parameters...'
nvec = [3];[1,2,3,4,5,7,10,50,100,1000];number of vectors
iter = 2500D;iterations
;making sure to use radians
proj = findgen(180) * (!pi / 180)
print, '...iterate for an ensemble of:'
;------------------------------------------------------------------------------
FOR i=0, n_elements(nvec)-1 DO BEGIN
  print, '......', nvec[i], ' vectors'
  print, '.............generate random angles'
  ;vecx = fltarr(nvec[i],iter) & vecy = fltarr(nvec[i],iter)
  ;vecang = fltarr(nvec[i],iter)
  rndang = imsl_random(nvec[i]*iter, /sphere, parameters=3)
  ;coordinates should be in a column of (x,y,z)...just remove the z-coord
  rndang = rndang[*,0:1]
  ;now find the angle and magnitude of the vector in x-y plane
  rndmag = fltarr(nvec[i]*iter)
  FOR j=0, n_elements(rndmag)-1 DO BEGIN
    rndmag[j] = sqrt((rndang[j,0]^2)+(rndang[j,1]^2));5;
  ENDFOR
  ;print, rndmag
  ;ang = asin(rndang[*,1] / rndang[*,0]) & print, ang & rndang = temporary(ang)      ;radians
  ang = atan(rndang[*,1], rndang[*,0]) & print, ang & rndang = temporary(ang)      ;radians
  
  ;from here the code is pretty similar, but needs to factor the magnitude in!
  ;reform the random angle to how the rest of the data will look for ease...
  print, '.............quickly reform the data'
  vecang = fltarr(nvec[i],iter);,n_elements(proj))
  FOR j=0L, n_elements(rndang)-1 DO vecang[j] = rndang[j]
  rndang = 0
  vecmag = fltarr(nvec[i],iter);,n_elements(proj))
  FOR j=0L, n_elements(rndmag)-1 DO vecmag[j] = rndmag[j]
  rndmag = 0
  print, '.............take vector projections'
  int = fltarr(nvec[i],iter,n_elements(proj))
  ;now take the projection of each vector onto every line!
  FOR j=0L, n_elements(proj)-1 DO BEGIN
    FOR k=0L, nvec[i]-1 DO BEGIN
      FOR l=0L, iter-1 DO BEGIN
        ;projection in 2D of a vector onto a line is r*cos(\theta)
        ;int[k,l,j] = cos(abs(vecang[k,l] - proj[j]))^2
        int[k,l,j] = vecmag[k,l]^2 * cos(abs(vecang[k,l] - proj[j]))^2
      ENDFOR
    ENDFOR
  ENDFOR
  
  
  ;now fit each experiment to get an M for each iteration of vector count
  ;M = (I_max -I_min) / (I_max + I_min)
  ;first need to sum up all the contributions of the individual vectors
  print, '.............fit an M value'
  m_pre = fltarr(n_elements(proj), iter)
  FOR j=0L, n_elements(proj)-1 DO BEGIN
    FOR k=0L, iter-1 DO BEGIN
      m_pre[j,k] = total(int[*,k,j])
    ENDFOR
  ENDFOR
  m_fin = fltarr(iter)
  FOR j=0L, iter-1 DO BEGIN
    imax = max(m_pre[*,j]) & imin = min(m_pre[*,j])
    m = (imax - imin) / (imax + imin)
    m_fin[j] = m
  ENDFOR
  m = 0
  
  ;need to print a file for thsi vector count...
  print, '.........print a file'
  fileparts = strcompress('/home/administrator/Desktop/Mdist_' + string(nvec[i]) + $
              'vec.txt', /remove_all)
  openw, lun, fileparts, /get_lun
  FOR j=0L, iter-1 DO printf, lun, m_fin[j]
  close, lun & free_lun, lun
  
ENDFOR

ENDPROGRAM:
print, 'END'
END
