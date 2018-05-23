function ovrcircle,r1,r2,d

RR = r1
r  = r2
if r1 lt r2 then begin
   RR = r2
   r  = r1
endif

; no overlap
if d gt (r+RR) then begin
  return,0
endif

;smaller circle totally contained
if (d+r) lt RR then begin
   ar=!pi*r^2
   return,ar
endif


cos1 = (d^2+r1^2-r2^2) / (2*d*r1 ) 
cos2 = (d^2+r2^2-r1^2) / (2*d*r2 )

ang1 = 2*acos( cos1 )
ang2 = 2*acos( cos2 )


ar1 = 0.5 * r1^2 *( ang1 - sin(ang1) )
ar2 = 0.5 * r2^2 *( ang2 - sin(ang2) )

ar  = ar1 + ar2
;STOP


return,ar
end