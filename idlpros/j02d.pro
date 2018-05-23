;------------------------------
; Bessel J_0(kr) specrtum of
; Del^2 phyd
;
;   for use after .run rdtoyp3
;------------------------------

nx=n_elements(x)
ny=n_elements(y)

r=fltarr( nx , ny )



for j=0,ny-1 do begin
FOR i=0,nx-1 do begin
  
  r(i,j)=sqrt(  x(i)^2 + y(j)^2 )

endfor
endfor




ks=0.1+findgen(701)/10.
aa=ks*0.


for i=0,700 do begin
  k = ks(i)
  j0r = beselj( k*reform(r)/4000. , 0 )

  j02d=reform( j0r, nx,ny )

  aa(i) = total( j02d* d2ph )

endfor

end
     


