pro sixcharge,f=f,x=x,y=y,gx=gx,gy=gy
; 2D Potential for 6-tuplet of charges arranged
;
;        -1     +2     -1
;
;      ---------------------
;
;        -1     +2     -1
;




; charged particles
q0={q:0.,x:0.,y:0.}

; ns is nmber of sixtuplets of charges
; each sixtuplet is further out and weaker
ns=5
qs=replicate( q0 , ns*6 )




for n=0,ns-1 do begin 
qs(n*6+0).x=-0.5 & qs(n*6+0).y=-0.5-n & qs(n*6+0).q=-1./(n+1.) 
qs(n*6+1).x= 0.  & qs(n*6+1).y=-0.5-n & qs(n*6+1).q= 2./(n+1.)  
qs(n*6+2).x= 0.5 & qs(n*6+2).y=-0.5-n & qs(n*6+2).q=-1./(n+1.)  

qs(n*6+3).x=-0.5 & qs(n*6+3).y= 0.5+n & qs(n*6+3).q=-1./(n+1.) 
qs(n*6+4).x= 0.  & qs(n*6+4).y= 0.5+n & qs(n*6+4).q= 2./(n+1.) 
qs(n*6+5).x= 0.5 & qs(n*6+5).y= 0.5+n & qs(n*6+5).q=-1./(n+1.) 

endfor

nx=n_elements(x)
ny=n_elements(y)

f = fltarr( nx, ny )

for n=0,ns*6-1    do begin
for j=0,ny-1 do begin
for i=0,nx-1 do begin

    rr = sqrt ( ( x(i) - qs(n).x )^2 +  ( y(j) - qs(n).y )^2 )  
   
    if rr lt 0.1 then rr=0.1

    f(i,j)  = f(i,j) + qs(n).q * alog( rr ) ; ( 1. / rr ) 

    gx(i,j) = gx(i,j) + qs(n).q * x(i)/rr^2

    gy(i,j) = gy(i,j) + qs(n).q * y(j)/rr^2


endfor
endfor
  print, " one charge "
endfor


return
end
