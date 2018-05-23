function piecequad,xdata,ydata,xi 

yi = fltarr(n_elements(xi)); % initialize
N = n_elements(xdata)-1;
i = 1 ;
for j=2,N  do begin
   ; determine polynomial on interval [xdata(j-2), xdata(j)]
    x1 = xdata(j-2); 
    x2 = xdata(j-1); 
    x3 = xdata(j);
    f1 = ydata(j-2); 
    f2 = ydata(j-1); 
    f3 = ydata(j);

    dx=x2-x1

    AA = f1
   
    ;QQ =  (1./dx)*( 1-dx/2.)/(dx^2/3.-1)
   
    ;BB = (f2-f1) / ( dx + QQ )
    ;CC = BB*QQ


    CC = ( f3 -f1 - (f2-f1)*(x3-x1)/(x2-x1) ) / ( (x3-x1)*(x3-x2) )
    BB = ( f2 -f1 - CC*(x2-x1)*(x2-x1) )/(x2-x1) 


    ; set function values for any xi(i) that lie in this interval:
    oo = where( xi ge x1 and xi le x2 ) 
    if min(oo) gt -1 then begin
          yi(oo) = AA + BB*(xi(oo)-x1) + CC*(xi(oo)-x1)*(xi(oo)-x1) ;
    endif


endfor

return,yi
end
