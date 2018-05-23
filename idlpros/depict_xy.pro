pro depict_xy,pce,initialize=initialize,xr=xr,yr=yr,maxoverlap=maxoverlap,rain=rain,level=level,slices=slices


a=pce.a
r=sqrt(a/!pi)

lp=n_elements(a)

phi=findgen(200)*!pi/100.

if keyword_set(initialize) then begin 
   oo=where( a eq max(a) )
   l=oo(0)
    x0=pce.x(l)
    y0=pce.y(l)

    ax=r(l)*cos(phi)
    ay=r(l)*sin(phi)
    plot,x0+ax,y0+ay,ps=3,xr=xr,yr=yr
endif


case 1 of 
keyword_set(maxoverlap): begin
   oo=max( pce.a, l )

    x0=pce.x(l)
    y0=pce.y(l)

    ax=r(l)*cos(phi)
    ay=r(l)*sin(phi)


       oplot,x0+ax,y0+ay
end

keyword_set(rain): begin
   oo=max( pce.a, l )

    x0=pce.x(l)
    y0=pce.y(l)

    ax=r(l)*cos(phi)
    ay=r(l)*sin(phi)


       ;oplot,x0+ax,y0+ay,colo=255-fix( -pce.rain_gauge*86400. )
       polyfill,x0+ax,y0+ay,colo=fix( -pce.rain_gauge*86400. )
end

keyword_set(level): begin
     
    l=fix(level)

    x0=pce.x(l)
    y0=pce.y(l)

    ax=r(l)*cos(phi)
    ay=r(l)*sin(phi)

       
       oplot,x0+ax,y0+ay
       xyouts,x0,y0,strtrim( string(pce.ipc) , 2) 
       
end

keyword_set(slices): begin
  for l=0,lp-1,1 do begin

    x0=pce.x(l)
    y0=pce.y(l)

    ax=r(l)*cos(phi)
    ay=r(l)*sin(phi)


       oplot,x0+ax,y0+ay
       ;;;xyouts,x0,y0,strtrim( string(pce%ip) , 2) )


  endfor 
end
endcase

;if not keyword_set(maxoverlap) then begin
;  for l=0,lp-1,5 do begin;
;
;    x0=pce.x(l)
;    y0=pce.y(l);
;
;    ax=r(l)*cos(phi)
;    ay=r(l)*sin(phi)
;
;
;       oplot,x0+ax,y0+ay
;
;  endfor 
;
;endif else begin
;
;   oo=max( pce.a, l )
;
;    x0=pce.x(l)
;    y0=pce.y(l)
;
;    ax=r(l)*cos(phi)
;    ay=r(l)*sin(phi);
;
;
;       oplot,x0+ax,y0+ay
;
;endelse


return
end
