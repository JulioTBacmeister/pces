pro depict_xyz,pce,initialize=initialize,range=range,noshafts=noshafts

a=pce.a
r=sqrt(a/!pi)

lp=n_elements(a)

phi=findgen(200)*!pi/100.


if not keyword_set(range) then begin
         ; limits for small volume around
         ; individual PCE 
     maxr=  sqrt( max( pce.a ) / !pi )
     oo=where( a eq max(a) )
     l=oo(0)
     x0=pce.x(l)
     y0=pce.y(l)
     minx = x0 - 1.2*r(l)
     maxx = x0 + 1.2*r(l)
     miny = y0 - 1.2*r(l)
     maxy = y0 + 1.2*r(l)
     minz = -1000.
     maxz = 20000.
endif else begin

    minx = range.x(0) & maxx = range.x(1)
    miny = range.y(0) & maxy = range.y(1)
    minz = range.z(0) & maxz = range.z(1)

endelse





if keyword_set(initialize) then begin 

   ; Set x data scale to range
   !X.S = [-(minx), 1.]/(maxx - (minx))  
   ; Set y data scale to range
   !Y.S = [-(miny), 1.]/(maxy - (miny))  
   ; Set z data scale to range
   !Z.S = [-(minz), 1.]/(maxz - (minz))  

endif


  for l=0,lp-1 do begin

    x0=pce.x(l)
    y0=pce.y(l)

    ax=r(l)*cos(phi)
    ay=r(l)*sin(phi)
    az=ay*0+pce.zo(l)

   pccolo=0
   if (pce.type_of_ce eq 3) then pccolo=2

   draw=1
   if (pce.type_of_ce eq 3) and keyword_set(noshafts) then draw=0

   if draw eq 1 then begin
        if (pce.zo(l) lt 15000.) and ( ( pce.ql(l)*1000. gt 0.5) or ( pce.qr(l)*1000. gt 0.5) or (l eq lp-1) ) then plots,x0+ax,y0+ay,az,/t3d,/data,colo=pccolo
         ;if (pce.zo(l) lt 15000.) and ( pce.a(l) gt 0.002*pce.a0 ) then plots,x0+ax,y0+ay,az,/t3d,/data
   endif

  endfor 


return
end
