 close,1
 openr,1,'../fort.311',/f77_u
 close,2
 openr,2,'../fort.313',/f77_u

plot_pos,2,2,/ref,/flip,pos=pos
    if keyword_set(animation) then window,xs=1100,ys=1000,re=2


first_read=1

more:
  nx=0l&ny=0l&nz=0l
 readu,1,nx,ny,nz,pcetime

 print,"Time == ",pcetime

 if first_read eq 1 then begin
 celli=fltarr( nx,ny,nz )
 cellu=fltarr( nx,ny,nz )
 celld=fltarr( nx,ny,nz )
 dphyd=fltarr( nx,ny,nz ) 
 phyd=fltarr( nx,ny,nz+1)
 d2phyd=fltarr( nx,ny,nz+1)
 d2ph=fltarr( nx,ny,nz+1)
 dpdz=fltarr( nx,ny,nz)
 th3=fltarr( nx,ny,nz)
 u3=fltarr( nx,ny,nz)
 v3=fltarr( nx,ny,nz)
 w3=fltarr( nx,ny,nz)
 ztmg=fltarr( nx,ny,nz)
 e2ij=fltarr( nx,ny,nz)
 ;ztz3=fltarr( nx,ny,nz)
 ur3=fltarr( nx,ny,nz)
 pnh=fltarr( nx,ny,nz+1)
 pnh2=fltarr( nx,ny,nz+1)
 p2=fltarr( nx,ny,nz+1)
 p2x=fltarr( nx,ny,nz+1)
 p2y=fltarr( nx,ny,nz+1)
 p2r=fltarr( nx,ny,nz+1)

 ur2=fltarr( nx,nz)
 zax=fltarr( nx,nz)&zax2=fltarr( nx,nz)
 ducz=fltarr( nx,nz)
 dke=fltarr( nx,nz)
 div2=fltarr( nx,nz)
 exz2=fltarr( nx,nz)

 xx= fltarr(nx)&yy= fltarr(ny)&zz= fltarr(nz)

 rhol=fltarr(nz)&thbck=fltarr(nz)
 rhd_u=fltarr(nz)&rhd_s=fltarr(nz)
 thd_u=fltarr(nz)&thd_s=fltarr(nz)
 qd_u=fltarr(nz) &qd_s=fltarr(nz)
 qld_u=fltarr(nz) &qld_s=fltarr(nz)
 qid_u=fltarr(nz) &qid_s=fltarr(nz)
 qrd_u=fltarr(nz) &qrd_s=fltarr(nz)
 qsd_u=fltarr(nz) &qsd_s=fltarr(nz)
 qhd_u=fltarr(nz) &qhd_s=fltarr(nz)

 ppc = fltarr( nz+1, 2 )
 ppcz = fltarr( nz , 2 )
 pcetheta = fltarr(nz)
 pcea = fltarr(nz)
 pcew = fltarr(nz)


    first_read=0
 endif

 
 readu,1,rhol,thbck,pcetheta,pcew,pcea
 readu,1,rhd_u,thd_u,qd_u,qld_u,qid_u,qrd_u,qsd_u,qhd_u
 readu,1,rhd_s,thd_s,qd_s,qld_s,qid_s,qrd_s,qsd_s,qhd_s
 
 readu,1,ppc,ppcz

 readu,1,celli    
 readu,1,cellu    
 readu,1,celld
 readu,1,dphyd,th3
 readu,1,phyd
 readu,1,pnh
 readu,1,pnh2
 ;readu,1,w3
 ;readu,1,ur3,u3,v3
 ;readu,1,ztx3,zty3,ztz3
 readu,1,d2ph
 readu,1,p2
 readu,1,p2x
 readu,1,p2y
 readu,1,p2r
 ;readu,1,ducz,dke,ur2,zax,div2,exz2

 pdyn = pnh-pnh2


 readu,2,nx,ny,nz,pcetime
 readu,2,xx,yy,zz
 readu,2,u3
 readu,2,v3
 readu,2,w3
 readu,2,e2ij
 readu,2,ztmg
 readu,2,ur3
 


 d2ph(*,*,0:nz-1) = d2ph(*,*,0:nz-1)-ur3(*,*,0:nz-1)

 phydx=fltarr(nz+1)
 pnhx=fltarr(nz+1)
 phydx0=fltarr(nz+1)
 pnhx0=fltarr(nz+1)
 phydx1=fltarr(nz+1)
 pnhx1=fltarr(nz+1)


aa=total( total( celli , 1) , 1)
armax=max( aa , larmax )

x=findgen(nx)*1000.
y=findgen(ny)*1000.

delz0 = 300.

zo=(nz-.5-findgen(nz))*delz0
ze=(nz-findgen(nz+1))*delz0

for l=0,nz-1 do begin
    dpdz(*,*,l)=(1./rhol(l))*( pnh(*,*,l) - pnh(*,*,l+1) ) / ( ze(l) - ze(l+1) )
endfor

 for l=0,nz-1 do begin
     pnhx(l+1)  = total( total( pnh(*,*,l+1)*cellu(*,*,l),1),1) / (total( total( cellu(*,*,l),1),1)+0.001) 
     phydx(l+1) = total( total( phyd(*,*,l+1)*cellu(*,*,l),1),1) / (total( total( cellu(*,*,l),1),1)+0.001)
  endfor
 for l=0,nz-1 do begin
     pnhx0(l+1)  = total( total( pnh(*,*,l+1)*celli(*,*,larmax),1),1) / (total( total( celli(*,*,larmax),1),1)+0.001) 
     phydx0(l+1) = total( total( phyd(*,*,l+1)*celli(*,*,larmax),1),1) / (total( total( celli(*,*,larmax),1),1)+0.001)
  endfor

 for l=0,nz-1 do begin
     pnhx1(l+1)  = total( total( pnh(*,*,l+1)*cellu(*,*,larmax),1),1) / (total( total( cellu(*,*,larmax),1),1)+0.001) 
     phydx1(l+1) = total( total( phyd(*,*,l+1)*cellu(*,*,larmax),1),1) / (total( total( cellu(*,*,larmax),1),1)+0.001)
  endfor

 for l=0,nz do begin
 for j=1,ny-2 do begin
 for i=1,nx-2 do begin
     d2phyd(i,j,l) = ( phyd(i,j+1,l)+phyd(i,j-1,l)-2.*phyd(i,j,l) )+ ( phyd(i+1,j,l)+phyd(i-1,j,l)-2.*phyd(i,j,l) )
  endfor
  endfor
  endfor


 w2=  reform( w3(*,49,*) )

  zax2(*,*)=0.
  for l=1,nz-2 do begin
  for j=1,ny-2 do begin
      
      zax2(j,l) = (ur2(j,l+1)-ur2(j,l-1))/(zo(l+1)-zo(l-1)) - (w2(j+1,l)-w2(j-1,l))/(x(j+1)-x(j-1))

  endfor
  endfor
        
if keyword_set(go_to) then begin

   if pcetime lt go_to then goto,more

endif 


pnlv=(findgen(30)-15)*1000

    if not keyword_set(animation) then begin 
       stop
    endif else begin

      if not keyword_set(set) then begin
         print," spcify SET=1 or SET=2 "
             STOP
      endif 


      
      if (set eq 1) then begin
         contour,reform( dphyd(*,49,*) ),x,zo,lev=(findgen(30)-15)*.1,tit='dphyd',pos=pos(*,0),/fill
         contour,reform( phyd(*,49,*) ),x,ze,lev=(findgen(30)-15)*1000,tit='dphyd',pos=pos(*,0),/noer
  
        contour,reform( w3(*,49,*) ),x,zo,lev=(findgen(30)-15)*10,tit='W',pos=pos(*,1),/noer,/fill
        contour,reform( pnh(*,49,*) ),x,ze,lev=pnlv,tit='W',pos=pos(*,1),/noer,c_line=(pnlv lt 0)
 
        contour,reform( d2ph(*,49,*) ),x,ze,lev=(findgen(30)-15)*0.001,tit='D2PH',pos=pos(*,2),/noer,/fill
 
        ;;contour, ur2,x,zo,lev=(findgen(30)-15)*10,tit='U!ir!n',pos=pos(*,2),/noer,/fill
        
      endif
      if (set eq 2) then begin
        wlv=(findgen(30)-15)*10
       !p.charsize=1.5
          ;contour,reform( w3(*,49,*) ),x,zo,lev=(findgen(30)-15)*10,tit='W',pos=pos(*,0),/fill
          ;contour,reform( dphyd(*,49,*) ),x,zo,lev=(findgen(30)-15)*.05,tit='dphyd',pos=pos(*,0),/fill
         contour,reform( th3(*,49,*) ),x,zo,lev=(findgen(30)-15)*5.,tit='Theta*',pos=pos(*,0),/fill
         ;contour,reform( w3(*,49,*) ),x,zo,lev=wlv,pos=pos(*,0),/foll,/noer,c_line=(wlv lt 0)
 
       contour,reform( pnh(*,49,*) ),x,ze,lev=pnlv,tit='P',pos=pos(*,1),/noer,/fill
       contour,reform( pdyn(*,49,*) ),x,ze,lev=pnlv,tit='P!idyn!n',pos=pos(*,2),/noer,/fill
       contour,reform( pnh2(*,49,*) ),x,ze,lev=pnlv,tit='P!ibuoy!n',pos=pos(*,3),/noer,/fill
 
        ;;contour,reform( d2ph(*,49,*) ),x,ze,lev=(findgen(30)-15)*0.001,tit='D2PH',pos=pos(*,2),/noer,/fill
 
        ;;contour, ur2,x,zo,lev=(findgen(30)-15)*10,tit='U!ir!n',pos=pos(*,2),/noer,/fill
        
      endif

        xyouts,/norm,.5,.5,align=.5,'time = '+string(pcetime)

     ; wait,0.1
;;STOP

      if keyword_set(stop_at) then begin
         if pcetime ge stop_at then STOP
      endif




      if keyword_set(go_to) then STOP

      if keyword_set(stops) then STOP
    endelse
goto,more
 close,1
 end
