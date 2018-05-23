pro check_energy_budget,qq=qq,be1=be1,be2=be2,ke1i=ke1i,pe1i=pe1i


    physparms,radius_e=radius_e  $
             ,grav=grav  $
             ,cpd=cp $
             ,rkap=rkap $
             ,ALHL=alhl $
             ,alhs=alhs $
             ,alhf=alhf $
             ,rgas=rgas





s=size(qq.a)&nt=s(3)&nz=s(1)&np=s(2)

epsi= 0.61
;epsi= 0.0

th00= reform( qq(0).thbck(nz-1,0) )

be1 = qq.th(*,1)/qq.thbck(*,0) + epsi *( qq.q(*,1)-qq.q(*,0) ) $      
;be1 = qq.th(*,1)/th00(0) + epsi *( qq.q(*,1)-qq.q(*,0) ) $      
   - ( qq.ql(*,1)-qq.ql(*,0) ) $
   - ( qq.qi(*,1)-qq.qi(*,0) ) $
   - ( qq.qr(*,1)-qq.qr(*,0) ) $
   - ( qq.qs(*,1)-qq.qs(*,0) ) $
   - ( qq.qh(*,1)-qq.qh(*,0) ) 

be2 = qq.th(*,2)/qq.thbck(*,0) + epsi *( qq.q(*,2)-qq.q(*,0) ) $      
   - ( qq.ql(*,2)-qq.ql(*,0) ) $
   - ( qq.qi(*,2)-qq.qi(*,0) ) $
   - ( qq.qr(*,2)-qq.qr(*,0) ) $
   - ( qq.qs(*,2)-qq.qs(*,0) ) $
   - ( qq.qh(*,2)-qq.qh(*,0) ) 


thbckz = qq.thbck(*,0)*0.

for l=1,nz-2 do begin

    thbckz(l,*) = ( qq.thbck(l-1,0) - qq.thbck(l+1,0) ) $
                / ( qq.zo(l-1,0) - qq.zo(l+1,0) )

endfor

l=0    
    thbckz(l,*) = ( qq.thbck(l,0) - qq.thbck(l+1,0) ) $
                / ( qq.zo(l,0) - qq.zo(l+1,0) )
l=nz-1
    thbckz(l,*) = ( qq.thbck(l-1,0) - qq.thbck(l,0) ) $
                / ( qq.zo(l-1,0) - qq.zo(l,0) )

H_th=100000.

pe1 = grav * be1*be1 *  qq.thbck(*,0) /thbckz / 2.
pe2 = grav * be2*be2 *  qq.thbck(*,0) /thbckz / 2.
;pe1 = grav * be1*be1 *  th00(0) /thbckz / 2.
;pe2 = grav * be2*be2 *  th00(0) /thbckz / 2.
;pe1 = grav * be1*be1 *  H_th / 2.
;pe2 = grav * be2*be2 *  H_th / 2.

ke1 = qq.w(*,1)^2 /2.
ke2 = qq.w(*,2)^2 /2.

rho = qq.mass(*,1)


ke1i = total( ke1 * qq.a(*,1) * rho , 1 )
ke2i = total( ke2 * qq.a(*,2) * rho , 1 )
pe1i = total( pe1 * qq.a(*,1) * rho , 1 )
pe2i = total( pe2 * qq.a(*,2) * rho , 1 )


temp=qq.th*0
ple = qq.pe(*,0)
plo = 0.5*(ple(1:nz,*)+ple(0:nz-1,*))
qst = qq.q*0.

for ic = 0,np-1 do begin
    temp(*,ic,*) = (( plo/1000.)^rkap) * ( qq.th(*,ic)+qq.thbck(*,0) )
    qst(*,ic,*)  = qsat(  temp(*,ic,*) , plo )
endfor

zo=qq(0).zo(*,0)

h   = cp*temp + alhl*qq.q + grav*qq.zo
hst = cp*temp + alhl*qst  + grav*qq.zo
sen = h(*,1,0)-alhl*qq(0).q(*,1) 
sen2= h(*,1,0)-alhl*qq(0).q(*,1) 

dh=-(hst(*,1,0)-h(99,1,0))
dh(where(dh lt 0))=0.
cape_init = total( rho*qq(0).a(*,1)*dh)


tpr2=fltarr(nz,nz)
tpr=sen*0.
cape2=tpr
for ll=nz-1,0,-1 do begin

tpr=sen*0
tpr(LL)=temp(LL,1,0)
qpr=qq(0).q(*,1)*0
qpr(LL)=qq(0).q(LL,1)

for l = LL-1,0,-1 do begin
   tpr(l) = (sen(l+1)-grav*qq(0).zo(L,1))/cp
   qpr(l) = qpr(l+1)
   qst0   = qsat(  temp(L,1,0) , plo(L) )
   if qpr(l) gt qst0 then begin
      dq=qpr(l)-qst0
      tpr(l) = tpr(l)+(ALHL/CP)*dq
      qpr(l) = qpr(L)-dq
   endif
   sen(l) = cp*tpr(l)+grav*qq(0).zo(L,1)
endfor

tpr2(*,LL)=tpr

hpr=cp*tpr+alhl*qpr+grav*qq.zo(*,1)

for l = ll,0,-1 do begin
    if ( (tpr(L)-temp(L,0,0) ) gt 0. ) then begin
    cape2(ll) = cape2(ll) + grav*(tpr(L)-temp(L,0,0) )* ( qq(0).ze(L)-qq(0).ze(L+1) ) $
            / temp(L,0,0)
    endif
endfor
endfor

cape2i = total( rho*qq(0).a(*,1)*cape2 )


STOP

return
end
