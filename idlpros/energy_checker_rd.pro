f='../fort.114'
close,1
;write(114) pce%time,ke0,ke1,ke2,pe0,pe1,kesrc1,pesrc1,pesrc2

  timex = 0.d
  ke0x = 0.d
  ke1x = 0.d
  ke2x = 0.d
  pe0x = 0.d
  peix = 0.d
  pefx = 0.d
  keix = 0.d
  kefx = 0.d
  pe1x = 0.d
  pe11x = 0.d
  tesrcbx = 0.d
  kesrcbx = 0.d
  pesrcbx = 0.d
  pesrcxx = 0.d
  pesrctx = 0.d
  kesrctx = 0.d
  pesrcvx = 0.d
  kesrcvx = 0.d
  
time   = dblarr(10000L)

kei    = dblarr(10000L)
kef    = dblarr(10000L)
pei    = dblarr(10000L)
pef    = dblarr(10000L)


ke0    = dblarr(10000L)
ke1    = dblarr(10000L)
ke2    = dblarr(10000L)
pe0    = dblarr(10000L)
pe1    = dblarr(10000L)
pe11    = dblarr(10000L)
pesrcb = dblarr(10000L)
pesrcx = dblarr(10000L)
kesrcb = dblarr(10000L)
tesrcb = dblarr(10000L)
pesrct = dblarr(10000L)
kesrct = dblarr(10000L)
pesrcv = dblarr(10000L)
kesrcv = dblarr(10000L)

ird=0
openr,1,/f77_u,f

while not eof(1) do begin

  readu,1,timex,ke0x,ke1x,ke2x,pe0x,pe1x,kesrcbx,pesrcbx,pesrcxx,pe11x,kesrctx,pesrctx,keix,peix,kefx,pefx,tesrcbx,kesrcvx,pesrcvx
  time(ird)   =  timex
  ke0(ird)    =  ke0x
  ke1(ird)    =  ke1x
  ke2(ird)    =  ke2x
  pe0(ird)    =  pe0x
  pe1(ird)    =  pe1x
  pe11(ird)   =  pe11x
  pesrcb(ird) =  pesrcbx
  pesrcx(ird) =  pesrcxx
  kesrcb(ird) =  kesrcbx
  pesrct(ird) =  pesrctx
  kesrct(ird) =  kesrctx

  pei(ird)    =  peix
  pef(ird)    =  pefx
  kei(ird)    =  keix
  kef(ird)    =  kefx

  pesrcv(ird) =  pesrcvx
  kesrcv(ird) =  kesrcvx
  tesrcb(ird) =  tesrcbx

  
  ird=ird+1

endwhile
close,1

  time=time(0:ird-1)  
  ke0=ke0(0:ird-1)  
  ke1=ke1(0:ird-1) 
  ke2=ke2(0:ird-1) 
  pe0=pe0(0:ird-1) 
  pe1=pe1(0:ird-1)  
  pe11=pe11(0:ird-1) 
  pesrcb=pesrcb(0:ird-1)
  pesrcx=pesrcx(0:ird-1)
  kesrcb=kesrcb(0:ird-1)
  tesrcb=tesrcb(0:ird-1)
  pesrct=pesrct(0:ird-1)
  kesrct=kesrct(0:ird-1)
  pesrcv=pesrcv(0:ird-1)
  kesrcv=kesrcv(0:ird-1)


  kei=kei(0:ird-1)  
  kef=kef(0:ird-1) 
  pei=pei(0:ird-1)  
  pef=pef(0:ird-1) 


  te0=ke0+pe0
  te2=ke2+pe1
  tei=kei+pei
  tef=kef+pef

  kerhs = kesrcv+kesrcb+kesrct  
  perhs = pesrcv+pesrcb+pesrct  

  ; ??/??/??
  ; tef-tei and -(perhs+kerhs) should be identical
  ; 04/25/17
  ; ke0 and kei are b4b id as are pe0 and pei
  ; just based on stored w0 and "beta" (kei,pei)
  ; vs calculated at stage 1 (ke0,pe0)
  

STOP

end
