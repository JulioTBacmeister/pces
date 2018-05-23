close,1&openr,1,/f77_u,'../fort.511'
  lp=100
  ir=0
np=700
  ttwx=dblarr(lp,np)
  rhsx=dblarr(lp,np)
  krhx=dblarr(lp,np)
  dkex=dblarr(lp,np)
  klh1x=dblarr(lp,np)
  klh2x=dblarr(lp,np)
  klh3x=dblarr(lp,np)
  flxx=dblarr(lp+1,np)
  plxx=dblarr(lp+1,np)
  a0x=dblarr(lp,np)
  a1x=dblarr(lp,np)
  w0x=dblarr(lp,np)
  w1x=dblarr(lp,np)
  rhox=dblarr(lp,np)

  ttw=dblarr(lp)
  rhs=dblarr(lp)
  krh=dblarr(lp)
  dke=dblarr(lp)
  klh1=dblarr(lp)
  klh2=dblarr(lp)
  klh3=dblarr(lp)
  flx=dblarr(lp+1)
  plx=dblarr(lp+1)
  rho=dblarr(lp)
  a0=dblarr(lp)
  a1=dblarr(lp)
  w0=dblarr(lp)
  w1=dblarr(lp)

  
  while ( not eof(1) and ir lt np-1) do begin
     readu,1,ttw,rhs,krh,dke,klh1,klh2,klh3,flx,plx,rho,w0,w1,a0,a1
     ttwx(*,ir)=ttw
     rhsx(*,ir)=rhs
     krhx(*,ir)=krh
     dkex(*,ir)=dke
     klh1x(*,ir)=klh1
     klh2x(*,ir)=klh2
     klh3x(*,ir)=klh3
     flxx(*,ir)=flx
     plxx(*,ir)=plx
     a0x(*,ir)=a0
     a1x(*,ir)=a1
     w0x(*,ir)=w0
     w1x(*,ir)=w1
     rhox(*,ir)=rho
     ir=ir+1
  endwhile

; this looks positive definite
; plot,total( (flxx(1:100,*)-flxx(0:99,*))*w0x,1)-total(0.5*(plxx(1:100,*)-plxx(0:99,*))*w0x^2,1),total(dkex,1)-total(klh2x,1)

  
  end
  
