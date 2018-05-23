pro simple_rd_1,filen=filen,pces=pcex,nt=nt,ipce=ipce
;-----------------------------
;  $Id: simple_rd_1.pro,v 1.16 2007/12/21 22:58:12 bacmj Exp $
;-----------------------------

init = 1
np = 10L
lp = 100
ipc=0L
l=0L
type_of_ce=0L
my_env    =0L
my_dndrft =0L
my_updrft =0L
my_prshft =0L
npopu     =0L
config='                '

if not keyword_set(filen) then filen='fort.112'

close,1
openr,1,filen,/f77_u

ird=0
irdx=0
while ird lt nt do begin
readu,1,np
for ip = 0,np-1 do begin
;!read in one PCE

     readu,1,ipc,l
     readu,1,type_of_ce
     readu,1,my_env
     readu,1,my_dndrft
     readu,1,my_prshft
     readu,1,my_updrft
     readu,1,npopu

     readu,1,time
     readu,1,age
     readu,1,birthtime
     readu,1,deathtime
     readu,1,config
     readu,1,maxw00
     readu,1,maxwtf
     readu,1,status
     readu,1,rain_gauge
     readu,1,snow_gauge
     readu,1,hail_gauge
     readu,1,xmean
     readu,1,ymean
     readu,1,a0
     readu,1,ash

     if (ipc eq ipce) then begin
     ;print, " constructing ",ipc,age,time

     ;! set up and read profile arrays
      x  = fltarr(l)     &   readu,1,x 
      y  = fltarr(l)     &   readu,1,y 
      z  = fltarr(l+1)   &   readu,1,z
      u  = fltarr(l)     &   readu,1,u 
      v  = fltarr(l)     &   readu,1,v 
      w  = fltarr(l)     &   readu,1,w 
      ur = fltarr(l)     &   readu,1,ur
      a  = fltarr(l)     &   readu,1,a
      q  = fltarr(l)     &   readu,1,q
      ql = fltarr(l)     &   readu,1,ql
      qr = fltarr(l)     &   readu,1,qr
      th = fltarr(l)     &   readu,1,th

      e  = fltarr(l)     &   readu,1,e
      pp = fltarr(l+1)   &   readu,1,pp
      phyd= fltarr(l+1)  &   readu,1,phyd
      atf= fltarr(l)     &   readu,1,atf
      aex= fltarr(l)     &   readu,1,aex
      ahal= fltarr(l)    &   readu,1,ahal
      ftf = fltarr(l)    &   readu,1,ftf

      ubck =fltarr(l)    &   readu,1,ubck
      vbck =fltarr(l)    &   readu,1,vbck
      wbck =fltarr(l+1)  &   readu,1,wbck
      thbck=fltarr(l)    &   readu,1,thbck
      qbck=fltarr(l)     &   readu,1,qbck
      qlbck=fltarr(l)    &   readu,1,qlbck
      qibck=fltarr(l)    &   readu,1,qibck


      thbc2=fltarr(l)    &   readu,1,thbc2
      uentr=fltarr(l)    &   readu,1,uentr
      tke=fltarr(l)      &   readu,1,tke
      rh =fltarr(l)      &   readu,1,rh
      revap=fltarr(l)    &   readu,1,revap
      pfz=fltarr(l)      &   readu,1,pfz
      byc=fltarr(l)      &   readu,1,byc
      bycp=fltarr(l)     &   readu,1,bycp
      ipentr=fltarr(l,5) &   readu,1,ipentr
      ipmerg=fltarr(50)  &   readu,1,ipmerg

      zo= ( z(1:l)+z(0:l-1) )/2.0


 
         pce   =   {    ipc    : ipc    ,  $    
                        l      : l      ,  $    
                        type_of_ce : type_of_ce ,  $    
                        my_env     : my_env     ,  $    
                        my_updrft  : my_updrft  ,  $    
                        my_dndrft  : my_dndrft  ,  $    
                        my_prshft  : my_prshft  ,  $    
                        npop       : npopu      ,  $    
                        rain_gauge : rain_gauge ,  $    
                        snow_gauge : snow_gauge ,  $    
                        hail_gauge : hail_gauge ,  $    
                        time   : time   ,  $    
                        age    : age    ,  $    
                        config : config ,  $    
                        maxw00 : maxw00 ,  $    
                        maxwtf : maxwtf ,  $    
                        status : status ,  $    
                        xmean  : xmean  ,  $    
                        ymean  : ymean  ,  $    
                        a0     : a0     ,  $    
                        x      : x      ,  $    
                        y      : y      ,  $    
                        z      : z      ,  $    
                        zo     : zo     ,  $    
                        u      : u      ,  $    
                        v      : v      ,  $    
                        w      : w      ,  $    
                        ur     : ur     ,  $    
                        e      : e      ,  $    
                        pp     : pp     ,  $    
                        phyd   : phyd   ,  $    
                        atf    : atf    ,  $    
                        aex    : aex    ,  $    
                        ahal   : ahal   ,  $    
                        a      : a      ,  $    
                        q      : q      ,  $    
                        ql     : ql     ,  $    
                        qr     : qr     ,  $    
                        th     : th     ,  $
                        ue     : uentr  ,  $
                        tke    : tke    ,  $
                        rh     : rh     ,  $
                        revap  : revap  ,  $
                        pfz    : pfz    ,  $
                        byc    : byc    ,  $
                        bycp   : bycp   ,  $

                        qbck   : qbck   ,  $    
                        qlbck  : qlbck  ,  $    
                        qibck  : qibck  ,  $    
                        ubck   : ubck   ,  $    
                        vbck   : vbck   ,  $    
                        thbck  : thbck  ,  $    
                        thbc2  : thbc2  ,  $    
                        wbck   : wbck      $    
                                          }    

      if init eq 1 then begin
            pcex=replicate(pce,nt)
            init=0
      endif
      pcex[irdx]=pce
      irdx=irdx+1
  endif else begin
  
      ;print, " skipping ",ipc

     readu,1  ;ipc,l
     readu,1  ;type_of_ce
     readu,1  ;my_env
     readu,1  ;my_dndrft
     readu,1  ;my_prshft
     readu,1  ;my_updrft
     readu,1  ;npopu

     readu,1  ;time
     readu,1  ;age
     readu,1  ;birthtime
     readu,1  ;deathtime
     readu,1  ;config
     readu,1  ;maxw00
     readu,1  ;maxwtf
     readu,1  ;status
     readu,1  ;rain_gauge
     readu,1  ;snow_gauge
     readu,1  ;hail_gauge
     readu,1  ;xmean
     readu,1  ;ymean
     readu,1  ;a0
     readu,1 

     ; set up and read profile arrays
       readu,1  ; ,x 
      readu,1  ;y 
       readu,1  ;z
       readu,1  ;u 
       readu,1  ;v 
     readu,1  ;w 
      readu,1  ;ur
        readu,1  ;a
        readu,1  ;q
        readu,1  ;ql
     readu,1  ;qr
        readu,1  ;th

      readu,1  ;e
       readu,1  ;pp
        readu,1  ;phyd
        readu,1  ;atf
      readu,1  ;aex
       readu,1  ;ahal
       readu,1  ;ftf

       readu,1  ;ubck
       readu,1  ;vbck
       readu,1  ;wbck
       readu,1  ;thbck
       readu,1  ;qbck
       readu,1  ;qlbck
        readu,1  ;qibck


      readu,1  ;thbc2
       readu,1  ;uentr
        readu,1  ;tke
      readu,1  ;rh
        readu,1  ;revap
        readu,1  ;pfz
       readu,1  ;byc
        readu,1  ;bycp
       readu,1  ;ipentr
       readu,1  ;ipmerg

  endelse

endfor

ird=ird+1

endwhile

return
end
