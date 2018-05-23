pro time_snap,filen=filen,pcex=pcex,nskip=nskip,oldbck=oldbck,sizeoffile=filesize,grpsize=grpsize
;-----------------------------
;  $Id: time_snap.pro,v 1.10 2008/06/09 20:50:11 bacmj Exp $
;-----------------------------

init = 1
np = 10L
lp = 100
lm = 72
i0=0l & j0=0l
ipc=0L
l=0L
type_of_ce=0L
my_env    =0L
my_dndrft =0L
my_updrft =0L
my_prshft =0L
npopu     =0L
config='                '

tke_flag=0L

    time = 0.d
    age = 0.d
    birthtime = 0.d
    deathtime = 0.d
    config = 0.d
    maxw00 = 0.d
    maxwtf = 0.d
    status = 0.d
    rain_gauge = 0.d
    snow_gauge = 0.d
    hail_gauge = 0.d
    xmean = 0.d
    ymean = 0.d
    a0 = 0.d
    ash = 0.d
    shtop = 0.d


              ; --- scalars
      qtbken = 0.d
      qtbkde = 0.d
      shextop = 0.d
      shexbot = 0.d





if not keyword_set(filen) then filen='../fort.112'

close,1
openr,1,filen,/f77_u

ird=0
irdx=0
reads=-1


reread:

readu,1,np

for ip = 0,np-1 do begin

         ;!read in one PCE to initialize
     readu,1,ipc,l
     readu,1,i0,j0
     readu,1,type_of_ce,tke_flag
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
     readu,1,shtop

     ;! set up and read profile arrays
      x  = dblarr(l)     &   readu,1,x 
      y  = dblarr(l)     &   readu,1,y 
      z  = dblarr(l+1)   &   readu,1,z
      pe = dblarr(l+1)   &   readu,1,pe
      u  = dblarr(l)     &   readu,1,u 
      v  = dblarr(l)     &   readu,1,v 
      w  = dblarr(l)     &   readu,1,w 
      ur = dblarr(l)     &   readu,1,ur
      a  = dblarr(l)     &   readu,1,a
      q  = dblarr(l)     &   readu,1,q
      ql = dblarr(l)     &   readu,1,ql
      qi = dblarr(l)     &   readu,1,qi
      qr = dblarr(l)     &   readu,1,qr
      qs = dblarr(l)     &   readu,1,qs
      qh = dblarr(l)     &   readu,1,qh
      th = dblarr(l)     &   readu,1,th

      abck= dblarr(l)    ; &   readu,1,abck
      thbck=dblarr(l)    &   readu,1,thbck
      ubck =dblarr(l)    ; &   readu,1,ubck
      vbck =dblarr(l)    ; &   readu,1,vbck
      wbck =dblarr(l)    ; &   readu,1,wbck
      qbck =dblarr(l)    ; &   readu,1,qbck
      qlbck =dblarr(l)   ; &   readu,1,qlbck
      qibck =dblarr(l)   ; &   readu,1,qibck
      qrbck =dblarr(l)   ; &   readu,1,qrbck
      qsbck =dblarr(l)   ; &   readu,1,qsbck
      qhbck =dblarr(l)   ; &   readu,1,qhbck

;;; PCEDIAGS
              ; --- profiles
      mutt      = dblarr(l)        &   readu,1,mutt
      pracs     = dblarr(l)        &   readu,1,pracs
      tmuph     = dblarr(l)        &   readu,1,tmuph
      qmuph     = dblarr(l)        &   readu,1,qmuph
      qvsl      = dblarr(l)        &   readu,1,qvsl
      qvsi      = dblarr(l)        &   readu,1,qvsi

      dhn       = dblarr(l)        &   readu,1,dhn
      e         = dblarr(l)        &   readu,1,e
      mass      = dblarr(l)        &   readu,1,mass
      rhoa      = dblarr(l)        &   readu,1,rhoa
      rhol      = dblarr(l)        &   readu,1,rhol
      ;dzrw      = dblarr(l)        &   readu,1,dzrw

      pp        = dblarr(l+1)      &   readu,1,pp
      phyd      = dblarr(l+1)      &   readu,1,phyd
      rcdot     = dblarr(l)        &   readu,1,rcdot


              ; --- scalars
      readu,1,qtbken
      readu,1,qtbkde
      readu,1,shextop
      readu,1,shexbot

      rfall     = dblarr(l)        &   readu,1,rfall
      sfall     = dblarr(l)        &   readu,1,sfall
      hfall     = dblarr(l)        &   readu,1,hfall

      if tke_flag eq 1 then begin
         tke      = dblarr(l)        &   readu,1,tke
      endif


      zo= ( z(1:l)+z(0:l-1) )/2.0

         pce   =   {    ipc    : ipc    ,  $    
                        l      : l      ,  $    
                        i      : i0      ,  $    
                        j      : j0      ,  $    
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
                        birthtime   : birthtime   ,  $    
                        deathtime   : deathtime   ,  $    
                        age    : age    ,  $    
                        config : config ,  $    
                        maxw00 : maxw00 ,  $    
                        maxwtf : maxwtf ,  $    
                        status : status ,  $    
                        xmean  : xmean  ,  $    
                        ymean  : ymean  ,  $    
                        a0     : a0     ,  $    
                        ash    : ash    ,  $    
                        shtop  : shtop  ,  $
                        x      : x      ,  $    
                        y      : y      ,  $    
                        z      : z      ,  $    
                        zo     : zo     ,  $    
                        pe     : pe     ,  $    
                        u      : u      ,  $    
                        v      : v      ,  $    
                        w      : w      ,  $    
                        ur     : ur     ,  $    
                        a      : a      ,  $    
                        q      : q      ,  $    
                        ql     : ql     ,  $    
                        qi     : qi     ,  $    
                        qr     : qr     ,  $    
                        qs     : qs     ,  $    
                        qh     : qh     ,  $    
                        th     : th     ,  $
                        abck   : abck   ,  $    
                        thbck  : thbck  ,  $    
                        wbck   : wbck   ,  $
                        ubck   : ubck   ,  $    
                        vbck   : vbck   ,  $    
                        qbck   : qbck   ,  $    
                        qlbck  : qlbck  ,  $    
                        qibck  : qibck  ,  $    
                        qrbck  : qrbck  ,  $    
                        qsbck  : qsbck  ,  $    
                        qhbck  : qhbck  ,  $    
                        mutt   : mutt   ,  $    
                        pracs  : pracs  ,  $    
                        tmuph  : tmuph  ,  $    
                        qmuph  : qmuph  ,  $    
                        qvsl   : qvsl   ,  $    
                        qvsi   : qvsi   ,  $    
                        dhn    : dhn    ,  $
                        mass   : mass   ,  $
                        rhoa   : rhoa   ,  $
                        rhol   : rhol   ,  $
                        ;;dzrw   : dzrw   ,  $
                        pp     : pp     ,  $
                        phyd   : phyd   ,  $
                        e      : e      ,  $
                        rfall  : rfall  ,  $
                        sfall  : sfall  ,  $
                        hfall  : hfall  ,  $
                        qtbken : qtbken ,  $
                        qtbkde : qtbkde ,  $
                        rcdot  : rcdot  ,  $
                        shextop : shextop,  $
                        shexbot : shexbot   $
                                          }    

                        if tke_flag eq 1 then begin
                           pce = create_struct( pce, "tke", tke )
                        endif



      if ip eq 0 then begin
            pcex=replicate(pce,np)
            init=0
      endif

      pcex[ip]=pce

      irdx=irdx+1

   endfor




point_lun,-1,nnn
ngoto=long64( nnn )
fli  =fstat(1)
filesize = fli.size
grpsize  = ngoto



; repoint data file and back to top

reads=reads+1
ngoto=nskip*ngoto

point_lun,1,ngoto

if reads ne 1 then goto, reread


return
end
