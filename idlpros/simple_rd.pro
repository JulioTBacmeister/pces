pro simple_rd,filen=filen,pces=pces,nt=nt

init = 1
np = 10L
lp = 100
ipc=0L
l=0L
my_env=0L
config='      '

close,1
openr,1,filen,/f77_u

ird=0
while ird lt nt do begin
readu,1,np
for ip = 0,np-1 do begin
;!read in one PCE
     readu,1,ipc,l,my_env
     readu,1,time
     readu,1,age
     readu,1,config
     readu,1,maxw00
     readu,1,maxwtf
     readu,1,status

     ;! set up and read profile arrays
      x = fltarr(l)     &   readu,1,x 
      y = fltarr(l)     &   readu,1,y 
      z = fltarr(l+1)   &   readu,1,z
      u = fltarr(l)     &   readu,1,u 
      v = fltarr(l)     &   readu,1,v 
      w = fltarr(l)     &   readu,1,w 
      ur= fltarr(l)     &   readu,1,ur
      a = fltarr(l)     &   readu,1,a
      q = fltarr(l)     &   readu,1,q
      ql= fltarr(l)     &   readu,1,ql
      qr= fltarr(l)     &   readu,1,qr
      th= fltarr(l)     &   readu,1,th
      e = fltarr(l)     &   readu,1,e
      pp= fltarr(l)     &   readu,1,pp

      zo= ( z(1:l)+z(0:l-1) )/2.0

 
         pce   =   {    ipc    : ipc    ,  $    
                        l      : l      ,  $    
                        my_env : my_env ,  $    
                        time   : time   ,  $    
                        age    : age    ,  $    
                        config : config ,  $    
                        maxw00 : maxw00 ,  $    
                        maxwtf : maxwtf ,  $    
                        status : status ,  $    
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
                        a      : a      ,  $    
                        q      : q      ,  $    
                        ql     : ql     ,  $    
                        qr     : qr     ,  $    
                        th     : th        }    

      if init eq 1 then begin
            pces=replicate(pce,np,nt)
            init=0
      endif
      pces[ip,ird]=pce
  
endfor

ird=ird+1

endwhile

return
end
