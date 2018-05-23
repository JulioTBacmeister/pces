

close,1
openr,1,/f77,'fort.201'


last_time=-999.

while not eof(1) do begin

icl=0l & nocl=0l



      pce_ip = 0l
      pce_type=0l
      pce_time=0.

      readu,1,pce_ip
      readu,1,pce_type
      readu,1,pce_time

if pce_time ne last_time then begin
   new_time_level=1
   gr={time:pce_time}
   grx=gr
endif else begin
   new_time_level=0
   gr={time:pce_time}
endelse



   gr = create_struct( gr , 'ip'   , pce_ip )
   gr = create_struct( gr , 'type' , pce_type )

print,"  pce_ip   = ",pce_ip
print,"  pce_type = ",pce_type
print,"  pce_time = ",pce_time

      numitouch=0l & lp=0l

      readu,1,numitouch,lp
print,"   num I touch ",numitouch,"  LP ",lp
    
     touch_ips=lonarr(numitouch) 

     readu,1,touch_ips
   gr = create_struct( gr , 'touch' , touch_ips )

print," who I touch ",touch_ips
      efxs=fltarr( lp, numitouch )

      readu,1,efxs
 
  gr = create_struct( gr , 'efxs' , efxs )

print,"                .... "

if new_time_level eq 0 then begin
   str='pce_'+strtrim(string(pce_ip),2)
   grx = create_struct( grx , str , gr )
endif

last_time = pce_time

endwhile

end
