close,1
openr,1,/f77_u,'../fort.911'

for i=0,198 do begin

ipc=0L&ncl=0L&lp=0L
readu,1,ipc,ncl,lp,pctime
e=fltarr(LP)&mute=fltarr(NCL+1,LP)
readu,1,e,mute

endfor

close,1

end
