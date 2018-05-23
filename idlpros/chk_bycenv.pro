close,1
openr,1,/f77,'fort.115'

ncl=0l&lp=0l
ipc=0l
for i=0,10000 do begin


readu,1,ipc
readu,1,ncl,lp
readu,1,time,age

xxp = fltarr( lp )

xxm = fltarr( ncl,lp)
wgt = fltarr( ncl, lp )

readu,1,xxp
readu,1,wgt
readu,1,xxm

STOP
endfor






end
