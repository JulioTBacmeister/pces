pro kei_rd,kei=kei,dkei=dkei,fxti=fxti,dfti=dfti,time=time,dqaw=dqaw $
          ,d2kei=d2kei,d3kei=d3kei,d4kei=d4kei
close,1
openr,1,'../fort.411',/f77_unf

time=dblarr(10000L,3)
kei=dblarr(10000L,3)
dkei=dblarr(10000L,3)
d2kei=dblarr(10000L,3)
d3kei=dblarr(10000L,3)
d4kei=dblarr(10000L,3)
fxti=dblarr(10000L,3)
dfti=dblarr(10000L,3)
dqaw=dblarr(10000L,3)
ipc=0L & pci=0L
kei0=0.d
dkei0=0.d
d2kei0=0.d
d3kei0=0.d
d4kei0=0.d
fxti0=0.d
dfti0=0.d
time0=0.d
dqaw0=0.d
ird=0L
while not eof(1) do begin
   readu,1,ipc,pci,time0
   readu,1,kei0
   readu,1,dkei0
   readu,1,fxti0
   readu,1,DFti0
   readu,1,dqaw0
   readu,1,d2kei0
   readu,1,d3kei0
   readu,1,d4kei0
   kei(ird,ipc-1)=kei0
   dkei(ird,ipc-1)=dkei0
   d2kei(ird,ipc-1)=d2kei0
   d3kei(ird,ipc-1)=d3kei0
   d4kei(ird,ipc-1)=d4kei0
   fxti(ird,ipc-1)=fxti0
   dfti(ird,ipc-1)=dfti0
   dqaw(ird,ipc-1)=dqaw0
   time(ird,ipc-1)=time0
   ird=ird+1
endwhile
close,1

kei=kei(0:ird-1,*)
dkei=dkei(0:ird-1,*)
d2kei=d2kei(0:ird-1,*)
d3kei=d3kei(0:ird-1,*)
d4kei=d4kei(0:ird-1,*)
fxti=fxti(0:ird-1,*)
dfti=dfti(0:ird-1,*)
dqaw=dqaw(0:ird-1,*)
time=time(0:ird-1,*)

return
end
