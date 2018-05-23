

close,1

if not keyword_set(fn) then fn='../fort.114'
lm=100

print,' reading FN=',fn


if not keyword_set(skipto) then skipto=0

openr,1,/f77,fn

time=0.

w0=fltarr(lm)
w1=fltarr(lm)
p=fltarr(lm+1)
pz=fltarr(lm)
byc=fltarr(lm)
a=fltarr(lm)
rhl=fltarr(lm)

ze=lm-findgen(lm+1)
zo=lm-findgen(lm)+.5

moreread:

readu,1, time , a0
readu,1, W0
readu,1, A
readu,1, RHL
readu,1, P
readu,1, PZ
readu,1, BYC
readu,1, W1 

if keyword_set(plots) then begin

plot_pos,3,2,/ref,/flip,pos=pos


   if not keyword_set(win) then win=0
   window,re=2,win

    plot,w1,zo,pos=pos(*,0),tit='W'
    oplot,w0,zo,ps=-4

    plot,pz,zo,/noer,tit='PZ',pos=pos(*,1)

    plot,p,ze,/noer,tit='P',pos=pos(*,2),xr=[-10,10]*1e3,/xst

    plot,rhl,zo,/noer,tit='RH',pos=pos(*,3)

    plot,a/a0,zo,/noer,tit='A',pos=pos(*,4),xr=[0,4]
    plot,byc*rhl,zo,/noer,tit='BYC*RHL',pos=pos(*,5)

     xyouts,/norm,.2,.5,fn
     xyouts,/norm,.5,.5,'time= '+string(time)    
         wait,.01


endif

 

print,'Time =',time

if time gt skipto then STOP

goto,moreread


end
