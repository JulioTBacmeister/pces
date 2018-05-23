im=0l&lm=0l&jm=0l&np=0l&lp=0l
close,1
 openr,1,/f77_u,'fort.222'

readmore:

 readu,1,im,jm,lm


thb=fltarr(im,jm,lm)&thc=fltarr(im,jm,lm)&thg=fltarr(im,jm,lm)&ab=fltarr(im,jm,lm)&ac=fltarr(im,jm,lm)
readu,1,ab,ac   
readu,1,thg,thb,thc
readu,1,np,lp   
th1=fltarr(lp)&th2=fltarr(lp)
a1=fltarr(lp)&a2=fltarr(lp)
readu,1,a1,a2
readu,1,th1,th2

wc_m=fltarr( lm , lp )
wm_c=fltarr( lp , lm )

ppeb=fltarr( im,jm,lm+1) & ppe1=fltarr(lp+1)
pplb=fltarr( im,jm,lm) & ppl1=fltarr(lp)
kcm=lonarr(lm+1)

readu,1,wc_m
readu,1,wm_c
readu,1,kcm
readu,1,ppeb,ppe1
readu,1,pplb,ppl1

dpb=ppeb(5,5,1:lm)-ppeb(5,5,0:lm-1)
dp1=ppeb(1:lp)-ppeb(0:lp-1)


STOP

goto,readmore


end
