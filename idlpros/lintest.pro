pro lintest,p=p,notes=notes,fname=fname

common psmngr_names,ps$,jpg$,publish
common psmngr_dimens,dims
common psmngr_state,state

if not keyword_set(fname) then fname='anidlplot'

a0=max(p(0).a)
rad=sqrt(a0/!pi)
lm=n_elements(p(0).a)
dt=p(1).time-p(0).time
dz=p(0).z(0:lm-1)-p(0).z(1:lm)

af=abs(fft(p(0).th,-1))
im=min(where(af eq max(af)))
lz=(max(p.zo)-min(p.zo))/im

hh=2.4*rad & dd=lz

omg2=(0.012)^2*(2.*!pi/hh)^2 /((2.*!pi/hh)^2 + (2.*!pi/dd)^2) 

period=2.*!pi/sqrt(omg2) 

publish={dir:'/home/bacmj/Documents/pce_notes/',fname:fname}

psmngr,/ps,/in
loadct,10

plot,p(0).th,p(0).zo,pos=[.1,.5,.2,.8],ytit='Z(m)',tit='!6init !7h!6',xticks=2,xtit='K'
plot,p(0).a,p(0).zo,pos=[.22,.5,.32,.8],ycharsi=.001,tit='!6init !7a!6',xticks=1,xtit='m!e2!n',/noer, $
       xr=[0.,a0*1.2]

contour,transpose(p.th),p.time,p(0).zo,pos=[.4,.5,.7,.8],/noer,xtit='Time (s)',tit='!7h!6(t,z)',lev=(findgen(30)-15.)*1.,/fill
contour,transpose(p.th),p.time,p(0).zo,pos=[.4,.5,.7,.8],/noer,xtit='Time (s)',tit='!7h!6(t,z)',lev=(findgen(30)-15.)*1.

plot,dz,p(0).zo,pos=[.775,.5,.875,.8],xticks=1,/noer,xtit='m',tit='dz',ycharsi=.001

xyouts,/norm,.2,.3,  "A0 =    "+string(a0)
xyouts,/norm,.2,.275,"Radius ="+string(rad) 
xyouts,/norm,.2,.250,"T-step ="+string(dt)
xyouts,/norm,.2,.225,"N levs ="+string(lm)

;xyouts,/norm,.3,.150,"theor. Period assuming L=2.4r and N=0.012 :"+string(period)+' s'
xyouts,/norm,.3,.125,"Configuration  : "+strtrim(p(0).config,2)

if keyword_set(notes) then begin
   xyouts,/norm,.3,.09,notes,size=.75
endif
  
psmngr,/ps,/cl,/di
return
end
