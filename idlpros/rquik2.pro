pro rquik2,filen=filen,dats=dats  

restart:

lp=100 ;35
nt=1000 ;000

x=fltarr(lp)
y=fltarr(lp)
Z=fltarr(lp)
Ze=fltarr(lp+1)
q=fltarr(lp)
r=fltarr(lp)
ql=fltarr(lp)
qe=fltarr(lp)
th=fltarr(lp)
thd=fltarr(lp)
the=fltarr(lp)
te=fltarr(lp)
byc=fltarr(lp)
pfz=fltarr(lp)
w=fltarr(lp+1)
wsub=fltarr(lp+1)
dthdt_e=fltarr(lp)
dqdt_e=fltarr(lp)
wl=fltarr(lp)
ef=fltarr(lp)
ur=fltarr(lp)
p=fltarr(lp+1)

rx1=fltarr(lp)
tar0=fltarr(lp)
tar1=fltarr(lp)
adot=fltarr(lp)

qst=fltarr(lp)
qst_p=fltarr(lp)
qst_t=fltarr(lp)

qq=fltarr(lp,nt)
qql=fltarr(lp,nt)
qqe=fltarr(lp,nt)
rr=fltarr(lp,nt)
bby=fltarr(lp,nt)
eef=fltarr(lp,nt)
uur=fltarr(lp,nt)
ww=fltarr(lp+1,nt)
wwsub=fltarr(lp+1,nt)
ffsub=fltarr(lp+1,nt)
wwl=fltarr(lp,nt)
pp=fltarr(lp+1,nt)
tth=fltarr(lp,nt)
tthe=fltarr(lp,nt)
tte=fltarr(lp,nt)
tthd=fltarr(lp,nt)
ppfz=fltarr(lp,nt)

ddthdt_e=fltarr(lp,nt)
ddqdt_e=fltarr(lp,nt)

qqst=fltarr(lp,nt)
qqst_p=fltarr(lp,nt)
qqst_t=fltarr(lp,nt)

rrx1=fltarr(lp,nt)
aadot=fltarr(lp,nt)

age=fltarr(nt)
mass=age
ccqmass0=age
ccqmass1=age
eeqmass0=age
eeqmass1=age
time=age
rradius=age
wtf=age
w00=age
stat=age
config='                     '

ll=1l
ipc=1l
close,1
openr,1,filen,/f77_u 
for n=0,nt-1 do begin 

readu,1,ipc,ll,age0,time0
readu,1,config
;!readu,1,cqmass0
;!readu,1,cqmass1
;!readu,1,eqmass0
;!readu,1,eqmass1
;!readu,1,radius
readu,1,maxw00
readu,1,maxwtf
readu,1,status
readu,1,x
readu,1,y 
readu,1,z
readu,1,wl
readu,1,ur
readu,1,r
readu,1,q
readu,1,ql
;!readu,1,thd
;!readu,1,te
;!readu,1,ef
;!readu,1,wsub
;!readu,1,fsub
;!readu,1,tar0
;!readu,1,ADOT
;!readu,1,p
;!readu,1,qst
;!readu,1,the
;!readu,1,qe
readu,1,rx1
;!readu,1,wsub
;!readu,1,dthdt_e
;!readu,1,dqdt_e


rr(*,n)=r
pp(*,n)=p
bby(*,n)=byc
eef(*,n)=ef
uur(*,n)=ur
ww(*,n)=w
wwsub(*,n)=wsub
;ffsub(*,n)=fsub
wwl(*,n)=wl
qq(*,n)=q
qqe(*,n)=qe
qql(*,n)=ql
tth(*,n)=th
tthd(*,n)=thd
tthe(*,n)=the
tte(*,n)=te
ppfz(*,n)=pfz
qqst(*,n)=qst
qqst_t(*,n)=qst_t
qqst_p(*,n)=qst_p
rrx1(*,n)=rx1
aadot(*,n)=adot
ddthdt_e(*,n)=dthdt_e
ddqdt_e(*,n)=dqdt_e

age(n)=age0
time(n)=time0
;mass(n)=mass0
;qmass(n)=qmass0
;!ccqmass1(n) = cqmass1
;!ccqmass0(n) = cqmass0
;!eeqmass1(n) = eqmass1
;!eeqmass0(n) = eqmass0
;!rradius(n) = radius
wtf(n) = maxwtf
w00(n) = maxw00
stat(n) = STATUS

endfor 
close,1

;goto, restart

ttht=tthd+tthe

dats = { name: filen , config:config, Z:z, ze:ze, age:age,   $
           radius:rradius, wtf:wtf , w00:w00 , status:stat , $
           rr:rr, $
           qq:qq, $
           qql:qql, $
           qqe:qqe, $
           uur:uur, $
           eef:eef, $
           wwl:wwl, $
           pp:pp, $
           rx1:rrx1, $
           adot:aadot, $
           qsat:qqst, $
           tte:tte, $
           tthe:tthe, $
           tthd:tthd, $
           ttht:ttht, dthdt_e:ddthdt_e,dqdt_e:ddqdt_e,w_e:wwsub }




;STOP


end
