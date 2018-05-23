noplot=1

if not keyword_set(window_initd) then ini
window_initd=1
runthru=1
stdcolo
;f='/home/bacmj/SCMrun/arm_scmx_A_scm_physdev_Eros-beta7p14-fin_PCES/fort.112'

if not keyword_set(filename) then begin 
   f='../fort.112'
endif else begin
   f=filename
endelse
if not keyword_set(olddiags) then olddiags=0

time_snap,f=f,p=p,nsk=0,old=olddiags,size=fsize,grp=gsize
nreads=long( fsize / gsize )


;range={x:[-120000.,20000.],y:[-120000.,20000.],z:[-1000.,100000.]}

if keyword_set(only2pces) then begin
range={x:[40000.,80000.],y:[20000.,60000.],z:[-1000.,100000.]}
endif else begin
range={x:[0000.,100000.]+450000.,y:[00000.,100000.]+450000.,z:[-1000.,100000.]}
endelse


;range={x:[-60000.,-20000.],y:[-60000.,-20000.],z:[-1000.,100000.]}
;range={x:[-120000.,-60000.],y:[-120000.,-60000.],z:[-1000.,100000.]}



T3D, /RESET, ROT=[-45, 6, 0], PERS=4.

p0=p
zo=p0(0).zo

npces=n_elements(p0.ipc)

;plot,p.xmean,p.ymean,ps=1,xr=[-10000,2000.]*11,yr=[-10000,2000.]*11,syms=.5,/xst,/yst

depict_xyz,p(0),/in,ran=range
q  = sel(p) ; { rain_gauge:p.rain_gauge, x:p.x, y:p.y, a: p.a, time: p.time }
qq = replicate( q , nreads )

iq  = 0
for n=0,nreads-1 do begin

    time_snap,f=f,p=p,nsk=n,old=olddiags

       qq[iq]= sel(p) ; { rain_gauge:p.rain_gauge, x:p.x, y:p.y, a: p.a, time: p.time }
       iq    =  iq + 1

    if not keyword_set( noplot ) then begin
       erase

       for i=0,npces-1 do begin 
           depict_xyz,p(i),noshafts=noshafts
       endfor

       if not keyword_set(runthru) then stop

       wait,0.05
    endif

endfor

qq=qq[0:iq-1]

qp=transpose( qq.qr+qq.qh+qq.qs, [2,0,1] )
qc=transpose( qq.ql+qq.qi , [2,0,1] )
qx=qc+qp
qv=transpose( qq.q , [2,0,1] )
ql=transpose( qq.ql , [2,0,1] )
qi=transpose( qq.qi , [2,0,1] )
qr=transpose( qq.qr , [2,0,1] )
qs=transpose( qq.qs , [2,0,1] )
qh=transpose( qq.qh , [2,0,1] )
qvb=transpose( qq.qbck(*,0) )

pcf=transpose(qq.rain_gauge)

th=transpose( qq.th , [2,0,1] )
tht=transpose( qq.th+qq.thbck , [2,0,1] )
thb=transpose( qq.thbck , [2,0,1] )

a=transpose( qq.a , [2,0,1] )
w=transpose( qq.w , [2,0,1] )
e=transpose( qq.e , [2,0,1] )

atot=transpose( qq.a(*,0)+qq.a(*,1)+qq.abck(*,0) )

;;qq=create_struct( qq, 'a0',a0 )

thbz=thb*0.
thbz(*,1:98,*)=(thb(*,2:99,*)-thb(*,0:97,*))

aqoo=qq.a(*,0)/p(0).a0
boo=where( aqoo gt 1.0 )
if min(boo) gt -1 then aqoo(boo)=1.

;thoo=transpose( qq.th(*,0)*qq.a(*,0)/p(0).a0 )
;woo=transpose( qq.w(*,0)*qq.a(*,0)/p(0).a0 )
thoo=transpose( qq.th(*,0)*aqoo )
woo=transpose( qq.w(*,0)*aqoo )
poo=transpose( (qq.phyd(1:100,0)+qq.pp(1:100,0))*aqoo )


r = sqrt( a / !pi )

ue= -e/(2*!pi*r)

time=qq.time(0)
z=p(0).zo
ze=p(0).z
a0=p(0).a0

;for l=1,98 do begin
;   thbz(*,l,*)=thbz(*,l,*)/(z(l+1)-z(l-1))
;endfor

;thbz(*,0,*)=thbz(*,1,*)
;thbz(*,99,*)=thbz(*,98,*)

;disp=th/thbz


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;---------------------------------------------------------------------
; FYI: The following command displays a time-height cross section of 
; updraft area for zero-th PCE
;
; IDL> contour,transpose(qq.a(*,0)),qq.time(0),p0(0).zo

print,'IDL> contour,transpose(qq.a(*,1))/qq(0).a0(1),qq.time(1),zo,lev=findgen(20)*.125'
contour,transpose(qq.a(*,1))/qq(0).a0(1),qq.time(1),zo,lev=findgen(20)*.125


end
