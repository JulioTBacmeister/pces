pro check_water_budget,qq=qq

s=size(qq.a)&nt=s(3)&nz=s(1)&np=s(2)

a0=qq(0).a0(0)

qtb  = qq.qbck + qq.qrbck + qq.qhbck + qq.qsbck + qq.qibck + qq.qlbck 
qtbm = total( qq.abck * qq.mass * qtb , 1 )

qtc  = qq.q + qq.qr + qq.qh + qq.qs + qq.qi + qq.ql 
qtcm = total( qq.a * qq.mass * qtc , 1 )
qlcm = total( qq.a * qq.mass * qq.ql, 1 )
qrcm = total( qq.a * qq.mass * qq.qr, 1 )
qvcm = total( qq.a * qq.mass * qq.q , 1 )

ymx=max( qtbm(0,*)+qtcm(0,*)+qtcm(1,*)+qq.rain_gauge(0)+qq.rain_gauge(1) )
yrn=[.8,1.2]*ymx

plot,qtbm(0,*)+qtcm(0,*)+qtcm(1,*)+qq.rain_gauge(0)+qq.rain_gauge(1),yr=yrn,/yst

oplot,qtbm(0,*),line=2

qbudg= qtbm(0,*)+qtcm(0,*)+qtcm(1,*)+qq.rain_gauge(0)+qq.rain_gauge(1)

STOP

return
end
