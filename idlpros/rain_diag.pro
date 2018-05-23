rho = grho(p)

epsl=0.1
itm=500
a=(qq(itm).a(*,0))
q=(qq(itm).q(*,0)) + (qq(itm).ql(*,0))
qbk=(qq(itm).qbck(*,0))

w=qq(itm).w(*,0)
s=size( p.z )&nl=s(1)
dz=p(0).z(0:nl-2)-p(0).z(1:nl-1)

intgl=0.
for l=0,nl-3 do begin

    intgl=intgl + $
    epsl*2.*sqrt(!pi) * rho(l) * qbk(l) * w(l)* sqrt(a(l))*dz(l)

endfor

sflx = rho(nl-3)*a(nl-3)*q(nl-3)*w(nl-3)
qxflx=rho*a *q *w

dtm=qq(itm+1).time(0)-qq(itm).time(0)
dqvu=   total( qq(itm).q(*,0) * qq(itm).a(*,0) * rho * dz )
dqvu= ( total( qq(itm+1).q(*,0) * qq(itm+1).a(*,0) * rho * dz ) - dqvu )/dtm
dqlu=   total( qq(itm).ql(*,0) * qq(itm).a(*,0) * rho * dz )
dqlu= ( total( qq(itm+1).ql(*,0) * qq(itm+1).a(*,0) * rho * dz ) - dqlu )/dtm
dqru=   total( qq(itm).qr(*,0) * qq(itm).a(*,0) * rho * dz )
dqru= ( total( qq(itm+1).qr(*,0) * qq(itm+1).a(*,0) * rho * dz ) - dqru )/dtm


dqvs=   total( qq(itm).q(*,1) * qq(itm).a(*,1) * rho * dz )
dqvs= ( total( qq(itm+1).q(*,1) * qq(itm+1).a(*,1) * rho * dz ) - dqvs )/dtm
dqls=   total( qq(itm).ql(*,1) * qq(itm).a(*,1) * rho * dz )
dqls= ( total( qq(itm+1).ql(*,1) * qq(itm+1).a(*,1) * rho * dz ) - dqls )/dtm
dqrs=   total( qq(itm).qr(*,1) * qq(itm).a(*,1) * rho * dz )
dqrs= ( total( qq(itm+1).qr(*,1) * qq(itm+1).a(*,1) * rho * dz ) - dqrs )/dtm



t1=intgl+sflx

t2=total( qq(itm).rain_gauge )

t3=t2+dqvu+dqlu+dqru+0*dqvs+dqls+dqrs


print," Rain totals "
print,t1
print,t2
print,t3

print," Relative error"

print,2*(t3-t1)/(t3+t1)

end