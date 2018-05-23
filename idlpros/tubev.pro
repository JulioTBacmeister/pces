function tubev,qq,theta=theta,w=w,press=press

aqoo=qq.a(*,0)/qq(0).a0(0)
boo=where( aqoo gt 1.0 )
if min(boo) gt -1 then aqoo(boo)=1.

;thoo=transpose( qq.th(*,0)*qq.a(*,0)/p(0).a0 )
;woo=transpose( qq.w(*,0)*qq.a(*,0)/p(0).a0 )
thoo=transpose( qq.th(*,0)*aqoo )
woo=transpose( qq.w(*,0)*aqoo )
poo=transpose( (qq.phyd(1:100,0)+qq.pp(1:100,0))*aqoo )



if keyword_set(theta) then foo=thoo
if keyword_set(w) then foo=woo

return,foo
end
