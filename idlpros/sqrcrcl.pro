pro sqrcrcl,rct=rct,crc=crc

;rct=[ys,xw,yn,xe]
;crc=[xc,yc,rad]


ys=rct(0)&xw=rct(1)&yn=rct(2)&xe=rct(3)
xc=crc(0)&yc=crc(1)&rad=crc(2)

x0=(xw+xe)/2. & y0=(ys+yn)/2.0

; 4 distances


dnw = sqrt ((xc-xw)^2 + (yc-yn)^2 )
dsw = sqrt ((xc-xw)^2 + (yc-ys)^2 )
dne = sqrt ((xc-xe)^2 + (yc-yn)^2 )
dse = sqrt ((xc-xe)^2 + (yc-ys)^2 )

; cyclonic (NH) ordering
vec=[ dsw, dse , dne, dnw ]


qin = vec le rad
cin = (xc ge xw) and (xc le xe) and (yc le yn) and (yc ge ys)

if cin eq 0 then begin ; circle center is NOT inside rectangle

case total(vec) of 

0: begin ; all corners outside
   over=0.
end

1: begin ; One corner inside
end

2: begin ; Two

3: begin
end

4: begin
   over=rcarea
end


endif 



STOP
return
end
