pro profile_ani,qq=qq,a0=a0,z=z


a=transpose( qq.a , [2,0,1] )

time=qq.time(0)

s=size(a)

for n=0,s(1)-1 do begin

    ;;if n eq 0 then plot,a(n,*,0),z,zr=[0,10*a0]

    ;;oplot,a(n,*,0),z
    plot,a(n,*,0),z,xr=[0,10*a0],/xst,/yst

    oplot,a(n,*,1),z,colo=2

    wait,.05
   


endfor


return
end
