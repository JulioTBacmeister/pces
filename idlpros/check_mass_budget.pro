pro check_mass_budget,qq=qq

s=size(qq.a)&nt=s(3)&nz=s(1)&np=s(2)

a0=qq(0).a0(0)

mass_env = total( qq(1:*).a(*,0)*qq(1:*).mass(*,1) , 1 )
mass_upd = total( qq(1:*).a(*,1)*qq(1:*).mass(*,1) , 1 )
mass_shf = total( qq(1:*).a(*,2)*qq(1:*).mass(*,1) , 1 )

mass_tot = mass_env + mass_upd + mass_shf

STOP

return
end
