function grho,p
grav=9.81


pe=p.pe * 100. ; convert from hPa to Pa
ze=p.z

s=size(pe)
;;STOP
case s(0) of 
1: begin
   rho=( pe(0:s(1)-2)-pe(1:s(1)-1) )/( ze(0:s(1)-2)-ze(1:s(1)-1) )
end
2: begin
   rho=( pe(0:s(1)-2,*)-pe(1:s(1)-1,*) )/( ze(0:s(1)-2,*)-ze(1:s(1)-1,*) )
end
3: begin
   rho=( pe(0:s(1)-2,*,*)-pe(1:s(1)-1,*,*) )/( ze(0:s(1)-2,*,*)-ze(1:s(1)-1,*,*) )
end

endcase
rho=-rho / grav




return,rho
end