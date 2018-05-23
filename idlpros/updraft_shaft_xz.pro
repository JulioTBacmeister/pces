pro updraft_shaft_xz,qq=qq,itm=itm

   ; stdcolo
zo=qq(itm).zo(*,1)
    
plot,qq(itm).x(*,1)+sqrt(qq(itm).a(*,1)/!pi),zo   ,xr=[-30,30]*1000. 
oplot,qq(itm).x(*,1)-sqrt(qq(itm).a(*,1)/!pi),zo 

oplot,qq(itm).x(*,2)+sqrt(qq(itm).a(*,2)/!pi),zo,colo=2 
oplot,qq(itm).x(*,2)-sqrt(qq(itm).a(*,2)/!pi),zo,colo=2


return
end
