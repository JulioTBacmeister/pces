;print,' give a time '
;read, ' ... ', tink


;oo=min( abs(time-tink), ioo )

 diag1=transpose(qq.rhol(*,0)*qq.w(*,0))
 diag2=diag1*0
  
  for l=1,98 do diag2(*,l)=(diag1(*,l+1)-diag1(*,l-1))/(zo(l+1)-zo(l-1))

  diag2 = a(*,*,0)*diag2



  diag3 = transpose( qq.rhol(*,0)*qq.e(*,0) )
  oss=where( a(*,*,0)/a0 gt .75  )

  plot,diag3(oss)/diag2(oss),ps=2,yr=[0,1.5],/yst


end
