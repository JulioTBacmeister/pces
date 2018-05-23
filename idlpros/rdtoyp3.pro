


;********************* 
        

nx=0l
ny=0l
nz=0l
;**********************************************************
  close,1
  openr,1,/f77,'../fort.211'
  readu,1,nx,ny,nz

  x    = dblarr(nx  )
  y    = dblarr(ny  )
  z    = dblarr(nz  )
  ze   = dblarr(nz+1  )

  readu,1,x,y,z,ze

  b    = dblarr(nx,ny,nz  )
  p    = dblarr(nx,ny,nz+1)
  ph   = dblarr(nx,ny,nz+1)
  d2ph = dblarr(nx,ny,nz+1)
  p1d  = dblarr(nx,ny,nz+1)

  readu,1,b
  readu,1,p
  readu,1,ph
  readu,1,d2ph
  readu,1,p1d
 
 close,1

end



