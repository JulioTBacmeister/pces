pro pp_phyd,dph=dph,pp=p,z=z,a=ar,pm=pm,ps=ps

lp=n_elements(ar)
p=fltarr(lp+1)
pm=fltarr(lp+1,lp+1)

     for L=1,LP-1 do begin
         ph0  = double( dph(L) )
         zc   = double( (Z(L)+Z(L-1))*0.5 )
         arl  = double( AR(L-1) )
         rad0 = double( sqrt( ARL / !PI ) ) 
         alpha0 = 2.40482555769577d0/rad0

         for LL=0,LP-1 do begin
            if ( z(LL) lt zc ) then begin
               p1d = 0.5d*ph0*exp(-alpha0*zc) * ( exp(-alpha0*z(LL) ) +  exp( alpha0*z(LL) ) ) - ph0
            endif else begin
               p1d = 0.5d*ph0*( exp(-alpha0*zc) -  exp( alpha0*zc) ) * exp( -alpha0*z(LL) )
            endelse
            wll   = MIN( [arl/(1.0d0*AR(LL)) , 1.0d0] )
            p(LL) = p(LL) + wll * p1d
            pm(LL,L) =  wll * p1d
         ENDFOR
      endfor
      p(LP)=p(LP-1)


      ps=pm
  
      for l=1,lp do begin
        ps(*,L)=ps(*,L-1)+pm(*,L)
      endfor


return
end