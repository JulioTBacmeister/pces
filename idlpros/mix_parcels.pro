pro mix_parcels,parcl=parcl,env=env,efrac=efrac,nitr=nitr,tv=tv,qlv=qlv,qv=qv

; EXAMPLE CALL 
;
;   mix_parcels,parcl={t:300.,q:0.03},env={t:299,q:0.01},efrac=0.09,nitr=100,tv=tv,qlv=qlv,qv=qv
;

common misc_parms, alhl,alhs,cp,rgas,airmw,h2omw,grav,rkap,EARTH_RADIUs $
                 ,tice
misc_parms


t0=parcl.t
q0=parcl.q
ql0=q0 - qsat( t0, 1000. )


te=env.t
qe=env.q


tv=[t0]&qv=[q0]&qlv=[ql0]

tx = efrac*te + (1.-efrac)*t0
qx = efrac*qe + (1.-efrac)*q0


for i=0,nitr do begin

    qlx = qx - qsat( tx , 1000. )
    dqst  = qsat( tx+1.,1000. ) - qsat(tx, 1000.)

    qlx = max( [qlx , 0. ] )

    dql = qlx - ql0
  
    tx  = tx + (1.-efrac)* dql * alhl/cp /(1.+dqst*alhl/cp ) 
    ;;tx  = tx + dql * alhl/cp /(1.+dqst*alhl/cp ) 

    t0  = tx
    ql0 = qlx
    q0  = qx


    tv = [ tv, t0 ]
    qlv= [ qlv , ql0 ]
    qv= [ qv , q0 ]

endfor





return
end
