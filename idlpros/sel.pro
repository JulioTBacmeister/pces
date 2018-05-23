function sel,p

q ={ rain_gauge:p.rain_gauge, i:p.i, j:p.j, x:p.x, y:p.y, a: p.a, time: p.time   $
    , type_of_ce:p.type_of_ce , my_prshft:p.my_prshft , my_updrft:p.my_updrft   $
    , age:p.age, status:p.status,maxw00:p.maxw00,maxwtf:p.maxwtf,ash:p.ash,shtop:p.shtop     $
    , thbck:p.thbck, qbck:p.qbck, abck:p.abck, ubck:p.ubck, vbck:p.vbck, wbck:p.wbck, q:p.q        $
    , ql:p.ql, qi:p.qi, qr: p.qr, qs:p.qs, qh:p.qh        $
    , qlbck:p.qlbck, qibck:p.qibck, qrbck: p.qrbck, qsbck:p.qsbck, qhbck:p.qhbck        $
    , xmean:p.xmean, ymean:p.ymean,th:p.th,e:p.e,mutt:p.mutt,tmuph:p.tmuph,qmuph:p.qmuph     $ 
    , pracs:p.pracs, rhoa:p.rhoa, rhol:p.rhol, mass:p.mass , pp:p.pp, phyd:p.phyd $ 
    , rfall:p.rfall ,  sfall:p.sfall ,  hfall:p.hfall , rcdot:p.rcdot $
    , dhn:p.dhn, ze:p.z , zo:p.zo , pe:p.pe, qtbken:p.qtbken, qtbkde:p.qtbkde $
    , shextop:p.shextop, shexbot:p.shexbot,qvsl:p.qvsl,qvsi:p.qvsi $
    , w:p.w, u:p.u , v:p.v , ur:p.ur, rad: sqrt(p.a/!pi) , a0:p.a0}


if tagset(p,'tke') then q=create_struct(q,'tke',p.tke)

return,q
end
