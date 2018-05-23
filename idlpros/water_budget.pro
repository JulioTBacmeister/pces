pro water_budget,qq=qq0

qq=qq0(5:*)


dz=transpose( qq.ze(0:99,0)-qq.ze(1:100)  )

z=(qq(0).ze(0:99)+qq(0).ze(1:100))/2.
time=qq.time(0)

mass=transpose( qq.mass(*,0) )

aru =transpose( qq.a(*,0) )

ars =transpose( qq.a(*,1) )

arb =transpose( qq.abck(*,0) )

rhau=transpose( qq.rhoa(*,0) )

rhas=transpose( qq.rhoa(*,1) )



qtu = transpose(  qq.q(*,0) + qq.ql(*,0) + qq.qi(*,0) + qq.qr(*,0) + qq.qs(*,0) + qq.qh(*,0) )
qtb = transpose(  qq.qbck(*,0) + qq.qlbck(*,0) + qq.qibck(*,0) + qq.qrbck(*,0) + qq.qsbck(*,0) + qq.qhbck(*,0) )
qts = transpose(  qq.q(*,1) + qq.ql(*,1) + qq.qi(*,1) + qq.qr(*,1) + qq.qs(*,1) + qq.qh(*,1) )

qvu = transpose(  qq.q(*,0) )
qvb = transpose(  qq.qbck(*,0) )
qvs = transpose(  qq.q(*,1) )

qcu = transpose(  qq.ql(*,0) + qq.qi(*,0)  )
qcb = transpose(  qq.qlbck(*,0) + qq.qibck(*,0) )
qcs = transpose(  qq.ql(*,1) + qq.qi(*,1) )

qpu = transpose(  qq.qr(*,0) + qq.qs(*,0) + qq.qh(*,0) )
qpb = transpose(  qq.qrbck(*,0) + qq.qsbck(*,0) + qq.qhbck(*,0) )
qps = transpose(  qq.qr(*,1) + qq.qs(*,1) + qq.qh(*,1) )


mqtu1 = total (  aru  * mass * qtu , 2 )
mqtu2 = total (  rhau * dz   * qtu , 2 )

mqts1 = total (  ars  * mass * qts , 2 )
mqts2 = total (  rhas * dz   * qts , 2 )

mqtb1 = total (  arb  * mass * qtb , 2 )


mqvu1 = total (  aru  * mass * qvu , 2 )
mqvs1 = total (  ars  * mass * qvs , 2 )
mqvb1 = total (  arb  * mass * qvb , 2 )

mqcu1 = total (  aru  * mass * qcu , 2 )
mqcs1 = total (  ars  * mass * qcs , 2 )
mqcb1 = total (  arb  * mass * qcb , 2 )

mqpu1 = total (  aru  * mass * qpu , 2 )
mqps1 = total (  ars  * mass * qps , 2 )
mqpb1 = total (  arb  * mass * qpb , 2 )


mqtt1 = mqtu1 + mqts1 + mqtb1
;mqtt1(0)=mqtt1(1)
;mqvb1(0)=mqvb1(1)
n=n_elements(mqtt1)
dmqtt = (mqtt1(0:n-2)-mqtt1(1:n-1))/(time(1)-time(0))

dmqvb = (mqvb1(0:n-2)-mqvb1(1:n-1))/(time(1)-time(0))
dmqvu = (mqvu1(0:n-2)-mqvu1(1:n-1))/(time(1)-time(0))
dmqvs = (mqvs1(0:n-2)-mqvs1(1:n-1))/(time(1)-time(0))

dmqcb = (mqcb1(0:n-2)-mqcb1(1:n-1))/(time(1)-time(0))
dmqcu = (mqcu1(0:n-2)-mqcu1(1:n-1))/(time(1)-time(0))
dmqcs = (mqcs1(0:n-2)-mqcs1(1:n-1))/(time(1)-time(0))

dmqpb = (mqpb1(0:n-2)-mqpb1(1:n-1))/(time(1)-time(0))
dmqpu = (mqpu1(0:n-2)-mqpu1(1:n-1))/(time(1)-time(0))
dmqps = (mqps1(0:n-2)-mqps1(1:n-1))/(time(1)-time(0))

dmqpct=dmqpu+dmqps+dmqpb + dmqcu+dmqcs+dmqcb
dmqvt =dmqvu+dmqvs+dmqvb 

prec  = qq.rain_gauge(0)+qq.rain_gauge(1)


STOP


return
end
