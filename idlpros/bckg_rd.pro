pro bckg_rd,nskip=nskip,atime=atime,bckg=bckg,fname=fname

if not keyword_set(fname) then fname='fort.111'

close,1
openr,1,/f77_u,fname

   ;WRITE( LUN )  TheTimeNow
   ;WRITE( LUN )  BCKG%IMS,BCKG%JMS,BCKG%LMS
   ;WRITE( LUN )  BCKG%AREA
   ;WRITE( LUN )  BCKG%XX
   ;WRITE( LUN )  BCKG%YY
   ;WRITE( LUN )  BCKG%ZGE
   ;WRITE( LUN )  BCKG%U
   ;WRITE( LUN )  BCKG%V
   ;WRITE( LUN )  BCKG%W
   ;WRITE( LUN )  BCKG%Q
   ;WRITE( LUN )  BCKG%TH


   READU,1,  TheTimeNow
ims=0L&jms=0L&lms=0L
   READU,1,  IMS,JMS,LMS


area = fltarr(ims,jms)
xx   = fltarr(ims,jms)
yy   = fltarr(ims,jms)
zge  = fltarr(ims,jms,lms+1)
u    = fltarr(ims,jms,lms)
v    = fltarr(ims,jms,lms)
w    = fltarr(ims,jms,lms+1)
q    = fltarr(ims,jms,lms)
th   = fltarr(ims,jms,lms)

aa   = fltarr(ims,jms,lms)
ppe  = fltarr(ims,jms,lms+1)




   READU,1,  AREA
   READU,1,  XX
   READU,1,  YY
   READU,1,  ZGE
   READU,1,  U
   READU,1,  V
   READU,1,  W
   READU,1,  Q
   READU,1,  TH
   READU,1,  AA
   READU,1,  PPE


   if nskip gt 0 then begin
 for n=1,nskip do begin

   READU,1,  TheTimeNow
   READU,1,  IMS,JMS,LMS
   READU,1,  AREA
   READU,1,  XX
   READU,1,  YY
   READU,1,  ZGE
   READU,1,  U
   READU,1,  V
   READU,1,  W
   READU,1,  Q
   READU,1,  TH
   READU,1,  AA
   READU,1,  PPE

 endfor 
   endif


close,1

bckg={ time:thetimenow, im:ims,jm:jms,lm:lms $
     , area:area,xx:xx,yy:yy,zge:zge,u:u,v:v,q:q,w:w,th:th,aa:aa,ppe:ppe}



return
end


