pro calc_byc,qq=qq,by1=by1,by2=by2


    physparms,radius_e=radius_e  $
             ,grav=grav  $
             ,cpd=cp $
             ,rkap=rkap $
             ,ALHL=alhl $
             ,alhs=alhs $
             ,alhf=alhf $
             ,rgas=rgas





s=size(qq.a)&nt=s(3)&nz=s(1)&np=s(2)

epsi= 0.61
;epsi= 0.0

th00= reform( qq(0).thbck(nz-1,0) )

be1 = qq.th(*,1)/qq.thbck(*,0) + epsi *( qq.q(*,1)-qq.q(*,0) ) $      
   - ( qq.ql(*,1)-qq.ql(*,0) ) $
   - ( qq.qi(*,1)-qq.qi(*,0) ) $
   - ( qq.qr(*,1)-qq.qr(*,0) ) $
   - ( qq.qs(*,1)-qq.qs(*,0) ) $
   - ( qq.qh(*,1)-qq.qh(*,0) ) 

be2 = qq.th(*,2)/qq.thbck(*,0) + epsi *( qq.q(*,2)-qq.q(*,0) ) $      
   - ( qq.ql(*,2)-qq.ql(*,0) ) $
   - ( qq.qi(*,2)-qq.qi(*,0) ) $
   - ( qq.qr(*,2)-qq.qr(*,0) ) $
   - ( qq.qs(*,2)-qq.qs(*,0) ) $
   - ( qq.qh(*,2)-qq.qh(*,0) ) 

by1 = { tot:be1 }
by2 = { tot:be2 }

;; Theta 
be1 = qq.th(*,1)/qq.thbck(*,0)

be2 = qq.th(*,2)/qq.thbck(*,0) 

by1 = create_struct( by1, 'th' , be1 )
by2 = create_struct( by2, 'th' , be2 )

;; QV
be1 = epsi *( qq.q(*,1)-qq.q(*,0) )

be2 = epsi *( qq.q(*,2)-qq.q(*,0) ) 

by1 = create_struct( by1, 'qv' , be1 )
by2 = create_struct( by2, 'qv' , be2 )

;; QC
be1 = $      
   - ( qq.ql(*,1)-qq.ql(*,0) ) $
   - ( qq.qi(*,1)-qq.qi(*,0) ) $
   - ( qq.qr(*,1)-qq.qr(*,0) ) $
   - ( qq.qs(*,1)-qq.qs(*,0) ) $
   - ( qq.qh(*,1)-qq.qh(*,0) ) 

be2 = $      
   - ( qq.ql(*,2)-qq.ql(*,0) ) $
   - ( qq.qi(*,2)-qq.qi(*,0) ) $
   - ( qq.qr(*,2)-qq.qr(*,0) ) $
   - ( qq.qs(*,2)-qq.qs(*,0) ) $
   - ( qq.qh(*,2)-qq.qh(*,0) ) 



by1 = create_struct( by1, 'qc' , be1 )
by2 = create_struct( by2, 'qc' , be2 )



return
end
