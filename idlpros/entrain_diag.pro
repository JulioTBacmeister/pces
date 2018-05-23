;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;*APRIL 12 2012 
;-------------------------------------------------
; This routine makes a diagnostic test for a 
; prescribed dynamics plume at steady-state (see notes
; form April 11,2012). Plume is assumed to have a constant
; vertical velocity W_0 upto a given depth and a rule for
; entrainment like U_ent = epsilon * W. Diagnostic tests shape
; of plume vs entrainment like
;
;          d(rho*a)/dz = C * rho * sqrt(a) 
;  
; where C=2*epsilon*sqrt(PI)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;*APRIL 13 2012 (~~~~ Friday the 13th ~~~~~ oooohhh!!!! ~~~~~)
;---------------------------------------------------
; More tests using prescribed dynamics plume. Shut off uphys. 
; (w/ #define SKIPUPH in ce1.F90) 
; So q is a good conservative tracer.  Under these circumstances
; you should get a relation
;
;   d(q)/dz = C * ( q_bck - q )/sqrt(a) 
;
; where C is same as before.  In fact you do. See fig in pces_dev.docx
; (current path shiya::/Users/juliob/Documents/PCES/)
;
;   IDL> plot,dq/( (qbk-q)/sqrt(a) ),yr=[0,1]
;   IDL> oplot,drha/(rho*sqrt(a)),ps=-4      
;
; This effectively overplots two estimates of coefficient "C". Runs in plot
; used epsilon=0.1, so C=0.3545

rho = grho(p)

itm=300
a=(qq(itm).a(*,0))
q=(qq(itm).q(*,0)) + (qq(itm).ql(*,0))
qbk=(qq(itm).qbck(*,0))

rha=rho*a

drha=rha*0
drha(1:98)=(rha(2:99)-rha(0:97))/(zo(2:99)-zo(0:97))
dq=q*0
dq(1:98)=(q(2:99)-q(0:97))/(zo(2:99)-zo(0:97))


drha(0:50)=0.     

;;plot,drha/rh*sqrt(a)              

;;cee=epsilon*2.*sqrt(!pi)


end