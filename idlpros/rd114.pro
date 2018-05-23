
close,1
openr,1,'../fort.114',/f77_u
lm=100
byc=fltarr(lm)
pfz=fltarr(lm)
ss=fltarr(lm)
pl=fltarr(lm)
trm3=fltarr(lm)
time=0.
for i=0,ntimes do begin
readu,1,time
readu,1,byc
readu,1,pfz
readu,1,ss
readu,1,pl
readu,1,trm3
endfor

close,1
end
