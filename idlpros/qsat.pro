function qsat, t , press, celsius=celsius , kelvin=kelvin
;!!!!! PRESSURES IN MILLIBARS !!!!!!!!

AIRMW  = 28.97
H2OMW  = 18.01

EPSILON = H2OMW/AIRMW

T0=T

if max([t]) gt 30. and not keyword_set(celsius) then begin
;   print," ASSUMING T IS IN KELVIN ... SUBTRACTING 273.16 "  
   T0=T-273.16
endif

esat = 6.112*exp( 17.67*T0 / (T0+243.5) )
qsat = EPSILON*( esat / (press - (1.-EPSILON)*esat) )

return, qsat
end
