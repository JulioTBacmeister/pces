function analytic_p,lambda=lb,z=z

h=max(z)

b  = 1./( 1. - exp( 2.*lb*h ) )
a  = -b * exp( 2.*lb*h )


p  = a * exp( -lb * z ) + b * exp( lb*z )

return,p
end
