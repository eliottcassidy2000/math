#!/usr/bin/env python3
"""Numerical algebra smoke tests, not checks of zeta zeros or RH."""
import cmath, math, json

def require(ok):
 if not ok: raise RuntimeError('Smoke test failed')
b=10.; rho=.5+14.25j
lam=cmath.exp((rho-1)*math.log(b)); c=-(cmath.exp(rho*math.log(b))-1)/((b-1)*rho)
errors=[]
for m in range(1,8):
 direct=-(cmath.exp(rho*(m+1)*math.log(b))-cmath.exp(rho*m*math.log(b)))/(rho*(b-1)*b**m)
 errors.append(abs(direct-c*lam**m))
require(max(errors)<1e-12)
rho_alias=rho+2j*math.pi/math.log(b)
require(abs(cmath.exp((rho_alias-1)*math.log(b))-lam)<1e-12)
require((1/rho+1/rho_alias).real>0)
z=.2+.1j
partial=sum(c*lam**m*z**m for m in range(1,60)); rational=c*lam*z/(1-lam*z)
require(abs(partial-rational)<1e-12)
trivial=[]
for m in range(1,5):
 direct=-(math.log1p(-b**(-2*m-2))-math.log1p(-b**(-2*m)))/(2*(b-1)*b**m)
 expanded=sum((b**(-2*j)-1)/(2*j*(b-1))*b**(-(2*j+1)*m) for j in range(1,30))
 trivial.append(abs(direct-expanded));require(abs(direct-expanded)<1e-18)
print(json.dumps({'status':'Numerical algebra smoke tests only; synthetic rho is not asserted to be a zero','max_single_mode_error':max(errors),'geometric_identity_error':abs(partial-rational),'trivial_correction_errors':trivial},indent=2,sort_keys=True))
