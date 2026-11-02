#!/usr/bin/env python3
"""Exact rational endpoints and symbolic constants of the n>=748 argument."""
from fractions import Fraction as F
from math import factorial
import json

def need(b,msg):
    if not b:raise RuntimeError(msg)
def kappa(n):
    r=n//2
    return F(4,r)*(1-F(3,4)**r-(F(1,4)**r if n%2==0 else 0))
def main():
    n=748;c=F(2499,10000);eta=F(1,55);C=F(181,8);cap=F(2,187)
    exp4=sum(F(4**j,factorial(j)) for j in range(8));need(exp4>50,'exponential Taylor bound')
    need(F(49,200)>2*cap,'V>2 bootstrap')
    heavy=lambda n:c*F(n-1,8)-F(4*n,n-1)-F(n,110)-F(n+2,60)
    need(heavy(n)==F(1736491,59760000)>0,'heavy endpoint')
    derivative=c/8-F(1,110)-F(1,60)
    need(derivative==F(14467,2640000)>0,'uniform increasing heavy residual')
    need(n-1>=224,'variance upper bound')
    need(16000*n*3**50<7*4**50,'uniform exponential endpoint')
    need(F(1,n)<F(1,60),'n epsilon decreases from log(4/3)>1/4')
    light=c-4*eta-F(1,4000)-4/C
    need(light==F(5063,39820000)>0,'light count bracket')
    energy=eta-C*eta*eta;need(energy==F(259,24200),'quadratic minimum')
    need(energy-cap==F(3,411400)>0 and eta>cap,'strict final contradiction')
    # Independent algebraic square completion, at every rational coefficient.
    xstar=1-2*C*eta
    need(0<xstar<1,'minimum is in the actual split interval')
    for x in [F(i,100) for i in range(101)]+[xstar]:
        need(eta*x+(1-x)**2/(4*C)==energy+(x-xstar)**2/(4*C),'additive energy identity')
    # Finite controls corroborate the uniform proof but are not its basis.
    for t in range(748,2001):
        need(heavy(t)>0 and kappa(t)<cap and t//2>=374,'finite threshold controls')
    lower=F(141421356,10**8);upper=F(141421357,10**8)
    need(lower*lower<2<upper*upper,'rational sqrt2 enclosure')
    ceiling_low=(3-2*upper)/16;ceiling_high=(3-2*lower)/16
    need(kappa(746)>ceiling_high and F(4,374)<ceiling_low,'ideal two-bin method boundary')
    need(1055**2-8*373**2==-7,'near-threshold Pell residue')
    print(json.dumps({'status':'passed','threshold':748,'new_orders':[748,1499],
        'remaining_open_orders':[10,747],'heavy_endpoint':str(heavy(n)),
        'heavy_derivative_floor':str(derivative),'light_bracket_surplus':str(light),
        'energy_floor':str(energy),'kappa_cap':str(cap),'contradiction_surplus':str(energy-cap),
        'finite_uniform_controls':1253,'ideal_method_ceiling_interval':[str(ceiling_low),str(ceiling_high)],
        'scope':'All-n extension comes from the written monotonicity and matching proofs; finite controls alone do not prove it. No all-order spectrum conclusion.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
