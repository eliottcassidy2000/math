#!/usr/bin/env python3
"""Exact endpoints for the mass-paid tail proof; checks remain active with -O."""
from fractions import Fraction as F
from math import factorial
import json

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def main():
    n=514;R=257;a=F(3,4);c=F(249999,1000000);d=F(2491,10000)
    t=d/8;lam=d/64;h=d/16
    need(sum(F(4**j,factorial(j)) for j in range(8))>50,'bootstrap exponential')
    need(F(49,200)>F(4,R) and F(49*(n-1),1600)>4,'variance bootstrap')
    need(F(1,4)*(1-F(1,50**4))>c,'nonlinear c floor')
    heavy=lambda n:c*F(n-1,2)-F(16*n,n-1)-d*F(n,4)-16
    need(heavy(n)==F(86393731,1026000000)>0,'integer heavy count17')
    need(c/2-d/4==F(125449,2000000)>0,'heavy monotonicity')
    Q=lambda r:F(r)+F(33,2)+F(8,r)
    big=Q(R)*a**32/(R-32)
    need(big==F(5789452743369435231,47408132269433547653120),'large-m rational')
    need(big<c-d and Q(R)*a**(R-1)<c-d,'both convex endpoints')
    small=lambda r:c*F(r,2)-8-F(8,2*r-1)-F(31,2)-Q(r)*a**17
    need(small(R)==F(463842995492257303903,70781597908992000000)>0,'small-m paid tail')
    need(c/2>a**17,'small-m monotonicity')
    need(t/(4*lam)==2,'exact mass coefficient')
    need(h-F(4,R)==F(187,41120000)>0 and t>F(4,R),'final contradiction')
    # Verify every admissible m on a substantial finite range, independent
    # of the producer's split at32. Uniformity is supplied by the written proof.
    controls=0
    for r in range(257,601):
        for m in range(17,r):
            min_u=max(F(0),c*F(r,4)-4-F(4,2*r-1)-F(m,4))
            need(Q(r)*a**m < (c-d)*(r-m)+2*min_u,'unsplit tail control')
            controls+=1
    # Hostiles: throwing away actual light mass loses this certificate;
    # the chosen constants cannot be promoted to order512.
    need(Q(R)*a**17>(c-d)*(R-17),'U-free shortcut must fail')
    need(h<F(4,256)*(1-a**256-F(1,4)**256),'n512 not certified')
    print(json.dumps({'status':'passed','threshold':514,'remaining_open_orders':[10,513],
        'heavy_endpoint_slack':str(heavy(n)),'small_tail_slack':str(small(R)),
        'large_tail_endpoint':str(big),'tail_allowance':str(c-d),
        'dual_coefficient':str(h),'final_surplus':str(h-F(4,R)),
        'unsplit_finite_controls':controls,'hostiles':['discarding U','same constants at512'],
        'scope':'Exact endpoints plus analytic monotonicity prove all n>=514; the finite controls alone do not.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
