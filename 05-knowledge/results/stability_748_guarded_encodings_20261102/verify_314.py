#!/usr/bin/env python3
"""Exact combined-budget checks, with active hostile controls under -O."""
from fractions import Fraction as F
from itertools import product
from math import factorial
import json

def need(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    R=157;c=F(249999,1000000);d=F(32,157);t=d/8;lam=d/64;h=t-4*lam;a=F(3,4)
    Q=lambda r:F(r)+F(33,2)+F(8,r)
    g=lambda r:(c-4*t)*r-16-F(8,r)
    tail=lambda r,m:(c-d)*(r-m)+c*F(r,2)-8-F(8,2*r-1)-F(m,2)-Q(r)*a**m
    need(sum(F(4**j,factorial(j)) for j in range(8))>50,'exp4 bootstrap')
    need(F(49,200)>t and F(49*313,1600)>4,'variance bootstrap')
    need(F(1,4)*(1-F(1,50**4))>c,'uniform c')
    need(1-4*t==F(141,157)>0,'positive slot divisor')
    need(g(R)-8*(1-4*t)==F(2225351,157000000)>0,'heavy count9')
    need(c-4*t>0,'heavy count monotone')
    need(c-d==F(7249843,157000000)>0,'tail budget')
    large=Q(R)*a**13/(R-13)
    need(large<c-d and Q(R)*a**(R-1)<c-d,'convex tail endpoints')
    need(tail(R,9)==F(363349507136067,402563072000000)>0,'concave endpoint9')
    need(tail(R,12)==F(175134030288958737,25764036608000000)>0,'concave endpoint12')
    need(F(3,2)*c-d-a**9>0,'small-tail uniform monotonicity')
    need(t==F(4,157) and lam==F(1,314) and h==F(2,157),'closed final constants')
    # Both expressions are multiaffine in these eight independent symbols;
    # equality on all Boolean vertices checks every coefficient exactly.
    for k,V,B,r,L,e,D,U in product([F(0),F(1)],repeat=8):
        before=(k*V-B)*(r-L)-4*(1-e)*(V-D-U)
        after=(k*r-4+4*e)*V+4*(1-e)*D-B*r-L*(k*V-B)+4*(1-e)*U
        need(before==after,'retained -BR identity')
    # The direct quadratic-minimum audit from the514 argument, now at d=32/157.
    for i in range(101):
        w=t*F(i,100);S=2*t-w
        residual=w*S/2+S*S/4-t*S+lam*(d-4*w)-t*w/4
        need(residual==w*(t-w)/4>=0,'positive U credit')
    controls=0
    for r in range(157,501):
        for m in range(9,r):
            min_u=max(F(0),c*F(r,4)-4-F(4,2*r-1)-F(m,4))
            need(Q(r)*a**m<(c-d)*(r-m)+2*min_u,'unsplit combined tail')
            controls+=1
        if lam*r<1:need((1-lam*r)*F(4,r)<=h,'uniform final contradiction')
    # At312 the selected constants do not justify the required count; and
    # the scalar m8 cell violates the tail bound. This is not a tournament.
    need(g(156)<8*(1-4*t),'312 count hostile')
    need(tail(156,8)<0,'312 relaxed tail hostile')
    # Both new credits matter: without -BR, the final coefficient h is too
    # small; without 2U the m9 geometric tail cannot be absorbed.
    need(h<F(4,R)*(1-a**R),'baseline omission hostile')
    need(Q(R)*a**9>(c-d)*(R-9),'mass omission hostile')
    print(json.dumps({'status':'passed','threshold':314,'remaining_open_orders':[10,313],
        'source_commits':['291b4d975784','9c7eb0b13e','e9dd62edb6'],
        'slot_endpoint_slack':str(g(R)-8*(1-4*t)),
        'small_tail_slacks':[str(tail(R,9)),str(tail(R,12))],
        'large_tail_endpoint':str(large),'tail_allowance':str(c-d),
        'baseline_polynomial_vertices':256,'dual_minima':101,
        'unsplit_finite_controls':controls,
        'hostiles':['same constants at312','omit -BR','omit actual U'],
        'scope':'Written monotonic/convex proof establishes all n>=314; finite controls corroborate it. No independent-agent audit claimed.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
