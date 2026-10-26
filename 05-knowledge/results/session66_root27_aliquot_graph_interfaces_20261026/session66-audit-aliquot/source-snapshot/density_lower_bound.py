#!/usr/bin/env python3
"""Certified Chen--Zhao Theorem1 evaluation, using only integer arithmetic.
For each nonnegative rational summand t, retain floor(10^24 t).
Sum/D <= g_M < (sum+max(1,number_of_positive_terms))/D.
This is a lower-density bound, not a natural-density existence result.
"""
import json
from pathlib import Path
from math import prod
D=10**24
BASE=[(2,6),(3,5),(5,4),(7,3),(11,2),(13,1),(17,1),(19,1),(23,1),(29,1),(31,1),(37,1),(41,1)]

def bound(factorization):
    # State d, sigma(2d), numerator/denominator of phi(M/d)/(M/d).
    states=[(1,1 if any(p==2 for p,e in factorization) else 3,1,1)]
    for p,e in factorization:
        nxt=[]
        for d,s,w,z in states:
            for k in range(e+1):
                # 2d exponent for p=2 is k+1, not k.
                ek=k+(p==2)
                sd=(p**(ek+1)-1)//(p-1)
                nw,nz=(w*(p-1),z*p) if k<e else (w,z)
                nxt.append((d*p**k,s*sd,nw,nz))
        states=nxt
    total=positive=0
    for d,s,w,z in states:
        sm=s-2*d
        if sm>2*d:
            # w/z * [1/(2d)-1/s(2d)]
            num=w*(sm-2*d); den=z*2*d*sm
            total+=(D*num)//den
            positive+=1
    return {'M_factorization':factorization,'M':prod(p**e for p,e in factorization),
        'divisor_count':len(states),'positive_terms':positive,
        'common_denominator':D,'lower_numerator':total,'strict_upper_numerator':total+max(1,positive),
        'lower_decimal':f'{total//D}.{total%D:024d}',
        'strict_upper_decimal':f'{(total+max(1,positive))//D}.{(total+max(1,positive))%D:024d}',
        'claim':'liminf U(X)/X >= g_M >= lower_numerator/common_denominator'}

if __name__=='__main__':
    zero_control=bound([])
    if zero_control['lower_numerator']!=0 or zero_control['strict_upper_numerator']!=1:
        raise RuntimeError('M=1 exact-zero enclosure guard')
    odd_control=bound([(3,2)])
    if not (odd_control['lower_numerator']*126 <= D < odd_control['strict_upper_numerator']*126):
        raise RuntimeError('odd M=9 control must enclose 1/126')
    out={'method':'Chen-Zhao2011 Theorem1, directed integer rounding of each rational term',
        'base':bound(BASE),'enlarged':bound(BASE+[(43,1)])}
    if out['base']['lower_numerator']*10_000_000 <= 602757*D: raise RuntimeError('published bound not reproduced')
    if out['enlarged']['lower_numerator']<=out['base']['strict_upper_numerator']: raise RuntimeError('no certified gain')
    path=Path(__file__).resolve().parents[1]/'results'/'density_lower_bound.json'
    path.write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
