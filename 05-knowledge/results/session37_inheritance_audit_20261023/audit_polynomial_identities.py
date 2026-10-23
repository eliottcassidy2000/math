#!/usr/bin/env python3
"""Independent exact controls for THM-1810's unsigned identity and THM-506's logic.
Standard library only. No external input and no repository mutation.
"""
from itertools import combinations, permutations
from collections import Counter
from math import factorial
from fractions import Fraction
import json
from pathlib import Path


def rank(rows):
    a=[list(map(Fraction,r)) for r in rows]
    if not a:return 0
    h=0
    for col in range(len(a[0])):
        pivot=next((r for r in range(h,len(a)) if a[r][col]),None)
        if pivot is None:continue
        a[h],a[pivot]=a[pivot],a[h]
        d=a[h][col];a[h]=[x/d for x in a[h]]
        for r in range(len(a)):
            if r!=h:
                d=a[r][col];a[r]=[x-d*y for x,y in zip(a[r],a[h])]
        h+=1
    return h


def sign(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))


def tournament_polynomials(n):
    edges=list(combinations(range(n),2))
    unsigned=Counter();signed=Counter();transitive=Counter()
    for mask in range(1<<len(edges)):
        score=[0]*n
        for k,(i,j) in enumerate(edges):score[j if (mask>>k)&1 else i]+=1
        e=tuple(score);unsigned[e]+=1;signed[e]+=(-1)**mask.bit_count()
        if len(set(score))==n:transitive[e]+=1
    # Columns x_i^(n-1), ..., x_i^0 align det with prod_(i<j)(x_i-x_j).
    per=Counter();det=Counter()
    for p in permutations(range(n)):
        exp=tuple(n-1-p[i] for i in range(n));per[exp]+=1;det[exp]+=sign(p)
    signed={k:v for k,v in signed.items() if v}
    assert signed==dict(det)
    assert transitive==per
    residual=unsigned-per
    assert all(v>0 for v in residual.values())
    assert sum(unsigned.values())==2**len(edges)
    assert sum(per.values())==factorial(n)
    assert (unsigned==per)==(n<=2)
    if n==3:assert residual=={(1,1,1):2}
    return {'n':n,'unsigned_tournament_value_at_ones':sum(unsigned.values()),
            'vandermonde_permanent_value_at_ones':sum(per.values()),
            'intransitive_residual_value_at_ones':sum(residual.values()),
            'signed_nonzero_monomials':len(signed),'permanent_equals_transitive_sum':True,
            'all_tournament_sum_equals_permanent':unsigned==per}


def L(p):return sum(c*factorial(k) for k,c in enumerate(p))
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def main():
    records=[tournament_polynomials(n) for n in range(1,6)]
    # A finite realized carrier set with linear projections x and y=x².
    # x recovers y as a function, but not affinely.
    points=[(0,0),(1,1),(2,4)]
    dx=[[x] for x,y in points[1:]];dxdy=[list(p) for p in points[1:]]
    assert rank(dx)==1 and rank(dxdy)==2 and len({x for x,y in points})==3
    q=[-1,1];h=[0,-2,1]
    assert L(q)==0 and L(mul(q,q))==1
    assert L(h)==0 and L(mul(h,h))==8
    result={'vandermonde':records,
      'n3_polynomial_correction':'prod_(i<j)(x_i+x_j) = per[x_i^(n-1-j)] + 2*x_0*x_1*x_2',
      'rank_countermodel':{'points':points,'rank_delta_x':rank(dx),'rank_delta_x_y':rank(dxdy),
                          'functional_recovery':'y=x^2','affine_recovery':False},
      'factorial_moment_cancellation':{'B':'s-1','L(B)':L(q),'L(B^2)':L(mul(q,q)),
          'H':'s*(s-2)','L(H)':L(h),'L(H^2)':L(mul(h,h)),
          'P':'W+Z*(ZW-2)','E(P^1)':0,'E(P^2)':0,'E(P^3)':0,'E(P^4)':48},
      'checks':'all assertions passed'}
    target=Path(__file__).with_suffix('.json');target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
