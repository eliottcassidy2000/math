#!/usr/bin/env python3
"""Exhaustive normal-form tests of the monomial coefficient-precision clock.
Uses direct ring multiplication and first hitting times, not derivative bounds.
All checks remain active with Python -O.
"""
from collections import Counter
from itertools import product
from math import prod
import json

checks = Counter()

def check(ok,lane,detail=None):
    checks[lane] += 1
    if not ok:
        raise RuntimeError(f'{lane}: {detail}')

def vp(x,p):
    if x == 0:
        return 10**9
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v

def hensel(x,e,p,a):
    t=x%p
    for j in range(1,a):
        t += ((pow(t,e,p**(j+1))-t)//p**j % p)*p**j
    return t

class MonomialRing:
    def __init__(self,p,weights):
        self.p=p
        self.weights={g:w for g,w in weights.items() if w>0}
        self.basis=sorted(self.weights,key=lambda g:(sum(g),g))
        self.index={g:i for i,g in enumerate(self.basis)}
        self.moduli=[p**self.weights[g] for g in self.basis]
        self.one=(1,)+(0,)*(len(self.basis)-1)
        self.mul_table=[]
        for g in self.basis:
            row=[]
            for h in self.basis:
                t=tuple(x+y for x,y in zip(g,h))
                row.append(self.index.get(t))
            self.mul_table.append(row)
        # Every divisor must exist with at least the same coefficient precision.
        for g,w in self.weights.items():
            for h in product(*(range(x+1) for x in g)):
                check(self.weights.get(h,0)>=w,'monomial_weight_validity',(g,h,w))
    def mul(self,x,y):
        z=[0]*len(self.basis)
        for i,a in enumerate(x):
            if a:
                for j,b in enumerate(y):
                    k=self.mul_table[i][j]
                    if b and k is not None:
                        z[k]=(z[k]+a*b)%self.moduli[k]
        return tuple(z)
    def power(self,x,e):
        y=self.one
        while e:
            if e&1:
                y=self.mul(y,x)
            x=self.mul(x,x)
            e>>=1
        return y
    def clocks(self,e):
        p=self.p
        s=vp(e,p)
        a=self.weights[self.basis[0]]
        H=0
        direct=0
        for g,w in self.weights.items():
            pos=[x for x in g if x]
            if pos:
                direct=max(direct,w+min(vp(x,p) for x in pos))
                if len(pos)==1:
                    x=pos[0]
                    h=vp(x,p)
                    if x==p**h:
                        H=max(H,h+w)
        check(H==direct,'axis_reduction_vs_all_monomials')
        return (max(a-1,H)+s-1)//s,H

cases=[]
seen=set()
def run_case(p,e,weights,label,cap=7000):
    key=(p,e,tuple(sorted((g,w) for g,w in weights.items() if w>0)))
    if key in seen:
        return
    seen.add(key)
    size=p**sum(weights.values())
    if size>cap:
        return
    ring=MonomialRing(p,weights)
    expected,H=ring.clocks(e)
    a=ring.weights[ring.basis[0]]
    worst=0
    fixed=0
    for x in product(*(range(m) for m in ring.moduli)):
        target=(hensel(x[0],e,p,a),)+(0,)*(len(ring.basis)-1)
        y=x
        t=0
        while y!=target:
            y=ring.power(y,e)
            t+=1
            if t>expected:
                raise RuntimeError(f'upper bound: {label,p,e,weights,x,expected}')
        worst=max(worst,t)
        fixed += t==0
        check(ring.power(x,e**expected)==target,'exhaustive_normal_form_retraction')
    check(worst==expected,'exact_axis_clock',(label,p,e,weights,worst,expected))
    check(fixed==p,'constant_teichmuller_fixed_points')
    cases.append({'label':label,'p':p,'exponent':e,'size':size,'basis_size':len(ring.basis),
                  'a':a,'H':H,'clock':worst,
                  'weights':{','.join(map(str,g)):w for g,w in sorted(ring.weights.items())}})

# Every monotone univariate profile of length <=5 and constant precision <=3
# at p=3, subject only to the explicit cardinality cap. Exponents 3 and 9.
for a in range(1,4):
    for b in range(1,6):
        for tail in product(range(1,a+1),repeat=b-1):
            ws=(a,)+tail
            if any(ws[j]<ws[j+1] for j in range(b-1)):
                continue
            weights={(j,):w for j,w in enumerate(ws)}
            for e in [3,9]:
                run_case(3,e,weights,'all_small_univariate_profiles')

# A second residue characteristic and a non-prime-power exponent.
for a in range(1,3):
    for b in range(1,5):
        for tail in product(range(1,a+1),repeat=b-1):
            ws=(a,)+tail
            if any(ws[j]<ws[j+1] for j in range(b-1)):
                continue
            for e in [5,45]:
                run_case(5,e,{(j,):w for j,w in enumerate(ws)},'p5_univariate_profiles')

# Every 2-variable square-zero-axis profile at p=3, precision <=3, cap 7000.
# Crossing terms may survive at any legal precision or be killed entirely.
for a in range(1,4):
    for u in range(1,a+1):
        for v in range(1,a+1):
            for z in range(0,min(u,v)+1):
                ws={(0,0):a,(1,0):u,(0,1):v,(1,1):z}
                for e in [3,9]:
                    run_case(3,e,ws,'all_small_bivariate_profiles')

# Named n=561 profiles, including a nonflat example in the theorem note.
for ws in [(3,1,1,1),(3,2,1,1),(2,1,1,1,1),(3,3),(2,2,2,2)]:
    run_case(3,561,{(j,):w for j,w in enumerate(ws)},'carmichael_exponent_profile')
# A three-variable example: many mixed monomials, all coefficient precision one.
run_case(3,561,{g:1 for g in product(range(2),repeat=3)},'trivariate_flat_box')

# Boundary control: odd-prime restriction is essential for the pure arithmetic
# sharpness term when no variable survives. Exponent6 modulo8 stabilizes in 1.
ring=MonomialRing(2,{(0,):3})
formal_bound,_=ring.clocks(6)
actual=max(0 if x==hensel(x,6,2,3) else 1 for x in range(8))
for x in range(8):
    check(ring.power((x,),6)==(hensel(x,6,2,3),),'dyadic_boundary_actual')
check(actual==1 and formal_bound==2,'dyadic_sharpness_hostile')

print(json.dumps({'status':'PASS','total_checks':sum(checks.values()),
                 'checks':dict(sorted(checks.items())),
                 'case_count':len(cases),
                 'enumerated_elements':sum(x['size'] for x in cases),
                 'universe':'All profiles specified in script; cardinality cap 7000 per ring',
                 'dyadic_boundary':{'p':2,'a':3,'exponent':6,'claimed_odd_formula':2,'actual':1},
                 'cases':cases},indent=2,sort_keys=True))
