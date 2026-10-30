#!/usr/bin/env python3
"""Exact E12 obstruction-clock bands; no E12 family closure is claimed."""
from fractions import Fraction as F
import json


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def floor(q):
    return q.numerator // q.denominator


def ceil(q):
    return -floor(-q)


def residue_interval(lo, hi, r):
    # Strict lower, weak upper endpoint.
    first = r + 8 * (floor((lo-r)/8)+1)
    last = r + 8 * floor((hi-r)/8)
    return None if first > last else (first, last)


def intervals(b):
    out=[]
    for r in range(8):
        f=F(5*r % 8,8)
        for name,lo,hi in [('early',F(0),D*(f-b*gamma)),
                           ('late',5*D*f,D*(1+f-b*gamma))]:
            iv=residue_interval(lo,hi,r)
            if iv:
                out.append(dict(residue=r,band=name,first=iv[0],last=iv[1],count=(iv[1]-iv[0])//8+1))
    return out


gamma=F(5878847,19531250)
lam11=F(12207031,19531250)
lam12=F(5,8)*(1-F(1,5**12))
Delta=lam12-lam11
D=8*5**10
T=2*5**10
need(Delta==F(1,2*5**11), 'slope gap')
need(lam11==F(5,8)-F(1,D), 'eighth decomposition')
need(lam12==F(5,8)-F(1,5*D), 'next eighth decomposition')
families=[(a,b) for b in range(1,4) for a in range(5) if F(a,5)+b*gamma<1]
need(families==[(0,1),(1,1),(2,1),(3,1),(0,2),(1,2),(0,3)], 'families')

# Independent finite exhaustive proof of the rational floor identity over
# every possible integer clock: subtract 5k+floor(5r/8), then partition each
# residue progression at every possible floor jump of both affine residuals.
# This visits all cells, NOT merely a finite selection of clock samples.
cells=0
for b in (1,2,3):
    upper=ceil((1-b*gamma)/Delta)-1
    bands=intervals(b)
    for r in range(8):
        klo=max(0,ceil(F(1-r,8)))
        khi=floor(F(upper-r,8))
        if klo>khi:
            continue
        cuts={klo,khi+1}
        f=F(5*r % 8,8)
        # Residual is intercept - rate*k; exact-floor changes just after
        # it hits integer j, i.e. at floor((intercept-j)/rate)+1.
        for intercept,rate in [(f-F(r,D)-b*gamma,F(8,D)),
                               (f-F(r,5*D),F(8,5*D))]:
            low=floor(intercept-rate*khi)-1
            high=ceil(intercept-rate*klo)+1
            for j in range(low,high+1):
                cut=floor((intercept-j)/rate)+1
                if klo<cut<=khi:
                    cuts.add(cut)
        for band in bands:
            if band['residue']==r:
                cuts.update(k for k in [(band['first']-r)//8,(band['last']-r)//8+1] if klo<k<=khi)
        points=sorted(cuts)
        for left,right in zip(points,points[1:]):
            if left>=right:
                continue
            vals=[]
            for k in (left,right-1):
                A=8*k+r
                exact=floor(lam11*A-b*gamma)==floor(lam12*A)
                derived=any(z['residue']==r and z['first']<=A<=z['last'] for z in bands)
                need(exact==derived, f'band mismatch {b=} {A=}')
                vals.append(exact)
            need(vals[0]==vals[1], 'nonconstant cell')
            cells+=1
    for band in bands:
        need(b*gamma+Delta*band['last']<1, 'strict coercivity')

out={'status':'PASS','gamma':str(gamma),'lambda12':str(lam12),'delta':str(Delta),
     'D':D,'T11':T,'exhaustive_floor_cells':cells,'families':[],
     'all_layer_content':'OPEN','E12':'OPEN','universal_Collatz':'OPEN'}
for a,b in families:
    minimum=a*T+3*b
    clipped=[]
    for z in intervals(b):
        lo=max(z['first'], minimum)
        first=lo+(z['residue']-lo)%8
        if first<=z['last']:
            clipped.append(dict(z,first=first,count=(z['last']-first)//8+1))
    out['families'].append({'a':a,'b':b,'total_clock_bands':clipped,
        'candidate_clocks':sum(z['count'] for z in clipped),
        'max_total_clock':max(z['last'] for z in clipped),
        'max_nonone_clock_sum':max(z['last'] for z in clipped)-a*T})
print(json.dumps(out,indent=2))
