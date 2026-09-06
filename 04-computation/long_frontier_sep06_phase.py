#!/usr/bin/env python3
"""Exact fixed-shape phase family; no census over the 384061 scales."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm, prod
from pathlib import Path
import hashlib
import json
import subprocess

Q = 91**6
PRIMES = (2,3,11,17,23,29)
P = prod(PRIMES)
V = (40341259,287243635,542783995,807423715,14321146945)
U = tuple(sorted(tuple(P//p for p in PRIMES)+(P,9*P)))
K = max(U)
D_MIN = min(gcd(a,b) for a,b in combinations(V,2))
LOW = F(Q,D_MIN)
HIGH = F(9*P)
WORD = {1:(2,1),2:(1,1),3:(2,3),4:(3,1)}
gates = 0


def check(test,label):
    global gates
    gates += 1
    if not test:
        raise RuntimeError(label)


def atlas(a,b):
    x = (a+b)//gcd(a,b)
    if x > 356:
        return False
    p = 2
    while p*p <= x:
        e = 0
        while x % p == 0:
            x //= p
            e += 1
        if e and (p % 3 != 2 or e > 2):
            return False
        p += 1
    return x == 1 or x % 3 == 2


def graph(row):
    remaining = set(range(13))
    groups = []
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        visited, pending = {start}, [start]
        while pending:
            i = pending.pop()
            for j in list(remaining):
                if atlas(row[i],row[j]):
                    remaining.remove(j)
                    visited.add(j)
                    pending.append(j)
        groups.append(sorted(visited))
    return groups


def ceildiv(a,b):
    return -((-a)//b)


def support(a,b,y):
    a,b = sorted((a,b))
    D = gcd(a,b)
    aa,bb = a//D,b//D
    delta = gcd(D,y)
    c,x = D//delta,y//delta
    check(bb <= Q,'internal pair height')
    r = (pow(aa,-1,bb)*x) % bb
    s = (x-aa*r)//bb
    low = max(ceildiv(-Q-r,bb),ceildiv(s-Q,aa))
    high = min((Q-r)//bb,(s+Q)//aa)
    return c,x,aa+bb,c<=Q and low<=high


def clearance(row,x):
    return min(min((n*x)%1,1-(n*x)%1) for n in row)


def coprime_count(low,high,primes):
    total = 0
    for mask in range(1<<len(primes)):
        d = prod(p for i,p in enumerate(primes) if mask>>i&1)
        total += (-1)**mask.bit_count()*(high//d-(low-1)//d)
    return total


def residue_count(low,high,r):
    # Independent count by 64 divisor conditions inside one class modulo5.
    total = 0
    for mask in range(1<<len(PRIMES)):
        d = prod(p for i,p in enumerate(PRIMES) if mask>>i&1)
        residue = r*pow(d,-1,5)%5
        count = (high//d-residue)//5-((low-1)//d-residue)//5
        total += (-1)**mask.bit_count()*count
    return total


def load_profiles():
    rel = '05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.json'
    root = Path(__file__).resolve().parent.parent
    file = root/rel
    raw = file.read_bytes() if file.exists() else subprocess.check_output(
        ['git','show','HEAD:'+rel],cwd=root)
    check(hashlib.sha256(raw).hexdigest() ==
          '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f',
          'inherited profile pin')
    return {int(d):{(c,tuple(gs)) for c,gs in z['profiles']}
            for d,z in json.loads(raw)['levels'].items()}


def main():
    q = (337,343,349,355)
    C = prod(q)
    check(V == tuple(sorted((C,)+tuple((356-x)*(C//x) for x in q))), 'five-star identity')
    check(P == 748374 and K == 6735366, 'larger shape constants')
    check(gcd(*V) == gcd(*U) == 1 and 1 not in V+U, 'primitive unitless shapes')
    check(all(v%2 for v in V), 'odd smaller shape')
    check(all(gcd(v,u)==1 for v in V for u in U), 'disjoint prime supports')
    check(all(gcd(u,5)==1 for u in U), 'larger fifth units')
    for x in q:
        check(atlas(C,(356-x)*(C//x)), 'five-star actual edge')
    for p,r in ((2,3),(3,17),(17,29),(2,23),(23,11)):
        check(atlas(P//p,P//r), 'cofactor actual edge')
    check(atlas(P//3,P) and atlas(P,9*P), 'two multiplier actual edges')
    max_sum = max((u+w)//gcd(u,w) for u,w in combinations(U,2))
    check(D_MIN == 115591 and min(V) == 349*D_MIN, 'exact opposing coefficient scale')
    check(max_sum == 262 < 349, 'all forward full supports below entry threshold')
    check(LOW*min(V)>355*K, 'no cross atlas edge throughout family')
    check(HIGH*sum(V)+sum(U)<Q*Q, 'entire family inside physical box')

    lower = Q//D_MIN+1
    upper = 9*P-1
    count = coprime_count(lower,upper,PRIMES+(5,))
    counts = {r:residue_count(lower,upper,r) for r in range(1,5)}
    check(count == sum(counts.values()) == 384061, 'exact scale count by two formulas')
    check((lower,upper)==(4912747,6735365), 'exact integer scale interval')

    # Whole continuous scale interval, with no scan of admissible scales.
    phase_bounds = {}
    for r,(k,delta) in WORD.items():
        check((2*r*k+delta)%10 == 5, 'half-grid identity in residue class')
        bounds = []
        for u in U:
            a = F((k*u)%5,5)
            lo = a+F(delta*u,10)/HIGH
            hi = a+F(delta*u,10)/LOW
            check(0<lo<=hi<1, 'affine image keeps its integer branch')
            bound = min(lo,1-hi)
            check(bound>F(1,6), 'whole-interval strict sixth clearance')
            bounds.append(bound)
        phase_bounds[r]=min(bounds)
    bound = min(phase_bounds.values())
    check(bound==F(2932887917,16555954870)>F(1,6), 'uniform certified family margin')

    # Complete scalar consumer comparison; delta=1 throughout by coprimality.
    native = []
    dual = []
    for u,w in combinations(U,2):
        D = gcd(u,w)
        aa,bb = u//D,w//D
        radius = Q*(aa+bb)-(aa-1)*(bb-1)
        for v in V:
            check(gcd(D,v)==1, 'forward cancellation absent')
            score = F(5*radius,63*K*v)
            check(score<1, 'all140 forward native gates fail')
            native.append(score)
            if w==K:
                check(30*D<=Q and 63*lcm(D,v)>5*Q, 'all35 endpoint gates fail in width')
    for a,b in combinations(V,2):
        D = gcd(a,b)
        for u in U:
            check(gcd(D,u)==1, 'dual cancellation absent')
            score=F(5*Q,63*K*D)
            check(score<1, 'all80 dual native gates fail in width')
            dual.append(score)
    check(5*HIGH<63*K and 8<42*max(V), 'both whole-arc gates fail for entire family')
    check(HIGH<=K, 'inherited seven-clock scale gate fails throughout')
    A23={1,3,5,6,10,11,12,13,17,18,20,22}
    check(not all(v%23 in A23 for v in V), 'incoming23 residue-class hypothesis absent')

    profiles=load_profiles()
    for p in PRIMES:
        check((p,(1,)*6) in profiles[6], 'only nontrivial symbolic profile retained')
    controls=[]
    for r in range(1,5):
        first=next(t for t in range(lower,lower+500) if t%5==r and gcd(t,5*P)==1)
        last=next(t for t in range(upper,upper-500,-1) if t%5==r and gcd(t,5*P)==1)
        for t in (first,last):
            row=tuple(t*v for v in V)+U
            check(gcd(*row)==1 and len(set(row))==13, 'actual physical primitive distinct row')
            check(sum(row)<Q*Q, 'actual physical box')
            check(graph(row)==[list(range(5)),list(range(5,13))], 'actual graph exactly5+8')
            check(all(gcd(t*v,u)==1 for v in V for u in U), 'actual cross coprimality')
            support_counts=[]
            for same,other,kind in ((tuple(t*v for v in V),U,'coefficient'),
                                   (U,tuple(t*v for v in V),'support')):
                n=0
                for a,b in combinations(same,2):
                    for y in other:
                        c,x,width,positive=support(a,b,y)
                        check(not positive, 'literal complete mixed-support test')
                        check(c>Q if kind=='coefficient' else x>Q*width,
                              'independent entry obstruction')
                        n+=1
                support_counts.append(n)
            check(support_counts==[80,140], 'both full support orientations')
            maxima=[]
            for tails in range(1,7):
                maximum=1
                for indices in combinations(range(13),13-tails):
                    occupied=set(indices)
                    c=gcd(*(row[i] for i in occupied))
                    gs=tuple(sorted(gcd(c,row[i]) for i in range(13) if i not in occupied))
                    check((c,gs) in profiles[tails], 'all inherited body profiles')
                    maximum=max(maximum,c)
                maxima.append(maximum)
            check(maxima==[1,1,1,1,1,29], 'profile maxima')
            k,delta=WORD[r]
            phase=F(k,5)+F(delta,10*t)
            margin=clearance(row,phase)
            check((t*phase)%1==F(1,2), 'literal smaller half phase')
            check(margin>bound>F(1,6), 'literal whole-row strict target')
            controls.append(dict(t=t,r=r,phase=str(phase),clearance=str(margin),
                                 physical_sum=sum(row),mixed_supports=support_counts))

    # A member hostile to the unrepaired nearest-fifth rule, not to LRC.
    hostile=4912753
    check(LOW<hostile<HIGH and gcd(hostile,5*P)==1 and hostile%5==3, 'hostile is in family')
    old_phase=F(4,5)+F(1,10*hostile)
    row=tuple(hostile*v for v in V)+U
    old_margin=clearance(row,old_phase)
    check(old_margin==F(309014,4912753)<F(1,14), 'nearest fifth lift fails actual target')
    repaired=F(2,5)+F(3,10*hostile)
    check(clearance(row,repaired)==F(1740589,9825506)>F(1,6), 'farther lift repairs the same row')

    manifest=dict(Q=Q,P=P,V=V,U=U,lower=lower,upper=upper,number_of_scales=count,
                  scale_counts_by_residue=counts,phase_word=WORD,
                  phase_bounds={r:str(v) for r,v in phase_bounds.items()},
                  uniform_clearance_bound=str(bound),
                  best_forward_native_score=str(max(native)),best_dual_native_score=str(max(dual)),
                  controls=controls,hostile=dict(t=hostile,phase=str(old_phase),
                  clearance=str(old_margin),repaired_phase=str(repaired)))
    print('PROVED fixed-shape actual-entry family, strict clearance>1/6; LRC14 remains OPEN')
    print(json.dumps(manifest,sort_keys=True,separators=(',',':')))
    print('EXPLICIT_GATES',gates)
    print('SEMANTIC_SHA256',hashlib.sha256(json.dumps(manifest,sort_keys=True,
                                separators=(',',':')).encode()).hexdigest())


if __name__=='__main__':
    main()
