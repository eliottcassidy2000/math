"""Small exact controls for the central reversible-core construction."""
from collections import Counter
from fractions import Fraction
from math import gcd
from pathlib import Path
import json

def factors(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out

def vp(a,p,e):
    v=0
    while v<e and a%p==0:
        v+=1
        a//=p
    return v

examples=[]
counts=Counter()
for n in range(3,102,2):
    fs=factors(n)
    p=min(fs)
    k=(n-1)//2
    block=[[a*b%n for b in (k,k+1)] for a in (k,k+1)]
    d=pow(4,-1,n)
    assert block==[[d,n-d],[n-d,d]]
    assert list(map(sum,block))==[n,n]
    first=next((r for r in range(k) if gcd(k-r,n)>1),k)
    assert 2*first==p-1
    for r in range(k):
        assert gcd(k-r,n)==gcd(k+1+r,n)==gcd(2*r+1,n)
    units=[a for a in range(n) if gcd(a,n)==1]
    for a in range(n):
        g=gcd(a,n); m=n//g
        fibers=Counter(a*x%n for x in range(n))
        assert len(fibers)==m and set(fibers.values())=={g}
        for x in range(n):
            y=a*x%n
            u=(pow(a//g,-1,m)*(y//g))%m if m>1 else 0
            t=x//m
            assert x==u+t*m
            counts['decoder_cases']+=1
        orbit={a*u%n for u in units}
        assert orbit=={x for x in range(n) if gcd(x,n)==g}
        for b in range(n):
            for q,e in fs.items():
                assert vp(a*b,q,e)==min(e,vp(a,q,e)+vp(b,q,e))
            counts['monoid_pairs']+=1
    for q,e in fs.items():
        observed=sum(Fraction(vp(a,q,e),n-1) for a in range(1,n))
        formula=sum(Fraction(n//q**j-1,n-1) for j in range(1,e+1))
        assert observed==formula
    counts['odd_moduli']+=1
    if n in (5,9,11,15,25,45):
        examples.append({'n':n,'central_block':block,'first_bad_radius_or_boundary':first,
                         'unit_core_width':2*first,'nonzero_defect_types':len({gcd(a,n) for a in range(1,n)})})
# Independent partition refinement of the zero-recognizing two-sheet automaton.
# Begin with the zero/nonzero observable; repeatedly retain successor classes.
for n in range(5,102,2):
    if n%3==0:
        continue
    classes=[int(x!=0) for x in range(n)]
    while True:
        signatures=[(classes[x],classes[(3*x+1)%n],classes[(3*x-1)%n]) for x in range(n)]
        names={sig:i for i,sig in enumerate(dict.fromkeys(signatures))}
        refined=[names[sig] for sig in signatures]
        if len(set(refined))==len(set(classes)):
            classes=refined
            break
        classes=refined
    assert len(set(classes))==n
    counts['affine_minimal_automata']+=1
assert {(3*x%9,x%3) for x in (0,3,6)}=={(0,0)}
assert gcd(1,5)==gcd(3,5) and gcd(3*1+1,5)!=gcd(3*3+1,5)
report={'status':'FINITE-EXACT','universe':'All odd n from 3 through 101; every a,b,x in Z/n for indicated tests',
        'counts':dict(counts),'examples':examples,
        'positive_controls':['prime n=5,11 have complete reversible central interval'],
        'hostile_controls':['n=9,a=3: low-digit sidecar cannot recover input',
                            'n=5: valuation quotient does not process 3x+1'],
        'all_assertions_passed':True}
Path(__file__).with_name('checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
