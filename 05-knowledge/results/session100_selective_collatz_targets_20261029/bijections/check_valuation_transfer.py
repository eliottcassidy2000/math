#!/usr/bin/env python3
"""Finite exact controls; no asserts, all checks remain active under -O."""
from fractions import Fraction as F
from itertools import product
from math import gcd
import json

checks=0
def require(ok, label):
    global checks
    checks+=1
    if not ok: raise RuntimeError(label)

def v2(n):
    if not n: raise ValueError('infinite valuation')
    n=abs(n)
    return (n & -n).bit_length()-1

def carry(b,w):
    D,C=0,0
    for a in w:
        C=b*C+2**D
        D+=a
    return D,C

def residue(b,w):
    D,C=carry(b,w); M=2**(D+1)
    return ((2**D-C)*pow(b**len(w),-1,M))%M

def rational_residue(x,M):
    return (x.numerator*pow(x.denominator,-1,M))%M

def inverse(b,w,y):
    for a in reversed(w): y=(2**a*y-1)/b
    return y

def itinerary(b,x,J):
    out=[]
    for _ in range(J):
        z=b*x+1
        a=v2(z.numerator)-v2(z.denominator)
        require(a>=1,'odd domain')
        out.append(a); x=z/2**a
    return tuple(out),x

def periodic(b,w):
    D,C=carry(b,w)
    return F(C,2**D-b**len(w))

words=[w for J in range(1,5) for w in product(range(1,6),repeat=J)]
for b in (3,5,7):
    for w in words:
        D,C=carry(b,w); M=2**(D+1); r=residue(b,w)
        require(r%2==1 and 0<r<M,'odd least representative')
        z=inverse(b,w,F(1))
        require(rational_residue(z,M)==r,'independent backwards residue')
        for k in (0,1,3):
            n=F(r+k*M)
            actual,end=itinerary(b,n,len(w))
            require(actual==w,'exact forward cylinder word')
            require(end==(b**len(w)*n+C)/2**D,'affine identity')
        p=periodic(b,w)
        actual,end=itinerary(b,p,len(w))
        require(actual==w and end==p,'periodic exact guard')
        rr=[residue(b,w[:j]) for j in range(1,len(w)+1)]
        require(rr==sorted(rr),'monotone prefix residues')
        if p.denominator>1:
            q=p.denominator
            for a in (1,2,5,11):
                z=inverse(b,(a,),p)
                require(z.denominator%q==0,'denominator persistence')
        for c in (3,5,7):
            q=periodic(c,w)
            require(itinerary(c,q,len(w))[0]==w,'parameter transfer')

# Finite source prefixes and reconstruction once modulus exceeds n.
for b in (3,5,7):
    for n in range(1,501,2):
        w,_=itinerary(b,F(n),20)
        for j in range(1,21):
            D=sum(w[:j]); r=residue(b,w[:j])
            require(r==n%(2**(D+1)),'integer-source residue')
            if 2**(D+1)>n: require(r==n,'integer stabilization')

hostile=(1,1,5)
require(periodic(5,hostile)==13,'5 cycle source')
require(periodic(3,hostile)==F(19,101),'positive rational transport')
rotations=[]
for i in range(3):
    w=hostile[i:]+hostile[:i]
    p=periodic(3,w)
    require(p.denominator==101 and p>0,'rotated rational denominator')
    rotations.append(str(p))
    for J in range(5):
        for pre in product(range(1,5),repeat=J):
            x=inverse(3,pre,p)
            require(x.denominator%101==0,'arbitrary prehistory denominator')
            require(itinerary(3,x,J)[1]==p,'prehistory guard')

residues=[]
for k in range(1,13):
    w=hostile*k; M=2**(sum(w)+1); r=residue(3,w)
    require(r==rational_residue(F(19,101),M),'repeated-word rational lift')
    require(101*r>=M-19,'quantitative denominator certificate')
    residues.append(r)
require(periodic(3,(2,))==1 and periodic(5,(2,))==-1,'positive-to-negative control')
require(periodic(3,(1,))==-1,'negative integer excluded by boundedness')
for k in range(1,21):
    require(residue(3,(1,)*k)==2**(k+1)-1,'negative residues diverge')

# Isometry on periodic sample words, measured exactly as rational differences.
samples=[w for w in words if len(w)<=2]
for u in samples:
    for w in samples:
        x=periodic(3,u)-periodic(3,w)
        y=periodic(5,u)-periodic(5,w)
        require((x==0)==(y==0),'same infinite coding equality')
        if x:
            require(v2(x.numerator)-v2(x.denominator)==v2(y.numerator)-v2(y.denominator),'parameter isometry')

print(json.dumps({'status':'PASS','checks':checks,'word_universe':{'lengths':[1,2,3,4],'digits':[1,2,3,4,5],'count':len(words),'parameters':[3,5,7]},'integer_sources':'odd 1..499, 20 valuations each','preperiod_universe':'length 0..4, digits 1..4, all three hostile rotations','transported_3_cycle':rotations,'hostile_3_residues_first_12_repetitions':residues,'no_claim':'No all-source convergence or aperiodic exclusion.'},sort_keys=True,indent=2))
