#!/usr/bin/env python3
"""Exact hostile/positive controls for marked compression, palettes and clocks.

These are our own test mechanisms, NOT an implementation of the Vohra paper.
"""
import hashlib
import itertools as it
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[2]


def need(ok,msg):
    if not ok:
        raise AssertionError(msg)


def sat_controls():
    # The eight strict clauses on three variables forbid one assignment each.
    assignments=list(it.product([0,1],repeat=3))
    checks=0
    for mask in range(256):
        forbidden=[a for i,a in enumerate(assignments) if mask>>i&1]
        H=[sum(a==b for b in forbidden) for a in assignments]
        m=len(forbidden)
        values=[]
        for h in H:
            val=F(1)
            for j in range(1,m+1):val*=1-F(h,j)
            need(val==int(h==0),'zero-energy polynomial')
            values.append(val)
        need(sum(values)==8-m,'literal truth table versus projector count')
        # Moments recover the count, but the moment oracle is extra input.
        z=s.Symbol('z')
        poly=s.Poly(s.prod(1-z/s.Integer(j) for j in range(1,m+1)),z)
        recovered=sum(F(c)*sum(h**i[0] for h in H) for i,c in poly.terms())
        need(recovered==sum(values),'all required moments retain zero-energy mass')
        checks+=1
    # p>m prevents energy wrapping, but does not prevent count cancellation.
    need(7%7==0 and 7>1,'one strict clause has seven solutions, invisible count mod7')
    # Seven copies of every clause: genuinely UNSAT, but H=7 becomes zero in F7.
    need(all(7*sum(a==b for b in assignments)==7 for a in assignments),'energy wrap hostile')
    dims=[]
    for r in range(1,7):
        # r clauses on disjoint variable triples: all 2^r violation patterns occur.
        patterns=list(it.product([0,1],repeat=r))
        H=[sum(p) for p in patterns]
        cyclic=s.Matrix([[h**j for j in range(r+1)] for h in H]).rank()
        joint=s.Matrix([[int(p==q) for q in patterns] for p in patterns]).rank()
        need(cyclic==r+1 and joint==2**r,'single-generator versus joint-generator closure')
        dims.append(dict(clauses=r,hamiltonian_cyclic_dimension=cyclic,joint_dimension=joint))
    dense=[]
    for r in range(1,6):
        # Seven strict clauses per triple force that triple to 111.
        # Their zero-energy projector is a single assignment indicator.
        n=3*r;N=1<<n
        H=[sum(((a>>(3*j))&7)!=7 for j in range(r)) for a in range(N)]
        f=[int(h==0) for h in H]
        need(sum(f)==1,'linear-size strict 3CNF with a unique assignment')
        step=1
        while step<N:
            for lo in range(0,N,2*step):
                for j in range(lo,lo+step):
                    a,c=f[j],f[j+step];f[j],f[j+step]=a+c,a-c
            step*=2
        support=sum(c!=0 for c in f)
        need(support==N,'zero-energy projector has full Walsh support')
        dense.append(dict(variables=n,clauses=7*r,cyclic_dimension=r+1,walsh_support=support))
    return dict(all_three_variable_clause_subsets=checks,
                finite_field_hostiles=['count mod7','energy wrap mod7'],changing_generators=dims,
                dense_projector_family=dense)


def palette_controls():
    # Fixed singleton c=5, two pairs and a second singleton on the other colours.
    colours=set(range(5));options=[]
    for p in it.combinations(colours,2):
        for q in it.combinations(colours-set(p),2):
            singleton=next(iter(colours-set(p)-set(q)))
            options.append((frozenset(p),frozenset(q),singleton))
    need(len(options)==30,'thirty palettes after a safe singleton')
    blocker_masks=[]
    for cls in [0,1]:
        for pair in it.combinations(range(5),2):
            mask=sum(1<<i for i,o in enumerate(options) if o[cls]==set(pair))
            need(mask.bit_count()==3,'one blocker removes exactly three options')
            blocker_masks.append(mask)
    # Exhaust all 9-element distinct blocker banks (smaller banks extend to one).
    minimum=30
    banks=0
    for bank in it.combinations(blocker_masks,9):
        union=0
        for mask in bank:union|=mask
        minimum=min(minimum,30-union.bit_count());banks+=1
    need(minimum>=3,'nine blockers cannot kill every safe-singleton palette')
    # Independent partition count of C5 and C4 joined at one port.
    edges=[(i,(i+1)%5) for i in range(5)]+[(0,5),(5,6),(6,7),(7,0)]
    counts=Counter()
    def extend(word):
        if len(word)==8:
            counts[max(word)+1]+=1;return
        i=len(word)
        for c in range(min(6,max(word)+2)):
            if all(word[a]!=c for a,b in edges if b==i and a<i) and all(word[b]!=c for a,b in edges if a==i and b<i):
                extend(word+[c])
    extend([0])
    need(dict(counts)=={3:30,4:180,5:247,6:113},'unchanged marked-port relation')
    return dict(options=30,exhaustive_nine_blocker_banks=banks,
                minimum_surviving_options=minimum,boundary_patterns=dict(counts),
                at_most_five_patterns=sum(v for k,v in counts.items() if k<=5))


def clock_controls():
    # Ambient degree d contains every excursion out of degree<d needed here.
    count=0
    for d in range(1,21):
        E=s.zeros(d+1);t=s.zeros(d+1)
        for j in range(1,d+1):E[j-1,j]=j
        for j in range(d):t[j+1,j]=1
        P=s.diag(*([1]*d+[0]));Q=s.eye(d+1)-P
        e=(P*E*P)[:d,:d];z=(P*t*P)[:d,:d]
        leak=(P*E*Q*t*P-P*t*Q*E*P)[:d,:d]
        need(e*z-z*e==s.eye(d)-leak,'compression retains exact boundary commutator')
        need(s.trace(leak)==d,'all dimension appears in leakage trace')
        need(leak.rank()==1 and leak[d-1,d-1]==d,'rank-one leakage has size d')
        count+=1
    return dict(exact_clock_compressions=count)


if __name__=='__main__':
    result=dict(status='PROVED elementary mechanisms with FINITE-EXACT controls; not a paper implementation',
                sat=sat_controls(),palettes=palette_controls(),clock=clock_controls())
    out=ROOT/'05-knowledge/results/compression_context_controls_20261015g.json'
    out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print('sha256',hashlib.sha256(out.read_bytes()).hexdigest())
