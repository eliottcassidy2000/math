#!/usr/bin/env python3
"""Exact local-probability verification and finite rational construction bounds."""
from itertools import permutations, product, combinations
from fractions import Fraction
from math import factorial
from pathlib import Path
import json
from verify_maxima import adjacency, require
HERE=Path(__file__).resolve().parent
A=adjacency(4,2,11)
require([sum(r) for r in A]==[3]*8,'regular local frame')
# A random image is implemented by pulling back each queried signed vertex.
frames=[]
for p in permutations(range(4)):
    for bits in range(16):
        f=[2*p[i]+(s^(bits>>i&1)) for i in range(4) for s in range(2)]
        frames.append([[A[f[i]][f[j]] for j in range(8)] for i in range(8)])
require(len(frames)==384,'signed permutation count')
counts={}
for r,target in [(1,Fraction(1,2)),(2,Fraction(5,16)),(3,Fraction(3,16))]:
    checked=0
    for inds in permutations(range(4),r+1):
        for signs in product(range(2),repeat=r+1):
            vs=[2*i+s for i,s in zip(inds,signs)]
            actual=Fraction(sum(all(F[u][v] for u,v in zip(vs,vs[1:])) for F in frames),384)
            require(actual==target,f'path probability {r}');checked+=1
    counts[f'{r}_edge_path']={'probability':str(target),'configurations_checked':checked}
checked=0
for inds in permutations(range(4)):
    for signs in product(range(2),repeat=4):
        vs=[2*i+s for i,s in zip(inds,signs)]
        actual=Fraction(sum(F[vs[0]][vs[1]]*F[vs[2]][vs[3]] for F in frames),384)
        require(actual==Fraction(1,4),'disjoint arc pair');checked+=1
counts['disjoint_edge_pair']={'probability':'1/4','configurations_checked':checked}
# An explicit (16,4,1) affine-plane design, independently checked.
def mul(a,b):
    ans=0
    while b:
        if b&1:ans^=a
        b>>=1;a<<=1
        if a&4:a^=7
    return ans
points=list(product(range(4),repeat=2));directions=[(1,t) for t in range(4)]+[(0,1)]
blocks=set()
for d in directions:
    for a in points:
        blocks.add(tuple(sorted((a[0]^mul(t,d[0]))*4+(a[1]^mul(t,d[1])) for t in range(4))))
require(len(blocks)==20,'affine block count')
pairs=[tuple(sorted(e)) for B in blocks for e in combinations(B,2)]
require(len(pairs)==120 and len(set(pairs))==120,'affine pair partition')
triples={R for B in blocks for R in combinations(B,3)}
d2=sum(not(set(U)&set(V)) for U,V in combinations(triples,2))
require(d2==80*(80-3*16+8)//2,'exact disjoint triple pair formula')
finite=[]
for k in [16,28,100,1000,10000,100000]:
    v=k-(k-4)%12;q=v*(v-1)//12;total=Fraction(0);terms=[]
    for r in range(min(12,k//3)+1):
        factors=[4*q-3*j*(k-1) for j in range(r)]
        if any(x<0 for x in factors):term=Fraction(0)
        else:
            prod=1
            for x in factors:prod*=x
            den=1
            for i in range(2*r):den*=k-i
            term=Fraction(3**r*prod,2**r*factorial(r)*den)
        terms.append(str(term));total+=term
    finite.append({'k':k,'v':v,'K4_blocks':q,'truncation_R':len(terms)-1,'ratio_lower_bound':str(total),'decimal':float(total)})
# Verify the rational termwise estimates used in the explicit O(1/k) bounds.
quantitative_checks=0
for k in range(4,257):
    v=k-(k-4)%12;q=v*(v-1)//12
    for r in range(k//3+6):
        greedy=Fraction(0)
        if r<=k//3:
            factors=[4*q-3*j*(k-1) for j in range(r)]
            if all(x>=0 for x in factors):
                prod=1;den=1
                for x in factors:prod*=x
                for j in range(2*r):den*=k-j
                greedy=Fraction(3**r*prod,2**r*factorial(r)*den)
        ar=Fraction(1,2**r*factorial(r))
        rhs=ar*(1-Fraction(23*r, k)-Fraction(9*r*(r-1),2*k))
        require(greedy>=rhs,'uniform quantitative termwise bound')
        if k%12==4:
            rhs2=ar*(1-Fraction(r,k)-Fraction(9*r*(r-1),2*k))
            require(greedy>=rhs2,'admissible-order quantitative termwise bound')
        quantitative_checks+=1
output={'local_frame':{'A':2,'b':11,'outdegrees':[3]*8,'signed_permutations':384},
        'local_probabilities':counts,'affine16_design':{'blocks':[list(B) for B in sorted(blocks)],'all_120_pairs_once':True,'triples':80,'disjoint_triple_pairs':d2},'finite_lower_bounds':finite,'quantitative_termwise_checks':quantitative_checks}
(HERE/'local_design_verification.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))
