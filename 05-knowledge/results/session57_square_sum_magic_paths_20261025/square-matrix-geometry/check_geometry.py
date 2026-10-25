#!/usr/bin/env python3
"""Exact symbolic and finite-field controls; includes a cited known K4,4 certificate."""
import itertools, json
from math import isqrt
from collections import Counter
from fractions import Fraction
from pathlib import Path
import sympy as sp

def rankmod(A,p):
    A=[[int(x)%p for x in row] for row in A]
    if not A: return 0
    h=0
    for c in range(len(A[0])):
        z=next((i for i in range(h,len(A)) if A[i][c]),None)
        if z is None: continue
        A[h],A[z]=A[z],A[h]
        inv=pow(A[h][c],-1,p)
        A[h]=[(x*inv)%p for x in A[h]]
        for i in range(len(A)):
            if i!=h and A[i][c]:
                fac=A[i][c]; A[i]=[(x-fac*y)%p for x,y in zip(A[i],A[h])]
        h+=1
        if h==len(A): break
    return h

def data(r,s):
    L=[]
    for i in range(r):
        for j in range(s):
            v=[0]*(r+s-1);v[i]=1
            if j:v[r+j-1]=1
            L.append(v)
    C=[]
    for i in range(1,r):
        for j in range(1,s):
            v=[0]*(r*s)
            for k,c in [(i*s+j,1),(0,1),(i*s,-1),(j,-1)]:v[k]+=c
            C.append(v)
    return L,C

out={'symbolic':[], 'finite_controls':[]}
for r,s in [(2,2),(2,3),(3,3),(3,4),(4,4)]:
    L,C=data(r,s);N=r*s;m=(r-1)*(s-1)
    assert sp.Matrix(L).rank()==r+s-1
    assert sp.Matrix(C).rank()==m
    assert sp.Matrix(C)*sp.Matrix(L)==sp.zeros(m,r+s-1)
    assert len({tuple(v) for v in L})==N
    for I in itertools.combinations(range(r),2):
        for J in itertools.combinations(range(s),2):
            assert sp.Matrix([L[i*s+j] for i in I for j in J]).rank()==3
    out['symbolic'].append(dict(r=r,s=s,codimension=m,dimension=r+s-2,
        canonical_twist=2*m-N,degree=2**m,projective_cover_degree=2**(N-1),
        singular_dimension=None if (r,s)==(2,2) else r+s-5))
    for p in [3,5]:
        # Exhaust all additive bases over F3; F5 exhaustive up to 3x3,
        # deterministic samples for larger dimensions.
        n=r+s-1
        bases=itertools.product(range(p),repeat=n)
        exhaustive=(p==3 or n<=5)
        if not exhaustive:bases=itertools.islice(bases,0,p**n,37)
        tested=singular=nonzero_bases=0
        for v in bases:
            z=[sum(x*y for x,y in zip(row,v))%p for row in L]
            if not any(z):continue
            nonzero_bases+=1
            zero={k for k,x in enumerate(z) if not x}
            rank_drop=rankmod([[x if k not in zero else 0 for k,x in enumerate(row)] for row in C],p)<m
            four_zero=any(all(i*s+j in zero for i in I for j in J)
                for I in itertools.combinations(range(r),2) for J in itertools.combinations(range(s),2))
            assert rank_drop==four_zero,(r,s,p,v)
            tested+=1;singular+=rank_drop
        out['finite_controls'].append(dict(r=r,s=s,p=p,exhaustive=exhaustive,
            nonzero_bases_tested=tested,singular_supports_with_multiplicity=singular))

C,x,y=sp.symbols('C x y'); S=[0,x,y,x+y]
forms=[sp.expand(C+t-s) for s in S for t in S]
mult=Counter(forms)
assert len(mult)==9 and sorted(mult.values())==[1,1,1,1,2,2,2,2,4]
coeff=sp.Matrix([[f.coeff(C),f.coeff(x),f.coeff(y)] for f in forms])
assert coeff.rank()==3
L44,C44=data(4,4)
assert sp.Matrix(C44)*coeff==sp.zeros(9,3)
magic=sp.Matrix([[C+x,C-x-y,C+y],[C-x+y,C,C+x-y],[C-y,C+x+y,C-x]])
sums=[sum(magic.row(i)) for i in range(3)]+[sum(magic.col(j)) for j in range(3)]
sums += [sum(magic[i,i] for i in range(3)),sum(magic[i,2-i] for i in range(3))]
assert all(sp.expand(a-3*C)==0 for a in sums)
assert set(magic)==set(forms)
# Coefficient image of the magic parameterization has rank 3, and magic
# row/column/diagonal equal-sum equations have kernel of dimension 3.
w=sp.symbols('w0:9');W=sp.Matrix(3,3,w)
lines=[sum(W.row(i)) for i in range(3)]+[sum(W.col(j)) for j in range(3)]
lines += [sum(W[i,i] for i in range(3)),sum(W[i,2-i] for i in range(3))]
constraints=sp.Matrix([[sp.expand(line-lines[0]).coeff(t) for t in w] for line in lines[1:]])
assert constraints.rank()==6
out['decoration']={'base_rank':3,'base_codimension':4,'surface_dimension':2,
 'distinct_forms':9,'sign_components':2**(16-9),'magic_constraints_rank':6,
 'multiplicities':{str(k):v for k,v in sorted(mult.items(),key=lambda kv:str(kv[0]))}}
a=[Fraction(1),Fraction(49)];b=[Fraction(0),Fraction(15)];c=Fraction(1,3);D=3
A=[D*D*(t-c) for t in a];B=[D*D*(t+c) for t in b]
Q=[[3,12],[21,24]]
assert len(set(A+B))==4 and all(t>0 and t.denominator==1 for t in A+B)
assert all(A[i]+B[j]==Q[i][j]**2 for i in range(2) for j in range(2))
out['gauge_control']={'row_vertices':list(map(int,A)),'column_vertices':list(map(int,B)),'roots':Q}
# Independently verify the numerical example Choudhry records from Guy.
KA=[324,54756,119716,264196];KB=[0,79200,227205,1258560]
KQ=[[isqrt(a+b) for b in KB] for a in KA]
assert all(KQ[i][j]**2==KA[i]+KB[j] for i in range(4) for j in range(4))
KAp=[a-1 for a in KA];KBp=[b+1 for b in KB]
assert len(set(KAp+KBp))==8 and min(KAp+KBp)>0
assert all(KQ[i][j]**2==KAp[i]+KBp[j] for i in range(4) for j in range(4))
parallelogram_counts=[sum(v[0]+v[3]==v[1]+v[2] for v in itertools.permutations(Z)) for Z in [KA,KB]]
reflection_count=sum(len({KA[i]+v[i] for i in range(4)})==1 for v in itertools.permutations(KB))
assert parallelogram_counts==[0,0] and reflection_count==0
out['known_K44_certificate']={'source':'https://arxiv.org/html/2508.07806v1',
 'attribution':'Numerical example recorded by Choudhry from Guy; Jean Lagrange attribution applies to the accompanying parametric family',
 'original_A':KA,'original_B':KB,'positive_A':KAp,'positive_B':KBp,'roots':KQ,
 'verified_cross_squares':16,'parallelogram_permutations_tested_per_part':24,
 'parallelogram_ordering_counts':parallelogram_counts,'reflected_matchings_tested':24,
 'constant_reflected_matchings':reflection_count,'conclusion':'U44(Q) nonempty; this point is not MSOS-decorated'}
out['status']='PASS: all assertions; geometric/arithmetic scope is stated in geometry.md'
Path(__file__).with_name('check_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
