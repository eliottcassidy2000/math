#!/usr/bin/env python3
"""Independent audit; no imports from producers, no removable assertions."""
from fractions import Fraction as F
from itertools import product
from hashlib import sha256
from pathlib import Path
import json

counts = {}
def check(condition, label):
    if not condition:
        raise RuntimeError(label)
    counts[label] = counts.get(label, 0) + 1

def step(n, b=3):
    x = b*n+1
    a = 0
    while x % 2 == 0:
        x //= 2
        a += 1
    return x, a

def trace(n, b, k):
    vals=[]; word=[]
    for _ in range(k):
        n,a=step(n,b); vals.append(n); word.append(a)
    return tuple(word),vals

def ceildiv(x,y):
    return -((-x)//y)

# Exhaust every residue, including excluded residues, rather than checking only
# sources selected by the producer's modular inverse.
for b in (3,5,9):
    for k in range(1,4):
        for word in product((1,2,3), repeat=k):
            A,B,C=1,1,0
            prefixes=[]
            for a in word:
                A,B,C=b*A,(2**a)*B,b*C+B
                prefixes.append((A,B,C))
            M=2*B
            residues=[n for n in range(1,M,2) if trace(n,b,k)[0]==word]
            check(len(residues)==1,'unique_exact_residue')
            r=residues[0]; t=(A*r+C)//B
            check(r==((B-C)*pow(A,-1,M))%M,'affine_residue_formula')
            check((A*r+C)%B==0 and t%2==1,'integer_odd_endpoint')
            for q in (0,1,3,10):
                got,ends=trace(r+M*q,b,k)
                check(got==word and ends[-1]==t+2*A*q,'periodic_cylinder_exactness')
            # Intersect ALL chronological source inequalities at once in final
            # word coordinates, independently of the recursive filter.
            for L,U in ((0,0),(0,11),(2,9),(8,8)):
                low,high=L,U
                for ap,bp,cp in prefixes:
                    d=(ap-bp)*M; h=(bp-ap)*r-cp
                    if d>0: low=max(low,ceildiv(h,d))
                    elif d<0: high=min(high,h//d)
                    elif h>0: high=low-1
                literal=[q for q in range(L,U+1)
                         if all(y>=r+M*q for y in trace(r+M*q,b,k)[1])]
                check(literal==list(range(low,high+1)),'all_prefix_interval_cut')
            # Enumerate next congruence rather than invert it to construct it.
            for a in range(1,6):
                P=2**a
                good=[q for q in range(P) if step(t+2*A*q,b)[1]==a]
                check(len(good)==1,'next_exact_guard_unique')
                s=good[0]
                check(s==((P//2-(b*t+1)//2)*pow(b*A,-1,P))%P,'next_guard_inverse_formula')
                for L,U in ((0,12),(1,1),(7,19)):
                    lp=max(0,ceildiv(L-s,P)); up=(U-s)//P
                    actual=[q for q in range(L,U+1) if step(t+2*A*q,b)[1]==a]
                    check(actual==[s+P*v for v in range(lp,up+1)],'interval_pullback')

# Literal finite source partitions with rational generating masses at nonempty
# prefixes. No geometric formula used to form the children.
for b in (3,5,9):
    for source_word in ((),(1,),(2,),(1,2)):
        sources=[n for n in range(3,258,2) if trace(n,b,len(source_word))[0]==source_word]
        groups={}
        for n in sources:
            endpoint=trace(n,b,len(source_word))[1][-1] if source_word else n
            a=step(endpoint,b)[1]
            groups.setdefault(a,[]).append(n)
        z=F(2,3)
        check(sum((sum((z**n for n in g),F(0)) for g in groups.values()),F(0))==sum((z**n for n in sources),F(0)), 'literal_source_mass_conservation')
        for g in groups.values():
            if len(g)>1:
                M=g[1]-g[0]
                check(g==list(range(g[0],g[-1]+1,M)),'child_progression')
                check(sum((z**n for n in g),F(0))==z**g[0]*(1-z**(M*len(g)))/(1-z**M),'geometric_mass')

# Check carry-sensitive fixed points using explicit rational branch composition.
def close(word,b=3):
    def apply(x):
        for a in word: x=(b*x+1)/2**a
        return x
    intercept=apply(F(0)); slope=apply(F(1))-intercept
    x=intercept/(1-slope)
    check(apply(x)==x,'rational_reclosure_identity')
    return x
check([close(w) for w in ((1,1,5),(1,3,3),(1,1,5,1),(1,3,3,1))]==[F(19,101),F(31,101),F(37,35),F(221,175)],'denominator_hostile')
x=13
for k in range(60):
    check(x>=13,'five_cycle_no_descent')
    x,a=step(x,5)
    check(a==(1,1,5)[k%3],'five_cycle_exact_word')
check(x==13,'five_cycle_return')

# 2x2 rational substochastic kernels (including genuinely randomized kernels),
# analytical inverse independent of producer Gaussian elimination.
def mul(A,v):
    return [sum((a*x for a,x in zip(row,v)),F(0)) for row in A]
def res(A,z,w):
    a,b=1-z*A[0][0],-z*A[0][1]
    c,d=-z*A[1][0],1-z*A[1][1]
    det=a*d-b*c
    return [(d*w[0]-b*w[1])/det,(-c*w[0]+a*w[1])/det]
columns=[(F(i,2),F(j,2)) for i in range(3) for j in range(3-i)]
kernels=[[[c0[0],c1[0]],[c0[1],c1[1]]] for c0,c1 in product(columns,repeat=2)]
w=[F(1,3),F(1,2)] # remaining 1/6 initially terminal
for Q in kernels:
    for z in (F(1,4),F(4,5)):
        g=res(Q,z,w)
        check(all(x>=0 for x in g),'rational_resolvent_positive')
        check(g==[x+z*y for x,y in zip(w,mul(Q,g))],'resolvent_equation')
        p=F(2,5)
        P=[[(1-p)*int(i==j)+p*Q[i][j] for j in range(2)] for i in range(2)]
        zz=z*p/(1-z*(1-p))
        check((1-z)*sum(res(P,z,w))==(1-zz)*sum(res(Q,zz,w)),'lazy_abel_identity')
    for R in kernels:
        alpha=F(1,3)
        L=[[(1-alpha)*Q[i][j]+alpha*R[i][j] for j in range(2)] for i in range(2)]
        q,l=w[:],w[:]
        for k in range(6):
            check(all(a >= (1-alpha)**k*b for a,b in zip(l,q)),'random_kernel_domination')
            q,l=mul(Q,q),mul(L,l)
# Initial-terminal convention and exact Abel identity on all 2-state maps.
for dest in product((0,1,2),repeat=2):
    Q=[[F(dest[j]==i+1) for j in range(2)] for i in range(2)]
    times=[]
    for start in (1,2):
        seen=set(); t=0; x=start
        while x and x not in seen:
            seen.add(x); x=dest[x-1]; t+=1
        times.append(None if x else t)
    z=F(3,4)
    d=(1-z)*sum(res(Q,z,w))
    check(d==sum(atom*(1-(z**t if t is not None else 0)) for atom,t in zip(w,times)),'abel_hitting_identity_with_terminal_mass')

# Construct inverse witnesses by search over all exponents 1..5, not producer's
# residue-case implementation. Test rank over GF(2), not rational elimination.
xs=[1]
for j in range(1,101):
    y=xs[-1]
    options=[((2**a*y-1)//3,a) for a in range(1,6) if (2**a*y-1)%3==0 and ((2**a*y-1)//3)%3 and ((2**a*y-1)//3)>y]
    check(bool(options),'inverse_witness_exists')
    x,a=options[0]; xs.append(x)
    check(step(x)==(y,a),'inverse_exact_valuation')
    check(3**j*x<=32**j,'inverse_size_bound')
    v=x
    for i in range(j):
        check(v>1,'inverse_no_premature_hit')
        v=step(v)[0]
    check(v==1,'inverse_exact_hitting_time')
N=100
rows=[0]*N
for j,x in enumerate(xs[1:]):
    for i in range(N):
        if x!=1: rows[i] |= 1<<j
        if x!=1: x=step(x)[0]
rank=0
for col in range(N):
    pivot=next((i for i in range(rank,N) if rows[i]>>col&1),None)
    if pivot is None: continue
    rows[rank],rows[pivot]=rows[pivot],rows[rank]
    for i in range(N):
        if i!=rank and rows[i]>>col&1: rows[i]^=rows[rank]
    rank+=1
check(rank==N,'literal_survival_rank_GF2')
for k in range(1,65):
    n=2**(k+1)-1
    for i in range(k):
        check(n==3**i*2**(k+1-i)-1 and n>1,'norm_one_witness')
        n,a=step(n)
        check(a==1,'norm_one_exact_guard')
    check(n>1,'norm_one_final_survival')

print(json.dumps({'status':'PASS','checks':counts,'total_checks':sum(counts.values()),'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'inverse_rank_GF2':rank,'scope':'finite independent controls only; proof audit in REPORT.md'},indent=2,sort_keys=True))
