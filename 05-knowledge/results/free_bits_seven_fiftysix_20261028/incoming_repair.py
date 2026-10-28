#!/usr/bin/env python3
"""Exact hostile controls and repairs for incoming THM-4840/HYP-9438."""
from collections import Counter
from itertools import combinations, product, permutations
from math import isqrt, prod
import json
checks = 0
def ck(ok, msg):
    global checks
    checks += 1
    if not ok: raise RuntimeError(msg)
def rank(xs):
    basis = {}
    for x in xs:
        while x:
            bit = x.bit_length()
            if bit not in basis: basis[bit] = x; break
            x ^= basis[bit]
    return len(basis)
def det(A):
    n=len(A); ans=0
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        ans += sign*prod(A[i][p[i]] for i in range(n))
    return ans
# Both sixteenth-coordinate constructions, over all K6 edge subsets.
edges=list(combinations(range(6),2)); weights_degree=Counter(); weights_extended=Counter()
for mask in range(1<<15):
    d=0
    for i,(u,v) in enumerate(edges):
        if mask>>i&1: d ^= (1<<u)|(1<<v)
    if d not in (0,63): continue
    w=mask.bit_count(); degree_bit=int(d==63); parity_bit=w%2
    weights_degree[w+degree_bit]+=1; weights_extended[w+parity_bit]+=1
    ck((w+parity_bit)%2==0,'correct extension must be even')
ck(sum(weights_degree.values())==2048,'degree extension dimension')
ck(min(w for w in weights_degree if w)>0,'nonzero code')
ck(min(w for w in weights_degree if w)==3,'triangle survives wrong extension')
ck(min(w for w in weights_extended if w)==4,'repaired extended distance')
# A regular binary column matroid whose particular matrix is not TU-signable.
A=[[1,1,1,0],[1,1,0,1],[1,0,1,1]]
cols=[sum(A[i][j]<<i for i in range(3)) for j in range(4)]
ck(cols==[7,3,5,6] and rank(cols)==3 and cols[1]^cols[2]^cols[3]==0,'triangle plus coloop')
ck(set(cols+[1,2,4])==set(range(1,8)),'identity extension is Fano')
slots=[(i,j) for i in range(3) for j in range(4) if A[i][j]]
possible=0
for signs in product((-1,1),repeat=len(slots)):
    B=[[0]*4 for _ in range(3)]
    for (i,j),s in zip(slots,signs): B[i][j]=s
    tu=all(abs(det([[B[i][j] for j in cs] for i in rs]))<=1
           for k in (2,3) for rs in combinations(range(3),k) for cs in combinations(range(4),k))
    possible += tu
    ck(not tu,'unexpected TU signing')
# Explicit cographic K3,3 columns from four fundamental cycles.
K=[(i,j) for i in range(3) for j in range(3,6)]
tree={0,1,2,3,6}; chords=[i for i in range(9) if i not in tree]
cycle_rows=[]
for chord in chords:
    start,end=K[chord]; stack=[(start,0,{start})]; found=None
    while stack:
        u,path,seen=stack.pop()
        if u==end: found=path|(1<<chord); break
        for e in tree:
            a,b=K[e]
            if u not in (a,b): continue
            v=b if u==a else a
            if v not in seen: stack.append((v,path|(1<<e),seen|{v}))
    if found is None: raise RuntimeError('tree path missing')
    cycle_rows.append(found)
co_cols=[sum(((row>>i)&1)<<j for j,row in enumerate(cycle_rows)) for i in range(9)]
prime_table=[None,7,2,5,19,89,59,29,251,13,83,41,139,5669,419,839]
bank=[prime_table[c] for c in co_cols]
ck(len(set(co_cols))==9 and rank(co_cols)==4,'simple rank-four cographic bank')
for p in bank: ck(all(p%d for d in range(2,isqrt(p)+1)),'prime')
for mask in range(1<<9):
    xs=[co_cols[i] for i in range(9) if mask>>i&1]
    # Graphic rank of complementary edges = 6-components.
    parent=list(range(6))
    def find(x):
        while parent[x]!=x: x=parent[x]
        return x
    for i,(u,v) in enumerate(K):
        if not(mask>>i&1): parent[find(u)]=find(v)
    comps=len({find(i) for i in range(6)})
    ck(rank(xs)==mask.bit_count()+1-comps,'cographic rank identity')
    syndrome=0; sigma=1
    for i in range(9):
        if mask>>i&1: syndrome^=co_cols[i]; sigma*=bank[i]+1
    ck((syndrome==0)==(isqrt(sigma)**2==sigma),'actual divisor-sum square test')
# Ideal evaluations and the failed arbitrary-half claim at norm three.
ck((2+1)%3==0 and (1+1+1)%3==0,'norm-three root')
ck((3+4)%7==0 and (4*4+4+1)%7==0,'norm-seven root is four')
ck((2+3)%5==0 and (3*3+1)%5==0,'norm-five root is three')
ck({1,(-1)%3,1}=={1,2},'mixed unit half has opposite residue pair')
# Product-weight Fourier assembly at all syndromes, including signed weights.
for r in (2,3,4):
    n=(1<<r)-1; weights=[(i+1,2*i-3) for i in range(n)]
    direct=[0]*(1<<r)
    for mask in range(1<<n):
        syndrome=0; weight=1
        for i,(f0,f1) in enumerate(weights):
            b=(mask>>i)&1; weight*=f1 if b else f0
            if b: syndrome^=i+1
        direct[syndrome]+=weight
    products=[prod(f0+(-1)**((b&(i+1)).bit_count())*f1 for i,(f0,f1) in enumerate(weights)) for b in range(1<<r)]
    for syndrome in range(1<<r):
        numerator=sum((-1)**((b&syndrome).bit_count())*value for b,value in enumerate(products))
        ck(numerator==(1<<r)*direct[syndrome],'exact weighted Fourier query')
result={'checks':checks,'degree_bit_extension_weights':dict(sorted(weights_degree.items())),
        'correct_even_extension_weights':dict(sorted(weights_extended.items())),
        'regular_matrix_signings_tested':2**len(slots),'tu_signings':possible,
        'K33_edges':K,'K33_cographic_columns':co_cols,'K33_prime_bank':bank,
        'K33_max_prime':max(bank),'K33_all_subset_checks':512}
print(json.dumps(result,indent=2,sort_keys=True))
