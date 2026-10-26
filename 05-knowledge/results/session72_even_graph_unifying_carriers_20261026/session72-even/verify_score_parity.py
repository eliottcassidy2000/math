#!/usr/bin/env python3
"""Independent score-parity specialization: exact finite checks only."""
from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path
import ast,json
CHECKS=Counter()
def check(v,tag):
    CHECKS[tag]+=1
    if not v: raise RuntimeError(tag)
def choose(n,k): return comb(n,k) if 0<=k<=n else 0
def kraw(n,a,r): return sum((-1)**j*choose(r,j)*choose(n-r,a-j) for j in range(a+1))
def rank_mod(a,p):
    a=[[v%p for v in row] for row in a]; r=0
    for j in range(len(a[0])):
        i=next((i for i in range(r,len(a)) if a[i][j]),None)
        if i is None: continue
        a[r],a[i]=a[i],a[r]; z=pow(a[r][j],p-2,p)
        a[r]=[v*z%p for v in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][j]:
                z=a[i][j]; a[i]=[(v-z*w)%p for v,w in zip(a[i],a[r])]
        r+=1
    return r

def main():
    rows=[]
    for n in range(65):
        m=choose(n,2); weights=list(range(m%2,n+1,2))
        atoms=[a for a in range(n//2+1) if not (2*a==n and a%2)]
        check(len(weights)==len(atoms),'cut_admissibility_count')
        A=[]
        for r in weights:
            A.append([choose(r,2) if s==r-2 else choose(n-r,2) if s==r+2 else r*(n-r) if s==r else 0 for s in weights])
        K=[[kraw(n,a,r) for a in atoms] for r in weights]
        eigen=[m-2*a*(n-a) for a in atoms]
        for i,r in enumerate(weights):
            check(sum(A[i])==m,'row_sum')
            for j,a in enumerate(atoms):
                check(sum(A[i][k]*K[k][j] for k in range(len(weights)))==eigen[j]*K[i][j],'krawtchouk_eigenvector')
            for j,s in enumerate(weights):
                check(choose(n,r)*A[i][j]==choose(n,s)*A[j][i],'stationary_balance')
        check(rank_mod(K,1000000007)==len(weights),'exact_modular_rank_certificate')
        check(len(set(eigen))==len(eigen),'simple_spectrum')
        if n%4==2:
            check(all(kraw(n,n//2,r)==0 for r in weights),'excluded_balanced_mode_zero')
        gaps=[m-z for z in eigen]
        if n>=3: check(min(x for x in gaps if x)==2*(n-1),'gap_twice_n_minus_one')
        else: check(gaps==[0],'small_boundary')
        if n<=6:
            edges=list(combinations(range(n),2)); fibres=Counter(); cycles=0
            p0=sum(((n-1-v)%2)<<v for v in range(n))
            for mask in range(1<<m):
                score=[0]*n; boundary=0
                for i,(u,v) in enumerate(edges):
                    score[v if mask>>i&1 else u]+=1
                    if mask>>i&1: boundary ^= (1<<u)|(1<<v)
                p=sum((s%2)<<v for v,s in enumerate(score))
                check(p==p0^boundary,'literal_tournament_boundary_map')
                check(p.bit_count()%2==m%2,'literal_score_parity')
                fibres[p]+=1
                if boundary==0: cycles+=1
            possible=[p for p in range(1<<n) if p.bit_count()%2==m%2]
            check(set(fibres)==set(possible),'onto_parity_affine_space')
            check(set(fibres.values())=={cycles},'cycle_space_fibre_size')
            for p in possible:
                r=p.bit_count(); out=Counter((p^(1<<u)^(1<<v)).bit_count() for u,v in edges)
                expected={s:A[weights.index(r)][j] for j,s in enumerate(weights) if A[weights.index(r)][j]}
                check(dict(out)==expected,'literal_pair_flip_rates')
        rows.append({'n':n,'types':len(weights),'score_weights':weights,'cut_sizes':atoms,'laplacian_spectrum':gaps})
    tree=ast.parse(Path(__file__).read_text())
    check(not any(isinstance(a,ast.Assert) for a in ast.walk(tree)),'no_asserts')
    check(not any(isinstance(a,ast.Constant) and isinstance(a.value,float) for a in ast.walk(tree)),'no_floats')
    print(json.dumps({'status':'PASS','checks':dict(sorted(CHECKS.items())),'total_checks':sum(CHECKS.values()),'literal_tournament_range':'0..6','tridiagonal_matrix_range':'0..64','rows':rows},indent=2,sort_keys=True))
if __name__=='__main__': main()
