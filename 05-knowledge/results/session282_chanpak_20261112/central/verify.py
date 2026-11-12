#!/usr/bin/env python3
"""Exact central-convolution controls; no check is disabled by python -O."""
from pathlib import Path
from itertools import product,combinations
from math import prod,gcd,lcm
from functools import reduce
import importlib.util,json
ROOT=Path(__file__).resolve().parents[4]
p=ROOT/'05-knowledge/results/session275_exact_selection_20261111/switches/carrier.py'
s=importlib.util.spec_from_file_location('carrier',p); c=importlib.util.module_from_spec(s);s.loader.exec_module(c)
def ck(ok,label):
    if not ok:raise RuntimeError(label)
def literal_hp(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not(mask>>v&1):continue
            previous=mask^(1<<v)
            if previous:
                dp[mask][v]=sum(dp[previous][u] for u in range(n) if previous>>u&1 and T[u]>>v&1)
    return sum(dp[-1])
def carry_select(h,E,z):
    D=prod(h);w=[D//a for a in h]
    residues=[z*pow(b,-1,a)%a for a,b in zip(h,w)]
    rem=z-sum(b*r for b,r in zip(w,residues))
    ck(rem%D==0,'integral absolute carry')
    wanted=rem//D
    sets=[[(d-r)//a for d in ds if (d-r)%a==0] for a,r,ds in zip(h,residues,E)]
    # Retain a witness, not only a Boolean hit.
    reachable={0:()}
    for ks in sets:
        reachable={old+k:word+(k,) for old,word in reachable.items() for k in ks}
    if wanted not in reachable:return None
    return tuple(r+a*k for a,r,k in zip(h,residues,reachable[wanted]))
A=(56,13,25,48,34,6); B=(52,13,24,49,34,6)
pa=c.profile(A);pb=c.profile(B)
ck(pa['H']==41 and pb['H']==43 and pa['h']==15,'seed values')
ck(all(pa[k]==pb[k] for k in pa if k!='H'),'physical central pair')
physical=[]
for X,Y in product((A,B),repeat=2):
    T=c.fuse(X,Y);H=literal_hp(T);M=c.mul(c.matrix(c.profile(X)),c.matrix(c.profile(Y)))
    ck(H==c.from_matrix(M)['H'],'physical product count')
    physical.append(H)
ck(sorted(set(physical))==list(range(min(physical),min(physical)+61,30)),'mesh15 progression')
# Session282 carry-loop pair, also independently replayed by wildcard lane.
X=(594,2020,2001,1991,1994,2013,1152,513,193,320,897)
Y=(1542,2020,2000,1991,1995,2013,1153,513,193,320,896)
px=c.profile(X);py=c.profile(Y)
ck(px['h']==169 and px['H']==2747 and py['H']==2749,'new central gate')
ck(all(px[k]==py[k] for k in px if k!='H'),'new whole packet')
coprime_counts=[]
for left,right in product((A,B),(X,Y)):
    H=literal_hp(c.fuse(left,right));coprime_counts.append(H)
ck(coprime_counts==[474053,474083,474391,474421],'literal coprime composition')
ck(gcd(169,15)==1 and 474055 not in coprime_counts,'primitive not interval')
for N in range(474051,474424,2):
    hit=carry_select((15,169),[(0,1),(0,1)],(N-474053)//2)
    ck((hit is not None)==(N in coprime_counts),'physical absolute target oracle')
# Exact integer arithmetic controls; no claim these synthetic banks are physical.
universes=0;targets=0
for hs in ((3,5),(3,7),(5,7),(3,5,7),(3,5,7,11)):
    for stretch in (1,2,3):
        ms=[stretch*(a-1) for a in hs]
        E=[list(range(m+1)) for m in ms];D=prod(hs);w=[D//a for a in hs]
        support={sum(x*b for x,b in zip(ds,w)) for ds in product(*E)}
        S=sum(m*b for m,b in zip(ms,w));C=(len(hs)-1)*D-sum(w)+1
        for z in range(-2,S+3):
            witness=carry_select(hs,E,z)
            ck((witness is not None)==(z in support),'absolute selector')
            if witness is not None:ck(sum(x*b for x,b in zip(witness,w))==z,'selected witness')
            if C<=z<=S-C:ck(z in support,'conductor interior')
            targets+=1
        universes+=1
for hs in ((3,5),(3,5,7)):
    E=[[-a,-1,0,1,a+2] for a in hs];D=prod(hs);w=[D//a for a in hs]
    support={sum(x*b for x,b in zip(ds,w)) for ds in product(*E)}
    for z in range(min(support)-1,max(support)+2):
        ck((carry_select(hs,E,z)is not None)==(z in support),'signed sparse digit bank')
        targets+=1
ck(carry_select((3,5),[range(3),range(5)],1)is None,'residue-only hostile')
valuation_checks=0
for r in range(1,5):
    for hs in product((1,3,5,9,15),repeat=r):
        D=prod(hs);w=[D//h for h in hs];G=reduce(gcd,w)
        ck(G==D//lcm(*hs),'prime ownership gcd')
        ck((G==1)==all(gcd(a,b)==1 for a,b in combinations(hs,2)),'primitive iff disjoint primes')
        valuation_checks+=1
print(json.dumps(dict(status='PASS',physical_pair_profiles=[pa,pb],literal_two_factor_counts=physical,physical_coprime_counts=coprime_counts,
 arithmetic_complete_bank_universes=universes,arithmetic_targets_tested=targets,
 prime_ownership_tuples=valuation_checks,arithmetic_not_physical=True),indent=2,sort_keys=True))
