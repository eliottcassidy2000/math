#!/usr/bin/env python3
"""
death-star-2026-07-20-S61d (HYP claim below) -- the owner's synthesis:
(A) CONFIRM boxeph's 3/8 mass at K8: E[eps over S_n], eps(s)=(-1)^{sum_{even cycles} len/2}.
(B) BABAI-CAMERON Remark 7.4: at n=1 mod 4 each switching class has a UNIQUE all-even-out-degree
    member (=> every automorphism fixes it => fixed-point-free-automorphism count = 0);
    at n=3 mod 4 NO even member exists (C(n,2) odd). The odd-n mod-4 dichotomy.
(C) 3,7 (=3 mod4, Paley TOURNAMENTS, no even member) vs 5,9 (=1 mod4, Paley GRAPHS, even member).
"""
from itertools import permutations, combinations
from fractions import Fraction as Fr
from math import comb, factorial

# ---------- (A) E[eps over S_n] = sum_{partitions} eps(lambda)/z_lambda ----------
def partitions(n, mx=None):
    if mx is None: mx=n
    if n==0: yield []; return
    for k in range(min(n,mx),0,-1):
        for rest in partitions(n-k,k): yield [k]+rest
def z_lambda(lam):
    from collections import Counter
    c=Counter(lam); z=1
    for part,m in c.items(): z*= (part**m)*factorial(m)
    return z
def eps_of(lam):
    s=sum(p//2 for p in lam if p%2==0)
    return (-1)**s
print("=== (A) E[eps over S_n]  (mass at the complete-core K_n class); confirm 3/8 at n=8 ===")
seq=[]
for n in range(2,11):
    E=sum(Fr(eps_of(lam), z_lambda(lam)) for lam in partitions(n))
    seq.append(E)
    fplus=(1+E)/2
    print(f"  n={n}: E[eps]={E}  f+=(1+E)/2={fplus}  mass|E|={abs(E)}  {'<-- 3/8 !' if abs(E)==Fr(3,8) else ''}")
print(f"  sequence n=2..10: {[str(e) for e in seq]}  (quarter-law {{0,1/2}} breaks to 3/8 at n=8,9)")

# ---------- (B)+(C) even-member dichotomy by n mod 4 ----------
print("\n=== (B) all-even-out-degree member per switching class (Babai-Cameron) ===")
def switch_v(bits,v,n,idx):
    for u in range(n):
        if u!=v: bits^=(1<<idx[(min(u,v),max(u,v))])
    return bits
def outdegs(bits,n,edges):
    od=[0]*n
    for k,(i,j) in enumerate(edges):
        if (bits>>k)&1: od[i]+=1
        else: od[j]+=1
    return od
def all_even(bits,n,edges): return all(d%2==0 for d in outdegs(bits,n,edges))

for n in [3,5]:                      # exhaustive: n=3 (=3 mod4), n=5 (=1 mod4)
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(edges); idx={e:k for k,e in enumerate(edges)}
    seen=set(); per_class=[]
    for start in range(1<<m):
        if start in seen: continue
        cls=set([start]); fr=[start]
        while fr:
            b=fr.pop()
            for v in range(n):
                nb=switch_v(b,v,n,idx)
                if nb not in cls: cls.add(nb); fr.append(nb)
        seen|=cls
        per_class.append(sum(1 for b in cls if all_even(b,n,edges)))
    from collections import Counter
    print(f"  n={n} ({n%4} mod 4): C(n,2)={comb(n,2)} ({'even' if comb(n,2)%2==0 else 'ODD'}); "
          f"even-members-per-class distribution = {dict(Counter(per_class))}  "
          f"({'UNIQUE even member' if set(per_class)=={1} else 'NO even member' if set(per_class)=={0} else 'mixed'})")

print("\n=== (C) the mod-4 dichotomy for ODD n (C(n,2) parity) ===")
for n in [3,5,7,9,11,13]:
    c=comb(n,2)
    paley = "Paley TOURNAMENT (q=3 mod4)" if n%4==3 else "Paley GRAPH self-compl (q=1 mod4)"
    member = "NO even member" if c%2==1 else "UNIQUE even member -> BC count 0"
    print(f"  n={n} = {n%4} mod4: C(n,2)={c} {'odd ' if c%2 else 'even'} => {member:32s} | {paley}")
print("\n  => 3,7,11 (=3 mod4): odd C(n,2), NO even member, Paley TOURNAMENTS -- the {7,21}-gap side.")
print("     5,9,13 (=1 mod4): even C(n,2), UNIQUE even member (BC=0), Paley GRAPHS/self-complementary.")
print("     This is why primes 3 and 7 pattern together (mod 4) and 5 sits with 1,9.")
