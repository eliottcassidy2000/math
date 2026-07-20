import itertools, random
import numpy as np
from collections import defaultdict
def scores(S,n): return [int(sum(1 for j in range(n) if j!=i and S[i,j]==1)) for i in range(n)]
def switch(S,n,W):
    T=S.copy()
    for i in range(n):
        for j in range(n):
            if i!=j and (((W>>i)&1)^((W>>j)&1)): T[i,j]=-S[i,j]
    return T
def rnd(n,rng):
    S=np.zeros((n,n),dtype=np.int64)
    for i in range(n):
        for j in range(i+1,n):
            v=1 if rng.random()<.5 else -1
            S[i,j]=v; S[j,i]=-v
    return S
print("THE PARITY ARGUMENT, WORKED OUT:")
print("  switching at W changes s_v by |W^c|-2a (v in W) or |W|-2b (v not in W),")
print("  so the score-PARITY vector flips on W iff |W^c| odd, off W iff |W| odd.")
print("  n ODD: exactly one of |W|,|W^c| is odd, and W ~ W^c, so WLOG |W| even and the")
print("         parity vector flips exactly on W.  W ranges over the 2^(n-1) EVEN subsets,")
print("         and the class has 2^(n-1) members => PARITY VECTOR IS A COMPLETE INVARIANT.")
print("  #odd scores = C(n,2) mod 2, so the reachable coset is:")
print("         n=1 mod 4 -> EVEN-weight coset, contains 0        -> unique ALL-EVEN member")
print("         n=3 mod 4 -> ODD-weight coset, contains all-ones  -> unique ALL-ODD member")
print("  BOTH anchors are permutation-invariant, so BOTH kill Babai-Cameron.")
print("  n EVEN: |W|,|W^c| have the SAME parity -> either NO vertex flips or ALL do,")
print("         so only 2 parity vectors occur and the argument COLLAPSES.\n")
rng=random.Random(7)
for n in [3,5,7,9,11]:
    dist_e=defaultdict(int); dist_o=defaultdict(int)
    trials = 40 if n>7 else 200
    for _ in range(trials):
        S=rnd(n,rng); e=o=0
        for W in range(1<<n):
            T=switch(S,n,W); sc=scores(T,n)
            if all(s%2==0 for s in sc): e+=1
            if all(s%2==1 for s in sc): o+=1
        dist_e[e//2]+=1; dist_o[o//2]+=1
    print(f"  n={n:2d} ({n%4} mod 4): all-EVEN members/class {dict(sorted(dist_e.items()))}"
          f"   all-ODD members/class {dict(sorted(dist_o.items()))}")
print()
print("  n EVEN control (argument should fail):")
for n in [4,6,8]:
    dist_e=defaultdict(int)
    for _ in range(60):
        S=rnd(n,rng); e=0
        for W in range(1<<n):
            if all(s%2==0 for s in scores(switch(S,n,W),n)): e+=1
        dist_e[e//2]+=1
    print(f"  n={n:2d}: all-EVEN members per class {dict(sorted(dist_e.items()))}"
          f"  <- NOT always 1 => no canonical anchor")
