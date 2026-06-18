"""
Part 4: INDEPENDENTLY re-enumerate the k=2 finite hard core (all speeds <= 62)
and verify: count, worst M, 0 sets below 1/14.

Definitions (from prompt case-split):
  - A covering 13-set S: |S|=13 distinct positive ints, primitive (gcd 1),
    contains a multiple of every q in 2..14 (THM-523 covering condition).
  - k = #{v in S : v > 13}.
  - Vmin = min S, Vmax = max S.
  - S3 = k>=2 AND Vmax >= 13*Vmin.
  - k=2 slice: exactly 2 speeds > 13, so 11 speeds in {1..13}.
  - hard core: all speeds <= 62.

So a k=2 hard-core set = (11-subset of {1..13}) U {two speeds in 14..62},
primitive, covering, with Vmax >= 13*Vmin.

Vmin: with 11 speeds in {1..13}, Vmin = min small speed (likely 1). Then
13*Vmin could be 13 -> Vmax>=13*1=13 is automatic if a large speed exists. But
if Vmin>1 (e.g. small part excludes 1), Vmin>=2 -> 13*Vmin>=26. We must apply
the S3 condition Vmax>=13*Vmin and the S2-exclusion (NOT clustered). Let's just
enumerate ALL k=2 covering primitive 13-sets with speeds<=62 and classify by S3.

We then M-check every S3 one. We report:
  - total k=2 covering primitive sets, speeds<=62
  - how many are S3 (Vmax>=13*Vmin)
  - worst M among S3, count below 1/14
We compare count to claimed 4865.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, time

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x * t) for x in S)
        if v > b: b = v
    return b

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def is_primitive(S):
    return reduce(gcd, S) == 1

print("="*78)
print("PART 4: independent re-enumeration of k=2 finite hard core (speeds<=62).")
print("="*78)
t0=time.time()
# 11-subset of {1..13}, then 2 large speeds in 14..62
small_parts = [list(c) for c in itertools.combinations(range(1,14), 11)]
large_pairs = list(itertools.combinations(range(14,63), 2))
print(f"  #small parts (11-subset of 1..13) = {len(small_parts)}")
print(f"  #large pairs (2-subset of 14..62) = {len(large_pairs)}")

total_cov=0
s3_sets=[]
for P in small_parts:
    Vmin = min(P)
    for (a,b) in large_pairs:
        S = P + [a,b]
        # primitive
        if reduce(gcd,S)!=1: continue
        if not is_covering(S): continue
        total_cov += 1
        Vmax = b
        # k: speeds>13. small part all <=13, a,b>=14 so k=2 exactly. good.
        # S3 condition: Vmax >= 13*Vmin
        if Vmax >= 13*Vmin:
            s3_sets.append(S)
print(f"  total k=2 covering primitive sets (speeds<=62): {total_cov}  ({time.time()-t0:.0f}s)")
print(f"  of those, S3 (Vmax>=13*Vmin): {len(s3_sets)}")
print(f"  claimed hard-core count: 4865")

# M-check the S3 ones
worst=None; worst_arg=None; below=0; below_list=[]
n=0
for S in s3_sets:
    n+=1
    m=Mval(S)
    if worst is None or m<worst: worst=m; worst_arg=S
    if m < F(1,14):
        below+=1; below_list.append((S,m))
    if n % 1000==0:
        print(f"    ...{n}/{len(s3_sets)} checked, worst so far {worst} ({time.time()-t0:.0f}s)")
print(f"  WORST M among S3 hard core = {worst} = {float(worst):.6f}  at S={worst_arg}")
print(f"  # with M < 1/14 : {below}")
for S,m in below_list[:20]:
    print(f"    COUNTEREXAMPLE S={S} M={m}={float(m):.6f}")
print(f"  ALL k=2 S3 hard-core sets have M>=1/14 ? {below==0}")
print(f"  ({time.time()-t0:.0f}s)")
print("DONE PART 4.")
