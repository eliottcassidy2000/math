"""
Part 6: BROAD adversarial hunt -- ignore the partition, directly search for ANY
covering primitive S3 set with M < 1/14. Two regimes:

 (A) k=2 with Vmax UP TO some large bound (beyond 62): verify the >=63 branch's
     promise empirically -- no k=2 covering S3 set with large Vmax has M<1/14, and
     in fact C fires via Vmax. We test C explicitly (W(S\{Vmax}) > 1/(7 Vmax)).
 (B) k>=3 sets: the self-declared OPEN regime. Hunt for M<1/14 across structured
     adversarial families (small cluster spread, comparable large speeds, varied
     offsets), AND test whether C fires (to confirm the gap is 'unproved' not 'false').

We use exact Mval and exact W. We also report, for k>=3 sets, how often the
single-tooth criterion C fails (the prompt claimed ~2.8%).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, time, random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h / u) % 1; b = (c + h / u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

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

def C_fires(S):
    S=sorted(S)
    for v in S:
        A=[x for x in S if x!=v]
        if Wwidth(A) > F(1,7*v):
            return True
    return False

def C_fires_via_max(S):
    S=sorted(S); v=S[-1]; A=S[:-1]
    return Wwidth(A) > F(1,7*v)

random.seed(99)
print("="*78)
print("PART 6A: k=2 covering S3 sets with Vmax in 63..LARGE -- confirm C fires via")
print("  Vmax and M>=1/14 (spot-check the >=63 branch beyond the hard core).")
print("="*78)
t0=time.time()
small_parts = [list(c) for c in itertools.combinations(range(1,14), 11)]
fail_M=[]; fail_C=[]; n=0
# sample large pairs with Vmax >= 63
for trial in range(6000):
    P = random.choice(small_parts)
    Vmin=min(P)
    a = random.randint(14, 400)
    b = random.randint(max(a+1,63), 600)
    S = sorted(P+[a,b])
    if reduce(gcd,S)!=1: continue
    if not is_covering(S): continue
    if b < 13*Vmin: continue   # must be S3
    n+=1
    m=Mval(S)
    if m < F(1,14): fail_M.append((S,m))
    if not C_fires_via_max(S): fail_C.append(S)
print(f"  k=2 S3 (Vmax>=63) sampled: {n}")
print(f"  # with M<1/14: {len(fail_M)}")
for S,m in fail_M[:10]: print(f"    M-FAIL S={S} M={m}")
print(f"  # where C does NOT fire via Vmax: {len(fail_C)}")
for S in fail_C[:10]: print(f"    C-via-max FAIL S={S}")
print(f"  ({time.time()-t0:.0f}s)")

print("\n"+"="*78)
print("PART 6B: k>=3 covering S3 sets -- hunt M<1/14, measure C-failure rate.")
print("="*78)
t1=time.time()
failM3=[]; cfail=0; n3=0
# structured: small part subset of 1..13 (size 13-k) + cluster of k large speeds
for trial in range(40000):
    k = random.randint(3, 6)
    ssize = 13-k
    if ssize < 1: continue
    P = sorted(random.sample(range(1,14), ssize))
    # cluster of k large speeds near some V0, with small spread (adversarial)
    V0 = random.randint(14, 250)
    spread = random.choice([5,10,14,20,30,45,60])
    L = set()
    while len(L) < k:
        L.add(V0 + random.randint(0, spread))
    L = sorted(L)
    S = sorted(set(P) | set(L))
    if len(S)!=13: continue
    if reduce(gcd,S)!=1: continue
    if not is_covering(S): continue
    Vmin=min(S); Vmax=max(S)
    if not (sum(1 for v in S if v>13)>=2 and Vmax>=13*Vmin): continue  # S3, k>=2
    if sum(1 for v in S if v>13) < 3: continue  # ensure k>=3
    n3+=1
    m=Mval(S)
    if m < F(1,14): failM3.append((S,m))
    if not C_fires(S): cfail+=1
print(f"  k>=3 S3 sampled: {n3}")
print(f"  # with M<1/14: {len(failM3)}")
for S,m in failM3[:15]: print(f"    *** M-FAIL S={S} M={m}={float(m):.5f}")
print(f"  # where criterion C does NOT fire (any v): {cfail}  ({100*cfail/max(n3,1):.2f}%)")
print(f"  ({time.time()-t1:.0f}s)")
print("DONE PART 6.")
