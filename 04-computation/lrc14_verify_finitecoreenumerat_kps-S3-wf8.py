"""
Part 8: targeted k>=3 'comparable large speeds' adversarial probe + the SHARPEST
possible counterexample hunt. Two goals:
 (1) CONFIRM the k>=3 gap is real: find k>=3 covering S3 sets where the single-tooth
     criterion C(S) (via v=Vmax, the dropped-max approach) FAILS, i.e.
     W(S\{Vmax}) <= 1/(7 Vmax). This shows the k=2 proof technique does NOT extend.
     (Even if M(S)>=1/14 still holds, the PROOF route breaks -> k>=3 unproved.)
 (2) Push hard for an actual M<1/14 counterexample in k>=3 with two/three nearly
     equal large speeds (the Moire/cluster-collapse regime).

We build sets: small part (subset of 1..13) + a TIGHT cluster of k large speeds
all within a small window. Covering enforced. We test:
   - C_via_max: W(S minus Vmax) > 1/(7 Vmax)
   - C_any: exists v with W(S minus v) > 1/(7 v)
   - M(S) >= 1/14
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

def C_via_max(S):
    S=sorted(S); v=S[-1]; A=S[:-1]
    return Wwidth(A) > F(1,7*v)

def C_any(S):
    S=sorted(S)
    for v in S:
        A=[x for x in S if x!=v]
        if Wwidth(A) > F(1,7*v): return True
    return False

random.seed(7)
print("="*78)
print("PART 8: k>=3 comparable-large-speeds probe. Drop-max C-failure & M-hunt.")
print("="*78)
t0=time.time()
cmax_fail=0; cany_fail=0; mfail=[]; n=0
examples_cmaxfail=[]
for trial in range(60000):
    k = random.randint(3,5)
    ssize = 13-k
    P = sorted(random.sample(range(1,14), ssize))
    # tight cluster: k large speeds within window of width <= 45, near V0
    V0 = random.randint(20, 300)
    win = random.choice([14,20,28,35,45])
    L=set()
    while len(L)<k:
        L.add(V0+random.randint(0,win))
    L=sorted(L)
    if max(L)<2*V0//1 and (max(L)-min(L))> 0:  # keep them comparable (ratio<2)
        pass
    S=sorted(set(P)|set(L))
    if len(S)!=13: continue
    if reduce(gcd,S)!=1: continue
    if not is_covering(S): continue
    Vmin=min(S); Vmax=max(S)
    klarge=sum(1 for v in S if v>13)
    if klarge<3: continue
    if Vmax < 13*Vmin: continue  # S3
    n+=1
    if not C_via_max(S):
        cmax_fail+=1
        if len(examples_cmaxfail)<8: examples_cmaxfail.append(S)
    if not C_any(S):
        cany_fail+=1
    m=Mval(S)
    if m < F(1,14): mfail.append((S,m))
print(f"  k>=3 S3 (comparable cluster) sampled: {n}  ({time.time()-t0:.0f}s)")
print(f"  drop-max criterion C-via-Vmax FAILS: {cmax_fail}  ({100*cmax_fail/max(n,1):.2f}%)")
for S in examples_cmaxfail:
    print(f"    C-via-max FAIL: S={S}  W(S\\max)={Wwidth(sorted(S)[:-1])}  thr=1/{7*max(S)}  M={Mval(S)}")
print(f"  criterion C (ANY v) fails: {cany_fail}  ({100*cany_fail/max(n,1):.2f}%)")
print(f"  *** actual M<1/14 counterexamples: {len(mfail)}")
for S,m in mfail[:15]:
    print(f"      M-FAIL S={S} M={m}={float(m):.5f}")
print(f"  ({time.time()-t0:.0f}s)")
print("DONE PART 8.")
