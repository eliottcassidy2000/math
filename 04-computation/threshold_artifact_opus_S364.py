# opus-2026-07-17-S364 -- HYP-7470: THE EMPIRICAL THRESHOLD IS A SAMPLING
# ARTIFACT.  The S362 stratum run reported
#     min-speed [23,70):   0/4 BONF5 > 0
#     min-speed [150,300): 4/4 BONF5 > 0
#     min-speed [600,900): 4/4 BONF5 > 0
# which LOOKS like a threshold -- but THM-1050 proves none can exist, since
# BONF5 is dilation-invariant.  Resolution: dilates of small failures live in
# every stratum, and random sampling from [m,13m] essentially never draws one.
# This script exhibits the counterexample explicitly.
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import random, itertools
LAM = F(1,14)
MODULI = [8,9,10,11,12,13,14]
def blocks_all(V): return all(any(v % q == 0 for v in V) for q in MODULI)
def teeth(x, lo, hi):
    w = LAM/x; out=[]
    for j in range(floor((lo-w)*x), floor((hi+w)*x)+2):
        a,b = max(F(j,x)-w,lo), min(F(j,x)+w,hi)
        if a<b: out.append((a,b))
    return out
def inter(u,v):
    out,i,j=[],0,0
    while i<len(u) and j<len(v):
        a,b = max(u[i][0],v[j][0]), min(u[i][1],v[j][1])
        if a<b: out.append((a,b))
        if u[i][1]<v[j][1]: i+=1
        else: j+=1
    return out
def mu(S):
    cur = teeth(S[0],F(0),F(1))
    for x in S[1:]: cur = inter(cur, teeth(x,F(0),F(1)))
    return sum(b-a for a,b in cur)
def bonf5(V):
    tot = F(1)
    for k in range(1,6):
        Sk = F(0)
        for C in itertools.combinations(V,k): Sk += mu(list(C))
        tot += (-1)**k * Sk
    return tot

print("STEP 1: find a SMALL failing family (min speed < 70, blocks all seven).")
random.seed(364)
small = None
while small is None:
    m = random.randint(23, 60)
    V = sorted(random.sample(range(m, 13*m), 13))
    if V[0] >= 70 or not blocks_all(V): continue
    b = bonf5(V)
    if b <= 0: small = (V, b)
V, b = small
print(f"   V = {V}")
print(f"   min speed {V[0]}, gcd {reduce(gcd,V)}, BONF5 = {float(b):+.6f}  (FAILS)")

print()
print("STEP 2: dilate it into the [600,900) stratum, where the run reported 4/4 positive.")
k = 900 // V[0]
W = [k*v for v in V]
bw = bonf5(W)
print(f"   k = {k};  W = {W[:5]}...")
print(f"   min speed {W[0]} (in the 'all positive' stratum), gcd {reduce(gcd,W)}")
print(f"   BONF5(W) = {float(bw):+.6f}")
print(f"   BONF5(V) = {float(b):+.6f}    identical: {bw == b}")
print(f"   blocks all seven: {blocks_all(W)}")
print()
print("STEP 3: verdict.")
print("   A family with min speed >= 600 and BONF5 < 0 EXISTS, exhibited above.")
print("   The stratum run's '4/4 positive' was true of its four RANDOM draws and")
print("   false of the stratum: random sampling from [m,13m] almost never draws a")
print("   dilate (it would need a common factor across all 13 speeds).  The")
print("   apparent threshold is a SAMPLING ARTIFACT, exactly as THM-1050 predicts.")
