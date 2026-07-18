#!/usr/bin/env python3
"""
Does COMPACT covering => M >= 1/13 ?   (boxeph-2026-07-18-S85)
=============================================================
The deep well (M=14/183, the global covering-min, THM-724) is NON-compact (rho=15).
Every compact (rho<13) primitive covering family seen so far has M >= 1/13 > 1/14.
CONJECTURE: primitive covering with rho = v_max/v_2nd < 13  =>  M >= 1/13.
If true, this CLOSES the compact residual (1/13 > 1/14) with a clean FLOOR.
Hunt hard for a counterexample: random, structured (dilated-AP+killer), and
adversarial near-tight compact covering families.  Report the minimum M and any
family with M < 1/13.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, itertools

def exact_M(V):
    if len(V) == 1: return F(1,2)
    best = F(0); qs = set([13,14])
    for i in range(len(V)):
        for j in range(i, len(V)):
            s = V[i]+V[j]
            for d in range(1, s+1):
                if s % d == 0: qs.add(d)
    for q in qs:
        for a in range(1, q):
            if gcd(a,q)==1:
                m = min(min((v*a)%q, q-(v*a)%q) for v in V)
                c = F(m,q)
                if c>best: best=c
    return best

def cover(V,n=14): return all(any(v%b==0 for v in V) for b in range(2,n+1))
def prim(V): return reduce(gcd,V)==1
def rho(V):
    Vs=sorted(V); return F(Vs[-1],Vs[-2])

THRESH = F(1,13)
below = []   # families with M < 1/13
minM = F(1); minV = None
count = 0

def consider(V):
    global minM, minV, count
    V = sorted(set(V))
    if len(V)!=13: return
    if not prim(V) or not cover(V): return
    if rho(V) >= 13: return   # compact only
    count += 1
    M = exact_M(V)
    if M < minM:
        minM = M; minV = V
    if M < THRESH:
        below.append((M, list(V)))

# --- 1. random compact covering ---
random.seed(20)
for _ in range(200000):
    V = random.sample(range(1,70), 13)
    consider(V)
    if count > 120000: break

# --- 2. structured: d*{1..12} u {killer}, near-dilated-tight ---
for d in range(1,10):
    ap = [d*i for i in range(1,13)]
    for w in range(1, 500):
        consider(ap+[w])

# --- 3. d*{1..13} with one element swapped (make primitive, stay covering) ---
for d in range(2,8):
    ap = [d*i for i in range(1,14)]
    for idx in range(13):
        for repl in range(1, 120):
            consider(ap[:idx]+ap[idx+1:]+[repl])

# --- 4. adversarial: perturb the tight AP {1..13} toward covering+primitive ---
tight = list(range(1,14))
for a in range(1,14):
    for b in range(14, 120):
        consider(tight[:a-1]+tight[a:]+[b])  # replace element a by b
# double replacement
for a1,a2 in itertools.combinations(range(1,14),2):
    for b in range(14, 60):
        base = [x for x in tight if x!=a1 and x!=a2]
        consider(base+[b, b+13])

# --- 5. scaled tight + coprime shifts ---
for d in [2,3,4,5,6]:
    for shift in range(1, 40):
        V = [d*i for i in range(1,13)] + [d*13 + shift]
        consider(V)
        V2 = [d*i for i in range(1,13)] + [shift]
        consider(V2)

print("="*90)
print("HUNT: compact (rho<13) primitive covering families with M < 1/13")
print("="*90)
print(f"tested {count} compact primitive covering families")
print(f"MINIMUM M found: {minM} = {float(minM):.6f}")
print(f"  1/14={float(F(1,14)):.6f}  1/13={float(F(1,13)):.6f}  14/183={float(F(14,183)):.6f}")
print(f"families with M < 1/13: {len(below)}")
if below:
    below.sort(key=lambda t:t[0])
    print("  SMALLEST (potential floor-breakers):")
    for M,V in below[:15]:
        print(f"    M={M} ({float(M):.6f})  rho={float(rho(V)):.2f}  V={V}")
    print("\n  => CONJECTURE 'compact covering => M>=1/13' is FALSE; floor is lower.")
else:
    print("  NONE. => CONJECTURE 'compact covering => M >= 1/13' SURVIVES the hunt.")
    print(f"  The minimizer sits at M={minM}, V={minV}")
