#!/usr/bin/env python3
"""
LRC(14) covering-bound investigation v2 (kind-pasteur-2026-06-27). FAST.

Float-screen M(S) (candidate method), exact-verify only the minimizers.

Probes:
  P1. Reduction sanity: 14-free => M>=1/14 at t=1/14 (trivial half).
  P2. The covering MARGIN: min M over covering families. Is it > 1/14?
      (mac-mini-S59: 'the bound is strictly weaker than the census'.)
  P3. Optimized-gamma 14-lift sieve on |M14|<=6: does a Q-safe gamma + surviving
      lift always certify M>=1/14?  Tabulate forbidden lifts by R-speed type.
"""
from fractions import Fraction as F
from itertools import combinations
import random, math

def nrm_f(x):
    x = x - math.floor(x)
    return min(x, 1.0 - x)

def M_float(speeds):
    """Float M via the complete crossing/peak candidate set (fast screen)."""
    speeds = sorted(set(speeds))
    n = len(speeds)
    cand = set()
    for v in speeds:
        for k in range(v):                       # peaks (2k+1)/(2v) in (0,1/2]
            cand.add((2 * k + 1, 2 * v))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = speeds[i], speeds[j]
            for D in (a + b, b - a):
                if D <= 0:
                    continue
                for k in range(1, (D // 2) + 1):
                    cand.add((k, D))
    best = 0.0
    bt = None
    for num, den in cand:
        t = num / den
        if t <= 0 or t > 0.5:
            continue
        d = min(nrm_f(v * t) for v in speeds)
        if d > best:
            best = d; bt = (num, den)
    return best, bt

def nrm(x):
    x = x - int(x)
    if x < 0: x += 1
    return min(x, 1 - x)

def M_exact(speeds):
    speeds = sorted(set(speeds)); n = len(speeds); cand = set()
    for v in speeds:
        for k in range(v):
            cand.add(F(2 * k + 1, 2 * v))
    for i in range(n):
        for j in range(i + 1, n):
            a, b = speeds[i], speeds[j]
            for D in (a + b, b - a):
                if D <= 0: continue
                for k in range(1, (D // 2) + 1):
                    cand.add(F(k, D))
    best = F(0)
    for t in cand:
        if 0 < t <= F(1, 2):
            d = min(nrm(v * t) for v in speeds)
            if d > best: best = d
    return best

def is_cov(S): return any(v % 14 == 0 for v in S)
ONE14 = 1.0 / 14

print("=" * 70)
print("LRC(14) COVERING MARGIN + LIFT SIEVE  (kps S31 cont., v2)")
print("=" * 70)

# ---- P1 ----
print("\n[P1] 14-free => M>=1/14 at t=1/14 (trivial reduction half).")
random.seed(1); bad = tested = 0
for _ in range(4000):
    S = random.sample(range(1, 70), 13)
    if is_cov(S): continue
    tested += 1
    if min(nrm(v * F(1, 14)) for v in S) < F(1, 14): bad += 1
print(f"   {tested} 14-free sets; t=1/14 fails on {bad} (expect 0).")

# ---- P2: covering margin ----
print("\n[P2] Covering MARGIN: min M over covering families (one mult of 14 present).")
def scan(gen, label, exact_top=3):
    rows = []
    for S in gen:
        if len(set(S)) != 13 or not is_cov(S): continue
        m, bt = M_float(S)
        rows.append((m, tuple(sorted(S))))
    rows.sort()
    print(f"   {label}: {len(rows)} sets; smallest M (float):")
    out = []
    for m, S in rows[:exact_top]:
        me = M_exact(list(S))
        out.append((me, S))
        print(f"        M={me}={float(me):.5f}  S={S}")
    return out

# A: 13 consecutive containing a multiple of 14
A = scan((list(range(a, a + 13)) for a in range(1, 160) if is_cov(list(range(a, a + 13)))),
         "consecutive-13 (covering)")
# B: AP {1..13}, replace one entry by a multiple of 14
def apswap():
    base = list(range(1, 14))
    for d in range(13):
        for mult in (14, 28, 42, 56, 70, 84):
            yield base[:d] + base[d+1:] + [mult]
B = scan(apswap(), "AP{1..13} one->mult14")
# C: random 12-subsets of {1..26}\{14} plus 14  (bounded comparable core)
def Csample(k=6000):
    pool = [v for v in range(1, 27) if v != 14]; rng = random.Random(5)
    for _ in range(k): yield rng.sample(pool, 12) + [14]
C = scan(Csample(), "rand 12-subset {1..26}\\{14} + {14}")
# D: dilations / structured: AP times small factor plus a 14 (look for near-tight covering)
def Dsample(k=6000):
    rng = random.Random(9)
    for _ in range(k):
        base = rng.sample(range(1, 30), 12)
        yield base + [14 * rng.randint(1, 4)]
D = scan(Dsample(), "rand 12 from {1..29} + small mult14")

allmins = [r[0][0] for r in (A, B, C, D) if r]
gmin = min(allmins)
print(f"\n   GLOBAL min M over covering sets tried: {gmin} = {float(gmin):.5f}")
print(f"   1/14={ONE14:.5f}  1/13={1/13:.5f}  1/12={1/12:.5f}  1/8={1/8:.5f}")
if gmin > F(1, 14):
    print(f"   => covering min is STRICTLY > 1/14.  Margin = {gmin - F(1,14)} = {float(gmin-F(1,14)):.5f}")
    print(f"      >= 1/13 ? {gmin >= F(1,13)}    >= 1/12 ? {gmin >= F(1,12)}")
else:
    print("   => a covering set reaches 1/14 (bound is tight on covering sets).")

# ---- P3: optimized-gamma 14-lift sieve on |M14|<=6 ----
print("\n[P3] Optimized-gamma 14-lift sieve (|M14|<=6 residual).")
def g14(b): return math.gcd(b, 14)
def sieve(S, gden=2800):
    Q = [v // 14 for v in S if v % 14 == 0]
    R = [v for v in S if v % 14 != 0]
    tc = {1: 0, 2: 0, 7: 0}
    for b in R: tc[g14(b)] += 1
    for gi in range(1, gden):
        gamma = F(gi, gden)
        if any(nrm(q * gamma) < F(1, 14) for q in Q): continue
        for m in range(14):
            t = (gamma + m) / 14
            if all(nrm(b * t) >= F(1, 14) for b in R):
                return True, gamma, m, tc
    return False, None, None, tc

tests = []
for r in (A, B, C, D):
    if r: tests.append(r[0][1])
tests += [tuple(sorted([14] + list(range(1, 13)))),
          tuple(sorted([14, 28] + list(range(1, 12)))),
          tuple(sorted([14, 28, 42] + list(range(1, 11)))),
          tuple(range(2, 15))]
fails = 0
for S in sorted(set(tests)):
    nM = sum(1 for v in S if v % 14 == 0)
    ok, g, m, tc = sieve(list(S))
    me = M_exact(list(S))
    if not ok: fails += 1
    print(f"   |M14|={nM} Rtype(1/2/7)={tc[1]}/{tc[2]}/{tc[7]} M={float(me):.4f} "
          f"{'OK g='+str(g)+' m='+str(m) if ok else 'SIEVE-FAILS'}  S={S}")
print(f"\n   sieve failures: {fails}/{len(set(tests))} "
      "(fail != M<1/14; only means this certificate didn't fire)")
print("=" * 70); print("DONE")
