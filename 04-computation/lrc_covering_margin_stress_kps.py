#!/usr/bin/env python3
"""
LRC(14) covering-bound STRESS TEST (kind-pasteur-2026-06-27).

Two conjectures emerged from v2:
  (C1) covering => M >= 1/13   [margin, min attained uniquely at {1..12,14}]
  (C2) optimized-gamma 14-lift sieve fires on every |M14|<=6 covering set

Stress both to find the breaking point (if any).
  S1: hunt covering sets with M < 1/13 (GW-type jumps, dilations, wide ranges).
  S2: optimized-gamma lift sieve on many random |M14|<=6 sets, esp. multi-7-type R.
  S3: the aliasing explanation: for LRC(13)-tight R (={1..12}*d), which multiples
      of 14 stay >= 1/13?  (14 == 1 mod 13 aliasing, fails iff 13 | k.)
"""
from fractions import Fraction as F
from itertools import combinations
import random, math

def nrm(x):
    x = x - int(x)
    if x < 0: x += 1
    return min(x, 1 - x)
def nrm_f(x):
    x = x - math.floor(x); return min(x, 1.0 - x)

def M_float(speeds):
    speeds = sorted(set(speeds)); n = len(speeds); cand = set()
    for v in speeds:
        for k in range(v): cand.add((2*k+1, 2*v))
    for i in range(n):
        for j in range(i+1, n):
            a, b = speeds[i], speeds[j]
            for D in (a+b, b-a):
                if D <= 0: continue
                for k in range(1, D//2 + 1): cand.add((k, D))
    best = 0.0
    for num, den in cand:
        t = num/den
        if 0 < t <= 0.5:
            d = min(nrm_f(v*t) for v in speeds)
            if d > best: best = d
    return best

def M_exact(speeds):
    speeds = sorted(set(speeds)); n = len(speeds); cand = set()
    for v in speeds:
        for k in range(v): cand.add(F(2*k+1, 2*v))
    for i in range(n):
        for j in range(i+1, n):
            a, b = speeds[i], speeds[j]
            for D in (a+b, b-a):
                if D <= 0: continue
                for k in range(1, D//2 + 1): cand.add(F(k, D))
    best = F(0)
    for t in cand:
        if 0 < t <= F(1,2):
            d = min(nrm(v*t) for v in speeds)
            if d > best: best = d
    return best

def is_cov(S): return any(v % 14 == 0 for v in S)
T13 = 1.0/13

print("="*70); print("LRC(14) COVERING STRESS  (kps S31 cont.)"); print("="*70)

# ---- S1: hunt covering M < 1/13 ----
print("\n[S1] Hunt covering 13-sets with M < 1/13 = %.5f (would refute C1)." % T13)
random.seed(2)
below = []
lo = 1.0
lo_set = None
N = 0
def consider(S):
    global lo, lo_set, N
    S = tuple(sorted(set(S)))
    if len(S) != 13 or not is_cov(S): return
    N += 1
    m = M_float(list(S))
    if m < lo - 1e-12:
        lo = m; lo_set = S
    if m < T13 - 1e-9:
        below.append((m, S))

# (a) GW-type: {1..11,13, X} with a jump X, made covering by forcing a mult of 14
for X in range(14, 60):
    consider(list(range(1,12)) + [13, X])         # GW family (X may be mult of 14)
# (b) AP{1..13} dilated by d, then make covering by replacing one entry with nearest mult of 14
for d in range(1, 6):
    base = [d*i for i in range(1,14)]
    for j in range(13):
        for mult in range(14, 14*40, 14):
            S = base[:j] + base[j+1:] + [mult]
            consider(S)
# (c) random wide-range covering
for _ in range(60000):
    k = random.choice([14,28,42,56,70,84,98])
    rest = random.sample([v for v in range(1, 90) if v % 14 != 0], 12)
    consider(rest + [k])
# (d) near-LRC13-tight cores (consecutive-ish) + a multiple of 14
for a in range(1, 30):
    for L in range(10, 13):
        core = list(range(a, a+L))
        for mult in range(14, 14*30, 14):
            if mult in core: continue
            need = 13 - len(core) - 1
            if need < 0: continue
            extra = [v for v in range(a+L, a+L+need+3) if v % 14 != 0][:need]
            S = core + extra + [mult]
            consider(S)

print(f"   examined {N} covering sets.")
print(f"   global min M (float) = {lo:.6f}  at {lo_set}")
if below:
    print(f"   *** {len(below)} sets with M < 1/13 FOUND (C1 refuted). Examples:")
    for m, S in sorted(below)[:6]:
        print(f"       M={M_exact(list(S))}={m:.5f}  {S}")
else:
    me = M_exact(list(lo_set))
    print(f"   NO covering set with M < 1/13 found.  min is exact = {me} = {float(me):.6f}")
    print(f"   => C1 ('covering => M >= 1/13') holds on this search; min attained at {lo_set}")

# ---- S2: optimized-gamma lift sieve failure boundary ----
print("\n[S2] Optimized-gamma 14-lift sieve on random |M14|<=6 sets (find failures).")
def g14(b): return math.gcd(b,14)
def sieve(S, gden=4200):
    Q = [v//14 for v in S if v % 14 == 0]
    R = [v for v in S if v % 14 != 0]
    tc = {1:0,2:0,7:0}
    for b in R: tc[g14(b)] += 1
    for gi in range(1, gden):
        gamma = F(gi, gden)
        if any(nrm(q*gamma) < F(1,14) for q in Q): continue
        for m in range(14):
            t = (gamma+m)/14
            if all(nrm(b*t) >= F(1,14) for b in R): return True, tc
    return False, tc

random.seed(3)
fails = []
by_n7 = {}   # failure rate by number of 7-type R-speeds
tested = 0
for _ in range(4000):
    nM = random.randint(1, 6)
    Q = random.sample([14*j for j in range(1, 7)], nM)
    R = random.sample([v for v in range(1, 70) if v % 14 != 0], 13 - nM)
    S = sorted(set(Q + R))
    if len(S) != 13: continue
    tested += 1
    ok, tc = sieve(S)
    by_n7.setdefault(tc[7], [0,0])
    by_n7[tc[7]][0] += 1
    if not ok:
        by_n7[tc[7]][1] += 1
        m = M_exact(S)
        fails.append((tc, float(m), tuple(S)))
print(f"   tested {tested} random |M14|<=6 covering sets; sieve failures: {len(fails)}")
print("   failure rate by # of 7-type R-speeds (gcd(b,14)=7):")
for n7 in sorted(by_n7):
    tot, fl = by_n7[n7]
    print(f"      n7={n7}: {fl}/{tot} failed ({100*fl/tot:.1f}%)")
if fails:
    print("   example sieve-failures (M shown -- note M may still be >=1/14):")
    for tc, m, S in fails[:6]:
        print(f"      Rtype1/2/7={tc[1]}/{tc[2]}/{tc[7]}  M={m:.4f}  S={S}")

# ---- S3: the aliasing explanation ----
print("\n[S3] Aliasing: R = {1..12}*d (LRC13-tight, M=1/13 at t=k/(13d)).")
print("     A multiple 14*j stays >= 1/13 unless it aliases to 0 mod 13.")
for d in (1, 2, 3):
    R = [d*i for i in range(1, 13)]
    mR = M_exact(R)
    print(f"   d={d}: R={R}, M(R)={mR}={float(mR):.5f}")
    safe = []
    unsafe = []
    for j in range(1, 14):
        mult = 14*j
        if mult in R: continue
        m = M_exact(R + [mult])
        (safe if m >= F(1,13) else unsafe).append((j, mult, float(m)))
    print(f"      14*j with M(R+14j) >= 1/13: j in {[j for j,_,_ in safe]}")
    if unsafe:
        print(f"      14*j DROPPING below 1/13: {[(j,mult,round(mm,4)) for j,mult,mm in unsafe]}")
        print(f"         (these j have 13 | (14j*k) at the optimum -> 13|j since 14==1 mod13)")
    else:
        print(f"      none drop below 1/13 for j=1..13.")
print("="*70); print("DONE")
