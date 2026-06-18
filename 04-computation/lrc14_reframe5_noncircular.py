#!/usr/bin/env python3
"""
REFRAME 5 part (2): THE non-circularity test, done carefully.

Established structure (mac-mini-S5):
  N=14 level config S = ({1..13}\{q}) u {w}, w=lcm(q,14).
  When w binds, the optimum is a SUM-crossing of (flank, w):
     tau* = j_tau/D,  D = flank + w,  gap = M(S) = j/D  (j unreduced crossing index),
  and M(S)>=1/14  <=>  D <= 14*j   (tautology: this IS the goal restated).

THE QUESTION (non-circular leverage): can we LOWER-BOUND j, or UPPER-BOUND the dip
  dip = M(A) - j/D
WITHOUT assuming M(S)>=1/14, using only arithmetic of (q,14,flank,w)?

FROZEN-CORE idea (C5): the core A=({1..13}\{q}) has gap M(A) at core-time(s) tau_A,
chosen WITHOUT reference to w. If at SOME core-optimal tau_A the parked runner w is
already at distance >= 1/14, then tau_A witnesses M(S) >= 1/14 directly -- the dip is a
red herring, and the bound is NON-CIRCULAR (it uses M(A) >= 1/14, the core's gap by
induction, plus a single arithmetic distance check ||w*tau_A|| >= 1/14).
"""
from fractions import Fraction as F
from math import gcd
import itertools


def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r


def g(S, t):
    return min(nrm(v * t) for v in S)


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


def Mall(S):
    """return (M, list_of_all_argmaxes)."""
    b = F(0); ats = []
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; ats = [t]
        elif v == b:
            ats.append(t)
    return b, ats


def lcm(a, b):
    return a * b // gcd(a, b)


print("=" * 78)
print("REFRAME 5 (2): FROZEN-CORE non-circular test")
print("=" * 78)
print("""
KEY IDEA (C5): the core A=({1..13}\\{q}) has gap M(A) attained at core-time(s) tau_A,
chosen WITHOUT reference to w. If at SOME core-optimal tau_A the parked runner w is
already >= 1/14 away, then tau_A witnesses M(S) >= 1/14 directly -- the dip is a red
herring, and the bound is NON-CIRCULAR (uses only M(A) >= 1/14 + one distance check).
""")

print(f"{'q':>2} {'M(A)':>6} {'#tauA':>6} {'best min(M(A),||w*tauA||)':>26} "
      f"{'clears?':>8} {'M(S)':>8} {'>=1/14':>7}")
allclear = True
for q in range(2, 14):
    A = [v for v in range(1, 14) if v != q]
    MA, tausA = Mall(A)
    w = lcm(q, 14)
    S = A + [w]
    MS, _ = Mall(S)
    best = F(0)
    for t in tausA:
        full = min(MA, nrm(w * t))
        if full > best:
            best = full
    clears = best >= F(1, 14)
    if not clears:
        allclear = False
    print(f"{q:>2} {str(MA):>6} {len(tausA):>6} {str(best):>26} "
          f"{str(clears):>8} {str(MS):>8} {str(MS>=F(1,14)):>7}")

print(f"\n  => every N=14-level single-drop core time already clears w (>=1/14)? {allclear}")

print("\n" + "-" * 78)
print("STRESS: frozen-core 1/14-witness across resonant multiples and multi-drop cores")
print("-" * 78)

fails = []; tested = 0
for q in range(2, 14):
    A = [v for v in range(1, 14) if v != q]
    MA, tausA = Mall(A)
    if MA < F(1, 14):
        continue
    L = lcm(q, 14)
    for k in range(1, 13):
        w = L * k
        if w in A:
            continue
        tested += 1
        best = max(min(MA, nrm(w * t)) for t in tausA)
        if best < F(1, 14):
            S = A + [w]
            MS, _ = Mall(S)
            fails.append((q, k, w, str(best), str(MS), MS >= F(1, 14)))

print(f"  single-drop resonant tests: {tested}, frozen-core failures: {len(fails)}")
for f in fails[:25]:
    print("   FAIL", f)

tested2 = 0; fails2 = []
for q1, q2 in itertools.combinations(range(2, 14), 2):
    A = [v for v in range(1, 14) if v not in (q1, q2)]
    MA, tausA = Mall(A)
    if MA < F(1, 14):
        continue
    for k in range(1, 4):
        L = lcm(lcm(q1, q2), 14)
        w = L * k
        if w in A:
            continue
        tested2 += 1
        best = max(min(MA, nrm(w * t)) for t in tausA)
        if best < F(1, 14):
            S = A + [w]
            MS, _ = Mall(S)
            fails2.append((q1, q2, k, w, str(best), str(MS), MS >= F(1, 14)))

print(f"  two-drop resonant tests: {tested2}, frozen-core failures: {len(fails2)}")
for f in fails2[:30]:
    print("   FAIL", f)

print("\nDONE.")
