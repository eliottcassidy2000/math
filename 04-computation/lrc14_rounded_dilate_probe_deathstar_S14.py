#!/usr/bin/env python3
"""
lrc14_rounded_dilate_probe_deathstar_S14.py — death-star-2026-07-12-S14

MISTAKE-101 discipline applied to my own THM-720 addendum: is the adversarial large-diameter
DC floor really CONSTANT 1/13, or does the j>=7 stratum (rounded REAL-ratio dilates of the
tight AP, v_i = round(lambda*i)) dip below 1/13 toward 1/14?

Probe: lambda = Lnum/Lden rational non-integer in [500, 2000]; v_i = round(lambda*i), i=1..13;
keep primitive DC hits; compute exact M (THM-668 pair-sum enumeration). Track anything with
M < 1/13. Also compute the off-lattice count j at the natural scale round(lambda) and the
2D-profile reach2 prediction via the pure sub-family.

Also: adversarial OFFSET probe at fixed L: hill-climb b in [-6,6]^13 (j>=7 patterns) minimizing
exact M of {L*i + b_i}, to see whether designed offsets can push M below 1/13 (u-obstruction).
"""
from fractions import Fraction
from math import gcd
from functools import reduce
import random, sys, time

def exact_M(v):
    vv = tuple(sorted(v))
    n = len(vv)
    rulers = sorted(set(vv[i] + vv[j] for i in range(n) for j in range(i, n)))
    bn, bd = 0, 1
    bq, bp = None, None
    for q in rulers:
        thr = (bn * q) // bd
        half = q >> 1
        for p in range(1, half + 1):
            m = q
            for x in vv:
                r = (x * p) % q
                if r > half:
                    r = q - r
                if r < m:
                    m = r
                    if m <= thr:
                        break
            if m > thr and m * bd > bn * q:
                bn, bd, bq, bp = m, q, q, p
                g = gcd(bn, bd)
                bn //= g; bd //= g
                thr = (bn * q) // bd
    return Fraction(bn, bd), bq, bp

def is_dc(v):
    return all(any(x % d == 0 for x in v) for d in range(2, 15))

def offlattice_j(v, L):
    return sum(1 for x in v if x - L * ((x + L // 2) // L) != 0)

ONE13 = Fraction(1, 13)
ONE14 = Fraction(1, 14)

print("=" * 78)
print("PROBE 1: rounded real-ratio dilates v_i = round(lambda*i), DC + primitive hits")
print("=" * 78)
rng = random.Random(14)
hits = 0
below_13 = []
t0 = time.time()
tried = 0
while hits < 60 and time.time() - t0 < 360:
    tried += 1
    Lden = rng.randint(2, 40)
    Lnum = rng.randint(500 * Lden, 2000 * Lden)
    if Lnum % Lden == 0:
        continue
    lam = Fraction(Lnum, Lden)
    v = sorted(set(int(round(float(lam * i))) for i in range(1, 14)))
    if len(v) != 13 or reduce(gcd, v) != 1 or not is_dc(v):
        continue
    hits += 1
    M, q, p = exact_M(v)
    L0 = int(round(float(lam)))
    j0 = offlattice_j(v, L0)
    flag = ""
    if M < ONE13:
        below_13.append((float(lam), v, M))
        flag = "  <<< BELOW 1/13"
    if M < ONE14:
        flag = "  <<<<<< BELOW 1/14 (!!)"
    print(f"  lam={float(lam):9.3f}  j@L0={j0:2d}  M = {str(M):>10} = {float(M):.6f}{flag}")
    sys.stdout.flush()
print(f"\n  DC hit rate: {hits}/{tried} tried; below 1/13: {len(below_13)}; "
      f"min M seen: {min((float(m) for _,_,m in below_13), default=float('nan')) if below_13 else 'none below'}")
if below_13:
    lam, v, M = min(below_13, key=lambda x: x[2])
    print(f"  WORST: lam={lam:.4f} v={v} M={M}={float(M):.6f}")

print()
print("=" * 78)
print("PROBE 2: adversarial offsets at L=1560 — hill-climb b (|b_i|<=6) minimizing exact M")
print("=" * 78)
L = 1560
def make(b):
    return sorted(L * i + b[i - 1] for i in range(1, 14))
def eval_b(b):
    v = make(b)
    if len(set(v)) != 13 or reduce(gcd, v) != 1:
        return None
    return exact_M(v)[0]

best_overall = None
for restart in range(4):
    rngh = random.Random(100 + restart)
    b = [rngh.randint(-6, 6) for _ in range(13)]
    if all(x == 0 for x in b):
        b[12] = 1
    cur = eval_b(b)
    while cur is None:
        b = [rngh.randint(-6, 6) for _ in range(13)]
        cur = eval_b(b)
    improved = True
    evals = 0
    while improved and evals < 260:
        improved = False
        for i in rngh.sample(range(13), 13):
            for d in (-2, -1, 1, 2):
                nb = list(b); nb[i] += d
                if abs(nb[i]) > 6:
                    continue
                m = eval_b(nb)
                evals += 1
                if m is not None and m < cur:
                    b, cur = nb, m
                    improved = True
                    break
            if improved:
                break
    j = sum(1 for x in b if x != 0)
    print(f"  restart {restart}: min M = {cur} = {float(cur):.6f}  at b={b}  (j={j})"
          + ("  <<< BELOW 1/13" if cur < ONE13 else ""))
    sys.stdout.flush()
    if best_overall is None or cur < best_overall[0]:
        best_overall = (cur, b)

M0, b0 = best_overall
print(f"\n  OVERALL adversarial min at L={L}: M = {M0} = {float(M0):.6f}, b = {b0}")
print(f"  vs 1/13 = {float(ONE13):.6f}: {'BELOW — the compressed floor is NOT 1/13 globally' if M0 < ONE13 else 'still >= 1/13 (u-obstruction not found by this climb)'}")
print(f"  vs 1/14 = {float(ONE14):.6f}: {'BELOW (!!)' if M0 < ONE14 else 'above'}")
print("\nDONE.")
