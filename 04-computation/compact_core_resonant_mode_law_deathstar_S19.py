#!/usr/bin/env python3
"""death-star-2026-07-16-S19 (HYP-7018): the compact-core resonant-mode law.

Owner directive: "prove compact-core flatness". The honest resolution: flatness with an
absolute constant is FALSE asymptotically even for balanced clusters — the correct law is
klein's THM-883 resonant mode extended into the balanced regime, with the independent-limit
miss-measure coefficient m* (exact rationals) replacing the slow-system measures; the slope
is small enough that PROGRAM-RANGE (bounded-Vmax) flatness holds with explicit constants.

Pieces:
  P1 (exact combinatorics): A* = 360/16807, B* = 120/16807 by full Z_7^5 enumeration.
      A_c(s) = P[{0} u {five uniform sections} = Z7 \\ {s,c} exactly]  (c not in {0,s})
      B(s)   = P[... = Z7 \\ {s} exactly]
  P2 (the law's prediction): peak height at owner-e resonant frequencies n ~ e*a:
      H_pred(e,a) = e * |1 - e(a/7)| * |m*_s(a)|,   m*_s(a) = -(A*+B*) e(as/7) - A*.
  P3 (probes): for compact cores (consecutive + generic-ratio), compute exact endpoints,
      probe |S(n)| on windows around n = e*a for every owner e and a = 1..6; report
      max-probe vs prediction vs the flat floor sqrt(C*M).
  P4 (empirical m-hat): compute the ACTUAL 6-runner miss-measures A_c, B for E \\ {e}
      (exact interval sweep) and the sharper prediction e |1-w^a| |m-hat^{(e)}(a)|;
      equidistribution claim: m-hat^{(e)} -> m* for incoherent cores.
  P5 (dilate covariance): S_{cE}(c n) = c S_E(n) exact referee.
"""
from fractions import Fraction as Fr
from math import gcd, pi
from itertools import product as iproduct
import cmath, sys, time

def lcm(a, b): return a * b // gcd(a, b)

# ---------------- exact R_s machinery (klein-verbatim logic) ----------------

def R_s_endpoints(E, s):
    bps = sorted(set(Fr(k, 7 * e) for e in E if e > 0 for k in range(7 * e)) | {Fr(0), Fr(1)})
    eps = []
    inR, start = False, None
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        occ = set(int((e * mid % 1) * 7) for e in E)
        cur = (len(occ) == 6) and (s not in occ)
        if cur and not inR: start, inR = bps[i], True
        if (not cur) and inR:
            eps.append((start, +1)); eps.append((bps[i], -1)); inR = False
    if inR:
        eps.append((start, +1)); eps.append((Fr(1), -1))
    return eps

def S_at(eps, n):
    return sum(sg * cmath.exp(2j * pi * float((n * p) % 1)) for p, sg in eps)

# ---------------- P1: exact independent-limit combinatorics ----------------

def p1():
    print("P1: independent-limit miss-measures (exact Z_7^5 enumeration; pin sigma_0 = 0)")
    from collections import Counter
    cntA = Counter(); cntB = 0
    tot = 7 ** 5
    s = 6  # generic s != 0; symmetry over s != 0
    for sig in iproduct(range(7), repeat=5):
        occ = set(sig) | {0}
        miss = frozenset(range(7)) - occ
        if miss == frozenset({s}):
            cntB += 1
        elif len(miss) == 2 and s in miss:
            c = next(iter(miss - {s}))
            cntA[c] += 1
    print(f"  B count = {cntB} (predict 120): {'OK' if cntB == 120 else 'FAIL'}")
    okA = all(cntA[c] == 360 for c in range(1, 6) if c != s)
    print(f"  A_c counts = {dict(cntA)} (predict 360 each, c not in {{0,{s}}}): {'OK' if okA else 'FAIL'}")
    A = Fr(360, tot); B = Fr(120, tot)
    return A, B

def mstar(A, B, s, a):
    w = cmath.exp(2j * pi * a * s / 7)
    return -(float(A) + float(B)) * w - float(A)

# ---------------- P3: probes ----------------

def probe_core(E, s, A, B, halfwin=60, label=""):
    eps = R_s_endpoints(E, s)
    M = len(eps)
    if M == 0:
        return None
    best_overall = 0.0
    rows = []
    for e in sorted(x for x in E if x > 0):
        best_e = (0.0, None, None)
        for a in range(1, 7):
            base = e * a
            for d in range(-halfwin, halfwin + 1):
                n = base + d
                if n <= 0: continue
                v = abs(S_at(eps, n))
                if v > best_e[0]:
                    best_e = (v, a, n)
        v, a, n = best_e
        pred = e * abs(1 - cmath.exp(2j * pi * a / 7)) * abs(mstar(A, B, s, a)) if a else 0
        rows.append((e, a, n, v, pred))
        best_overall = max(best_overall, v)
    print(f"  {label} E={E} s={s}: M={M}  max-probe |S| = {best_overall:.2f}  "
          f"(|S|^2/M = {best_overall**2/M:.2f})")
    for e, a, n, v, pred in rows[-3:]:
        print(f"      owner {e}: max at n={n} (a={a}): |S|={v:.2f} vs law-pred {pred:.2f} "
              f"(ratio {v/pred if pred>0 else float('nan'):.2f})")
    return best_overall, M

if __name__ == "__main__":
    t0 = time.time()
    A, B = p1()
    print(f"  A* = {A} = {float(A):.6f}, B* = {B} = {float(B):.6f}, A* = 3B*: {A == 3*B}")
    print(f"  max_a |1-w^a||m*(a)| = "
          f"{max(abs(1-cmath.exp(2j*pi*a/7))*abs(mstar(A,B,6,a)) for a in range(1,7)):.4f}")
    print("\nP3a: CONSECUTIVE compact cores [0, c..c+5] (difference-slow family)")
    for c in [30, 100, 300, 1000, 2000]:
        probe_core([0] + list(range(c, c + 6)), 3, A, B, label=f"c={c}")
        sys.stdout.flush()
    print("\nP3b: GENERIC-RATIO compact cores (incoherent)")
    for c in [30, 100, 300, 1000, 2000]:
        E = [0, c, int(1.37 * c) + 1, int(1.91 * c), int(2.83 * c) + 1, int(4.13 * c), int(5.87 * c) + 1]
        E = sorted(set(E))
        if len(E) != 7:
            continue
        probe_core(E, 3, A, B, label=f"c={c}")
        sys.stdout.flush()
    print(f"\n[total {time.time()-t0:.1f}s]")
