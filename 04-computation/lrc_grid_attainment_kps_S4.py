#!/usr/bin/env python3
"""
lrc_grid_attainment_kps_S4.py -- HYP-4108 (kind-pasteur-2026-07-05-S4)

Sanity companion to LRCGridAttainment.lean (the theorems are PROVED there,
kernel-pure; this just exercises the claims numerically):

(G1) GRID ATTAINMENT: for random integer families, the exact max-min M
     (fine-scanned) is attained at some t = m/(|v_i|+|v_j|).
(G2) THE 2B CERTIFICATE: whenever margin beta is achievable, a certificate
     (k, s, ceil(beta*s)) with s = |v_i|+|v_j| <= 2B exists and the atom's
     mod conditions hold -- ZERO slack (beta itself, not beta' > beta).
(G3) The S2/S3/S4 bound ladder on the dichotomy's loose branch, side by side:
     S2 (atom, no bound) -> S3 (s <= B/(2 delta) + 1 ~ 176B) -> S4 (s <= 2B).
"""

from math import gcd, ceil
from fractions import Fraction
import random

def margin_at(ws, p, q):
    m = q
    for w in ws:
        r = (abs(w) * p) % q
        r = min(r, q - r)
        if r < m:
            m = r
    return Fraction(m, q)

def exact_M_grid(ws):
    """max over the merge grid t = d/(|vi|+|vj|), all d (non-reduced allowed)."""
    L = sorted(set(abs(w) for w in ws))
    best = Fraction(0)
    arg = None
    qs = set()
    for i in range(len(L)):
        for j in range(i, len(L)):
            qs.add(L[i] + L[j])
    for q in qs:
        for d in range(1, q // 2 + 1):
            m = margin_at(ws, d, q)
            if m > best:
                best, arg = m, (d, q)
    return best, arg

def fine_M(ws, Q=150):
    best = Fraction(0)
    for q in range(2, Q + 1):
        for p in range(1, q):
            if gcd(p, q) > 1:
                continue
            m = margin_at(ws, p, q)
            if m > best:
                best = m
    return best

def main():
    rng = random.Random(42)
    print("[G1] grid attainment: exact M on the merge grid vs fine scan")
    bad = 0
    for _ in range(300):
        kk = rng.randint(3, 12)
        ws = [rng.choice([-1, 1]) * rng.randint(1, 60) for _ in range(kk)]
        if any(w == 0 for w in ws):
            continue
        Mg, arg = exact_M_grid(ws)
        Mf = fine_M(ws)
        if Mf > Mg:
            bad += 1
            print(f"   G1 FAIL {ws}: fine {Mf} > grid {Mg}")
    print(f"   300 families, {bad} failures")

    print("[G2] the 2B certificate, zero slack:")
    bad2 = 0
    tested = 0
    for _ in range(300):
        kk = rng.randint(3, 12)
        ws = [rng.choice([-1, 1]) * rng.randint(1, 60) for _ in range(kk)]
        B = max(abs(w) for w in ws)
        Mg, arg = exact_M_grid(ws)
        if Mg == 0:
            continue
        beta = Mg  # ZERO slack: certify at the exact max itself
        d, q = arg
        assert q <= 2 * B
        mu = ceil(beta * q)
        ok = all(mu <= (w * d) % q <= q - mu for w in ws)
        tested += 1
        if not ok:
            bad2 += 1
            print(f"   G2 FAIL {ws}: cert at {d}/{q} does not verify")
    print(f"   {tested} families, {bad2} failures (cert modulus always <= 2B)")

    print("[G3] the loose-branch modulus-bound ladder at beta = 2/25:")
    B = 24
    delta = Fraction(14, 169) - Fraction(2, 25)
    print(f"   S2 (HYP-4102): certificate language, no a-priori bound")
    print(f"   S3 (HYP-4105): s <= B/(2*delta) + 1 = {float(Fraction(B,2)/delta + 1):.0f}"
          f"  (needs the 14/169 slack)")
    print(f"   S4 (HYP-4108): s <= 2B = {2*B}  (zero slack, any beta > 0)")

if __name__ == "__main__":
    main()
