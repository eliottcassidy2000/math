#!/usr/bin/env python3
"""death-star-2026-07-17-S39 (HYP-7185): referee for LRCMomentCertificates.lean.

(1) pointwise capped certificate 30(C3+C5) <= 27-27c+28C2+33C4 on c in 0..6 (exact);
(2) the summation floor 30*B5 >= 3(q-1) - 3S1 + 2S2 - 3S4 on capped families,
    verified against direct B5 on random (v, q) with bandCount <= 6 everywhere;
(3) the wall witness: masses (450,702,1248,1)/2401 on c=(0,1,3,13) match the
    equidistributed moments and carry ledger -342/2401 (exact rational identities);
(4) gap accounting: 40x -> 1.08x."""
import random
from fractions import Fraction as Fr
from math import comb, gcd

def in_band(vi, q, p):
    return q <= 14 * ((vi * p) % q) <= 13 * q

def referee_pointwise():
    ok = all(30 * (comb(c, 3) + comb(c, 5))
             <= 27 - 27 * c + 28 * comb(c, 2) + 33 * comb(c, 4) for c in range(7))
    tight = [c for c in range(7) if 30 * (comb(c, 3) + comb(c, 5))
             == 27 - 27 * c + 28 * comb(c, 2) + 33 * comb(c, 4)]
    print(f"pointwise cert: {'PASS' if ok else 'FAIL'}; tight at {tight}")

def referee_floor(trials=6000, seed=39):
    rnd = random.Random(seed)
    ok = True
    n = 0
    while n < 200:
        q = rnd.randint(30, 120)
        v = [rnd.randint(1, 10**5) for _ in range(13)]
        counts = [sum(1 for vi in v if not in_band(vi, q, p)) for p in range(1, q)]
        if max(counts) > 6:
            continue
        n += 1
        S = [sum(comb(c, d) for c in counts) for d in range(6)]
        b5 = sum((-1) ** d * S[d] for d in range(6))
        if not 3 * (q - 1) - 3 * S[1] + 2 * S[2] - 3 * S[4] <= 30 * b5:
            ok = False
            print(f"  FAIL floor at q={q}")
        if n > trials:
            break
    print(f"capped floor: {'PASS' if ok else 'FAIL'} ({n} capped families found)")

def referee_wall():
    h = [Fr(450, 2401), Fr(702, 2401), Fr(1248, 2401), Fr(1, 2401)]
    cs = [0, 1, 3, 13]
    mom = lambda f: sum(hh * f(c) for hh, c in zip(h, cs))
    ok = (mom(lambda c: 1) == 1 and mom(lambda c: c) == Fr(13, 7)
          and mom(lambda c: comb(c, 2)) == Fr(78, 49)
          and mom(lambda c: comb(c, 4)) == Fr(715, 2401))
    led = mom(lambda c: sum((-1) ** d * comb(c, d) for d in range(6)))
    print(f"wall witness: moments {'MATCH' if ok else 'FAIL'}; "
          f"ledger = {led} = {float(led):.4f} (want -342/2401 = {float(Fr(-342,2401)):.4f})")

def gap_accounting():
    E = Fr(2052, 16807)
    trivial = 40  # THM-944's per-triple ceiling, ~40(q-1)
    capped_short = Fr(675, 2401) / 30  # 22.5/2401
    print(f"gap: trivial ceiling ~{trivial}(q-1); capped-certificate shortfall "
          f"{float(capped_short):.5f}(q-1) vs equilibrium {float(E):.5f}(q-1): "
          f"ratio {(float(capped_short)+float(E))/float(E):.3f}x (was ~{trivial/float(E):.0f}x)")

if __name__ == "__main__":
    referee_pointwise()
    referee_floor()
    referee_wall()
    gap_accounting()
