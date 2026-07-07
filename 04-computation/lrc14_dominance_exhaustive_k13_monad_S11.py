#!/usr/bin/env python3
r"""
lrc14_dominance_exhaustive_k13_monad_S11.py   (monad-explorer-2026-07-07-S11, HYP-5157, part 9)

(A) EXHAUSTIVE layer-cake dominance at k=8: M_m(E) <= M_m(AP_8) for all m = 2..7,
    over ALL translation-normalized 8-subsets of {0..N}, N <= 22 (pure counting, no
    integrals; dilation-reduction applied).  This is the finite core of the Sigma_3
    AP-extremality lemma (the weights 1/q' are decreasing, so dominance of cumulative
    counts + Abel summation => Sigma_3(E)'s theta^2-part <= AP's).

(B) k=13 CONSTANTS for the same program on the k=13 leg (the tail-diameter residual
    diam >= 76): universal part = 1 - 13/7 + C(13,2)/49 = 36/49; exact E[V], E[V^2],
    PZ at AP_13; B3 sanity.  The k=13 bar is m_P = 0.0565 (MISTAKE-123: k=13 unaffected).

(C) The dichotomy numbers: how large can Sigma_3(E) be for a diam >= 76 13-shape that
    is NOT a dilated AP?  (probe: dilated-AP-with-defect, two-scale, random spread)
"""
import sys, random
from fractions import Fraction as F
from itertools import combinations
from math import gcd

exec(open('/home/bigo/math/04-computation/lrc14_cubic_moment_gate_monad_S11.py')
     .read().split("if __name__")[0])
_src = open('/home/bigo/math/04-computation/lrc14_window_correlation_calculus_monad_S11.py').read()
_body = _src[:_src.rfind('if __name__')]
_body = _body[_body.index('THETA = F(1, 7)'):]
exec(_body)
THETA = F(1, 7)

def M_counts_vec(E, mmax=7):
    cnt = [0] * (mmax + 1)
    Es = sorted(E)
    for a, b, c in combinations(Es, 3):
        p, q = b - a, c - a
        g = gcd(p, q)
        qr = q // g
        if qr <= mmax:
            cnt[qr] += 1
    # cumulative
    for m in range(3, mmax + 1):
        cnt[m] += cnt[m - 1]
    return cnt[2:]

if __name__ == "__main__":
    print("=" * 100)
    print("PART 9a -- EXHAUSTIVE M_m DOMINANCE, k=8, all 0-anchored 8-subsets of {0..N}")
    print("=" * 100)
    ap = M_counts_vec(range(8))
    print(f"  AP_8 M_m (m=2..7): {ap}")
    total = viol = 0
    worst_by_m = [0] * 6
    argmax_by_m = [None] * 6
    for N in range(7, 23):
        # subsets containing 0 and N (diameter exactly N), translation-normalized
        for mid in combinations(range(1, N), 6):
            E = (0,) + mid + (N,)
            # dilation reduce
            g = 0
            for e in E[1:]:
                g = gcd(g, e)
            if g > 1:
                continue  # counted at smaller N
            c = M_counts_vec(E)
            total += 1
            bad = False
            for i in range(6):
                if c[i] > ap[i]:
                    bad = True
                if c[i] > worst_by_m[i] and tuple(E) != tuple(range(8)):
                    worst_by_m[i] = c[i]
                    argmax_by_m[i] = E
            if bad:
                viol += 1
                if viol <= 5:
                    print(f"  VIOLATION: {E}  M = {c}  vs AP {ap}")
        sys.stdout.flush()
    print(f"  exhaustive over diam <= 22 primitives: {total} shapes, {viol} violations")
    print(f"  runner-up counts per m: {worst_by_m}")
    for i in range(6):
        print(f"    m={i+2}: max non-AP M_m = {worst_by_m[i]} (AP: {ap[i]}) at {argmax_by_m[i]}")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 9b -- k=13 CONSTANTS (leg bar m_P = 14249/252252 = 0.056487)")
    print("=" * 100)
    k = 13
    univ13 = 1 - k * THETA + comb(k, 2) * THETA ** 2
    print(f"  universal part 1 - 13/7 + 78/49 = {univ13} = {float(univ13):.6f}")
    for name, E in [("AP_13", list(range(13))),
                    ("GW", [0,1,2,3,4,5,6,7,8,9,10,12,23]),
                    ("record-diam43ish", [0,2,4,5,6,7,8,9,10,11,12,14,43]),
                    ("dilAP7+defect (diam 85)", [0,7,14,21,28,35,42,49,56,63,70,77,85]),
                    ("spread-random diam 100", [0,9,17,26,38,47,55,64,71,80,88,95,100])]:
        aV, MV, m3V, vmaxV = excess_moments(E, [THETA])
        pz = aV[0] ** 2 / MV[0][0]
        s3 = sigma_m(E, 3)
        b3 = univ13 - s3
        at3 = 1 - atom_zero_bound_3mom(aV[0], MV[0][0], m3V, vmaxV[0])
        print(f"  {name:24s} E[V]={float(aV[0]):.5f} E[V^2]={float(MV[0][0]):.5f} "
              f"PZ={float(pz):.4f} AT3={float(at3):.4f} | Sigma_3={float(s3):.4f} "
              f"B3={float(b3):+.4f}  (bar 0.0565)")
        sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART 9c -- SPREAD-SIDE HEADROOM at k=13: PZ vs bar on spread shapes (the diam>=76 residual)")
    print("=" * 100)
    rng = random.Random(9137)
    worst_pz = (None, None)
    for trial in range(12):
        E = sorted(rng.sample(range(1, 150), 12))
        E = [0] + [e for e in E]
        E = sorted(set(E))[:13]
        if len(E) < 13 or (max(E) - min(E)) < 76:
            continue
        aV, MV, _, _ = excess_moments(E, [THETA])
        pz = float(aV[0] ** 2 / MV[0][0])
        if worst_pz[0] is None or pz < worst_pz[0]:
            worst_pz = (pz, E)
        print(f"  diam {max(E)-min(E):3d}: PZ = {pz:.4f}  {E}")
        sys.stdout.flush()
    print(f"  worst random spread PZ = {worst_pz[0]:.4f} (bar 0.0565) at {worst_pz[1]}")
