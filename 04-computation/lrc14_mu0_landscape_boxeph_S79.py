#!/usr/bin/env python3
"""
THE MU_0 LANDSCAPE ON THE TRAPPED-CORE BOUNDARY + THE LIPSCHITZ CONVERTER
(boxeph-2026-07-17-S79)

THM-984 reduces LRC(14) to: mu_0 > 0 on trapped cores.  This referee:

(A) THE LIPSCHITZ CONVERTER (one line, proved):  at the argmax t_0 the margin
    M = min_i ||v_i t_0|| exceeds 1/14 by M - 1/14, slopes are <= v_max, so
        mu_0 >= 2 (M - 1/14) / v_max.
    ANY quantitative M-floor becomes a mu_0-floor, hence (THM-984 bridge) a
    census certificate at an explicit modulus.  Refereed exactly here
    (exact M via kink+crossing candidates; exact mu_0 via sweep).

(B) THE LANDSCAPE: exact (M, mu_0) on controls (GW tight; deep well) and on
    randomized trapped-core-boundary families (covering-proxy, spread > 13,
    max >= 23, chain-dense, small-relation-carrying): does any non-tight
    family show mu_0 = 0?  (No -- mu_0 = 0 iff M = 1/14 exactly, and the
    converter quantifies.)  Converter tightness ratios mapped; extremal
    (minimum-mu_0) families listed for the fee-ledger/witness routes.
"""
import sys, random
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, '04-computation')

FT = Fr(1, 14)

def dist_int(x):
    f = x - int(x) if x >= 0 else x - (int(x) - 1)
    f = f % 1
    return min(f, 1 - f)

def exact_M(fam):
    """max over t of min_i ||v_i t||: attained at a kink (k/(2v)) or a
    crossing (k/(v_i +- v_j)); evaluate exactly."""
    cands = set()
    for v in fam:
        for k in range(2 * v + 1):
            cands.add(Fr(k, 2 * v))
    for i in range(len(fam)):
        for j in range(i + 1, len(fam)):
            for s in (fam[i] + fam[j], abs(fam[i] - fam[j])):
                if s:
                    for k in range(s + 1):
                        cands.add(Fr(k, s))
    best, arg = Fr(0), Fr(0)
    for t in cands:
        m = min(dist_int(v * t) for v in fam)
        if m > best:
            best, arg = m, t
    return best, arg

def exact_mu0(fam):
    """measure{t in [0,1) : all ||v_i t|| > 1/14} by exact sweep."""
    bps = sorted(set([Fr(0), Fr(1)] +
                     [Fr(m, v) + s * FT / v for v in fam
                      for m in range(v + 1) for s in (-1, 1)
                      if 0 <= Fr(m, v) + s * FT / v <= 1]))
    tot = Fr(0)
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        if all(dist_int(v * mid) > FT for v in fam):
            tot += bps[i + 1] - bps[i]
    return tot

def covering_proxy(fam):
    return all(any(v % q == 0 for v in fam) for q in range(2, 15))

def chain_dense(fam):
    s = sorted(fam)
    return any(Fr(s[i + 1], s[i]) < Fr(7, 3) for i in range(len(s) - 1))

def small_relation(fam):
    """exists eps in {-1,0,1}^13, support <= 4, sum eps v = 0 (nontrivial)."""
    from itertools import combinations, product
    for k in (2, 3, 4):
        for idx in combinations(range(len(fam)), k):
            for signs in product((-1, 1), repeat=k - 1):
                if fam[idx[0]] + sum(s * fam[idx[j + 1]]
                                     for j, s in enumerate(signs)) == 0:
                    return True
    return False

if __name__ == "__main__":
    print("THE MU_0 LANDSCAPE + LIPSCHITZ CONVERTER (boxeph S79)")
    print("=" * 76)
    controls = [("GW tight", [1,2,3,4,5,6,7,8,9,10,11,13,24]),
                ("deep well w=14", [1,2,3,4,5,6,7,8,9,10,11,12,14]),
                ("deep well w=25", [1,2,3,4,5,6,7,8,9,10,11,12,25]),
                ("doubling", [3,6,12,24,48,96,192,384,768,1536,23,46,92])]
    rows = []
    for name, fam in controls:
        M, t0 = exact_M(fam)
        mu0 = exact_mu0(fam)
        vmax = max(fam)
        conv = 2 * (M - FT) / vmax if M > FT else Fr(0)
        ok = mu0 >= conv
        ratio = float(mu0 / conv) if conv > 0 else float('inf')
        rows.append((name, fam, M, mu0, conv, ratio))
        print(f"  [{name}] M = {M} ({float(M):.5f}); mu0 = {float(mu0):.6f}; "
              f"converter floor = {float(conv):.6f} (ratio x{ratio:.1f}); "
              f"law holds: {ok}")
        assert ok
    print()
    print("  randomized trapped-core-boundary sample (covering-proxy + spread>13")
    print("  + max>=23 + chain-dense + small-relation):")
    rng = random.Random(79)
    found = 0
    zero_mu = []
    worst = None
    while found < 12:
        base = sorted(rng.sample(range(1, 13), 9))
        big = sorted(rng.sample(range(23, 60), 4))
        fam = sorted(set(base + big))
        if len(fam) != 13:
            continue
        if not (covering_proxy(fam) and chain_dense(fam)
                and max(fam) / min(fam) > 13 and small_relation(fam)):
            continue
        found += 1
        M, t0 = exact_M(fam)
        mu0 = exact_mu0(fam)
        vmax = max(fam)
        conv = 2 * (M - FT) / vmax if M > FT else Fr(0)
        assert mu0 >= conv
        if mu0 == 0:
            zero_mu.append(fam)
        if worst is None or mu0 < worst[0]:
            worst = (mu0, fam, M)
        print(f"    {fam}: M = {float(M):.5f}, mu0 = {float(mu0):.6f}, "
              f"conv floor {float(conv):.6f}"
              + ("  <-- MU0 = 0!" if mu0 == 0 else ""))
    print()
    print(f"  zero-mu0 families in sample: {zero_mu if zero_mu else 'NONE'}")
    print(f"  extremal (min mu0): {worst[1]} with mu0 = {float(worst[0]):.6f}, "
          f"M = {float(worst[2]):.6f}")
    print("=" * 76)
    print("done")
