"""
lrc14_tail_true_min_opus_S155.py   (opus-2026-07-08-S155)

THE TRUE k=11 prim-diam>=25 tail D3-MINIMUM -- correcting klein-S186/LEM-009.

FINDING (this session): the claimed tail minimizer (block+outlier {0..9,25}, D3=0.4587) is NOT the
minimum.  The DILATED-BLOCK-with-interior-point (0,3,6,8,9,12,15,18,21,24,27) = 3*block_10 u {8}
has D3 = 0.452986 < 0.4587 (prim-diam 27, primitive), the tail analog of the EXHAUSTIVE minimizer
2*block_10 u {9} = (0,2,4,6,8,9,10,12,14,16,18), D3 = 0.4356 at prim-diam 18.  D3 is
dilation-invariant while klein's fixed-window "cluster size" is NOT -- so the cluster-monotonicity
proof device (D3 >= D3_{c(E)}, D3_c decreasing) FAILS: these shapes have small window-cluster but
low D3.  The GLOBAL extremal conclusion must be re-derived.

This script finds the true tail infimum (thorough float search seeded by the low-D3 structures +
broad random + perturbed exhaustive-minimizers), exact-verifies the lowest, and checks it clears
bar (the closure question).  METHOD: fast numpy float D3 (matches exact to ~5e-7), exact-verify
lowest candidates.
"""
import sys, random, itertools
import numpy as np
from fractions import Fraction as Fr
import math

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_d3_exact_verify_klein_S184 import D3 as D3_exact, BAR

K = 11
TH = 1.0 / 7
Mval = 6.0 / 7
NG = 9000
_X = (np.arange(NG) + 0.5) / NG


def D3f(E):
    A = np.array(sorted(E), float)
    ph = (A[:, None] * _X[None, :]) % 1.0
    ph.sort(axis=0)
    g = np.empty_like(ph)
    g[:-1] = ph[1:] - ph[:-1]; g[-1] = 1 - ph[-1] + ph[0]
    W = np.clip(g - TH, 0, None).sum(0)
    m1, m2, m3 = W.mean(), (W * W).mean(), (W * W * W).mean()
    den = m2 - m3 / Mval
    return m1 / Mval + (m1 - m2 / Mval) ** 2 / den if den > 0 else m1 / Mval


def prim(E):
    E = sorted(E); g = 0
    for i in range(1, len(E)):
        g = math.gcd(g, E[i] - E[0])
    return g == 1, (E[-1] - E[0]) // max(g, 1)


def keep(E):
    E = tuple(sorted(set(E)))
    if len(E) != K:
        return None
    p, d = prim(E)
    if not p or d < 25:
        return None
    return E


def main():
    print("=" * 96)
    print("TRUE k=11 prim-diam>=25 tail D3-minimum (correcting klein tail-min 0.4587)")
    print(f"  bar = {float(BAR):.6f}")
    print("=" * 96)
    pool = {}

    def add(E):
        E = keep(E)
        if E and E not in pool:
            pool[E] = D3f(E)

    r = random.Random(11)
    # 1) dilated block_10 + interior/exterior point, all d, all positions
    for d in range(2, 12):
        blk = [d * i for i in range(10)]
        for p in range(1, d * 9 + 80):
            add(blk + [p])
    # 2) dilated block_c (c=8,9) + points
    for d in range(2, 10):
        for c in (8, 9):
            blk = [d * i for i in range(c)]
            for _ in range(1500):
                add(blk + r.sample(range(1, d * (c - 1) + 90), K - c))
    # 3) perturbations of the exhaustive low-D3 minimizers, dilated
    seeds = [(0, 2, 4, 6, 8, 9, 10, 12, 14, 16, 18), (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 17),
             (0, 2, 4, 6, 8, 9, 10, 12, 14, 16, 17), (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]
    for s in seeds:
        for d in range(2, 8):
            base = [d * x for x in s]
            for _ in range(2000):
                E = base[:]
                j = r.randrange(len(E)); E[j] += r.choice([-2, -1, 1, 2])
                add(E)
    # 4) broad random
    for _ in range(30000):
        D = r.randint(25, 70)
        add([0] + r.sample(range(1, D), K - 2) + [D])
    # 5) two-scale (block + dilated tail)
    for _ in range(8000):
        d = r.choice([2, 3, 4])
        a = r.randint(5, 8)
        E = [d * i for i in range(a)] + r.sample(range(d * a, d * a + 90), K - a)
        add(E)

    print(f"  pooled {len(pool)} distinct primitive tail shapes")
    ranked = sorted(pool.items(), key=lambda kv: kv[1])
    print("\n  lowest 15 by float-D3 (then EXACT-verified):")
    print(f"    {'shape':46s} {'floatD3':>9s} {'exactD3':>9s} {'margin':>9s}")
    truemin = None
    for E, vf in ranked[:15]:
        ve = float(D3_exact(E))
        if truemin is None or ve < truemin[0]:
            truemin = (ve, E)
        _, pd = prim(E)
        print(f"    {str(E):46s} {vf:9.6f} {ve:9.6f} {ve-float(BAR):+9.6f}  pd={pd}")
    print("\n" + "=" * 96)
    print(f"TRUE tail min (exact, over search) = {truemin[0]:.6f} at {truemin[1]}")
    print(f"  vs klein's claimed 0.458714 (block+outlier) -- CORRECTED DOWN")
    print(f"  vs exhaustive compact min 0.435613 (prim-diam 18)")
    print(f"  margin over bar = +{truemin[0]-float(BAR):.6f}  => closure {'HOLDS' if truemin[0]>=float(BAR) else 'FAILS'}")
    print("=" * 96)


if __name__ == "__main__":
    main()
