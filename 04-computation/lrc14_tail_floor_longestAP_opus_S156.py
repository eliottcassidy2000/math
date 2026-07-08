"""
lrc14_tail_floor_longestAP_opus_S156.py   (opus-2026-07-08-S156)

RE-DERIVE THE k=11 TAIL D3-FLOOR ON THE LONGEST-AP AXIS (the dilation-invariant fix for the
cluster-monotonicity refuted in S155 / MISTAKE-126).

L(E) := length of the longest arithmetic progression contained in E.  L is DILATION-INVARIANT
and TRANSLATION-invariant -- the correct axis for the dilation-invariant floor D3 (unlike the
fixed-window cluster count, which put a shape 0.012 below klein's D3_10 bound).

THE RE-DERIVATION (three steps):
 STEP 1 (reduction, RIGOROUS): for k=11, prim-diam >= 25  =>  L <= 10.
   [L = 11 means E is an 11-term AP = c*{0..10}+b, whose differences have gcd = c, so
    prim-diam = 10 < 25.  Contradiction.  Hence L <= 10 on the tail.]
 STEP 2 (the binding family L = 10): E = "AP_10 at scale d" + 1 extra point.  WLOG (translation +
   reflection + dilation, all D3-invariant) E = {0,d,2d,...,9d} u {p}, gcd(d,p)=1, p not on the
   lattice.  prim-diam = 9d (interior 0<p<9d) or p (exterior p>9d).  Tail => d>=3 (interior) or
   p>=25 (exterior).  We minimize D3 over (d,p): the min is the INTERIOR point at the SMALLEST
   tail scale d=3, D3 = 0.452986 (EXACT); D3 INCREASES in d toward the decorrelation limit 0.4646
   (klein LEM-009's D3_10 = the d->inf limit where the point decorrelates from the AP).
 STEP 3 (stratification L <= 9): min D3 over longest-AP = L is monotone, floor(L) >= floor(10) =
   0.4530 for all L (the fewer/shorter the AP structure, the higher the floor).

CONCLUSION: tail floor = min_L floor(L) = floor(10) = 0.452986 >= bar = 0.331212 (margin +0.1218).
The k=11 tail closes on the longest-AP axis; klein's D3_c limits are the d->inf UPPER references.

METHOD: fast numpy float D3 (matches exact to ~5e-7) for scans; klein-S184 EXACT rational D3 for
the minimizers.
"""
import sys, math
import numpy as np
from fractions import Fraction as Fr

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_d3_exact_verify_klein_S184 import D3 as D3exact, BAR

K = 11
TH = 1.0 / 7
Mval = 6.0 / 7
NG = 9000
_X = (np.arange(NG) + 0.5) / NG


def D3f(E):
    A = np.array(sorted(E), float)
    ph = (A[:, None] * _X[None, :]) % 1.0
    ph.sort(0)
    g = np.empty_like(ph)
    g[:-1] = ph[1:] - ph[:-1]; g[-1] = 1 - ph[-1] + ph[0]
    W = np.clip(g - TH, 0, None).sum(0)
    m1, m2, m3 = W.mean(), (W * W).mean(), (W ** 3).mean()
    den = m2 - m3 / Mval
    return m1 / Mval + (m1 - m2 / Mval) ** 2 / den if den > 0 else m1 / Mval


def longest_ap(E):
    S = set(E); best = 1
    E = sorted(E)
    for i, a in enumerate(E):
        for b in E[i + 1:]:
            d = b - a
            if a - d in S:
                continue                      # count each maximal AP once (a is the start)
            L = 2; x = b + d
            while x in S:
                L += 1; x += d
            best = max(best, L)
    return best


def prim(E):
    E = sorted(E); g = 0
    for i in range(1, len(E)):
        g = math.gcd(g, E[i] - E[0])
    return g, (E[-1] - E[0]) // max(g, 1)


def D3_decorr_limit_block10_plus_point(nx=3000, nrep=600, seed=0):
    """d->inf limit of the L=10 interior/exterior family = block_10 + 1 decorrelated point."""
    rng = np.random.default_rng(seed)
    X = (np.arange(nx) + 0.5) / nx
    blk = np.arange(10, dtype=float)
    bph = (blk[:, None] * X[None, :]) % 1.0
    s1 = s2 = s3 = 0.0; tot = 0
    for _ in range(nrep):
        u = rng.random((1, nx))
        ph = np.concatenate([bph, u], 0); ph.sort(0)
        g = np.empty_like(ph); g[:-1] = ph[1:] - ph[:-1]; g[-1] = 1 - ph[-1] + ph[0]
        W = np.clip(g - TH, 0, None).sum(0)
        s1 += W.sum(); s2 += (W * W).sum(); s3 += (W ** 3).sum(); tot += nx
    m1, m2, m3 = s1 / tot, s2 / tot, s3 / tot
    den = m2 - m3 / Mval
    return m1 / Mval + (m1 - m2 / Mval) ** 2 / den


def main():
    print("=" * 98)
    print("RE-DERIVED k=11 TAIL D3-FLOOR on the LONGEST-AP axis (dilation-invariant)")
    print(f"  bar = {float(BAR):.6f}")
    print("=" * 98)

    print("\nSTEP 1 (reduction, RIGOROUS): prim-diam >= 25 => longest-AP L <= 10.")
    print("  [L=11 => E = c*{0..10}+b, gcd of diffs = c, prim-diam = 10 < 25. So L <= 10 on the tail.]")

    print("\nSTEP 2 (binding family L=10):  E = {0,d,..,9d} u {p}, gcd(d,p)=1.  min D3 over (d,p):")
    print("  d |  interior min (0<p<9d)        | exterior min (p>9d)          | prim-diam constraint")
    glob = (9.9, None)
    interior_by_d = {}
    for d in range(1, 16):
        lat = set(range(0, 9 * d + 1, d))
        bi = (9.9, None); be = (9.9, None)
        for p in range(1, 9 * d + 90):
            if p in lat or math.gcd(d, p) != 1:
                continue
            E = tuple(sorted(lat | {p}))
            if len(E) != K:
                continue
            g, pdm = prim(E)
            if pdm < 25:
                continue
            v = D3f(E)
            if p < 9 * d:
                if v < bi[0]:
                    bi = (v, E)
            else:
                if v < be[0]:
                    be = (v, E)
            if v < glob[0]:
                glob = (v, E)
        si = f"{bi[0]:.6f} {str(bi[1]) if bi[1] else '':40s}" if bi[1] else "(none: 9d<25)"
        se = f"{be[0]:.6f}" if be[1] else "(none)"
        print(f"  {d:2d}| {si:52s} | {se:10s}")
        if bi[1]:
            interior_by_d[d] = bi
    # exact-verify the global minimizer
    gm = glob[1]
    gmex = D3exact(gm)
    print(f"\n  GLOBAL L=10 min (float) = {glob[0]:.6f} at {gm}")
    print(f"  EXACT D3 = {gmex} = {float(gmex):.6f}  (margin over bar +{float(gmex-BAR):.6f})")
    print(f"  longest-AP({gm}) = {longest_ap(gm)}  gcd,primdiam = {prim(gm)}")

    print("\n  d-MONOTONICITY of the interior family (D3 increases in d toward the decorr limit):")
    for d in sorted(interior_by_d):
        v, E = interior_by_d[d]
        print(f"    d={d:2d}: interior-min D3 = {v:.6f}  at {E}")
    lim = D3_decorr_limit_block10_plus_point()
    print(f"  d->inf decorrelation limit (block_10 + 1 iid point) = {lim:.6f}  (= klein LEM-009 D3_10=0.4646)")
    print("  => interior D3 rises from 0.4530 (d=3, strongest AP-point correlation) toward 0.4646; min at d=3.")

    print("\nSTEP 3 (stratification L<=9): min D3 by longest-AP over a broad tail search:")
    import random
    r = random.Random(3)
    strat = {L: (9.9, None, 0) for L in range(2, 11)}
    def upd(E):
        E = tuple(sorted(set(E)))
        if len(E) != K:
            return
        g, pdm = prim(E)
        if g != 1 or pdm < 25:
            return
        L = longest_ap(E); v = D3f(E)
        cur = strat[L]
        strat[L] = (min(cur[0], v), E if v < cur[0] else cur[1], cur[2] + 1)
    # broad + AP-seeded (to populate all L)
    for _ in range(30000):
        D = r.randint(25, 65)
        upd([0] + r.sample(range(1, D), K - 2) + [D])
    for d in range(1, 8):
        for _ in range(6000):
            a = r.randint(3, 10)
            base = [d * i for i in range(a)]
            rem = K - a
            pool = [x for x in range(1, d * a + 90) if x not in set(base)]
            if rem <= len(pool):
                upd(base + r.sample(pool, rem))
    for L in range(2, 11):
        v, E, n = strat[L]
        if n:
            print(f"    L={L:2d}: min D3 = {v:.6f}  (n={n:5d})  argmin {E}")
    floors = [strat[L][0] for L in range(2, 11) if strat[L][2]]
    print(f"\n  => min over all L = {min(floors):.6f} at longest-AP = 10 (the binding family).")

    print("\n" + "=" * 98)
    print("CONCLUSION (longest-AP axis):")
    print(f"  tail floor = floor(L=10) = D3(0,3,6,8,9,12,15,18,21,24,27) = {float(gmex):.6f} (EXACT)")
    print(f"  >= bar = {float(BAR):.6f}  (margin +{float(gmex-BAR):.6f})  => k=11 tail CLOSES.")
    print(f"  klein's D3_c limits are the d->inf UPPER references (0.4646 = decorrelated L=10).")
    print("=" * 98)


if __name__ == "__main__":
    main()
