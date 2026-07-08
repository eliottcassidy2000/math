"""
lrc14_cluster_monotonicity_opus_S155.py   (opus-2026-07-08-S155)

WORK THE CLUSTER-MONOTONICITY STEP (klein-S186 / LEM-009) for the k=11 covering tail.

klein reduced the prim-diam >= 25 tail to: every such shape has D3 >= D3_{c(E)} >= D3_10 = 0.4646,
c(E) = max points in a length-9 window, D3_c the c-block+(11-c)-decorrelated-outlier limit
(D3_10/9/8/7 = 0.4646/0.5235/0.6021/0.6785, decreasing).  The make-or-break for the whole closure
is the EXTREMAL claim: is the 10-block+outlier the GLOBAL D3-minimizer over ALL prim-diam >= 25
shapes (klein: min D3 = 0.4587 at {0..9,25}), or can some structure go lower?

STRUCTURAL SIMPLIFICATION (opus): c = 10 (10 integer points in a span-9 window) FORCES the block
{j..j+9} (only 10 integers fit in [j,j+9]).  So the densest cluster is UNIQUELY the block.  The
extremal claim splits: (i) c=10 = block + 1 outlier, all diam -> min D3 = 0.4587 [klein LEM-009];
(ii) c <= 9 -> D3 > 0.4587 [the monotonicity, must survive a hard search].

METHOD: fast numpy float D3 via W(x) = sum_i (circular_gap_i - 1/7)_+ (the uncovered measure), on
a fine x-grid, for the SEARCH; klein-S184 EXACT rational D3 to verify any near-minimal candidate.
"""
import sys, random
import numpy as np
from fractions import Fraction as Fr
from collections import Counter, defaultdict

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_d3_exact_verify_klein_S184 import D3 as D3_exact, BAR, moments

K = 11
TH = 1.0 / 7
Mval = 6.0 / 7
BARf = float(BAR)
NGRID = 9000
_X = (np.arange(NGRID) + 0.5) / NGRID          # midpoint grid, avoids x=0 alias


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def is_primitive(E):
    E = sorted(E); g = 0
    for i in range(1, len(E)):
        g = gcd(g, E[i] - E[0])
    return g == 1


def prim_diam(E):
    E = sorted(E); g = 0
    for i in range(1, len(E)):
        g = gcd(g, E[i] - E[0])
    return (E[-1] - E[0]) // max(g, 1)


def cluster_size(E, W=9):
    E = sorted(E); best = 1
    for i in range(len(E)):
        best = max(best, sum(1 for e in E if E[i] <= e <= E[i] + W))
    return best


def R2(E):
    c = Counter()
    for i in range(len(E)):
        for j in range(len(E)):
            if i != j:
                c[E[i] - E[j]] += 1
    return sum(v * v for v in c.values())


def D3_float(E):
    """Fast float D3 via W(x)=sum(gap-theta)_+ on the x-grid."""
    Earr = np.array(sorted(E), dtype=np.float64)
    ph = (Earr[:, None] * _X[None, :]) % 1.0          # (k, N)
    ph.sort(axis=0)
    gaps = np.empty_like(ph)
    gaps[:-1] = ph[1:] - ph[:-1]
    gaps[-1] = 1.0 - ph[-1] + ph[0]
    W = np.clip(gaps - TH, 0.0, None).sum(axis=0)      # uncovered measure per x
    m1 = W.mean(); m2 = (W * W).mean(); m3 = (W * W * W).mean()
    den = m2 - m3 / Mval
    if den <= 0:
        return m1 / Mval
    return m1 / Mval + (m1 - m2 / Mval) ** 2 / den


def main():
    print("=" * 100)
    print("CLUSTER-MONOTONICITY / GLOBAL EXTREMALITY hard-test for the k=11 prim-diam>=25 tail")
    print(f"  bar = {BARf:.6f};  float-D3 grid N={NGRID}")
    print("=" * 100)

    claim = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 25)
    d3c_ex = D3_exact(claim)
    d3c_fl = D3_float(claim)
    print(f"\n(A) claimed global min {claim}:")
    print(f"    EXACT  D3 = {float(d3c_ex):.6f} (margin {float(d3c_ex-BAR):+.6f}); float D3 = {d3c_fl:.6f}"
          f"  [float err {abs(d3c_fl-float(d3c_ex)):.2e}]  c={cluster_size(claim)} R2={R2(claim)}")
    THRESH = float(d3c_ex)

    best = (d3c_fl, claim, 10, 590)
    near_thresh = []          # shapes with float-D3 within 0.01 of THRESH (exact-verify later)

    def upd(E):
        nonlocal best
        E = tuple(sorted(set(E)))
        if len(E) != K or not is_primitive(E) or prim_diam(E) < 25:
            return
        d3 = D3_float(E)
        if d3 < best[0]:
            best = (d3, E, cluster_size(E), R2(E))
        if d3 < THRESH + 0.01:
            near_thresh.append((d3, E))

    # (B) broad random search
    print("\n(B) broad random tail search (diam 25..70) ...")
    rng = random.Random(20260708)
    ntot = 0
    for _ in range(22000):
        D = rng.randint(25, 70)
        E = [0] + rng.sample(range(1, D), K - 2) + [D]
        Es = sorted(set(E))
        if len(Es) != K or not is_primitive(Es):
            continue
        ntot += 1
        upd(Es)
    print(f"    scanned {ntot} primitive shapes; random min float-D3 = {best[0]:.6f} at {best[1]} (c={best[2]},R2={best[3]})")

    # (C) structured adversaries
    print("\n(C) STRUCTURED adversaries (most likely to beat block+outlier):")

    def scan(name, gen):
        b = None; n = 0
        for E in gen:
            Es = sorted(set(E))
            if len(Es) != K or not is_primitive(Es) or prim_diam(Es) < 25:
                continue
            n += 1
            d3 = D3_float(Es)
            upd(Es)
            if b is None or d3 < b[0]:
                b = (d3, tuple(Es), cluster_size(Es), R2(Es))
        if b:
            print(f"    {name:40s}: min D3={b[0]:.6f} at {b[1]} (c={b[2]},R2={b[3]})  [n={n}]")

    r = random.Random(1)
    for c in (10, 9, 8, 7, 6):
        def bo(c=c):
            for _ in range(2500):
                outs = sorted(r.sample(range(c + 1, 100), K - c))
                yield list(range(c)) + outs
        scan(f"{c}-block + {K-c} outliers", bo())

    scan("block_10 + outlier D (D=25..400)", ([*range(10), D] for D in range(25, 401)))

    for (a, b) in [(5, 5), (6, 4), (6, 5), (7, 4), (5, 4), (8, 3), (7, 3), (4, 4)]:
        def tb(a=a, b=b):
            for _ in range(2500):
                off = r.randint(a + 2, 45)
                blk = list(range(a)) + [off + i for i in range(b)]
                rem = K - a - b
                pool = [x for x in range(90) if x not in set(blk)]
                outs = sorted(r.sample(pool, rem)) if rem > 0 else []
                yield sorted(blk + outs)
        scan(f"two blocks {a}+{b} + {max(0,K-a-b)} out", tb())

    for d in (2, 3, 4, 5, 6):
        def na(d=d):
            base = [i * d for i in range(K)]
            for _ in range(2000):
                E = base[:]
                j = r.randrange(1, K); E[j] += r.choice([-1, 1])
                yield E
        scan(f"near-AP d={d} (perturbed prim)", na())

    # AP with a few points removed + replaced far (near-extremal energy)
    def ap_variants():
        for d in (1, 2, 3):
            for _ in range(1500):
                base = [i * d for i in range(K + 2)]
                drop = sorted(r.sample(range(1, K + 1), 2))
                E = [x for idx, x in enumerate(base) if idx not in drop]
                E = E[:K - 1] + [base[-1] + r.randint(1, 40)]
                yield E
    scan("AP-minus-2 + far point", ap_variants())

    # (D) axis test
    print("\n(D) AXIS TEST: min float-D3 stratified by cluster size c:")
    strat = defaultdict(lambda: [9.0, None, 0])
    rr = random.Random(999)
    for _ in range(18000):
        D = rr.randint(25, 60)
        E = sorted(set([0] + rr.sample(range(1, D), K - 2) + [D]))
        if len(E) != K or not is_primitive(E):
            continue
        d3 = D3_float(E); c = cluster_size(E)
        s = strat[c]; s[2] += 1
        if d3 < s[0]:
            s[0] = d3; s[1] = tuple(E)
    for c in sorted(strat):
        s = strat[c]
        print(f"    c={c:2d}: min D3={s[0]:.6f}  (n={s[2]:5d})  argmin {s[1]}")

    # verify near-threshold candidates EXACTLY
    print("\n(E) EXACT verification of all float-candidates within 0.01 of the claimed min:")
    seen = set()
    cands = sorted(set(E for _, E in near_thresh))
    below = []
    for E in cands:
        if E in seen:
            continue
        seen.add(E)
        d3e = D3_exact(E)
        tag = ""
        if d3e < BAR:
            tag = "  *** BELOW BAR ***"
        if d3e < d3c_ex:
            tag += "  <<< BEATS CLAIMED MIN"
            below.append((float(d3e), E))
        if float(d3e) < THRESH + 0.003:   # only print the closest
            print(f"    {E}: exact D3={float(d3e):.6f} (margin {float(d3e-BAR):+.6f}){tag}")

    print("\n" + "=" * 100)
    print(f"GLOBAL min float-D3 = {best[0]:.6f} at {best[1]} (c={best[2]}, R2={best[3]})")
    print(f"claimed min block+outlier {claim}: exact D3 = {float(d3c_ex):.6f}")
    if below:
        print(f"!!! {len(below)} shapes BEAT the claimed min in EXACT arithmetic:")
        for v in sorted(below)[:10]:
            print(f"      D3={v[0]:.6f}  {v[1]}")
    else:
        print(f"NO shape beats the claimed min (exact-verified). block+outlier {claim} SURVIVES.")
        print(f"  => every prim-diam>=25 shape searched has D3 >= {float(d3c_ex):.6f} >= bar (+{float(d3c_ex-BAR):.4f}).")
    print("=" * 100)


if __name__ == "__main__":
    main()
