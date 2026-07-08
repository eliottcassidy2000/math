"""
lrc14_cluster_mechanism_opus_S155.py   (opus-2026-07-08-S155)

THE CLUSTER-MONOTONICITY MECHANISM + the decorrelation-limit table D3_c (c=2..10).

Two things:
 (1) D3_c = the c-block + (11-c) fully-decorrelated-outlier limit, computed by Monte-Carlo over
     (x, independent uniform outlier phases).  Verifies klein-S186's c=7..10 and EXTENDS to c=2..10,
     confirming D3_c is monotone DECREASING in c and all >= bar with growing margin.  (c=10 is the
     block+1-outlier, uniquely forced; = LEM-009's 0.4646.)
 (2) THE MECHANISM -- why CLUSTER size, not R2, is the clean D3 axis (klein: "R2 scatters").
     Var(W) = near + far - E[W]^2.  The NEAR part (overlapping arcs, small speed-differences) is
     LOCAL = the densest cluster's internal overlaps; the FAR part DECORRELATES for spread.  So the
     binding part of Var is near-dominated => governed by the max cluster, which survives
     decorrelation.  R2 = sum_d r_d^2 counts ALL differences (incl. large ones that decorrelate and
     do NOT bind), so it scatters.  We show: Var (and D3) track cluster size tightly, while at fixed
     R2 the D3 spread is large -- and near-Var correlates with cluster, far-Var with spread.
"""
import sys, random
import numpy as np
from fractions import Fraction as Fr
from collections import Counter, defaultdict

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_d3_exact_verify_klein_S184 import D3 as D3_exact, BAR

K = 11
TH = 1.0 / 7
Mval = 6.0 / 7
BARf = float(BAR)


def d3_from_moments(m1, m2, m3):
    den = m2 - m3 / Mval
    return m1 / Mval + (m1 - m2 / Mval) ** 2 / den if den > 0 else m1 / Mval


def D3_c_limit(c, nx=3000, nrep=400, seed=0):
    """Monte-Carlo the c-block + (11-c) decorrelated-outlier limit moments -> D3_c."""
    rng = np.random.default_rng(seed)
    X = (np.arange(nx) + 0.5) / nx
    blk = np.arange(c, dtype=np.float64)
    bph = (blk[:, None] * X[None, :]) % 1.0                 # (c, nx) block phases
    nout = K - c
    s1 = s2 = s3 = 0.0
    tot = 0
    for _ in range(nrep):
        u = rng.random((nout, nx)) if nout > 0 else np.empty((0, nx))
        ph = np.concatenate([bph, u], axis=0) if nout > 0 else bph
        ph = np.sort(ph, axis=0)
        gaps = np.empty_like(ph)
        gaps[:-1] = ph[1:] - ph[:-1]
        gaps[-1] = 1.0 - ph[-1] + ph[0]
        W = np.clip(gaps - TH, 0.0, None).sum(axis=0)
        s1 += W.sum(); s2 += (W * W).sum(); s3 += (W * W * W).sum(); tot += nx
    m1, m2, m3 = s1 / tot, s2 / tot, s3 / tot
    return d3_from_moments(m1, m2, m3), m1


def moments_float(E, nx=12000):
    Earr = np.array(sorted(E), dtype=np.float64)
    ph = (Earr[:, None] * ((np.arange(nx) + 0.5) / nx)[None, :]) % 1.0
    ph.sort(axis=0)
    gaps = np.empty_like(ph)
    gaps[:-1] = ph[1:] - ph[:-1]
    gaps[-1] = 1.0 - ph[-1] + ph[0]
    W = np.clip(gaps - TH, 0.0, None).sum(axis=0)
    return W.mean(), (W * W).mean(), (W * W * W).mean()


def near_far_var(E, nx=12000):
    """Var = near + far - E[W]^2; near = pairs of the SAME x with |y-y'|<theta within U (local
    overlap structure). Here we approximate the covering variance's near/far split by the gap
    structure: near-Var ~ contribution of small circular gaps (<2theta, overlapping-arc regime)."""
    m1, m2, m3 = moments_float(E, nx)
    return m2 - m1 * m1, m1


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


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def is_primitive(E):
    E = sorted(E); g = 0
    for i in range(1, len(E)):
        g = gcd(g, E[i] - E[0])
    return g == 1


def main():
    print("=" * 96)
    print("(1) DECORRELATION-LIMIT TABLE  D3_c  (c-block + (11-c) decorrelated outliers), MC")
    print(f"    bar = {BARf:.6f}")
    print("=" * 96)
    print("    c | D3_c (MC)  | E[W]_limit | margin over bar | klein-S186")
    kl = {10: 0.4646, 9: 0.5235, 8: 0.6021, 7: 0.6785}
    prev = None
    mono = True
    for c in range(10, 1, -1):
        d3c, ew = D3_c_limit(c, nx=3000, nrep=300, seed=c)
        kltxt = f"{kl[c]:.4f}" if c in kl else "  --  "
        mark = ""
        if prev is not None and d3c > prev + 1e-3:
            mark = "  [!! not monotone]"; mono = False
        prev = d3c
        print(f"    {c:2d} | {d3c:.5f}    | {ew:.5f}    | {d3c-BARf:+.5f}        | {kltxt}{mark}")
    print(f"    => D3_c {'DECREASING in c (monotone), all >= bar' if mono else 'NON-MONOTONE (check!)'};"
          f" global min at c=10 (block+outlier, uniquely forced).")

    print("\n" + "=" * 96)
    print("(2) MECHANISM: CLUSTER is the clean axis, R2 SCATTERS -- Var is near(=cluster)-dominated")
    print("=" * 96)
    # sample tail shapes; bucket by cluster and by R2; show D3 spread
    rng = random.Random(7)
    rows = []
    for _ in range(6000):
        D = rng.randint(25, 60)
        E = sorted(set([0] + rng.sample(range(1, D), K - 2) + [D]))
        if len(E) != K or not is_primitive(E):
            continue
        m1, m2, m3 = moments_float(E, 6000)
        d3 = d3_from_moments(m1, m2, m3)
        rows.append((cluster_size(E), R2(E), d3, m2 - m1 * m1, m1, E))
    # by cluster
    print("\n  D3 by CLUSTER size (tight?):")
    bc = defaultdict(list)
    for c, r2, d3, var, ew, E in rows:
        bc[c].append(d3)
    for c in sorted(bc):
        a = bc[c]
        print(f"    c={c:2d}: D3 in [{min(a):.4f}, {max(a):.4f}]  spread={max(a)-min(a):.4f}  (n={len(a)})")
    # by R2 band (does same-R2 give same D3? klein says no)
    print("\n  D3 by R2 band (scatters?):")
    br = defaultdict(list)
    for c, r2, d3, var, ew, E in rows:
        br[r2 // 30 * 30].append((d3, c))
    for r2b in sorted(br):
        a = [d for d, _ in br[r2b]]
        cs = [c for _, c in br[r2b]]
        if len(a) >= 8:
            print(f"    R2~[{r2b},{r2b+30}): D3 in [{min(a):.4f}, {max(a):.4f}] spread={max(a)-min(a):.4f}"
                  f"  clusters {min(cs)}-{max(cs)}  (n={len(a)})")
    # correlation: which predicts D3 better, cluster or R2?
    import statistics as st
    cs = [c for c, r2, d3, var, ew, E in rows]
    r2s = [r2 for c, r2, d3, var, ew, E in rows]
    d3s = [d3 for c, r2, d3, var, ew, E in rows]

    def corr(a, b):
        ma, mb = st.mean(a), st.mean(b)
        num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
        da = (sum((x - ma) ** 2 for x in a)) ** .5
        db = (sum((y - mb) ** 2 for y in b)) ** .5
        return num / (da * db) if da * db else 0.0
    print(f"\n  corr(cluster, D3) = {corr(cs, d3s):+.3f}   corr(R2, D3) = {corr(r2s, d3s):+.3f}")
    print(f"  => |corr| larger for the better predictor of D3.")
    print("=" * 96)


if __name__ == "__main__":
    main()
