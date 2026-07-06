#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S20 -- HYP-4482: THE AP-UNIQUENESS as a FREIMAN-TYPE RIGIDITY.

opus-S112: safe(S,beta) = theta-sum over the relation lattice L(S)={a: sum a_i v_i=0}.
safe=0 (the AP, unique tiler) <=> MAXIMAL cancellation <=> RICHEST relation lattice
at low height <=> maximal ADDITIVE ENERGY / minimal DOUBLING <=> the AP (Freiman).

This tests the Freiman frame:
 (1) additive energy E(S)=#{(i,j,k,l): v_i+v_j=v_k+v_l} and doubling |S+S| vs
     safe(2/25): does safe->0 track max E / min |S+S| (the AP)?
 (2) among primitive families, is minimal-doubling (|S+S|=2*12-1=23) UNIQUELY the
     AP => safe=0 <=> minimal doubling <=> AP (Freiman rigidity)?
 (3) the low-height relation count (rank of the theta-sum's dominant terms).
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def safe_measure(W, rho):
    ivals = []
    for v in W:
        r = F(rho) / v
        for m in range(v):
            c = F(m, v)
            ivals.append((c - r, c + r))
    pts = set([F(0), F(1)])
    for lo, hi in ivals:
        for x in (lo, hi):
            xr = x - int(x)
            if xr < 0:
                xr += 1
            pts.add(xr)
    spts = sorted(pts)
    def in_danger(t):
        for v in W:
            y = v * t
            y = y - int(y)
            d = min(y, 1 - y)
            if d < rho:
                return True
        return False
    safe = F(0)
    for i in range(len(spts) - 1):
        a, b = spts[i], spts[i + 1]
        if not in_danger((a + b) / 2):
            safe += (b - a)
    return safe

def additive_energy(W):
    """#{(i,j,k,l): v_i+v_j = v_k+v_l} via sum-frequency."""
    freq = {}
    for i in range(len(W)):
        for j in range(len(W)):
            s = W[i] + W[j]
            freq[s] = freq.get(s, 0) + 1
    return sum(c * c for c in freq.values())

def doubling(W):
    """|S+S| = #distinct pairwise sums (incl i=j)."""
    return len({W[i] + W[j] for i in range(len(W)) for j in range(len(W))})

def low_relations(W, K=1):
    """count nonzero integer relations sum a_i v_i = 0 with |a_i| <= K (the
    dominant theta-sum terms).  K=1: relations with coeffs in {-1,0,1}."""
    n = len(W)
    cnt = 0
    for a in itertools.product(range(-K, K + 1), repeat=n):
        if any(a) and sum(a[i] * W[i] for i in range(n)) == 0:
            cnt += 1
    return cnt

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

RHO = F(2, 25)

def main():
    print("=" * 84)
    print("THE FREIMAN FRAME: safe(2/25) vs additive energy E / doubling |S+S|")
    print("=" * 84)
    AP = tuple(range(1, 13))
    print(f"  AP {{1..12}}: safe={float(safe_measure(AP,RHO)):.5f}, E={additive_energy(AP)}, "
          f"|S+S|={doubling(AP)} (min = 2*12-1 = 23)")
    random.seed(20)
    rows = []
    for _ in range(2500):
        W = primitive(tuple(sorted(random.sample(range(1, 45), 12))))
        if len(set(W)) != 12:
            continue
        s = safe_measure(W, RHO)
        rows.append((s, additive_energy(W), doubling(W), W))
    rows.sort()  # by safe ascending
    print(f"\n  {len(rows)} non-AP families.  Do the LOWEST-safe families have MAX E / MIN |S+S|?")
    print(f"  {'safe':>10} {'E':>7} {'|S+S|':>7}   family")
    for s, E, d, W in rows[:6]:
        print(f"  {float(s):>10.5f} {E:>7} {d:>7}   {list(W)}")
    print(f"  ... (highest-safe end) ...")
    for s, E, d, W in rows[-3:]:
        print(f"  {float(s):>10.5f} {E:>7} {d:>7}   {list(W)}")
    # correlation
    import statistics
    safes = [float(s) for s, E, d, W in rows]
    Es = [E for s, E, d, W in rows]
    ds = [d for s, E, d, W in rows]
    def corr(x, y):
        mx, my = statistics.mean(x), statistics.mean(y)
        num = sum((a - mx) * (b - my) for a, b in zip(x, y))
        den = (sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y)) ** 0.5
        return num / den if den else 0
    print(f"\n  corr(safe, additive energy E) = {corr(safes, Es):+.3f} "
          f"(expect NEGATIVE: low safe <-> high E)")
    print(f"  corr(safe, doubling |S+S|)    = {corr(safes, ds):+.3f} "
          f"(expect POSITIVE: low safe <-> low doubling)")
    # (2) minimal-doubling uniqueness
    print()
    print("  MINIMAL-DOUBLING among primitive 12-families (|S+S| = 23 = the AP's):")
    mind = [(s, W) for s, E, d, W in rows if d == 23]
    print(f"    primitive families with |S+S|=23: {len(mind)} "
          f"({'the AP is UNIQUELY minimal-doubling => Freiman' if all(W==AP for s,W in mind) else 'others exist'})")
    for s, W in mind[:5]:
        print(f"      |S+S|=23: {list(W)} safe={float(s):.5f}")
    print()
    print(f"  AP low-height relations (coeffs in {{-1,0,1}}): {low_relations(AP,1)} "
          f"(the theta-sum's dominant terms; richest => max cancellation => safe=0)")

if __name__ == "__main__":
    main()
