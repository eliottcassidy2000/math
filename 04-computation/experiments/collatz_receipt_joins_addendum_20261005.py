#!/usr/bin/env python3
"""
Addendum to collatz_receipt_joins_census_20261005.py (opus, 2026-10-05).

(a) product joins broken down by the size of the join point z = x y relative to a, b, ab (hub / mid / upward);
(b) synchronous two-step fusion U^2(ab) = U^2(a) U^2(b) searched to a, b <= 3001 (none expected);
(c) the Applegate--Lagarias multiplier classes (Tables 1-2 of arXiv math/0411140) replayed as parametrized
    receipts: for each class s mod 2^k with a multiplier m inserted at the first odd point y = U(x), the
    orbit endpoint and the moving multiplier label m*y(t); the full insertion defect also includes
    -[y(t)] and the wild-certificate boundary. The plain classes are checked too;
(d) the exact odd density of the multiplier classes (the classes where the universal receipt has a composite defect).
"""
import sys, math, time
from fractions import Fraction

def U(n):
    n = 3 * n + 1
    while n % 2 == 0:
        n //= 2
    return n

def orbit(n):
    out = [n]
    while n != 1:
        n = U(n); out.append(n)
    return out

def main():
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    t0 = time.time()
    # ---- (a)
    Amax = 301
    odds = list(range(3, Amax + 1, 2)); orb = {a: orbit(a) for a in odds}
    pairs = 0; hub = 0; mid = 0; up = 0; mid_ex = []; up_ex = []
    zhist = {}
    for ia, a in enumerate(odds):
        Oa = orb[a]
        for b in odds[ia:]:
            Ob = orb[b]; ab = a * b
            Oab = orbit(ab); pos = {z: i for i, z in enumerate(Oab) if i >= 1 and z != 1}
            pairs += 1
            h = m = u = False
            for j, x in enumerate(Oa):
                if x == 1: break
                for l, y in enumerate(Ob):
                    if y == 1: break
                    z = x * y
                    if z in pos:
                        zhist[z] = zhist.get(z, 0) + 1
                        if z <= max(a, b): h = True
                        elif z < ab:
                            m = True
                            if len(mid_ex) < 12: mid_ex.append((a, b, pos[z], j, l, x, y, z))
                        else:
                            u = True
                            if len(up_ex) < 8: up_ex.append((a, b, pos[z], j, l, x, y, z))
            hub += h; mid += m; up += u
    top = sorted(zhist.items(), key=lambda kv: -kv[1])[:12]
    P(f"(a) pairs {pairs}: product joins with z <= max(a,b) (hub): {hub} ({100*hub/pairs:.1f}%); max(a,b) < z < ab (mid): {mid} ({100*mid/pairs:.1f}%); z >= ab (upward): {up} ({100*up/pairs:.1f}%)")
    P(f"    most frequent join points z (count over all (pair, hit)): {top}")
    P(f"    mid examples (a,b,i,j,l,x,y,z): {mid_ex}")
    P(f"    upward examples: {up_ex}")
    P(f"    [{time.time()-t0:.0f}s]")
    # ---- (b)
    Bmax = 3001
    hits2 = []
    Us = {a: U(a) for a in range(3, Bmax + 1, 2)}
    U2 = {a: U(Us[a]) for a in Us}
    for a in range(3, Bmax + 1, 2):
        ua, u2a = Us[a], U2[a]
        if ua == 1 or u2a == 1: continue
        for b in range(a, Bmax + 1, 2):
            ub, u2b = Us[b], U2[b]
            if ub == 1 or u2b == 1: continue
            ab = a * b
            z = U(U(ab))
            if z == u2a * u2b:
                hits2.append((a, b, z))
    P(f"(b) synchronous two-step fusions U^2(ab) = U^2(a) U^2(b), a <= b <= {Bmax}, all > 1: {len(hits2)} {hits2[:10]}")
    P(f"    [{time.time()-t0:.0f}s]")
    # ---- (c) Applegate--Lagarias classes (Table 1-2, arXiv math/0411140): (residue, modulus, multiplier or 1, where)
    # format: (s, 2^k, [(position_in_T_path, multiplier)...]) -- we replay with the accelerated map: multiply at the
    # listed odd point.  Plain classes: multiplier list empty.
    classes = [
        (1, 4, []), (3, 16, []), (11, 32, []), (27, 128, [(1, 13)]), (91, 256, [(0, 25)]), (219, 256, []),
        (59, 128, []), (123, 256, [(0, 7)]), (251, 256, [(0, 23)]), (7, 64, [(0, 5)]), (39, 128, [(0, 35)]),
        (103, 512, [(2, 13)]), (359, 512, [(0, 35)]), (231, 256, [(0, 5)]), (23, 32, []), (15, 128, []),
        (79, 256, []), (207, 256, [(0, 5), (1, 5)]), (47, 128, [(0, 13)]), (111, 128, [(0, 11)]), (31, 64, [(0, 11)]),
        (63, 128, [(0, 11)]), (127, 256, [(0, 43)]), (255, 512, [(0, 43)]), (511, 1024, [(1, 29)]),
        (1023, 2048, [(0, 11), (1, 11)]), (2047, 4096, [(0, 11), (1, 11)]),
    ]
    P("(c) Applegate--Lagarias classes replayed (accelerated map; multiply by m at the listed odd index, then iterate until below (76/79) x or 12 odd steps):")
    dens_plain = Fraction(0); dens_mult = Fraction(0); ok_all = True
    for (s, M, mults) in classes:
        if s % 2 == 0: continue
        d = Fraction(2, M)
        if mults: dens_mult += d
        else: dens_plain += d
        worst = 0.0; steps_max = 0; good = True; defect_labels = []
        for t in range(0, 64):
            x = s + M * t
            if x == 1: continue
            y = x; k = 0; path = [x]
            mult_at = dict(mults)
            while True:
                if k in mult_at:
                    y = y * mult_at[k]
                    if t < 2: defect_labels.append(y)
                y = U(y); k += 1; path.append(y)
                if 79 * y < 76 * x or k >= 14: break
            r = Fraction(y, x)
            worst = max(worst, r)
            steps_max = max(steps_max, k)
            if r >= 1: good = False
        ok_all &= good
        P(f"    class {s} mod {M}: multipliers {mults or 'none'}; worst ratio over t<64: {float(worst):.4f}; odd steps <= {steps_max}; descent for all t<64: {good}" + (f"; moving multiplier labels (t=0,1): {defect_labels}" if mults else ""))
    P(f"    all listed classes descend (t < 64): {ok_all}")
    P(f"(d) odd density of plain-descent classes: {dens_plain} = {float(dens_plain):.4f}; multiplier classes: {dens_mult} = {float(dens_mult):.4f}; uncovered -1 mod 4096: {Fraction(2,4096)} = {2/4096:.4f}; total {float(dens_plain+dens_mult+Fraction(2,4096)):.4f}")
    # Endpoint and all-prefix coefficient events must not be conflated.
    from math import comb
    endpoint12 = Fraction(sum(comb(19, i) for i in range(12, 20)), 2 ** 19)
    counts = {0: 1}
    for depth in range(1, 13):
        new = {}
        bound = (3**depth).bit_length() - 1
        for old, count in counts.items():
            for total in range(old + 1, bound + 1):
                new[total] = new.get(total, 0) + count
        counts = new
    prefix12 = sum((Fraction(count, 2**total) for total, count in counts.items()), Fraction())
    if endpoint12 != Fraction(11773, 65536) or prefix12 != Fraction(427, 8192):
        raise ValueError("endpoint/prefix control")
    P(f"    comparison: terminal coefficient survival at 12 odd steps = {endpoint12} = {float(endpoint12):.4f}; survival at every prefix = {prefix12} = {float(prefix12):.4f}")
    P(f"    class counts: {len(classes)} total; {sum(not m for _, _, m in classes)} plain; {sum(bool(m) for _, _, m in classes)} multiplier")
    P(f"    [{time.time()-t0:.0f}s total]")
    with open(__file__.replace(".py", ".out"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
