#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_add_energy_law_refine_macmini.py  (mac-mini 2026-06-21, THREAD D, lead ii)

REFINE the additive-energy law (HYP-2735). The honest gap: AE alone has Kendall
concordance ~0.72 with measS7 and FAILS to separate equal-AE sets (consec 0.345 vs
odd-AP 0.309). The IE fingerprint (lrc14_consec_vs_oddAP_IE) showed the consec
advantage is carried by the EVEN inclusion-exclusion terms S_2, S_4, S_6 (esp. j=4)
-- exactly the even-Krawtchouk band the Delsarte dual lives on (HYP-2726).

NEW HYPOTHESIS (HYP-2737 candidate):
  measS7(E) is governed by the SIGNED EVEN-IE functional
      Phi(E) := S_0 - S_2 + S_4 - S_6   (the even part of the IE alternating sum)
  better than by additive energy. Equivalently, since p0 = sum (-1)^j S_j and the
  ODD part contributes the "easy" coupon floor, the CORRECTION corr = p0 - coupon is
  dominated by the even-band imbalance.

TESTS:
  (1) compute Phi(E) and even-only / odd-only IE partial sums; rank concordance with p0.
  (2) test the cleaner candidate: does the EVEN-Krawtchouk moment  E[K_2(N)]  (N = #empty
      sectors, K_2 the binary Krawtchouk) order p0 monotonically?  (K_2 is convex, so
      E[K_2(N)] is Schur-convex -> consec-max would follow from a convex-order argument).
  (3) does (AE, then Phi) two-key order beat (AE, then span)?  -- the real tie-break.

If a SINGLE even-band scalar orders p0 monotonically, that's a non-finite-check route
to consec-max (the convex-order / Schur-convexity argument, no exhaustive atlas).
"""
import itertools
from fractions import Fraction as Fr
from math import gcd, comb
from collections import defaultdict, Counter

P = 7
def sector(yf): return int(P * yf)

def breakpoints(E):
    bp = {Fr(0), Fr(1)}
    for e in E:
        if e == 0: continue
        for t in range(0, P * e):
            bp.add(Fr(t, P * e))
    return sorted(bp)

def empty_count_law(E):
    """dict: N (#empty sectors) -> Lebesgue measure."""
    E = [int(e) for e in E if int(e) != 0]
    xs = breakpoints(E); law = defaultdict(lambda: Fr(0))
    for a, b in zip(xs, xs[1:]):
        mid = (a + b) / 2
        N = P - len(set(sector((e * mid) % 1) for e in E))
        law[N] += (b - a)
    return law

def S_terms(E):
    law = empty_count_law(E)
    S = [Fr(0)] * (P + 1)
    for N, m in law.items():
        for j in range(N + 1):
            S[j] += comb(N, j) * m
    return S

def measS7(E):
    S = S_terms(E)
    return sum((-1)**j * S[j] for j in range(P + 1))

def add_energy(E):
    r = Counter()
    for a in E:
        for b in E: r[a + b] += 1
    return sum(v * v for v in r.values())

# binary Krawtchouk K_j(x; m), m = P-1 = 6 (the dual space dimension in HYP-2726)
def krawtchouk(j, x, m=6):
    return sum((-1)**i * comb(x, i) * comb(m - x, j - i) for i in range(j + 1))

def Kmoment(E, j):
    """E_v[ K_j(N) ] over the empty-count law."""
    law = empty_count_law(E)
    return sum(krawtchouk(j, N) * m for N, m in law.items())

def concordance(rows, keyfn, dirn=+1):
    conc = disc = 0
    n = len(rows)
    for i in range(n):
        for jj in range(i+1, n):
            ti, tj = rows[i]["p0"], rows[jj]["p0"]
            if ti == tj: continue
            si, sj = keyfn(rows[i]), keyfn(rows[jj])
            if si == sj: continue
            if (si - sj) * (ti - tj) * dirn > 0: conc += 1
            else: disc += 1
    return conc / (conc + disc) if (conc + disc) else float('nan')

def main():
    print("#"*80)
    print("# ADDITIVE-ENERGY LAW REFINEMENT (THREAD D, lead ii)")
    print("#"*80)

    # population: primitive 8-sets min=1 max<=14
    rows = []
    for combo in itertools.combinations(range(2, 15), 7):
        E = (1,) + combo
        g = 0
        for e in E: g = gcd(g, e)
        if g != 1: continue
        S = S_terms(E)
        p0 = sum((-1)**j * S[j] for j in range(P + 1))
        rows.append({
            "E": E, "p0": p0,
            "AE": add_energy(E),
            "even_IE": S[0] - S[2] + S[4] - S[6],
            "odd_IE": S[1] - S[3] + S[5] - S[7],
            "S2": S[2], "S4": S[4],
            "K2": Kmoment(E, 2), "K4": Kmoment(E, 4), "K6": Kmoment(E, 6),
            "span": max(E) - min(E),
        })
    print(f"\npopulation size = {len(rows)}")

    print("\n=== Kendall concordance with measS7 (1.0 = perfect monotone order) ===")
    cands = [
        ("AE",      lambda r: r["AE"],       +1),
        ("even_IE", lambda r: r["even_IE"],  +1),
        ("odd_IE",  lambda r: r["odd_IE"],   +1),
        ("S2",      lambda r: r["S2"],       -1),
        ("S4",      lambda r: r["S4"],       +1),
        ("K2(N)",   lambda r: r["K2"],       +1),
        ("K4(N)",   lambda r: r["K4"],       +1),
        ("K6(N)",   lambda r: r["K6"],       +1),
        ("-K2(N)",  lambda r: r["K2"],       -1),
        ("span",    lambda r: r["span"],     -1),
    ]
    out = []
    for name, fn, d in cands:
        out.append((concordance(rows, fn, d), name, d))
    for c, name, d in sorted(out, reverse=True):
        print(f"  {name:9s} dir={d:+d}  concordance={c:.4f}")

    # two-key orderings: does (AE, even_IE) reproduce the p0 order's TOP better than (AE,span)?
    print("\n=== two-key tie-break: predicted #1 by (AE desc, KEY) ===")
    true_top = sorted(rows, key=lambda r: -r["p0"])[0]
    print(f"  TRUE #1: p0={float(true_top['p0']):.5f} E={true_top['E']}")
    for keyname, keyfn in [("even_IE desc", lambda r: -float(r["even_IE"])),
                            ("K2 desc",      lambda r: -float(r["K2"])),
                            ("span asc",     lambda r:  r["span"]),
                            ("S4 desc",      lambda r: -float(r["S4"]))]:
        pred = sorted(rows, key=lambda r: (-r["AE"], keyfn(r)))[0]
        print(f"  (AE, {keyname:13s}) -> #1 E={pred['E']} p0={float(pred['p0']):.5f}  "
              f"{'MATCH consec' if pred['E']==tuple(range(1,9)) else 'MISS'}")

    # the equal-AE separation: among the max-AE sets, which scalar separates consec from odd-AP?
    print("\n=== equal-AE class (max AE) ranked by p0; even_IE & K2 alongside ===")
    maxAE = max(r["AE"] for r in rows)
    cls = sorted([r for r in rows if r["AE"] == maxAE], key=lambda r: -r["p0"])
    print(f"  max AE = {maxAE}, #sets = {len(cls)}")
    for r in cls:
        print(f"    p0={float(r['p0']):.5f} even_IE={float(r['even_IE']):.4f} "
              f"K2={float(r['K2']):.4f} span={r['span']} E={r['E']}")

    print("\nDONE.")

if __name__ == "__main__":
    main()
