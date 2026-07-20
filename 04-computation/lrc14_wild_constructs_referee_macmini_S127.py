#!/usr/bin/env python3
"""Wild-constructs referee — HYP-8010 / 8015 / 8025 (+8030 note)
(mac-mini-2026-07-19-S127; owner: hypothesize wildly, combine ideas).

(A) HYP-8010 — THE LONELINESS TRAJECTORY.  Tournament = opus-S407's
    observer-inclusive half-turn tournament on {0} ∪ W (14 vertices):
    i → j  iff  ((w_j − w_i)·t mod 1) < 1/2.
    It is piecewise-constant in t, changing only on the DIFFERENCE lattice
    {k/d, (2k+1)/(2d) : d = pairwise differences}.  S407 verified
    c3(t*) ∈ {108,112} for near-floor families (r = −0.926 vs M).
    THE ADVERSARIAL QUESTION S407's endpoint cannot answer: what is the
    c3(t) PROFILE over the whole circle?  If c3 ≈ 112 on most of the
    circle, the endpoint correlation is a base rate, not steering; if the
    near-regular region is small and t* sits inside it, the bridge is a
    genuine trajectory phenomenon.  We compute: measure of {t : c3(t) =
    max}, c3 at t*, quantile of c3(t*) in the profile.

(B) HYP-8015 — THE DYADIC DESCENT CHAIN.  σ(W) = evens(W)/2.  THEOREM
    (one line, verified here): M(W) ≤ M(σW) (the even sub-family's
    clearance at t equals the halved family's at 2t; dropping runners
    raises the max-min).  The chain (M(σ^j W))_j vs the AP ceiling
    reference 1/⌈14/2^j⌉ = a new near-floor invariant; deviations locate
    the 2-adic depth of each family's defect.

(C) HYP-8025 — DUTY RIGIDITY + DEFICIT LAPLACIAN.
    (a) duty table q ∈ {2..14} → carriers in W; load count = ∏ #carriers
        (0 if a duty is uncovered = non-covering there); exclusive count =
        permanent of the duty×carrier incidence (injective assignments).
        Wild claim under test: extremal families are duty-degenerate
        (stacked duties force permanent collapse).
    (b) deficit Laplacian: weights δ(a,b) = 1/(2(b²−a²)) on mixed-parity
        pairs (HYP-7960), 0 on both-odd; Fiedler value λ₂ census vs M.

(D) HYP-8030 — CB filtration: no computation here (conjecture + in-range
    citations); see the INDEX entry.
"""

from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd
import sys

sys.path.insert(0, "04-computation")
from lrc14_dyadic_tower_ladders_macmini_S124 import exact_M

ONE14 = F(1, 14)


# ------------------------------------------------------------------ (A)
def c3_profile(W):
    """c3(t) of the observer-inclusive half-turn tournament, exact
    piecewise census over [0,1)."""
    V = [0] + sorted(W)
    n = len(V)
    diffs = sorted({abs(a - b) for a, b in combinations(V, 2)})
    pts = set()
    for d in diffs:
        for k in range(2 * d):
            pts.add(F(k, 2 * d))     # k/d and half-shifted (2k+1)/(2d)
    pts = sorted(pts)
    pts.append(F(1))
    segs = []                        # (length, c3)
    for a, b in zip(pts, pts[1:]):
        if a == b:
            continue
        t = (a + b) / 2
        out = [0] * n
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                fr = ((V[j] - V[i]) * t) % 1
                if fr < F(1, 2):
                    out[i] += 1
        c3 = (n * (n - 1) * (n - 2)) // 6 - sum(o * (o - 1) // 2 for o in out)
        segs.append((b - a, c3))
    return segs


def analyze_trajectory(name, W, Mval=None, tstar=None):
    if tstar is None:
        Mval, tstar, _, _ = exact_M(sorted(W))
    segs = c3_profile(W)
    total = sum(l for l, _ in segs)
    assert total == 1
    cmax = max(c for _, c in segs)
    meas_max = sum(l for l, c in segs if c == cmax)
    # c3 at t*: find the segment containing t* (t* may be a breakpoint —
    # evaluate just left and right)
    eps = F(1, 10 ** 9)
    for probe, side in ((tstar - eps, "left"), (tstar + eps, "right")):
        pass
    def c3_at(t):
        V = [0] + sorted(W)
        n = len(V)
        out = [0] * n
        for i in range(n):
            for j in range(n):
                if i != j and ((V[j] - V[i]) * t) % 1 < F(1, 2):
                    out[i] += 1
        return (n*(n-1)*(n-2))//6 - sum(o*(o-1)//2 for o in out)
    c3s = c3_at(tstar + eps)
    below = sum(l for l, c in segs if c < c3s)
    print(f"  {name:<22} M={str(Mval):>8} c3(t*+)={c3s:>3}  max c3={cmax:>3} "
          f"on measure {float(meas_max):.4f}  quantile(c3(t*))={float(below):.4f}")
    return c3s, cmax, meas_max


# ------------------------------------------------------------------ (B)
def sigma(W):
    return sorted(w // 2 for w in W if w % 2 == 0)


def descent_chain(name, W):
    chain = []
    cur = sorted(W)
    n0 = 14
    j = 0
    while cur:
        Mv, t, _, _ = exact_M(cur) if len(cur) > 1 else (F(1, 2), F(1, 2), 0, 0)
        ref = F(1, -(-n0 // 2 ** j))       # 1/ceil(14/2^j)
        chain.append((j, len(cur), Mv, ref, Mv == ref))
        nxt = sigma(cur)
        if nxt == cur or not nxt:
            break
        # verify the descent inequality at this step
        cur = nxt
        j += 1
    # inequality check across chain
    ok = all(chain[i][2] <= chain[i + 1][2] for i in range(len(chain) - 1))
    marks = " ".join(f"{str(m)}{'*' if hit else ''}" for _, _, m, _, hit in chain)
    print(f"  {name:<22} chain: {marks}   descent-ineq OK: {ok}")
    assert ok
    return chain


# ------------------------------------------------------------------ (C)
def duty_table(W):
    tab = {}
    for q in range(2, 15):
        tab[q] = [w for w in W if w % q == 0]
    return tab


def permanent01(rows):
    """Permanent of a small 0-1 matrix given as list of column-sets."""
    n = len(rows)
    if n == 0:
        return 1
    total = 0
    cols = rows[0]
    for c in cols:
        rest = [r - {c} for r in rows[1:]]
        if any(not r for r in rest):
            # zero row → contributes 0 unless handled; recursion handles
            pass
        total += permanent01(rest)
    return total


def duty_census(name, W):
    tab = duty_table(W)
    covered = [q for q in range(2, 15) if tab[q]]
    load = 1
    for q in covered:
        load *= len(tab[q])
    # exclusive permanent over the covered duties with ≥1 carrier,
    # restricted to duties 7..14 (small multiples; duties ≤6 have many
    # carriers and swamp the permanent) — the STACKING signal lives high.
    high = [set(tab[q]) for q in range(7, 15) if tab[q]]
    perm = permanent01(high) if high else 0
    print(f"  {name:<22} covered q∈[2,14]:{len(covered):>2}  load={load:<12} "
          f"high-duty(7..14) exclusive assignments={perm}")
    return perm


def fiedler(name, W):
    Ws = sorted(W)
    n = len(Ws)
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a, b = Ws[i], Ws[j]
            g = gcd(a, b)
            aa, bb = a // g, b // g
            w = 0.0 if (aa % 2 and bb % 2) else 1.0 / (2 * (bb * bb - aa * aa))
            L[i][i] += w
            L[j][j] += w
            L[i][j] -= w
            L[j][i] -= w
    # smallest nonzero eigenvalue by power iteration on pseudo-inverse-ish:
    # simple Jacobi eigen via numpy-free — use characteristic through
    # deflation: cheap route: numpy if available
    try:
        import numpy as np
        ev = sorted(np.linalg.eigvalsh(np.array(L)))
        lam2 = ev[1]
    except Exception:
        lam2 = float("nan")
    print(f"  {name:<22} deficit-Laplacian lambda2 = {lam2:.3e}")
    return lam2


def main():
    fams = [
        ("AP (tight)",        list(range(1, 14))),
        ("GW (tight)",        list(range(1, 12)) + [13, 24]),
        ("deep well",         list(range(1, 13)) + [182]),
        ("F3 {1..11,13,36}",  list(range(1, 12)) + [13, 36]),
        ("K2 {1..12,26}",     list(range(1, 13)) + [26]),
        ("{1..12,14}",        list(range(1, 13)) + [14]),
        ("translate {2..14}", list(range(2, 15))),
        ("BF38 band",         [3, 5, 7, 8, 9, 10, 11, 12, 13, 15, 21, 35]),
        ("resistant S125",    [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 108]),
        ("powers of 2",       [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]),
    ]

    print("== (A) HYP-8010: trajectory profiles (observer half-turn c3(t)) ==")
    print("   [question: is c3(t*) high because the whole profile is high?]")
    for name, W in fams[:7]:
        analyze_trajectory(name, W)

    print("\n== (B) HYP-8015: dyadic descent chains (M(σ^j W)); * = equals 1/⌈14/2^j⌉ ==")
    for name, W in fams:
        descent_chain(name, W)

    print("\n== (C) HYP-8025a: duty rigidity census ==")
    for name, W in fams:
        duty_census(name, W)

    print("\n== (C) HYP-8025b: deficit Laplacian Fiedler census ==")
    for name, W in fams:
        fiedler(name, W)

    print("\nDone; interpretations to the INDEX entries.")


if __name__ == "__main__":
    main()
