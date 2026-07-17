#!/usr/bin/env python3
"""
THE MULTI-PARENT HUNTER AND THE c = 9/10 CONSECUTIVE CROSSINGS
(boxeph-2026-07-17-S73)

THE MULTI-PARENT HUNTER (LEM-045): for ANY ancestor sets P(i) subset {0..i-1},
    mu(Union_{i<n} A_i) + sum_{i>=1} mu(A_i cap Union_{j in P(i)} A_j)
        <= sum_{i<n} mu(A_i)
-- same leaf-peeling disjointification (Union_{P(n)} A_j subset S); interpolates
union bound (P empty) -> tree-Hunter (|P| = 1) -> EXACT identity (P = all).

Consequence: two-parent credits mu(A_i cap (A_{i-1} u A_{i-2})) =
pair + pair - triple revive the consecutive route past LEM-044's pair-only
boundary: c = 9 needs total credit > 2/7 -- measured here exactly.

Referee: (1) inequality exact on random families/ancestor-sets; (2) the
c = 9 and c = 10 consecutive crossings with 2-/3-parent credits (exact).
"""
import sys, random
from fractions import Fraction as Fr
sys.path.insert(0, '04-computation')
from lrc14_pair_overlap_law_boxeph_S69 import mu_brute, danger_breaks, in_danger

def mu_inter_union(block, i, parents):
    """mu(A_i cap Union_{j in parents} A_j) exact by sweep."""
    speeds = [block[i]] + [block[j] for j in parents]
    bps = sorted(set([Fr(0), Fr(1)] + [x for v in speeds for x in danger_breaks(v)]))
    tot = Fr(0)
    for t in range(len(bps) - 1):
        mid = (bps[t] + bps[t + 1]) / 2
        if in_danger(block[i], mid) and any(in_danger(block[j], mid) for j in parents):
            tot += bps[t + 1] - bps[t]
    return tot

if __name__ == "__main__":
    print("MULTI-PARENT HUNTER + c = 9/10 CONSECUTIVE (boxeph S73)")
    print("=" * 74)
    rng = random.Random(3)
    for trial in range(25):
        c = rng.randint(3, 6)
        fam = sorted(rng.sample(range(2, 40), c))
        P = {i: sorted(rng.sample(range(i), rng.randint(0, min(i, 2))))
             for i in range(1, c)}
        lhs = mu_brute(fam, want_all=False) + sum(
            mu_inter_union(fam, i, P[i]) for i in range(1, c) if P[i])
        rhs = Fr(c, 7)
        assert lhs <= rhs, (fam, P)
    print("  multi-parent inequality exact on 25 random family/ancestor sets")
    print()
    for c, npar in [(9, 2), (10, 2), (10, 3), (11, 3)]:
        ok_all = True
        worst = None
        for v in (1, 5, 20, 50):
            blk = [v + i for i in range(c)]
            credits = sum(mu_inter_union(blk, i,
                          list(range(max(0, i - npar), i)))
                          for i in range(1, c))
            need = Fr(c - 7, 7)
            margin = credits - need
            good = 1 - mu_brute(blk, want_all=False)
            ok = margin > 0 and good >= margin
            ok_all = ok_all and ok
            if worst is None or margin < worst[0]:
                worst = (margin, v)
            if v == 50:
                print(f"  c={c} ({npar}-parent) v=50: credits "
                      f"{float(credits):.5f} vs need {float(need):.5f}; "
                      f"margin {float(margin):+.5f}; good {float(good):.5f}; "
                      f"all v tested cross: {ok_all}; worst margin "
                      f"{float(worst[0]):+.5f} at v={worst[1]}")
    print("=" * 74)
    print("done")
