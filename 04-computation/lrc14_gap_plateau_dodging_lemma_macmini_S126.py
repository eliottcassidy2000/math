#!/usr/bin/env python3
"""Gap-level plateau-dodging lemma — constants referee (mac-mini-S126).

MECHANISM CREDIT: this is kind-pasteur THM-735's simultaneous Bonferroni
multi-peel, run at level 3/41 (gap exclusion) instead of level 0 (LRC
closure), powered by the SETTLED small-n floors instead of a computed body
mass.  Checked before filing that canon has the 1/14-level version
(THM-733/735/738) but — per grep — not the 3/41-level closure statement.

LEMMA (two-far form).  Let S ⊂ {1..13}, |S| = 11.  By settled LRC(12),
M(S) ≥ 1/12.  Let t* be a maximizer; clear_S has slopes ≤ 13, so the
plateau I = component of {clear_S ≥ 3/41} at t* has length
    |I| ≥ ℓ := 2·(1/12 − 3/41)/13 = 5/3198.
For any integers x, y ≥ 2 the danger set {||xt|| < 3/41} meets I in measure
≤ (x|I|+1)·6/(41x) = (6/41)|I| + 6/(41x).  If
    (12/41)|I| + (6/41)(1/x + 1/y) < |I|   ⟸   1/x + 1/y ≤ 29ℓ/6 = 145/19188,
some t ∈ I clears every runner at ≥ 3/41, so M(S∪{x,y}) ≥ 3/41 — NOT in
the gap.  Numerically 2/265 < 145/19188 ⟺ 2·19188 < 265·145 ⟺ 38376 <
38425 ✓: **both fars ≥ 265 ⟹ M ≥ 3/41, for every complement S, no duty
or rung hypotheses.**

MIXED STRIP.  For fixed x, S∪{x} has 12 speeds: M ≥ 1/13 (settled LRC(13));
slope ≤ max(13,x); plateau ℓ₁ ≥ 2(1/13 − 3/41)/max(13,x) = (4/533)/x for
x ≥ 14.  One remaining far y needs (6/41)ℓ₁ + 6/(41y) < ℓ₁ ⟺
y > 6/(35·ℓ₁) = 3198x/70·(6/(35·4/533)=…) — computed exactly below; comes
to y ≥ 23x (integer form y ≥ 23x suffices).  So the ALL-HEIGHTS two-far
gap closure = [S125 sweep, fars ≤ 300] + [both ≥ 265 lemma] + [mixed
strip x ≤ 264, 300 < y < 23x: FINITE, sweepable].

STRATA LADDER.  k fars, |S| = 13−k, M(S) ≥ 1/(14−k) settled for k ≥ 1:
all-big bound B_k from (6k/41)ℓ_k + (6/41)Σ1/x_i < ℓ_k needs 6k < 41 —
DIES AT k = 7 (the apex-7 / MISTAKE-122 / THM-735-j≤6 wall, on cue).

This referee: (a) exact-arithmetic verification of every constant above;
(b) empirical spot-check: 200 random two-far families with fars ≥ 265
(no filters), assert exact M ≥ 3/41 via the S124 verifier; (c) the k-ladder
constants table k = 1..6.
"""

from fractions import Fraction as F
import random
import sys

sys.path.insert(0, "04-computation")
from lrc14_dyadic_tower_ladders_macmini_S124 import exact_M

GAP_HI = F(3, 41)


def main():
    # (a) exact constants
    ell = 2 * (F(1, 12) - F(3, 41)) / 13
    assert ell == F(5, 3198), ell
    thresh = F(29, 41) * ell / (F(6, 41))          # 1/x+1/y bound = 29ℓ/6
    assert thresh == F(145, 19188), thresh
    assert F(2, 265) < thresh and F(2, 264) > thresh
    print(f"(a) two-far: ℓ = {ell}, 1/x+1/y bound = {thresh}; "
          f"both ≥ 265 works (2/265 < {thresh} < 2/264)  EXACT")

    ell1 = 2 * (F(1, 13) - F(3, 41))               # over slope x: ℓ₁ = this/x
    assert ell1 == F(4, 533)
    ybound = 6 / (35 * ell1)                       # y > this·x
    print(f"(a) mixed strip: ℓ₁·x = {ell1}, y > {ybound}·x = "
          f"{float(ybound):.2f}x  ⟹ y ≥ 23x suffices (23 > {float(ybound):.2f})")
    assert 23 > ybound

    print("(a) k-ladder (k fars, settled floor 1/(14-k), slope 13):")
    for k in range(1, 8):
        if 6 * k >= 41:
            print(f"    k={k}: 6k/41 = {F(6*k,41)} ≥ 1 — WALL (apex-7)")
            continue
        ellk = 2 * (F(1, 14 - k) - GAP_HI) / 13
        bound = (F(41 - 6 * k, 41) * ellk) / F(6, 41)   # Σ1/x_i < this
        Bk = int(k / bound) + 1
        print(f"    k={k}: ℓ_k = {ellk}, Σ1/x bound = {bound}, "
              f"all-big B_k = {Bk}")

    # (b) empirical spot-check, NO duty/rung filters, fars ≥ 265
    rng = random.Random(20260719)
    from itertools import combinations
    bad = 0
    for trial in range(200):
        removed = rng.sample(range(1, 14), 2)
        S = [u for u in range(1, 14) if u not in removed]
        x = rng.randint(265, 2000)
        y = rng.randint(265, 2000)
        while y == x:
            y = rng.randint(265, 2000)
        W = sorted(set(S) | {x, y})
        M, t, act, pairs = exact_M(W)
        if M < GAP_HI:
            bad += 1
            print(f"    !! LEMMA VIOLATED: {W}: M = {M}")
    print(f"(b) 200 random both-big (≥265, ≤2000) two-far families: "
          f"{200-bad}/200 have M ≥ 3/41  {'OK' if bad == 0 else 'FAIL'}")
    assert bad == 0

    print("\nLemma constants verified exactly; empirical check clean.")
    print("Status: DERIVED (proof sketch above, THM-735 mechanism at level")
    print("3/41 + settled floors); promotion needs the plateau-interval and")
    print("edge-tooth steps written out — one careful session (or kps).")


if __name__ == "__main__":
    main()
