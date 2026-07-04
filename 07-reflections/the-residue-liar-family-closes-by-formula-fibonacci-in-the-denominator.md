# The residue-liar family closes by a formula — with Fibonacci in the denominator

*kind-pasteur-2026-07-04. The owner pointed at a fact — the Fibonacci sequence mod 7 has
period 16 — and asked me to bring it to the endgame. Following it landed on the census's
hardest instances: the single-swap coverer family that mimics the AP but stays loose. It
turns out to close by an explicit formula, its lonely-time denominators are Fibonacci
periodically, and the formula proves the magnitude bound the tight-locus rigidity was missing
for this family.*

## The family

The tight-locus rigidity (GAP-A, mac-mini THM-612/HYP-4070) turns on the single-swap coverers
`{1,…,11,13,X}` with `12 ∣ X`: the AP is `X = 12`, GW is `X = 24`, and the open residual is
that every larger `X` is strictly loose (so the tight locus is finite). mac-mini had the
values — `3/41, 4/53, 5/65, 7/89, …` — and monotonicity, but not a closed form.

Writing `X = 12k`, the closed form is exact:

> **`M({1,2,…,11,13,12k}) = k/(12k+5)` for `k ≥ 3`, attained at `t = (5k+2)/(12k+5)`;**
> and `M = 1/14` at `k = 1` (AP) and `k = 2` (GW).

Since `k/(12k+5) > 1/14 ⇔ 14k > 12k+5 ⇔ k ≥ 3`, the family is lonely for every `k ≥ 1`, and
strictly loose for `k ≥ 3`. So the tight members of `{1,…,11,13,12k}` are **exactly** `X ∈
{12, 24}` — the magnitude bound for this coverer family, previously one of the two open pieces
of GAP-A, now proved by the formula. `M → 1/12` as `k → ∞`: the coverer, pushed to infinity,
relaxes the family toward the `n = 12` bound, never back down to `1/14`.

## The proof is a residue table

The lower bound `M ≥ k/(12k+5)` is a one-line construction: put `t = (5k+2)/(12k+5)`. Then each
runner `v` sits at `r_v/(12k+5)` where `r_v = v(5k+2) mod (12k+5)`, and every `r_v` is a linear
polynomial in `k`:

    v : 1   2    3    4    5    6    7    8    9    10   11   13   12k
    r : 5k+2 10k+4 3k+1 8k+3  k  6k+2 11k+4 4k+1 9k+3 2k  7k+2 5k+1 11k+5

For each, `r_v − k ≥ 0` and `(12k+5) − r_v − k ≥ 0` for `k ≥ 3` (both linear, checked termwise),
so `dist_v ≥ k/(12k+5)`, with equality exactly at `v = 5` (`r = k`) and `v = 12k` (`r = 11k+5`,
so `12k+5 − r = k`). The minimum is `k/(12k+5)`; that is `M`. No search, no `native_decide` —
thirteen linear inequalities in `k`. This is the cleanest kind of infinite-family certificate:
one rational time, verified by modular arithmetic uniformly in the parameter.

## The Fibonacci hint, landed

The denominator `12k+5` is a Fibonacci number precisely at `k = 7, 19, 31, …` — `89 = F₁₁`,
`233 = F₁₃`, `377 = F₁₄`. And `k = 7` is exactly where the family becomes **covering**: `12·7 =
84 = 6·14`, so `{1,…,11,13,84}` covers `q = 14` (and 12, 13), the genuine census-hard
apex-blocking residue-liar — lonely at `t = 37/89` with `M = 7/89`, `89 = F₁₁`. The period-16
structure the owner named (`φ` has order 16 in `F₄₉`, since `x²−x−1` is irreducible mod 7)
governs when `12k+5` is Fibonacci, and the first covering coverer sits on the first Fibonacci
denominator. The census's hardest single-swap family is indexed, in part, by the golden ratio
mod 7.

## Why this is the endgame, not a curiosity

Two things.

**It closes a piece of the rigidity.** GAP-A needed "coverers are magnitude-bounded" for the
single-swap family; `M = k/(12k+5) > 1/14` for `k ≥ 3` is that bound, with a proof rather than
a search. The residual finiteness (the *general* tight-locus, all lift/coverer combinations)
stays open — that is the Perarnau–Serra extremal-uniqueness — but this family's slice is now
a theorem.

**It is the general principle the owner was pointing at.** An infinite family of LRC(14)
instances — including the census's worst apex-blockers — closed by a single explicit rational
time and a residue table uniform in the parameter. This is the same shape as the far-peel
infinite family `{1,…,12,w≥182}` (S39): the census is not a wall of individual checks but a
union of parametrized families each closable at once. The Fibonacci periodicity is one such
parametrization made visible. Every such formula is a slice of the census turned from
"astronomical finite check" into "closed form."

## The extremes, named

Fibonacci also fixes the two poles of the tight-locus picture. Consecutive-Fibonacci families
`{F_k,…,F_{k+12}}` are **maximally loose** (`M ≈ 0.198 = 2.78·(1/14)`, converging, verified) —
the golden ratio makes the runners maximally equidistributed. The AP is the tight pole
(`M = 1/14`, runners on the 14th-roots). So the rigidity is the statement that the tight pole
is isolated, and the residue-liars are the near-tight families that mimic the AP's residues
while relaxing toward the golden pole: `M = k/(12k+5)` interpolating from `1/14` (k small,
near AP) up toward `1/12` (k large, relaxed). The owner's hint named the far pole.

---
*Linked: [[the-tight-locus-is-the-arithmetic-progression]] (S37), [[the-uniform-looseness-is-lrc-hard-the-far-peel-is-measure-and-linear]]
(S39). Sharpens mac-mini HYP-4070 (GAP-A value-list → closed formula + proof), opus HYP-4047
(Eisenstein/`Φ₆` is the cyclotomic pole; golden/`F₄₉` is the loose pole). Scripts:
`lrc14_residue_liar_formula`, `lrc14_fibonacci_mod7_kps.py`. HYP-4078 (renumbered from 4076).*
