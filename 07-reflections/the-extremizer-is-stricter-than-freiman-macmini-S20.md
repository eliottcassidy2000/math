---
source: mac-mini-2026-07-06-S20
status: synthesis + verified refinement of the AP-uniqueness frontier
tags:
  - lonely-runner
  - density-floor
  - freiman
  - additive-energy
  - doubling
  - relation-lattice
  - theta-sum
  - residue-pinning
---

# The LRC extremizer is STRICTER than the Freiman optimum

opus-S112 gave `safe(S,β)` an exact form: a theta-sum over the relation lattice
`L(S) = {a : Σ a_i v_i = 0}`, `safe = Σ_{a∈L} ∏_i ĥ_{a_i}(β)`. `safe = 0` (the
AP, unique tiler) ⟺ maximal cancellation. The natural additive-combinatorics
reading — richest lattice = most relations = max additive energy = min doubling
= AP (Freiman) — is HALF right, and the half that fails is the interesting part.

## What the Freiman frame gets right

- **The AP is uniquely minimal-doubling.** `|S+S| ≥ 2n−1` for any n-set (sort;
  `a_1+a_1 < … < a_1+a_n < a_2+a_n < … < a_n+a_n` are 2n−1 distinct sums), with
  equality iff S is an arithmetic progression. Among PRIMITIVE 12-families, only
  the AP {1,…,12} attains `|S+S| = 23` (verified: 0 non-AP primitives).
- **The AP maximizes additive energy** (E = 1156) and has ~9910 low-height
  relations (coeffs in {−1,0,1}) — the dominant theta-sum terms.

## What it gets wrong (the refinement)

**`safe = 0` correlates only WEAKLY with additive energy / doubling over random
families** (corr ≈ −0.03, +0.07). The reason: the theta-sum is a SIGNED
cancellation, not a count. Many relations do not imply cancellation to 0 — the
signs must conspire. So "most relations" is necessary, not sufficient.

**Sharper: `safe = 0` is STRICTER than minimal doubling.** Shifted APs
`{a, a+d, …, a+11d}` with `a ≠ d` are minimal-doubling (`|S+S| = 23`) yet have
`safe > 0`:

| family | \|S+S\| | M | safe(2/25) |
|---|---|---|---|
| dilated AP c·{1,…,12} | 23 | **1/13** | **0** |
| {2,…,13} (shift) | 23 | 2/15 | 0.062 |
| {1,3,5,…,23} (d=2) | 23 | 1/2 | 0.112 |
| {3,…,14} (shift) | 23 | 3/17 | 0.103 |

So Freiman-optimality (min doubling) is achieved by EVERY AP, but LRC-extremality
(`safe = 0`, `M = 1/13`) only by the **dilated AP c·{1,…,12}** — first term =
common difference.

## Why: the extra rigidity is residue-pinning (a ≡ d mod n)

A shifted AP `{a+id : i=0..11}` reduces mod 13 to a complete-minus-one residue
system; it MISSES the residue 0 (⟺ residues = {1,…,12}, the tight structure) iff
`a ≡ d (mod 13)`. The dilated AP c·{1,…,12} has a = d = c, so `a ≡ d` and its
residues are exactly {1,…,12} — residue_pinning_13. A generic shifted AP has
`a ≢ d`, hence INCLUDES a multiple of 13, breaking the tight structure.

So the LRC extremizer sits at the intersection of TWO minimalities:
- **additive** (Freiman): minimal doubling `|S+S| = 2n−1` ⟹ AP;
- **multiplicative/residue** (residue_pinning): complete unit residues mod 13
  ⟹ `a ≡ d`.

Their intersection is the single ray `c·{1,…,12}`. This is the sum-product
coincidence (opus-S107) made exact: **the AP is the unique family that is BOTH
Freiman-minimal AND residue-complete.** Neither alone suffices; loneliness needs
both. The theta-sum's full cancellation is the analytic shadow of this double
minimality — the signed relations cancel only when the additive AP-structure and
the mod-13 completeness align, which pins first-term = difference.

## Consequence for the (U)/(G) proof

The (U) rigidity (`safe = 0 ⟹ dilated AP`) factors:
1. `safe = 0 ⟹ additive minimality` (`|S+S| = 2n−1`) — the theta-sum /
   additive-energy step (OPEN; the signed-cancellation-forces-min-doubling
   lemma);
2. `|S+S| = 2n−1 ⟹ AP` — CLASSICAL (formalizable; the minimal-sumset theorem);
3. `+ residue-pinning ⟹ a ≡ d ⟹ dilated AP` — GREEN (`residue_pinning_13`).

So the open heart of (U) is step 1 alone — "full theta-sum cancellation forces
minimal doubling" — a clean additive-combinatorics target, sharper than raw
extremizer-uniqueness, with steps 2 and 3 either classical or already green.

## Formalized

`LRCMinimalSumset.lean` (this session): the classical `|S+S| ≥ 2n−1` for a
finite set of integers (the Freiman base, step 2's lower bound) — self-contained,
the additive-combinatorics anchor of the frame.

-> HYP-4482, HYP-4446/opus-S112 (theta-sum), HYP-4472/S19 (compactness),
HYP-4452/S17 (safe), residue_pinning_13, HYP-4396/opus-S107 (sum-product).

## S21 addendum: the equality-case core formalized, and the circularity honestly flagged

**Formalized (green, kernel-pure).** `LRCMinimalSumset.lean` now has both:
- `two_mul_card_sub_one_le`: `|S+S| ≥ 2|S|−1` (the bound);
- `sumset_eq_translates`: at equality `|S+S| = 2|S|−1`, `S+S = (min S + S) ∪ (S + max S)`
  EXACTLY — the structural core of "minimal doubling ⟹ AP" (step 2).

The full AP conclusion (translate structure ⟹ arithmetic progression) is the
classical remainder — a fiddly ordering induction, left for a dedicated pass.

**The circularity, flagged honestly.** The (U)-factoring
`safe=0 ⟹ [step1] min-doubling ⟹ [step2] AP ⟹ [step3] dilated AP` is
STRUCTURALLY clean but NOT an independent route to (U): step 1
(`safe=0 ⟹ min-doubling`) is TRUE (all safe=0 families are dilated APs, which
are min-doubling `|S+S|=23`, verified c=1..7) but its ONLY known proof routes
through (U) itself (`safe=0 ⟹ dilated AP ⟹ min-doubling`). An INDEPENDENT proof
of step 1 — "full theta-sum cancellation at 2/25 forces minimal doubling",
without assuming (U) — is the genuinely open additive-combinatorics target the
factoring isolates. So the value of the factoring is diagnostic (it names the
one open link cleanly) + the formalized classical steps (2 green), not a closure.
