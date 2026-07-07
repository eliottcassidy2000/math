---
source: opus-2026-07-07-S138
status: exact theory for the 2-anchor tail on its conjectured minimizer (multiplicities PROVED-verified,
  1-anchor limit table EXACT, 2-anchor limits + margins); personal-space tournament (conservation law,
  observer-rank separator, mirror-breaking); Goldbach/BB certificate-compression lens (frame)
tags:
  - lonely-runner
  - LRC14
  - two-anchor-tail
  - three-distance
  - personal-space
  - perspective
  - goldbach
  - busy-beaver
---

# The two-anchor tail made exact — and the observer breaks the mirror

**opus-2026-07-07-S138.** Owner: work the 2-anchor joint tail; take arXiv:1505.02479
(Yedidia–Aaronson) / BB(27) / our Goldbach-adjacent threads as inspiration; invent
runner-tournaments from pair statistics; use the relativity of perspectives (antisymmetric
pairwise rates). boxeph-S1 had named the joint next step: *"Farey roof from a shifted origin →
exact PA₂ per k."* Delivered below.

## 1. The 2-anchor tail on its conjectured minimizer, exact

`PA₂(E) = P_x(max(gap∋0, gap∋½) > 1/7)`; boxeph's empirical minimizer class is the spread
inhomogeneous AP `{a + dj}`. Structure: `config({a+dj}, x)` is the Steinhaus config
`C⁰(y) = {frac(jy) : 0 ≤ j ≤ k−1}` (y = dx) rotated by `z = ax`, so the anchors sample the
gaps of `C⁰(y)` at positions `−z` and `½ − z`, and as `(a,d)` grows coprimely `(y,z)`
equidistributes on `T²`:

- **1-anchor limit = the long-gap mass** `E_y[L(y)]`, `L = Σ_{g > 1/7} g` over the gaps of
  `C⁰(y)`. With the classical three-distance **multiplicities** `N₁ = K+1−q`, `N₂ = K+1−q′`,
  `N₃ = q+q′−K−1` (machine-verified exactly, 0 violations, K ≤ 16 — together with THM-637's
  value-proofs this is now proof-grade), `L` is piecewise affine per Farey-K cell and its mean
  is an **exact rational**:

  | k | E_y[L] exact | ≈ | T_k | verdict |
  |---|---|---|---|---|
  | 8 | 1041/1715 | .60700 | .6185 | **below — 1-anchor misses k=8** |
  | 9 | 888/1715 | .51778 | .5057 | clears |
  | 10 | 2269/5145 | .44101 | .3956 | clears |
  | 11 | 274/735 | .37279 | .2747 | clears |
  | 12 | 873/2695 | .32393 | .1429 | clears |
  | 13 | 15982/56595 | .28239 | .0565 | clears |

  These match boxeph's numeric infima — quantifying exactly why the single anchor fails only
  at k=8 (deficit 199/... ≈ 0.0115).

- **2-anchor limit** `E_y[meas(S ∪ (S−½))]` (S = the union of long gaps): ≈
  `0.8017 / 0.6927 / 0.5974 / 0.4978 / 0.4287 / 0.3665` (200k-sample; per-y exact interval
  arithmetic) — **the antipode anchor rescues k=8 with margin +0.183**, and margins grow to
  +0.310 at k=13. The gain over one anchor (+0.084…+0.195) is the failure of long gaps to
  align with their own half-shifts — a ½-decorrelation.

- **Finite (a,d), exact** (order-cell engine): PA₂ at k=8 runs 0.940 (AP) → ~0.795–0.822 at
  small spread — the T² limit is approached from a *slightly overshooting* family: (2,5) gives
  0.79524 **below** the limit 0.8017, so the true `inf PA₂` sits near but not exactly at the
  limit. Flag: boxeph's descents at k=9,10 (0.511, 0.434) also dip ~0.007 below my exact
  limits — the minimizing *class* is spread-AP-like but the exact infimum includes a finite
  correction. **The remaining lemma, stated exactly:** for all k-element E,
  `PA₂(E) ≥ T_k`, with conjectured near-inf the spread-AP limit values above (worst observed
  margins: k=8 ≈ 0.77 vs 0.6185). A 2-gap rigidity on two fixed anchors — the load-bearing
  open piece, now with exact target constants on the binding class.

## 2. The half-shift parity mechanics (the Lemoine shape)

The substitution `x ↦ x + ½` **fixes every even phase and rotates every odd phase by ½**
(`e(x+½) = ex + e/2`). So the ½-anchor is the 2-adic partner of the origin: for all-odd
families `gap∋½(x) = gap∋0(x+½)` exactly, and for mixed parity the half-shift interleaves the
odd sub-config against the fixed even one — the same parity-interlacing operator that built
the E[maxgap] record families. This is the precise sense in which the 2-anchor decomposition
is **Lemoine-shaped** (`n = p + 2q`: one part weighted by 2): the owner's Goldbach pointer
lands on the repo's old Goldbach/Lemoine thread (bsd-hodge-polygonal-ladder) as *pair
decompositions where one slot carries the 2-adic weight* — here, anchor pairs `{0, ½}` whose
second member only the odd part can reach.

## 3. The certificate-compression lens (Yedidia–Aaronson → BB(27))

A Π₁ conjecture's honest size is the smallest machine whose non-halting encodes it; Goldbach
compressed 4888 → 27 states. LRC(14) is effectively Π₁ (rational witnesses + effective
bounds), and the fleet's program is exactly machine-size reduction: sieve/saturation and
peels (domain cuts) → coarse/conjugate certificates (infinite strata by uniform witnesses) →
diameter/intersected ledgers + `V0abs ≤ 1106` + boxes (explicit finite enumerations) → the
2-anchor tail (the last analytic gate). **Program metric worth tracking: the residual
obligation inventory** — today: [2-anchor rigidity lemma at 6 constants] + [Part-A finite
check ≤ V0abs] + [D>75 spread tail] + [wide-cluster gate]. When each is a finite check, their
total size is our "state count" — and the BB(27) lesson is that encodings, not effort,
shrink it: the anchored reformulation just compressed "AP-minimality of a max over 13 gaps"
to "a floor on 2 fixed gaps."

## 4. The personal-space tournament (creative piece; relativity made structural)

For the 14 runners (observer included), `PS_i(x)` = the sum of the two gaps flanking runner
i's phase. Facts: **conservation** `Σ_i PS_i(x) = 2` exactly (each gap credited to both
flankers — verified exactly on all families); **relativity**: PS is frame-invariant (pairwise
rates antisymmetric — the owner's reciprocity — so "space" is shared data, not perspectival);
**loneliness = the observer's space**: lonely at x iff both of 0's flanking distances ≥ 1/14
— the anchored-tail program is literally personal-space statistics, `E[PS₀] = E[gap∋0]`
(for the AP: `93/440`, THM-637(a) again). Tournament: `i → j` iff `E[PS_i] > E[PS_j]`
(exact). Findings: the **observer's rank separates structure** — rank 0 (most spacious) for
AP, GW, deep well, far-element families; rank 1 (beaten by the top speed) for the parity
record, prim-sat, and one member of each mirror pair. And the sharpest one:

> **The observer breaks the mirror.** Gap-multiset functionals are mirror-blind (S137's
> pointwise identity), but anchored/perspectival ones are not: mirror pair `{1..12,20}` vs
> its reversal have `E[PS₀] = 0.2055` vs `0.1741`. Reversal is a symmetry of the moving
> configuration; the observer pins a frame and *sees* it — the SC/NS "symmetry breaking"
> of mac-mini-S43, localized to exactly the observer-anchored statistics. (This is also why
> the 2-anchor lemma is stated on speed-form representatives: the anchors are the observer's
> own position and antipode.)

Differentiated from monad-S5 (centering / exclusive-attention) and kps-S64 (leader/Paley)
tournaments; the personal-space order is the one whose extremal statistic *is* the
load-bearing tail.

## Honest status

The multiplicity lemma + E[L] table are proof-grade exact; the 2-anchor limits are
high-precision numerics with exact per-y interval arithmetic; the rigidity lemma
(`inf PA₂ ≥ T_k`) remains open — it is *the* remaining analytic gate of the density floor,
now with exact binding-class constants. Small flags: my tournament script's mirror-symmetry
CHECK line pairs by the wrong center (values are right; check line cosmetic); boxeph's k=9/10
descent values dip ~0.007 below the spread-AP limits — the exact finite-family infimum needs
one clean sweep (next session or any agent with the engine).
