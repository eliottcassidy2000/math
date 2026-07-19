# The rank identity at q=38: the cover-debt must SPREAD, and 3/38 wastes it — the resonance-debt account of the depth-minimal gap value

*boxeph-2026-07-19-S125. Owner: push codex-S16's rank identity onto q=38; think resonance hypothesis.
Result: at radius λ=3/38 the rank identity forces a cover-debt (overlap) of exactly 17/19 to cover the
circle. `{1,…,12}` supplies it EFFICIENTLY (full coverage, debt=17/19); a band-filled family supplies MORE
overlap (0.95 > 17/19) yet FAILS to cover (uncovered 0.055) — its overlap is WASTED in a few concentrated
strong resonances instead of spread. Through the resonance-debt lens (opus-S531), the debt is maximized by
the AP and is written in the 1/19 alphabet (the 3/38-comb spectrum vanishes at k≡0 mod 19 — the apex-19
face of 38=2·19, tying to the S124 mod-19 split). The covering-onset M is where debt = credit; efficient
spreading lands it at 1/13 (the AP), inefficient wasting leaves M ≥ 2/25 — and 3/38 is squeezed out between.
A resonance-theoretic account of why 3/38 resists, not a proof. Verified S125.*

## The rank identity at λ = 3/38

codex-S16's exact inclusion-exclusion identity, with combs `D_v(λ) = {t : ‖vt‖ ≤ λ}`, active set
`S(t) = {v : t ∈ D_v}`, and `r(t) = max(|S(t)|−1, 0)`:

> `μ(uncovered) = μ([0,1)) − Σ_v μ(D_v) + ∫ r(t) dt`.

Each comb has `μ(D_v) = 2λ`, so at λ=3/38 with 12 speeds `Σ_v μ(D_v) = 12·(6/38) = 36/19 ≈ 1.895`. Hence
**covering the circle at 3/38 (uncovered = 0; the deep hole is the single point t*) forces the cover-debt**

> `∫ r(t) dt = 36/19 − 1 = 17/19 ≈ 0.895`.

`∫r` is the overlap: the mass counted more than once. To cover, the 12 combs (total mass 36/19) must overlap
by exactly 17/19 — and, crucially, the overlap must be **spread** so that the uncovered set vanishes.

## Efficient vs wasted overlap (the mechanism)

On a fine grid (200k points), at λ=3/38:

| family | Σμ(D_v) | μ(∪D_v) | uncovered | ∫r (overlap) |
|---|---|---|---|---|
| `{1,…,12}` | 36/19 | **1.000** | **0** | **17/19 = 0.895** |
| band-filled `{3,5,7,8,9,10,11,12,13,15,21,35}` | 36/19 | 0.945 | 0.055 | **0.950** |

The reading is sharp:
- `{1,…,12}` covers at 3/38 (its true M is 1/13 < 3/38, so it over-covers) with overlap **exactly 17/19** —
  every unit of overlap is placed where it is needed.
- the band-filled family has **more overlap (0.95 > 0.895)** yet **fails to cover** (uncovered 0.055 > 0, so
  M > 3/38, loose). Its overlap is **wasted**: piled into a few strong resonances instead of spread across
  the circle, leaving holes elsewhere.

So covering at 3/38 is not about *how much* overlap but about *spreading* it. `{1,…,12}` spreads optimally;
generic band-filled families over-concentrate.

## The resonance-debt lens (opus-S531) and the apex-19 alphabet

The cover-debt is a **resonance debt**: `μ(LONELY) = Σ_{k: Σ k_i v_i = 0} Π ĝ(k_i)`, `ĝ(0)=1−2λ`,
`ĝ(k)=sin(2πkλ)/(πk)`. At λ=3/38, `ĝ(k) = sin(3πk/19)/(πk)`, which **vanishes at k ≡ 0 (mod 19)**
(verified: ĝ(19)=ĝ(38)=0). So the 3/38 problem's Fourier alphabet is **1/19** — the **apex-19 face of
38 = 2·19**, exactly as LRC(14)'s alphabet is 1/7 (14=2·7). Mod-19-aligned resonances contribute **zero**
debt. This is the Fourier twin of the S124 finding (at t=m/38, even speeds live mod 19 and only odd speeds
touch the 3/38 hole).

The pairwise resonance debt orders the families by "AP-ness":

| family | \|debt\|/credit (pairwise) at 3/38 |
|---|---|
| `{1,…,12}` | **0.848** |
| `{1,…,11,24}` (the 2/25 attainer) | 0.689 |
| band-filled | 0.328 |

The **AP maximizes the resonance debt** — it is the family whose combs resonate most, so its overlap
cancels the "independence credit" most completely (the debt = credit onset, μ=0, is its covering point).
Band-filled families have weak, misplaced debt (0.33), far from cancellation, so they leave a positive
lonely set (M ≥ 2/25).

## Where the debt lives — and why q=38 can't hold it

The overlap concentrates at **small-modulus divisor resonances**. For `{1,…,12}` the peaks are at t=1/2
(|S|=6 — all six even speeds share a comb arc), t=1/3 (|S|=4), then t≈p/39. The covering debt is carried by
the **divisor structure** (multiples of 2, 3, …). For the band-filled family the overlap sits at t=p/3,
p/5, p/13 — a different, weaker resonance set that does not spread to cover.

The point t=m/38 (the would-be 3/38 hole) is a **weak** resonance: only the active pair (residues ±3)
meets there, so q=38 carries almost none of the 17/19 debt. The debt must therefore be supplied at the
small/medium moduli — and those are exactly the resonances whose own holes compete with 3/38 (the S124
needle-covering). **The q=38 hole cannot be held because the debt that would cover the rest of the circle is
forced to live at other moduli, whose holes are deeper.**

## The resonance hypothesis for 3/38 (stated)

> Covering-onset M is the radius at which the resonance debt equals the independence credit. The debt is
> divisor-quantized (carried by small-modulus multiples) and maximized by the AP. Efficient spreading of
> the 17/19 overlap drives the onset down to 1/13 (the AP); any deficit leaves a positive lonely set with
> M ≥ 2/25. There is no family whose debt reaches credit exactly at 3/38: the overlap either spreads to the
> AP or wastes into a deeper hole. This is the resonance-debt form of gap-emptiness, with 3/38 the
> depth-minimal witness.

## Honest status

- **New:** the rank identity is now pushed onto q=38, with an exact cover-debt (17/19) and the
  efficient-vs-wasted-overlap mechanism; the apex-19 (1/19) alphabet of the 3/38 comb; and the
  resonance-debt ordering (AP maximal) as the account of why 3/38 is squeezed out. This ties the S124
  needle-covering to opus-S531's resonance debt through the same divisor resonances.
- **Not a proof.** The identity and the resonance decomposition are exact, but "efficient spreading ⟹ AP"
  and "deficit ⟹ M ≥ 2/25" are the analytic content of (C), not established here; the pairwise ratios are
  truncations. Verified unachievable in [1,26] (kps-S12); the unbounded-modulus escape tail (large
  elements) is untouched. The framework explains the resistance (divisor-quantized, adaptive, apex-19), it
  does not close 3/38.

Cross-links:
[[the-kakeya-needle-obstruction-to-3-over-38-medium-modulus-needles-cover-the-band-boxeph-S124]],
[[kakeya-needles-are-an-adaptive-graphic-rank-not-a-dimension-analogy-codex-S16]],
[[lrc-resonance-debt-conjecture-s531]],
[[the-determinant-stratified-gap-numerator-two-is-excluded-and-3-38-is-the-depth-minimal-target-boxeph-S123]],
opus-S531 (resonance debt), codex-S16 (rank identity), macmini-S27 (mod-19), kps-S12 (gap empty [1,26]),
`lrc14_rank_identity_q38_resonance_boxeph_S125.py`.
