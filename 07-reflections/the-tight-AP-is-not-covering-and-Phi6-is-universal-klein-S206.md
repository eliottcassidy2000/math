---
source: klein-2026-07-09-S206
status: THREE findings, one of them potentially load-bearing for the whole covering programme.
  (1) The "13-comb Eisenstein resonance" (t*=14/183, 183=Phi6(14)) is NOT special to n=14: `n^3 ≡ -1`,
  `n(n-1) ≡ -1 (mod Phi6(n))` and `CF(n/Phi6(n)) = [0; n-1, n]` are IDENTITIES for every n (from
  `x^3+1=(x+1)Phi6(x)`). What IS n-dependent is the covering cushion `n^2/Phi6(n) → 1` — 25.0% at n=4,
  7.10% at n=14 — which is why LRC gets monotonically harder, with no n=14-specific coincidence.
  (2) VERIFIED EXHAUSTIVELY (n=4..7) that the TIGHT locus (`M=1/n` exactly) is entirely NON-COVERING:
  `#tight ∩ covering = 0`; tight sets are `{AP, GW}`. The covering-min is `2/(2n-1)` (stable under cap
  increase), cushion `2n/(2n-1)`.
  (3) THE CONSEQUENCE: `{1,…,13}` is NOT covering. Its co-offset cluster `{0,…,12}` at `Vmax=13` — the
  cluster with NO good period, which killed FIVE fleet routes (MISTAKE-127 arc-count, -128 c<D3, -129
  smooth-mean, kps kissing, -130 widest-arc) — never belongs to the covering case at all: THM-523
  dispatches it with the explicit witness `τ=1/14`. VERIFIED: EVERY primitive covering 13-set has a
  STRICT good period (966 exhaustive at speeds ≤18, 400 random ≤60, adversarial hill-climb at caps
  20/30/50/100): min margin `7·maxgap/Vmax = 1.2353`, never ≤ 1.
tags: [lrc14, covering-set, thm-523, phi6, tight-locus, good-period, prove-or-disprove, cross-n]
---

# The tight AP is not covering — and Φ₆ is universal

**klein-2026-07-09-S206.** Owner: work the bounded corner; investigate the uniqueness of the 13-comb
Eisenstein resonance and its relation to LRC(14)'s hardness across runner counts; see how LRC(14) must be
proved or disproved.

## 1. The Φ₆ resonance is universal, and the *cushion* is what depends on n

`Φ₆(x) = x²−x+1` (the 6th cyclotomic polynomial); `Φ₆(14) = 183`. The "13-comb at `t* = 14/183`" rests on
`13·14 ≡ −1 (mod 183)` and `14³ ≡ −1 (mod 183)` (so `14` is a primitive 6th root — the Eisenstein flavour).
Both are **identities in n**, immediate from `x³+1 = (x+1)(x²−x+1)`:

> `n³ ≡ −1` and `n(n−1) ≡ −1 (mod Φ₆(n))`, and `contfrac(n/Φ₆(n)) = [0; n−1, n]`  — for **every** n.

(Verified n=3..20.) The comb spacing `n−1` is simply the first continued-fraction denominator of `n/Φ₆(n)`.
So the `14/183` resonance is the `n=14` member of a universal family; **nothing about it singles out 14.**

What *does* depend on n is the **covering cushion**. The proved/conjectured covering floor is `n/Φ₆(n)`
(THM-523/HYP-3778), against the LRC target `1/n`, so the cushion ratio is

> `ρ(n) = n²/Φ₆(n) = n²/(n²−n+1) → 1`:  `1.2308` (n=4), `1.1395` (n=7), **`1.0710` (n=14)**, → 1.

**That monotone shrink — not any n=14 coincidence — is the hardness.** LRC is proved for n ≤ 13; n = 14 is
simply the first n where a 7.1% cushion defeats current methods.

## 2. The tight locus is entirely non-covering (exhaustive)

THM-523: a 13-set omitting a multiple of some `q ∈ {2,…,14}` has the explicit lonely witness `τ = 1/q`, so
`M ≥ 1/q ≥ 1/14`. Hence LRC(14) reduces to **covering** sets (a multiple of *every* `q ∈ {2,…,14}`).
Exhaustively, for n = 4,5,6,7 (primitive, speeds ≤ 20–26), computing `M(S)` **exactly** (the local maxima of
the piecewise-linear `min_i‖v_i t‖` sit at `t = p/(v_i+v_j)`):

| n | #tight | #tight ∩ covering | tight sets | covering-min | cushion |
|---|---|---|---|---|---|
| 4 | 1 | **0** | {1,2,3} | 2/7 | 8/7 |
| 5 | 2 | **0** | {1,2,3,4}, {1,3,4,7} | 2/9 | 10/9 |
| 6 | 2 | **0** | {1,2,3,4,5}, {1,3,4,5,9} | 2/11 | 12/11 |
| 7 | 1 | **0** | {1,2,3,4,5,6} | 2/13 | 14/13 |

The tight locus is exactly `{AP, GW}` (Goddyn–Wong), confirming THM-612. The covering-min is `2/(2n−1)`,
**stable** under cap increase (n=4: cap 20→70; n=5: 22→55). At n=14 the covering-min descends with the speed
cap (`2/23` at ≤17, `1/12` at 18–20, the repo's `7/89` at larger speeds) but never approaches `1/14`.

## 3. The consequence: the counterexample that killed five routes is out of scope

`{1,…,13}` contains no multiple of 14 ⟹ **it is NOT covering**, and THM-523 dispatches it with `τ = 1/14`
(attaining equality `M = 1/14`). But its co-offset cluster is `E = {13−v} = {0,…,12}` at `Vmax = 13` — *the*
cluster klein-S201 showed has **no good period** (`margin = 7·maxgap/Vmax = 0.538 < 1`), and which broke, in
order: the arc-count pigeonhole (MISTAKE-127), `c < D3` (MISTAKE-128), opus's smooth grid-mean (MISTAKE-129),
kps's kissing certificate, and mac-mini's widest-arc pigeonhole (MISTAKE-130).

**That cluster never belongs to the covering case.** Restricting to clusters arising from *covering* velocity
sets, a strict good period always exists, comfortably:

> **966 exhaustive primitive covering 13-sets (speeds ≤ 18): 0 without a strict good period, min margin
> `1.2353`.  400 random (≤60): 0 failures, min `1.9091`.  Adversarial hill-climb minimising the margin at
> caps 20/30/50/100: minima `1.47 / 1.40 / 1.64 / 1.99` — never ≤ 1.**

(The worst case sits at small speeds — `(1,2,3,4,7,8,9,10,11,12,13,14,17)`, margin `1.2353` — i.e. inside the
exhaustive range, so this is not a random-sampling artefact of the kind MISTAKE-126/127 warned about.)

**Claim to reconcile with canon (not asserted over it):** the covering-case cluster analysis in the repo
quantifies over *primitive clusters `E` with `spread ≤ Vmax`*. The covering-derived family is strictly
smaller — `E` must satisfy, relative to `Vmax`, `∀q∈{2..14} ∃e∈E : q | (Vmax−e)`. Dropping that constraint
admits the tight-AP pathology. If THM-527's reduction does not already carry the covering constraint, then
LEM-013's adversarial band, the `7`-structured obstruction, and the whole "existence is a MAX not a mean"
saga were being fought on clusters LRC(14) never needs. **This should be checked before more effort goes into
the hard clusters.** (Filed as a hypothesis, not a canon override.)

## 4. So: how must LRC(14) be proved — or disproved?

> **LRC(14) ⟺ every primitive covering 13-set has `M(S) ≥ 1/14`.**

- **The equality locus is non-covering.** `M = 1/14` is attained only by `{AP, GW}`, all non-covering. They are
  certified by an *explicit witness* `τ = 1/q` — elementary, exact, with equality. **No margin, averaging, or
  discrepancy argument can ever certify them**, because there is no margin to find. Every fleet route that
  died did so by being pointed at this locus.
- **The covering branch has a strict cushion**, so margin arguments belong *there*: covering-min `≥ 7/89`
  empirically (10% above `1/14`), conjectured floor `n/Φ₆(n) = 14/183` (7.10% above; HYP-3778), gap
  `7/89 − 14/183 = 35/16287`. Uniform looseness remains HYP-2566 — **that is the real open problem.**
- **To disprove**, one must exhibit a covering 13-set with `M < 1/14`. The covering-min descends with the
  speed cap but stays far above; the conjectured floor `n/Φ₆(n)` is exactly the Φ₆ resonance of §1.

So the proof shape is forced: **[explicit `τ=1/q` on the non-covering branch, equality allowed] + [strict
margin on the covering branch, cushion `n²/Φ₆(n)`]**. And the shrinking cushion `→ 1` says the same method
cannot survive large n — LRC in general needs something else.

## 5. The bounded corner (honest)

klein-S205's drift-absorbed embed needs `Vmax > 1.41·spread`, i.e. every speed `> 0.29·Vmax`. Covering sets
contain small speeds, so `spread ≈ Vmax` and **the embed does not fire on them**. The embed treats all 13
runners as one cluster; the general case needs the `P ∪ L` split with a slow-runner part `R` (as
`scale_separation_phase` has). That extension — gap-centred fast phase **plus** a 1-Lipschitz slow-safe `R` —
is the concrete next Lean step for `hembed`, and the corner is the covering case itself, where good periods
are comfortable (margin ≥ 1.2353) and only the *realization* is missing.

Files: `lrc_phi6_resonance_across_n_klein_S206`, `lrc_covering_min_across_n_klein_S206`,
`lrc14_covering_min_descent_klein_S206`, `lrc14_covering_sets_have_good_periods_klein_S206`,
`lrc14_covering_margin_adversarial_klein_S206` (+outs). Touches THM-523, THM-612, HYP-2566, HYP-3778,
MISTAKE-127/128/129/130, klein-S201/S205.
