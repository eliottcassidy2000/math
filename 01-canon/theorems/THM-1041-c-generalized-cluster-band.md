---
id: THM-1041
title: THE c-GENERALIZED CLUSTER BAND (reach more than doubled, and SHARP) + THE SMALL-MODULUS WITNESS CRITERION — the wide-cluster case. (I) RESIDUE LAW: for q = v₁ + c the killer residues are exactly eᵢ = |δᵢ − c| (δᵢ = vᵢ − v₁), since vᵢ < q ⟹ q − vᵢ = c − δᵢ and vᵢ > q ⟹ vᵢ − q = δᵢ − c; THM-1032 is the single case c = M. (II) LETTING c FLOAT: if some integer c has dist(c, δᵢ) ∈ [μ, M] for every i, and M ≤ 12μ, then q = v₁ + c with a = ⌈q/(14μ)⌉ certifies min‖v·t‖ ≥ 1/14. For two killers such a c exists **iff D ≤ M − μ or 2μ ≤ D ≤ 2M**, and when M ≥ 3μ (true for every core ⊆ {1,…,12}) these merge to **D ≤ 2M** — the explicit choice is the CLUSTER MIDPOINT c = ⌈D/2⌉. (III) SHARP: at D = 2M+1 no c exists (c ≤ M forces |c − D| ≥ M+1), and the measured first failure is exactly D = 2M+1 for all 12 cores (25 at M=12, 23 at M=11). So the explicit reach goes from THM-1032's D ≤ M − μ = 10–11 to **D ≤ 2M = 22–24**, verified 16,874/16,874 exhaustive with ZERO failures. (IV) WIDENED CAP: admitting residues up to U with U < 13μ extends to D ≤ 2U; at μ = 2 the cap U = 22 covers D ≤ 36 fully. (V) THE μ=1 DEGENERACY DIAGNOSED: the band certificate fails on exactly 21 of 7880 covering wide clusters with μ = 1 and 0 of 715 with μ ≥ 2 — because e_min = 1 forces the admissible interval to length q/168 sitting inside (1,2), holding no integer, while the true witnesses (a = 4,5,7) work through WRAPAROUND, which the no-wrap band forbids. (VI) THE SMALL-MODULUS WITNESS CRITERION, which tolerates wraparound: if some q and a satisfy ‖v·a/q‖ ≥ ⌈q/14⌉/q for every speed, then M ≥ 1/14; searching only q ≤ 28 certifies 8589/8595 covering wide clusters and q ≤ 40 certifies **all 8595**, including all 943 families of the hard initial-segment core {1,…,11}
status: (I)-(IV) PROVED (residue arithmetic + the interval-existence argument of THM-1032; the c-existence characterization and its sharpness are elementary and are confirmed exhaustively 16,874/16,874 with the first failure at exactly 2M+1). (V) exact census. (VI) the IMPLICATION is proved and trivial; the EXISTENCE of a good small modulus is VERIFIED (8595/8595 at q ≤ 40, median 12 valid (q,a) pairs) but NOT proved, and the honest reason is recorded below: a union bound over 13 speeds costs 13/7 > 1, the same overshoot that blocks LRC(14) itself
source: kind-pasteur-2026-07-18-S128 (cont.58; owner: work the wide cluster case)
depends_on:
  - THM-1032         # the c = M case this generalizes
  - THM-1018 (II)    # the band lemma both consume
related:
  - THM-523          # t = 1/q; the q ≤ 14, threshold-1 case of (VI)
  - THM-1007         # single-killer + lacunary chains
  - THM-1026         # opus: the 13·(2λ) = 13/7 union-bound overshoot named in (VI)
script: 04-computation/c_generalized_band_kps_S128c58.py, small_modulus_avoid_kps_S128c58.py, wide_mu_split_kps_S128c58.py, initial_segment_core_kps_S128c58.py, wide_misses_kps_S128c58.py (+ .out)
---

# THM-1041 — the c-generalized cluster band, and the wide-cluster census

## (I) The residue law

Let V = P ∪ K, core P ⊆ [μ, M], killers K = {v₁ < … < v_r} all > 13M, offsets
δᵢ = vᵢ − v₁ (so δ₁ = 0, δ_r = D). For **q = v₁ + c** with 0 < c and D − c < q/2:

> **eᵢ = |δᵢ − c|.**

*Proof.* If δᵢ < c then vᵢ < q, so vᵢ mod q = vᵢ and q − vᵢ = c − δᵢ; since v₁ > 13M the
least-absolute residue is c − δᵢ. If δᵢ > c then vᵢ > q and vᵢ − q = δᵢ − c, which is the
least-absolute residue. ∎

THM-1032 is the single choice c = M, giving eᵢ = M − δᵢ and requiring δᵢ ≤ M − μ. Nothing
forces c = M.

## (II) Letting c float

> **If some integer c satisfies dist(c, δᵢ) ∈ [μ, M] for every i, and M ≤ 12μ, then
> q = v₁ + c and a = ⌈q/(14μ)⌉ give min_{v∈V} ‖v·t‖ ≥ 1/14.**

Every residue then lies in [μ, M], so e_min = μ and e_max = M exactly as in THM-1032, and
the integrality argument there applies verbatim (q ≥ 14M + 1 and M ≤ 12μ).

**When does such a c exist?** For r = 2 killers, Δ = {0, D}, the requirement is
c ∈ ([−M,−μ] ∪ [μ,M]) ∩ ([D−M, D−μ] ∪ [D+μ, D+M]), which is non-empty **iff**

  D ≤ M − μ  (take c = M — THM-1032's branch),  or  2μ ≤ D ≤ 2M.

On the second branch the explicit choice is the **cluster midpoint**

> **c = ⌈D/2⌉**,  giving  |c| = ⌈D/2⌉ ∈ [μ, M]  and  |c − D| = ⌊D/2⌋ ∈ [μ, M].

When M ≥ 3μ the two branches overlap (2μ ≤ M − μ) and merge into a single interval
**D ≤ 2M**. Every core P ⊆ {1,…,12} with |P| = 11 has (μ,M) ∈ {(1,12),(2,12),(1,11)}, so
M ≥ 3μ always holds here and the clean statement is:

> **Explicit reach: every spread D ≤ 2M.**

## (III) Sharpness

At D = 2M + 1 no c exists: c ∈ [μ,M] forces |c − D| = D − c ≥ M + 1 > M, and c ∈ [−M,−μ]
forces |c − D| = D + |c| > M. Measured first failure, per core:

| core drop | 1 | 2–11 | 12 |
|---|---|---|---|
| (μ, M) | (2,12) | (1,12) | (1,11) |
| 2M | 24 | 24 | 22 |
| first failing D | **25** | **25** | **23** |

Exhaustive over all 12 cores × D ≤ 2M × 59 values of v₁: **16,874 / 16,874 certified,
zero failures.** Against THM-1032's D ≤ M − μ this is a reach of 22–24 in place of 10–11.

## (IV) Widening the cap

Admitting residues up to U rather than M (still e_min = μ) needs U < 13μ and
q ≥ 14Uμ/(13μ − U); the reach becomes D ≤ 2U. At μ = 2 the cap U = 22 is affordable
(q ≥ 154, and q = v₁ + c > 156) and covers D ≤ 36 fully, 260/260. At μ = 1 there is no
room at all — 13μ = 13 ≤ M + 1 — which is the same degeneracy as (V).

## (V) The μ = 1 no-wrap degeneracy

Census of covering wide clusters (core an 11-subset of {1,…,12}, killers 13a and 14b both
> 13M, spread > M − μ), band-certified by modulus scan:

| μ | tested | band-certified | rate |
|---|---|---|---|
| 1 | 7880 | 7859 | 0.9973 |
| 2 | 715 | 715 | **1.0000** |

All 21 failures have μ = 1, and the cause is visible: at the working modulus the *ratio*
condition passes comfortably (e_max/e_min ≤ 9 ≤ 13), but the admissible interval
[q/14, 13q/(14 e_max)] evaluates to e.g. [1.071, 1.990], [1.143, 1.857], [1.214, 1.973] —
**always inside (1,2), never containing an integer.** With e_min = 1 the interval length is
q(13/e_max − 1)/14, which at e_max ≈ 8 and q ≈ 16 is ≈ 0.7. The genuine witnesses use
a = 4, 5, 7 — far outside the band — and work because e·a *wraps* past q and lands back in
the safe zone. The no-wrap hypothesis of THM-1018(II) is what fails, not the mathematics.

## (VI) The small-modulus witness criterion

Dropping the no-wrap requirement gives the natural generalization:

> **If some q and some a satisfy la(v·a mod q, q) ≥ ⌈q/14⌉ for every speed v, then
> t = a/q gives min‖v·t‖ ≥ ⌈q/14⌉/q ≥ 1/14.**

This contains both earlier tools: THM-523's t = 1/q is q ≤ 14 with threshold 1 and a = 1;
THM-1018(II)'s band is the sub-case where no product wraps. For q ≤ 28 the threshold is 2,
so the criterion reads simply "v·a mod q ∉ {0, 1, q−1} for all v".

Census results on the 8595 covering wide clusters:

| search range | certified |
|---|---|
| q ≤ 28 | 8589 / 8595 |
| q ≤ 40 | **8595 / 8595** |

The 6 residual-at-28 families all have core **{1,…,11}** — the initial segment, i.e. the AP
extremizer — and are witnessed at t = 3/37 or 3/38, giving M ≈ 0.081 against the threshold
0.0714. A full sweep of that stratum (943 covering families with core {1,…,11}) is
**943/943 certified**. First-successful-modulus distribution is spread across q = 15…28
with no modulus dominating, and the slack is comfortable: over 3364 families the number of
valid (q,a) pairs has min 0, p10 = 6, median 12, max 30.

## Honest status of the existence question

(VI)'s implication is trivial; its *existence* half is verified, not proved, and the
obstruction is structural rather than a missing lemma. A union bound over the 13 speeds
costs each speed a bad-fraction of about 2⌈q/14⌉/q ≈ 1/7, for a total of 13/7 = 1.857 > 1.
That is precisely the overshoot opus named in THM-1026 as the reason pairs cannot close 13
runners. So the small-modulus existence step is not a soft counting exercise: **it is
LRC(14) itself in miniature**, and any proof of it will need the same structure the main
problem needs. Recording this rather than presenting the census as a proof.

What *is* proved outright is (I)–(IV): the explicit construction, now reaching D ≤ 2M and
sharp there.

## The gap map for the clustered regime

| spread D | tool | status |
|---|---|---|
| D ≤ 2M (= 22–24) | q = v₁ + ⌈D/2⌉, explicit | **PROVED** (here) |
| 2M < D ≤ ~2·22 at μ=2 | widened cap U | **PROVED** (here, μ ≥ 2 only) |
| middle | small-modulus witness | verified 8595/8595, not proved |
| D > 12·v₁ | killers themselves lacunary | **PROVED** (THM-1007) |

The middle band is the whole of the remaining wide-cluster problem, and it is now sharply
delimited on both sides.

## Named next
- Close the middle band. The union bound cannot do it (13/7); what might: the killers
  contribute only 2 of the 13 residues and the core is a *fixed* subset of {1,…,12}, so the
  criterion depends on the family only through (core, two killer residues) — a finite
  pattern space per modulus. A covering-system argument over a coprime modulus set
  (15,16,17,19,23,…) is the shape to try, since CRT makes the bad conditions independent.
- Lean: (I) and (II) are short additions to `LRCClusterBand` — the residue law is a
  two-case `Nat` argument and the midpoint choice c = ⌈D/2⌉ discharges by `omega`.
