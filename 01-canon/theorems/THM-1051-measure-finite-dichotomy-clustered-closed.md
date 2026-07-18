---
id: THM-1051
title: THE MEASURE/FINITE DICHOTOMY — the two-killer clustered case of "covering ⟹ M > 1/14" is CLOSED, and the non-covering lemma it replaces is REFUTED. (0) REFUTATION of my own THM-1041 named-next: no fixed modulus range can work, because k₁ = lcm(15,…,Q) with k₂ = 2k₁ is a LEGAL trapped-core-shape family (covering, compressed, gap family, verified legal for Q ≥ 30) every one of whose moduli in [15,Q] is killing — a killer divisible by q is unsafe at q for every multiplier. So "the killing sets never cover" is FALSE. (I) But defeating every modulus forces the killers to EXPLODE — the minimum max(k₁,k₂) defeating all of [15,Q] runs 271 (Q=18) → 1274 (Q=20) → >4000 (Q=24) — and huge killers are exactly what a MEASURE argument wants. (II) THE CORE-SAFE SET: S(P) = {t : ‖pt‖ ≥ 1/14 ∀p ∈ P} is an exact finite union of rational intervals; for the twelve cores P ⊆ {1,…,12} with |P| = 11 it has 12–20 components, total measure 0.051–0.127, largest component 0.0077–0.042. (III) THE MEASURE HORN: inside an interval of length L on which core and k₁ are safe, k₂'s unsafe set has measure ≤ L/7 + 2/(7k₂), so a good t survives whenever **k₂ > 1/(3L)**. Removing any single killer k₁ < 874 from S(P) exactly, the worst surviving component over ALL twelve cores and ALL k₁ is L = 0.00098184 (attained at k₁ = 873, uniformly across cores), giving a threshold of **339.5**. (IV) THE FINITE HORN: all covering families with BOTH killers < 874 — there are 41,986 — are certified by the small-modulus criterion at q ≤ 40, **41,986 / 41,986, zero uncertified**. (V) THE HORNS OVERLAP: 339.5 < 874, so every two-killer clustered covering family is certified by one horn or the other, with no gap and a factor-2.6 margin. The adversarial lcm-killers of (0), which defeat every modulus, are certified by the measure horn at the explicit t = 556667/7280000
status: PROVED FINITE-EXACT for two killers and any core P ⊆ {1,…,12} with |P| = 11. The uniform split has three cases: both killers below 874 use the 41,986-family exact finite horn; one below and one at least 874 uses the exact all-core/small-killer component scan, whose minimum is L=2/2037; both at least 874 use the explicit core-safe interval of length 1/1092 and an elementary two-comb measure bound. This proof does not use THM-1061's sampled scaling claim. The finite banks have one implementation and are not kernel-internalized
source: kind-pasteur-2026-07-18-S128 (cont.59; owner: work the non-covering lemma on the composite-killer extremals)
depends_on:
  - THM-1041         # the criterion, and the named-next this refutes and replaces
  - THM-1032         # the explicit-modulus construction, now subsumed for r = 2
related:
  - THM-933          # block gluing — this is its two-block case made explicit and elementary
  - THM-1007         # single-killer + lacunary chains
  - THM-1026         # opus: the 13/7 overshoot that blocks the pure counting route
script: 04-computation/noncover_adversarial_kps_S128c59.py, scale_dichotomy_kps_S128c59.py, perturbation_closure_kps_S128c59.py, mixed_killer_gap_kps_S128c59.py, full_safe_set_kps_S128c59.py (+ .out)
---

# THM-1051 — the measure/finite dichotomy

## (0) The non-covering lemma is false

THM-1041's named-next asked for a lemma saying the killing sets of distinct moduli cannot
union to everything. **It is false, and the counterexample is immediate.** A killer
divisible by q is unsafe at q for *every* multiplier a, so q is killing. Take

> k₁ = lcm(15, 16, …, Q),  k₂ = 2k₁.

Then every q ∈ [15,Q] divides k₁, so every modulus in the range is killing and the
small-modulus criterion has no witness at all. And the family is **legal**: with core
{2,…,12} it is covering (13 and 14 both divide k₁ once Q ≥ 30), compressed (k₂/k₁ = 2), a
gap family (k₁ ≫ 13·12), and has 13 distinct speeds. Verified legal at Q = 30, 40, 60,
killing at 16/16, 26/26, 46/46 moduli respectively.

Since Q is arbitrary, **no fixed modulus range can suffice.** Recording this against my own
proposal rather than letting it stand.

## (I) But the refutation carries its own cure

Defeating every modulus is expensive. The minimum max(k₁,k₂) over legal families defeating
all of [15,Q]:

| Q | 18 | 20 | 22 | 24 |
|---|---|---|---|---|
| min max(k₁,k₂) | 271 | 1274 | 1274 | > 4000 |

Killers that beat the modulus search are *large*, and large killers are precisely what a
measure argument wants. That is the dichotomy.

## (II) The core-safe set, exactly

S(P) = {t ∈ [0,1] : ‖pt‖ ≥ 1/14 for all p ∈ P} is a finite union of intervals with rational
endpoints (breakpoints j/p ± 1/(14p)). Computed exactly for all twelve cores:

| core drop | components | measure | largest component |
|---|---|---|---|
| 1 | 12 | 0.1055 | 0.04167 |
| 6 | 14 | 0.0514 | 0.00765 |
| 7 | 16 | 0.1274 | 0.02083 |
| 11 | 20 | 0.1060 | 0.01190 |

A cheaper explicit sub-interval also exists and is worth recording: since 13 is prime and
exceeds 12, ‖p/13‖ ≥ 1/13 for every p ≤ 12, so ‖pt‖ ≥ 1/13 − 12|t − 1/13| ≥ 1/14 on
**I = [1/13 − 1/2184, 1/13 + 1/2184]**, length 1/1092. Verified exact: the minimum of
min_p‖pt‖ over all twelve cores and the endpoints of I is exactly 1/14. But I alone is not
enough — k₁ = 156 = 12·13 makes 156t an integer at t = 1/13 and its unsafe set swallows all
of I. Using the full S(P) removes that obstruction, which is why (III) uses S(P).

## (III) The measure horn

> On an interval of length L where the core and k₁ are safe, the unsafe set of a further
> killer k₂ is a union of intervals of length 1/(7k₂) spaced 1/k₂, so its measure there is
> at most L/7 + 2/(7k₂). A safe t survives whenever L/7 + 2/(7k₂) < L, i.e.
> **k₂ > 1/(3L)**.

Removing each small killer k₁ ∈ (13·max P, 874) from S(P) exactly and taking the largest
surviving component, the worst case over all twelve cores and all such k₁ is

> **L = 0.00098184** (attained at k₁ = 873, the same value for every core),
> **threshold 1/(3L) = 339.5**.

## (IV) The finite horn

Every covering family with both killers below 874 is a finite list. Checked exhaustively
with the small-modulus criterion at q ≤ 40, via bitmask intersection:

> **41,986 families, 41,986 certified, 0 uncertified.**

## (V) The horns overlap — no gap

Write the ordered killers as `k₁<k₂`.  The uniform proof has three cases:

- if `k₂<874`, the finite horn certifies it (IV);
- if `k₁<874≤k₂`, the exact removal scan gives a surviving component of
  length at least `L=2/2037`, and `k₂≥874>1/(3L)=679/2`, so (III)
  certifies it;
- if `874≤k₁<k₂`, use the explicit core-safe interval `I` of length
  `ell=1/1092`.  Each killer removes at most `ell/7+2/(7kᵢ)`, hence together
  they remove at most

  `2ell/7 + 4/(7*874) < ell`.

  A point of `I` survives both combs.

These cases are exhaustive.  Direct spot-checks additionally confirm the
two measure regimes: 220/220 random families with both killers ≥ 874, and 294/294 *constructed*
covering mixed families (small k₁, large k₂ up to 5·10⁵). The adversarial lcm-killers of
(0) — the families designed to defeat every modulus — are certified by the measure horn at
the explicit witness **t = 556667/7280000**.

## Why this is the right shape

The two horns are the two regimes the problem actually has, and the refutation in (0) is
what reveals it: a family can only escape the arithmetic (modulus) argument by making its
killers divisible by everything in range, and that same act makes them enormous, which
hands the family to the analytic (measure) argument. The failure mode of one tool is the
hypothesis of the other. This is THM-933's block gluing specialised to two blocks and made
elementary — no discrepancy machinery, just an interval length and a union bound.

Also worth recording: M({2,…,12}) = **1/7** exactly, at t = 1/14 (the binding speed is 2,
giving 2/14). And M({2,…,12} ∪ {V, 2V}) = 1/7 for every V from 157 to 1600 — the
scale-separated families sit at exactly **twice** the 1/14 threshold, not near it.

## Why no scaling conjecture is needed

The mixed case uses the exact bounded removal scan only when `k₁<874`.
When `k₁≥874`, the proof changes to the simultaneous two-comb union bound on
the fixed interval `I`; it never extrapolates a sampled component ratio.
Thus THM-1061(III), withdrawn by MISTAKE-163, is not an input to this theorem.

## Named next
- **Run the finite horn for r = 3,…,6 killers.** The measure horn already generalises:
  r killers need k_min > 2r/(L(7−r)), giving thresholds 1638 (r=3), 2912 (r=4), 5460 (r=5),
  13104 (r=6), all finite. The finite horn for r = 3 is the next computation, and it is
  larger but still bounded. At r ≥ 7 the union bound dies (7 − r ≤ 0) and the core has ≤ 6
  speeds — a genuinely different regime.
- Lean: (III) is a measure computation over a `Finset` of rational intervals; the union
  bound is elementary. (IV) is a bounded decision procedure and is the natural
  `decide`-style target, though 41,986 families will want a reflected checker rather than
  kernel reduction.
