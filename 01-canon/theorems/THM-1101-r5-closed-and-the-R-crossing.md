---
id: THM-1101
title: r=5 bottom-window R crossing and finite horn below 235; uniform r=5 remains open after the sampled all-scale inference is withdrawn
status: PARTIAL / AUDITED BY MISTAKE-164 — the explicit bottom-window R>1 row and the bounded finite-horn computation are retained, but neither fixed bottom windows nor one translated decay ray proves an all-scale tail. The claimed uniform r=5 closure is withdrawn. THM-1097 independently proves uniform r=4 by a sharp discrepancy theorem and a guarded full-complement bank
source: kind-pasteur-2026-07-18-S128 (cont.63; owner: run the r=5 finite horn and settle r=4's R over all triples)
depends_on:
  - THM-1081    # the R-ladder and the r=4 closure this settles and extends
related:
  - THM-1051, THM-1061, THM-1071, THM-1097, MISTAKE-164
script: 04-computation/r4_R_alltriples_kps_S128c63.py, r5_R_crossing_kps_S128c63.py, r5_finite_horn_kps_S128c63.py, r5_finite_horn_fixed_kps_S128c63.py, r5_all_high_proof_split_gap_fraction_referee_codex_S66.py (+ .out)
---

# THM-1101 — bounded r=5 horn and a sampled R crossing

> **Audit correction (codex-S73; MISTAKE-164).**  The computations below are
> fixed-window data.  They do not prove that every omitted four-killer tuple
> has `R<1`, because the scripts contain no reduction to their windows and the
> decay table follows only one translated tuple.  The former uniform `r=5`
> conclusion is withdrawn.  Uniform `r=4` is now proved independently by
> THM-1097, whose analytic guards plus finite bank cover every tuple.

## (I) r=4's bounded R scan, superseded by THM-1097

cont.62 scanned only killers in [lo, lo+55). Widening to **all triples with k₃ ≤ lo+100
across all 220 cores** returns the same answer:

> **max R = 0.98453**, at core [1,2,3,5,6,7,8,9,11], killers (150,156,158), T = 155.6.

And the decay beyond the window, with k₁,k₂ pinned at the worst pair:

| k₃ | 158 | 208 | 308 | 558 | 1158 | 3158 | 9158 |
|---|---|---|---|---|---|---|---|
| R | **0.98453** | 0.54096 | 0.38889 | 0.38889 | 0.30736 | 0.23252 | 0.20453 |

R falls away monotonically, through the 7/18 = 0.38889 asymptotic plateau and below. The
scan maximum lies at the bottom.  The final sentence formerly inferred that
the same holds for all triples; that inference was not proved.

THM-1097 now supplies the honest all-scale statement for `r=4`, using the
sharper final-comb target `L>1/(7k3)` rather than extrapolating this `R` scan.

## (II) An explicit r=5 obstruction to the older component target

> **max R = 1.28495**, at core [1,2,4,5,7,9,11,12], killers (158,160,162,164), T = 210.7.

An independent exact-rational replay sharpens this row to

```text
L = 41/25920,
T = 1/(3L) = 8640/41,
R = T/164 = 2160/1681 > 1,
7*164*L = 11767/6480 > 1.
```

Thus it genuinely refutes the older `L>1/(3k4)` target, while still passing
the sharper `L>1/(7k4)` target that would suffice for the fifth killer.

The ladder predicted a crossing at r=5 and there it is. The measure horn genuinely fails,
and for the first time in this arc the finite horn is **mandatory**, not a check.

Inside the scanned 22-point bottom windows, **1011** quadruples have `R>=1`
and the largest fourth killer is 178.  Translating the displayed worst tuple
puts its own `R` below one.  This does not prove that the infinite failure
region is confined to that bank.

### An exact all-high row missed by both branches

The omission is witnessed far beyond the finite horn by

```text
P=(1,2,4,5,7,9,11,12),
(k1,k2,k3,k4)=(294,298,299,303),
k5=320.
```

Exact endpoint subtraction after the first four killers gives

```text
N=198,
mu=14258767904/152794649007,
L=431/415716,
T_mass=N/(6mu)=5042223417231/14258767904,
T_comp=1/(3L)=138572/431,
T=min(T_mass,T_comp)=138572/431=321.512...,
R=T/303=138572/130593>1.
```

All five killers exceed `235`, while `k5=320<T`; therefore neither the
below-235 finite horn nor the asserted `k5>T` measure tail applies.  This is
not a failure of LRC: the thirteen-speed family has the explicit lonely
witness `(q,a)=(22,7)`.  It is in the intended covering branch for
`q=2,...,14`: the core covers the ordinary small moduli, `294` supplies 14,
`299` supplies 13, and `320` supplies 8 and 10.

The gap is recurrent rather than an isolated mistranslation.  On the same
core, the toothpick template `(K,K+4,K+5,K+9)` has `R>=1` at 27 scanned
values through `K=363`.  For example `(238,242,243,247)` has
`T=278.512...`.  The one translated decay ray therefore cannot be a
monotonicity principle for independently moving combs.

## (III) The bound correction

My first r=5 run set the finite-horn bound at **max k₄ + 20 = 198**. That is wrong. The
fifth killer is certified by the measure horn only once it exceeds **T**, not once it
exceeds the largest removed killer — and max T over the failing region is **210.7**. So
KB = 198 left a real gap at k₅ ∈ [198, 211], and the run that reported 11,702,422
quintuples with zero failures did not actually close r=5.

Rerun at **KB = max T + 25 = 235**. Recording the error rather than quietly replacing the
number, because the first run looked exactly as convincing as the second.

## (IV) The bounded r=5 finite horn

> **263,708,305 quintuples** passing the covering-necessary condition, over all 495
> eight-speed cores, **ZERO uncertified.**

This computation has declared scope `k_i<235`.  The former closure sentence
combined it with the sampled `max T=210.7`; because that maximum is not
uniform, tuples outside the finite box remain open.  A separately structured
exact referee for the large bounded census is also desirable.

The raw below-235 bank has `17,484,487,020` ordered five-subsets before the
covering prune.  Merely extending the blind cutoff past the exact gap would
already require `B>=322`, whose corresponding raw bank has
`583,797,584,496` rows, and still would give no theorem for the tail.

The all-high row is independently replayable from exact fractions at

```text
04-computation/r5_all_high_proof_split_gap_fraction_referee_codex_S66.py
05-knowledge/results/r5_all_high_proof_split_gap_fraction_referee_codex_S66.out
source SHA-256  a71ca9fee18ee5b4e999697f71ad48ee059ce582ec89b90a81c9ca0e8f08bab7
output SHA-256  0de9c790841527872d8faa5484af83fa034af7425a6df2a5272540c41be5ee7f
```

Normal and optimized Python runs byte-match the frozen output.  The replay
retains exact endpoint owners and also verifies the post-fifth `(22,7)`
witness component, so the proof-gap diagnosis does not depend on decimals.

## (V) The measured ladder

> **0.51852 (r=2) → 0.73375 (r=3) → 0.98453 (r=4) → 1.28495 (r=5)**

The displayed values are useful bottom-window telemetry, not a complete
all-scale ladder.  THM-1097 does reveal a rigorous phase transition for the
coarse sharp-discrepancy mass/component mechanism: its asymptotic inequality
works through three removals and fails to force the fourth-removal case.
Failure of that sufficient inequality is not failure of the conjecture.

More explicitly, if a core interval has length `ell`, four exact removals
combined only with the one-comb mass and incidence bounds give the sufficient
condition

```text
14 ell k4 - 7 ell(k1+k2+k3)
 - 6(k4/k1+k4/k2+k4/k3) - 45 > 0.                    (1)
```

Along near-equal scales this tends like `-7 ell K-63`, so (1) gets worse as
the toothpicks dilate.  THM-1097's coarse mass/component guard cannot simply
be iterated once more.  A uniform proof has to recover overlap or the
endpoint-owner recurrence that the separate mass and cut counts discard.

## Named next
- r=5: prove a uniform four-comb endpoint/self-similarity theorem or a valid
  reduction of every large tuple to a finite guard bank.  Additional rays do
  not discharge this quantifier.
- Independently replay the below-235 horn without floating-point pruning.
- r=6 should not be promoted until the r=5 all-scale bridge is real.
- r ≥ 7: the union bound underlying the measure horn dies outright (7 − r ≤ 0) and the core
  has ≤ 6 speeds — a different regime, not an extension of this one.
- Keep the finite horns as bounded verification and state their boxes.
