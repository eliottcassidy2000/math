---
id: HYP-2770
title: The affine CJJ route must retain gap-word or relation-code marginals after full-residue localization
status: OPEN proof target; exact k=8 scout
source: codex-2026-06-21-S74
depends_on:
  - HYP-2749
  - HYP-2754
  - HYP-2747
  - HYP-2744
  - HYP-2738
  - HYP-2726
  - THM-531
  - THM-534
related:
  - HYP-2758
  - HYP-2762
  - HYP-2756
  - HYP-2751
  - HYP-2746
  - HYP-2727
  - HYP-2723
  - THM-126
  - THM-134
  - OPEN-Q-108
external:
  - https://arxiv.org/abs/2211.01248
  - https://arxiv.org/abs/2112.09221
---

# HYP-2770: Affine CJJ Needs Gap-Word / Relation Data

## Claim

The CJJ linear-code hierarchy is still the right template for the LRC(14)
consec-max proof, but the affine replacement for "linear code" cannot be a
residue-only quotient.  After HYP-2749 reduces the problem to the
full-`Z/7`-residue stratum, `Aff(1,F_7)` is too transitive: residue-count,
affine-pair, affine-triple, and affine-quad profiles are all identical for the
entire bounded k=8 full-residue stratum tested below.

The useful hierarchy must therefore factor translation+dilation while retaining
one of the following non-residue channels:

```text
integer gap / generated AP word
  or
individual relation-code co-missing marginals P(A subset M).
```

This refines HYP-2754 rather than refuting it.  The affine group is the right
symmetry to factor; residue-only affine moments are just not the final
certificate.

## Exact Scout

Script:

```text
04-computation/lrc14_affine_full_residue_obligation_codex_s74.py
```

Stored output:

```text
05-knowledge/results/lrc14_affine_full_residue_obligation_codex_s74.out
```

External CJJ reading used:

- `arXiv:2112.09221`, *A Complete Linear Programming Hierarchy for Linear
  Codes*, introduces the Delsarte-generalizing LP hierarchy and its complete
  linear-code target.
- `arXiv:2211.01248`, *Exact Completeness of LP Hierarchies for Linear Codes*,
  proves exact recovery at low enough level and clarifies the role of
  pseudoprobabilities / linear-code structure.

Affine orbit facts over `F_7`:

```text
r=0: 1 orbit
r=1: 1 orbit
r=2: 1 orbit
r=3: 2 orbits, reps (0,1,3), (0,1,2)
r=4: 2 orbits, reps (0,1,2,4), (0,1,2,3)
r=5: 1 orbit
r=6: 1 orbit
r=7: 1 orbit
```

So level `<=2` residue-only affine data must collapse.

For the k=8 full-residue stratum inside `[0,15]`:

```text
total shapes = 528
max measS7 = 481/1470
argmax shapes = (0,1,2,3,4,5,6,7) and (0,2,4,6,8,10,12,14)

quotient class counts:
  residue_count          1 class, AP-class size 528
  affine_pair            1 class, AP-class size 528
  affine_triple_profile  1 class, AP-class size 528
  affine_quad_profile    1 class, AP-class size 528
  integer_gap_word     528 classes, AP-class size 1
```

The second argmax is the dilation `2*AP`, as predicted by THM-531.  Thus the
right statement is not "consec is unique" but "the AP dilation orbit is unique"
inside the full-residue stratum.  However, all 528 shapes are indistinguishable
by residue-only affine profiles through degree 4, so the proof cannot stop at
that quotient.

Incoming KPS `affine_symmetric_lp_kpswf6.py` is consistent with this: AP
maximizes the distinct-residue-pair count and the signed `L_y` value for bounded
k=8,9,10 scans, but `256/432/400` shapes share AP's affine pair signature in
the k=8/9/10 bounded banks.  That signature localizes the safe full-residue
regime; it does not pin the AP orbit.

Late-pulled S19/LAYER-2 work strengthens the interpretation rather than
changing it.  The new multiplicative-symmetry proof records
`W_a(cE)=W_{ca}(E)`, so the `Z_7^*` dilation action is a real symmetry spine
for the consec/AP layer.  After the S20 correction (HYP-2762), do not read this
as saying Paley/QR is the large-p tournament H-extremizer; the additive
interval/AP is the H-side driver, while `Z_7^*` remains a useful symmetry
feature.  In HYP-2770 language, this supports P1: quotient by translation and
dilation first.  The simultaneous
Freiman/anchored-profile output says the translation-invariant leg profile is
not sufficient, but consec is maximal inside the minimal anchored leg-sum
layer; the extra second-profile/gap word is exactly the non-residue channel
P2 retains.

The later mac-mini Route-3 Huffer-Shepp port is another warning against a
lossy quotient.  The reflection identity `W_a(-E)=W_{7-a}(E)` and monotonicity
`W_a(E) <= W_a(E union -E)` both port cleanly, but consec is not a per-cell
maximizer: k=8..10 have shapes that beat it on individual `W_a` cells while
losing in the aggregate `measS7` sum.  Thus the LRC wall is genuinely
aggregate.  HYP-2770 should keep aggregate relation/moment channels, not try
to prove six independent per-cell inequalities.

Incoming KPS c14-lift work gives an external-method reason for the same
architecture, with an important correction: the state-of-the-art polynomial
method stalls because `14=2*7` is composite, but the paper's `c=14` lift is not
a literal multiplicative CRT factorization.  The corrected reading is sharper:
the missing step is certifying the canonical/consec tuple when the Fermat-field
argument no longer applies.  HYP-2770 is a compatible moment-hierarchy version
of that analytic replacement: use the seven-sector/full-residue quotient to
localize the canonical orbit, then retain the gap/relation data needed after
that quotient collapses.

The latest KPS decorrelated-wide note strengthens P3/P4 from the far-runner
side.  For decorrelated far elements, the cover probability is the closed-form
linear functional `sum_t P_t^(r) p_t(B)` of the base missed-depth distribution;
this is exactly the THM-534 moment-dual shape.  The new finite check proves the
decorrelated wide maximum for k=8..12 is always the `r=1` single-far/consec-base
plateau `Q(k-1) < cap_k`.  Thus the decorrelated main term is no longer a
mystery in this route; the remaining wide residual is the signed resonance
error, where HYP-2757 replaces the uniform weights by finite curve widths.
The next KPS S24 update tightens this further: against the moment-dual baseline,
the commensurable k=9 atlas has signed error at most about `0.01211`, only about
9 percent of the `Q(k-1)`-to-`cap_k` margin.  So the wide side is effectively
verified closed in the current atlas; the formal proof target is now a loose
joint-discrepancy/resonance-error bound, not a sharp constant hunt.

## Tournament Analysis

Vertices are proof lenses, not runners or arcs.  Pair observable:

```text
(pins_consec, avoids_circularity, preserves_gap_word, formal_ready, reusable)
```

The switch is lexicographic comparison.  The resulting tournament is transitive:

```text
affine_gap_word
  > relation_pair_marginals
  > full_residue_localization
  > signed_delsarte_dual
  > theta_prime_ceiling
  > signed_tanner_dessin
  > boolean_mobius_integrality
  > affine_residue_pairs
  > raw_tanner_expansion

score histogram = {0:1, 1:1, ..., 8:1}
directed 3-cycles = 0
SCC sizes = nine singleton SCCs
Hamiltonian path count = 1
```

## Proof Obligation

HYP-2770 proposes the next rigorous target:

```text
P0. Use HYP-2749 to discard non-full-residue shapes.
P1. Normalize the remaining shapes by translation and dilation on the
    seven-sector/full-residue quotient; do not claim a literal CRT shortcut for
    the paper's additive `c=14` lift.
P2. Retain integer gap/generated-word data, or equivalently the individual
    relation-code co-missing marginals.
P3. Prove the AP dilation orbit maximizes the level-2 relation-code bound on
    this affine-normalized full-residue quotient.
P4. Splice that extremality into the signed THM-534 Delsarte dual, using the
    KPS decorrelated-wide closure and tiny signed-error atlas to isolate the
    remaining wide obligation as a loose resonant-error bound.
```

The theta-prime / tournament-Paley bridge remains valuable as a shared ceiling:
it says the LRC consec problem and tournament Paley extremality are both
aggregate theta-ceiling problems.  It is not, by itself, a proof of the AP
orbit extremality.
