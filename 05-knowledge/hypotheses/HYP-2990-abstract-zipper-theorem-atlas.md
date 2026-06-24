---
id: HYP-2990
title: Abstract zipper theorem atlas
status: EXPLORATORY PROOF-TECHNOLOGY / cross-domain theorem templates; not a proof
source: codex-2026-06-24-S165
related:
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2984
  - HYP-2969
  - HYP-2945
  - HYP-2942
  - HYP-2934
  - HYP-2932
  - HYP-2887
  - HYP-2664
  - HYP-2620
  - HYP-2618
  - HYP-1815
  - THM-354
  - THM-408
  - OPEN-Q-108
results:
  - 04-computation/abstract_zipper_theorem_atlas_codex_s165.py
  - 05-knowledge/results/abstract_zipper_theorem_atlas_codex_s165.out
---

# HYP-2990: Abstract Zipper Theorem Atlas

HYP-2990 records the cross-domain version of the HYP-2987 zipper target.  A
zipper theorem is not a claim that two scalar summaries are similar.  It is a
controlled-kernel theorem:

```text
a quotient may forget a coordinate only when the target predicate is constant
on fibers, the coordinate is reconstructible from retained labels, a dual
certificate annihilates the forgotten direction, or the direction is routed to
a named residual sector.
```

The useful abstraction is two toothed carriers and a slider.  The left teeth
and right teeth each certify local predicates.  The slider identifies their
packet fibers.  The theorem is valid only if the slider keeps the labels needed
to transport certificates across the teeth.

## Computation

Script:

```text
04-computation/abstract_zipper_theorem_atlas_codex_s165.py
```

Stored output:

```text
05-knowledge/results/abstract_zipper_theorem_atlas_codex_s165.out
```

The script compares eleven carrier candidates, including a negative control:

```text
lrc14_certificate_handoff_zipper
kernel_tope_smoothing_zipper
octahedral_hodge_current_zipper
c27_unital_pair_completion_zipper
farey_product_k33_zipper
boundary_moment_multichart_zipper
shell1_root_packet_mouth_zipper
unit_distance_spine_ear_zipper
ocf_activity_coimage_zipper
good_cut_scc_support_zipper
raw_scalar_shadow
```

Tournament vertices are proof carriers, not runners or arcs.  The pairwise
observable is majority comparison over:

```text
predicate retention
fiber labels
two-sided transfer
declared kernel
finite checkability
family compression
formalizability
anti-scalar guardrail
cross-domain signal
```

## Tournament Readout

The S165 fingerprint is:

```text
score_hist={0: 1, 1: 1, 2: 1, 4: 2, 5: 1, 6: 1, 7: 2, 9: 1, 10: 1}
directed_3cycles=4
SCC_sizes=[1, 1, 1, 6, 1, 1]
Hamiltonian_path_count=15
```

The top carriers are:

```text
lrc14_certificate_handoff_zipper  score 10
kernel_tope_smoothing_zipper      score  9
octahedral_hodge_current_zipper   score  7
c27_unital_pair_completion_zipper score  7
farey_product_k33_zipper          score  6
```

The nontrivial six-carrier SCC is the important middle layer.  Octahedral
currents, C27/unital pair completion, Farey/K33 incidence, boundary moments,
unit-distance ears, and good-cut/SCC support all dominate one another in
different retained-label directions.  That means the next theorem should not
linearize them prematurely.  It should define typed handoffs.

## Zipper Templates

**Labelled packet zipper lemma.**  If two carriers `L` and `R` map to the same
labelled packet base `P`, the target predicate is constant on `P`-fibers, and
every kernel coordinate is reconstructed, annihilated, or routed to a named
residual sector, then local `L`- and `R`-certificates glue to a certificate over
the fiber product `L x_P R`.

**No-free-slider lemma.**  A proposed zipper fails exactly where the slider
forgets a load-bearing label.  A failure certificate should name one of four
defects: predicate mixing on a fiber, unreconstructed coordinate, unannihilated
kernel, or unnamed residual sector.

**Alternating gauntlet theorem.**  Run left teeth and right teeth alternately.
If every failed tooth emits a smaller named packet or a boundary-moment /
state-lift debt, and the debt order is well founded, the gauntlet terminates in
a certificate, a known equality atom, or a forbidden lift.

**Harmonic residual theorem.**  The last bucket cannot be an anonymous
exception.  It must be a named representation, homology, cocircuit, curl,
Johnson-harmonic, or state-lift sector with an explicit predicate.

## Cross-Domain Carriers

The past-topic connections are now more than analogies:

- HYP-2987 is the main labelled-packet zipper: source packets zip to Fejer,
  Ramanujan, endpoint, twist, moment, analytic-sieve, and state-lift carriers.
- HYP-2984/HYP-2985/HYP-2986 form a local zipper between open topes, boundary
  cocircuits, kernel support radii, and admissible smoothing exits.
- HYP-2887 is a divergence/curl zipper on the octahedral `L(K4)` support.
- HYP-2942 is a C27/unital pair-completion zipper; global claims fail if the
  q=3 branch label is forgotten.
- HYP-2932/HYP-2934/HYP-2945 give the Farey/K33 incidence zipper; product data
  is useful only after exact `M`, denominator shell, factor fiber, and minor
  labels are retained.
- HYP-2969 gives the boundary-moment multi-chart zipper; one all-covered
  denominator chart is not an obstruction.
- HYP-1815/HYP-2664 give a shell-1/root-packet zipper; carry gates must precede
  root-packet or comb scalarization.
- HYP-2620/THM-408 give the unit-distance endpoint-ear zipper; graph traceability
  alone forgets the unit-cocyclic geometry that can still block an extension.
- HYP-2618 gives the OCF activity/coimage zipper; `H(T)=I(Omega(T),2)` should
  be transported through compatible packet addresses, not read as raw noise
  stability.
- THM-354/INV-237 give the good-cut/SCC support zipper; the base-path interval
  coordinate becomes `g_P(T)=n-#SCC(T)` only after the base path and component
  boundary order are declared.

## LRC14 Transfer

The strongest next LRC14 move is a typed packet theorem, not another scalar
ranker.  Fejer intervals, Ramanujan exact-period projectors, endpoint bridges,
tope/cocircuit walls, boundary moments, and state-lift obligations should be
declared as projections of one labelled packet object.

The finite obstruction target becomes:

```text
F7 = the named harmonic/state-lift residual sector left after every zipper
tooth either certifies a strict interval or emits AP/Goddyn-Wong equality.
```

This strengthens HYP-2987's O6 obligation.  `F7` must state which predicate it
preserves and what information it destroys.  Otherwise it is not a theorem
residual; it is just an untyped failure bucket.

## Next Work

1. Formalize the no-free-slider lemma as a checklist for every LRC14 quotient:
   invariant, reconstructible, annihilated, or named residual.
2. Give HYP-2987's `F7` a concrete harmonic/state-lift definition.
3. Build a family-compression prototype where one zipper certificate covers a
   whole K33 or petal packet family rather than a named row.
4. Test whether the good-cut/SCC support residue can become a packet-coordinate
   for LRC14 state-lift construction.
5. Revisit the octahedral current route with the zipper language: divergence,
   curl, and tail labels should be separate teeth, not scalar current mass.
