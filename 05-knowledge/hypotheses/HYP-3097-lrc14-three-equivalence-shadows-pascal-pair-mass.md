---
id: HYP-3097
title: LRC14 three equivalence shadows and Pascal pair-mass invariant
status: SYNTHESIS / invariant proposal; not a proof
source: codex-2026-06-27-S257
tangent: T1175
related:
  - HYP-2187
  - HYP-2998
  - HYP-2999
  - HYP-3000
  - HYP-3003
  - HYP-3085
  - HYP-3087
  - HYP-3089
  - HYP-3090
  - HYP-3091
  - HYP-3092
  - HYP-3093
  - HYP-3094
  - HYP-3096
  - THM-563
  - THM-576
  - OPEN-Q-108
---

# HYP-3097: LRC14 Three Equivalence Shadows And Pascal Pair-Mass Invariant

## Claim

LRC14 should track three different shadows before any scalar quotient is allowed
to stand in for the object:

```text
equidistribution:     measure / Weyl / density shadow
equinumerosity:       cardinal / Pascal / Burnside shadow
equidecomposability:  retained scissors fiber, sidecars, and obstruction data
```

This is the Pascal/pair-mass companion to HYP-3091's lonely-set fiber,
HYP-3092's verified pair-normalized cap mass, and HYP-3093's broad
equivalence-triad audit: HYP-3091 records the three-sameness invariant
`Phi(S)`, HYP-3092 verifies the cap-side pair mass, HYP-3093 names the
forgetting-cost protocol, and this note asks where the row-14 Pascal counts,
cap defects, and Farey pair-mass coordinates live inside that fiber.

HYP-2187 already separated equinumerosity from equidecomposability for
tournaments: equal counts or equal `H` values are only volume shadows unless the
`beta1`, odd-cycle, strong-component, or packet fiber survives.  The same
discipline now looks useful for LRC14.  Equidistribution proves that a large or
committed speed removes the expected fraction of a safe set; equinumerosity
counts how many bases or packet rows exist; equidecomposability asks whether two
rows with the same count or measure decompose into the same endpoint, sector,
pairwise-avoidance, level-7, and Farey-address pieces.

The proposed invariant is therefore not a number.  It is a packet:

```text
LRC_scissors_packet(S) =
  (measure shadow,
   cardinal/Pascal shadow,
   component and largest-arc data,
   sector-pair avoidance matrix,
   Pascal row or pair-mass address,
   Farey root (p,q), p+q, p*q,
   level-7 lift status,
   endpoint-owner / route sidecars,
   named destroyed coordinate and terminal exit).
```

Two rows are "equidecomposable" for a proof route only if the route's predicate
is constant after this packet is reattached, or if the lost coordinate is
reconstructed, dual-annihilated, descended, an AP/GW boundary atom, or named
residual debt.

## The 1001, 2002, 3003, 4004 Signal

The numbers in the prompt expose a pair-normalized Pascal structure.

```text
C(14,2) = 91

1001 = 11 * 91 = C(14,4) = C(14,10)
2002 = 22 * 91 = C(14,5) = C(14,9)
3003 = 33 * 91 = C(14,6) = C(14,8)
4004 = 44 * 91 = 2*C(14,5) = C(14,4)+C(14,6)
```

Thus `1001,2002,3003` are literal row-14 Pascal entries, while `4004` is the
affine pair-mass completion.  It is not itself a row-14 entry, but it is the
next value in the `11,22,33,44` multiples of the 14-runner pair count.

This matters because the same numbers already occur in the LRC14 work:

- THM-563's all-bounded-base finite check counts bases by binomial slices:
  `3003,3432,3003,2002,1001,364`.
- HYP-3090/HYP-3092 record the cap law
  `cap_k = C(k+1,2)/C(14,2)` for `k=10,11,12,13`.
- The two binding deviations are
  `cap_9 = 45/91 - 1/4004` and
  `cap_8 = 36/91 - 1081/76440`.
- HYP-3003's unit-excess Farey route has first product-wall packet
  `p/q=3/41`, with additive lane `p+q=44`, so `4004=44*C(14,2)` is also the
  pair-apex normalization of that first Kpq/K33-facing additive lane.

The conservative reading is: the denominator `91` is the pair apex, the row-14
Pascal entries count cardinal shadows of hidden obligation subsets, and the
factor `11,22,33,44` is a decomposition coordinate measuring how many pair
slices the shadow occupies.  The `k=9` defect `1/4004` is therefore not just a
small rational error.  It is a one-unit defect in the affine pair-mass
completion.

## Invariant Shape

The hidden invariant to test is:

```text
pair_mass_rank(E) = 91-normalized Pascal/Farey shadow
sector_scissors(E) = sector-pair avoidance matrix plus endpoint-owner data
cap_defect(E) = triangular_pair_cap(E) - actual_cap_or_p0_bound(E)
```

For `k>=10`, HYP-3090/HYP-3092 say the cap shadow is pure triangular pair
mass.  For `k=8,9`, the defect is the obstruction.  The new question is
whether those defects are exactly the Dehn-like scissors part of the LRC
object: invisible to equinumerosity, almost invisible to equidistribution, but
visible in the sector-pair decomposition and route sidecars.

## Relation To Past Work

This synthesis uses four existing threads.

1. **HYP-2187:** cardinal equality is not equidecomposition.  In LRC terms,
   equal caps or equal base counts need not preserve the sector and endpoint
   scissors fiber.

2. **HYP-2998/HYP-2999/HYP-3000/HYP-3003:** Fibonacci/Pascal rows and Farey
   `p+q`, `p*q` lanes are retained packet fields.  The row sum is a shadow;
   the row vector and carry/product fiber are proof data.

3. **HYP-3085/HYP-3090/HYP-3092/THM-576:** the cap side is pairwise
   avoidance.  The sector-pair matrix is the first real scissors decomposition
   of the cap shadow.

4. **HYP-3087/HYP-3089/HYP-3096:** the hyperoperation and polynomial-method
   ledgers say that `(p,q)`, `p+q`, `p*q`, CRT lift status, direct lonely
   components, and
   finite bad-denominator budgets must survive before an LRC witness route is
   legal.

## Proof Program

The next finite audit should add these fields to HYP-2963-style packet rows:

```text
pascal_pair_mass_unit
pascal_row14_shadow
triangular_cap_shadow
cap_defect_numerator
sector_pair_scissors_signature
same_cap_non_equidecomp_class
farey_additive_lane_mod_91
farey_product_lane_mod_91
level7_lift_status
destroyed_coordinate
```

Then test:

1. Rows with the same triangular cap shadow but different `sector_pair` matrix
   split into different route outcomes.
2. The `k=8,9` cap defects are carried by a small number of sector-pair
   scissors classes, not by raw base count.
3. The `4004=44*91` defect aligns with the first `p+q=44` Kpq/K33 product-wall
   lane only after exact `(p,q)` and endpoint owners are retained.
4. Equal base counts such as the two `3003` slices in THM-563 can be separated
   by decomposition data: base count alone is equinumerosity, not proof
   equivalence.

## S258 Back-And-Forth With The Witness Ledger

`04-computation/lrc14_observer_gluing_ledger_codex_s258.py` performs the first
joint audit with HYP-3096.  It does not merely print cap numerology; it attaches
pair/Pascal scissors fields to exact direct lonely-set rows:

```text
H7_pair_shadow = C(count_7_divisible,2) / C(14,2)
even_pair_shadow = C(count_even,2) / C(14,2)
mod7_residue_counts and mod14_residue_counts
Farey lanes p+q mod 91 and p*q mod 91 when an apex/root is named
```

The sample is deliberately small, but it already separates live rows that a
scalar residual label would merge.  The seven live observer-glue rows have one
terminal label, yet the sample has `5` distinct mod-7 scissors signatures.
The apex family `{1,...,12,14V}` keeps the same mod-7/mod-14 signature as `V`
changes while its direct largest arc falls from `3/637` at `V=13` to `3/9800`
at `V=200`; the divisor-loaded rows keep a different scissors signature and
shatter further to `1/82320` at `B=8`.  Thus:

```text
same terminal residual status
  != same component topology
  != same pair-scissors packet
  != same legal quotient for a proof route.
```

This confirms the invariant's role as a gluing sidecar.  Pair/Pascal mass
does not by itself prove the witness route; it tells the proof when the direct
arc chart, normalized arc chart, and moment/branch charts are still observing
the same packet.

Incoming kps-S31ag adds a compatible warning from the tournament side:
`lrc14_coverage_vs_tournament_H_kps.py` shows the coarse mod-14 winding
tournament degenerates for every `k>=8`, because `8` residues force an
antipodal pair `{r,r+7}` and hence a permanent tie.  Thus the coverage
extremality cannot be imported from coarse H-maximization at exactly the
binding rows.  In this HYP's language, coarse H is only another cardinal/count
shadow unless the fine-scale mod-`p` or packet scissors sidecar is retained.
This strengthens the S258 guardrail: observer vertices should be sector-pair
and fine-scale proof packets, not coarse mod-14 winding tournaments.

## Tournament Analysis

Vertices are invariant shadows and proof packets, not runners:

```text
sector_pair_scissors_packet
direct_lonely_component_packet
farey_root_sum_product_packet
pascal_pair_mass_rank
triangular_cap_shadow
node3_equidistribution_shadow
bounded_base_equinumerosity
raw_scalar_cap
```

Pairwise observable:

```text
(LRC predicate preserved,
 destroyed coordinate named,
 defect visible,
 finite-check compatibility,
 route/terminal-exit usefulness)
```

Gauge: coordinate majority.  Tie Hamiltonian path is the listed order.  The
declared fingerprint is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

The ranking is a guardrail, not a theorem.  It says the cap scalar is legal
only after the sector-pair scissors packet, component packet, and Farey address
have either survived or been paid for.

## Assumption Challenge

This pass explicitly considered alternate tournament vertices:

```text
runners, speed sets, danger arcs, sector boundaries, sector-pair events,
Pascal row entries, pair-mass units, Farey packets, exact denominator grids,
level-7 lift events, endpoint owners, components of L(S), proof obligations.
```

Chosen vertices are invariant/proof packets because the LRC predicate
`M(S)>=1/14` is not preserved by raw equinumerosity or raw equidistribution
alone.  The quotient preserves only a measure or count unless the sidecars are
reattached.  It destroys endpoint ownership, sector-pair correlations,
component topology, CRT lift status, and Farey product/sum fibers.

Challenged assumption: `4004,3003,2002,1001` should be read only as Pascal
triangle numerology.  The stronger reading is that Pascal counts are cardinal
shadows of a pair-apex scissors object, and `4004` is the affine completion
where the LRC cap defect is visible.
