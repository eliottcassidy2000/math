---
id: HYP-3042
title: LRC14 owner-strip filtration after residual repairs
status: SYNTHESIS / necessary-condition sharpening; not a proof
source: codex-2026-06-26-S205
tangent: T1123
script: 04-computation/lrc14_owner_strip_spectral_sequence_codex_s205.py
result: 05-knowledge/results/lrc14_owner_strip_spectral_sequence_codex_s205.out
related:
  - HYP-3041
  - HYP-3038
  - HYP-3037
  - HYP-3036
  - HYP-3035
  - HYP-3031
  - HYP-3018
  - HYP-2997
  - HYP-2963
  - THM-572
  - LTI-190
  - LTI-189
  - LTT-088
  - LTT-087
  - T1122
  - OPEN-Q-108
---

# HYP-3042: LRC14 Owner-Strip Filtration

## Claim

The recent residual stack suggests a short filtration rather than another
scalar invariant:

```text
raw shadow
  -> status gate
  -> primitive-period deck
  -> Haar/drop-add zeta
  -> endpoint-owner strip current
  -> labelled route certificate
```

The hidden detail is that `endpoint-owner strip` is not a new ad hoc label.
It is the common boundary current already visible in older work:

- HYP-2997 split the Haar `zeta` channel from the endpoint-owner boundary
  cocycle.
- HYP-3018 retained active bottleneck endpoint owners and residue sums as a
  normal-fan support.
- HYP-3035 found that all `15` coarse ET+unit residual first repairs are
  `owner_strip` repairs.
- HYP-3036 showed primitive-period decks schedule direct q-witness residuals
  but leave the boundary/covering layer as a separate page.
- HYP-3041 showed that even an owner-strip collision can be repaired by a
  hidden primitive clock: the AP-tail rows with the same mod-14 owner strip
  split after reattaching `q13_puncture_bit` / `ap_tail_certificate_kind`.
- HYP-3038 showed the `q=23` petal/covering diagonal needs nonzero
  Haar/drop-add zeta first, and then endpoint-owner strips split the two
  diagonal routes even though the coarse endpoint word is the same `B18Z6`.

Thus the sharpened necessary condition is:

```text
Any LRC14 counterexample surviving the current protected residual stack must
simultaneously have
  no positive primitive safe mass at q <= 13,
  no AP-tail q=13 puncture / reciprocal fixed-point clock,
  no useful drop/add Haar-zeta square that opens or descends,
  and no distinguishing endpoint-owner strip current
before it can be treated as genuine named F7/THM-572 debt.
```

## Computation

Script:

```text
04-computation/lrc14_owner_strip_spectral_sequence_codex_s205.py
```

Stored output:

```text
05-knowledge/results/lrc14_owner_strip_spectral_sequence_codex_s205.out
```

The script is a connection-mining table over HYP-3041, HYP-3038, HYP-3037,
HYP-3036, HYP-3035, HYP-3031, HYP-3018, and HYP-2997.  It does not recompute
the full packet bank.  It records which coordinate each thread had already
named and then builds a tournament over filtration pages.

## Definitions

An `endpoint-owner strip current` is the boundary-owner multiset that survives
after a coarse endpoint word has forgotten identity.  In S201, both q=23
diagonal rows have the coarse endpoint-current word `B18Z6`, but the external
owner multisets differ:

```text
petal diagonal:    12:26x6, 6:20x4
covering diagonal: 2:16x6
```

So `B18Z6` is only a scalar shadow.  The strip current must retain which
external speed/residue owns the boundary facets and with what multiplicity.

The filtration pages are:

```text
E0 raw_shadow
E1 status_gate
E2 primitive_period_deck
E3 haar_zeta_grid
E4 endpoint_owner_strip
E5 labelled_route_certificate
```

They should not be read as a formal spectral sequence yet.  The useful analogy
is that each page kills one class of quotient ambiguity and passes only the
remaining residual current to the next page.

HYP-3041 is the warning example for `E2`: a mod-14 owner strip can still collide
when it forgot `m mod 13`; the quotient becomes theorem-safe only after
`q13_puncture_bit` or `ap_tail_certificate_kind` is reattached.

## Tournament Analysis

Vertices are filtration layers / proof carriers, not runners.

Pairwise observable:

```text
status, route, primitive_period, haar_zeta, owner_strip,
topology, family_transfer, compression, low_cost
```

Switch:

```text
majority of retained proof coordinates
tie path =
labelled_route_certificate > endpoint_owner_strip > haar_zeta_grid
> primitive_period_deck > status_gate > raw_shadow
```

Output fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1]
hamiltonian_path_count=1
path=endpoint_owner_strip>labelled_route_certificate>haar_zeta_grid>
     primitive_period_deck>status_gate>raw_shadow
```

The path ranks `endpoint_owner_strip` before labelled route because the gauge
rewards non-route legal splitting and compression.  Route labels are complete
certificates, but they are the last page rather than the first proof object.

## Next Pull

Add these sidecars to a cached HYP-2963 packet ledger:

```text
primitive_safe_deck_2_13
drop_add_square_id
exact_M_zeta
endpoint_owner_strip_current
owner_strip_page
first_surviving_filtration_page
```

Then search for residual packets whose first surviving page is beyond
`endpoint_owner_strip`.  Those are the candidates that deserve named F7,
THM-572, harmonic, or state-lift attention.

## Assumption Challenge

Alternate vertices considered: runners, residual fibers, denominators, drop
pairs, add pairs, endpoint walls, safe bars, Haar squares, primitive periods,
cocycle channels, normal-fan supports, and proof obligations.

The chosen vertices are filtration pages because the preserved predicate is
strict-open status plus route schedulability after quotienting.  This destroys
raw runner identity, exact route labels before the last page, and scalar
endpoint words unless the owner multiset is retained.
