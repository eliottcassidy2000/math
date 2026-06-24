---
id: HYP-2954
title: LRC14 missing-picture bridge
status: PROOF-INTERFACE / functorial bridge target, not a proof
source: codex-2026-06-24-S149
related:
  - HYP-2953
  - HYP-2952
  - HYP-2951
  - HYP-2950
  - HYP-2949
  - HYP-2948
  - HYP-2947
  - HYP-2946
  - HYP-2944
  - HYP-2942
  - HYP-2940
  - HYP-2937
  - HYP-2928
  - HYP-2910
  - HYP-2909
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_missing_picture_bridge_codex_s149.py
  - 05-knowledge/results/lrc14_missing_picture_bridge_codex_s149.out
---

# HYP-2954: LRC14 Missing-Picture Bridge

## Claim

The current LRC14 proof history is no longer missing a scalar invariant.  It is missing a functor:

```text
primitive reduced residual
  -> exact M/Farey branch
  -> Haar-Baire open-or-boundary front
  -> C27/unital/K33 owner address
  -> discharge, AP/GW boundary atom, covering strictness, or TournamentStateLift.
```

The functor has to preserve the LRC14 witness predicate: exact scale, strict-open Haar mass, closed endpoint debt, C27 hole/double owner labels, q=3 unital chart position, and K33/state-lift obligation.  It may forget raw runner names only after those labels have been recorded.

## Computation

`04-computation/lrc14_missing_picture_bridge_codex_s149.py` reuses the exact kernels from S124 and S146:

- S124 supplies exact `M(S)`, argmax denominators, and `q_threshold`.
- S146 supplies strict-safe interval fronts, Haar mass, endpoint boundary owners, and bounded AP-neighborhood scans.

The named-row audit separates the low frontier into proof routes:

| row | exact readout | route |
| --- | --- | --- |
| AP | `M=1/14`, `q=14`, strict Haar mass `0` | boundary atom |
| GW `12->24` | `M=1/14`, `q=14`, strict Haar mass `0` | boundary atom |
| near `12->36` | `M=3/41`, strict mass `1/1260`, packet `K33` | K33 state-lift obligation |
| petal `10->20` | `M=2/27`, strict mass `1/980`, packet `P10` | unit-petal discharge |
| petal `13->26` | `M=2/27`, strict mass `1/182`, packet `P13` | unit-petal discharge |
| splice `(10,12)->(20,24)` | `M=2/27`, packet `P10+GW` | unit-petal discharge |
| splice `(10,12)->(20,36)` | `M=2/27`, packet `P10+K33` | K33 state-lift obligation |
| residue liar `12->26` | same AP residues but `q=12`, `M=1/12` | q-witness loose |
| magnitude liar `12->96` | `q=14` but `M=8/101` | magnitude-aware open front |
| covering comb `12->84` | `14`-covering, `M=7/89` | covering strictness |

The reused S146 local database remains the hard boundary calibration:

- one-swap AP neighborhood with added value `<=160`: only AP and GW are boundary-only;
- two-swap AP neighborhood with added values `14..40`: no boundary-only rows;
- the smallest two-swap positive rows are exactly the S138 splices.

Rebase integration with the incoming S148/S149 work strengthens the same picture from adjacent sides.  HYP-2950's gauntlet found no below-threshold row among named adversaries, AP single-swaps through `v<=360`, or AP two-swaps through added values `<=42`, and the low `M<=2/27` packet remains AP, GW, K33 near-miss, unit petals, `P10+GW`, and `P10+K33`.  HYP-2952's derived-boundary tournament classes show that AP, GW, the near/K33 row, petals, and splices can share the same transitive six-unit apex-pressure class, so the apex tournament is a necessary pressure filter, not the missing invariant.  HYP-2953's source-spectrum pullback names the same need in source-fiber language: the functor here is the computable packet ledger whose pullback should preserve Farey source, Haar/Baire boundary status, and C27/unital/K33 labels.  That is exactly why the bridge functor must retain exact scale and owner labels.

## Tournament Analysis

Vertices are proof interfaces/quotients, not runners.  Pairwise observable: branch retention, boundary retention, C27/unital retention, K33/state-lift fit, finite checkability, theorem endpoint strength, and anti-scalarization guard.  The switch is componentwise score sum, with ties broken by the declared Hamiltonian path.

The proof-interface tournament is transitive:

```text
C27_transfer
> exact_M_Farey
> unital_pair_chart
> K33_incidence
> Haar_Baire_front
> TournamentStateLift
> boundary_owner_skeleton
> affine_depth
> raw_scalar_M
> raw_apex_tournament
```

This ranking is not a proof, but it states the missing-interface order.  Raw apex tournaments are magnitude-blind: AP and residue-liars can share residue data while `q_threshold` separates them.  Raw scalar `M` is also insufficient because it does not say which owner packet produced the mass.  The bridge must carry the tuple.

## Missing Theorem

Every primitive reduced LRC14 residual not discharged by q-threshold or Haar-open strictness either:

1. is AP/GW boundary-only after C27/unital labels are retained;
2. constructs the HYP-2908 / THM-572 `TournamentStateLift`;
3. lies in the covering/shell branch and is strict by comb-margin or shell-collapse.

Equivalently: the unresolved arrow is not "find a better count."  It is the quotient-preservation theorem sending any bad residual into a typed packet category whose endpoint is already known or sharply isolated.

## Assumption Challenge

Tournament vertices need not be runners, arcs, or an isomorphism class of the raw phase tournament.  In this bridge the plausible vertices are:

- strict-safe interval fronts;
- boundary owner pairs;
- C27 hole/double transfers;
- q=3 unital chart points and blocks;
- K33 incidence obligations;
- covering comb teeth;
- exact-period packet states;
- support-six relation packets;
- proof obligations in the HYP-2908 state lift.

The quotient preserves the LRC14 witness predicate only if it records exact scale, boundary/open status, and owner labels.  It destroys raw speed identities and must therefore not be used before the `q_threshold` and C27/Farey labels are attached.
