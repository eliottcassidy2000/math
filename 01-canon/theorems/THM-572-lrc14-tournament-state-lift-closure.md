---
id: THM-572
title: LRC(14) tournament-state lift closure -- any bad atom that lifts to H=7 is impossible
status: FORMALIZED in Lean; conditional on constructing the state lift
source: codex-2026-06-22-S128
depends_on:
  - THM-343   # H(T) != 7 for all tournaments
  - HYP-2908  # LRC14 forbidden-H7 state-lift target
related:
  - HYP-2907
  - HYP-2910
  - HYP-2924
  - HYP-2930
  - THM-079
  - THM-201
  - OPEN-Q-108
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCTournamentStateLift.lean
  - 04-computation/lean/TournamentH7/TournamentH7/Verify.lean
---

# THM-572 -- LRC(14) tournament-state lift closure

## Statement

Let `Bad` be any proposed LRC14 bad-atom predicate.  If every witness of
`Bad` constructs a tournament-state packet whose activity-two packet value is
`7` and agrees with the tournament Hamiltonian-path count `H`, then `Bad` is
empty.

Equivalently, if a proposed bad atom constructs a tournament `T` with

```text
H(T) = 7,
```

then the bad atom is impossible.

## Lean formalization

The Lean module `LRCTournamentStateLift.lean` exposes three reusable closure
forms.

```lean
theorem no_tournament_state_lift
    (L : LonelyRunner.TournamentStateLift) : False

theorem not_bad_of_tournament_state_lift {Bad : Prop}
    (hLift : Bad -> LonelyRunner.TournamentStateLift) : Not Bad

theorem not_bad_of_H_eq_seven_lift {Bad : Prop}
    (hLift : Bad -> Exists fun n => Exists fun T : Tournament n => H T = 7) :
    Not Bad
```

The proof is a direct application of the already formalized THM-343 theorem
`Tournament.H_ne_seven`.

## Proof sketch

For a `TournamentStateLift`, the fields give

```text
packetValue = H(T)
packetValue = 7.
```

Thus `H(T)=7`, contradicting THM-343.

## Proof role in LRC(14)

This theorem does **not** construct the missing LRC-to-tournament state lift.
It pins down the endpoint of the tournament proof route:

```text
LRC14 bad atom
  -> TournamentStateLift
  -> H(T)=7
  -> contradiction by THM-343.
```

After THM-568/HYP-2929 and THM-571, the best target for the first arrow is the
remaining `|M14|<=6` scale-separated / finite-core covering strictness branch,
or any shell-height `h>1` atom left by apex-shell collapse.  Such an atom must
be lifted into complete tournament-conflict data, not arbitrary digraph data;
HYP-2907 remains the guardrail showing that loose digraph categories can
realize value `7`.

## Assumption challenge

The theorem deliberately keeps `Bad` abstract.  This avoids assuming that the
right tournament vertices are runners.  Candidate vertices for the missing
lift remain sector-state words, wall-crossing packets, cover-arc packets,
exact-period packets, support-six relation packets, or finite proof
obligations.  The quotient must preserve the activity-two packet value and its
agreement with `H`; raw tournament isomorphism alone is insufficient by
HYP-2924.
