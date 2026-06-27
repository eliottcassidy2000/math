# LRC14 Counterexample Families And Sporadics

Created: 2026-06-24T15:24:18Z
Author: codex-2026-06-24-S151
Hypothesis: HYP-2961

I pushed the counterexample hunt into a family/sporadic classifier.  The exact
gate is now:

`TRUE-COUNTEREXAMPLE-CANDIDATE = qdiv>14 && exact_threshold_status == covered`.

Everything outside that gate has a discharge reason: q-clock witness, AP/GW
boundary tightness, q14 positive migration, unit-petal/GW strip, K33/state-lift
packet, single-14-tail comb, or unlabelled covering repair with positive mass.

The S151 audit covered AP, one AP swap through `add<=420`, two swaps through
`add<=60`, and three swaps through `add<=30`.  Exact `qdiv>=14` rows audited:
`68368`.  True counterexample candidates: `0`.

Family counts:

- AP/GW boundary: `2`.
- K33 state-lift: `338`.
- mixed K33/petal lift: `2`.
- q14 positive migration: `41906`.
- single-14-tail comb: `235`.
- unit-petal/GW strip: `9632`.
- unlabelled covering repair: `16253`.

The last bucket is the interesting one.  It is the sporadic reservoir, not the
sporadic counterexample bucket.  Examples like `drop(4,6)->add(19,42)` and
`drop(2,6)->add(17,42)` are unlabelled distributed-cover repairs, but the exact
interval classifier still finds positive safe mass.  So the next theorem should
try to convert their divisor-cover skeleton directly into a witness interval.

Tournament readout: vertices are classifier gates/families rather than runners.
The observable is "which quotient preserves the strict-counterexample predicate
and adds the strongest discharge label?"  The resulting tournament is
transitive with path:

`qclock-excluded > q14-boundary-tight > q14-open-migration > covering-comb-family > unit-petal-family > K33-state-lift-family > unlabelled-covering-repair > true-covered-sporadic`.

This is not the proof.  It is a way to make future proof attempts honest: any
proposed counterexample now has to state which gate it escapes, whether it is a
family member, or why it is the first genuine sporadic.
