---
id: THM-1153
title: The literal covering no-Lonely13 dominance premise is refuted by the doubled tight AP
status: PROVED — the positive family W=2*{1,...,13} covers every modulus 2,...,14, has no 1/13-lonely time because dilation preserves the exact AP maximum 1/14, and has largest/second-largest ratio 26/24<13.  Thus the Lean proposition INVcov is false.  Its conditional consumer theorems remain kernel-valid, while ResidualINV and residualINV_iff_LRC14 remain the honest exact interface.  A dedicated Lean module checks the counterexample and `not_INVcov` without sorry
source: codex-2026-07-18-S67 continuation, counterexample identified by concurrent S74 frontier audit
depends_on: [THM-358]
related: [THM-1131, THM-1149, MISTAKE-170]
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCINVcovCounterexample.lean
---

# THM-1153 — `INVcov` is false

The formal proposition introduced in the corrected M-split was

```text
INVcov:
  positive(v) and Covering(2..14,v) and no Lonely13(v)
  => some speed 13-dominates every other speed.        (1)
```

Adding `Covering(2..14)` repaired the earlier missing-modulus error, but it
did not make (1) true.  The obstruction is dilation.

## Exact counterexample

Let

```text
W=2*{1,...,13}={2,4,6,...,26}.                         (2)
```

All speeds are positive and distinct.

### Covering

For every `2<=q<=13`, the speed `2q` belongs to `W`, so `q|2q`.  Modulus
`14` divides the speed `14`.  Hence `W` is literally `Covering(2..14)` in
the sense of `LRCSieveDispatch.lean`.

### No `Lonely13` time

The initial segment `{1,...,13}` has exact lonely maximum `1/14`: the unit
time `1/14` attains it, and Dirichlet's 14-point pigeonhole gives the matching
upper bound at every time.  Multiplying every speed by two only reparametrizes
time, so

```text
M(W)=M({1,...,13})=1/14<1/13.                          (3)
```

Thus no time is `1/13`-lonely for `W`.

### No dominance

The two largest speeds are `26` and `24`, and

```text
13*24>26.                                              (4)
```

No speed can 13-dominate all the others.  Equations (2)--(4) refute (1).

## What remains formally sound

The theorems

```text
LRC14_of_INVcov
LRC14_finset_of_INVcov
residualINV_of_INVcov
```

are valid implications from a false premise.  Lean did not prove a false
conclusion; the mathematical mistake was calling their premise an open
supplier.  They are retained as historical conditional lemmas and labelled
accordingly.

The exact interface

```text
ResidualINV:
  positive(v) and Covering(v) and no Lonely14(v)
  => 13-fold dominance
```

is unaffected.  Under the cited AP bridge,
`residualINV_iff_LRC14` proves that it is equivalent to the working LRC(14)
statement.  It is therefore a faithful counterexample language, but not a
noncircular reduction.

A potentially useful stronger supplier must normalize dilation first.  A
primitive version of (1) is not refuted by (2), but it is not enough merely to
write `gcd=1`: one must prove that gcd reduction preserves absence of a lonely
time and then **rederive** Covering for the reduced family.  Covering itself is
not dilation-invariant—`{1,...,13}` misses modulus 14 whereas (2) covers it.
The existing `exists_lonely_smul` theorem supplies the first ingredient; the
complete primitive supplier/consumer bridge remains to be stated and proved.

## Structural lesson and Tournament Analysis

The speed-order tournament of `{1,...,13}` is identical to that of its
dilation (2): both are transitive, have score histogram `{0,...,12}`, no
directed cycles, thirteen singleton SCCs, and one Hamiltonian path.  Yet
Covering changes and the proposed inverse premise changes truth value.  Thus
the runner tournament destroys precisely the gcd/divisor-fibre data at issue.

Candidate vertices challenged here are runners, speed ranks, residue classes,
divisor obligations, gcd fibres, and deletion crowns.  The faithful quotient
for this correction must retain the common dilation factor together with the
`q=2,...,14` divisibility obligations.  A tournament on ranks or residues
cannot distinguish the counterexample from its non-Covering primitive core.

## Formal verification

`LRCINVcovCounterexample.lean` defines (2), proves positivity and Covering,
uses the exact AP tightness theorem plus dilation invariance to rule out a
`Lonely13` time, proves failure of every dominance choice, and concludes
`not_INVcov`.  The module is audited for axioms and contains no `sorry`.
