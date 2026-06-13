---
id: HYP-1992
status: OPEN
source: codex-2026-06-01-S520
related:
  - HYP-1991
  - HYP-1990
  - HYP-1989
  - HYP-1988
  - HYP-1986
  - HYP-1987
  - HYP-1981
  - HYP-1982
  - HYP-1839
  - HYP-1840
  - HYP-1952
  - HYP-1985
  - THM-381
  - THM-382
  - THM-383
  - THM-384
  - THM-385
  - THM-386
  - THM-377
  - THM-380
---

# HYP-1992: The n=18 battlefield inherits the observer-source gate/tightness split

## Statement

THM-382, THM-383, HYP-1981, THM-384/HYP-1986, HYP-1987,
THM-385/HYP-1988, and THM-386 turn the n=18 gate battlefield into a
threshold-decorated observer-source target problem:

```text
no 18-gate branch:     all unit points a/18 remain observer-source targets;
has 18-gate branch:    unit source targets are killed, but endpoint debt is exported.
```

Thus an `n=18` open-cover counterexample must contain a speed divisible by
`18`, while proof work in that branch should focus on descendant endpoint
source layers, especially the `9`-ladder, `18`-gate, and `36`-double-gate
layers.

The older gate/tightness lesson is still present: the no-gate branch is
dismissed by unit witnesses.  HYP-1981 sharpens what those witnesses mean:
they are marked tournament source targets for the stationary observer.
THM-384 sharpens the local target further to the two adjacent observer gaps,
and HYP-1987 says the global target lives in the tiny arc-confined source menu
inside `A000568(17)`.  THM-385/HYP-1988 adds the adjacent repair ledger:
observer indegree is exactly blocker count, so the gate branch should next be
stratified by almost-source and higher blocker layers.
THM-386 adds the directed-flow constraint: true source-gap entry is an
`LS -> LL` gap race rather than an average gap-sum phenomenon.

## Evidence

`lrc_n18_gate_battlefield_s520.py` audits the unit skeleton, canonical rows,
targeted one-step repairs, pure gate replacements, HYP-1981 observer-source
fingerprints, and a representative branch tournament.

The residue action is exact:

```text
unit points: 1/18, 5/18, 7/18, 11/18, 13/18, 17/18
residue 0 mod 18 covers all six;
every other residue tested covers none.
```

Canonical rows show the split:

```text
initial: boundary_only, gate=0, unit_safe=6, unprotected=6
initial 8->18: positive_gap, gate=1, unit_safe=0
initial 17->18: positive_gap, gate=1, unit_safe=0
```

The S383 lpd/gate ladder chain is especially rigid:

```text
9-ladder skip 8:
  gap/th=1/176, unprotected=176, first=11/162

18-gate ladder skip 8:
  gap/th=1/352, unprotected=352, first=19/324

36-double-gate skip 8:
  gap/th=1/704, unprotected=704, first=37/648
```

So the dyadic gate lift halves the normalized gap, doubles unprotected
endpoint debt, and doubles the first descendant denominator.

The targeted one-step repair scan around coordinates

```text
6, 8, 9, 12, 16, 17
```

found only positive-gap repairs.  In the no-gate half, the best rows keep
all six unit witnesses safe and unprotected:

```text
swap 16->48: gate=0, unit_safe=6, unit_unprotected=6, gap/th=0.020833
swap 12->24: gate=0, unit_safe=6, unit_unprotected=6, gap/th=0.029762
```

In the gate half, the best rows kill the unit witnesses but keep positive
gaps and descendant debt:

```text
swap 8->108: gate=1, unit_safe=0, gap/th=0.021044, first=37/306
swap 17->18: gate=1, unit_safe=0, gap/th=0.033333, first=19/324
```

The branch Tournament Analysis is transitive:

```text
H=1, c3=0, SCCs=(1,1,1,1,1,1,1,1,1)
top rows:
  double-gate 36* skip 8
  gate ladder 18* skip 8
  lpd ladder 9* skip 8
  swap 8->108
```

This says the sampled branch pressure is an ordered endpoint-debt ladder, not
a cyclic repair mechanism.

The observer-source fingerprints supply the HYP-1981 reading.  With observer
edges oriented by

```text
observer -> runner  iff  ||v*t|| >= 1/18,
```

the representative witness times all have observer outdegree `17/17`.  The
initial tight row reaches source targets only at boundary unit points, so
THM-383's tie-wall compactification is essential.  Gate rows destroy those
unit source targets but still create open source targets at descendant
endpoint denominators, exactly the threshold-decorated behavior THM-382 says
must be retained.

## Predictions

1. The no-`18`-multiple branch cannot contain an open-cover counterexample.
   It is dismissed by the six unit observer-source targets.
2. The `18`-gate branch should be attacked by an endpoint-debt certificate,
   not by speed-set search.
3. The `9 -> 18 -> 36` chain should obey a dyadic debt law:

```text
gap/th halves,
unprotected endpoints double,
first descendant denominator doubles.
```

4. Any credible n=18 counterexample search must first produce a marked
   tournament walk avoiding every observer-source target while also carrying a
   nonpeeling endpoint-pressure core inside this exported debt.  The sampled
   branch tournament found no cyclic repair shape.

## Sources

- `04-computation/lrc_n18_gate_battlefield_s520.py`
- `05-knowledge/results/lrc_n18_gate_battlefield_s520.out`
- `07-reflections/lrc-n18-observer-source-gate-battlefield-s520.md`
- HYP-1991
- HYP-1990
- HYP-1989
- HYP-1988
- HYP-1986
- HYP-1987
- HYP-1981
- HYP-1982
- HYP-1839
- HYP-1840
- HYP-1952
- THM-381
- THM-382
- THM-383
- THM-384
- THM-385
- THM-386
- THM-380
