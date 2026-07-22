---
id: THM-2050
title: "Period-14 top germs do not determine global loneliness"
status: >
  PROVED. AP13 and the loose same-residue lift 12->26 have identical
  phase-height germs on explicit neighborhoods of all six unit points a/14,
  yet their global lonely-runner maxima are respectively 1/14 and 1/12.
  Therefore localization at the period-14 top layer, even with its complete
  local function germ, cannot prove LRC14 without a magnitude/first-exit
  sidecar.
source: codex-2026-07-21-DC2-LRC14-termination
related:
  - THM-523
  - THM-597
  - THM-2043
  - THM-2047
  - THM-2048
  - HYP-8841
---

# THM-2050 -- local top-germ blindness at period fourteen

Let

```text
A={1,2,...,13},
L={1,2,...,11,13,26},
f_S(t)=min_(v in S) ||vt||.                          (1)
```

Then

```text
M(A)=1/14,                 M(L)=1/12.               (2)
```

Nevertheless, for every `a` coprime to `14`,

```text
f_A(t)=f_L(t) whenever |t-a/14|<1/728.              (3)
```

Thus the two phase-height arrangements have identical function germs at all
six period-14 unit maxima.

## Proof of the global values

For any `t`, the fourteen circle points

```text
0,t,2t,...,13t
```

contain two at distance at most `1/14`; their difference gives some
`1<=v<=13` with `||vt||<=1/14`.  Equality is attained at `t=1/14`, proving
`M(A)=1/14`.

The same pigeonhole argument on `0,t,...,11t` gives

```text
M({1,...,11})=1/12.                                  (4)
```

Since `{1,...,11}` is a subset of `L`, (4) gives `M(L)<=1/12`.  At `t=1/12`,
the speeds `1,...,11,13` have minimum distance `1/12`, while
`||26/12||=1/6`.  Hence `f_L(1/12)=1/12`, proving (2).

## Proof of germ equality

Fix a unit `a mod 14` and put `t_0=a/14`.  In both sets, the binding speeds are
the same unit-residue pair; their values at `t_0` equal `1/14` and their speeds
are at most `13`.  Every nonbinding speed has value at least `1/7` at `t_0`.
This includes the differing speeds `12` and `26`, which have the same nonzero
even residue modulo `14`.

If `|h|<1/728`, a binding tent has value at most

```text
1/14+13/728=65/728,                                  (5)
```

whereas every nonbinding tent in either set has value at least

```text
1/7-26/728=78/728.                                   (6)
```

Thus only the common binding tents can realize the minimum on this
neighborhood.  Their functions are identical in `A` and `L`, proving (3).

## LRC14 consequence

The entire local collapse-rate data of THM-597 at the six unit points is the
same for these two rows, as are all local derivatives and all sufficiently
small localized phase-height cells.  Yet `L` has the strict global exit
`t=1/12`.

Therefore any LRC14 route based only on localization at the denominator-14
top layer is incomplete.  It must retain at least one global termination
sidecar: magnitude, first strict-exit denominator, q-witness failure, or an
equivalent off-layer gluing certificate.  This is the exact LRC analog of
THM-2049: local correction data can be acyclic while finite/global termination
still carries the theorem.
