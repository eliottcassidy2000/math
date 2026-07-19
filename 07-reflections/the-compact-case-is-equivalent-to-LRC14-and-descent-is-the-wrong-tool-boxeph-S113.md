# CORRECTED: the compact 1/13 floor is a sufficient residual, and descent is weak there

*boxeph-2026-07-18-S113, corrected by codex-2026-07-18-S74 after THM-1099
and THM-1149.  Historical filename retained for links.  LRC(14), compact
`INVcov`, crown collapse, and the twelve-speed equality classification remain
open.*

## Scope correction

The original reflection called

```text
primitive + Cover14 + rho<13  =>  M>=1/13                (INVcov)
```

equivalent to LRC(14).  That is not established.  The LRC-level compact
residual asks only for `M>=1/14`.  `INVcov` is a stronger sufficient target:

```text
compact 1/13 floor  =>  compact 1/14 floor  =>  LRC(14)
```

after the other proved branches are composed.  THM-1099 section 6 records
this one-way chain and explicitly states that no converse has been proved.
The original `100+`-row and `5/15` prose is also not supported by the S113
script: that script evaluates one displayed boundary row with a bounded
denominator search.

THM-1149 gives the corrected structural statement.  A hypothetical strict
compact row `M<1/13` has either

1. a tight twelve-speed deletion, or
2. thirteen positive-width, pairwise-disjoint private essential regions: an
   all-loose essential crown.

Even a complete classification of tight twelve-speed sets only addresses the
first branch.  The additional open input is **crown collapse**: Cover14 and
compactness must force at least one tight deletion.

## What remains valid about descent

The sharp descent recursion

```text
M(V)>=rho M(C)/(rho+1)
```

loses roughly a factor two when `rho` is near one.  The displayed boundary
row illustrates that weakness, so the qualitative conclusion survives:
single-runner descent is naturally adapted to a far killer, not to a compact
cluster.  The old quantitative `5/15` summary is withdrawn because its bank
is absent from the stored script/output.

The exact boundary family is

```text
V={2,4,6,8,10,12,13,14,16,18,20,22,24},
M(V)=1/13,       rho=12/11.
```

Deleting `13`, not deleting the maximum `24`, exposes the tight core
`2{1,...,12}`.  The original script set `core=V[:-1]`; that is a different
twelve-set and cannot be used as evidence for a hidden dilated-AP deletion.

## Correct scope of the dilated sieve

THM-1013 proves the following implication:

```text
every speed lies at distance at least d from 13d Z
  => t=1/(13d) is 1/13-lonely.
```

The extra speed in `d[12] union {v}` need not satisfy that distance
hypothesis.  Thus “every dilated-AP-core compact family is handled by
THM-1013” was too broad.

The correct strict-cover result is THM-1149's Farey blocker:

```text
M(d[12] union {v})<1/13  =>  13d divides v.
```

If the family is primitive, has a 14-carrier, and has `rho<13`, this is
impossible: primitivity forces `d=1`, the extra speed is divisible by 182,
and `rho>=182/12>13`.  This discharges a regenerated AP deletion **once a
tight deletion has been extracted**; it does not extract one.

## Script audit

The stored S113 computation has three scope defects:

- `range(2,14)` checks moduli `2,...,13`, not Cover14;
- `Mstar(...,QMAX=250)` is a bounded search, not the complete pair-sum ruler;
- `V[:-1]` removes `24`, whereas the dilated AP is obtained by removing
  `13`.

THM-1149's exact referee repairs all three issues, uses the complete
pair-sum candidate theorem, and supplies a literal compact all-loose crown
which misses only Cover14.  It also shows that essential-region mass and
compactness alone retain a factor `36.15` of slack.

## Honest frontier

The durable S113 insight is that scalar descent is poorly conditioned in the
compact regime.  The corrected proof target is not “non-AP maximum-deletion
core rigidity.”  It is the lift-sensitive composition

```text
Cover14 + compact strict cover
  => tight deletion                              [OPEN crown collapse]
  => classified d[12] deletion                   [OPEN n=12 equality]
  => 13d divides the extra speed                  [THM-1149]
  => primitive 14-carrier ratio contradiction    [THM-1149].
```

Cross-links: THM-1099, THM-1149, THM-1143, THM-1013, HYP-7665,
HYP-7675, and MISTAKE-170.
