---
id: THM-383
name: lrc-half-turn-boundary-compactification-audit
status: PROVED (exact finite computation, bounded scope)
date: 2026-06-01
session: codex-2026-06-01-S513
depends_on:
  - THM-373
  - THM-374
---

# THM-383: Tie-wall compactification expands the bounded half-turn class menu

## Statement

In the bounded S512 circular-menu audit

```text
04-computation/lrc_iso_class_constraint_s512.py
```

the half-turn runner clock was sampled for primitive speed systems in the
following finite windows:

```text
total n=3, moving speeds <= 16
total n=4, moving speeds <= 10
total n=5, moving speeds <= 8
```

Open-cell samples use midpoints between consecutive event walls.  The
wall-compactified samples add the exact event-wall times, with the script's
fixed Hamiltonian tie path completing boundary tournaments.

The exact bounded counts are:

```text
total n   open phase classes   wall-compactified classes   A000568(n)
   3              2                         2                   2
   4              2                         4                   4
   5              4                        11                  12
```

The open-cell class details are:

```text
total n=3:
  H=1   score=(0,1,2)       c3=0
  H=3   score=(1,1,1)       c3=1

total n=4:
  H=1   score=(0,1,2,3)     c3=0
  H=5   score=(1,1,2,2)     c3=2

total n=5:
  H=1   score=(0,1,2,3,4)   c3=0
  H=9   score=(1,1,2,3,3)   c3=3
  H=11  score=(1,2,2,2,3)   c3=4
  H=15  score=(2,2,2,2,2)   c3=5
```

Thus, in this bounded audit, the open half-turn runner clock lands in a small
circular-menu subset of A000568 classes, while equality walls can add boundary
classes under the fixed tie path.

## Scope

This theorem records the exact bounded S512 audit.  It is not a full
classification of all circular half-turn tournaments at every `n`, and it is
not a claim that total `n=5` can reach all `12` A000568 classes after
wall-compactification.  The certified statement is the finite menu expansion
shown above.

## Proof

For each bounded speed system, the script forms the exact half-turn and LRC
event set.  For open cells it samples the rational midpoint of each adjacent
wall interval.  For the compactified menu it additionally samples all event
walls.  At each sampled time, it constructs the half-turn phase tournament,
using the fixed tie path on boundary comparisons, and canonicalizes the
tournament under relabeling.

The stored output reports the counts and class details in the statement.
Because the enumeration is finite and all event times are exact rational
values, this proves the bounded audit claim.

## Verification Record

Stored output:

```text
05-knowledge/results/lrc_iso_class_constraint_s512.out
```

The relevant blocks are:

```text
Circular half-turn menu probe for total n=3
Circular half-turn menu probe for total n=4
Circular half-turn menu probe for total n=5
```

The A000568 base values used in the comparison are computed by the script's
odd-partition Burnside formula and cross-checked by direct enumeration through
`n<=5`.

## Significance

This formalizes the boundary warning behind HYP-1982.  The open half-turn
clock is a restricted circular menu, matching the S24/S26 view of half-turn
tournaments.  LRC, however, often uses equality witnesses.  Once equality
walls are included and ties are completed, the class image can expand sharply:
from `2` to `4` classes at total `n=4`, and from `4` to `11` classes at total
`n=5` in the bounded audit.

Future A000568-style LRC proofs therefore need a boundary compactification or
an explicit rule excluding boundary-only witnesses.

## Related

- THM-373: runner phase clock wall decomposition.
- THM-374: transitive half-turn tournaments and semicircles.
- THM-382: threshold-decorated bounded fiber separation.
- HYP-1951: runner tournament clock circular menu.
- HYP-1982: threshold-decorated class fiber.
