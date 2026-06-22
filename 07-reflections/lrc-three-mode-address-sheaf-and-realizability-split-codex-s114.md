---
source: codex-2026-06-22-S114
status: reflection / proof-order correction
tags: [lrc14, half-tiling, mobius, legendre, eisenstein, exact-period, finite-atlas, weyl]
---

# The three modes are an address sheaf, and the finite route splits from the analytic route

This is the HYP-2902 refinement of HYP-2901's lcm denominator-wall guardrail.

The useful correction is that the three recurrences should not be read as
three competing scalar identities.  They are local charts on the same
half-tiling address object.

Mobius is always present as the full inclusion-exclusion skeleton.  Eisenstein
is the even/pronic half-mode.  Legendre is the odd/square half-mode.  The same
letters `A..G` mean different subtournament sizes in each chart.

For Legendre odd size `N`, the corrected slots are:

```text
A,B : N-1
C,D : N-2
E,F : N-3
G   : N-4
```

The corners are `A,D,B`; the edges are `A+D-E`, `B+D-F`, `A+B-C`; and the
center is `A+B+D-C-E-F+G`.  In scalar cardinality, `C-D` cancels because both
are size `N-2`.  Geometrically they are not interchangeable: `D` is a corner,
while `C` is the subtraction on the `AB` edge.  This is exactly the kind of
label that LRC proof attempts keep rediscovering after a scalar quotient fails.

At LRC14, the top half-mode is even:

```text
h(14)=42=7*6.
```

The Eisenstein recurrence samples children of sizes `13` and `12`.  The apex
`7` is the pronic fold parameter, not one of those child sizes.  That means
the operator composition has two size coordinates:

```text
N=14 -> local recurrence children 13,12
N=14 -> fold/apex coordinate 7
```

The Legendre chart then applies on the exposed odd coordinate, with its own
local sizes.  This is the cleanest way to make the owner's "Mobius always,
Legendre odd, Eisenstein even" correction precise without flattening the
chart structure.

The same anti-scalarization rule explains the finite-basis failure.  For

```text
S_X={1,...,11,13,lcm(2,...,X)}
```

all denominators `D<=X` are killed by the last speed.  This proves the needed
unbounded-denominator obstruction.  But the first witness is not generally
`nextprime(X)`: exact scans give `X=14 -> D=41`, `X=24 -> D=53`, `X=41 -> D=53`.
The first opening is controlled by the AP-core safe-denominator ladder after
the divisor-loaded tail has killed the small chart.

So the proof split should now be stated carefully:

```text
Node 2: finite/algebraic AP-hull or three-gap majorization, with labels.
Node 3: analytic/exact-period floor after divisor-killed packets are removed.
```

A fixed finite residue basis cannot be the proof.  A finite atlas can still
handle coherent low packets, but an analytic Weyl/L2/equidistribution input is
irreducible for divisor-loaded tails.  The three-mode address sheaf is the
bridge: it says which labels have to survive long enough for the finite and
analytic branches to meet.
