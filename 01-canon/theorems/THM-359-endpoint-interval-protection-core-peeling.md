---
id: THM-359
name: endpoint-interval-protection-core-peeling
status: PROVED
date: 2026-05-30
session: codex-2026-05-30-S362
depends_on:
  - THM-357
  - HYP-1811
---

# THM-359: Endpoint-Interval Protection Core Peeling

## Statement

Let `E` be a finite set of endpoints and `I` a finite set of intervals.  Suppose
we have two incidence maps:

```text
B : I -> P(E)      boundary endpoints of an interval
P : E -> P(I)      intervals that strictly protect an endpoint
```

Call a pair `(E', I')`, with `E' subset E` and `I' subset I`, a
**protection core** when:

1. every endpoint in `E'` is protected by at least one interval in `I'`; and
2. every interval in `I'` has all of its boundary endpoints in `E'`.

That is,

```text
for all e in E',  P(e) cap I' != empty,
for all i in I',  B(i) subset E'.
```

Start from `(E_0, I_0)=(E,I)` and repeatedly remove:

```text
D_m = { e in E_m : P(e) cap I_m = empty },
E_{m+1} = E_m \ D_m,
I_{m+1} = I_m \ { i in I_m : B(i) cap D_m != empty }.
```

When `D_m` is empty, stop.  The terminal pair is the largest protection core:
it contains every protection core.  In particular, if the terminal pair is
empty, then no nonempty protection core exists.

## Proof

The construction is monotone decreasing in the finite product poset

```text
P(E) x P(I).
```

Therefore it terminates.

Let `(C_E, C_I)` be any protection core.  We prove by induction that

```text
C_E subset E_m and C_I subset I_m
```

for every stage `m`.  This is true at `m=0`.  If it holds at stage `m`, then no
endpoint of `C_E` lies in `D_m`, because every `e in C_E` has some protector in
`C_I`, and `C_I subset I_m`.  Hence `C_E subset E_{m+1}`.

Likewise, no interval of `C_I` is removed at stage `m`, because every
`i in C_I` has `B(i) subset C_E`, and `C_E` is disjoint from `D_m`.  Hence
`C_I subset I_{m+1}`.  This completes the induction.

At termination, `D_m` is empty, so every remaining endpoint has a protector
among the remaining intervals.  Any remaining interval has never had a boundary
endpoint removed, so its boundary is contained in the remaining endpoint set.
Thus the terminal pair is itself a protection core.

Since every core survives every peeling step, and the terminal pair is a core,
it is the largest core.

## LRC Specialization

For a Lonely Runner speed set `V`, take `I` to be the pulled-back forbidden
intervals and `E` to be their distinct rational endpoints in

```text
Q(V) = (k+1)lcm(V).
```

The map `B` sends an interval to its two endpoints.  The map `P` sends an
endpoint to all open intervals that strictly contain it.  Empty terminal core
is therefore a finite certificate that no all-protected endpoint/interval
subsystem survives this peeling notion.

## Related

- THM-357: Lonely Runner endpoint-protection trichotomy.
- THM-358: initial-segment unit skeleton.
- HYP-1811: no all-protected endpoint core.
- HYP-1813: Bohr-boundary descent.
- `04-computation/lonely_runner_bohr_descent_s362.py`.
