# A single compatibility edge can create the entire recurrence head

**Status: SYNTHESIS COMPANION TO
[THM-3305](../01-canon/theorems/THM-3305-single-edge-visible-zero-mode-and-rank-two-resolvent-update.md).**
The exhaustive finite statement is now proved and independently audited in
that theorem.  The primary exact companion is
[`gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.py`](../04-computation/gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.py),
with frozen
[`output`](../05-knowledge/results/gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.out).
An independent reconstruction is frozen in a separate
[`audit`](../04-computation/gmc_single_edge_zero_mode_resolvent_update_independent_audit_20260803.py)
and
[`audit output`](../05-knowledge/results/gmc_single_edge_zero_mode_resolvent_update_independent_audit_20260803.out).
It extends the static relation-walk analysis of
[THM-3288](../01-canon/theorems/THM-3288-maximizing-witness-lifted-walk-rational-series.md)
and the
[selector half-transfer scout](selector-paired-half-transfer-observer-quotient-codex-20260803.md).
The relations are symmetric compatibility edges, not chronology or a
tournament, and no FC, GMC, LRC, or positivity conclusion follows.

## Inheritance and exact question

The closest proved mechanism is THM-3288's finite transfer matrix: the selected
eleven-edge tree has `23` active `(row,witness)` vertices, `40` directed
compatibility arrows, and an all-ones walk series of minimal order `14`.  The
canonical hostile is still the row path `2-17-16`, whose two edges require
disjoint middle fibres.  The corrected near miss is the selector-only quotient:
the sign bit gives parity but not an equitable two-state dynamics.  The
least-used sidecars are the active-start vector and the full-core zero-mode
head.

The new bounded question is deliberately local.  The full core contains eleven
row edges outside the selected tree.  Add each one separately, rebuild its
literal maximizing-witness relation, and ask:

1. did the active vertex set change;
2. what are the exact full and parity Krylov orders; and
3. is the all-ones vector orthogonal to the adjacency kernel?

This is an exhaustive `11`-case census inside the already frozen core, not a
random edge experiment.

## The unique visible one-edge zero mode

The selected-tree control has

```text
active vertices 23, arrows 40, order 14,
even/odd half-orders 7/7, kernel dimension 9, zero mass 0.
```

Among all eleven one-edge additions, exactly

```text
row edge (7,18),     witness relation Q-7 <-> Q+1                 (1)
```

creates a nonzero all-ones projection onto the adjacency kernel.  Both endpoint
vertices already exist, and `(1)` contributes exactly two directed matrix
entries.  Thus the new adjacency is literally a rank-two update with no change
in the active-start vector.  Its exact signature is

```text
active vertices 23, arrows 42, order 15,
even/odd half-orders 8/7, kernel dimension 9, zero mass 3/5.       (2)
```

The other ten one-edge additions have zero visible kernel mass.  Three of them
create one new active vertex:

```text
(10,16) creates (10,Q-8),
(11,13) creates (13,Q-8),
(16,21) creates (21,Q-8).                                       (3)
```

They therefore require the active-start sidecar even though none creates the
head in `(2)`.  At the opposite extreme, adding `(19,22)` creates no vertex and
reduces the observer order from `14` to `12`.  Hence exact scalar recurrence
order is neither monotone under edge addition nor determined by the numbers of
vertices and arrows.

## Closed form and the head vector

For the selected tree plus `(7,18)`, the walk sequence starts

```text
23, 42, 116, 252, 692, 1602, 4404, 10518,
28966, 70344, 194026, 475722, 1313706, 3239490, 8953208, ...      (4)
```

Its minimal prefix-valid characteristic polynomial is

```text
z (z^2-2) (z^4-5z^2+3)
  (z^8-14z^6+62z^4-94z^2+30).                                  (5)
```

The exact Hankel rank is `15`, while its rational tail denominator has degree
`14`:

```text
(1-2x^2)(1-5x^2+3x^4)
 (1-14x^2+62x^4-94x^6+30x^8).                                  (6)
```

The numerator is

```text
23 + 42x - 367x^2 - 630x^3 + 2235x^4 + 3576x^5
-6528x^6 - 9600x^7 + 9436x^8 + 12456x^9
-6284x^10 - 7032x^11 + 1668x^12 + 1368x^13 - 108x^14.           (7)
```

As in the full core, the order/denominator discrepancy is one explicit
zero-mode head.  If

```text
p(t)=(t-2)(t^2-5t+3)(t^4-14t^3+62t^2-94t+30),
h=p(A^2)1,
```

then

```text
A h=0,       1^T h=-108,       h^T h=19440,
p(0)=-180,   1^T pi_0(1)=3/5.                                  (8)
```

The support is unusually small: `h=108` at `(18,Q+1)` and `h=-36` at

```text
(18,Q+2), (18,Q+4), (18,Q+5),
(22,Q+1), (22,Q+2), (22,Q+5).                                  (9)
```

There is a literal incidence explanation.  The central vertex `(18,Q+1)` has
the two neighbours

```text
(7,Q-7), (19,Q-1).
```

The other three row-18 vertices in `(9)` see only `(19,Q-1)`, while the three
row-22 vertices see only `(7,Q-7)`.  Therefore their neighbour-incidence
vectors obey

```text
3 N(18,Q+1)
 = N(18,Q+2)+N(18,Q+4)+N(18,Q+5)
   +N(22,Q+1)+N(22,Q+2)+N(22,Q+5).                     (10)
```

The primitive kernel coefficients in `(10)` sum to `-3`, so this is not an
all-ones-invisible difference of twins.  Multiplying them by `36` gives the
head in `(8)--(9)`.  The new edge is load-bearing because it supplies the
second central neighbour that completes `(10)`.

The degree-seven even recurrence fails only at its first possible index, by
`-108`, and holds thereafter.  Thus one compatibility edge creates the same
kind of one-step boundary debt that the entire full core displays, but with an
exactly localized support.

## A `2 x 2` update computes the new infinite sequence

Let `B` be the selected-tree adjacency, let `U` select the two vertices in
`(1)`, and put

```text
C=[0 1; 1 0],       A=B+U C U^T,
R=(I-xB)^(-1),      s=U^T R 1,       K=U^T R U.                  (11)
```

The scalar Woodbury identity gives the whole change in the generating
function from a two-dimensional calculation:

```text
G_A(x)-G_B(x)=x s^T (C-xK)^(-1) s.                              (12)
```

The scout verifies `(12)` symbolically over `Q(x)`, not from a truncated
series.  Every entry of `s` and `K` has the common local denominator

```text
d_old=1-13x^2+54x^4-77x^6+24x^8,
```

and the update introduces

```text
d_new=1-14x^2+62x^4-94x^6+30x^8.                              (13)
```

The correction denominator is exactly `d_old*d_new`; the two common factors
in `(6)` pass through unchanged.  This is the promised efficient closed-form
operation: a node-preserving relation edge can be compiled by a `2 x 2`
resolvent update while retaining the head.  The node-creating cases in `(3)`
instead require a bordered resolvent together with the changed start vector.

## Synthesis and scope

This exact example sharpens the boundary lesson shared with THM-3302's Lucas
doubling.  A short scalar update can alter only a small local factor and still
create a new boundary atom.  The safe compilation payload is therefore

```text
(rational observer, active-start vector, local resolvent, zero-mode head),  (14)
```

not the denominator alone.  Equation `(14)` is a procedural bridge, not an
identification of the graph and Lucas constructions.

The next finite test is clear: process node-preserving core edges by repeated
low-rank updates, process the three node births by bordered Schur complements,
and compare the resulting local secular factors with the final full-core head
`87/113`.  This would give an exact incremental compiler for these static walk
sequences.  It would still count compatibility walks, not physical response
time.

Reproduce with

```text
python3 04-computation/gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.py
python3 -O 04-computation/gmc_single_edge_zero_mode_resolvent_update_partial_scout_20260803.py
```

Both executions byte-match the stored output.  LF-normalized SHA-256 hashes:

```text
script 5704ee8d6641a12e50d6a6179ccf58f911ba782cfe369ac70a32db51887565b2
output b52db8c18d3182fb79ed5a230d9f75775155f438c7f32e974ab1ff3fcb15e556.
```

The independent audit imports or executes neither the primary scout nor the
selector half-transfer scout.  It rebuilds all eleven additions from the
pinned THM-3287 relation constructor and independently recovers `(1)--(13)`.
Its hashes are

```text
audit        f7f7fdfb32e26eff7d64d6e654916111d3b227122eaf3d69ba7b7501dd4cdb43
audit output 2fca320ccbf92fcf47a4a5eda9ed7c22c07a02843599c50eef4a29bdffd63acf.
```
