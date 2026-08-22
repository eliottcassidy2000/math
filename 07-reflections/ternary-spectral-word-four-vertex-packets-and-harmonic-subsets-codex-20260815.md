# Two four-vertex packets produce a ternary order-22 spectral word

**Research reflection / synthesis, not a truth source.**
The exact recurrence package is routed to audit-pending THM-3484.  The graph
determinant interpretation uses proved THM-3482.

## The compression mechanism

The private-support graph has two `K4` packets and one forced bridge:

```text
edges:    6 + 6 + 1 = 13,
tree degrees: 3 + 3 + 1 = 7.
```

The private-count gradient turns this into a degree-seven determinant on each
of the three `k mod 3` lanes.  A generic period-three degree-seven word needs
order `24`, because each of the three cubic Fourier colours may have degree
seven.  Here all lanes begin with the same term

```text
-16384 k^7.
```

The nontrivial cubic colours annihilate that shared top term.  Their degrees
drop to six, while the invariant colour remains degree seven.  The exact
dimension count is therefore

```text
8 + 7 + 7 = 22.
```

This is not a fitted recurrence.  The minimal characteristic polynomial is

```text
(x-1)^8 (x^2+x+1)^7,
```

and the first `22x22` Hankel determinant is nonzero.

## What the tournament language preserves

Each four-vertex packet has six pairwise coactivities.  Calling it a
bidirected tournament is useful only at that level: it preserves which pairs
can meet.  The recurrence needs more data:

- the owner-order orientation;
- the private-count potential;
- the two sixteen-tree sums;
- the ternary residue state.

Thus “tournament of size four” explains the six edge coordinates, while the
tree polynomial explains why only three edge factors survive in each spanning
tree.  The order-22 recurrence belongs to the weighted packet, not to the
unweighted tournament quotient.

## Subsets of the harmonic series

The three lanes partition the natural numbers into periodic subsets.  Each
therefore carries exactly one third of the logarithmic harmonic mass:

```text
sum_(k<=N, k==r mod3) 1/k = (1/3) log N + O(1).
```

This makes the spectral word a precise three-colouring of a harmonic
subseries.  The same phrase is dangerous for a general subset: every subset
does select a harmonic subseries, but that subseries may converge, diverge
irregularly, or lack a logarithmic coefficient.  Here periodicity—not merely
sethood—creates the `1/3` invariant.

## New transfer questions

The useful transplant pattern is now explicit:

```text
source:       a periodic family of packet weights,
representation: cubic Fourier colours,
invariant:    degree of each colour,
operation:    interlacing the three lanes,
target:       minimal linear recurrence,
lost data:    individual sheet positions and absolute graph flux.
```

Two questions are worth carrying back to the larger repo.

1. Does the THM-2334 relation-current transform have several residue colours
   with a common leading asymptotic, so that its apparent recurrence or Prony
   rank drops by the same Fourier-cancellation mechanism?
2. In the factorial/Jacobian lanes, do degree-graded monoid families share a
   leading term across residue classes, making their interlaced word have a
   smaller annihilator than the naive period-times-degree bound?

The present result does not answer either question, but it supplies the exact
test: Fourier-decompose the residue polynomials first, then count degrees by
colour.  A common top coefficient is not cosmetic; it removes one Jordan rung
at every nontrivial character.
