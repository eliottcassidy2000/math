# Berggren wall phases and Fibonacci projective phases share one signed C4 class

**Status: VERIFIED-EXACT ELEMENTARY CANDIDATE; unnumbered.**  The exact
companion is
`04-computation/berggren_fibonacci_signed_c4_cospan_20260816.py`.  This note
identifies one lawful common quotient of two previously separate four-state
systems.  It does not identify their states before quotienting, transfer
their densities, or construct an LRC/JC current.

## 1. The two four-cycles

For the two descendant-angle walls of THM-3334, the inverse Berggren map has
period-four itineraries

```text
alpha: B A B B,             beta: B C B B.              (1)
```

The affine CDF update reverses order on branch `B` and preserves it on
branches `A,C`.  Hence both cycles carry the sign word

```text
epsilon_wall=(-,+,-,-).                                  (2)
```

Independently, put

```text
x0=[1:0], x1=[0:1], x2=[1:1], x3=[1:2] in P^1(F3).
```

The Fibonacci matrix

```text
G=[0 1;1 1]
```

acts by the directed cycle

```text
x0 -> x1 -> x2 -> x3 -> x0.                              (3)
```

With the displayed representatives, the determinants of successive pairs
have signs

```text
epsilon_Fib=(+,-,+,+).                                   (4)
```

Both products are `-1`.

Here `(3)` is the **raw Fibonacci recurrence modulo three**.  At indices
where consecutive Fibonacci numbers are both odd, THM-3339 applies a
content-two parity normalization before interpreting the pair as a primitive
Berggren node.  That normalization is not slope-preserving (MISTAKE-418).
The cospan below deliberately retains only the raw projective recurrence; it
does not assert that `(3)` is the canonical normalized angle-state path.

## 2. The exact switching equivalence

Changing the representative at a projective vertex by `-1` multiplies the
two incident determinant signs by `-1`.  This is the usual vertex-sign gauge
on edge cochains.  Apply the alternating section

```text
s=(+,-,+,-).                                             (5)
```

Then, edge by edge around `(3)`,

```text
epsilon_Fib(i) s(i)s(i+1)=epsilon_wall(i).               (6)
```

In binary notation, the negative-edge indicators are

```text
c_wall=(1,0,1,1),          c_Fib=(0,1,0,0).              (7)
```

Their difference is the coboundary of the alternating vertex bit.  Both have
odd cycle sum and therefore represent the unique nonzero class

```text
[c_wall]=[c_Fib] !=0 in H^1(C4;F2)=F2.                  (8)
```

The exact companion exhausts all sixteen vertex sections.  Global sign is
the kernel, so the nonzero class has eight representatives; `(2)` and `(4)`
belong to that same orbit.

Thus the common object is precisely

```text
(pointed directed C4, nonzero F2 edge class).             (9)
```

It is stronger than an abstract cardinality-four coincidence and weaker than
an identification of the two dynamics.

## 3. Why four vertices naturally expose six pairs

Four vertices have six unordered pairs.  The successor law observes four of
them:

```text
01,12,23,30.
```

The two antipodal pairs

```text
02,13                                                   (10)
```

are absent from the recurrence.  Keeping only forward transitions gives an
oriented `C4` with two missing edges.  Keeping both forward and inverse
transitions gives four bidirected pairs and two missing pairs.  Both are
honest directed graphs of the kind suggested by “tournaments with missing
or both-way edges”; neither is a tournament.

There are four ways to orient `(10)` and complete the forward `C4` to a
tournament.  None is invariant under the phase rotation.  The half-turn fixes
each antipodal pair setwise while swapping its endpoints, which would force
an invariant arrow to point both ways.

On the Fibonacci side, a choice of all four projective representatives does
orient all six determinant pairs and hence supplies one tournament section.
Changing a representative switches every incident edge.  The result is a
switching class of eight sections, not a canonical tournament.  On the wall
side, no diagonal observable is present at all.  Consequently the size-six
object is best read as

```text
four observed transition pairs + two gauge-dependent/unobserved diagonals,
```

not as six new vertices.

## 4. What the cospan preserves and destroys

Point phase `i` in either wall orbit and map it to `x_i` in `(3)`.  After the
alternating section `(5)`, this gives the cospan

```text
alpha/beta wall phase --> signed pointed C4 <-- Fibonacci line mod 3. (11)
```

It preserves:

- the successor cycle;
- the vertex-sign gauge;
- the parity of orientation-reversing edges; and
- the nonzero cohomology class `(8)`.

It destroys:

- the quadratic wall values in `Q(sqrt(145))`;
- whether the positive wall branch is `A` or `C`;
- the affine CDF offsets `0` or `2/3` and contraction scale `1/3`;
- the exact Fibonacci fraction, height, and integral Farey frame;
- the q=15 owner/block interpretation of the projective cycle; and
- every full-tree branch word outside the four-state quotient.

The missing affine offsets are decisive.  The sign holonomy and scale explain
why the wall return denominator is `3^4+1`; the offsets distinguish `32/81`
from `50/81` and therefore produce the actual `16/41,9/41,16/41` split.
The Fibonacci projective cycle does not contain those offsets, so `(11)`
cannot transfer the angle densities.

## 5. Subsets of the naturals and harmonic series

The two four-state systems induce different subsets of the natural numbers.
The full ternary Berggren action on `P^1(F3)` partitions heap addresses into
four state fibres of density `1/4`; any union of `j` fibres has harmonic
coefficient `j/4`.  The angle-wall automaton instead partitions nodes into
three angle languages with coefficients `16/41,9/41,16/41`.

Equation `(11)` maps automaton **phases**, not tree nodes, so it does not map
one support subset to the other.  What survives is the edge cohomology class,
not the counting measure.  This makes the general harmonic-series slogan
precise: every subset of the natural numbers gives a harmonic subseries, but
only additional regularity or periodicity supplies a density and a controlled
logarithmic coefficient.  The carrier and enumeration cannot be omitted.

## 6. Relation to LRC and the D5-map program

This is a complete small example of how to write an explicit word-to-`H^1`
map:

```text
vertices = four reduction phases,
edges    = lawful successor transitions,
cochain  = order-reversal bit,
gauge    = vertex sign changes,
class    = odd cycle sum.                                (12)
```

It is useful as a template for the requested LRC word-current to JC-flux
map because every type in `(12)` is explicit.  It does not supply that map.
THM-3537 now supplies a third, genuinely Jacobian-side leg: the determinant
line of its quartic old-`L` inertia orbit has the same nonzero class.  In the
one-cut gauge `(1,0,0,0)`, the vertex switch `(0,0,0,1)` gives the wall word
`(1,0,1,1)`.  The accompanying quadratic orbit carries a second nonzero bit,
so the total old-`L` determinant class is `1+1=0`; this cancellation is data
that the isolated signed `C4` quotient forgets.  The explicit audit and full
loss ledger are in
`keller-inertia-and-berggren-wall-share-one-signed-c4-class-codex-20260816.md`.

The provisional THM-3534 response cospan lives over a large odd finite field,
requires quotienting an endpoint line, and lacks the physical same-copy
closure edge.  There is no coefficient-compatible identification with the
binary class `(8)`, and no characteristic-zero Jacobian flux follows.

Likewise RXTX's `4 x 4` block output has six off-diagonal positions, but no
map from those products to `(10)` or `(12)` is known.  The shared denominator
`41` remains an exact recurrence echo, not a tensor transfer.

## 7. Reproduction and scope

Run

```text
python -B 04-computation/berggren_fibonacci_signed_c4_cospan_20260816.py
python -B -O 04-computation/berggren_fibonacci_signed_c4_cospan_20260816.py
```

Normal and optimized runs match the stored transcript.  The LF-normalized
script hash and semantic hash are

```text
44f8eea13cac2b7c4e05c69ec15b67c828ef9ca4571389c9aa479733bf133e18
93a9f158e462686b4393163e46fcb211f607282eaaf3b8fd20a9f214156a7f9e.
```

No canonical tournament, angle-density transfer, Farey lift, physical LRC
current, D5 class, Jacobian flux, or matrix-multiplication statement is
proved.
