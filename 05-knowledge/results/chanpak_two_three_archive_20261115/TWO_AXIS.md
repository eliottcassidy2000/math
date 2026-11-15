# Two genuinely independent physical slot switches

The [binary switch](CHAN_HISTORY.md) has equal H and Q₃ at its
two settings, and changes only Q₂ between 1,943,040 and 1,946,880.
A different pair from the [Cayley reader session](../chanpak_bott_cayley_periods_20261115/CAYLEY.md) does the opposite:
H=495 and Q₂=43,200 are equal, while Q₃ is 1,058,805 or 1,063,125.
These were previously studied as a collision of low-slot readers.
Putting the two constructions together gives an exact product of
two independently variable scales on *one physical family*.

Ordinally concatenate D_i (the binary switch setting), E_j
(the Cayley choice, with j=0 the lower Q₃), and a common transitive
tail of s vertices. Optional-slot profiles multiply. Thus the
resulting module has constant h=1155·495 and independent
coordinates

    Q₂ = 43200·Q₂(D_i)·2^s,
    Q₃ = 308859075·Q₃(E_j)·3^s.

Substitute this module into the
[three physical reader states](../chanpak_physical_routes_20261115/PORT_ROUTES.md).
Their readout is 4Q₃+9Q₂+c h, c∈{−3,−1,+1}.
This gives a literal 2×2×3 box of twelve strong tournaments,
with signed steps

    change i:  9·43200·3840 · 2^s
    change j:  4·308859075·4320 · 3^s
    change c by 2:  2·1155·495.

The three scopes are disjoint: the 11 internal D arcs, four
internal E arcs (codes 76→12 and 88→408 inverted for the increasing
Q₃ direction), and the two single exterior flips, one at a time.
No claim about the intermediate counts of a bundled edit is made.
All cross mixed differences between different coordinates vanish
exactly. The step scales are separated: the j-step exceeds the
i-step plus twice the c-step, and the i-step exceeds twice the
c-step. Thus all twelve count values are distinct even at s=0.
Their three step magnitudes have gcd 1350 for every s≥0, by
factorization (the powers 2 and 3 already exceed their
valuations in the constant step). Even this expanded box does
not provide a coprime pair of return shifts.

With m copies of D and k copies of E, choosing i of the upper-D
and j of the upper-E versions gives the product of two exact
sweeps, not just a union of one-dimensional experiments. All D
choices change Q₂ but not Q₃; all E choices change Q₃ but not Q₂.
The common powers 2^s and 3^s are still driven by the same padding
clock, but the multiplicative settings inside each coordinate
now really are independent. Within a single coordinate the steps
remain geometric, so neither coordinate has a unit additive
return. Different order of the same number of switches gives no
new counts. The [growth cone](SLOT_ENTROPY.md) and
[frozen-grade sieve](FROZEN_GRADES.md) explain two separate obstacles
to turning this box into the prescribed exact route.

The new idea is taking the product of two *different kernel
relations*. An actual module can give equality of a selected
set of readouts without retaining every other coordinate. Combining
relations with complementary defects creates independent degrees of
freedom. The remaining requirement is to keep a small absolute
step repeatable after the move, not merely usable in a finite
preallocated box.
