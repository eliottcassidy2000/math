# Incoming Session321 changes the retained arithmetic target

Fetched main tip: `0eab57c432`, 2026-11-13. These sources arrived after this
packet's initial derivations at baseline `68d0d86bef`. They were read before
publication; the correction and overlaps below are part of the final result.
Their exact source blobs are pinned in `SOURCE_PINS.json`.

## The shared metagraph presentation

Incoming [WEDGE_WIRING.md](../session321_selective_returns_and_fano_wiring_20261113/WEDGE_WIRING.md)
independently gives the same cut-difference integration and exact spanning-path
guard as our [native note](native/RESULT.md), including the same complementary
order-four witness. This is concurrent independent recovery, not priority
assigned to this packet. It additionally states the useful necessary score
identity `out_T(v)+out_U(v)=n-1+[v=start]-[v=end]` in the successful alignment.
Our additions are the uniform native jump weight, the explicit ordered-carrier
phase and dart-weight guard, and the independent complete small-order replay.
Neither presentation eliminates the need to enumerate permutation alignment
orbits at general order.

The earlier packet's full native merged matrices count both members of a
black reflection pair. Incoming Fano degrees use the fundamental black graph,
where that pair is counted once. For a chiral source the former is `2H/Aut`
and the latter `H/Aut`; these are compatible conventions. They must not be
mixed in a stochastic normalization or an incidence recurrence.

## The arithmetic correction is stronger than equal nullity

Incoming [BSD_PARITY.md](../session321_selective_returns_and_fano_wiring_20261113/BSD_PARITY.md)
shows that the whole labelled matrix `Phi(M)=M^2+M+I` is unchanged by
converse. The first bridge draft had proved only its transpose relation,
then incorrectly said that named vectors in its kernel require an orientation
choice. This false caveat was repaired before promotion; see
[CORRECTION.md](bridges/CORRECTION.md) and the dated entry in
[MISTAKES.md](../../../01-canon/MISTAKES.md).

For row-zero tournament M, `M^t=M+J+I`, `MJ=0`, `JM=(n-1)J` and `J^2=nJ`.
These identities make Phi symmetric and exactly converse-invariant. Its
diagonal and row sums are all one. A labelled wedge field therefore determines
Phi and its embedded binary kernel, as well as nullity. On this kernel every x has even coordinate sum, so changing orientation
changes the doubled Monsky coordinates `(x,(M+I)x)` by `(0,(n mod2)x)`.
The lift is therefore invariant at even n and needs the orientation at odd n;
the actual marked prime realization remains essential. The new target is larger; it
is not full BSD data.

Incoming [BSD_WIRING_STENCIL.md](../session321_selective_returns_and_fano_wiring_20261113/BSD_WIRING_STENCIL.md)
then supplies an exact transport constraint. If L is the binary Laplacian of
the marked Hamilton path, native transport gives

    M' = M + J + (n mod2)I + L,
    Phi(M') + Phi(M) = ML + LM + L^2 + L.

The last matrix is alternating and kills the all-ones vector. At interior
spine positions separated by at least three, it is the stated four-neighbour
stencil. It is a necessary constraint on a candidate alignment, not a
sufficient substitute for the wedge path guard. The independently written
[check_stencil.py](incoming/check_stencil.py) checks full converse equality
on every labelled tournament through6 and the complete identity on every
fixed-spine word through6, with a triangle's named kernel as a hostile to the
old caveat. All gates are explicit exceptions and both execution modes agree.

Incoming selective transitive ports preserve only the parity of their length
for this first-descent kernel dimension, while their Hamilton counts grow.
That strong equal-H ray and our all-order equal-full-Q ordinal family are
complementary failures: the former retains strong connectivity but not equal
Q, while the latter retains equal Q and scores but uses nonstrong suffixes.
Neither should silently inherit the other's stronger qualifier.

## Three Fano carriers, three distinct retained maps

1. The recovered 63 minimum Paley-7 blockers map by a type-aware decoder to
   nonzero `F8 direct-sum F8`. This is equivariant, but native deleted-arc XOR
   does not become vector addition without transporting the operation.
2. Our [thirty labelled Fano systems](bridges/FANO_ORBITS.md) carry the
   exponent G21 action. Two are fixed and the other28 are four natural
   seven-point orbits. Their exact map retains a type as well as a stabilizer
   centre; the centre alone is four-to-one.
3. Incoming [FANO_METAGRAPH.md](../session321_selective_returns_and_fano_wiring_20261113/FANO_METAGRAPH.md)
   studies 128 independent internal triangle twists of one 21-vertex
   tournament, yielding12 ordinary and6 merged classes. The seven canonical
   triangles have a Paley majority quotient, but their single dissent arcs
   retain Frobenius ports. The quotient converse conflicts with those forced
   ports, explaining why it does not lift. This is a different 128-state
   carrier from the thirty Fano line systems.

The incoming 8-to-1 majority quotient refines the earlier packet's negative
claim that evaluation on F8 is not a tournament quotient: majority is indeed
well defined here, but it is not uniform orientation across a pair of fibres.
The dissent port is the missing coordinate. The family has no internal native
edges because every source is regular and a native neighbour has different
degrees at its two endpoints. This finite independent set does not refute
all-order black-core connectivity through vertices outside the family.

The independent [Fano incoming audit](incoming/FANO_REVIEW.md) retains the
forced-map proof and replays all128 twist orientations, their canonical
partition, majority and mask actions. The large order21 Hamilton DP values
are inherited finite computations and are not claimed as independently
rerun in this packet.
