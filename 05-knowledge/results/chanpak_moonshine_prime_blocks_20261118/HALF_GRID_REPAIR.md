# Half-grid reflection and a continuation-preserving repair at actual misses

**Status:** all-order replacement law for shared-root ordinal arms;
literal protected unit moves, a two-digit composition, and a scoped exact
comparison against the frozen decoder bank. Chan–Pak remains open.

## What the half-integer analogy actually transfers

The earlier [half-grid selector](../chanpak_half_grid_selector_20261118/README.md)
and the earlier midpoint-chamber work
used two integer-translation classes in (1/2)Z. This is an additive group,
not an ordinary ring; the transported product x⋆y=2xy stays in it.
The coordinate H/2 of an odd count is a proper half-integer and a +2
repair is a unit translation. Negation preserves both cosets and changes
sign. A reflection about an integer center preserves them too, whereas
translation by 1/2 exchanges them. Thus parity and sign are distinct
pieces of information. Negative signs here denote increments or
antisymmetric reads, not negative path counts.

For [Riemann's completed function](https://dlmf.nist.gov/25.4),
s=1/2+w turns its reflection into w→−w. An even function of w can
still have zeros away from the symmetry line. The analogy is to
use coordinates which separate information fixed by reflection from
information reversed by it, not to infer RH from reflection or
Chan–Pak from RH. The older
[quasi-Riemann/Collatz archive](../collatz_integers_as_irrationality_measures_20261014.md)
makes the analogous distinction between a bound strong enough for
a particular threshold and a bound that settles all scales.

For a tournament arm, taking converse exchanges the weighted starts
and ends, and leaves H fixed. Set M=(U+Ubar)/2, D=(U−Ubar)/2,
and use the same splitting for the parity-signed reads V,Vbar.
Converse becomes (M,D)→(M,−D). In an order exchange of two arms,
the cross term changes by
```
U_A Ubar_B − U_B Ubar_A = 2(D_A M_B − M_A D_B)
```
plus the analogous signed expression. This is an exact two-directional
read, not a claim about zeta zeros. To change H while preserving *all*
these coordinates, we need replacement arms, not merely converse.

## Protected replacement through arbitrary word contexts

An arm consists of a tournament with a marked root. Its signature is
(g,B,U,V,Ubar,Vbar), where g counts Hamiltonian paths of the deletion,
B counts ordered two-slot covers of it allowing empty slots,
U sums paths on nonempty deletion subsets whose starts are beaten
by the root, and Ubar uses ends beating the root. The weight of
each path is the Hamiltonian count of the complementary subset.
V and Vbar have the same summands with sign (−1)^subset_size.

In the [shared-root arm product](ALIGNED_ARMS.md), identify roots
and orient every deletion edge from the earlier arm X to the later
arm Y. The previous transfer reads
```
g' = g_X g_Y,    B' = B_X B_Y,
U' = B_Y U_X + g_X U_Y,    V' = g_X V_Y,
Ubar' = B_X Ubar_Y + g_Y Ubar_X,    Vbar' = g_Y Vbar_X,
H' = g_Y H_X + g_X H_Y − g_X g_Y + U_X Ubar_Y + V_X Vbar_Y.
```
The weighted reads follow by placing a selected path in the ordered
deletions and putting the other vertices in the complementary slot.
The signed sum of optional splits of a nonempty deletion is zero:
the two-cover convolution evaluated at (−1)+1=0 gives this by
f_Q(−1)=(−1)^|Q| f_Q(1), the
[odd-cycle expansion](https://arxiv.org/abs/2307.05569).
For H, count spanning packings by odd directed cycles and singletons.
A cycle through the root either stays in one arm or uses a nonempty
path in both, with even combined order away from the root. Its
weight two cancels the denominator two of the parity projector.
That gives the two cross terms; the two confined terms double-count
the packings with the root singleton, giving −g_X g_Y.

It follows immediately that replacing an arm by another with exactly
the same signature changes nothing the arm product reads except H.
A change ΔH in arm i is multiplied by the product of the g_j for
j≠i, regardless of the surrounding word. With several independent
replacements, these changes *add*, with no mixed interaction.
These are substitutions of explicit tournaments, not necessarily
single-arc edits or transitions between arbitrary representatives
of a scalar count.

Here are two fixed templates. Masks list outgoing arcs, indexed
from zero. The common signature of each pair can be checked by
subset DP.

| H | outgoing masks | root | common signature |
|---:|---|---:|---|
| 41 | (126,28,56,112,96,66,6) | 3 | (15,240,63,−9,99,−15) |
| 43 | (126,28,40,112,100,66,6) | 4 | (15,240,63,−9,99,−15) |
| 65 | (62,124,24,112,96,68,5) | 2 | (17,224,72,−6,135,−11) |
| 69 | (16,65,3,71,14,31,53) | 1 | (17,224,72,−6,135,−11) |

In particular the first pair supplies **±2g** with any exterior
arm word whose deletion count is g. This includes arbitrary
transitive deletions (g=1), so increasing the exterior size does
not spend the unit step. The last pair supplies ±4g. Putting both
pairs next to each other makes independent digits of size 17 and
30 **on the half-grid**, with literal H=10021,10055,10081,10115.
These four also have identical retained-root continuation
(e,t,d)=(165,15447,255). In general a common arm signature fixes
(t,d)=(g+U,g); if ordinary endpoint e differs, an ordinal sink
sets e=0, doubles t, and leaves H and d unchanged. All later words
of the [retained-root decoder](ANS_INVERSE.md) preserve the exact
signed difference.

## Place the repair against a fixed bank

Append a right arm whose deletion is a transitive tape and whose
incontacts encode I=Σ 2^j in increasing tape order. The endpoint
reads of this arm are Ubar=I and Vbar=−(I mod 2), by the binomial
sum over predecessors of each endpoint. Its standalone count
is 1+C(I), where
```
C(I) = sum 2^max(1,j−i−1) over i<j, bit i=0, bit j=1.
```
Choose an even number of interior vertices between the two
contacts of the unique possible root cycle and apply the weight
two from the odd-cycle expansion. This is also the inherited
binary carry checksum: C(2J+1)=C(J), C(2J)=C(J)+J+(J mod 2).
For the 41/43 pair the two *absolute* outputs are
```
H_-(I)=41+63I+15C(I)+9(I mod 2),
H_+(I)=H_-(I)+2.
```
The destination mod 3 first chooses the orientation (2 for minus,
1 for plus, and no match for 0), then the formula checks the
absolute target. The nonnegative checksum gives the exact search
bound I≤(N−41)/63; the same code decides failure in this lane.
Positive I is encoded without leading high zero bits, which can
only pad counts. The empty context can supply 41/43 separately.

Compare with the unchanged 11-source bank from
[the prequential test](AMM12592_SOURCE_BANK.md). Count a repair
only when exactly one endpoint was supplied by that frozen bank,
both lie in the tested block, and the absent endpoint is
constructed by the paired replacement. Results:

| Target block | Absent odds reached by a protected step | Of them prime |
|---|---:|---:|
| [4096,16384) | 37 | 10 |
| [16384,65536) | 121 | 41 |
| [65536,262144) | 552 | 152 |

Example I=892 gives **66349 → 66347**. The original bank supplies
the former but not the latter. The paired templates re-encode
the former in the correct interface first. The two physical
graphs have the identical (e,t,d)=(15,76107,15), confirmed by
direct subset counting, so **no reset is needed in this
example**. We do not claim that every old 66349 representative
admits this edit. This is a post-training structural repair of
the already frozen bank, not a second blind holdout.

## The local-to-global obligation

A shared-root word with pairwise coprime deletion counts g_i
weights local half-steps by P/g_i, where P=product g_i. If each
local palette contains all residues mod g_i, the Chinese
remainder theorem would select independent physical digits to
match a target mod P. Proof: the ith summand is the only one
nonzero mod g_i and P/g_i is a unit there. These actual binary
palettes do not supply complete residue sets. Even if enlarged
palettes do, the prescribed *absolute* equality requires a physical
carry or a bound on the remaining error, as with a local-to-global
representation approach to Goldbach. Congruences or aggregate
word capacity do not supply that step. Repeating the same g
instead of choosing pairwise coprime g_i creates a common
divisor in the steps and cannot simulate independent residues.

The new exact step is a durable repair, not a selector for every
target. There are at most 2X/63+O(1) totals through X in the
41/43 right-tape lane, and arbitrary-context replacements still
need a proof of simultaneous center and shift placement. The
stronger demand is to build palettes large enough to match the
target phase and independently move the absolute baseline, both
while retaining continuation.

Replay: `python3 reflection_replay.py` with the archived
`benchmark.py` and `fixtures.py` in the same packet.
