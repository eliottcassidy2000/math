# Spanning trees, physical growth, and exact Hamilton-count flow

2026-11-13. **PROVED scoped / FINITE-EXACT controls. Chan–Pak and ordinary
Collatz remain OPEN.** The explicitly enlarged Collatz reachability problem
below is solved. The user clarified that every genuine reverse Collatz
edge should also be permitted; the earlier extra odd inverse is retained.

## A complete descent for the enlarged graph

Let B contain ordinary Collatz, every genuine reverse edge, and the extra
move n→(n−1)/3 at odd n>1 with n≡1 mod3. The last move has an even quotient
and is not a genuine Collatz inverse. Define, for n>1,

| Input | Smaller macro-parent D(n) | Native edges required |
|---|---:|---:|
| even n | n/2 |1|
| odd n≡0 mod3 | n/3 |7|
| odd n≡1 mod3 | (n−1)/3 |1|
| odd n≡2 mod3 | (2n−1)/3 |2|

The decisive multiple-of-three route is

    n→3n+1→6n+2→12n+4→4n+1→4n/3→2n/3→n/3.

Each edge has its literal parity/divisibility guard. Since0<D(n)≤2n/3,
every positive integer reaches1. Its route has peak≤12n+4 and length at
most7 ceil(log_(3/2)n). The macro-parent tree on[1,N] grows by attaching
one leaf at each new integer. Expanding its edges requires the stated
larger label budget; the native prefix[1,N] need not be connected.
There is also a native-edge spanning arborescence, obtained by choosing
successors that decrease the shortest distance to1.

Genuine reverse edges **without** the extra odd-to-even inverse preserve
the ordinary basin, so their reachability question is exactly ordinary
Collatz. The incoming literal-relaxation peak obstruction uses a smaller
move alphabet and does not apply to B. See [the proof](collatz/PROOF.md),
[old-work recovery](collatz/HISTORY.md), and [incoming integration](INCOMING.md).

## Exact growth of tournament occurrences

For an ordered word W of r nontrivial strong components, with total order
m, singleton gaps form the lattice N^(r+1). At order n there are exactly
binom(n−m+r,r) ordinary classes of that word. The fixed-H rational series
and eventual polynomial are inherited. The [new growth treatment](multiplicity/PROOF.md)
adds a canonical collision-free parent and permanent named ports.

The canonical padding forest becomes asymptotically unary: its branching
fraction at fixed H=h is r_h/n+O_h(n^−2). Counts with k ordered permanent
ports have degree r_h+k. However indegree-one unit-ready rooted counts
eventually stabilize at a constant. At H9 and order14 there are65 ordinary
classes,730 rooted classes and7870 ordered-two-port classes, but only one
unit-ready rooted class. A large fibre need not supply a new useful interface.

For complete labelled rows there is an exact birth flow. Let B_T(y) count
the2^n newborn cuts of a source T by their resulting Hamilton total. Then

    L_(n+1)(k)=[y^k] sum_[T],|T|=n (n!/|Aut(T)|) B_T(y).

The source coefficient needs the actual endpoint/two-path response.
Two order-five H3 sources, with the same core word and automorphism order,
have different next-count support. Relocating one named singleton across
C3 preserves H3, yet the same named newborn cut gives17 before preparation
and29 afterward. See [the exact flow and literal preparation](multiplicity/BIRTH_FLOW.md).

The [general relocation law](multiplicity/NEUTRAL_GAP.md) works for every
nonempty core and every binary tail, giving an exact four-coordinate
response and a two-phase arithmetic inverse for the count difference.
When the newborn beats the moved singleton, that singleton has indegree
one and defect S_J, the number of core paths starting in its selected
entry set. Thus S_J=1 supplies a unit-ready output port for every tail.
In the literal example the H29 output is ready at the moved singleton.
This handoff is essential: the newborn itself is not ready there.
The ready output family is the inherited selected-entry funnel language;
the operational gain is its explicit same-name, same-star preparation and
port handoff, not an additional arithmetic coverage family.

The preparation preserves even the full polynomial of Hamilton counts on
every named induced subset. That polynomial still cannot predict the next
fixed-star attachment; the missing information is the endpoint response.

## The tree/path matrix bridge retains the right information

Ordinary A^k counts walks, with repeated vertices allowed. Hamilton paths
on n vertices have n−1 edges; Hamilton cycles have n edges. Neither the
exponent nor a diagonal/off-diagonal test by itself enforces occupancy.
For instance TT_n has one Hamilton path but A^n=0.

With commuting x_v satisfying x_v²=0 and X=diag(x_v),

    H(G)=[x_V] 1ᵀX(AX)^(n−1)1.

The known determinant quotient explains tournament parity from
A+Aᵀ=J−I. Marking each edge of an out-arborescence by its **tail** variable
and discarding repeated powers instead extracts its Hamilton-path layer.
Marking heads would keep every tree and would not detect branching.
This gives a precise tree/path correspondence, with exponential subset
information retained. [The matrix proof](matrix/PROOF.md) includes endpoint
cofactors, Boolean coefficient extraction, cycles, and the inherited exact
directed-cut birth law. Irving–Omar's existing parity treatment is credited.

Outside tournaments, inherited one-missing-pair witnesses already attain
7 and21. Completing their missing pair changes7 to9 in either orientation,
and21 to23 or41. Thus the tournament condition constrains which repairs
are available; ordinary digraph walk counts do not explain the exclusions.

## What still separates this from a Chan–Pak proof

Let S_n be the Hamilton counts attained at order n. Source padding gives
S_n⊆S_(n+1). Therefore the assertion that the missing odd values tend
pointwise to{7,21} means precisely

    for every odd h outside{7,21}, some n has N_n(h)>0.

That is the remaining conjecture, not an input justified by growing total
mass or fixed-count multiplicity. For each finite bound L, one must supply
actual witnesses for every permitted h≤L, with a proof for arbitrary L.

The macro-parent gives an [arithmetic tree of all desired odd targets](TARGET_TREE.md)
using four guarded affine rules and five finite bypasses. Its tournament
lift remains CONDITIONAL. The tree contains13→27, impossible under any
pure-birth lift at every order. It also hits a demonstrated all-order
unit-ready hole at159252697. Thus a complete lift needs certified changes
of old arcs, cores or interfaces; neither a fair enumeration nor neutral
padding alone supplies it.

## Reproduction and provenance

[REVIEW.md](REVIEW.md) records exact universes, independent proof reviews,
commands and scope boundaries. The complete small inventories are inherited
and credited; this session does not claim a new tournament census or a
larger verified spectrum. [BOARD.md](BOARD.md) retains the live mechanisms;
[INCOMING.md](INCOMING.md) distinguishes concurrent results and graph models.
