# The self-converse spectrum genuinely excludes 63 (recovered)

The exclusion and mirrored-component reduction already appear in the
[Session 24 multiplicity study](../../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md).
This is an independent short replay plus a blue-source interpretation.

Let S be the set of all tournament Hamilton counts, and let C be
the counts of strong self-converse tournaments together with 1.
Then the exact self-converse spectrum is

    S_SC = {g² c : g in S, c in C}.

The strongly connected components of a tournament are totally ordered:
arcs between successive components all go forward. Each Hamilton path
concatenates Hamilton paths of the components, so H is multiplicative.
Converse reverses their order and takes each to its own converse.
A self-converse tournament therefore has mirrored pairs of components
and at most one center, which must itself be self-converse. Conversely
the ordinal stack X -> B -> X^op is self-converse for any center
B of the indicated type (including empty). This proves the formula.

In particular there is **no self-converse 63-path tournament at any
size**, though 63 is an actual tournament count. Indeed the only
possible square factors are 1 and 9, and dividing by 9 leaves the
forbidden 7. The center with count 63 must be strong. Busch's sharp
[minimum for strong tournaments](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
is already 75 at order 9 and increases thereafter. The verified
birth atlas rules out order <=7 (also obtained independently among
blue tilings by this replay), and the even self-converse mod-4 law
(archive THM-791) rules out order 8.

More fully, the stand-alone `holes.py` generates only blue fixed-frame
orientations through order 8, keeps the strong ones and computes
their Hamilton counts by subset DP. This enumerates the **complete**
potential centers below 75: any self-converse tournament has an
anti-involution (take a Sylow 2-subgroup in its automorphisms and
anti-automorphisms; the ordinary subgroup has odd order). Reverse
each directed Hamilton order and apply the anti-involution.
The action is an involution on an odd set, by Redei, so a fixed
Hamilton order realizes it as standard reversal. Thus every
self-converse center has a representative among the enumerated frames.
Only squared side factors 1,9,25 can occur below 75: the odd square
49 would require a count-7 side. This gives the exact all-order
self-converse omissions below 75:

    7, 19, 21, 23, 35, 39, 43, 47, 53, 55, 59, 63.

The ordinary spectrum in this range misses only 7 and 21 (archive
THM-4744). Thus restricting the *witness itself* to the self-converse
vertices would disprove a stronger, incorrect formulation. It does
not obstruct computation via a blue encoding or using blue sources
followed by a symmetry-breaking output.

## An exact two-arc exit from blue at every even size >=6

The prior [two-contact power construction](https://github.com/eliottcassidy2000/math/tree/395537b8edc50a89a42bb7742e61c20828620ab3/05-knowledge/results/chanpak_exponent_dimension_20261117/apex)
gives H=4^(k+1)+1 and H=4^(k+1)-1 at order 2k+4, k>=1.
Its upper tag word is 1^(2k+2)0. In the order starting at the apex
then following its transitive spine this is precisely a transitive
tournament with the extreme arc reversed. Thus its frame reversal
is an anti-involution. At the lower endpoint the even-order mod-4
obstruction certifies that the tournament is chiral. The two flips are
literal chords of the common Hamilton framing, not a native blue
metagraph line. At k=2 this gives a specific blue 65-path source
whose two-contact edit reaches one of the ordinary 63-path birth
classes. This exemplifies a viable weaker goal: choose a blue
seed and a *controlled exit*, rather than demanding the target itself
be blue. No theorem yet gives such a seed and exit at every target.
