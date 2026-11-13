# Reading nodes, line color and A000568 without conflating the weights

Use forward frame labels 0,...,n−1. The n−1 path arcs i→i+1
are frozen; the binom(n−1,2) other arcs are frame bits, indexed
lexicographically by i<j, j−i≥2. A bit 1 denotes i→j. Call a
framed instance a **presentation**, the orbit under free-bit
complement a **line instance**, and an ordinary isomorphism
class folded with its converse a **merged node**. Reflection
(i,j)↦(n−1−j,n−1−i) is the presentation of converse; the line
is blue precisely when its two presentations are individually
fixed by reflection. These definitions agree with THM-643/796.

Three distinct weights answer different questions:

1. **Presentation mass.** On an ordinary class [T] it is
   H(T)/|Aut(T)|: the automorphism action on ordered spanning
   paths is free. On a merged node the factor is doubled
   exactly when the class is not self-converse.
2. **Line incidence.** The blue or black half-edge count at a
   node is the number of presentations of that color in the
   node, for n≥3. A loop uses two presentations of that node.
   Different literal lines with the same node pair remain
   different instances. At n=1,2 the complement is trivial,
   so the presentation is a singleton orbit, not two half-edges.
3. **Class mass.** An ordinary class counts once in A000568
   regardless of its presentation mass. A merged node counts
   once if blue is incident and twice if not (including a blue
   loop as incidence).

The last statement is the inherited THM-643 color detector:
blue gives a self-converse tournament, and on the fibre of a
self-converse class the reflection involution fixes a
presentation because its size H(T)/|Aut(T)| is odd. Rédei gives
odd H(T); an involutory tournament automorphism would reverse
some two-cycle of names, which is impossible, so Cauchy's
theorem gives odd automorphism order. Thus, once the colored
merged graph is constructed,

    A000568(n) = 2 · (# merged nodes) − (# blue-incident nodes).

This is not an asymptotic approximation or a shortcut for
constructing the quotient. The associated unbiasing from raw
presentations is Σ_t |Aut(T_t)|/H(T_t). Those weights are
exact, but computing them need not be faster than direct
canonical isomorphism classification. It would also be
incorrect to count every blue line as a different self-
converse class.

## End-parent classes are not sufficient data

Two end deletions share a contiguous central overlap; their
actual frame bits must agree there. There are then exactly two
upper presentations, distinguished by the spanning corner arc.
Both lifts are blue iff the ordered end cards are mutual
reflections (then the corner is mirror-fixed); otherwise both
are black. In particular color does not follow from the two
merged end nodes.

Here are explicit hostile controls (rows are outgoing bitmasks).
At n=4 take the free masks 5 and 7. The rows are respectively
(6,12,8,1) and (14,12,8,0). In both, both end cards are transitive
3-tournaments and both upper frames are blue. The second
upper tournament is transitive; the first has a directed
triangle, so their merged nodes differ. At n≤3 both end cards
are trivial; already n=3 shows the two possible corner choices
can give the transitive tournament or the triangle. Order 4 is
the first case here with nonempty lower free-bit sets.

At n=5 take masks 0 and 1. Their rows are (2,4,9,19,7) and
(6,4,8,19,7). The right end card is unchanged and the left end
cards are isomorphic after allowing converse; permutation
canonicalization gives key 4 for all four lower cards (the key
is the minimum of the ordinary forward-pair mask over all
vertex orderings and their bitwise converses). The upper
canonical keys are 10 and 8 respectively. Upper mask 0 is blue;
mask 1 is black since reflection sends chord (0,2) to (2,4).
This witnesses loss of even color at the merged-pair projection.
These are exact finite controls supplementing the inherited
THM-796 non-Markovian weighted node spans.

The useful invariant is therefore the **shortcut-phase card**,
not the quotient node or its weighted degree. The reconstruction
and sharp minimal-card theorems in PROOF.md identify how little
framed information is sufficient to lift the line system across
size. Recovering the class mergers, rather than the lines, still
requires actual class information (e.g. an isomorphism algorithm).
