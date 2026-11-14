# Exact profile multiplicities in the triangular component family

**Proved in this restricted class.** Point the added apex and suppose every
SCC of the old core is a singleton or a directed triangle. Distinct
pointed isomorphism classes can share the same entire two-path response.
The multiplicity is exactly a dimer-placement coefficient.

Use the inherited [contact decoder](../chanpak_f1_physical_carry_20261114/SUPPORT_DECODER.md)
and the ordinal laws in [ORDINAL_SELECTOR](ORDINAL_SELECTOR.md). For h=3^r,
f/h=2^w, the normalized coordinate P/h is a fixed-width bit word of width
w: each singleton contributes its contact bit, and each triangle uses
two bits giving the number c=0,1,2,3 of its contacted vertices. The
two normalized bases are f_i/h_i=2 and 4, and p_i=P_i/h_i=c.
For Q/h the local digit is the complementary contact number (1−c
or 3−c). Subtract the integer formed by reversing and complementing
the P bits. A mixed triangle contributes +2^a for c=1 and −2^a for
c=2, with a its starting bit position measured from the left. All
other blocks contribute 0. The nonzero positions are nonadjacent.
Such signed expansions are unique: at the lowest differing position,
presence on only one side contradicts parity, and opposite signs
on both sides contradict the next modulus 4 because their next
digits must both vanish. Thus h,f,P,Q recover every mixed triangle.

Normalize all seven profile coordinates by h. A homogeneous triangle
has the same normalized coordinates as two consecutive equal-bit
singletons: on contacts all zero or all one the normalized pair
total is 4 and the nonempty pair-end/start count is respectively 0
or 2. Direct starts/ends and the normalized single-restriction pair
counts also agree. The two-block laws retain equality in arbitrary
ordinal contexts. Consequently no residual variation remains at
fixed h,f,P,Q in this family, despite multi-class fibers.

Let m be the number of mixed triangles detected. Keep their slots
as barriers and split the remaining slots into maximal constant-bit
runs of lengths L. Each of the r−m remaining triangles occupies a
domino within one such run. In a length-L run there are
binom(L−j,j) ways to place j nonoverlapping dominoes, by subtracting
one from each successive left-end separation. Hence the precise
number of old marked isomorphism classes in this profile is

    [z^(r−m)] product over L of (sum_j binom(L−j,j) z^j).

The SCC order is canonical and all same-cardinality subsets on a
triangle are equivalent by its 3 rotations. Thus placements are
neither label counts nor merely a surrogate statistic. For example
the all-zero width-six word with r=2 gives binom(4,2)=6 distinct
marked classes, all with (h,f,P,Q,S,T,K)=(9,576,0,567,0,9,9).
An unpointed class after adjoining the apex may identify additional
objects by moving the apex, so the formula cannot be inserted
as an unpointed A000568 coefficient.
