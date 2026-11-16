# A triangle budget for induced repair

Let H(T) count Hamiltonian paths (H(empty)=1) and τ(T) count directed
triangles. Use the Grinberg–Stanley expansion of H as the sum of 2 to
the number of cycles in a vertex-disjoint packing of directed odd cycles
(empty packing included). In particular H is induced-monotone, H≥1+2τ,
and adding x gives the exact difference
$$
r=(H(T+x)-H(T))/2
  =\sum_{C\text{ directed odd cycle through }x}H(T-V(C)).
$$
As usual cycles are oriented cyclic orderings, not counted with a marked
starting point. Write O=N⁺(x), I=N⁻(x) in T. If r>0 there is a triangle
through x: any return to x must contain some O-to-I arc. Let there be
t such crossing arcs aᵢ→bᵢ. Fix any one, a→b, and put q=H(T−{a,b}).

**Triangle-budget theorem.** If r>0, then
$$
                         r\ge\tau(T)+(q-1)/2.
$$
This is true for every rooted triangle, so one may maximize q.

For each aᵢ→bᵢ let dᵢ count the predecessors of aᵢ in O and eᵢ the
successors of bᵢ in I. There are at least binomial(dᵢ+eᵢ,2) distinct
rooted 5-cycles with precisely this one O-to-I transition: take two
predecessors in their arc order, two successors in their arc order,
or one of each, and concatenate through aᵢ→bᵢ. All mixed-side old
triangles have exactly one O-to-I arc and are bounded in number by
Σ(dᵢ+eᵢ)≤t+Σ binomial(dᵢ+eᵢ,2), using s≤1+binomial(s,2).
Internal triangles in O through a give *further* rooted 5-cycles:
u→v→a→u becomes x→u→v→a→b→x; u is not a predecessor of a.
Reversing the argument for internal triangles in I through b yields
further cycles with three I vertices instead of three O vertices.
These images are disjoint from every previously charged 5-cycle
and from one another. The remaining old triangles avoid a,b, so
number at most (q−1)/2. In the exact rooted sum, all the rooted
triangles and charged 5-cycles have weight at least one, and the
chosen root triangle has extra weight q−1. Subtracting the charged
categories proves the inequality. Keeping the triangle count of
the remainder exact gives the stronger clause
`r ≥ τ(T)+q−1−τ(T−{a,b})`. Thus forbid the conjunction
`x→a and b→x` whenever this threshold or q exceeds the target
budget r. These are two-literal constraints for a single port.

## Consequences for an absolute window

For any induced extension T⊂U with H(U)=H(T)+2r and r>0, insert the
extra vertices one at a time and consider the first strict increase.
It shows τ(T)≤r. More precisely, at each step the triangle increment
is at most its half-gap (and is zero at a neutral step). Thus
τ(U)≤τ(T)+r≤2r. Every vertex in a nonsingleton strong component
belongs to a directed triangle: otherwise its outneighbors cannot
reach its inneighbors within that strong component. All other
vertices are singleton strong components. Consequently the whole
*pair*, after removing vertices not in directed triangles of U
from both graphs, fits on at most 6r vertices, without changing
either count (ordinal-sum counts multiply). An empty trimmed source
means count one and can use the empty convention. In particular a
positive bounded additive repair cannot keep a large witness as an
induced subtournament. The rough numerical bound for the old count
is (3r)! by counting only the support of old triangles. This is an
all-orders result, not a restriction to the known ribbon recipes.

**A growing forbidden band.** In fact H(T)≤6^τ(T). Encode a
Hamilton path by the relative ordering of each of its directed
triangles. Call uRv when u→v is in no directed triangle. In
outneighborhood notation this is exactly N⁺(v)⊂N⁺(u), with
strict inclusion, so R is a strict partial order. Two permutations
with the same code agree on all incomparable pairs of R.
They are therefore related by swaps of adjacent R-comparable
vertices: bring the desired first vertex leftward, crossing
only comparable vertices, then induct.

Rewrite such adjacent swaps only toward the order prescribed
by R. This terminates by inversion count in any linear extension
of R. Disjoint rewrites commute; overlapping rewrites are a
three-term descending chain and both routes reach its ascending
ordering. Hence it is confluent (or induct on remaining inversion
count using this local diamond). A Hamiltonian path is already
irreducible: every adjacent pair follows its tournament arc.
Thus no two paths have the same triangle code. The 6^τ bound follows. If b is the number of components in
the graph of directed triangles joined when they share an arc,
the sharper bound is 6^b·3^(τ−b): introduce the triangles of
each component in spanning-tree order, and each subsequent
triangle has one pair already ordered, leaving at most three
options. Equivalently the code determines an acyclic orientation
of the graph of all triangle-supported pairs.

In particular any literal induced extension which increases H
must add at least 2⌈log₆ H⌉. The forbidden band below a prescribed
large target grows, even if we are free to add arbitrarily many
vertices. This is a quantitative reason to look for cancellation
between two reader choices, not a small positive debit from one.

For r=1 the old τ is at most one, so its count is 1 or 3. Both
1→3 and 3→5 occur. This classification also holds with arbitrarily
many added vertices. For single insertions the r=2 and r=3 lists
are respectively {1→5,5→9} and {3→9,5→11,9→15}.
Here is a reduction of the latter two lists to tiny local facts.
A strong SCC on s≥3 vertices has at least s−2 triangles: Moon's
vertex-pancyclicity supplies a strong induced (s−1)-vertex set
for s≥4, and the omitted vertex lies in an additional triangle.
The strong classes of orders 3 and 4 have (τ,H)=(1,3),(2,5).
At order 5 the τ=3 classes all have H=9 (the accompanying
standalone labelled check also gives the other τ fibres). Thus,
when τ≤3, possible old counts are 1,3,5,9,15,27; with τ≤2
they are 1,3,5,9 only.

If τ=r, the budget forces q=1 for every root triangle. For old
H27 and τ≤3 the old graph is three triangular SCCs, and two removed
vertices cannot kill all three. For old H15 and τ≤3 it is a
three-vertex and four-vertex SCC, and for old H9 and τ≤2 it is two
three-vertex SCCs. In either two-block case a root triangle must
meet both blocks. A mixed orientation from x to a strong block
would produce a root triangle fully in that block, so the
orientations from x to each block are uniform. A root triangle
meeting both forces all cross pairs to make root triangles, 12
or 9, exceeding the budget.

The remaining 3→7 and 1→7 options are excluded by the known H7
exclusion (or the strong order bound of Busch: only a single strong
factor is possible for prime 7, and its order is at most 4).
All lower-list positive cases have literal witnesses in check_core.py's
enumeration, independent of any saved isomorphism atlas. For a
general induced extension with total r=2 or 3, concatenate
smaller single-step moves, allowing neutral insertions. The only
additional chain for r=2 is 1→3→5 and for r=3 is 3→5→9; their
endpoints are already in the displayed lists. Thus all three
small-gap classifications hold with any number of insertions.
For N≥17 a tournament of count N cannot contain any induced
subtournament with counts N−2, N−4 or N−6.

## Import from the Boolean Pythagorean proof: a backbone test

[Heule, Kullmann and Marek, §§5–6](https://arxiv.org/html/1605.00723)
eliminate degree-one constraints iteratively and split the remaining
problem into cubes. The final 7825 vertex contributes seven triples,
or fourteen clauses; 6494 of 7825 numbers occur in triples. They
identify opposite forced colors at the new vertex from two triples
in a backbone of a restricted satisfiable instance, not simply from
the number of new clauses.

For any *partial* tournament, put a hyperedge on a triple if it
can be directed cyclically in some completion. Its components
are the finest fixed partition within which all actual SCCs lie
in every completion. In an actual tournament the triangle
components equal the SCCs: an external vertex has uniform orientation
to a directed triangle in a different component (otherwise it
makes a cross triangle), and this propagates along overlapping
triangles. Uniform orientation between all two components follows;
the quotient is triangle-free, hence transitive. Thus arcs between
possible-triangle components are irrelevant to H, not just
uninfluential in one candidate, and one can factor an exact target
across independent component completion problems.

A one-arc hereditary blocker is an especially small cube.
If completing its missing arc produces counts N−2 and N+2, with
N>5, then **no** supertournament containing that partial template
can have count N: its low side cannot increase by two, and its
high side cannot decrease. Here rows are outgoing adjacency
bitmasks in zero-based order. Unassign the indicated arc:

| N | Rows | Arc | Completed counts |
|--:|---|---|---|
| 7 | (22,28,24,17,0) | 2–4 | 5,9 |
| 21 | (22,60,56,49,32,1) | 0–4 | 19,23 |
| 35 | (14,52,56,50,33,1) | 1–4 | 33,37 |

The last row demonstrates why this is a hereditary *template*
obstruction and not a global new forbidden value. To prove coverage
using such a search, one must still prove a target-dependent structural
statement: for example, complete the cubes of a new replaceable-reader
grammar, not merely rule out cubes for the fixed induced-extension
grammar.

Sources: [Chan and Pak, §6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf);
[Moon, vertex-pancyclicity](https://doi.org/10.4153/CMB-1966-038-7).

The optional [weighted cycle-pack splitter](packing_certificate.py)
reproduces both known exclusions without an isomorphism atlas.
Fix a spanning directed cycle in the strong factor (Camion);
unassigned arcs remain free. Sum weights of forced packs for a
lower bound and all still-compatible packs for an upper bound.
Split on one free arc until each branch excludes the target.
For H7 one need only check strong orders 3–4 by Busch, as above.
For H21 excluding factor 7 leaves only one nontrivial strong
factor of count 21, of order at most 6 by the same bound.
The script recalculates every interval and verifies the branch
volumes add to the whole directed-cycle cube. At order 6 it
uses 216 leaves instead of 512 full orientations in that cube.
These counts certify finite kernels; unlike the vanishing
reflection cross-term, the splitter is not a universal inverse.
