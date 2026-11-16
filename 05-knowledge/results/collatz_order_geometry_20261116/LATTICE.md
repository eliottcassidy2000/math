# From the mod-nine hexagon to a genuinely hyperbolic tree

## Inheritance and the two meanings of "hexagonal"

The archive [Excursion forests: flat stars, boundary closure and unbounded tiling addresses](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/forest_20260926_tilings.md) explains
the 3 versus 5 transition on regular maps. For p-gons and q at a vertex, the
flat boundary is (p-2)(q-2)=4, giving (3,6),(4,4),(6,3).
The integer points below it are (3,3),(3,4),(4,3),(3,5),(5,3), the five
spherical cases. In duality language the flat side is one fixed point and
one pair; the spherical side retains one fixed point and gains a second
pair. The fixed point changes from the square tiling to the tetrahedron.
This is not a direct splitting map of individual lattices to solids.
The same archive has an important counterexample to local-to-global
reasoning: a flat 3.8.24 vertex forces alternation of 8 and 24 around a
triangle, which cannot close.

The [Collatz coherence note](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_coherence_20260930_trichotomy_hexagon_motifs.md) has another exact hexagon, the six units mod nine in doubling
order (1,2,4,8,7,5), with negation acting by the half-turn. That is a
finite address space, not a geometric tiling of all the Collatz edges.

## A two-child rule visible in the hexagon

Use f(n)=n/2 on evens and f(n)=3n+1 on odds. Let m be odd and prime to 3.
The inverse accelerated odd children are u_a=(2^a m-1)/3 over a≥1 with
2^a m≡1 mod 3. Exactly one parity of a is admissible. In its first three
exponents, either 1,3,5 or 2,4,6, the values 2^a m modulo 9 run through
1,4,7 in a cyclic order. This is the two-step action on the mod-nine
hexagon. Exactly one choice yields u_a divisible by 3, leaving exactly
two odd unit children. This gives a canonical left/right rule by taking
them in increasing a. All three are positive, and rooting the construction at 5 keeps it
away from the familiar accelerated loop at 1.

The resulting addresses at depth k give precisely 2^k distinct odd
integers. Indeed the accelerated forward Collatz map recovers the
parent uniquely, and 5 is not periodic, as its ordinary trajectory is
5,16,8,4,2,1 followed by the usual cycle. Each child arc expands to an
ordinary path of length a+1≤7. Thus the unconditional basin of 1
already has at least 2^k distinct nodes at distance at most 5+7k from
1. This follows from a local hexagon without making any assumption
that this basin exhausts the positive integers.

## Geometric consequence: impossible in a flat lattice, easy in H²

No injection of the basin into any fixed d-dimensional uniformly
discrete Euclidean set can send all tree edges to pairs at bounded
distance. If points have separation δ and edges have length at most L,
then k binary generations fit in the Euclidean ball of radius 7kL about
the image of 5. Disjoint δ/2 balls give

    2^k ≤ (1 + 14 k L/δ)^d,

which contradicts exponential growth. Even for finite depth k the
largest edge needs L≥δ(2^(k/d)-1)/(14 k). Square, triangular and
honeycomb lattices all have polynomial volume and so all fail the
same bounded-stretch test. This does not forbid drawings that use
arbitrarily long edges, unlabeled quotients, or multiple nodes over
one lattice point.

In contrast the inverse tree of the basin has at most two children
per vertex, after removing the outgoing cycle edge of 1. Assign a bit
to each child, with a fixed convention if there is one. For a vertex
of depth d with address w, take x to be the midpoint of the dyadic
interval for w and y=2^(-d) in the upper-half-plane model. Different
depths have hyperbolic separation at least ln 2, because the
length of any curve is at least the absolute change in ln y.
At the same depth, distinct midpoints are at horizontal distance
at least y, hence have distance at least arcosh(3/2). Parent to child
has differences in x and y of magnitude y/4 and y/2, giving

    cosh(distance)=1+((y/4)^2+(y/2)^2)/(2y(y/2))=21/16.

So the whole inverse basin has a uniformly separated placement with
constant edge length in the hyperbolic plane. This is a placement,
not a claim that the tree is the vertex-edge graph of a regular tiling
or that every arithmetic label is a geometric coordinate. The
conceptual transfer works because the local mod-nine operation
certifies branching; branching, rather than the word hexagon, picks
out the global metric geometry.

## What 3-to-5 counts, and how to vary it

Put D=(p−2)(q−2). Its integer shell is exactly the ordered factor pairs
(a,b) with ab=D, by (p,q)=(a+2,b+2). Thus the number of marked regular map
types on the shell is the divisor count τ(D). The 3 flat types are τ(4),
while the 5 spherical types are τ(1)+τ(2)+τ(3). Swapping p,q exchanges
factor coordinates; on shell D it fixes precisely one type if D is a
square, none otherwise. All-order generating functions (absolute
convergence for Re(s)>1) for marked vs up-to-duality shell counts are

    Σ τ(D)D^(−s) = ζ(s)^2,
    Σ #dual-orbits(D) D^(−s) = (ζ(s)^2+ζ(2s))/2.

The tiny 3/5 comparison is therefore a boundary slice through a divisor
lattice, not a direct geometric map of three plane lattices to five solids.
Even the curvature magnitude retains more information than the shell:
the excess of the (2,p,q) triangle divided by π on S² (signed in other
curvatures) is (4−D)/(2pq), and pq=D+2(a+b)+4. For a fixed nonflat shell
it is largest in absolute value for the most balanced factor pair and
smallest at (1,D), since the denominator is increasing in a+b.

### Change the observation: pass from faces/vertices to flags

Start with the (2,p,q) reflection triangulation, using a vertex for
each chamber and an edge across its sides. This is the chamber dual
map, and its three incident face lengths at every vertex are (4,2p,2q),
because around a corner with angle π/m there are 2m chambers.
The reflection criterion over S², E² and H² is classical, e.g.
[Michael Davis, Example 2.59](https://people.math.osu.edu/davis.12/IGAP.pdf).
Forgetting the colors that distinguish the last two slots collapses the
five spherical types to three figures (4,6,6),(4,6,8),(4,6,10).
The three flat types collapse to two (4,6,12),(4,8,8). The regular
3-to-5 count therefore becomes 2-to-3 when observed in unoriented chamber
figures. This is just a different count of objects, not a new Platonic
classification.

Conversely, consider any cubic map with three globally fixed, distinct
face sizes at every vertex and whose face boundaries are ordinary cycles
(as in an edge-to-edge tiling). Walking around a given face, the sizes
of its neighbors along edges must alternate between the other two
types. Thus every face size is even. Write the three sizes (2a,2b,2c)
with 2≤a<b<c. Nonnegative angle charge gives 1/a+1/b+1/c≥1.
If a≥3, the sum is at most 1/3+1/4+1/5<1; if a=2,b≥4, it is at most
1/2+1/4+1/5<1. Hence a=2,b=3 and c=4,5,6. The chamber construction
realizes all three. In particular, the exact recurring band is
(4,6,10),(4,6,12),(4,6,14), positive, flat, negative. This separates
an arithmetic vertex-charge criterion from the combinatorial closure
of neighbor words; local flat stars with odd sizes may fail closure.

## A literal 3-to-5 split in the admissible Collatz growth

Keep all three least-parity children, but tag the one divisible by 3
as terminal in this odd compressed model. It really has no odd
Syracuse preimage: 2^a m−1 never vanishes mod 3 for m divisible by 3.
The two other children are live. Expanding a live leaf therefore
replaces one live leaf with two live leaves and one terminal leaf.
Starting from the live 5, after m such choices the frontier has
m+1 live and m terminal leaves. Its total goes 1,3,5,7,...,
irrespective of the order of expansion or batching; the first
3-to-5 step is literally replacing one of the live leaves by its
triple. There is no claim that these counters are the Schläfli
types: the Platonic 3-to-5 is a shell count, and this one is a
frontier count.

The number of legal chronological expansion histories of length m
is m!: immediately before the i-th expansion there are i live leaves.
But the count of different rooted, ordered partial shapes (equivalently
different terminal marked frontiers in this fixed integer construction)
is the Catalan number C_m. To see this, suppress the forced terminal
legs and keep the two live slots ordered by exponent. Then shaped
histories are full ordered binary trees with m internal vertices,
and the usual root-splitting recurrence is
C_0=1, C_m=Σ_{j=0}^{m−1} C_j C_{m−1−j}.
The actual inverse acceleration has unique parents, so these shapes
are distinguished by their labeled frontier, not just by the
unlabeled abstract shape. For a fixed shape T, its chronological
multiplicity is m!/∏_v size(T_v), where size counts internal vertices;
this follows by shuffling the two child histories after their parent.
Neither formula says the compressed tree covers every odd integer:
it deliberately discards larger exponents and uses 5 as its root.
