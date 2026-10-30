# Two cyclic parents conceal a collective phase of unbounded arity

**PROVED:** the group, its exact two-closure, common tournament banks, and the first arity detecting its missing phase. **FINITE-EXACT:** the accompanying checker. This extends the composite-chart hostile in [the preceding intersection packet](../../layer_lifts_cyclic_intersections_20261030/INTERSECTIONS.md), rather than extending its prime-degree theorem falsely. No unrestricted tournament enumeration or major-conjecture closure follows.

## C1. An explicit group between two cyclic parents and their pairwise closure

Let r,m≥3 be odd, n=rm, and label vertices by x=i+rj with i∈Z/r and j∈Z/m. Put

    a(x)=x+1,
    v(x)=x+r if x=0 mod r, and v(x)=x otherwise,
    b=v a v^(-1),   G=<a,b>,   g=gcd(r,m).

The two cyclic parents <a>,<b> are distinct. Let W=C_m wr C_r act by arbitrary independent translations inside the r fibres, followed by a cyclic permutation of the fibres. Write its elements as

    (i,j) -> (i+t, j+z_i),    t∈Z/r, z∈(Z/m)^r.

Use representatives 0≤t<r. The map

    chi(t,z)=sum_i z_i−t mod g

is a surjective homomorphism W→Z/g. Wrapping t modulo r changes this expression by a multiple of r, hence by zero modulo g. Then

    G=ker chi,       |G|=r m^r/g,       |W|=r m^r.        (C1)

**Proof.** The element a has outer shift1 and exactly one internal carry1, so chi(a)=0. Conjugation preserves chi, giving chi(b)=0. The commutator b a^(-1)=v a v^(-1)a^(-1) translates two adjacent fibres by opposite amounts. Its a-conjugates generate the entire augmentation subgroup of internal translations with sum zero in Z/m. Moreover a^r translates every fibre by1. These generate precisely the internal vectors whose sum belongs to r(Z/m), equivalently whose sum is zero modulo g. Every possible outer shift is supplied by a. This proves equality with the kernel, including the order. The two parents are distinct because b's internal carry pattern differs from that of any power of a with outer shift1; alternatively a power of a giving b would have to have the same increment at every named vertex, whereas b has unequal increments.

## C2. Pairwise data forget exactly the collective congruence

For any ordered tuple of vertices that omits at least one fibre, G and W have the same orbit. Indeed an arbitrary W transporter can have its chi value corrected by a translation in an unused target fibre, without changing any tuple entry. If a tuple meets every fibre, its W pointwise stabilizer is trivial: one point fixes the outer shift, and a point in each fibre fixes every internal translation. Its W-orbit consequently splits into exactly g G-orbits.

Therefore, if g>1:

- G and W have identical orbits on all ordered t-tuples for every t<r;
- they have different orbits already at t=r;
- on a named transversal with one coordinate in each fibre, the extra invariant is the sum of internal coordinates modulo g after aligning the ordered fibre labels.

When g=1 there is no missing phase: G=W.

In particular the **two-closure** of G, the full permutation group preserving each directed pair orbit, is exactly W. To see the upper bound, one within-fibre directed difference-one relation consists of r disjoint directed m-cycles, so its connected components identify the fibres. The cross-fibre difference-one relation forces a cyclic translation of those fibres. Within each fibre the directed m-cycle relation forces an internal translation. These are exactly W. The lower bound follows from equality of pair orbits, since r≥3.

Taking r=m through unbounded odd values shows that no fixed tuple arity detects every collective phase in this family. This is a statement about these explicit orbit consumers, not about all possible succinct descriptions: chi itself needs only a short congruence formula. It is the same distinction between small symbolic state and a faithful bounded-arity relational observer that matters in the recent return-observer bounds.

## C3. Exact common banks and native-chart intersections

The common G-invariant tournaments are exactly lexicographic products

    S[R],

where S is a circulant tournament on r named fibre labels and R is a circulant tournament on m named internal labels, with the same R in every fibre.

**Proof.** The augmentation translations act transitively on ordered cross-fibre pairs with a prescribed ordered fibre pair: two internal displacements can be chosen independently and a third fibre absorbs their sum. Within a fibre they preserve and realize every internal translation. The outer a action identifies equal internal differences in all fibres and equal cyclic differences of fibres. Thus these are precisely the directed pair orbits. Assigning orientations independently to each antipodal pair gives S and R. Conversely every such product is fixed by W and hence G.

There are `(r−1)/2` outer and `(m−1)/2` inner orientation coordinates. Consequently

    |F_<a> ∩ F_<b>| = 2^((r+m−2)/2).                   (C2)

The natural Hamilton path P=(0,1,...,n−1) is legal exactly when the outer difference-one orientation points forward. No inner orientation is pinned. Applying the SAME legal path complement to both parent banks translates their common restricted affine flat injectively. Their resulting native charts therefore intersect in exactly

    2^((r+m−4)/2)                                      (C3)

tournaments, not merely at least that many. This recovers the earlier r=3 lower-bound family and supplies its exact intersection plus a two-parameter extension. At r=m=3 the chart intersection has2 points; at r=3,m=9 it has16. These are named charts, not graph-isomorphism counts.

Their common homogeneous direction has dimension `(r+m−2)/2`; the shared path removes exactly one outer coordinate. The large overlap is not repaired by carrying that same path's endpoints or ports. The unobserved inner directions are actual independent tournament coordinates.

## C4. Exact marking counts reveal the arity threshold

Let A_W(t) and A_G(t) count orbits of ordered t-tuples, with repetitions. For t≥1,

    A_W(t) = (1/r) sum_(j=1..min(r,t)) (r)_j S(t,j) m^(t−j),

where S(t,j) is a Stirling number and (r)_j is falling factorial. The block pattern has `(r)_j S(t,j)` choices. Each occupied fibre leaves `m^(number_of_entries−1)` relative internal choices, and the free outer cyclic translation divides by r.

For t≥r,

    A_G(t)−A_W(t) = (g−1)(r−1)! S(t,r) m^(t−r).         (C4)

For t<r the difference is zero; A_G(0)=A_W(0)=1. Formula(C4) follows because tuples meeting every fibre are a free W-set and every one of their W-orbits splits into g G-orbits. For r=m=3 the first positive-t counts are W:1,5,29 and G:1,5,33. The pairwise invariant is identical, while the three-point consumer detects the missing global phase.

## C5. Why this matters, and what it cannot do

A tournament automorphism group is two-closed: preserving its directed edges is a binary-relation condition. Thus, when g>1, **no tournament on this named vertex set has full automorphism group G**. Every G-invariant tournament has the larger W symmetry. A specification by two cyclic generators therefore does not by itself identify the full graph symmetry. Conversely an r-point relation recording chi can distinguish G from W; changing the consumer changes the faithful symmetry group.

This sharpens two existing boundaries for enumeration: common-parent charts at composite order can overlap in large affine flats, and pairwise symmetry can miss a global phase even when every lower-arity marking is retained. For an open-problem reduction, one must decide whether the target is a pair relation, a named r-port constraint, a positive existence predicate or an exact orbit count before quotienting. These results supply an explicit test family, not an inference about LRC, Jacobian or Collatz truth.
