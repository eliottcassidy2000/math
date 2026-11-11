# Complete minimum-defect families for every indegree at least three

**PROVED IN NOTE from inherited positive exposure capacities; independently
reviewable elementary proof.** Literal finite controls are in `minimal_defect.json`.
This is an all-order classification of one extremal stratum, not a source-return
or Chan–Pak theorem. The path-count formula for a cyclic collar is a recovered
run decomposition; the new conclusion is that the minimum-defect stratum
exhausts exactly these collars.

## 1. Positive defect and multiplicative exposure

For a strong rooted tournament (T,r), write I=N^-(r), d=|I|>=2,
R=T-({r} union I), U=T-r, and b=F_r(2)-H(T). The inherited identity is

    b = sum_(a in I) s_U(a) + sum_{a,c subset I} c_ac(U),
    c_ac(U)=W_ac(U)+W_ca(U).                                (1)

Here W_ac counts ordered covers by two nonempty paths, the first ending at a,
the second starting at c. The first sum is q-s>=0. These definitions and the
rooted identity are reviewed in [the fixed-defect proof](../defect/PROOF.md).

[THM-4167, exposure-capacity deletion support](../../../../01-canon/theorems/THM-4167-tournament-exposure-capacity-deletion-support-moment-and-parity-holonomy.md)
expresses c_ac as a positive sum of 2^(number of cycles) over tagged odd-cycle
packings. The distinguished cycle uses a new ear with local segment
c->ear->a or a->ear->c; the tag retains which ear orientation is used.

A tagged packing on I+ear and an ordinary odd-cycle packing on the disjoint
vertex set R have a disjoint union in U+ear. The map is injective, keeps its
tag and multiplies the weights. Consequently

    c_ac(U) >= H(R)c_ac(I).

Summing over pairs and over the ordered two-path covers of I gives

    b >= H(R)(Q_I(2)-2H(I)) >= H(R)(2^d-2).                (2)

The last sum is over all nonempty proper subsets J of I, with term
H(I[J])H(I\J)>=1. The multiplicative bound is a direct corollary of the
inherited tagged positivity, not a claim that a single W_ac is monotone.

## 2. Exact equality classification

**Theorem.** For every d>=3, a strong rooted tournament of indegree d has
b>=2^d-2. Equality holds exactly for the cyclic ordered blocks

    r -> TT_k -> I -> r,       k>=1,                       (3)

where I is either TT_3 or C3 when d=3, and I=TT_d when d>=4.
All arrows between the three displayed blocks are uniform. This is a
classification up to rooted ordinary isomorphism, with every k retained.

**Proof.** Equality in (2) forces H(R)=1, so R is a transitive chain, and

    q-s=0,   Q_I(2)-2H(I)=2^d-2,
    c_ac(U)=c_ac(I) for every unordered pair {a,c} in I.    (4)

The positive tagged interpretation says in particular that there can be no
new tagged odd-cycle packing using a vertex of R.

Suppose a vertex x in R has mixed contacts to I. Since d>=3, choose three
inports with two on one side of x and one on the other. If b,c both beat x
and x beats a, order b,c so b->c. Then

    b -> c -> x -> a

is a directed path on four vertices with inport endpoints. If a beats x and
x beats b,c, order b,c so b->c, and use

    a -> x -> b -> c.

Close this path by the corresponding tagged two-edge ear. It becomes an odd
five-cycle using x, contributing positive mass to the endpoint capacity that
was absent on I. This contradicts (4). Thus every x in R either beats ALL I
or loses to ALL I.

If an earlier chain vertex x loses to all I and a later y beats all I, then
for any two distinct inports a,c,

    a -> x -> y -> c

again closes to a new tagged five-cycle. Hence the chain contacts have form
all-ones followed by all-zeros. The last vertex of R cannot have all-zero
contacts: it would be a sink of T. Also R is nonempty because otherwise r is
a sink. Strongness therefore forces every chain vertex to beat all I.
This proves the cyclic block form (3).

Equality in the proper-subset sum in (4) makes every proper induced subset
of I transitive. At d3 this imposes no restriction: all one/two-vertex subsets
are transitive, so I is TT3 or C3. At d>=4, every triple is proper, so I has no
directed triangle and is transitive.

Conversely every graph in (3) is strong. Every path of U=TT_k->I starts in
the chain, giving q=s=H(I). In a W_ac cover, the second path starts in I and
cannot return to the chain. Thus every chain vertex belongs, in its unique
order, to the first path before its I segment. Removing the chain is a
bijection with W_ac(I). Hence b=Q_I(2)-2H(I)=2^d-2. QED.

The boundary d2 is essential. There are only two inports, so the mixed-contact
five-cycle argument has no third inport. The inherited five ready-tail families
give additional minimum-defect sources at b2. The theorem does not erase that
exceptional pentagon; it explains why the corresponding freedom disappears
starting at three inports.

## 3. Counts, symmetry and composition

The root is intrinsic in the marked problem, I is its in-neighbourhood, and
the residual chain has a unique order. Consequently (3) has one rooted type
per k when I=TT_d. At d3 there are two, distinguished by the actual inport
tournament. Their rooted automorphism groups are respectively trivial and C3.
Therefore the minimum-defect rooted-class series are

    d3: 2 z^5/(1-z),
    d>=4: z^(d+2)/(1-z).                                 (5)

The broader cyclic collar with any fixed inport tournament I has

    H(T)=(Q_I(2)-H(I))2^k+3H(I)-Q_I(2),
    b=Q_I(2)-2H(I).                                       (6)

For the path count, a Hamilton word with the root internal is uniquely
R_1,I_1,r,R_2,I_2: R_2 and I_1 are nonempty; R_1,I_2 may be empty. There
are 2^k-1 choices of the nonempty second chain subsequence, and the weighted
I split contributes Q_I(2)-H(I). The root-first and root-last cases contribute
H(I) each. This proves (6) without fitting a recurrence.

For transitive I, (6) becomes the symmetric Mersenne-product expression

    H=(2^d-1)(2^k-1)+2.                                  (7)

Converse exchanges d and k and therefore exchanges the two endpoint defect
values 2^d-2 and 2^k-2 while preserving H. This is an actual operation on
the supplied block source. It does not allow arbitrary factor choices in (7).
For I=C3, H=9*2^k-3, with defect6 at the displayed root.

## 4. The H141 test and exact stopping reason

The C3 collar with k4 has order8 and H141. In the gauge r0, inports1,2,3,
chain4,5,6,7, its rows are

    (240,5,9,3,238,206,142,14).

Its automorphism group has order3, inherited from the inport triangle.
The complete256 one-vertex attachment scan has neither425 nor427. Thus
minimum root defect6 does not automatically supply either ternary branch.
This is a one-birth failure only; no all-order induced-extension obstruction
is asserted here. The source may be tested with the older sharp frontier
method or with genuine replacement operations.

The connection contract is explicit: rooted source -> its inport tournament
and residual chain preserves the full graph exactly inside the equality
stratum. Formula (7) forgets d,k and the root, so a scalar factorization
requires those blocks as a construction certificate. The different d2 and
d3 boundaries are the hostile control for any attempted fixed-defect shortcut.
