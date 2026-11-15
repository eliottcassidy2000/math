# Move the contact, then test the reset

The transfer below distinguishes a genuinely invertible change of reader
from a rank-one reset. It rules out a finite alphabet on both hard
factorization lanes by a weighted, all-order argument. An explicit
triangle block breaks the algebra behind the no-go; it is a working
way to enrich the alphabet, not a coverage claim.

## State and construction

For a tournament T with exposed vertex u let g=H(T). Let s,d be,
respectively, the total and odd-minus-even sum of nonempty paths
starting at u, each weighted by the Hamilton count of the complement.
Build T→R as an ordinal sum, then add a new exposed vertex v. Of T it
beats just u; of R it beats a prescribed marked subset. Write b=H(R),
B=H(R+v), Q for the ordered two-cover count (empty parts allowed), s_R,e_R
for total and signed starts at the marked vertices,
and t,f for total and signed ends at the unmarked vertices. Then

    G = B g + (b+t)s + (f−b)d,
    s' = Qs + (b+s_R)g,
    d' = (b−e_R)g.                                  (1)

Indeed on T alone the new vertex contributes s−d to g: an odd cycle
through v uses exactly a nonempty even path that starts at u, with
weight two. The odd-cycle packing expansion separates the crossing
cycles into a product of total and signed path counts, giving (1)'s
first line. For paths starting at v, the singleton contributes bg.
A longer path either first visits u, followed optionally by a path
in R, or first visits a mark of R. The former contributes Qs, and its
signed sum cancels because the unrestricted optional path in nonempty
R has signed weight zero. For this cancellation use the two-color
version of the odd-cycle packing expansion: each block in a packing
is odd, and assigning it to either path with a sign on the first
color yields opposite weights. The latter contributes gs_R and
reverses the sign.

With an empty R the exact transformation is E, and singleton R marked
0 or 1 gives A or B respectively (the letter B now denotes an operation):

    E(g,s,d) = (g+s−d, g+s, g−d)
    A(g,s,d) = (g+2s, g+2s, g)
    B(g,s,d) = (g+s−d, 2g+2s, 0).

Each E adds one vertex; each A/B adds two, with literal arcs as above.
E is invertible over the integers, with inverse
(s+d−g,g−d,s−g), determinant one, and characteristic polynomial
x³−x²−x−1. Thus the alphabet really adds a full-rank operation to
the two singleton choices, not another nominal label on a reset.

Every physical exposed state satisfies s≥g≥d≥0. To see this, for
any nonempty spine the sum of all signed weighted starts is its
Hamilton count, by unrestricted optional-path cancellation and the
empty-path term. Each individual signed weighted start is nonnegative
by induction: at u it equals H(T−u) minus the sum of signed starts at
outneighbors in T−u, a sub-sum of nonnegative values whose total is
H(T−u). Use signed start 1 for the singleton. Removing u and
re-exposing it gives s−g equal to exactly that nonnegative sub-sum.
The bound d≤g follows as well (or from the total of the nonnegative
per-vertex signed starts).

After one E/A/B move the closed cone 1≤s/g≤4, 0≤d/g≤1 holds and
remains invariant. For E use (g+s)/(g+s−d)≤2; A has new s=G, and B
has new s at most 4G. The first coordinate never decreases;
if s>d it strictly increases under every move. Once s>d holds it
persists. Nothing in this reduction assumes the original seed strong.

## Exact resets, a rewrite, and prime/7p word forms

Words below are chronological. Matrix multiplication yields

    BE = EA,
    AA(g,s,d) = (g+2s)(3,3,1),
    BB(g,s,d) = (3g+3s−d)(1,2,0).

The rewrite BE→EA removes one B, so normalization terminates, and
preserves the whole numerical state, not just the scalar count. It
need not preserve the tournament. The normalized word still produces
a physical tournament by direct assembly.

AA makes the count divisible by 3 and gives a common factor, hence
excludes prime or 7p outputs except possibly p=3. A BB reset makes a
common integer factor C=3g+3s−d≥5g. It is itself attained by the physical
prefix ending BB. A nonempty suffix acts homogeneously and strictly
increases the first coordinate from the representative (1,2,0). For a
prime, BB must therefore be terminal if present. For 7p and p≠3, if
an internal BB appears then C=p (7 is impossible and C=7p leaves no
room for a strict increase), and the remaining suffix realizes a
dimensionless multiplier 7. The suffix has bounded length: writing
W=g+s, on the cone each move increases W by at least 3/2, and if the
final multiplier is 7 then W≤35, starting at W=3. Taking the first BB,
every prime or 7p output therefore comes from a word avoiding BE,AA,BB,
followed by one of finitely many bounded tails (including the empty
tail). The bounded-tail reduction does not assume smoothness,
equidistribution, or uniqueness of reachability.

## Four-letter contraction and prime-lane sparsity

For a word j of length four avoiding BE, AA and BB internally, let
r_j=W(jx)/W(x). There are 24 literal words. For fixed W=1 their
reciprocal gains are convex functions on the cone slice. Its four
vertices correspond to (g,s,d)=(1,a,b) for a∈{1,4}, b∈{0,1}, scaled
to W=1. The maximum of the sum is thus at a vertex. Exact arithmetic
in the replay verifies the four strict upper bounds, each less than
7/8 (approximate values 0.769862, 0.872636, 0.710647, 0.745025).

Each single step increases W by at most 4 on the cone, so r_j≤256.
For u=49/50, the total gain weight is bounded by
(7/8)256^(1/50)<1, equivalently 256·7^50<8^50.
Group a normalized word into four-letter blocks after its first move.
Summing W_end^(−u) over all whole blocks is a geometric series,
regardless of which projective state a previous block left. Ignore
cross-block exclusions for an upper bound. The first letter, leftover
0–3 letters, and bounded tails have only finitely many options and
cannot decrease g. Since W≤5g on the cone, the number of possible
prefixes ending at g≤X is O_seed(X^(49/50)); the same bounds the prime
and 7p outputs. This works for any fixed finite seed list. In
particular the grammar cannot realize all primes, by the ordinary
prime-counting asymptotic. The theorem does not apply to arbitrary
changing exposed vertices, larger blocks, or a seed bank growing
with the target.

## Bypass the reset with a cyclic block

In (1), choose R to be a directed triangle and mark one or two of
its vertices. The two choices, C and D, give

    C = [[5,9,−1],[6,12,0],[2,0,0]],
    D = [[5,6,−2],[9,12,0],[1,0,0]].

Here b=3,Q=12,B=5, each vertex contributes 3 to the unsigned starts or
ends and 1 to the signed starts or ends. Both matrices have determinant
24, so neither ever collapses projective states or uniformly forces a
common factor in the manner of AA or BB. Starting from the exposed singleton,
C gives the prime 13; EC gives 23; EAC gives 69, in chronological
order. Independent subset DP in the replay checks each assembled
graph (including strong examples) and all words of length at most
three. A monotone-state census finds 141 distinct counts up to 1000
from the singleton with E/A/B and 267 with C/D added. The improved
alphabet still misses 11,19,33 there; this is a scoped finite result
(the first coordinate is nondecreasing). Thus removing the reset
increases the supply and escapes specific holes, but it is not
a universal inverse compiler. One next test is to allow **two**
old contacts so normal forms no longer follow a single spine of
exposures; for the stronger enlarged one-contact alphabet one also
needs a new entropy/overlap argument rather than the E/A/B bound.

## An unbounded inverse grammar that keeps the missing coordinate

After an A move the state is (x,x,y), with x,y odd and
z=(x−y)/2 a positive integer. As a further redesign, return to this
same two-dimensional surface by applying one of EA, BA, DA.
In coordinates (x,z) the three physical returns are

    EA: (5x+2z, 2x)
    BA: (9x+2z, 4x)
    DA: (51x+4z, 21x).

This gives an exact arithmetic decoder of arbitrary depth on the
**two-coordinate** states. Since x is odd, v₂(z') is respectively
1, 2, 0, so terminal moves are disjoint. Divide z' by 2,4,21 as
appropriate to recover x and solve the first equation for z; reject
if integral positivity and the physical precondition fail.
Values of v₂(z')≥3 reject outright except at a given seed.
Every accepted inverse strictly decreases x, so from any finite
chosen seed bank the algorithm terminates and either gives a
literal construction or excludes that precise profile from the
return grammar. For a single seed its orbit is a free rooted
ternary tree with no joins: decoding separates letters, and the
root cannot recur by strict growth. The recurrence recovers an
absolute profile rather than a residue. It still needs a rule
choosing z from a prescribed scalar x, which this decoder does not
pretend to supply. There is an elementary global bottleneck even
before that selection: each digit multiplies x by at least 5, so
this three-return family alone reaches at most
O(X^(log₅3)) possible profiles through x≤X.

## The exact state stops being sufficient at two contacts

Writing rows as outgoing bit masks, compare the strong order-five
tournaments (6,12,24,17,3) and (6,12,16,5,11), each exposed at
vertex 4 (zero-based). Both have the same state (15,21,3) and agree
on every E/A/B/C/D continuation. In the first graph every vertex
has (s,d)=(21,3). In the second the other vertices have
(15,3) once and (19,3) three times. If the new apex beats 4 and
one distinct old vertex, and all other old vertices beat it, the
new count is H+s_4−d_4+s_j−d_j−K_{4j}−K_{j4}, where K_{ij}
is twice the weighted count of even paths starting at i and ending
at j. An apex cycle must return through an **unmarked** endpoint:
adding the single-root increments without the two cross terms
would incorrectly predict 51 versus {45,49}. The adjusted palettes,
checked by literal subset DP, are {41,43} and {33,43}, respectively.
The extra contact opens new branches and is not encoded in a
single-reader state. A valid universal construction pursuing that
route needs not just each individual root profile, but its joint
endpoint interaction with the other contact.
## A single cyclic move still cannot supply the prime lane

The two enrichments C and D individually fail at all depths, even
though they defeat the rank-one reset argument. For either alphabet
EABC or EABD the conclusion for prime and 7p outputs strengthens to
the scoped bound O(X^(999/1000)) for every finite seed menu.
The same cone is invariant: for C, s'/G=(6g+12s)/(5g+9s−d) is in
[1,4] and 0<=2g/G<=1, and for D use (9g+12s)/(5g+6s−2d)
and 0<=g/G<=1. Each step increases W by a factor in [3/2,21].
The reset normal forms still apply because C and D are homogeneous
and increase g strictly on the reset representative and the cone.

Use eight-letter blocks instead of four. In either alphabet there
are 16897 admissible words. At the four cone corners, each sum of
reciprocal W gains is strictly less than 19/20. The companion
`replay_cyclic_pressure.py` uses integral **upper** rounding of each
term to a multiple of 10^-9, so this is a rational inequality, not
a numerical estimate at the exponent 1. Convexity and
`21^8 * 19^1000 < 20^1000` then give strict contraction at exponent
999/1000. The proof says choosing just one of the new cyclic moves
is insufficient: any successful construction in this direction
must genuinely mix at least two different cyclic types or add
another mechanism. This does not exclude the joint EABCD alphabet.
