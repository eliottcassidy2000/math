# Triangle inequality, curvature, and the capacity of a branching graph

**Status:** PROVED elementary estimates and explicit construction; FINITE-EXACT
replays; CITED standard hyperbolic metric/area inputs. No priority claim.
This packet does not prove Collatz convergence or construct a hyperbolic
hexagonal tessellation. The new gap construction concerns an abstract binary
tree, not every native Collatz edge. Independent review is recorded separately.

## 1. What the triangle inequality does—and does not—say

Let G be a connected graph with intrinsic unweighted graph distance d_G, and
let f:V(G) -> X map its vertices into a metric space. Suppose each graph edge
uv satisfies d_X(f(u),f(v)) <= L. For any shortest graph path, the triangle
inequality gives

    d_X(f(u),f(v)) <= L d_G(u,v),
    f(B_G(o,r)) subset B_X(f(o),Lr).                       (1)

For a directed predecessor set, the same containment follows from the
available directed paths to o, without asserting it is a full undirected ball.
One may instead bound drawn route lengths: endpoint distance is no greater
than route length, so (1) still follows. Crossing or noncrossing is irrelevant
to this necessary bound.

There is **no reverse inequality** in (1). It is an upper bound on ambient
radius, not a lower bound on it. An ambient geodesic may leave the drawn graph
and make a very large shortcut. Even an injective, properly embedded tree
whose edges are ambient geodesics of the same length need not satisfy a
coarse inverse distance bound. Section 4 gives an explicit hostile example.

Capacity needs one further input: distinct source vertices must occupy enough
separate target positions. For a lattice graph this means injectivity or a
uniform fibre bound; for a continuum it requires uniform separation of
occupied positions. Mere injectivity into the Euclidean plane gives no such
bound. Negative curvature does not relax the triangle inequality. It changes
how many separated positions can lie inside the ball forced by (1).

## 2. A quantitative curvature–branching bound

Write H^2_kappa for the complete simply connected plane of **constant**
curvature -kappa^2, with kappa>0. In curvature -1 its upper-half-plane metric
and distance are

    ds^2=(dx^2+dy^2)/y^2,
    cosh d(z,w)=1+|z-w|^2/(2 Im(z) Im(w)).                (2)

The area of a radius-R disk is 2 pi (cosh R-1). These are standard inputs:
[Anna Felikson, Geometry III/V, §6.2, Theorem 6.22](https://www.maths.dur.ac.uk/users/anna.felikson/Geometry/Geometry23/notes.pdf)
and [Stephen M. Buckley, Nonpositive Curvature and Complex Analysis, disk-area discussion before §2.9](https://archive.maths.nuim.ie/staff/sbuckley/Papers/nonpos.pdf).
Rescaling lengths by 1/kappa gives

    A_kappa(R)=2 pi/kappa^2 [cosh(kappa R)-1].           (3)

**Packing theorem.** Suppose an N-vertex source set has paths of at most r
edges to o, every image edge has length <=L, each occupied image position
has at most K source vertices, and different occupied positions are at least
delta>0 apart. Then

    (N/K)[cosh(kappa delta/2)-1]
        <= cosh(kappa Lr+kappa delta/2)-1.              (4)

Indeed, at least N/K different positions lie in radius Lr by (1). Their open
radius-delta/2 disks are pairwise disjoint: an intersection would contradict
the triangle inequality and separation. Another application of the same
inequality puts these disks in radius Lr+delta/2. Adding their areas proves
(4). Disk boundaries have zero area, so equality at separation delta causes
no difficulty. The analogous Euclidean estimate is

    N <= K (1+2Lr/delta)^2.                            (5)

For r>0, (4) gives the explicit necessary edge-scale bound

    L >= [acosh(1+(N/K)(cosh(kappa delta/2)-1))
                    -kappa delta/2]/(kappa r).        (6)

This is the exact inversion of this packing bound, not an assertion that
optimal disk packings attain equality. With no edge assumption, replace Lr
by the actual outer radius R of the occupied positions.

In particular, a full rooted b-ary tree (b>=2) has
N_r=(b^(r+1)-1)/(b-1) vertices through depth r. Taking logarithmic growth in
(4), with fixed finite K,L and fixed delta>0, proves

    kappa L >= log b.                                 (7)

More generally, source growth exponent h=limsup_r log N_r/r requires
h<=kappa L. For any separated placement of the depth-r ball, its outer
radius necessarily has liminf R_r/r >= log b/kappa. Formula (5) rules out
bounded L in the Euclidean plane. The curvature input here is quantitative
constant-curvature area growth; the phrase “negative curvature” alone is
not a substitute for (3). The threshold (7) is sharp as an infimum for b=2,
even with equal-length geodesic edges, properness, noncrossing, and a
quasi-isometric distance bound: see Section 5. The critical equality case
L=log 2/kappa is **not** settled here.

## 3. Native Collatz witnesses and the stopped-tree condition

Use the ordinary map C(n)=n/2 for even n and 3n+1 for odd n. Its predecessors
of m are 2m and, exactly when m=4 mod 6, (m-1)/3. Define

    I_r={n>=1: the first forward C-hit of 5 is at time <=r}.

The first-hit time, not a cyclicly reused time or an odd-only clock, is the
relevant distance upper bound. Every selected path is stopped at 5. The
forward orbit 5,16,8,4,2,1,4,... never returns to 5, so inverse generation
has no repeated depth. Exact predecessor enumeration and forward replay give

    |I_21|=704;  |B_honeycomb(21)|=694;  excess=10,
    |I_26|=2263; |B_triangular(26)|=2107; excess=156.    (8)

Thus the numerical difference 704-694 is **10, not 20**. These are the first
radius-count violations for this root within the checked prefixes, not a
claim of a globally smallest obstruction. The complete 2263 labelled states
and first-hit clocks are replayed against the inherited pinned TSV.

The honeycomb vertices are (a,b,s), s in {0,1}. A black (a,b,0) is adjacent
to white (a,b,1),(a-1,b,1),(a,b-1,1). The triangular lattice on Z^2 has steps
±(1,0), ±(0,1), ±(1,-1). Their exact ball sizes are

    H(r)=1+3r(r+1)/2;  T(r)=1+3r(r+1).                 (9)

For the triangular graph, a radius-r ball is |a|,|b|,|a+b|<=r and shell r
has 6r points. In the honeycomb, black vertices within 2k steps form that
triangular radius-k ball. White vertices within 2k+1 steps are its union
with shifts (-1,0),(0,-1); equivalently
-k-1<=a,b<=k and -k-1<=a+b<=k. This region has 3(k+1)^2 points, by removing
the two triangular corners of a 2(k+1)-by-2(k+1) integer box. At radius 2k
there are 3k^2 white vertices. These counts give (9) at both parities.
Independent BFS is an additional finite check, not the proof of (9).

An injective edge-preserving map would send I_r into the corresponding
radius-r lattice ball by (1), contradicting (8). This obstruction occurs
inside a certified convergent, cycle-free stopped subtree; it does not
rely on the 1->4->2->1 parity obstruction or on unknown Collatz behavior.

**General stopped-tree lemma.** For any function C on any set and root o,
the vertices whose first hit of o is finite, with the outgoing edge of o
removed, form a rooted tree. Every other vertex has one parent and its
nonnegative first-hit time decreases by one. Every vertex reaches o, so the
graph is connected. An undirected finite cycle would force each of its
vertices to use its unique outgoing edge on that cycle, hence would give
a directed cycle, contradicting decreasing first-hit time. The lemma
remains true when o is periodic in the un-stopped graph. Merely writing
“inverse tree” while retaining the root's outgoing cycle does not suffice.

For a finite rooted tree of distinct odd states with actual odd-map edges
u->q and a(u)=v_2(3u+1), expand the edge into

    u -> 2^a q -> ... -> 2q -> q.

If A(q) is the set of incoming exponents at q, distinct children have
distinct exponents. The exact union of native states has

    |V_odd| + sum_(q with children) max A(q)             (10)

vertices, not the count from independent subdivisions. The even chains for
different q are disjoint by unique odd part; chains entering one q are
nested. For a binary odd tree the independent-subdivision excess is
sum_q min A(q). The SF5 root-5 tree through odd depth 8 has 511 odd nodes,
2001 native union states, and 2811 independently subdivided states:
810 duplicate even states must be identified. These numbers are rechecked
here, together with all union clocks and the tree edge count N-1.

For a uniform asymptotic lower bound, recover the stronger prior five-step
supply rather than replace it by the older eight-step SF5 bound. For n not
divisible by 3 choose k according to

    n mod9: 1 2 4 5 7 8
    k:      2 1 2 3 4 1.

The endpoints F0(n)=32n and F1(n)=(16n-2^(4-k))/3 are distinct positive
integers prime to 3, both exceed n, and C^5(Fi(n))=n. The second reverse
route has k doublings, one legal odd predecessor, then 4-k doublings;
2^k n is 4 or 7 mod9. Distinct parents have disjoint five-step fibres.
Iterating gives 2^h distinct endpoints at time 5h, and the strict increase
also prevents cross-level collisions. For root 5, the 2^(h+1)-1 selected
endpoints through macro-depth h lie in I_(5h). Native intermediate tails
may merge; no independent-subdivision assertion is being used.

Consequently, if native images are delta-separated and every native edge
has length <=L in H^2_kappa, then

    kappa L >= (log 2)/5.                             (11)

A fixed fibre bound K does not change this exponent. For a honeycomb target
with integer edge bound L and fibre bound K, the finite inequality is

    2^(h+1)-1 <= K[1+3(5hL)(5hL+1)/2].                (12)

For h>0 it implies

    L >= {sqrt[1+(8/3)((2^(h+1)-1)/K-1)]-1}/(10h).

The analogous Euclidean lower bound is
L>=delta(sqrt[(2^(h+1)-1)/K]-1)/(10h). Thus unbounded route lengths or
unbounded fibres are forced in a fixed Euclidean lattice. In (11), meeting
the necessary curvature threshold is not a sufficient construction of the
native Collatz graph, and says nothing about universal forward convergence.

## 4. The inherited dyadic example: exact radial efficiency, no reverse bound

The prior construction places (k,j) at

    z_(k,j)=((j+1/2)2^(-k),2^(-k)), 0<=j<2^k.

Every edge has length ell=acosh(21/16); distinct vertices are at least ell
apart, with equality exactly for parent-child pairs. Geodesic edges give a
proper noncrossing embedding. The pinned prior proof supplies these all-level
claims; the checker here independently replays their rational certificates.

For the root z_(0,0), (2) gives

    cosh d(root,z_(k,j))
      =[( (j+1/2)2^(-k)-1/2 )^2+1+2^(-2k)]/(2*2^(-k)).

The farthest level-k points are the two outermost centers, so the exact
outer radius of all levels <=k is

    R_k=acosh[(5*2^k-2+5*2^(-k))/8].                  (13)

The expression increases with k>=0. Moreover

    k log2 <= R_k < k log2+log(5/4),
    R_k-k log2 -> log(5/4).                          (14)

The lower bound drops the horizontal term in (2). For the upper bound,
cosh(k log2+log(5/4)) minus the bracket in (13) is
1/4-(9/40)2^(-k)>0. Finally acosh t=log(2t)+o(1) as t tends to infinity,
which proves the limit. Thus the radial rate log2 is the smallest possible
for a uniformly separated binary-tree placement in curvature -1, although
the edge-path upper bound ell*k is larger.

Nevertheless, for k>=1 take j=2^(k-1)-1 and j+1=2^(k-1). Their tree LCA is
the root, giving d_T=2k, but their ambient distance is always acosh(3/2).
Thus no constants a>0,B<infinity satisfy d_H>=a d_T-B for all pairs.
Their intrinsic distance **along the embedded tree** is 2k*ell. This does
not conflict with the ambient triangle inequality: that inequality says
acosh(3/2)<=2k*ell, not the reverse. The embedding is separated and radially
efficient but not quasi-isometric. “Negative curvature fixes capacity”
therefore must not be inflated into “every natural drawing preserves tree
metrics.”

## 5. A strict sibling gap supplies the missing reverse estimate

Here is a bounded constructive repair, still for the abstract binary tree.
Choose 0<q<1/2 and 0<W<=1. Start with an interval I_empty=[0,W]. Replace
each interval [a,a+Wq^k] by its leftmost and rightmost subintervals of
length Wq^(k+1). For an address w of length k let x_w be its interval center
and place

    z_w=(x_w,q^k).

Let s=log(1/q). A child of (x,y) has coordinates
(x ± W(1-q)y/2, qy). Equation (2) proves that every edge has length

    ell(q,W)=acosh[1+(1+W^2/4)(1-q)^2/(2q)].          (15)

**Explicit quasi-isometric embedding.** For all vertices u,v,

    s d_T(u,v)-C <= d_H(z_u,z_v)
                   <= s d_T(u,v)+log(W^2+2),          (16)
    C=-2 log[W(1-2q)] > 0.

To prove this, let their LCA have depth h and their depths be h+a,h+b.
Then d_T=a+b. Both horizontal centers lie in the LCA interval, so
|dx|<=Wq^h, while y_u,y_v<=q^h. Equation (2), in the equivalent form
cosh d=(dx^2+y_u^2+y_v^2)/(2y_uy_v), yields

    cosh d <= (W^2+2)q^(-(a+b))/2.

Since d<=log(2 cosh d), this proves the upper estimate. If one vertex is
the LCA, vertical separation gives d>=s(a+b). Otherwise they lie in
different child intervals, separated by the **positive** gap
W(1-2q)q^h. The identity

    d=2 asinh(sqrt(dx^2+dy^2)/(2 sqrt(y_u y_v)))

and asinh t>=log(2t), t>0, give the lower estimate in (16). This is a
quasi-isometric **embedding** of the tree; no coarse density in all H^2 is
claimed. Its extra ingredient is relative subtree separation, not a
reverse triangle inequality.

For completeness, different levels are separated by at least s. At the
same level k>=1, neighboring interval centers differ by at least
W(1-q)q^(k-1), attained by siblings. To see the minimum, if their LCA has
depth h<k, the two closest level-k centers in its different children differ
by W[(1-2q)q^h+q^k], minimized at h=k-1. Thus every pair of distinct
vertices is separated by the positive constant

    delta(q,W)=min{s, 2 asinh[W(1-q)/(2q)]}.           (17)

Geodesic edges are noncrossing and proper. Put beta=W(1-q)/2. For an edge
from (x,y) to (x+beta*y,qy), its Euclidean circle center on the real axis is

    c=x+[(beta^2+q^2-1)/(2 beta)]y < x.

The strict sign holds since W<=1 and q<1/2. Both endpoints are to the right
of c, so the geodesic segment has monotonically decreasing height, staying
inside its endpoint horizontal interval and band [qy,y]. The left-child
case is its reflection. Different parents at depth k have center separation
at least W(1-q)q^(k-1), exceeding the sum W(1-q)q^k of their outgoing half
widths. Their horizontal edge intervals are disjoint. Different levels have
disjoint open height bands; band boundaries are attained only at vertices.
This proves noncrossing. A compact set has positive minimum height, so it
meets only finitely many levels and edges, proving properness.

At q=1/2 the sibling gap vanishes, C in (16) diverges, and Section 4 supplies
actual bounded-distance cross-root pairs. This is a genuine failure boundary,
not merely a bad estimate. If q tends to 1/2 from below and W tends to zero,
(15) tends to acosh(5/4)=log2. Every fixed parameter choice still has
positive separation and all the properties above. Therefore the necessary
binary threshold ell>=log2 is sharp as an **infimum**. Separation is not
uniform across this limiting family; (17) tends to zero as W tends to zero.
For small W, same-level siblings are closer than the equal edge length ell.
Thus this repaired family need not retain the inherited dyadic property
“minimum-distance pairs are exactly tree edges.” It trades that property
for a coarse lower-distance bound and an arbitrarily near-optimal edge scale.

More explicitly, for any desired lambda>log2, choose s<lambda sufficiently
close to lambda with s>log2, let q=e^(-s), and set

    W=2 sqrt[(cosh lambda-1)/(cosh s-1)-1].

Then 0<W<=1 for s close enough to lambda and (15) equals lambda exactly.
Rescaling to curvature -kappa^2 proves existence for every prescribed equal
edge length L with kappa L>log2. No attainment claim is made at equality.

## 6. A separate, genuine role for 20 in curvature

This is not the numerical gap in (8). In a finite triangulation of a closed
surface, counting incidences gives 3F=2E and sum_v deg(v)=2E. Euler's identity
V-E+F=chi therefore gives the combinatorial curvature sum

    sum_v (6-deg(v))=6 chi.                            (18)

If the surface is a sphere and every vertex has degree 5, then V=12, E=30,
F=20. Here 20 is forced by an actual curvature budget. The triangular
lattice has degree 6 and zero local defect; a closed regular triangulation
of degree >=7 has chi<0. A honeycomb is a hexagonal-face, degree-3 graph,
not a triangulation; its dual is the triangular lattice. These facts
explain a structural occurrence of 20 but give no new Collatz consequence.

## 7. Scope, inheritance, and reproduction

Closest inherited mechanisms: the constant-edge dyadic hyperbolic tree;
the finite radius-21/26 Collatz obstruction; the five-native-step binary
supply; and exact shared-tail source unions. The hostile inherited example
is the bounded-distance pair across the root split. The repaired near miss
is “proper separated drawing implies quasi-isometry”: it does not; a
uniform relative sibling gap supplies the missing estimate. The least-used
sidecar here is the actual native first-hit clock and source union, not the
unweighted odd-tree depth.

The small concept board is: intrinsic graph metric / ambient metric /
separation-and-capacity / curvature scale / first-hit and shared-tail clocks.
The connection map is the vertex placement f. It preserves source identity
when injective and sends short source paths to short target routes; it need
not preserve ambient lower distances or arithmetic forward behavior. The
cheapest decisive controls are the cross-root pair, a vanished sibling gap,
and literal native unions rather than independently subdivided macro edges.

Run `python3 verify.py` and `python3 -O verify.py` from any working directory.
The checker is stdlib-only, resolves pinned inputs relative to its own file,
and uses explicit exceptions rather than assert. All finite algebraic
geometry checks use rational arithmetic; printed metric examples use floats
and are labelled illustrative. The all-level arguments are the proofs above,
not extrapolations from finite tests. `SOURCE_PINS.json` records exact Git
inputs and external metric sources. `PUBLIC_FILES.json` and `SHA256SUMS`
seal the portable payload and its file selector; the checksum file excludes
only itself.
