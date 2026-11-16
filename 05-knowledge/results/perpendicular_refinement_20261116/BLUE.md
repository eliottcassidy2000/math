# A blue hub with a converse-symmetric five-line refinement

## Two different three-point quotients

In the [archived tiling gauge](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-643-blue-black-line-parity-structure.md),
declare 0→1→⋯→n−1 the distinguished Hamilton path. Nonconsecutive pairs
(i,j) are the frame coordinates; their bits record arcs backward from that
order. Reflection σ sends (i,j) to (n−1−j,n−1−i) and induces converse.
Global bit flip τ gives the physical line {t,τt}. These involutions commute.
A line is blue when both frame endpoints are individually fixed by σ.
In particular a blue line cannot end at a chiral tournament.

Take a blue frame u, a reflected coordinate pair e,e′, and a distinct
reflection-fixed coordinate k. In the three-cube select only
u, u+e, u+e′, u+e+k, u+e′+k. In the real cube, the converse direction
(e−e′)/√2 and the fixed split direction k are genuinely orthogonal.
Forgetting k collapses two parallel bars, leaves the center alone, and
gives a signed three-point quotient. Merging by σ gives a different
three-point quotient: it keeps the split level, forgets the sign, and has
orbit sizes 1,2,2. These operations are not interchangeable. In the full
cube, dropping k would also split the center into two preimages. The five
points are a selected subfamily, not the whole cube.

## An actual connected star at every order divisible by six

Let n=6m and M=2m. Define the hub frame C by partitioning the vertices
into residue classes A,B,D mod 3, orienting A→B→D→A, and using decreasing
order inside each class. This is C₃[TT_M,TT_M,TT_M]; the consecutive
Hamilton path uses cyclic neighbor classes. Set u=τC. The reflection is
an anti-automorphism of C and fixes its frame. Use coordinates
e=(0,3), e′=(n−4,n−1), and k=(3(m−1)+1,3m+1).

Reflection swaps e and e′ and fixes k. These coordinates are adjacent
internal arcs in distinct hub modules. Toggling any subset in C merely
changes the transitive orders by disjoint adjacent swaps. Thus the five
τ-lines share the same hub class: one blue line to u and four black
lines to the arm classes. These are actual line-multigraph edges, not
merely adjacencies in a cube diagram.

For distinctness, u has outdegree 2M−1−floor(i/3) at each nonendpoint i,
unique maximum 2M at the first vertex, and unique minimum M−1 at the last.
The τ-flipped cyclic directions supply M outneighbors outside the
module, with one loss at the predecessor and one gain at the successor
of the Hamilton path; exactly one correction is absent at each end.
Toggling e removes the unique max, while e′ removes the unique min; k
only swaps two consecutive degrees in the middle. The hub has neither
extreme, so the min/max tests separate it, the blue leaf, and both arms.

For m≥3 the two heights on the e side are distinguished by the
outdegree-(2M−m) induced subtournament. Before the k toggle it is row
(3(m−1),3(m−1)+1,3(m−1)+2), which is transitive. After toggling k its
middle vertex is replaced by 3m+1, making a directed triangle. The
first and last row vertices remain linked forward; the new vertex
beats the first and loses to the last since both pairs are nonconsecutive.
Converse handles the other arm. At m=1 the heights have Hamilton counts
37 versus 33. At m=2 the multiset of sorted outneighbor degrees at
outdegree-7 vertices distinguishes them: signature (3,4,5,5,6,6,7)
occurs twice before k and once afterward. The replay checks these low
cases directly and controls the row formula through m=20.

At order six order the ten free pairs lexicographically, lowest bit
first. The five leaf masks are 68,70,324,102,356. Their complements
955,953,699,921,667 are all isomorphic to the common hub. Its Hamilton
count is 45; the five leaf counts are 33,37,37,33,33.
In the unmerged class metagraph these five physical lines form a
selected star, not necessarily the entire set of incident lines.
In the converse-merged multigraph the pairs become parallel, so the
support star has three edges with multiplicities 1,2,2.

## Why five classes need not already be five incident lines

At order five, a minimal frame witness has u=0, e=(0,2), e′=(2,4),
k=(1,3), and the five masks 0,1,32,9,40. They have distinct sorted
outdegrees and Hamilton counts 1,3,3,5,5, where distinctness uses
degrees, not counts. Their complement ends do not all represent one
hub, so these lines do not form a star. There are only four classes at
order four: this is minimal for the five-class frame statement, not
a claim of global minimality for all ways to select incident lines.

Run `python3 verify_blue.py` for literal permutation checks at order
six and direct larger controls using the module and degree proofs.
