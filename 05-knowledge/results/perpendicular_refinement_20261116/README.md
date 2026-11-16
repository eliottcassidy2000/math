# One fixed point, two converse pairs

The common operation is an **equivariant refinement**, not an identity
between the underlying objects. A three-point observer with coordinates
x=−1,0,1 has involution x↦−x. Keep the fixed center alone but refine each
of the two exchanged sides by the same second coordinate t. The resulting
cross has points (0,0),(±1,±1), and the involution sends (x,t)↦(−x,t).
So the two perpendicular rails have three properties worth checking:
the center is fixed, the wings are paired across the same t, and forgetting
t recovers the original three points. The count of fixed points matters
as much as the total count. The aligned split adds one even and one odd
dimension on permutation representations; the dimensions change from
(2,1) to (3,2).

## p²qr: when the visible square has no vertical lift

The [archived balance](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/arithmetic_braids_20260917_divisors.md)
has F=10, U=3, S=7. The new 3,5,7 are counts of *states* of those
ten divisors, not different values of F. For d=p^a q^b r^c, take
x=b−c. The center x=0 consists of p,p²,qr,pqr and is pointwise
fixed by Q, the exchange of q and r. On the wings retain whether
a=1 or a∈{0,2}, i.e. t=2|a−1|−1. The aligned pairs are
{pq}↔{pr} and {q,p²q}↔{r,p²r}.
The complement C:d↦p²qr/d and Q both flip x and preserve t.

This five-state chart is forced if both maps must descend, the center
must stay one class, and the two wings must each have two classes:
Q identifies their partitions, and complement forces invariance under
a↦2−a. The single-element class must be {1}. A vertical swap on the
five observed states does **not** lift to a divisor permutation, because
the two classes on a wing contain one and two elements. It is also
not the same quotient as the five complement orbits, on which Q has
three fixed points and just one exchanged pair. See [DIVISOR.md](DIVISOR.md)
and the archived [Klein dictionary](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/black_detours_closed_blue_lines_klein_dictionary_20261020.md).

Keeping the full signed a−1 on each wing yields seven points. Then Q
flips the horizontal coordinate only, C flips both, and CQ flips
only the vertical. Seven is minimal for seeing the faithful
Klein action with the center frozen at p²qr. More generally for
p^e qr with e≥1, a five-state chart of this kind can show the faithful
Klein action precisely when e is odd: dividing each wing at e/2 makes
C exchange the two vertical cells while Q preserves them. If e is
even the wing's fixed midpoint a=e/2 must belong to a fixed cell
under a↦e−a, so with only two cells neither can be exchanged.
A three-cell wing (below, at, and above e/2) restores faithfulness.
This is a parity obstruction about the observable chart, not a
statement that every other divisor observable is nonfaithful.

## A literal blue-centered 1+2+2 star at order six

Use the repo's [THM-643 sigma–tau frame](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-643-blue-black-line-parity-structure.md).
Fix the Hamilton path 0→1→2→3→4→5. Put bit 1 on a backward
non-path arc; list the free pairs lexicographically. The involution
sigma sends (i,j) to (5−j,5−i), inducing converse on classes.
The line flip tau complements the ten bits. A line is blue exactly
when its literal mask is sigma-fixed.

The common hub is the cyclic composition C₃[TT₂,TT₂,TT₂] on
A={0,3}, B={1,4}, C={2,5}, with A→B→C→A and both vertices in
each module ordered decreasing. This is frame mask 955. Start with
its line partner 68 and use coordinates e=(0,3), e'=(2,5),
k=(1,4). Sigma exchanges e,e' and fixes k. The five masks are:

| Leaf mask | Canonical class code | H | Role |
|---:|---:|---:|---|
| 68 | 145 | 33 | fixed center, blue line |
| 70 | 152 | 37 | first converse wing |
| 324 | 281 | 37 | its converse |
| 102 | 156 | 33 | second converse wing |
| 356 | 282 | 33 | its converse |

Canonical code is the minimum backward-arc bitmask over permutations,
using every pair in lexicographic order. The complements of all five
masks are in class 408 (H=45), since e,e',k are the three internal
module pairs, and flipping an internal TT₂ arc changes only its
labeling. It follows that these are five genuine lines on a common
hub, with one blue leaf and four distinct black leaves in two
converse pairs. The line multigraph retains five edges; converse
merging leaves three support edges with multiplicities 1,2,2.
The little star is a selected part of the metagraph, not the whole
neighborhood of the hub or an induced five-vertex pattern.

This star persists at **every order n=6m**, not just six. Use the hub
C₃[TT_(2m),TT_(2m),TT_(2m)] on the three residue classes, its complement
frame, and the three internal-pair toggles (0,3), (n−4,n−1),
(3m−2,3m+1). The first two toggle the unique maximum or minimum degree
of the central leaf. For m≥3 the middle toggle changes the relevant
equal-degree row in each wing from a transitive triangle to a cycle.
This separates both heights; m=1,2 have literal small certificates.
See [the uniform construction and proof](BLUE.md).

There is an actual vertical operation here, unlike on the divisors:
toggle k iff the e and e' bits differ. The guard is unchanged by
the toggle; sigma swaps the two guard bits and tau complements both.
Thus the operation is an involution of the framed cube commuting
with both sigma and tau. It fixes the center and swaps the heights
on each wing. On the selected five *lines* it acts as the second
reflection of the cross. The full cube operation is not claimed
to descend to an operation on all unframed tournament classes.

## The Collatz fork clock has exactly the same visible cross

Recover the [weighted binary skeleton and live hexagon](https://github.com/eliottcassidy2000/math/blob/9dec4d1efcbb1389b41b1b4116154b39950ed050/05-knowledge/results/collatz_order_geometry_20261116/README.md).
Use the shortcut inverse branches D(x)=2x and O(x)=(2x−1)/3
when x≡2 mod3, away from the root cut. Fork sites mod 9 have
types 2,5,8. Type 5 has a sterile odd child. After discarding the
sterile branches, live forks (types 2 and 8) on the double branch
make the mod-27 clock 2→8→20→26→11→17→2.
The double-edge lengths are 2,4,2,4,2,4.

Now reflect residues by R(x)=19−x mod27. On this hexagon it sends
index i to 5−i, reverses the clock, and preserves its undirected
edge lengths. It exchanges the live type-2 and type-8 sheets
and preserves the sterile residues {5,14,23}, fixing 23.
Keep the sterile sheet as one point. Refine each live sheet by
whether the next live fork along the odd branch requires an extra
two doubles (because its first fork is sterile). The delayed
live residues are exactly 11 and 8, exchanged by R. The no-delay
pairs are 2↔17 and 20↔26. Thus the 3→5 cross occurs in the
clock itself. Its one-versus-two fibers again prevent the
vertical reflection from becoming a literal bijection of phases.

This is not a symmetry of the whole Collatz dynamics. For example
the full odd-branch destinations on the six phases are two sets
A={2,11,20}, B={8,17,26}. R exchanges A,B, but the reflected
sources 2 and 17 both have destination set A. The delay bit
is retained; the next digit is not. A geometric picture that
omitted that distinction would mistakenly turn a reversible clock
into a reversible tree rule.

## The fork grammar is an exact distance error code

[DISTANCE.md](DISTANCE.md) gives the sharp bounds and exact correction.
A weight-e edge takes the outer 3^(-e) fraction of a parent's interval,
and the hyperbolic upper-half-plane coordinates are midpoint and
interval length. Use the [standard upper-half-plane distance formula](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Analysis/Complex/UpperHalfPlane/Metric.html).
For two fork vertices put a,b for their weighted lengths after the
common ancestor, A=3^(-a), B=3^(-b), u for their midpoint difference
in that ancestor's unit interval, and D=a+b. If Q=u²+A²+B²,
```
d_H − D ln 3 = ln[(Q + sqrt(Q²−4 A² B²))/2].
```
The error is a logarithm of one local expression, not a sum of
roundoff along the path. For incomparable vertices it is between
2 ln(sibling gap) and zero. For ancestor pairs it lies in
[0,ln(5/4)]. For an incomparable pair whose fork lengths are
(e,f), the four immediate lower bounds come from gaps
```
(e,f)     (2,2)  (2,4)  (1,4)  (3,4)
gap        7/9   71/81  53/81  77/81
```
The local bound depends only on the fork type; the exact error
uses endpoint prefixes as well. Reading K weighted digits after
the common ancestor truncates the correction with error at most
(2/gap)·3^(−K) for an incomparable pair (provided the branches
have been identified; short suffixes can be read exactly).
Indeed each missing ternary tail changes its midpoint by at most
half of 3^(−K); the derivative in u is 2u/sqrt(Q²−4A²B²)≤2/u.
The log of the gap is the midpoint of the locally certified
correction interval.

The worst gap 53/81 is just enough. The entire signed error
interval has width
```
ln[(5/4)(81/53)²] = ln(32805/11236) < ln 3.
```
Therefore for all fork vertices the following is exact, even
without reading their ancestor or fork type:
```
D = floor[(d_H + 2 ln(81/53))/ln 3]
  = floor[(d_H + ln(12/5))/ln 3].
```
You can compare rational cosh(d_H) to rational bin boundaries and
avoid floating point entirely. The simpler shift ln(12/5) follows
from (53/81)²>5/12 and exp(d_H)/3^D<5/4 for every finite pair.
If comparability is already given, floor(d_H/ln 3) suffices for
ancestor pairs and ceil(d_H/ln 3) for incomparable pairs.
There is an even shorter rational
decoder for distinct vertices: compute V=2cosh(d_H). The formula
V=3^D Q and the elementary bounds gap²≤Q≤5/4<3(53/81)² imply
```
(53/81)²·3^D ≤ V < (53/81)²·3^(D+1).
```
Here Q≤5/4 in the ancestor case follows from
1+B²+(1−B)²/4≤5/4 for B≤1/3. Equality of vertices can be handled
separately as D=0. Only integer or rational comparisons occur. Or center
the signed metric interval before
rounding and tolerate a distance measurement error of up to
(but not including) half ln(33708/32805), about 0.0136.
Also the sign of the error now recovers a combinatorial fact.
For distinct comparable vertices the midpoint difference is strictly
positive and the vertical-distance lower bound is strict; for
incomparable vertices Q<1 and the upper bound is strict. So the same
single hyperbolic distance reports whether the vertices are comparable
once D is decoded. Two distances to the root then give the ancestor
direction when comparable and the weighted depth of their common
ancestor via (D(root,u)+D(root,v)−D(u,v))/2 in all cases. In
particular the induced hyperbolic metric on the marked skeleton
determines its weighted edges: a pair is adjacent exactly when no
third skeleton vertex splits the decoded distance additively.
The earlier lossless expansion then reconstructs the rooted basin.

This explains why the **pairs** of lengths matter, not just the
alphabet {1,2,3,4}. For the same outer-thirds embedding, a single
common shift-and-floor decoder works whenever the worst sibling gap
exceeds sqrt(5/12). All pairs with both lengths ≥2 clear that
threshold. If one length is 1, its partner must be ≥4. The unphysical
fork (1,3) would already break this common decoder if both branches
are full: the lower gap bound and the upper ancestor bound are
attained in limits along extreme endpoint rays. The actual Collatz
grammar excludes precisely these short partners. The full tree
before compression has a separate universal decoder using outer
fifths; it is not claimed to satisfy the thirds bound. For unweighted binary
outer-edge placement, five is the smallest integer contraction base
with this universal one-bin guarantee: base b has sibling gap
1−2/b, and for b=3,4 its sharp global error width
ln[(5/4)/(1−2/b)²] exceeds ln b. Base 2 has no sibling gap.

## Exact controls and scope

Run `python3 check.py`: the arithmetic equivariance and Collatz
hexagon, literal order-six tournament class/HP certificate, and
33,331 pairwise rational distance tests (weighted thirds and
unweighted fifths). Optional independent replays:
`python3 check_divisor.py`, `python3 check_distance.py`, and
`python3 verify_blue.py` (uniform star controls to order 120).

The shared theme is an auxiliary coordinate that is invariant under
the named reflection but still resolves an ambiguity. In the tournament
it gives an actual controlled move; in the divisor and Collatz
observers it is a quotient coordinate, not a new automorphism upstairs.
In the metric it determines the local error budget. None of the
constructions asserts convergence of all integers or full Chan–Pak
coverage.
