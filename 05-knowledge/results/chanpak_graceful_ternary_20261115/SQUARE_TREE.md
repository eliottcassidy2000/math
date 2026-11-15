# A square-sum absorber tree with an exact tournament count inverse

**Proved:** a branch-preserving embedding of a full rooted ternary tree
into the *absolute* three-count windows of the literal Chan–Pak bank,
with an exact inverse on its image. An optional fourth branch gives
a full quaternary tree. This is not a new cofinite family: its image
has explicitly sublinear counting growth.

## From the Berggren shell to a fixed aperture

The [October 19 norm-one tree and square-sum absorber construction](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/square_sum_berggren_absorber_bijection_20261019.md)
uses positive solutions a²+b²=c²+1. Given an edge of the square-sum graph
with endpoints x,y and S=x+y square, put

    u=a²S−x,  v=b²S−y,  u+v=c²S.

Then the three links u–x–y–v and the shortcut u–v are all
square-sum edges. The all-parity coefficient-shell forest is a
labeled Berggren-tree isomorphism via the common matrix addresses.
It does not follow that the image of arbitrary absorber operations
is a free tree on host edges.

Here deliberately hold c=18. The shell has the four ordered
nonboundary points

    (a,b)=(6,17),(10,15),(15,10),(17,6),

all with a²+b²=325. Start with the square-sum edge (1,3) of sum 4.
Use the first three coefficient choices for a ternary alphabet, or
all four for a quaternary alphabet. If a new pair (u,v) has sum S',
its parent sum is S'/324, and

    (a²−1) S'/324 < u < a² S'/324.

The integer inverse is exact: the ceiling of u/(S'/324) identifies
a², the alphabet identifies b, and then x=a²S'/324−u,
y=b²S'/324−v. Reject unless both are positive and the same logic
eventually reaches the prescribed ordered root. Distinct last letters
have disjoint open intervals. Thus even both *absolute labels*, not
just abstract matrix words, give a branch-preserving isomorphism with
the full rooted tree. This is a free iterated system of square-sum
absorber maps, **not** a claim that those three maps are the original
three Berggren matrices on the coefficient shell.

## Both absolute labels can be read through the count

Now use the [literal independent-address tournament bank](https://github.com/eliottcassidy2000/math/blob/2aa5cd874debced736c3461a65bd544f41d17f95/05-knowledge/results/chanpak_target_address_20261115/INVERSE.md)
at each host pair (x,y), with counts

    H_j(x,y)=432xy+129x+258y+15[x odd]+30[y odd]+45(p(x)+p(y))+135+2j,
    j=0,1,2,  p(0)=0, p(2u)=p(u)+u+[u odd], p(2u+1)=p(u).

All three choices are actual strong tournaments; they sit at spacing
+2 without selecting different addresses. The bound 0≤p(x)≤x follows
from the recurrence. A positive-depth address has S=4·324^d and
t=min(x,y) in the interval S/10<t<S/3: its smaller coefficient is
6 or 10, and a coordinate's normalized offset is strictly less
than 1/324, so 35/324<t/S<100/324.

For any fixed S≥1296, on the larger test wedge 1≤t≤floor(S/3) the
count windows are strictly sorted first by t and then by orientation
(the smaller count has (x,y)=(S−t,t)). This is an exact inequality,
not a monotonicity hypothesis for p:

- For fixed t the swap difference is
  129(S−2t)+15([S−t odd]−[t odd]). It exceeds 4, because
  S−2t≥S/3. The unordered checksum is identical in the two choices.
- The lowest product-plus-linear value at t+1 minus the highest at
  t is 303(S−2t−1). The total checksum change can cost at most
  45S, the parity flags at most 45, the window width at most 4.
  The resulting lower bound 56S−352 is positive for S≥1296.

Consequently different host vertices at the same depth have disjoint
absolute windows. For all positive depths each count lies strictly
between 35S² and 700S²: on [1/10,1/3] the normalized product is at
least 9/100, and the upper bound 108S²+303S+184≤700S² follows
from the same checksum bound. The next depth multiplies S² by 324²,
so these bands are disjoint, and the root is also outside them.
Each vertex of the tree therefore has a distinct literal tournament
window. The transport does not assert that its parent edge is an
arc flip or a local tournament insertion.

## Target-first recognition

Given N, select the unique possible depth from the separated bands.
The depth-zero triple can be tested directly. For positive depth,
j is determined by 2j≡N mod 3. Binary-search the literal count with
this j in the ordered wedge; either it misses or it recovers a
unique (x,y). Then run the integer absorber inverse. Accept precisely
if all recovered letters are in the alphabet and the final ordered
pair is (1,3). This certifies exact membership in the family without
enumerating levels: there are O(log N) depth and wedge tests and
O(log N) descent steps, all on O(log N)-bit integers. A count
evaluation uses O(log N) elementary integer operations.

For alphabet size m∈{3,4}, the number of attained counts ≤X has order

    Θ(X^(log(m)/(4 log(18)))),

by the disjoint depth bands and exact m^d vertices at level d.
So a perfect ternary inverse and exact +2 control are compatible
with an extremely sparse absolute spectrum. The missing ingredient
for Chan–Pak is an aperture that can vary with the target without
losing separation, or a way of coupling different apertures so their
absolute images fill gaps. Infinite cardinality or tree
isomorphism alone cannot supply this.

## Why mixing the apertures is not automatic

For a generic absorber write q=c², A=a², and use coordinates (S,x).
Its action is M_(q,A):(S,x)↦(qS,AS−x). Two such actions commute iff

    A/(q+1)=A'/(q'+1),

by comparing the lower coordinate of the two two-step composites.
Since q+1=a²+b² and coefficients are positive, this is equivalent
to equality of the ordered leg ratios a/b. In particular the
different symmetric Pell absorbers (5,5,7) and (29,29,41) commute.
Their two word orders give exactly the same absolute host pair from
every seed, not only the same hypotenuse or residue. Thus the prior
Berggren address on *coefficient triples* does not make the
variable-aperture operations on *hosts* a free tree. The fixed
aperture repairs the clash by replacing a numerical resemblance
with disjoint inverse intervals.

Replay: run `python verify_square.py`.

## A support barrier even when all square-sum apertures and seeds mix

Allow **every** positive address pair x,y with S=x+y a square, not
merely the preceding tree. Then the literal bank has at most
O(X^(3/4)) distinct counts ≤X in this whole class.

Indeed H_j≥432xy. For a fixed sum S, the smaller address t has
xy=t(S−t)≥tS/2, so the number of ordered pairs whose count could
be ≤X is O(min(S,X/S)). With S=c², sum c² up to c=X^(1/4) and
X/c² thereafter. Both sums are O(X^(3/4)). Multiplying by the
three center choices does not change the bound. This applies
to all square-sum absorber trees at once, whatever their roots,
coefficient sector or word-weight choices. It is a support bound,
not a statement that their algebraic side lacks exact inverses.
A decoder retaining the square-sum graph as a transport device
would have to change how its two labels are read, rather than
just enlarging the palette of square sums.

More generally, if admissible sums form a set A with
#(A∩[1,Y])=O(Y^alpha), 0≤alpha<1, then the same argument with
dyadic bands gives O(X^((1+alpha)/2)) counts. For bounded address ratio, ordinary
zero density of A also gives o(X): then S=O(sqrt(X)), so there are
o(sqrt(X)) sums and O(sqrt(X)) pairs per sum. Without that ratio
restriction the counting bound alone does not give o(X) from zero
density. Above sqrt(X), a dyadic band near Y contributes at most
X times its local density, and there are O(log X) such bands.
For instance, density on the order of 1/log Y in each band does
not make the sum tend to zero after division by X. This says
nothing about additional collisions in the unrestricted readout.
