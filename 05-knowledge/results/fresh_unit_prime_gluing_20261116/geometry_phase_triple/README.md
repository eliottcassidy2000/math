# Rotating the contacts: geometry, square-sum endpoints, and a new exact window

The prescribed 6015/6017/6019 continuation and Chan–Pak remain open.
There is a three-block construction for a different target, a reusable
formula for its phases, and a finite reduction for the prescribed triple.

## Geometry: keep the face permutation

The [5–6–7 flag-carrier session](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/geometry_567_flag_carriers_20261004.md)
works with edge reversal α and local cyclic rotation σ on directed edge ends.
Faces are the cycles of σ∘α. For a cubic cellular map on an orientable
closed surface, Euler's equation gives
```
6χ = Σ_faces (6 − face length).
```
Thus with only pentagons, hexagons and heptagons, F5−F7 = 12(1−genus).
The hexagon contributes zero, and pentagons and heptagons have equal
opposite charges even in a mixed map. The equivelar versions are spherical,
Euclidean and hyperbolic respectively (the last for the regular
{7,3} metric; see
[the classical tessellation criterion](https://www.d.umn.edu/~ddunham/dunbrid00.pdf)).
The archived hostile example is the same Heawood graph with the same
Paley order-21 action: equal point/line rotation signs give 7 hexagons
and genus one; opposite signs give 3 fourteen-gons and genus three.
A feasible charge balance is not a map; the local rotations matter.

## Eight and nine: the cut of a reflection orbit

The archived [THM-4510](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md)
can be read directly on the cyclic set Z/16. Take R(x)=−x and S(x)=9−x.
Their product S∘R is translation by 9, so it is transitive. The
loop-erased two-matching graph is connected, of maximum degree two,
with ends exactly the two R-fixed points 0 and 8. There are no S-fixed
points and no common edges. Removing 0 exposes its S-neighbor 9.
The three integer sums are 9, 16 and 25, the squares from (3,4,5):
8 is a surviving reflection fixed point and 9 is an exposed contact.
These are two different roles, not a curvature pair.

More generally, replace 16 and 9 by even B and odd A with 0<A<B and
gcd(A,B)=1. The three target sums A,B,A+B give a path on [1,B−1]
from B/2 to A, and B attaches at A. For the same three-sum template,
B+1 has just one neighbor A−1. It extends the unique path on [1,B]
exactly when A−1=B/2, because deleting that leaf from any new
Hamilton path would have to leave the old path. At 9 and 16 this hits
the other end exactly. Vertex 18 attaches to 7 instead; the three
leaves 16,17,18 preclude a Hamilton path. The full square-sum graphs at 15–18 also contain the chord 1–3 of
sum 4. It changes no leaves. At 15–17 the other edges are the
spanning path just described; using the chord would require dropping
one path edge incident with both chord endpoints, which is impossible.
Thus these paths are still unique up to reversal in the full graph. In the
Pythagorean subfamily A=a², B=b²=4k², the second-attachment test
is the simultaneous Pell condition a²−2k²=c²−6k²=1. This is an
exact arithmetic translation, separate from the tournament construction.

## Geometry of a fixed tournament deletion

Let G be the root-deletion tournament, g=H(G), and let x_i in {0,1}
mark root outcontacts. For i≠j let M_ij count even-order paths in G
from i to j with weight H(G minus the selected path). Write W_ij=M_ij+M_ji,
and r_i=Σ_j(M_ij−M_ji). The even-path cut identity gives
```
H−g = 2Σ_ij M_ij x_i(1−x_j)
    = Σ_{i<j} W_ij(x_i−x_j)^2 + r·x .
```
E and O are affine-linear functions of x by first successor.
If L is the symmetric weighted Laplacian and Lz=−r/2, the latter
formula is H=g−zᵀLz+(x−z)ᵀL(x−z). Such z exists on the zero-sum
subspace because r sums to zero and W_ij≥1 (use the direct arc).
A fixed continuation is a Boolean cube slice by two affine equations;
the requested H is an exact ellipsoid level on that slice. Continuous
intersection does not imply a Boolean point. This is a Euclidean
geometry statement, not a claim that tournaments inherit the sign
of tessellation curvature.

### The three cyclic blocks

If G has a rotation preserving all arcs and cycling each of several
orbits, the coefficients of E and O are constant on each orbit.
Thus independently rotating a fixed contact subset inside each orbit
keeps (E,O,g) fixed. The energy still sees relative phases:
the difference across orbits i,j is a function of p_j−p_i, since
simultaneously rotating both preserves W. Intra-orbit terms and the
drift are constant. After removing simultaneous rotation, there are
two relative phases for three cyclic blocks. This is the analogue of
keeping the relative rotation in the Heawood figure.

Use outgoing row masks
```
G = (74,148,289,22,37,11,190,317,123)
```
on indices 0..8. The rotation simultaneously increments indices mod
three within the blocks (0,1,2), (3,4,5), (6,7,8). Attach root 9:
select one contact at index i in the first block, omit index 3+j in
the second block, and select only 6 in the third. The resulting table is
```
          j=0   j=1   j=2
 i=0      6875  6815  6923
 i=1      7287  7175  7201
 i=2      7289  7311  7285
```
Every entry has (E,O,g)=(5586,6357,1521), with 771 Hamilton paths
ending at the root. Hence 7285,7287,7289 share the same physical deletion
and the entire reader signature. Their contact masks are 92,114,116.
From phase (2,2), move (−1,+1) then (+1,0) mod 3: both steps give +2.
Reversing order takes the intermediate count 6923; the mixed difference
is 364. The phase moves commute as relabelings of contact sets but their
deposits cannot be assigned independently. This is an exact example of
why keeping the rotation matters.

An arbitrary length m≥1 of near-sinks (all nonroot old vertices beat
the new block, and the new transitive block beats the root) gives
```
H = 10422·2^m − 2387 + d,       d=0,2,4,
E = 10422·2^(m−1),              O = E+1521,        g=1521.
```
Here the inherited [terminal suffix law](../real_boundary_signatures/COMPACTIFICATION.md)
reads Z=E+O−g=10422 and I=H−2O+2g=−2387+d. This is an unbounded family
of exact consecutive windows, not overlap at every prescribed baseline;
it cannot produce the height-three target 6015.

## What remains at 6015 / 6017 / 6019

There is now an exact finite decision rule for the whole triple with no
assumed ceiling on the input tournament size. Write X=E+O. For existence
of a common continuation, it suffices and is necessary that their cores
share a sink key (g,oddpart X), since even a common positive O−E can be
turned into a common zero by adding ordinary sinks. Deleting ordinal
singletons outside the root SCC preserves H,g and oddpart X. A suffix
singleton halves X and a prefix singleton leaves it unchanged. Root
singleton cases have g=H, so cannot be used by all three. For each
factorization of H into SCC counts, use [Busch's bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
on the individual nontrivial components. Thus the reduced orders
are at most 19, 17 and 17 respectively (not an enumeration through those
orders). For example the largest bound for 6015 is attained by the
factor list 3,5,401 with order caps 3,4,12.

For a proposed sink key one can sharpen by treating the unknown
root SCC separately. If A→C→B and the root is in C, then H and g of C
scale by H(A)H(B); X scales by H(A)P₂(B), with P₂ the ordered two-slot
cover count including empty slots. These formulas follow by splitting
each path at the ordinal seams; choice of the two B slots flips the
selected parity in aggregate when B is nonempty (by the endpoint identity
a root in C ends no spanning path of the total). In particular:
- To close the existing key (729,12217), the missing 6019 necessarily
  has one nontrivial component: gcd(6019,729)=1.
- At (1185,8051), a 6015 witness cannot have a nontrivial prefix
  (gcd of all three numbers is 1). The only possible nontrivial suffix
  has H=5, P₂=32. Suffix H=3 and H=15 are ruled out by their odd P₂
  factors: for strong components respectively {3} and {15,27},
  and a decomposable H=15 has a triangle factor. Hence the alternative
  to a single core is exactly the reduced key H=1203, g=237,
  oddpart X=8051, with root SCC order ≤14. The small suffix table is
  implemented by suffix_spectrum.cpp. **The reduced key is impossible
  at all orders:** [complete finite certificate and reduction](REDUCED_237.md).
  This eliminates the nontrivial suffix route at (1185,8051), but does
  not eliminate a single strong 6015 core or a different key.
- There is an additional congruence gate in the **deletion**, not an
  ordinal suffix of the full rooted host. Appending one deletion component
  of path count h and two-slot count R sends the previous (g,X) to
  (g h, R X+g(h−R+L)), where L is the sum of selected start weights in
  this component for the contacts there. This formula partitions two-slot
  selections by whether their chosen path visits the previous deletion,
  starts in the new component, or is empty. If δ is the gcd of h and
  all the single-contact start weights, then δ also divides R because
  their sum is R−h. Induction shows the product of the δ's divides
  gcd(g,X). In particular every triangle deletion SCC contributes a
  factor 3 (h=3, R=12, each start weight 3). Neither given mantissa
  is divisible by 3, so neither key can have even one literal triangle
  deletion SCC. In the reduced g=237=3·79 case this leaves a single
  nontrivial deletion SCC of count 237, rather than a 3-by-79 split.
  It does not eliminate that remaining component.
- Every root-preserving tournament automorphism acts freely on the
  objects counted by X: an ordered selected root path and an ordered
  complementary Hamilton path fix every vertex if the object is fixed.
  Its action on the deletion Hamilton paths is also free. Thus the group
  order divides gcd(g,X), forcing rooted rigidity for both known keys.
  Any unmarked deletion symmetry must be broken by the two-coloring of
  vertices into outgoing and incoming contacts.
- Count-neutral induced extension cannot move a witness to another
  sink key. Indeed any insertion at positive cut weight strictly
  increases H. Zero cut weight means every incontact beats every
  outcontact (otherwise the offending two-vertex path has positive
  weight), so the added vertex is an ordinal singleton at an SCC cut.
  Iterating proves the assertion for arbitrary induced neutral extension.

The all-contact experiments are scoped. The factor sampler (one million
seeded draws) yielded four pair keys but no triple among its ordinal
schemes. The rotational scan (150,000 seeded draws with repetitions)
yielded 20 sink keys for the target values but no pair. The exact
same-deletion chart scan enumerates all orientations subject to the
displayed 2- or 3-cycle-block template at deletion orders 6–11 (cycles chosen with the same reference orientation); its
first same-signature, consecutive triple is the one above. It does
not characterize all tournaments of those orders. The separate
diag_windows.cpp enumerates all labeled invariant deletions at orders 9
and 10 under the exact simultaneous cycle (and a fixed tenth point),
allowing both orientations of each internal cycle. With g≤6019 the
complete center sets are {7287,8981} at order nine and
{7287,8981,18459,21071,21609,24221} at order ten. Thus 6017 is absent as
a same-deletion center in those symmetry classes, not merely missed at
one parameter choice. The unrestricted larger-order triple is not excluded. 

The proof obligation for Chan–Pak is still a *transversal at prescribed
absolute levels*: sufficiently many actual deletion shells, with adjustable
centers and weights, must cover all admissible targets, not merely give
long runs whose base cannot be chosen. The relative-rotation construction
supplies exact local deposits; it does not discharge that placement
obligation.

The independent unrestricted seeded searches (orders 9/10/11, plus
a focused run against the two known sink keys) combined to 146 sink
keys containing two of the old targets but zero containing all three.
This count is empirical and not needed for the theorem or any exclusion.

## Replay

- `python3 check_phase.py`: independent literal path DP, nine phase profiles,
  square-sum endpoints, and the first near-sink window.
- `c++ -O3 -std=c++17 same_deletion.cpp -o /tmp/same && /tmp/same 11`:
  the complete scoped chart scan.
- `c++ -O3 -std=c++17 suffix_spectrum.cpp -o /tmp/suffix && /tmp/suffix`:
  tiny exhaustive factor table.
- `c++ -O3 -std=c++17 factor_search.cpp -o /tmp/factors && /tmp/factors 1000000`:
  seeded factor search.
- `c++ -O3 -std=c++17 rotation_search.cpp -o /tmp/rotate && /tmp/rotate 150000`:
  seeded rotational search.
- `c++ -O3 -std=c++17 diag_windows.cpp -o /tmp/diag && /tmp/diag 9 && /tmp/diag 10`:
  exact diagonal C3 center sets with g≤6019.
