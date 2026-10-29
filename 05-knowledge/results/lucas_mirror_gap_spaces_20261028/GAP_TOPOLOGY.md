# Threshold gaps realize an independence complex, with the marked-target boundary

**PROVED** elementary equivariant homotopy equivalence below. Cycle homotopy
types are **CITED**, not new. Rational controls are **FINITE-EXACT** and the
audit is same-agent. This is a static configuration-space theorem, not an
LRC trajectory theorem or a proof of LRC(14).

## 1. A general threshold-space theorem

Let K be a nonempty simplicial complex on a subset of[n], with[n] not a face.
Let a_i>0 and sum a_i=1. In the standard closed simplex Delta define

    B_K(a)={g in Delta : S_a(g)={i:g_i>=a_i} is a face of K}.

The weak inequality is intentional. Define its geometric realization by

    |K|={u in Delta : support(u) is a face of K}.

**Theorem G1.** B_K(a) is homotopy equivalent to |K|. The same is true when
all gap coordinates g_i are required to be strictly positive. The maps and
homotopies are equivariant for every coordinate-permutation group preserving
K and the vector a.

**Explicit maps.** For g in B_K(a), put

    e_i=max(g_i-a_i,0), E=sum e_i,       r(g)=e/E.

E>0: otherwise g_i<=a_i for all i and equal sums force g=a, whose safe set
is[n], excluded by hypothesis. The support of r(g) lies inside S_a(g), so
r(g) belongs to |K|. For u in |K| put

    s(u)=(a+u)/2.

Its safe set is {i:u_i>=a_i}, a subset of support(u), so s(u) belongs to
B_K(a), with strictly positive coordinates.

**Homotopies.** The segment

    H(g,t)=(1-t)g+t*s(r(g))

stays in B_K(a). If i is outside S_a(g), then g_i<a_i, r(g)_i=0, and
s(r(g))_i=a_i/2<a_i. Thus the safe set along the segment is a subset of
the original face S_a(g). The other homotopy is the segment in |K| from
u to r(s(u)); its support is contained in support(u), hence stays inside
that simplex of K. These prove s r homotopic to id and r s homotopic to id.
Continuity and the claimed equivariance follow from the formulas. The first
homotopy also retains strict positivity if it was present initially. QED.

This is a homotopy equivalence, not an equality of safe masks or a claim that
s is a literal inverse. In particular, zero-excess coordinates at equality
may disappear under r. The proof treats that boundary explicitly.

## 2. Apply the theorem to the inherited LRC gap criterion

For n>=3 points in a fixed cyclic order on the unit circle, let g_i be their
clockwise gaps, summing to1. A vertex is at distance at least1/n from all
other vertices iff both incident gaps are at least1/n. This is inherited
[THM-370, LRC two-neighbor cycle independence](../../../01-canon/theorems/THM-370-lrc-two-neighbor-cycle-independence.md).

Consequently the static space B_n of configurations in which **no vertex is
lonely** has safe sets that are nonempty independent sets of C_n. The empty
set is impossible because every gap being less than1/n contradicts the sum.
Every nonempty independent S actually occurs. For k=|S| choose

    g_i=(n+k)/(2nk) if i in S,       g_i=1/(2n) otherwise.

These positive rational gaps sum to1 and have exactly S as safe set.
Theorem G1 gives the stronger continuous statement

    B_n is homotopy equivalent to Ind(C_n),

equivariantly under the dihedral permutations of the gap positions.

Kozlov's classical cycle calculation, reproved as Theorem4.3 in
[Barmak, *Star clusters in independence complexes of graphs*](https://arxiv.org/pdf/1007.0418),
therefore gives

    B_(3k)   ~ S^(k-1) wedge S^(k-1),
    B_(3k+1) ~ S^(k-1),
    B_(3k+2) ~ S^k.

In particular B5~S1 and B14~S4. These concern static, cyclically ordered
all-vertex failure configurations; they do not assert anything about every
integer-speed trajectory entering or leaving them.

At n=5 the ten nonempty safe patterns are precisely the nonzero Lucas words.
Under single-bit adjacency their graph is C10, the outer graph of the golden
eleven-point realization. The barycentric subdivision of Ind(C5)=C5 is also
this C10: its vertices are the five singleton faces and five two-element
faces, and its edges are their inclusions. This gives an actual map, not
merely two objects of size10. The zero Lucas word supplies the cone point
for a cubical enlargement but has no gap-source counterpart at threshold1/5.

## 3. A specified runner changes the topology completely

Keep one distinguished runner, with incident gaps indexed0 and1 after naming
the cyclic frame. Its failure condition is only

    not(g_0>=1/n and g_1>=1/n).

The allowed safe faces form K_* = Ind(K2 disjoint (n-2) isolates). This is
S0 joined with an(n-3)-simplex. For n>=3 it is a cone and is contractible.
Theorem G1 proves that this marked-runner failure space is contractible too.

Thus the S4 in the fourteen-gap all-runner problem cannot be used as though
it were the failure topology for one chosen runner. The first failed
implication would be replacing a union of distinct success predicates by a
single marked predicate. This is a target-loss obstruction, not a refutation
of topological approaches to LRC.

The dihedral action here is a **spatial relabelling of gaps**. It differs
from time reflection t->1-t for a fixed velocity list. The latter's exact
trace obstruction is inherited from
[odd-sign spectra](../odd_sign_spectra_trijections_20261028/LRC_UNIT_DISTANCE.md)
and is strengthened by the concurrent
[reflection persistence theorem](../session91_packing_golden_carriers_20261029/geometry/LRC_REFLECTION_AND_GOLDEN_GAP.md):
mixed-parity time-reflection traces vanish across the whole positive-threshold
persistence module. One cannot substitute one reflection action for the other.

## 4. Why the six-period sign appears, and what it does not measure

For n>=3 the independent-set polynomial of C_n is

    I_n(t)=trace([[1,t],[1,0]]^n).

This follows by multiplying the two-state cyclic transfer matrix, with a
factor t for each occupied position. At t=1 its trace is L_n. At t=-1 its
eigenvalues are exp(+-i pi/3), giving

    I_n(-1)=2 cos(n pi/3),       reduced_chi(B_n)=-I_n(-1).

This proves a genuine six-period relation between Lucas digit parity and the
Euler characteristic of the gap space. The same recurrence gives L_n even
exactly when3 divides n. Therefore odd word length, odd number of words,
odd number of occupied positions, and odd tournament automorphism order are
different conditions. In particular W3 has four words.

The cubical complex made of coordinate cubes for independent supports is
contractible: every cube contains zero and radial contraction stays in it.
Its radial simplex section sum x_i=1 is Ind(C_n). A contractible cubical
carrier can thus have a spherical or wedge-of-spheres link. Adding or deleting
the common zero point changes the relevant topology; dimension alone records
none of that distinction.

## 5. Productive next obligation

An arithmetic LRC application must now specify the map from a fixed labelled
integer-velocity instance into this gap space through collision times,
retain the marked runner while cyclic order changes, and prove an obstruction
to remaining in its own failure set. The equivalence above alone supplies
none of those three assertions. A concrete next experiment is a relative
chain complex of the different marked success sets with collision transport,
tested first on the known AP13 boundary and the mixed-parity empty/point/
interval controls. A nonzero global Euler characteristic is insufficient.

The late incoming [joint-port/germ packet](../joint_port_germs_20261028/README.md)
sharpens this next target: its positive threshold-growth detector distinguishes
isolated witnesses from emptiness, but the rows{1,3,6} and{2,3,6} at1/3 have
matching local germs and owner types while adding speed7 kills only the first.
Their surviving phases differ. That is a concrete reason for retaining actual
phase in the proposed chain model, beyond merely asking for a positive statistic.

`check_bridges.py` tests330 rational cyclic-gap sources and their selected
homotopy points; additional tests use nonuniform thresholds, equality cases,
zero gaps and the marked-target complex. Those checks supplement the
all-parameter proof of G1. They do not establish a velocity-bound or convergence.

## 6. Fixed topology can hide exponentially growing return information

The concurrent [guarded-flip theorem](../session91_packing_golden_carriers_20261029/operators/GUARDED_FLIP_LEAKAGE.md)
computes the forbidden-return Gram D=PB(1-P)BP, where B sums named bit flips
and P projects onto independent supports. Its rank depends on the legal
operator, not just on the homotopy type of Ind(H).

Here is an explicit all-size separation, also proved directly. Take the
constraint graph H=K_(1,m), m>=1. Its independence complex is the disjoint
union of a point(the centre) and the simplex of all leaves; hence it is
homotopy equivalent to S0 for every m. Theorem G1 gives the same fixed
homotopy type for all associated threshold spaces.

Nevertheless

    rank_Q D = 2^m-1.

Indeed the empty state contributes zero. The singleton exchange block is
the signless Laplacian of the star, of rank m: its quadratic form is the
sum over leaves of(x_c+x_i)², with a one-dimensional kernel. Every leaf subset
of size at least2 has one forbidden insertion(the centre), with no other
legal predecessor of that forbidden state. It therefore contributes an
independent diagonal pin. There are2^m-m-1 such subsets. Summing the ranks
gives2^m-1. The checker constructs the full native Gram and row-reduces it
over Q for1<=m<=5.

This is a uniform obstruction to bounding an exact real linear return
interface by the homotopy type of the threshold space. The topology stays
S0 while the required Gram-factor rank grows exponentially. It does not
preclude succinct symbolic formulas or nonlinear encodings. In particular,
the S4 description at fourteen gaps does not supply a small source-weighted
Gram controller merely by being a low-dimensional homotopy model.
