# Bounded insertion gaps are finite; consecutive attained counts need reencoding

**Status: PROVED from cited OCF and Moon vertex-pancyclicity; FINITE-EXACT all-order gain-four classification. The unrestricted odd Hamiltonian-count spectrum remains OPEN. No external priority claim.**

Write H(T) for the number of directed Hamiltonian vertex words of a finite
tournament, with H(empty)=1. An induced extension changes no old arc. The
results below permit arbitrary added vertices and arbitrary orientations of
their new arcs; they do not assume a fixed block grammar, bounded feedback
set, or fixed source bank.

## 1. Inheritance and the precise new boundary

The inherited exact cut law is F=h+2K, with K a nonnegative weighted cut
capacity. Session124's `session124_tournament_cut_response_20261031/README.md`
already proves that a positive gain2K forces every crossing pair deletion to
have H at most K, and bounds the parent's triangle packing. Its equality
criterion identifies singleton-SCC insertion as the only zero-gain case.
It does not assert that every small numerical cut capacity is realizable.

The present consequence is stronger in a different coordinate: **fixing a
positive gain bounds the old Hamiltonian count itself, at all tournament
orders**, and hence prevents an eventual additive induction even after
arbitrary neutral padding or multi-vertex induced extension. The exact gains
two and four have only two possible endpoint pairs each.

Other load-bearing inherited reductions are not rediscovered here:

* `THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`
  already reduces global completeness to the strong p and7p lanes, with
  finite carry atoms49,63,343.
* `session124_million_path_spectrum_20261031/PROOF.md` already proves the
  arbitrary-block p-singular valuation bound p^max(0,r-2) dividing F.
* `session124_digital_constructor_inverse_20261031/PROOF.md` already gives
  exact inverse decision procedures, not universal positive answers.
* `square_sum_spectrum_overnight_20261031l/README.md` corrects the cyclic
  graph-product completeness claim. No such graph-transfer claim is used.

## 2. Two cited mechanisms and the local SCC budget

The odd-cycle collection formula (OCF) gives

    H(T) = sum over vertex-disjoint odd directed cycle collections Gamma
                 2^|Gamma|,

including the empty collection. Cycles have length at least three and are
identified by rotation, not by their vertex support. Dividing the collections
according to their unique cycle containing x gives the exact deletion law

    H(T)-H(T-x) = 2 sum_(odd directed cycles D containing x) H(T-V(D)).  (1)

This is inherited `THM-070-claim-a-from-ocf.md`. A current primary source for
OCF is Irving--Omar, *Revisiting the Redei--Berge Symmetric Functions via Matrix
Algebra*, [Corollary20](https://arxiv.org/html/2412.10572v3#S3.SS4), attributed
there to Grinberg--Stanley. Path reversal removes the complement convention.
Every summand is a positive integer, and induced deletion is monotone.

Moon, *On Subtournaments of a Tournament*, Canadian Mathematical Bulletin9
(1966), [Theorem1, p297](https://doi.org/10.4153/CMB-1966-038-7), says that
**each specified vertex** of a strong tournament of order m lies on directed
cycles of all lengths3,...,m. This is vertex-pancyclicity, not merely the
existence of one cycle of each length somewhere. The theorem number is1;
Theorem2 in that paper is a different strong-subtournament counting result.

Suppose an actual one-vertex deletion has

    F=H(T), h=H(T-x), F-h=2K>0.

Let C be the strongly connected component containing x, put m=|C|, and let a
be the product of H over every other SCC. SCCs are totally ordered by
dominance, so Hamiltonian paths concatenate and

    F=a H(C), h=a H(C-x), K=a K_C,
    K_C=(H(C)-H(C-x))/2.                                  (2)

If C were a singleton the gain would be zero. Hence m>=3 and K_C>=1.
In particular a is positive odd, a divides K, and a<=K. By Moon, x belongs
to at least floor((m-1)/2) different odd cycles, one at each odd length.
Applying (1) inside C, each such cycle contributes at least one, so

    K_C >= floor((m-1)/2),
    m <= 2K_C+2 = 2K/a+2.                                 (3)

This is an **affected-SCC bound**. The total tournament order is not bounded:
source/sink padding may be arbitrarily long, and other SCCs remain present.
The scalar factor a retains their entire multiplicative contribution.

### A canonical actual-pair kernel, not only a scalar bound

Let core(T) retain precisely the vertices in nonsingleton SCCs. If the outside
core is nonempty, then a>=3 and its order is at most a. Indeed every strong
block B has H(B)>=|B|, since a Hamilton cycle supplies |B| distinct rotated
Hamilton paths. For several outside blocks, their H values are at least3,
and their sum is at most their product a. If the outside core is empty,
a=1 and its order is zero. Therefore

    |core(T)| <= 2K+2.                                    (4a)

For a=1 this is (3). For a>=3, use K>=a and

    a+2K/a+2 <= 2K+2,

because 2K(1-1/a)-a >= a-2 >=1.

The same bound holds for **every** induced pair A subset T with positive
gain2K. Some x in V(T) minus V(A) must belong to core(T), since removing
only singleton SCCs preserves H. For that vertex, (1) and Moon give a positive
decrement, while monotonicity gives

    0 < H(T)-H(T-x) <= H(T)-H(A) = 2K.

Apply (4a) to this single deletion. Consequently the simultaneous reduction

    (A,T) -> (A intersect core(T), core(T))                 (4b)

preserves both Hamiltonian counts, the induced-inclusion relation, and a
strictly positive gain, and its larger tournament has order at most2K+2.
Every vertex discarded in (4b) was a singleton SCC of T and remains a
singleton SCC in every induced subgraph containing it, so deleting it from
A also preserves H. The padding sidecar consists of the removed singleton
vertices' positions among the ordered SCCs and their A-membership flags;
retain labels as well when reconstructing the original marked objects.

The order bound2K+2 is sufficient, not claimed sharp for every K.

This gives an exact logical finite decision/reconstruction procedure for a proposed
induced count pair with gap2K: enumerate every tournament T0 of order at
most2K+2 and every proper vertex subset A0, testing both path counts.
Acceptance gives an actual pair. Conversely every actual pair reduces to
one in this finite universe. No separate unproved spectrum-membership oracle
for the outside factor a is needed. The brute-force universe has size
2^{O(K^2)}; no efficient algorithm in log K is claimed.

This is not transitive closure of numerical one-vertex count pairs: two
edges with a common scalar H need not have compatible intermediate witnesses.
The kernel retains the actual T0 and the actual induced subset A0.

## 3. Uniform endpoint bound, including arbitrary induced extensions

For every tournament U of order r>=2,

    H(U) <= r!/2,                                         (4)

because reversing a vertex word pairs the r! permutations, and a word and
its reverse cannot both be directed paths in a tournament. Since |C-x|>=2,
(2)--(4) imply

    h <= a (2K/a+1)!/2 <= (2K+1)!/2.                      (5)

For the last inequality, a=1 is immediate. If a>=2, then K>=2 and
2K/a+1<=K+1<=2K. The factorial ratio
(2K+1)!/(2K/a+1)! contains the factor2K+1, which is at least a.
Thus the displayed inequality holds without a hidden growth assumption.
Consequently F<= (2K+1)!/2 + 2K.

**Theorem (bounded positive induced gap).** If A is an induced subtournament
of T and H(T)-H(A)=2K>0, then

    H(A) <= (2K+1)!/2,
    H(T) <= (2K+1)!/2+2K.                                 (6)

**Proof.** Add the omitted vertices in any order. By (1), each one-vertex
increment is nonnegative and even. Let the first strictly positive increment
be2k. Its old count is still H(A), and 1<=k<=K. Apply (5) at that step and
use factorial monotonicity. The second inequality follows by adding2K.
Zero-gain steps before or after that step cause no difficulty. QED.

Alternatively use the canonical pair (4b). Its smaller tournament has at
most2K+1 vertices, because its count is strictly smaller and it is a proper
induced subset. Apply (4), treating its possible orders0 and1 separately
with H=1. This gives the same bound and retains the entire actual pair.

Thus there are only finitely many raw-count endpoint pairs for any fixed
positive induced increment. This is not just an obstruction to a particular
seed bank. In gain notation delta=2K, the old endpoint bound is

    2 H(A) <= (delta+1)!.

As H(A) tends to infinity, this implies

    delta >= (1-o(1)) log H(A) / log log H(A).              (7)

For example log H(A)<= (delta+1)log(delta+1) suffices to invert the bound;
if delta is already at least log H(A), (7) is automatic. Otherwise
log(delta+1)<=log log H(A)+o(1), giving the claimed asymptotic. Natural
logarithms are used. No asymptotic assertion is needed for the exact (6).

**Nested-ray consequence.** Every chain of tournaments under induced
containment has only O(X loglog X/log X) distinct H values at most X, hence
zero density among the positive odd integers. Below sqrt(X) use the trivial
sqrt(X) bound. Above sqrt(X), successive distinct counts have gaps bounded
below by a constant times log X/loglog X by (7), giving the stated estimate.
Finite chains or rays whose H stays bounded are harmless. A finite union of
such chains also has density zero.

An equivalent local separation statement is useful: if
N>(2K+1)!/2, then any tournaments with **distinct** H values in
[N,N+2K] are pairwise incomparable under induced containment. Apply (6)
to a hypothetical comparable pair, whose half-gap is at most K. Equal H
values are excluded from this assertion because neutral padding nests them.

This does **not** bound the union of arbitrarily many branching chains.
Nor does it forbid arc flips, deleting/replacing old vertices, changing the
old induced tournament, or rebuilding an unrelated witness for the same h.
Those operations are precisely plausible ways around the obstruction.

## 4. Sharp classifications at gains two and four

**Gain two, all orders and all induced extensions:**

    (H(A),H(T)) is either (1,3) or (3,5).                  (8)

For K=1, (6) gives H(A)<=3. Its positive oddness leaves1 or3. Both occur:
a directed triangle has a transitive two-vertex deletion; the strong
four-vertex tournament with rows(8,1,3,6) has H5 and a triangle deletion.
This argument does not rely on the finite census.

In particular no witness of any attained odd h>=5 is an induced subtournament
of a witness for h+2. Merely proving h and h+2 separately attainable is not a
compatible insertion certificate.

**Gain four, all orders and all induced extensions:**

    (H(A),H(T)) is either (1,5) or (5,9).                  (9)

For one-vertex gain4, the external factor a is an odd divisor of2, hence1.
By (3) the affected strong SCC has at most six vertices. The complete finite
strong-tournament check in this packet finds only these two pairs. Both have
retained actual witnesses. For an arbitrary induced extension of total gain4,
the positive step sizes are either one gain4 or two gain2 steps. The latter
can only be1->3->5 by (8), giving the already listed pair(1,5). This proves
(9) at all orders from a proved cutoff, not extrapolation from small n.

The producer covers all33,867 labelled tournaments of orders1..6 by binary
arc assignments. Its separate strong presentation consumer fixes a directed
Hamilton cycle and enumerates all remaining arcs at orders3..6. Moon's theorem
(or Camion) guarantees that every strong class has such a presentation.
The latter counts all paths and all one-vertex deletions by literal
permutations, without importing the producer's DP or an old atlas.

The all-induced-subset audit also gives exact maximal canonical core orders:
4 at gain2, and5 at gain4. Thus the bound6 at gain4 is not sharp. At gain6,
C3 ordinal S4 gives a core of order7 and count15 with an induced count9; no
complete sharp gain6 classification is claimed here.

## 5. Nonvacuous infinite adjacent pairs, from actual source states

For every k>=6 there exist actual tournaments with both

    H=2^k+1 and H=2^k+3.                                  (10)

Nevertheless no tournament realizing the first value can be an induced
subtournament of any tournament realizing the second, by (8).

Here is a compact source for the second family. Take the parent rows

    (20,1,2,7,14)

on vertices0,...,4. An outgoing row is a bit mask. Adjoin apex x=5 with
out-neighbors{0,4}, so its cut mask is17. The full six-vertex rows are

    (20,33,34,39,14,17).

Literal path/cover counting gives the actual marked state

    (F,q,g,h)=(33,32,2,9),

where h counts parent paths, F counts full paths, g counts parent paths
ending in the apex cut, and q counts ordered two-path covers(P,Q), Q nonempty,
whose Q starts in that cut. P may be empty. This is an adjacency witness,
not a formally positive state assumed realizable.

Append a new vertex z dominated by all parent vertices and satisfying z->x.
Keep x marked. The inherited positive singleton action is

    (F,q,g,h) -> (F+q+g,2q,0,h).

For completeness, a full enlarged path either ends at z or contains z
immediately before x. The former case has F-(h-g) paths; splitting at z,x
in the latter gives q+h. Their sum is F+q+g. Each parent two-path cover
extends by placing z at the end of either path, giving2q; all parent paths
now end at z, outside the cut, giving g'=0 and h'=h.

After r>=1 zero steps the full count is

    2^r q+F+g-q = 2^(r+5)+3.

The result has6+r=k+1 vertices. The positive-tail condition matters:
the unextended source count33 is not2^5+3=35.

For the first family, take the single parent vertex a and apex x->a.
Its state is(F,q,g,h)=(1,1,1,1); k>=1 zero steps give2^k+1 on k+2 vertices.
`family_check.py` independently counts the source paths/covers literally,
reconstructs every tail edge and checks subset-DP counts at k=6,...,12.
The recurrences prove all k, not the finite tail tests.

The formulas are used as infinite witnesses for the new non-nesting
consequence, not as a claim of priority for these elementary count families.
They do not advance the existing finite global coverage ceiling.

## 6. Exact remaining proof obligation

The failed induction would infer an actual cut with K=1 merely from an
existing h witness. Equations (6) and (8) prohibit this at every h>=5, even
if a different h witness or arbitrarily many new vertices are chosen.
A bounded finite menu of positive additive increments also fails at all
sufficiently large h. The obstruction is an actual induced-embedding
predicate, not a numerical congruence obstruction.

A viable global construction must therefore use at least one of:

* unbounded increments, with an exact target-dependent source and cut;
* branching families whose union is not a finite union of nested rays;
* controlled reencoding or arc changes in an old witness; or
* an arithmetic factor descent supplying actual prime and7-prime atoms.

None of these alternatives is proved universally successful here. Global
attainment of all positive odd values other than7,21 remains OPEN.

## 7. Cheap hostile to an overstrong linear refinement

The small n<=6 data might suggest H(T)<=5K for every positive deletion
gain2K. That extrapolation is false. The strong seven-vertex tournament
with rows(14,4,8,18,39,79,31) has H27, while deleting vertex0 leaves H17.
Here K=5 and27>5K. Equivalently, the proposed uniform deletion ratio
H(T-x)/H(T)<=3/5 fails. The literal checker retains this witness; no claim
is made that every possible linear bound is refuted. The factorial theorem
above does not depend on this exploratory guess.
