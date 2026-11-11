# A low-ratio attachment cannot scramble neutral padding

**PROVED-IN-NOTE:** singleton-padding inversion bound, exact interval capacity,
the low-ratio normal form, and isolated-triangle rigidity below. **FINITE-EXACT:** companion enumeration.
**INHERITED:** positive odd-cycle expansion, attachment cut quadratic,
submodularity, and neutrality classification. **OPEN:** the universal two
branches and Chan–Pak. No canon ID or literature-priority claim.

## 1. Inheritance and the change of carrier

The closest mechanism is the positive odd-cycle formula (THM-002) and its
odd-path capacity descendant THM-4177. The usual quadratic attachment
formula and strict curvature are already in THM-4097, THM-4114, THM-4118,
and the explicitly credited section C of THM-5127. We do not rediscover them.
The hostile is THM-4114's order-five H9 parent whose extension image has
interior gaps despite full interaction rank. The corrected near miss is
MISTAKE-756: a uniformly contacted module can be reentered. The least-used
sidecar here is the **ordered list of singleton SCCs**, retained together
with the orientation-resolved even-vertex path capacity.

Concept board: (i) branch h -> 3h+2/4; (ii) positive directed capacities;
(iii) arbitrary neutral singleton padding; (iv) inversion partitions;
(v) switching Fourier boundary. The new pull is from (ii) to (iii):
nonnegative atom costs turn a bounded count ratio into a very small binary
language. Unlike a Fourier moment argument, no cancellation or support
inference is involved. Anchor: branch search; niche: ordinal padding;
wildcard: partition-number compression of binary attachment words.

## 2. The inherited exact positive carrier

Let B be a nonempty tournament, h=H(B), and let a new vertex v beat exactly
P subset V(B). For distinct i,j define the nonnegative integer

    c_B(i,j) = sum_Q H(B - V(Q)),

where Q runs over directed simple paths from i to j on a positive **even
number of vertices**. The empty tournament has H=1. Then

    H(B+v_P) = h + 2 sum_(i in P, j not in P) c_B(i,j).       (1)

This is the directed form behind the inherited symmetric curvature
mu(i,j)=2(c_B(i,j)+c_B(j,i)). To see (1) directly from the inherited odd-cycle
formula, separate cycle collections not using v (total h) from those using
v. Such a collection has a unique odd cycle through v. Remove v to obtain
an even-vertex path i -> ... -> j, where v beats i and j beats v; its cycle
weight is 2 and all disjoint remaining cycles sum to H(B-V(Q)). This is a
bijection and preserves every weight. Thus all terms in (1) are positive.

In particular if i->j is an arc, its one-edge path contributes

    c_B(i,j) >= H(B-{i,j}).                                (2)

This proof uses the odd-cycle identity as an inherited theorem, not as a
newly verified external claim. Formula (1) is not the switching energy:
switching can have signed interactions and higher boundary degrees;
attachment here is an oriented nonnegative cut with the source fixed.

## 3. New exact interval-capacity formula

Suppose u,w are singleton strong components in that order. Write

    B = L ▷ {u} ▷ M ▷ {w} ▷ R,

allowing any of L,M,R to be empty. Put m=|M|. Every path from u to w lies
entirely in {u} union V(M) union {w}; there are no paths from w to u.
Choosing its internal set S gives

    c_B(u,w)
      = H(L)H(R) sum_(S subset V(M), |S| even)
                         H(M[S]) H(M-S),                  (3)
    c_B(w,u)=0.

Indeed u followed by any Hamiltonian path of M[S], followed by w, is
exactly one admissible path. When S is empty use the one-edge path u->w.
Deleting this path leaves the displayed ordered product. This proof does
not assume that a path visits a nontrivial SCC only once after arbitrary
new arcs have been added; it uses paths in the **unextended** ordinal sum.

For each Hamiltonian path of M, cut after each even number of vertices,
including zero and including m when m is even. There are floor(m/2)+1
such cuts. Their ordered path pairs contribute to the sum in (3), and
concatenation recovers the original path and cut uniquely. Therefore

    c_B(u,w) >= (floor(m/2)+1) h.                          (4)

This remains true for empty M. In particular every forward pair of
singleton SCCs costs at least h; if two or more vertices lie between them,
it costs at least 2h. No unknown strong-tournament minimum is used.

## 4. New normal form below five times the source count

List all singleton SCCs in source order as z_1,...,z_s, and put
b_i=1 if z_i belongs to P. Let inv(b) be the number of pairs i<j with
b_i=1 and b_j=0. From (1), (2), and H(B-{z_i,z_j})=h,

    H(B+v_P) >= (1+2 inv(b)) h.                           (5)

**Theorem.** If H(B+v_P)<5h, then b has at most one inversion. If it has one,
the two vertices forming it are **adjacent singleton SCCs in the full
SCC order**, with no vertex or nontrivial component between them.

Proof. Two inversions contribute at least 2h to the cut in (1), contradicting
the strict ratio. If the unique inverted pair has m>=2 intervening vertices,
(4) gives the same contradiction. If m=1, that intervening vertex is itself
a singleton SCC. Its bit is either 0, creating an inversion with the first
endpoint, or 1, creating an inversion with the second. Either yields a
second inversion. Thus m=0. QED.

Equivalently, the singleton word is one of

    0^a 1^b,       or       0^a 1 0 1^b.                  (6)

In the second language the displayed 10 pair must be adjacent in the full
source, not merely in the singleton subsequence. These are necessary
conditions only. The nontrivial-SCC attachment bits still matter.

For the branch targets N=3h+2 and 3h+4, h>=3 implies N<5h, so the theorem
applies to both branches. In the one-inversion case the adjacent pair has
capacity exactly h by (3). The **entire remaining positive cut budget** is
therefore exactly 1 for 3h+2 and exactly 2 for 3h+4:

    sum_(i in P,j not in P, (i,j) != inverted pair) c_B(i,j)
       = 1 or 2.                                         (7)

This last condition is an exact, strong search filter. It is not a proof
that the residual budget can be realized. It says any atom of capacity
larger than 1 or 2 must not cross that cut.

### Stronger rigidity: an inversion below five times is exactly an isolated triangle

**Theorem.** Under N=H(B+v_P)<5h, if the singleton word has an inversion,
then the whole child is L ▷ C3 ▷ R, where
B=L ▷ {u} ▷ {w} ▷ R and the unique inversion is u in P, w not in P.
In particular N=3h. Conversely this isolated triangle always has count 3h.

Proof. The preceding theorem makes u,w adjacent. In fact for i in L,

    c_B(i,w) = H(R) U_L(i),
    U_L(i) = sum_(nonempty directed paths Q starting at i in L) H(L-V(Q)).

Every path from i to w stays in L until possibly visiting u. Its even
vertex parity forces exactly one of omitting or including u, according as
|Q| is odd or even. This proves the exact identity. Dually for j in R,

    c_B(u,j) = H(L) V_R(j),
    V_R(j) = sum_(nonempty directed paths Q ending at j in R) H(R-V(Q)).

Suppose some i in L is also in P. We claim U_L(i)>=H(L), and therefore
c_B(i,w)>=h. For every Hamiltonian path of L, split it
just before i into a prefix A (possibly empty) and a nonempty suffix Q
starting at i. If |Q| is odd, append w to Q. If |Q| is even, append u,w.
In either case this is an even-vertex directed path from i to w. Its
complement contains A, R, and possibly the singleton u, with count
H(A)H(R). More explicitly, summing the capacity terms over all suffix
sets Q and their paths includes every ordered pair (path on the prefix
set, path Q) obtained by this splitting, with weight H(R). Distinct
original paths give distinct ordered pairs, since concatenation recovers
the original. Therefore c_B(i,w)>=H(L)H(R)=h.
The already inverted pair has c_B(u,w)=h, so the cut would be at least 2h,
contradicting N<5h. Thus P contains no vertex of L. By reversing all arcs
and exchanging the roles of left/right, if any j in R were outside P then
c_B(u,j)>=h; this is likewise impossible. Hence the new vertex loses to
all L and beats all R. The child is exactly the asserted ordinal triangle,
with path count H(L)*3*H(R)=3h. QED.

**Branch consequence.** For every h>=3, a 3h+2 or 3h+4 attachment has
**no singleton inversion at all**. Its singleton word is 0^a 1^b.
The candidate count for these branches improves from (8) below to

    (s+1) 2^c,                                             (8a)

including s=0. The tiny residual budgets in (7) are therefore impossible,
not new existential subproblems. The earlier normal form remains useful
for other targets, and the rigidity theorem identifies its entire
one-inversion case exactly. Examining whether the residual budget can connect a core gives this
stronger classification.
The boundary inequality U_L(i)>=H(L) is the same path-cut injection as
Session250 `constructive/PROOF.md`, Lemma4; the new consequence is the
parity-resolved exact capacity and isolated-triangle rigidity.

### Corollary: one strong core needs only one active flank

Suppose B=TT_a ▷ C ▷ TT_b with C a nontrivial strong tournament, and the
child has count N<5H(C). If the SCC of the new vertex contains C, it cannot
also contain singleton vertices on **both** sides of C. In particular, if
the whole child is strong then a=0 or b=0.

Proof. If the child is strong and a,b>0, the new vertex must beat the
leftmost singleton: otherwise that vertex remains a source. Similarly the
rightmost singleton must beat the new vertex, or it remains a sink. These
two singleton bits are a forward inversion spanning C, contradicting the
normal form. More generally the SCC containing the new vertex is a union
of a consecutive interval of old SCCs: old SCCs cannot split upon adding a
vertex, and any two joined old SCCs force the intermediate ones to join.
Restrict to this interval. If it includes C and both flanks, the same
source/sink argument applies to its extreme singleton vertices, again
contradicting the original capacity bound. QED.

For N=3h+2 or 3h+4 with odd h=H(C)>=3, this gives an exact one-sided
reduction for this source grammar. Indeed gcd(h,N)=1. If C lay outside the
new vertex's SCC, its factor h would divide N, impossible. Thus C lies
inside, and every component outside is a singleton removable without
changing either count. The remaining witness is strong and has padding
on at most one side of C. This does not say that arbitrary source
realizers have only one nontrivial SCC.

This corollary connects the fixed-core compiler to the inversion theorem;
its proof uses only that theorem and the inherited multiplicativity of H
over ordered SCCs.

## 5. Algorithmic consequence and a finite-state missing lemma

Fix B, let c be the number of vertices in nontrivial SCCs, and let a(B)
count adjacent singleton-SCC pairs in the **full** source. For s>=1, there
are exactly s+1 zero-inversion singleton words and a(B) allowable
one-inversion words. Hence every attachment with child count <5h is among

    (s+1+a(B)) 2^c <= 2s 2^c                              (8)

candidates, instead of 2^(s+c). For s=0 retain the usual 2^c candidates.
This bound is linear in the number of neutral singleton vertices; only
the nontrivial-core vertex count remains in the exponent. It does not
bound the number of different sources.

More generally (5) gives inv(b)<=floor((N-h)/(2h)). For a fixed inversion
budget r, binary words of length s with at most r inversions number at most

    (s+1) sum_(k=0..r) p(k),                              (9)

where p(k) is the partition number. For fixed zero count, a word corresponds
to a partition inside the zero-by-one rectangle, of size inv(b); discard
the rectangle constraint and sum over at most s+1 zero counts. This is a
standard binary-inversion/partition correspondence, applied here as a new
padding compression consequence, not claimed as new partition theory.

**Exact remaining realization obligation.** For each realizable odd h>=95
and each e in {1,2}, find some source B with H(B)=h, a choice of its
nontrivial-SCC bits, and one of (6) such that the cut in (1) is h+e. The
sources may differ between e=1 and e=2. This is equivalent to the original
branch assertion, with a proved restricted padding language. By the stronger rigidity theorem, the unique-inversion state is impossible
for either branch: only the monotone state with full h+e budget remains. The finite language does not
close the unbounded collection of strong source components or prove the
required residual feasibility. No bounded universal transfer matrix is
asserted.

Connection contract: source is labelled padded tournaments with marked
attachment; target is their inversion-bounded singleton word plus exact
capacity cut. The map retains singleton order and all capacity entries.
It preserves the exact child count (1) and the low-ratio feasibility test.
Keeping only the word loses the strong-core attachment and residual
capacity; (7) or the full matrix is the mandatory sidecar. The cheapest
positive next experiment is separate searches for monotone and one-defect
witnesses for the known h>=95 branches, using only the monotone state for either branch, rather than claiming
interval support from submodularity.

## 6. Failure boundary and exact tests

The strict ratio is essential: TT4 with ascending orientation and beaten
set consisting of its second vertex has two singleton inversions and
child H=5, source H=1. The one-inversion state is not vacuous for nontrivial
h: an ordinal C3 together with an adjacent singleton pair, with the new
vertex completing a triangle on that pair and meeting the C3 uniformly,
has child 3h=9. These examples are constructed directly, not inferred from
moments. For arbitrary large h, placing a 3-cycle after any source gives
the same 3h equality mechanism whenever h is realizable.

Run from the repository root:

    python3 -B 05-knowledge/results/session251_creative_completion_20261111/geometry/check_padding.py

The script exhausts all labelled sources of orders 1 through 5 and all
attachments; there are no isomorphism or strongness filters. Integer subset
DP counts the child independently of path-capacity enumeration. It checks
(1), the new interval lower bound (4), inversion bound (5), the full-order
adjacency conclusion, branch cases, isolated-triangle rigidity, the one-sided-core corollary,
and the strict-threshold hostile. The
saved output is the exact finite scope, not an all-order computational proof.
No inherited large census is rerun. This is verification of new consequences.

## 7. Source map

- `01-canon/theorems/THM-002-ocf.md`: inherited odd-cycle expansion.
- THM-4177, `01-canon/theorems/THM-4177-root-split-johnson-current-and-source-sink-boundary-census.md`:
  inherited positive even-vertex path capacities and orientation sidecar.
- THM-1900, insertion-response calculus: inherited neutral insertions are
  condensation down-sets; repaired proof warns against deleting an apex
  from a path and assuming the result stays a path.
- THM-5127, heat-flow/vertex-cube theorem: explicitly credits the earlier
  attachment quadratic, submodularity, and strict pair curvature.
- `05-knowledge/results/session250_sixteen_20261110/geometry/SWITCHING_BOUNDARY_RECURRENCE.md`
  and `ISING_INVERSE_INTERFACE.md`: switching boundary carrier, contrasted
  with the positive directed attachment carrier here.
- `05-knowledge/results/decomposition16_odd_spectrum_20261110/README.md`:
  inherited threshold95 branch reduction and permission for different
  source realizers. This note strengthens padding geometry, not that
  unresolved existence quantifier.

No new general research-method card is claimed from this single bridge.
