# Which lattice and spectral identifications preserve Hamilton paths?

**PROVED:** the effective root-frame realization criterion below and the
stated algebraic consequences of conference switching.
**FINITE-EXACT:** the twelve-core family and retained controls.
The spectrum assertion “all positive odds except7,21 occur” remains **OPEN**.

## Recovering the actual integral bridge

The relevant prior source is
[Session36, METRIC_CURL_RECOVERY, sections1-5](../session36_information_curl_contexts_20261023/lattice/METRIC_CURL_RECOVERY.md),
with its hostile and lineage in
[RECOVERY_AND_ROUNDS](../session36_information_curl_contexts_20261023/lattice/RECOVERY_AND_ROUNDS.md).
For a tournament with skew matrix S, retain the integral augmentation lattice

    L=A_(n-1)={x in Z^n:sum x_i=0},
    g(x,y)=x^t y,       omega(x,y)=x^t S y.

In the labelled star basis b_i=e_i-e_r, r=n-1, its Gram matrix is I+J and
its alternating matrix is

    K_ij=S_ij-S_ir+S_jr.                               (8)

The prior fixed-label fibre theorem proves that identical K permits
exactly cyclic rotations of the ordered strong-component list. Such
rotations preserve H, since ordered joins multiply Hamilton counts.
In particular the map is injective on strongly connected labelled
tournaments. With integral metric isometries retained, it classifies
strong tournaments up to isomorphism. This is a genuine H-preserving
structural identification, not a numerical resemblance of spectra.

### An explicit inverse: tournament feasibility is 2-SAT in a star frame

**Theorem R1.** A skew integer matrix K of size n-1 with zero diagonal
is realizable as (8) iff there are signs t_i in {-1,+1} satisfying

    K_ij+t_i-t_j in {-1,+1} for every i!=j.            (9)

The decision and a realization use a 2-CNF formula of quadratic size,
and hence a quadratic-time strongly-connected-component algorithm.

**Proof.** If a realization exists put t_i=S_ir, giving (9).
Conversely define S_ir=t_i,S_ri=-t_i, and
S_ij=K_ij+t_i-t_j. These entries give a tournament and recover K exactly.

Each pair admits only four sign assignments. Forbid each invalid assignment
by a two-literal clause. Entries outside {-3,-1,1,3} are immediately
impossible. The implication graph has2(n-1) literals and O(n^2) arcs;
the usual mutual-reachability criterion and component order produce an
assignment or an inconsistency. QED.

The literal constraints have a useful interpretation:
K_ij=3 forces t_i=-1,t_j=+1; K_ij=-3 forces the reverse signs.
K_ij=1 forbids(+1,-1), and K_ij=-1 forbids(-1,+1).

For example the correctly typed skew form

    K=((0,3,3),(-3,0,3),(-3,-3,0))

is infeasible: its01 pair forces t_1=+1 and its12 pair forces t_1=-1.
Thus arbitrary integral alternating forms with plausible entry sizes do
not automatically represent tournaments. The code checks the inverse and
H preservation on all1099 labelled tournaments through5 and on the prime
twelve controls.

The new contribution here is the explicit feasible-form compiler attached
to the recovered fibre theorem. No literature priority is claimed.

## The coarser conference identification loses H

[Session35 CONFERENCE_COLLAPSE](../session35_dimensions_dualities_20261023/zeta/CONFERENCE_COLLAPSE.md)
already proves the general conference algebra and the order-eight H split.
Session36 already records order-twelve masks0,1 with H95095 and196815 and
the same selected genus-three coefficient660. Those discoveries are
credited here. The present calculation extends to the full specified
2048-member orbit, its prime strata, full named response, and simple-path
support sizes.

For the matrix S in [D4](PROOF.md), exact multiplication gives S^2=-11I.
Every diagonal sign switch S_X=D_X S D_X has the same identity. Hence

    I-S_X^2=12I,
    det(tI-S_X)=(t^2+11)^6,
    1^t(tI-S_X)^(-1)1=12t/(t^2+11).                   (10)

Thus the entire positive Gram lift, scalar theta/Epstein data attached to
it, and rooted skew resolvent are unchanged. This is inherited conference
algebra, not a new spectral theorem.

Put A_X=(J-I+S_X)/2. The determinant lemma and resolvent identity give

    det(tI-A_X)=t(t-5)(t^2+t+3)^5,
    1^t(tI-A_X)^(-1)1=(12t+6)/(t(t-5)).                (11)

Expanding at infinity proves that the total directed-walk counts are

    W_0=12,       W_k=66*5^(k-1) for k>=1.             (12)

Equations(10)-(12) hold across the full switching orbit, while H takes
four distinct values. In particular masks1 and3 are both strong and prime,
with H196815 and236447. Their complete scalar linear walk observers are
identical but their Hamilton counts differ. An orthogonal or spectral
identification here is not an isomorphism of the named tournament
adjacency relation.

The exact support-size profiles locate the first lost information:

| Switch mask | Simple paths on1..5 vertices | On6 vertices | On12 vertices |
|---:|---|---:|---:|
|0|12,66,330,1485,5940|20460|95095|
|1|12,66,330,1485,5940|20700|196815|
|3|12,66,330,1485,5940|20796|236447|
|11|12,66,330,1485,5940|20844|256023|

The full2048-member certificate has exactly one value at each size1..5 and
four at every size6..12. Walks permit repeat visits; these counts retain
the actual visited set. This is the concrete coordinate lost by the
linear-system analogy.

## Root-shell observations versus the full named cubic

Session36's general conference theorem makes the **complete two-root shell**
distribution of metric and alternating pairings identical at fixed order.
This does not identify the whole genus-two decorated theta series.

Its selected genus-three coefficient has
Gamma=((2,-1,0),(-1,2,0),(0,0,2)),
and alternating entries Omega_12=-3,Omega_13=-4,Omega_23=2.
Independent root-tuple and five-vertex-incidence implementations give660
for each of masks0,1,3,11. Consequently even the three strong prime
controls remain H-distinct after adding that specific scalar observation.
No assertion is made that their full genus-three data agree.

The full named cubic two-port response is a different object. For example,
at(a,c)=(e_0,0), mask1 and its converse give691765 and997475;
mask3 and its converse give874729 and1341291.
These are actual response differences, not root-shell coefficients.
The common number “three” does not identify the two observation spaces.
The full integral alternating form remains the H-preserving sidecar;
its recoverability from two-port data at prime order12 remains unproved.

## Retained three-state inflation and its limits

Incoming [Session186 G1-G3](../session186_proof_frontiers_20261104/gate/PROOF.md)
and the [homogeneous-pair packet](../homogeneous_pair_response_20261103/PROOF.md)
retain the three strong six-vertex seeds

    Q0=(16,5,33,39,46,3), Q1=(16,5,49,39,42,3),
    Q2=(18,4,49,39,42,3).

Inflating name5 by TT_r gives three strong tournaments linked by single
old-edge reversals, with the three-point odd interval

    H(T_j(r))=4*3^r+9*2^r-3+2j,    j=0,1,2, r>=1.

The genuine3^r and2^r terms come from the exact occupied-run functional

    H(Q_j[M])=(27+2j)H(M)+21q_2(M)+4q_3(M).

It retains ordered path-cover numbers q_k, yet forgets their named endpoints.
Converse preserves each q_k. Nonuniform two-port attachments restore the
lost information: for consecutive names a=5,b=6,w=7 and r>=3,

    C_bw-C_aw+C_wa-C_wb=8*3^(r-2)+18*2^(r-2)-6+4j>0.

Either whole converse or internal chain reversal negates this contrast.
Thus these exponential intervals give exact growing separators, not a
full-response collision. The prior packet retains complete named cubic
responses for all three states at orders7 and12. Its five-run hostile also
prevents transferring the scalar three-run cap to arbitrary port probes.

At each r there are only three attained values in this family; the formula
does not supply an inverse for an arbitrary odd target or an indefinitely
reusable operation on all source states. That missing coverage assertion
cannot be inferred from the simultaneous exponential terms.

## Implication for the missing spectrum question

The odd-cycle conflict-graph identity **THM-002**,
[ocf](../../../01-canon/theorems/THM-002-ocf.md), is another exact bridge:

    H(T)=I(Omega(T),2).

Its target is a restricted image of conflict graphs of actual tournaments.
Arbitrary independence polynomials are not enough: complete graphs K_3
and K_10 have I(2)=7 and21, precisely values excluded in the tournament
image. Realizability and cycle-overlap support cannot be omitted.

The integral root-frame compiler above solves feasibility of a proposed
form. It does not prescribe its Hamilton count. A block direct sum of two
nonempty star-frame forms generally introduces forbidden zero cross-entries;
it is not a free additive tournament construction. The known ordered-join
operation gives multiplication, which still needs new strong prime atoms.

The current exact reduction is **THM-4094**,
[hamiltonian-matching-deficit-and-two-prime-lane-completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md):
with its proved finite carry atoms, completeness is equivalent to strong
realizations of every odd prime p!=7 and every7p with odd prime p!=3.
The incoming three-state exponential intervals and finite edit repairs do
not supply both infinite lanes.

A route through the lattice would need an H-controlled operation on
feasible integral forms producing those targets, with a proof that the
2-SAT constraints stay satisfiable. A route through the walk model needs
a visited-set state, not only its spectral transfer matrix. These are
precise remaining obligations; neither was discharged in this session.

