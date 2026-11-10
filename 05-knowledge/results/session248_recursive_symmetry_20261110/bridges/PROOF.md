# Recursive symmetry as a vertex-support spectrum

**Status: PROVED elementary new operator formulas and deductions; inherited
Walsh orbit-basis theorem is CITED as source-reported, not re-audited.**
The finite check concerns the new transition operator only. No external
literature priority claim, no new canon identifier, and no claim to solve
the existing tiling-line metagraph's spectrum.

## 0. Inheritance, discovery boundary, and live concepts

Closest mechanism: THM-4721, section 5, supplies signed graph orbit sums as
a basis of tournament isomorphism invariants:
`01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`.
It is the source-reported theorem W, not a new result here. Its underlying
monomial action is THM-4531, theorem G:
`01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md`.

A much earlier complementary source is
`01-canon/theorems/THM-071-walsh-fourier-diagonalization.md`:
its order-5/order-7 Hamiltonian-count decomposition organizes **edge degree**.
Its own general-n conjectural clauses are not imported. The new operator
instead organizes **vertex support**, giving an all-order invariant spectrum
without needing any general Hamiltonian-count decomposition.

Canonical hostile and least-used sidecar:
`05-knowledge/results/overnight_hexagon_sep05_boolean_consumer.md`, sections
1–2, separates multiplicity walks, simple support graphs, and stationary
weights. We retain transition multiplicities and orbit volume throughout.
That source studies Eulerian graphs (all degrees even), whereas the Fourier
partners here are **sign-even graphs** (all automorphisms have positive
orientation sign). These meanings of “even” must not be conflated.

Corrected near miss: THM-4531's failure of a natural set bijection does not
prevent an equivariant **linear** Fourier transport. Our map is a map of
function spaces and Markov operators, not a map pairing tournament types
with graph types and preserving their automorphism groups.

Concept board:
1. Recursive birth of invariants by vertex support.
2. Signed-even graph Fourier basis.
3. Random vertex forgetting and fair extension.
4. Edge-degree noise versus vertex-support noise.
5. Orbit-volume weighting versus Boolean adjacency.
6. Individual symmetry groups versus statistical relabelling invariance.

The new content is the exact all-order resampling spectrum, its recursive
birth multiplicities, its deletion/extension interpretation, and its
commuting two-parameter noise calculus. The basis and equinumerosity are
explicitly inherited. These mechanisms are sufficiently elementary that
no claim of global mathematical novelty is intended.

## 1. Objects and the signed basis

Let T_n be labelled tournaments on [n], encoded by independent coordinates
x_ij in {−1,+1}, i<j, with +1 meaning i→j. Write μ_n for the uniform law,
and I_n for the real vector space of functions constant on tournament
isomorphism classes. Its dimension t_n is the number of tournament types;
t_0=t_1=t_2=1. For an ordinary simple graph F on [n],

    χ_F(T) = product over ij in E(F) of x_ij(T).

Relabelling permutes these Walsh monomials with signs, because an arc whose
endpoints change order changes its coordinate sign. A graph F is sign-even
if every graph automorphism reverses an even number of its ordered edges.
For such F the appropriately signed sum Φ_F over distinct labelled copies
is nonzero and invariant. A choice of orientation convention fixes its
overall sign. For a non-sign-even graph the invariant projection is zero.

CITED theorem W says that one Φ_F for each sign-even graph type on [n]
forms an orthogonal basis of I_n. The elementary reason, recalled to fix
scope, is that the Walsh basis decomposes into signed permutation orbits;
an orbit contributes a one-dimensional invariant line exactly when its
stabilizer acts trivially on that line. This is linear invariant theory,
not an identification of the vertex automorphism groups of F and T.

Define v(F) as the number of nonisolated vertices and e(F)=|E(F)|. The
empty graph has v=e=0. Neither parameter changes within an orbit.

## 2. The retention-resampling metagraph

For 0≤k≤n, define K_(n,k) by the following concrete random move:

1. Choose a uniformly random k-subset S of [n].
2. Retain every arc with both endpoints in S.
3. Independently choose a fair orientation for every other pair.

The output again has n vertices; the retained induced tournament is the
cross-size interface. Relabelling commutes with this experiment, so it
induces a weighted directed kernel on tournament types. It is reversible,
not generally uniform over types. Explicitly its stationary probability is

    π_n([T]) = (n! / |Aut(T)|) / 2^binom(n,2).

This is a metagraph in the precise sense of a weighted transition graph
on isomorphism classes. It is NOT asserted to be the repository's tiling
line metagraph, and is not an arbitrary orientation of a tie-filled graph.

### Theorem B1 — exact support spectrum (PROVED)

For every labelled graph F,

    K_(n,k) χ_F = λ_v(n,k) χ_F,
    λ_v(n,k) = binom(k,v)/binom(n,v),  v=v(F),

where λ_0=1 and binom(k,v)=0 for v>k. Therefore the same formula holds for
Φ_F. The full spectrum on I_n is obtained by taking one copy of λ_v for
each sign-even graph type of nonisolated support v.

**Proof.** Condition on S. If an edge of F is not entirely retained, its
new independent fair sign occurs once in the monomial, making its
conditional expectation zero. Otherwise χ_F is unchanged. The latter
event says precisely V_nonisolated(F)⊆S and has probability
binom(n−v,k−v)/binom(n,k)=binom(k,v)/binom(n,v). Signed orbit sums preserve
v and so inherit the eigenvalue. This proves the statement without
assuming any formula for tournament enumeration. Each fixed-S operation
is conditional expectation in L²(μ_n), hence self-adjoint; its average is
self-adjoint and equivariant. Lumping by isomorphism gives the asserted
orbit-volume reversibility. ∎

### Theorem B2 — recursive birth dimensions (PROVED)

Let c_v be the number of sign-even graph types on exactly v vertices with
no isolated vertices; set c_0=1. Then

    t_n = sum_(v=0)^n c_v,
    c_v = t_v − t_(v−1) for v≥1.

Consequently, for every 0≤k≤n,

    rank(K_(n,k) on I_n) = t_k,
    dim ker(K_(n,k) on I_n) = t_n − t_k.

Its eigenvalue list with multiplicities is

    1 once;  λ_v(n,k) with multiplicity t_v−t_(v−1), 1≤v≤n.

Coincident eigenvalues must of course be combined; all v>k contribute to
one zero eigenspace.

**Proof.** Adding isolated vertices preserves sign-evenness. For a graph
with at least one edge its nonisolated vertices are characteristic; its
automorphism group is the product of that of the nonisolated graph and
the symmetric group on isolates. Permuting isolates reverses no graph
edges. Removing all isolates therefore gives a unique sign-even
isolate-free type, proving the sum by theorem W. The empty graph is handled
separately. The rank and kernel follow because λ_v>0 exactly for v≤k. ∎

Interpretation: a genuinely new invariant direction appears at size v
exactly when a sign-even Fourier partner first uses all v vertices. This
is a precise recursive emergence of **invariant information**, not a
claim that every new invariant direction creates a new nontrivial
vertex-automorphism group.

## 3. Exact gap and its symmetry acceleration

### Corollary B3 (PROVED)

For n≥3 and 0≤k<n, the nontrivial spectral gap of the type kernel is

    gap(K_(n,k)|I_n) = 1 − binom(k,3)/binom(n,3).

In particular the one-vertex-refresh kernel K_(n,n−1) has gap 3/n.

**Proof.** No nonempty sign-even graph has support 1 or 2: support 1 is
impossible, while the single edge has a sign-reversing endpoint swap.
The three-vertex path P3 is sign-even: its leaf swap reverses both edges.
It spans a nonconstant invariant direction, so c_3=1. For fixed k<n the
positive values λ_v strictly decrease with v, since
λ_(v+1)/λ_v=(k−v)/(n−v)<1. Thus λ_3 is the largest nonconstant eigenvalue;
if k<3 every such eigenvalue is zero. Finally
binom(n−1,3)/binom(n,3)=(n−3)/n. ∎

For comparison, on **all labelled functions**, the first nonconstant
Walsh monomials are single edges with v=2. The labelled one-vertex-refresh
gap is 2/n. Passing to isomorphism-invariant observations removes this
entire slow sector and gives a factor 3/2 increase in the exact gap.
This is a clean mathematical effect of symmetry reduction on relaxation,
not a claim about total-variation cutoff times.

At k=2 the type kernel is rank one even though a labelled edge survives.
An unoriented two-vertex isomorphism type carries no information, so this
apparently retained information disappears entirely after quotienting.
This is the smallest hostile test for confusing labelled retention with
retention of invariant information.

## 4. A genuine size-changing map and its loss ledger

Define D_(n,k):I_k→I_n by averaging a k-vertex invariant over all induced
k-subsets. Define R_(k,n):I_n→I_k by extending a given k-tournament with
n−k labelled vertices and fair independent remaining arcs, and taking
expected value. Relabelling invariance makes these definitions independent
of the chosen identification of a subset with [k]. Then

    K_(n,k)|I_n = D_(n,k) R_(k,n).

For an isolate-free sign-even graph F on v≤k vertices, define a coherent
signed embedding sum

    A_(F,n)(T) = sum over injections f:V(F)→[n]
                    product over oriented edges a→b of F of x_(f(a),f(b)),

using antisymmetric coordinates x_ba=−x_ab. Fix the reference orientation
of F once. The sum is an invariant polynomial. Sign-evenness ensures
that its equal labelled edge-support terms do not cancel, so it is a
nonzero constant multiple of Φ_(F plus isolates). Write

    U_(F,n) = A_(F,n)/(n)_v,

where (n)_v is the falling factorial. Direct counting gives

    D_(n,k) U_(F,k) = U_(F,n),
    R_(k,n) U_(F,n) = ((k)_v/(n)_v) U_(F,k).

The first equation counts the chance that a particular injection is
contained in the random subset; the second kills every monomial touching
a new vertex. These are compatible normalized coordinates as size grows,
not merely a dimension count. Unnormalized sums obey
R_(k,n) A_(F,n)=A_(F,k).

The complete image of D_(n,k) is exactly the span of the support-v≤k
modes. Thus small induced-subtournament averages retain all and only those
modes. They lose every support-v>k coefficient, even when the original
statistic has low edge degree compared with other statistics.

Map/loss ledger:
- Source: real functions of tournament isomorphism type at order n.
- Target: expected statistics of a retained order-k tournament.
- Map: R_(k,n), with D_(n,k) the return map.
- Preserved: all support≤k modes up to explicit nonzero scalar factors.
- Destroyed: all support>k modes exactly, not approximately.
- Restoration sidecar: the omitted signed-even-graph Fourier coefficients.
- Cheapest decisive test: at n=4,k=3, eigenvalues 1,1/4,0,0; two new
  size-four invariant directions disappear while the P3 mode survives.
- Scope boundary: retaining these coefficients does not pair the graph
  automorphism group with the tournament automorphism group.

The maps R are expectation operators, not deterministic recovery of an
individual tournament from one sampled card. Likewise equality of all
small-support expectations need not imply isomorphism.

### Cross-lane interpretation: a down-up chain

The [birth analysis](../birth/PROOF.md) supplies the type-level
birth probabilities B_j(T,U)=a(T,U)/2^j and uniform-vertex-deletion
probabilities D_(j+1)(U,T)=d(U,T)/(j+1), together with their orbit-volume
balance law. Successively delete uniformly to size k, then add fair new
vertices to size n. The retained k-set is uniform; all subsequently
exposed edges are independent fair orientations. This is exactly K_(n,k)
after forgetting labels. Thus the spectral decomposition here diagonalizes
the down-up composition built from that lane's birth/death kernels.

In particular, at k=n−1 the birth lane's finite combinatorial transition
counts have spectrum 1 and 1−v/n with multiplicities t_v−t_(v−1). The
balance law explains the correct inner product; the signed-even graph
basis explains every relaxation scale. Neither ingredient alone supplies
both pieces.

### Corollary B4 — the graded birth-death commutation relation (PROVED)

Use row-stochastic chronological multiplication in this paragraph, so
B_n D_(n+1) means birth from n to n+1 followed by deletion back to n.
For n≥1,

    (n+1) B_n D_(n+1) − n D_n B_(n−1) = I.

Indeed the deletion selects the newborn with probability 1/(n+1), leaving
the original type unchanged. With probability n/(n+1) it selects a
uniform old vertex; the newborn is then a fair replacement of that vertex,
which is exactly D_n B_(n−1)=K_(n,n−1). This proves the identity of kernels,
including transition multiplicities, not merely an identity of spectra.
For n=1 the size-zero birth/death kernels have the evident unique-state
interpretation.

Therefore the birth-then-delete operator on order n has eigenvalues
1−v/(n+1), with multiplicity c_v for 0≤v≤n. It has no zero eigenvalues,
whereas delete-then-birth at order n+1 has exactly c_(n+1) zero modes.
The nonzero spectra coincide across the adjacent sizes, and all newly
born invariant directions account for the discrepancy. This is a precise
graded operator relation tying recursive birth to support relaxation.
It is not asserted to be a new abstract theory of differential posets.

## 5. Two-axis noise and a commuting calculus

Let N_ρ independently reverse each arc with probability (1−ρ)/2,
0≤ρ≤1. The standard Walsh calculation gives

    N_ρ χ_F = ρ^e(F) χ_F.

Combining with B1 gives a new operational bigrading of the invariant space:

    N_ρ K_(n,k) Φ_F = ρ^e(F) [binom(k,v(F))/binom(n,v(F))] Φ_F.

All retention kernels K_(n,k), and every N_ρ, commute. Hence edge-degree
and vertex-support relaxation are simultaneously measurable; edge noise
alone conflates graphs with equal e but different v, and vertex refresh
alone conflates equal-v modes with different e. The common eigenspaces of the entire two-parameter family
are indexed by (e,v), with dimension equal to the number of sign-even
graph types with those parameters and no isolates.

For an explicit control, P3 has (e,v)=(2,3), P5 has (4,5), and the
disjoint union of two P3's has (4,6); each is sign-even. For P5 the
reflection reverses four edges. For two P3's, component automorphisms
have positive sign and exchanging the components contributes the same
orientation-parity twice. Thus edge noise cannot distinguish the last
two modes, but vertex-retention noise does for n≥6 and 5≤k<n. This uses
no C3-module or wreath-product tournament construction.

A further exact consequence is the support-energy decay formula. If
f−E f = sum_(F nonempty) a_F ψ_F in an orthonormal signed-orbit basis,
then after r independent retention steps,

    ||K_(n,k)^r(f−E f)||² = sum_F a_F² λ_v(F)^(2r).

For one-vertex refresh λ_v=1−v/n. This is a complete hierarchy of
relaxation rates indexed by the number of vertices an invariant needs,
not by the number of arcs its Fourier monomial uses.

## 6. Targeted finite check of the new operator

`check_resampling.py` independently constructs tournament types by direct
permutation canonicalization for n=3,4,5. For every k<n it enumerates the
actual retained subsets and all fair arc completions, constructs the
rational type transition matrix, and checks its characteristic polynomial,
row sums, and detailed balance with labelled-orbit volumes. It does not
construct the Fourier basis or audit the old even-graph theorem.

Reproduce from repository root:

    python3 05-knowledge/results/session248_recursive_symmetry_20261110/bridges/check_resampling.py

The recorded output `check_resampling.out` is FINITE-EXACT for these
universes. Examples:

- n=3,k=2: spectrum 1,0.
- n=4,k=3: spectrum 1,1/4,0,0.
- n=5,k=4: spectrum 1,2/5,1/5,1/5, and eight zeroes.

Controls include k=0 (complete reset), k=2 (retained labelled edge but
complete type reset), n−1 (predicted gap 3/n), and exact reversibility.
The finite checks support the new formulas; the all-order proof is B1–B3.

## 7. What this opens, and what remains open

PROVED: the recursive invariant dimension jump t_n−t_(n−1) is precisely
the nullity of one-vertex-refresh on type functions. Every old direction
survives with an explicitly known rate, and every genuinely n-supported
direction dies in one step. This is an exact separation of inherited and
newly born invariant information as size grows.

OPEN / NEXT: compare the actual tiling-line transition operator with this
support grading. Its commutator with K_(n,n−1) would measure how strongly
line adjacency mixes different birth scales. A zero commutator would mean
that the line dynamics respects recursive support; a nonzero block would
identify a precise coupling between support sizes. No such commutation
or transfer is claimed here: the labelled sampling mechanism and its
multiplicities must first be specified for that distinct metagraph.

OPEN / NEXT: classify which signed-even support modes detect nontrivial
tournament automorphisms efficiently. The full basis determines the type,
but the spectral theorem alone gives no small-support certificate for
rigidity or a prescribed symmetry group. This distinction keeps the
statistical-symmetry bridge from silently becoming a false group theorem.
