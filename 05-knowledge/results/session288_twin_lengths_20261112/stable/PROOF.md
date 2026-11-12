# Stable lengths as a simplicial complex and a finite union of exact cylinders

**PROVED structural statements; FINITE-EXACT source-bank compilation.**
The global stable-port conjecture and Chan–Pak remain OPEN.

## 1. Provenance and the precise question

Use the selected-set interface in `five_stress_connection_atlas_20261112/SEEDS.md`:
`delta=H-W(O)`, where W counts ordered two-nonempty-path spanning covers with
first end and second start selected. U counts the same covers without the
first-end restriction, plus Hamilton paths starting selected. A unit seed has
delta=1. Its inherited physical tail has count `1+2^t U`, t>=1.

The finite atlas proves joint stability for its whole exact stable-port bank
through order8, and its stable-twin contraction lemma is inherited. The
canonical maximal-transitive-block partition is already proved in
`session197/converse-clone-lengths/PROOF.md`, sections1–2. The path-run
substitution mechanism is THM-3121, `path-cover-walk-content-substitution-kernel`.
We do not rename these mechanisms as new discoveries. What follows organizes
the *admissible free length coordinates*, removes overlaps between inherited
families, and derives exact all-size descendant and physical-count bounds.

## 2. Stable faces and their obstruction ideal

Fix a unit marked tournament `(T,O)` and a set P of unselected ports. For
`S subset P`, call S a stable face if independently replacing its vertices by
nonempty transitive modules keeps delta=1 at every positive length tuple.
The empty set is a face, and setting an arbitrary coordinate to1 proves that
every subset of a face is a face. Thus the stable subsets form a simplicial
complex K. This definition does not assume P transitive.

For an **acyclic** P the following exact finite compiler applies. Put
`k=|T|-|P|` and `D=k+2`. In an ordered cover with at most two paths, compress
maximal port runs. Between successive fixed exterior vertices, the walk in
P moves strictly forward in its acyclic order. Consequently each port has
at most k+2 runs, independent of all module lengths. The standard lift factor
for r runs in a chain of length m is `F_r(m)=r! S(m,r)`. Therefore

    delta(m)=sum_(1<=r_i<=D) a_r product_i F_(r_i)(m_i).       (1)

Zero-run terms vanish on positive lengths. `tensor.py` compiles all compressed
skeletons directly: each fixed vertex is visited once, port visits are counted,
and an optional single cut is allowed only between two selected vertices.
Uncut words contribute +1; cut words contribute -1. There are no port self
steps, because such steps belong inside a maximal run. Finite termination is
justified by the preceding acyclic run bound, not by a sampled maximum.

Convert (1) to positive exponential bases1,...,D, then use the anchored basis

    phi_b(m)=b^m-b (b=2,...,D),   phi_0(m)=1.

The functions are linearly independent on positive integers by the
Vandermonde determinant (equivalently the triangular Stirling evaluation).
Write uniquely

    delta(m)-1=sum_alpha c_alpha product_(i:alpha_i>0) phi_(alpha_i)(m_i).

The constant is zero because the original state is a unit seed. Setting all
coordinates outside S to1 kills exactly the terms whose supports are not
contained in S. Linear independence proves the exact iff:

> S is stable iff no nonzero anchored coefficient has support contained in S.

Thus inclusion-minimal nonzero supports are the minimal nonfaces of K. Their
squarefree monomials generate its Stanley–Reisner ideal. This is an exact
map from a physical signed path-cover tensor to a combinatorial obstruction
ideal. It discards coefficients and hence does not classify isolated resonant
length tuples. The full tensor is the required sidecar for those tuples.

**Boundary:** for cyclic P the number of runs can be unbounded, so (1) is not
a general finite tensor theorem. K is still defined, but this compiler has no
claim to decide all its faces. Also absence of singleton obstructions does not
algebraically imply absence of mixed obstructions: the formal function
`1+(2^a-2)(3^b-3)` has two stable singleton faces and a forbidden pair. This
formal hostile is **not** asserted to be realizable by a tournament. Proving
that tournament geometry excludes such behavior is exactly the remaining
stable-port problem, not a consequence of tensor algebra.

## 3. Canonical marked twin blocks and cylinder counting

Partition the unselected vertices into maximal unselected transitive modules;
keep every selected vertex as a singleton. This is the marked version of the
inherited canonical transitive partition. Intersecting transitive modules have
a transitive modular union, which proves uniqueness. Equivalently, connect
unselected twin pairs; its components are the same chains. Here twins means
identical orientations to every vertex outside the pair. Each vertex has at
most one successor twin and one predecessor twin, and these chains cannot
cycle. This gives a direct computation of the canonical blocks.

Positive inflation changes only the corresponding block length and never
merges distinct canonical marked blocks. Indeed, contracting an alleged new
larger transitive module would give a forbidden larger old module; selected
singletons remain separating marks. Consequently marked isomorphism is exactly
weighted-quotient isomorphism. The weights and selected status are essential.

Take a jointly stable set of unselected ports. Multiple stable ports in one
canonical chain represent **one** clock, its total length, not independent
new structure. Contracting the stable part of such a chain to any positive
length preserves the seed, not only lengths at least the original size:
the finite exponential response is constant on the infinite product tail,
and repeated univariate exponential independence extends that identity to
all positive chain lengths. This argument applies simultaneously because the
original stable face is acyclic. Frozen chain lengths must not be discarded.

Each source therefore gives a cylinder in the positive-integer weight lattice:
some block lengths are arbitrary positive integers and all other lengths are
fixed. Canonical quotient automorphisms are retained. Intersections of two
cylinders are either empty (conflicting fixed lengths) or another cylinder.
Finite inclusion–exclusion computes the exact union at every order. A cylinder
with d free coordinates and total minimum order b has series

    z^b/(1-z)^d.                                         (2)

For a quotient automorphism group G, include every orbit translate of each
cylinder before taking the union, and divide the weighted-lattice count by
|G|. This division is valid here because every resulting marked unit seed
has trivial marked automorphism group: the group acts freely on Hamilton
paths and the selected two-path covers, so its order divides delta=1. In the
particular compiled bank below all 88 quotient groups are already trivial.

More generally, if a chain-reduced rigid marked chart has acyclic port set P
and stable complex K, its freely inflatable face union has series
`z^|T| sum_(S in K) (z/(1-z))^|S|`, where S is the set of coordinates whose
length is strictly greater than1. This is a direct physical realization of
the face polynomial. It does not count isolated nonfree resonances.

## 4. Exact all-size descendant classification of the inherited bank

**FINITE-EXACT INPUT:** the atlas's complete strong marked unit-seed census
through order8, including its exact stable-port coefficient certificates.
We do not independently rerun that enormous census. The source hash is
recorded in `cones.json` and `normalforms.json`.

`cones.py` canonically compresses every one of the 161 states with at least
one stable port, expands the finite quotient symmetry action, removes contained
cylinders, and computes their exact intersections. It then **recomputes the
full compressed-run defect tensor for every retained cylinder** on its actual
reduced graph, rather than assuming that contraction preserved a multiport
identity. All 98 tensors are identically1. There are88 quotient types,
79 one-clock cylinders and19 two-clock cylinders. The union series is

    C(z)=(z^3+4z^5+8z^6+20z^7+46z^8)/(1-z)
         +(5z^7+14z^8)/(1-z)^2.                         (3)

The first counts, n=1,...,10, are

    0,0,1,1,5,13,38,103,122,141.

In particular the union contains **exactly19n-49 marked isomorphism classes
at every n>=8**. This is an all-size result for the specified bank's
independent stable-clone descendants, not a prediction of the full sequence
of all unit seeds or all stable-clone births. It says nothing about descendants
under additional graph-changing operations.

A literal two-clock witness has rows `(14,28,8,16,37,79,31)`, selected mask3,
and ports2,6. Its full tensor is only `a_(1,1)=1`. Its two independent positive
lengths give exactly n-6 distinct marked seeds at every order n>=7. Neither
port can be contracted with an unselected twin in its seven-vertex base.

## 5. Stable reduction: exact finite receipt, conditional global confluence

The inherited stable-twin contraction reduces order while preserving a unit
seed and stability of the contracted port. If two reductions overlap in one
chain, they have the same contracted marked graph up to isomorphism. If they
act in different chains whose two directions are jointly stable, their
contractions commute and each direction remains stable afterward, by the
positive-length identity in section3.

Therefore on any reduction-closed class in which every pair of individually
stable directions is jointly stable, stable-twin contraction is confluent:
termination plus these local diamonds gives a unique terminal marked class.
This is **CONDITIONAL** globally. The atlas's global stable-port hypothesis
would supply the condition, but this session does not prove that hypothesis.

The exact inherited bank through8 has161 states,63 distinct reduction arcs,
and a unique terminal class for every state. There are **no nontrivial local
diamonds** in that finite quotient DAG, so the finite calculation does not
stress-test genuinely competing stable directions. The full terminal-order
histograms appear in `normalforms.json`; a class born at7 has two distinct
order8 stable-clone descendants, so counting one descendant per birth would
already be wrong.

## 6. Why this finite bank cannot cover the Chan–Pak tail

This is an application of the inherited fixed-template polylogarithmic
obstruction in `session175_tournament_spectrum_20261104/clone_obstruction/PROOF.md`,
also explicitly credited in Session181. The new specialization retains U and
the additional tail clock, and inserts the exact79/19 cylinder counts.

**PROVED for every finite strong stable-clone bank.** Let d be the number of
free twin-block lengths in one cylinder. Every vertex of a strong tournament
lies on a directed triangle: otherwise all its in-neighbors beat all its
out-neighbors and there is no path back from the latter to the former. Hence
an unselected transitive block of length m lies in an induced
`C3[TT_m,1,1]`. The odd-cycle formula gives H>=1+2^m: choose an odd nonempty
subset of the chain to obtain2^(m-1) directed odd cycles and their weight2.
The extra vertices cannot remove those odd-cycle contributions.

For a unit seed, U>=W=H-1 since every cover taxed by W is included in U.
Consequently each free block length satisfies U>=2^m. If a literal tail
output `1+2^t U` is at most X, then each m and t is at most
`L=floor(log2(X-1))`. One d-clock cylinder supplies at most L^(d+1) such
outputs. Collisions can only reduce this upper bound. Thus the98-cylinder
bank supplies at most

    79 L^2 +19 L^3                                      (4)

counts at most X. In particular its physical tail outputs have density zero.
This concerns **actual counts bounded in height**, not odd parts of slopes
of arbitrarily high seeds. It is stronger than a bounded sample missing values,
but does not rule out infinitely many new cores, new operations between
cylinders, or Chan–Pak itself.

## 7. Targeted hostile exploration and remaining frontier

`explore.py` chooses12,000 distinct deterministic one-vertex extensions of
inherited order6–8 stable seeds (seed288), with the old selection retained or
the new vertex additionally selected. There are268 unit seeds;245 have at
least two ports passing single-length2 and3 screens. No pair of such ports
fails the mixed(2,2) test. These are merely screens, **not** all-length
stability certificates and not a complete order9 census. Normal/optimized
receipts are supplied. The two known false length2 positives are retained as
exact tensor controls elsewhere in this lane.

This does not establish closure under adaptive moves that become stable only
after a special inflation. Such newly unlocked directions need their own
coefficient test; (3) classifies the original cylinder union only.

The constructive task is now sharper: free twin lengths form controlled
cylinders, but their finite union has too little height-bounded supply.
A useful new move must change the frozen length constraints/core, realize
nonfree resonances systematically, or compose different charts with a
retained endpoint interface. The obstruction ideal separates the free face
problem from that resonance problem. Neither is solved for every tournament.
