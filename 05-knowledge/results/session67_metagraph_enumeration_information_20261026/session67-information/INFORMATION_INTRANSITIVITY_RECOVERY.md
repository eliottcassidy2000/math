# Tournaments as information: recovered lineage and a lossless family enumeration

**Session67 — 2026-10-26.** Source checkout read-only, pinned at
`b4234f728c0516b5fc521c001d69e26f3d2508d5`.

**Status.** Historical claims are distinguished below from their repairs.
Elementary identities, the weighted-core classification, and the realization-orbit
formula are **PROVED** here. Recovery of the earlier C3/interval-reversal theorem
is **CITED**, not a new result. The supplied independent programs are
**FINITE-EXACT** through order five. This note does not claim a new globally novel
classification, efficient all-order enumeration, or a solution of LRC.

## 1. What was actually recovered

The requested relation-first and information-first work predates the current
session by months. `SOURCE_MANIFEST.json` records exact repository paths, fixed
commit, blob object, SHA-256, local immutable copies, and last-change commits for
every item below. Sparse checkout absence did not mean absence from the repository:
these files were recovered with read-only `git show PIN:PATH`.

| Source | Historical date/evidence | Relevant recovered content | Status for this session |
|---|---|---|---|
| S01–S02: `04-computation/information_theory_deep.py`, `05-knowledge/results/information_theory_deep.out` | Header self-identifies `opus-2026-03-13-S67j`; latest file commit is `2476ce4fbe27a5d5ff644c860ba11a48734f881d`, June 12 | Pairwise comparisons as a noisy channel; path entropy; information lost by scores | Provenance recovered. Channel interpretations need the repairs in §3; header date is not independently authenticated creation time. |
| S03–S04: `04-computation/ocf_information_theory.py`, `05-knowledge/results/ocf_information_theory.out` | Header self-identifies March 13, S67k; same June 12 archival commit | Odd-cycle coefficients as information channels; loss from collapsing them to one Hamiltonian count | Useful quotient viewpoint, not an injective binary code: coefficients can carry into other powers of two. |
| S05: `05-knowledge/hypotheses/HYP-3564-relations-not-things-the-relational-reframe.md` | June 29; `05df17e60a66c95187cbf06538c41bff8ef57237` | Relation, rather than its endpoints, is the object; intransitivity is relational content | REFRAME, not proof. “Only quadratic invariant” does not mean “only invariant of any degree.” No LRC conclusion is inherited. |
| S06: `05-knowledge/hypotheses/HYP-3595-rib-sc-side-even-and-silent-flip-eulerian.md` | June 29; `b9d29da1ba69a9e1a4cdd9a1b61674146c035ed6` | Already warns that the cut/cycle cardinality split is labelled-only and stabilizers do not decouple | This warning takes precedence over the next day's simpler slogan. Its GF(2) cut coordinates must not be called integer scores. |
| S07: `05-knowledge/hypotheses/HYP-3602-tournaments-are-intransitivity-among-n-things.md` | June 30; `ea7a9474cd4e4a791b448bccc0f6131728de54bc` | Explicit relation-first ontology, triangles as irreducible obstructions to ranking, “order × shape” proposal | Correct core intuition; literal coboundary/cohomology and unlabelled-product statements need repair. |
| S08: `01-canon/theorems/THM-1382-map-graph-star-flips-nested-cut-cycle.md` | July 20; `107cb4957769c83699d0b308aec13530fd99525e` | Star-orbit invariants belong to a path-marked tiling; they do not descend to tournament classes | PROVED negative boundary; relevant warning against counting gauge-dependent shadows as classes. |
| S09: `01-canon/theorems/THM-1979-tournament-space-is-a-spectrum-single-point-to-continuum.md` | Corrected July 21; `98c31416f4956765b7ec4cb9dd8d4e4406ce3185` | Score variance determines triangle count but not score shape or the structural fibre | CORRECTED source, used rather than the earlier monotonic-richness story. |
| S10: `01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md` | October 21 source label; committed `894dfee3681f0992f1c1c7f5372793caf4224419` | Exact score/curl split and metagraph line currents; compressed port data do not determine the graph | Its elementary algebra is independently proved below. The reported 48,000 order-six L3-compatible transports remain inherited finite evidence, not rerun here. |
| S11–S13: Session36 files `context/CURL-AND-CLAN-COMPLETION.md`, `context/REALIZATION_PHASES.md`, `lattice/METRIC_CURL_RECOVERY.md` under `05-knowledge/results/session36_information_curl_contexts_20261023/` | October 23; `5390aea9cbe0a23328eef89eaf95f1d85f299d01`, `564fbdc91c129583cf2f25e29fec06044e48bc06` | Exact unsigned-curl response quotient; XOR/NAE realization; signed-curl fibre and missing cut | Existing in-repository results, recovered rather than claimed anew. |

Connected saved-conversation searches for “intransitivity”, “binary
intransitivity”, and “tournaments information” returned no older saved ChatGPT
matches; the room search returned recent October messages. That is a search
scope result, not proof no older conversation exists. The exact repository
objects above supply the historical evidence.

## 2. Binary intransitivity is intrinsic; linear curl is a different object

Let `A_ij=1` mean `i -> j`, with `A_ji=1-A_ij` for distinct vertices, and
`S_ij=A_ij-A_ji` in `{+1,-1}`. There are no ties. Relabelling is the only
identification in A000568; converse is NOT silently identified.

For three distinct vertices let

```
c_ijk = S_ij + S_jk + S_ki,
chi_ijk = A_ij A_jk A_ki + A_ik A_kj A_ji.
```

Then `chi` is the binary cyclic-triple indicator and

```
chi_ijk = (c_ijk^2 - 1)/8.
```

**Proof.** Around a cyclic triangle all three cyclic-direction signs agree, so
`c=+3` or `-3`. Otherwise precisely one or two agree and `c=+1` or `-1`.
The displayed formula distinguishes those two cases exactly. QED.

Consequently a tournament is transitive iff every `chi` vanishes. To see the
nontrivial direction, choose a shortest directed cycle; if its length exceeds
three, its chord between the first and third vertices either shortens it or
makes a directed triangle. Thus no triangle implies acyclicity and a total order.
Equivalently there exists an injective real ranking `p` with
`S_ij=sign(p_i-p_j)`.

**Important correction.** This does NOT assert `S_ij=p_i-p_j`. Already the
transitive order `0 -> 1 -> 2`, `0 -> 2` has
`c_012=1+1-1=1`, not zero. The real ±1 flow of a total order is not a gradient.
On the full clique complex of `K_n`, first cohomology is zero; on the graph
alone the cycle space exists, but it is not an invariant-dependent “cohomology
class of binary intransitivity.” Those are different coefficient systems and
complexes. This is the precise repair to HYP-3602, not a rejection of its
relation-first viewpoint. [HodgeRank primary source](https://arxiv.org/abs/0811.1067)
likewise treats real edge-flow decomposition, not equality of a sign ranking
with an exact additive potential.

Let `r_i=sum_j S_ij=2s_i-(n-1)`. Summing curl over the third vertex gives the
lossless labelled formula

```
n S_ij = r_i-r_j + sum_k c_ijk.                    (2.1)
```

Repeated-index curls are zero. Indeed the right-hand sum of curls is
`n S_ij+r_j-r_i`. Thus **scores plus the full signed, addressed triangle-curl
field recover every original binary relation**. Scores alone, unsigned curls,
multisets of curls, and curl energy need not.

The least-squares projection of `S` onto gradients has potential `p_i=r_i/n`.
With one term per unordered edge,

```
E_grad = (n^2-1)/3 - 8 c3/n,
E_residual = (n-1)(n-2)/6 + 8 c3/n.                (2.2)
```

For `n>=1` these follow from
`sum_i binom(s_i,2)=binom(n,3)-c3` and orthogonal projection on the complete
graph. The positive baseline at `c3=0` is indispensable. Scalar energy records
only `c3`; it is not a full classification.

A sharp GF(2) warning: in sorted edge order `01,02,12`, the transitive bit word
`111` becomes the cyclic word `010` after switching the cut at vertex 1. Both
have odd edge parity. Thus a switching-invariant cycle shadow can merge the
transitive and intransitive relations. The equality of labelled cardinalities

```
2^binom(n,2) = 2^(n-1) 2^binom(n-1,2)
```

is not an intrinsic unlabelled product. Over GF(2), cut and cycle spaces of
`K_n` are not even complementary when `n` is even: their intersection has
dimension `n-2`. A spanning-tree coordinate split can always be chosen, but it
is a gauge and its action must be carried along.

## 3. Three different entropies: the March information calculation repaired

Write `h(T)` for Hamiltonian-path count, reserving `H_Sh` for Shannon entropy.

1. `log2 h(T)` is entropy of a uniformly selected directed Hamiltonian path
   **conditional on T**. Such a path only satisfies consecutive comparisons;
   it is not a total ranking satisfying all pairwise comparisons.
2. `H_Sh(h(T))` is entropy of the scalar path count under a specified tournament
   distribution. It is a statistic of T, not the entropy in item 1.
3. `H_Sh([T])` is entropy of the unlabelled class under a specified source.

For a uniform labelled tournament on n vertices, `m=binom(n,2)` independent
binary comparisons give `H_Sh(T)=m`. Orbit-stabilizer gives

```
Pr([T]=C) = n! / (2^m |Aut(C)|),
H_Sh([T]) = m-log2(n!)+ E log2 |Aut(T)|.            (3.1)
```

**Proof.** Class C has `n!/|Aut(C)|` labelled realizations. Substitute that
probability into `-sum_C p_C log2 p_C`. QED.

A *uniform class* source instead has entropy `log2 a_n`, `a_n=A000568(n)`.
These are unequal in general: at n=5 the two entropies are respectively
`3.395179984853084...` and `3.584962500721156...`. A worst-case fixed-length
class index uses `ceil(log2 a_n)` bits; this says nothing about how quickly
that index can be constructed.

The March script's channel surrogate `log2(n!)-E log2 h(T)` is not generally
`I(ranking;T)`. At comparison accuracy 1/2 the channel output is independent of
the latent ranking, so the true mutual information is exactly zero. At n=3
the surrogate is `log2 6-(1/4)log2 3 = 2.188721875540867...`, already a strict
counterexample. The related printed claim that Rédei's path theorem makes
every tournament strongly connected is also false: a transitive tournament
is an immediate hostile. None of those historical interpretations is used in
the enumeration.

The OCF expression `h=1+2 alpha_1+4 alpha_2+...` is an exact integer statistic
when its cycle coefficients are correctly defined, but it is not a unique
binary expansion of the coefficient vector. At order six the old output
already records h=9 both from `(alpha_1,alpha_2)=(4,0)` and `(2,1)`.
Coefficients and the full incidence of their cycles are separate information.

## 4. The repaired “order × intransitivity” grammar is a weighted-core orbit

A **module** M is a nonempty vertex set every exterior vertex sees uniformly.
A transitive module induces a total order. A **pair-free** tournament Q has no
two-vertex module. This term must not be confused with a modular-prime
(indecomposable) tournament: larger nontransitive modules may remain.

### Theorem 4.1 — unique maximal-transitive-module coordinates

Every nonempty tournament has a unique partition into maximal transitive
modules. Its quotient Q is pair-free, and the tournament is

```
T = Q[TT_(l_1), ..., TT_(l_q)],     l_i >= 1.        (4.1)
```

Two such presentations are isomorphic exactly when their quotients are
isomorphic by a map carrying the attached positive sizes. For a fixed
representative Q, sizes must therefore be taken modulo `Aut(Q)`.

**Proof.** If two transitive modules overlap, their union is a module: an
exterior vertex has the same relation to both through a shared vertex.
Their union is transitive: any triangle not contained entirely in one module
has two vertices in one of the modules and its third vertex outside, which
cannot be cyclic because of homogeneous exterior adjacency. Thus intersecting
maximal transitive modules coincide. Singletons ensure they cover V.

If two quotient vertices formed a module, the union of the corresponding two
transitive blocks would again be a transitive module, contradicting maximality.
Conversely, given a pair-free Q, a transitive module meeting multiple displayed
blocks could be enlarged by the overlapping-union argument to contain every
block it meets. Its support would be a nontrivial transitive module of Q.
Two consecutive vertices in that module's order would then form a two-vertex
module of Q, a contradiction. Hence the displayed blocks are exactly the
maximal ones. An isomorphism must permute this intrinsic partition, yielding
the weighted quotient isomorphism; conversely such a map lifts using the
unique isomorphisms of the transitive factors. QED.

For each unlabelled pair-free Q, Burnside on positive size vectors supplies
its all-size ordinary generating function

```
F_Q(x) = (1/|Aut Q|) sum_(g in Aut Q)
          product_(cycles C of g) x^|C|/(1-x^|C|),
A(x) = 1 + sum_[pair-free Q, nonempty] F_Q(x).       (4.2)
```

A vector fixed by g has one positive integer weight per cycle, contributing
that cycle's length times its weight. That proves each factor; Theorem 4.1
proves exhaustiveness and disjointness of the total sum. For a coefficient
of degree n only q<=n is relevant, so the formula is well-defined formally.
This does not eliminate the hard task of enumerating all pair-free cores;
it identifies exactly the missing family bank rather than hiding it in a
scalar statistic.

**Information consequence.** Let X be the unlabelled tournament class and
let `K(X)=([Q],[l]_Aut(Q))`. The theorem says K is bijective, hence for EVERY
probability distribution on the finite order-n classes,

```
H_Sh(X | K(X)) = 0,
H_Sh(X) = H_Sh([Q]) + H_Sh([l]_Aut(Q) | [Q]).        (4.3)
```

There is no independence assertion. Core identity and sizes are coupled by
both total order n and automorphism action. This is the precise replacement
for multiplying two unlabelled shape counts. The core still carries the full
binary relation between blocks, not only how many cyclic triples it has.

### Same core and present symmetry do not determine future enumeration weights

The sibling Session67 recursion note (`session67-recursion/REPORT.md`, §3)
identifies the exact child-family transfer

```
K_n(T,G) = sum_(b: T+b in G) |Aut(T+b)| / ((n+1)|Aut(T)|),
```

where b ranges over all binary one-vertex attachments. This is an exact
unrooting weight, not a probability. A transfer matrix based solely on a
parent family must give the same K for every parent in that family if it is
to work for arbitrary input weights on parent classes.

Let `n=3k-1`, k>=2, and take the two strong tournaments

```
T_A = C3[TT_(3k-3),1,1],
T_B = C3[TT_(k-1),TT_k,TT_k].
```

Both have the same pair-free core C3 and trivial automorphism group: a
nonconstant size triple is stabilized by no nonidentity cyclic rotation.
Nevertheless for the target family of core order three,

```
K_n(T_A, core3) = (3k+2)/(3k),
K_n(T_B, core3) = (5k+2)/(3k),                    (4.4)
```

whose difference is always 2/3.

**Proof.** A core-three-preserving extension inserts the new vertex into
one of the existing transitive blocks. A block of size l offers l+1 distinct
insertion positions. Every resulting child of T_A is rigid, so the total
weighted count is `(3k-2)+2+2=3k+2`. For T_B, the k positions in the shorter
block produce the equal-size cyclic blowup with automorphism order three;
the other two blocks contribute k+1 rigid positions each. The weighted total
is `3k+2(k+1)=5k+2`. Divide by n+1=3k. No other core-three child is possible:
deleting the new vertex from its maximal-transitive-module presentation
must leave all three nonempty blocks and recovers the parent's intrinsic
partition. QED.

For k=2 both parents have eight raw core-three extensions, but their
symmetry-weighted sums are eight and twelve. The independent program
`verify_growth_information.py` checks all 32 attachments to each order-five
parent by direct module extraction and permutation automorphism tests,
yielding K=4/3 versus 2.

The missing information is **the size vector's position relative to a
future symmetry**, not merely the current automorphism order. This validates
keeping `[l]_Aut(Q)` in (4.3). It refutes a universal coarse transfer valid
for arbitrary parent weights; it does not rule out analytic global averaging,
aggregate scalar identities, or all useful family recurrences.

## 5. A second exact enumeration directly from binary intransitivity

The complete binary intransitivity field is the labelled 3-uniform hypergraph
`H=C3(T)` of cyclic triples. It still forgets information. The classical
[C3-structure theorem](https://perso.ens-lyon.fr/stephan.thomasse/liste/bilt.pdf)
of Boussaïri, Ille, Lopez and Thomassé (2005 manuscript, Theorem 2) identifies
its labelled realization fibre with transformations by interval reversals;
for an indecomposable realization the only alternatives are it and its
converse. These are CITED results. Session36 proved that the same quotient
is the full ordinary commutative walk series and complete vertexwise
substitution Hamiltonian response, not a complete tournament invariant.

For a candidate H, introduce edge bits `x_ij=1 iff i->j`, i<j. On i<j<k the
cyclic-direction bits are `(x_ij,x_jk,1-x_ik)`. Thus require

```
positive triple in H:     x_ij XOR x_jk=0,
                         x_ij XOR x_ik=1;
negative triple not H:    NAE(x_ij,x_jk,1-x_ik).     (5.1)
```

Positive XOR components reduce to phase bits; the signed not-all-equal
clauses are indispensable. For example H containing precisely `014` and
`023` on five vertices leaves six consistent positive phase bits but has no
realization once the negative triples are enforced. This is an inherited
Session36 hostile, independently replayed here.

Let `R(H)` be the solution set of (5.1), with its natural relabelling action
by the hypergraph automorphism group. Then the complete relation-first
enumeration is

```
a_n = sum_[realizable unlabelled H on n vertices] |R(H)/Aut(H)|
    = sum_[H] (1/|Aut(H)|) sum_(g in Aut(H)) |Fix_R(H)(g)|.  (5.2)
```

**Proof.** Formula (5.1) specifies exactly all labelled tournament realizations
of H, with no duplicates. An isomorphism between two such realizations
preserves their cyclic-triple set, so it belongs to `Aut(H)`. Conversely an
`Aut(H)` map carrying one edge solution to another is a tournament
isomorphism. Every tournament has one hypergraph-isomorphism type. Sum the
resulting disjoint realization orbits and apply Burnside. QED.

The group action on directed bits is a **signed** edge permutation: relabelling
can reverse the increasing-index reference orientation. On consistent parity
components it becomes an affine signed permutation of the phase coordinates.
Discarding these signs or the stabilizers invalidates the formula. Converse
is not an extra group generator unless anti-isomorphism is expressly the target.

The binary C3 hypergraph is not the GF(2) even-graph shadow from the old
cut/cycle count. Those quotient constructions have different fibres.

The smallest visible unlabelled loss is order four: `1 -> C3` and `C3 -> 1`
have the same labelled cyclic triple and h=3 but are nonisomorphic (source
versus sink). Thus three unlabelled C3-shapes lift with multiplicities
`1,2,1` to the four tournament classes. At order five nine shapes lift with
multiplicities `1,3,2,1,1,1,1,1,1` to twelve classes. Session36 also supplies
a strong order-seven nonisomorphic, nonconverse pair with identical C3 fields:
`C3[1 -> C3, TT_2, 1]` versus `C3[(1 -> C3)^op, TT_2, 1]`.

Equation (5.2) is not claimed to outperform direct tournament enumeration.
It makes the user's “binary intransitivity in the relation itself” into an
exact finite family index with its information-restoring sidecar.

## 6. Metagraph connection contract and boundaries

| Source -> target | Preserved | Lost / necessary sidecar | Decisive test |
|---|---|---|---|
| T -> weighted pair-free core | Entire isomorphism class | Nothing if quotient orientation, sizes, and `Aut(Q)` action are kept | Reconstruct, canonicalize, and compare all n<=5 classes |
| T -> C3 hypergraph | Every binary triangle-intransitivity event; full substitution response | Interval orientation, component order, and generally converse; restore one realization orbit under `Aut(H)` | Order-four source/sink hostile, plus formula (5.2) |
| T -> signed addressed curl | All arc data modulo the known SCC cut-rotation fibre | Lost cut position; scores restore all arcs by (2.1) | Every labelled edge through n=5 |
| Path-marked T -> tiling mask | T and its selected Hamiltonian-path gauge | Forgetting the path requires correct orbit/fibre multiplicities | THM-1382 shows star quotient is not class quotient |
| T -> scores/c3/electrical energy | Score field or exact scalar triangle statistic; line-flux constraints | Placement and relational data; no all-order reconstruction | THM-4739's order-six ambiguity and the class fibres above |

The metagraph's vertices are tournament isomorphism classes, while a merged
metagraph may additionally identify converse. These are different target
sets. Bare blue/black port degrees, line energies, or unsigned C3 data cannot
by themselves supply A000568. Either (4.2)'s full weighted core or (5.2)'s
realization orbit restores enough information; the same canonical class can
then be routed into whatever accurately defined metagraph family is used.

**Characteristic-two guardrail.** Equality generated by local moves is an
incidence relation. Do not replace its quotient by the kernel of a mod-two
Laplacian without proof: for the connected graph K_(2,3), the incidence
matrix has rank four over GF(2), while `D D^t` has rank two. The Laplacian
therefore has three-dimensional kernel although the component-constant
space is one-dimensional. The extra vectors are not extra classes.
`verify_growth_information.py` checks this hostile directly; no positive
real-energy argument transfers automatically to characteristic two.

A single tournament has `binom(n,2)` binary choices before quotienting.
A fixed Hamiltonian path removes n-1 free choices and leaves
`binom(n-1,2)` tiling coordinates, but not one tiling per class. There are
`h(T)/|Aut(T)|` tilings per unmerged class (automorphisms act freely on ordered
Hamiltonian paths). Consequently presentation entropy, class entropy, and
path entropy should not be interchanged. This is another precise reason the
old count factor is not an unlabelled family factorization.

## 7. Independent finite certificates and stopping scope

Run in the deliverable directory:

```
python3 verify_information.py
python3 -O verify_information.py
python3 verify_core_information.py
python3 -O verify_core_information.py
python3 verify_growth_information.py
python3 -O verify_growth_information.py
```

All three programs use explicit `require` checks, not assertions disabled by `-O`.
Normal and optimized output files agree byte-for-byte. They use only the
Python standard library and do not import prior computation code.

| n | Labelled tournaments | A000568 | Unlabelled realizable C3 shapes | Possible C3 inputs tested |
|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 1 | 1 |
| 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 1 | 1 | 1 |
| 3 | 8 | 2 | 2 | 2 |
| 4 | 64 | 4 | 3 | 16 |
| 5 | 1024 | 12 | 9 | 1024 |

`verification.json` records exact realization fibres, automorphism orders,
Burnside fixed-point sums, score/curl reconstruction checks, entropy checks,
and hostiles. Of 1,024 possible C3 hypergraphs at order five, 253 have
consistent positive parity equations and 45 of those still have no realization.
Thus even full positive-cycle data require the negative-transitivity clauses.

`core_information_verification.json` independently extracts every transitive
module, verifies the maximal-module partition and pair-free quotient, and
checks that weighted-core orbits and full tournament classes coincide. At
n=5 the class counts by core order are `{1:1, 3:2, 4:4, 5:5}`.

No all-order proof is inferred from these finite counts. All-order claims
are the explicit algebraic and orbit proofs above; historical numerical
claims not listed as replayed remain attributed to their exact source.
The primary C3 theorem was directly read in its author PDF, including its
statement, Corollary 1 and proof; `PRIMARY_SOURCE_VERIFICATION.json` records
the inspected passages. It was not accepted solely on Session36 attribution.
The source checkout was never modified. Publication/integration is delegated
to the parent session's authorized checkpoints.
