# Closed counts, operation responses, and the common label boundary

Status: **PROVED in the stated scopes; independent
[audit PASS](../six_color/INFORMATION_AUDIT.md).** The accompanying
exhaustive computation is FINITE-EXACT, separately from the all-order proofs.
No new literature-priority claim, H algorithm, metagraph enumeration speedup,
or Hadwiger–Nelson bound is asserted.

## 1. Inheritance and the live comparison

The older `07-reflections/information-geometry-of-tournaments.md` and
`07-reflections/directed-graph-metagraph.md` supplied questions, not dependencies.
The operative recovered mechanisms are:

- THM-076, `01-canon/theorems/THM-076-walsh-ocf-factorization.md`, gives signed
  even-path Walsh coefficients. THM-071,
  `01-canon/theorems/THM-071-walsh-fourier-diagonalization.md`, identifies the
  degree-two support with adjacent edge pairs. Its stated general-scope limits
  are not silently removed here: the needed degree-two formula is proved
  directly below for every order, without OCF.
- THM-4739,
  `01-canon/theorems/THM-4739-kirchhoff-and-hodge-laws-of-the-tournament-metagraph-intransitivity-currents-second-order-determination-and-the-consecutive-start-port-ledgers.md`,
  separates scores, full signed triangle circulation, and compressed port
  statistics. Its finite reconstruction claims do not identify the full tensor
  with per-vertex cyclic-triangle counts.
- [Session36 metric curl recovery](../../session36_information_curl_contexts_20261023/lattice/METRIC_CURL_RECOVERY.md),
  §§2,8, already proves that a full labelled curl fibre consists of the cyclic
  rotations of the ordered strongly connected components, has size equal to
  their number, and is resolved by one named vertex score. Its
  [realizability decoder](../../session36_information_curl_contexts_20261023/root/CURL_REALIZABILITY.md)
  and [unlabelled necklace transform](../../session36_information_curl_contexts_20261023/root/UNLABELLED_CURL_NECKLACES.md)
  are recovered results, not discoveries here.
- THM-4731,
  `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`,
  retains marked Hamiltonian paths when complementing non-path arcs. The
  [Session66 tc scope correction](../../session66_root27_aliquot_graph_interfaces_20261026/session66-tc-scope/REPORT.md)
  distinguishes H from H/|Aut| and an induced boundary action from an inflated
  gadget automorphism group. No forbidden-H result is exported to tc here.

Anchor: predict operations on a labelled relation. Niche: recover the exact
quotient used in complete tournament counting. Wildcard: combine two lossy
records through their common labels. The live board is signed curl / complete
flip response / marked path incidence / separate versus joint isomorphism
quotients / scalar statistical projections. The canonical hostile is a closed
count preserved by a quotient whose next marked operation is not preserved.

## 2. The relation and its coordinates

Fix V=[n], E=binom(V,2), m=|E|. For e={i,j}, i<j, write x_e=+1 for i→j and
x_e=−1 for j→i. This is a complete antisymmetric binary relation; transitivity
is an additional property, not part of its definition. Write f_ij=−f_ji for
the same intrinsic relation and

    K_ijk=f_ij+f_jk+f_ki.

For an oriented triple K is ±3 on a directed cycle and ±1 on a transitive
triple. Thus binary intransitivity is in the relation itself. Changing labels
or the reference orientation of an edge changes coordinates, not that relation.
Physically reversing all arcs sends f to −f; this is a different operation.

For a named subset F⊆E let T⊕F reverse exactly those arcs. Define the complete
operation response

    R_T(F)=H(T⊕F),        F⊆E.

Every named arc is available, including any previously chosen base-path arcs.
This is the FULL arc cube, not the tiling cube with a Hamiltonian path held
fixed. MISTAKES.md's tiling-complement correction remains load-bearing.

## 3. Complete flip responses lose exactly global converse

**I1.** For n≥3, R_T=R_U iff U=T or U=T^op, with labels fixed.

For two distinct adjacent edges e,f of K_n, let ε_ef be the product of their
reference signs when traversing their two-edge path. Reversing that traversal
changes both signs and leaves ε unchanged. Put

    κ_n=2(n−2)!/2^(n−1) > 0.

The degree-two Walsh coefficient of H is

    Hhat({e,f})=ε_ef κ_n          if e,f meet;
                 0               if e,f are disjoint.

Proof of the coefficient. Expand the indicator of every directed Hamiltonian
vertex word as 2^(−n+1) times the product of its n−1 arc factors (1+σ_e x_e).
A word contributes to this coefficient only if it uses both selected edges.
For a two-edge path there are 2(n−2)! such words; both traversals have the same
sign ε. For disjoint edges, reversing the traversal of just one two-vertex
block pairs the words with opposite contributions. This proof works for both
odd and even n.

For the response, orthogonality of Boolean characters gives

    2^(−m) sum_(F⊆E) (−1)^|F∩{e,f}| R_T(F)
       = Hhat({e,f}) x_e x_f.                         (1)

Consequently equal responses force x_e x_f=y_e y_f whenever e and f meet.
The ratios y_e/x_e must be constant on the connected line graph L(K_n).
All arcs therefore agree, or all are reversed. Conversely H(T)=H(T^op) by
reversing every Hamiltonian vertex word, and this equality commutes with every
fixed flip F. This proves both directions. For n=2 H is constant and the two
labelled tournaments are exactly one converse pair; n=1 has one tournament.

**I2. Sharp binary compression of this response.** Choose any spanning tree
Q of L(K_n). Divide each left side of (1), for {e,f}∈E(Q), by ε_ef κ_n.
The resulting m−1 signs are x_e x_f. Choose the sign of one edge e0 and
propagate along Q to recover every x_e. The two choices produce T and T^op.
There are 2^(m−1) converse classes, and all sign assignments on Q occur.
Thus m−1 binary features are necessary and sufficient among fixed-length
binary encodings; one additional oriented arc restores the full relation.

This is a statement about the size of the retained exact observations, NOT
about the number of H queries. A correlation in (1) can require the entire
2^m response table. There is no fast H algorithm or efficient oracle sampling
bound here. The reference orientation and chosen tree are known coordinate
data; their replacement transforms the decoder equivariantly.

## 4. Two complementary losses, and the missing alignment

**I3.** For n≥3 the pair (K_T,R_T), with its common vertex and arc labels,
determines T exactly. Indeed I1 leaves only T and T^op. Their curl tensors
are K and −K, which differ since every triple value is a nonzero odd integer.
The spanning-tree correlations may replace R in this assertion. At n=2 both
curl tensors are empty and both responses constant: the lower boundary fails.

Under simultaneous relabelling, the complete pair classifies tournament
isomorphism types. Hence its orbit count is exactly A000568. The response
alone instead classifies tournaments up to isomorphism AND converse. If a_n
is the tournament count and s_n the self-converse count, its orbit count is
(a_n+s_n)/2. The inherited full-curl quotient has the different SCC-necklace
count. These equalities specify objects and quotients; none is an enumeration
acceleration or permission to divide a labelled count by n! uniformly.

**I4. Taking the two isomorphism quotients separately destroys I3, minimally
at n=4.** Let C be a labelled directed triangle and v a singleton. Set

    T=C→v,               U=v→C,

with C internally unchanged. Reversing the uniform cut cyclically rotates
the SCC list, so K_T=K_U even with fixed labels. Since C is self-converse,
T and U are converse-isomorphic; their response objects are isomorphic.
But their score multisets are respectively (0,2,2,2) and (1,1,1,3), so T
and U are not isomorphic. The curl isomorphism and response isomorphism
cannot be chosen to be the SAME vertex permutation. Their alignment is the
missing sidecar. All tournament types through n=3 are self-converse, so a
response isomorphism already distinguishes their types; four is minimal.

This is a concrete full-boundary lesson: two records can be jointly faithful
while their separately unmarked versions are not. It concerns the actual
common relation, not merely two matching counts.

## 5. Marked operations explain why closed H is insufficient

The same-curl n=3 pair

    T: 0→1→2 and 0→2;       U: 1→2→0 and 1→0

has K_012=1 and H=1 in both cases. Reverse the SAME named pair {0,1}.
The first result is transitive with H=1; the second is a directed triangle
with H=3. Thus curl determines closed H but is not a sufficient state for
named arc interventions. This is minimal, since H is constant for n≤2.
The correct repair for arbitrary interventions is the complete response
(or its m−1 correlation code) together with one oriented arc, or the aligned
I3 pair. A single orientation bit added to curl alone need not resolve its
multiple cuts; curl does not already retain the marked wiring.

Keeping just the H-values of SINGLE arc flips is also insufficient. Even
after adjoining full curl, eight collision fibres remain at n=4 and 184 at
n=6. There are none at n=3,5 in this census; no parity theorem is inferred.
The proof of I1 needs the complete operation response, not just its radius-one
neighbourhood or the value of its quadratic component at one input.

**I5. A marked valid Hamiltonian path repairs the curl loss for the actual
tiling operation.** The map (T,P)↦(K_T,P), where P is a directed Hamiltonian
vertex word of T, is injective. Every such word traverses the ordered SCCs
contiguously and in their order. Among the rotations in the inherited curl
fibre, exactly one places the component containing the first vertex of P
first. That reconstructs T. In fact the first vertex alone suffices when
accompanied by the PROMISE that it starts a valid Hamiltonian path; an
arbitrarily marked vertex does not suffice.

Consequently (K_T,P) can reconstruct the complement-line endpoint which
reverses E\E(P), and can support the true marked metagraph incidence. A scalar
H, tc, or port count does not retain P or its compatibility. This is a recovery
statement, not an efficient computation of the whole metagraph.

## 6. Historical triage with exact hostile controls

The all-order implication H(T)=H(T^op) ⇒ zero odd Walsh coefficients is valid
in the full arc cube. The following additional claims in the March reflection
do not follow and have decisive small controls:

1. **Fisher maximum at β=0 is false.** For the finite Gibbs law proportional
   to exp(βH), I(β)=Var_β(H), and I'(0) is the third central moment under the
   uniform law. At n=3, H=1 has multiplicity 6 and H=3 multiplicity 2. Therefore
   I(β)=12 exp(2β)/(3+exp(2β))², maximized at β=log(3)/2, with I'(0)=3/4.
   At the reflection's own n=5 scope the exact derivative is 165/32>0.
   Finite partition functions are analytic; this says nothing about an
   unproved asymptotic phase transition statement.
2. **Score conditioning is not the Walsh degree-zero-plus-two projection.**
   At n=5 the latter equals 3c3. On the transitive score fibre it is 0, whereas
   E[H|score]=1. Thus degree-four content has a nonzero conditional mean on
   this fibre. Exhaustively, Var(H)=285/16, Var(E[H|score])=1935/112, and the
   degree-two Walsh energy is 135/8. The two explained fractions are 129/133
   and 18/19, respectively, not equal. Score conditioning is an orthogonal
   projection onto all score-measurable functions, a different subspace.
3. **Betti_0=1 does not imply simple connectivity or vanishing higher
   persistence.** The n=3 H=1 sublevel consists of the six transitive vertices
   of Q3, excluding its two opposite cyclic vertices. Its induced graph is a
   six-cycle, with Betti_0=1 and Betti_1=1. The induced cubical subcomplex has
   no square faces, so the same counterexample applies there.
4. **The cube's Euler characteristic is not generally 0.** A filled m-cube
   has Σ_d (−1)^d binom(m,d)2^(m−d)=(2−1)^m=1. Its one-skeleton has Euler
   characteristic 2^m−m2^(m−1); these are different carriers. No discrete
   Hessian signature alone supplies a Morse complex or its critical indices.

The historical text remains provenance. Parent integration owns the maintained
warning and MISTAKES lineage; this packet changes no historical or canon file.

## 7. Reproduction, controls, and loss ledger

Run from the repository root:

```powershell
python 05-knowledge/results/hn_metagraph_information_20261026/information/check.py
python -O 05-knowledge/results/hn_metagraph_information_20261026/information/check.py
```

[check.py](check.py) uses subset-DP H, independently checked against all vertex
permutations through n=5. It enumerates every labelled tournament through n=6
without filters; verifies all quadratic coefficients by integer FWHT, every
tree-decoder fibre, and every aligned joint invariant; directly groups full
response tables through n=4; checks shifted correlation sums exhaustively
through n=4 and at four shifts in each larger order; verifies the path-start
selector in every curl fibre; and enumerates all isomorphism quotients through
n=5. Positive controls include every converse
pair; hostile controls include the n=3 intervention and n=4 alignment failure,
the radius-one collisions, the six-cycle sublevel, and the Fisher/OCR
corrections. Normal and optimized
receipts are identical: [receipt.json](receipt.json),
[receipt_optimized.json](receipt_optimized.json). The current run has 667,959
explicit checks over 33,867 labelled tournaments.

| Map | Exactly retained or lost | Repair / boundary |
|---|---|---|
| T→full signed curl | SCC interiors and circular order; loses initial cut | One named score, or a valid path-start promise (recovered/derived) |
| T→complete labelled flip response | Relation up to global converse | One oriented arc; m−1 binary correlation features are sharp |
| T→aligned (curl,response) | Full labelled relation, n≥3 | Common label carrier is essential; n=2 fails |
| T→separate isomorphism classes of those records | Can lose their relative identification | C3→v versus v→C3 is the first hostile |
| (T,P)→(curl,P) | Full marked tiling-line state | Retain path validity and its vertex word |
| H→score conditional projection / degree-two projection | Different scalar summaries | Neither summary is full marked incidence |

The useful common framework is an operation-relative state: specify which
named attachments or changes must still be evaluable, then retain a shared
boundary relation sufficient for them. These exact tournament maps support
that principle without identifying a tournament with a Euclidean colouring
instance or claiming a new chromatic bound.
