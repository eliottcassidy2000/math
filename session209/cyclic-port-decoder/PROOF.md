# Exact binary query complexity on the leaf atom grammar

**Status: PROVED specialization and synthesis; FINITE-EXACT independent controls.**
For a known dominance-join decomposition containing s independently oriented
C3 or cyclic-C5 factors, exactly s two-port queries, with each answer retained
only modulo 4, recover all s orientation bits. No collection of one-port
queries recovers any of them. The independent choices are distinct labelled
presentations, not distinct tournament isomorphism classes.

This connects Sessions14/15's rigid leaf grammar to Session 27's marked
responses. The general invisibility and two-port visibility results already
exist in later repository work and are **recovered dependencies**, not new
claims here. The small explicit decoder and its exact query lower bound are
the specialization proved below. No literature-priority claim is made.

## 1. Recovery and freshness

The exact source paths and blob hashes are recorded in `SOURCES.json`, pinned
to repository commit `286f5f4ddbdf25a74ecbfe74e8ef9a735a142393`.

- Session 14: `paired_black_leaf_reduction_20261020/BLACK_QUOTIENT_LEAF_CLASSIFICATION.md`
  proves that nonself-converse paired-black leaves are nonpalindromic
  singleton/C3 SCC words, modulo reversal. Mixed leaves are rigid wrappers
  around a strong self-converse deficit-two center. These are leaves in the
  paired-black **quotient**, not degree-one vertices in the labelled native
  complement graph. The quotient forgets neighbor incidence.
- Session 15: `deficit_two_centers_20261020/README.md` and its complete census
  prove that the center is uniquely S4, with H=5. Hence mixed leaves have the
  form A -> S4 -> A^op, where A is a singleton/C3 word. The Session 14 center
  question is not still open.
- Session 27: `session27_repunits_lattices_20261022/session27-response/proof.md`
  defines arbitrary labelled transitive-base attachment responses, proves
  two-port tomography of ordered two-path covers, and records the opposite
  C3 responses 11/13. Its finite module-flip statement is historical.
- Later `proof_principles_20261103/continuation/PROOF.md`, C1-C3, proves
  one-port invisibility for endpoint-balanced modules in arbitrary uniform
  ambient quotients. In particular every vertex-transitive module qualifies.
- Later `two_port_beyond_nine_20261103/algebra/PROOF.md`, A1-A3, proves that
  reversals of disjoint odd modules are distinguished by two-port responses,
  even modulo 8. This is stronger in ambient scope than the ordered-sum
  recovery below. Universal arbitrary-core two-port reconstruction remains
  open in that source.
- `two_port_unit_distance_20261103/reconstruction/PROOF.md` already gives
  the C5 response pair 89/91. These local counts are independently rechecked,
  not claimed as discoveries.

All paths above are under `05-knowledge/results/`. The inheritance board is
leaf SCC words / labelled probes / module symmetry / odd multiplicative
amplification / spectrum threshold. The corrected near miss is treating
Session 27's historical open statement as current. The canonical hostile is
a cyclic triple that ceases to be a module after a nonuniform attachment.

## 2. Definitions and an elementary invisibility proof

Let H(T) be the number of directed Hamilton vertex words. A port mask b has
bit u equal to1 precisely when old vertex u points to the added vertex x.
For two masks b,c, add x,y with x->y and the indicated old-to-new incidences.
Write F_T(b,c) for the resulting Hamilton count. Every experiment adds fresh
ports to the original core T and then discards them. Keeping earlier ports
would change the core and is not a sequence of one-port observations here.

A module M has uniform incidence to every exterior vertex. Define C3 and C5
on Z/nZ by i->j iff j-i lies in {1,...,(n-1)/2}. Both are vertex-transitive,
and their converses have the same labelled vertex set.

**Lemma 1 (recovered special case).** Reversing a C3 or C5 module preserves
H and every one-port response, in any ambient tournament.

Proof. Every binary coloring of a C3 has two equally colored vertices.
Their transposition is an antiautomorphism preserving the coloring. For C5,
complement the marked set if needed so it has size at most two. Empty sets
admit every reflection; {a} is fixed by i -> 2a-i; {a,b} is fixed by
 i -> a+b-i. These reflections reverse every C5 arc and preserve the marked
set. Extend the selected reflection by the identity on the exterior and the
port. Uniform module incidences and preserved port colors make it an actual
isomorphism of the two extended tournaments. With no port any reflection
works. Thus the counts agree. QED.

This is a short direct proof, not a replacement for the stronger all-run
endpoint-balanced theorem. Its color argument does not automatically extend
to every odd cyclic tournament: at order 7 some binary subsets have no
reflection. That failure of this proof is not a failure of count equality.

For s disjoint such modules, reverse them successively. All2^s labelled
variants share the entire one-port table. The isomorphism used may depend
on the queried mask; therefore this does not identify the labelled cores.
The variants are also mutually isomorphic as unmarked tournaments, so this
is not a multiplicity claim in the converse-merged metagraph.

## 3. Local two-port constants and dominance amplification

Exact direct enumeration gives

    F_C3(1,2)=11,       F_(C3^op)(1,2)=13,
    F_C5(1,2)=89,       F_(C5^op)(1,2)=91.               (1)

These inherited constants use the mask convention above. `verify.py`
independently checks each by subset endpoint DP and by all 5! or 7! vertex
permutations, respectively. Thus (1) is a finite exact ingredient, not an
assumed pattern in larger cycles.

Consider a known ordered decomposition

    T_epsilon = B_1 -> ... -> B_m,

where selected blocks are C3 or C5 with known local labels, but each is
independently replaced by its converse according to an unknown bit epsilon_j.
All other blocks are fixed arbitrary nonempty tournaments. The locations,
local label maps, and unflipped orientations are known. The only unknowns
are these s bits; this is not an arbitrary reconstruction problem.

H multiplies under dominance joins, since a Hamilton word must visit the
blocks in their order and traverse each block in one interval. For selected
block j, put

    q_j = product_(i != j) H(B_i).

This is a known positive odd integer, independent of all orientation bits.
Oddness follows from tournament Hamilton-path parity; alternatively it is
direct for the leaf atoms and fixed factors used in the controls.

Place two ports inside the jth position: every earlier block points to both
ports, both ports point to every later block, and their internal masks on
B_j are 1 and 2. This is a legal marked transitive-base query. Its resulting
count is exactly

    F_j = q_j (a_j + 2 epsilon_j),
    a_j=11 for C3, and a_j=89 for C5.                  (2)

Indeed the extended tournament is still an ordered sum, now with the local
five- or seven-vertex extension as its jth block. No other unknown block
orientation affects its factor H.

Since q_j is odd,

    F_j - q_j a_j = 2 epsilon_j mod4.                 (3)

The division implicit in extracting epsilon_j is only division of the known
even integer representative 0 or2, not division by a nonunit in Z/4Z.
Equivalently, epsilon_j=0 iff F_j=q_j a_j mod4. Equation(3) is one binary
observation of exactly the jth unknown bit, with no interference.

## 4. Exact oracle complexity

**Theorem 2.** For the family above with s>=1:

1. Arbitrarily many zero-/one-port queries, even with exact integer answers,
   cannot distinguish any of its 2^s labelled variants.
2. The s nonadaptive localized two-port queries(2), retaining only H modulo 4,
   reconstruct the orientation vector.
3. Any deterministic strategy using responses only modulo 4 requires at
   least s queries in the worst case, even if it may adapt queries or use
   more than two ports in each query.

The first claim is Lemma 1. The second is(3). For the third, every queried
extension is still a tournament, so H is odd. Each response modulo 4 therefore
has only two possible values, 1 or3. A depth-r adaptive decision tree has
at most 2^r leaves; separating2^s candidates needs r>=s. Thus the presented
algorithm is optimal for this declared oracle. The same lower bound holds
for a zero-error randomized strategy with a hard worst-case query cap.
No lower bound is asserted for queries returning unrestricted integers.

For s=0 no queries are needed. Computing q_j, selecting the known masks,
and evaluating the oracle are separate costs. An oracle query count is not
an algorithm for cheaply computing H, and the theorem is not a general
complexity lower bound for tournament reconstruction.

## 5. The complete one-port fibres over leaf presentations

The previous decoder assumes the decomposition is known. On the entire
Sessions 14/15 leaf family, the full one-port table itself supplies the missing
sidecar, up to exactly the triangle bits.

**Theorem 3 (PROVED with the FINITE-EXACT order-four ingredient).** Let T be any labelled tournament representing a paired-black
quotient leaf, and let r be its TOTAL number of triangle SCCs. Among all
labelled tournaments on the same vertex set, its full one-port response
fibre consists of exactly the2^r independent reversals of those triangles.
Thus after observing the full one-port table, exactly r additional two-port
queries modulo 4 are necessary and sufficient to reconstruct T.

Proof. The recovered C4 of `proof_principles_20261103/continuation/PROOF.md`
proves that equality of complete one-port tables forces equality of the
ordered SCC partition, as actual labelled subsets. It uses the positive
weighted-cut insertion formula: differing arcs form a positive circulation,
and therefore cannot cross a strong-component boundary. Its conclusion is
stronger than equality of the list of SCC sizes.

By Sessions14/15, T has only singleton/C3 SCCs, with possibly one strong
four-vertex S4 center. On three fixed labels the two strong tournaments
are precisely the two orientations of a directed cycle. Singleton blocks
have no choice. If the four-vertex center occurs, its 24 labelled strong
orientations have 24 distinct full one-port tables. This last assertion is
FINITE-EXACT: `verify.py` sweeps all 64 four-vertex tournaments, tests strong
connectivity by every-start reachability, and checks all 16 port masks for
each of the 24 strong cores. No classification of larger strong cores is
being inferred from this check.

To read the center table from the ambient table, localize a single port
between the preceding and following SCCs. The response is its local value
times 3^r, independently of all triangle orientation choices. Hence equal
ambient tables force equal center tables and the center orientation is fixed.
Conversely all independent triangle flips preserve the full ambient table
by Lemma 1. These are exactly 2^r distinct fixed-label possibilities.
The decoder and lower bound now follow from Theorem 2. QED.

This is a sharp SECOND-STAGE query result. Obtaining the entire first-stage
one-port table may require 2^n queries; no efficient total-query reconstruction
algorithm is asserted. The scope is complete fibres over known leaf
presentations, not all fibres of arbitrary tournaments.

## 6. The concrete Sessions 14/15 and numeric15 connections

A nonpalindromic rigid word containing r triangles and t singletons has

    order = 3r+t,       H = 3^r,       |Aut| = 3^r.

When it represents a paired-black quotient leaf, all its 2^r fixed-label
triangle orientation variants are indistinguishable by every one-port
experiment. Exactly r two-port H-mod4 queries recover those bits.
Their common unmarked leaf node and its complement wiring are not being
reconstructed by this experiment.

For a mixed leaf A -> S4 -> A^op, suppose A has r triangles and t singletons.
Then

    order = 6r+2t+4,       H = 5*9^r,       |Aut| = 9^r.

There are 2r triangle blocks, so the same second-stage task takes
exactly 2r queries (none when r=0). Independently choosing their labelled
orientations does not destroy self-converseness as an isomorphism-class
property: each triangle is self-converse up to relabelling. The quotient
leaf remains the same. These formulas are SCC consequences of the inherited
classification, not a new proof of its completeness.

For the numeric15 lane, k copies of C5 have H=15^k on 5k vertices. Their k
orientation bits likewise take exactly k localized H-mod4 queries. This
is a measurement statement, not an efficient realization bound. Session 208's
independently verified order14 atom has H=15^6, much smaller than the30
vertices of six separate C5 blocks. It has a different strong-component
structure, and our decoder does not assume that this atom splits into C5s.
The Sessions 14/15 leaf spectrum consists only of powers3^r (r>=1) and
5*9^r (r>=0), so no value 15^k with k>=1 occurs at a paired-black leaf.

**Threshold consequence and limit.** The bridge explains precisely which
marked orientation information a scalar/one-port observer drops, and how to
recover it on the leaf grammar. It supplies neither a new Hamilton-count
value nor an all-order construction for every sufficiently large odd value.
The eventual threshold problem remains open. It also supplies no count of
nonleaf bridges, no native-metagraph connectivity theorem, and no assertion
that raw H or query complexity determines graph wiring.

## 7. Hostiles and reproducibility

- Module hypothesis: extend C3 and C3^op by mask 1. Both four-vertex cores
  have H=5, but a further port of mask 10 gives 11 and 13. The old triple is
  now nonuniform to the first port, hence no longer a module.
- Known decomposition: localizing a query inside one SCC requires its
  location and labels. Without that sidecar, the theorem is not a decoder
  for an unknown arbitrary tournament.
- Isomorphism quotient: all chirality variants above are isomorphic;
 2^s cannot be read as a number of metagraph class nodes.
- Finite checks support the proof rather than extrapolate it. The verifier
  checks every ambient quotient on up to five vertices with a C3/C5 module,
  every port mask, and explicit isomorphisms; these are bounded checks.

Run `python3 verify.py` and `python3 -O verify.py`. Both regenerate
`RESULT.json` with the same content. The code uses explicit exceptions,
not removable assertions. The quotient controls enumerate 677,320
isomorphisms in total. Local constants are checked by independent
permutation enumeration, and the localized decoder is checked by full
extension DP on all variants of three stated SCC words. All 24 labelled strong S4
centers are separately checked for one-port faithfulness. Source documents
are read-only dependencies; the verifier needs no repository checkout,
network, optional package, or original scratch path.
