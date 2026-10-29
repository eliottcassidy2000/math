# Session102 independent audit

**Status: VERIFIED proof audit; FINITE-EXACT controls.** Source pin
`3b6e3ef82382da222fc59aa02914eea72d5a105b`. No publication performed.

## Scope and inherited hazards

Anchor: signed modular recursion and signed-cut/root normalization.
Niche: native-path versus switching operation types.
Wildcard: sharp order-seven two-path switching construction.
Concept board: strong-module isotypes, cycle indices, switching cut rank,
root stabilizers, native guards, jointly ordered run ports.

Read the source AGENTS, START-HERE, research protocol and active guardrails.
Source was clean at its specified pin. Targeted processor startup was unavailable
because `.machine-id` is absent; no identity was invented. This is not a failed
mathematical test. No claim of classifying all prime tournaments is supported.

## Proof audit verdicts

1. **Signed modular G1–G4: PASS.** G1's denominator is the positive compatible
   quotient stabilizer; a quotient cycle of length l closes with sign epsilon^l.
   A fixed composite has exactly a^(l−1) lift choices. Arbitrary noncanonical
   partitions would count only partition-preserving maps; the canonical qualifier
   is load-bearing. Negative permutations have only even cycles of 2-adic
   valuation one and at most one fixed vertex, so only one recursive branch can
   remain negative. Complement-converse attachment orbits use alternating subset
   indicators and hence the even-cycle evaluation. The odd unpointed descent
   is valid because negative involutions are positively conjugate and their fixed
   vertices comprise exactly one ordinary automorphism orbit.
2. **Switching/root group consequences: PASS.** In root-source gauge, a positive
   switching permutation fixing the root has constant switching signs; an
   extended negative one has the opposite common sign off the root. Thus its
   rooted restriction is respectively an automorphism or anti-automorphism of
   the deleted tournament. Tournament automorphism groups have odd order (an
   involution would transpose a pair and reverse its arc). A Sylow 2-subgroup
   therefore acts freely in the positive case; in the extended case point
   stabilizers have order at most two. Orbit sizes give divisibility by n and
   2n, respectively. This concerns permutation groups, not sign-vector lifts.
3. **Rank-one cuts: PASS.** All 2x2 cross-sign products being positive is equivalent
   to proportional cross rows. Multiplication by each row's reference-column
   sign makes all rows identical. This gives a module after switching, not an
   ordinary module before switching. With an outside root, source normalization
   simultaneously gives the module in the root-deleted tournament.
4. **Native operation and SCC proofs: PASS.** Named-path guards forbid any
   nontrivial switch preserving the same path. One native move is not a switch
   from order three, and not converse-plus-switch from order four. Two native
   moves toggle the symmetric difference of their path edge sets. Cut degrees
   are at most four; order eight would require degree four everywhere, violated
   at a path endpoint. The supplied order-seven K4,3 decomposition is legal and
   sharp, including 448 black-black fillings. The general even/odd word bounds
   follow from degree≤2r and endpoint degree≤2r−1. SCC blue counting requires
   dual-compatible opposite component types and retains the center separately.
   The boundary-port defect destroys each old nonsingleton child module, not
   necessarily every possible new module.

5. **Pfaffian parity/multiport proof: PASS.** The rank-one cross two-form squares
   to zero, leaving zero or one cross pair according to parity. Zero-based odd
   cofactors give the positive last-root expansion. Eliminating an even interior
   yields E+U^T A^{-1}U (plus sign); the mixed root column is
   y+v(u^T A^{-1}x). For an odd interior, moving pivot a to the end incurs
   (-1)^a; eliminating its invertible even deletion retains gamma=z+b^T C^{-1}V
   and the full mixed port skew form. Every even tournament Pfaffian is an odd
   integer and hence invertible, so the odd tournament rank is exactly n−1.
   The scalar-energy hostile has Pfaffians 1 and 3 without invalidly resetting
   the root gauge. These are standard algebraic interface identities, not a
   classification or a claim of planar matchgate realization.

No proof corrections remain in these audited lanes. One implementation defect
was found and repaired by the grammar producer: Python `assert` validation
would disappear under `-O`; all load-bearing gates now use explicit exceptions.

## Independent controls (no producer imports)

- `check_signed_cuts_independent.py`: all labelled tournaments n=2..5; 31,668
  proper cuts and 16,932 normalized switching checks. Small nonmodule rank-one
  hostile: descending transitive order4, cut {0,3}|{1,2}. Operation masks tested
  through n=8, including every small exception.
- `check_switch_groups_independent.py`: every permutation of every labelled
  tournament n=3..5; signed root-stabilizer identity and Sylow divisibilities.
  The root-source+C3 witness has exactly A4 and S4 as its two permutation groups.
- `check_historical_hostiles_independent.py`: independently constructs the exact
  literal-converse decorated-C3 hostile. Unrestricted degree/adjacency backtracking
  finds full groups 81/27; all 4,096 named subset-H entries agree; full H=462915;
  67 labelled cyclic triples agree; all literal point stabilizers agree; binary
  attachment orbit counts are 176/512. Equality for every nonempty tuple follows
  by intersection of these point stabilizers. Independently constructs the
  all-strong source C3[A,D,D], A=(10,12,9,16,7), D=(2,4,9,3), checks strongness
  of children, source and partners, and obtains H=2,899,675 versus 2,906,091 for
  the two recorded run orders by subset-endpoint dynamic programming.
- `check_two_move_boundary_independent.py`: uses directed arc sets rather than
  producer bit matrices, verifies both native guards and final switching
  equality for all 512 within-part fillings; 448 black-black and 64 blue-blue.

- `check_pfaffian_independent.py`: a fractional two-by-two skew-elimination
  Pfaffian oracle, independent of the producer's recursive matching sum;
  300 mixed-root cases, 120 even multiport cases, and 640 odd pivot cases,
  with fixed seed 20261029. Normal and optimized outputs agree.
- `check_projective_group_independent.py`: generates degree-eight permutation
  groups from translation, negative inversion and nonsquare dilation rather
  than enumerating projective matrices. Full signed-permutation enumeration
  identifies orders 168/336. Centralizer counts give class sizes
  1,21,24,24,42,56; the only identity-containing class-union sizes dividing 168
  are 1 and 168, certifying simplicity of the positive group. This is a finite
  certificate for the classical PSL2(7) example, not a new group classification.

All six independent scripts were run normally and with `python3 -O`; matching receipt
pairs are byte-identical. Producer grammar, metagraph, switching and Pfaffian scripts were also
replayed normally/optimized, with byte-identical outputs retained here.
The computation is finite evidence; all-order statements rest on proofs above.

## Reproduction

From the session directory, for each `audit/check_*_independent.py` run Python
normally and with `-O`; compare JSON stdout. Producer replay receipts are
`grammar-producer-{replay,optimized}.json` and
`metagraph-producer-{replay,optimized}.json`, with the analogous switching and
Pfaffian pairs. Final metagraph audit additionally checks its complete order-seven
extremal skeleton argument: fixing P forces the bipartition, and the residual
K4,3 path has exactly two orientations. The labelled ordered counts are
5,160,960 total, 645,120 blue-blue, and 4,515,840 black-black.

## Remaining boundaries

No full prime classification, class-support metagraph connectivity, arbitrary
run-list minimal carrier, unrestricted unmarked extension recurrence, or
interchangeability of converse/switching/native signs follows. Fixed-label
switching word bounds do not cover arbitrary final relabelling. The 81/27
hostile and all-strong run-order obstruction remain distinct and exact.


## Historical attribution guard

The switching point-stabilizer/Sylow divisibility proof is an elementary
reconstruction. Existing source THM-4531/THM-1460 and the classical
Babai–Cameron switching-group classification take precedence; the new report
must not present a weak divisibility corollary as a new group classification.

Final switching companion text audit passes. Its two primality hostiles are
also independently tested: cyclic5 is ordinarily prime with ten nontrivial
signed-cut sides; source+Paley7 has its seven-vertex ordinary module and no
nontrivial signed cut. The classical attribution was checked directly against
[Babai–Cameron's primary PDF](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v7i1r38/pdf):
the introduction explicitly attributes the full abstract realization criterion
to Theorem5.2 and the prescribed-action semiregularity criterion to Theorem5.1.
The report correctly keeps these stronger classical results separate from its
reconstructed point-stabilizer/divisibility proof.

## Final Cayley/syndrome/parity extension audit

**PASS.** This section updates the earlier provenance boundary: order-eight
path-generator injectivity and the resulting two-native-step obstruction are
inherited from `tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md`,
§8 Theorem3, not new session theorems. The complete order-seven legal skeleton,
its colour/syndrome action, and the explicit representative-level lift loss
are retained as the scoped refinements.

The Cayley graph is the existential aggregate across every tournament in a
labelled switching class. Tree-gauge normalization gives one representative
per directed path, not permission to switch representatives during a native
word. Its generator differences span the even-edge space image; cut parity
then proves the stated component and bipartiteness classification. The
transitive-source no-lift examples correctly exclude all competing paths
when generator injectivity applies. No native connectivity follows.

The nine-bit syndrome has three disjoint two-bit supports. Its eight fibres
are all64-element sets; the blue condition is precisely the zero fibre.
Exactly four permutations preserve the unordered pair of undirected paths.
The explicit cross-cut corrections form a group action because the cut is
preserved and the reversal flag is a homomorphism. This realized V4 acts
trivially on all three syndrome coordinates: no Fano/GL3 promotion is valid.

The final edge-parity condition is necessary:
`r*binom(n−1,2) = epsilon*binom(n,2)+s(n−s) (mod2)`.
Its mod4 table is correct: residues0,1,2,3 give respectively r=s, no extra
condition, epsilon=s, and r=epsilon. In particular n=2 mod4 forbids reaching
an odd-size pure cut by any native word, independently of the size bound.

Two additional standalone independent checkers were added:

- `check_native_quotient_independent.py` generates legal edges by fixing an
  ordered directed Hamilton path and all free nonpath orientations, then
  runs union-find. This differs from producer path-backtracking. At n=3 it
  verifies degree sequence six leaves/two degree-three centers and two
  four-vertex components, hence two K1,3 graphs. At n=4,5,6 it obtains
  respectively96,3,840,368,640 edges and component sizes64;
  512+512;16,384+16,384. These are complete finite labelled censuses.
- `check_syndrome_independent.py` uses directed arc sets to verify all512
  fillings, all7! possible skeleton permutations, the corrected V4 action
  and its composition law, the eight equal fibres, and the parity table.

Both new independent checkers pass normally/optimized with identical
receipts. Final Cayley producer replay passes750,053 explicit checks in
both modes; syndrome producer replay passes14,850 checks. The root synthesis
correctly distinguishes inherited results, new scoped refinements, and a
proposed organizing framework. The even-invertible-interior qualifier on scalar isotropic-port closure is
now applied in the final README, avoiding an invalid extension to odd residual
modes. Final synthesis audit passes with no pending proof or scope corrections.
The final local audit has eight independent scripts and six producer replay
pairs (14 normal/optimized byte-identical pairs in total). Publication may
retain each independent script and a single receipt; duplicate optimized
outputs and producer replay copies are local verification records, not needed
as additional mathematical artifacts.
