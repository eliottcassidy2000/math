# Inheritance, connections, and information retained

## Portfolio and live board

Anchor: one-vertex tournament growth, with a full unmarked group theorem
and a concrete A000568 enumeration component.

Niche: the growing-palette limit of the old stabilizer recurrences.
The old fixed-palette formulas were recovered before making claims about
new work; the limit concerns the residual group, not just its order.

Wildcard: rigid forbidden-tournament families and the Erdős–Hajnal
interface. This lane overtook the initial analogy when the Paley-eleven
deletion test produced a usable all-size guard.

The live board is **regular deletion / decorated action / fixed subspaces /
named ports / weighted transitive sets / source-sensitive limits**.
The LRC(14) anchor remains OPEN: no result here supplies its missing weighted
target-point bound.

Inheritance pass:

- Closest proved mechanism: weighted stabilizer-necklace recurrence S1 in
  [STABILIZERS.md](../stabilizer_necklace_kernels_20261028/STABILIZERS.md),
  with its ordinary versus fixed-converse-centralizer distinction.
- Canonical hostile: the twelve-vertex 81/27 pair; its ordinary child
  decorations differ although its labelled cyclic triples agree.
- Corrected near miss: the global form of
  [THM-4803](../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md), as repaired in
  [Session68 CLASSIFICATION](../session68_tournament_composition_classification_20261026/session68-classification/CLASSIFICATION.md).
  Local group formulas survive; equal cyclic triples do not determine
  the global ordinary automorphism group.
- Least-used relevant sidecar recovered: the actual point stabilizer
  and the two roots of a clone extension. These are cheaper than a
  complete child deck, and regularity makes them sufficient here.
- An independent hostile added this session: a rigid extension of
  C3[C3] creates R5. H-freeness cannot be inferred from rigidity.

## Source recovery and what changed

| Source | Recovered mechanism | Use and boundary in this packet |
|---|---|---|
| [Session68 classification](../session68_tournament_composition_classification_20261026/session68-classification/CLASSIFICATION.md) | Canonical strong-module tree, actual quotient action, stabilizer of ordinary child types | Supplies full wreath actions. R1/R2 concern arbitrary regular parents, not only a visible wreath subgroup. |
| [Stabilizer necklace kernels](../stabilizer_necklace_kernels_20261028/STABILIZERS.md) | Cost and stabilizer polynomial; conditional future stabilizers | Rooted polynomial is inherited. R2 subtracts the point-stabilizer orbit to give unmarked children. |
| [Session24 palette dynamics](../exception_palette_dynamics_20261021/session24-palette-dynamics/PALETTE_DYNAMICS.md) | Ordinary free-color orbit polynomials for C3,R5,P7; fixed-Q rarity | Recurrences are not rediscovered. Growing Q gives a new scoped limit and error bound in R4. |
| [Session107 colored recurrence](../session107_quasirandom_arithmetic_20261029/colored/COLORED_RECURRENCE.md) | Named colors, odd-cycle Burnside sum, DU=2UD+Q2^n I; normalized Weyl relation | Full A000568 already has an exact global count formula. R2 supplies local duplicate-free augmentation in the regular-card sector, not a replacement global scalar formula. |
| [Session108 palette symmetry](../session108_sigma_symmetry_20261030/tournaments/PALETTE_SYMMETRY.md) | Palette permutations are an additional quotient, with joint fixed data | Random ports here are fixed individually. Allowing their permutation changes the group and requires the palette sidecar. |
| [Growing arithmetic attachments](../growing_attachment_families_20261029/README.md) | Explicit finite-window arithmetic attachment bounds | An asymptotic IID palette for a fixed bank does not automatically justify Q growing with n. The simultaneous error remains a separate obligation. |
| [Session109 Pisano symmetry](../session109_pisano_symmetry_20261030/README.md) | Ramified lift and source clock retain an extra field-valued coordinate | A color parity can be enough for this specific stabilizer query; this does not authorize discarding a phase in arithmetic returns. |
| [Wedge/descent/168](../wedge_descent_168_20261029/THEORY.md) | Fixed-port nilpotence and its re-marking failure; literal versus switching action | R1 obtains an intrinsic fixed port by regularity. It does not extend a fixed-port descent argument to arbitrary re-marking. |
| [Tournament open-problem routes](../tournament_open_problem_routes_20261029/PROOFS.md) | Exact weighted transitive-set substitution and the failure of a scalar Hamiltonian proxy | R6 preserves a forbidden-pattern guard and constructs rigid examples. A lower bound uniform over all allowed prime atoms is still absent. |
| [Native 81/27 control](../overnight_rank_return_20261029/atom_precision.py) | Explicit labelled pair and native automorphism tests | The new verifier reconstructs the pair without importing that code, then exhausts its regular binary masks. |

Source hashes in the manifest pin the inherited Git versions, including
historically future-dated folder names. Folder dates are provenance labels,
not a reason to override the status hierarchy.

## A common decomposition, three different calculations

The same canonical tree supports three distinct operations:

1. **Groups:** internal automorphisms form a direct product, extended by
   the quotient permutations preserving actual child types.
2. **Counts:** fixed-point counts use the cycle index; trivial stabilizers
   require intersections of fixed-color subspaces, or equivalent subgroup
   information. Unmarking additionally requires the clone correction.
3. **Extremal sets:** maximum transitive size is a maximum of sums of
   child weights over transitive quotient subsets.

These operations share the ports and the LCA tree, but their numerical
outputs are not interchangeable. R6 is a particularly explicit test:
the stabilizer can be made trivial while the third operation retains
the same power-law obstruction. R5's hostile pair shows that recording
cyclic triples loses information needed even for the first operation.

The random-port threshold is an explanation by codimension. For a
nonidentity permutation, x_i=x_(g(i)) forces b-c(g) independent equalities.
There are linearly many bottom blocks, so a defect of probability Q^-d
becomes visible at n/Q^d of constant order. Higher-node coincidences
require longer blocks and vanish in that regime. The resulting small
cyclic factors are attached to independent bottom locations.
This is why an order-21 base leaves order-three defects at its threshold.

The explicit rigid constructions instead use a parity character on words.
For a base cut S, the total character sum is (b-2|S|)^k, hence the two
weights (b^k +/- (b-2|S|)^k)/2. A trivial base set stabilizer and unequal
weights keep the two child orbit types distinct at every node. The
coincidence of the imbalance 3^k in the R5 and P7 constructions has this
specific algebraic cause: 5-2=7-4=3.

## Connection ledger

| Source -> target / map | Preserved predicate | Information lost / required sidecar | Cheapest decisive test |
|---|---|---|---|
| Regular T + subset S -> unmarked child T+v | Full literal automorphism group equals G_S | Root may have two regular-card candidates; retain oriented clone edge or apply R2 correction | All cuts of every regular parent through seven, plus C3[C3] |
| Parent group -> Burnside extension count | Number of unmarked children with a regular card | Scalar group order loses action on subsets and vertices | Compare native child classes with subset-orbit minus vertex-orbit counts |
| r fixed ports -> Q=2^r named colors | Pointwise port stabilizer, all incidence bits | Permutable ports need a palette action; arithmetic rows need their feasibility relations | Q=2,3,4 base coloring exhaustion |
| Recursive action -> bottom Poisson defects | Bounded event probabilities in stated scaling | Full tails, weighted moments, exact higher-node action disappear | Exact recurrence versus rational union bound through depth four |
| P7 word -> parity coloring | A constructive trivial stabilizer; signed weight imbalance | Other colorings and arithmetic phase data are not retained | Native groups at child orders eight and fifty |
| One new vertex -> forbidden-pattern test | H-freeness when H and all H-v are sufficiently large prime tournaments | H-free old core alone is insufficient | All modules of P11 and its cards; rigid R5-creating ten-vertex hostile |
| Substitution tree -> transitive size | Exact weighted maximum, not just group order | Unbounded prime-atom weight inequalities remain | Explicit all-zero product transitive set attains the upper bound |
| Probability bound -> LRC/Collatz target | No unconditional transfer asserted | Actual target weights, source values, legal words, clocks, and return coverage required | The existing nonroot 5x+1 cycle and target-point weighted controls must remain live |

## Why the 168 symmetry does not conflict

For the literal tournament P7 plus a source, regular-parent theorem R1
fixes the source and leaves Aut(P7) of order 21. A switching-class
automorphism is allowed to change the orientation representative by
switching; its order-168 group acts on a different object. The missing
coordinate is the selected representative/root. Quotienting it away can
increase apparent symmetry, whereas adding a fixed port intersects
stabilizers. This is a mathematical distinction between operations, not
a numerical analogy between 7,21,168.

## Primary literature and novelty discipline

- **CITED:** [Melody Chan, The distinguishing number of the direct product
  and wreath product action](https://arxiv.org/abs/math/0601414).
  Distinguishing-color orbit methods for wreath actions are established
  literature. The inherited repository recurrences belong to this framework.
- **CITED:** [Alon, Pach, Solymosi, Ramsey-type theorems with forbidden
  subgraphs](https://web.math.princeton.edu/~nalon/PDFS/aps4.pdf), Theorem1.2.
  The graph and tournament Erdős–Hajnal formulations are equivalent.
  We use the transitive-subtournament formulation and a scoped substitution
  construction; historical open-status remarks elsewhere in that old paper
  are not being treated as current facts.
- **CITED through the audited repository structural route:** the canonical
  decomposition sources linked in Session68. The full decomposition theorem
  is not reproved here. The new standalone regular-parent proof does not
  depend on it.

The elementary proofs in this packet certify their statements; they do not
establish that these statements were previously unknown. Neither a sharper
random-color limit nor a rigid extremal family supplies an Erdős–Hajnal
lower bound for arbitrary forbidden tournaments.
