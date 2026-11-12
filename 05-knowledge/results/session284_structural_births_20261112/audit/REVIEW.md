# Independent Session284 review

**Verdict: PASS, with minor terminology/link repairs applied.**

Scope: mathematical review of the new architecture, births, extension,
wildcard and cross-family bridge notes; standalone implementation controls.
This is not a re-audit of all inherited repository theorems, a full formal
proof, or a census at orders50/200. No all-tournament finite-core classification
or Chan–Pak conclusion is certified.

## 1. Mathematical review

- **Architecture:** the canonicalized extension recursion is exhaustive;
  Burnside counts marked children, and deletion-vertex orbits give the exact
  forgetting map. The W-parent lower bound follows from the explicit rigidity
  lemma plus inherited THM-4169. The information bound is for injective
  fixed-length addresses, not a prohibition on concise algorithms. Repaired
  naming from 'canonical augmentation' to extension-and-canonicalization
  avoids implying an isomorph-free construction-path algorithm. A doubled
  slash in a reference was removed. No numerical statement changed.
- **Births:** all modules of a marked transitive-deletion tournament are the
  stated intervals; endpoint conditions are necessary and sufficient. The
  full module count handles the transitive and empty-word boundaries.
  Alternating selected subsequences classify all induced prime subsets.
  The shifted-index bijection establishes the binomial inventory; the
  two-state matrix gives its recurrence. Rigidity and explicit even-order
  prime extension arguments are sound. The bounded-prime-seed obstruction
  does not rule out a finite parametric family rule.
- **Extension:** the signature/run data are lossless when the old subset is
  retained. Positive run lengths leave skeleton core rank unchanged. The
  newborn-middle-vertex hostile correctly prevents a false run-count formula.
  Fixed-template Burnside includes possible action kernels. The weighted
  root-forgetting identity keeps actual child automorphism orders. The
  walk/Stirling semantics correctly orders the internal path-cover blocks.
  The labelled random-growth proof uses independent forest edge differences;
  newborn-edge expectation and total-variation coupling are valid. The
  Poisson limit has its necessary m=n-b→infinity regime. External-signature
  onset k+2^k is not confused with ordinary VC onset2^k.
- **Wildcard:** canonical strong-component/C3 decomposition gives the exact
  hereditary class. Ordered minimax birth depth is necessary as well as
  sufficient; one cannot use an unordered Kraft capacity. The added greedy
  dyadic allocator is equivalent to ordered prefix-code feasibility: interval
  refinement and unary-node suppression preserve the needed inequalities,
  and earliest aligned placement is coordinatewise optimal. The cyclic orbit
  recurrence separates nonconstant and constant triples and counts full
  stabilizers only for the ternary tower. Pointwise stabilizer occupancy is
  simultaneously sharp at every level. Its free action on path groups with
  fixed old neighbours of newborns proves the divisibility bound without
  assuming those actions are automorphisms of the extended graph. Sharpness
  of the group floor is not promoted to Hamilton-count gcd equality.
  General odd-prime towers explicitly retain only a subgroup.
- **Bridge:** the complete intersection criterion r(b)≤1 follows from the
  W5 obstruction and the normal form. The padded-word change-point bijection
  and total-size generating-function shift are correct. Counts are marked.

## 2. Independent controls

`independent.py` imports no producer code. It uses literal adjacency matrices,
permutation automorphism tests and explicit subset/partition searches.

- 33,864 fixed-parent batch cases: every labelled parent of orders1–4 with
  two newborns, checking survival distribution and expected rank for each
  parent individually (not only a rank-pooled ensemble).
- W5 has trivial automorphism group; exactly20 prime attachment patterns
  form14 unmarked prime child classes. The forgetting fibres are eight
  singletons and six pairs.
- An explicit81-element permutation group of the nine-tower verifies
  pointwise-stabilizer floors on all512 vertex subsets.
- Direct weighted root-forgetting checks for(parent order,batch size)
  (1,2),(2,2),(3,1),(3,2) return the expected class totals2,4,4,12.
- Actual membership in the ordinal/C3 class is tested by recursive literal
  bipartition/tripartition search on all127 apex words through length6,
  independently of the descent formula.

Normal and optimized runs agree. All six documented commands were replayed in normal and optimized modes.
Producer scripts were separately replayed
against their saved receipts; wildcard's final source includes the literal
9! permutation census and 21,844 greedy-versus-split-DP cost lists. Their stated finite control universes are not used
as evidence for unproved all-order extrapolation: the notes supply the
all-order proofs.

## 3. Inherited boundaries

Canonical chain decomposition, THM-4169 prime-parent classification,
THM-3121 substitution semantics, prior modular species, and ordinary VC
onset are inherited dependencies. This review checks their new uses and
small controls, not their entire original proof/computational archives.

The complete **marked** grammar is not an efficient complete unmarked atlas.
The complete family grammars are not all tournament families. Exact object
generation is not absolute Hamilton-count selection. These distinctions
remain explicit in the deliverables and must survive the user summary.

`SOURCE_HASHES.json` pins the reviewed source/receipt bytes after final
presentation repairs. Its absence means the final snapshot is not yet made.
