# Adversarial review: double-triple symmetry boundary

20 October 2026. **Verdict: PASS, no mathematical blocker found.**

The frozen producer packet is
`session13-census-extension/double-triple-boundary`, with inventory SHA-256
`d49efa234ef37af5b7eed7df6ba1e5a4ebfe7679050d025384fad105a872a67c`.
All twelve inventory entries passed byte-size and SHA-256 verification; the
packet's own `SHA256SUMS.txt` also passed. No producer or repository file was
edited. This review and its new computation live separately.

## 1. Claims accepted and scope retained

The displayed proofs establish equality of full tournament and globally
sign-even graph types for the two specified embeddings on `5+3+3` with zero
or one fixed point: synchronized `C5×C3` and independent `C5×C3²`.
At degree thirteen, each predicate has exact deficit six, accounted for by
four globally even eligible-C5-copy fusion losses and two relatively even
eligible-copy orbits on globally odd graphs.

The regular-tournament module inequality `3|A|≤d` and its directed-triangle
equality decomposition are valid as stated. The single-cycle corollary
requires proper support; the packet correctly excludes a full cycle.

This is **not** a complete degree-thirteen subgroup census, a universal
higher-degree failure result, a claim about all embeddings of the abstract
groups, or a global bijection preserving every prime simultaneously.
The packet states these boundaries correctly.

The original degree-eleven/twelve counts and degree-thirteen cube totals
remain **FINITE-EXACT** computational results. The exact deficit-six proof
does not require the degree-thirteen cube totals.

## 2. Predicate-restricted matching audit

For a fixed single-five-cycle subgroup K,
`N(K)/K = C4×S_(n−5)`. The inherited relative theorem preserves the conjugacy
class of a Sylow-3 subgroup of the quotient stabilizer, hence every specified
3-subgroup containment predicate. Its use here does not require ambient
pronormality of K.

The lift from quotient containment to the original embedded H is valid:
`N(K)∩Aut(X)` contains K. If its quotient contains a conjugate of U, it
therefore contains the **entire preimage** of that conjugate. The preimage
of the displayed residual double-three-cycle U is precisely the required
K×C3, and that of independent triple rotations is K×C3². Conversely either
specified extension gives this quotient containment. Quotient conjugators
lift to normalizer conjugators.

Full objects with a marked eligible K-copy correspond to the relevant
normalizer orbits. Forgetting the mark counts conjugacy classes of
**H-extendible ambient copies of the single-cycle K**, not all abstract C5
subgroups and not all single-cycle copies irrespective of extension.
Automorphism conjugacy preserves this eligibility predicate.

With source multiplicity one, the relative count is

    R = tau_T(H) = tau_E(H) + F_E^K + Z^K.

The definition of Z uses relatively even marked-copy orbits, not bare
globally odd graph types. In the two actual escape types there is one such
orbit each. The resulting equation is sound. The packet separately records
the H-based defects; it does not confuse synchronized `Z^H=258` with
restricted `Z^K=2`.

## 3. Graph geometry and exhaustive anomaly reduction

The C5 support in a relatively even graph must be independent twins.
There are only four circulants on five vertices. Reflection has negative
edge-reversal sign on either five-cycle; multiplier two has negative sign
on K5. Extending these permutations by the identity outside the support
produces no external sign contribution. This gate was also independently
tested from the definition of sign on individual edges.

Every internal permutation of an independent-twin class is sign-even.
If an automorphism preserves the marked support's maximal twin class, an
internal permutation can conjugate its image of K back to K. It follows
that the class stabilizer is sign-even whenever the marked graph is
relatively even. Internal permutations also conjugate all single-five-cycle
subgroups in that class.

Consequently fusion or sign escape requires at least two maximal
independent-twin classes of size at least five. There are at most two at
degrees eleven through thirteen. This is the key exhaustiveness statement
that justifies a small graph computation; it is not supplied merely by
enumerating the producer's selected quotient patterns.

At eleven and twelve:

- Weights 5+6 and 5+7 have only the weight-five class eligible. A copy in the
  larger class leaves fewer than three unused points there and at most one
  three-cycle in the other class.
- Two eligible classes require 6+6, where the graph is K_(6,6) and the class
  swap is sign-even.
- Equal five-classes leave fewer than three residual vertices, so an
  extending double-three-cycle cannot fit; equal six-classes have an even
  swap. Thus sign escape is absent.

At thirteen all large-class weights are covered by
5+5, 5+6, 5+7, 5+8, 6+6, and 6+7, with the remaining vertices arbitrary.
For 5+5 the extra cycle must use all three residual vertices. Their
adjacency to each large class is uniform and their induced graph is empty
or K3. K3 supplies an odd transposition fixing the marked support, so only
the independent triple can survive relative evenness. The 5+6 and 5+7
patterns still have one eligible class and no possible class swap.

After merging equal open neighborhoods, the surviving anomalies are
exactly:

| Full graph | Global sign | Eligible-copy orbits | Contribution |
|---|---|---:|---:|
| K_(5,3) disjoint union I5 | even | 2 | fusion 1 |
| K_(5,8) | even | 2 | fusion 1 |
| K_(6,1) disjoint union I6 | even | 2 | fusion 1 |
| K_(6,7) | even | 2 | fusion 1 |
| K_(5,5) disjoint union I3 | odd | 1, relatively even | escape 1 |
| K_(5,5,3) | odd | 1, relatively even | escape 1 |

The two equal-five-class swaps reverse 25 edges. The equal-six-class swap
reverses 36 and contributes no escape. Unequal classes, or equal-weight
classes distinguished by degree, cannot be exchanged. The internal
independent-class permutations are even, so the four fusion examples are
globally even. Their two triple rotations may be chosen independently:
the same six anomalies apply to both predicates.

## 4. Independent computation: larger universe and direct permutation tests

`independent_residual_audit.py` imports no producer code. Rather than
assuming the four selected weighted quotients, it enumerates **every**
weighted presentation `(a,b,1,...,1)` with `5≤a≤b` and total degree
11, 12, or 13, and every quotient edge setting. The residual vertices are
unrestricted, including edges among them. Each presentation is expanded to
an ordinary adjacency matrix. Maximal independent-twin classes are then
recovered by identical full adjacency rows and canonicalized under all
quotient permutations.

The computation checks both subgroup predicates by actual vertex
permutations:

- enumerate pairs of disjoint triples outside a representative five-support;
- test both relative orientations of their synchronized product;
- test each triple rotation separately for the independent predicate.

It tests quotient automorphism lifts on the expanded adjacency matrix and
calculates sign directly as the parity of reversed **individual graph
edges**, not by the producer's weighted sign formula. All omitted internal
symmetric-group generators are directly checked to preserve adjacency and
have sign +1.

| Degree | Raw presentations | Distinct full graph types | Fusion | Escape |
|---:|---:|---:|---:|---:|
| 11 | 10 | 6 | 0 | 0 |
| 12 | 76 | 28 | 0 | 0 |
| 13 | 1,108 | 148 | 4 | 2 |

Both predicates agree on this candidate universe. All six anomaly types
agree exactly with the proof and producer audit. The independent run has
3,341 explicit runtime checks, all retained under optimized Python.
Normal and `-O` output files are byte-identical.

Hostile and boundary controls include:

- K_(5,5)⊔I2 admits a single triple off K but no double-triple extension.
- K_(5,5)⊔I3 gives the first corresponding escape.
- K_(5,5)⊔K3 admits the extensions but has no relatively even eligible
  marking, demonstrating why the residual true-twin triangle is excluded.
- K_(5,6) and K_(5,7) each have two large independent-twin classes but only
  one eligible class.
- K_(6,6) has two eligible classes joined by an even automorphism.
- Nonempty C5 circulants have the claimed odd normalizer elements.
- Directed-triangle blowups of regular cyclic modules of sizes
  3, 5, 7, and 9 attain the regular-module bound.

This computation is a complete check of the anomaly candidate universe,
not all graphs or all degree-thirteen subgroups. Its completeness for the
anomaly question rests on the independently checked two-large-twin-class
necessity above. It does not rederive the general Walsh matching theorem
or independently compute the absolute tournament totals.

## 5. Regular-module proof and independent source-side strengthening

Let a regular d-vertex tournament have a proper nontrivial module A of
size m. Its internal outdegree is constant, hence `(m−1)/2`, and the
numbers of outside vertices beating and losing to A both equal
`h=(d−m)/2`. For the set L beating A, the sum of its outdegrees outside A
is `h((d−1)/2−m)`. This is at least its internal arc count `h(h−1)/2`,
giving `d≥3m`. Equality says no arc goes from L to the opposite outside
set U, yielding `L→A→U→L`, with all three induced tournaments regular.
No equality condition or strictness issue is missing.

The original source Sylow proof is valid. There is additionally a shorter
independent proof that avoids all inherited small-degree image
classification and strengthens the source statement to degree fourteen:

**Corollary.** If an odd subgroup G≤S_n, `11≤n≤14`, contains a single
5-cycle and an element of cycle type `3+3+1^(n−6)`, then
`v_5(|G|)=1`. In particular this holds for either specified H.

**Proof.** Suppose 25 divides |G|. Since `v_5(n!)=2`, a Sylow-5 subgroup
of G is an ambient Sylow subgroup and contains two disjoint single
5-cycles, with supports A and B. A lies inside one G-orbit Ω. The induced
odd action on Ω admits an invariant tournament; transitivity makes it
regular. The orbit size is odd and at most thirteen. If A were a proper
subset of Ω, the single-cycle module bound would require `|Ω|≥15`.
Therefore A itself is a five-vertex G-orbit; the same applies to B.

An element of order three preserves each of these G-orbits. A nontrivial
restriction to either five-orbit would be a single 3-cycle on a regular
five-tournament, contradicting `5≥3·3`. Hence it fixes both A and B
pointwise. Only `n−10≤4` vertices remain, insufficient for two disjoint
three-cycles. Contradiction. The given single five-cycle ensures the
5-valuation is exactly one. QED.

For a tournament automorphism group one can use the induced tournament
on each orbit directly. Sylow conjugacy now gives the required one
eligible-K-copy orbit per full source type. This strengthening is a
review addendum, **not a modification of the frozen proof**, and does not
assert graph equality at fourteen.

## 6. Replays, receipts, and reproducibility

Both producer scripts were replayed normally and under `python3 -O`.
All four outputs match their corresponding frozen outputs byte for byte.
The cube replay confirms:

- synchronized: `R=1024, tau_T=768, tau_E=762, F_T^H=256, F_E^H=4, Z^H=258`;
- independent: `R=256, tau_T=256, tau_E=250, F_T^H=0, F_E^H=4, Z^H=2`.

These producer replays share the pinned cube/isomorphism/sign engine and
are reproducibility checks, not a second independent algorithm for those
absolute totals.

From this review directory:

```sh
python3 independent_residual_audit.py > /tmp/residual-review.json
cmp /tmp/residual-review.json independent_residual_audit.json
python3 -O independent_residual_audit.py > /tmp/residual-review-O.json
cmp /tmp/residual-review-O.json independent_residual_audit.optimized.json
sha256sum -c SHA256SUMS.txt
```

`review_receipts.json` records hashes of every frozen producer file and
each replay/new audit artifact. Particularly:

- frozen descent proof:
  `018fe2d3fec520e6e7cf71e1a9b57aa6e25b3064e78c2bac69c44bcf1175b674`;
- frozen degree-thirteen proof:
  `a6b78702b976254e90aae0ef68572cb37c3a3747932ce7b6852e7af4b51c1353`;
- frozen regular-module proof:
  `cf8354eb674efe169c80cbd6888a6e4f840a596c7d75401b8327a8451b18cb48`;
- inherited relative theorem:
  `c33157c00a065b2e470d59e67fcbe21e91084fb799343943186cd398fe136684`;
- pinned extension engine:
  `73a502ec5fc0bc299ac07e99706da9ea41a6882ae692431ee267fb627de36bc0`;
- pinned fixed-cube engine:
  `46cdad3ae5c4008dc912e83122383bfd6abd7f3a761688bc934e4705412d1684`.

No correction is required for the audited claims. The source-side
simplification is optional. The mathematical next boundary remains the
higher-degree eligible-copy geometry or a genuinely complete subsequent
subgroup census, not an extrapolation of these selected counts.
