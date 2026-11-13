# Audit, scope, and correction lineage

**Verdict: PASS for the explicitly scoped statements and finite universes.**
Status promotion is bound by reversible before/after edits in PROMOTIONS.json.
No general graph algorithm, classical theorem or inherited census count is
claimed as a new discovery. No theorem ID was reserved.

## Complete census and independently reconstructed incidence

The closure producer constructs all relabeling orbits of every labelled
orientation through n7 (2^21 orientations at the top order) using a C++
adjacent-transposition orbit map. It checks every representative's unrestricted
canonical minimum by all vertex permutations. It constructs the full ordinary
and merged matrices in two ways: fixed-spine words, and physical Hamilton paths
of ordinary representatives divided by the free automorphism action. Both
ordinary and merged blue matrices agree as well. Ordered covers are independently
computed by permutation break counts and subset convolution. There are no
strongness, height or symmetry filters. Normal and optimized runs give the
same data, with 67,869 always-active gates.

The root [independent checker](audit/independent.py) imports no producer code.
It builds relabeling orbits directly by all vertex permutations, then enumerates
literal physical paths and covers. It matches every ordinary/merged colored
entry through n6, verifies full-Q refinement through n7 and independently
certifies weighted merged rigidity. It has 12,305 gates and matching normal/-O
receipts. This is a second complete census through n6, not through n7. The
pairing producer supplies an additional complete through-six census and six
independently reconstructed full physical-path rows at seven.

## Minimum measurements

[The closure review](audit_closure/REPORT.md) checks the finite partition
universal property, all four exact measurement palettes, every named-anchor
witness and packing certificate. Its brute-force subset search independently
establishes minima through n6. Its differential small-family oracle and graph
partition controls give 16,097 gates, identical under -O.

The root [second anchor solver](audit/anchors.py) independently reconstructs
all 184 equal-Q collision pairs at n7. It uses fresh residual antichains at
every search state, rather than the producer's global mask reduction. It rejects
budget eight for all four palettes and independently checks nine-column
witnesses. States visited: A 38,757; P 36,606; colored A 38,620; colored P
35,943. Its 493 gates include 480 brute-force small controls.
[Its independent source audit](audit_closure/ROOT_ANCHOR_SUPPLEMENT.md) validates
every pruning and branching rule, including 24,775 additional differential
controls. The full normal and optimized root receipts are byte-identical and bind the actual n7 exhaustion; the review
records the distinction between source review and execution evidence.

The coarsest stable partition is discrete through7. This justifies reducing
closure to separation for these orders; it is not assumed in the measurement
optimization. Costs are named columns, with two outputs for a colored column.
An arbitrary scalar coding the whole class bank is not admitted as compression.

## Pairings, transport, and root proofs

[The signed-lift review](closure/PAIRINGS_REVIEW.md) audits reconstruction,
source normalization, the support-component kernel, guarded base lifts,
representation sectors and cycle indices. Its separate checker has 1,260,321
always-active checks, with identical normal/-O receipts. These include full
operator intertwining, a separately implemented union-find and weighted graph
refinement, kernel actions and small literal Burnside/cycle checks. This tests
all n7 matrix entries from the frozen complete census, not a second generation
of every order-seven orbit.

[The bridge review](audit_bridges/REPORT.md) passes B1–B6 with 3,504 independent
gates, including small graph holonomy, nonabelian triangle controls, separate
finite-site plus calculations, tensor sections and S5/A5 perfect-quotient tests.
The bridge producer passes 5,432 gates. The primary Stacks, Bott, ABS and Quillen
sources were checked to the scopes actually cited. Neményi attribution retains
the prior packet's source limitations; an unread erratum is not imported.

[The root review](audit_bridges/ROOT_REVIEW.md) passes O1–O3 and C1–C3 after
two scope repairs, with 3,605 exact controls and identical normal/-O receipts.
The elementary algebra, derivative and Clifford claims have written general
proofs; their finite controls do not substitute for those proofs.

[The black-core root review](closure/BLACK_REVIEW.md) separately verifies
connected-fibre contraction, the self-converse attachment criterion, and the
order-five converse hostile. Its producer adds 69,050 normal/-O gates. The
status-only edit is bound in BLACK_PROMOTION.json; it leaves the primary
promotion manifest unchanged.

## Corrections before promotion

1. O1 initially required the starting space merely to contain the initial
   partition indicators. Extra initial observers can force a finer algebra.
   The repaired hypothesis makes that space exactly the initial block-constant
   algebra. Three points with identity dynamics and the extra observer (0,1,2)
   give the minimal hostile used by the independent review.
2. O3's consequence is expressly a lower bound for strongly lumpable partition
   quotients refining the specified observer. Selected scalar statistics can
   have a different linear realization.
3. The ordinary random-walk plus block is P_Y in the natural fibre-constant
   basis, and M^(1/2) P_Y M^(-1/2) in the orthonormal sum basis. The signed-lift
   draft omitted this basis distinction; the code and every count were unchanged.
4. A pre-audit bridge discussion was tightened to preserve **flat** transport:
   -1 holonomy can remain nontrivial while the underlying complex line, or
   sixteen-fold real bundle, is explicitly trivial. No KO/KU persistence claim
   is made. The final source already includes both explicit trivializations.

These mechanisms are recorded in the correction lineage, including the entry
in 01-canon/MISTAKES.md. No accepted theorem is retracted.

## Reproduce without changing the truth surface

From the repository root, the following producers/checkers have matching normal
and optimized receipts. Use Python with -B to avoid incidental bytecode files.
The closure producer compiles its local C++17 helper in ignored `_scratch`.

- `closure/check.py`: complete matrices; optional `--output data.optimized.json`.
- `closure/analyze.py`: exact anchor minima and Boolean sidecars.
- `pairings/check.py`: signed kernels, valid lifts and independent small census.
- `bridges/check_bridges.py`: finite holonomy/site/period controls.
- `audit/independent.py`: independent orbit, cover, matrix and closure checks.
- `audit/anchors.py`: independent n7 minimum lower bound.
- `closure/check_signed_review.py`: independent signed-operator consequences.
- `audit_closure/check.py`: finite solver and universal-property audit.
- `audit_closure/check_root_anchors.py`: source-audit solver differential controls.
- `audit_bridges/check_independent.py`: independent bridge controls.
- `audit_bridges/check_root.py`: independent observer and Clifford controls.

All paths above are relative to this packet. Some audit receipts bind the
reviewed prepromotion text; PROMOTIONS.json retains the exact reversible relation
and postpromotion rebindings are explicitly recorded where needed. The closed
manifest excludes generated `_scratch`, bytecode caches and its own digest.
SOURCE_PINS.json binds inherited local references to the Git baseline rather
than silently following future edits. No Lean result is claimed.

Frozen producer JSON retains its generation-time candidate status as provenance;
the proof promotion and this audit determine the current scoped status.
Repository check: `python3 agents/check_docs.py` passed the bounded startup,
frontmatter, correction-sentinel and router checks.

The final [joint-shadow review](audit/JOINT_SHADOW_REVIEW.md) passes the exact
all-k corollary combining the incoming 21-vertex core with full-Q ordinal
multiplicativity. The root replayed the core certificate normally and with -O.
SHADOW_PROMOTION.json records its status-only promotion; no core Hamilton-path
census or all-order native-incidence claim is made.

The [F8 action review](audit/F8_REVIEW.md) independently checks the algebraic
map and replays all 5,320 finite gates in both Python modes, including the
failed seven-point quotient. F8_PROMOTION.json retains the status-only edit.
