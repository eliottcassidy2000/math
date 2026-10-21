# Independent review: balanced deletion groups and cut ports

2026-10-21. **PASS at the exact hashes in `AUDIT_PINS.json`.**

No unresolved mathematical or execution blocker remains. The reviewed main
note and separate all-depth addendum are sound within their stated labelled
and family-specific scopes. No repository file or Git state was modified.

## Proof review

- The inherited two-return decoder really makes deletion of the two
  centre-flanking vertices a Hamiltonian-path surjection. The centre itself
  is retained. Infinity means the full exceptional shore, not an ordinary
  finite subset; the endpoint descriptions are correct.
- The independent left/right new-bit toggles give exactly the fibres as
  V4 orbits. The two size-two stabilizers are distinct coordinate lines.
  Their multiplicities and the free/fixed strata in (3) and (4) are correct.
  A single globally identified stabilizer line cannot represent both types
  equivariantly. This does not exclude a more carefully marked arithmetic
  bridge, or a bridge on a restricted set of fibres.
- Rigidity and unique anti-reflection are inherited correctly. Reflection
  swaps the axes, giving a faithful dihedral group of order eight. This is
  different from the inherited simultaneous-complement/swap Klein four.
  The blue diagonal, its non-invariance under individual toggles, and the
  black/blue deletion direction are correct.
- The defect identity has the right sign: retained old path edges equal
  E(p) minus the bypass edges, hence E(p) symmetric-difference C. Converse
  commutes with deletion, so K has the same correction. Equality without
  correction is a statement about labelled tournaments, not isomorphism
  classes. The block-count proof gives the general r bound.
- The new two-seam reconstruction is valid at all finite and infinity
  boundaries, including START and END. Both infinity bits are physically
  forced to 1. Occupancy is pre = (epsilon=1 or eta=0) and post =
  (epsilon=0 or eta=1). The last scope sentence correctly says reconstruction
  does not automatically give a legal one-move labelled parent neighbour.
- The displayed T3 to T2 witness has the stated two literal paths, common
  parent, cut sets, black status and score spectra. It is at the smallest
  admissible parent parameter k=2. The two target score multisets prove
  nonisomorphism and nonisomorphism even after converse-merging.
- The all-order census is proved, not inferred from finite values. The
  identity |C|=2-e-a and all endpoint/adjacency overlap counts are correct.
  They give N0=2q, N1=2q^2-4 and N2=2q^2-6q+5.
- The all-depth addendum is also correct. Consecutive balanced insertions
  still fill only the same two seams, so the correction bound stays two,
  independently of d. The six disjoint boundary rows prove all three A1
  counts with D=2^d, including their reduction to the one-step case.

## Independent controls

`independent_literal_review.py` imports no producer or repository code.
It builds tournaments as Boolean adjacency matrices and directly enumerates
all directed Hamiltonian words by bitset DFS. Its group operations move
literal inserted vertices between runs. Its primary cut observable is the
XOR difference between the two actual target tournaments, checked against
retained-word gaps and then the two candidate seams.

- All 5,709 source Hamiltonian paths at orders 7,9,11,13,15 are obtained
  by literal DFS, not a coordinate generator. All orbit, stabilizer, V4,
  reflection, blue and eight-permutation checks pass.
- Their defect counts, for parent k=2,3,4,5,6 respectively, are
  (6,14,5), (10,46,25), (18,158,113), (34,574,481), (66,2174,1985).
- All 7,457 path/deletion cases for k=2..4 and d=1..3 pass the all-depth
  literal target identity, A1 counts and all nine oriented boundary strata.
- The general defect identity passes 174,720 Hamiltonian-path/deletion
  cases across all 1,024 labelled tournaments of order five, retaining
  only deletions that leave a directed Hamiltonian path.
- Normal and optimized independent outputs are byte-identical. Frozen
  copies of both producer checkers also replay normally and with -O with
  byte-identical outputs. The main producer retains 91,344 composition
  checks. All independent validation uses exceptions, not assertions.

## Labelled versus unlabelled neighbour countercontrol

For the displayed P00 child hostile, its deleted target is not exactly the
labelled R or K target of any of the nine parent T2 Hamiltonian paths.
However, the stronger claim after isomorphism is FALSE: parent path
(3,4,2,0,1) has a transitive R target with score multiset (0,1,2,3,4), and
its K target is transitive as well. Therefore the displayed target is an
unlabelled/converse-merged neighbour. The retained independent checker
asserts both facts. This does not weaken the actual obstruction: the two
child paths in one fibre have distinct deleted target classes.

## Changes resolved during review

1. Added explicit producer assertions for the defect census and the
   standard singleton-path correction.
2. Replaced ambiguous central-deletion wording by centre-flanking pair
   deletion, consistently retaining the centre.
3. Narrowed “repair of the incidence map” to deleted-target reconstruction
   and explicitly retained the labelled-neighbour scope.

The relevant inherited decoder, earlier 4/2/1 deletion and scalar laws,
arbitrary interior thinning, endpoint-growth results, and the two cited
label/permutation corrections were checked against the cited local sources.
No theorem-ID claim, unrestricted tournament classification, global
arithmetic identification or literature-priority claim is certified here.

## Reproduction

Run from this directory:

    python3 independent_literal_review.py > independent.normal.json
    python3 -O independent_literal_review.py > independent.optimized.json
    cmp independent.normal.json independent.optimized.json
    python3 frozen/check_balanced_group_ports.py > producer.normal.json
    python3 -O frozen/check_balanced_group_ports.py > producer.optimized.json
    cmp producer.normal.json producer.optimized.json
    python3 frozen/check_all_depth_ports.py > depth_producer.normal.json
    python3 -O frozen/check_all_depth_ports.py > depth_producer.optimized.json
    cmp depth_producer.normal.json depth_producer.optimized.json

The precise source and review-artifact hashes are retained in
`AUDIT_PINS.json` and `SHA256SUMS`.
