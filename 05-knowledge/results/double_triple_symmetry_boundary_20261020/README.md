# Double-triple symmetry: descent through twelve and exact failure at thirteen

20 October 2026. **PROVED structural results; FINITE-EXACT selected counts.**

For the two specified actions

- H_d=C5×C3 on 5+3+3+1^(n−11), with the C3 generator synchronized across the two triples;
- H_i=C5×C3² on the same orbit partition, with independent triple cycles,

full tournament types and globally sign-even graph types agree at n=11,12, and differ by exactly six at n=13.

| embedding | n=11: tau_T=tau_E | n=12: tau_T=tau_E | n=13: tau_T | n=13: tau_E |
|---|---:|---:|---:|---:|
| H_d, synchronized | 12 | 96 | 768 | 762 |
| H_i, independent | 4 | 32 | 256 | 250 |

The degree-eleven/twelve numerical rows are inherited from the [complete bounded census](../mixed_prime_census_extension_20261020/README.md). Degree thirteen concerns **only these two explicitly embedded subgroups**. It is not a complete subgroup census.

## Mechanism and proofs

1. [Predicate-restricted descent at n=11,12](DOUBLE_TRIPLE_DESCENT.md). Localize the relative Walsh matching at a single C5, retaining the demand for two disjoint triple cycles. Count only H-extendible copies of C5. At these degrees eligible copies fuse uniquely and no odd graph automorphism escapes the C5 normalizer.
2. [Sharp degree-thirteen boundary](SHARP_N13_BOUNDARY.md). Exactly four globally even graph types have two eligible-copy classes, and exactly two relatively even graph types have an escaping odd automorphism. The source has no C5-copy fusion, giving the exact defect 4+2=6 for each embedding. The six anomaly types are proved exhaustive within this mechanism, not discovered by assuming a subgroup list at thirteen.
3. [Regular-module bound](THREE_CYCLE_SUPPORT_LEMMA.md). Every proper nontrivial homogeneous module A of a regular d-vertex tournament has 3|A|≤d. Equality forces three regular modules arranged as a directed triangle. A single-cycle automorphism support is such a module. This elementary geometric statement excludes the unwanted transitive source case at thirteen, without a degree-thirteen group classification.

Graph sign throughout means the edge-reversal character, not ordinary vertex-permutation parity. The synchronized group has order fifteen; the independent group has order forty-five but is not cyclic C45. Their embeddings and their H-based versus C5-based fusion defects must not be conflated.

## Reproduce the controls

From this directory, in the published sibling layout:

```sh
python3 weighted_copy_audit.py > /tmp/weighted-copy-audit.json
cmp /tmp/weighted-copy-audit.json weighted_copy_audit.json
python3 -O weighted_copy_audit.py > /tmp/weighted-copy-audit.optimized.json
cmp /tmp/weighted-copy-audit.optimized.json weighted_copy_audit.optimized.json
python3 double_triple_n13_probe.py > /tmp/double-triple-n13.json
cmp /tmp/double-triple-n13.json double_triple_n13_probe.json
python3 -O double_triple_n13_probe.py > /tmp/double-triple-n13.optimized.json
cmp /tmp/double-triple-n13.optimized.json double_triple_n13_probe.optimized.json
sha256sum -c SHA256SUMS.txt
```

The weighted-quotient audit is independent of the cube/sign/isomorphism engine. The selected full-cube calculation uses the pinned engine from the predecessor census extension. See [VERIFICATION.md](VERIFICATION.md) for precise review and independence limits.

## Scope and next target

The bounded descent and its degree-thirteen exact defect are proved. A full degree-thirteen/fourteen subgroup census, a classification of defects in every higher degree, and a general object-relative descent theorem remain open here. A useful next target is to track the eligible independent-twin-class geometry as further fixed vertices are added, preserving the extension predicate rather than reverting to unrestricted C5 copies.

The [inherited relative fixed-cube theorem](../pronormal_quotient_bijections_20261020/PRONORMAL_QUOTIENT_THEOREM.md), Theorem 1, supplies the matching. The [degree-at-most-ten census](../odd_mixed_prime_census_20261020/README.md) and its [degree-eleven/twelve extension](../mixed_prime_census_extension_20261020/README.md) supply the bounded orbit-image classifications. The regular-module argument is included in full; no literature-priority claim is made.
