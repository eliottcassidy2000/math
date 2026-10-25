# Independent range-signature algebra audit

Status: PASS for mathematics and finite computations. Final producer assertion-mode guard verified; independent normal/-O replay passes.

Audited packet: /tmp/session56-overlap-signatures. No producer edits, deletions, remote writes, or changes to the research checkout were made.

## Mathematical verdict

The family idempotence iff meet-closure proof is correct, including the empty family. The compression lemma correctly uses finite minimal elements, and the antichain idempotence iff empty/singleton proof does not incorrectly assume the first new intersection is itself minimal. The nonnegative-integer histogram classification follows from total mass and is correctly scoped away from signed/probability coefficients. The distinction between independent equal-characteristic factors, extension fields, and Z/25Z is correct.

The overlapping query explanation correctly preserves order: P*Q*Q*R need not equal P*Q*R. The disjoint sparse-table highest differing bit chooses a single block midpoint and exact, ordered, nonoverlapping partition. The segment-tree two-sided accumulators preserve order. Both correctness claims require well-formed, nonempty input arrays and valid interval indices; empty input/invalid queries are not audited API guarantees. The duplicate-component guard is sound for well-formed leaves and rejects the concrete overlap. It is not a full validation API for arbitrarily forged Signature objects. Pairwise coprimality is checked for each test sequence, not inferred from unique IDs.

The complexity claims distinguish aggregate-operation counts from mask-word, integer coefficient, eager provenance, output, and local enumeration costs. O(n^2 K) eager disjoint-table witness references is a valid bound: at a scale of length h, total interval lengths are O(nh), and summing scales gives O(n^2). No bit-complexity or constant-size signature claim is being smuggled into RMQ asymptotics.

## Independent finite checks

check.py uses explicit exceptions throughout, so its verification remains live under Python optimization. It does not call the producer main verifier. It independently counts all arrays from top-row/center line equations, counts full multiplicities, independently performs incremental CRT, and checks local residues, square membership, all eight line sums, and exact 36-bit collision masks.

Normal and optimized runs both pass:
- 31 interval histograms across 20 distinct moduli, including composite prime powers.
- 134 global CRT witnesses, all emitted with entries, exact collision mask, and square roots in global-witnesses.json.
- 88,560 noncommutative tuple-concatenation interval queries at lengths 1 through 80.
- All 256 three-bit families, 65,536 compression pairs, and 256 nonnegative-integer histogram probes.
- Duplicate-leaf and actual overlapping-query provenance rejection.
- Z/30Z true count 112 versus duplicated-F5 count 784; mask zero appears only in the latter. True minimal masks are 0x1000880 and 0x80040200. Z/25Z remains Parker.

The independent enumeration follows the same line-equation strategy as the inherited independent enumerator but is separately written, counts multiplicities rather than only supports, and does not call either inherited enumeration. Data-structure code is imported from the producer and tested against independently constructed outputs, rather than copied into an alternative implementation.

## Source review

Primary paper: Michael A. Bender and Martin Farach-Colton, The LCA Problem Revisited, LATIN 2000, LNCS 1776, pp.88-94. Independently read section 3 at printed p.91 (PDF page 4) via https://cadmo.ethz.ch/education/lectures/HS18/SAADS/papers/lca.pdf on 2026-10-25. It supports overlapping power-of-two RMQ blocks and O(n log n) preprocessing/O(1) query in its model. No claim of novelty for disjoint sparse tables or segment trees, and no attribution of the family algebra to that paper. This audit paraphrases only; no third-party text or PDF is republished.

The frozen Session55 verifier hash matches SOURCE_PINS.txt exactly. The inherited theorem original was subsequently resolved at /workspace/scratch/07d5532f2fce/session55-crt-collision-wiring/THEOREM.md, and its SHA-256 72945f751dd6a1f79b899b202233bbf27dc78fefd0a3420f217979334b956450 exactly matches the recorded source pin.

## Reproduction and limitations

python3 /tmp/session56-audit-range-signatures/check.py
python3 -O /tmp/session56-audit-range-signatures/check.py

An optional first positional argument supplies the producer packet directory when relocated.

The producer README explicitly supports only unoptimized execution. Its original verifier nevertheless prints PASS with assertions disabled, including removal of certificate equality checks. The producer subsequently added an explicit rejection in main(). Its optimized CLI now fails with RuntimeError as intended; normal --verify succeeds. Import under optimization remains supported for the independent checker. No claim is made that the producer itself passes optimized verification. Final producer verify.py SHA-256: fff525318dfc1e181ce38bb6908273416be0650975eaa5546a7fa50818a5b5fb.

Startup was read-only. Targeted processor inspection was blocked by absent .machine-id; no identity was invented. The research checkout remained clean at the final status inspection. Parent handles publication and remote checkpoints.
