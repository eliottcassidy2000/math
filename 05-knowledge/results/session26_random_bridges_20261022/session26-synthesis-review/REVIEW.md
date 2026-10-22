# Final Session 26 synthesis and exploration audit

Verdict: PASS. Reviewed 2026-10-22. No remaining blocking correction identified. This is a source-bound factual and mathematical synthesis audit, not a novelty assessment or a replay of every exhaustive census. No Git or repository writes were performed.

## Final source binding

- SYNTHESIS.md: SHA256 `0ec0d0e03f402bfaea30c5d0135ff7e09ad838bbe3db87c9a35950e11d083802`.
- EXPLORATION.md: SHA256 `cd559cf6596eecbf3a21455a542fe0ab4727ad57f735626d19f96868e7a3ad55`.

These are the final files in the root handoff directory. Their original versions were read completely; the three requested precision edits were subsequently verified. All 24 Markdown links resolve to supplied source, review, or manifest files under the specified publication layout: these two files at packet root, lane folders as siblings below that root. This is a local target/layout check, not a claim that remote publication or user access has already been verified.

## Cross-checks and resolved wording

The feedback inequalities and fixed-k asymptotics agree with the fixed-feedback proof and independent audit. Bounds are correctly scoped to actual path-complement edges, both colors, and fixed k. The upper inequality is not advertised as sharp. The generic class and tiling formulas use different measures and preserve their multiplicities.

The one-apex class, self-converse, blue-incidence, and tiling formulas match the full blue note. Independent arithmetic reproduces all entries of the n=8,9,10 boundary table. The last-wiring independent audit confirms the exact six black edges, two path components, double edge, 47 isolates, eight touched nodes, no loops, and all distinct endpoint labels. The 12-incidence deduction at n=8 is correct. Induced bridges are not promoted to full-graph bridges.

The cut/Ferrers formulas, strong extrema, H=13 exclusion, product density, mean and variance, CDF rate, bounded-H main terms, reverse-complement bound, and missing-odd lower density match the classification/spectral proof and its review. Padding and strong-core restrictions remain explicit. The electrical mean is unconditional and is not confused with the strong-word mean. The binary-product approximation is normalized at the relevant scale; its source gives the explicit error (m/2+1)2^(-m).

The attachment transform, profile recovery, signed-channel control, counter state bound, and fixed-k arithmetic complexity match their proofs and independent audits. The observer, quadratic-reader, and runner sections preserve the important distinction between deterministic memory, linear rank, arithmetic reachability, and signed normalized moments. Runner-count recurrence obstruction is not generalized to nonlinear reconstruction.

Three precision repairs were requested and confirmed:
1. EXPLORATION now describes the four-dimensional rational Hankel-row span, rather than suggesting it is a much larger memory dimension.
2. EXPLORATION's exact lcm denominator-hostile family explicitly has integer Q>=14, matching its proved source domain.
3. SYNTHESIS specifies the finite abelian reachable group, binary affine update, invertible linear update, and homomorphic output hypotheses needed by the constrained-reader theorem.

The archival-repair summary matches the corrections packet, including the exact 1/8 survivor and the absence of a bounded rational census implication. The preservation statements concerning incoming S50/S51 work are publication bookkeeping supplied by the parent, not new independently replayed mathematics. No such replay is claimed in the synthesis.

## Independent electrical check

For each j<i and binary entries, 2 b_j(1-b_i)=(b_j-b_i)^2+b_j-b_i. Multiplication by w_(i-j) and summation gives precisely b^T L b+d dot b, with the displayed sign of d. Thus the new identity is correct for every binary word, including the empty word with zero matrices.

The field sum is zero. Complement leaves the quadratic energy unchanged and negates the field; reversal does the same by distance invariance and antisymmetry of d. Reverse-complement therefore preserves H. The two complement sum/difference identities follow immediately. Words 10 and 01 are a valid hostile to energy alone.

For m>=1 the total conductance is sum_(d=1)^(m-1)(m-d)w_d=2^m-2. Independent fair bits cross each pair with probability one half, while the field has mean zero. Hence E H=1+(2^m-2)/4=2^(m-2)+1/2, exactly the positive-index k=1 attachment mean. The displayed m>=1 restriction matters: the empty word has H=1, not 3/4.

The persisted independent check exhausts 1,022 words of lengths 1 through 9. For each it computes the cut, Laplacian energy and field, complement identities, and a separate adjacency/subset-DP Hamiltonian count. Every identity and unconditional average passes. No producer recurrence is imported.

## Provenance and controls

The four sampling manifests contain exactly 103 paths, all distinct, in sizes 24,24,24,31. All 24 root archive copies match the SHA256 values in ROOT_SAMPLE.json, in the same order. This verifies local copied-file integrity and manifest consistency; it does not independently refetch the historical commit or certify every assertion in the sampled archive.

Run check.py or python -O check.py from this review directory. Normal and optimized outputs match byte-for-byte. check.json records source hashes, sample/hash/link totals, direct Hamiltonian controls, and exact boundary rows. Guards are explicit exceptions, not optimization-sensitive assertions. The checker uses absolute source locations for this audit workspace; a relocated publication must adjust those input paths before rerunning.

Scope exclusions in both final documents are appropriate: no complete Hamiltonian spectrum, unrestricted tournament count asymptotic, all-metagraph wiring or bridge classification, LRC(14), Collatz, Hadwiger–Nelson, or native 153N ledger completion is asserted.

## Portable publication replay addendum

A separate `check_portable.py` now supports the published layout without altering the historical `check.py`, `check.json`, or `check-optimized.json`. Run `python3 session26-synthesis-review/check_portable.py` from any working directory when this review folder sits directly beneath packet root. Alternatively pass the packet root as its sole positional argument. Python `-O` is supported.

The portable checker verifies the exact final SYNTHESIS/EXPLORATION hashes, all 103 manifest paths and their uniqueness, the four lane sizes against SAMPLING_SUMMARY.json, the 24 root manifest path/hash records, all 24 links strictly relative to packet root, the same 1,022 electrical/direct-Hamiltonian controls, and the exact boundary arithmetic against the displayed synthesis table. It imports no producer code.

Archived sample copies intentionally do not ship with the packet. Consequently portable output explicitly reports archived-copy verification as NOT_RERUN and explains that omission. Manifest consistency is not substituted for fresh archive-byte or remote-commit verification. The historical output records the original local verification of all 24 available copied files; that historical claim remains unchanged and is not attributed to the portable replay.

The portable checker was run both normally and under `-O` against a temporary packet-equivalent layout outside the repository, with symlinks to the frozen root documents/manifests and lane directories, and no root archived sample copies. Both runs passed and produced byte-identical `portable.json` and `portable-optimized.json`. This establishes portability across the packaging change, not remote publication verification. No mathematics, source claims, original checker, or original outputs were changed by this follow-up.
