# Independent audit: transport-seeded generator and native metagraph

Session105, 2026-10-29. **PASS within stated finite scope; no correctness defect found.** No publication performed by this audit.

## Independent oracle and results

`check_independent.py` is an adjacency-matrix/permutation oracle, not a copy of either subject's verification helpers. Subject imports are used to obtain outputs and probe guards. Native canonicalization minimizes over all permutations with orbit memoization; generator tests enumerate complete isomorphism sets directly. The two subjects use opposite bit conventions, explicitly converted here.

Normal and `python3 -O` runs agree in all mathematical fields, with **1,649 explicit checks** each. `receipt.json` and `receipt.optimized.json` record the results (only the optimization-mode flag differs).

* Every raw labelled parent of orders 0..4: 1,1,2,8,64 matrices. Direct permutation isomorphism enumeration confirms every attaining map, all transported automorphisms, the exact normalized cache entry, and all binary-mask orbit representatives/literal children. Pointed bank sizes are 1,2,4,12,48.
* Generated ordinary child banks through order 5 match the complete labelled universe without duplicates: 1,1,2,4,12 classes. For **every pointed state**, independently recomputed complete child automorphism groups and every vertex-indexed rooted transporter set match. This catches retaining only one attaining map or an incomplete group even where class counts happen to match.
* Native M and blue B matrices independently match every entry through order 5. Independent Hamilton-permutation enumeration verifies W=a M, H=a h, symmetry, automorphism-mass detailed balance, and exact last-vertex strip recurrence.
* Independent n4→5 words 4/18 have equal ordered source/native-target class keys but unequal extension histograms. Words 2/6 retain the obstruction after converse-merging each endpoint and forgetting endpoint order. No ordered pair obstruction occurs at n2 or n3.
* Production canonicalization bound, partial-row guard, and native `require(False,...)` reject under normal and optimized Python.

The inherited 81/27 hostile is **replayed, not newly independently derived**: all 4,096 named subset-H entries agree, literal point stabilizers agree, yet attachment orbit counts are 176 versus 512 and transporter 0→4 sizes are 27 versus 0. Inherited partial/absent-port rejection and misleading algebraic Q identities are also replayed. These results are credited to the existing Session87/88 machinery.

## Full subject replays

`generator_replay.json` and its optimized counterpart execute the subject through child order 7. All six output hashes agree between old/new preparations at every level and between Python modes: states, Q, all root transports, automorphism groups, parent alignment data, ordered literal bank. The inherited independent audit exhausts labelled tournaments through six. Reproduced request/search pairs:

|child order|old requests|new requests|old searches|new searches|
|---:|---:|---:|---:|---:|
|1|3|1|1|1|
|2|5|2|1|1|
|3|9|4|2|1|
|4|29|11|4|3|
|5|110|42|11|8|
|6|652|256|56|46|
|7|6289|2649|576|525|

Native replay files agree byte-for-byte between normal and optimized Python: **23,038 checks**, full matrices through n6, and exact strip recurrence through n5→6. This final source includes the proposal-clock checks. Earlier audit replays had 19,586 checks before §7; the underlying native M/B matrices are unchanged. No discrepancy persists on the same source version.

## Proof review

The all-order torsor identity is correct with maps C→P: fixing a0, the maps a0^-1∘a give exactly Aut(C), bijectively. Canonical degree cells ensure sorted automorphisms have the inherited canonicalizer's attaining order. Increasing-mask full-orbit traversal produces exactly the old keys. Complete attaining sets are essential, and the cache remains keyed by literal matrices. This is an exact preparation optimization of unchanged Session88 adaptive full-port output, not a new canonical augmentation algorithm.

The native proof correctly separates symmetric path-orbit M, labelled-witness W=a M, and Markov P=M/h. Orbit-size detailed balance for W and stationary h-mass for P are distinct, correctly stated. The compression J*CJ and its [-1,1] spectral bound are valid; repeated P uses fresh path sampling, not repeated application of one involution. The SCC specialization preserves the actual child words and bridge ports. Strip intertwining is exact arcwise. The obstruction rejects strong lumpability of the stated endpoint-pair partition, **not every possible recurrence or larger sufficient state**.

This is the Hamilton-path native graph, not the separate up/down attachment operator, switching graph, or fundamental paired-black graph. No connectivity or scalar A000568 recurrence follows from these finite checks.

## Benchmark scope

The timed `build` covers full preparation, subgroup closure, Q, transports and groups; full output equality runs after the timed interval. Requests eliminated here are mostly cache hits, not searches. The documented ~1.05x illustrative timing must not be replaced by the 6289/2649 request ratio or described as an all-order runtime gain. Adaptive classified-root work is unchanged. The inherited Session87 smaller ordinary-class-only deliverable is explicitly not the timing comparator. The documented scope and caveats are honest.

## Reproduce

```sh
python3 -B audit_native/check_independent.py > audit_native/receipt.json
python3 -B -O audit_native/check_independent.py > audit_native/receipt.optimized.json
python3 -B generator/check_preparation.py --max-order 7 --repeats 1 > audit_native/generator_replay.json
python3 -B -O generator/check_preparation.py --max-order 7 --repeats 1 > audit_native/generator_replay.optimized.json
python3 -B metagraph/check_weighted_transport.py > audit_native/native_replay.json
python3 -B -O metagraph/check_weighted_transport.py > audit_native/native_replay.optimized.json
```

Run from the packet directory. The exact pinned Session87/88 source bytes are bundled under `vendor/`; the generator prefers that root, and the independent script passes it explicitly. No external checkout is required. Finite oracle claims stop at the bounds above; all-order statements rest on the displayed proofs, not extrapolation from these runs.

## Final §7 clock bridge audit

The fixed-permutation guarded move is an involution on the entire labelled universe (identity on invalid guards), so uniform labelled stationarity is sound. Reversal pairs bound H by n!/2, making R half-lazy. Reversible half-lazy R is positive semidefinite, as is the refresh kernel. The shared mu proportional to 1/Aut makes the Dirichlet-form mixture argument valid; it proves both absolute and right gap bounds 3 alpha/n for n>=3. This does not use or prove native-only connectedness. The conditioned actual-path kernel instead has stationary mass H/Aut and must not be silently substituted into the mixture.

The up/down proof defines its class transition matrix column-stochastically. In the native row-kernel convention its transpose is used, leaving spectrum and reversibility unchanged. This orientation clarification is now explicit in §7 rather than treating unlike matrix conventions as identical.

Independent added checks construct refresh entries by directly resampling each labelled vertex's incident edges, without importing the up/down checker. Through n5, exact rational checks verify both kernels' row sums, detailed balance, weighted positive semidefiniteness, proposal half-laziness, the n3 matrix, and stationary acceptance exactly 2^{−(n−1)}. At alpha=1/4,1/2,1, rational Schur elimination proves the claimed mixed gap bound on the entire small class space via positive semidefiniteness of the weighted matrix (1−3 alpha/n)I−L+(3 alpha/n)Pi. These finite checks supplement, not replace, the all-order form proof.

The sparse proposal acceptance caveat is essential and correctly retained: the guaranteed gap comes from the refresh fraction; this is not a native-only efficiency result or a uniform sampler of isomorphism classes. Final subject/proof SHA256 values are recorded in both independent receipts.

## Portable final packet replay

The final generator loader and four exact inherited sources are self-contained under `vendor/`; `SOURCE_MANIFEST.json` records their provenance, hashes and sizes. The audit script locates its packet with `Path(__file__).resolve().parents[1]` and supplies that packet's vendor root explicitly. It has no workspace-name or sister-checkout dependency.

After the final loader update, the entire packet was copied beneath a fresh temporary directory and both independent Python modes were run with working directory `/tmp`. All 1,649 checks passed and mathematical receipts agreed. The standalone generator was separately run without a baseline override from that isolated copy through child order 7. All six structural output hashes at every order exactly match the historical pre-vendoring full replay. The older timing replay files remain historical measurements, not timings of the final isolated run. `portability.json` records this verification; independent receipts now pin the final generator, proofs and vendor manifest.
