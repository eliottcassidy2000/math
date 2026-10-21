# Incoming cycle-count polynomial module: bounded independent review

Status: mathematical normalization, integer implementations, signed typing, and scoped comparison PASS. Incoming work is attributed to the `tournament_tail_kernels_20261021b` packet at commit `21ea39dea84ee7f229d1f92bcbf8387e13e29cd5`, not claimed as this session's invention. The packet compares against our degree-local scalar source from `622e4451384c6d6ec6ee56bba0f71dc50732512d`.

## Mathematical audit

For actual odd cycles a_i, set t to their number, e=sum_i(a_i-1)/2+sum_{i<j}gcd(a_i,a_j), g=e-binomial(t,2), and S_j=sum_i(gcd(a_i,j)-1). Subtracting the Gaussian increment kt+binomial(k,2) from the old insertion exponent gives exactly (ell-1)k²/2+k S_ell. This is independent of the retained polynomial exponent t, so aggregation by degree and residual future-response key is legitimate. The polynomial coefficient of y^t retains the information required to recover the original sum. The final linear functional G(y^t)=2^binomial(t,2) is not multiplicative (G(y²)=2 but G(y)²=1).

The raw normalized weight 2^g/z has g>=0 and z dividing M!, where M is actual mass. Thus N! and committed-degree d! both clear denominators when M<=d<=N. The global coefficient ladder divides by ell*k exactly on every raw-history contribution and hence on merged signed sums. The local multiplier (d+ell*k)!/(d! ell^k k!) is an integer by labelled cycle selection. Excess then multiplies by the integer factorial ratio d'!/(d+ell*k)! and signed binomial coefficient. Final divisions use N! and d! respectively; the unnormalized minimal odd denominator cannot simply be recycled after removing the Gaussian factor.

Core excess s changes degree by ell*s and carries (-1)^s binomial(k+s-1,s), but does not create actual cycles. Only k changes t and future responses. At the terminal length, G(y^k P) is common to all excess choices and may be calculated once before distributing signed excess. This is a valid linearity optimization, not an interpretation of signed intermediate aggregates as nonnegative class counts. No-fixed coefficients remain rational masses (degree 3 is 2/3).

Hybrid conversion is valid only with fixed points, since R_1 supplies t. A full-degree state is settled directly. Affordable future lengths only decrease, and residual responses fit the existing bounded digit representation; polynomial exponents remain essential sidecars. The result is an exact all-finite-order change of state, not an all-order polynomial runtime bound, classification, representative construction, automorphism recovery, or metagraph reconstruction.

## Independent controls

The checker requires an explicit repository root and does not infer a checkout location, use Git discovery, or fall back to legacy paths. From this workspace, run:

    python -B session22-incoming-enumeration-review/check.py --repo-root math
    python -B -O session22-incoming-enumeration-review/check.py --repo-root math

Both return identical 25,763-check PASS receipts. No assertion is load-bearing and the checker does not invoke the producer's controls. Tests cover direct unmerged rational partitions and explicit signed substitution in all four modes through 24, global and local polynomial kernels, the frozen degree-local scalar, full-payload versus terminal folding, four hybrid switch lengths, normalized raw-history arithmetic, local transfer divisibility, nonmultiplicativity, rational observer typing, benchmark source pins, median calculation and all timed coefficient hashes against the saved full vectors.

The incoming README reports an independent parent proof/API audit PASS. There is no separate independent audit artifact in this five-file packet; this review does not treat that statement as our own execution evidence. The existing degree-local independent audit in `session22-enumeration-review/README.md` was also read. It establishes the scalar normalization and cautions about unequal telemetry in the older comparison.

The declared upstream verifier was also rerun after its inherited source dependency was materialized from the sparse checkout. Normal and optimized execution each pass all 79,565 checks and reproduce the frozen receipt exactly; see the identical `upstream-normal.json` and `upstream-optimized.json` receipts.

## Current engine selection and benchmark comparability

Stored incoming measurements are three alternating-order end-to-end samples on a shared Windows CPython 3.13.14 host, including setup and output conversion. They compare identical complete coefficient vectors. They are not isolated attribution to a single algebraic identity, and must not be pooled with this session's differently hosted earlier benchmarks.

| Engine | Ordinary through 200 median seconds | Core through 100 median seconds |
|---|---:|---:|
| Previous common-odd scalar | 4.645314000001235 | 1.0660102999991068 |
| Our degree-local scalar | 2.1967132999998285 | 0.8516732999996748 |
| Incoming global polynomial | 3.4817127999995137 | 0.24856429999999818 |
| Incoming local polynomial | 3.2261822000000393 | 0.3150282000005973 |
| Scalar prefix, polynomial unit tail | 2.1777284999989206 | 0.3222315000002709 |

Practical scoped selection: retain our degree-local scalar as the ordinary baseline; hybrid's approximately 0.86% lower median is effectively tied in these three shared-host samples, not a demonstrated ordinary improvement. For signed core through 100, select incoming global polynomial: approximately 3.426 times faster than our scalar in these samples. There is no universal fastest-engine conclusion.

Ordinary through 200 structural transitions: 2,994,814 to 1,515,929; peak profile keys: 96,836 to 63,299. Polynomial peak coefficient cells are 96,842, and global polynomial performs 2,355,764 coefficient-ladder divisions. The structural count reduction is not a corresponding reduction in total scalar arithmetic. Core through 100 structural transitions: 1,697,975 to 289,706; polynomial peak coefficient cells 5,493 versus 5,480 scalar keys.

Separate traced-memory observations, not timing runs: ordinary through 100 peaks are 1,449,988 bytes for our scalar and 2,109,552 for incoming global polynomial. Core through 60 peaks are 4,428,468 and 315,624 bytes respectively; previous common-odd scalar is smaller still at 278,020 bytes. Thus neither fewer keys nor the observed core speed advantage establishes universal memory dominance.

## Source pins and scope

All hashes are LF-normalized SHA-256 and were verified locally against benchmark metadata:

- Incoming check.py: 6c72647b829565fa97429c188ca0f58da5bd35fc4d617efb4977c0c69c8a2272
- Frozen our incoming_local.py: b76167ef82dda0b134681abce33e159d5dd80d15e92edb9e31ab9b490f6939ad
- Incoming receipt.json: d45245f0b1ae5321beb3152c37ce8e5246a355eb65d2349cf31031df9c11c5d2
- Incoming benchmarks.json: 2bb79e980b5f11b825db6ab3c591cbc1b45cbd0163f46ab8eb2cc2b9784f1ca9
- Previous scalar check.py: 8ef8429ff53764dc54af36ac45eb1e2729b7cf8e67c0b3d1309b8e9a308fcdbd

No expensive benchmarks were repeated. No source optimization, Git mutation, or remote write was undertaken by this reviewer. This is a bounded audit and integration recommendation, not a new optimization branch.

## Final portable checker pin

Final checker SHA-256: `2881b920373247bdc87063e10c7700d6bde26b38b4697aecf671f384241006ee`. Both commands above were rerun after adding the required `--repo-root`; the resulting 25,763-check normal and optimized receipts are identical. In another checkout, replace `math` with its repository root, and invoke this checker by its actual path. Incoming source and benchmark hashes are unchanged.
