# Hybrid palette-compression audit addendum

**Correctness: PASS. Performance: workload- and path-planning-dependent; no uniform speedup certified.**

This addendum is separate from the earlier frozen audit. All files covered by the original `SHA256SUMS` were rechecked and remain unchanged. Producer files were read, not modified.

## Exact tensor indexing

The auditor independently checks metadata alignment, edge order, tensor-axis order in the einsum specification, per-axis palette lists, and every entry of the dense coefficient tensor. Re-expanding each compressed tensor recovers all **31,396 original dense entries exactly**, from **12,970 stored compressed entries**. Unsupported indices are exactly zero. This is an exact storage reduction of the full tensor boxes, distinct from the previously certified 553 nonzero coefficients and 263 quotient contraction monomials.

The axis order is preserved by `np.ix_(*palettes)`: the jth axis selects exactly the palette for the jth core edge. The corresponding matrix operand selects that same ordered palette. A loop uses its matrix diagonal, consistent with the single vertex index in the contraction specification.

Four cores have all four channels on every edge and take the fallback path. These include K4 and the doubled square, whose specialized original routines are retained. They also include four parallel edges and the three-vertex multigraph with multiplicities (2,2,1). The audit covers every branch.

## Independent execution checks

1. For all 31 cores, the auditor brute-force sums over all vertex labels using four completely arbitrary dense integer matrices of order 3. No skew-symmetry, even/odd parity, or rational relationship among matrices is imposed. Every result matches the hybrid contraction exactly. This is a direct axis/operand test, rather than merely a comparison on special tournament inputs.
2. Four additional fixtures are used: a weighted skew matrix of order 3, a transitive tournament of order 4, and almost-regular deleted cyclic tournaments of orders 6 and 10. At z=1 and z=0.37+0.19i, all **248 core comparisons** agree with the original five-basis engine. Maximum absolute difference: **2.220446049250313 × 10⁻¹⁵**.
3. Real and complex output agreement supports the basis implementation as well as the compressed contractions. The exact identities and independent table transformations remain those proved in the main audit.

## Timing anomaly and controlled comparison

An initial warm benchmark unexpectedly made the old four-basis engine tens to hundreds of times slower than the hybrid. Those ratios must not be advertised as a palette-compression speedup.

Profiling a fresh process showed that original `CoreB.value_b` uses NumPy's `optimize='greedy'` default intermediate-size limit and caches the resulting path without an n key. Fresh planning at n=10 already chose an all-at-once contraction for some cores. Thus this particular slowdown is not simply caused by reusing a smaller-n path. The hybrid uses `optimize=('greedy', 2**24)` and keys its main path cache by matrix order. Existing smaller-fixture plans can incidentally perform much better than fresh plans under the default cap; comparisons depend on planning history as well as representation.

To isolate this issue, `benchmark_equal_cap_addendum.py` replans both engines at each n with the same greedy intermediate cap, 2**24. It updates runtime path caches only and preserves the specialized K4/doubled-square implementations. After warmup it alternates evaluation order over seven repeats under one-thread BLAS settings. Median old/hybrid time ratios were:

- n=10: 0.903, so the hybrid was slightly slower.
- n=20: 0.972, again slightly slower.
- n=40: 1.104, so the hybrid was approximately 10% faster.

These are local empirical results, not a guaranteed speedup. They include basis-matrix construction and inversions, exclude cold compilation and planning, and concern weight two on deleted cyclic fixtures. The algebraic correctness and storage reduction are independent of these timing outcomes. The producer's different modest warm ratios need not be numerically reproduced because the paths and history differ; neither benchmark justifies a uniform speedup theorem.

## Files and reproduction

- `audit_hybrid_addendum.py`: indexing, brute-force contraction, original-engine comparison, and initial path-history-sensitive timings.
- `hybrid-addendum-receipt.json`: exact counts, error receipts, and initial timings, retained for transparency.
- `profile_hybrid_addendum.py` and `hybrid-profile-receipt.json`: fresh-process per-core path profiling.
- `benchmark_equal_cap_addendum.py` and `equal-cap-addendum-receipt.json`: controlled replanning comparison.
- `HYBRID_INPUT_SHA256SUMS`: pins the implementation and inputs audited here.
- `HYBRID_SHA256SUMS`: hashes only this addendum's new files.

Run the Python scripts with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1`. The initial unconstrained benchmark can take several minutes because it intentionally preserves the original unfavorable path choices. The equal-cap benchmark is the appropriate controlled performance comparison.

No uniform asymptotic remainder or all-orders palette theorem follows from this implementation audit.
