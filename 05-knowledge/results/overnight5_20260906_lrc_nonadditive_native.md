# Independent native confirmation of the nonadditive ceiling

**Status: FINITE-EXACT independent confirmation.** The all-height theorem is
already proved by the [combined fourth-checkpoint corollary](overnight4_20260906_lrc_parityfree_native.md#5-combined-nonadditive-ceiling-and-exact-old-threshold-classification),
using incoming THM-4437/4441 and the independently audited norm-four theorem.
This larger native run is additional evidence, not its missing dependency.

The exact universe is every primitive sorted distinct positive ternary-unit
triple `1<=a<b<c<=535`, excluding `a+b=c`. The implementation retains the
literal six-sheet interval engine independently audited in the fourth
checkpoint, without using raw carrier congruences or roof formulas to
compute a projection. All endpoints use the exact denominator `42abc`;
comparisons use widened integer products. Nonprimitive tails are outside
this census, although the all-height theorem separately descends by Haar
covering.

The [frozen output](overnight5_20260906_lrc_nonadditive_native_h535.out) gives:

| Quantity | Exact value |
|---|---:|
| Primitive unit triples | 6,514,176 |
| Excluded additive triples | 10,908 |
| Complete nonadditive rows | 6,503,268 |
| Empty dictionaries | 44,018 |
| Literal sheet contacts | 638,334,708 |
| Maximum selected minimum | 11/140 |
| Maximum physical mass | 11/140 |
| Equality triple for either maximum | (2,11,20) |

At the extremizer the three projections are
`(131/1540,11/140,3/35)`. Outside both norm-four and norm-five relations,
all 6,462,442 rows have every projection strictly below `11/140`; the
largest individual projection in this head is `6/77`, first at `(7,16,22)`.
This preserves incoming THM-4437's non-strict individual boundary.
The exact row-serialization FNV1a64 checksum is `17918696151674804143`.

The [cutoff/universe referee](overnight5_20260906_lrc_nonadditive_reduction_audit.md)
independently proves the generic `c>=535` reduction from the pinned
coefficient box, and reconstructs all universe counts by divisor/Mobius
methods. Its 2,661 gates do not independently recompute these masses.
The [norm-five reconstruction](overnight5_20260906_lrc_norm5_profiles.md)
and its [independent profile referee](overnight5_20260906_lrc_norm5_audit.md)
provide compact complete bases and analytic tails for that exceptional
sector. The current shorter all-height proof uses the stronger incoming
generic theorem and does not need the new large head.

The unchanged native engine previously agreed row by row with the raw H63
table. This fifth H535 run was not paired with a fresh full raw engine;
incoming THM-4437 separately has an independently matched H611 head. The
frozen native output's final sentence records its pre-integration status;
the scope stated here is current. Root corrected the copied header and
normalized line endings when filing, without changing executable logic.

Reproduce from the repository root, writing the executable outside the repo:

```bash
g++ -std=c++17 -O3 -DNDEBUG 04-computation/overnight5_20260906_lrc_nonadditive_native.cpp -o C:/w/overnight5_norm5/nonadditive_native.exe
C:/w/overnight5_norm5/nonadditive_native.exe 535
python3 04-computation/overnight5_20260906_lrc_nonadditive_reduction_audit.py
```

The fifth manifest pins raw LF hashes. This is a local tail statistic;
arbitrary body entry, synchronization, and LRC(14) remain open.
