# Scale-three Hamming-six depth-three/full-run benchmark

Date: 2026-07-15.  This is a planning benchmark, not a census or an emptiness
claim.  It used the optimized executable built from
`lrc13_scale_three_hamming_six_depth_two_scout_codex_S16.cpp`, source SHA-256
`4c6144e14d12a4badf734d26844a26980dadd29f5553e1c7927bb23a74d66ca9`.

The exact depth-two ledger supplied light, median, and heavy contexts within
each stratum.  The command was

```bash
/usr/bin/time -l /tmp/thm862-d2-final-O3 \
  --depth 3 --context-start INDEX --context-limit 1
```

`--depth 3` deliberately materializes every depth-three child and does not
activate the THM-857 full-tooth/streaming gates; those gates are activated by
the engine only for `--depth 6`.  User time is the useful comparison because
other repository computations distorted wall time during the heavier pilots.

| stratum | load | context | depth 2 | depth 3 | user s | wall s | peak RSS |
|---|---|---:|---:|---:|---:|---:|---:|
| `1^2 3^4` | light | 103 | 3,011 | 97,386 | 1.29 | 2.42 | 12,550,144 |
| `1^2 3^4` | median | 1,287 | 8,797 | 480,983 | 10.12 | 29.48 | 32,358,400 |
| `1^2 3^4` | heavy | 1,184 | 26,093 | 2,453,403 | 79.49 | 387.05 | 74,350,592 |
| `1 3^5` | light | 1,444 | 4,530 | 179,258 | 2.74 | 4.10 | 17,563,648 |
| `1 3^5` | median | 1,120 | 8,240 | 435,801 | 8.47 | 19.30 | 29,179,904 |
| `1 3^5` | heavy | 1,057 | 22,800 | 1,994,613 | 59.23 | 303.63 | 72,728,576 |
| `3^6` | light | 1,452 | 4,600 | 182,026 | 2.72 | 5.20 | 18,661,376 |
| `3^6` | median | 395 | 9,828 | 566,682 | 12.34 | 50.34 | 33,488,896 |
| `3^6` | heavy | 895 | 33,133 | 3,497,310 | 118.96 | 492.69 | 98,631,680 |

The depth-three/depth-two ratios range from `32.34` to `105.55`.  Scaling the
three stratum totals by the light, median, and heavy ratios respectively gives
nonrigorous all-context projections of `568.6M`, `823.4M`, and `1.426B`
depth-three logical nodes.  Scaling the measured user cost in the same way
gives approximately `2.33`, `4.72`, and `12.65` CPU-hours for the raw,
uncached `--depth 3` implementation.  These are workload projections only.

One bounded on/off check at context 103 gave the same `97,386` depth-three
nodes with and without `--no-early-cap-gate`, as expected: bounded depth-three
mode intentionally leaves the gates inactive.

The full-depth gates were then exercised on one light context per stratum:

| stratum | context | depth 3 | full-tooth at d3 | stream at d3 | depth 4 | later nodes | user s | peak RSS |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| `1^2 3^4` | 103 | 97,386 | 87,184 | 9,629 | 1,143 | none | 0.21 | 8,699,904 |
| `1 3^5` | 1,444 | 179,258 | 163,966 | 14,568 | 1,671 | none | 0.34 | 11,452,416 |
| `3^6` | 1,452 | 182,026 | 168,177 | 12,895 | 2,341 | one d5 node | 0.34 | 9,076,736 |

All three full-depth pilots had zero covering and loose terminals.  The gates
certified more than 99% of depth-three children immediately, and every
remaining branch died by depth five in these samples.  Extrapolating the
light gated cost over the three raw-node projections suggests roughly 18--46
CPU-minutes, but the light samples do not bound heavy contexts.  A cautious
planning range for the present full-depth engine is therefore 0.5--3
CPU-hours, not a theorem.

## Recommendation

Build a geometry-batched depth-three engine before an all-context run.  The
exact depth-two computation has `14,992,263` logical lanes but only `4,307,561`
literal `(R,x1,x2)` geometries, a factor `14,992,263/4,307,561 = 3.48045...`.
The parent components, complete-tooth test, and intersection with a proposed
third speed depend only on that geometry (and the third speed); candidate-ray
availability and the least later speed remain lane-specific.  Thus the cache
can eliminate repeated literal interval work while retaining every labelled
language and its shortcut ancestry.  It is an execution cache, not a
proof-state quotient.
