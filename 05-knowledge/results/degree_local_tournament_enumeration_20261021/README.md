# Degree-local tournament enumeration checkpoint

New Session22 arithmetic and executable evidence, independently reviewed. This is not a restoration of the lost Session21 source or timings.

- [Proof, implementation and producer evidence](session22-enumeration/README.md)
- [Independent review and telemetry-matched ablation](session22-enumeration-review/README.md)

The kernel uses degree-local odd-factorial normalization with integral precomputed transfers. Mathematical response states and transitions are unchanged. Shared-host telemetry-matched median speedups were 1.83–2.16 at ordinary bounds 80–200 and 1.33 at core bound 100. These are finite measurements, not universally-fastest or asymptotic runtime claims. Signed core aggregates and rational no-fixed observers retain their original meaning; individual tournament classes and metagraph adjacency are not recovered.

## Reproduce

From this packet directory:

```sh
python3 session22-enumeration/check.py --receipt /tmp/degree-local-producer.json
python3 -O session22-enumeration/check.py --receipt /tmp/degree-local-producer-opt.json
python3 session22-enumeration-review/check.py
python3 -O session22-enumeration-review/check.py
python3 session22-enumeration-review/benchmark.py
python3 session22-enumeration-review/benchmark_extra.py
```

The sibling directory names are deliberately preserved so that independent scripts run without path edits. Existing receipts preserve their original observations; reruns may have different timing samples. Both benchmark scripts regenerate the included telemetry-stripped baseline locally.

## Immutable source provenance

Baseline is byte-identical to `05-knowledge/results/tournament_enumeration_20261021a/check.py` at commit `0bc8da01d9e5e2a72d2c3f410aa07d5245e61967`, inherited from `881a9903995884c02a03ffb45c68e9174f8b6c58`. Its SHA-256 is `8ef8429ff53764dc54af36ac45eb1e2729b7cf8e67c0b3d1309b8e9a308fcdbd`. The reviewed new kernel SHA-256 is `b76167ef82dda0b134681abce33e159d5dd80d15e92edb9e31ab9b490f6939ad`.

The top-level manifest pins every included file except itself. Nested producer/reviewer files are copied without rewriting; their historical commands mentioning a workspace root now refer to this packet root.
