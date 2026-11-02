# Square-guard common blocks and the first cubic response

- **Recovered:** the exact order-25 polynomial and two-block swap already
  proved in `joint_incidence_clocked_frames_20261102/PATH_KERNEL.md`.
- **New proof refinement:** canonical common-block normal form; opposite
  connector conflicts are exactly reversed adjacency between singleton
  blocks. Source-preserving directed block blow-up preserves all Hamilton
  paths, not merely the displayed paths.
- **New finite continuation:** order 26 first permits compatible triples,
  giving `2^291(17179857952+10240z+928z^2+64z^3)`; complete guarded spectra
  are checked through order 30. This counts labelled full tournaments with
  a retained square-sum edge guard, not unrestricted tournament H.

Read `PROOF.md`. `normal.json` contains every order-25/26 path and every
nontrivial compatible subset, plus all finite spectra and union profiles.

## Reproduce anywhere

Python 3.8+ standard library only; run from any working directory:

```sh
python3 /path/to/square_paths/verify.py > /tmp/square-normal.json
python3 -O /path/to/square_paths/verify.py > /tmp/square-optimized.json
cmp /tmp/square-normal.json /tmp/square-optimized.json
cmp /tmp/square-normal.json /path/to/square_paths/normal.json
```

No repository access, network, downloaded dataset, or all-tournament
enumeration is used. The first algorithm enumerates compatible subsets
of all square paths, retaining actual directed-arc unions. The second
counts the active orientation cubes directly, assigning all other guard
orientations response zero. They share the transparent all-start DFS
path bank; the separate auditor independently reconstructs that bank.
All checks are explicit exceptions and survive `-O`.

`SOURCE_PINS.json` authenticates the three small inherited source snapshots
under `sources/`; they are provenance only, not runtime dependencies.
`PUBLIC_FILES.json` is the portable release selector; `SHA256SUMS` checks
its payload. Scratch probes are intentionally excluded.
