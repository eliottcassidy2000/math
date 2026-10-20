# Hollow-frame blue extrema

Read `HOLLOW_FRAME_BLUE_EXTREMA.md`.

Main results:

- Exact labelled external-frame mean `E h=k!`.
- Gluck's regular-subset theorem gives `max B=max h`.
- Hence `max B>=k!+1` for `k>=2`, ruling out every fixed exponential cap.
- An explicit elementary family independently has `B>=3(m!)^3`.
- Complete independent censuses give maxima `1,3,9,39,195` for orders
  `3,5,7,9,11`.
- The order-eleven maximum is attained by exactly eight self-converse
  classes, one auxiliary `Q3` packet.

Scope and status are separated in the main manuscript. Exact variances
`0,1,6,75,1500` are finite data, not an all-orders formula.

Artifacts:

- `enumerate_frames.cpp`: complete subcube-sum census.
- `replay_subset_dp.cpp`: independent pointwise subset-DP replay.
- `verify_extrema.py`: inherited-class and witness checks.
- `build_witness.py`: explicit full-word witness and maximum-orbit checks.
- `exact_k*.json`, `replay_k*.json`, `verification.json`: matching outputs.
- `h_k*.bin`: complete raw uint16 counts, little-endian on the recorded host.
- `n11_witness.json`: all 195 symmetric words, adjacency and group certificate.
- `PUBLICATION_MANIFEST.json`, `SHA256SUMS`: frozen artifact inventory.

Reproduction commands are relative to this packet's directory. The
inherited-class census is optional; use `verify_extrema.py --inherited PATH`
to select it explicitly. All core checks run without a sibling census.

No repository or remote changes were made by this worker. All work is staged
in this directory for the parent to review and publish.
