# Portable GSMA / norm / sparse-unit packet

Read `INHERITANCE_AND_KERNEL.md` first. It separates recovered source material,
proposed GSMA architecture, proved algebra, cited finiteness, and the finite
inverse-power census at orders 3, 4, 5.

* `SOURCE_PINS.json`: 10 full-path repository sources pinned to commit
  `522535a7b6f519fd0c1e4304a2c5b1c909fce4c9`, with blob and SHA-256 identities.
* `check_unit_kernel.py`: standard-library-only exact checker; no repository imports.
* `unit_kernel.normal.json`, `unit_kernel.optimized.json`: identical 15,544-check receipts.
* `RELOCATION_VERIFICATION.json`: same execution after relocation.
* `PRIMARY_SOURCE.json`: inspected primary source for Thue–Mahler terminology and finiteness.
* `SHA256SUMS` and `verify_release.py`: portable file integrity and executable replay.

Run `python3 -B verify_release.py` after copying this folder. No network,
external packages, source repository, or original absolute path is needed.

The `k=4,5` sparse-unit tables are complete only for `0 <= n <= 500`.
They are not full Thue-equation solution lists. The cubic global classification
and CF-305 certification belong to the companion Tribonacci packet.
