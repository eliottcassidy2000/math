# Mixed-(2,2), four-common-neighbor exact exclusion

Entry point: `mixed22_c4_exact_exclusion_20261019.md`.

The completed packet checks all 1,165,395 degree-valid canonical cores:
950,128 analytic exclusions and 215,267 residual cores checked order by order.
`independent_full_replay_receipt_20261019.json` records the terminal audit.
The four compressed residual streams are partitioned by index modulo four.
Each is independently hash-checked; the verifier checks their exact union.

Reproduce with Python 3.12 (tested on 3.12.14), standard library only:

```sh
python -O verify_mixed22_c4_exact_exclusion_20261019.py
```

For four-way replay, use the commands in the theorem note, then run its
`--assemble` command. All checks remain active under `-O`. Replay regenerates
receipt files with new elapsed times; compare their input hashes and exact
counts rather than the runtime-dependent receipt bytes.

The complete archive is self-contained. The two smaller delivery archives
are an alternative to it: extract both into the same parent directory.

```sh
tar -xzf mixed22_c4_exact_exclusion_20261019_delivery_part1of2.tar.gz
tar -xzf mixed22_c4_exact_exclusion_20261019_delivery_part2of2.tar.gz
cd mixed22_c4_exact_exclusion_20261019
python -O verify_mixed22_c4_exact_exclusion_20261019.py
```

The manifest covers the assembled directory. Earlier analytic-stage notes
and receipts are retained as provenance; the full theorem and terminal replay
receipt give the completed scope. Previous c<=2/c=3 packets are linked only
for the combined branch corollary and are not runtime dependencies here.
