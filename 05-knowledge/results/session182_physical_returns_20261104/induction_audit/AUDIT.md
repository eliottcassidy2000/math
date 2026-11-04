# Independent audit: ordinal-component multiport construction

PASS on mathematical claims, with scope as stated in the source proof.

1. Terminal-SCC criterion is necessary and sufficient for full strongness, including arbitrary nonstrong suffixes and empty port sets. The converse routes through the actual terminal port, root, and nonempty prefix.
2. All seven ordered-sum formulas are correct. In particular the four joint-endpoint placements contribute jH, hJ, eS and gF; omitting gF would be incorrect. Empty component port sets are handled.
3. The all-clock +4 construction has the stated two actual signatures, negative intercepts −5 and−11, defects6 and12, and physical old-name/edge retention except the expressly changed root–1 edge. It is an infinite one-branch construction, not iterable two-branch closure or a native guarded program.
4. The bit recurrence and leading-nonport shifts preserve the stated joint defect; the final port bit is precisely the full strongness condition.

The author verifier ran successfully both normally and under Python -O, producing identical receipts:1380 ordinal checks,1098 strongness checks,1098 attachment/defect checks,62 binary signature controls and four full +4 example pairs. Receipts are `normal.json` and `optimized.json`.

An independent verifier `audit.py` enumerates two-slot subset assignments directly for every port bitmask on TT_m through m11, giving4094 recurrence checks. This algorithm uses no author signature routine. It independently verifies all-strong physical retention schedules for eight prefix lengths: add terminal port4, add nonport3 before4, then reverse root–1. The first two steps preserve every old incidence, and the last changes exactly root–1. This establishes an adjacency-level strong schedule, not native Hamilton-complement legality. Receipt: `independent.json`.

No forbidden file, startup helper or installation used. Source hashes are retained in `source_sha256.txt`.

## Latest strengthened schedule re-audit

Re-read the updated Section5 paragraph after the author incorporated the all-strong order add4, insert3, reverse root–1. PASS: this is exactly the schedule independently checked in `audit.py` for eight prefix lengths. At each birth the previously existing induced graph is unchanged, terminal4 remains an actual port, and the final reversal does not remove terminal4. The stated all-order strongness follows from the proved terminal-SCC criterion. The distinction from native guarded legality remains explicit. The author verifier's SHA256 is unchanged from the successful normal/-O runs, so no verifier rerun is required.
