# Independent audit: event, work and peak-clock optima

**FINITE-EXACT / PASS.** Two separate new checkers work from the already independently authenticated complete legal fibres. They do not import the producer's search, optimizer or image checker.

`audit_optimization.py` builds the complete pairwise state graph using literal physical-path clock distances. Independent minimax Dijkstra reproduces the original threshold-scan bottlenecks. Breadth-first search obtains minimum event counts, while integer Dijkstra obtains minimum cumulative owner-clock work. Both are evaluated at the exact slot-optimal and clock-optimal thresholds.

| Circuit and buffers | Minimum slot peak | Minimum events at that slot peak | Minimum total clock work | Minimum clock peak |
|---|---:|---:|---:|---:|
| 16 with 3015 | 6 | 4 | 50 | 13 |
| 16 with 1394 | 8 | 3 | 40 | 16 |
| 16 with 1394,2905 | 6 | 4 | 48 | 12 |
| 18, original cost43, with 1394 | 6 | 4 | 55 | 17 |
| 18, original cost37, with 2071 | 6 | 4 | 49 | 12 |

For each row a path achieves both the minimum event count and minimum total work at the optimal slot peak. In the final row, requiring the optimal clock peak twelve instead forces at least five events and total work at least fifty-seven. Both bounds are attained. Therefore the four-event, work49, peak13 route and five-event, work57, peak12 route are distinct optimal trade-offs; neither statistic may be silently substituted for another.

The optimizer passes **55 explicit checks** under both ordinary and optimized Python with byte-identical receipts. Full attaining paths and step weights are retained in `optimization.normal.json`.

## Native validation of all three new variants

`audit_variants.py` checks each new state against the complete authenticated fibre, then independently recomputes all 3,528 physical current rows for each authentic generator. It also checks literal target margins, every fixed physical-owner clock, the full four-port resource orbits and simultaneous normalized-resource injection. The underlying actual matrices, K/converse closures and source-private owner symbols are tied by hash to their independent primary audit.

- 16 with buffer1394: 65280 → 116276 → 86071 → 255. Slot supports 8,6,6; clock costs 16,12,12; three events and total work40. It improves the original work44 route with the same event count and peaks.
- Cost43 eighteen-source circuit with buffer1394: 261632 → 437825 → 274217 → 279403 → 511. Four six-slot events and work55 instead of59, with clock peak17.
- Cost37 eighteen-source circuit with buffer2071: 261632 → 453144 → 428849 → 360765 → 360700 → 511. Slot supports 6,6,6,4,6; clock costs 12,12,12,9,12; five events, work57 and clock peak12. The cost-three source2835 is changed inside the four-slot event.

All **nine newly supplied full Laplacian image vectors** are multiplied independently over the original 438,216 edge rows, checking every coordinate of `Lz=charge+256correction`. Together with the unchanged twenty previously checked vectors, they cover the exact source charge of every new performed step; the common target charge remains covered by the prior audit. Hash checks establish that the old image vectors are the same bytes previously independently multiplied, rather than trusting new assertions about them.

The variant audit passes **32,447 explicit checks** under both ordinary and optimized Python, with byte-identical receipts. Its new-state resource-reentry counts are 3,5,9 respectively. All buffers return exactly, but none of these routes is strictly temporally nonreturning. Their stronger no-reacquisition interpretation remains blocked by the primitive-atom proof in ROUTE_AUDIT.md.

## Scope and reproduction

All optimization minima are in their specified complete pinned legal menus. The global full28 invariant proves the separate six-slot lower bound, but it does not prove all these event or total-work minima in a larger menu. Actual resource availability and a compatible unchanged remainder remain hypotheses; no full-atlas rewiring or all-size native supply is asserted.

Run `audit_optimization.py` with `--packet`, `--routes` (the independent primary receipt), and `--out`. Run `audit_variants.py` with `--repo`, `--packet`, `--variants`, `--images` (the new nine-vector file), `--base-images` followed by the three original vector files, `--routes`, `--optimization` (the independent optimizer receipt), and `--out`. Every input file is SHA256-pinned in the resulting receipt. Repeat with `python -O` and compare bytes.
