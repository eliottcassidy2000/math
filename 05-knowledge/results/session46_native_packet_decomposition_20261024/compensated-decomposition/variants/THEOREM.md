# Exact event, total-work, and peak-clock trade-offs

These are optimizations of the complete finite guarded fibers in the parent certificate, not claims about a larger native supply. All parent hypotheses, resource-availability conditions, and the strict-temporal-nonreturn warning remain in force.

Fix the smallest possible maximum changed-slot support for the stated menu. Connect two legal states when their Hamming distance is within that support. Give an edge its actual sum of changed owner-clock costs. Breadth-first search determines minimum move count; exact integer Dijkstra determines minimum total clock work. A second threshold search determines the smallest maximum single-event clock cost. The checker retains attaining paths and all Dijkstra distance labels.

| Original circuit and buffers | Sharp slot budget | Minimum events | Minimum total clock work | Peak clocks on a minimum-work route | Minimum peak clocks |
|---|---:|---:|---:|---:|---:|
| 16; buffer 3015 | 6 | 4 | 50 | 13 | 13 |
| 16; buffer 1394 | 8 | 3 | 40 | 16 | 16 |
| 16; buffers 1394,2905 | 6 | 4 | 48 | 12 | 12 |
| 18, original cost43; buffer1394 | 6 | 4 | 55 | 17 | 17 |
| 18, original cost37; buffer2071 | 6 | 4 | 49 | 13 | 12 |

In every row one path attains both minimum event count and minimum total work. The original direct atoms cost 32,43,37 clocks respectively, but require 16,18,18 simultaneous changed slots. Reducing atomic support is therefore not a claim of reduced total clock work.

Three additional paths are retained in `../packet_variants.json`:

- The 16-slot circuit with buffer1394 has minimum-work route 65280 -> 116276 -> 86071 -> 255. Its slot supports are 8,6,6 and its clock costs are 16,12,12: total40, improving the original displayed 44-clock route without adding events.
- The cost43 18-slot circuit with buffer1394 has minimum-work route 261632 -> 437825 -> 274217 -> 279403 -> 511. Its total work is55 rather than the first displayed route's59. The checker retains its exact per-step support and clock costs.
- The cost37 18-slot circuit with buffer2071 has peak-clock-optimal route 261632 -> 453144 -> 428849 -> 360765 -> 360700 -> 511. Its slot supports are 6,6,6,4,6 and its clock costs are 12,12,12,9,12, totaling57. The cost-three source2835 changes during the four-slot move, rather than during a six-slot move.

For the last row the four-event minimum-work route costs49 total and peaks at13, whereas peak12 requires at least five events and at least57 total clocks. These latter two lower bounds are exact in this 38-state menu under the sharp six-slot budget. Thus peak coordination, number of events, and cumulative work genuinely conflict; they are not interchangeable statistics.

Nine new full Laplacian image vectors certify the new step charges. Other steps reuse already certified image equations, not resource capacity. Each performed move names its precise image witness. All physical matrices, target margins, shared/private clocks, literal K/converse lifts, selected X mod8 phases, and normalized-resource injection continue to be checked independently.

Every path restores its compensators and consequently reoccupies some released resource. The exact reentry ledger is retained. Under strict temporal nonreturn these paths are disallowed, and the parent's primitive-atom obstruction applies instead.
