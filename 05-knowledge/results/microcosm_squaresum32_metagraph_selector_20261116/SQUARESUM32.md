# The first 32 nodes of the square-sum graph: paths from 25, the first cycle at 32, and the analogy

Owner: *"Consider how this relates analogously to ideas around the first 32 nodes of the sum-square graph."* The
square-sum graph on `1..n` joins `i, j` when `i + j` is a perfect square (THM-4505, THM-4860 for the prime version).
Counts by depth-first search ([receipts/hamcyc.out](receipts/hamcyc.out), `hamcyc.c`, under a second):

| n | 15 | 16 | 17 | 18–22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 | 33 | 34 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Hamiltonian paths | 1 | 1 | 1 | 0 | 3 | 0 | 10 | 12 | 35 | 52 | 19 | 20 | 349 | 392 | 669 | 4041 |
| Hamiltonian cycles | – | – | – | – | – | – | – | – | – | 0 | 0 | 0 | 0 | **1** | 1 | 11 |

So the square-sum graph has a Hamiltonian path for every `n ≥ 25` in this range after the holes 18–22 and 24, and its
**first Hamiltonian cycle appears exactly at 32**, uniquely (M-4 confirmed; the path counts agree with the published
sequence). The analogy with the count spectrum, typed:

| | square-sum graph on `1..n` | odd Hamiltonian-path counts |
|---|---|---|
| microcosm | paths only at 15, 16, 17, 23, then from 25; cycles only from 32 | counts 1, 3, 5, then 9–19, then from 23; the separators 7 and 21 |
| the two separators | the gaps 18–22 and 24 (six values) | 7 and 21 (two values) |
| the 32 | `1..32` is the first node set with a Hamiltonian cycle (one cycle) | the first 32 odd numbers end at 63 = the atom through both holes; they sum to 32² |
| macrocosm | conjectured: paths for all `n ≥ 25`, cycles for all `n ≥ 32` (verified far; not proved here) | conjectured: every odd `≥ 23` is a count (Chan–Pak; verified to 10⁸ by the banks) |

Both are "finitely many small exceptions, then everything" statements whose macrocosm is verified but not proved;
nothing structural joins them beyond the shape, and the shared number 32 is `2^5` on one side (the first cycle) and
`2^5` odd numbers on the other. **Type: ANALOGY.**
