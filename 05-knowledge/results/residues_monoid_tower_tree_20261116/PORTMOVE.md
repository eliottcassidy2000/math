# P107: the root-moving port move on the ready seeds (FINITE-EXACT, complementing wave 46)

Session klein-2026-11-16-monoid. Tool `scripts/portmove.c` (all rooted tournaments of an order, or a given seed; every
extension by one or two vertices with arbitrary contacts; every new root); scripts `seeds.py`, `seeds2.py`,
`onparent.py`; receipts `portmove_3..6.out`, `seeds*.out`, `onparent.out`, `iterate_children.out`, `seedruns.out`.

**Coordinates (wave 45 §6).** For a root `x`: `r = E`, `s = O` (the even and odd root-start path sums), `g = H(T − x)`,
charge `a = 2r + 1 − H`, defect `b = s − r − g`. A root is ready when `b = 0`. The defect is `b = −(D + g)` in the
notation of CONTINUATION.md (path_splice_rows_continuation_20261116), so the law `D + g ≥ 0` found there is `b ≤ 0`
for every rooted tournament to order 8, and the ready roots are exactly the roots with `D = −g` — the class that
carried every literal digit hit. P107 asks for an operation that may move the root and sends a ready `(r, a)` to a
ready `(2r + 1, a)`, i.e. `H′ = 2H + 1 + a`, `a′ = a`, `b′ = 0`.

**Already resolved (wave 46 §1, PROVED / FINITE-EXACT).** No single module, root moving or not, is a binary port
move; in one-parent growth no `{0, 2}` binary cylinder exists to order 8. Not repeated.

**This test.** Arbitrary contacts, one or two added vertices, every new root:

| base | ready roots | extensions with `H′ = 2H + 1 + a` | (extension, new root) pairs ready and same charge |
|---|---|---|---|
| all classes of order 3 | 4 | 20 | 21, all from `H = 1`, charge 0 (`1 → 3`) |
| order 4 | 7 | 112 | 27, all `H = 1`, charge 0 |
| order 5 | 17 | 255 | 33, all `H = 1`, charge 0 |
| order 6 | 64 | 568 | 39 from `H = 1`, charge 0; **3 from `H = 15`, charge 4** |

The three nontrivial moves are from the order-6 base with gentourng string `111011111111101`, root 4, profile
`(15, 9, 10, 1)`, charge 4: two added vertices (`cx = 57, cy = 0` or its relabellings) give `(35, 19, 24, 5)` at the
same root, charge 4, ready — `r: 9 ↦ 19 = 2·9 + 1`. They are one-offs: from the child no extension by one or two
vertices even reaches the next target `75 = 2·35 + 1 + 4` (`iterate_children.out`). So no binary cylinder starts there.

**On wave 45's ready seeds 13, 23, 41, 59, 75, 95.** The caterpillar automaton (`seeds.py`, verified exactly against
the rooted profile on all 511 words of length ≤ 8) has, among words of length ≤ 13, a ready root at the end only for
75 (`SLLSLLSL`, charge −24); the one-parent trees (`onparent.py`) have ready roots at all six counts at orders 5–8,
with many charges (13 at charge 0; 23 at −4, 0, 2, 6, 8; 41 at −12 … 26; 59 at −18 … 32; 75 at −14 … 40; 95 at −30 …
48). For all 83 representatives (one per count and charge, orders 5–8; two added vertices at orders ≤ 7, one at order 8; every new root) the search finds **no** ready same-charge move: 51 seeds have no extension of the target count at all, the other 32 reach the count `2H + 1 + a` (5 to 21 extensions each) but never ready at the same charge (receipt `seedruns.out`). Pre-registered P-6 and P-9 **confirmed**.

**Reading.** The ready class is `D = −g`; appending at the root freezes the charge and doubles `r`; the port move
`r ↦ 2r + 1` exists sporadically (charge 4, one base, one step) but not as a reusable letter — consistent with wave
45's charge principle (a reusable digit needs a shared conserved charge) and with the previous packet's finding that
synchronized digits exist only at count 5. The honest statement of P107 after this test: **no root-moving port move
of this kind exists on the ready seeds with up to two added vertices** (83 seed representatives, every count and charge, 0 hits).
