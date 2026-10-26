# The reverse valuation potential on the production ordinary A000568 graph, and fixed-precision pruning

2026-10-25, klein. **IMPLEMENTED and FINITE-EXACT** for `N <= 100`; **PROVED** pruning statement.
Source construction: [REVERSE_POTENTIAL](../../port_residual_integral_20261025/counting/REVERSE_POTENTIAL.md)
(proved there; implementation declared OPEN). Engine formulas (I5)–(I7) with `e = 0`, one channel,
from [counting/README](../../port_residual_integral_20261025/counting/README.md) and its `check.py`
(`bound`, `intrinsic`, `transfer`, defect update `D'(d) = D(d) + k(d - gcd(ell,d))`). No cycle-index
engine is imported; the graph is rebuilt from those formulas.

## 1. What was implemented

`reverse_potential.py` builds the structurally reachable, terminal-reaching graph with vertices
`(ell, S, D)` (stage, remaining degree, retained defect tuple), edges labelled by the odd factor
`T_e = transfer(N-S, ell, k)` and the cost `J_e` of (I7). A reverse pass computes
`V(v) = min_e (J_e + V(w))`, `B = V(initial)`, reduced costs `J'_e = J_e + V(w) - V(v) >= 0`, and the
transformed recurrence `Z` modulo `2^K` from `Z_initial = 1`. The terminal identity
`Z_terminal = O(N) a(N) / 2^B (mod 2^K)` is checked against the exact partition oracle.

| N | vertices | edges | `B` | `floor(N/2)` | `v_2(a(N))` | identity mod `2^20` |
|---:|---:|---:|---:|---:|---:|---|
| 10 | 20 | 28 | 5 | 5 | 6 | yes |
| 20 | 80 | 142 | 10 | 10 | 12 | yes |
| 30 | 213 | 472 | 15 | 15 | 17 | yes |
| 40 | 471 | 1,229 | 20 | 20 | 26 | yes |
| 50 | 927 | 2,732 | 25 | 25 | 26 | yes |
| 100 | 13,445 | 50,044 | 50 | 50 | — | yes (`K = 1..64`) |

`B = floor(N/2)` holds at every tested `N >= 3` (the (R3) identification), and `B < v_2(a(N))` at
even `N = 8, 10, 12, 20, 30, 40, 50` (carries among minimum-cost terms), exactly as REVERSE_POTENTIAL
warns. The graph at `N = 100` has the 13,445 states and 50,044 contributions that the packet reports
for its `D`-gauge response graph, so this is the production graph, not a toy. The naive hostile
"subtract `B` at the initial edge" occurs on the real graph at every `N`: one initial edge has
`J_e < B` (receipt field `naive_initial_negative_edges`).

## 2. Fixed-precision pruning (PROVED)

Modulo `2^K`, an edge with reduced cost `J'_e >= K` contributes `T_e 2^{J'_e} Z_v ≡ 0`. Hence:

**Proposition.** Let `G_K` be the subgraph of edges with `J'_e < K`, restricted to vertices that
are forward-reachable from the initial vertex and backward-reachable from the terminal inside it.
Running the transformed recurrence on `G_K` modulo `2^K` returns the same `Z_terminal mod 2^K`.
*Proof.* Deleting an edge of reduced cost `>= K` changes no `Z_w mod 2^K`; a vertex not on an
initial-to-terminal path of `G_K` cannot contribute to `Z_terminal`. QED.

Measured live subgraphs (`reverse_prune.py`, all results correct against the exact oracle):

| N | full V / E | K=1 | K=4 | K=8 | K=16 | K=32 | K=64 |
|---:|---|---|---|---|---|---|---|
| 60 | 1,721 / 5,499 | 150 / 156 | 443 / 529 | 797 / 1,113 | 1,143 / 2,041 | 1,420 / 3,048 | 1,590 / 3,980 |
| 80 | 5,120 / 18,109 | 360 / 374 | 713 / 831 | 1,791 / 2,483 | 2,920 / 5,292 | 4,041 / 9,072 | 4,668 / 12,322 |
| 100 | 13,445 / 50,044 | 550 / 568 | 1,361 / 1,577 | 3,769 / 5,270 | 6,857 / 12,495 | 10,087 / 23,248 | 11,994 / 32,761 |

So sixteen normalized binary digits of `a(100)` need 51% of the states and 25% of the transitions of
the unrestricted recurrence, eight digits 28% / 10.5%, and the parity of `a(N)/2^{floor(N/2)}` needs
a graph with fewer than six states per unit of `N` (`K = 1`: 550 states at `N = 100`). This is a
smaller **unrestricted** recurrence for normalized digits, in the packet's own sense: no tournament
filter is imposed; only arithmetically invisible transitions are removed. The reverse pass over the
full graph is still needed once to compute `V`; the pruned graph can then be reused for every
smaller `K`.

## 3. Hostiles and limits

- Structural reachability, not nonzero residues, defines the graph: a live state may have `Z ≡ 0`.
- `G_K` is not valid for `K' > K`; the baseline `B` is the minimum summand valuation, not `v_2(a(N))`.
- The signed clone-excess recurrence has its own graph and sign bookkeeping; it is not implemented.
- No asymptotic bound on `|E(G_K)|` is proved; the table is measurement. The exact oracle
  (`ordinary_direct`) enumerates odd partitions and is the cost bottleneck above `N ≈ 100`, not the
  reverse pass.
