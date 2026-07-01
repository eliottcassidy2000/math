---
id: HYP-3803
title: THE FLIP-RANK of tournament iso classes + the BALANCED-MAX-CUT shape + a PHASE TRANSITION at n=6. In the tiling/cube model a labeled tournament = a bit per edge (i<j); a REALIZING SUBCUBE of dimension k = k "free" edges + a fixed orientation of the other C(n,2)-k edges whose 2^k completions meet EVERY iso class. FLIP-RANK rho(n) = min such k (owner's n=4 "flip only 2 arcs"). VERIFIED (exhaustive n=4,5,6): rho(n) = 1,2,4,7 for n=3,4,5,6. The information floor LB=ceil(log2|G_n|)=1,2,4,6 is TIGHT for n<=5 but NOT at n=6 (rho(6)=7=LB+1, exhaustive: no 6-edge subcube realizes G_6). THE SHAPE of the minimal config for n<=5 = FIX A BALANCED VERTEX-BIPARTITION CUT (the complete bipartite K_{a,b}, a=ceil(n/2), b=floor(n/2)) and FLIP THE TWO SIDES (the within-part sub-tournaments); free-count f(n)=C(a,2)+C(b,2). STRIKING: f(n)=LB EXACTLY for n=3..7 (1,2,4,6,9), then f(n)<LB for n>=8 (12<13, ...), so the bipartite config becomes information-theoretically IMPOSSIBLE at n>=8. But it BREAKS EARLIER, at n=6: the balanced-bipartite config REALIZES for n=4 (fix K_{2,2}, free the matching {(0,1),(2,3)}) and n=5 (fix K_{3,2}, free triangle{0,1,2}+edge{3,4}), but FAILS for EVERY split of n=6 (3+3 f=6, 2+4 f=7, 1+5 f=10 all fail for all cut orientations). So the max-cut shape works for n<=5 and dies at n=6; the n=6 achiever (rho=7) is a NON-bipartite connected 7-edge config (degrees 2,2,2,2,3,3). Redundancy 2^rho-|G_n| = 0,0,4,72 for n=3,4,5,6 (n=3,4 are EXACT bijections since |G|=2,4 are powers of 2). Related notions defined: transversal-rank (exact bijection, only n=3,4), the cut-orientation "mixing" condition, the base-path-anchored (naive, k=C(n-1,2)) vs free
status: MIXED (exhaustive small-n verification + clean shape + open general-n). VERIFIED EXACTLY: rho(3,4,5,6)=1,2,4,7 (n=4,5 exhaustive over all subcubes; n=6 k=6 EXHAUSTIVE over all C(15,6) free-sets x fixed orientations => none realizes => rho(6)>=7, and a 7-subcube realizes => rho(6)=7). Balanced-bipartite realizes n=4,5 (specific cut orientations), fails ALL splits n=6. f(n)=C(ceil(n/2),2)+C(floor(n/2),2)=LB for n=3..7 (arithmetic fact), f(n)<LB for n>=8. HONEST: rho(n) for n>=7 NOT computed (n=7 canonicalization = 2^21 x 5040 too slow here); the "shape at n>=6" (non-bipartite) not yet characterized beyond one example; the n=6 transition is verified, its general-n continuation is conjectural. A structural discovery, not tied to a proof of LRC.
source: klein-2026-07-01-S71
depends_on:
  - HYP-3802   # S70: tiling model on the atoms; circulant=>difference-striped (this studies the general subcube shape)
related:
  - HYP-3801   # loop-map dictionary (the cube Q_m and its S_n action)
  - THM-002    # OCF; the iso classes G_n being coded here
external: A000568 (# tournaments up to iso); the tournament cube Q_{C(n,2)} with S_n action; max-cut / vertex bipartition; critical/defining sets in combinatorial design
results:
  - 04-computation/flip_rank_realizing_subcube_klein.py
  - 04-computation/flip_rank_bipartite_shape_klein.py
  - 04-computation/flip_rank_settle_n6_klein.py
  - 05-knowledge/results/flip_rank_realizing_subcube_klein.out
  - 05-knowledge/results/flip_rank_bipartite_shape_klein.out
  - 05-knowledge/results/flip_rank_settle_n6_klein.out
---

# HYP-3803 — the flip-rank, the balanced-max-cut shape, and the n=6 transition

## Definition
A labeled tournament on `n` vertices = a bit per edge `(i<j)` (`1: i->j`). A **realizing subcube** of
dimension `k` = a set of `k` "free" edges + a fixed orientation of the other `C(n,2)-k` edges, such that
the `2^k` completions include a representative of every iso class (`G_n`, A000568). The **flip-rank**
`rho(n)` = the minimum such `k`. (This makes precise the owner's n=4 observation: all 4 iso classes are
reachable by flipping only 2 arcs, if the 4 fixed arcs are configured right. The naive Hamiltonian-path
model uses `k = C(n-1,2)`; the flip-rank is the true minimum.)

## Verified results (exhaustive small n)
| `n` | `\|G_n\|` | `LB=ceil(log2)` | **`rho(n)`** | redundancy `2^rho-\|G\|` |
|---|---|---|---|---|
| 3 | 2 | 1 | **1** | 0 (exact bijection) |
| 4 | 4 | 2 | **2** | 0 (exact bijection) |
| 5 | 12 | 4 | **4** | 4 |
| 6 | 56 | 6 | **7** | 72 |
`rho(6)=7` is **exhaustively verified**: no 6-edge subcube (any free-set, any fixed orientation) realizes
all 56 classes, and a 7-edge one does. So **the information floor `LB=ceil(log2|G_n|)` is tight for
`n<=5` but breaks at `n=6`** (`rho(6)=7 > 6=LB`).

## The shape: fix a balanced max-cut, flip the two sides
For `n<=5` the minimal realizing subcube has a clean shape: **fix the complete bipartite cut `K_{a,b}` of a
balanced vertex bipartition `A|B` (`a=ceil(n/2)`, `b=floor(n/2)`), and free the within-part edges** (the two
sub-tournaments on `A` and `B`). Free-count `f(n) = C(a,2)+C(b,2)`:
- **n=4**: parts `{0,1}|{2,3}`, fix `K_{2,2}` (4 arcs = the cut = a 4-cycle), free the matching
  `{(0,1),(2,3)}` (2 arcs). This IS the owner's "flip 2 arcs, special config on the 4 fixed" — the config
  rule is *the fixed arcs form the bipartite cut*.
- **n=5**: parts `{0,1,2}|{3,4}`, fix `K_{3,2}` (6 arcs), free `triangle{0,1,2} + edge{3,4}` (4 arcs).

**Arithmetic coincidence**: `f(n) = C(ceil(n/2),2)+C(floor(n/2),2)` equals `LB=ceil(log2|G_n|)` **exactly for
`n=3..7`** (`1,2,4,6,9`), then `f(n) < LB` for `n>=8` (`12<13`, `16<18`, ...) — so for `n>=8` the bipartite
config has too few free edges to encode all classes (information-theoretically impossible).

## The phase transition at n=6
The bipartite shape **works for `n<=5` and dies at `n=6`** — earlier than the `n>=8` information barrier:
- balanced-bipartite REALIZES for `n=4,5` (verified, for specific "mixing" cut orientations);
- for `n=6` it FAILS for **every** split (`3+3` `f=6`, `2+4` `f=7`, `1+5` `f=10`) and **every** cut
  orientation — the two triangles plus any fixed `K_{3,3}` bridge produce too many isomorphic `T_6`'s to
  cover all 56 classes. The `n=6` achiever (`rho=7`) is a **non-bipartite connected 7-edge** config
  (degree sequence `2,2,2,2,3,3`).
So there are two thresholds: **`n=6`** (the max-cut shape stops realizing, `rho` exceeds `LB`) and **`n=8`**
(the max-cut free-count drops below `LB`). The shape's failure precedes its impossibility.

## Related notions to study (defined here)
- **Transversal-rank**: min subcube that BIJECTS with `G_n` (`2^k=|G_n|`) — only possible when `|G_n|` is a
  power of 2, i.e., `n=3` (`2`) and `n=4` (`4`); never again (`|G_n|` odd for `n>=5`: `12` is even but not a
  power of 2... `|G_5|=12`, so no). So exact "flip-codes" exist only at `n<=4`.
- **Cut-orientation mixing**: which orientations of the fixed `K_{a,b}` cut realize? (Not all do; the
  "totally-A-beats-B" orientation fails — need a mixing cross-tournament.) Characterizing the realizing
  cut orientations is open.
- **Base-path-anchored flip-rank**: the minimum when the fixed arcs must contain a Hamiltonian path (the
  naive model's constraint). This is `>= rho(n)` and `<= C(n-1,2)`.
- **Rainbow number**: the MAX subcube dimension with all `2^k` completions in DISTINCT classes (an
  independent set in the "same-class" graph on the cube) — the dual extremal quantity.
- **The cut/cycle link**: the balanced-cut config fixes a *vertex-bipartition cut*; the base path fixes an
  *edge-cut* (score hierarchy, cut space) with the tiles as cycle space. Whether the two "cut" notions
  align (GF(2) cut/cycle, per CLAUDE.md) is worth checking.

## Net
The flip-rank `rho(n) = 1,2,4,7` for `n=3..6`; the minimal realizing shape is the **balanced max-cut** (fix
the cut, flip the two sides) for `n<=5`, matching the log-2 floor; it undergoes a **phase transition at
n=6** where it fails entirely and `rho` exceeds the floor by 1. The bipartite free-count `f(n)` meets the
floor for `n<=7` and falls below it for `n>=8`. A structural discovery about how iso classes sit in the
tournament cube; the general-`n` shape (post-transition) and `rho(n)` for `n>=7` are open.
