# The blue/black flip-lines are an EVEN GRAPH + a T-JOIN — the grid reflection is the complement, and SC/NS tiling-parity is the T-join boundary

*opus-2026-07-01-S18. The owner asked to describe and constrain, as precisely as possible, the process by which
blue/black flip-lines (tournament-tiling-explorer) are assigned between pure-blue / pure-black / mixed merged
nodes so each node gets its correct tiling count — with new metrics, reframings, recursion, and corrections to
the conjecture.*

> **CONVERGENCE (concurrent):** mac-mini-S83 (HYP-3808) independently derived the identical parity
> decomposition — black = even graph, blue = T-join/all-odd on SC, τ-parity = SC/NS, the same closed forms
> (#blue = 2^{(m+⌊(n−1)/2⌋)/2−1}), the same eligibility rules, and the same n=6 "self-loops-only-on-mixed"
> correction. This reflection agrees with all of it. **My unique contribution is the GENERATING MECHANISM:
> the grid reflection R = the tournament complement** — which *proves* the parities (R pairs black lines →
> even; grid-sym = R-fixed = SC) and *characterizes the self-loops* (black self-loop ⟺ R(t)=flip(t)), plus
> the Z/2-chain reframe (boundary = SC) and the realization-degeneracy metric. I record it as an extension of
> HYP-3808, not a duplicate. The process resolves into two objects (even graph + T-join) with one generator
> (complement).

## The setup, made precise
A **line** = a flip-pair `{t, flip(t)}` where `flip` reverses all `m = C(n-1,2)` tiles (the d=m waggly layer).
Each tiling is an endpoint of exactly one line, so the merged-metagraph pairing identity (THM-346 bucket
balance, kind-pasteur) reads, per node `v`:
```
    tau(v) = 2 * (self-loop lines at v) + (edge lines incident to v),    sum_v tau(v) = 2^m.
```
A line is **blue** iff its tiling is grid-symmetric (`isGridSym`: invariant under (x,y)→(n−y+1,n−x+1)), **black**
otherwise. Nodes are **pure-blue** (all tilings grid-sym), **pure-black** (none), or **mixed** (both). Verified
`tau = 2*self + edges` universally.

## The generating fact: the grid reflection R IS the complement
Let `R` = the grid reflection on tiles. **`class(R(t)) = complement(class(t))` for every tiling** (verified
exactly: 8/8, 64/64, 1024/1024 at n=4,5,6). Two immediate consequences:
- **grid-symmetric ⟺ complement-fixed ⟺ SC.** A grid-sym tiling has `t = R(t)`, so `class(t) =
  complement(class(t))`, i.e. self-complementary. Hence **pure-blue and mixed nodes are exactly the SC nodes;
  pure-black are exactly the NS nodes.** Blue lines live entirely in the SC world; black lines touch NS.
- **R fixes merged nodes** (`mn(R(t)) = mn(t)`, since the merged node is complement-invariant) and **commutes
  with flip**. So R permutes the lines, preserving color and the merged-node pair each line joins.

## The two objects the process decomposes into
### BLACK = an even graph (Eulerian)
For a black line `t` (non-grid-sym, `R(t) ≠ t`), R sends it to a *different* black line `{R(t), flip(R(t))}`
between the *same* merged pair — **unless** `R(t) = flip(t)`, in which case the line is R-fixed. But `R(t) =
flip(t)` means `class(flip(t)) = class(R(t)) = complement(class(t))`, so `flip(t)` lands in `t`'s merged node —
i.e. the line is a **self-loop**. Therefore:
- **Black EDGES pair up under R ⇒ every node has EVEN black-edge-degree ⇒ the black edge graph is an EVEN
  GRAPH.** (Verified: connected, cycle-rank E−V+1 = 1, 14, 425 at n=4,5,6, non-bipartite for n≥5.)
- **Black SELF-LOOPS are exactly the R-fixed black lines** — the tilings whose grid-reflection equals their
  tile-complement (`R(t)=flip(t)`). This is *why* self-loops don't spoil the even-degree count.

### BLUE = a T-join with T = the SC nodes
Grid-sym tilings are R-fixed, so they are not paired away; each SC class carries an odd grid-symmetric core.
Result: **every SC node has ODD blue-edge-degree, every NS node has zero — the blue edge graph is a T-JOIN with
T = all SC nodes.** For a T-join to exist, `|T| = |SC|` must be even — **verified even at every n** (2,8,12,88
for n=4..7). Blue self-loops occur only on **mixed** nodes (the owner's "only mixed" conjecture is correct here).

### The τ-parity law (a corollary)
Since `tau ≡ blue-deg + black-deg (mod 2)`, black-deg is even, and blue-deg is `[v is SC]`:
**`tau(v)` is ODD for SC nodes and EVEN for NS nodes.** (NS-merged glues C and Cᵒᵖ, `tau = 2·tilings(C)`; SC is
a single class with an odd core.) This is the T-join boundary condition wearing a tiling-count hat.

## The eligibility rules (the process, formalized — and the conjecture corrected)
- **Blue line** = an edge between two SC nodes with ≥1 **mixed** endpoint (mixed–mixed or mixed–pure-blue;
  **never** pure-blue–pure-blue), OR a blue self-loop **only on a mixed** node. Blue never touches pure-black.
- **Black line** = an edge between two nodes of {mixed, pure-black} (**any** of mixed–mixed, mixed–pure-black,
  **pure-black–pure-black**), OR a black self-loop on a **mixed or pure-black** node. Black never touches pure-blue.
- **Corrections to the conjecture** (both first appear at n=6): black self-loops **do** occur on pure-black
  nodes, and **pure-black–pure-black** black edges exist (290 of them at n=6). So "self-loops only on mixed"
  and "pure-black connects only to mixed" hold for BLUE but **break for BLACK** — the black even-graph is a full
  even graph on {mixed ∪ pure-black}, not a pure-black↔mixed bipartite one.

## Metrics (the illuminating quantities)
| n | (B, M, K) | SC=B+M | line total 2^{m−1} | blue 2^{(m+f)/2−1} | black | black cycle-rank |
|---|---|---|---|---|---|---|
| 4 | (1,1,1) | 2 | 4 | 2 | 2 | 1 |
| 5 | (3,5,2) | 8 | 32 | 8 | 24 | 14 |
| 6 | (2,10,22) | 12 | 512 | 32 | 480 | 425 |
| 7 | (4,84,184) | 88 | 16384 | 128 | 16256 | — |

`f = ⌊(n−1)/2⌋` = # anti-diagonal (R-fixed) tiles; blue lines `= 2^{(m+f)/2 − 1}`. Further metrics worth
tracking: **the interface size M(n)** (mixed nodes — the *only* nodes in both objects, coupling the SC/blue and
NS/black worlds); the **black even-graph girth / odd-cycle spectrum** (non-bipartite ⇒ odd cycles, tying to
perfection/odd-holes); the **blue T-join excess** over a perfect matching on SC (`blue-edges − |SC|/2`); the
**self-loop split** (blue-self on M, black-self on M∪K = the R-fixed core `R(t)=flip(t)`).

## Reframings (to make the process even more precise)
1. **Z/2 chain view.** The black graph is an element of the **cycle space** (Eulerian, chargeless); the blue
   graph is a chain with **boundary = the SC nodes**. The tiling-count parity is literally that boundary. So
   the whole flip-line assignment is a Z/2 chain realizing boundary = SC, split by the R-involution into an
   R-symmetric part (black, even) and an R-fixed part (blue, T-join).
2. **The realization/degeneracy question.** The parities + category-support are *necessary* conditions. Are they
   *sufficient* to pin the exact edges, or is the true flip-line structure one of many valid realizations? The
   count of valid realizations (vs the true one) is a new **degeneracy metric** measuring how much the parity
   structure alone determines the metagraph — a concrete next computation.
3. **The complement as the single generator.** Everything above flows from `R = complement`: the SC/NS split,
   the even/odd of black/blue, the self-loop cores. The process is "one involution (complement) acting on the
   flip-pairing," not three independent rules.

## Recursion and the concrete next target
- **Recursion.** `K(n)` = NS-merged = 1,2,22,184 (= (A000568−SC)/2); `SC = B+M` = 2,8,12,88; the black
  even-graph is connected with cycle-rank `2^{m−1} − blue − (nodes−1)` growing quadratically. Open: a clean
  Mode-A (vertex-insertion) recursion for the black even-graph and the blue T-join separately — does inserting
  a vertex *suspend* the previous even-graph (cone) or glue along the interface M?
- **Toward the flip-rank exhaustion.** The τ-parity (SC classes odd tiling-count, NS even) is a rigid parity
  invariant of the tiling distribution, and the flip-rank obstruction (the Paley heptagon) is an **SC node**
  (odd tiling-count, in the blue T-join). This suggests a **parity/T-join pruning** for the covering search:
  any covering sub-structure must respect the SC-odd/NS-even parities and the even-graph structure of the
  black layer — a lever the raw arc-hypercube search ignores. Making that pruning concrete (does the T-join
  boundary obstruct low-dimensional covers of the SC classes?) is the bridge back to k(7).

## Status
- **Established (verified n≤6, n=7 category counts):** R = complement; grid-sym ⟺ SC; BLACK = even graph;
  BLUE = T-join (T=SC, |SC| even); τ-parity (SC odd / NS even); the eligibility rules.
- **Corrected:** black self-loops occur on pure-black, and pure-black–pure-black black edges exist (n≥6) — the
  "only mixed" intuition holds only for blue.
- **Proposed (new):** the Z/2 chain reframe (boundary = SC); the realization-degeneracy metric; the black
  even-graph odd-cycle spectrum; a Mode-A recursion; T-join/parity pruning as a bridge to the flip-rank.
- **Honest scope:** structural (verified small n + the R=complement generator, which proves the parities in
  general); the recursion and the flip-rank bridge are proposed, not yet proved.

Related: THM-346 (bucket balance / half-line conservation — the base identity), HYP-3805 (flip-rank; the Paley
obstruction is an SC/T-join node), the even-graph first-class mandate (E_n), the spine/ribs/sea SC-NS
decomposition, MISTAKE-033 (flip tiling ≠ Tᵒᵖ — but R = the true complement). HYP-3808 (this). Scripts:
04-computation/mmg_{blueblack_parity,grid_reflection_is_complement,evengraph_tjoin_decomp,category_recursion_n7}_opus_20260701.py.
