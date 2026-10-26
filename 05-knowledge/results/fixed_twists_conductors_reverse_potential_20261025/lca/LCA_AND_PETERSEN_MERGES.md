# Least common ancestors, the Petersen family, and what merges exactly

2026-10-25, klein. Inheritance survey (two read-only subagents over the navigation layer, the canon,
Session56 and the Petersen/Kuratowski files) plus three exact merges. **PROVED** where marked;
suggestions are marked as such. Files are cited by path; nothing here is a new bridge between
Petersen graphs and tournaments beyond what the cited theorems state.

## 1. Inheritance: Petersen, K5, K3,3

- Kohl's Tait-coloured Collatz graph has no `K5`, `K3,3` or Petersen minor; its Kempe `<a,c>` chains
  have half-length `k = v_2(y+1)`; Althöfer's union graph acquires `K3,3` at `N = 52`, `K5` at `68`,
  a Petersen-family minor at `92` and Petersen at `104`
  ([procgen_kuratowski](../../procgen_kuratowski_20260925_tait_kempe_triples.md), PROVED / FINITE-EXACT).
- The Petersen family is the `Delta-Y` class of `K6`; a linear `K6` has three linked triangle pairs
  iff its Gale tournament is `C3[TT2,TT2,TT2]`
  ([THM-4529](../../../01-canon/theorems/THM-4529-petersen-and-heawood-families-in-paley-coordinates-anti-circulant-all-odd-tournaments.md), PROVED).
- Petersen is the disjoint-arc-support graph (THM-261); at `n = 5` the degree-two Walsh support of `H`
  is `L(K5)`, the complement of Petersen (THM-071); `Omega = K5` is the first obstruction to planar
  `Omega(T)` (THM-519/520); the metagraph `G_7/Z_2` has an explicit `K6` minor and `h >= 12`
  (deathstar S59x/S59y reflections), hence is not linklessly embeddable; `G_5/Z_2` contains no
  Petersen subgraph (petersen_metagraph_s287.out).
- Square-sum side: a magic square of squares exists iff a square-sum `K3,3` with AP shores exists
  (Session57); two hourglasses glue along a spanning tree of `K3,3` (Session58).
- MISTAKE-643 corrected the del Pezzo index and the line rule in the Petersen/quintic lane.

## 2. Inheritance: ancestors

- [Session56](../../session56_rmq_lca_connections_20261025/README.md): RMQ = LCA in a stable
  Cartesian tree; in a tournament's canonical strong-module hierarchy the least common module of two
  vertices determines the quotient where their arc is read; low-digit-first `p`-adic tries turn
  equality depth into LCA, `depth = v_p(x-y)`; the minimum number of same-label hierarchies covering a
  cluster family is the chromatic number of its crossing graph.
- Ports: the canonical parent forest `p(T) = T - x` ([ports_owners_f1 anchor §4](../../ports_owners_f1_20261025/anchor.md));
  an LCA exists only inside one `alpha`-column and is a common induced subtournament.
- Berggren/Stern–Brocot ancestry (THM-3345, THM-3497, THM-4057/4059/4061); Collatz inverse tree
  (collatz_bugs: `f(X,Y) = h(X) - h(LCA)`).

## 3. Exact merges

**M1 (PROVED, [ANTIHYDRA A3](../antihydra/ANTIHYDRA.md)).** For the Antihydra odd map and for the
Syracuse maps, the `k`-cylinder of a valuation sequence is a residue class modulo
`2^{1+v_1+...+v_k}`; so two odd 2-adic integers share their first `k` valuations iff their trie LCA
depth is at least `1 + v_1 + ... + v_k`. Valuation partial sums are LCA depths, and the Antihydra
barrier is a depth condition on the cylinders of `27`. Kohl's Kempe half-length `v_2(y+1)` is the
same object: the LCA depth of `y` and `-1`.

**M2 (PROVED, `check_lca_merges.py`).** The Gale tournament `C3[TT2,TT2,TT2]` of a linked `K6`
(THM-4529 Theorem L) is CLONE_INVOLUTION's canonical hostile; its gap action is the regular `C3`
on three gaps, its clone-fibre polynomial is `1 + t + t^2 + t^3`, and its depth triple is
`(1; 1 mod 3; 3 mod 2)` — the smallest instance of THM-4764's modular depth attached to a member of
the Petersen family. Over `Q` the depth is one (one rotation cycle), so the three clone pairs of the
linked configuration cancel at first order only.

**M3 (PROVED, `check_lca_merges.py`).** The crossing graph of the ten 2-subsets of a 5-set is
`L(K5)`, the complement of Petersen, with chromatic number five; by Session56's hierarchy-cover
theorem the ten pairs need exactly five laminar hierarchies. (The Petersen graph itself is the
*compatibility* graph of the same family: disjoint pairs are laminar.)

## 4. What does not merge (recorded stopping reasons)

- The canonical parent forest's `alpha`-grading cannot carry clique minors: one arc flip changes
  `alpha` by at most one, so contracting `alpha`-columns yields a path (subagent observation, not a
  theorem of the repo; recorded as a suggestion).
- `K5` and `K3,3` are automatic in every tournament's underlying complete graph; the meaningful
  analogue is Kim–Seymour strong minors (decoder_minors, CITED), and pair modules are not strong
  contractions. No Petersen-to-tournament map beyond THM-4529/THM-261/THM-071 is claimed.
- Spanning-tree LCA on `K3,3` would only index Session59's fundamental-cycle equations.
