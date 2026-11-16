---
id: THM-5259
title: "Accelerated Collatz T(n) = n/2 (even), (3n+1)/2 (odd); inverse tree with children 2n and, for n = 2 mod 3, (2n-1)/3. (A) Prune the multiples of 3 and contract the unary nodes: every remaining node has exactly two children, every edge length (original edges) lies in {1, 2, 3, 4}, ray edges have length 2 or 4, branch edges 1, 2, 3, 4, determined by the fork n modulo 27 (2, 20: (2, 2); 11: (2, 4); 8: (4, 3); 17, 26: (4, 1); 5, 14, 23: sterile branch, contracted); no two consecutive forks on a ray are sterile; lengths with the left/right structure recover every label; without lengths the tree is the complete binary tree. (B) Tower laws: rows row(3j) = 2j, row(3j+1) = row(3j+2) = 2j+1; the halving arrow 2m -> m spans 2(floor((2m-1)/3) - floor((m-1)/3)) rows, +2 per +3 in m; the chain 3 <- 6 <- 12 <- ... climbs 2, 4, 8, 16 rows; the odd arrows fall into three families by n mod 6 with first members 1 -> 2, 3 -> 5, 5 -> 8, each next member 4 rows higher and 2 rows longer."
status: >
  PROVED (klein-2026-11-16-monoid; residue analysis modulo 9 and 27; not independently audited). FINITE-EXACT: the tree to
  depth 18 (524,287 nodes) and the tower laws to 10^5 (scripts/collatz_tower_tree.py). The statements are the owner's
  description of 2026-11-16 made exact. CITED: COLLATZ_ROWS.md (path_splice_rows_continuation_20261116), THM-5257, wave 46
  of chanpak_waves_20261113 (the numeric-order Collatz tournament reverses the rays).
session: klein-2026-11-16-monoid
date: 2026-11-16
source: 05-knowledge/results/residues_monoid_tower_tree_20261116/TOWER_TREE.md
scripts:
  - 05-knowledge/results/residues_monoid_tower_tree_20261116/scripts/collatz_tower_tree.py (+ receipts/collatz_tower_tree.out, tree_mod27.out)
related:
  - 05-knowledge/results/path_splice_rows_continuation_20261116/COLLATZ_ROWS.md
  - 01-canon/theorems/THM-5257-tree-tournaments-built-one-node-at-a-time-transitive-iff-depth-first-the-3-cycles-are-ancestor-descendant-interloper-triples-and-the-collatz-breadth-first-counts.md
---
# THM-5259 — the pruned accelerated tree and the tower

Proofs in the packet's TOWER_TREE.md §1–2: along a ray the residues alternate 2, 1 modulo 3; a fork `n ≡ 2 (mod 3)`
has odd child `k = (2n − 1)/3 ≡ 0, 1, 2 (mod 3)` for `n ≡ 5, 2, 8 (mod 9)`; `4n ≡ 5 (mod 9)` iff `n ≡ 8 (mod 9)` and
`n ≡ 5 ⟹ 4n ≡ 2 (mod 9)`, so a ray skips at most one sterile fork; the branch's first fork is `k` or `2k` and may skip
one sterile fork; the sixteen cases modulo 27 give the table. The tower laws are two-line computations with floors.
