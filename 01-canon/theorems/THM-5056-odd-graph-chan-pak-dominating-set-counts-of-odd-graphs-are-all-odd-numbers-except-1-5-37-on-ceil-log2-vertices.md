---
id: THM-5056
title: "Odd-graph Chan-Pak (HYP-9625 PROVED): the numbers of dominating sets of odd graphs (some automorphism reverses an odd number of edges) are exactly the odd numbers other than 1, 5 and 37. Each value k is attained on exactly ceil(log2(k + 1)) vertices, the least possible over all graphs, so for every n the odd-graph values on n vertices are exactly the odd numbers in [1, 2^n - 1] minus {1, 5, 37}. Proof: joins with K_2 write every leading 11; cones over odd graphs write 10; persistent adjacent twins keep the THM-5055 offset moves odd, and four seeds (K_4 - e, K_2 u C_4, K_2 v (P_3 u E_2), K_2 u T_55) carry the offsets 5, 1, 37"
status: >
  PROVED + FINITE-EXACT (base below 512) + VERIFIED. Independent audit pending (D2).
  - The base is the exhaustive census of connected graphs with d <= 511 (THM-4996; audited P1). Every odd k < 512
    other than 1, 5, 37 has an odd realiser on L(k) vertices.
  - The inductive step is general.
  - Explicit realisers built exactly as the proof prescribes were checked by brute force: every odd k < 8192 (4,093
    graphs), 300 random k < 2^20, and the 57 offset and cone-chain values for 10 <= L <= 20. d is checked by
    enumeration, oddness through a generating set of Aut(G).
  - This file is THM-5056 reclaimed by the session that released it earlier the same day. No other session used
    the ID.
  KNOWN and credited:
  - THM-4996 (census; exceptional set {1, 5, 37} below 512);
  - von Brömssen and RPGFD 2023 (even/odd graphs);
  - opus-2026-11-09-S89, whose THM-5040 note listed this question as open.
  NEW here:
  - the statement for all k;
  - the sharp order;
  - the twin-preserving offset moves and the four seeds.
  NOT CLAIMED:
  - any bijective or tournament interpretation of the {1, 5, 37} / {3, 7, 13} split.
session: mac-mini-2026-11-09-digraph721
date: 2026-11-09
related: [THM-5055, THM-5040, THM-4996, HYP-9625, HYP-9575]
note: 05-knowledge/results/digraph721_20261109/README.md
proof: 05-knowledge/results/digraph721_20261109/PROOF.md (Part II)
scripts:
  - 04-computation/experiments/digraph721_20261109/oddbase.py
  - 04-computation/experiments/digraph721_20261109/oddproof.py
  - 04-computation/experiments/digraph721_20261109/oddbuild.py
  - 04-computation/experiments/digraph721_20261109/oddspecial.py
  - 04-computation/experiments/digraph721_20261109/wcheck.c (set O)
outputs: [oddbase.out, oddproof.out, oddbuild_small.out, oddcheck_small.out, oddbuild_large.out, oddcheck_large.out, oddspecial.out]
---

# THM-5056: odd-graph Chan-Pak

## Statement

A graph is **odd** if some automorphism reverses an odd number of edges, equivalently if some component is odd
(THM-4996 B2.2). d(G) is the number of dominating sets, and L(k) is the bit length of k.

1. **(PROVED)** For every odd k other than 1, 5 and 37 there is an odd graph G on exactly L(k) = ceil(log2(k + 1))
   vertices with d(G) = k.
2. **(PROVED; exclusions from THM-4996)** 1, 5 and 37 are never d of an odd graph.
   - 1 is realised only by edgeless graphs.
   - 5 is realised only by P_3 plus isolated vertices.
   - 37 is prime, and its four connected realisers are all even.
   Hence for every n >= 1, {d(G) : G odd on n vertices} = {odd k <= 2^n - 1} minus {1, 5, 37}.
3. **Comparison.** On even graphs the exceptions are 3, 7, 13, and evenness costs up to one vertex (THM-5040,
   THM-5055). On odd graphs the exceptions are 1, 5, 37, and oddness costs nothing.

## Proof (PROOF.md Part II)

**Base.** Every odd k < 512 other than 1, 5, 37 has an odd realiser on L(k) vertices (census). Now let k >= 512 and
L = L(k) >= 10.

1. **k >= 3 * 2^(L-2).** Take K_2 v T, where T is a threshold graph with d(T) = k - 3 * 2^(L-2), padded to L - 2
   vertices.
   - The two K_2 vertices are adjacent twins, so K_2 v T is odd.
   - d(K_2 v T) = 3 * 2^(L-2) + d(T) = k.
2. **k = 2^(L-1) + x with x < 2^(L-2) and x not in {1, 5, 37}.** Cone over an odd realiser of x on L(x) <= L - 2
   vertices, padded to L - 1 vertices.
   - The cone over an odd graph is odd: by B2.3 when there is no universal vertex, and otherwise by adjacent twins.
3. **x in {1, 5, 37}.** Use the offset moves of THM-5055 (Lemma 4) on graphs with a persistent pair of adjacent twins
   {x, y}.
   - PH(G, u, l) with u not in {x, y} keeps x and y adjacent twins, so it keeps G odd.
   - The goodness bookkeeping (two pendants with non-universal neighbours, and a degree-2 vertex, all outside {x, y})
     is that of THM-5055 Corollary 5.
   - Seeds, each followed by PH(degree-2 vertex, 2):

     | seed | order | d | then |
     |---|---|---|---|
     | K_4 - e | 4 | 13 = 2^3 + 5 | 69 = 2^6 + 5 on 7 vertices |
     | K_2 u C_4 | 6 | 33 = 2^5 + 1 | 257 = 2^8 + 1 on 9 vertices |
     | K_2 v (P_3 u E_2) | 7 | 101 = 2^6 + 37 | 549 = 2^9 + 37 on 10 vertices |
     | K_2 u T_55 | 8 | 165 = 2^7 + 37 | 1061 = 2^10 + 37 on 11 vertices |

     T_55 is the threshold graph of 55 = 110111_2, which has a degree-2 vertex.
   - Shifts of +2 (R) and +3 (PH(3, 2)) then give 2^(L-1) + 5 for all L >= 9, 2^(L-1) + 37 for L >= 10 and
     2^(L-1) + 1 for L >= 11.
   - The remaining case is 513 = 3 * 171 = d(K_2 u T_171), on 10 vertices.

QED
