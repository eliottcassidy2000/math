---
id: THM-5042
title: "Second, independent proof of THM-5056 (mac-mini, odd-graph Chan-Pak; first in canon), as the odd-graph companion of THM-5040: the dominating-set counts of odd graphs (von Brömssen) on at most n vertices are exactly the odd numbers up to 2^n - 1 except 1, 5 and 37. Oddness costs no vertex: every odd k other than 1, 5, 37 is d(G) for an odd graph on ceil(log2(k + 1)) vertices. With THM-5040 this settles the even/odd split of the dominating-set spectrum: exceptional sets {3, 7, 13} (even) and {1, 5, 37} (odd)"
status: >
  PROVED, with a computer-assisted base and an independent audit pending. This is a second, independent proof of
  THM-5056 (mac-mini-2026-11-09-digraph721, in canon about an hour earlier; found without knowledge of it). The two
  proofs share the cone step and differ in the leading-11 step and the seeds:
  - THM-5056 uses K_2-joins and offset moves from K_4 - e, K_2 u C_4, K_2 v (P_3 u E_2) and K_2 u T_55;
  - this proof uses cones over cones and SR from six census seeds at orders 8 and 9.
  - Base: the exhaustive census of THM-4996 (3) (all graphs on <= 9 vertices; d <= 511 at all orders), plus six
    explicit seed graphs on 8 and 9 vertices, each re-verified by brute force.
  - The inductive step is general.
  - Checks: odd_companion.py recomputes the exact odd-graph sets for n <= 9 from the census, re-verifies the seeds, and
    builds and brute-forces 65 witnesses at n = 10, 11 (odd_companion.out).
session: opus-2026-11-09-S89
date: 2026-11-09
related: [THM-5056, HYP-9625, THM-5040, THM-5055, THM-4996, HYP-9575, THM-4531, HYP-9610]
script: 04-computation/experiments/egcp_20261109/odd_companion.py
---

# THM-5042 - Odd-graph companion (second proof of THM-5056)

> Priority: [THM-5056](THM-5056-odd-graph-chan-pak-dominating-set-counts-of-odd-graphs-are-all-odd-numbers-except-1-5-37-on-ceil-log2-vertices.md) (mac-mini) states and proves the same theorem first; HYP-9625 is the hypothesis it settled. This file records an independent proof, cross-cited, in the same way that THM-5055 is cross-cited as an independent second proof of THM-5040.

Notation as in [THM-5040](THM-5040-sharp-even-graph-chan-pak.md). G is **odd** if some automorphism reverses an odd
number of edges. Let E_odd(n) = {d(G) : G odd with at most n vertices}. Padding keeps oddness: the odd automorphism
extends.

**Theorem.** For every n >= 2, E_odd(n) = {odd k <= 2^n - 1} minus {1, 5, 37}.

**Consequences.**
- {d(G) : G odd} is the set of odd numbers other than 1, 5, 37.
- The least order of an odd graph with d = k is ceil(log2(k + 1)), the same as over all graphs (THM-4996 (1); S.
  Alikhani, Graphs Combin. 29 (2013) 1175-1181, for completeness). Even graphs need one extra vertex on the top bands
  and at X(n) (THM-5040).

## Proof

**Exclusions.**
- d <= 2^n - 1 for every graph on n vertices.
- 1 is realised only by edgeless graphs, and 5 only by P_3 plus isolated vertices. Both are even.
- 37 is prime. Its realisers are a connected graph plus isolated vertices, and the four connected realisers are even.
- These facts come from the exhaustive census of THM-4996 (3), which is complete for d <= 511 at every order by its
  Lemma B3.

**Lemma A (cone over a cone).** For any graph G on n - 2 vertices, K_1 v (K_1 v G) is odd, with
d = 2^(n-1) + 2^(n-2) + d(G).

*Proof.*
- The two apexes are universal and adjacent, so they are adjacent twins. Their transposition reverses only the edge
  between them (labelling the apexes first), so the graph is odd (B2.1).
- The value comes from applying the cone formula twice (THM-4996 Lemma B1). QED

**Lemma B (padded cone).** If G is odd on at most n - 2 vertices, then K_1 v (G u E_j), with j >= 1 chosen so that the
inner graph has n - 1 vertices, is odd.

*Proof.* The inner graph has an isolated vertex, so it has no universal vertex. Apply THM-4996 B2.3: such a cone is even
iff the inner graph is. QED

**Lemma C (odd seeds).** For every m >= 8 and c in {1, 5, 37}, there is an odd graph on m vertices with
d = 2^(m-1) + c that has no universal vertex and has a leaf q fixed by some odd automorphism.

*Proof.*
- **Base, m = 8 and 9.** The seeds below are census graphs, hex vertex pairs, with q the marked leaf:

  | (m, d) | q | edges |
  | --- | --- | --- |
  | (8, 129) | 3 | 0102030406071516172757 |
  | (8, 133) | 3 | 01020304051617262767 |
  | (8, 165) | 5 | 010203041526273637464767 |
  | (9, 257) | 4 | 0102030405071618273868 |
  | (9, 261) | 6 | 0102030405071627384858 |
  | (9, 293) | 6 | 010203040507081617182838485878 |

- **Step.** Apply the leaf stable recursion SR of THM-5040 Lemma 7 at q. It maps d to 3*2^(m-1) + d, so 2^(m-1) + c
  goes to 2^(m+1) + c. The new graph has the new leaf u and no universal vertex.
- An odd automorphism g fixing q extends to SR(G) by fixing the hub and u. Its sign is unchanged: labelling the hub
  first, the hub edges keep their orientation. So SR(G) is odd, and g fixes the leaf u. The property is inherited.
  QED

**Induction.** The base n <= 9 is the census (odd_companion.py checks the formula for 2 <= n <= 9). For n >= 10, the
odd values not already in E_odd(n - 1) are k = 2^(n-1) + x with x odd in [1, 2^(n-1) - 1]:
- **x >= 2^(n-2) + 1.** Lemma A, with G any graph on n - 2 vertices with d(G) = x - 2^(n-2) <= 2^(n-2) - 1 (threshold
  graphs, THM-4996 (1)).
- **x <= 2^(n-2) - 1 and x not in {1, 5, 37}.** Lemma B, with G odd on <= n - 2 vertices (induction hypothesis).
- **x in {1, 5, 37}.** Lemma C. Note that n >= 10 forces 2^(n-2) > 37.

So E_odd(n) contains every odd k <= 2^n - 1 outside {1, 5, 37}; the exclusions give the reverse inclusion. QED

## Remarks

- **The even/odd split.** Over all graphs, d takes every odd value (Alikhani 2013) on ceil(log2(k + 1)) vertices.
  - Splitting by von Brömssen parity, the odd side keeps the minimal orders and loses only {1, 5, 37}.
  - The even side loses {3, 7, 13} and pays one extra vertex on the top bands (M(n - 1), 2^(n - 1)) and at X(n - 1).
  - The extra cost sits on the even side because two universal vertices, or a K_2 or K_3 component, force oddness.
    The cheapest graphs for the top values (cones over cones, K_n) are odd.
- **Tournaments: a typed ANALOGY ([HYP-9610](../../05-knowledge/hypotheses/HYP-9610-even-graph-domination-no-universal-and-dictionary.md)).**
  - Tournaments are equinumerous with even graphs (RPGFD).
  - Their path counts H are odd (Rédei), and the only known exceptions are {7, 21} (Chan-Pak).
  - No map carries H to d (THM-4996 (4), THM-4531). The analogy is between exceptional-set phenomena, not values.
