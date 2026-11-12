---
id: HYP-9677
title: "Block semigroup conjecture for the Chan-Pak problem. For a binary word w of length n let P be the value of its letters 1 (first letter = top bit), Q the value of its letters 0 (first letter = unit bit), g the one-apex function, and M_w = [[g(P), Q+1, -[Q even]], [P+1, 2^n, 0], [[P even], 0, -[n = 0]]]. The first coordinates of the vectors M_(w_k)...M_(w_0)(1,0,0)^T are the numbers of Hamiltonian paths of the apex-path tournaments (iterated joins of one-apex tournaments at their apexes; THM-5148). (A) Every odd number greater than 2647 is such a first coordinate. (B) Already with the 63 words of length at most 5: every odd number greater than 10391 is the first coordinate of gamma (1,0,0)^T for a product gamma of the matrices M_w, |w| <= 5. Either statement implies the Chan-Pak conjecture."
status: >
  OPEN. FINITE-EXACT evidence (THM-5148 (D), (E); exhaustive closures, pre-registered): (A) holds for all odd
  numbers up to 2^18, with exactly 29 exceptions below (7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247,
  251, 277, 299, 421, 449, 503, 533, 569, 583, 637, 821, 1067, 1339, 1861, 2647); (B) holds for all odd numbers up
  to 2^17, with exactly 49 exceptions below, the largest 10391. For words of length at most 4 the exceptions do
  not stop in range (127 below 2^17, the largest 130247); for length at most 3 there are 1355.
  Implication to Chan-Pak: PROVED given the KNOWN fact that every odd number below 80,557 except 7 and 21 is a path
  count (the exceptions of (A) and (B) other than 7 and 21 are counts of other tournaments).
  No proof strategy is established. The statement is a local-to-global question for an orbit of a finitely
  generated semigroup of integer matrices. Related families in canon with finite completeness: bank x one-apex
  joins to 2^35 (THM-4885 (D)), two-apex words below 2^23 (opus-2026-11-11-S93, CONCURRENT), strong span-2 words
  (THM-4880), two-core language (Session270); none is proved complete.
  AUDIT: pending at this commit (PRE-AUDIT DRAFT).
session: kind-pasteur-2026-11-11-S25
source: 05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/README.md
related:
  - 01-canon/theorems/THM-5148-one-apex-chain-semigroup-matrix-law-and-its-incompleteness-in-range.md
  - 01-canon/theorems/THM-4866-single-cycle-floor-bottom-band-uniquely-hamiltonian-one-arc-lorentzian-joins-and-crevasses.md
  - 01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md
---

## Statement

`g(0) = 1`, `g(2v+1) = g(v)`, `g(2v) = g(v) + v + (v mod 2)`. For a word `w` in `{0, 1}^n` (`n >= 0`):
`P(w) = sum of 2^(n-j)` over the positions `j` with `w_j = 1`, `Q(w) = sum of 2^(j-1)` over the positions with
`w_j = 0`,

    M_w = [[ g(P), Q + 1, -[Q even] ], [ P + 1, 2^n, 0 ], [ [P even], 0, -[n = 0] ]].

Examples: the empty word gives `[[1, 1, -1], [1, 1, 0], [1, 0, -1]]`; `w = 0`: `[[1, 2, 0], [1, 2, 0], [1, 0, 0]]`;
`w = 1`: `[[1, 1, -1], [2, 2, 0], [0, 0, 0]]`; `w = 10`: `[[3, 3, -1], [3, 4, 0], [1, 0, 0]]`.

**(A)** Every odd `N > 2647` is the first coordinate of `M_(w_k) ... M_(w_0) (1, 0, 0)^T` for some words
`w_0, ..., w_k`.

**(B)** Every odd `N > 10391` is such a first coordinate with all `w_j` of length at most 5.

## What the vectors are

`M_(w_k) ... M_(w_0) (1, 0, 0)^T = (A, W, d)`: `A` is the number of Hamiltonian paths of the tournament
`Y_0 => Y_1 => ... => Y_k` (ordinal sum of one-apex tournaments with apex words `w_0, ..., w_k`, the arc between
consecutive apexes reversed), `W` is the sum of `H(T - P)` over the directed paths `P` starting at the last apex,
and `d` the same sum with the sign `(-1)^(|P| + 1)` (THM-5148 (D)). One has `0 <= d <= A <= W` and every step
satisfies `A' >= W`.

## Evidence and limits

- Exhaustive closures (every state with `A` below the bound is generated): THM-5148 (D), (E).
- The miss sets shrink with the word length: 1355, 127, 49 odd numbers below `2^17` for lengths at most 3, 4, 5;
  for length 4 the largest miss in range is 130247, so the range tested does not show whether length 4 suffices
  eventually.
- The number of live states grows faster than the range (8.2 million for the 63 blocks below `2^17`), so the
  multiplicity grows; a fixed finite bank does not have this property (THM-4885).
- Nothing here bounds the exceptional set. A proof would have to control the distribution of the first
  coordinates modulo powers of 2 and in short intervals.

## What would refute it

An odd number above the stated bounds that is not reached. The closures are exhaustive, so a refutation in range
is excluded; beyond the range the cost of the closure grows like the bound to a power about 1.35.
