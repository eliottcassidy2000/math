---
id: HYP-9677
title: "Bounded-run root-change words. Every odd number greater than 10391 is the value of an O/Z/J root-change word in which every run of spine letters (O, Z) has length at most 5; equivalently, it is the first coordinate of gamma (1,0,0)^T for a product gamma of the 63 matrices M_w = [[g(P), Q+1, -[Q even]], [P+1, 2^n, 0], [[P even], 0, -[n = 0]]], w a binary word of length n <= 5 (THM-5148 (D)). With the KNOWN finite base this implies the Chan-Pak conjecture. A variant, inside the uniquely Hamiltonian grammar, of HYP-9640 (B); part (A) of the pre-audit draft (every odd N > 2647 with words of any length) is WITHDRAWN as a restatement of HYP-9494 (D)."
status: >
  OPEN. FINITE-EXACT evidence (THM-5148 (E); exhaustive closures, pre-registered): holds for all odd numbers up to
  2^19, with exactly 49 exceptions below 10391. For run length at most 4 the exceptions do not stop in range (127
  below 2^17, the largest 130247); for run length at most 3 there are 1355.
  WITHDRAWN: part (A) of the pre-audit draft ('every odd N > 2647 is a value of the family with words of any
  length') is the statement of HYP-9494 (D) (kind-pasteur-2026-11-07-S14), FINITE-EXACT there to 10^9; the draft
  filed it again with evidence to 2^18 (audit M1, M4; MISTAKE-768).
  Relation to canon: HYP-9640 (B) (kind-pasteur-2026-11-09-S18) is a finitely generated statement of the same kind
  with other generators (56 two-port block matrices of order at most 4, every odd v >= 493, FINITE-EXACT to
  2 * 10^5; a 13-matrix subset, v >= 961); HYP-9494 (D) bounds the number of root changes (at most 7 above 2647),
  not the run length. No implication between these statements and this one is known.
  Implication to Chan-Pak: PROVED given the KNOWN fact that every odd number below 80,557 except 7 and 21 is a path
  count.
  No proof strategy is established. The statement is a local-to-global question for an orbit of a finitely
  generated semigroup of integer matrices (as HYP-9640).
  AUDIT: the numbers below 2^16 and 2^17 were reproduced by the independent audit of the packet as a by-product;
  the statement was not audited (05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/audit/AUDIT.md,
  section 6.6).
session: kind-pasteur-2026-11-11-S25
source: 05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/README.md
related:
  - 01-canon/theorems/THM-5148-one-apex-chain-semigroup-matrix-law-and-its-incompleteness-in-range.md
  - 05-knowledge/hypotheses/HYP-9494-root-change-carmichael-hierarchy-two-root-change-coverage-c1-depth-three-adic-factor.md
  - 05-knowledge/hypotheses/HYP-9640-chan-pak-as-a-matrix-semigroup-orbit-and-one-defect-digraph-spectra.md
  - 01-canon/theorems/THM-4904-root-change-word-fibers-carmichael-chamber-parity-law-sqrt-density.md
  - 01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md
---

## Statement

`g(0) = 1`, `g(2v+1) = g(v)`, `g(2v) = g(v) + v + (v mod 2)`. For a word `w` in `{0, 1}^n` (`n >= 0`):
`P(w) = sum of 2^(n-j)` over the positions `j` with `w_j = 1`, `Q(w) = sum of 2^(j-1)` over the positions with
`w_j = 0`,

    M_w = [[ g(P), Q + 1, -[Q even] ], [ P + 1, 2^n, 0 ], [ [P even], 0, -[n = 0] ]].

Examples: the empty word gives `[[1, 1, -1], [1, 1, 0], [1, 0, -1]]`; `w = 0`: `[[1, 2, 0], [1, 2, 0], [1, 0, 0]]`;
`w = 1`: `[[1, 1, -1], [2, 2, 0], [0, 0, 0]]`; `w = 10`: `[[3, 3, -1], [3, 4, 0], [1, 0, 0]]`.

**Conjecture.** Every odd `N > 10391` is the first coordinate of `M_(w_k) ... M_(w_0) (1, 0, 0)^T` for some words
`w_0, ..., w_k` of length at most 5.

## What the vectors are

`M_(w_k) ... M_(w_0) (1, 0, 0)^T = (A, W, d)`: `A` is the number of Hamiltonian paths of the tournament
`Y_0 => Y_1 => ... => Y_k` (ordinal sum of one-apex tournaments with apex words `w_0, ..., w_k`, the arc between
consecutive apexes reversed), `W` is the sum of `H(T - P)` over the directed paths `P` starting at the last apex,
and `d` the same sum with the sign `(-1)^(|P| + 1)`. In the O/Z/J grammar of Sessions 197-200 (THM-4904) the
empty word is the letter `J` (a root change), a letter 1 is `O`, a letter 0 is `Z`, and `(D, E, U) = (A, d, W)`;
the conjecture is about the words whose runs of `O`, `Z` have length at most 5. The strong members of the grammar
are the uniquely Hamiltonian strong tournaments (Session197).

## Evidence and limits

- Exhaustive closures (every value of the sub-language below the bound is generated): THM-5148 (E).
- The miss sets shrink with the run length: 1355, 127, 49 odd numbers below `2^17` for lengths at most 3, 4, 5;
  for length 4 the largest miss in range is 130247, so the range tested does not show whether length 4 suffices
  eventually.
- The number of live states grows faster than the range (57 million for the 63 blocks below `2^19`).
- Nothing here bounds the exceptional set. A proof would have to control the distribution of the first
  coordinates modulo powers of 2 and in short intervals.

## What would refute it

An odd number above 10391 that is not reached. The closures are exhaustive, so a refutation below `2^19` is
excluded; beyond the range the cost of the closure grows like the bound to a power about 1.35.
