---
id: THM-5202
title: "Bottleneck factorization: in bounded-span ribbon tournaments the powers of the Moon-Busch floor blocks have rank-one transfer, so H(p u s) = L(p) R(s) for integers L, R whenever a floor run u sits inside the word; prime counts must avoid every bulk floor run, and the cheapest bottleneck-free growth is ((9 + sqrt 113)/2)^(1/4) = 1.7700005 per vertex against 5^(1/3) = 1.7099759 for the floor"
status: >
  (A) PROVED (elementary, probably classical). (B) PROVED (exact integer linear algebra on the minimal ribbon
  representations of THM-4880). (C) PROVED from (B). (D) FINITE-EXACT (periodic words; ranks modulo two primes, exact
  sympy on the named blocks). (E) FINITE-EXACT within the generated near-floor family. The consequence for the
  near-floor conjecture NF2 is HYP-9742 (CONJECTURE). Chan-Pak OPEN.
session: mac-mini-2026-11-15-bottlenecks (wave 28 of chanpak_waves_20261113)
date: 2026-11-15
source: 05-knowledge/results/chanpak_waves_20261113/README.md (wave 28)
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/w28_*.py, w28_*.c, w28_receipts.txt
related:
  - THM-4880 (ribbon transfer, minimal representations), THM-4881 (block 1-.-2 triples), THM-4744 (Busch floor)
  - HYP-9742 (NF2 fails for large primes), HYP-9732 (crossing lemma), THM-5192 (carrier forms)
  - 7-21 repository, notes/05 D9 (NF2) and the concurrent notes/06 E3 (p-adic degeneration), credited below
---

# THM-5202 - bottleneck factorization and the prime premium

Notation. f(n) = min{3^a 5^b : 2a + 3b = n - 1} is the Moon-Busch floor (minimum H of a strong n-tournament). A ribbon
word of span w lists letters a_q ⊆ {1..w}; vertex q has back arcs q → q - s for s in a_q, and all other arcs point
forward. H(word) = fin · M_{a_n} ⋯ M_{a_1} · init in the minimal integer representation of THM-4880: dimension 14 for
w = 2, 62 for w = 3. Letters are written as digit strings ("0" = empty, "12" = {1,2}).

## (A) Coalescence inequality (PROVED)

Let T = A ∪ B with A ∩ B = {v}, and suppose every arc between A - v and B - v points from A - v to B - v.

- *Unique decomposition.* Every Hamiltonian path of T is α β v γ ζ, uniquely, with α, γ inside A - v and β, ζ inside
  B - v (any of them may be empty). A path can pass from B - v back to A - v only through v.
- *Injection.* The map (α v γ, β v ζ) ↦ α β v γ ζ embeds pairs of Hamiltonian paths of A and B. Hence

      H(T) = H(A) H(B) + Extra,   Extra ≥ 0,

  where Extra counts the paths with (α ≠ ∅, β ≠ ∅, v → last α) or (γ ≠ ∅, ζ ≠ ∅, first ζ → v).
- *Exact cases.* Two 3-cycles glued at a vertex give 9 = 3·3. A 3-cycle followed by the strong 4-tournament in
  ribbon form 0.0.2.0.2.2 gives 17 = 3·5 + 2.

## (B) Floor powers have rank one (PROVED, exact)

- *Span 2.* The 3-unit 1.2 has characteristic polynomial x^13 (x - 3). Its powers have ranks 4, 3, 2, 1, 1, …, and
  M_{1.2}^4 = c r^T with

      c = (15, 15, 12, 14, 18, 32, 20, 12, -24, 24, 12, 0, 0, 0),
      r = (1, 2, 3, 0, 0, 0, 0, 0, 2, 2, 0, 2, 1, 0).

  The unit 0.2 also has characteristic polynomial x^13 (x - 3).
- *Span 3.* The 5-units 0.1.3, 0.3.1 and 0.1.123 have power ranks 3, 2, 1, so their cubes have rank one. The units
  1.2.3, 1.12.3 and 0.12.3 reach rank one at the fourth or third power. In each case the trace, and the eigenvalue, is 5.
- *Span 2 is irreducible on single letters* (THM-4880 (G)). The only rank-one factors of length ≤ 3 there are
  non-strong cuts such as 0.0 and 0.1. The bottlenecks above come from *powers* of strong floor units, whose
  remaining spectrum is nilpotent.

## (C) Factorization through a floor run (PROVED)

If rank M_u = 1, write M_u = c r^T with c a primitive integer vector; then r is integral. For all words p and s,

    H(p u s) = L(p) · R(s),   L(p) = r · y(p),   R(s) = fin · M_s c,

both integers. For example H(p (1.2)^j s) = 3^(j-4) L(p) R(s) for j ≥ 4. A replay on 300 random span-2 words matched
300 of 300; one instance is 2,130,237 = 753 × 2,829.

**Corollary (prime avoidance).**
- If H(p u s) is prime, then |L(p)| = 1 or |R(s)| = 1.
- In a strong word, L and R grow with the lengths of p and s.
- So a ribbon word with prime count contains no floor run, of the lengths in (B), outside a bounded distance from its
  two ends.
- More generally, every bulk floor run contributes an exact factor 3^j or 5^j times a factorization L·R. This is the
  mechanism behind the concurrent session's observation (7-21 notes/06 E3) that the 3- and 5-adic valuations of
  near-floor counts grow linearly with the order.

## (D) The premium of bottleneck-free words (FINITE-EXACT)

The test is on periodic words whose bulk is strongly connected (the middle third lies in one strong component). Among
them, every block of growth below 1.7700 per letter is a bottleneck: rank M^dim = 1 modulo the primes 33,554,393 and
33,554,383. The families checked are:
- span 2, all periods ≤ 10: 591 classes;
- span 3, all periods ≤ 6: 408 classes;
- span 3, B^j s with B a floor 5-unit, j ≤ 2 and |s| ≤ 4: 1,926 words.

The cheapest bottleneck-free blocks are 0.2.1.2 and 0.12.2.2 (and the span-3 analogue 0.13.12.3). Their characteristic
polynomial is x^12 (x² - 9x - 8), exactly. So

    ρ* = ((9 + √113)/2)^(1/4) = 1.7700005,   ρ*/5^(1/3) = 1.035104.

The second root -0.81507 is what lets these counts avoid exact factorization. The lowest-growth strong blocks of all
are bottlenecks with integer eigenvalues: 3 per 2 letters, 5 per 3, 243 per 10, 85 per 8, 17 per 5, and so on.

## (E) Near-floor census (FINITE-EXACT within the generated family)

The program `w28_nearfloor_ribbon_dfs.c` runs a DFS over strong ribbon words with the prefix prune
H(prefix of length j) < 4 f(j). This prune is a heuristic restriction: it is complete only if near-floor words have
near-floor prefixes.

- *Span 3 to order 14* (1.09 × 10^9 words at order 14). The windows [5f/3, 3f) miss only:
  - 63 at order 7;
  - 637 at order 12;
  - 1067, 1109 and 1339 at order 13;
  - 1879, 1903, 1943, 1949, 1993, 2051 and 2089 at order 14.

  The first three lines agree with the 7-21 session's independent degree-1 generator.
- *Every miss is prime or has a prime factor ≥ 13*: 7²·13, 11·97, 13·103, 11·173, 29·67, 7·293.
- *Smallest prime count over f(k)*: 1.178, 1.187, 1.336, 1.249, 1.389, 1.410, 1.423 at k = 8, …, 14. Span 2 and
  span 3 agree there. Span 2 continues 1.492, 1.511, 1.536, 1.713 at k = 15, …, 18.
- *Growth of that ratio*: 3.2% per vertex from k = 10 to k = 18, against the predicted premium of 3.51% from (D).

## Credit and scope

- (A) is the elementary series-composition argument and is surely classical in some form.
- The concurrent 7-21 session (notes/06, in progress during this session) found by statistics that NF2(c) degenerates
  for every fixed c. This theorem supplies the exact algebraic reason and the explicit constant.
- Nothing here bounds non-ribbon tournaments. The 7-21 census shows the near-floor classes have span ≤ 4.
