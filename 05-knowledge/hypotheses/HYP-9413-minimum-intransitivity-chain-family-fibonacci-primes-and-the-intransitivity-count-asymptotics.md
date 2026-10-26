---
id: HYP-9413
title: "After THM-4803: (a) the strong n-tournaments with the minimum n-2 cyclic triples are exactly the 2^(n-4) chain-reversed transitive tournaments (reverse the arcs of a chain 1 = a_0 < ... < a_r = n in TT_n, interior chain set S in {3..n-2}), pairwise non-isomorphic, and the prime ones are the F(n-4) chains with all gaps <= 2; (b) (C3, score function) determines a tournament iff no prime node of its modular tree is size-balanced; (c) the intransitivity counts satisfy c(n)/A000568(n) -> 1/2 and l(n)/2^C(n,2) = 1/2 - (C(n,2) + n) 2^(-n) + O(n^3 4^(-n)); (d) collect the order-10 prime cycle-index types and predict A000568(11), c(11) without enumerating 903,753,248 classes"
status: >
  OPEN (klein-2026-10-26). (a) FINITE-EXACT for n = 5..9 by construction and canonical forms (chain_family.py; the n = 10
  counts 64 and 8 agree); the chain tournaments are strong with n - 2 cyclic triples by a two-line count (proved in
  FAMILIES.md §E), and a gap >= 3 forces a module (proved), so the open parts are the converse characterizations:
  every extremal strong tournament is a chain tournament (an induction through Moon's strong-subtournament theorem
  looks feasible: the deleted vertex lies on exactly one cyclic triple), and gap <= 2 suffices for primeness.
  (b) FINITE-EXACT n = 3..7 (balanced_check.py; the balanced direction is proved, the no-cancellation direction is not).
  (c) the limit 1/2 for both ratios is PROVED (almost all tournaments are prime; 1/N <= 1/2); the second-order term is a
  heuristic from the two-modules and sink/source events (observed 1/2 - l/2^m = 0.0502 at n = 10 against 0.0537).
  (d) a task, not a claim.
session: klein-2026-10-26 (codex/session-tournament-twist-unify)
source: 05-knowledge/results/a000568_intransitivity_families_20261026/FAMILIES.md (§E) and INTRANSITIVITY.md (§3, §5)
scripts:
  - 05-knowledge/results/a000568_intransitivity_families_20261026/scripts/chain_family.py (+ chain-normal.txt / chain-optimized.txt)
  - 05-knowledge/results/a000568_intransitivity_families_20261026/scripts/balanced_check.py (+ balanced-normal.txt / balanced-optimized.txt)
related:
  - 01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md
  - 01-canon/theorems/THM-462-the-cubic-spectrum-of-tournaments-is-gap-free.md
  - 01-canon/theorems/THM-866-axis-level-completeness.md
---

# HYP-9413 — the minimum-intransitivity chain family and the count asymptotics

(a) Moon's bound `c3 >= n - 2` for strong tournaments is attained by exactly `2^{n-4}` classes
for `5 <= n <= 10`, all of them chain-reversed transitive tournaments; among them exactly `F_{n-4}`
are prime (1, 1, 2, 3, 5, 8), the chains whose consecutive elements differ by at most 2 (compositions
of `n - 5` into parts 1 and 2). The level `c3 = n - 2` of the metagraph axis (THM-866) is therefore
a binary cube `2^{n-4}` on the strong side, with a Fibonacci sub-family of primes.

(b) The score function plus the cyclic-triple hypergraph determine the tournament unless a prime
node of the modular tree is balanced with respect to block sizes; exact for `n <= 7`.

(c) `c(n)/A000568(n)`: 0.49 at `n = 10`, limit `1/2` proved; `l(n)/2^{C(n,2)}`: 0.450 at `n = 10`,
limit `1/2` proved, second-order term conjectured.

(d) The 23 cycle-index types through `n = 9` (`data/prime_cycle_types_n9.tsv`) plus the order-10
types (not collected; `tourn_enum.c --primes` with the hash registry, ~10 minutes on 8 cores)
would give `A000568(11)` and `c(11)` by the species equations; `A000568(11) = 903753248` is the
check.
