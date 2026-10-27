---
id: HYP-9413
title: "After THM-4803: (a) the strong n-tournaments with the minimum n-2 cyclic triples are exactly the 2^(n-4) chain-reversed transitive tournaments (reverse the arcs of a chain 1 = a_0 < ... < a_r = n in TT_n, interior chain set S in {3..n-2}), pairwise non-isomorphic, and the prime ones are the F(n-4) chains with all gaps <= 2; (b) (C3, score function) determines a tournament iff no prime node of its modular tree is size-balanced; (c) the intransitivity counts satisfy c(n)/A000568(n) -> 1/2 and l(n)/2^C(n,2) = 1/2 - (C(n,2) + n) 2^(-n) + O(n^3 4^(-n)); (d) collect order-10 cycle indices and identify the remaining order-11 prime inputs"
status: >
  PARTIALLY CLOSED / CORRECTED (klein-2026-10-26; Session67 proof of (b)). (a) FINITE-EXACT for n = 5..9 by construction and canonical forms (chain_family.py; the n = 10
  counts 64 and 8 agree); the chain tournaments are strong with n - 2 cyclic triples by a two-line count (proved in
  FAMILIES.md §E), and a gap >= 3 forces a module (proved), so the open parts are the converse characterizations:
  every extremal strong tournament is a chain tournament (an induction through Moon's strong-subtournament theorem
  looks feasible: the deleted vertex lies on exactly one cyclic triple), and gap <= 2 suffices for primeness.
  (b) PROVED using cited modular-fibre structure: module-summed addressed scores prevent cross-level cancellation; exact fibre 2^b. Independent full labelled checks n<=5 supplement the proof.
  (c) the limit 1/2 for both ratios is PROVED (almost all tournaments are prime; 1/N <= 1/2); the second-order term is a
  heuristic from the two-modules and sink/source events (observed 1/2 - l/2^m = 0.0502 at n = 10 against 0.0537).
  (d) a task, not a claim; order-11 totals also require P11 (and self-converse P11 for C3), not just order-10 cycle indices.
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
node of the modular tree is balanced with respect to block sizes. **PROVED at
every order in Session67, using the cited shared modular tree:** the exact
labelled fibre is `2^b`. For each module M, addressed score sums recover
`D(M)=(sum d(v)-binom(|M|,2))/|M|`; child-minus-parent D separates every
quotient contribution, preventing cancellation across tree levels. The
[proof and bounded checks](../results/session67_metagraph_enumeration_information_20261026/session67-incoming/INCOMING_INTEGRATION.md) close (b), not the other open clauses.

(c) `c(n)/A000568(n)`: 0.49 at `n = 10`, limit `1/2` proved; `l(n)/2^{C(n,2)}`: 0.450 at `n = 10`,
limit `1/2` proved, second-order term conjectured.

(d) The 23 cycle-index types through `n = 9` (`data/prime_cycle_types_n9.tsv`) plus the order-10
types (not collected; `tourn_enum.c --primes` with the hash registry, ~10 minutes on 8 cores)
give the lower-order contributions to `A000568(11)` and `c(11)`. The first
coefficient still has the additive unknown P11, and C3 needs
`(P11+P11,self-converse)/2`. A known A11=903753248 can instead infer P11;
order-10 cycle indices alone do not predict it.

> **Update 2026-10-27 (klein, THM-4814).** Item (d) is resolved in its honest form: the order-10 cycle-index types (nine) were collected; the species identity cannot output A000568(11) because the order-11 primes are their own term, but it determines the non-prime side exactly (19,288,658 non-strong, 87,588,161 strong non-prime classes), and with OEIS A000568(11) = 903,753,248 and A002785(11) = 279,968 this gives A280772(11) = 796,876,429, 259,307 self-converse primes, c(11) = 447,114,461 and l(11) = 16,901,802,400,214,078. See `05-knowledge/results/order10_prime_types_n11_ledger_spectral_blackcore_20261027/README.md` §2.
