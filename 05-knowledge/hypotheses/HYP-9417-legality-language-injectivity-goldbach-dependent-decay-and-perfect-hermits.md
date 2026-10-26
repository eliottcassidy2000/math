---
id: HYP-9417
title: "After THM-4806: (1) the Goldbach-dependent share of the odd aliquot values in the decade [10^d, 10^(d+1)) is c / ln(5 x 10^d) with c in [1.3, 1.6] (pre-registered 0.0812 at d = 7 CONFIRMED), so almost all odd values have a non-semiprime source (density-zero dependence); next decade predicted in [0.070, 0.078], the 2 x pq pattern share in [0.185, 0.200], the 2 x prime class share in [0.112, 0.125]; (2) the seventh even perfect number 137438691328 is not a hermit: a companion with at least five odd prime factors exists; (3) the PSD leakage rank of boundary_support, recomputed with named guard operators, reproduces 2^k - d_t on the session70 fibre. (The all-order injectivity of the named legality language, listed here before the audit, is PROVED: session70, THM-4807, THM-4806 R'; MISTAKE-677.)"
status: "NUMERICAL / OPEN: (1) four of five pre-registered 10^8 predictions confirmed (decade share 0.0812 in [0.074, 0.083]; count = S57 + 1; 2 x prime class 0.1278 in [0.125, 0.145]; enrichment 1.55 in [1.45, 1.75]) and one refuted (2 x pq pattern share 0.1990 against [0.200, 0.240]); the 10^9 predictions are untested; (2) OPEN (no companion with at most two odd prime factors, and none among the enumerated larger shapes); (3) untested. The former item (1), all-order injectivity, was already proved and audited on main (session70's guards lane; THM-4807) when this file was first written; it is removed (MISTAKE-677)."
session: opus-2026-10-26-S61
source: 05-knowledge/results/legality_language_goldbach_free_perfect_companions_20261026.md
related:
  - 01-canon/theorems/THM-4806-named-legality-language-single-pair-lemma-goldbach-free-shapes-and-row-companions-of-even-perfect-numbers.md
  - 01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md
  - 05-knowledge/hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md
  - 05-knowledge/hypotheses/HYP-9415-aliquot-table-predictions-twin-value-sources-self-only-perfect-numbers-mersenne-abc-hits.md
  - 05-knowledge/results/session70_attachment_guard_alignment_20261026/README.md
---

# HYP-9417: after THM-4806

1. **Goldbach-free coverage in density.** The dependent share of the decade `[10^d, 10^(d+1))` is `c / ln(5 x 10^d)` with `c` tending to a constant in `[1.3, 1.6]` (`1.40, 1.39, 1.44` for `d = 5, 6, 7`), so the Goldbach-dependent odd values have density zero: for almost every odd `v`, coverage does not depend on the Goldbach antecedent. Pre-registered for `10^9` (untested): decade share in `[0.070, 0.078]`; the `2 x pq` share of the dependent set in `[0.185, 0.200]`; the dependent share within `v - 1 = 2 x prime` in `[0.112, 0.125]`.
2. **The seventh perfect number.** `137438691328 = 2^18 (2^19 - 1)` has no row companion (`q_a` composite for every proper divisor `a` of 18) and no companion among the enumerated shapes (THM-4806 P, scope in its status); conjecture: a companion with five or more odd prime factors exists, so 28 remains the only perfect hermit (HYP-9415 (4)).
3. **Leakage rank.** The PSD leakage identity of the strong-atom HISTORY (`I - K^2 = [(I - E)UE]^*[(I - E)UE]`) and the boundary_support Prop. 0 rank, evaluated with the session70 guard operators instead of the complement, give the missing rank `2^k - d_t`.

Removed after the audit: the all-order injectivity of `T -> L(T)` (proved and audited on main before this file was written: session70's guards lane; THM-4807; second proof THM-4806 R').
