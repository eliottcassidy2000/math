---
id: HYP-9375
title: "The Legendre tournament on primes congruent to 3 modulo 4 is quasirandom: with p -> q iff (p/q) = +1 (a tournament by quadratic reciprocity), the induced k-vertex isomorphism classes among the primes of P_3 below X occur with the random-tournament frequencies as X grows. FINITE-EXACT: triangle density 0.25339 on the first 120 primes (exact) and 0.2497 on the first 20,000 (sampled), four-vertex frequencies within 0.002 of (3/8,1/8,1/8,3/8); all four 4-vertex classes occur; the tournament is not regular (scores 47..68 around 59.5 at 120 primes)."
status: >
  PROVED (2026-10-29), using CITED Heath-Brown quadratic large sieve Theorem1
  and the prime number theorem in the fixed progression3 mod4. The signed
  matrix has operator norm O_epsilon(X^(1/2+epsilon)); a direct telescoping
  counting lemma gives every fixed induced class the random frequency.
  See tournament_open_problem_routes_20261029/PROOFS.md P1 for the full
  derivation and normalization. Same-agent audit; no literature-priority claim.
  Earlier finite evidence remains below. The numeral four of quadratic
  reciprocity's residue table has no map to A000568(4)=4.
session: klein-2026-10-25 (codex/session-tournament-twist-unify)
source: 05-knowledge/results/analogy_reduction_atlas_20261025/QR_TOURNAMENT.md
scripts:
  - 05-knowledge/results/analogy_reduction_atlas_20261025/qr_tournament.py (+ qr-normal.jsonl / qr-optimized.jsonl)
related:
  - 01-canon/theorems/THM-640-the-paley-bridge.md
  - 01-canon/theorems/THM-307-burnside-qr-decomposition.md
---

# HYP-9375 — The Legendre tournament on P_3 is quasirandom

**Closed:** [P1, prime-prefix quasirandomness](../results/tournament_open_problem_routes_20261029/PROOFS.md).
The proof executes the original large-sieve suggestion and includes the
81/27 frequency ratio and arbitrary subsets above the square-root scale.
This is an averaged asymptotic statement, not fixed-root port recovery.
See the original source for its definition and historical finite census.
