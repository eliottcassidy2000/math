---
id: HYP-9375
title: "The Legendre tournament on primes congruent to 3 modulo 4 is quasirandom: with p -> q iff (p/q) = +1 (a tournament by quadratic reciprocity), the induced k-vertex isomorphism classes among the primes of P_3 below X occur with the random-tournament frequencies as X grows. FINITE-EXACT: triangle density 0.25339 on the first 120 primes (exact) and 0.2497 on the first 20,000 (sampled), four-vertex frequencies within 0.002 of (3/8,1/8,1/8,3/8); all four 4-vertex classes occur; the tournament is not regular (scores 47..68 around 59.5 at 120 primes)."
status: >
  OPEN (klein-2026-10-25). The definition is exact (quadratic reciprocity: (p/q) = -(q/p) for p, q = 3 mod 4);
  the evidence is a finite census, receipts in the source packet. Likely provable by large-sieve estimates for
  products of Legendre symbols over primes; not attempted. Recorded with the hostile that the numeral four of
  quadratic reciprocity's residue table has no map to A000568(4) = 4.
session: klein-2026-10-25 (codex/session-tournament-twist-unify)
source: 05-knowledge/results/analogy_reduction_atlas_20261025/QR_TOURNAMENT.md
scripts:
  - 05-knowledge/results/analogy_reduction_atlas_20261025/qr_tournament.py (+ qr-normal.jsonl / qr-optimized.jsonl)
related:
  - 01-canon/theorems/THM-640-the-paley-bridge.md
  - 01-canon/theorems/THM-307-burnside-qr-decomposition.md
---

# HYP-9375 — The Legendre tournament on P_3 is quasirandom

See the source note for the definition, the census table, and the salvaged exact content of the
"reciprocity as a tournament" reading. Cheapest decisive test: a `k = 5` census with a large-sieve
error term, or a biased `k`-class persisting at 10^5 primes.
