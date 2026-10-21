---
id: HYP-9358
title: "Open laws for Camion pseudoprimes (THM-4746): (a) the Camion-Wieferich primes (s(p) = 2^((p-1)/2) mod p^2, equivalently p^2 a Camion pseudoprime) are infinite in number but of density zero, like Wieferich primes (only 3 and 13 below 3000); (b) the sporadic zero classes {r >= 4 : p | s(r) - c(r)} are Poisson with mean (p-4)/p; (c) the Camion pseudoprimes up to x are 5x/(6 log x)(1 + o(1)), all but o(x/log x) of them of the form 2p or 3p; (d) modulo 9, s(m) = 2^((m-1)/2) for every odd m not = 1 mod 6 (or m = 1), and fails for every m = 1 mod 6 with m >= 7; (e) the quasi-periodicity of THM-4746 holds for every class of labelled structures whose components are linearly ordered and whose total count is a finite combination of exponentials modulo p"
status: >
  OPEN. FINITE-EXACT evidence: (a) p <= 3000 (camion_list_20000.out) and the agreement of the two forms for
  p <= 1500 (camion_verify.out); (b) the 75 primes 7 <= p <= 400: observed 31, 25, 13, 4, 2, 0 primes with
  0..5 sporadic classes against the Poisson prediction 29.1, 27.5, 13.0, 4.1, 1.0, 0.2 (camion_sporadic.out);
  (c) up to 20000: 2145 squarefree Camion pseudoprimes, of which 857 are 3p and 1228 are 2p, leaving 60 exotic
  ones (37 odd, 23 even); (d) odd m <= 700; (e) proved for tournaments (THM-4746) and strong semicomplete digraphs
  (THM-4746 (7)); a counterexample would be a linearly ordered decomposition whose total is not quasi-periodic.
session: mac-mini-2026-10-21-fanopsp
source: 05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md
related:
  - 01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md
scripts:
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_sporadic.py (+ .out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_list.c (+ camion_list_20000.out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_verify.c (+ camion_verify.out)
---

# HYP-9358 — open laws for Camion pseudoprimes

## Statements and evidence

**(a) Camion–Wieferich primes.**

* The *Camion quotient* `q_C(p) = (s(p) − 2^((p−1)/2))/p mod p` is the analogue of the Fermat quotient.
* Its zeros below 3000 are `p = 3` (`D(3) = 0` exactly) and `p = 13` (`13² | D(13)`).
* The random model predicts about `Σ 1/p ≈ log log x` zeros up to `x`. That is infinitely many, but very sparse: the same heuristic as for the Wieferich primes 1093 and 3511.
* Curiosity: `3511 | D(15)`.

**(b) Sporadic classes.** Besides the universal classes `{1, 2, 3}`, a prime `p` has the classes `r ∈ [4, p−1]` with `p | D(r)`.

* The counts over 75 primes fit Poisson with the model mean `(p−4)/p ≈ 0.945` (observed mean 0.947).
* So `D(r) mod p` behaves like a random residue for `r ≥ 4`. The universal classes are exactly the orders at which every strong tournament is a circulant.

**(c) Counting.** The families `2p` and `3p` already give `~5x/(6 log x)`.

* The exotic Camion pseudoprimes are the Carmichael numbers, the Chernick-type `3b(3b−2)`, mixed `{1,2,3}`-residue products, and products using sporadic classes (`1239, 3765, 8401, 12019, 14317` below 20000).
* They are conjectured to be `O(x^θ)` for some `θ < 1`.

**(d) Modulo 9.** The bad odd classes are exactly `m ≡ 1 (mod 6)`, `m ≥ 7`.

* `D` is not periodic modulo `p²` with any period dividing `2p(p−1)`. So a proof needs a mod-`p²` refinement of the group-ring argument, where binomial coefficients are no longer governed by Lucas' theorem alone.
* If (d) holds, the factor `9` never obstructs: every odd `n` divisible by `9` passes modulo 9, because `n ≡ 3 (mod 6)`.

**(e) General principle.**

* Positive instances: tournaments, and semicomplete digraphs with multiplier `(3/p)`.
* Hostile control: connected graphs, whose components are unordered, have no quasi-period.
* Candidate next instances:
  * tournaments with a marked vertex;
  * ordered set partitions weighted by `a^C(block,2)`;
  * `k`-hypertournaments if they have a linear condensation.

## What would refute

* (a): a proof that the Camion quotient vanishes only finitely often.
* (b): a structured excess of sporadic classes, e.g. a residue `r ≥ 4` that is a zero class for a positive proportion of primes.
* (d): an odd `m ≡ 3, 5 (mod 6)` with `D(m) ≢ 0 (mod 9)`.
