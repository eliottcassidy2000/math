---
id: HYP-9358
title: "Open laws for Camion pseudoprimes (THM-4746). (a) The Camion-Wieferich primes, with s(p) = 2^((p-1)/2) (mod p^2), equivalently (THM-4746 (7)) p^2 a Camion pseudoprime, are infinite in number but of density zero, like Wieferich primes (only 3 and 13 below 10000). (b) The sporadic zero classes {r >= 4 : p | s(r) - c(r)} are Poisson with mean (p-4)/p. (c) The Camion pseudoprimes up to x number 5x/(6 log x)(1 + o(1)), all but o(x/log x) of them of the form 2p or 3p. (d') For every odd prime p and k >= 1, s(m + phi(p^k)) = (2/p) s(m) (mod p^k) for all m >= 2k - 1 (PROVED for every k >= 1 by the independently audited Session25 continuation). (e) Thresholds: for each class of labelled structures obtained from 2^C(m,2) by Hurwitz inversion and shifts, determine the exact index from which the mod-p law holds."
status: >
  PARTIALLY PROVED: (d') is PROVED for all odd primes p and k >= 1 by the Session25 continuation.
  Parts (a), (b), (c), and the exact-threshold classification in (e) remain OPEN; (a) gained the PROVED tower-height
  statement THM-4746 (9) and the open sub-question (a') h(p) <= 2 (mac-mini continuation, 2026-10-21). FINITE-EXACT evidence:
  - (a): p <= 10000 (camion_wieferich_10000.out), extended to p <= 10^5 (THM-4750 (5)): the Camion-Wieferich primes
    are exactly 3, 13, 22409; the equivalence with "p^2 passes" is PROVED by Lemma B of the audit, THM-4746 (7).
  - (b): the 75 primes 7 <= p <= 400. Observed 31, 25, 13, 4, 2, 0 primes with 0..5 sporadic classes, against the
    Poisson prediction 29.1, 27.5, 13.0, 4.1, 1.0, 0.2 (camion_sporadic.out; reproduced exactly by the audit).
    Extended to the 1226 primes 7 <= p <= 10^4 (THM-4750 (4), camion_irregular_10000.out): 423, 473, 250, 67, 11, 2,
    mean 1.002 and variance 0.875 (slightly underdispersed). No odd sporadic class: 0.610 (Siegel-type e^(-1/2) = 0.607);
    no even one: 0.581; neither: 0.345 (e^(-1) = 0.368).
  - (c): up to 20000 there are 2145 squarefree Camion pseudoprimes. Of these 857 are 3p and 1228 are 2p, leaving
    60 exotic ones (37 odd, 23 even).
  - (d'): PROVED for every odd prime p and k >= 1, from m >= 2k - 1, by the Session25 filtered-Hurwitz proof.
    The original k = 3, 4 finite evidence is preserved in audit a17b; k = 1 is THM-4746 (1), and k = 2 is the incoming audit's Lemma B.

    The former (d), the mod-9 law, is now PROVED (THM-4746 (7)).
  - (e): the thresholds known so far are m >= 1 for strong tournaments and strong semicomplete digraphs, and m >= 2
    for connected graphs (mod p). Mod p^2 they are m >= 3 for strong tournaments; for connected graphs the audit's a17
    observed m >= 3 as well.
  AUDITED 2026-10-21 (MISTAKE-651). The pre-audit (e) used connected graphs as a hostile control; they obey the law.
  The pre-audit (d) is now a theorem.
session: mac-mini-2026-10-21-fanopsp
source: 05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md
related:
  - 05-knowledge/results/session25_spectrum_wiring_20261021/session25-congruence-review/ADDITIVE_PRIME_POWER_PROOF.md (all-prime-power additive closure)
  - 05-knowledge/results/session25_spectrum_wiring_20261021/session25-additive-review/REVIEW.md (independent proof and corollary audit)
  - 01-canon/theorems/THM-4746-camion-pseudoprimes-strong-tournament-counts-are-quasi-periodic-mod-p-korselt-criterion-carmichael-numbers-and-the-triangle.md
  - 05-knowledge/results/tournament_pseudoprimes_audit_20261021.md (Lemmas B and C)
  - 01-canon/theorems/THM-4750-the-camion-function-strong-tournament-counts-interpolate-p-adically-at-half-speed-fermat-quotient-curvature-wieferich-collapse-second-order-zeros-in-pairs-and-the-camion-limit.md (continuation: Camion function, tower heights, 22409)
  - 05-knowledge/hypotheses/HYP-9361-camion-function-open-laws-split-inert-equidistribution-tower-heights-uniformity-k-speed.md
scripts:
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_sporadic.py (+ .out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_list.c (+ camion_list_20000.out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_wieferich.c (+ camion_wieferich_10000.out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/camion_verify.c (+ camion_verify.out)
  - 04-computation/experiments/tournament_pseudoprimes_20261021/audit/a17_mod_p2_quasiperiod.py, a17b_mod_pk_quasiperiod.py
---

# HYP-9358 — open laws for Camion pseudoprimes

## Statements and evidence

**(a) Camion–Wieferich primes.**

* The *Camion quotient* `q_C(p) = (s(p) − 2^((p−1)/2))/p mod p` is the analogue of the Fermat quotient.
* Its zeros below 10000 are `p = 3` (`D(3) = 0` exactly) and `p = 13` (`13² ∥ D(13)`).
* `q_C` shows no relation to the Euler quotient of 2 (scratch comparison for `p ≤ 139`). So Camion–Wieferich primes are not Wieferich primes in disguise.
* The random model predicts about `log log x` zeros up to `x`. That is infinitely many, but very sparse: the same heuristic as for the Wieferich primes 1093 and 3511.
* Curiosity: `3511 | D(15)`.
* **Tower height (THM-4746 (9), PROVED).** The exponents `a` with `p^a` a Camion pseudoprime form an initial segment `{1, …, h(p)}`, and `h(p)` is the frozen value of `v_p(D(p^a))`. `h(p) ≥ 2` iff `p` is Camion–Wieferich; `h(3) = h(13) = 2` exactly; `h(p) = 1` for every other odd `p ≤ 10^4`. Open sub-question (a′), the analogue of the absence of Wieferich primes of order 3: is `h(p) ≤ 2` for every `p`?
* **Continuation (THM-4750 (5)).** `h(p) = v_p(Λ_p)` for the Camion limit `Λ_p = lim ε^a D(p^a)`, and `Λ_p mod p^K` is computable from `s(n)`, `n ≤ (2K−1)p`.
  * The search to `10^5` found one new Camion–Wieferich prime, **22409**, with `h(22409) = 2`. So every prime `p ≤ 10^5` has `h(p) ≤ 2`.
  * For `p ≥ 5`, `h(p) ≥ 3` needs the class-1 Camion zero to satisfy `ℓ ≡ p (mod p²)`, a second `1/p` event.
  * Neither Wieferich prime (1093, 3511) is Camion–Wieferich. Wieferich primes instead make the Camion curvature vanish (THM-4750 (3)).
  * The search to `10^5` is the end of the feasible range for `O(p²)` methods. The expected number of further Camion–Wieferich primes in `(10^5, 10^6]` is about `0.2`.

**(b) Sporadic classes.** Besides the universal classes `{1, 2, 3}`, a prime `p` has the classes `r ∈ [4, p−1]` with `p | D(r)`.

* The counts over 75 primes fit Poisson with the model mean `(p−4)/p ≈ 0.945` (observed mean 0.947).
* So `D(r) mod p` behaves like a random residue for `r ≥ 4`. The universal classes are exactly the orders at which every strong tournament is a circulant.

**(c) Counting.** The families `2p` and `3p` already give `~5x/(6 log x)`.

* The exotic Camion pseudoprimes are:
  * the Carmichael numbers;
  * the Chernick-type `3b(3b−2)`;
  * mixed `{1, 2, 3}`-residue products;
  * products using sporadic classes: odd `1239, 3765, 8401, 12019, 14317` and even `186, 530, 1898, 3290, 7202, 8282, 17066` below 20000.
* They are conjectured to be `O(x^θ)` for some `θ < 1`.

**(d′) Higher prime powers — PROVED.** The [Session25 filtered-Hurwitz/Newton proof](../results/session25_spectrum_wiring_20261021/session25-congruence-review/ADDITIVE_PRIME_POWER_PROOF.md), with [independent audit](../results/session25_spectrum_wiring_20261021/session25-additive-review/REVIEW.md), proves `s(m + phi(p^k)) ≡ (2/p)s(m) (mod p^k)` and the same law for `D`, for every odd prime `p`, `k ≥ 1`, and `m ≥ 2k − 1`. This extends the incoming audit's square-case Lemma B; it does not claim that the sufficient threshold is always minimal. The full prime-power CRT criterion in the same proof removes the squarefree/cube-free restriction for the aggregate discrepancy.

**(e) Exact thresholds.** The proved sufficient thresholds are `m ≥ 2k − 1` for the strong count (Session25; THM-4746 §7) and `n ≥ 2k` for connected graphs (one shift). Let `τ(p,k)` be the least `m_0` such that `s(m + φ(p^k)) ≡ (2/p)s(m) (mod p^k)` for every `m ≥ m_0`. Since the law holds from `2k − 1`, `τ(p,k) − 1` is the largest failing index below it.

* **Data (FINITE-EXACT, `modpk/pk_period_*.out`).**

| `p` | `k = 2` | `k = 3` | `k = 4` | `k = 5` | `k = 6` | `k = 7..10` |
|---|---|---|---|---|---|---|
| 3 | 3 | 3 | 5 | 7 | 9 | 9, 9, 11, 13 |
| 5 | 3 | 4 | 5 | 7 | 8 | — |
| 7 | 2 | 5 | 6 | 7 | — | — |
| 11 | 3 | 5 | 7 | — | — | — |
| 13 | 2 | 5 | 7 | — | — | — |
| 17, 19, 23 | 3 | 5 | — | — | — | — |
| 29 … 97 | 3 | — | — | — | — | — |

  The sufficient value is `2k − 1 = 3, 5, 7, 9, 11, 13–19`.

* **Conjecture (e1), generic sharpness.** For each fixed `k`, `τ(p,k) = 2k − 1` for all primes `p` outside a set of density zero.
  * The exceptions at `k = 2` (`p = 7, 13`, where `τ = 2`) are accidental zeros of `s(2 + φ(p²)) mod p²`.
  * Heuristically, an exception at level `k` has probability about `1/p`.
* **Conjecture (e2), small primes.** For `p < 2k − 1`, `τ(p,k) < 2k − 1`. Every tested case fits: `p = 3` with `k ≥ 3`, `p = 5` with `k ≥ 4`, `p = 7` with `k = 5`. At `p = 2k − 1` both behaviours occur: `τ(3,2) = 3`, `τ(5,3) = 4`, `τ(7,4) = 6`. The reason is presumably extra `p`-divisibility of the delta coefficients `β_(0,m)`. Mod `p^k`, `s(m + φ(p^k)) − (2/p)s(m) ≡ (2/p)·β_(0,m)`, where `β_(0,m)` is the `δ_m` coefficient of `e − S` in the filtered algebra. A plausible source of the extra divisibility is the odd double factorials `(2j−1)!!` in the Newton coefficients of `C(C(m,2), j)`. This mechanism is HEURISTIC; no formula for `τ(3,k)` is proposed. THM-4750 §2 confirms the double-factorial mechanism for the Mahler speed: the extra divisibility at `p = 5, 7` sits exactly at `i ∈ {2k−1, 2k}` with `p | (2k−1)!!`. The threshold statement (e2) itself remains open.
* **Connected graphs.** The failures in 15 tested pairs stop at `n ≤ 2k − 2`, one below the proved `2k` (`modpk/cg_pk.out`).
* **Hostile control.** Rooted counts such as `m·s(m)` are genuinely not `(p−1)`-quasi-periodic. Their quasi-period is `p(p−1)` (`camion_controls.out`).

## What would refute

* (a): a proof that the Camion quotient vanishes only finitely often.
* (b): a structured excess of sporadic classes, e.g. a residue `r ≥ 4` that is a zero class for a positive proportion of primes.
* The former (d′) falsification target is closed by the proof above; the remaining open parts retain their stated scope.
* (a′): a prime with `p³ | D(p²)`, i.e. `h(p) ≥ 3`.
* (e1): a positive proportion of primes with `τ(p,k) < 2k − 1` at some fixed `k`.
