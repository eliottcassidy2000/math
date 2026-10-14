---
id: THM-4633
title: "Collatz cycles as rational approximations of theta = log2 3 with an admissible-minimum budget: a Syracuse cycle with k odd elements, total exponent K and minimum m satisfies 2^K/3^k = prod (1 + 1/(3 n_i)), hence 0 < K/k - theta < 1/(3 m ln 2); every reduced fraction p/q > theta has the budget m_max(p/q) = q/(3(p ln 2 - q ln 3)) (about a_(j+1) q_j^2/(3 ln 2) at a convergent), and a cycle with K/k = p/q has m < m_max; if k^2 <= (3/2) m ln 2 then K/k is a convergent (Legendre); with the verified bound m >= 2^71 the convergents p_1/q_1 = 2/1 (budget 1.16, realised by the trivial cycle at m = 1), 8/5, 65/41, 485/306, 24727/15601, 125743/79335, 301994/190537, 17087915/10781274, 272500658/171928773, 630138897/397573379 and 10439860591/6586818670 are excluded and the first survivor is 217976794617/137528045312, so every nontrivial cycle has more than 4.95 10^10 odd elements (Legendre route) and at least 137,528,045,312 (Hercher's Corollary 29, the same convergent); an effective exponent |theta - p/q| > c q^(-mu_0) gives k > (3 c m ln 2)^(1/mu_0), so the owner's mu(theta) = 2 would force k of order sqrt(m), which the computed continued fraction (partial quotients at most 55 through index 40) already delivers in the verified window"
status: >
  PROVED (statements 1-3 and 6: elementary; statement 5 first clause via Legendre's theorem, the table
  and Barina's verification). FINITE-EXACT: the convergent table (120-digit arithmetic, exact integers
  p_j, q_j, sign of p_j ln 2 - q_j ln 3, budgets) in log23_convergents.py. CITED: the product identity
  and the approximation inequality are classical (Steiner 1977; Eliahou 1993, whose method this is;
  Lagarias's survey); the verification to 2^71 (Barina, J. Supercomput. 81 (2025) art. 810); the
  odd-element bound 137,528,045,312 (Hercher, J. Integer Seq. 26 (2023), Corollary 29, with corrigendum;
  his finer intermediate-fraction analysis replaces the Legendre threshold), which equals q_23 of this
  table; the effective-measure literature as marked. NEW here: the budget function m_max as the
  dictionary 'convergent -> admissible minimum', the explicit exclusion ledger at 2^71, the Legendre
  route as a self-contained proof of a 4.95 10^10 bound, and the irrationality-measure parametrisation
  (statement 6) with its comparison against the actual continued fraction. Not independently audited yet.
session: mac-mini-2026-10-14-irrmeasure
source: 05-knowledge/results/collatz_irrationality_measures_20261014.md
scripts:
  - 04-computation/experiments/irrmeasure_20261014/log23_convergents.py (+ .out)
related:
  - THM-4634 (parity-number dictionary), THM-4591 (integer T-cycles to period 301,993; the convergent 301994/190537 is its period bound), THM-4484 (free and sporadic cycles), THM-4516 (perfect-power clocks)
  - HYP-9260 / HYP-9261 (convergents and semiconvergents of theta as resonance frequencies; MISTAKE-602 on semiconvergents)
  - HYP-9263 (mu = 2 for logarithms)
  - S. Eliahou, The 3x+1 problem: new lower bounds on nontrivial cycle lengths, Discrete Math. 118 (1993) 45-56
  - C. Hercher, There are no Collatz m-cycles with m <= 91, J. Integer Seq. 26 (2023), Article 23.3.5 (and corrigendum); arXiv:2201.00406
  - D. Barina, Improved verification limit for the convergence of the Collatz conjecture, J. Supercomput. 81 (2025) art. 810 (2^71)
  - R. P. Steiner, A theorem on the Syracuse problem, Proc. 7th Manitoba Conf. (1977); J. Simons, B. de Weger, Theoretical and computational bounds for m-cycles of the 3n+1 problem, Acta Arith. 117 (2005)
  - G. Rhin, Approximants de Pade et mesures effectives d'irrationalite, Progr. Math. 71 (1987) (mu(gamma) < 8.616 for nonzero gamma in Q log 2 + Q log 3)
---

# THM-4633 — cycles as rational approximations of log₂3

## Setting

`U(n) = (3n + 1)/2^(v_2(3n+1))` is the Syracuse map on odd positive integers. A cycle `n_0 → n_1 → ... →
n_(k-1) → n_0` has `k` odd elements, exponents `e_i = v_2(3n_i + 1)`, total exponent `K = Σ e_i` (the period of
the corresponding `T`-cycle, THM-4634), and minimum `m = min n_i`. `θ = log_2 3 = 1.5849625...`, with
continued fraction `[1; 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 15, 1, 9, 2, 5, 7, 1, 1, 4, 8, 1, 11, 1,
20, 2, 1, 10, 1, 4, 1, ...]` and convergents `p_j/q_j` (`p_0/q_0 = 1/1`, `p_1/q_1 = 2/1`, `3/2`, `8/5`, `19/12`,
`65/41`, `84/53`, `485/306`, `1054/665`, `24727/15601`, ...; the odd-indexed ones lie above `θ`).

## Statement

1. **Product identity and approximation (KNOWN).** `2^(e_i) n_(i+1) = 3 n_i + 1` multiplied around the cycle
   gives `2^K/3^k = Π_i (1 + 1/(3 n_i))`, hence

       0 < K ln 2 - k ln 3 < k/(3m)      and      0 < K/k - θ < 1/(3 m ln 2),

   with equality `2^2/3 = 1 + 1/3` exactly for the trivial cycle (`k = 1`, `K = 2`, `m = 1`).
2. **Budget.** For a reduced fraction `p/q > θ` put `m_max(p/q) := q/(3 (p ln 2 - q ln 3))`. Every cycle whose
   ratio `K/k` reduces to `p/q` (any common multiple `K = dp`, `k = dq`) has `m < m_max(p/q)`; at a
   convergent, `m_max(p_j/q_j) = a_(j+1) q_j^2/(3 ln 2) · (1 + O(1/a_(j+1) q_j))`.
3. **Legendre.** If `k^2 ≤ (3/2) m ln 2`, then `|θ - K/k| < 1/(2k^2)` and `K/k` in lowest terms is a convergent of
   `θ`, necessarily one with `p_j/q_j > θ`.
4. **Exclusion ledger at `2^71` (FINITE-EXACT).** For the convergents above `θ`:

   | j | p_j (= K) | q_j (= k) | budget m_max | status vs 2^71 ≈ 2.36·10^21 |
   |---|---|---|---|---|
   | 1 | 2 | 1 | 1.159 | realised by the trivial cycle (m = 1 ≤ 1.159) |
   | 3 | 8 | 5 | 31.98 | excluded |
   | 5 | 65 | 41 | 1192 | excluded |
   | 7 | 485 | 306 | 99,781 | excluded |
   | 9 | 24727 | 15601 | 2.86·10^8 | excluded |
   | 11 | 125743 | 79335 | 7.22·10^9 | excluded |
   | 13 | 301994 | 190537 | 9.85·10^11 | excluded |
   | 15 | 17087915 | 10781274 | 2.94·10^14 | excluded |
   | 17 | 272500658 | 171928773 | 3.20·10^16 | excluded |
   | 19 | 630138897 | 397573379 | 1.25·10^18 | excluded |
   | 21 | 10439860591 | 6586818670 | 2.17·10^20 | excluded |
   | 23 | 217976794617 | 137528045312 | 5.10·10^22 | alive |
   | 25 | 8573543875303 | 5409303924479 | 2.73·10^25 | alive |

   (the even-indexed convergents lie below `θ` and cannot be cycle ratios, since `2^K > 3^k`).
5. **Cycle length.** Let `m ≥ 2^71` for every nontrivial cycle (Barina 2025: every `n < 2^71` reaches 1).
   Then a nontrivial cycle cannot have `k ≤ sqrt((3/2)·2^71 ln 2) = 4.95·10^10` odd elements: by 3 its
   ratio would be a convergent, by 4 a surviving one, so `k ≥ q_23 = 1.375·10^11`, a contradiction. Hence
   **every nontrivial Syracuse cycle has more than `4.95·10^10` odd elements** (PROVED here), and
   **at least `137,528,045,312`** (Hercher, Corollary 29; CITED), which is `q_23`, the first surviving
   convergent denominator — the exclusion ledger is the mechanism behind Hercher's number.
6. **Irrationality-measure parametrisation.** If `|θ - p/q| > c q^(-μ_0)` for all `p/q` (an effective
   irrationality exponent `μ_0` with constant `c`), then every cycle with minimum `≥ M_0` has

       k > (3 c M_0 ln 2)^(1/μ_0).

   With `M_0 = 2^71` and `c = 1`: `μ_0 = 2` gives `k > 7.0·10^10`; `μ_0 = 3`: `1.7·10^7`; `μ_0 = 5.117` (the
   exponent that a three-term linear-form bound `|p + q_1 ln 2 + q_2 ln 3| > H^(-4.117)` would give for `θ`,
   reported in the literature search and NOT verified here): `1.7·10^4`; `μ_0 = 8.616` (Rhin's measure for
   the numbers `a ln 2 + b ln 3`, a different statement, quoted for scale): `329`. The exact continued
   fraction does better than any of these: `m_max(p_j/q_j) ≍ a_(j+1) q_j^2`, and the partial quotients
   through index 40 are at most 55, so in the verified window `θ` behaves exactly like a number of
   exponent 2 — the owner's hypothesis `μ(θ) = 2` (HYP-9263) is what the convergent ledger realises
   numerically, and it would add nothing effective beyond the ledger.

## Proofs

**(1).** `2^(e_i) n_(i+1) = 3n_i + 1 = 3 n_i (1 + 1/(3n_i))`; the product over the cycle telescopes the
`n_i`. Taking logarithms, `K ln 2 - k ln 3 = Σ ln(1 + 1/(3n_i)) ∈ (0, k ln(1 + 1/(3m))] ⊂ (0, k/(3m))` for
`m ≥ 1`, with equality in the first inclusion iff all `n_i = m`, i.e. the trivial cycle. Dividing by `k ln 2`
gives the second form.

**(2).** If `K/k = p/q` with `K = dp`, `k = dq`, then `K ln 2 - k ln 3 = d(p ln 2 - q ln 3) < dq/(3m)`, so
`m < q/(3(p ln 2 - q ln 3))`. At a convergent, `|q_j θ - p_j| = 1/(q_(j+1) + q_j θ_(j+1)')` lies between
`1/(2 q_(j+1))` and `1/q_(j+1)` with `q_(j+1) = a_(j+1) q_j + q_(j-1)`, giving the stated size.

**(3).** Legendre: a reduced `p/q` with `|θ - p/q| < 1/(2q^2)` is a convergent. The hypothesis gives
`1/(3 m ln 2) ≤ 1/(2k^2)`, and `K/k > θ` selects the convergents above `θ`.

**(4).** `log23_convergents.py` computes `θ` to 120 digits, the partial quotients (the integers listed; the
first 40 are reliable far beyond `q_37 ≈ 4.6·10^18`), the exact integer convergents, the sign of
`p_j ln 2 - q_j ln 3` and the budgets. (For reference, `p_13/q_13 = 301994/190537` is the convergent behind
THM-4591's period bound `301,993`, and `p_15/q_15 = 17087915/10781274` is Eliahou's 1993 bound from the
verification limit of his day.)

**(5).** If `k ≤ 4.95·10^10 = sqrt((3/2)·2^71 ln 2)` then `k^2 ≤ (3/2) m ln 2` for `m ≥ 2^71`, so by (3) `K/k`
reduces to a convergent `p_j/q_j > θ` with `q_j ≤ k ≤ 4.95·10^10`; the ledger (4) shows every such
convergent has `m_max < 2^71 ≤ m`, contradicting (2). Hercher's Corollary 29 reaches `q_23` itself by
controlling the intermediate fractions between consecutive convergents (his `m`-cycle machinery with the
`2^71` input); the agreement `q_23 = 137,528,045,312` is the consistency check.

**(6).** Combine `c k^(-μ_0) < |θ - K/k| < 1/(3 m ln 2)`.  ∎

## Reading

* **The budget is the dictionary.** Each convergent above `θ` is a "constructed approximation"
  `p_j/q_j` with a price: a cycle realising it must have all its elements below `m_max ≍ a_(j+1) q_j^2`.
  The trivial cycle is the convergent `2/1` realised at the only integer `m = 1 ≤ 1.159`; every later
  convergent up to `j = 21` is paid for only by elements below `2^71`, which verification forbids. In the
  owner's words, `1 → 2 → 1` is the only "abnormality" because `2/1` is the only convergent whose budget
  admits an integer below the verified limit.
* **What irrationality measures can and cannot do.** An effective exponent `μ_0` turns the verified limit
  `M_0` into the odd-length bound `(3 c M_0 ln 2)^(1/μ_0)`; the best conceivable, `μ_0 = 2`, gives
  `k ≳ sqrt(M_0)`, and the actual continued fraction already achieves that scale (the ledger). No
  Diophantine statement about `θ` alone excludes cycles with `k ≫ sqrt(M_0)`, because for such `k` the
  interval `(θ, θ + 1/(3 M_0 ln 2))` contains rationals with denominator `k` for trivial reasons; this is
  why the `m`-cycle theorems (Simons–de Weger, Hercher) need the cycle's internal structure and not only `θ`.
* **Divergence is untouched.** The dictionary of THM-4634 places divergent trajectories at irrational
  parity numbers; `θ`'s Diophantine properties say nothing about them.
