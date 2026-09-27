---
id: HYP-9164
title: "Bounded shadow error: for a divergent Syracuse orbit m_l with real shadow xi 3^l/2^(d_l) = m_l + eta_l, the error eta_l = sum_(k>=0) 2^(d_(l+k)-d_l)/3^(k+1) (the real Bernstein value of the tail word, obeying 2^(v_(l+1)) eta_(l+1) = 3 eta_l - 1) is unbounded: for every divergent orbit on either sheet sup_l eta_l = infinity, i.e. every divergent orbit has unbounded valuations or unboundedly deep dips relative to the current point (the positive cycles of the 3x-1 sheet have bounded error, so the statement concerns divergent orbits only). Trivial bound: sup eta_l >= 5/3, sharp over words."
status: >
  OPEN intermediate statement, implied by the no-divergence half of the Collatz
  conjecture on the respective sheet (on the plus sheet a convergent orbit has eta = infinity; on the minus sheet
  the positive cycles 1, 5-7, 17-...-91 and the orbits entering them have
  bounded error equal to the cycle point, so the statement concerns divergent
  orbits only; a first draft said 'no positive integer orbit', corrected by the
  audit). PROVED context: the recursion, eta >= 1,
  the budget 2^v <= 3 eta - 1 and the dip bound m_(l+k) >= m_l/(3 eta_l); for
  sup eta < B the valuations are <= log_2(3B - 1) and for B = 2 the word is the
  itinerary of f(eta) = (3 eta - 1)/2 on [1, 5/3), (3 eta - 1)/4 on [5/3, 2)
  (eight periodic itineraries with p <= 10, none with a positive plus-sheet
  integer point; the itinerary (1) is the minus-sheet fixed point 1). The sub-case 'no divergent orbit has bounded valuations' is the
  Mahler-shaped (bounded-cut) case. Not implied by any density statement.
source: opus-2026-09-27 session collatz-shadow-flp-20260927
related: [THM-4476 (Bernstein series, two places), HYP-9161 (shell revisits), S13 note collatz_directions_20260926.md (E_inf), the Mahler frontier]
verification: 04-computation/experiments/collatz_shadow_error_20260927.py -> .out; collatz_shadow_error_20260927_itineraries.py -> .out (490 sampled itineraries, none with a positive plus-sheet integer point below 2^41; at least three v = 1 steps after every v = 2); audit collatz_shadow_error_20260927_audit.py -> .out
---

# HYP-9164 -- the shadow error of a divergent orbit is unbounded

See `05-knowledge/results/collatz_shadow_error_flp_20260927.md`, sections 1-3.

**Statement.** Let `n` be a positive odd integer whose Syracuse orbit is
divergent, so that `eta_l = sum_(k>=0) 2^(d_(l+k) - d_l)/3^(k+1) < infinity` for
every `l` (THM-4476). Then `sup_l eta_l = infinity`. Same for the `3x-1` sheet.

**Why it is the right intermediate.** `eta_l >= 1` always; `2^(v_(l+1)) <= 3 eta_l
- 1`; `m_(l+k) >= m_l/(3 eta_l)`. So a bounded error means bounded valuations
and no deep dips: the orbit would climb at a controlled pace using only
halvings of bounded depth, and its word would be the itinerary of an interval
map with finitely many rational-slope branches. The statement is pointwise,
weaker than the conjecture, and is what a Flatto-Lagarias-Pollington-type
argument (the fractional parts of a real shadow cannot be confined) would
have to establish first, in the form "the error cannot be confined to
`[1, B)`". The trivial case `B = 5/3` holds (all errors below `5/3` force the
word `1^infinity`).

**Evidence and non-evidence.** The eight periodic `f`-itineraries with `p <= 10`
have no positive plus-sheet integer point (`(1)` is the minus-sheet fixed point
`1`, an integer itinerary with `eta = 1`); along orbits reaching `1` the
truncated errors are unbounded in the sense that they grow with the depth of
the coming fall (`332` at `m = 3077` on the orbit of `27`); along the `-5`
shadow families the error is `5(1 - (8/9)^k)`. None of this decides the
hypothesis.

**What would refute it.** A divergent orbit whose valuations are bounded and
whose value never falls below a fixed fraction of any earlier value.
