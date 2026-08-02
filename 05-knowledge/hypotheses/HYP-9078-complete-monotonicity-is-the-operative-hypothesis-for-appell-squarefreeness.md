---
id: HYP-9078
title: "Complete monotonicity is the operative hypothesis for Appell squarefreeness (and hence for SlotSFC_1(2))"
status: >
  CONJECTURE / EVIDENCE ONLY -- NOT PROVED. THM-3021 showed that SlotSFC_1(2) at
  p = 0 is equivalent to squarefreeness of the Appell sequence
  Phi_n(z) = int (z-u)^n dnu(u), and killed the two obvious hypotheses:
  positivity of the coefficients (H2 counterexample) and positivity of the
  measure nu (H6 counterexample, Phi_4 a perfect square for a 3-atom positive
  measure). This file proposes the hypothesis that survives.
  CONJECTURE (CM). If nu has a COMPLETELY MONOTONE density on (0,inf), then
  every Phi_n is squarefree.
  Motivation: the SFC measure dnu = (1/s) u^{1/s-1} e^{-u^{1/s}} du IS
  completely monotone for every s >= 1 (u^{1/s-1} has exponent <= 0; e^{-u^{1/s}}
  is CM because u^{1/s} is a Bernstein function for 0 < 1/s <= 1; CM is closed
  under products), whereas an atomic measure is not.
  Equivalent form, via Bernstein: CM density <=> dnu = (int e^{-au} dsigma(a))du,
  and then Phi_n(z) = int a^{-n-1} Psi_n(a z) dsigma(a) with
  Psi_n(y) = (-1)^n n! e_n(-y), e_n the truncated exponential -- so the
  conjecture says a POSITIVE MIXTURE OF DILATES of one fixed squarefree
  polynomial is squarefree.
  EVIDENCE: exhaustive numerical search over finite exponential mixtures
  (2 atoms in sigma: exactly determined system, n = 3..7, 6000 restarts each;
  3 atoms: underdetermined, n = 4..6, 6000 restarts each) found ZERO
  positive-weight multiple roots, against hundreds of solutions for atomic nu
  at the same search effort. If CONJECTURE (CM) is true, SlotSFC_1(2) at p = 0
  follows for every s and every window.
source: death-star-2026-07-31-coinC2
depends_on:
  - THM-3021
  - THM-3020
related:
  - THM-3019
  - THM-2836
script: 04-computation/appell_complete_monotonicity_hyp9078.py
output: 05-knowledge/results/appell_complete_monotonicity_hyp9078.out
---

# HYP-9078 -- complete monotonicity is the operative hypothesis

**Notation correction (MISTAKE-350).** `SlotSFC_1(2)` is the two-monomial
restriction inside ambient `SFC(1)`, not the original ambient `SFC(2)`.

## 1. Where this comes from

THM-3021 established

```text
SlotSFC_1(2) at window k, support {0,s}   <=>   Phi_{k+2} squarefree,
Phi_n(z) = int_0^inf (z-u)^n dnu(u),    Phi_n' = n Phi_{n-1},
```

with `dnu(u) = (1/s) u^{1/s-1} e^{-u^{1/s}} du` the pushforward of `e^{-t}dt`
along `u = t^s`, and then **killed the two natural hypotheses**:

* positivity of the coefficients of `A` -- false, `A = (lam+5)(lam+9)`;
* positivity of the measure `nu` -- false, `nu = delta_0 + (81/16)delta_1 +
  (1/16)delta_3` gives `Phi_4(z) = (7z^2-12z+9)^2/8`, a perfect square.

So the operative property is finer than positivity. This file names a
candidate and tests it.

## 2. The conjecture

**CONJECTURE (CM).** If `nu` has a completely monotone density on `(0,inf)`,
then `Phi_n(z) = int (z-u)^n dnu(u)` is squarefree for every `n`.

**The SFC measure qualifies.** `dnu/du = (1/s) u^{1/s-1} e^{-u^{1/s}}` is a
product of two completely monotone functions: `u^{1/s-1}` (a power with
exponent `<= 0`) and `e^{-u^{1/s}}` (completely monotone because `u -> u^{1/s}`
is a Bernstein function for `0 < 1/s <= 1`), and CM is closed under products.
An atomic measure has no density at all, which is exactly why H6 escapes.

**Equivalent mixture form.** By Bernstein's theorem a CM density is a mixture
of exponentials, `dnu(u) = (int_0^inf e^{-au} dsigma(a)) du` with `sigma >= 0`.
Substituting and rescaling `w = au`,

```text
Phi_n(z) = int_0^inf a^{-n-1} Psi_n(a z) dsigma(a),
Psi_n(y) = int_0^inf (y-w)^n e^{-w} dw = (-1)^n n! e_n(-y),
```

with `e_n` the truncated exponential. `Psi_n` is squarefree (that is the
`s = 1` case, THM-3019 S4). **So (CM) says: a positive mixture of positive
dilates of one fixed squarefree polynomial is squarefree.** Dilation by `a > 0`
preserves the arguments of the roots and scales only their moduli, so the
constituent polynomials all have their roots on the same finite union of rays
-- which is presumably the structure a proof should exploit.

## 3. Evidence

Search for a positive-weight multiple root, `Phi_n(z_0) = Phi_{n-1}(z_0) = 0`,
by Newton solve from random restarts:

```text
sigma with 2 atoms   (q = (1,q2) > 0, a = (1,a2) > 0):
    4 real equations in 4 real unknowns -- EXACTLY DETERMINED
    n = 3,4,5,6,7,  6000 restarts each   ->  0 solutions
sigma with 3 atoms   (q = (1,q2,q3) > 0, a = (1,a2,a3) > 0):
    underdetermined, so strictly easier to solve
    n = 4,5,6,      6000 restarts each   ->  0 solutions
```

For comparison, the *same* solver and effort on 3-atom **nu** (H6) returned
33 / 69 / 117 valid positive-weight solutions at `n = 4/5/6`. The contrast is
the evidence: the phenomenon is abundant without CM and absent with it.

## 4. What it would buy, and what is not claimed

If (CM) holds then SlotSFC_1(2) is proved for `p = 0`, every `s >= 1`, every window
`k`, since the SFC measure is CM. That is the whole `p = 0` family at once,
which is strictly more than THM-3020 (K5) achieved (`s = 1` only).

**Not claimed:** (CM) is a conjecture. It is supported only by the search
above, which is finite and numerical and proves nothing. In particular no
claim is made that CM is *necessary*; some weaker property may suffice. The
`p >= 1` families are untouched -- there `I_m = K(pm, m)` and the Appell
identity `Phi_n' = n Phi_{n-1}` is not available in the same form, so the
reformulation of THM-3021 section 5 does not directly apply.

Referee: `appell_complete_monotonicity_hyp9078.py` reproduces both searches
and the sanity check that a single exponential returns `Psi_n` with root
separations bounded away from zero.
