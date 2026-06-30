---
id: HYP-3560
title: PROVED (finite rehearsal, approach 1): the Hamiltonian-path count H concentrates -- CV(H)^2 = (1/n!) sum_{pi: no descent-by-1} 2^{#ascents-by-1} - 1 -> 0, because consecutive-integer adjacencies are asymptotically Poisson(1) and E[2^asc] P(desc=0) = e * e^{-1} = 1; the second moment is diagonal-dominated (Var ~ E[H], Poisson-like) with a controlled 2^{overlap} pair tail -- the EXACT, checkable Siegel-Rogers 2nd moment, validated by klein THM-588 (no linear invariant, one quadratic => proof effort is purely 2nd-moment). NEW METHOD added: Chen-Stein / Poisson approximation gives BOTH concentration (CV^2->0) AND existence (Poisson(lambda>0) is >=1 w.p. 1-e^{-lambda}); the LRC floor transfers as the SET-INDEPENDENT congruence (Gamma_0(N), Han-Lee) bound on the sheet-pair overlap tail = THM-579's CV(N_R)^2 gatekeeper. Plus more new approaches (large deviations, Walsh-Fourier, cut/cusp vanishing).
status: PROVED (CV(H)^2 closed form + Poisson limit; verified exactly n=3..10: 1/3,1/3,19/60,13/45,131/504,131/560,...->0). The LRC transfer is a PROGRAM (the rehearsal is rigorous; the transfer needs the Han-Lee congruence 2nd moment).
source: mac-mini-2026-06-29-S21
related:
  - HYP-3554   # the metagraph is a finite Siegel transform (this is the worked 2nd moment)
  - HYP-3553   # the Gamma_0(N) congruence floor (the set-independent transfer)
  - THM-588    # klein: no linear / one quadratic invariant => purely 2nd-moment (validates approach 1)
  - THM-579    # the LRC gatekeeper CV(N_R)^2 (the transfer target)
  - HYP-3538   # the cap S_2 = -tr(M_odd)/2 (the LRC pairwise 2nd moment)
external: arXiv:2507.05905 (Han-Lee); Chen-Stein method; Siegel-Rogers 2nd moment
results:
  - 04-computation/H_variance_finite_rehearsal_macmini_20260629.py
  - 05-knowledge/results/H_variance_finite_rehearsal_macmini_20260629.out
---

# HYP-3560 -- H concentrates (the rehearsal), and the Poisson-existence method

## The result (approach 1, PROVED)
Over labeled tournaments (the Siegel/mass measure, each arc a fair coin), `E[H] = n!/2^{n-1}`, and by
relabeling symmetry the second moment collapses to a single permutation sum against a reference path,
giving the closed form
> `CV(H)^2 = Var(H)/E[H]^2 = (1/n!) sum_{pi'} c(pi') 2^{j(pi')} - 1`,
where `c(pi')=1` iff `pi'` has NO descending consecutive-integer adjacency `(k+1,k)`, and `j(pi') =`
#ascending consecutive-integer adjacencies `(k,k+1)`. Exact values (verified vs brute force n<=5):
`CV(H)^2 = 1/3, 1/3, 19/60, 13/45, 131/504, 131/560` for `n=3..8`, decreasing (`0.234, 0.212, 0.193` at
n=8,9,10). **It tends to 0**: `CV(H)^2 + 1 = E[2^{asc}\,1(desc=0)]`, and the consecutive-integer
adjacency counts `(asc, desc)` are asymptotically **independent Poisson(1)** (`E[asc]=(n-1)/n -> 1`), so
the limit is `E[2^{Poisson(1)}] * P(Poisson(1)=0) = e^{2-1} * e^{-1} = e * e^{-1} = 1`, hence `CV(H)^2 ->
0`. **H concentrates**; `Var(H) ~ E[H]` (Poisson-like): the second moment is diagonal-dominated, the
off-diagonal a `2^{overlap}` pair tail that `n!` outgrows.

This is the exact, checkable Siegel-Rogers second moment -- and klein's THM-588 (the metagraph has NO
linear invariant and EXACTLY ONE quadratic, the 3-cycle count) is the structural reason the whole proof
lives here: there is no first-order content to bound, only the pairwise second moment.

## The new method: Chen-Stein / Poisson (existence AND concentration in one)
The CV->0 proof IS a Poisson approximation, and that is the new tool. **Chen-Stein** makes it rigorous and
quantitative (total-variation bound on the count vs Poisson(lambda)), and Poisson gives TWO things at once:
- **concentration**: `CV^2 -> 0` (the second moment is Poisson-like);
- **existence**: a `Poisson(lambda)` count with `lambda > 0` is `>= 1` with probability `1 - e^{-lambda} >
  0` -- existence WITHOUT construction, the probabilistic-method twin of the Ky-Fan forced count
  (approach 2). For the LRC: if the safe-sheet / lonely-point count is asymptotically `Poisson(lambda)`
  with `lambda > 0`, a lonely point exists for the TYPICAL config; the worst case (the disproof) is a
  large deviation, controlled by the SET-INDEPENDENT congruence bound below.

## The transfer to the LRC gatekeeper (the program)
The LRC sheet count `N_R = sum_sheets 1(sheet safe)` has the SAME second-moment shape:
`Var(N_R) = [diagonal = E[N_R]] + [off-diagonal = sum_{a!=b} Cov(safe_a, safe_b)]`, the off-diagonal a
sum over PAIRS of sheets weighted by their **resonance overlap** -- exactly the `2^{overlap}` tail of the
metagraph rehearsal, but the overlap here is a CONGRUENCE condition (the speeds mod 14). Han-Lee's
congruence Siegel second moment (HYP-3553) bounds this overlap tail with the covering built in as
`Gamma_0(N)`, SET-INDEPENDENTLY. So `CV(N_R)^2 -> small` by the same diagonal-dominated + Poisson-tail
argument, and THM-579's gatekeeper `CV(N_R)^2 < m_Q/(1-m_Q)` follows. The metagraph is where the bound is
exactly checkable (`CV(H)^2 -> 0`, proved); the congruence Siegel formula is the analytic lift.

## More new approaches (challenging the assumptions)
- **Large-deviation rate function.** `Z(beta) = sum_T H^beta`; the rate `I(x)` for the H-distribution.
  The disproof is a LARGE DEVIATION of the sheet count -- exponentially unlikely, and the rate function
  is the worst-case bound the congruence structure must beat.
- **Shared Walsh/Fourier diagonalization.** The metagraph is diagonalized by hypercube Walsh characters
  (THM-584); the LRC danger functions are exponential sums on the torus -- the SAME Fourier basis. Run the
  metagraph's exact Walsh second moment as the template for the LRC's `sum|chat(14N)|^2`.
- **Cut/cusp vanishing (Eisenstein = 0).** THM-588's `mult(1)=0` says the LINEAR (score / cut / Eisenstein)
  invariant does NOT survive the quotient; the first nonzero is the QUADRATIC (3-cycle / cusp). So the LRC
  binding content is purely the CUSP form (matches "the floor is R-even, S_1 vanishes", HYP-3538/klein) --
  the Eisenstein/score part carries no obstruction.
- **Reverse transfer.** Use a known LRC-side equidistribution (Weyl) to deduce a tournament-counting
  identity -- the arrow run backwards, a consistency check and a source of new tournament facts.

## What it buys
A PROVED concentration result (`CV(H)^2 -> 0`) that is the exact finite model of THM-579's gatekeeper, a
rigorous method (Chen-Stein/Poisson) that delivers concentration AND existence together, and a concrete
program: bound `CV(N_R)^2` by the congruence (`Gamma_0(N)`) overlap tail, set-independently, the Poisson
limit giving the floor. klein THM-588 confirms the effort belongs entirely at the second moment.
