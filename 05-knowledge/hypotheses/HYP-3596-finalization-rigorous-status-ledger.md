---
id: HYP-3596
title: FINALIZATION / honest status ledger of the klein metagraph->cusp arc -- what is PROVED (THM-584/587/588/589/590), what is exactly VERIFIED (the counts, the modular data, the binding floor constant), what is CONJECTURAL REDUCTION (rho_j = the Z_7-core gap; the Gamma_0(14) Eisenstein (+) f_14 decomposition; the uniform floor R'>=3/pi^2), and the explicit statement that LRC(14) itself is OPEN and NOT proved by this work
status: FINALIZATION (a status accounting, per CLAUDE.md 'mark VERIFIED or CONJECTURE; never claim proved without the proof'). The PROVED/VERIFIED items are rigorous; the reduction and LRC(14) are explicitly OPEN.
source: klein-2026-06-29-S15
depends_on:
  - THM-590    # the rigorous apex cyclotomic gap (the self-contained core)
related:
  - HYP-3593   # the truth across frames (the one cyclotomic atom)
  - HYP-3587   # genus = local-global gap; the f_14 reduction (conjectural)
  - HYP-3575   # mac-mini: rho_j = the Z_7 Gram gap (the reduction)
  - HYP-3415   # the LRC14 critical path (kps)
---

# HYP-3596 — finalization: the honest rigorous status

This finalizes the klein arc (metagraph -> signed cycle index -> floor -> finite cyclotomic min -> cusps ->
genus -> the one obstruction atom) by separating, precisely, the rigorous from the conjectural. **LRC(14) is
NOT proved by this work**; what is established is a rigorous structural framework + constants + a (conditional)
reduction to one cusp-form bound.

## A. PROVED (canon theorems, with proofs)

- **THM-584** complement = the antipodal map of the arc-hypercube `Q_{C(n,2)}`; the iso-class metagraph is its
  `S_n`-quotient; eigenvalues `d-2k`; `R`-even/`R`-odd = even/odd hypercube level; mod-4 eigenvalue law.
- **THM-587** the per-level signed cycle index `P_n(x) = (1/n!) sum_sigma prod_cyc (1 + s_c x^{ell_c})`;
  `P_n(1) = A000568`, `P_n(-1) = SC(n)` (Lefschetz/antipodal Euler number).
- **THM-588** the arc-flip metagraph Laplacian spectrum `{2k}`; algebraic connectivity `= 4` for all `n>=3`
  (`mult(1)=0` proved; cyclicity is the level-2 invariant); walk gap `= 4/C(n,2)`; Fiedler = 3-cycle count.
- **THM-589** `CV(H)^2` even-run closed form (odd overlap cancels by orientation parity); `CV(H)^2 ~ 2/n`.
- **THM-590** the apex-7 cyclotomic gap: `g(O)=0` iff `O in {empty,Z_7}`; else `g(O) >= 4cos^2(3pi/7)`,
  equality iff `|O| in {2,5}` (the doublet). The minimal nonzero apex obstruction is `4cos^2(3pi/7)`.

## B. EXACTLY VERIFIED (finite computations; rigorous as stated, scope noted)

- The three "even-graph" counts (HYP-3591/3592): A000088 (all), A000568 (tournaments = Royle-even = odd-order
  Burnside = `P_n(1)`), A002854 (Eulerian/degree-even = cycle space); the sandwich `Euler <= tourn <= all`;
  `BLUE(n)=2^{e}`, `e = k^2`(odd)/`k(k-1)`(even); **blue-classes = SC exactly** (n<=6); pure-black = NS.
- The modular data (HYP-3586): `genus(X_0(2p)) = 0,0,1,2,2` (`p=3,5,7,11,13`), jumping `0->1` at `N=14`;
  the 4 cusps of `X_0(14)` form the Klein group under Atkin-Lehner `W(14)=(Z/2)^2`; `nu_2 = 0` iff apex
  `= 3 mod 4` (Paley). (Standard modular formulas.)
- `CV(N_R)^2` is set-dependent and unbounded as `m_R -> 0` (HYP-3554) -- the variance route cannot be uniform.
- The exact binding floor constant (HYP-3593): `inf_{scan} R' = 114382/332563 = 0.343941` at `R={1..13}\{7}`,
  `Q={1,2}`, which **clears** the clean bound `1/(2 zeta(2)) = 3/pi^2 = 0.303964` by `+0.040`.
  SCOPE: this is the minimum over a 942-covering scan, not a proof of the infimum over the full family.
  CORRECTION (HYP-3597): `inf R'` over the FULL (infinite) covering family is actually `0` -- the lonely
  MEASURE `m_S` vanishes as `S` grows (verified `0.344->0.253->0.173->0.107->0`), but this does NOT break
  LRC (LRC is EXISTENCE, not positive measure; the over-large `S` trivially cover and are outside the
  relevant family). So `0.344` is a scan-slice value, NOT a family infimum. The PROVABLE positive floor is
  the FINITE-family minimum `4cos^2(3pi/7)` (THM-590), reached by the descent's finitization -- the measure
  infimum over the infinite family is the wrong object.

## C. CONJECTURAL REDUCTION (NOT proved; the program)

- **The reduction `rho_j = g(O_j)`**: that the LRC per-level decorrelation equals the apex gap of the
  2-adic-descended `Z_7`-core (mac-mini HYP-3575/3576; THM-580 descent). Conditional.
- **The automorphic decomposition**: that the `Gamma_0(14)` second moment splits Eisenstein (bulk) `(+)`
  cusp form `f_14`, with the obstruction = the genus-1 cusp form's apex-cusp coefficient (HYP-3587, mac-mini
  HYP-3553). Conjectural.
- **The uniform floor `R' >= 3/pi^2`** over ALL coverings (HYP-3550/3553; verified on scans only).
- **"hardness = genus"** (HYP-3586): a structural reframe on 2 relevant data points (`N=6,14`).

## D. OPEN

- **LRC(14) itself** -- a case of the lonely runner conjecture -- is OPEN and is NOT proved here. The
  critical path (kps HYP-3415) = [q-witness + LRC<=13 induction] + [the covering floor `R'>0` uniformly]; the
  floor is the open piece, and within it the open content is exactly **C** above.

## The honest final statement

> This arc rigorously establishes the metagraph/antipodal/signed-cycle-index framework (THM-584/587/588/589),
> the apex-7 cyclotomic gap `4cos^2(3pi/7)` (THM-590), the modular geometry of `X_0(14)` (genus 1, cusps =
> Klein), and the exact binding floor constant. It REDUCES the LRC(14) covering floor -- CONDITIONALLY on the
> descent-to-`Z_7`-core reduction -- to (i) the positivity of a single cyclotomic atom, which is now a
> theorem (`4cos^2(3pi/7) > 0`, THM-590), and (ii) bounding the genus-1 cusp-form component against the
> `3/pi^2` Eisenstein bulk. It does NOT prove LRC(14). The remaining mathematics is the conjectural
> reduction (C) -- one structural step (the reduction) and one analytic step (the cusp-form bound).

So: the *self-contained truth* of the arc -- "`4cos^2(3pi/7) > 0` is the minimal apex obstruction, and it is
positive" -- is now PROVED (THM-590). The bridge from it to LRC(14) is the explicitly-open program (C).
