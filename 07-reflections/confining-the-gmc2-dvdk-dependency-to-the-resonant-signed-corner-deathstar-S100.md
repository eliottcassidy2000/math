# Confining the GMC(2) DvdK dependency: it is elementary except in the resonant-signed corner

**death-star-2026-07-21-S100** (HYP-8877). Owner: creative math (not Lean) — sharpen targets / bypass the
GMC(2) dependency on Duistermaat–van der Kallen. Result: DvdK is **not** needed for generic GMC(2); its one hard
use confines to a single corner. This does not eliminate DvdK (the corner is genuine one-variable DvdK), but it
isolates it precisely and shows two large regimes are DvdK-free by elementary means.

## The setup
GMC(2)/THM-2022 needs, for the lowest balanced face `f_F = Σ_{i∈F} c_i u^{q_i}` (distinct charges `q_i`
straddling 0), that some power has a nonzero constant term: `∃ m0, CT(f_F^{m0}) ≠ 0`. That is the DvdK input.
Crucially, `E[P^m] = Σ_{balanced channels r} (multinomial(m,r)·A(r)!) · c^r` has **positive integer weights**;
the *only* source of cancellation — hence the *only* reason DvdK is nontrivial — is the sign/phase of the
coefficient monomials `c^r`.

## Two independent elementary reductions (verified)
**Reduction 1 — the lowest face is generically an EDGE (2 monomials), where DvdK is a binomial.**
The face comes from the LP `min Σ a_i x_i` s.t. `x ≥ 0, Σx=1, Σq_i x_i = 0` — **two** equality constraints
(mass, charge). A basic optimum therefore has `≤ 2` nonzero `x_i`, i.e. is supported on an **edge** `{u^{-a}, u^{b}}`
(`a,b>0`, straddling), and generically the tilted-height-`δ` set is exactly those two monomials (the rest strictly
higher). For an edge, DvdK is elementary: with `g=gcd(a,b)`, `m0=(a+b)/g`, `k=b/g`,
```
CT(f_edge^{m0}) = C(m0, k) · c_-^k · c_+^{m0-k} ≠ 0.
```
Verified for `(-a,b) ∈ {(-1,1),(-1,2),(-2,3),(-3,5),(-2,7),(-5,4)}` — first nonzero CT at `m0=(a+b)/g`, always a
nonzero binomial. A `≥3`-monomial face needs `≥3` tilted heights `a_i−λq_i` concurrent at the straddling `λ*` — a
**codimension-≥1 resonance** on the Newton support.

**Reduction 2 — positive real coefficients give `CT(f^m) > 0` for free.**
Since the channel weights are positive, `c_i > 0 ⟹` every term of `CT(f^m)` is positive, so `CT(f^m) > 0` with no
cancellation (a nonnegative walk-return / recurrence). Verified: charges `(-1,0,1)` with `(1,1,1)` give the
central trinomial `1,3,7,19,51,141,393 > 0` (A002426), and `(2,1,3)` gives `1,13,37,289,… > 0`. So the DvdK input
is automatic whenever the descent point can be taken positive-real.

## The hard corner (where DvdK genuinely lives)
DvdK is nontrivial **iff** the lowest face is a **resonance** (`≥3` distinct charges, a degenerate/non-generic
support) **and** the coefficients are **signed/complex with cancellation**. Verified example: charges `(-2,-1,1,2)`
with signed coeffs `(1,1,-1,-1)` gives `CT(f^{1..8}) = 0,-4,0,36,0,-370,0,4004` — low powers vanish, DvdK
guarantees the eventual nonzero (here at `m=2`, consistent with ESV `≤ m+n`). This corner is exactly the
S89–S91 charge-**resonance** / confluent case, whose tied-core weights are the **central trinomial = a
free-probability moment** (S90) and whose generating function `W(t)=CT(1/(1−t f_F))` is the algebraic
transfer/zeta resolvent of Monsky's proof (the tournament-zeta lens, S99). So the residual hard input is precisely
"the resonant multi-clock face" in the S99 scale-then-clock picture — the edge case is a single 2-clock (period
`(a+b)/g`), the resonance is the genuine multi-clock.

## Honest scope and value
- **Not a full bypass.** The resonant-signed corner is genuine one-variable DvdK and still requires it (or Monsky,
  ≈person-months per the S95 roadmap). No claim that GMC(2) is unconditionally DvdK-free.
- **A verified confinement / sharpening.** GMC(2) is DvdK-free (elementary binomial or positivity) for (i) every
  support whose lowest balanced face is an edge — the generic case — and (ii) every point reducible to positive-real
  coefficients. The hard input lives only in the resonant ∩ signed corner, which is codimension-≥1 in support space
  and is exactly the central-trinomial/transfer-operator object already studied (S90, S99).
- **Formalization payoff.** A Lean GMC(2) could discharge the generic edge case elementarily (binomial CT) and cite
  DvdK only for the resonant face — splitting the one non-Mathlib input off a codimension-≥1 stratum rather than the
  whole problem.

Cross-links: THM-2022 (the GMC2 proof), THM-1630 (DvdK), S90 (central trinomial = tied-core weights),
S99/HYP-8876 (scale-then-clock; edge = 2-clock, resonance = multi-clock), S95 (DvdK roadmap), memory
`nc2-gmc2-lean-formalization-state`. Script `04-computation/dvdk_confinement_deathstar_S100.py` (+ `.out`). HYP-8877.
