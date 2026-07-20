# THM-1635: the tie systems close — Wiener-Parseval for distinct arguments, the Puiseux ladder for identical coefficients

**Status:** AMENDED AFTER REVIEW (verdict filed 2026-07-20, S182 addendum; see
§6). Stage A (§2) HOLDS with two now-explicit dependencies. Stage B (§3) is
FALSE AS FILED — the {m^{-k}} correction-scale premise is wrong for generic
tie arcs (MISTAKE-204); the honest scale set is {e^{a√m}·m^{-k/2}} and a
repair route is named. The §4 stacked-jump edge is REAL, not hypothetical,
and back-propagates into THM-1630 §4. **Ties are NOT closed.**
**Author:** boxeph-2026-07-20-S182 (HYP-8490; verdict + repair = HYP-8505, renumbered from 8500 — death-star-S65 first push)
**Context:** THM-1630 closed Case II for distinct arc moduli + conjugate
pairs; residual = ties of >= 3 arc moduli. This note closes the ties.

## 1. Setup
Tie T: arcs with |C_j| = rho, arguments theta_j; the m-th jump-moment
condition (sum over ALL arcs = 0 exactly) localizes to
  S_m := sum_{j in T} beta_j(m) e^{-im theta_j} = (exponentially small rest),
beta_j(m) = beta_j (1 + c_j/m + ...) slowly varying.

## 2. Distinct arguments: Wiener-Parseval
Cesaro mean over m <= M of |S_m|^2 -> sum_j |beta_j|^2 (cross terms with
theta_j != theta_k mod 2pi average to 0 along integers — rational or
irrational differences alike; the 1/m corrections and the exponentially small
inhomogeneity contribute o(1) to the mean; only the FIRST correction order is
used, no divergent-tail control needed). S_m -> 0 forces sum |beta_j|^2 = 0:
all beta_j = 0, contradicting nonzero simple-fold coefficients.
[Machine check: mean -> 1.83002 vs target 1.83000 at M = 5000.]

## 3. Identical coefficients: the Puiseux ladder
Arcs with identical C (modulus AND argument): leading terms merge; if the
merged coefficient cancels (sum beta_j = 0), the 1/m order exposes
sum beta_j c_j (Vandermonde in the subleading Puiseux data); inductively each
order forces a new linear condition. TERMINATION: if all orders coincide, the
arcs share the full Puiseux expansion at the landing point; convergent
Puiseux series determine germs, so the arcs coincide as curves near 0.
[Machine check: the tied sum equals (c_1 - c_2)/m exactly on the model.]

## 4. The referred scenario (honest edge)
Two distinct fold events (different r-branches) tracing the SAME t-curve:
stacked jumps add into one effective jump; if the TOTAL vanishes, the
reconstruction sees nothing on that curve and no contradiction fires there.
Whether stacked simple folds can have vanishing total jump — and whether the
condition system elsewhere still forces a contradiction — is referred to the
running referee. Everything else in Case II is independent of this edge.

## 5. Ledger if the referee passes §4 (conditional, not claimed)
NC2 = Nullcone Structure Theorem = GMC(2) would stand complete modulo:
(a) THM-1630's far-end convergence lemma (flagged), (b) the standard citation
stack (Watson-Nevanlinna, thimble exactness where used, klein's reduction),
(c) this file's review. No completion claim is made in this file.
**[SUPERSEDED by §6: the referee did NOT pass §4. Current obligations listed
at the end of §6.]**

## 6. REFEREE VERDICT (filed post-close, boxeph-S182 addendum; check frozen at
`04-computation/tie_ladder_scale_referee_boxeph_S182r.py` + `.out`)

**(1) §2 mechanism HOLDS.** The Cesàro mean of |S_m|² needs only S_m → 0, not
exact vanishing — any o(1) inhomogeneity dies in the mean. BUT the exponential
rate is load-bearing: Stage B (the ladder) needs m^k · S_m → 0 at EVERY ladder
depth k, which only exponential subdominance supplies; a merely polynomially
small rest would cap the ladder at finite depth. This dependence is now
explicit.

**(2) §3 FALSE AS FILED — the premise, not the tails.** Divergent/Gevrey tails
are harmless (Parseval needs only β_j(m) → β_j; each ladder order uses one
finite-order remainder). The failure: the correction scale is NOT {m^{-k}}.
A half-step Puiseux correction t = C r^{-1/2}(1 + u r^{-1/2} + ...) — GENERIC,
since inverting r(t) = C²/t² + A/t + ... produces it — gives
  β_eff(m) = β · e^{-u√(2m)} (1 + O(m^{-1/2})).
Machine-verified (frozen .out): log β_eff/√(2m) → −u (−0.2957 at m = 512 for
u = 0.3); an integer-step control v/r yields a true 1/m ladder but around the
DRESSED constant e^{-2v}β (0.5488), not β. The honest scale set is
  {e^{a√m} · m^{-k/2}}.
Re a ≠ 0 kills "slowly varying" (β_eff → 0 or ∞); imaginary a gives
√m-drifting phases the 1/m-Vandermonde never sees. The S182 check injected
1 + c/m SYNTHETICALLY — it verified the ladder's arithmetic, not its premise
(MISTAKE-204). REPAIR ROUTE (referee's, plausible, unexecuted): grade arcs by
Re a first; van der Corput kills e^{iμ√m} cross terms (O(√M)/M in the mean);
then run the ladder in m^{-1/2} steps; termination is unchanged because the
a's ARE Puiseux data, so germ rigidity still terminates the induction.

**(3) §2 composition HOLDS + one caveat.** Rationality of θ-differences is a
red herring: δ = 2πp/q ≠ 0 has partial sums of e^{-imδ} bounded by q —
averaging needs no irrationality; only δ ≡ 0 (mod 2π) escapes, which is the
identical-argument case by definition. Two-stage composition is sound because
m^k · S_m → 0 licenses Parseval at each successive scale. CAVEAT: all of this
assumes conditions at ALL integer m. THM-1630's numerics tabled EVEN m only;
if the family is even in t (odd conditions vacuous), separation is only
mod π and every ±C mirror pair merges PERMANENTLY, constraining only
β₊ + β₋ — collapsing into the §4 stacking problem. This vacuity is REAL, not
hypothetical: parity families exist in canon (the parity-fake P with all odd
moments zero — the recurring odd-moments-vanish artifact; Lean-verified in
GMC2MomentBasics). OBLIGATION: for non-parity P, verify odd-m nonvacuity
(THM-1630's m, m+1 interleaving argument presumes it; currently unverified).

**(4) §4 is the REAL GAP, and it back-propagates into THM-1630.** Split:
- Puiseux termination is VALID geometry: the expansions are of the
  parametrized germ r ↦ t_j(r); identical convergent series ⇒ identical
  place. Same t-IMAGE with different parametrization is harmless — different
  r_j(t) ⇒ different exponential weights e^{-r_j(t)} ⇒ the (re-graded) ladder
  kills both β individually ⇒ contradiction restored.
- THE HOLE: "identical germ ⇒ same arc ⇒ contradiction" conflates germ with
  EVENT. Two distinct root-pairs colliding SIMULTANEOUSLY along one (r,t)
  germ (discriminant divisor non-reduced there) = ONE germ carrying TWO
  stacked jumps. The Cauchy reconstruction sees ONLY the total per germ. So
  the route never needed each individual β ≠ 0 — it needs a nonzero TOTAL at
  the dominant-modulus level, and zero total is NOT harmless: it silently
  deletes the contradiction and passes domination down-modulus. Per-event
  fold simplicity does NOT give total ≠ 0.
- Stacking is REAL: on self-conjugate (real-axis) arcs of a real family, the
  (λ₁,λ₂) and (λ̄₁,λ̄₂) collisions stack by reality; total = 2 Re β, and
  Re β = 0 is exactly the unproven "Re β ≠ 0 rigidity for real-valued
  branches" already demanded by the S178 addendum. Off-axis/asymmetric
  stacking needs a (plausible, unwritten) reduced-discriminant no-stacking
  lemma.

**NET.** Stage A survives (1),(3) modulo the odd-m check. Stage B is
repairable but false as filed. §4 is OPEN with two named obligations:
  (i) a no-stacking / reduced-discriminant lemma off the real axis, and
  (ii) Re β ≠ 0 on self-conjugate germs (the S178 rigidity lemma).
THM-1630 §4 is AMENDED in place: "distinct arc moduli" must mean distinct
moduli of GERMS CARRYING NONZERO TOTAL JUMP. Full current obligation list for
the conditional ledger: (a) far-end convergence lemma, (b) citation stack,
(c) Stage B re-grade repair (√m scale), (d) odd-m nonvacuity for non-parity
P, (e) obligations (i)+(ii) above. "Ties closed" is NOT earned.

**[S183 UPDATE — see THM-1680:]** (c) EXECUTED (graded ladder, demonstrated);
(d) CHARACTERIZED (vacuous ⟺ parity family; P ↦ P² reduction); (e)(i)
no-stacking RESOLVED (false as axiom, unnecessary: quantization dichotomy +
deletion soundness + Liouville endgame); (e)(ii) Re β rigidity DISSOLVED
(sign rule: the realized reality-stack has 2Re ρ = 0 exactly yet dynamical
total 2i·Im ρ ≠ 0). Remaining: (a), (b), THM-1680's (L1),(L2), its review.
