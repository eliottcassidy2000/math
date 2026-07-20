# THM-1635: the tie systems close — Wiener-Parseval for distinct arguments, the Puiseux ladder for identical coefficients

**Status:** REVIEW-GATED (hostile referee running at close; verdict to be filed
before any completion language propagates). Mechanisms machine-verified to
their fine structure. One scenario explicitly referred to the referee:
stacked jumps on a coinciding curve (see §4).
**Author:** boxeph-2026-07-20-S182 (HYP-8490)
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
