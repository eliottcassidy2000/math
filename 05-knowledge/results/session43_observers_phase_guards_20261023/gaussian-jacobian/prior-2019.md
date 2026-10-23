---
id: THM-2019
title: "NC2 for affine charge-height supports and their common radial multiples"
status: >
  PROVED.  On one affine line in (charge,total degree), every balanced
  length-m word has the same Wick factorial.  For a common nonzero radial
  multiplier B(ZW), the moment factors exactly into a toral constant term and
  one ordinary exponential moment.  Duistermaat--van der Kallen plus EMP then
  rule out every two-sided or charge-zero nullcone point.  The proof includes
  the rational-intercept/integrality subsequence and arbitrary complex
  coefficients.  Exact independent verification is stored for two genuinely
  nonhomogeneous weighted supports.
source: codex-2026-07-21-NC2-followup
depends_on:
  - THM-1630  # DvdK's one-variable Laurent nullcone theorem
  - THM-1510  # EMP: eventual nonvanishing of Laplace moments of a polynomial
  - THM-1645  # charge/radial moment dictionary
related:
  - THM-2014
  - THM-2017
  - THM-2022
  - HYP-8765
script: 04-computation/gmc2_affine_height_supports_codex_20260721.py
---

# THM-2019 — the affine height lock

Put (s=ZW), and let

\[
 Q(Z,W)=\sum_{i=1}^k c_i Z^{a_i}W^{b_i}\ne0,
 \qquad q_i=a_i-b_i,
 \qquad h_i=a_i+b_i .                         \tag{1}
\]

Assume that the support obeys one affine charge-height law

\[
 h_i=\lambda q_i+\delta\qquad(1\le i\le k).   \tag{2}
\]

Equivalently, (Q) is weighted-homogeneous for the possibly mixed weights
(1-\lambda,1+\lambda):

\[
 (1-\lambda)a_i+(1+\lambda)b_i=\delta .       \tag{3}
\]

Let (B\in\mathbb C[s]\setminus\{0\}) be arbitrary and set

\[
 P(Z,W)=B(ZW)Q(Z,W).                           \tag{4}
\]

> **Theorem.** If (mathbb E[P^m]=0) for every (m\ge1), then all
> charges occurring in (Q), and hence in (P), are strictly positive or
> all are strictly negative.  Conversely either strict one-sided condition
> makes every moment vanish.  Thus NC2 holds for every affine-height support
> and for every one of its common nonzero radial multiples.

For each fixed nonzero (B), its coefficients and the coefficients (c_i) are
arbitrary complex numbers. The actual support is taken after deleting zero
(c_i); such deletion merely passes to a subset of the same affine line.

## 1. Exact separation of the toral and radial addresses

Define the Laurent polynomial retaining only the charge address,

\[
 A(u)=\sum_i c_i u^{q_i}.                     \tag{5}
\]

There is no collision hidden in (5): on the affine line (2), (q_i)
determines (h_i), and ((q_i,h_i)) determines ((a_i,b_i)).  Equal
monomials are combined before writing (1).

Consider a word with multiplicity vector (r=(r_1,\ldots,r_k)),
(sum_i r_i=m).  If it is charge-balanced, then

\[
 \sum_i q_i r_i=0,
 \qquad
 \sum_i h_i r_i=\lambda\sum_iq_ir_i+\delta m=\delta m. \tag{6}
\]

Its total (Z)- and (W)-degrees are therefore both

\[
 N_m={\delta m\over2}.                         \tag{7}
\]

In particular, whenever a balanced word of length (m) exists, (N_m) is
automatically a nonnegative integer.  Every such word has the *same* Wick
factor (N_m!).  Consequently the charge-zero part and the Gaussian moment
factor exactly:

\[
 \operatorname{CT}_{\rm charge}(Q^m)
   =s^{N_m}\operatorname{CT}_u(A^m),
 \qquad
 \boxed{\;
 \mathbb E[P^m]
   =\operatorname{CT}_u(A^m)\,
      L\!\left(s^{N_m}B(s)^m\right),
 \;}                                                       \tag{8}
\]

where (L(s^n)=n!). Formula (8) is asserted at reachable lengths, for which
(N_m) is an integer. If no balanced word of length (m) exists, then
(
\operatorname{CT}_u(A^m)=\mathbb E[P^m]=0
\)
directly, and no fractional power of (s) is formed. This is the
missing-address repair in its cleanest form: charge decides which words
survive, while (2) makes their height address constant.

## 2. The rational-intercept subsequence

Suppose first that (A) has a positive and a negative exponent.  Choosing
support points (q_+>0>q_-), equation (2) gives

\[
 \delta={q_+h_-+|q_-|h_+\over q_++|q_-|}>0.    \tag{9}
\]

Thus (delta\in\mathbb Q_{>0}).  Choose an integer (ell\ge1) for which

\[
 D={\ell\delta\over2}\in\mathbb Z_{>0},
 \qquad
 H(s)=s^D B(s)^\ell\ne0.                      \tag{10}
\]

At (m=\ell n), (8) becomes the exact product

\[
 \mathbb E[P^{\ell n}]
   =\operatorname{CT}_u\!\left((A^\ell)^n\right)L(H^n).     \tag{11}
\]

The Laurent polynomial (A^\ell) is still two-sided: its lowest and
highest exponents are (ell q_{\min}<0<ell q_{\max}), with nonzero extreme
coefficients.  Duistermaat--van der Kallen's Theorem 2 (THM-1630) gives

\[
 \limsup_{n\to\infty}
 \left|\operatorname{CT}_u((A^\ell)^n)\right|^{1/n}>0,       \tag{12}
\]

so the toral factor is nonzero for arbitrarily large (n).

The exponential moment theorem (EMP, THM-1510) says that for every nonzero
polynomial (H), (L(H^n)\ne0) for all sufficiently large (n).  (For
positive degree this is its leading-factorial asymptotic; for a nonzero
constant it is immediate.)  Choose a sufficiently large (n) from the
infinite set supplied by (12).  Both factors in (11) are nonzero, contrary
to the nullcone assumption.  Hence a nullcone member cannot straddle charge
zero.

This subsequence is essential bookkeeping, not decoration.  The intercept
(delta) need not be an integer; only (ell\delta/2) must be integral.

## 3. Neutral charge is also impossible

It remains to exclude support contained in (q\ge0) (or (q\le0)) but
containing (q=0).  On (2), the charge-zero point has (h=\delta=2a), so
(delta/2\in\mathbb Z_{\ge0}).  Let (c_0\ne0) be its coefficient.  If
all other charges are nonnegative, the only charge-zero words use that term,
and (8) reads

\[
 \mathbb E[P^m]
   =c_0^m L\!\left((s^{\delta/2}B(s))^m\right).              \tag{13}
\]

EMP makes (13) nonzero for all sufficiently large (m).  The nonpositive
case is identical.  Therefore the only nullcone supports are strictly
one-sided.  Conversely, a strict one-sided support has no balanced word in
any positive power, proving the reverse implication and the theorem.  \(\square\)

## 4. Pair radicals on this many-circuit stratum

Fix (B\ne0) and an affine-height support (S), and let

\[
 I_S=\langle\mathbb E[P],\mathbb E[P^2],\ldots\rangle
       \subset\mathbb C[c_1,\ldots,c_k].                    \tag{14}
\]

The theorem identifies (V(I_S)) exactly with the union of the two strict
one-sided coordinate subspaces.  Hence the Nullstellensatz gives

\[
 c_i c_j\in\sqrt{I_S}\quad(q_i>0>q_j),
 \qquad
 c_t\in\sqrt{I_S}\quad(q_t=0).                              \tag{15}
\]

Thus HYP-8765's desired pair-radical conclusion is proved here even with
arbitrarily many primitive circuits and arbitrary complex phases.  No
effective claim such as the conjectural ((k-1)R(S)) cutoff is made;
Noetherianity only guarantees some finite set of moment levels.

## 5. Strictly beyond ordinary homogeneity

Ordinary homogeneous polynomials are the special case (lambda=0),
(B=1).  The theorem is substantially larger.  For example,

\[
 Q=aZ^6+bZ^4W+cZ^2W^2+dW^3                       \tag{16}
\]

lies on (a_{\rm exp}+2b_{\rm exp}=6), not on a constant-total-degree
shell.  Its charges and heights are

\[
 (q,h)=(6,6),(3,5),(0,4),(-3,3),
 \qquad h={q\over3}+4.                            \tag{17}
\]

For every nonzero (B(s)),

\[
 \mathbb E[(B(s)Q)^m]
 =\operatorname{CT}(au^6+bu^3+c+du^{-3})^m\,
   L(s^{2m}B(s)^m).                               \tag{18}
\]

A rational-intercept example is the weighted-degree-eight support
(a_{\rm exp}+2b_{\rm exp}=8).  Its five points have

\[
 q=(8,5,2,-1,-4),qquad h={q\over3}+{16\over3}.    \tag{19}
\]

All balanced lengths are multiples of three, and with (ell=3,D=8),

\[
 \mathbb E[(B(s)Q)^{3n}]
 =\operatorname{CT}(A^{3n})L\!\left((s^8B(s)^3)^n\right).  \tag{20}
\]

The exact verification script checks (18) directly through (m=8) with
(B=1-2s+3s^2), and checks (20), including all zero off-subsequence
moments, through (m=6).

## 6. What this says about the remaining many-circuit problem

For a chosen shear (lambda), attach to each support monomial its affine
intercept

\[
 \Delta_i=h_i-\lambda q_i.                       \tag{21}
\]

The theorem closes **affine-height rank one**, and also the special
parallel-layer stacks produced by one common radial factor (B(s)).  This
is the precise transfer from the repository's address-coordinate work:
once all channels share one radial address, the scalar moment splits into a
toral address and an EMP address.

What remains in HYP-8765 is not merely “more affine lines.”  It is the case
where different charge sectors carry *incompatible* radial factors, so no
common (B(s)^m) can be extracted.  There the factorial-Hankel/resultant
tower is genuinely needed to reconstruct the lost layer address.  THM-2019
therefore supplies both a broad closure and a sharp diagnostic for the
residual.

## Novelty and scope

The ordinary homogeneous case is classical and already covered by the
published Gaussian-moment results.  The weighted-homogeneous affine-line
factorization, especially with an arbitrary common radial multiplier, was
not previously stated in this repository; no claim of external literature
novelty is made.  The proof is an exact composition of known DvdK and EMP
theorems.  It does **not** itself cover independent radial polynomials attached
to different charges. Those cases, and full NC2, are now proved separately by
THM-2022's lowest-face Frobenius congruence; the present factorization remains
a sharper description on its stratum.
