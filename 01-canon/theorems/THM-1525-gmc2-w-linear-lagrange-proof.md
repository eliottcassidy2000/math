# THM-1525: GMC(2) holds on the W-linear class — the Lagrange-inversion proof

**Status:** VERIFIED (complete proof; master identity machine-checked, each branch elementary)
**Author:** boxeph-2026-07-20-S167 (HYP-8360)
**Context:** GMC(N) false for N >= 3 (owner-supplied counterexample, verified S165);
GMC(1) and homogeneous-GMC(2) known; unrestricted GMC(2) open. This settles the
class P = A(Z) + W B(Z) (equivalently Z-degree <= 1 by symmetry), strictly
containing the S166 radial family.

## Statement
Let Z = X+iY, W = X-iY for independent standard real Gaussians X, Y, and let
P = A(Z) + W B(Z) with A, B in C[s]. If E[P^m] = 0 for all m >= 1, then either
P = beta W (B constant, A = 0) or P is a pure positive-charge polynomial
(A(0) = 0 and val B >= 2). In all cases, for every fixed Q in C[Z, W],
E[Q P^m] = 0 for all sufficiently large m. GMC(2) holds on this class.

## Proof
**(0) Master identity (directed Wick => unicyclic vacuum diagrams; verified
coefficient-wise to m = 7 against brute Wick on random instances):**
  E[e^{tP}] = e^{t A(sigma)} / (1 - 2t B'(sigma)),   sigma = 2t B(sigma),
and the insertion law E[Z^alpha e^{tP}] = sigma(t)^alpha * E[e^{tP}].
(W pairs only into Z: each vertex has at most one out-arrow, so connected
vacuum diagrams are trees + one cycle; trees give A(sigma), the cycle gives
the log-det denominator. Exact, not semiclassical.)

**(1) Branch B(0) != 0.** sigma is an invertible change of parameter;
nullcone <=> exp(sigma A / 2B) = (B - sigma B')/B in C(sigma)-series. The
right side is rational; the classical logarithmic-derivative lemma (y'/y of a
rational function has only simple poles with integer residues — never a
nonzero polynomial or rational with the wrong poles) forces A = 0. Then
nullcone <=> [s^m] B(s)^m = 0 for all m, and induction (b_1 = 0, then
[s^m]B^m = m beta^{m-1} b_m) forces B = const. P = beta W: safe (E[Q W^m] = 0
once m exceeds the Z-degree of Q).

**(2) Branch val B = 1.** Valuation-exact extraction: [s^k](B^k A^{m-k}) =
b_1^k A(0)^{m-k}, so E[P^m] = sum_k C(m,k) 2^k k! b_1^k A(0)^{m-k} depends
only on (A(0), b_1); m = 1, 2 give a_0 = -2 b_1 and then 4 b_1^2 = 0. The
nullcone is EMPTY here.

**(3) Branch val B >= 2 (or B = 0).** [s^k](B^k A^{m-k}) = 0 for k >= 1, so
E[P^m] = A(0)^m and the nullcone forces A(0) = 0. Then every monomial of P
has U(1)-charge >= +1 (Z-terms: charge >= 1; W Z^j with j >= 2: charge j-1 >= 1),
so P^m has charge >= m and E[Q P^m] = 0 as soon as m exceeds the maximal
negative charge of Q. Safe by charge positivity. QED

## Why N = 3 escapes and where the wall is
The third Gaussian twists the Lagrange sum by E[U^k] = (1/2)_k — the
fiber-fraction weights — enabling the (1+x)^{-1/2} collapse (the S165
counterexample). In N = 2, W-degree >= 2 vertices have several out-arrows:
vacuum diagrams of arbitrary loop order, generically DIVERGENT formal series
(e.g. E[e^{t Z^2 W^2}] has coefficients 4^m (2m)!/m!). The remaining wall for
full GMC(2) is exactly this resurgent regime; the W-linear case is the
complete convergent-unicyclic stratum.

Scripts: 04-computation/gmc2_wlinear_theorem_boxeph_S167.py (+ frozen out).
