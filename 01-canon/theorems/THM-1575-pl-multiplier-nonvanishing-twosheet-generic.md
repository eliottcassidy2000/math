# THM-1575: Picard-Lefschetz multiplier nonvanishing — the two-sheet and distinct-rate cases

**Status:** VERIFIED (two-sheet and pairwise-distinct-rate cases; resonant strata
remain OPEN and are precisely delimited)
**Author:** boxeph-2026-07-20-S172 (HYP-8400)
**Context:** By Global Rigidity (S171), a nullcone P satisfies A(s) = E[e^{-sP}] == 1
exactly. The Nullcone Structure Theorem needs: mixed P => A != 1.

## Statement
Let P in C[Z, W] be mixed (charges of both signs), and suppose the large-s
structure of A(s) on a maximal sector decomposes by steepest descent as
   A(s) = sum_j c_j(s) e^{rho_j(s)} + (flat part),
with finitely many nonconstant rate functions rho_j (algebraic in s, of
distinct leading behavior) and prefactors c_j algebraic times powers of s.
Then:
(a) [TWO-SHEET CASE] If there are at most two exponential sheets (rho, -rho),
    then A != 1. Hence no mixed P in two-sheet position is in the nullcone.
(b) [DISTINCT-RATE CASE] If the nonzero rates rho_j are pairwise distinct as
    germs at infinity, then A != 1.
(c) The remaining case — RESONANT STRATA, where two or more distinct saddles
    share the same exact rate function — is open and is the entire residual
    of the Nullcone Structure Theorem.

## Proof
(b) [(a) is a special case] Suppose A == 1. Then
   sum_j c_j(s) e^{rho_j(s)} = 1 - (flat part) =: F(s),
where F is flat-dominated (subexponential). The functions e^{rho_j} with
pairwise distinct nonconstant rates are linearly independent over the
differential field K of algebraic-with-s-powers prefactors: a nontrivial
relation sum c_j e^{rho_j} = F with c_j in K, F subexponential, evaluated
along a ray where one rate strictly dominates (exists by distinctness:
compare leading exponents/coefficients; pick the maximal rho_{j0} on a
subsector) gives |c_{j0} e^{rho_{j0}}| >> |all other terms| — contradiction
unless c_{j0} = 0; induct downward. So all c_j = 0.
But the c_j are LOCAL saddle data: for a nondegenerate saddle z_j of the
phase (s P(z, zbar) + |z|^2/2, in the rotated/Lefschetz-thimble presentation),
c_j = (Gaussian fluctuation determinant)^{-1/2} x (nonzero Jacobian factors),
which is NONZERO; degenerate saddles contribute Airy-type prefactors, also
not identically zero as germs. Contradiction. QED (a, b)

Citations: steepest-descent/thimble decomposition for exponential integrals
of polynomial phase (classical; Pham, Delabaere-Howls for the two-real-
dimensional oscillatory setting); linear independence of exponentials over
Hardy-type differential fields (classical).

## The exhibit and the checked instances
- P = Z + W (simplest mixed): A(s) = e^{2s^2} EXACTLY (moments C(2k,k)2^k k!,
  machine-verified): one visible sheet, multiplier 1 != 0.
- P = Z^2 + Z + W: A = 1 + 2s^2 + O(s^3) (E2 = 4): sheet at second order.
- The S171 parity fake 2X(1+iY): dies at m = 4; consistent — its two sheets
  are RESONANT through m = 3 only, then split.

## What remains (the precise residual)
Resonant strata: coefficient loci where distinct thimbles carry identical
rate functions to all orders. On such strata the independence argument fails
and cancellation c_{j1} e^{rho} + c_{j2} e^{rho} = (c_{j1}+c_{j2}) e^{rho}
could vanish. Residual conjecture (the last stone): on resonant strata the
SUM of same-rate multipliers is still nonzero for mixed P — equivalently the
vanishing-cycle classes over a common critical value never sum to zero in
homology with the Gaussian orientation. Attack routes: (i) monodromy — a
loop in coefficient space permutes same-rate thimbles, and a vanishing sum
is monodromy-invariant: compute the monodromy representation (irreducibility
would force the sum nonzero unless it is the invariant class, computable);
(ii) positivity at special rays (real-positive P-directions make same-rate
contributions positively aligned); (iii) exact elimination stratum by
stratum (S171 pattern — every resonant stratum met so far died by m <= 4).
