# GVC(3) claimed counterexample (Lambda = Delta^6): exact verification of the finite instances

- **Agent**: boxeph (subagent, multifront worktree)
- **Date**: 2026-08-03
- **Script**: `04-computation/gvc3_delta6_counterexample_verification_boxeph.py`
- **Output**: `05-knowledge/results/gvc3_delta6_counterexample_verification_boxeph.out` (30 PASS, 0 FAIL, exit 0, wall 59.6 s, sympy 1.9, exact integer/rational arithmetic only)
- **Provenance**: UNCONFIRMED. The user-supplied citation (arXiv 2606.17854) resolves to an unrelated paper. Everything below is judged on the arithmetic alone.

## Setup

Variables x, y, t. Operator **Delta f = 4 f_xy + f_tt**.

- rho = t^2 + x*y
- A = rho + x^2
- C = y*rho^2 - 2*x*t^2*rho - x^3*t^2   (homogeneous, degree 5)
- P = A*C^2   (homogeneous, degree 12)
- Q = x^2

Claimed counterexample to the Generalized Vanishing Conjecture GVC(3) for
Lambda = Delta^6: Lambda^m(P^m) = 0 for all m >= 1, yet Lambda^m(Q*P^m) != 0
for all m >= 1 — contradicting GVC's conclusion "Lambda^m(Q*P^m) = 0 for m >> 0".

## Verified — FINITE-EXACT

All of the following were checked by full expansion in sympy (no floats, no numpy):

1. **deg(P) = 12** and **P has exactly 23 terms** fully expanded; P is homogeneous. — FINITE-EXACT
2. **Key identity**: `x*C == rho^3 - t^2*(rho + x^2)^2` exactly. Consequence (also checked): `Q*P = x^2*A*C^2 = A*(rho^3 - t^2*A^2)^2`. — FINITE-EXACT
3. **Vanishing**: Delta^6(P) == 0; Delta^12(P^2) == 0; Delta^18(P^3) == 0.
   Timings: 1.2 s / 9.1 s / 28.9 s (P^3 is degree 36, 168 terms — cheap). — FINITE-EXACT for m = 1, 2, 3
4. **Nonvanishing with exact constants**, for m = 1, 2: since Q*P^m is homogeneous of degree 12m+2 (checked) and each Delta drops degree by exactly 2, Delta^(6m+1)(Q*P^m) is a constant (checked: no free symbols). Exact values:
   - m = 1: Delta^7(Q*P) = **697 426 329 600** = 2^9 * 7! * 2! * 15!!/5!!  (matches claimed formula)
   - m = 2: Delta^13(Q*P^2) = **4 424 683 459 217 616 116 121 600 000** = 2^17 * 13! * 4! * 27!!/9!!  (matches)
   - m = 0 edge case **FAILS as claimed**: Delta(x^2) = 0 while the formula would give 2*3!! = 6; the closed form is for m >= 1 only. — FINITE-EXACT for m = 1, 2
5. **Consistency**: Delta^(6m)(Q*P^m) != 0 for m = 1, 2 (as forced by 4: one more Delta of zero is zero). Explicitly Delta^6(Q*P) = 2372198400*(15*t^2 + 17*x^2 + 66*x*y). — FINITE-EXACT

## CITED-UNVERIFIED — what remains for a genuine GVC(3) refutation

The refutation of GVC(3) needs two **all-m** statements, of which only finite instances are verified here:

- **(V)** Delta^(6m)(P^m) = 0 for **all** m >= 1. Verified m = 1, 2, 3 only. The claimed source's induction (unlocatable; cited arXiv id resolves to an unrelated paper) is required for all m. — CITED-UNVERIFIED
- **(N)** Delta^(6m+1)(Q*P^m) = 2^(8m+1) * (6m+1)! * (2m)! * (12m+3)!! / (4m+1)!! (hence != 0) for **all** m >= 1. Verified m = 1, 2 only. This is the claimed closed form to be recorded; nonvanishing for infinitely many m is what contradicts "for m >> 0". — CITED-UNVERIFIED

Conditional conclusion: **if** (V) and (N) hold for all m (finite instances above are consistent and exact), then GVC(3) fails for Lambda = Delta^6, P = A*C^2, Q = x^2. The finite checks alone refute only a vanishing threshold m0 <= 2 and prove Lambda is "admissible-so-far" up to m = 3. Status of the counterexample as a whole: **OPEN pending (V), (N)**; every individually checked instance: **FINITE-EXACT**.

## Structural observations (mechanism probes)

- **Delta is the Laplacian of rho.** The Gram matrix of q = rho = t^2 + xy in (x,y,t) is M = [[0,1/2,0],[1/2,0,0],[0,0,1]]; the dual form q*(xi) = xi^T M^{-1} xi = 4 xi_x xi_y + xi_t^2 is exactly the symbol of Delta (checked). So rho is the Fourier dual of the symbol, and Delta = q*(d) is the Laplace operator of the nondegenerate (signature (2,1)) quadratic form rho.
- **Radial law** (n = 3): Delta(rho^k) = 2k(2k+1) rho^(k-1) (checked k = 1..4), hence Delta^k(rho^k) = 2^k * k! * (2k+1)!! (checked k = 6, 7, 13).
- **Fischer/apolarity reading.** For homogeneous f of degree 2k, Delta^k(f) is (up to the normalization Delta^k(rho^k)) the Fischer pairing <f, rho^k>. So claim (V) says **P^m is Fischer-orthogonal to rho^(6m)**, and claim (N) says <Q*P^m, rho^(6m+1)> != 0; in normalized form the claim-4 constant divided by Delta^(6m+1)(rho^(6m+1)) is exactly **2^(2m) (2m)! / (4m+1)!!** (= 8/15 at m=1, 128/315 at m=2; checked). Any all-m proof plausibly computes this normalized pairing.
- **The mechanism is NOT Hessian nilpotency of C**: Delta(C) != 0 (C is not harmonic), q*(grad C) != 0 (gradient not isotropic), Delta(C^2), Delta(C^3) != 0. Also Delta(P) != 0 and Delta(A) = 6.
- **Sharpness**: Delta^5(P) = -24710400*(24 t^2 + 8 x^2 - 12 x y - 21 y^2) != 0 and Delta^11(P^2) != 0, so the Delta-vanishing order of P^m is exactly 6m at m = 1, 2 — the exponent 6 in Lambda = Delta^6 is tight, not slack.
- **Where the induction likely lives**: by the key identity, x*C = rho^3 - t^2*A^2, so Q*P = A*(rho^3 - t^2*A^2)^2 lies in the subring generated by rho, A, t (equivalently rho, x^2, t^2); more generally x^(2m)*P^m = A^m*(rho^3 - t^2*A^2)^(2m). The double factorials in the constant are the fingerprints of the radial law Delta(rho^k) = 2k(2k+1) rho^(k-1) applied down a rho-power ladder. Also q*(grad A) = 4*(A + x^2) (checked), so A generates a closed little system under Delta and the dual pairing.
- **m = 0 anomaly** is structural: the formula's m >= 1 restriction reflects that Q = x^2 alone is Delta-harmonic (x^2 contains no y, t^2-mixables), and only multiplication by P^m couples it to rho-powers.

## Reproduction

```
cd <worktree>
python3 04-computation/gvc3_delta6_counterexample_verification_boxeph.py
```

Exit 0 iff all claimed-instance checks pass; exploratory probes are reported as `[PROBE] ...: YES/NO` and never fail the run.
