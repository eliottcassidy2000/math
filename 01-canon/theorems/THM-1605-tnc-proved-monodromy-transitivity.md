# THM-1605: THE TORAL NULLCONE CONJECTURE, PROVED — monodromy transitivity kills Pi = ct

**Status:** VERIFIED — ADVERSARIALLY REVIEWED (S176): the original step (3)
(Puiseux-DFT local lemma) contained an OVERCLAIM (MISTAKE-202: equal products
over disjoint subsets do NOT force equal subsets; e.g. I = {0,2} in Z_4 has
S_1 = S_3 = 0 automatically) — found in self-review BEFORE the fleet pass and
REPLACED by the simpler ORBIT-PRODUCT proof below, which two hostile referees
then confirmed (one pinhole — the c = 0 exclusion — patched; two exposition
debts paid: the Rouche cluster definition and the r0 != 0 irreducibility
role). Classical inputs: Gauss lemma/irreducibility => transitive monodromy;
permanence of identities under continuation; Vieta.
**Author:** boxeph-2026-07-20-S175 (HYP-8440)
**Consequence (fleet chain):** with klein's Gamma bridge (S351: TNC => NC2 =>
GMC(2)) and the assembled strata theorems: NC2, GMC(2), AND the 2-D NULLCONE
STRUCTURE THEOREM (nullcone = one-sided exactly) are COMPLETE.

## Statement
Let R be a polynomial with R(0) = r0 != 0, deg R = d = M + N, M >= 1, N >= 1.
Let u_1(t), ..., u_M(t) be the M roots of u^M = t R(u) tending to 0, and
Pi(t) = prod u_i(t). Then Pi(t) = c t exactly (some constant c) is IMPOSSIBLE.
Equivalently (THM-1550): CT(Lambda^m) = 0 for all m >= 1 with Lambda =
u^{-M} R(u) forces N = 0: the toral nullcone is one-sided. TNC holds.

## Proof
**(1) Fiber reformulation.** u^M = tR(u) iff H(u) = t where H(u) := u^M / R(u),
a rational map of degree d = M+N on P^1 (zero of order M at u = 0 since
r0 != 0; pole of order N at u = infinity since deg R = M + N). For small
t != 0 the fiber H^{-1}(t) consists of d simple points: M near 0 (the
0-cluster S_0) and N near infinity (the infinity-cluster). Pi(t) is the
product over the 0-cluster.

**(2) THE ORBIT-PRODUCT PROOF (replaces the original steps 2-4).**
Cluster definition (Rouche): choose delta, eps > 0 with no root of Phi_t on
|u| = delta for 0 < |t| < eps and no critical value of H in the punctured
disk; the root count in |u| < delta is continuous and integer, hence
constantly M: the 0-cluster S0 and Pi are single-valued analytic on the
punctured disk. Fix a basepoint t_b there.
PERMANENCE: the hypothesis gives p := prod_{v in S0} u_v == c t as germs at
t_b. Continuation along any loop gamma (in the base = punctured disk of
noncritical values, which also excludes 0) is a ring homomorphism on germs
sending p to prod_{v in gamma(S0)} u_v and fixing ct; hence EVERY monodromy
image S' = gamma(S0) satisfies prod_{v in S'} u_v == c t.
TRANSITIVITY: Phi = u^M - t R(u) is linear in t with gcd(u^M, R) = 1
(exactly r0 != 0), hence irreducible in C[u, t], hence irreducible over
C(t) (Gauss): the cover is connected and the monodromy group G is
transitive on the d labels. [If r0 = 0 the curve is reducible and the claim
is false — Pi == 0; r0 != 0 is essential and used exactly here.]
ORBIT COUNT: let O = {gamma(S0)} (r >= 1 distinct M-subsets, each of product
ct). Equivariance (v in S iff gamma v in gamma S) makes eta_v = #{S in O :
v in S} G-invariant, hence constant eta = rM/d by transitivity. Therefore
  prod_{S in O} prod_{v in S} u_v = (prod_{all v} u_v)^eta = ((-1)^d r0/r_d)^eta,
CONSTANT in t by Vieta (deg_u Phi_t = d exactly for t != 0 since r_d != 0 and
M < d; constant term -t r0: the t's cancel). But the left side equals (ct)^r.
C != 0 EXCLUSION (the referee's pinhole, two independent patches): (A)
Phi_t(0) = -t r0 != 0, so no fiber point vanishes and Pi(t) != 0 — Pi = 0*t
is impossible; (B) Puiseux: Pi(t) = (-1)^{M-1} r0 t (1 + o(1)), forcing
c = (-1)^{M-1} r0 != 0 (machine-verified: Pi/t -> -2 = -r0 on the test case).
Hence (ct)^r is nonconstant while equal to a constant. Contradiction. QED

**(2') [SUPERSEDED — kept for the record] Single-valuedness under continuation.** Suppose Pi(t) = ct. The germ
p(t) := prod_{u in S(t)} u (S(t) = the continued cluster) satisfies p = ct
near 0. Analytic continuation along any path in the complement of the finite
critical-value set preserves identities: after continuation along a loop
gamma, the continued function equals prod over the monodromy image S' =
gamma(S_0), and the identity gives prod_{S'} u = ct as germs at the SAME
basepoint. Hence for EVERY monodromy image S' of the 0-cluster:
prod_{u in S'} u(t) = ct, so for any two images: prod_{S \ S'} u = prod_{S' \ S} u
identically — equal products over DISJOINT subsets of fiber points.

**(3) The local lemma (Puiseux + discrete Fourier).** Let gamma be a small
loop around a critical value t* not in {0, infinity} (fiber points there
avoid 0 and infinity: H(0) = 0 and H(infinity) = 0 for N >= 1, so u = 0 or
infinity would force t* = 0). Its local monodromy is a product of disjoint
cycles; on a k-cycle the branches are w_i(t) = W(zeta^i mu), mu = (t-t*)^{1/k},
zeta = e^{2 pi i/k}, W(0) = w* in C^*, W'(0) != 0. Suppose the identity of (2)
forces prod_{i in I} W(zeta^i mu) = prod_{i in J} W(zeta^i mu) for disjoint
I, J within one cycle (ramification classes separate by Puiseux exponents, so
multi-cycle differences split per cycle). Taking log and expanding: the mu^j
coefficient is a triangular combination whose leading new term is a nonzero
multiple of sum_{i in I} zeta^{ij} - sum_{i in J} zeta^{ij}; inductively all
character sums of I and J agree (j = 0: |I| = |J|; j = 1..k-1 from the
expansion). Fourier inversion on Z_k gives indicator_I = indicator_J, so
I = J — contradicting disjointness unless I = J = empty. Hence NO local
monodromy element moves the cluster: gamma(S_0) = S_0 for every loop.
(Machine checks: no disjoint equal-character-sum pairs exist for k <= 8;
mixing monodromy demonstrably occurs for concrete R — see frozen out.)

**(4) Transitivity.** The cover H : P^1 -> P^1 is connected (P^1 is
irreducible), so its monodromy group is TRANSITIVE on the d fiber points.
But (3) shows the group preserves the proper nonempty subset S_0 (N >= 1
makes it proper). Contradiction. QED

## Cross-checks
- M = 2, all N: the independent INVOLUTION proof: u_2 = c u_1 / R(u_1)
  forces psi(x) = cx/R(x) to satisfy R(x) R(psi(x)) = c^2; the polynomial
  identity has LHS degree d^2 (top coefficient r0 r_d^d != 0, no cancellation
  possible against lower terms) vs RHS degree d^2 - d: impossible for d >= 2.
  (Degree-gap machine-checked d = 3, 4, 5.)
- All prior closures — N = 1 all M (S173), (2,2), (2,3) by hand, (2,4), (3,3)
  by gauged elimination, 25/25 Hensel deviations — are corollaries/corroboration.

## Consequences (the fleet chain, now complete)
TNC (this) => NC2 (klein's Gamma bridge, S351: factorial moment weights make
the growing-degree top term dominate the radial average) => GMC(2) and the
2-D NULLCONE STRUCTURE THEOREM (the GMC(2) nullcone = exactly the one-sided
polynomials; strata assembled across THM-1500..1580, all agents). The
Gaussian Moment Conjecture in two variables is resolved: TRUE, in maximal
contrast with GMC(N >= 3) FALSE — the charge-radius lock made quantitative.
Formalization scoping (from S173) now applies to a COMPLETE proof chain.
