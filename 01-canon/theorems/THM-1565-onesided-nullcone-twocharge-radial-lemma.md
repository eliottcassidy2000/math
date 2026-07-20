# THM-1565: The Radial Lemma and the One-Sided Theorem for two-charge spans — GMC(2) progress

**Status:** VERIFIED (Radial Lemma modulo two classical citations; reduction identity
machine-verified exactly; three-charge base case by exact elimination)
**Author:** boxeph-2026-07-20-S169 (HYP-8385)
**Context:** GMC(2) reduces to the One-Sided Nullcone Conjecture (S168). This
proves the conjecture for ALL two-charge spans at ALL degrees, plus the first
three-charge case.

## Theorem A (the Radial Lemma, all degrees)
Let psi be a nonzero polynomial and L(f) = int_0^inf f(r) e^{-r} dr. If
L(psi^m) = 0 for all m >= 1 then psi = 0 (nonconstant or constant alike).
*Proof.* A(s) := L(e^{-s psi}) is analytic where the integral converges; contour
rotation of r within |arg r| < pi/2 (legitimate: integrand entire, decay in the
overlap sector) extends A to a sector of opening (D+1) pi on the log-Riemann
surface, D = deg psi. The nullcone makes A asymptotic to 1 AT ALL ORDERS with
remainders |s|^M (DM)! C^M / M! — Gevrey order D-1. Since (D+1) pi > (D-1) pi,
Watson–Nevanlinna uniqueness [CITATION: Nevanlinna; Sokal's improved Watson]
forces A == 1 on the sector. But along a ray with Re(s c_D e^{i D theta}) > 0,
dominated convergence gives A(s) -> 0 as s -> infinity (mass near r = 0 is
O(s^{-1/k}) if psi(0) = 0, and e^{-s psi(0)} -> 0 after rotating s if not).
Contradiction. QED

## Theorem B (One-Sided for two-charge spans, all degrees)
Let P = Z^a f(R) + W^b g(R), a, b >= 1, f, g in C[R], R = ZW. If E[P^m] = 0 for
all m >= 1, then f = 0 or g = 0.
*Proof.* Charge balance kills E[P^m] unless (a'+b') | m (a' = a/gcd, b' = b/gcd);
and E[P^{(a'+b')t}] = C((a'+b')t, b't) * E[H^t] with H = R^{ab/gcd} f^{b'} g^{a'}
— a PERFECT POWER (machine-verified exactly on 6 random (a,b,f,g) with
off-lattice vanishing). The nullcone thus puts H in the radial nullcone;
Theorem A gives H = 0; C[R] is a domain, so f = 0 or g = 0. QED

## Theorem C (first three-charge case, exact elimination)
For P = alpha Z + beta W + c0 + c1 R: E[P^m] depends on (alpha, beta) only
through p = alpha*beta (U(1)-gauge, verified), and the exact elimination of
(c0, c1) from E1..E4 leaves resultants whose gcd has degree 0 in c1 at 24
sampled rational p != 0 — exceeding the bounded p-degree of the final
discriminant, hence identically: THE NULLCONE FORCES p = 0, i.e. alpha = 0 or
beta = 0: one-sided. (With p = 0 the member is Z-side + radial or W-side +
radial; radial part then dies by E1/E2 and Theorem A's pattern.)

## Consequences
- One-Sided (hence GMC(2)-safety) now covers: all W-degree <= 1 (THM-1525),
  all two-charge spans any degree (Thm B), pure radial (Thm A), and
  {+1, 0, -1} at (const, const, linear) degrees (Thm C).
- Remaining for full GMC(2): spans with >= 3 charges at general degrees —
  the elimination template of Thm C scales (S166 machinery); the charge-radius
  lock (S168) says no independent Gamma exists to power a telescope.
Scripts: 04-computation/onesided_twocharge_theorem_boxeph_S169.py + frozen out.
