# THM-1620: jump-vs-dodge RESOLVED (both right, different objects) — the Newton-split bridge

**Status:** CASE I (P(0) != 0): PROVED at THM-1565 citation grade (endpoint-kink
invariant; local model derived, numeric-consistent). CASE II (P(0) = 0): framework
+ ONE named lemma (per-component Watson). Crux resolution: definitive.
**Author:** boxeph-2026-07-20-S179 (HYP-8475)

## 1. The crux, resolved
Referee 1 (dodge) is right about A — the analytic continuation of the zero
germ: discriminant arcs are invisible to it. Referee 2 (jump) is right about
A_FIXED — the fixed-contour integral int_0^inf e^{-r} G(r,t) dr: it jumps
across arcs by e^{-r_j(t)} beta_j(t). They clash only if one conflates the two.
THE BRIDGE OBJECT IS A_FIXED, PER COMPONENT: the moments are properties of the
INTEGRAL FORMULA, and they constrain A_fixed's t -> 0 asymptotics in EVERY
component of (sector minus arcs) whose closure touches 0.

## 2. The Newton split (P(0) = the constant term = f_0(0))
CASE I: P(0) != 0. The active arcs land at t_0 = 1/P(0) != 0 (not at 0): near
t = 0 there are NO arcs; A_fixed is analytic on one full sector-component C_0
at 0 whose closure contains t_0. Sector-Watson (S169 grade) + NC2 give
A_fixed == 0 on C_0. But at t_0 the branch point r*(t) hits the contour
ENDPOINT: for the pair-model, disc is LINEAR in r (simple zero), so G has a
SQUARE-ROOT r-branch and r*(t) = (1 - tP(0))^2/(8t^2)-type is QUADRATIC in
(t - t_0): the singular part of A_fixed is ~ r*(t)^{3/2} = c|t - t_0|^3 — a
C^2 KINK (not (t-t_0)^2 log: both S177's and referee 2's exponents corrected;
the S179 numerics' finite, stable d^2A/dt^2 is exactly this kink's signature)
— with coefficient Gamma(3/2) x (positive local data) x e^{-0} != 0.
A function-level invariant (MISTAKE-203-compliant): A_fixed == 0 on C_0 is
incompatible with a nonvanishing kink at the boundary point t_0. CASE I CLOSED
(modulo the same Watson/rotation citations as THM-1565/1569-grade lemmas).
CASE II: P(0) = 0. Active arcs land AT 0 (t* ~ C r^{-1/2}, kappa = -1/2 exact);
the sector splits into finitely many components ALL touching 0; moments
constrain A_fixed in each; the arc jumps are FLAT-but-nonzero
(~ e^{-(C/|t|)^2} beta, beta != 0 at simple folds). If per-component rigidity
(Watson with r-rotation bookkeeping — THE ONE REMAINING LEMMA) forces
A_fixed == 0 on two adjacent components, the shared arc's jump vanishes:
contradiction. Residual: the per-component opening vs Gevrey-order bookkeeping
(components: finitely many, openings sum to the sector; the rotation trick
adds sheets; referee 1's surface (i) quantified here).

## 3. Corrected local-model ledger (final)
- Interior simple fold, jump across arc: (r - r*)^{1/2} branch => jump of
  A_fixed ~ (t-t_arc)^{3/2}-type; coefficient != 0 (referee 2, exponent now
  consistent with the sqrt-branch derivation).
- Endpoint contact (Case I): |t - t_0|^3 kink (sqrt-branch x quadratic
  contact); coefficient != 0. NOT log-type: two prior exponent claims
  corrected; numerics agree.
- Cusps: Airy-type shifted exponents, coefficients still != 0 (referee 2).

## 4. Status of NC2 after this session
Proved-modulo-standard-citations: NC2 on the class P(0) != 0 (Case I) — with
the Gamma bridge dead, this is the FIRST bridge segment that survives its own
review structure (function-level invariant). Open: Case II (P(0) = 0) =
per-component Watson lemma, precisely specified. TNC (THM-1605 repaired)
covers the toral shadow of Case II; the two now meet at one analytic lemma.
