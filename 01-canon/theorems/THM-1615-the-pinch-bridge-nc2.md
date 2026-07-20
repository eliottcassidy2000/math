# THM-1615: THE PINCH BRIDGE — a rebuilt TNC => NC2 route that localizes instead of dominating

**Status:** REFUTED AS SKETCHED (S178 adversarial pass, MISTAKE-203). Referee
verdict: (d) is FATAL — crossing a fold curve of the REPRESENTATION is a
Stokes crossing, not a singularity of A: the r-contour deforms around the
migrating branch point; genuine singularities of the continued A occur ONLY at
TRAPPED PINCHES (branch-point collisions: a finite discriminant-of-discriminant
set) or ENDPOINT CONTACT r_j(t) = 0 — which the sweep itself protects
(|t*(0+)| = infinity => contact only as |t| -> infinity). A generic ray MISSES
the finite trapped set: the ray-genericity defense deletes exactly the
singularities the argument needed. (c) is a GAP: the THM-1565 Gevrey/rotation
citations do not transfer (integrand not entire in r; |t Lambda_r| < 1 fails
for r >~ |t|^{-2/w}; uniformity dead; repair needs an unproven |kappa_j| <= sigma
for every Puiseux branch). (b)'s sweep was verified only on a zero-radial
instance; with f_0(0) != 0, t*(r) -> 1/f_0(0) FINITE as r -> 0 — the stated
universal sweep is false, and ironically that is the only configuration where
the endpoint model could fire. SALVAGE CONDITION (the real remaining content):
prove ONE trapped, non-cancelling pinch (or an endpoint branch with finite
t*(0)) is reachable on the principal sheet — a Landau-set analysis. What
survives: lemma (a) (R-H critical values at every r) and the reduction; kept
below for the record.
**Author:** boxeph-2026-07-20-S177 (HYP-8455)
**Context:** CASE-gamma-bridge-domination-step + THM-1585: klein's domination
step is FALSE (top-term share -> 0). This bridge replaces global comparison
with a LOCAL singularity that no radial averaging can cancel.

## The route
NC2 <=> A(t) := integral_0^infty e^{-r} CT_u log(1 - t Lambda_r(u)) dr == 0
(klein's reduction [1], referee-confirmed), Lambda_r = the locked Laurent
polynomial P(sqrt(2r)e^{i theta}) in u = e^{i theta}.

(a) CRITICAL VALUES EXIST AT EVERY FIXED r (PROVED): for mixed Lambda_r the
cover H_r(u) = t has degree d >= 2 with cluster ramification only (M-1)+(N-1)
= d-2 of the required 2d-2 (Riemann-Hurwitz): at least d ramification points
lie in C^*, so finite nonzero critical values t*(r) exist — always, not
generically. [One-sided P: CT log == 0, no curves, the bridge correctly
does not fire — matches the true nullcone members.]

(b) THE SWEEP (verified numerically + scaling): |t*(r)| -> infinity as r -> 0
(the lock kills charged coefficients: weak coupling) and -> 0 algebraically
as r -> infinity (coefficients grow: |t*| ~ r^kappa, kappa < 0; measured
kappa ~ -0.64 on the test instance, monotone 3.33 -> 0.02). The fold curves
cross EVERY modulus.

(c) WATSON + ROTATION (THM-1565 pattern): the NC2 hypothesis makes A
asymptotically 0 at all orders with Gevrey bounds; r-contour rotation widens
the sector past the Gevrey threshold; Watson-Nevanlinna: A == 0 EXACTLY on a
wide sector.

(d) FIRST PINCH (the localization that replaces domination): continue t along
a generic ray in the sector; by (b) the deformed r-contour is eventually
pinched or the singularity hits the endpoint r = 0. Local models: fold —
integral of e^{-r} log(t - t_c - a(r-r_c)^2) contributes a sqrt(t-t_c) term
with coefficient ~ pi e^{-r_c}/sqrt(a) != 0; endpoint — (t-t_0)log(t-t_0)
with coefficient e^{-0} != 0. First-crossing simplicity holds for a dense set
of rays (genericity in the RAY, not in P — legitimate). A nonvanishing LOCAL
singular term contradicts A == 0. UNLIKE the dead bridge, no global sum is
compared: averaging cannot cancel a local branch singularity.

## Consistency checks
- Pure radial P: CT log = log(1 - t f_0(r)): endpoint model reproduces the
  Radial Lemma (THM-1565) — an independent second proof.
- One-sided P: F == 0 identically: no contradiction fires. Exactly the
  true nullcone.
- Sweep instance: Lambda_r = sqrt(2r)(u + 1/u) + 2r u^2: table in the frozen out.

## What review should attack
(i) uniformity of the Gevrey bounds under r-rotation with the log integrand;
(ii) the pinch-vs-dodge dichotomy for the deformed contour (can the r-contour
dodge forever within the damped region?); (iii) simultaneous pinches with
canceling phases on non-generic rays (handled by ray-genericity — verify);
(iv) the r -> 0 endpoint degeneration of the cover (M(r) cluster-size jumps).
