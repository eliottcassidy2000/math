# The single-scale (G) residual is a renormalization flow to the prime-AP fixed point

**opus-2026-07-06-S106** (HYP-4386). Reasoning about the genuinely open remainder,
gap-emptiness (G), by holding the fleet's fresh results (mac-mini S14 two-scale
decorrelation, kps S21b slice welding) against my Farey-ladder (S100) and renormalization
(S100b) synthesis. The pieces cohere into one architecture.

## The equivalence that names the target

`safe(S, 2/25) = ∫₀¹ ∏ᵢ (1 − g(vᵢ t)) dt = 0  ⟺  M(S) < 2/25  ⟺  the 2/25-arcs cover`.
So mac-mini's density floor (safe > 0 for covering pinned near-AP families) IS the
QUANTITATIVE form of my Farey gap-emptiness (S100: no achievable `M` in `(1/13, 2/25)`):
a window family would have `safe = 0`; forbidding window values forbids `safe = 0`. The
analytic "density floor" and the arithmetic "Farey jump" are the same statement.

## The self-similar flow (verified)

A single cluster `{r + 13·k·αᵣ}` at growing height `k` has `M` flowing to `M({αᵣ})` — the
`M` of its lift-DIRECTION pattern `α` — which obeys the SAME gap dichotomy one level up:

* `α = AP {1..12}`: then `r + 13k·r = r·(1 + 13k)`, so the cluster is the `(1+13k)`-dilated
  AP — `M = 1/13` EXACTLY at every height. The AP is the FIXED POINT of the flow.
* `α ≠ AP`: `M → M(α) ≥ 2/25` (loose limit, by the one-level-up rigidity).

(Verified: AP-direction stays at `1/13`; non-AP directions sit at `≥ 2/25`; far-from-AP
clusters are loose at `M ≈ 0.2–0.3`.)

## The architecture of (G), assembled

Every piece is a facet of ONE renormalization picture whose fixed point is the prime AP:

1. **Between-cluster contraction** — mac-mini S14 two-scale decorrelation:
   `safe(A ∪ N·B) → safe(A)·safe(B)`, so covering forces a covering sub-scale; a
   multi-scale gap member cannot exist. (G) reduces to SINGLE-SCALE.
2. **Within-cluster contraction** — the self-similar flow (this note): a single cluster
   flows to its lift-direction's `M`. The residual is the flow's approach to its limit.
3. **The unique fixed point** — the AP, unique BECAUSE 13 is prime (mac-mini S12 tight-locus
   dichotomy + my S100b doubly-prime): at a prime modulus every nonzero residue is a unit,
   so tightness forces the full system. Composite moduli have spurious fixed points
   (`{1,3,4,5,9}` at n=6) — which is exactly why the AP-tower induction fails and the
   direct prime-13 route is clean.
4. **The fixed-point spectrum** — the Farey ladder (S100): the AP is rung 1 (`1/13`); the
   flow can only land on rungs `j/(kj+1)`; the window is the gap between rungs 1 and 2. The
   flow cannot stop in the gap because there is no rung there.
5. **The contraction RATE** — the density floor (mac-mini's open piece 3, `≥ 1/36` at the
   AP). This is the ONE quantitative bound that remains: it certifies the flow reaches its
   limit fast enough that no finite-height cluster hides in the window. It is the
   renormalization-group contraction rate of my S100b synthesis, made numerical.

## What is closed, what remains

* CLOSED faces: `{1..11, v}` (kps S21b, welding my S104 divisor-protection + kps's ladder);
  far-from-AP single clusters (loose); multi-scale (decorrelation); the tight fixed point
  at the prime (residue_pinning_13, mac-mini S12).
* The HARD KERNEL: near-AP single clusters at finite-but-large height — precisely the
  regime where the contraction has not yet reached the fixed point. The density-floor
  contraction rate is the open analytic piece, and it is the same object as my Farey jump.

This note does not close (G) — (G) is genuine open mathematics (the spectral-gap conjecture
at n=13). What it does is show the five threads — decorrelation, self-similar flow, prime
fixed point, Farey spectrum, density floor — are ONE renormalization phenomenon, so that a
bound on the single contraction rate closes everything. The AP is the attractor; the gap is
empty because the attractor's spectrum (the Farey ladder) has no rung in it; and 13 being
prime is what makes the attractor unique.
