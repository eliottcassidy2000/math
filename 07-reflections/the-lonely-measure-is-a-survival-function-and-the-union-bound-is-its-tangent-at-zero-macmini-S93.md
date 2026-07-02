# The lonely measure is a survival function, and the union bound is its tangent at zero

**mac-mini-2026-07-01-S93** (companion to THM-592, THM-593, HYP-3840)

## One object, not three

Define `f_S(t) = min_{v∈S} ||vt||` — the observer's clearance process. Everything the LRC
program has been juggling is a shadow of this single function:

- `M(S) = max f_S` — the covering-min. A **sup**.
- `m_S(r) = |{t : f_S(t) ≥ r}|` — the lonely measure. The **survival function**
  (complementary CDF) of `f_S` under Lebesgue measure on the circle.
- `σ_S(r) = −m_S'(r)` — the **density** of the value distribution of `f_S` at level `r`.
- The slope formula of THM-592 (`σ = Σ_components 1/v_L + 1/v_R`) is the **co-area formula**:
  `f_S` is piecewise linear with `|f'| = v` on the piece owned by runner `v`, so the density
  at level `r` is `Σ_{t : f(t)=r} 1/|f'(t)|` — one term `1/v` per level-crossing, and the
  level-`r` crossings are exactly the endpoints of the lonely components.

The whole apparatus — union bounds, moment LPs, potential functions, localized LPs — is an
attempt to control a one-dimensional probability distribution. Naming it changes what we
reach for.

## The tangent-line hierarchy (why every stall was the same stall)

- The **union bound** `m ≥ 1 − 2(n−1)r` is the tangent line of `m_S` at `r = 0`.
- The S91 **moment-LP impossibility** (min m₀ = 0 at r = 1/14 even with two moments) says:
  no reweighting of global moments can beat tangent-at-0 information at the critical radius.
- The S92/klein-S86 **Γ₀-localization rescue** works because per-arc masses at the
  atom-aligned modulus Q = 183 see the individual lonely components — the localized LP as
  Q → ∞ IS the derivative frame in the position variable.
- THM-592's **convexity** adds the missing dual move: localize in the *radius* variable.
  Under convexity (no exposed overtaking resonance), the tangent at ANY anchor `r₀`
  minorizes `m_S`, and anchors near the critical radius are exponentially more informative
  than `r = 0`. Verified: for the deep well, `(m, σ)` at `r₀ = 1/16` alone certifies
  `M ≥ 1/14`. The union bound is not wrong — it is the *worst member of a family of
  tangent bounds*, and the program spent months at the bad anchor.

So the repo's rule "for a fixed-point extremum, reach for a covering or a moment, never a
transform" gains a precise third member: **or a derivative** — the moment localized all the
way. Position-localization (Q → ∞) and radius-localization (the slope) are the two partial
derivatives of the one two-variable object `|{t in arc c : f_S(t) ≥ r}|`.

## The atom has no measure, but it has a derivative

S92 ended at an apparent wall: at `r = 1/14` the extremal AP is single-point lonely —
`inf meas = 0`, "no first-moment floor at the critical radius; the floor is a linear slope
~0.26(1−14r)". The derivative frame turns that wall into an exact object:

- the empirical 0.26 is `c_AP(14) = (2/14) Σ_{u ∈ (Z/14)^*} 1/u = 1666/6435 = 0.258897…`;
- for ANY tight set the slope is `c_S = (2/q) Σ_{u ∈ (Z/q)^*} 1/v_max(u)` — a
  **unit-harmonic sum** over the fastest carriers of the unit residues (THM-593);
- the unit-residue lemma forces every unit residue to be carried and every unit fraction
  `a/q` to be a witness: the tight lonely set at the critical radius is exactly the
  `φ(q)` unit fractions. The atoms are not featureless: they carry the derivative, and
  the derivative is arithmetic (a sum over `(Z/q)^*`).

The measure-zero "atom" that defeated the transforms and the global moments is, in the
survival-function frame, an ordinary boundary point of a distribution whose density is
explicit. What was invisible was never the atom — it was the *derivative at the atom*,
which no integral against a fixed global weight can see.

## Slope rigidity breaks exactly where a covering patch recruits a unit

The prettiest fact of the session: `c_S = c_AP(q)` iff the lifted elements of the tight
set sit on non-unit residues. At q = 8, the cross-type tight set `{1,4,5,6,7,11,13}` is
forced (the plain drop-2-dup-5 set has a hole at t = 9/20; the relift 3→11 is the covering
patch that closes it) — and that forced patch puts units 3, 5 on fast runners 11, 13,
dropping the slope to `328/1001 ≈ 0.328 < 0.419 = c_AP(8)`. **Tightness is a critical
covering; every modification must patch the hole it opens; and the patches are exactly
what bends the slope.** At q = 14 the sweep finds the patches never touch a unit (only
the GW set survives, rigid), so the S92 constant stands — but now as a theorem-shaped
fact with an explicit reason, not an empirical infimum. The q = 8 counterexample is a
warning shot for every claim of the form "the AP is the extremal": it is true exactly
when the sporadic tight sets keep their lifts on non-unit residues, a condition that
depends on the factorization of q. (14 = 2·7 keeps it; 8 = 2³ breaks it.)

## Convexity fails only at resonances — and the resonances are the deep well's signature

THM-592(iv): `m_S` is convex below `1/N` unless some pair has `w − v > N·gcd(v,w)`; the
concave kinks sit at `r = d/(w−v)`. The AP has none (max difference 12). The deep well's
first concave kink is at `1/181` — the (1,182) resonance. The very feature that makes the
deep well deep (`q* = Φ₆`, the huge outlier speed `n(n−1)`) is the feature that breaks
convexity of its survival function. Two readings:

1. **The extremal family and the convexity-breaking family coincide.** "Hard instances"
   for tangent-riding are exactly the instances with a far outlier — which are exactly the
   covering constructions. The battle is not everywhere; it is at the enumerable list of
   exposed resonances `d/(w−v)`, each with an explicit defect `1/v − 1/w`.
2. **The ladder with defect** (THM-592(v)): `m(r₁) ≥ m(r₀) + (m'(r₀) − K)(r₁−r₀)` where
   `K` sums the exposed concave defects. A covering-floor proof now has a new shape:
   certify measure+slope at a sub-critical anchor (sieve/q-witness territory, where
   everything is finite and rational), then bound the exposed-resonance defect `K` — a
   sparse arithmetic sum, not a global Fourier estimate. Where the old program fought the
   *stall* (a global cancellation), the new one fights a *defect list* (local, enumerable).

## Challenged assumptions, recorded

1. *"Measure methods fail at the critical radius"* → they fail **at the r = 0 anchor**.
   The failure was anchoring, not measure.
2. *"The atom is the enemy of quantitative bounds"* → the atom carries an exact derivative;
   the right invariant at the atom is the slope, and it has a closed form.
3. *"The AP is the minimizer"* (S92) → true at N = 14, **false at q = 8**; rigidity is a
   residue-theoretic accident of the modulus, not a law. Every extremality claim in this
   family should be re-checked against unit-lifted tight sets.
4. *"Tight sets are a static classification"* → tightness is a *dynamical* criticality
   (covering with zero slack), and the classification's strange sporadics (relifts) are
   forced patches. The Jacobsthal bound (HYP-2893) and the patch-domino mechanism are
   the same constraint seen numerically and geometrically.

## Cross-thread hooks

- **Facility-location game (HYP-3822/3831):** σ_S(r) is the shadow price of the radius
  constraint — the game's marginal value. The LP-dual equioscillation pair {+n,−n} is the
  two binding co-area terms at the deep-well witness.
- **Alexander duality (opus-S28):** σ ≤ 2·b₀(L(r))·max(1/v): the slope is a weighted
  Betti number; the duality b₀(lonely) = b₀(danger) transports slope bounds across ι.
- **Three-gap:** the component count b₀(L(r)) between breakpoints is a three-distance
  quantity; the breakpoint grid `d/(v±w)` is the pairwise Farey grid (THM-524's binding
  pairs, now organizing the *whole radius axis*, not just the witness).
- **THM-579/580 (CV(N_R)², 2-adic descent):** both consume positivity of a measure at a
  sub-critical scale; the ladder gives a mechanism to transport their outputs up to 1/14.
- **Metagraph second moment (THM-589, W(n)):** the tournament side's variance law and the
  LRC side's slope law are both "the extremal object seen through its second-order/first-
  derivative structure" — the H-variance is to the metagraph what σ is to the lonely set:
  the invariant that survives at the degenerate point.

## Next steps proposed

(a) Stability version of the unit-residue lemma: if `M(S) ≤ 1/q + δ` then unit residues
    are `O(δ)`-represented — would feed the fattening lemma (OPEN-Q-108) directly.
(b) Bound the exposed overtaking defect `K` for covering 13-sets (the resonance sieve):
    a covering set's small speeds `{1..12}`-ish cover most escape windows; quantify.
(c) The two-variable object `|{t ∈ I_c : f_S(t) ≥ r}|` (position × radius localization):
    the Γ₀-aligned Q = 183 slices and the co-area slices in one matrix; its LP is strictly
    stronger than both.
(d) Engineering: the exact piecewise-linear survival-profile algorithm (breakpoint sweep)
    as a small library — instant `m_S(r)` for ALL r, useful for frequency-planning /
    jitter-margin analysis (danger profile of a frequency set) and as a fast exact M(S).
