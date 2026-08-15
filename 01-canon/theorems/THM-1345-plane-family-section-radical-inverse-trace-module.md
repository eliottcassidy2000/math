# THM-1345 — The plane family and its section, the radical inverse, and the trace module of the Jacobian counterexample

**Status:** VERIFIED-EXACT (rational-function identities + numeric radical-inverse
verification + overdetermined exact trace fits; slots marked pending).

> **NON-RADICAL SUCCESSOR (THM-3438 / MISTAKE-396).**  The statements below
> about the sporadic cubic and its iterates remain exact, but “every known
> counterexample is radical-invertible” is superseded.  The weighted degree-five
> Keller map has geometric monodromy `S_5`, so it is an explicit non-radical
> inverse.  An `A_5` realization remains interesting but is no longer the first
> Abel--Ruffini-obstructed Keller map.
**Author:** kind-pasteur-2026-07-20-S128c101 (HYP-8150). Input: the owner's two
explicit fiber families. Builds on THM-1310 (fiber geometry), THM-1335
(trisection modulus, trace-polynomiality), mac-mini THM-1340 (engine
trichotomy), klein-S325 (Smith rule).

Setting: F the JC₃ counterexample; targets (a,b,g); L the Jelonek/resolvent
quartic; N(x) = L·x³ + (4−3bg)x − 2g the fiber cubic (depressed).

## (1) The owner's families, verified and placed

**Family 1** (r-continuum): F(0,0,−1/(4r²)) = F(±r, ∓3/(2r), 13/(2r²)) =
(−1/(4r²), 0, 0) — verified identically in r. This is the ℂ*-flow collision
continuum (c97/opus-S414); in the present frame it is the v = 0 slice of:

**Family 2** (the whole plane {g = 0}): with u the target's first coordinate,
s² = v² − 16u, x± = ±2/s:

    F(0, v, u − 4v²) = F(x±, v/4 − 3/(2x±), 13/(2x±²) − 3v/(4x±)) = (u, v, 0).

Verified as EXACT rational-function identities via the parametrization
u = (v² − t²)/16, s = t (both branches + the x=0 branch; a first attempt with
modular s²-reduction produced a spurious False on the x₋ branch — negative
powers of s escape `subs`; the rational parametrization is the clean method).
Placement in the THM-1310 frame:
- **s IS the resolvent:** −L|_{g=0} = b² − 16a = v² − 16u = s². Family 2 is the
  explicit rational splitting of the fibers over the double cover
  s² = −L restricted to the plane.
- **The excluded parabola v² = 16u is exactly Jelonek ∩ {g=0}** (L₀ = 0).
- The 13/2 z-numerator is universal on the plane (z = (2−3s)/x² with the
  collision decode 2 − 3s* = 13/2).
- **σ(a, b) = (0, b, a − 4b²) is a POLYNOMIAL SECTION of F over {g = 0}**:
  F∘σ = (a, b, 0) identically. The non-invertible F admits an honest polynomial
  section over a coordinate plane — the x = 0 leaf. (Over {g=0} the fiber cubic
  factors x·(L₀x² + 4): monodromy degenerates S₃ → C₂ over the plane, which is
  why the rational branch exists.)

## (2) THE RADICAL INVERSE

The fiber cubic is depressed and S₃ is solvable, so Cardano applies globally:
for a target Q = (a,b,g) off the walls, with p = (4−3bg)/L, q = −2g/L,

    x_j = ω^j·∛(−q/2 + √(q²/4 + p³/27)) − p/(3·ω^j·∛(…)),   j = 0,1,2,
    y_j = s(x_j)/x_j  (the degree-1 subresultant),   z_j  linear (z-affinity),

is a **complete explicit 3-valued radical inverse of F**. Verified numerically:
6 random targets × 3 branches, worst error 7.8e−9 (typical 1e−14). The radical
tower is: one square root (the discriminant — the √(−L) resolvent class of
THM-1310) followed by one cube root — Cardano = the trisection of THM-1335 made
algebraic; over {g=0} it degenerates to Family 2's rational-plus-quadratic form.

**Consequence for the salvage catalog / realization program (corrected
dichotomy):** `F`, its conjugates, and the towers `F^{∘m}` are
**RADICAL-INVERTIBLE** because their iterated-`S_3` monodromy is solvable.
THM-3438's weighted degree-five map instead has full `S_5` monodromy and is
the first explicit **Abel--Ruffini-obstructed** Keller inverse in the repo.
“Radical-invertible or not” remains a genuine coarse invariant, while `A_5`
is now a separate alternating-monodromy realization problem rather than the
threshold witness.

## (3) The trace module (fit method, 18+ exact samples each, overdetermined)

Beyond THM-1335's coordinate traces, the following are POLYNOMIAL:

    Tr(xy)   = −3            (constant! ⟺ Tr(u) = 0 cross-check)
    Tr(y²)   = −18a + 9b²/4
    Tr(xz)   = −9b²g/4 − 3b/2
    Tr(yz)   = −81a²bg² − 243a²g/2 + 135ab²g/4 − 21ab + 33b⁴g/16 − 21b³/8
    Tr(xyz)  = 81a²g²/2 − 27abg/2 + 78a − 15b³g/8 − 21b²/4
    Tr(x²z²) = −81a²bg³ − 297a²g² + 54ab²g² + 69abg − 338a + 33b⁴g²/16 − 5b³g/4 + 89b²/4
    Tr(z)    = −81a²g²/2 + 27abg − 51a + 15b³g/8 − 3b²/4   (independent re-derivation of THM-1335(4) ✓)

and, completing the set (deg-8 fits, 34 samples each, with the g=0 anchors
COMPUTED INDEPENDENTLY from the owner's Family 2 both checking True):

    Tr(y³z) = 2187a⁴g³/4 + 729a³bg²/2 + 1782a³g − 1863a²b³g²/8 − 3807a²b²g/8
              + 432a²b + 945ab⁴g/16 − 111ab³/2 + 129b⁶g/64 − 69b⁵/32
    Tr(y⁴)  = −81a²bg + 162a² − 27ab² + 33b⁴/16

while pure x-powers pole: Tr(x²) = 2(3bg−4)/L, Tr(x³) = 6g/L. The mechanism of
the poles is the master identity itself (on fibers x³ = [(3bg−4)x + 2g]/L: pure
x-powers reduce with 1/L; and a structural constraint: trace poles can only sit
on the Jelonek set {L=0} — Family 2 exhibits finite fibers over {g=0} — so no
g-poles are possible, as observed). The trace-polynomial module is LARGE: all
coordinates and every mixed monomial tested; only the pure-x cone poles.

**Tower transitivity ⟹ THE F∘F CENTROID IS POLYNOMIAL** (Tr_{F∘F} = Tr_F∘Tr_F):
Tr_{F∘F}(x) = 0, Tr_{F∘F}(y) = Tr_F(3y/2) = 9b/4, Tr_{F∘F}(z) = Tr_F(ζ) — a
polynomial since {y², xyz, x²z², y³z} are all trace-polynomial. **The
polynomial-centroid law (THM-1335's sharper-JC candidate) survives its first
non-trivial test and propagates through the composition monoid** — with NO
degree-9 computation (backlog c99(b) resolved).

## Files
- `04-computation/jacobian_radical_inverse_kps_S128c101.py` + `.out`
  (+ the focused y³z/y⁴ deg-8 fit run in the same .out family)
