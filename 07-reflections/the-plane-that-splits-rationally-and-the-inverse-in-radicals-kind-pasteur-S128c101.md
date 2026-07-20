# The plane that splits rationally, the inverse in radicals, and the module that will not pole

**kind-pasteur-2026-07-20-S128c101** (HYP-8150, THM-1345) · owner: *"work your
next steps and consider possibly the below"* + two explicit fiber families.

## 1. What the owner's families are, in the frame the fleet built

Family 1 (the r-continuum) is the ℂ*-flow collision family — the v = 0 slice of
Family 2, which is the real gift: an explicit rational splitting of EVERY fiber
over the plane {g = 0}. Placed in THM-1310's coordinates, three facts snap
together:

- the owner's s satisfies s² = v² − 16u = **−L restricted to the plane**: the
  splitting variable IS the quadratic resolvent. Family 2 says: over {g=0} the
  S₃-cover trivializes to the double cover of its own resolvent, rationally.
- the excluded parabola v² = 16u is **Jelonek ∩ {g=0}** — the owner's
  non-degeneracy condition is exactly "stay off the asymptotic variety."
- the x = 0 branch is a **polynomial section** σ(a,b) = (0, b, a − 4b²),
  F∘σ = id — a non-invertible étale map with an honest polynomial section over
  a coordinate plane. (Mechanism: over {g=0} the depressed cubic factors
  x(L₀x² + 4); the monodromy drops S₃ → C₂ over the plane.)

All three branches verified as exact rational-function identities (via
u = (v²−t²)/16, s = t — after a first modular-reduction attempt false-flagged
the x₋ branch: `subs(s²→…)` misses negative powers; lesson filed).

## 2. The inverse in radicals — and the Abel–Ruffini slice of the realization program

The fiber cubic is depressed; S₃ is solvable; Cardano is global: THM-1345 §2
writes the **complete 3-valued radical inverse** of F (one square root — the
√(−L) resolvent — then one cube root; verified numerically to 1e−14 across
branches). The counterexample is *non-invertible but radical-invertible*: its
failure of injectivity is exactly as tame as a cubic equation.

This slices the realization program (T1549) along the solvable/unsolvable line:
every currently-known counterexample (F, conjugates, towers F^{∘m}) has
solvable monodromy, hence fibers in radicals. Klein's Smith rule allows A₅ at
degree 5: an icosahedral Keller map would be the first **Abel–Ruffini-obstructed**
example — non-injective and non-radical — the polynomial incarnation of quintic
unsolvability, exactly as F is the polynomial incarnation of trisection
impossibility (THM-1335) and as Pisano-60 is the A₅-side clock (HYP-8145). The
ladder of classical impossibilities is becoming the ladder of Keller maps.

## 3. The trace module refuses to pole

The fit method (18+ exact rational samples per monomial, overdetermined) keeps
returning polynomials: Tr(xy) = −3 (a constant — dual to the trace laws
Tr(x) = Tr(u) = 0), Tr(y²), Tr(xz), Tr(yz), Tr(xyz), Tr(x²z²) all polynomial —
only the pure x-powers pole (Tr(x²), Tr(x³) with denominator L), and the master
identity explains why (x³ reduces with 1/L on fibers). The trace-polynomial
module is far larger than "the coordinates": it looks like everything OUTSIDE
the pure-x cone. Tower transitivity (Tr_{F∘F} = Tr_F∘Tr_F) turns the F∘F
centroid question into four module-membership checks — three confirmed, the
y³z fit at degree 8 recorded in the results file. Open shaped question: is the
module exactly ℂ[x,y,z] minus the pure-x directions — i.e., generated over the
polynomial centroid by {y, z}-mixing — and WHY (a z-affinity/depression proof
would explain THM-1335(4) rather than just compute it).

## 4. Honest ledger

- The x₋ false-flag was mine (mechanical), caught by the hand-check of F₃ = 0
  on the branch before touching the owner's formula — verify the instrument
  before doubting the input.
- Radical inverse is numeric-verified (1e−14 typical; one 7.8e−9 near-wall
  target); the symbolic identity is Cardano + the subresultant, both exact
  objects — a Lean-ready composition.
- Trace fits are overdetermined exact fits, not symbolic derivations; the
  y³z/y⁴ slots and the F∘F verdict are recorded in the .out as they land.
