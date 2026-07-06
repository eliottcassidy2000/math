# The density floor is a Cohn–Elkies problem: opus's theta-sum IS Poisson summation over the relation lattice

*mac-mini-2026-07-06-S24 (HYP-4532). Owner: work even more creative refinements,
search for connections. The density floor (the sole open piece) reduces — per
opus-S114 — to a Beurling–Selberg majorant (the "Selberg tail bound"). This note
makes the connection that names the right machinery: opus's safe theta-sum
(HYP-4446) is *literally Poisson summation over the relation lattice*, which is the
**Cohn–Elkies / Viazovska** framework — the state of the art for exactly this shape
(Poisson + positivity + lattices). Every adjacent term is 0 prior repo mentions
(Poisson summation, Cohn–Elkies, Viazovska, lattice theta, quadrature, uncertainty).*

## opus-S114's clarity, and where it points

opus-S114 (GREEN `LRCLonelyOpen`): `safe(S,β) > 0 ⟺ M(S) > β`, so `safe(·,2/25)=0`
for *both* the AP (`M=1/13`) and any gap member (`M∈(1/13,2/25)`) — the theta /
Fekete / Paley / Freiman pictures are **faithful reformulations, not reductions**.
The genuine open case is bounded/single-cluster, needing **(i)** a height *upper*
bound (`max ≤ ?`; S113's Farey wall gives only the lower bracket `max ≥ (3k+2)/2`)
or **(ii)** the **Selberg–Beurling majorant tail bound**. Both n-specific.

This note is about route (ii) — and it identifies the framework where such bounds
are actually constructed.

## The connection: safe = Poisson summation over the relation lattice

The safe measure is an integral over a *line* in the torus:
`safe(S,β) = ∫₀¹ F(v₁t, …, vₙt) dt`, `F(x) = ∏ᵢ(1 − g_β(xᵢ))`. The pushforward of
Lebesgue under `t ↦ (v₁t,…,vₙt)` is a measure on `Tⁿ` supported on that line, and
**Poisson summation** for it reads

> `∫ F dμ_S = Σ_{a : Σ aᵢvᵢ = 0} F̂(a) = Σ_{a ∈ L(S)} ∏ᵢ f̂(aᵢ)`,

with `L(S) = {a ∈ ℤⁿ : Σ aᵢvᵢ = 0}` the **relation lattice** (rank `n−1`) and
`f̂(0)=1−2β`, `f̂(m)=−ĝ_β(m)`. **This is exactly opus's theta identity (HYP-4446)**
— now named: the density-floor theta-sum *is* Poisson summation of a fixed
test-function's transform over the relation lattice. The density floor —
`Σ_{a∈L(S)} ∏f̂(aᵢ) > 0` for non-AP `L(S)` — is a **positivity of a Poisson sum over
a lattice**, the precise shape of the **Cohn–Elkies linear-programming bound**.

## What the Cohn–Elkies / Viazovska framework buys

In Cohn–Elkies, one *chooses* the test function to certify a packing/covering
bound; the sharp bound comes from a "magic function" whose Poisson sum over the
extremal lattice telescopes — constructed by Viazovska from **modular forms** (the
dim-8/24 solutions). Transporting the dictionary:

- **the test function is the Beurling–Selberg majorant** `g⁺ ≥ g_β` (my S22): a
  *suboptimal* Cohn–Elkies function. Its band-limit `N ~ 2k²` (must carry
  opus-S113's width) is the "resolution" the certificate needs.
- **the sharp certificate is a magic function** — the width-`N` majorant whose
  Poisson sum is provably positive over every non-AP `L(S)` and zero over `L(AP)`.
  This is exactly route (ii), and the Cohn–Elkies/Viazovska machinery is where such
  functions are *built*, not just posited.
- **the extremal lattice is `L(AP)`** — the relation lattice of `{1,…,n−1}`,
  `{a : Σ aᵢ·i = 0}`. Its theta function (with the arc weights) hits `0`; the floor
  is the *isolation* of `L(AP)` as the unique such lattice among relation lattices
  of covering families. The natural home of that theta's modularity is the
  project's own **`X₀(14)`** (cusps `{1,2,7,14}`, the `14=2·7` apex) — the modular
  curve where a magic function for `n=13` would live.

So the floor is not merely "another reformulation": it is the reformulation into
the *one framework where Poisson-sum positivity problems are solved* — and it names
the missing object (a modular magic function of resolution `N~2k²` for `L(AP)`),
rather than leaving it a generic "tail bound."

## The dual reading: LP certificate ⇒ height bound (route i from route ii)

Cohn–Elkies is a *linear program*; its dual is a *configuration* (a lattice point).
A Cohn–Elkies certificate that `Σ_{L(S)} ∏f̂ > 0` fails only if `L(S)` admits a
short vector pattern matching `L(AP)`'s — i.e. `S` is a generalized AP with the
AP's short relations. But S113's Farey wall forces `max ≥ (3k+2)/2` and my lever
forces `q ≤ 2·max`, so a certificate-evading `S` is a *bounded-ratio generalized AP
of height ≥ (3k+2)/2*. The Cohn–Elkies dual therefore feeds route (i): the
test-function positivity, where it fails, pins `S` to a finite family of generalized
APs — which is exactly the **height upper bound** opus-S114 flagged as missing. The
two routes are LP-primal and LP-dual of the same certificate.

## Adjacent connections surfaced (all fresh)

- **Quadrature / Gauss nodes:** the AP at `t=1/n` *is* the `n`-point Gauss quadrature
  node set on the circle (roots of unity); `safe` is a quadrature-error functional,
  minimized (to 0) at the optimal nodes. Non-optimal nodes ⇒ positive error ⇒
  `safe>0`. The Beurling–Selberg majorant is the standard quadrature-error control.
- **Uncertainty principle (Bourgain–Clozel–Kahane / ±eigenfunctions):** the majorant
  `g⁺` with `ĝ⁺ ≥ 0` on a band and sign control is a `+1`-eigenfunction-type object;
  the `N~2k²` resolution is an uncertainty budget (space `β` vs frequency `N`).
- **Selberg's original use:** the Beurling–Selberg majorant was built for the
  *explicit formula* (bounding `Σ over zeros`) — the same Poisson-sum-positivity
  template. The LRC floor is that template with the relation lattice in place of
  the zeros.

## Net

- **The density floor = Cohn–Elkies positivity of a Poisson sum over the relation
  lattice `L(S)`** — opus's theta-sum, named as its true classical object.
- The **Beurling–Selberg majorant is a suboptimal Cohn–Elkies test function**; the
  sharp one is a **modular magic function** of resolution `N~2k²` for `L(AP)`, whose
  natural home is `X₀(14)`. This is the concrete construction target for route (ii).
- The **LP dual** turns the certificate into route (i): where positivity fails, `S`
  is pinned to a finite set of generalized APs — the missing **height upper bound**.
  The two open routes are one LP.

## Pointers

- opus HYP-4466 (safe⟺M>β, reformulation-not-reduction), HYP-4446 (theta-sum),
  HYP-4456 (structure×width, Farey wall); mac-mini HYP-4512/4522 (Beurling–Selberg
  N~2k², convergence), HYP-4432 (q≤2max), HYP-4482 (n=7); kps HYP-4467 (harmonic).
- New leads (0 prior mentions): Poisson summation / Cohn–Elkies LP / Viazovska magic
  functions / lattice theta / Gauss quadrature / Beurling–Selberg for the explicit
  formula. The relation lattice `L(AP)` and its theta on `X₀(14)` are the objects to build.
