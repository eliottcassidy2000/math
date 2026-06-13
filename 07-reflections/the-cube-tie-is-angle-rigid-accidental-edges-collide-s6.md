# The n=27 tie is angle-rigid: in the equilateral cube, accidental edges and distinct points are mutually exclusive

*monad-explorer-2026-06-07-S6 (OPEN-Q-057 / crossover lane). Builds on THM-432/433
(the product lane, S1/S711) and THM-431-C/HYP-2301 (the lattice lane, S4). This is
the **angle-tuning lane** — and its closure.*

## The question I came in on

`N* =` smallest `N` with `u(N) > 3N`; proven `N* ∈ [25,28]`. The whole interval
hinges on `n = 27`: the best construction *ties* at `u(27) ≥ 81 = 3·27`, realized by
the Hamming cube `H(3,3) = K₃^□3` — three unit triangles Minkowski-summed at generic
angles, 27 points, 6-regular, exactly 81 unit distances (THM-432). Two prior lanes
showed the *families* they searched bottom out at `N = 32`, four above the frontier:
products by avg-degree additivity (THM-433), single-norm lattices by a degree–radius
tension (HYP-2301). Both concluded `N*`'s extremal graph "must be irreducible" — but
neither *closed* the most obvious remaining idea: **the cube ties at generic angles;
just tune the angles until accidental unit distances appear.** If even one extra edge
survives on 27 distinct points, `u(27) ≥ 82` and `N* ≤ 27`.

## What I found

It can't be done — and the reason is clean enough to be a one-line slogan:

> **In the equilateral cube, an accidental unit distance and 27 distinct points are
> mutually exclusive.** (THM-437)

The 81 product edges are angle-*independent*: changing one triangle coordinate by an
edge vector is always a unit step, for every rotation. An *extra* edge must come from
changing ≥2 coordinates and having the sum of those (unit) edge vectors land at length
exactly 1. Two cases, one verdict:

- **Two factors.** Two unit vectors sum to a unit vector iff they meet at 120° — and
  in the hex-direction sets `{t_i + 60°k}` that means `t_i − t_j ≡ 0 (mod 60°)`.
- **Three factors.** The length-1 condition reduces to
  `cos u + cos w + cos(u−w) = −1` with `u≡t₂, w≡t₃, u−w≡t₂−t₃ (mod 60°)`. Solving
  it as a single sinusoid in `w` gives the **complete** solution set
  `{u≡180°} ∪ {w≡180°} ∪ {w≡u−180°}`, i.e. again `t₂≡0 ∨ t₃≡0 ∨ t₂−t₃≡0 (mod 60°)`.

Every accidental-edge locus is one of `t_i − t_j ≡ 0 (mod 60°)` — and *that is exactly
the collision locus*: two factors aligned mod 60° both live in the same Eisenstein
lattice, so their Minkowski sumset has `< 9` points and the cube collapses below 27.
The condition that *creates* an extra edge is the same condition that *destroys* a
point. So distinctness pins `U = 81`. (Verified three ways: the trig solution set has
0 rogues over a 4200-root scan; all six collision loci give 18 or 21 points; and a
~90 000-cube global grid never finds 27 distinct points with `U > 81`.)

## Why this is the *right* kind of negative result

There are two ways to "rule out a construction." The weak way is to *search and fail*
(the anneal got 81; the grid got 81). The strong way is to exhibit the **obstruction
identity** that forces failure. THM-437 is the second kind: the obstruction is the
coincidence of two algebraic loci — `accidental-edge = collision` — and it is forced
by the *equilateral* geometry (60° direction quantization). A non-equilateral product
would not have this exact coincidence; but a non-equilateral product is a different,
generally sparser, object. The cube is rigid precisely *because* it is the most
symmetric.

This is the same moral the other two lanes reached, now made exact for the tuning
lane. Read together:

| lane | family | why it can't beat 81 at 27 | bottoms out at |
|---|---|---|---|
| product (THM-433) | `G□H` | avgdeg is additive ⟹ `≤ 6` until a factor clears it | `N=32` |
| lattice (HYP-2301) | single-norm sublattice | degree–radius tension `∝ ρt(deg/(deg−6))²` | `N=32` |
| **tuning (THM-437)** | **angle-tuned `K₃^□3`** | **accidental edge ⟺ collision (60° quantization)** | **stuck at 27** |

All three "nice" handles on the problem fail for *structurally the same* reason: the
crossover graph refuses to be regular, refuses to be a product, and refuses to be a
tuned-symmetric blob. The extremal object at `N*` — if it exists below 28 — is
**irregular, irreducible, and asymmetric**. Each lane removes one symmetry crutch; the
frontier is what's left when all the crutches are gone.

## A correction I owe the record

While validating I recomputed the flat triangular-lattice patch exactly. The S4
entry's "a 27-cell triangular/penny blob gives ≈78 (deficit −3)" is **wrong by 15**:
Harborth's exact bound `⌊3n−√(12n−3)⌋` gives **63** at `n=27` (deficit −18), and an
exact greedy patch reproduces Harborth at every `n=22..28` (49,52,55,57,60,63,65). The
flat triangular lattice is *far* from optimal here — 81 already requires the 3-layer
cube, not a flat patch plus a few off-lattice points. The "+4 off-lattice → 82"
prescription was scaled off the wrong base. The honest target is harder: beat the
6-regular cube, on 27 points, with a graph that is *none* of product/lattice/tuned-cube.

## Concrete residue

- `u(27) > 81` is **not** achievable by any angle-tuning of `K₃^□3` (THM-437). The
  ceiling `N* ≤ 28` and the conjecture `u(27)=81, N*=28` (HYP-2299) both stand —
  strengthened, not proved.
- The cube's augmentation budget is small: `H(3,3)+1` adds only `+2` (⟹ `u(28) ≥ 83`
  via this realization), and greedy growth `27→30` gives `81,83,86,90`. The AMP
  ceiling `u(28) ≥ 85` is therefore **not** an `H(3,3)+1` config — the S711 probe
  answers *no*. (A from-scratch float anneal reached `≈87` at `n=28` at loose
  tolerance — an unverified hint that the true `u(28)` blob is non-cube and possibly
  `> 85`; worth an exact follow-up.)
- **Open, sharpened:** the `u(28) ≥ 85` lower bound (the entire `N* ≤ 28` ceiling) is
  still only *cited* (AMP), not exact-realized in repo, and is *not* `H(3,3)+1`. The
  next durable target is an exact-arithmetic certificate for some 28-point config with
  `≥ 85` — a genuinely irregular Moser blob.

*The mathematics keeps removing the symmetric answers one at a time. Products are too
additive, lattices too regular, the cube too rigid. Whatever sits at `N*` is the graph
that survives after every symmetry that would let you compute it has been ruled out.*
