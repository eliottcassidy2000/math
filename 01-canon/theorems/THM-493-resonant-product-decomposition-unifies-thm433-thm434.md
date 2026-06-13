---
id: THM-493
name: The resonant-product decomposition — the Moser lattice IS the Minkowski
      product at a resonant angle; the non-product 3N-crossover is a product PLUS
      the THM-434 transverse bonus
status: PROVED (from THM-434) + VERIFIED exact-integer (every distance² over ℚ(√3,√11))
date: 2026-06-13
session: monad-explorer-2026-06-13
depends_on:
  - THM-434   # #units(L_t) = 12 + r_E(t): the THREE families of unit vectors of L_t
  - THM-433   # avgdeg additive under the GENERIC-angle product; N* is "non-product"
  - THM-431   # u(21)=57; N* ∈ [25,28]; u(28) ≥ 85 (Engel/Schade Moser lattice)
  - HYP-2460  # this theorem's hypothesis stub
external:
  - "Erdős, AMM 1946 (the Minkowski-sum / product unit-distance construction)."
  - "Engel et al., arXiv:2406.15317 (the Moser lattice; u(28) ≥ 85)."
---

# THM-493: the resonant-product decomposition

THM-433 studied the **generic-angle** Minkowski product (average degree is additive;
the 3N crossover `N*` is "non-product"). THM-434 counted the **unit vectors** of the
Moser-ladder lattice `L_t = ℤ[ζ₆] ⊕ ω_t·ℤ[ζ₆]`. This theorem shows they are **two
faces of one operation**: the Moser lattice IS the Minkowski product of the
triangular lattice with a copy of itself rotated by the **resonant (Moser) angle**,
and THM-434's *transverse* unit vectors are exactly the **extra edges** that the
resonant angle creates beyond the Cartesian product. The "non-product" crossover of
THM-433 is therefore a **product at the resonant angle**, and the resonance bonus is
what does the crossing.

## Setup

`ℤ[ζ₆]` = Eisenstein integers (the triangular lattice), `ζ₆ = (1+i√3)/2`.
Fix `t ≥ 1` with `d = 4t−1` **not** three times a square (so `ω_t ∉ ℚ(√−3)` and
`L_t` is a genuine rank-4 lattice), and let
```
   ω_t = ((2t−1) + i√d) / (2t),   |ω_t| = 1   (the t-th Moser angle, cos = (2t−1)/2t).
```
For finite `G, H ⊂ ℤ[ζ₆]`, the **Minkowski product at angle `ω_t`** is the planar
point set
```
   P = G ⊞_t H := { g + ω_t·h : g ∈ G, h ∈ H } ⊂ L_t.
```
Write `e(·)` for the number of unit edges (unit distances) of a point set, and for
`α ∈ ℤ[ζ₆]` let `m_α(G) = #{(g,g') ∈ G² : g − g' = α}` (ordered) and
`N(α) = α·ᾱ` (the Eisenstein norm `= x²+xy+y²` for `α = x+yζ₆`).

## Statement

**(0) `P` has exactly `|G|·|H|` distinct points** (the map `(g,h) ↦ g + ω_t h` is
injective).

**(1) [PROVED] Exact unit-distance count.**
```
   U(P)  =  e(G)·|H|  +  |G|·e(H)  +  Δ_t(G,H),

   Δ_t(G,H)  =  (1/2) · Σ_{α : N(α)=t}  m_α(G)·m_α(H)
             =  #{ unordered {(g,h),(g',h')} : g−g' = α, h−h' = −α, N(α) = t }.
```
The first two terms are the Cartesian-product edges (G-edges at fixed `h`; H-edges
at fixed `g`) — the **generic-angle count** of THM-433. The third term `Δ_t` is the
**resonance bonus**: diagonal edges that exist *only* at the resonant angle, indexed
by the `r_E(t)` transverse unit vectors of THM-434.

**(2) [PROVED] Generic vs resonant.** At a generic relative angle `θ` (so `e^{iθ}`
is transcendental over `ℚ(√−3)`) the product realizes the Cartesian product `G □ H`
with **no** diagonal edges: `U = e(G)|H| + |G|e(H)`, `Δ = 0` (THM-433). The bonus
`Δ_t > 0` requires (i) the resonant angle `ω_t`, AND (ii) both factors to contain a
matching `√t`-separation (`m_α(G), m_α(H) > 0` for some `N(α)=t`). For `t` **not**
Loeschian (`r_E(t)=0`, e.g. `t=2`: the `√7` rung) there are no transverse vectors,
so `Δ_t ≡ 0` for *every* `G,H` — resonance at that angle buys nothing.

## Proof

**(0) Injectivity.** If `g + ω_t h = g' + ω_t h'` then `g − g' = ω_t (h' − h)`. The
left side lies in `ℤ[ζ₆] ⊂ ℚ(√−3)`. If `h ≠ h'`, the right side is `ω_t` times a
nonzero Eisenstein integer; since `ω_t ∈ ℚ(√−d)` with `√−d ∉ ℚ(√−3)` (the rank-4
hypothesis), `ω_t·(nonzero) ∉ ℚ(√−3)` — contradiction. So `h = h'`, then `g = g'`. ∎

**(1) Decomposition.** Two product points `(g,h) ≠ (g',h')` are at distance 1 iff
the lattice vector
```
   v = (g − g') + ω_t·(h − h')  ∈  L_t        satisfies   |v| = 1,
```
i.e. `v` is a **unit vector of `L_t`**. THM-434 proves the unit vectors of `L_t` are
exactly three disjoint families (`z = a + ω_t b`, `a,b ∈ ℤ[ζ₆]`):

| family | shape | forces | edge type | count |
|---|---|---|---|---|
| triangular rosette | `ζ₆^j` (`b=0`) | `h−h'=0`, `g−g'=ζ₆^j` (unit) | G-edge at fixed `h` | `e(G)·|H|` |
| `ω_t` rosette | `ζ₆^j·ω_t` (`a=0`) | `g−g'=0`, `h−h'=ζ₆^j` (unit) | H-edge at fixed `g` | `|G|·e(H)` |
| **transverse** | `α(1−ω_t)`, `N(α)=t` | `g−g'=α`, `h−h'=−α` | diagonal | `Δ_t(G,H)` |

(The transverse identity uses `|1−ω_t|² = 1/t` exactly, THM-434 Step 0, so
`|α(1−ω_t)|² = N(α)/t = 1 ⟺ N(α)=t`.) The three `(g−g',\,h−h')`-signatures are
mutually exclusive, and by THM-434 they **exhaust** the unit vectors, so the unit
edges of `P` partition into the three types and `U(P)` is their sum. The transverse
count: each undirected diagonal edge `{(g,h),(g',h')}` has `g−g'=α, h−h'=−α` for a
unique `α` with `N(α)=t`; summing `m_α(G)·m_α(H)` over all such `α` counts each
undirected edge twice (once for `α`, once for `−α`), hence the `1/2`. ∎

**(2)** At a generic angle `L_t` is replaced by a dense (non-lattice) set; the only
*forced* unit distances are the Cartesian-product edges (Erdős 1946), and the
transverse family does not exist (`|（g−g')+e^{iθ}(h−h')| = 1` with both parts
nonzero would impose an algebraic relation on the transcendental `e^{iθ}`). So
`Δ = 0`. ∎

## Corollary — the u(28) ≥ 85 crossover is the resonant product `W₇ ⊞₃ R`

Take `G = W₇` (the Eisenstein rosette / 6-wheel: hub + 6 unit neighbours, 7 pts,
`e=12`), `H = R` (the unit rhombus `{0,1,ζ₆,1+ζ₆}` = two glued unit triangles, 4
pts, `e=5`), and the **Moser angle** `t=3` (`ω₃ = (5+i√11)/6`). Then `P = W₇ ⊞₃ R`
is a **28-point** unit-distance graph (Corollary 0: distinct) with
```
   U  =  e(W₇)·4  +  7·e(R)  +  Δ₃  =  48  +  35  +  2  =  85   >   84 = 3·28.
```
(Exact-integer recount over `ℚ(√3,√11)`, 0 coincident points:
`resonant_product_bonus_monad.py`.) **The same product graph has only `83 = P(28)`
unit distances at a generic angle** (the THM-433 product cap, `< 84`); the resonant
angle adds exactly the `Δ₃ = 2` transverse edges, and **those 2 edges are the entire
crossing of `3N`** (`83 < 84 < 85`). This makes THM-433's "the crossover is
non-product/irreducible" precise and constructive:

> The first construction to beat `3N` is the *same product* that realises the
> product cap, evaluated at the **resonant angle** — and the resonance bonus
> (THM-434's transverse vectors) is exactly what carries it across the line.

So `u(28) ≥ 85` (THM-431's cited ceiling for `N*`) is reproduced **in-repo as an
explicit resonant product** with an exact integer count, independent of Engel's
search. (It lands in the Moser lattice `L₃` because `W₇, R ⊂ ℤ[ζ₆]` and `ω₃` is the
Moser angle — consistent with "AMP/Engel extremals live in the Moser lattice".)

## The n = 27 frontier (consistent with `u(27)=81`, `N*=28`)

A curated exact search (`resonant_product_Nstar_search_monad.py`, all rungs
`t ∈ {3,4,7,9,12,13,16}`, factor patches = densest hex prefixes / parallelograms /
triangles) finds the best **two-factor resonant product** for `n=24..27`:
```
   n    24  25  26  27
   3n   72  75  78  81
   best 68  72  61  75      (all ≤ 3n)
```
None beats `3N` before `n=28` (where `4×7 → 85` does). The pure cube `K₃^□3 = 81`
ties at 27 with **zero** bonus — `K₃` is a *unit* triangle (only `N=1` displacements),
so no rung gives it a transverse partner. To beat 81 at 27 a resonant product would
need a factorization (`3·9`) whose 3-point factor carries a `√t`-pair *and* whose
9-point factor matches enough of them; the search shows the bonus never covers the
product deficit. This is **consistent with** (does not prove) the HYP-2299 conjecture
`u(27)=81, N*=28`. The clean two-tie pattern of THM-433 (products tie 3N at 27, 30)
survives the resonant family: `n=30` reaches `90 = 3·30` (generic 87 + bonus 3) — a
**tie**, not a beat.

## Relation to canon

- **Unifies** THM-433 and THM-434 with no contradiction (no court case): THM-433 is
  the `Δ=0` (generic) case; THM-434 supplies the `Δ`-indexing unit vectors; this
  theorem is the bridge. Sharpens HYP-2299's "non-product crossover" to
  "product-at-resonant-angle."
- **Re-derives** `u(28) ≥ 85` constructively (THM-431(B) had it cited/annealed).
- **Engineering:** a clean exact-arithmetic "resonant product" primitive — pick two
  triangular patches with matching `√t`-spectra, glue at `ω_t`, gain
  `(1/2)Σ_{N(α)=t} m_α(G)m_α(H)` extra unit distances. The bonus is a **correlation
  of the two factors' norm-`t` displacement spectra** — a designable quantity.
- **Resonance with the project's web:** this is another instance of "everything is
  the triangle" — the bonus lives on the `√t = √(4t−1)`-ladder of CM directions
  (THM-434), and the Moser lattice is `triangular ⊗ √−11`, the geometry-side mirror
  of the LRC two-tower `clock × shell` (THM-427) flagged in the S4 Moser reflection.
