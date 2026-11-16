# The dimension-3 Keller map as a fold family after an inversion, and what the plane lacks

Session `mac-mini-2026-11-15-tapes`, 2026-11-15. Status of the conjecture **in this repository**: the Jacobian
conjecture is false in every dimension `n >= 3` (THM-1300, an owner-supplied map, verified exactly there and
again here); **the planar conjecture JC₂ is OPEN, and nothing below proves or refutes it.** Not independently
audited. All identities are checked by `scripts/jacobian_levels.py` (`out/jacobian_levels.txt`).

## 1. The map, re-verified

    F = ( u³z + y²u(4 + 3xy),   y + 3xu²z + 3xy²(4 + 3xy),   2x - 3x²y - x³z ),      u = 1 + xy.

Exact computation: `det JF = -2`, degrees (7, 6, 4), and
`F(-2, 1/3, -2) = F(-4, 1/3, 0) = F(6, -1/6, 1/4) = (0, 1/3, -24)`. So `F` is a Keller map of `ℂ³` that is not
injective. (The collision of THM-1300, at `(-1/4, 0, 0)`, is a different fibre.)

## 2. Anatomy: a zero cancels a pole

Write `l = 2 - 3xy - x²z`, so that `F₃ = x·l`. Then

    F = Ψ ∘ τ,     Ψ(s, w, c) = ( s² + sw - cs³,  4s + 2w - 3cs²,  c ),     τ(x, y, z) = ( (1 + xy)/x,  1/x,  x·l ).

* `Ψ` is a polynomial map with `det JΨ = 2w`. It is the root cover of a cubic family: if `(a, b, c) = Ψ(s, w, c)`
  then `c s³ - 2s² + b s - 2a = 0` and `2w = b - 4s + 3cs²` is the derivative of that cubic at its root `s`. So `Ψ`
  is 3-to-1 and its critical divisor `w = 0` is where the root is double.
* `τ` is birational with `det Jτ = -x = -1/w`; it is an isomorphism of `{x ≠ 0}` onto `{w ≠ 0}`, with inverse
  `(s, w, c) ↦ (1/w, s - w, 5w² - 3sw - cw³)`.
* Hence `det JF = 2w · (-1/w) = -2`: **the simple zero of the fold family along its critical divisor cancels the
  simple pole of the inversion.** The image of `τ` misses exactly the critical divisor of `Ψ`, and `F` extends
  polynomially over the inverted plane `x = 0`.

The fold family is one member of a general shape: for every polynomial `r(s, c)`,

    Ψ_r(s, w, c) = ( sw + r - s ∂r/∂s,   2w - 2 ∂r/∂s,   c )      has      det JΨ_r = 2w,

and THM-1300 is `r = (cs³ - 2s²)/2`. (The tangent-sweep maps of THM-4625 and the fibre cubic of THM-4621 are
the same objects in other coordinates; the factorisation and the reading "zero against pole" are what is added.)

## 3. Level surfaces

* For `c ≠ 0` the level surface `{F₃ = c}` is the plane with coordinates `(s, w)` minus the line `w = 0`: a
  cylinder `ℂ × ℂ*`. On it `(F₁, F₂)` is the planar polynomial map `Φ_c = (s² + sw - cs³, 4s + 2w - 3cs²)`, of
  Jacobian `2w`, with its critical line deleted. In the original coordinates `(x, y)` the Jacobian is `2x⁻³`.
* `{F₃ = 0}` is the disjoint union of the plane `{x = 0}`, on which `F = (z + 4y², y, 0)` is an automorphism onto
  `{c = 0}`, and the cylinder `{l = 0}`, on which `F = ((v+1)(v+2)/x², 2(2v+3)/x, 0)`, `v = xy`, is 2-to-1.

So every collision of `F` happens on a cylinder, and the only level-surface component that is a plane carries an
injective map. This is forced if JC₂ holds, because of the (two-line) equivalence:

> JC₂ ⟺ every Keller map of `ℂ³` is injective on every component of a level set of one of its coordinates that is
> isomorphic to `ℂ²`.

(`⟸`: pad a planar Keller map by `z`. `⟹`: the restriction to such a component is an étale polynomial self-map of
the plane.) After THM-1300 the planar conjecture is exactly the statement that **planes** cannot play the role
that cylinders play here.

## 4. What the plane lacks

The planar analogue of section 2 would be `G = ψ_r ∘ τ₂` with `ψ_r(s, w) = (sw + r - s r', 2w - 2r')`,
`r ∈ ℂ[s]`, and `τ₂(x, y) = (γ x³ y + φ(x), 1/x)`, `φ` a Laurent polynomial: this is the only shape of `τ₂` with
`det Jτ₂ = γx = γ/w`, so that `det JG = 2γ`.

**Proposition (PROVED, for this ansatz only).** `G` is never a polynomial map.

*Proof.* The second component is `2/x - 2r'(γx³y + φ(x))`. On `y = 0` the pole of `r'(φ(x))` must be simple, so
`deg r = 2` and `φ = 1/(ρx) +` polynomial, `ρ` the leading coefficient of `r'`. Then the first component is
`1/(2ρx²) +` polynomial. ∎

In dimension 3 the cubic is allowed because **its leading coefficient is itself a coordinate**: `c = x·l` vanishes
on the inverted plane, so `c s³ = l u³/x²` has the pole order of a quadratic term; and the same coordinate supplies
`∂c/∂z = -x³`, which turns the slice Jacobian `2x⁻³` into a constant. In the plane the coefficients of the fold
family are constants and the pole count fails at the first step. One fibre dimension is rigid; two leave room for
a coefficient that is a function.

## 5. The planar work this touches

* **THM-2045**: the factorized coordinate `x(a - b x^r q^s)` has no planar Jacobian mate. `F₃ = x·l` is the
  spatial member of that family (reducible zero fibre, a plane and a cylinder), and it has two mates.
* **THM-2699 – THM-2722** classify when projections of polynomial graphs are planar Keller maps. `Ψ` is such a
  projection (the graph of a cubic family projected along its root); in the plane its critical curve cannot be
  removed, in space it is sent to infinity by `τ`.
* **THM-4727** (torus reduction): the quotient by the weights (1, -1, -2) is `(v, t) ↦ (lA, l²B)` with Jacobian
  `2l²`; in the coordinates above the torus acts by `(s, w, c) ↦ (s/λ, w/λ, λc)`, and `l = cw`.
* Reductions that raise the dimension (to cubic homogeneous or symmetric form) now land in statements that are
  false in general, so they can no longer carry the planar case as they stand; JC₂ still implies the Dixmier
  conjecture for the first Weyl algebra, and the tools that remain are two-dimensional (Newton polygons, the
  degree lattice, the planar strata of THM-2102 – THM-2129).

## 6. Two questions

1. Is there a Keller map of `ℂ³` with a coordinate all of whose level surfaces are planes and which is not
   injective? By section 3 this is JC₂ itself when the level surfaces form a trivial family; for a non-trivial
   family of planes it is a genuinely spatial question.
2. Classify the pairs `(r, τ)` for which `Ψ_r ∘ τ` is polynomial, with `τ` an inversion along a hyperplane and
   `c` a coordinate vanishing on it (HYP-9832 (c)); the known corpus (THM-4625) should be exactly the solutions
   with `c` of the form `x·l`.

**Analogy with the tournament problem of this session (ANALOGY only).** There, one tape gives a count that is
affine in one binary word and a bank of bounded capacity; a second tape makes the leading coefficient a free
product `432AB` (README, section 4). Here, one fibre dimension forces a pole and a second one lets a coefficient
vary. In both, the jump is from a constant to a free second parameter.
