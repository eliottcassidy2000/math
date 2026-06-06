---
source: claude-2026-06-06-S687
status: VERIFIED (W6 χ=3 = roots of z^7−z; Moser spindle χ=4 explicit; cos θ=5/6, e^{iθ} root of 3z²−5z+3, disc −11) + the field-tower conjecture
tags: [hadwiger-nelson, moser-spindle, FTA, roots-of-unity, cyclotomic, eisenstein, heegner, mahler-measure, niven, chromatic-plane, field-tower, z7-z]
---

# Hadwiger–Nelson on the FTA root-locus: the χ=3→χ=4 jump is "leave the cyclotomic field"

The seed: *a degree-`n` polynomial has `n+1` coefficients and `n` roots (with
multiplicity) in the complex plane; the constant term is the "+1"; the roots map
to `n` points of ℂ.* But **ℂ = ℝ² = the Hadwiger–Nelson plane**. So the
fundamental-theorem-of-algebra root-map `ℂ^{n+1} → ℂ^n` lands inside the very
plane HN colors. Following that seed all the way down gives a clean, new
field-theoretic grading of the plane's chromatic obstruction. (Companion to
HYP-2275, which ran the same seed toward LRC: the covering-depth PGF, `p_0` =
constant term, worry-set ⟺ `z=0` is a root. This is the HN face.)

## The χ=3 floor is the cyclotomic locus: `z^7 − z`

The smallest Eisenstein unit-distance gadget is the **unit hexagon + its center**:
the six 6th-roots of unity (adjacent chord `2 sin 30° = 1`) plus the origin (radius
`1` to all six). Its unit-distance graph is the **wheel `W₆`** (hub + 6-cycle),
`χ = 3` (the even cycle wants 2, the hub forces a 3rd) — the local reason the
triangular lattice is not 2-colorable.

Those 7 points are **exactly the roots of**
```
z^7 − z  =  z · (z^6 − 1).
```
This *is* the seed's "`n` roots + the constant": `z^6 − 1` (the hexagon) times `z`
(the center/origin). And the same polynomial has a second life: **over `𝔽₇`,
`z^7 − z = ∏_{a∈𝔽₇}(z − a)`** (Fermat's little theorem), so its roots are *all of
`𝔽₇`* — the 7 colors of the hexagonal 7-coloring. One polynomial, two readings:
over ℂ it is the χ=3 *gadget*, over `𝔽₇` it is the 7-coloring *palette*. And
`7 = 6 + 1 = Φ₃(2) = |PG(2,2)|`, tying into the repo's `Φ₃(2^k)`/projective-plane
account of the forbidden tournament values (HYP-1043/1104/1419). The "needle" the
repo flagged — *are these 7s the same?* — gets a crisp partial answer here: they
are the **two evaluations of `z^7 − z`**, one geometric (roots in ℂ), one
arithmetic (roots = `𝔽₇`).

## χ=4 needs a second number field: the Moser spindle in ℚ(√−3, √−11)

The Eisenstein lattice `ℤ[ζ₆]` (the triangular lattice) is 3-colorable, so **every
finite unit-distance graph drawn inside `ℤ[ζ₆]` has `χ ≤ 3`**. Staying on the
roots-of-unity lattice caps you at the floor. To force `χ ≥ 4` you must leave it.

I built the **Moser spindle** explicitly and verified `χ = 4` (7 vertices, 11 unit
edges, exact backtracking): two `ℤ[ζ₆]` unit-rhombi `{0, 1, ζ₆, 1+ζ₆}` sharing the
vertex `0`, the second rotated about `0` by the spindle angle `θ`, with a unit edge
joining the two far vertices `1+ζ₆` and `e^{iθ}(1+ζ₆)` (both at distance `√3` from
`0`, so `2√3 sin(θ/2) = 1`). The arithmetic of that rotation is the whole story:

- `cos θ = 1 − 2 sin²(θ/2) = 1 − 1/6 = **5/6**`. By **Niven's theorem**, a rational
  cosine forces a rational angle only for `cos ∈ {0, ±½, ±1}`; `5/6` is not among
  them, so `θ/π` is **irrational** — `e^{iθ}` is **not a root of unity**.
- `e^{iθ} = (5 + i√11)/6` is a root of `**3z² − 5z + 3**` (discriminant `−11`,
  Mahler measure `3 ≠ 1`, hence a *non-cyclotomic* minimal polynomial — leading
  coefficient `3`, the prime-3 once more, but now as a *scaling*, not a lattice).
- Therefore the rhombi live in `ℚ(√−3)` (Eisenstein) and the rotation in
  `ℚ(√−11)`; the spindle is realized in the **biquadratic field `ℚ(√−3, √−11)`**.
  Both `−3` and `−11` are **Heegner numbers** (and `11` is on the repo's Heegner
  tail, HYP-2226).

The verifier confirms `rhombus A ⊂ ℤ[ζ₆]` while the **rotated far vertex leaves the
lattice** — the concrete sense in which `χ ≥ 4` escapes the cyclotomic locus.

## The new statement (and a conjecture)

> **On the FTA root-locus, the plane's chromatic obstruction is graded by the
> number field of the algebraic points.** Mahler-measure-1 / cyclotomic points
> (roots of unity; the Eisenstein field `ℚ(√−3)`; the lattice `ℤ[ζ₆]`) force only
> `χ ≥ 3`. Forcing `χ ≥ 4` requires a **non-cyclotomic** rotation — here a point of
> `ℚ(√−11)` with Mahler measure 3 — i.e. a **second imaginary-quadratic field**.
> The "incommensurate rotation" the repo kept flagging (S699g's 33.56°) is exactly
> *this algebraic point is not a root of unity*, and its field is `ℚ(√−11)`.

[CONJECTURE] `χ(ℝ²)` is bounded below by **how many pairwise-incommensurate
imaginary-quadratic rotations** one unit-distance graph can force. de Grey's
`χ ≥ 5` graph — a Minkowski-sum of many rotated Moser/spindle copies — would then
be a **tower of several such fields**, and the chromatic number of the plane would
be a statement about the rank of a multiplicatively-generated rotation group over
`ℂ^×`. (Speculative; the honest verified content is the `ℚ(√−3,√−11)` realization
of the single spindle and the `χ≤3` ceiling on `ℤ[ζ₆]`.)

## The FTA/spectral bridge (why "roots of unity" is the right object)

The user's seed and the repo's spectral HN bound are the same thing: for a
circulant distance graph the **eigenvalues are the connection-set polynomial
evaluated at roots of unity** (the FTA locus). The hexagon `C₆(1)` has eigenvalues
`(z + 1/z)|_{6th roots} = {±2, ±1}`, Hoffman `χ ≥ 1 − 2/(−2) = 2` (tight,
bipartite). As `m → ∞` the `6m`-th roots fill the unit circle and the polynomial
becomes the **Bessel `J₀`** symbol — the plane's `λ_min ≈ −0.4028`, `χ(ℝ²) ≥ 3.48`
floor (S699g / HYP-2264). So "evaluate a polynomial on the FTA root-locus" *is* the
spectral chromatic bound, discretely and in the continuum.

## Where this sits in the web

- **HYP-2275** — the FTA duality, LRC face (`p_0` = constant term, worry-set = `z=0`
  root). This is its HN twin.
- **HYP-2264 / S699g** — the HN/UD/LRC distance-Cayley unification, the `J₀` floor,
  and the "incommensurate rotation = Diophantine resonance" flag — now given its
  arithmetic identity `ℚ(√−11)`.
- **HYP-1043/1104/1419** — `{7,21} = Φ₃(2),Φ₃(4) = |PG(2,2)|,|PG(2,4)|`; the line-5143
  caution ("cyclotomic carrier echoes, not scalar identities") is respected — `z^7−z`
  is a *shared polynomial*, not a claim that the 7s are one number.
- **Cl₂(π/3)** (S599u) — the Mahler-measure / log-sin thread; the spindle's rotation
  having Mahler measure 3 is the same "Mahler measure grades the geometry" idea.

## Next

1. **Field of the next gadget.** Find the Golomb graph / a `χ=4` graph with a
   *different* rotation field, and the smallest `χ=5` configuration's field set —
   test the "tower of imaginary-quadratic fields" conjecture concretely (which
   `√−d` appear? are they forced distinct?).
2. **Which `−d` are reachable?** `cos θ = 5/6 → √−11`. For a unit edge between two
   norm-`N` Eisenstein points the rotation has `cos θ = 1 − 1/(2N)`, giving
   `√−(4N²−... )`. Map the function `N ↦` the rotation field; do Heegner `−d`
   dominate the small gadgets?
3. **Spectral via roots-of-unity.** Optimize a connection-set polynomial's minimum
   over roots of unity (a finite Delsarte LP) to push the circulant Hoffman bound
   above 3.48 — the discrete approach to the plane's measurable `χ` from the FTA
   side.
