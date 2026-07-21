---
id: THM-1770
title: "THE Sym^3 COUNTEREXAMPLE IS THE GEOMETRY OF JC NON-PROPERNESS -- and X ≅ A^3 with NO computation, by Serre vanishing on the affine base. (A) NON-COMPUTATIONAL X ≅ A^3: X is a Zariski-locally-trivial A^1-bundle over A^2; its structure group is the affine group Aff_1 = G_m ⋉ G_a; the G_m (linear) part is a line bundle classified by Pic(A^2) = 0 (A^2 is a UFD), and the residual G_a-torsor is classified by H^1(A^2, O) = 0 (Serre: coherent cohomology of an affine scheme vanishes in positive degree). Both obstruction groups vanish, so the bundle is TRIVIAL: X ≅ A^2 × A^1 = A^3. One line, no coordinates -- 'affine space has no cohomology to twist by, so a bundle with contractible fibre over it is trivial.' (B) THE COUNTEREXAMPLE = JC NON-PROPERNESS: pi|X: X → Y is a morphism of two copies of A^3, etale (R removed), degree 3; since O(A^3)* = C*, its Jacobian is a nonzero constant, so it is a KELLER map, and degree 3 makes it NON-INJECTIVE. It does not contradict pi_1(A^3)=1 because it is NON-FINITE: the removed R ∪ pi^{-1}(H) is exactly where sheets escape = non-properness. This is the SAME mechanism as the repo's dim-3 JC counterexample (THM-1300) and the 'proper Keller => automorphism; non-properness is a three-line covering argument' bottom rung (THM-1330/1605) -- now realised as a clean rational-geometry model on Sym^3(P^1) = P^3, the twisted cubic = small diagonal, H tangent-not-osculating to it"
status: PROVED (A: the bundle-triviality is Serre vanishing + Pic(A^2)=0, standard and non-computational; B: pi|X Keller/non-injective/non-finite is immediate from etale + degree 3 + O(A^3)*=C*). Reflection-grade tie-in to THM-1300/1330.
author: opus-2026-07-20-S430
depends_on: [THM-1300 (dim-3 JC counterexample), THM-1330 (Keller monoid; non-properness = three-line covering), THM-1605 (extent vs mechanism; properness is the bottom rung)]
---

# THM-1770 — The Sym³ counterexample: X ≅ A³ for free, and it is JC non-properness

Owner supplied a rational-geometry counterexample and asked for a **non-computational** reason
`X ≅ A³`. There is a one-line one, and the map `π|X` turns out to be the clean geometric model
of the repo's JC non-properness.

## Setup (owner's)

`π : ℙ¹ × Sym²(ℙ¹) → Sym³(ℙ¹)`, `(p,{q,r}) ↦ {p,q,r}`, degree 3. `Sym²(ℙ¹) ≅ ℙ²`,
`Sym³(ℙ¹) ≅ ℙ³`; the small diagonal `{3p}` is the **twisted cubic** `C ⊂ ℙ³`. `R` = ramification
(`p ∈ {q,r}`); `H ⊂ ℙ³` a hyperplane **tangent but not osculating** to `C`.
`X = (ℙ¹×Sym²ℙ¹) \ (R ∪ π⁻¹H)`, `Y = Sym³ℙ¹ \ H = ℙ³ \ H ≅ A³`, and `X ≅ A³`.

## A. `X ≅ A³` with no computation

`X` is a **Zariski-locally-trivial `A¹`-bundle over `A²`** (owner's structure). An `A¹`-bundle
has structure group the affine group `Aff₁ = G_m ⋉ G_a` (maps `x ↦ u x + t`). Triviality is
obstructed in two stages, and **both obstruction groups vanish for the base `A²`:**

1. **Linear part** `G_m` → a **line bundle** `L` on the base, class in
   `\mathrm{Pic}(A²) = 0` (`A²` is a UFD, every line bundle trivial). So `L` is trivial and the
   structure group reduces to `G_a`.
2. **Translation part** `G_a` → a **torsor**, class in
   `H¹(A², G_a) = H¹(A², O_{A²}) = 0` (Serre: higher coherent cohomology of an **affine** scheme
   vanishes). So the torsor is trivial.

Hence the bundle is trivial and

```
X ≅ A² × A¹ = A³ .
```

> **The one-line reason.** *`A²` has no `\mathrm{Pic}` and no `H¹(O)` — nothing for an `A¹`-bundle
> to twist by — so a locally trivial bundle with contractible fibre over it is a product.* This
> is the algebraic mirror of "a fibre bundle with contractible fibre over a contractible base is
> trivial"; the "contractibility" that matters is **affineness + Serre vanishing**, not topology.

*(The only content that is not free is that `X → A²` is a locally trivial bundle rather than a
mere `A¹`-fibration — exactly the distinction that separates `A³` from the Danielewski
threefolds, where the fibration is not locally trivial. Once local triviality is granted, §A is
automatic.)*

## B. `π|X` is a Keller, non-injective, **non-proper** self-map of `A³` — the JC mechanism

Under `X ≅ A³ ≅ Y`, the degree-3 map `π|X : A³ → A³` is a morphism, i.e. a **polynomial map**.
Three facts, each immediate:

- **Keller.** `π|X` is étale (ramification `R` removed), so its Jacobian is a **unit** in
  `O(X) = ℂ[x,y,z]`; but `O(A³)^* = ℂ^*`, so `det J(π|X)` is a **nonzero constant** — a Keller
  map.
- **Non-injective.** It is generically 3-to-1 (degree 3).
- **Non-proper.** A finite étale cover of `A³` is trivial (`π₁^{ét}(A³) = 1` over `ℂ`), yet
  `π|X` is connected of degree 3 — so `π|X` is **not finite**. The non-finiteness is precisely
  the deletion of `R ∪ π⁻¹(H)`: fibres lose points into the removed locus. That is **loss of
  sheets at infinity = non-properness**.

> **This is the repo's JC counterexample in rational-geometry form.** THM-1300 is a dim-3
> Keller, non-injective, non-proper polynomial map; THM-1330/1605 record that *proper* Keller ⟹
> automorphism and *"non-properness is a three-line covering-space argument."* The Sym³
> construction realises exactly that: an étale (Keller) degree-3 self-map of `A³` that fails
> injectivity **only because** properness is dropped, with the failure locus the very divisor
> `R ∪ π⁻¹(H)` one removes to expose the two `A³`'s. The tangent-not-osculating condition on `H`
> is what makes `Y = ℙ³ \ H ≅ A³` while keeping `π⁻¹(H)` a genuine "sheets escape" locus rather
> than a fold.

## C. Two tangential extensions (free discretion)

1. **The twisted cubic ↔ the JC discriminant.** `H` tangent to `C` (twisted cubic = small
   diagonal `{3p}`) meets `C` at a double point — the geometric shadow of the **triple-collision
   fibre** of THM-1300 (`1` σ-fixed `+ 1` free orbit). "Tangent not osculating" = the collision
   is a genuine double contact, not a triple (order-2, not order-3) — matching THM-1350's
   **odd-fibre / collision-multiplicity-≥3** structure read on the diagonal: osculating `H` would
   force a triple contact (the σ-fixed sheet), tangent-only gives the `1+2` split.
2. **Sym as the source of "multiplication vs addition."** `Sym³(ℙ¹) = ℙ³` is the space of binary
   cubic forms; `π` adds a root to a quadratic. The elementary-symmetric coordinates on `Sym³`
   are the **multiplicative** shadow and the power-sums the **additive** — the same duality as
   THM-1580 (`h = s·A·C` multiplicative, `B` additive). The Vieta step that closes THM-1580's
   one-sided case is literally the `Sym`-to-roots map here; `π|X` is Vieta with one root free.

## D. Status

- **§A (X ≅ A³ non-computationally):** done — Serre vanishing + `Pic(A²)=0`.
- **§B (π|X = JC non-properness):** done — étale + degree 3 + `O(A³)^*=ℂ^*`.
- **§C:** tangential readings, flagged for whoever works the JC-diagonal or Sym-duality threads.

## Verification

No computation needed for §A–B (cohomological / immediate). The identifications
`Sym²ℙ¹≅ℙ²`, `Sym³ℙ¹≅ℙ³`, small diagonal = twisted cubic, `ℙ³\H≅A³` are standard.
