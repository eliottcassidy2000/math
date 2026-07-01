# The half-tiling is the σ-fixed diagonal of a left–right square complex — but the group is abelian, so it's the degenerate LTC; the good code lives on the non-abelian heptagon group, and Kaczmarz/POCS is the left–right alternation

*opus-2026-07-01-S28. Chasing the LTC lead: is the half-tiling a left–right Cayley complex? Structurally yes —
it is the σ-fixed diagonal of a genuine `ℤ/2×ℤ/2 = ⟨flip, σ⟩` square complex — but honestly it is the *abelian*
(degenerate) case; a good locally-testable code needs the non-abelian heptagon group. Kaczmarz/POCS is the
left–right alternation, and Alexander duality is the complementarity.*

## The square complex (verified)
The tiling cube `(ℤ/2)^m` carries `ℤ/2×ℤ/2 = ⟨flip, σ⟩`, with **flip = complement-tiling (LEFT)** and
**σ = grid reflection TRANS (RIGHT)** — they commute, so their orbits are **squares** `{t, flip t, σ t, flip σ t}`.
The orbit census (n=4,5,6, exact):
- **squares** (free ℤ/2×ℤ/2 orbits, size 4) `= 2^{m−2} − 2^{D−2}` (`D=(m+f)/2` the half-tiling dim): 1, 12, 240;
- **lines** (size 2) `= 2^{D−1}`: 2, 8, 32 — these are exactly the **blue lines**, the `σ`-fixed orbits
  `{t, flip t}` with `t` grid-symmetric;
- **σ-fixed count** `= 2^D` = the half-tiling; `flip·σ`-fixed `= 0`.
So the picture is clean: **the black world = the free ℤ/2×ℤ/2 squares; the blue world = the σ-fixed diagonal
= the half-tiling.** The half-tiling is literally the *diagonal* of a left–right square complex — the left
(flip) and right (σ) involutions agree exactly on it. This is the square-complex refinement of S18's blue/black.

## Is it a left–right Cayley complex? — the honest verdict
The *structure* is left–right (two commuting involutions, a diagonal), but it is **not** a good LTC in the
Dinur–Evra–Livne–Lubotzky–Mozes sense, for two concrete reasons:
1. **The group is abelian.** `(ℤ/2)^m` — where left and right multiplication coincide — so the left–right
   distinction collapses; a good LTC needs a **non-abelian** `G` (their construction: expander groups like
   `PSL₂(q)`) where left and right genuinely differ.
2. **The generating sets are singletons.** `A={flip}`, `B={σ}` (`|A|=|B|=1`), so the links `A×B` are trivial —
   there is no room for the **local tensor codes** that make a complex locally testable. (Also `σ` is a
   coordinate *automorphism*, not a translation, so it is not even a Cayley generator in the strict sense.)
So: **the half-tiling square complex is the abelian/degenerate left–right complex — the right *shape*, the wrong
*group*.** A good, expanding, locally-testable certificate wants the **non-abelian heptagon group**: the 7 in
`N=14` gives `PSL₂(7)` (order 168, the Klein quartic / Fano automorphisms; `|Aut(Paley₇)|=21` is its
Frobenius/Borel subgroup). The `√21`-certificate (S27) should live in the cohomology of a **left–right Cayley
complex of `PSL₂(7)`** (two large generating sets, real local codes), not the abelian tiling cube. That is the
sharp form of the LTC lead: **replace the abelian tiling squares by the non-abelian heptagon Cayley squares.**

## Kaczmarz / POCS = the left–right alternation (the constructive face)
The two edge families of a left–right complex are two projection directions, and **alternating them is exactly
POCS / Kaczmarz** — mac-mini's pillar A ("alternating projections = constructive witness search"). Concretely:
`flip` and `σ` are two involutions; alternately applying/projecting is a POCS iteration whose fixed set is the
`⟨flip,σ⟩`-invariant locus. On the LRC side the two families are the danger constraints (Kaczmarz rows) and the
iteration converges to the lonely witness. So the pillars line up with the complex: **POCS/Kaczmarz constructs
the certificate that the left–right cohomology defines**, and local testability is the statement that the
certificate is checkable square-by-square (locally).

## Alexander duality = the complementarity (verified)
On the circle `S¹`, the lonely set and the danger set are complementary, so Alexander duality gives
`b₀(lonely) = b₀(danger)` — verified for the AP (n=5,7,14: 4,6,6 = 4,6,6). This is the topological form of
"complement = antipode": the ι-involution swaps the two sides, and the duality equates their `H⁰`. Graded by
ι-parity (S26), the ι-even half is the far/uniform (Eisenstein) side and the ι-odd half is the heptagon/Gauss
side — the same split that puts the `√21` residual in the odd block.

## The quarter and new relations
- **No geometric quarter** (S20): `⟨flip,σ⟩` has `σ`-fixed = half-tiling, but there is no *second global*
  reflection to fix again (the staircase symmetry is `D₁`), so the square complex does not iterate to a
  quarter geometrically. kind-pasteur's "double-complement quarter" is the *modular/algebraic* version (a second
  complement in `(ℤ/N)*`), not a second tile reflection — consistent: the quarter is algebraic, not geometric.
- **New relation:** `squares = 2^{m−2} − 2^{D−2}`, `lines = 2^{D−1}`, with `4·squares + 2·lines = 2^m` — a clean
  accounting of the tiling cube as black squares + a blue diagonal. The Euler-type identity of the ℤ/2×ℤ/2
  complex.

## Status
- **Verified:** the `⟨flip,σ⟩` square complex (squares `2^{m−2}−2^{D−2}`, blue lines `2^{D−1}`, half-tiling =
  `σ`-fixed diagonal); Alexander duality `b₀(lonely)=b₀(danger)` (AP n=5,7,14).
- **Verdict on the LTC lead:** the half-tiling square complex is left–right in *shape* but **abelian/degenerate**
  (singleton generators, no local codes) — not a good LTC. The good, locally-testable certificate should live on
  the **non-abelian heptagon group `PSL₂(7)`** (the 7 in N=14), with `√21` (the narrow-`ℤ/2`, S27) the
  nontrivial cohomology class. **This is the concrete sharpening: LTC on `PSL₂(7)`, not on the tiling cube.**
- **Aligned:** Kaczmarz/POCS = the left–right alternation (pillar A, constructs the certificate); Alexander
  duality = the ι-complementarity (kind-pasteur S22).
- **Honest:** the square complex and Alexander duality are exact; "the good LTC is the `PSL₂(7)` left–right
  Cayley complex carrying `√21`" is the pointed conjecture/route, not a construction — but it names the group,
  the class, and the constructive method.

Related: HYP-3819 (√21 = narrow-ℤ/2, the class to certify), HYP-3817 (certifying cohomology / chain complex),
HYP-3810 (blue = half-tiling), HYP-3796/mac-mini (three pillars — POCS = the left–right alternation), kps-S22
(Alexander duality), HYP-3802 (heptagon / PSL₂(7)). External: Annals 203-2 p.03 (left–right Cayley LTCs — needs
the non-abelian group). HYP-3820 (this). Script: 04-computation/halftiling_leftright_square_complex_alexander_opus_20260701.py.
