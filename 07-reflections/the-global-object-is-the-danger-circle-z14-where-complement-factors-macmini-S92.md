# The global object that breaks arc-transitivity is the danger circle Z/14 — where complement factors into the Atkin-Lehner V₄

*mac-mini-2026-07-13-S92. Owner asked me to find the global object that breaks arc-transitivity — the
structure that could carry the 2·7 Atkin-Lehner V₄ that the Sₙ-symmetric tournament metagraph cannot
(S91). I found it, and it is not exotic: it is the **danger circle `Z/14 = Z/2 × Z/7`** itself. On it,
the complement involution **factors** `W₁₄ = W₂·W₇` — the exact factorization Sₙ forbids — because the
CRT decomposition canonically distinguishes the 2-part from the 7-part. This closes the S89→S92 arc:
the V₄ lives on the runner side (S91) precisely because the runner side carries the level-14 structure,
and imposing that level is what breaks the transitivity.*

---

## The obstruction, restated (S91)

Complement is the only involution on the tournament iso-class metagraph, and it is **irreducible**: a
factoring `W₁₄ = W₂·W₇` would split the `C(n,2)` arcs into two Sₙ-invariant blocks, but **Sₙ is
transitive on the arcs**, so the only invariant blocks are `∅` and everything. There is no canonical
"2-part" and "7-part" of the arcs. The `14 = 2·7` factorization is arithmetic, and the tournament's
symmetry washes it out. So the object that carries the V₄ must **break arc-transitivity** — impose a
structure under which some arcs are canonically distinguished.

## The object: the circle Z/14 with its CRT decomposition

The runners live on the circle at the tight time; the loneliness threshold is `1/14`, so the natural
circle is `Z/14`. And **`Z/14 = Z/2 × Z/7` by CRT**. That decomposition is exactly the missing
canonical 2-vs-7 split. On it the Atkin-Lehner V₄ is elementary and real (verified):

> `V₄ = { id, W₂ = (x ↦ x+7), W₇ = (x ↦ 7−x), W₁₄ = (x ↦ −x) }`, and
> **`W₂ ∘ W₇ = W₁₄` — complement factors.**

- `W₂` is the order-2 **translation** by `7` — the 2-part, the 2-adic descent `Z/14 → Z/7` (THM-580).
- `W₇` is the **reflection** `7−x` — the 7-part.
- `W₁₄ = −x` is complement/inversion — their product, with fixed points `{0,7}` (the 2-torsion).

The point is not that these particular maps are the arithmetic Atkin-Lehner operators (they are not —
see the honest bound). The point is that **on `Z/14` the complement is reducible**: the CRT gives two
commuting order-2 factors, so the V₄ *exists*. On the arc-transitive metagraph it does not. The circle's
arithmetic breaks the transitivity the tournament's symmetry imposes.

## The combinatorial shadows of the two factors

Both factors have a concrete home already in the corpus:

- **The 7-part = the apex-7 Paley tournament.** The QR/Paley tournament on `Z_7` (connection set the
  quadratic residues `{1,2,4}`) has scores all `3` (regular), is self-complementary (`7 ≡ 3 mod 4`),
  and has automorphism group the Frobenius `Z_7 ⋊ Z_3` of order `21` — which **breaks `S_7`**
  (order 5040). Its number of directed 3-cycles — its OCF odd-cycle count — is exactly
  > **`c₃(Paley_7) = C(7,3) − 7·C(3,2) = 35 − 21 = 14 = 2·7`.**
  The level `14` appears combinatorially as the odd-cycle count of the apex-7 tournament, and *both*
  primes sit inside `Z_7`'s arithmetic: the `7` is the additive `Z_7`, the `2` is the Legendre quotient
  `(Z_7)*/QR` (= QR-vs-nonQR = the complement itself).
- **The 2-part = the 2-adic descent** (THM-580, the `u = 2t` circle cover `Z/14 → Z/7`), which is
  exactly `W₂` on the circle.

So the tournament sees the 7-part as the Paley `c₃ = 14`, and the 2-part as the descent — and the
circle `Z/14` glues them by CRT into the V₄.

## The clean picture: tournament : runner :: X(1) : X₀(14)

This is a level-structure statement. The **bare tournament metagraph** (tournaments mod `Sₙ`) is the
analog of `X(1) = SL₂(ℤ)\ℍ` — bare objects, full symmetry, no arithmetic, no cusp form. The **runner
configuration on `Z/14`** adds a **level-14 structure** (the phases on the circle), the analog of
`X₀(14) = Γ₀(14)\ℍ`. The Atkin-Lehner V₄ acts on the level, i.e. on the `2·7` — and *imposing the
level is exactly what breaks the arc-transitivity*. This is why (S91) the V₄ generators W₂ (2-adic) and
W₇ (apex-7) are "on the runner side": the runner side *is* the leveled object. The bulk/cusp split, the
`f₁₄` obstruction, the whole modular thread — they belong to level 14, and the tournament was level 1
all along.

## The honest bound (what remains)

I found the **group**, not yet the **arithmetic**. On `Z/14` the V₄ acts by translations/reflections,
which are **fixed-point-free** (`W₂`, `W₇` free; `W₁₄` fixes the 2-torsion `{0,7}`). But the genuine
Atkin-Lehner involutions on `X₀(14)` have **class-number fixed points** — `W₂` fixes 4 CM points
(discriminants among `−4,−8,−7,−56`; klein-S59). The circle gives the V₄ *group* and the *level*; the
CM fixed-point counts require the **moduli lift** to `X₀(14)` (elliptic curves with the level-14
structure), which the flat circle cannot see. So:

> **Found:** the object that breaks arc-transitivity = the danger circle `Z/14 = Z/2 × Z/7`, on which
> complement factors into the Atkin-Lehner V₄. **Remaining:** lift the circle's V₄-group to the moduli
> `X₀(14)`'s fixed-point arithmetic — the functorial bridge, which is the genuinely modular content.

Note also (S91, resolved here): on the *tiling model* — the base-path staircase, another arc-transitivity
breaker (arcs split into cut = base path ⊕ cycle = tiles) — complement is **fixed-point-free** (the
antipode of the tile-hypercube), matching the AL regular action, unlike the iso-class complement that
fixes SC. The tiling model and the circle are the two natural breakers; the circle is the one that
carries the `2·7` *arithmetic*, the tiling model the one that carries the fixed-point-free *shape*.

## Coda — the arc closes

S89 asked how tournaments relate to the last bit; S90 found the odd-graph→cusp transfer provably blind;
S91 traced the blindness to one root — Sₙ-transitivity on arcs makes complement irreducible; and S92
finds the cure: the danger circle `Z/14`, whose CRT arithmetic supplies the `2·7` split the tournament
lacks. The V₄ was never going to be found by a second *tournament* involution — there is none. It is
found by *adding the level* — the circle `Z/14` — which is what the runner side is and the metagraph is
not. The tournament is the bulk at level 1; the cusp form is the arithmetic at level 14; and the circle
is the bridge object, group-complete and arithmetic-incomplete. The last inch of that bridge — the
moduli lift — is the same modular content the covering-min's Dedekind/`f₁₄` residual has been pointing
at all along.

---

*Cross-links: the obstruction — HYP-6565/S91 (complement irreducible under Sₙ), THM-584 (complement =
antipodal). The circle V₄ — verified `W₂·W₇=W₁₄` on `Z/14`. The 7-part — the apex-7 Paley tournament
(`c₃=14`, HYP-3733, klein-S10 apex cusp `d=7`). The 2-part — THM-580 (2-adic descent = `W₂`). The
moduli — X₀(14)=14a, klein-S59 (W₂ 4 CM points), HYP-3768/f₁₄ (the cusp-form residual). The level-
structure frame — the bare metagraph = X(1), the runner circle = X₀(14). My computation: HYP-6575,
`04-computation/global_object_arc_transitivity_macmini_S92.py`.*
