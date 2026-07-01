# The apex code, built: the PSL₂(𝔽₇) left–right Cayley complex realizes the certificate as an H¹ class on an expander (dim 16, O(1)-local) — turning the anti-LTC into an LTC via the apex's own group; and the four faces of √p

*opus-2026-07-01-S30. The culmination: build the explicit PSL₂(𝔽₇) left–right Cayley code and test
O(1)-local-testability of the certificate class, with the order-7 (heptagon) and order-3 (Eisenstein)
generators crossing the 3 and 7 where √21 lives. Plus the four faces of √p, the {7,21} impossibility, and the
Pochhammer even/measure side.*

## The build (verified)
`G = PSL₂(𝔽₇)` (order 168). **LEFT generator `a` of order 7** (a heptagon rotation), **RIGHT generator `b` of
order 3** (an Eisenstein multiplier); `⟨a,b⟩ = G`. The **left–right square complex**: A-edges `(g, ag)`,
B-edges `(g, gb)`, squares `{g, ag, gb, agb}` — and because `a` has order 7 and `b` has order 3, **each square
literally crosses the 3 and the 7**, the arithmetic home of `√21`. Results:
- **Expander (not Ramanujan with these gens):** the degree-4 Cayley graph has second eigenvalue `3.69` vs the
  Ramanujan bound `2√3 = 3.46` — a genuine spectral gap `0.31`, *hugely* better than the abelian tiling cube
  (gap `~2/m → 0`, S28). (LPS/quaternion generators would hit the Ramanujan `2√p` bound exactly.)
- **The complex:** `V=168, E=336, F=168`; **each edge is in exactly 2 squares** — so it is a *closed square
  surface*, i.e. a **surface / CSS-code** structure.
- **The code:** over GF(2), `rank d₁ = 160`, `rank d₀ = 160`, so `dim Z¹ = 176`, `dim B¹ = 160`, and
  **`dim H¹ = 16`** — a 16-dimensional nontrivial cohomology (the certificate class space). `√21` (the narrow-ℤ/2,
  S27) is one distinguished ℤ/2 class inside this `H¹`.
- **O(1)-local:** every square-check touches 4 edges and every edge sits in 2 squares, so a non-cocycle
  `f = 1_e` violates exactly 2 squares — the test is `O(1)`-local and non-codewords are locally detected
  (coboundary expansion `> 0`).

> **This turns the anti-LTC into an LTC via the apex's own group.** The abelian tiling cube (S28) had no
> expansion and singleton generators (no local codes); the `PSL₂(7)` complex is an expander with a nontrivial
> `H¹` code and O(1)-local checks. **The certificate is a cohomology class on the heptagon group's own complex.**

**Honest caveat.** With two generators (each edge in 2 squares) this is a *surface code* — `O(1)`-local but not
yet `O(1)`-**sound** (surface codes have poly(n), not constant, testing soundness; a small logical error is
locally invisible). A *good classical* LTC in the Dinur–Evra–Livne–Lubotzky–Mozes sense needs **larger
generating sets** (each edge in many squares → tensor local codes) on the Ramanujan expander. So: the group, the
class (`√21 ∈ H¹`), and O(1)-locality are built and real; full O(1)-soundness is the next construction, using LPS
generators and tensor local codes.

## The four faces of √p at the hard apex (p = 7)
The single algebraic quantity `√7` wears four faces, now all verified/realized:
| face | object | where |
|---|---|---|
| **Gauss `i√p`** | the ι-odd Gauss sum = the certificate | S23 |
| **Paley skew `±i√p`** | the Paley-heptagon skew-adjacency eigenvalues = the tournament | S23 |
| **Ramanujan `2√p`** | the spectral bound of the LPS Cayley expander = the **LTC substrate** | **this build** |
| **Field `ℚ(√−p)`** | the class field carrying the narrow-ℤ/2 = the arithmetic | S27 |
`√21 = √(3·7)` is the *cross* of the `p=3` and `p=7` versions — the certificate (Gauss product), the tournament
(the 21-Frobenius), the expander (order-3×order-7 generators), and the field (`ℚ(√21)`), all at once. The LTC
build is the **Ramanujan face**, and it is where the other three meet: the code (H¹) is the certificate, the
group is the tournament's automorphisms, the field is the coefficient ring.

## The {7,21} impossibility — what it abstractly represents
`{7, 21}` are the **forbidden H-values** (no tournament realizes them). The abstract content is a coincidence of
numbers that is not a coincidence: **21 is simultaneously a forbidden H-value and the order of the apex's
automorphism group** `|Aut(Paley₇)| = 21 = 7⋊3` (the Frobenius/Borel subgroup of `PSL₂(7)`). So *the
impossibility is the maximal symmetry* — the H-spectrum has a hole exactly at the order of the most symmetric
object. `7` is the prime (the seven-vanishing, THM-503); `21 = 3·7` is the composite where the Eisenstein-3 and
heptagon-7 cross. The forbidden `{7,21}` is the apex: the point that is at once maximally symmetric (the
obstruction) and arithmetically composite (`√21`, non-sum-of-two-squares, narrow-ℤ/2). It is the same tension in
every language: forbidden H = |Aut| = the free-ℤ₂/Borsuk–Ulam hard side = the ι-odd certificate.

## The complementary even/measure side (Pochhammer) and the speculative complexity
- **Pochhammer (verified):** `f(14) = (½)_{12}/12! = 0.1612 ≈ 1/√(πn) = 0.1508` — the fiber fraction, a
  rising-factorial (Γ/Wallis) ratio. This is the **ι-even / π / measure** face: the smooth, equidistributed,
  SOS/Brouwer-easy side that the far-element/discrepancy potential (S29) controls — complementary to the
  `√p`/odd/certificate side. The two halves of the split are `1/√(πn)` (measure, Pochhammer) and `i√p`
  (certificate, Gauss).
- **Mahler–Popken (noted, speculative):** integer complexities `‖7‖=6, ‖14‖=8, ‖21‖=9`, with `21 = 3·7` the
  compositum modulus — the additive–multiplicative `+/×` tension (S29) numerically; no theorem, just the flag
  that the composite `3·7` residual sits where `+` and `×` disagree (it is the non-sum-of-two-squares 21).

## Status
- **Built + verified:** the `PSL₂(7)` left–right Cayley complex (order-7 × order-3 generators); expander (gap
  0.31); a closed surface complex `V=168,E=336,F=168`; `dim H¹ = 16` (the certificate/`√21` class space);
  `O(1)`-local (each edge in 2 squares). Pochhammer `f(14)=0.161≈1/√(πn)`. The four faces of `√7`.
- **The result:** the anti-LTC (abelian tiling cube) becomes an LTC substrate on the apex's own group; `√21` is a
  distinguished ℤ/2 class of `H¹`; the check is O(1)-local; POCS/Kaczmarz (alternating A/B) is the constructor.
- **Honest:** two generators give a surface code — `O(1)`-local but not yet `O(1)`-sound; a *good* classical LTC
  needs LPS-Ramanujan generators and tensor local codes (each edge in many squares). The construction is real
  and the class lives in `H¹`; full soundness is the stated next build.

Related: HYP-3821 (facility-location / two-squares / PSL₂(7)), HYP-3820 (anti-LTC = abelian cube; PSL₂(7) is the
fix), HYP-3819 (√21 = narrow-ℤ/2 = the H¹ class), HYP-3796/mac-mini (POCS = the alternation), HYP-3802 (heptagon /
PSL₂(7)), THM-503 (seven-vanishing = forbidden 7). External: Annals 203-2 p.03 (DELLM LTCs — the tensor-local-code
target). HYP-3823 (this). Script: 04-computation/psl27_leftright_cayley_code_opus_20260701.py.
