---
id: THM-4562
title: "Slice lemma for Keller maps (a locally nilpotent dual field splits C^n = X_j x A^1 and reduces injectivity to a Keller map X_j -> A^(n-1), i.e. JC_n with a slice = cancellation in dimension n-1 plus JC on X_j), and THM-1300's Jacobian counterexample is fully non-slice: no dual field is locally nilpotent and no component is a coordinate or stable coordinate (fibre point counts over F_q)"
status: "PROVED (elementary: slice theorem for locally nilpotent derivations; uniqueness of derivations on separable algebraic extensions; spreading out and point counting). The fibre counts are PROVED for every q with gcd(q, 6) = 1 (derived in general by the audit; brute force at q = 5, 7, 11, 13, 17 here and q = 25, 49, 121, 125 in the audit). Found by the session's reader of openai/math #047 (affine cancellation). INDEPENDENTLY AUDITED 2026-10-06 (audit A: PASS WITH CORRECTIONS, applied: the lemma's (iff) needed 'F_j algebraically independent over R'; MISTAKE-579)."
session: mac-mini-2026-10-06-oaimath2
source: 05-knowledge/results/oai2_openai_math_second_reading_20261006.md
related:
  - THM-1300 (the map: u = 1 + xy, F = (u^3 z + y^2 u (4 + 3xy), y + 3x u^2 z + 3x y^2 (4 + 3xy), 2x - 3x^2 y - x^3 z); Alpoge 2026)
  - THM-3605 / THM-3561 (the Russell cylinder: #047's base threefold R = {xy = z(z+1)} x A^1 is THM-3605's Y_1 x A^1 under B = 4z, C = x, Y = 16y)
  - HYP-9216 (D(X) = A_4?)
  - Fujita (1979), Miyanishi-Sugie (1980): cancellation for A^2
---

# THM-4562 — slices, cancellation, and the non-slice counterexample

## Slice lemma (PROVED)

**Setting.** `F = (F_1, …, F_n) : C^n → C^n` is a Keller map (`det JF ∈ C^×`). Its dual fields `V_j = Σ_k ((JF^T)^(−1))_(jk) ∂_k` are polynomial, commute, and satisfy `V_j F_i = δ_(ij)`. (For THM-1300 these are the images of the `∂_j` under the Dixmier transfer.)

**Lemma.**
1. `V_j` is locally nilpotent **iff** `C[x] = R[F_j]` with `F_j` algebraically independent over `R` (that is, `C[x] = R^[1]` in the variable `F_j`), for some subring `R` containing every `F_i` (`i ≠ j`).
   * (⇒) This is the slice theorem: `F_j` is a slice of `V_j`, so `C[x] = ker(V_j)[F_j]`. Moreover `ker V_j ≅ C[x]/(F_j) = O(X_j)`, with `X_j = F_j^(−1)(0)`.
   * (⇐) `∂/∂F_j` over `R` (well defined because `F_j` is transcendental over `R`) agrees with `V_j` on all the `F_i`. A derivation of `C(x)` is determined on the separable algebraic extension `C(x)/C(F)`, so `V_j = ∂/∂F_j`, which is locally nilpotent.
2. In that case `C^n ≅ X_j × A¹`, and `F` becomes `(p, t) ↦ (G(p), t)` with `G = (F_i)_(i≠j) : X_j → A^(n−1)` étale. So **`F` is injective iff `G` is.**

**Consequences.**
* `n = 3`. Surface cancellation (Fujita; Miyanishi–Sugie) gives `X_j ≅ A²`. So a JC₃ counterexample with a locally nilpotent dual field would be a JC₂ counterexample.
* `n = 5` (CONDITIONAL on openai/math #047). #047's exotic `X` (`X × A¹ ≅ A⁵`, `X ≇ A⁴`) admits no injective Keller map to `A⁴`.
  * An étale injective map would be an open immersion.
  * The units are `C^×`, and Hartogs then forces `X ≅ A⁴`.
  * So any Keller map `X → A⁴` would give a slice-type JC₅ counterexample.

## THM-1300 is fully non-slice (PROVED)

**Fibre counts** over `F_q`, `gcd(q, 6) = 1`:

| component | weight | zero fibre | fibre over `c ≠ 0` |
|---|---|---|---|
| `F_1` | −2 | `2q² − 2q + 1` | `q² − q + 1` |
| `F_2` | −1 | `q² − q + 1` | `q² + 1` |
| `F_3` | +1 | `2q² − q` | `q² − q` |

These are PROVED for every `q` with `gcd(q, 6) = 1`.
* `F_3 = x(2 − 3xy − x²z)`.
* `F_1` and `F_2` are linear in `z`, with coefficients `u³` and `3xu²`.
* The exceptional loci `u = 0` and `x = 0` give the corrections.

They were also checked by brute force at `q = 5, 7, 11, 13, 17` (this session and the reader) and at `q = 25, 49, 121, 125` (audit).

* **No `V_j` is locally nilpotent.** Otherwise all fibres of `F_j` would be isomorphic to `X_j` (lemma, part 2). Spreading the isomorphism over a finitely generated ring and reducing mod almost every `p` would make the zero-fibre and nonzero-fibre counts equal. They differ for every `q`.
* **No component is a coordinate or a stable coordinate.** The fibres of a (stable) coordinate are (stably) `A²`, so they would have `q²` points mod almost every `p`. None of the counts is `q²`.
  * By Dutta–Lahiri (*On residual and stable coordinates*, J. Pure Appl. Algebra 225 (2021) 106707, Thm 3.4), a residual coordinate over `k[p]` is a stable coordinate. So no component is residual over a coordinate line either.
* **Fibre geometry (PROVED; re-derived by the audit):** `F_1^(−1)(c) ≅ A² ∖ {xy = −1}` and `F_3^(−1)(c) ≅ C^× × C` for `c ≠ 0`.

So THM-1300's counterexample is the opposite extreme from #047's stable coordinate `H`: its non-injectivity cannot be pushed down a dimension by the slice lemma. That fits JC₂ being open.

## Same threefold (PROVED, a change of variables)

#047's base `R = C[x, y, z, u]/(xy − z(z + 1))` becomes `Y_1 × A¹` with `Y_1 = {CY = B(B + 4)}` under `B = 4z`, `C = x`, `Y = 16y`. This is THM-3605's Russell cylinder, `≅ Y_2 × A¹` over THM-3561's surface.

So the stabilized THM-3561 map is an étale, non-injective map `A³ → Spec R`, and #047 principalizes its ideal on `SL_2 × A¹` over the same threefold. Any deeper link is OPEN.
