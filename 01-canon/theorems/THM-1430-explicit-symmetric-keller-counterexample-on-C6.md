---
id: THM-1430
title: "AN EXPLICIT SYMMETRIC-CASE KELLER COUNTEREXAMPLE ON ℂ⁶, and the answer to whether the classical reductions are witness-effective. (0) THEY ARE. Both Bass–Connell–Wright and de Bondt–van den Essen are compositions with explicit polynomial AUTOMORPHISMS plus stabilisation, F′ = E₁ ∘ F^[m] ∘ E₂, so collisions transport in both directions with explicit sections — they are NOT truth-only equivalences. (Note BCW's Theorem 2.1(a) is *phrased* one-directionally; cite the construction, not the sentence.) (I) THE OBJECT: applying the dBvdE symmetric reduction DIRECTLY to the THM-1300 counterexample — skipping BCW, which the symmetric step does not require — gives G = T⁻¹∘Φ∘T on ℂ⁶ with Φ(x,y) = (F(x), JF(x)ᵀy) and T(u,v) = (u+iv, (u−iv)/2). VERIFIED INDEPENDENTLY: JG is SYMMETRIC (0 asymmetric entries), so on simply-connected ℂ⁶ the map G−z is a GRADIENT; det JG ≡ 4 = (det JF)²; and the three collision points transport by z_a = (a/2, −ia/2) to three DISTINCT points sharing the image (−1/8,0,0,i/8,0,0). (II) det is a PROOF not a computation: JΦ = [[JF,0],[∗,JFᵀ]] is block lower-triangular, so det JΦ = (det JF)², and T is linear. No 6×6 symbolic determinant is needed — an earlier run died inside simplify() attempting exactly that. (III) IT EXPLAINS opus-S422's NEGATIVE: the naive potential ⟨y,H(x)⟩ genuinely is not Hessian-nilpotent; the T-twist is the entire content of the reduction, since Hess P ~ −S·N with S·N block-triangular. (IV) HONEST LIMIT: Hess P is NOT nilpotent here (P is inhomogeneous because BCW was skipped), so this is NOT a Zhao vanishing-conjecture witness — that needs the homogeneous quartic, i.e. the BCW chain first"
status: >
  (0) Established from the construction shapes, with the literature's own phrasing
  caveat flagged.  Corroborated by precedent: Hubbers ran the same chain on Pinchuk's
  map (dim 2 → 203 → 1999) with maximum fibre size 2 preserved throughout.
  (I) VERIFIED-EXACT AND INDEPENDENTLY.  The construction was proposed by a research
  agent; I did not take it on trust.  Every claim below was recomputed from the
  definitions in exact arithmetic over ℚ(i): symmetry of JG by polynomial identity
  (expand, no simplify), the determinant by block triangularity plus four random
  spot-checks, and the collision by direct substitution.
  (II) PROVED.
  (III) Explanatory, and consistent with the S422 output it accounts for.
  (IV) VERIFIED at a random rational point: Hess(P)^k ≠ 0 for k ≤ 7.
  NOT VERIFIED BY ME: the agent's dimension counts for the full BCW→dBvdE chain
  (N = 79 cubic homogeneous, then 158 with a quartic P of 1660 monomials).  Those are
  reported as the agent's computation and are upper bounds from a naive greedy split.
  This does not resolve JC₂, DC₁ or DC₂, and it is not a VC witness.
source: kind-pasteur-2026-07-20-S128c112 (owner: pin the BCW/dBvdE constructions from the texts and run the transport; think of running equivalences witness-effectively versus equivalences being truth only)
depends_on:
  - THM-1300    # the counterexample being transported
related: [THM-1375, THM-1325]
script: 04-computation/symmetric_transport_verify_kps_S128c112.py (+ .out)
---

# THM-1430 — the symmetric transport, executed and checked

## 0. Witness-effective, not truth-only

The owner's distinction is the right one to ask, and the answer is favourable. Both
classical reductions have the shape

> `F′ = E₁ ∘ F^[m] ∘ E₂`,  `E₁, E₂` explicit polynomial **automorphisms**, `F^[m]` a
> stabilisation.

Injectivity therefore transports **in both directions**, constructively, with explicit
sections. There is no proof by contradiction and no genericity anywhere in the
construction. One caveat worth carrying: **BCW's Theorem 2.1(a) is *phrased*
one-directionally** ("if `F̃` is invertible then so also is `F`"), which reads truth-only.
The construction is manifestly bidirectional. Cite the construction, not the sentence.

Precedent settles it empirically: Hubbers ran exactly this chain on Pinchuk's degree-25
map, `dim 2 → 203` (Yagzhev cubic homogeneous) `→ 1999` (Drużkowski cubic linear), with
maximum fibre size 2 preserved throughout — a literal two-point collision surviving into
dimension 203.

## I. The object

The symmetric step does **not** require homogeneity, so it can be applied to `F` directly.
With `Φ(x,y) = (F(x), JF(x)ᵀ y)` on `ℂ⁶` and the twist `T(u,v) = (u + iv, (u − iv)/2)`,

> **`G := T⁻¹ ∘ Φ ∘ T`**,  `T⁻¹(x,y) = ((x+2y)/2, −i(x−2y)/2)`.

Recomputed from the definitions, not taken on trust:

| check | result |
|---|---|
| `T⁻¹ ∘ T = id` | exact |
| **`JG` symmetric** | **0 asymmetric entries** — so `G − z` is a **gradient** on simply-connected `ℂ⁶` |
| `det JG` | **≡ 4 = (det JF)²**, four random spot-checks |
| collisions | `z_a = (a/2, −ia/2)`; three **distinct** points, common image `(−1/8, 0, 0, i/8, 0, 0)` |

```
z₁ = (0, 0, −1/8, 0, 0, i/8)
z₂ = (1/2, −3/4, 13/4, −i/2, 3i/4, −13i/4)
z₃ = (−1/2, 3/4, 13/4, i/2, −3i/4, −13i/4)
```

> **`G` is an explicit non-injective Keller map on `ℂ⁶` with symmetric Jacobian** — a
> counterexample in the symmetric / gradient category, written down.

Field note: the main-diagonal-symmetric form needs `√−1`, so `P` and the witness
coordinates live over `ℚ(i)`, not `ℚ`. For a field-independent version one uses de Bondt's
anti-diagonal variant.

## II. The determinant is a proof

`JΦ` is **block lower-triangular**:

> `JΦ = [[ JF , 0 ], [ ∗ , JFᵀ ]]` ⟹ `det JΦ = det JF · det JFᵀ = (det JF)²`,

and `T` is linear, so `det JG = det JΦ = 4`. **No 6×6 symbolic determinant is required**,
and none was taken: a first attempt at this script died inside `simplify()` on exactly that
object. The spot-checks confirm the identity rather than establish it.

## III. Why the earlier in-repo attempt failed

opus-S422 recorded that the naive potential `⟨y, H(x)⟩` is **not** Hessian-nilpotent, and
that negative is correct — but incomplete. The missing ingredient is the linear twist `T`:
after conjugation, `Hess P` is similar to `−S·N` with

> `S·N = [[ JH , 0 ], [ M , JHᵀ ]]` block lower-triangular,

so `Hess P` is nilpotent **iff** `JH` is. The `S`-twist is the entire content of the
reduction. S422's other negative — that the naive *substitution* gadget breaks Keller-ness
— is likewise correct about the wrong gadget: BCW subtracts a full **product** and adjoins
**two** variables, and then `det J` is preserved exactly.

**And the repo was one transpose-inverse away from this since klein-S323.** That session
built the symplectic cotangent lift `Φ(x,ξ) = (F, JF⁻ᵀξ)`, which preserves `ω`. The dBvdE
object is the *same doubling* with `JFᵀ` in place of `JF⁻ᵀ` — and it is the transpose, not
the inverse-transpose, that is a gradient map.

## IV. What this is not

`Hess P` is **not** nilpotent here — verified, `Hess(P)^k ≠ 0` for `k ≤ 7` at a random
rational point. That is expected: skipping BCW leaves `P` inhomogeneous, while Zhao's
vanishing conjecture needs a **homogeneous quartic** with nilpotent Hessian. So this is a
symmetric-case counterexample and **not** a VC witness.

The VC witness needs the BCW chain first. Reported from the agent's computation and **not
independently verified here**: `N = 79` (cubic homogeneous) then `N = 158`, with `P` a
quartic of 1660 monomials over `ℚ(i)` — and both are *upper bounds* from a naive greedy
factorisation, likely reducible by exploiting that `F` is written in powers of `u = 1+xy`.

One further point the agent makes that is worth recording because it is cheap and rigorous:
*certifying* that the transported `P` violates VC needs no large computation, since the
chain is composition-with-automorphisms and THM-1300 §3 already proves `F`'s formal inverse
is non-polynomial. Exhibiting a **specific** `m` with `Δᵐ P^{m+1} ≠ 0` is the expensive
part, and Zhao's bound (`m > (3/2)(3^{n−2} − 1)`, astronomical at `n = 158`) gives no help.

## Named next

- Run the VC certification on the `ℂ⁶` object first as a dry run: 198 monomials, and
  S422's `lap`/`lappow` machinery already exists.
- Then the BCW chain, exploiting `u = 1+xy` to cut the gadget count below the naive 18.
