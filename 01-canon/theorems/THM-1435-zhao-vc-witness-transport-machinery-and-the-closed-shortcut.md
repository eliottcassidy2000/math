---
id: THM-1435
title: "THE ZHAO VC WITNESS — the transport machinery is BUILT AND CONTROL-VALIDATED, and the shortcut is CLOSED. (A) The de Bondt conjugation is proved here as a symbolic MATRIX IDENTITY, not cited: for N = [[A,0],[B,Aᵀ]], T N T⁻¹ is symmetric iff B is symmetric, with T = (1/√2)[[I,iI],[iI,I]]. (B) The cotangent lift G(x,y) = (H(x), JH(x)ᵀy) has exactly that shape, with B_ij = Σ_k (∂²H_k/∂x_i∂x_j) y_k SYMMETRIC because mixed partials commute — that one line is the whole mechanism. (C) Implemented end to end and validated with TWO positive controls (P homogeneous quartic, Hess P nilpotent, Δ^m(P^m) = 0 for m = 1,2,3) AND a negative control (non-nilpotent JH ⟹ Hess P not nilpotent and Δ^m(P^m) ≠ 0) — so the test has power. (D) CORRECTION TO death-star-S61: its parenthetical 'Yagzhev-normalize G = L⁻¹∘F = X + H (JH nilpotent)' is FALSE. Verified exactly: JH³ ≠ 0, and e₁, e₂, e₃ are each nonzero while e₁+e₂+e₃ = 0 identically — the latter is *all* that det(I+JH) = 1 says. Separate vanishing is forced only by HOMOGENEITY, which Alpöge's H does not have. So Bass–Connell–Wright is load-bearing, not cosmetic, and it is the entire remaining content of the transport. (E) The VC-WITNESS DIMENSION is defined and bracketed: 5 ≤ vcwd ≤ ~20."
status: >
  (A) PROVED -- a symbolic matrix identity in generic A and symmetric B, checked at
      n = 1,2,3 (not an example: the entries are free symbols).
  (B) PROVED -- one line (equality of mixed partials), and verified structurally by the
      implementation at every control.
  (C) IMPLEMENTED + VALIDATED.  Two positive controls and one negative control, all
      behaving as required.  This is machinery, not a result.
  (D) VERIFIED-EXACT, and it CORRECTS a fleet colleague's reflection.  det JF = -2 and the
      triple collision are independently re-verified here (third independent check in the
      repo).
  (E) DEFINITION + BRACKET.  The lower bound 5 is cited (VC proven for cubic-homogeneous
      maps in dim <= 4: Wright n=3, Hubbers n=4) and NOT re-derived here.  The upper bound
      ~20 is death-star-S61's count, which this file does not re-derive either.
  THE WITNESS IS NOT PRODUCED.  What is produced is everything downstream of the BCW step,
  plus a proof that the BCW step cannot be skipped.  Stated as a partial result on purpose.
source: klein-2026-07-20-S337 (owner: think in terms of VC-witness dimension and explicit witness to Zhao's vanishing conjecture)
attribution: >
  THE COUNTEREXAMPLE IS NOT OURS.  The Keller map F is Levent Alpoge's, with co-credit to
  Akhil Mathew, obtained with Claude Fable 5 and announced 19 July 2026 (THM-1300's
  attribution blocks, mac-mini-S127/S129).  Every corollary -- Dixmier, Zhao's vanishing
  conjecture, the image conjecture, Mathieu subspaces -- is a corollary of Alpoge's theorem
  and is being chased by the whole field.  NO PRIORITY CLAIM IS AVAILABLE OR MADE HERE.
depends_on:
  - THM-1300  # the counterexample + its attribution record
related:
  - THM-1370  # the weighted symmetry (1,-1,-2), externally confirmed
  - "death-star-2026-07-20-S61 reflection: the dimension count M ~ 20 and the program"
script: 04-computation/zhao_vc_witness_transport_klein_S337.py (+ .out)
---

# THM-1435 — the Zhao VC witness: machinery built, shortcut closed

## 0. The target, and what is at stake

**Zhao's Vanishing Conjecture (VC).** For `P` homogeneous with `Δ^m(P^m) = 0` for all `m ≥ 1`,
one has `Δ^m(P^{m+1}) = 0` for `m ≫ 0`. Zhao's Prop. 1.2 makes the hypothesis structural:

> `P` homogeneous is **Hessian nilpotent** ⟺ `Δ^m(P^m) = 0` for all `m ≥ 1`.

So a *witness* to VC's failure is an explicit homogeneous `P` with **nilpotent Hessian** such
that `x + ∇P` is **not injective** — the non-injectivity forces `Δ^m(P^{m+1}) ≠ 0` for
infinitely many `m` via Zhao's inversion formula. Crucially the witness is **a collision, not
a vanishing pattern**: two points, checkable in finite time. No such `P` is known in any
dimension.

## 1. Independent exact verification of the input (third in this repo)

`u = 1+xy`, `F = (u³z + y²u(4+3xy), y + 3xu²z + 3xy²(4+3xy), 2x − 3x²y − x³z)`:

- `det JF = −2`, constant. ✔
- `F(0,0,−¼) = F(1,−3/2,13/2) = F(−1,3/2,13/2) = (−¼, 0, 0)`. ✔
- component degrees `(7, 6, 4)`; `JF` is **not** symmetric, so the de Bondt doubling is forced.

## 2. The conjugation lemma, proved rather than cited

**Lemma.** Let `T = (1/√2)[[I, iI],[iI, I]]`. For `N = [[A, 0],[B, Aᵀ]]`,

```text
T N T⁻¹ = ½ [[ A + Aᵀ + iB ,  −iA + B + iAᵀ ],
             [ iA + B − iAᵀ ,  A + Aᵀ − iB  ]]
```

which is **symmetric if and only if `B` is symmetric**. Conjugation preserves nilpotency, and
`N` is nilpotent iff `A` is.

Verified as an identity in *free* symbols (generic `A`, generic symmetric `B`) at `n = 1,2,3`.
The `i` is why the construction lives over `ℂ` and not `ℝ`.

## 3. The cotangent lift — the whole mechanism in one line

For `H` homogeneous of degree `d` on `ℂⁿ` with `JH` nilpotent, define on `ℂ^{2n}`

```text
G(x, y) = ( H(x),  JH(x)ᵀ y ),        JG = [[ JH, 0 ], [ B, JHᵀ ]],
B_ij = ∂/∂x_j ( Σ_k ∂H_k/∂x_i · y_k ) = Σ_k (∂²H_k / ∂x_i ∂x_j) y_k.
```

> **`B` is symmetric because mixed partials commute.** That is the entire reason the de Bondt
> hypothesis is met, and it is why the *cotangent* lift is the right lift.

`G` is homogeneous of degree `d`; `JG` is block lower-triangular with nilpotent diagonal
blocks, hence nilpotent. By §2, `G' = T∘G∘T⁻¹` has symmetric Jacobian, so `G' = ∇P` with `P`
homogeneous of degree `d+1` (recovered by Euler: `P = (Σ_k w_k G'_k)/(d+1)`), and
`Hess P = JG'` is **nilpotent**. Therefore `Δ^m(P^m) = 0` for all `m` **by construction** —
half of VC holds for free.

**Collision transport.** If `x₁ ≠ x₂` with `x₁ + H(x₁) = x₂ + H(x₂)`, then for any `w`, setting
`y_a = (I + JH(x_a)ᵀ)⁻¹ w` (invertible, since `det(I + JH) = 1`) makes `(x₁,y₁) ≠ (x₂,y₂)` collide
under `Id + G`, hence — `T` being a linear bijection — under `Id + ∇P`. One collision in, one
collision out.

## 4. Implemented, with the negative control that gives the test power

| run | `JH` nilpotent? | `Hess P` nilpotent? | `Δ^m(P^m)`, `m = 1,2,3` |
|---|---|---|---|
| control A: `H = (0, x₀³)` | yes | **yes** | `0, 0, 0` |
| control B: `H = (x₁³, 0)` | yes | **yes** | `0, 0, 0` |
| **negative**: `H = (x₀³, x₁³)` | **no** | **no** | **nonzero** |

Both positive controls produce `P` homogeneous of degree 4 with `∇P = G'` verified and
`Hess P` nilpotent. The negative control fails both tests, which is what makes the positive
results informative rather than vacuous.

## 5. The shortcut is CLOSED — correcting death-star-S61

death-star-S61 writes: *"Yagzhev-normalize `G = L⁻¹∘F = X + H` (**JH nilpotent**; `L` the
antidiagonal linear part, `det JG = 1`, the triple collision transports)."* The `L` and
`det JG = 1` are right; **the nilpotency parenthetical is false.**

With `L = [[0,0,1],[0,1,0],[2,0,0]]` (`det L = −2`) and `H = L⁻¹F − X`:

```text
JH³ ≠ 0.
e₁ = tr JH  = x(2x²y³ + 12x²yz + 60xy² + 9xz + 48y)/2          ≠ 0
e₂          = −x(6x⁵y⁴z + 18x⁴y⁵ + … + 48y)/2                   ≠ 0
e₃ = det JH = 3x²(2x⁴y⁴z + 6x³y⁵ + … + 48y²)/2                  ≠ 0
e₁ + e₂ + e₃ = 0        ← and this is ALL that det(I + JH) = 1 says.
```

**Why.** `det(I + JH) = 1 + e₁ + e₂ + e₃`, so the Keller condition is the *single* equation
`Σ e_i = 0`. Separate vanishing — i.e. nilpotency — is forced only when `H` is **homogeneous**,
because then `e_i` is homogeneous of degree `i(d−1)` and the summands cannot cancel across
degrees. Alpöge's `H` is not homogeneous, the cancellation is available, and it is exactly what
occurs.

> **Consequence.** The cotangent lift requires a homogeneous `H` with nilpotent `JH`.
> Alpöge's normalised `H` has **neither**. So Bass–Connell–Wright degree-reduction *and*
> homogenisation is not a bookkeeping step that a dimension count can gloss — it is
> load-bearing, and it is the **entire remaining content** of the transport.

This does not contradict death-star's dimension estimate (`M ≈ 20`), which is about the *cost*
of BCW. It corrects the claim that the step before it is already done.

## 6. The VC-witness dimension

**Definition.** `vcwd` := the least `M` for which there exists `P ∈ ℂ[w₁,…,w_M]` homogeneous
with nilpotent Hessian such that `x + ∇P` is not injective — the least dimension carrying an
explicit witness to the failure of Zhao's VC.

```text
5  ≤  vcwd  ≤  ~20        (upper bound: death-star-S61's transport count, ≤ 34)
```

The lower bound is cited, not re-derived: VC is proven for cubic-homogeneous maps in dimension
`≤ 4` (Wright, `n=3`; Hubbers, `n=4`), so no witness exists below 5. Both bounds are soft — the
lower one is a theorem about a restricted class, the upper one an unexecuted construction.
**Narrowing this bracket is the concrete open problem this file leaves.**

*On the name.* Vapnik–Chervonenkis dimension is the largest configuration a class can shatter;
`vcwd` is the smallest configuration that certifies a failure. Both measure *how big a
configuration certification costs* — but the shared initials are a pun, not a theorem, and
nothing here transports between the two notions.

## 7. What this is, honestly

The witness is **not** produced. What is produced is (i) every step downstream of BCW, built
and control-validated, so that a cubic-homogeneous nilpotent counterexample can be converted to
an explicit `P` mechanically; (ii) a proof that BCW cannot be skipped; and (iii) a correction to
the one claim that suggested it could. Combined with THM-1300's attribution record: none of the
underlying mathematics is ours, and the corollary is being chased by the field right now.

*Files: `04-computation/zhao_vc_witness_transport_klein_S337.py` (+ `.out`).*
