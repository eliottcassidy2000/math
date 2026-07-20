---
id: THM-1395
title: "THE INDEX-1 COLLAPSE IS A PROPERTY OF THE GROUP, NOT OF THE TORUS — (ℤ/2)ⁿ escapes it and restores codimension n. (0) FIRST, the obvious escape is CLOSED: ind is monotone under ℤ/2-maps, so any invariant subspace A ⊆ T^k on which the involution is free has ind(A) ≤ ind(T^k) = 1 — in particular T^k contains no invariant free S^m for m ≥ 2, and no subspace restriction can revive Borsuk–Ulam. (I) But THM-1385 is a statement about ℤ/2, i.e. about ONE character. Let G = (ℤ/2)ⁿ act on T^k by half-translations h₁..hₙ independent over 𝔽₂ (n ≤ k). Every hypothesis of THM-1385 still holds — the action is free, the quotient is the aspherical T^k, Γ is torsion-free — yet the Fadell–Husseini index does NOT collapse: it is ker(𝔽₂[x₁..xₙ] → Λ(e₁..e_k)), which contains every xᵢ² (that IS the ℤ/2 collapse) but NOT the product x₁⋯xₙ, since an exterior algebra kills squares and not products of distinct generators. Hence there is no G-map T^k → S(V) for V the sum of the n distinct characters, so every G-equivariant f: T^k → V has a ZERO — codimension n, not 1. (II) The reachable codimension is exactly k, the 𝔽₂-rank of the 2-torsion, verified by exact rank computation. (III) WHY mac-mini's witness escapes and why distinct characters cannot: (cos 2πt₁, sin 2πt₁) carries the SAME character on both coordinates, so its target is really S¹ and an equivariant map exists; with distinct characters the target is a genuine (S^{n-1}, free) and the product class obstructs. Verified 120/120 certified zeros with distinct characters, against a VALIDATED control that does find zero-free maps when the characters are repeated"
status: >
  (0) PROVED — one line from monotonicity of the index, and it closes the first escape
  anyone would try.
  (I) PROVED modulo standard Fadell–Husseini index facts (Index_G(S(V)) = (x_1⋯x_n) for
  a sum of distinct nontrivial characters; monotonicity of the index under G-maps).  The
  computation of Index_G(T^k) is exact: H*(T^k/G;𝔽₂) is an exterior algebra and the
  classifying map sends x_i to independent degree-1 classes.
  (II) VERIFIED-EXACT by 𝔽₂ rank computation over all half-period subgroups, k ≤ 4.
  (III) VERIFIED: 120/120 random (ℤ/2)²-equivariant maps T² → ℝ² with distinct characters
  carry a zero certified by a non-zero winding number on a grid cell.  INSTRUMENT
  VALIDATED against the repeated-character control, which does produce zero-free maps
  (1/120 by random search, plus the explicit witness) — a search that could not find the
  escape where the escape exists would have proved nothing (MISTAKE-196).
  This does NOT resolve LRC(14) or any named open problem.  It reopens a route THM-1385
  closed, and says exactly how far that route can go.
source: kind-pasteur-2026-07-20-S128c107 (owner: one involution that is free and carries an odd map; the k-torus of the resonance lattice needs the ℤ/2-index form)
depends_on:
  - THM-1385    # mac-mini: free involution on an aspherical space has ind = 1
  - THM-1380    # opus: freeness and oddness sit on different involutions
related: [THM-581, THM-584]
script: 04-computation/elementary_abelian_escape_kps_S128c107.py (+ .out)
---

# THM-1395 — the collapse is about the group

## 0. The obvious escape, closed

The first thing one tries after THM-1385 is to restrict to a smaller, non-aspherical
invariant piece of the torus. It cannot work:

> **Corollary.** If `A ⊆ T^k` is invariant and the involution acts freely on `A`, then
> `ind(A) ≤ ind(T^k) = 1`.
> *Proof.* The inclusion `A ↪ T^k` is a ℤ/2-map and the index is monotone. ∎

In particular **`T^k` contains no invariant free `S^m` with `m ≥ 2`.** Recording this
because it is the natural next move and it is dead.

## I. The escape that works

THM-1385 proves `ind = 1` for a free ℤ/2-action on an aspherical space, and the proof is
correct. But `ind` here is the ℤ/2- (Conner–Floyd/Yang) index — a statement about **one**
character. Take instead

> `G = (ℤ/2)ⁿ` acting on `T^k` by translation by half-periods `h₁,…,hₙ`
> independent over `𝔽₂` (so `n ≤ k`).

Every hypothesis of THM-1385 still holds: the action is free, the quotient `T^k/G ≅ T^k`
is aspherical, `Γ` is torsion-free. What changes is the invariant. The Fadell–Husseini
index is the ideal

```
Index_G(T^k) = ker( 𝔽₂[x₁,…,xₙ] ⟶ H*(T^k/G; 𝔽₂) = Λ(e₁,…,e_k) ),   x_i ↦ w_i
```

with the `w_i` independent degree-1 classes. An exterior algebra kills **squares** —
`w_i² = 0`, so `x_i² ∈ Index_G`, and *that is exactly the ℤ/2 collapse* — but it does not
kill **products of distinct generators**: `w₁⋯wₙ ≠ 0`. Since
`Index_G(S(V)) = (x₁⋯xₙ)` for `V = ⊕ᵢ ℝ_{χᵢ}` the sum of the `n` distinct nontrivial
characters, we get `Index_G(S(V)) ⊄ Index_G(T^k)`, hence no `G`-map `T^k → S(V)`:

> **Theorem.** Every `G`-equivariant `f : T^k → V` has a zero.
> Concretely: if `f = (f₁,…,fₙ)` with `fᵢ(t + h_j) = (−1)^{δᵢⱼ} fᵢ(t)`, then `f` vanishes
> somewhere. **Codimension `n`, not 1.**

## II. How far it goes

`w₁⋯wₙ ≠ 0` in `Λ(𝔽₂^k)` exactly when the `hᵢ` are independent, so the reachable
codimension is precisely `k`, the `𝔽₂`-rank of the 2-torsion `(T^k)[2] ≅ (ℤ/2)^k`.
Verified by exact `𝔽₂` rank computation: independent sets of size `n` exist for `n ≤ k`
and not for `n = k+1`, at `k = 2,3,4`.

## III. Why the sharpness witness escapes, and why distinct characters cannot

mac-mini's zero-free map `f(t) = (cos 2πt₁, sin 2πt₁) : T^k → ℝ²` is odd for
`h = (½,0,…,0)`. Reproduced here at `k = 1,2,3,5,8` (`|f| ≡ 1`, oddness error `~10⁻¹⁵`).
It escapes for a specific reason: **both coordinates carry the same character**, so the
target is really the circle `S¹` with the antipodal action, and `T^k → S¹` equivariant
maps exist in abundance (any suitable linear character).

With `n` **distinct** characters that route is unavailable, and the product class
obstructs. Tested at `n = k = 2` with `h₁ = (½,0)`, `h₂ = (0,½)`: `f₁` supported on
Fourier modes with `m₁` odd and `m₂` even, `f₂` the mirror.

| character pattern | maps sampled | zero-free found |
|---|---|---|
| **distinct** `(1,0), (0,1)` | 120 | **0** — every one has a zero certified by non-zero winding |
| **repeated** `(1,0), (1,0)` | 120 | 1, plus the explicit witness |

The control matters more than the result. A search that finds no zero-free map proves
nothing unless it *can* find one where one exists; here it does, on the repeated-character
family, so the empty distinct-character column is evidence rather than an instrument
artifact.

## What this means for the resonance torus

The owner's caveat — `T^k ≠ S^k` blocks plain BU, use the index — is right, and
THM-1385's answer of `1` is right *for ℤ/2*. The correction is that the collapse is not a
property of the torus at all. The torus has plenty of `𝔽₂`-cohomology; what it lacks is
**squares**. Any argument that only ever squares one class sees a 1-dimensional object;
an argument that multiplies `k` distinct classes sees all `k` dimensions.

So the design rule for a Borsuk–Ulam attack on the resonance lattice: **do not look for a
single involution with a large index — look for a large elementary abelian 2-group of
half-period translations and give each coordinate of the test map a different character.**

## Named next

- Identify the resonance lattice's actual half-period group: its `𝔽₂`-rank is the ceiling
  on any BU-style codimension, and it is computable directly from the speed set.
- The odd-`p` analogue: `H*(T^k; 𝔽_p)` is exterior on degree-1 classes tensor a polynomial
  part, so `(ℤ/p)ⁿ` may behave differently — worth one computation before assuming it
  mirrors `p = 2`.
- Whether any *non-free* action does better via the Fadell–Husseini index of the fixed-point
  set — the one direction neither THM-1385 nor this theorem touches.
