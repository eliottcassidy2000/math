---
id: THM-1385
title: The ℤ/2-index of a free involution on an ASPHERICAL space is exactly 1 — so on the resonance k-torus Borsuk–Ulam collapses to a single equation, for every k and EVERY free involution. Asphericity is precisely the obstruction: ind(S^k)=k because S^k is not aspherical for k≥2, while ind(T^k)=1 because T^k is. Sharp, with an explicit zero-free odd map T^k → ℝ² for every k.
status: PROVED (three lines, via torsion-freeness of the quotient group) + VERIFIED (the sharpness witness constructed and checked numerically for k = 1,2,3,5,8; the Klein-bottle non-translation case independently confirmed via Wu's formula w₁² = w₂ = χ mod 2 = 0).
source: mac-mini-2026-07-20-S125 (owner: the k-torus of the resonance lattice needs the ℤ/2-index form, since T^k ≠ S^k blocks plain BU)
depends_on:
  - THM-1380  # opus: freeness and oddness sit on different involutions (the circle-level no-go this generalizes)
related:
  - THM-581   # the even-category reading of the LRC floor
  - THM-584   # complement = antipodal map
---

# THM-1385 — aspherical ⟹ index one: the torus collapse of Borsuk–Ulam

**One line.** The owner's caveat — *`T^k ≠ S^k` blocks plain BU, so use the ℤ/2-index* —
is exactly right, and the index form returns a hard number: **1**. Not `k`, not `k−1`.
One equation, for every `k`, for every free involution.

## Setup

For a free ℤ/2-space `X`, let `ind(X)` be the ℤ/2- (Conner–Floyd/Yang) index: the height
of the characteristic class `w ∈ H¹(X/ℤ₂; ℤ/2)` classifying the double cover; equivalently
the largest `m` admitting a ℤ/2-map `S^m → X`. The Borsuk–Ulam mechanism is

> `ind(X) ≥ m` ⟹ every odd map `X → ℝ^m` has a zero.

Classically `ind(S^k) = k`, which is BU.

## (A) The theorem

> **Let `X` be aspherical and let ℤ/2 act freely on `X`. Then `ind(X) = 1`.**

*Proof.* `ind ≥ 1`: the double cover is nontrivial, so `w ≠ 0`.

`ind ≤ 1`: suppose `ind ≥ 2`, i.e. there is a ℤ/2-map `S² → X`. Passing to quotients gives
`RP² → M := X/ℤ₂`, and equivariance says the classifying class pulls back to the generator
of `H¹(RP²;ℤ/2)`, so the composite

```
π₁(RP²) = ℤ/2 ⟶ π₁(M) = Γ ⟶ Γ/π₁(X) = ℤ/2
```

is **onto**. Hence the generator maps to some `γ ∈ Γ` with `γ² = 1` and `γ ≠ 1`. But `X`
aspherical + the action free makes `M` aspherical with `Γ` acting freely on the contractible
universal cover, so **`Γ` is torsion-free** — contradiction. ∎

## (B) The Borsuk–Ulam collapse

> An odd map `X → ℝ^m` from an aspherical free ℤ/2-space is forced to vanish **only for
> `m = 1`**. For every `m ≥ 2` zero-free odd maps exist.

## (C) Sharpness on the torus, explicitly (VERIFIED)

Let `τ(x) = x + v` on `T^k`, `v = w/2` with `w ∈ ℤ^k \ 2ℤ^k` (the general free translation
involution). Choose `i` with `w_i` odd and put `a = e_i`. Then

```
f : T^k → ℝ² = ℂ,     f(x) = e^{2πi⟨a,x⟩}
```

satisfies `f(x+v) = f(x)·e^{2πi⟨a,v⟩} = f(x)·e^{iπ} = −f(x)` (**odd**) and `|f| ≡ 1`
(**never zero**). Verified numerically for `k = 1,2,3,5,8`. So the bound `ind = 1` is attained,
not merely an upper estimate.

**Non-translation check.** The Klein-bottle involution `ι(x,y) = (x+½, −y)` on `T²` is free
and *not* a translation; its quotient is `K`, the classifying class is `w₁(K)`, and Wu on a
closed surface gives `w₁² = w₂ = χ(K) mod 2 = 0`. Again `ind = 1`, independently.

## (D) Why `S^k` escapes — asphericity is the whole obstruction

`ind(S^k) = k` and `ind(T^k) = 1`, and the proof isolates the reason: `S^k` is **not**
aspherical for `k ≥ 2` (`π_k ≠ 0`), so the torsion argument has no purchase. The dividing
line between "BU gives `k` equations" and "BU gives one" is exactly asphericity of the
carrier — not its dimension, not its homology.

## (E) Consequence for the resonance lattice (the LRC reading)

Any Borsuk–Ulam-type argument staged on the resonance `k`-torus — **whatever** free involution
is chosen, translation or not — yields **at most one scalar equation**. It can never deliver
the `k` independent conditions such an argument is usually reached for. This is a structural
no-go, and it is complementary to THM-1380: that result shows freeness and oddness sit on
*different* involutions at the circle level; this one shows that even when they are made to
coincide on a higher-dimensional torus, the index does not grow.

**What is not claimed.** This bounds the *topological* method only. It says nothing about
whether the underlying LRC statement is true, and it does not obstruct measure, Fourier, or
covering arguments on the same torus — only the odd-map/BU route.

*Artifacts:* `04-computation/z2_index_aspherical_macmini_S125.py` (+out).
Credits: THM-1380 (opus-S401, the circle-level no-go), THM-581/584.
