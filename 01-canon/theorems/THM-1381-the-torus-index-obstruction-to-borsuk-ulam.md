---
id: THM-1381
title: THE TORUS INDEX OBSTRUCTION — a free translation involution on T^k has ℤ/2-index exactly 1 for EVERY k (explicit nonvanishing equivariant map f(x)=e^{2πix₁}); and for integer speeds the LRC resonance orbit is a closed circle, so the index is 1 for a second, independent reason. Hence Borsuk–Ulam-type arguments can force ONE constraint, never n. The free involution that distinguishes LRC from tournaments sits on a space too flat to exploit it.
status: PROVED (the witness map is explicit and verified to machine precision at k = 1,2,3,5,12; the cohomological argument is two lines; the circle-orbit fact is immediate for integer speeds). A NEGATIVE about a method family, not about LRC(14).
source: klein-2026-07-20-S332
depends_on:
  - THM-1043  # scale coordinates
related:
  - THM-1042  # the component-length obstruction (additive certificates)
  - THM-127   # dihedral anti-automorphism (the tournament-side involution, which is NOT free)
---

# THM-1381 — the torus index obstruction

klein-S322 found that the LRC involution `s ↦ −s` on `(ℤ/qℤ)*` is **free** for every `q ≥ 3`, while the
tournament complement `T ↦ T^op` has fixed points (the self-complementary classes). Freeness is exactly
the hypothesis Borsuk–Ulam-type arguments want, so it is natural to ask what equivariant topology buys on
the LRC side. This theorem answers: **nothing beyond one constraint**, and says precisely why.

## 1. The index is 1, for every k

**Theorem.** Let `c ∈ T^k` be a nonzero 2-torsion element and let `ℤ/2` act on `T^k` by `x ↦ x + c`.
This action is free, and its **ℤ/2-index is exactly 1**, independently of `k`.

*Proof.* Freeness gives index `≥ 1`. For the upper bound, write `c = (c₁,…,c_k)` with some `c_j = ½`;
relabel so `c₁ = ½`. Define

```text
f : T^k → ℝ²,     f(x) = (cos 2πx₁, sin 2πx₁).
```

Then `f(x + c) = (cos(2πx₁ + π), sin(2πx₁ + π)) = −f(x)`, so `f` is equivariant, and `|f| ≡ 1`, so `f`
never vanishes. An equivariant map into `ℝ² ∖ {0}` exists, hence index `≤ 1`. ∎

*(Verified numerically at `k = 1,2,3,5,12`: `max|f(x+c)+f(x)| ≈ 1.9e−15`, `min|f| = 1`.)*

The cohomological reading is the same fact: `H*(T^k; ℤ/2)` is an **exterior** algebra on degree-1
generators, so every degree-1 class squares to zero; in particular `w₁² = 0` for the classifying class of
any free translation action, and the index cannot exceed 1. **Dimension buys nothing.**

Contrast: on `S^k` the antipodal action has index `k`, and no equivariant map `S^k → ℝ^k ∖ {0}` exists —
that *is* Borsuk–Ulam. The torus fails precisely where the sphere succeeds, and the failure is in the ring
structure, not the dimension.

## 2. A second, independent reason, special to LRC

For **integer** speeds `v₁,…,v_n`, the map `t ↦ (v₁t, …, v_nt) mod 1` is invariant under `t ↦ t+1`, so it
descends to `S¹ → T^n` and its image is a **closed circle** — a 1-dimensional subtorus. The LRC problem
therefore lives on a 1-parameter space no matter how many speeds there are. Whatever free involution one
puts on it, the index is 1.

So the obstruction is doubly robust: it survives both the "work on the big resonance torus `T^{n−1}`"
route (§1) and the "work on the actual orbit" route (§2).

## 3. Consequence for the method family

An equivariant map `X → ℝ^m` is guaranteed a zero only when `m ≤ index(X)`. With `index = 1`:

> **Borsuk–Ulam-type arguments on the LRC resonance torus can force one constraint, never `n` of them.**

They cannot reach LRC(n) for `n ≥ 3`. This is the fourth delimited family in the map of what does not
work, alongside pairwise invariants (klein-S324), alternating truncations (klein-S325), and
additive/proportional certificates (THM-1042).

## 4. The synthesis this closes

S322 posed a puzzle: the LRC side has a *free* involution and the tournament side does not, yet the
tournament side is where the structural machinery lives. The resolution:

- **Tournaments**: involution with fixed points ⟹ a spine (SC classes) to build on, but no BU.
- **LRC**: involution free ⟹ BU is *applicable*, but the space is a torus ⟹ index 1 ⟹ BU is *empty*.

**Freeness is necessary but not sufficient; the space must also be curved enough to carry index.** The
very feature that distinguishes LRC from tournaments is defused by the flatness of the space it acts on.
That is why the fibration of S322 transports concepts and not proofs in *either* direction.

## 5. Scope

This bounds a method, not the problem. It says nothing about `M` for any family and contradicts nothing.
It does predict that any future proposal routed through equivariant topology on the resonance torus —
Borsuk–Ulam, Yang index, Fadell–Husseini, or the ℤ/2-index form — will yield at most one constraint, and
should be checked against that ceiling before being developed.

*Files: `04-computation/lrc_torus_index_klein_S332.py` (+ .out).*
