---
theorem_id: THM-131
title: D_{14} irreducible representation decomposition of H(T_7) = 189
status: PROVED (computational)
proved_by: opus-2026-03-12-S58
date: 2026-03-12
related_theorems: [THM-126, THM-127]
tags: [paley, dihedral, representation-theory, hamiltonian-paths]
---

## Statement

The 189 Hamiltonian paths of the Paley tournament T_7 form a permutation
representation of D_{14} = ⟨r, s | r^7=s^2=1, srs=r^{-1}⟩, which decomposes
into irreducible representations as:

```
V_{HP} ≅ 18·ρ₀ ⊕ 9·ρ₁ ⊕ 27·ρ₂ ⊕ 27·ρ₃ ⊕ 27·ρ₄
```

where:
- ρ₀ = trivial representation (dim 1)
- ρ₁ = sign representation (dim 1): trivial on rotations, -1 on reflections
- ρ₂, ρ₃, ρ₄ = 2-dimensional representations (k=1,2,3)

Dimension check: 18·1 + 9·1 + 27·2 + 27·2 + 27·2 = 18 + 9 + 162 = 189 ✓

## Orbit Structure

**Z_7 orbits:** 27 orbits of size 7 (= H/p = 189/7)
- 9 self-paired orbits (s maps orbit to itself)
- 9 swapped pairs (s exchanges two orbits)

**D_{14} orbits:** 18 (by Burnside: (189 + 0 + 0 + 0 + 7·9)/14 = 252/14 = 18)

**Character on conjugacy classes:**

| Conjugacy class | Size | χ_V |
|----------------|------|-----|
| e | 1 | 189 |
| r, r^6 | 2 | 0 |
| r², r^5 | 2 | 0 |
| r³, r^4 | 2 | 0 |
| s, sr, ..., sr^6 | 7 | 9 |

All non-identity rotations fix 0 HPs (since p=7 is prime, no HP is periodic).
Each reflection fixes exactly 9 HPs (= number of self-paired orbits × 1 fixed
point per self-paired orbit via rotation shift).

## Proof

The D_{14} action on HPs is:
- Rotation r: path (v_0,...,v_6) ↦ (v_0+1,...,v_6+1) mod 7
- Reflection s: path (v_0,...,v_6) ↦ ((-v_6) mod 7,...,(-v_0) mod 7)

The reflection maps T_7 to T_7^op (THM-127), and since T_7 ≅ T_7^op, the
reflected-reversed path is also an HP of T_7.

Character computation:
- χ(e) = 189 (identity fixes all)
- χ(r^j) = 0 for j=1,...,6 (no HP fixed by non-trivial rotation, since p prime)
- χ(s) = 9 (computed by enumeration: exactly 9 HPs satisfy s(P) = P)

Multiplicity formula m_i = (1/|G|) Σ_C |C|·χ_V(C)·χ_i(C)*:
- m(ρ₀) = (189 + 63)/14 = 252/14 = 18
- m(ρ₁) = (189 - 63)/14 = 126/14 = 9
- m(ρ_k) = (378 + 0)/14 = 378/14 = 27 for k=1,2,3

## Interpretation

**18 D_{14}-invariant HP combinations:** These are HP counts/combinations that
are symmetric under both rotation AND reflection. They correspond to the 18
orbits of D_{14} acting on HPs.

**9 sign representations:** Rotation-invariant combinations that flip sign under
reflection. These detect the "chirality" of HP structure.

**Equal 2-dim multiplicities (27 each):** All three 2-dimensional irreps appear
with the same multiplicity. This is a consequence of the Galois action on the
character table: the 2-dim irreps are permuted by Gal(Q(ω_7)/Q), and since
H(T_7) is rational, all must appear equally.

## Generalization Conjecture

For Paley T_p (p ≡ 3 mod 4 prime), the HP permutation representation under
D_{2p} decomposes as:

```
V_{HP} ≅ a·ρ₀ ⊕ b·ρ₁ ⊕ c·(ρ₂ ⊕ ... ⊕ ρ_{(p-1)/2})
```

where a + b = H/p (number of Z_p orbits) and all 2-dim multiplicities are
equal (= c = H/p by Galois rationality). Then:

a + b + 2c·(p-1)/2 = H
a + b = H/p
c(p-1) = H - H/p = H(p-1)/p

So c = H/p. And a = (H + H/p · (# fixed points of s)) / (2p) ... needs more
careful analysis.

For T_7: a=18, b=9, c=27. Check: 18+9 = 27 = H/p ✓. c = 27 = H/p ✓.

## Scripts

Script: `04-computation/d14_irrep_decomposition.py`
Output: `05-knowledge/results/d14_irrep_decomposition.out`
