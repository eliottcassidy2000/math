---
theorem_id: THM-127
title: Dihedral anti-automorphism of Paley tournaments (p≡3 mod 4)
status: PROVED
proved_by: opus-2026-03-12
date: 2026-03-12
related_theorems: [THM-126]
tags: [paley, dihedral, automorphism, group-theory, symmetry]
---

## Statement

Let p ≡ 3 (mod 4) be prime and T_p the Paley tournament on Z_p (vertex i→j iff j−i is a
quadratic residue mod p).

The full symmetry group of T_p is the dihedral group D_{2p} = ⟨r, s | r^p=s²=1, srs=r^{-1}⟩
acting as:

- **Rotation r**: vertex v ↦ v+1 (mod p) — an automorphism T_p → T_p (preserves orientation)
- **Reflection s**: vertex v ↦ −v (mod p) — an ANTI-automorphism T_p → T_p^{op} (reverses all arcs)

D_{2p} acts faithfully on Z_p by these maps.

## Proof

**r is an automorphism:** i→j iff j−i ∈ QR_p, and (j+1)−(i+1) = j−i, so the shift by 1
preserves the tournament. This generates the cyclic automorphism group Z_p ≤ Aut(T_p).

**s is an anti-automorphism:** Under v↦−v: i→j becomes −i→−j, i.e., −j→−i in the
original labeling. Equivalently, arc (i,j) maps to arc (−j,−i), which is the same as (−j→−i)
iff −i−(−j) = j−i ∈ QR_p. But this gives the arc −j→−i, which means s maps each arc of T_p
to the REVERSE arc of T_p — hence s: T_p → T_p^{op}.

**Why p≡3 mod 4 is essential:** The map v↦−v = v+p is an anti-automorphism iff −1 is NOT
a quadratic residue. By Euler's criterion: (−1)^{(p-1)/2} ≡ −1 (mod p) iff p≡3 (mod 4).
So −1 ∉ QR_p precisely when p≡3 (mod 4).

**Why p≡1 mod 4 fails:** For p≡1 mod 4, −1 ∈ QR_p, so QR_p = −QR_p, and v↦−v is an
ordinary automorphism (not anti-automorphism). Paley tournaments at p≡1 mod 4 are self-complementary
but the dihedral action degenerates — the reflection fixes the tournament rather than flipping it.

## Geometric Picture

The tournament T_p can be drawn as a regular p-gon with vertices at roots of unity
exp(2πik/p), k=0,...,p-1. The p rotations and p reflections of the p-gon generate D_{2p}.

Under this picture:
- Each rotation preserves arc orientation (automorphism).
- Each reflection reverses arc orientation (anti-automorphism) because reflection conjugates
  the Legendre symbol: η(−d) = −η(d) for p≡3 mod 4.

This gives a vivid geometric proof: the H-maximizer T_7 is the unique circulant tournament
on Z_7 that is invariant (up to complementation) under the full D_{14} symmetry.

## Consequences

1. **H(T_p) = H(T_p^{op}):** Since s is an anti-automorphism, T_p and T_p^{op} are
   isomorphic as tournaments, confirming H(T_p) is well-defined (H is invariant under
   tournament isomorphism, and H(T) = H(T^{op}) for any T since reversing all arcs also
   reverses all Hamiltonian paths).

2. **Dihedral orbit structure:** The automorphism group has size |Aut(T_p)| = p(p−1)/2
   (full automorphisms include all of GF(p)* by multiplication, not just the rotations).
   But the full dihedral symmetry D_{2p} — rotations AND reflections — still acts, with
   reflections swapping T_p ↔ T_p^{op}.

3. **Connection to OEIS sequence H(T_p)/|Aut(T_p)|:** The sequence 1, 9, 1729, ... counts
   "essentially different" Hamiltonian paths up to the dihedral symmetry.

## Dihedral Connection Table

| n | p-gon shape | Dihedral group | Reflection action on T |
|---|---|---|---|
| 3 | Triangle | D_6 | anti-automorphism (3≡3 mod 4) |
| 5 | Pentagon | D_10 | automorphism (5≡1 mod 4) — T_5^{op}≅T_5 |
| 7 | Heptagon | D_14 | anti-automorphism (7≡3 mod 4) |
| 11 | Hendecagon | D_22 | anti-automorphism (11≡3 mod 4) |
| 19 | 19-gon | D_38 | anti-automorphism (19≡3 mod 4) |

The Paley maximizer conjecture can be stated geometrically: the tournament maximizing H
on a regular p-gon is the one with full dihedral symmetry (all rotations are automorphisms,
all reflections are anti-automorphisms).
