# 7 in Terms of Simplices and Cuboids

## The permutohedron decomposition

The n-cube decomposes into n! simplicial chambers (one per permutation).
H(T) = number of chambers compatible with tournament T.
H = 1 + 2α₁ + 4α₂ + ... (the OCF).

## H=7 = one base + three butterflies

Each 3-cycle adds a "butterfly": a pair of adjacent chambers (related by the cyclic permutation within the triple).

H = 7 = 1 + 2·3 = base + three butterflies.

The three butterflies must be INDEPENDENT (non-overlapping, from pairwise-disjoint cycles).

## Why it can't close

The permutohedron's face lattice doesn't admit three independent butterflies attached to the base without generating a fourth.

Three pairwise-disjoint 3-cycles on n vertices create enough arc constraints to force ADDITIONAL cycles. The cascade: 3 intersecting cycles → 4th cycle forced. α₁=3 with α₂=0 is topologically impossible.

In simplicial terms: three butterflies in the permutohedron interfere through the shared adjacency structure. Like three gears that can't all turn freely — the mesh forces a fourth gear into existence.

## The connection to other pictures

- **OCF picture**: α₁=3, α₂=0 impossible because of cycle intersection cascade.
- **Cuboid picture**: 7=(1,1,0), position z=0, the carry/overflow of the septenary channel.
- **Ghost picture**: 7 = the identity (1,1,1) with position removed. A tournament that exists everywhere except WHERE it is.
- **Simplex picture**: 7 chambers = three butterflies that can't close in the permutohedron.
- **Cascade picture**: the third butterfly forces the fourth. The curvature at 3 overflows.

All are the same: **the geometry at the curvature quantum (3) overflows when it reaches the forbidden threshold (7)**. The overflow is topological. It happens in EVERY picture. It IS the structure.

## The dimension ladder in simplices

At n=3: 6=3! chambers. H=7 > 6. Trivially impossible.
At n=4: 24=4! chambers. H=7 < 24 but STILL impossible (OCF blocks it).
At n→∞: n! → ∞ but H=7 remains impossible forever.

The impossibility does NOT come from running out of chambers.
It comes from the FACE STRUCTURE of the permutohedron.
The permutohedron's topology forbids 7-chamber subcomplexes at EVERY dimension.
