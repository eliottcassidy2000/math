# The odd-cycle packing identity does not commute with isomorphism

There are two distinct ways to put path or cycle data on an unlabelled
tournament: count presentations with stabilizer weights, or quotient the
decorated objects themselves. The distinction matters in the metagraph.
This note uses the usual odd-cycle packing formula: H(T) is the number
of vertex-disjoint packs of directed odd cycles of lengths ≥3, with
each chosen cycle assigned one of two colors. Uncovered vertices do
not carry a color.

## A scaffold Burnside rule

An oriented scaffold F prescribes some arcs and may distinguish
vertices or components by marks. Its oriented completions are
tournaments containing exactly those prescribed orientations on
the prescribed pairs (all other pairs are free). Isomorphism of a
completion is required to preserve the scaffold and its marks.
Let G be the automorphism group of that scaffold. A permutation
can preserve a tournament only when it has odd order: an even cycle
in its permutation decomposes at half its length into some swapped
pair, which cannot preserve that pair's orientation. Conversely an
odd-order permutation is free of all pair-reversing powers, so it
permits one independent orientation per unordered-pair orbit not
already prescribed by F. Thus the exact number of decorated classes is
$$
 rac1{|G|}\sum_{g\in G,\, |g|\ {m odd}}
          2^{E(\lambda(g))-e_F(g)},
$$
where e_F(g) is the number of prescribed-arc orbits. For an odd
partition λ of the vertices into permutation cycle lengths,
$$
 E(\lambda)=\sum_i (\lambda_i-1)/2+
                         \sum_{i<j}\gcd(\lambda_i,\lambda_j).
$$
The first sum counts distance classes within each odd orbit. The
second counts the coupled rotations across two vertex cycles.
At F empty this is the standard unmarked partition sum for
A000568. For a fixed full directed Hamilton path its scaffold
group is trivial, giving exactly its metagraph presentation count.
For a Hamilton cycle its scaffold group is cyclic; now the same
calculation recovers the [cycle necklace counts](CYCLE_NECKLACES.md)
with their variable numbers of cuts. The scaffold formula does
not confuse marked orbit counts with the plain number of classes.

## Closed partition formula for colored packing classes

Give a packing type by r uncovered vertices and m_{ℓ,c} directed
cycles of each odd length ℓ≥3 and color c=0,1. Equal-color cycles of
the same length are *not* individually distinguished. Its scaffold
group is S_r times the product of the wreath products
C_ℓ wr S_(m_{ℓ,c}). Its prescribed arcs are in equivariant bijection
with the covered vertices by taking the tail, so e_F(g) is exactly
the number of permutation cycles on the covered vertices.

Here is the explicit partition sum without enumerating G. Choose an
odd partition λ_0 of r, with weight 1/z_(λ_0), as in A000568.
For each (ℓ,c) choose an odd partition α of m_(ℓ,c), with weight
1/z_α. Its parts a are the lengths of cycles of the permutation
of the equal-length equal-color components. For each such part
choose b|ℓ, with weight φ(b)/ℓ, and contribute ℓ/b copies of
vertex-cycle length ab to λ. Concatenate these copies and λ_0.
If u is the number of contributed parts (excluding λ_0), sum
the product of these weights times 2^(E(λ)−u).

Indeed after a round through a permuted cycle of a components,
their internal rotations combine into one element of C_ℓ.
This element is uniform under the wreath-product average;
there are φ(b) elements of order b. It splits the aℓ vertices
into ℓ/b cycles of length ab. Since ℓ is odd, the resulting
vertex lengths are odd precisely when a is odd. This proves
the formula, including repeated cycle parts and identical
colors. The unmarked enumeration and the metagraph's one-ring
enumeration occur as the two extreme scaffold choices.

## Consequence: the orbit-level path/packing correspondence fails exactly at symmetry

For a given tournament T let G=Aut(T). Its H(T) Hamiltonian paths
have free G-action, because an ordered full path fixes every
vertex. The colored packings have exactly the same raw
cardinality by the packing identity, but their orbit count is
$$
 {H(T)\over |G|}+{1\over|G|}\sum_{1\ne g\in G}
          |\mathrm{Fix}_{\rm colored\ packs}(g)|.
$$
The second term is zero exactly when T is rigid: every g fixes
the empty packing. For the directed triangle, one path orbit
and three packing orbits are already different. Thus a
path-to-odd-cycle-packing bijection at raw cardinality cannot
be equivariant in any nonrigid tournament. The missing
stabilizers are not an artifact of the choice of cycle roots.

This explains why swapping a path formula for a packing formula
after passing to isomorphism classes can silently change the
enumerated object. It also gives a positive replacement: either
keep stabilizer weights (giving equal raw orbit *mass*), or count
packing orbits by the colored scaffold sum above. This is a
refinement for decorated classes, not a new closed formula for
the joint unmarked statistic (H,C), nor evidence of path-spectrum
completeness.

Run `python replay_packing.py` for independent literal orbit tests.
