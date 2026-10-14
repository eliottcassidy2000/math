# Fermat primes, Platonic curvature, and the cyclotomic tournament norm

Research synthesis, 14 October 2026. Status: proved elementary lemmas and classical implications, with an exact finite census and explicitly open transfer targets. No claim of historical novelty. Canon consulted at snapshot `0f7b1b200b1b96a6ce61bbe304bb92c603515bbf`; no canonical theorem is changed.

## 1. The structural connection, including the literal 27

There are two distinct inputs:

- The user's literal list is `3,5,27,257,65537`, with `27=3^3`.
- Replacing 27 by 17 gives the five **known** Fermat primes `3,5,17,257,65537`. It is not known that this is the complete list. The five convex Platonic solids, by contrast, are a complete classification. [Brakke's explicit polygon constructions](https://www.kenbrakke.com/RulerAndCompass/big-gon.html).

The useful common structure is the group

`Gal(Q(zeta_p)/Q) = F_p^*`.

It controls both the quadratic construction of the regular prime-sided polygon and the multiplier symmetries of a cyclic tournament. At a Fermat prime it is a 2-group: that makes the polygon constructible, while forcing every nontrivial multiplier to break any tournament orientation. The two consequences come from the same group, rather than the two lists having length five.

Literal 27 is valuable as the contrasting case. A regular 27-gon is not straightedge-and-compass constructible, since its cyclotomic degree is `phi(27)=18`. But the field `F_27`, whose additive group is `C_3^3`, supports a Paley tournament with 1,053 symmetries and only two distinct nonzero skew eigenvalues. **The field F_27 and the ring Z/27Z are different objects.**

## 2. Orientation-count lemma and a characterization of Fermat primes

Let `G` act faithfully on a finite vertex set. A `G`-invariant tournament exists if and only if `|G|` is odd.

Proof. If `|G|` is even, Cauchy's theorem supplies a nonidentity involution, which swaps some two vertices and reverses their directed edge. Conversely, if an ordered-pair orbit met its reversed orbit, an element would exchange two vertices and have an even cycle, hence even order. For odd `|G|`, therefore, the ordered-pair orbits come in distinct reversal pairs. Choose one orbit in each pair to orient all edges. If there are `r` ordered-pair orbits off the diagonal, the number of invariant labelled orientations is `2^(r/2)`.

For an odd finite field `F_q` and `H <= F_q^*`, apply this to the affine group of translations and multipliers. Equivalently choose a skew connection set `S` with

`S disjoint (-S), S union (-S) = F_q^*`.

An `H`-invariant choice exists exactly when `|H|` is odd. In that case multiplication acts freely on `F_q^*`, and the exact number of labelled connection sets is

`2^((q-1)/(2|H|))`.

These are choices on the fixed labelled field, not isomorphism classes.

### Prime-order consequence

For an odd prime `p`, the following are equivalent:

1. `p` is a Fermat prime.
2. The regular `p`-gon is constructible with straightedge and compass.
3. Every vertex-transitive tournament on `p` vertices has automorphism group exactly `C_p`.
4. The multiplier action `F_p^*` on all skew connection sets is free.

For (1) and (2), the cyclotomic Galois group is cyclic of order `p-1`; a quadratic tower exists exactly when that order is a power of two. If `2^m+1` is prime, then `m` has no odd divisor greater than one, so `m` is a power of two. See the square-root-tower and constructibility theorems in [Elman, algebra text, sections 52 and 60](https://www.math.ucla.edu/~rse/algebra_book.pdf).

For (3), Burnside's prime-degree theorem says a transitive group of prime degree is either doubly transitive or a subgroup of `AGL(1,p)`. The doubly transitive possibility cannot preserve a tournament. Thus its automorphism group lies in `C_p semidirect F_p^*`. A tournament automorphism group has odd order. If `p-1` is a 2-power, only the translation subgroup remains. The classical group-theoretic input is [Müller, *Permutation groups of prime degree, a quick proof of Burnside's theorem*](https://arxiv.org/abs/math/0310200). Conversely, when `p-1` has an odd factor greater than one, the orientation-count lemma constructs a tournament with a nontrivial odd multiplier subgroup. This proves the converses in (3) and (4) as well.

### Exact isomorphism-class counts

At a Fermat prime, every prime-order vertex-transitive tournament is cyclic and has precisely its translation group as automorphisms. An isomorphism between two such tournaments conjugates this group to itself, hence is affine. Multiplier orbits are free and translations do not change the connection set. Therefore the number of isomorphism classes of vertex-transitive tournaments is exactly

`2^((p-1)/2)/(p-1)`.

The counts at the five known Fermat primes are respectively

`1, 1, 16, 2^120, 2^32752`.

This counts vertex-transitive tournament isomorphism classes, not all tournaments on `p` vertices and not labelled tournaments. Reversal introduces no additional quotient: the multiplier `-1` already identifies each cyclic tournament with its opposite.

## 3. The cyclotomic norm carrier: an exact link to the corpus's open bound

Let `p` be any odd prime, let `S` be skew as above, and let `T_S` have arcs `i -> j` exactly when `j-i in S`. Put

`K=A-A^T`, `alpha_S=sum_(s in S) zeta_p^s`, `lambda_S=1+2 alpha_S`.

### Proposition A: the determinant is an algebraic norm

The Fourier eigenvalues of `K` are zero and the `p-1` Galois conjugates of `lambda_S`. Indeed, on the Fourier vector `(zeta_p^(bi))_i` the eigenvalue is

`sum_(s in S) zeta_p^(bs) - sum_(s in -S) zeta_p^(bs)`

`=1+2 sum_(s in S) zeta_p^(bs)`.

Consequently the corpus's invariant is exactly

`disc(T_S) = det(I+K)/2^(p-1) = Norm_(Q(zeta_p)/Q)(1+alpha_S)`.

The right side is a positive integer. This proof is valid for every odd prime, without a Fermat hypothesis.

### Proposition B: the Galois stabilizer is precisely the multiplier stabilizer

For `a in F_p^*`, `sigma_a(lambda_S)=lambda_S` if and only if `aS=S`. To prove the nontrivial implication, subtract their coefficient polynomials. The difference has degree at most `p-1` and constant coefficient zero. If it vanishes at `zeta_p`, it is a multiple of `Phi_p=1+x+...+x^(p-1)`; its zero constant coefficient forces that multiple to vanish. Thus the connection sets agree.

At a Fermat prime, Proposition B and the free multiplier action imply

`Q(lambda_S)=Q(zeta_p)`.

Since `S` is skew, complex conjugation sends `lambda_S` to `-lambda_S`. If `sigma_a(lambda_S^2)=lambda_S^2`, characteristic zero gives `sigma_a(lambda_S)=lambda_S` or `-lambda_S`. The first case forces `a=1`; the second, after composing with conjugation, forces `a=-1`. Therefore

`Q(lambda_S^2)=Q(zeta_p+zeta_p^(-1))`.

In particular, every such cyclic tournament has all `p-1` nonzero skew eigenvalues distinct, and its squared eigenvalues generate the complete real cyclotomic subfield. The nonzero factor of its skew characteristic polynomial is irreducible over `Q` and even. Constructibility of its eigenvalues does not imply small algebraic degree.

### The remaining combinatorial inequality

The canonical [LEM-003](../../01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md) proves that automorphisms act freely on directed Hamiltonian paths. Hence at a Fermat prime

`H(T_S)/p = number of Hamiltonian-path orbits under translations`.

Every vertex-transitive tournament is regular, so the inverse response in [THM-1950](../../01-canon/theorems/THM-1950-h-ge-disc-reduced-to-strongly-connected.md) is `s=p`. Its still-open strong-base inequality, restricted to this sector, becomes the precise conjecture

`number of Hamiltonian-path orbits >= Norm(1+alpha_S)`.

This is a count-versus-norm statement with both sides fully identified. An orbit injection, or a positive expansion of the norm into bounded path-orbit contributions, would be real progress. Neither the norm identity nor constructibility supplies that injection. The factor `p` must not be exported to arbitrary strong tournaments: the known order-seven example with `H=31, disc=6, s=2/3` refutes that extension.

## 4. A finite-exact census and a useful failed simplification

The companion standard-library verifier enumerates all skew connection sets for `p=3,5,17`, reduces by the exact multiplier action, counts Hamiltonian paths by fixed-start subset dynamic programming, and computes determinants by integer Bareiss elimination. It independently computes the norm as the determinant of multiplication by `1+sum_(s in S)x^s` on `Z[x]/Phi_p`. All eighteen classes pass the count-versus-norm inequality.

At `p=17`, the unique class maximizing `disc` also minimizes `H/(17 disc)`:

- `S={3,9,10,11,12,13,15,16}`
- `H=13,551,423,911`
- `disc=24,617`
- `H/(17 disc)=797,142,583/24,617`, approximately `32,381.7924`

The carousel has `H=13,689,269,499` and `disc=1`; its Hamiltonian count matches LEM-004's recorded value. The broad safety margin at 17 is finite evidence, not an all-order proof.

There is a genuine arithmetic distinction from Paley flatness. For a regular tournament on `n=4t+1` vertices, set `B=KK^T`. If `i -> j` and the vertices have `c` common out-neighbours, then

`B_ij=4c-n+2 = 1 (mod 4)`.

The flat spectral target would have off-diagonal entries `-1`, which is impossible. Thus

`||B-nI+J||_F^2 >= 4n(n-1)`.

This floor follows term by term from `|B_ij+1|>=2`. Flatness is available for doubly regular tournaments at orders `3 mod 4`, not `1 mod 4`.

However, **minimizing this variance does not identify the hardest class for the bound**. At 17 the minimum is 2,176, attained by

`S={1,4,9,10,11,12,14,15}`, `disc=14,621`,

whereas the minimum-ratio/max-disc class above has variance 3,264. This exact counterexample rules out replacing the full norm by one spectral-variance statistic.

## 5. What the literal 27 supplies

Construct `F_27=F_3[t]/(t^3-t-1)`. The cubic has no root in `F_3` and is irreducible. Direct an edge `x -> y` when `y-x` is a nonzero square. Since `27=3 mod 4`, negation interchanges squares and nonsquares; this is a tournament.

Its explicitly verified automorphisms are

`x -> a x^(3^j)+b`, with `a` a nonzero square, `b in F_27`, `j=0,1,2`.

There are `13*27*3=1,053` distinct maps. They form its full automorphism group by the classical Paley-tournament theorem; the companion code independently verifies all 1,053 maps, without using that classification. The classical spectral identity is `KK^T=27I-J`; both facts are described in [Jones, sections 2 and 9.7](https://arxiv.org/html/1702.00285).

Thus its nonzero skew eigenvalues are `+3 sqrt(-3)` and `-3 sqrt(-3)`, each with multiplicity 13, and

`disc=28^13/2^26=7^13=96,889,010,407`, `s=27`.

The open strong-base target is `H >= 27*7^13=2,616,003,280,989`. The script does not count its Hamiltonian paths and does not claim to settle this instance. It certifies the matrix and the symmetry group. By LEM-003, `1053 | H`.

This contrast is structural: Fermat prime cyclic tournaments have trivial multiplier stabilizer and maximal cyclotomic spectral degree; the literal prime-power example has large odd stabilizers and highly repeated spectrum in a quadratic field.

## 6. Why there are five Platonic solids, and where the comparison stops

For a regular spherical polyhedron with `p` sides per face and `q` faces at a vertex, Euler's formula gives

`1/p+1/q-1/2=1/E > 0`.

With `p,q>=3`, the complete ordered list is

`(3,3), (3,4), (4,3), (3,5), (5,3)`.

This is why there are exactly five. The zero boundary contains `(3,6),(4,4),(6,3)`; negative values lead to the hyperbolic regime. Swapping `p,q` gives duality, so the five solids have only three rotation-group types: `A4,S4,A5`. A direct modern discussion of the Euler/Diophantine link is [Khudaverdian–Mkrtchyan](https://arxiv.org/abs/1604.06062).

The triangular, square, and pentagonal faces are all constructible. But constructibility of the regular 17-gon cannot make it a Platonic face: even the smallest allowed valence `q=3` gives `1/17+1/3<1/2`. The curvature inequality and the Galois 2-group condition impose different constraints.

Every Platonic rotation group has even order, hence cannot act faithfully as automorphisms of a tournament. Allowing a group element to reverse *all* arcs does not rescue these groups: the orientation-preserving kernel has odd order and index at most two, whereas `12,24,60` are divisible by four. Odd dihedral polygon groups do pass this weaker test: translations preserve any cyclic tournament and `x -> -x` reverses it.

A reciprocal-sum comparison with deficient/perfect/abundant integers must likewise specify its map. If a sphere is given cone orders equal to all divisors `d>1` of `n`, its orbifold Euler characteristic is

`chi_orb=2-tau(n)+sigma(n)/n`.

A perfect number therefore gives `chi_orb=4-tau(n)`, not generally zero. Perfect 6 gives `(2,3,6)` and zero curvature; perfect 28 gives `(2,4,7,14,28)` and Euler characteristic `-2`. The reciprocal-sum bridge is exact, but the three arithmetic signs are not the three curvature signs without a fixed number of cone points.

## 7. A precise HN seam transfer target, without changing the predicate silently

The HN seam construction supplies an undirected-graph involution induced by `sqrt(11) -> -sqrt(11)`. Its ability to preserve unit edges is compatible with the tournament obstruction: undirected graphs permit even-order automorphisms. A tournament-anti-automorphism involution can fix at most one vertex, since two fixed vertices would give an edge required to reverse itself.

A less restrictive candidate is an orientation **up to switching**. Let `K` be a skew signed adjacency matrix of an oriented graph, and let `P` represent a graph automorphism `tau`. Ask whether

`P K P^T = D K D`

for a diagonal sign matrix `D=diag(s_v)`. Define on each undirected edge the direction-independent ratio `r_ij=(P K P^T)_ij/K_ij`. Such a `D` exists if and only if the product of `r_e` around every graph cycle is `+1`.

Proof. Necessity telescopes the vertex signs around the cycle. For sufficiency, fix a root in each component and define `s_v` as the product of edge ratios on a root-to-v path; the cycle condition makes it path independent. This is a finite cycle-space/coboundary test.

Switching preserves `det(I+K)` because `I+DKD=D(I+K)D`. It does not preserve Hamiltonian-path counts or the unlabelled inverse response: switching the cyclic triangle gives a transitive tournament, changing `H` from 3 to 1 while `disc=1` remains fixed. The response can be retained only with its source vector tracked: `v^T(I+K)^(-1)v` is preserved under `(K,v)->(DKD,Dv)`.

Thus the test should be applied to an intrinsically specified orientation and a clearly chosen invariant. An arbitrary vertex-order orientation tests only that choice. Moreover, a sparse oriented unit-distance graph is not a tournament, so neither THM-1950's path inequality nor a coloring obstruction automatically transfers.

## 8. Sharp next questions

1. **Fermat-sector orbit/norm conjecture:** for every Fermat prime and skew `S`, is `H(T_S)/p >= Norm(1+alpha_S)`? This is a precise restriction of the corpus's open base, now supplied with an exact common field carrier; exhaustive at 3, 5, 17. The missing object is a count-preserving positive or injective comparison.
2. **Extremizer question, low-confidence:** does a maximum-norm class minimize `H/(p disc)` at every Fermat prime? It holds at the small classified orders, including all sixteen classes at 17. A variance-only replacement already fails there. This question should be attacked via the full cyclotomic norm, not inferred from one spectral moment.
3. **HN switching test:** after specifying a geometrically meaningful seam orientation, does its Galois involution give a trivial cycle-space obstruction? If so, preserve the determinant and the signed source vector explicitly, then identify a separate invariant that still detects the target coloring obstruction. Without the last step, the construction is a symmetry transfer, not progress on the HN chromatic predicate.

## Reproduction

`python3 04-computation/experiments/fermat_platonic_orientation_check.py`

The frozen JSON records all eighteen class representatives, exact Hamiltonian counts, determinant/norm agreement, variance values, and the F_27 certificates. Hamiltonian DP is independently checked against literal permutations at orders 3 and 5. The census scope is vertex-transitive tournaments at orders 3, 5, and 17; there is no claimed exhaustive census at 257, 65537, or 27, and no new solution of an open major problem.
