# Unlabelled curl necklaces, primitive periods, and connected even graphs

Status: PROVED counting consequences of the companion curl classification; CITED equinumerosity theorem as explicitly identified below. No new natural graph bijection is claimed.

## 1. Necklace transform of A000568

Write A(z)=Σ_(n≥0)a_n z^n for unlabelled tournaments (a_0=1) and S(z) for nonempty strong tournament isomorphism types, with a singleton treated as strong. Unique ordered SCC decomposition gives A=1/(1−S) as an ordinary formal series.

The companion theorem says that isometry classes of the standard augmentation root lattice together with its tournament alternating form correspond exactly to cyclic words in strong tournament types, graded by total vertex count. A lattice isometry is ± a coordinate permutation; its global minus enters an alternating pairing twice and cancels. Thus this does not identify converse tournaments unless the appropriate cyclic strong-word is already isomorphic. The n=1,2 cases are checked directly.

Let c_n denote these nonempty unlabelled curl classes, C(z)=Σ c_n z^n. Burnside for rotations of a word, with each letter weighted by its vertex count, gives

C(z)=Σ_(d≥1) φ(d)/d · log A(z^d).

For completeness, length k rotations of cycle length d|k have k/d letter cycles. Their fixed-word series is S(z^d)^(k/d), and there are φ(d) such rotations. Averaging by k and writing k=dj yields Σ_d φ(d)/d Σ_j S(z^d)^j/j, which is the displayed formula.

It is invertible:

log A(z)=Σ_(d≥1) b(d)/d · C(z^d),  b(d)=∏_(p|d)(1−p), b(1)=1.

Indeed b is the Dirichlet inverse of φ: φ=id*μ, so φ^(−1)=1*(μ·id), whose divisor sum is the displayed product. These formulas transform exact counts; they are not a speed improvement without an independently faster algorithm for curl orbits.

At n=4, c_4=3 whereas a_4=4 and connected ordinary graph types number6. The labelled connected-graph identity must NOT be carried over unchanged to unlabelled objects.

## 2. Primitive curl classes recover the older Lyndon/even-graph bridge

A cyclic strong-component word is primitive if it is not a repetition of a shorter cyclic word. Let p_n count primitive words of total vertex weight n. Every cyclic word is a unique r-th power of a primitive word, up to cyclic rotation, so

c_n=Σ_(d|n)p_d,  p_n=Σ_(d|n)μ(d)c_(n/d).

Recover the older repository bridge THM-4721, Lemma E and Theorem A (`01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`). Here 'even graph' has the Royle–Praeger–Glasby–Freedman–Devillers meaning: every automorphism reverses an even number of edges relative to a chosen orientation. It does NOT mean every vertex has even degree.

Their published theorem says that unlabelled even graphs and tournaments have equal counts. Primary source: https://arxiv.org/html/2204.01947v2 , Theorem1.1; journal J. Algebraic Combin.57 (2023),515–524, https://doi.org/10.1007/s10801-022-01197-0 . Only this counting theorem is imported.

An even graph is a disjoint union of connected even graphs. Necessity follows by extending an odd component automorphism by the identity. For sufficiency choose the same orientation on isomorphic component copies. A permutation of those copies preserves the chosen orientations, and the remaining automorphisms inside components reverse even numbers of edges. Thus every global automorphism has even reversal parity.

If e_n counts connected even graph types, their multiset decomposition gives

A(z)=∏_(n≥1)(1−z^n)^(−e_n).

Substitution in the necklace formula yields

C(z)=Σ_(d,k,n≥1) φ(d)e_n z^(dkn)/(dk)
    =Σ_(m,n≥1)e_n z^(mn),

because Σ_(d|m)φ(d)=m. Therefore

c_n=Σ_(d|n)e_d, and hence p_n=e_n.

So primitive cyclic lists of strong tournament types are equinumerous with connected even graphs. This weighted-Lyndon/connected-even count equality and an order-dependent type codec were ALREADY proved in `05-knowledge/results/tournament_even_bijection_20261011c.md`, sections2–4, and reused in `tournament_atoms_abc_and_plane_20261011c.md`. Each primitive cyclic word has one Lyndon representative after choosing an alphabet order. The Session36 addition is the actual metric-curl-fibre interpretation, its necklace transform and its stabilizer/weighted-response consequences, not a new discovery of the older equinumerosity or codec. The count equality does not furnish a canonical or relabelling-equivariant map. Historical no-natural-bijection results are not bypassed.

## 3. Symmetry records the lost initial gap

For a fixed labelled tournament with k strong components, its curl stabilizer inside S_n may permute those blocks only by a cyclic rotation. Let r be the number of rotations preserving their cyclic isomorphism-type word, equivalently its repetition multiplicity. The kernel of this action is the product of the internal strong-component automorphism groups. Every type-preserving rotation is realizable by choosing isomorphisms between the respective blocks, so

|Stab_(S_n)(curl T)| = r · ∏_i |Aut(S_i)|.

Each tournament automorphism group has odd order: an even-order group would contain an involution, whose transposed vertex pair reverses its own arc. Consequently the curl stabilizer has even order exactly when r is even. Primitive cyclic words have no extra rotational stabilizer, but odd nontrivial repetition can also leave the stabilizer odd.

This is a statement about the permutation stabilizer of labelled curl, not the full metric-lattice automorphism group: for n≥3 the latter also contains the global central −I acting trivially on the alternating form.

The fundamental distinction is now visible: a tournament orders its strong blocks along a line; full curl retains their oriented circle and forgets a starting gap. The extra symmetry is exactly the periodicity made possible by forgetting that gap.

## 4. Hamiltonian-spectrum refinement and its necessary Adams operation

For a fixed real parameter s, define the ordinary weighted isotype series

A(z,s)=1+Σ_(T nonempty tournament types) H(T)^s z^|T|,
C(z,s)=Σ_(curl types [T]) H(T)^s z^|T|.

SCC multiplicativity again gives A=1/(1−S) for the weighted strong-type series. A rotation whose letter cycles have length d repeats every selected block d times; its Hamiltonian weight is therefore raised to the d-th power. Weighted Burnside consequently gives

C(z,s)=Σ_(d≥1) φ(d)/d · log A(z^d, d s).

The simultaneous replacement s→ds is essential. Using only z→z^d would be wrong once a repeated nontrivial strong component appears. The first such possibility is two copies of C3, at order6, with H=3²=9 rather than3.

The same identity can retain the ENTIRE Hamiltonian-count frequency distribution by working in the rational monoid algebra with basis symbols [h], multiplication [h][k]=[hk], and Adams operation ψ_d([h])=[h^d]. Then the formula is Σ_d φ(d)/d · ψ_d(log A(z^d)). The weighted inverse uses the same coefficients b(d)/d as above and the same Adams operations, since ψ_d ψ_e=ψ_de.

This is an exact transform of multiplicities, not a characterization of which H-values are realizable. It supplies a cross-check on future enumeration tables and prevents confusing path-count weights with an additive marking variable.
