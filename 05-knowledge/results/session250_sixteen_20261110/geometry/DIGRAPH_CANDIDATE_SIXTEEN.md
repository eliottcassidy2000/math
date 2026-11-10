# A literal 1-5-4-5-1 reconstruction from three-vertex digraphs

**Status: PROVED-IN-NOTE reconstruction; FINITE-EXACT symmetry table. NOT recovered historical provenance.** The user's original referent has not been established. This is a mathematically natural candidate, separate from the Clebsch switch construction.

## 1. The exact five bins

Take loopless simple directed graphs on three vertices, allowing both arcs on a pair, modulo vertex relabeling. Let e be the number of arcs. Their arc-count polynomial is

D_3(z)=[(1+z)^6+3(1+z^2)^3+2(1+z^3)^2]/6
      =1+z+4z^2+4z^3+4z^4+z^5+z^6.

Proof: identity, transposition, and3-cycle of S3 have arc-position orbit lengths respectively1^6,2^3,3^2. A fixed graph includes either all or none of each orbit. Average their generating polynomials.

The complementary five bins are

empty(e=0) / sparse(e=1 or2) / half-full(e=3) / dense(e=4 or5) / full(e=6),

with cardinalities **1,5,4,5,1**. This partition retains whether the graph is empty/full, and otherwise the sign of e−3. It deliberately forgets the one-versus-two-arc distinction and its complementary four-versus-five-arc distinction. Thus it is not the fine arc grading.

## 2. Complement and transpose are different symmetries

Complement C toggles all six arc positions. Transpose R reverses every arc. They commute with each other and with relabeling, giving a Klein four action on the16 types. Complement swaps the outer bins and sparse/dense bins; transpose preserves every bin and preserves H by reversing every directed Hamiltonian path.

Exact enumeration of all64 labelled graphs, canonized under all six vertex permutations, gives

Fix(identity)=16, Fix(R)=10, Fix(C)=4, Fix(CR)=2.

Therefore there are13 transpose-orbits,10 complement-orbits,9 complement-transpose-orbits and8 Klein-four-orbits. **All four half-full types are individually self-complementary.** They are: the transitive tournament, the directed3-cycle, a digon with one outgoing third-vertex arc, and its transpose. The first two are self-transpose; transpose exchanges the latter pair. In the sparse two-arc grade, the four types are outgoing star, incoming star, directed path, and digon plus isolated vertex. Their complements form the four-arc grade.

These assertions are class statements, not labelled fixed-point statements: no labelled digraph equals its bitwise complement. They retain the necessary relabeling witness implicitly via canonicalization.

The16 classes have H multiplicities {0:5,1:4,2:3,3:2,4:1,6:1}, obtainable directly from the supplied rows. H is not determined by the five bins: the half-full bin contains H=1 and H=3. Complement is not an H symmetry: empty/full give0/6. Thus neither coarse bins nor Klein orbits can replace the Hamiltonian-count sidecar.

## 3. All-order edge grading and twisted traces

Let X_n={(i,j):i!=j}, N=n(n−1). For pi in S_n and epsilon in{0,1}, let q_(pi,epsilon) act by relabeling and then, if epsilon=1, swapping the two coordinates. Let its cycle lengths on X_n be L_1,...,L_r. Define

F_(pi,epsilon,0)(z)=product_a(1+z^(L_a)),

F_(pi,epsilon,1)(z)=0 if any L_a is odd,
                     2^r z^(N/2) otherwise.

These are exact arc-count polynomials of labelled digraphs fixed by relabeling pi, transpose^epsilon and complement^delta together. Proof: without complement the membership bit is constant around each q-cycle; with complement it must alternate, possible precisely for even cycle lengths, with two choices and exactly half the positions occupied.

Consequences:

D_n(z)=1/n! sum_pi F_(pi,0,0)(z),

D_n(z)=z^N D_n(1/z).

The ordinary digraph-class sequence is D_n(1). Adjoining transpose gives edge-graded class polynomial (1/(2n!))sum_(pi,epsilon) F_(pi,epsilon,0)(z). Adjoining both transpose and complement gives the unweighted Klein quotient

K_n=(1/(4n!))sum_(pi,epsilon,delta) F_(pi,epsilon,delta)(1).

**Grading boundary:** the analogous average at general z is not an ordinary arc-graded quotient polynomial, since complement does not preserve arc count. It is a twisted fixed-point trace; using it as an orbit generating polynomial would be a type error. One may instead grade complement-orbits by min(e,N−e), but that is a different enumerator.

These formulas re-express the inherited pair-state Burnside mechanism with the actual arc grading and the commuting complement/transpose actions. They are not claimed as new classical enumeration formulas. Their use here is to make the candidate16 partition structurally meaningful and prevent conflating class-level and labelled fixed points.

## 4. True marking is not differentiation

Let M_n(z) count isomorphism classes (D,a), where a is a distinguished PRESENT ordered arc, graded by e(D). If pi has f fixed vertices, it fixes b=f(f−1) ordered arc positions. All other arc cycles have lengths greater than1. Burnside gives

M_n(z)=1/n! sum_pi b z (1+z)^(b−1) product_(arc cycles L>1)(1+z^L),

where the summand is zero when b=0. The marked arc must be fixed individually; it is not enough that it lies in a longer cycle of pi.

At n=3 this is simply M_3(z)=z(1+z)^5, with coefficient total32. In contrast zD_3'(z) has total48. The first counts Aut(D)-orbits of present arcs in each D; the second weights each unmarked class by its raw number of arcs. The five-bin reconstruction forgets precisely the action data needed to compare these objects.

For Hamiltonian paths the corresponding loss is milder only because an entire ordered spanning path fixes every vertex. The path-marked class mass is sum_[D] H(D)/|Aut(D)|, not sum_[D] H(D). This distinction links the candidate16 reconstruction to A000568-style class enumeration without claiming the digraph category or its H-spectrum is the tournament category.

## 5. Reproduction and boundaries

Run `python3 05-knowledge/results/session250_sixteen_20261110/geometry/digraph_candidate.py`; the matching JSON lists all16 canonical masks, arc counts, H, automorphism orders and three symmetry images. Bit positions are explicitly ordered pairs in lexicographic order. Enumeration is all64 labelled graphs, with no filters.

Inherited provenance: `05-knowledge/results/digraph_metagraph_s20dm.out` already records the fine grade1,1,4,4,4,1,1. The new candidate reconstruction packages that grading into the requested five bins and makes its involutions, loss ledger and marked-count distinction explicit. The later `digraph_matrix_boundary_20261110` packet contains the broader pair-state/matrix boundary and marked-path framework; it remains the larger inheritance source.

Open: whether this was the user's intended historical decomposition; whether its coarse bins have an additional intrinsic characterization in that original work; and any inference from these digraph classes to the complete tournament Hamiltonian spectrum. None is asserted.
