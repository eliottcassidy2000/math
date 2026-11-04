# A conversed child cannot survive full named-response equality

**Proved integration lemma, pending independent review.** This is a direct general corollary of the already reviewed Session188 compressed exposed-anchor theorem, not a new parity mechanism. It applies to arbitrary nonempty children, including singletons. Combining it with Session189 transitive-child quartic reconstruction removes that theorem's converse ambiguity under a full-response-collision hypothesis.

## 1. Pinned inputs and terminology

Use the physical observer F_T(a,c), with x->y, a_v=[v->x], and c_v=[v->y]. A child is a fixed physical module block. The two queried tournaments have the same blocks and names; no relabeling is allowed.

The precise prior input is [S188 returns_next/UNIFORM.md](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/returns_next/UNIFORM.md):

- U2 permits a partition into disjoint module blocks and remaining individual singleton names, with any chosen uncovered singleton anchor v. Identifying attachment bits within the module blocks preserves the whole-converse parity obstruction.
- U3 gives the explicit compressed battery below, at most q+2 actual queries for q blocks including singleton blocks.
- U4 lifts its contrast2 modulo4 through every fixed ambient by giving both ports the original core-root exterior incidences.

The unbounded whole-converse/ambient separation underneath it is [incoming E3--E4](https://github.com/eliottcassidy2000/math/blob/76b3b1a77d124100f03026420f6b821533cdec15/05-knowledge/results/two_port_parity_separation_20261104/PROOF.md). The two-old-name boundary, where U2's E1 identity is not asserted, is handled there directly.

The profile substitution used below is exact over the integers. If both ports see a child uniformly, then that child remains a module in the augmented graph. Replacing its interior by an equal-complete-ordered-cover-profile interior preserves every ordered-cover coordinate of the augmented core, and hence preserves each queried Hamilton count after substitution into any fixed ambient. Uniform means separately constant on that child for each port; the two ports may have different directions.

## 2. General conversed-child lemma

Let

    T=Q[D_1,...,D_k],      U=Q^op[D_1',...,D_k']

on the same nonempty physical blocks, with total core order at least2. Suppose every corresponding pair of children has identical complete ordered-path-cover profiles. Assume that for at least one index z,

    D_z'=D_z^op

as literal named tournaments. The quotient Q need not be prime. Embed either core at the same root of any fixed ambient R, leaving all exterior names, arcs, and root incidences fixed.

**Lemma.** A physical two-port query separates the two ambient tournaments by2 modulo4. It may be chosen nonuniform only on D_z and uniform on every other child. A source-dependent battery of at most

    |D_z|+k+1

actual queries works in every fixed ambient, with the prescribed exterior attachment bits. This count equals q+2 for q=|D_z|+k-1 compressed blocks. The selected individual query need not be the same across ambients.

**Proof.** Write A=D_z. Keep each physical name of A as a separate singleton block, and keep each D_i, i!=z, as one compressed module block. This is a legal partition for U2: every non-singleton block is an actual module of T, and any chosen v in A is an uncovered singleton. No condition on the internal orientation of A is needed; it need not be transitive, prime, or strongly connected. The other names of A are merely singleton blocks too.

Apply U2--U3 to the whole-core-converse comparison T versus T^op. It gives an actual query constant on every D_i, i!=z, with contrast2 modulo4. Give both ports the original core-root incidence to each exterior name. U4/E4 gives the same nonzero residue through the fixed ambient, multiplying the core contrast by odd H(R).

Now evaluate the actual U at that same physical query. Its A interior already equals the A interior of T^op. Its quotient is exactly Q^op, also that of T^op. For i!=z, the profiles of D_i', D_i and D_i^op agree. Because the ports are uniform on D_i, restore D_i' to D_i^op without changing the actual queried count. This is an integer equality, not a congruence and not a claim that arbitrary endpoint marginals agree. Perform the disjoint replacements one at a time. They transform U into T^op inside the ambient. Consequently

    F_ambient(U)(a,c)=F_ambient(T^op)(a,c),
    F_ambient(T)(a,c)-F_ambient(U)(a,c)=2 modulo4.

For ambient substitution, restoration holds in every augmented-core ordered-cover coordinate, not merely its local Hamilton coordinate. This justifies arbitrary fixed exteriors. QED.

A singleton A is explicitly allowed: A=A^op and its unique name is the uncovered anchor. The lemma then says that a profile-preserving quotient converse with a singleton child is always detected. At total core order2 the literal endpoint query from E3 replaces the E1 battery and still respects the compression.

## 3. Complete physical battery

For total core order at least3, let B be the q-block partition just described, and choose v in A. Let d be the number of other blocks C with v->C; this is a block count, not the physical outdegree.

Inside the core set c=e_v. For every other compressed block C include

    a=e_C          if C->v,
    a=e_C+e_v      if v->C.

Here e_C has value1 at every physical name of C, not at an abstract replacement vertex. Also include a=e_v when d is odd, and a=0 when q-1+d is odd. Include the endpoint query (a,c)=(e_v,0). Other core attachment bits are zero.

For every exterior name t, in every one of these queries set

    a_t=c_t=[t->core].

The physical edge remains x->y. Let b_v=(S_v(T)+E_v(T)) modulo2. If b_v=1, the endpoint query has core contrast2 modulo4. If b_v=0, the sum of the c=e_v battery contrasts is2 modulo4, so at least one member has that residue. Uniform compression was already proved in U2 before this battery was formed. This is why the anchor may lie inside A even though A as a whole is a module: A is deliberately not one of the compressed blocks.

There is no truncation of the original U comparison at degree3. The actual target and whole-core-converse queried counts are first proved equal by integer profile restoration. Only the resulting whole-converse comparison invokes its already proved parity identity. Arbitrary original-module comparisons may still have genuine quartics.

## 4. Original-module and transitive-child consequences

Suppose U is an original-module XOR target of T and F_T=F_U. The S188 hereditary theorem gives equal complete profiles on every fixed subset. At a highest changed original strong-tree node M, the source quotient is prime and the target quotient is its converse. Other highest changed nodes can be restored because they are disjoint and both ports imitate M on the exterior.

Apply the lemma at M. It follows that **no child of M can have target interior exactly equal to its source whole converse**. This assertion applies to every child shape and size, not only transitive ones. In particular M can have no singleton child under the collision hypothesis.

Now suppose a source child A is transitive and |A|>=5. Session189's internal-distance theorem first forces its target to be either A or A^op. The present lemma excludes the A^op branch. Therefore

    U[A]=T[A]

literally, with the same physical order.

If every child of such a highest changed prime node is transitive of size at least five, the only surviving candidate at that node is

    Q[A_1,...,A_k]  versus  Q^op[A_1,...,A_k],

with all child interiors unchanged. Thus the remaining candidate is the quotient-only converse, rather than an independent reversal cube of child interiors. This does **not** prove that the quotient-only comparison can be hidden or visible for every prime quotient. No full-response collision is exhibited.

## 5. Attribution and boundaries

The general lemma is already implicit in U2--U4 plus profile substitution. S188's [pair-root proof](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/prime_pair_root_family/PROOF.md), Section3, explicitly uses its pair-child special case; [ordered_pair_ambient/PROOF.md](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/ordered_pair_ambient/PROOF.md), Section3, supplies the fixed-ambient pair integration. The TT3 degree-tomography root-family proof uses the same final whole-converse restoration after identifying a reversed target triple.

The new conclusion here is the general explicit formulation and its integration with the newly proved transitive distance reconstruction. It is not a claim of a newly discovered exposed-anchor mechanism, arbitrary-prime closure, or all-pairs orbit reconstruction.
