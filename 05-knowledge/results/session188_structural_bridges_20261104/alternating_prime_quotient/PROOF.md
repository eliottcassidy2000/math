# Pair-child separation in alternating-hub prime quotients

**PROVED:** an unbounded second strong-prime quotient family satisfies the marked-pair separation criterion. With a pair child at any quotient vertex and every other child nontrivial, quotient converse with all named child interiors unchanged is visible to an actual two-port query. Other child sizes and orientations are arbitrary. The proof is elementary and does not use a finite inference of primality.

## 1. The quotient family

For r>=2 let W_r have chain names 0,...,2r-1 with every edge i->j for i<j, and a hub h=2r. The hub beats the even chain names and loses to the odd chain names.

This tournament is strong: the cycle

    0->1->...->2r-1->h->0

visits all names. It is also prime, meaning it has no proper module of size at least two.

To prove primeness, a module containing two chain vertices i<j must contain every chain vertex strictly between them: such a vertex sees i and j differently. Thus it contains the whole interval from i to j, including both parities, so the hub also distinguishes two of its members and must lie in the module.

A module containing the hub and one chain vertex i must contain the chain endpoints 0 and 2r-1. If i>0, vertex0 beats i but loses to the hub, so it must belong to the module; if i=0 this is already true. If i<2r-1, the last chain vertex loses to i but beats the hub, so it also must belong; if i=2r-1 it is already present. The two endpoints force every intermediate chain name into the module. Hence the module is all of W_r.

Any set of at least two vertices contains either two chain vertices or the hub and a chain vertex. These cases prove primeness. The same elementary argument also works at r=1, where the quotient is a directed triangle; the main statement below keeps its small-size boundary explicit.

## 2. What deleting a quotient vertex exposes

Fix z.

- If z is the hub, the remainder is the transitive chain.
- If z=0, the new first chain vertex1 is a source in the remainder: it beats every later chain name and the hub.
- If z=2r-1, the new last chain vertex2r-2 is a sink in the remainder: every earlier chain name and the hub beat it.
- If 1<=z<=2r-2, the two chain neighbors z-1,z+1 are a homogeneous pair in W_r-z, oriented z-1->z+1. No remaining chain name lies between them, and every other chain name sees both uniformly. They have the same parity, so the hub also sees both uniformly.

Thus deleting any quotient vertex exposes either a global ordered cut or an ordered two-vertex module. These assertions use the physical labels displayed above and hold for every r.

## 3. Inflating arbitrary child tournaments

Let T=W_r[A_i] and U=W_r^op[A_i] retain all child interiors and physical names. Suppose A_z={p,q} is a pair oriented p->q. The other children are arbitrary nonempty tournaments.

The marked-pair coefficient identity gives, for suitable disjoint unions A,B of exterior child blocks,

    [a_q a_A c_p c_B](F_T-F_U)=C_BA(R)-C_AB(R),
    R=T-A_z.

If M=A->B is a module of R, the independently proved ordered-module inequality is

    C_BA-C_AB >= H(R/M)
       [q_2(A)H(B)+H(A)q_2(B)+q_2(A)q_2(B)].        (1)

It is strictly positive when |M|>=3. Each coefficient is the signed sum of sixteen actual named queries, so a positive coefficient guarantees an actual response separator. No cubic or quartic term is dropped.

Apply Section2. If z is the hub or a chain endpoint, take A,B across the displayed global cut. Since r>=2 and every remaining child is nonempty, the exterior has at least four physical names. Thus |M|=|R|>=4 and (1) is positive, without any further child-size assumption.

If z is an interior chain vertex, take

    A=A_(z-1),  B=A_(z+1).

Their union is an ordered module in R. It suffices that |A_(z-1)|+|A_(z+1)|>=3, equivalently that at least one neighboring child is nontrivial.

**Theorem.** For every r>=2, these explicit conditions separate quotient converse in T and U. In particular, if all children are nontrivial and at least one is a pair, choosing that pair as A_z always works, wherever it occurs in the quotient.

The remaining children have no restrictions on size, internal tournament type, or ordered-cover profile. The conclusion is an exact integer distinction, not necessarily a residue-two distinction.

## 4. Small cases and the limit of the criterion

At r=1 the quotient is C3 and all three possible deletions give a two-vertex ordered quotient. The coefficient is positive if the two exterior children contain at least three physical names in total. If they are both singletons, the ordered-cover skew is zero. This is only a zero boundary of this particular quartic: it is not an exact response collision, and the exposed-anchor theorem already supplies other separators.

For an interior deletion at any larger r, if both neighboring children are singletons, the displayed ordered module has only two physical names and (1) gives zero. No failure of the complete response is inferred. The theorem does not silently replace an arbitrary larger selected child by a pair.

The result concerns the whole substitution and its quotient converse. Arbitrary strong outer ambient lifting is not claimed. Its dependencies are the exact grouped quartic in `../marked_pair_quotient/PROOF.md` and the cover injection in `../ordered_module_skew/PROOF.md`.

## 5. Verification

The checker tests the strong/prime definitions and every deletion structure through r=7, including r=1. It retains full sixteen-query controls for all root positions at r=1,2 with nontrivial children, larger quotient examples with arbitrary remaining sizes, and the zero quartic boundary at r=1 with singleton exterior children. The unbounded assertions follow from Sections1–3; finite checks validate indexing and implementation only.
