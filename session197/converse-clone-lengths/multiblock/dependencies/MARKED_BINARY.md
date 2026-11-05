# Complete marked self-converse binary unit-clock rigidity

Session 197. This packet classifies one precise physical cloning mechanism at every order. It does not classify unrestricted Hamilton-path values. The relevant self-converse symmetry must fix the cloned root.

## 1. Statement

Let a nonempty tournament T have an anti-involution tau fixing a root r. Replace r by a transitive m-vertex module, with its old exterior incidences, and denote its Hamilton-path count by H_m. Suppose that, for real A>0,

    H_m = A*2^m + 1

holds for every sufficiently large positive integer m.

Then A=4^j for a unique integer j>=0, and the identity holds for every positive m. The rooted source has the following complete description:

- Its unique nontrivial strong component is the tournament on 2j+3 naturally ordered vertices with every arc forward except the last-to-first arc. The root is its middle vertex.
- A transitive prefix and suffix of the same length surround this component, with all arcs forward between these three blocks.

Conversely every such rooted source works. Thus every nonnegative j and padding length is realized. In particular the strong source for each amplitude is unique up to rooted isomorphism.

For every positive m, its nontrivial strong component has 2j+m+2 vertices and count

    H_m = 2^(m+2j)+1.

More generally, no root-fixed self-converse source of crossing matching rank one can give 7p for any prime p at any positive clone length, even without a unit offset. For the classified unit-offset family, seven never divides 2^d+1. If its count is prime, m+2j must be a power of two; this does not assert that all such Fermat numbers are prime.

## 2. Matching rank one and its exact symmetric form

The inherited root-run expansion writes H_m as a finite sum of integer exponentials. Its largest base is 1+nu, where nu is the maximum matching size in the directed crossing bipartite graph from the root out-neighbors O to its in-neighbors I. Its leading coefficient is positive. Comparing eventual growth with A*2^m+1 forces nu=1 and makes A an integer. Linear independence of the finitely many exponential bases makes the identity exact for every positive m.

An anti-involution fixing r cannot fix another vertex: it would reverse the arc between two fixed vertices. Thus it pairs O with I. Label

    O={o_1,...,o_t}, I={i_1,...,i_t}, tau(o_u)=i_u.

The tournament induced on O is some tournament S; the induced tournament on I is its converse under this pairing. Let C_(u,v)=1 mean o_u -> i_v. Applying tau shows C_(u,v)=C_(v,u), so C is symmetric.

A symmetric zero-one matrix of matching rank one has exactly one nonzero entry, and that entry is diagonal. An off-diagonal one would give two disjoint matching edges (u,v),(v,u). Two diagonal ones also give a size-two matching. At least one entry is nonzero because the rank is one.

Consequently there is a distinguished a in S such that every crossing arc goes I -> O except the single arc o_a -> i_a. The root lies between I and O: I -> r -> O. This proves the claimed form of all rank-one sources, without assuming strongness.

## 3. A square formula for every rank-one source

For an arbitrary nonempty S with distinguished vertex a, define

    h=H(S),
    q=sum_P H(S-V(P)),

where P ranges over all nonempty directed simple paths in S ending at a. The empty tournament has H=1. Equivalently q counts a cover of S by two labelled paths, the first nonempty and ending at a, the second possibly empty. In particular q is a positive integer.

**Rank-one formula.** For the source of Section 2 and every m>=1,

    H_m = h^2 + q^2*2^m.                                (1)

Use the inherited odd-cycle-collection identity: H is the sum of 2^(number of cycles) over collections of vertex-disjoint directed odd cycles, including the empty collection. Every cycle that is not wholly internal to I or O must use the single backward crossing arc o_a -> i_a. Therefore a cycle collection has at most one such external cycle.

Collections without an external cycle contribute H(I)H(O)=h^2. The clone chain is acyclic and contributes no internal cycle.

Remove the backward arc from an external cycle. The resulting directed path starts at i_a and ends at o_a. It consists of:

1. A nonempty path in I starting at i_a.
2. An arbitrary subset of the m clone vertices, in their transitive order, possibly empty.
3. A nonempty path in O ending at o_a.

The first and third paths cannot leave and return to their own sides: every remaining inter-block arc is forward in I -> clone chain -> O. Conversely every such choice of odd total length gives exactly one directed cycle, counted without a cyclic starting-point multiplicity because its unique backward arc identifies the decomposition. There is one invalid even-length boundary choice: both side paths are singletons and the clone subset is empty, so closing would require the absent arc i_a -> o_a and produce a forbidden two-cycle. Its total length is two, so the odd-parity filter below never includes it. Every other empty-clone choice uses a different, present I-to-O joining arc.

Fix the I and O paths. The external cycle is odd exactly when the selected clone-subset cardinality has one prescribed parity. For m>=1 there are exactly 2^(m-1) subsets of either parity. The external cycle has odd-cycle weight two. The remaining vertices of I and O independently contribute H of their induced tournaments. Summing over the two paths gives

    2 * 2^(m-1) * q * q = 2^m q^2.

The two path sums agree by applying tau and reversing a path. This proves (1), including mixed or empty residual induced tournaments. No formula for deletion at m=0 is claimed: the parity-balancing step needs a nonempty clone chain.

### An all-order 7p obstruction before imposing unit offset

Formula (1) also excludes 7p for every rank-one source in Section 2. The nonzero squares modulo seven are 1,2,4, and every 2^m is itself a nonzero square modulo seven. Thus if seven divides h^2+2^m q^2, neither h nor q can be nonzero modulo seven: otherwise their squared ratio would be the nonsquare -2^m. Hence seven divides both h and q, and 49 divides H_m. If H_m=7p with p prime, necessarily p=7 and H_m=49; but positive h,q divisible by seven give H_m>=49(1+2^m)>=147. This is impossible.

This obstruction includes arbitrary nontransitive S and weak sources. It does not use the global impossibility of Hamilton-path counts 7 or 21. The two elementary norm forms are h^2+(2^(m/2)q)^2 for even m and h^2+2(2^((m-1)/2)q)^2 for odd m.

## 4. Unit offset forces the complete physical classification

Compare (1) with A*2^m+1. Then A=q^2 and h^2=1, so H(S)=1. A tournament has H=1 exactly when it is transitive: the forward order gives the unique path in a transitive tournament, while a nontransitive tournament has a directed triangle and hence contributes at least the empty collection and that triangle in the odd-cycle identity.

Order S transitively and let j be the number of vertices preceding a. A path ending at a consists of a followed in reverse description by any subset of its j predecessors, in increasing order. There are 2^j choices, and every residual induced tournament has H=1. Thus q=2^j and A=4^j.

Write t=|S|. The t-j-1 vertices after a in O form a transitive terminal block; their paired I vertices form a transitive initial block of the same size. The remaining core has natural order

    i_a, i_(a-1), ..., i_1, r, o_1, ..., o_(a-1), o_a,

where the notation refers to the transitive rank of a. All its arcs are forward except o_a -> i_a. It is strong: the natural consecutive order followed by that backward arc is a directed Hamilton cycle. Its order is 2j+3 and r is the middle vertex. These core and padding blocks are its exact strong-component decomposition, apart from the singleton components inside each padding block.

Replacing r by a nonempty transitive chain retains the same backward arc and changes the core order to 2j+m+2. The square formula gives the displayed count. This establishes uniqueness and the converse construction.

For the prime restriction, if d=m+2j has an odd divisor v>1, write d=uv. Then 2^d+1 is divisible by 2^u+1, with a complementary factor greater than one. Hence d must be a power of two. Modulo seven the powers of two are 1,2,4, none equal to minus one, so no output is divisible by seven, including the exceptional possibilities H=21 or H=49.

## 5. Sources, boundaries and controls

The matching/maximal-run theorem and positive transitive-clone expansion are inherited from S196 spectral moves and the already reviewed Session197 `../positive-clock/PROOF.md` (SHA256 `b725136569ee84ddd46c5ef4f0b45b40e216eec9c314ffb604609e16a1e96f6f`). The odd-cycle-collection identity is the same inherited identity used and traced in `../PROOF.md`. The single-backward-arc count is also consistent with the reviewed one-seam family in `../../session197-lowcycle/SELF_CONVERSE_SEAMS.md` (SHA256 `eb47d9a4f5eb86c7f34480711c5306eecda081f262baf9fea1537952252a6995`). The present all-order classification follows from the symmetric matching-rank-one form and square identity proved here. No literature-priority claim is made.

The source can be weak; equal transitive padding is explicitly included. The root-fixed anti-involution is essential to the argument and is not weakened to unmarked self-converseness. The source of order one has matching rank zero and constant H_m=1, so it is excluded by A>0. There is no empty-source or m=0 claim.

`verify.py` checks all symmetric crossing matrices through size five, every physical rank-one source obtained from a labelled S through order four and every distinguished a, direct Hamilton counts for three positive clone lengths, the lifted anti-involutions, and the complete unit-offset classification. Its checks use explicit exceptions and are replayed normally and under Python -O. Finite tests are controls for the universal proof, not substitutes for it.
