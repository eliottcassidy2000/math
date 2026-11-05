# Unbounded erased-deletion ambiguity at a fixed order

Session197. A separate corollary of the reviewed ordinary one-seam law and weak-source decomposition. The frozen weak-oracle theorem is not changed.

For a rooted tournament write (h,e,R) for its Hamilton-path, root-ending, and root-starting ordered-two-cover responses, and d=H(T-root). The gadget M_L in `ROOTED.md` has

    h'=g(L)h-(L mod2)e+LR,  e'=h-[L=1]e,
    R'=Q(L)h+B(L)R.

For L=2^j, j>=1, its constants are g=Q=1 and B=L, so on a balanced input R=h it gives

    (h,e,h) -> ((2^j+1)h, h, (2^j+1)h).             (1)

It adds exactly j+1 vertices. Deleting the new root leaves the old input followed by a transitive tail, hence d'=h. These are actual gadgets, not a formal matrix construction.

**Theorem.** For every m>=1, put

    P_m = product_(j=1,...,m) (2^j+1),
    n_m = 2 + m(m+3)/2.

There are at least m rooted tournaments of order n_m, each with exactly one Hamilton path starting at its root and the same weak triple

    (h,e,R)=(P_m,0,2P_m),                            (2)

whose deletion counts are the m distinct integers

    d_j=P_m/(2^j+1),  j=1,...,m.                     (3)

Their initial strong components all have exactly one Hamilton cycle, and each graph has just one further singleton SCC.

**Proof.** Start from the singleton, whose profile is(1,1,1). For each j=1,...,m apply M_(2^j) once, in any order, choosing j to be the last operation for the j-th witness. Every step after the first has a balanced input by(1); the first step does as well. The physical rooted start-one grammar guarantees that the resulting nontrivial core is strong, has start count one, and belongs to the C=1 class. Repeatedly applying(1) gives core h=R=P_m, while its last root-ending and deletion responses are both P_m/(2^j+1).

Append one universal sink. The weak-source decomposition in `WEAK_START_ONE_ORACLE.md` shows that this preserves h and d, makes e zero, doubles R, and preserves the start count. Thus all witnesses have(2) and deletion response(3). Each gadget adds j+1 vertices, regardless of its position in the list, so their common order is1+sum(j+1)+1=n_m. The factors2^j+1 are distinct, and P_m is positive, so the m deletion counts are distinct. Different deletion counts also ensure distinct rooted isomorphism classes. QED.

This rules out any uniform finite bound on the number of deletion responses compatible with a weak triple, even when the order is known. It does not assert that the displayed set is exhaustive, nor that more than one extra coordinate is necessary: retaining d itself always distinguishes these examples. No unrestricted linear-rank claim follows.

The verifier constructs every arc for the m choices at m=1,...,4 and directly counts each weak graph, its root deletion, its core cycle count and its start response. It also checks the exact gadget response/order identities for m through50, without claiming direct graph-DP counts at those larger orders. The all-order claim rests on(1) and the physical grammar, not those finite tests.

Dependencies:
- ROOTED.md, SHA256 ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e.
- WEAK_START_ONE_ORACLE.md, SHA256 55048003e63b68ef21a51d424117145a4afaeb4bdfb43b5946d19524b6006ec0.
