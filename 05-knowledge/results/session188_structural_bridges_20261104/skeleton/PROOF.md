# Rooted skeleton transform and the three-state gate

Let Q be a finite tournament with distinguished vertex v, and let Q[TT_r] replace v by a transitive tournament on r labelled vertices, with every external arc inherited from v. Permit r=0, which deletes v, and use H(empty)=1 for the empty path. Define B_k as the number of words consisting of all vertices of Q-v, each once, and k indistinguishable placeholders, subject to: no placeholders are consecutive; every consecutive pair of exterior vertices is a directed edge of Q; and every exterior/placeholder adjacency is oriented as if the placeholder were v. Set B_0=H(Q-v). There are finitely many B_k, and they are nonnegative integers.

For any tournament A, let q_k(A) count ordered partitions of its vertex set into k nonempty directed paths, counting the order of the paths and each path's internal orientation. Set q_0(empty)=1. The same contraction argument first gives H(Q[A])=sum_k B_k q_k(A). This expresses the full scalar response to uniform module substitution in terms of its ordered path-cover profile. It is the familiar path-cover substitution identity; the special case below gives a compact exponential basis.

**Theorem.** For every r>=0,

    H(Q[TT_r]) = sum_k B_k k! S(r,k),
    sum_r H(Q[TT_r]) z^r/r! = sum_k B_k (exp(z)-1)^k,
    B_k = (1/k!) sum_{r=0}^k s(k,r) H(r),

Here S(r,k) is a Stirling number of the second kind, and s(k,r) is a signed Stirling number of the first kind. Thus the signed Stirling transform of H is nonnegative and divisible by k!; if d is the maximum k with B_k>0, the sequence is an integer linear combination of 1^r,...,d^r for r>=1, with a genuine d^r term of coefficient B_d. Consequently its ordinary generating function is rational and its recurrence has characteristic dividing product_{j=1}^d (X-j), after r>=1. The r=0 term may contain a separate 0^r contribution.

**Proof.** Given a directed Hamiltonian path in Q[TT_r], contract each maximal consecutive run of module vertices to a placeholder. The result is a skeleton counted by B_k. Conversely every skeleton and every ordered partition of the r module labels into k nonempty blocks uniquely determine a Hamiltonian path: list each block in the fixed transitive order. The number of ordered partitions is k! S(r,k), including the empty case. This establishes the formula. The exponential generating function is the standard exponential formula for ordered nonempty sets, and its polynomial in exp(z)-1 yields the remaining assertions. The inverse formula follows by substituting z=log(1+t) in the exponential generating function, or by the inversion identity for the two Stirling transforms. The leading exp(dz) coefficient is B_d, hence cannot cancel. QED.

**Three-state application.** With row bitmasks in vertex order 0,...,5, let the three roots be

    (16,5,33,39,46,3), (16,5,49,39,42,3), (18,4,49,39,42,3),

and v=5. Direct skeleton enumeration gives B=(9+2j,27+2j,21,4), j=0,1,2. Therefore

    F_j(z) = -1 + (-3+2j)e^z + 9e^(2z) + 4e^(3z),
    H_j(r) = 4*3^r+9*2^r-3+2j   (r>=1),
    H_j(0) = 9+2j.

All three states share their 2- and 3-entry skeletons; only B_0 and B_1 change, by 2 per state. The genuine 3^r contribution comes from the four three-entry skeletons. The fixed 2 spacing persists for all r>=0, including the deleted-root case. The cubic ordinary recurrence applies for r>=1. The full named cubic response is a different, richer invariant; this scalar formula neither replaces it nor implies a reusable return for arbitrary starting H.

`verify.py` independently enumerates skeletons and checks each inflated tournament using a subset dynamic program for r=0,...,6, with explicit runtime checks surviving Python -O. It confirms the data, not the universal theorem by itself. The transform is an elementary path-cover/skeleton identity and no novelty claim is made for it; the coefficients give a compact structural synthesis of the already recorded gate.

**Finite certification of an all-order spacing.** Let Q and Q' be rooted tournaments whose skeleton profiles vanish above D, and let c be an integer. If H(Q[TT_r])-H(Q'[TT_r])=c for r=1,...,D, then this equality holds for every r>=1. Indeed the difference minus c is a linear combination of 1^r,...,D^r; the D-by-D evaluation matrix at r=1,...,D is a Vandermonde matrix multiplied by a nonzero diagonal, and is nonsingular over the rationals. A direct sufficient-and-necessary coefficient criterion is B_k(Q)=B_k(Q') for k>=2 and B_1(Q)-B_1(Q')=c; B_0 is unrestricted because it only adjusts r=0. For r=0 spacing also to equal c, require B_0(Q)-B_0(Q')=c. These are scalar spacing tests, not tests of full named response or of closure under future gadget operations.
