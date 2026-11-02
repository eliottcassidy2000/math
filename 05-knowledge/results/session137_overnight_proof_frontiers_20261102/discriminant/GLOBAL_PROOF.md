# A vertex-cycle induction for the Hamiltonian determinant bound

**PROVED, independently audited; exact computational controls supplement the all-order proof.** This manuscript concerns nonempty tournaments. The full tournament *value-spectrum* conjecture is a different problem and remains open.

## The target and inherited ingredients

For a loopless tournament with adjacency matrix A on n>=1 vertices, let H count directed Hamiltonian vertex words. Put

    u(T)=det(I+A),   v(T)=(-1)^n det A,
    K=A-A^T,        D(T)=det(I+K)/2^(n-1),
    s(T)=1^T(I+K)^(-1)1.

The conclusion is

    H(T) >= u(T)+|v(T)| = max(1,s(T)) D(T) >= D(T).       (1)

The stronger middle quantity is the strong-component base left open in `THM-1950-h-ge-disc-reduced-to-strongly-connected.md` and identified with rooted parity energies in `THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average.md`. The proof below does not assume either open inequality. It uses the inherited odd-cycle collection formula (OCF), ordinary determinant expansion, and elementary skew-matrix algebra. OCF remains credited to Grinberg–Stanley and its cited matrix-algebra presentations; no claim of literature priority is made for this bound without a separate source review.

For the empty tournament define H=u=v=1. Inequality (1) is FALSE for the empty tournament and is never asserted there. The weaker H>=u is valid at that base.

## 1. Algebraic positivity, with no path-count hypothesis

Let M=I+K and J=11^T. Tournament complementation gives

    M+J=2(I+A),     M-J=-2A^T.

Since K is real skew, det(M)>0. If z=M^(-1)1, then

    s=1^Tz=z^TMz=||z||²>0.

The rank-one determinant lemma therefore gives

    u=D(1+s)/2,    v=D(1-s)/2,
    u+v=D>0,      u-v=sD>0.                              (2)

Thus u>|v| for every NONEMPTY induced tournament. At the empty tournament, u+v=2 and u-v=0, still nonnegative. This positivity is established algebraically before the induction; it is not inferred from the desired H bound.

## 2. Exact vertex-cycle recurrences

Let T be the old tournament, add a vertex x, and let S be its out-neighborhood among old vertices. Let P range over all simple directed paths in T whose first vertex lies in S and last vertex lies outside S. Because these endpoints belong to disjoint sets, P has at least two vertices. Let E and O denote paths with an even and odd number of vertices, respectively, and write R_P=T-V(P).

The three exact identities are

    u(T+x) = u(T) + sum_E u(R_P) - sum_O u(R_P),
    v(T+x) = -sum_(E union O) v(R_P),
    H(T+x) = H(T) + 2 sum_E H(R_P).                      (3)

**First identity.** In det(I+A), a permutation is a collection of disjoint directed cycles of length at least two, with unused vertices contributing diagonal1. Separate collections avoiding x from those containing one cycle through x. Deleting x from that cycle yields exactly P. If P has ell vertices, its new cycle has length ell+1 and determinant sign(-1)^ell. The remaining collection contributes u(R_P). No cyclic-rotation multiplier occurs because x is distinguished.

**Second identity.** The signed spanning-cycle-cover quantity v(T)=(-1)^n det A assigns weight(-1) to each cycle of a spanning cover. The new distinguished cycle therefore contributes a factor-1 regardless of its length, and the remaining spanning cover contributes v(R_P). This also includes an empty remainder, where v=1.

**Third identity.** OCF assigns weight2 per odd cycle in a vertex-disjoint odd-cycle collection. Collections avoiding x contribute H(T). A cycle through x is odd exactly when P has even vertex count, and the remaining collection contributes H(R_P). This proves the last identity.

## 3. Induction and both margins

Induct on the number of vertices for the WEAKER assertion H>=u, starting at the empty tournament. Suppose it is known for every smaller tournament. Subtract (3). For either sign epsilon in{+1,-1},

    H(T+x)-u(T+x)+epsilon v(T+x)
      = H(T)-u(T)
        +sum_E [2(H(R_P)-u(R_P)) + u(R_P)-epsilon v(R_P)]
        +sum_O [u(R_P)-epsilon v(R_P)].                 (4)

Every summand is nonnegative: induction handles H-u, and (2) handles u-epsilon v, including the explicitly checked empty remainder. Therefore both margins in (4) are nonnegative. In particular H(T+x)>=u(T+x)+|v(T+x)|>=u(T+x), proving the induction assertion at the next size and simultaneously proving (1) for every nonempty tournament.

There is no assumption about strongness, feedback number, parity of n, or a finite computational cutoff. The earlier two-feedback-vertex finite-kernel route is valid as a specialization but is unnecessary for the global proof.

## 4. Equality and strictness in the ordinary tournament case

A transitive tournament has H=u=1 and v=0, so equality holds throughout. For C3, (H,u,v,D)=(3,2,-1,1), hence H=u+|v| but H>D.

Every nontransitive tournament contains a directed triangle. If its order is at least four, choose x on that triangle. The other two triangle vertices form an even path P in (3), and its remainder is nonempty. Equation (2) makes u(R_P)-epsilon v(R_P) strictly positive for BOTH signs epsilon. Thus both margins in (4) are strictly positive.

Consequently

    H=u+|v|  iff  T is transitive or is C3;
    H=D      iff  T is transitive.                      (5)

These equality statements are for nonempty unweighted tournaments only.

## 5. The exact proper-cycle defect

For an ordinary tournament let C_all count all directed cycles (length at least three, modulo cyclic rotation), and let C_proper count those not spanning all vertices. Both count different cyclic orders on the same support separately.

First the same induction gives

    delta(T):=H(T)-u(T) >= C_all(T).                     (6)

Indeed C_all(T+x)=C_all(T)+sum_P 1. The recurrence for delta has an even-path summand 2(H(R_P)-u(R_P))+u(R_P) and an odd-path summand u(R_P). Every u(R_P) is a positive integer, including u(empty)=1, so each new cycle contributes at least one.

In either signed margin (4), delta(T) counts every cycle avoiding x. Every path with nonempty remainder contributes at least one because u(R_P)-epsilon v(R_P) is a positive integer. These paths correspond precisely to proper cycles through x. Empty-remainder terms are nonnegative. Consequently

    H(T) >= u(T)+|v(T)|+C_proper(T).                    (7)

If T contains two vertex-disjoint directed triangles, choose x in one and use the even path consisting of its other two vertices. Its remainder contains the other triangle, so delta(R_P)>=1. Its summand therefore contributes at least three instead of the one already allocated to that proper cycle. Both signed margins exceed C_proper by at least two.

Conversely, if no two directed triangles are vertex-disjoint, no two directed cycles are vertex-disjoint: each induced tournament on a directed cycle contains a triangle. Let o,e count odd/even cycles, and c_n count Hamiltonian cycles. Since there are no multiple-cycle collections, determinant expansion and OCF give

    u=1+o-e,  v=-c_n,  H=1+2o.

Thus H-u-|v|=o+e-c_n=C_proper exactly. We obtain the sharp classification

    equality in (7) iff no two directed triangles are vertex-disjoint;
    two disjoint triangles imply H-u-|v|-C_proper >=2.  (8)

This strengthening is unweighted; the step u(R)>=1 uses integrality and must not be transferred to fractional matrices.

## 6. Weighted extension and its equality cases

Let A have zero diagonal and real entries satisfying0<=A_ij<=1 and A_ij+A_ji=1 for i!=j. Define H(A) as the sum of products of arc weights over Hamiltonian vertex words. It equals the expectation of H in the independent random tournament with edge probabilities A_ij. Taking expectations in OCF gives its weighted version: vertex-disjoint odd cycles use distinct unordered edges, so every collection's expectation is its product weight.

The identities (2) still hold. In (3), each path now carries the nonnegative weight of its old arcs and the two incident arcs through x. Paths of ONE old vertex are allowed in the determinant recurrences: they represent weighted two-cycles x-i-x. These paths belong to O, hence do not occur in weighted OCF, which only uses odd cycles of length at least three. They contribute only the beneficial nonnegative O terms in (4). All longer paths are treated exactly as before. The same induction therefore proves

    H(A) >= det(I+A)+|det A|                             (9)

for every nonempty weighted tournament matrix. For n>=3, any fractional pair {x,i} gives a positive one-vertex-path weight A_xi A_ix. Its remainder has n-2>=1 vertices, so its odd-path summand in each signed margin is strictly positive by (2). Hence fractional edges force strict H>u+|v|. Combining the deterministic classification gives the complete equality rule:

- n=1 or n=2: H=u+|v| for every weighted tournament matrix;
- n>=3: H=u+|v| exactly for deterministic transitive tournaments or the deterministic C3.

In order two, H=1, u=1-p(1-p), v=-p(1-p), and D=1-2p(1-p). Therefore H=D holds, in every nonempty order, exactly for deterministic transitive matrices. These low-order exceptions explain why the ordinary equality rule cannot be copied without adjustment.

The proof uses weighted OCF and weighted determinant recurrences directly. It does NOT infer the determinant inequality by Jensen or by taking expectations of a deterministic determinant: products A_ij A_ji vanish in deterministic samples but need not vanish in their mean matrix.

## 7. Verification and status boundaries

The universal argument is (2)–(4), not extrapolation from a census. Independent checks separately verify the determinant signs, path-endpoint orientation, every empty-complement contribution, weighted one-vertex path terms, and both margins. Literal path counting and independent exact determinants verify these small hostile examples; separate modular-Gaussian and Bareiss implementations audit the identities.

This result closes the repository's H>=disc target and its stronger P bound. It does not prove that every odd integer except7/21 occurs as H, does not settle LRC14, and gives no efficient exact Hamiltonian-path counter. The ordinary determinant remains a lower bound rather than a complete count invariant.

## Primary-source attribution

The OCF used in (3) is [Corollary20 of Irving–Omar, *Revisiting the Rédei–Berge Symmetric Functions via Matrix Algebra*, EJC32(4) (2025), #P4.43](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v32i4p43/pdf/). That corollary credits Grinberg–Stanley, Theorem1.39 and Lemma6.5, in [arXiv:2307.05569](https://arxiv.org/abs/2307.05569). The determinant recurrences and induction are proved explicitly here. A targeted source check found no matching inequality in the sources inspected; this is not a comprehensive novelty or priority claim.
