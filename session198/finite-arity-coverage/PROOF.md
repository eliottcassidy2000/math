# A finite-arity arithmetic coverage problem with physical witnesses

Session198. **PROVED:** internal sink runs can be absorbed into a gadget index; an explicit finite-arity arithmetic family has physical C1/C2 witnesses with at most floor(log2 n)+10 vertices. **FINITE-EXACT:** this restricted family covers every odd n from3 through100,001 except the eight values below; six literal repairs leave only7 and21. **OPEN:** the stated finite-arity number-theoretic range lemma. The finite check does not establish that lemma or an all-target theorem.

## 1. Explicit arithmetic functions and matrices

For L>=1 define three integer functions by

    (g(1),Q(1),B(1))=(1,1,1),
    (g(2l),Q(2l),B(2l))=(g(l),Q(l),2B(l)),
    (g(2l-1),Q(2l-1),B(2l-1))
      =(g(l)+l-(l mod2),Q(l)+B(l),2B(l)), l>=2.

In particular B(L) is the least power of two at least L. Define the fixed-size integer map

    T_L(D,E,U)=(g(L)D-(L mod2)E+LU,
                 D-[L=1]E,
                 Q(L)D+B(L)U).                           (1)

These are the complete physical source gadget maps projected onto (D,E,U), as proved in the reviewed source grammar. The omitted deletion count is not needed by any map(1). Every finite composition in(1), beginning at(1,1,1), therefore has an actual source graph. No numerical-cone test is being substituted for physical membership.

The other explicit source family is

    w(L)=(g(L),L mod2,L,1), L>=2.                       (2)

It is the converse of the actual one-apex gadget G(L), rooted at its apex. The old root-ending count of G(L) is one, so its converse is start-one; the other entries in(2) are the reviewed one-apex profile identities.

The inherited proofs are the source grammar at commit6f3f93a3915ef2a7eba9a611f67a1bb75bd392c9, the complete C2 bilinear theorem frozen with it, and the physical one-apex family in session198/source-pair-modular at commitc71ff2d7a3bc2bd6121457a35d122a92f0c2bb32. Frozen copies accompany the release.

## 2. Sink absorption is an exact physical-profile identity

Let O append a universal sink, retaining the root. For r>=1,

    O^r(D,E,U,d)=(D,0,2^r U,d).

The complete indexed map M_L additionally outputs d'=D. The binary recurrences give

    g(2^r L)=g(L), Q(2^r L)=Q(L), B(2^r L)=2^r B(L).

Since 2^r L is even and at least two, direct substitution proves

    M_L O^r(v)=M_(2^r L)(v), r>=1.                  (3)

This is equality of all four response coordinates of two actual constructions. A graph-isomorphism claim is not needed. It also holds at the singleton boundary. The r=0 case is simply the identity M_L=M_L and is treated separately.

Consequently every source profile built using at most four M moves and any number of O moves has a presentation

    O^r M_(L_k)...M_(L_1)(singleton), 0<=k<=4,r>=0, (4)

with no internal O moves. Repeatedly absorb each sink run immediately before the next gadget using(3). We do not claim that every source, or every C1 scalar, has such a bounded-gadget presentation.

For a source v=(D,E,U,d), its C1 closure has count F(v)=D-E+U. If r=0 this is the first coordinate of T_1(v). If r>=1,

    F(O^r(v))=D+2^r U,

which is the first coordinate of T_(2^r)(v). Thus the exact scalar family furnished by(4) and C1 closure is

    A = { first coordinate of
          T_(2^r) T_(L_k)...T_(L_1)(1,1,1)
          : 0<=k<=4, r>=0, L_i>=1 } minus {1}.       (5)

Every member of A is realized by a strong C1 graph. The excluded value1 is the harmless singleton/edge boundary. Formula(5) has at most four free gadget integers and one nonnegative exponent, not an unbounded word or an implicit physical-membership condition.

## 3. A second family involving only two integer arguments

For v=(D,E,U,d), w=(F,G,V,f), put

    N(v,w)=(D-E+U)F-(D-E+d)G+(D+3U-d)V-(U+E)f,
    M(v,w)=(D+2U)F-(U+d)G+(2D-E+6U-d)V-(U+E)f.

The complete C2 theorem realizes these counts on every pair of nontrivial actual sources. Define the pure two-integer arithmetic set

    P = { N(w(L),w(K)), M(w(L),w(K)) : L,K>=2 },      (6)

using the explicit vector(2). Thus every member of P has an actual strong C2 witness.

## 4. The remaining finite-arity number-theoretic lemma

Here is a sufficient, explicit arithmetic statement with no unrestricted source grammar in its quantifiers:

> **Finite-arity coverage lemma (OPEN).** Every odd integer n>1861 belongs to A union P, where A and P are exactly the finite-parameter expressions(5)-(6).

The proposed threshold1861 comes from finite evidence, not a proof of eventual coverage. The assertion remains open even though both families are physical, every parameter choice is allowed, and their union passes the check below.

**Conditional consequence.** If this arithmetic lemma holds, every odd integer n>=3 other than7 and21 is the Hamilton-path count of a strong tournament with at most three Hamilton cycles. For n>1861, use(5) or(6). For n<=1861, the finite certificate below supplies the appropriate recipe, with six explicit repairs. This is a precise sufficient reduction; it is not an equivalence claim between all C1/C2 tournaments and the bounded family(5)-(6).

## 5. A uniform graph-order bound for every represented target

The gadget M_L adds 1+ceil(log2 L) vertices. Its U output is at least L times the old U, because Q>0 and B>=L. Along the k arbitrary gadgets in(5), starting with U=1, therefore

    U >= product_i L_i.

For the final power-of-two gadget, the target n is at least 2^r U, including r=0 since D>=E. Hence

    2^r product_i L_i <= n.                         (7)

The direct physical graph in(5) has exactly

    k+2+r+sum_i ceil(log2 L_i)

vertices. Since k<=4, equation(7) bounds this by log2 n+2k+2<=log2 n+10, and thus by floor(log2 n)+10.

For(6), the one-apex source w(L) has 1+ceil(log2 L) vertices. The N construction adds one more vertex beyond its two sources, and M adds two. The reviewed lower bound n>=3LK+1 implies the same, looser floor(log2 n)+10 bound for both. The six fixed repairs also satisfy this bound by literal inspection. Consequently the conditional all-target result above would have a uniform logarithmic-order physical witness, with leading coefficient one and additive constant ten.

This size bound is unconditional for every parameter tuple already represented by(5)-(6); only existence of such a tuple for every large odd n is conditional.

## 6. Exact bounded verification, not extrapolation

A complete target-bounded enumeration of(5) uses only its at most four prefix maps. Put F=D-E+U. The reviewed monotonicity theorem gives strict increase of F under each M_L and the bound F(M_L(v))>=LU. Thus a prefix relevant to a final target n<=X has F<=X, and only L<=floor(X/U) need be tried. The terminal exponent r is tested exactly while D+2^r U<=X, with r=0 handled by F. Duplicate projected triples can be merged at their first, hence smallest, gadget depth because every future T_L and every terminal readout depends only on that triple.

For(6), symmetry permits L<=K, and the exact physical lower bound n>=3LK+1 makes the loop K<=floor((X-1)/(3L)) complete. The same weaker bound safely covers M as well as N.

The standalone verifier reconstructs every selected positive recipe with an independent bit-by-bit evaluation of g,Q,B. It checks the sink-absorption identity on a linear basis, verifies the order bound for every selected witness, and directly counts paths and cycles in each repair graph. The independently encoded normalized enumeration agrees with the earlier unnormalized O/M depth-limited probe; only the normalized standalone enumeration is part of this release.

At X=100,001, the odd values missing from A union P are exactly

    7,21,67,231,247,533,821,1861.

The six values other than7 and21 have the supplied literal repairs: C3 for67 and C2 for the other five. They are unchanged witnesses from the reviewed session198/c1-monotone packet. Hence every odd integer from3 through100,001 other than7 and21 is physically realized. This statement does not assert that7 or21 are impossible in unrestricted higher-cycle tournaments.

NORMAL.json and OPTIMIZED.json record the profile and certificate counts and must agree byte-for-byte. REPAIRS.json contains the six actual adjacency lists. The support hash pins the deterministically reconstructed positive certificates without a large support-table file. No cofinite or unbounded coverage claim is inferred from the stabilization of the missing list.
