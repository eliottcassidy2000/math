# Twelve physical full-query states detect all first-order collision motion

**Status: PROVED linear-algebra theorem in every dimension at least three; FINITE-EXACT, reproducible physical nine-coordinate bank with sharp cardinality and radius.** This is a concrete continuation of the planar-collision and tournament-interface research. It does not construct a Keller map, physical tournament deformations, or source-closed ternary offset returns. No literature-priority claim is made.

## 1. Source, idea, lever, obstacle

The historical source is [THM-4411](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/01-canon/theorems/THM-4411-first-order-collision-transgression-seminormal-tradeoff.md): on a smooth surface, tangent-relation periods exactly detect first-order collision compatibility. [The later quadratic-sidecar theorem](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/planar_jc48_sep06_collision.md) proves that in dimension three the exact blind quotient consists of quadrics through the tangent directions. [S166's physical interpolation banks](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/session166_planar_metagraph_20261103/collision_banks/NOTE.md) repair that quotient with six literal three-state guard directions.

The new lever is to retain **all three raw guard states underlying a fixed-port observer**, rather than the single three-state guard. [S164's guard algebra, sections 4–5](https://github.com/eliottcassidy2000/math/blob/2374ba81efc7ad585d723015630b546a3cf01ff1/05-knowledge/results/session164_emergent_20261103/profile_algebras/GUARD_ALGEBRA.md) already proves the exact nine-coordinate physical transfer. The present result supplies a sharp collision-detection bank for this full interface. Three dimensions do not merely become nine copies of the old test: the shared word acts on all three query blocks, and the full radius-three bank has a new three-dimensional blind quotient.

The obstacle remains physical source closure. Collision velocities below are arbitrary affine-space germs whose **tangents** are actual state vectors. They are not perturbations of the tournament or new allowed attachments. Completeness of their period test cannot manufacture ternary +2/+4 source returns.

## 2. The exact collision quotient in any dimension

Work over a characteristic-zero field. Let V have dimension n≥3 and let nonzero labelled t_1,…,t_m span V. Define

    N = direct_sum_i V / k t_i,
    R = {lambda : sum_i lambda_i t_i = 0},
    delta(w) = ([w],…,[w]),
    P([z_i])(lambda) = sum_i lambda_i t_i wedge z_i.

A branch system f_i(u,epsilon)=u t_i+epsilon z_i, epsilon²=0, preserves one common target to first order, up to branch reparametrization, precisely when [z_i] belongs to delta(V).

Define

    S_t = {A: V → wedge² V linear : t_i wedge A(t_i)=0 for every i},
    K_0 = {A_w(v)=v wedge w : w in V}.

**Lemma 1.** There is a natural isomorphism

    ker P / delta(V)  =  S_t / K_0.                    (1)

**Proof.** Vanishing of P says that t_i↦t_i wedge z_i respects every linear relation, so extends uniquely to a linear A on V. For a nonzero t, the kernel of wedge multiplication t:wedge²V→wedge³V is t wedge V. Hence the displayed condition on A is exactly what permits recovery of each z_i modulo k t_i. The common-motion image corresponds exactly to A_w. These constructions are inverse. ∎

This avoids assuming that the special three-dimensional scalar-quadric description survives unchanged in higher dimensions. The relevant quadratic response is now the vector-valued polynomial v↦v wedge A(v).

## 3. Sharp n+3-direction theorem

**Theorem 2.** Among banks of directions spanning an n-dimensional V, n≥3, exactly n+3 directions can suffice for ker P=delta(V), and fewer never suffice. An explicit complete bank is

    e_1,…,e_n,
    u=(1,…,1),
    a=(a_1,…,a_n),
    b=(a_1²,…,a_n²),                                (2)

where the a_i are distinct and nonzero.

The spanning hypothesis matters. Nonspanning banks can have other, smaller completeness behavior; this theorem makes no universal lower-bound claim for them.

**Lower bound.** Hom(V,wedge²V) has dimension n binom(n,2). Each t condition has codimension at most binom(n−1,2), because it is the condition A(t)∈t wedge V. The common-motion subspace has dimension n. For m≤n+2,

    dim(S_t/K_0)
       ≥ n binom(n,2) − (n+2) binom(n−1,2) − n
       = n−2 > 0.                                   (3)

**Sufficiency.** The coordinate-axis conditions force

    A(e_i)=e_i wedge w_i,
    w_i = sum_(j≠i) c_ij e_j.

For distinct i,j,k, the coefficient of e_i wedge e_j wedge e_k in t wedge A(t) is

    t_i t_j (c_jk−c_ik)
      + t_i t_k (c_ij−c_kj)
      + t_j t_k (c_ki−c_ji).                         (4)

For the three extra directions in (2), divide (4) by t_i t_j t_k. The resulting equations on the three parenthesized differences have coefficient rows

    (1,1,1),
    (a_k^(-1),a_j^(-1),a_i^(-1)),
    (a_k^(-2),a_j^(-2),a_i^(-2)).

Their Vandermonde determinant is nonzero. Every difference is therefore zero. For each fixed j, all c_ij with i≠j are equal, so there is a single w with A(e_i)=e_i wedge w for all i. Thus S_t=K_0. ∎

The theorem is elementary exterior algebra. It is a theorem about arbitrary spanning directions, not about their realization by a specified physical language.

## 4. The physical nine-coordinate carrier

Use temporal word order: a word is processed from left to right. A guarded tournament X has guard q with unique in-neighbor r and raw state

    v(X,q)=(s,u,d),
    s=# Hamilton paths starting at q,
    u=# Hamilton paths not starting at q,
    d=# Hamilton paths of X−q starting at r.

The inherited one-vertex moves are

    C(s,u,d)=(u,s+u+d,s),
    D(s,u,d)=(s+d,2u,0).

Each adds a with q→a and a→X−q. C marks a; D retains q. Both preserve strength. Their matrices are

    C=[[0,1,0],[1,1,1],[1,0,0]],
    D=[[1,0,1],[0,2,0],[0,0,0]].

Keep a distinct observer port p. A query adds a vertex beating only p and marks that new vertex as the next observer port; the guard remains q. Let X_0=X, X_1 and X_2 be the zero-, one- and two-query graphs. Retain the entire matrix

    V_X = [v(X_0,q) | v(X_1,q) | v(X_2,q)].           (5)

S164 proves that each C/D attachment commutes with these queries, so a word W acts by V_X↦M_W V_X. This is an exact nine-coordinate state, not the old three-coordinate observer itself.

Take the directed triangle 0→1→2→0, guard q=0 and port p=1. The query adjacency row masks are

    X_0 = (2,4,1),
    X_1 = (10,4,9,2),
    X_2 = (26,20,25,2,8).

Their raw state columns are (1,2,0), (2,3,1), (3,6,2), giving

    V_0=[[1,2,3],[2,3,6],[0,1,2]],     det V_0=−2.     (6)

Thus M↦vec(MV_0) is an invertible nine-dimensional linear map. Matrix-word coordinates and literal nine-query-state coordinates have exactly the same collision-completeness behavior. Row-major vectorization is used throughout.

## 5. A simultaneously cardinality- and radius-optimal physical bank

**Theorem 3, finite-exact.** The following twelve temporal words give a complete, spanning collision bank in the nine-dimensional state (5):

    CDDD, DD, CDCD, CCD, CDC, DDCC,
    DCD, DCDC, CDD, CCC, CCCC, C.                     (7)

Every word uses at most four moves. The base tournaments have at most seven vertices; the two queried graphs used to read a state have at most nine. All these graphs are literal and strong.

The bank has minimum possible cardinality among spanning nine-dimensional collision banks, by Theorem 2. It also has minimum possible maximum word length in this triangle-seeded C/D full-query language.

### Positive finite certificate

The first nine listed physical vectors are a basis with determinant 2048. In that basis, the coordinate-axis conditions leave 9·8=72 coefficients c_ij. The constraints supplied by the last three directions have rank 63 modulo the prime 1009. This is a sound characteristic-zero lower bound. Common motion supplies a nine-dimensional kernel over Q, giving the matching upper bound. Therefore the kernel is exactly common motion.

As an independent formulation, the checker directly constructs the original wedge-period map using a basis of the three tangent relations. Its raw-velocity matrix is 108×108 and has rank 87 modulo 1009. There are always 12 tangential directions and nine common-motion directions in its kernel; these are independent because the tangents span V. This again proves its exact characteristic-zero rank is 87 and that there is no further collision ambiguity.

### The radius-three obstruction

Use the following nine word matrices as a coordinate basis:

    empty, C, D, CC, CD, DC, DD, CCD, CDD.             (8)

Their row-major vectorizations have determinant 16. Let e_0,…,e_8 denote this basis, transported to physical coordinates by M↦vec(MV_0), and write x=(x_0,…,x_8). Every one of the fifteen words of length at most three satisfies x wedge A_j(x)=0 for all three maps

    A_1(x)=(x_4−x_5) e_4 wedge e_5,
    A_2(x)=x_7 e_7 wedge (−e_0−e_1−e_2+2e_3+e_4+e_5+e_6),
    A_3(x)=x_8 e_8 wedge ( e_0+e_1+e_2−2e_4−2e_5−e_6). (9)

These identities are checked over Q without modular reduction. None of the three is common motion, and their classes are independent: all vanish on e_0 and e_1, so a common-motion combination would have w=0; then their respective e_4, e_7 and e_8 responses force all coefficients zero.

The radius-three normalized constraint matrix has rank 60 modulo 1009. The nine common motions and the three exact maps (9) give a twelve-dimensional rational kernel. Consequently rank is exactly 60 and the blind quotient has exactly dimension three. No spanning subbank can be complete, since a complete subbank would force every map satisfying the full bank to be common motion. Radius two has span only seven. This proves the sharp radius four claim.

## 6. What this changes for the ternary research

1. A full fixed-port guard interface really is nine-dimensional, already from the smallest triangle seed. No abstract enlargement was needed to obtain spanning physical tangent directions.
2. The planar collision mechanism now yields an exact, optimal test bank for that interface. The radius-three failure has explicit invisible motions, not just a failed rank search.
3. A new candidate first-order response formula on the full interface can be tested for collision compatibility using (7), provided its polynomial/linear dependence and interpretation are proved separately.
4. The bank detects compatibility of arbitrary ambient velocities on the twelve selected branches; extending that conclusion to every physical state requires a separate theorem. It does **not** convert an algebraic common motion into a physical attachment, identify source histories, establish affine ternary offsets, or prove any Hamilton-spectrum or Jacobian conjecture.

The next genuine construction question is whether an offset gadget has a lawful response on this full carrier whose target state lies in a reusable physical source bank. The present theorem makes one necessary linear test exact; it does not replace that source-return obligation.

## 7. Reproduction

Run:

    python3 check_collision_bank.py > collision_results.json
    python3 -O check_collision_bank.py > /tmp/collision_results_optimized.json
    cmp collision_results.json /tmp/collision_results_optimized.json

The standard-library checker uses always-active failures. Its 186 checks include generic banks in dimensions 3–9, independent literal Hamilton-path subset DP on every selected/base word and its two query graphs, strength and unique-in-neighbor checks, the seed determinant, both physical rank certificates, the exact radius-three hostiles, their independence, and the radius-two span. Normal and optimized output agree byte-for-byte. Modulo-1009 ranks are used only as characteristic-zero lower bounds; every needed upper bound is supplied by explicit rational kernel vectors or common motion.
