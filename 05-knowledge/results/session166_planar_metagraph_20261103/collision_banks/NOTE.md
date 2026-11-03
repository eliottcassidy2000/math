# Sharp physical interpolation banks from collision geometry

**Status: PROVED, with exact finite controls.** A bridge from the planar-Jacobian collision research identifies the correct quadratic test for the three-state tournament guard language. It leads to interpolation banks that are simultaneously optimal in cardinality and maximum word length. The bridge constructs explicit map germs from physical tournament states; it does not identify tournament operations with Keller maps or claim any progress on polynomial invertibility.

## 1. Inherited mechanisms and the connection

The source collision theorem is [the earlier quadratic-sidecar theorem](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/planar_jc48_sep06_collision.md), following [THM-4411, first-order collision transgression](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/01-canon/theorems/THM-4411-first-order-collision-transgression-seminormal-tradeoff.md). In dimension two, scalar tangent-relation periods detect common first-order motion. In dimension three, even retaining all target two-forms loses exactly the quadrics through the tangent directions. Arbitrarily many directions on one conic remain blind.

The tournament input is [S165 guard-state geometry](https://github.com/eliottcassidy2000/math/blob/6b24421039fd2ffa820e9b56d107e458d75498ed/05-knowledge/results/session165_geometry_representation_20261103/guard_geometry/NOTE.md). A tournament X with marked vertex q having a unique in-neighbor r has state

    v=(s,u,d),
    s = number of Hamilton paths of X starting at q,
    u = number of Hamilton paths of X not starting at q,
    d = number of Hamilton paths of X−q starting at r.

Start from a directed triangle v0=(1,2,0). Add a vertex a with q→a and a→X−q. Operation C moves the guard to a; D keeps q. Both preserve strongness. Their column matrices are

    C = [[0,1,0],[1,1,1],[1,0,0]],
    D = [[1,0,1],[0,2,0],[0,0,0]].

**Word convention:** words below are temporal. Thus D^j C^k gives state C^k D^j v0. Length is the number of added vertices.

The connection map is explicit: send each physical guarded tournament to its vector v, then use v as the tangent of the branch germ f(u)=u v in affine three-space. This preserves the guard-count vector and literal realizability of the indexing tournament. It forgets its individual paths and graph isomorphism type. Branch deformations are arbitrary affine-space germs, not physical deformations of tournaments. The collision theorem then asks for a quadratic-unisolvent bank of physical vectors. Reset planes supply precisely such banks, and their mechanism extends to every degree.

## 2. Physical reset planes

Put L={d=0}, L_k=C^k L, and let

    ell_k(v) = (third coordinate of C^(−k) v).

Since det C=1, ell_k is an integral linear form. The planes L_k are distinct. Indeed, if k>j and a,b>0 then C^(k−j)(a,b,0) has strictly positive third coordinate. In particular every positive point C^k(a,b,0) lies outside L_j for j<k. This gives both plane separation and the nonvanishing factors needed by interpolation.

Two literal point families on L are

    D^j v0 = (1,2^(j+1),0),                              (1)
    D^n C D^m v0 = A_m(1,2^n,0),
    A_m=1+2^(m+1),  m≥0, n≥1.                           (2)

The first varies along the affine line s=1. The second varies along distinct radial lines u=2^n s, with distinct positive scales A_m≠1. These are actual positive words, not linear combinations of tournament states.

## 3. A minimal bank for homogeneous degree delta

For an integer delta≥0, use

    H_delta = {D^j C^k : j,k≥0, j+k≤delta}.             (3)

**Theorem A.** Evaluation on the physical states in H_delta is an isomorphism from the space of homogeneous polynomials of degree delta in three variables over Q to Q^(binom(delta+2,2)). Every test uses at most delta moves, hence at most delta+3 vertices. Both the number of states and the maximum word length are optimal within this C/D language.

**Proof of sufficiency.** Let a homogeneous P of degree delta vanish on the bank. On L_0, it has delta+1 distinct projective zeros (1,2^(j+1),0), so P|L_0=0 and ell_0 divides P. Suppose ell_0...ell_(k−1) has been removed. At the delta−k+1 points on L_k, these previous factors are nonzero by positivity. Pullback by C^k therefore gives a homogeneous binary form of degree delta−k with that many distinct slopes. It vanishes identically and ell_k divides the residual. Continuing gives delta+1 distinct linear factors in degree delta, so P=0. Dimension equality then gives an evaluation isomorphism and in particular proves that no two bank states coincide. ∎

**Optimality.** A vector space of dimension binom(delta+2,2) cannot be tested injectively by fewer scalar evaluations. For radius, every word of length at most R lies in one of L_0,...,L_R: take the last D, or the seed if there is no D, and then its final C-run. Consequently the nonzero homogeneous polynomial ell_0...ell_R vanishes on every such state. When delta≥1, all words of length at most delta−1 therefore fail to detect a degree-delta polynomial. The bank (3) attains the required delta moves. For delta=0 the single seed suffices.

The same bank detects a homogeneous polynomial whose degree is known only to be at most delta: multiplying a hypothetical vanishing polynomial of degree e by s^(delta−e) reduces to the stated theorem.

## 4. A minimal bank for all polynomials of bounded total degree

Augment (3) by

    B_delta = H_delta union
      {D^m C D^n C^k : m,k≥0, n≥1, m+n+k≤delta}.        (4)

**Theorem B.** For delta≥0, B_delta is an interpolation bank for all polynomials of total degree at most delta over Q. Its cardinality is exactly binom(delta+3,3), the dimension lower bound. For delta≥1 its maximum word length is delta+1, which is also optimal within the C/D language. The corresponding tournaments have at most delta+4 vertices. For delta=0, the bank is just the triangle.

### A two-dimensional triangular lemma

For r≥0, consider the following points (a,b) in an affine plane:

    (1,2^(j+1)),                       0≤j≤r;
    A_m(1,2^n),                       1≤n≤r, 0≤m≤r−n.

There are (r+1)+r(r+1)/2=binom(r+2,2) points, and they interpolate all total-degree-at-most-r polynomials.

To prove this, a vanishing polynomial restricts to zero at r+1 points on a=1, hence has factor a−1. The residual has degree at most r−1. On b=2a there are r distinct nonzero scales A_m; they avoid a=1. Its restriction vanishes identically, so b−2a is another factor. On b=4a there are r−1 scales, avoiding the previous lines, and so on. The resulting product of r+1 distinct affine-line factors exceeds the original degree. This proves the lemma, including r=0.

### Proof of Theorem B

On each L_k the part of (4) uses the two-dimensional lemma with r=delta−k, pulled forward by C^k. If P vanishes on B_delta, the k=0 lemma gives ell_0|P. At k=1 all old factors are nonzero on the test nodes, so the lemma applies to the residual. Continue as in Theorem A. This forces delta+1 linear plane factors and gives P=0.

The count is

    sum_(k=0)^delta binom(delta−k+2,2)=binom(delta+3,3).

Words in H_delta cost at most delta. Additional words cost m+1+n+k≤delta+1. The proof of injectivity also verifies distinctness, so the count does not accidentally count duplicated physical states.

For sharp radius, delta≥1, define the affine form

    a_delta(v)=(1,0,−1) C^(−(delta−1)) v − 1

and the nonzero polynomial

    Q_delta(v) = a_delta(v) product_(j=0)^(delta−2) ell_j(v),  (5)

with the empty product interpreted as 1. Its total degree is delta. Any word of length at most delta with a final C-run of length at most delta−2 lands on one of the displayed planes. The remaining possible states are

    C^(delta−1)v0, C^(delta−1)Dv0, C^delta v0.

After applying C^(−(delta−1)), these become v0,Dv0,Cv0, all satisfying s−d=1. Thus a_delta kills all of them. Formula (5) vanishes on the entire radius-delta word ball, not just on a selected bank. No bank of that radius can interpolate degree-delta polynomials, and (4) attains delta+1. ∎

**Improvement over S165.** The prior homogeneous bank used (delta+1)^2 tests and up to 2delta moves; the prior general bank used (delta+1)^3 tests and up to 3delta+2 moves. Here both test counts are the exact vector-space dimensions, and the radii delta and delta+1 are proved optimal for this literal language. This strengthens the physical polynomial-identity test, not physical state reachability or odd Hamilton-count coverage.

## 5. The sharp six-marked-state collision test

At delta=2, Theorem A gives these six words and directions:

    empty : (1,2,0)
    D     : (1,4,0)
    DD    : (1,8,0)
    C     : (2,3,1)
    DC    : (4,5,1)
    CC    : (3,6,2).

All are literal strong guarded tournaments with at most five vertices. These are six distinct marked states; C and D at order four have isomorphic unmarked tournaments. In the quadratic basis (s²,u²,d²,su,sd,ud), their evaluation matrix is

    [[1, 4,0, 2,0, 0],
     [1,16,0, 4,0, 0],
     [1,64,0, 8,0, 0],
     [4, 9,1, 6,2, 3],
     [16,25,1,20,4,5],
     [9,36,4,18,6,12]],                  determinant = 192. (6)

For these tangents t_i, take arbitrary velocities n_i in Q³ and form branch germs

    f_i(u,epsilon)=u t_i+epsilon n_i, epsilon²=0.

They admit one common target motion after independent branch reparametrizations if and only if

    sum_i lambda_i (t_i wedge n_i)=0                 (7)

for the following three rows lambda:

    (2,−3,1,0,0,0),
    (−3,1,0,−1,1,0),
    (2,−1,0,−2,0,1).

Each row is an actual relation sum lambda_i t_i=0, and they form the full relation space. Thus only nine scalar wedge-coordinate tests are required; their linear map on eighteen velocity coordinates has rank nine. Its kernel consists exactly of common motions plus six independent tangential reparametrizations.

### Why the collision theorem applies

For spanning nonzero t_i in V=Q³, write N=direct_sum_i V/Qt_i and let delta:V→N be common motion. The condition (7) means the assignment t_i↦t_i wedge n_i extends to a linear A:V→wedge²V. A chosen volume form turns A into the bilinear form b(u,v)=vol(u wedge A(v)). Its diagonal q(v)=b(v,v) vanishes at every t_i. Modulo common motion, precisely the symmetric quadratic part remains. Hence

    ker(periods)/delta(V) is isomorphic to I_2(t),

the quadrics vanishing on the t_i. Conversely every such quadric gives A by polarization, and A(t_i) belongs to t_i wedge V because q(t_i)=0, so it arises from some n_i. This proves the stated quotient rather than assuming that a rank test alone means common motion. Equation (6) makes I_2(t)=0.

Six directions are necessary for this kind of complete three-dimensional scalar-period test: fewer evaluations cannot kill all six quadratic coefficients. The bank also has the shortest possible maximum C/D word length by Theorem A.

### A retained hostile

Drop the final CC direction. The other five lie on the reducible conic

    q(s,u,d)=d(s−u+d)=0.

Thus their scalar-relation periods are not complete. An explicit invisible non-common motion uses

    n_1=(0,−1/2,0), n_2=(0,−3/2,0), n_3=(0,−7/2,0),
    n_4=(0,1/4,1/4), n_5=(0,1/8,1/8).

Indeed t_i cross n_i=M t_i with

    M=[[0,0,1/2],[0,0,−1/2],[1/2,−1/2,1]].

All tangent relations therefore kill the wedge data. The first three branches force a putative common velocity to (1/2,1/2,0), which is incompatible with the fourth. Adding the sixth direction detects the missing quadric, since q(3,6,2)=−2.

## 6. Boundaries and reproduction

The general theorems are characteristic-zero statements. The small determinant 192 already warns that reducing the six-point certificate modulo 3 changes the conclusion. The full algebra generated by C and D at odd primes is a different question from this chosen interpolation bank. No modular universality follows here.

All directions and test states are physically realized by the specified literal guard operations. Polynomial expressions evaluated on them are still required to be independently established as polynomial response formulas. The branch germs in Section 5 are actual affine-space map germs, but their deformations are not asserted to preserve a tournament, descend to a planar source, terminate polynomially, or supply a Keller pair.

Run `python3 check.py` and `python3 -O check.py`. The standalone standard-library checker verifies exact ranks/determinants, literal tournament counts, bank cardinalities and size bounds, finite polynomial-rank controls, radius-lower-bound witnesses, and both the complete and five-direction hostile collision systems. Always-active checks and deterministic JSON output make optimized replay meaningful. The interpolation and optimality proofs establish all degrees; finite controls do not replace those proofs. No claim of literature priority is made.
