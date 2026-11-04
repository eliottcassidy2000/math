# Named rank reconstruction and persistent quartic clocks in the three-state gate

**PROVED:** explicit all-order singleton-port response clocks, reconstruction of the six seeded orientations with all named transitive-module orders, one fixed all-order whole-converse query, and a positive all-order pair-reversal quartic coefficient. **FINITE-EXACT:** complete named response-difference polynomials and their actual interpolation counts for all three states at r=1,...,7. This is a restricted-family reconstruction theorem, not universal tournament reconstruction or exclusion of the prime twelve-core residual.

## 1. Physical objects, source mechanisms, and what is new

Use the S185 seeds, outgoing masks on fixed names 0,...,5:

    Q0=(16,5,33,39,46,3),
    Q1=(16,5,49,39,42,3),
    Q2=(18,4,49,39,42,3).

Replace name 5 by the transitive module on names 5,...,r+4, with that increasing order, obtaining T_j(r). Every seed and inflation is strong. All exterior incidences retain the physical root type: the module beats exactly exterior names 0,1. The inherited scalar law is

    H_j(r)=4*3^r+9*2^r-3+2j.

The two gate changes are the physical exterior edges {2,4} and {0,1}. A port query adjoins x->y; a_i=1 means i->x, c_i=1 means i->y. F_T counts Hamilton paths in this actual tournament. All comparisons retain names.

Closest mechanisms: S185's occupied-gap law; S186 G1's twin contrast; the incoming homogeneous-pair packet's calibrated singleton-query contrast; S186 structure's deletion-cover quartic identity; and T4's whole-converse cubic theorem. The previous named twin separator is credited, not claimed anew. The new objects here are the complete rank clock, its restricted reconstruction consequence, and the all-order nonzero quartic clock. The scalar-context hostile and prime twelve-core residual remain untouched.

Write E_v for the number of Hamilton paths ending at v, and C_vw for the number of ordered spanning two-path covers with first path ending at v and second starting at w, including directions absent as actual edges. Put

    S_v(T)=E_v(T)+sum_w C_vw(T).

This counts permutations in which all consecutive old arcs except the outgoing adjacency of v must be directed correctly. If v is last, there is no ignored adjacency.

## 2. A general named ranking theorem for transitive modules

**Theorem 1.** For every tournament T and every named v,

    F_T(e_v,0)=H(T)+2S_v(T).                                      (1)

If a->b is a homogeneous pair and Q contracts it to z, then

    S_b(T)-S_a(T)=S_z(Q)>=H(Q)>0.                                (2)

Consequently, if a prescribed vertex set is a transitive module in T, sorting its singleton-port values F_T(e_v,0) in increasing order recovers its entire directed order, in any fixed exterior tournament.

**Proof of (1).** Here y beats every old vertex, x beats every old vertex except v, and x->y. If x is first, y must immediately follow it; this contributes H(T). Otherwise x must follow v. Either y is first, or it immediately follows x. Each choice is a bijection with an old permutation counted by S_v: insert x after v and place y first, or insert the block x,y after v. The cases are disjoint and exhaustive. This includes v at the end. Hence (1).

**Proof of (2).** Compare the wildcard-outgoing permutation interpretations of S_b and S_a by swapping names a,b. Words where a,b are not adjacent cancel in pairs, since their exterior incidences agree. Of the adjacent cases, words a,b counted in S_a correspond to words b,a counted in S_b, with the internal arc ignored in each. The remaining words are a,b counted in S_b; words b,a counted in S_a are impossible. Contract the remaining block a,b to z. Its outgoing adjacency remains ignored, giving exactly S_z(Q). Every actual Hamilton path of Q is among these words, and every tournament has a Hamilton path. This proves (2). Consecutive vertices in a transitive module form homogeneous pairs; applying (2) successively gives strict increase along the module. QED.

The theorem assumes a prescribed transitive module. It is not a test recognizing arbitrary modules from a score ranking, and it makes no claim that singleton probes reconstruct arbitrary tournaments. It strengthens the usable interface from one pair's sign to the total order on an arbitrarily large named module.

## 3. Closed response clocks for every gate module position

For k=1,...,r, let v_k=4+k, the kth module vertex. Define

    (B_0,C_0,D_0)=(44,-5,-4),
    (B_1,C_1,D_1)=(42,-1,-3),
    (B_2,C_2,D_2)=(42, 3,-3).

**Theorem 2.** For every j=0,1,2 and r>=1,

    S_(v_k)(T_j(r))
      =18*3^(r-k)*4^(k-1)
       +B_j*2^(r-k)*3^(k-1)
       +C_j*2^(k-1)+D_j*1[k=r].                              (3)

Together with (1), this gives each actual named singleton-port count, not just a difference or an endpoint marginal. The response clocks include a genuine 4^(k-1) term that is absent from the scalar H clock.

**Finite skeleton certificate for the boundary.** Keep the first module vertex v_1, and replace the other r-1 module vertices by runs inserted into permutations of the five exterior names and v_1. Ignore only the outgoing adjacency of v_1. The numbers of valid skeletons with k occupied module gaps, k=0,1,2,3, are

    state 0: (53,137,98,18),
    state 1: (56,137,96,18),
    state 2: (60,141,96,18).                                (4)

Here a gap may be occupied iff its preceding name beats the remaining module (or is absent) and the remaining module beats its following name (or is absent). The ignored outgoing adjacency is relevant only for an unoccupied gap after v_1. Every backward tested exterior adjacency must occupy its gap. Thus for a skeleton order p, if B is the mandatory-gap set and G the permissible-gap set, its coefficient is binom(|G|-|B|,k-|B|) when B is a subset of G, and zero otherwise. This is a complete 720-permutation integer certificate per state, reproduced in verify.py. There are at most three module runs, because each nonterminal run must precede one of the distinct names 0,1. The separate ignored adjacency at v_1 does not change that fact.

Write q_k(m)=k!S(m,k), with q_0(0)=1, q_0(m)=0 for m>0. Ordered nonempty run assignments of the remaining linearly ordered vertices contribute q_k(r-1). Therefore (4) proves, for all r,

    f_j(r):=S_(v_1)(T_j(r))
       =18*3^(r-1)+B_j*2^(r-1)+C_j+D_j*1[r=1].            (5)

For successive positions, contraction in (2) gives

    S_j(r,k+1)-S_j(r,k)=S_j(r-1,k).                         (6)

The initial values (5) and recurrence (6) uniquely give

    S_j(r,k)=sum_(i=0)^(k-1) binom(k-1,i) f_j(r-i).

Applying the binomial theorem to each exponential proves (3), including the endpoint correction because r-i=1 occurs in this sum only when k=r and i=r-1. QED.

For example state 0 at r=4 has module scores (833,1166,1636,2296). Its actual singleton-port counts are H_0(4)+twice these numbers. Each adjacent difference is twice the corresponding contracted score, not an extrapolated numerical fit.

## 4. One fixed global-converse probe and restricted reconstruction

Set A_j(r)=S_0(T_j(r)) and Bop_j(r)=S_0(T_j(r)^op). Direct occupied-gap certificates at the exterior marked name 0 give the following coefficients of (q_0(r),q_1(r),q_2(r),q_3(r)):

    j     original              whole converse
    0     (15,45,39,8)          (9,27,21,4)
    1     (27,75,63,14)         (12,30,21,4)
    2     (18,54,51,14)         (18,39,24,4).

These enumerate the 120 permutations of the other five seed names, with outgoing adjacency of 0 ignored. The same mandatory/permissible gap formula proves the all-order result. Since r>=1, the q_0 term contributes zero. Expanding yields

    A_0=8*3^r+15*2^r-9,     Bop_0=4*3^r+9*2^r-3,
    A_1=14*3^r+21*2^r-9,    Bop_1=4*3^r+9*2^r,
    A_2=14*3^r+9*2^r-6,     Bop_2=4*3^r+12*2^r+3.

**Theorem 3.** The same physical query (a,c)=(e_0,0) separates every T_j(r) from its whole converse, with positive difference

    j=0: 2*(4*3^r+6*2^r-6),
    j=1: 2*(10*3^r+12*2^r-9),
    j=2: 2*(10*3^r-3*2^r-9).                              (7)

Positivity is immediate for j=0,1; for j=2, 10*3^r-3*2^r-9>=7*3^r-9>0 for r>=1. Formula (1) and converse invariance of H prove (7). This covers r=1,2 as well as the larger sizes handled by the old three-module-name contrast.

**Restricted reconstruction corollary.** Fix the five exterior names and the module's named vertex set. Consider all tournaments obtained by choosing j, choosing either the original exterior orientation or its whole converse, and choosing any of the r! transitive orders on the module. This family has 6r! distinct named tournaments. The following r+2 actual queries reconstruct the entire tournament:

1. F(0,0)=H determines j from H_j(r), with r determined by the given order n=r+5.
2. F(e_0,0) chooses between the two exterior orientations: its values are H_j+2A_j and H_j+2Bop_j, distinct by (7).
3. The r values F(e_v,0) for module names v recover their transitive order by Theorem 1.

Changing the named transitive order leaves exterior S_0 unchanged, since it merely relabels vertices inside a uniform module. In the whole-converse family the actual module order is still recovered by increasing singleton score, irrespective of the reversed exterior. Thus no unproved calibration for the opposite exterior is needed. All members are strong, since both seed orientations are strong. This is a strict finite-description family at every r, not a universal reconstruction claim or an optimal-query lower bound.

## 5. A positive quartic coefficient at every gate size

Let T'_j(r) differ from T_j(r) only by reversal of the consecutive pair edge 5->6, for r>=2. This is a physical homogeneous-pair reversal, not whole converse.

**Theorem 4.** For every j=0,1,2 and every r>=2,

    [a_0 a_3 c_0 c_6](F_(T_j(r))-F_(T'_j(r)))
       =2^(r-2)+2-1[r=2] >0.                              (8)

Consequently this entire internal-reversal response difference has degree exactly four at every size. Both tournaments are strong; T' simply exchanges the first two ranks of the transitive module.

**Proof.** Contract {5,6} to z in Q. The inherited exact contact formula gives

    [a_i a_j c_i c_v](F_T-F_T')=C_zj(Q-i)-C_jz(Q-i)

for twins u->v and distinct exterior i,j. Take u=5,v=6,i=0,j=3. The quotient module has r-1 vertices and z is its first member. Delete exterior 0. In C_z3, joining the two path words gives a word containing the consecutive block z,3, whose intervening arc is ignored. For C_3z use block 3,z instead. Keep these blocks adjacent while inserting the remaining r-2 module vertices into all other permissible gaps. The exact run coefficients are

    state 0: C_z3 profile (3,6,2), C_3z profile (1,2,1),
    state 1: C_z3 profile (6,11,4), C_3z profile (4,7,3),
    state 2: C_z3 profile (6,11,4), C_3z profile (4,7,3).

The common difference is (2,4,1), indexed by k=0,1,2. These are exact finite skeleton counts from permutations of the four remaining exterior names and z, subject to the frozen named adjacency. The gap between the frozen names is forbidden; the boundary arc itself is untested. The general gap formula reproduces every entry. Hence the coefficient is

    2*q_0(r-2)+4*q_1(r-2)+q_2(r-2)
       =2^(r-2)+2-1[r=2].

Every full response has bidegree at most (2,2), since a Hamilton word has at most two incidences at either port. Its total degree is at most four. The positive displayed coefficient proves equality. QED.

This is stronger than the old finite r=2 hostile. It does not say low-degree probes fail to distinguish the pair: Theorem 1 already distinguishes it. It says they cannot reconstruct its complete difference by pretending it is cubic. Whole converse remains covered by T4's degree-at-most-three theorem, and (7) makes that cubic-or-lower difference nonzero at every size.

## 6. Complete named data and independent checks

`compile.cpp` evaluates literal Hamilton subset DP on the actual port extensions. For each state and r=1,...,7 it retains every physical count on the weight-at-most-two-per-port battery for whole converse, and additionally pair reversal for r>=2. `full_response.json` contains 39 full comparison packets, including physical rows, the entire battery of paired actual counts, every nonzero named multilinear coefficient, degree, and four high-weight controls. No named coordinate is sorted or replaced by a marginal. Entries not retained as nonzero coefficients are zero by interpolation, not dropped unknowns.

The square battery has (1+n+binom(n,2))^2 entries. Boolean inversion gives the complete response-difference polynomial because bidegree is at most (2,2). Global-converse packets independently show the predicted zero quartic and pure-quadratic terms, so the same data imply the smaller T4 battery of 1+2n+n^3 suffices. Pair packets retain the quartic part, including (8). Complete named finite data are provided; the all-order claims rely on the proofs and finite skeleton certificates above, not extrapolation from these packets.

`verify.py` independently computes every S coordinate from endpoint and cover DP for all states and r=1,...,7; checks all clock formulas and recurrence; enumerates all finite gap certificates; checks the quartic by its literal 16-query finite difference at every r; inverts every retained battery again; and independently re-evaluates one dense physical query pair in every packet in Python. Normal and `python -O` replays use explicit exception gates and must agree byte for byte. A compiler is already available; no installation is needed.

Reproduction:

    g++ -O3 -std=c++17 compile.cpp -o /tmp/s187_gate_compile
    /tmp/s187_gate_compile > full_response.json
    python verify.py > receipt.json
    python -O verify.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

`base.py` and `literal.py` retain the separately written S186 endpoint/cover and literal-extension implementations. Exploratory files are not proof dependencies. No theorem about arbitrary prime twelve-cores, universal full-response injectivity, an exact response collision, an unrestricted spectrum, or a reusable ternary gate is claimed.
