# Quartic deletion tomography and an exterior edge-derivative theorem

**PROVED:** a complete named quartic normal form; a sharp five-old-vertex threshold for arbitrary quartic response contrasts; and an all-order exponential quartic for twin reversal over a transitive exterior, independent of every root incidence. **FINITE-EXACT:** exhaustive normal-form controls through order five, stated larger controls, and complete literal response tables for one strong hostile. **INCOMING PROVED:** commit `76b3b1a` closes every single nontrivial hidden module converse by a separate parity theorem, including the former strong prime-core residual. **OPEN:** exact full-response collisions between arbitrary unrelated tournaments and universal named reconstruction.

This extends the S186 twin-contact and S187 active-SCC reductions. It does not infer equality of full responses from selected marginals, or discard quartics for a partial module reversal. No literature-priority claim is made.

## 1. Observer and cover data

All names are physical. Let F_T(a,c) count Hamilton paths after adding x->y, where a_v=1 means v->x and c_v=1 means v->y. Its unique Boolean multilinear polynomial has bidegree at most (2,2), because each port has at most two path incidences. A quartic coordinate therefore has exactly two distinct a-names and two distinct c-names. A physical name may occur in both ports.

For distinct j,k in a tournament R, let C_jk(R) count ordered covers (P1,P2) by two nonempty directed paths, with P1 ending at j and P2 starting at k. The paths cover R, and the boundary arc j->k is **not** tested.

For four distinct names u,v,w,z, let D_(uv;wz)(T) count ordered covers (P1,P2,P3) by three nonempty directed paths, with P1 ending at u, P2 starting at v and ending at w, and P3 starting at z. Neither boundary arc u->v nor w->z is tested. Thus the middle path contains at least two names. These are specified covers, not products of endpoint marginals without disjoint support.

## 2. Complete quartic normal form

**Theorem 1.** Every named quartic of F_T is given by exactly one of the following three cases.

1. Two shared names:

       [a_i a_j c_i c_j] F_T = 0.

2. Exactly one shared name, with i,j,k distinct:

       [a_i a_j c_i c_k] F_T
           = C_jk(T-i) + C_kj(T-i).                         (1)

3. Disjoint pairs A={i,j} and B={k,l}:

       [a_i a_j c_k c_l] F_T
           = sum over (u,v)=(i,j),(j,i)
                     and (w,z)=(k,l),(l,k)
               (D_(uv;wz)(T) + D_(wz;uv)(T)).               (2)

In particular every quartic coefficient is nonnegative. The complete top-degree part is symmetric under exchanging the ports, and is unchanged by replacing T by its whole converse. Neither assertion applies to arbitrary lower-degree coordinates.

**Proof.** A term contributing degree four must have both ports internal to its Hamilton word, with two old neighbors each. Every incoming port factor is a named variable; every outgoing port factor is one minus a named variable. Taking the degree-four part selects the negative term of each of the two outgoing factors, giving positive sign. A word using x->y cannot contribute degree four because that port-port arc leaves at most two old-port incidences.

If both ports have the same two old neighbors i,j, those four vertices would form a closed four-cycle in the undirected path-incidence graph. This is impossible in a Hamilton path. This proves case 1.

In case 2 the only possible consecutive blocks are

       j,x,i,y,k       and       k,y,i,x,j.

Removing i,x,y from the first leaves exactly a two-path cover of T-i counted by C_jk; the second leaves one counted by C_kj. Inserting the block is the inverse map, including singleton cover paths. Their degree-four signs are positive. This proves (1). In particular **all old edges incident to i disappear from this coordinate**.

In case 3 the port neighborhoods are two disjoint consecutive blocks. Either block can come first, and each has two orientations. Removing the ports gives precisely the eight three-cover classes in (2), with the two gap arcs untested. Conversely each cover uniquely reinserts x and y. This proves the formula.

Exchanging ports just exchanges the two unordered marked pairs and preserves these counts. Reversing every old path exchanges C_jk with C_kj and sends each three-cover class to another one among the same eight. Hence the entire quartic cancels under whole converse. QED.

For n old names, the full bidegree-(2,2) quartic basis still has binomial(n,2)^2 named coordinates, including zeros. The normal form has n*binomial(n-1,2) shared-name capacity entries and 3*binomial(n,4) disjoint-pair entries; the remaining coordinates are recovered by the exact named symmetries above. This is a fixed-coordinate identity, not coefficient sorting or an unlabeled quotient.

## 3. Consequence for the hidden-module question

**Exact deletion-capacity obstruction.** Full named response equality F_T=F_U implies

       C_jk(T-i)+C_kj(T-i)
         = C_jk(U-i)+C_kj(U-i)

for every three distinct physical names i,j,k, as well as equality of every eight-term three-cover sum in (2).

Apply this with T=Q[M] and U=Q[M^op], with the exterior and all names fixed. These are necessary conditions on any hidden module converse in **any** quotient. They keep precisely the degree-four information omitted by assuming a partial reversal is a whole converse. They do not establish exact descent to M, and their logical independence from the complete degree-at-most-three response data is not asserted. Existing residue tests, twin-sign separation, and S187 active-SCC localization remain applicable. The fresh incoming E1-E4 parity theorem at commit `76b3b1a` now excludes every single nontrivial module converse, including the previously open strong prime residual. Our quartic conditions remain valid, but are no longer an open-route necessity for that particular exclusion. See `INTEGRATION.md`.

## 4. Sharp universal order threshold

**Theorem 2.** For every pair of tournaments T,U on the same n<=4 old names, F_T-F_U has total degree at most three. Five is the smallest old order at which an arbitrary named-response difference can have degree four, even when both tournaments are strong.

**Proof of the lower bound.** There are no quartics for n<=2. For n=3, any allowed shared-name coefficient deletes to a two-vertex tournament; each of C_jk and C_kj equals 1, so the coefficient is 2. For n=4, a shared-name coefficient deletes to the three-vertex set {j,k,l}. Its capacity is

       1[l->j]+1[k->l]+1[l->k]+1[j->l] = 2.

For a disjoint coefficient at n=4, a three-cover has two singleton end paths and a middle path on two vertices. Across the eight orientation/order choices in (2), each direction of each of the four cross-pair edges occurs exactly once. Opposite directions pair up, so exactly four choices are valid. Thus every disjoint coefficient is 4. All quartics are therefore independent of the old orientation at these orders, and cancel in every comparison. The five-vertex strong witness below proves sharpness. QED.

This is stronger in scope than the S186 finite minimum-order census restricted to homogeneous-pair reversals. It does not extend the whole-converse cubic theorem to higher-order arbitrary comparisons.

## 5. Twin reversal sees every exterior edge derivative

Let T contain homogeneous pair u->v, and let T' reverse only that edge. Write Q for its contraction to z and R=Q-z for the exterior tournament. For distinct exterior names i,j:

**Theorem 3.**

       [a_v a_i c_u c_j](F_T-F_T')
           = C_ji(R)-C_ij(R).                              (3)

The right side is independent of every incidence between z and R.

**Proof.** Terms using neither internal orientation cancel. A surviving word in T contains u,v consecutively; in T' it contains v,u. For the specified quartic, x has old neighbors v,i and y has neighbors u,j. In T the only relevant block is

       j,y,u,v,x,i,

which deletes to a two-cover of R counted by C_ji. In T' the block is

       i,x,v,u,y,j,

which deletes to a cover counted by C_ij. Each quartic sign is positive. Any alternate orientation of either port block would make the required u-v adjacency impossible. Thus the difference is (3). The block contains both twins and both ports, so no root-to-exterior arc is tested at its two deleted boundaries. QED.

For clarity, C_ij(R) and C_ji(R) are independent of the existing edge between i and j, because they place its ends in separate paths. If R^+ forces i->j and R^- forces j->i while keeping every other edge, then

       H(R^+)-H(R^-)=C_ij(R)-C_ji(R).

Indeed paths avoiding that edge cancel, and paths using it contract to the indicated two-covers. Hence (3) is exactly the **negative exterior edge-flip derivative**. Any exterior edge whose reversal changes H certifies a genuine degree-four twin response. Vanishing of all these derivatives is only necessary for this family of quartics to vanish; it is not claimed sufficient for the entire quartic to vanish.

## 6. An all-order exponential quartic, with strong realizations

Let the m>=3 exterior names form the transitive order b1->...->bm, with arbitrary common incidences to u,v. Then

       C_b1,bm(R)=0,           C_bm,b1(R)=2^(m-2).

For the first count, a path ending at the first name must be its singleton and a path starting at the last name must also be its singleton, leaving an uncovered name when m>=3. For the second, b_m belongs to the first path and b_1 to the second. Each of the other m-2 names may be assigned independently to either path, and its transitive path order is unique.

Therefore, with n=m+2 old names,

       [a_v a_b1 c_u c_bm](F_T-F_T')
           = 2^(m-2) = 2^(n-4) > 0.                       (4)

This holds for every one of the 2^m possible common root-incidence patterns, not just strong quotients. For strong examples, make the root beat b1 and make every other exterior name beat the root. The quotient contains the Hamilton cycle z,b1,...,bm,z. Replacing z by either orientation of a nonempty pair preserves strongness, as paths can leave and return through the quotient. Thus (4) supplies strong examples at every old order n>=5.

The first example, on names 0,...,4 with u=0, v=4 and b1,b2,b3=1,2,3, is

       Q  = (2,12,9,1),
       T  = (18,12,25,17,2),
       T' = (2,12,25,17,3).

Both are strong, and [a_4 a_1 c_0 c_3](F_T-F_T')=2. The complete partial-reversal difference, the complete whole-converse difference, and all sixteen literal finite-difference terms are retained. This is a new convenient family and labeling, not a claim that a five-vertex quartic hostile was previously unknown. The lower-degree inherited twin sign also separates these tournaments; no necessity of quartic data for first detection is claimed.

## 7. Reproduction and evidence boundary

Run the standalone standard-library checker:

    python verify.py > receipt.json
    python -O verify.py > receipt.optimized.json

Every check raises an exception on failure, including optimized execution. The checker compares the full polynomial obtained by literal Hamilton-word expansion against the separate deletion endpoint DP and old three-cover enumeration. It enumerates every labeled tournament on 1,...,5 old names: 1,099 cases, retaining each adjacency and a digest of all named quartic coordinates, including zeros. It checks whole-converse quartic cancellation throughout that census.

For the exterior derivative identity it checks every labeled fixed-root quotient through order four, plus 20 explicitly retained larger quotients. It verifies literal sixteen-query certificates for every root-incidence pattern at m=3 and the strong family at m=4,...,7. Finally it evaluates all 1,024 actual attachment queries for the five-vertex hostile pair by independent Hamilton subset DP and compares the complete difference polynomial, not merely one coefficient. Further literal high-weight controls check full responses across the small census.

The unbounded results are the bijections and transitive-cover count above. The finite cases do not prove universal reconstruction or existence of a full-response collision. The all-order single-module exclusion is inherited from the separate incoming parity theorem, not inferred from these finite controls.
