# Transitive modules lift an old-arc repair to C1-to-C3 replacement families

**RESERVED / proof and finite coefficient certificate under independent audit.**

This supplies explicit unbounded families of physical rebuilding operations. It does not prove an eventual spectrum threshold or a recurrence for arbitrary represented integers. The resulting counts may have other known realizations; no novelty of the scalar values is asserted.

## Convention and operation

H is the number of directed Hamiltonian paths. C is the number of directed Hamiltonian cycles modulo rotation. Every row below is an outgoing bitmask on vertices numbered from zero.

Given a core Q and a mask M on its vertices, Q[M;k] adds a transitive module x 1 to ... to xk. Every xi beats exactly the core vertices selected by M and loses to the others. The module vertices have identical incidences to the core. Empty modules (k=0) are allowed.

## The two C1-to-C3 families

The seven-vertex parent cores, eight-vertex intermediate child cores, and masks are:

- B: P=(32,5,1,39,47,70,31), T=(0,5,129,39,47,199,159,27). Parent mask 5, child mask 5. The only changed old arc is 0 to 5; vertex 7 is new, beating exactly core vertices 0,1,3,4. It beats every module vertex.
- C: P=(2,4,41,3,15,91,31), T=(0,5,41,131,15,219,159,23). Parent mask 27, child mask 155. The only changed old arc is 0 to 1; vertex 7 is new, beating exactly core vertices 0,1,2,4. It loses to every module vertex.

For each k, construct the parent P[M;k]. Reverse the indicated old arc and add vertex 7, with the stated incidences to the module. In both families old vertex 0 is now a universal sink. Delete it. The result R_k has the same order as the parent, and all arcs between retained old vertices are unchanged. Equivalently, this is replacement of one old vertex by one new vertex; the intermediate reversal is incident to the discarded vertex.

For family B, for every k>=0:

    H(P_k) = 20*2^k+5,       C(P_k)=1,
    H(R_k) = 40*2^k+11,      C(R_k)=3.

For family C, k=0 gives H(P_0)=25 and H(R_0)=51. For every k>=1:

    H(P_k) = 27*2^k-1,       C(P_k)=1,
    H(R_k) = 54*2^k-1,       C(R_k)=3.

Thus H(R_k)=2H(P_k)+1 in both families, with no change in vertex count. In particular both the parent and the replacement are strong: a positive Hamiltonian cycle passes through every vertex. The weak intermediate child is permitted and has the same path count as its sink deletion.

## A strong-parent, strong-child one-reversal/one-birth family

A third family A uses:

    P=(2,24,3,37,13,23), M_P=15;
    T=(2,8,67,101,15,87,19), M_T=79.

Reverse the old arc 1 to 4 and add vertex 6. It beats core vertices 0,1,4 and loses to every module vertex. No old vertex is removed. Both cores are strong, and each module vertex has an in-neighbor and an out-neighbor in its core, so every inflated tournament is strong.

At k=0 the counts are H(P_0)=25, H(T_0)=51. For k>=1:

    H(P_k)=8*3^k+17*2^k+1,
    H(T_k)=16*3^k+34*2^k+3=2H(P_k)+1.

At every k>=0 the cycle counts are C(P_k)=2^k+1 and C(T_k)=2^(k+1)+1. This family preserves all old vertices and changes precisely one old arc, but its cycle counts are not bounded by 3 as k grows. No minimality assertion beyond the separately audited k=0 ancestry obstruction is made.

## Exact module run expansion

Fix a core Q of order n and a mask M. For each permutation of the core vertices, consider its n+1 linear gaps, including the two end gaps. A gap may be empty exactly when its neighboring core vertices form a forward edge, with either boundary always allowed. A gap may be occupied by a nonempty module run exactly when its left core endpoint, if any, beats M and M beats its right core endpoint, if any.

Multiply (e_g+f_g z) over all gaps, where e_g and f_g are these zero-one permissions, and sum over the core permutations. Write the resulting polynomial as a_0+a_1 z+...+a_(n+1) z^(n+1).

Every Hamiltonian path in Q[M;k] has a unique core order and unique set of occupied gaps. For r occupied gaps, assign the k labeled, transitively ordered module vertices to those r gaps, using each gap. Each gap then has its unique increasing internal order. The number of assignments is r! S(k,r), where S denotes a Stirling number of the second kind. Conversely each such assignment produces exactly one legal path. Therefore

    H(Q[M;k]) = sum_r a_r * r! S(k,r).

Here 0!S(0,0)=1 and S(k,0)=0 for k>0. In particular the constant skeleton coefficient only governs k=0; it must not be extrapolated as a constant contribution for positive k.

For cycles, fix core vertex 0 as the first vertex and permute the other core vertices. Use the n cyclic gaps, including the final-to-first gap, with the same permissions. This removes rotational multiplicity. The identical argument gives the cycle expansion with coefficients c_r. Fixing a core start is valid even when module runs are present because the core is nonempty.

## Finite coefficient certificate

The complete skeleton coefficient lists, with trailing zero coefficients omitted, are:

- A parent paths: (25,59,41,8); child paths: (51,119,82,16).
- A parent cycles: (2,3,1); child cycles: (3,5,2).
- B parent paths: (25,45,20); child paths, before and after deleting 0: (51,91,40).
- B parent cycles: (1,1); intermediate child cycles: (0); replacement cycles: (3,3).
- C parent paths: (25,53,27); child paths, before and after deleting 0: (51,107,54).
- C parent cycles: (1,1); intermediate child cycles: (0); replacement cycles: (3,3).

These are finite exact identities, obtained by every core permutation and every permitted gap choice, not a fit to finitely many module sizes. verify.py recomputes all coefficients literally, validates every tournament and the exact edge changes, and independently checks path and cycle counts by subset DP for k=0,...,4. The reduced-child mask is obtained by deleting old vertex 0 and shifting the remaining mask bits down by one.

In each path pair, the child coefficient vector minus twice the parent vector is (1,1,0,...). This contributes exactly 1 for every k: at k=0 only r=0 contributes, and for k>0 only r=1 contributes. Thus doubling plus 1 holds at all module sizes. The displayed formulas follow from

    1!S(k,1)=1,
    2!S(k,2)=2^k-2,
    3!S(k,3)=3^k-3*2^k+3        (k>=1).

The cycle vectors immediately give C1-to-C3 in B and C. Universal-sink deletion preserves H because every Hamiltonian path ends at that sink, and deleting it is bijective with Hamiltonian paths of the remaining tournament.

## Scope, inheritance, and reproduction

The starting H25-to-H51 one-edit/birth repair was found while pursuing the independently audited obstruction that no H51 tournament contains any induced H25 subtournament. The new operation changes ancestry by discarding one vertex (B,C) or changing an old arc (A). It therefore does not conflict with that obstruction. The module-run expansion is elementary and self-contained; the all-k identities here need no Busch bound or exhaustive tournament bank.

This is a positive constructor in three sparse families. Deletion is an explicit part of families B and C, so no minimum old-arc edit count is asserted when deletions are allowed. Family B is not closed under this rewrite: its child count is six larger than the next parent count 20*2^(k+1)+5. It does not show that an arbitrary Hn parent admits either binary child, nor does it establish a complete residue class or all sufficiently large odd counts. In B, only n=20*2^k+5 is covered; in C, only n=25 or 27*2^k-1 (k>=1) is covered.

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json
