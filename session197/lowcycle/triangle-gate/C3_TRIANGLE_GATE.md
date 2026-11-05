# A physical triangle gate and a two-parameter three-cycle family

Session197. **All-order identities and a constructive family, not a complete scalar classification.** The ordinary one-seam calculus and binary gadget function g are inherited from `ROOTED.md`. The new ingredient here is an explicit marked triangle replacement at an outdegree-one vertex, including its two-cover response. The two-cover case split was independently derived by the low-cycle and ribbon workers before comparison. No literature-priority claim is made.

For a rooted graph (T,a), use h=H(T), e_a for paths ending at a, s_a for paths starting at a, and R_a for ordered spanning path covers (P,Q) with P possibly empty and Q starting at a. A singleton path is allowed; H(empty)=1.

## 1. General triangle replacement, with exact sidecars

Suppose a has exactly one out-neighbor b in a tournament T of order at least two. Replace a by a directed triangle, each of whose vertices has a's original incidence to every outside vertex. Retain one triangle vertex as the new root x. Write T_triangle for this actual graph. Then

    H(T_triangle) = 3h + 3R_b − 3e_b,
    e_x(T_triangle) = e_a + R_b − e_b,
    R_x(T_triangle) = 2h + R_a + 4R_b − 2e_b,
    H(T_triangle−x) = h + R_b − e_b.                 (1)

If T is strong, then C(T_triangle)=3C(T). If |T|>=3, then s_x(T_triangle)=3s_a(T). In particular replacing the outdegree-one vertex of a strong C1 tournament produces a strong C3 graph with the displayed returning response.

**Proof.** Put S=T−a and define

    d0=H(S), f=e_b(S), K=R_b(S), k0=H(S−b).

All vertices of S except b beat a, while a beats b. Direct insertion of a into an old path or a two-cover gives

    e_a=d0−f, h=d0−f+K, R_a=d0+K,
    e_b(T)=k0, R_b(T)=2K−k0.                        (2)

For the last identity, in a cover whose reserved path starts at b, a must be the last vertex of one of the two paths. Placing it on the unreserved path contributes K. Placing it on the reserved path contributes K−k0, excluding the old singleton reserved path b.

Every exit from the triangle module to S must enter b. Hence at most one such exit can occur in a Hamilton path or in a spanning two-cover.

A Hamilton path traversing the module in one chunk gives 3h possibilities. For two nonempty chunks, S splits into an arbitrary preceding path and a path starting at b which cannot be the singleton b. There are K−k0 such splits and six ordered two-path partitions of a directed triangle. Thus H(T_triangle)=3h+6(K−k0), giving (1). Deleting x leaves a transitive two-vertex module, with one single-chunk path and two ordered nonempty two-path partitions. This gives H(T_triangle−x)=h+2(K−k0).

A path ending outside the module must use one module chunk, hence there are 3(h−e_a) such paths. The triangle's cyclic automorphism makes the three module endpoint counts equal. Subtraction from H(T_triangle) and division by three gives the stated e_x.

For R_x, let Q be the path starting at x and P the other path. There are three disjoint cases.

- No module-to-S exit is used. Then Q lies entirely in the module and all S vertices lie on P, possibly followed by a module tail. There is one whole-module path starting at x and two proper partitions into P-tail and Q. The contribution is d0+2(d0−f)=3d0−2f.
- Q uses the exit. Partition the three module vertices into the P-tail, the nonempty Q-prefix starting at x, and the Q-tail. There are seven such ordered path partitions. If the S portion of Q is singleton b, its module tail must be empty; only three partitions remain. The contribution is 7K−4k0. Empty old P is allowed and causes no extra exclusion.
- P uses the exit. Its module middle chunk is nonempty, its module tail may be empty, and Q is a module-only path starting at x. There are four such partitions, reduced to two if P's S-suffix is singleton b. The contribution is 4K−2k0. Its old S-prefix may be empty.

The total is 3d0−2f+11K−6k0, which equals the R_x expression in (1) by (2). These cases include T of order two. The verifier also literally enumerates the seven/three and four/two module partitions.

For a Hamilton cycle, the module must form one chunk because every exit goes through b. Contracting that chunk gives a Hamilton cycle of T, and there are three possible internal directed Hamilton paths of the triangle. This proves the cycle law for strong T. A Hamilton path starting at x has one S-path starting at b, with a module prefix and optional module tail. When |S|>=2, the S-path cannot end at b; exactly three module partitions are available. Since s_a(T)=s_b(S), this gives s_x=3s_a. The size-two boundary is deliberately excluded from this last law. QED.

## 2. Literal fixed-two-port implementation

Start with fixed ports a->b. For u>=2, write u−1 in binary from most significant bit to least. For each bit append a vertex with mask1 if the bit is1, and mask3 if it is0. A mask records the ports beaten by the newborn; every older nonport beats it. The two original ports are always retained. Call the resulting graph U_u.

Writing the new vertices x_1,...,x_m, they are transitive in that order, all beat a, and b->x_i exactly when the corresponding binary bit is1. The leading bit is1, so a,b,x_1,...,x_m,a is a Hamilton cycle. Deleting a leaves a transitive base with one apex b. Set

    D=g(u), p=u mod2.

Here g is the ordinary binary gadget function

    g(1)=1,
    g(2l)=g(l),
    g(2l−1)=g(l)+l−(l mod2)  for l>=2.             (3)

Indeed reversing the one-apex graph's arcs and reversing the base order gives exactly the allowed gadget indexed by u in `ROOTED.md`. The old b-ending count is p, and its b-starting two-cover count is 1 plus the binary integer u−1, namely u. Therefore

    H(U_u)=D+u−p, e_a(U_u)=D−p, R_a(U_u)=D+u,
    e_b(U_u)=1, R_b(U_u)=2u−1, s_a(U_u)=C(U_u)=1.  (4)

Append masks2,3. The old a and these two newborns form a directed triangle, each beating b and losing to all x_i. Thus this is literally the replacement in Section1. At the old, still-exposed a, the seed profile (h,e,R,d) is

    h=3D+9u−3p−6,
    e=D+2u−p−2,
    R=3D+11u−2p−6,
    d=D+3u−p−2,                                   (5)

with s_a=C=3. This is an actual fixed-port gate, not a profile inferred from equal scalar counts.

## 3. Arbitrary second binary phase

After the gate, append an arbitrary word over masks0,1 which is empty or ends in1. Encode it by

    v=1+sum_i z_i 2^(i−1),

where z_i is its i-th mask, so v>=1 and v=1 means the empty suffix. Every integer v>=2 gives exactly one such suffix.

For any rooted graph with more than one vertex, let d=H(T−a). Keeping old a exposed, a mask0 birth O and mask1 birth Z have the exact four-coordinate laws

    O(h,e,R,d)=(h,0,2R,d),
    Z(h,e,R,d)=(h−e+R,d,2R−d,d).                   (6)

These inherited old-root one-seam sidecars also follow by splitting paths and covers at the newborn. Both preserve s_a. An O output is weak; a Z output has cycle count s_a of its input, since its only reverse exit is newborn->a. Therefore every empty or final-Z suffix in this construction is strong with exactly three cycles. Intermediate weak states do not invalidate the final conclusion.

For an arbitrary input profile, the scalar suffix response is

    h_out=h+(g(v)−1)d+(v−1)(R−d)
              +(1−t)(d−e),  t=v mod2.             (7)

For v=1 this is h. For a first mask1, write v=2l and substitute Z followed by the suffix encoded by l. For a first mask0, write v=2l−1 with l>=2 and substitute O followed by that suffix. Equations (3),(6) give (7) in both cases, proving it for every word, including the empty suffix.

Substituting (5), the family has the closed formula

    E(u,v) = [g(u)+3u−p−2]g(v)
             +[2g(u)+8u−p−4]v − u(1+t) − p,       (8)

for every u>=2,v>=1. Its order is

    4+bit_length(u−1)+bit_length(v−1).

Every graph has feedback number at most two, since removing the two fixed ports leaves the transitive birth order. No claim that (8) exhausts C3/FVS<=2 graphs is made.

## 4. Exact family inverse and useful boundaries

For all allowed parameters,

    E(u,v)>=7uv+1.                                (9)

For u>=4, using g(u),g(v)>=1 and p,t<=1 gives E−7uv>=(u−3)(v+1)>=2. At u=2, g(2)=1,p=0 gives E=5g(v)+14v−2(1+t)>=14v+1. At u=3, g(3)=3,p=1 gives E=9g(v)+25v−3(1+t)−1>=25v+2>21v.

Thus to decide whether a particular positive integer X belongs to this family, test precisely the finitely many u>=2,v>=1 with uv<=(X−1)/7 using (3),(8). The construction then gives a graph for every successful pair. This is a complete inverse for this named family, not a cofinite theorem for the full low-cycle spectrum.

Examples are

    E(4,2)=67, with fixed-port word11231 and order7;
    E(8,7)=533, with word11123011 and order10;
    E(6,31)=1861, with word1312301111 and order12.

The H533 graph is the separately audited `FEEDBACK533_C3.md` witness after relabeling. Together with the independent negative theorem there, it demonstrates why allowing a third cycle changes the fixed-port scalar class.

The exact family image through5001 has 1581 distinct values. The bound (9) makes this image calculation complete for (8). In particular821 is not in this family. This last absence is not asserted for all fixed-port C3 graphs or all tournaments; it is an explicit remaining obstruction to using this family alone as an arithmetic compiler.

## 5. Verification

Run `python3 verify_c3_triangle_gate.py` and `python3 -O verify_c3_triangle_gate.py`. The standard-library verifier checks all1390 outdegree-one incidences in labelled tournaments through order5, all module-partition constants, exact physical family graphs and profiles, the three examples, the exact image/inverse bound, and normal/optimized receipt equality. The two-cover counts use subset start/end DP on literal graphs, not formula fitting. The proof is all-order; these controls are separate finite evidence.
