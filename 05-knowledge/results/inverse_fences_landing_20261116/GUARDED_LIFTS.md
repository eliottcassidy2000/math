# Physical rays for both guarded signed inverse moves

Status: all-order within explicit marked graph families; not a full tournament spectrum proof. When substituting the arithmetic incoming branches in the half-count chart, both targets are 4g−1: for g odd the signed parent is 3g; for g even it is 3g−1. The physical problem is not to divide an arbitrary tournament by three.

## The reusable two-apex counter

For a 0/1 word w, put its vertices in a transitive spine and give a distinct apex x the outgoing edges x→i exactly at 1 bits. A word beginning with 1 and ending with 0 gives a strong tournament. Write h(w) for the Hamilton count, A(w) for the binary integer, and Q(w)=1+Σ_(i≥0) (1−w_i)2^i with i counted from the left. With p the last bit (or 0 for empty), the exact update is

    1: (h,A,p) ↦ (h,2A+1,1);
    0: (h,A,p) ↦ (h+A+p,2A,0);

starting at (1,0,0). One explanation is the odd directed cycle formula: all these cycles use the apex, and a new zero, cinched onto the end of the spine, admits an apex-led increasing fragment from an earlier outgoing position. The number of odd cycles newly closed is (A+p)/2, by the subset parity of intermediate spine positions. A new one closes none.

For two strong words u,v, direct every cross edge u→v except the apex bridge, directed backwards y→x. Call this R(u,v). Its count is

    H(R(u,v)) = h(u)h(v)+(A(u)+1)Q(v)−1.

This is the two-apex bridge in the [preceding proof](../chanpak_target_ribbon_audit_20261116/TWO_APEX.md). Briefly, an odd cycle mixing blocks must use both apices and the reversed bridge; it blocks all original cycles. The two independent choices of increasing rooted fragments number A(u)+1 and Q(v), each with signed even-minus-odd balance −1. This gives ((A(u)+1)Q(v)−1)/2 mixed odd cycles, adding twice that number to the ordinal product h(u)h(v).

For any fixed word W starting and ending in 1, let A=A(W) and B=h(W)−A+1. For s≥1 a tail of s zeros gives

    h(W 0^s)=A·2^s+B,   A(W 0^s)+1=A·2^s+1.

Consequently for any fixed right core v with j=h(v), Q=Q(v), the exact dyadic pencil is

    H(R(W 0^s,v))=(j+Q)A·2^s+jB+Q−1.

These are all real directed graphs, and the input includes which part of its spine is adjustable.

## Plus branch: 3g ↦ 4g−1 with g odd

Take W=1001, so A=9 and B=−3; take s≥2. The source is the ordinal sum of a directed triangle and R(W 0^s,100). The target is R(W 0^(s−1),100000). The two right words have (j,Q)=(5,7) and (33,63), respectively. Therefore

    d=324·2^s−27,    H_target=432·2^s−37=(4d−3)/3.

The source is 3 mod 6. To carry it out, remove the triangle component, take one zero out of the left tail, and add three zeros to the right tail. This is a fixed-size surgery, not an induced extension. At s=2 this gives 1269→1691.

This solution is unique inside the following aligned-pencil stencil: fixed same W on both sides, right words 10^m and 10^m' with m,m'≥1, source preceded by an ordinal triangle, and the target left-tail length shifted by a fixed integer. Indeed for 10^m one has j=2^m+1,Q=2^(m+1)−1. After dividing the source by the triangle factor, matching the constant terms of H_target=4H_source_core−1 forces

    (B+2)(2^m'−2^(m+2))=3B−7.

Thus B+2 divides 13. From h(W)≤(A+1)/2 one has B odd and ≤1. (The general bound h(w)≤A(w)+1 for words starting in 1 is inductive: track the nonnegative D=A+1−h, noting that a 1 increases D by A+1 and a 0 subtracts 1 only if the preceding step was a 1.) The remaining B choices −1,−3,−15 give differences −10,16,4 respectively. A difference −10 is impossible between powers of two (its odd part is 5, incompatible with 2^r−1), and +4 is impossible with the subtracted power 2^(m+2)≥8. Thus B=−3 and m=2,m'=5. Now A−h(W)=4 forces A≤9; testing the possible odd binary integers 1,3,5,7,9 gives only W=1001. The leading terms then force shift −1.

## Minus branch: 3g−1 ↦ 4g−1 with g even

Take W=11 or 101, so A=3 or 5 and B=−1. For s≥1 use source R(W 0^s,1110) and target R(W 0^(s+2),10). For a right word 1^m0, j=Q=2^m+1. Hence

    d=18A·2^s−1,    H_target=24A·2^s−1=(4d+1)/3.

The source is 5 mod 6. To carry this out, move two vertices from the right fan to the left zero tail, with the prescribed incidences and no change in total order. First examples are 107→143 and 179→239.

In the aligned stencil with fixed W and rights 1^m0 and 1^m'0, the count anchor is (B+1)j−1. The inverse map fixes −1 formally, so matching anchors gives (B+1)(3j'−4j)=0. Both j,j' are odd and B=−1 is forced. Matching slopes up to a power of two now asks j=3j'. The equation 2^m+1=3(2^m'+1) gives m'=1,m=3 (if m'≥2, the right side minus one has 2-adic valuation 1 and is larger than 2). The tail shift is then +2, and A−h(W)=2 forces A≤5, giving exactly the two W above.

## Why concatenation is not automatic

A fixed dyadic stencil has count functions a·2^s+b. In any affine identity of two such functions after shifting s by a fixed integer, both signed inverse lifts (slopes 4/3) decrease the valuation v₃(a) by exactly one. The even lift (slope 2) leaves it unchanged. Consequently no finite bank of these aligned identities can support arbitrarily long words with more than max v₃(a)−min v₃(a) sign inverses, no matter how many even steps are interspersed. This concerns infinitely valid pencil identities, not individual isolated coincidences. The minus anchors are −1 and the plus anchors differ; matching absolute leading profiles, not just scalar congruence, is necessary. A true ternary refill, a nonaligned state-dependent reader, or unboundedly many templates are necessary in this mode.

## Arbitrary right cores: the exact ports and their obstruction

The same single backward bridge R(U,V) works for any nonempty tournaments U,V, with roots x,y. For a root r, define a weighted end port

    P_end(T,r)=Σ_{S∋r} end_r(T[S]) H(T\S),

and P_start by reversing paths. Here end and start count Hamilton paths with the given endpoint, and H(empty)=1. The signed version of P_end inserts (−1)^|S| into the sum and is exactly −start_r(T); the signed P_start is −end_r(T). The bridge count therefore reads

    H(R(U,V)) = H(U)H(V) + P_start(U,x)P_end(V,y) − end_x(U)start_y(V).

Indeed a mixed odd cycle decomposes uniquely into a rooted outgoing fragment in U and an incoming one in V, whose total vertex count is odd. The choice of either odd/even type yields the product-minus-signed-product formula, and the remaining cycle packs on the two sides factor by H on their complements. The only disallowed rooted pair (both roots alone) has even total size, so contributes zero. Old packs still give H(U)H(V).

Proof of the signed port identity: work in the commutative ring with x_v²=0, put X=diag(x_v) and u=1ᵀX. The row G_A(x)=u(I−AX)^−1 enumerates all nonempty simple paths by their last vertex (a repeat vanishes in the ring); F_A=1+G_A1 enumerates them without an endpoint. Reversal identifies F_{Aᵀ}=F_A. The tournament identity Aᵀ=J−I−A gives I−AᵀX=(I+AX−JX)(I+X). Premultiplying its inverse by u, the factor (I+X)^−1 disappears because uX=0. The remaining matrix has rank-one modification JX=1u, so

    G_{Aᵀ}(x)=u(I+AX)^−1/[1−u(I+AX)^−1 1]=−G_A(−x)/F_A(−x).

Summing and adding 1 gives F_{Aᵀ}=1/F_A(−x). The full squarefree coefficient in G_{Aᵀ}(x)=−G_A(−x)F_A(x) is exactly the asserted signed identity.

The legacy negative-anchor condition P_end(V,y)=H(V) and start_y(V)=1 is rigid. Indeed P_end enumerates all ordered pairs of complementary directed covering paths where the first ends at y and the second may be empty, while H(V) counts just those whose seam is valid. If P_end=H, every Hamilton path of V−y must start at an out-neighbor of y; otherwise the singleton first part y gives an invalid seam. Thus start_y(V)=H(V−y), which has value 1 only for a transitive deletion. For a strong V the one-apex equality classification then gives just the right fan 1^m0. In particular one cannot freely replace this same two-port handoff by an arbitrary ternary reservoir.

The exact wider guard is **compensation**, not equality of each constituent: with W=11 or 101 on the left as above, negative anchor −1 holds precisely when

    P_end(V,y)−H(V) = start_y(V)−1,

and then

    H(R(W 0^s,V)) = A(W)(H(V)+P_end(V,y)) 2^s −1.

This is a strict extension of the fan lane. For example, outgoing masks (6,8,2,5) with root y=0 give (H,P_end,start)=(5,6,2), hence the odd coefficient 11 in H+P_end. The remaining physical question is to grow H+P_end by exact factors of 3 repeatedly while maintaining this compensated port.

## Three consecutive plus lifts with strong witnesses

There is a larger positive lane; in particular the triangle-ordinal input of the first construction is not necessary. Fix W=11 and define four strong right cores V_i whose root is vertex 0. Row masks below use bit j for the edge i→j.

| i | Rows of V_i | (H, P_end, start_0) | Left tail | R count |
|---:|---|---|---|---|
| 0 | (56,37,1,6,14,28) | (37,44,4) | s | 243·2^s+3 |
| 1 | (44,41,2,4,15,28) | (23,31,5) | s+1 | 324·2^s+3 |
| 2 | (16,9,3,5,14) | (15,21,3) | s+2 | 432·2^s+3 |
| 3 | (50,40,3,5,14,28) | (43,53,7) | s+1 | 576·2^s+3 |

In each row the anchor is −H+P_end−start=3, and the coefficient from the arbitrary-core bridge is 3(H+P_end) before shifting the tail. Thus for every s≥1 the four displayed counts give three successive applications of d↦(4d−3)/3, each with legal parent 3 mod6. At s=1 the chain is 489→651→867→1155. The same tabulated right cores and shifts work for W=101, scaling all leading coefficients by 5/3. These are actual common templates for successive physical moves; no separate witness search is needed at each s. Both sides of each bridge are strong, and a backwards arc plus the remaining forwards cross arcs make the full bridge strong. One changes the right core and at most one zero in the left tail at each step. The checker independently verifies all port profiles and strongness. This still consumes one factor of 3 from the odd scale at every step; the examples do not supply ternary refill.

## Four consecutive minus lifts with strong witnesses

Here the same compensated anchor persists through four negative steps. Again the right root is vertex 0 and all masks are ordinary outgoing rows.

| i | Rows of V_i | (H, P_end, start_0) | W | Left tail | R count |
|---:|---|---|---|---|---|
| 0 | (92,81,34,6,12,27,60) | (109,134,26) | 101 | s | 1215·2^s−1 |
| 1 | (52,49,2,7,12,28) | (37,44,8) | 101 | s+2 | 1620·2^s−1 |
| 2 | (44,17,2,6,13,30) | (25,29,5) | 101 | s+3 | 2160·2^s−1 |
| 3 | (116,97,34,7,14,24,60) | (131,157,27) | 101 | s+1 | 2880·2^s−1 |
| 4 | (22,16,2,7,12) | (9,11,3) | 11 | s+6 | 3840·2^s−1 |

Every anchor is −1. Each parent is 5 mod 6 for s≥1. Each step d↦(4d+1)/3 trades one factor of 3 from the odd coefficient for a factor of 4 from the tail. In particular s=1 gives 2429→3239→4319→5759→7679. The right cores are all strong; the watcher can retain the first s existing zero-tail vertices at every step, since only 0 to 6 additional zeros are specified in the finite recipes. The checker independently computes all profiles, directed reachability, and one full literal chain.

### A three-step ternary reservoir, not yet an unbounded refill

The opposite-port even/odd sums r,s in [the near-sink proof](INVERSE.md) give P=r+s and E=s−r for the converse of the chosen right module. Thus the compensated condition P−H=E−1 is exactly **2r=H−1**, the one-hit (not necessarily repeatable) near-sink condition. It is a physical bridge between the two inverse architectures. The right-core coefficient K=H+r+s has the following distinct witnessed values. Use the converse of each listed rooted tournament as the right module; masks are outgoing.

| K | H | (r,s) | Root | Rows before taking converse |
|---|---:|---|---:|---|
| 2 | 1 | (0,1) | 0 | (0) |
| 6 | 3 | (1,2) | 0 | (2,4,1) |
| 18 | 9 | (4,5) | 1 | (2,4,25,19,3) |
| 54 | 25 | (12,17) | 2 | (4,33,34,7,15,25) |

With left word 101 0^s, the general bridge formula gives d=5·2^s K−1 for s≥1. Moving 54→18→6→2 and adding two zeros to the left at each step gives **three consecutive exact minus inverse moves**, since (4d+1)/3 replaces K by K/3 and s by s+2. At s=1 the path is 539→719→959→1279. This uses genuinely different right profiles and permits the final singleton right core (the bridge proof remains valid with empty complements and paths). It is not an unbounded generator; a physical construction for higher successive powers of 3 in K would be a direct way to continue.
