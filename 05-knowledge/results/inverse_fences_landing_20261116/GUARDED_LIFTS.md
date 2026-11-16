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
