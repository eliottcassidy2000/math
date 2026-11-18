# Complete low-knob classification and a target-first inverse

**PROVED for the two-source periodic unit-gain family only.** This identifies all low-address choices, their exact continuations, and why fixed-continuation center freedom stops at two positions.

At level j>=1 put r=4j+2, B=2^r=4*16^j, u=(7*16^j-2)/5 and C=(299*16^j+36-10j)/5. A free position is any zero bit of u below r, together with r+1; the distinguished bit r is not free. For any subset of these positions let P be its binary address. Use

    WL=B+P+4BV, WR=u+P+4BV.

Let c_i denote the reduction in gap constant upon including free bit i. The exact difference is

    HL-HR = C - sum_(i selected)c_i - 2V.

## All response families, not a finite census

As a multiset, the free responses are

    {4,6,6} union {a_k,b_k:0<=k<=j-2} union {a_(j-1)},
    a_k=104*16^k-2, b_k=208*16^k-2.

There are exactly 2j+2 free bits. Proof: at j=1 direct evaluation gives {4,6,6,102}. Extending the mask by four places retains every old free response at its shifted position (the new negative edits at positions1 and2 cancel the factor16 change in the weighted-start term). The two new free positions 0 and3 have responses a_j and b_(j-1), respectively. For position0 this follows from c_0=-18+40(B-u)=104*16^j-2; for position3, c_3=-4+5(B-u)=208*16^(j-1)-2. Here the old low bit is zero for j>=1, so there is no extra adjacency term. This proves the induction.

Their total is

    S_j=(39*16^j-14-20j)/5,
    C-S_j=52*16^j+10+2j.

Every response is even. Therefore every one of the 2^(2j+2) low patterns has an exact nonnegative collision continuation

    V=(C-sum c_i)/2,
    26*16^j+5+j <= V <= C/2.

This gives the full physical collision bank without searching any high word. For a requested even gap g the same statement holds with V=(C-g-sum c_i)/2 whenever g<=C-S_j.

## The independent center knob is maximal at fixed V and gap

After the three small responses 4,6,6, every response in sorted order exceeds the sum of all preceding responses. Indeed, the sum before a_k is

    16+sum_(l<k)(a_l+b_l)=(104*16^k-24-20k)/5 < a_k,

and b_k>a_k plus that preceding sum. Consequently equality of two response subset sums forces agreement at every large response. Among {4,6,6}, the only ambiguity is exchanging which 6 is selected. Thus there are at most two low patterns with the same V and gap; the previously proved center swap exhausts them. More local searching cannot produce a third independent center in this fixed mask family. Changing V, the mask, sources or alphabet remains possible.

## A prescribed-height decoder

For gap zero every supplied center H lies between 400B^2 and 1000B^2. To see this, use 6.5B<=V<7.5B, 0<=B+P<4B and K(W)<=W. The ordinary count 9+23W-5K(W)-(W mod2) is at least 18W and at most 9+23W. The lower bound exceeds 468B^2; the upper is below 1000B^2 for B>=64.

As j increases, B multiplies by16. These certified height bands are disjoint. Given an absolute target N, first select its sole possible level (if any) using these bounds. Then enumerate the 2^(2j+2) low patterns, solve the forced V, and compare one exact count with N. On a hit return both physical words and sources. This is a complete decision algorithm for this entire infinite family, using O(N^(1/4) log N) arithmetic operations with the direct binary checksum evaluator. No maximum tournament size is chosen in advance.

The same bounds imply at most O(X^(1/4)) distinct centers through X. This family alone therefore cannot prove Chan–Pak. The new information is a physical target-first selector with a completely classified continuation bank, and an exact explanation of which freedom is still insufficient. A productive extension must change the response spectrum or add genuinely new source states, rather than mistake the finite two-choice center ambiguity for a freely adjustable interval.

Replay: `python3 finite_edit_recontact/low_knob_replay.py` checks the response classification, total, all low patterns through j=5, and center inverse on every generated height through j=4. It uses the literal source controls already independently checked in the companion scripts.
