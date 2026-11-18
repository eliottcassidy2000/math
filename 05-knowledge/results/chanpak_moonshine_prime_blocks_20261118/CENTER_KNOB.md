# A center move that preserves the exact prescribed difference

**PROVED physical construction.** This adds an independent center choice to the preceding unit-gain difference controller. It does not yet place that center at every prescribed integer.

Retain its literal sources L=(9,1,23,5), R=(9,1,25,5), r=4j+2 with j>=1, u=(7*16^j-2)/5 and gap constant C_j. Let a=r-3 and b=r+1. Both positions are free in the edit mask: u has a zero at a, and the b position is outside u and the distinguished r position.

For an arbitrary high address V>=0 define two pairs of addresses:

    lower-center: WL=2^r+2^a+2^(r+2)V, WR=u+2^a+2^(r+2)V;
    upper-center: WL=2^r+2^b+2^(r+2)V, WR=u+2^b+2^(r+2)V.

**Theorem.** Both pairs have exactly the same signed count difference

    HL-HR=C_j-6-2V.

Moving from the lower-center pair to the upper-center pair increases BOTH counts by exactly

    Delta=345*2^(r-3)+75V+10-5(V mod2)>0.

Thus choosing a requested even gap g<=C_j-6 fixes V=(C_j-6-g)/2, while the bit swap still provides two different absolute placements with that exact gap. No source count is inferred from arithmetic: both members are supplied by the two literal order-five sources and binary words.

## Proof of equal gap response

For a free bit i the reduction of the gap constant when that bit is switched on is

    c_i=2*2^i+5 sum_(j in edit support) z_j(2^(|i-j|-1)+[|i-j|=1]).

The r+1 bit has c_(r+1)=A+5=6, since A=1 and z_r=1. The recursive mask near its upper end has the same free-bit response c_(r-3)=6. This can also be checked without trusting a pattern: u_(j+1)=16u_j+6 transports an old free bit i to i+4. The new low bits 1 and 2 are negative edits. Their contribution to c_(i+4) is -5(2^(i+2)+2^(i+1))=-30*2^i, exactly cancelling the change from 2*2^i to 2*2^(i+4). Hence old free-bit responses survive the four-bit extension. At the base r=6 the mask is (0,-1,-1,0,-1,0,+1), so c_3=16+5(-2-2-2+4)=6. Transporting that base response gives i=r-3 at every subsequent level. Both free bits therefore reduce the difference by exactly 6.

## Proof of the common center shift

It suffices to compare L's two counts, because their gaps agree. For the lower address, the low bits r-3 and r contribute 4 internally and 18V against the high tail. For the upper address, bits r and r+1 contribute 2 internally and 3V against the high tail, plus the adjacency term V mod2. Thus

    K(upper)-K(lower)=-15V-2+(V mod2).

The address rises by 15*2^(r-3), and neither address has bit zero. Substituting in H=9+23W-5K(W)-(W mod2) proves the stated Delta. Every term and sign is exact.

## Interpretation and remaining target

The common high address is a gap controller; the equal-response free-bit swap is a gap-neutral center controller. This is stronger than finding two unrelated equal gaps after enumeration: the independence is an all-order identity with a repeatable mask recursion. It also identifies the remaining challenge honestly: the center controller currently has only two positions, separated by a substantial amount. For universal placement we need a sufficiently rich bank of equal-response bit patterns whose common shifts fill the required interval, or source changes that align this interval with N. Coprimality or merely having two placements is not such a theorem.
