# The automatic even-order transfer does not extend to modulo sixteen

The modulo-eight theorem cannot be upgraded to a universal odd-unit modulo-sixteen transfer merely from even core order and a nonzero Hamilton residue. Here is an exact six-name example, with outgoing adjacency masks:

    A=(48,17,59,19,32,10),
    B=(50,16,59,19,0,26).

Their ordered-cover profiles, at k=1,...,6, are

    q(A)=(11,178,828,1752,1800,720),
    q(B)=( 3, 90,612,1608,1800,720).

Thus H(A)-H(B)=8, nonzero modulo sixteen. In the directed-triangle ambient context, replace one root by A or B. Its exact substitution formula is

    H(C3[A,1,1])=3H(A)+q_2(A).

Consequently the two ambient counts are 211 and 99, with difference 112, which is zero modulo sixteen. An odd multiple of eight could not vanish modulo sixteen. Hence the general even-order transfer statement is genuinely limited to modulo eight unless further profile information is controlled.

This is not an exact Hamilton collision: 211 differs from 99. It is not a full named-response collision, and no claim is made that this pair arises from one original-module XOR family. It is a sharp boundary for the general ambient congruence theorem, whose hypothesis permits arbitrary same-size even module pairs.

`sharpness.py` independently enumerates every old word and every legal cut pattern through the binomial cut formula to recover both complete six-name profiles. It also counts both eight-name ambient Hamilton paths by literal permutation enumeration. Normal and optimized receipts agree.
