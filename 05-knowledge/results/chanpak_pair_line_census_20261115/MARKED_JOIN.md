# A physical reader algebra with a central fork

**Proved:** exact ordinal compiler, a literal primitive two-option
fork that survives arbitrary surrounding marked hosts, and a
prime-divisibility budget for repeated blocks. Not an absolute
selector.

For a nonempty marked tournament B write h=H(B), q=P₂(B)+2h,
s/e for spanning paths starting/ending in the marked set, and:
M counts ordered two-slot covers with a marked start of slot 2;
E counts the corresponding marked end of slot 1;
W counts both requirements together. The unmentioned slot may
be empty for M or E. If A is a nonempty source, form B⇒A then add
an apex that beats exactly the marked vertices of B. The result
has H=(h−W)H(A)+M Q(A), Q(A)=P₂(A)+2H(A), by the marked-closer law.

## Ordinal compiler

For marked A⇒B, with the union of the marks, h and q multiply, and

    s = s_A h_B;                       e = h_A e_B;
    M = h_A M_B + M_A q_B;             E = E_A h_B + q_A E_B;
    W = h_A W_B + W_A h_B + M_A E_B + e_A s_B.

For M and E divide at which component contains the observed
endpoint. For W consider in turn both endpoints in B, both in A,
end of slot 1 in B and start of slot 2 in A, or the other way
around. The third case separates into two independent two-slot
covers, hence M_A E_B. The last forces each slot to consist of
a single component's spanning path, hence e_A s_B. This is why
dropping the Hamilton endpoints gives a wrong continuation state.

All these identities are represented at once by the integral matrix

    [ h  M  e  W ]
    [ 0  q  0  E ]
    [ 0  0  0  s ]
    [ 0  0  0  h ]

which multiplies on ordinal join. This lifts the old Heisenberg
analogy to a physical matrix compiler: the central cross term is
the usual triangular product, alongside an additional endpoint
term with a zero middle diagonal. The missing inverse operation
is real; one must realize a difference between two physical words,
rather than insert formal group inverses as building blocks.

## Same continuation, different intercept

The two order-6 marked host orientations (code, mark mask)
(257,51), (8737,29) have identical

    (h,q,M,E,s,e)=(17,224,135,135,11,11)

but W=76 and W=72 respectively. Codes use lexicographic i<j
with bit 1 meaning i→j. No class assumptions are needed:
marked_join.py checks literal arcs and performs exact subset DP.
An earlier same-host fork at order 5 has common
(15,120,42,42,6,6) with W=10 and 8, but shares bad prime factors.

The entire difference of the order-6 matrices is −4 times the
matrix unit in the top-right corner. This unit propagates
through arbitrary surrounding marked hosts by their ordinary H
counts only. Therefore a join of k order-6 blocks chosen
independently from the pair gives exactly k+1 distinct intercepts
with common slope and with spacing 4·17^(k−1)·H(A) on any common
source A. Every switch leaves every other profile coordinate
unchanged, even after further marked-host joins. The baseline
still grows, and the spacing grows too; no unit-width interval
is implied. At one block and a singleton source, the two
literal outputs are 211 and 215.

For k≥3 and a transitive source of t≥1 vertices the resulting
counts are all coprime to 17. Modulo 17 every diagonal h vanishes,
so the only contribution to M_k is M·q^(k−1), and the only
one to W_k is M·E·q^(k−2) (the zero diagonal kills the endpoint
route of length at least three). Thus the count modulo 17 is
135·224^(k−2)·(224·2^t−135), which cannot vanish because
3·2^t=16 mod17 would require 2^t=11 mod17, outside the subgroup
of powers of two.

## A prime budget for nonprimitive host joins

For each block i and prime p let γ_i=v_p(gcd(h_i,q_i)).
In the matrix product, any nonzero entry is a sum of products
of entries along weakly increasing index paths. A nonzero path
can take at most two strict steps: 1→2→4 or 1→3→4 for two
steps. Every diagonal step contributes h_i or q_i (the third
diagonal never contributes). Thus every entry is divisible
by p^(sum γ_i − largest γ_i − second largest γ_i), with absent
largest terms treated as zero. In particular the marked closer
inherits this divisor, independent of the source. The attractive
order-5 unit-half-step fork repeated k≥3 therefore forces a
factor 15^(k−2), whereas the primitive order-6 pair does not.
A low-denominator progression or a large window that ignores
the host's gcd budget can give no prime coverage at all.
