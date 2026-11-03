# Strong tournament instructions with an encoded integer state

Status: exact constructions and conditional state transitions, not a proof that every allowed odd integer is a Hamilton-path count. The strong wrapper and path-cover identities are from [S159's positive profile operator](https://github.com/eliottcassidy2000/math/blob/7589463e54e20585aff6a2e68019ee4a4af582cd/05-knowledge/results/session159_constructions_20261103/profile_operator/NOTE.md); the one-apex recurrences are from [S158's binary bank](https://github.com/eliottcassidy2000/math/blob/7589463e54e20585aff6a2e68019ee4a4af582cd/05-knowledge/results/session158_creative_20261103/binary_bank/NOTE.md). No priority or literature novelty claim is made.

Write h(X) for the number of directed Hamilton paths and F_X(r) for the number of ordered r-slot path covers, empty slots permitted. Let T_n be the transitive tournament on n vertices. Then F_Tn(r)=r^n. For nonempty X, let G(X) be the directed triangle substitution X→a→b→X. It is strongly connected and

    h(G(X))=h(X)+F_X(2),
    F_G(X)(r)=r² F_X(r)+r[F_X(r+1)−F_X(r−1)].          (1)

The latter holds as a polynomial identity, and F_(X⊕Y)(r)=F_X(r)F_Y(r) for an ordinal sum.

## A decoded integer instruction set

For n>=1 set E(n)=G(T_n), and set E(0) to a singleton. All E(n) are strong (taking a singleton to be strong), and (1) gives

    h(E(0))=1,       h(E(n))=1+2^n for n>=1.             (2)

Thus the decoded state is n=0 when h=1, and n=ν₂(h−1) otherwise. The map E followed by this decoder covers every n in N_0 injectively. For n>=1, replacing the transitive block T_n in its wrapper by T_(2n) or T_(2n+1) gives two exact digit instructions on actual witnesses. If x=h(E(n)), their count transitions are

    digit 0: x → 1+(x−1)²,     digit 1: x → 1+2(x−1)².  (3)

Start at E(1) with count 3 and process the remaining binary digits. For example, 3→5 or 9; from 5, the successors are 17 and 33. The zero state has a separate entry instruction to E(1); the digit instructions above require n>=1. The representation uses n+2 vertices, so its size is linear in the decoded integer, not necessarily in its bit length. These rules realize n as a decoded label, not h=n and not every odd integer as h.

## Ordered rational instructions on the decoded state

Take any finite ordered list of positive rational fractions a_i/b_i in lowest terms. At state E(n), n>=1, find the first fraction with b_i dividing n, and replace T_n by T_(a_i n/b_i). If none is applicable, halt. This exactly simulates the ordinary first-applicable FRACTRAN rule on the positive integer n. On the count x=1+2^n, the corresponding partial transition is

    x → 1+((x−1)^(1/b_i))^a_i.                           (4)

The b_i-th root is integral exactly when b_i divides n. The guard concerns n, not divisibility of x by b_i. The rule is an explicit state-level construction; no efficiency or complexity claim follows.

For example, with ordered fractions (3/2, 5/3), input n=12 follows

    12 → 18 → 27 → 45 → 75 → 125  (halt).

The corresponding counts are 1+2 raised to each displayed exponent. The first fraction takes priority when both fractions are applicable. This finite example is not a claim that the two-fraction program is universal, and replacing blocks is not a bounded-size local graph rewrite.

## A genuine count automaton and an obstruction to scalar state

S158's strong one-apex bank has F(n)=n+A(n), n>=1, with A(1)=2 and

    A(2n)=A(n)+n,    A(2n+1)=A(n)+1−(n mod 2).

An actual witness is obtained from the transitive base whose binary orientation word is the ordinary expansion of n followed by 0; an apex beats exactly the base vertices with bit 1. In state (h,n,t)=(F(n),n,n mod 2), the exact updates are

    digit 0: (h,n,t) → (h+2n,       2n,   0),
    digit 1: (h,n,t) → (h+n+2−t,     2n+1, 1).             (5)

Start at (3,1,1). The extra state cannot in general be dropped: F(2)=F(3)=5, whereas their digit-0 successors have counts F(4)=9 and F(6)=11. Thus a deterministic digit-0 rule on the scalar count alone cannot reproduce this bank. Distinct witnesses with the same count need not react alike to the same encoded instruction. The count outputs are exact tournament counts, but their range is only this restricted bank.

## An exact strong affine lane

Formula (1) gives

    F_E(n)(r)=r^(n+2)+r[(r+1)^n−(r−1)^n],
    c_n=F_E(n)(2)=2^(n+2)+2(3^n−1).

For n>=1,k>=0 define W(n,k)=G(E(n)⊕T_k), with T_0 empty. This is strong, and the ordinal-sum profile identity yields

    h(W(n,k))=1+2^n+2^k c_n.                              (6)

Hence, with n retained, k→k+1 acts as h→2h−(1+2^n). For n=1 the lane begins 15,27,51,99; for n=2 it begins 37,69,133. These are exact counts of individual strong tournaments, but the conditional affine lane is not an unrestricted addition instruction.

## Complementary dyadic tiles and the missing assembler

Every odd y>=3 can be written uniquely as y=2m+1, m=2^k+r, 0<=r<2^k. The boundary r=0 is h(E(k+1)); the interior is not automatically realized. More generally the bits of m give the exact arithmetic identity

    y=1+Σ_{i: bit_i(m)=1} (h(E(i+1))−1).                  (7)

This is a sum over a tuple of witnesses, not a count formula for one tournament. For example 7=1+(3−1)+(5−1) and 21=1+(5−1)+(17−1); none of these arithmetic identities supplies a tournament of count 7 or 21. A putative universal count-level assembler must establish how to realize such sums, and cannot be inferred from the existence of separate witnesses. Similarly, encoded universality in (2)–(4) and modular universality of other banks do not establish exact odd-count coverage.

## Independent finite controls

Run `python3 check.py` and `python3 -O check.py`. The checker constructs literal directed graphs and counts Hamilton paths independently by subset/endpoint dynamic programming, checks the profile formula and W(n,k) for small inputs, and separately verifies actual one-apex digit witnesses. It checks (3), first-applicable FRACTRAN transitions, integral-root guards, and dyadic arithmetic. These finite tests control examples; the displayed all-order statements follow from the proofs above and their cited input identities.
