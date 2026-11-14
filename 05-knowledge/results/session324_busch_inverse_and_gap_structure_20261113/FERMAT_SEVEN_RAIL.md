# Seven times every noninitial Fermat number is a strong count

**Theorem.** Every integer 7(2^s+1) with s≥2 is the Hamilton path count of a strong tournament. For s=1 it is 21, the known exclusion. In particular all 7F_r with r≥1 occur, with no Fermat primality assumption.

## A single repeatable rail from exponent 10

Use the [cut-one one-return-star formula](../session317_exact_anchor_atlas_20261113/ONE_RETURN_STAR.md). Start with the ordinal join I→U→W and add apex v; it is dominated by U, beats a specified set S in I, and has a specified cut in W. Use

| block | outgoing masks | apex out-cut |
|---|---|---:|
| I | (4,1,10,67,15,31,55) | 81 |
| W | (8,1,3,6) | 13 |

In I, let α count ordered two-path covers in which the second slot is nonempty and starts in S, and λ count those where the first slot is nonempty and ends in S as well. In W, let β count ordered two-slot covers in which the first slot is empty or ends at a vertex beating v. Put d=H(W+v). Literal subset DP yields

    H(I)=29, α=256, λ=74;
    H(W)=5,  d=13,  β=14.

The exact star law is

    H = αβ Q_U(2) + [H(I)d−H(W)λ]H(U)
      = 3584 Q_U(2) + 7 H(U).

For clarity, delete v and keep the pre- and post-v paths. Both meeting U gives αβ[F₂(U)]. U entirely before v gives αβ H(U) from cases with a nonempty I suffix, plus H(I)d H(U) when that suffix is absent. U entirely after v gives [αβ−H(W)λ]H(U), the subtraction enforcing the I-side incoming guard when there is no W prefix. Their sum is the displayed formula. The quotient is strong: I and W are strong, the apex enters I through S, the ordinal links reach U and W, U beats the apex, and W can reach its unique in-cut vertex 1, which beats the apex.

Take U=TT_m for m≥1. Then H(U)=1 and Q_U(2)=2^m, giving exactly 7(2^(m+9)+1). Adding the next transitive vertex at the end of U gives the next exponent with all old arcs unchanged. This is a literal fixed-core infinite rail; the offset +7 is an exact cancellation 29·13−5·74, not a missing count used as a homogeneous component.

## The finite initial bridge

The following outgoing masks give all remaining s=2,…,9; strongness and path counts are replayed directly.

| s | H | outgoing masks |
|---:|---:|---|
| 2 | 35 | (54,64,90,83,2,94,17) |
| 3 | 63 | (238,220,216,16,225,206,8,72) |
| 4 | 119 | (132,73,10,65,111,79,5,126) |
| 5 | 231 | (150,176,194,215,164,77,147,32) |
| 6 | 455 | (48,173,41,209,38,200,151,21) |
| 7 | 903 | (440,273,483,486,460,146,291,322,32) |
| 8 | 1799 | (450,16,107,163,109,387,42,342,94) |
| 9 | 3591 | (282,680,699,32,682,897,447,777,542,73) |

The family supplies an explicit infinite prime-cofactor test family whenever 2^s+1 is prime, and also realizes every composite member. It neither proves infinitely many prime cofactors nor covers 7p for arbitrary primes.
