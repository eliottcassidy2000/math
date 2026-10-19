# An independent odd-prime one-subdivision construction

2026-10-19. **PROVED**. No external priority claim. This companion gives a second arithmetic mechanism independent of the power-of-two square-scaling proof.

Fix an odd prime p for which 2 is a quadratic nonresidue. On n>=1 original vertices indexed by i=0,...,n-1, put q_i=pi+1. Let L be the least odd power of p exceeding 4n+2. For every edge i<j define

    A_i=4Lq_i,
    d=q_j-q_i,
    r=Ld-1,
    X_ij=r^2-A_i.

The intended sums are r^2 and (r+2)^2. Since q_i+q_j is 2 modulo p, the p-adic valuation of A_i+A_j is the odd exponent of L. Thus the original family is independent. Every X_ij is 1 modulo p, so new-new sums are 2 modulo p and are nonsquares.

Suppose X_ij+A_k=s^2, with s>=0, and put t=s-r. Then

    2L(dt-2(q_k-q_i))=t(2-t).

Because p is odd, at most one of t,t-2 is divisible by p. Consequently p^a=L divides one of the two factors. On the other hand,

    |t|=4L|q_k-q_i|/(s+r)
        <=4Lp(n-1)/(Lp-1)<4n.

As L>4n+2, divisibility forces t=0 or2. The square equation therefore has exactly the intended endpoints k=i,j.

All new labels exceed every original label: d>=p, so

    X_ij >= (Lp-1)^2-4Lpn > 4Lpn >= max_k A_k,

using Lp>12n+6. For a fixed d, different i yield different X_ij. Distinct positive differences are distinct positive multiples of p; the minimum squared-root gap is 3L^2p^2-2Lp, which exceeds the maximum possible difference 4Lpn between subtracted original labels. Hence all labels are distinct.

Minimality of the odd exponent gives L<=p^2(4n+2)<=6p^2n. Therefore all labels are at most 36p^6 n^4. For p=3 the sharper direct estimate is 26244n^4. Thus every fixed eligible prime gives an absolute quartic construction, with a constant depending on p; p=3 is the smallest eligible prime.

The odd-prime proof uses valuation parity to suppress old-old edges, a fixed quadratic nonresidue to suppress new-new edges, and a prime-power divisibility condition to select endpoints. It does not preserve arbitrary numerical seed labels and is superseded for that purpose by the square-scaled pair-extension theorem.

The independent checker covers p in {3,5,11,13,19,29} and complete graphs through order twenty-four, verifying the full intended adjacency relation. These controls supplement the all-n proof and do not claim an exhaustive prime classification.
