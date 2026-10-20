# Exact consequences once D(z)=z^4 is verified

Status: algebraic consequences of the complete deficit-two proof; no independent bridge-classification assertion.

Let a_0=1, a_n=0 for n<0, and a_n=a_(n-1)+a_(n-3) for n>=1. These count words in rigid strong atoms of sizes1 and3.

With D(z)=z^4, the complete paired-black leaf series is

    L(z)=(2z^4-z^5-z^7)/[(1-z-z^3)(1-z^2-z^6)]
        =(2z^4-z^5-z^7)/(1-z-z^2+z^5-z^6+z^7+z^9).

The mixed-leaf series is

    M(z)=z^4/(1-z^2-z^6).

Therefore, for all r>=0,

    L_(2r)=(a_(2r)-a_r)/2+a_(r-2),
    L_(2r+1)=(a_(2r+1)-a_r-a_(r-1))/2,
    M_(2r)=a_(r-2),
    M_(2r+1)=0.

The complete total-leaf recurrence, using zero coefficients at negative indices, is

    L_n=L_(n-1)+L_(n-2)-L_(n-5)+L_(n-6)-L_(n-7)-L_(n-9)

for every n>=8. It fails at n=7 by the explicit numerator term, so n=8 is the exact homogeneous onset with this negative-index convention. As a conventional order-nine recurrence, it holds from n=9 onward. Mixed leaves satisfy M_n=M_(n-2)+M_(n-6) for every n>4, with the only nonzero forcing term M_4=1.

An extra exact cancellation occurs in even order. For r>=5,

    a_r-2a_(r-2)=a_(r-7),

by two applications of the recurrence (equivalently t^7+2t^2-1 is divisible by1-t-t^3). Thus

    L_(2r)=(a_(2r)-a_(r-7))/2,   r>=5.

In particular, total leaf counts are exactly half the rigid-word count at n=10 and12, and strictly below half at every even n>=14 and every odd n>=3. This sharpens the earlier asymptotic L_n/a_n->1/2 without changing it.

The coefficients n=0,...,30 are given in leaf_series_coefficients.json and independently checked by direct convolution and both parity formulas. These enumerate leaf vertices, not bridge edges; the two leaves at n=4 are the endpoints of the same edge.

The total-leaf denominator above is reduced, so the minimal constant-coefficient eventual recurrence over the rationals has order nine. Indeed, after removing the harmless factor z^4 from the numerator, write C=2-z-z^3. For P=1-z-z^3 we have C=P+1, hence gcd(C,P)=1. Modulo C, Q=1-z^2-z^6 reduces to -3+4z-2z^2. Modulo the corresponding quadratic 2z^2-4z+3, C reduces to5-7z/2. Its only possible common root is z=10/7, where the quadratic equals67/49, not zero. Also Q(0)=P(0)=1. Thus gcd(2z^4-z^5-z^7,PQ)=1. Likewise the mixed-leaf recurrence has minimal order six because gcd(z^4,1-z^2-z^6)=1.
