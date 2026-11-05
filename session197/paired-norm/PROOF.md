# Paired modules, quadratic norms and a parity obstruction for 7p

Session 197. This is a separate consequence of the reviewed three-module identity. It gives explicit self-converse strong graphs, an all-order 7p obstruction for one parity, and a 25-vertex example for the other. It does not claim all prime counts, all 7p counts, or infinitely many prime quotients.

## 1. The physical paired construction

Use the nine-vertex quotient Q with outgoing rows

    [30,64,10,66,174,335,21,367,95]

and anti-involution

    tau=[0,7,8,5,6,3,4,1,2].

For a nonempty tournament U and integer m>=1, replace vertex0 by Tr_m, vertex2 by U and vertex8 by U^op, its converse. All other quotient vertices remain singletons. Call the actual result P_m(U).

Q is strong, and nonempty substitution preserves strongness. The anti-involution reverses the Tr_m chain, exchanges each vertex of U with its corresponding vertex in U^op, and exchanges the other paired quotient vertices. It reverses every arc. Thus P_m(U) is self-converse and strong, with order

    6+m+2|U|.

Write h=H(U) and b=B_U(2), where B is the labelled-bin path-cover polynomial. Reversing each bin path shows that U and U^op have equal B profiles. The reviewed exact three-module formula therefore gives

    H(P_m(U)) = h^2 + 18 b^2 3^m.                           (1)

This is an actual graph identity for arbitrary nonempty U, not a prescription of arbitrary h,b. The bin-label swap acts freely on two-bin covers of a nonempty tournament, so b is a positive even integer. Also h is positive and odd.

## 2. Quadratic norm forms and elementary congruences

For even m=2r, r>=1, put y=3^(r+1)b. For odd m=2r+1, r>=0, use the same y. Equation(1) becomes

    even m: H=h^2+2y^2;
    odd m:  H=h^2+6y^2.                                    (2)

These are the elementary norms of h+y sqrt(-2) and h+y sqrt(-6), respectively. No unique factorization or class-number theorem is used.

Because b is even and h is odd, every count in (1) is one modulo eight. Modulo three it equals h^2. Consequently any prime output is one modulo24; an output divisible by three is larger than three and hence composite. The minimum possible output is217.

If H=7p with p prime, then p is seven modulo24. Indeed p>3 because H>=217, and reducing H modulo eight and three gives respectively p=7 modulo8 and p=1 modulo3. This is a necessary condition only.

## 3. Even root multiplicity can never produce 7p

**Theorem.** For every nonempty U and every even m>=2, H(P_m(U)) is not 7p for any prime p.

If seven divides h^2+2y^2, then h^2=-2y^2 modulo seven. The nonzero quadratic residues modulo seven are1,2,4, while -2 is5. Hence y=0 and h=0 modulo seven. It follows that49 divides H.

If H were7p with p prime, this would force p=7 and H=49. But y>=3^2*2=18, so H>=1+2*18^2=649. This contradiction proves the theorem. It includes the possible quotient primes3 and7 rather than silently excluding them by a congruence.

The argument is the inert-prime norm obstruction at seven, expressed entirely by elementary modular arithmetic. It applies to arbitrary paired input modules, not only transitive ones.

## 4. Odd multiplicity gives an explicit quotient norm certificate

For odd m, suppose seven divides H=h^2+6y^2. Then h=+y or -y modulo seven. In the first case,

    H/7 = ((h+6y)/7)^2 + 6((y-h)/7)^2.

In the second case,

    H/7 = ((h-6y)/7)^2 + 6((h+y)/7)^2.                      (3)

The displayed coordinates are integers by the respective congruence, and the identities follow by expansion. They are odd because h is odd, y is even and seven is odd. Thus every prime quotient in the odd-m 7p channel comes with an explicit representation p=x^2+6z^2 with x,z odd. This is not a converse representation theorem or a claim that every such prime is a physical output.

For the original source m=1,U=K1, one has h=1,y=6 and p=31=5^2+6*1^2.

## 5. A cyclotomic subfamily and a new explicit 7p graph

When U=Tr_s, h=1 and b=2^s. Thus

    H(P_m(Tr_s))=1+2^(2s+1)3^(m+2).

Taking any odd integer q>=3, set s=(q-1)/2 and m=q-2. This gives a self-converse strong graph of order2q+3 with

    H=6^q+1.

Its quotient by seven is the elementary alternating sum

    (6^q+1)/7 = 6^(q-1)-6^(q-2)+...-6+1.

If this quotient is prime, q must be prime. For composite odd q=ab with a,b>=3, factor

    6^q+1=(6^a+1)[(6^a)^(b-1)-(6^a)^(b-2)+...-(6^a)+1].

Both factors exceed seven, which is impossible for7 times a prime. A prime q is not sufficient: q=5 gives7*11*101, and q=7 gives7^2*29*197.

For q=11, choose m=9 and U=Tr_5. The resulting25-vertex graph has

    H=6^11+1=362797057=7*51828151,

and51828151 is prime. The verifier tests every integer divisor from2 through7199, the exact integer square-root bound. With h=1,y=7776, (3) gives the additional exact certificate

    51828151=6665^2+6*1111^2.

`witness.json` contains every outgoing adjacency mask, the complete physical anti-involution, the block labels and arithmetic certificate. Its Hamilton count follows from the proved all-order identity (1). We do not pretend that a25-vertex endpoint DP was run. Independent small physical counts control the construction, while direct checks verify all25-vertex orientations, the anti-involution and strongness.

## 6. What the paired unit-offset specialization can and cannot do

For fixed U, (1) is an exact ternary family with amplitude18 b^2 and offset h^2. It has offset one exactly when H(U)=1, equivalently when U is transitive. In that unit-offset subclass,

    amplitude=18*4^s=72*4^(s-1), s>=1.

Thus arbitrary paired modules do not automatically create arbitrary unit-offset amplitudes. This is a scope statement for this specific construction, consistent with the separate sharp amplitude72 theorem and its unique equality core.

## 7. Sources and reproducibility

The substantive inherited input is `../../session197-odd-cycle-path-polynomials/THREE_MODULE_CLOCK.md`, frozen SHA256 `25ac0883cb93cd55b1109160bdb253794413e3aa3b56682bccf28f1dd6542941`. Its independent audit is `../../session197-bridges/THREE_MODULE_REVIEW.md`, SHA256 `59f31c66a2d23f61cf2d2f9d9246ca2358a183ff14c12e7070c2dfb6ff759d1c`. That packet proves H(Q[T,U,V])=18 B_T(3)B_U(2)B_V(2)+H(T)H(U)H(V). The present paired specialization, parity-sensitive norm obstruction and quotient certificates are consequences proved here. No literature-priority claim is made.

Run `verify.py` normally and under `python -O`; stdout receipts must agree. It is standalone standard-library code. It directly counts sixteen small physical graphs by Hamilton subset DP, checks every lifted anti-involution, validates2,048 integer norm controls, builds the explicit25-vertex graph and tests its primality certificate. Every check uses an explicit exception, not a removable assertion. The large graph count is identified in the receipt as an evaluation of the proved three-module identity.
