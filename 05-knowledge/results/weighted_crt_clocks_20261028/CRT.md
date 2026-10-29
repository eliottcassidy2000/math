# Odd Lucas readers beyond76: relative phase is an actual coordinate

**PROVED** C1-C2; **PROVED** finite-state minimization principle C3;
**FINITE-EXACT** census applications. [Script](check_crt.py), [JSON](check_crt.json).

## C1. A uniform odd-period ring chart

Let n>=3 be odd, O=Z[phi], phi^2=phi+1, and R_n=O/(phi^n-1).
Write F_j,L_j for Fibonacci and Lucas numbers. Set

    A=F_(n-1)-1, B=F_n,  phi^n-1=A+B phi,
    M=L_n                 if 3 does not divide n,
    M=L_n/4               if 3 divides n.

Then M is odd, gcd(B,M)=1, and with r=-A/B modulo M,

    R_n = Z/M                              (3 does not divide n),
    R_n = F4 x Z/M                         (3 divides n),           (C1)
    phi -> r,                  respectively phi -> (omega,r).

The ring isomorphisms include multiplicative structure, not just cardinality.
Factor M into its prime powers to obtain the full CRT chart; a prime-power
factor is not in general a field. The convention Z/1 is harmless at n=3.

Proof. Multiplication by A+B phi has matrix [[A,B],[B,A+B]], determinant
Norm(phi^n-1)=-L_n, so |R_n|=L_n. A common odd prime divisor of A and B
would give phi^n=1 in O/p but Norm(phi^n)=-1, impossible in odd characteristic.
Modulo2, phi has order3 in F4. If3 divides odd n, B is2 mod4 (Fibonacci
period modulo4 is6), hence gcd(A,B)=2; otherwise it is1. Dividing A,B by2
in the first case gives odd norm modulo2, since the second coordinate is odd;
thus v_2(L_n)=2. The identity L_n^2-5F_n^2=-4 shows gcd(B,M)=1.

The relation A+B phi=0 therefore forces phi=r in the odd quotient; its norm
gives r^2-r-1=0 there. Evaluation onto Z/M is surjective. When3 divides n,
evaluation also maps onto F4; the coprime integer moduli give surjectivity
onto their product. Equality of cardinalities proves(C1).

The inherited all-odd reader sends cyclic no-11 words to sum_i w_i phi^i;
its bijectivity is credited to
[THM-4592, current all-period addendum](../../../01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md) and the
[76-state application](../lucas_prime_reflection_towers_20261028/COMPOSITE_76.md).
The present script independently checks bijectivity for every state in its
declared seven-length universe. It does not treat that census as the proof
of the inherited all-size reader theorem.

Rotation tau multiplies by phi. Word reflection R reverses exponents.
On F4 it is exactly Frobenius a->a^2. On the odd factors it generally
requires the full word decoder. In particular phi^-1=phi-1 does not satisfy
the same quadratic polynomial: f(phi^-1)=-2phi^-1. Thus reflection is not
silently promoted to a ring automorphism. Arithmetic translation T:x->x+1
is evaluated in the ring and then decoded back to a legal word.

## C2. Compatible phases, including zeros and shorter periods

Choose one representative in each orbit of multiplication by the local
root in every prime-power factor (and in F4, when present). For a local
state, let d_i be its orbit period and a_i in Z/d_i its phase. Include the
zero state with period1. Since phi^n=1, each d_i divides n.

For a fixed tuple of local orbit types, synchronized rotation adds1 to all
phases. Put D=lcm_i d_i. The number of global rotation orbits is exactly

    product_i d_i / D.                                  (C2)

Every diagonal orbit has length D, proving the count. A tuple of local
phases can simultaneously be moved to the chosen representatives iff

    a_i = a_j mod gcd(d_i,d_j) for every i,j.            (C3)

This is precisely generalized CRT compatibility. In general retain the
relative-phase class in (product Z/d_i)/diagonal(Z/D), followed by one
global phase. Their sizes multiply back to product d_i: this is a faithful
chart organization, not free compression of the full reader.

Choose a seed per global rotation orbit. Since R tau R=tau^-1, reflection
sends a chart phase i to j-i in the reflected seed's chart. The chart type
and seed shift j are indispensable. Arithmetic addition generally changes
chart type and relative phase; it must be evaluated and recharted.

| n | L_n | CRT factors | independent local orbit types | global rotation orbits |
|---|---:|---|---:|---:|
|9|76|4,19|6|10|
|11|199|199|19|19|
|13|521|521|41|41|
|15|1364|4,11,31|18|94|
|21|24476|4,29,211|110|1170|
|23|64079|139,461|147|2787|
|25|167761|11,101,151|105|6713|

At n=15 the local nonzero periods are3,5,15. A fully nonzero tuple therefore
needs15 relative-phase classes for each tuple of local orbit types. The
global periods have counts `{1:1,3:1,5:2,15:90}`. The excess94-18=76 counts
omitted global orbits in this specific projection; it is not an identification
of that set with R_9. The agreement of integers has a stated counting source.

The independent necklace formula agrees:

    number of rotation orbits = (1/n) sum_(d|n) EulerPhi(d) L_(n/d).

It follows by counting words fixed by rotation, including period1. Prime
Lucas values collapse(C1) to a field; the recent
[prime-degree reflection theorem](../session93_reversal_lucas_clocks_20261029/lucas_groups/LUCAS_GAP_THEOREM.md)
then applies under its explicit primality hypothesis. Its A/S parity rule
depends on n mod12; it cannot be carried to these composite degrees by
forgetting their invariant CRT factors.

## C3. Future attachments can destroy the proposed quotient

For a finite state set X, observation o, and named operations f_1,...,f_s,
define P_0 by equality of o and refine using

    key_(t+1)(x)=(P_t(x), P_t(f_1 x),...,P_t(f_s x)).     (C4)

Induction shows P_t agrees exactly when all operation words of length<=t
give the same observation. A stable partition is the coarsest congruence
refining the observation, hence the smallest deterministic future carrier.
This proof applies to noninvertible operations too. Our operations R,T
are permutations, but no unobserved rotation generator is added to the menu.

Apply(C4) to one CRT projection, allowing only reflection R and addition T
of1. Every listed odd projection gives the full reader partition:

| n | observed modulus | smallest future carrier | separating depth |
|---|---|---:|---|
|9|19|76|2|
|15|11;31|1364|6;4|
|21|29;211|24476|12;7|
|23|139;461|64079|12;10|
|25|11;101;151|167761|34;18;18|

Depth is the first refinement whose classes are all singletons, not the
length of one universal distinguishing word. The script exhausts every
reader state and keeps the full refinement counts in JSON. Projection to
F4, when present, remains exactly a four-state carrier: both Frobenius and
addition1 are already defined there. At prime degrees199 and521 the only
factor observation is identity, so faithfulness there is a trivial control.

**OPEN:** every nontrivial odd CRT-factor observation forces the full reader
for all odd periods under this operation menu. The finite positive results
do not establish it, especially at odd prime powers not exhausted here.
Nor is the minimal quotient of one specified observer a universal minimal
carrier for every arithmetic or tournament operation.
