# An exact closed trace carrier for the quoted quadratic units

**PROVED elementary**, local results Q1-Q3; **FINITE-EXACT** controls.
Same-agent audit. This specializes and extends the incoming
[Session93 golden trace-clock theorem](../session93_reversal_lucas_clocks_20261029/arithmetic/TRACE_CLOCK_CARRIERS.md)
to the actual units in the user's formulas. Smith normal form, quadratic
recurrences and trace-pairing discriminants are classical mechanisms.
No external-priority claim is made.

## Q1. Specify the source, clocks and conjugation

For an integer a>=1 let

    O = Z[u],       u²-a u-1=0,       Delta=a²+4,
    bar(u)=a-u,     Tr(r+s u)=2r+a s,
    N(r+s u)=r²+a r s-s².

Work on all m² states of O/mO, m>=1. Fix positive clock exponents
e_1,...,e_t and write b_i=u^e_i=A_i+B_i u. The allowed operations are
multiplication by any b_i or its inverse, and conjugation. Observe Tr(x).
Every operation acts on the SAME source state. Arbitrarily independently
chosen states would not have these joint observations.

**Theorem.** The exact smallest behavioral carrier is the vector

    C(x)=(Tr x, Tr(b_1 x),...,Tr(b_t x)) mod m.

Put g=gcd(B_1,...,B_t), delta=gcd(2,a). Every fiber has size

    d = gcd(m,delta) gcd(m,Delta*g/delta),

and there are exactly m²/d carrier states. This quotient is closed under
all the stated clocks, their inverses, and conjugation.

**Proof of completeness and closure.** The Z-module spanned by 1 and the
b_i is Z+g Z u, by Bezout on their B_i coefficients. It is a subring:
(g u)²=g²+a g²u belongs to it. Since N(b_i)=(-1)^e_i,
b_i^-1=(-1)^e_i(Tr(b_i)-b_i) belongs to the same subring. So every future
clock word is in the span already observed. Conversely time zero and each
single clock are allowed, forcing every listed observation.

Conjugation is also closed: Tr(bar x)=Tr x and

    Tr(b_i bar x)=Tr(b_i) Tr(x)-Tr(b_i x).

Clock transitions are integral linear combinations of the listed entries;
they can be obtained by the same Bezout expression and multiplication in
Z+gZu. No inverse of a nonunit modulo m is needed to update the quotient.

The observation matrix has first row (2,a), and subsequent rows

    (2A_i+aB_i, aA_i+(a²+2)B_i).

Its entry gcd is delta. Its two-by-two minors with the first row are
Delta B_i. A minor between two other rows is, up to sign,
Delta(A_i B_j-A_j B_i). Thus the gcd of all minors is Delta g.
The Smith invariant factors are (delta,Delta*g/delta). Over Z/m their
kernel has the displayed cardinality d, proving the exact minimum. QED.

This counts all ring states, not only units, a particular Pell orbit,
or norm-one states. Those restrictions can change fiber sizes.

## Q2. A single clock has a two-observation recurrence and a retained norm

For b=A+B u and y_k=Tr(b^k x),

    y_(k+2)=Tr(b)y_(k+1)-N(b)y_k.

Furthermore

    y_1²-Tr(b)y_0 y_1+N(b)y_0² = -Delta B² N(x).

Both identities follow directly from b²-Tr(b)b+N(b)=0 and the conjugate
pair. These are exact over Z and modulo every m. Recovering N(x) by division
requires Delta B² to be invertible; the displayed identity itself does not.

The primitive clock b=u has Smith factors (delta,Delta/delta). Decimation
to a large exponent introduces the additional divisor B. Thus a beautiful
scalar recurrence can coexist with a much larger lost phase. This is a
computable loss, rather than a generic warning about information.

## Q3. The 22,37,58 units and their exact missing coordinates

The [Ramanujan formulas](RAMANUJAN.md) use the following real units and even
exponents. Their imaginary CM fields are separate objects.

| n | u | a | exponent e | (A,B) for u^e | Smith factors |
|---|---|---:|---:|---|---|
|22|1+sqrt(2)|2|12|(5741,13860)|(2,55440)|
|37|6+sqrt(37)|12|6|(21169,255780)|(2,18927720)|
|58|(5+sqrt(29))/2|5|12|(13741001,71351280)|(1,2069187120)|

At37 the real order Z[6+sqrt(37)]=Z[sqrt(37)] has discriminant148 and
conductor2 in the maximal real order. Replacing148 by37 would change the
even-modulus trace kernel. At22 and58 the displayed orders are maximal.

The determinant factorizations identify concrete appearances of the small
primes, rather than matching numerals after the fact:

    a=2, e=12: Delta B = 2^5 * 3^2 * 5 * 7 * 11;
    a=12,e=6:  Delta B = 2^4 * 3^2 * 5 * 7^2 * 29 * 37;
    a=5, e=12: Delta B = 2^4 * 3^4 * 5 * 7 * 11^2 * 13 * 29.

For example, in the58 unit ring modulo13, u^12=1. The entire sampled
trace history is constant and has only13 possible values on169 states.
Adding the primitive clock Tr(u x) recovers all169, because29 is a unit
modulo13. In the22 unit ring modulo11, u^12=-1: the sampled history only
alternates the initial trace. At the ramified modulus29 in the58 example,
even the primitive trace pairing has a29-element kernel. That loss cannot
be repaired merely by adding another power clock.

Modulo8, the22 sampled carrier has four cells on64 source states; each
fiber has16 states. A universal 'one missing bit' description is therefore
false for this ring-state menu. Over the positive real Fricke branch the
two-valued choice in [RAMANUJAN.md](RAMANUJAN.md) is correct; the source
space and operation are different.

This gives an explicit arithmetic comparison for the tournament objective.
The determinant ideal computes exactly which joint observations recover
phase here. The tournament counterpart is an actual intersection of endpoint
permutation images and, for stronger operations, simultaneous transporters.
The two mechanisms are not identified as one group or one state space.

Reproduction: [unit_clocks.py](unit_clocks.py), [unit_clocks.json](unit_clocks.json).
The full finite universe is a=1,2,5,12; m=1..30; exponent menus
(1),(2),(6),(12),(4,6),(6,12),(5,7), plus named controls at8,13,29,37.
Every ring state is checked; kernel and image counts are measured directly,
and both algebraic identities are verified before reducing modulo m.
