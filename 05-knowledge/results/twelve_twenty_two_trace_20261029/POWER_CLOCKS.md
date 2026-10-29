# Integral power clocks and the first mixed-observation boundary

**Status: PROVED scoped**, except explicitly finite tables and cited inputs.
The algebraic lemmas are elementary; no priority claim is made. The
Session104 counter and October28 quadratic observer are inherited.

## P1. The unit-clock gcd law

Let A be a commutative ring free of finite rank over Z and let u be a unit.
For positive a,b,

\[
\mathbb Z[u^a,u^b]=\mathbb Z[u^{\gcd(a,b)}].
\]

Multiplication by any unit v has a monic integer characteristic polynomial
with constant term plus or minus1. Cayley--Hamilton expresses v inverse
as an integer polynomial in v. Bezout's identity therefore puts
u^gcd(a,b) in the left-hand ring. The reverse inclusion is immediate.
This is a ring identity, not a positive-word reachability assertion.

For A=Z[u]/(u^3-2u^2-u+1), direct reduction gives

\[
u=u^4-5u^2+2,\qquad
d=3d^2-d^4-1\quad(d=u-1).
\]

Hence Z[u^2]=Z[u] and Z[d^2]=Z[d]. Applied to the matrix C of the README,
this proves integral recovery from the squared clock. It does **not** say
that the map sending an arbitrary state to its square is linear or an
automorphism; nor is u mapped to u^2 an endomorphism of this fixed algebra.

## P2. A rank-three sum lemma

Let A be free of rank3 over Z. If unital subalgebras B1,B2 generate A as
a ring, then their **additive** sum B1+B2 is all of A.

To prove this, reduce modulo any prime p. If the images have a sum of
dimension at most2, either one is contained in the other or both have
dimension1. In either case their sum is already a proper subalgebra,
contradicting generation of A/pA. Thus the sum surjects modulo every
prime. The finitely generated quotient A/(B1+B2) is zero: a nonzero free
part or nonzero finite part would survive modulo some prime.

The argument permits ramification and nilpotents. It does not require a
field or a semisimple algebra. Cayley--Hamilton in rank3 also says that
Z[v] is additively spanned by1,v,v^2 for every v in A. Consequently, if
Z[u^gcd(a,b)]=A, the five elements

\[
1,\ u^a,\ u^{2a},\ u^b,\ u^{2b}
\]

span A over Z. For an arbitrary common-source linear readout ell, their
observations recover every ell(u^j x) by integer linear combinations.
This recovers the primitive observer, which need not recover x itself.

## P3. Two observers on the same cubic clock

For the actual safe-word counter use column coordinates (B0,A,B1),
initial state (1,0,0), and ell=(1,1,1). The primitive observation matrix is

\[
O_{\rm count}=\begin{pmatrix}1&1&1\\3&2&2\\7&4&5\end{pmatrix},
\qquad\det O_{\rm count}=-1.
\]

Thus it reconstructs an arbitrary integer state modulo every m. A
literal independent enumeration of the six-state, five-letter grammar
agrees with ell C^n x for n=0 through7. The first counts are
1,3,7,16,36,81,182,409. Their interpretation and the distinction between
safe counting and legality are inherited from
[Session104](../session104_descent_exterior_20261029/README.md).

For comparison use the arithmetic trace of A, with source basis1,u,u^2.
Here Tr(1)=3, Tr(u)=2, Tr(u^2)=6, and the primitive trace Gram matrix has
Smith factors1,7,7 and determinant49. Its full clock-observation carrier
modulo m has cardinality

\[
m^3/\gcd(m,7)^2.
\]

At7, f(u)=(u-3)^3. Write u=3+epsilon in
F7[epsilon]/(epsilon^3). Multiplication trace is three times the constant
coefficient, so the two-dimensional ideal (epsilon) is invisible even
after **every** primitive clock power. Frobenius/Galois transport cannot
make a trace see this ideal. This loss is a property of the readout,
absent from O_count, rather than a property of C alone.

The exact three-reading decimations have the following Smith data.
`Index` is [A:Z[u^e]], computed as a lattice index.

| e | Index | Counting observer | Arithmetic trace observer |
|---|---:|---|---|
|1|1|(1,1,1)|(1,7,7)|
|2|1|(1,1,1)|(1,7,7)|
|12|2739659 = 13^3*29*43|(1,13,210743)|(1,91,1475201)|
|22|3282660883657 = 881*5347*696851|(1,1,3282660883657)|(1,7,22978626185599)|

The last two indices are coprime. This is an alternative explicit proof
that the two full separate-clock banks span A. P1--P2 prove the same
fact structurally: gcd(12,22)=2 and Z[u^2]=A.

## P4. Four readings by44, three by48

**PROVED finite optimization, with exact integer certificates.** In the
prescribed bank (0,12,24,22,44), no three power vectors span A over Z.
There are three successful four-element subbanks; one is (0,12,22,44).
The saved HNF transformation in `clocks.json` reconstructs1,u,u^2 by
integer combinations of these four vectors. Thus the actual counting
observer is injective over every modulus. The trace observer recovers
exactly its primitive carrier, with its existing7-primary kernel.

Now allow all times in the semigroup {12a+22b:a,b>=0}. The nontrivial
successive lattice indices up to the first complete horizon are

| Horizon |22|24|34|36|44|
|---|---:|---:|---:|---:|---:|
|Index of all available power vectors|542607|13|13|13|1|

The lower bound has a direct reason. Modulo13, u^12=1. Before time44
at most one22-step is possible; every available vector lies in the span
of1,u^22, of dimension2. At time44 the missing squared22-step becomes
available. The integer HNF proves that no other prime obstruction remains.

If only three readings are allowed, the first complete bank is
(44,46,48). The determinant of its power vectors is -1. Multiplication
by the unit u^44 reduces this to the basis (1,u^2,u^4), whose determinant
is1 up to sign. Exhaustion of all triples of semigroup times at most48
shows that this is the sole unit-determinant triple in that range. This
finite exhaustion certifies the earliest three-reading horizon; it is
not a claim that four observations remain necessary with unbounded time.
Two scalar linear readings cannot reconstruct a rank-three source over Q.

**Operation guard.** The recovered integer combinations can have negative
coefficients. They do not assert a positive-word implementation of a
single C-step. For example modulo3, C has order26 and C^2 has order13;
all12/22 words stay in the even subgroup, although their observations
linearly reconstruct the full counting state. State identifiability,
legal trajectory reachability and physical port realization are distinct.

## P5. An all-r13-adic precision recurrence

Put u^12=1+13A0. Exact reduction gives that
(1,A0,A0^2) has determinant of magnitude29*43, a13-adic unit. For r>=0,
put s=r+1 and write

\[
u^{12\cdot13^r}=1+13^s A_r.
\]

The binomial theorem gives A_(r+1)=A_r modulo13. Thus each
(1,A_r,A_r^2) is a basis over Z13. The three power rows1,b,b^2, with
b=1+13^s A_r, reduce integrally to

\[
1,\quad13^s A_r,\quad13^{2s}A_r^2.
\]

The **13-primary Smith exponents are exactly (0,s,2s)**. Both the
counting observer and the arithmetic trace pairing are invertible at13,
so the law holds for both readouts. At modulus13^t the logarithm to base13
of the visible state count is

\[
E_r(t)=t+(t-r-1)_+ +(t-2r-2)_+,
\]

and, for t>=0,

\[
E_r(t+1)-E_r(t)=1+\mathbf1_{t\ge r+1}+\mathbf1_{t\ge2r+2}.
\]

This is an all-r proof, checked at r=0 through8 using modular arithmetic
at precision2r+3. For every r, gcd(12*13^r,22)=2. P1--P2 therefore prove
that adding readings at22 and44 restores the **whole primitive observer
over Z**, not merely its13-primary part. Increasing precision memory
and inexpensive complementary recovery coexist in this particular algebra.

## P6. Ramanujan22: the quadratic hostile control

The inherited quadratic trace analysis concerns v=1+sqrt(2), satisfying
v^2=2v+1. Ramanujan's g22=sqrt(v) yields the integer trace approximation
64 Tr(v^12)+24=2508952; see [the original paper, section5](https://ramanujan.sirinudi.org/Volumes/published/ram06.pdf)
and [the prior trace carrier](../cursor_cm_phase_20261028/TRACE_CARRIER.md).

Unlike the cubic counter,

\[
\mathbb Z[v^{12},v^{22}]=\mathbb Z[v^2]
=\mathbb Z+2\mathbb Z v\ne\mathbb Z[v].
\]

The index is2. The primitive trace observer has Smith factors(2,4),
whereas all even-clock observations have factors(2,8). At modulus8,
z=v-1=sqrt(2) has zero trace against every even power, but
Tr(vz)=4. The two clocks retain an **additional dyadic loss relative
to the primitive observer**; this is not a claim that the complete source
requires just one bit. At11, v^12=-1 in F121, so the twelfth-clock trace
has only11 states out of121; adding22 restores the primitive observer
there, since the remaining index2 is invertible.

More generally, for a quadratic norm-minus-one unit with v^2=a v+1,
the squared-power order has index |a| in the declared order Z[v], for
a nonzero integer a. The cubic index1 recovery and the Pell index2 loss
are instances of an integer-basis question, not consequences of equal
field degrees or matching periods.

## P7. Rank four: separate clocks miss a mixed amplitude

Use the actual two quadratic factors of Ramanujan's g130:
a=2+sqrt(5), b=(3+sqrt(13))/2, so g130^2=ab. In
R=Z[a,b], use the rank-four basis(1,b,a,ab). The relations are
a^2=4a+1 and b^2=3b+1. The two separate clock algebras generate R as a
ring but their additive sum has rank3. Thus P2 fails in rank4.

Let z=(2a-4)(2b-3)=2sqrt(65). Field trace factors across the two quadratic
extensions. For **every j>=0**,

\[
\operatorname{Tr}(z a^j)=\operatorname{Tr}(z b^j)=0,
\qquad \operatorname{Tr}(zab)=260.
\]

For a pure a-reading the second factor has trace zero, and conversely.
The mixed reading uses Tr((2a-4)a)=20 and Tr((2b-3)b)=13. Its product260
is nonzero. In this basis the trace Gram matrix is

\[
\begin{pmatrix}
4&6&8&12\\6&22&12&44\\8&12&36&54\\12&44&54&198
\end{pmatrix},
\]

with Smith factors(2,2,130,130). Modulo3 its first three rows have rank3
and all four have rank4: a mixed reading raises the visible state count
from27 to81. The source, operation and witness are explicit; the entries
12,22,44 are not an identification with any other objects of those sizes.

The inherited g130 warning concerned neglected conjugate terms in an
analytic near-integer estimate. Here the same biquadratic source gives a
different, exact observation experiment. Both require joint data, but the
predicate and observation language must be stated separately. This example
does not prove a joint-port theorem for arbitrary tournaments.
