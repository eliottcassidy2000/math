# Golden values with an exact binary source: three registers and guarded receipts

2026-10-07. **Status: PROVED** elementary integral-ring, source-guard and
receipt-composition statements; **FINITE-EXACT** for the declared controls.
The result supplies a faithful arithmetic interface. It adds no universal
Collatz coverage and does not infer a prescribed source's fate from density.

## 1. Inheritance, types, and the useful new coordinate

The anchor is a supplied integer's actual guard. The niche is a joint
golden/binary digit reader. The wildcard is the specific list
105, 223, 233, **332**, 425. The concept board is source polynomial / golden
value / binary value / carry / clock / authenticated common future.

The closest mechanisms are the exact pair reader and unique carry in
[golden_digit_carry_20261003](golden_digit_carry_20261003.md), the distinct
integer-modulus and principal-ideal clocks in
[golden_prime_clocks_20261003](golden_prime_clocks_20261003.md), and the
ordered guarded interface in
[golden_primes_carries_route_compiler_20261003](golden_primes_carries_route_compiler_20261003.md).
The hostile is the same word read under different radices. The repaired
near miss is treating a golden shift as division by two. The least-used
sidecar is the carry **evaluated at the binary radix**, rather than its
whole polynomial or just its golden value.

[THM-4528, Collatz parity in base phi](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md)
reads an **orbit itinerary**, not the ordinary binary expansion of its
source. Its signed boundary repairs are retained in
[collatz_golden_carriers_20261004](collatz_golden_carriers_20261004.md).
Our input here is instead a finite polynomial P(z) whose coefficients are
the actual source's binary digits. These types must not be identified.
Likewise the literal positional base-phi expansion of the decimal integer n
has value n, unlike the golden reading of its ordinary binary digits.

The current [THM-4590, equidistributed Collatz classes](../../01-canon/theorems/THM-4590-collatz-classes-are-equidistributed-residues-windows-slowly-varying-density.md)
concerns arrangements of saturated classes and explicitly does not prove
individual coverage. None of its statements is needed below. Its referenced
`golden_collatz_resonance_20261007.md` was absent from the synced checkout;
no assertion is inherited from that missing document.

## 2. An integral CRT carrier, with three exact coordinates

Write f(z)=z^2-z-1 and O=Z[phi], where phi^2=phi+1. The two evaluations are

    P(phi)=A+B phi,       P(2)=N.

Because f(2)=1, the ideals (f) and (z-2) are comaximal in Z[z]. Indeed
f(z)-1=(z-2)(z+1). Their intersection is their product, so the ordinary
Chinese remainder argument, with no division of integers, gives

    Z[z]/((z-2)f(z))  ~=  O x Z.                         (1)

The explicit inverse sends (A+B phi,N) to

    A+Bz+D f(z),       D=N-A-2B.                         (2)

Thus (A,B,D) is a three-integer carrier; N=A+2B+D is recovered exactly.
It retains both evaluations, not the whole input polynomial or its parse.
Its kernel is exactly ((z-2)f(z)). The polynomial f represents the central
idempotent (0,1) in O x Z; 1-f represents (1,0).

Appending an integer digit d by P -> zP+d gives

    (A,B,D) -> (B+d, A+B, 2D+B).                         (3)

For binary d this is a streaming reader of the ordinary binary source.
With d=0 the matrix is

    J = [[0,1,0],[1,1,0],[0,1,2]],
    det J=-2,   char_J(z)=(z-2)(z^2-z-1).

Its three eigenvalues are 2, phi, and -phi^(-1). In particular the discarded
source mode grows like 2^k: it is not a bounded carry. The minimal polynomial
has degree three. More intrinsically, any integral linear quotient retaining
both evaluations has kernel contained in the kernel of (1), hence has rank
at least three. This is minimality in the stated linear, arbitrary-integer-
coefficient model; it is not a lower bound on every nonlinear coding scheme.

Addition is componentwise. If s=(A,B,D), t=(C,E,F), their product is

    A'=AC+BE,
    B'=AE+BC+BE,
    D'=BE+D(C+2E)+F(A+2B)+DF.                            (4)

This follows by multiplying independently in O and Z. It proves
associativity and preserves both input values. The term BE is required:
the state of z is (0,1,0), while that of z^2 is (1,1,1), not (1,1,0).
The latter represents z+1, with binary value 3 instead of 4.

This is multiplication of **digit polynomials**, which can have overlapping
digits. It is not the assertion that golden evaluation of the canonical
binary expansion is multiplicative. For example the square of P_3=1+z has
state (2,3,1), whereas P_9=1+z^3 has state (2,2,3). Both binary values are 9,
but their golden values differ by phi. Binary normalization has its own
unique carry, a multiple of z-2; here
(1+z)^2-(1+z^3)=-(z-2)z(z+1). It preserves the Z component while changing
the O component. Golden normalization by multiples of f does the reverse.
The joint carrier makes both losses explicit; a chosen normal form still
needs its corresponding carry certificate.

Given all three coordinates modulo 2^K, every source predicate modulo 2^K
is recoverable. The same statement holds for any modulus. This preserves
an already proved guard; it does not prove the guard eventually applies.

## 3. The requested 223/233 pair detects the missing coordinate exactly

Let P_n(z) be the ordinary binary polynomial of n. Directly,

    P_233(z)-P_223(z)=(z^2-z-1)(z^3+z).                 (5)

Consequently both golden values equal 18+28 phi, but their D registers
are respectively 159 and 149. At z=2 the missing carry is 2^3+2=10.
Their actual first odd Collatz valuations are different:

    v2(3*223+1)=1,       v2(3*233+1)=2.

This extends to an infinite exact family. For every t>=0, the lower eight
bits of 256t+223 and 256t+233 remain those in (5), while their common higher
polynomial is z^8 P_t(z). Hence they have identical golden values and
their D registers differ by 10. The former's next odd step grows; the
latter's next odd step descends. No function of the golden value alone can
classify this first-step guard, even on this one explicit family.

The smaller collision 100_phi=011_phi identifies binary sources 4 and 3.
It rules out a source-residue decoder from golden values modulo **every**
integer m>1, since their binary difference is 1. Among odd sources, 7 and 9
also have equal golden values and different next valuations.

The precision bill is real. For K>=2, the two odd binary sources

    n=2^(K+1)+1,       m=3*2^(K-1)+1

have equal golden values, but their D values differ by 2^(K-1). Thus
retaining D only modulo 2^(K-1) can lose a K-bit source guard even when
the golden value is known exactly. For K=1, sources 3 and 4 give the parity
hostile. This claim is about arbitrary source residues, not an assertion
that every individual guard needs all K bits.

The ordinary binary reading is not a primality test either: 3 and 4, or
7 and 9, supply prime/composite collisions. This does not claim that the
literal positional base-phi values 3 and 4 are equal; they are not.

## 4. Two clocks, not one: 105, mod 18, and mod 19

Let pi_phi(m) be the order of multiplication by phi on O/(m). If m is odd,
the exact period of the joint zero-padding operator is

    lcm(pi_phi(m), ord_m(2)).                            (6)

This is immediate in the product coordinates (A+B phi,N): zero padding
multiplies the two factors by phi and 2 respectively. The state of 1 attains
the full period, so (6) is exact, not just a period bound.

| Modulus | Golden period | Binary period | Joint period |
|---:|---:|---:|---:|
| 3 | 8 | 2 | 8 |
| 5 | 20 | 4 | 20 |
| 7 | 16 | 3 | 48 |
| 11 | 10 | 10 | 10 |
| 19 | 18 | 18 | 18 |
| 105 | 80 | 12 | 240 |
| 223 | 448 | 37 | 16576 |
| 233 | 52 | 29 | 1508 |
| 425 | 900 | 40 | 1800 |

The table's finite orders are checked by exhaustive return to the identity;
the general formula is proved above. The golden identities and the first
five small-prime mechanisms are inherited from the clock note. At 105,
phi^80=1 modulo 105 but 2^80=46 modulo 105. A register that keeps only the
golden period 80 therefore changes the actual binary source guard. The
joint period is 240. At 19 both orders are 18, but this agreement of periods
does not identify the two operators or supply Collatz legality.

If m=2^K q with q odd, the joint operator is not invertible when K>0.
Its maximal preperiod is K and its exact eventual period is

    lcm(pi_phi(m), ord_q(2)), with ord_1(2)=1.            (7)

The binary 2-power component vanishes after K paddings, whereas the golden
component is always invertible. The state of 1 has exactly this preperiod
and period. Thus mod 18 the pair is (preperiod,period)=(1,24); mod 332 it is
(2,6888), although phi alone has period 168. A periodic golden clock must
not be advertised as reversible storage of the missing source bits.

## 5. Exact polynomial receipts compose through the carrier

Let U(n)=(3n+1)/2^v2(3n+1) on positive odd integers. For a word
w=(a_1,...,a_r) of positive integers, put A=sum a_i and define its ordered
polynomial carry recursively by

    C_empty(z)=0,
    C_(w,a)(z)=3 C_w(z)+z^A.

Then C_w(2) is the familiar affine carry B_w. A receipt consists of positive
odd binary sources n,m, the ordered word w, and a polynomial R satisfying

    3^r P_n(z)+C_w(z)-z^A P_m(z)=(z-2)R(z).             (8)

**PROVED iff:** such an integer R exists exactly when w is the actual
valuation prefix taking n to m. Necessity follows from its scalar endpoint
identity and monic division by z-2. For sufficiency, evaluate at 2.
Reduce this identity modulo 3 to invert the last step. The remaining
odd multiplier and powers of 2 are units modulo 3, forcing the last inverse
numerator divisible by 3; induction does the same at every earlier inverse.
Starting with positive odd m, an integral inverse (2^a m-1)/3 is positive
and odd. Therefore every intermediate endpoint is positive odd and each
valuation is exactly the declared a. This includes the empty word.

Evaluation at phi gives the compatible golden equation

    3^r P_n(phi)+C_w(phi)-phi^A P_m(phi)=(phi-2)R(phi).

But phi-2=-phi^(-2) is a unit in O. Golden equality alone cannot enforce
the binary endpoint identity or its valuations. For example replacing
source 7 by 9 in its one-step receipt to 11 preserves all golden values
in the equation while making the actual arithmetic endpoint false.

For consecutive words u,v with a matching actual intermediate source,

    C_uv=3^len(v) C_u+z^A_u C_v,
    R_uv=3^len(v) R_u+z^A_u R_v.                        (9)

Expanding (8) proves both formulas. Thus polynomial receipts can be stored
and composed without normalizing their golden digits. They preserve order,
source, target, valuation guard, and the binary/golden evaluation pair.
The implementation authenticates the shared binary polynomial, rejecting a
seam that merely has the same golden value. A ROOT receipt has the extra
chronological requirement that no proper prefix reaches 1; this is checked
separately. An empty ROOT receipt is valid only for exact integer source 1.

The state triple alone is not a full receipt: it forgets the ordered word
and the discarded polynomial multiple of (z-2)f. Retain w and R for the
advertised polynomial verification. This is lossless for its specified
values and receipt interfaces, not for arbitrary input history.

## 6. Arithmetic roles, including the distinct source 332

The prior [223/233/322 package](collatz_prime_partition_223_20261007.md)
already proves their common-future families and their scope. The current
332 is 4*83, not the Lucas number 322=L_12=89+233. Exact odd routes give

    332 ->166 ->83 ->125 ->47 ->71 ->107 ->161,

where the first two arrows are ordinary halvings and the remaining arrows
are accelerated odd steps. It therefore joins the already certified
161/322 route. In the Terras convention T(odd n)=(3n+1)/2,

    T^33(332)=425,

whereas the prior clocks are T^7(223)=T^15(233)=T^25(322)=425. The time
sidecar changes when 322 is replaced by 332. The shared destination does
not identify these clocks or turn an arbitrary equal golden value into
an orbit relation.

| Source | Ordinary-binary golden pair (A,B) | D | Initial even halvings |
|---:|---:|---:|---:|
| 105 | (10,15) | 65 | 0 |
| 223 | (18,28) | 149 | 0 |
| 233 | (18,28) | 159 | 0 |
| 322 | (18,30) | 244 | 1 |
| 332 | (20,32) | 248 | 2 |
| 425 | (26,41) | 317 | 0 |

The literal positional expansion
105=1001010101.0101001001_phi and its factor/carry identity are inherited
from the digit note; that word is not the ordinary binary word 1101001.
The inherited roles 233=F_13 and the 223/233 prime subgroup orders remain
separate from the source-value collision (5).

The genuine prior common-future family is

    223+62208t,   233+65536t  ->  425+118098t,  t>=0.

Its authenticated smaller-child proof transports to the larger source.
It is **not** the equal-golden-value family 256t+223,256t+233 from section 3;
their moduli and preserved predicates differ. Both families keep their
exact source labels. No claim of added stopping coverage is made.

| Source -> target | Map / preserved predicate | Lost data / required sidecar |
|---|---|---|
| Integer digit polynomial -> three registers | Both evaluations and all binary residue guards | Original polynomial/history; retain polynomial receipt if required |
| Three registers -> golden pair | Quotient by the binary mode | N and its guards; restore D at the required precision |
| Odd-word receipt -> golden equation | Correct evaluated identity | Binary source/valuation; retain (8), positivity and odd endpoints |
| Verified receipts -> concatenated receipt | Ordered common future | First-hit condition is extra, not a scalar endpoint property |
| Joint modular clock -> padded input | Finite register state after (6) or (7) | Exact height, discarded initial binary bits, and eventual coverage |

## 7. Reproduction and scope

Run the [script](../../04-computation/experiments/golden_guard_arithmetic_20261007.py)
normally and with `python -B -O`; the saved output has the same stem.
The finite universe is every signed coefficient word of lengths 0..5 over
{-1,0,1}; all binary polynomial products of lengths 0..3; sources 0..4095
with guard precisions 1..10; 128 checks of the named infinite collision;
the displayed modular clocks; every prefix of lengths 1..6 at all 512 odd
sources below 1024, with all tested midpoint compositions; the six named
first-hit ROOT routes; and ten typed or false-receipt hostiles.
All 52,488 checks are active in normal and optimized Python.

The positive output is an exact composable guard representation and
polynomial verifier. Neither golden periodicity, the number 9/4, an
algebraic repeller, a Diophantine exceptional list, nor saturated-class
equidistribution by itself supplies a missing source receipt. No external
paper theorem is re-audited or used as a new convergence premise here.
