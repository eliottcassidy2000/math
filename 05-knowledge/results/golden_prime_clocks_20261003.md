# Golden-ratio digit clocks: prime quotients, 105, and source-guard registers

2026-10-03 (America/Denver).

**Status:** PROVED elementary ring and clock statements below; FINITE-EXACT
for the stated script universes; CITED for the historical first-nonflat
cyclotomic result. These clocks preserve modular registers. They give no
primality criterion for arbitrary digit words and no universal Collatz coverage.

**Inheritance.** The closest mechanism is the Lucas/monodromy quotient in
[the September 30 circulant note, section 2](collatz_circulant_20260930_circulants_lucas_cubic_monotile.md).
The current guard mechanism is the exact phase language in
[compressed route guards, sections 1–2](collatz_fourier_guards_20261003.md),
which inherits [THM-4501, recursive motif families](../../01-canon/theorems/THM-4501-collatz-recursive-motif-families-and-frequency.md).
The hostile example is a nonzero nilpotent modulo 5: nonzero does not permit
division. The corrected near miss is the signed norm of phi^4−1: it is −5,
while the quotient has cardinality 5. The least-used sidecar is the second
coordinate of a+b phi, together with the distinction between the ideals (p)
and (phi^n−1). The companion
[digit and carry reader](golden_digit_carry_20261003.md) owns the exact
positional/Zeckendorf conversion and the literal decimal-105 expansion.

The live concepts are digit carries, ring norms, prime splitting, finite
shift clocks, and exact source guards. The transfer developed here is a
finite register for the same digit reader, rather than an identification
of these objects with a Collatz orbit.

## 1. One exact two-coordinate reader

Let O=Z[phi], phi^2=phi+1, and represent a+b phi by (a,b). Then

    (a,b)(c,d) = (ac+bd, ad+bc+bd),
    phi = (0,1),       phi^(-1) = (-1,1),
    conjugate(a,b) = (a+b,-b),
    N(a,b) = a^2+ab-b^2.

Multiplication by phi on coefficient columns is

    M = [[0,1],[1,1]],       det M = −1.

Consequently M is invertible modulo every positive integer, including 105.
The affine reader for a digit d is (A,B) -> (B+d,A+B). Negative positional
exponents need only powers of the exact unit phi^(-1). In particular,
`11_phi = 100_phi` is the carry relation 1+phi=phi^2; it is not an equality
between the decimal integers 11 and 100.

For any u=(a,b), the multiplication matrix is
[[a,b],[b,a+b]], of determinant N(u). Thus

    u is a unit in O/(m)  iff  gcd(N(u),m)=1.

The sufficiency follows by multiplying the conjugate by the inverse of
N(u) modulo m. Necessity follows by taking determinants of an inverse
multiplication map. For m=1 the zero ring has its single identity; the API
uses the consistent period-one convention.

**Hostile.** In O/(5), eta=phi−3 is nonzero but eta^2=0. In O/(105), the
nonzero element phi−3 has norm 5, and its product with (63,84) is zero.
The shift M is invertible, but arbitrary digit values are not. Any algorithm
that divides by a represented value must retain the norm/unit guard.

## 2. The two identities are trace identities, with different prime clocks

Define Fibonacci numbers F_0=0, F_1=1 and Lucas numbers L_0=2, L_1=1, each
with the same two-term addition recurrence. Direct induction gives

    phi^n = F_(n−1)+F_n phi,
    M^n = [[F_(n−1),F_n],[F_n,F_(n+1)]],
    phi^(2n)+(-1)^n = L_n phi^n,
    N(phi^n−1) = 1+(-1)^n−L_n.

The third identity also follows from conjugation phi'=-phi^(-1). At n=4,5
it is exactly

    phi^8+1 = 7 phi^4,       phi^10 = 1+11 phi^5.

The signed norm at n=4 is −5, not +5; the finite quotient by the nonzero
principal ideal (phi^4−1) has cardinality |−5|=5. Cardinality is the absolute
determinant of its multiplication matrix.

For an odd prime p≠5 the discriminant of X^2−X−1 is 5.

* If 5 is a square modulo p, O/(p) is F_p × F_p and phi^(p−1)=1.
  Its order divides p−1; equality is not asserted.
* If 5 is a nonsquare modulo p, O/(p) is F_(p^2). Frobenius exchanges the
  roots, so phi^p=1−phi=−phi^(-1). Therefore phi^(p+1)=−1: its order divides
  2(p+1) and does not divide p+1.
* Modulo 2 the irreducible quadratic gives F_4 and phi has order 3.
* Modulo 5, write phi=3+eta=3(1+2eta), eta^2=0. The first factor has order
  4 and the second order 5, so phi has order 20.

The exact small clocks are as follows. The residue-ring type is part of
the statement, not an optional decoration.

| Modulus | Ring | Order of phi | Decisive identity |
| --- | --- | ---: | --- |
| 3 | F_9 | 8 | phi^4=−1 |
| 5 | F_5[eta]/(eta^2) | 20 | phi=3(1+2eta) |
| 7 | F_49 | 16 | phi^8=−1 |
| 11 | F_11 × F_11 | 10 | phi maps to (4,8) |
| 105 | product of the 3,5,7 rows | 80 | lcm(8,20,16) |

For 3 and 7, the displayed half-turn forces the stated power-of-two order.
For 11, phi^10=1 while phi^2 and phi^5 are not 1, excluding the proper
divisors needed to establish order 10. The Chinese remainder theorem
applied to both coordinates proves the 105 row.

There are 8·20·48=7680 units modulo 105. The 80 powers of phi visit only
1/96 of them. A cyclic shift clock is not an enumeration of all units.

**Two quotient clocks at the same prime.** The inherited monodromy quotient
O/(phi^5−1) is F_11 with phi mapping to 4, of order 5. Indeed phi^5−1=2+5phi
forces phi=4 modulo 11, and the quotient has size 11. The integer-modulus
quotient O/(11) retains both roots (4,8) and has shift order 10. Forgetting
which ideal defined the quotient loses precisely this branch information.

## 3. Exact schedules for binary/ternary guard registers

Let pi(m) denote the order of M modulo m. For a,p>=1,

    pi(2^a) = 3·2^(a−1),       pi(3^p) = 8·3^(p−1).

Here is a direct lifting proof, without an assumption about other primes.
We have M^3−I=2M and M^6−I=4M^3. The latter has entrywise 2-adic valuation
exactly 2 and, after division by 4, is I modulo 2. If
M^(3·2^r)=I+2^(r+1)U for r>=1 and U is I modulo 2, squaring gives

    M^(3·2^(r+1))−I = 2^(r+2)(U+2^r U^2).

The bracket is still I modulo 2, so the exact valuation rises by one.
Together with the order 3 modulo 2, this proves the first formula.

Similarly,

    M^8−I = 3 [[4,7],[7,11]].

Cubing I+3^s B, s>=1, gives I+3^(s+1)B' with B'=B modulo 3. The initial
B is nonzero modulo 3. Its exact entrywise 3-adic valuation therefore rises
by one at each cubing. The order 8 modulo 3 proves the second formula.

CRT now gives, allowing a,p>=0 and using pi(1)=1,

    T(a,p) = lcm( a==0 ? 1 : 3·2^(a−1),
                  p==0 ? 1 : 8·3^(p−1) ).

When both exponents are positive this simplifies to

    T(a,p) = 2^max(a−1,3) · 3^max(1,p−1).

For example, modulo 2^8·3^3=6912 the exact clock is 1152. The script exports
`guard_clock(a,p)` and signed `phi_power(k,modulus)`, so a Laurent reader
can correct its recorded fractional length without a modular division risk.
No blanket prime-power lifting law is claimed for primes other than 2 and 3.

**Typed compiler connection.** Read a digit word to its pair modulo
2^a3^p. Appending T(a,p) zero digits multiplies the state by M^T=I. This
preserves both modular coordinates and hence every predicate on that pair.
It does not preserve exact size, integrality as an ordinary rational integer,
or a prescribed first-hit route. If the represented object is an integer n,
its state is (n,0); after a clock shift the residue remains (n,0), but the
exact value will generally no longer be a rational integer.

| Source | Target and map | Preserved predicate | Loss and required sidecar |
| --- | --- | --- | --- |
| Exact Laurent digit word | Pair in O/(2^a3^p), by Horner and phi^(-R) | Every source congruence in the retained pair | Exact value, magnitude and rational integrality; keep the exact pair and R |
| Shift of digit positions | Clock phase modulo T(a,p) | Modular reader state after a full clock | Actual position; keep signed position or length |
| Certified route source | Its binary/ternary guard registers | The stated finite congruence guard | A root exception and the certified tail are not periodic state alone |

The exact valuation bit matters. A route with division budget A commonly
requires a source guard modulo 2^(A+1), not only 2^A; use a=A+1 when that is
the guard's actual modulus. For example n=1 and n=5 are identical modulo 4,
but v2(3n+1) is respectively 2 and 4. Similarly the isolated (m,t)=(1,0)
first-hit exception in the inherited phase compiler cannot be removed by
deleting an entire periodic class. A finite register organizes a proved
guard; it does not prove that every input eventually reaches that guard.

## 4. A separate meaning of 105: the first nonflat cyclotomic polynomial

This section treats 105 as an index, separately from the modulus in section 2.
Paul Garrett's primary
[note on the 105th cyclotomic polynomial](https://www-users.cse.umn.edu/~garrett/m/algebra/notes_2023-24/105th_cyclotomic_poly.pdf)
identifies 105=3·5·7 as the first index with a coefficient outside {−1,0,1}.
Its degree is 48 and its coefficients at degrees 7 and 41 are −2. The
script independently computes every Phi_n for 1<=n<=105 by exact monic
division and verifies these claims, including the entire lower range.

A second exact check uses

    Phi_105(X) = ((X^105−1)(X^3−1)(X^5−1)(X^7−1))
                 /((X^35−1)(X^21−1)(X^15−1)(X−1)).

Evaluating in the same two-coordinate reader gives

    Phi_105(phi) = 5224431949 + 8453308464 phi,
    N(Phi_105(phi)) = 16271615641 = 21211·767131.

This is an exact arithmetic intersection with the base-phi reader, not a
claim that the coefficient word is binary. In fact its two −2 coefficients
are the precise obstruction to that description. The value's norm is
composite; the two listed factors are prime, checked by exact trial division.

For an odd prime index q, Phi_q(X)=1+X+...+X^(q−1) is a binary repunit.
Since phi−1=phi^(-1),

    Phi_q(phi) = phi(phi^q−1),       N(Phi_q(phi)) = L_q.

Thus its principal ideal is the inherited ideal (phi^q−1): the binary
repunit and monodromy quotient are exactly connected by a unit. Prime
index does not make this an algebraic prime. The concrete hostile is

    Phi_23(phi) = phi^12 (11+2phi)(21+phi),
    N(11+2phi)=139,      N(21+phi)=461,
    L_23=64079=139·461.

Both factors are nonunits in O, so this is a proper factorization. Conversely
composite norm alone is not a universal disproof of algebraic primality:
Phi_3(phi)=2phi^2 generates (2), whose quotient is the field F_4, even though
its norm is 4. Norm a rational prime is sufficient for a principal element
to be algebraically prime; fieldness of the quotient is the exact criterion.

These distinctions separate three questions: a prime integer index, a binary
positional digit word, and a prime element of O. The 23 example breaks their
putative equivalence while preserving the exact repunit-to-quotient map.

## 5. Reproduction and boundaries

Run from the repository root:

    python 04-computation/experiments/golden_prime_clocks_20261003.py

Saved output: [golden_prime_clocks_20261003.out](golden_prime_clocks_20261003.out).
The script uses integer arithmetic only and explicit failure checks, including
under `python -O`. Its independent paths are direct Fibonacci-register
iteration versus pair/matrix powers; brute inverse existence versus the norm
criterion; and recursive cyclotomic division versus the squarefree product.

The finite universe is 95 primes <=499, all moduli 1..128, all 89440 pairs
modulo 1..64 (650 with an independent brute inverse test), all 117 mixed
moduli with 0<=a<=12 and 0<=p<=8, 20 direct prime-power cycles, and all
cyclotomic indices 1..105. Positive controls include the literal Lucas
identities and exact 105 clock. Hostile controls include nilpotents, the
signed norm, a lost valuation bit, exact magnitude loss, prime-index
factorization, and the norm-4 prime ideal.

The useful next operation is to attach these small reversible modular
registers to the exact digit/carry reader and evaluate an already-proved
source guard there. The unresolved part remains the production or coverage
of certified routes, not the modular clock arithmetic.
