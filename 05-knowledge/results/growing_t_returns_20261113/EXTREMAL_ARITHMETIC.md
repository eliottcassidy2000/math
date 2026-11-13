# A finite near-square intersection forced by the optimal graph geometry

**PROVED, scoped:** factorization and residue restrictions for the exact
strong unit-ready minimum mu(M), and its complete intersection with the
square-plus-one primes. These are consequences of the independently
reviewed [extremal theorem](inverse/UNIT_READY_MINIMUM.md). They do not
restrict the primes represented by other tournament families.

The inherited [aliquot extremal interface](../uniform_return_metagraph_20261111/arithmetic/PROOF.md)
connects primality with A135244. The actual new source here is the optimal
binary-tail arrangement, whose count has an exact two-parity formula.
The relevant cutoff mechanism is eventual residue stabilization, not an
assumption that every finite collection must share a universal cause.

## 1. Two parity families, with explicit factors

Put r=M-4. For odd r=2q+1, q>=2, set x=2^q. Then

    mu(M)=6x^2+5x-1=(6x-1)(x+1).                    (1)

Thus every odd-order minimum with M>=9 is composite. Both factors exceed1.
The underlying optimal tournament is still strong: factorization of its
Hamilton count is not a decomposition of the graph into ordinal factors.

For even r=2q, q>=2, set x=2^(q-1). Then

    mu(M)=12x^2+7x-1.                               (2)

If q is odd, this is0 modulo3 and exceeds3. Hence every minimum at
M=2 modulo4 with M>=10 is composite. Consequently a prime minimum at
M>=8 is possible only when M is divisible by4. This condition is necessary,
not sufficient; no infinitude of such primes is asserted.

## 2. The square-plus-one cutoff is complete

For M>=11, the powers x in(1),(2) are divisible by8, except for the
already separated M10 case. Both formulas then give

    mu(M)=7 mod8.

A square plus1 is1,2, or5 modulo8. Therefore no mu(M), M>=11, is even
of square-plus-one form. The remaining finite head is

    M:     4  5  6  7  8   9   10
    mu(M): 5  9 17 33 61 115  219.

Only5=2^2+1 and17=4^2+1 are square-plus-one values in this head, and
both are prime. Thus the complete set-theoretic intersection is

    {mu(M): M>=4} intersect A002496 = {5,17}.         (3)

[A002496](https://oeis.org/A002496) concerns *all* primes n^2+1; its
infinitude remains unproved. Equation(3) classifies only this optimal
unit-ready tournament slice. It gives a concrete explanation for its
cutoff: powers of two eventually eliminate every permitted square residue.
No prime search or statistical extrapolation is used in the proof.

## 3. What transfers to enumeration

The sharp order minimum and induced forbidden cores reduce searches for
unit-ready representatives of a prescribed Hamilton count. They do not
remove those graphs from [A000568](https://oeis.org/A000568), which counts
all unlabelled tournaments. No new term of A000568 is claimed.

For A135244, use the inherited extremal equivalence

    A135244(p+1)=p^2 iff p is prime.

It holds for integers p>=2, with A135244 defined as the largest source of
the specified aliquot sum, or0 when there is no source. Applied to p=mu(M),
(1) and the divisible-by3 part of(2) rule out equality at all odd M>=9 and
all M=2 mod4, M>=10. This is an exact rejection rule for that derived
subsequence, not a new general aliquot inversion algorithm. The two
near-square intersections in(3) give the familiar equalities
A135244(6)=25 and A135244(18)=289.

The first75 minimum counts and checks through order1000 are recorded in
`extremal_arithmetic.json`. The proof, rather than those finite checks,
establishes the infinite statements. A targeted OEIS web search for the
displayed minimum prefix found no match in the returned results; this is
not a priority claim or a claim that the sequence is absent from OEIS.
