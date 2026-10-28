# Divisor boxes and the exact overlap at twenty

All results below are **PROVED** by the displayed elementary arguments,
except the explicitly **CITED** Ljunggren input in section4. The finite
checks support, rather than replace, the all-order arguments. The
near-perfect family is classical; the contribution here is its precise
application to the user's proposed overlap and its merger with current
repository work.

## 1 Why cubes and fifth powers have composite counterparts

For `n=product_i p_i^a_i`, every divisor has a unique exponent vector
`(e_i)` with `0<=e_i<=a_i`. Thus

```text
tau(n)=product_i(a_i+1),
P_a(X)=product_i(1+X_i+...+X_i^a_i)-product_i X_i^a_i,
s(n)=P_a(p_1,...,p_r).
```

P has exactly tau(n)-1 monomials, each with coefficient 1. It records
actual proper divisors, before evaluation forgets their identities.
Factoring the small target values of tau gives the complete shape lists:

- tau=4: exponent list(3) or(1,1), hence p³ or pq.
- tau=5: list(4), hence p⁴.
- tau=6: list(5) or(2,1), hence p⁵ or p²q.
- tau=8: list(7),(3,1),(1,1,1).
- tau=9: list(8) or(2,2).

Distinct variables correspond to distinct primes. Repeated substitution
of one prime into two distinct variables would double-count divisors.
This is why scalar regrouping cannot silently replace the divisor box.
In particular20 and 12 share the polynomial
`1+X+X²+Y+XY`, at `(X,Y)=(2,5)` and `(2,3)` respectively.

Prime powers are always deficient: for p>=2 and a>=1,
`s(p^a)=(p^a-1)/(p-1)<p^a`. If p<q are primes, then
`Delta(pq)=1+p+q-pq=2-(p-1)(q-1)<=0`, with equality only at(2,3).

## 2 All nondeficient numbers with five proper divisors

Exactly five proper divisors means tau(n)=6. The prime-power shape p⁵
is deficient by section1. In the other shape p²q, p and q are distinct.
If both are odd, then

```text
sigma(n)/n < [p/(p-1)] [q/(q-1)] <= (3/2)(5/4)=15/8<2,
```

so n is deficient. The two even cases are:

1. `n=4q`, q odd prime. Here `s(n)=7+3q`, so
   `Delta(n)=7-q`. Nondeficiency requires q in{3,5,7}, giving12,20,28.
2. `n=2p²`, p odd prime. Here `s(n)=3+3p+p²`, so
   `Delta(n)=3+3p-p²`. It is 3 at p=3 and negative for p>=5,
   giving18 as the only nondeficient case.

These cases are exhaustive. The corresponding excesses are 4,2,0,3.
In particular20 is the only integer with tau(n)=6 and Delta(n)=2.
Equivalently, it is the only integer with five proper divisors whose
nontrivial proper divisors sum to n+1. QED.

## 3 Every two-prime source with excess two

**Theorem.** Let n have exactly two distinct prime factors. Then
`Delta(n)=2` iff

```text
n=2^a q,   a>=2,   q=2^(a+1)-3 prime.
```

**Proof.** If both primes are odd, the same 15/8 upper bound on the
abundancy excludes Delta(n)>0, regardless of their positive exponents.
Thus `n=2^a q^b` with a,b>=1 and q odd prime. Put `X=2^(a+1)`.
The equation sigma(n)=2n+2 is

```text
(X-1)(q^(b+1)-1)/(q-1)=X q^b+2.
```

Multiplication by q-1 and rearrangement give

```text
(X-q)q^b=X+2q-3,
X=q+3/(1+q+...+q^(b-1)).
```

If b>=2 the denominator is at least q+1>=4. The right side then lies
strictly between q and q+1, contradicting integrality of X. Hence b=1
and X=q+3. The case a=1 gives q=1, so a>=2. Conversely the prescribed
prime q and the divisor formula give Delta(n)=2 exactly. QED.

This is the classical classification in Proposition 2.3 of
[Li–Ma–Zhang (2013)](https://publi.math.unideb.hu/paper/1759/download/).
The incoming [Session83 analysis](../session83_aliquot_divisor_grammar_20261028/FIXED_EXCESS_AND_FIBERS.md)
also gives it using the factorized identity
`Delta(2^a*q^b)+1=(2^(a+1)-q)(1+q+...+q^(b-1))`.
These even sources have redundant divisor 2 in Pollack–Shevelev's
terminology. The proof here does not import the stronger
classification of every possible redundant divisor.
The least example is 20; an integer with three distinct primes is at
least 30 and prime powers are deficient, so20 is also the least positive
integer of any factorization type with excess 2. The statement that all
such integers have two prime factors is false:650 has excess 2 and
factorization2*5²*13.

For every even n with excess 2, `sigma(n)=2(n+1)` has 2-adic valuation one.
It cannot be a square. This includes the complete two-prime family,
so none of its own aliquot edges is square-sum. QED.

**Corollary connecting the aliquot and fourth-power lanes.** In the
two-prime family, the odd prime factor q is 5 modulo 8, since a>=2.
By the exact Eisenstein residue calculation below it needs at least five
fourth-power summands in that ring. This is a necessary representation
obstruction on a specified prime factor, not on the whole source n.
The converse fails at q=37: it is 5 modulo 8 but q+3 is not a power of two.

## 4 The inherited one thousand and one characterization

The square-repunit theorem of Ljunggren, in the formulation of
[Bugeaud–Mignotte–Roy](https://msp.org/pjm/2000/193-2/pjm-v193-n2-p.pdf),
implies `sigma(ell³)` is square for prime ell iff ell=7.
This is **CITED**, not derived from a finite prime search.

For distinct primes p<q<r, the equation s(pqr)=343 is

```text
1+p+q+r+pq+pr+qr=343,
(p+q+1)(p+r+1)=343+p(p+1).
```

The left sum is at least 1+3p+3p², so p<11. The case p=2 gives even s,
whereas343 is odd. For p=3,5,7 the factor-pair targets are 355,373,399.
The smaller factor must be at least p+q+1 for a prime q>p. For355 its
possible smaller factors1,5 are too small; for 373 only 1 is possible.
For399 the smaller factors are 1,3,7,19, and only 19 is large enough.
Its partner21 yields q=11,r=13. Thus the unique source is 1001.

This reproduces [Session82's proof](../session82_decimal_aliquot_square_sum_20261028/ANCHOR_CHARACTERIZATION.md)
and preserves its restricted scope. Complete inverse checking uses a
general bound: if s(n)=t>=2, n is composite; for its least prime divisor p,
the distinct proper divisors 1 and n/p give
`t>=1+n/p>=1+sqrt(n)`. Hence `n<=(t-1)²`.
At t=343 this is 116964. The complete fibre has 19 sources, only one of
which is squarefree with three distinct primes. QED.

## 5 The ring calculation behind the residue budget

With `omega²=-1-omega`, exact expansion gives

```text
(a+b*omega)^4 = (a^4-6a²b²+4ab³)
              +(4a³b-6a²b²+b^4)*omega.
```

Modulo8 this equals0,1,omega,omega² according as `(a,b) mod2` is
`(0,0),(1,0),(0,1),(1,1)`. This follows from an odd integer's square
being1 modulo 8, its fourth power being1 modulo 16, and the displayed
even coefficients. It proves the residue palette for every ring element,
not only for a bounded search of representatives.

For k<=5 summands, let r,s,t count the palette values 1,omega,omega².
Being a rational integer modulo 8 requires `s-t=0 mod8`; as
`|s-t|<=k<8`, this is exactly s=t. The resulting integer coefficient
is r-t, with `r+2t<=k`. The attainable integer coefficients form the
whole interval `[-floor(k/2),k]`: nonnegative values use t=0 and
negative value-j uses r=0,t=j. Modulo8, k=4 misses precisely5, whereas
k=5 reaches every residue. Taking the first k that reaches each residue
gives the table in README. This is a complete modular calculation;
its sufficiency does not assert a global representation of each integer.

## 6 Connection contracts and stopping boundaries

| Source and map | Preserved information | Lost information and required coordinate | Decisive control |
|---|---|---|---|
| Eisenstein coefficients `(a,b)` to parity coordinates `(2a-b,b)` | Lattice addition and exact norm after its metric is transported | Square-grid metric is not the transported norm; retain the factor3 | Six units on the norm-one circle |
| Exponent box to integer weights, then sum | The divisor sum | Individual divisors and exponent shape; retain prime labels and P_a | 12 and 20 share shape, but have excess 4 and 2 |
| Proper-divisor sum to excess `s(n)-n` | The chosen near-perfect equality | Prime support and term count; retain n or factorization | 20 and 650 share excess 2 |
| Two-prime excess-two source to its odd prime q | q=5 modulo 8, hence exclusion from four Eisenstein fourth powers | The congruence alone loses the power-of-two guard q+3 | 13 is a parameter,37 is not |
| Four-term count vectors to integer or Eisenstein residues | The specified representation or congruence predicate | Budget or integral lifting; retain counts and the permitted ring | 37 has a residue obstruction;50 has an integer budget obstruction repaired in L |
| Source n to its inverse-fibre label s(n) | Its aliquot target | Source compatibility and dilation carry | 27 and 35 have target 13 but different doubled outputs |
| Aliquot edge to endpoint sum sigma(n) | Square-sum predicate when n is retained | Arbitrary square adjacency need not be an aliquot edge | 343–57 versus 1001–1001000 |

The concept board changed in three concrete ways. The clarification about20
replaced a proposed partition by simultaneous shape and excess coordinates.
The full two-prime calculation produced a complete family and the 650
boundary. The unit-residue calculation connected the fourth-power budget
to the Eisenstein relation `1+omega+omega²=0`, while the explicit50 witness
located its repair in the first grid. These changes supply actual maps;
they do not infer a universal identification from repeated13/37 counts.

The research cards used are “keep observer, recurrence class, and finite
head” and the source/target type ledger in META-PATTERNS. No new universal
meta-pattern or general classification of all near-perfect numbers is
promoted from these examples.
