# Independent audit of the complete square-class lift

**Verdict: the proposed labelled-family classification passes.** This adds a
strictly broader hypothesis than fixed-character or separate-stratum-character
covariance: the symbol need only be invariant under unit squares. It neither
classifies arbitrary tournaments nor gives an isomorphism census at order231.
No producer implementation is imported.

## 1. Existence and exact count

Let N>1 be odd, N=product p^a_p. Consider full signs on nonzero residues satisfying
sigma(-d)=-sigma(d) and sigma(u^2 d)=sigma(d) for every unit u.
On the stratum of additive order m, reduction from U(N) onto U(m) is surjective,
and consequently its square subgroup maps onto U(m)^2. The square-action orbits
are precisely the classes of

    U(m)/U(m)^2 = (F_2)^omega(m).

If some prime p dividing N is1 mod4, negation is a square on the prime stratum
S_p. Square invariance and oddness then contradict each other. This proves the
necessity that every prime factor be3 mod4, without assuming anything about the
other strata.

Conversely assume every p dividing N is3 mod4. On a stratum with r distinct prime
factors, negation acts on its square-class cube by adding the all-ones vector.
That action has no fixed point and pairs its2^r classes into2^(r-1) pairs.
A full odd sign chooses one sign independently on each pair. Thus this stratum
has2^(2^(r-1)) possibilities. Different strata are disjoint and independent.
The total count is2^C, where

    C = sum_{m|N,m>1} 2^(omega(m)-1)
      = [product_{p|N}(1+2a_p)-1]/2.

The product identity follows by expanding over the chosen prime-power exponents:
an omitted prime contributes1 and any positive exponent contributes a factor2.
The m=1 term is removed before division by two. This is a labelled-symbol count,
not the number of tournament isomorphism classes.

The separate-stratum-character family has exponent

    L = sum_{m|N,m>1} omega(m)
      = sum_{p|N} a_p product_{q!=p}(a_q+1).

The inequalities2^(r-1)>=r are equalities exactly at r=1,2 and strict at every
r>=3. Therefore C=L exactly when N has at most two distinct prime factors.
If N has three or more, a divisor containing three primes already yields a
strict contribution. The smallest admissible three-prime modulus is
3*7*11=231, with C=13 and L=12:8192 square-invariant symbols versus4096 locally
character-covariant symbols.

The unit stabilizer remains exactly U(N)^2 for any of these full symbols.
On every prime stratum its odd sign must distinguish the two square classes,
forcing a preserving unit to be a square at each prime. Those conditions suffice
by the assumed square invariance. The all-local-nonsquare units reverse every
stratum, since they add the all-ones vector on every square-class cube. This is
only a statement about unit multipliers, not the full automorphism group.

## 2. Complete three-bit antipodal classification

The relevant involution on F_2^3 is x->x+(1,1,1), explicitly **not** x->-x,
which would be the identity. Write the sign as (-1)^f(x), where

    f(x+111)=f(x)+1.

There are sixteen such truth tables, from four independent antipodal pairs.
An independent algebraic-normal-form calculation classifies all of them.
The derivative of the cubic monomial has nonzero quadratic part, so the cubic
coefficient must vanish. The three quadratic coefficients must coincide, say b,
because the linear coefficients of the derivative vanish. Its constant term gives

    coefficient(x1)+coefficient(x2)+coefficient(x3)+b = 1.

Thus there are exactly two cases:

- b=0: any odd-weight linear coefficient vector, plus an arbitrary constant.
  These are the eight signed odd characters (affine Boolean functions).
- b=1: the quadratic part x1x2+x1x3+x2x3, plus any even-weight linear part and
  an arbitrary constant. These are the other eight functions.

The quadratic part is the Boolean majority function. If s_i=(-1)^x_i, its sign is

    majority(s1,s2,s3) = (s1+s2+s3-s1s2s3)/2.

Multiplying this sign by an even-degree character and a global sign gives the
eight quadratic cases. All these Boolean tables have weight four, consistent
with the separate odd-weight firewall; they are not examples of odd truth-table
weight. The literal checker enumerates all256 Boolean tables and recovers exactly
the stated8+8 split.

## 3. Literal non-character tournament at N=231

On the top stratum U(231), define

    sigma(d)=majority((d/3),(d/7),(d/11)).

On each proper additive-order stratum S_m, choose the Legendre character at the
smallest prime dividing m, evaluated on its normalized unit coordinate.
This defines a sign at every nonzero residue. Each local Legendre value flips
under negation, so the majority sign and the proper-stratum signs are odd.
All unit squares preserve all local Legendre values. Hence the construction is
a full square-invariant circulant tournament symbol.

It is not a character even on U(231). Choose unit square classes(-,+,+) and
(+,-,+). Each has majority value+1, but their product has class(-,-,+) and value-1.
Since sigma(1)=1, no global sign adjustment repairs multiplicativity. The checker
records literal residues realizing these classes, not merely formal class names.

The checker tests all unit-square covariance equations and all unordered vertex
pairs; the exact finite counts appear in the JSON receipt. No expensive isomorphism
census is performed, and no claim is made that this example improves a Ramsey
bound or any other major-conjecture frontier.

## 4. Reproduction

Run `python3 check_square_lift.py` and `python3 -O check_square_lift.py`.
`square_normal.json` and `square_optimized.json` agree byte for byte. The checker
uses actual square orbits on every odd modulus through301, verifies the existence
criterion and counting exponents, exhausts the three-bit antipodal tables, and
constructs the literal231-vertex symbol. Uniform statements follow from the
proofs above rather than extrapolation of these finite controls.
