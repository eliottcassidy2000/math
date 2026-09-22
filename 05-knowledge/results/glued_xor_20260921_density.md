# Squarefree density, an exact lattice coupling, and two signed triangles

**2026-09-21. Status: PROVED elementary maps, counts and obstructions;
FINITE-EXACT controls; inherited squarefree and divisor-sign densities.**
No novelty claim. No Collatz convergence, cycle classification, or orbit
equidistribution follows.

There is a precise model of two three-cycles exchanged by sign: rotate
primitive integer lattice points through a hexagon. It also carries the
density `6/pi^2`. However, the sixfold multiplicity cancels from the
normalized count, and a fourfold model gives the identical finite density.
The actual common mechanism with squarefree integers is one excluded
residue out of `p^2` at every prime, connected by an explicit CRT bijection.

## 1. Inheritance, scope and the live board

The inherited inputs are the
[proper-divisor means and sign law](catalan_elliptic_20260921_statistics.md),
the [signed Collatz sheet and carry audit](catalan_elliptic_20260921_catalan.md),
and [squarefree residue law SF1](arithmetic_braids2_20260917_squarefree_symmetry.md).
The trace-one matrix behind the hexagonal action also appears in
[the rational-cycle lift](arithmetic_seams_20260921_dynamics.md); its use as
a linear rotation is inherited, while the explicit sieve coupling and
matched square/hexagon counting are the present comparison.

Closest mechanism: a finite CRT product preserves labelled local
divisibility. Canonical hostile: sign-forgetting at Collatz magnitude3
does not preserve the next magnitude. Corrected near miss: doubling a
sample space doubles its numerator and denominator, not its probability.
Least-used sidecar: the relationship between an ordinary integer height
box and a compatible family of residue coordinates.

| Lane / object | Exact retained predicate | Cheap hostile or lost coordinate |
|---|---|---|
| Anchor: squarefree sieve | `p^2` does not divide n | does not refer to a Collatz parameter |
| Niche: CRT digit pairing | p divides both coordinates iff p^2 divides n | addition, multiplication and ordinary height |
| Wildcard: hexagonal rotation | gcd and sixfold orbit | orbit indicator stays constant |
| Square comparison | the same normalized primitive counts | fourfold instead of sixfold symmetry |
| Signed sampling | absolute-value predicate | normalization cancels the extra sheet |
| Collatz dynamics | guarded forward arrows | a distinct map and an unproved orbit law |

The classical visible-point interpretation is also stated in
[Pomerance's author-hosted Euler-phi lecture, slides8--9](https://math.dartmouth.edu/~carlp/PDF/phitalk2.pdf).
All maps and count formulas used below are proved explicitly; no theorem
about the dynamics of visible-point sets is imported from that source.

## 2. An exact finite bijection behind the common Euler product

Let P be any finite nonempty set of primes and `M=prod_(p in P)p`.
For each residue `n mod M^2`, write its p-component uniquely as

```text
n mod p^2 = a_p+p b_p,       0<=a_p,b_p<p.
```

Use CRT separately to assemble `a mod M` and `b mod M` from the digits
`a_p` and `b_p`. This defines a bijection

```text
Phi_P: Z/M^2 Z -> (Z/M Z)^2.                             (DX1)
```

Indeed the digits recover `n mod p^2`, and CRT then recovers `n mod M^2`.
More importantly, for every `p in P`,

```text
p^2|n iff a_p=b_p=0 iff p|a and p|b.                    (DX2)
```

Consequently the image of the finite squarefree sieve is precisely the
finite primitive-pair sieve:

```text
#{n mod M^2 : p^2 does not divide n for every p in P}
 =#{(a,b) mod M : gcd(a,b,M)=1}
 =prod_(p in P)(p^2-1).
```

Uniform sampling on either finite space gives exactly
`prod_(p in P)(1-p^(-2))`. This is a predicate-preserving bijection, not
merely an agreement between two constants.

For `P={2,3}`, the local counts are3 out of4 and8 out of9. Their product
is24 out of36, namely2/3. There are **24**, not six, primitive pairs
modulo6. The remaining primes contribute

```text
prod_(p>=5)(1-p^(-2))=9/pi^2,
(1-1/4)(1-1/9)=2/3,
(2/3)*(9/pi^2)=6/pi^2.                                 (DX3)
```

### What the bijection does not preserve

At p=2, `Phi(1)=(1,0)` and `Phi(2)=(0,1)`. Addition in the target is
componentwise modulo2, so `Phi(1)+Phi(1)=(0,0)`, whereas `Phi(2)=(0,1)`.
Also `Phi(2*2 mod4)=(0,0)` but `Phi(2)^2=(0,1)`. Thus this map is neither
an additive nor multiplicative homomorphism.

It preserves only the declared primes. At `M=6`,

```text
Phi_{2,3}(25)=(1,2),    Phi_{2,3}(29)=(5,0).
```

The first input is not squarefree although its output has gcd1. The
second input is squarefree although its output has gcd5. Both obey (DX2):
prime5 was not in the finite sieve. This is why a finite CRT certificate
must not silently become a global integer bijection.

There is a useful precise completion. At depth h, split the base-p digits
of `n mod p^(2h)` into their even and odd positions, producing two numbers
modulo `p^h`. CRT gives compatible bijections

```text
Z/M^(2h) Z <-> (Z/M^h Z)^2.
```

Passing to all depths and primes gives a measure-preserving homeomorphism
of the underlying profinite probability spaces `Zhat <-> Zhat^2`.
This means exactly that compatible digit cylinders are carried bijectively
to cylinders of the same measure; it is not a group or ring isomorphism.
The all-prime squarefree condition corresponds to a pair not jointly
divisible by any prime in this completed space.

The completion does not preserve ordinary integer representatives. Already
`n=2` has second digit-coordinate1 modulo2 and0 modulo every odd prime.
An ordinary integer with all the latter congruences would have to be zero,
contradicting the former. Thus this explicit completed map does not send
every ordinary integer to an ordinary integer pair. It supplies no
archimedean height, ordering, or Collatz-time transport.

## 3. Passing from finite sieves to the two spatial densities

Let mu denote the Mobius function. The squarefree count in `1,...,X` is

```text
Q(X)=sum_(d<=sqrt X)mu(d) floor(X/d^2)
    =X/zeta(2)+O(sqrt X).
```

This is inherited from the statistics note. For ordered positive lattice
pairs in the square `1<=a,b<=N`, write `V(N)=#{gcd(a,b)=1}`. The analogous
finite identity is

```text
V(N)=sum_(d<=N)mu(d) floor(N/d)^2
    =N^2/zeta(2)+O(N log N).                            (DX4)
```

To prove the first equality, sum `sum_(d|gcd(a,b))mu(d)` over pairs.
For the second, replacing each squared floor costs `O(N/d)` and the sum
of those errors is `O(N log N)`; the reciprocal-square tail beyond N
costs `O(N)`. Absolute convergence gives
`sum mu(d)/d^2=prod_p(1-p^(-2))=1/zeta(2)=6/pi^2`.

These estimates justify the limiting densities without assuming that
all growing-prime events are independently sampled on a finite height box.
The finite height populations are not exactly identified by (DX1):
`Q(9)=6`, while `V(3)=7`. Agreement of their limiting normalized counts
does not make their finite boxes equal.

## 4. A genuine two-triangle model, and a fourfold hostile with the same counts

Consider the integer linear map

```text
R(a,b)=(-b,a+b),
R^2(a,b)=(-a-b,a),     R^3(a,b)=(-a,-b),     R^6=I.     (DX5)
```

Its determinant is1. Therefore gcd is preserved, as is
`a^2+ab+b^2`. It also preserves the hexagonal height
`h(a,b)=max(|a|,|b|,|a+b|)`.

Every nonzero integer vector has exact R-period6. A smaller period would
give a fixed vector of R, R^2 or R^3. The matrices `R-I`, `R^2-I` and
`R^3-I` have determinants1,3 and4 respectively, so their rational kernels
are zero. Each six-orbit splits into two exact three-cycles under R^2,
exchanged by the central sign `R^3=-I`. For example,

```text
R-orbit: (1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1).
R^2 triangle: (1,0),(-1,1),(0,-1);
its negative: (-1,0),(1,-1),(0,1).
```

This realizes the suggested two signed copies of a three-cycle on a
specified carrier, with a specified map. It is not yet a Collatz carrier.

Define `H_N={v in Z^2:0<h(v)<=N}` and let `P_H(N)` count its primitive
vectors, meaning gcd1. The shell `h(v)=r` is the disjoint union of the six
R-images of

```text
{(r-t,t):0<=t<r}.
```

This is a half-open side of the hexagon, so each corner is counted once.
Its r elements contain `phi(r)` primitive vectors because
`gcd(r-t,t)=gcd(r,t)`; use `phi(1)=1`. Summing shells proves

```text
|H_N|=6 sum_(r<=N)r=3N(N+1),
P_H(N)=6 sum_(r<=N)phi(r).                              (DX6)
```

The sixfold factor appears in both numerator and denominator. It cancels
from `P_H(N)/|H_N|`. Splitting every orbit into its two triangles likewise
changes neither the primitive fraction nor its limit.

There is a decisive alternative symmetry. Let
`S_N={v in Z^2:0<max(|a|,|b|)<=N}` and rotate by
`J(a,b)=(-b,a)`. This is a free order-four action on nonzero vectors and
also preserves gcd. Its total count is `4N(N+1)`. Primitive vectors with
both coordinates nonzero contribute `4V(N)`; the four unit-axis points
contribute4. Since `V(N)=2 sum_(r<=N)phi(r)-1`,

```text
P_S(N)=8 sum_(r<=N)phi(r),
P_H(N)/|H_N| = P_S(N)/|S_N|
             =2 sum_(r<=N)phi(r)/(N(N+1))
             =(V(N)+1)/(N(N+1)) -> 6/pi^2.              (DX7)
```

Thus fourfold and sixfold models have the **same exact primitive fraction
for every positive integer N**. The numerator6 in the final zeta value
cannot be inferred from sixfold orbit multiplicity alone. The arithmetic
primitive predicate and the normalization determine the count.

Two boundary distinctions are essential. First, freeness in (DX5) is an
integer-lattice statement: modulo2 the primitive residues form one
three-cycle; modulo3 they form one two-cycle and one six-cycle. Modulo6
there are four six-cycles. Reducing a lattice action can collapse periods.
Second, even in the genuine integer six-cycle model,

```text
1_{gcd=1}(R^j v)=1_{gcd=1}(v) for every j.
```

Every orbit's primitive time average is therefore exactly0 or1, while
the spatial limiting proportion is `6/pi^2`. The primitive orbit through
`(1,0)` and the nonprimitive orbit through `(2,0)` witness both values.
An actual two-triangle structure does not imply the required ergodicity.

## 5. Signing, residue conditioning and Collatz parameter variation

For uniform `n in {-X,...,-1,1,...,X}`, absolute squarefreeness has exactly
the same finite probability as on `1,...,X`:

```text
2Q(X)/(2X)=Q(X)/X.
```

Negation is a genuine second sheet, but it cannot double a normalized
density. The inherited conditioned limits are

```text
P(squarefree | odd)=8/pi^2,
P(squarefree | n=1 mod6)=P(squarefree | n=5 mod6)=9/pi^2,
P(squarefree | n=3 mod6)=6/pi^2.                         (DX8)
```

Oddness removes the local square obstruction at2. Rows1 and5 also remove
it at3; in row3 exactly two of the three lifts modulo9 avoid9. Negation
exchanges rows1 and5 and preserves row3. Those facts explain the constants
through explicit local conditioning, without counting dynamical cycles.

The divisor defect `D=F-S-U` from the prior statistics note is negative
exactly on the squarefree set outside a zero-density exceptional set.
Consequently (DX8) also gives its row-conditioned negative-sign densities.
This inherited consequence concerns uniform sampling of integers.

Squarefreeness itself does not depend on the Collatz additive parameter.
For an explicit dynamic contrast, use `U_b(n)=(3n+b)/2^v2(3n+b)` on
nonzero odd integers where the numerator is nonzero. A fixed point obeys

```text
(2^k-3)n=b,            k>=1.
```

For b=1 the complete fixed-point list is `-1,1`; for b=5 it is `-5,5,1`,
with exponents respectively `(1,2)` and `(1,2,3)`. Completeness follows
by listing divisors `2^k-3` of b; since `2^k-3>=-1`, no further signed
divisor is possible in either case. The uniform squarefree density is
the same for both maps. Thus that density alone cannot recover even the
number of fixed points in this family, much less classify its longer
cycles. This does not purport to rule out every more elaborate proposed
connection; it identifies the missing information in a density-only claim.

## 6. Connection contract, exact controls and stopping point

| Source -> target | Map | Preserves | Does not supply |
|---|---|---|---|
| Finite squarefree sieve -> primitive residue pairs | prime-digit splitting and CRT | each declared p-square exclusion | primes outside the sieve; ordinary height; ring operations |
| Completed squarefree sieve -> completed primitive pairs | interleaving all local digits | cylinder measures and all local exclusions | ordinary integer representatives |
| Hexagonal lattice -> two triangles | powers R^2 with sign R^3 | gcd and exact finite periods | a Collatz semiconjugacy or spatial-to-time equality |
| Hexagon -> square comparison | equal totient-count ratios | exact primitive proportion | sixfold symmetry as a unique explanation |
| Positive -> signed sample | adjoining a sign | absolute-value predicate | change of normalized density |

Run from the repository root:

```text
python 04-computation/experiments/glued_xor_20260921_density.py
python -O 04-computation/experiments/glued_xor_20260921_density.py --output C:/tmp/glued_xor_density_optimized.json
```

The [standard-library script](../../04-computation/experiments/glued_xor_20260921_density.py)
and [JSON](../../04-computation/experiments/glued_xor_20260921_density.json)
check93,026 residues across eight complete digit/CRT universes, including
depths1,2,3. They check every vector in the square and hexagon at
`N=1,2,3,5,10,25,50,100`, every indicated orbit and gcd, and every
hexagonal shell at those scales. Local rotation periods are enumerated
modulo2,3,5,6,30. The finite squarefree census uses all integers through
600,000 and its three odd residue rows. Positive and hostile controls
include the nonhomomorphism, omitted-prime, height-box, sign-normalization
and parameter-change examples above. Explicit checks survive Python `-O`.

The genuine survivor is an exact square-exclusion/primitive-pair coupling
and a geometric realization of two signed three-cycles. The obstruction
is equally exact: symmetry multiplicity cancels, the CRT map loses
ordinary height, and individual rotation orbits have time average0 or1.
A Collatz interpretation must provide another explicit map preserving
legal steps and the needed trajectory quantity; those obligations remain
OPEN here.
