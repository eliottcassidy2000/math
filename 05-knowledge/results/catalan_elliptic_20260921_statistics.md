# Proper-divisor statistics: the squarefree constant controls a mean and a sign law

**2026-09-21. Status: PROVED exact identities and density consequences below;
CITED classical prime asymptotics and normal-order theorem; FINITE-EXACT
controls.** No novelty claim. The probability space is explicitly uniform
integers in `1,...,X`, with integer `X>=1`; no Collatz orbit measure is assumed.

The useful new consequence of the inherited divisor classification is

```text
density{F-S-U<0}=6/pi^2,
density{F-S-U=0}=0,
density{F-S-U>0}=1-6/pi^2.                                (ST1)
```

At the same time `E_X(F-S-U)~(1-6/pi^2)log X>0`. Thus the limiting
majority sign is negative while the mean is eventually positive. Both
statements have exact proofs; neither can replace a trajectory argument.

## Inheritance and live board

The closest proved mechanism is the exponent box versus Boolean subcube
in [the divisor-balance classification, DB1--DB3](arithmetic_braids_20260917_divisors.md),
ultimately attached to
[THM-2422, operation fibres and twin-center ancestry](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md).
That classification and the
[squarefree residue densities, SF1](arithmetic_braids2_20260917_squarefree_symmetry.md)
are inherited, not rediscovered claims. The canonical hostile is `p^2`:
removing a repeated-factor diagonal destroys a genuine divisor witness.
The corrected near miss here is treating `omega`, `Omega`, and the user's
proper-prime count as the same function. The least-used sidecar is the
sampling measure and its endpoint corrections.

| Lane / concept | Preserved quantity | Hostile or missing coordinate |
|---|---|---|
| Anchor: divisor box | exact F,S,U and their defect | n=1, a prime, p^2, p^3 |
| Niche: incidence summation | exact uniform mean | one fixed integer is not random |
| Wildcard: sign of the defect | squarefreeness outside a sparse exception set | majority sign need not equal mean sign |
| Euler product | local exclusion of p^2 | fixed finite CRT independence is not uniform growing-sieve independence |
| Collatz residue rows | conditioned spatial density | no invariant orbit measure or descent inequality supplied |

The source-to-target maps are divisor incidence followed by summation,
and prime exponents followed by the sign of the defect. The first forgets
the order in which integers are sampled; the second forgets magnitudes and
most exponent data. An orbit conclusion needs those coordinates restored.

## 1. Definitions, with both endpoint corrections retained

For every positive integer n, use the original proper, nontrivial sets:

```text
F(n)=#{d|n:1<d<n},
S(n)=#{d|n:1<d<n, d squarefree},
U(n)=#{p|n:p prime, p<n}.
```

Let `tau(n)` count all positive divisors, `omega(n)` distinct prime
divisors, `Omega(n)` prime divisors with multiplicity, and `mu` be the
Mobius function. Then, including n=1,

```text
F(n)=tau(n)-2+[n=1],
S(n)=2^omega(n)-1-mu(n)^2+[n=1],
U(n)=omega(n)-[n prime].                                 (ST2)
```

Thus `F(1)=S(1)=U(1)=0`, and all three also vanish at a prime.
The phrase "number of prime factors" is otherwise ambiguous: for example
at `n=8`, `(U,omega,Omega)=(1,1,3)`. At a prime these are `(0,1,1)`.
The paper of Hardy--Ramanujan uses its own letters `f,F` for the two
prime-factor counts; its `F` must not be substituted for the user's `F`.

Write `E_X f=X^(-1)sum_(1<=n<=X)f(n)`, `Q(X)=sum_(n<=X)mu(n)^2`, and
`pi(X)=#{p<=X:p prime}`. Double-counting the pairs `(d,n)` with `d|n`
gives exact finite formulas:

```text
X E_X F = sum_(d<=X) floor(X/d)-2X+1,
X E_X S = sum_(d<=X) mu(d)^2 floor(X/d)-X-Q(X)+1,
X E_X U = sum_(p<=X) floor(X/p)-pi(X),
X E_X Omega = sum_(p^a<=X,a>=1) floor(X/p^a).              (ST3)
```

These are identities for each integer X, not asymptotic probabilities.
The `+1` and `-pi(X)` terms matter for literal finite comparisons.

## 2. What the logarithms actually describe

**CITED prime input.** The prime number theorem implies
`P_X(n prime)=pi(X)/X~1/log X`. It is a statement about the uniform prefix,
not a probability attached to a fixed known integer or an arbitrary short
interval. A primary proof is
[Zagier, Newman's short proof of the prime number theorem (1997), pp.705--708](https://people.mpim-bonn.mpg.de/zagier/files/doi/10.2307/2975232/fulltext.pdf).

**CITED prime-harmonic input.** We use
`sum_(p<=X)1/p=log log X+B_1+O(1/log X)` and
`pi(X)=O(X/log X)`. Theorem2.7(d), Corollary2.6 and equation(2.22) of
[Montgomery--Vaughan, author-hosted Chapter2](https://personal.science.psu.edu/rcv4/personal/Publications/MNTI/06.0_pp_35_75_The_elementary_theory_of_arithmetic_functions.pdf)
give these estimates and the corresponding mean of omega. Consequently
the third identity of (ST3) yields

```text
E_X U = log log X+B_1+O(1/log X).                         (ST4)
```

Counting multiplicity changes the constant, not this leading scale:

```text
E_X Omega = log log X+B_1+sum_p 1/(p(p-1))+o(1).
```

Indeed the additional incidence sum is over `p^a`, `a>=2`. Its limiting
mean is the convergent sum `sum_p sum_(a>=2)p^(-a)`; the floor errors are
`O(sqrt(X)log X)/X`, and the omitted reciprocal tail tends to zero.

For the full divisor count a direct hyperbola decomposition gives, with
`r=floor(sqrt X)`,

```text
sum_(n<=X)tau(n)=2 sum_(d<=r)floor(X/d)-r^2
               =X log X+(2 gamma-1)X+O(sqrt X).
```

The second equality follows by replacing the floors and using
`sum_(d<=r)1/d=log r+gamma+O(1/r)`. Hence

```text
E_X F = log X+2 gamma-3+O(X^(-1/2)).                     (ST5)
```

Thus the two proposed scales belong to different counts: the mean of F
is logarithmic, while the mean of U is doubly logarithmic. All logarithms
here are natural.

**CITED normal-order distinction.** Hardy--Ramanujan's
[original 1917 paper, sections I--II](https://ramanujan.sirinudi.org/Volumes/published/ram35.html)
shows that both omega and Omega have normal order `log log n`. Together
with `2^omega(n)<=tau(n)<=2^Omega(n)`, this gives
`tau(n)=(log n)^(log 2+o(1))` outside a set of density zero. The same
logarithmic exponent applies to F and S by (ST2). This is smaller than
their arithmetic-mean order `log X`; replacing a mean by a typical value
would lose the contribution of the larger values.

## 3. Why the squarefree density also multiplies the mean of S

**PROVED here from elementary convolution.** The two local identities

```text
mu(n)^2=sum_(r^2|n)mu(r),
2^omega(n)=sum_(d|n)mu(d)^2                              (ST6)
```

follow respectively by excluding exponents at least2 and by choosing a
subset of the prime support. The first yields

```text
Q(X)=sum_(r<=sqrt X)mu(r)floor(X/r^2)
    =c X+O(sqrt X),
c=sum_(r>=1)mu(r)/r^2=prod_p(1-p^(-2))=1/zeta(2)=6/pi^2.
```

Absolute convergence bounds the tail by `O(X^(-1/2))`; the floor errors
number at most `sqrt X`. The same local square exclusion appears in the
absolutely convergent Dirichlet series, for `Re(s)>1`:

```text
sum mu(n)^2/n^s = zeta(s)/zeta(2s),
sum 2^omega(n)/n^s = zeta(s)^2/zeta(2s).                 (ST7)
```

The second series has one extra divisor-incidence factor `zeta(s)`.
This explains the common constant without relying on an unproved
independence model. A proof of its two-term mean avoids any Tauberian
assumption: writing `T(y)=sum_(n<=y)tau(n)`, (ST6) gives

```text
sum_(n<=X)2^omega(n)=sum_(r<=sqrt X)mu(r)T(X/r^2)
                   =c X log X+C X+O(sqrt X log X),
C=(2 gamma-1)/zeta(2)-2 zeta'(2)/zeta(2)^2.              (ST8)
```

Insert the hyperbola estimate for T. The summed errors are bounded by
`sqrt X sum_(r<=sqrt X)1/r`. Extending the two main absolutely convergent
series to infinity costs `O(sqrt X log X)`. Finally
`sum mu(r)log(r)/r^2=zeta'(2)/zeta(2)^2`, by differentiating the absolutely
convergent series for `1/zeta(s)` at2. Combining (ST2), (ST8) and Q gives

```text
E_X S=c log X+C-1-c+O(log X/sqrt X),
E_X S/E_X F -> c,             E_X U/E_X F -> 0,
E_X(F-S-U)=(1-c)log X-log log X+O(1).                    (ST9)
```

The ratio here is a ratio of means. It is not the mean of S/F, which
also requires a convention at primes and at1 where F=0. The Euler product
relates two specified sums; it does not make divisor events on one
integer independent.

## 4. The sign of the divisor defect detects squarefreeness almost everywhere

Put `D(n)=F(n)-S(n)-U(n)`. The inherited exact classification says:

```text
D=0 iff n=1, p, p^3, or p^2qr;
D<0 iff n is squarefree composite, p^2, or p^2q,
```

where letters in each mixed pattern denote distinct primes. For clarity,
the decisive general inequality is elementary. At a nonsquarefree n with
`r=omega(n)`,

```text
D(n)=tau(n)-2^r-r-1 >= 2^(r-1)-r-1 > 0   if r>=4.      (ST10)
```

At a squarefree composite, D=-r. Thus the symmetric difference of the
sets `{D<0}` and `{n squarefree}` is contained in `{omega(n)<=3}`.

**PROVED sparse-exception lemma.** For every fixed K, the set
`{n:omega(n)<=K}` has natural density zero. Here is a finite-CRT proof,
which does not require the normal-order theorem. For a finite prime set
P, let `h_P(n)=#{p in P:p|n}`. On a full residue system modulo the product
of P, CRT gives exactly

```text
E z^h_P=prod_(p in P)(1-1/p+z/p).
```

For fixed `0<z<1`, the indicator of `omega<=K` is bounded above by
`z^(-K) z^h_P`, since `h_P<=omega`. Periodic averaging therefore gives

```text
upper_density{omega<=K}
 <= z^(-K) prod_(p in P)(1-(1-z)/p).
```

As P exhausts the primes, the product tends to zero because the sum of
prime reciprocals diverges. This last fact already follows from the
prime-harmonic estimate used in (ST4); no quantitative rate is needed.
The order of limits is fixed P, then X to infinity, then enlarge P.

Applying the lemma with K=3 and using Q(X) proves all three statements
in (ST1). Equality in the user's `F=S+U` equation is therefore a
zero-density event, despite its three infinite prime-exponent families.

There is no contradiction between the majority-negative sign law and
the positive mean (ST9): D is unbounded. The positive values have enough
total weight to outweigh the more numerous negative values. For a
positive/hostile pair, `210=2*3*5*7` has D=-4, while
`420=2^2*3*5*7` has D=3. These are fixed examples, not a rate estimate.

## 5. The odd Collatz rows change the constant, not the proof obligation

The inherited squarefree row law has an immediate new consequence for D.
Relative to uniform positive integers `n<=X` in the indicated row,

| Row | Limiting P(D<0) | Limiting P(D=0) | Limiting P(D>0) |
|---|---:|---:|---:|
| 1 modulo6 | 9/pi^2 | 0 | 1-9/pi^2 |
| 3 modulo6 | 6/pi^2 | 0 | 1-6/pi^2 |
| 5 modulo6 | 9/pi^2 | 0 | 1-9/pi^2 |

For rows1 and5, neither2 nor3 divides n, leaving the product
`prod_(p>=5)(1-p^(-2))=9/pi^2`. In row3 exactly two thirds of its
lifts modulo9 avoid9, giving the factor2/3. A fixed cutoff CRT sieve
followed by the square-divisor tail proves these densities; this is SF1
in the inherited note. The zero-density exceptional set in section4
remains zero-density relative to each row because each row has positive
ambient density. Consequently the sign law transports to these rows.

For signed sampling on `{-X,...,-1,1,...,X}`, define each function at n
using `|n|`. The unconditional finite distribution is then exactly the
positive distribution. Negation interchanges rows1 and5 and preserves
row3. This signed reflection does not produce forward-orbit sampling.

The missing map to Collatz is explicit: uniform spatial density does not
identify frequencies along a single deterministic orbit. Even a proposed
orbit frequency law would still need a legal-step relation and a descent
quantity. The existing
[simultaneous-squarefree growth theorem](arithmetic_braids2_20260917_squarefree_symmetry.md)
already supplies arbitrarily long increasing prefixes whose every odd
node is squarefree. Nothing in the present sign or mean law contradicts
those prefixes or upgrades them to infinite trajectories.

## 6. Exact controls and stopping boundary

Run from the repository root:

```text
python 04-computation/experiments/catalan_elliptic_20260921_statistics.py
python -O 04-computation/experiments/catalan_elliptic_20260921_statistics.py --output C:/tmp/catalan_elliptic_statistics_optimized.json
```

The [standard-library script](../../04-computation/experiments/catalan_elliptic_20260921_statistics.py)
and [JSON](../../04-computation/experiments/catalan_elliptic_20260921_statistics.json)
use every integer through1,000,000 without filters. All tests use explicit
exceptions and remain active under `-O`. Independent direct proper-divisor
incidence and trial-factor profiles are checked through10,000. Eleven
prefixes, including X=1,2,3,4,8, verify every endpoint in (ST3), both
squarefree convolution paths, and the hyperbola formula. The sign
classification is checked on the full million; the three residue rows and
finite-prime CRT generating polynomials are counted exactly. Decimal
displays are labelled as such and are not certified asymptotic error bounds.

The exact finite results illustrate how slowly the equality families can
thin out; density zero does not mean negligible at a particular cutoff.

| X | Count D<0 | Count D=0 | Count D>0 | E_X D (exact terminating decimal) |
|---:|---:|---:|---:|---:|
| 10,000 | 5,639 | 2,175 | 2,186 | 0.2811 |
| 100,000 | 56,889 | 18,953 | 24,158 | 0.92155 |
| 1,000,000 | 573,459 | 165,862 | 260,679 | 1.617065 |

At X=1,000,000, the exact means of `(F,S,U)` are
`(11.970035,7.577760,2.775210)`. There are607,926 squarefree inputs.
These counts are independent checks at stated scales, not the proof of
the limits. Normal and optimized runs produce byte-identical JSON.

The stopping result is an exact comparison of spatial statistics:
squarefree exclusion determines the leading coefficient of E S and,
outside a zero-density set, the sign of F-S-U. The open transport question
is whether any arithmetic dynamics preserves a measure and a load-bearing
descent predicate compatible with these observables. An Euler product
alone does not provide that transport.
