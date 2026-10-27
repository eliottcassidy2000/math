# Exact Liouville and Möbius transforms

The definitions, divisor inversion and Dirichlet series are **CITED**
classical arithmetic: [DLMF 27.2](https://dlmf.nist.gov/27.2),
[27.4](https://dlmf.nist.gov/27.4), [27.5](https://dlmf.nist.gov/27.5).
The proofs below derive the applications explicitly. The square lift was
already used in [prime power bias and the Collatz collector](../prime_power_bias_collatz_collector_20261014.md),
section 1. No novelty is claimed for those identities or necklace inversion.

## A1 Prime powers explain the distinction

Write `Omega(n)=sum_p v_p(n)`, `omega(n)=#{p:p|n}` and
`lambda(n)=(-1)^Omega(n)`. Möbius satisfies `mu(n)=(-1)^omega(n)` on
squarefree n and zero otherwise. Both functions take value 1 at n=1.
Thus `mu(n)=lambda(n) mu(n)^2` pointwise. Liouville is completely
multiplicative; Möbius is multiplicative only for coprime arguments.

For Dirichlet convolution `(f*g)(n)=sum_(d|n) f(d)g(n/d)`, let 1 denote
the constant sequence and epsilon the unit supported at 1. Let s be the
indicator of positive squares and t(d²)=mu(d), with t zero off squares.
Then

```text
mu * 1 = epsilon                 lambda * mu² = epsilon
lambda * 1 = s                   lambda = mu * s
s * t = epsilon                  mu = lambda * t.
```

**Proof.** At a prime p, put z for its local variable. The generating
series for 1, mu, lambda, mu² and s are respectively
`1/(1-z)`, `1-z`, `1/(1+z)`, `1+z`, `1/(1-z²)`.
The identities follow by multiplication, including
`1/(1+z)=(1-z)/(1-z²)`. Multiplicativity completes the proof. QED.

This is an incidence-algebra distinction. The divisor interval of
`p^a` is a chain; its Möbius number from bottom to top is -1 for a=1
and zero for a>=2. Products of these chains give the arithmetic mu.
Liouville instead assigns the parity of the full exponent rank, including
every repeated step. A square contributes an even number of steps to
lambda and annihilates mu. They are different inverses: mu inverts 1,
whereas lambda inverts the squarefree indicator mu².

## A2 Dirichlet series and summatory functions

For `Re(s)>1`, absolute convergence gives

```text
sum mu(n)/n^s      = 1/zeta(s)
sum lambda(n)/n^s  = zeta(2s)/zeta(s)
sum mu(n)^2/n^s    = zeta(s)/zeta(2s).
```

These have different analytic behavior at the square scale. We use no
analytic continuation or global sign conjecture to prove the results here.
If `M(x)=sum_(n<=x) mu(n)` and `L(x)=sum_(n<=x) lambda(n)`, the finite
identities are

```text
L(x) = sum_(d<=sqrt(x)) M(floor(x/d²))
M(x) = sum_(d<=sqrt(x)) mu(d) L(floor(x/d²)).
```

For any fixed alpha>1/2, an `O(x^alpha)` bound for either sum implies
the same bound for the other, because `sum d^(-2alpha)` converges.
At alpha=1/2 the absolute-value argument loses a logarithm. In particular
the two families of bounds `O_epsilon(x^(1/2+epsilon))` are equivalent;
this restates a classical equivalence and proves neither bound.
The earlier prime-power-bias note already distinguishes Pólya-type bias
from Mertens-type cancellation. A finite prefix or the identities alone
does not settle their global sign behavior.

## A3 Lambert sums and ordinary Euler products

All identities in this section are formal power-series identities over Z;
each coefficient involves finitely many factors. The unweighted Lambert
transform reads the divisor convolution:

```text
sum_(n>=1) mu(n) q^n/(1-q^n)     = q
sum_(n>=1) lambda(n) q^n/(1-q^n) = sum_(d>=1) q^(d²).
```

The recovered [Euler ghost note](../../../07-reflections/euler-product-ghosts-and-zeta-side-cancellation-gates.md)
instead uses the weighted transform. Put

```text
P_b(q) = product_(n>=1) (1-q^n)^b(n)
G_b(n) = sum_(d|n) d b(d),      -q (log P_b)' = sum G_b(n)q^n.
```

Its finite observation that mu and lambda have dense ghosts has the
following all-order proof:

```text
G_mu(n) = product_(p|n) (1-p)
G_lambda(n) = product_(p^a || n) [1-(-p)^(a+1)]/(1+p).
```

Both are nonzero for every n>=1. Their signs are `mu(rad(n))` and
`lambda(n)` respectively, including the empty product at 1. This follows
by multiplicativity and the finite local geometric sum. The square lift
also persists at the product level:

```text
P_lambda(q) = product_(d>=1) P_mu(q^(d²))
P_mu(q) = product_(d>=1) P_lambda(q^(d²))^mu(d).
```

Indeed the exponent of `(1-q^n)` in the first product is
`sum_(d²|n) mu(n/d²)=lambda(n)`; the inverse is A1. This completes the
exact Dirichlet/product comparison requested by that reflection. It does
not identify ordinary coefficients with Dirichlet coefficients or infer
coefficient sparsity from ghost density. Through degree80, direct product
multiplication and the independent ghost recurrence agree. At degrees
11,13 the ordinary coefficients are `(1,1)` for P_mu and `(-2,2)` for
P_lambda; those are coefficients, not the values of mu or lambda.

## A4 What multiplication contexts recover

For an arithmetic observer f, define the response of m as the function
`k -> f(mk)` on positive integers. Equality of responses is the equivalence
relation that must survive multiplication by every context k.

**Proposition.** Lambda responses are equal iff the initial lambda values
are equal. Mu responses have one class containing all nonsquarefree
integers and a separate class for every squarefree integer.

**Proof.** Complete multiplicativity gives the lambda statement, with k=1
proving necessity. A nonsquarefree m has mu(mk)=0 for every k. For squarefree
m, the response at 1 is nonzero, and for every prime p,
`mu(mp)=0 iff p|m`. Hence the response recovers all prime divisors and
therefore m. Conversely an identical integer has an identical response. QED.

More generally, for squarefree m,k,
`mu(mk)=mu(m)mu(k)` when their prime supports are disjoint and is zero
when they overlap. This supplies an exact finite-dimensional contrast.
Fix a set P of r distinct primes and index subsets A,B by
`n_A=product_(p in A) p`. The two composition matrices are

```text
M[A,B] = mu(n_A n_B) = (-1)^(|A|+|B|) [A intersect B is empty]
L[A,B] = lambda(n_A n_B) = (-1)^(|A|+|B|).
```

L is a nonzero outer product, of rank one. M is the r-fold tensor power
of `[[1,-1],[-1,0]]`, with integral inverse
`M^(-1)[A,B]=(-1)^r [A union B=P]`. Thus it has rank `2^r` over Q and
over every field. The inverse follows from the local inverse
`[[0,-1],[-1,-1]]`. This is the Boolean disjointness coordinate retained
by mu; it is an exact operation-response classifier, not a classifier of
arbitrary tournament isomorphism types.

**Actual tournament contexts.** For order-join `T -> U`, every Hamiltonian
path traverses T before U, so `H(T -> U)=H(T)H(U)`. The recovered five-vertex
classes T11,T13 realize the prime alphabet `{11,13}` and all products
formed by joining them. In the subset order empty,{11},{13},{11,13},

```text
M = [ 1 -1 -1  1 ]       L has rank 1.
    [-1  0  1  0 ]
    [-1  1  0  0 ]
    [ 1  0  0  0 ]
```

Here `mu(H(T11 -> T11))=0` and `mu(H(T13 -> T11))=1`, while all
lambda responses of T11 and T13 agree. The first distinguishing displayed
join has ten vertices. No claim that ten is minimal among all possible
contexts or operations is needed. The alphabet must consist of attained
prime H-values; all primes are not assumed attainable.

For the symmetry orders 81 and27, the roles reverse: lambda separates,
while both mu responses are identically zero. Applying lambda to a group
order is compatible with finite exact sequences by order multiplication.
It is not the permutation sign of its automorphisms: all tournament
automorphisms are even permutations, including those in these examples.

## A5 Primitive periods and the Frobenius connection

Let W(n)=q^n count words of length n on q letters and F(n) those of least
period n. Every word repeats a unique primitive word, so
`W=1*F` and `F=mu*W`. Primitive necklaces number `F(n)/n`.
This is the classical formula inherited by the
[carousel and normal-basis dictionary](../module_shattering_lines_carousel_finite_field_20261026.md),
Theorem2: binary rotation corresponds to Frobenius after choosing a normal
basis; period equals field degree. The choice of basis is a necessary
coordinate, not a canonical identification of a tournament with a field element.

Using lambda in place of mu gives exactly

```text
(lambda * W)(n) = (lambda * 1 * F)(n)
                 = sum_(d²|n) F(n/d²).
```

The error is therefore the positive contribution of primitive words at
proper square-divisor scales. On squarefree n the replacement happens to
work; at nonsquarefree n it fails for q>=1. At n=4,q=2, it gives14 instead
of12 primitive words, and cannot count free four-orbits.

For q=2 and prime periods11,13, the primitive necklace counts are
`(2^11-2)/11=186` and `(2^13-2)/13=630`. Both arithmetic functions agree
on every divisor of these periods. At n=p² the correct primitive count
is `2^(p²)-2^p`; lambda adds F(1)=2. For p=11,13 its numerator therefore
has remainder2 modulo p², so cannot count necklaces of exact period p².
This is a decisive hostile to treating agreement at prime periods as an
interchangeable inversion theory. Direct binary-word enumeration verifies
all periods1..13; the prime-square tests use exact integer arithmetic.
