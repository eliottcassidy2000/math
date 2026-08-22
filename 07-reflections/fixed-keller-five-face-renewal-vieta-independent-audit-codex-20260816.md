# The renewal hostile closes by taking the norm before reducing the root label

**Status: exact independent hostile audit; subsequently incorporated into
proved THM-3522.**  This reflection records the disjoint Vieta derivation and
exact companion.  It did not itself alter canonical status; the theorem owner
later combined it with the primary proof and the other independent audit to
promote THM-3522.

## Inheritance and the apparent obstruction

THM-3506 transports the first three faces of a complete packet `A(e,m)`
through

```text
Q=L^e N(P),       e'=7e-2m,       m'=3e-2m,
```

provided `Q` is polynomial.  THM-3513 obtains the remaining top-`z` and
minimum-`gamma` faces for the fixed `J -> G` step from a singleton that is
root-independent branch by branch.

For a general packet, put

```text
a=2e-4m/3,       b=2e-2m/3,
h=a-3b=-4e+2m/3.                                (1)
```

The exponent `h` need not be divisible by three: it is `-1018` for `G` and
`-6386` for `R_5`.  Thus `q^h` is not individually constant on the three
residual roots.  That is a genuine hostile to the branchwise proof, but not
to the norm.

## The shared endpoint

For any monomial of `P`, the packet inequalities give

```text
delta_6=gamma-k >= -10e+8m/3,
delta_8=gamma-3k >= -14e+4m.                      (2)
```

Equality in either bound requires both minimum `gamma` and maximum `k`.
The complete minimum-`gamma` face is

```text
z^e(27x^2z+y^3)^(e-2m/3),
```

and its maximum-`k` endpoint is exactly the complete singleton

```text
x^a z^b.                                           (3)
```

Therefore the same nonzero coefficient `c` controls both hybrid limits,
with no possible equal-weight cancellation.

## Vieta closes both norms

In THM-3513's top-`z` scaling, each branch contributes

```text
c(-C)^b q^h s^(10e-8m/3),
27A^2C q^3-2=0.
```

Although `q^h` is not root-independent,

```text
product(q)=2/(27A^2C).
```

After multiplication by
`L^e ~(27A^2C^2)^e s^(6e)`, the complete leading term is

```text
c_z' A^(10e-4m/3) C^(12e-8m/3),
```

which is precisely
`c_z' x^(2e'-4m'/3)z^(2e'-2m'/3)`.

In the minimum-`gamma` scaling, each branch contributes

```text
c(-C)^b q^h t^(-14e+4m),
Dq^3-3Bq-2=0,       D=27A^2C+B^3.
```

Now `product(q)=2/D`; multiplying by
`L^e ~(CD)^e t^(-8e)` gives

```text
c_gamma' C^(7e-2m)D^(5e-2m/3),
```

at weight `-50e+12m=-8e'+2m'`.  This is exactly the complete face

```text
c_gamma' z^e'(27x^2z+y^3)^(e'-2m'/3).                (4)
```

Polynomiality of `Q` remains load-bearing: it is what converts these generic
hybrid asymptotics into complete polynomial initial forms.  The argument does
not prove polynomiality of the next norm.

## Scalar coherence

Retaining the nonmonic Vieta factors gives more than support:

```text
c_gamma'=(-1)^(3b) 2^h c^3,
c_z'=27^(e'-2m'/3)c_gamma'.                           (5)
```

For `J -> G`, (5) independently reproduces THM-3513's exact values

```text
c_gamma(G)=3^513/2^117,
c_z(G)=3^1128/2^117.
```

It then predicts the previously unexpanded renewal scalars

```text
c_gamma(R_5)=3^3384/2^1369,
c_z(R_5)=3^7251/2^1369,

c_gamma(R_6)=3^21753/2^10493,
c_z(R_6)=3^46008/2^10493.                             (6)
```

Together with the already-proved polynomiality gates, the candidate lemma
would give complete packets

```text
A(1699,615) for R_5,
A(10663,3867) for R_6.
```

It does not prove that `R_7` is polynomial, that either polynomial is an
image equation or prime, a fifth nonproperness component, degree-`243`
separability, an unconditional all-level orbit, arbitrary-map closure,
`JC(2)`, or LRC.

## Exact companion

```text
04-computation/keller_five_face_renewal_vieta_independent_audit_20260816.py
05-knowledge/results/keller_five_face_renewal_vieta_independent_audit_20260816.out
```

The companion checks `15,250` admissible packets with `1<=e<=300`; in
`10,167` rows, `h` is nonzero modulo three, so branchwise closure is
unavailable.  Three split primes independently verify the product-of-powers
Vieta identity for the `J`, `G`, and `R_5` input states.  A hostile that drops
the nonmonic leading coefficient fails the required `A,C,D` exponents in
every nonzero packet row.  Normal and optimized outputs are byte-identical;
the semantic digest is

```text
7468b48f87fab27aacf77af954ea5b63189eb8ba22ec1402b9df57d6d1378375.
```
