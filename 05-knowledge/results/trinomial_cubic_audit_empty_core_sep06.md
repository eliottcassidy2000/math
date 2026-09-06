# Independent cubic trace-certificate audit

**Status: PROVED ANALYTICALLY + FINITE-EXACT identity certificate; audit PASS.**
This independently verifies the full uniform sign claim in the
[cubic-family proof](trinomial_cubic_empty_core_sep06.md): for every real
`g>=11`, the normalized anchored doubled residue has negative value at
all three first-row roots. Its first-return interpretation applies only
to integer `g>=11` with `gcd(g,21)=1` and nonzero complex coefficients
on `(-21,2g-21,3g-21)`.

The [independent source](../../04-computation/trinomial_cubic_audit_empty_core_sep06.py)
uses only Python's standard library, exact fractions, the companion
recurrence, and direct determinant expansion. It reads the four literal
candidate coefficient arrays and three denominators from the hash-pinned
primary source using `ast.literal_eval`; it does not import or execute
that producer, call a computer algebra system, or use its Newton-sum
implementation. These arrays are certificate data to verify, not
assumed identities. The [output](trinomial_cubic_audit_empty_core_sep06.out)
records 205 explicit gates and the separate identity/interpretation scopes.

## 1. Independent reduction and the finite identity argument

Put `u=g-7`, `K=(2g)_14/128`, and let `q_j` denote the coefficient of
`tau^j` in `Q_g/K`. Work in the quotient by

```text
P_g=72tau³+504u*tau²+84u(u-1)tau+u(u-1)(u-2).
```

Multiplication by `tau` uses the recurrence

```text
tau³=-7u*tau²-(7/6)u(u-1)tau-u(u-1)(u-2)/72.       (1)
```

The lower-carry inverse is

```text
tau^-1=-[72tau²+504u*tau+84u(u-1)]/[u(u-1)(u-2)].  (2)
```

All divisions in `(2)` are valid on `g>=11`. The apparently rational
parameter denominator cancels from the only inverse term because

```text
q_0/[u(u-1)(u-2)]
 =1024(2g-15)(2g-17)(2g-19)(2g-20)/21!,

q_j=128 product_(k=14)^(20-j)(2g-k)/[(21-3j)!(2j)!],
                                                    j=1,...,7. (3)
```

The empty product for `j=7` equals one. These identities follow directly
by cancelling the fourteen common falling-factorial factors and then
`(2g-14)(2g-16)(2g-18)=8u(u-1)(u-2)`.
The independent recurrence constructs the complete residue as

```text
R=q_0*tau^-1 + sum_(j=1)^7 q_j*tau^(j-1) mod P_g.
```

Assign weight one to both `g` and `tau`. Every term of the recurrence
replacement `(1)` has weight at most three. The inverse contribution
in `(2)-(3)` has weight at most six. For `j>=1`, the coefficient `q_j`
has degree `7-j`, so its term also has weight at most six. Therefore
`R=r_0+r_1 tau+r_2 tau²`, where

```text
deg_g r_k <= 6-k,      k=0,1,2.                    (4)
```

This proves polynomiality as well as degree bounds before any finite
specialization is used.

For the monic cubic the coefficient of `tau^(3-j)` has parameter
degree at most `j`. Its Newton identity, or induction using `(1)`,
gives `deg_g Tr(tau^n)<=n`. Consequently the weighted trace entry
`H_ij=Tr(R*tau^(i+j))` has degree at most `6+i+j`.
For its leading `k` by `k` determinant, each permutation term has
degree at most

```text
6k + sum_(i=0)^(k-1)i + sum_(i=0)^(k-1)i = k²+5k.
```

The three signed-minor bounds are thus `6,14,24`. The displayed
candidate factorizations have those same bounds. The independent
source evaluates both sides exactly at the **25 distinct integers
`g=11,...,35`**. Their differences have degree at most 24 and have
25 roots, so they vanish identically. This is an exact polynomial
identity certificate with proved degree bounds, not an empirical
parameter census. The residue coefficient identities, of degree at
most six, are verified by the same samples.

At each parameter the audit constructs multiplication traces by their
three diagonal entries in the basis `(1,tau,tau²)`, rather than by the
producer's Newton-sum method. It then expands each determinant directly.
The independently evaluated standard cubic discriminant has degree at
most six and agrees at all 25 points with
`15552(g-8)(g-7)² F3(g)`.

## 2. Sign conclusion and exact scope

All 38 literal coefficients of the four factors shifted by `g=s+11`
are positive, as are the three denominators. The verified identities
therefore prove all three signed leading minors strictly positive for
every real `g>=11`. The discriminant identity proves three distinct
real roots. Since `P_g` has positive coefficients and positive constant
term, all three roots are strictly negative.

The quotient trace matrix is
`V^T diag(R(lambda_1),R(lambda_2),R(lambda_3)) V`, with invertible real
Vandermonde matrix `V`. The three signed leading minors make it negative
definite, hence each `R(lambda_i)<0`. Multiplying by `K>0` gives the
actual sign `tau^-1 Q_g(tau)<0` at each first root. The carry is essential.
A trace and determinant without the middle minor would not suffice
for this three-root conclusion.

For integer `g` coprime to 21, the balance equation
`g(2y+3z)=21m` forces `g|m`. Its complete rows at masses `g` and `2g`
are respectively the four and eight displayed channels of the primary
proof. Their negative-charge counts are positive for `g>=11`.
Thus the noncancellation gives the exact first nonzero moment `g` or
`2g`, with both possibilities attained. The positive factor `K` scales
the `k`th signed minor by `K^k`; the three inherited unnormalized values
at `g=11` are reproduced exactly.

Four named admissible controls `g=11,13,16,22` independently enumerate
all multiplicity vectors by their original charge equation, including
weights and all earlier masses. The nonprimitive control `g=12` has
first support return four, confirming that the gcd hypothesis cannot
be dropped from the first-return interpretation. Formal mass-row and
sign identities remain valid there.

The separate [sector companion](trinomial_cubic_sector_empty_core_sep06.md)
proves the exact gauge and root overlap with THM-4440. At least two
first roots per support, and all three for admissible `g>=22`, lie
outside its real-rooted-core sector. That comparison is independently
proved and is not needed for the signed-minor identities.

The primary proof's complete rows, monomial anchor, inverse carry,
normalization, discriminant, inertia implication, primitive domain,
both attainment statements, and distinction from THM-4440 have all
passed the final text audit. No mathematical correction is needed.

## 3. Reproduction and frozen evidence

```sh
python3 -B 04-computation/trinomial_cubic_audit_empty_core_sep06.py
python3 -B -O 04-computation/trinomial_cubic_audit_empty_core_sep06.py
```

Normal and optimized output are byte-identical. Semantic digest:
`5df58bf57f9eaeee1eaf3e11669b6f49516f1670b9935077269b9fbd821c169c`.
Raw SHA-256:

```text
independent source
06dcd740dc0f3f133451ede0de2f258c22f34e7ed28393c77fa14baceeb2fa85
independent output
319af5bd209616715c379528d99587a7d537567a3cfee29848787326db348614
primary source pinned for literal certificate data
7ebdffa7cd529cc1d641dc319954c23a594cff5b4e98ae538e21ed18849b77f5
```

No general higher-channel negative-definiteness or two-rung theorem is
claimed. The next decisive object is a complete larger-dimensional
trace form, retaining every leading minor and the exact carry.
