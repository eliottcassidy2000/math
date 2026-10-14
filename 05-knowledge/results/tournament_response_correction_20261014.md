# Tournament response domains and Vandermonde typing: narrow correction

Status: PROVED algebraic correction + FINITE-EXACT minimality audit.
Date: 2026-10-14.
Scope: the domain claims surrounding THM-1950 and representation-language
claims surrounding THM-1805. Neither theorem's surviving exact identities are
retracted. The all-order Hamiltonian/Pfaffian bound remains OPEN.

## Inheritance and correction search

The closest proved mechanism is THM-1950's Schur-complement law, and the
canonical hostile is THM-4051's strong order-seven representative `299`, the
recorded minimum-response class. The retained coordinate is the signed
inverse response, not merely a determinant. The prior MISTAKE-491 repaired
a sampled order-seven extremum, not the super-multiplicativity assertion.
Targeted searches of the correction ledger and canonical files found no
existing repair of the claims corrected here. Both target files were
unchanged between audit snapshot `bd2d0d6660b3889179e226172a9572351bd122bc`
and publication-worktree base `200076ca71ea3abd4e6e6fb6bc8f1c2fd55092d2`.

## 1. Exact hostile and the two different minimality questions

Use `A[i][j]=1` iff `i` beats `j`, on vertices `0,...,6`. Mask `299` means
bit one orients `i->j` in lexicographic order on pairs `i<j`. Its adjacency
matrix is

```text
0110101
0000100
0100000
1110000
0011000
1111100
0111110
```

This tournament `C` is strongly connected. With `K=A-A^T`, exact calculations
give

```text
H(C)=31,  det(I+K)=384,  disc(C)=384/64=6,
(I+K)^(-1)1=(-1/6,1/2,1/3,1/6,1/3,-1/3,-1/6),
s(C)=2/3.
```

Both the sum and squared norm of the displayed inverse vector are `2/3`.
Let `v` be a singleton, with `disc(v)=s(v)=1`. The exact composition law
therefore gives

```text
disc(C join v)=6*(1+2/3)/2=5 < 6=disc(C)*disc(v),
H(C join v)=31,  s(C join v)=1.
```

These are precisely its two strongly connected components. Thus universal
super-multiplicativity over SCCs is false. For another control,
`disc(C join C)=26<36`, `s(C join C)=12/13`, and `H(C join C)=961`.

The independent exact labelled census through order six finds minimum
strong responses `3,3,2,1` at orders `3,4,5,6`, with strong counts
`2,24,544,22320`. Singletons have response one and no order-two tournament
is strong. Every SCC of a reducible tournament on at most seven vertices
therefore has response at least one. Expanding the multi-component
correction as an average of even-cardinality products shows it is at least
one. Hence **eight is the minimum order for failure relative to the product
over all SCCs**.

A distinct question allows two arbitrary, possibly reducible, blocks.
The directed triangle has `(H,disc,s)=(3,1,3)`. Consequently

```text
A=C3 join C3: (H,disc,s)=(9,5,3/5),
A join v:    (H,disc,s)=(9,4,1).
```

This is a two-block failure already at order seven, but not an SCC-product
failure: its SCCs are two triangles and a singleton. The census checks that
all tournaments through order five have `s>=1`, so no smaller pair of
nonempty arbitrary blocks can fail. **Seven is the minimum order for the
two-arbitrary-block claim**. These universes must not be conflated.

## 2. Exact repair and an always-defined multiplicative pair

For any two nonempty tournaments, retain the exact laws

```text
d12=d1*d2*(1+s1*s2)/2,
s12=(s1+s2)/(1+s1*s2).
```

The first is super-multiplicative exactly when `s1*s2>=1`, equal exactly
when the product is one, and sub-multiplicative when it is below one.
Nonnegativity alone supplies no choice of inequality direction.

Define

```text
Dplus=d*(1+s)/2 > 0,
Dminus=d*(1-s)/2.
```

Direct substitution proves `Dplus12=Dplus1*Dplus2` and
`Dminus12=Dminus1*Dminus2`. Conversely, `d=Dplus+Dminus` and
`d*s=Dplus-Dminus`. Thus this is an invertible coordinate change on the
response pair, not an analogy or a claimed new invariant theorem.
The second coordinate may be negative; it records the missing sign.

The global coordinate `q=Dminus/Dplus=(1-s)/(1+s)` satisfies
`q12=q1*q2`, including the absorbing case `s=1`, where `q=0`.
For the reduction's envelope,

```text
P=max(1,s)*d=Dplus+abs(Dminus).
```

Writing factor coordinates `(a,b)` and `(c,e)`, with `a,c>0`, gives

```text
P12=ac+abs(be) <= (a+abs(b))*(c+abs(e))=P1*P2.
```

The difference is `a*abs(e)+c*abs(b)`. Equality in this envelope law holds
exactly when both factor responses are one. Therefore the conditional
induction from `H(C)>=P(C)` on strong components to `H(T)>=P(T)>=disc(T)`
is unchanged. It did not require super-multiplicativity of `disc`.
The conjectured all-order equality characterization `H=disc iff transitive`
is not established by the weak strong-base inequality alone.

## 3. Rational composition is not a global real SL2 action

The map `y -> (s+y)/(1+s*y)` is represented by
`[[1,s],[s,1]]`, whose determinant is `1-s²`. At `s=1` the map is constant;
for `s>1` its determinant remains negative under every real nonzero scalar
rescaling. Real `arctanh(s)` is available only for `0<=s<1`, whereas regular
tournaments have `s=n`. Thus the universal real SL2/rapidity description is
incorrect. The signed coordinate `q` above is valid on the entire domain.

Nor does reducibility imply `s<1`: a singleton factor gives exactly one,
and `C3 join C` gives `11/9`. The exact boundary is

```text
s12-1=-(s1-1)*(s2-1)/(1+s1*s2).
```

## 4. THM-1805: preserve the identity, repair the representation

The all-order polynomial identity

```text
prod_(i<j)(x_i-x_j)=sum_T (-1)^(upsets(T))*x^(score(T))
```

and its score-wise cancellation are unaffected. Alternation is under the
permutation group on root variables. For `f(t)=a*prod_i(t-x_i)`,
`Disc(f)=a^(2n-2)*Delta(x)^2`; a quadratic has
`Disc(f)=a²*(x1-x2)²`, not `x1-x2`.

The correctly typed binary-form connection uses ordered linear factors
`ell_i=a_i X+b_i Y` and SL2 brackets `[ij]=a_i*b_j-a_j*b_i`. On the affine
slice `ell_i=X-x_iY`, the product of pair brackets is `Delta(x)`. This is
not an unspecified SLn action on the root coordinates or a single
n-vector determinant bracket.

The score vector is a multidegree. An actual one-parameter specialization
`x_i -> t^(w_i)*x_i` assigns charge `sum_i w_i*score(T)_i`. A transfer to the
Gaussian functional still needs polynomial coefficients, Wick weights,
and its target predicate preserved. The multigrading alone does not supply
that transfer. The classical stable cycle-choice involution remains cited;
this audit has not reconstructed it or claimed novelty for the identity.

## Reproduction and independent checks

```bash
python3 04-computation/experiments/tournament_response_correction_20261014.py
python3 -O 04-computation/experiments/tournament_response_correction_20261014.py
```

Both runs byte-match the frozen `.out`. The standalone verifier uses only
the standard library. Named witnesses use rational elimination checked
against independent integer Bareiss determinants; Hamiltonian subset DP
is independently checked against literal permutations through order eight.
The finite census computes response via the augmented determinant
`det([[I+K,1],[-1^T,1]])/det(I+K)-1`, avoiding the inverse-solve route.
It also checks the multiplicative pair, signed coordinate, and envelope on
positive, negative, and absorbing controls. No floating-point value is
used to decide an identity or inequality. No all-order bound is inferred
from the finite census.
