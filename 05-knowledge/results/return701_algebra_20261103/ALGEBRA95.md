# Exact rational classification of the persistent-port algebra

**PROVED by explicit rational matrix certificates and their consequences.** Local labels A1--A6.
The operation language and actual-source27-dimensional linear observer are
those of [JE4](../joint_port_geometry_20261103/JOINT_EFFECTS.md). Every generator
keeps the original second port; C3p/C5p use three live ports internally.
There is no assertion of closure under arbitrary native complement.

## A1 Complete algebra and its radical

Let A be the unital Q-algebra generated on that observer by C3p,C5p,N,J0,...,J3.
Then

```
dim_Q A = 95,
dim J=69, dim J^2=18, dim J^3=4, J^4=0,
A/J = Q^4 direct-sum M2(Q) direct-sum M3(Q)^2,
```

where J is the Jacobson radical. On the faithful27-dimensional source
module V, the dimensions of V,VJ,VJ^2,VJ^3,VJ^4 are

`27,21,7,2,0`.

This includes the extension data: [multiplication.json](multiplication.json)
is a complete multiplication table on a specified95-element rational
basis, with1647 nonzero basis products, the identity, radical powers and
all named generators. The concrete matrices give a faithful representation.
Thus the classification is not merely a list of radical dimensions and a
semisimple quotient; the table specifies the entire Q-algebra.

## A2 Rational upper and lower certificates

Use the physical capped-cover source space of dimension27. Row-reduce its
27-by46 basis B, choose pivot columns P so B[:,P]=I, and set
`A_L=B M_L[:,P]` for the46-state cap operators M_L. This coordinate choice
avoids the large fractions introduced by the previous future-query basis.
The exact conjugacy `A_L C=C A_L(old)` is checked for every generator.

The95 matrix words previously independent modulo1009 are reconstructed
over Q. Their flattened rows have rank95. Row reduction gives a basis R
with95 pivot matrix entries. Exact multiplication on the right by each
generator stays in span(R), proving an upper bound95. Because R is a word
basis containing the identity in its span, this is exactly A.

Take G_ij=tr(R_i R_j). Its exact rank is26, with69-dimensional kernel J0.
Cyclicity of trace and algebra closure imply this kernel is a two-sided
ideal: tr((ax)y)=tr(xya) and tr((xa)y)=tr(xay).
The literal module filtration shows J0^4=0 and J0^3 is nonzero. Therefore
J0 is a nilpotent ideal, without needing an imported trace-radical theorem.

## A3 The quotient splits over Q

The quotient A/J0 has dimension26 and its center has dimension7.
A specified central element has regular characteristic polynomial

`(x-2522)(x-472)(x-197)^9(x-14)^9(x+1)(x+7)^4(x+18)`.

Every central root is rational and distinct. Adapt a basis to the module
flag VJ^3 subset VJ^2 subset VJ subset V. Its four diagonal layers have
dimensions2,5,14,6. On these layers J0 acts as zero. Polynomial projectors
for the seven central roots select invariant submodules. Choosing a
smallest nonzero projected submodule for each root gives dimensions

`1,2,1,3,3,1,1`.

These supply explicit representations into four Q factors, one M2(Q),
and two M3(Q). The images of the95 basis matrices have rank26 and kill
exactly J0. The target product has dimension4+4+9+9=26, so the map is
onto the whole product and induces an isomorphism. In particular the
degree-three simple factors are split matrix algebras; rational center
alone would not have proved that. The quotient is semisimple, and the
nilpotent ideal J0 is therefore exactly the Jacobson radical J.

The sparse table is obtained from actual basis-matrix products at their
pivot coordinates. It reproduces all seven independently computed right
generator actions. Closure plus the faithful matrix embedding determines
associativity and every product; it is not an arbitrary fitted table.

## A4 All-word expansion has at most three layer changes

In the adapted basis every generator is block triangular on the four
radical layers. An entry of a word A1...AL is the sum over layer-index
paths. A path may remain on a layer for arbitrarily many steps, but it
can change layers strictly at most three times. Between those changes it
uses the corresponding diagonal matrix products; at each change it uses
the actual off-diagonal block of the generator at that time.

This is an exact noncommutative expansion, with word order and jump times
retained. For a weighted constant operator D+B, where D is block diagonal
and B strictly block triangular, its resolvent is the finite expression

`(I-z(D+B))^-1 = sum_(j=0)^3 [(I-zD)^-1 zB]^j (I-zD)^-1`.

Indeed any product with four strictly descending layer changes is zero,
even with block-diagonal factors between them. This gives rational weighted
recurrences retaining the transient terms. It does **not** mean all
information is forgotten after four ordinary operation steps.

## A5 Formal clocks do not manufacture legal words

[radical_clocks.py](radical_clocks.py) supplies an explicit M in J with
M^4=0, M^3 nonzero, entries in Z[1/2], and primitive integer numerator.
Its cube stays nonzero modulo3. For U=I+M,

`U^t=I+tM+binom(t,2)M^2+binom(t,3)M^3`.

The exact order modulo p^e is p^e for every prime p>=5 and is3^(e+1)
for p=3, for every e>=1. For p>=5, the truncated logarithm
`M-M^2/2+M^3/3` is p-primitive: it is M times an invertible polynomial
with constant term1. Exponential/logarithm, or the same binomial valuations,
give exact valuation e for U^(p^e)-I. At3 the coefficient
binom(3^t,3) has valuation t-1 while the other nonconstant coefficients
have valuation t, and M^3 has a3-unit entry. This proves the extra factor3.
The initial choice of M had cube zero modulo3; that choice would not
support this clock claim. The supplied matrix is part of the theorem.

These are **formal algebraic clocks**. C3p and C5p have rank26; N and every
J_m have rank23. Every nonempty legal positive word is singular. U is
invertible and nonidentity, so it cannot be any such word. The algebra
allows linear combinations of operations; the legal construction language
does not. A future claimed arithmetic controller must preserve this guard.

## A6 A uniform recurrence of order at most 24

For any fixed A in the operation algebra, let chi_i(x) be the characteristic
polynomial of its image in simple component i. Index the components by
the seven central eigenvalues in increasing order: -18,-7,-1,14,197,472,2522.
Their matrix sizes are 1,2,1,3,3,1,1. Then the explicit polynomial

`P_A(x)=chi_1 * chi_2^3 * chi_3^3 * chi_4^3 * chi_5 * chi_6 * chi_7`

has degree 24 and satisfies `P_A(A)=0` on the entire 27-dimensional
observer, including the radical contributions.

Proof: on each of the four semisimple layers, take the product of chi_i
over the different simple types present on that layer. Multiplicity of
one type does not require another factor, because the layer is a direct
sum of copies of its simple module. These four layer polynomials have
degrees 2,4,12,6. By Cayley--Hamilton each kills its own layer. Applying
them from the top layer down therefore kills V. Their evaluations commute
because all are polynomials in the same A, so their product annihilates
V regardless of order. The component incidences in semisimple.json give
the displayed exponents 1,3,3,3,1,1,1. This proves the formula for every A,
not just a generic linear combination of generators.

In particular, if A=sum w_L A_L is the weighted sum of the seven named
operators, every response `c_n=v A^n ell` obeys this constant-coefficient
recurrence. This is the weighted sum over length-n attachment words with
the specified initial source observer and terminal query. Coefficients
may be indeterminates: the argument works over their rational-function
field and gives polynomial identities after clearing denominators.
The order can drop for a particular weight or response; 24 is a uniform
upper bound, not a claim of minimality. Changing the generator at every
time still uses the ordered layer-path expansion A4 rather than this
constant-operator recurrence.

For the positive weights `(1,2,3,4,5,6,7)` on
`(C3p,C5p,N,J0,J1,J2,J3)`, the seven factors are, in the same order,

```text
x-436
x^2-1022x-7704
x-61
x^3-10810x^2+207633x+100913724
x^3-78426x^2+50513552x+17906890752
x+25
x+100
```

Raising them to exponents `(1,3,3,3,1,1,1)` gives the displayed
degree-24 annihilator. Exact Horner evaluation on the full matrix is
zero, so this example checks the radical terms too. The following
standalone check, executed in this session, reproduces that evaluation
from the pinned certificates (run as Python from the repository root):

```python
from pathlib import Path
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
p = Path('05-knowledge/results/return701_algebra_20261103')
a = json.loads((p/'algebra_exact.json').read_text())
b = json.loads((p/'semisimple.json').read_text())
w = {k: i+1 for i, k in enumerate(a['generators'])}
A = sum((w[k]*s.Matrix(v) for k, v in a['generators'].items()), s.zeros(27))
x = s.Symbol('x')
P = s.Integer(1)
for c, e in zip(b['components'], (1,3,3,3,1,1,1)):
    C = sum((w[k]*s.Matrix(v) for k, v in c['generator_matrices'].items()),
            s.zeros(c['dimension']))
    P *= C.charpoly(x).as_expr()**e
D = DM.from_Matrix(A).to_field()
I, Z = DM.eye(27, D.domain), DM.zeros((27,27), D.domain)
for c in s.Poly(P, x).all_coeffs():
    Z = Z*D + I.scalarmul(D.domain.convert(c))
if not Z.is_zero_matrix:
    raise RuntimeError('Degree-24 annihilator failed')
print('PASS exact degree', s.degree(P,x))
```
