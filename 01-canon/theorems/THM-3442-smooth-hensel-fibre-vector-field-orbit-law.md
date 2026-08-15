---
id: THM-3442
title: "Smooth Hensel-fibre vector-field orbit law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A
  congruence-to-identity automorphism of a smooth
  p-adic scheme has free equal-length orbits on every Hensel fibre exactly
  when its first-carry vector field has no zero on the relevant special-fibre
  points.  The first-carry criterion is sharp, including the p=2 depth
  boundary.
source: root-smooth-hensel-vector-field-orbits-2026-08-15
audit: independent global construction and proof reconstruction; nonlinear conjugate-translation control; normal/optimized/stored replay; hash, AST, ID, and documentation gates clean
depends_on: []
related:
  - THM-3439-near-identity-grassmannian-hensel-orbit-law
  - THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers
script: 04-computation/smooth_hensel_vector_field_orbits_thm3442.py
output: 05-knowledge/results/smooth_hensel_vector_field_orbits_thm3442.out
script_sha256: 26edd280a13f8a9e9cff84d0b1d480fc839dfc2134d0223b49ab4db5902a455c
output_sha256: 2e95fc184e1c2ec732f67d190a4333fd5f7e8d315e5dbc7b2f3e0e52e9be9abc
semantic_sha256: c72801a6c4fcb534b9e71ab0bb23da9163fc1051c45895b9bcd65241fd0a1c49
hash_basis: LF-normalized bytes
---

# THM-3442 -- smooth Hensel-fibre vector-field orbit law

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The proof and exact companion have passed an independent derivation and
replay audit.

## 1. Exact statement

Let `p` be prime and let `c>=1`.  Assume

```text
p is odd, or c>=2.                                      (1)
```

Let `X` be a smooth finite-type `Z_p`-scheme of pure relative dimension
`d>=1`, with `X(F_p)` nonempty.  Let `g` be a `Z_p`-automorphism whose
reduction modulo `p^c` is scheme-theoretically the identity.

The first infinitesimal displacement of the graph of `g` from the diagonal,
divided by `p^c` and reduced modulo `p`, is a vector field

```text
delta_g in H^0(X_(F_p),T_X).                            (2)
```

On an invariant affine formal chart this is the derivation

```text
delta_g(f)=((g^*f-f)/p^c) mod p.                        (3)
```

Flatness of `X/Z_p` makes division unique, and the product rule modulo
`p^(c+1)` makes `(3)` a derivation.  Equivalently, `(2)` is defined from the
first neighbourhood of the diagonal; this formulation is independent of a
chart.

For `a>=c`, reduction gives

```text
pi_(a,c):X(Z/p^a Z)->X(Z/p^c Z).                        (4)
```

If `delta_g` is nonzero at every point of `X(F_p)`, then every fibre of `(4)`
has

```text
fibre size       =p^(d(a-c)),
every orbit size =p^(a-c),
number of orbits =p^((d-1)(a-c)).                       (5)
```

For `a>c`, the action is transitive on every fibre exactly when `d=1`.

The zero condition is exact.  If `delta_g(xbar)=0`, choose any lift
`x_c` of `xbar`.  Every one of the `p^d` points over `x_c` at level `c+1` is
fixed by `g`.  Since smoothness lifts every `F_p` point, the following are
equivalent:

```text
delta_g has no zero on X(F_p);
g is free on every nonempty depth-c fibre at level c+1;
g is free on every nonempty depth-c fibre at every level a>c.              (6)
```

The statement is local on the smooth formal completion.  It therefore also
applies verbatim to a smooth `p`-adic formal scheme, locally topologically of
finite presentation and of pure relative dimension `d`, with finite
residue-point set and a congruence-to-identity formal automorphism.

## 2. Inheritance and connection ledger

THM-3439 is the linear homogeneous-space instance: on a Grassmannian, the
vector field induced by a matrix `E` vanishes at a plane exactly when that
plane is `E`-invariant.  The present proof isolates the actual mechanism—the
first-carry vector field—and also incorporates the surviving `p=2,c>=2`
range.  It does not use THM-3439 as a dependency.

The canonical hostile is a vector field with one zero: free orbits on other
fibres do not repair the fixed fibre above that zero.  The corrected near miss
is to say merely that a group element is “near identity.”  One needs an
actual automorphism, scheme-theoretic congruence modulo `p^c`, and the induced
derivation `(2)`.  The least-used sidecar is the orbit-bank label when
`d>1`.

| field | exact content |
|---|---|
| source | a smooth depth-`c` Hensel fibre with first-carry vector field |
| target | equal cyclic orbit packets |
| map | iteration of the congruence-to-identity automorphism |
| preserved | the base point modulo `p^c` and calibrated cyclic exponent |
| destroyed by one exponent | `d-1` transverse tangent coordinates |
| required sidecar | one of `p^((d-1)(a-c))` orbit labels |
| cheapest positive | translation on `A^d`, `x_1 -> x_1+p^c` |
| cheapest hostile | dilation on `A^1`, whose vector field `x partial_x` vanishes at zero |

## 3. Smooth fibre count

For each square-zero lifting step

```text
Z/p^(r+1) -> Z/p^r,
```

formal smoothness says that the lifts of a point form a torsor under its
`d`-dimensional tangent space over `F_p`.  There are therefore exactly `p^d`
lifts at every step.  Iterating from `c` to `a` proves the first formula in
`(5)` and shows that every special-fibre point lifts to every level.

Because `g=id mod p^c`, it preserves each fibre of `(4)`.

## 4. First-carry iteration

In local etale coordinates, write

```text
g^*=I+p^c D.                                             (7)
```

For `t>=0`, condition `(1)` and the binomial valuation give

```text
g^(p^t)*=I+p^(c+t)D_t,       D_t mod p=delta_g.          (8)
```

For odd `p`, the degree-two and higher binomial terms lie one layer deeper
starting at `c=1`.  For `p=2`, the same is true starting at `c=2`; after the
first squaring the depth only increases.  This proves `(8)` uniformly in the
declared range.

Take

```text
0<j<p^(a-c),       j=p^t q,       p does not divide q.   (9)
```

Then

```text
g^j*=I+q p^(c+t)D_t       mod p^(c+t+1).                (10)
```

If `x_a` were fixed by `g^j`, evaluation of `(10)` at `x_a`, followed by
reduction modulo `p`, would give

```text
q delta_g(xbar)=0
```

in the tangent space at its reduction.  Thus a shorter fixed orbit forces a
zero of the vector field.  Conversely, `(8)` with `t=a-c` gives

```text
g^(p^(a-c))=id mod p^a.                                 (11)
```

When `delta_g` has no zero, `(9)--(11)` prove that every orbit has exactly
`p^(a-c)` points.  Division into the fibre count proves `(5)`.

## 5. The first-lift translation and necessity

Fix `xbar in X(F_p)` and a lift `x_c`.  The `p^d` lifts to level `c+1` form
an affine torsor for `T_(xbar)X`.  In local etale coordinates, substituting a
lift `x_c+p^c v` into `(7)` shows that `g` acts on this torsor as

```text
v |-> v+delta_g(xbar).                                   (12)
```

Terms involving both the lift displacement and `p^cD` are divisible by
`p^(2c)` and vanish modulo `p^(c+1)`.  If the vector field vanishes, `(12)`
is the identity on the whole first lifting fibre.  If it is nonzero,
translation by it has no fixed point.  This proves the necessity and both
directions of `(6)`.

## 6. Examples and sharp boundaries

1. **Affine translation.**  On `A^d`, let

   ```text
   g(x_1,...,x_d)=(x_1+p^c,x_2,...,x_d).
   ```

   Then `delta_g=partial_(x_1)` is nowhere zero.  Each fibre splits exactly
   as `(5)`; for `d=1` it is one cyclic torsor.
2. **A mixed zero set.**  On `A^2`, the triangular automorphism

   ```text
   g(x,y)=(x+p^c y^2,y)
   ```

   has `delta_g=y^2 partial_x`.  Fibres reducing to `y!=0` are free, while
   every first lift above `y=0` is fixed.  This makes the pointwise
   quantifier in `(6)` sharp.
3. **The two-adic boundary.**  At `p=2,c=1`, take a `2x2` first carry with
   polynomial `T^2+T+1`.  On `P^1(Z/8)`, the square of the near-identity
   element is scalar, giving two 2-cycles rather than a 4-cycle.  At `c=2`,
   the extra quadratic carry is one layer deeper and the expected 4-cycle at
   level four returns.
4. **Smoothness.**  Without formal smoothness, lift fibres need not have
   constant cardinality or be tangent-space torsors.  No singular or
   nonreduced analogue is asserted.
5. **No physical transport.**  A Hensel orbit is not an LRC owner, current,
   or safe-time interval.  The theorem gives no LRC(14), `JC(2)`, or
   boundary-response consequence.

## 7. Grassmannian specialization

Let `X=Gr_k(Z_p^n)` and let `g` be induced by

```text
U=I+p^cE.
```

At a special-fibre plane `W`, the value of `delta_g` is the class of

```text
E|_W in Hom(W,F_p^n/W).
```

It vanishes exactly when `E(W) subseteq W`.  Thus `(5)--(6)` recover
THM-3439 for odd `p` and extend its conclusion to `p=2,c>=2`.  This is an
exact specialization, not a map from Grassmannian orbit labels to the
applications that motivated THM-3439.

## 8. Exact companion

The standard-library companion enumerates every fibre point and orbit for
affine translations in dimensions one through three, a nonlinear triangular
shear with both zero and nonzero vector-field fibres, multiplicative dilation
on `G_m` and its excluded zero, and projective `p=2` controls at `c=1,2`.
It checks the exact point count, orbit length/count, first-lift translation,
and normal/optimized equality using integer arithmetic only.

Reproduce with

```text
python3 -B 04-computation/smooth_hensel_vector_field_orbits_thm3442.py
python3 -B -O 04-computation/smooth_hensel_vector_field_orbits_thm3442.py
```

QED.
