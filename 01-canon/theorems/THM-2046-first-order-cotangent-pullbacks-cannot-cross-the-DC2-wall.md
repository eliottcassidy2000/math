---
id: THM-2046
title: "First-order cotangent pullbacks cannot cross the DC(2) wall"
status: >
  PROVED. In every rank, a Weyl-generator ansatz whose position images are
  multiplication polynomials and whose momentum images have differential order
  at most one forces the position map to have constant nonzero Jacobian. The
  first-order coefficients are its inverse-transpose Jacobian and the scalar
  terms form a flat connection. Consequently no such A2 ansatz can have first
  position x(a-b x^r q^s), by THM-2045. The exact A3 certificate therefore does
  not descend to A2 inside its own filtered pullback class.
source: codex-2026-07-21-DC2-filtered-pullback-wall
related:
  - THM-1300
  - THM-2044
  - THM-2045
  - HYP-8802
  - HYP-8803
script: 04-computation/dc2_ore_descent_codex_20260721.py
external:
  - https://github.com/techno-optimist/erdos-frontier-atlas/tree/main/certificates/dixmier-conjecture
---

# THM-2046 -- the first-order cotangent-pullback wall

Use the Weyl convention

```text
[d_i,x_j]=delta_ij,       [x_i,x_j]=[d_i,d_j]=0.
```

Let `A_n=D(C[x_1,...,x_n])`.  Suppose an assignment of its generators has

```text
Phi(x_j)=P_j(x),
Phi(d_i)=D_i=sum_k b_ik(x)d_k+c_i(x),                (1)
```

where every `P_j,b_ik,c_i` is a polynomial.  If the images satisfy the Weyl
relations, then, writing `J=(partial_k P_j)_(j,k)` and `B=(b_ik)`,

```text
B J^T=I,                    det J in C^*.             (2)
```

In particular `P=(P_1,...,P_n)` is a Keller map and necessarily

```text
B=(J^T)^(-1).                                       (3)
```

If `delta_i=sum_k b_ik partial_k`, then the `delta_i` commute, and the final
Weyl relations are exactly the flatness equations

```text
delta_i(c_j)-delta_j(c_i)=0.                         (4)
```

Thus (1) is not a larger quantum escape class surrounding the cotangent
pullback: its mixed commutators force the classical Keller map back into the
problem.

## Proof

Multiplication operators commute, and a scalar term has zero commutator with a
multiplication operator.  Hence

```text
[D_i,P_j]=sum_k b_ik partial_k(P_j).                 (5)
```

The `n^2` mixed Weyl relations say precisely `B J^T=I`.  Taking determinants
in the polynomial ring gives

```text
det(B) det(J)=1.                                    (6)
```

The only units of `C[x_1,...,x_n]` are the nonzero constants, so both
determinants are constant and nonzero.  This proves (2)-(3).

For completeness, (3) already forces the vector-field parts to commute.
Indeed `[delta_i,delta_j](P_k)=0` for every `k`, because each `delta_i(P_k)` is
constant.  If `v` is the coefficient row of `[delta_i,delta_j]`, this says
`v J^T=0`.  The matrix `J^T` is invertible over the polynomial ring, so `v=0`.
Expanding `[D_i,D_j]` now leaves only (4).  This proves the theorem.

## Rank-two corollary for the suspension coordinate

Let

```text
R=x(a-b x^r q^s),       a*b != 0, r,s positive.      (7)
```

There is no endomorphism in the class (1) of `A_2=D(C[x,q])` with
`Phi(x_1)=R`.  Otherwise its second multiplication-position image `Q(x,q)`
would obey

```text
R_x Q_q-R_q Q_x=c in C^*.                           (8)
```

After replacing `Q` by `Q/c`, this contradicts THM-2045.  This includes the
coordinate `R=x(2-3xq)` of THM-2044.

Notice the quantifiers.  The corollary excludes all *multiplication-position,
first-order-momentum* lifts with this `R`; it does not exclude a genuinely
nonfiltered Weyl endomorphism.  It therefore narrows DC(2), rather than proving
it.

## Why the exact A3 certificate stops here

The linked certificate was replayed exactly at external commit `a42a35a`.
It verifies a map

```text
Phi:A_3 -> A_3,
Phi(x_i)=F_i,
Phi(d_i)=sum_k ((JF)^T)^(-1)_ik d_k,                (9)
```

including all fifteen defining Weyl relations and a three-point collision.
That is precisely the class (1).  Its own terminal verdict is `DC_3` (and
higher), and its README explicitly says the object says nothing about `DC_1`
or `DC_2`.  Equation (2) explains structurally why deleting one pair is not an
indexing formality: the three multiplication positions of (9) would have to
become a planar Keller pair.

THM-2044 crosses the wall in the only currently viable direction.  Relative to
the source momenta `(p,z)`, its four classical output symbols have exact
differential orders

```text
ord(R,T,D,S)=(0,2,5,3).                             (10)
```

Indeed `ell=L+g` has order one, while the sheared polynomials have

```text
deg_ell(T)=2,       LC_ell(T)=x^3(3xq-2)/3,
deg_ell(S)=3,       LC_ell(S)=-x^4(3xq-2)/54,       (11)
```

and (4) of THM-2044 gives `D` a nonzero `ell^5` term.  Thus the Poisson
counterexample does not nearly belong to (1): the stabilization turns one of
the former base coordinates into a momentum direction, and the connection
then reaches order five.  HYP-8802 and HYP-8803 study exactly this nonfiltered
region.

## Scope guardrail

This theorem refutes the inference

```text
explicit A3 cotangent pullback + rank-two Poisson suspension
    => DC(2) counterexample by direct substitution.                (false)
```

It does not refute the target `DC(2) is false`.  It proves that any witness
extending the present `R` must use a momentum-dependent second position,
higher-order momentum images, or both.  THM-2044 already supplies the correct
principal symbols in that forbidden-to-filtered sector; polynomial closure of
their exact Weyl commutators remains the live gate.
