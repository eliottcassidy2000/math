---
id: THM-3922
title: "Affine-plane open boundary basis and class-group obstruction"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED. If a normal integral affine
  surface X contains a dense open U isomorphic to A2, then the prime
  divisorial components of X minus U form a Z-basis of Cl(X). In particular
  Cl(X) is torsion-free; if U is proper then the divisorial boundary is
  nonempty and Cl(X) has positive rank. Consequently the finite normal
  completion of a planar Keller map of generic degree greater than one
  cannot have finite class group or any class-group torsion. This is a
  necessary completion gate, not a proof of JC(2).
source: jc_degree6_one_place / post-THM-3920 boundary-localization lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_zero_debt_lift, 2026-08-23). The full
  Weil-localization sequence was rederived: restriction identifies both
  unit groups with k*, so the boundary map is injective, while Cl(A2)=0
  makes it surjective. Normal Hartogs plus affineness really forces U=X
  when there is no divisorial boundary. Finite normalization, Zariski Main,
  and triviality of connected finite etale covers of A2 verify the
  generic-degree-greater-than-one Keller corollary. The punctured-plane and
  third-Veronese hostiles correctly isolate the unit, normality, and
  actual-completion scope. No repair was required.
related:
  - THM-3578-zariski-main-boundary-rank-and-sheet-debt
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3924-decic-cubic-index-five-ramification-class-obstruction
---

# THM-3922 -- an affine-plane boundary is a class-group basis

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.** Let `k` be an algebraically
closed field of characteristic zero. Let `X` be a normal integral affine
surface of finite type over `k`, and suppose that it contains a dense open
subscheme

```text
j:U isomorphic to A2_k  ->  X.                                  (1)
```

Let

```text
D_1,...,D_s                                                   (2)
```

be the prime divisors that are irreducible components of `X minus U`.
Codimension-two components of the complement are deliberately absent from
this list. Then the boundary-class map is an isomorphism

```text
direct_sum_(i=1)^s Z[D_i]  --isomorphic-->  Cl(X).             (3)
```

Consequently:

1. `Cl(X)` is free abelian, and its rank is exactly `s`;
2. every `D_i` is primitive and the boundary classes have no nontrivial
   integral relation;
3. if `U` is a proper open subset of `X`, then `s>=1` and `Cl(X)` has
   positive rank; and
4. if `Cl(X)` is finite, or merely has nonzero torsion, then `X` cannot be
   a proper normal affine completion of an affine plane.

The same statements hold over an arbitrary field for which the displayed
open is the ordinary affine plane; algebraic closure and characteristic zero
enter only the Keller-map corollary below.

## 1. The localization sequence has neither kernel nor cokernel

Restriction gives an injection

```text
Gamma(X,O_X)^* -> Gamma(U,O_U)^*=k^*.                       (4)
```

The image already contains `k^*`, so

```text
Gamma(X,O_X)^*=Gamma(U,O_U)^*=k^*.                          (5)
```

For the normal integral scheme `X`, localization of Weil divisor class
groups along `j` gives the exact sequence

```text
Gamma(X,O_X)^* -> Gamma(U,O_U)^*
  -> direct_sum_i Z[D_i] -> Cl(X) -> Cl(U) -> 0.            (6)
```

The middle unit map records the boundary divisor of a rational function
that is invertible on `U`. Equation `(5)` makes its quotient zero. Meanwhile

```text
Cl(U)=Cl(A2_k)=0.                                           (7)
```

Thus the boundary map in `(6)` is both injective and surjective, proving
`(3)`. Notice that this is stronger than a rank count: it excludes torsion
and makes the actual prime boundary classes a basis.

Equivalently, any putative relation

```text
div_X(g)=sum_i n_i D_i                                    (8)
```

would make `g|_U` a unit. Since the only units on `A2` are scalars, every
`n_i` must vanish. Conversely, any Weil class on `X` restricts trivially to
the factorial open `U`, so it is represented by a boundary divisor. These
are the two exact directions of `(3)`.

## 2. A proper affine extension cannot hide in codimension two

Suppose `s=0`. Then `X minus U` has codimension at least two. Normal Hartogs,
or equivalently the height-one valuation intersection description of a
Krull domain, gives

```text
Gamma(X,O_X) --isomorphic--> Gamma(U,O_U)=k[x,y].           (9)
```

Both `X` and `U` are affine, and the restriction map in `(9)` is the ring
map induced by the open immersion `(1)`. Hence `(1)` is an isomorphism:

```text
s=0  implies  U=X isomorphic to A2.                        (10)
```

Therefore a proper inclusion in `(1)` has `s>=1`. Combining this with `(3)`
proves assertions 3 and 4.

## 3. Finite normal completions of Keller planes

Let

```text
f=(A,C):A2_(x,y) -> A2_(A,C)                              (11)
```

be a dominant planar Keller map, let `K=k(x,y)`, and let

```text
X=Spec S,
S=the integral closure of k[A,C] in K.                    (12)
```

The normalization is finite. Since `f` is etale and therefore quasi-finite,
normalization-form Zariski Main factors it as

```text
A2_(x,y) --open--> X --finite--> A2_(A,C).                (13)
```

If the first arrow were an isomorphism, `(13)` would be a connected finite
etale cover of `A2`. In characteristic zero such a cover is trivial, so the
generic degree would be one. Hence the finite completion attached to any
hypothetical Keller map of generic degree greater than one has a proper
affine-plane open. Applying Sections 1--2 yields the necessary invoice

```text
Cl(X)=Z^s,                s=# prime divisorial boundary >=1. (14)
```

In particular:

```text
finite Cl(X), any torsion in Cl(X), or a relation among the omitted
prime divisors  =>  no nontrivial Keller affine-plane atlas.        (15)
```

This complements THM-3578. That theorem lower-bounds the number of boundary
components by class-group rank for a general quasi-finite plane chart;
here the special source `A2`, including its scalar units and factoriality,
upgrades the inequality to the exact integral basis `(3)`.
Equivalently, THM-3811's deleted-divisor unit criterion supplies the
injective half of `(3)`, while `Cl(A2)=0` supplies the surjective half and
normal Hartogs supplies the properness conclusion.

## 4. Hostile boundaries and exact scope

The unit hypothesis in `(6)` is load-bearing, although it is automatic from
`U isomorphic to A2`. For the open

```text
G_m times A1 = D(x) subset A2,                              (16)
```

the omitted line is principal and `Cl(A2)=0`: the nonconstant unit `x` on
the open is exactly the kernel element that prevents `(3)`. Thus one cannot
replace the affine-plane open by an arbitrary factorial open.

The third-Veronese surface of THM-3801 and THM-3920 gives the sharp torsion
hostile:

```text
X=Spec k[t^3,t^2w,tw^2,w^3],               Cl(X)=Z/3.       (17)
```

It cannot contain a dense affine-plane open. In the explicit Chebyshev
attempt, deleting the total-ramification divisor `t=0` makes its torsion
relation visible as the forbidden unit `t^3`. This recovers the constant-row
failure in THM-3920 without using the special Chebyshev formula beyond the
identification of its finite completion.

Normality is essential to the Weil localization and Hartogs steps, and
affineness is essential when `(9)` is turned into the isomorphism `(10)`.
The obstruction must be applied to the **actual finite normal cubic
completion** `X`, not transferred from a quadratic resolvent, a projective
compactification, or another birational surface. A torsion class on one of
those related objects does not by itself trigger `(15)`.

Finally, `(14)` is only a necessary atlas invoice. A normal rational affine
surface with free positive class group does not automatically contain an
affine plane, and a basis of candidate boundary classes does not prove
effectivity, simultaneous contractibility, or etaleness of a proposed open.
Thus THM-3922 sharply prunes the counterexample search but does not prove
`JC(2)`. **QED.**
