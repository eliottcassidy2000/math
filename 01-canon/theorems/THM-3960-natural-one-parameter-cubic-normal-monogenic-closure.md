---
id: THM-3960
title: "Natural one-parameter cubic normal-monogenic closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For arbitrary
  C,E in k[t], the natural cubic F=T^3-3PT-(E+CP) is integral exactly when
  C^3+27E is nonzero. In that case its hypersurface has only finitely many
  singular points and is therefore normal. Its global monogenic different
  becomes a forbidden nonconstant unit on every putative Keller
  affine-plane open. If C^3+27E vanishes identically, F is a linear times
  quadratic product; the linear component is degree one and the quadratic
  affine-plane component is excluded by its ramification-line unit. For an
  irreducible hidden cubic E+C h^2-2h^3, the branch polynomial
  (E+CP)^2-4P^3 is additionally primitive and squarefree, giving an
  independent codimension-one maximal-order proof of normality. Together
  with THM-3956 and THM-3958 this closes the entire natural one-parameter
  globally monogenic family, not JC(2).
source: jc-degree6-one-place / post-THM-3956--3958 culmination, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-extra-debt-local, final_merge_referee,
  and audit_boundary_forest_3951, 2026-08-24). The three audits independently
  reconstructed the rational-root/integrality dichotomy and P-degree
  argument; finite-singular-locus R1+S2 normality, including the E=0
  doubled-branch hostile; the hidden-irreducible branch's primitivity,
  squarefreeness, and discriminant-index maximality; the normalization-form
  Zariski Main and global-different forbidden-unit contradictions; and both
  reduced components when B=0. They checked exhaustion by THM-3956/3958 and
  the exact function-field scope. Normal and optimized runs LF-normalize to
  the frozen 42-gate output; raw and semantic hashes and documentation checks
  pass, with no repair required.
depends_on:
  - THM-3956-split-hidden-cubic-integrality-and-repeated-root-trichotomy
  - THM-3958-one-hidden-root-principal-different-and-pure-power-boundary
related:
  - THM-3862-russell-finite-completion-nonmonogenic-branch-contract
  - THM-3847-two-place-cubic-deformation-monogenic-unit-debt
script: 04-computation/jc2_natural_one_parameter_cubic_normal_monogenic_thm3960.py
output: 05-knowledge/results/jc2_natural_one_parameter_cubic_normal_monogenic_thm3960.out
script_sha256: 86155a2c9427ff2ae4266a11d99b1c41d05b4c0fdbf88d5dc2d8a0bf18356e11
output_sha256: 853caab7ac647851bae242170257d6400c730a1e466554b2f84fbe019e37a784
semantic_sha256: e351375ce01977071e544aa926685b89d0d36d955a790a4ab8d1e2f6c23f5db1
hash_basis: raw LF bytes
---

# THM-3960 -- the full natural one-parameter cubic is closed

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Put

```text
R=k[P,t],
F(T)=T^3-3PT-(E+CP),                  C,E in k[t],          (1)
A=R[T]/(F),                           X=Spec A,
pi:X -> A2_(P,t),
G(h)=E+C h^2-2h^3.                                            (2)
```

The theorem closes this complete **natural one-parameter monogenic
grammar**. It does not assert that a hypothetical Keller normalization must
have the presentation `(1)`; escaping global monogenicity remains mandatory.

## 1. Exact integral-versus-reducible dichotomy

Define

```text
B=C^3+27E in k[t].                                           (3)
```

Suppose first that `F` is reducible over `Frac(R)`. As a cubic it has a root
there. Since `F` is monic, the root is integral over the normal UFD `R`, so
it lies in `R`. If its `P`-degree were positive, the terms `T^3` and `3PT`
at that root would have incompatible top `P`-degrees. Hence the root lies in
`k[t]`. Comparing the constant and linear `P` rows gives

```text
T=-C/3,                         B=0.                       (4)
```

Conversely, if `B=0`, put `x=C/3`, so `E=-x^3`. Then

```text
F=(T+x)(T^2-xT+x^2-3P).                                  (5)
```

Therefore

```text
F irreducible  iff  B!=0 as a polynomial.                 (6)
```

In the first case `A` is a domain and is finite free of rank three over
`R`, with basis `(1,T,T^2)`.

## 2. Every integral natural cubic is automatically normal

Assume `B!=0`. The hypersurface domain `A` is `S2`. Its partial derivatives
are

```text
F_T=3(T^2-P),
F_P=-(3T+C),
F_t=-(E'+C'P).                                             (7)
```

At a singular point, the first two equations force

```text
T=-C/3,                         P=C^2/9.                  (8)
```

Substitution in `F` gives

```text
F=-(C^3+27E)/27=-B/27.                                    (9)
```

Thus every singular point lies over a root of the nonzero one-variable
polynomial `B`. There are finitely many such parameters, `(8)` gives at most
one point over each, and the last equation in `(7)` only removes candidates.
The singular locus is finite, hence has codimension two. Therefore `A` is
`R1+S2` and Serre's criterion proves

```text
B!=0  implies  A is a normal domain.                       (10)
```

This argument deliberately includes nonreduced branch divisors. For example,
`E=0,C=t` gives

```text
F=T^3-3PT-tP,
H=P^2(t^2-4P),                                             (11)
```

whose branch has a doubled component, while `B=t^3!=0` and `(10)` still
proves that the surface itself is normal. Reducedness of the target branch
is sufficient for normality below, but not necessary.

## 3. Irreducible hidden cubic: the branch is primitive and squarefree

Now impose the sharper condition that the hidden cubic `G` in `(2)` is
irreducible over `k(t)`. Then

```text
E!=0,                             B!=0.                    (12)
```

Indeed `E=0` gives `G=h^2(C-2h)`, while `B=0` gives the exact factorization

```text
G=-(3h-C)^2(6h+C)/27.                                    (13)
```

The target branch polynomial and the order discriminant are

```text
H=(E+CP)^2-4P^3,
disc_T(F)=-27H.                                           (14)
```

As a cubic in `P`, its discriminant is

```text
disc_P(H)=-16E^3(C^3+27E)=-16E^3B.                       (15)
```

By `(12)`, this is nonzero in `k(t)`, so `H` is squarefree over `k(t)`.
Its leading `P^3` coefficient is the unit `-4`; hence its coefficients are
primitive and no vertical polynomial in `k[t]` can divide it. A repeated
irreducible factor in `k[t,P]` would therefore remain repeated after passage
to `k(t)[P]`, contradicting `(15)`. Thus

```text
H is reduced in k[t,P] and has no vertical component.      (16)
```

This gives an independent codimension-one proof of `(10)`. Let `p` be a
height-one prime of `R`. Away from `H`, the order `A_p` is etale. If `p`
divides `H`, reducedness gives

```text
ord_p(disc(A_p/R_p))=1.                                   (17)
```

Let `Abar_p` be the maximal `R_p`-order in `Frac(A)`. The DVR discriminant-
index formula is

```text
ord_p disc(A_p)=ord_p disc(Abar_p)
                 +2 length_(R_p)(Abar_p/A_p).             (18)
```

The rightmost term is a nonnegative even integer, while the maximal-order
discriminant exponent is nonnegative. Equation `(17)` forces the index
length to be zero. Hence `A_p=Abar_p`; every height-one local ring `A_q`
over `p` is a DVR. This proves `R1`, while the hypersurface gives `S2`.
Thus the squarefree-discriminant route independently proves normality in the
irreducible-hidden-cubic case, including the edge `C=0` and arbitrary common
factors of `C,E` (which cannot become vertical factors because of `-4P^3`).

## 4. The global different is a forbidden unit

Assume `B!=0`, so `(10)` identifies `X` with the finite normalization of
`A2_(P,t)` in its cubic function field. Suppose there were source
coordinates giving a same-function-field Keller map. Normalization-form
Zariski Main would identify the source with an open immersion

```text
j:U isomorphic to A2 -> X                                (19)
```

on which `pi` is etale. The monogenic presentation gives

```text
Omega_(A/R)=A dT/(delta dT),
delta=F_T=3(T^2-P).                                       (20)
```

Etaleness on `U` makes `delta|_U` a unit. But every unit on `A2` is a scalar.
The element `delta` cannot be scalar in `A`: its expression in the free
`R`-basis `(1,T,T^2)` has nonzero `T^2` coefficient `3`. This contradiction
proves

```text
B!=0  implies  no same-function-field nontrivial Keller A2 chart. (21)
```

This is a direct proof, not a circular invocation of THM-3862. It also
explains that theorem's relevance: a genuine Keller finite completion of
degree greater than one cannot be globally monogenic.

## 5. No lower-degree component survives when `B=0`

It remains to interpret the reducible boundary `(5)` honestly. Its reduced
normalization is the disjoint union

```text
T=-x(t):                         k[P,t],
T^2-xT+x^2-3P=0:                k[t,w],       w=2T-x,     (22)
```

where on the quadratic component

```text
P=x^2/4+w^2/12,
3(4P-x^2)=w^2,                  partial P/partial w=w/6.   (23)
```

The linear component maps isomorphically to the target and is the trivial
degree-one sheet. The quadratic component is an affine plane, but its
degree-two map ramifies on `w=0`. Any etale affine-plane open must delete
this line, after which `w` is a forbidden nonconstant unit. This includes
`x=0`, where `F=T(T^2-3P)`. Hence the reducible boundary supplies no
nontrivial Keller component.

Combining `(21)--(23)` proves the universal decision

```text
for every C,E in k[t], the natural family (1) has no
same-function-field nontrivial planar Keller chart.                 (24)
```

## 6. Relation to the structural split theorems and exact escape

THM-3956 records the much finer boundary triangle, repeated-root,
normalization, and class anatomy when `G` splits. THM-3958 records the exact
graph/residual principal-different and pure-power equality cases when `G`
has one root and an irreducible quadratic cofactor. The present theorem adds
the no-root squarefree branch calculation and shows that the simpler
normal-monogenic mechanism closes all factorization types at once.

The result does **not** prove JC(2). It proves that a counterexample cannot
have finite normalization equal to the natural globally monogenic order
`R[T]/(F)` in `(1)`. In fact, when `B!=0`, `(10)` says this order is already
the maximal finite order in its function field, so there is no proper
overorder escape with the same base and field. A viable construction must
leave this one-parameter linear-in-`P` depressed-cubic **function-field
grammar**, for example by starting with one of the genuinely nonmonogenic
S3 fields/orders studied elsewhere in the repo. Those constructions are not
overorders of `(1)` and remain relevant.
**QED.**
