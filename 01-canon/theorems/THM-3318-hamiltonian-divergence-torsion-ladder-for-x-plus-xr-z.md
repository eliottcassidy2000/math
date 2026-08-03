---
id: THM-3318
title: "Hamiltonian divergence torsion ladder for x plus lambda x^r z"
status: >
  PROVED + GENERIC-SYMBOLIC + FINITE-EXACT HOSTILE-AUDITED.  For
  P=x+lambda*x^r*z over a characteristic-zero field, r>=2 and lambda!=0,
  the canonical divergence class is nonzero with exact K[P]-annihilator
  (P^r); its first P-multiple is the unit-response class, whose annihilator is
  (P^(r-1)).  Both become exact after deleting P=0, but their Laurent
  primitives have unavoidable poles, so no polynomial Jacobian mate exists.
source: root/creative-jacobian-lrc/2026-08-03
depends_on: []
related:
  - THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus
scripts:
  - 04-computation/jacobian_divergence_family_agent_exact.py
outputs:
  - 05-knowledge/results/jacobian_divergence_family_agent_exact.out
---

# THM-3318 -- Hamiltonian divergence torsion ladder for `x+lambda x^r z`

**PROVED + GENERIC-SYMBOLIC + FINITE-EXACT HOSTILE-AUDITED.**

## Statement

Let `K` be a characteristic-zero field, `R=K[x,z]`, `r>=2`, and
`lambda in K*`.  Put

```text
P=x+lambda x^r z,
D_P(f)=P_x f_z-P_z f_x,
C_P=R/D_P(R).                                             (1)
```

The quotient is a `K[P]`-module because `D_P(F(P)f)=F(P)D_P(f)`.  For any
polynomial Bezout row

```text
A P_x+B P_z=1,                                            (2)
```

the class

```text
mu(P)=[A_x+B_z] in C_P                                   (3)
```

is independent of the row.  For this family,

```text
mu=[lambda r(r+1)x^(r-2)z],              theta=[1],       (4)

Ann_(K[P])(mu)=(P^r),
Ann_(K[P])(theta)=(P^(r-1)),
P mu=-(r-1)theta.                                        (5)
```

Consequently

```text
K[P]mu ~= K[T]/(T^r),
K[P]theta=P K[P]mu ~= K[T]/(T^(r-1)).                    (6)
```

In particular `mu!=0`, `theta!=0`, and no `Q in R` satisfies
`Jac(P,Q)=1`.

## Canonical representative

One exact Bezout row is

```text
A=1-lambda r x^(r-1)z,
B=lambda r^2 x^(r-2)z^2.                                (7)
```

Direct multiplication and differentiation give

```text
A P_x+B P_z=1,
A_x+B_z=lambda r(r+1)x^(r-2)z,                          (8)
```

proving `(4)` as a representative.

If `(A',B')` is another row, the difference is

```text
(A'-A,B'-B)=(hP_z,-hP_x)                                (9)
```

for some `h in R`, since `(P_x,P_z)` is unimodular.  Its divergence is
`-D_P(h)`, so `(3)` is well-defined.  Moreover `mu=0` is equivalent to a
polynomial mate: a primitive of the divergence turns `(7)` into a closed
Bezout row, which integrates to `(Q_z,-Q_x)`; the converse follows from the
mate's closed row.

## Localization and complete response fibres

After inverting `x`,

```text
S=R[x^(-1)]=K[P,x,x^(-1)],
z=(P-x)/(lambda x^r),
D_P=-lambda x^r partial_x                              (10)
```

when `P` is held fixed.  A Laurent expansion in `x` shows

```text
ker(D_P:S->S)=K[P].                                     (11)
```

Define

```text
h_0=r z/x-lambda^(-1)x^(-r),
Q_0=x^(1-r)/(lambda(r-1)).                              (12)
```

Then

```text
D_P(h_0)=lambda r(r+1)x^(r-2)z,
D_P(Q_0)=1.                                             (13)
```

Equations `(11)--(13)` classify all localized solutions:

```text
{h in S:D_P(h)=A_x+B_z}=h_0+K[P],
{Q in S:D_P(Q)=1}=Q_0+K[P].                             (14)
```

The poles in `(12)` cannot be cancelled by a polynomial in `P`, so neither
fibre meets `R`.

## Exact annihilators

Put `v=lambda x^(r-1)z`, so `P=x(1+v)`.  The Laurent primitives yield
polynomial killing witnesses:

```text
P^(r-1)Q_0=(1/(lambda(r-1)))(1+v)^(r-1) in R,
P^r h_0=(1/lambda)(1+v)^r(rv-1) in R.                   (15)
```

Thus the powers displayed in `(5)` annihilate `theta` and `mu`.

They are sharp.  Let a nonzero `F(T) in K[T]` have `T`-adic order `m`, so
`F(T)=T^mG(T)` with `G(0)!=0`.  By `(14)`, a polynomial primitive of
`F(P)` or `F(P)(A_x+B_z)` would have to be, respectively,

```text
F(P)Q_0+H(P),             F(P)h_0+H(P),       H in K[T]. (16)
```

Along the divisor `x=0`, their first terms have exact orders

```text
m+1-r,                    m-r.                           (17)
```

If the displayed order is negative, its leading coefficient is a nonzero
multiple of `G(0)` and cannot be cancelled by `H(P)`.  Conversely `(15)`
handles every `m` at or above the threshold.  This proves both annihilator
equalities in `(5)`.

Finally, with

```text
g_0=(r-1)z+lambda r x^(r-1)z^2,                         (18)
```

one has

```text
P h_0+(r-1)Q_0=g_0,
D_P(g_0)=P(A_x+B_z)+(r-1).                              (19)
```

Passing to `C_P` gives the bridge `Pmu=-(r-1)theta` and hence `(6)`.

## Special fibre and Kummer sidecar

Both classes vanish after inverting `P`.  Indeed `P=x(1+v)`, so `x` becomes
a unit on `P!=0` and the primitives `(12)` apply.  The obstruction is
therefore supported at the special fibre as a polynomial-extension failure;
it is not a generic-fibre class that remains nonexact.

The Laurent mate has a separate finite-cover sidecar.  Writing `p=P`,
`q=Q_0`, and `w=x^(-1)` gives

```text
w^(r-1)=lambda(r-1)q.                                   (20)
```

After inverting `q`, this is a finite etale Kummer cover of degree `r-1`,
with geometric deck group `mu_(r-1)`.  For `r=2` the deck is trivial although
the polynomial extension still fails; therefore the deck is not the cause of
the obstruction.  At `r=3` it is a `C_2` only at the group level.  There is no
carrier map from this noncritical localization to THM-3306's exceptional
quadratic at a critical subresultant base locus.

## Verification and boundary

The exact companion checks `(7)--(8)`, `(12)--(13)`, `(15)`, `(19)`, and
`(20)` symbolically with formal `r,lambda`.  It tests `r=2,...,9` at four
nonzero scalar values, both bracket signs, sharp pole orders, killing powers,
and Kummer factor degrees.  For `r=2,...,8`, an exact bounded matrix bank
confirms that adjoining the divergence representative raises the image rank
by one while an actual image control does not.  The matrix bank is diagnostic;
the localization-and-pole proof is uniform in degree.  Normal and optimized
runs equal the frozen transcript.

```text
script sha256 (LF)  fc95ebad308c9d953c5a331a6c545c273907665cbc24bb0b32639775492f3731
output sha256 (LF)  004b82f379ae9726cb569f561349f89f23b41737945d1a581bbaa5e35675a22a
```

THM-2045 already proves no mate for the containing weighted family.  The new
content here is the exact canonical cokernel block, its annihilators, and its
one-step relation.  The family is not Keller: it has no polynomial mate and
therefore supplies neither a counterexample nor a new case of `JC(2)`.
