---
id: THM-3982
title: "Polynomial-shear submersion rational exactness and two-color image"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every
  height n>=2 and every polynomial h, A_h=x+h(x^n t) is an affine
  submersion. Its Hamiltonian has a nonzero rational invariant value exactly
  when deg(h)<=1. For nonconstant h this is the generic inverse-branch
  residue criterion delta^n(u)=0, which forces u to be a polynomial in h
  and hence h to be affine-linear. The source-polynomial and two-color
  completion invariant images are classified exactly in all three degree
  regimes. In particular no polynomial shear has a constant-Jacobian mate
  in B_n, although the constant shears do have the source mate t.
source: jc-degree6-one-place / post-THM-3978 polynomial-shear classification, 2026-08-24
audit: >
  PASS (root / jc-cohn3709, 2026-08-24). The audit independently rederived
  the Hamiltonian equation in (A,u), the inverse-branch residue formula,
  and the rational-exactness criterion, including the fact that infinity
  contributes no missing pole. It checked the delta-integration argument,
  its cancellation-free nonlinear leading term, the constant-shear Laurent
  intersection and exact q=n row, and the affine translation to THM-3978
  with the plus sign at the second color. Normal, optimized, and frozen
  outputs byte-match at CHECKS=217; all hashes agree.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
  - THM-3978-linear-seam-submersion-rational-mate-pole-obstruction
related:
  - THM-3974-height-tower-few-weight-darboux-support-obstruction
  - THM-3979-two-color-formal-cusp-darboux-lifting
script: 04-computation/jc2_polynomial_shear_submersion_rational_exactness_thm3982.py
output: 05-knowledge/results/jc2_polynomial_shear_submersion_rational_exactness_thm3982.out
script_sha256: 6685e693249d23d49dfd65b8b5abb5f39444ce338057f29e6184f6508b864983
output_sha256: 4982833deef5987838b2d1df618bcafb9e5222020c3904a88105353ca4f3df33
semantic_sha256: d6296d65a17dfae753d00943d1553503b0a7009d1ab6996e855d57dccac0b916
hash_basis: raw LF bytes
---

# THM-3982 -- polynomial shears are submersions but never completion coordinates

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over
an algebraically closed field `k` of characteristic zero. Fix `n>=2`, put

```text
u=x^n t,                  z=1+u,
p=x^-n u(u+1),            y=x^(-n-1)u^2(u+1),
B_n=k[x,u,p,y] subset k[x,t],                              (1)
```

and let `h in k[u]`. The polynomial shear

```text
A=A_h=x+h(u) in B_n                                      (2)
```

has no affine critical point for **every** `h`. Nevertheless it never has a
mate in `B_n` with constant nonzero Jacobian. More precisely, write

```text
I_rat(A)=J(A,k(x,t)) intersect k(A),
I_R(A)=J(A,R) intersect k[A]                             (3)
```

for the rational invariant image and the invariant image of the Hamiltonian
`J(A,-)` on a `k[A]`-algebra `R`, respectively. The complete table is as
follows.

If `h=beta` is constant, then

```text
I_rat(A)=k(A),
I_(k[x,t])(A)=k[A],
I_(B_n)(A)=(A-beta)^n k[A].                              (4)
```

If `h=alpha u+beta` with `alpha!=0`, then

```text
I_rat(A)=k(A),
I_(k[x,t])(A)=(A-beta)^(n-1) k[A],
I_(B_n)(A)=((A-beta)(A-beta+alpha))^(n-1) k[A].          (5)
```

Finally, if `deg(h)>=2`, then already

```text
I_rat(A)=0,                                              (6)
```

and hence both polynomial invariant images are zero. Thus `(4)--(6)`
separate three different debts:

```text
nonlinear h:       rational exactness fails;
affine alpha!=0:   rational exactness survives, polynomiality fails;
constant h:        the source mate survives, completion regularity fails. (7)
```

## 1. Every polynomial shear is a submersion

In the source coordinates,

```text
A_t=x^n h'(u),
A_x=1+n x^(n-1)t h'(u),
x A_x-n t A_t=x.                                       (8)
```

If both derivatives vanished, the last identity would give `x=0`. Since
`n>=2`, the second row then gives `A_x=1`, a contradiction. Equivalently,
when `x!=0`, a zero of `A_t` is a zero of `h'(u)` and again forces `A_x=1`.
This proves the submersion assertion uniformly, including the critical values
of the one-variable polynomial `h`.

## 2. Rational mates are generic inverse-branch quadratures

The rational fields satisfy

```text
k(x,t)=k(x,u)=k(A,u),                  x=A-h(u).          (9)
```

Because `dx wedge du=x^n dx wedge dt`, differentiation at fixed `A` gives

```text
J_(x,t)(A,Q)=x^n (partial Q/partial u)|A.                (10)
```

Consequently, for nonzero `r(A) in k(A)`, the equation

```text
J(A,Q)=r(A)                                             (11)
```

has a rational solution exactly when the one-form

```text
omega_A=du/(A-h(u))^n                                   (12)
```

is rationally exact over `k(A)`. Multiplication by `r(A)` neither creates
nor removes exactness.

Assume first that `h` is nonconstant. On `k(u)` put

```text
delta=(1/h'(u)) d/du,                  delta(h)=1.        (13)
```

Let `rho` be any generic inverse branch, so `h(rho)=A`. With local parameter
`w=h(u)-A`, one has `du=(delta u)dw`, and Taylor expansion in `w` gives

```text
Res_(u=rho) omega_A
  =(-1)^n (delta^n u)(rho)/(n-1)!.                       (14)
```

A rational differential on the projective line has a rational primitive if
and only if all of its geometric residues vanish. A primitive constructed
over a finite splitting field descends to `k(A)(u)` by averaging its Galois
conjugates. The generic branch embeddings are injective, so `(14)` vanishes
on every conjugate branch exactly when

```text
delta^n u=0 in k(u).                                    (15)
```

This is the exact residue gate; it is not merely a necessary condition.

## 3. Iterated integration forces an affine shear

The constants of the derivation `delta` on `k(u)` are exactly `k`, and
`delta(h)=1`. If `(15)` holds, repeated integration therefore gives

```text
u=P(h(u))                  for some P in k[s],
deg(P)<=n-1.                                               (16)
```

Indeed, `delta^(n-1)u` is constant; subtract its corresponding multiple of
`h^(n-1)`, and descend one derivative at a time. Since `u` is nonconstant,
`P` is nonconstant, and polynomial degrees in `(16)` give

```text
1=deg(P) deg(h).                                         (17)
```

Thus a nonconstant `h` passes the rational exactness gate if and only if
`deg(h)=1`.

There is also a cancellation-free check at infinity. If
`h=a_d u^d+lower` with `d>=2`, then

```text
delta^n u
 = [product_(j=1)^(n-1)(1-jd)]/(d a_d)^n
     u^(1-nd)(1+O(u^-1)),                                (18)
```

whose leading coefficient is nonzero in characteristic zero. This directly
rules out a hidden exceptional nonlinear polynomial.

Conversely, the two surviving rational primitives are explicit:

```text
h=beta:                 Q_0=t=u/(A-beta)^n,

h=alpha u+beta:         Q_0=(A-h(u))^(1-n)/(alpha(n-1))
                            =x^(1-n)/(alpha(n-1)).        (19)
```

They satisfy `J(A,Q_0)=1`. Multiplying them by an arbitrary element of
`k(A)` proves the first rows of `(4)--(5)`. If `deg(h)>=2` and a nonzero
`r(A)` occurred in the rational invariant image, `Q/r(A)` would be a
rational constant mate, contradicting `(14)--(18)`. This proves `(6)`.

## 4. Source-polynomial images

For `h=beta`, the mate `t` in `(19)` is already a source polynomial, and

```text
J(A,r(A)t)=r(A).                                        (20)
```

Hence the source image is all of `k[A]`.

For `h=alpha u+beta`, translate the invariant coordinate by putting

```text
A_0=A-beta=x+alpha u.                                   (21)
```

THM-3978 applies with `c=alpha` and gives

```text
J(A,k[x,t]) intersect k[A]
  =A_0^(n-1)k[A_0]=(A-beta)^(n-1)k[A].                  (22)
```

The generator is sharp. If

```text
V=A_0/x=1+alpha x^(n-1)t,
Q_src=(V^(n-1)-1)/(alpha(n-1)),                         (23)
```

then `Q_src in k[x,t]` and `J(A,Q_src)=A_0^(n-1)`.
Together with `(6)`, this proves the middle rows of the table.

## 5. The constant shear owes one full height

It remains to prove the new constant-shear completion row. Put `h=beta`, so
`A-beta=x`. If `J(A,Q)=r(A)`, integration of `(10)` gives

```text
Q=x^-n u r(beta+x)+H(x),                 H in k(x).      (24)
```

Expand

```text
r(beta+x)=sum_(j>=0)c_j x^j.                            (25)
```

Every element of `B_n` has a finite weight decomposition in
`k[x,x^-1,u]`. Thus if `(24)` lies in `B_n`, then `H` is Laurent in `x`;
here one uses
`k(x) intersect k[x,x^-1,u]=k[x,x^-1]`.
For `0<=j<n`, the weight `j-n=-q` has `1<=q<=n` and coefficient

```text
c_j u+h_(j-n).                                          (26)
```

The exact negative module of THM-3973 is

```text
(B_n)_(-q)
 =x^-q u^ceil(q/n)(u+1)^ceil(q/(n+1)) k[u]
 =x^-q u(u+1)k[u]                    (1<=q<=n).          (27)
```

The linear polynomial in `(26)` is divisible by `u(u+1)` only when both
coefficients vanish. Hence `c_0=...=c_(n-1)=0`, or equivalently

```text
(A-beta)^n divides r(A).                                (28)
```

Conversely, for `r(A)=(A-beta)^n s(A)`, the element

```text
Q=u s(A) in B_n                                         (29)
```

has Jacobian `r(A)`. This proves the last row of `(4)`, including its exact
exponent `n` rather than `n-1`.

## 6. The affine shear owes both colors

For `alpha!=0`, the translated coordinate `(21)` is exactly the linear seam
of THM-3978. Its two colors are

```text
u=0:                  A=beta,
u=-1 along D:         A=beta-alpha.                     (30)
```

The corresponding completion charges are therefore

```text
(A-beta)^(n-1),             (A-beta+alpha)^(n-1).       (31)
```

THM-3978 proves that these charges are jointly necessary and sufficient,
which gives the last row of `(5)`. An exact generator is obtained from
`(23)` by

```text
Q_B=(A_0+alpha)^(n-1)(V^(n-1)-1)/(alpha(n-1)),
J(A,Q_B)=(A_0(A_0+alpha))^(n-1).                        (32)
```

This also freezes the plus sign in the second factor.

## 7. Consequence and scope

All three completion image ideals in `(4)--(6)` omit `1`. Hence

```text
there is no Q in B_n with J(A_h,Q)=1
for any n>=2 and any polynomial h.                      (33)
```

The mechanism changes sharply with the degree of `h`: generic inverse-branch
residues, source polynomial poles, and the second completion color are three
different obstructions. The result closes the full polynomial-shear first-
coordinate family against companions of arbitrary support in `B_n`. It does
not classify arbitrary first coordinates, formal power-series shears,
nonpolynomial `h`, arbitrary Darboux pairs on `B_n`, or unrestricted `JC(2)`.
**QED.**

## Reproduction

```bash
python3 04-computation/jc2_polynomial_shear_submersion_rational_exactness_thm3982.py
python3 -O 04-computation/jc2_polynomial_shear_submersion_rational_exactness_thm3982.py
python3 agents/check_docs.py
```
