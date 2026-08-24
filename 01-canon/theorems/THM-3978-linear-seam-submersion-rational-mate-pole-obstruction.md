---
id: THM-3978
title: "Linear-seam submersion and rational-mate pole obstruction"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
  HOSTILE AUDIT REQUIRED. For every height n>=2 and c!=0, the two-weight
  element A_c=x+c(z-1) of the exact-volume completion B_n has no affine
  critical point. Nevertheless all of its rational constant-Jacobian mates
  have a compulsory pole, so none is polynomial. More sharply, its
  Hamiltonian image on the source polynomial ring meets k[A_c] in
  A_c^(n-1)k[A_c], while on B_n it meets k[A_c] in
  (A_c(A_c+c))^(n-1)k[A_c]. The second factor is the exact footprint of the
  second DPD color. This closes one named two-by-arbitrary first-coordinate
  family only; no unrestricted Darboux pair or JC(2) conclusion is claimed.
source: jc-cohn3709 / post-THM-3973 critical-free linear-seam lane, 2026-08-24
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
related:
  - THM-3974-height-tower-few-weight-darboux-support-obstruction
  - THM-3975-danielewski-one-arm-modification-cubic-control-and-hyperelliptic-no-mate
script: 04-computation/jc2_linear_seam_submersion_mate_ideal_thm3978.py
output: 05-knowledge/results/jc2_linear_seam_submersion_mate_ideal_thm3978.out
script_sha256: dcd7f328c96f410a1bbe8cf2f51dcb65df073a1cf564149c24ad8d265e095f9f
output_sha256: fb9a74e639b464baa6490d0b09e664baa2fbe32ff47a30d4f1edae05822c22aa
semantic_sha256: 402e6801c873e6e6abfb488f813ed82cc2b4250fd160558bd2b0e008970cf21b
hash_basis: raw LF bytes
---

# THM-3978 -- a submersion can still owe both completion colors

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
HOSTILE AUDIT REQUIRED.** Work over an algebraically closed field `k` of
characteristic zero. For `n>=2`, use the height tower

```text
u=x^n t,                 z=1+u,
p=x^-n u(u+1),           y=x^(-n-1)u^2(u+1),
B_n=k[x,u,p,y] subset k[x,t]                              (1)
```

of THM-3973. Fix `c in k^*` and put

```text
A=A_c=x+cu=x+c x^n t.                                    (2)
```

Then `A` has no critical point on the source affine plane. Despite that
positive control, it has no polynomial constant-Jacobian mate. The exact
statements are

```text
J(A,k[x,t]) intersect k[A] = A^(n-1)k[A],                (3)

J(A,B_n) intersect k[A]
  = (A(A+c))^(n-1)k[A].                                  (4)
```

Here the left side means the image of the Hamiltonian derivation
`J(A,-)` intersected with the displayed invariant subring. In particular,
neither image contains `1`.

Every rational mate of `A` is visible explicitly:

```text
J(A,Q)=1
iff
Q=x^(1-n)/(c(n-1))+H(A),             H in k(A).           (5)
```

The theorem therefore separates three obligations which should not be
conflated: absence of affine critical points, polynomial exactness along a
fibre, and regularity at both colors of the chosen completion. It closes the
named first coordinate `(2)` against companions of arbitrary support, not
all two-weight first coordinates and not `JC(2)`.

## 1. The first coordinate is a genuine submersion

In source coordinates,

```text
A_t=cx^n,
A_x=1+cnx^(n-1)t.                                       (6)
```

If `A_t=0`, then `x=0`; because `n>=2`, the second row of `(6)` is then
`A_x=1`. Thus `dA` never vanishes. This is not a near-miss caused by an
ordinary affine critical point.

The special zero fibre already records the subtler issue. Write

```text
V=A/x=1+c x^(n-1)t.                                     (7)
```

Then

```text
A^-1(0)=V(x) disjoint union V(V)
        isomorphic to A1 disjoint union G_m.             (8)
```

The two components are disjoint because `V=1` on `x=0`. The second component
is where the rational primitive below acquires its unavoidable pole.

## 2. Integrating the Jacobian equation exactly

Pass temporarily to the coordinates `(x,u)`. Since

```text
dx wedge du=x^n dx wedge dt,
```

one has, for every rational `Q`,

```text
J_(x,t)(A,Q)=x^n(Q_u-cQ_x).                              (9)
```

Let

```text
D=partial_u-c partial_x.
```

The change of variables `k(x,u)=k(A,x)` gives

```text
D(A)=0,                   D(x)=-c,                      (10)
ker(D)=k(A),
D[x^(1-n)/(c(n-1))]=x^-n.                               (11)
```

Equations `(9)--(11)` prove `(5)`. More generally, if

```text
J(A,Q)=r(A),                 r in k[A],                  (12)
```

then every rational solution is

```text
Q=r(A)x^(1-n)/(c(n-1))+H(A).                            (13)
```

This is a complete rational classification, not a bounded ansatz.

## 3. The first debt is exact polynomiality

Suppose the solution in `(13)` lies in `k[x,t]`, and set

```text
q_0(x)=Q(x,0) in k[x].                                  (14)
```

At `t=0` one has `A=x`, so `(13)` determines `H` and gives the unique
normal form

```text
Q=q_0(A)
  +r(A)(x^(1-n)-A^(1-n))/(c(n-1)).                      (15)
```

With `V` as in `(7)`, put `S(V)=1+V+...+V^(n-2)`. Then

```text
(x^(1-n)-A^(1-n))/(c(n-1))
  =t S(V)/((n-1)V^(n-1)).                               (16)
```

The polynomial `V=1+c x^(n-1)t` is prime. At its generic point, `x` and
`t` are units and `S(V)=1`. Hence the second term of `(15)` is regular
along `V=0` exactly when

```text
ord_(A=0)(r)>=n-1.                                      (17)
```

The regular summand `q_0(A)` cannot cancel a pole there. Conversely, when
`r=A^(n-1)s(A)`, the expression

```text
Q=s(A)(V^(n-1)-1)/(c(n-1))                              (18)
```

is polynomial and has Jacobian `r(A)`. This proves `(3)`. In particular,
the submersion `(2)` has rational mates but no polynomial constant mate.

## 4. The completion charges the second color

It remains to decide when the polynomial solution `(15)` belongs to `B_n`.
Once `(17)` holds, both `q_0(A)` and `r(A)A^(1-n)` are polynomials in `A`
and hence lie in `B_n`. Therefore

```text
Q in B_n
iff
x^(1-n)r(A) in B_n.                                     (19)
```

Taylor expansion at the weight-zero part `cu` gives

```text
x^(1-n)r(x+cu)
 =sum_(j>=0) x^(j+1-n) r^[j](cu)/j!.                    (20)
```

For `0<=j<=n-2`, put `q=n-1-j`. Then `1<=q<=n-1`, and the exact
homogeneous-piece formula of THM-3973 specializes to

```text
(B_n)_(-q)
 =x^-q u(u+1)k[u].                                      (21)
```

Thus every negative row in `(20)` is regular precisely when

```text
u(u+1) divides r^[j](cu),          0<=j<=n-2.           (22)
```

At `u=0` and `u=-1`, condition `(22)` says respectively

```text
ord_0(r)>=n-1,                 ord_(-c)(r)>=n-1.         (23)
```

These conditions are also sufficient, because every nonnegative row belongs
to `B_n`. Hence `(19)--(23)` prove the exact two-color ideal `(4)`.

The first factor in `(4)` is the source-polynomial pole payment already seen
in `(16)`. The second is new completion data: the color `u=-1`, which is the
boundary `z=0`, maps under `(2)` to the target value `A=-c`. The two DPD
colors have become two target factors with the same multiplicity.

## 5. Sharp positive controls

For

```text
r(A)=(A(A+c))^(n-1)h(A),
```

an explicit companion realizing `(4)` is

```text
Q=(A+c)^(n-1)h(A)(V^(n-1)-1)/(c(n-1)).                  (24)
```

Its apparent powers of `x` satisfy exactly the two-color module conditions
in `(21)`. At height two this becomes especially concrete:

```text
n=2,
A=x+cu,
Q=u+cxp,
J(A,Q)=A(A+c).                                          (25)
```

Thus the obstruction is sharp as an image-ideal calculation. It is not a
claim that the Hamiltonian derivation has no polynomial values, only that
every invariant value pays the exact two-color product and therefore cannot
be a nonzero constant.

## 6. Exact companion and scope

The assertion-free companion checks `(6)`, `(5)`, `(24)`, every negative
homogeneous row of `(24)`, both one-color near misses, and `(25)` for heights
`2<=n<=9`. Ordinary and optimized Python runs agree byte for byte:

```bash
python3 04-computation/jc2_linear_seam_submersion_mate_ideal_thm3978.py
python3 -O 04-computation/jc2_linear_seam_submersion_mate_ideal_thm3978.py
```

The proof of `(3)--(4)` is symbolic for every `n>=2`; the finite replay is an
independent hostile check, not a finite substitute for the proof. A general
two-weight first coordinate can have a different Hamiltonian foliation.
Unrestricted Darboux pairs on `B_n` and the planar Jacobian conjecture remain
open.
