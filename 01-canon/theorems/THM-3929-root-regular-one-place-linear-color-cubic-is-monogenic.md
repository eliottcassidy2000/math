---
id: THM-3929
title: "Root-regular one-place linear-color cubics force a scalar-unit index value"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the
  binary-cubic grammar Phi=a(A)U^3+C U^2V+c(A)UV^2+d(A)V^3, suppose an
  irreducible repeated-root component normalizes to A1, its affine repeated
  root is regular on that A1, and its projection to the A-line has degree
  three. Then the repeated root is integral over k[A], so its monic cubic
  forces a|c,d. Consequently a unit coefficient ideal forces a in k*, and
  Phi itself represents a scalar unit: the associated cubic order is
  monogenic. A separate centered-Mobius calculation also closes the first
  finite-pole seam: a degree-one raw root w=1/u forces q=0 and then the
  opposite endpoint coefficient is a scalar unit. This closes the
  root-regular and centered-Mobius one-place linear-color grammars, not
  homogeneous root maps of degree at least two. The latter are exactly the
  live THM-3927 escape mechanism.
source: jc_zero_debt_lift / post-THM-3927 one-place compression boundary, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-cohn3709, 2026-08-23). The audit
  reconstructed the degree-three integrality/minimal-polynomial argument,
  checked that root regularity is an explicit hypothesis rather than a
  consequence of one-place geometry, and verified the coordinate corollary.
  It independently derived the centered w=1/u seam, including primitivity,
  both Euclidean/resultant cancellations, the p^2=3q partial-cancellation
  hostile, and the q=0 scalar endpoint. The theorem remains explicitly
  scoped away from degree-at-least-two homogeneous root maps. The 39-gate
  companion byte-matches in normal and optimized mode, its frozen output and
  raw hashes pass, and documentation checks pass.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
related:
  - THM-3927-unit-ideal-rational-sextic-affine-address-cap-two-place-boundary
  - THM-3925-fivefold-conic-contact-torus-sextic-one-place-fold
script: 04-computation/jc2_root_regular_one_place_linear_color_cubic_thm3929.py
output: 05-knowledge/results/jc2_root_regular_one_place_linear_color_cubic_thm3929.out
script_sha256: 863c809aad0dc905a5143ab34444fab08dc5445a48a05f9431545f2f936a86d1
output_sha256: 983f42b62fb761c48d070479a7c8f6c4041ebc4a0cbbb82a6b1d5bbaa10c0f62
semantic_sha256: 38b0bfb97e7063f80f78b2bd14fa84092c4aeaf1fca55dd447a2b49966bfe0f6
hash_basis: raw LF bytes
---

# THM-3929 -- a regular repeated root makes one-place compression monogenic

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Put

```text
R=k[A,C]
```

and consider a binary cubic

```text
Phi(U,V)=a(A)U^3+C U^2V+c(A)UV^2+d(A)V^3.                 (1)
```

Here `a,c,d in k[A]` and `a` is nonzero. The intended application is when
`Phi` is the binary index form of a finite-flat cubic algebra, but the main
divisibility statement is purely about `(1)`.

## 1. The exact root-regular hypotheses

Use the affine root chart `V!=0`, with repeated-root coordinate `t=U/V`.
Let `I` be an irreducible component of

```text
Phi(t,1)=0,                 partial_t Phi(t,1)=0.           (2)
```

Assume the following.

1. `I` maps generically birationally onto an irreducible discriminant
   component `Gamma`.
2. The normalization is a one-place rational affine curve

   ```text
   I^nu = Gamma^nu = A1_s.                                 (3)
   ```

3. The chosen affine repeated root has no pole at a finite point of `(3)`:

   ```text
   t in k[s].                                               (4)
   ```

4. The generic polynomial obtained below is irreducible of degree three over
   `k(A)`. Equivalently, the projection `I^nu -> A1_A` has degree three.

Hypothesis `(4)` is the **root-regular** condition. In homogeneous language,
the repeated-root map

```text
I^nu -> P1_[U:V]                                           (5)
```

does not hit `[1:0]` over a finite point of `A1_s`. It is not a consequence
of rationality or of a one-place target curve. This distinction is the
load-bearing scope boundary of the theorem.

## 2. The incidence cubic

Write

```text
f(t)=a t^3+C t^2+c t+d.
```

The exact elimination identity

```text
2f-t f'=-a t^3+c t+2d                                    (6)
```

shows that the repeated root satisfies

```text
M_0(T)=a(A)T^3-c(A)T-2d(A)=0.                             (7)
```

By Hypothesis 4, `M_0/a` is the monic minimal polynomial of `t` over
`k(A)`, and

```text
[k(s):k(A)]=3.                                             (8)
```

Both target coordinates pull back regularly to `(3)`, so `A,C in k[s]`.
For a nonconstant polynomial `A(s)`, the field degree in `(8)` equals
`deg_s A`; hence

```text
deg_s A=3.                                                 (9)
```

In particular, after dividing by its nonzero leading coefficient, the
identity `A=A(s)` is a monic equation for `s` over `k[A]`. Thus

```text
k[s] is integral over k[A].                               (10)
```

The root-regular condition `(4)` now does the decisive work: `t` belongs to
the integral ring `k[s]`, so it is integral over `k[A]`.

## 3. Integrality turns coefficient primitivity into an actual unit value

Let `M(T)` be the monic minimal polynomial of `t` over `k(A)`. Since `t` is
integral over `k[A]`, every coefficient of `M` is integral over `k[A]`.
Those coefficients lie in `k(A)`, and `k[A]` is integrally closed. Therefore

```text
M(T) in k[A][T].                                          (11)
```

But `(7)` and irreducibility give

```text
M(T)=T^3-(c/a)T-(2d/a).                                   (12)
```

Comparison of `(11)` and `(12)` proves

```text
a divides c,                 a divides d       in k[A].    (13)
```

Consequently the coefficient ideal of `(1)` is not merely constrained; it
collapses exactly to

```text
(a,C,c,d)=(a,C) in k[A,C].                                (14)
```

If the coefficient ideal is the unit ideal, specializing `C=0` in `(14)`
shows `(a)=k[A]`. Hence

```text
a in k*.                                                   (15)
```

The binary cubic then represents the scalar unit `a` at `(U,V)=(1,0)`.
If `(1)` is a binary index form, THM-3801's index-form criterion gives

```text
the associated cubic algebra is globally monogenic.       (16)
```

Thus a hypothetical constant-unit degree-three Keller completion cannot
enter this grammar: THM-3801 requires its finite cubic completion to be
nonmonogenic.

## 4. Coordinate subcase: the missing quadratic and linear terms

There is a useful transparent specialization of the integrality argument.
Suppose the repeated root itself is the normalization coordinate, `t=s`.
Write

```text
A(s)=alpha s^3+beta s^2+gamma s+delta,      alpha!=0.      (17)
```

Substitute

```text
s^3=(A-beta s^2-gamma s-delta)/alpha
```

into `(7)`. The free `k[A]`-basis `(1,s,s^2)` gives successively

```text
beta=0,
c=-(gamma/alpha)a,
d=a(A-delta)/(2alpha).                                    (18)
```

The derivative equation in `(2)` now reads

```text
C(s)=-a(A(s))(3s^2-gamma/alpha)/(2s).                     (19)
```

If `gamma!=0`, polynomiality at `s=0` forces `a(delta)=0`; then `(14)` has
the common zero `(A,C)=(delta,0)` and is not the unit ideal. Under the unit
coefficient hypothesis, `(15)` therefore forces

```text
gamma=0,
A=alpha s^3+delta,
c=0,
d=a(A-delta)/(2alpha),
C=-(3a/2)s.                                               (20)
```

This is the literal `q=0`, constant-`a` monogenic fold. The coordinate
calculation is a corollary; the proof in Section 3 also covers regular
repeated roots such as `t=s^2` that generate `k(s)` together with `A=s^3`.

## 5. The first finite-pole seam also collapses

One finite-pole family can be handled without root regularity. Suppose the
normalization coordinate is `u`, and after affine normalization

```text
A=u^3+p u^2+q u+r,                 x=A-r.                  (21)
```

Put `w=1/u`. Its equation over `k(A)` is

```text
xw^3-qw^2-pw-1=0.                                        (22)
```

The trace-centered root

```text
t=w-q/(3x)                                                (23)
```

has the exact primitive incidence numerator

```text
27x^3t^3-9x(3px+q^2)t
 -(27x^2+9pqx+2q^3)=0.                                   (24)
```

Assume that `(23)` is the repeated root of a linear-color binary cubic
`(1)`, that its coefficient ideal is the unit ideal, and that the pulled-back
color `C(u)` is polynomial. These are exact normalized hypotheses; no claim
is made that every degree-one homogeneous root map can be put in this form
while preserving `(1)`.

When `q!=0`, the three coefficients in `(24)` are primitive in `k[x]`.
Therefore any polynomial incidence triple defining the same monic cubic is
`h(x)` times this triple for some `h in k[x]`. Its common coefficient content
is `h`, so the unit coefficient ideal forces `h in k*`.

The derivative equation gives, up to this harmless scalar,

```text
C(u)= -27u^2(u^2+pu+q)^2 E(u)
       / (2(3u^2+3pu+2q)),                                (25)

E(u)=9u^4+21pu^3+(12p^2+12q)u^2+15pqu+5q^2.              (26)
```

Set `D=3u^2+3pu+2q`. Exact Euclidean division and a resultant give

```text
E mod D = q(pu+q),
Res_u(D,u(u^2+pu+q))=2q^3.                                (27)
```

Thus for `q!=0`, the denominator `D` is coprime to the displayed prefactor
in `(25)`, but cannot divide `E`; its remainder is a nonzero polynomial of
degree at most one. At the apparent exceptional seam `p^2=3q`,

```text
Res_u(D,E)=-27q^3(p^2-3q),                                (28)
```

so one root can coincide, but `(27)` shows the other still does not cancel.
Polynomiality of `C(u)` therefore forces

```text
q=0.                                                       (29)
```

For `q=0`, divide `(24)` by its content `27x^2`. The primitive incidence is

```text
xt^3-pt-1=0.                                               (30)
```

Hence the actual coefficient triple is

```text
(a,c,d)=h(x)(x,p,1/2).                                    (31)
```

It is primitive, so the same coefficient-ideal argument forces `h in k*`.
Now `d=h/2` is a scalar unit and `Phi(0,1)=d`; THM-3801 again makes the cubic
order monogenic. Directly,

```text
t=1/u,        x=u^2(u+p),        C=-h u(3u+4p)/2,          (32)
```

so the surviving `q=0` row really is polynomial.

This closes the centered-Mobius finite-pole row. It does not close root maps
of degree at least two.

## 6. Why the omitted higher homogeneous-root case is real

If `(5)` hits `[1:0]` over a finite normalization point, the affine ratio
`t=U/V` has a finite pole and no longer belongs to `k[s]`. The implication

```text
t regular on A1_s  =>  t integral over k[A]               (33)
```

then fails at its first step, so `(13)` cannot be inferred.

The exact THM-3927 sidecar illustrates this mechanism. With normalization
coordinate `u`, its repeated root contains

```text
t=-4(u+2)/((u-1)(u+1)),                                   (34)
```

and hence has finite root-at-infinity points `u=1,-1`. Multiplying the
color coordinate by `A` cancels a different finite pole and produces a
polynomial one-place parametrization, but it does not make `(34)` integral;
the scaled coefficient ideal instead acquires the common factor address
`(A,B)=(0,0)`. The companion checks these identities directly.

Therefore this theorem does **not** close arbitrary `b=C` cubics, homogeneous
repeated-root maps of degree at least two, or JC(2). The next live object is
the finite pole divisor of `(5)`: one must either eliminate it without losing
the unit coefficient ideal or prove that its divisor payments force
monogenicity by a different argument.

## Reproduction

```bash
python3 04-computation/jc2_root_regular_one_place_linear_color_cubic_thm3929.py
python3 -O 04-computation/jc2_root_regular_one_place_linear_color_cubic_thm3929.py
```

Both runs must byte-match
`05-knowledge/results/jc2_root_regular_one_place_linear_color_cubic_thm3929.out`.
