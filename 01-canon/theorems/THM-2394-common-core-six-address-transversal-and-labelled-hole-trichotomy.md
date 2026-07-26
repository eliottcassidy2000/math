---
id: THM-2394
title: "Common-core six-address transversal and labelled-hole trichotomy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING. Assume
  THM-2393's sole no-clean residual M=1 and
  (C_1,C_2,c_1,c_2)=h(1,13,13,169). On every generic fibre where
  q_*,C_3,c_3 are safe, the guard and four lower q words form an exact
  six-address transversal: their multiplicity K is 0 at one root and 1
  at the other six. Writing A,B,C for the singleton roots of
  D_h,D_(13h),D_(169h), the hole is either B or A=C. There are exactly
  three labelled cases. The hole is B on base mass at least 3972/8281,
  and the generic type with hole B and distinct double root C has mass
  at least 3335/8281. Hence one fixed ordered hole/double pair has mass
  at least 3335/347802 and carries every nonzero septimal colour with
  exact charged phase. The address transition is carry-corrected:
  b(y)=kappa_13(y)-a({13y}) and
  c(y)=a({169y})-kappa_169(y)=kappa_13(y)-b({13y})
  modulo seven. The carry-free rule is false already at y=1/10. This is
  a labelled carry reduction, not a row exclusion, ledger decrement, or
  proof of LRC(14).
source: codex-2026-07-26-common-core-transversal
depends_on:
  - THM-2388-thirteen-root-multiplicity-reflection-and-blocker-caged-toothpick-law
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
  - THM-2393-c3-safe-double-fibre-capacity-and-common-core-residual
related:
  - THM-2232-same-core-signed-eigen-markov-dual-exclusion
  - THM-2257-depth-three-common-core-169-image-sieve-exclusion
  - THM-2362-target-shift-successor-role-jet-and-inverse-root-sheet-anchor
script: 04-computation/lrc14_common_core_transversal_thm2394.py
output: 05-knowledge/results/lrc14_common_core_transversal_thm2394.out
script_sha256: 74f75a65f88c0b604bbbacf34f906bdf883a3eb80e1f271ca3e918e4bb767cf6
output_sha256: 90da48855d57c5f3cecc6cf94bd16acbb6f27d54b4120242d29e7d7a944e5d3f
hash_basis: working-tree bytes (LF)
---

# THM-2394 -- the common-core fibre is a labelled dipole

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

THM-2393 leaves one exact boundary:

```text
M=1,

(C_1,C_2,c_1,c_2)=h(1,13,13,169),       gcd(h,91)=1.       (1)
```

This theorem retains the labels which an unlabelled capacity count
discarded. The result is stronger than “one low-cage address exists”:
the six guard/lower-`q` incidences form a transversal, and on a large
subcarrier the missing address and the unique double address are a
fixed ordered pair.

The price of following that pair under multiplication by thirteen is
an explicit base-digit carry. Omitting it gives a false transition law.

## 1. The high-safe seven-root fibre

Retain THM-2393's notation

```text
K=1_(E_H)+sum_(q_i != q_*) 1_(D_(q_i)),

X_0
 =D_(q_*)^c intersection D_(C_3)^c intersection D_(c_3)^c.
                                                               (2)
```

Here `K` contains the guard word and the four lower ordinary words; the
top word is absent on `X_0`. Equation (11) of THM-2393 gives

```text
mu(X_0)=396/637.                                      (3)
```

Disintegrate over

```text
x_r=(y+r)/7,                         r in F_7.         (4)
```

Every high condition in `X_0` is constant on (4). Away from the
finitely many strict endpoints, put

```text
k_r=K(x_r).                                           (5)
```

The guard has two roots and each lower ordinary word has one, so

```text
sum_r k_r=2+4=6.                                      (6)
```

In the common-core gauge (1), let `A,B,C` be the singleton root words

```text
A_r=1_(D_h)(x_r),
B_r=1_(D_(13h))(x_r),
C_r=1_(D_(169h))(x_r).                               (7)
```

The same word `B` has two roles:

```text
quotient blocker C_2=13h,
actual blocker c_1=13h.                              (8)
```

Assume the no-clean alternative `delta=0`. On almost every fibre in
`X_0`, the scalar cover, the THM-2388 collision cage, and the definition
of a clean hole give respectively

```text
k_r+B_r+C_r>=1,                                      (9a)

k_r>=2  implies  A_r+B_r>=1,                         (9b)

k_r=0   implies  A_r+B_r>=1.                         (9c)
```

Equation (9c) uses `C_3`-safety: the quotient-blocker cage reduces from
`D_(C_1) union D_(C_2) union D_(C_3)` to `A union B`.

## 2. Six-address transversal lemma

> **Lemma.** Let `k_r` be seven nonnegative integers of sum six, and let
> `A,B,C` be singleton words on `F_7`. If (9a)--(9c) hold, then
>
> ```text
> k_r in {0,1},
> ```
>
> with exactly one hole `d`. Moreover precisely one of the following
> labelled cases holds:
>
> ```text
> I.    d=B=C;
>
> II.   d=B,       C!=B;
>
> III.  d=A=C,     B!=d.                              (10)
> ```

**Proof.** Put

```text
p_r=k_r+B_r+C_r.
```

By (6) and the singleton counts,

```text
sum_r p_r=8.
```

Equation (9a) says all seven `p_r` are positive. Hence

```text
sum_r(p_r-1)=1,                                      (11)
```

so the full word has exactly one double root and no higher
multiplicity.

Suppose `k` has a collision. Equation (11) then forces one root with
`k_r=2`, no actual blocker there, two `k`-holes, and no other
collision. By (9b) the collision root is `A` rather than `B`. The two
holes must be the distinct singleton roots `B,C`. Equation (9c) forces
the `C`-hole also to be `A`. Thus the same root is both the `k`
collision `A` and the actual blocker `C`, contradicting (11).

Therefore `k_r<=1`. Its sum is six, so there is one hole `d`.
Equations (9a) and (9c) give

```text
d in (B union C) intersection (A union B)
  =B union (A intersection C).                       (12)
```

If `d=B`, split according to `B=C` or `B!=C`, giving I and II. If
`d!=B`, equation (12) gives III. These cases are disjoint and complete.
QED.

The exact abstract enumeration contains

```text
type I:       49 labelled states,
type II:     294 labelled states,
type III:     42 labelled states,
total:       385.                                    (13)
```

Its only seven `k`-words are the complements of the seven singleton
holes. Thus the lemma is not an averaging statement.

## 3. Type II occupies a large carrier

Write `a(y),b(y),c(y)` for the root addresses of `A,B,C`. THM-2393's
common-core occupancy law has a direct labelled interpretation. Since

```text
U_(1,1)=B union (A intersection C),                  (14)
```

its two-root event is exactly

```text
a(y)=c(y)!=b(y).
```

Equation (37) of THM-2393 therefore gives

```text
mu{a=c!=b}=24/169.                                   (15)
```

The event `b=c` is also exact. Root disintegration and the common
dilation identity give

```text
mu_y{b=c}
 =7 mu_x(D_(13h) intersection D_(169h))
 =7 mu_x(D_h intersection D_(13h))
 =7/91
 =1/13.                                              (16)
```

Case III is contained in (15). Hence on a subset of `X_0` of mass at
least

```text
396/637-24/169
 =3972/8281,                                         (17)
```

the unique `K`-hole is the middle root `B`. Case I is contained in
`{b=c}`. Consequently case II alone has base mass at least

```text
3972/8281-1/13
 =3335/8281.                                         (18)
```

There are seven possible hole addresses, so one fixed middle-hole
address occurs on mass at least

```text
3972/(7*8281)=3972/57967.                            (19)
```

There are `7*6=42` ordered distinct pairs `(b,c)`. Thus one fixed
type-II hole/double pair occurs on a base cell `Y_(b,c)` with

```text
rho:=mu(Y_(b,c))
 >=3335/(42*8281)
 =3335/347802.                                       (20)
```

No independence of `X_0` from the address events is used in
(17)--(20); only subtraction of their exact global masses is used.

## 4. Exact septimal dipole and phase

On a type-II cell, the lemma says pointwise

```text
K_y=1-1_({b}),

K_y+B_y+C_y=1+1_({c}),                b!=c.           (21)
```

Let `zeta_7=exp(2 pi i/7)` and use the normalized root DFT. For every
nonzero `ell in F_7`,

```text
Khat_y(ell)=-(1/7)zeta_7^(-ell b),

Jhat_y(ell):=
  Fourier[K_y+B_y+C_y-1](ell)
  =(1/7)zeta_7^(-ell c).                             (22)
```

Therefore every nonzero septimal colour survives, and the charged
product retains its exact orientation:

```text
Khat_y(ell) conjugate(Jhat_y(ell))
 =-(1/49)zeta_7^(ell(c-b)).                          (23)
```

On the fixed cell (20), the integrated coefficient and charged-product
magnitudes are at least

```text
rho/7 >=3335/2434614,

rho/49>=3335/17042298.                               (24)
```

Pointwise Parseval and fourth-power sums are

```text
sum_(ell!=0)|Khat_y(ell)|^2=6/49,

sum_(ell!=0)|Khat_y(ell)|^4=6/2401.                  (25)
```

This is a literal labelled hole/double dipole, not an uncoloured energy
statement.

## 5. The transition has a compulsory carry

For a seven-unit speed `v`, let `r_v(y)` be the unique generic address
such that

```text
(y+r_v(y))/7 in D_v.
```

Put

```text
T(y)={13y},

kappa_m(y)=floor(my) mod 7.                          (26)
```

Comparing

```text
13v(y+r)/7
```

with the root expression at base `T(y)` gives the exact law

```text
r_(13v)(y)
 =kappa_13(y)-r_v(Ty)                     mod 7.     (27)
```

Applying (27) once and twice to the common core yields

```text
b(y)=kappa_13(y)-a(Ty),

c(y)=a(T^2y)-kappa_169(y)
    =kappa_13(y)-b(Ty)                     mod 7.     (28)
```

The carries are load-bearing. With `h=1` and `y=1/10`,

```text
(a(y),b(y),c(y))=(0,1,4),                            (29)
```

whereas deleting the carries would predict

```text
(a(y),-a(Ty),a(T^2y))=(0,0,6).                      (30)
```

Thus the binary address and the carry-free sign reversal are both
insufficient state spaces. The smallest honest successor state retains
at least

```text
(hole address, double address, kappa_13, owner type). (31)
```

## 6. Scope and next target

The theorem supplies a positive, fixed-phase septimal dipole inside the
counterfactual no-clean branch. It does **not** yet supply:

- invariance of `Y_(b,c)` under `T`;
- a target-`13` colour or terminal endpoint reference;
- an identification of the dipole with THM-2232's signed residual
  eigenline; or
- the `169`-image inclusion required by THM-2257.

The sharp next test is a carry-labelled successor coupling: prove that a
positive portion of (20) retains its owner type after `T`, or construct a
lawful hostile showing that all of it can escape. THM-2362's
successor-role jet is the closest proved mechanism, but its target-shift
observable has not yet been identified with (22).

No thirteen-adic row is removed, the scalar ledger remains `165`, and
LRC(14) remains open.

## 7. Exact companion

The dependency-free exact companion:

- exhausts all `924*7^3=316,932` abstract incidence states before
  filtering and obtains exactly the `385` states in (13);
- verifies that every survivor is a six-address transversal and checks
  the trichotomy;
- reconstructs the full physical address histogram
  `1,24,12,132` over denominator `169`;
- proves (15)--(20) with exact rational arithmetic;
- checks (27)--(28) on every strict endpoint cell and records the
  carry-free hostile (29)--(30);
- verifies every Fourier magnitude and quantitative floor in
  (22)--(25); and
- runs every assertion under ordinary and optimized Python.

Run

```bash
python3 04-computation/lrc14_common_core_transversal_thm2394.py
python3 -O 04-computation/lrc14_common_core_transversal_thm2394.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_common_core_transversal_thm2394.out
```

after LF normalization.
