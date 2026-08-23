---
id: THM-3829
title: "Misaligned second-row R2Z2 10/7 profiles do not enter the cubic pseudo-plane Darboux packet"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the first
  canonical second-row extension of THM-3821, suppose the nonzero r^2z^2
  profile X in the first coordinate has its necessarily proportional mate
  Y=lambda X, but the preceding rz^2 profiles are misaligned:
  L=T-lambda S!=0.  Then the pair is not Darboux.  The r6 bucket gives a
  10/7 tower, the r5 bucket integrates S, and a cancellation-safe r4
  valuation forces M=e v^2 R with R a unit on the zero set of v.  The r3
  bucket then has an unavoidable first local term at every nonzero v-root
  and an odd 2d+1 origin term otherwise.  Combined with the independently
  audited aligned closure THM-3828, this excludes every X!=0 pair in the
  fixed r^2z^2 second-row ansatz.  The one-sided orientation X=0,Y!=0 and
  higher canonical slots remain OPEN.
source: jc_zero_debt_lift / cubic-pseudoplane second-row profile lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn-boundary, 2026-08-23).
  The audit checked completeness of the S integral; every nonzero-root and
  origin comparison in the cancellation-safe r4 source, including zero
  orders and the tied payment; the global M=e v^2R factorization; the full
  non-p and non-g order floors in every zero/nonzero branch; the p/g
  resonance bounds; polynomiality of all terminal remainders; and both
  final local coefficients.  Normal and optimized runs byte-match the
  frozen transcript and the raw hashes agree.  The deterministic companion
  has 58
  active gates checking the Poisson Casimir, monic canonical reduction,
  arm and top Wronskian buckets, exact 10/7 tower and r5 integration, the
  full pretyped r4 source polynomial and both local valuation families,
  the r4 root payment, the zero/nonzero N and H branches, the p and g
  lower-bound differential blocks at nonzero roots and the origin, all
  four N/H terminal specializations, polynomiality of every terminal
  remainder, and both final leading coefficients.  Normal and optimized
  runs byte-match the frozen transcript.  No finite-field inference is
  used.
depends_on:
  - THM-3828-proportional-second-row-r2z2-profile-nonentry
related:
  - THM-3821-cubic-pseudoplane-rz2-odd-ladder-terminal-riccati-gate
  - THM-3814-nodal-rz-kummer-profile-degree-gate
  - THM-3811-nodal-arm-bezout-law-for-cubic-pseudoplane-darboux-pairs
script: 04-computation/jc2_cubic_pseudoplane_misaligned_second_row_thm3829.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_misaligned_second_row_thm3829.out
script_sha256: 27b15411a51ad8e6facf1dfa39c724c278a5a9811d7898c88a529935881fed67
output_sha256: d728949c7f9f331a80f64b446e4b6c3f84bb5f5899cac1c68a46a4edb88549be
semantic_sha256: d045d251b2f0016d9909c329a37a8d30c8e4b7c50d0fef84e9000f53012fe8c6
hash_basis: raw LF bytes
---

# THM-3829 -- the misaligned 10/7 row is also empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be an
algebraically closed field of characteristic zero and put

```text
B=k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary profiles

```text
f,g,h,kappa,p,q,S,T,X,Y in k[e]
```

consider the first canonical second-row extension of THM-3821:

```text
A=e^2-z/3+r g+z^2 f+rz p+rz^2 S+r^2z^2 X,
C=e^3-e-ez/2+r h+z^2 kappa+rz q+rz^2 T+r^2z^2 Y.      (2)
```

Suppose `X!=0`.  If

```text
Y=lambda X,                 L:=T-lambda S !=0          (3)
```

for some `lambda in k`, then

```text
{A,C}!=1.                                                (4)
```

Consequently, THM-3828 and this theorem together exclude every Darboux pair
with `X!=0` in (2): the top bucket first forces `Y=lambda X`; THM-3828 closes
`L=0`, and (4) closes `L!=0`.

The statement is only about the displayed fixed ansatz.  It does **not**
exclude the one-sided top orientation `X=0,Y!=0`, any further canonical slot,
or a planar Jacobian counterexample by another profile.

## 1. Exact canonical buckets

Reduce `{A,C}-1` by the monic relation

```text
z^3=r^2e+r.
```

The result has `z`-degree at most two, with no quotient loss.  The pure `z`
and top `r^7` buckets are

```text
[z]   = (36e^2f-24e kappa-12f+1)/2,
[r^7] = 30e^2(XY'-YX').                                (5)
```

Thus a vanishing top bucket with `X!=0` makes `(Y/X)'=0` in `k(e)`, hence
`Y=lambda X`.  Introduce the remaining target-direction differences

```text
M=kappa-lambda f,   N=q-lambda p,   H=h-lambda g,
L=T-lambda S.                                            (6)
```

The arm becomes

```text
D f=2eM-1/12,             D=3e^2-2lambda e-1.            (7)
```

In particular `M!=0`: otherwise the nonconstant quadratic `D` would divide
the unit `-1/12` in `k[e]`.

The next three useful buckets are

```text
[r^6]   =-3e(7eLX'-10eXL'-2LX),
[r^4z^2]=15e(2XN'-NX'),
[r^4z]  =3(10eXH'-3eHX'-2XH).                            (8)
```

Every equation below is an identity in the polynomial domain `k[e]`; the
displayed harmless powers of `e`, `v`, and nonzero scalar constants may
therefore be cancelled as polynomials.  We never evaluate and divide at a
root.

## 2. The 10/7 tower and the exact `S` integral

Because `X,L` are nonzero, the first equation in (8) gives

```text
7 ord_pi(X)-10 ord_pi(L) = 2                 if pi=e,
7 ord_pi(X)-10 ord_pi(L) = 0                 if pi!=e.   (9)
```

The nonnegative solutions are

```text
X=alpha e^6v^10,             L=beta e^4v^7,              (10)
```

where `alpha,beta in k*` and `v in k[e]\{0}`.  Conversely (10) solves the
`r^6` equation exactly.  This is the `10/7` tower left OPEN by THM-3828.

After (10), the full `r^5` equation integrates without dividing by `M`:

```text
S=(10alpha/(7beta))e^2v^3 M
  +(2alpha/7)e^5v^10+gamma e^4v^7,          gamma in k. (11)
```

The exact companion substitutes (11) into the unreduced bucket and obtains
zero identically.  If `S_1,S_2` are two solutions, their difference `W`
satisfies `-LW'+WL'=0`, so `(W/L)'=0` in `k(e)` and `W` is a scalar multiple
of `L`; its scalar is absorbed into `gamma`.  Thus no integration constant or
exceptional `M=0` branch has been discarded.

## 3. Cancellation-safe `r^4` typing

Substitute (10)--(11) into `[r^4]` and remove its nonzero common polynomial
factor `-3e^3v^2/(7beta)`.  Call the remaining source `P_4`.  The complete
source is frozen in the companion.  Its only two blocks capable of being
first at a zero of `v` are

```text
30alpha M(4eM v'-evM'+2Mv),
49beta^2 e^2v^4(-4efv'+evf'-2fv).                       (12)
```

All other terms have one of the strictly later forms

```text
e^6v^15,  e^7v^14v',
e^3Mv^8,  e^4Mv^7v',  e^4v^8M',
e^2Mv^5,  e^3Mv^4v',  e^3v^5M'.                         (13)
```

This separation prevents a hidden cancellation assumption.

Let `rho!=0` be a root of `v`, write

```text
v=(e-rho)^m u,       M=(e-rho)^s w,
m>=1,                u(rho)w(rho)!=0.                   (14)
```

The first block in (12) has order `2s+m-1` and leading coefficient

```text
30alpha rho u(rho)w(rho)^2(4m-s).                       (15)
```

If `s>0`, the arm (7) makes `D(rho)f(rho)=-1/12`; hence both `D(rho)` and
`f(rho)` are nonzero.  The second block then has order `5m-1` and leading
coefficient

```text
-196beta^2 rho^3 f(rho)m u(rho)^5.                      (16)
```

If `s<2m`, (15) is uniquely first; its only resonance `s=4m` lies outside
that range.  If `s>2m`, then `s>0` and (16) is uniquely first.  Therefore

```text
s=2m.                                                    (17)
```

The origin must be handled separately because the explicit `e` powers in
(12) change the orders.  Put

```text
v=e^d u,       M=e^s w,       d,s>=0,       u(0)w(0)!=0. (18)
```

Equation (7) fixes `f(0)=1/12`.  The two orders and leading coefficients are

```text
M^2 block:  2s+d,     30alpha u(0)w(0)^2(4d-s+2),
f block:    5d+2,    -98beta^2 f(0)u(0)^5(2d+1).        (19)
```

If `s<2d+1`, the first block is unique and its resonance `s=4d+2` lies above
the range.  If `s>2d+1`, the second is unique.  Hence

```text
s=2d+1.                                                  (20)
```

Equations (17) and (20), at every irreducible factor of `v`, give the global
factorization

```text
M=e v^2 R,              R in k[e],       gcd(R,v)=1.     (21)
```

At the tied order, the leading coefficient of `P_4=0` also gives, uniformly
at every zero `xi` of `v` (including the origin),

```text
15alpha R(xi)^2=49beta^2 f(xi).                          (22)
```

We need only its consequence: `R` and `f` are units in every local ring at a
zero of `v`.  This also follows for `f` directly from (7): after (21),

```text
D f=2e^2v^2R-1/12.                                      (23)
```

## 4. The `N` and `H` branches

The last two equations in (8) are honest zero/nonzero dichotomies:

```text
N=0       or       N=delta e^3v^5,    delta in k*,
H=0       or       H=theta e^2v^3,    theta in k*.       (24)
```

They follow by taking prime valuations exactly as in (9).  No division by
`N` or `H` is used in the zero branches.

The next bucket is needed only to keep lower profiles from arriving earlier
than the terminal obstruction.  In `[r^3z^2]`, in either branch of `N`, the
entire differential block containing `p` is

```text
735beta^2 e^3p v^4v'
-147beta^2 e^3v^5p'
+441beta^2 e^2p v^5.                                   (25)
```

At a nonzero root (14), if `t=ord_rho(p)<3m`, (25) is uniquely first among
the terms of the source: it has order `t+5m-1` and coefficient

```text
147beta^2 rho^3u(rho)^5 p_t(5m-t).                      (26)
```

Its resonance `t=5m` lies above the assumed range; every non-`p` term has
order at least `8m-1`.  At the origin, if `t=ord_0(p)<3d+2`, the corresponding
order and coefficient are

```text
t+5d+2,       147beta^2u(0)^5 p_t(5d-t+3),              (27)
```

while every non-`p` term has order at least `8d+4`.  Again the resonance is
outside the range.  Therefore

```text
p=e^2v^3 U,                 U in k[e].                  (28)
```

If `N!=0`, the `[r^3z]` bucket similarly has the entire `g` block

```text
147beta^2 e^4g v^6v'
-49beta^2 e^4v^7g'
+98beta^2 e^3g v^7.                                    (29)
```

Using (21), (24), and (28), its other terms have order at least `8m-1` at a
nonzero `v`-root and `8d+4` at the origin.  If
`t=ord_rho(g)<m`, the leading coefficient of (29) is

```text
49beta^2 rho^4u(rho)^7 g_t(3m-t),                       (30)
```

whose resonance lies above the range.  At the origin the corresponding
coefficient is

```text
49beta^2u(0)^7 g_t(3d-t+2)                              (31)
```

under `t<d+1`.  Consequently

```text
g=e v V                  when N!=0.                     (32)
```

When `N=0`, every `g` term disappears from the decisive bucket below, so no
claim like (32) is needed or made.

## 5. The terminal root obstruction

Substitute (10), (11), (21), (24), (28), and, only when `N!=0`, (32) into
the full `[r^3]` bucket.  In each of the four branches

```text
(N=0 or N!=0) x (H=0 or H!=0),                          (33)
```

the bucket factors exactly as

```text
[r^3]=-3e^2v Q/(7beta),                                  (34)

Q= -56beta eRf v'
   +28beta eRv f'
   -28beta efv R'
   -28beta Rfv
   +e^2v^4 Psi,                                         (35)
```

where `Psi` is a polynomial in the displayed profiles and their derivatives.
The exact companion independently constructs all four `Psi` and checks that
their denominators are one.  Thus (35) is not a Laurent truncation.

Suppose first that `v` has a nonzero root `rho` of multiplicity `m`.  By
(21) and (23), `R(rho)f(rho)!=0`.  In (35), the first term is uniquely first,
of order `m-1`, with coefficient

```text
-56beta rho R(rho)f(rho)m u(rho),                       (36)
```

which is nonzero.  This contradicts `Q=0`.

Therefore `v` has no nonzero root.  Since `k` is algebraically closed,

```text
v=c e^d,                  c in k*,       d>=0.           (37)
```

At the origin, the first and fourth terms in (35) are the only terms of
order `d`; together they have coefficient

```text
-28beta R(0)f(0)c(2d+1).                                (38)
```

Here `R(0)!=0` by (21), `f(0)=1/12` by (19), and `2d+1` cannot vanish in
characteristic zero.  Thus (38) is nonzero, another contradiction.

This proves (4).

## 6. Exact boundary and design consequence

The proof closes the entire previously OPEN `L!=0` tower of THM-3828,
including all four `N/H` zero seams.  Together, the two theorems show:

```text
fixed ansatz (2) + X!=0  ==>  no Darboux pair.           (39)
```

The first genuinely live orientation in this profile grammar is now the
one-sided branch

```text
X=0,                 Y!=0,                               (40)
```

or else a higher canonical slot whose top bucket changes (5).  Neither is
covered here.  In particular (39) is an ansatz obstruction, not a theorem
against planar Jacobian counterexamples in general.

## 7. Reproducibility

Run

```bash
python3 04-computation/jc2_cubic_pseudoplane_misaligned_second_row_thm3829.py
python3 -O 04-computation/jc2_cubic_pseudoplane_misaligned_second_row_thm3829.py
```

Both executions must byte-match

```text
05-knowledge/results/jc2_cubic_pseudoplane_misaligned_second_row_thm3829.out
```

The companion works over exact SymPy expressions.  There is no random or
finite-field step, no disabled Python `assert`, and no inference from sampled
orders.  The local-order formulas are symbolic consequences of the frozen
full bucket identities.
