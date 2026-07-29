---
id: THM-2929
title: "Effective high-gap cubic-null holonomy cutoff"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  For every integer R>=13, the fixed rational box
  [217/100,221/100] x [245/100,249/100] contains exactly one positive
  cubic-null point for the support (0,1,2,R), and its normalized
  quartic endpoint is positive.  The certificate is finite-exact for
  13<=R<=19 and uniform by factorial-ratio contraction for R>=20.
  At R=12 one Bernstein coefficient in this particular certificate is
  negative; this is not a nonexistence or optimal-cutoff assertion.
source: codex-gmc-uniform-width-extension-2026-07-29
depends_on:
  - THM-2914-eventual-high-gap-cubic-null-positive-holonomy-branch
related:
  - THM-2872-four-slot-shared-multipole-quartic-norm-and-response-secant-reduction
  - THM-2879-all-shift-cubic-null-endpoint-holonomy-exit
  - THM-2910-nonconsecutive-cubic-null-endpoint-holonomy-sign-reversal
  - THM-2914-eventual-high-gap-cubic-null-positive-holonomy-branch
script: 04-computation/gmc_effective_high_gap_cubic_null_cutoff_thm2929.py
output: 05-knowledge/results/gmc_effective_high_gap_cubic_null_cutoff_thm2929.out
script_sha256: 458b28af435c21be8d55b572ee1f9b931977878b775b15c6694048eede24a0c2
output_sha256: fc14210f7d192de465e4b7fa50756df256c020caa2fd1ee1934a6003a60cd728
hash_basis: LF-normalized bytes
---

# THM-2929 -- effective high-gap cubic-null holonomy cutoff

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

This is the effective fixed-box continuation promised, but not quantified,
by THM-2914.

## 1. Statement

Let

```text
L(s^k)=k!,                    f_j=s^j/j!,
d_0=f_1-f_0,                  d_1=f_2-f_1,
C_R=f_R-f_2,                  T_k(R)=(kR)!/(R!)^k,
epsilon_R=T_3(R)^(-1/3).                                  (1)
```

For real `(xi,eta)`, put

```text
U_R=d_0+epsilon_R xi C_R,
V_R=d_1+epsilon_R eta C_R,                               (2)
```

and write the binary moments of `U_R+zV_R` as

```text
q_R(z)=L((U_R+zV_R)^2)
      =g_0+2g_1z+g_2z^2,

c_R(z)=L((U_R+zV_R)^3)
      =t_0+3t_1z+3t_2z^2+t_3z^3.                       (3)
```

For every integer `R>=13`, the rational box

```text
B=[217/100,221/100] x [245/100,249/100]                 (4)
```

contains exactly one point `(xi_R,eta_R)` at which `q_R` divides
`c_R`.  Both coordinates are positive.  If `J_R` denotes THM-2872's
quartic endpoint determinant and

```text
H_R=epsilon_R^4 T_4(R)>0,                               (5)
```

then at this point

```text
J_R/H_R > 19                         for 13<=R<=19,
J_R/H_R > 7599/400                   for R>=20.          (6)
```

In particular, this positive cubic-null point is not quartic-null.
The uniqueness assertion is only uniqueness inside the fixed box
`B`; it is not global branch uniqueness.

## 2. Exact divisibility map

The two cubic remainder invariants used in THM-2914 are

```text
F_(R,1)=3t_1g_0g_2-t_3g_0^2-2t_0g_1g_2,
F_(R,2)=3t_2g_0g_2-2t_3g_1g_0-t_0g_2^2.               (7)
```

The polynomials `U_R,V_R` are linearly independent: the `f_0`
coefficient occurs only in `U_R`.  Moreover

```text
L(P^2)=integral_0^infinity P(s)^2 exp(-s) ds>0          (8)
```

for every nonzero real polynomial `P`.  The quadratic Gram matrix is
therefore positive definite.  Consequently

```text
F_R=(F_(R,1),F_(R,2))=0
iff q_R divides c_R.                                    (9)
```

The limiting map from THM-2914 is

```text
F_(infinity,1)=-eta^3+6eta xi^2-4xi^3-14,
F_(infinity,2)=-2(eta^3-3eta^2xi+2xi^3+4).            (10)
```

Set

```text
c=(219/100,247/100),

       [ 18870800/285129003   -6982600/285129003 ]
A =    [ 13965200/285129003     1635200/95043001 ].   (11)
```

Exact multiplication gives

```text
A D F_infinity(c)=I.                                    (12)
```

We work with the equivalent maps

```text
G_R=A F_R,                    G_infinity=A F_infinity. (13)
```

The positive limiting selector is also exact.  If `u=eta^3`, then

```text
p(u)=2u^3+18u^2-729u+54,                               (14)
```

and its unique root in

```text
15069319/10^6 < u < 15069320/10^6                     (15)
```

gives

```text
eta=u^(1/3),
xi=eta(2u^2+21u-603)/189.                              (16)
```

Sturm counts give exactly one root of `p` in each of
`(-25,-24)`, `(0,1)`, and the interval `(15)`.  The root in `(0,1)`
has negative `xi`; `(15)--(16)` give the unique positive-quadrant
limiting zero, and exact rational bounds place it strictly inside
`B`.

## 3. The finite block `13<=R<=19`

For every `12<=R<=19`, the companion encloses the positive algebraic
number `epsilon_R` between adjacent multiples of `10^-18`.  Each
enclosure is verified by cubing and comparing with `1/T_3(R)`; no
floating-point comparison is truth-bearing.  It then reduces all
expressions modulo

```text
epsilon_R^3=1/T_3(R)                                  (17)
```

and converts their restrictions to `B` to exact Bernstein form.

For each `13<=R<=19`, every Bernstein coefficient of the four oriented
faces

```text
-G_(R,1)(217/100,eta),
 G_(R,1)(221/100,eta),
-G_(R,2)(xi,245/100),
 G_(R,2)(xi,249/100)                                  (18)
```

is strictly greater than `1/160`.  Poincare--Miranda therefore gives
a zero of `G_R` in `B`.

Let

```text
Phi_R=id-G_R.                                           (19)
```

On all of `B`, the exact Bernstein row-sum bound is

```text
||D Phi_R||_infinity < 1/5.                            (20)
```

Thus `Phi_R` is a contraction on the convex box.  Any two zeros of
`G_R` would be fixed points of `Phi_R`, so `(20)` makes the
Poincare--Miranda zero unique in `B`.

If

```text
L((U_R+zV_R)^4)
 =A_0+4A_1z+6A_2z^2+4A_3z^3+A_4z^4,
```

the endpoint evaluated by the companion is

```text
J_R=(2A_1g_0-A_0g_1)g_2^2
    -(2A_3g_2-A_4g_1)g_0^2.                            (21)
```

Every two-variable Bernstein coefficient of `J_R/H_R` on `B` is
strictly greater than `19`.  This proves `(6)` for the finite block,
and in particular at its unique zero.

The finite certificate contains 160 face coefficients, including the
`R=12` boundary audit, 560 contraction coefficients, and 252 endpoint
coefficients.  Their exact canonical-expression digest is

```text
3f760cd9dff45ec3064b323623ff0876cbb641f9e637703509fd7d5e5f5af0c7.
                                                                    (22)
```

## 4. The uniform tail `R>=20`

THM-2914 expands each error as finitely many primitive factorial-ratio
atoms.  The companion regenerates this expansion literally rather
than trusting a copied list.  After duplicates are removed, there are

```text
9 quadratic types,       22 cubic types,
45 normalized-quartic types.                           (23)
```

For an atom with `j` copies of the high monomial and offset `a`, put

```text
m_(j,a)(R)
 =product_(h=1)^j (jR+a+h)/(R+1)^j,

r_3(R)=3(3R+1)(3R+2)/(R+1)^2,
r_4(R)=4(4R+1)(4R+2)(4R+3)/(R+1)^3.                  (24)
```

Here `r_k(R)=T_k(R+1)/T_k(R)`.  If `k` is the number of high factors
in an atom, its successive absolute ratio is less than `1/2` once

```text
r_3(R)^k > 8m_(j,a)(R)^3                              (25)
```

for a quadratic or cubic atom, or

```text
r_4(R)^3 > 8r_3(R)^(4-k)m_(j,a)(R)^3                 (26)
```

for a normalized-quartic atom.  For every one of the `76` types, the
numerator of the corresponding difference in `(25)` or `(26)`,
after writing `R=20+n`, has strictly positive integer coefficients.
Thus every primitive error atom contracts by a factor strictly less
than `1/2` from one width to the next for all `R>=20`.  The exact
coefficient digest is

```text
91ec49b07cf15bb869c168be9079a4fbbfce0b13069684d4c56382756b2fdfba.
                                                                    (27)
```

At `R=20`, the exact enclosure

```text
1200603/10^15 < epsilon_20 < 1200604/10^15            (28)
```

and coefficientwise absolute envelopes on `B` give

```text
||G_20-G_infinity||_infinity                  < 1/700,
||D(G_20-G_infinity)||_(infinity,row)          < 1/900,
||J_20/H_20-J_infinity||_infinity             < 1/400. (29)
```

Because each underlying atom subsequently shrinks by more than a
factor two, the same bounds hold for every `R>=20`.  Their exact
envelope digest is

```text
0b2493cd5f6f0593244866311d22eb2c765822fb3ff4acfdc8d765b2877bc561.
                                                                    (30)
```

Exact Bernstein minimization of the limiting expressions gives

```text
oriented face floor
 =13337106413/712822507500 > 1/54,

||D(id-G_infinity)||_(infinity,row)
 =305733792/2376075025 < 13/100,

J_infinity floor
 =1929107681/100000000 > 19.                          (31)
```

Combining `(29)` and `(31)` yields, throughout `B`,

```text
oriented face floor for G_R       > 323/18900 > 0,
||D Phi_R||_(infinity,row)         < 59/450 < 1,
J_R/H_R                            > 7599/400 > 0.      (32)
```

Poincare--Miranda, followed by the same contraction uniqueness
argument as in Section 3, proves the theorem for every `R>=20`.

## 5. The `R=12` boundary

The companion also computes the exact `R=12` data using

```text
6660271211996/10^18
 < epsilon_12
 < 6660271211997/10^18.                               (33)
```

On the lower-`eta` oriented face, the degree-four Bernstein coefficient
with zero-based slot index `3` has upper endpoint strictly below
`-1/200`.  Therefore the coefficientwise positive-face certificate
used above does not extend to `R=12`.

This is only a failure of this fixed box and this Bernstein certificate.
It does **not** prove that the cubic-null branch is absent at `R=12`,
that its endpoint has changed sign, or that `13` is the optimal cutoff
for any broader formulation.

## 6. Meaning and scope

THM-2914 proved eventual local continuation but left its threshold
existential.  The present theorem makes that one lane effective:

```text
support=(0,1,2,R),       fixed box B,       R_0=13.    (34)
```

It proves neither global uniqueness nor an assertion about arbitrary
four-slot supports.  In particular, it does not close the general
four-slot SFC/GMC problem, classify other sign branches, or bypass
the channel-collision obstruction.  The mechanism is the combination
of:

1. exact finite Bernstein certificates for the initial widths;
2. a limiting Poincare--Miranda box with a contraction uniqueness
   estimate; and
3. a coefficientwise-positive factorial-ratio ladder that propagates
   one exact tail envelope to every later width.

## 7. Exact companion

The companion hash-pins and imports the independently audited THM-2914
script, regenerates the moment tensors from `(1)--(3)`, and verifies:

1. the limiting selector and all three real eliminant roots;
2. every exact algebraic-number face, derivative, and endpoint
   certificate for `12<=R<=19`;
3. all `76` primitive tail-ratio inequalities;
4. the exact `R=20` value, derivative, and endpoint envelopes; and
5. all displayed rational constants and certificate digests.

Reproduction:

```text
python3 04-computation/gmc_effective_high_gap_cubic_null_cutoff_thm2929.py
python3 -O 04-computation/gmc_effective_high_gap_cubic_null_cutoff_thm2929.py
```

Both executions agree byte-for-byte with

```text
05-knowledge/results/gmc_effective_high_gap_cubic_null_cutoff_thm2929.out.
```

The transcript is a compact evidence digest; the mathematical
certificate consists of the exact constructions and explicit checks
in the companion together with the argument above.
