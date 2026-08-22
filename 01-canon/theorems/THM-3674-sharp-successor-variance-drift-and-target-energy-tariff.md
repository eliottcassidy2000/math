---
id: THM-3674
title: "Sharp successor variance, drift, and target-energy tariff"
status: >
  PROVED; PENDING INDEPENDENT HOSTILE AUDIT.  The THM-2365 drift splits
  exactly into p^-2 times the variance of its successor marginal plus the
  nonzero-deep, nonzero-target energy.  For arbitrary tensors the sharp
  variance constant is p^-2.  The lawful diagonal-plane zero improves it
  sharply to 1/(p(p-1)), even under nonnegativity.  At p=13, one distinct
  three-site successor defect Delta forces drift at least |Delta|^2/158184
  and deep-target energy at least |Delta|^2/2056392.  This is a quantitative
  transfer theorem, not a covering-row nonvanishing result or LRC(14) proof.
source: kps-s194 / Averroes marginal decomposition, 2026-08-21
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
related:
  - THM-3665-lrc-support-minimal-three-twist-target-detector
  - THM-3670-lrc-successor-mass-transfer-and-thirteen-number-gate
  - THM-3672-lrc-successor-mass-all-packet-positive-control
---

# THM-3674 -- successor variation pays a sharp target-energy tariff

**PROVED; PENDING INDEPENDENT HOSTILE AUDIT.**  This theorem keeps all
normalizations from THM-2365.  Its first identity is an exact orthogonal
decomposition, not an estimate.

## 1. Drift, successor marginal, and deep-target energy

Let `p` be prime and let

```text
H:F_p^3 -> C,

B(a,b,h)
 =p^-3 sum_(r,s,t)H(r,s,t)zeta^(ar+bs+ht),          (1)
```

where `zeta=exp(2*pi*i/p)`.  Let `P` be THM-2365's orthogonal projection
onto the functions of `r-t`, and put

```text
D=p^-3 sum_(r,s,t)|H(r,s,t)-(PH)(r,s,t)|^2.         (2)
```

Define the successor marginal and its normalized variance by

```text
S(s,t)=sum_r H(r,s,t),
S_bar=p^-2 sum_(s,t)S(s,t),
Var(S)=p^-2 sum_(s,t)|S(s,t)-S_bar|^2.              (3)
```

THM-2365's target coordinate of `B(a,b,h)` is `(b,a+h)`.  Reparameterize
`q=a+h` and define the nonzero-deep, nonzero-target energy

```text
E_dt=sum_((b,q)!=(0,0), a!=0)|B(a,b,q-a)|^2.        (4)
```

Then the exact identity is

```text
D=Var(S)/p^2+E_dt.                                  (5)
```

Indeed, with the same plus-sign convention as (1),

```text
S_hat(b,q)
 =p^-2 sum_(s,t)S(s,t)zeta^(bs+qt)
 =p B(0,b,q).                                       (6)
```

Parseval gives

```text
Var(S)=p^2 sum_((b,q)!=(0,0))|B(0,b,q)|^2.          (7)
```

Meanwhile THM-2365 (17), after `h=q-a`, is

```text
D=sum_((b,q)!=(0,0), a in F_p)|B(a,b,q-a)|^2.      (8)
```

Separating the `a=0` slice and using (7) proves (5).  In particular

```text
D>=Var(S)/p^2.                                      (9)
```

The constant `p^-2` is sharp for arbitrary `H`.

## 2. The lawful diagonal zero gives the sharp improvement

Now impose the exact THM-2365 condition

```text
H(t,s,t)=0 for every s,t.                           (10)
```

Its target-null identity is

```text
sum_(a in F_p)B(a,b,q-a)=0                          (11)
```

for every `(b,q)`.  On a nonzero target line write
`z_a=B(a,b,q-a)`.  Then

```text
|z_0|^2
 =|sum_(a!=0)z_a|^2
 <=(p-1)sum_(a!=0)|z_a|^2.                         (12)
```

Summing (12), and using (5) and (7), yields the sharp pair

```text
D>=Var(S)/(p(p-1)),

E_dt>=Var(S)/(p^2(p-1)).                           (13)
```

It also recovers THM-2365 (19): `E_dt>=D/p`.

For completeness, the first inequality in (13) has a direct time-domain
proof.  Put `Q=(I-P)H` and `u=r-t`.  Condition (10) gives `Q(t,s,t)=0`.
For each `(s,t)`, Cauchy--Schwarz on the remaining `p-1` values of `u`
gives

```text
|S(s,t)-S_bar|^2
 <=(p-1)sum_r|Q(r,s,t)|^2.                          (14)
```

Summing and applying the normalizations in (2)--(3) gives (13).

## 3. Equality and nonnegative sharpness

For arbitrary `H`, equality in (9) holds exactly when

```text
Q(r,s,t)=g(s,t),              sum_(s,t)g(s,t)=0.    (15)
```

Under (10), equality throughout (13) holds exactly when

```text
Q(t+u,s,t)=0       for u=0,
Q(t+u,s,t)=g(s,t)  for u!=0,
sum_(s,t)g(s,t)=0.                                  (16)
```

These extremizers can be nonnegative.  Choose a real mean-zero `g` and
`M>=max|g|`, then set

```text
H(t+u,s,t)=0          for u=0,
H(t+u,s,t)=M+g(s,t)   for u!=0.                     (17)
```

Its projection is zero on `u=0` and `M` on `u!=0`; hence (16) holds and
every bound in (13) is attained.  Positivity alone cannot improve the
constants.

## 4. One three-site defect gives an explicit tariff

Let `x,a,b in F_p^2` and put

```text
Delta=S(x)+S(x+a)-2S(x+b),

w=delta_x+delta_(x+a)-2delta_(x+b),
kappa=sum_y|w(y)|^2.                                (18)
```

Since `sum_y w(y)=0`, Cauchy--Schwarz gives the sharp estimate

```text
Var(S)>=|Delta|^2/(kappa p^2).                      (19)
```

Equality holds exactly when `S-S_bar` is a scalar multiple of `w`.  If the
three sites are distinct, `kappa=6`; composing (19) with (13) gives

```text
D>=|Delta|^2/(6p^3(p-1)),

E_dt>=|Delta|^2/(6p^4(p-1)).                       (20)
```

At `p=13`, this is

```text
Var(S)>=|Delta|^2/1014,

D>=|Delta|^2/158184,

E_dt>=|Delta|^2/2056392.                            (21)
```

The distinct-site hypothesis is load-bearing.  Directly combining repeated
sites in (18) gives

```text
kappa=6  if a,b!=0 and a!=b,
kappa=8  if a=0 and b!=0,
kappa=2  if b=0 and a!=0,
kappa=2  if a=b!=0,
kappa=0  if a=b=0.                                  (22)
```

Taking `g` proportional to the combined mask `w` in (16)--(17) attains
(19)--(20), so all displayed constants are sharp simultaneously.

For THM-3670, `a=(1,0)` and `b=(0,1)`, so the three sites are distinct.
Any exact nonzero rational successor defect therefore pays the tariff (21).
The theorem does not prove that such a defect exists on a hypothetical
covering row, retain a preselected frequency triangle, or prove LRC(14).

**QED.**
