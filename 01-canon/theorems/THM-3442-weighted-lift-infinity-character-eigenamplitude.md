---
id: THM-3442
title: "The weighted-lift infinity torsor carries a nonzero fundamental character amplitude"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / INDEPENDENT AUDIT REQUIRED.
  On the cyclic Q=infinity branch torsor of THM-3440, the normalized source
  coordinate x has leading branch profile j -> const*zeta_n^j.  Its fundamental
  Fourier coefficient is therefore a formal unit.  At n=91, CRT identifies
  this with a nonzero (1,1) line after the explicit factor-root normalization.
  This supplies a canonical Keller-side eigen-amplitude after branch and root
  markings, but not an LRC amplitude intertwiner, bispectrum identity, physical
  current, or LRC(14) consequence.
source: codex2 Puiseux reconstruction-current derivation, 2026-08-15
depends_on:
  - THM-3438-weighted-lift-keller-degree-spectrum
  - THM-3440-weighted-lift-cyclic-infinity-torsor-and-7x13-character-grid
related:
  - THM-2512
  - THM-3431-valuative-persistence-multiset-and-lrc-jc-boundary
script: 04-computation/jc_weighted_lift_infinity_character_thm3442.py
output: 05-knowledge/results/jc_weighted_lift_infinity_character_thm3442.out
---

# THM-3442 -- the weighted-lift infinity torsor carries a nonzero fundamental character amplitude

**RESERVED / PROVISIONAL PROOF CANDIDATE / INDEPENDENT AUDIT REQUIRED.**

## 1. Exact local expansion

Use a lawful weighted lift from THM-3438.  Let

```text
d=deg(p),             n=d+1,
p(w)=a_d w^d+lower,   a_d!=0,
R(w)=integral_0^w p=(a_d/n)w^n+lower.                 (1)
```

Fix a generic `P=P_0` and the target line `(A,B,C)=(Q,P_0,1)`.  Put

```text
t=Q^(-1),             s^n=t,
w=s^(-1)W.                                               (2)
```

The inverse equation `R(w)-P_0w+cQ=0`, after multiplication by `t`, becomes

```text
(a_d/n)W^n+c+terms divisible by s=0.                    (3)
```

Choose `rho` with

```text
(a_d/n)rho^n+c=0.                                       (4)
```

Over a field containing `rho` and `zeta_n`, Hensel's lemma gives `n` branches

```text
W_j(s)=rho zeta_n^j+O(s),
w_j(s)=s^(-1)W_j(s),             j in Z/n.              (5)
```

These are the marked cyclic branches of THM-3440.

The weighted reconstruction has

```text
gamma=(P_0-p(w))/c,          x=1/gamma                  (6)
```

on `C=1`.  Substitution of `(1)` and `(5)` gives

```text
gamma_j(s)=-(a_d/c)rho^(n-1)zeta_n^(j(n-1))s^(-(n-1))
           +O(s^(-(n-2))),                              (7)

X_j(s):=s^(-(n-1))x_j(s)
       =-(c/a_d)rho^(-(n-1))zeta_n^j+O(s).              (8)
```

The exponent in `(8)` is `+j` because `-(n-1)=1 mod n`.

## 2. Exact Fourier nonvanishing

With the marked branch ordering, define

```text
Xhat(k;s)=(1/n) sum_(j mod n) X_j(s) zeta_n^(-kj).      (9)
```

Orthogonality and `(8)` give

```text
Xhat(1;0)=-(c/a_d)rho^(-(n-1))!=0,
Xhat(k;0)=0                         for k!=1.            (10)
```

Therefore `Xhat(1;s)` is a unit of the formal power-series ring.  In
particular it is identically nonzero as an algebraic Puiseux germ.  This is a
structural consequence of the weighted reconstruction, not a numerical
sample and not a genericity guess.

The source coordinate `x` is thus a natural Keller-side local current: after
the unique ramification scaling `s^(-(n-1))`, its leading branch profile spans
the fundamental `C_n` character line.

## 3. The degree-91 `(1,1)` line

For THM-3438's degree-91 seed, `a_d=-1`, `c=1`, and `(4)` is

```text
rho^91=91.                                               (11)
```

Let `tau` be the oriented `C_91` inertia generator.  In the CRT coordinates
of THM-3440, use the factor generators

```text
(1,0)=tau^78,             (0,1)=tau^14.                 (12)
```

Normalize the corresponding primitive roots by

```text
eta_7=zeta_91^78,         eta_13=zeta_91^14.            (13)
```

Since every branch label has the CRT reconstruction

```text
j=78(j mod 7)+14(j mod 13) mod 91,                      (14)
```

one has

```text
zeta_91^j=eta_7^(j mod 7) eta_13^(j mod 13).            (15)
```

Thus the fundamental character in `(8)--(10)` is exactly the `(1,1)` tensor
character relative to the explicit factor-root normalization `(12)--(13)`,
and

```text
Xhat_91((1,1);0)=rho^(-90)!=0.                          (16)
```

This upgrades THM-3440's common character **carrier** to one explicit nonzero
Keller-side eigen-amplitude.

## 4. Typed bridge and missing sidecar

After choosing the inertia orientation, `rho`, and the factor roots `(13)`,
the JC line in `(16)` is one-dimensional.  THM-2512's LRC `(1,1)` coefficient
also lives in a one-dimensional `7 tensor 13` eigenspace.  A linear map sending
the marked LRC basis to the normalized germ `(8)` is therefore explicit.

What is not canonical is the amplitude comparison: changing `rho` or the
inertia generator rescales/permutates the JC basis, while the LRC coefficient
has its own ancestry, centering, and response normalization.  A theorem-level
D5 map still needs an intertwiner proving that those normalizations are
compatible.  Dimension matching and separate nonvanishing do not provide it.

| field | exact content |
|---|---|
| source | weighted inverse branches and reconstructed coordinate `x` |
| target | a nonzero fundamental Fourier/Puiseux eigenline |
| map | ramification normalization `(2),(8)` followed by branch Fourier transform `(9)` |
| preserved | cyclic character, CRT factor address, valuation, and formal nonvanishing |
| destroyed | finite-target values, LRC ancestry/centering, positivity, physical time, and response normalization |
| needed sidecar | a gauge-compatible amplitude intertwiner between the LRC response and `X` |
| cheapest hostile | reverse `tau`: mode `1` becomes mode `-1`, so unmarked `(1,1)` equality fails |

## 5. Boundaries

- The result is a formal germ over target infinity, not a finite Jelonek
  component or a claim about every target value.
- The nonzero object is one Fourier coefficient of a reconstructed source
  coordinate.  It is not a bispectrum and does not prove THM-2512's LRC
  contraction nonzero.
- The `(1,1)` label uses the explicit CRT root normalization `(13)`.  With
  conventional roots `zeta_7=zeta_91^13`, `zeta_13=zeta_91^7`, the same
  character receives a unit-rescaled dual label; this is gauge, not a
  contradiction.
- No additive map to THM-3437's characteristic-zero Tor tower is asserted.
- Nothing here proves LRC(14), `JC(2)`, or a physical D5 current.

## 6. Exact companion

The companion must verify the Newton/Hensel leading equations, exponent
conversion in `(7)--(10)`, the degree-91 seed coefficient, CRT identities
`(12)--(15)`, Fourier orthogonality, and generator/root hostile controls.
Until that replay and an independent audit are merged, this file remains
provisional.
