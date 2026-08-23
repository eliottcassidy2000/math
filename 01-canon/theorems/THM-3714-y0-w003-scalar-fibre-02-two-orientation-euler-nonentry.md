---
id: THM-3714
title: "All-scale W003 scalar-fibre-02 two-orientation Euler nonentry"
status: >
  PROVED + VERIFIED-EXACT.  On the W003 ray, every post-parity placement whose
  scalar fibre is 02+11 is Darboux-empty in the y=0 collision ring: one
  orientation exists at n=2 and both exist for n>=3.  The upper rows force the
  hidden weight-minus-two coefficient into a common Euler charge q+2=delta p,
  so the scalar row factors H'K+delta HK'.  The n=2 boundary instead forces
  K|H and factors J'K+2JK'.  This closes the complete named scalar-fibre pair,
  not all of W003, the full 3x4 cell, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  SELF-AUDITED.  The exact companion reconstructs the post-parity
  actual-support census on the named fibre, checks both generic integrations
  in all three gcd residue branches, verifies the charge and nonlinear-wedge
  identities, and checks the exceptional n=2 peel.  Normal and optimized runs
  byte-match the frozen transcript.  Independent hostile audit remains open.
depends_on:
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
  - THM-3712-y0-w007-complete-all-scale-arm-transport-nonentry
script: 04-computation/jacobian_y0_w003_scalar02_euler_thm3714.py
output: 05-knowledge/results/jacobian_y0_w003_scalar02_euler_thm3714.out
script_sha256: b91a594f3734aac39b9324dbde72be8a24bdfc78c0d5772996185411aea9f4fa
output_sha256: 8f36643788358195faf7070b3d50f89e45d458ad3aa342a3965d9f52cc3349d9
hash_basis: LF-normalized bytes
---

# THM-3714 -- the W003 scalar fibre 02 is empty

**PROVED + VERIFIED-EXACT.**  Work over `C` in the THM-3696 collision ring.
Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

## 1. Exact placement scope

The W003 ray has gaps and fibres

```text
(X,Y;U,V,W)=(n,2n;n,n,n),                    n>=1,

00; 01+10; 02+11; 03+12+20; 13+21; 22; 23.            (2)
```

Restrict to scalar fibre `02+11`, deduplicate actual supports, and apply the
inherited THM-3613 parity gate.  There is no placement at `n=1`, one at
`n=2`, and exactly the following two orientations for `n>=3`:

| family | `wt(P)` | `wt(Q)` | live scales |
|---|---|---|---|
| A | `(-n-2,-2,2n-2)` | `(1-n,1,n+1,2n+1)` | `n>=2` |
| B | `(1-n,1,2n+1)` | `(-n-2,-2,n-2,2n-2)` | `n>=3` |

The companion reconstructs this census directly.  It counts an actual
support pair, not a choice of arm anchor.

## 2. Family A for n>=3

Put

```text
delta=gcd(n-1,3),       alpha=(n+2)/delta,
beta=(n-1)/delta.                                       (3)
```

The negative singleton `00` and positive singletons `22,23` give

```text
f_0=aH^alpha,                   g_0=b_0H^beta,
f_2=d_0K^(2n-2),                g_2=e_0K^(n+1),
g_3=t_0K^(2n+1).                                      (4)
```

The positive gcd is one because `gcd(n+1,2n+1)=1`.  Write the hidden arm
coefficient `f_1=M` and its mate `g_1=L`.  The upper double `13+21=0`
integrates completely to

```text
L=lambda K+rho MK^3,
rho=(2n+1)t_0/[2(n-1)d_0].                            (5)
```

The triple `03+12+20=0` then gives

```text
M=kappa K^-2
  -(2n+1)a t_0/[(n+1)e_0] H^alpha K^n
  +2(n-1)b_0d_0/[(n+1)e_0] H^beta K^(n-3).           (6)
```

Equation `(6)` is an identity in `C(b)`.  For `n>=3`, its last two terms are
polynomials.  Since `K` is nonconstant, polynomiality of `M` forces
`kappa=0`.

Every remaining monomial `H^pK^q` of `M` obeys the exact charge law

```text
q+2=delta p.                                           (7)
```

Consequently, with

```text
E_delta(H,K)=H'K+delta HK',                            (8)
```

one has term by term

```text
W_(-2,1)(H^pK^q,K)=pH^(p-1)K^q E_delta.               (9)
```

The nonlinear part of `(5)` introduces no new charge:

```text
W_(-2,1)(M,MK^3)=3MK^2 W_(-2,1)(M,K).                (10)
```

The other scalar address also factors `E_delta`, because
`delta alpha=n+2`.  Thus the scalar row `02+11` is

```text
W_(-n-2,n+1)(f_0,g_2)+W_(-2,1)(M,L)
       =E_delta Phi_A                                  (11)
```

for the explicit polynomial `Phi_A` checked by the companion.  It cannot be
one because Section 5 proves `E_delta` is a nonunit.

## 3. Family B for n>=3

Keep `(3)` and now write

```text
f_0=aH^beta,                    g_0=b_0H^alpha,
f_2=d_0K^(2n+1),                g_2=e_0K^(n-2),
g_3=t_0K^(2n-2).                                      (12)
```

The three positive weights again have gcd one: any common divisor divides
both `3=(2n+1)-(2n-2)` and `2=(2n-2)-2(n-2)`.  Let `M=g_1` be the hidden
weight-minus-two arm coefficient and `L=f_1` its weight-one mate.  The upper
double and triple integrate to

```text
L=lambda K+rho MK^3,
rho=(2n+1)d_0/[2(n-1)t_0],                            (13)

M= 2(n-1)b_0t_0/[(n-2)e_0] H^alpha K^n
  -4(n-1)^2a t_0^2/[(2n+1)(n-2)d_0e_0]
       H^beta K^(n-3).                                (14)
```

As in A, the omitted homogeneous term before polynomiality is
`kappa K^-2`, hence zero.  Both monomials in `(14)` obey `(7)`.  Now

```text
W_(1,-2)(L,M)=-W_(-2,1)(M,L),                         (15)
```

so the scalar row again has an exact factorization

```text
W_(1-n,n-2)(f_0,g_2)+W_(1,-2)(L,M)
       =E_delta Phi_B.                                (16)
```

This closes B at every live scale.

## 4. The exceptional family A at n=2

At `n=2`, `(6)` has one negative extra exponent.  The complete polynomiality
identity is

```text
K^2M=kappa-5a t_0/(3e_0) H^4K^4
             +2b_0d_0/(3e_0) HK.                    (17)
```

Evaluate at a root of the nonconstant polynomial `K`: `kappa=0`.  Divide by
one `K`; polynomiality then forces `K|H`.  Put `H=KJ`.  The arm coefficient
and mate become

```text
M=-5a t_0/(3e_0) K^6J^4+2b_0d_0/(3e_0)J,
L=lambda K+5t_0/(2d_0)MK^3.                          (18)
```

Direct substitution in the upper double, triple, and scalar row gives

```text
scalar row=(J'K+2JK') Phi_2                            (19)
```

with `Phi_2` polynomial.  The exact companion freezes the full identity and
all signs.

## 5. The Euler factors are nonunits

Negative coefficient membership in `(4)` or `(12)` and squarefreeness give
`h|H`.  Positive coefficient membership and the gcd-one power form give
`b|K`.  Hence

```text
deg E_delta=deg H+deg K-1>=2,                          (20)
```

with nonzero leading multiplier `deg H+delta deg K`.

In the exceptional branch `H=KJ`, `h|H` gives `deg H>=2` and `b|K` makes
`K` nonconstant.  Therefore

```text
deg(J'K+2JK')=deg H-1>=1,                             (21)
```

with nonzero leading multiplier `deg J+2deg K`.  No displayed Euler factor
is a unit.

## 6. Conclusion and scope

Equations `(11),(16),(19)` close every retained W003 placement whose scalar
fibre is `02+11`, in both orientations and at every scale.  Output exchange
gives the transposed conclusion.  This does not close the other W003 scalar
fibres, W004--W006, arbitrary `3 x 4` supports, general quartic C3 data, or
`JC(2)`.

Run

```bash
python3 -B 04-computation/jacobian_y0_w003_scalar02_euler_thm3714.py
python3 -B -O 04-computation/jacobian_y0_w003_scalar02_euler_thm3714.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
