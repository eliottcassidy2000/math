---
id: THM-3724
title: "W004 scalar-13 Kummer-twist nonentry"
status: >
  PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.  Every W004
  placement in the named all-scale family with weights
  P=(1-3n,1-2n,1), Q=(-2n-2,-n-2,-2,2n-2) and scalar fibre 13+22 is
  Darboux-empty in the y=0 collision ring.  The end rows create a Kummer
  power twist; the two middle rows are differentially incompatible on both
  its zero and nonzero sheets.  This is not proved canon until audit
  promotion, and it does not close all of W004, the full 3x4 cell, general
  quartic C3 data, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The exact companion checks the W004 fibre
  word, singleton arithmetic in every residue modulo 8, both transport
  identities, the zero-sheet Cramer compatibility, the nonzero Kummer-sheet
  determinant and three-term compatibility numerator, and the n=1 boundary.
  Normal and optimized runs byte-match the frozen transcript.  Independent
  hostile audit remains open.
depends_on:
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3722-y0-w004-scalar-03-norm-twist-nonentry
script: 04-computation/jacobian_y0_w004_scalar13_kummer_twist_thm3724.py
output: 05-knowledge/results/jacobian_y0_w004_scalar13_kummer_twist_thm3724.out
script_sha256: 78b0e627dbaedf438961cd4593c639709cb6a2e8a8fef08dc6e00ceb3203b742
output_sha256: 181bfafb7a6719e0bb3e3ee6cce2b1741d83cd6bd3382989f1cbd9a6717555ce
hash_basis: LF-normalized bytes
---

# THM-3724 -- a complete W004 Kummer-twist family is empty

**PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.**  Work over
`C` in the THM-3696 collision ring.  Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

All coefficient functions lie in `C[b]`; primes mean `d/db`.

## 1. Exact family and singleton bases

The W004 ray has supports `(0,n,3n)` and `(0,n,2n,4n)`, hence fibres

```text
00; 01+10; 02+11; 12+20; 03+21; 13+22; 23.            (2)
```

We close the following actual-support placement at every `n>=1`:

```text
wt(P)=(1-3n,1-2n,1),
wt(Q)=(-2n-2,-n-2,-2,2n-2),
scalar fibre 13+22.                                    (3)
```

First suppose `n>=2`.  Set

```text
delta=gcd(3n-1,2n+2)=gcd(n-3,8),
alpha=(3n-1)/delta,             beta=(2n+2)/delta,
R=2n-1,                         T=2n-2,
m=2alpha-beta=4(n-1)/delta.                            (4)
```

The singleton rows `00` and `23` give nonzero constants `a,c,d,t` and
nonconstant `H,K in C[b]` with

```text
f0=aH^alpha,          g0=cH^beta,
f2=dK,                g3=tK^T.                         (5)
```

THM-3696 membership gives `h|H` and `b|K`.  In particular both bases are
nonconstant.  Put

```text
U=HK^delta,           Pi=K^R f1,
L=g1,                 N=g2,              Z=K^2N.       (6)
```

## 2. The two end transports create a Kummer sheet

The upper zero row `03+21` is exactly

```text
(K^(n+2)L)'=(atT/d)(U^alpha)'.                         (7)
```

At any root of the nonconstant polynomial `K`, both sides before
differentiation vanish after integration.  Thus the integration constant is
zero and, with `rho=atT/d`,

```text
L=rho K^(-(n+2))U^alpha=rho H^alpha K^(R-2).           (8)
```

The lowest zero row `01+10` now reduces to

```text
delta U Pi'-R Pi U'=D U^m U',
D=a rho alpha(R-2)/(c beta).                           (9)
```

Since `delta*m-R=R-2`, put

```text
A=a rho alpha/(c beta),               F=Pi-AU^m.       (10)
```

Equation `(9)` becomes

```text
delta U F'-R F U'=0,
hence                         F^delta=kappa U^R         (11)
```

for some `kappa in C`; `kappa=0` includes `F=0`.  This is the Kummer
analogue of THM-3722's quadratic norm sheet.  It is not itself a
contradiction: both its zero and nonzero sheets must be tested against the
two middle rows.

## 3. Exact middle-row system

After substituting `(6)` and `(8)`, the rows `12+20` and `02+11` are

```text
R Pi Z'-2Pi'Z=cd beta U^(beta-1)U',                    (12)

a alpha(delta U Z'-2U'Z)
 +rho(R alpha Pi U'-(n+2)U Pi')=0.                    (13)
```

These two equations remember both algebraic value and derivative data.  The
rest of the proof shows that no sheet of `(11)` can satisfy both.

## 4. The zero Kummer sheet is differentially incompatible

Let `kappa=0`, so `Pi=AU^m`, and put

```text
e=beta-m=2(3-n)/delta,                V=Z'/U'.          (14)
```

Treating `(12)--(13)` as a linear system for `(V,Z)`, its determinant is

```text
2a alpha A(R-2)U^m !=0.                                (15)
```

Cramer's rule gives

```text
Z= -[R rho A(n-3)/(2a alpha delta)]U^m
   -[delta cd beta/(2A(R-2))]U^e,

V= -[rho A m(n-3)/(a alpha delta)]U^(m-1)
   -[cd beta/(A(R-2))]U^(e-1).                         (16)
```

But `V=Z'/U'` requires `dZ/dU=V`, whereas `(16)` gives

```text
dZ/dU-V=
 [cd beta(n-2)/(A(R-2))]U^(e-1)
 -[rho A m(n-3)(R-2)/(2a alpha delta)]U^(m-1).         (17)
```

The two exponents are distinct: `e=m` would imply `6n=10`.  At `n=2` the
second coefficient in `(17)` is nonzero; at `n=3` the first is nonzero; and
at every `n>=4` both are nonzero.  Since a nonconstant polynomial `U` is
transcendental over `C`, distinct Laurent powers of `U` are linearly
independent.  Thus `(17)` cannot vanish.

## 5. Every nonzero Kummer sheet is also incompatible

Now let `kappa!=0`.  The integer `delta` divides `8`, while `R` is odd, so
`gcd(delta,R)=1`.  Unique factorization in `C[b]`, followed by a harmless
scalar rescaling over `C`, turns `(11)` into

```text
U=v^delta,                F=Bv^R,              B!=0,   (18)
```

for a nonconstant polynomial `v`.  Since `delta*m=4n-4`,

```text
Pi=Av^(4n-4)+Bv^R.                                     (19)
```

Put `Y=Z'/v'`.  In the variables `(Y,Z)`, equations `(12)--(13)` have
determinant

```text
2a alpha delta A(2n-3)v^(4n-4) !=0.                   (20)
```

Solving algebraically for `(Y,Z)` and then imposing the necessary derivative
compatibility `dZ/dv=Y` gives the exact rational function

```text
dZ/dv-Y = N(v)/D(v),                                  (21)

D(v)=2ac rho(n+1)(2n-3)(3n-1)v^(2n+5),

N(v)=
 -3Bc rho^2(n-2)(n+1)(2n-3)^2(2n-1)v^(4n+3)
 -2a rho^3(n-3)(n-1)(2n-3)^2(3n-1)v^(6n)
 +8c^3d(n-2)(n+1)^3v^10.                              (22)
```

At `n=2`, only the middle term of `N` survives, and it is nonzero.  At
`n=3`, the first and third terms survive with distinct powers.  At every
`n>=4`, all three coefficients are nonzero, and the powers `4n+3`, `6n`,
and `10` are pairwise distinct.  Injectivity of `C[s] -> C[b]`, `s |-> v`,
therefore gives `N(v)!=0`.  This contradicts `(21)` on every nonzero Kummer
sheet.

## 6. The n=1 boundary

When `n=1`, the top singleton `23` has weights `(1,0)`.  Its zero bracket is

```text
W_(1,0)(f2,g3)=-f2 g3'=0.                              (23)
```

The domain and actual support give `g3'=0`; hence the weight-zero piece `g3`
is a literal scalar.  Delete it without changing any bracket.  The remaining
pair has at most three graded pieces in each output and is impossible by the
all-degree `3 x 3` theorem THM-3592.

## 7. Scope and reproduction

Sections 2--6 close every scale of the named W004 placement `(3)`.  The
mechanism is stronger than a scalar-row obstruction: every candidate fails
the derivative compatibility of its zero rows before reaching `13+22`.
Other W004 placements, W005--W006, arbitrary `3 x 4` supports, unrestricted
quartic C3/cofactor data, and `JC(2)` remain open.

Run

```bash
python3 -B 04-computation/jacobian_y0_w004_scalar13_kummer_twist_thm3724.py
python3 -B -O 04-computation/jacobian_y0_w004_scalar13_kummer_twist_thm3724.py
```

Both commands must agree byte for byte with the frozen transcript.  Audit
promotion is required before this proof candidate enters proved canon.
