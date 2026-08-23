---
id: THM-3717
title: "W003 scalar-triple square-law Euler nonentry"
status: >
  PROVED + VERIFIED-EXACT.  The complete post-parity W003 placement family
  whose scalar fibre is the triple 03+12+20 is Darboux-empty in the y=0
  collision ring at every scale n>=2.  An upper double produces a square law;
  the lower double then gives either a nonzero Euler factor or an impossible
  constant polynomial product.  The n=2 square boundary forces an H^2 arm.
  This closes residual family C, not residual family D, all of W003, the full
  3x4 cell, general quartic C3 data, or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  SELF-AUDITED.  The exact companion checks the W003 fibre word, the upper
  integration, the middle square primitive, its logarithmic derivative, and
  the lower Euler factor in every residue class modulo 5.  It separately
  verifies the n=3 and n=4 coefficient boundaries and the n=2 square-arm
  scalar factor.  Normal and optimized runs byte-match the frozen transcript.
  Independent hostile audit remains open.
depends_on:
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
  - THM-3715-y0-w003-odd-lowest-scalar-two-orientation-half-charge-nonentry
related:
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
  - THM-3714-y0-w003-scalar-fibre-02-two-orientation-euler-nonentry
script: 04-computation/jacobian_y0_w003_scalar_triple_thm3717.py
output: 05-knowledge/results/jacobian_y0_w003_scalar_triple_thm3717.out
script_sha256: 8edf85fbbc02602fabef10edba28c591d275265d8ab94957571c8f97f815e639
output_sha256: 4fa2f1024a83706b3003fbe208f77f1cd48b03d3c7d6bf13f7e896d64ee6510e
hash_basis: LF-normalized bytes
---

# THM-3717 -- the W003 scalar triple is empty

**PROVED + VERIFIED-EXACT.**  Work over `C` in the THM-3696 collision ring.
Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

All coefficient functions below lie in `C[b]`; primes mean `d/db`.

## 1. Exact placement and singleton normalization

THM-3715 leaves the following family as one of the two exact W003 residuals:

```text
wt(P)=(-n-2,-2,2n-2),
wt(Q)=(1-2n,1-n,1,n+1),                     n>=2.       (2)
```

Its seven W003 fibres are

```text
00; 01+10; 02+11; 03+12+20; 13+21; 22; 23,             (3)
```

and `03+12+20` is the scalar fibre.  All seven displayed coefficients are
active: otherwise the support drops to an inherited smaller-support gate.

Set

```text
delta=gcd(n+2,2n-1)=gcd(n+2,5),
alpha=(n+2)/delta,             beta=(2n-1)/delta.       (4)
```

The singleton rows and the UFD common-power law give nonzero constants
`a,c,d,e,t` and nonconstant `H,K in C[b]` such that

```text
f0=a H^alpha,       g0=c H^beta,
f2=d K^(2n-2),      g2=e K,             g3=t K^(n+1).
                                                               (5)
```

THM-3696 membership gives `h|H` and `b|K`.  In particular, for every positive
integer `q`,

```text
E_q=H'K+qHK'                                             (6)
```

is nonzero: its leading coefficient is
`(deg H+q deg K) lc(H)lc(K)`.  It is in fact a nonunit of degree
`deg H+deg K-1>=2`.

Write `M=f1` at weight `-2` and `L=g1` at weight `1-n`.

## 2. The upper double fixes the hidden transport

The zero row `13+21` is

```text
0 = W_(-2,n+1)(M,tK^(n+1))
    +W_(2n-2,1-n)(dK^(2n-2),L)

  = K^(n-1)[(n+1)t(K^2M)'-2(n-1)d(K^(n-1)L)'].          (7)
```

Hence

```text
K^(n-1)L = lambda + rho K^2M,
rho=(n+1)t/[2(n-1)d].                                   (8)
```

At a complex root of the nonconstant polynomial `K`, both nonconstant terms
in `(8)` vanish.  Thus `lambda=0`.  For `n>=3`, put `Y=L/rho`; then

```text
M=K^(n-3)Y,             L=rho Y,
T=K^(n-1)Y=K^2M.                                        (9)
```

## 3. The middle double is an exact square law

Using `(4)` and `(9)`, the zero row `02+11` satisfies

```text
K^(n+1)[W_(-n-2,1)(aH^alpha,eK)
        +W_(-2,1-n)(K^(n-3)Y,rho Y)]

 = [ae H^alpha K^(n+2) + rho(3-n)T^2/2]'.              (10)
```

At `n=3`, the hidden bracket in `(10)` is identically zero, whereas the outer
bracket is

```text
ae E_5 != 0.                                             (11)
```

So the family is empty at `n=3`.

Now suppose `n>3`.  Integrating `(10)` gives

```text
T^2 = mu + [2ae/(rho(n-3))] H^alpha K^(n+2).            (12)
```

Because `T=K^(n-1)Y`, evaluation at a root of `K` forces `mu=0`.  Logarithmic
differentiation in the fraction field gives

```text
T'/T = (alpha/2)(H'/H+delta K'/K).                      (13)
```

This is the useful invariant hidden by the original coefficient equations:
the middle row does not merely determine a coefficient; it synchronizes the
`H` and `K` charges of every remaining row.

## 4. The lower double is impossible

Substitute `(9),(13)` in the zero row `01+10`.  Direct collection gives

```text
0 = W_(-n-2,1-n)(aH^alpha,rho Y)
    +W_(-2,1-2n)(K^(n-3)Y,cH^beta)

  = E_delta T/2 [a rho alpha(4-n)H^(alpha-1)K^(-n)
                 -c beta(n-2)H^(beta-1)K^(-3)].         (14)
```

The factors `E_delta` and `T` are nonzero.  At `n=4`, the first term in the
bracket vanishes and the second does not, a contradiction.

For `n>4`, cancellation in the fraction field turns `(14)` into

```text
H^(beta-alpha) K^(n-3)
  = [a rho alpha(4-n)]/[c beta(n-2)].                   (15)
```

But

```text
beta-alpha=(n-3)/delta>0.                               (16)
```

The left side of `(15)` is a nonconstant polynomial while the right side is a
nonzero scalar.  This closes every `n>4`.

## 5. The n=2 square boundary

At `n=2`, `(8)` instead reads

```text
L=rho KM,                 rho=3t/(2d).                  (17)
```

Put `S=K^2M`.  The middle double has the exact primitive

```text
K^3[W_(-4,1)(aH^4,eK)+W_(-2,-1)(M,rho KM)]
  = [ae H^4K^4 + rho S^2/2]'.                           (18)
```

As above, integration and root evaluation give

```text
S^2 = -(2ae/rho)H^4K^4.                                 (19)
```

Unique factorization over `C[b]` therefore gives a nonzero scalar `sigma`
with

```text
M=sigma H^2,             L=rho sigma H^2K,
sigma^2=-2ae/rho.                                         (20)
```

The scalar triple now factors exactly as

```text
W_(-4,3)(aH^4,tK^3)+W_(-2,1)(sigma H^2,eK)
 +W_(2,-3)(dK^2,cH^3)

 = 2H E_1(6atH^2K^2-3cdHK+e sigma).                    (21)
```

It is divisible by the nonunit `H`, so it cannot equal one.  This closes
`n=2`.

## 6. Consequence and exact scope

Equations `(11),(14),(15),(21)` close the complete W003 placement family C
from THM-3715:

```text
P=(-n-2,-2,2n-2), Q=(1-2n,1-n,1,n+1), n>=2,
scalar fibre 03+12+20.                                  (22)
```

The only W003 family not yet closed is its family D:

```text
P=(1-3n,1-2n,1), Q=(-n-2,-2,n-2,2n-2), n>=3,
scalar fibre 13+21.                                     (23)
```

Thus this theorem does not yet close W003, arbitrary `3 x 4` supports,
general quartic C3 data, or `JC(2)`.

Run

```bash
python3 -B 04-computation/jacobian_y0_w003_scalar_triple_thm3717.py
python3 -B -O 04-computation/jacobian_y0_w003_scalar_triple_thm3717.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
