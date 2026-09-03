---
id: THM-4362
title: "Tetrahedral diagonal annihilator and sharp depth threshold"
status: >
  PROVED FINITE DEPTH-MODULE THEOREM RELATIVE TO THM-4308 + VERIFIED-EXACT
  + INDEPENDENTLY AUDITED. For every m>=7 and d>=0, the primitive alternating
  tetrahedral functional L_m on the coefficient diagonal (n,2n-10),
  5<=n<=m, annihilates pi_m(P_d) if and only if d<=m-7. The positive side is
  a fourth finite-difference identity applied generator by generator; at
  d=m-6 an explicit retained monomial has value +/-1, and nesting gives
  failure at every greater depth. Its m=9 and m=10 instances are exactly the
  functionals used in THM-4358 and THM-4361. No bracket, all-row membership,
  seam-entry, Keller-pair, JC(2), or DC(2) conclusion is asserted.
source: root + clean-room finite-difference referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4308-source-normal-bracket-hasse-truncation-through-row-eight
related:
  - THM-4358-source-normal-s4339-row-ten-delayed-depth-extinction
  - THM-4361-source-normal-beta-zero-row-ten-joint-depth-extinction
  - THM-3417-exchangeable-krawtchouk-shell-inversion-and-radial-high-minor-converse
mistake_firewall:
  - MISTAKE-522
primary_script: 04-computation/tetrahedral_diagonal_depth_annihilator_thm4362.py
primary_output: 05-knowledge/results/tetrahedral_diagonal_depth_annihilator_thm4362.out
primary_script_sha256: 0d2d1aa8a21460c02850edcab4c0d2be769a1b1c2aff925613459f94a73663c7
primary_output_sha256: 3de17aff44166a5ae4841624a6b23547089d11a676d6d9e628a8f6f5055fd9da
independent_referee_script: 04-computation/tetrahedral_diagonal_depth_annihilator_independent_referee_thm4362.py
independent_referee_output: 05-knowledge/results/tetrahedral_diagonal_depth_annihilator_independent_referee_thm4362.out
independent_referee_script_sha256: 2d11c433436dc93159f8d867617cd72aeece88ec6826be8970f0df0ed15f99d2
independent_referee_output_sha256: d65b7cc6f0f2d0b2b14bd4a1b49d3effbf7d0990a0d63b1a15d56c131138f6dd
hash_basis: raw LF bytes
audit: >
  PASS. The primary enumerates 85,281 exact projected columns in a hostile
  finite range and checks explicit boundary generators. An import-free
  clean-room implementation separately proves and tests the support,
  endpoint, partial-binomial, sign, nesting, and figurate identities, with
  2,639,852 grouped exact checks. Normal/-O/frozen LF streams match.
---

# THM-4362 -- Tetrahedral diagonal annihilator and sharp depth threshold

**PROVED FINITE DEPTH-MODULE THEOREM RELATIVE TO THM-4308 + VERIFIED-EXACT
+ INDEPENDENTLY AUDITED. THIS IS A SHARP STATEMENT ABOUT `pi_m(P_d)`. IT
DOES NOT ASSERT A BRACKET SOLUTION, INFINITE DEPTH MEMBERSHIP, AN ALL-ROW
LIFT, SEAM ENTRY, A KELLER PAIR, `JC(2)`, OR `DC(2)`.**

## 1. Statement and inheritance

Work over the algebraically closed characteristic-zero field and the source-
normal coordinates of THM-4308. For a projected polynomial `h`, write

```text
h_(n,r)=[x^r t^n]h.
```

Here `pi_m` retains exactly the coefficient rows `0<=n<=m`.

For every integer `m>=7`, define

```text
L_m(h)=sum_(n=5)^m
       (-1)^(n-5) C(m+3-n,3) h_(n,2n-10).               (1)
```

The last coefficient is `+/-C(3,3)=+/-1`, so `L_m` is primitive. The sharp
theorem is

```text
L_m vanishes identically on pi_m(P_d)  iff  d<=m-7       (2)
```

for every integer `d>=0`.

The inheritance pass was:

- closest proved mechanism: THM-4308's exact monomial spanning and source-
  normal expansion laws;
- canonical hostile: an annihilator observed in one finite matrix need not
  persist one depth farther;
- corrected near miss: a rank/nullspace pattern becomes a theorem only after
  every generator and both truncation endpoints are controlled;
- least-used sidecar: the three omitted rows carrying zero consecutive
  tetrahedral weights.

The live concept board was

```text
depth generator | coefficient diagonal | binomial packet length
truncation buffer | finite-difference order | unit hostile | nesting.       (3)
```

## 2. One generator meets all or none of the diagonal

THM-4308 proves

```text
P_d=span{x^a u^b p^c y^e: a,b,c,e>=0, a+b<=d},          (4)
```

and each generator expands as

```text
x^a u^b p^c y^e
=x^(a+2b+e)t^(b+c+2e)(1+x^2t)^(c+e).                   (5)
```

Set

```text
r0=a+2b+e,               n0=b+c+2e,               N=c+e.
```

The term indexed by `0<=k<=N` in `(5)` is

```text
C(N,k)x^(r0+2k)t^(n0+k).                                (6)
```

It meets the diagonal `r=2n-10` exactly when

```text
r0+2k=2(n0+k)-10  iff  a=2c+3e-10.                     (7)
```

The condition is independent of `k`. Thus a generator misses every
coordinate in `(1)`, or its complete binomial packet lies on that diagonal.

Assume `(7)`. Since `a>=0`,

```text
2c+3e>=10,
N=c+e>=4,                         n0>=c+2e>=5.          (8)
```

The final row of the packet is

```text
nmax=n0+N=b+2c+3e=b+a+10.                               (9)
```

These three elementary bounds contain the whole mechanism: at least four
binomial steps are present, none lies below row five, and depth controls the
upper endpoint.

## 3. Fourth-difference annihilation

Suppose `a+b<=d<=m-7`. Equation `(9)` gives

```text
nmax<=m+3.                                               (10)
```

The only packet terms above the retained row `m` can be at `m+1,m+2,m+3`.
Their continued weights in `(1)` are

```text
C(2,3)=C(1,3)=C(0,3)=0.                                 (11)
```

Hence extending the sum from the retained terms to all `k=0,...,N` changes
nothing. Put `R=m+3-n0`. By `(9)--(10)`, `R>=N`. Ordinary coefficient
extraction, with no negative-top binomial convention, gives

```text
L_m(generator)
=(-1)^(n0-5) sum_(k=0)^N (-1)^k C(N,k)C(R-k,3)
=(-1)^(n0-5)[s^3]
  (1+s)^R(1-(1+s)^(-1))^N
=(-1)^(n0-5)[s^3]s^N(1+s)^(R-N)
=0,                                                        (12)
```

because `N>=4`. Linearity and `(4)` prove the forward implication in `(2)`.
Equivalently, the tetrahedral weight is a cubic polynomial in the backward
row address, and a packet of binomial order `N>=4` (at least five terms)
takes its fourth finite difference.

## 4. The one-depth-beyond unit hostile

Let

```text
d0=m-6,                 epsilon=d0 mod 2,
a=d0, b=0, e=epsilon,   c=(10+d0-3epsilon)/2.            (13)
```

For `m>=7`, these are nonnegative integers, `a+b=d0`, and `(7)` holds. Here
`N>=5` and

```text
n0+N=d0+10=m+4.                                        (14)
```

The retained terms are precisely `0<=k<=N-4`. With `R=N-1`, their unsigned
partial difference is

```text
S_N=sum_(k=0)^(N-4)(-1)^k C(N,k)C(N-1-k,3).
```

The next three ordinary binomial coefficients vanish. Therefore

```text
S_N=[s^3]sum_(k=0)^(N-1)(-1)^k C(N,k)(1+s)^(N-1-k)
   =[s^3](s^N-(-1)^N)/(1+s)
   =(-1)^N.                                              (15)
```

Combining `(14)--(15)`,

```text
L_m(x^a p^c y^e)=(-1)^(n0-5+N)
                 =(-1)^(m-1)=(-1)^(d0+7).              (16)
```

This is a unit, not merely a nonzero rational with possible characteristic
or denominator issues. For every `d>=d0`, the same monomial belongs to
`P_d` by `(4)`, so failure persists by nesting. This proves the reverse
implication and the sharpness of `(2)`.

## 5. The row-nine and row-ten consumers

At `m=9`, `(1)` is

```text
35h_(5,0)-20h_(6,2)+10h_(7,4)-4h_(8,6)+h_(9,8).        (17)
```

The threshold is `d<=2`; on `h=A`, `(17)` is exactly THM-4358's primitive
`P_2` functional `H_A`.

At `m=10`, `(1)` is

```text
56h_(5,0)-35h_(6,2)+20h_(7,4)-10h_(8,6)
 +4h_(9,8)-h_(10,10).                                  (18)
```

The threshold is `d<=3`; on `h=C`, `(18)` is exactly THM-4361's primitive
`P_3` functional `H_C`. The present theorem proves that the two rows
annihilate their whole projected depth modules. It does not compute their
values on a bracket-selected jet; those extinction arguments remain the
content of THM-4358 and THM-4361.

## 6. Triangular numbers and natural addresses

For `r>=3`,

```text
C(r,3)=sum_(q=1)^(r-2) T(q),              T(q)=q(q+1)/2. (19)
```

Thus the absolute weights in `(1)` are cumulative triangular numbers. On
the integer extension,

```text
T(-z)=T(z-1),              T(z+2)-T(z-2)=2+4z.          (20)
```

Equation `(20)` is the centered first-layer identity; the proof of `(2)` is
the corresponding fourth-difference cancellation for the cumulative cubic
weights. The natural discrete address is the backward row distance
`m-n` together with tetrahedral degree three (annihilating difference order
four). Replacing these by an arbitrary ordinal label would forget the
truncation endpoint that determines
the sharp depth threshold.

The comparison with THM-3417's Krawtchouk inversion stops at signed binomial
transforms: the source, diagonal, and depth endpoint in `(7)--(10)` are
indispensable here. No tournament or LRC statement follows from that shared
syntax.

## 7. Audit and scope

The 86,508-check primary verifies the finite-difference core, 85,281 exact
projected columns for `7<=m<=18`, 31 sharp hostiles, both downstream
stencils, and the figurate identities. The independent 2,639,852-check
implementation separately rebuilds the columns, partial-sum hostile, signs,
endpoints, nesting, and a larger diagonal/order test bank. Both normal and
optimized streams byte-match their frozen LF outputs.

The argument after `(4)--(5)` is integral and would work over any nonzero
commutative coefficient ring with that module presentation. The identification
of the actual `P_d` with `(4)`, however, is imported only in THM-4308's
algebraically closed characteristic-zero setting; this theorem does not
silently enlarge that geometric scope.

Reproduce from the repository root:

```text
python3 -B 04-computation/tetrahedral_diagonal_depth_annihilator_thm4362.py
python3 -B -O 04-computation/tetrahedral_diagonal_depth_annihilator_thm4362.py
python3 -B 04-computation/tetrahedral_diagonal_depth_annihilator_independent_referee_thm4362.py
python3 -B -O 04-computation/tetrahedral_diagonal_depth_annihilator_independent_referee_thm4362.py
```
