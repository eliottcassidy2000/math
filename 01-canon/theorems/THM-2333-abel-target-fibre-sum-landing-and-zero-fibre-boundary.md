---
id: THM-2333
title: "Abel target-fibre sum landing and the zero-fibre boundary"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. On every
  positive strict shallow-owner word stratum covered by THM-2327, the
  zero mode of the extra word factor and the two nonzero bare endpoint
  coefficients form a nonzero convolution slice of the marked
  word/deepest-comb/bare current. Under THM-2309's target quotient, the
  base-index Abel expansion of this slice decomposes into 169
  target-vector fibre sums, every fibre sum has an Abel limit, and at
  least one limit is nonzero. The surviving fibre may be zero. This is
  sharp at the finite-group level: exact rational full-support endpoint
  weights can put 169 nonzero term pairs in every fibre while cancellation
  leaves only the zero fibre. The theorem does not retain the nonconstant
  word mode, force a nonzero target gain, give an all-91-unit or bounded
  visible survivor, or transport terminal phase. No scalar row is
  excluded and LRC(14) remains open.
source: codex-2026-07-25-abel-target-fibre-sum
depends_on:
  - THM-2302-same-label-expiration-dichotomy-and-pure-terminal-shell-no-go
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2327-two-colour-marked-unit-c3-triangle
related:
  - THM-2319-crt-unit-bispectrum-needle-and-mixed-polarization-no-go
  - THM-2329-boundary-triple-rerooting-and-transverse-gain-obstruction
  - THM-2331-two-sided-septimal-address-embedding-in-marked-current
script: 04-computation/lrc14_abel_target_fibre_sum_thm2333.py
output: 05-knowledge/results/lrc14_abel_target_fibre_sum_thm2333.out
script_sha256: db9a490bd74e7905db26e2d59e591c980c81c7ac533975cd8c638d74cdb3a6ac
output_sha256: 2f5e64410a82476fedb6b54303b67554698ba3fe5791cc321f69cabf6434fb59
hash_basis: working-tree bytes (LF)
---

# THM-2333 -- one target fibre survives, but it may be zero

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2331 embeds every prescribed target address as a nonzero term of the
marked current's Abel expansion. It deliberately does not say that terms
in one target fibre have a nonzero sum. The first aggregate statement is
available, but its sharp boundary is the zero target fibre.

The operation is finite-character projection:

```text
expanded harmonic pair (u,v)
  -> exact relation r=u+m e_3-v
  -> target quotient q=pi(r) in K_13/L_13
  -> sum all terms with the same q.                 (1)
```

There are only `169` fibres. The sum over all of them is a nonzero
zero-word-mode convolution slice, so one fibre survives. Nothing in this
argument makes that fibre nonzero in the target plane.

## 1. A nonzero zero-word-mode slice

Use THM-2327 on a positive strict shallow-owner word stratum. Its exact
mixed triangle is

```text
(1_(E_Q))_hat(X)
(1_(D_(c_3)))_hat(m c_3)
conjugate((1_E)_hat(Y)) !=0,

Y=X+m c_3,                       gcd(m,91)=1,       (2)
```

where

```text
(1_E)_hat(X)(1_E)_hat(Y)!=0.                        (3)
```

Write

```text
W=T^(-(lambda_j+1))Q,
E_Q=E intersection W,

omega=(1_W)_hat(0)=measure(W)=measure(Q)>0.         (4)
```

In the convolution expansion of `(1_(E_Q))_hat(X)`, retain only the
`h=0` Fourier mode of `1_W`. This defines the separate slice

```text
C_0
 =omega (1_E)_hat(X)
   (1_(D_(c_3)))_hat(m c_3)
   conjugate((1_E)_hat(Y))
 !=0.                                               (5)
```

Equation (5) is not the full current (2). Its nonzero word modes may
cancel it. The point of (5) is different: both endpoint functions are
now the same nine-factor exclusive-owner rectangle, so their relation
terms carry a common target quotient.

The strict shallow grade gives

```text
13 divides X,Y,c_3.                                 (6)
```

Thus every endpoint harmonic address below reduces into

```text
K_13=(w mod 13)^perp.                               (7)
```

## 2. Base-index Abel expansion

Write the exclusive rectangle, up to null endpoints, as

```text
1_E(t)=prod_(i=0)^8 I_i(w_i t).                     (8)
```

THM-2302 gives eight danger/safe factors and one wider guard factor.
Write

```text
a(u)=prod_i (I_i)_hat(u_i),                         (9)
```

and base-Poisson smooth each factor:

```text
(I_i)_hat(k) -> rho^|k| (I_i)_hat(k),
                         0<rho<1.                  (10)
```

The resulting absolutely convergent rectangle coefficient is

```text
F_N(rho)
 =sum_(u.w=N) a(u) rho^(||u||_1).                  (11)
```

Composition with each integer speed preserves Haar `L^1` convergence, so

```text
lim_(rho->1-)F_N(rho)=(1_E)_hat(N).                 (12)
```

Smooth the base deepest-comb factor by `rho^|m|` and put

```text
C_0(rho)
 =omega rho^|m| d_hat(m)
    F_X(rho) conjugate(F_Y(rho)),                   (13)

d_hat(m)=sin(pi*m/7)/(pi*m)!=0.
```

Then

```text
lim_(rho->1-)C_0(rho)=C_0!=0.                       (14)
```

This regularizes the slice (5), not the full word current (2).

## 3. The 169 target-fibre sums

Use THM-2309's exact decomposition

```text
K_13=L_13 direct_sum span(e_a,e_b),

G:=K_13/L_13 isomorphic to F_13^2,       |G|=169,  (15)
```

and write `pi:K_13->G` for the quotient map. For `q in G`, define

```text
A_q(rho)
 =omega rho^|m| d_hat(m)
   sum_(
     u.w=X, v.w=Y,
     pi(u+m e_3-v)=q
   )
     a(u) conjugate(a(v))
     rho^(||u||_1+||v||_1).                         (16)
```

Every sum is absolutely convergent. The relation inside (16) is exact:

```text
(u+m e_3-v).w=X+m c_3-Y=0.                         (17)
```

Equations (6)--(7) make each quotient label in (16) well-defined.
Partitioning all pairs by their label gives

```text
sum_(q in G) A_q(rho)=C_0(rho).                     (18)
```

It remains to justify passing to the Abel boundary separately in every
fibre; (18) alone would not do so for infinitely many terms.

## 4. Finite-character projection gives every limit

Let `G^` be the `169` characters of `G`. For `chi in G^`, define the
twisted rectangle coefficient

```text
F_N^chi(rho)
 =sum_(u.w=N)
    a(u)rho^(||u||_1) chi(pi(u)).                  (19)
```

This is well-defined because `13|N` makes `u mod 13` lie in `K_13`.
Extend `chi after pi` from `K_13` to a linear character on `F_13^9`.
If its coordinate exponents are `tau_i`, then (19) is the ordinary
frequency-`N` coefficient of

```text
prod_i I_(i,rho)(w_i t+tau_i/13),                  (20)
```

up to the harmless Fourier-sign convention. A different extension changes
the exponent by a multiple of `w`; it is trivial on `u.w=N=0 mod 13`.

Every shifted factor in (20) converges in `L^1`, and the factors are
uniformly bounded Poisson averages. Hence

```text
F_N^chi(1-):=lim_(rho->1-)F_N^chi(rho)              (21)
```

exists.

Finite-character orthogonality applied to (16) gives

```text
A_q(rho)
 =1/169 sum_(chi in G^)
    conjugate(chi(q)) chi(pi(m e_3))
    omega rho^|m| d_hat(m)
    F_X^chi(rho) conjugate(F_Y^chi(rho)).           (22)
```

Therefore every limit

```text
A_q(1-)=lim_(rho->1-)A_q(rho)                      (23)
```

exists, and the finite sum commutes with the limit:

```text
sum_(q in G) A_q(1-)=C_0!=0.                       (24)
```

At least one target fibre satisfies

```text
A_q(1-)!=0.                                         (25)
```

This is a cancellation-surviving aggregate relation-address incidence
theorem. It does not yet prove `q!=0`.

## 5. A sharp full-support zero-fibre hostile

The zero-fibre possibility is not a cosmetic caveat. It survives exact
factorization and complete termwise support at the quotient level.

Let `G=F_13^2`, let `N=|G|=169`, and fix any `p in G`. In the rational
group algebra define

```text
U=delta_0+1_G,

V_0=delta_0-(1/(N+1))1_G,
V(x)=V_0(x-p).                                      (26)
```

Thus

```text
U(0)=2,             U(x)=1 for x!=0,

V(p)=N/(N+1),       V(x)=-1/(N+1) for x!=p.        (27)
```

Every value of `U` and `V` is nonzero. Define the fibre correlation

```text
H(q)=sum_(u in G)U(u)V(u+p-q).                      (28)
```

Every `q` contains all `169` nonzero term pairs. Direct calculation gives

```text
H(0)=1,
H(q)=0 for every q!=0.                              (29)
```

Indeed, at `q=0` the numerator is

```text
2N-(N-1)=N+1,
```

while for `q!=0` it is

```text
-2+N-(N-2)=0.                                      (30)
```

Equivalently, `V_0` is the convolution inverse of `U`.

This is an exact algebraic hostile, not a claim that the rational weights
in (26) are produced by a canonical LRC rectangle. It proves that the
formal inputs

```text
nonzero total
 + factorized endpoint weights
 + every quotient fibre termwise occupied
```

do not force a nonzero target vector. A future proof must use a semantic
word mode, target-specific positivity, phase, or another property absent
from arbitrary group-algebra weights.

## 6. Exact remaining boundary

The progression is now:

```text
THM-2325:
  every nonzero target vector has exact all-91-unit addresses;

THM-2331:
  every such address occurs as a nonzero fixed-rho term;

THM-2333:
  some target fibre has a nonzero Abel aggregate;

still open:
  prove that a surviving fibre is nonzero in G and carries the semantic
  function role of the positive word, retain an all-91-unit or visible
  address inside that aggregate, and transport terminal phase.        (31)
```

The zero-word-mode operation intentionally replaces the word by the scalar
`omega`. It retains positive word mass but destroys which target word was
realized. The hostile (26)--(30) shows why that destroyed coordinate is
exactly what the next theorem must restore.

No scalar profile is excluded. The exact ledger remains `165`, and
LRC(14) remains open.

## 7. Exact companion

The companion verifies the `169`-character orthogonality distribution,
the exact owner-packet quotient dimension, every rational hostile fibre
sum, all `169^2` nonzero hostile term pairs, and the convolution-inverse
identity. Every load-bearing check raises explicitly under ordinary and
optimized Python.

Reproduce with

```bash
python3 04-computation/lrc14_abel_target_fibre_sum_thm2333.py
python3 -O 04-computation/lrc14_abel_target_fibre_sum_thm2333.py
```

Both transcripts must match

```text
05-knowledge/results/lrc14_abel_target_fibre_sum_thm2333.out
```

byte-for-byte after LF normalization.
