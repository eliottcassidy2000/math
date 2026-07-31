---
id: THM-2991
title: "PF-infinity arbitrarily delayed global Newton-ratio return"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED + REPAIRED
source: codex-gmc-pf-infinity-delayed-return-2026-07-30
audit: >
  The promotion audit correctly verified the original two-cluster directional
  turn, but a post-promotion scope audit caught that this did not return below
  the leading edge.  MISTAKE-335 records the distinction.  The repaired
  three-cluster proof was independently rederived, including the complete
  leading ladder, central growth, reciprocal last-edge limit, strict gap,
  PF-infinity/Hurwitz/ULC typing, and exact controls through n=60.
depends_on: []
related:
  - THM-2982-first-gap-wall-stripped-norm-core-strict-ulc-through-thirty-four
  - THM-2989-first-gap-wall-stripped-all-width-leading-edge-positivity
script: 04-computation/gmc_pf_infinity_arbitrarily_delayed_newton_ratio_return_thm2991.py
output: 05-knowledge/results/gmc_pf_infinity_arbitrarily_delayed_newton_ratio_return_thm2991.out
script_sha256: 6b018ef8aef3cc7444ad95ac37275894ba76cba30819c05711d25907d23d7496
output_sha256: 1e5d43275521216768b78215a23e59f69bc962c1fca0b80ce1684afff4600b9c
hash_basis: LF-normalized bytes
---

# THM-2991 -- PF-infinity arbitrarily delayed global Newton-ratio return

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED + REPAIRED.**

## 1. Statement

For every integer `K>=1` there are integers `n>=max(2,K)` and `C>3` such
that

```text
P_(n,C)(x)=(x+1)^n(x+3)(x+C)^n                         (1)
```

has only negative real roots, is PF-infinity, and has strictly ULC
coefficients.  Put `d=2n+1` and write

```text
P_(n,C)(x)=sum_(k=0)^d e_k x^(d-k),
b_k=e_k/binom(d,k),
R_k=b_k^2/(b_(k-1)b_(k+1)),       1<=k<d.              (2)
```

The parameters may be chosen so that

```text
R_1<R_2<...<R_n                 and                 R_(2n)<R_1.  (3)
```

Thus at least `K` leading circuits can improve strictly while a later
circuit is genuinely tighter than the leading edge.  In particular, no
fixed leading Newton-ratio prefix, even together with PF-infinity, Hurwitz
stability, and strict ULC, implies the global no-return property
`R_j>=R_1` for every deeper circuit.

This is a structural no-go, not a counterexample to the first-gap norm-core
family.  A family-specific no-return theorem may use its wall invoice,
response geometry, or another continuation sidecar absent from `(1)`.

## 2. An arbitrarily long leading ladder

The `e_k` in `(2)` are the elementary symmetric functions of `n` copies of
`1`, one copy of `3`, and `n` copies of `C`.  For `0<=k<=n`,

```text
b_k(C)=c_k C^k+O(C^(k-1)),
c_k=binom(n,k)/binom(d,k).                              (4)
```

For `1<=k<n`, therefore,

```text
lim_(C->infinity) R_k(C)
 =c_k^2/(c_(k-1)c_(k+1))
 =(n-k+1)(d-k)/((d-k+1)(n-k))
 =1+(n+1)/((d-k+1)(n-k))=:L_k.                         (5)
```

The denominator in the last expression strictly decreases with `k`, so

```text
L_1<L_2<...<L_(n-1).                                   (6)
```

At `k=n`, the three relevant raw leading terms are

```text
e_(n-1)~n C^(n-1),    e_n~C^n,    e_(n+1)~(n+3)C^n.    (7)
```

The adjacent binomial quotients give

```text
lim_(C->infinity) R_n(C)/C=1/((n+3)(n+2))>0.           (8)
```

Hence `R_n(C)` tends to infinity.  For each fixed `n`, `(5)--(8)` imply
all the finitely many strict inequalities `R_1<...<R_n` once `C` is large.

## 3. The opposite edge returns below the leading edge

For any positive root multiset `r=(r_1,...,r_d)`, complementation gives

```text
e_(d-j)(r)=e_d(r)e_j(r^(-1)),
R_(d-j)(r)=R_j(r^(-1)).                                 (9)
```

Apply `(9)` with `j=1`.  As `C` tends to infinity, the reciprocal multiset
in `(1)` tends to `n` copies of `1`, one copy of `t=1/3`, and `n` zeros.
Directly from its first two elementary symmetric functions,

```text
lim R_(2n)(C)=2(n+t)^2/(d(n-1+2t)),                    (10)
lim R_1(C)   =2n^2/(d(n-1)).                           (11)
```

Their difference is

```text
lim (R_1-R_(2n))
 =2t(2n-t(n-1))/(d(n-1)(n-1+2t))>0,                   (12)
```

because `n>=2` and `t=1/3`.  The ratios are rational functions of `C`
with positive denominators for `C>0`.  The strict limiting ladder and gap
therefore hold simultaneously for every sufficiently large real `C`, and
hence for a sufficiently large integer `C>3`.  Taking
`n>=max(2,K)` proves `(3)`.

The certified return in `(3)` occurs at the last circuit; the theorem does
not claim this is the first circuit below `R_1`.  It only needs the sharp
global fact that some arbitrarily deep circuit beats the leading edge.

## 4. PF-infinity, Hurwitz, and strict ULC do not repair it

Every factor `(x+r)`, `r>0`, has a two-term PF-infinity coefficient sequence;
convolution preserves total nonnegativity of the Toeplitz matrix.  Therefore
`P_(n,C)` is PF-infinity.  Its roots are literally `-1,-3,-C`, so it is
Hurwitz stable.

Newton's inequalities say that the binomially normalized elementary
symmetric sequence `b_k` is log-concave.  Equality at any internal circuit
for positive roots forces all roots to be equal.  Since `1,3,C` are not all
equal, every circuit is strict.  Thus `(3)` occurs in the strict ULC cone,
not at a coefficient zero or stability boundary.

The smallest frozen control is

```text
P_(2,20)=x^5+45x^4+607x^3+2283x^2+2920x+1200,
(b_0,...,b_5)=(1,9,607/10,2283/10,584,1200),
(R_1,...,R_4)=(810/607,368449/205470,
               5212089/3544880,42632/34245),           (13)
```

so `R_1<R_2` and `R_4<R_1` exactly.

## 5. Correction boundary and preserved information

The first promoted version used the reciprocal two-cluster family
`(x+1)^n(x+B)^n`.  It rigorously proved an arbitrarily delayed *directional*
turn `R_(n+1)<R_n`, but reciprocal symmetry also gives
`R_(n+1)=R_(n-1)>R_1` when `n>2`.  It therefore did not refute the stronger
leading-edge no-return property.  MISTAKE-335 records the failed implication;
Sections 2--3 supply the repaired three-scale global return.

The exact connection is now

```text
source:     a finite leading list R_1,...,R_K;
target:     the global Newton-circuit ratio path;
preserved:  positivity, PF-infinity, Hurwitz stability, strict ULC,
            and an arbitrarily long improving prefix;
destroyed:  the reciprocal tail scale and the unobserved last edge;
sidecar:    a family-specific continuation law excluding a later edge return.
```

Consequently THM-2989's all-width leading-edge sign, even if its encoded
wall invoice is proved unconditionally, is not by itself an all-width ULC
theorem.  The missing no-return statement must use more than a fixed finite
top jet or a generic stability class.

## 6. Exact evidence

The standalone companion verifies `(4)--(13)` by exact rational arithmetic,
including concrete three-scale witnesses through `n=60`, the reciprocal
last-edge identity, every strict ULC circuit in those witnesses, and the
degree-five control.  It also freezes the original two-cluster family as a
directional-only boundary.  Normal and optimized transcripts are
LF-identical (13 lines, 844 LF bytes).  Reproduce with

```text
python 04-computation/gmc_pf_infinity_arbitrarily_delayed_newton_ratio_return_thm2991.py --output .scratch/thm2991.global.normal.out
python -O 04-computation/gmc_pf_infinity_arbitrarily_delayed_newton_ratio_return_thm2991.py --output .scratch/thm2991.global.opt.out
```

Frozen LF hashes are

```text
script  6b018ef8aef3cc7444ad95ac37275894ba76cba30819c05711d25907d23d7496
output  1e5d43275521216768b78215a23e59f69bc962c1fca0b80ce1684afff4600b9c
```

The repaired proof was independently rederived before canonization.  The
audit checked the reciprocal indexing, the central normalization, the strict
limit gap, the family scope, and the distinction between a directional turn
and a return below `R_1`; no defect remained.

**QED.**
