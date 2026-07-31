---
id: THM-2991
title: "PF-infinity arbitrarily delayed Newton-ratio return"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED
source: codex-gmc-pf-infinity-delayed-return-2026-07-30
audit: >
  An independent hostile audit rederived the complete asymptotic ladder,
  central growth constant, reciprocal indexing, PF-infinity/Hurwitz typing,
  and strict Newton equality boundary.  It extended exact convolution
  controls through n=60 and replayed normal, optimized, and stored output.
  Its sole evidence request--directly checking the central leading-coefficient
  quotient--is installed before promotion.
depends_on: []
related:
  - THM-2982-first-gap-wall-stripped-norm-core-strict-ulc-through-thirty-four
  - THM-2989-first-gap-wall-stripped-all-width-leading-edge-positivity
script: 04-computation/gmc_pf_infinity_arbitrarily_delayed_newton_ratio_return_thm2991.py
output: 05-knowledge/results/gmc_pf_infinity_arbitrarily_delayed_newton_ratio_return_thm2991.out
script_sha256: 29f2bd92b13ec180badc334f47792c5c69a26865da9e2b1e1fc8208c3be9d054
output_sha256: 15a83465e78783e127453925715d61b78fa715a5b27d35555b5c142710101b2a
hash_basis: LF-normalized bytes
---

# THM-2991 -- PF-infinity arbitrarily delayed Newton-ratio return

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

For every integer `K>=1` there are integers `n>=max(2,K)` and `B>1` such
that

```text
P_(n,B)(x)=(x+1)^n(x+B)^n                              (1)
```

has only negative real roots, is PF-infinity, and has strictly ULC
coefficients.  Write

```text
P_(n,B)(x)=sum_(k=0)^(2n) e_k x^(2n-k),
b_k=e_k/binom(2n,k),
R_k=b_k^2/(b_(k-1)b_(k+1)),       1<=k<2n.             (2)
```

The parameters may be chosen so that

```text
R_1<R_2<...<R_n,
R_(2n-k)=R_k,
R_(n+1)=R_(n-1)<R_n.                                  (3)
```

Thus the first backward step occurs at circuit `n+1`, after an arbitrarily
long strictly improving leading prefix.  In particular, no fixed number of
leading Newton circuits, even together with PF-infinity, Hurwitz stability,
and strict ULC, implies the global **no-return** property.

This is a structural no-go, not a counterexample to the first-gap norm-core
family.  A family-specific no-return theorem may use its wall invoice,
response geometry, or another continuation sidecar not present in `(1)`.

## 2. Leading asymptotics

The `e_k` in `(2)` are the elementary symmetric functions of `n` copies of
`1` and `n` copies of `B`.  For `0<=k<=n`,

```text
b_k(B)=c_k B^k+O(B^(k-1)),
c_k=binom(n,k)/binom(2n,k).                             (4)
```

For `1<=k<n`, therefore,

```text
lim_(B->infinity) R_k(B)
 =c_k^2/(c_(k-1)c_(k+1))
 =(n-k+1)(2n-k)/((2n-k+1)(n-k))
 =1+n/((n-k)(2n-k+1))=:L_k.                            (5)
```

The denominator in the last expression strictly decreases with `k`, so

```text
L_1<L_2<...<L_(n-1).                                   (6)
```

At the central circuit the two adjacent leading terms are instead

```text
b_(n-1)~ n B^(n-1)/binom(2n,n-1),
b_n    ~   B^n/binom(2n,n),
b_(n+1)~ n B^n/binom(2n,n+1).                          (7)
```

Using `binom(2n,n-1)=binom(2n,n+1)` gives

```text
lim_(B->infinity) R_n(B)/B=1/(n+1)^2.                  (8)
```

Hence `R_n(B)` tends to infinity.  Equations `(5)--(8)` show that, for each
fixed `n`, all finitely many strict inequalities

```text
R_1(B)<...<R_n(B)                                       (9)
```

hold for every sufficiently large real `B`.  The ratios are rational
functions with rational coefficients and positive denominators on `B>0`, so
one may choose an integer `B>1`.  Taking `n>=max(2,K)` supplies the requested
arbitrarily long prefix.

## 3. Reciprocal symmetry forces the return

Complementing a chosen subset of the multiset of roots gives the exact
identity

```text
e_(2n-k)=B^(n-k)e_k.                                   (10)
```

Since `binom(2n,2n-k)=binom(2n,k)`,

```text
b_(2n-k)=B^(n-k)b_k.                                   (11)
```

The powers of `B` cancel in each three-term quotient, proving

```text
R_(2n-k)=R_k.                                          (12)
```

In particular, `(9)` and `(12)` give

```text
R_(n+1)=R_(n-1)<R_n.                                   (13)
```

Thus the later return is not a numerical accident: it is forced by the
same reciprocal symmetry that permits the long initial improvement.

## 4. PF-infinity, Hurwitz, and strict ULC do not repair it

Every factor `(x+r)`, `r>0`, has a two-term PF-infinity coefficient sequence;
convolution preserves total nonnegativity of the Toeplitz matrix.  Therefore
`P_(n,B)` is PF-infinity.  Its roots are literally `-1` and `-B`, so it is
Hurwitz stable.

Newton's inequalities say that the binomially normalized elementary
symmetric sequence `b_k` is log-concave.  Equality at any internal circuit
for positive roots forces all roots to be equal.  Since `B>1`, every circuit
is strict.  Thus the return in `(13)` occurs inside the strict ULC cone, not
at a zero, repeated-coefficient, or stability boundary.

The smallest exact control is

```text
P_(2,2)=x^4+6x^3+13x^2+12x+4,
(b_0,...,b_4)=(1,3/2,13/6,3,4),
(R_1,R_2,R_3)=(27/26,169/162,27/26).                   (14)
```

It is already PF-infinity and strict ULC, yet its middle improvement returns
immediately on the opposite side.

## 5. Preserved and destroyed information

The exact connection is

```text
source:     a finite leading list R_1,...,R_K;
target:     the global Newton-circuit ratio path;
preserved:  positivity, PF-infinity, Hurwitz stability, strict ULC,
            and an arbitrarily long improving prefix;
destroyed:  reciprocal position and the unobserved coefficient tail;
sidecar:    a family-specific continuation law excluding a later turn.
```

Consequently THM-2989's all-width leading-edge sign, even if its encoded
wall invoice is later proved unconditionally, is not by itself an all-width
ULC theorem.  The missing no-return statement must use more than a fixed
finite top jet or a generic stability class.

## 6. Exact evidence

The standalone companion verifies `(4)--(14)` by exact rational arithmetic,
including concrete witnesses through `n=40`, the reciprocal identities at
every coefficient, and the degree-four hostile.  Normal and optimized
transcripts are LF-identical (12 lines, 705 bytes).  Reproduce with

```text
python 04-computation/gmc_pf_infinity_arbitrarily_delayed_newton_ratio_return_thm2991.py --output .scratch/thm2991.normal.out
python -O 04-computation/gmc_pf_infinity_arbitrarily_delayed_newton_ratio_return_thm2991.py --output .scratch/thm2991.opt.out
```

Frozen LF hashes are

```text
script  29f2bd92b13ec180badc334f47792c5c69a26865da9e2b1e1fc8208c3be9d054
output  15a83465e78783e127453925715d61b78fa715a5b27d35555b5c142710101b2a
```

The independent audit rederived every proof identity, extended the exact
coefficient-convolution controls through `n=60`, and replayed normal,
optimized, and stored transcripts.  After its requested direct central
leading-coefficient quotient check was added, no defect remained.

**QED.**
