---
id: THM-3180
title: "Six-exit factorial Newton--Euclidean closure through r=2198"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For the resonant exact-support quadratic factorial pair, the THM-3152
  first Euclidean flag and fixed prime bank through 47 close every six-exit
  residual with 2001<=d<=2200.  There are exactly 79 such rows; two exact
  implementations agree on every progressive degree-set trace, and no prime
  beyond 29 is actually needed.  Together with THM-3161 and the inherited
  prime, prime-power-predecessor, two-step-prime, three-step-prime,
  four-step-prime, and five-step-prime exits, every exact three-term
  quadratic factorial window beginning at 1<=r<=2198 contains a nonzero
  moment.  This is bounded and does not prove an all-height result or FC(3).
audit: >
  The primary companion imports THM-3152's pinned exact recurrence,
  zero-root-aware denominator-capacity barcode, Fraction lower hull, and
  first full Euclidean remainder.  A self-contained independent companion
  reconstructs primality and prime-power filters, every polynomial row and
  remainder, coefficient valuations, and determinant lower hulls.  Both
  replay all 79 six-exit residuals in four deterministic chunks.  Their four
  chunk trace digests, global semantic trace digest, killer-prime histogram,
  and d=2009 hostile trace agree exactly.  Planted factors v+1 and v survive.
  Both scripts use explicit RuntimeError checks, contain no assert statement,
  and have byte-identical normal and optimized transcripts.
source: root/frontier-synthesis/factorial-six-exit/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3131-prime-resonance-newton-slope-separation
  - THM-3142-prime-power-predecessor-newton-separation-and-composite-window-census
  - THM-3143-two-step-prime-resonance-euclidean-newton-separation
  - THM-3146-three-step-prime-resonance-full-euclidean-newton-separation
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
  - THM-3153-four-step-prime-resonance-second-euclidean-newton-separation
  - THM-3161-factorial-newton-euclidean-closure-through-r1998
  - THM-3170-five-step-prime-resonance-euclidean-newton-holotopy
related:
  - THM-3176-six-step-prime-resonance-third-euclidean-newton-separation
script: 04-computation/factorial_six_exit_first_flag_2200_thm3180.py
output: 05-knowledge/results/factorial_six_exit_first_flag_2200_thm3180.out
script_sha256: 50b4691c72d5d0b942eafdb922273ee749127f9725bd3728a9b422bfc494e69e
output_sha256: a7f5d4514b7ecab2da1fee1d892ca2d46664e0c41f59d49563dcd930c4109634
independent_script: 04-computation/factorial_six_exit_first_flag_2200_independent_audit_thm3180.py
independent_output: 05-knowledge/results/factorial_six_exit_first_flag_2200_independent_audit_thm3180.out
independent_script_sha256: 5d00729ae39cdf1d085eca73d4612d4960a0788dd746f945ab2e8936ebbd466d
independent_output_sha256: c244bea4a36e0663f421892b45eed53ecf6f7a3a421003b112542700cadd3ae0
hash_basis: LF-normalized bytes
---

# THM-3180 -- six-exit factorial Newton--Euclidean closure through r=2198

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**

## 1. Statement and scope

Let

```text
L(t^m)=m!,                    q(t)=a+bt+ct^2,               (1)
```

with `abc!=0`.  For every integer `1<=r<=2198`, the three moments

```text
L(q^r),                 L(q^(r+1)),                 L(q^(r+2)) (2)
```

cannot all vanish.

This is a bounded exact-support `{0,1,2}` theorem.  It does not cover a
missing coefficient, translated or arbitrary support, all of `SFC(1)`,
`SFC(3)`, or the three-variable Factorial Conjecture `FC(3)`.

## 2. Inherited range, exits, and exact residual universe

THM-3161 already proves `(2)` for `1<=r<=1998`, equivalently `3<=d<=2000`
under

```text
d=r+2.                                                        (3)
```

Fix `1999<=r<=2198`, so `2001<=d<=2200`.  THM-3124 forces a hypothetical
bad window to the resonance `b/a=-1/d` and reduces it to a common root of

```text
P=A_(d-2)^(d),          Q=A_(d-1)^(d),
A_n^(d)(v)=L((d-t+vt^2)^n).                                 (4)
```

The uniform inherited exits in the new range are

```text
d prime:                 THM-3131,
d-1 a prime power:       THM-3142,
d-2 prime:               THM-3143,
d-3 prime:               THM-3146,
d-4 prime:               THM-3153,
d-5 prime:               THM-3170.                          (5)
```

The prime hypotheses in the offset theorems are odd; every prime appearing
in the last four lines of `(5)` is greater than two in the present range.
Call `d` a **six-exit residual** when none of `(5)` applies.  Exact
primality and factorization give the following progressive survivor counts:

```text
universe          d     d-1   d-2   d-3   d-4   d-5
3..1000           831   642   513   390   301   211
1001..2000        865   725   617   511   426   341
3..2000          1696  1367  1130   901   727   552
2001..2200        176   149   130   111    95    79.        (6)
```

Only the last row is a new computational proof obligation.  Its exact
ordered residual-list digest is

```text
8d93a3f8b5fb2877e35872d96b2dbce0d3a219265e0718471e5913b54236fad5.
                                                                    (7)
```

The deterministic work chunks `[2001,2050]`, `[2051,2100]`,
`[2101,2150]`, `[2151,2200]` contain respectively `15,16,18,30`
six-exit residuals.

## 3. First Euclidean flag and finite closure

Put `n=d-2`.  THM-3152 proves that the integral first full Euclidean row

```text
R_1=(2n-1)(Q-2(n+1)(2n+1)vP)+2d(n+1)P                     (8)
```

has degree at most `n-1` and

```text
gcd_Q(P,Q)=gcd_Q(P,R_1).                                    (9)
```

For every six-exit residual in `(6)`, the exact degrees are

```text
(deg P,deg Q,deg R_1)=(d-2,d-1,d-3).                       (10)
```

At every prime in the fixed bank

```text
S={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47},              (11)
```

the companions construct THM-3152's zero-root-aware necessary
common-factor degree barcode `D_p(P,Q,R_1)`.  For all 79 rows,

```text
intersection_(p in S) D_p(P,Q,R_1) intersection Z_(>0)
 =empty.                                                     (12)
```

By the proved degree-barcode lemma, `(12)` excludes every nonconstant
rational factor common to `P,Q,R_1`, hence every nonconstant common factor
of `P,Q`.  Since rational polynomials with a common complex root have a
nonconstant rational gcd, no six-exit residual supports a bad window.
Together with `(5)` this closes `2001<=d<=2200`; THM-3161 supplies
`3<=d<=2000`.  Thus `(2)` holds through `r=d-2=2198`.  QED.

The full semantic trace records the surviving positive degrees after every
successive prime, not merely final emptiness.  Both implementations give

```text
[2001,2050] f1278e4b6f1855e75edcff6e02f0d6c2b00a1ce82a2abf15492d7c3d2a0a6a0d
[2051,2100] 5d33c39afdfa96a0e8d0ecc3c34a4ebbcc0ea644b3d22a2fc8927689e87a7281
[2101,2150] 7884bb2d58460b834be608fa68ea6acfadc76c489e9806f4fbf0ce4aaeb73341
[2151,2200] 30617545e693dc35b7226deb83d6e5e93d81bce89d578f3360d642f5eb7eb2c6
global       650960f5ec5129fbb170360331dec10843172f81584c954c55216c5f9c166052.
                                                                    (13)
```

The first killing-prime histogram is

```text
2:8, 3:5, 5:20, 7:12, 11:10, 13:1, 17:11, 19:5, 23:5, 29:2. (14)
```

Thus no prime beyond 29 is actually needed, although the audited fixed bank
is retained through 47.

## 4. Independent audit and controls

The primary artifact imports the exact THM-3152 companion with pinned hash

```text
f804d3996abe4df981dbf7db877af4aeca9218df64b0ac382af876a3cdca15a0.
                                                                    (15)
```

The independent artifact imports none of that mathematical engine.  It
separately implements trial-division primality and factorization, the
integer moment recurrence, the full first remainder, `p`-adic valuations,
and a determinant cross-product lower hull.  Both replay the entire promoted
range in the same deterministic chunks and agree on `(7),(13),(14)`.

For the hostile row `d=2009`, `p=2` retains the 15 degrees

```text
128,256,...,1920;
```

`p=3` retains the same set and `p=5` kills it.  Both implementations
reproduce that full trace.  Planted common factors `v+1` and `v` retain
degree one, guarding against an overstrong observer.  All load-bearing
checks use explicit `require` calls raising `RuntimeError`; neither script
contains an `assert` statement.  Normal and `python -O` executions compare
byte-for-byte with the stored transcripts.

## 5. THM-3176 and the exact first unaudited case

THM-3176 proves the additional uniform exit

```text
d-6 prime.                                                   (16)
```

It is related, not a dependency of this theorem: the computation in
Section 3 closes the larger 79-row six-exit universe, including the 13 rows
in that universe for which `d-6` is prime.  If `(16)` is applied first, only
66 seven-exit residuals remain in `2001<=d<=2200`.

The first height outside the proved interval is `r=2199`, with `d=2201`.
It is also the first post-`2200` seven-exit residual, not merely the next
integer:

```text
d   =2201=31*71,
d-1 =2200=2^3*5^2*11                 (not a prime power),
d-2 =2199=3*733,
d-3 =2198=2*7*157,
d-4 =2197=13^3,
d-5 =2196=2^2*3^2*61,
d-6 =2195=5*439.                                           (17)
```

Thus none of the seven uniform exits applies at `r=2199`, and the finite
first-flag computation stops immediately before it.  Equation `(17)` is an
exact boundary invoice, not a common-root witness or counterexample.

## 6. Reproduction

Run from the repository root:

```text
python3 04-computation/factorial_six_exit_first_flag_2200_thm3180.py
python3 -O 04-computation/factorial_six_exit_first_flag_2200_thm3180.py
python3 04-computation/factorial_six_exit_first_flag_2200_independent_audit_thm3180.py
python3 -O 04-computation/factorial_six_exit_first_flag_2200_independent_audit_thm3180.py
```

Compare LF-normalized bytes with the two declared stored outputs.  The
primary companion resolves its pinned THM-3152 dependency relative to its
own file, so it is independent of the caller's current directory.
