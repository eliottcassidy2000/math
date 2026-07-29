---
id: THM-2964
title: "General pure factorial-moment resonance ladder"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For every order k>=2 and width M, the pure denominator-cleared
  factorial-moment coefficient has an interior negative-integer root -r
  exactly when M=(k+1)d+1, r=kd+1, and (k-1)d is even.  Every such root is
  simple.  This proves the quadratic, cubic, and quartic resonance ladders
  seen in THM-2960 for all widths and all normalized four-slot interiors.
  It classifies one coefficient only: no global cofactor divisibility,
  matrix-bar exclusion, or GMC nonvanishing is claimed.
source: codex-gmc-general-pure-resonance-2026-07-29
depends_on: []
related:
  - THM-2925-general-width-terminal-pole-cancellation
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
  - THM-2960-local-smith-jet-fitting-barcode-and-negative-depth-chamber-atlas
script: 04-computation/gmc_general_pure_factorial_moment_resonance_ladder_thm2964.py
output: 05-knowledge/results/gmc_general_pure_factorial_moment_resonance_ladder_thm2964.out
script_sha256: 488783af2b96c2c21f1751fa09e457ac94682d1c9ce4f8433d4fd4a4483de9a7
output_sha256: d2cccee97efa77e4d8b903ac6545447626cd8538291ac13b1369704843495044
hash_basis: LF-normalized bytes
---

# THM-2964 -- general pure factorial-moment resonance ladder

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

The local Smith atlas in THM-2960 found three parallel wall laws:

```text
q200: M=3d+1, r=2d+1, d even;
c300: M=4d+1, r=3d+1;
f400: M=5d+1, r=4d+1, d even.                         (1)
```

The quartic members `(M,r)=(11,9),(21,17)` initially looked like unrelated
matrix anomalies.  They are instead the first two instances of one
all-order analytic ladder.

## 1. Universal pure coefficient

Let \(k\ge2\), \(M\ge3\), and use the rising-factorial notation

```text
(x)_a=x(x+1)...(x+a-1),        (x)_0=1.
```

After eliminating the coefficient at exponent \(M\), the pure coefficient
of the order-\(k\) factorial-moment form is, up to the standard positive
normalization,

```text
          k
S_(k,M)(n)=sum (-1)^j binom(k,j)
          j=0

                         (kn+1)_(jM)
                    * -----------------.               (2)
                          (n+1)_M^j
```

Indeed, the \(j\)-th term chooses the eliminated exponent \(M\) in exactly
\(j\) of the \(k\) tensor positions.  Formula `(2)` therefore does not
depend on the other interior exponents in a normalized support
`(0,a,b,M)`.

Put

```text
D_(k,M)(n)
 =prod_(s=1)^(M-1)(n+s)^(k-1) * (n+M)^(k-2),

P_(k,M)(n)=D_(k,M)(n) S_(k,M)(n).                     (3)
```

The apparent poles of `(2)` lie at the negative integer walls.  At
\(-s\), \(s<M\), every summand has pole order at most \(k-1\).  At the
terminal wall \(-M\), the order-\(k-1\) parts of the last two summands
cancel, leaving order at most \(k-2\).  Thus `(3)` is a polynomial.  An
equivalent exact integer construction, used by the companion, is

```text
P_(k,M)
 =
 sum_j (-1)^j binom(k,j)(kn+1)_(jM)(n+1)_M^(k-j)
 ----------------------------------------------------. (4)
                   (n+1)_M(n+M)
```

The division in `(4)` is exact.

## 2. The resonance theorem

For every integer wall

```text
1<=r<M,
```

the following are equivalent:

1. \(P_{k,M}(-r)=0\);
2. there is an integer \(d\ge1\) such that

```text
M=(k+1)d+1,        r=kd+1,        (k-1)d is even.      (5)
```

Every root in `(5)` has multiplicity exactly one.

For odd \(k\), the parity condition is automatic.  For even \(k\), it says
that \(d\) is even.  Hence `(5)` specializes exactly to the three laws in
`(1)`.

## 3. Laurent proof of the zero classification

Fix \(1\le r<M\), set

```text
z=n+r,                    d=M-r,
```

and write

```text
(n+1)_M=z A(1+O(z)),
A=(-1)^(r-1)(r-1)!d!.                                  (6)
```

Denote the \(j\)-th summand in `(2)` by \(T_j\).  Its numerator contains
the factor \(z\) exactly when

```text
kr<=jM.                                                  (7)
```

Consequently, if \(kr\le(k-1)M\), only \(T_k\) reaches pole order
\(k-1\).  Its leading coefficient is nonzero, so \(P_{k,M}(-r)\ne0\).

Suppose \(kr>(k-1)M\).  Exactly \(T_{k-1}\) and \(T_k\) now reach pole
order \(k-1\).  The magnitudes of their leading Laurent coefficients are
equal precisely when

```text
(r-1)!d!
---------- = (kd)!,                                     (8)
(M-kd-1)!
```

or, equivalently,

```text
binom(r-1,(k-1)d)=binom(kd,(k-1)d).                     (9)
```

Because \(d\ge1\), the lower binomial argument in `(9)` is positive.
For fixed positive lower argument, the binomial coefficient is strictly
increasing in its upper argument.  Therefore `(9)` holds exactly when

```text
r-1=kd,
```

which is equivalent to the first two equations in `(5)`.

At this magnitude equality, the leading coefficient from \(T_k\) is
negative.  The leading coefficient from \(T_{k-1}\) is positive exactly
when

```text
(k-1)d=0 mod 2.                                         (10)
```

Thus cancellation occurs exactly under all three conditions in `(5)`.
This proves the zero classification.

## 4. Why every resonance is simple

Assume `(5)`, and let

```text
H_a=sum_(i=1)^a 1/i.
```

Relative to the common positive order-\(-(k-1)\) coefficient of
\(T_{k-1}\), the next Laurent coefficient of
\(T_{k-2}+T_{k-1}+T_k\) is

```text
(k+1)(H_d-H_(kd))

 +(-1)^d (k-1)/2 * (kd)!(d!)^2/((k+2)d+1)!.            (11)
```

All earlier summands have lower pole order and do not enter `(11)`.
If \(d\) is odd, both terms in `(11)` are negative.  If \(d\) is even,
the second term is positive but has magnitude less than \((k-1)/2\),
whereas

```text
(k+1)(H_(kd)-H_d)
 >= (k+1)(k-1)/k
 >  (k-1)/2.                                            (12)
```

Therefore `(11)` never vanishes.  Multiplication by the factor
\(z^{k-1}\) in `(3)` leaves a nonzero linear term, so every resonance
root is simple.

## 5. Exact regression and controls

The standalone companion reconstructs `(4)` without importing any earlier
GMC script.  It checks every one of the

```text
3,564 walls with 2<=k<=7 and 3<=M<=35.                 (13)
```

All \(27\) roots agree with `(5)` and have valuation one.  The frozen
semantic digest is

```text
8bdc888577f7903f86f82b8fde177da3c02e22d405fce693ffb909b766cc011a.
```

Normal and optimized executions byte-match the stored transcript.

The finite computation is a hostile implementation check, not the proof:
Sections 3--4 establish the statement for every \(k,M\).

## 6. Scope and first failed implication

The exact connection is:

```text
source:      one pure factorial-moment coefficient;
target:      its negative-integer wall divisor;
preserved:   exact width, order, root, parity, and multiplicity;
lost:        the rest of the Macaulay matrix and every other coefficient;
sidecar:     local Smith jets or a global resultant/cofactor argument.
```

A root of \(P_{k,M}\) can account for a local Smith bar, as the quartic
roots do at widths `11` and `21`.  It does **not** imply that the whole
polynomial \(P_{k,M}\) divides a chosen cofactor, and absence of a pure
root does not exclude a matrix-level bar.  THM-2960's `(M,r)=(12,5)`
wall, where `q200,c300,f400` are all units, is the sharp current hostile.
No cofactor nonvanishing, arbitrary-width GMC closure, or fixed-prime rank
claim follows from this theorem alone.

**QED (candidate; independent hostile audit pending).**
