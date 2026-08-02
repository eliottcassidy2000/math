---
id: THM-3148
title: "Fixed-offset Frobenius endpoint-resultant classification"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.
  The factorial functional modulo p is a Frobenius endpoint projector.  At
  every fixed offset d=p+s it identifies the complete height-zero residual
  of A_(p+a) with a fixed degree-a polynomial.  For the resonant pair this
  descends the large window r=p+s-2 to the base window r=s-2, and its fixed
  resultant exactly classifies the exceptional primes on the unit-root face.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-3131-prime-resonance-newton-slope-separation
  - THM-3138-left-factorial-resonance-newton-slope-separation
  - THM-3143-two-step-prime-resonance-euclidean-newton-separation
  - THM-3146-three-step-prime-resonance-full-euclidean-newton-separation
script: 04-computation/factorial_fixed_offset_frobenius_endpoint_resultant_thm3148.py
output: 05-knowledge/results/factorial_fixed_offset_frobenius_endpoint_resultant_thm3148.out
script_sha256: 014cc689c1f5794008e37232e000405f463a4f6ee24ae4bc998c437e0919f84c
output_sha256: 86b6a987641ca64080ee5dcbfe7f71ce05a3f095f009e8925976770e6df4918b
hash_basis: LF-normalized bytes
---

# THM-3148 -- fixed-offset Frobenius endpoint-resultant classification

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

The prime-resonance theorems THM-3131, THM-3138, THM-3143, and THM-3146
computed their low Newton faces one offset at a time.  The common mechanism is
an exact Frobenius descent: modulo `p`, the factorial functional forgets every
nonconstant term of a `p`-th power.  Thus a window of size about `p` has the
same height-zero residual as a fixed small window.

## 1. The factorial Frobenius projector

Let

```text
L(t^k)=k!,
```

extended coefficientwise from `Z[v][t]` to `Z[v]`.  Let `p` be prime and let
`P,G in Z[v][t]`, with the constant term `P(0)=d` an integer independent of
`v`.  Then

```text
L(P(t)^p G(t)) = d L(G(t))                         (mod p). (1)
```

There is no degree hypothesis on `G` in `(1)`.

### Proof

In characteristic `p`,

```text
P(t)^p=d^p+sum_(i>=1) P_i(v)^p t^(ip).                     (2)
```

Every monomial contributed by the displayed sum has `t`-degree at least
`p`, and its factorial is therefore zero modulo `p`.  Only the constant term
survives.  Fermat gives `d^p=d (mod p)`, proving `(1)`.

## 2. Fixed-offset residual polynomial

For integers `s>=1`, `a>=0`, put

```text
A_(p+a)^[p+s](v)=L((p+s-t+v t^2)^(p+a)),
F_(a,s)(v)=L((s-t+v t^2)^a).                               (3)
```

Applying `(1)` with `P=G=p+s-t+v t^2` gives the coefficientwise polynomial
congruence

```text
A_(p+a)^[p+s](v) = s F_(a,s)(v)                    (mod p). (4)
```

Explicitly,

```text
[v^j]F_(a,s)
=binom(a,j) sum_(ell=0)^(a-j)
 binom(a-j,ell)s^(a-j-ell)(-1)^ell(2j+ell)!,               (5)
```

and the leading coefficient is `(2a)!`.  Consequently, if

```text
p>2a,                         p does not divide s,          (6)
```

then `sF_(a,s) mod p` has degree exactly `a` and is the complete
untruncated height-zero residual of the large polynomial in `(3)`: every
coefficient above degree `a` has positive `p`-valuation, while the nonzero
coefficients of `sF_(a,s)` have valuation zero.

The qualifier *height-zero residual* is load-bearing.  Equation `(4)` says
nothing about the positive-height faces of the Newton polygon.

## 3. Resonant windows descend to smaller resonant windows

Now take `s>=2` and the quadratic-moment resonance

```text
d=p+s,                  r=d-2=p+s-2.                       (7)
```

The pair whose coprimality is needed is

```text
A_r^[d],                         A_(r+1)^[d].               (8)
```

Set `r_0=s-2`.  Under

```text
p>2(s-1),                                                     (9)
```

their complete height-zero residual pair is

```text
sF_(r_0,s),                       sF_(r_0+1,s).              (10)
```

Define the primitive fixed resultant

```text
delta_s=Res_v(F_(s-2,s),F_(s-1,s)).                         (11)
```

Since both degrees in `(10)` are preserved, the two residuals have a common
affine root over the algebraic closure of `F_p` if and only if

```text
p divides delta_s.                                          (12)
```

That common root, when it exists, is automatically nonzero.  Indeed, put

```text
I_a(s)=F_(a,s)(0)=L((s-t)^a).
```

One integration by parts gives

```text
I_(a+1)(s)=s^(a+1)-(a+1)I_a(s).                             (13)
```

If two consecutive constants vanished modulo `p`, then `(13)` would force
`p|s`, contrary to `(9)`.  Thus `(12)` classifies genuine unit-root,
height-zero coincidences rather than a spurious root at the origin.

The descent is literal.  If `P_(n,r_0)(u)` denotes the normalized resonant
polynomial from THM-3124, then

```text
F_(n,s)(v)=s^n P_(n,r_0)(v/s).                              (14)
```

Hence

```text
delta_s != 0
iff the base resonant pair at window r_0=s-2 is coprime over C. (15)
```

THM-3124 proves the right side for `0<=r_0<=200`.  Therefore

```text
delta_s != 0 for every 2<=s<=202,                            (16)
```

and each such fixed offset has only finitely many exceptional primes on its
height-zero face.  No assertion that `delta_s!=0` for all `s` is smuggled
into the theorem: by `(15)`, that assertion is exactly the still-open
unbounded quadratic moment problem in smaller coordinates.

When `p|delta_s`, the ordinary subresultant/Euclidean chain over `F_p`
recovers the exact residual gcd degree and polynomial.  This is a finite
calculation depending only on `s`, not on the large degree `p+s`.

## 4. The first exceptional spectra are genuinely nonlocal

The exact primitive resultants through `s=6` are

```text
s  base r_0   factorization of delta_s
2      0      1
3      1      2^2 * 29
4      2      2^11 * 3^2 * 7 * 4547
5      3      2^16 * 3^7 * 5^2 * 11^2 * 20747 * 249721
6      4      2^47 * 3^15 * 5^4 * 7 * 139 * 3767
              * 12041 * 807241.                            (17)
```

The resultant of the actual pair `(10)` is

```text
Res(sF_(s-2,s),sF_(s-1,s))=s^(2s-3) delta_s.               (18)
```

The outer power of `s` creates no extra exception under `(9)`.

At `s=3`,

```text
F_(1,3)=2v+2,                  F_(2,3)=24v^2+5,             (19)
```

and the unique eligible exception is `p=29`, with common residual factor
`v+1`.  This is exactly the endpoint residue `87=3*29` encountered in
THM-3146.

At `s=4`,

```text
F_(2,4)=24v^2+4v+10,
F_(3,4)=720v^3-72v^2+24v+34,                               (20)

F_(3,4)=(30v-8)F_(2,4)+(456-976v).                         (21)
```

The eligible exceptional prime `4547` has exact gcd

```text
v+1304,                         root v=3243 in F_4547.      (22)
```

Thus exceptional primes are not confined to the old finite window census.
The prime `61`, although it divides the leading coefficient `976` in the
next Euclidean quotient, does **not** divide `delta_4`; its modular gcd is
`1`.  This separates quotient-denominator walls from actual residual-root
exceptions.

For every eligible exceptional prime in `(17)`, the companion computes a
linear, nonzero residual gcd.  Higher offsets may have larger gcd degree;
the theorem claims the subresultant classification, not universal linearity.

## 5. What the descent preserves and destroys

The map in `(4)` preserves:

- the full reduction of the large moment polynomial modulo `p`;
- every height-zero coefficient and its affine residual roots;
- the exact fixed-offset resultant and subresultant obstruction.

It destroys:

- all coefficients of positive `p`-valuation;
- every positive Newton slope and its residual polynomial;
- a lift of a modular residual root to a common algebraic root over `Q`.

Therefore `p not dividing delta_s` excludes a common *unit* root of `(8)`
over `Q_p`, but does not by itself prove that `(8)` is coprime.  A complete
argument still has to separate every common positive slope, as THM-3146 did
at `s=3`.  Conversely, `p|delta_s` is an obstruction to this particular
low-face certificate, not a counterexample to the Gaussian Moment
Conjecture.

No claim is made here for arbitrary radial polynomial coefficients, other
supports, all of NC2, or the Gaussian Moment Conjecture in two dimensions.

## 6. Exact companion

Run

```text
python 04-computation/factorial_fixed_offset_frobenius_endpoint_resultant_thm3148.py
python -O 04-computation/factorial_fixed_offset_frobenius_endpoint_resultant_thm3148.py
```

and compare byte-for-byte with

```text
05-knowledge/results/factorial_fixed_offset_frobenius_endpoint_resultant_thm3148.out.
```

The companion uses only integer and rational arithmetic.  It checks `1,400`
coefficientwise Frobenius congruences, independent bivariate expansions, the
scaling and zero-root recurrences, Bareiss resultants, exact factorizations,
all eligible modular gcds through `s=6`, and the hostile `p=61` denominator
wall.

**QED (candidate pending independent hostile audit).**
