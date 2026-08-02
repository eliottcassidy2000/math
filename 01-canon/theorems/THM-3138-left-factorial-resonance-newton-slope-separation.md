---
id: THM-3138
title: "Left-factorial resonance Newton-slope separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every odd
  prime p, no exact quadratic has three zero factorial moments beginning at
  r=p-1.  At the composite resonance d=r+2=p+1, the first polynomial has an
  unconditional midpoint unit: all slopes to its left are nonpositive and
  its final slope is 2/(p-1).  The second polynomial has sole slope 2/p, so
  they are coprime.  The left-factorial residue !p controls whether the left
  part is horizontal, but Kurepa's conjecture is not needed for separation.
audit: >
  An independent hostile audit rederived the multinomial coefficients,
  derangement/left-factorial congruence, unconditional midpoint unit, the
  lower-hull vertex and convex-slope argument even under a hypothetically
  raised or missing constant endpoint, the complete A_p polygon, p=3, and
  common-factor base change.  The companion scans all 1,228 odd primes below
  10,000, reconstructs selected full integer polygons, exhausts 90,484
  synthetic raised-endpoint valuation words, and replays normal, optimized,
  and stored transcripts exactly.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-3131-prime-resonance-newton-slope-separation
  - THM-3142-prime-power-predecessor-newton-separation-and-composite-window-census
script: 04-computation/factorial_left_factorial_resonance_newton_slopes_thm3138.py
output: 05-knowledge/results/factorial_left_factorial_resonance_newton_slopes_thm3138.out
script_sha256: 227eab165bca3f8c4bf68f698188cde72a44cca0d557d96f7364c429d889e1f6
output_sha256: f305691769975b2b479edf2f0b4b0a1803545d3af764c9ae443789937450f03f
hash_basis: LF-normalized bytes
---

# THM-3138 -- left-factorial resonance Newton-slope separation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(t^k)=k!,                         q(t)=a+bt+ct^2,             (1)
```

with `abc!=0`.  For every odd prime `p`, put

```text
r=p-1,                             d=r+2=p+1.                 (2)
```

Then

```text
L(q^r), L(q^(r+1)), L(q^(r+2))                                (3)
```

cannot all vanish.  This is a uniform all-prime theorem, not a finite
census and not a consequence conditional on Kurepa's conjecture.

## 1. Resonant integral pair

THM-3124 proves that `(3)` forces

```text
b/a=-1/(r+2)=-1/d.                                             (4)
```

After division by `a`, write `u=c/a`.  It remains to exclude a common root of

```text
P_(p-1)(u)=L((1-t/d+u t^2)^(p-1)),
P_p(u)    =L((1-t/d+u t^2)^p).                                (5)
```

The invertible change `v=d u`, followed by multiplication by `d^n`, gives
the integral polynomials

```text
A_n(v)=d^n P_n(v/d)=L((d-t+v t^2)^n) in Z[v].                 (6)
```

Writing `A_n(v)=sum_(j=0)^n a_(n,j)v^j`, direct expansion gives

```text
a_(n,j)=binom(n,j) sum_(ell=0)^(n-j)
  binom(n-j,ell)d^(n-j-ell)(-1)^ell(2j+ell)!.                 (7)
```

Since `d=p+1==1 (mod p)`, the neighboring prime `p=d-1` sees a coefficient
geometry unavailable at primes dividing the composite resonance `d`.

## 2. Midpoint rigidity and the first polygon

Put `m=(p-1)/2`.  At the midpoint `j=m`, every summand of `(7)` with
`ell>=1` contains `(p-1+ell)!` and hence vanishes modulo `p`.  Wilson's
theorem and `binom(p-1,m)==(-1)^m (mod p)` give

```text
a_(p-1,m) == binom(p-1,m)(p-1)!
            == (-1)^(m+1) !=0                         (mod p). (8)
```

Thus `(m,0)` is a Newton-polygon point.  Every earlier nonzero coefficient
has nonnegative integral valuation.  Moreover `(m,0)` cannot be bypassed by
a chord to a later point: every later coefficient has valuation at least one,
so such a chord has strictly positive height at `m`.  Hence `(m,0)` is a
lower-hull vertex.  The convex slopes of the hull to its left are all
nonpositive, because their right endpoint has the minimum possible height
zero.

For `j>m`, every factorial in `(7)` is divisible by `p`, while

```text
a_(p-1,p-1)=(2p-2)!,             v_p((2p-2)!)=1.              (9)
```

Every intervening valuation is at least one and lies strictly above the
segment from `(m,0)` to `(p-1,1)`.  Consequently

```text
slopes(NP_p(A_(p-1)))
    are <=0 to the left of m, followed by exactly 2/(p-1).    (10)
```

Zero coefficients are simply omitted from the polygon and do not change
this argument.

## 3. The left-factorial endpoint is informative, not load-bearing

The constant coefficient reduces to a classical residue:

```text
a_(p-1,0) == L((1-t)^(p-1)) == D_(p-1)              (mod p), (11)
```

where `D_n` is the derangement number.  Indeed

```text
D_(p-1)=sum_(ell=0)^(p-1)
  binom(p-1,ell)(-1)^ell ell!
       == sum_(ell=0)^(p-1) ell! = !p                 (mod p), (12)
```

because `binom(p-1,ell)==(-1)^ell (mod p)`.  If Kurepa's conjectural
nonvanishing `!p!=0 (mod p)` holds, the left portion in `(10)` is exactly the
horizontal segment `(0,0)--(m,0)`, so the full slope set is
`{0,2/(p-1)}`.

But this endpoint unit is not needed.  At a hypothetical Kurepa prime the
left endpoint could rise, producing one or more negative slopes into the
unconditional midpoint unit.  Those slopes remain separated from the
positive slope of the second polynomial below.  This is why the initially
tempting conditional theorem strengthens to an unconditional one.

In fact `a_(p-1,0)>0` over the integers, since it is the integral of the
nonzero nonnegative function `exp(-t)(d-t)^(p-1)` on `[0,infinity)`.  Even
without using that positivity, an exact zero constant would contribute only
a factor `v`; it cannot be common because the next polynomial has unit
constant coefficient.

## 4. The second polygon

For `A_p`, the constant coefficient is a unit.  In `(7)` the `ell=0` term is
`d^p==1 (mod p)`, while every other term contains either `binom(p,ell)` or
`p!`.

For every `0<j<p`, the outer coefficient `binom(p,j)` supplies one factor of
`p`.  If `j>=m+1`, every factorial `(2j+ell)!` supplies a second factor.
Finally

```text
a_(p,p)=(2p)!,                         v_p((2p)!)=2.          (13)
```

Thus for `1<=j<=m` the valuation is at least one, strictly above the line of
height `2j/p`; for `m+1<=j<p` it is at least two, again strictly above that
line.  Therefore

```text
NP_p(A_p) : (0,0)--(p,2),             sole slope 2/p.        (14)
```

This positive slope is unequal to `2/(p-1)` and cannot equal any of the
nonpositive prefix slopes in `(10)`.

## 5. Common-factor exclusion

The Newton product lemma says that over a valued field the lower slopes of a
product, counted with horizontal multiplicity, are the multiset union of the
factor slopes.  Equivalently, after splitting over an algebraic closure of
`Q_p`, a nonzero root `alpha` contributes slope `-v_p(alpha)`.

If `A_(p-1)` and `A_p` had a nonconstant common factor over `Q`, base change
to `Q_p` would preserve a nonconstant common factor.  Any nonzero common root
would contribute a common finite slope, forbidden by `(10)` and `(14)`; zero
is not a root of `A_p`.  The pair is therefore coprime over `Q`, hence has no
common complex root.  By the invertible scaling `(6)`, neither do the pair in
`(5)`, contradicting `(3)`.  QED.

## 6. Exact controls and arithmetic boundary

The companion scans all `1,228` odd primes below `10,000`.  It finds no
vanishing left-factorial residue and checks independently, by the derangement
recurrence, every congruence `(12)`.  The semantic residue digest is

```text
8f44042786c94d77af31433cac2739dc88967074eb42250e3f3d791bff03ed5a. (15)
```

For `p=3,5,7,11,101,211` it reconstructs the complete integer polynomials
from `(7)`, computes their lower hulls, and gets `(10)` with a horizontal
prefix and `(14)`.  A repeated-product bivariate constructor gives an
independent coefficient path for `p<=11`.  Finally, `90,484` synthetic
valuation words exhaust prefix heights in `{0,1,2,3}` and later heights in
`{1,2,3}` through midpoint `m=5`; every raised-endpoint hull has only
nonpositive prefix slopes and the final slope `1/m`, while avoiding `2/p`.
This finite synthetic control illustrates the general convex proof rather
than replacing it.

If Kurepa nonvanishing failed, only the exact horizontal-prefix description
would fail.  The midpoint, disjointness, coprimality, and factorial-moment
conclusion survive unchanged.  Thus Kurepa is an arithmetic classifier of
the polygon chamber, not a proof dependency.

Run

```text
python3 04-computation/factorial_left_factorial_resonance_newton_slopes_thm3138.py
python3 -O 04-computation/factorial_left_factorial_resonance_newton_slopes_thm3138.py
```

and compare byte-for-byte with the declared output.

## 7. Connection contract and scope

The source object is THM-3124's resonant pair `(5)`.  The map is the integral
rescaling `(6)` followed by the `p`-adic Newton polygon at the neighboring
prime `p=d-1`.  It preserves common algebraic factors and forgets root phase
and residue data.  The preserved predicate is coprimality; the load-bearing
sidecars are the unconditional midpoint unit `(8)` and leading valuations
`(9),(13)`.  The left-factorial endpoint `(12)` refines the chamber but is not
load-bearing.

Combining this theorem with THM-3131 and THM-3124, any still-open bad exact
quadratic window must satisfy

```text
r>=201,       d=r+2 composite,       and d-1 composite.       (16)
```

The surprising move is to inspect a prime adjacent to, rather than dividing,
the composite resonance.  This is an exact `{0,1,2}` / `SFC(1)`
factorial-moment theorem.  It does not settle arbitrary-support `SFC(3)`,
`GMC(2)`, `NC(2)`, or `LRC(14)`.

THM-3142 subsequently extends the exclusion conclusion to every prime-power
predecessor `d-1=p^k` and closes the finite exact residual through `d=1000`;
the left-factorial chamber analysis here remains the finer `k=1`, odd-prime
statement.

For context on the still-open left-factorial problem, see
Andrejic--Bostan--Tatarevic, *Improved algorithms for left factorial
residues*, arXiv:1904.09196 (2019).  That external computation is not a
dependency of any claim here.

## 8. Independent hostile audit

The audit independently rebuilt `(7)`, `(8)`, and `(12)`, checked every
valuation threshold and the `p=3` boundary, and first accepted the stronger
conditional exact-polygon statement.  It then attacked the Kurepa unit and
confirmed the unconditional repair: `(m,0)` cannot be bypassed, convexity
forces all prefix slopes nonpositive, and a missing zero root cannot be
shared with the unit-constant second polynomial.  It also rechecked the
Newton product/base-change implication and exact transcript hashes.  No
conditional arithmetic claim is used in the proof.

**End of proof.**
