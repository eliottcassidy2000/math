---
id: THM-3894
title: "Cusp residual all-degree gauge--Kummer parity filtration"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  Iterating the THM-3886 leading symbol gives
  an all-degree parity passport on the THM-3884 equality seam.  Every
  even-degree survivor aligns through the middle gauge depth and then must
  enter one Kummer symbol q_0=xu^2.  Every odd-degree survivor aligns one
  layer farther and must have square q_0.  This is associated-graded only:
  no gauge-peeling invariance or existence is claimed.  An independent
  vertical-gauge corollary proves that a Kummer passport cannot be completed
  with all remaining data in k[x].
source: jc_quartic_c3_construct / post-THM-3886 iterated leading-symbol filtration, 2026-08-23
audit: >
  SELF-AUDITED proof candidate.  The exact companion verifies the universal
  component cancellation, all degree comparisons and parity crossings, the
  Kummer and square normal forms, and the vertical-gauge quartic
  discriminant factorization.  It freezes the L-adic reduction and Bezout
  coprimality used in the vertical no-go.  Normal and optimized runs must
  byte-match the frozen transcript.  Independent audit must recheck the
  induction indices, exhaustion of competing residual terms, the UFD
  square implications, and the even-in-y reduction in the vertical slice.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3884-cusp-residual-total-degree-leading-gauge-filtration
  - THM-3886-cusp-residual-equality-seam-second-layer-trichotomy
related:
  - THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure
script: 04-computation/jc2_cusp_residual_all_degree_gauge_kummer_parity_thm3894.py
output: 05-knowledge/results/jc2_cusp_residual_all_degree_gauge_kummer_parity_thm3894.out
script_sha256: a3e21a4c4e684848ef5f4fbf3758442f1cd01644120243a0900ad040de010b71
output_sha256: 0d7cb91145ce7cabd39852cdfc55ab72084aa4454092a48212f0259ee2292d6a
semantic_sha256: d8395de15e06826f7214f2ba22ef8c194543d56ff967145cc42126a4448af300
hash_basis: raw LF bytes
---

# THM-3894 -- the equality seam has one parity passport

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  Use the THM-3881 notation

```text
a=x+1,                 L=9x+4,
K=K_2+K_1+K_0,
K_2=y^2-15x^2,         K_1=-15x,             K_0=-4,      (1)
P=aL^2,
r=aT+Kf,               A=KT+aPf,
B=Pf^2-T^2.                                                (2)
```

Let `S(T,f)` be the residual from THM-3881.  Suppose it is a square and the
THM-3884 equality seam occurs:

```text
deg f=n>=1,            deg T=n+1,
f_n=xq_0,              T_(n+1)=-K_2q_0,                  (3)
```

where `q_0` is nonzero homogeneous of degree `n-1`.

## 1. The iterated gauge symbol

Put `q_(-2)=q_(-1)=0`.  Say that `(T,f)` is *aligned to depth `j>=1`* if
there are homogeneous polynomials

```text
deg q_i=n-1-i,                         0<=i<=j-1,          (4)
```

such that, for `0<=i<=j-1`,

```text
f_(n-i)=q_(i-1)+xq_i,
T_(n+1-i)=-K_2q_i+15xq_(i-1)+4q_(i-2).                   (5)
```

Depth one is exactly `(3)`.  Let the next unaligned pieces be

```text
s=f_(n-j),                t=T_(n-j+1),                    (6)
```

and define the homogeneous polynomial of degree `n-j+2`

```text
R_j=K_2(s-q_(j-1))
       +x(t-15xq_(j-1)-4q_(j-2)).                         (7)
```

Then `R_j` is exactly the highest remaining homogeneous part of `r`.
Indeed, the already aligned pieces cancel coefficient by coefficient in
`aT+Kf`; at degree `n-j+2`, substituting `(5)` leaves `(7)`.

The two other universal leading forms do not change with `j`:

```text
A_(n+4)=81q_0x^5,             B_(2n+3)=81q_0^2x^5.       (8)
```

Consequently the only two possible top contributions to the residual are

```text
3r^2B:  degree 4n+7-2j,   form 243q_0^2x^5R_j^2,
8AB:    degree 3n+7,      form 52488q_0^3x^10.           (9)
```

All remaining terms have degrees at most

```text
3n+6-2j,        2n+6,        n+5,        4,              (10)
```

and hence never tie either maximum in its stated regime.

## 2. The three regimes

### 2.1. Before the middle: another gauge jet

If `n>2j`, the first line of `(9)` strictly dominates.  If `R_j` were
nonzero, its `x`-valuation would be

```text
5+2v_x(q_0)+2v_x(R_j),                                    (11)
```

which is odd.  A nonzero highest form of a square is itself a square, so
`R_j=0`.  Since `gcd(x,K_2)=1`, reduction of `(7)` modulo `x` gives

```text
s=q_(j-1)+xq_j,
t=-K_2q_j+15xq_(j-1)+4q_(j-2),                           (12)
```

for a unique homogeneous `q_j` of degree `n-j-1`.  Thus alignment advances
by one full gauge jet.

### 2.2. At the middle: the Kummer symbol

If `n=2j`, the two lines of `(9)` tie and the complete top form is

```text
243q_0^2x^5(R_j^2+216q_0x^5).                            (13)
```

Its total degree `3n+7` is odd.  Therefore a square residual forces `(13)`
to vanish.  Unique factorization, followed by absorption of a scalar square,
gives the exact normal form

```text
q_0=xu^2,                R_j=epsilon*x^3u,
deg u=j-1,               epsilon^2=-216.                  (14)
```

Both signs in `(14)` are genuine at this filtration level.  The theorem does
not assert that either lifts.

### 2.3. Beyond the middle: the square profile

If `n<2j`, the second line of `(9)` strictly dominates.  Its degree is
`3n+7`.  Hence a square residual requires `n` odd; when `n` is odd, unique
factorization of

```text
52488q_0^3x^10                                             (15)
```

requires `q_0` itself to be a square up to a scalar, and the scalar is a
square because `k` is algebraically closed.

Iterating from depth one yields the promised parity passport.

* If `n=2m`, alignment is forced through depth `m`; at the middle the only
  possible leading symbol is `(14)` with `deg u=m-1`.
* If `n=2m+1`, alignment is forced through depth `m+1`; at the first layer
  beyond the middle, `q_0` must be a square.

Thus every survivor has `ceil(n/2)` aligned leading gauge jets before its
terminal parity symbol.  Later `q_i` are allowed to be zero; only `q_0` is
required to be nonzero.

The boundaries `n=1,2` are empty by THM-3886.  For `n>=3`, these conclusions
are necessary associated-graded conditions, not existence statements.

## 3. A vertical completion of the Kummer passport is impossible

There is one useful all-degree closure inside the remaining even passport.
Let `q,w in k[x]`, with `q!=0`, and consider the full gauge plus one residual
coordinate

```text
f=aq,                       T=-Kq+w.                      (16)
```

Then `S(T,f)` is not a square in `k[x,y]`.

To prove this, temporarily regard `K` as an independent coordinate over
`k(x)`.  Substitution in the exact residual makes `S` a quartic in `K`, and
direct elimination gives

```text
Disc_K(S)=128L^2q^8(3a^2q+2)^2
  *(-32L^4q+72L^2aqw^2+27w^4)
  *(-4L^2a^2q-2L^2+a^3qw^2)^3.                           (17)
```

As a polynomial in `y`, `S` has degree eight.  If it were a square, its
square root would have degree four and, because `S(x,-y)=S(x,y)`, would be
even in `y`.  It would therefore be a square already in `k[x,K]`, so its
quartic discriminant would vanish.

None of the three nonconstant factors in `(17)` can vanish identically.
The first and third specialize at `a=0` to `2` and `-2L^2`, respectively.
For the middle factor, evaluation at `L=0` first forces `L|w`; write `w=Lv`.
After canceling `L^4`, its vanishing becomes

```text
8q(9av^2-4)=-27v^4.                                      (18)
```

But

```text
gcd(v,9av^2-4)=1,                                         (19)
```

so `(18)` forces `9av^2-4` to be a unit.  This is impossible for nonzero
`v`, since its degree is `1+2deg v`; while `v=0` forces `q=0`.  This proves
the claim.  In particular, the full-gauge-plus-one-coordinate completion of
the first open even passport cannot remain vertical.  A completion leaving
the special form `(16)`, including general `y`-dependent lower data, remains
open.

## 4. Scope

The result is a filtration theorem.  It does not authorize subtracting the
aligned gauge polynomial: THM-3881 shows that such a subtraction changes the
residual by a `Delta`-debt.  It also does not prove that either parity
passport lifts, classify the `deg T>=deg f+2` region, produce a Keller atlas,
or settle JC(2).  What it does prove is that the equality seam has only two
all-degree terminal symbols, and that the most economical vertical
completion of the even symbol is empty.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_cusp_residual_all_degree_gauge_kummer_parity_thm3894.py
python3 -O 04-computation/jc2_cusp_residual_all_degree_gauge_kummer_parity_thm3894.py
```

Both streams must byte-match
`05-knowledge/results/jc2_cusp_residual_all_degree_gauge_kummer_parity_thm3894.out`.
