---
id: THM-3474
title: "Binary-submask factor degrees and all-height prime-power reset families"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  For every even N, the
  complete 2-adic lower Newton polygon of A_N^(N+1) has one edge for each
  set bit of N; its necessary factor-degree barcode is exactly the set of
  binary submasks of N.  Intersecting this with THM-3467's complete p^k
  reset gives an exact bit-test coprimality compiler.  In particular, for
  every k>=1 it closes d=2^s p^k+1 for every power of two 2<=2^s<p.
  These are infinite
  sparse families of exact quadratic three-moment windows, not SFC or FC.
audit: >
  The binary polygon is a closed-form corollary of THM-3161's proved digital
  skeleton, not a new coefficient-unit theorem.  Fraction-hull and
  self-contained determinant-hull engines agree with all predicted hulls
  and submask barcodes for 60 even N through 120.  They also agree on 15
  cross-place reset controls, all with empty positive intersections, retain
  planted factors v and v+1, and reproduce the sharp odd-N=5 hostile.
source: root/factorial-jacobian-alternation/2026-08-15
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3161-factorial-newton-euclidean-closure-through-r1998
  - THM-3467-factorial-seven-exit-newton-barcode-extension-after-r2498
related:
  - THM-3142-prime-power-predecessor-newton-separation-and-composite-window-census
  - THM-3152-multi-place-newton-degree-barcode-and-euclidean-flag-census
script: 04-computation/factorial_binary_submask_reset_families_thm3474.py
output: 05-knowledge/results/factorial_binary_submask_reset_families_thm3474.out
script_sha256: ffe323328fc0fcf41ca5166ecf851a9b1b23e19b11772d6e4cdf18e4d04fedac
output_sha256: edb29cc1e2b982adeaac7b283f725ace8ef62c229b10208010cfb2f7a53b96cd
semantic_sha256: 49272719898eec5b506e435140c0ce6aac1ce0b07ce275fff648df55b18b61f8
hash_basis: raw bytes
---

# THM-3474 -- binary submasks and all-height reset families

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Complete binary polygon

Use the exact factorial moments

```text
L(x^r)=r!,
A_n^(d)(v)=L((d-x+v x^2)^n).                               (1)
```

Let `N` be a positive even integer and put

```text
G_N(v)=A_N^(N+1)(v).                                       (2)
```

Write the set bits of `N` as

```text
N=sum_(q=1)^s 2^(e_q),            1<=e_1<...<e_s,
J_0=0,                            J_q=sum_(i=1)^q 2^(e_i).  (3)
```

Then the complete lower `2`-adic Newton polygon of `G_N` has vertices

```text
(J_q,2J_q-q),                       0<=q<=s.                (4)
```

Its `q`-th maximal edge has

```text
slope       (2^(e_q+1)-1)/2^(e_q),
capacity    2^(e_q),
denominator 2^(e_q).                                         (5)
```

The constant coefficient is nonzero, so the coordinate-root capacity is
zero.  Consequently THM-3152's exact local degree barcode is

```text
D_2(G_N)={r: 0<=r<=N and r & ~N = 0}.                       (6)
```

Thus every rational factor of `G_N` has degree equal to a binary submask of
`N`.  Equation (6) is a necessary degree address.  It neither constructs a
factor for every submask nor makes a nonempty address sufficient.

## 2. Proof from the inherited digital skeleton

THM-3161 proves, for every prime dividing `N`, that the actual Newton polygon
is the lower hull of its digital weight.  At `p=2`, that weight simplifies to

```text
w_N(j)=nu_2 C(N,j)+nu_2((2j)!)
      =2j+s_2(N-j)-s_2(N),                  0<=j<=N,        (7)
```

where `s_2` is binary digit sum.  THM-3161 also proves exact coefficient
valuation at every vertex of this formal lower hull.  The new step is to
solve that hull in closed form.

At `j=J_q`, equation (7) gives

```text
w_N(J_q)=2J_q-q.                                           (8)
```

Between consecutive prefix sums, write

```text
j=J_(q-1)+u,                      0<u<2^(e_q).              (9)
```

The binary supports in

```text
N-j=(N-J_q)+(2^(e_q)-u)                                   (10)
```

are disjoint.  Hence

```text
w_N(j)=2j-q+s_2(2^(e_q)-u).                                (11)
```

The chord between the two points in (8) has height

```text
2j-q+1-u/2^(e_q).                                         (12)
```

Subtracting (12) from (11) gives

```text
s_2(2^(e_q)-u)-1+u/2^(e_q)>0.                              (13)
```

Thus every interior point lies strictly above the chord.  The slopes in (5)
strictly increase with `q`, so (4) is the entire convex lower hull.  The
exactness clause inherited from THM-3161 passes it to the actual coefficient
polygon.

Each numerator in (5) is odd.  Its reduced denominator equals its horizontal
capacity, so a local factor may use either zero or the entire `2^(e_q)` block.
Summing independently over the edges gives precisely the submasks in (6).
QED.

For completeness, the coefficient-unit mechanism behind the inherited
exactness is especially transparent here.  If

```text
[v^j]G_N=C(N,j)(2j)! Z_j,
```

then at `j=J_q` the `ell=0` term of `Z_j` is odd; the `ell=1` term has the
even factor `N-J_q` when present, and every `ell>=2` rising product contains
an even factor.  Hence `Z_(J_q)==1 (mod 2)`.

## 3. Cross-place bit compiler

Let `p>=5` be prime, `k>=1`,

```text
H=p^k,               m=(p-1)/2,               N=aH,       (14)
```

and suppose `a` is even with

```text
2<=a<p,                    T=min(a-1,m).                    (15)
```

Put `d=N+1`, `F=A_(N-1)^(d)`, `G=A_N^(d)`, and let `E` be the first full
Euclidean row.  THM-3467 gives the complete common `p`-adic barcode

```text
D_p(F,G,E)={0,H,2H,...,TH}.                                (16)
```

Since `N` is even, the present theorem applies to the same row `G`.  Every
positive rational common-factor degree must therefore be

```text
r=tH,       1<=t<=T,       and       (tH) & ~(aH)=0.       (17)
```

This yields an exact, finite bit compiler:

> If no `t` in the displayed range satisfies the submask test in (17), then
> `gcd_Q(F,G)=1`.

Indeed, a common factor divides `G` and hence obeys (6), while the Euclidean
flag preserves it and hence forces (16).  Empty intersection excludes every
positive degree by THM-3152.  QED.

## 4. Infinite power-of-two families

Suppose

```text
a=2^ell,                         ell>=1.                    (18)
```

Since `H` is odd, `aH` is the binary expansion of `H` shifted left by `ell`.
Every binary submask of `aH` is divisible by `a`.  But for `1<=t<a`,

```text
nu_2(tH)=nu_2(t)<ell,                                      (19)
```

so `tH` is not divisible by `a` and cannot be such a submask.

For every `2<=a<p`, one has `T=min(a-1,m)<a`, so (19) excludes every
candidate in (17).  Therefore

```text
gcd_Q(A_(a p^k-1)^(a p^k+1),A_(a p^k)^(a p^k+1))=1        (20)
```

for every `k>=1` whenever `a` is a power of two smaller than `p`.

By THM-3124's exact resonance reduction, if

```text
q(x)=A+Bx+Cx^2,                         ABC!=0,             (21)
```

then at least one of

```text
L(q^(a p^k-1)),          L(q^(a p^k)),
L(q^(a p^k+1))                                             (22)
```

is nonzero.  This proves infinitely many all-height-in-`k` windows for each
admissible pair `(p,a)`.  Examples include

```text
d=2*5^k+1,        d=4*5^k+1,        d=4*7^k+1,             (23)
```

as well as every other power-of-two multiplier below `p`.  In particular,
the whole family `d=4*5^k+1` is structural; `d=2501` is its `k=4` member.

These are sparse exact-support quadratic windows.  They do not give a
contiguous all-height theorem, arbitrary one-variable support, `SFC(3)`, or
the multivariable Factorial Conjecture.

## 5. Hostiles and failure boundaries

- **Parity is load-bearing.**  Formally applying (4) to odd `N=5` predicts
  `((0,0),(1,1),(5,8))`, but the actual polygon of `A_5^(6)` is
  `((0,3),(5,8))`; the vertex-unit mechanism fails.
- **Reset completeness is load-bearing.**  The compiler uses THM-3467's
  complete tail separation; the first common face alone would not suffice.
- **A surviving bit is inconclusive.**  Some non-power-of-two multipliers
  pass (17); that is only an observer address, never a constructed factor.
- **The two places act together.**  Neither a binary submask condition nor a
  `p^k` divisibility condition alone proves coprimality.

Power-of-two `a` is a clean sufficient family, not a necessary condition for
the compiler to close a row.  Testing (17) directly gives additional isolated
closures.

## 6. Exact verification

The companion hash-pins two different inherited engines: THM-3152's
Fraction lower hull and THM-3180's self-contained determinant lower hull.  It
checks (4)--(6) for all `60` even integers `2<=N<=120`, with semantic digest

```text
888b5d72e72bfeae0ef53140420e728fc5c7a78bfab60e6ac29e694a0ee60edd. (24)
```

It then checks `15` power-of-two reset cells, using both engines at both
places and requiring empty positive intersections.  Their digest is

```text
3512b41ae23e168e5c500a989f37b6120e2dee7670c700b1f40f32b2fe3c45e0. (25)
```

Planted factors `v+1` and `v` retain degree one.  The odd `N=5` hostile is
reconstructed independently.  Run

```text
python3 04-computation/factorial_binary_submask_reset_families_thm3474.py
python3 -O 04-computation/factorial_binary_submask_reset_families_thm3474.py
```

and compare raw bytes with the declared output.  The script contains no
Python `assert`, floating-point, or random truth gate.

**QED.**
