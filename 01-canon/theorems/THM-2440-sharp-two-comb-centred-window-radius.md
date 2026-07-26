---
id: THM-2440
title: "Sharp two-comb centred-window radius"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE
  AUDIT PENDING. If two radius-1/14 integer danger combs cover
  almost every point of {||ny||<rho} for rho>1/14, then both speeds
  are multiples of n and rho<=15/182. Equality is possible only for
  the scaled pair {n,13n}, which covers the open window almost
  everywhere but misses its two internal seam points. Thus no two
  combs cover {||ny||<8/91}. This is a standalone sharp geometry
  theorem and does not remove a new LRC(14) row.
source: codex-2026-07-26-two-comb-centred-radius
depends_on: []
related:
  - THM-1094-exact-two-comb-component-theorem
  - THM-1147-exact-two-comb-gap-law
  - THM-2434-endpoint-cage-gcd-invoice-for-guard-top-tilings
  - THM-2436-punctured-ninety-one-stalk-mixed-mode-and-repeated-step-closure
script: 04-computation/lrc14_sharp_two_comb_centred_radius_thm2440.py
output: 05-knowledge/results/lrc14_sharp_two_comb_centred_radius_thm2440.out
script_sha256: ff6045d45842d2daa5ffb7c97481be7b93495bb303e33b0cc377167c2d5f54fb
output_sha256: 774b94744cc4260be6295ac623d59eda7fdf682c70fcf3428f4997705049085f
hash_basis: working-tree bytes (LF)
---

# THM-2440 -- the seam-safe two-comb radius is `15/182`

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE
AUDIT PENDING.**

For a positive integer `u`, write

```text
D_u={y in R/Z: ||uy||<1/14}.                              (1)
```

The open endpoint convention is load-bearing. Two teeth may meet at
a point which neither danger set contains, so ordinary connected
components of `D_A union D_B` do not describe almost-everywhere
coverage. The theorem below permits those seams exactly.

## 1. The sharp theorem

Let `n,A,B` be positive integers and let `rho>1/14`. If

```text
{y:||ny||<rho} subset_ae D_A union D_B,                  (2)
```

then

```text
n|A,                         n|B,                        (3)

rho<=15/182.                                             (4)
```

If equality holds in (4), then, up to exchanging `A,B`,

```text
(A,B)=(n,13n).                                           (5)
```

Conversely, the pair in (5) covers

```text
{y:||ny||<15/182}
```

almost everywhere. Its only misses inside each pulled-back component
are the two seams corresponding to `ny=+-13/182` (hence `2n`
points on the original circle).

In particular,

```text
{y:||ny||<8/91} not subset_ae D_A union D_B              (6)
```

for every two positive integer speeds `A,B`, since
`8/91=16/182>15/182`.

## 2. Centre masks force the common pullback

Every centre `j/n` of the left side of (2) belongs to

```text
closed(D_A) union closed(D_B).                            (7)
```

Indeed, otherwise a neighborhood of that centre would have positive
uncovered measure.

For one speed `u`, put

```text
d=gcd(u,n),                     N=n/d.                    (8)
```

The values `uj/n` run through the `N`th roots, each with multiplicity
`d`. Hence the exact number of centres in `closed(D_u)` is

```text
d(2 floor(N/14)+1).                                      (9)
```

If `N>=3`, this is at most `n/3`. If `N=2`, it is `n/2`;
the hit set is the unique index-two subgroup of `Z/nZ`. Therefore two
proper masks cannot cover all `n` centres:

- two masks with `N>=3` cover at most `2n/3`;
- one `N=2` mask and one `N>=3` mask cover at most `5n/6`;
- two `N=2` masks coincide and cover only `n/2`.

Thus at least one of `A,B` is divisible by `n`. Say

```text
A=an.                                                     (10)
```

Coverage itself forces the other divisibility. First, (2) and the
measure bound `mu(D_A union D_B)<=2/7` give `rho<=1/7`.
For

```text
1/(14a)<t<min(rho,13/(14a)),                              (11)
```

the `n` points

```text
y_j=(j+t)/n,                     j in Z/nZ,                (12)
```

all lie in the window and all avoid `D_A`. If `n` did not divide
`B`, the values `B y_j` would contain `N>=2` equally spaced circle
points. They cannot all lie in the single arc `||x||<1/14`.
Consequently at least one `y_j` avoids `D_B` for every `t` in the
positive-length interval (11). Integrating over `t` and pigeonholing
the finite labels `j` produces a positive-measure violation of (2).
Therefore `n|B`, proving (3).

This argument uses open intervals only away from finitely many
endpoints. It does not turn almost-everywhere coverage into literal
endpoint coverage.

## 3. The reduced seam chain

Write

```text
A=an,                         B=bn
```

and exchange the two speeds so that `a<=b`. Pullback by `y -> ny`
preserves Haar measure, so (2) becomes

```text
{x:||x||<rho} subset_ae D_a union D_b.                   (13)
```

The positive central tooth of `D_a` ends at

```text
r_a=1/(14a),                                               (14)
```

and its next positive tooth begins at `13/(14a)`.
The central tooth of `D_b` is contained in the central tooth of
`D_a`.

If `b<13a`, the first noncentral `D_b` tooth begins at
`13/(14b)>r_a`. There is a positive gap, so (13) cannot pass
`r_a` even almost everywhere:

```text
rho<=1/(14a).                                             (15)
```

If `b>=13a`, any `D_b` tooth which touches or overlaps the central
`D_a` tooth has length `1/(7b)` and therefore ends no later than

```text
r_a+1/(7b).                                               (16)
```

The next `D_b` tooth is separated by a positive gap, while (16) is
strictly before the next `D_a` tooth. Thus no alternating tooth chain
can extend farther, and

```text
rho<=1/(14a)+1/(7b)
    <=1/(14a)+1/(91a)
     =15/(182a)
    <=15/182.                                             (17)
```

Equality throughout (17) forces `a=1,b=13`. For this pair,

```text
D_1 covers             |x|<13/182,

D_13 covers            13/182<|x|<15/182                 (18)
```

near zero. The points `x=+-13/182` are the two open seams.
This proves the converse, the equality statement, and (6).

## 4. Exact companion and controls

Run

```text
python 04-computation/lrc14_sharp_two_comb_centred_radius_thm2440.py
python -O 04-computation/lrc14_sharp_two_comb_centred_radius_thm2440.py
```

The dependency-free `Fraction` companion:

- checks the exact centre count (9) and every proper-mask inequality
  through denominator `1000`;
- constructs the actual almost-everywhere centred radius of
  `D_a union D_b` by merging exact tooth intervals, with touching
  endpoints deliberately treated as seams;
- checks all `1<=a<=b<=120`, finding the unique maximum
  `15/182` at `(1,13)`;
- applies the proved centre/divisibility gate and then the independent
  exact tooth merge to all `856,800` triples
  `1<=A<B,n<=120`, finding no cover of the target radius `8/91`;
- finds exactly the nine bounded equality triples
  `(A,B)=(n,13n)`, `1<=n<=9`; and
- retains THM-2434's hostile seam
  `D_3,D_11` at `5/14`: the endpoint is missed although its two
  punctured sides have different covering owners.

Normal and optimized transcripts must reproduce

```text
05-knowledge/results/lrc14_sharp_two_comb_centred_radius_thm2440.out
```

byte-for-byte.

## 5. Scope

THM-2440 is a standalone sharp one-dimensional geometry theorem. It
repairs the tempting but false inference that an almost-everywhere
open cover catches every endpoint: only closed centre masks are used,
and the equality example itself has two missed seams.

The `8/91` corollary would close the proposed exceptional step-one
two-blocker route, but hostile-audited THM-2436 already empties the
entire deep-`c_3` branch by a stronger punctured-stalk atlas. Thus
THM-2440 supplies a reusable sharp sidecar and an equality
classification; it does not remove another valuation shape or scalar
row, prove owner-conditioned survival in the live `c_3<=M` graft, or
prove LRC(14).

QED conditional on independent audit.
