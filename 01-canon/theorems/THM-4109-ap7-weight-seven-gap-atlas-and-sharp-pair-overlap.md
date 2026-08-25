---
id: THM-4109
title: "AP7 weight-seven gap atlas and sharp selected-pair overlap"
status: >
  PROVED + PROVED RELATIVE TO THM-2072/4092/4101 + VERIFIED-EXACT +
  INDEPENDENTLY VERIFIED-EXACT. On the exact AP7 safe interval, a four-outlier
  weight-seven bank with one even and three odd speeds is absorbable whenever
  two odd speeds differ by 4 and every outlier is at least 197. Gap 8 is
  absorbable from height 232. These are the first exact uniform heights for
  this selected-pair second-moment certificate; the previous heights have
  explicit negative-margin hostiles. Gaps 4,8,12 are the only whole-interval
  resonances, gaps 2,6,10 have zero overlap, and no other even gap improves
  height 264 within this certificate. This is neither arbitrary weight-seven
  absorption nor LRC(14).
source: codex-arithmetic-boundary-breakthrough-20260825
depends_on:
  - THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate
  - THM-4092-parity-weighted-antipodal-k-comb-density-compiler
  - THM-4101-ap7-weight-seven-gap-four-second-moment-absorption
related:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4079-lrc14-antipodal-outlier-absorption-and-adaptive-clock
  - THM-4098-weight-seven-antipodal-scale-escape-and-missing-parity-rows
  - THM-4100-residual-component-three-outlier-lrc-compiler
  - MISTAKE-402
script: 04-computation/lrc_ap7_even_gap_overlap_atlas_thm4109.py
output: 05-knowledge/results/lrc_ap7_even_gap_overlap_atlas_thm4109.out
independent_audit_script: 04-computation/lrc_ap7_even_gap_overlap_atlas_thm4109_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_ap7_even_gap_overlap_atlas_thm4109_independent_audit.out
script_sha256: 505784e664c703cbcee7f1f66d579e3d203ae88d9970a9dde2598b4655719d56
output_sha256: 391df0ccf5aaf4b5584ac359f97a53930e6006024bb9686a829b93a2f4156952
independent_audit_script_sha256: ba940638a2d8a0667c1ccfb1da045f98707cbbaa5b43dd5e874ebccc6b9b6add
independent_audit_output_sha256: f1a163991a220164effedf9e92880da0b1e19b29b08685daf0a4300cf08764c9
hash_basis: raw LF bytes
audit: >
  PASS. The primary uses merged exact comb teeth, wall quadrature for the
  triangular mean, all 490 odd residue classes modulo 980 for each resonant
  gap, positive quadratic ray gates, and exact height-263 hostiles. The
  independent path uses literal theta-wall cells, a closed antiderivative,
  explicit left/right packet endpoints, remote derivation plus backward and
  forward checks of every residue law, and independently searched hostiles.
  It checks 11,760 endpoint/quasipolynomial gates. After CRLF-to-LF
  normalization, normal and optimized runs match both frozen outputs exactly.
---

# THM-4109 -- the AP7 weight-seven selected-pair gap atlas

Put `delta=1/14` and retain THM-4092's exact two-phase AP7 safe interval

```text
D={1,...,7},             J=[4/35,13/98],
L=|J|=9/490,             J subset G_D^+-.                (1)
```

For a speed `v`, write

```text
U_v={theta:min(||v theta||,||v(theta+1/2)||)<delta},
omega(v)=1 if v is even, and 2 if v is odd.              (2)
```

This theorem refines THM-4101 in two directions: it computes the actual
selected odd-pair overlap rather than a uniform lower envelope, and it
classifies every even gap that could improve the old height `264` inside the
same second-moment certificate. A failed margin below is only a certificate
hostile unless coverage is separately checked.

## 1. The selected-pair second-moment criterion

Let the four distinct outliers be

```text
V={e,w,u,u+g},
e even,                 u,w,u+g odd,                    (3)
```

so `g` is positive and even. Put

```text
M_g(u)=mu(J intersect U_u intersect U_(u+g)).            (4)
```

THM-4092's one-comb estimates give

```text
sum_(v in V) mu(J intersect U_v)
 <=L+T(e,w,u,g),                                        (5)

T(e,w,u,g)
 =6/(49e)+(5/49)(1/w+1/u+1/(u+g)).                     (6)
```

The four-comb pointwise inequality from THM-4101 is

```text
1_(F>0)<=F-(1/2)binom(F,2),       0<=F<=4.              (7)
```

Retaining only the selected pair in the nonnegative second factorial moment
therefore gives

```text
mu(J intersect union_(v in V)U_v)
 <=L+T(e,w,u,g)-(1/2)M_g(u).                            (8)
```

Consequently

```text
boxed: M_g(u)/2>T(e,w,u,g)                              (9)
```

is a strict survivor certificate for `D union V`.

## 2. A moving-target formula and convergence theorem

For odd `u`, put

```text
n=2u,            A=(-1/7,1/7) in R/Z,
x=n theta,       s=2g theta.                            (10)
```

The two selected danger conditions are

```text
x in A,          x+s in A.                              (11)
```

Thus their instantaneous phase-intersection length is the triangular wave

```text
H(s)=max(0,2/7-||s||).                                  (12)
```

Define its exact average on `J` by

```text
Q_g=int_J H(2g theta)d theta
   =(1/(2g))int_(2gJ)H(s)ds.                            (13)
```

> **Lemma 2.1 (uniform moving-target error).** For every positive even `g`
> and positive odd `u`,
>
> ```text
> Q_g-(5gL/4+2/7)/u <= M_g(u)
>                    <= Q_g+(5gL/4+1)/u.                (14)
> ```
>
> In particular `M_g(u)->Q_g` for fixed `g` as `u->infinity` through odd
> integers.

### Proof

Partition `J` into complete `1/n`-periods for `x=n theta` and at most two
boundary fragments. On a complete period `I` with midpoint `c`, the shift
`2g theta` moves by circular distance at most `g/n` from its midpoint value.
Intersecting all moving target arcs over `I` can remove at most `2g/n` of
phase length from `H(2gc)`. Pullback under `x` divides phase length by `n`,
so this costs at most `2g/n^2` on `I`.

Symmetrically, every phase in a moving intersection lies in the midpoint
intersection after enlarging its shifted target arc by `g/n` at each
endpoint. That enlargement gains at most `2g/n` of phase length, hence at
most `2g/n^2` after pullback. Thus the same complete-period error controls
both sides.

The function `theta -> H(2g theta)` is `2g`-Lipschitz. Midpoint quadrature on
an interval of length `1/n` therefore costs at most

```text
int_I 2g|theta-c|d theta=g/(2n^2).                      (15)
```

There are at most `nL` complete periods. Their total lower or upper error is
at most `5gL/(2n)`. The two discarded fragments have total length below
`2/n`; their `H`-mass is at most `4/(7n)`, while their literal overlap mass
is at most `2/n`. Substitute `n=2u` to obtain the two sides of `(14)`. QED.

## 3. Complete resonance classification

The limiting overlap `Q_g` is positive exactly when `2gJ` has
positive-length intersection with an interval

```text
(k-2/7,k+2/7),            k in Z.                       (16)
```

For the seven even gaps through `14`, direct rational endpoint comparison
shows

```text
Q_g=0 exactly for g in {2,6,10}.                        (17)
```

For `g>=16`, the interval `2gJ` has length

```text
2gL=9g/245>3/7,                                         (18)
```

longer than every complementary gap of the support of `H`, so `Q_g>0`.
This proves `(17)` for all even gaps.

A stronger resonance occurs when `2gJ` contains one integer and lies wholly
inside its support interval in `(16)`. Whole-interval containment forces
`9g/245<4/7`, hence `g<=14`; exact inspection then gives only

```text
(g,k)=(4,1),(8,2),(12,3).                              (19)
```

All three resonance integers pull back to the same internal split
`theta=1/8`.

## 4. The exact mod-980 overlap law

Let `g=4h` be one of `4,8,12`, put `k=h`, `n=2u`, and `m=n+2g`. Before the
split point, the `p`-tooth for `n` meets the `(p+k)`-tooth for `m`; after the
split the opposite endpoint pair is active. The unclipped components are

```text
I_p^-=((p+k-1/7)/m,(p+1/7)/n)       on theta<=1/8,
I_p^+=((p-1/7)/n,(p+k+1/7)/m)       on theta>=1/8.       (20)
```

On each side the live indices `p` form one consecutive integer interval;
only the two extreme components can be clipped by the endpoints of `J`.
Every switch index is a floor or ceiling of an affine form at one of

```text
4/35, 1/8, 13/98, 1/7.                                 (21)
```

Their denominators divide `1960` in `n=2u`, hence all fractional parts and
branch choices repeat under `u -> u+980`. On a fixed odd residue
`r mod 980`, each switch index is affine in the ray parameter. Summing the
arithmetic progressions of component lengths in `(20)` first gives a
degree-at-most-two numerator over `u(u+g)`. Lemma 2.1 gives
`M_g(u)-Q_g=O(1/u)` on the same ray, so the quadratic coefficient after
subtracting `Q_g` must vanish. Therefore:

> **Lemma 4.1 (exact residue law).** For `g in {4,8,12}` and each odd
> `r mod 980`, there are rational `alpha_(g,r),beta_(g,r)` such that every
> positive odd `u congruent r mod 980` satisfies
>
> ```text
> u(u+g)(M_g(u)-Q_g)=alpha_(g,r)u+beta_(g,r).            (22)
> ```

This is an infinite quasipolynomial identity, not an interpolation claim.
The endpoint derivation above proves its form. The primary certificate
constructs all `490` residue pairs for each gap and checks two second
differences plus a held-out fourth point. The independent audit derives each
pair from remote samples, checks backward and forward samples, and compares
literal theta-wall cells with the endpoint sum on both a low and a remote
representative of every residue class.

Exact integration of `(13)` gives

```text
Q_4 =2187/480200,
Q_8 = 927/240100,
Q_12=1521/480200.                                      (23)
```

## 5. Sharp uniform heights for this certificate

Fix a lower height `B` for all four outliers. The largest reciprocal taxes
in `(6)` occur at the least allowed even speed and least available odd
speeds. There are two exhaustive placements:

1. the selected pair starts at the least odd speed, and the third odd speed
   is the least distinct admissible one; or
2. the selected pair starts later, and the third odd speed is the least odd
   speed.

In the second case, substitute `(22)` into `(9)` and multiply by `u(u+g)`.
On each of the 490 residue rays the result is a quadratic in `u`. Its leading
coefficient, value, and derivative at the first admissible point are all
strictly positive at the certified heights. Hence every later point on the
ray is positive. The first case is checked directly. The exact outcomes are

| gap | first height | least boundary proof gate | gate margin |
|---:|---:|---|---:|
| 4 | **197** | `{197,198,203,207}` | `68161/17847663372` |
| 8 | **232** | `{232,233,235,241}` | `11027/15001161644` |
| 12 | `268` | `{268,269,271,283}` | `20945/2764475878` |

At the previous heights the selected-pair margins are nonpositive:

| gap | previous height | certificate hostile | margin |
|---:|---:|---|---:|
| 4 | 196 | `{196,197,203,207}` | `-28297/11357603964` |
| 8 | 231 | `{231,232,233,241}` | `-500165/73729113612` |
| 12 | 267 | `{267,268,269,279}` | `-54/1604047` |

Therefore the heights `197,232,268` are exact first uniform heights for
this selected-pair certificate. In particular, gap four strengthens
THM-4101's sufficient height `264` to `197`, and gap eight supplies a new
family from height `232`. The negative rows prove sharpness of the
certificate, not failure of actual absorption.

## 6. No other even gap improves 264 in this certificate

At height `263`, retain the least third odd and even speeds `263,264`, and
let the selected pair start tend to infinity. Lemma 2.1 makes the limiting
fixed tax

```text
C_263=5/(49*263)+6/(49*264)=69/81004.                   (24)
```

Thus any even gap capable of a uniform height at most `263` must survive the
asymptotic comparison `Q_g/2>=C_263`. Exact integration for `g<202` has no
equality case and leaves only

```text
g in {4,8,12,16,20,34,38,42,46}.                       (25)
```

For `g>=202`, a complete-period area plus at most one residual triangle gives

```text
Q_g<=4L/49+2/(49g)
   <=4L/49+2/(49*202)<2C_263,                            (26)
```

so no larger gap qualifies. Finally, exact negative-margin banks at height
`263` eliminate every candidate in `(25)` except `4,8`:

| gap | certificate hostile | margin |
|---:|---|---:|
| 12 | `{263,264,265,275}` | `-37921/751312100` |
| 16 | `{263,264,265,279}` | `-334935/931627004` |
| 20 | `{263,264,265,283}` | `-4901605/8504852972` |
| 34 | `{263,264,265,297}` | `-734207/1014271335` |
| 38 | `{263,264,267,305}` | `-32696563/46175925180` |
| 42 | `{263,264,267,309}` | `-11066945/15593837028` |
| 46 | `{263,264,267,313}` | `-17389205/23693548494` |

A failed bank at height `263` belongs to every lower-height universe, so
these rows prove that no other even gap improves height `264` through the
criterion `(9)`.

## 7. Endpoint clock, hostiles, and precise scope

Strict inequality in `(9)` leaves a nonempty closed survivor in `J` after
removing the four open danger combs. An endpoint of a survivor component is
rational with an even presentation whose denominator is at most
`14 max(V)`. THM-2072 therefore supplies the associated endpoint clock and,
with the inherited owner atlas, closes the corresponding dyadic two-tail
seam through every common dilation.

The canonical bank

```text
{8,9,11,13}                                             (27)
```

still covers all of `J`. Its gap-two selected pairs have zero overlap, and
its gap-four selected-pair margin is negative. This is the hostile boundary
between an actual cover and the new sufficient families. By contrast, the
negative rows in Sections 5--6 are only certificate hostiles unless a
separate union computation proves coverage.

THM-4100 is orthogonal: it absorbs three arbitrary-parity outliers over AP8
from height `93`; the present theorem absorbs four weight-seven outliers over
AP7 by retaining an odd-pair overlap. Neither theorem extracts its core from
an arbitrary residual. General unscaled weight-seven absorption, physical
entry from the LRC(14) reduction, and LRC(14) itself remain **OPEN**.

**QED.**
