# Independent referee report: proposed THM-4088

**Verdict:** **PASS**, subject to the two scope conditions below.  The result
is elementary but useful: it upgrades the exact interval-margin mechanism
already present in THM-613/HYP-3114 to an arithmetic-type density statement,
and it makes the scalar order-tournament loss exact.

## Required scope conditions

1. Say **scalar strict-order tournament**: vertices carry distinct real
   labels and `i -> j` means `x_i < x_j` (or the globally reversed
   convention).  This conclusion does not apply to an intrinsically
   arithmetic orientation with extra edge labels.
2. Say **strict / positive-margin witness**.  The density conclusion is false
   for a boundary-only witness with `F_S(t)=b`.  Arithmetic type is also a
   statement in the chosen time coordinate; arbitrary real rescaling can
   change it.

## Inheritance pass

- **Closest proved mechanism:** the proved first part of
  THM-613, `margin-measure-slope-bridge`, already establishes that
  `F_S(t)=min_{v in S} ||vt||` is `B=max S` Lipschitz and fattens a maximizing
  witness into an interval.
- **Canonical hostile:** THM-358, `lrc-initial-segment-unit-skeleton`, says
  that for `S={1,...,n-1}` at threshold `1/n`, the entire closed safe set is
  the finite rational unit skeleton `{a/n:gcd(a,n)=1}`.  For `n=14` there is
  no strict witness and no irrational equality witness.
- **Corrected near miss:** MISTAKE-170 / THM-1158,
  `invcov-refuted-by-doubled-ap`, shows that `{1,...,13}` and its common
  dilation have the same transitive speed-order tournament while the
  auxiliary Covering predicate changes.  Rank order is not a divisor or
  arithmetic carrier.
- **Least-used relevant sidecar:** HYP-3114 retains the witness interval,
  margin, maximum speed, continued-fraction entry and exceptional
  approximants.  THM-4057 supplies Stern--Brocot interval addresses.  The new
  statement needs only the interval/margin part; CF and arithmetic-height
  data matter only for a prescribed denominator or named number.
- **Method card:** “Tournament Analysis must preserve content.”  Here the
  unadorned scalar-order tournament provably preserves only rank and destroys
  precisely the claimed arithmetic data.

## Canon-ready theorem

Let `V` be a finite label set.

### A. Exact image and fibres of scalar order tournaments

For pairwise distinct real labels `(x_v)_(v in V)`, orient

```text
u -> v  iff  x_u < x_v.                                  (1)
```

Then the resulting tournament is transitive.  Conversely, every transitive
tournament on `V` arises in this way.

More strongly, fix a transitive tournament `T` and assign to every vertex one
of the following arithmetic types:

```text
rational;
real algebraic of an arbitrarily prescribed exact degree d >= 2;
transcendental.                                          (2)
```

There are pairwise distinct labels of exactly those prescribed types whose
order tournament is `T`.  Hence the fibre of the scalar-order map over every
transitive tournament meets every finite prescribed arithmetic-type stratum.
In particular, the unlabelled orientation cannot determine whether any
vertex is rational, quadratic irrational, or transcendental.

### B. Positive-margin LRC witness density

Let `S` be a nonempty finite set of positive integers, put

```text
F_S(t)=min_(v in S) ||v t||,    B=max S,
W_b^>(S)={t in R/Z : F_S(t)>b}.                           (3)
```

If `t_0 in W_b^>(S)` and

```text
eta=F_S(t_0)-b>0,                                        (4)
```

then

```text
d_R/Z(t,t_0)<eta/B  =>  t in W_b^>(S).                   (5)
```

Consequently `W_b^>(S)` is open, and if it is nonempty then each of

```text
W_b^>(S) intersect Q/Z,
W_b^>(S) intersect {quadratic irrationals mod Z},
W_b^>(S) intersect {transcendentals mod Z}               (6)
```

is dense in `W_b^>(S)` in its relative topology.  The same holds with
“quadratic” replaced by “real algebraic of exact degree `d`” for any fixed
`d>=2`.

For LRC(14), take `b=1/14`.  This says that any already-proved strict witness
has representatives of all three arithmetic types arbitrarily nearby.  It
does not prove that such a strict witness exists for an open row.

### C. Sharp boundary

The hypothesis `eta>0` cannot be weakened to `F_S(t_0)>=b`.  At

```text
S={1,...,13},    b=1/14,                                 (7)
```

THM-358 gives

```text
{t:F_S(t)>=1/14}={a/14:gcd(a,14)=1},
{t:F_S(t)>1/14}=empty.                                   (8)
```

Thus the equality set is finite and rational.  Also the strict radius in
(5) is sharp from Lipschitz information alone: for `S={1}`, `b=1/4`, and
`t_0=1/2`, one has `eta/B=1/4`, while both endpoints `1/4,3/4` have value
exactly `b`, not strictly above it.

## Proof

For A, (1) inherits transitivity from `<`.  A transitive tournament has a
unique topological total order, so any strictly increasing labels realize it.
For every fixed `d>=2`, the real algebraic numbers of exact degree `d` are
dense: if `alpha=2^(1/d)`, then `q+r alpha` has exact degree `d` for
`q,r in Q`, `r!=0`, and such values can be placed in any open interval.
Rationals are dense.  Transcendentals are dense because the algebraic numbers
are countable while every open interval is uncountable.  Choose one point of
the prescribed type in each of disjoint open intervals ordered according to
`T`.

For B, distance to a fixed closed subset is 1-Lipschitz, hence

```text
| ||v t||-||v u|| | <= |v| d_R/Z(t,u).                   (9)
```

The minimum of finitely many functions with common Lipschitz constant `B` is
`B`-Lipschitz.  Therefore

```text
F_S(t) >= F_S(t_0)-B d_R/Z(t,t_0)>F_S(t_0)-eta=b,         (10)
```

which proves (5).  Each arithmetic class in (6), and each exact algebraic
degree class, is dense in the circle because the corresponding class is
dense in every lifted real interval.  Intersecting those dense classes with
the open set `W_b^>(S)` proves relative density.  C is THM-358 plus the
displayed one-speed calculation.

## Optional Stern--Brocot corollary

This is worth one short corollary, not a headline theorem.  Every nonempty
open arc in `(0,1)` contains a complete finite Stern--Brocot interval: take an
irrational point in the arc and follow its nested Stern--Brocot addresses;
the adjacent-boundary denominators tend to infinity, so the interval widths
`1/(bd)` tend to zero.  Hence every positive-margin LRC witness arc contains
a Stern--Brocot interval all of whose interior is still strict, and that
interval contains rational, exact-degree algebraic, and transcendental
points.  The finite address therefore classifies a location interval, not an
arithmetic type.  Do not say that one finite continued-fraction word decides
rationality, periodicity, or transcendence.

## Connection ledger

| source | target/map | preserves | destroys | required sidecar | hostile |
|---|---|---|---|---|---|
| distinct real labels | scalar order tournament by pairwise `<` | rank permutation | gaps, field degree, height, arithmetic type | numerical labels/type data | THM-1158 common dilation |
| strict LRC time `t_0` | open ball of radius `eta/B` | `F_S>b` | exact time and Diophantine quality | `(eta,B)` and fixed normalization | THM-358 AP13 boundary |
| open witness arc | contained Stern--Brocot interval | strict witness predicate | chosen representative/type | finite address plus margin | boundary-only unit skeleton |

## Status and dependency recommendation

Promote THM-4088 as **PROVED (elementary)** with dependencies

```text
depends_on: [THM-358, THM-613]
related: [THM-1158, THM-4057, HYP-3114]
```

The proof can be self-contained, but these links record inheritance and the
sharp boundary.  Do not claim progress on the existence part of LRC(14), a
classification of rational/irrational/transcendental numbers by tournaments,
or invariance of arithmetic type under time rescaling.

## Independent exact probe

Reproduce from the session worktree with

```text
python -B .scratch/lrc_density_referee_20260825/lrc_density_exact_probe.py
python -O -B .scratch/lrc_density_referee_20260825/lrc_density_exact_probe.py
```

The probe uses exact `Fraction` arithmetic.  It checks the complete rational
grid equality skeleton for `S={1,...,n-1}`, `2<=n<=14`; exhausts the circle
Lipschitz inequality over every nonempty subset of `{1,...,7}` on the
`1/84` grid; tests every induced positive-margin implication there; and
checks the sharp one-speed endpoint hostile.  It is a hostile control only;
the proof above is the theorem.

Normal and optimized replays agree exactly:

```text
THM-4088 independent exact probe
tight-initial-segment grid values: 6240
circle-Lipschitz pairs: 896112
strict-margin implications: 22935
open-radius endpoint hostile: PASS
PASS
```

SHA-256 after the agreeing replays:

```text
lrc_density_exact_probe.py
  304ef9b7bc69b83c600064988a2b3bdc953e129d5e496200f9e16cf31b7e516b
```
