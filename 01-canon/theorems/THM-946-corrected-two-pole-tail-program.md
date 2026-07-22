---
id: THM-946
title: THE CORRECTED TWO-POLE TAIL PROGRAM — a rigorous discrete-convolution bound, the refuted logarithm-free claim, and the still-conditional B5 exhaustion route
status: PARTIAL / TARGET BYPASSED BY THM-2051. The corrected two-pole convolution lemma below is proved with an explicit constant 64. The former C0=8 bound is refuted exactly. T3 still needs a punctured near-pole congruence estimate; T4 and T5 still need affine resonance-strip and resonance-slab estimates as absolute relation-mass bounds. However, THM-2051 now proves the coarse alternative “support-at-most-five relation of height at most 2^20 OR continuous B5 > 0” by whole-product Fejer--BV approximation, without these estimates. The small-relation branch remains unclassified.
source: kind-pasteur-2026-07-17-S128 cont.38, corrected by codex-S48 audit
depends_on:
  - THM-935 (relation-mass identity and the proved support-2 tail)
  - THM-930 (leverage identity frame)
related:
  - THM-944 (the first-pushed B5 race scoreboard; this page was renumbered after collision)
  - THM-945/947 (moment wall and arc wire)
  - HYP-7189 (audit-corrected tail program)
  - THM-2051 (whole-product bypass for the coarse universal alternative)
---

# THM-946 — corrected two-pole tail program

> **2026-07-21 update (THM-2051).** The estimates on this page remain valid
> targets if one wants termwise absolute control or a much smaller relation
> horizon. They are no longer a prerequisite for the coarse universal
> small-relation-or-BONF5 alternative: THM-2051 keeps each signed centered
> product intact, kills its bounded Fejer part by character orthogonality, and
> absorbs every higher relation collectively in one `L1` error.

This page supersedes the collided page formerly named
`THM-944-ts-tail-lemma-universal-exhaustion.md`.  The old page made two distinct
overclaims:

1. its two-pole estimate omitted a necessary `log(2 + Delta)` factor and is false;
2. even after correcting that atom, its one-paragraph `T4/T5` nesting does not prove
   the required affine-coset estimates.

The rigorous result is the corrected convolution lemma below.  The rest of the page
records exactly what remains before it can feed THM-935.

## I. Corrected discrete two-pole convolution lemma — proved

Let `A,B >= 1`, let `alpha,beta` be real, and put

```text
Delta = |alpha - beta|,
S(A,B;alpha,beta)
  = sum_(m in Z) 1 /
      (max(1,A|m-alpha|) max(1,B|m-beta|)).
```

Then

```text
S(A,B;alpha,beta)
  <= 64(1 + log(2+Delta)) / (A B (1+Delta))
     + 6/(1+A Delta) + 6/(1+B Delta).                 (1)
```

Coprimality is not needed for this analytic estimate.  It enters later only when an
integer relation line is reduced to a primitive one-parameter progression.

### Proof

Write

```text
u_m = |m-alpha|,    v_m = |m-beta|.
```

Call `m` regular when both distances are at least `1/2`.  At a regular point,

```text
1 / (max(1,A u_m) max(1,B v_m)) <= 1/(A B u_m v_m).  (2)
```

We also use the elementary inequality, valid for every `z >= 0`,

```text
1/max(1,z/2) <= 3/(1+z).                              (3)
```

Indeed, for `z <= 2` the right side is at least one, while for `z >= 2` the
claim is `2/z <= 3/(1+z)`.

### Small separation: `0 <= Delta <= 1`

There are at most two nonregular integers.  Since

```text
Delta <= u_m + v_m,
```

at least one of the two distances is at least `Delta/2`.  Applying (3) to the
corresponding factor and summing over the at most two exceptional integers gives

```text
exceptional contribution <= 6/(1+A Delta) + 6/(1+B Delta).   (4)
```

For regular integers with `u_m < 2`, there are at most five choices and each
`1/(u_m v_m)` is at most four, for a contribution at most `20`.  If `u_m >= 2`,
the reverse triangle inequality gives

```text
v_m >= u_m - Delta >= u_m/2,
```

so `1/(u_m v_m) <= 2/u_m^2`.  Unit shells contain at most four relevant lattice
points and `sum_(k>=2) 1/k^2 <= 1`, so this part is at most `8`.  Thus the regular
reciprocal sum is at most `28`.  Since

```text
(1+log(2+Delta))/(1+Delta) >= 1/2,
```

the first term of (1), together with (2), covers the regular contribution.

### Large separation: `Delta >= 1`

Assume `alpha <= beta`.  The two exceptional half-unit neighborhoods are disjoint
and contain at most one integer apiece.  Near `alpha` one has
`v_m >= Delta/2`; near `beta` one has `u_m >= Delta/2`.  By (3), their total is at
most

```text
3/(1+A Delta) + 3/(1+B Delta).                        (5)
```

For regular integers between the poles, `u_m+v_m=Delta`, hence

```text
1/(u_m v_m) = (1/u_m + 1/v_m)/Delta.
```

The distances form unit-spaced sequences beginning at least at `1/2`, and harmonic
sum comparison gives

```text
sum_between 1/(u_m v_m) <= (4+2 log(2 Delta))/Delta.  (6)
```

On either exterior ray the distances have the form `r+k` and `r+k+Delta`, with
`r >= 1/2`.  Monotone integral comparison gives

```text
sum_(k>=0) 1/((r+k)(r+k+Delta))
  <= 1/(r(r+Delta))
     + integral_0^infty dx/((r+x)(r+x+Delta))
  <= (2+log(1+2 Delta))/Delta.                        (7)
```

Combining (6) and the two copies of (7), the regular reciprocal sum is at most

```text
(8+2 log(2 Delta)+2 log(1+2 Delta))/Delta
  <= 16(1+log(2+Delta))/(1+Delta).                    (8)
```

Equations (2), (5), and (8), with the harmless relaxation `16 -> 64` and `3 -> 6`,
prove (1).  QED.

## II. Exact refutation of the former lemma and sharpness of the logarithm

Take

```text
A=B=1, alpha=0, beta=N, N=2^22.
```

The middle terms alone give

```text
S(1,1;0,N) >= sum_(m=1)^(N-1) 1/(m(N-m))
            = 2 H_(N-1)/N
            >= 2 log(N)/N
            = 44 log(2)/N.                            (9)
```

The old claimed `C0=8` right side is

```text
8(3+log(2))/(N+1).                                    (10)
```

The elementary tangent bound

```text
1/x >= 4/3 - 4x/9   on [1,2]
```

follows from `(2x-3)^2/(9x) >= 0`; integration gives `log(2) >= 2/3`.
After cross multiplication, (9) minus (10) has the sign of

```text
(36 log(2)-24)N + 44 log(2),
```

which is strictly positive.  Thus the old lemma is false.

In fact the full sum in this example is

```text
S(1,1;0,N) = 4 H_(N-1)/N + 2/N + 2/N^2
           ~ 4 log(N)/N.                              (11)
```

Therefore the `log(2+Delta)/(A B(1+Delta))` behavior in (1) is necessary, not a
proof artifact.

## III. What the corrected lemma does and does not give for `T3`

For positive speeds `a,b,c`, put

```text
g  = gcd(a,b),          g0 = gcd(g,c),
d  = g/g0,              a' = a/g,
b' = b/g,               c' = c/g0.
```

Solvability forces `h3=d t`, and the reduced slice is

```text
a' h1 + b' h2 = -c' t.
```

If `(x_t,y_t)` is one solution, all solutions are

```text
(x_t+m b', y_t-m a'),
```

and the gauge-invariant pole separation is

```text
Delta_t = c'|t|/(a'b').                               (12)
```

Changing the particular solution translates both pole centers by the same integer,
so it neither changes (12) nor the sum.  The distinction between the original outer
coordinate `h3=d t` and the reduced parameter `t` is mandatory.

The main term of (1) becomes

```text
64(1+log(2+c'|t|/(a'b'))) / (a'b' + c'|t|).           (13)
```

After the outer `1/(d|t|)` weight, the range `|t|>H` in (13) has the expected size

```text
O((1+log(2+Vmax H))/H).
```

Under `H`-dissociation, primitive pair relations imply `H<Vmax`; thus this part is
`O((1+log Vmax)/H)`.  The corrected logarithm does not by itself raise the expected
final exponent: after the small-`t` harmonic sum, the plausible target remains

```text
T3(H) = O((1+log Vmax)^2/H).                           (14)
```

But (14) is not proved here.  The two near-pole terms in (1) are too coarse to sum
pointwise.  In the full-support relation mass the exact pole points correspond to a
zero coordinate and must be deleted; the neighboring points occur only in specific
integer congruence classes.  A proof still needs a **punctured near-pole congruence
lemma** which sums those classes globally.  The old phrase “`1/|t|`-compatible forms”
did not supply this lemma.  The committed numerical program computes only positive
box truncations, hence lower approximations to the infinite absolute tail; it cannot
certify an upper bound without an omitted-tail estimate.

## IV. Why `T4` and `T5` are not mechanical nesting

At support four, after fixing two outer coordinates `(u,t)`, the innermost two-pole
separation is, up to primitive gcd reduction,

```text
Delta_(u,t) = |c u + d t|/(a b).                      (15)
```

It does not grow with `|t|` uniformly: it vanishes on arbitrarily large points of the
resonance line `c u+d t=0`.  The missing result is a uniform affine-coset estimate
which sums both the resonance strip around that line and its complement with explicit
constants.

At support five the corresponding obstruction is a resonance slab around

```text
c u + d t + e r = 0.                                  (16)
```

This is a rank-two plane of near-cancellation, not another copy of the homogeneous
`T3` statement.  The required strip and slab estimates, plus their gcd, zero-coordinate,
and constant bookkeeping, remain open.

## V. The honest conditional exhaustion statement

Let

```text
L = 1+log(2+Vmax).
```

Assume one eventually proves, for every relevant support,

```text
Ts(H) <= Cs L^(s-1)/H,    s=3,4,5.                    (17)
```

Then THM-935 gives

```text
|B5-2052/16807|
 <= [312/343 + (6864/49)C3 L^2
              + (1430/7)C4 L^3
              + 1287 C5 L^4] / H.                   (18)
```

Consequently, if `H >= K L^4` and

```text
K > (16807/2052)
      (312/343 + (6864/49)C3 + (1430/7)C4 + 1287 C5),
```

then absence of a support-at-most-five relation of height at most `H` implies
`B5>0`.  This would prove only the analytic alternative

```text
small support relation  OR  B5 > 0.                   (19)
```

The first side of (19) is not a completed LRC classification.  Ratio, cascade,
covering, and near-AP work provide routes for parts of it, but the structured branch
remains an explicit proof obligation.  Thus neither (17) nor even its future proof
would by itself close LRC(14).

## VI. Tournament analysis and challenged vertex assumption

Runner vertices are not the natural quotient for this tail problem.  The preserved
object is the affine resonance arrangement of slice cosets and pole windows.  The
pair observable is the signed order of the two pole centers; changing a particular
solution shifts both centers by the same integer, so this edge and its absolute
separation are gauge invariant.  With only two poles this is a single edge and has no
nontrivial score histogram, cycle, SCC, or Hamiltonian-path fingerprint.

For (15)--(16), the honest higher object is closer to an **oriented matroid**: sign
cells and circuits record resonance lines and planes.  Treating runners as tournament
vertices destroys the affine offsets that control the estimate.  A future elimination-
order tournament may be useful for choosing which coordinate pair to leave innermost,
but no canonical antisymmetric observable or LRC-preservation theorem is presently
known.  Tournament Analysis is therefore not forced into the current scalar script;
the challenged assumption and the alternate coset/pole-window vertices are recorded
explicitly.

## Exact remaining tasks

1. Prove the punctured near-pole congruence sum needed for `T3`, with constants.
2. Prove the affine resonance-strip estimate for `T4`.
3. Prove the affine resonance-slab estimate for `T5`.
4. Insert the resulting `C3,C4,C5` into (18).
5. Separately solve the structured small-relation branch of (19).
