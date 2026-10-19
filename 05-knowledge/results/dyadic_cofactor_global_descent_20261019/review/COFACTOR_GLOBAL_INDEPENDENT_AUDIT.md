# Independent audit of the dyadic cofactor-map classification

Date: 19 October 2026.

## Verdict

**Pass. No mathematical gap found.** The proof establishes the stated global
classification for the auxiliary map on all positive even integers. Its
odd-prime-guarded aliquot consequences follow. It makes no claim about aliquot
orbits after leaving that domain, and no general aliquot or Collatz convergence
conclusion is justified or needed.

Audited source:

`COFACTOR_MAP_GLOBAL_CLASSIFICATION.md` (reviewed before its status-only publication update; see `SOURCE_PROVENANCE.json` for both hashes)

SHA-256:

`1572502c1ca5151562d91a19b36b662eabf27336477472c1f245f14d6659bab8`

## Independent proof review

Write `n=cp`, with `c` its exact power-of-two factor and `p` positive odd, and
put `q=p+1`. The defining expression is positive and even, and subtraction gives
`F(n)-n=2c-q`. Thus immediate descent, fixed points, and abundance are correctly
partitioned by `q>2c`, `q=2c`, and `q<2c`.

### Non-power-of-two q

If `q=dk` with `d=2^v2(q)` and `k>=3` odd, abundance forces `d<c`:
otherwise `q>=3c>2c`. Both factors are powers of two, so `c/d` is even.
Consequently

    F(n)/d = (c/d)(q+1)-k

is odd; the claimed cofactor is exact. Applying the next step with that exact
cofactor yields

    d(F²(n)-n) = (c+d)(2d-q-1)+q.

The estimate `2d-q-1<=-d-1`, followed by `q<2c`, gives a strict negative upper
bound. No primality assumption is used here.

### Power-of-two starts

For `C>=8` a power of two, the first image `3C-2` is exactly divisible by 2,
as is the second image `3C/2+2`. This proves all three consecutive formulas,
not merely formal affine substitutions. The third image is below `C` exactly
when `C>16`, so `C>=32` covers the stated nonexceptional power starts.

The restriction `C>=8` matters: at `C=4` the second image has cofactor 8,
not 2. The candidate does not apply the formula there.

### q=d a power of two, with p>=3

Here `d>=4`, and abundance implies `c>=d`. These cases are exhaustive:

1. `c=d`;
2. `d<c<d²`;
3. `c=d²`;
4. `c>d²`.

For the strict middle case, set `e=c/d`. Then `2<=e<=d/2`, and both `e`
and `d/e` are powers of two with `d/e` even. The exact cofactors at the first
and second images are respectively `d` and `e`. Independent expansion yields

    F³(n)-n = -d(d-e)+e+d-d/e
              <= -d²/2+3d/2-2 < 0.

The final bound is already `-4` at `d=4` and decreases thereafter.

For `c>d²`, let `h=c/d²>=2`. The first two image cofactors are both exactly
`d`, because their displayed quotients are odd. The third-step difference is

    -h(d-1)+3d-2.

For `h>=4` this is at most `-d+2<0`. Powers of two leave only `h=2`.
In that case the third image is `n+d` with exact cofactor `d`; the genuine
fourth-step difference is

    -(2d-1)(d-2)<0.

### Resonances

For `c=d`, the first image is exactly `d²`. For `d>=8`, the power-start block
therefore gives

    F⁴(n)-n = -d²/4+d+4 <= -4.

For `c=d²`, the first image has exact cofactor `d`, and the second image is
exactly `d³`. For `d>=8`, the same block gives

    F⁵(n)-n = -d³/4+d²+4 <= -60.

The remaining `d=4` cases are correctly separated: `n=12` is exceptional;
`n=48` has the genuine chain

    48, 76, 64, 190, 98, 52, 46.

All first five images of 48 are larger than 48. Therefore six is a genuinely
necessary uniform bound for this descent statement, not an artifact of a
loose estimate.

### Global argument

The starting-value exceptions `{2,4,8,12,16}` and their displayed trajectories
are correct. Every nonfixed, nonexceptional block start has a strictly smaller
actual iterate after at most six steps. Repeating such blocks produces strictly
decreasing positive even block starts. Well-ordering forces termination at a
fixed point or declared exception. The finite exceptional trajectories then
give exactly the two stated nontrivial cycles.

The argument does not assume boundedness, eventual periodicity, or a prior
classification of cycles. It proves those properties. Nor does it require
every individual step to descend.

### Native aliquot domain

For `p` an odd prime, multiplicativity gives

    sigma(cp)-cp = (2c-1)(p+1)-cp = F(cp).

Agreement holds at each retained-domain input. Both nontrivial auxiliary-map
cycles contain a power of two, with odd part 1. Thus no aliquot orbit remaining
in the odd-prime class forever can enter either cycle. Its auxiliary fixed
point has `p=2c-1` prime, hence is an even perfect number.

This also excludes every nontrivial aliquot cycle contained entirely in the
class, with no restriction on cycle length. It does not exclude aliquot cycles
that leave this class. In particular, `F(120)=120` while the true aliquot sum of
120 is 240; the composite guard is essential.

## Independent executable controls

Script: `audit_cofactor_global.py`.

Recorded output: `cofactor_global_audit_results.json`.

The script uses direct integer evaluation of `F` independently of the proof's
branch certificates. All arithmetic is exact.

- Sixteen independently reconstructed symbolic identities simplify to zero.
- Exhaustive test of all 2,000,000 positive even starts through 4,000,000:
  every claimed branch-specific endpoint descends, except the ten fixed
  points and five stated exceptions.
- Independently measured first-descent counts: 1,996,001 at step 1; 3,875 at
  step 2; 90 at step 3; 13 at step 4; five at step 5; one at step 6.
  The sole six-step start in this range is 48.
- 150,558 additional structured and large-integer checks. These include
  pairs of dyadic exponents through 180, explicit valuation assertions in
  every transition regime, power starts through exponent 600, and 20,000
  reproducibly seeded arbitrary odd-part checks with cofactor exponents and
  odd-part bit lengths through 2,048.
- A separate actual-orbit cycle scan, without using descent cases, from every
  positive even start through 500,000 visits 250,011 distinct integers and
  finds only the two asserted nontrivial cycles and the encountered fixed
  points. Largest visited integer: 786,430.
- Independent divisor enumeration agrees with `F` on all 3,869 odd-prime-
  guarded starts through 30,000. Explicit negative controls distinguish
  `F` from the aliquot map at 8, 16, and 120.

These finite checks are supplementary controls. The complete result rests on
the exhaustive symbolic case split and well-ordering argument reviewed above.

## Scope and presentation

No required repair was identified. The source appropriately distinguishes the
auxiliary map from the genuine aliquot map and expressly leaves historical
priority unestablished. This audit does not establish novelty or priority.
