---
id: THM-2546
title: "Integral-coordinate dichotomy of the sporadic Keller map and the scope of the parity lens"
status: >
  PROVED for (1) the y-integrality (exact symbolic y-eliminant
  2y^3 - 3b y^2 + 18a y + (27a^2c - 18ab + b^3), constant lead 2), (2)
  the gauge-hostile Tr(y) = 3b/2 killing the naive stratum-wide
  depressed-law conjecture, and (3) the pointwise parity rigidity
  lemma; VERIFIED for the z-integrality (five integer-target exact
  z-eliminant cores, leads in {1,2,4,8}, exact symbolic lead OPEN).
  RELEASES HYP-9030 test (ii) as stated (pointwise odd-real-count
  witnesses cannot exist for even-degree wild maps) and records the
  redirected parity program.  Nothing here excludes G1 or any stratum.
source: opus-2026-07-27
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
related:
  - HYP-9030-keller-degree-semigroup (test (ii))
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
script: 04-computation/keller_integral_coordinate_parity_opus_20260727.py
output: 05-knowledge/results/keller_integral_coordinate_parity_opus_20260727.out
---

# THM-2546 -- one coordinate carries all the escape

## (1) The integral-coordinate dichotomy [y PROVED, z VERIFIED]

For the sporadic Keller map `F` (THM-1300/2473) with target `(a,b,c)`:

```text
y-eliminant core:  2 y^3 - 3b y^2 + 18a y + (27 a^2 c - 18 a b + b^3),
```

with CONSTANT leading coefficient `2 = |det J_F|`: the coordinate `y`
is integral over `Q[a,b,c]` -- no fiber sequence ever escapes in `y`.
The z-eliminant cores at integer targets have leads in `{1,2,4,8}`
(2-adic units up to content) and 2-power denominators in the symmetric
functions: `z` is integral as well (VERIFIED; symbolic lead OPEN).
Since the x-eliminant lead is the nonconstant irreducible
`L = 27a^2c^2 - 18abc + 16a + b^3c - b^2` (THM-2473(2)):

> **Monogenic escape law.** All Jelonek escape of `F` is carried by the
> single non-integral coordinate `x`; every escaping fiber point has
> `|x| -> infinity` with `y, z` bounded, and the even escape parity of
> THM-2473(3) governs the entire boundary behavior.

## (2) The depressed law is an x-gauge fact [PROVED hostile]

`Tr(x) = 0` identically (THM-2473(2)) but `Tr(y) = 3b/2 != 0`: the
"depressed trisection shape" is a property of the escape coordinate,
not of the map. This kills the naive conjecture (raised in this
session's point-cap analysis) that the graph polynomial of a point-cap
Keller map is depressed in every (or any a-priori-chosen) primitive
coordinate; any depressed law must FIRST select the escape coordinate.

## (3) Pointwise parity rigidity and the release of HYP-9030 test (ii)

**Lemma (trivial, load-bearing).** Let `G` be a real polynomial map
`R^n -> R^n` with everywhere-nonzero constant Jacobian, of complex
degree `d`. At every value whose complex fiber is FULL (`d` points),
complex conjugation is an involution of the fiber fixing exactly the
real points, so `#real = d (mod 2)`.

Consequences: (i) an even-degree wild map has an EVEN real count at
every full-fiber value -- a pointwise odd-real-count witness CANNOT
exist, so HYP-9030 test (ii) ("does det J < 0 force an odd-real-count
value on the point-cap stratum") is structurally vacuous as stated:
at regular values it is *equivalent to* the conjecture it was meant to
attack, and odd counts occur only on escape strata where the
congruence fails (witness: THM-2473's empty-fiber curve gives `0` real
preimages for the odd-degree `F`). **Test (ii) is RELEASED as stated.**

**The redirected parity program** (what survives, in decreasing
concreteness): (a) escape-parity laws per non-integral coordinate --
for `F` the escape is even BECAUSE the escape coordinate's eliminant is
depressed; for a hypothetical degree-4 point-cap witness the analogous
question is the subleading coefficient of its escape-coordinate
eliminant on `Z(lead)`; (b) the sigma-fixed-locus mechanism (S128c97):
an equivariance with odd fixed-fiber forces odd counts; (c) integrated
sphere degree: at any value `v` over which the map is proper, the
normalized boundary map `S^2_R -> S^2` has degree `-k(v)` for large
`R`, so a leading-form computation forcing odd degree at one proper
value would contradict (i) -- this is the only pointwise-flavored
route left, and it must go through properness-over-a-value, not
regularity.

## Loss ledger

The dichotomy is proved for the sporadic `F` only; no claim is made
for general Keller maps or the point-cap stratum. The z-integrality is
sampled (5 exact integer-target eliminants), not symbolic. The release
of test (ii) does not weaken HYP-9030's atom conjecture; it removes
one of its three proposed tests and replaces it with the redirected
program above.
