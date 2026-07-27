---
id: HYP-9033
title: "The discriminant tower: Jelonek divisors as classification data for the Keller monoid, and the genus axis of collapse vs rigidity"
status: >
  OPEN synthesis with one PROVED identity (the corrected saturated cuspidal
  law 27 c^2 L = S^2 - T^3, exact; MISTAKE-287), several VERIFIED
  inputs (THM-2473/2546 and referee), five falsifiable predictions
  (P1-P2 under computation in flight), and a typed dictionary entry
  for HYP-9031.  Nothing here is a JC(2) or classification claim.
source: opus-2026-07-28 (from the uniform disc law found by the
  THM-2546 referee run: disc = -4 (square)^2 L for all three
  coordinate cubics of the sporadic Keller map)
related:
  - HYP-9030-keller-degree-semigroup (the monoid; atoms {1,3} conjectured)
  - THM-2473 / THM-2546 (branch tower; integral-coordinate dichotomy)
  - THM-1330 (necessary monoid atlas -- gains an invariant here)
  - HYP-9031 (the D5 dictionary -- gains the genus axis)
  - THM-2465 (W1/W2 members; the -(det J)^2 law test target)
  - THM-2566 (two-chart saturated cusp atlas and parasitic-plane ledger)
---

# HYP-9033 -- the discriminant tower and the genus axis

## 1. The saturated cuspidal law [PROVED, exact; corrected by MISTAKE-287]

For the sporadic Keller map `F` (THM-1300/2473), with
`T = 4 - 3bc`, `S = 27ac^2 - 9bc + 8`:

```text
27 c^2 L = S^2 - T^3        (exact polynomial identity)
```

so the Jelonek quartic `L` is the ELLIPTIC DISCRIMINANT (difference of
square and cube) of the trisection pencil `L x^3 + T x - 2c`.  However,
the **raw global pullback** of the cuspidal cubic under
`Phi_c:(a,b,c) -> (S,T)` is

```text
Phi_c^* V(S^2-T^3) = V(c^2L),
underlying set       = V(c) union V(L),
effective divisor    = 2V(c) + V(L).                (C)
```

Thus `Z(L)` is exactly the pullback only on `c!=0`; globally it is the
`c`-saturated pullback.  The lost factor is geometrically real:
`Phi_c` collapses the honest affine target plane `c=0` to the smooth cusp
point `(8,4)`.  For example `(a,b,c)=(1,0,0)` maps to `(8,4)` but has
`L=16`, so it is not on the Jelonek set.  This corrects the former global
pullback sentence (MISTAKE-287); THM-2566 gives the exact two-chart atlas.

On the localized chart `c!=0`, the stratification matches THM-2473 exactly:

- off `Z(L)`: full fibers (3) -- off the cusp curve on this chart;
- `Z(L) \ E`: drop to 1, survivor `x = 2c/T` = the LINEAR REMNANT of
  the degenerated pencil -- smooth points of the cusp curve;
- the empty-fiber curve `E` (where `T = 0`, hence `S = 0`): the
  preimage of THE CUSP POINT `(S,T) = (0,0)` itself -- total escape.

The last assertion is global, not merely localized: `T=S=0` forces
`bc=4/3` and `a=4/(27c^2)`, so `c` is automatically nonzero and the
cusp-point preimage is exactly `E`.

The trisection anatomy `4T^3 - 3T -+ 1 = (T -+ 1)(2T +- 1)^2`
(HYP-9030's branch identity) is this same cuspidal geometry at the
Chebyshev point. The uniform coordinate law (THM-2546 + referee):
`disc_xi = -4 (square)^2 L` for ALL THREE coordinates, with
`-4 = -(det J)^2` and integral-coordinate leads `{2, 8} = {|det J|,
|det J|^3}` -- the whole eliminant package is `det J`-graded.

## 2. The discriminant tower [prediction P1, computation in flight]

For a composite `G = F_out o F_in`: `S_G = S_(F_out) u F_out(S_(F_in))`.
Prediction: the odd part of `disc` of `G`'s escape-coordinate eliminant
factors as the product of the tower's Jelonek factors
(`L * L_2` for `F o F`, with `L_2` cutting `F(Z(L))`), each to odd
multiplicity, times projection-collision squares; and the composite's
integral-coordinate leads are powers of `det J_G = 4`. **The monoid
grading of HYP-9030 is mirrored by a divisor tower: each composition
appends one Jelonek factor.**

## 3. Two crevasse invariants for atomhood

The atom-degree floor program (HYP-9030's Busch-crevasse analogue)
gains a second invariant:

1. **Monodromy primitivity** (THM-2473: `F o F` is maximally
   imprimitive; a primitive composite-degree monodromy = new atom).
2. **Jelonek irreducibility**: `F`'s `L` is irreducible (atom-
   consistent); a `k`-fold composite generically has `>= k` Jelonek
   components. Prediction P4: grade-`3^k` members have Jelonek sets
   with generically exactly `k` components; an irreducible Jelonek set
   at composite grade would be strong atom evidence (one-way test:
   components can coincide).

Classification data proposal for THM-1330's atlas: the triple
`(fiber degree 3^k; Jelonek divisor with its component tower, mod
Aut(target); monodromy subgroup of the iterated wreath with its block
system)`. `L` is invariant under source-tame conjugation and covariant
under target automorphisms, so the divisor class is a genuine invariant
of the monoid element modulo the tame moves that HYP-9030 says act
within a grade.

## 4. The genus axis [new HYP-9031 dictionary entry]

The cuspidal cubic is RATIONAL (genus 0) -- and its rational
parametrization is realized by the map's own exceptional curve
`t -> (4/(27t^2), 4/(3t), t)` sitting over the cusp. **The JC >= 3
counterexample family exists because its discriminant geometry is
genus 0: escape/collision can walk along a rational curve.** By
contrast, every degree-18/22 JC(2)-side plane closure succeeded by
restoring a root coordinate and forcing GENUS >= 2 on the resulting
curve (THM-2463/2468/2469/2470/2472/2475/2476/2480): rigidity exactly
where the discriminant-side geometry has no rational walk.

> **Genus axis (dictionary row for HYP-9031):** collapse happens where
> the branch/discriminant geometry is rational; rigidity is proved
> where it has genus >= 2. JC(3) fell on a genus-0 cusp; JC(2)'s
> closures manufacture genus >= 2; LRC(14)'s spectral side has no
> curve at all (Diophantine rigidity) -- the three frontiers sit at
> genus 0 / genus >= 2 / no-geometry on one axis.

Prediction P5: remaining JC(2)-adjacent closures (higher mixed strata)
will continue to run through genus >= 2 of restored discriminant-side
curves; a stratum whose restored curve is forced rational would be the
place to HUNT for a JC(2) counterexample instead.

## 5. Predictions ledger

- P1 (in flight): `S_(F o F) = Z(L) u Z(L_2)` and the disc odd part is
  `L * L_2` on slices.
- P2 (in flight): W1 obeys `disc = -(det J)^2 (square)^2 L_W1` with
  `L_W1 = L o T1^(-1)`; leads `{|det J|, |det J|^3}` after suitable
  normalization; W2 likewise in its own normal form if explicit.
- P3: any new atom candidate must simultaneously break monodromy
  imprimitivity AND Jelonek reducibility expectations.
- P4: component count of the Jelonek divisor = the grade exponent `k`,
  generically.
- P5: the genus axis governs the remaining JC(2) strata.

## Loss ledger

The saturated cuspidal law is proved for the sporadic `F` only; the raw
pullback has the extra affine plane `c=0` described in (C).  The tower and
`-(det J)^2` laws are predictions under test; "generic component
count" needs a precise genericity notion before any classification
claim; none of this closes JC(2), classifies the monoid, or excludes
G1 -- it supplies invariants and a hunting map.
