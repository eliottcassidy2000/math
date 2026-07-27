---
id: HYP-9033
title: "The discriminant tower: Jelonek divisors as classification data for the Keller monoid, and the genus axis of collapse vs rigidity"
status: >
  OPEN synthesis with a PROVED fixed-map core: the corrected saturated
  cuspidal law 27 c^2 L = S^2 - T^3 (MISTAKE-287) and THM-2570's global
  normalization/conductor theorem.  THM-2576 proves the set-level composition
  law and the two irreducible components of S_(F o F); its discriminant
  multiplicities remain open.  The synthesis retains several VERIFIED inputs,
  four wider predictions, and a typed dictionary entry for HYP-9031.  Nothing
  here is a JC(2) or classification claim.
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
  - THM-2570 (global normalization, conductor, and cusp-cylinder structure)
  - THM-2576 (image divisor and two-component nonproperness law for F o F)
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

## 2. The discriminant tower [set law and grade two PROVED; multiplicities OPEN]

THM-2576 proves for dominant polynomial maps, directly from escape sequences,

```text
S_(F_out o F_in)=S_(F_out) union F_out(S_(F_in)).
```

For the fixed sporadic map it computes the irreducible degree-25 polynomial
`H` cutting `closure(F(V(L)))` and proves

```text
S_(F o F)=V(LH),
```

with exactly the two distinct irreducible components `V(L)` and `V(H)`.
The raw inverse-coordinate resultant is `-a^8 c^18 S^8 H`; the first three
factors are chart artifacts and say nothing about composite multiplicity.

The remaining P1 prediction is that the odd part of the composite
escape-coordinate discriminant is `L*H` times a square, and that its integral
coordinate leads are powers of `det J_(F o F)=4`.  At higher grades the exact
set law gives the image tower, but distinctness and irreducibility of its
successive closures remain open.

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

For this fixed map, THM-2570 now makes the genus-0 statement exact rather
than metaphorical: on `c!=0` the whole Jelonek component has coordinate ring
`Q[c,c^-1,theta^2,theta^3]`, its normalization is
`Q[c,c^-1,theta]`, and the missing exponent `theta^1` is precisely the
conductor defect supported on the empty-fibre curve.  The global normalizing
coordinate `lambda=(2-theta)/c` continues through the smooth `c=0` boundary.
This verifies the proposed axis for the sporadic map only; it does not prove
that genus controls a general Keller map or any remaining planar stratum.

Prediction P5: remaining JC(2)-adjacent closures (higher mixed strata)
will continue to run through genus >= 2 of restored discriminant-side
curves; a stratum whose restored curve is forced rational would be the
place to HUNT for a JC(2) counterexample instead.

## 5. Predictions ledger

- P1 (set part proved at grade two): `S_(F o F)=V(LH)` with two irreducible
  components (THM-2576).  The disc odd part `L*H` times a square remains open.
- P2 (in flight): W1 obeys `disc = -(det J)^2 (square)^2 L_W1` with
  `L_W1 = L o T1^(-1)`; leads `{|det J|, |det J|^3}` after suitable
  normalization; W2 likewise in its own normal form if explicit.
- P3: any new atom candidate must simultaneously break monodromy
  imprimitivity AND Jelonek reducibility expectations.
- P4: component count of the Jelonek divisor = the grade exponent `k`,
  generically.
- P5: the genus axis governs the remaining JC(2) strata.

## Scorecard (same session, keller_disc_tower_opus_20260728)

P1 RESOLVED, half-refuted by a BETTER law: L_2 = equation of F(Z(L))
computed and PROVED (irreducible, degree 25, degrees (14,21,12), 361
terms); S_(F o F) = Z(L) u Z(L_2) CONFIRMED with exact fiber-drop
witnesses (9/7/3/0). But the disc odd part at level 2 is L_2 ALONE --
L appears to EVEN multiplicity 8: the middle level's -L class cancels
against the norm law prod_(fiber) L = L_2/(64 L) (a conductor-
discriminant mechanism). Corrected tower law: **odd part of disc at
level k = L_k alone; the level-k escape-coordinate lead IS L_k (up to
64L normalization); the depressed (trace-zero) structure persists at
depth 2** (xi^8 coefficient identically zero -- induction candidate
for all F^k). P2 CONFIRMED on BOTH tame conjugates W1 AND W2:
disc = -4 (square)^2 L_W with L_W = L o T^(-1) (full target
covariance) and lead pattern {L_W, -|det J|, |det J|^3} -- the law is
a gauge-covariant of the ATOM. The composite-lead prediction is
REFUTED as posed and replaced: integral-coordinate leads iterate by
CUBING (8 -> 512 = 2^9; lead_(k+1) = lead_k^3), tied to |det J_F| = 2
per-atom, not to det of the composite. P4 consistent (2 components at
grade 9). Next: level 3 (is the degree-27 disc odd part L_3 alone,
with L_3 the norm-transfer of L_2?).

## Loss ledger

The saturated cuspidal law is proved for the sporadic `F` only; the raw
pullback has the extra affine plane `c=0` described in (C).  The set-level
composition law is general and the two grade-two components are proved, but
the discriminant multiplicities, `-(det J)^2` law, and higher component count
remain predictions.  "Generic component count" needs a precise genericity
notion before any classification claim; none of this closes JC(2), classifies
the monoid, or excludes G1 -- it supplies invariants and a hunting map.
