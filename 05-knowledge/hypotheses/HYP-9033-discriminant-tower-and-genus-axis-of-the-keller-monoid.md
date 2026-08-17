---
id: HYP-9033
title: "The discriminant tower: Jelonek divisors as classification data for the Keller monoid, and the genus axis of collapse vs rigidity"
status: >
  OPEN synthesis with a substantially PROVED fixed-map core.  The saturated
  cusp law and normalization are proved; THM-2582/3495/3498/3525/3526 give
  the first six literal-coordinate discriminant square classes.
  THM-3528 proves raw complete packets at every level, THM-3529 proves every
  complete packet is a finite-sheet unit, and THM-3530 proves every raw rung
  is an absolute image prime with generic image degree one and exactly n
  reduced Jelonek components for F^n.  THM-3531 proves the intrinsic
  all-level class `[(-1)^nP_(n-1)]`; THM-3533 gives newest-prime
  normalization exponent one; and THM-3535 proves full wreath monodromy plus
  every nonzero constant linear form primitive at all depths.  THM-3532 proves
  full covariance on the fixed polynomial-conjugacy orbit and only one-step
  covariance for independent source/target changes.  Constant-form order
  indices, old-prime multiplicities, wider tame-equivalence covariance, and
  monoid-wide classification remain open.
  Nothing here is a JC(2) or arbitrary-map classification claim.
source: opus-2026-07-28 (from the uniform disc law found by the
  THM-2546 referee run: disc = -4 (square)^2 L for all three
  coordinate cubics of the sporadic Keller map)
related:
  - THM-3438-weighted-lift-keller-degree-spectrum (S_n atoms in every grade)
  - HYP-9030-keller-degree-semigroup (old atom spectrum refuted)
  - THM-2473 / THM-2546 (branch tower; integral-coordinate dichotomy)
  - THM-1330 (necessary monoid atlas -- gains an invariant here)
  - HYP-9031 (the D5 dictionary -- gains the genus axis)
  - THM-2465 (W1/W2 members; the -(det J)^2 law test target)
  - THM-2566 (two-chart saturated cusp atlas and parasitic-plane ledger)
  - THM-2570 (global normalization, conductor, and cusp-cylinder structure)
  - THM-2576 (image divisor and two-component nonproperness law for F o F)
  - THM-3495 (third prime image, three components, and degree-27 square class)
  - THM-3498 (fourth old-boundary valuation and degree-81 square-class gate)
  - THM-3504 (fourth prime image and four-component nonproperness set)
  - THM-3506 (five-face transform and the next old-boundary valuation)
  - THM-3521 (fixed-R5 finite sheet and following old-boundary valuation)
  - THM-3522 (fixed-chart complete-packet renewal propagation)
  - THM-3528 (all-level raw complete packets and defect identity)
  - THM-3529 (all complete packets are finite-sheet units)
  - THM-3530 (all-level fixed raw image primes and component count)
  - THM-3531 (intrinsic all-level discriminant square class)
  - THM-3532 (two-sided one-step and fixed conjugacy-orbit covariance)
  - THM-3533 (newest-prime normalization exponent one and index-square law)
  - THM-3535 (full wreath monodromy and all-level constant-linear primitivity)
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

Over the generic target field this has the precise square-class consequence

```text
[disc_x]=[disc_r]=[disc_z]=[-L].
```

Thus the three cubics share one quadratic sign-resolvent/Kummer extension;
they do not supply three independent binary classes.  This synchronizes only
the sign quotient of their root monodromy, not their cubic splitting fields,
labelled roots, or boundary-effective sections.  The exact scope and the
infinite-family classification boundary are separated in
`07-reflections/jc-three-cubics-one-kummer-class-and-family-classification-boundary-codex-20260815.md`.

## 2. The discriminant tower [set law and fixed grades two through four PROVED]

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
With the standard Sylvester/root-product convention corrected by MISTAKE-290,
the raw inverse-coordinate resultant is `a^8 c^18 S^8 H`; the first three
factors are chart artifacts and say nothing about composite multiplicity.

The concurrent exact inverse-section computation
`keller_disc_tower_opus_20260728.py` independently reconstructs the same `H`
(identical 361-term coefficient-ledger hash).  THM-2582 proves the resulting
function-field norm identity

```text
product_(q in F^-1(t)) L(q) = H(t)/(64L(t)).
```

On three exact rational `(b,c)` slices, the degree-nine escape-coordinate
discriminant contains `H` to odd multiplicity one and `L` to even multiplicity
eight.  Thus the former prediction “odd part `L*H`” is **REFUTED on all three
slices**.  THM-2582 proves the global survivor: the odd square-free divisor is
`H` alone.  Its all-degree block lemma shows that the cross-block resultant is
alternating precisely for odd block degree and contributes the outer
discriminant class; together with `Norm(L)=H/(64L)`, this cancels `L` in
parity.  This determines square class, not exact positive multiplicity.
Likewise the former “powers of
`det J_(F o F)=4`” lead prediction is refuted for this normalization; the
observed canonical integral leads are `8` and `512`, obtained by cubing the
atom-level leads.

THM-3495 proves the next fixed-map step.  With primitive absolutely
irreducible `J=2^35L^7N(H)`,

```text
closure(F(V(H)))=V(J),
S_(F^3)=V(LHJ),
[Delta_3]=[-2J].
```

The three factors are pairwise distinct and the restriction
`V(H)->V(J)` has generic degree one.

THM-3498 proves the next discriminant gate without constructing the next
image equation.  The exposed face

```text
in(J)=-2^58 3^51 13^8 79^4 313^2 x^43(3xz-2y)^15
```

and a finite-sheet unit control give `v_L(N(J))=-43`, hence
`G=L^43N(J)` lies in `Q[a,b,c]` and is coprime to `L`.  A squarefree
degree-81 good-reduction fibre then makes the odd-block recursion lawful and
gives

```text
[Delta_4]=[2G].
```

The constant class `[2]` is retained.  THM-3504 closes the image question
without a global expansion: on `b=c=1`, a proved degree bound `542` and an
exact squarefree `F_1009` specialization force the only possible norm
exponent to be one.  The same specialization is coprime to the old factors.
Consequently

```text
closure(F(V(J)))=V(G),
S_(F^4)=V(LHJ G),
gcd(G,LHJ)=1,
```

`G` is absolutely irreducible, and `V(J)->V(G)` has generic degree one.
THM-3506 then derives the hidden five-face transform directly from the inverse
chart.  The first untested pair is `(271,99)`, not `(259,87)`; a finite-sheet
unit gives

```text
v_L(N(G))=-271,
R_5=L^271N(G) in Q[a,b,c],
gcd(R_5,L)=1.
```

It also proves the exposed face
`R_5 ~ x^1699(3xz-2y)^615` in the relevant initial form.  At that theorem
state the fifth image role, integral normalizations, global degree ledgers,
and degree-243 separability were open.
THM-3521 closes the next finite-sheet question by proving

```text
v_L(N(R_5))=-1699,
R_6=L^1699N(R_5) in Q[a,b,c],
gcd(R_6,L)=1.
```

It also transports the next top face with pair `(10663,3867)`.  THM-3522
then proves the general fixed-chart implication

```text
complete A(e,m) + polynomial L^eN(P)
  => complete A(7e-2m,3e-2m),
```

and consequently gives complete packets `A(1699,615)` for `R_5` and
`A(10663,3867)` for `R_6`.  THM-3523/3527 pay the next two finite-sheet
gates, while THM-3525/3526 prove the degree-243/729 separability gates and
classes `[-2R_5]`, `[2R_6]`.  THM-3528 then proves raw polynomial complete
packets at every level.  THM-3529 identifies the finite divisor as `V(F^*L)`
and proves its order is zero on every complete packet.  Finally THM-3530
uses the primitive top-face exponents to force every norm-image exponent to
one, proving all raw rungs are absolute image primes and

```text
S_(F^n)=union_(j=0)^(n-1) V(P_j)
```

with exactly `n` reduced components.

THM-3531 now closes the intrinsic square-class recursion at every depth.  If
`K_n/K_0` is the degree-`3^n` function-field extension and `P_j` denotes the
raw rung, then

```text
d(K_n/K_0)=[(-1)^nP_(n-1)].
```

This is basis-independent and makes the newest raw prime the sole odd
horizontal discriminant divisor.  THM-3533 strengthens this at the newest
prime: its normalization-discriminant multiplicity is exactly one, while an
integral primitive power order has multiplicity `1+2i` for local index length
`i`.  THM-3535 supplies the separate coordinate gate: a supported
transposition forces full iterated wreath monodromy, and every nonzero constant
linear form is primitive at every depth.  Its local power-order index remains
open.  THM-3532 transports the raw packet, prime, intrinsic class, newest
different, Jelonek tower, and each primitive observation `ell` as
`ell o phi^(-1)` under honest conjugacy.  For affine `phi` these transported
observations remain affine-linear (the constant is irrelevant to generated
fields); for nonlinear `phi` they need not be standard linear forms.  A
general left/right tame equivalence transports one edge but inserts a middle
automorphism on iteration, so no broader family theorem follows.

## 3. Two crevasse invariants for atomhood (spectrum repaired)

THM-3438 supersedes the old atom-degree floor: weighted lifts give a primitive
`S_n` atom in every grade `n>=3`, and numerical product grades also contain
composites.  The two invariants below now distinguish maps *within* a mixed
grade rather than predict which atom degrees exist:

1. **Monodromy primitivity** (THM-2473: `F o F` is maximally
   imprimitive; THM-3438 realizes primitive `S_n` monodromy at every grade).
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

- P1 (fixed component, intrinsic square-class, and constant-linear towers
  proved all levels):
  `S_(F o F)=V(LH)` with two irreducible components (THM-2576), while the
  degree-nine eliminant has discriminant square class `[H]` (THM-2582).
  THM-3495 gives `S_(F^3)=V(LHJ)` and degree-27 class `[-2J]`.  Odd part
  `L*H` at grade two is refuted; exact positive multiplicities remain open.
  At grade four THM-3498 proves the square-class/localization gate
  `[Delta_4]=[2G]`, while THM-3504 proves `S_(F^4)=V(LHJ G)` with four
  pairwise-distinct prime components and generic degree one on the newest
  image restriction.  THM-3525/3526 add `[-2R_5]` and `[2R_6]` at depths
  five and six.  THM-3528/3529 give all-level complete packets and old-`L`
  units, and THM-3530 proves the newest raw polynomial is an absolute image
  prime of generic degree one at every depth.  THM-3531 proves the intrinsic
  class `[(-1)^nP_(n-1)]` at every depth; THM-3533 gives newest-prime
  normalization exponent one; and THM-3535 proves every nonzero constant
  linear form primitive there, so all corresponding coordinate eliminants
  have that class.  Their local indices and exact old-prime multiplicities
  remain open.
- P2 (PROVED fixed-orbit covariance; VERIFIED-EXACT controls): W1 and W2 are
  target postcompositions and obey the one-step law
  `disc=-(det J)^2(square)^2L_W`, with `L_W=L o T^(-1)`, but fail the naive
  second-iterate transport.  THM-3532 proves all-level transport for the
  honest conjugates `T o F o T^(-1)` and, more generally, every polynomial
  conjugate of the fixed map.  Covariance for arbitrary tame-equivalent
  Keller maps remains open.
- P3 (reframed): within a mixed grade, primitivity proves atomness; Jelonek
  reducibility is a second one-way classifier, not an atom-degree floor.
- P4: component count of the Jelonek divisor = the grade exponent `k`,
  generically.  THM-3530 proves this for the fixed raw tower at every `k`.
  Any family-generic or monoid-wide form remains open.
- P5: the genus axis governs the remaining JC(2) strata.

## Concurrent scorecard (typed after audit)

The set law and image equation are canonized in THM-2576 and independently
audited coefficient-by-coefficient.  THM-2582 now canonizes the exact norm
identity and the global composite discriminant square class `[H]`; THM-3495
adds the third prime `J`, the three-component set, and `[-2J]`; THM-3498 adds
the fourth old-`L` valuation, polynomial localization, degree-81 genericity
gate, and `[2G]`; THM-3504 adds the fourth prime `G`, image multiplicity one,
and the four-component set; THM-3506 adds the conditional five-face matrix
`(e,m)->(7e-2m,3e-2m)`, the exact pair `(271,99)`, and the next old-`L`
valuation/localization; THM-3513 closes the two renewal faces of `G`; and
THM-3521 adds `v_L(N(R_5))=-1699`, polynomial `R_6`, and its exposed top
pair; THM-3522 adds conditional fixed-chart renewal and the complete packets
of `R_5` and `R_6`; and THM-3523 adds `v_L(N(R_6))=-10663`, polynomial,
`L`-coprime `R_7`, and `A(66907,24255)`.  THM-3528 promotes raw packet
polynomiality to all levels, THM-3529 promotes finite-sheet units to all
complete packets, and THM-3530 promotes the fixed raw image-prime/component
tower to all levels.  THM-3531 promotes the intrinsic discriminant class,
THM-3533 fixes the newest normalization exponent at one, and THM-3535
promotes full wreath monodromy and every nonzero constant-linear primitive
view to all levels.  THM-3532 proves full fixed-conjugacy-orbit covariance
with a sharp one-step/tower boundary.  Constant-view indices, exact old-prime
multiplicities, and wider covariance remain open.
Trace-zero persistence at depth two and the `8 -> 512` cubing pattern are
VERIFIED for the fixed construction only.  The W1/W2 identities are exact for
those two target postcompositions, not a classification theorem.  Their honest
conjugates are covered by THM-3532.  P4 is now proved for the fixed raw tower
at every depth.  The next decisive tests are constant-view indices, old-prime
discriminant multiplicities, and transport beyond the fixed conjugacy orbit.

## Loss ledger

The saturated cuspidal law is proved for the sporadic `F` only; the raw
pullback has the extra affine plane `c=0` described in (C).  The set-level
composition law and odd-block square-class lemma are general.  THM-3530
proves the fixed raw component tower and its generic-degree-one image maps at
all depths, while THM-3531 proves its intrinsic discriminant square class,
THM-3533 proves newest-prime normalization multiplicity one, and THM-3535
proves all-level constant-linear primitivity.  THM-3532 carries these data
around the fixed polynomial-conjugacy orbit, provided the five-weight chart,
raw scalar normalization, and transported observation travel too; none is a
family-generic law.  Constant-view indices, exact old-prime multiplicities,
the monoid-wide `-(det J)^2` law, and covariance between distinct tame-
equivalence classes remain outside proved canon.
"Generic component count" still needs a precise family parameter space before
any classification claim; none of this closes JC(2), classifies maps within a
numerical grade, or controls the weighted G1 witness -- it supplies invariants
and a hunting map.
