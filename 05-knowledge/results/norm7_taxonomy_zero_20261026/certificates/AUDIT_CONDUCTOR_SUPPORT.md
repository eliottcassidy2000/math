# Independent audit of odd-prime unramifiedness

2026-10-26. **PASS.** No blocking gap found in the stated full-carrier theorem.
This audit does not identify a pure constituent, rank, conductor exponent,
or individual Euler factor.

Reviewed `../arithmetic/CONDUCTOR_SUPPORT.md`, SHA256
`2aca53f7c50f441758787f674630e3a47ecfaffe158c1cd3d607316979494b6d`.
The reviewed status was “Proof candidate”; a later promotion of that status
alone is outside this hash but does not change the proof audited here.

A subsequent scope-only paragraph was separately checked against the
inherited `normseven_arithmetic_f1_20261026/reduction/RESULT.md`, R2. It
distinguishes the free configuration carrier from the norm-seven hourglass
arrangement, whose triples136 and245 appear in characteristic3 and destroy
freeness of the four marked quotient kernels. The amended candidate has
SHA256 `e1f6ef43a434a72560078dc66f3851f56b7b6edccbc1288dff93c8e09e440bb3`.
This addition correctly narrows transport scope; the PASS remains unchanged.

## 1. Carrier and integral configuration strata

The literal cover is finite etale over the locus with the required
inequalities and nonzero radicand. Removing a zero of any positive-exponent
edge does not change its character sum: that summand was zero. Even edge
multiplicities do not change this definedness condition. Repeated factors
after a merger remain present; the argument does not replace allowed
collisions by an injectivity condition.

Every admissible equality partition has two distinguished blocks at0 and1.
Its other blocks are pairwise distinct variables away from those values.
Thus it is precisely M_(0,r+3), with infinity as the third fixed marking,
over the same integral base. Edges between the fixed blocks contribute only
units1 or−1. No arithmetic constant with an odd prime factor is introduced.

The finite union of strata with at most k total blocks is closed in U:
it is the intersection with U of a finite union of diagonal closed loci.
Its successive difference is the disjoint union of the exact-k strata.
This supplies the closed filtration needed for the final localization step;
the stratification is more than a partition of generic rational points.

## 2. The proper smooth model and the parity blowups

I independently read Knudsen's Theorem2.7 and its proof on printed page179
of the [author-hosted collected PDF](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1976b--DetDiv%2BKnudsenSequel-GT.pdf).
The theorem states properness, smoothness over Z and relative normal
crossings for the stable-curve boundary. The proof starts from
M_(0,3)=Spec Z and uses the universal pointed curve. In genus zero each
stable labelled curve has no automorphism, and the universal-curve induction
gives schemes. A compatible collection of stable label partitions specifies
a tree; its corresponding intersection is a product of the moduli spaces
at the vertices. Incompatible collections have empty intersection. This
explains the *simple*, globally labelled boundary needed by the argument.

The source PDF was also inspected through the local rendered page
`/tmp/knudsen-page-58.png` (printed page179). PDF SHA256:
`56bb67830d2d16a1c3c6e6b9ccd35e406d374f2e36669b938c10f78d45b0b0f2`;
page-image SHA256:
`5d8196eb0061c98646286fd562b46ded99a1c741bbb0b4716648ca1df81d18e8`.

Every coordinate difference is a unit on the integral open stratum, not
merely on its characteristic-zero fibre. Its divisor on the proper model
therefore has no vertical component. Boundary zeros and poles both enter
the integer multiplicity list; negative orders are correctly reduced modulo2
using the floor formula.

For a pair of odd boundary components, the center is their whole smooth
intersection, of relative codimension2. In the two blowup charts,
`x1=t,x2=t*v` and `x1=t*v,x2=t`, a monomial radicand acquires exceptional
order m1+m2. This is even. A third SNC component has coordinate x3 and
does not contain the generic point of that center. Its order is therefore
not added to the exceptional order, even when a triple intersection occurs.
The two selected strict transforms are disjoint. A strict-transform
intersection maps into the old intersection, so previously disjoint odd
components cannot start intersecting. There are no newly odd exceptional
components. These facts prove termination after at most one operation per
original odd pair, without needing a resolution theorem in mixed
characteristic.

The same coordinate charts preserve a relative SNC total boundary and a
smooth proper base model. With all odd components disjoint, write
div(f)=2A+B. The indicated isomorphism L²=O(B), for L=O(−A), must be the one
induced by this **particular f**; it retains its unit and orientation sign.
Under that understood choice, the double-cover algebra has exactly the
given generic radicand. Locally it is z²=u or z²=u*x1. Etale extraction of
the unit square root is legitimate because2 is invertible. Both charts
are smooth over the base.

In the branched chart, the reduced full boundary is locally
`z*x2*...*xk=0`; in the unbranched chart it is an etale pullback of the old
boundary. Hence each labelled reduced boundary inverse image and each of
their intersections is smooth proper, even when disconnected. Taking the
entire boundary, not just the branch divisor, is necessary and is done in
the candidate. The covers with square radicand and r=0 cause no exception.

## 3. The integral compact-support argument

The augmented boundary complex is exact on geometric stalks: at a point on
k boundary components it is the augmented simplex complex, and on the open
it is the identity of the constant sheaf. Pushforwards of the constant
sheaves on the closed intersections are understood in this expression.
The complex is finite and computes j_!Lambda before any passage to the
generic fibre.

I read the exact primary statements:

- [Stacks Proposition63.15.2](https://stacks.math.columbia.edu/tag/0GLG)
  applies to each smooth proper boundary intersection with finite constant
  coefficients Z/ell^a, since ell is inverted on the chosen base.
- [Stacks Lemma18.43.5](https://stacks.math.columbia.edu/tag/093U)
  makes these finite locally constant modules a weak Serre subcategory:
  kernels, cokernels and extensions occurring in the spectral sequence
  and localization long exact sequences remain in that category.
- [Stacks Theorem59.91.11 and Lemma59.91.12](https://stacks.math.columbia.edu/tag/095S)
  identify the geometric-fibre cohomology for the proper maps and bounded
  torsion complexes used here.

All their hypotheses apply. The boundary spectral sequence is a bounded
spectral sequence over S_ell, so its abutment cohomology sheaves are finite
locally constant. Compact-support localization for the finite closed
collision filtration then has the same property. This may be implemented
in a common proper compactification; extension by zero and proper base
change retain the required geometric-fibre identifications. The common
compactification need not itself be smooth because it is used for
localization and proper base change, not for proper-smooth base change.

The inverse systems at ell^a have finite cohomology groups on geometric
fibres, hence satisfy the needed Mittag-Leffler condition. Their limits
and the resulting Q_ell representations inherit trivial inertia at each
odd p different from ell. The minus projector is taken only over Q_ell,
so ell=2 is permitted. The same construction without the cover handles U.

This reasoning does **not** infer unramifiedness merely from unramified
graded generic representations. Each filtration and extension is realized
by locally constant sheaves over the integral base. That distinction closes
the delicate extension gap identified in the audit request.

## 4. Precise outcome and hostile boundaries

The proved conclusion is unramifiedness of every full-carrier H_c^i at odd
p≠ell, and consequently of its actual Galois-stable subquotients. If a
compatible pure constituent is later constructed inside these carriers,
one may choose ell≠p for each odd p to conclude no odd conductor support.
This does not construct that constituent or determine its rank, purity,
irreducibility, local Euler factor, or conductor exponent at2.

The proof correctly fails after arbitrary constant twists: y²=3 is ramified
at3. It also correctly retains infinity in the parity divisor, and does not
allow a new fixed rational marker to collide with0 or1 at an odd prime.
Those hypotheses are material. The raw finite-field trace remains an
alternating sum of cohomological traces, rather than one specified pure
representation. These boundaries should remain visible after promotion.
