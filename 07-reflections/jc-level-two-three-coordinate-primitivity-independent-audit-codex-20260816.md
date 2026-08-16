# Level-two three-coordinate primitivity: independent audit reflection

**Historical research reflection, 2026-08-16.**  The proof source is
[THM-3508](../01-canon/theorems/THM-3508-level-two-sporadic-keller-three-coordinate-primitive-discriminant-square-class.md),
not this note.

## Inheritance pass

The anchor was the missing level-two statement that all three source
coordinates of the fixed sporadic Keller square generate its generic
degree-nine function field.

- **Closest proved mechanism:** THM-3494's trace-form identity says primitive
  views of one separable extension have discriminants differing by squares.
- **Canonical hostile:** a coordinate of the intermediate point `F(x,y,z)`
  sees only the degree-three block quotient.
- **Corrected near miss:** MISTAKE-413 forbids replacing a full discriminant
  square class by its odd divisor after silently deleting a constant unit.
- **Least-used sidecar:** the unsquared determinant of a coordinate power
  basis, computed in a regular representation rather than inferred from a
  resultant factorization.

The inherited theorem THM-2582 already supplied the exact `x`-class `[H]`.
The only genuinely missing implication was therefore primitivity of `y,z`;
once that gate was paid, trace-form congruence transported the class.

## Portfolio and concept board

The live lanes were:

```text
Anchor:   generic x/y/z primitivity for F^2;
Niche:    legality of a disconnected but etale specialization;
Wildcard: a direct split finite-field fibre with no inverse formulas.
```

The concept board stabilized at five objects.

| object / representation | predicate | operation | lost coordinate / obstruction | cheapest decisive test |
|---|---|---|---|---|
| generic cubic-over-cubic tower | coordinate generates degree 9 | power-basis determinant | specialization denominators | one lawful nonzero determinant |
| target `(1,1,1)` rank-3 etale algebra | inverse tower is defined | exact regular representation | connectedness | unit norms plus map replay |
| three source power bases | rank 9 and common trace class | basis change | unsquared volume | exact FLINT rank/determinant |
| intermediate first coordinate | expected rank 3 | embed in rank 9 | inner block label | minimal versus characteristic polynomial |
| direct `F_41` split fibre | nine distinct coordinate values | exhaustive map evaluation | nonsplit geometric points | Vandermonde determinants |

Each new result changed the other lanes.  The factorization
`25t^3+t-2=(5t-2)(5t^2+2t+1)` initially looked hostile to the rational
specialization, but clarified the actual logical need: determinant
nonvanishing in a finite etale algebra does not require a number field.  The
rank-three intermediate coordinate then made the primitivity gate visible,
while the direct finite-field split fibre checked the same distinction without
using the inverse section.

## Exact outcome

At target `(1,1,1)`, the clean-room python-flint implementation reconstructed
both inverse stages from the map equations and verified all denominators as
units.  Relative to the cubic-over-cubic tower basis, the power matrices for
`x,y,z` have ranks `(9,9,9)` and nonzero exact rational determinants.  Their
minimal polynomials equal their degree-nine characteristic polynomials.

The intermediate first coordinate instead has rank `3`, zero nine-column
determinant, and characteristic polynomial equal to the cube of its cubic
minimal polynomial.  This is the promised hostile: omitting primitivity turns
a degree-nine coordinate resultant into a repeated cubic with zero
discriminant.

For the three primitive coordinates, exact discriminant ratios equal squared
basis-volume ratios.  Each specialized discriminant divided by
`H(1,1,1)=951326441195` is separately an exact rational square.  Generically,
the nonzero specialization determinants prove primitivity, and THM-2582 then
supplies the full class `[H]`.

The orthogonal direct control enumerated the original `F^2` equations on
`F_41^3`.  Target `(13,0,11)` has nine split source points; all three source
coordinates take nine distinct values, with Vandermonde determinants
`(1,14,12)`, while the intermediate first coordinate has three values of
multiplicity three.

## Connection contract

| field | audited content |
|---|---|
| source | a primitive coordinate power basis in the degree-nine tower |
| target | the common class `[H]` in `K*/K*2` |
| map | trace Gram congruence under basis change |
| preserved | permutation sign and the full discriminant square class |
| destroyed | labelled sheets, exact multiplicities, projection collisions, unsquared volume |
| sidecar | exact basis determinant and saturated degree-nine eliminant |
| hostile | the degree-three intermediate coordinate |
| decisive test | nonzero `9x9` determinant at one lawful specialization |

Eliminant normalization is safe only after the degree and saturation gates:
for degree nine, `Disc(cM)=c^16 Disc(M)`, and `c^16` is a square.  An
unsaturated resultant can carry chart factors; a nonprimitive coordinate can
produce a repeated lower-degree polynomial.  Scalar parity repairs neither.

## Method and remaining frontier

The session used the existing META-PATTERNS cards **Use redundant paths as
detectors** and **Audit saturation and basis covariance before naming a
lattice bridge**.  No new method card is proposed: this is strong evidence in
one Keller thread, while the repository already has the relevant general
cards.

Direct progress is THM-3508.  The niche result is the explicit explanation of
why a disconnected etale fibre is lawful for generic determinant
nonvanishing.  The wildcard became a second exact implementation and a compact
split-fibre witness.  The honest remaining frontier is level three and above:
no computation here proves that every later source coordinate remains
primitive, that every raw eliminant is globally saturated, or that the newest
image prime alone controls every coordinate class.
