# Norm-seven arithmetic and F1: twists, exceptional primes and boundary lifts

2026-10-26. **PROVED** results with independent proof audits and separate
**FINITE-EXACT** controls. Claims retain their stated fields and models.
No new theorem ID is reserved. The inherited THM-4773 correction is included.

- [Fixed-twist central boundary and resolution](boundary/RESULT.md): the
  complete central divisor has zero, eight or infinitely many rational
  points; its eight rational nodes need not lift to the smooth resolution.
- [Gale arithmetic](gale/RESULT.md): rational quadratic twists of the outer
  double-plane K3 are classified by square classes; NS has Galois representation
  Q^18 plus the quadratic twist character. The Gale minus-one twist is
  genuinely nonisomorphic over Q. Nonsquare constant twists lose the free
  rational section in the marked elliptic fibration.
- [Marked F1 Laurent lifts](f1/RESULT.md): all-rank determinant criterion,
  one compensating root character, and the stricter affine-cone obstruction.
- [Reduction at three](reduction/RESULT.md): inertia-corrected missing-sheet
  law; both norm-seven pairs lose freeness at the extra characteristic-three
  triples, while normalization and its conductor survive.

## Inheritance and concept board

Anchor: arithmetic content of the complementary norm-seven quotient maps.
Niche: prescribed-twist boundary points and exceptional conics. Wildcard:
root-order congruence lifts and the different roles of primes three and seven.

Closest proved mechanisms: complementary Fano codes, the full character
algebra with conductor, outer K3 Picard rank19, and fixed-twist conic tests.
Canonical hostile: rational coarse quotient points with a nonsplit residual
torsor. Corrected near misses: a presentation determinant mistaken for an
intrinsic obstruction, and open-central none-or-infinite promoted to the
complete divisor. Least-used sidecars: inertia intersections with each kernel,
Galois action on divisor classes, and arithmetic of exceptional conics.

| Live concept | Operation / predicate | Information to retain / decisive test |
|---|---|---|
| Complementary norm-seven covers | Normalize the coarse product; recover the field | Inertia and residual torsor; characteristic-three triple136 |
| Gale minus-one twist | Rational descent of the outer double-plane K3 | Galois divisor classes; geometric rank19 and anti-invariant section |
| Central boundary | Fixed-twist rational lifting and resolution | Closed triple ports and tangent conics; finite-eight witness |
| F1 torus coordinates | Integral lift of marked root characters | Root order, determinant class, units and positive cone |
| Twist counting | Same-source quotient product and normalization | Full joint distribution; the mean alone erases missing sheets |

The incoming/current packets were recovered before derivation:
`normseven_f1_20261025`, `normseven_boundary_20261025`,
`twist_conductor_reverse_20261025`,
`fixed_twists_conductors_reverse_potential_20261025`, and Session63 twist
mass extend the original Gale/Session60/Session61 sources substantially.
Their Laurent repair, conductor atlas, stack identity, conic criterion,
elliptic decomposition and twist means are inherited, not new results here.

The board changes together: the determinant criterion concerns ambient
coordinates and cannot erase a residual torsor; characteristic-three
inertia changes its dimension; a rational point of the normal source may
still fail the exceptional-conic test; geometric isomorphism of K3 twists
does not erase their Galois divisor data. These are maps with distinct
domains, not numerical identifications of six or seven.

Method cards used from `00-navigation/META-PATTERNS.md`: search the statement
before the method; classify the target over the exceptional stratum before
closing by density; retain power classes before classifying genus; controlled
forgetting requires a sidecar. No new method card is claimed from this one
session. No LRC, integer-square hourglass, or plane-colouring claim follows.

## Audits, reproduction and remaining frontier

Proof reviews and independent computational paths are recorded in
[AUDIT.md](AUDIT.md). All scripts use Python's standard library. From this
directory run `python3 -B verify.py --replay` to check the file manifest,
paired receipts and replay all eight producer/auditor scripts in both modes.
The complete manifest and source pins retain the final exact bytes.

The work extends current main through `983768efde84`; a later fetch before
publication found no incoming commits to integrate. The original shared
checkout's unrelated files were preserved in place.

Open targets remain the interior rational/integral hourglass problem,
nonconstant K3 twists, and an arithmetic atlas of resolved quotient maps
retaining every exceptional conic. Laurent coordinate repairs do not imply
whole-surface torification. The determinant theorem is a classical congruence
mechanism applied to these marked root algebras; the Gaussian norm criterion
is classical too. No historical-priority claim is made for either.
