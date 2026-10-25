# Independent review of the root-deletion, branch-cube and Walsh packet

2026-10-25. **PASS** for the stated characteristic-zero claims D1–D4 in
[DELETION_AND_CUBE.md](../root_deletion/DELETION_AND_CUBE.md), and for the
scope boundaries of its canonical/euler formulas. Reviewer: independent
Fano/code lane; no producer implementation was imported by the checker.
The exact producer proof and script bytes reviewed are pinned in
[root_independent_receipt.json](root_independent_receipt.json).

## Proof review

1. **Finite cover versus resolved quotient.** The forgotten root is a
   square root of a section of `O(2)` on the six-root surface. The other
   six roots cannot all vanish at a projective point, so the projection
   is everywhere defined, and its algebra is `O + O(-H)`. The branch
   misses every pre-existing six-root node, since the omitted line
   passing through such a triple would make four concurrent lines.
   An incident triple instead gives a smooth six-root chart and a
   reduced two-branch node. The other four nonzero coordinates give
   `2^4/2=8` branch nodes. These two node types are correctly separated.

2. **Local resolution and fixed locus.** In `w^2=st`, substitution
   `t=s b^2,w=s b` is the normalized pullback after blowing up the base.
   It is smooth and is the standard A1 resolution. The lifted deck
   involution sends `b` to `-b` on the exceptional rational curve. Its
   two fixed points lie on the strict branch components; the exceptional
   curve is not pointwise fixed. The quotient is `Bl_(8t)(K3)`, and the
   smooth-cover branch class is `2(sigma^*H-sum E)`. Thus the text does
   not incorrectly call the map to the *minimal* K3 finite.

3. **Components and genera.** Restricting six forms to the omitted line
   creates exactly `t` pairs and `6-t` distinct zero points. The square
   classes of root ratios have rank `5-t`, independently of constants
   over the algebraically closed base. Their valuation vectors give
   that rank directly. Consequently there are `2^t` components, each
   of degree `2^(5-t)` over the line, with order-two inertia at every
   one of the `6-t` zeros. Riemann–Hurwitz yields genera `17,5,1,0`.
   At paired zeros the two branches differ in one sign bit, so they
   belong to different components. Hence each component is smooth;
   no self-node has been silently included in its genus.

4. **Intersection multiplicities and class relations.** Every paired
   zero contributes eight nodes, partitioned among `2^(t-1)` cube edges.
   This gives edge multiplicity `2^(4-t)`. Adjunction on the K3 gives
   the claimed diagonal entries. The matrix spectra and rank are
   correct for all `t=0,1,2,3`. The singleton null modes are genuine
   linear equivalences: the two factors `x_a ± c x_b` have product
   proportional to the omitted line's pullback and are each a
   hyperplane section, so they cut out the two half-cubes with
   multiplicity one. These `t` relations are independent. Matrix rank
   bounds the divisor-span rank below, while these relations bound it
   above, giving exactly `2^t-t`. No integral saturation is inferred.

5. **Canonical and topological checks.** On the blown-up K3,
   `K_S=sum E` and `M=H-sum E`, so `K_S+M=H`. This proves `K_Y^2=16`
   without ignoring blowups. The invariant canonical summand is the
   pulled-back K3 form and the anti-invariant summand has dimension six.
   The local differential confirms the sign convention; the global
   coordinate action still needs its determinant factor. The fixed-locus
   Euler characteristic is the sum of the normalized branch components,
   `16t-32`, and `2(24+8t)-(16t-32)=80`. The arithmetic genus of the
   union remains 17. All four genus cases and both hourglass deletion
   types satisfy these equations.

6. **Seven-support interpretation.** The seven nonempty subsets index
   seven nonconstant characters of the conic cube, not seven curves or
   seven independent classes. The text retains the constant positive
   mode and identifies the three singleton modes as actual relations.
   This is the required loss accounting. It makes no false claim that
   the monoid `P^2` itself is a K3 or that support incidence alone
   recovers the curve configuration.

## Independent controls

The producer's normal and optimized checkers were replayed and agree
byte for byte. The independent checker
[check_root_independent.py](check_root_independent.py) verifies:

- all 252 labelled deletions by a direct determinant census;
- genera from independent rational Hurwitz arithmetic;
- rational Gaussian matrix ranks and a separately generated Walsh basis;
- all eight explicit central conics and all their intersections over
  `F_41`, using `i=9`, distinct from the producer's `F_17` control;
- all 5,040 Fano markings against the common cube-translation criterion.

Normal and `-O` independent outputs are identical. The finite-field
calculation is a control of explicit conic equations, not a replacement
for the characteristic-zero geometry. It finds 42 points on each conic,
24 distinct two-component intersections, and exactly the doubled cube.
The checker includes no assertion statements disabled by optimization.

## Shared-action addendum

The Fano packet's section 6a has also been checked. For central label `c`
and opposite-pair differences `d_1,d_2,d_3`, their sum is `c`, and all
are nonzero. A bad triple means one difference equals `c`, equivalent
to equality of the other two. In the free case the three differences
are distinct with nonzero sum and therefore form a basis of `F2^3`.
Thus free action on the seven-root resolution corresponds exactly to
regular translation on the eight *component labels* of the central K3.
The latter statement is not freeness on the whole K3. The nonfree
strata have ranks two and one as stated. The character correspondence
from Walsh supports to canonical labels is an explicit invertible
linear change of the character basis; it does not equate divisor
classes with differential forms.

## Limits retained

No fourfold intersections, characteristic two, nonclosed-field splitting,
or full integral Picard lattice is included. The source's reference to
incoming Session60 accurately acknowledges earlier/concurrent K3
deletion and canonical-character mechanisms. No proof that seven
indexed K3 quotients are pairwise nonisomorphic is present or needed.
No mathematical defect or correction to inherited canon was found.

Reproduce from the repository root:

```
python -B 05-knowledge/results/six_seven_geometry_20261025/audits/check_root_independent.py 05-knowledge/results/six_seven_geometry_20261025/audits/root_independent_receipt.json
python -B -O 05-knowledge/results/six_seven_geometry_20261025/audits/check_root_independent.py 05-knowledge/results/six_seven_geometry_20261025/audits/root_independent_optimized.json
```
