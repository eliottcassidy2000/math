# Surgical correction guidance for THM-4727

This is replacement guidance, not a source edit. Preserve unaffected mathematical content and historical provenance. Target only `01-canon/theorems/THM-4727-mixed-weight-torus-reduction-of-keller-maps-ramification-exponent-twisted-planar-jacobian-equation-and-the-dimension-profile.md`.

## Banner immediately after the main heading

> **S144 scope correction (2026-11-03).** The Gao–Gallagher tangent-sweep **seed family** has affine orbit scale. It is not the whole known counterexample corpus. Reorder the fixed map as `G=(F_3,F_2,F_1)` to make source and target weights agree. Its iterates are equivariant Keller counterexamples with orbit scales of degrees `1,7,37,187,...`; the scale-zero divisor of the mth iterate has at least m pairwise disjoint smooth components. Thus existence with `deg l>=2` is settled by composition. Classification of primitive solutions and of arbitrary Keller maps remains open. See the S144 composition proof packet. Here `{x_1=0}` is an **invariant coordinate hyperplane**, not a pointwise fixed plane; the fixed locus of this nonzero-weight action is the origin.

Replace “PROVED (pre-audit)” in the summary with “PROVED IN THE STATED SCOPES; S144 CORRECTED for the all-known-family and higher-degree-existence assertions.” The new universal iterate statements have a separate elementary proof and exact companion; do not label every legacy assertion reaudited merely because of this correction.

## Title/status repairs

Replace the title clause “so the whole known dimension-three Keller counterexample corpus is the W = 3, deg l = 1 stratum of (R)” by “so this tangent-sweep seed family lies in the W = 3, deg l = 1 stratum of (R); compositions already leave that stratum”.

In prose and title, replace “fixed plane”, “fixed hyperplane”, and “fixed-plane targets” by the corresponding “invariant coordinate plane”, “invariant coordinate hyperplane”, and “coordinate-hyperplane targets”. These are invariant subsets, not fixed loci.

## Section 2 heading

`## 2. The tangent-sweep seeds have affine orbit scale; their iterates do not`

## Replace the overbroad family sentence

“So every member of the displayed Gao–Gallagher tangent-sweep family, including THM-1300 and Gao's degree-four map, is `(1,-1,-2)`-equivariant with affine orbit scale. This does not exhaust known compositions or classify all counterexamples.”

Keep the displayed quotient formula and fixed-map examples with their explicit scopes.

## Replace the open-question paragraph

“Classifying the `(1,-1,-2)`-equivariant Keller maps reduces to equation (T) together with its polynomial lifting conditions. The displayed tangent-sweep seeds have affine `l`. Higher-degree `l` already exists by composition: put `G=(F_3,F_2,F_1)`, `l=2-3v-t`, `A=v+3(1+v)^2t+3v^2(4+3v)`, `B=(1+v)^3t+v^2(1+v)(4+3v)`, and `H=(lA,l^2B)`. Then `G^2` has determinant `4`, geometric degree `9`, and scale `l_2=l(2-3lA-l^2B)` of degree `7`. More generally `l_m=product_(j<m)(l composed with H^j)`, with degree `(3*5^(m-1)-1)/2`; its zero divisor has at least m disjoint smooth components. These cannot be removed by equivariant two-sided polynomial automorphisms, which have constant orbit scale and induce quotient automorphisms. Even primitive maps with arbitrary-degree nonlinear scale arise by the equivariant source shear `(x,y,z) -> (x,y+x^(2k-1)z^k,z)`, which changes the scale to `2-3v-t-3t^k` but preserves the equivalence class. The remaining questions concern primitive maps inequivalent under equivariant polynomial source/target automorphisms to the known affine-scale seeds, and a classification of nonsweep solutions. No arbitrary-map classification follows.”

## Strengthen the necessary weight gate after Corollary R'

“The derivative at the origin gives a stronger necessary condition than the equality of total weights: since all weights are nonzero, `F(0)=0`, and the invertible `JF(0)` intertwines the two torus representations. Therefore their entire weight multisets agree. Equality `W=W'` alone is not sufficient for an equivariant Keller map.”

## Additional boundary issues for the independent auditor

The universal claim that a coordinate-hyperplane fiber has the full generic degree must be restricted if that hyperplane is itself a nonproperness component. The formula `Jac(Fbar)=c l^(W-1)` only implies vanishing on `l=0` when `W>1`; when `W=1`, it remains nonzero. Neither issue is resolved by the tower alone. Have the independent boundary audit own any further correction rather than silently broadening this patch.
