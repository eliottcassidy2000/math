# Exact one-apex certificates: two independent nonzero shifts

**Proved:** for disjoint bundles at one apex, zero-, one-, and two-bundle
readings determine every Hamilton count and certify gap preservation and
additivity. **Finite-exact:** the displayed physical examples. No claim
of unbounded supply, arbitrary absolute placement, or plane-coloring bounds.

## What to retain from the unit-distance work

The named-port relation in
`boundary_palette_compression_20261011.md` distinguishes boundary supports
that look identical after counting colors. `boundary_context_weights_20261011d.md`
goes further: equal Boolean extension support does not mean equal extension
multiplicities. Here the analogous carrier is the weighted table indexed
by the two vertices where a path touches the apex, not the total two-cover
count. `boundary_fibre_atlas_20261018d.md` also distinguishes local phase
tests from exact compatibility. No geometric embedding is being transferred.

## Certify every setting

Fix the deletion tournament and its initial apex spokes. Write x_i=1
when the apex beats i. Let E_i and S_i count Hamilton paths on the deletion
tournament ending and starting at i. Let B_ij count ordered covers by two
nonempty paths, the first ending at i and the second starting at j.
The inherited quadratic (THM-4097/4114) is

    H(x) = sum_i ((1-x_i) E_i + x_i S_i)
           + sum_(i≠j) (1-x_i) x_j B_ij.

Remove the apex to prove this: it was last, first, or the joining arc
pair between the two disjoint paths. In particular a path depends on
at most two apex spokes. Disjoint bundles assign each spoke to at most
one input bit. Each of the two tournament counts is thus a multilinear
quadratic in the bundle bits, and their gap is also quadratic.

The constant, linear, and quadratic coefficients are recovered by
finite differences on the empty setting, singleton settings and pair
settings. This proves (for any number of bundles) that the gap is
constant iff these readings all have the original gap, and that the
baseline shifts add iff each pair mixed coefficient is zero. Once the
gap is constant, testing additivity on one side suffices for both.
The certificate requires disjoint bundles, or another reason each
spoke depends on at most one bit; arbitrary overlapping edits do not
give this degree argument.

For U,V disjoint, write eta_i=1−2x_i. The mixed coefficient is

    J(U,V) = −sum_(i in U,j in V) (B_ij+B_ji) eta_i eta_j.

The nonnegative weights have the same sign contribution across the
entire rectangle when each bundle is internally homogeneous in its
starting orientation. If any cross weight is positive, such bundles
cannot be additive. Mixed bundles can have exact cancellation. This
sign formula was independently replayed on 1,098 apex completions and
51,984 disjoint-bundle checks in the companion weighted-boundary code.

## Physical examples: each pair is (lower witness, upper witness)

Settings are empty, U, V, U+V in that order.

| Order | Four pairs | Result |
|---|---|---|
| 5 | (1,3), (11,13), (3,5), (9,9) | Both single tests pass for the gap; the joint test is essential |
| 7 | (9,11), (37,39), (17,19), (45,47) | Independent shifts +28 and +8 while preserving +2 |
| 8 | (173,175), (171,173), (71,73), (69,71) | Independent shifts −2 and −102; two exact three-count windows |

The last row yields the windows {69,71,73} and {171,173,175}.
These are not all intervening targets. The seven-vertex row already
uses genuinely nonzero independent controls, but its increments are
too large to create one wide gap-free interval.

For reconstruction, flip each apex spoke to i indicated by bit i
in the mask (apex is the last label). Rows are outgoing bitmasks.

* Order 5: literal bases [0,1,3,7,15] and [10,0,3,6,15]; U=3,V=4.
* Order 7: shared seed [38,64,34,103,111,2,37]; initial masks 30,31;
  U=3,V=36.
* Order 8: shared seed [50,224,123,179,130,16,57,101]; initial masks
  3,112; U=48,V=78.

Run `python3 apex_bundle_replay.py`. It verifies every spoke orientation
around each of the six bases by comparing literal subset DP with the
named-port formula, and checks the displayed corners.
