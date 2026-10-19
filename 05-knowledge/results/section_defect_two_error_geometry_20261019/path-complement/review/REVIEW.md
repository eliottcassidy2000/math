# Independent audit: path-complement addendum

Date: 2026-10-19. Verdict: **PASS. No mathematical correction required.**

The signed-complement construction removes the same-direction normalization exception in the frozen two-extra-error packet. Consequently, disjoint mixed-cube supports and homogeneous incidence dimension four are sufficient for the two-double-cycle exclusion at p = 29, without an additional sparsity condition. The numerical bounds remain **89 and 445**. The rank-six hostile and overlapping supports are not excluded by this criterion.

## Frozen scope and provenance

The original audited input hashes are retained below; published presentation hashes appear in `review-manifest.json`.

Frozen README SHA-256:

`546196c69c31bdb368ad562d6c2a9fca734cd6b211d0829236e9f4db71c205c3`

The four publication files are `README.md`, `verify_path_complement.py`, `verify_path_complement.out`, and `path-complement-exact-results.json`. Their initial hashes, together with the publication-list hash, are in this directory's `input-sha256.txt`. The original packet remains separately frozen at README hash `2ea9c2a21fec044d5424443998ef6b6710c281a1dc214d456372e448772e1960` and is reviewed in `../../review/REVIEW.md`.

## 1. The possible obstruction is localized correctly

For two same-direction dipoles on distinct e-cycles, a nonexceptional component d has correction supports on L1, L1−d, L2, L2−d. A shared column occurs exactly when the nonzero quotient difference of the lines is ±d. There can be only one shared column for a fixed component: two coincidences would require 2d to vanish in the quotient.

The listed quotient projections are correct: [1,9], [28,26], and [20,3]. In each pair, both entries are nonzero and unequal up to sign. Therefore at most one of the two nonexceptional components can have a shared column. This argument uses the specified directions at the specified prime, not an asserted general-position principle.

Each original segment has at most 14 points. Hence exceeding the 27-position column bound requires both lengths to be 14 and the two segments on the unique shared column to be disjoint. All other configurations already meet the column lemma. Distinct exceptional cycles are essential here and are guaranteed by the remaining 89-failure distribution.

## 2. The signed complement repairs every affected component at once

Replacing S1 by S1−1_(L1) gives minus the indicator of the 15 complementary positions. Its e-derivative is unchanged, since the whole e-line is e-invariant. Its mixed cube is unchanged because the triple derivative includes D_e. Both equalities are exact over the integers.

On the crowded shared column, the other 14-point segment lies inside the first segment's 15-point complement. Thus the union shrinks from 28 to 15 positions. All other columns of that component are unshared, and the other nonexceptional component has no shared columns at all. Their support lengths are therefore at most 15. In the exceptional component the derivative itself is unchanged and occupies at most four points.

The resulting support is always in at most four columns, with at most 27 positions per column and 29 > 2·4. Any cancellation after multiplication by the coefficients can only shrink the support, so the argument does not need assumptions about their signs or nonvanishing.

The revised potential keeps the same integer mixed derivative, all total sums, and all individual d-cycle sums. The previous binary-normalization argument consequently applies. The local-motion identity also remains exact with the revised potential, placing the same impossible seam triple in the original carry coset. The construction does not discard the global period obstruction or replace it with compact-stress sufficiency.

## 3. The reduced universe and all 24 repairs check exactly

Translation puts P1 at zero. Any potentially crowded geometry has two paths of length 14 and P2 = sign·d+h·e. The unequal quotient projections make the choice of nonexceptional d and sign unique; h is then unique. Thus the stated 3·2·2·29 = 348 cases cover the reduced universe without omission or duplication.

An independent dense-grid implementation checks all 348 geometries pointwise, without importing submitted code. It verifies every before/after column certificate, every correction's individual cycle sums, the precise signed-complement support, and exact preservation of both the dipole and mixed cube. All submitted records agree exactly.

Exactly 24 cases require repair, those with h = 14 or 15 for each of the twelve direction/sign choices. Each repaired nonexceptional component has at most 15 positions per column. The remaining 324 cases already meet the original bounds.

## 4. A second representation checks the full short-path universe

The independent audit also extends beyond the reduced universe. In the coordinate basis (e,f), with f the next listed direction, it represents each e-column as a 29-bit cyclic set. Every relative second origin has coordinates (s,t), with t nonzero because the exceptional cycles are distinct. Both short lengths range independently from one through fourteen.

This covers exactly

3 × 28 × 29 × 14 × 14 = **477,456 ordered geometries up to translation of P1**.

Every correction support is rebuilt from the coordinate shifts of the two original masks. Only 24 geometries have a 28-position column, and every one is repaired by the same signed complement. All 477,456 final geometries satisfy the normalization bounds; all 24 repaired ones satisfy the stronger 15-position bound.

This is a support enumeration, not a section search or rank enumeration. It proves no further exclusion of rank-six or overlapping-cube geometries. Its purpose is to independently stress-test the claimed absence of any remaining same-direction sparsity exception.

As a deliberate hostile control, replacing the path with the *positive* complement reverses the endpoint dipole and is rejected by the exact-preservation check. The negative sign in the addendum is necessary and correct.

## Reproduction and limits

Run `python3 -B independent_controls.py` in this directory. The implementation uses explicit exception checks, so they remain enabled under `-O`; its normal and optimized JSON/stdout receipts match byte for byte.

The unchanged submitted verifier is replayed in `submitted-replay/`. Both its JSON and stdout reproduce the published candidates byte for byte. **The submitted verifier uses assertions and must be run without `-O`.**

The addendum does not establish a section with 89 failures, a 90/450 lower bound, a chromatic-number consequence, or a theorem for other primes/direction triples. Its conditional exclusion claim and stated limitations are accurate.
