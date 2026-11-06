# Independent audit: actual dyadic-class spanning addendum

**Verdict: PASS for the stated product-of-independent-classes result.** The addendum validly extends the earlier formal obstruction to actual identities holding on a full Cartesian product of two dyadic residue-class tails. It does not exclude coupled indices or thinner subsets. The only issue found is a non-substantive terminology point: the homogenized vectors have **linear span/rank three**; the original planar points have affine dimension two, not three.

This is a separate addendum. The existing main and digit audit reports and their independent verifiers/receipts are unchanged.

## Frozen inputs

The exact originals are preserved under `frozen-inputs/`:

- DYADIC-SPAN.md: `efb7ad8eefe15b68bd3eb7862ec585b03ef88f71a11a36822c2d83bc5bf4fbc7`.
- check_dyadic_span.py: `10f9902931ad1568df127c3313c019fb028bd6484dba7719d0f04591bb4fd9e0`.
- DYADIC-normal.json and DYADIC-optimized.json: both `413b7fcc8f0890dea35a911ded655ed584e5496f8da975fbf7112d9ff2ca7b32`.

No original input was modified, and no GitHub or Git-index operation was performed.

## Independent proof review

### 1. Exact binary-block normal form

Let a word of length k have bits epsilon_i in {0,1}, with O=0 and Z=1, indexed i=0,...,k-1. Direct expansion gives

- A=2^k,
- r=sum_i epsilon_i*2^(k-1-i),
- a=sum_i epsilon_i*2^i,
- b=-epsilon_0 for k>0, and b=0 for the empty word,
- c=-sum_{j<i} epsilon_i*epsilon_j*2^(i-1-j) - sum_{i>=1} epsilon_i*epsilon_(i-1).

Thus `u'=A*h-r` and `g(u')=g(h)+a*h+b*(h mod 2)+c`, exactly as claimed. The current index has fixed parity after the first letter, so only an initial Z can contribute the original h parity. This explains why b is exactly 0 or -1 rather than an arbitrary accumulated parity coefficient.

The map from binary words to r in {0,...,2^k-1} is bijective. Every dyadic residue class is represented, and the constructed image contains its entire sufficiently large positive tail. All chosen h are at least two, so the intended unit-profile index restriction is respected.

### 2. Actual vectors span every tail

For m>=1, the recurrence gives

`g(2^m)=g(2^(m+1))=1`, and `g(2^(m+2)-1)=2^(m+1)+1`.

For the three h values in the proof, the determinant of the base row vectors `(g(h)+b*(h mod 2),h,1)` is

`-2^m*(2^(m+1)+b)`.

Applying the invertible block map multiplies this by A=2^k. Hence the actual transformed determinant has the exact closed form

`-2^(k+m)*(2^(m+1)+b)`.

It is nonzero for every m>=1 and b in {0,-1}. Increasing m makes all three image indices arbitrarily large in the selected dyadic class. This proves the claimed linear span over Q for **actual** g-values, including the empty binary block and every residue class.

The proof's displayed equation (2) and its linear-independence argument are correct. For strict terminology, change the section heading “A three-dimensional affine span in every tail” to “Three linearly independent vectors in every tail” or “Rank three in every tail.” The receipt's phrase “affine span3” should be understood as this homogenized rank-three statement; it is not an additional affine-dimension assertion.

### 3. From a full product identity to coefficient equality

After fixing parities, each expression is bilinear in the three coordinates `(g(u),u,1)` and `(g(v),v,1)`. The additive constant is represented by the coefficient of the product of the final coordinates. Separate O/Z/sibling compositions have fixed affine formulas on dyadic residue subclasses, and their branch/parity choices can be fixed by such refinement.

If `X(u)^T C Y(v)=0` for all sufficiently large independent indices in a product of two selected classes, choose three of the proved independent X vectors and three independent Y vectors, all beyond the required threshold. Writing these vectors as the rows of invertible matrices P and R gives `P C R^T=0`, hence C=0.

This reasoning genuinely needs the Cartesian-product hypothesis: all nine selected index pairs must satisfy the identity. A coupled-index locus, a finite set, or another thinner subset does not supply that implication. The addendum explicitly retains this boundary.

The earlier formal coefficient obstruction then applies exactly as stated, including the empty central-word cases and exchange of the indices. No additional target-occupancy or spectrum conclusion follows.

## Independent exact computation

`check_dyadic_independent.py` does not import the author checker's affine or determinant functions. It uses the closed bit-sum formulas above, iterative recurrence evaluation, its own determinant calculations, and full 9-by-9 product-evaluation matrices.

Checks passed:

- All 2,047 O/Z blocks of length at most ten, with complete coverage of every residue class at each length.
- 202,653 exact block-profile controls, using every seed h from 2 through 100.
- 40,940 actual span determinants, using m=1,...,20 for every tested block. Every determinant agrees with the closed formula above.
- 225 full 9-by-9 product-evaluation determinants, covering every ordered pair of the 15 residue classes from moduli 1,2,4,8. Each equals `det(P)^3*det(R)^3` and is nonzero. These explicitly test the independent-product inference as well as the separate three-vector spans.

Both independent Python runs are byte-identical. Normal and optimized author replays are byte-identical to the original supplied receipts. The all-class and all-tail statements rest on the algebraic proof, not these finite test ranges.

## Reproduction

```bash
python3 session201-prefix-repair-audit/check_dyadic_independent.py > session201-prefix-repair-audit/dyadic-independent-normal.json
python3 -O session201-prefix-repair-audit/check_dyadic_independent.py > session201-prefix-repair-audit/dyadic-independent-optimized.json
cmp session201-prefix-repair-audit/dyadic-independent-normal.json session201-prefix-repair-audit/dyadic-independent-optimized.json
python3 -B session201-prefix-repair-audit/frozen-inputs/check_dyadic_span.py > session201-prefix-repair-audit/dyadic-author-replay-normal.json
python3 -B -O session201-prefix-repair-audit/frozen-inputs/check_dyadic_span.py > session201-prefix-repair-audit/dyadic-author-replay-optimized.json
cmp session201-prefix-repair-audit/frozen-inputs/DYADIC-normal.json session201-prefix-repair-audit/dyadic-author-replay-normal.json
cmp session201-prefix-repair-audit/dyadic-author-replay-normal.json session201-prefix-repair-audit/dyadic-author-replay-optimized.json
```

All reproduction commands and comparisons passed.

## Verified terminology-only revision

After the proof audit, the author made the recommended terminology correction. Exact diffs confirm only these changes:

1. The section heading now says “Three linearly independent homogenized vectors in every tail.”
2. The checker receipt's scope string now says “actual homogenized linear span3.”
3. The regenerated receipts record that new scope text and the corresponding checker hash. All numerical fields are unchanged.

No proof formula, algorithm, test universe, or main/digit packet source changed. The corrected files are preserved separately under `terminology-inputs/`, leaving the original frozen inputs intact.

Verified revised hashes:

- DYADIC-SPAN.md: `af63b6279c32d73e0bd41abd803cc4dbb126375688d35bac4b53d305114ecc9a`.
- check_dyadic_span.py: `a121e03109243afe3218a1ac56fd2cb87a2b9b90c951b8a36d2179042533108d`.
- Both DYADIC receipts: `dd30a0decff4e51eef6c13f74a31b4c984c8475ceca47e990f1797babd172149`.

Normal and optimized replays of the corrected checker were independently run and are byte-identical to the revised author receipts. They are saved as `dyadic-terminology-replay-normal.json` and `dyadic-terminology-replay-optimized.json`. **PASS carries forward to this exact terminology-corrected version without qualification.**
