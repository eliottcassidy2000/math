# Independent audit of the golden-language normalization

**PROVED dependency audit; FINITE-EXACT independent controls.** The all-order
normalization theorem in [arithmetic_interior.md, section 6](arithmetic_interior.md)
passes an independent review of its path choices, endpoint cases, activity
conditions, and converse quotient. The proof uses the existing four-K interior
identity, the unrestricted prefix interchange, and the separate sink-context
identity in [moves.md, section 10](moves.md). The finite controls below supplement
that proof; they are not an extrapolation from a bounded graph census.

## 1. The simultaneous path-choice obligation

The important strengthening in the interior theorem is compatibility, rather
than just individual existence of black paths.

* Its third path is always black. The endpoint score condition would force
  `k>=m+2`; at `k>=m+3`, reflecting the W interval contradicts `L -> 1`, while
  the remaining equality case contradicts the two arcs `3 -> 0` and `1 -> R`.
* For its second path, a first-SCC vertex of W with score different from `m+2`
  supplies a black starting choice. If every such vertex has score `m+2`, a
  hypothetical blue path sends an A vertex to a first-SCC W vertex whose
  original W score equals its post-K A score, at most m. This is incompatible
  with that common score. The first-SCC lower bound ensures the reflected
  vertex really is in the first SCC and is not an endpoint of the first path.
* The source and target cannot both be self-converse: the replacement flips
  three distinct self-converse SCC-type parities, whereas a self-converse SCC
  word has at most one odd such coordinate.
* If the source is self-converse, a mirrored R4 path can be varied independently
  of W's chosen starting vertex. When that R4 is W's first component, the
  already-proved small-first-component case instead makes the second move
  black for every choice. If the target is self-converse, its black path fixes
  the context paths, but its reflected SCC form forces the first W score to
  be `m+1`, `m+3`, or at least `m+4`, never `m+2`.

Thus one pair of context paths makes all four moves black simultaneously. No
unjustified independent-choice inference remains.

## 2. Boundary and word arguments

If the ordinary right-end target `A -> C3 -> 1` is pure blue, it is a rigid
SCC palindrome. Its initial component must therefore be a singleton. The
separate boundary ledger applies to its converse and reaches `A -> 1 -> C3`;
this target is nonself-converse because its first and last SCC types differ.
This includes the smallest context `A=1`. Erasing an empty context from the
interior ledger is neither needed nor asserted. The isolated R4 has its own
two black path presentations, each reaching the converse-merged `1 -> C3` type.

For an interior last R4, its two replacement words cannot both be palindromes:
their disagreement set is exactly two adjacent unequal letters. Reflection
would have to preserve this set, which would put unequal letters at the middle
pair of a palindrome. If another R4 remains, activity is automatic.

After dissolution, `r=n-3f` singletons and `t=f` triangles remain. Activity
requires and is possible exactly when both inventories are positive. A rigid
word is active precisely when it is not palindromic. Its inversion number I
changes to `rt-I` under converse. Below the midpoint, bubble downward; above
it, bubble upward. Those routes stay nonpalindromic. A nonpalindrome at the
midpoint has an adjacent inversion and moves immediately below the midpoint.
The two sorted extremes are converse-equivalent. This proves one ambient
black component contains every active word of the fixed inventory.

The theorem does **not** assert that the induced terminal-word graph is
connected, that f stays fixed at intermediate tournaments, or that different
f values label different connected components. It preserves n throughout;
f is preserved at the authenticated rewrite endpoints.

## 3. Independent exact controls

[arithmetic_normalization_audit.py](arithmetic_normalization_audit.py) imports
no producer code. It implements adjacency and K moves independently, checks
all path arcs and final labelled adjacency matrices, and verifies all four
endpoint score excesses in every standard-path sink context of orders 1–7.
The universe has **1,101 presentations**, including the collapsed singleton
port. Every case gives `(-1,4-n,n-3,1)`, so all four paths are black.

The separate symbolic control enumerates **all 142,129 ordinary words** on
`{1,C3,R4}` of vertex weight at most 24. It normalizes all **141,726 active
words** across **84 nonempty (n,f) strata**. It exercises 1,612 cases needing
the alternate interior target, 645 right-boundary repairs, 86 left-boundary
repairs, and 3,934 nonpalindromic midpoint starts. Every one of its 2,038,158
adjacent rigid swaps has active endpoints. All words of each inventory have
the same converse-merged normal form. These word tests use the proved macro
theorems; they do not enumerate every intermediate tournament of every route.

Reproduce from the repository root:

```text
python 05-knowledge/results/scc_basins_20261020d/arithmetic_normalization_audit.py
python -O 05-knowledge/results/scc_basins_20261020d/arithmetic_normalization_audit.py
```

Normal and optimized runs passed with byte-identical JSON receipts. Their
receipt SHA256 is
`f585c8c35acad229a73aaa7db208188cf9d32c81ac254f693d871cbe11cc122d`.
The [receipt](arithmetic_normalization_audit.json) includes the complete
by-order census, branch counts, direct R4 certificates, and source digest.
