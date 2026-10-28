# Recovered March square codec: exact repair and its information ledger

Session 89, 2026-10-28. **PROVED elementary repair**, with separate finite exhaustive checks. All source observations refer to commit `03408f84ac43fcdcfb0887e5d6bc1fc8792415cf`; the two live codecs receive only the narrow repairs described below; saved historical outputs remain untouched.

## The source really is old, but its implementation is not lossless

The March 24 script `04-computation/tiling_square_s325.py` proposes a square made from consecutive fixed-Hamiltonian-path tilings. Its committed `05-knowledge/results/tiling_square_s325.out` actually reports `Filled: 7 / 9 cells` at side three. The code's `build_square` maps smaller-tiling tiles `(i,j)` to `(i,j-2)`, retaining them only when `i<j-2`. Thus every smaller-tiling range-two tile is discarded. It identifies smaller-tiling masks `000` and `101` at size four. `04-computation/new_sequences_s327.py`, lines 198–211, repeats this same mapping in its rank/determinant census.

The historical script also retracts its own earlier claim that the square diagonal is the insertion layer. With its lower-triangle map, the diagonal consists of range-two tiles. The insertion layer is the last row. These are different subsets.

## Corrected all-size map

Let `m>=1`. A fixed forward-path tiling of size `m+2` has free bits `a_ij`, `0<=i<j<=m+1`, `j-i>=2`; a second of size `m+1` has free bits `b_ij`, `0<=i<j<=m`, `j-i>=2`. Define an `m×m` matrix by

```
M[j-2,i] = a_ij,          0<=i<j<=m+1, j-i>=2;
M[i,j-1] = b_ij,          0<=i<j<=m,   j-i>=2.
```

The first index set maps bijectively to `r>=c`, the second to `r<c`. Their inverse maps are respectively `(r,c)->(c,r+2)` and `(r,c)->(r,c+1)`. They are disjoint and exhaustive, so this is a coordinate permutation, not merely a cardinality argument. Consequently it preserves full bit data, Hamming distance, the uniform measure, and independent bit flips. It works over every alphabet, and over a field is a linear coordinate isomorphism.

**Source / target / map / preserved / lost:** source is two *marked-path* tilings of orders `m+2,m+1`; target is an `m×m` array; map is displayed above; all coordinates and cube adjacency are preserved. No labelled information is lost. Ordinary tournament isomorphism, natural matrix conjugation, and vertex deletion are not asserted to correspond.

THM-830 supplies the further gap-tournament bijection: the tiling bit on `(i,j)` becomes an edge on gap vertices `(i,j-1)`. Thus the sources are, by another explicit coordinate map, unrestricted labelled tournaments on `m+1,m` gap vertices. They are **not** obtained by forgetting the fixed path or by taking an induced subtournament of the original tiling tournament. This explains the two-order shift between the historical tiling terminology and the user's square-bit-matrix observation.

## Exact defect of the historical code

The old larger-tiling half is correct. Its smaller half discards exactly `m-1` range-two bits and maps the rest bijectively to the upper triangle excluding the last column. Therefore its image is

```
{ M : M[i,m-1]=0 for 0<=i<m-1 },
```

of dimension `m²-m+1`, and every image matrix has exactly `2^(m-1)` preimages. At `m=3`, its 512 input pairs yield only 128 matrices, each four times. For arbitrary coefficient alphabets the same missing-coordinate description holds. Over any field its matrices have block form `[[B,0],[v,d]]`, hence `det(M)=d det(B)`. The old rank/determinant distributions therefore cannot be read as distributions over all binary square matrices. This is a map-specific correction, not a refutation of the triangular-number identity.

## What this repair does not establish

A coordinate bijection transports an action only when the action is transported too: `g⋆M = F(g F^{-1}(M))`. It does not establish that `⋆` is matrix conjugation. Likewise Hamming-isometry is a labelled statement; quotient Hamming distances require the same specified group action on both sides. Bit rank over `F_2`, ordinary real rank, and integer determinant are three different observables and must be named separately.

The independent checker reconstructs the exact historical coordinate rule, tests the missing-coordinate/fibre formula, checks all corrected round trips through side four, and checks cube-edge Hamming preservation through side three. It does not execute either exploratory historical script or import its producer code.


## Minimal source repairs and live-code regression

`repairs/tiling_square_s325.patch` corrects the illustrative coordinate and `build_square`; `repairs/new_sequences_s327.patch` corrects the repeated inline coordinate. The latter changes exactly one executable constant, independently checked by comparing ASTs after undoing that edit. Both preserve all unrelated exploratory code and add explicit warnings about historical rank/determinant data. No saved numerical output is rewritten or regenerated.

The public `repairs/audit_patched_codec.py` uses repository-relative live sources by default, with `--repo-root`, `--s325-source`, and `--s327-source` overrides. It extracts the actual S325 function and the actual S327 inline constructor by AST, without executing either exploratory script. Both pass complete side-1–4 universes plus coordinate-basis checks through side 12 against the inverse-cell oracle: **399,112 always-active checks**, identical normally and under `python -O`. Each original buggy source is rejected at side two, and both patches replay byte-for-byte. See `repairs/README.md` for scope and commands.
