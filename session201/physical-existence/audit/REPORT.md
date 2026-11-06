# Independent audit: adjacent-cut inverse and even diagonal repair

**PASS, without mathematical correction.** The coefficient-three adjacent-cut inequality, tighter full/balanced/joint inverses, specified-word sharpness, even-index diagonal witnesses, and four-profile distinction are valid as stated. The sharper inverse preserves existence and minimum precursor J costs; it does not establish universal occupancy. The optimal coefficient among all possible representations of a target is not proved.

This report covers the frozen main proof and verifier. Any later asymmetric-budget or letter-specific corollary requires a separate input record and review.

## Frozen inputs

All source and dependency hashes are in `INPUTS.json`. Principal inputs:

- `PROOF.md`: `741bbb44bfc26315bc0aedc4356bab99c3fe5328444e06265cf129d6988c2ea0`.
- `verify.py`: `61e9c7254cfac277b737a05f98155c1b95a29dbfd70091862d0cdaba33badf52`.
- Both author receipts: `d4955d9b23c73ed4f6ca0572a425ac167210db1b64f547d3cc389b06ad164325`.

The physical-adjoint, synchronized-inverse and Carmichael-fiber dependencies match their earlier frozen audited versions. Their physical-source, common-word and budget conventions remain load-bearing.

## 1. The matrix inequality and cut orientation

The four raw images of the nonnegative coordinate basis are `(1,0,1,0)`, `(0,0,1,0)`, `(2,1,2,1)` and `(1,0,2,1)`. Independent raw matrix multiplication recovers all three displayed A and B rows and residual matrices. Every residual entry is nonnegative. Thus multiplication by two nonnegative coordinate vectors proves

    U(c v) U(sigma(c) w) <= 3[K(c v,w)-d_v d_w]

throughout the stated cone, in particular on actual physical sources. The subtraction uses the deletion coordinates of v and w before crossing, not those of the two transformed endpoints. No arbitrary cone point is promoted to a physical source.

At consecutive cuts the pair is exactly `(v,sigma(c)w)` followed by `(cv,w)`. Left cover strictly increases and right cover strictly decreases. At the first crossing of these covers, the two candidate maxima are respectively B=U(sigma(c)w) and A=U(cv). Therefore the smaller candidate maximum squared is at most AB. Since physical d_v,d_w are positive integers, AB<=3(n-1). The argument does not combine unrelated representations or independent boundaries.

Endpoint cases are complete. The full-C1 empty word gives t,t at height two and target three. Otherwise full words cross. For balanced/joint words the right terminal seed has height three, so every nonempty word crosses; the empty word gives heights two and three, target five, and satisfies the bound directly. The separate balanced singleton value three is handled independently.

The joint auxiliary vector moves along the identical suffix and remains coordinatewise bounded by its physical boundary by the inherited synchronized theorem. The chosen cut only repartitions a word, so the two stored J costs add to its original cost. This proves the all-order bound for a supplied witness, not just an empirical smaller inventory.

## 2. Exact inverse and target branches

Sufficiency is inherited literal-word reconstruction; necessity is the tighter cut theorem. Minimum-cost versions are exact because a minimum-cost word can be cut without increasing its cost and every pair reconstructs a word with the sum of its two costs. In the joint version the complete pair (w,q) retains its own minimum cost. Deduplicating by w alone remains invalid.

For a Mersenne branch n=m+2(2^k-1)h, the strict joint inequality m<2h gives n>2^k m. Replacing m by floor(n/2^k) therefore gives the stated sufficient numerical cutoff on each unweighted physical factor. The finite k bound follows from m>=5 and h>=3. The weighted boundary itself need not be small. No new parent is implied for a target that was previously absent.

## 3. Sharpness has exactly its declared scope

For the displayed family, direct iteration gives w=(2d-1,d-1,2d,d), with d=2^(k+1)+1. The left source is `(3,0,3*2^k,1)` when k>=1; at k=0 it is `(3,1,3,1)`. The proof explicitly retains this endpoint exception, so it does not substitute a nonphysical k=0 profile.

After crossing O, the paired scalar is 3d^2 and adjacent candidate maxima are 3(d-1) and 3d. Their product divided by `n-d_vd_w` equals `9(d-1)/(3d-1)`, tending to three. This proves sharpness of the one-letter coefficient on actual physical sources.

For the specified word `Z O^(k+1) J O^(k+1)`, all other cuts have a maximum no smaller than the two neighboring candidates, so its exact best cut height is 3(d-1). The ratio of its square to n-1 tends to three. Removing its final O and changing the right seed to b produces exactly the preceding full-word cuts: dagger of the removed O acts as Z and takes t to b. Hence the same fixed-word obstruction applies to balanced and synchronized cuts. None of this rules out a different representation of the same target with a smaller cutoff.

## 4. Even diagonal family and distinct profiles

For p>=0,q>=1, the middle profile after O^p Z^q has exactly the stated coordinates. The final O^p scales only U and clears E when p>0; at p=0 it is the identity. Because G is independent of E, both cases yield

    G(t.O^p Z^q O^p)=3+2^q(2^(2p+2)-1).

Setting p=(a-2)/2 and q=a+2 is integral for even a>=2 and gives the required diagonal value `(2^(a+1)-1)^2+2`. The zero-J precursor and closing OJ are literal physical words. Their graph order is 2a+4; the original diagonal closing has order 2a+3.

The four F recipes all have the same target. A is the original diagonal; B and C are the inherited O/J twins; D is exactly the dagger of B, which preserves F by the adjoint. Their numerical distinctness is fully checked, not inferred from different words. C alone has deletion greater than one. B and A have different endpoint coordinates. The remaining D-versus-A/B comparisons use the displayed U formulas. For r=2^p, the first difference is positive; the second is -14 at r=1 and positive for r>=2. The signs are never zero. Thus all four numerical profiles are distinct, without asserting graph nonisomorphism.

The a=2 boundary correctly includes a D recipe ending in O rather than Z; its U formula and deletion distinction remain valid. The proof does not assert that its construction repairs odd a. The whole spectrum and general nonempty-fiber conjecture remain open.

## 5. Independent reproducibility

Author runs at cap 5,001 were replayed under ordinary Python and Python -O. Both outputs exactly match both frozen receipts.

The separate `independent.py` imports no author modules. It uses raw matrix generators, increasing-height bucket inventories instead of the author's heap, independent exact-cost image maps, and direct DFS enumeration of directed paths rather than subset dynamic programming for graph controls. Its normal and optimized receipts are byte-identical.

The complete finite comparison independently recovers:

- 126,910 physical source profiles through U=5,000.
- Smaller height cutoff 122: 514 physical sources and 398 synchronized states.
- 66,644 accepted target-cutoff pair presentations.
- 29,322 joint values, 2,471 full C1 targets, and 2,458 balanced targets through 5,001.
- Exact minimum J costs identical in direct, full-pair, balanced, joint and Mersenne-branch calculations.

The direct cutoff is justified by F>=U+1, not by sampled word length. Every orbit transition raises U, so the increasing-height search includes all profiles at its declared bound. The exact target-dependent cutoff is applied to each accepted pair, not only to the largest target in the batch.

Further independent controls:

- All 48 residual matrix entries and all A/B rows verified by raw-basis calculation.
- 19,682 full/balanced word-cut tests through length eight, with the exact chosen adjacent cut and boundary cases checked.
- 9,841 synchronized word-cut tests retaining both readouts and the coordinatewise sidecar bound.
- Sharpness families checked for k=0..32, including the k=0 endpoint and equality of shifted balanced cuts.
- 169 palindromic-family parameter checks and 80 diagonal numerical profiles across even a=2..40.
- Ten literal graph DFS controls: all four profile closings and the balanced closing at each a=2 and a=4. They have respectively 51 and 963 Hamilton paths and exactly one Hamilton cycle; graph orders also match.
- Unsynchronized pair (15,13), equal-w/different-q states, and the no-k=1 hostiles 443 and 617 remain excluded in the correct scopes.

From this audit directory:

    python3 independent.py > check-normal.json
    python3 -O independent.py > check-optimized.json
    cmp check-normal.json check-optimized.json
    cmp check-normal.json INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

The all-order claims follow from the matrix positivity, cut argument and explicit symbolic family, not from finite agreement. The square-root coefficient is only proved sharp under the specified restriction, and universal coverage is not concluded.
