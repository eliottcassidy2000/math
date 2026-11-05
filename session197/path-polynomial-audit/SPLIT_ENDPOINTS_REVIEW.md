# Independent audit: positive split-endpoint closure

Status: **PASS, no correction required.** This review concerns one distinct endpoint pair retained inside one intact input module, and the separate coincident-endpoint boundary. It is not a review of arbitrary multiport connectivity closure.

## Frozen inputs

Files in `session197-odd-cycle-path-polynomials/`:

- `SPLIT_ENDPOINTS.md`: SHA256 `3229f2a573ff7858439cf5e69ba5e0dd3d323c8ec31c283d4d7749cc7a173f1f`.
- `verify_split_endpoints.py`: `0cb0fb43f84448d050d90b37e87b06a3a24383482f826808b533cb842a26a26d`.
- `split_normal.stdout` and `split_optimized.stdout`: `0f9cf8abc0aaf0de72942fcaab39ced59431ddba86ad80b3eadd059367663914`.
- `verify_matching_amplitude.py`: `e062734acf04e628a5fee3ed15d57123a22bcb6f7a2f1058edd1972f0566773c`.
- `MATCHING_AMPLITUDE.md`: `443097af71c4d1c364c4bc8ad5d417453df6aea71d734ccb5700b74b93bf48ee`.
- `CONSEQUENCES.md`: `5616a1d747edad5828f48fb36efb615cdb01d60eaf770bef9e00029d81b74e1d`.

The remaining named marked/unmarked dependencies are bound in the adjacent `MARKED_CONTEXTS_REVIEW.md`; the marked proof's final clarified hash is `885af59f322f7965bcc98b6359bff5a69bb8bec33ca22cd003eb8c0a994eefea`.

## Proof checks

The linked state reserves one bin and the split state reserves two. Both use the same variable x for the remaining ordinary bins, so the different shifts J(x+1) and K(x+2) are correct.

If the endpoints lie in one maximal module run of the reserved output path, that run is the entire path: no vertices can precede its prescribed start or follow its prescribed end. Choosing it leaves the ordinary context operator applied to the input linked profile. If the endpoints lie in different runs, the first and last runs are exactly the two marked input fragments. All other fragments are unmarked, and whole-module uniformity makes their recombination depend only on their number. These disjoint exhaustive cases prove the triangular linear law.

The reverse-anchor transitive test v->Tr_m->u has j=0 and k=x^m. Its output linked path gives the connecting operator on every monomial. Cutting free-chain runs in this open path and in disjoint odd cycles produces exactly the positive skeleton factors claimed. There is no path-label or cycle-rotation overcount. The inherited positive nabla and crossing-cycle kernels therefore prove formal-basis positivity, not merely positivity on physical profiles.

For u=v the only permitted path is the singleton, and the complement really is Q[r:T-u], including the empty module convention. No split term is appropriate. For W, the only split connection uses both outside vertices in the order b,a, so B_W=I. The two j=0 counterexamples are actual graph inputs with different output linked profiles and demonstrate why j alone cannot close.

The matching bound consumes one extra O-to-I exterior segment on the open path. Selecting one crossing edge from each disjoint segment gives 1+t_0+sum t_i<=nu. Conversely a maximum matching can be ordered along the open path, supplying nu-1 free-chain runs. Every maximal forward-shift contribution has positive leading shift coefficient, so it cannot cancel. The nu=0 case consequently has no connecting block.

When N=2nu, reaching the maximal shift leaves exactly two exterior vertices per segment. Hence each segment is one perfect-matching edge, all exterior vertices are used, and the remaining matching edges form crossing permutation cycles. For a fixed matching, selecting and ordering j>=1 open-path edges and permuting the remainder gives sum binom(nu,j) j! x^(overline(nu-j)) = nu (x+1)^(overline(nu-1)). Multiplication by the number of perfect matchings proves the claimed coefficient. Boundary nu=1 gives the constant one for the triangle.

## Source replay and independent controls

Normal and optimized source runs reproduced the frozen receipt byte for byte.

A separate literal-path/bin-assignment checker reproduces every full low-order polynomial control without importing source routines:

- 1,508 split context polynomials, checked at 5,604 integer interpolation points;
- 841 coincident-endpoint polynomials, checked at 3,857 points;
- 820 W endpoint-pair polynomial systems, checked at 4,044 points.

The checker reuses only the independent literal enumerator from the preceding marked audit. It counts directed path words by permutations, not the source rooted endpoint DP.

A second independent checker verifies all four higher-matching templates. It enumerates literal u-to-v paths, counts the complement bins through odd-cycle collections, and solves the numeric shift Vandermonde system by exact Gaussian elimination. It therefore uses neither the source endpoint DP nor its centered-moment/Lagrange extraction. The top-shift polynomials are respectively 1, 2(x+1), 4(x+1), and 3(x+1)(x+2), with independently counted perfect-matching multiplicities 1,1,2,1. Optimized and normal independent receipts agree.

Independent artifact hashes:

- `check_split_endpoints.py`: `90feaff6ad80bb3436dd8a807c64b5a34892add62cd8d509ff21cee76e4d38f7`.
- `check_split_top.py`: `d7d27ce049b3a1feb9d602d68dd1699278adc18dcb4ad22f3eda3c85291bb86e`.
- `split-independent.json`: `feaa17d3620b82e683f8947b13deff1f8b843569340cfa75c5fd1e140c29bef7`.
- `split-top-independent.json`: `e640606476f20051d0b430340e51d211aa49852270b98d1a0fa1b4dfcaaf2150`.
- Shared independent dependency `check_marked_contexts.py`: `d9446de257ba232509fc7aa0032d3c2da2dc19ddfea95fb5edcd5bda00863620`.

The all-order claims and their implementation pass at the frozen hashes above.
