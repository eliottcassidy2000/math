# Permutation-orbit refinement of the port-word cone

2026-10-19 UTC. The five-seam theorem improves from balance 24/55 to **236/525 = 9/20 − 1/2100**. This is a valid outer-cone certificate for actual uniform exterior-poset laws; sharpness over realizable contexts is not asserted.

## 1. Actual boundary laws have a symmetry constraint

Let U be a designated set of maximal elements of any finite exterior poset. For a subset S of U, consider extensions whose maximal terminal run of U consists precisely of S. All permutations of S occur equally often: the members of S are an antichain of maximal elements, so any permutation of those final positions is valid and preserves the preceding non-U element (if present). This is a bijection fixing the rest of the extension. Dually, initial runs from a designated set of minimal elements have equal counts under every permutation of their underlying subset.

This argument applies to the entire exterior poset, even when that exterior was built using earlier seams. It therefore constrains its true uniform extension law, not just the unweighted internal law of its terminal core.

For two maximal ports {6,8}, terminal-profile counts have the form

    (x_empty,x_6,x_8,t,t)

in profile order (empty,6,8,68,86). The two nontrivial orderings must be grouped as one generator. For minimal ports {0,3}, initial counts similarly have the form (y_empty,y_0,y_3,u,u). In the present core words all maximal elements of a nonempty left exterior are its two outgoing ports, so x_empty=0; an empty exterior is handled separately by the all-ones weight.

For each five-by-five shuffle kernel K, valid left generators are its rows K[1], K[2], K[3]+K[4], together with the all-ones empty-context vector. Valid right generators are columns 0,1,2 and their sum column 3+column 4, restricted to the four nonempty core suffix profiles. Arbitrary nonnegative combinations contain every actual left and right weighting.

## 2. Exact primitive generators and deficits

Removing duplicates up to positive scaling gives these eight left generators, in core prefix order (empty,0,3,03,30):

    (1,1,1,1,1), (1,1,2,1,2),
    (1,2,1,2,1), (1,2,2,3,3),
    (2,3,3,3,3), (2,3,3,4,4),
    (2,3,6,4,8), (2,6,3,8,4).

There are six right generators, in core suffix order (6,8,68,86):

    (1,1,1,1), (1,2,2,1), (1,3,3,1),
    (2,1,1,2), (2,2,3,3), (3,6,8,4).

Before primitive-ray normalization there are seven distinct right vectors: the additional vector (3,3,3,3) repeats the all-ones ray, and (2,6,6,2) is normalized to (1,3,3,1). Thus the refined certificate still has only 8·6=48 generator products per core.

Using the exact core profile-count matrices N and marked matrices C for event 5<7 gives:

- B8 probability range: [236/525,18/35];
- B9 probability range: [1095/2372,600/1153].

Both have balance at least 236/525. All 48 lower deficits `525 lCr − 236 lNr` and upper deficits `289 lNr − 525 lCr` are nonnegative. With the displayed primitive generators, their minima are respectively 0 and 3720 for B8, and 6641 and 15797 for B9.

The minima of the probabilities occur at left weight (1,2,1,2,1) and right weight (1,3,3,1): a factor two for a core extension starting at 0, and a factor three for one ending at 8. For B8 this gives numerator/denominator 236/525; for B9 it gives 1095/2372. Individual rays in this outer cone need not be realizable by actual exteriors, so these extrema are cone extrema, not a proved optimum over actual block-word contexts.

The bilinear nonnegative-combination argument now transports this stronger bound to every block in every word allowed by the base five-seam theorem. Its degree, width, endpoint-witness and conditional bridgeless assertions do not change.

## 3. Reproduction and scope

`verify_port_word_refinement_20261019.py` imports the unchanged standard-library base verifier, recomputes its shuffle kernels and core profile tables, constructs the permutation-orbit generator sums, normalizes the rays, and checks all 48 exact rational comparisons per core. It also checks the core profile equalities and that designated endpoints remain maximal/minimal in two mixed exteriors. Exact arrays and all individual deficits are saved in `port_word_refined_cone_20261019.json`.

The all-context claim rests on the permutation bijection and cone proof. The computations are fixed exact controls, not a census or an all-length numerical experiment.


Reproduce from the repository root with `python3 05-knowledge/experiments/verify_port_word_refinement_20261019.py`. The three standard-library verifiers are colocated for imports; published JSON receipts are in `05-knowledge/results`.
