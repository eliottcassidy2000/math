# A constructive eventual square-sum cycle theorem with the prescribed seed

## Theorem

For every integer N ≥ 128453, there is a Hamiltonian cycle on the vertex set {1,…,N} in which every adjacent pair, including the closing pair, has square sum, and in which 482, 3362, 6242 occur consecutively in that order after choosing the cycle's orientation.

This is a computer-assisted constructive theorem proved by the finite certificates identified below and the induction in this document. The threshold is sufficient; no optimality claim is made. The ordinary interval square-sum Hamiltonicity theorem is credited to Robert Gerbicz (2018). No historical novelty or priority claim is made for the prescribed-seed result here.

## 1. Nice pairs and the three families

A nice pair at index n consists of square-sum Hamiltonian paths S0 on [1,n] and S1 on [1,n+1], such that every common vertex has the same positional parity in the two paths. Positions are one-based. The extra vertex n+1 necessarily has positional parity n+1: subtracting the common odd/even-position counts proves this directly.

The families used here impose these additional conditions:

- B: both paths start 1,3; even-length paths end 8 and odd-length paths end 10.
- Aplus: both paths start 69,127; even-length paths end 10 and odd-length paths end 8.
- Aminus: both paths start 69,127; even-length paths end 8 and odd-length paths end 10.

Both A families also contain the edge 2–7 in each path, with 2 at an odd position and 7 at an even position. The direction of this edge may differ between the two paths.

## 2. Affine blocks and exact coverage

For offsets −24 ≤ c ≤ 24, put

T(c,i)[j] = 49 S_i[j] + (−1)^j c.

Every internal square-sum edge remains a square-sum edge: its sum is multiplied by 49, which is 7². For a common vertex x, the two input paths place x at the same parity, so its lifted representative agrees in both possible blocks. As c ranges over all 49 offsets, the representatives of x fill exactly [49x−24,49x+24].

To obtain a path of length M = 49n+d, where 24 ≤ d ≤ 73, choose

i(c)=1 precisely when 49+(−1)^(n+1)c ≤ d.

These chosen blocks partition [25,M]; adding the singleton connectors 1,…,24 partitions [1,M]. Reversing any block changes no coverage. A finite glue word specifies every chosen block, its orientation, and every singleton exactly once. Numerical endpoint sums check every edge between atoms.

For paired outputs at residues d and d+1, equal singleton position parities and equal original-position shifts for each block certify the nice-pair condition. This certificate depends only on n modulo2, rather than n itself or the internal order of either input path.

## 3. Why a small repair is necessary, and why it is valid

The unique block whose input choice changes between the outputs has offset c*=48−d for even n, or c*=d−48 for odd n. Its newly added vertex has input-position parity n+1 and output-position parity 49n+d+1. Thus its shared shift must be d modulo2.

Consequently, if the changed block's shorter-input terminal port has no external neighbor except the required terminal singleton of the shorter output, the whole-block model is impossible: both would have output-position parity n+d. At parity0 residue48, the changed block c=0 has port392, whose only external neighbor is singleton8. This obstructs B and Aminus. At parity0 residue50, the changed block c=−2 has port390, whose only external neighbor is singleton10. This obstructs Aplus. The corresponding dual states are obstructed as well. These are obstructions to an intermediate architecture, not to the theorem.

Every B input has the guaranteed initial edge1–3. We permit this edge to be cut in one declared block c, producing:

- a head singleton 49−c;
- a tail, the lift of S_i with its first vertex removed, whose initial endpoint is147+c and terminal endpoint is49·end(S_i)+(−1)^(n+i)c.

The pieces remain disjoint and partition exactly the original block. The head and tail retain one shared original-position shift. If a forward tail starts at parity p, its shift is p; if reversed, its shift is p+(n+i) modulo2. The head must have parity1+shift. This is the same parity certificate for every common lifted vertex, including the separated head.

The frozen table needs only these three repairs at input parity0:

- B, residue48: cut c=0.
- Aminus, residue48: cut c=0.
- Aplus, residue50: cut c=−2.

All other input-parity0 transitions use intact blocks. For the first repair, the cut makes392–49 available; for the third, it makes390–51 available. Both sums are441. The full words, not just these local connections, have been checked.

## 4. The finite transition theorem

The attached typed table contains every family, input parity, and residue24,…,72: 3·2·49 = 294 states. It contains288 whole-block states and6 single-prefix-split states. The universal checker independently reconstructs endpoints, checks exact atom coverage, verifies every numerical connector, checks the family prefix and terminal, checks the protected2–7 edge and parities for both A families, and compares complete parity profiles across both words.

The parity-1 states are also explicit. They are obtained from the corresponding parity-0 state at residue96−d by exchanging the order of the output words and complementing each selected input index. Offsets, orientations and head atoms remain unchanged. Complementing input parity and selected index preserves every block-length parity and numerical endpoint, including those of a split tail. The complete-table checker verifies this exact dual relation in addition to replaying each state separately.

Therefore every B_n produces B_(49n+d), Aplus_(49n+d), and Aminus_(49n+d), for every24 ≤ d ≤72. These are universal implications for any admissible input pair, not just the concrete test pairs used in expansion controls.

## 5. Induction to all B and A indices

The independently checked base bundle supplies B_n for every53 ≤ n ≤2620, including both members of each nice pair. Since49·53+24=2621, strong induction now gives B_n for every n≥53: whenever n≥2621, write

m = floor((n−24)/49),  d = n−49m.

Then24 ≤ d ≤72 and53 ≤ m<n. Apply the appropriate B-to-B certificate to B_m.

For every n≥2621, the same decomposition and a B-to-A certificate produce both protected A families at index n. These transformations are one-time uses of B pairs; they do not require any shared parity coloring between the A families.

## 6. Final lift with the prescribed seed

Given N≥128453, write N=49n+d with24 ≤ d ≤72 using n=floor((N−24)/49). Then n≥2621. Put j=0 if d≤56 and j=1 if d≥57. Choose whichever A family makes its member S_j end in10: Aplus if n+j is even, and Aminus if n+j is odd.

Lift this pair using the same forced block choices. In every chosen block cut the protected edge2–7, producing two segments. Its two cut ports are98−c and343+c because the original position of2 is odd and that of7 is even. The two possible edge directions in each member give four orientation states. Together with two input parities and49 residues, this is a finite table of392 cases.

The frozen final table provides, in every case, a glue word using all98 segments and all24 singletons, beginning at1 and ending at8, and containing the connector482–3362. Closing the edge8–1 gives a square-sum Hamiltonian cycle.

Indeed, the block ending at482 has c=−8(−1)^(n+j), and its forced choice is j. The block at c=19 begins49·69−19=3362 followed by49·127+19=6242. The cut2–7 does not disturb that prefix. Thus the connector482–3362 is immediately followed by3362–6242, or the whole seed appears in the reverse orientation. Orienting the cycle appropriately gives the stated consecutive triple.

Finally,49·2621+24=128453. This proves the theorem.

## 7. Exact evidence and reproduction

Paths below are relative to this `seeded_eventual_cycle_theorem_20261019` directory. The theorem and its two dependency bundles are siblings under `05-knowledge/results/`.

1. Base bundle: `../square_sum_B_base_interval_20261019/`. Its `B_base_chunks.json` indexes11 compressed root-byte files and their SHA-256 digests. `verify_B_base_interval_20261019.py` exhaustively checks2568 nice pairs,6866832 vertices and6861696 square-root steps, with17173589 explicit checks. Its normal and optimized receipts pass.
2. Pair-transition bundle: `certified-pair-transitions/`, indexed and hashed in `certified_pair_transition_manifest.json`. `verify_pair_transition.py` is the independent universal replay checker; its concrete60/61 expansion input is frozen beside it as `low_prefix_pair_test.json`. `verify_complete_transition_table.py` checks all294 states, every witness hash and all147 duality relations. Both normal and optimized receipts pass. Discovery scripts and timeout files are not proof dependencies.
3. Final bundle: `../seeded_final_split_table_20261019/`. `verify_final_split_table_20261019.py` checks all392 conditional states with234025 explicit checks. The certificate file `final_split_glue_certificates_20261019.json` has SHA-256 `3a035b14c8caa9f1293f66bbb0c9e77e9f1e933b7e5a07e523e2a5e7f59dbffe`.
4. Constructive expansion: `construct_seeded_cycles.py` decodes B bases, recursively constructs B pairs, applies either typed pair transition, selects the A state and final cut orientations, expands the final segments, and directly checks every vertex and square edge in the resulting cycle. Its independent controls pass414343 checks in normal and optimized modes. They construct N=128453,132063,132128, exercising both repaired A transitions, and additionally expand the repaired B pair at2694/2695. The three compressed cycles are included; each canonical presentation starts482,3362,6242.

To replay the transition bundle, run `python verify_complete_transition_table.py` and `python -O verify_complete_transition_table.py` from its directory, and compare the resulting receipts. The base and final bundles have their own similarly named checkers and existing normal/optimized receipts. To reconstruct the three concrete cycles, run `python construct_seeded_cycles.py` with the listed sibling bundles retained in the relative directory layout; passing any other integer N≥128453 requests that cycle instead.

All correctness checks use explicit failures, not Python assertions removable by optimization. The proof uses the finite replay results plus the mathematical coverage, parity, induction and final-lift arguments above; it does not rely on the behavior or trustworthiness of the discovery MILP solver.
