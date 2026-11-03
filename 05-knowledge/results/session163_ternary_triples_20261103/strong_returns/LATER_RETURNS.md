# Later source returns: two scoped obstructions

**Status.** The first observation is an all-word finite-state certificate for a specific 65-atom language. The second is an explicit 60-graph test of one kind of nonuniform strong repair. Neither excludes other source languages, source re-markings, multiple cross-edge reversals or arbitrary strong tournaments. In particular, the tournament path-count spectrum and the closed ternary source bank remain open. No novelty relative to the wider literature is asserted.

Let `(h,x,y)` be the physical observer defined in [the strong-return note](NOTE.md), and t=(2,−1,0). The two returned strong sources there have states

    R2=(285,123,225),  R4=(291,120,225).

Their exact full-state ternary targets are

    target source   3R                3R+t              3R+2t
    R2              (855,369,675)     (857,368,675)     (859,367,675)
    R4              (873,360,675)     (875,359,675)     (877,358,675).

These full-state targets ensure the exact correction t required by the inherited mixed-prefix observer; a tournament with the desired first coordinate but different observer state is not excluded.

## The inherited 65-atom source language

Start at the empty-source state `(1,0,0)` and multiply the 65 matrices from [S156's pinned `bank65.json`](https://github.com/eliottcassidy2000/math/blob/3f58f67ceef35168766783b0fa6a837a7153aeb7/05-knowledge/results/session156_spectrum_20261103/returns/bank65.json) in chronological order. This is a literal single-bridge source language; its construction and matrix calibration are inherited from [S156](https://github.com/eliottcassidy2000/math/blob/3f58f67ceef35168766783b0fa6a837a7153aeb7/05-knowledge/results/session156_spectrum_20261103/returns/NOTE.md). The copied bank file has SHA256 `26fad58d82f2428be3bba1e7b72e93e6494983eb8e33d1a3dba1f10c2c62ca20`.

All matrices have nonnegative entries, with upper-left entry at least 1. Every reachable state is nonnegative, so its first coordinate h cannot decrease. Breadth-first closure from `(1,0,0)`, discarding states with h>877, terminates at exactly **7,372** states and contains none of the six targets above. Because h never decreases, any word ending at one of those targets would remain under 877 at every prefix. Therefore **no finite word in this 65-atom language realizes any of the six exact source states**. This is not merely a bound on word length. The closure count and the h-monotonicity are both checked in `later_bank_check.py`; independent calibration of the inherited matrices is in the S156 reference.

A test of reduced arithmetic alone can miss this issue: for each modulus m=2,…,13, all six target residues occur in the 65-matrix orbit of `(1,0,0)` modulo m. This statement is a finite modular computation only, not an assertion about every modulus. Thus those congruence tests cannot establish the six exact identities; they do not make the exact exclusion global.

For any fixed common suffix W in `{b,q2}` applied to both sides, the obstruction also rules out a predecessor word V in this alphabet satisfying `v(VW)=(3R_i+j t)P_W`: both atom matrices are nonsingular, with determinants 1 and 5, so cancellation would give the already excluded `v(V)=3R_i+j t`. This applies only to a common suffix and the named source alphabet. It does not rule out a graph with a different suffix representation, a different port, or an external construction.

## One nonuniform crossing edge

The universally available marked-source 0 branch gives the nonstrong ordinal products C3⊕R2 and C3⊕R4, with counts 855 and 873 and the original source port retained. For each product reverse exactly one of its 30 cross edges, leaving all internal edges unchanged. Every resulting graph is strong: its two internally strong blocks now have edges in both directions.

The [odd-cycle-packing identity](https://github.com/eliottcassidy2000/math/blob/7589463e54e20585aff6a2e68019ee4a4af582cd/05-knowledge/results/session159_constructions_20261103/cycle_profiles/NOTE.md) gives a general lower bound. If the reversed cross edge is v→u, where v lies in X and u lies in C3, then

    H(new) ≥ 3H(X)+2H(X−v).

Indeed every old packing remains, and the new triangle u→w→v→u, with w the out-neighbor of u in C3, can be combined with every odd-cycle packing of X−v. By the same identity their total weight is 2H(X−v). Direct nine-vertex path enumeration gives

    H(R2−v): 163, 175, 95, 15, 27, 123, 87, 27, 81, 165
    H(R4−v):  97,  97, 97,159, 27,  15, 81, 51, 81, 159.

Thus the symbolic lower bound alone gives a jump of at least 30 in either case, excluding +2 and +4.

There are only ten distinct results by cyclic symmetry of C3. Listed by the R vertex label 0,…,9, exact Hamilton-path count changes are

    R2: 2162, 1000, 760, 1080, 2790, 4426, 3670, 4470, 10572, 12570
    R4: 1668, 1668, 1668, 822, 846, 1524, 4770, 6288, 10050, 11508.

Thus every one-edge nonuniform strong repair overshoots the desired +2 and +4 count changes, by a minimum jump of 760 or 822. This excludes exactly this construction, not all nonuniform repair gadgets. The general lower bound above does not yield the displayed exact jumps.

## Reproduction

Place `bank65.json` beside `later_bank_check.py` and run both `python3 later_bank_check.py` and `python3 -O later_bank_check.py`. The program verifies the pinned checksum, h-monotonicity, full bounded reachable-state closure, the six exclusions and the modular orbits. Run `g++ -std=c++17 -O2 -Wall -Wextra crossing_check.cpp -o crossing_check && ./crossing_check`; this independently constructs all 60 labelled crossing reversals, checks strong connectivity, and uses endpoint/subset dynamic programming to verify every exact count and all 20 deletion counts. Assertions remain active under optimized Python and C++ builds.
