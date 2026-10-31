# All-size square-sum cycles containing an explicit unbounded-length path

## Theorem

For t≥0 put Q_t=49^t, C_t=(Q_t−1)/12, and

    H_t=(Q_t−C_t, 3Q_t+C_t, 6Q_t−C_t, 19Q_t+C_t).

Let P_l be H_0,H_1,…,H_(l−1), followed by the first two entries of H_l. Then P_l has4l+2 vertices and the P_l are nested prefixes of one simple infinite square-sum ray.

For every l≥1 and every integer

    N ≥ T_l=(107·49^l−1)/2,

there is a square-sum Hamiltonian cycle on[1,N] whose oriented traversal begins with the **original, unscaled** path P_l. It ends at8 for even N and15 for odd N, closing to1 by9=8+1 or16=15+1.

Examples:

- P1=(1,3,6,19,45,151), for every N≥2621.
- P2=(1,3,6,19,45,151,290,935,2201,7403), for every N≥128453.
- P3 extends P2 by14206,45819,107845,362751, for every N≥6294221.

Every finite prefix of this explicit infinite path therefore extends to a Hamiltonian cycle on every sufficiently large initial interval. This is **not** arbitrary prescribed-path embedding. Ordinary square-sum Hamiltonicity is credited to Gerbicz2018, and no historical-priority claim is made here.

## Minimal dependency boundary

This package uses only:

1. the existing certified B base pairs at indices53 through2620, stored in11 root-byte chunks in the repository-relative directory `05-knowledge/results/square_sum_B_base_interval_20261019/`;
2. the new98 paired B transitions in `paired/`;
3. the new98 single-output terminal-cycle transitions in `cycles/`.

`BASE_DEPENDENCIES.json` freezes the base manifest and every compressed/decoded chunk hash. The data remain in their existing directory; this package does not duplicate them. **No old B transition, A transition,392-state final table, external C implementation or unlicensed source code is needed for this core theorem.** Earlier2303-copy and depth/width corollaries are separate results with additional dependencies and are not silently included in this package.

## Proof

A B_n pair consists of Hamiltonian square-sum paths on[1,n] and[1,n+1], both starting1,3, with equal one-based positional parity on every common vertex. An even-length member ends8; an odd-length member ends10. The11 base chunks directly certify such pairs for every53≤n≤2620.

For offset c∈[−24,24], the affine block maps an input vertex v in position j to49v+(−1)^j c. Internal square sums multiply by49. Input compatibility makes common representatives independent of the chosen input member. Choosing the longer member exactly when49+(−1)^(n+1)c≤d partitions[25,49n+d]; singleton1,…,24 complete the interval. The extra input vertex n+1 has positional parity n+1 by subtraction of the shared parity counts, which proves this selector rule.

A prefix split replaces a block by its head49−c and tail starting147+c. Exact ownership is unchanged. The paired checker verifies every atom, forced input, numerical connector, complete traversal, endpoint and parity profile. In particular, every separated head and tail has the same original-to-output parity shift, and the complete profiles agree between outputs. Thus the new paired table maps every weak B input pair to another B pair, for every input parity and every residue24 through72. Each output starts with

    singleton1, singleton3, singleton6, singleton19, forward intact block4.

Define B(n) deterministically: use the frozen base pair when53≤n≤2620; otherwise set m=floor((n−24)/49), d=n−49m, recursively build B(m), and apply the **new** paired state. For n≥2621,53≤m<n and24≤d≤72. Strong induction proves a B pair for every n≥53 without any old transition table.

The stronger invariant is that B(n) has common prefix P_l whenever n≥T_l. It holds for l=0 because P0=(1,3). Since T_(l+1)=49T_l+24, the quotient m is at least T_l whenever n≥T_(l+1). The fixed five-atom start then transforms P_l into

    P_(l+1)=(1,3,6,19) concatenated with (49p_j+(−1)^j4)_j.

This is exactly the stated closed family because C_(t+1)=49C_t+4. Its internal H_t edge sums are4Q_t,9Q_t,25Q_t and its next-block boundary sum is64Q_t. Positivity and simplicity follow recursively: the four new labels lie in[1,24], while each lifted previous vertex lies in its own disjoint band above24. Nesting follows from the recurrence.

The new terminal consumer uses the same forced blocks and interval partition, begins with the same five atoms, and supplies **one cycle word** for every input parity/residue. It verifies exact symbolic length49n+d, all ownership, all connector squares, and the closing edge from8 or15 to1. Output paired-parity and head/tail-shift obligations are intentionally absent: this is a terminal consumer and its output is never reused as a compatible pair. Input compatibility remains essential for coverage.

Given N≥T_l,l≥1, set m=floor((N−24)/49),d=N−49m. Then m≥T_(l−1). The compatible B(m) construction supplies P_(l−1), and the terminal converter produces a cycle beginning with P_l. This proves all sizes, not merely descendants of a particular finite base.

## Verification and reproduction

Everything here uses Python's standard library. Correctness checks use explicit failures rather than removable assertions.

If this package is installed beside the existing base directory, run:

    python verify_core.py > core.normal.json
    python -O verify_core.py > core.optimized.json
    python construct_cycle.py 2621

Otherwise give the base directory explicitly:

    python verify_core.py --base-dir /path/to/math/05-knowledge/results/square_sum_B_base_interval_20261019
    python construct_cycle.py 6294221 --base-dir /path/to/the/base/directory

`verify_core.py` verifies the frozen own-file hashes, all2568 base pairs, all98 new paired states,49 exact paired dualities and all98 terminal states. It then constructs and stream-checks cycles at2621,2622,128453,128454 and6294221, including their complete unscaled prefixes. `construct_cycle.py` stores no large cycle; it returns a deterministic decimal-lines stream hash. Construction and full verification require O(N) time and memory.

The new tables have separately passed independent checkers with no producer imports. Independent audits cover the multi-split paired semantics, all2568 bases, the terminal consumer without inappropriate output-pair requirements, direct expansions at both input parities, hostile mutations and the dependency-reduced induction. The included `audit_independent.py` replays the minimal core without producer imports:106,801,611 checks in each mode, byte-identical receipts. See `INDEPENDENT_REPRODUCE.md` for its portable commands.

The constructor is intentionally separate from discovery. No MILP solver or search log is a proof dependency. The original frozen certificate bytes are preserved in this publication package.
