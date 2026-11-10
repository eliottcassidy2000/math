# Independent audit: exact cyclic orbit costs

Verdict: **PASS after one verification correction**, applied by producer. No mathematical correction was required. This audits `../cyclic_cost/PROOF.md`, including arbitrary cyclic orbit profiles, all prime-action multiplicities, sharp sevenfold equality, and the small-q corollary. Chan–Pak completeness is not claimed.

## Proof reconstruction and hostile checks

1. For a faithful C_m action, every orbit is the regular action of a cyclic quotient C_d. Its restriction is therefore a circulant of odd order d. A vertex-transitive tournament is strong unless it is a singleton: automorphisms preserve the linear order of its strong components, contradicting vertex transitivity if there is more than one component. The lcm of the orbit lengths is m, not their sum or product.
2. OCF restricted to internal cycle packings gives the product of actual orbit-block H counts. All excluded contributions have strictly positive weight. Equality is precisely absence of mixed directed odd cycles. The assertion is stronger than simply multiplying automorphism orders, and does not conflate H/m with H/|Aut(T)|.
3. The equality proof is valid. Given v with mixed orientation into a strong block B, take a directed path from an x with v->x to a y with y->v. The first orientation change supplies a mixed triangle. Hence every outside vertex sees B uniformly. For two strong blocks, applying this to both forces their mutual orientation to be uniform. The quotient tournament then has no triangles and is transitive. Conversely an ordinal sum forbids any cycle crossing blocks. This establishes both directions, including singleton blocks.
4. Replacing each actual block count by mu_d preserves equality exactly when every block minimizes. Ordinally concatenating minimizing circulants, with simultaneous translations, realizes every specified profile whose lcm is m. This audits the full arbitrary-composite cyclic theorem, not just the p=7 specialization.
5. For a nontrivial prime action, all nonfixed orbits have size p. Thus the s-orbit formula and equality description follow for every prime p and every s>=1. For p=7 the independently reproduced minimum175 yields H/7>=175^s/7. Equality H/7=25 forces s=1; equality for a fixed larger s means an ordinal sequence of minimizing seven-blocks and singletons. Larger full automorphism groups do not invalidate the specified diagonal subgroup action.
6. The small-q corollary uses **inherited** THM-4812(c), `01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md`. Its exact threshold is23 and its exceptions are K1,C3,R5,P7. The producer uses full-Aut(T) orbits to invoke it; those induced tournaments really are vertex-transitive. Restriction of the specified cyclic A embeds it in the product of cyclic projections A_i, giving |A|<=product|A_i| and the displayed inequality in the correct direction. Every local quotient H_i/|A_i| is an integer >=1. A P7 orbit forces a factor at least189/7=27 because its cyclic subgroup orders are1,3,7. Removing it leaves a product of C3/C5 groups, hence exponent15. A cyclic subgroup then has order dividing15. All four possibilities are attained by ordinal K1/C3/R5 combinations. This audit did **not** rerun the large vertex-transitive classification underlying THM-4812(c); that dependency remains explicit.

## Independent finite verification

`check.py` here imports no producer code and counts Hamilton paths by literal recursive DFS rather than subset DP. It constructs invariant completions from directed pair orbits rather than from skew connection sets.

- All30 cyclic tournaments for orbit sizes3,5,7,9 were checked. Minima are3,15,175,3159 respectively; the order9 control is not needed by the sharp-seven theorem.
- The order7 histogram is six175 and two189. Independently verified that the six minimizing connection sets are exactly the nonzero affine multiplier images of {1,2,3}.
- All208 invariant tournaments across profiles (3,3), (3,1,1), (5,1), (3,3,1), (3,5), and (1,1,1) satisfy the product lower bound, both equality equivalences, and divisibility by the lcm. Equality counts are8,12,8,24,16,6 respectively. These include multiple nontrivial prime orbits, singleton placements, and the composite C15 action with orbit profile(3,5).
- Literal enumeration of the21 affine P7 transformations yields element orders exactly{1,3,7}, independently checking the local cyclic-group obstruction used in the small-q argument. Fullness of that automorphism group is inherited, not reclassified here.
- Both ordinary and `python3 -O` runs pass and their JSON outputs are byte-identical. `--inject-failure` replaces expected mu7 by176 and is rejected with a nonzero exit under both modes.

Reproduce from repository root:

    python3 05-knowledge/results/session238_planar_chanpak_20261110/cyclic_cost_audit/check.py
    python3 -O 05-knowledge/results/session238_planar_chanpak_20261110/cyclic_cost_audit/check.py

Append `--inject-failure` for the expected-failure controls. Output files are `results.json`, `results.optimized.json`, `failure.out`, and `failure.optimized.out` in this directory.

## Correction and frozen target

The initial producer checker used Python assertions, which disappear under optimization. This was reported before finalization; the producer replaced them with explicit RuntimeError gates and added its own corruption control. The proof itself was sound. No new all-order input is inferred from these finite tests.

Audited target SHA256:

- PROOF.md: `5c046ae29aac62934e914f7246871f7d796fa4e39e5f3016a877d64a6e6ca3a9`
- check.py: `d0bc0ce3633ebcb9ad581e7e9154bf9d5aff4425a91a0e351f90ab13c303943c`
