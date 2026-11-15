# Changing core and chronology without losing the exact count

2026-11-14. Continuation of the [Collatz address and 63 packet](../chanpak_collatz_exact_addresses_20261114/README.md). This session studies which physical changes can genuinely add missing outcomes to an exact count construction. **Chan–Pak completeness remains open.** The statements below concern explicit families or the indicated chronological parameters, not a census of the full tournament spectrum.

## The cleanest constructive result: change the boundary, not the three sizes

[One-outlet transfer](ONE_OUTLET.md) lets an arbitrary tournament replace one vertex that beats only a specified exterior outlet. Its full count is a linear expression in the replacement's Hamilton count and its ordered two-cover count. With transitive replacements, four literal hosts have the resulting arms

    19·2^m−13,    18·2^m−5,    9·2^m−3,    9·2^m−1.

The first reaches 63. Changing the outside host and its width gives all three consecutive odds **9·2^m−5, 9·2^m−3, 9·2^m−1** for every m≥2, including **67, 69, 71**. This solves exact placement on those arms. The improvement is genuinely geometric: arbitrary substitution of **three** tournaments into a directed triangle (not just three transitive ones) misses 63, 67 and 71 at every order, while hitting 69. The first arm crosses the reindexed lower arm of the triple at 139 only: their exact difference is 2^m−8 when comparing 19·2^m−13 with 9·2^(m+1)−5. A single crossing is not an indefinitely overlapping corridor. A fifth marked host gives an exact [shifted-trit switch](SHIFTED_TRIT_RETURN.md): it sends the entire lower arm h=9·2^m−5 to 3h+2 without preserving the forbidden triangle bank, for example 67→203. Output and input arms differ mod 9, so the switch is not closed under repetition.

The same note exposes a norm obstruction: synchronizing the two core ratios in a singleton triangle wrap forces a scaled x²+2y². Its unique positive representation gives a target-first search for primes p and 7p and an explicit upper bound on core order. That search is exact within the synchronization hypothesis, but the compatible prime targets that pass are sparse even among primes.

## What two Hamiltonian paths actually buy

A path budget asks how many full Hamiltonian paths are needed to cover all backward arcs against one common rank.

- [Exact rectangular core](BALANCED_TWO_RUN.md): for C3[1,TT_r,TT_r] this budget is r−floor(√r) for r≥2. There is a full rectangular formula, an exact target inverse for its Mersenne-product counts, and a logarithmic fixed-budget spectrum in this family. One vertex deletion suffices to make the tournament transitive, so feedback vertex number and chronological path budget measure genuinely different things.
- [Induced hull](INDUCED_HULL.md): a tournament embeds as an induced subtournament of some k-path tournament iff it has a rank in which every backward indegree and outdegree is ≤k. Edge coloring splits backward arcs into disjoint chains, and private bridge vertices complete the chains into full paths over the same augmented vertex set. That completion need not preserve H.
- [Two-path changing-core strip](READOUT_HIERARCHY.md): a cyclic matching construction preserves κ=2 while replacing a middle part by a chosen κ≤2 core. Its exact count uses the first three ordered cover statistics. One core gives **18,893=7·2699** and admits every positive 7-adic valuation along an explicitly liftable sequence. In contrast to that local success, even allowing arbitrary ordinal sums of one-apex blocks as candidate cores gives only (X^{\log_3 2}\exp(O(\sqrt{\log X}))) different counts up to X in this strip. Neither infinitely many prime cofactors nor exact coverage follows.

## Next exact obligation

The one-outlet route gives the equation N=α+β·2^m, but the two coefficients must be realized by the **same marked exterior graph**. The three-cover route offers more degrees of freedom, but the successful 7-adic lift chooses only the valuation, not the remaining absolute cofactor. A viable continuation needs a theorem supplying joint marked profiles adapted to N, or a further contact mode that changes the absolute baseline without undoing readout readiness. Do not infer such a theorem from residue coverage, a repeated dyadic center, or induced embeddability.

## Replay

Use Python 3, with or without `-O`:

- `replay_two_run.py`: 1,764 literal path covers and independent small permutation path counts.
- `replay_induced_hull.py`: literal common-rank supertournaments and colored-chain certificates.
- `replay_one_outlet.py`: direct subset path counts, outlet arms and the order-cutoff-controlled exact spectrum of all cyclic triple wraps through count 98.
- `replay_readout_hierarchy.py`: strip counts, rank/cover certificates and finite controls for the exact 7-adic isometry.

Replays are finite checks accompanying the all-order arguments, not substitutes for them. See [Chan–Pak §6](https://doi.org/10.1016/J.TCS.2024.114776) for the conjecture.
