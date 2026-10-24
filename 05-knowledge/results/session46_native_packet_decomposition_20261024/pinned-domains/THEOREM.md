# Complete common-clock domains for the expanded source supply

## Exact finite result

Starting with the 28 actual source options inherited from Session45 for targets
U=1965 and V=2234, retain every indexed owner clock common to the two saved
paths of each source. All owners at unfixed clocks must still form a permutation.

The complete pinned Hamiltonian domain has exactly two paths at each of 27
sources. Source 1923 alone has three paths; its additional path has target 1559.
There are 57 paths in total, with 57 distinct normalized black F/R resources and
228 distinct closed ports. Every source automorphism group is trivial, checked
on the original canonical matrices; hidden lifts are not additional capacity.

Exactly 56 paths have target U or V, one of each at every source. Consequently,
when the complete assignment retains the original target margins supported only
on U and V, these are genuinely complete binary source domains, not arbitrary
selected submenus. The extra path to 1559 cannot occur: its target count would
be positive where the required count is zero, and actual nonnegative assignment
counts cannot cancel that contribution.

This distinction matters. Raw pins alone do not give a binary domain at 1923.
The target-margin guard, not a silently discarded path, supplies that restriction.

## Retained data

The receipt saves every canonical source matrix, all fixed clock/owner pairs,
full path, normalized word, exact target, complete four-word resource and the
actual owner frame. The declared common roles are first owner, last owner and
the two owners of the saved indexed internal arc. Every other fixed owner clock
is retained as a private coordinate. Identifying these declared roles does not
assert that different source tournaments amalgamate into one common host.

The checker reconstructs actual K partners, proves the literal involution,
checks reflected converse and its owner-clock reversal, and verifies all original
fixed coordinates for every path. No quotient or kernel label substitutes for
an actual presentation. In this particular supply the source automorphism lift
is trivial; no all-size triviality assertion is made.

Both original targets have equal paired-charge coordinates, so their common
choice can be used by the separately audited dyadic row criterion. This domain
certificate alone is not a cancellation certificate for a collection of slots.
It does not assert every ordinary individual-port harmonic phase is unchanged.

## Proof and reproduction

For each source, fix the common clock positions and enumerate all permutations
of exactly the remaining owners. This is the full domain because a Hamiltonian
path uses each owner once. Test all seven directed consecutive arcs against the
original matrix. Compute the normalized word and use the original owner map
for its source and complemented target. Direct four-word orbit intersections
establish disjointness, not a comparison of raw lift labels.

All numerical conclusions follow by these exact finite tests. An independent
audit enumerates all 8! owner permutations per source rather than merely the
unfixed positions. The source files are pinned by SHA256 in the receipt.

Run with the original repository and inherited supply file:

    python3 -B check.py --repo REPO --supply REPO/05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json --out normal.json
    python3 -O -B check.py --repo REPO --supply REPO/05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json --out optimized.json

The development invocation may point --supply to the original separate packet.
The supplied inventory has a stable logical provenance key, so the saved
receipt is identical under either invocation and after relocation.

Normal and optimized runs pass 795 checks identically. This is an expanded
native supply certificate for selected packets, not a full-atlas rewiring.
